import csv
import json
from collections import defaultdict

# INPUTS
SERVICES_CSV = "Oservices_data.csv"
PATTERNS_CSV = "OWR-patterns-consolidated-sequential.csv"
SUPPLY_CSV   = "WR_supply.csv"

# OUTPUTS
OUT_DOWN_EVENTS_CSV = "1-o-event-ids_DOWN.csv"
OUT_UP_EVENTS_CSV   = "1-o-event-ids_UP.csv"
OUT_JSON = "milp_preprocessed.json"
OUT_ROW_DUMP = "constraints_dump_rowwise.txt"

DEP_BASE = 20000
ARR_BASE = 100

# UTILS

def read_csv_dict(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append({k.strip(): (v.strip() if v else "") for k, v in r.items()})
    return rows

def read_csv_rows(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.reader(f))

def safe_int(x, default=0):
    try:
        return int(float(x))
    except:
        return default

def normalize_station(s):
    return (s or "").strip()

def dir_key(d):
    return (d or "").strip().upper()

def type_key(t):
    return (t or "").strip().lower()


# SECTION 1: READ INPUTS
services = read_csv_dict(SERVICES_CSV)
patterns_raw = read_csv_rows(PATTERNS_CSV)
supply_raw = read_csv_rows(SUPPLY_CSV)

# ============================================================
# SECTION 2: EXTRACT STATION SEQUENCE FROM SUPPLY
# ============================================================
def extract_station_sequence_from_supply(supply_rows):
    station_seq = []
    for row in supply_rows:
        if len(row) >= 2 and row[1].strip().lower() == "station":
            for s in row[2:]:
                s = normalize_station(s)
                if not s:
                    break
                if "total" in s.lower():
                    break
                station_seq.append(s)
            break
    return station_seq

STATION_SEQ_DOWN = extract_station_sequence_from_supply(supply_raw)
STATION_SEQ_UP = list(reversed(STATION_SEQ_DOWN))

if not STATION_SEQ_DOWN:
    raise RuntimeError("Station sequence could not be extracted from WR_supply.csv")

# ============================================================
# SECTION 3: PARSE PATTERNS
# ============================================================
def parse_patterns(pattern_rows):
    header = pattern_rows[0]
    idx_pattern_id = header.index("Pattern_ID")

    patterns = {}
    for row in pattern_rows[1:]:
        if len(row) <= idx_pattern_id:
            continue
        pat_id = row[idx_pattern_id].strip()
        if not pat_id:
            continue

        segs = []
        for i, col in enumerate(header):
            if col.startswith("Major_Segment_"):
                seg = row[i].strip() if i < len(row) else ""
                time_col = col.replace("Major_Segment_", "Time_")
                if time_col in header:
                    j = header.index(time_col)
                    dt = row[j].strip() if j < len(row) else ""
                else:
                    dt = ""

                if seg and "-" in seg:
                    a, b = seg.split("-", 1)
                    segs.append((normalize_station(a), normalize_station(b), safe_int(dt, 0)))

        patterns[pat_id] = segs
    return patterns

PATTERNS = parse_patterns(patterns_raw)

def build_station_path_from_pattern(pat_segments, direction):
    if not pat_segments:
        return [], []

    # Build physical path exactly as in patterns (DRD->CCG for UP, CCG->DRD for DOWN)
    stations = [pat_segments[0][0]]
    travel = []
    for a, b, t in pat_segments:
        stations.append(b)
        travel.append(t)

    return stations, travel


# ============================================================
# SECTION 4: HEADERS 
# ============================================================
def make_header(station_seq):
    header = ["SrNum", "Time", "Type", "Dir", "PatNum", "From", "To"]
    for st in station_seq:
        header.append(st + "a")  # arrival
        header.append(st + "d")  # departure
    return header

DOWN_HEADER = make_header(STATION_SEQ_DOWN)
UP_HEADER   = make_header(STATION_SEQ_UP)

down_rows = [DOWN_HEADER]
up_rows   = [UP_HEADER]

# ============================================================
# SECTION 5: MAIN PREPROCESSING
# ============================================================
next_dep_id = DEP_BASE
next_arr_id = ARR_BASE

# turnaround station -> list of "dep_arr" strings
# turnaround primitive sets: terminal -> Dir -> {dep:[], arr:[]}
# turnaround primitive sets: terminal -> type -> Dir -> {dep:[], arr:[]}
turnaround = defaultdict(lambda: {
    "fast": {
        "UP":   {"dep": [], "arr": []},
        "DOWN": {"dep": [], "arr": []}
    },
    "slow": {
        "UP":   {"dep": [], "arr": []},
        "DOWN": {"dep": [], "arr": []}
    }
})


# distribution: route_key -> list of (start_time, dep_event_id)
route_dep_events = defaultdict(list)

# headway dict: Dir -> Type -> Station -> [dep_ids]
headway = {
    "UP": {"fast": defaultdict(list), "slow": defaultdict(list)},
    "DOWN": {"fast": defaultdict(list), "slow": defaultdict(list)}
}

milp_json = {
    "station_sequence": {"DOWN": STATION_SEQ_DOWN, "UP": STATION_SEQ_UP},
    "services": {},          # SrNum -> dict
    "turnaround": {},        # terminal_station -> list of "dep_arr"
    "distribution": {},      # route_key -> list of "dep_dep"
    "headway": {}            # to be filled
}

def init_service_row(srnum, start_time, tp, d, pat, fr, to, header):
    row = {h: "" for h in header}
    row["SrNum"] = str(srnum)
    row["Time"] = str(start_time)
    row["Type"] = tp
    row["Dir"] = d
    row["PatNum"] = pat
    row["From"] = fr
    row["To"] = to
    return row

for s in services:
    srnum = safe_int(s.get("SrNum", ""), 0)
    start_time = safe_int(s.get("Time", ""), 0)
    d = dir_key(s.get("Dir", ""))
    tp = type_key(s.get("Type", ""))
    pat = s.get("PatNum", "").strip()
    fr = normalize_station(s.get("From", ""))
    to = normalize_station(s.get("To", ""))

    segs = PATTERNS.get(pat, [])
    stations, travel_times = build_station_path_from_pattern(segs, d)

    if not stations:
        stations = [fr, to] if fr and to and fr != to else [fr]
        travel_times = [0] * max(0, len(stations)-1)

    # Fix endpoints if mismatch
    if fr and stations and stations[0] != fr:
        stations[0] = fr
    if to and stations and stations[-1] != to:
        stations[-1] = to

    # Create events
    offset = 0
    segments = {}  # "id1_id2": dt

    # dep at first station
    dep0 = next_dep_id
    next_dep_id += 1

    # headway add
    if d in headway and tp in headway[d]:
        headway[d][tp][stations[0]].append(dep0)

    prev_event_id = dep0

    # distribution
    route_key = f"{stations[0]}-{stations[-1]}-{tp}-{d}"
    route_dep_events[route_key].append((start_time, dep0))

    # prepare service CSV row
    if d == "DOWN":
        row = init_service_row(srnum, start_time, tp, d, pat, stations[0], stations[-1], DOWN_HEADER)
        row[stations[0] + "d"] = str(dep0)  # first departure
    else:
        row = init_service_row(srnum, start_time, tp, d, pat, stations[0], stations[-1], UP_HEADER)
        row[stations[0] + "d"] = str(dep0)

    last_arr = None

    for i in range(1, len(stations)):
        dt = travel_times[i-1] if i-1 < len(travel_times) else 0
        offset += dt

        # arrival
        arr_id = next_arr_id + 1
        next_arr_id = arr_id
        segments[f"{prev_event_id}_{arr_id}"] = dt
        last_arr = arr_id

        # fill arrival column
        row[stations[i] + "a"] = str(arr_id)

        prev_event_id = arr_id

        # intermediate departure
        if i < len(stations) - 1:
            dep_id = next_dep_id
            next_dep_id += 1
            segments[f"{arr_id}_{dep_id}"] = 0
            prev_event_id = dep_id

            # headway add
            if d in headway and tp in headway[d]:
                headway[d][tp][stations[i]].append(dep_id)

            # fill dep column
            row[stations[i] + "d"] = str(dep_id)

    # turnaround store at terminal
    # turnaround primitive sets at terminal (NO PAIRS)
    
    # true terminal depends on direction
    # true terminal depends on direction
    # ---------------------------------------------------
# TURNAROUND PRIMITIVES (NO PAIRING, PURE TERMINALS)
# ---------------------------------------------------

    origin = stations[0]
    dest   = stations[-1]

    # departure belongs to ORIGIN terminal (by type + direction)
    turnaround[origin][tp][d]["dep"].append(dep0)

    # arrival belongs to DESTINATION terminal (by type + direction)
    if last_arr is not None:
        turnaround[dest][tp][d]["arr"].append(last_arr)




    # add row to correct csv
    if d == "DOWN":
        down_rows.append([row[h] for h in DOWN_HEADER])
    elif d == "UP":
        up_rows.append([row[h] for h in UP_HEADER])

    # store service in JSON
    milp_json["services"][str(srnum)] = {
        "start_time": start_time,
        "Dir": d,
        "Type": tp,
        "PatNum": pat,
        "From": stations[0],
        "To": stations[-1],
        "segments": segments,
        "first_dep": dep0,
        "last_arr": last_arr
    }

# ============================================================
# SECTION 6: TURNAROUND JSON
# ============================================================
# ============================================================
# SECTION 6: TURNAROUND JSON (primitive sets only)
# ============================================================
milp_json["turnaround"] = {st: dict(v) for st, v in turnaround.items()}


# ============================================================
# SECTION 7: DISTRIBUTION JSON
# ============================================================
# ============================================================
# SECTION 7: DISTRIBUTION JSON (ordered first departures only)
# ============================================================
distribution_map = {}
for route_key, items in route_dep_events.items():
    items.sort(key=lambda x: x[0])  # sort by start time
    # store ONLY ordered departure ids (no pairing)
    dep_ids = [eid for _, eid in items]
    if len(dep_ids) >= 2:
        distribution_map[route_key] = dep_ids

milp_json["distribution"] = distribution_map


# ============================================================
# SECTION 8: HEADWAY JSON (convert defaultdict to normal dict)
# ============================================================
def normal_dict(x):
    if isinstance(x, defaultdict):
        x = dict(x)
    if isinstance(x, dict):
        return {k: normal_dict(v) for k, v in x.items()}
    return x

milp_json["headway"] = normal_dict(headway)

# ============================================================
# SECTION 9: WRITE OUTPUT FILES
# ============================================================
with open(OUT_DOWN_EVENTS_CSV, "w", newline="", encoding="utf-8") as f:
    csv.writer(f).writerows(down_rows)

with open(OUT_UP_EVENTS_CSV, "w", newline="", encoding="utf-8") as f:
    csv.writer(f).writerows(up_rows)

with open(OUT_JSON, "w", encoding="utf-8") as f:
    json.dump(milp_json, f, indent=2)

# ============================================================
# SECTION 10: WRITE ROW-WISE READABLE DUMP (for easy reading)
# ============================================================
with open(OUT_ROW_DUMP, "w", encoding="utf-8") as f:
    f.write("===== SERVICES SEGMENTS (Row-wise) =====\n")
    for sr in sorted(milp_json["services"].keys(), key=lambda z: int(z)):
        segs = milp_json["services"][sr]["segments"]
        seg_str = ", ".join([f"{k}={v}" for k, v in segs.items()])
        f.write(f"SERVICE {sr}: {seg_str}\n")

    f.write("\n===== TURNAROUND (terminal_station -> dep_arr) =====\n")
    for st in sorted(milp_json["turnaround"].keys()):
        f.write(f"TURNAROUND {st}: {', '.join(milp_json['turnaround'][st])}\n")

    f.write("\n===== DISTRIBUTION (route -> ordered dep ids) =====\n")
    for rk in sorted(milp_json["distribution"].keys()):
        ids = milp_json["distribution"][rk]
        f.write(f"DISTRIBUTION {rk}: {', '.join(str(x) for x in ids)}\n")

    f.write("\n===== HEADWAY (Dir Type Station -> dep_ids) =====\n")
    for d in ("UP", "DOWN"):
        for tp in ("fast", "slow"):
            if d not in milp_json["headway"] or tp not in milp_json["headway"][d]:
                continue
            for st in sorted(milp_json["headway"][d][tp].keys()):
                ids = milp_json["headway"][d][tp][st]
                if ids:
                    f.write(f"HEADWAY {d} {tp} {st}: {', '.join(str(x) for x in ids)}\n")

print("Created:")
print("  ", OUT_DOWN_EVENTS_CSV)
print("  ", OUT_UP_EVENTS_CSV)
print("  ", OUT_JSON)
print("  ", OUT_ROW_DUMP)
