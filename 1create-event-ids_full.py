import csv
import json
from collections import defaultdict

class Utils:

    @staticmethod
    def Safeint(x, default=0):
        try:
            return int(float(x))
        except:
            return default

    @staticmethod
    def Normalizestation(s):
        return (s or "").strip()

    @staticmethod
    def Dirkey(d):
        return (d or "").strip().upper()

    @staticmethod
    def Typekey(t):
        return (t or "").strip().lower()

class Event:
    def __init__(self, eventid, station, eventtype):
        self.id = eventid
        self.station = station
        self.type = eventtype  # "arrival" or "departure"


class Station:
    def __init__(self, name):
        self.name = name
        self.arrivals = []
        self.departures = []

    def addevent(self, event):
        if event.type == "arrival":
            self.arrivals.append(event)
        else:
            self.departures.append(event)


class Service:
    def __init__(self, srnum, starttime, direction, servicetype, patternid):
        self.srnum = srnum
        self.starttime = starttime
        self.direction = direction
        self.servicetype = servicetype
        self.patternid = patternid

        self.events = []
        self.segments = {}
        self.Firstdep = None
        self.Lastarr = None

    def addevent(self, event):
        self.events.append(event)

    def addsegment(self, fromevent, toevent, traveltime):
        key = f"{fromevent.id}_{toevent.id}"
        self.segments[key] = traveltime

    def Todict(self):   # Converts the Service object into a dictionary.
        return {
            "start_time": self.starttime,
            "Dir": self.direction,
            "Type": self.servicetype,
            "PatNum": self.patternid,
            "segments": self.segments,
            "first_dep": self.Firstdep,
            "last_arr": self.Lastarr
        }

class MILPCSVPreprocessor:

    Departbase = 20000
    Arrbase = 100

    Outputup = "1-o-event-ids_DOWN2.csv"
    Outputdown = "1-o-event-ids_UP2.csv"
    Outputjson = "milp_preprocessed2.json"
    Outputdump = "constraints_dump_rowwise2.txt"

    def __init__(self, servicescsv, patternscsv, supplycsv):

        self.servicescsv = servicescsv
        self.patternscsv = patternscsv
        self.supplycsv = supplycsv

        self.nextdepid = self.Departbase
        self.nextarrid = self.Arrbase

        self.stations = {}
        self.services = []
        self.patterns = {}

        self.Stationseqdown = []
        self.Stationsequp = []

        # MILP structures
        self.turnaround = defaultdict(lambda: {
            "fast": {
                "UP": {"dep": [], "arr": []},
                "DOWN": {"dep": [], "arr": []}
            },
            "slow": {
                "UP": {"dep": [], "arr": []},
                "DOWN": {"dep": [], "arr": []}
            }
        })

        self.Routedepevents = defaultdict(list)

        self.headway = {
            "UP": {"fast": defaultdict(list), "slow": defaultdict(list)},
            "DOWN": {"fast": defaultdict(list), "slow": defaultdict(list)}
        }


    def readcsvdict(self, path):
        with open(path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            return [{k.strip(): (v.strip() if v else "") for k, v in r.items()} for r in reader]

    def readcsvrows(self, path):
        with open(path, newline="", encoding="utf-8") as f:
            return list(csv.reader(f))


    def loadsupply(self):
        rows = self.readcsvrows(self.supplycsv)

        for row in rows:
            if len(row) >= 2 and row[1].strip().lower() == "station":
                for s in row[2:]:
                    s = Utils.Normalizestation(s)
                    if not s or "total" in s.lower():
                        break
                    self.Stationseqdown.append(s)
                break

        self.Stationsequp = list(reversed(self.Stationseqdown))

        for st in self.Stationseqdown:
            self.stations[st] = Station(st)

    def loadpatterns(self):
        rows = self.readcsvrows(self.patternscsv)
        header = rows[0]
        idx_pattern_id = header.index("Pattern_ID")

        for row in rows[1:]:
            pat_id = row[idx_pattern_id].strip()
            if not pat_id:
                continue

            segs = []
            for i, col in enumerate(header):
                if col.startswith("Major_Segment_"):
                    seg = row[i].strip()
                    time_col = col.replace("Major_Segment_", "Time_")
                    dt = 0
                    if time_col in header:
                        j = header.index(time_col)
                        dt = Utils.Safeint(row[j])
                    if "-" in seg:
                        a, b = seg.split("-", 1)
                        segs.append((a.strip(), b.strip(), dt))

            self.patterns[pat_id] = segs


    def processservices(self):

        rawservices = self.readcsvdict(self.servicescsv)

        for s in rawservices:

            srnum = Utils.Safeint(s.get("SrNum"))
            starttime = Utils.Safeint(s.get("Time"))
            d = Utils.Dirkey(s.get("Dir"))
            tp = Utils.Typekey(s.get("Type"))
            pat = s.get("PatNum", "").strip()

            pattern = self.patterns.get(pat)
            if not pattern:
                continue

            service = Service(srnum, starttime, d, tp, pat)

            stations = [pattern[0][0]] + [b for _, b, _ in pattern]
            traveltimes = [t for _, _, t in pattern]

            # FIRST DEPARTURE
            firststation = self.stations[stations[0]]
            depevent = Event(self.nextdepid, firststation, "departure")
            self.nextdepid += 1

            service.Firstdep = depevent.id
            service.addevent(depevent)
            firststation.addevent(depevent)

            # HEADWAY add
            self.headway[d][tp][stations[0]].append(depevent.id)

            # DISTRIBUTION add
            route_key = f"{stations[0]}-{stations[-1]}-{tp}-{d}"
            self.Routedepevents[route_key].append((starttime, depevent.id))

            Previousevent = depevent

            # Traverse pattern
            for i in range(1, len(stations)):

                stationobj = self.stations[stations[i]]

                # ARRIVAL
                arrevent = Event(self.nextarrid + 1, stationobj, "arrival")
                self.nextarrid = arrevent.id

                service.addevent(arrevent)
                stationobj.addevent(arrevent)
                service.addsegment(Previousevent, arrevent, traveltimes[i-1])

                Previousevent = arrevent
                service.Lastarr = arrevent.id

                # INTERMEDIATE DEPARTURE
                if i < len(stations) - 1:
                    depevent = Event(self.nextdepid, stationobj, "departure")
                    self.nextdepid += 1

                    service.addevent(depevent)
                    stationobj.addevent(depevent)
                    service.addsegment(Previousevent, depevent, 0)

                    # HEADWAY add
                    self.headway[d][tp][stations[i]].append(depevent.id)

                    Previousevent = depevent

            # TURNAROUND primitive sets
            origin = stations[0]
            dest = stations[-1]

            self.turnaround[origin][tp][d]["dep"].append(service.Firstdep)
            if service.Lastarr is not None:
                self.turnaround[dest][tp][d]["arr"].append(service.Lastarr)

            self.services.append(service)

    def exportcsv(self):

        def makeheader(stations):
            h = ["SrNum", "Time", "Type", "Dir", "PatNum", "From", "To"]
            for st in stations:
                h += [st + "a", st + "d"]
            return h

        downheader = makeheader(self.Stationseqdown)
        upheader = makeheader(self.Stationsequp)

        downrows = [downheader]
        uprows = [upheader]

        for service in self.services:

            header = downheader if service.direction == "DOWN" else upheader
            row = {h: "" for h in header}

            row.update({
                "SrNum": service.srnum,
                "Time": service.starttime,
                "Type": service.servicetype,
                "Dir": service.direction,
                "PatNum": service.patternid,
                "From": service.events[0].station.name,
                "To": service.events[-1].station.name
            })

            for e in service.events:
                if e.type == "arrival":
                    row[e.station.name + "a"] = e.id
                else:
                    row[e.station.name + "d"] = e.id

            if service.direction == "DOWN":
                downrows.append([row[h] for h in downheader])
            else:
                uprows.append([row[h] for h in upheader])

        with open(self.Outputdown, "w", newline="", encoding="utf-8") as f:
            csv.writer(f).writerows(downrows)

        with open(self.Outputup, "w", newline="", encoding="utf-8") as f:
            csv.writer(f).writerows(uprows)


    def exportjsonanddump(self):

        # DISTRIBUTION
        distributionmap = {}
        for routekey, items in self.Routedepevents.items():
            items.sort(key=lambda x: x[0])
            depids = [eid for _, eid in items]
            if len(depids) >= 2:
                distributionmap[routekey] = depids

        # convert defaultdict to dict
        def normaldict(x):
            if isinstance(x, defaultdict):
                x = dict(x)
            if isinstance(x, dict):
                return {k: normaldict(v) for k, v in x.items()}
            return x

        milp_json = {
            "station_sequence": {
                "DOWN": self.Stationseqdown,
                "UP": self.Stationsequp
            },
            "services": {str(s.srnum): s.Todict() for s in self.services},
            "turnaround": dict(self.turnaround),
            "distribution": distributionmap,
            "headway": normaldict(self.headway)
        }

        with open(self.Outputjson, "w", encoding="utf-8") as f:
            json.dump(milp_json, f, indent=2)

        with open(self.Outputdump, "w", encoding="utf-8") as f:
            for s in self.services:
                f.write(f"SERVICE {s.srnum}: {s.segments}\n")

    # ========================================================

    def run(self):
        self.loadsupply()
        self.loadpatterns()
        self.processservices()
        self.exportcsv()
        self.exportjsonanddump()

        print("Created:")
        print("  ", self.Outputdown)
        print("  ", self.Outputup)
        print("  ", self.Outputjson)
        print("  ", self.Outputdump)


# ============================================================

if __name__ == "__main__":
    MILPCSVPreprocessor(
        "Oservices_data.csv",
        "OWR-patterns-consolidated-sequential.csv",
        "WR_supply.csv"
    ).run()
