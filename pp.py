import json, os
import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# PARAMETERS
# ============================================================

class Param:

    BIG_M = 2000
    MIN_HEADWAY = 3.0

    TURN_LB = 5.0
    TURN_UB = 10

    RAKE_WEIGHT = 1000000.0
    DELTA_WEIGHT = 1.0

    TIME_LIMIT = 600
    OUTPUT_DIR = "outputs_json"

    EVENT_COLS = [
        'CCGa','CCGd','MCTa','MCTd','DDRa','DDRd','BAa','BAd','ADHa','ADHd',
        'GMNa','GMNd','BVIa','BVId','BYRa','BYRd','BSRa','BSRd','VRa','VRd','DRDa','DRDd'
    ]

# ============================================================
# JSON LOADER
# ============================================================

class JsonLoader:

    def __init__(self,fname):

        with open(fname) as f:
            self.J=json.load(f)

    def services(self):

        return [
            v for v in self.J["services"].values()
            if v.get("first_dep") is not None and v.get("last_arr") is not None
        ]

    def events(self):

        ev=set()

        for v in self.J["services"].values():

            if v.get("first_dep"): ev.add(int(v["first_dep"]))
            if v.get("last_arr"): ev.add(int(v["last_arr"]))

            for k in v.get("segments",{}):

                a,b=k.split("_")
                ev.add(int(a))
                ev.add(int(b))

        return sorted(ev)

    def turnaround_pairs(self):

        P=[]
        T=self.J["turnaround"]

        for st in T:

            for typ in ["fast","slow"]:

                up_arr=T[st][typ]["UP"]["arr"]
                down_dep=T[st][typ]["DOWN"]["dep"]

                down_arr=T[st][typ]["DOWN"]["arr"]
                up_dep=T[st][typ]["UP"]["dep"]

                for a in up_arr:
                    if a is None: continue
                    for d in down_dep:
                        if d is None: continue
                        P.append((int(a),int(d)))

                for a in down_arr:
                    if a is None: continue
                    for d in up_dep:
                        if d is None: continue
                        P.append((int(a),int(d)))

        return P

    def headway_pairs(self):

        P=[]

        for d in self.J["headway"]:
            for typ in self.J["headway"][d]:
                for st,ids in self.J["headway"][d][typ].items():
                    for i in range(len(ids)):
                        for j in range(i+1,len(ids)):
                            P.append((int(ids[i]),int(ids[j])))

        return P

    def distribution(self):

        return self.J["distribution"]

# ============================================================
# MILP MODEL
# ============================================================

class Scheduler:

    def __init__(self,loader):

        self.loader=loader
        self.services=loader.services()
        self.events=loader.events()
        self.turn=loader.turnaround_pairs()
        self.head=loader.headway_pairs()
        self.dist=loader.distribution()

        self.m=pyo.ConcreteModel()

    def build(self):

        m=self.m

        m.EVENTS=pyo.Set(initialize=self.events)
        m.LINKS=pyo.Set(initialize=self.turn,dimen=2)
        m.HW=pyo.Set(initialize=self.head,dimen=2)

        m.t=pyo.Var(m.EVENTS,within=pyo.NonNegativeReals)
        m.delta=pyo.Var(within=pyo.NonNegativeReals)

        m.X=pyo.Var(m.LINKS,within=pyo.Binary)
        m.p=pyo.Var(m.HW,within=pyo.Binary)

        ARR=set()
        DEP=set()

        for s in self.services:

            DEP.add(int(s["first_dep"]))
            ARR.add(int(s["last_arr"]))

        m.source=pyo.Var(DEP,within=pyo.Binary)
        m.sink=pyo.Var(ARR,within=pyo.Binary)

        m.trav=pyo.ConstraintList()
        m.headway=pyo.ConstraintList()
        m.rake=pyo.ConstraintList()
        m.distc=pyo.ConstraintList()
        m.start = pyo.ConstraintList() 

        # traversal

        for s in self.services:
    
            for k,tt in s["segments"].items():

                a,b=k.split("_")

                a=int(a)
                b=int(b)

                m.trav.add(m.t[b]-m.t[a]==tt)
                m.trav.add(m.t[b]>=m.t[a])

        for s in self.services:
            m.start.add(m.t[int(s["first_dep"])] >= 480)
            m.start.add(m.t[int(s["first_dep"])] <= 660)

        # headway

        for i,j in m.HW:

            m.headway.add(m.t[j]-m.t[i]+(Param.BIG_M+Param.MIN_HEADWAY)*m.p[i,j]>=Param.MIN_HEADWAY)
            m.headway.add(m.t[j]-m.t[i]+(Param.BIG_M+Param.MIN_HEADWAY)*m.p[i,j]<=Param.BIG_M)

        # rake conservation

        for j in ARR:

            m.rake.add(sum(m.X[i,j] for i,jj in self.turn if jj==j)+m.sink[j]==1)

        for i in DEP:

            m.rake.add(sum(m.X[ii,j] for ii,j in self.turn if ii==i)+m.source[i]==1)

        m.rake.add(sum(m.source[i] for i in DEP)==sum(m.sink[j] for j in ARR))

        # turnaround

        for a,d in self.turn:

            m.rake.add(m.t[d]-m.t[a]>=Param.TURN_LB-Param.BIG_M*(1-m.X[a,d]))
            m.rake.add(m.t[d]-m.t[a]<=Param.TURN_UB+Param.BIG_M*(1-m.X[a,d]))
        
        # distribution constraints
        for ids in self.dist.values():
            if len(ids) > 1:

                ideal = 60.0 / len(ids)

                for u in range(len(ids)-1):

                    a = int(ids[u])
                    b = int(ids[u+1])

                    m.distc.add(m.t[b] - m.t[a] >= ideal - m.delta)
                    m.distc.add(m.t[b] - m.t[a] <= ideal + m.delta)

        # objective

        m.obj=pyo.Objective(

            expr=
            Param.RAKE_WEIGHT*(sum(m.source[i] for i in DEP)+sum(m.sink[j] for j in ARR))
            + Param.DELTA_WEIGHT*m.delta,

            sense=pyo.minimize)

        return m

# ============================================================
# SOLVER
# ============================================================

class Optimizer:

    def __init__(self,m): self.m=m

    def solve(self):

        s=SolverFactory("gurobi")

        s.options["TimeLimit"]=Param.TIME_LIMIT
        s.options["MIPGap"]=0.001
        s.options["Threads"]=0

        return s.solve(self.m,tee=True)

    def times(self):

        return {i:float(self.m.t[i].value) for i in self.m.EVENTS}

    def rake_links(self):

        return [(i,j) for (i,j) in self.m.LINKS if self.m.X[i,j].value>0.5]

    def rake_count(self):

        return int(sum(self.m.source[i].value for i in self.m.source))

# ============================================================
# OUTPUT FUNCTIONS
# ============================================================

def save_rake_links(active_links,model):

    rows=[]

    for dep,arr in active_links:

        rows.append({
            "dep_event":dep,
            "arr_event":arr,
            "turnaround_time":model.t[dep].value-model.t[arr].value
        })

    pd.DataFrame(rows).to_csv("rake_linkages.csv",index=False)

def replace_csv(times, up="1-o-event-ids_UP2.csv", down="1-o-event-ids_DOWN2.csv"):
    
    os.makedirs(Param.OUTPUT_DIR, exist_ok=True)

    def run(infile, outfile):

        df = pd.read_csv(infile)

        # detect valid event columns present in file
        valid_cols = [c for c in Param.EVENT_COLS if c in df.columns]

        for idx, row in df.iterrows():

            is_up = ("Dir" in df.columns and str(row.get("Dir", "")).upper() == "UP")

            cols = valid_cols[::-1] if is_up else valid_cols

            for c in cols:

                v = row[c]

                if pd.notna(v):

                    event_id = int(v)

                    if event_id in times:

                        df.at[idx, c] = round(times[event_id], 2)

                    else:

                        df.at[idx, c] = np.nan

        # compute optimized start time safely
        df["Optimized_Start_Time"] = df[valid_cols].min(axis=1)

        # sort by optimized start
        df.sort_values("Optimized_Start_Time", inplace=True)

        df.to_csv(outfile, index=False)

        print("Saved", outfile)

    run(up, f"{Param.OUTPUT_DIR}/optimized_UP.csv")
    run(down, f"{Param.OUTPUT_DIR}/optimized_DOWN.csv")

# ============================================================
# PLOT
# ============================================================

def plot_time_vs_station_with_linkages(
        up_csv,
        down_csv,
        active_links,
        event_times,
        event_term_type):

    df_up = pd.read_csv(up_csv)
    df_dn = pd.read_csv(down_csv)

    # define station order locally (avoid global dependency)
    STATION_ORDER = ["CCG","MCT","DDR","BA","ADH","GMN","BVI","BYR","BSR","VR","DRD"]
    station_idx = {s:i for i,s in enumerate(STATION_ORDER)}

    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(18,6),sharey=True)

    # ----------------------------
    # draw services
    # ----------------------------
    def draw_services(ax, df, color):

        if df.empty:
            return

        for _,r in df.iterrows():

            xs=[]
            ys=[]

            for s in STATION_ORDER:

                for suf in ["a","d"]:

                    c=f"{s}{suf}"

                    if c in df.columns and pd.notna(r[c]):

                        xs.append(float(r[c]))
                        ys.append(station_idx[s])

            if len(xs)>1:
                ax.plot(xs,ys,color,alpha=.7)

    # FAST services
    draw_services(ax1,
                  df_up[df_up["Type"].str.lower()=="fast"],
                  "r-")

    draw_services(ax1,
                  df_dn[df_dn["Type"].str.lower()=="fast"],
                  "b-")

    ax1.set_title("FAST Track: UP + DOWN with Solver Linkages",
                  weight="bold")

    # SLOW services
    draw_services(ax2,
                  df_up[df_up["Type"].str.lower()=="slow"],
                  "r-")

    draw_services(ax2,
                  df_dn[df_dn["Type"].str.lower()=="slow"],
                  "b-")

    ax2.set_title("SLOW Track: UP + DOWN with Solver Linkages",
                  weight="bold")

    # ----------------------------
    # overlay rake link arrows
    # ----------------------------
    for dep, arr in active_links:

        if dep not in event_times or arr not in event_times:
            continue

        dep_info = event_term_type.get(dep)
        arr_info = event_term_type.get(arr)

        if dep_info is None or arr_info is None:
            continue

        dep_station, dep_type = dep_info
        arr_station, arr_type = arr_info

        # must be same terminal + same type
        if dep_station != arr_station or dep_type != arr_type:
            continue

        if dep_station not in station_idx:
            continue

        y = station_idx[dep_station]

        x_arr = event_times[arr]
        x_dep = event_times[dep]

        ax = ax1 if dep_type=="fast" else ax2

        ax.annotate(
            "",
            xy=(x_dep, y),
            xytext=(x_arr, y),
            arrowprops=dict(
                arrowstyle="->",
                linestyle="--",
                color="black",
                lw=1.5
            )
        )

    # ----------------------------
    # axis formatting
    # ----------------------------
    for ax in (ax1,ax2):

        ax.set_yticks(range(len(STATION_ORDER)))
        ax.set_yticklabels(STATION_ORDER)
        ax.set_xlabel("Time (minutes)")
        ax.grid(True,alpha=.3)

    ax1.set_ylabel("Station")

    plt.tight_layout()

    outfile = "time_vs_station_solver_linkages_fast_slow_split.png"

    plt.savefig(outfile,dpi=200)
    plt.close()

    print("Saved",outfile)
def build_event_terminal_type_map(loader):
    
    M={}

    T=loader.J["turnaround"]

    for st in T:

        for typ in ["fast","slow"]:

            for d in ["UP","DOWN"]:

                for a in T[st][typ][d]["arr"]:
                    if a is not None:
                        M[int(a)] = (st,typ)

                for dep in T[st][typ][d]["dep"]:
                    if dep is not None:
                        M[int(dep)] = (st,typ)

    return M
# ============================================================
# MAIN
# ============================================================

def main():
    
    print("Loading JSON...")

    loader = JsonLoader("milp_preprocessed2.json")

    print("Building MILP model...")

    scheduler = Scheduler(loader)

    model = scheduler.build()

    print("Building terminal-type map...")

    event_term_type = build_event_terminal_type_map(loader)

    print("Solving MILP using Gurobi...")

    opt = Optimizer(model)

    res = opt.solve()

    print("Solver status:", res.solver.termination_condition)

    if res.solver.termination_condition != pyo.TerminationCondition.optimal:

        print("Model not optimal. Exiting.")
        return

    print("Extracting solution...")

    times = opt.times()

    active_links = opt.rake_links()

    print("Saving rake linkages...")

    save_rake_links(active_links, model)

    print("Generating optimized CSV files...")

    replace_csv(
        times,
        up="1-o-event-ids_UP2.csv",
        down="1-o-event-ids_DOWN2.csv"
    )

    print("Generating plot...")

    plot_time_vs_station_with_linkages(

        f"{Param.OUTPUT_DIR}/optimized_UP.csv",

        f"{Param.OUTPUT_DIR}/optimized_DOWN.csv",

        active_links,

        times,

        event_term_type
    )

    rake_count = opt.rake_count()

    print("\n===================================")
    print("Minimum number of rakes required =", rake_count)
    print("Total rake linkages =", len(active_links))
    print("===================================")

    print("Done.")

if __name__=="__main__":
    main()
