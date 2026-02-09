import json, os
from dataclasses import dataclass
import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class Param:
    BIG_M = 2000
    MIN_HEADWAY = 3.0

    TURN_LB = 5.0
    TURN_UB = 10.0

    RAKE_WEIGHT = 10000.0
    LINK_REWARD = 2000.0
    DELTA_WEIGHT = 1.0

    TIME_LIMIT = 600
    OUTPUT_DIR = "outputs_json"

    EVENT_COLS = [
        'CCGa','CCGd','MCTa','MCTd','DDRa','DDRd','BAa','BAd','ADHa','ADHd',
        'GMNa','GMNd','BVIa','BVId','BYRa','BYRd','BSRa','BSRd','VRa','VRd',
        'DRDa','DRDd'
    ]

@dataclass
class Service:
    start: float
    typ: str
    direction: str
    frm: str
    to: str
    segments: dict
    first_dep: int
    last_arr: int

@dataclass
class Station:
    code: str
    index: int
    arrivals: set = None
    departures: set = None

    def __post_init__(self):
        self.arrivals = self.arrivals or set()
        self.departures = self.departures or set()

    def add_arrival(self, e):
        self.arrivals.add(e)

    def add_departure(self, e):
        self.departures.add(e)


STATION_ORDER = ["CCG","MCT","DDR","BA","ADH","GMN","BVI","BYR","BSR","VR","DRD"]

class JsonLoader:

    def __init__(self,fname):
        with open(fname) as f:
            self.J=json.load(f)

    def services(self):
        out=[]
        for v in self.J["services"].values():
            out.append(Service(
                start=v["start_time"],
                typ=v["Type"],
                direction=v["Dir"],
                frm=v["From"],
                to=v["To"],
                segments=v["segments"],
                first_dep=int(v["first_dep"]),
                last_arr=int(v["last_arr"])
            ))
        return out

    def events(self):
        ev=set()
        for s in self.J["services"].values():
            ev.add(int(s["first_dep"]))
            ev.add(int(s["last_arr"]))
            for k in s["segments"]:
                a,b=k.split("_")
                ev.add(int(a)); ev.add(int(b))
        return sorted(ev)

    # Build station objects
    def stations(self):

        stations = {
            s: Station(s, i)
            for i, s in enumerate(STATION_ORDER)
        }

        for v in self.J["services"].values():

            if v["first_dep"] is not None:
                stations[v["From"]].add_departure(int(v["first_dep"]))

            if v["last_arr"] is not None:
                stations[v["To"]].add_arrival(int(v["last_arr"]))

        return stations

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
                    for d in down_dep:
                        P.append((int(d),int(a)))

                for a in down_arr:
                    for d in up_dep:
                        P.append((int(d),int(a)))

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

class Scheduler:

    def __init__(self,loader):
        self.services=loader.services()
        self.events=loader.events()
        self.turn=loader.turnaround_pairs()
        self.head=loader.headway_pairs()
        self.dist=loader.distribution()
        self.stations = loader.stations()   # ⭐ NEW
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
            if s.last_arr is not None:
                ARR.add(s.last_arr)
            DEP.add(s.first_dep)

        m.source=pyo.Var(DEP,within=pyo.Binary)
        m.sink=pyo.Var(ARR,within=pyo.Binary)

        m.trav=pyo.ConstraintList()
        m.headway=pyo.ConstraintList()
        m.rake=pyo.ConstraintList()
        m.distc=pyo.ConstraintList()
        m.start=pyo.ConstraintList()

        # traversal
        for s in self.services:
            for k,tt in s.segments.items():
                a,b=k.split("_")
                a=int(a); b=int(b)
                m.trav.add(m.t[b]-m.t[a]==tt)
                m.trav.add(m.t[b]>=m.t[a])

        # headway
        for i,j in m.HW:
            m.headway.add(m.t[j]-m.t[i]+(Param.BIG_M+Param.MIN_HEADWAY)*m.p[i,j]>=Param.MIN_HEADWAY)
            m.headway.add(m.t[j]-m.t[i]+(Param.BIG_M+Param.MIN_HEADWAY)*m.p[i,j]<=Param.BIG_M)

        # rake conservation
        for j in ARR:
            m.rake.add(sum(m.X[i,j] for i,jj in self.turn if jj==j) + m.sink[j] == 1)

        for i in DEP:
            m.rake.add(sum(m.X[ii,j] for ii,j in self.turn if ii==i) + m.source[i] == 1)

        m.rake.add(sum(m.source[i] for i in DEP)==sum(m.sink[j] for j in ARR))

        # turnaround
        for a,d in self.turn:
            m.rake.add(m.t[d]-m.t[a]>=Param.TURN_LB-Param.BIG_M*(1-m.X[a,d]))
            m.rake.add(m.t[d]-m.t[a]<=Param.TURN_UB+Param.BIG_M*(1-m.X[a,d]))

        # distribution
        for ids in self.dist.values():
            if len(ids)>1:
                ideal=60.0/(len(ids))
                for u in range(len(ids)-1):
                    a=int(ids[u]); b=int(ids[u+1])
                    m.distc.add(m.t[b]-m.t[a]>=ideal-m.delta)
                    m.distc.add(m.t[b]-m.t[a]<=ideal+m.delta)

        # start windows
        first=True
        for s in sorted(self.services,key=lambda x:x.start):
            if first:
                m.start.add(m.t[s.first_dep]==480)
                first=False
            else:
                m.start.add(m.t[s.first_dep]>=s.start)
                m.start.add(m.t[s.first_dep]<=s.start+19)

        # objective
        m.obj=pyo.Objective(
            expr=Param.RAKE_WEIGHT*(sum(m.source[i] for i in DEP)+sum(m.sink[j] for j in ARR))
            + Param.DELTA_WEIGHT*m.delta,
            sense=pyo.minimize)

        return m

class Optimizer:
    def __init__(self,m): self.m=m

    def solve(self):
        s=SolverFactory("cbc")
        s.options["seconds"]=Param.TIME_LIMIT
        return s.solve(self.m,tee=True)

    def times(self):
        return {i:float(self.m.t[i].value) for i in self.m.EVENTS if self.m.t[i].value is not None}


def two_csv(times,up="1-o-event-ids_UP.csv",down="1-o-event-ids_DOWN.csv"):
    os.makedirs(Param.OUTPUT_DIR,exist_ok=True)

    def run(f,out):
        df=pd.read_csv(f)

        for idx,row in df.iterrows():

            is_up = ("Dir" in row and str(row["Dir"]).upper()=="UP")
            cols = Param.EVENT_COLS[::-1] if is_up else Param.EVENT_COLS

            for c in cols:
                if c in df.columns:
                    v=row[c]
                    if pd.notna(v):
                        df.at[idx,c]=round(times.get(int(v),np.nan),2)

        df["Optimized_Start_Time"]=df[Param.EVENT_COLS].min(axis=1)
        df.sort_values("Optimized_Start_Time").to_csv(out,index=False)
        print("Saved",out)

    run(up,f"{Param.OUTPUT_DIR}/optimized_UP.csv")
    run(down,f"{Param.OUTPUT_DIR}/optimized_DOWN.csv")


def plot_time_vs_station_with_linkages(up_csv, down_csv, active_links, event_times, stations):

    df_up=pd.read_csv(up_csv)
    df_dn=pd.read_csv(down_csv)

    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(18,6),sharey=True)

    def draw_services(ax,df,color):
        for _,r in df.iterrows():
            xs=[]; ys=[]
            for s in STATION_ORDER:
                for suf in ["a","d"]:
                    c=f"{s}{suf}"
                    if c in df.columns and pd.notna(r[c]):
                        xs.append(r[c])
                        ys.append(stations[s].index)
            ax.plot(xs,ys,color,alpha=.7)

    draw_services(ax1, df_up[df_up["Type"].str.lower()=="fast"], "r-")
    draw_services(ax1, df_dn[df_dn["Type"].str.lower()=="fast"], "b-")
    ax1.set_title("FAST Track")

    draw_services(ax2, df_up[df_up["Type"].str.lower()=="slow"], "r-")
    draw_services(ax2, df_dn[df_dn["Type"].str.lower()=="slow"], "b-")
    ax2.set_title("SLOW Track")

    # rake link arrows
    for dep,arr in active_links:

        if dep not in event_times or arr not in event_times:
            continue

        st=None
        for s in stations.values():
            if dep in s.departures or dep in s.arrivals:
                st=s
                break

        if st is None:
            continue

        y=st.index
        x1=event_times[arr]
        x2=event_times[dep]

        for ax in (ax1,ax2):
            ax.annotate("",xy=(x2,y),xytext=(x1,y),
                        arrowprops=dict(arrowstyle="->",linestyle="--",color="k"))

    for ax in (ax1,ax2):
        ax.set_yticks(range(len(STATION_ORDER)))
        ax.set_yticklabels(STATION_ORDER)
        ax.set_xlabel("Time (minutes)")
        ax.grid(True,alpha=.3)

    ax1.set_ylabel("Station")
    plt.tight_layout()
    plt.savefig("time_vs_station_solver_linkages.png",dpi=200)
    plt.close()
    print("Saved time_vs_station_solver_linkages.png")

def main():

    os.makedirs(Param.OUTPUT_DIR,exist_ok=True)

    loader=JsonLoader("milp_preprocessed.json")
    scheduler=Scheduler(loader)
    model=scheduler.build()
    stations=scheduler.stations

    res=Optimizer(model).solve()

    m=model

    active_links=[(i,j) for (i,j) in m.LINKS if m.X[i,j].value>0.5]

    rows=[{
        "dep_event":dep,
        "arr_event":arr,
        "turnaround_time":m.t[dep].value-m.t[arr].value
    } for dep,arr in active_links]

    pd.DataFrame(rows).to_csv("rake_linkages.csv",index=False)

    if res.solver.termination_condition!=pyo.TerminationCondition.optimal:
        return

    times=Optimizer(model).times()

    two_csv(times)

    plot_time_vs_station_with_linkages(
        f"{Param.OUTPUT_DIR}/optimized_UP.csv",
        f"{Param.OUTPUT_DIR}/optimized_DOWN.csv",
        active_links,
        times,
        stations
    )

    print("Solved. Events:",len(times))


if __name__=="__main__":
    main()
