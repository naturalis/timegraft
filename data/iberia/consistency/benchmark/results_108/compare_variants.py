"""The table of results_108/README.md (#108, step 1): the variants of ages.timetree_time compared on the backbone
metrics of the gate, from builds in /tmp/g108/VARIANT/{clean,noise_s1..s5} (paths of the machine it ran on).

    PYTHONPATH=benchmark:src python compare_variants.py [VARIANT ...] > variants.tsv
"""
import json, csv, statistics, inspect, collections, sys
import benchlib
cls=[c for n,c in inspect.getmembers(benchlib, inspect.isclass) if hasattr(c,"stem_age")][0]
def factor(a,b): return max(a,b)/min(a,b)
def q(x,p): x=sorted(x); return x[min(len(x)-1,int(p*len(x)))]
V=sys.argv[1:] or ["raw","downweight","replace","adjusted"]
out=[]
ref=cls("/tmp/g108/raw/clean/backbone.tre")
for v in V:
    D=f"/tmp/g108/{v}/clean"
    m=json.load(open(D+"/metadata.json"))["ages"]["final"]
    bb=cls(D+"/backbone.tre"); moves=[]; stems=[]
    for s in range(1,6):
        nb=cls(f"/tmp/g108/{v}/noise_s{s}/backbone.tre"); a,r=nb.node_ages(),bb.node_ages()
        moves+=[factor(a[k],r[k]) for k in a if k in r and r[k]>0 and a[k]>0]
        stems+=[factor(nb.stem_age(g),bb.stem_age(g)) for g in nb.leaf if g in bb.leaf]
    e=[factor(float(r["age"]),float(r["timetree_age"])) for r in csv.DictReader(open(f"/tmp/g108/{v}/holdout.tsv"),delimiter="\t") if r["method"]=="fit"]
    a,ra=bb.node_ages(),ref.node_ages()
    f=[factor(a[k],ra[k]) for k in a if k in ra and ra[k]>0]
    st=[factor(bb.stem_age(g),ref.stem_age(g)) for g in bb.leaf]
    out.append(dict(variant=v, calibrations=m["n"], severe_2a=m["severe"], rms_2b=m["weighted_rms"], negative_2c=bb.negative,
        holdout_median_4a=round(statistics.median(e),3), holdout_p90_4b=round(q(e,0.9),3),
        noise_p90_5a=round(q(moves,0.9),3), noise_over_1_5_5b=round(sum(x>1.5 for x in moves)/len(moves),3), noise_stem_p90=round(q(stems,0.9),3),
        nodes_moved_10pc=sum(x>1.1 for x in f), nodes_moved_1_5=sum(x>1.5 for x in f), stems_moved_10pc=sum(x>1.1 for x in st),
        shrimp_stem=round(bb.stem_age("Eusergestes"),1)))
w=csv.DictWriter(sys.stdout,fieldnames=list(out[0]),delimiter="\t"); w.writeheader(); [w.writerow(r) for r in out]
