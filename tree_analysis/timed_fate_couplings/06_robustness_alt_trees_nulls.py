#!/usr/bin/env python3
"""Dated coupling hierarchy under alternative trees, nulls and cell-type sets.

Clades and statistic exactly as c7: every lineage crossing t defines a clade of its
sampled descendants, requiring >=3 tips of the blastomere; a pair is called at a slice
when z=(obs-exp)/sd >= 3 and obs >= 5 in BOTH blastomeres; coupling depth = earliest
slice still called. The null is by explicit permutation (c7's closed form assumes
exchangeable tips, which a stratified null breaks).

env:  RUN_TREE  npz (default: placement tree)
      NULL_MODE global | traj | clone
      DROP_BLOOD 1 to exclude haematopoietic cell types
      NPERM, TAG
"""
import numpy as np, csv, os, sys, json, time
from collections import defaultdict

import os as _os
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))

def _results(*parts):
    """Intermediate/output dir for this analysis step, override with DTT_RESULTS."""
    base = _os.environ.get("DTT_RESULTS", _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "out"))
    p = _os.path.join(base, *parts) if parts else base
    _os.makedirs(_os.path.dirname(p) if _os.path.splitext(p)[1] else p, exist_ok=True)
    return p
import sys as _sys
_sys.path.insert(0, _os.path.join(_REPO, "tools"))
from tree_io import load_tree

def _support(name, env=None, hint=None):
    """Resolve an input under support_data/, overridable by env var."""
    if env:
        v = _os.environ.get(env)
        if v:
            return v
    p = _os.path.join(_REPO, "support_data", name)
    if not _os.path.exists(p) and hint:
        raise FileNotFoundError(f"{p} not found. {hint}")
    return p


FP = _results() + _os.sep
TREE=os.environ.get("RUN_TREE",_os.environ.get("DTT_TREE", _os.path.join(_REPO, "support_data", "merged_full_placed.nwk")))
NULL=os.environ.get("NULL_MODE","global"); NPERM=int(os.environ.get("NPERM","20"))
DROP_BLOOD=os.environ.get("DROP_BLOOD","")=="1"; TAG=os.environ.get("TAG","run")
MIN_CELLS,MIN_CLADE,Z_THR,MIN_OBS=100,3,3.0,5
TIMES=[round(7.0+0.25*i,2) for i in range(26)]
BLOOD={"White_blood_cells","Definitive_erythroid","Primitive_erythroid","Megakaryocytes",
       "Mast_cells","B_cells","T_cells"}
rng=np.random.default_rng(0); t0=time.time()

z=load_tree(TREE)
par=z["parent"].astype(np.int64); tm=np.asarray(z["time"],float)
isleaf=np.asarray(z["is_leaf"],bool); names=[str(v) for v in z["names"]]
N=len(par); kids=[[] for _ in range(N)]
for v in range(1,N): kids[par[v]].append(v)
side=np.full(N,-1,np.int8)
for s,c in enumerate(kids[0]):
    st=[c]
    while st: u=st.pop(); side[u]=s; st.extend(kids[u])
ct_of,tj_of={},{}
with open(_os.path.join(_REPO, "support_data", "cell_metadata.v8.txt.gz")) as f:
    r=csv.reader(f,delimiter="\t"); h=next(r)
    a,b,c=h.index("cell_id"),h.index("celltype"),h.index("major_trajectory")
    for x in r: ct_of[x[a]]=x[b]; tj_of[x[a]]=x[c]
leaves=np.flatnonzero(isleaf)
cnt=defaultdict(lambda:[0,0])
for v in leaves:
    t=ct_of.get(names[v],"")
    if t and side[v]>=0: cnt[t][side[v]]+=1
types=sorted(t for t,(x,y) in cnt.items() if x>=MIN_CELLS and y>=MIN_CELLS)
if DROP_BLOOD:
    types=[t for t in types if tj_of.get(next(names[v] for v in leaves if ct_of.get(names[v],"")==t),"") not in BLOOD]
tix={t:i for i,t in enumerate(types)}; T=len(types)
ct=np.full(N,-1,np.int16)
for v in leaves: ct[v]=tix.get(ct_of.get(names[v],""),-1)

if NULL=="traj":
    us=sorted({tj_of.get(names[v],"") for v in leaves}); um={t:i for i,t in enumerate(us)}
    STRAT=np.array([um.get(tj_of.get(names[v],""),-1) if isleaf[v] else -1 for v in range(N)],np.int64)
elif NULL=="clone":
    STRAT=np.full(N,-1,np.int64)
    CLONE_T=float(os.environ.get("CLONE_T","6.0"))
    cross=[v for v in range(1,N) if tm[par[v]]<CLONE_T<=tm[v]]
    for i,v in enumerate(cross):
        st=[v]
        while st: u=st.pop(); STRAT[u]=i; st.extend(kids[u])
else:
    STRAT=np.zeros(N,np.int64)
KEEP_N=int(os.environ.get("TIP_KEEP_N","0"))
if KEEP_N:
    _elig=np.array([v for v in leaves if ct[v]>=0])
    _drop=rng.choice(_elig, size=max(0,len(_elig)-KEEP_N), replace=False)
    ct[_drop]=-1
    print(f"subsampled tips: {len(_elig):,} -> {KEEP_N:,}", flush=True)
print(f"[{TAG}] tree={os.path.basename(TREE)} null={NULL} clone_T={os.environ.get('CLONE_T','-')} "
      f"types={T} nperm={NPERM}",flush=True)

def cooc(lab,cid,idx,ncl):
    P=np.zeros((ncl,T),np.float32); P[cid[idx],lab[idx]]=1.0
    return P.T@P

depth={}
for t in TIMES:
    cross=[v for v in range(1,N) if tm[par[v]]<t<=tm[v]]
    cid=np.full(N,-1,np.int64)
    for k,v in enumerate(cross):
        st=[v]
        while st: u=st.pop(); cid[u]=k; st.extend(kids[u])
    ncl=len(cross); sig=np.ones((T,T),bool)
    for s in (0,1):
        m=isleaf&(side==s)&(cid>=0)&(ct>=0)
        sz=np.bincount(cid[m],minlength=ncl)
        keep=np.flatnonzero(sz>=MIN_CLADE)
        idx=np.flatnonzero(m&np.isin(cid,keep))
        obs=cooc(ct,cid,idx,ncl)
        by=defaultdict(list)
        for v in idx: by[int(STRAT[v])].append(v)
        by={k:np.array(v) for k,v in by.items()}
        acc=np.zeros((T,T)); acc2=np.zeros((T,T))
        for _ in range(NPERM):
            p=ct.copy()
            for _,ix in by.items():
                vv=p[ix].copy(); rng.shuffle(vv); p[ix]=vv
            cc=cooc(p,cid,idx,ncl); acc+=cc; acc2+=cc*cc
        mu=acc/NPERM; sd=np.sqrt(np.maximum(acc2/NPERM-mu**2,1e-9))
        sig&=((obs-mu)/sd>=Z_THR)&(obs>=MIN_OBS)
    np.fill_diagonal(sig,False)
    for i,j in zip(*np.where(np.triu(sig,1))):
        depth.setdefault((types[i],types[j]),t)
    print(f"   E{t:<6.2f} clades {ncl:>8,}  called {int(np.triu(sig,1).sum()):>5}  "
          f"cumulative pairs {len(depth):>4}  [{time.time()-t0:.0f}s]",flush=True)

out=f"{FP}/fig6_v8/coupling_{TAG}.csv"
with open(out,"w",newline="") as f:
    w=csv.writer(f); w.writerow(["celltype_1","celltype_2","coupling_depth_E"])
    for (a,b),d in sorted(depth.items(),key=lambda kv:kv[1]): w.writerow([a,b,d])
cov=sorted({x for p in depth for x in p})
print(f"\n[{TAG}] coupled pairs {len(depth)}   cell types covered {len(cov)}/{T}   "
      f"earliest E{min(depth.values()) if depth else float('nan')}   -> {out}")
