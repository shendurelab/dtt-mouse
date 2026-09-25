#!/usr/bin/env python3
"""Timed-clade couplings under a trajectory-matched null (preliminary).

Clades as in c7: every lineage crossing time t, clade = its sampled descendants,
requiring >=3 tips of the blastomere. Statistic as in c7: obs = clades containing
both types; called when z=(obs-exp)/sd >= 3 and obs >= 5 in BOTH blastomeres.

Null is by explicit permutation (c7 uses a closed form that assumes exchangeable
tips, which a stratified null breaks):
   global : cell-type labels shuffled within blastomere
   traj   : shuffled within (blastomere, major trajectory)
"""
import numpy as np, csv, sys
from collections import defaultdict

import os as _os
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))

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


FP="/Users/jay.shendure/Dropbox/claude/current/final_push"
V8="/Users/jay.shendure/Dropbox/claude/mouse_sprint/tape_pipeline/figures/v8"
NPERM=int(sys.argv[1]) if len(sys.argv)>1 else 10
SLICES=[float(x) for x in (sys.argv[2:] or ["9.0","11.0","13.0"])]
MIN_CLADE, Z_THR, MIN_OBS = 3, 3.0, 5
rng=np.random.default_rng(0)

z=np.load(_support("mergedtree_dttpq_v8.npz", "DTT_MERGED_TREE_NPZ", 'Generate it with: python3 tools/make_merged_tree_npz.py (derived from support_data/merged_full_placed.nwk).'),
          allow_pickle=True)
par=z["parent"].astype(np.int64); tm=np.asarray(z["time"],float)
isleaf=np.asarray(z["is_leaf"],bool); names=[str(v) for v in z["names"]]
N=len(par); kids=[[] for _ in range(N)]
for v in range(1,N): kids[par[v]].append(v)

ct_of,tj_of={},{}
with open(f"{FP}/data/cell_metadata.v8.txt") as f:
    r=csv.reader(f,delimiter="\t"); h=next(r)
    a,b,c=h.index("cell_id"),h.index("celltype"),h.index("major_trajectory")
    for x in r: ct_of[x[a]]=x[b]; tj_of[x[a]]=x[c]

leaves=np.flatnonzero(isleaf)
# universe: >=100 cells in each blastomere subtree, as in c1/c7
side=np.full(N,-1,np.int8)
for s,c0 in enumerate(kids[0]):
    st=[c0]
    while st: u=st.pop(); side[u]=s; st.extend(kids[u])
cnt=defaultdict(lambda:[0,0])
for v in leaves:
    t=ct_of.get(names[v],"")
    if t and side[v]>=0: cnt[t][side[v]]+=1
types=sorted(t for t,(x,y) in cnt.items() if x>=100 and y>=100)
tix={t:i for i,t in enumerate(types)}; T=len(types)
ct=np.full(N,-1,np.int16)
for v in leaves: ct[v]=tix.get(ct_of.get(names[v],""),-1)
tjs=sorted({tj_of.get(names[v],"") for v in leaves})
tjx={t:i for i,t in enumerate(tjs)}
tjc=np.full(N,-1,np.int16)
for v in leaves: tjc[v]=tjx.get(tj_of.get(names[v],""),-1)
print(f"universe: {T} cell types")

def clade_of(t):
    """clade id per leaf at time t, per blastomere; -1 if clade too small"""
    cross=[v for v in range(1,N) if tm[par[v]]<t<=tm[v]]
    cid=np.full(N,-1,np.int64)
    for k,v in enumerate(cross):
        st=[v]
        while st:
            u=st.pop(); cid[u]=k; st.extend(kids[u])
    return cid,len(cross)

def cooc(lab,cid,ncl,mask):
    """clades x types presence -> co-occurrence counts (T x T)"""
    sel=mask&(lab>=0)&(cid>=0)
    P=np.zeros((ncl,T),np.float32)
    P[cid[sel],lab[sel]]=1.0
    keep=P.sum(1)>=0            # all clades; size filter applied via mask upstream
    return P.T@P

def run_slice(t):
    cid,ncl=clade_of(t)
    out={}
    for s in (0,1):
        m=(side==s)&isleaf
        # clade must hold >=MIN_CLADE tips of this blastomere
        sz=np.bincount(cid[m&(cid>=0)],minlength=ncl)
        ok=np.isin(cid,np.flatnonzero(sz>=MIN_CLADE))
        mask=m&ok
        obs=cooc(ct,cid,ncl,mask)
        nulls={}
        for kind in ("global","traj"):
            acc=np.zeros((ncl if False else T,T),np.float64); acc2=np.zeros_like(acc)
            by=defaultdict(list)
            for v in np.flatnonzero(mask):
                if ct[v]>=0:
                    by[(int(tjc[v]) if kind=="traj" else 0)].append(v)
            for _ in range(NPERM):
                p=ct.copy()
                for _,idx in by.items():
                    idx=np.array(idx); vv=p[idx].copy(); rng.shuffle(vv); p[idx]=vv
                c=cooc(p,cid,ncl,mask); acc+=c; acc2+=c*c
            mu=acc/NPERM; sd=np.sqrt(np.maximum(acc2/NPERM-mu**2,1e-9))
            nulls[kind]=(mu,sd)
        out[s]=(obs,nulls)
    return out,ncl

cp=[(r["celltype_1"],r["celltype_2"]) for r in
    csv.DictReader(open(f"{V8}/coupling_depth/figSXa_coupling_depth_AB.csv"))
    if r["coupling_depth_both_E"]]
cp=[(a,b) for a,b in cp if a in tix and b in tix]
tjof={}
for v in leaves:
    t=ct_of.get(names[v],"")
    if t in tix and t not in tjof: tjof[t]=tj_of.get(names[v],"")

print(f"{'slice':>7}{'clades':>9}   {'null':<8}{'called (all 88x88)':>20}{'of the 129 coupled':>20}"
      f"{'same-traj':>11}{'cross-traj':>11}")
for t in SLICES:
    res,ncl=run_slice(t)
    for kind in ("global","traj"):
        sig=np.ones((T,T),bool)
        for s in (0,1):
            obs,nulls=res[s]; mu,sd=nulls[kind]
            zz=(obs-mu)/sd
            sig &= (zz>=Z_THR)&(obs>=MIN_OBS)
        np.fill_diagonal(sig,False)
        tot=int(sig.sum()//2)
        hit=[(a,b) for a,b in cp if sig[tix[a],tix[b]]]
        ns=sum(1 for a,b in hit if tjof[a]==tjof[b]); nc=len(hit)-ns
        print(f"E{t:<6.2f}{ncl:>9,}   {kind:<8}{tot:>20,}{len(hit):>20}{ns:>11}{nc:>11}")
