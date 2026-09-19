#!/usr/bin/env python3
"""STEP 2 -- permutation test for recurrent progenitor -> post-mitotic sibling pairings,
swept over the evidence-threshold grid.

UNIT OF ANALYSIS = the post-mitotic CELL. Each post-mitotic cell is counted at most once:
among its qualifying non-post-mitotic tree-siblings it is assigned the NEAREST one
(fewest discordant DNA-typewriter sites; ties -> most shared edits -> deterministic pair
order). A cell in a terminal polytomy therefore contributes exactly one candidate division,
not C(k,2) of them.

QUALIFYING PAIR = tree-sibling pair (cherry or terminal polytomy; see h1_pairs.py) with
    MRCA age >= AMIN      ("the ancestor is dated late")
    discordant sites <= DMAX  ("few edits separate them")
Both criteria are label-independent, so they define a fixed eligible pair set that is
identical for observed and permuted labels.

NULL = tip-label permutation: cell-type labels are shuffled across ALL annotated tips of the
same blastomere with the tree held fixed (the null used in the preceding "tree siblings share
fate" subsection), and the ENTIRE procedure above -- heterotypic test, one-use-per-cell
nearest-partner assignment, counting -- is recomputed. Labels are always permuted WITHIN
blastomere, so each independently reconstructed half-embryo keeps its own cell-type
composition.

The PRIMARY test is over the whole tree (A+B counts pooled, null pooled the same way), since
that is the quantity of interest: is this post-mitotic type recurrently paired with this
progenitor above chance. Because A (B1) and B (B2) are separately reconstructed half-embryos,
the per-blastomere tests are also reported, and their agreement is used as validation rather
than as a gate -- gating on FDR<1% in both halves independently mostly rejects on power in the
smaller half B (e.g. otic epithelium -> otic sensory neurons is 28 obs vs 1.3 exp, q=4e-16 in
A but has only 3 obs in B), which discards real signal.

STATISTICS per (progenitor -> post-mitotic) pairing, for A, B and pooled:
    obs   observed # of captured divisions (= # of distinct post-mitotic cells)
    exp   mean of the permuted counts;  vperm = variance of the permuted counts
    fold  (obs + 1) / (exp + 1)            <- pseudocount, as in the subsection above
    p     upper-tail probability of obs under the permutation null, evaluated parametrically
          because an empirical p-value floors at 1/(NPERM+1) and could never clear BH FDR<1%
          across thousands of pairings. Poisson(exp+1) when the null is not overdispersed;
          negative binomial matched to (exp, vperm) when it is. The pseudocount in the mean
          keeps pairings resting on one or two cells from reaching significance.
    q     Benjamini-Hochberg across all pairings with obs >= MINOBS in that stratum.

CAPTURED DIVISION (primary) = pooled q < FDR and pooled fold >= FOLD_MIN.
  tier "both"  = additionally q < FDR and fold >= FOLD_MIN in each blastomere independently
                 (the strictest tier; what the earlier backbone analysis reported)
  tier "concordant" = additionally enriched in the same direction (fold > 1) in both halves

Outputs: out/h2b_traj_sweep.csv (grid summary), out/h2b_traj_pairings_<tag>.csv (per-pairing detail for
every grid cell), out/h2b_traj_raw.npz (obs/exp/var matrices).  Deterministic given NPERM/SEED.
"""
import numpy as np, os, sys, time, csv
from scipy.stats import poisson, nbinom, spearmanr, pearsonr

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NPERM = int(os.environ.get("NPERM", 500))
SEED  = int(os.environ.get("SEED", 0))
PSEUDO, MINOBS, FDR, FOLD_MIN = 1.0, 3, 0.01, 3.0
AMINS = [(0.0, "none"), (12.0, "E12"), (12.5, "E12.5"), (13.0, "E13")]
DMAXS = [2, 1, 0]
log = lambda *a: print(*a, file=sys.stderr, flush=True)
t0 = time.time()

d = np.load(f"{HERE}/out/h1_pairs.npz", allow_pickle=True)
I, J = d["i"].astype(np.int64), d["j"].astype(np.int64)
AGE, DISC, SHAR, KG = d["age"], d["disc"], d["shared"], d["kgroup"]
pair_blast, leaf_ct, leaf_blast = d["pair_blast"], d["leaf_ct"].astype(np.int64), d["leaf_blast"]
types = list(d["types"]); is_pm_t = d["type_is_pm"]
NT = len(types)
log(f"{len(I):,} tree-sibling pairs; {NT} cell types ({int(is_pm_t.sum())} post-mitotic); "
    f"NPERM={NPERM} seed={SEED}")

BLASTS = [("B1", "A"), ("B2", "B")]


def select_counts(ii, jj, dd, ss, pm_all, lab):
    """Count captured divisions per (progenitor type, post-mitotic type).

    One post-mitotic cell -> at most one pairing: its nearest qualifying sibling.
    Returns a flat (NT*NT) count vector indexed prog_type*NT + pm_type.
    """
    a, b = pm_all[ii], pm_all[jj]
    het = a ^ b                                     # exactly one member post-mitotic
    if not het.any():
        return np.zeros(NT * NT, np.int64)
    ai = a[het]
    pmleaf = np.where(ai, ii[het], jj[het])
    prleaf = np.where(ai, jj[het], ii[het])
    dh = dd[het].astype(np.int64); sh = ss[het].astype(np.int64)
    # rank key: group by post-mitotic cell, then fewest discordant, then most shared edits.
    # Remaining ties fall to original pair order via a stable sort -> deterministic.
    key = pmleaf * 8192 + dh * 128 + (66 - sh)
    o = np.argsort(key, kind="stable")
    pl = pmleaf[o]
    first = np.r_[True, pl[1:] != pl[:-1]]          # keep the best row per post-mitotic cell
    sel = o[first]
    return np.bincount(lab[prleaf[sel]] * NT + lab[pmleaf[sel]], minlength=NT * NT)


def bh(p):
    m = len(p)
    if m == 0: return p
    o = np.argsort(p); q = np.empty(m)
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


# map each cell type to its major trajectory (types are nested within trajectories)
import csv as _csv
_tj = {}
with open("/Users/jay.shendure/Dropbox/claude/current/final_push/data/cell_metadata.v8.txt") as _f:
    _r = _csv.reader(_f, delimiter="\t"); _h = next(_r)
    _a, _b = _h.index("celltype"), _h.index("major_trajectory")
    for _x in _r: _tj.setdefault(_x[_a], _x[_b])
_names = [str(t) for t in d["types"]]
_uniq = sorted({_tj.get(t, "") for t in _names})
TJ_OF_TYPE = np.array([_uniq.index(_tj.get(t, "")) for t in _names], np.int64)
log(f"trajectory-matched null: {len(_uniq)} trajectories over {len(_names)} cell types")

cells = [(am, al, dm) for am, al in AMINS for dm in DMAXS]
STRATA = ["B1", "B2", "pool"]
sub = {}          # (blast, cellkey) -> eligible pair arrays
pool_ix = {}      # blast -> tip indices in the permutation pool
for bl, blname in BLASTS:
    pmask = pair_blast == bl
    pool_ix[bl] = np.flatnonzero((leaf_blast == bl) & (leaf_ct >= 0))
    # TRAJECTORY-MATCHED NULL: split the pool by the major trajectory of each tip's
    # cell type, so permuted labels can never leave their trajectory.
    _sub = {}
    for _v in pool_ix[bl]:
        _sub.setdefault(int(TJ_OF_TYPE[leaf_ct[_v]]), []).append(_v)
    pool_ix[bl] = [np.array(v) for v in _sub.values()]
    log(f"\n=== blastomere {blname} ({bl}): {int(pmask.sum()):,} sibling pairs, "
        f"{sum(len(x) for x in pool_ix[bl]):,} annotated tips in "
        f"{len(pool_ix[bl])} trajectory strata ===")
    for am, al, dm in cells:
        m = pmask & (AGE >= am) & (DISC >= 0) & (DISC <= dm)
        sub[(bl, (al, dm))] = (I[m], J[m], DISC[m], SHAR[m], int(m.sum()))
        if am == 0.0 or dm == 1:
            log(f"   age>={al:<6} disc<={dm}: {int(m.sum()):,} eligible sibling pairs")

lab0 = leaf_ct.copy()
pm_all0 = is_pm_t[lab0]
obs = {}
for st in STRATA:
    for al, dm in [(a[1], b) for a in AMINS for b in DMAXS]:
        if st == "pool":
            obs[(st, (al, dm))] = obs[("B1", (al, dm))] + obs[("B2", (al, dm))]
        else:
            v = sub[(st, (al, dm))]
            obs[(st, (al, dm))] = select_counts(v[0], v[1], v[2], v[3], pm_all0, lab0)
log("\nobserved captured divisions (pooled A+B), by threshold cell:")
for al, dm in [(a[1], b) for a in AMINS for b in DMAXS]:
    v = obs[("pool", (al, dm))]
    log(f"   age>={al:<6} disc<={dm}: {int(v.sum()):,} post-mitotic cells in "
        f"{int((v>0).sum()):,} distinct pairings")

S = {(st, k): np.zeros(NT * NT) for st in STRATA for k in [(a[1], b) for a in AMINS for b in DMAXS]}
S2 = {k: np.zeros(NT * NT) for k in S}
rng = np.random.default_rng(SEED)
labp = lab0.copy()
for p in range(NPERM):
    for bl, _ in BLASTS:
        for _ix in pool_ix[bl]:
            labp[_ix] = rng.permutation(lab0[_ix])
    pm_all = is_pm_t[labp]
    for k in [(a[1], b) for a in AMINS for b in DMAXS]:
        cs = {}
        for bl, _ in BLASTS:
            v = sub[(bl, k)]
            cs[bl] = select_counts(v[0], v[1], v[2], v[3], pm_all, labp).astype(np.float64)
        cs["pool"] = cs["B1"] + cs["B2"]
        for st in STRATA:
            S[(st, k)] += cs[st]; S2[(st, k)] += cs[st] * cs[st]
    if (p + 1) % 25 == 0:
        log(f"   perm {p+1}/{NPERM} [{time.time()-t0:.0f}s]")

res = {}
for st in STRATA:
    for k in [(a[1], b) for a in AMINS for b in DMAXS]:
        e = S[(st, k)] / NPERM
        var = np.maximum(S2[(st, k)] / NPERM - e * e, 0.0)
        np_ = (sub[("B1", k)][4] + sub[("B2", k)][4]) if st == "pool" else sub[(st, k)][4]
        res[(st, k)] = dict(obs=obs[(st, k)], exp=e, var=var, npairs=np_)

np.savez_compressed(f"{HERE}/out/h2b_traj_raw.npz",
    **{f"{st}|{al}|{dm}|{w}": res[(st, (al, dm))][w]
       for st in STRATA for al, dm in [(a[1], b) for a in AMINS for b in DMAXS]
       for w in ("obs", "exp", "var")},
    types=np.array(types, object), NT=NT, nperm=NPERM)

# ---- per-pairing statistics + BH, per blastomere -------------------------------------------
def stats_table(st, key):
    r = res[(st, key)]
    obs, e, var = r["obs"], r["exp"], r["var"]
    idx = np.flatnonzero(obs >= MINOBS)
    if len(idx) == 0:
        return {}, 0.0
    o = obs[idx].astype(float); ee = e[idx] + PSEUDO; vv = var[idx] + PSEUDO
    disp = float(np.median(vv / ee))
    pv = np.where(vv > ee * 1.05,
                  nbinom.sf(o - 1, np.maximum(ee ** 2 / np.maximum(vv - ee, 1e-9), 1e-9),
                            np.clip(ee / np.maximum(vv, 1e-9), 1e-12, 1 - 1e-12)),
                  poisson.sf(o - 1, ee))
    q = bh(pv)
    out = {}
    for n, ii in enumerate(idx):
        out[(ii // NT, ii % NT)] = dict(obs=int(obs[ii]), exp=float(e[ii]),
            fold=(obs[ii] + PSEUDO) / (e[ii] + PSEUDO), p=float(pv[n]), q=float(q[n]))
    return out, disp


gl = {}
for fn, ok in [(f"/Users/jay.shendure/Dropbox/claude/penultimate_clade_k_analysis/germ_layer_map_validated.csv",
                lambda r: r["status"] in ("data-backed", "tree-resolved")),
               (f"/Users/jay.shendure/Dropbox/claude/penultimate_clade_k_analysis/germ_layer_map_v6.csv",
                lambda r: r["confidence"] == "clear")]:
    if not os.path.exists(fn): continue
    for r in csv.DictReader(open(fn)):
        if r["celltype"] not in gl and r["germ_layer"] and ok(r): gl[r["celltype"]] = r["germ_layer"]

rows = []
robust = {}
for am, al, dm in cells:
    tP, dispP = stats_table("pool", (al, dm))
    tA, dispA = stats_table("B1", (al, dm))
    tB, dispB = stats_table("B2", (al, dm))
    # PRIMARY: pooled FDR<1% and pooled fold>=FOLD_MIN
    core = sorted([k for k, v in tP.items() if v["q"] < FDR and v["fold"] >= FOLD_MIN],
                  key=lambda k: -tP[k]["fold"])
    # replication tiers
    conc = [k for k in core if tA.get(k, {}).get("fold", 0) > 1 and tB.get(k, {}).get("fold", 0) > 1]
    both = [k for k in core if k in tA and k in tB and tA[k]["q"] < FDR and tB[k]["q"] < FDR
            and tA[k]["fold"] >= FOLD_MIN and tB[k]["fold"] >= FOLD_MIN]
    common = set(tA) & set(tB)
    if common:
        fA = np.array([tA[k]["fold"] for k in common]); fB = np.array([tB[k]["fold"] for k in common])
        rho = float(spearmanr(fA, fB).correlation); rlog = float(pearsonr(np.log2(fA), np.log2(fB))[0])
    else:
        rho = rlog = float("nan")
    nsame = sum(1 for k in core if gl.get(types[k[0]]) and gl.get(types[k[0]]) == gl.get(types[k[1]]))
    rows.append(dict(age_min=al, disc_max=dm,
        eligible_pairs=res[("pool", (al, dm))]["npairs"],
        cells_paired=int(res[("pool", (al, dm))]["obs"].sum()),
        tested_pool=len(tP), captured=len(core), captured_concordant=len(conc),
        captured_both_blastomeres=len(both), within_germlayer=nsame,
        dispersion_pool=round(dispP, 2), rho_AB=round(rho, 3), r_log2_AB=round(rlog, 3)))
    robust[(al, dm)] = core
    tag = f"age{al}_disc{dm}"
    with open(f"{HERE}/out/h2b_traj_pairings_{tag}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["progenitor", "postmitotic", "prog_germlayer", "pm_germlayer", "same_germlayer",
                    "obs", "exp", "fold", "q", "obsA", "foldA", "qA", "obsB", "foldB", "qB",
                    "captured", "concordant_AB", "sig_both_blastomeres"])
        for k in sorted(tP, key=lambda k: -tP[k]["fold"]):
            gi, pi = k; g1, g2 = gl.get(types[gi], ""), gl.get(types[pi], "")
            a, b = tA.get(k), tB.get(k)
            w.writerow([types[gi], types[pi], g1, g2, int(bool(g1) and g1 == g2),
                tP[k]["obs"], round(tP[k]["exp"], 2), round(tP[k]["fold"], 2), f"{tP[k]['q']:.3e}",
                a["obs"] if a else "", round(a["fold"], 2) if a else "", f"{a['q']:.3e}" if a else "",
                b["obs"] if b else "", round(b["fold"], 2) if b else "", f"{b['q']:.3e}" if b else "",
                int(k in set(core)), int(k in set(conc)), int(k in set(both))])

with open(f"{HERE}/out/h2b_traj_sweep.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader()
    for r in rows: w.writerow(r)

print("\n=== THRESHOLD SWEEP ===")
print("captured = pooled permutation FDR<1% and fold>=3  |  +conc = also enriched in both halves")
print("  +both = also FDR<1% & fold>=3 in each half independently (strictest)\n")
print(f"{'age':>7} {'disc':>5} {'elig pairs':>11} {'cells':>8} {'tested':>7} "
      f"{'CAPTURED':>9} {'+conc':>6} {'+both':>6} {'sameGL':>7} {'rho A/B':>8} {'disp':>6}")
for r in rows:
    print(f"{r['age_min']:>7} {r['disc_max']:>5} {r['eligible_pairs']:>11,} {r['cells_paired']:>8,} "
          f"{r['tested_pool']:>7} {r['captured']:>9} {r['captured_concordant']:>6} "
          f"{r['captured_both_blastomeres']:>6} {r['within_germlayer']:>7} {r['rho_AB']:>8} "
          f"{r['dispersion_pool']:>6}")

allc = [set(v) for v in robust.values()]
inter = set.intersection(*allc) if allc else set()
union = set.union(*allc) if allc else set()
print(f"\nthreshold-robust: {len(inter)} pairings captured in ALL {len(cells)} grid cells; "
      f"{len(union)} in at least one")
cnt = {}
for v in robust.values():
    for k in v: cnt[k] = cnt.get(k, 0) + 1
print(f"\npairings by threshold robustness:")
for k, c in sorted(cnt.items(), key=lambda x: (-x[1], -robust[(AMINS[0][1], DMAXS[0])].count(x[0]))):
    g1, g2 = gl.get(types[k[0]], "?"), gl.get(types[k[1]], "?")
    flag = "same" if (g1 == g2 and g1 != "?") else f"{g1[:4]}>{g2[:4]}"
    print(f"   {c:>2}/{len(cells)}  [{flag:>9}]  {types[k[0]]} -> {types[k[1]]}")
print(f"\nwrote out/h2b_traj_sweep.csv, out/h2b_traj_pairings_*.csv, out/h2b_traj_raw.npz [{time.time()-t0:.0f}s]")
