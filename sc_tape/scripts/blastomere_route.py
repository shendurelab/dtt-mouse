#!/usr/bin/env python3
"""
v4 blastomere routing stage  --  split the single-cell TAPE consensus matrix into
two per-blastomere (B1 / B2) tree-input matrices, cleaning the first (B1-vs-B2)
cell division by ROUTING ON THE FOUNDER SIGNATURE instead of letting NJ distance
decide it.

WHY
---
On the merged E13.5 NJ tree, ~4% of cells are misassigned at the very first
(blastomere) split: low-coverage cells that either lost their blastomere-defining
tapes (pulled into the low-signal clade by shared missingness) or picked up a
low-support ambient "bleed" edit at a founder locus. They are not doublets and
not homoplasy; the B1/B2 founder edits are robust and non-circular. Because the
two blastomere trees are built separately downstream, we can fix this BEFORE NJ
by routing each cell to its blastomere from the founder edits directly.

WHAT THIS STAGE DOES  (calls/filters otherwise UNCHANGED)
---------------------------------------------------------
  1. Derive the B1/B2 DEFINING SITES: site-1 founder alleles that are near-fixed
     on one side of the first division and enriched vs the other. Derived from
     CONFIDENT (well-covered) cells via a rarity-weighted PCA bipartition, so
     shared missingness never drives the partition. (Non-circular: re-deriving
     from the clean clades is identical.)
  2. ROUTE every cell by a majority vote over the defining sites -> B1 or B2.
  3. MASK only the cross-blastomere defining calls: any founder locus whose call
     carries the OPPOSITE blastomere's allele is set to missing (provably bleed).
     This cleans the vote and stops the stray edit distorting within-tree placement.
     No other (within-set) masking is done.
  4. SET ASIDE only the genuinely ambiguous remainder (no founder recovered, a
     tie, or bleed across >=2 tapes). Everything confidently-routing is KEPT --
     including low-coverage-but-clean cells -- so sampling is not biased.

OUTPUTS  (into --outdir)
------------------------
  <out>.B1_tape_consensus.tsv   B1 cells, tree-input format for parse_tape_consensus.R
  <out>.B2_tape_consensus.tsv   B2 cells,  "
  <out>.set_aside.tsv           the ambiguous cells + reason (handed to colleagues)
  <out>.routing_labels.tsv      per-cell: blastomere, votes, masked tape, reason
  <out>.defining_sites.tsv      the derived B1/B2 founder alleles (auditable)
  <out>.routing_report.txt      counts + parameters

Downstream is UNCHANGED: run the committed R distance/NJ/dating pipeline once per
blastomere, e.g.
  DTT_TSV=<out>.B1_tape_consensus.tsv DTT_OUT=.../B1.dist.gz Rscript run_full_distance.R
  DTT_TSV=<out>.B2_tape_consensus.tsv DTT_OUT=.../B2.dist.gz Rscript run_full_distance.R

Usage:
  python3 blastomere_route.py <consensus.tsv[.gz]> [--outdir DIR] [--out PREFIX]
Input is the v3 single-cell consensus matrix (cell + 11 integration dash-chain
columns + n_loci/n_doublet_loci/mean_dominance), e.g. tables/v3/e3.cell_consensus.tsv.gz.
"""
from __future__ import annotations
import argparse, csv, gzip, os, sys
from collections import Counter
import numpy as np

# --- parameters (single source of truth for this stage) --------------------
QC_COLS       = {"n_loci", "n_doublet_loci", "mean_dominance"}
N_SITES       = 6
FOUNDER_SITE  = 0          # site index used as the founder locus (site 1)
FIX_THRESH    = 0.85       # a defining allele is near-fixed (>= this) on its side
DISC_MIN      = 0.30       # ...and enriched by >= this vs the other side
CONF_MIN_LOCI = 8          # "confident" cell (used only to DERIVE the split axis)
MIN_REC       = 20         # min recovered cells per side to call a locus discriminating
MASK_TOKEN    = "NA"       # tree-input missing token (parse_tape_consensus.R -> all-missing tape)


# --- IO --------------------------------------------------------------------
def load_consensus(path):
    """Read the v3 consensus TSV(.gz). Returns (cells, barcodes, tokens, nloci, qc).
    tokens[i, g] = 6-tuple of site strings ('U' unedited) or None (integration not
    observed). qc[i] = (n_loci, n_doublet_loci, mean_dominance) strings (passthrough)."""
    op = gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)
    rows = list(csv.reader(op, delimiter="\t"))
    op.close()
    hdr = rows[0]
    bcs = [c for c in hdr[1:] if c not in QC_COLS]
    bidx = [hdr.index(b) for b in bcs]
    qc_idx = {c: hdr.index(c) for c in QC_COLS if c in hdr}
    n, G = len(rows) - 1, len(bcs)
    cells = [r[0] for r in rows[1:]]
    tok = np.empty((n, G), dtype=object)
    nloci = np.zeros(n, int)
    qc = np.empty((n, 3), dtype=object)
    for i, r in enumerate(rows[1:]):
        for j, bi in enumerate(bidx):
            v = r[bi]
            if v == "":
                tok[i, j] = None
            elif v == "-":
                tok[i, j] = ("U",) * N_SITES
            else:
                e = v.split("-")
                tok[i, j] = tuple(e) + ("U",) * (N_SITES - len(e))
        nloci[i] = int(r[qc_idx["n_loci"]]) if "n_loci" in qc_idx else \
                   sum(1 for j in range(G) if tok[i, j] is not None)
        qc[i] = [r[qc_idx[c]] if c in qc_idx else "" for c in
                 ("n_loci", "n_doublet_loci", "mean_dominance")]
    return cells, bcs, tok, nloci, qc


def encode_founder(tok):
    """Encode the founder-site allele per (cell, integration): 0 unedited, -1 missing,
    >=1 allele id. Returns (code, vocab) where vocab[g] = {allele_str: id}."""
    n, G = tok.shape
    vocab = [{} for _ in range(G)]
    code = np.full((n, G), -1, np.int16)
    for g in range(G):
        vv = vocab[g]
        for i in range(n):
            t = tok[i, g]
            if t is None:
                continue
            a = t[FOUNDER_SITE]
            if a == "U":
                code[i, g] = 0
            else:
                if a not in vv:
                    vv[a] = len(vv) + 1
                code[i, g] = vv[a]
    return code, vocab


# --- derive defining sites -------------------------------------------------
def rarity_weights(code, A):
    """-log(global recovered-allele fraction) per (integration, allele). Rare founder
    alleles carry more evidence."""
    G = code.shape[1]
    W = np.zeros((G, A))
    for g in range(G):
        rec = code[:, g][code[:, g] >= 0]
        if len(rec) == 0:
            continue
        frac = np.bincount(rec, minlength=A) / len(rec)
        W[g] = np.where(frac > 0, -np.log(np.where(frac > 0, frac, 1.0)), 0.0)
    return W


def confident_axis(code, W, conf, A):
    """Leading PCA axis of confident cells over sqrt(rarity) one-hot founder edits
    (row L2-normalized), then an Otsu (balance-weighted) threshold -> boolean split.
    Deterministic (no randomness). Returns the tentative split mask over `conf`."""
    G = code.shape[1]
    p = G * A
    F = np.zeros((len(conf), p), np.float32)
    for ri, i in enumerate(conf):
        for g in range(G):
            e = code[i, g]
            if e >= 1 and W[g, e] > 0:
                F[ri, g * A + e] = np.sqrt(W[g, e])
        nrm = np.linalg.norm(F[ri])
        if nrm > 0:
            F[ri] /= nrm
    F = F[:, (F != 0).any(0)]
    mu = F.mean(0)
    Fc = F - mu
    _, V = np.linalg.eigh(Fc.T @ Fc)
    proj = Fc @ V[:, -1]
    if proj.mean() < 0:
        proj = -proj
    # Otsu / balance-weighted between-class variance threshold
    x = np.sort(proj); m = len(x); cs = np.cumsum(x); tot = cs[-1]
    cand = np.unique(np.clip((np.linspace(0.02, 0.98, 385) * m).astype(int), 1, m - 1))
    best, bt = -np.inf, cand[len(cand) // 2]
    for t in cand:
        f = t / m
        sc = f * (1 - f) * ((tot - cs[t - 1]) / (m - t) - cs[t - 1] / t) ** 2
        if sc > best:
            best, bt = sc, t
    thr = 0.5 * (x[bt - 1] + x[bt])
    return proj > thr


def derive_defining(code, vocab, nloci, A):
    """Derive per-side defining founder alleles. Returns (defs, sideL_barcodes) where
    defs = list of (g, allele_id, side in {'L','R'}, fL, fR). L/R are provisional; the
    caller orients them to B1/B2."""
    W = rarity_weights(code, A)
    conf = np.where(nloci >= CONF_MIN_LOCI)[0]
    if len(conf) < 2 * MIN_REC:
        raise SystemExit(f"too few confident cells ({len(conf)}) to derive the split")
    right = confident_axis(code, W, conf, A)
    L0, R0 = conf[~right], conf[right]
    defs = []
    for g in range(code.shape[1]):
        cl = code[L0, g]; cr = code[R0, g]
        recL = cl[cl >= 0]; recR = cr[cr >= 0]
        if len(recL) < MIN_REC or len(recR) < MIN_REC:
            continue
        fL = np.bincount(recL, minlength=A) / len(recL)
        fR = np.bincount(recR, minlength=A) / len(recR)
        for a in range(1, A):
            if fL[a] >= FIX_THRESH and (fL[a] - fR[a]) >= DISC_MIN:
                defs.append((g, a, "L", float(fL[a]), float(fR[a])))
            if fR[a] >= FIX_THRESH and (fR[a] - fL[a]) >= DISC_MIN:
                defs.append((g, a, "R", float(fL[a]), float(fR[a])))
    return defs, (L0, R0)


# --- route + mask + set aside ----------------------------------------------
def route(code, defs):
    """Vote each cell over the defining alleles and apply the route/mask/set-aside
    rule. Returns a dict of per-cell arrays.

    Rule (per cell):
      w = winning-side votes, l = losing-side votes.
      ROUTE to argmax side iff  w >= 1 and w > l and l <= 1.
        - if l == 1: MASK the single opposite-side founder integration (set NA).
      SET ASIDE iff  w == 0 (no founder)  or  w == l (tie)  or  l >= 2 (multi-bleed).
    """
    n, G = code.shape
    Ldef = [(g, a) for g, a, s, *_ in defs if s == "L"]
    Rdef = [(g, a) for g, a, s, *_ in defs if s == "R"]
    vL = np.zeros(n, int); vR = np.zeros(n, int)
    # also remember, per cell, which integration carries an opposite-side allele
    L_int = {g: a for g, a in Ldef}       # (only one defining allele per (g,side) in practice)
    R_int = {g: a for g, a in Rdef}
    for g, a in Ldef:
        vL += (code[:, g] == a)
    for g, a in Rdef:
        vR += (code[:, g] == a)
    w = np.maximum(vL, vR); l = np.minimum(vL, vR)
    winnerL = vL > vR
    routed = (w >= 1) & (w > l) & (l <= 1)
    none_ = (w == 0)
    tie = (w == l) & (w > 0)
    multibleed = (l >= 2) & (w > l)
    aside = none_ | tie | multibleed
    # blastomere label 'L'/'R' for routed cells; masked integration index (or -1)
    mask_g = np.full(n, -1, int)
    for i in np.where(routed & (l == 1))[0]:
        # the losing side's defining integration this cell carries
        losing = R_int if winnerL[i] else L_int
        for g, a in losing.items():
            if code[i, g] == a:
                mask_g[i] = g
                break
    reason = np.array(["routed"] * n, dtype=object)
    reason[none_] = "set_aside:no_founder"
    reason[tie] = "set_aside:tie"
    reason[multibleed] = "set_aside:multibleed"
    return dict(vL=vL, vR=vR, routed=routed, winnerL=winnerL, aside=aside,
                none_=none_, tie=tie, multibleed=multibleed, mask_g=mask_g, reason=reason)


def orient(routed, winnerL):
    """Orient provisional L/R to B1/B2. Convention: B1 = the LARGER routed group
    (matches the tree pipeline's blast-1 = larger root subclade). Returns a function
    cell_index -> 'B1'/'B2' for routed cells."""
    nL = int((routed & winnerL).sum()); nR = int((routed & ~winnerL).sum())
    L_is_B1 = nL >= nR
    return L_is_B1


# --- emit ------------------------------------------------------------------
def tape_token(t):
    return MASK_TOKEN if t is None else "|".join(t)


def cell_pass_qc(qc_row, min_n_loci, max_doublet_loci, min_mean_dom):
    """pass_qc for a cell from its passthrough QC (n_loci, n_doublet_loci,
    mean_dominance). A cell passes iff it clears ALL active thresholds. Defaults
    (min_n_loci=0, max_doublet_loci=inf, min_mean_dom=0) => every cell passes,
    matching the pre-filter behaviour."""
    nl = int(qc_row[0]); nd = int(qc_row[1]); md = float(qc_row[2])
    return int(nl >= min_n_loci and nd <= max_doublet_loci and md >= min_mean_dom)


def write_matrix(path, cells, bcs, tok, qc, keep_idx, mask_g,
                 min_n_loci=0, max_doublet_loci=10**9, min_mean_dom=0.0):
    """Write a per-blastomere tree-input TSV: cell_id + 11 |-token tapes + QC + pass_qc.
    For a kept cell, the masked founder integration (mask_g) is written NA (bleed).
    pass_qc is set from the min_tapes / doublet / dominance thresholds (Tweak B); the
    cell stays in the matrix either way (parse_tape_consensus.R selects pass_qc==1).
    Returns (n_written, n_pass)."""
    n_pass = 0
    with open(path, "w") as f:
        f.write("cell_id\t" + "\t".join(bcs) +
                "\tn_loci\tn_doublet_loci\tmean_dominance\tpass_qc\n")
        for i in keep_idx:
            toks = []
            for g in range(len(bcs)):
                toks.append(MASK_TOKEN if g == mask_g[i] else tape_token(tok[i, g]))
            pq = cell_pass_qc(qc[i], min_n_loci, max_doublet_loci, min_mean_dom)
            n_pass += pq
            f.write(cells[i] + "\t" + "\t".join(toks) + "\t" +
                    "\t".join(map(str, qc[i])) + f"\t{pq}\n")
    return len(keep_idx), n_pass


def main():
    ap = argparse.ArgumentParser(description="v4 blastomere routing stage")
    ap.add_argument("consensus", help="v3 cell_consensus TSV(.gz)")
    ap.add_argument("--outdir", default=".", help="output directory")
    ap.add_argument("--out", default="e3", help="output prefix (default e3)")
    # Tweak B: pre-ship cell QC filters -> pass_qc (cells stay in the matrix; the
    # downstream distance step selects pass_qc==1). Defaults = off (all pass).
    ap.add_argument("--min-n-loci", type=int, default=0,
                    help="min tapes recovered (n_loci) to pass_qc [default 0 = off]")
    ap.add_argument("--max-doublet-loci", type=int, default=10**9,
                    help="max n_doublet_loci to pass_qc [default off]")
    ap.add_argument("--min-mean-dominance", type=float, default=0.0,
                    help="min mean_dominance to pass_qc [default 0.0 = off]")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    OP = lambda s: os.path.join(args.outdir, f"{args.out}.{s}")

    cells, bcs, tok, nloci, qc = load_consensus(args.consensus)
    n, G = tok.shape
    print(f"[route] loaded {n:,} cells x {G} integrations", flush=True)

    code, vocab = encode_founder(tok)
    A = max((len(v) for v in vocab), default=1) + 1
    defs, _ = derive_defining(code, vocab, nloci, A)
    R = route(code, defs)
    L_is_B1 = orient(R["routed"], R["winnerL"])

    # per-cell blastomere: B1/B2 or ""(set aside)
    is_b1 = R["routed"] & (R["winnerL"] == L_is_B1)
    is_b2 = R["routed"] & (R["winnerL"] != L_is_B1)
    b1_idx = np.where(is_b1)[0]; b2_idx = np.where(is_b2)[0]

    # ---- defining sites table (oriented) ----
    inv = [{v: k for k, v in vocab[g].items()} for g in range(G)]
    def side_label(s):  # provisional L/R -> B1/B2
        return "B1" if (s == "L") == L_is_B1 else "B2"
    with open(OP("defining_sites.tsv"), "w") as f:
        f.write("integration\tsite\tblastomere\tallele\tfrac_own_side\tfrac_other_side\n")
        for g, a, s, fL, fR in sorted(defs, key=lambda d: (d[0], d[2])):
            own, oth = (fL, fR) if s == "L" else (fR, fL)
            f.write(f"{bcs[g]}\tsite{FOUNDER_SITE+1}\t{side_label(s)}\t{inv[g][a]}\t"
                    f"{own:.3f}\t{oth:.3f}\n")

    # ---- per-cell routing labels ----
    with open(OP("routing_labels.tsv"), "w") as f:
        f.write("cell\tblastomere\tvotes_B1\tvotes_B2\tmasked_integration\treason\n")
        for i in range(n):
            bl = "B1" if is_b1[i] else "B2" if is_b2[i] else ""
            vb1, vb2 = (R["vL"][i], R["vR"][i]) if L_is_B1 else (R["vR"][i], R["vL"][i])
            mg = bcs[R["mask_g"][i]] if R["mask_g"][i] >= 0 else ""
            f.write(f"{cells[i]}\t{bl}\t{vb1}\t{vb2}\t{mg}\t{R['reason'][i]}\n")

    # ---- set-aside table ----
    with open(OP("set_aside.tsv"), "w") as f:
        f.write("cell\treason\tn_loci\tvotes_B1\tvotes_B2\n")
        for i in np.where(R["aside"])[0]:
            vb1, vb2 = (R["vL"][i], R["vR"][i]) if L_is_B1 else (R["vR"][i], R["vL"][i])
            f.write(f"{cells[i]}\t{R['reason'][i]}\t{nloci[i]}\t{vb1}\t{vb2}\n")

    # ---- the two matrices (with Tweak B pass_qc) ----
    qc_kw = dict(min_n_loci=args.min_n_loci, max_doublet_loci=args.max_doublet_loci,
                 min_mean_dom=args.min_mean_dominance)
    b1_tot, b1_pass = write_matrix(OP("B1_tape_consensus.tsv"), cells, bcs, tok, qc,
                                   b1_idx, R["mask_g"], **qc_kw)
    b2_tot, b2_pass = write_matrix(OP("B2_tape_consensus.tsv"), cells, bcs, tok, qc,
                                   b2_idx, R["mask_g"], **qc_kw)

    # ---- report ----
    n_mask = int((R["routed"] & (R["mask_g"] >= 0)).sum())
    disc_tapes = sorted(set(g for g, *_ in defs))
    lines = [
        "v4 blastomere routing report",
        "============================",
        f"input                : {args.consensus}",
        f"cells                : {n:,}",
        f"integrations         : {G}",
        f"discriminating tapes : {len(disc_tapes)}  ({[bcs[g] for g in disc_tapes]})",
        f"defining alleles     : {len(defs)}",
        f"orientation          : B1 = {'larger' if True else ''} routed group "
        f"(provisional {'L' if L_is_B1 else 'R'} side)",
        "",
        f"routed B1            : {len(b1_idx):,}",
        f"routed B2            : {len(b2_idx):,}",
        f"total routed         : {int(R['routed'].sum()):,} "
        f"({R['routed'].mean()*100:.1f}%)",
        f"masked (1 bleed->NA) : {n_mask:,}",
        "",
        f"set aside total      : {int(R['aside'].sum()):,} "
        f"({R['aside'].mean()*100:.2f}%)",
        f"  no founder (w=0)   : {int(R['none_'].sum()):,}",
        f"  tie                : {int(R['tie'].sum()):,}",
        f"  multibleed (l>=2)  : {int(R['multibleed'].sum()):,}",
        "",
        f"params: FIX_THRESH={FIX_THRESH} DISC_MIN={DISC_MIN} "
        f"CONF_MIN_LOCI={CONF_MIN_LOCI} FOUNDER_SITE=site{FOUNDER_SITE+1}",
        "",
        "Tweak B pre-ship QC filter (pass_qc; cells retained in matrix):",
        f"  thresholds         : n_loci>={args.min_n_loci}, "
        f"n_doublet_loci<={args.max_doublet_loci if args.max_doublet_loci < 10**9 else 'inf'}, "
        f"mean_dominance>={args.min_mean_dominance}",
        f"  B1 pass_qc=1       : {b1_pass:,} / {b1_tot:,} "
        f"({b1_pass*100/max(1,b1_tot):.1f}%)",
        f"  B2 pass_qc=1       : {b2_pass:,} / {b2_tot:,} "
        f"({b2_pass*100/max(1,b2_tot):.1f}%)",
        f"  total pass_qc=1    : {b1_pass+b2_pass:,} / {b1_tot+b2_tot:,} "
        f"({(b1_pass+b2_pass)*100/max(1,b1_tot+b2_tot):.1f}%)",
    ]
    report = "\n".join(lines) + "\n"
    with open(OP("routing_report.txt"), "w") as f:
        f.write(report)
    print(report)
    print(f"[route] wrote outputs to {args.outdir}/ (prefix '{args.out}')")


if __name__ == "__main__":
    main()
