#!/usr/bin/env python3
"""Stage 02 -- parse reads to per-barcode genotype patterns.

Usage:  python3 scripts/02_parse.py DTTz_3_S3 [--limit N]

For every read: extract + whitelist-correct the TAPE-BC, resolve the 6-site
array with the robust parser, drop ordering violations, and keep VALID reads
(BC in whitelist, all sites resolved, insertions in vocab). Writes:
  data/<tag>.patterns_full.pkl        {barcode: {pattern tuple: read count}}
  data/<tag>.persite_editrate.tsv     overall per-site edit rate
  data/<tag>.perbc_persite_editrate.tsv
  data/<tag>.parse_summary.tsv        read funnel
Requires the raw fastq under RAW_DIR and the refs from stage 01.
"""
import sys
from pathlib import Path
from collections import Counter, defaultdict
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config
from tape.io import iter_reads, load_whitelist, load_vocab, save_pickle
from tape.barcode import extract_bc, correct
from tape.parser import walk_robust, pattern_of

N = config.N_SITES


def main(tag: str, limit=None):
    emb = config.get_embryo(tag)
    if not emb.fastq_path.exists():
        sys.exit(f"raw fastq not found: {emb.fastq_path} (see README).")
    if not emb.whitelist_tsv.exists():
        sys.exit(f"missing {emb.whitelist_tsv}; run stage 01 first.")
    wl = load_whitelist(emb.whitelist_tsv)
    vocab = load_vocab(emb.vocab_tsv)
    print(f"[02] parsing {emb.label} ({tag}); {len(wl)} whitelist BCs, {len(vocab)} vocab")

    site_ed = [0] * N
    site_res = [0] * N
    perbc_ed = {w: [0] * N for w in wl}
    perbc_res = {w: [0] * N for w in wl}
    patterns = defaultdict(Counter)
    tot = bc_ok = bc_wl = err = res_all = valid = 0

    for seq in iter_reads(emb.fastq_path, limit):
        tot += 1
        bc, cs = extract_bc(seq)
        if bc is None or len(bc) != config.BC_LEN:
            continue
        bc_ok += 1
        corr = correct(bc, wl)
        if corr is not None:
            bc_wl += 1
        state, e, ra = walk_robust(seq, cs)
        if e:
            err += 1
            continue
        if ra:
            res_all += 1
        for i in range(N):
            s = state[i]
            if s is None:
                continue
            site_res[i] += 1
            edited = (s != "U")
            if edited:
                site_ed[i] += 1
            if corr is not None:
                perbc_res[corr][i] += 1
                if edited:
                    perbc_ed[corr][i] += 1
        # valid read -> contribute a genotype pattern
        if corr is not None and ra:
            if all((s == "U") or (s[0] == "E" and s[1] in vocab) for s in state):
                valid += 1
                patterns[corr.decode()][pattern_of(state)] += 1

    save_pickle({k: dict(v) for k, v in patterns.items()}, emb.out("patterns_full.pkl"))

    with open(emb.out("persite_editrate.tsv"), "w") as f:
        f.write("site\tn_resolved\tn_edited\tedit_rate\n")
        for i in range(N):
            r = site_ed[i] / site_res[i] if site_res[i] else 0
            f.write(f"{i+1}\t{site_res[i]}\t{site_ed[i]}\t{r:.4f}\n")
    with open(emb.out("perbc_persite_editrate.tsv"), "w") as f:
        f.write("tapebc\tsite\tn_resolved\tn_edited\tedit_rate\n")
        for w in wl:
            for i in range(N):
                res, ed = perbc_res[w][i], perbc_ed[w][i]
                f.write(f"{w.decode()}\t{i+1}\t{res}\t{ed}\t{(ed/res if res else 0):.4f}\n")
    with open(emb.out("parse_summary.tsv"), "w") as f:
        f.write("metric\tcount\tpct_of_total\n")
        for name, v in [("total_reads", tot), ("bc_extractable", bc_ok),
                        ("bc_in_whitelist", bc_wl), ("ordering_errors", err),
                        ("all_sites_resolved", res_all), ("valid_reads", valid)]:
            f.write(f"{name}\t{v}\t{v/max(1,tot)*100:.2f}\n")

    site6 = site_ed[5] / site_res[5] * 100 if site_res[5] else 0
    print(f"[02] reads={tot:,}  bc_wl={bc_wl/max(1,tot)*100:.1f}%  valid={valid/max(1,tot)*100:.1f}%  "
          f"site6 edit={site6:.1f}%  barcodes={len(patterns)}")
    print(f"[02] wrote {emb.out('patterns_full.pkl')} + edit-rate/summary tsvs")


if __name__ == "__main__":
    args = sys.argv[1:]
    limit = None
    if "--limit" in args:
        i = args.index("--limit")
        limit = int(args[i + 1])
        args = args[:i] + args[i + 2:]
    if len(args) != 1:
        sys.exit(__doc__)
    main(args[0], limit)
