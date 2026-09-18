#!/usr/bin/env python3
"""Emit tidy CSVs for the bulk-TAPE summary figure panels from existing stage
outputs (no re-parsing). Run after stages 02 + 04 for the embryos of interest.

Usage:  python3 scripts/bulk_figure_tables.py DTTz_3_S3 DTTz_2_S2

Writes to data/tables/:
  panel1_integration_summary.csv     per-embryo: n_barcodes, n_integrations, n_multicopy
  panel1_barcode_copies.csv          per-barcode: embryo, tapebc, reads, copies
  panel2_persite_by_embryo.csv       per-embryo x site: n_resolved, n_edited, edit_rate
  panel3_persite_by_integration.csv  per-barcode x site edit rate (each embryo)
"""
import sys, csv
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config

TAGS = sys.argv[1:] or ["DTTz_3_S3", "DTTz_2_S2"]
OUT = config.TABLES_DIR
OUT.mkdir(parents=True, exist_ok=True)


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


# ---- panel 1: integration barcode + copy-number summary (from stage 04) ----
p1_bc, p1_sum = [], []
for tag in TAGS:
    emb = config.get_embryo(tag)
    dm = emb.out("depth_model.tsv")
    if not dm.exists():
        print(f"[skip {tag}] missing {dm.name} (run stage 04)"); continue
    rows = read_tsv(dm)
    n_bc = len(rows)
    n_int = sum(int(r["copies"]) for r in rows)
    n_multi = sum(1 for r in rows if int(r["copies"]) > 1)
    # headline editing rates (site 1 and site 6) from the per-site table
    site_rate = {}
    ps = emb.out("persite_editrate.tsv")
    if ps.exists():
        site_rate = {r["site"]: r["edit_rate"] for r in read_tsv(ps)}
    p1_sum.append({"embryo": emb.label, "tag": tag, "n_barcodes": n_bc,
                   "n_integrations": n_int, "n_multicopy_barcodes": n_multi,
                   "site1_edit_rate": site_rate.get("1", ""),
                   "site6_edit_rate": site_rate.get("6", "")})
    for r in rows:
        p1_bc.append({"embryo": emb.label, "tapebc": r["tapebc"], "reads": r["reads"],
                      "copies": r["copies"], "cells": r["cells"]})

# ---- panel 2: per-site edit rate per embryo (from stage 02) ----
p2 = []
for tag in TAGS:
    emb = config.get_embryo(tag)
    ps = emb.out("persite_editrate.tsv")
    if not ps.exists():
        continue
    for r in read_tsv(ps):
        p2.append({"embryo": emb.label, "site": r["site"], "n_resolved": r["n_resolved"],
                   "n_edited": r["n_edited"], "edit_rate": r["edit_rate"]})

# ---- panel 3: per-integration per-site edit rate (from stage 02) ----
p3 = []
for tag in TAGS:
    emb = config.get_embryo(tag)
    pb = emb.out("perbc_persite_editrate.tsv")
    if not pb.exists():
        continue
    for r in read_tsv(pb):
        p3.append({"embryo": emb.label, "tapebc": r["tapebc"], "site": r["site"],
                   "edit_rate": r["edit_rate"]})


def write_csv(name, rows, cols):
    path = OUT / name
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
    print(f"wrote {path}  ({len(rows)} rows)")


write_csv("panel1_integration_summary.csv", p1_sum,
          ["embryo", "tag", "n_barcodes", "n_integrations", "n_multicopy_barcodes",
           "site1_edit_rate", "site6_edit_rate"])
write_csv("panel1_barcode_copies.csv", p1_bc, ["embryo", "tapebc", "reads", "copies", "cells"])
write_csv("panel2_persite_by_embryo.csv", p2,
          ["embryo", "site", "n_resolved", "n_edited", "edit_rate"])
write_csv("panel3_persite_by_integration.csv", p3, ["embryo", "tapebc", "site", "edit_rate"])
