"""Estimate genomic copy number per TAPE integration barcode.

In bulk, a barcode's total read abundance is proportional to (its genomic copy
number) x (number of input cells) x (per-molecule depth). Because every barcode
is sampled from the same cells at the same depth, per-barcode abundance is a
direct copy-number readout: single-copy integrations share one abundance level,
and a barcode carried as 2/3/... genomic copies (e.g. a concatemeric integration)
sits at the corresponding integer multiple.

This matters for two things:
  - the true number of INTEGRATIONS = sum of copies, which exceeds the number of
    distinct valid barcodes whenever any barcode is multi-copy;
  - cell counts: a multi-copy barcode's reads/lambda over-counts cells by its copy
    number, so per-cell estimates should divide by it.
"""
from __future__ import annotations
import numpy as np


def estimate_copies(per_bc_reads: dict):
    """Estimate integer genomic copy number per barcode from per-barcode total reads.

    Args:
      per_bc_reads: {barcode: total_reads}
    Returns:
      (copies, single_copy_unit, total_integrations)
        copies             = {barcode: int copies}
        single_copy_unit   = estimated reads for one genomic copy
        total_integrations = sum of copies over all barcodes
    """
    if not per_bc_reads:
        return {}, 0.0, 0
    reads = {b: float(r) for b, r in per_bc_reads.items()}
    # Single-copy unit = the abundance level of a single integration. Most barcodes
    # are single-copy, so start from the overall median, then refine to the median
    # of the barcodes that round to one copy (robust when several are multi-copy).
    unit = float(np.median(list(reads.values())))
    for _ in range(3):
        singles = [r for r in reads.values() if max(1, round(r / unit)) == 1]
        if not singles:
            break
        new = float(np.median(singles))
        if abs(new - unit) < 1e-6:
            break
        unit = new
    copies = {b: max(1, int(round(r / unit))) for b, r in reads.items()}
    return copies, unit, sum(copies.values())
