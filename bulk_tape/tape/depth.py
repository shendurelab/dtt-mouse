"""Read-depth model: turn per-pattern read counts into genomic equivalents (cells).

If each genomic equivalent (one integration in one cell) yields ~lambda reads,
then a clean pattern's reads/lambda estimates how many cells carry it. lambda is
the mode of the per-pattern read-count distribution in log space -- the "single
genomic equivalent" hump, below the 2x/3x/... clonal multiplets.
"""
from __future__ import annotations
import numpy as np
from scipy.stats import gaussian_kde

from config import LAMBDA_MIN_READS, LAMBDA_KDE_BW


def estimate_lambda(counts) -> float:
    """Estimate lambda from an iterable of per-pattern read counts.

    Uses a Gaussian KDE in log10 space over patterns with >= LAMBDA_MIN_READS
    reads (the canonical method behind the original cell_data.json). Falls back
    to a smoothed histogram mode when too few patterns clear the floor.
    """
    vals = np.asarray([x for x in counts if x >= LAMBDA_MIN_READS], float)
    if vals.size >= 5:
        v = np.log10(vals)
        if v.max() > v.min():
            kde = gaussian_kde(v, bw_method=LAMBDA_KDE_BW)
            xs = np.linspace(v.min(), v.max(), 400)
            return float(10 ** xs[np.argmax(kde(xs))])
    return _hist_mode(counts)


def _hist_mode(counts) -> float:
    """Fallback lambda: peak of a smoothed 40-bin log10 histogram over all counts."""
    vals = np.asarray(list(counts), float)
    vals = vals[vals > 0]
    if vals.size == 0:
        return 1.0
    v = np.log10(vals)
    if v.max() == v.min():
        return float(vals[0])
    edges = np.linspace(v.min(), v.max(), 40)
    h, e = np.histogram(v, bins=edges)
    hs = np.convolve(h, np.ones(3) / 3, mode="same")
    centers = (e[:-1] + e[1:]) / 2
    return float(10 ** centers[np.argmax(hs)])


def cells_per_pattern(counts: dict, lam: float, min_cells: int = 1) -> dict:
    """Map {pattern: reads} -> {pattern: n_cells} with n_cells = round(reads/lam),
    keeping patterns with >= min_cells."""
    out = {p: int(round(c / lam)) for p, c in counts.items()}
    return {p: n for p, n in out.items() if n >= min_cells}


def genomic_equivalents(counts, lam: float) -> int:
    """Total genomic equivalents in a barcode = sum of round(reads/lambda)."""
    return int(np.round(np.asarray(list(counts), float) / lam).sum())
