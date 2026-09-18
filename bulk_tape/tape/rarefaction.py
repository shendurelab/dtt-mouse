"""Rarefaction / extrapolation of lineage-genotype diversity.

Given the abundances of distinct lineage genotypes at one integration (cells per
genotype, from the depth model), how many distinct genotypes are recovered as a
function of the number of cells sampled? We interpolate below the observed sample
size with the Hurlbert expectation E[S(m)] and extrapolate above it with the
Chao1 / iNEXT estimator, giving an accumulation curve and an asymptotic richness
(Chao1) against which the fraction captured is read off.

Abundances are integer counts (cells per genotype). Reference: Colwell et al.
2012 (iNEXT interpolation/extrapolation); Chao1 (Chao 1984).
"""
from __future__ import annotations
import math
import numpy as np


def chao1(counts):
    """Return (S_obs, S_chao1, f1, f2, f0) for a list of integer abundances.
    f1/f2 = # genotypes seen once/twice; f0 = estimated undetected genotypes."""
    counts = [c for c in counts if c > 0]
    n = sum(counts)
    S = len(counts)
    f1 = sum(1 for c in counts if c == 1)
    f2 = sum(1 for c in counts if c == 2)
    if n == 0:
        return 0, 0.0, 0, 0, 0.0
    f0 = ((n - 1) / n) * (f1 * f1 / (2 * f2) if f2 > 0 else f1 * (f1 - 1) / 2)
    return S, S + f0, f1, f2, f0


def _lcomb(n, k):
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def rarefy(counts, mgrid):
    """Expected number of distinct genotypes at each sample size m in mgrid.

    m <= n: exact Hurlbert interpolation E[S(m)] = S - sum_i C(n-c_i, m)/C(n, m).
    m  > n: Chao1 extrapolation S_obs + f0*(1 - (1 - f1/(n*f0+f1))^(m-n)).
    """
    counts = [c for c in counts if c > 0]
    n = sum(counts)
    S = len(counts)
    _, _, f1, f2, f0 = chao1(counts)
    out = []
    for m in mgrid:
        if m <= n:
            ln_den = _lcomb(n, m)
            s = 0.0
            for c in counts:
                lt = _lcomb(n - c, m)
                if lt > -math.inf:
                    s += math.exp(lt - ln_den)
            out.append(S - s)
        else:
            denom = n * f0 + f1
            inc = f0 * (1 - (1 - f1 / denom) ** (m - n)) if (denom > 0 and f0 > 0) else 0.0
            out.append(S + inc)
    return np.array(out)


def accumulation(counts, extrap: float, n_interp: int, n_extrap: int):
    """Build the sample-size grid and its rarefaction/extrapolation curve.

    Returns (mgrid, curve, n_obs) where n_obs is the observed number of cells."""
    n = sum(c for c in counts if c > 0)
    if n <= 0:
        return np.array([0]), np.array([0.0]), 0
    hi = max(n, int(round(n * extrap)))
    mgrid = np.unique(np.concatenate([
        np.linspace(1, n, n_interp),
        np.linspace(n, hi, n_extrap),
    ]).astype(int))
    return mgrid, rarefy(counts, mgrid), n
