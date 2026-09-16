#!/usr/bin/env python3
"""Per-blastomere clock rate: DTT (edit-distance) units per day.

The placer computes a query's pendant in DTT/edit units (dtt_lengths.py). To place
that pendant on the DATED backbone we divide by the per-side strict-clock rate LSD2
fitted when it dated that side:

    pendant_days = pendant_DTT / rate            (rate is DTT-units per day)

We take the rate LSD2 REPORTED (authoritative), not a value we recompute, because the
constrained (min-age) LSD2 adjusts branch times non-uniformly, so a naive
sum(DTT)/sum(time) slope would disagree with LSD2's fitted rate.

The rate is set via env RATE -- run_side.sh hardcodes it per side from the real v6
dating run's lsd2_final_rate (read once from that run's LSD2 output; not regenerated
at run time, since 3_date_tree does not produce a machine-readable provenance file).
"""
import os
import sys


def get_rate(side):
    """Resolve the per-side clock rate (DTT/day) and a provenance string.

    Aborts if RATE is unset so we never silently date with a missing/zero rate.
    """
    r = os.environ.get("RATE")
    if not r:
        sys.exit(f"[clockrate] no rate for {side}: set env RATE (DTT/day, from that "
                 f"side's 3_date_tree LSD2 run -- see run_side.sh).")
    return float(r), f"env RATE={r}"


if __name__ == "__main__":
    side = os.environ.get("SIDE", "B1")
    rate, src = get_rate(side)
    print(f"{side}\trate={rate:.6f} DTT/day\t1_DTT={1.0/rate:.4f} days\tsource: {src}")
