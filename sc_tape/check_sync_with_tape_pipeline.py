#!/usr/bin/env python3
"""Guard against drift between this release and the working pipeline.

The modules in tape/ and the single-cell scripts are copies of the working
tree (mouse_sprint/tape_pipeline/ and mouse_sprint/src/0-blastomere-route/).
Copies drift. This asserts they have not.

  python3 check_sync_with_tape_pipeline.py [--pipeline ../tape_pipeline] \
                                           [--routing ../src/0-blastomere-route]

config.py is deliberately NOT compared -- this release carries a single-cell-only
config with the bulk-only knobs removed; the constants that both must share are
compared by value instead. blastomere_route.py is compared against the COMMIT
THAT PRODUCED THE PUBLISHED MATRICES (6ee057f), not against the working HEAD,
which replaced the multibleed rule with a stricter founder-conflict veto -- see
README, "Provenance".
"""
import argparse
import hashlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SHARED_MODULES = ["barcode.py", "denoise.py", "editchain.py", "io.py",
                  "merge.py", "parser.py", "sc.py", "__init__.py"]
SHARED_SCRIPTS = ["sc_consensus.py", "sc_consensus_par.py", "00_merge_pairs.py",
                  "sc_recovery.py", "sc_umi_stats.py"]
ROUTING_COMMIT = "6ee057f"
# constants that must agree with the working config.py, by value
SHARED_CONSTANTS = ["BC_RIGHT", "BC_RIGHT_SPACER", "BC_LEN", "N_SITES", "GAP_MIN",
                    "GAP_MAX", "ANCHOR", "TERM", "TERM_LONG", "MONO", "LEADIN",
                    "BC_CORRECT_MAX", "COLLAPSE_RATIO", "MIN_READS", "CROSS_BC_DEPTH",
                    "DROPOUT_FOLD_RATIO", "BRANCH_MIN_FRAC", "BRANCH_MIN_READS",
                    "UMI_LEN", "UMI_MERGE_DIST", "MIN_READS_PER_UMI",
                    "SC_MIN_MOLECULES", "SC_DOMINANCE_MIN",
                    "FASTP_OVERLAP_LEN_REQUIRE", "FASTP_OVERLAP_DIFF_LIMIT",
                    "FASTP_OVERLAP_DIFF_PCT_LIMIT"]


def md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pipeline", default=str(ROOT.parent / "tape_pipeline"))
    ap.add_argument("--routing", default=str(ROOT.parent / "src" / "0-blastomere-route"))
    a = ap.parse_args()
    pipe, routing = Path(a.pipeline), Path(a.routing)
    if not pipe.exists():
        sys.exit(f"working pipeline not found at {pipe} -- pass --pipeline. "
                 f"(Standalone checkouts have nothing to compare against; skip this check.)")

    bad = []
    for m in SHARED_MODULES:
        mine, theirs = ROOT / "tape" / m, pipe / "tape" / m
        if not theirs.exists():
            bad.append(f"MISSING upstream: tape/{m}")
        elif md5(mine) != md5(theirs):
            bad.append(f"DRIFT: tape/{m}")
    for s in SHARED_SCRIPTS:
        mine, theirs = ROOT / "scripts" / s, pipe / "scripts" / s
        if not theirs.exists():
            bad.append(f"MISSING upstream: scripts/{s}")
        elif md5(mine) != md5(theirs):
            bad.append(f"DRIFT: scripts/{s}")

    # blastomere_route.py against the published commit, via git
    if routing.exists():
        try:
            blob = subprocess.run(
                ["git", "-C", str(routing), "show",
                 f"{ROUTING_COMMIT}:src/0-blastomere-route/blastomere_route.py"],
                capture_output=True, check=True).stdout
            if hashlib.md5(blob).hexdigest() != md5(ROOT / "scripts" / "blastomere_route.py"):
                bad.append(f"DRIFT: scripts/blastomere_route.py vs {ROUTING_COMMIT}")
        except subprocess.CalledProcessError:
            bad.append(f"could not read {ROUTING_COMMIT} -- is the repo shallow?")

    # constants, by value
    sys.path.insert(0, str(ROOT))
    import config as mine_cfg
    sys.path.insert(0, str(pipe))
    for k in list(sys.modules):
        if k == "config":
            del sys.modules[k]
    sys.path.remove(str(ROOT))
    import config as their_cfg
    for c in SHARED_CONSTANTS:
        if not hasattr(their_cfg, c):
            bad.append(f"MISSING upstream constant: {c}")
        elif getattr(mine_cfg, c) != getattr(their_cfg, c):
            bad.append(f"CONSTANT DRIFT: {c} = {getattr(mine_cfg, c)!r} here, "
                       f"{getattr(their_cfg, c)!r} upstream")

    if bad:
        print("\n".join(bad))
        print(f"\n{len(bad)} problem(s). Re-sync the two copies before shipping.")
        sys.exit(1)
    print(f"in sync: {len(SHARED_MODULES)} modules, {len(SHARED_SCRIPTS)} scripts, "
          f"blastomere_route.py @ {ROUTING_COMMIT}, {len(SHARED_CONSTANTS)} constants")


if __name__ == "__main__":
    main()
