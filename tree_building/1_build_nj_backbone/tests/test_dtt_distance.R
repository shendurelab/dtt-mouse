# test_dtt_distance.R
# Run from the repo root (tree_building/):  Rscript 1_build_nj_backbone/tests/test_dtt_distance.R
#
# Two cells, three tapes each (tapeX, tapeY, tapeZ). The per-tape distance is
#   depth1 + depth2 - 2 * shared_edit_prefix.

source("1_build_nj_backbone/dtt_distance.R")

# ---- tiny test harness ------------------------------------------------------
passed <- 0; failed <- 0
run_test <- function(name, got, want, tol = 1e-9) {
  ok <- !is.na(got) && abs(got - want) < tol
  if (ok) { passed <<- passed + 1; cat(sprintf("PASS  %-40s got %.4f\n", name, got)) }
  else    { failed <<- failed + 1; cat(sprintf("FAIL  %-40s got %s, want %.4f\n",
                                               name, format(got), want)) }
}

# ---- shared definitions -----------------------------------------------------
# cellA: same in both tests.
cellA <- list(
  tapeX = c("CTT", "ATA", "CAC", "None", "None", "None"),  # depth 3
  tapeY = c("CTT", "None", "None", "None", "None", "None"),  # depth 1
  tapeZ = c("GTA", "CCC", "CCC", "None", "None", "None")   # depth 3
)

# Per-tape distances of cellA vs cellB below (hand-computed):
#   tapeX: depth 3 & 3, prefix 2 (CTT,ATA match; CAC vs GGG differ) -> 3+3-2*2 = 2
#   tapeY: depth 1 & 1, prefix 1 (identical)                        -> 1+1-2*1 = 0
#   tapeZ: depth 3 & 3, prefix 1 (GTA match; CCC vs TTT differ)     -> 3+3-2*1 = 4

# ============================================================================
# TEST 1: all three tapes present in both cells.
#   distance = mean(2, 0, 4) = 2
# ============================================================================
cellB <- list(
  tapeX = c("CTT", "ATA", "GGG", "None", "None", "None"),  # depth 3
  tapeY = c("CTT", "None", "None", "None", "None", "None"),  # depth 1
  tapeZ = c("GTA", "TTT", "AAA", "None", "None", "None")   # depth 3
)
run_test("test1: all 3 tapes present", dtt_distance(cellA, cellB), (2 + 0 + 4) / 3)

# ============================================================================
# TEST 2: tapeZ undetected in the second cell (all sites "None").
#   Only tapeX and tapeY are shared -> distance = mean(2, 0) = 1
#   (Dropping the most-divergent tape lowers the distance AND the denominator
#    goes 3 -> 2; the missing integration is excluded, not penalized.)
# ============================================================================
cellB_missingZ <- list(
  tapeX = c("CTT", "ATA", "GGG", "None", "None", "None"),  # depth 3
  tapeY = c("CTT", "None", "None", "None", "None", "None"),  # depth 1
  tapeZ = c("None", "None", "None", "None", "None", "None")  # undetected
)
run_test("test2: tapeZ missing in cell B", dtt_distance(cellA, cellB_missingZ), (2 + 0) / 2)

# ---- summary ----------------------------------------------------------------
cat(sprintf("\n%d passed, %d failed\n", passed, failed))
if (failed > 0) quit(status = 1)
