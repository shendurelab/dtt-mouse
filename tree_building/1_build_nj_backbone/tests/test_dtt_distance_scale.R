# test_dtt_distance_scale.R
# Prove that dtt_distance_matrix_fast() computes the EXACT same quantity as the
# reference dtt_distance_matrix() in dtt_distance.R. "Exact" means: identical NA
# pattern and a maximum absolute difference of exactly 0 (not a tolerance) -- both
# implementations reduce a pair to a mean of small integers, so they must agree
# to the bit.
#
# Run from the repo root (tree_building/), against the v6 data this repo ships:
#   DATA_VERSION=v6 DTT_TSV=processed_data/e3v5v6.B1_tape_consensus.ge7_founderok.tsv.gz \
#     Rscript 1_build_nj_backbone/tests/test_dtt_distance_scale.R
# (DATA_VERSION=v1/v2 also work if you have that delivery's data on disk --
#  DTT_TSV then defaults to that version's single consensus file.)

source("lib/paths.R")                              # DATA_VERSION, data_path
source("1_build_nj_backbone/dtt_distance.R")          # reference: dtt_distance_matrix
source("1_build_nj_backbone/dtt_distance_scale.R")    # fast:      dtt_distance_matrix_fast
source("1_build_nj_backbone/parse_tape_consensus.R")  # parse_cells (TSV -> dtt cells)

# DTT_TSV lets this test point at a per-side file directly (needed for v4/v5/v6,
# which have no single canonical consensus file -- consensus_tsv_path() only
# resolves a default for v1/v2/v3). Sys.getenv(var, unset) evaluates `unset`
# eagerly even when `var` IS set, so only call consensus_tsv_path() when needed.
TSV <- Sys.getenv("DTT_TSV")
if (!nzchar(TSV)) TSV <- consensus_tsv_path(".")
N_ROWS <- nrow(read.delim(TSV, colClasses = "character", check.names = FALSE))

# ---- test harness -----------------------------------------------------------
passed <- 0; failed <- 0

# Compare the reference and fast matrices for one set of cells. Asserts the NA
# pattern matches and the largest absolute value difference is exactly 0.
check_equiv <- function(name, cells) {
  ref  <- dtt_distance_matrix(cells)
  fast <- dtt_distance_matrix_fast(cells)
  na_ok   <- all(is.na(ref) == is.na(fast))                 # 1. same NA pattern
  dim_ok  <- identical(dimnames(ref), dimnames(fast))       # 2. same labels/order
  maxdiff <- max(abs(ref - fast), na.rm = TRUE)             # 3. largest difference
  if (!is.finite(maxdiff)) maxdiff <- 0                     #    (all-NA -> no diff)
  ok <- na_ok && dim_ok && maxdiff == 0
  if (ok) { passed <<- passed + 1; cat(sprintf("PASS  %-34s maxdiff %.1g\n", name, maxdiff)) }
  else    { failed <<- failed + 1
            cat(sprintf("FAIL  %-34s na_ok=%s dim_ok=%s maxdiff=%s\n",
                        name, na_ok, dim_ok, format(maxdiff))) }
}

# Assert one reference distance equals a hand-computed value (sanity-checks that
# the reference itself behaves as we believe before we trust it as ground truth).
check_value <- function(name, got, want) {
  ok <- (is.na(got) && is.na(want)) || (!is.na(got) && !is.na(want) && got == want)
  if (ok) { passed <<- passed + 1; cat(sprintf("PASS  %-34s = %s\n", name, format(got))) }
  else    { failed <<- failed + 1; cat(sprintf("FAIL  %-34s got %s want %s\n",
                                               name, format(got), format(want))) }
}

E <- "ETY"; N <- "None"  # unedited placeholder; missing site

# ============================================================================
# 1. The existing unit cases from test_dtt_distance.R
# ============================================================================
cellA <- list(
  tapeX = c("CTT", "ATA", "CAC", N, N, N),  # depth 3
  tapeY = c("CTT", N, N, N, N, N),          # depth 1
  tapeZ = c("GTA", "CCC", "CCC", N, N, N))  # depth 3
cellB <- list(
  tapeX = c("CTT", "ATA", "GGG", N, N, N),  # prefix 2 -> dist 2
  tapeY = c("CTT", N, N, N, N, N),          # prefix 1 -> dist 0
  tapeZ = c("GTA", "TTT", "AAA", N, N, N))  # prefix 1 -> dist 4
cellB_missingZ <- list(
  tapeX = c("CTT", "ATA", "GGG", N, N, N),
  tapeY = c("CTT", N, N, N, N, N),
  tapeZ = c(N, N, N, N, N, N))              # tapeZ undetected
check_equiv("unit: all 3 tapes present",  list(a = cellA, b = cellB))
check_equiv("unit: tapeZ missing in b",   list(a = cellA, b = cellB_missingZ))
check_value("unit: mean(2,0,4)/3",        dtt_distance(cellA, cellB),         (2 + 0 + 4) / 3)
check_value("unit: mean(2,0)/2",          dtt_distance(cellA, cellB_missingZ), (2 + 0) / 2)

# ============================================================================
# 2. Adversarial cases that exercise each subtle rule (each expects exact match)
# ============================================================================
# (a) no shared recovered tape -> NA. t1 recovered only in c1, t2 only in c2.
no_share <- list(
  c1 = list(t1 = c("AAG", E, E, E, E, E), t2 = c(N, N, N, N, N, N)),
  c2 = list(t1 = c(N, N, N, N, N, N),     t2 = c("AAG", E, E, E, E, E)))
check_equiv("adv: no shared tape",        no_share)
check_value("adv: no shared tape -> NA",  dtt_distance(no_share$c1, no_share$c2), NA_real_)

# (b) all-unedited shared tape -> depth 0 -> distance 0.
all_u <- list(c1 = list(t1 = rep(E, 6)), c2 = list(t1 = rep(E, 6)))
check_equiv("adv: all-unedited tape",     all_u)
check_value("adv: all-unedited -> 0",     dtt_distance(all_u$c1, all_u$c2), 0)

# (c) non-contiguous edits: depth counts ALL edited sites, not a prefix run.
noncontig <- list(
  c1 = list(t1 = c("AAG", E, "GCC", E, E, E)),   # depth 2
  c2 = list(t1 = c("AAG", E, "CAC", E, E, E)))   # depth 2, prefix 1 -> dist 2
check_equiv("adv: non-contiguous edits",  noncontig)
check_value("adv: non-contiguous -> 2",   dtt_distance(noncontig$c1, noncontig$c2), 2)

# (d) prefix matches, diverges, then a later site re-matches: re-match must NOT count.
rematch <- list(
  c1 = list(t1 = c("AAG", "GCC", "CAC", "TTT", E, E)),  # depth 4
  c2 = list(t1 = c("AAG", "GCC", "ACT", "TTT", E, E)))  # depth 4, prefix 2 -> dist 4
check_equiv("adv: re-match after divergence", rematch)
check_value("adv: re-match -> 4",         dtt_distance(rematch$c1, rematch$c2), 4)

# (e) variable-width edit tokens compared by full identity (GATG != GAT).
varwidth <- list(
  c1 = list(t1 = c("GATG", E, E, E, E, E)),  # depth 1
  c2 = list(t1 = c("GAT",  E, E, E, E, E)))  # depth 1, differ at site1 -> dist 2
check_equiv("adv: variable-width tokens", varwidth)
check_value("adv: GATG vs GAT -> 2",      dtt_distance(varwidth$c1, varwidth$c2), 2)

# (f) a cell with zero recovered tapes -> its diagonal is NA.
zero_rec <- list(
  dead  = list(t1 = c(N, N, N, N, N, N)),
  alive = list(t1 = c("AAG", E, E, E, E, E)))
check_equiv("adv: zero-recovered cell",   zero_rec)
m_zero <- dtt_distance_matrix(zero_rec)
check_value("adv: dead diagonal -> NA",   m_zero["dead", "dead"], NA_real_)

# ============================================================================
# 3. Real e3 subsamples (independent cross-check on actual data)
# ============================================================================
for (n in c(50L, 200L, 500L)) {
  set.seed(n)                                   # reproducible subsample per size
  rows <- sort(sample.int(N_ROWS, n))           # N_ROWS data rows in the TSV
  cells <- parse_cells(TSV, rows = rows)
  check_equiv(sprintf("real e3: %d random cells", n), cells)
}

# ---- summary ----------------------------------------------------------------
cat(sprintf("\n%d passed, %d failed\n", passed, failed))
if (failed > 0) quit(status = 1)
