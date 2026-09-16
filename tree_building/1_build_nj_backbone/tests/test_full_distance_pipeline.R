# test_full_distance_pipeline.R
# End-to-end check of the full-distance pipeline on a small subset of QC-pass
# cells: build the matrix with the fast engine, impute NA -> median, write the
# lower-triangular PHYLIP as plain text with write_phylip_lower() and gzip it
# with pigz (the same local-scratch-then-compress path run_full_distance.R
# uses), read it back and confirm it round-trips the matrix, then (if the
# DecentTree binary is present) build a tree from it and confirm DecentTree
# ingests the format and returns all tips. This is the same code path as
# run_full_distance.R, just size-capped.
#
# Run from the repo root (tree_building/), against the v6 data this repo ships:
#   DATA_VERSION=v6 QC_PASS_MODE=all \
#     DTT_TSV=processed_data/e3v5v6.B1_tape_consensus.ge7_founderok.tsv.gz \
#     Rscript 1_build_nj_backbone/tests/test_full_distance_pipeline.R
# (DATA_VERSION=v1/v2 also work if you have that delivery's data on disk, with
#  QC_PASS_MODE=version -- see below.)
# Optional environment variables:
#   DATA_VERSION : v1, v2, or v6 (required, no default; see lib/paths.R)
#   DTT_TSV     : input consensus TSV (required for v4/v5/v6, which have no
#                 single canonical file; defaults to that DATA_VERSION's
#                 consensus file for v1/v2/v3)
#   QC_PASS_MODE : version (default; apply DATA_VERSION's QC rule via
#                  qc_pass_cell_ids) or all (trust every cell_id in DTT_TSV --
#                  for already-filtered inputs like this repo's
#                  ge7_founderok.tsv.gz files, which need DTT_TSV set)
#   DTT_TEST_N  : number of QC-pass cells to test on (default 300)
#   DECENTTREE  : path to the decenttree binary
#                 (default decenttree/build/decenttree; skipped if absent)

source("lib/paths.R")                                 # DATA_VERSION, data_path
source("lib/qc.R")                                     # qc_pass_cell_ids
source("1_build_nj_backbone/parse_tape_consensus.R")  # parse_cells
source("1_build_nj_backbone/dtt_distance_scale.R")    # dtt_distance_matrix_fast
Rcpp::sourceCpp("1_build_nj_backbone/dtt_phylip_writer.cpp")  # write_phylip_lower (plain text)
if (!nzchar(Sys.which("pigz"))) stop("pigz not found on PATH")

QC_PASS_MODE <- Sys.getenv("QC_PASS_MODE", "version")  # version | all
if (!QC_PASS_MODE %in% c("version", "all")) {
  stop(sprintf("Unknown QC_PASS_MODE '%s' (known: version, all)", QC_PASS_MODE))
}

# write the plain PHYLIP then gzip it with pigz, mirroring run_full_distance.R.
write_and_gzip <- function(M, labels, med) {
  plain <- tempfile(fileext = ".phy")
  on.exit(unlink(plain), add = TRUE)
  n_imp <- write_phylip_lower(M, labels, na_value = med, path = plain, ndigits = 6)
  status <- system2("pigz", c("-f", shQuote(plain)))
  if (status != 0) stop("pigz compression failed (exit status ", status, ")")
  list(n_imp = n_imp, gz = paste0(plain, ".gz"))
}

# Sys.getenv(var, unset) evaluates `unset` eagerly even when `var` IS set, so
# only call consensus_tsv_path() -- which errors for DATA_VERSION=v4/v5/v6 --
# when DTT_TSV is actually unset.
TSV  <- Sys.getenv("DTT_TSV")
if (!nzchar(TSV)) TSV <- consensus_tsv_path(".")
N    <- as.integer(Sys.getenv("DTT_TEST_N", "300"))
BIN  <- Sys.getenv("DECENTTREE", "decenttree/build/decenttree")
NWK  <- tempfile(fileext = ".nwk")

# ---- tiny test harness ------------------------------------------------------
passed <- 0; failed <- 0
ok <- function(name, cond, extra = "") {
  if (isTRUE(cond)) { passed <<- passed + 1; cat(sprintf("PASS  %-34s %s\n", name, extra)) }
  else              { failed <<- failed + 1; cat(sprintf("FAIL  %-34s %s\n", name, extra)) }
}

# 1. build the matrix for the first N QC-pass cells (same path as the pipeline).
#    QC_PASS_MODE=all trusts every cell_id in DTT_TSV as-is -- for inputs that
#    are already a pre-filtered subset (this repo's ge7_founderok.tsv.gz files),
#    same convention as run_full_distance.R.
all_ids <- read.delim(TSV, colClasses = "character", check.names = FALSE)$cell_id
if (QC_PASS_MODE == "all") {
  idx <- head(seq_along(all_ids), N)
} else {
  pass_ids <- qc_pass_cell_ids(".", DATA_VERSION, tsv_path = TSV)
  idx   <- head(which(all_ids %in% pass_ids), N)
}
cells <- parse_cells(TSV, rows = idx)
M     <- dtt_distance_matrix_fast(cells)

# 1b. verify the fully-unedited outgroup (make_unedited_cell) that
#     run_full_distance.R appends for rooting: its distance to a cell equals that
#     cell's mean edit-depth, it is never NA, and its self-distance is 0.
mean_depth <- function(cell) {                        # mean #edits over recovered tapes
  rec <- Filter(function(s) !all(s %in% c("None", "?")), cell)
  mean(vapply(rec, function(s) sum(!(s %in% c("None", "?", "ETY"))), numeric(1)))
}
cells_og <- c(cells, setNames(list(make_unedited_cell(names(cells[[1]]))), "UNEDITED_OUTGROUP"))
Mo  <- dtt_distance_matrix_fast(cells_og)
exp <- vapply(cells, mean_depth, numeric(1))          # expected outgroup distances
got <- Mo["UNEDITED_OUTGROUP", names(exp)]
ok("outgroup dist = cell mean depth", max(abs(got - exp)) <= 1e-9,
   sprintf("maxdiff %.1e", max(abs(got - exp))))
ok("outgroup never NA",        all(!is.na(Mo["UNEDITED_OUTGROUP", ])))
ok("outgroup self-distance 0", Mo["UNEDITED_OUTGROUP", "UNEDITED_OUTGROUP"] == 0)

# 2. median over finite off-diagonal pairs; impute NA -> median.
d   <- M[upper.tri(M)]
med <- median(d[!is.na(d)])
M[is.na(M)] <- med
cat(sprintf("Tested on %d cells; median = %.6f; NA pairs = %d\n",
            N, med, sum(is.na(d))))

# 3. write the plain lower-triangular PHYLIP file, then gzip it with pigz.
w <- write_and_gzip(M, rownames(M), med)
n_imp <- w$n_imp; GZ <- w$gz
ok("writer imputed-count matches", n_imp == sum(is.na(d)),
   sprintf("(%d)", as.integer(n_imp)))

# 4. read it back and reconstruct the full symmetric matrix.
ln    <- readLines(gzfile(GZ))
n     <- as.integer(ln[1])
rows  <- strsplit(trimws(ln[2:(n + 1)]), "[[:space:]]+")
back  <- matrix(0, n, n)
for (r in 2:n) back[r, 1:(r - 1)] <- as.numeric(rows[[r]][-1])
back[upper.tri(back)] <- t(back)[upper.tri(back)]
labs_file <- vapply(rows, `[`, character(1), 1)

ok("header count correct",        n == N)
ok("labels match matrix order",   identical(rownames(M), labs_file))
ok("matrix round-trips (<=1e-6)", max(abs(M - back)) <= 1e-6,
   sprintf("maxdiff %.1e", max(abs(M - back))))

# 4b. deterministically exercise NA imputation (a 300-cell subset may have no NA
#     pairs, so force two lower-triangle entries to NA and confirm they are
#     written as the median and counted).
Msyn <- M
Msyn[2, 1] <- NA; Msyn[1, 2] <- NA   # lower entry (2,1)
Msyn[3, 1] <- NA; Msyn[1, 3] <- NA   # lower entry (3,1)
w2  <- write_and_gzip(Msyn, rownames(Msyn), med)
ni2 <- w2$n_imp
r2  <- strsplit(trimws(readLines(gzfile(w2$gz))[2:(N + 1)]), "[[:space:]]+")
v21 <- as.numeric(r2[[2]][2]); v31 <- as.numeric(r2[[3]][2])
ok("synthetic NA imputed to median",
   ni2 == 2 && abs(v21 - med) <= 1e-6 && abs(v31 - med) <= 1e-6,
   sprintf("(count=%d)", as.integer(ni2)))

# 5. if DecentTree is available, build a tree and check it ingests the format.
if (nzchar(Sys.which(BIN)) || file.exists(BIN)) {
  status <- system2(BIN, c("-in", GZ, "-t", "NJ-R", "-out", NWK, "-nt", "1", "-no-banner"),
                    stdout = FALSE, stderr = FALSE)
  if (requireNamespace("ape", quietly = TRUE) && file.exists(NWK)) {
    tr <- ape::read.tree(NWK)
    ok("DecentTree returned all tips",   ape::Ntip(tr) == N, sprintf("(%d tips)", ape::Ntip(tr)))
    ok("DecentTree tip labels match",    setequal(tr$tip.label, rownames(M)))
  } else {
    cat("SKIP  DecentTree tree check (ape missing or no output)\n")
  }
} else {
  cat(sprintf("SKIP  DecentTree run (binary not found at '%s'; set DECENTTREE=)\n", BIN))
}

cat(sprintf("\n%d passed, %d failed\n", passed, failed))
if (failed > 0) quit(status = 1)
