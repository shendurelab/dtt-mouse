# run_full_distance.R
# Build the full cell-cell DTT distance matrix over a set of QC-pass cells,
# append a synthetic founder-genotype outgroup for rooting, then build the NJ
# tree.
#
# Two output paths, chosen by which env var is set:
#   DTT_NJ_INMEM : hand the matrix straight to decenttree's RapidNJ in-process
#                  (production v6 path -- see run_b1.sh/run_b2.sh). No PHYLIP
#                  file, no intermediate parse.
#   DTT_OUT      : write a gzipped lower-triangular PHYLIP file instead, for
#                  decenttree's CLI to read separately (older path, kept for
#                  smaller/legacy DATA_VERSIONs).
#
# Key environment variables:
#   DATA_VERSION       : data delivery to use (required; see lib/paths.R)
#   DTT_TSV            : input consensus TSV (required for v6, which
#                        have no single canonical default)
#   DTT_NTHREADS       : kernel threads (default = all detected cores)
#   SYNTHROOT_SIDE     : "B1" or "B2" -- which blastomere's founder genotype
#                        the rooting outgroup is built from (required)
#   DEFINING_SITES_TSV : defining-sites table for that outgroup (default
#                        processed_data/e3v5v6.defining_sites.tsv)
#   QC_PASS_MODE       : "version" (apply DATA_VERSION's QC rule, default) or
#                        "all" (trust every cell_id in DTT_TSV as already
#                        QC-pass, e.g. this repo's ge7_founderok.tsv.gz files)
#   DTT_NJ_INMEM       : output .nwk path -- selects the in-memory NJ path
#   DTT_NJ_PRECISION   : decenttree float precision for that path (default 6)
#   DTT_NJ_NDIGITS     : decimal places to round distances to first (default 0)
#   DTT_DIR            : decenttree source dir (default <repo>/decenttree)
#   DTT_OUT            : output .dist.gz path for the PHYLIP path (default
#                        results/1-dist-matrix/<DATA_VERSION>/synthroot/...)
#   DTT_SCRATCH        : node-local scratch dir for that path's plain-text +
#                        pigz intermediate (default $TMPDIR or /tmp)
#
# e.g. (this repo's actual v6 production invocation -- see run_b1.sh/run_b2.sh):
#   DATA_VERSION=v6 SYNTHROOT_SIDE=B1 QC_PASS_MODE=all \
#     DTT_TSV=processed_data/e3v5v6.B1_tape_consensus.ge7_founderok.tsv.gz \
#     DEFINING_SITES_TSV=processed_data/e3v5v6.defining_sites.tsv \
#     DTT_NJ_INMEM=results/B1/e3v5v6_nj_ge7.nwk DTT_NJ_NDIGITS=0 \
#     Rscript 1_build_nj_backbone/run_full_distance.R

# Resolve this script's own directory so it sources its siblings and locates the
# repo's data/results dirs regardless of the working directory it is launched in.
.this_file <- local({
  a <- commandArgs(FALSE)
  m <- grep("^--file=", a, value = TRUE)               # set by Rscript
  if (length(m)) sub("^--file=", "", m[[1]]) else sys.frame(1)$ofile  # or source()
})
SCRIPT_DIR <- dirname(normalizePath(.this_file))        # <repo>/1_build_nj_backbone
REPO_DIR   <- dirname(SCRIPT_DIR)                       # <repo> (tree_building/)

source(file.path(REPO_DIR, "lib/paths.R"))            # DATA_VERSION, data_path, results_path, write_provenance
source(file.path(REPO_DIR, "lib/qc.R"))                # qc_pass_cell_ids
source(file.path(SCRIPT_DIR, "parse_tape_consensus.R"))  # parse_cells
source(file.path(SCRIPT_DIR, "dtt_distance_scale.R"))    # dtt_distance_matrix_fast
source(file.path(SCRIPT_DIR, "make_synthetic_root_outgroup.R"))  # make_synthetic_root_cell
Rcpp::sourceCpp(file.path(SCRIPT_DIR, "dtt_phylip_writer.cpp"))  # write_phylip_lower (plain text)

# ---- configuration ----------------------------------------------------------
SYNTHROOT_SIDE <- Sys.getenv("SYNTHROOT_SIDE")
if (!SYNTHROOT_SIDE %in% c("B1", "B2")) {
  stop(sprintf("SYNTHROOT_SIDE must be 'B1' or 'B2' (got '%s')", SYNTHROOT_SIDE))
}
DEFINING_SITES_TSV <- Sys.getenv("DEFINING_SITES_TSV", file.path(REPO_DIR, "processed_data/e3v5v6.defining_sites.tsv"))

QC_PASS_MODE <- Sys.getenv("QC_PASS_MODE", "version")  # version | all
if (!QC_PASS_MODE %in% c("version", "all")) {
  stop(sprintf("Unknown QC_PASS_MODE '%s' (known: version, all)", QC_PASS_MODE))
}

# Sys.getenv(var, unset) evaluates `unset` EAGERLY even when `var` IS set (R's
# usual lazy-argument evaluation does not apply to this base primitive) -- so
# calling consensus_tsv_path(REPO_DIR) unconditionally here would throw for
# DATA_VERSION=v4 (which has no single consensus file, by design; see
# lib/paths.R) even on a call that always overrides DTT_TSV anyway. Only
# resolve the default when DTT_TSV is actually unset.
TSV <- Sys.getenv("DTT_TSV")
if (!nzchar(TSV)) TSV <- consensus_tsv_path(REPO_DIR)
# "synthroot" subfolder keeps this path distinct from older, non-synthroot outputs.
DEFAULT_OUT <- results_path(REPO_DIR, "1-dist-matrix", "synthroot", "e3_distance_qcpass.dist.gz")
OUT    <- Sys.getenv("DTT_OUT", DEFAULT_OUT)
SCRATCH <- Sys.getenv("DTT_SCRATCH", Sys.getenv("TMPDIR", "/tmp"))
MAXC   <- as.integer(Sys.getenv("DTT_MAX_CELLS", "0"))
NTHR   <- as.integer(Sys.getenv("DTT_NTHREADS", "0"))
# DTT_NJ_INMEM (optional): if set to an output newick path, use the in-memory
# NJ path (below) instead of the PHYLIP write. Default "" -> PHYLIP path.
NJ_INMEM_OUT <- Sys.getenv("DTT_NJ_INMEM")
if (!nzchar(NJ_INMEM_OUT) && !nzchar(Sys.which("pigz"))) stop("pigz not found on PATH")

# decenttree_nj_inmem()/write_phylip_lower() both take a required `na_value` to
# impute NA distances with. For this repo's v6 data no NA can ever occur:
# dtt_distance() (dtt_distance.R) returns NA only when two cells share no
# recovered integration, but the upstream QC filter requires n_loci >= 7 of 11
# for every cell (lib/qc.R), so any two cells recovered sets overlap in at
# least 7 + 7 - 11 = 3 integrations (and the synthetic root outgroup below
# always has all 11 recovered). So this value is passed but can never be used.
NA_VALUE_UNUSED <- 0

# 1. cap kernel threads to the allocation (default = all cores if unset).
if (!is.na(NTHR) && NTHR > 0) RcppParallel::setThreadOptions(numThreads = NTHR)

# 2. positional indices of QC-pass rows. QC_PASS_MODE="version" applies this
#    DATA_VERSION's standard rule (lib/qc.R); "all" trusts DTT_TSV's own
#    row set as-is (for inputs that are already a pre-vetted subset).
all_ids <- read.delim(TSV, colClasses = "character", check.names = FALSE)$cell_id
if (QC_PASS_MODE == "all") {
  qc_idx <- seq_along(all_ids)
  cat("QC_PASS_MODE=all: using every cell_id in DTT_TSV, no additional filtering.\n")
} else {
  pass_ids <- qc_pass_cell_ids(REPO_DIR, DATA_VERSION, tsv_path = TSV)
  qc_idx  <- which(all_ids %in% pass_ids)
}
if (!is.na(MAXC) && MAXC > 0) qc_idx <- head(qc_idx, MAXC)   # quick-test cap
cat(sprintf("QC-pass cells to process: %d\n", length(qc_idx)))

# 3. parse those cells into the dtt representation.
cells <- parse_cells(TSV, rows = qc_idx)

# 4. append the synthetic root outgroup (SYNTHROOT_SIDE's founder genotype),
#    for rooting the tree downstream -- see make_synthetic_root_outgroup.R.
barcodes       <- names(cells[[1]])
defining_sites <- read.delim(DEFINING_SITES_TSV, colClasses = "character")
root_id        <- paste0("SYNTHETIC_ROOT_", SYNTHROOT_SIDE)
cells[[root_id]] <- make_synthetic_root_cell(barcodes, defining_sites, SYNTHROOT_SIDE)
cat(sprintf("%d cells total (incl. %s).\n", length(cells), root_id))

# 5. compute the full distance matrix (the big-memory step -- O(n^2) floats).
cat("Computing full distance matrix...\n")
t_mat <- system.time(M <- dtt_distance_matrix_fast(cells))["elapsed"]
cat(sprintf("  built %d x %d matrix in %.1f s\n", nrow(M), ncol(M), t_mat))
labels <- rownames(M)
rm(cells)                                  # free the parsed list before writing

# verify the n_loci>=7 no-NA guarantee (see NA_VALUE_UNUSED above) actually
# holds for this input, once, before either output path substitutes it in.
n_na <- sum(is.na(M))
if (n_na > 0) {
  stop(sprintf("FATAL: %d distance pairs are NA -- the n_loci>=7 QC guarantee (every pair shares >=3 integrations) does not hold for this input.", n_na))
}

# 5b. OPTIONAL in-memory NJ (DTT_NJ_INMEM): hand M straight to decenttree's
#     RapidNJ, skipping the text PHYLIP round-trip -- for memory reasons at
#     full scale (see run_b1.sh/run_b2.sh for sizing).
if (nzchar(NJ_INMEM_OUT)) {
  cat("DTT_NJ_INMEM set -> in-memory NJ-R handoff (no PHYLIP write, no decenttree parse).\n")
  DTT_DIR <- Sys.getenv("DTT_DIR", file.path(REPO_DIR, "decenttree"))

  if (Sys.info()[["sysname"]] == "Darwin") {
    omp_cxxflags <- ""
    omp_libs      <- ""
  } else {
    omp_cxxflags <- "-fopenmp"
    omp_libs      <- "-fopenmp"
  }
  Sys.setenv(PKG_CXXFLAGS = paste0("-I", normalizePath(DTT_DIR),
                                   " -I", normalizePath(file.path(DTT_DIR, "build")),
                                   " -O2 ", omp_cxxflags))
  Sys.setenv(PKG_LIBS = omp_libs)
  Rcpp::sourceCpp(file.path(SCRIPT_DIR, "nj_inmemory.cpp"))
  prec     <- as.integer(Sys.getenv("DTT_NJ_PRECISION", "6"))
  # ndigits: distance rounding before the NJ. Default 0 = FULL PRECISION 
  ndig     <- as.integer(Sys.getenv("DTT_NJ_NDIGITS", "0"))
  nthreads <- if (!is.na(NTHR) && NTHR > 0) NTHR else parallel::detectCores()
  dir.create(dirname(NJ_INMEM_OUT), showWarnings = FALSE, recursive = TRUE)
  cat(sprintf("In-memory NJ-R (nthreads=%d, precision=%d, ndigits=%d) -> %s\n",
              nthreads, prec, ndig, NJ_INMEM_OUT))
  t_nj <- system.time(
    decenttree_nj_inmem(M, labels, na_value = NA_VALUE_UNUSED, nthreads = nthreads,
                        precision = prec, out_path = NJ_INMEM_OUT, ndigits = ndig)
  )["elapsed"]
  cat(sprintf("  NJ built in %.1f s\n", t_nj))
  # same silent-tip-drop guard as the CLI driver: every input taxon must appear.
  nwk      <- readChar(NJ_INMEM_OUT, file.info(NJ_INMEM_OUT)$size, useBytes = TRUE)
  n_actual <- length(gregexpr("[(,][^,():]+:", nwk, perl = TRUE)[[1]])
  if (n_actual != length(labels)) {
    stop(sprintf("FATAL: in-memory NJ tree has %d tips, expected %d (decenttree dropped tips).",
                 n_actual, length(labels)))
  }
  cat(sprintf("In-memory NJ tree: tip count verified (%d == %d). Done.\n", n_actual, length(labels)))
  write_provenance(dirname(NJ_INMEM_OUT), REPO_DIR, inputs = TSV,
                   extra = list(n_cells = length(labels), synthroot_side = SYNTHROOT_SIDE,
                                qc_pass_mode = QC_PASS_MODE, nj_algo = "NJ-R", nj_mode = "in-memory"))
  quit(save = "no", status = 0)
}

# 6. write the lower-triangular PHYLIP.
dir.create(SCRATCH, showWarnings = FALSE, recursive = TRUE)
plain      <- tempfile("dtt_dist_", tmpdir = SCRATCH, fileext = ".phy")
scratch_gz <- paste0(plain, ".gz")
on.exit(unlink(c(plain, scratch_gz)), add = TRUE)


cat(sprintf("Formatting -> %s (local scratch)\n", plain))
t_fmt <- system.time(
  n_imp <- write_phylip_lower(M, labels, na_value = NA_VALUE_UNUSED, path = plain, ndigits = 4)
)["elapsed"]
cat(sprintf("  formatted in %.1f s\n", t_fmt))
stopifnot("n_na check above guarantees this" = n_imp == 0)

pigz_threads <- if (!is.na(NTHR) && NTHR > 0) NTHR else parallel::detectCores()
cat(sprintf("Compressing with pigz -p %d...\n", pigz_threads))
t_gz <- system.time(
  status <- system2("pigz", c("-p", pigz_threads, "-f", shQuote(plain)))
)["elapsed"]
if (status != 0) stop("pigz compression failed (exit status ", status, ")")
cat(sprintf("  compressed in %.1f s\n", t_gz))

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Copying -> %s\n", OUT))
if (!file.copy(scratch_gz, OUT, overwrite = TRUE))
  stop("failed to copy compressed output to ", OUT)

# 7. report.
cat(sprintf("Done. %d cells; %.0f NA pairs imputed (expected 0); file %.2f GB.\n",
            length(labels), n_imp,
            file.info(OUT)$size / 2^30))

# 8. provenance sidecar: trace this output back to the exact commit + inputs.
write_provenance(dirname(OUT), REPO_DIR, inputs = TSV,
                 extra = list(n_cells = length(labels), n_imputed = n_imp,
                              max_cells_cap = MAXC, synthroot_side = SYNTHROOT_SIDE,
                              qc_pass_mode = QC_PASS_MODE))
