# lib/paths.R
# Shared data-version config for the e3 pipeline (stages 1-4). Sourced by every
# stage script so the data-directory mapping and the results-path convention
# live in exactly one place instead of being copy-pasted per script.
#
# DATA_VERSION selects which delivery of the raw data a run uses (e.g. the
# tape-consensus TSV changed shape between v1 and v2 -- see lib/qc.R). It
# must be set explicitly (no silent default), and it namespaces every output
# path so v1 and v2 results/provenance never collide or get mixed up:
#   results/<stage>/<DATA_VERSION>/...
#
# Usage (from a stage script, after resolving REPO_DIR):
#   source(file.path(REPO_DIR, "lib/paths.R"))
#   TSV <- data_path(REPO_DIR, "e3_tape_consensus.tsv")
#   OUT <- results_path(REPO_DIR, "1-dist-matrix", "e3_distance_qcpass.dist.gz")

DATA_VERSION <- Sys.getenv("DATA_VERSION", "")
if (!nzchar(DATA_VERSION)) {
  stop("DATA_VERSION must be set (e.g. DATA_VERSION=v1 or DATA_VERSION=v2) -- ",
       "there is no default, so a result's version is always explicit.")
}

# version -> data directory, relative to the repo root. Centralizes the v1/v2
# naming inconsistency (v1-e3_tape vs v2-e3-tape) so no other script has to
# know about it.
.DATA_DIRS <- c(
  v1 = "data/nobackup/v1-e3_tape",
  v2 = "data/nobackup/v2-e3-tape",
  v3 = "data/nobackup/v3-e3-tape",
  # v4 = blastomere-routed B1/B2 split of v3's genotype calls (masks 7,275
  # cross-bleed tape calls to NA at the routing step -- a real genotype-call
  # change, not just a subset, so it gets its own version rather than being
  # folded into v3). Lives in tape_pipeline/tables/v4 (git-tracked already;
  # unlike v1-v3 there is no data/nobackup copy to keep in sync).
  v4 = "tape_pipeline/tables/v4",
  # v5 = e3v5, a fresh STANDALONE single-cell call for embryo #3 from the seq4
  # run (NOT pooled with v3's e3 matrix). Like v4 it ships as two per-blastomere
  # B1/B2 tape-consensus matrices with a pass_qc column (git-tracked in
  # tape_pipeline/tables/v5; no data/nobackup copy). See lib/qc.R for its
  # QC rule and tape_pipeline/tables/v5/README.md for provenance.
  v5 = "tape_pipeline/tables/v5",
  # v6 = e3v5v6, the seq4 e3v5 call REPROCESSED with the v6 tape-consensus
  # delivery (processed_data/e3v5v6.{B1,B2}_tape_consensus.tsv.gz). Same
  # two-per-blastomere shape as v4/v5 (11 integration barcodes + n_loci/
  # n_doublet_loci/mean_dominance/pass_qc), a genuinely different delivery from
  # v5 -- so its results/provenance MUST be namespaced v6, not v5. This is the
  # high-quality backbone rebuild input (pass_qc & n_loci>=7 & founder-consistent;
  # see lib/qc.R for its QC rule and processed_data/README.md).
  v6 = "processed_data"
)
if (!DATA_VERSION %in% names(.DATA_DIRS)) {
  stop(sprintf("Unknown DATA_VERSION '%s' (known: %s)",
               DATA_VERSION, paste(names(.DATA_DIRS), collapse = ", ")))
}

# Path to a file inside this run's data directory.
#   repo : repo root (e.g. REPO_DIR)
#   ...  : path components joined onto the data directory
data_path <- function(repo, ...) file.path(repo, .DATA_DIRS[[DATA_VERSION]], ...)

# The tape-consensus file itself changed name AND compression between
# versions (v1: e3_tape_consensus.tsv; v2: e3_tape_consensus.tsv.gz), so
# callers must go through this instead of hardcoding the filename -- R's
# read.delim()/read.table() auto-detect gzip from the file content once the
# path is correct, so nothing else needs to change at the call site.
#
# v4 deliberately has NO entry here: it ships as TWO parallel matrices
# (e3.B1_tape_consensus.tsv.gz, e3.B2_tape_consensus.tsv.gz), not one, so
# there is no single canonical "the" consensus file to default to -- callers
# working with v4 must always pass DTT_TSV explicitly (see
# src/13-blastomere-route/README.md).
.CONSENSUS_FILE <- c(v1 = "e3_tape_consensus.tsv", v2 = "e3_tape_consensus.tsv.gz", v3 = "e3_tape_consensus.tsv.gz")
consensus_tsv_path <- function(repo) data_path(repo, .CONSENSUS_FILE[[DATA_VERSION]])

# The canonical cell-metadata file changed name between versions (v1's
# delivery calls it cell_metadata.v2.txt; v2's calls it cell_metadata.v3.txt;
# v5's calls it cell_metadata.v4.txt; v6's calls it cell_metadata.v6.txt --
# the numberings are unrelated), so callers must go through this instead of
# hardcoding the filename.
.METADATA_FILE <- c(v1 = "cell_metadata.v2.txt", v2 = "cell_metadata.v3.txt",
                    v5 = "cell_metadata.v4.txt", v6 = "cell_metadata.v6.txt")
# v5/v6's metadata are big single-cell tables like v1-v3's, so they live in the
# gitignored data/nobackup tree -- NOT alongside their git-tracked per-blastomere
# tape tables (.DATA_DIRS[["v5"]] = tape_pipeline/tables/v5, etc). So their
# directories are resolved separately from data_path() here; other versions
# keep the old behavior. v6's delivery (1,753,895 cells) covers all 1,708,269
# v6-called cells, including the seq5 batch that v5's metadata predates.
.METADATA_DIRS <- c(v5 = "data/nobackup/v5-e3-tape", v6 = "data/nobackup/v6")
metadata_path <- function(repo) {
  dir <- if (DATA_VERSION %in% names(.METADATA_DIRS)) .METADATA_DIRS[[DATA_VERSION]]
         else .DATA_DIRS[[DATA_VERSION]]
  file.path(repo, dir, .METADATA_FILE[[DATA_VERSION]])
}

# Path to a file inside this run's version-namespaced results directory for a
# given stage, e.g. results_path(repo, "1-dist-matrix", "x.dist.gz") ->
# <repo>/results/1-dist-matrix/<DATA_VERSION>/x.dist.gz
results_path <- function(repo, stage, ...) {
  file.path(repo, "results", stage, DATA_VERSION, ...)
}

# Write a provenance sidecar next to a stage's output, so a result can be
# traced back to the exact code, data, and parameters that produced it (the
# directory name only tells you the data VERSION, not the commit or inputs).
#   out_dir : directory to write provenance.json into (created if needed)
#   repo    : repo root, for `git -C <repo> ...`
#   inputs  : character vector of input file paths to checksum
#   extra   : named list of additional fields (e.g. key env vars, seeds)
#   tag     : distinguishes the sidecar filename when multiple pipeline steps
#             write into the SAME directory (e.g. stage 4's n<N>_s<seed>/, which
#             02/03/04 and run_sciphy.sh all populate) -- otherwise each step's
#             provenance.json would silently overwrite the previous step's.
#             Default "provenance" is fine when a step owns its output dir alone.
write_provenance <- function(out_dir, repo, inputs = character(0), extra = list(),
                             tag = "provenance") {
  git_commit <- tryCatch(
    system2("git", c("-C", repo, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_)
  git_dirty <- tryCatch(
    length(system2("git", c("-C", repo, "status", "--porcelain"), stdout = TRUE, stderr = FALSE)) > 0,
    error = function(e) NA)

  input_info <- lapply(inputs, function(p) {
    list(path = p, exists = file.exists(p),
         md5 = if (file.exists(p)) unname(tools::md5sum(p)) else NA_character_)
  })
  names(input_info) <- inputs

  info <- c(
    list(
      timestamp    = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      git_commit   = git_commit,
      git_dirty    = git_dirty,
      data_version = DATA_VERSION,
      r_version    = R.version.string,
      hostname     = Sys.info()[["nodename"]],
      sge_job_id   = Sys.getenv("JOB_ID", NA),
      inputs       = input_info
    ),
    extra
  )

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(info, file.path(out_dir, paste0(tag, ".json")),
                         auto_unbox = TRUE, null = "null", na = "null")
  } else {
    # jsonlite not available -- fall back to a flat, still-diffable text dump
    # rather than silently skipping provenance.
    dump <- utils::capture.output(utils::str(info))
    writeLines(dump, file.path(out_dir, paste0(tag, ".txt")))
  }
}
