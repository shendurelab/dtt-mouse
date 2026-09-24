# lib/paths.R
# Data-version config for the e3v8 pipeline (stage 1, 1_build_nj_backbone).
# dtt-mouse ships only the e3v8 delivery, so DATA_VERSION is always "v8"; this
# file centralizes where its files live so no stage script hardcodes them, and
# namespaces every result under results/<stage>/v8/... .
#
# Usage (from a stage script, after resolving REPO_DIR):
#   source(file.path(REPO_DIR, "lib/paths.R"))
#   OUT <- results_path(REPO_DIR, "1-dist-matrix", "e3_distance_qcpass.dist.gz")

DATA_VERSION <- Sys.getenv("DATA_VERSION", "")
if (!nzchar(DATA_VERSION)) {
  stop("DATA_VERSION must be set (v8) -- there is no default, so a result's ",
       "version is always explicit.")
}

.DATA_DIRS <- c(
  # v8 = e3v8: the tape-consensus delivery this repo ships. The RAW per-side
  # consensus (support_data/e3v8.{B1,B2}_tape_consensus.tsv.gz, repo root) and
  # the derived, founder-filtered ge7_founderok.tsv.gz backbone input
  # (tree_building/processed_data/, from 00_filter_highqual_consensus.R) live
  # in different directories -- see consensus_tsv_path() below.
  v8 = "processed_data"
)
if (!DATA_VERSION %in% names(.DATA_DIRS)) {
  stop(sprintf("Unknown DATA_VERSION '%s' (known: %s)",
               DATA_VERSION, paste(names(.DATA_DIRS), collapse = ", ")))
}

# Path to a file inside this run's data directory.
#   repo : repo root (e.g. REPO_DIR)
#   ...  : path components joined onto the data directory
data_path <- function(repo, ...) file.path(repo, .DATA_DIRS[[DATA_VERSION]], ...)

# v8 ships as TWO per-blastomere consensus matrices (B1/B2), not one, so there
# is no single canonical "the" consensus file to default to -- callers must
# always pass DTT_TSV/tsv_path explicitly (support_data/e3v8.{B1,B2}_tape_
# consensus.tsv.gz for the full callset, or processed_data/e3v8.{B1,B2}_tape_
# consensus.ge7_founderok.tsv.gz for the founder-filtered backbone input).
consensus_tsv_path <- function(repo) {
  stop("consensus_tsv_path: v8 has no single consensus file -- pass DTT_TSV/",
       "tsv_path explicitly.")
}

# v8's cell metadata lives in the repo root's support_data/, not under a
# per-version data directory.
metadata_path <- function(repo) file.path(repo, "..", "support_data", "cell_metadata.v8.txt.gz")

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
#             write into the SAME directory -- otherwise each step's
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
