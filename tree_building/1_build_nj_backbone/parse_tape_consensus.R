# parse_tape_consensus.R
# Read the e3 TAPE consensus TSV into the cell representation that
# dtt_distance.R expects: a named list of cells, where each cell is a named
# list of integrations (tapes) and each tape is a length-6 character vector of
# site states. The list name is the integration id (barcode).
#
# Token mapping from the TSV to dtt_distance.R conventions:
#   - a whole-cell integration value of NA (integration NOT observed)
#                                         -> a tape of six "None" (all missing,
#                                            so tape_recovered() == FALSE)
#   - within an observed tape "s1|...|s6", a site "U" (UNEDITED)
#                                         -> "ETY" (the unedited placeholder)
#   - any other site token (an edit, e.g. "AAG", "GATG")  -> kept verbatim
#
# Run (example):
#   source("lib/paths.R")   # DATA_VERSION must be set, e.g. DATA_VERSION=v1
#   source("1_build_nj_backbone/parse_tape_consensus.R")
#   cells <- parse_cells(consensus_tsv_path("."), rows = 1:200)

MISSING_TAPE <- rep("None", 6)  # six missing sites == integration not recovered

# Convert one TSV integration cell into a length-6 tape vector.
parse_tape <- function(value) {
  # 1. an unobserved integration is read as NA -> all-missing tape
  if (is.na(value)) return(MISSING_TAPE)
  # 2. split the "s1|s2|...|s6" string into its six site tokens
  sites <- strsplit(value, "|", fixed = TRUE)[[1]]
  # 3. guard the data contract: every observed tape has exactly six sites
  if (length(sites) != 6L) {
    stop(sprintf("expected 6 sites, got %d in tape '%s'", length(sites), value))
  }
  # 4. map the unedited placeholder "U" to dtt's "ETY"; edits pass through
  sites[sites == "U"] <- "ETY"
  sites
}

# Parse the consensus TSV into a named list of cells (the dtt representation).
#   tsv_path : path to e3_tape_consensus.tsv
#   rows     : optional integer vector to keep only some data rows (for tests)
parse_cells <- function(tsv_path, rows = NULL) {
  # 1. read every column as character so site tokens are never coerced
  df <- read.delim(tsv_path, colClasses = "character", check.names = FALSE)
  # 2. optionally subset to the requested data rows
  if (!is.null(rows)) df <- df[rows, , drop = FALSE]
  # 3. the 11 integration barcodes are columns 2..12 (col 1 is cell_id)
  barcodes <- colnames(df)[2:12]
  # 4. build one cell (named list of 11 tapes) per row
  cells <- lapply(seq_len(nrow(df)), function(r) {
    tapes <- lapply(barcodes, function(b) parse_tape(df[[b]][r]))
    setNames(tapes, barcodes)
  })
  # 5. name each cell by its cell_id
  setNames(cells, df[["cell_id"]])
}

# A fully-unedited cell: every integration recovered but all six sites unedited
# ("ETY", depth 0). Used as an outgroup to root the tree -- it is the ancestral
# genotype, so its distance to any cell equals that cell's mean edit-depth.
#   barcodes : the integration ids to give the cell (use names() of a real cell
#              so the outgroup shares all integrations with every other cell)
make_unedited_cell <- function(barcodes) {
  setNames(rep(list(rep("ETY", 6L)), length(barcodes)), barcodes)
}
