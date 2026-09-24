# 00_filter_highqual_consensus.R
# =============================================================================
# Produce the HIGH-QUALITY backbone input for the v6 NJ rebuild: filter a v6
# per-side tape-consensus TSV down to the cells passing ALL of
#   (1) pass_qc == 1
#   (2) n_loci >= MIN_LOCI            (default 7)
#   (3) FOUNDER-CONSISTENT            (NOT divergent from this side's founder)
# and write the surviving rows back out in the SAME column format the NJ build
# (src/1-dist-mat/run_full_distance.R) reads, so it can be fed as DTT_TSV with
# QC_PASS_MODE=all (filtering already applied).
#
# The founder-consistency gate is a byte-for-byte reuse of the rule in
# src/16-placement/04_prep_v6_queries.R (itself the audited copy of
# src/8-support/15_divergent_cell_filter.R): parse the canonical founder config
# D_b from the supp-table-1 genotypes (B1 -> A_founder_genotype, B2 ->
# B_founder_genotype), then DROP a cell if at ANY defining (barcode, site) it
# OBSERVED a non-"?" token that DISAGREES with the founder token (a divergent
# edit, or an unedited "0" where the founder is edited). Whole-missing tapes
# ("?") never cause a drop.
#
# This script does NOT touch pars_ancestral_state or any reconstruction code --
# it only changes which INPUT cells the backbone is built from.
#
# Usage (env vars; defaults shown are dtt-mouse's e3v8 case):
#   SIDE=B1 MIN_LOCI=7 Rscript tree_building/1_build_nj_backbone/00_filter_highqual_consensus.R
#
# Expected survivor counts (fail-fast checkpoint, step 6), e3v8:
#   MIN_LOCI=7 (the NJ backbone input) -> B1 410,925 ; B2 244,776
#   MIN_LOCI=4 (the placement QUERY input -- pass_qc==1 never has n_loci<4, so
#     this is "no additional floor beyond pass_qc"; B1 750,417 ;
#     B2 530,861. 
# =============================================================================

suppressMessages({ library(data.table) })

# ---- 0. config -------------------------------------------------------------
SIDE        <- Sys.getenv("SIDE", "B1")
stopifnot("SIDE must be B1 or B2" = SIDE %in% c("B1", "B2"))
MIN_LOCI    <- as.integer(Sys.getenv("MIN_LOCI", "7"))
QUERY_TSV   <- Sys.getenv("QUERY_TSV",
                 sprintf("support_data/e3v8.%s_tape_consensus.tsv.gz", SIDE))
FOUNDER_CSV <- Sys.getenv("FOUNDER_CSV",
                 "tree_building/processed_data/e3v8.supp_table1_founder_genotypes.wide.csv")
OUT         <- Sys.getenv("OUT",
                 if (MIN_LOCI == 4L)
                   sprintf("tree_building/processed_data/e3v8.%s_tape_consensus.founderok.tsv.gz", SIDE)
                 else
                   sprintf("tree_building/processed_data/e3v8.%s_tape_consensus.ge%d_founderok.tsv.gz",
                           SIDE, MIN_LOCI))
SIDE_KEY    <- if (SIDE == "B1") "A" else "B"       # founder table uses A/B for B1/B2
SITE_COLS   <- paste0("Site", 1:6)
# expected founder-consistent survivor counts, for the fail-fast checkpoint.
# dtt-mouse is e3v8-only, so the default IS the e3v8 target (unlike
# mouse_sprint's script, which defaults to v6 and needs EXPECT_KEPT_B1/B2
# overridden for v8); SKIP_EXPECT_KEPT=1 still available for a first
# exploratory run against different data.
EXPECT_KEPT <- (if (MIN_LOCI == 7L) c(B1 = 410925L, B2 = 244776L)
               else if (MIN_LOCI == 4L) c(B1 = 750417L, B2 = 530861L)
               else c(B1 = NA, B2 = NA))
if (nzchar(Sys.getenv("EXPECT_KEPT_B1", ""))) EXPECT_KEPT["B1"] <- as.integer(Sys.getenv("EXPECT_KEPT_B1"))
if (nzchar(Sys.getenv("EXPECT_KEPT_B2", ""))) EXPECT_KEPT["B2"] <- as.integer(Sys.getenv("EXPECT_KEPT_B2"))
if (Sys.getenv("SKIP_EXPECT_KEPT", "0") == "1") EXPECT_KEPT[] <- NA
for (f in c(QUERY_TSV, FOUNDER_CSV)) stopifnot("missing input" = file.exists(f))
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)

# fread's default na.strings turns the literal "NA" (a whole-missing barcode
# cell) into <NA> -- exactly what 04_prep_v6_queries.R relies on for the gate,
# so we read the SAME way to reproduce its counts, and write back with na="NA"
# (step 5) to restore the original file format byte-for-byte.
read_tsv <- function(p) {
  if (grepl("\\.gz$", p)) fread(cmd = sprintf("zcat < %s", shQuote(p)), colClasses = "character")
  else fread(p, colClasses = "character")
}

# ---- 1. load the full per-side consensus -----------------------------------
dt      <- read_tsv(QUERY_TSV)
bcs     <- colnames(dt)[2:12]                       # the 11 integration barcodes
n_loci  <- as.integer(dt$n_loci)
pass_qc <- as.integer(dt$pass_qc)
cat(sprintf("[1] %s: %d rows loaded from %s\n", SIDE, nrow(dt), QUERY_TSV))

# ---- 2. quality band: pass_qc==1 & n_loci>=MIN_LOCI ------------------------
qc_ok <- pass_qc == 1L & !is.na(n_loci) & n_loci >= MIN_LOCI
cat(sprintf("[2] pass_qc==1 & n_loci>=%d: %d cells\n", MIN_LOCI, sum(qc_ok)))

# ---- 3. parse the canonical founder config D_b (identical to 04_/15_) ------
ft <- read.csv(FOUNDER_CSV, comment.char = "#", stringsAsFactors = FALSE,
               colClasses = "character")
parse_founder <- function(integrations, genotypes) {
  rows <- list()
  for (i in seq_along(integrations)) {
    g <- genotypes[i]
    if (is.na(g) || g == "(unedited)" || g == "") next
    toks <- strsplit(g, "-", fixed = TRUE)[[1]]
    for (k in seq_along(toks)) rows[[length(rows) + 1L]] <-
      data.frame(barcode = integrations[i], site = SITE_COLS[k], founder_token = toks[k],
                 stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}
Db <- parse_founder(ft$integration,
        if (SIDE_KEY == "A") ft$A_founder_genotype else ft$B_founder_genotype)
cat(sprintf("[3] founder config D_%s (%s): %d defining (barcode,site) positions (expect %d)\n",
            SIDE_KEY, SIDE, nrow(Db), if (SIDE_KEY == "A") 23L else 19L))
stopifnot("|D_b| != expected 23/19" =
            nrow(Db) == (if (SIDE_KEY == "A") 23L else 19L))

# ---- 4. FOUNDER-CONSISTENCY GATE (identical rule to 04_prep_v6_queries.R) --
# recode a barcode column to its 6 site tokens (U->"0", missing/NA/""->"?"),
# then flag a cell if at ANY defining (barcode,site) it observed a NON-"?" token
# that DISAGREES with the founder (divergent edit OR unedited "0").
recode_site <- function(col, k) {
  parts <- tstrsplit(col, "|", fixed = TRUE)
  v <- if (length(parts) >= k) parts[[k]] else rep(NA_character_, length(col))
  v[is.na(v) | v == ""] <- "?"; v[v == "U"] <- "0"; v
}
bad_founder <- rep(FALSE, nrow(dt))
for (r in seq_len(nrow(Db))) {
  g <- Db$barcode[r]; k <- as.integer(sub("Site", "", Db$site[r])); tok <- Db$founder_token[r]
  if (!(g %in% bcs)) next
  val <- recode_site(dt[[g]], k)
  bad_founder <- bad_founder | (val != "?" & val != tok)
}
cat(sprintf("[4] founder gate: %d of the %d QC-band cells are founder-divergent (dropped)\n",
            sum(qc_ok & bad_founder), sum(qc_ok)))

# ---- 5. keep = quality band AND founder-consistent; write filtered TSV -----
keep <- qc_ok & !bad_founder
out  <- dt[keep]
fwrite(out, OUT, sep = "\t", quote = FALSE, na = "NA")   # na="NA" restores the input format
cat(sprintf("[5] wrote %d %s cells -> %s\n", nrow(out), SIDE, OUT))

# ---- 6. fail-fast checkpoint: row count must match the precomputed target --
if (!is.na(EXPECT_KEPT[[SIDE]])) {
  cat(sprintf("[6] survivors = %d ; expected (founder-consistent) = %d\n",
              nrow(out), EXPECT_KEPT[[SIDE]]))
  stopifnot("survivor count != expected founder-consistent target -- founder filter WRONG" =
              nrow(out) == EXPECT_KEPT[[SIDE]])
  cat("[6] OK: founder-filtered survivor count matches counts.tsv.\n")
}
