# 00a_reshape_v8_founder_genotypes.R
# v8's Supplementary Table 1 ships long (one row per integration x blastomere)
# instead of v6's wide (one row per integration, A/B genotype side by side).
# Reshape it to the wide schema so 00_filter_highqual_consensus.R's
# parse_founder() can read it unchanged.
#
# Run from the repo root:
#   Rscript tree_building/1_build_nj_backbone/00a_reshape_v8_founder_genotypes.R

IN  <- Sys.getenv("IN", "support_data/e3v8.supp_table1_founder_genotypes.csv")
OUT <- Sys.getenv("OUT", "tree_building/processed_data/e3v8.supp_table1_founder_genotypes.wide.csv")

# some header lines are quoted ("# ..."), so a leading comment.char doesn't
# catch them -- strip any line starting with # or "# before parsing
lines <- readLines(IN)
lines <- lines[!grepl('^"?#', lines)]
ft <- read.csv(text = lines, stringsAsFactors = FALSE, colClasses = "character")
stopifnot(all(c("integration", "blastomere", "n_defining_sites", "defining_prefix") %in% names(ft)))

# "none"/0 sites -> "(unedited)", matching v6's wide-format sentinel
genotype <- ifelse(ft$defining_prefix == "none", "(unedited)", ft$defining_prefix)

a <- ft[ft$blastomere == "A", c("integration", "defining_prefix", "n_defining_sites")]
b <- ft[ft$blastomere == "B", c("integration", "defining_prefix", "n_defining_sites")]
names(a) <- c("integration", "A_founder_genotype", "A_n_edits")
names(b) <- c("integration", "B_founder_genotype", "B_n_edits")
a$A_founder_genotype[a$A_founder_genotype == "none"] <- "(unedited)"
b$B_founder_genotype[b$B_founder_genotype == "none"] <- "(unedited)"

wide <- merge(a, b, by = "integration", all = TRUE)
stopifnot("integration present in only one blastomere's rows" = !anyNA(wide))

writeLines(c(
  "# Supplementary Table 1 (v8, e3v8): reshaped from e3v8.supp_table1_founder_genotypes.csv",
  "# (long format) into v6's wide format. Genotype content is unchanged (see v8 README.md:",
  "# byte-identical to v6, A=23/B=19 edits).",
  paste0("# total founder edits  A: ", sum(as.integer(wide$A_n_edits)),
         "   B: ", sum(as.integer(wide$B_n_edits)))
), OUT)
write.table(wide, OUT, sep = ",", row.names = FALSE, quote = FALSE, append = TRUE)
cat(sprintf("wrote %d integrations -> %s\n", nrow(wide), OUT))
