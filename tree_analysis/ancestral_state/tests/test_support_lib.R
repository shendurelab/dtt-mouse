#!/usr/bin/env Rscript
# test_support_lib.R
# Oracle test for 01_parsimony_lib.R + 02_support_lib.R, in the same
# stopifnot-oracle style used by this repo's other test_*.R files.
#
# Part 1 replays toy_example/ -- a small hand-worked tree+tape example --
# through the REAL production path (reconstruct_support()) and asserts every
# node's per-tape edit count and total support against the expected values.
# Part 2 checks the missing-tape ("?") conservative-handling edge cases
# directly against pars_ancestral_state()/count_edits().
#
# Run from tree_analysis/:
#   Rscript ancestral_state/tests/test_support_lib.R

source("ancestral_state/01_parsimony_lib.R")
source("ancestral_state/02_support_lib.R")
source("ancestral_state/03_tape_tip_state_mapping.R")

pass <- function(msg) cat("PASS:", msg, "\n")

# ---- Part 1: toy_example/ end-to-end, matches the hand-worked oracle -------

tree <- ape::read.tree("ancestral_state/tests/toy_example/toy_tree.nwk")
tip_states_long <- read.csv("ancestral_state/tests/toy_example/toy_tip_states.csv",
                             colClasses = "character")
# tape2's Site2 is always "0" (inert padding, not part of the traced
# example): support.jpeg's tape2 is single-site, but reconstruct_support()
# takes one site_cols list for every tape in a call, so tape2 gets a second,
# always-unedited column just so both tapes can share that one argument.

result <- reconstruct_support(tree, tip_states_long,
                               barcodes = c("tape1", "tape2"),
                               site_cols = c("Site1", "Site2"),
                               mc.cores = 1)
support_df <- result$support_df

node_of <- function(label) {
  if (label %in% tree$tip.label) return(match(label, tree$tip.label))
  ape::getMRCA(tree, strsplit(label, "\\+")[[1]])
}
row_for <- function(label) support_df[support_df$node == node_of(label), ]

check_row <- function(label, exp_tape1, exp_tape2, exp_support) {
  r <- row_for(label)
  ok <- identical(unname(r$tape1), exp_tape1) &&
        identical(unname(r$tape2), exp_tape2) &&
        identical(unname(r$support), exp_support)
  cat(sprintf("  [%s] node %-8s tape1=%s tape2=%s support=%s (expected %s,%s,%s)\n",
              ifelse(ok, "OK", "FAIL"), label, r$tape1, r$tape2, r$support,
              exp_tape1, exp_tape2, exp_support))
  if (!ok) stop("support.jpeg oracle mismatch at node ", label)
}

check_row("t1",       0, 1, 1)
check_row("t2",       0, 1, 1)
check_row("t3",       1, 1, 2)
check_row("t1+t2",    1, 0, 1)   # MRCA(t1,t2) -- "node 4" in support.jpeg
check_row("t1+t2+t3", 1, 0, 1)   # root -- "node 5" in support.jpeg

pass("toy_example reproduces support.jpeg's support column (1,1,2,1,1) exactly")

# ---- Part 2: missing-tape ("?") conservative handling -----------------------

# One child entirely missing: the known sibling's state passes through
# unchanged, and BOTH children get 0 edits charged on this branch (we don't
# know whether the known state's edits happened above or below this split).
missing_one <- pars_ancestral_state(c("AB", "0"), c("?", "?"))
stopifnot(identical(missing_one, c("AB", "0")))
stopifnot(identical(count_edits(missing_one, c("AB", "0")), 0L))
stopifnot(identical(count_edits(missing_one, c("?", "?")), 0L))
pass("one child all-\"?\": sibling passes through, 0 edits charged either side")

# Both children entirely missing: parent is deferred upward as all-"?", not
# forced to "0" -- never manufacture an edit out of pure ambiguity.
missing_both <- pars_ancestral_state(c("?", "?"), c("?", "?"))
stopifnot(identical(missing_both, c("?", "?")))
stopifnot(identical(count_edits(rep("0", 2), missing_both), 0L))
pass("both children all-\"?\": parent stays \"?\", contributes 0 edits upward")


# ---- Part 3: sentinel recoding (03_tape_tip_state_mapping.R) ---------------

raw <- data.frame(cell_id = c("c1", "c2"), barcode = c("bc1", "bc1"),
                   Site1 = c("ETY", "None"), Site2 = c("AAG", "None"),
                   stringsAsFactors = FALSE)
recoded <- recode_tape_tokens(raw, site_cols = c("Site1", "Site2"))
stopifnot(identical(recoded$Site1, c("0", "?")))
stopifnot(identical(recoded$Site2, c("AAG", "?")))
pass("recode_tape_tokens(): ETY -> \"0\", None -> \"?\", edit tokens kept verbatim")

cat("\nAll support_lib.R tests passed.\n")
