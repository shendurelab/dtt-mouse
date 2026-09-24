# fig_s12_placement_accuracy_analysis.R
# Fig. S12 A-D: held-out-cell placement accuracy. Drop real backbone cells,
# mask their tape genotype down to a randomly-drawn 4-8 retained tapes
# (simulating incomplete recovery), place them back onto the pruned backbone
# via DTT-ancestor nearest-neighbour search, then score how far (tree edges,
# and days to MRCA) the placement lands from the cell's true original
# position -- vs. the best anchor actually available in the pruned backbone
# (the floor no placement algorithm can beat). 
#
#
# PILOT SCALE: n_reps=1, n_cells=2000 (1000/side), vs. full scale in paper
# n_reps=5, n_cells=100000.
#
# Run from the repo root:
#   Rscript tree_analysis/placement_accuracy/fig_s12_placement_accuracy_analysis.R

library(ape)
library(data.table)
library(castor)

# Define parameters
n_reps = 5
n_cells = 100000
METRIC = "dtt"   # dna typewriter distance, same as in tree reconstruction 

set.seed(1)

SIDES <- list(
  B1 = list(
    tree      = "tree_building/results/3-dated-tree/perside_B1_minB2h_lineage_constrained.nwk",
    consensus = "support_data/e3v8.B1_tape_consensus.tsv.gz"),
  B2 = list(
    tree      = "tree_building/results/3-dated-tree/perside_B2_minB2h_lineage_constrained.nwk",
    consensus = "support_data/e3v8.B2_tape_consensus.tsv.gz"))
for (s in SIDES) stopifnot(file.exists(s$tree), file.exists(s$consensus))

# Load backbone tree
trees <- list()
for (SIDE in names(SIDES)) {
  message("[load] ", SIDE, " ...")
  tr <- read.tree(SIDES[[SIDE]]$tree)
  tip_depths <- node.depth.edgelength(tr)[seq_along(tr$tip.label)]
  message(sprintf("[load] %s: %d tips, tip depth range %.5f - %.5f days",
                  SIDE, Ntip(tr), min(tip_depths), max(tip_depths)))
  trees[[SIDE]] <- tr
}


# Drop cells at random from backbone tree
n_side <- n_cells / length(SIDES)

drop_sets <- list()      # drop_sets[[side]][[rep]] = character vector of dropped cell_ids
pruned_trees <- list()   # pruned_trees[[side]][[rep]] = backbone tree with those tips removed
for (SIDE in names(SIDES)) {
  drop_sets[[SIDE]] <- list(); pruned_trees[[SIDE]] <- list()
  for (rep in seq_len(n_reps)) {
    dropped <- sample(trees[[SIDE]]$tip.label, n_side)
    drop_sets[[SIDE]][[rep]] <- dropped
    pruned_trees[[SIDE]][[rep]] <- drop.tip(trees[[SIDE]], dropped)
    message(sprintf("[drop] %s rep %d: dropped %d -> pruned backbone %d tips",
                    SIDE, rep, length(dropped), Ntip(pruned_trees[[SIDE]][[rep]])))
  }
}


# Subsample the genotype of the cells to 4-6 recovered tapes

# each side's consensus, once; TAPE columns = everything but cell_id + QC cols
QC_COLS <- c("n_loci", "n_doublet_loci", "mean_dominance", "pass_qc")
consensus <- list(); tape_cols <- list()
for (SIDE in names(SIDES)) {
  dt <- fread(cmd = sprintf("zcat < %s", shQuote(SIDES[[SIDE]]$consensus)), colClasses = "character")
  setkey(dt, cell_id)
  consensus[[SIDE]] <- dt
  tape_cols[[SIDE]] <- setdiff(colnames(dt), c("cell_id", QC_COLS))
  message(sprintf("[genotype] %s: %d cells, %d TAPE columns", SIDE, nrow(dt), length(tape_cols[[SIDE]])))
}

# a tape is absent if blank, or "NA" -- as text OR as a real NA (fread converts
# the literal text "NA" to an actual NA, which plain `%in% c("NA")` misses)
is_absent <- function(tape_values) is.na(tape_values) | tape_values %in% c("", "NA")

n_recovered <- function(tape_values) sum(!is_absent(tape_values))

# per-tape non-recovery rate, in tape_cols order -- over the backbone tips only,
# not the whole consensus
tape_dropout_rate <- list()
for (SIDE in names(SIDES)) {
  tips <- trees[[SIDE]]$tip.label
  stopifnot(all(tips %in% consensus[[SIDE]]$cell_id))   # a missing tip would join as an all-NA row
  dt <- consensus[[SIDE]][.(tips)]
  tape_dropout_rate[[SIDE]] <- vapply(tape_cols[[SIDE]],
                                      function(col) mean(is_absent(dt[[col]])), numeric(1))
  message(sprintf("[dropout] %s: non-recovery rate over %d backbone tips: %s", SIDE, nrow(dt),
                  paste(sprintf("%s=%.3f", tape_cols[[SIDE]], tape_dropout_rate[[SIDE]]), collapse = " ")))
}

# blank recovered tapes in one cell's row until only `target` remain, drawing
# them in proportion to how often each tape fails to be recovered population-wide
# (sample() renormalises prob, so the raw rates work as relative weights)
mask_to_target <- function(tape_values, target, dropout_weight) {
  recovered <- which(!is_absent(tape_values))
  n_blank <- length(recovered) - target
  # guard: sample(x, 0, prob=) still validates prob against seq_len(x) when x is scalar
  if (n_blank > 0)
    tape_values[sample(recovered, n_blank, prob = dropout_weight[recovered])] <- ""
  tape_values
}

# retained-tape target for one cell: uniform over 4-8, capped by what the cell
# actually has. Drawn per cell over its own achievable range rather than
# pmin(sample(4:8), present), which would dump every over-ambitious draw onto the
# top reachable bin. 
draw_targets <- function(present) {
  lo <- pmin(4, present); hi <- pmin(8, present)
  lo + floor(runif(length(present)) * (hi - lo + 1))
}

masked_genotypes <- list()   # [[side]][[rep]] = data.table, TAPE cols masked + retained_tapes/present_tapes
for (SIDE in names(SIDES)) {
  cols <- tape_cols[[SIDE]]
  masked_genotypes[[SIDE]] <- list()
  for (rep in seq_len(n_reps)) {
    sub <- copy(consensus[[SIDE]][.(drop_sets[[SIDE]][[rep]])])
    tape_matrix  <- as.matrix(sub[, ..cols])
    present      <- apply(tape_matrix, 1, n_recovered)
    target       <- draw_targets(present)
    masked       <- t(mapply(mask_to_target, asplit(tape_matrix, 1), target,
                             MoreArgs = list(dropout_weight = tape_dropout_rate[[SIDE]])))
    sub[, (cols) := as.data.table(masked)]
    sub[, retained_tapes := target]
    # pre-masking coverage: only cells with present_tapes >= 8 can reach the top
    # retained_tapes bin, so that bin is drawn from better-covered cells
    sub[, present_tapes := present]
    masked_genotypes[[SIDE]][[rep]] <- sub
    message(sprintf("[mask] %s rep %d: %d cells, present_tapes %d-%d (mean %.2f), retained_tapes %d-%d",
                    SIDE, rep, nrow(sub), min(present), max(present), mean(present),
                    min(target), max(target)))
  }
}


# Run placement of the dropped cells on the backbone tree

# dtt-mouse's own placement engine (tree_building/4_full_tree_placement/),
# unmodified; the "anchor universe" it sees is just the pruned backbone +
# this rep's masked cells, not the full production query set
BESTMATCH_PY  <- "tree_building/4_full_tree_placement/bestmatch.py"
DUMP_TSV_PY   <- "tree_analysis/placement_accuracy/dump_bestmatch_tsv.py"
ACC_DIR <- "tree_analysis/placement_accuracy/results"
col_order <- colnames(consensus[["B1"]])   # same TSV column layout for both sides

placements <- list()   # [[side]][[rep]] = data.table(cell_id, matched_id, mb_score, side, rep)

for (SIDE in names(SIDES)) {
  placements[[SIDE]] <- list()
  for (rep in seq_len(n_reps)) {
    # rep_dir: the masked backbone/consensus -- identical regardless of METRIC, so
    # re-running with a different METRIC never needs to rebuild these
    rep_dir <- file.path(ACC_DIR, SIDE, sprintf("rep%d", rep))
    dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)

    backbone_nwk <- file.path(rep_dir, "backbone.nwk")
    write.tree(pruned_trees[[SIDE]][[rep]], backbone_nwk)

    consensus_tsv <- file.path(rep_dir, "consensus.tsv.gz")
    backbone_rows <- consensus[[SIDE]][.(pruned_trees[[SIDE]][[rep]]$tip.label)]
    fwrite(rbindlist(list(backbone_rows, masked_genotypes[[SIDE]][[rep]][, ..col_order])),
           consensus_tsv, sep = "\t", compress = "gzip")

    # out_dir: METRIC-specific bestmatch.py output -- unlinked every run since
    # bestmatch.py reuses cellids_*.npy if it finds one already on disk
    out_dir <- file.path(rep_dir, METRIC)
    unlink(out_dir, recursive = TRUE)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    Sys.setenv(SIDE = SIDE, CONSENSUS_TSV = consensus_tsv, BACKBONE_NWK = backbone_nwk,
               METRIC = METRIC, OUTDIR = out_dir)
    system2("python3", BESTMATCH_PY,
            stdout = file.path(out_dir, "bestmatch.log"), stderr = file.path(out_dir, "bestmatch.log"))

    result_tsv <- file.path(out_dir, "placements.tsv")
    system2("python3", c(DUMP_TSV_PY, out_dir, SIDE, result_tsv))
    placements[[SIDE]][[rep]] <- fread(result_tsv)
    placements[[SIDE]][[rep]][, `:=`(side = SIDE, rep_id = rep, metric = METRIC)]   # NOT `rep` -- would shadow the loop variable inside data.table's `[`
    message(sprintf("[place] %s rep %d (%s): %d placements -> %s",
                    SIDE, rep, METRIC, nrow(placements[[SIDE]][[rep]]), result_tsv))
  }
}


# Score placement in terms of distance of the placed cell to the true location
# of the dropped cell in the backbone tree. Keep track of retained TAPE number,
# edit depth, cell type of that cell.

celltypes <- fread(cmd = "zcat < support_data/cell_metadata.v8.txt.gz",
                   select = c("cell_id", "celltype", "major_trajectory"))

# total edited sites a query carries: summed over its RETAINED (non-blanked)
# tapes, 6 sites each. Scales with tape count by construction, so it tracks
# retained_tapes closely (median 22 edits at 4 tapes, 42 at 8)
total_edit_depth <- function(tape_values) {
  retained <- tape_values[!is_absent(tape_values)]
  edited <- vapply(strsplit(retained, "|", fixed = TRUE),
                   function(s) sum(!(s %in% c("", "NA", "U", "-"))), numeric(1))
  sum(edited)
}

# for every node of tr, a REPRESENTATIVE tip within its clade that's in `survivor_tips`
# (NA if none); postorder guarantees a node's children are resolved before the node
representative_survivor <- function(tr, survivor_tips) {
  n_tip <- Ntip(tr)
  rep_tip <- c(ifelse(tr$tip.label %in% survivor_tips, tr$tip.label, NA_character_),
              rep(NA_character_, tr$Nnode))
  edge <- reorder(tr, "postorder")$edge
  for (i in seq_len(nrow(edge))) {
    p <- edge[i, 1]; ch <- edge[i, 2]
    if (is.na(rep_tip[p]) && !is.na(rep_tip[ch])) rep_tip[p] <- rep_tip[ch]
  }
  rep_tip
}

# a dropped tip's TRUE anchor is a surviving tip from its sibling clade in the
# original tree T; if that whole clade was also dropped, climb to the next
# ancestor and check ITS other children, and so on -- NA only if the tip's
# entire side of the tree (up to the root) was wiped out
true_anchor_tip <- function(tr, dropped_tips, survivor_tips) {
  parent_of <- integer(Ntip(tr) + tr$Nnode)
  parent_of[tr$edge[, 2]] <- tr$edge[, 1]
  children_of <- split(tr$edge[, 2], tr$edge[, 1])
  rep_tip <- representative_survivor(tr, survivor_tips)

  dropped_idx <- match(dropped_tips, tr$tip.label)
  vapply(dropped_idx, function(xi) {
    node <- xi
    repeat {
      p <- parent_of[node]
      if (p == 0) return(NA_character_)              # reached the root, no survivor anywhere
      survs <- na.omit(rep_tip[setdiff(children_of[[as.character(p)]], node)])
      if (length(survs)) return(survs[1])
      node <- p                                       # this whole level was also dropped; climb
    }
  }, character(1))
}

scored <- list()
for (SIDE in names(SIDES)) {
  ntip <- Ntip(trees[[SIDE]])
  tips_per_node <- count_tips_per_node(trees[[SIDE]])   # indexed by node - ntip
  for (rep in seq_len(n_reps)) {
    masked <- masked_genotypes[[SIDE]][[rep]]
    # `..` only resolves a bare symbol -- masked[, ..tape_cols[[SIDE]]] silently
    # returns the column NAMES as data instead of selecting the columns
    cols <- tape_cols[[SIDE]]
    total_edits <- apply(as.matrix(masked[, ..cols]), 1, total_edit_depth)
    stopifnot(length(total_edits) == nrow(masked))
    anchor_tip <- true_anchor_tip(trees[[SIDE]], drop_sets[[SIDE]][[rep]],
                                  pruned_trees[[SIDE]][[rep]]$tip.label)
    stopifnot(all(na.omit(anchor_tip) %in% pruned_trees[[SIDE]][[rep]]$tip.label))

    sc <- copy(placements[[SIDE]][[rep]])   # copy() -- data.table's := would otherwise mutate `placements` too
    sc[, retained_tapes := masked$retained_tapes[match(cell_id, masked$cell_id)]]
    sc[, present_tapes := masked$present_tapes[match(cell_id, masked$cell_id)]]
    sc[, total_edits := total_edits[match(cell_id, masked$cell_id)]]
    sc[, celltype := celltypes$celltype[match(cell_id, celltypes$cell_id)]]
    sc[, major_trajectory := celltypes$major_trajectory[match(cell_id, celltypes$cell_id)]]
    sc[, true_anchor_tip := anchor_tip[match(cell_id, drop_sets[[SIDE]][[rep]])]]
    sc[, true_anchor_survived := !is.na(true_anchor_tip)]

    # how far the true anchor itself sits from the query -- a distant witness (deep
    # in the sibling clade, close relatives all dropped) carries little genotype
    # signal, so placement can't be expected to land close even if it "succeeds"
    sc[true_anchor_survived == TRUE, true_anchor_dist := get_pairwise_distances(
      trees[[SIDE]], A = cell_id, B = true_anchor_tip, as_edge_counts = TRUE)]

    # true topological distance, IN THE REAL TREE, between the dropped cell's
    # actual tip and the tip it was placed next to (NA if unplaceable)
    sc[matched_id != "", node_dist := get_pairwise_distances(
      trees[[SIDE]], A = cell_id, B = matched_id, as_edge_counts = TRUE)]

    # same, but on the PRUNED backbone (using true_anchor_tip -- cell_id itself isn't
    # a tip there) -- drop.tip already collapses every edge whose other branch was
    # entirely dropped cells, so this is node_dist with exactly the invisible-to-the-
    # algorithm edges removed. Undefined if either endpoint is missing.
    sc[true_anchor_survived == TRUE & matched_id != "", node_dist_visible := get_pairwise_distances(
      pruned_trees[[SIDE]][[rep]], A = true_anchor_tip, B = matched_id, as_edge_counts = TRUE)]

    # same idea in TIME (days) rather than edges: total developmental time both
    # lineages spent apart since their MRCA (tree is dated/ultrametric, so this
    # is 2x how long ago the split happened) -- a placement can be many edges
    # away yet recently diverged, or few edges away but split very early
    sc[matched_id != "", mrca_divergence_days := get_pairwise_distances(
      trees[[SIDE]], A = cell_id, B = matched_id, as_edge_counts = FALSE)]
    sc[true_anchor_survived == TRUE, true_anchor_divergence_days := get_pairwise_distances(
      trees[[SIDE]], A = cell_id, B = true_anchor_tip, as_edge_counts = FALSE)]

    # tips under MRCA(query, match): the set of cells the placement can't tell the
    # query apart from. Reads node_dist as a neighbourhood size without assuming
    # the tree is balanced -- a d-edge error spans 2^(d/2) tips only if it is.
    sc[matched_id != "", mrca_clade_size := tips_per_node[get_pairwise_mrcas(
      trees[[SIDE]], A = cell_id, B = matched_id, check_input = FALSE) - ntip]]
    sc[true_anchor_survived == TRUE, true_anchor_clade_size := tips_per_node[get_pairwise_mrcas(
      trees[[SIDE]], A = cell_id, B = true_anchor_tip, check_input = FALSE) - ntip]]

    scored[[length(scored) + 1]] <- sc
    message(sprintf("[score] %s rep %d: median MRCA clade %.0f tips (anchor floor %.0f)",
                    SIDE, rep, median(sc$mrca_clade_size, na.rm = TRUE),
                    median(sc$true_anchor_clade_size, na.rm = TRUE)))
    message(sprintf("[score] %s rep %d: %d/%d placed, median node_dist %.1f (%.2f days), exact %.1f%%, true anchor survived %.1f%% (median dist to it %.1f, %.2f days)",
                    SIDE, rep, sum(sc$matched_id != ""), nrow(sc),
                    median(sc$node_dist, na.rm = TRUE), median(sc$mrca_divergence_days, na.rm = TRUE),
                    100 * mean(sc$node_dist == 0, na.rm = TRUE),
                    100 * mean(sc$true_anchor_survived), median(sc$true_anchor_dist, na.rm = TRUE),
                    median(sc$true_anchor_divergence_days, na.rm = TRUE)))
  }
}
scored <- rbindlist(scored)

# save both artifacts, tagged by METRIC
placements_flat <- rbindlist(unlist(placements, recursive = FALSE))
fwrite(placements_flat, file.path(ACC_DIR, sprintf("placements_%s.tsv", METRIC)), sep = "\t")
fwrite(scored, file.path(ACC_DIR, sprintf("scored_%s.tsv", METRIC)), sep = "\t")
message(sprintf("[save] wrote placements_%s.tsv and scored_%s.tsv to %s", METRIC, METRIC, ACC_DIR))

# staged copy for making_figures/Figure-S12.R
FIG_CSV <- "figures_data/fig_s12_placement_accuracy_scored_dtt.csv"
fwrite(scored, FIG_CSV)
message("wrote ", FIG_CSV)
