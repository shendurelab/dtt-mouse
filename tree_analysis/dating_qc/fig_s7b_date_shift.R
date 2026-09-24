#!/usr/bin/env Rscript
# fig_s7b_date_shift.R
# Fig. S7B: per-internal-node date shift under constrained re-dating, vs. the
# node's unconstrained date. 
#
# node date (E) = node.depth.edgelength() on a merged time tree whose root is
#                 the day-0 zygote and whose tips are E13.5.
# shift (days)  = new (constrained) node date - old (unconstrained) node date.
#
#
# Run from the repo root:
#   Rscript tree_analysis/dating_qc/fig_s7b_date_shift.R

suppressPackageStartupMessages({ library(ape); library(castor) })
source("tree_analysis/ancestral_state/01_parsimony_lib.R")   # .get_root()/.build_kids()/.internal_postorder()

UNC <- "tree_building/results/3-dated-tree/merged_minB2h_unconstrained.nwk"
CON <- "tree_building/results/3-dated-tree/merged_minB2h_lineage_constrained.nwk"
stopifnot(file.exists(UNC), file.exists(CON))
OUT_CSV <- "figures_data/fig_s7_dating_panel_C_date_shift.csv"

u <- read.tree(UNC); c <- read.tree(CON)
stopifnot(identical(u$tip.label, c$tip.label))

# ---- match every internal node of c to its counterpart in u ---------------
# Representative-tip trick: for every internal node of c, pick one descendant
# tip per child (from c's own structure only). The MRCA in u of two children's
# representative tips IS the node's counterpart in u, regardless of which tip
# is picked -- so every node's match can be queried in one batched call.
ntip      <- Ntip(u)
n_total_c <- ntip + c$Nnode

root_c <- .get_root(c)
kids_c <- .build_kids(c, n_total_c)
post_c <- .internal_postorder(c, root_c, ntip)

rep_tip_c <- integer(n_total_c)
rep_tip_c[seq_len(ntip)] <- seq_len(ntip)   # tips: same label order in both trees
for (node in post_c) rep_tip_c[node] <- rep_tip_c[kids_c[[node]][1]]

# Binary nodes (near-universal) go through the one batched call; polytomy
# nodes (rare -- LSD2 collapsing a handful of near-zero branches) get one
# castor::get_mrca_of_set() call each.
is_binary    <- lengths(kids_c[post_c]) == 2L
binary_nodes <- post_c[is_binary]
poly_nodes   <- post_c[!is_binary]

mrca_id_u <- integer(n_total_c)
mrca_id_u[seq_len(ntip)] <- seq_len(ntip)

A <- vapply(kids_c[binary_nodes], function(k) rep_tip_c[k[1]], integer(1))
B <- vapply(kids_c[binary_nodes], function(k) rep_tip_c[k[2]], integer(1))
mrca_id_u[binary_nodes] <- get_pairwise_mrcas(u, A, B)
for (node in poly_nodes) mrca_id_u[node] <- get_mrca_of_set(u, rep_tip_c[kids_c[[node]]])

cat(sprintf("node matching: %d binary (1 batched castor call), %d polytomy (individual calls)\n",
            length(binary_nodes), length(poly_nodes)))

# ---- per-internal-node old / new / shift -----------------------------------
ir    <- (ntip + 1):(ntip + c$Nnode)
old   <- node.depth.edgelength(u)[mrca_id_u[ir]]  # E-day, unconstrained counterpart
new   <- node.depth.edgelength(c)[ir]              # E-day, constrained
shift <- new - old

# ---- blastomere: larger root-child clade = B1, smaller = B2 ----------------
# Same "larger = A/B1" convention as fig_2jk_analysis.R; labels kept as B1/B2
# here (not relabelled to Blastomere A/B) since that's what Figure-S7.R's
# Fig. S7B section already expects in this column.
suppressPackageStartupMessages(library(phangorn))
root <- setdiff(u$edge[, 1], u$edge[, 2])
kids <- u$edge[u$edge[, 1] == root, 2]
stopifnot(length(kids) == 2)
csize <- vapply(kids, function(k) length(phangorn::Descendants(u, k, type = "tips")[[1]]), integer(1))
kids <- kids[order(-csize)]   # [B1_node, B2_node]

desc <- phangorn::Descendants(u, kids, type = "all")
is_b1 <- logical(ntip + u$Nnode); is_b1[c(kids[1], desc[[1]])] <- TRUE
is_b2 <- logical(ntip + u$Nnode); is_b2[c(kids[2], desc[[2]])] <- TRUE
# is_b1/is_b2 are sized/indexed for u's own node space; ir indexes c's node
# space, so look membership up via each c-node's matched u-node (mrca_id_u).
blast <- ifelse(is_b1[mrca_id_u[ir]], "B1", ifelse(is_b2[mrca_id_u[ir]], "B2", NA_character_))

d <- data.frame(old = old, shift = shift, blast = factor(blast, levels = c("B1", "B2")))
d <- d[!is.na(d$blast), ]   # drops only the root/stem node(s) under neither clade
cat(sprintf("internal nodes: %d (B1=%d, B2=%d)   shift range: [%+.3f, %+.3f] d\n",
            nrow(d), sum(d$blast == "B1"), sum(d$blast == "B2"), min(shift), max(shift)))

write.csv(d, OUT_CSV, row.names = FALSE)
cat("wrote", OUT_CSV, "\n")
