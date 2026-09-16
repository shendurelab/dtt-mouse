

#######################################################################
### Figure 5. Lineage predicts fate at every scale, from sibling pairs 
### to a dated hierarchy of cell-type couplings


#############################################################################
### Fig. 5A: Heatmap of the 37 progenitor & post-mitotic heterotypic pairings

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/figures/v8/heterotypic"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/fig6B_captured_divisions.csv"))

plot_dat <- dat %>%
    mutate(log2_fold = log2(fold))

prog_levels <- c(
    "Border-associated macrophages",
    "Definitive early erythroblasts (CD36-)",
    "Hematopoietic stem cells (Cd34+)",
    "Hematopoietic stem cells (Mpo+)",
    "Kupffer cells",
    "Myelinating Schwann cells",
    "Myelinating Schwann cells (Tgfb2+)",
    "Neural crest (PNS glia)",
    "Astrocytes",
    "Dorsal telencephalon",
    "Eye field",
    "Hindbrain",
    "Hypothalamus",
    "Intermediate neuronal progenitors",
    "Naive retinal progenitor cells",
    "Olfactory bulb cells",
    "Retinal progenitor cells",
    "Spinal cord dorsal progenitors",
    "Spinal cord/r7/r8",
    "Telencephalon",
    "Branchial arch epithelium",
    "Olfactory epithelial cells",
    "Otic epithelial cells",
    "Pre-epidermal keratinocytes"
)

post_levels <- c(
    "Granulocytes",
    "Megakaryocytes",
    "Osteoclasts",
    "Primitive erythroid cells",
    "Neural crest (PNS neurons)",
    "Parasympathetic neurons",
    "Sympathetic neurons",
    "Cajal-Retzius cells",
    "Cerebellar Purkinje cells",
    "Cranial motor neurons",
    "Deep-layer neurons",
    "GABAergic cortical interneurons",
    "GABAergic neurons",
    "Retinal ganglion cells",
    "Spinal cord motor neurons",
    "Granular keratinocytes",
    "Olfactory sensory neurons",
    "Otic sensory neurons"
)

plot_dat <- plot_dat %>%
    mutate(progenitor = factor(progenitor, levels = rev(prog_levels)),  # rev so first is on top
           postmitotic = factor(postmitotic, levels = post_levels))

p = ggplot(plot_dat, aes(x = postmitotic, y = progenitor, fill = log2_fold)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = n_cells, color = log2_fold > 3),
              size = 1.5) +
    scale_fill_viridis_c(option = "magma", name = expression(log[2]~fold~enrichment)) +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = "none") +
    coord_fixed() +
    labs(x = "post-mitotic derivative",
         y = "progenitor / regional state",
         title = "Captured differentiation divisions") +
    theme_minimal(base_size = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))

ggsave(paste0(save_path, "/Fig5/Fig5A.pdf"), p, height = 9, width = 4.5)


####################################################
### Fig. 5C: Timed fate couplings between cell types

library(ggplot2)
library(dplyr)
library(forcats)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/figures/v8/coupling_depth"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/figXa_coupling_carpet.csv"))

plot_dat <- dat |>
    mutate(
        pair = paste(celltype_1, celltype_2, sep = "  |  "),
        # row_in_figure already encodes the coupling-depth ordering (1 = top)
        pair = fct_reorder(pair, row_in_figure, .desc = TRUE)
    )

# anchor the diverging palette at 0 so black always means log2 enr = 0,
# regardless of the asymmetric range (~ -1.5 to +2.5)
lim  <- max(abs(range(plot_dat$log2_enrichment_mean, na.rm = TRUE)))
brks <- c(-lim, -lim/2, 0, lim/2, lim)
pal  <- c("#00E5FF", "#1E5FCC", "#000000", "#C0661A", "#FFD400")

p = ggplot(plot_dat, aes(clade_ancestor_age_E, pair, fill = log2_enrichment_mean)) +
    geom_raster() +
    geom_point(
        data = \(d) filter(d, significant_both_blastomeres == 1),
        colour = "white", size = 0.35, stroke = 0, shape = 16
    ) +
    scale_fill_gradientn(
        colours = pal,
        values  = scales::rescale(brks),
        limits  = c(-lim, lim),
        name    = expression(log[2]~enrichment~(mean~of~A*","~B))
    ) +
    scale_x_continuous(
        breaks = 7:13, labels = paste0("E", 7:13), expand = c(0, 0)
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
        title    = "Cell-type coupling across developmental time  (n = 129 pairs)",
        subtitle = "clades defined by the age of their dated ancestor; rows ordered by coupling depth; white dots = significant in both blastomeres",
        x = "developmental time of the clade ancestor  (earlier to the left)",
        y = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
        axis.text.y        = element_text(size = 5.5),
        axis.ticks.y       = element_blank(),
        panel.grid         = element_blank(),
        plot.title.position = "plot",
        legend.key.height  = unit(1.2, "cm"),
        legend.key.width   = unit(0.35, "cm")
    )

ggsave(paste0(save_path, "/Fig5/Fig5C.pdf"), p, height = 9, width = 9)

write.table(unique(plot_dat$pair), paste0(save_path, "/Fig5/Fig5C_rownames.txt"), row.names=F, col.names=F, sep="\t", quote=F)


######################################################
### Fig. 5D: A draft hierarchy of timed fate couplings

library(ggplot2)
library(ggdendro)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/figures/v8/coupling_depth"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

INK  <- "#1f2328"
INK2 <- "#4a5158"
EDGE <- "#7f868d"
DATE <- "#8a3324"

# ---- 1. Read data --------------------------------------------------------
depth_mat <- read.csv(file.path(data_path, "figXb_coupling_depth_matrix.csv"),
                      row.names = 1, check.names = FALSE)
labels   <- rownames(depth_mat)
n_leaves <- length(labels)

hier <- read.csv(file.path(data_path, "figXb_coupling_hierarchy.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
hier <- hier[order(hier$merge_order), ]
hier$member_list <- strsplit(hier$members, "; ", fixed = TRUE)

# ---- 2. Reconstruct hclust merge matrix from member lists ---------------
# clusters: named list; name is cluster id (-k for leaf k, +i for merge i)
clusters <- setNames(as.list(labels), as.character(-seq_len(n_leaves)))

merge_left   <- integer(nrow(hier))
merge_right  <- integer(nrow(hier))
merge_height <- numeric(nrow(hier))

for (i in seq_len(nrow(hier))) {
    target <- hier$member_list[[i]]
    cluster_ids <- names(clusters)
    K <- length(clusters)
    found <- NULL
    for (a in seq_len(K - 1)) {
        for (b in seq(a + 1, K)) {
            u <- c(clusters[[a]], clusters[[b]])
            if (length(u) == length(target) && setequal(u, target)) {
                found <- c(a, b); break
            }
        }
        if (!is.null(found)) break
    }
    stopifnot(!is.null(found))
    merge_left[i]   <- as.integer(cluster_ids[found[1]])
    merge_right[i]  <- as.integer(cluster_ids[found[2]])
    merge_height[i] <- hier$branch_point_age_E[i]
    clusters <- clusters[-found]
    clusters[[as.character(i)]] <- target
}
stopifnot(!is.unsorted(merge_height))

# ---- 3. Derive a valid leaf order by post-order traversal ---------------
merge_mat <- cbind(merge_left, merge_right)
get_leaves <- function(id) {
    if (id < 0) return(-id)
    c(get_leaves(merge_mat[id, 1]), get_leaves(merge_mat[id, 2]))
}
leaf_order <- get_leaves(nrow(merge_mat))

hc <- structure(
    list(merge = merge_mat, height = merge_height, order = leaf_order,
         labels = labels, method = "average", call = match.call(),
         dist.method = "user"),
    class = "hclust"
)

# ---- 4. Rank-transform heights (ties at 13.5 would otherwise collapse) ---
real_heights   <- hc$height
hc_rank        <- hc
hc_rank$height <- rank(hc$height, ties.method = "first")

dd  <- dendro_data(as.dendrogram(hc_rank), type = "rectangle")
seg <- dd$segments
lab <- dd$labels
lab$label <- as.character(lab$label)

seg <- transform(seg,
                 plot_x    = -y,    plot_y    = x,
                 plot_xend = -yend, plot_yend = xend)

# ---- 5. Branch-point date labels ----------------------------------------
top_segs <- subset(seg, y == yend & x != xend)
tops <- do.call(rbind, lapply(split(top_segs, top_segs$y), function(g) {
    data.frame(rank     = round(g$y[1]),
               leaf_mid = mean(range(c(g$x, g$xend))),
               stringsAsFactors = FALSE)
}))
tops$age <- real_heights[tops$rank]
tops     <- subset(tops, age <= 13.25)     # 51 branch points

# ---- 6. Group brackets (right side) -------------------------------------
groups <- list(
    "Otic"               = c("Otic sensory neurons", "Otic epithelial cells"),
    "Gut tube"           = c("Pancreatic acinar cells", "Midgut/Hindgut epithelial cells"),
    "Keratinocytes"      = c("Pre-epidermal keratinocytes", "Granular keratinocytes",
                             "Branchial arch epithelium", "Lung progenitor cells"),
    "Olfactory"          = c("Olfactory epithelial cells", "Olfactory bulb cells",
                             "Olfactory sensory neurons"),
    "Haematopoiesis"     = c("Microglia", "Definitive erythroblasts (CD36+)", "Mast cells",
                             "Monocytes", "Hematopoietic stem cells (Mpo+)", "Kupffer cells",
                             "Megakaryocytes", "Border-associated macrophages",
                             "Hematopoietic stem cells (Cd34+)", "Primitive erythroid cells"),
    "Retina / eye field" = c("Retinal progenitor cells", "Retinal pigment cells", "Eye field",
                             "Retinal ganglion cells", "Naive retinal progenitor cells"),
    "Neural crest / PNS" = c("Neural crest (PNS neurons)", "Dorsal root ganglion neurons",
                             "Neural crest (PNS glia)", "Myelinating Schwann cells",
                             "Sympathetic neurons", "Parasympathetic neurons",
                             "Myelinating Schwann cells (Tgfb2+)", "Cranial motor neurons",
                             "Melanocyte cells"),
    "Kidney"             = c("Nephron progenitors", "Metanephric mesenchyme", "Ureteric bud")
)

grp_df <- do.call(rbind, lapply(names(groups), function(nm) {
    ys <- lab$x[match(groups[[nm]], lab$label)]
    ys <- ys[!is.na(ys)]
    if (!length(ys)) return(NULL)
    data.frame(name = nm, y_pos = mean(range(ys)),          # midpoint, not mean
               y_min = min(ys), y_max = max(ys),
               stringsAsFactors = FALSE)
}))

# ---- 7. Draw -------------------------------------------------------------
depth_max    <- max(seg$y)
label_x      <- 0.4                    # leaf label column (right side)
group_text_x <- -(depth_max + 6)       # group label column (left of root)

p <- ggplot() +
    geom_segment(data = seg,
                 aes(x = plot_x, y = plot_y, xend = plot_xend, yend = plot_yend),
                 colour = EDGE, linewidth = 0.32, lineend = "round") +
    
    geom_text(data = tops,
              aes(x = -rank, y = leaf_mid, label = sprintf("E%.2f", age)),
              hjust = 1, vjust = -0.4, size = 2.0, colour = DATE,
              nudge_x = -0.15) +
    
    geom_text(data = lab,
              aes(x = label_x, y = x, label = label),
              hjust = 0, size = 2.2, colour = INK) +
    
    # left-side group labels (no brackets)
    geom_text(data = grp_df,
              aes(x = group_text_x, y = y_pos, label = name),
              hjust = 0, size = 2.5, colour = INK, fontface = "italic") +
    
    scale_x_continuous(expand = expansion(add = c(8, 10))) +
    scale_y_continuous(expand = expansion(add = c(0.4, 0.4))) +
    coord_cartesian(clip = "off") +
    
    labs(
        title    = "Cell types grouped by the depth at which their lineage coupling remains detectable",
        subtitle = "dates are coupling depths in E-days, not divergence times; groups closing earliest stay co-restricted deepest",
        caption  = "branch points labelled with the developmental time at which the group's coupling is still detectable (horizontal distance not to scale)",
        x = NULL, y = NULL
    ) +
    theme_void(base_size = 9) +
    theme(
        plot.title    = element_text(size = 11,  hjust = 0.5, colour = INK,  margin = margin(b = 2)),
        plot.subtitle = element_text(size = 9.5, hjust = 0.5, colour = INK2, margin = margin(b = 12)),
        plot.caption  = element_text(size = 9,   hjust = 0.5, colour = INK2, margin = margin(t = 10)),
        plot.margin   = margin(12, 14, 12, 14)
    )

ggsave(paste0(save_path, "/Fig5/Fig5D.pdf"), p, height = 10, width = 4)

