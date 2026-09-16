# test dates of tips wrt dating sensitivity analyses
#
#------------------------------------------------------------------------------#
## libs
library(ggplot2)
library(castor)
library(reshape2)
library(wesanderson)

#------------------------------------------------------------------------------#
## figure settings

out_path = "figures/sensitivity_dating/node_time_shifts.pdf"
out_path2 = "figures/sensitivity_dating/node_time_shifts.png"

base_size <- 10

#------------------------------------------------------------------------------#
## Load trees
# load trees
tree_files = c("tree_building/results/3-dated-tree/merged_minB2h_lineage_constrained.nwk",
               "tree_building/results/3-dated-tree-sensitivity/merged/merged_time_tree_ge7_attempt1_minage_sourced_minB2h_l0.01_sidefrac63-37.nwk",
               "tree_building/results/3-dated-tree-sensitivity-5842/merged/merged_time_tree_ge7_attempt1_minage_sourced_minB2h_l0.01_sidefrac58-42.nwk")

trees = lapply(tree_files, function(x) {
  ape::read.tree(x)
})

names(trees) = c("main", "sensitivity_63-37", "sensitivity_58-42")

# The trees have the exact same nodes, but different node dates.
# For each node in the sensitivity trees, quantify the difference in node time
# (increasing from the root) to the main tree node. 

root_depths = lapply(trees, castor::get_all_distances_to_root)

#------------------------------------------------------------------------------#
## Assign blastomere labels

Ntip = length(trees$main$tip.label)
root_id = setdiff(trees$main$edge[,1], trees$main$edge[,2])   # root = node with no parent
children = trees$main$edge[trees$main$edge[,1] == root_id, 2]
stopifnot(length(children) == 2)   # root must be bifurcating for a well-defined A/B split

sub1 = castor::get_subtree_at_node(trees$main, children[1] - Ntip)
sub2 = castor::get_subtree_at_node(trees$main, children[2] - Ntip)

if (length(sub1$new2old_tip) >= length(sub2$new2old_tip)) {
  bigger = sub1; smaller = sub2
} else {
  bigger = sub2; smaller = sub1
}

blastomere = rep(NA_character_, Ntip + trees$main$Nnode)                                      
blastomere[bigger$new2old_tip]         = "A"
blastomere[Ntip + bigger$new2old_node] = "A"
blastomere[smaller$new2old_tip]          = "B"
blastomere[Ntip + smaller$new2old_node]  = "B"

#------------------------------------------------------------------------------#
# compute differences in node times
diff_63_37 = root_depths[["sensitivity_63-37"]] - root_depths$main                          
diff_58_42 = root_depths[["sensitivity_58-42"]] - root_depths$main

diffs = data.frame(original_depth = root_depths$main,
                   blastomere = blastomere,
                   diff_63_37 = diff_63_37,
                   diff_58_42 = diff_58_42)

#exclude tip depths
diffs = diffs[diffs$original_depth < 13.49,]

#exclude root depth because it is fixed and gives NA values to blastomeres
diffs = diffs[diffs$original_depth > 1.49,]

# reshape + rescale to hours
diffs_long <- diffs %>%
  pivot_longer(cols = starts_with("diff_"), names_to = "sensitivity", values_to = "depth_diff") %>%
  mutate(depth_diff_hours = depth_diff * 24)   # days -> hours


diffs_long <- diffs_long %>%
  mutate(depth_bin = cut(original_depth, breaks = breaks, labels = day_labels,
                         include.lowest = TRUE, right = FALSE))

summary_df <- diffs_long %>%
  mutate(depth_bin = as.numeric(as.character(depth_bin))) %>%
  group_by(sensitivity, depth_bin) %>%
  summarise(
    mean_diff_hours = mean(depth_diff_hours),
    se_diff_hours   = sd(depth_diff_hours) / sqrt(n()),
    mean_abs_diff_hours = mean(abs(depth_diff_hours)),
    n = n(),
    .groups = "drop"
  )

cols = c("#7B4FA0", "#1B9AAA") #wesanderson::wes_palette("Zissou1", n = 5, type = "discrete")[c(1,5)]
cols

p = ggplot(diffs_long, aes(x = original_depth, y = depth_diff_hours, color = sensitivity)) +
  facet_wrap(~blastomere)+
  geom_point(size = 0.5, alpha = 0.05) +
  #geom_line(data = summary_df, aes(x = depth_bin, y = mean_diff_hours, color = sensitivity),
  #          linewidth = 1, inherit.aes = FALSE) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-3, 3), linetype = "dashed", color = "black") +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "grey") +
  scale_color_manual(values = cols) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  labs(x = "Original node time [days]", y = "Node time difference [hours]",
       color = "Sensitivity analysis") +
  scale_x_continuous(breaks = seq(1.5, 13.5, by = 1)) +
  scale_y_continuous(breaks = c(-12, -6, -3, -1, 0, 1, 3, 6, 12)) +
  theme_classic(base_size = base_size, base_family = "sans")+
  theme(legend.position = "top") 

p
ggsave(out_path, p, height = 4, width = 8, dpi = 300)
ggsave(out_path2, p, height = 4, width = 8, dpi = 300)
