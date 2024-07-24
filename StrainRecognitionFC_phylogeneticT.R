# 1. Clear WS, set WD, load libraries ----

rm(list = ls())
setwd("/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM")

library(ape)
library(phangorn)
library(DECIPHER)
library(ggplot2)
library(ggtree)
library(cowplot)
library(tidyverse)
library(gridExtra)
library(grid)
library(gtable)

seed <- 777
set.seed(seed)

source(file = "/Projects1/Fabian/paper_theme_fab.R")


# 2. Load data ----

sequences <- readDNAStringSet("/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/sequences_14strain.fasta")


# 3. Construct phylogenetic tree ----

aligned_sequences <- DECIPHER::AlignSeqs(DNAStringSet(sequences)) # align sequences

aligned_sequences_ape <- as.phyDat(as(aligned_sequences, "matrix")) # Convert alignment to suitable format for ape

## Construct tree

# Calculate distance matrix
distance_matrix <- dist.ml(aligned_sequences_ape)

# Construct tree using neighbor-joining
nj_tree <- nj(distance_matrix)

# Construct tree using maximum likelihood
fit_ml <- pml(nj_tree, data = aligned_sequences_ape)
fit_ml <- optim.pml(fit_ml, model = "HKY")

ml_tree <- fit_ml$tree


# 4. Visualization tree ----

plot(nj_tree)
plot(ml_tree)

plot_nj_tree <- ggtree(nj_tree, layout = "roundrect", linewidth = 0.8) +
  geom_tiplab(size = 6, fontface = "italic") +
  xlim_tree(xlim = c(0, 0.85))+
  theme_tree()
print(plot_nj_tree)

plot_ml_tree <- ggtree(ml_tree) +
  geom_tiplab() +
  theme_tree()
print(plot_ml_tree)
  

# 5. Calculate mean distance for each mock/co-culture ----

distance_mx <- as.matrix(distance_matrix)
distance_df_full <- as.data.frame(distance_mx)
distance_mx[upper.tri(distance_mx)] <- NA
distance_df <- as.data.frame(distance_mx)
distance_df[distance_df==0] <- NA

# Mock 1: So, Fn
distance_m1 <- distance_df[c(4, 11), c(4, 11)]
mean_m1 <- mean(unlist(distance_m1), na.rm = TRUE)
sd_m1 <- sd(unlist(distance_m1), na.rm = TRUE)

# Mock 2: So, Fn, Pg
distance_m2 <- distance_df[c(4, 5, 11), c(4, 5, 11)]
mean_m2 <- mean(unlist(distance_m2), na.rm = TRUE)
sd_m2 <- sd(unlist(distance_m2), na.rm = TRUE)

# Mock 3: So, Fn, Pg, Vp
distance_m3 <- distance_df[c(4, 5, 7, 11), c(4, 5, 7, 11)]
mean_m3 <- mean(unlist(distance_m3), na.rm = TRUE)
sd_m3 <- sd(unlist(distance_m3), na.rm = TRUE)

# Mock 4: So, Ssal, Ssan, Smi, Sg, Vp, Av, An
distance_m4 <- distance_df[c(2, 3, 7, 8, 9, 11, 12, 13), c(2, 3, 7, 8, 9, 11, 12, 13)]
mean_m4 <- mean(unlist(distance_m4), na.rm = TRUE)
sd_m4 <- sd(unlist(distance_m4), na.rm = TRUE)

# Mock 5: So, Ssal, Ssan, Smi, Sg, Smu, Ssob
distance_m5 <- distance_df[c(8:14), c(8:14)]
mean_m5 <- mean(unlist(distance_m5), na.rm = TRUE)
sd_m5 <- sd(unlist(distance_m5), na.rm = TRUE)

# Mock 6: Aa, Fn, Pg, Pi, Smu, Ssob
distance_m6 <- distance_df[c(1, 4, 5, 6, 10, 14), c(1, 4, 5, 6, 10, 14)]
mean_m6 <- mean(unlist(distance_m6), na.rm = TRUE)
sd_m6 <- sd(unlist(distance_m6), na.rm = TRUE)

# Mock 7: So, Ssal, Ssan, Smi, Sg, Vp, Av, An, Aa, Fn, Pg, Pi, Smu, Ssob
distance_m7 <- distance_df
mean_m7 <- mean(unlist(distance_m7), na.rm = TRUE)
sd_m7 <- sd(unlist(distance_m7), na.rm = TRUE)

# Join data in data frame
mean_distance_mocks <- data.frame(mock = c("Mock 1", "Mock 2", "Mock 3", "Mock 4", "Mock 5", "Mock 6", "Mock 7", "Mock 8", "Mock 9"),
                                  mean_dist = c(mean_m1, mean_m2, mean_m3, mean_m4, mean_m5, mean_m6, mean_m7, mean_m1, mean_m1),
                                  sd_dist = c(sd_m1, sd_m2, sd_m3, sd_m4, sd_m5, sd_m6, sd_m7, sd_m1, sd_m1))

plot_mean_dist_mocks <- ggplot(data = mean_distance_mocks, aes(x = mock, y = mean_dist)) +
  geom_point(size = 7) +
  geom_errorbar(aes(ymin = mean_dist-sd_dist, ymax = mean_dist+sd_dist), width = 0.1, size = 1) +
  labs(x = NULL, y = "Mean phylogenetic distance") +
  paper_theme_fab +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(limits = c(0, 0.42), breaks = seq(0, 0.4, by = 0.1)) +
  scale_x_discrete(labels = c("Mock 1" = "Mock 1/\nCo-culture 1", "Mock 2" = "Mock 2/\nCo-culture 2", "Mock 3" = "Mock 3/\nCo-culture 3", "Mock 4" = "Mock 4", "Mock 5" = "Mock 5", "Mock 6" = "Mock 6", "Mock 7" = "Mock 7", "Mock 8" = "Mock 8", "Mock 9" = "Mock 9"))
print(plot_mean_dist_mocks)

# Create boxplot
distance_m1_flat <- unlist(distance_m1)
distance_m1_flat <- distance_m1_flat[!is.na(distance_m1_flat)]
distance_m1_flat <- unname(distance_m1_flat)

distance_m2_flat <- unlist(distance_m2)
distance_m2_flat <- distance_m2_flat[!is.na(distance_m2_flat)]
distance_m2_flat <- unname(distance_m2_flat)

distance_m3_flat <- unlist(distance_m3)
distance_m3_flat <- distance_m3_flat[!is.na(distance_m3_flat)]
distance_m3_flat <- unname(distance_m3_flat)

distance_m4_flat <- unlist(distance_m4)
distance_m4_flat <- distance_m4_flat[!is.na(distance_m4_flat)]
distance_m4_flat <- unname(distance_m4_flat)

distance_m5_flat <- unlist(distance_m5)
distance_m5_flat <- distance_m5_flat[!is.na(distance_m5_flat)]
distance_m5_flat <- unname(distance_m5_flat)

distance_m6_flat <- unlist(distance_m6)
distance_m6_flat <- distance_m6_flat[!is.na(distance_m6_flat)]
distance_m6_flat <- unname(distance_m6_flat)

distance_m7_flat <- unlist(distance_m7)
distance_m7_flat <- distance_m7_flat[!is.na(distance_m7_flat)]
distance_m7_flat <- unname(distance_m7_flat)

distance_box <- data.frame(mock = c(rep("Mock 1", length(distance_m1_flat)),
                                    rep("Mock 2", length(distance_m2_flat)),
                                    rep("Mock 3", length(distance_m3_flat)),
                                    rep("Mock 4", length(distance_m4_flat)),
                                    rep("Mock 5", length(distance_m5_flat)),
                                    rep("Mock 6", length(distance_m6_flat)),
                                    rep("Mock 7", length(distance_m7_flat)),
                                    rep("Mock 8", length(distance_m1_flat)),
                                    rep("Mock 9", length(distance_m1_flat))),
                           distance = c(distance_m1_flat, distance_m2_flat, distance_m3_flat, distance_m4_flat, distance_m5_flat, distance_m6_flat, distance_m7_flat, distance_m1_flat, distance_m1_flat))

plot_dist_box <- ggplot(data = distance_box, aes(x = mock, y = distance)) +
  geom_boxplot(outlier.size = 3) +
  geom_point(data = aggregate(distance ~ mock, data = distance_box, mean),
             aes(x = mock, y = distance), color = "darkred", size = 7) +
  labs(x = NULL, y = "Phylogenetic distance") +
  paper_theme_fab +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_discrete(labels = c("Mock 1" = "Mock 1\nCo-culture 1", "Mock 2" = "Mock 2\nCo-culture 2", "Mock 3" = "Mock 3\nCo-culture 3", "Mock 4" = "Mock 4", "Mock 5" = "Mock 5", "Mock 6" = "Mock 6", "Mock 7" = "Mock 7", "Mock 8" = "Mock 8", "Mock 9" = "Mock 9"))
print(plot_dist_box)

# 6. Combined plots ----

# Create data frame with members of each mock/co-culture
community_members <- data.frame(community = c("Mock 1/8/9\nCo-culture 1", "Mock 1/8/9\nCo-culture 1", "Mock 1/8/9\nCo-culture 1",
                                              "Mock 2\nCo-culture 2", "Mock 2\nCo-culture 2", "Mock 2\nCo-culture 2", "Mock 2\nCo-culture 2",
                                              "Mock 3\nCo-culture 3", "Mock 3\nCo-culture 3", "Mock 3\nCo-culture 3", "Mock 3\nCo-culture 3", "Mock 3\nCo-culture 3",
                                              "Mock 4\n ", "Mock 4\n ", "Mock 4\n ", "Mock 4\n ", "Mock 4\n ", "Mock 4\n ", "Mock 4\n ", "Mock 4\n ", "Mock 4\n ",
                                              "Mock 5\n ", "Mock 5\n ", "Mock 5\n ", "Mock 5\n ", "Mock 5\n ", "Mock 5\n ", "Mock 5\n ", "Mock 5\n ",
                                              "Mock 6\n ", "Mock 6\n ", "Mock 6\n ", "Mock 6\n ", "Mock 6\n ", "Mock 6\n ", "Mock 6\n ",
                                              "Mock 7\n ", "Mock 7\n ", "Mock 7\n ", "Mock 7\n ", "Mock 7\n ", "Mock 7\n ", "Mock 7\n ",
                                              "Mock 7\n ", "Mock 7\n ", "Mock 7\n ", "Mock 7\n ", "Mock 7\n ", "Mock 7\n ", "Mock 7\n "),
                                member = c("F. nucleatum", "S. oralis", "                                          ",
                                           "F. nucleatum", "P. gingivalis", "S. oralis", "                                          ",
                                           "F. nucleatum", "P. gingivalis", "S. oralis", "V. parvula", "                                          ",
                                           "A. naeslundii", "A. viscosus", "S. gordonii", "S. mitis", "S. oralis", "S. salivarius", "S. sanguinis", "V. parvula", "                                          ",
                                           "S. gordonii", "S. mitis", "S. mutans", "S. oralis", "S. salivarius", "S. sanguinis", "S. sobrinus", "                                          ",
                                           "A. actinomycetemcomitans", "F. nucleatum", "P. gingivalis", "P. intermedia", "S. mutans", "S. sobrinus", "                                          ",
                                           "A. actinomycetemcomitans","A. naeslundii", "A. viscosus", "F. nucleatum", "P. gingivalis", "P. intermedia", "S. gordonii",
                                           "S. mitis", "S. mutans", "S. oralis", "S. salivarius", "S. sanguinis", "S. sobrinus", "V. parvula"))

community_members_wide <- community_members %>% 
  dplyr::group_by(community) %>% 
  dplyr::mutate(row = row_number()) %>% 
  pivot_wider(names_from = community, values_from = member) %>% 
  dplyr::select(-row)
community_members_wide[is.na(community_members_wide)] <- ""

table_theme <- ttheme_minimal(core = list(fg_params = list(fontface = 3, cex = 1.2)),
                              colhead = list(fg_params = list(fontface = 2, cex = 1.6, parse = TRUE))
                              )

separators <- replicate(ncol(community_members_wide) - 1,
                        segmentsGrob(x1 = unit(0, "npc"), gp = gpar(lty = 1)),
                        simplify = FALSE)

community_grob <- gridExtra::tableGrob(community_members_wide, theme = table_theme, rows = NULL)

community_table <- gtable::gtable_add_grob(community_grob, grobs = separators,
                                           t = 2, b = nrow(community_grob), l = seq_len(ncol(community_grob)-1)+1)

community_table <- gtable::gtable_add_grob(community_table, grobs = segmentsGrob(x0 = unit(0, "npc"),
                                                                                 y0 = unit(0, "npc"),
                                                                                 x1 = unit(1, "npc"),
                                                                                 y1 = unit(0, "npc"),
                                                                                 gp = gpar(lwd = 5)),
                                           t = 1, b = 1, l = 1, r = ncol(community_table))

p_table_community <- grid.arrange(community_table)

plot_phylogenetic <- cowplot::plot_grid(plot_nj_tree, plot_dist_box, labels = c("B", "C"), label_size = 34, hjust = -0.15, ncol = 2, nrow = 1)
print(plot_phylogenetic)

plot_phylo_table <- cowplot::plot_grid(p_table_community, NULL, plot_phylogenetic, labels = c("A", "", ""), label_size = 34, hjust = -0.15, nrow = 3, ncol = 1, rel_heights = c(0.6, 0.05, 1))
print(plot_phylo_table)
