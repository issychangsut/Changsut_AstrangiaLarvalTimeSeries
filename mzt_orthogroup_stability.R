# ============================================================
# SCRIPT 1: mzt_orthogroup_stability.R
# ============================================================
# Place at: ~/Larval_ortholog2026/mzt_orthogroup_stability.R
# Run with: Rscript ~/Larval_ortholog2026/mzt_orthogroup_stability.R
# ============================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

broad_og  <- "~/Larval_ortholog2026/Proteomes_broad/OrthoFinder/Results_broad_run/Orthogroups/Orthogroups.tsv"
cnid_og   <- "~/Larval_ortholog2026/Proteomes_cnidaria/OrthoFinder/Results_cnidaria_run/Orthogroups/Orthogroups.tsv"
dtw_file  <- "~/Larval_ortholog2026/DTW_clusters_with_annotations.tsv"
out_dir   <- "~/Larval_ortholog2026/MZT_stability_outputs"
dir.create(out_dir, showWarnings = FALSE)

dtw <- read_tsv(dtw_file, col_types = cols(.default = "c")) %>%
  select(gene_id, cluster)
cat("DTW genes loaded:", nrow(dtw), "\n")
cat("Clusters:", unique(dtw$cluster), "\n")

load_og_long <- function(path, run_label) {
  og <- read_tsv(path, col_types = cols(.default = "c"))
  apoc_col <- grep("apoculata|poculata", colnames(og), ignore.case = TRUE, value = TRUE)[1]
  cat(run_label, "- Astrangia column:", apoc_col, "\n")
  og %>%
    select(Orthogroup, all_of(apoc_col)) %>%
    rename(genes = all_of(apoc_col)) %>%
    separate_rows(genes, sep = ",\\s*") %>%
    filter(!is.na(genes) & genes != "") %>%
    mutate(genes = sub("^protein\\|evm\\.model\\.", "", genes)) %>%
    mutate(run = run_label)
}

cat("Loading broad run...\n")
broad_long <- load_og_long(broad_og, "broad_10tax")
cat("Loading cnidaria run...\n")
cnid_long  <- load_og_long(cnid_og,  "cnidaria_7tax")

broad_mapped <- dtw %>%
  left_join(broad_long %>% select(gene_id = genes, broad_OG = Orthogroup), by = "gene_id")
cnid_mapped <- dtw %>%
  left_join(cnid_long %>% select(gene_id = genes, cnid_OG = Orthogroup), by = "gene_id")

full_map <- broad_mapped %>%
  left_join(cnid_mapped, by = c("gene_id", "cluster"), relationship = "many-to-many")

write_tsv(full_map, file.path(out_dir, "DTW_genes_orthogroup_assignments.tsv"))
cat("Mapped", sum(!is.na(full_map$broad_OG)), "genes to broad orthogroups\n")
cat("Mapped", sum(!is.na(full_map$cnid_OG)),  "genes to cnidaria orthogroups\n")

summary_broad <- full_map %>%
  filter(!is.na(broad_OG)) %>%
  group_by(cluster) %>%
  summarise(n_genes_mapped = n(), n_orthogroups = n_distinct(broad_OG), .groups = "drop") %>%
  mutate(run = "broad_10tax")

summary_cnid <- full_map %>%
  filter(!is.na(cnid_OG)) %>%
  group_by(cluster) %>%
  summarise(n_genes_mapped = n(), n_orthogroups = n_distinct(cnid_OG), .groups = "drop") %>%
  mutate(run = "cnidaria_7tax")

summary_all <- bind_rows(summary_broad, summary_cnid)
write_tsv(summary_all, file.path(out_dir, "cluster_orthogroup_summary.tsv"))
print(summary_all)

stability <- full_map %>%
  filter(!is.na(broad_OG) & !is.na(cnid_OG)) %>%
  group_by(cluster, broad_OG) %>%
  summarise(
    n_genes    = n(),
    n_cnid_OGs = n_distinct(cnid_OG),
    stability  = ifelse(n_distinct(cnid_OG) == 1, "stable", "split"),
    .groups    = "drop"
  )

write_tsv(stability, file.path(out_dir, "orthogroup_stability_broad_vs_cnidaria.tsv"))
cat("\nStability summary:\n")
print(table(stability$cluster, stability$stability))

breadth_colors <- c(stable = "#2ecc71", split = "#e74c3c")
p2 <- stability %>%
  count(cluster, stability) %>%
  ggplot(aes(x = paste("Cluster", cluster), y = n, fill = stability)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = breadth_colors) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  labs(title = "Orthogroup Stability: Broad vs Cnidaria Run",
       x = "DTW Cluster", y = "Proportion", fill = "Stability")
ggsave(file.path(out_dir, "cluster_stability.png"), p2, width = 7, height = 5, dpi = 150)

cat("\nDone. Outputs in:", out_dir, "\n")
