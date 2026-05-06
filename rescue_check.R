# ============================================================
# SCRIPT 2: rescue_check.R
# ============================================================
# Place at: ~/Larval_ortholog2026/rescue_check.R
# Run with: Rscript ~/Larval_ortholog2026/rescue_check.R
# ============================================================

library(dplyr); library(readr); library(tidyr)

dtw <- read_tsv("~/Larval_ortholog2026/DTW_clusters_with_annotations.tsv",
                col_types = cols(.default = "c")) %>% select(gene_id, cluster)

broad_og <- read_tsv("~/Larval_ortholog2026/Proteomes_broad/OrthoFinder/Results_broad_run/Orthogroups/Orthogroups.tsv",
                     col_types = cols(.default = "c")) %>%
  select(Orthogroup, Astrangia_poculata) %>%
  separate_rows(Astrangia_poculata, sep = ",\\s*") %>%
  filter(!is.na(Astrangia_poculata) & Astrangia_poculata != "") %>%
  mutate(gene_id = sub("^protein\\|evm\\.model\\.", "", Astrangia_poculata))

cnid_og <- read_tsv("~/Larval_ortholog2026/Proteomes_cnidaria/OrthoFinder/Results_cnidaria_run/Orthogroups/Orthogroups.tsv",
                    col_types = cols(.default = "c")) %>%
  select(Orthogroup, Astrangia_poculata) %>%
  separate_rows(Astrangia_poculata, sep = ",\\s*") %>%
  filter(!is.na(Astrangia_poculata) & Astrangia_poculata != "") %>%
  mutate(gene_id = sub("^protein\\|evm\\.model\\.", "", Astrangia_poculata))

unmapped_broad <- dtw %>% anti_join(broad_og, by = "gene_id")
unmapped_cnid  <- dtw %>% anti_join(cnid_og,  by = "gene_id")
rescued        <- unmapped_broad %>% inner_join(cnid_og, by = "gene_id")
truly_orphan   <- unmapped_broad %>% anti_join(cnid_og,  by = "gene_id")

cat("Unmapped in broad run:", nrow(unmapped_broad), "\n")
cat("Unmapped in cnidaria run:", nrow(unmapped_cnid), "\n")
cat("Rescued by cnidaria run (cnidarian-specific):", nrow(rescued), "\n")
cat("  Cluster 1:", sum(rescued$cluster == "1"), "\n")
cat("  Cluster 2 (ZGA):", sum(rescued$cluster == "2"), "\n")
cat("Truly orphan (Astrangia-specific):", nrow(truly_orphan), "\n")
cat("  Cluster 1:", sum(truly_orphan$cluster == "1"), "\n")
cat("  Cluster 2 (ZGA):", sum(truly_orphan$cluster == "2"), "\n")

write_tsv(rescued,      "~/Larval_ortholog2026/MZT_stability_outputs/cnidaria_rescued_ZGA_genes.tsv")
write_tsv(truly_orphan, "~/Larval_ortholog2026/MZT_stability_outputs/astrangia_specific_ZGA_genes.tsv")
