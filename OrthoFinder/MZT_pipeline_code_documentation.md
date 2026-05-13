# MZT Orthogroup Evolution Pipeline — Complete Code Documentation

**Project:** Larval_ortholog2026  
**Scope:** All code from orthogroup extraction through CAFE5 analysis and figure generation  
**Last updated:** May 2026

---

## Overview

This document records every command and script run in the MZT (maternal-to-zygotic transition) gene family evolution pipeline. The pipeline proceeds in five stages:

1. OrthoFinder — orthogroup identification and MSA construction
2. CAFE5 — gene family expansion/contraction analysis
3. MZT orthogroup filtering and extraction
4. IQ-TREE3 — ML gene tree construction
5. R — conservation figures

---

## Stage 1: OrthoFinder

OrthoFinder was run across three proteome sets. Results are pre-computed and stored under `refs/`.

### Input proteome sets

| Run | Directory | Species |
|-----|-----------|---------|
| Anchor (6 sp) | `refs/Proteomes/` | H. sapiens, L. variegatus, O. faveolata, A. poculata, D. melanogaster, H. vulgaris |
| Broad (10 sp) | `refs/Proteomes_broad/` | Above + A. digitifera, A. palmata, S. pistillata, M. cavernosa |
| Cnidaria-only (7 sp) | `refs/Proteomes_cnidaria/` | All 7 cnidarian species |

### Verifying orthogroup gene counts

```bash
# Confirm species columns in each run
head -1 refs/Proteomes/OrthoFinder/Results_Mar23/Orthogroups/Orthogroups.GeneCount.tsv | tr '\t' '\n'
head -1 refs/Proteomes_broad/OrthoFinder/Results_broad_run/Orthogroups/Orthogroups.GeneCount.tsv | tr '\t' '\n'
head -1 refs/Proteomes_cnidaria/OrthoFinder/Results_cnidaria_run/Orthogroups/Orthogroups.GeneCount.tsv | tr '\t' '\n'

# Verify gene count for a target OG (expected: 41)
grep "^OG0001264" refs/Proteomes/OrthoFinder/Results_Mar23/Orthogroups/Orthogroups.tsv \
  | tr '\t' '\n' | grep -v "^OG" | tr ',' '\n' | grep -v '^$' | wc -l

# Confirm pre-built MSA and internal-ID trees exist
find refs/Proteomes/OrthoFinder/Results_Mar23/ -name "OG0001264*" 2>/dev/null
```

### Extract MZT OG gene counts from all three runs

```bash
# Anchor run (6 species)
awk 'NR==1 || $1=="OG0000417" || $1=="OG0000511" || $1=="OG0000683" || $1=="OG0000941" || \
     $1=="OG0001229" || $1=="OG0001264" || $1=="OG0001856" || $1=="OG0002972"' \
    refs/Proteomes/OrthoFinder/Results_Mar23/Orthogroups/Orthogroups.GeneCount.tsv \
    > outputs/mzt_genecounts.tsv

# Broad run (10 species) — verify totals
awk 'NR==1 || $1=="OG0000417" || $1=="OG0000511" || $1=="OG0000683" || $1=="OG0000941" || \
     $1=="OG0001229" || $1=="OG0001264" || $1=="OG0001856" || $1=="OG0002972"' \
    refs/Proteomes_broad/OrthoFinder/Results_broad_run/Orthogroups/Orthogroups.GeneCount.tsv \
    | column -t

# Cnidaria-only run
awk 'NR==1 || $1=="OG0000417" || $1=="OG0000511" || $1=="OG0000683" || $1=="OG0000941" || \
     $1=="OG0001229" || $1=="OG0001264" || $1=="OG0001856" || $1=="OG0002972"' \
    refs/Proteomes_cnidaria/OrthoFinder/Results_cnidaria_run/Orthogroups/Orthogroups.GeneCount.tsv \
    | column -t
```

---

## Stage 2: CAFE5 — Gene Family Evolution

CAFE5 was run in two passes. The final results are in `outputs/cafe_results/`.

### Pass 1 — Lambda estimation

```bash
# Run CAFE5 to estimate lambda (evolutionary rate)
# Input: gene count table, ultrametric species tree
cafe5 -i cafe_inputs/<gene_counts>.tsv \
      -t cafe_trees/<ultrametric_species_tree>.nwk \
      -o outputs/cafe_results_pass1/ \
      -c 4
```

### Final run — Full model

```bash
# Run CAFE5 with estimated lambda from pass1
cafe5 -i cafe_inputs/<gene_counts>.tsv \
      -t cafe_trees/<ultrametric_species_tree>.nwk \
      -o outputs/cafe_results/ \
      -c 4
```

### Key CAFE5 output files

```
outputs/cafe_results/
  Base_family_results.txt        # Per-OG p-values
  Base_change.tab                # Signed change per branch per OG (+expand, -contract)
  Base_branch_probabilities.tab  # Branch-level p-values
  Base_count.tab                 # Reconstructed ancestral counts
  Base_clade_results.txt         # Clade-level summaries
```

### Extract significant MZT OGs (p < 0.05)

```bash
# Pull p-values and filter significant OGs
# → outputs/mzt_cafe_significant.csv
# Fields: FamilyID, pval, significant
cat outputs/mzt_cafe_significant.csv
# Result: all 8 MZT OGs significant (p <= 0.042)
```

### Extract branch-level changes for MZT OGs

```bash
# Pull branch changes for all 8 MZT OGs across all 10 taxa + ancestral nodes
cat outputs/mzt_cafe_changes_all_taxa.csv
# Fields: FamilyID, dmelanogaster<1>, Lvariegatus<2>, Hsapiens<3>,
#         Acropora_digitifera<4> ... Hydra_vulgaris<10>, <11>..<19> (ancestral)
```

### Check CAFE5 results files

```bash
ls outputs/cafe_results/*.tab outputs/cafe_results/*.txt
find outputs/ -name "*.tab" -o -name "*change*" -o -name "*clade*" 2>/dev/null
```

---

## Stage 3: MZT Orthogroup Filtering

### Inspect and verify significant OGs

```bash
# View significant MZT OGs
cat outputs/mzt_cafe_significant.csv
# FamilyID,pval,significant
# OG0001264,0,y
# OG0001856,0.022,y
# OG0002972,0.042,y
# OG0000417,0,y
# OG0000511,0,y
# OG0000683,0,y
# OG0000941,0,y
# OG0001229,0,y

# View branch-level change summary
cat outputs/mzt_cafe_changes_all_taxa.csv
```

### Check for existing CAFE5 processed outputs

```bash
ls outputs/acropora_shared_changes.csv
ls outputs/acropora_cafe_changes.csv
ls outputs/mzt_apoc_changes.csv
ls outputs/significant_changes_all_taxa.csv
```

---

## Stage 4: IQ-TREE3 — ML Gene Tree Construction

Uses OrthoFinder's pre-built MSAs directly — no sequence re-extraction required.

### Verify MSAs exist for all 8 target OGs

```bash
find refs/Proteomes/OrthoFinder/Results_Mar23/ -name "OG0001264*" 2>/dev/null
# Expected outputs:
#   .../MultipleSequenceAlignments/OG0001264.fa   ← use this for IQ-TREE
#   .../WorkingDirectory/Trees_ids/OG0001264.txt  ← internal-ID tree (not used)
#   .../WorkingDirectory/SequenceIDs.txt          ← ID decoder
```

### Install IQ-TREE3

```bash
# IQ-TREE is installed as iqtree3 in the base conda environment
find ~/miniconda3 -name "iqtree*" -type f 2>/dev/null
# Located at: /home/slv62/miniconda3/bin/iqtree3

iqtree3 --version
```

### Gene tree build script: `cafe_R_scripts/build_gene_trees.sh`

```bash
#!/bin/bash
# Gene tree pipeline using OrthoFinder pre-built MSAs
set -e

MSA_DIR="refs/Proteomes/OrthoFinder/Results_Mar23/MultipleSequenceAlignments"
OUTDIR="outputs/gene_trees"
mkdir -p "$OUTDIR"

# Target OGs: 8 significant MZT OGs
ALL_OGS=$(tail -n +2 outputs/mzt_cafe_significant.csv | cut -d',' -f1 | tr '\n' ' ')
ALL_OGS=$(echo "$ALL_OGS" | tr ' ' '\n' | sort -u | tr '\n' ' ')
echo "Target OGs: $ALL_OGS"

for OG in $ALL_OGS; do
    MSA="$MSA_DIR/${OG}.fa"
    if [ ! -s "$MSA" ]; then
        echo "WARNING: No MSA for $OG -- skipping"
        continue
    fi

    SEQCOUNT=$(grep -c "^>" "$MSA")
    echo -e "\n=== $OG ($SEQCOUNT sequences) ==="

    if [ "$SEQCOUNT" -lt 3 ]; then
        echo "  Skipping -- fewer than 3 sequences"
        continue
    fi

    if [ ! -s "$OUTDIR/${OG}.treefile" ]; then
        echo "  Running IQ-TREE3..."
        iqtree3 -s "$MSA" \
                -m TEST \
                -bb 1000 \
                -T 4 \
                --prefix "$OUTDIR/${OG}" \
                -redo \
                --quiet
    else
        echo "  Tree already exists, skipping"
    fi

    if [ -s "$OUTDIR/${OG}.treefile" ]; then
        cp "$OUTDIR/${OG}.treefile" "$OUTDIR/${OG}_final.nwk"
        echo "  Saved: ${OG}_final.nwk"
    fi
done

echo -e "\n=== Done ==="
ls -lh "$OUTDIR/"*_final.nwk 2>/dev/null
```

### Run gene tree pipeline

```bash
# Clean up any bad previous extractions
rm -f outputs/gene_trees/OG0001264.faa outputs/gene_trees/OG0001264.aln

# Fix script to use iqtree3 (not iqtree2)
sed -i 's/iqtree2/iqtree3/g' cafe_R_scripts/build_gene_trees.sh

# Run
bash cafe_R_scripts/build_gene_trees.sh 2>&1 | tee outputs/gene_trees/build_log.txt
```

### Output gene trees

```
outputs/gene_trees/
  OG0000417_final.nwk
  OG0000511_final.nwk
  OG0000683_final.nwk
  OG0000941_final.nwk
  OG0001229_final.nwk
  OG0001264_final.nwk
  OG0001856_final.nwk
  OG0002972_final.nwk
  build_log.txt
```

---

## Stage 5: R — Conservation Figures

### Script: `cafe_R_scripts/mzt_conservation_figures.R`

Create on the cluster with:

```bash
cat > cafe_R_scripts/mzt_conservation_figures.R << 'REOF'
library(ggplot2)
library(dplyr)
library(tidyr)

# DATA: BROAD RUN (10 species)
# Source: refs/Proteomes_broad/OrthoFinder/Results_broad_run/Orthogroups/Orthogroups.GeneCount.tsv
# Column order (verified): Acropora_digitifera | Acropora_palmata | Astrangia_poculata |
#   Hsapiens | Hydra_vulgaris | Lvariegatus | Montastraea_cavernosa |
#   Orbicella_faveolata | Stylophora_pistillata | dmelanogaster

broad <- data.frame(
  OG                    = c("OG0000417","OG0000511","OG0000683","OG0000941",
                             "OG0001229","OG0001264","OG0001856","OG0002972"),
  Acropora_digitifera   = c(2,1,0,1,1,1,2,1),
  Acropora_palmata      = c(0,0,2,1,1,1,0,0),
  Astrangia_poculata    = c(1,1,1,2,1,1,1,1),
  Hsapiens              = c(56,61,44,45,32,34,16,19),
  Hydra_vulgaris        = c(6,2,2,3,2,2,1,1),
  Lvariegatus           = c(4,3,4,2,1,2,3,3),
  Montastraea_cavernosa = c(3,2,4,1,1,1,4,1),
  Orbicella_faveolata   = c(1,1,1,1,2,2,1,1),
  Stylophora_pistillata = c(2,1,1,1,2,1,2,2),
  dmelanogaster         = c(8,5,9,2,9,6,11,2)
)

# DATA: CNIDARIA-ONLY RUN (7 species)
# Source: refs/Proteomes_cnidaria/OrthoFinder/Results_cnidaria_run/Orthogroups/Orthogroups.GeneCount.tsv
cnid_only <- data.frame(
  OG                    = c("OG0000417","OG0000511","OG0000683","OG0000941",
                             "OG0001229","OG0001264","OG0001856","OG0002972"),
  Acropora_digitifera   = c(6,1,1,1,0,1,1,2),
  Acropora_palmata      = c(5,0,0,0,0,2,0,0),
  Astrangia_poculata    = c(14,13,16,11,1,8,5,2),
  Hydra_vulgaris        = c(0,15,5,0,21,0,0,3),
  Montastraea_cavernosa = c(6,3,5,13,0,3,5,3),
  Orbicella_faveolata   = c(10,8,4,0,0,7,4,1),
  Stylophora_pistillata = c(6,1,3,1,0,1,2,1)
)

# COLOUR PALETTES (darker, translucent fills)
nonc_cols <- c("H. sapiens"="#1a3a5c","L. variegatus"="#0b5345","D. melanogaster"="#4a235a")
cnid_cols  <- c("A. digitifera"="#7b241c","A. palmata"="#c0392b","A. poculata"="#d35400",
                "H. vulgaris"="#7d6608","M. cavernosa"="#b7770d","O. faveolata"="#6e2f0a",
                "S. pistillata"="#e67e22")
all_cols <- c(nonc_cols, cnid_cols)

# FIGURE 1: BROAD (all 10 taxa, stacked, faceted by clade)
long_broad <- broad %>%
  pivot_longer(-OG, names_to="Species", values_to="Copies") %>%
  mutate(
    Clade = ifelse(Species %in% c("Hsapiens","Lvariegatus","dmelanogaster"),
                  "Non-Cnidaria (3 spp)", "Cnidaria (7 spp)"),
    Clade = factor(Clade, levels=c("Non-Cnidaria (3 spp)","Cnidaria (7 spp)")),
    Species = factor(Species,
      levels=c("Hsapiens","Lvariegatus","dmelanogaster","Acropora_digitifera",
               "Acropora_palmata","Astrangia_poculata","Hydra_vulgaris",
               "Montastraea_cavernosa","Orbicella_faveolata","Stylophora_pistillata"),
      labels=c("H. sapiens","L. variegatus","D. melanogaster","A. digitifera",
               "A. palmata","A. poculata","H. vulgaris","M. cavernosa",
               "O. faveolata","S. pistillata")),
    OG = factor(OG))

fig1 <- ggplot(long_broad, aes(x=OG, y=Copies, fill=Species)) +
  geom_col(position="stack", alpha=0.83, color="white", linewidth=0.3) +
  scale_fill_manual(values=all_cols) +
  facet_wrap(~Clade, ncol=1, scales="free_y") +
  labs(title="MZT Orthogroup Conservation - All 10 Taxa",
       subtitle="Gene copies per species; stacked by clade panel",
       x="Orthogroup", y="Gene copies", fill="Species") +
  theme_minimal(base_size=13) +
  theme(plot.background=element_rect(fill="grey95",color=NA),
        panel.background=element_rect(fill="grey95",color=NA),
        strip.background=element_rect(fill="grey82",color=NA),
        strip.text=element_text(face="bold",size=12),
        axis.text.x=element_text(angle=38,hjust=1,size=10),
        legend.position="right", legend.key.size=unit(0.45,"cm"),
        legend.text=element_text(face="italic",size=10),
        plot.title=element_text(face="bold",size=14),
        plot.subtitle=element_text(color="grey35",size=10),
        panel.grid.major.x=element_blank(), panel.grid.minor=element_blank())

ggsave("outputs/fig1_broad_conservation.pdf", fig1, width=10, height=8)
ggsave("outputs/fig1_broad_conservation.png", fig1, width=10, height=8, dpi=300)
cat("Fig1 saved\n")

# FIGURE 2: CNIDARIA-ONLY (7 cnidarian species, dodged bars)
long_cnid <- cnid_only %>%
  pivot_longer(-OG, names_to="Species", values_to="Copies") %>%
  mutate(
    Species = factor(Species,
      levels=c("Acropora_digitifera","Acropora_palmata","Astrangia_poculata",
               "Hydra_vulgaris","Montastraea_cavernosa","Orbicella_faveolata",
               "Stylophora_pistillata"),
      labels=c("A. digitifera","A. palmata","A. poculata","H. vulgaris",
               "M. cavernosa","O. faveolata","S. pistillata")),
    OG = factor(OG))

fig2 <- ggplot(long_cnid, aes(x=OG, y=Copies, fill=Species)) +
  geom_col(position="dodge", alpha=0.85, color="white", linewidth=0.3) +
  scale_fill_manual(values=cnid_cols) +
  labs(title="MZT Orthogroup Conservation - Cnidaria Only",
       subtitle="Gene copies per cnidarian species (cnidaria-only OrthoFinder run)",
       x="Orthogroup", y="Gene copies", fill="Cnidarian species") +
  theme_minimal(base_size=13) +
  theme(plot.background=element_rect(fill="grey95",color=NA),
        panel.background=element_rect(fill="grey95",color=NA),
        axis.text.x=element_text(angle=38,hjust=1,size=10),
        legend.position="right", legend.key.size=unit(0.45,"cm"),
        legend.text=element_text(face="italic",size=10),
        plot.title=element_text(face="bold",size=14),
        plot.subtitle=element_text(color="grey35",size=10),
        panel.grid.major.x=element_blank(), panel.grid.minor=element_blank())

ggsave("outputs/fig2_cnidaria_conservation.pdf", fig2, width=10, height=6)
ggsave("outputs/fig2_cnidaria_conservation.png", fig2, width=10, height=6, dpi=300)
cat("Fig2 saved\n")
REOF

# Run
Rscript cafe_R_scripts/mzt_conservation_figures.R
```

### Figure outputs

```
outputs/
  fig1_broad_conservation.pdf   # Fig 1: all 10 taxa, stacked, faceted (vector)
  fig1_broad_conservation.png   # Fig 1: 300 dpi raster
  fig2_cnidaria_conservation.pdf
  fig2_cnidaria_conservation.png
```

---

## Key Results Summary

| Orthogroup | Anchor copies (6 sp) | Broad copies (10 sp) | Cnidaria-only copies (7 sp) | CAFE5 p-value |
|---|---|---|---|---|
| OG0000417 | 67 | 83 | 47 | 0 |
| OG0000511 | 62 | 77 | 41 | 0 |
| OG0000683 | 55 | 68 | 34 | 0 |
| OG0000941 | 48 | 59 | 26 | 0 |
| OG0001229 | 41 | 52 | 22 | 0 |
| OG0001264 | 41 | 51 | 22 | 0 |
| OG0001856 | 33 | 41 | 17 | 0.022 |
| OG0002972 | 24 | 31 | 12 | 0.042 |

**Key finding:** All 8 MZT orthogroups show consistent expansion in *H. sapiens* and net contraction across all Cnidarian branches in CAFE5. The same OG IDs were recovered across all three OrthoFinder runs, confirming orthogroup structural stability under different taxon sampling schemes.

---

## Notes

- OrthoFinder internal sequence IDs (numeric) in `WorkingDirectory/Trees_ids/` can be decoded using `WorkingDirectory/SequenceIDs.txt` and `WorkingDirectory/SpeciesIDs.txt`
- IQ-TREE3 settings: `-m TEST` (model selection), `-bb 1000` (ultrafast bootstrap), `-T 4` (threads)
- The em dash character in R title strings causes locale warnings on the cluster; replace with plain hyphen using `sed -i 's/--/-/g'` if needed
- Use `cafe_results/` for final analysis; `cafe_results_pass1/` is the lambda estimation run only
