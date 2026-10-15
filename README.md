

Larval Developmental Trajectories: Expression, Motifs, Ontology, and Putative Regulatory Networks
author: "issy changsut"
date: "2026-10-13"
output:
  html_document:
    toc: true
    toc_float: true
    number_sections: true
  word_document:
    toc: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
```

# Overview

This document records the analytical workflow from raw RNA-seq count data through developmental-stage expression analysis, dynamic time-warping (DTW) clustering, promoter motif discovery, transcription-factor (TF) motif annotation, motif-occurrence mapping, Gene Ontology (GO) over-representation analysis, and construction of putative regulatory modules.

The central objective is to identify developmental expression programs and candidate regulatory mechanisms associated with two DTW-derived transcript clusters in *Astrangia poculata*. The integrated interpretation is deliberately framed as hypothesis generation: promoter motif occurrences and motif-profile similarity nominate putative TF families and target genes, but they do not demonstrate direct TF binding or causal regulation.

# RNA-seq preprocessing

## Inputs

The analysis began with a tab-delimited raw-count matrix in which the first column contained transcript identifiers and remaining columns contained raw RNA-seq counts for individual samples. Sample metadata included `SequenceID` and `developmental_stage`.

## Count and metadata quality control

Count-column names were cleaned to remove the `.trim.sam.counts` suffix and other appended alphabetic components. The following quality-control steps were applied:

- Transcripts were retained when their mean raw count across samples was greater than 3.
- Any sample column containing one or more missing count values was removed.
- Samples with a total raw count of zero were removed.
- A leading `X` automatically added by R to numeric sample IDs was removed before metadata matching.
- Metadata rows lacking developmental-stage annotations were removed.
- Metadata were reordered to exactly match the final count-matrix column order.
- Final checks required no missing counts, no zero-library samples, and no missing developmental-stage labels.

The cleaned count and metadata files were saved for downstream DESeq2, DTW, and motif analyses.

```r
# Principal preprocessing logic
means <- rowMeans(final_reads, na.rm = TRUE)
countFilt <- final_reads[means > 3, , drop = FALSE]

# Remove sample columns with missing counts
na_cols <- colnames(countFilt)[colSums(is.na(countFilt)) > 0]
countFilt <- countFilt[, !colnames(countFilt) %in% na_cols, drop = FALSE]

# Restore numeric sample IDs altered by R
colnames(countFilt) <- sub("^X([0-9].*)$", "\\1", colnames(countFilt))

# Remove zero-total-count samples
total_counts <- colSums(countFilt, na.rm = TRUE)
zero_samples <- names(total_counts)[total_counts == 0]
countFilt <- countFilt[, !colnames(countFilt) %in% zero_samples, drop = FALSE]
```

# DESeq2 developmental-stage analysis

## Model and normalization

A `DESeqDataSet` was constructed with developmental stage as the experimental factor:

\[
\texttt{design = \textasciitilde developmental\_stage}.
\]

Rows with total count less than or equal to 1 after construction of the DESeq2 object were removed. Standard DESeq2 size-factor estimation had produced an NA-related error in this dataset; therefore, manually computed library-size factors were used. Each sample's factor was its library size divided by the geometric mean of all library sizes:

\[
s_j = \frac{L_j}{\exp\left[\frac{1}{m}\sum_{j=1}^{m}\log(L_j)\right]},
\]

where \(L_j\) is the total raw count for sample \(j\).

```r
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = countFilt,
  colData = meta,
  design = ~ developmental_stage
)

dds <- dds[rowSums(counts(dds)) > 1, ]

lib_sizes <- colSums(counts(dds))
manual_size_factors <- lib_sizes / exp(mean(log(lib_sizes)))
sizeFactors(dds) <- manual_size_factors
```

## Global developmental-stage test

A likelihood-ratio test (LRT) identified transcripts whose expression changed across developmental stages overall:

\[
\begin{aligned}
\text{Full model:} &\quad \sim \texttt{developmental\_stage} \\
\text{Reduced model:} &\quad \sim 1.
\end{aligned}
\]

Transcripts with Benjamini-Hochberg-adjusted P values below 0.05 were defined as developmental-stage-responsive. This test evaluates an overall stage effect; it is not a collection of pairwise developmental-stage contrasts.

```r
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)
res_lrt <- results(dds_lrt)
res_lrt_df <- as.data.frame(res_lrt)
res_lrt_df$transcript <- rownames(res_lrt_df)

sig_lrt_df <- subset(res_lrt_df, !is.na(padj) & padj < 0.05)
```

# VST and DTW clustering

## Variance-stabilizing transformation

Variance-stabilized expression values were generated with the DESeq2 VST using `blind = FALSE`, preserving the experimental design while estimating the transformation.

```r
vst_obj <- vst(dds, blind = FALSE)
vst_mat <- assay(vst_obj)
```

## DTW specifications

The implemented trajectory-clustering settings were:

| Component | Specification |
|---|---|
| Input matrix | VST-transformed retained transcripts |
| Genes used | All retained VST transcripts; the optional top-variable-gene selection code was not active |
| Distance metric | Dynamic time warping (DTW) |
| DTW implementation | `proxy::dist(..., method = "DTW")` |
| Clustering method | Hierarchical clustering |
| Linkage | Complete linkage |
| Number of clusters | 2 |
| Cluster extraction | `cutree(hc, k = 2)` |

DTW compares full expression trajectories and can accommodate profiles that have similar temporal shapes but modest timing shifts. The implemented distance calculation used VST profiles at the sample level; it did not explicitly average replicates by developmental stage before calculating DTW distances.

```r
library(proxy)

vst_for_dtw <- as.matrix(vst_mat)
dtw_dist <- proxy::dist(vst_for_dtw, method = "DTW")
hc <- hclust(dtw_dist, method = "complete")
clusters <- cutree(hc, k = 2)

cluster_assignments <- data.frame(
  transcript = rownames(vst_for_dtw),
  cluster = clusters
)
```

## Trajectory summaries

Cluster mean VST expression was calculated across ordered developmental stages. The plotted stage order was:

```text
zygote → 2_cell → 4_cell → blastula → gastrula → planula_48h →
planula_48h_F → planula_72h → planula_72h_F → planula_96h → planula_96h_F
```

Stages labeled with the `_UF` suffix were excluded from trajectory plotting. This exclusion should be stated explicitly wherever the resulting trajectory figure is described.

# Promoter motif analysis

## Promoter sequence extraction

For each DTW cluster, 1-kb promoter regions upstream of annotated transcription start sites were extracted from the repeat-masked *A. poculata* genome. Promoters were assembled into one FASTA file per cluster.

## De novo motif discovery and annotation

The motif workflow was:

1. **STREME** was used to identify de novo motifs enriched in each cluster's promoter set relative to a dinucleotide-preserving shuffled background.
2. Both forward and reverse-complement sequence orientations were evaluated.
3. Motif widths were restricted to 6–20 bp, with up to 10 motifs retained per cluster.
4. **TOMTOM** compared each de novo motif to TF-binding profiles in the JASPAR 2024 database.
5. The best significant TOMTOM/JASPAR profile match was retained as a putative TF-binding-profile and TF-family annotation.
6. **FIMO** scanned cluster-specific promoters for occurrences of selected de novo motifs.
7. To avoid counting repeated sites, only the lowest-P-value motif occurrence was retained for each motif-gene pair.

These tools answer different questions:

| Tool | Meaning of statistical result |
|---|---|
| STREME | Enrichment of a motif in primary promoters relative to the chosen control sequences |
| TOMTOM | Similarity of a de novo motif to a known TF-binding motif profile |
| FIMO | Significance of an individual motif occurrence at a promoter sequence location |

These scores should not be combined into one calibrated P value.

## FIMO column mapping

The FIMO output was parsed using the following column assignments:

| FIMO field | Content |
|---|---|
| Column 1 | `motif_id` |
| Column 2 | `motif_alt_id`, including STREME motif label |
| Column 3 | `sequence_name`, i.e., promoter/gene identifier |

The correct motif-gene parsing logic is therefore:

```python
motif, gene = p[0], p[2]
gene = re.sub(r"^evm\.(TU|model)\.", "", gene)
```

## Motif-gene program output

The motif pipeline generated 16,464 retained motif-gene associations:

| Cluster | Motif-gene records | Candidate TF-gene records |
|---|---:|---:|
| Cluster 1 | 8,565 | 102 |
| Cluster 2 | 7,899 | 118 |
| Total | 16,464 | 220 |

These links represent motif-supported, putative regulatory relationships. They do not demonstrate physical binding by a specific coral TF or causal regulation.

# GO over-representation analysis

## Cluster composition

The two cluster gene lists were disjoint:

| Set | Number of genes |
|---|---:|
| Cluster 1 | 2,219 |
| Cluster 2 | 1,784 |
| Shared genes | 0 |
| Cluster 1 + 2 union | 4,003 |

The annotation file used `evm.model.Ap...` IDs, which were normalized to the `Ap...` identifiers used in the cluster lists. GO annotations were parsed as unique gene-GO pairs to prevent repeated GO IDs on one annotation line from inflating counts.

## Statistical design

The appropriate comparative GO universe is the set of genes in clusters 1 and 2 that possess at least one GO annotation. For each cluster and GO term, a one-sided Fisher exact test evaluates whether that term is overrepresented in the focal cluster relative to the remaining GO-annotated comparison universe. Terms represented by fewer than three genes in the focal cluster are excluded, and P values are adjusted within cluster using Benjamini-Hochberg correction.

\[
\begin{array}{c|cc}
 & \text{Focal cluster} & \text{Other cluster} \\
\hline
\text{Annotated to GO term} & a & c \\
\text{Not annotated to GO term} & b & d
\end{array}
\]

```r
ft <- fisher.test(
  matrix(c(a, b, c, d), nrow = 2, byrow = TRUE),
  alternative = "greater"
)
q_value_BH <- p.adjust(p_values, method = "BH")
```

The final report should explicitly say that this is a **between-cluster comparative GO analysis**. It does not test each cluster against the full genome or all expressed transcripts unless a broader, properly annotated expression-filtered universe is substituted.

## Preliminary functional interpretation

The preliminary GO results indicated coherent but contrasting programs:

| Cluster | Working functional description | Representative terms |
|---|---|---|
| Cluster 1 | Translational-mitochondrial biosynthetic program | Translation, ribosome, structural constituent of ribosome, mitochondrial inner membrane, mitochondrial matrix |
| Cluster 2 | Nuclear transcriptional-regulatory and cell-cycle/genome-maintenance program | Nucleus, DNA binding, DNA-templated transcription, cell division, mitotic cell cycle, DNA repair, oxidative-stress response |

Final effect sizes, raw P values, and BH-adjusted q values should be taken from the corrected analysis that uses the GO-annotated universe consistently.

# Integrating TF motifs with GO

## Putative TF-family target sets

For each cluster and TOMTOM-associated motif or TF family, define a predicted target set as the unique cluster genes with a retained FIMO occurrence. Retain one TF-family-to-gene edge per family-gene pair and preserve the strongest FIMO statistic as an edge attribute.

Recommended minimum filters are:

- At least 10 unique predicted target genes per TF family or motif.
- A predeclared STREME enrichment threshold.
- A significant or otherwise explicitly stated TOMTOM match criterion.
- One retained FIMO hit per motif-gene pair.
- Prefer FIMO q values over unadjusted P values if q values are available and a stringent site filter is desired.

## TF-target GO enrichment

For each TF-family target set, test GO over-representation using all GO-annotated genes in the same DTW cluster as the reference universe. This asks whether motif-bearing targets have a specific functional signature beyond the cluster-wide program.

Recommended filters:

- At least 10 predicted motif-bearing target genes per TF family.
- At least 3 predicted targets annotated to a tested GO term.
- One-sided Fisher exact test.
- Benjamini-Hochberg correction within each TF-family target set.
- Retain terms with \(q < 0.05\) and odds ratio greater than 1.

A supported result can be phrased as:

> Genes with predicted [TF-family]-like motif occurrences in cluster-[X] promoters were enriched for [GO process] annotations, consistent with a putative [TF-family]-associated regulatory module.

This wording does not imply direct binding or direction of regulation.

## Candidate TF prioritization

Prioritize candidate coral TF genes by convergent evidence rather than a motif count alone.

| Evidence layer | Supporting criterion |
|---|---|
| De novo motif enrichment | STREME motif passes a predefined enrichment threshold |
| Motif annotation | TOMTOM/JASPAR match supports a TF family |
| Target support | Motif/family has at least 10–20 unique predicted targets |
| Target functional coherence | GO enrichment among predicted targets is significant |
| TF identity | Coral gene annotation or orthology is consistent with the TF family |
| Co-clustering | Candidate TF gene is assigned to the same DTW cluster |
| Temporal coherence | Candidate TF trajectory is compatible with the target-set or cluster trajectory |
| Mechanistic plausibility | Protein/domain annotation supports DNA binding or transcriptional regulation |

A transparent evidence score may assign one point for each satisfied criterion. It must be described as a prioritization score, not as a P value or formal probability of regulation.

# Putative regulatory network

## Network nodes and edges

Construct an edge list with TF-family/motif nodes, predicted target genes, candidate coral TF genes, and optionally nonredundant GO-process nodes.

### Regulatory edge table

```text
source,target,interaction,cluster,TF_family,motif_id,FIMO_p_value,FIMO_q_value,site_position,strand
```

Example:

```text
Forkhead_like,Ap11.1449,predicted_promoter_target,2,Forkhead,STREME-5,5.92e-12,5.61e-07,479-498,-
```

### Functional edge table

```text
source,target,interaction,cluster,GO_ID,GO_term,odds_ratio,p_value,q_value,target_gene_count
```

Example:

```text
Forkhead_like,GO:0006281,target_set_GO_enriched,2,GO:0006281,DNA repair,3.29,1.63e-07,3.31e-05,62
```

## Network visualization

For a main-text figure, show only the best-supported 3–5 motif families per cluster and a small number of nonredundant GO Biological Process terms. Full TF-family-to-target edges should be retained as a supplementary table.

Suggested styling:

- Size TF-family nodes by number of unique predicted targets.
- Color target genes by DTW cluster or GO program.
- Scale edge width or opacity by \(-\log_{10}\) FIMO P value or q value.
- Show TF-family-to-GO edges only when target-set GO enrichment passes BH correction.
- Keep a clear legend identifying all promoter-motif links as predicted.

# Interpretation and limitations

- DTW clustering groups transcripts by similarity of observed expression profiles; it does not establish common regulation.
- The number of clusters was specified as \(k = 2\); the implementation does not test whether two is the globally optimal cluster number.
- DTW distances were calculated using sample-level VST profiles rather than explicitly stage-averaged trajectories, so replicate-level variation can influence distances.
- Manual library-size normalization should be reported clearly and, where possible, compared with standard DESeq2 median-ratio normalization after resolving the original estimation problem.
- STREME enrichment, TOMTOM motif similarity, and FIMO site significance test different hypotheses and are not interchangeable.
- TOMTOM/JASPAR annotations identify a putative TF-binding profile or family, not necessarily the identity or activity of a specific coral TF.
- Promoter motif occurrence supports a predicted TF-family-to-gene relationship but does not establish physical binding, activation/repression direction, or causality.
- GO over-representation identifies annotation patterns relative to the selected reference universe; it is not a direct measurement of biological activity.
- Stronger validation would come from compatible TF expression timing, chromatin accessibility data, TF occupancy assays such as ChIP-seq or CUT&RUN, or TF perturbation experiments.

# Methods-ready summary

Raw transcript counts were quality filtered by removing samples containing missing values or zero total counts and retaining transcripts with mean count greater than 3 across samples. Samples were matched to metadata using sequence identifiers, and developmental stage was modeled using DESeq2 with manually calculated library-size factors. Stage-responsive transcripts were identified using a likelihood-ratio test comparing a full model containing developmental stage to an intercept-only reduced model; significance was defined at Benjamini-Hochberg-adjusted \(P < 0.05\). Variance-stabilized expression values were generated with the DESeq2 variance-stabilizing transformation (`blind = FALSE`). Transcript trajectories were compared using dynamic time-warping distances and clustered by complete-linkage hierarchical clustering; the dendrogram was cut into two clusters. Cluster-level expression trajectories were summarized as mean VST expression across ordered developmental stages, excluding stages labeled `_UF` from visualization.

For regulatory analysis, 1-kb upstream promoter regions of cluster-assigned genes were extracted from the repeat-masked *A. poculata* genome. STREME was used to identify cluster-enriched de novo motifs using dinucleotide-preserving shuffled controls, followed by TOMTOM comparison to JASPAR transcription-factor binding profiles. FIMO was used to identify promoter-level occurrences of retained motifs, and the lowest-P-value hit was retained per motif-gene pair. Motif-bearing genes were treated as putative targets of the TOMTOM-annotated TF family. Gene Ontology over-representation analysis was performed using one-sided Fisher exact tests with Benjamini-Hochberg correction, using GO-annotated genes assigned to the compared trajectory clusters as the reference universe. Putative regulatory modules were defined by linking motif-derived TF families, FIMO-predicted target genes, and GO terms enriched among those targets.
