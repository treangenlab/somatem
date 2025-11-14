# Plot Interpretation Guide

## 1 - Jaccard Heatmap

- **What it shows:** Pairwise dissimilarity (1 - Jaccard similarity) between unitigs, assembly, and all bins based on k-mer hash overlap.
- **What to look for:** Dark colors (low values) indicate high similarity between sequences; light colors suggest sequences are distinct or have little overlap.

## MDS on 1 - Jaccard

- **What it shows:** 2D projection showing relationships between unitigs, assembly, and bins in similarity space.
- **What to look for:** Closely positioned points share more k-mer content; distant points are compositionally different or represent distinct sequences.

## Top Bins by Explained Hashes

- **What it shows:** Absolute count of unitig k-mer hashes explained by each bin (top contributors).
- **What to look for:** Bins with high counts contain sequences well-represented in the unitig set; low counts may indicate incomplete or low-quality bins.

## Unitigs Partition

- **What it shows:** Fraction of unitig hashes found only in assembly (A only), only in bins (B only), in both, or in neither (unexplained).
- **What to look for:** High "both" fraction indicates good agreement between assembly and binning; high "unexplained" suggests missing sequences in both outputs.

## Cumulative Explained vs #Bins

- **What it shows:** Greedy cumulative coverage: how much of the unitig content is explained as you add bins in order of their contribution.
- **What to look for:** A steep initial rise shows a few bins explain most unitigs (good); a slow rise suggests content is spread across many bins or missing.

## Key Metrics Table

- **What it shows:** Summary statistics including unitig/assembly/bin set sizes, overlap fractions, AUC (area under cumulative curve), and PAM score.
- **What to look for:** High "frac_unitigs_in_A" and "frac_unitigs_in_B" indicate good recovery; PAM score near 1.0 means excellent overall quality and balance.

## Sankey Diagram

- **What it shows:** Flow of unitig hashes into four categories: A only, A∩B (both), B only, and unexplained.
- **What to look for:** Thick flows to "A∩B" indicate strong concordance; thick flows to "unexplained" suggest sequences missing from both assembly and bins.

## Per-Bin Explained (%)

- **What it shows:** Percentage of total unitig hashes explained by each of the top bins.
- **What to look for:** Bins with high percentages are major contributors to the unitig pool; uneven distribution may indicate dominant populations or fragmented binning.

## Parameters Table

- **What it shows:** Analysis parameters including k-mer size, scaling factor, random seed, number of bins, and whether RocksDB was used.
- **What to look for:** These settings affect sensitivity and memory usage; lower 'scaled' values provide finer resolution but require more memory.

---

## Overall Interpretation

**PAM (Pigeon Appraisal Metric):** A composite score (0-1) combining explained fraction (60%), AUC (30%), and balance between A-only vs B-only (10%). Higher PAM indicates better assembly-binning concordance.

**Ideal Results:** High overlap in both assembly and bins, minimal unexplained fraction, steep cumulative curve, and PAM > 0.8.
