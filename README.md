# Promoter Hybrid (Most Optimal)
The Promoter program is a bioinformatics tool that uses the Smith-Waterman-Gotoh algorithm for local DNA sequence alignment. It identifies homologous genes and predicts promoter regions in genomic data using BLOSUM62 scoring and Sigma70 consensus patterns, which is ideal for genomic research and DNA analysis.

## Hybrid / Most Optimal Version
In this optimised branch, the following changes were made:

  - **Sequential.java** was modified to introduce the GeneProcessor class for efficient parallel gene processing.
  - The original **SmithWatermanGotoh** algorithm was retained, as it provided better performance without parallelisation overhead.

