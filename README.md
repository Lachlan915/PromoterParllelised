# Promoter Hybrid (Most Optimal)
The Promoter program is a bioinformatics tool that uses the Smith-Waterman-Gotoh algorithm for local DNA sequence alignment. It identifies homologous genes and predicts promoter regions in genomic data using BLOSUM62 scoring and Sigma70 consensus patterns, which is ideal for genomic research and DNA analysis.

## Hybrid / Most Optimal Version
In this optimised branch, the following changes were made:

  - **Sequential.java** was modified to introduce the *GeneProcessor* class for efficient parallel gene processing.
  - The original **SmithWatermanGotoh** algorithm was retained, as it provided better performance without parallelisation overhead.

## Optimisation Overview
  1. **Gene Processing in Sequential.java:** The *GeneProcessor* class was introduced to parallelise the processing of gene files. Each file is handled concurrently across multiple CPU cores, significantly reducing the total processing time. This modification allowed the program to benefit from parallel execution while avoiding unnecessary overhead.
  2. **Original SmithWatermanGotoh Algorithm:** After testing, the original sequential Smith-Waterman-Gotoh algorithm was found to be more efficient than its parallelised version due to the reduced overhead associated with thread management and task splitting.

## Performance Improvement
With this hybrid approach, the computational time was reduced to around 32 seconds, achieving a speedup of approximately 3.75 times compared to the original sequential implementation. This version combines efficient parallel file processing with the optimal use of the original SmithWatermanGotoh algorithm for the best overall performance.
