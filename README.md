# Promoter Parallelised
The Promoter program is a bioinformatics tool that uses the Smith-Waterman-Gotoh algorithm for local DNA sequence alignment. It identifies homologous genes and predicts promoter regions in genomic data using BLOSUM62 scoring and Sigma70 consensus patterns, which is ideal for genomic research and DNA analysis.

## Full Parallelised
In this parallelized version of the Promoter program, two main files have been modified:

  - **SmithWatermanGotoh.java**
  - **Sequential.java**

### Parallelisation Overview
  1. **SmithWatermanGotoh Algorithm:** The scoring matrix computation has been parallelized to improve efficiency. A threshold of 1500 is used to determine when to split a task for parallel execution. This optimization helps distribute workload across multiple threads, using available CPU cores.
  2. **Gene Processing in Sequential.java:** Previously, gene files were processed sequentially, one at a time. Gene files are divided and processed simultaneously across different CPU cores in the parallelised version. This approach significantly reduces the overall processing time by leveraging concurrent file handling.

### Performance Improvement
With these changes, the computational time was reduced by approximately 1 minute and 15 seconds, achieving a speedup of about 1.6 times. However, further testing showed that using the original SmithWatermanGotoh algorithm was even faster due to reduced overhead, indicating that parallelization is not always the optimal solution.
