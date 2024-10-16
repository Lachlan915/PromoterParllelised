# Promoter Parllelised
The Promoter program is a bioinformatics tool that uses the Smith-Waterman-Gotoh algorithm for local DNA sequence alignment. It identifies homologous genes and predicts promoter regions in genomic data using BLOSUM62 scoring and Sigma70 consensus patterns, ideal for genomic research and DNA analysis.

## Full Parallelised
In this version of the Promoter program, the following files have been changed:
  <br/>**- SmithWatermanGotoh.java**
  <br/>**- Sequential.java**
<br/>In the **SmithWatermanGotoh** algorithm, the computation of creating the scoring matrix was parallelised to help with efficiency. The threadshold is placed at 1500. This threadshold is used to determine if a particular task needs to split. In the **Sequential** class, it was found that the gene files were read one at a time and what better way than parallelising the class for the files to be split up and be worked on diffrent CPU cores. 

<br/>All in all this roughly reduced the computational time around 1 minute and 15 seconds so about 1.6 times faster. After that an even faster method was discovered, by using the original **SmithWatermanGotoh** algorithm. This is due to the fact that there was too much over head with the parallelised version of the **SmithWatermanGotoh** algorithm.
