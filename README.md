# Genome-scale-technologies

Script is destined to work on single-end reads originating from the same
strand of a single chromosome. Typical parameters of input dataset:
- number of reads: 1000,
- read length: 80bp,
- average percentage of mismatches: ≤ 5%,
- average coverage: ≥ 5×.


### Training data and evaluation

Information about data available in the training/reads directory:
- reads1.fasta – 1000 reads containing ∼1% mismatches,
- reads2.fasta – 1000 reads containing 1-3% mismatches,
- reads3.fasta – 1000 reads containing 3-5% mismatches.

To run the script evaluating assembly bowtie2 program and Python module pysam are required.

## Approach

Parameters (k-mer size and minimal number of k-mer occurrences)
We ran the script for K from 10 to 25 inclusive for training data and after comparing all the results we chose the most optimal
k-mer size (18) and threshold for k-mer filtering (<= 1 occurrences).
For error correction, the total number of occurrences of each 18-mer in the input reads was calculated. 
Then, for those that occur only once, a similar one (Hamming distance = 1) but with more occurrences was found and replaced with it.

#### Contig assembly
Using corrected reads we build DeBruijn weighted graph, where each node represents a 17-mer present in the reads
and each directed edge represents a 18-mer created from joining two connected 17-mers. 
Weights on the edges correspond to the total number of occurrences of the 18-mer in reads.
We iteratively merge nodes (using greedy approach) that are connected by the edge with highest weight, creating new contig. 
After merging we delete edges that previously started at one node or ended at second node and jump to the next edge. 
We ignore loops (edges starting and ending at the same node). The merging process is finished when there are no more edges in the graph.

### Usage
```assembly.py input_reads.fasta output_contigs.fasta```


##### Authors: Izabela Fedorczyk, Roksana Malinowska, Weronika Trawińska
