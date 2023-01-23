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

### Approach

one weighted
edge per distinct k-mer








##### Authors: Izabela Fedorczyk, Roksana Malinowska, Weronika Trawińska
