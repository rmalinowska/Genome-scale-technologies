#! /usr/bin/env python3
# authors: Izabela Fedorczyk, Roksana Malinowska, Weronika Trawińska

'''
Parameters (k-mer size and minimal number of k-mer occurrences)
We ran the script for K from 10 to 25 inclusive for training data and after comparing all the results we chose the most optimal
k-mer size (18) and threshold for k-mer filtering (<= 1 occurrences).

For error correction, the total number of occurrences of each 18-mer in the input reads was calculated. 
Then, for those that occur only once, a similar one (Hamming distance = 1) but with more occurrences was found and replaced with it.

Contig assembly
Using corrected reads we build DeBruijn weighted graph, where each node represents a 17-mer present in the reads
and each directed edge represents a 18-mer created from joining two connected 17-mers. 
Weights on the edges correspond to the total number of occurrences of the 18-mer in reads.

We iteratively merge nodes (using greedy approach) that are connected by the edge with highest weight, creating new contig. 
After merging we delete edges that previously started at one node or ended at second node and jump to the next edge. 
We ignore loops (edges starting and ending at the same node). The merging process is finished when there are no more edges in the graph.
'''

from Bio import SeqIO

import argparse
import copy

parser = argparse.ArgumentParser()
parser.add_argument('input_file', type = str)
parser.add_argument('output_file', type = str)
args = parser.parse_args()


input_file = args.input_file
output_file = args.output_file
K = 18

def import_reads(filename):
    reads = []
    data = SeqIO.parse(filename, "fasta")
    for record in data:
        reads.append(record.seq)
    return reads

def kmerHist(reads, k):
    ''' Return k-mer histogram and average # k-mer occurrences '''
    kmerhist = {}
    for read in reads:
        for kmer in [ read[i:i+k] for i in range(len(read)-(k-1)) ]:
            kmerhist[kmer] = kmerhist.get(kmer, 0) + 1
    return kmerhist

def neighbors1mm(kmer, alpha):
    ''' Generate all neighbors at Hamming distance 1 from kmer '''
    neighbors = []
    for j in range(len(kmer)-1, -1, -1):
        oldc = kmer[j]
        for c in alpha:
            if c == oldc: continue
            neighbors.append(kmer[:j] + c + kmer[j+1:])
    return neighbors

def correct1mm(read, k, kmerhist, alpha, thresh):
    ''' Return an error-corrected version of read.  k = k-mer length.
        kmerhist is kmer count map.  alpha is alphabet.  thresh is
        count threshold above which k-mer is considered correct. '''
    # Iterate over k-mers in read
    for i in range(len(read)-(k-1)):
        kmer = read[i:i+k]
        # If k-mer is infrequent...
        if kmerhist.get(kmer, 0) <= thresh:
            # Look for a frequent neighbor
            for newkmer in neighbors1mm(kmer, alpha):
                if kmerhist.get(newkmer, 0) > thresh:
                    # Found a frequent neighbor; replace old kmer
                    # with neighbor
                    read = read[:i] + newkmer + read[i+k:]
                    break
    # Return possibly-corrected read
    return read

def correct_reads(reads, k, alpha, thresh):
    """ Correct all reads iteratively. """
    corrected = []
    khist = kmerHist(reads, k)
    for read in reads:
        corrected.append(correct1mm(read, k, khist, alpha, thresh))
    return corrected

class DeBruijnGraph:
    ''' De Bruijn directed multigraph built from a collection of
        strings. User supplies strings and k-mer length k.  Nodes
        are k-1-mers.  An Edge corresponds to the k-mer that joins
        a left k-1-mer to a right k-1-mer. '''

    @staticmethod
    def chop(st, k):
        ''' Chop string into k-mers of given length '''
        for i in range(len(st)-(k-1)):
            yield (st[i:i+k], st[i:i+k-1], st[i+1:i+k])

    def __init__(self, strIter, k):
        ''' Build de Bruijn multigraph given string iterator and k-mer
            length k. Nodes in the graph are implemented as strings in self.nodes set and edges as tuples
            (left, right):weight, where (left, right) is a directed edge from left to right node and
            weight is number of single edges that would occure between those edges.'''
        self.nodes = set()
        self.edges = {}
        for st in strIter:
            for _, left, right in self.chop(st, k):
                self.nodes.add(left)
                self.nodes.add(right)
                if (left, right) in self.edges.keys():
                    self.edges[(left, right)] += 1
                else:
                    self.edges[(left, right)] = 1


def contigs(graph):

    contigs = list(graph.nodes)
    edges = graph.edges

    edges = list(dict(sorted(edges.items(), key = lambda item: item[1])).items())
    while edges != []:
        to_merge = edges.pop()

        merged = to_merge[0][0] + to_merge[0][1][K-2:]
        next_edges = []

        contigs.remove(to_merge[0][0])
        if to_merge[0][1] in contigs:
            contigs.remove(to_merge[0][1])
        contigs.append(merged)
        for e in edges:
            if e[0][0] == to_merge[0][1]:
                new = (merged, e[0][1])
                weight = e[1]
                next_edges.append((new, weight))

            elif e[0][1] == to_merge[0][0]:
                new = (e[0][0], merged)
                weight = e[1]
                next_edges.append((new, weight))

            elif e[0][0] != to_merge[0][0] and e[0][1] != to_merge[0][1]:
                next_edges.append(e)

        edges = next_edges
    return contigs

FILE = input_file

THRESHOLD = 1


def filter_contigs(contigs, minlen):
    filtered = []
    for contig in contigs:
        if len(contig) >= minlen:
            filtered.append(contig)
    return filtered

def main():
    reads = import_reads(FILE)
    corrected_reads = correct_reads(reads, K, "ACTG", THRESHOLD)
    deBruijnG = DeBruijnGraph(corrected_reads, K)
    CONTIGS = contigs(deBruijnG)
    CONTIGS = filter_contigs(CONTIGS, 300)
    with open(output_file, 'w+') as f:
        for i, con in enumerate(CONTIGS):
            f.write(f">read_{i}\n")
            f.write(str(con)+"\n")

if __name__ == "__main__":
    main()
