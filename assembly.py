#! /usr/bin/env python3

from Bio import SeqIO

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('input_file', type = str)
parser.add_argument('output_file', type = str)
parser.add_argument('K', type = int)
args = parser.parse_args()


input_file = args.input_file
output_file = args.output_file
K = args.K

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

    def contigs(self):
        def values(self):
            maximal = 0
            self.values = {}
            for i, j in self.edges.items():
                self.values.setdefault(j, []).append(i)
                if j > maximal:
                    maximal = j

            return maximal

        while self.edges != {}:
            maximal = values(self)
            to_merge = self.values[maximal][0]
            del self.edges[to_merge]

            merged = to_merge[0] + to_merge[1][K-2:]
            edges_copy = list(self.edges.keys())

            if to_merge[0] == to_merge[1]:
                print("TAKIE SAME")

            for i in edges_copy:
                if i[1] == to_merge[0]:
                    new = (i[0], merged)
                    k = self.edges.pop(i)
                    self.edges[new] = k

                if i[0] == to_merge[1]:
                    new = (merged, i[1])
                    try:
                        k = self.edges.pop(i)
                        self.edges[new] = k
                    except KeyError:
                        self.edges[new] = k

            self.nodes.discard(to_merge[0])
            self.nodes.discard(to_merge[1])
            self.nodes.add(merged)

# FILE = "c:/Users/roksa/Desktop/sem7/TWSG2/projekt2/training/reads/reads1.fasta"
FILE = input_file
# K = 18
K = K
THRESHOLD = 1

def main():
    reads = import_reads(FILE)
    corrected_reads = correct_reads(reads, K, "ACTG", THRESHOLD)
    deBruijnG = DeBruijnGraph(corrected_reads, K)
    # deBruijnG = DeBruijnGraph(reads, K)


    deBruijnG.contigs()
    # print(deBruijnG.nodes)


    with open(output_file, 'w+') as f:
        for i, con in enumerate(deBruijnG.nodes):
            f.write(f">read_{i}")
            f.write("\n")
            f.write(str(con))
            f.write('\n')



if __name__ == "__main__":
    main()
