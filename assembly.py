#! /usr/bin/python

from Bio import SeqIO

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
        print(k,"- mery dla", read,"\n\n", [ read[i:i+k] for i in range(len(read)-(k-1)) ], '\n')
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
            _, left, right = self.chop(st, k)
            self.nodes.add(left)
            self.nodes.add(right)
            if (left, right) in self.edges.keys():
                self.edges[(left, right)] += 1
            else:
                self.edges[(left, right)] = 1
