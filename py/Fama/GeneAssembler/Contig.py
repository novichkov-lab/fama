#!/usr/bin/python
from collections import defaultdict

class Contig:
    def __init__(self, contig_id = "", sequence = ""):
        self.contig_id = contig_id
        self.sequence = sequence
        self.read_count = defaultdict(int) # Count of reads for each sample
        self.read_segments = defaultdict(int) # Length of all aligned segments in each sample
        self.reads = []
        self.genes = {}

    def update_coverage(self, sample, segment_length):
        if sample in self.read_count:
            self.read_count[sample] += 1
        else:
            self.read_count[sample] = 1
        if sample in self.read_segments:
            self.read_segments[sample] += segment_length
        else:
            self.read_segments[sample] = segment_length

    def get_coverage(self, sample = None):
        if self.sequence == '':
            return 0.0
        if sample:
            if sample in self.read_count and sample in self.read_segments:
                return self.read_segments[sample] / len(self.sequence)
            else:
                return 0.0
        else:
            return sum(self.read_segments.values()) / len(self.sequence)

    def add_gene(self, gene):
        self.genes[gene.gene_id] = gene

    def get_rpkm(self, sample_readcount, sample = None):
        if sample:
            if sample in self.read_count:
                return self.read_count[sample] * 1000000000 / len(self.sequence) / sample_readcount
            else:
                return 0.0
        else:
            return len(self.reads) * 1000000000 / len(self.sequence) / sample_readcount

    def get_rpm(self, sample_readcount, sample = None):
        if sample:
            if sample in self.read_count:
                return self.read_count[sample] * 1000000 / sample_readcount
            else:
                return 0.0
        else:
            return len(self.reads) * 1000000 / sample_readcount

    def get_read_count(self, sample = None):
        if sample:
            if sample in self.read_count:
                return self.read_count[sample]
            else:
                return 0.0
        else:
            return len(self.reads)
