"""This module describes Contig class"""
from collections import defaultdict


class Contig:
    """Contig objects stores data about an assembled contig: sequence,
    number of mapped reads, list of aligned reads, sizes of alignment etc.

    Attributes:
        contig_id (str): contig identifier
        sequence (str): contig sequence
        read_count (:obj:defaultdict[str, int]): key is sample identifier,
            value is count of reads for the sample
        read_segments (:obj:defaultdict[str, int]):  key is sample identifier,
            value is length of all aligned segments for the sample
        reads (list of str): identifiers of reads mapped to the contig
        genes (dict[str,:obj:Gene]): genes predicted in the contig

    """
    def __init__(self, contig_id='', sequence=''):
        """
        Args:
            contig_id (str): contig identifier
            sequence (str): contig sequence
        """
        self.contig_id = contig_id
        self.sequence = sequence
        self.read_count = defaultdict(int)  # Count of reads for each sample
        self.read_segments = defaultdict(int)  # Length of all aligned segments in each sample
        self.reads = []
        self.genes = {}

    def update_coverage(self, sample, segment_length):
        """Updates read_count and read_segments attributes"""
        if sample in self.read_count:
            self.read_count[sample] += 1
        else:
            self.read_count[sample] = 1
        if sample in self.read_segments:
            self.read_segments[sample] += segment_length
        else:
            self.read_segments[sample] = segment_length

    def get_coverage(self, sample=None):
        """Returns read coverage for a given sample or total coverage"""
        result = 0.0
        if self.sequence == '':
            pass
        elif sample:
            if sample in self.read_count and sample in self.read_segments:
                result = self.read_segments[sample] / len(self.sequence)
        else:
            result = sum(self.read_segments.values()) / len(self.sequence)
        return result

    def add_gene(self, gene):
        """Adds Gene object to dictionary of genes"""
        self.genes[gene.gene_id] = gene

    def get_rpkm(self, sample_readcount, sample=None):
        """Returns RPKM score for a given sample or all samples"""
        result = 0.0
        if sample in self.read_count:
            result = self.read_count[sample] * 1000000000 / len(self.sequence) / sample_readcount
        elif sample is None:
            result = len(self.reads) * 1000000000 / len(self.sequence) / sample_readcount
        return result

    def get_rpm(self, sample_readcount, sample=None):
        """Returns RPM score for a given sample or all samples"""
        result = 0.0
        if sample in self.read_count:
            result = self.read_count[sample] * 1000000 / sample_readcount
        elif sample is None:
            result = len(self.reads) * 1000000 / sample_readcount
        return result

    def get_read_count(self, sample=None):
        """Returns raw read count for a given sample or all samples"""
        result = 0.0
        if sample in self.read_count:
            result = self.read_count[sample]
        elif sample is None:
            result = len(self.reads)
        return result
