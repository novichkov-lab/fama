"""Describes Gene class"""
#!/usr/bin/python
from collections import defaultdict
from lib.utils.const import STATUS_CAND, STATUS_GOOD, STATUS_BAD

class Gene:
    """Gene object stores data for a gene predicted in assembled contig

    Attributes:
        contig_id (str): contig identifier
        gene_id (str): gene identifier
        protein_sequence (str): protein sequence
        start (int): start position
        end (int): end position
        strand (str): strand identifier, either '-1' or '1'
        functions (defaultdict[str:float]): key is function identifier,
            value is a score (currently, RPKM)
        hit_list (:obj:DiamondHitList): list of DIAMOND hits
        status = STATUS_CAND
        status (str): can have one of three possible values defined
            in STATUS_CAND, STATUS_GOOD, STATUS_BAD constants
        taxonomy (str): NCBI Taxonomy ID of Lowest Common Ancestor
    """
    def __init__(self, contig_id='', gene_id='', sequence=None, start=None, end=None, strand=None):
        """Args:
            contig_id (str): contig identifier
            gene_id (str): gene identifier
            sequence (str): protein sequence
            start (int): start position
            end (int): end position
            strand (str): strand identifier, either '-1' or '1'
        """
        self.contig_id = contig_id
        self.gene_id = gene_id
        self.protein_sequence = sequence
        self.start = start
        self.end = end
        self.strand = strand
        self.functions = defaultdict(float)
        self.hit_list = None
        self.status = STATUS_CAND
        self.taxonomy = None

    def set_functions(self, functions):
        """Appends scores from functions argument to scores in self.functions
        attribute

        Args:
            functions (dict[str,float]): key is function identifier, value
                is a score (currently, RPKM)
        """
        for function in functions:
            self.functions[function] += functions[function]

    def set_status(self, status):
        """Changes gene status. STATUS_CAND (default value) can be
        cahnged to STATUS_GOOD or STATUS_BAD. STATUS_BAD can be changed
        to STATUS_CAND or STATUS_GOOD. STATUS_GOOD can be changed only
        to STATUS_CAND.

        Args:
            status (str): can have one of three possible values defined
                in STATUS_CAND, STATUS_GOOD, STATUS_BAD constants
        """
        if self.status == STATUS_CAND or self.status == STATUS_BAD:
            self.status = status
        elif self.status == STATUS_GOOD:
            if status != STATUS_BAD:
                self.status = status
