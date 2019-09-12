"""Describes AnnotatedRead class"""
from collections import defaultdict
from lib.utils.const import STATUS_CAND, STATUS_GOOD, STATUS_BAD

class AnnotatedRead:
    """AnnotatedRead stores sequence read data (sequence, identifier, quality
    etc.) together with functional and taxonomic annotations.

    Attributes:
        read_id (str): read identifier. A part of first line in FASTQ entry
            up to the first space symbol
        read_id_line (str): For FASTQ entry, this is the entire first line
        sequence (str): For FASTQ entry, this is the entire second line
        quality (str): For FASTQ entry, this is the entire fourth line
        self.line3 (str): For FASTQ entry, this is the entire third line
        self.hit_list (:obj:DiamondHitList): DiamondHitList with all hits for the read
        self.status (str): Can be STATUS_CAND,STATUS_GOOD,STATUS_BAD
        self.functions (:obj:'defaultdict' of float): dictionary with function
            identifiers as keys and RPKM scores as values
        self.pe_id (str): For FASTQ entry, this is the entire first line of paired end
        self.pe_sequence (str): For FASTQ entry, this is the entire second line of paired end
        self.pe_quality (str): For FASTQ entry, this is the entire fourth line of paired end
        self.pe_line3 (str): For FASTQ entry, this is the entire third line of paired end
        self.taxonomy (str): NCBI Taxonomy ID set by LCA algorithm
    """
    def __init__(self, read_id=None):
        """
        Args:
            read_id(str, optional): read identifier
        """
        self.read_id = read_id
        # Please note that read ID may not contain
        # entire ID line from FASTQ file. It contains only part of the
        # line before the first space symbol.
        self.read_id_line = None  # 1st FASTQ line
        self.sequence = None      # 2nd FASTQ line
        self.quality = None       # 4th FASTQ line
        self.line3 = None         # 3rd FASTQ line
        self.hit_list = None      # hit_list is a DiamondHitList object
        self.status = STATUS_CAND
        self.functions = defaultdict(float) # functions dictionary
        self.pe_id = None          # 1st PE FASTQ line
        self.pe_sequence = None    # 2nd PE FASTQ line
        self.pe_quality = None     # 4th PE FASTQ line
        self.pe_line3 = None       # 3rd PE FASTQ line
        self.taxonomy = None       # NCBI Taxonomy ID set by LCA algorithm

    # Functional annotations
    def set_functions(self, functions):
        """Adds up scores from functions argument to functions attribute

        Args:
            functions (:obj:'dict' of float): dictionary with function
                identifiers as keys and RPKM scores as values
        """
        for function in functions:
            self.functions[function] += functions[function]

    def append_functions(self, functions):
        """Adds up scores from functions argument to functions attribute

        Args:
            functions (:obj:'dict' of float): dictionary with function
                identifiers as keys and RPKM scores as values
        """
        for function in functions:
            self.functions[function] += functions[function]

    def set_status(self, status):
        """Sets read status"""
        if status in set([STATUS_CAND, STATUS_GOOD, STATUS_BAD]):
            self.status = status
        else:
            raise ValueError('Unknown read status: ' + status)

    def show_hits(self):
        """Prints hits from hit_list attribute"""
        self.hit_list.print_hits()

    def __str__(self):
        return 'Annotated read: ' + self.status + '\n'+ str(self.hit_list) \
            + '\t' + ','.join(self.functions)
