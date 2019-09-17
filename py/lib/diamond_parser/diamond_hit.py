"""Defines DiamondHit class"""


class DiamondHit(object):
    """DiamondHit object stores data of one DIAMOND hit.

    DIAMOND hit is a pair of sequences, one is query and other is subject,
    which have some similarity. This similarity is described by amino acid
    identity percent, e-value and bit-score. The area of similarity (or
    alignment) is described by start and end coordinates in both sequences.
    Length of alignment is a distance from the first position to the last
    position of the alignemnt in the query sequence.

    Attributes:
        query_id (str): query sequence identifier (usually, query is a sequence read or a protein)
        subject_id (str): subject sequence identifier (subject is always a reference protein)
        identity (float): amino acid identity %
        length (int): length of alignment
        mismatch (int): number of mismatches
        s_len (int): length of subject sequence
        q_start (int): first position of alignment in the query sequence
        q_end (int): last position of alignment in the query sequence
        s_start (int): first position of alignment in the subject sequence
        s_end (int): last position of alignment in the subject sequence
        evalue (float): e-value of the hit
        bitscore (float): bit-score of the hit
        functions (list): list of function identifiers

    """

    def __init__(self, query_id="", subject_id="", identity=0.0,
                 length=0, mismatch=0, s_len=0, q_start=0,
                 q_end=0, s_start=0, s_end=0, evalue=0.0,
                 bitscore=0.0):
        """Args:
            query_id (str): query sequence identifier (usually, query is a read or a protein)
            subject_id (str): subject sequence identifier (subject is always a reference protein)
            identity (float): amino acid identity %
            length (int): length of alignment
            mismatch (int): number of mismatches
            s_len (int): length of subject sequence
            q_start (int): first position of alignment in the query sequence
            q_end (int): last position of alignment in the query sequence
            s_start (int): first position of alignment in the subject sequence
            s_end (int): last position of alignment in the subject sequence
            evalue (float): e-value of the hit
            bitscore (float): bit-score of the hit

        """
        self.query_id = query_id
        self.subject_id = subject_id
        self.identity = identity
        self.length = length
        self.mismatch = mismatch
        self.s_len = s_len
        self.q_start = q_start
        self.q_end = q_end
        self.s_start = s_start
        self.s_end = s_end
        self.evalue = evalue
        self.bitscore = bitscore
        self.functions = []

    def create_hit(self, tabular_output_fields):
        """Fills DiamondHit attributes with values from DIAMOND output in
        tabular format

        Args:
            tabular_output_fields (list): list of values for 10 fields of
                DIAMOND output in tabular format

        Raises:
            TypeError: field value cannot be converted to int or float
            IndexError: field list is shorter than expected
        """
        try:
            tabular_output_fields[2] = float(tabular_output_fields[2])
            tabular_output_fields[3] = int(tabular_output_fields[3])
            tabular_output_fields[4] = int(tabular_output_fields[4])
            tabular_output_fields[5] = int(tabular_output_fields[5])
            tabular_output_fields[6] = int(tabular_output_fields[6])
            tabular_output_fields[7] = int(tabular_output_fields[7])
            tabular_output_fields[8] = int(tabular_output_fields[8])
            tabular_output_fields[9] = int(tabular_output_fields[9])
            tabular_output_fields[10] = float(tabular_output_fields[10])
            tabular_output_fields[11] = float(tabular_output_fields[11])
        except TypeError as detail:
            print('Value parsing failed for read ' + tabular_output_fields[0])
            print(str(detail))
            raise
        except IndexError as detail:
            print('Wrong number of fields in tabular output for read ' + tabular_output_fields[0])
            print(str(detail))
            raise
        self.query_id = tabular_output_fields[0]
        self.subject_id = tabular_output_fields[1]
        self.identity = tabular_output_fields[2]
        self.length = tabular_output_fields[3]
        self.mismatch = tabular_output_fields[4]
        self.s_len = tabular_output_fields[5]
        self.q_start = tabular_output_fields[6]
        self.q_end = tabular_output_fields[7]
        self.s_start = tabular_output_fields[8]
        self.s_end = tabular_output_fields[9]
        self.evalue = tabular_output_fields[10]
        self.bitscore = tabular_output_fields[11]

    def import_hit(self, entry_tokens):
        """Fills DiamondHit attributes with values from DIAMOND output in
        tabular format, with the last element of the list containing
        pipe-separated function identifiers

        Args:
        entry_tokens (list): list of values for 10 fields of DIAMOND output
            in tabular format, with 11th element storing identifiers
            of functions separated by pipe symbol

        Raises:
            IndexError: field list is shorter than expected
        """
        self.create_hit(entry_tokens[:-1])
        try:
            self.functions = entry_tokens[12].split('|')
        except IndexError:
            print('Unable to parse function list:', entry_tokens)

    def annotate_hit(self, ref_data):
        """Fills functions with list of functions assigned to the subject sequence

        """
        if self.subject_id:
            self.functions = ref_data.lookup_protein_function(self.subject_id)

    def __str__(self):
        """Returns string representation of DiamondHit, tab-separated"""
        return "\t".join((self.query_id, self.subject_id, str(self.identity),
                          str(self.length), str(self.mismatch), str(self.s_len),
                          str(self.q_start), str(self.q_end), str(self.s_start),
                          str(self.s_end), str(self.evalue), str(self.bitscore),
                          "|".join(self.functions)))

    def print_original_hit(self):
        """Returns string representation of DIAMOND hit in tabular format"""
        return "\t".join((self.query_id, self.subject_id, str(self.identity),
                          str(self.length), str(self.mismatch), str(self.s_len),
                          str(self.q_start), str(self.q_end), str(self.s_start),
                          str(self.s_end), str(self.evalue), str(self.bitscore)))
