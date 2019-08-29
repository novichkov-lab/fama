import re

class DiamondHit:
    def __init__(self, query_id = "", subject_id = "", identity = 0.0,\
                length = 0, mismatch = 0, s_len = 0, q_start = 0, \
                q_end = 0, s_start = 0, s_end = 0, evalue = 0.0, \
                bitscore = 0.0):
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
        # Takes list of tabular output fields
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
        except ValueError as detail:
            print ('Value parsing failed for read ' + tabular_output_fields[0])
            print (str(detail))
            raise
        except IndexError as e:
            print ('Wrong number of fields in tabular output for read ' + tabular_output_fields[0])
            print (str(e))
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
        # Takes list of tabular output fields
        try:
            entry_tokens[2] = float(entry_tokens[2])
            entry_tokens[3] = int(entry_tokens[3])
            entry_tokens[4] = int(entry_tokens[4])
            entry_tokens[5] = int(entry_tokens[5])
            entry_tokens[6] = int(entry_tokens[6])
            entry_tokens[7] = int(entry_tokens[7])
            entry_tokens[8] = int(entry_tokens[8])
            entry_tokens[9] = int(entry_tokens[9])
            entry_tokens[10] = float(entry_tokens[10])
            entry_tokens[11] = float(entry_tokens[11])
        except ValueError as detail:
            print ('Value parsing failed for read ' + entry_tokens[0])
            print (str(detail))
            raise
        except IndexError as e:
            print ('Wrong number of fields in tabular output for read ' + entry_tokens[0])
            print (str(e))
            raise
        self.query_id = entry_tokens[0]
        self.subject_id = entry_tokens[1]
        self.identity = entry_tokens[2]
        self.length = entry_tokens[3]
        self.mismatch = entry_tokens[4]
        self.s_len = entry_tokens[5]
        self.q_start = entry_tokens[6]
        self.q_end = entry_tokens[7]
        self.s_start = entry_tokens[8]
        self.s_end = entry_tokens[9]
        self.evalue = entry_tokens[10]
        self.bitscore = entry_tokens[11]
        try:
            self.functions = entry_tokens[12].split('|')
        except IndexError:
            print ('Unable to parse function list:',entry_tokens)
    
    def annotate_hit(self, ref_data):
        if self.subject_id:
            self.functions = ref_data.lookup_protein_function(self.subject_id)

    def __str__(self):
        return "\t".join((self.query_id, self.subject_id, str(self.identity), \
                          str(self.length), str(self.mismatch), str(self.s_len), str(self.q_start), \
                          str(self.q_end), str(self.s_start), str(self.s_end), str(self.evalue), str(self.bitscore), "|".join(self.functions)))

    def print_original_hit(self):
        return "\t".join((self.query_id, self.subject_id, str(self.identity), \
                          str(self.length), str(self.mismatch), str(self.s_len), str(self.q_start), \
                          str(self.q_end), str(self.s_start), str(self.s_end), str(self.evalue), str(self.bitscore)))
