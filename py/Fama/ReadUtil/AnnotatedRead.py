from collections import defaultdict
from Fama.DiamondParser.DiamondHit import DiamondHit
from Fama.DiamondParser.DiamondHitList import DiamondHitList

class AnnotatedRead:
    def __init__(self, read_id = None):
        self.read_id = read_id  
        # Please note that read ID may not contain
        # entire ID line from FASTQ file. It contains only part of the 
        # line before the first space symbol. 
        self.read_id_line = None  # 1st FASTQ line
        self.sequence = None    # 2nd FASTQ line
        self.quality = None     # 4th FASTQ line
        self.line3 = None       # 3rd FASTQ line
        self.hit_list = None    # hit_list is a DiamondHitList object
        self.status = 'unaccounted'
        self.functions = defaultdict(float)     # functions dictionary key is function ID and value is RPKM score
        self.pe_id = None          # 1st FASTQ line
        self.pe_sequence = None    # 2nd FASTQ line
        self.pe_quality = None     # 4th FASTQ line
        self.pe_line3 = None       # 3rd FASTQ line
        self.taxonomy = None       # NCBI Taxonomy ID set by LCA algorithm
        
    #~ def get_read_id(self):
        #~ return self.read_id

    #~ def set_read_id(self,read_id):
        #~ self.read_id = read_id
        
    #~ def set_hit_list(self,hit_list):
        #~ self.hit_list = hit_list
        
    #~ def get_hit_list(self):
        #~ return self.hit_list
        
    #~ # Sequence data
    #~ def set_read_id_line(self,line):
        #~ self.read_id_line = line
        
    #~ def get_read_id_line(self):
        #~ return self.read_id_line
        
    #~ def set_sequence(self, seq):
        #~ self.sequence = seq
        
    #~ def get_sequence(self):
        #~ return self.sequence
        
    #~ def set_quality(self, quality):
        #~ self.quality = quality
    
    #~ def get_quality(self):
        #~ return self.quality
        
    #~ def set_line3(self,line3):
        #~ self.line3 = line3
        
    #~ def get_line3(self):
        #~ return self.line3

    #~ # Paired-end data
    #~ def set_pe_id(self, pe_id):
        #~ self.pe_id = pe_id
        
    #~ def set_pe_sequence(self, seq):
        #~ self.pe_sequence = seq
        
    #~ def set_pe_quality(self, quality):
        #~ self.pe_quality = quality
    
    #~ def set_pe_line3(self,line3):
        #~ self.pe_line3 = line3

    # Function data
    def set_functions(self, functions):
        for function in functions:
            self.functions[function] += functions[function]

    def append_functions(self, functions):
        for function in functions:
            self.functions[function] += functions[function]

    #~ def get_functions(self):
        #~ return self.functions

    def set_status(self, status):
        if status in ['unaccounted', 'nofunction', 'function']:
            self.status = status
        else:
            raise ValueError('Unknown read status: ' + status)
        
    #~ def get_status(self):
        #~ return self.status

    def show_hits(self):
        self.hit_list.print_hits()
            
    def __str__(self):
        return 'Annotated read: ' + self.status + '\n'+ str(self.hit_list) + '\t' + ','.join(self.functions)
