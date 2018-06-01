from collections import Counter
from lib.DiamondParser.DiamondHit import DiamondHit
from lib.DiamondParser.DiamondHitList import DiamondHitList

class AnnotatedRead:
    def __init__(self, read_id):
        self.read_id = read_id  # 1st FASTQ line
        self.sequence = None    # 2nd FASTQ line
        self.quality = None     # 4th FASTQ line
        self.line3 = None       # 3rd FASTQ line
        self.hit_list = None    # hit_list is a DiamondHitList object
        self.status = 'unaccounted'
        self.functions = Counter()     # functions dictionary key is function ID and value is RPKM score
        
    def get_read_id(self):
        return self.read_id
        
    def set_hit_list(self,hit_list):
        self.hit_list = hit_list
        
    def get_hit_list(self):
        return self.hit_list
        
    def set_sequence(self, seq):
        self.sequence = seq
        
    def get_sequence(self):
        return self.sequence
        
    def set_quality(self, quality):
        self.quality = quality
    
    def get_quality(self):
        return self.quality
        
    def set_line3(self,line3):
        self.line3 = line3
        
    def get_line3(self):
        return self.line3
    
    def set_functions(self, functions):
        for function in functions:
            self.functions[function] += functions[function]

    def append_functions(self, functions):
        for function in functions:
            self.functions[function] += functions[function]

    def get_functions(self):
        return self.functions

    def set_status(self, status):
        if self.status == 'unaccounted' or self.status == 'nofunction':
            self.status = status
        elif self.status == 'function':
            if status != 'nofunction':
                self.status = status
        
    def get_status(self):
        return self.status

    def show_hits(self):
        self.hit_list.print_hits()
            
    def __str__(self):
        return 'Annotated read: ' + self.status + '\n'+ str(self.hit_list) + '\t' + ','.join(self.functions)
