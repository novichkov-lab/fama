#!/usr/bin/python
from collections import defaultdict

class Gene:
    def __init__(self, contig_id = "", gene_id = "", sequence = None):
        self.contig_id = contig_id
        self.gene_id = gene_id
        self.protein_sequence = sequence
        self.tax_id = None
        self.functions = defaultdict(float)     # functions dictionary key is function ID and value is RPKM score
        self.hit_list = None
        self.status = 'unaccounted'


    def set_protein_sequence(self, sequence):
        self.protein_sequence = sequence

    # Function data
    def set_functions(self, functions):
        for function in functions:
            self.functions[function] += functions[function]

    def set_hit_list(self,hit_list):
        self.hit_list = hit_list

    def set_status(self, status):
        if self.status == 'unaccounted' or self.status == 'nofunction':
            self.status = status
        elif self.status == 'function':
            if status != 'nofunction':
                self.status = status
        
    def get_status(self):
        return self.status

    def set_taxonomy_id(self, tax_id):
        self.tax_id = tax_id

    def get_taxonomy_id(self):
        return self.tax_id
        
