#!/usr/bin/python
from collections import defaultdict

class Gene:
    def __init__(self, contig_id = "", gene_id = "", sequence = None, start = None, end = None, strand = None):
        self.contig_id = contig_id
        self.gene_id = gene_id
        self.protein_sequence = sequence
        self.start = start
        self.end = end
        self.strand = strand
        self.tax_id = None # Obsolete: Taxonomy from UniProt search
        self.functions = defaultdict(float)     # functions dictionary key is function ID and value is RPKM score
        self.hit_list = None
        self.uniref_hit = None
        self.status = 'unaccounted'
        self.taxonomy = None       # NCBI Taxonomy ID set by LCA algorithm


    # Function data
    def set_functions(self, functions):
        for function in functions:
            self.functions[function] += functions[function]

    def set_uniref_hit(self,hit):
        self.uniref_hit = hit

    def set_status(self, status):
        if self.status == 'unaccounted' or self.status == 'nofunction':
            self.status = status
        elif self.status == 'function':
            if status != 'nofunction':
                self.status = status
        
        
