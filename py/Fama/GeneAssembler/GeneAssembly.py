import os, csv
from subprocess import Popen, PIPE, CalledProcessError
from collections import defaultdict, Counter
from Fama.GeneAssembler.Contig import Contig
from Fama.GeneAssembler.Gene import Gene
from Fama.DiamondParser.hit_utils import autovivify

class GeneAssembly:
    def __init__(self):
        self.reads = autovivify(2) # Stores read IDs and sample IDs: self.reads[function_id][read_id] = sample_id
        self.contigs = autovivify(2, Contig)# Stores Contig objects: self.contigs[function_id][contig_id] = Contig
    
    def calculate_mean_coverage(self, sample=None):
        contig_count = 0
        coverage_sum = 0.0
        for function_id in self.contigs.keys():
            for contig_id in self.contigs[function_id].keys():
                contig_count += 1
                contig_coverage = 0.0
                if sample:
                    contig_coverage = self.contigs[function_id][contig_id].get_coverage(sample)
                else:
                    contig_coverage = self.contigs[function_id][contig_id].get_coverage()
                coverage_sum += contig_coverage
        return coverage_sum/contig_count
    
