"""Describes GeneAssembly class"""
from lib.gene_assembler.contig import Contig
from lib.utils.utils import autovivify

class GeneAssembly(object):
    """GeneAssembly stores input reads and resulting contigs

    Attributes:
        reads (defaultdict[str,defaultdict[str,str]]): reads mapped to a
            certain function. Outer key is function identifier
            in function-specific assemblies (or 'Coassembly' for co-assembly),
            inner key is read idenrifier, value is sample identifier
        contigs (defaultdict[str,defaultdict[str,:obj:Contig]]): assembled
            contigs. Outer key is function identifier in function-specific
            assemblies (or 'Coassembly'), inner key is read idenrifier,
            value is Contig object
    """
    def __init__(self):
        self.reads = autovivify(2)
        self.contigs = autovivify(2, Contig)

    def calculate_average_coverage(self, sample=None):
        """Returns average read coverage of the assembly for a given sample
        or for all reads

        Args:
            sample (str, optional): sample identifier
        """
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
