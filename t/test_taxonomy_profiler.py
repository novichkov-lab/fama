#!/usr/bin/python
import os, csv, operator
import unittest
from collections import Counter,defaultdict

from context import Fama
from Fama.DiamondParser.DiamondParser import DiamondParser
from Fama.DiamondParser.hit_utils import cleanup_protein_id
from Fama.ReferenceLibrary.TaxonomyData import TaxonomyData
from Fama.OutputUtil.Report import generate_report

data_dir = 'data'
config_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'config.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project.ini')
sample = 'test_sample'
end = 'pe1'

class DiamondParserTest(unittest.TestCase):

    def setUp(self):
        self.parser = DiamondParser(config_file=config_path, project_file=project_path, sample=sample, end=end)

    def test_1_load_taxdata(self):
        tax_data = TaxonomyData(self.parser.config)
        tax_data.load_taxdata(self.parser.config)
        self.assertEqual(len(tax_data.names), 1721025)
        self.assertEqual(len(tax_data.nodes), 1721025)

    def test_2_get_taxonomy_profile(self):
        tax_data = TaxonomyData(self.parser.config)
        tax_data.load_taxdata(self.parser.config)
        self.parser.parse_background_output()
        for read in self.parser.reads:
            print(self.parser.reads[read].get_status())
        generate_report(self.parser)
        tax_stats = Counter()
        identity_stats = defaultdict(float)
        rpkm_stats = defaultdict(float)
        for read in self.parser.reads.keys():
            if self.parser.reads[read].get_status() == 'function,besthit':
                hits = self.parser.reads[read].get_hit_list().get_hits()
                for hit in hits:
                    protein_taxid = self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                    tax_stats[protein_taxid] += 1
                    identity_stats[protein_taxid] += hit.get_identity()
                if len(hits) == 1:
                    read_functions = self.parser.reads[read].get_functions()
                    for function in read_functions:
                        rpkm_stats[self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hits[0].get_subject_id()))] += read_functions[function]
                else:
                    read_functions = self.parser.reads[read].get_functions()
                    protein_taxids = {}
                    for hit in hits:
                        hit_taxid = self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                        hit_functions = hit.get_functions()
                        for hit_function in hit_functions:
                            protein_taxids[hit_taxid] = hit_function
                    for taxid in protein_taxids:
                        if protein_taxids[taxid] in read_functions:
                            rpkm_stats[taxid] += read_functions[protein_taxids[taxid]]

        print(tax_stats)
        counts_per_rank, identity_per_rank, rpkm_per_rank = tax_data.get_taxonomy_profile(tax_stats, identity_stats, rpkm_stats)
        print (counts_per_rank)
        self.assertEqual(len(counts_per_rank), 7)
        self.assertEqual(counts_per_rank['superkingdom']['Bacteria'], 7)


    def tearDown(self):
        self.parser = None

if __name__=='__main__':
    unittest.main()
