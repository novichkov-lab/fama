#!/usr/bin/python
import os, csv, operator
import unittest
from context import lib
from collections import Counter

from lib.DiamondParser.DiamondParser import DiamondParser
from lib.DiamondParser.DiamondParser import cleanup_protein_id
from lib.ReferenceLibrary.TaxonomyData import TaxonomyData


data_dir = 'data'
config_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'config.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project.ini')
sample = 'test_sample'
end = 'pe1'

class DiamondParserTest(unittest.TestCase):

    def setUp(self):
        self.parser = DiamondParser(config_path, project_path, sample, end)

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
        from lib.OutputUtil.Report import generate_report
        generate_report(self.parser)
        tax_stats = Counter()
        identity_stats = Counter()
        for read in self.parser.reads.keys():
            if self.parser.reads[read].get_status() == 'function,besthit':
                print('37: ',read)
                for hit in self.parser.reads[read].get_hit_list().get_hits():
                    tax_stats[self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))] += 1
                    identity_stats[self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))] += hit.get_identity()
        
        print(tax_stats)
        counts_per_rank, identity_per_rank = tax_data.get_taxonomy_profile(tax_stats, identity_stats)
        print (counts_per_rank)
        self.assertEqual(len(counts_per_rank), 7)
        self.assertEqual(counts_per_rank['superkingdom']['Bacteria'], 6)


    def tearDown(self):
        self.parser = None

if __name__=='__main__':
    unittest.main()
