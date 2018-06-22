#!/usr/bin/python
import os, csv, operator
import unittest
from collections import Counter,defaultdict

from context import Fama
from Fama.DiamondParser.DiamondParser import DiamondParser
from Fama.OutputUtil.JSONUtil import import_annotated_reads
from Fama.DiamondParser.hit_utils import cleanup_protein_id
from Fama.ReferenceLibrary.TaxonomyData import TaxonomyData
from Fama.OutputUtil.Report import generate_report
from Fama.TaxonomyProfile import TaxonomyProfile
from Fama.OutputUtil.KronaXMLWriter import generate_taxonomy_chart
from Fama.Project import Project

data_dir = 'data'
config_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'config.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project.ini')
#sample = 'test_sample'
end = 'pe1'

project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_Hans_sulfate_t.ini')
sample = 'HL1A'

class DiamondParserTest(unittest.TestCase):

    def setUp(self):
        self.parser = DiamondParser(config_file=config_path, project_file=project_path, sample=sample, end=end)

    @unittest.skip("for faster testing")
    def test_1_load_taxdata(self):
        tax_data = TaxonomyData(self.parser.config)
        tax_data.load_taxdata(self.parser.config)
        self.assertEqual(len(tax_data.names), 1721025)
        self.assertEqual(len(tax_data.nodes), 1721025)

    @unittest.skip("for faster testing")
    def test_2_get_taxonomy_profile(self):
        tax_data = TaxonomyData(self.parser.config)
        tax_data.load_taxdata(self.parser.config)
        self.parser.reads = import_annotated_reads(os.path.join(self.parser.project.get_project_dir(self.parser.sample), self.parser.sample + '_' + self.parser.end + '_' + self.parser.project.get_reads_json_name()))

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

    @unittest.skip("for faster testing")
    def test_3_build_taxonomy_profile(self):
        tax_data = TaxonomyData(self.parser.config)
        tax_data.load_taxdata(self.parser.config)
        self.parser.reads = import_annotated_reads(os.path.join(self.parser.project.get_project_dir(self.parser.sample), self.parser.sample + '_' + self.parser.end + '_' + self.parser.project.get_reads_json_name()))
        scores = defaultdict(lambda : defaultdict(float))
        print(sample, end, 'read count', str(len(self.parser.reads)))
        for read in self.parser.reads.keys():
            if self.parser.reads[read].get_status() == 'function,besthit' or self.parser.reads[read].get_status() == 'function':
                hits = self.parser.reads[read].get_hit_list().get_hits()
                for hit in hits:
                    protein_taxid = self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                    scores[protein_taxid]['count'] += 1.0
                    scores[protein_taxid]['identity'] += hit.get_identity()
                if len(hits) == 1:
                    read_functions = self.parser.reads[read].get_functions()
                    for function in read_functions:
                        scores[self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hits[0].get_subject_id()))]['rpkm'] += read_functions[function]
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
                            scores[taxid]['rpkm'] += read_functions[protein_taxids[taxid]]

        print(sample, end, 'tax id count', str(len(scores)))
        tax_profile = TaxonomyProfile()
        #print(scores)
        tax_profile.build_taxonomy_profile(tax_data,scores)
        print ('Root children:',','.join(tax_profile.tree.data['1'].children))
        if '2157' in tax_profile.tree.data:
            print ('Archaea found') 
        if '10239' in tax_profile.tree.data:
            print ('Viruses found') 
        generate_taxonomy_chart(tax_profile, self.parser.sample, 'test.xml')
        #print(tax_profile.print_taxonomy_profile())
        #for taxid in tax_profile.tree.data:
        #    print(taxid, tax_profile.tree.data[taxid].name, tax_profile.tree.data[taxid].rank, tax_profile.tree.data[taxid].get_attribute('score'))
        
        #self.assertEqual(len(tax_profile.tree.data), 47)
        self.assertTrue(tax_profile.tree.data)
    
    def test_4_build_taxonomy_profile(self):
        tax_data = TaxonomyData(self.parser.config)
        tax_data.load_taxdata(self.parser.config)
        
        project = Project(config_file=config_path, project_file=project_path)
        project.load_functional_profile()

        for sample in sorted(project.samples.keys()):
            for end in sorted(project.samples[sample].keys()):

                reads = project.samples[sample][end]
                #print(sample, end, 'read count', str(len(reads)))
                scores = defaultdict(lambda : defaultdict(float))
                for read in reads:
                    if reads[read].get_status() == 'function,besthit' or reads[read].get_status() == 'function':
                        hits = reads[read].get_hit_list().get_hits()
                        for hit in hits:
                            protein_taxid = project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                            scores[protein_taxid]['count'] += 1.0
                            scores[protein_taxid]['identity'] += hit.get_identity()
                        if len(hits) == 1:
                            read_functions = reads[read].get_functions()
                            for function in read_functions:
                                scores[project.ref_data.lookup_protein_tax(cleanup_protein_id(hits[0].get_subject_id()))]['rpkm'] += read_functions[function]
                        else:
                            read_functions = reads[read].get_functions()
                            protein_taxids = {}
                            for hit in hits:
                                hit_taxid = project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                                hit_functions = hit.get_functions()
                                for hit_function in hit_functions:
                                    protein_taxids[hit_taxid] = hit_function
                            for taxid in protein_taxids:
                                if protein_taxids[taxid] in read_functions:
                                    scores[taxid]['rpkm'] += read_functions[protein_taxids[taxid]]
                print(sample, end, 'tax id count', str(len(scores)))
                tax_profile = TaxonomyProfile()
                outfile = os.path.join(project.options.get_project_dir(sample), project.options.get_output_subdir(sample), sample + '_' + end + '_' + 'taxonomy_profile.xml')
                tax_profile.build_taxonomy_profile(tax_data, scores)
                #print ('Root children:',','.join(tax_profile.tree.data['1'].children))
                #if '2157' in tax_profile.tree.data:
                #    print ('Archaea found') 
                #if '10239' in tax_profile.tree.data:
                #    print ('Viruses found') 
                #print(tax_profile.print_taxonomy_profile())
                generate_taxonomy_chart(tax_profile, sample, outfile)
                self.assertTrue(tax_profile.tree.data)


    def tearDown(self):
        self.parser = None

if __name__=='__main__':
    unittest.main()
