#!/usr/bin/python3
import os, csv, operator
import unittest
import json
from collections import Counter

from context import lib

from lib.utils.const import ENDS

from lib.project.project import Project
from lib.project.sample import Sample
from lib.project.project_options import ProjectOptions

from lib.diamond_parser.diamond_hit import DiamondHit
from lib.diamond_parser.diamond_hit_list import DiamondHitList
from lib.diamond_parser.diamond_parser import DiamondParser
from lib.diamond_parser.hit_utils import compare_protein_hits_lca

from lib.sequences.annotated_read import AnnotatedRead
from lib.output.json_util import import_annotated_reads,export_annotated_reads

from lib.protein_functional_pipeline import functional_profiling_pipeline, parse_background_output


data_dir = 'data'
config_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'config_test.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'p', 'project_protein_ht1_isolates_nitrogen10.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'p', 'project_protein_enigma_nonpseudomonas_universal1.4.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'p', 'project_protein_isolates_universal1.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'p', 'project_protein_enigma_isolates_nitrogen10.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'p', 'project_protein_cami_ref_universal_v1.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'projects', 'protein_project.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'p', 'project_protein_test.ini')
sample = 'test_sample1'
#sample = 'GW821-FHT04F04'
end = 'pe1'

class DiamondParserTest(unittest.TestCase):

    def setUp(self):
        #self.parser = DiamondParser(config_file=config_path, project_file=project_path, sample=sample, end=end)
        self.project = Project(config_file=config_path, project_file=project_path)
        for sample_id in self.project.list_samples():
            sample = Sample(sample_id=sample_id)
            sample.load_sample(self.project.options)
            self.project.samples[sample_id] = sample
        
    def test_1_test_pipeline(self):
        #parser = functional_profiling_pipeline(config_file=config_path, project_file=project_path, sample=sample, end=end)
        pass

    def test_2_load_proteins(self):
        self.project.import_reads_json(sample, ENDS)
        protein = 'D16-4706_contig_11213_8'
        self.assertTrue(protein in self.project.samples[sample].reads[end])
        print(self.project.samples[sample].reads[end][protein].taxonomy)
        

    def test_3_protein_taxonomy(self):
        self.project.import_reads_json(sample, ENDS)
        protein = 'D16-4706_contig_11213_7'
        print('D16-4706_contig_11213_7 taxonomy')
        print(self.project.samples[sample].reads[end][protein].taxonomy)
        
        parser = DiamondParser(config=self.project.config,
                       options=self.project.options,
                       taxonomy_data=self.project.taxonomy_data,
                       ref_data=self.project.ref_data,
                       sample=self.project.samples[sample],
                       end=end)
        parser.parse_reference_output()
        print(str(parser.reads[protein]))
        
#        parse_background_output(parser)
        hit_line = 'D16-4706_contig_11213_7|4|257	fig|408672.3.peg.2637	63.0	254	94	256	1	254	2	255	1.1e-97	362.1'
        hit = DiamondHit()
        hit.create_hit(tabular_output_fields=hit_line.split('\t'))
        hit_list = DiamondHitList('D16-4706_contig_11213_7|4|257')
        hit_list.add_hit(hit)
        hit_list.annotate_hits(self.project.ref_data)
        hit_list.filter_list_by_identity(self.project.ref_data)
        print('hit_list')
        print(hit_list)
        
        compare_protein_hits_lca(parser.reads[protein], 4, 257, hit_list, 0.03, 1.0, 1.0, self.project.taxonomy_data, self.project.ref_data)
        print(parser.reads[protein].taxonomy)
        self.assertEqual(parser.reads[protein].taxonomy, '408672')
        

    def test_4_sample_taxonomy(self):
        with open('samples_taxonomy2.txt', 'w') as outfile:
            for sample_id in self.project.list_samples():
                self.project.import_reads_json(sample_id, ENDS)
                taxonomy_ids = []
                for protein_id, protein in self.project.samples[sample_id].reads[end].items():
                    taxonomy_ids.append(protein.taxonomy)
                lca_taxonomy = self.project.taxonomy_data.get_lca(taxonomy_ids)
                outfile.write('\t'.join([sample_id, lca_taxonomy, self.project.taxonomy_data.get_name(lca_taxonomy)]) + '\n')
        
        #~ self.project.import_reads_json(sample, ENDS)
        #~ for protein_id, protein in self.project.samples[sample].reads[end].items():
            #~ taxonomy_ids.append(protein.taxonomy)
        #~ lca_taxonomy = self.project.taxonomy_data.get_lca(taxonomy_ids)
        #~ print(sample, lca_taxonomy, self.project.taxonomy_data.get_name(lca_taxonomy))
        #~ self.assertEqual(lca_taxonomy, '28216')


    def test_5_export_proteins(self):
        with open('proteins.faa', 'w') as outfile:
            for sample_id in self.project.list_samples():
                self.project.import_reads_json(sample_id, ENDS)
                for protein_id, protein in self.project.samples[sample_id].reads[end].items():
                    if protein.status == 'function':
                        outfile.write('>' + protein_id + '|' + sample_id + '|' + ';'.join(protein.functions.keys()) +  '|' + self.project.taxonomy_data.get_name(protein.taxonomy) + '\n')
                        outfile.write(protein.sequence + '\n\n')

    def test_6_export_protein_table(self):
        out_file = os.path.join(self.project.options.work_dir, 'proteins.list.txt')
        with open(out_file, 'w') as outfile:
            for sample_id in self.project.list_samples():
                self.project.import_reads_json(sample_id, ENDS)
                for protein_id, protein in self.project.samples[sample_id].reads[end].items():
                    if protein.status == 'function':
                        protein_length = len(protein.sequence)
                        ref_length = int(protein.hit_list.hits[0].s_len)
                        outfile.write('\t'.join([
                            sample_id,
                            protein_id,
                            ';'.join(sorted(protein.functions.keys())),
                            '{0:.4f}'.format(protein_length/ref_length),
                            protein.taxonomy,
                            self.project.taxonomy_data.data[protein.taxonomy]['name']
                            ]) + '\n')
    def tearDown(self):
        self.parser = None


if __name__=='__main__':
    unittest.main()
