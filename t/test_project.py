#!/usr/bin/python3
import os, csv, operator
import unittest
import json
from collections import Counter,defaultdict

from context import Fama
from Fama.Project import Project
from Fama.Sample import Sample
from Fama.ProjectUtil.ProjectOptions import ProjectOptions

config_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'config.ini')
#config_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'config_rpL6_singleDB.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_EB271_nitrogen_t.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_Hans_sulfate_t.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW306_cazy_t.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW306_nitrogen_t.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_EB271_sulfate_t.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW3062M_rpl6_testdb.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW306_nitrogen9.ini')

def autovivify(levels=1, final=dict):
    return (defaultdict(final) if levels < 2 else
            defaultdict(lambda: autovivify(levels - 1, final)))

class ProjectTest(unittest.TestCase):

    def setUp(self):
        self.project = Project(config_file=config_path, project_file=project_path)
        for sample_id in self.project.list_samples():
            sample = Sample(sample_id=sample_id)
            sample.load_sample(self.project.options)
            self.project.samples[sample_id] = sample
    
    @unittest.skip("for faster testing")
    def test_project_options(self):
        print ('Print list of samples1')
        options = ProjectOptions(project_path)
        print (options.parser.sections())
        self.assertEqual(len(options.parser.sections()), 6)
        self.assertEqual(options.parser.sections()[0], 'sample1')
    
    @unittest.skip("for faster testing")
    def test_list_samples(self):
        print ('Print list of samples2')
        samples = self.project.list_samples()
        print (samples)
        self.assertEqual(len(samples), 6)
        self.assertEqual(samples[0], 'sample1')
        
    @unittest.skip("for faster testing")
    def test_check_project(self):
        print ('Print problems found in test project: ')
        self.project.check_project()

    @unittest.skip("for faster testing")
    def test_check_config(self):
        print ('Print problems found in test project: ')
        self.project.check_project()
        self.assertEqual(self.project.config.get_biscore_range_cutoff(self.project.options.get_collection()), 0.2)
        self.assertEqual(self.project.config.get_identity_cutoff(self.project.options.get_collection()), 40.0)
        ranks_cutoff = self.project.config.get_ranks_cutoffs(self.project.options.get_collection())
        
        print(ranks_cutoff)
        print(ranks_cutoff['species'])
    
    @unittest.skip("for faster testing")
    def test_load_project(self):
        print ('Load project from INI file')
        self.project = None
        self.project = Project(config_file=config_path, project_file=project_path)
        self.project.load_project()
        self.assertEqual(len(self.project.samples), 6)

    @unittest.skip("for faster testing")
    def test_collect_fragment_stats(self):
        print ('Load project from JSON')
        #self.project.load_functional_profile()
        
        sample_stats = autovivify(2, int)
        for sample_id in self.project.list_samples():
            if not self.project.samples[sample_id].is_paired_end:
                continue
            self.project.import_reads_json(sample_id, self.project.ENDS)
            print('Mapping data loaded for sample', sample_id)
            both_ends_mapped_reads = {}
            pe1_multiple_functions = {}
            pe2_multiple_functions = {}
            sample_stats['reads_pe1_total'][sample_id] = len(self.project.samples[sample_id].reads['pe1'])
            sample_stats['reads_pe2_total'][sample_id] = len(self.project.samples[sample_id].reads['pe2'])
            pe1_reads = self.project.samples[sample_id].reads['pe1']
            pe2_reads = self.project.samples[sample_id].reads['pe2']
            for read in pe1_reads:
                if pe1_reads[read].get_status() == 'function,besthit' or pe1_reads[read].get_status() == 'function':
                    if read in pe2_reads and (pe2_reads[read].get_status() == 'function,besthit' or pe2_reads[read].get_status() == 'function'):
                        sample_stats['both ends mapped'][sample_id] += 1
                        both_ends_mapped_reads[read] = 1
                    else:
                        sample_stats['pe1 mapped only'][sample_id] += 1
            for read in pe2_reads:
                if read not in both_ends_mapped_reads and (pe2_reads[read].get_status() == 'function,besthit' or pe2_reads[read].get_status() == 'function'):
                        sample_stats['pe2 mapped only'][sample_id] += 1

            for read in pe1_reads:
                if pe1_reads[read].get_status() == 'function,besthit' or pe1_reads[read].get_status() == 'function':
                    sample_stats['reads_pe1_mapped'][sample_id] += 1
                    if len(pe1_reads[read].get_functions()) == 1:
                        sample_stats['pe1 single function'][sample_id] += 1
                    elif len(pe1_reads[read].get_functions()) > 1:
                        sample_stats['pe1 multiple functions'][sample_id] += 1
                        pe1_multiple_functions[read] = 1

            for read in pe2_reads:
                if pe2_reads[read].get_status() == 'function,besthit' or pe2_reads[read].get_status() == 'function':
                    sample_stats['reads_pe2_mapped'][sample_id] += 1
                    if len(pe2_reads[read].get_functions()) == 1:
                        sample_stats['pe2 single function'][sample_id] += 1
                    elif len(pe2_reads[read].get_functions()) > 1:
                        sample_stats['pe2 multiple functions'][sample_id] += 1
                        pe2_multiple_functions[read] = 1
            
            for read in both_ends_mapped_reads:
                if len(pe1_reads[read].get_functions()) == 1:
                    sample_stats['pe1 single function, both ends mapped'][sample_id] += 1
                elif len(pe1_reads[read].get_functions()) > 1:
                    sample_stats['pe1 multiple functions, both ends mapped'][sample_id] += 1
                if len(pe2_reads[read].get_functions()) == 1:
                    sample_stats['pe2 single function, both ends mapped'][sample_id] += 1
                elif len(pe2_reads[read].get_functions()) > 1:
                    sample_stats['pe2 multiple functions, both ends mapped'][sample_id] += 1
                
            for read in pe1_reads:
                if pe1_reads[read].get_status() == 'function,besthit' or pe1_reads[read].get_status() == 'function':
                    hits = pe1_reads[read].get_hit_list()
                    if len(hits.get_hits()) == 1:
                        sample_stats['pe1 single hit'][sample_id] += 1
                    elif len(hits.get_hits()) > 1:
                        sample_stats['pe1 multiple hits'][sample_id] += 1

            for read in pe2_reads:
                if pe2_reads[read].get_status() == 'function,besthit' or pe2_reads[read].get_status() == 'function':
                    hits = pe2_reads[read].get_hit_list()
                    if len(hits.get_hits()) == 1:
                        sample_stats['pe2 single hit'][sample_id] += 1
                    elif len(hits.get_hits()) > 1:
                        sample_stats['pe2 multiple hits'][sample_id] += 1

            for read in pe1_multiple_functions:
                if pe1_reads[read].get_status() == 'function,besthit' or pe1_reads[read].get_status() == 'function':
                    hits = pe1_reads[read].get_hit_list()
                    if len(hits.get_hits()) == 1:
                        sample_stats['pe1 multiple functions, single hit'][sample_id] += 1
                    elif len(hits.get_hits()) > 1:
                        sample_stats['pe1 multiple functions, multiple hits'][sample_id] += 1
                
            for read in pe2_multiple_functions:
                if pe2_reads[read].get_status() == 'function,besthit' or pe2_reads[read].get_status() == 'function':
                    hits = pe2_reads[read].get_hit_list()
                    if len(hits.get_hits()) == 1:
                        sample_stats['pe2 multiple functions, single hit'][sample_id] += 1
                    elif len(hits.get_hits()) > 1:
                        sample_stats['pe2 multiple functions, multiple hits'][sample_id] += 1
            
            for read in both_ends_mapped_reads:
                hits = pe1_reads[read].get_hit_list()
                if len(hits.get_hits()) == 1:
                    sample_stats['pe1 single hit, both ends mapped'][sample_id] += 1
                elif len(hits.get_hits()) > 1:
                    sample_stats['pe1 multiple hits, both ends mapped'][sample_id] += 1
                hits = pe2_reads[read].get_hit_list()
                if len(hits.get_hits()) == 1:
                    sample_stats['pe2 single hit, both ends mapped'][sample_id] += 1
                elif len(hits.get_hits()) > 1:
                    sample_stats['pe2 multiple hits, both ends mapped'][sample_id] += 1
            
            self.project.samples[sample_id].reads = None
                
        with open('outfile.tsv', 'w') as of:
            for sample_id in self.project.list_samples():
                of.write('\t' + sample_id)
            of.write('\n')
            for item in sorted(sample_stats.keys()):
                of.write(item)
                for sample_id in self.project.list_samples():
                    if sample_id in sample_stats[item]:
                        of.write('\t' + str(sample_stats[item][sample_id]))
                    else:
                        of.write('\t0')
                of.write('\n')
            of.close()
                
        self.assertEqual(len(self.project.samples), 6)

    @unittest.skip("for faster testing")
    def test_top_size(self):
        print ('Load reads from JSON')
        sample_id = 'sample3'
        self.project.import_reads_json(sample_id,['pe1',])
        
        tsvfile = '/mnt/data3/FEBA/4703/nitrogen_v7.1_fama/sample3_pe1_bgr_tabular_output.txt'
        outfile = 'top_hit_count.txt'

        current_query_id = None
        top_size = 0
        identity_cutoff = 50.0
        length_cutoff = 15
        bitscore_range_cutoff = 0.97
        bitscore_cutoff = 0.0
        print ('Identity cutoff: ', identity_cutoff, ', Length cutoff: ', length_cutoff)

        with open(outfile, 'w') as of:
            with open(tsvfile, 'r', newline='') as f:
                tsvin = csv.reader(f, delimiter='\t')
                for row in tsvin:
                    if current_query_id is None:
                        current_query_id = row[0]
                        bitscore_cutoff = float(row[11]) * bitscore_range_cutoff
                    
                    # filtering by identity and length
                    if float(row[2]) < identity_cutoff:
                        continue # skip this line
                    if float(row[3]) < length_cutoff:
                        continue # skip this line

                    if row[0] != current_query_id:
                        read_id = current_query_id.split('|')[0]
                        of.write(read_id + '\t' + str(top_size) + '\t' + self.project.samples[sample_id].reads['pe1'][read_id].get_status() + '\n')
                        current_query_id = row[0]
                        top_size = 0
                        bitscore_cutoff = float(row[11]) * bitscore_range_cutoff

                    if float(row[11]) >= bitscore_cutoff:
                        top_size += 1
                f.closed
            of.closed

        
    def tearDown(self):
        self.project = None
        


if __name__=='__main__':
    unittest.main()
