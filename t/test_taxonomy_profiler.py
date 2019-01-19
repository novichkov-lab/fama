#!/usr/bin/python3
import os, csv, operator
import unittest
from collections import Counter,defaultdict

from context import Fama
from Fama.utils import autovivify,cleanup_protein_id,sanitize_file_name

from Fama.DiamondParser.DiamondParser import DiamondParser
from Fama.OutputUtil.JSONUtil import import_annotated_reads
from Fama.OutputUtil.Report import generate_fastq_report
from Fama.TaxonomyProfile import TaxonomyProfile
from Fama.OutputUtil.KronaXMLWriter import generate_taxonomy_chart
from Fama.OutputUtil.KronaXMLWriter import generate_loose_taxonomy_chart
from Fama.OutputUtil.KronaXMLWriter import generate_taxonomy_series_chart
from Fama.Project import Project

data_dir = 'data'
config_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'config.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project.ini')
sample_id = 'test_sample'
end = 'pe1'


#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_Hans_sulfate_t.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_EB271_nitrogen_t.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW306_cazy_t.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW3062M_universal1.ini')
sample_id = 'sample6'

class TaxonomyProfilingTest(unittest.TestCase):

    def setUp(self):
        self.project = Project(config_file=config_path, project_file=project_path)
        self.project.load_project()
        self.parser = DiamondParser(config = self.project.config, 
                            options=self.project.options, 
                            taxonomy_data=self.project.taxonomy_data,
                            ref_data=self.project.ref_data,
                            sample=self.project.samples[sample_id], 
                            end=end)

#    @unittest.skip("for faster testing")
    def test_1_load_taxdata(self):
        self.assertEqual(len(self.project.taxonomy_data.names), 1770528)
        self.assertEqual(len(self.project.taxonomy_data.nodes), 1770528)

#    @unittest.skip("for faster testing")
    def test_2_get_taxonomy_profile(self):
        self.parser.reads = import_annotated_reads(os.path.join(self.parser.options.get_project_dir(self.parser.sample.sample_id), self.parser.sample.sample_id + '_' + self.parser.end + '_' + self.parser.options.get_reads_json_name()))

        for read in self.parser.reads:
            print(self.parser.reads[read].get_status())
        generate_fastq_report(self.parser)
        tax_stats = Counter()
        identity_stats = defaultdict(float)
        rpkm_stats = defaultdict(float)
        for read in self.parser.reads.keys():
            if self.parser.reads[read].get_status() == 'function':
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
        counts_per_rank, identity_per_rank, rpkm_per_rank = self.project.taxonomy_data.get_taxonomy_profile(tax_stats, identity_stats, rpkm_stats)
        print (counts_per_rank)
        self.assertEqual(len(counts_per_rank), 7)
        self.assertEqual(counts_per_rank['superkingdom']['Bacteria'], 7)

#    @unittest.skip("for faster testing")
    def test_3_build_taxonomy_profile(self):
        self.parser.reads = import_annotated_reads(os.path.join(self.parser.options.get_project_dir(self.parser.sample.sample_id), self.parser.sample.sample_id + '_' + self.parser.end + '_' + self.parser.options.get_reads_json_name()))
        scores = defaultdict(lambda : defaultdict(float))
        print(sample_id, end, 'read count', str(len(self.parser.reads)))
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

        print(sample_id, end, 'tax id count', str(len(scores)))
        tax_profile = TaxonomyProfile()
        #print(scores)
        tax_profile.build_taxonomy_profile(self.project.taxonomy_data,scores)
        print ('Root children:',','.join(tax_profile.tree.data['1'].children))
        if '2157' in tax_profile.tree.data:
            print ('Archaea found') 
        if '10239' in tax_profile.tree.data:
            print ('Viruses found') 
        generate_taxonomy_chart(tax_profile, self.parser.sample.sample_id, 'test.xml')
        #print(tax_profile.print_taxonomy_profile())
        #for taxid in tax_profile.tree.data:
        #    print(taxid, tax_profile.tree.data[taxid].name, tax_profile.tree.data[taxid].rank, tax_profile.tree.data[taxid].get_attribute('score'))
        
        #self.assertEqual(len(tax_profile.tree.data), 47)
        self.assertTrue(tax_profile.tree.data)
    
#    @unittest.skip("for faster testing")
    def test_4_build_taxonomy_profile(self):
        
        for sample in self.project.list_samples():
            self.project.import_reads_json(sample,self.project.ENDS)
            for end in self.project.ENDS:
                
                #print(sample, end, 'read count', str(len(reads)))
                scores = defaultdict(lambda : defaultdict(float))
                for read_id,read in self.project.samples[sample_id].reads[end].items():
                    
                    if read.get_status() == 'function,besthit' or read.get_status() == 'function':
                        hits = read.get_hit_list().get_hits()
                        for hit in hits:
                            protein_taxid = self.project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                            scores[protein_taxid]['count'] += 1.0
                            scores[protein_taxid]['identity'] += hit.get_identity()
                        if len(hits) == 1:
                            read_functions = read.get_functions()
                            for function in read_functions:
                                scores[self.project.ref_data.lookup_protein_tax(cleanup_protein_id(hits[0].get_subject_id()))]['rpkm'] += read_functions[function]
                        else:
                            read_functions = read.get_functions()
                            protein_taxids = {}
                            for hit in hits:
                                hit_taxid = self.project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                                hit_functions = hit.get_functions()
                                for hit_function in hit_functions:
                                    protein_taxids[hit_taxid] = hit_function
                            for taxid in protein_taxids:
                                if protein_taxids[taxid] in read_functions:
                                    scores[taxid]['rpkm'] += read_functions[protein_taxids[taxid]]
                print(sample, end, 'tax id count', str(len(scores)))
                tax_profile = TaxonomyProfile()
                outfile = os.path.join(self.project.options.get_project_dir(sample), self.project.options.get_output_subdir(sample), sample + '_' + end + '_' + 'taxonomy_profile.xml')
                tax_profile.build_taxonomy_profile(self.project.taxonomy_data, scores)
                #print ('Root children:',','.join(tax_profile.tree.data['1'].children))
                #if '2157' in tax_profile.tree.data:
                #    print ('Archaea found') 
                #if '10239' in tax_profile.tree.data:
                #    print ('Viruses found') 
                #print(tax_profile.print_taxonomy_profile())
                generate_taxonomy_chart(tax_profile, sample, outfile)
                self.assertTrue(tax_profile.tree.data)
                
            self.project.samples[sample].reads[end] = None

#    @unittest.skip("for faster testing")
    def test_5_build_multifile_taxonomy_profile(self):
        scores = defaultdict(lambda : defaultdict(dict))

        for sample in self.project.list_samples():
            self.project.import_reads_json(sample,self.project.ENDS)
            for end in self.project.ENDS:
                scaling_factor = 1.0
                if end == 'pe1':
                    scaling_factor = self.project.options.get_fastq1_readcount(sample)/(self.project.options.get_fastq1_readcount(sample) + self.project.options.get_fastq2_readcount(sample))
                elif end == 'pe2':
                    scaling_factor = self.project.options.get_fastq2_readcount(sample)/(self.project.options.get_fastq1_readcount(sample) + self.project.options.get_fastq2_readcount(sample))
                else:
                    raise Exception('Unknown end identifier')

                multiple_hits = 0
                read_count = 0
                
                for read_id in self.project.samples[sample].reads[end]:
                    read = self.project.samples[sample].reads[end][read_id]
                    if read.get_status() == 'function,besthit' or read.get_status() == 'function':
                        
                        read_functions = read.get_functions()
                        hits = read.get_hit_list().get_hits()
                        if len(hits) >1:
                            multiple_hits += 1
                        read_count += 1
                        # Count hits and identity
                        for hit in hits:
                            hit_taxid = self.project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                            if sample in scores[hit_taxid]:
                                scores[hit_taxid][sample]['count'] += 1.0/len(hits)
                                scores[hit_taxid][sample]['hit_count'] += 1.0
                                scores[hit_taxid][sample]['identity'] += hit.get_identity()
                            else:
                                scores[hit_taxid][sample]['count'] = 1.0/len(hits)
                                scores[hit_taxid][sample]['hit_count'] = 1.0 # we need hit count to calculate average identity
                                scores[hit_taxid][sample]['identity'] = hit.get_identity()
                                # Initialize 'rpkm' here
                                scores[hit_taxid][sample]['rpkm'] = 0.0

                            
                        # Count RPKM
                        # If we have only one hit, all RPKM scores of the read would be assigned to the tax id of the hit
                        # If we have more than one hit, RPKM scores would be equally divided between tax ids of the hits, with regard to functional assignments of the hits
                        function_taxids = defaultdict(list)
                        # First, collect all taxonomy IDs for each function assigned to the hit
                        for function in read_functions:
                            for hit in hits:
                                hit_taxid = self.project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                                hit_functions = hit.get_functions()
                                for hit_function in hit.get_functions():
                                    if hit_function == function:
                                        function_taxids[function].append(hit_taxid)
                        # Second, for each tax ID add its share of the RPKM score assigned to hit's function
                        for function in function_taxids:
                            tax_count = len(function_taxids[function])
                            for hit_taxid in function_taxids[function]:
                                scores[hit_taxid][sample]['rpkm'] += read_functions[function] * scaling_factor/tax_count
                print(sample, end, 'read count', str(read_count))
                print ('Reads with multiple hits: ', multiple_hits)
                
                self.project.samples[sample].reads[end] = None
                
        # Now, we have all taxonomy ids, from each samples, in a single nested dictionary. 
        print('tax id count', str(len(scores)))

        tax_profile = TaxonomyProfile()
        outfile = os.path.join(self.project.options.get_work_dir(), self.project.options.get_name() + '_taxonomy_profile.xml')
        outfile = outfile.replace(' ', '_')
        outfile = outfile.replace("'", "")

        tax_profile.build_functional_taxonomy_profile(self.project.taxonomy_data, scores)

        generate_taxonomy_series_chart(tax_profile, sorted(self.project.samples.keys()), outfile)
        self.assertTrue(tax_profile.tree.data)

#    @unittest.skip("for faster testing")
    def test_6_build_project_taxonomy_profile_1function(self):
        
        target_function = 'GH9'
        
        
        scores = defaultdict(lambda : defaultdict(dict))
        
        for sample in self.project.list_samples():
            self.project.import_reads_json(sample,self.project.ENDS)
            for end in self.project.ENDS:
                
                
                scaling_factor = 1.0
                if end == 'pe1':
                    scaling_factor = self.project.options.get_fastq1_readcount(sample)/(self.project.options.get_fastq1_readcount(sample) + self.project.options.get_fastq2_readcount(sample))
                elif end == 'pe2':
                    scaling_factor = self.project.options.get_fastq2_readcount(sample)/(self.project.options.get_fastq1_readcount(sample) + self.project.options.get_fastq2_readcount(sample))
                else:
                    raise Exception('Unknown end identifier')

                multiple_hits = 0
                read_count = 0
                
                for read_id,read in self.project.samples[sample].reads[end].items():
                    #read = self.project.samples[sample].reads[end][read_id]
                    if read.get_status() == 'function,besthit' or read.get_status() == 'function':
                        
                        read_functions = read.get_functions()
                        
                        if target_function not in read_functions:
                            continue
                            
                        hits = read.get_hit_list().get_hits()
                        if len(hits) >1:
                            multiple_hits += 1
                        read_count += 1
                        # Count hits and identity
                        for hit in hits:
                            hit_taxid = self.project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                            if sample in scores[hit_taxid]:
                                scores[hit_taxid][sample]['count'] += 1.0/len(hits)
                                scores[hit_taxid][sample]['hit_count'] += 1.0
                                scores[hit_taxid][sample]['identity'] += hit.get_identity()
                            else:
                                scores[hit_taxid][sample]['count'] = 1.0/len(hits)
                                scores[hit_taxid][sample]['hit_count'] = 1.0 # we need hit count to calculate average identity
                                scores[hit_taxid][sample]['identity'] = hit.get_identity()
                                # Initialize 'rpkm' here
                                scores[hit_taxid][sample]['rpkm'] = 0.0

                            
                        # Count RPKM
                        # If we have only one hit, all RPKM scores of the read would be assigned to the tax id of the hit
                        # If we have more than one hit, RPKM scores would be equally divided between tax ids of the hits, with regard to functional assignments of the hits
                        function_taxids = defaultdict(list)
                        # First, collect all taxonomy IDs for each function assigned to the hit
                        for function in read_functions:
                            if function != target_function:
                                continue
                            for hit in hits:
                                hit_taxid = self.project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                                hit_functions = hit.get_functions()
                                for hit_function in hit.get_functions():
                                    if hit_function == function:
                                        function_taxids[function].append(hit_taxid)
                        # Second, for each tax ID add its share of the RPKM score assigned to hit's function
                        for function in function_taxids:
                            if function != target_function:
                                continue
                            tax_count = len(function_taxids[function])
                            for hit_taxid in function_taxids[function]:
                                scores[hit_taxid][sample]['rpkm'] += read_functions[function] * scaling_factor/tax_count
                print(sample, end, 'read count', str(read_count))
                print ('Reads with multiple hits: ', multiple_hits)
                
                self.project.samples[sample].reads[end] = None
                
        # Now, we have all taxonomy ids, from each samples, in a single nested dictionary. 
        print('tax id count', str(len(scores)))

        tax_profile = TaxonomyProfile()
        outfile = os.path.join(self.project.options.get_work_dir(), self.project.options.get_name() + '_' + target_function + '_taxonomy_profile.xml')
        outfile = outfile.replace(' ', '_')
        outfile = outfile.replace("'", "")

        tax_profile.build_functional_taxonomy_profile(self.project.taxonomy_data, scores)

        generate_taxonomy_series_chart(tax_profile, sorted(self.project.list_samples()), outfile)
        self.assertTrue(tax_profile.tree.data)


#    @unittest.skip("for faster testing")
    def test_7_build_loose_taxonomy_profile(self):
        metrics = 'rpkm'
        for sample in self.project.list_samples():
            if sample != sample_id:
                continue
            self.project.import_reads_json(sample,self.project.ENDS)
            
            norm_factor = 1000000/self.project.options.get_fastq1_readcount(sample_id)

            for end in self.project.ENDS:
                
                #print(sample, end, 'read count', str(len(reads)))
                scores = autovivify(2,float)
                for read_id,read in self.project.samples[sample_id].reads[end].items():
                    
                    if read.get_status() != 'function':
                        continue
                    
                    if read.taxonomy is None:
                        print ('No taxonomy ID assigned to ' + read_id)
                       # raise ValueError('No taxonomy ID assigned to ' + read_id)
                    else:
                        scores[read.taxonomy]['count'] += 1
                        scores[read.taxonomy][metrics] += norm_factor * sum(list(read.get_functions().values()))
                        scores[read.taxonomy]['hit_count'] += 1
                        scores[read.taxonomy]['identity'] += max([hit.get_identity() for hit in read.get_hit_list().get_hits()])

                print(sample, end, 'tax id count', str(len(scores)))
                tax_profile = TaxonomyProfile()
                outfile = os.path.join(self.project.options.get_project_dir(sample), self.project.options.get_output_subdir(sample), sample + '_' + end + '_' + metrics + '_taxonomy_profile.xml')
                tax_profile.build_taxonomy_profile(self.project.taxonomy_data, scores)
                #print ('Root children:',','.join(tax_profile.tree.data['1'].children))
                #if '2157' in tax_profile.tree.data:
                #    print ('Archaea found') 
                #if '10239' in tax_profile.tree.data:
                #    print ('Viruses found') 
                #print(tax_profile.print_taxonomy_profile())
                generate_loose_taxonomy_chart(tax_profile, sample, outfile, score=metrics)
                self.assertTrue(tax_profile.tree.data)
                
            self.project.samples[sample].reads[end] = None

#    @unittest.skip("for faster testing")
    def test_8_build_loose_project_taxonomy_profile(self):
        metrics = 'rpkm'
        scores = autovivify(3,float)
        for sample in self.project.list_samples():
            self.project.import_reads_json(sample,self.project.ENDS)
            
            norm_factor = 1000000/self.project.options.get_fastq1_readcount(sample)

            for end in self.project.ENDS:
                
                #print(sample, end, 'read count', str(len(reads)))
                for read_id,read in self.project.samples[sample].reads[end].items():
                    
                    if read.get_status() != 'function':
                        continue
                    
                    if read.taxonomy is None:
                        print ('No taxonomy ID assigned to ' + read_id)
                       # raise ValueError('No taxonomy ID assigned to ' + read_id)
                    else:
                        scores[read.taxonomy][sample]['count'] += 1
                        scores[read.taxonomy][sample][metrics] += norm_factor * sum(list(read.get_functions().values()))
                        scores[read.taxonomy][sample]['hit_count'] += 1
                        scores[read.taxonomy][sample]['identity'] += max([hit.get_identity() for hit in read.get_hit_list().get_hits()])

                print(sample, end, 'tax id count', str(len(scores)))
                
            self.project.samples[sample].reads[end] = None


        tax_profile = TaxonomyProfile()
        outfile = sanitize_file_name(os.path.join(self.project.options.get_work_dir(), self.project.options.get_name() + '_' + metrics + '_taxonomy_profile.xml'))

        tax_profile.build_functional_taxonomy_profile(self.project.taxonomy_data, scores)
        generate_taxonomy_series_chart(tax_profile, sorted(self.project.samples.keys()), outfile, score=metrics)
        self.assertTrue(tax_profile.tree.data)




    def tearDown(self):
        self.parser = None


if __name__=='__main__':
    unittest.main()
