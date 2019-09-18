#!/usr/bin/python3
import os
import unittest
from collections import Counter, defaultdict

from context import lib

from lib.utils.const import ENDS, STATUS_GOOD
from lib.utils.utils import autovivify, cleanup_protein_id, sanitize_file_name
from lib.diamond_parser.diamond_parser import DiamondParser
from lib.output.json_util import import_annotated_reads
from lib.output.report import generate_fastq_report, get_function_taxonomy_scores, get_function_scores
from lib.taxonomy.taxonomy_profile import TaxonomyProfile
from lib.output.krona_xml_writer import make_taxonomy_chart, make_taxonomy_series_chart, make_function_taxonomy_chart
from lib.project.project import Project

data_dir = 'data'
config_path = os.path.join(
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    'py',
    'config.ini'
)
#~ config_path = os.path.join(
    #~ os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    #~ 'py',
    #~ 'config_rpL6_singleDB.ini'
#~ )
project_path = os.path.join(
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    'py',
    'project.ini'
)
sample_id = 'test_sample'
end = 'pe1'


#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_Hans_sulfate_t.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_EB271_nitrogen_t.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW306_cazy_t.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW3062M_universal1.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'p', 'project_FW3062M_rpl6.ini')
#sample_id = 'FW306-4701'

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
    def test_2_get_hits_taxonomy_profile(self):
        self.parser.reads = import_annotated_reads(os.path.join(self.parser.options.get_project_dir(self.parser.sample.sample_id), self.parser.sample.sample_id + '_' + self.parser.end + '_' + self.parser.options.reads_json_name))
        generate_fastq_report(self.parser)
        tax_stats = Counter()
        identity_stats = defaultdict(float)
        rpkm_stats = defaultdict(float)
        for read_id, read in self.parser.reads.items():
            print(read.status, read.read_id, read.taxonomy, self.project.taxonomy_data.get_name(read.taxonomy))
            if read.status == STATUS_GOOD:
                hits = read.hit_list.hits
                for hit in read.hit_list.hits:
                    protein_taxid = self.parser.ref_data.lookup_protein_tax(hit.subject_id)
                    tax_stats[protein_taxid] += 1
                    identity_stats[protein_taxid] += hit.identity
                if len(hits) == 1:
                    for function in read.functions:
                        rpkm_stats[self.parser.ref_data.lookup_protein_tax(hits[0].subject_id)] += read.functions[function]
                else:
                    protein_taxids = {}
                    for hit in read.hit_list.hits:
                        hit_taxid = self.parser.ref_data.lookup_protein_tax(hit.subject_id)
                        hit_functions = hit.functions
                        for hit_function in hit_functions:
                            protein_taxids[hit_taxid] = hit_function
                    for taxid in protein_taxids:
                        if protein_taxids[taxid] in read.functions:
                            rpkm_stats[taxid] += read.functions[protein_taxids[taxid]]

        print(tax_stats)
        counts_per_rank, identity_per_rank, rpkm_per_rank = self.project.taxonomy_data.get_taxonomy_profile(tax_stats, identity_stats, rpkm_stats)
        print(counts_per_rank)
        print([self.project.taxonomy_data.get_name(taxid) + ':' + str(rpkm_stats[taxid]) for taxid in rpkm_stats])
        self.assertEqual(len(counts_per_rank), 7)
        self.assertEqual(counts_per_rank['superkingdom']['Bacteria'], 8)


#    @unittest.skip("for faster testing")
    def test_3_build_taxonomy_profile(self):
        self.parser.reads = import_annotated_reads(os.path.join(self.parser.options.get_project_dir(self.parser.sample.sample_id), self.parser.sample.sample_id + '_' + self.parser.end + '_' + self.parser.options.reads_json_name))
        scores = defaultdict(lambda : defaultdict(float))
        print(sample_id, end, 'read count', str(len(self.parser.reads)))
        for read_id, read in self.parser.reads.keys():
            if read.status == STATUS_GOOD:
                for hit in read.hit_list.hits:
                    protein_taxid = self.parser.ref_data.lookup_protein_tax(hit.subject_id)
                    scores[protein_taxid]['count'] += 1.0
                    scores[protein_taxid]['identity'] += hit.identity
                if len(hits) == 1:
                    for function in read.functions:
                        scores[self.parser.ref_data.lookup_protein_tax(hits[0].subject_id)]['rpkm'] += read.functions[function]
                else:
                    read_functions = self.parser.reads[read].functions
                    protein_taxids = {}
                    for hit in read.hit_list.hits:
                        hit_taxid = self.parser.ref_data.lookup_protein_tax(hit.subject_id)
                        for hit_function in hit.functions:
                            protein_taxids[hit_taxid] = hit_function
                    for taxid in protein_taxids:
                        if protein_taxids[taxid] in read.functions:
                            scores[taxid]['rpkm'] += read.functions[protein_taxids[taxid]]

        print(sample_id, end, 'tax id count', str(len(scores)))
        tax_profile = TaxonomyProfile()
        print(scores)
        tax_profile.make_taxonomy_profile(self.project.taxonomy_data,scores)
        print ('Root children:',','.join(tax_profile.tree.data['1'].children))
        if '2157' in tax_profile.tree.data:
            print ('Archaea found') 
        if '10239' in tax_profile.tree.data:
            print ('Viruses found') 
        make_taxonomy_chart(tax_profile, self.parser.sample.sample_id, 'test.xml', self.project.config.krona_path)
        print(tax_profile.print_taxonomy_profile())
        #for taxid in tax_profile.tree.data:
        #    print(taxid, tax_profile.tree.data[taxid].name, tax_profile.tree.data[taxid].rank, tax_profile.tree.data[taxid].get_attribute('score'))
        
        #self.assertEqual(len(tax_profile.tree.data), 47)
        self.assertTrue(tax_profile.tree.data)

    
#    @unittest.skip("for faster testing")
    def test_4_build_taxonomy_profile(self):
        
        for sample in self.project.list_samples():
            self.project.import_reads_json(sample,ENDS)
        metrics = 'efpkg'
        scores_function_taxonomy = get_function_taxonomy_scores(
            self.project, sample_id, metrics=metrics
        )
        print(scores_function_taxonomy)
            #~ for end in ENDS:
                
                #~ #print(sample, end, 'read count', str(len(reads)))
                #~ scores = defaultdict(lambda : defaultdict(float))
                #~ for read_id,read in self.project.samples[sample_id].reads[end].items():
                    
                    #~ if read.status == 'function':
                        #~ hits = read.hit_list.hits
                        #~ for hit in hits:
                            #~ protein_taxid = self.project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.subject_id))
                            #~ scores[protein_taxid]['count'] += 1.0
                            #~ scores[protein_taxid]['identity'] += hit.identity
                        #~ if len(hits) == 1:
                            #~ read_functions = read.functions
                            #~ for function in read_functions:
                                #~ scores[self.project.ref_data.lookup_protein_tax(cleanup_protein_id(hits[0].subject_id))]['rpkm'] += read_functions[function]
                        #~ else:
                            #~ read_functions = read.functions
                            #~ protein_taxids = {}
                            #~ for hit in hits:
                                #~ hit_taxid = self.project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.subject_id))
                                #~ hit_functions = hit.functions
                                #~ for hit_function in hit_functions:
                                    #~ protein_taxids[hit_taxid] = hit_function
                            #~ for taxid in protein_taxids:
                                #~ if protein_taxids[taxid] in read_functions:
                                    #~ scores[taxid]['rpkm'] += read_functions[protein_taxids[taxid]]
                #~ print(sample, end, 'tax id count', str(len(scores)))
        scores = defaultdict(lambda : defaultdict(float))
        for tax_id in scores_function_taxonomy:
            for func_id in scores_function_taxonomy[tax_id]:
                for key, val in scores_function_taxonomy[tax_id][func_id][sample_id].items():
                    scores[tax_id][key] += val
        print('Scores',str(scores))
        tax_profile = TaxonomyProfile()
        outfile = os.path.join(
            self.project.options.get_project_dir(sample_id),
            self.project.options.get_output_subdir(sample_id),
            sample_id + '_' + end + '_' + 'taxonomy_profile.xml'
        )
        tax_profile.make_taxonomy_profile(self.project.taxonomy_data, scores)
        print ('Root children:',','.join(tax_profile.tree.data['1'].children))
        if '2157' in tax_profile.tree.data:
            print ('Archaea found') 
        if '10239' in tax_profile.tree.data:
            print ('Viruses found') 
        print(tax_profile.stringify_taxonomy_profile())
        make_taxonomy_chart(
            tax_profile, sample_id, outfile, self.project.config.krona_path, score=metrics
        )
        self.assertTrue(tax_profile.tree.data)
                
        #~ self.project.samples[sample].reads[end] = None


#    @unittest.skip("for faster testing")
    def test_5_build_sample_function_taxonomy_profile(self):
        for sample in self.project.list_samples():
            self.project.import_reads_json(sample,ENDS)
        scores = defaultdict(lambda : defaultdict(dict))
        metric = 'efpkg'
        all_scores = get_function_taxonomy_scores(self.project, sample_id, metric=metric)
        
        functions = set()
        for tax_id in all_scores:
            for func_id in all_scores[tax_id]:
                functions.add(func_id)
                scores[tax_id][func_id] = all_scores[tax_id][func_id][sample_id]
        #~ for sample in self.project.list_samples():
            #~ self.project.import_reads_json(sample,ENDS)
            #~ for end in ENDS:
                #~ scaling_factor = 1.0
                #~ if end == 'pe1':
                    #~ scaling_factor = self.project.options.get_fastq1_readcount(sample)/(
                        #~ self.project.options.get_fastq1_readcount(sample)
                        #~ + self.project.options.get_fastq2_readcount(sample)
                    #~ )
                #~ elif end == 'pe2':
                    #~ scaling_factor = self.project.options.get_fastq2_readcount(sample)/(
                        #~ self.project.options.get_fastq1_readcount(sample)
                        #~ + self.project.options.get_fastq2_readcount(sample)
                    #~ )
                #~ else:
                    #~ raise Exception('Unknown end identifier')

                #~ multiple_hits = 0
                #~ read_count = 0
                
                #~ for read_id in self.project.samples[sample].reads[end]:
                    #~ read = self.project.samples[sample].reads[end][read_id]
                    #~ if read.status == 'function':
                        
                        #~ read_functions = read.functions
                        #~ hits = read.hit_list.hits
                        #~ if len(hits) >1:
                            #~ multiple_hits += 1
                        #~ read_count += 1
                        #~ # Count hits and identity
                        #~ for hit in hits:
                            #~ hit_taxid = self.project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.subject_id))
                            #~ if sample in scores[hit_taxid]:
                                #~ scores[hit_taxid][sample]['count'] += 1.0/len(hits)
                                #~ scores[hit_taxid][sample]['hit_count'] += 1.0
                                #~ scores[hit_taxid][sample]['identity'] += hit.identity
                            #~ else:
                                #~ scores[hit_taxid][sample]['count'] = 1.0/len(hits)
                                #~ scores[hit_taxid][sample]['hit_count'] = 1.0 # we need hit count to calculate average identity
                                #~ scores[hit_taxid][sample]['identity'] = hit.identity
                                #~ # Initialize 'rpkm' here
                                #~ scores[hit_taxid][sample]['rpkm'] = 0.0

                            
                        #~ # Count RPKM
                        #~ # If we have only one hit, all RPKM scores of the read would be assigned to the tax id of the hit
                        #~ # If we have more than one hit, RPKM scores would be equally divided between tax ids of the hits, with regard to functional assignments of the hits
                        #~ function_taxids = defaultdict(list)
                        #~ # First, collect all taxonomy IDs for each function assigned to the hit
                        #~ for function in read_functions:
                            #~ for hit in hits:
                                #~ hit_taxid = self.project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.subject_id))
                                #~ for hit_function in hit.functions:
                                    #~ if hit_function == function:
                                        #~ function_taxids[function].append(hit_taxid)
                        #~ # Second, for each tax ID add its share of the RPKM score assigned to hit's function
                        #~ for function in function_taxids:
                            #~ tax_count = len(function_taxids[function])
                            #~ for hit_taxid in function_taxids[function]:
                                #~ scores[hit_taxid][sample]['rpkm'] += read_functions[function] * scaling_factor/tax_count
                #~ print(sample, end, 'read count', str(read_count))
                #~ print ('Reads with multiple hits: ', multiple_hits)
                
                #~ self.project.samples[sample].reads[end] = None
                
        # Now, we have all taxonomy ids, from each samples, in a single nested dictionary. 
        print('tax id count', str(len(scores)))

        
        tax_profile = TaxonomyProfile()
        outfile = os.path.join(
            self.project.options.work_dir,
            self.project.options.project_name + '_taxonomy_profile.xml'
            )
        outfile = outfile.replace(' ', '_')
        outfile = outfile.replace("'", "")

        tax_profile.make_function_taxonomy_profile(
            self.project.taxonomy_data,
            scores
            )
        print(tax_profile.stringify_taxonomy_profile())
        make_function_taxonomy_chart(
#        make_taxonomy_series_chart(
            tax_profile, sorted(list(functions)),  # sorted(self.project.samples.keys()),
            outfile, self.project.config.krona_path
            )
        self.assertTrue(tax_profile.tree.data)


#    @unittest.skip("for faster testing")
    def test_6_build_project_taxonomy_profile_1function(self):
        
        target_function = 'GH9'
        
        
        scores = defaultdict(lambda : defaultdict(dict))
        
        for sample in self.project.list_samples():
            self.project.import_reads_json(sample,ENDS)
            for end in ENDS:
                
                
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
                    if read.status == 'function':
                        
                        read_functions = read.functions
                        
                        if target_function not in read_functions:
                            continue
                            
                        hits = read.hit_list.hits
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
        outfile = os.path.join(self.project.options.work_dir, self.project.options.project_name + '_' + target_function + '_taxonomy_profile.xml')
        outfile = outfile.replace(' ', '_')
        outfile = outfile.replace("'", "")

        tax_profile.make_function_taxonomy_profile(self.project.taxonomy_data, scores)

        make_taxonomy_series_chart(tax_profile, sorted(self.project.list_samples()), outfile, self.project.config.krona_path)
        self.assertTrue(tax_profile.tree.data)


#    @unittest.skip("for faster testing")
    def test_7_build_loose_taxonomy_profile(self):
        metrics = 'rpkm'
        for sample in self.project.list_samples():
            if sample != sample_id:
                continue
            self.project.import_reads_json(sample, ENDS)
            
            norm_factor = 1000000/self.project.options.get_fastq1_readcount(sample_id)

            for end in ENDS:
                
                #print(sample, end, 'read count', str(len(reads)))
                scores = autovivify(2,float)
                for read_id,read in self.project.samples[sample_id].reads[end].items():
                    
                    if read.status != 'function':
                        continue
                    
                    if read.taxonomy is None:
                        print ('No taxonomy ID assigned to ' + read_id)
                       # raise ValueError('No taxonomy ID assigned to ' + read_id)
                    else:
                        scores[read.taxonomy]['count'] += 1
                        scores[read.taxonomy][metrics] += norm_factor * sum(list(read.functions.values()))
                        scores[read.taxonomy]['hit_count'] += 1
                        scores[read.taxonomy]['identity'] += max([hit.identity for hit in read.hit_list.hits])

                print(sample, end, 'tax id count', str(len(scores)))
                tax_profile = TaxonomyProfile()
                outfile = os.path.join(self.project.options.get_project_dir(sample), self.project.options.get_output_subdir(sample), sample + '_' + end + '_' + metrics + '_taxonomy_profile.xml')
                tax_profile.make_taxonomy_profile(self.project.taxonomy_data, scores)
                #print ('Root children:',','.join(tax_profile.tree.data['1'].children))
                #if '2157' in tax_profile.tree.data:
                #    print ('Archaea found') 
                #if '10239' in tax_profile.tree.data:
                #    print ('Viruses found') 
                #print(tax_profile.print_taxonomy_profile())
                make_taxonomy_chart(tax_profile, sample, outfile, self.project.config.krona_path, score=metrics)
                self.assertTrue(tax_profile.tree.data)
                
            self.project.samples[sample].reads[end] = None


#    @unittest.skip("for faster testing")
    def test_8_build_loose_project_taxonomy_profile(self):
        metrics = 'rpkm'
        scores = autovivify(3,float)
        for sample in self.project.list_samples():
            self.project.import_reads_json(sample, ENDS)
            
            norm_factor = 1000000/self.project.options.get_fastq1_readcount(sample)

            for end in ENDS:
                
                #print(sample, end, 'read count', str(len(reads)))
                for read_id,read in self.project.samples[sample].reads[end].items():
                    
                    if read.status != 'function':
                        continue
                    
                    if read.taxonomy is None:
                        print ('No taxonomy ID assigned to ' + read_id)
                       # raise ValueError('No taxonomy ID assigned to ' + read_id)
                    else:
                        scores[read.taxonomy][sample]['count'] += 1
                        scores[read.taxonomy][sample][metrics] += norm_factor * sum(list(read.functions.values()))
                        scores[read.taxonomy][sample]['hit_count'] += 1
                        scores[read.taxonomy][sample]['identity'] += max([hit.identity for hit in read.hit_list.hits])

                print(sample, end, 'tax id count', str(len(scores)))
                
            self.project.samples[sample].reads[end] = None


        tax_profile = TaxonomyProfile()
        outfile = sanitize_file_name(os.path.join(self.project.options.work_dir, self.project.options.project_name + '_' + metrics + '_taxonomy_profile.xml'))

        tax_profile.make_function_taxonomy_profile(self.project.taxonomy_data, scores)
        make_taxonomy_series_chart(tax_profile, sorted(self.project.samples.keys()), outfile, self.project.config.krona_path, score=metrics)
        self.assertTrue(tax_profile.tree.data)

    def tearDown(self):
        self.parser = None


if __name__ == '__main__':
    unittest.main()
