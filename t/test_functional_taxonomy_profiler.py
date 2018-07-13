#!/usr/bin/python
import os, csv, operator
import unittest
from collections import Counter,defaultdict

import xlsxwriter
import pandas as pd

from context import Fama
from Fama.DiamondParser.DiamondParser import DiamondParser
from Fama.OutputUtil.JSONUtil import import_annotated_reads
from Fama.DiamondParser.hit_utils import cleanup_protein_id
from Fama.ReferenceLibrary.TaxonomyData import TaxonomyData
from Fama.OutputUtil.Report import generate_report
from Fama.TaxonomyProfile import TaxonomyProfile
from Fama.OutputUtil.KronaXMLWriter import generate_functional_taxonomy_chart
from Fama.Project import Project

data_dir = 'data'
config_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'config.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project.ini')
#sample = 'test_sample'
end = 'pe1'

#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_Hans_nitrogen_t.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_EB271_nitrogen_t.ini')
sample = 'HL1G'

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
    
    
    def test_4_build_functional_taxonomy_profile(self):
        tax_data = TaxonomyData(self.parser.config)
        tax_data.load_taxdata(self.parser.config)
        
        project = Project(config_file=config_path, project_file=project_path)
        project.load_functional_profile()
        outfile = project.options.get_name() + '_functions_taxonomy.xlsx'
        #outfile = 'test_functions_taxonomy.xlsx'
        outfile = outfile.replace(' ', '_')
        outfile = outfile.replace("'", "")

        writer = pd.ExcelWriter(outfile, engine='xlsxwriter')

        for sample in sorted(project.samples.keys()):
            scores = defaultdict(lambda : defaultdict(dict))
            function_list = set()


            for end in sorted(project.samples[sample].keys()):

                scaling_factor = 1.0
                if end == 'pe1':
                    scaling_factor = project.options.get_fastq1_readcount(sample)/(project.options.get_fastq1_readcount(sample) + project.options.get_fastq2_readcount(sample))
                elif end == 'pe2':
                    scaling_factor = project.options.get_fastq2_readcount(sample)/(project.options.get_fastq1_readcount(sample) + project.options.get_fastq2_readcount(sample))
                else:
                    raise Exception('Unknown end identifier')

                reads = project.samples[sample][end]
#        self.parser.reads = import_annotated_reads(os.path.join(self.parser.project.get_project_dir(self.parser.sample), self.parser.sample + '_' + self.parser.end + '_' + self.parser.project.get_reads_json_name()))
#        reads = self.parser.reads


                multiple_hits = 0
                read_count = 0
                
                for read in reads:
                    if reads[read].get_status() == 'function,besthit' or reads[read].get_status() == 'function':
                        
                        read_functions = reads[read].get_functions()

                        hits = reads[read].get_hit_list().get_hits()
                        if len(hits) >1:
                            multiple_hits += 1
                        read_count += 1

                        # Count RPKM
                        # If we have only one hit, all RPKM scores of the read would be assigned to the tax id of the hit
                        # If we have more than one hit, RPKM scores would be equally divided between tax ids of the hits, with regard to functional assignments of the hits
                        function_taxids = defaultdict(lambda : defaultdict(float))
                        # First, collect all taxonomy IDs for each function assigned to the hit
                        tax_functions_count = 0.0
                        
                        for read_function in read_functions:
                            for hit in hits:
                                hit_taxid = self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))

                                hit_functions = hit.get_functions()
                                for hit_function in hit_functions:
                                    if hit_function == read_function:
                                        tax_functions_count += 1.0
                                        # Add function to non-redundant list of functions
                                        function_list.add(read_function)
                                        
                                        function_taxids[read_function][hit_taxid] += 1.0
                                        # Count identity and hits (per function) here
                                        if read_function in scores[hit_taxid]:
                                            scores[hit_taxid][read_function]['hit_count'] += 1.0
                                            scores[hit_taxid][read_function]['identity'] += hit.get_identity()
                                        else:
                                            scores[hit_taxid][read_function]['hit_count'] = 1.0 # hit count required for average identity calculation
                                            scores[hit_taxid][read_function]['identity'] = hit.get_identity()
                                            # Initialize 'count' and 'rpkm' here
                                            scores[hit_taxid][read_function]['count'] = 0.0
                                            scores[hit_taxid][read_function]['rpkm'] = 0.0

                                        
                        # Second, for each tax ID add its share of read count and RPKM score assigned to hit's function
                        for function in function_taxids:
                            tax_count = len(function_taxids[function])
                            for hit_taxid in function_taxids[function]:
                                if hit_taxid in scores and function in scores[hit_taxid]:
                                    scores[hit_taxid][function]['rpkm'] += read_functions[function] * scaling_factor / tax_count
                                    scores[hit_taxid][function]['count'] += function_taxids[function][hit_taxid] / tax_functions_count


#                print(sample, end, 'read count', str(read_count))
#                print ('Reads with multiple hits: ', multiple_hits)


        #for read in reads:
            #if reads[read].get_status() == 'function,besthit' or reads[read].get_status() == 'function':
                #read_functions = reads[read].get_functions()
                #hits = reads[read].get_hit_list().get_hits()
                #for function in read_functions:
                    #for hit in hits:
                        #protein_taxid = self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                        #for hit_function in hit.get_functions():
                            #if hit_function == function:
                                #function_list.add(function)
                                #if function in scores[protein_taxid]:
                                    #scores[protein_taxid][function]['count'] += 1.0
                                    #scores[protein_taxid][function]['identity'] += hit.get_identity()
                                    #scores[protein_taxid][function]['rpkm'] += read_functions[function]
                                #else:
                                    #scores[protein_taxid][function]['count'] = 1.0
                                    #scores[protein_taxid][function]['identity'] = hit.get_identity()
                                    #scores[protein_taxid][function]['rpkm'] = read_functions[function]

            tax_profile = TaxonomyProfile()
            
            #outfile = os.path.join(self.parser.project.get_project_dir(sample), self.parser.project.get_output_subdir(sample), sample + '_' + end + '_' + 'functional_taxonomy_profile.xml')
            outfile = os.path.join(sample + '_' + 'functional_taxonomy_profile.xml')
            tax_profile.build_functional_taxonomy_profile(tax_data, scores)
            #print(tax_profile.print_functional_taxonomy_profile())
            generate_functional_taxonomy_chart(tax_profile, sorted(function_list), outfile)
            #print(tax_profile.print_functional_taxonomy_table())
            df = tax_profile.convert_function_taxonomic_profile_into_df()
        
            #print(df)
            df.to_excel(writer, sheet_name=sample)
            workbook  = writer.book
            worksheet = writer.sheets[sample]
            superkingdom_format = workbook.add_format({'bg_color': '#FF6666'})
            phylum_format = workbook.add_format({'bg_color': '#FF9900'})
            class_format = workbook.add_format({'bg_color': '#FFCC99'})
            order_format = workbook.add_format({'bg_color': '#FFFFCC'})
            family_format = workbook.add_format({'bg_color': '#99FFCC'})
            genus_format = workbook.add_format({'bg_color': '#99FFFF'})
            worksheet.conditional_format('C4:C1048560', {'type':     'text',
                                   'criteria': 'containing',
                                   'value':    'superkingdom',
                                   'format':   superkingdom_format})
            worksheet.conditional_format('C4:C1048560', {'type':     'text',
                                   'criteria': 'containing',
                                   'value':    'phylum',
                                   'format':   phylum_format})
            worksheet.conditional_format('C4:C1048560', {'type':     'text',
                                   'criteria': 'containing',
                                   'value':    'class',
                                   'format':   class_format})
            worksheet.conditional_format('C4:C1048560', {'type':     'text',
                                   'criteria': 'containing',
                                   'value':    'order',
                                   'format':   order_format})
            worksheet.conditional_format('C4:C1048560', {'type':     'text',
                                   'criteria': 'containing',
                                   'value':    'family',
                                   'format':   family_format})
            worksheet.conditional_format('C4:C1048560', {'type':     'text',
                                   'criteria': 'containing',
                                   'value':    'genus',
                                   'format':   genus_format})

            
            worksheet.set_column(1, 1, 30)
            worksheet.set_column(2, 2, 15)





        writer.save()
        
        self.assertTrue(tax_profile.tree.data)
            
        #for sample in sorted(project.samples.keys()):
            #for end in sorted(project.samples[sample].keys()):

                #reads = project.samples[sample][end]
                ##print(sample, end, 'read count', str(len(reads)))
                #scores = defaultdict(lambda : defaultdict(float))
                #for read in reads:
                    #if reads[read].get_status() == 'function,besthit' or reads[read].get_status() == 'function':
                        #hits = reads[read].get_hit_list().get_hits()
                        #for hit in hits:
                            #protein_taxid = project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                            #scores[protein_taxid]['count'] += 1.0
                            #scores[protein_taxid]['identity'] += hit.get_identity()
                        #if len(hits) == 1:
                            #read_functions = reads[read].get_functions()
                            #for function in read_functions:
                                #scores[project.ref_data.lookup_protein_tax(cleanup_protein_id(hits[0].get_subject_id()))]['rpkm'] += read_functions[function]
                        #else:
                            #read_functions = reads[read].get_functions()
                            #protein_taxids = {}
                            #for hit in hits:
                                #hit_taxid = project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                                #hit_functions = hit.get_functions()
                                #for hit_function in hit_functions:
                                    #protein_taxids[hit_taxid] = hit_function
                            #for taxid in protein_taxids:
                                #if protein_taxids[taxid] in read_functions:
                                    #scores[taxid]['rpkm'] += read_functions[protein_taxids[taxid]]
                #print(sample, end, 'tax id count', str(len(scores)))
                #tax_profile = TaxonomyProfile()
                #outfile = os.path.join(project.options.get_project_dir(sample), project.options.get_output_subdir(sample), sample + '_' + end + '_' + 'taxonomy_profile.xml')
                #tax_profile.build_taxonomy_profile(tax_data, scores)
                ##print ('Root children:',','.join(tax_profile.tree.data['1'].children))
                ##if '2157' in tax_profile.tree.data:
                ##    print ('Archaea found') 
                ##if '10239' in tax_profile.tree.data:
                ##    print ('Viruses found') 
                ##print(tax_profile.print_taxonomy_profile())
                #generate_taxonomy_chart(tax_profile, sample, outfile)
                #self.assertTrue(tax_profile.tree.data)


    def tearDown(self):
        self.parser = None

if __name__=='__main__':
    unittest.main()
