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

#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_EB271_nitrogen8_t.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_fw106_fw301_sulfate.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW306_cazy1_t.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_B1F_cazy_t.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'protein_project_EB271_assembly.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW306_nitrogen8test_t.ini')
sample = 'sample1'

ENDS = ('pe1', 'pe2')

def autovivify(levels=1, final=dict):
    return (defaultdict(final) if levels < 2 else
            defaultdict(lambda: autovivify(levels - 1, final)))

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
    
    @unittest.skip("for faster testing")
    def test_4_build_functional_taxonomy_profile(self):
        tax_data = TaxonomyData(self.parser.config)
        tax_data.load_taxdata(self.parser.config)
        
        project = Project(config_file=config_path, project_file=project_path)
        #project.load_functional_profile()
        outfile = project.options.get_name() + '_functions_taxonomy.xlsx'
        #outfile = 'test_functions_taxonomy.xlsx'
        outfile = outfile.replace(' ', '_')
        outfile = outfile.replace("'", "")

        writer = pd.ExcelWriter(outfile, engine='xlsxwriter')

        for sample in sorted(project.list_samples()):
            scores = defaultdict(lambda : defaultdict(dict))
            function_list = set()

            for end in ENDS:
                project.load_annotated_reads(sample, end)

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
                
                for read_id in project.samples[sample][end]:
                    read = project.samples[sample][end][read_id]
                    if read.get_status() == 'function,besthit' or read.get_status() == 'function':
                        
                        read_functions = read.get_functions()

                        hits = read.get_hit_list().get_hits()
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
                project.samples[sample][end] = None
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

    @unittest.skip("temporary excluded")
    def test_5_build_fpkm_taxonomy_profile(self):
        tax_data = TaxonomyData(self.parser.config)
        tax_data.load_taxdata(self.parser.config)
        
        project = Project(config_file=config_path, project_file=project_path)
        project.load_functional_profile()
        outfile = project.options.get_name() + '_fpkm_functions_taxonomy.xlsx'
        #outfile = 'test_functions_taxonomy.xlsx'
        outfile = outfile.replace(' ', '_')
        outfile = outfile.replace("'", "")
        writer = pd.ExcelWriter(outfile, engine='xlsxwriter')
        ENDS = ('pe1', 'pe2')
        fpkm_scores = autovivify(2, float)
        all_functions_list = set()
        for sample in sorted(project.samples.keys()):
            scores = autovivify(3, float)
            reads_processed = set()
            function_list = set()

            reads_pe1 = project.samples[sample]['pe1']
            reads_pe2 = {}
            if 'pe2' in project.samples[sample]:
                reads_pe2 = project.samples[sample]['pe2']
                
            
            for read in reads_pe1:
                if reads_pe1[read].get_status() == 'function,besthit' or reads_pe1[read].get_status() == 'function':
                    reads_processed.add(read)
                                            
                    if read in reads_pe2: 
                        if reads_pe2[read].get_status() == 'function,besthit' or reads_pe2[read].get_status() == 'function':
                            # Both ends are mapped
                            fpkm = defaultdict(float) # Stores FPKM scores for the current read
                            read_functions = set() # List of functions assigned to the current read
                            if self.have_similar_functions(reads_pe1[read], reads_pe2[read]): # Ends have overlapping functions, count them together

                                read1_functions = reads_pe1[read].get_functions()
                                read2_functions = reads_pe2[read].get_functions()

                                read_functions.update(read1_functions.keys())
                                read_functions.update(read2_functions.keys())
                                # We filled list of functions. Let's calculate FPKM score
                                
                                for read_function in read_functions:
                                    if read_function in read1_functions and read_function in read2_functions: # Calculate average score
                                        fpkm[read_function] += (read1_functions[read_function] + read2_functions[read_function]) / 2
                                    elif read_function in read1_functions:
                                        fpkm[read_function] += read1_functions[read_function]
                                    elif read_function in read2_functions:
                                        fpkm[read_function] += read2_functions[read_function]
                            else:
                                read_functions = reads_pe1[read].get_functions()
                                for read_function in read_functions:
                                    fpkm[read_function] += read_functions[read_function]
                                read_functions = reads_pe2[read].get_functions()
                                for read_function in read_functions:
                                    fpkm[read_function] += read_functions[read_function]
                            
                            # Add FPKM score of the read to FPKM score of the sample
                            for function in fpkm:
                                all_functions_list.add(function)
                                fpkm_scores[sample][function] += fpkm[function]
                            
                            # Build FPKM-based functional taxonomic profile
                            tax_functions_count = 0.0 # Total number of functions assigned to all hits for the read
                            function_taxids = autovivify(2, float) # for each function, store list of tax_ids and count of hits for each tax_id
                            max_bitscore = -1.0
                            hits1 = reads_pe1[read].get_hit_list().get_hits()
                            hits2 = reads_pe2[read].get_hit_list().get_hits()
                            for read_function in read_functions:
                                # Find max bitscore for the function
                                if read_function in read1_functions and reads_pe1[read].get_status() == 'function,besthit':
                                    for hit in hits1:
                                        for hit_function in hit.get_functions():
                                            if hit_function == read_function and hit.get_bitscore() > max_bitscore:
                                                function_list.add(read_function)
                                                max_bitscore = hit.get_bitscore()
                                if read_function in read2_functions and reads_pe2[read].get_status() == 'function,besthit':
                                    for hit in hits2:
                                        for hit_function in hit.get_functions():
                                            if hit_function == read_function and hit.get_bitscore() > max_bitscore:
                                                function_list.add(read_function)
                                                max_bitscore = hit.get_bitscore()

                                # Ignore all hits for this function with bitscore < cutoff. 

                                if max_bitscore > 0.0:
                                    min_bitscore = max_bitscore * (1 - project.config.get_biscore_range_cutoff(project.options.get_collection()))
                                    for hit in (hit for hit in hits1 if reads_pe1[read].get_status() == 'function,besthit'):
                                        hit_taxid = self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                                        for hit_function in (hit_function for hit_function in hit.get_functions() if hit.get_bitscore() > min_bitscore and hit_function == read_function):
                                            # Count hits
                                            tax_functions_count += 1.0
                                            # Count hits per function per taxid
                                            function_taxids[read_function][hit_taxid] += 1.0
                                            # Count identity and hits (per function) 
                                            scores[hit_taxid][read_function]['hit_count'] += 1.0 
                                            scores[hit_taxid][read_function]['identity'] += hit.get_identity()

                                    for hit in (hit for hit in hits2 if reads_pe2[read].get_status() == 'function,besthit'):
                                        hit_taxid = self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                                        for hit_function in (hit_function for hit_function in hit.get_functions() if hit.get_bitscore() > min_bitscore and hit_function == read_function):
                                            # Count hits
                                            tax_functions_count += 1.0
                                            # Count hits per function per taxid
                                            function_taxids[read_function][hit_taxid] += 1.0
                                            # Count identity and hits (per function) here
                                            scores[hit_taxid][read_function]['hit_count'] += 1.0
                                            scores[hit_taxid][read_function]['identity'] += hit.get_identity()
                                
                            # Second, for each tax ID add its share of read count and FPKM score assigned to hit's function
                            for function in function_taxids:
                                if function in fpkm:
                                    tax_count = sum(function_taxids[function].values()) # How many times this function was counted in all tax_ids
                                    for hit_taxid in function_taxids[function]:
                                        scores[hit_taxid][function]['fpkm'] += fpkm[function] * function_taxids[function][hit_taxid] / tax_count
                                        scores[hit_taxid][function]['count'] += function_taxids[function][hit_taxid] / tax_functions_count
                                
                    else: #Only end1 is mapped
                        read_functions = reads_pe1[read].get_functions()
                        hits = reads_pe1[read].get_hit_list().get_hits()
                        # Count FPKM
                        for function in read_functions:
                            all_functions_list.add(function)
                            fpkm_scores[sample][function] += read_functions[function]

                        # First, collect all taxonomy IDs for each function assigned to the hit
                        tax_functions_count = 0.0
                        function_taxids = autovivify(2, float) # for each function, store list of tax_ids and count of hits for each tax_id
                        
                        if reads_pe1[read].get_status() == 'function,besthit':
                            for read_function in read_functions:
                                for hit in hits:
                                    hit_taxid = self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                                    for hit_function in hit.get_functions():
                                        if hit_function == read_function:
                                            # Count hits
                                            tax_functions_count += 1.0
                                            # Add function to non-redundant list of functions
                                            function_list.add(read_function)
                                            # Count hits per function per taxid
                                            function_taxids[read_function][hit_taxid] += 1.0
                                            # Count identity and hits (per function) here
                                            scores[hit_taxid][read_function]['hit_count'] += 1.0
                                            scores[hit_taxid][read_function]['identity'] += hit.get_identity()
                                            
                            # Second, for each tax ID add its share of read count and RPKM score assigned to hit's function
                            for function in function_taxids:
                                    tax_count = sum(function_taxids[function].values()) # How many times this function was counted in all tax_ids
                                    for hit_taxid in function_taxids[function]:
                                        scores[hit_taxid][function]['fpkm'] += read_functions[function] * function_taxids[function][hit_taxid] / tax_count
                                        scores[hit_taxid][function]['count'] += function_taxids[function][hit_taxid] / tax_functions_count

            for read in reads_pe2:
                if read not in reads_processed: #Only end2 was mapped
                    if reads_pe2[read].get_status() == 'function,besthit' or reads_pe2[read].get_status() == 'function':
                        read_functions = reads_pe2[read].get_functions()
                        hits = reads_pe2[read].get_hit_list().get_hits()
                        # Count FPKM
                        for function in read_functions:
                            all_functions_list.add(function)
                            fpkm_scores[sample][function] += read_functions[function]

                        tax_functions_count = 0.0
                        function_taxids = autovivify(2, float) # for each function, store list of tax_ids and count of hits for each tax_id

                        if reads_pe2[read].get_status() == 'function,besthit':
                            for read_function in read_functions:
                                for hit in hits:
                                    hit_taxid = self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                                    for hit_function in hit.get_functions():
                                        if hit_function == read_function:
                                            # Count hits
                                            tax_functions_count += 1.0
                                            # Add function to non-redundant list of functions
                                            function_list.add(read_function)
                                            # Count hits per function per taxid
                                            function_taxids[read_function][hit_taxid] += 1.0
                                            # Count identity and hits (per function) here
                                            scores[hit_taxid][read_function]['hit_count'] += 1.0
                                            scores[hit_taxid][read_function]['identity'] += hit.get_identity()

                            # Second, for each tax ID add its share of read count and RPKM score assigned to hit's function
                            for function in function_taxids:
                                tax_count = sum(function_taxids[function].values()) # How many times this function was counted in all tax_ids
                                for hit_taxid in function_taxids[function]:
                                    scores[hit_taxid][function]['fpkm'] += read_functions[function] * function_taxids[function][hit_taxid] / tax_count
                                    scores[hit_taxid][function]['count'] += function_taxids[function][hit_taxid] / tax_functions_count

            tax_profile = TaxonomyProfile()
            
            #outfile = os.path.join(self.parser.project.get_project_dir(sample), self.parser.project.get_output_subdir(sample), sample + '_' + end + '_' + 'functional_taxonomy_profile.xml')
            outfile = os.path.join(sample + '_fpkm_functional_taxonomy_profile.xml')
            tax_profile.build_functional_taxonomy_profile(tax_data, scores)
            #print(tax_profile.print_functional_taxonomy_profile())
            generate_functional_taxonomy_chart(tax_profile, sorted(function_list), outfile, 'fpkm')
            #print(tax_profile.print_functional_taxonomy_table())
            df = tax_profile.convert_function_taxonomic_profile_into_df('fpkm')
        
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
        
        # Add FPKM worksheet
        bold = workbook.add_format({'bold': True})
        scores_worksheet = workbook.add_worksheet('Functions FPKM')
        row = 0
        col = 0
        scores_worksheet.write(row, col, 'Function', bold)
        for sample in sorted(project.samples.keys()):
            col += 1
            scores_worksheet.write(row, col, sample, bold)
        col += 1
        scores_worksheet.write(row, col, 'Definition', bold)
        for function in sorted(all_functions_list):
            row += 1
            col = 0
            scores_worksheet.write(row, col, function, bold)
            for sample in sorted(project.samples.keys()):
                col += 1
                if function in fpkm_scores[sample]:
                    scores_worksheet.write(row, col, '{0:.3f}'.format(fpkm_scores[sample][function]))
                else:
                    scores_worksheet.write(row, col, 0.0)
            col += 1
            scores_worksheet.write(row, col, project.ref_data.lookup_function_name(function))
        # adjust column width
        scores_worksheet.set_column(0, 0, 10)
        scores_worksheet.set_column(col, col, 50)

        writer.save()
        
        self.assertTrue(tax_profile.tree.data)

    @unittest.skip("for faster testing")
    def test_6_lazy_build_fpkm_taxonomy_profile(self):
        tax_data = TaxonomyData(self.parser.config)
        tax_data.load_taxdata(self.parser.config)
        
        project = Project(config_file=config_path, project_file=project_path)
        #project.load_functional_profile()
        outfile = project.options.get_name() + '_fpkm_functions_taxonomy.xlsx'
        #outfile = 'test_functions_taxonomy.xlsx'
        outfile = outfile.replace(' ', '_')
        outfile = outfile.replace("'", "")
        writer = pd.ExcelWriter(outfile, engine='xlsxwriter')
        ENDS = ('pe1', 'pe2')
        fpkm_scores = autovivify(2, float)
        all_functions_list = set()
        for sample in project.list_samples():
            for end in ENDS:
                project.load_annotated_reads(sample, end)

            scores = autovivify(3, float)
            reads_processed = set()
            function_list = set()

            reads_pe1 = project.samples[sample]['pe1']
            reads_pe2 = {}
            if 'pe2' in project.samples[sample]:
                reads_pe2 = project.samples[sample]['pe2']
                
            
            for read in reads_pe1:
                if reads_pe1[read].get_status() == 'function,besthit' or reads_pe1[read].get_status() == 'function':
                    reads_processed.add(read)
                                            
                    if read in reads_pe2: 
                        if reads_pe2[read].get_status() == 'function,besthit' or reads_pe2[read].get_status() == 'function':
                            # Both ends are mapped
                            fpkm = defaultdict(float) # Stores FPKM scores for the current read
                            read_functions = set() # List of functions assigned to the current read
                            if self.have_similar_functions(reads_pe1[read], reads_pe2[read]): # Ends have overlapping functions, count them together

                                read1_functions = reads_pe1[read].get_functions()
                                read2_functions = reads_pe2[read].get_functions()

                                read_functions.update(read1_functions.keys())
                                read_functions.update(read2_functions.keys())
                                # We filled list of functions. Let's calculate FPKM score
                                
                                for read_function in read_functions:
                                    if read_function in read1_functions and read_function in read2_functions: # Calculate average score
                                        fpkm[read_function] += (read1_functions[read_function] + read2_functions[read_function]) / 2
                                    elif read_function in read1_functions:
                                        fpkm[read_function] += read1_functions[read_function]
                                    elif read_function in read2_functions:
                                        fpkm[read_function] += read2_functions[read_function]
                            else:
                                read_functions = reads_pe1[read].get_functions()
                                for read_function in read_functions:
                                    fpkm[read_function] += read_functions[read_function]
                                read_functions = reads_pe2[read].get_functions()
                                for read_function in read_functions:
                                    fpkm[read_function] += read_functions[read_function]
                            
                            # Add FPKM score of the read to FPKM score of the sample
                            for function in fpkm:
                                all_functions_list.add(function)
                                fpkm_scores[sample][function] += fpkm[function]
                            
                            # Build FPKM-based functional taxonomic profile
                            tax_functions_count = 0.0 # Total number of functions assigned to all hits for the read
                            function_taxids = autovivify(2, float) # for each function, store list of tax_ids and count of hits for each tax_id
                            max_bitscore = -1.0
                            hits1 = reads_pe1[read].get_hit_list().get_hits()
                            hits2 = reads_pe2[read].get_hit_list().get_hits()
                            for read_function in read_functions:
                                # Find max bitscore for the function
                                if read_function in read1_functions and reads_pe1[read].get_status() == 'function,besthit':
                                    for hit in hits1:
                                        for hit_function in hit.get_functions():
                                            if hit_function == read_function and hit.get_bitscore() > max_bitscore:
                                                function_list.add(read_function)
                                                max_bitscore = hit.get_bitscore()
                                if read_function in read2_functions and reads_pe2[read].get_status() == 'function,besthit':
                                    for hit in hits2:
                                        for hit_function in hit.get_functions():
                                            if hit_function == read_function and hit.get_bitscore() > max_bitscore:
                                                function_list.add(read_function)
                                                max_bitscore = hit.get_bitscore()

                                # Ignore all hits for this function with bitscore < cutoff. 

                                if max_bitscore > 0.0:
                                    min_bitscore = max_bitscore * (1 - project.config.get_biscore_range_cutoff(project.options.get_collection()))
                                    for hit in (hit for hit in hits1 if reads_pe1[read].get_status() == 'function,besthit'):
                                        hit_taxid = self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                                        for hit_function in (hit_function for hit_function in hit.get_functions() if hit.get_bitscore() > min_bitscore and hit_function == read_function):
                                            # Count hits
                                            tax_functions_count += 1.0
                                            # Count hits per function per taxid
                                            function_taxids[read_function][hit_taxid] += 1.0
                                            # Count identity and hits (per function) 
                                            scores[hit_taxid][read_function]['hit_count'] += 1.0 
                                            scores[hit_taxid][read_function]['identity'] += hit.get_identity()

                                    for hit in (hit for hit in hits2 if reads_pe2[read].get_status() == 'function,besthit'):
                                        hit_taxid = self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                                        for hit_function in (hit_function for hit_function in hit.get_functions() if hit.get_bitscore() > min_bitscore and hit_function == read_function):
                                            # Count hits
                                            tax_functions_count += 1.0
                                            # Count hits per function per taxid
                                            function_taxids[read_function][hit_taxid] += 1.0
                                            # Count identity and hits (per function) here
                                            scores[hit_taxid][read_function]['hit_count'] += 1.0
                                            scores[hit_taxid][read_function]['identity'] += hit.get_identity()
                                
                            # Second, for each tax ID add its share of read count and FPKM score assigned to hit's function
                            for function in function_taxids:
                                tax_count = sum(function_taxids[function].values()) # How many times this function was counted in all tax_ids
                                for hit_taxid in function_taxids[function]:
                                    scores[hit_taxid][function]['fpkm'] += fpkm[read_function] * function_taxids[function][hit_taxid] / tax_count
                                    scores[hit_taxid][function]['count'] += function_taxids[function][hit_taxid] / tax_functions_count
                                
                    else: #Only end1 is mapped
                        read_functions = reads_pe1[read].get_functions()
                        hits = reads_pe1[read].get_hit_list().get_hits()
                        # Count FPKM
                        for function in read_functions:
                            all_functions_list.add(function)
                            fpkm_scores[sample][function] += read_functions[function]

                        # First, collect all taxonomy IDs for each function assigned to the hit
                        tax_functions_count = 0.0
                        function_taxids = autovivify(2, float) # for each function, store list of tax_ids and count of hits for each tax_id
                        
                        if reads_pe1[read].get_status() == 'function,besthit':
                            for read_function in read_functions:
                                for hit in hits:
                                    hit_taxid = self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                                    for hit_function in hit.get_functions():
                                        if hit_function == read_function:
                                            # Count hits
                                            tax_functions_count += 1.0
                                            # Add function to non-redundant list of functions
                                            function_list.add(read_function)
                                            # Count hits per function per taxid
                                            function_taxids[read_function][hit_taxid] += 1.0
                                            # Count identity and hits (per function) here
                                            scores[hit_taxid][read_function]['hit_count'] += 1.0
                                            scores[hit_taxid][read_function]['identity'] += hit.get_identity()
                                            
                            # Second, for each tax ID add its share of read count and RPKM score assigned to hit's function
                            for function in function_taxids:
                                tax_count = sum(function_taxids[function].values()) # How many times this function was counted in all tax_ids
                                for hit_taxid in function_taxids[function]:
                                    scores[hit_taxid][function]['fpkm'] += read_functions[function] * function_taxids[function][hit_taxid] / tax_count
                                    scores[hit_taxid][function]['count'] += function_taxids[function][hit_taxid] / tax_functions_count

            for read in reads_pe2:
                if read not in reads_processed: #Only end2 was mapped
                    if reads_pe2[read].get_status() == 'function,besthit' or reads_pe2[read].get_status() == 'function':
                        read_functions = reads_pe2[read].get_functions()
                        hits = reads_pe2[read].get_hit_list().get_hits()
                        # Count FPKM
                        for function in read_functions:
                            all_functions_list.add(function)
                            fpkm_scores[sample][function] += read_functions[function]

                        tax_functions_count = 0.0
                        function_taxids = autovivify(2, float) # for each function, store list of tax_ids and count of hits for each tax_id

                        if reads_pe2[read].get_status() == 'function,besthit':
                            for read_function in read_functions:
                                for hit in hits:
                                    hit_taxid = self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                                    for hit_function in hit.get_functions():
                                        if hit_function == read_function:
                                            # Count hits
                                            tax_functions_count += 1.0
                                            # Add function to non-redundant list of functions
                                            function_list.add(read_function)
                                            # Count hits per function per taxid
                                            function_taxids[read_function][hit_taxid] += 1.0
                                            # Count identity and hits (per function) here
                                            scores[hit_taxid][read_function]['hit_count'] += 1.0
                                            scores[hit_taxid][read_function]['identity'] += hit.get_identity()

                            # Second, for each tax ID add its share of read count and RPKM score assigned to hit's function
                            for function in function_taxids:
                                tax_count = sum(function_taxids[function].values()) # How many times this function was counted in all tax_ids
                                for hit_taxid in function_taxids[function]:
                                    scores[hit_taxid][function]['fpkm'] += read_functions[function] * function_taxids[function][hit_taxid] / tax_count
                                    scores[hit_taxid][function]['count'] += function_taxids[function][hit_taxid] / tax_functions_count

            tax_profile = TaxonomyProfile()
            
            #outfile = os.path.join(self.parser.project.get_project_dir(sample), self.parser.project.get_output_subdir(sample), sample + '_' + end + '_' + 'functional_taxonomy_profile.xml')
            outfile = os.path.join(sample + '_fpkm_functional_taxonomy_profile.xml')
            tax_profile.build_functional_taxonomy_profile(tax_data, scores)
            #print(tax_profile.print_functional_taxonomy_profile())
            
            
            
            
            
            #function_list = ['GH5_4','GH10','GH55', 'GH1', 'GH6', 'GH9']
            function_list = sorted(all_functions_list)
            
            
            
            
            generate_functional_taxonomy_chart(tax_profile, sorted(function_list), outfile, 'fpkm')
            #print(tax_profile.print_functional_taxonomy_table())
            df = tax_profile.convert_function_taxonomic_profile_into_df('fpkm')
        
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
            for end in ENDS:
                project.samples[sample][end] = None
        
        # Add FPKM worksheet
        bold = workbook.add_format({'bold': True})
        scores_worksheet = workbook.add_worksheet('Functions FPKM')
        row = 0
        col = 0
        scores_worksheet.write(row, col, 'Function', bold)
        for sample in sorted(project.samples.keys()):
            col += 1
            scores_worksheet.write(row, col, sample, bold)
        col += 1
        scores_worksheet.write(row, col, 'Definition', bold)
        for function in sorted(function_list):
            row += 1
            col = 0
            scores_worksheet.write(row, col, function, bold)
            for sample in sorted(project.samples.keys()):
                col += 1
                if function in fpkm_scores[sample]:
                    scores_worksheet.write(row, col, '{0:.3f}'.format(fpkm_scores[sample][function]))
                else:
                    scores_worksheet.write(row, col, 0.0)
            col += 1
            scores_worksheet.write(row, col, project.ref_data.lookup_function_name(function))
        # adjust column width
        scores_worksheet.set_column(0, 0, 10)
        scores_worksheet.set_column(col, col, 50)

        writer.save()
        
        self.assertTrue(tax_profile.tree.data)

#    @unittest.skip("for faster testing")
    def test_7_lazy_build_fpkm_naive_functions_xslx(self):
        tax_data = TaxonomyData(self.parser.config)
        tax_data.load_taxdata(self.parser.config)
        
        project = Project(config_file=config_path, project_file=project_path)
        outfile = project.options.get_name() + '_fpkm_functions_taxonomy.xlsx'
        #outfile = 'test_functions_taxonomy.xlsx'
        outfile = outfile.replace(' ', '_')
        outfile = outfile.replace("'", "")
        writer = pd.ExcelWriter(outfile, engine='xlsxwriter')
        ENDS = ('pe1', 'pe2')
        fpkm_scores = autovivify(2, float)
        all_functions_list = set()
        for sample in project.list_samples():
            for end in ENDS:
                project.load_annotated_reads(sample, end)

            scores = autovivify(3, float)
            reads_processed = set()
            function_list = set()

            for read_id in project.samples[sample]['pe1'].keys():
                read_pe1 = project.samples[sample]['pe1'][read_id]
                if read_pe1.get_status() == 'function,besthit' or read_pe1.get_status() == 'function':
                    reads_processed.add(read_id)
                                            
                    if 'pe2' in project.samples[sample] and read_id in project.samples[sample]['pe2'].keys(): 
                        read_pe2 = project.samples[sample]['pe2'][read_id]
                        if read_pe2.get_status() == 'function,besthit' or read_pe2.get_status() == 'function':
                            # Both ends are mapped
                            fpkm = defaultdict(float) # Stores FPKM scores for the current read
                            read_functions = set() # List of functions assigned to the current read
                            if self.have_similar_functions(read_pe1, read_pe2): # Ends have overlapping functions, count them together

                                read1_functions = read_pe1.get_functions()
                                read2_functions = read_pe2.get_functions()

                                read_functions.update(read1_functions.keys())
                                read_functions.update(read2_functions.keys())
                                # We filled list of functions. Let's calculate FPKM score
                                
                                for read_function in read_functions:
                                    if read_function in read1_functions and read_function in read2_functions: # Take higher score
                                        if read1_functions[read_function] > read2_functions[read_function]:
                                            fpkm[read_function] = read1_functions[read_function]
                                        else:
                                            fpkm[read_function] = read2_functions[read_function]
                                    elif read_function in read1_functions:
                                        fpkm[read_function] += read1_functions[read_function]
                                    elif read_function in read2_functions:
                                        fpkm[read_function] += read2_functions[read_function]
                            else:
                                read_functions = read_pe1.get_functions()
                                for read_function in read_functions:
                                    fpkm[read_function] += read_functions[read_function]
                                read_functions = read_pe2.get_functions()
                                for read_function in read_functions:
                                    fpkm[read_function] += read_functions[read_function]
                            
                            # Add FPKM score of the read to FPKM score of the sample
                            for function in fpkm:
                                all_functions_list.add(function)
                                fpkm_scores[sample][function] += fpkm[function]

                    else: #Only end1 is mapped
                        read_functions = read_pe1.get_functions()
                        # Count FPKM
                        for function in read_functions:
                            all_functions_list.add(function)
                            fpkm_scores[sample][function] += read_functions[function]
            for read_id in project.samples[sample]['pe2'].keys():
                if read_id not in reads_processed: #Only end2 was mapped
                    read_pe2 = project.samples[sample]['pe2'][read_id]
                    if read_pe2.get_status() == 'function,besthit' or read_pe2.get_status() == 'function':
                        read_functions = read_pe2.get_functions()
                        # Count FPKM
                        for function in read_functions:
                            all_functions_list.add(function)
                            fpkm_scores[sample][function] += read_functions[function]
            
            # Delete reads from memory
            for end in ENDS:
                project.samples[sample][end] = None

        workbook  = writer.book

        # Add FPKM worksheet
        bold = workbook.add_format({'bold': True})
        scores_worksheet = workbook.add_worksheet('Functions FPKM')
        row = 0
        col = 0
        scores_worksheet.write(row, col, 'Function', bold)
        for sample in sorted(project.list_samples()):
            col += 1
            scores_worksheet.write(row, col, sample, bold)
        col += 1
        scores_worksheet.write(row, col, 'Definition', bold)
        for function in sorted(all_functions_list):
            row += 1
            col = 0
            scores_worksheet.write(row, col, function, bold)
            for sample in sorted(project.list_samples()):
                col += 1
                if function in fpkm_scores[sample]:
                    scores_worksheet.write(row, col, fpkm_scores[sample][function])
                else:
                    scores_worksheet.write(row, col, 0.0)
            col += 1
            scores_worksheet.write(row, col, project.ref_data.lookup_function_name(function))
        # adjust column width
        scores_worksheet.set_column(0, 0, 10)
        scores_worksheet.set_column(col, col, 50)

        writer.save()
        
        self.assertTrue(fpkm_scores)

    @unittest.skip("for faster testing")
    def test_8_lazy_build_fpkm_taxonomy_profile(self):
        tax_data = TaxonomyData(self.parser.config)
        tax_data.load_taxdata(self.parser.config)
        
        project = Project(config_file=config_path, project_file=project_path)
        outfile = project.options.get_name() + '_fpkm_functions_taxonomy.xlsx'
        #outfile = 'test_functions_taxonomy.xlsx'
        outfile = outfile.replace(' ', '_')
        outfile = outfile.replace("'", "")
        writer = pd.ExcelWriter(outfile, engine='xlsxwriter')
        ENDS = ('pe1', 'pe2')
        fpkm_scores = autovivify(2, float)
        all_functions_list = set()
        for sample in project.list_samples():
            for end in ENDS:
                project.load_annotated_reads(sample, end)

            scores = autovivify(3, float)
            reads_processed = set()
            function_list = set()

            for read_id in project.samples[sample]['pe1'].keys():
                read_pe1 = project.samples[sample]['pe1'][read_id]
                if read_pe1.get_status() == 'function,besthit' or read_pe1.get_status() == 'function':
                    reads_processed.add(read_id)
                                            
                    if 'pe2' in project.samples[sample] and read_id in project.samples[sample]['pe2'].keys(): 
                        read_pe2 = project.samples[sample]['pe2'][read_id]
                        if read_pe2.get_status() == 'function,besthit' or read_pe2.get_status() == 'function':
                            # Both ends are mapped
                            fpkm = defaultdict(float) # Stores FPKM scores for the current read
                            read_functions = set() # List of functions assigned to the current read
                            if self.have_similar_functions(read_pe1, read_pe2): # Ends have overlapping functions, count them together

                                read1_functions = read_pe1.get_functions()
                                read2_functions = read_pe2.get_functions()

                                read_functions.update(read1_functions.keys())
                                read_functions.update(read2_functions.keys())
                                # We filled list of functions. Let's calculate FPKM score
                                
                                for read_function in read_functions:
                                    if read_function in read1_functions and read_function in read2_functions: # Take higher score
                                        if read1_functions[read_function] > read2_functions[read_function]:
                                            fpkm[read_function] = read1_functions[read_function]
                                        else:
                                            fpkm[read_function] = read2_functions[read_function]
                                    elif read_function in read1_functions:
                                        fpkm[read_function] += read1_functions[read_function]
                                    elif read_function in read2_functions:
                                        fpkm[read_function] += read2_functions[read_function]
                            else:
                                read_functions = read_pe1.get_functions()
                                for read_function in read_functions:
                                    fpkm[read_function] += read_functions[read_function]
                                read_functions = read_pe2.get_functions()
                                for read_function in read_functions:
                                    fpkm[read_function] += read_functions[read_function]
                            
                            # Add FPKM score of the read to FPKM score of the sample
                            for function in fpkm:
                                all_functions_list.add(function)
                                fpkm_scores[sample][function] += fpkm[function]

                            # Build FPKM-based functional taxonomic profile
                            function_maxbitscores = {}
                            hits1 = [hit for hit in read_pe1.get_hit_list().get_hits() if read_pe1.get_status() == 'function,besthit']
                            hits2 = [hit for hit in read_pe2.get_hit_list().get_hits() if read_pe2.get_status() == 'function,besthit']
                            # Find max. bitscore for each function
                            for hit in hits1:
                                for hit_function in hit.get_functions():
                                    if hit_function in read_functions:
                                        if hit_function in function_maxbitscores:
                                            if hit.get_bitscore() > function_maxbitscores[hit_function]:
                                                function_maxbitscores[hit_function] = hit.get_bitscore()
                                        else:
                                            function_maxbitscores[hit_function] = hit.get_bitscore()
                            for hit in hits2:
                                for hit_function in hit.get_functions():
                                    if hit_function in read_functions:
                                        if hit_function in function_maxbitscores:
                                            if hit.get_bitscore() > function_maxbitscores[hit_function]:
                                                function_maxbitscores[hit_function] = hit.get_bitscore()
                                        else:
                                            function_maxbitscores[hit_function] = hit.get_bitscore()
                            # Count hits, FPKM scores, fragments
                            for hit in hits1:
                                hit_taxid = self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                                for hit_function in hit.get_functions():
                                    if hit_function in read_pe1.get_functions():
                                        if hit_function in function_maxbitscores and hit.get_bitscore() == function_maxbitscores[hit_function]:
                                            scores[hit_taxid][hit_function]['hit_count'] += 1.0 
                                            scores[hit_taxid][hit_function]['identity'] += hit.get_identity()
                                            scores[hit_taxid][hit_function]['fpkm'] += read_pe1.get_functions()[hit_function]
                                            scores[hit_taxid][hit_function]['count'] += 1.0
                                            del function_maxbitscores[hit_function]

                            for hit in hits2:
                                hit_taxid = self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                                for hit_function in hit.get_functions():
                                    if hit_function in read_pe2.get_functions():
                                        if hit_function in function_maxbitscores and hit.get_bitscore() == function_maxbitscores[hit_function]:
                                            scores[hit_taxid][hit_function]['hit_count'] += 1.0 
                                            scores[hit_taxid][hit_function]['identity'] += hit.get_identity()
                                            scores[hit_taxid][hit_function]['fpkm'] += read_pe2.get_functions()[hit_function]
                                            scores[hit_taxid][hit_function]['count'] += 1.0
                                            del function_maxbitscores[hit_function]


                    else: #Only end1 is mapped
                        read_functions = read_pe1.get_functions()
                        # Count FPKM
                        for function in read_functions:
                            all_functions_list.add(function)
                            fpkm_scores[sample][function] += read_functions[function]


                        # Build FPKM-based functional taxonomic profile
                        function_maxbitscores = {}
                        hits1 = [hit for hit in read_pe1.get_hit_list().get_hits() if read_pe1.get_status() == 'function,besthit']
                        # Find max. bitscore for each function
                        for hit in hits1:
                            for hit_function in hit.get_functions():
                                if hit_function in read_functions:
                                    if hit_function in function_maxbitscores:
                                        if hit.get_bitscore() > function_maxbitscores[hit_function]:
                                            function_maxbitscores[hit_function] = hit.get_bitscore()
                                    else:
                                        function_maxbitscores[hit_function] = hit.get_bitscore()
                        # Count hits, FPKM scores, fragments
                        for hit in hits1:
                            hit_taxid = self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                            for hit_function in hit.get_functions():
                                if hit_function in read_functions:
                                    if hit_function in function_maxbitscores and hit.get_bitscore() == function_maxbitscores[hit_function]:
                                        scores[hit_taxid][hit_function]['hit_count'] += 1.0 
                                        scores[hit_taxid][hit_function]['identity'] += hit.get_identity()
                                        scores[hit_taxid][hit_function]['fpkm'] += read_pe1.get_functions()[hit_function]
                                        scores[hit_taxid][hit_function]['count'] += 1.0
                                        del function_maxbitscores[hit_function]

            for read_id in project.samples[sample]['pe2'].keys():
                if read_id not in reads_processed: #Only end2 was mapped
                    read_pe2 = project.samples[sample]['pe2'][read_id]
                    if read_pe2.get_status() == 'function,besthit' or read_pe2.get_status() == 'function':
                        read_functions = read_pe2.get_functions()
                        # Count FPKM
                        for function in read_functions:
                            all_functions_list.add(function)
                            fpkm_scores[sample][function] += read_functions[function]

                        # Build FPKM-based functional taxonomic profile
                        function_maxbitscores = {}
                        hits2 = [hit for hit in read_pe2.get_hit_list().get_hits() if read_pe2.get_status() == 'function,besthit']
                        # Find max. bitscore for each function
                        for hit in hits2:
                            for hit_function in hit.get_functions():
                                if hit_function in read_functions:
                                    if hit_function in function_maxbitscores:
                                        if hit.get_bitscore() > function_maxbitscores[hit_function]:
                                            function_maxbitscores[hit_function] = hit.get_bitscore()
                                    else:
                                        function_maxbitscores[hit_function] = hit.get_bitscore()
                        # Count hits, FPKM scores, fragments
                        for hit in hits2:
                            hit_taxid = self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                            for hit_function in hit.get_functions():
                                if hit_function in read_functions:
                                    if hit_function in function_maxbitscores and hit.get_bitscore() == function_maxbitscores[hit_function]:
                                        scores[hit_taxid][hit_function]['hit_count'] += 1.0 
                                        scores[hit_taxid][hit_function]['identity'] += hit.get_identity()
                                        scores[hit_taxid][hit_function]['fpkm'] += read_pe2.get_functions()[hit_function]
                                        scores[hit_taxid][hit_function]['count'] += 1.0
                                        del function_maxbitscores[hit_function]
            
            tax_profile = TaxonomyProfile()
            
            #outfile = os.path.join(self.parser.project.get_project_dir(sample), self.parser.project.get_output_subdir(sample), sample + '_' + end + '_' + 'functional_taxonomy_profile.xml')
            outfile = os.path.join(sample + '_fpkm_naive_functional_taxonomy_profile.xml')
            tax_profile.build_functional_taxonomy_profile(tax_data, scores)
            #print(tax_profile.print_functional_taxonomy_profile())
            
            function_list = sorted(all_functions_list)
            
            generate_functional_taxonomy_chart(tax_profile, sorted(function_list), outfile, 'fpkm')
            #print(tax_profile.print_functional_taxonomy_table())
            df = tax_profile.convert_function_taxonomic_profile_into_df('fpkm')
        
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

            # Delete reads from memory
            for end in ENDS:
                project.samples[sample][end] = None


        # Add FPKM worksheet
        bold = workbook.add_format({'bold': True})
        scores_worksheet = workbook.add_worksheet('Functions FPKM')
        row = 0
        col = 0
        scores_worksheet.write(row, col, 'Function', bold)
        for sample in sorted(project.list_samples()):
            col += 1
            scores_worksheet.write(row, col, sample, bold)
        col += 1
        scores_worksheet.write(row, col, 'Definition', bold)
        for function in sorted(all_functions_list):
            row += 1
            col = 0
            scores_worksheet.write(row, col, function, bold)
            for sample in sorted(project.list_samples()):
                col += 1
                if function in fpkm_scores[sample]:
                    scores_worksheet.write(row, col, fpkm_scores[sample][function])
                else:
                    scores_worksheet.write(row, col, 0.0)
            col += 1
            scores_worksheet.write(row, col, project.ref_data.lookup_function_name(function))
        # adjust column width
        scores_worksheet.set_column(0, 0, 10)
        scores_worksheet.set_column(col, col, 50)

        writer.save()
        
        self.assertTrue(fpkm_scores)

    def tearDown(self):
        self.parser = None
        
    def have_similar_functions(self, read1, read2):
        ret_val = False
        read1_functions = read1.get_functions()
        for function in read1.get_functions():
            if function in read1_functions:
                ret_val = True
        return ret_val

if __name__=='__main__':
    unittest.main()
