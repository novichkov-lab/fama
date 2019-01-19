#!/usr/bin/python3
import os, csv, operator
import unittest
from collections import Counter,defaultdict

import xlsxwriter
import pandas as pd

from context import Fama
from Fama.utils import autovivify,cleanup_protein_id,sanitize_file_name
from Fama.Project import Project
from Fama.DiamondParser.DiamondParser import DiamondParser
from Fama.TaxonomyProfile import TaxonomyProfile

from Fama.OutputUtil.JSONUtil import import_annotated_reads
from Fama.OutputUtil.Report import generate_fastq_report
from Fama.OutputUtil.KronaXMLWriter import generate_functional_taxonomy_chart,generate_taxonomy_series_chart

data_dir = 'data'
config_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'config.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project.ini')
sample_id = 'test_sample'
end = 'pe1'

#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_EB271_nitrogen8_t.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_fw106_fw301_sulfate.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW306_cazy1_t.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_B1F_cazy_t.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'protein_project_EB271_assembly.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW306_nitrogen8test_t.ini')
#sample_id = 'sample1'
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW3062M_universal1.ini')
sample_id = 'sample6'

class FunctionTaxonomyProfilingTest(unittest.TestCase):

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
    def test_1_build_functional_taxonomy_profile(self):
        
        project = Project(config_file=config_path, project_file=project_path)
        #project.load_functional_profile()
        outfile = self.project.options.get_name() + '_functions_taxonomy.xlsx'
        #outfile = 'test_functions_taxonomy.xlsx'
        outfile = outfile.replace(' ', '_')
        outfile = outfile.replace("'", "")

        writer = pd.ExcelWriter(outfile, engine='xlsxwriter')

        for sample in sorted(project.list_samples()):
            scores = defaultdict(lambda : defaultdict(dict))
            function_list = set()
            self.project.import_reads_json(sample,self.project.ENDS)
            for end in self.project.ENDS:

                scaling_factor = 1.0
                if end == 'pe1':
                    scaling_factor = self.project.options.get_fastq1_readcount(sample)/(self.project.options.get_fastq1_readcount(sample) + self.project.options.get_fastq2_readcount(sample))
                elif end == 'pe2':
                    scaling_factor = self.project.options.get_fastq2_readcount(sample)/(self.project.options.get_fastq1_readcount(sample) + self.project.options.get_fastq2_readcount(sample))
                else:
                    raise Exception('Unknown end identifier')

#        self.parser.reads = import_annotated_reads(os.path.join(self.parser.options.get_project_dir(self.parser.sample), self.parser.sample + '_' + self.parser.end + '_' + self.parser.options.get_reads_json_name()))
#        reads = self.parser.reads

                multiple_hits = 0
                read_count = 0
                
                for read_id,read in self.project.samples[sample].reads[end].items():
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
                self.project.samples[sample].reads[end] = None
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
            
            #outfile = os.path.join(self.parser.options.get_project_dir(sample), self.parser.options.get_output_subdir(sample), sample + '_' + end + '_' + 'functional_taxonomy_profile.xml')
            outfile = os.path.join(sample + '_' + 'functional_taxonomy_profile.xml')
            tax_profile.build_functional_taxonomy_profile(self.project.taxonomy_data, scores)
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
            

#    @unittest.skip("temporary excluded")
    def test_2_build_fpkm_taxonomy_profile(self):
        
        project = Project(config_file=config_path, project_file=project_path)
        outfile = self.project.options.get_name() + '_fpkm_functions_taxonomy.xlsx'
        #outfile = 'test_functions_taxonomy.xlsx'
        outfile = outfile.replace(' ', '_')
        outfile = outfile.replace("'", "")
        writer = pd.ExcelWriter(outfile, engine='xlsxwriter')
        fpkm_scores = autovivify(2, float)
        all_functions_list = set()
        for sample in sorted(self.project.list_samples()):
            self.project.import_reads_json(sample,self.project.ENDS)
            scores = autovivify(3, float)
            reads_processed = set()
            function_list = set()

            reads_pe1 = self.project.samples[sample].reads['pe1']
            reads_pe2 = {}
            if 'pe2' in self.project.samples[sample].reads:
                reads_pe2 = self.project.samples[sample].reads['pe2']
                
            
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
                                    min_bitscore = max_bitscore * (1 - self.project.config.get_biscore_range_cutoff(self.project.options.get_collection()))
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
            
            #outfile = os.path.join(self.parser.options.get_project_dir(sample), self.parser.options.get_output_subdir(sample), sample + '_' + end + '_' + 'functional_taxonomy_profile.xml')
            outfile = os.path.join(sample + '_fpkm_functional_taxonomy_profile.xml')
            tax_profile.build_functional_taxonomy_profile(self.project.taxonomy_data, scores)
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
        for sample in sorted(self.project.list_samples()):
            col += 1
            scores_worksheet.write(row, col, sample, bold)
        col += 1
        scores_worksheet.write(row, col, 'Definition', bold)
        for function in sorted(all_functions_list):
            row += 1
            col = 0
            scores_worksheet.write(row, col, function, bold)
            for sample in sorted(self.project.list_samples()):
                col += 1
                if function in fpkm_scores[sample]:
                    scores_worksheet.write(row, col, '{0:.3f}'.format(fpkm_scores[sample][function]))
                else:
                    scores_worksheet.write(row, col, 0.0)
            col += 1
            scores_worksheet.write(row, col, self.project.ref_data.lookup_function_name(function))
        # adjust column width
        scores_worksheet.set_column(0, 0, 10)
        scores_worksheet.set_column(col, col, 50)

        writer.save()
        
        self.assertTrue(tax_profile.tree.data)

#    @unittest.skip("for faster testing")
    def test_3_lazy_build_fpkm_taxonomy_profile(self):
        
        outfile = sanitize_file_name(self.project.options.get_name() + '_fpkm_functions_taxonomy.xlsx')
        #outfile = 'test_functions_taxonomy.xlsx'
        writer = pd.ExcelWriter(outfile, engine='xlsxwriter')
        fpkm_scores = autovivify(2, float)
        all_functions_list = set()
        for sample in self.project.list_samples():
            self.project.import_reads_json(sample,self.project.ENDS)

            scores = autovivify(3, float)
            reads_processed = set()
            function_list = set()

            reads_pe1 = self.project.samples[sample].reads['pe1']
            reads_pe2 = {}
            if 'pe2' in self.project.samples[sample].reads:
                reads_pe2 = self.project.samples[sample].reads['pe2']
            
            for read in reads_pe1:
                if reads_pe1[read].get_status() == 'function':
                    reads_processed.add(read)
                                            
                    if read in reads_pe2: 
                        if reads_pe2[read].get_status() == 'function':
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
                                if read_function in read1_functions and reads_pe1[read].get_status() == 'function':
                                    for hit in hits1:
                                        for hit_function in hit.get_functions():
                                            if hit_function == read_function and hit.get_bitscore() > max_bitscore:
                                                function_list.add(read_function)
                                                max_bitscore = hit.get_bitscore()
                                if read_function in read2_functions and reads_pe2[read].get_status() == 'function':
                                    for hit in hits2:
                                        for hit_function in hit.get_functions():
                                            if hit_function == read_function and hit.get_bitscore() > max_bitscore:
                                                function_list.add(read_function)
                                                max_bitscore = hit.get_bitscore()

                                # Ignore all hits for this function with bitscore < cutoff. 

                                if max_bitscore > 0.0:
                                    min_bitscore = max_bitscore * (1 - self.project.config.get_biscore_range_cutoff(self.project.options.get_collection()))
                                    for hit in (hit for hit in hits1 if reads_pe1[read].get_status() == 'function'):
                                        hit_taxid = self.parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                                        for hit_function in (hit_function for hit_function in hit.get_functions() if hit.get_bitscore() > min_bitscore and hit_function == read_function):
                                            # Count hits
                                            tax_functions_count += 1.0
                                            # Count hits per function per taxid
                                            function_taxids[read_function][hit_taxid] += 1.0
                                            # Count identity and hits (per function) 
                                            scores[hit_taxid][read_function]['hit_count'] += 1.0 
                                            scores[hit_taxid][read_function]['identity'] += hit.get_identity()

                                    for hit in (hit for hit in hits2 if reads_pe2[read].get_status() == 'function'):
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
                        
                        if reads_pe1[read].get_status() == 'function':
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
                    if reads_pe2[read].get_status() == 'function':
                        read_functions = reads_pe2[read].get_functions()
                        hits = reads_pe2[read].get_hit_list().get_hits()
                        # Count FPKM
                        for function in read_functions:
                            all_functions_list.add(function)
                            fpkm_scores[sample][function] += read_functions[function]

                        tax_functions_count = 0.0
                        function_taxids = autovivify(2, float) # for each function, store list of tax_ids and count of hits for each tax_id

                        if reads_pe2[read].get_status() == 'function':
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
            
            #outfile = os.path.join(self.parser.options.get_project_dir(sample), self.parser.options.get_output_subdir(sample), sample + '_' + end + '_' + 'functional_taxonomy_profile.xml')
            outfile = os.path.join(sample + '_fpkm_functional_taxonomy_profile.xml')
            tax_profile.build_functional_taxonomy_profile(self.project.taxonomy_data, scores)
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
            for end in self.project.ENDS:
                self.project.samples[sample].reads[end] = None
        
        # Add FPKM worksheet
        bold = workbook.add_format({'bold': True})
        scores_worksheet = workbook.add_worksheet('Functions FPKM')
        row = 0
        col = 0
        scores_worksheet.write(row, col, 'Function', bold)
        for sample in sorted(self.project.list_samples()):
            col += 1
            scores_worksheet.write(row, col, sample, bold)
        col += 1
        scores_worksheet.write(row, col, 'Definition', bold)
        for function in sorted(function_list):
            row += 1
            col = 0
            scores_worksheet.write(row, col, function, bold)
            for sample in sorted(self.project.list_samples()):
                col += 1
                if function in fpkm_scores[sample]:
                    scores_worksheet.write(row, col, '{0:.3f}'.format(fpkm_scores[sample][function]))
                else:
                    scores_worksheet.write(row, col, 0.0)
            col += 1
            scores_worksheet.write(row, col, self.project.ref_data.lookup_function_name(function))
        # adjust column width
        scores_worksheet.set_column(0, 0, 10)
        scores_worksheet.set_column(col, col, 50)

        writer.save()
        
        self.assertTrue(tax_profile.tree.data)

#    @unittest.skip("for faster testing")
    def test_4_lazy_build_fpkm_naive_functions_xslx(self):
        
        outfile = self.project.options.get_name() + '_fpkm_functions_taxonomy.xlsx'
        #outfile = 'test_functions_taxonomy.xlsx'
        outfile = outfile.replace(' ', '_')
        outfile = outfile.replace("'", "")
        writer = pd.ExcelWriter(outfile, engine='xlsxwriter')
        fpkm_scores = autovivify(2, float)
        all_functions_list = set()
        for sample in self.project.list_samples():
            self.project.import_reads_json(sample,self.project.ENDS)

            scores = autovivify(3, float)
            reads_processed = set()
            function_list = set()

            for read_id in self.project.samples[sample].reads['pe1'].keys():
                read_pe1 = self.project.samples[sample].reads['pe1'][read_id]
                if read_pe1.get_status() == 'function,besthit' or read_pe1.get_status() == 'function':
                    reads_processed.add(read_id)
                                            
                    if 'pe2' in self.project.samples[sample].reads and read_id in self.project.samples[sample].reads['pe2'].keys(): 
                        read_pe2 = self.project.samples[sample].reads['pe2'][read_id]
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
            for read_id in self.project.samples[sample].reads['pe2'].keys():
                if read_id not in reads_processed: #Only end2 was mapped
                    read_pe2 = self.project.samples[sample].reads['pe2'][read_id]
                    if read_pe2.get_status() == 'function,besthit' or read_pe2.get_status() == 'function':
                        read_functions = read_pe2.get_functions()
                        # Count FPKM
                        for function in read_functions:
                            all_functions_list.add(function)
                            fpkm_scores[sample][function] += read_functions[function]
            
            # Delete reads from memory
            for end in self.project.ENDS:
                self.project.samples[sample].reads[end] = None

        workbook  = writer.book

        # Add FPKM worksheet
        bold = workbook.add_format({'bold': True})
        scores_worksheet = workbook.add_worksheet('Functions FPKM')
        row = 0
        col = 0
        scores_worksheet.write(row, col, 'Function', bold)
        for sample in sorted(self.project.list_samples()):
            col += 1
            scores_worksheet.write(row, col, sample, bold)
        col += 1
        scores_worksheet.write(row, col, 'Definition', bold)
        for function in sorted(all_functions_list):
            row += 1
            col = 0
            scores_worksheet.write(row, col, function, bold)
            for sample in sorted(self.project.list_samples()):
                col += 1
                if function in fpkm_scores[sample]:
                    scores_worksheet.write(row, col, fpkm_scores[sample][function])
                else:
                    scores_worksheet.write(row, col, 0.0)
            col += 1
            scores_worksheet.write(row, col, self.project.ref_data.lookup_function_name(function))
        # adjust column width
        scores_worksheet.set_column(0, 0, 10)
        scores_worksheet.set_column(col, col, 50)

        writer.save()
        
        self.assertTrue(fpkm_scores)

#    @unittest.skip("for faster testing")
    def test_5_lazy_build_fpkm_taxonomy_profile(self):
        
        outfile = sanitize_file_name(self.project.options.get_name() + '_fpkm_functions_taxonomy.xlsx')
        writer = pd.ExcelWriter(outfile, engine='xlsxwriter')
        fpkm_scores = autovivify(2, float)
        all_functions_list = set()
        for sample in self.project.list_samples():
            self.project.import_reads_json(sample,self.project.ENDS)

            scores = autovivify(3, float)
            reads_processed = set()
            function_list = set()

            for read_id in self.project.samples[sample].reads['pe1'].keys():
                read_pe1 = self.project.samples[sample].reads['pe1'][read_id]
                if read_pe1.get_status() == 'function,besthit' or read_pe1.get_status() == 'function':
                    reads_processed.add(read_id)
                                            
                    if 'pe2' in self.project.samples[sample].reads and read_id in self.project.samples[sample].reads['pe2'].keys(): 
                        read_pe2 = self.project.samples[sample].reads['pe2'][read_id]
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
                            hits1 = [hit for hit in read_pe1.get_hit_list().get_hits() if read_pe1.get_status() == 'function']
                            hits2 = [hit for hit in read_pe2.get_hit_list().get_hits() if read_pe2.get_status() == 'function']
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

            for read_id in self.project.samples[sample].reads['pe2'].keys():
                if read_id not in reads_processed: #Only end2 was mapped
                    read_pe2 = self.project.samples[sample].reads['pe2'][read_id]
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
            
            #outfile = os.path.join(self.parser.options.get_project_dir(sample), self.parser.options.get_output_subdir(sample), sample + '_' + end + '_' + 'functional_taxonomy_profile.xml')
            outfile = os.path.join(sample + '_fpkm_naive_functional_taxonomy_profile.xml')
            tax_profile.build_functional_taxonomy_profile(self.project.taxonomy_data, scores)
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
            for end in self.project.ENDS:
                self.project.samples[sample].reads[end] = None


        # Add FPKM worksheet
        bold = workbook.add_format({'bold': True})
        scores_worksheet = workbook.add_worksheet('Functions FPKM')
        row = 0
        col = 0
        scores_worksheet.write(row, col, 'Function', bold)
        for sample in sorted(self.project.list_samples()):
            col += 1
            scores_worksheet.write(row, col, sample, bold)
        col += 1
        scores_worksheet.write(row, col, 'Definition', bold)
        for function in sorted(all_functions_list):
            row += 1
            col = 0
            scores_worksheet.write(row, col, function, bold)
            for sample in sorted(self.project.list_samples()):
                col += 1
                if function in fpkm_scores[sample]:
                    scores_worksheet.write(row, col, fpkm_scores[sample][function])
                else:
                    scores_worksheet.write(row, col, 0.0)
            col += 1
            scores_worksheet.write(row, col, self.project.ref_data.lookup_function_name(function))
        # adjust column width
        scores_worksheet.set_column(0, 0, 10)
        scores_worksheet.set_column(col, col, 50)

        writer.save()
        
        self.assertTrue(fpkm_scores)

#    @unittest.skip("for faster testing")
    def test_6_lazy_build_fragment_function_taxonomy_profile(self):
        
        outfile = sanitize_file_name(self.project.options.get_name() + '_fpkm_functions_taxonomy.xlsx')
        writer = pd.ExcelWriter(outfile, engine='xlsxwriter')
        metrics = 'fpkm'



        all_scores = autovivify(2, float) # For function comparative table
        all_functions_list = set()
        for sample_id in self.project.list_samples():
            self.project.import_reads_json(sample_id,self.project.ENDS)


            if metrics == 'fpkm':
                norm_factor = self.project.samples[sample_id].rpkm_scaling_factor
                if norm_factor == 0.0:
                    norm_factor = 1000000/project.options.get_fastq1_readcount(sample_id)
            elif metrics == 'fpkg':
                norm_factor = self.project.samples[sample_id].rpkg_scaling_factor
                if norm_factor == 0.0:
                    raise ValueError('FPKG scaling factor is missing')


            sample_scores = autovivify(3, float) # For taxonomy profiling
            reads_processed = set()
            function_list = set()

            for read_id in self.project.samples[sample_id].reads['pe1'].keys():
                read_pe1 = self.project.samples[sample_id].reads['pe1'][read_id]
                if read_pe1.get_status() != 'function':
                    continue
                reads_processed.add(read_id)
                                        
                if 'pe2' in self.project.samples[sample_id].reads and read_id in self.project.samples[sample_id].reads['pe2'].keys(): 
                    read_pe2 = self.project.samples[sample_id].reads['pe2'][read_id]
                    if read_pe2.get_status() == 'function':
                        # Both ends are mapped
                        fragment_taxonomy = self.project.taxonomy_data.get_lca([read_pe1.taxonomy, read_pe2.taxonomy])
                        
                        #fragment_scores = defaultdict(float) # Stores FPKM scores for the current read
                        fragment_functions = set() # List of functions assigned to the current read

                        read1_functions = read_pe1.get_functions()
                        read2_functions = read_pe2.get_functions()

                        fragment_functions.update(read1_functions.keys())
                        fragment_functions.update(read2_functions.keys())

                        # We filled list of functions. Let's calculate FPKM score
                            
                        for function in fragment_functions:
                            all_functions_list.add(function)
                            sample_scores[fragment_taxonomy][function]['count'] += 1
                            if function in read1_functions and function in read2_functions: # Take higher score
                                all_scores[sample_id][function] += max (norm_factor * read1_functions[function], norm_factor * read2_functions[function])
                                sample_scores[fragment_taxonomy][function][metrics] += max (norm_factor * read1_functions[function], norm_factor * read2_functions[function])
                            elif function in read1_functions:
                                all_scores[sample_id][function] += norm_factor * read1_functions[function]
                                sample_scores[fragment_taxonomy][function][metrics] += norm_factor * read1_functions[function]
                            elif function in read2_functions:
                                all_scores[sample_id][function] += norm_factor * read2_functions[function]
                                sample_scores[fragment_taxonomy][function][metrics] += norm_factor * read2_functions[function]

                        # Let's get hits data for calculation of identity% 
                        
                        function_maxbitscores = {}
                        hits1 = [hit for hit in read_pe1.get_hit_list().get_hits() if read_pe1.get_status() == 'function']
                        hits2 = [hit for hit in read_pe2.get_hit_list().get_hits() if read_pe2.get_status() == 'function']
                        # Find max. bitscore for each function
                        for hit in hits1:
                            for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                                if hit_function in function_maxbitscores:
                                    if hit.get_bitscore() > function_maxbitscores[hit_function]:
                                        function_maxbitscores[hit_function] = hit.get_bitscore()
                                else:
                                    function_maxbitscores[hit_function] = hit.get_bitscore()
                        for hit in hits2:
                            for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                                if hit_function in function_maxbitscores:
                                    if hit.get_bitscore() > function_maxbitscores[hit_function]:
                                        function_maxbitscores[hit_function] = hit.get_bitscore()
                                else:
                                    function_maxbitscores[hit_function] = hit.get_bitscore()
                        # Count hits and sum identity
                        for hit in hits1:
                            for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                                if hit_function in function_maxbitscores and hit.get_bitscore() == function_maxbitscores[hit_function]:
                                    sample_scores[fragment_taxonomy][hit_function]['hit_count'] += 1.0 
                                    sample_scores[fragment_taxonomy][hit_function]['identity'] += hit.get_identity()
                                    del function_maxbitscores[hit_function]

                        for hit in hits2:
                            for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                                if hit_function in function_maxbitscores and hit.get_bitscore() == function_maxbitscores[hit_function]:
                                    sample_scores[fragment_taxonomy][hit_function]['hit_count'] += 1.0 
                                    sample_scores[fragment_taxonomy][hit_function]['identity'] += hit.get_identity()
                                    del function_maxbitscores[hit_function]


                else: #Only end1 is mapped
                    fragment_functions = set()
                    fragment_taxonomy = read_pe1.taxonomy
                    read_functions = read_pe1.get_functions()
                    fragment_functions.update(read_functions.keys())
                    # Count FPKM
                    for function in fragment_functions:
                        all_functions_list.add(function)
                        all_scores[sample_id][function] += norm_factor * read_functions[function]
                        sample_scores[fragment_taxonomy][function]['count'] += 1
                        sample_scores[fragment_taxonomy][function][metrics] += norm_factor * read_functions[function]
                    # Build FPKM-based functional taxonomic profile
                    function_maxbitscores = {}
                    hits1 = [hit for hit in read_pe1.get_hit_list().get_hits() if read_pe1.get_status() == 'function']
                    # Find max. bitscore for each function
                    for hit in hits1:
                        for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                            if hit_function in function_maxbitscores:
                                if hit.get_bitscore() > function_maxbitscores[hit_function]:
                                    function_maxbitscores[hit_function] = hit.get_bitscore()
                            else:
                                function_maxbitscores[hit_function] = hit.get_bitscore()
                    # Count hits and sum identity
                    for hit in hits1:
                        for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                            if hit_function in function_maxbitscores and hit.get_bitscore() == function_maxbitscores[hit_function]:
                                sample_scores[fragment_taxonomy][hit_function]['hit_count'] += 1.0 
                                sample_scores[fragment_taxonomy][hit_function]['identity'] += hit.get_identity()
                                del function_maxbitscores[hit_function]

            for read_id in self.project.samples[sample_id].reads['pe2'].keys():
                if read_id in reads_processed: 
                    continue #Only end2 was mapped
                    
                read_pe2 = self.project.samples[sample_id].reads['pe2'][read_id]
                if read_pe2.get_status() != 'function':
                    continue

                fragment_functions = set()
                fragment_taxonomy = read_pe2.taxonomy
                read_functions = read_pe2.get_functions()
                fragment_functions.update(read_functions.keys())

                # Count FPKM
                for function in fragment_functions:
                    all_functions_list.add(function)
                    all_scores[sample_id][function] += norm_factor * read_functions[function]
                    sample_scores[fragment_taxonomy][function]['count'] += 1
                    sample_scores[fragment_taxonomy][function][metrics] += norm_factor * read_functions[function]
                # Build FPKM-based functional taxonomic profile
                function_maxbitscores = {}
                hits = [hit for hit in read_pe2.get_hit_list().get_hits() if read_pe2.get_status() == 'function']
                # Find max. bitscore for each function
                for hit in hits:
                    for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                        if hit_function in function_maxbitscores:
                            if hit.get_bitscore() > function_maxbitscores[hit_function]:
                                function_maxbitscores[hit_function] = hit.get_bitscore()
                        else:
                            function_maxbitscores[hit_function] = hit.get_bitscore()
                # Count hits and sum identity
                for hit in hits:
                    for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                        if hit_function in function_maxbitscores and hit.get_bitscore() == function_maxbitscores[hit_function]:
                            sample_scores[fragment_taxonomy][hit_function]['hit_count'] += 1.0 
                            sample_scores[fragment_taxonomy][hit_function]['identity'] += hit.get_identity()
                            del function_maxbitscores[hit_function]

            #~ for tax in sample_scores.keys():
                #~ for func in sample_scores[tax].keys():
                    #~ if not isinstance(sample_scores[tax][func]['identity'], float):
                        #~ print (tax, func)
            
            
            tax_profile = TaxonomyProfile()
            
            #outfile = os.path.join(self.parser.options.get_project_dir(sample), self.parser.options.get_output_subdir(sample), sample + '_' + end + '_' + 'functional_taxonomy_profile.xml')
#            outfile = os.path.join(sample + '_fpkm_naive_functional_taxonomy_profile.xml')
            outfile = sanitize_file_name(os.path.join(self.project.options.get_work_dir(), sample_id + '_' + metrics + '_fpkm_naive_functional_taxonomy_profile.xml'))
            tax_profile.build_functional_taxonomy_profile(self.project.taxonomy_data, sample_scores)
#            print(tax_profile.print_functional_taxonomy_profile(score=metrics))
            
            function_list = sorted(all_functions_list)
            generate_taxonomy_series_chart(tax_profile, function_list, outfile, score=metrics)

#            generate_functional_taxonomy_chart(tax_profile, function_list, outfile, 'fpkm')
#            print(tax_profile.print_functional_taxonomy_table())
            df = tax_profile.convert_function_taxonomic_profile_into_df(score=metrics)
        
            #print(df)
            df.to_excel(writer, sheet_name=sample_id)
            workbook  = writer.book
            worksheet = writer.sheets[sample_id]
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
            for end in self.project.ENDS:
                self.project.samples[sample_id].reads[end] = None


        # Add FPKM worksheet
        bold = workbook.add_format({'bold': True})
        scores_worksheet = workbook.add_worksheet('Functions FPKM')
        row = 0
        col = 0
        scores_worksheet.write(row, col, 'Function', bold)
        for sample in sorted(self.project.list_samples()):
            col += 1
            scores_worksheet.write(row, col, sample, bold)
        col += 1
        scores_worksheet.write(row, col, 'Definition', bold)
        for function in sorted(all_functions_list):
            row += 1
            col = 0
            scores_worksheet.write(row, col, function, bold)
            for sample in sorted(self.project.list_samples()):
                col += 1
                if function in all_scores[sample]:
                    scores_worksheet.write(row, col, all_scores[sample][function])
                else:
                    scores_worksheet.write(row, col, 0.0)
            col += 1
            scores_worksheet.write(row, col, self.project.ref_data.lookup_function_name(function))
        # adjust column width
        scores_worksheet.set_column(0, 0, 10)
        scores_worksheet.set_column(col, col, 50)

        writer.save()
        
        self.assertTrue(all_scores)

    def tearDown(self):
        self.parser = None
        
    def have_similar_functions(self, read1, read2):
        ret_val = False
        read1_functions = read1.get_functions()
        for function in read1.get_functions():
            if function in read1_functions:
                ret_val = True
        return ret_val

    def get_fragment_taxid(self, read1, read2):
        return self.project.taxonomy_data.get_lca([read1.taxonomy, read2.taxonomy])

if __name__=='__main__':
    unittest.main()
