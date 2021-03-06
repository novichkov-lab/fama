#!/usr/bin/python3
import os, csv, operator
import unittest
from collections import Counter,defaultdict

import xlsxwriter
import pandas as pd

from context import lib
from lib.utils.const import ENDS
from lib.utils.utils import autovivify,cleanup_protein_id,sanitize_file_name
from lib.project.project import Project
from lib.diamond_parser.diamond_parser import DiamondParser
from lib.taxonomy.taxonomy_profile import TaxonomyProfile

from lib.output.json_util import import_annotated_reads
from lib.output.report import generate_fastq_report, get_function_scores, get_function_taxonomy_scores
from lib.output.xlsx_util import make_function_sample_xlsx, make_func_tax_sample_xlsx, make_sample_tax_func_xlsx
from lib.output.krona_xml_writer import make_function_taxonomy_chart, make_taxonomy_series_chart

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
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW3062M_universal1.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW306_universal1_lca_test.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW306_nitrogen9_lca.ini')
#sample_id = 'sample6'

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
        
        outfile = 'test.xlsx'#sanitize_file_name(self.project.options.get_name() + '_fpkm_functions_taxonomy.xlsx')
        writer = pd.ExcelWriter(outfile, engine='xlsxwriter')
        metrics = 'fpkg'



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
                        fragment_lca = self.project.taxonomy_data.get_lca([read_pe1.taxonomy, read_pe2.taxonomy])
#                        if fragment_lca == '1':
#                            print ('Unknown LCA', read_id)
                        #fragment_scores = defaultdict(float) # Stores FPKM scores for the current read
                        fragment_functions = set() # List of functions assigned to the current read

                        read1_functions = read_pe1.get_functions()
                        read2_functions = read_pe2.get_functions()

                        fragment_functions.update(read1_functions.keys())
                        fragment_functions.update(read2_functions.keys())
                        
#                        if 'RP-S9' in fragment_functions:
#                            print(read_id, '\tboth ends\t', fragment_taxonomy)

                        # We filled list of functions. Let's calculate FPKM score
                            
                        for function in fragment_functions:
                            all_functions_list.add(function)
                            if function in read1_functions and function in read2_functions: # Take higher score
                                sample_scores[fragment_lca][function]['count'] += 1
                                all_scores[sample_id][function] += max (norm_factor * read1_functions[function], norm_factor * read2_functions[function])
                                sample_scores[fragment_lca][function][metrics] += max (norm_factor * read1_functions[function], norm_factor * read2_functions[function])
                            elif function in read1_functions:
                                sample_scores[read_pe1.taxonomy][function]['count'] += 1
                                #sample_scores[fragment_lca][function]['count'] += 1
                                all_scores[sample_id][function] += norm_factor * read1_functions[function]
                                sample_scores[read_pe1.taxonomy][function][metrics] += norm_factor * read1_functions[function]
                                #sample_scores[fragment_lca][function][metrics] += norm_factor * read1_functions[function]
                            elif function in read2_functions:
                                sample_scores[read_pe2.taxonomy][function]['count'] += 1
                                #sample_scores[fragment_lca][function]['count'] += 1
                                all_scores[sample_id][function] += norm_factor * read2_functions[function]
                                sample_scores[read_pe2.taxonomy][function][metrics] += norm_factor * read2_functions[function]
                                #sample_scores[fragment_lca][function][metrics] += norm_factor * read2_functions[function]

                        # Let's get hits data for calculation of identity% 
                        
                        function_maxbitscores = defaultdict(dict)
                        hits1 = [hit for hit in read_pe1.get_hit_list().get_hits()]
                        hits2 = [hit for hit in read_pe2.get_hit_list().get_hits()]
                        # Find max. bitscore for each function
                        for hit in hits1:
                            for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                                if hit_function in function_maxbitscores:
                                    if hit.get_bitscore() > function_maxbitscores[hit_function]['bitscore']:
                                        function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                        function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                        function_maxbitscores[hit_function]['taxonomy'] = read_pe1.taxonomy
                                else:
                                    function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                    function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                    function_maxbitscores[hit_function]['taxonomy'] = read_pe1.taxonomy
                        for hit in hits2:
                            for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                                if hit_function in function_maxbitscores:
                                    if hit.get_bitscore() > function_maxbitscores[hit_function]['bitscore']:
                                        function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                        function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                        function_maxbitscores[hit_function]['taxonomy'] = read_pe2.taxonomy
                                else:
                                    function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                    function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                    function_maxbitscores[hit_function]['taxonomy'] = read_pe2.taxonomy
                        # Count hits and sum identity
                        for function in fragment_functions:
                            if function in function_maxbitscores:
                                if function in read1_functions and function in read2_functions:
                                    #sample_scores[function_maxbitscores[function]['taxonomy']][function]['hit_count'] += 1.0 
                                    #sample_scores[function_maxbitscores[function]['taxonomy']][function]['identity'] += function_maxbitscores[hit_function]['identity'] 
                                    
                                    sample_scores[fragment_lca][function]['hit_count'] += 1.0 
                                    sample_scores[fragment_lca][function]['identity'] += function_maxbitscores[hit_function]['identity'] 
                                elif function in read1_functions:
                                    sample_scores[read_pe1.taxonomy][function]['hit_count'] += 1.0 
                                    sample_scores[read_pe1.taxonomy][function]['identity'] += function_maxbitscores[hit_function]['identity'] 
                                elif function in read2_functions:
                                    sample_scores[read_pe2.taxonomy][function]['hit_count'] += 1.0 
                                    sample_scores[read_pe2.taxonomy][function]['identity'] += function_maxbitscores[hit_function]['identity'] 
                                    
                            else:
                                print('Function',function,'not found in hits of read', read_id)
                                
                                
                        #~ for hit in hits1:
                            #~ for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                                #~ if hit_function in function_maxbitscores and hit.get_bitscore() == function_maxbitscores[hit_function]:
                                    #~ if hit_function in read1_functions and hit_function in read2_functions:
                                        #~ sample_scores[fragment_lca][hit_function]['hit_count'] += 1.0 
                                        #~ sample_scores[fragment_lca][hit_function]['identity'] += hit.get_identity()
                                    #~ else:
                                        #~ sample_scores[read_pe1.taxonomy][hit_function]['hit_count'] += 1.0 
                                        #~ sample_scores[read_pe1.taxonomy][hit_function]['identity'] += hit.get_identity()
                                    #~ del function_maxbitscores[hit_function]

                        #~ for hit in hits2:
                            #~ for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                                #~ if hit_function in function_maxbitscores and hit.get_bitscore() == function_maxbitscores[hit_function]:
                                    #~ sample_scores[read_pe2.taxonomy][hit_function]['hit_count'] += 1.0 
                                    #~ sample_scores[read_pe2.taxonomy][hit_function]['identity'] += hit.get_identity()
                                    #~ del function_maxbitscores[hit_function]


                else: #Only end1 is mapped
                    fragment_functions = set()
                    fragment_taxonomy = read_pe1.taxonomy
                    read_functions = read_pe1.get_functions()
                    fragment_functions.update(read_functions.keys())

#                    if 'RP-S9' in fragment_functions:
#                        print(read_id, '\tpe1\t', fragment_taxonomy)

                    # Count FPKM
                    for function in fragment_functions:
                        all_functions_list.add(function)
                        all_scores[sample_id][function] += norm_factor * read_functions[function]
                        sample_scores[fragment_taxonomy][function]['count'] += 1
                        sample_scores[fragment_taxonomy][function][metrics] += norm_factor * read_functions[function]
                    # Build FPKM-based functional taxonomic profile
                    function_maxbitscores = defaultdict(dict)
                    hits1 = [hit for hit in read_pe1.get_hit_list().get_hits()]
                    # Find max. bitscore for each function
                    for hit in hits1:
                        for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                            if hit_function in function_maxbitscores:
                                if hit.get_bitscore() > function_maxbitscores[hit_function]['bitscore']:
                                    function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                    function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                    function_maxbitscores[hit_function]['taxonomy'] = fragment_taxonomy
                            else:
                                function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                function_maxbitscores[hit_function]['taxonomy'] = fragment_taxonomy
                    # Count hits and sum identity
                    for function in fragment_functions:
                        if function in function_maxbitscores:
                            sample_scores[function_maxbitscores[function]['taxonomy']][function]['hit_count'] += 1.0 
                            sample_scores[function_maxbitscores[function]['taxonomy']][function]['identity'] += function_maxbitscores[hit_function]['identity']
                        else:
                            print('Function',function,'not found in hits of read', read_id)

                    #~ for hit in hits1:
                        #~ for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                            #~ if hit_function in function_maxbitscores and hit.get_bitscore() == function_maxbitscores[hit_function]:
                                #~ sample_scores[fragment_taxonomy][hit_function]['hit_count'] += 1.0 
                                #~ sample_scores[fragment_taxonomy][hit_function]['identity'] += hit.get_identity()
                                #~ del function_maxbitscores[hit_function]

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

#                if 'RP-S9' in fragment_functions:
#                    print(read_id, '\tpe2\t', fragment_taxonomy)

                # Count FPKM
                for function in fragment_functions:
                    all_functions_list.add(function)
                    all_scores[sample_id][function] += norm_factor * read_functions[function]
                    sample_scores[fragment_taxonomy][function]['count'] += 1
                    sample_scores[fragment_taxonomy][function][metrics] += norm_factor * read_functions[function]
                # Build FPKM-based functional taxonomic profile
                function_maxbitscores = defaultdict(dict)
                hits = [hit for hit in read_pe2.get_hit_list().get_hits()]
                # Find max. bitscore for each function
                for hit in hits:
                    for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                        if hit_function in function_maxbitscores:
                            if hit.get_bitscore() > function_maxbitscores[hit_function]['bitscore']:
                                function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                function_maxbitscores[hit_function]['taxonomy'] = fragment_taxonomy
                        else:
                            function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                            function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                            function_maxbitscores[hit_function]['taxonomy'] = fragment_taxonomy

                # Count hits and sum identity
                for function in fragment_functions:
                    if function in function_maxbitscores:
                        sample_scores[function_maxbitscores[function]['taxonomy']][function]['hit_count'] += 1.0 
                        sample_scores[function_maxbitscores[function]['taxonomy']][function]['identity'] += function_maxbitscores[hit_function]['identity']
                    else:
                        print('Function',function,'not found in hits of read', read_id)

                #~ for hit in hits:
                    #~ for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                        #~ if hit_function in function_maxbitscores and hit.get_bitscore() == function_maxbitscores[hit_function]:
                            #~ sample_scores[fragment_taxonomy][hit_function]['hit_count'] += 1.0 
                            #~ sample_scores[fragment_taxonomy][hit_function]['identity'] += hit.get_identity()
                            #~ del function_maxbitscores[hit_function]

            #~ for tax in sample_scores.keys():
                #~ for func in sample_scores[tax].keys():
                    #~ if not isinstance(sample_scores[tax][func]['identity'], float):
                        #~ print (tax, func)
            

            tax_profile = TaxonomyProfile()
            
            #outfile = os.path.join(self.parser.options.get_project_dir(sample), self.parser.options.get_output_subdir(sample), sample + '_' + end + '_' + 'functional_taxonomy_profile.xml')
#            outfile = os.path.join(sample + '_fpkm_naive_functional_taxonomy_profile.xml')
            outfile = sanitize_file_name(os.path.join(self.project.options.get_work_dir(), sample_id + '_' + metrics + '_naive_functional_taxonomy_profile.xml'))
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
            worksheet.conditional_format('C3:C1048560', {'type':     'text',
                                   'criteria': 'containing',
                                   'value':    'superkingdom',
                                   'format':   superkingdom_format})
            worksheet.conditional_format('C3:C1048560', {'type':     'text',
                                   'criteria': 'containing',
                                   'value':    'phylum',
                                   'format':   phylum_format})
            worksheet.conditional_format('C3:C1048560', {'type':     'text',
                                   'criteria': 'containing',
                                   'value':    'class',
                                   'format':   class_format})
            worksheet.conditional_format('C3:C1048560', {'type':     'text',
                                   'criteria': 'containing',
                                   'value':    'order',
                                   'format':   order_format})
            worksheet.conditional_format('C3:C1048560', {'type':     'text',
                                   'criteria': 'containing',
                                   'value':    'family',
                                   'format':   family_format})
            worksheet.conditional_format('C3:C1048560', {'type':     'text',
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

#    @unittest.skip("for faster testing")
    def test_7_build_fragment_function_lca_profile(self):
        self.project.import_reads_json(sample_id, self.project.ENDS)
        
        outfile = sanitize_file_name(self.project.options.get_name() + '_fpkm_functions.xlsx')
        metrics = 'fpkm'
        
        scores = get_function_scores(self.project,sample_id=sample_id,metrics='fpkm')
        
        generate_function_sample_xlsx(self.project, scores, metrics='fpkm', sample_id=sample_id)
        

        #~ tax_profile = TaxonomyProfile()
        #~ outfile = sanitize_file_name(os.path.join(self.project.options.get_work_dir(), sample_id + '_' + metrics + '_lca_functional_taxonomy_profile.xml'))
        #~ tax_profile.build_functional_taxonomy_profile(self.project.taxonomy_data, scores)
        #~ generate_lca_taxonomy_chart(tax_profile, sample=sample_id, outfile=outfile, score='fpkm')

        
        self.assertTrue(len(scores), 30)

    def test_8_build_fragment_function_taxonomy_lca_profile(self):
        #sample_id = 'sample1'
        self.project.import_reads_json(sample_id, ENDS)
        
        outfile = sanitize_file_name(self.project.options.project_name + '_fpkm_functions_taxonomy.xlsx')
        metric = 'efpkg'
        
        scores = get_function_taxonomy_scores(self.project,sample_id=sample_id,metric=metric)
        
        make_func_tax_sample_xlsx(self.project, scores, metric=metric, sample_id=sample_id)
        
        # Subsetting scores
        sample_scores = autovivify(3, float)
        for tax in scores.keys():
            for f in scores[tax].keys():
                if sample_id in scores[tax][f]:
                    for k,v in scores[tax][f][sample_id].items():
                        sample_scores[tax][f][k] = v



        tax_profile = TaxonomyProfile()
        outfile = sanitize_file_name(
            os.path.join(
                self.project.options.work_dir,
                sample_id + '_' + metric + '_lca_functional_taxonomy_profile.xml'
                )
            )
        tax_profile.make_function_taxonomy_profile(self.project.taxonomy_data, sample_scores)
        #generate_lca_taxonomy_chart(tax_profile, sample=sample_id, outfile=outfile, score=metrics)
        make_taxonomy_series_chart(
            tax_profile,
            sample_list=sorted(self.project.ref_data.functions_dict.keys()),
            outfile=outfile, krona_path=self.project.config.krona_path, metric=metric
            )
        
        self.assertTrue(sample_scores)

    def test_9_build_fragment_function_taxonomy_comparison_table (self):
        for sample_id in self.project.list_samples():
            self.project.import_reads_json(sample_id, ENDS)
        
        metrics = 'fpkg'
        scores = get_function_taxonomy_scores(self.project,sample_id=None,metrics=metrics)
        generate_sample_taxonomy_function_xlsx(self.project, scores, metrics=metrics, function_id=None)
        
        self.assertTrue(scores)

    def test_10_build_fragment_function_taxonomy_filtered_table (self):
        for sample_id in self.project.list_samples():
            self.project.import_reads_json(sample_id, self.project.ENDS)
        
        metrics = 'fpkg'
        taxon = 'phylum'
        scores = get_function_taxonomy_scores(self.project,sample_id=None,metrics=metrics)
        
        xlsxfile = 'test.xlsx'

        writer = pd.ExcelWriter(xlsxfile, engine='xlsxwriter')

        for function in sorted(self.project.ref_data.functions_dict.keys()):
        
            # Subsetting scores
            sample_scores = autovivify(3, float)
            for tax in scores.keys():
                if function in scores[tax].keys():
                    for s in self.project.list_samples():
                        if s in scores[tax][function]:
                            for k,v in scores[tax][function][s].items():
                                sample_scores[tax][s][k] = v
                        else:
                            sample_scores[tax][s][metrics] = 0.0
                    


            tax_profile = TaxonomyProfile()
            tax_profile.build_functional_taxonomy_profile(self.project.taxonomy_data, sample_scores)

            df = tax_profile.convert_taxonomic_profile_into_score_df(score=metrics)
            #print(df.index)
            #print(list(df))
            df_filtered = df[df[('','Rank')] == 'phylum']
            #print(df_filtered.head())
            
            df_filtered.to_excel(writer, sheet_name=function)
            workbook  = writer.book
            worksheet = writer.sheets[function]

            superkingdom_format = workbook.add_format({'bg_color': '#FF6666'})
            phylum_format = workbook.add_format({'bg_color': '#FF9900'})
            class_format = workbook.add_format({'bg_color': '#FFCC99'})
            order_format = workbook.add_format({'bg_color': '#FFFFCC'})
            family_format = workbook.add_format({'bg_color': '#99FFCC'})
        #    genus_format = workbook.add_format({'bg_color': '#99FFFF'})
            worksheet.conditional_format('C3:C1048560', {'type':     'text',
                                   'criteria': 'containing',
                                   'value':    'superkingdom',
                                   'format':   superkingdom_format})
            worksheet.conditional_format('C3:C1048560', {'type':     'text',
                                   'criteria': 'containing',
                                   'value':    'phylum',
                                   'format':   phylum_format})
            worksheet.conditional_format('C3:C1048560', {'type':     'text',
                                   'criteria': 'containing',
                                   'value':    'class',
                                   'format':   class_format})
            worksheet.conditional_format('C3:C1048560', {'type':     'text',
                                   'criteria': 'containing',
                                   'value':    'order',
                                   'format':   order_format})
            worksheet.conditional_format('C3:C1048560', {'type':     'text',
                                   'criteria': 'containing',
                                   'value':    'family',
                                   'format':   family_format})
            #~ worksheet.conditional_format('C3:C1048560', {'type':     'text',
                                   #~ 'criteria': 'containing',
                                   #~ 'value':    'genus',
                                   #~ 'format':   genus_format})
            worksheet.set_column(1, 1, 30)
            worksheet.set_column(2, 2, 15)

        writer.save()
        
        self.assertTrue(scores)

    def test_11_build_fragment_function_taxonomy_filtered_table (self):
        for sample_id in self.project.list_samples():
            self.project.import_reads_json(sample_id, self.project.ENDS)
        
        metrics = 'fpk'
        
        #rank = 'phylum'
        #ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus']
        scores = get_function_taxonomy_scores(self.project,sample_id=None,metrics=metrics)
        #for rank in ranks:
        #generate_sample_taxonomy_function_xlsx(self.project, scores, metrics=metrics, function_id='RP-S3', rank = None)
        generate_function_taxonomy_sample_xlsx(self.project, scores, metrics=metrics, sample_id=None, rank = None)
        
        self.assertTrue(scores)


    @unittest.skip("for faster testing")
    def test_2_get_hits_taxonomy_profile(self):
        self.parser.reads = import_annotated_reads(os.path.join(
            self.parser.options.get_project_dir(self.parser.sample.sample_id),
            self.parser.sample.sample_id + '_' + self.parser.end + '_' + self.parser.options.reads_json_name
        ))
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

    @unittest.skip("for faster testing")
    def test_3_make_taxonomy_profile(self):
        self.parser.reads = import_annotated_reads(os.path.join(
            self.parser.options.get_project_dir(self.parser.sample.sample_id),
            self.parser.sample.sample_id + '_' + self.parser.end + '_' + self.parser.options.reads_json_name
        ))
        scores = defaultdict(lambda : defaultdict(float))
        print(sample_id, end, 'read count', str(len(self.parser.reads)))
        for read_id, read in self.parser.reads.items():
            if read.status == STATUS_GOOD:
                for hit in read.hit_list.hits:
                    protein_taxid = self.parser.ref_data.lookup_protein_tax(hit.subject_id)
                    scores[protein_taxid]['count'] += 1.0
                    scores[protein_taxid]['identity'] += hit.identity
                if len(read.hit_list.hits) == 1:
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

    
    @unittest.skip("for faster testing")
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


    @unittest.skip("for faster testing")
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


    @unittest.skip("for faster testing")
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


    @unittest.skip("for faster testing")
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


    @unittest.skip("for faster testing")
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
