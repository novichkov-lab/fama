#!/usr/bin/python
import os,sys,argparse
from collections import defaultdict
import xlsxwriter
import pandas as pd

from Fama.utils import cleanup_protein_id
from Fama.Project import Project
from Fama.ReferenceLibrary.TaxonomyData import TaxonomyData
from Fama.TaxonomyProfile import TaxonomyProfile
from Fama.OutputUtil.KronaXMLWriter import generate_functional_taxonomy_chart

def get_args():
    desc = '''This program generates Krona chart for taxonomy profile of all functions.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--config', help='Path to config.ini')
    parser.add_argument('--project', help='Path to project.ini')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    if not args.config or not args.project:
        parser.print_help()
        sys.exit(1)
    return args

def create_functional_taxonomic_profile_singlefile(project,tax_data,sample,end,writer):
    reads = project.samples[sample][end]
    scores = defaultdict(lambda : defaultdict(dict))
    function_list = set()
    for read in reads:
        if reads[read].get_status() == 'function,besthit' or reads[read].get_status() == 'function':
            read_functions = reads[read].get_functions()
            hits = reads[read].get_hit_list().get_hits()
            protein_taxids = {}
            for function in read_functions:
                for hit in hits:
                    protein_taxid = project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                    for hit_function in hit.get_functions():
                        if hit_function == function:
                            function_list.add(function)
                            if function in scores[protein_taxid]:
                                scores[protein_taxid][function]['count'] += 1.0
                                scores[protein_taxid][function]['identity'] += hit.get_identity()
                                scores[protein_taxid][function]['rpkm'] += read_functions[function]
                            else:
                                scores[protein_taxid][function]['count'] = 1.0
                                scores[protein_taxid][function]['identity'] = hit.get_identity()
                                scores[protein_taxid][function]['rpkm'] = read_functions[function]
    tax_profile = TaxonomyProfile()
    outfile = os.path.join(project.options.get_project_dir(sample), project.options.get_output_subdir(sample), sample + '_' + end + '_' + 'functional_taxonomy_profile.xml')
    tax_profile.build_functional_taxonomy_profile(tax_data, scores)
    generate_functional_taxonomy_chart(tax_profile, sorted(function_list), outfile)
    df = tax_profile.convert_function_taxonomic_profile_into_df()

    #print(df)
    df.to_excel(writer, sheet_name=sample+'_'+end)
    workbook  = writer.book
    worksheet = writer.sheets[sample + '_' + end]
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

def create_functional_taxonomic_profile(project,tax_data,sample,writer):
    scores = defaultdict(lambda : defaultdict(dict))
    function_list = set()
    for end in sorted(project.samples[sample]):
        scaling_factor = 1.0
        if end == 'pe1':
            scaling_factor = project.options.get_fastq1_readcount(sample)/(project.options.get_fastq1_readcount(sample) + project.options.get_fastq2_readcount(sample))
        elif end == 'pe2':
            scaling_factor = project.options.get_fastq2_readcount(sample)/(project.options.get_fastq1_readcount(sample) + project.options.get_fastq2_readcount(sample))
        else:
            raise Exception('Unknown end identifier')


        reads = project.samples[sample][end]
        for read in reads:
            if reads[read].get_status() == 'function,besthit' or reads[read].get_status() == 'function':
                read_functions = reads[read].get_functions()
                hits = reads[read].get_hit_list().get_hits()
                # If we have only one hit, all RPKM scores of the read would be assigned to the tax id of the hit
                # If we have more than one hit, RPKM scores would be equally divided between tax ids of the hits, with regard to functional assignments of the hits
                function_taxids = defaultdict(lambda : defaultdict(float))
                # First, collect all taxonomy IDs for each function assigned to the hit
                tax_functions_count = 0.0
                for read_function in read_functions:
                    for hit in hits:
                        hit_taxid = project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
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

    tax_profile = TaxonomyProfile()
    outfile = os.path.join(project.options.get_project_dir(sample), project.options.get_output_subdir(sample), sample + '_functional_taxonomy_profile.xml')
    tax_profile.build_functional_taxonomy_profile(tax_data, scores)
    generate_functional_taxonomy_chart(tax_profile, sorted(function_list), outfile)
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


def main():
    args = get_args()
    
    project = Project(config_file=args.config, project_file=args.project)
    tax_data = TaxonomyData(project.config)
    tax_data.load_taxdata(project.config)
    project.load_functional_profile()
    
    xlsxfile = os.path.join(project.options.get_work_dir(), project.options.get_name() + '_functions_taxonomy.xlsx')
    xlsxfile = xlsxfile.replace(' ', '_')
    xlsxfile = xlsxfile.replace("'", "")
    xlsxfile = xlsxfile.replace('"', '')
    writer = pd.ExcelWriter(xlsxfile, engine='xlsxwriter')

    for sample in sorted(project.samples):
        create_functional_taxonomic_profile(project,tax_data,sample,writer)
    
    writer.save()

if __name__=='__main__':
    main()

