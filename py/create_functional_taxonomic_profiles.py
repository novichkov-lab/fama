#!/usr/bin/python
import os,sys,argparse
from collections import defaultdict
import xlsxwriter
import pandas as pd

from Fama.DiamondParser.hit_utils import cleanup_protein_id
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

def create_functional_taxonomic_profile(project,tax_data,sample,end,writer):
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
        for end in sorted(project.samples[sample]):
            create_functional_taxonomic_profile(project,tax_data,sample,end,writer)
    
    writer.save()

if __name__=='__main__':
    main()

