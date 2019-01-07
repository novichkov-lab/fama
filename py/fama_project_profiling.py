#!/usr/bin/python3
import os,sys,argparse
import pandas as pd
from Fama.Project import Project
from Fama.ReferenceLibrary.TaxonomyData import TaxonomyData
from Fama.OutputUtil.Report import create_functions_xlsx
from Fama.OutputUtil.Report import create_functions_markdown_document
from create_functional_taxonomic_profiles import create_functional_taxonomic_profile

def get_args():
    desc = '''This program runs functional profiling pipeline.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--config', help='Path to config.ini')
    parser.add_argument('--project', help='Path to project.ini')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return args

def main():
    args = get_args()
    project = Project(config_file=args.config, project_file=args.project)
    project.generate_functional_profile()

    if not os.path.isdir(project.options.get_work_dir()):
        os.mkdir(project.options.get_work_dir())
        
    #create_functions_xlsx(project)
    #create_functions_markdown_document(project)

    tax_data = project.taxonomy_data

    xlsxfile = os.path.join(project.options.get_work_dir(), project.options.get_name() + '_functions_taxonomy.xlsx')
    xlsxfile = xlsxfile.replace(' ', '_')
    xlsxfile = xlsxfile.replace("'", "")
    xlsxfile = xlsxfile.replace('"', '')
    writer = pd.ExcelWriter(xlsxfile, engine='xlsxwriter')

    ENDS = ['pe1','pe2']
    for sample in sorted(project.list_samples()):
        for end in ENDS:
            project.load_annotated_reads(sample, end)
        create_functional_taxonomic_profile(project,tax_data,sample,writer)
        for end in ENDS:
            project.samples[sample][end] = None

    writer.save()

    print('Done!')

if __name__=='__main__':
    main()

