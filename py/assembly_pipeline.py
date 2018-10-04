#!/usr/bin/python
import os,sys,argparse
from Fama.Project import Project
from Fama.GeneAssembler.GeneAssembler import GeneAssembler
from Fama.OutputUtil.JSONUtil import export_gene_assembly

def get_args():
    desc = '''This program runs gene co-assembly pipeline on a project 
    after read mapping.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--config', help='Path to config.ini')
    parser.add_argument('--project', help='Path to project.ini')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return args

def gene_assembly_pipeline(config_file, project_file):
    project = Project(config_file=config_file, project_file=project_file)
    assembler = GeneAssembler(project)
    assembler.assemble_contigs()
    assembler.map_genes()
    assembler.map_genes2uniprot()
    assembler.generate_output()
    export_gene_assembly(assembler.assembly,os.path.join(project.options.get_assembly_dir(), 'all_contigs_assembly.json'))

def main():
    args = get_args()
    gene_assembly_pipeline(config_file = args.config, project_file = args.project)

if __name__=='__main__':
    main()

