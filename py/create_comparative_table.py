#!/usr/bin/python
import os,sys,argparse
from subprocess import Popen, PIPE, CalledProcessError
from collections import Counter
from Fama.Project import Project
from Fama.OutputUtil.Report import generate_functions_scores_table
from Fama.OutputUtil.Report import generate_functions_scores_list
from Fama.OutputUtil.JSONUtil import export_annotated_reads

def get_args():
    desc = '''This program generates comparative table for a project.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--config', help='Path to config.ini')
    parser.add_argument('--project', help='Path to project.ini')
    parser.add_argument('--outfile', help='Output file name')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    if not args.config or not args.project or not args.outfile:
        parser.print_help()
        sys.exit(1)
    return args


def main():
    args = get_args()
    
    project = Project(config_file=args.config, project_file=args.project)
    #project.load_functional_profile()
    with open(args.outfile, 'w') as of:
        of.write(generate_functions_scores_table(project))
        of.close()

if __name__=='__main__':
    main()

