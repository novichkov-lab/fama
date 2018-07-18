#!/usr/bin/python
import os,sys,argparse
from subprocess import Popen, PIPE, CalledProcessError
from collections import Counter
from Fama.Project import Project
from Fama.OutputUtil.Report import create_functions_xlsx
from Fama.OutputUtil.Report import create_functions_markdown_document


def get_args():
    desc = '''This program generates table of all functions for all samples in a project.'''
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


def main():
    args = get_args()
    # load data
    project = Project(config_file=args.config, project_file=args.project)
    project.load_functional_profile()
    # generate output
    create_functions_xlsx(project)
    create_functions_markdown_document(project)

if __name__=='__main__':
    main()

