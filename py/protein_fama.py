#!/usr/bin/python
import sys,argparse
from Fama.fasta_protein_pipeline import functional_profiling_pipeline

def get_args():
    desc = '''This program  parses DIAMOND tabular output of sequence reads
    search against reference protein library.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--config', help='Path to config.ini')
    parser.add_argument('--project', help='Path to project.ini')
    parser.add_argument('--coverage', help='Path to coverage file')
    parser.add_argument('--sample', help='Sample ID')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return args

def main():
    args = get_args()
    functional_profiling_pipeline(config_file=args.config, project_file=args.project, sample=args.sample)
    print('Done!')

if __name__=='__main__':
    main()

