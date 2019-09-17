#!/usr/bin/python3
"""Runs Fama functional profiling pipeline"""
import sys
import argparse
from lib.utils.const import ENDS
from lib.functional_profiling_pipeline import fastq_pipeline
from lib.protein_functional_pipeline import protein_pipeline

def get_args():
    """Returns command-line arguments"""
    desc = '''This program runs functional profiling pipeline.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-c', dest='config', type=str, help='Path to config.ini')
    parser.add_argument('-p', dest='project', type=str, help='Path to project.ini')
    parser.add_argument('-s', dest='sample', type=str, default=None,
                        help='Sample ID (optional)')
    parser.add_argument('-e', dest='end', type=str, default=None,
                        help='paired end ID (expected values: pe1 or pe2). Optional')
    parser.add_argument('--prot', dest='prot', action='store_true',
                        help='Process protein sequences in FASTA format (default: False)')
    parser.set_defaults(prot=False)
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    if args.end is not None and args.end not in ENDS:
        print('End parameter must be either pe1 or pe2')
        parser.print_help()
        sys.exit(1)
    return args

def main():
    """Main function calling functional profiling module"""
    args = get_args()
    if args.prot:
        protein_pipeline(args)
    else:
        if args.sample is None:
            print('Running functional profiling for all samples in the project')
        else:
            print('Running functional profiling only for ', args.sample)
        fastq_pipeline(
            config_file=args.config, project_file=args.project,
            sample_identifier=args.sample, end_identifier=args.end)
    print('Done!')

if __name__ == '__main__':
    main()
