#!/usr/bin/python
import sys,argparse
from Fama.functional_profiling_pipeline import fastq_pipeline
from Fama.fasta_protein_pipeline import protein_pipeline

def get_args():
    desc = '''This program runs functional profiling pipeline.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-c', dest='config', type=str, help='Path to config.ini')
    parser.add_argument('-p', dest='project', type=str, help='Path to project.ini')
    parser.add_argument('-s', dest='sample', type=str, default=None,
                        help='Sample ID (optional)')
    parser.add_argument('-e', dest='end', type=str, default=None,
                        help='paired end ID (expected values: pe1 or pe2). Optional')
    parser.add_argument('--prot', dest='prot', default=False,
                        help='Process protein sequences in FASTA format (default: False)')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    if not args.end is None:
        if args.end != 'pe1' and args.end != 'pe2':
            print ('End parameter must be either pe1 or pe2')
            parser.print_help()
            sys.exit(1)
    return args

def main():
    args = get_args()
    if args.prot:
        protein_pipeline(args)
    else:
        fastq_pipeline(config_file=args.config, project_file=args.project, sample_identifier=args.sample, end_identifier=args.end)
    print('Done!')

if __name__=='__main__':
    main()

