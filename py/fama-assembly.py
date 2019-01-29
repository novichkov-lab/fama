#!/usr/bin/python3
import sys,argparse
from Fama.assembly_pipeline import assembly_pipeline

def get_args():
    desc = '''This program runs functional profiling pipeline.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-c', dest='config', type=str, help='Path to config.ini')
    parser.add_argument('-p', dest='project', type=str, help='Path to project.ini')
    parser.add_argument('-s', dest='sample', type=str, default=None,
                        help='Sample ID (optional)')
    parser.add_argument('-e', dest='end', type=str, default=None,
                        help='paired end ID (expected values: pe1 or pe2). Optional')
    parser.add_argument('--coassembly', dest='coassembly', action='store_true',
                        help='Coassemble reads for all functions (default: False)')
    parser.set_defaults(coassembly=False)
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
    assembly_pipeline(args)

if __name__=='__main__':
    main()

