#!/usr/bin/python3
import sys,argparse
from Fama.trimming_pipeline import trimming_pipeline

def get_args():
    desc = '''This program runs functional profiling pipeline.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-c', dest='config', type=str, help='Path to config.ini')
    parser.add_argument('-i', dest='input_file', type=str, help='Path to list of sequence files (tab-separated: sample ID, path to FASTA/first FASTQ, path to paired-end FASTQ)')
    parser.add_argument('-r', dest='collection', type=str, help='Reference collection ID')
    parser.add_argument('-p', dest='project_name', type=str, default='Fama project', help='Project name (optional)')
    parser.add_argument('--prot', dest='prot',  action='store_true',
                        help='Input files contain protein sequences (optional, default: False)')
    parser.set_defaults(prot=False)
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return args

def main():
    args = get_args()
    trimming_pipeline(config_file=args.config, sequence_list_file=args.input_file, project_name=args.project_name, collection=args.collection, is_protein=args.prot)
    print('Done!')

if __name__=='__main__':
    main()

