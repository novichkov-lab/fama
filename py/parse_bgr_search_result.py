#!/usr/bin/python
import sys,argparse
from lib.DiamondParser.DiamondParser import DiamondParser

def get_args():
    desc = '''This script parses DIAMOND tabular output of sequence reads
    search against reference protein library.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--config', help='Path to config.ini')
    parser.add_argument('--project', help='Path to project.ini')
    parser.add_argument('--sample', help='Sample ID')
    parser.add_argument('--end', help='paired end ID. Must be either pe1 or pe2')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return args

def main():
    args = get_args()

    parser = DiamondParser(args.config, args.project, args.sample, args.end)
    parser.parse_background_output()
    
    read_dict = parser.get_reads()
    print (len(read_dict))
    for read in read_dict.keys():
        print(read)
        print(read_dict[read].get_status())
        print(read_dict[read].get_functions())
        read_dict[read].show_hits()
        pass

if __name__=='__main__':
    main()

