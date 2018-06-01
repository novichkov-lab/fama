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
    if args.end != 'pe1' or args.end != 'pe2':
        parser.print_help()
        sys.exit(1)
    return args

def main():
    args = get_args()

    #parser = DiamondParser(config_file, collection)
    parser = DiamondParser(args.config, args.project, args.sample, args.end)
    #parser.parse_tabular_output(infile)
    parser.parse_reference_output()
    
    #Import sequence data for selected sequence reads
    print ('Reading FASTQ file')
    parser.import_fastq()
    
    print ('Exporting FASTQ ')
    parser.export_hit_fastq()
    parser.export_read_fastq()
    print ('Exporting hits')
    parser.export_hit_list()
    
    read_dict = parser.get_reads()
    print (len(read_dict))
    for read in read_dict.keys():
        read_dict[read].show_hits()

if __name__=='__main__':
    main()

