#!/usr/bin/python
import sys,argparse
from lib.DiamondParser.DiamondParser import DiamondParser

def get_args():
    desc = '''This script parses DIAMOND tabular output of sequence reads
    search against reference protein library.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--config', help='Path to config.ini')
    parser.add_argument('--infile', help='Input file (tabular DIAMOND output)')
    parser.add_argument('--fastq', help='FASTQ file')
    parser.add_argument('--col', help='Reference collection name')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return args

def main():
    args = get_args()

    config_file = '/mnt/data2/SEED/fama/py/config.ini'
    collection = 'nitrogen_v7'
    infile = '/mnt/data2/SEED/fama/t/data/test_diamond_output.txt'
    outdir = '/mnt/data2/SEED/fama/t/data/'

    #parser = DiamondParser(config_file, collection)
    parser = DiamondParser(args.config, args.col)
    #parser.parse_tabular_output(infile)
    parser.parse_tabular_output(args.infile)
    
    #Import sequence data for selected sequence reads
    print ('Reading FASTQ file')
    parser.import_fastq(args.fastq)
    
    print ('Exporting FASTQ ')
    parser.export_hit_fastq(outdir)
    print ('Exporting hits')
    parser.export_hit_list(outdir)
    
    read_dict = parser.get_reads()
    print (len(read_dict))
    for read in read_dict.keys():
        read_dict[read].show_hits()

if __name__=='__main__':
    main()

