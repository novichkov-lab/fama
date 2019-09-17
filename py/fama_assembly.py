#!/usr/bin/python3
"""Script starting Fama assembly pipeline"""
import sys
import argparse
from lib.assembly_pipeline import assembly_pipeline


def get_args():
    """Returns command-line arguments"""
    desc = '''This program runs gene assembly pipeline.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-c', dest='config', type=str, help='Path to config.ini')
    parser.add_argument('-p', dest='project', type=str, help='Path to project.ini')
    parser.add_argument('-a', dest='assembler', type=str, default='metaspades',
                        help='contig assembler (expected values: metaspades or megahit).\
                         Optional, default:metaspades')
    parser.add_argument('--coassembly', dest='coassembly', action='store_true',
                        help='Coassemble reads for all functions (default: False)')
    parser.set_defaults(coassembly=False)
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return args


def main():
    """Main function calling assembly module"""
    args = get_args()
    assembly_pipeline(args)

if __name__ == '__main__':
    main()
