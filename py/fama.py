#!/usr/bin/python
import os,sys,argparse
from subprocess import Popen, PIPE, CalledProcessError
from collections import Counter
from lib.DiamondParser.DiamondParser import DiamondParser

def get_args():
    desc = '''This program  parses DIAMOND tabular output of sequence reads
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
    if args.end != 'pe1' and args.end != 'pe2':
        print ('End parameter should be either pe1 or pe2')
        parser.print_help()
        sys.exit(1)
    return args

def run_ref_search(parser):
    print ('Starting DIAMOND')
    diamond_args = ['/usr/bin/diamond',
                    'blastx',
                    '--db',
                    parser.config.get_reference_diamond_db(parser.project.get_collection(parser.sample)),
                    '--query',
                    parser.project.get_fastq_path(parser.sample,parser.end),
                    '--out',
                    os.path.join(parser.project.get_project_dir(parser.sample), parser.sample + '_' + parser.end + '_'+ parser.project.get_ref_output_name()),
                    '--max-target-seqs',
                    '50',
                    '--evalue',
                    str(parser.config.get_evalue_cutoff(parser.project.get_collection(parser.sample))),
                    '--threads',
                    parser.config.get_threads(),
                    '--outfmt','6','qseqid','sseqid','pident','length','mismatch','slen','qstart','qend','sstart','send','evalue','bitscore'
                    ]

    with Popen(diamond_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)

    print ('DIAMOND finished')

def run_bgr_search(parser):
    print ('Starting DIAMOND')
    diamond_args = ['/usr/bin/diamond',
                    'blastx',
                    '--db',
                    parser.config.get_background_diamond_db(parser.project.get_collection(parser.sample)),
                    '--query',
                    os.path.join(parser.project.get_project_dir(parser.sample), parser.sample + '_' + parser.end + '_'+ parser.project.get_ref_hits_fastq_name()),
                    '--out',
                    os.path.join(parser.project.get_project_dir(parser.sample), parser.sample + '_' + parser.end + '_'+ parser.project.get_background_output_name()),
                    '--max-target-seqs',
                    '50',
                    '--evalue',
                    str(parser.config.get_evalue_cutoff(parser.project.get_collection(parser.sample))),
                    '--threads',
                    parser.config.get_threads(),
                    '--outfmt','6','qseqid','sseqid','pident','length','mismatch','slen','qstart','qend','sstart','send','evalue','bitscore'
                    ]

    with Popen(diamond_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)

    print ('DIAMOND finished')

def write_output(parser):
    outfile = os.path.join(parser.project.get_project_dir(parser.sample), parser.sample + '_' + parser.end + '_'+ parser.project. get_report_name())
    with open(outfile, 'w') as of:
        # Write general info
        of.write('\nRun info:\n')
        of.write('Sample ID:\t' + parser.sample + '\n')
        of.write('Paired end:\t' + parser.end + '\n')
        of.write('FASTQ file:\t' + parser.project.get_fastq_path(parser.sample, parser.end) + '\n')
        of.write('Total number of reads:\t' + str(parser.project.get_fastq1_readcount(parser.sample)) + '\n\n')

        # Write read statistics
        of.write('\nRead statistics:\n')
        read_stats = Counter()
        for read in parser.reads.keys():
            read_stats[parser.reads[read].get_status()] += 1
        for status in read_stats:
            of.write(status + '\t' + str(read_stats[status]) + '\n')

        # Write function scores
        of.write('\nFunction statistics:\n')
        func_stats = Counter()
        func_counts = Counter()
        for read in parser.reads.keys():
            functions = parser.reads[read].get_functions()
            for function in functions:
                func_stats[function] += functions[function]
                func_counts[function] += 1
        of.write('\nFunction\tRPKM score\tRead count\n')
        for function in sorted(func_stats.keys()):
            of.write(function + '\t' + parser.ref_data.lookup_function_name(function) 
                    + '\t' + str(func_stats[function]) + '\t' +
                    str(func_counts[function]) + '\n')
                    
        of.write('\nList of reads\n')
        # Write list of reads
        for read in sorted(parser.reads.keys()):
            of.write(read + '\t' + parser.reads[read].get_status() + '\n')
            of.write('\t' + str(parser.reads[read].get_functions()) + '\n')
            for hit in parser.reads[read].get_hit_list().get_hits():
                of.write('\t' + str(hit) + '\n')

            


def main():
    args = get_args()
    parser = DiamondParser(args.config, args.project, args.sample, args.end)

    # Search in reference database
    #run_ref_search(parser)
    
    # Process output of reference DB search
    #parser.parse_reference_output()
    
    #Import sequence data for selected sequence reads
    print ('Reading FASTQ file')
    #parser.import_fastq()
    
    print ('Exporting FASTQ ')
    #parser.export_hit_fastq()
    #parser.export_read_fastq()
    print ('Exporting hits')
    #parser.export_hit_list()
    
    # Search in background database
    #run_bgr_search(parser)

    # Process output of reference DB search
    parser.parse_background_output()
    
    # Generate output
    write_output(parser)

if __name__=='__main__':
    main()

