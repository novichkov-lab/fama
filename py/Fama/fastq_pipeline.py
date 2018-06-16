#!/usr/bin/python
import os,sys,argparse
from subprocess import Popen, PIPE, CalledProcessError
from collections import Counter
from Fama.DiamondParser.DiamondParser import DiamondParser
from Fama.OutputUtil.Report import generate_report
from Fama.OutputUtil.PdfReport import generate_pdf_report
from Fama.OutputUtil.KronaXMLWriter import generate_xml
from Fama.OutputUtil.JSONUtil import export_annotated_reads

# This program runs functional profiling for individual FASTQ file


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
                    str(parser.config.get_background_db_size(parser.project.get_collection(parser.sample)) 
                        * parser.config.get_evalue_cutoff(parser.project.get_collection(parser.sample))
                        / parser.config.get_reference_db_size(parser.project.get_collection(parser.sample))),
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


def functional_profiling_pipeline(config_file, project_file, sample, end):
    parser = DiamondParser(config_file=config_file, project_file=project_file, sample=sample, end=end)
    
    if not os.path.isdir(parser.project.get_project_dir(parser.sample)):
        os.mkdir(parser.project.get_project_dir(parser.sample))
    if not os.path.isdir(os.path.join(parser.project.get_project_dir(parser.sample),parser.project.get_output_subdir(parser.sample))):
        os.mkdir(os.path.join(parser.project.get_project_dir(parser.sample),parser.project.get_output_subdir(parser.sample)))

    # Search in reference database
    run_ref_search(parser)
    
    # Process output of reference DB search
    parser.parse_reference_output()
    
    #Import sequence data for selected sequence reads
    print ('Reading FASTQ file')
    parser.import_fastq()
    
    print ('Exporting FASTQ ')
    parser.export_hit_fastq()
    print ('Exporting hits')
    parser.export_hit_list()
    
    # Search in background database
    run_bgr_search(parser)

    # Process output of reference DB search
    parser.parse_background_output()

    parser.export_read_fastq()
    parser.export_paired_end_reads_fastq()
    export_annotated_reads(parser)
    
    # Generate output
    generate_report(parser)
    generate_pdf_report(parser)
    generate_xml(parser)


def main():
    
    print('This program is not intended to run directly')

if __name__=='__main__':
    main()

