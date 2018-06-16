#!/usr/bin/python
import os,sys,argparse,gzip
from subprocess import Popen, PIPE, CalledProcessError
from collections import defaultdict
from Fama.DiamondParser.DiamondParser import DiamondParser
from Fama.OutputUtil.Report import generate_report
from Fama.OutputUtil.PdfReport import generate_protein_pdf_report
from Fama.OutputUtil.KronaXMLWriter import generate_xml
from Fama.OutputUtil.JSONUtil import export_annotated_reads
from Fama.OutputUtil.JSONUtil import import_annotated_reads

# This program no longer runs functional profiling for individual samples


def get_args():
    desc = '''This program  parses DIAMOND tabular output of sequence reads
    search against reference protein library.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--config', help='Path to config.ini')
    parser.add_argument('--project', help='Path to project.ini')
    parser.add_argument('--sample', help='Sample ID')
    parser.add_argument('--end', help='paired end ID. Must be pe1 for single-end runs')
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
                    'blastp',
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
                    'blastp',
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

def import_protein_fasta(parser):
    fasta_file = parser.project.get_fastq_path(parser.sample,parser.end)
    sequence  = []
    current_id = None
    fh = None
    if fasta_file.endswith('.gz'):
        fh = gzip.open(fasta_file, 'rb')
    else:
        fh = open(fasta_file, 'rb')
    if fh:
        for line in fh:
            line = line.decode('utf8').rstrip('\n\r')
            if line.startswith('>'):
                if current_id:
                    seq_id = current_id[1:].split(' ')[0]
                    parser.reads[seq_id].set_read_id_line(current_id)
                    parser.reads[seq_id].set_sequence(''.join(sequence))
                sequence = []
                seq_id = line[1:].split(' ')[0]
                if seq_id in parser.reads:
                    current_id = line
                else: 
                    current_id = None
                    seq_id = None
            else:
                if current_id:
                    sequence.append(line)
        if current_id:
            parser.reads[seq_id].set_read_id_line(current_id)
            parser.reads[seq_id].set_sequence(''.join(sequence))
        fh.close()

def calculate_protein_coverage(parser, infile):
    proteins = {} # this dict stores lists of tuples [start position, end position] for all proteins
    coverage_values = defaultdict(dict) # this dict stores coverage values for all positions of interest
    
    for read in parser.reads:
        protein_id = parser.reads[read].get_read_id_line()
        protein_id_tokens = protein_id.split(' # ')
        contig_id = protein_id_tokens[0]
        contig_id = '_'.join(contig_id.split('_')[:-1])
        contig_id = contig_id[1:]
        if contig_id not in proteins:
            proteins[contig_id] = {}
        proteins[contig_id][read] = {}
        proteins[contig_id][read]['start'] = int(protein_id_tokens[1])
        proteins[contig_id][read]['end'] = int(protein_id_tokens[2])
   
    print ('Reading coverage file...')
    fh = None
    if infile.endswith('.gz'):
        fh = gzip.open(infile, 'rb')
    else:
        fh = open(infile, 'rb')
    if fh:
        for line in fh:
            line = line.decode('utf8').rstrip('\n\r')
            contig_cov, position, coverage = line.split('\t')
            if contig_cov in proteins.keys():
                for protein_id in proteins[contig_cov].keys():
                    start = proteins[contig_cov][protein_id]['start']
                    end = proteins[contig_cov][protein_id]['end']
                    if (int(position) > start -1) and int(position) <= end:
                        coverage_values[contig_cov][int(position)] = int(coverage)
        fh.close()
    print ('Calculating average coverage ...')
    
    for read in parser.reads:
        # first, calculate average count for protein
        protein_id = parser.reads[read].get_read_id_line()
        protein_id_tokens = protein_id[1:].split(' ')
        contig_id = '_'.join(protein_id_tokens[0].split('_')[:-1])
        if contig_id in coverage_values:
            i = proteins[contig_id][read]['start']
            cov_arr = []
            if i:
                while i <= proteins[contig_id][read]['end']:
                    if i in coverage_values[contig_id]:
                        cov_arr.append(coverage_values[contig_id][i])
                    i += 1
        coverage_avg = sum(cov_arr)/len(cov_arr)
        
        for function in parser.reads[read].get_functions():
            if coverage_avg:
                parser.reads[read].functions[function] = coverage_avg
            else:
                parser.reads[read].functions[function] = 0.0 # just in case we have no coverage data
    

    
def functional_profiling_pipeline(config_file, project_file, sample, end, coverage_file):
    parser = DiamondParser(config_file=config_file, project_file=project_file, sample=sample, end=end)
    
    if not os.path.isdir(parser.project.get_project_dir(parser.sample)):
        os.mkdir(parser.project.get_project_dir(parser.sample))
    if not os.path.isdir(os.path.join(parser.project.get_project_dir(parser.sample),parser.project.get_output_subdir(parser.sample))):
        os.mkdir(os.path.join(parser.project.get_project_dir(parser.sample),parser.project.get_output_subdir(parser.sample)))

    # Search in reference database
    #run_ref_search(parser)
    
    # Process output of reference DB search
    #parser.parse_reference_output()
    
    ##Import sequence data for selected sequence reads
    #print ('Reading FASTQ file')
    #import_protein_fasta(parser)
    
    #print ('Exporting FASTQ ')
    #parser.export_hit_fasta()
    #print ('Exporting hits')
    #parser.export_hit_list()
    
    ## Search in background database
    ##run_bgr_search(parser)

    ## Process output of reference DB search
    #parser.parse_background_output()
    
    if not parser.reads:
        print ('Import JSON file')
        parser.reads = import_annotated_reads(os.path.join(parser.project.get_project_dir(sample), sample + '_' + end + '_' + parser.project.get_reads_json_name()))
        
    #coverage_file = '/mnt/data2/FEBA/proteins/D16-4706_coverage.tsv.gz'
    print ('Calculating coverage values')
    #calculate_protein_coverage(parser, coverage_file)
    
    print('Exporting JSON')
    parser.export_read_fasta()
    export_annotated_reads(parser)
    
    # Generate output
    print('Generating reports')
    generate_report(parser)
    generate_protein_pdf_report(parser)
    generate_xml(parser)

def main():
    args = get_args()
    functional_profiling_pipeline(config_file=args.config, project_file=args.project, sample=args.sample, end=args.end)
    print('Done!')

if __name__=='__main__':
    main()

