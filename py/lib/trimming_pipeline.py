#!/usr/bin/python3
import os,gzip
from subprocess import Popen, PIPE, CalledProcessError
from collections import defaultdict

from lib.utils.const import ENDS
from lib.project.program_config import ProgramConfig

def trimming_pipeline(config_file, sequence_list_file, project_name, collection, is_protein, notrim):
    
    # read config, check collection
    config = ProgramConfig(config_file)
    if collection not in config.collections:
        raise ValueError('Collection ' + collection + ' not found in Fama config file ' + config_file)
    
    # read sequence list, check file paths
    input_files = defaultdict(dict)
    working_directory = os.path.dirname(os.path.realpath(sequence_list_file))
    print(working_directory)
    with open(sequence_list_file, 'r') as f:
        for line in f:
            sample_id, fwd_sequence_file, rev_sequence_file = line.rstrip('\n\r').split('\t')
            if os.path.exists(fwd_sequence_file):
                input_files[sample_id][ENDS[0]] = fwd_sequence_file
            elif os.path.exists(os.path.join(working_directory, fwd_sequence_file)):
                input_files[sample_id][ENDS[0]] = os.path.join(working_directory, fwd_sequence_file)
            else:
                raise OSError('Sequence file for sample ' +  sample_id + ' not found: ' + fwd_sequence_file)
            if rev_sequence_file == '':
                input_files[sample_id][ENDS[1]] = ''
            elif os.path.exists(rev_sequence_file):
                input_files[sample_id][ENDS[1]] = rev_sequence_file
            elif os.path.exists(os.path.join(working_directory, rev_sequence_file)):
                input_files[sample_id][ENDS[1]] = os.path.join(working_directory, rev_sequence_file)
            else:
                raise OSError('Sequence file for sample ' +  sample_id + ' not found: ' + rev_sequence_file)
    
    sequence_files = defaultdict(dict)
    # run Trimmomatic, if is_protein and notrim are False
    for sample_id in input_files:
        if is_protein or notrim:
            sequence_files[sample_id][ENDS[0]] = input_files[sample_id][ENDS[0]]
            sequence_files[sample_id][ENDS[1]] = input_files[sample_id][ENDS[1]]
        else:
            if is_fastq(input_files[sample_id][ENDS[0]]):
                fwd_outfile, rev_outfile = run_trimmomatic(input_files[sample_id][ENDS[0]], input_files[sample_id][ENDS[1]], sample_id, working_directory, config.threads)
                sequence_files[sample_id][ENDS[0]] = fwd_outfile
                sequence_files[sample_id][ENDS[1]] = rev_outfile
            else:
                sequence_files[sample_id][ENDS[0]] = input_files[sample_id][ENDS[0]]
                sequence_files[sample_id][ENDS[1]] = input_files[sample_id][ENDS[1]]
                

    # save project.ini
    
    project_ini_path = os.path.join(working_directory, 'project_' + collection + '.ini')

    working_directory = os.path.join(working_directory,collection)

    with open (project_ini_path, 'w') as of:
        # Boilerplate
        of.write('[DEFAULT]\n')
        of.write('project_name = \'' + project_name + '\'\n')
        of.write('collection = ' + collection + '\n')
        of.write('ref_output_name = ref_tabular_output.txt\n')
        of.write('background_output_name = bgr_tabular_output.txt\n')
        of.write('ref_hits_list_name = ref_hits.txt\n')
        of.write('ref_hits_fastq_name = ref_hits.fq\n')
        of.write('reads_fastq_name = reads.fq\n')
        of.write('pe_reads_fastq_name = reads_pe.fq\n')
        of.write('output_subdir = out\n')
        of.write('report_name = report.txt\n')
        of.write('xml_name = krona.xml\n')
        of.write('html_name = functional_profile.html\n')
        of.write('reads_json_name = reads.json\n')
        of.write('assembly_subdir = assembly\n')
        of.write('work_dir = ' + working_directory + '\n')
        
        #Sample sections
        replicate = 0
        for sample_id in sequence_files:
            of.write('\n[' + sample_id + ']\n')
            of.write('sample_id = ' + sample_id + '\n')
            of.write('fastq_pe1 = ' + sequence_files[sample_id][ENDS[0]] + '\n')
            of.write('fastq_pe2 = ' + sequence_files[sample_id][ENDS[1]] + '\n')
            of.write('sample_dir = ' + os.path.join(working_directory, sample_id) + '\n')
            of.write('replicate = ' + str(replicate) + '\n')
            replicate += 1


def run_trimmomatic(file1,file2,sample_id,working_directory,threads):
    print ('Starting Trimmomatic')
    if file2 == '':
        outfile = os.path.join(working_directory, sample_id + '_SE.fastq.gz'),
        trimmomatic_args = ['TrimmomaticSE',
                            '-threads', threads,
                            '-phred33',
                            file1,
                            outfile,
                            'ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-PE.fa:2:30:10',
                            'LEADING:3',
                            'TRAILING:3',
                            'SLIDINGWINDOW:4:14',
                            'MINLEN:50'
                            ]

        with Popen(trimmomatic_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
            for line in p.stdout:
                print(line, end='')
        if p.returncode != 0:
            raise CalledProcessError(p.returncode, p.args)
        print ('Trimmomatic finished')
        return outfile, ''
        
    else:
        trimmomatic_args = ['TrimmomaticPE',
                            '-threads', threads,
                            '-phred33',
                            file1,
                            file2,
                            '-baseout', os.path.join(working_directory, sample_id + '.fastq.gz'),
                            'ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-PE.fa:2:30:10',
                            'LEADING:3',
                            'TRAILING:3',
                            'SLIDINGWINDOW:4:14',
                            'MINLEN:50'
                            ]

        with Popen(trimmomatic_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
            for line in p.stdout:
                print(line, end='')
        if p.returncode != 0:
            raise CalledProcessError(p.returncode, p.args)
        print ('Trimmomatic finished')
        outfile1 = os.path.join(working_directory, sample_id + '_1P.fastq.gz')
        outfile2 = os.path.join(working_directory, sample_id + '_2P.fastq.gz')
        return outfile1, outfile2
    
def is_fastq(infile):
    first_symbol = ''
    if infile.endswith('.gz'):
        fh = gzip.open(infile, 'rb')
        first_symbol = fh.readline().decode('utf8')[0]
        fh.close()
    else:
        with open(infile,'r') as f:
            first_symbol = f.readline()[0]

    if first_symbol == '@':
        return True
    else:
        return False
    
