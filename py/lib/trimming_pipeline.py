#!/usr/bin/python3
"""Runs Fama trimming pipeline"""
import os
import gzip
from collections import defaultdict

from lib.utils.const import ENDS
from lib.utils.utils import run_external_program
from lib.project.program_config import ProgramConfig


def trimming_pipeline(config_file, sequence_list_file, project_name, collection,
                      is_protein, notrim):
    """Defines steps of the pipeline

    Args:
        config_file (str): path to program configuration ini file
        sequence_list_file (str): path to file with list of input files
        project_name (str): project name
        collection (str): Fama collection identifier
        is_protein (bool): True if input sequences are proteins, False if nucleotides
        notrim (bool): True if no trimming needed
    """
    # read config, check collection
    config = ProgramConfig(config_file)
    if collection not in config.collections:
        raise ValueError(
            'Collection ' + collection + ' not found in Fama config file ' + config_file
            )

    # read sequence list, check file paths
    input_files = defaultdict(dict)
    working_directory = os.path.dirname(os.path.realpath(sequence_list_file))
    print(working_directory)
    with open(sequence_list_file, 'r') as infile:
        for line in infile:
            sample_id, fwd_sequence_file, rev_sequence_file = line.rstrip('\n\r').split('\t')
            if os.path.exists(fwd_sequence_file):
                input_files[sample_id][ENDS[0]] = fwd_sequence_file
            elif os.path.exists(os.path.join(working_directory, fwd_sequence_file)):
                input_files[sample_id][ENDS[0]] = os.path.join(
                    working_directory, fwd_sequence_file
                    )
            else:
                raise OSError(
                    'Sequence file for sample ' + sample_id + ' not found: ' + fwd_sequence_file
                    )
            if rev_sequence_file == '':
                input_files[sample_id][ENDS[1]] = ''
            elif os.path.exists(rev_sequence_file):
                input_files[sample_id][ENDS[1]] = rev_sequence_file
            elif os.path.exists(os.path.join(working_directory, rev_sequence_file)):
                input_files[sample_id][ENDS[1]] = os.path.join(
                    working_directory, rev_sequence_file
                    )
            else:
                raise OSError(
                    'Sequence file for sample ' + sample_id + ' not found: ' + rev_sequence_file
                    )

    sequence_files = defaultdict(dict)
    # run Trimmomatic, if is_protein and notrim are False
    for sample_id in input_files:
        if is_protein or notrim:
            sequence_files[sample_id][ENDS[0]] = input_files[sample_id][ENDS[0]]
            sequence_files[sample_id][ENDS[1]] = input_files[sample_id][ENDS[1]]
        else:
            if is_fastq(input_files[sample_id][ENDS[0]]):
                fwd_outfile, rev_outfile = run_trimmomatic(
                    input_files[sample_id][ENDS[0]],
                    input_files[sample_id][ENDS[1]],
                    sample_id, working_directory, config.threads
                    )
                sequence_files[sample_id][ENDS[0]] = fwd_outfile
                sequence_files[sample_id][ENDS[1]] = rev_outfile
            else:
                sequence_files[sample_id][ENDS[0]] = input_files[sample_id][ENDS[0]]
                sequence_files[sample_id][ENDS[1]] = input_files[sample_id][ENDS[1]]

    # save project.ini
    project_ini_path = os.path.join(working_directory, 'project_' + collection + '.ini')

    working_directory = os.path.join(working_directory, collection)

    with open(project_ini_path, 'w') as outfile:
        # Boilerplate
        outfile.write('[DEFAULT]\n')
        outfile.write('project_name = \'' + project_name + '\'\n')
        outfile.write('collection = ' + collection + '\n')
        outfile.write('ref_output_name = ref_tabular_output.txt\n')
        outfile.write('background_output_name = bgr_tabular_output.txt\n')
        outfile.write('ref_hits_list_name = ref_hits.txt\n')
        outfile.write('ref_hits_fastq_name = ref_hits.fq\n')
        outfile.write('reads_fastq_name = reads.fq\n')
        outfile.write('pe_reads_fastq_name = reads_pe.fq\n')
        outfile.write('output_subdir = out\n')
        outfile.write('report_name = report.txt\n')
        outfile.write('xml_name = krona.xml\n')
        outfile.write('html_name = functional_profile.html\n')
        outfile.write('reads_json_name = reads.json\n')
        outfile.write('assembly_subdir = assembly\n')
        outfile.write('work_dir = ' + working_directory + '\n')

        # Sample sections
        replicate = 0
        for sample_id in sequence_files:
            outfile.write('\n[' + sample_id + ']\n')
            outfile.write('sample_id = ' + sample_id + '\n')
            outfile.write('fastq_pe1 = ' + sequence_files[sample_id][ENDS[0]] + '\n')
            outfile.write('fastq_pe2 = ' + sequence_files[sample_id][ENDS[1]] + '\n')
            outfile.write('sample_dir = ' + os.path.join(working_directory, sample_id) + '\n')
            outfile.write('replicate = ' + str(replicate) + '\n')
            replicate += 1


def run_trimmomatic(file1, file2, sample_id, working_directory, threads):
    """Runs trimmomatic for one or two files

   Args:
        file1 (str): path to input file 1 (FASTQ)
        file2 (str): path to input file 2 (FASTQ paired-end) or None
        sample_id (str): sample identifier
        working_directory: directory where Trimmomatic will write output files
        threads (str): number of threads

    Returns:
        outfile1 (str): path to FASTQ file with trimmed paired-end1 reads
        outfile2 (str): path to FASTQ file with trimmed paired-end2 reads
    """
    print('Starting Trimmomatic')
    outfile1 = ''
    outfile2 = ''
    if file2 == '':
        outfile1 = os.path.join(working_directory, sample_id + '_SE.fastq.gz')
        trimmomatic_args = ['TrimmomaticSE',
                            '-threads', threads,
                            '-phred33',
                            file1,
                            outfile1,
                            'ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-PE.fa:2:30:10',
                            'LEADING:3',
                            'TRAILING:3',
                            'SLIDINGWINDOW:4:14',
                            'MINLEN:50']
        run_external_program(trimmomatic_args)
        print('Trimmomatic finished')
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
                            'MINLEN:50']
        run_external_program(trimmomatic_args)
        print('Trimmomatic finished')
        outfile1 = os.path.join(working_directory, sample_id + '_1P.fastq.gz')
        outfile2 = os.path.join(working_directory, sample_id + '_2P.fastq.gz')
    return outfile1, outfile2


def is_fastq(infile):
    """Checks input file format

    Args:
        infile (str): path to input sequence file

    Returns:
        result(bool): True if infile is FASTQ file, otherwise returns False
    """
    first_symbol = ''
    if infile.endswith('.gz'):
        file_handle = gzip.open(infile, 'rb')
        first_symbol = file_handle.readline().decode('utf8')[0]
        file_handle.close()
    else:
        with open(infile, 'r') as file_handle:
            first_symbol = file_handle.readline()[0]
    result = False
    if first_symbol == '@':
        result = True
    return result
