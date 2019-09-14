"""Runs Fama functional profiling pipeline"""
import os
from subprocess import Popen, PIPE, CalledProcessError

from lib.utils.const import ENDS, STATUS_GOOD
from lib.project.project import Project
from lib.project.sample import Sample
from lib.diamond_parser.diamond_parser import DiamondParser
from lib.output.report import generate_fastq_report, generate_sample_report
from lib.output.pdf_report import generate_pdf_report
from lib.output.krona_xml_writer import generate_functions_chart
from lib.output.json_util import export_annotated_reads, export_sample
from lib.third_party.microbe_census import run_pipeline, report_results

def run_ref_search(parser, command):
    """Runs pre-selection DIAMOND search

    Args:
        parser (:obj:DiamondParser): parser object processing an input sequence file
        command (str): either 'blastx' or 'blastp' (see DIAMOND manual)
    """
    print('Starting DIAMOND')
    diamond_args = [parser.config.diamond_path,
                    command,
                    '--db',
                    parser.config.get_reference_diamond_db(parser.options\
                    .get_collection(parser.sample.sample_id)),
                    '--query',
                    parser.options.get_fastq_path(parser.sample.sample_id, parser.end),
                    '--out',
                    os.path.join(parser.options.get_project_dir(parser.sample.sample_id),
                                 parser.sample.sample_id + '_' + parser.end + '_' \
                                 + parser.options.ref_output_name),
                    '--max-target-seqs',
                    '50',
                    '--evalue',
                    str(parser.config.get_evalue_cutoff(parser.options\
                        .get_collection(parser.sample.sample_id))),
                    '--threads',
                    parser.config.threads,
                    '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length',
                    'mismatch', 'slen', 'qstart', 'qend', 'sstart', 'send',
                    'evalue', 'bitscore'
                   ]

    with Popen(diamond_args, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
        for line in proc.stdout:
            print(line, end='')
    if proc.returncode != 0:
        raise CalledProcessError(proc.returncode, proc.args)
    print('DIAMOND finished')

def run_bgr_search(parser, command):
    """Runs classification DIAMOND search

    Args:
        parser (:obj:DiamondParser): parser object processing an input sequence file
        command (str): either 'blastx' or 'blastp' (see DIAMOND manual)
    """
    print('Starting DIAMOND')
    diamond_args = [parser.config.diamond_path,
                    command,
                    '--db',
                    parser.config.get_background_diamond_db(parser.options\
                    .get_collection(parser.sample.sample_id)),
                    '--query',
                    os.path.join(parser.options.get_project_dir(parser.sample.sample_id),
                                 parser.sample.sample_id + '_' + parser.end + '_' \
                                 + parser.options.ref_hits_fastq_name),
                    '--out',
                    os.path.join(parser.options.get_project_dir(parser.sample.sample_id),
                                 parser.sample.sample_id + '_' + parser.end + '_' \
                                 + parser.options.background_output_name),
                    '--max-target-seqs',
                    '100',
                    '--evalue',
                    str(parser.config.get_background_db_size(parser.options\
                    .get_collection(parser.sample.sample_id)) \
                    * parser.config.get_evalue_cutoff(parser.options\
                    .get_collection(parser.sample.sample_id)) \
                    / parser.config.get_reference_db_size(parser.options\
                    .get_collection(parser.sample.sample_id))),
                    '--threads',
                    parser.config.threads,
                    '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                    'slen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
                   ]

    with Popen(diamond_args, stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
        for line in proc.stdout:
            print(line, end='')
    if proc.returncode != 0:
        raise CalledProcessError(proc.returncode, proc.args)
    print('DIAMOND finished')

def run_microbecensus(sample, config):
    """Runs MicrobeCensus

    Args:
        sample (:obj:Sample): sample analyzed
        config (:obj:ProgramConfig): program configuration object
    """
    args = {}
    if sample.is_paired_end:
        args['seqfiles'] = [sample.fastq_fwd_path, sample.fastq_rev_path]
    else:
        args['seqfiles'] = [sample.fastq_fwd_path]
    args['verbose'] = True
    args['diamond'] = config.diamond_path
    args['data_dir'] = config.microbecensus_datadir
    args['outfile'] = os.path.join(sample.work_directory, 'microbecensus.out.txt')
    args['threads'] = int(config.threads)
    args['no_equivs'] = True
    if sample.fastq_fwd_readcount < 1500000:
        # MicrobeCensus subsamples 2M reads by default, but sequence library
        # must have more reads as some reads are always discarded by filtering
        args['nreads'] = sample.fastq_fwd_readcount // 2
    elif sample.fastq_fwd_readcount < 3000000:
        args['nreads'] = sample.fastq_fwd_readcount - 1000000
    else:
        args['nreads'] = 2000000
    print(args)
    est_ags, args = run_pipeline(args)
    report_results(args, est_ags, None)


def fastq_pipeline(config_file, project_file, sample_identifier=None, end_identifier=None):
    """Functional profiling pipeline for entire project

    Args:
        config_file (str): program configuration file
        project_file (str): project parameters file
        sample_identifier (str): sample identifier
        end_identifier (str): end identifier
    """
    project = Project(config_file=config_file, project_file=project_file)

    for sample_id in project.list_samples():
        if sample_identifier and sample_identifier != sample_id:
            continue
        sample = Sample(sample_id)
        sample.load_sample(project.options)
        project.samples[sample_id] = sample
        for end in ENDS:
            if end_identifier is not None and end != end_identifier:
                continue
            if end == 'pe2' and not project.samples[sample_id].is_paired_end:
                continue
            project.samples[sample_id].reads[end] = \
                run_fastq_pipeline(project,
                                   sample=project.samples[sample_id],
                                   end_id=end)
        export_sample(project.samples[sample_id])
        # Generate output for the sample or delete sample from memory
        generate_sample_report(project, sample_id)
        project.options.set_sample_data(project.samples[sample_id])

    # Generate output for the project
    if sample_identifier is None:
        # Skip project report if the pipeline is running for only one sample
        project.generate_report()

    # Rename existing project file and save current version
    project.save_project_options()

def run_fastq_pipeline(project, sample, end_id):
    """Functional profiling pipeline for single FASTQ file processing

    Args:
        project (:obj:Project): current project
        sample (:obj:Sample): current sample
        end_id (str): end identifier
    """
    parser = DiamondParser(config=project.config,
                           options=project.options,
                           taxonomy_data=project.taxonomy_data,
                           ref_data=project.ref_data,
                           sample=sample,
                           end=end_id)

    if not os.path.isdir(project.options.get_project_dir(sample.sample_id)):
        os.makedirs(project.options.get_project_dir(sample.sample_id), exist_ok=True)
    if not os.path.isdir(os.path.join(project.options.get_project_dir(sample.sample_id),
                                      project.options.get_output_subdir(sample.sample_id))):
        os.mkdir(os.path.join(project.options.get_project_dir(sample.sample_id),
                              project.options.get_output_subdir(sample.sample_id)))

    # Search in reference database
    if not os.path.exists(os.path.join(parser.options.get_project_dir(parser.sample.sample_id),
                                       parser.sample.sample_id + '_' + parser.end + '_'\
                                       + parser.options.ref_output_name)):
        run_ref_search(parser, 'blastx')

    # Process output of reference DB search
    parser.parse_reference_output()

    #Import sequence data for selected sequence reads
    print('Reading FASTQ file')
    read_count, base_count = parser.import_fastq()

    if end_id == 'pe1':
        if sample.fastq_fwd_readcount == 0:
            sample.fastq_fwd_readcount = read_count
        if sample.fastq_fwd_basecount == 0:
            sample.fastq_fwd_basecount = base_count

    elif end_id == 'pe2':
        if sample.fastq_rev_readcount == 0:
            sample.fastq_rev_readcount = read_count
        if sample.fastq_rev_basecount == 0:
            sample.fastq_rev_basecount = base_count

    if sample.rpkg_scaling_factor == 0.0:
        sample.import_rpkg_scaling_factor()
    if sample.rpkg_scaling_factor == 0.0:
        run_microbecensus(sample=sample, config=project.config)
        sample.import_rpkg_scaling_factor()
    project.options.set_sample_data(sample)

    if not parser.reads:
        # No hits found
        return {}

    print('Exporting FASTQ ')
    parser.export_hit_fastq()
    print('Exporting hits')
    parser.export_hit_list()

    # Search in background database
    if not os.path.exists(os.path.join(parser.options.get_project_dir(parser.sample.sample_id),
                                       parser.sample.sample_id + '_' + parser.end + '_' \
                                       + parser.options.background_output_name)):
        run_bgr_search(parser, 'blastx')

    # Process output of background DB search
    parser.parse_background_output()

    parser.export_read_fastq()
    if sample.is_paired_end:
        parser.export_paired_end_reads_fastq()
    export_annotated_reads(parser)

    # Generate output
    generate_fastq_report(parser)
    generate_pdf_report(parser)
    generate_functions_chart(parser)

    #return parser.reads
    # Return only good reads
    return {read_id:read for (read_id, read) in parser.reads.items() if read.status == STATUS_GOOD}

def main():
    """Main function"""
    print('This program is not intended to run directly. Run fama.py instead.')

if __name__ == '__main__':
    main()
