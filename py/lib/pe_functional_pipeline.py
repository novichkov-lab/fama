"""Runs Fama functional profiling pipeline"""
import os
import gzip

from lib.utils.const import ENDS, STATUS_GOOD
from lib.se_functional_pipeline import run_fastq_pipeline
from lib.utils.utils import run_external_program
from lib.project.project import Project
from lib.project.sample import Sample
from lib.diamond_parser.diamond_parser import DiamondParser
from lib.output.report import generate_fastq_report, generate_sample_report
from lib.output.pdf_report import generate_pdf_report
from lib.output.krona_xml_writer import make_functions_chart
from lib.output.json_util import export_annotated_reads, export_sample
from lib.third_party.microbe_census import run_pipeline, report_results
from lib.diamond_parser.hit_utils import get_paired_end, parse_fastq_seqid

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
                    parser.config.get_reference_diamond_db(
                        parser.options.get_collection(parser.sample.sample_id)
                    ),
                    '--query',
                    parser.options.get_fastq_path(parser.sample.sample_id, parser.end),
                    '--out',
                    os.path.join(
                        parser.options.get_project_dir(parser.sample.sample_id),
                        parser.sample.sample_id + '_' + parser.end + '_'
                        + parser.options.ref_output_name
                    ),
                    '--max-target-seqs',
                    '50',
                    '--evalue',
                    str(parser.config.get_evalue_cutoff(
                        parser.options.get_collection(parser.sample.sample_id)
                    )),
                    '--threads',
                    parser.config.threads,
                    '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length',
                    'mismatch', 'slen', 'qstart', 'qend', 'sstart', 'send',
                    'evalue', 'bitscore']
    run_external_program(diamond_args)
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
                    parser.config.get_background_diamond_db(
                        parser.options.get_collection(parser.sample.sample_id)
                    ),
                    '--query',
                    os.path.join(
                        parser.options.get_project_dir(parser.sample.sample_id),
                        parser.sample.sample_id + '_' + parser.end + '_'
                        + parser.options.ref_hits_fastq_name
                    ),
                    '--out',
                    os.path.join(
                        parser.options.get_project_dir(parser.sample.sample_id),
                        parser.sample.sample_id + '_' + parser.end + '_'
                        + parser.options.background_output_name
                    ),
                    '--max-target-seqs',
                    '100',
                    '--evalue',
                    str(
                        parser.config.get_background_db_size(
                            parser.options.get_collection(parser.sample.sample_id)
                        ) * parser.config.get_evalue_cutoff(
                            parser.options.get_collection(parser.sample.sample_id)
                        ) / parser.config.get_reference_db_size(
                            parser.options.get_collection(parser.sample.sample_id)
                        )),
                    '--threads',
                    parser.config.threads,
                    '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                    'slen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    run_external_program(diamond_args)
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


def import_fastq_pe(parser1, parser2):
    """Reads uncompressed or gzipped FASTQ file, finds sequences of
    selected reads and stores them

    Returns:
        read_count (int): number of reads in the file
        base_count (int): total number of bases in all reads
    """
    fastq_file1 = parser1.options.get_fastq_path(parser1.sample.sample_id, parser1.end)
    line_counter = 0
    read_count1 = 0
    base_count1 = 0
    current_read = None
    infile_handle = None
    if fastq_file1.endswith('.gz'):
        infile_handle = gzip.open(fastq_file1, 'rb')
    else:
        infile_handle = open(fastq_file1, 'rb')
    for line in infile_handle:
        # count lines as each FASTQ entry has exactly four lines
        line_counter += 1
        if line_counter == 5:
            line_counter = 1
        line = line.decode('utf8').rstrip('\n\r')
        if line_counter == 1:
            read_count1 += 1
            (read_id, _) = parse_fastq_seqid(line)
            current_read = read_id
            if current_read in parser1.reads:
                parser1.reads[current_read].read_id_line = line
            if current_read in parser2.reads:
                parser2.reads[current_read].pe_id = line
        elif line_counter == 2:
            base_count1 += len(line)
            if current_read in parser1.reads:
                parser1.reads[current_read].sequence = line
            if current_read in parser2.reads:
                parser2.reads[current_read].pe_sequence = line
        elif line_counter == 3:
            if current_read in parser1.reads:
                parser1.reads[current_read].line3 = line
            if current_read in parser2.reads:
                parser2.reads[current_read].pe_line3 = line
        elif line_counter == 4:
            if current_read in parser1.reads:
                parser1.reads[current_read].quality = line
            if current_read in parser2.reads:
                parser2.reads[current_read].pe_quality = line
    infile_handle.close()

    fastq_file2 = parser1.options.get_fastq_path(parser2.sample.sample_id, parser2.end)
    line_counter = 0
    read_count2 = 0
    base_count2 = 0
    current_read = None
    if fastq_file2.endswith('.gz'):
        infile_handle = gzip.open(fastq_file2, 'rb')
    else:
        infile_handle = open(fastq_file2, 'rb')
    for line in infile_handle:
        # count lines as each FASTQ entry has exactly four lines
        line_counter += 1
        if line_counter == 5:
            line_counter = 1
        line = line.decode('utf8').rstrip('\n\r')
        if line_counter == 1:
            read_count2 += 1
            (read_id, _) = parse_fastq_seqid(line)
            current_read = read_id
            if current_read in parser1.reads:
                parser1.reads[current_read].pe_id = line
            if current_read in parser2.reads:
                parser2.reads[current_read].read_id_line = line
        elif line_counter == 2:
            base_count2 += len(line)
            if current_read in parser1.reads:
                parser1.reads[current_read].pe_sequence = line
            if current_read in parser2.reads:
                parser2.reads[current_read].sequence = line
        elif line_counter == 3:
            if current_read in parser1.reads:
                parser1.reads[current_read].pe_line3 = line
            if current_read in parser2.reads:
                parser2.reads[current_read].line3 = line
        elif line_counter == 4:
            if current_read in parser1.reads:
                parser1.reads[current_read].pe_quality = line
            if current_read in parser2.reads:
                parser2.reads[current_read].quality = line
    infile_handle.close()
    return (parser1, parser2, read_count1, read_count2, base_count1, base_count2)


def export_paired_end_reads_fastq(parser):
    """ For paired-end sequence reads, write paired-end reads for pre-selected 
    reads into a separate FASTQ file
    """
    fastq_file = parser.options.get_fastq_path(parser.sample.sample_id, get_paired_end(parser.end))
    outdir = parser.sample.work_directory
    read_ids = {}
    for read_id in sorted(parser.reads.keys()):
        read_ids[read_id] = read_id
    line_counter = 0
    fastq_outfile = os.path.join(outdir,
                                 parser.sample.sample_id + '_'
                                 + parser.end + '_'
                                 + parser.options.pe_reads_fastq_name + '.gz')
    with gzip.open(fastq_outfile, 'wt') as outfile:
        for read_id in sorted(parser.reads.keys()):
            outfile.write(parser.reads[read_id].pe_id + '\n')
            outfile.write(parser.reads[read_id].pe_sequence + '\n')
            outfile.write(parser.reads[read_id].pe_line3 + '\n')
            outfile.write(parser.reads[read_id].pe_quality + '\n')


def fastq_pe_pipeline(project, sample_identifier=None, end_identifier=None):
    """Functional profiling pipeline for entire project

    Args:
        project (:obj:Project): current project
        sample_identifier (str, optional): sample identifier
        end_identifier (str, optional): end identifier
    """
    for sample_id in project.list_samples():
        if sample_identifier and sample_identifier != sample_id:
            continue
        sample = Sample(sample_id)
        sample.load_sample(project.options)
        project.samples[sample_id] = sample
        if end_identifier:
            project.samples[sample_id].reads[end] = \
                run_fastq_pipeline(project,
                                   sample=project.samples[sample_id],
                                   end_id=end_identifier)
        else:
            project.samples[sample_id].reads = \
                run_pe_fastq_pipeline(project,
                                   sample=project.samples[sample_id])
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


def run_pe_fastq_pipeline(project, sample):
    """Functional profiling pipeline for single FASTQ file processing

    Args:
        project (:obj:Project): current project
        sample (:obj:Sample): current sample
    """
    result = {}
    parser1 = DiamondParser(config=project.config,
                           options=project.options,
                           taxonomy_data=project.taxonomy_data,
                           ref_data=project.ref_data,
                           sample=sample,
                           end=ENDS[0])
    parser2 = DiamondParser(config=project.config,
                           options=project.options,
                           taxonomy_data=project.taxonomy_data,
                           ref_data=project.ref_data,
                           sample=sample,
                           end=ENDS[1])

    if not os.path.isdir(project.options.get_project_dir(sample.sample_id)):
        os.makedirs(project.options.get_project_dir(sample.sample_id), exist_ok=True)
    if not os.path.isdir(os.path.join(project.options.get_project_dir(sample.sample_id),
                                      project.options.get_output_subdir(sample.sample_id))):
        os.mkdir(os.path.join(project.options.get_project_dir(sample.sample_id),
                              project.options.get_output_subdir(sample.sample_id)))

    # Search in reference database
    if not os.path.exists(
            os.path.join(
                parser1.options.get_project_dir(parser1.sample.sample_id),
                parser1.sample.sample_id + '_' + parser1.end + '_' + parser1.options.ref_output_name
            )
    ):
        run_ref_search(parser1, 'blastx')

    if not os.path.exists(
            os.path.join(
                parser2.options.get_project_dir(parser2.sample.sample_id),
                parser2.sample.sample_id + '_' + parser2.end + '_' + parser2.options.ref_output_name
            )
    ):
        run_ref_search(parser2, 'blastx')

    # Process output of reference DB search
    parser1.parse_reference_output()
    parser2.parse_reference_output()


    parser1_read_ids = set(parser1.reads.keys())
    parser2_read_ids = set(parser2.reads.keys())

    # Import sequence data for selected sequence reads
    print('Reading FASTQ file')
    (parser1, parser2, read_count1, read_count2, base_count1, base_count2) = import_fastq_pe(parser1, parser2)

    if sample.fastq_fwd_readcount == 0:
        sample.fastq_fwd_readcount = read_count1
    if sample.fastq_fwd_basecount == 0:
        sample.fastq_fwd_basecount = base_count1
    if sample.fastq_rev_readcount == 0:
        sample.fastq_rev_readcount = read_count2
    if sample.fastq_rev_basecount == 0:
        sample.fastq_rev_basecount = base_count2

    if sample.rpkg_scaling_factor == 0.0:
        sample.import_rpkg_scaling_factor()
    if sample.rpkg_scaling_factor == 0.0:
        run_microbecensus(sample=sample, config=project.config)
        sample.import_rpkg_scaling_factor()
    project.options.set_sample_data(sample)

    if parser1.reads:
        parser1.export_hit_fastq()
        print('Hits for forward end reads exported in FASTQ format')
        parser1.export_hit_list()
        print('List of hits fo forward end reads exported')
        if not os.path.exists(
                os.path.join(
                    parser1.options.get_project_dir(parser1.sample.sample_id),
                    parser1.sample.sample_id + '_' + parser1.end + '_'
                    + parser1.options.background_output_name
                )
        ):
            run_bgr_search(parser1, 'blastx')
        print('Classification DB search finished')
        parser1.parse_background_output()
        print('Classification DB search results imported')
        parser1.export_read_fastq()
        print('Classified forward end reads exported in FASTQ format')
        export_paired_end_reads_fastq(parser1)
        print('Paired reads for classified forward end reads exported')
        export_annotated_reads(parser1)
        print('Classified forward end reads exported in JSON format')
        generate_fastq_report(parser1)
        print('Text report for forward end reads created')
        generate_pdf_report(parser1)
        print('PDF report for forward end reads created')
        make_functions_chart(parser1)
        print('Krona chart for forward end reads created')
        result[ENDS[0]] = {read_id: read for (read_id, read) in parser1.reads.items() if read.status == STATUS_GOOD}
    else:
        # No hits found
        print('Pre-selection search did not find any hits for forward end reads')
        result[ENDS[0]] = {}

    if parser2.reads:
        parser2.export_hit_fastq()
        print('Hits for reverse end reads exported in FASTQ format')
        parser2.export_hit_list()
        print('List of hits for reverse end reads exported')
        if not os.path.exists(
                os.path.join(
                    parser2.options.get_project_dir(parser2.sample.sample_id),
                    parser2.sample.sample_id + '_' + parser2.end + '_'
                    + parser2.options.background_output_name
                )
        ):
            run_bgr_search(parser2, 'blastx')
        print('Classification DB search for reverse end reads finished')
        parser2.parse_background_output()
        print('Classification DB search results for reverse end reads imported')
        parser2.export_read_fastq()
        print('Classified reverse end reads exported in FASTQ format')
        export_paired_end_reads_fastq(parser2)
        print('Paired reads for classified reverse end reads exported')
        export_annotated_reads(parser2)
        print('Classified reverse end reads exported in JSON format')
        generate_fastq_report(parser2)
        print('Text report for reverse end reads created')
        generate_pdf_report(parser2)
        print('PDF report for reverse end reads created')
        make_functions_chart(parser2)
        print('Krona chart for reverse end reads created')
        result[ENDS[1]] = {read_id: read for (read_id, read) in parser2.reads.items() if read.status == STATUS_GOOD}
    else:
        # No hits found
        print('Pre-selection search did not find any hits for reverse end reads')
        result[ENDS[1]] = {}


    #~ print('Exporting FASTQ ')
    #~ parser2.export_hit_fastq()
    #~ print('Exporting hits')
    #~ parser2.export_hit_list()

    #~ # Search in background database
    #~ if not os.path.exists(
            #~ os.path.join(
                #~ parser2.options.get_project_dir(parser2.sample.sample_id),
                #~ parser2.sample.sample_id + '_' + parser2.end + '_'
                #~ + parser2.options.background_output_name
            #~ )
    #~ ):
        #~ run_bgr_search(parser2, 'blastx')

    #~ # Process output of background DB search
    #~ parser2.parse_background_output()

    #~ parser2.export_read_fastq()
    
    #~ # TODO: make new export_paired_end_reads_fastq
    #~ export_paired_end_reads_fastq(parser2)
    #~ export_annotated_reads(parser2)

    #~ # Generate output
    #~ generate_fastq_report(parser2)
    #~ generate_pdf_report(parser2)
    #~ make_functions_chart(parser2)

    return result

def main():
    """Main function"""
    print('This program is not intended to run directly. Run fama.py instead.')

if __name__ == '__main__':
    main()
