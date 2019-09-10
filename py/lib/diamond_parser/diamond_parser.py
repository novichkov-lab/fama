"""Describes DiamondParser class"""
import os
import csv
import gzip

from lib.utils.const import STATUS_GOOD
from lib.project.program_config import ProgramConfig
from lib.project.project_options import ProjectOptions
from lib.reference_library.reference_data import ReferenceData
from lib.reference_library.taxonomy_data import TaxonomyData
from lib.diamond_parser.diamond_hit import DiamondHit
from lib.diamond_parser.diamond_hit_list import DiamondHitList
from lib.sequences.annotated_read import AnnotatedRead
from lib.diamond_parser.hit_utils import compare_hits_erpk_lca, get_paired_end, parse_fastq_seqid

class DiamondParser(object):
    """DiamondParser performs various operations with DIAMOND input and
    output files: imports and exports sequences, imports DIAMOND search
    results in tabular formats, stores the results etc.
    Basically, this is a workhorse of Fama pipeline.

    Processing of each FASTQ or FASTA file requires its own instance of DiamondParser.

    During a run of Fama profiling, DiamondParser methods are called in
    the following order:
    1. parse_reference_output: reads first DIAMOND output file and finds
        candidate query sequences of interest
    2. import_fastq (or import_fasta): reads initial query file and finds
        sequences of interest.
    3. export_hit_fastq (or export_hit_fasta): creates a sequence file of hits found.
    4. export_hit_list : creates a table of hits selected for subsequent analysis
    5. parse_background_output: reads seconf DIAMOND output file and finds
        which query sequences are truly sequences of interest, which should
        be counted for functional or taxonomic profile
    6. export_paired_end_reads_fastq: for paired-end reads, finds paired
        reads for sequences of intereast and exports them

    Attributes:
        reads (:obj:'dict' of :obj:'AnnotatedRead'): dictionary of query
            sequences (sequence reads or proteins)
        sample (str): sample identifier
        end (str): end identifier, always 'pe1' or 'pe2'. Single-end
            reads and protein projects use only 'pe1'
        config (:obj:'ProgramConfig'): Fama configuration parameters
        options  (:obj:'ProjectOptions'): Fama project options
        collection (str): reference collection identifier
        ref_data (:obj:'ReferenceData'): reference dataset for the
            collection (list of functions, list of proteins etc.)
        taxonomy_data (:obj:'TaxonomyData'): NCBI taxonomy dataset for the collection

    Todo:
        * add paired-end processing in DiamondParser. For paired-end
        reads, let's process two DIAMOND outputs first, then read each
        of FASTQ files only once.
        But for that, reads dictionary must have two levels: outer level
        for read  identifiers, and inner level for end identifiers.

    """
    def __init__(self, sample=None, end='', config_file=None, project_file=None,
                 config=None, options=None, ref_data=None, taxonomy_data=None):
        """
        Args:
            sample (str): sample identifier
            end (str): end identifier, always 'pe1' or 'pe2'. Single-end
                reads and protein projects use only 'pe1'
            config_file (str): full path to program configuration ini file.
                Ignored, if config argument is set
            project_file (str): full path to project ini file. Ignored,
                if options argument is set
            config (:obj:'ProgramConfig'): Fama configuration parameters
            options  (:obj:'ProjectOptions'): Fama project options
            ref_data (:obj:'ReferenceData'): reference dataset for the
                collection (list of functions, list of proteins etc.)
            taxonomy_data (:obj:'TaxonomyData'): NCBI taxonomy dataset for the collection

        Raises:
            Exception: if collection identifier from project options
                does not match any collection in program config.

        """
        self.reads = {}
        self.sample = sample
        self.end = end
        self.config = config
        if not self.config:
            self.config = ProgramConfig(config_file)
        self.options = options
        if not self.options:
            self.options = ProjectOptions(project_file)
        collection = self.options.get_collection(self.sample.sample_id)
        if collection not in self.config.collections:
            raise Exception('Collection ' + collection
                            + ' not found. Available collections are: '
                            + (',').join(self.config.collections))
        self.collection = collection
        self.ref_data = ref_data
        if self.ref_data is None:
            self.ref_data = ReferenceData(self.config)
            self.ref_data.load_reference_data(self.collection)
        self.taxonomy_data = taxonomy_data
        if self.taxonomy_data is None:
            self.taxonomy_data = TaxonomyData(self.config)
            self.taxonomy_data.load_taxdata(self.config)

    def parse_reference_output(self):
        """Reads and processes DIAMOND tabular output of the first DIAMOND
        search.

        Note: this function finds query sequences similar to reference
        proteins. Since a query sequence may have more than one areas of
        similarity (for instance, in fusion proteins of two subunits or
        in multi-domain proteins), it will try to find as many such areas
        as possible.

        DIAMOND hits are filtered by two parameters: length of alignment
        and amino acid identity %.

        This function does not return anything. Instead, it populates
        'reads' dictionary with AnnotatedRead objects.

        """
        tsvfile = os.path.join(self.options.get_project_dir(self.sample.sample_id),
                               self.sample.sample_id + '_' + self.end
                               + '_'+ self.options.ref_output_name)
        current_sequence_read_id = ''
        hit_list = DiamondHitList(current_sequence_read_id)
        identity_cutoff = self.config.get_identity_cutoff(self.collection)
        length_cutoff = self.config.get_length_cutoff(self.collection)
        print('Identity cutoff: ', identity_cutoff, ', Length cutoff: ', length_cutoff)
        with open(tsvfile, 'r', newline='') as infile:
            tsvin = csv.reader(infile, delimiter='\t')
            for row in tsvin:
                hit = DiamondHit()
                (row[0], _) = parse_fastq_seqid(row[0])
                hit.create_hit(row)
                # filtering by identity and length
                if hit.identity < identity_cutoff:
                    continue # go to next hit
                if hit.length < length_cutoff:
                    continue # go to next hit

                if hit.query_id != current_sequence_read_id:
                    # when new query ID reached, process collected hits,
                    # then start over with new query identifier
                    # filtering: remove overlapping hits
                    hit_list.filter_list(self.config.get_overlap_cutoff(self.collection))
                    # if any hits left, assign function to hits and populate reads dictionary
                    if hit_list.hits_number != 0:
                        hit_list.annotate_hits(self.ref_data)
                        read = AnnotatedRead(current_sequence_read_id)
                        read.hit_list = hit_list
                        self.reads[current_sequence_read_id] = read
                    # start over
                    current_sequence_read_id = hit.query_id
                    hit_list = DiamondHitList(current_sequence_read_id)
                hit_list.add_hit(hit)
            # when EOF reached, process collected hits
            if hit_list.hits_number != 0:
                hit_list.filter_list(self.config.get_overlap_cutoff(self.collection))
                hit_list.annotate_hits(self.ref_data)
                read = AnnotatedRead(current_sequence_read_id)
                read.hit_list = hit_list
                self.reads[current_sequence_read_id] = read

    def parse_background_output(self):
        """Reads and processes DIAMOND tabular output of the second DIAMOND
        search.

        Note: this function takes existing list of hits and compares each
        of them with results of other similarity serach (against larger DB).
        For the comparison, it calls compare_hits_erpk_lca function, which
        in turn updates entries in the 'reads' dictionary.

        Raises:
            KeyError if read identifier not found in the 'reads' dictionary

        """
        if not self.reads:
            # Something went wrong and 'reads' dictionary is empty.
            # Let's try to import list of reads from file.
            self.reads = self.import_hit_list()

        tsvfile = os.path.join(self.sample.work_directory,
                               self.sample.sample_id + '_' + self.end
                               + '_'+ self.options.background_output_name)

        average_read_length = self.sample.get_avg_read_length(self.end)

        current_query_id = None
        hit_list = None
        identity_cutoff = self.config.get_identity_cutoff(self.collection)
        length_cutoff = self.config.get_length_cutoff(self.collection)
        bitscore_range_cutoff = self.config.get_biscore_range_cutoff(self.collection)
        print('Identity cutoff: ', identity_cutoff, ', Length cutoff: ', length_cutoff)

        with open(tsvfile, 'r', newline='') as infile:
            tsvin = csv.reader(infile, delimiter='\t')
            for row in tsvin:
                if current_query_id is None:
                    current_query_id = row[0]
                    hit_list = DiamondHitList(current_query_id)
                hit = DiamondHit()
                hit.create_hit(row)
                # filtering by identity and length
                if hit.identity < identity_cutoff:
                    continue # skip this line
                if hit.length < length_cutoff:
                    continue # skip this line

                # when new query ID reached, process collected hits,
                # then start over with new query identifier
                if hit.query_id != current_query_id:
                    # assign functions to selected hits
                    hit_list.annotate_hits(self.ref_data)
                    # extract initial read identifier from identifier of the hit
                    current_query_id_tokens = current_query_id.split('|')
                    read_id = '|'.join(current_query_id_tokens[:-2])
                    # compare list of hits from search in background DB
                    # with existing hit from the first similarity search
                    try:
                        compare_hits_erpk_lca(
                            self.reads[read_id],
                            int(current_query_id_tokens[-2]), #hit_start
                            int(current_query_id_tokens[-1]), #hit_end
                            hit_list, bitscore_range_cutoff, length_cutoff,
                            average_read_length, self.taxonomy_data, self.ref_data,
                            rank_cutoffs=self.config.get_ranks_cutoffs(
                                self.options.get_collection()
                                )
                            )
                    except KeyError:
                        print('Read not found: ', read_id)
                    # starting over
                    current_query_id = hit.query_id
                    hit_list = DiamondHitList(current_query_id)
                hit_list.add_hit(hit)
            # when EOF reached, process collected hits
            # assign functions to selected hits
            hit_list.annotate_hits(self.ref_data)
            # extract initial read identifier from identifier of the hit
            current_query_id_tokens = current_query_id.split('|')
            read_id = '|'.join(current_query_id_tokens[:-2])
            # compare list of hits from search in background DB with
            # existing hit from the first similarity search
            try:
                compare_hits_erpk_lca(
                    self.reads[read_id],
                    int(current_query_id_tokens[-2]), #hit_start
                    int(current_query_id_tokens[-1]), #hit_end
                    hit_list, bitscore_range_cutoff, length_cutoff,
                    average_read_length, self.taxonomy_data, self.ref_data,
                    rank_cutoffs=self.config.get_ranks_cutoffs(
                        self.options.get_collection()
                        )
                    )
            except KeyError:
                print('Read not found: ', read_id)

    def import_fastq(self):
        """Reads uncompressed or gzipped FASTQ file, finds sequences of
        selected reads and stores them

        Returns:
            read_count (int): number of reads in the file
            base_count (int): total number of bases in all reads
        """
        fastq_file = self.options.get_fastq_path(self.sample.sample_id, self.end)
        line_counter = 0
        read_count = 0
        base_count = 0
        current_read = None
        infile_handle = None
        if fastq_file.endswith('.gz'):
            infile_handle = gzip.open(fastq_file, 'rb')
        else:
            infile_handle = open(fastq_file, 'rb')
        for line in infile_handle:
            # count lines as each FASTQ entry has exactly four lines
            line_counter += 1
            if line_counter == 5:
                line_counter = 1
            line = line.decode('utf8').rstrip('\n\r')
            if line_counter == 1:
                read_count += 1
                (read_id, _) = parse_fastq_seqid(line)
                current_read = None
                if read_id in self.reads:
                    current_read = read_id
                    self.reads[current_read].read_id_line = line
            elif line_counter == 2:
                base_count += len(line)
                if current_read is None:
                    continue
                self.reads[current_read].sequence = line
            elif line_counter == 3:
                if current_read is None:
                    continue
                self.reads[current_read].line3 = line
            elif line_counter == 4:
                if current_read is None:
                    continue
                self.reads[current_read].quality = line
        infile_handle.close()
        return read_count, base_count

    def import_fasta(self):
        """Reads uncompressed or gzipped FASTA file, finds sequences of
        selected reads and stores them

        Returns:
            read_count (int): number of reads in the file
            base_count (int): total number of bases in all reads
        """
        fasta_file = self.options.get_fastq_path(self.sample.sample_id, self.end)
        sequence = []
        read_count = 0
        base_count = 0
        current_id = ''
        infile_handle = None
        if fasta_file.endswith('.gz'):
            infile_handle = gzip.open(fasta_file, 'rb')
        else:
            infile_handle = open(fasta_file, 'rb')
        for line in infile_handle:
            line = line.decode('utf8').rstrip('\n\r')
            if line.startswith('>'):
                read_count += 1
                if current_id != '':
                    current_id = current_id[1:]
                    self.reads[current_id].read_id_line = current_id
                    self.reads[current_id].sequence = ''.join(sequence)
                sequence = []
                seq_id = line[1:]
                if seq_id in self.reads:
                    current_id = line
                else:
                    current_id = ''
            else:
                base_count += len(line)
                if current_id == '':
                    continue
                sequence.append(line)
        if current_id != '':
            self.reads[seq_id].read_id_line = current_id
            self.reads[seq_id].sequence = ''.join(sequence)
        infile_handle.close()
        return read_count, base_count

    def export_read_fastq(self):
        """Exports sequence reads as gzipped FASTQ file"""
        outdir = self.sample.work_directory
        with gzip.open(os.path.join(outdir,
                                    self.sample.sample_id + '_'
                                    + self.end + '_'
                                    + self.options.reads_fastq_name + '.gz'), 'wt') as outfile:
            for read_id in sorted(self.reads.keys()):
                if self.reads[read_id].status == STATUS_GOOD:
                    outfile.write(self.reads[read_id].read_id_line + '\n')
                    outfile.write(self.reads[read_id].sequence + '\n')
                    outfile.write(self.reads[read_id].line3 + '\n')
                    outfile.write(self.reads[read_id].quality + '\n')

    def export_read_fasta(self):
        """Exports sequences as gzipped FASTA file"""
        outdir = self.sample.work_directory
        fastq_file = os.path.join(outdir,
                                  self.sample.sample_id + '_'
                                  + self.end + '_'
                                  + self.options.reads_fastq_name + '.gz')
        with gzip.open(fastq_file, 'wt') as outfile:
            for read_id in sorted(self.reads.keys()):
                if self.reads[read_id].status == STATUS_GOOD:
                    outfile.write(self.reads[read_id].read_id_line + '\n')
                    outfile.write(self.reads[read_id].sequence + '\n')

    def export_hit_fastq(self):
        """Exports sequences of DAIMOND hits as gzipped FASTQ file"""
        outdir = self.sample.work_directory
        with open(os.path.join(outdir,
                               self.sample.sample_id + '_'
                               + self.end + '_'
                               + self.options.ref_hits_fastq_name), 'w') as outfile:
            for read_id, read in self.reads.items():
                for hit in read.hit_list.hits:
                    start = hit.q_start
                    end = hit.q_end
                    outfile.write("@" + read.read_id + '|' + \
                        str(start) + '|' + str(end) + '\n')
                    if start < end:
                        # hit on + strand
                        start = start - 1
                        end = end
                    else:
                        # hit on - strand
                        temp_value = start
                        start = end - 1
                        end = temp_value
                    try:
                        outfile.write(read.sequence[start:end] + '\n')
                        outfile.write(read.line3 + '\n')
                        outfile.write(read.quality[start:end] + '\n')
                    except TypeError:
                        print('TypeError occurred while exporting ', read_id)

    def export_hit_fasta(self):
        """Exports sequences of DAIMOND hits as gzipped FASTA file"""
        outdir = self.sample.work_directory
        with open(os.path.join(outdir,
                               self.sample.sample_id + '_'
                               + self.end + '_'
                               + self.options.ref_hits_fastq_name), 'w') as outfile:
            for read_id, read in self.reads.items():
                for hit in read.hit_list.hits:
                    start = hit.q_start
                    end = hit.q_end
                    outfile.write(">" + read.read_id + '|' + \
                        str(start) + '|' + str(end) + '\n')
                    if start < end:
                        # hit on + strand
                        start = start - 1
                        end = end
                    else:
                        # hit on - strand
                        temp_value = start
                        start = end - 1
                        end = temp_value
                    try:
                        outfile.write(read.sequence[start:end] + '\n')
                    except TypeError:
                        print('TypeError occurred while exporting ', read_id)

    def export_hit_list(self):
        """Exports tab-separated table of DAIMOND hits"""
        outfile = os.path.join(self.sample.work_directory,
                               self.sample.sample_id + '_'
                               + self.end + '_'
                               + self.options.ref_hits_list_name)
        with open(outfile, 'w') as outfile:
            for _, read in self.reads.items():
                for hit in read.hit_list.hits:
                    outfile.write(str(hit) + '\n')

    def export_paired_end_reads_fastq(self):
        """ For paired-end sequence reads, reads FASTQ file,
        collects paired-end reads for pre-selected reads and writes them
        to a separate FASTQ file

        """
        fastq_file = self.options.get_fastq_path(self.sample.sample_id, get_paired_end(self.end))
        outdir = self.sample.work_directory
        read_ids = {}
        for read_id in sorted(self.reads.keys()):
            read_ids[read_id] = read_id
        line_counter = 0
        fastq_outfile = os.path.join(outdir,
                                     self.sample.sample_id + '_'
                                     + self.end + '_'
                                     + self.options.pe_reads_fastq_name + '.gz')
        with gzip.open(fastq_outfile, 'wt') as outfile:
            current_read = None
            infile_handle = None
            if fastq_file.endswith('.gz'):
                infile_handle = gzip.open(fastq_file, 'rb')
            else:
                infile_handle = open(fastq_file, 'rb')
            for line in infile_handle:
                line_counter += 1
                if line_counter == 5:
                    line_counter = 1
                line = line.decode('utf8').rstrip('\n\r')
                if line_counter == 1:
                    current_read = None
                    (read_id, _) = parse_fastq_seqid(line)
                    if read_id in read_ids:
                        if self.reads[read_id].status == STATUS_GOOD:
                            current_read = read_id
                            self.reads[current_read].pe_id = line
                            outfile.write(line + '\n')
                elif current_read is not None:
                    outfile.write(line + '\n')
                    if line_counter == 2:
                        self.reads[current_read].pe_sequence = line
                    elif line_counter == 3:
                        self.reads[current_read].pe_line3 = line
                    elif line_counter == 4:
                        self.reads[current_read].pe_quality = line
            infile_handle.close()

    def import_hit_list(self):
        """Imports tab-separated table of DAIMOND hits. Use for resuming
        analysis after Fama restart

        Returns:
            :obj:dict[str, :obj:AnnotatedRead]

        """
        infile = os.path.join(os.path.join(self.sample.work_directory,
                                           self.sample.sample_id + '_'
                                           + self.end + '_'
                                           + self.options.ref_hits_list_name))
        ret_val = {}
        hit_list = None
        current_read_id = None
        with open(infile, 'r', newline='') as infile:
            tsvin = csv.reader(infile, delimiter='\t')
            for row in tsvin:
                if current_read_id is None:
                    # initialize
                    current_read_id = row[0]
                    hit_list = DiamondHitList(current_read_id)
                elif current_read_id != row[0]:
                    ret_val[current_read_id] = AnnotatedRead(current_read_id)
                    ret_val[current_read_id].hit_list = hit_list
                    current_read_id = row[0]
                    hit_list = DiamondHitList(current_read_id)
                hit = DiamondHit()
                hit.import_hit(row)
                hit_list.add_hit(hit)
            ret_val[current_read_id] = AnnotatedRead(current_read_id)
            ret_val[current_read_id].hit_list = hit_list

        return ret_val
