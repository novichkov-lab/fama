"""This module describes GeneAssembler class and contains some
accessory functions for running gene-centric assembler
"""
import os
import csv
import shutil
from collections import Counter, defaultdict

from lib.utils.const import ENDS, ROOT_TAXONOMY_ID, STATUS_GOOD, STATUS_BAD
from lib.utils.utils import autovivify, run_external_program, run_external_program_ignoreerror
from lib.gene_assembler.contig import Contig
from lib.gene_assembler.gene import Gene
from lib.gene_assembler.gene_assembly import GeneAssembly
from lib.diamond_parser.diamond_hit_list import DiamondHitList
from lib.diamond_parser.diamond_hit import DiamondHit
from lib.output.json_util import export_gene_assembly
from lib.taxonomy.taxonomy_profile import TaxonomyProfile
from lib.output.krona_xml_writer import make_assembly_taxonomy_chart
from lib.output.report import generate_assembly_report
from lib.output.xlsx_util import make_assembly_xlsx
from lib.reference_library.taxonomy_data import TaxonomyData


class GeneAssembler(object):
    """GeneAssembler is a working horse of Fama assembly pipeline. It exports
    sequence reads, feeds external assembler with them, imports resulting
    contigs, maps reads to contigs with Bowtie, finds genes with Prodigal,
    assigns functions to the genes and sends gene assembly data to report generator

    Attributes:
        project (:obj:Project): Project instance storing sample data,
            reference data and reads for assembly
        assembler (str): external asembly program, valid values are 'metaspades'
            (default) and 'megahit'
        assembly (:obj:GeneAssembly): gene assembly with contigs and genes
        is_paired_end (bool): True for paired-end project, False for others
        assembly_dir (path): directory for assembly files
    """
    def __init__(self, project, assembler='metaspades'):
        """Args:
            project (:obj:Project): Project instance storing sample data,
                reference data and reads for assembly
        assembler (str): external asembly program, valid values are 'metaspades'
            (default) and 'megahit'
        """
        self.project = project
        self.assembler = assembler
        project.load_project()
        self.assembly = GeneAssembly()
        self.is_paired_end = None
        self.assembly_dir = self.project.options.assembly_dir
        #~ if os.path.exists(self.assembly_dir):
            #~ raise FileExistsError ('Assembly subdirectory already exists. Delete existing directory or change subdirectory name.')
        if not os.path.isdir(self.assembly_dir):
            os.mkdir(self.assembly_dir)
        if not os.path.isdir(os.path.join(self.assembly_dir, 'out')):
            os.mkdir(os.path.join(self.assembly_dir, 'out'))

    def export_reads(self, do_coassembly=True):
        """Exports annotated reads in FASTQ format.

        Args:
            do_coassembly (bool): if True, all reads are exported into
                Coassembly_1.fastq (and Coassembly_2.fastq for paired-end
                reads) files. If False, for each function a separate file
                (or pair of files) will be created.
        """
        # Delete all existing FASTQ files
        for filename in os.listdir(self.assembly_dir):
            if filename.endswith('.fastq'):
                os.remove(os.path.join(self.assembly_dir, filename))
        # Load reads, export reads in FASTQ format, remove reads from memory
        for sample_id in sorted(self.project.list_samples()):
            self.is_paired_end = self.project.samples[sample_id].is_paired_end
            self.project.import_reads_json(sample_id, ENDS)
            for end in ENDS:
                # print ('Loading mapped reads: ', sample, end)
                # self.project.load_annotated_reads(sample, end) # Lazy load
                for read_id in self.project.samples[sample_id].reads[end]:
                    read = self.project.samples[sample_id].reads[end][read_id]
                    if read.status != STATUS_GOOD:
                        continue
                    for function in read.functions:
                        if do_coassembly:
                            function = 'Coassembly'
                        if read_id in self.assembly.reads[function]:
                            continue
                        self.assembly.reads[function][read_id] = sample_id
                        fwd_outfile = os.path.join(self.assembly_dir, function + '_pe1.fastq')
                        rev_outfile = os.path.join(self.assembly_dir, function + '_pe2.fastq')
                        if self.is_paired_end:
                            if end == 'pe2':
                                fwd_outfile = os.path.join(self.assembly_dir,
                                                           function + '_pe2.fastq')
                                rev_outfile = os.path.join(self.assembly_dir,
                                                           function + '_pe1.fastq')
                            with open(rev_outfile, 'a') as rev_of:
                                rev_of.write(read.pe_id + '\n')
                                rev_of.write(read.pe_sequence + '\n')
                                rev_of.write(read.pe_line3 + '\n')
                                rev_of.write(read.pe_quality + '\n')
                        with open(fwd_outfile, 'a') as fwd_of:
                            fwd_of.write(read.read_id_line + '\n')
                            fwd_of.write(read.sequence + '\n')
                            fwd_of.write(read.line3 + '\n')
                            fwd_of.write(read.quality + '\n')
                # Delete reads from memory
                self.project.samples[sample_id].reads[end] = None

    def assemble_contigs(self):
        """Assembles contigs from annotated reads, a separate assembly for
        each of functions, runs read mapping, calculates read coverage
        """
        # Run Assembler ('megahit' for Megahit or 'metaSPAdes' for metaSPAdes)
        if self.assembler == 'megahit':
            run_assembler(sorted(self.assembly.reads.keys()),
                          self.project.config.megahit_path,
                          self.assembly_dir)
        elif self.assembler == 'metaspades':
            run_assembler(sorted(self.assembly.reads.keys()),
                          self.project.config.metaspades_path,
                          self.assembly_dir,
                          is_paired_end=self.is_paired_end)
        else:
            raise ValueError('Unknown assembler: ' + self.assembler)

        # Filter contigs by size
        self.filter_contigs_by_length()
        # Run Bowtie
        run_mapper_indexing(sorted(self.assembly.reads.keys()),
                            self.assembly_dir,
                            self.project.config.bowtie_indexer_path)
        run_mapper(sorted(self.assembly.reads.keys()),
                   self.assembly_dir,
                   self.project.config.bowtie_path,
                   is_paired_end=self.is_paired_end)
        self.import_contigs()
        self.import_read_mappings()

    def import_contigs(self):
        """Imports assembled contigs from filtered FASTA file"""
        for function in sorted(self.assembly.reads.keys()):
            contig_file = os.path.join(self.assembly_dir,
                                       function,
                                       'final.contigs.filtered.fa')
            if os.path.exists(contig_file):
                with open(contig_file, 'r') as infile:
                    current_id = None
                    sequence = ''
                    for line in infile:
                        line = line.rstrip('\n\r')
                        if line.startswith('>'):
                            if current_id is not None:
                                contig = Contig(contig_id=current_id, sequence=sequence)
                                self.assembly.contigs[function][current_id] = contig
                            if self.assembler == 'megahit':
                                line_tokens = line[1:].split(' ')
                                current_id = line_tokens[0]
                            elif self.assembler == 'metaspades':
                                current_id = line[1:]
                            else:
                                raise ValueError('Unknown assembler: ' + self.assembler)

                            sequence = ''
                        else:
                            sequence += line
                    if current_id is not None:
                        contig = Contig(contig_id=current_id, sequence=sequence)
                        self.assembly.contigs[function][current_id] = contig
            else:
                print('File ' + contig_file + ' does not exist.')

    def import_read_mappings(self):
        """Imports read mapping data from SAM file(s)"""
        for function in sorted(self.assembly.reads.keys()):
            sam_file = os.path.join(self.assembly_dir, function, 'contigs.sam')
            if os.path.exists(sam_file):
                with open(sam_file, 'r') as infile:
                    for line in infile:
                        if line.startswith('@'):
                            continue
                        line_tokens = line.split('\t')
                        if len(line_tokens) > 9:
                            read_id = line_tokens[0]
                            contig_id = line_tokens[2]
                            alignment_length = len(line_tokens[9])
                            if contig_id in self.assembly.contigs[function]:
                                self.assembly.contigs[function][contig_id].update_coverage(
                                    self.assembly.reads[function][read_id],
                                    alignment_length
                                    )
                                self.assembly.contigs[function][contig_id].reads.append(read_id)
            else:
                print('File ' + sam_file + ' does not exist.')

    def filter_contigs_by_length(self):
        """Filters list of contigs by length

        TODO:
            make contig_length_threshold a parameter in ProgramConfig or constant
        """
        contig_length_threshold = 300
        for function in self.assembly.reads.keys():
            contig_file = os.path.join(self.assembly_dir, function, 'final.contigs.fa')
            if not os.path.exists(contig_file):
                continue

            outfile = os.path.join(self.assembly_dir, function, 'final.contigs.filtered.fa')
            with open(outfile, 'w') as outfile:
                with open(contig_file, 'r') as infile:
                    current_id = None
                    sequence = []
                    for line in infile:
                        line = line.rstrip('\n\r')
                        if line.startswith('>'):
                            contig_sequence = ''.join(sequence)
                            if current_id and len(contig_sequence) >= contig_length_threshold:
                                outfile.write('\n'.join([current_id, contig_sequence, '']))
                            line_tokens = line.split(' ')
                            current_id = line_tokens[0]
                            sequence = []
                        else:
                            sequence.append(line)
                    contig_sequence = ''.join(sequence)
                    if len(contig_sequence) >= contig_length_threshold:
                        outfile.write('\n'.join([current_id, contig_sequence, '']))

    def parse_reference_output(self):
        """Reads and processes DIAMOND tabular output of the preselection
        DIAMOND search.

        Note: this function finds query sequences similar to reference
        proteins. Since a query sequence may have more than one areas of
        similarity (for instance, in fusion proteins of two subunits or
        in multi-domain proteins), it will try to find as many such areas
        as possible.

        DIAMOND hits are filtered by two parameters: length of alignment
        and amino acid identity %, which are defined in program config ini.
        """
        tsvfile = os.path.join(self.assembly_dir,
                               'all_contigs_' + self.project.options.ref_output_name)
        current_id = ''
        hit_list = DiamondHitList(current_id)
        identity_cutoff = self.project.config.get_identity_cutoff(
            self.project.options.get_collection())
        length_cutoff = self.project.config.get_length_cutoff(
            self.project.options.get_collection())
        print('Parse reference output: Identity cutoff: ',
              identity_cutoff,
              ', Length cutoff: ',
              length_cutoff)

        with open(tsvfile, 'r', newline='') as infile:
            tsvin = csv.reader(infile, delimiter='\t')
            for row in tsvin:
                hit = DiamondHit()
                hit.create_hit(row)
                # filtering by identity and length
                if hit.identity < identity_cutoff:
                    continue  # skip this line
                if hit.length < length_cutoff:
                    continue  # skip this line

                if hit.query_id != current_id:
                    # filter list for overlapping hits
                    hit_list.filter_list(self.project.config.get_overlap_cutoff(
                        self.project.options.get_collection()))
                    if hit_list.hits_number != 0:
                        # annotate_hits
                        hit_list.annotate_hits(self.project.ref_data)
                        function_id, contig_id, _ = parse_gene_id(current_id)
                        self.assembly.contigs[function_id][contig_id].\
                            genes[current_id].hit_list = hit_list

                    current_id = hit.query_id
                    hit_list = DiamondHitList(current_id)
                hit_list.add_hit(hit)
            hit_list.filter_list(
                self.project.config.get_overlap_cutoff(self.project.options.get_collection()))
            if hit_list.hits_number != 0:
                # annotate_hits
                hit_list.annotate_hits(self.project.ref_data)
                function_id, contig_id, _ = parse_gene_id(current_id)
                self.assembly.contigs[function_id][contig_id].genes[current_id].hit_list = \
                    hit_list

    def export_hit_fasta(self):
        """Exports hit sequences as gzipped FASTA file"""
        outfile = os.path.join(
            self.assembly_dir, 'all_contigs_' + self.project.options.ref_hits_fastq_name
        )

        with open(outfile, 'w') as outfile:
            for function in sorted(self.assembly.contigs.keys()):
                for contig_id in sorted(self.assembly.contigs[function].keys()):
                    for gene_id in self.assembly.contigs[function][contig_id].genes.keys():
                        gene = self.assembly.contigs[function][contig_id].genes[gene_id]
                        if not gene.hit_list:
                            continue
                        for hit in gene.hit_list.hits:
                            start = hit.q_start
                            end = hit.q_end
                            outfile.write('>' + '|'.join([gene_id, str(start), str(end)]) + '\n')
                            start = start - 1
                            try:
                                outfile.write(gene.protein_sequence[start:end] + '\n')
                            except TypeError:
                                print('TypeError occurred while exporting ', gene.gene_id)

    def parse_background_output(self):
        """Reads and processes DIAMOND tabular output of the classification
        DIAMOND search.

        Note: this function takes existing list of hits and compares each
        of them with results of new similarity serach (against classification DB).
        For the comparison, it calls compare_hits_lca function.

        """
        tsvfile = os.path.join(self.assembly_dir,
                               'all_contigs_' + self.project.options.background_output_name)
        current_query_id = None
        hit_list = None
        identity_cutoff = self.project.config.get_identity_cutoff(
            self.project.options.get_collection())
        length_cutoff = self.project.config.get_length_cutoff(
            self.project.options.get_collection())
        biscore_range_cutoff = self.project.config.get_biscore_range_cutoff(
            self.project.options.get_collection())
        print('Identity cutoff: ', identity_cutoff, ', Length cutoff: ', length_cutoff)

        average_coverage = self.assembly.calculate_average_coverage()

        with open(tsvfile, 'r', newline='') as infile:
            tsvin = csv.reader(infile, delimiter='\t')
            function_id = ''
            contig_id = ''
            gene_id = ''
            coverage = ''
            for row in tsvin:
                if current_query_id is None:
                    current_query_id = row[0]
                    hit_list = DiamondHitList(current_query_id)

                hit = DiamondHit()
                hit.create_hit(row)
                # filtering by identity and length
                if hit.identity < identity_cutoff or hit.length < length_cutoff:
                    continue  # skip this line

                if hit.query_id != current_query_id:
                    hit_list.annotate_hits(self.project.ref_data)
                    # compare list of hits from search in background DB with existing
                    # hit from search in reference DB
                    current_query_id_tokens = current_query_id.split('|')
                    function_id = current_query_id_tokens[0]
                    contig_id = '_'.join(current_query_id_tokens[1].split('_')[:-1])
                    gene_id = '|'.join(current_query_id_tokens[:-2])
                    coverage = self.assembly.contigs[function_id][contig_id].get_coverage()
                    try:
                        compare_hits_lca(
                            self.assembly.contigs[function_id][contig_id].genes[gene_id],
                            int(current_query_id_tokens[-2]),  # hit_start
                            int(current_query_id_tokens[-1]),  # hit_end
                            hit_list,
                            biscore_range_cutoff,
                            average_coverage,
                            coverage,
                            self.project.taxonomy_data,
                            self.project.ref_data,
                            rank_cutoffs=self.project.config.get_ranks_cutoffs(
                                self.project.options.get_collection()
                                )
                            )
                    except KeyError:
                        print(' '.join(['Gene not found:', gene_id, 'in', function_id, contig_id]))
                    current_query_id = hit.query_id
                    hit_list = DiamondHitList(current_query_id)
                hit_list.add_hit(hit)
            hit_list.annotate_hits(self.project.ref_data)
            current_query_id_tokens = current_query_id.split('|')
            function_id = current_query_id_tokens[0]
            contig_id = '_'.join(current_query_id_tokens[1].split('_')[:-1])
            gene_id = '|'.join(current_query_id_tokens[:-2])
            coverage = self.assembly.contigs[function_id][contig_id].get_coverage()
            try:
                compare_hits_lca(
                    self.assembly.contigs[function_id][contig_id].genes[gene_id],
                    int(current_query_id_tokens[-2]),  # hit_start
                    int(current_query_id_tokens[-1]),  # hit_end
                    hit_list,
                    biscore_range_cutoff,
                    average_coverage,
                    coverage,
                    self.project.taxonomy_data,
                    self.project.ref_data,
                    rank_cutoffs=self.project.config.get_ranks_cutoffs(
                        self.project.options.get_collection()
                    )
                )
            except KeyError:
                print(' '.join(['Gene not found:', gene_id, 'in', function_id, contig_id]))

    def predict_genes(self):
        """Filters contigs by coverage, runs Prodigal on remaining contigs,

        Todo:
            make contig_coverage_cutoff a parameter or a constant
        """
        # Filter contigs by coverage
        contig_coverage_cutoff = 3.0

        prodigal_infile = os.path.join(self.assembly_dir, 'all_contigs.fa')
        with open(prodigal_infile, 'w') as outfile:
            for function in sorted(self.assembly.contigs.keys()):
                for contig in sorted(self.assembly.contigs[function].keys()):
                    if self.assembly.contigs[function][
                            contig
                    ].get_coverage() >= contig_coverage_cutoff:
                        outfile.write('>' + function + '|' + contig + '\n')
                        outfile.write(self.assembly.contigs[function][contig].sequence + '\n')

        # Run Prodigal
        prodigal_outfile = os.path.join(self.assembly_dir, 'all_contigs.prodigal.out.faa')
        run_prodigal(prodigal_infile, prodigal_outfile, self.project.config.prodigal_path)

        with open(prodigal_outfile, 'r') as infile:
            current_id = None
            sequence = ''
            for line in infile:
                line = line.rstrip('\n\r')
                if line.startswith('>'):
                    if current_id:
                        line_tokens = current_id.split(' # ')
                        function_id, contig_id, _ = parse_gene_id(line_tokens[0])
                        gene = Gene(contig_id=contig_id,
                                    gene_id=line_tokens[0],
                                    sequence=sequence,
                                    start=line_tokens[1],
                                    end=line_tokens[2],
                                    strand=line_tokens[3])
                        self.assembly.contigs[function_id][contig_id].add_gene(gene)
                    line_tokens = line.split(' ')
                    current_id = line[1:]  # line_tokens[0][1:]
                    sequence = ''
                else:
                    sequence += line
            line_tokens = current_id.split(' # ')
            function_id, contig_id, _ = parse_gene_id(line_tokens[0])
            gene = Gene(contig_id=contig_id,
                        gene_id=line_tokens[0],
                        sequence=sequence,
                        start=line_tokens[1],
                        end=line_tokens[2],
                        strand=line_tokens[3])
            self.assembly.contigs[function_id][contig_id].add_gene(gene)

    def annotate_genes(self):
        """Runs pre-selection DIAMOND search, runs classification DIAMOND search,
        exports assembly in JSON format

        Todo:
            make contig_coverage_cutoff a parameter or a constant
        """
        # Search in reference database
        run_ref_search(self.project)

        # Process output of reference DB search
        self.parse_reference_output()
        export_gene_assembly(
            self.assembly, os.path.join(self.assembly_dir, 'all_contigs_assembly.json'))

        # Import sequence data for selected sequence reads
        print('Reading FASTQ file')
        self.export_hit_fasta()

        # Search in background database
        run_bgr_search(self.project)

        # Process output of reference DB search
        self.parse_background_output()

        print('Exporting JSON')
        export_gene_assembly(self.assembly,
                             os.path.join(self.assembly_dir,
                                          'all_contigs_assembly.json'))

    def generate_taxonomy_chart(self, taxonomy_data):
        '''
        Collects data about functional genes in assembly and
        creates one Krona chart for all functions

        Args:
            taxonomy_data (:obj:TaxonomyData): NCBI taxonomy data
        '''
        functions_list = set()
        genes = autovivify(2)  # genes[gene][function][parameter] = parameter_value
        scores = autovivify(2)  # scores[taxonomy ID][function][parameter] = parameter_value

        total_read_count = 0
        for sample in self.project.list_samples():
            total_read_count += self.project.options.get_fastq1_readcount(sample)

        for function in self.assembly.contigs:
            functions_list.add(function)
            for _, contig in self.assembly.contigs[function].items():
                for gene_id, gene in contig.genes.items():
                    if gene.status != STATUS_GOOD:
                        continue
                    taxonomy_id = gene.taxonomy  # Was get_taxonomy_id()
                    for hit in gene.hit_list.hits:
                        identity = hit.identity
                        for hit_function in hit.functions:
                            functions_list.add(hit_function)
                            if 'rpkm' in scores[taxonomy_id][hit_function]:
                                scores[taxonomy_id][hit_function]['rpkm'] += \
                                    contig.get_rpkm(total_read_count) * \
                                    len(gene.protein_sequence) * 3 / len(contig.sequence)
                            else:
                                scores[taxonomy_id][hit_function]['rpkm'] = \
                                    contig.get_rpkm(total_read_count) * \
                                    len(gene.protein_sequence) * 3 / len(contig.sequence)
                            if 'count' in scores[taxonomy_id][hit_function]:
                                scores[taxonomy_id][hit_function]['count'] += \
                                    contig.get_read_count() * len(gene.protein_sequence) * 3 / \
                                    len(contig.sequence)
                            else:
                                scores[taxonomy_id][hit_function]['count'] = \
                                    contig.get_read_count() * len(gene.protein_sequence) * 3 / \
                                    len(contig.sequence)
                            if 'hit_count' in scores[taxonomy_id][hit_function]:
                                scores[taxonomy_id][hit_function]['hit_count'] += 1
                            else:
                                scores[taxonomy_id][hit_function]['hit_count'] = 1
                            if 'identity' in scores[taxonomy_id][hit_function]:
                                scores[taxonomy_id][hit_function]['identity'] += \
                                    identity
                            else:
                                scores[taxonomy_id][hit_function]['identity'] = \
                                    identity
                            if 'genes' in scores[taxonomy_id][hit_function]:
                                scores[taxonomy_id][hit_function]['genes'] += gene_id + ' '
                            else:
                                scores[taxonomy_id][hit_function]['genes'] = gene_id + ' '

                            genes[gene_id][hit_function]['Length'] = \
                                str(len(gene.protein_sequence)) + 'aa'
                            genes[gene_id][hit_function]['Completeness'] = '{0:.0f}'.format(
                                len(gene.protein_sequence) * 100 / hit.s_len
                            )
                            genes[gene_id][hit_function]['identity'] = '{0:.1f}'.format(
                                identity
                            )
                            genes[gene_id][hit_function]['rpkm'] = '{0:.6f}'.format(
                                contig.get_rpkm(
                                    total_read_count
                                ) * len(gene.protein_sequence) * 3 / len(contig.sequence)
                            )
                            genes[gene_id][hit_function]['count'] = '{0:.0f}'.format(
                                contig.get_read_count() * len(
                                    gene.protein_sequence
                                ) * 3 / len(contig.sequence)
                            )
                            genes[gene_id][hit_function]['coverage'] = '{0:.1f}'.format(
                                contig.get_coverage()
                            )

        taxonomic_profile = TaxonomyProfile()
        taxonomic_profile.make_assembly_taxonomy_profile(taxonomy_data, scores)

        outfile = os.path.join(self.assembly_dir, 'assembly_taxonomic_profile.xml')
        make_assembly_taxonomy_chart(
            taxonomic_profile, genes, sorted(functions_list), outfile,
            self.project.config.krona_path, metric='rpkm'
            )

    def generate_function_taxonomy_charts(self, taxonomy_data):
        '''
        Generates series of Krona charts visualizing functions in assembly:
        one function per file, separate stats for each sample

        Args:
            taxonomy_data (:obj:TaxonomyData): NCBI taxonomy data
        '''
        functions_list = set()
        samples_list = sorted(self.project.list_samples())

        total_read_count = 0
        for sample in self.project.list_samples():
            total_read_count += self.project.options.get_fastq1_readcount(sample)

        # Make list of functions
        for function in self.assembly.contigs:
            for contig in self.assembly.contigs[function]:
                for gene_id, gene in self.assembly.contigs[function][contig].genes.items():
                    if gene.status == STATUS_GOOD:
                        for gene_function in gene.functions:
                            functions_list.add(gene_function)

        for function in sorted(functions_list):
            genes = autovivify(2)  # genes[gene][sample][parameter] = parameter_value
            scores = autovivify(2)  # scores[taxonomy ID][sample][parameter] = parameter_value
            outfile = os.path.join(self.assembly_dir, 'out', function + '_taxonomic_profile.xml')
            for assembly_function in self.assembly.contigs:
                for _, contig in self.assembly.contigs[assembly_function].items():
                    for gene_id, gene in contig.genes.items():
                        function_counted = False
                        if gene.status != STATUS_GOOD or function not in gene.functions:
                            continue
                        taxonomy_id = gene.taxonomy
                        if taxonomy_id not in scores:
                            for sample_id in samples_list:
                                scores[taxonomy_id][sample_id]['rpkm'] = 0.0
                                scores[taxonomy_id][sample_id]['count'] = 0
                                scores[taxonomy_id][sample_id]['hit_count'] = 0
                                scores[taxonomy_id][sample_id]['identity'] = 0.0
                                scores[taxonomy_id][sample_id]['genes'] = ''
                            scores[taxonomy_id]['All samples']['rpkm'] = 0.0
                            scores[taxonomy_id]['All samples']['count'] = 0
                            scores[taxonomy_id]['All samples']['hit_count'] = 0
                            scores[taxonomy_id]['All samples']['identity'] = 0.0
                            scores[taxonomy_id]['All samples']['genes'] = ''

                        for hit in gene.hit_list.hits:
                            identity = hit.identity
                            if function in hit.functions:
                                if function_counted:
                                    continue
                                for sample in samples_list:
                                    if sample in contig.read_count:
                                        scores[taxonomy_id][sample]['rpkm'] += contig.get_rpkm(
                                            self.project.options.get_fastq1_readcount(sample),
                                            sample
                                            ) * len(gene.protein_sequence) * 3 / len(
                                                contig.sequence
                                            )
                                        scores[taxonomy_id][sample]['count'] += \
                                            contig.get_read_count(sample) * \
                                            len(gene.protein_sequence) * 3 / \
                                            len(contig.sequence)
                                        scores[taxonomy_id][sample]['hit_count'] += 1
                                        scores[taxonomy_id][sample]['identity'] += identity
                                        scores[taxonomy_id][sample]['genes'] += gene_id + ' '

                                        genes[gene_id][sample]['Length'] = \
                                            str(len(gene.protein_sequence)) + 'aa'
                                        genes[gene_id][sample]['Completeness'] = '{0:.0f}'.format(
                                            len(gene.protein_sequence) * 100 / hit.s_len
                                        )
                                        genes[gene_id][sample]['identity'] = '{0:.1f}'.format(
                                            identity
                                        )
                                        genes[gene_id][sample]['rpkm'] = '{0:.7f}'.format(
                                            contig.get_rpkm(
                                                self.project.options.get_fastq1_readcount(
                                                    sample
                                                ),
                                                sample
                                            ) * len(gene.protein_sequence) * 3 / len(
                                                contig.sequence
                                            )
                                        )
                                        genes[gene_id][sample]['count'] = 3 * '{0:.0f}'.format(
                                            contig.get_read_count(sample) * len(
                                                gene.protein_sequence
                                            ) / len(contig.sequence)
                                        )
                                        genes[gene_id][sample]['coverage'] = '{0:.1f}'.format(
                                            contig.get_coverage(sample)
                                        )
                                scores[taxonomy_id]['All samples']['rpkm'] += \
                                    contig.get_rpkm(total_read_count) * \
                                    len(gene.protein_sequence) \
                                    * 3 / len(contig.sequence)
                                scores[taxonomy_id]['All samples']['count'] += \
                                    contig.get_read_count() * len(gene.protein_sequence) \
                                    * 3 / len(contig.sequence)
                                scores[taxonomy_id]['All samples']['hit_count'] += 1
                                scores[taxonomy_id]['All samples']['identity'] += identity
                                scores[taxonomy_id]['All samples']['genes'] += gene_id + ' '
                                genes[gene_id]['All samples']['Length'] = \
                                    str(len(gene.protein_sequence)) + 'aa'
                                genes[gene_id]['All samples']['Completeness'] = '{0:.0f}'.format(
                                    len(gene.protein_sequence) * 100 / hit.s_len
                                )
                                genes[gene_id]['All samples']['identity'] = \
                                    '{0:.1f}'.format(identity)
                                genes[gene_id]['All samples']['rpkm'] = '{0:.7f}'.format(
                                    contig.get_rpkm(total_read_count) * len(
                                        gene.protein_sequence
                                    ) * 3 / len(contig.sequence)
                                )
                                genes[gene_id]['All samples']['count'] = '{0:.0f}'.format(
                                    contig.get_read_count() * len(
                                        gene.protein_sequence
                                    ) * 3 / len(contig.sequence)
                                )
                                genes[gene_id]['All samples']['coverage'] = '{0:.1f}'.format(
                                    contig.get_coverage()
                                )
                                function_counted = True
            taxonomic_profile = TaxonomyProfile()
            taxonomic_profile.make_assembly_taxonomy_profile(taxonomy_data, scores)

            output_sample_ids = sorted(self.project.list_samples())
            output_sample_ids.append('All samples')
            make_assembly_taxonomy_chart(
                taxonomic_profile, genes, output_sample_ids, outfile,
                self.project.config.krona_path, metric='rpkm'
                )

    def write_sequences(self):
        """Exports gene and protein sequences in FASTA format"""
        genes = autovivify(2)  # genes[function][gene][parameter] = parameter_value

        for function in self.assembly.contigs:
            for contig in self.assembly.contigs[function]:
                for gene_id, gene in self.assembly.contigs[function][contig].genes.items():
                    if gene.status == STATUS_GOOD:
                        for hit in gene.hit_list.data:
                            taxonomy_id = gene.taxonomy
                            for hit_function in hit.functions:
                                start = gene.start
                                end = gene.end
                                strand = gene.strand
                                genes[hit_function][gene_id]['start'] = start
                                genes[hit_function][gene_id]['end'] = end
                                genes[hit_function][gene_id]['strand'] = strand
                                genes[hit_function][gene_id]['taxonomy'] = taxonomy_id
                                gene_sequence = self.assembly.contigs[function][contig].\
                                    sequence[int(start) - 1: int(end)]
                                if strand == '-1':
                                    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
                                    gene_sequence = ''.join(
                                        [complement[nucl] for nucl in reversed(gene_sequence)]
                                    )
                                genes[hit_function][gene_id]['sequence'] = gene_sequence
                                genes[hit_function][gene_id]['protein'] = gene.protein_sequence
                                genes[hit_function][gene_id]['aai'] = hit.identity
                                genes[hit_function][gene_id]['completeness'] = \
                                    len(gene.protein_sequence) * 100 / hit.s_len
        for function in genes:
            outfile = os.path.join(self.project.options.assembly_dir,
                                   'out',
                                   function + '_genes_Fama.fna')
            with open(outfile, 'w') as outfile:
                for gene_id in genes[function]:
                    lineage = self.project.taxonomy_data.get_taxonomy_lineage(
                        genes[function][gene_id]['taxonomy'])
                    outfile.write('>' + gene_id + '|' +
                                  genes[function][gene_id]['start'] + '|' +
                                  genes[function][gene_id]['end'] + '|' +
                                  genes[function][gene_id]['strand'] + '|' +
                                  lineage + '\n')  # '|'
                    outfile.write(genes[function][gene_id]['sequence'] + '\n')
            outfile = os.path.join(self.project.options.assembly_dir,
                                   'out',
                                   function + '_proteins_Fama.faa')
            with open(outfile, 'w') as outfile:
                for gene_id in genes[function]:
                    lineage = self.project.taxonomy_data.get_taxonomy_lineage(
                        genes[function][gene_id]['taxonomy'])
                    outfile.write('>' + gene_id + '|' +
                                  genes[function][gene_id]['start'] + '|' +
                                  genes[function][gene_id]['end'] + '|' +
                                  genes[function][gene_id]['strand'] + '|' +
                                  lineage + '\n')  # '|'
                    outfile.write(genes[function][gene_id]['protein'] + '\n')

    def generate_output(self):
        """Sends assembly data to Excel report generator, exports genes and
        proteins, calls methods for taxonomy chart generation
        """
        self.write_sequences()
        make_assembly_xlsx(self)
        self.generate_taxonomy_chart(self.project.taxonomy_data)
        self.generate_function_taxonomy_charts(self.project.taxonomy_data)
        generate_assembly_report(self)


def run_assembler(functions, assembler, output_dir, is_paired_end=True):
    """Fires up external assembler, either metaSPAdes or MEGAHIT"""
    if assembler.endswith('megahit'):
        run_megahit(functions, output_dir, assembler, is_paired_end)
    elif assembler.endswith('metaspades.py'):
        if is_paired_end:
            run_spades(functions, output_dir, assembler, is_paired_end)
        else:
            raise RuntimeError(
                'Current version of metaSPAdes does not support single-end libraries.'
                )


def run_megahit(functions, output_dir, assembler_command, is_paired_end=True):
    """Runs MEGAHIT assembler on exported reads"""
    print('Starting assembly')
    for function in functions:
        print('Run assembler for function', function)
        if is_paired_end:
            assembler_args = [assembler_command,
                              '-1',
                              os.path.join(output_dir, function + '_pe1.fastq'),
                              '-2',
                              os.path.join(output_dir, function + '_pe2.fastq'),
                              '-o',
                              os.path.join(output_dir, function)]
        else:
            assembler_args = [assembler_command,
                              '-r',
                              os.path.join(output_dir, function + '_pe1.fastq'),
                              '-o',
                              os.path.join(output_dir, function)]
        run_external_program(assembler_args)
        print('Assembler finished for function ', function)
    print('Assembly finished')


def run_spades(functions, output_dir, assembler_command, is_paired_end=True):
    """Runs metaSPAdes assembler on exported reads"""
    print('Starting metaSPAdes')
    tmp_dir = os.path.join(output_dir, 'tmp')
    for function in functions:
        print('Run metaSPAdes for function', function)
        assembler_args = [assembler_command,
                          '--meta',
                          '-t',
                          '12',
                          '-m',
                          '50',  # TODO: make a parameter
                          '-k',
                          '33,55,71,91,111',  # TODO: adjust automatically
                          '-o',
                          os.path.join(output_dir, function),
                          '--tmp-dir',
                          tmp_dir]
        if is_paired_end:
            assembler_args.extend(['-1',
                                   os.path.join(output_dir, function + '_pe1.fastq'),
                                   '-2',
                                   os.path.join(output_dir, function + '_pe2.fastq')])
        else:
            assembler_args.extend(['-s',
                                   os.path.join(output_dir, function + '_pe1.fastq')])

        run_external_program_ignoreerror(assembler_args)
        if os.path.exists(os.path.join(output_dir, function, 'contigs.fasta')):
            shutil.copyfile(os.path.join(output_dir, function, 'contigs.fasta'),
                            os.path.join(output_dir, function, 'final.contigs.fa'))
        print('Assembler finished for function ', function)
    print('metaSPAdes finished')


def run_mapper_indexing(functions, output_dir, mapper_command):
    """Runs Bowtie2 indexer on filtered contigs"""
    mapper_command = 'bowtie2-build'

    for function in functions:
        if not os.path.exists(os.path.join(output_dir, function, 'final.contigs.filtered.fa')):
            print('Contigs file for function', function, 'not found')
            continue
        print('Run indexing for function', function)
        if os.path.getsize(os.path.join(output_dir, function, 'final.contigs.filtered.fa')) > 0:
            if not os.path.exists(os.path.join(output_dir, function, 'index')):
                os.mkdir(os.path.join(output_dir, function, 'index'))
            mapper_args = [mapper_command,
                           '-f',
                           os.path.join(output_dir, function, 'final.contigs.filtered.fa'),
                           os.path.join(output_dir, function, 'index', 'index')]
            run_external_program(mapper_args)


def run_mapper(functions, output_dir, mapper_command, is_paired_end=True):
    """Runs Bowtie2 mapper on filtered contigs"""
    mapper_command = 'bowtie2'

    for function in functions:
        if not os.path.exists(os.path.join(output_dir, function, 'final.contigs.filtered.fa')):
            continue
        if os.path.getsize(os.path.join(output_dir, function, 'final.contigs.filtered.fa')) > 0:
            print('Run read mapping for function', function)
            if is_paired_end:
                mapper_args = [mapper_command,
                               '-q',
                               '--very-sensitive',
                               '--quiet',
                               '-x',
                               os.path.join(output_dir, function, 'index', 'index'),
                               '-1',
                               os.path.join(output_dir, function + '_pe1.fastq'),
                               '-2',
                               os.path.join(output_dir, function + '_pe2.fastq'),
                               '>' + os.path.join(output_dir, function, 'contigs.sam')]
            else:
                mapper_args = [mapper_command,
                               '-q',
                               '--very-sensitive',
                               '--quiet',
                               '-x',
                               os.path.join(output_dir, function, 'index', 'index'),
                               '-U',
                               os.path.join(output_dir, function + '_pe1.fastq'),
                               '>' + os.path.join(output_dir, function, 'contigs.sam')]
            run_external_program(mapper_args)


def run_prodigal(infile, outfile, prodigal_path):
    """Runs Prodigal gene prediction on filtered contigs"""
    print('Starting Prodigal')
    prodigal_args = [prodigal_path,
                     '-p',
                     'meta',
                     '-a',
                     outfile,
                     '-i',
                     infile,
                     '-o',
                     outfile+'prodigal.txt']
    run_external_program(prodigal_args)
    print('Prodigal finished')


def run_ref_search(project):
    """Runs DIAMOND pre-selection search on predicted genes"""
    print('Starting DIAMOND')
    diamond_args = [project.config.diamond_path,
                    'blastp',
                    '--db',
                    project.config.get_reference_diamond_db(project.options.get_collection()),
                    '--query',
                    os.path.join(project.options.assembly_dir,
                                 'all_contigs.prodigal.out.faa'),
                    '--out',
                    os.path.join(project.options.assembly_dir,
                                 'all_contigs_' + project.options.ref_output_name),
                    '--max-target-seqs',
                    '50',
                    '--evalue',
                    str(project.config.get_evalue_cutoff(project.options.get_collection())),
                    '--threads',
                    project.config.threads,
                    '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length',
                    'mismatch', 'slen', 'qstart', 'qend', 'sstart', 'send',
                    'evalue', 'bitscore']
    run_external_program(diamond_args)
    print('DIAMOND finished')


def run_bgr_search(project):
    """Runs DIAMOND classification search on predicted genes"""
    print('Starting DIAMOND')
    diamond_args = [project.config.diamond_path,
                    'blastp',
                    '--db',
                    project.config.get_background_diamond_db(project.options.get_collection()),
                    '--query',
                    os.path.join(
                        project.options.assembly_dir,
                        'all_contigs_' + project.options.ref_hits_fastq_name
                    ),
                    '--out',
                    os.path.join(
                        project.options.assembly_dir,
                        'all_contigs_' + project.options.background_output_name
                    ),
                    '--max-target-seqs',
                    '50',
                    '--evalue',
                    str(project.config.get_background_db_size(project.options.get_collection())
                        * project.config.get_evalue_cutoff(project.options.get_collection())
                        / project.config.get_reference_db_size(project.options.get_collection())),
                    '--threads',
                    project.config.threads,
                    '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length',
                    'mismatch', 'slen', 'qstart', 'qend', 'sstart', 'send',
                    'evalue', 'bitscore']
    run_external_program(diamond_args)
    print('DIAMOND finished')


def parse_gene_id(gene_id):
    """Extracts contig identifier and function identifier from gene identifier"""
    (function_id, gene) = gene_id.split('|')
    gene_id_tokens = gene.split('_')
    gene_id = gene_id_tokens[-1]
    contig_id = '_'.join(gene_id_tokens[:-1])
    return function_id, contig_id, gene_id


def get_abundance(function_fraction, average_coverage, coverage):
    """Calculates  relative abundance from contig coverage"""
    if function_fraction > 1.0:
        print('FUNCTION FRACTION TOO BIG!', function_fraction)
    ret_val = coverage * function_fraction/average_coverage
    return ret_val


def compare_hits_lca(gene, hit_start, hit_end, new_hit_list, bitscore_range_cutoff,
                     average_coverage, coverage, taxonomy_data, ref_data, rank_cutoffs=None):
    """Compares DiamondHit object assigned to a Gene object with list of new
    DiamondHit objects, assigns abundance to functions and taxonomy

    Args:
        gene (:obj:Gene): gene being under analysis
        hit_start (int): start position of known hit
        hit_end (int): end position of known hit
        new_hit_list (:obj:'DiamondHitList'): list of hit from new search
        bitscore_range_cutoff (float): lowest acceptaple bitscore
            (relative to top bit-score, default 0.97)
        average_coverage (float): average read coverage across all contigs
        coverage (float): read coverage of the contig
        taxonomy_data (:obj:TaxonomyData): taxonomic data
        ref_data (:obj:ReferenceData): functional reference data
        rank_cutoffs (:obj:dict[str, float]): key is taxonomy rank, value
            is of amino acid identity % threshold for this rank

    This function compares one hit assigned to an annotated read with a list
    of new hits. It looks through the hit list, finds hits with bitscore
    above cutoff and assigns their functions to the read. If any hits to
    functional protein of interest are found, the read gets status 'function'.
    Otherwise, it gets status 'nofunction'.

    This function does not return anything. It sets status of gene and
    assigns abundance to each function associated with the gene.

    """
    if rank_cutoffs is None:
        rank_cutoffs = {}
    # Find best hit
    for hit in gene.hit_list.hits:
        if hit.q_start == hit_start and hit.q_end == hit_end:
            best_bitscore = 0.0
            best_hit = None
            for new_hit in new_hit_list.hits:
                if new_hit.bitscore > best_bitscore:
                    best_hit = new_hit
                    best_bitscore = new_hit.bitscore
            # Set status of read
            if best_hit is not None:
                if '' in best_hit.functions:
                    gene.set_status(STATUS_BAD)
                    return
                else:
                    gene.set_status(STATUS_GOOD)
            else:
                gene.set_status(STATUS_BAD)
                return

            # Filter list of hits by bitscore
            bitscore_lower_cutoff = best_bitscore * (1.0 - bitscore_range_cutoff)
            new_hits = [
                new_hit for new_hit in new_hit_list.hits if new_hit.bitscore > bitscore_lower_cutoff
                ]

            if hit.subject_id not in [
                    new_hit.subject_id for new_hit in new_hits
            ] and hit.bitscore >= best_bitscore:
                new_hits.append(hit)

            # Collect taxonomy IDs of all hits for LCA inference
            taxonomy_ids = set()
            # If rank-specific AAI cutoffs are not set
            if not rank_cutoffs:
                taxonomy_ids = set([ref_data.lookup_protein_tax(h.subject_id) for h in new_hits])

            # If rank-specific AAI cutoffs were calculated for the reference dataset:
            else:
                for new_hit in new_hits:
                    subject_taxon_id = ref_data.lookup_protein_tax(new_hit.subject_id)
                    subject_rank = taxonomy_data.get_rank(subject_taxon_id)
                    while subject_taxon_id != ROOT_TAXONOMY_ID:
                        if subject_rank not in rank_cutoffs:
                            (subject_taxon_id, subject_rank) = \
                                taxonomy_data.get_upper_level_taxon(subject_taxon_id)
                        elif new_hit.identity < rank_cutoffs[subject_rank]:
                            (subject_taxon_id, subject_rank) = \
                                taxonomy_data.get_upper_level_taxon(subject_taxon_id)
                        else:
                            taxonomy_ids.add(subject_taxon_id)
                            break
            

            # Make non-redundant list of functions from hits after filtering
            new_functions = {}
            new_functions_counter = Counter()
            new_functions_dict = defaultdict(dict)
            # Find best hit for each function: only one hit with highest
            # bitscore to be reported for each function
            for new_hit in new_hits:
                for hit_func in new_hit.functions:
                    new_functions_counter[hit_func] += 1
                    if hit_func in new_functions_dict:
                        if new_hit.bitscore > new_functions_dict[hit_func]['bit_score']:
                            new_functions_dict[hit_func]['bit_score'] = new_hit.bitscore
                            new_functions_dict[hit_func]['hit'] = new_hit
                    else:
                        new_functions_dict[hit_func]['bit_score'] = new_hit.bitscore
                        new_functions_dict[hit_func]['hit'] = new_hit

            # If the most common function in new hits is unknown, set status STATUS_BAD and return
            if new_functions_counter.most_common(1)[0][0] == '':
                gene.set_status(STATUS_BAD)
                return

            # Calculate RPK scores for functions
            for function in new_functions_dict:
                if function == '':
                    continue
                new_functions[function] = get_abundance(1.0, average_coverage, coverage)

            gene.set_functions(new_functions)

            # Set new list of hits
            _hit_list = DiamondHitList(gene.gene_id)
            for new_function in new_functions_dict:
                if new_function == '':
                    continue
                good_hit = new_functions_dict[new_function]['hit']
                good_hit.query_id = gene.gene_id
                good_hit.annotate_hit(ref_data)
                _hit_list.add_hit(good_hit)
            gene.hit_list = _hit_list
            # Set read taxonomy ID
            gene.taxonomy = taxonomy_data.get_lca(taxonomy_ids)


def cleanup_read_id(read_id):
    """Removes end identifier from the end of read identifier"""
    result = read_id
    if read_id.endswith('.1'):
        result = read_id[:-2]
    elif read_id.endswith('.2'):
        result = read_id[:-2]
    return result
