import os, csv, operator, shutil
from subprocess import Popen, PIPE, CalledProcessError
from collections import Counter, defaultdict

from Fama.utils import autovivify,cleanup_protein_id
from Fama.GeneAssembler.Contig import Contig
from Fama.GeneAssembler.Gene import Gene
from Fama.GeneAssembler.GeneAssembly import GeneAssembly
from Fama.DiamondParser.DiamondHitList import DiamondHitList
from Fama.DiamondParser.DiamondHit import DiamondHit
from Fama.OutputUtil.JSONUtil import export_gene_assembly
from Fama.TaxonomyProfile import TaxonomyProfile
from Fama.OutputUtil.KronaXMLWriter import generate_assembly_taxonomy_chart
from Fama.OutputUtil.Report import generate_assembly_report
from Fama.OutputUtil.XlsxUtil import create_assembly_xlsx
from Fama.ReferenceLibrary.TaxonomyData import TaxonomyData
from Fama.ReferenceLibrary.UniprotData import UniprotData

from Fama.OutputUtil.JSONUtil import import_gene_assembly

class GeneAssembler:
    def __init__(self, project, assembler = 'metaspades'):
        self.project = project
        self.assembler = assembler
        project.load_project()
        self.assembly = GeneAssembly()
        self.is_paired_end = None
        self.assembly_dir = os.path.join(self.project.options.get_assembly_dir())
        if not os.path.isdir(self.assembly_dir):
            os.mkdir(self.assembly_dir)
        if not os.path.isdir(os.path.join(self.assembly_dir,'out')):
            os.mkdir(os.path.join(self.assembly_dir,'out'))
        self.uniprot = UniprotData(self.project.config)


    def assemble_contigs(self):
        # Export reads in FASTQ format
        for f in os.listdir(self.assembly_dir):
            if f.endswith('.fastq'):
                os.remove(os.path.join(self.assembly_dir, f))        

        for sample_id in sorted(self.project.list_samples()):
            self.is_paired_end = self.project.samples[sample_id].is_paired_end
            self.project.import_reads_json(sample_id, self.project.ENDS)
            for end in self.project.ENDS:
                #print ('Loading mapped reads: ', sample, end)
                #self.project.load_annotated_reads(sample, end) # Lazy load
                for read_id in self.project.samples[sample_id].reads[end]:
                    read = self.project.samples[sample_id].reads[end][read_id]
#                    print(read_id)
                    if read.get_status() == 'function,besthit' or read.get_status() == 'function':
                        for function in read.get_functions():
                            if read_id in self.assembly.reads[function]:
                                continue
                            self.assembly.reads[function][read_id] = sample_id
                            outfile1 = os.path.join(self.assembly_dir,function + '_pe1.fastq')
                            if self.is_paired_end:
                                outfile2 = os.path.join(self.assembly_dir,function + '_pe2.fastq')
                                if end == 'pe2':
                                    outfile1 = os.path.join(self.assembly_dir,function + '_pe2.fastq')
                                    outfile2 = os.path.join(self.assembly_dir,function + '_pe1.fastq')
                            
                            with open(outfile1, 'a') as of1:
                                of1.write(read.read_id_line + '\n')
                                of1.write(read.sequence + '\n')
                                of1.write(read.line3 + '\n')
                                of1.write(read.quality + '\n')
                                of1.closed
                            if self.is_paired_end:
                                if read.pe_id is None:
                                    print(read.read_id_line)
                                with open(outfile2, 'a') as of2:
                                    of2.write(read.pe_id + '\n')
                                    of2.write(read.pe_sequence + '\n')
                                    of2.write(read.pe_line3 + '\n')
                                    of2.write(read.pe_quality + '\n')
                                    of2.closed
                # Delete reads from memory
                self.project.samples[sample_id].reads[end] = None
        
        # Run Assembler ('megahit' for Megahit or 'metaSPAdes' for metaSPAdes)
        if self.assembler == 'megahit':
            run_assembler(sorted(self.assembly.reads.keys()), self.project.config.get_megahit_path(), self.assembly_dir)
        elif self.assembler == 'metaspades':
            run_assembler(sorted(self.assembly.reads.keys()), self.project.config.get_metaspades_path(), self.assembly_dir, is_paired_end = self.is_paired_end)
        else:
            raise ValueError('Unknown assembler: ' + self.assembler)
        
        # Filter contigs by size
        self.filter_contigs_by_length()

        # Run Bowtie
        
        run_mapper_indexing(sorted(self.assembly.reads.keys()), self.assembly_dir, self.project.config.get_bowtie_indexer_path())
        run_mapper(sorted(self.assembly.reads.keys()), self.assembly_dir, self.project.config.get_bowtie_path(), is_paired_end = self.is_paired_end)

        # Import contig sequences
        for function in sorted(self.assembly.reads.keys()):
            contig_file = os.path.join(self.assembly_dir,function,'final.contigs.filtered.fa')
            if os.path.exists(contig_file):
                with open (contig_file, 'r') as f:
                    current_id = None
                    sequence = ''
                    for line in f:
                        line = line.rstrip('\n\r')
                        if line.startswith('>'):
                            if current_id:
                                contig = Contig(contig_id=current_id,sequence=sequence)
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
                    if current_id:
                        contig = Contig(contig_id=current_id,sequence=sequence)
                        self.assembly.contigs[function][current_id] = contig
            else:
                print('File ' + contig_file + ' does not exist.')

        export_gene_assembly(self.assembly, os.path.join(self.assembly_dir, 'all_contigs_assembly.json'))

        # Import contig mapping data
        for function in sorted(self.assembly.reads.keys()):
            sam_file = os.path.join(self.assembly_dir,function,'contigs.sam')
            if os.path.exists(sam_file):
                with open (sam_file, 'r') as f:
                    for line in f:
                        if line.startswith('@'):
                            continue
                        line_tokens = line.split('\t')
                        if len(line_tokens) > 9:
                            read_id = line_tokens[0]
                            contig_id = line_tokens[2]
                            alignment_length = len(line_tokens[9])
                            if contig_id in self.assembly.contigs[function]:
                                self.assembly.contigs[function][contig_id].update_coverage(self.assembly.reads[function][read_id],alignment_length)
                                self.assembly.contigs[function][contig_id].reads.append(read_id)
                    f.closed
            else:
                print('File ' + sam_file + ' does not exist.')

    def coassemble_contigs(self):
        for f in os.listdir(self.assembly_dir):
            if f.endswith('.fastq'):
                os.remove(os.path.join(self.assembly_dir, f))        
        # Export reads in FASTQ format
        for sample_id in sorted(self.project.list_samples()):
            self.is_paired_end = self.project.samples[sample_id].is_paired_end
            self.project.import_reads_json(sample_id, self.project.ENDS)
            for end in self.project.ENDS:
#                print ('Loading mapped reads: ', sample, end)
#                self.project.load_annotated_reads(sample, end) # Lazy load
                for read_id in self.project.samples[sample_id].reads[end]:
                    read = self.project.samples[sample_id].reads[end][read_id]
                    if read.get_status() == 'function,besthit' or read.get_status() == 'function':
                        if read_id in self.assembly.reads['Coassembly']:
                            continue
                        self.assembly.reads['Coassembly'][read_id] = sample_id
                        outfile1 = os.path.join(self.assembly_dir,'Coassembly_pe1.fastq')
                        if self.is_paired_end:
                            outfile2 = os.path.join(self.assembly_dir,'Coassembly_pe2.fastq')
                            if end == 'pe2':
                                outfile1 = os.path.join(self.assembly_dir,'Coassembly_pe2.fastq')
                                outfile2 = os.path.join(self.assembly_dir,'Coassembly_pe1.fastq')
                        
                        with open(outfile1, 'a') as of1:
                            of1.write(read.read_id_line + '\n')
                            of1.write(read.sequence + '\n')
                            of1.write(read.line3 + '\n')
                            of1.write(read.quality + '\n')
                            of1.closed
                        if self.is_paired_end:
                            with open(outfile2, 'a') as of2:
                                of2.write(read.pe_id + '\n')
                                of2.write(read.pe_sequence + '\n')
                                of2.write(read.pe_line3 + '\n')
                                of2.write(read.pe_quality + '\n')
                                of2.closed
                # Delete reads from memory
                self.project.samples[sample_id].reads[end] = None
        
        # Run Assembler ('megahit' for Megahit or 'metaSPAdes' for metaSPAdes)
        if self.assembler == 'megahit':
            run_assembler(sorted(self.assembly.reads.keys()), self.project.config.get_megahit_path(), self.assembly_dir, is_paired_end = self.is_paired_end)
        elif self.assembler == 'metaspades':
            run_assembler(sorted(self.assembly.reads.keys()), self.project.config.get_metaspades_path(), self.assembly_dir, is_paired_end = self.is_paired_end)
        else:
            raise ValueError('Unknown assembler: ' + self.assembler)

        # Filter contigs by size
        self.filter_contigs_by_length()

        # Run Bowtie
        contig_file = os.path.join(self.assembly_dir,'Coassembly','final.contigs.filtered.fa')
        print('Run read mapping')
        run_mapper_indexing(sorted(self.assembly.reads.keys()), self.assembly_dir, self.project.config.get_bowtie_indexer_path())
        run_mapper(sorted(self.assembly.reads.keys()), self.assembly_dir, self.project.config.get_bowtie_path(), is_paired_end = self.is_paired_end)

        # Import contig sequences
        if os.path.exists(contig_file):
            with open (contig_file, 'r') as f:
                current_id = None
                sequence = ''
                for line in f:
                    line = line.rstrip('\n\r')
                    if line.startswith('>'):
                        if current_id:
                            contig = Contig(contig_id=current_id,sequence=sequence)
                            self.assembly.contigs['Coassembly'][current_id] = contig
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
                if current_id:
                    contig = Contig(contig_id=current_id,sequence=sequence)
                    self.assembly.contigs['Coassembly'][current_id] = contig
        else:
            print('File ' + contig_file + ' does not exist.')
        
        # Import contig mapping data
        sam_file = os.path.join(self.assembly_dir,'Coassembly','contigs.sam')
        if os.path.exists(sam_file):
            with open (sam_file, 'r') as f:
                for line in f:
                    if line.startswith('@'):
                        continue
                    line_tokens = line.split('\t')
                    if len(line_tokens) > 9:
                        read_id = cleanup_read_id(line_tokens[0])
                        contig_id = line_tokens[2]
                        alignment_length = len(line_tokens[9])
                        if contig_id in self.assembly.contigs['Coassembly']:
                            self.assembly.contigs['Coassembly'][contig_id].update_coverage(self.assembly.reads['Coassembly'][read_id],alignment_length)
                            self.assembly.contigs['Coassembly'][contig_id].reads.append(read_id)
                f.closed
        else:
            print('File ' + sam_file + ' does not exist.')

    def filter_contigs_by_length(self):
        
        contig_length_threshold = 300
        
        for function in self.assembly.reads.keys():
            contig_file = os.path.join(self.assembly_dir,function,'final.contigs.fa')
            if not os.path.exists(contig_file):
                continue

            outfile = os.path.join(self.assembly_dir,function,'final.contigs.filtered.fa')
            with open (outfile, 'w') as of:
                with open (contig_file, 'r') as f:
                    current_id = None
                    sequence = ''
                    for line in f:
                        line = line.rstrip('\n\r')
                        if line.startswith('>'):
                            if current_id:
                                if len(sequence) >= contig_length_threshold:
                                    of.write(current_id + '\n')
                                    of.write(sequence + '\n')
                            line_tokens = line.split(' ')
                            current_id = line_tokens[0]
                            sequence = ''
                        else:
                            sequence += line
                    if len(sequence) >= contig_length_threshold:
                        of.write(current_id + '\n')
                        of.write(sequence + '\n')
                of.closed

    def parse_reference_output(self):
        tsvfile = os.path.join(self.assembly_dir, 'all_contigs_' + self.project.options.get_ref_output_name())
        
        current_id = ''
        _hit_list = DiamondHitList(current_id)
        identity_cutoff = self.project.config.get_identity_cutoff(self.project.options.get_collection())
        length_cutoff = self.project.config.get_length_cutoff(self.project.options.get_collection())
        print ('Parse reference output: Identity cutoff: ', identity_cutoff, ', Length cutoff: ', length_cutoff)
        
        with open(tsvfile, 'r', newline='') as f:
            tsvin = csv.reader(f, delimiter='\t')
            for row in tsvin:
                hit = DiamondHit()
                hit.create_hit(row)
                # filtering by identity and length
                if hit.get_identity() < identity_cutoff:
                    continue # skip this line
                if hit.get_length() < length_cutoff:
                    continue # skip this line

                if hit.get_query_id() != current_id:
                    # filter list for overlapping hits
                    _hit_list.filter_list(self.project.config.get_overlap_cutoff(self.project.options.get_collection()))
                    if _hit_list.get_hits_number() != 0:
                        # annotate_hits
                        _hit_list.annotate_hits(self.project.ref_data)
                        function_id,contig_id,gene_id = parse_gene_id(current_id)
                        #print(current_id,function_id,contig_id,gene_id)
                        #print('Genes:',self.assembly.contigs[function_id][contig_id].genes.keys())
                        self.assembly.contigs[function_id][contig_id].genes[current_id].set_hit_list(_hit_list)

                    current_id = hit.get_query_id()
                    _hit_list = DiamondHitList(current_id)
                _hit_list.add_hit(hit)
            if _hit_list.get_hits_number() != 0:
                _hit_list.filter_list(self.project.config.get_overlap_cutoff(self.project.options.get_collection()))
                # annotate_hits
                _hit_list.annotate_hits(self.project.ref_data)
                function_id,contig_id,gene_id = parse_gene_id(current_id)
                self.assembly.contigs[function_id][contig_id].genes[current_id].set_hit_list(_hit_list)

    def export_hit_fasta(self):
        outfile = os.path.join(self.assembly_dir, 'all_contigs_'+ self.project.options.get_ref_hits_fastq_name())
        
        with open(outfile, 'w') as of:
            for function in sorted(self.assembly.contigs.keys()):
                for contig_id in sorted(self.assembly.contigs[function].keys()):
                    for gene_id in self.assembly.contigs[function][contig_id].genes.keys():
                        gene = self.assembly.contigs[function][contig_id].genes[gene_id]
                        if not gene.hit_list:
                            continue
                        for hit in gene.hit_list.get_hits():
                            start = hit.get_query_start()
                            end = hit.get_query_end()
                            of.write(">" + gene_id + '|' + str(start) + '|' + str(end) + '\n')
                            start = start - 1
                            try:
                                of.write(gene.protein_sequence[start:end] + '\n') 
                            except TypeError:
                                print ('TypeError occurred while exporting ', gene.gene_id)
            of.closed

    def parse_background_output(self):
        
        tsvfile = os.path.join(self.assembly_dir, 'all_contigs_'+ self.project.options.get_background_output_name())
        
        current_query_id = None
        _hit_list = None
        identity_cutoff = self.project.config.get_identity_cutoff(self.project.options.get_collection())
        length_cutoff = self.project.config.get_length_cutoff(self.project.options.get_collection())
        biscore_range_cutoff = self.project.config.get_biscore_range_cutoff(self.project.options.get_collection())
        print ('Identity cutoff: ', identity_cutoff, ', Length cutoff: ', length_cutoff)
        
        average_coverage = self.assembly.calculate_average_coverage()
        
        with open(tsvfile, 'r', newline='') as f:
            tsvin = csv.reader(f, delimiter='\t')
            for row in tsvin:
                if current_query_id is None:
                    current_query_id = row[0]
                    _hit_list = DiamondHitList(current_query_id)
                
                hit = DiamondHit()
                hit.create_hit(row)
                # filtering by identity and length
                if hit.get_identity() < identity_cutoff:
                    continue # skip this line
                if hit.get_length() < length_cutoff:
                    continue # skip this line

                if hit.get_query_id() != current_query_id:
                    _hit_list.annotate_hits(self.project.ref_data)
                    # compare list of hits from search in background DB with existing hit from search in reference DB
                    current_query_id_tokens = current_query_id.split('|')
                    hit_end = current_query_id_tokens[-1]
                    hit_start = current_query_id_tokens[-2]
                    function_id = current_query_id_tokens[0]
                    contig_tokens = current_query_id_tokens[1].split('_')
                    contig_id = '_'.join(contig_tokens[:-1])
                    gene_id = '|'.join(current_query_id_tokens[:-2])
                    hit_start= int(hit_start)
                    hit_end = int(hit_end)
                    if gene_id in self.assembly.contigs[function_id][contig_id].genes:
                        coverage = self.assembly.contigs[function_id][contig_id].get_coverage()
                        compare_hits_lca(self.assembly.contigs[function_id][contig_id].genes[gene_id], hit_start, hit_end, _hit_list, biscore_range_cutoff, average_coverage, coverage, self.project.taxonomy_data, self.project.ref_data, rank_cutoffs = self.project.config.get_ranks_cutoffs(self.project.options.get_collection()))  # here should be all the magic
                    else:
                        print ('Gene not found: ', gene_id, ' in ', function_id, contig_id)
                        raise TypeError
                    current_query_id = hit.get_query_id()
                    _hit_list = DiamondHitList(current_query_id)
                _hit_list.add_hit(hit)
            _hit_list.annotate_hits(self.project.ref_data)
            current_query_id_tokens = current_query_id.split('|')
            hit_end = current_query_id_tokens[-1]
            hit_start = current_query_id_tokens[-2]
            read_id = '|'.join(current_query_id_tokens[:-2])
            hit_start= int(hit_start)
            hit_end = int(hit_end)
            if gene_id in self.assembly.contigs[function_id][contig_id].genes:
                compare_hits_lca(self.assembly.contigs[function_id][contig_id].genes[gene_id], hit_start, hit_end, _hit_list, biscore_range_cutoff, average_coverage, coverage, self.project.taxonomy_data, self.project.ref_data, rank_cutoffs = self.project.config.get_ranks_cutoffs(self.project.options.get_collection()))  # here should be all the magic
            else:
                print ('Gene not found: ', gene_id)
                raise TypeError


    def parse_uniprot_output(self):
        tsvfile = os.path.join(self.assembly_dir, 'all_contigs_proteins.uniprot.diamondout.txt')
        current_id = ''
        identity_cutoff = self.project.config.get_identity_cutoff(self.project.options.get_collection())
        length_cutoff = self.project.config.get_length_cutoff(self.project.options.get_collection())
        print ('Parse Uniprot output: Identity cutoff: ', identity_cutoff, ', Length cutoff: ', length_cutoff)
        
        with open(tsvfile, 'r', newline='') as f:
            tsvin = csv.reader(f, delimiter='\t')
            for row in tsvin:
                hit = DiamondHit()
                hit.create_hit(row)
                # filtering by identity and length
                if hit.get_identity() < identity_cutoff:
                    continue # skip this line
                if hit.get_length() < length_cutoff:
                    continue # skip this line
                taxonomy_id = self.uniprot.get_uniprot_taxid(hit.get_subject_id())
                gene_id = hit.get_query_id()
                function_id, contig_id, _ = parse_gene_id(gene_id)
                
                if gene_id in self.assembly.contigs[function_id][contig_id].genes:
                    self.assembly.contigs[function_id][contig_id].genes[gene_id].set_uniref_hit(hit)
                    if taxonomy_id:
                        self.assembly.contigs[function_id][contig_id].genes[gene_id].set_taxonomy_id(taxonomy_id)
                else:
                    print ('Gene not found ', gene_id)
                

    def map_genes2uniprot(self):
        run_uniprot_search(self.project)
        self.parse_uniprot_output()        


    def map_genes(self):
        # Filter contigs by coverage
        contig_coverage_cutoff = 3.0
                
        prodigal_infile = os.path.join(self.assembly_dir,'all_contigs.fa')
        with open(prodigal_infile, 'w') as of:
            for function in sorted(self.assembly.contigs.keys()):
                for contig in sorted(self.assembly.contigs[function].keys()):
                    if self.assembly.contigs[function][contig].get_coverage() >= contig_coverage_cutoff:
                        of.write('>' + function + '|' + contig  + '\n')
                        of.write(self.assembly.contigs[function][contig].sequence + '\n')
            of.closed
        
        #~ with open('out.txt', 'w') as of:
            #~ for function in sorted(self.assembly.contigs.keys()):
                #~ for contig in sorted(self.assembly.contigs[function].keys()):
                    #~ of.write(function + '\t' + contig  + '\t' + 
                            #~ str(self.assembly.contigs[function][contig].get_coverage()) + '\t' +
                            #~ str(len(self.assembly.contigs[function][contig].sequence)) + '\t' +
                            #~ ','.join(self.assembly.contigs[function][contig].read_count.keys()) + '\t' +
                            #~ ','.join(str(x) for x in self.assembly.contigs[function][contig].read_count.values()) + '\t' +
                            #~ ','.join(str(x) for x in self.assembly.contigs[function][contig].read_segments.values()) + '\t' +
                             #~ '\n')
            #~ of.closed
        
        # Run Prodigal
        prodigal_outfile = os.path.join(self.assembly_dir,'all_contigs.prodigal.out.faa')
        if not os.path.exists(prodigal_outfile):
            run_prodigal(prodigal_infile, prodigal_outfile, self.project.config.get_prodigal_path())

        with open (prodigal_outfile, 'r') as f:
            current_id = None
            sequence = ''
            for line in f:
                line = line.rstrip('\n\r')
                if line.startswith('>'):
                    if current_id:
                        line_tokens = current_id.split(' # ')
                        function_id,contig_id,gene_id = parse_gene_id(line_tokens[0])
                        #~ print (line_tokens[0],function_id,contig_id,gene_id)
                        gene = Gene(contig_id=contig_id, gene_id=line_tokens[0], sequence=sequence, start=line_tokens[1], end=line_tokens[2], strand=line_tokens[3])
                        self.assembly.contigs[function_id][contig_id].add_gene(gene)
                    line_tokens = line.split(' ')
                    current_id = line[1:] #line_tokens[0][1:]
                    sequence = ''
                else:
                    sequence += line
            line_tokens = current_id.split(' # ')
            function_id,contig_id,gene_id = parse_gene_id(line_tokens[0])
            gene = Gene(contig_id=contig_id, gene_id=line_tokens[0], sequence=sequence, start=line_tokens[1], end=line_tokens[2], strand=line_tokens[3])
            self.assembly.contigs[function_id][contig_id].add_gene(gene)
            f.closed

        # Search in reference database
        if not os.path.exists(os.path.join(self.assembly_dir, 'all_contigs_' + self.project.options.get_ref_output_name())):
            run_ref_search(self.project)
        
        # Process output of reference DB search
        self.parse_reference_output()

        export_gene_assembly(self.assembly, os.path.join(self.assembly_dir, 'all_contigs_assembly.json'))
        
        ##Import sequence data for selected sequence reads
        print ('Reading FASTQ file')
        self.export_hit_fasta()
        
        # Search in background database
        if not os.path.exists(os.path.join(self.assembly_dir, 'all_contigs_'+ self.project.options.get_background_output_name())):
            run_bgr_search(self.project)

        # Process output of reference DB search
        self.parse_background_output()
        
        print('Exporting JSON')
        export_gene_assembly(self.assembly, os.path.join(self.assembly_dir, 'all_contigs_assembly.json'))
        
        # Generate output
        #print('Generating reports')
        #generate_assembly_report(self, taxonomy_data)
        
    def generate_taxonomy_chart(self, taxonomy_data):
        '''
        Collects data about functional genes in assembly and 
        creates one Krona chart for all functions
        '''
        functions_list = set()
        samples_list = sorted(self.project.list_samples())
        genes = autovivify(2) # genes[gene][function][parameter] = parameter_value 
        scores = autovivify(2) # scores[taxonomy ID][function][parameter] = parameter_value 
        

        total_read_count = 0
        for sample in self.project.list_samples():
            total_read_count += self.project.options.get_fastq1_readcount(sample)

        for function in self.assembly.contigs:
            functions_list.add(function)
            for contig in self.assembly.contigs[function]:
                for gene_id in self.assembly.contigs[function][contig].genes:
                    gene = self.assembly.contigs[function][contig].genes[gene_id]
                    if gene.get_status() == 'function':
                        taxonomy_id = gene.taxonomy # Was get_taxonomy_id()
                        for hit in gene.hit_list.get_hits():
                            identity = hit.get_identity()
                            #~ if not taxonomy_id:
                                # If UniRef-based taxonomy ID was not set, guess it from Fama hit
#                                print ('Taxonomy ID is missing for gene ', gene_id)
                                #~ taxonomy_id = self.project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                                #~ gene.set_taxonomy_id(taxonomy_id)
#                            else:
#                                print ('Taxonomy ID for gene ', gene_id, ' is ', taxonomy_id)
                            hit_functions = hit.get_functions()
                            for hit_function in hit_functions:
                                functions_list.add(hit_function)
                                if 'rpkm' in scores[taxonomy_id][hit_function]:
                                    scores[taxonomy_id][hit_function]['rpkm'] += self.assembly.contigs[function][contig].get_rpkm(total_read_count) * len(gene.protein_sequence) * 3 / len(self.assembly.contigs[function][contig].sequence)
                                else:
                                    scores[taxonomy_id][hit_function]['rpkm'] = self.assembly.contigs[function][contig].get_rpkm(total_read_count) * len(gene.protein_sequence) * 3 / len(self.assembly.contigs[function][contig].sequence)
                                if 'count' in scores[taxonomy_id][hit_function]:
                                    scores[taxonomy_id][hit_function]['count'] += self.assembly.contigs[function][contig].get_read_count() * len(gene.protein_sequence) * 3 / len(self.assembly.contigs[function][contig].sequence)
                                else:
                                    scores[taxonomy_id][hit_function]['count'] = self.assembly.contigs[function][contig].get_read_count() * len(gene.protein_sequence) * 3 / len(self.assembly.contigs[function][contig].sequence)
                                if 'hit_count' in scores[taxonomy_id][hit_function]:
                                    scores[taxonomy_id][hit_function]['hit_count'] += 1
                                else:
                                    scores[taxonomy_id][hit_function]['hit_count'] = 1
                                if gene.uniref_hit:
                                    if 'identity' in scores[taxonomy_id][hit_function]:
                                        scores[taxonomy_id][hit_function]['identity'] += gene.uniref_hit.identity
                                    else:
                                        scores[taxonomy_id][hit_function]['identity'] = gene.uniref_hit.identity
                                else:
                                    if 'identity' in scores[taxonomy_id][hit_function]:
                                        scores[taxonomy_id][hit_function]['identity'] += identity
                                    else:
                                        scores[taxonomy_id][hit_function]['identity'] = identity
                                if 'genes' in scores[taxonomy_id][hit_function]:
                                    scores[taxonomy_id][hit_function]['genes'] += gene_id + ' ' #+= gene_id + ', ' + str(len(gene.protein_sequence)) + ' aa(' + '{0:.0f}'.format(len(gene.protein_sequence) * 100 / hit.get_subject_length())+ '%)<br>'
                                else:
                                    scores[taxonomy_id][hit_function]['genes'] = gene_id + ' '# + ', ' + str(len(gene.protein_sequence)) + ' aa(' + '{0:.0f}'.format(len(gene.protein_sequence) * 100 / hit.get_subject_length())+ '%)<br>'
                                    
                                genes[gene_id][hit_function]['Length'] = str(len(gene.protein_sequence)) + 'aa'
                                genes[gene_id][hit_function]['Completeness'] = '{0:.0f}'.format(len(gene.protein_sequence) * 100 / hit.get_subject_length())
                                if gene.uniref_hit:
                                    genes[gene_id][hit_function]['Best hit'] = gene.uniref_hit.subject_id
                                    genes[gene_id][hit_function]['identity'] = '{0:.1f}'.format(gene.uniref_hit.identity)
                                else:
                                    genes[gene_id][hit_function]['identity'] = '{0:.1f}'.format(identity)
                                genes[gene_id][hit_function]['rpkm'] = '{0:.6f}'.format(self.assembly.contigs[function][contig].get_rpkm(total_read_count) * len(gene.protein_sequence) * 3 / len(self.assembly.contigs[function][contig].sequence))
                                genes[gene_id][hit_function]['count'] = '{0:.0f}'.format(self.assembly.contigs[function][contig].get_read_count() * len(gene.protein_sequence) * 3 / len(self.assembly.contigs[function][contig].sequence))
                                genes[gene_id][hit_function]['coverage'] = '{0:.1f}'.format(self.assembly.contigs[function][contig].get_coverage())
        
        taxonomic_profile = TaxonomyProfile()
        taxonomic_profile.build_assembly_taxonomic_profile(taxonomy_data, scores)
#        taxonomic_profile.print_taxonomy_profile()
        outfile = os.path.join(self.assembly_dir, 'assembly_taxonomic_profile.xml')
        generate_assembly_taxonomy_chart(taxonomic_profile, genes, sorted(functions_list), outfile, self.project.config.get_krona_path(), score='rpkm')

    def generate_function_taxonomy_charts(self, taxonomy_data):
        '''
        Generates series of Krona charts visualizing functions in assembly: 
        one function per file, separate stats for each sample
        '''
        functions_list = set()
        samples_list = sorted(self.project.list_samples())
        
        total_read_count = 0
        for sample in self.project.list_samples():
            total_read_count += self.project.options.get_fastq1_readcount(sample)

        # Make list of functions
        for function in self.assembly.contigs:
            for contig in self.assembly.contigs[function]:
                for gene_id in self.assembly.contigs[function][contig].genes:
                    gene_status = self.assembly.contigs[function][contig].genes[gene_id].get_status()
                    if gene_status == 'function':
                        for f in self.assembly.contigs[function][contig].genes[gene_id].functions.keys():
                            functions_list.add(f)
        
        for function in sorted(functions_list):
            genes = autovivify(2) # genes[gene][sample][parameter] = parameter_value 
            scores = autovivify(2) # scores[taxonomy ID][sample][parameter] = parameter_value 
            outfile = os.path.join(self.assembly_dir, function + '_taxonomic_profile.xml')
            for f in self.assembly.contigs:
                for contig in self.assembly.contigs[f]:
                    for gene_id in self.assembly.contigs[f][contig].genes:
                        function_counted = False
                        gene = self.assembly.contigs[f][contig].genes[gene_id]
                        if gene.get_status() == 'function' and function in gene.functions:
                            taxonomy_id = gene.taxonomy # Was get_taxonomy_id()
                            for hit in gene.hit_list.get_hits():
                                identity = hit.get_identity()
                                #~ if not taxonomy_id:
    #                                print ('Taxonomy ID is missing for gene ', gene_id)
                                    #~ taxonomy_id = self.project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                                    #~ gene.set_taxonomy_id(taxonomy_id)
    #                            else:
    #                                print ('Taxonomy ID for gene ', gene_id, ' is ', taxonomy_id)
                                if function in hit.get_functions():
                                    if function_counted:
                                        continue
                                    for sample in samples_list:
                                        if sample in self.assembly.contigs[f][contig].read_count:
                                            if 'rpkm' in scores[taxonomy_id][sample]:
                                                scores[taxonomy_id][sample]['rpkm'] += self.assembly.contigs[f][contig].get_rpkm(self.project.options.get_fastq1_readcount(sample), sample) * len(gene.protein_sequence) * 3 / len(self.assembly.contigs[f][contig].sequence)
                                            else:
                                                scores[taxonomy_id][sample]['rpkm'] = self.assembly.contigs[f][contig].get_rpkm(self.project.options.get_fastq1_readcount(sample), sample) * len(gene.protein_sequence) * 3 / len(self.assembly.contigs[f][contig].sequence)
                                            if 'count' in scores[taxonomy_id][sample]:
                                                scores[taxonomy_id][sample]['count'] += self.assembly.contigs[f][contig].get_read_count(sample) * len(gene.protein_sequence) * 3 / len(self.assembly.contigs[f][contig].sequence)
                                            else:
                                                scores[taxonomy_id][sample]['count'] = self.assembly.contigs[f][contig].get_read_count(sample) * len(gene.protein_sequence) * 3 / len(self.assembly.contigs[f][contig].sequence)
                                            if 'hit_count' in scores[taxonomy_id][sample]:
                                                scores[taxonomy_id][sample]['hit_count'] += 1
                                            else:
                                                scores[taxonomy_id][sample]['hit_count'] = 1
                                            if gene.uniref_hit:
                                                if 'identity' in scores[taxonomy_id][sample]:
                                                    scores[taxonomy_id][sample]['identity'] += gene.uniref_hit.identity
                                                else:
                                                    scores[taxonomy_id][sample]['identity'] = gene.uniref_hit.identity
                                            else:
                                                if 'identity' in scores[taxonomy_id][sample]:
                                                    scores[taxonomy_id][sample]['identity'] += identity
                                                else:
                                                    scores[taxonomy_id][sample]['identity'] = identity
                                            if 'genes' in scores[taxonomy_id][sample]:
                                                scores[taxonomy_id][sample]['genes'] += gene_id + ' ' #+= gene_id + ', ' + str(len(gene.protein_sequence)) + ' aa(' + '{0:.0f}'.format(len(gene.protein_sequence) * 100 / hit.get_subject_length())+ '%)<br>'
                                            else:
                                                scores[taxonomy_id][sample]['genes'] = gene_id + ' '# + ', ' + str(len(gene.protein_sequence)) + ' aa(' + '{0:.0f}'.format(len(gene.protein_sequence) * 100 / hit.get_subject_length())+ '%)<br>'

                                            genes[gene_id][sample]['Length'] = str(len(gene.protein_sequence)) + 'aa'
                                            genes[gene_id][sample]['Completeness'] = '{0:.0f}'.format(len(gene.protein_sequence) * 100 / hit.get_subject_length())
                                            if gene.uniref_hit:
                                                genes[gene_id][sample]['Best hit'] = gene.uniref_hit.subject_id
                                                genes[gene_id][sample]['identity'] = '{0:.1f}'.format(gene.uniref_hit.identity)
                                            else:
                                                genes[gene_id][sample]['identity'] = '{0:.1f}'.format(identity)
                                            genes[gene_id][sample]['rpkm'] = '{0:.7f}'.format(self.assembly.contigs[f][contig].get_rpkm(self.project.options.get_fastq1_readcount(sample), sample) * len(gene.protein_sequence) * 3 / len(self.assembly.contigs[f][contig].sequence))
                                            genes[gene_id][sample]['count'] = '{0:.0f}'.format(self.assembly.contigs[f][contig].get_read_count(sample) * len(gene.protein_sequence) * 3 / len(self.assembly.contigs[f][contig].sequence))
                                            genes[gene_id][sample]['coverage'] = '{0:.1f}'.format(self.assembly.contigs[f][contig].get_coverage(sample))

                                    if 'rpkm' in scores[taxonomy_id]['All samples']:
                                        scores[taxonomy_id]['All samples']['rpkm'] += self.assembly.contigs[f][contig].get_rpkm(total_read_count) * len(gene.protein_sequence) * 3 / len(self.assembly.contigs[f][contig].sequence)
                                    else:
                                        scores[taxonomy_id]['All samples']['rpkm'] = self.assembly.contigs[f][contig].get_rpkm(total_read_count) * len(gene.protein_sequence) * 3 / len(self.assembly.contigs[f][contig].sequence)
                                    if 'count' in scores[taxonomy_id]['All samples']:
                                        scores[taxonomy_id]['All samples']['count'] += self.assembly.contigs[f][contig].get_read_count() * len(gene.protein_sequence) * 3 / len(self.assembly.contigs[f][contig].sequence)
                                    else:
                                        scores[taxonomy_id]['All samples']['count'] = self.assembly.contigs[f][contig].get_read_count() * len(gene.protein_sequence) * 3 / len(self.assembly.contigs[f][contig].sequence)
                                    if 'hit_count' in scores[taxonomy_id]['All samples']:
                                        scores[taxonomy_id]['All samples']['hit_count'] += 1
                                    else:
                                        scores[taxonomy_id]['All samples']['hit_count'] = 1
                                    if gene.uniref_hit:
                                        if 'identity' in scores[taxonomy_id]['All samples']:
                                            scores[taxonomy_id]['All samples']['identity'] += gene.uniref_hit.identity
                                        else:
                                            scores[taxonomy_id]['All samples']['identity'] = gene.uniref_hit.identity
                                    else:
                                        if 'identity' in scores[taxonomy_id]['All samples']:
                                            scores[taxonomy_id]['All samples']['identity'] += identity
                                        else:
                                            scores[taxonomy_id]['All samples']['identity'] = identity
                                    if 'genes' in scores[taxonomy_id]['All samples']:
                                        scores[taxonomy_id]['All samples']['genes'] += gene_id + ' ' #+= gene_id + ', ' + str(len(gene.protein_sequence)) + ' aa(' + '{0:.0f}'.format(len(gene.protein_sequence) * 100 / hit.get_subject_length())+ '%)<br>'
                                    else:
                                        scores[taxonomy_id]['All samples']['genes'] = gene_id + ' '# + ', ' + str(len(gene.protein_sequence)) + ' aa(' + '{0:.0f}'.format(len(gene.protein_sequence) * 100 / hit.get_subject_length())+ '%)<br>'
                                    
                                    genes[gene_id]['All samples']['Length'] = str(len(gene.protein_sequence)) + 'aa'
                                    genes[gene_id]['All samples']['Completeness'] = '{0:.0f}'.format(len(gene.protein_sequence) * 100 / hit.get_subject_length())
                                    if gene.uniref_hit:
                                        genes[gene_id]['All samples']['identity'] = '{0:.1f}'.format(gene.uniref_hit.identity)
                                        genes[gene_id]['All samples']['Best hit'] = gene.uniref_hit.subject_id
                                    else:
                                        genes[gene_id]['All samples']['identity'] = '{0:.1f}'.format(identity)
                                    genes[gene_id]['All samples']['rpkm'] = '{0:.7f}'.format(self.assembly.contigs[f][contig].get_rpkm(total_read_count) * len(gene.protein_sequence) * 3 / len(self.assembly.contigs[f][contig].sequence))
                                    genes[gene_id]['All samples']['count'] = '{0:.0f}'.format(self.assembly.contigs[f][contig].get_read_count() * len(gene.protein_sequence) * 3 / len(self.assembly.contigs[f][contig].sequence))
                                    genes[gene_id]['All samples']['coverage'] = '{0:.1f}'.format(self.assembly.contigs[f][contig].get_coverage())
                                    function_counted = True
            
            taxonomic_profile = TaxonomyProfile()
            taxonomic_profile.build_assembly_taxonomic_profile(taxonomy_data, scores)
    #        taxonomic_profile.print_taxonomy_profile()
            output_sample_ids = sorted(self.project.list_samples())
            output_sample_ids.append('All samples')
            
            generate_assembly_taxonomy_chart(taxonomic_profile, genes, output_sample_ids, outfile, self.project.config.get_krona_path(), score='rpkm')
    
    def generate_output(self):
        taxonomy_data = TaxonomyData(self.project.config)
        taxonomy_data.load_taxdata(self.project.config)

        create_assembly_xlsx(self, taxonomy_data)
        self.generate_taxonomy_chart(taxonomy_data)
        self.generate_function_taxonomy_charts(taxonomy_data)

def run_assembler(functions, assembler, output_dir, is_paired_end = True):
    if assembler.endswith('megahit'):
        run_megahit(functions, output_dir, assembler, is_paired_end)
    elif assembler.endswith('metaspades.py'):
        if is_paired_end:
            run_spades(functions, output_dir, assembler, is_paired_end)
        else:
            raise ProgrammingError('Current version of metaSPAdes does not support single-end libraries.')
        

def run_megahit(functions, output_dir, assembler_command, is_paired_end = True):
    print ('Starting assembly')
    #assembler_command = 'megahit'
    for function in functions:
        print ('Run assembler for function', function)
        if is_paired_end:
            assembler_args = [assembler_command,
                            '-1',
                            os.path.join(output_dir,function + '_pe1.fastq'),
                            '-2',
                            os.path.join(output_dir,function + '_pe2.fastq'),
                            '-o',
                            os.path.join(output_dir,function)
                            ]
        else:
            assembler_args = [assembler_command,
                            '-r',
                            os.path.join(output_dir,function + '_pe1.fastq'),
                            '-o',
                            os.path.join(output_dir,function)
                            ]

        with Popen(assembler_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
            for line in p.stdout:
                print(line, end='')
        if p.returncode != 0:
            print ('Assembler finished with error for function ', function)
            #raise CalledProcessError(p.returncode, p.args)

    print ('Assembly finished')


def run_spades(functions, output_dir, assembler_command, is_paired_end = True):
    print ('Starting metaSPAdes')
    #assembler_command = 'metaspades.py'
    tmp_dir = os.path.join(output_dir, 'tmp')
    for function in functions:
        print ('Run metaSPAdes for function', function)
        assembler_args = [assembler_command,
                        '--meta',
                        '-t',
                        '12',
                        '-m',
                        '50', # TODO: make a parameter
                        '-k',
                        '33,55,71,91,111', # TODO: adjust automatically
                        '-o',
                        os.path.join(output_dir,function),
                        '--tmp-dir',
                        tmp_dir
                        ]
        if is_paired_end:
            assembler_args.extend([
                        '-1',
                        os.path.join(output_dir,function + '_pe1.fastq'),
                        '-2',
                        os.path.join(output_dir,function + '_pe2.fastq')
                        ])
        else:
            assembler_args.extend([
                        '-s', 
                        os.path.join(output_dir,function + '_pe1.fastq')
                        ])
                        

        with Popen(assembler_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
            for line in p.stdout:
                print(line, end='')
        if p.returncode != 0:
            print ('metaSpades finished with error for function ', function)
            #raise CalledProcessError(p.returncode, p.args)
        if os.path.exists(os.path.join(output_dir,function,'contigs.fasta')):
            shutil.copyfile(os.path.join(output_dir,function,'contigs.fasta'),os.path.join(output_dir,function,'final.contigs.fa'))
    print ('metaSPAdes finished')

def run_mapper_indexing(functions, output_dir, mapper_command):
    mapper_command = 'bowtie2-build'
    
    for function in functions:
        if not os.path.exists(os.path.join(output_dir,function, 'final.contigs.filtered.fa')):
            print ('Contigs file for function', function, 'not found')
            continue
        print ('Run indexing for function', function)
        if (os.path.getsize(os.path.join(output_dir,function, 'final.contigs.filtered.fa')) > 0):
            if not os.path.exists(os.path.join(output_dir,function, 'index')):
                os.mkdir(os.path.join(output_dir,function, 'index'))
            mapper_args = [mapper_command,
                            '-f',
                            os.path.join(output_dir,function, 'final.contigs.filtered.fa'),
                            os.path.join(output_dir,function, 'index', 'index')
                            ]

            with Popen(mapper_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
                for line in p.stdout:
                    print(line, end='')
            if p.returncode != 0:
                raise CalledProcessError(p.returncode, p.args)

def run_mapper(functions, output_dir, mapper_command, is_paired_end = True):
    mapper_command = 'bowtie2'
    
    for function in functions:
        if not os.path.exists(os.path.join(output_dir,function, 'final.contigs.filtered.fa')):
            continue
        if (os.path.getsize(os.path.join(output_dir,function, 'final.contigs.filtered.fa')) > 0):
            print ('Run read mapping for function', function)
            if is_paired_end:
                mapper_args = [mapper_command,
                                '-q',
                                '--very-sensitive',
                                '--quiet',
                                '-x',
                                os.path.join(output_dir,function, 'index', 'index'),
                                '-1',
                                os.path.join(output_dir,function + '_pe1.fastq'),
                                '-2',
                                os.path.join(output_dir,function + '_pe2.fastq'),
                                '>' + os.path.join(output_dir,function, 'contigs.sam')
                                ]
            else:
                mapper_args = [mapper_command,
                                '-q',
                                '--very-sensitive',
                                '--quiet',
                                '-x',
                                os.path.join(output_dir,function, 'index', 'index'),
                                '-U',
                                os.path.join(output_dir,function + '_pe1.fastq'),
                                '>' + os.path.join(output_dir,function, 'contigs.sam')
                                ]

            with Popen(mapper_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
                for line in p.stdout:
                    print(line, end='')
            if p.returncode != 0:
                raise CalledProcessError(p.returncode, p.args)

def run_prodigal(infile, outfile, prodigal_path):
    print ('Starting Prodigal')
    prodigal_args = [prodigal_path,
                    '-a',
                    outfile,
                    '-i',
                    infile,
                    '-o',
                    outfile+'prodigal.txt',
                    ]

    with Popen(prodigal_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)

    print ('Prodigal finished')

def run_ref_search(project):
    print ('Starting DIAMOND')
    diamond_args = [project.config.get_diamond_path(),
                    'blastp',
                    '--db',
                    project.config.get_reference_diamond_db(project.options.get_collection()),
                    '--query',
                    os.path.join(project.options.get_assembly_dir(), 'all_contigs.prodigal.out.faa'),
                    '--out',
                    os.path.join(project.options.get_assembly_dir(), 'all_contigs_' + project.options.get_ref_output_name()),
                    '--max-target-seqs',
                    '50',
                    '--evalue',
                    str(project.config.get_evalue_cutoff(project.options.get_collection())),
                    '--threads',
                    project.config.get_threads(),
                    '--outfmt','6','qseqid','sseqid','pident','length','mismatch','slen','qstart','qend','sstart','send','evalue','bitscore'
                    ]

    with Popen(diamond_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)

    print ('DIAMOND finished')

def run_bgr_search(project):
    print ('Starting DIAMOND')
    diamond_args = [project.config.get_diamond_path(),
                    'blastp',
                    '--db',
                    project.config.get_background_diamond_db(project.options.get_collection()),
                    '--query',
                    os.path.join(project.options.get_assembly_dir(), 'all_contigs_'+ project.options.get_ref_hits_fastq_name()),
                    '--out',
                    os.path.join(project.options.get_assembly_dir(), 'all_contigs_'+ project.options.get_background_output_name()),
                    '--max-target-seqs',
                    '50',
                    '--evalue',
                    str(project.config.get_background_db_size(project.options.get_collection()) 
                        * project.config.get_evalue_cutoff(project.options.get_collection())
                        / project.config.get_reference_db_size(project.options.get_collection())),
                    '--threads',
                    project.config.get_threads(),
                    '--outfmt','6','qseqid','sseqid','pident','length','mismatch','slen','qstart','qend','sstart','send','evalue','bitscore'
                    ]

    with Popen(diamond_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)

    print ('DIAMOND finished')

def run_uniprot_search(project):
    print ('Starting DIAMOND')
    outfile = os.path.join(project.options.get_assembly_dir(), 'all_contigs_proteins.uniprot.diamondout.txt')
    
    diamond_args = [project.config.get_diamond_path(),
                    'blastp',
                    '--db',
                    project.config.get_uniprot_diamond_db(),
                    '--query',
                    os.path.join(project.options.get_assembly_dir(), 'all_contigs.prodigal.out.faa'),
                    '--out',
                    outfile,
                    '--max-target-seqs',
                    '1',
                    '--evalue',
                    str(project.config.get_background_db_size(project.options.get_collection()) 
                        * project.config.get_evalue_cutoff(project.options.get_collection())
                        / project.config.get_reference_db_size(project.options.get_collection())),
                    '--threads',
                    project.config.get_threads(),
                    '--outfmt','6','qseqid','sseqid','pident','length','mismatch','slen','qstart','qend','sstart','send','evalue','bitscore'
                    ]

    if not os.path.exists(outfile):
        with Popen(diamond_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
            for line in p.stdout:
                print(line, end='')
        if p.returncode != 0:
            raise CalledProcessError(p.returncode, p.args)

    print ('DIAMOND finished')

def parse_gene_id(gene_id):
    (function_id, gene) = gene_id.split('|')
    gene_id_tokens = gene.split('_')
    gene_id = gene_id_tokens[-1]
    contig_id = '_'.join(gene_id_tokens[:-1])
    return function_id, contig_id, gene_id

def get_abundance(function_fraction, average_coverage, coverage):
    if function_fraction > 1.0:
        print('FUNCTION FRACTION TOO BIG!', function_fraction)
    #else:
    #    print('FUNCTION FRACTION ', function_fraction)
    ret_val = coverage * function_fraction/average_coverage
    return ret_val


def compare_hits_lca(gene, hit_start, hit_end, new_hit_list, bitscore_range_cutoff, average_coverage, coverage, taxonomy_data, ref_data, rank_cutoffs = {}):
    # This function compares hits assigned to an annotated read with functions
    # from a Diamond hit list. It looks through the hit list, finds 
    # hits with bitscore above cutoff and takes their functions.
    #
    # If there is one hit with the highest bit-score, read gets status 'function,besthit'
    # If there are several hits with the highest bit-score, read gets status 'function'
    # Otherwise, read gets status 'nofunction'
    #
    # hit_start and hit_end parameters are used for identification of hit for
    # comparison, since multiple hits can be associated with a read 
    #
    # This function does not return anything. It sets status of read and 
    # assigns RPKM score to each function of the read
    #
    # Find best hit
    
    for hit in gene.hit_list.get_hits():
        if hit.get_query_start() == hit_start and hit.get_query_end() == hit_end:
            best_bitscore = 0.0
            best_hit = None
            for new_hit in new_hit_list.get_hits():
                bitscore = new_hit.get_bitscore()
                if bitscore > best_bitscore:
                    best_hit = new_hit
                    best_bitscore = bitscore
            # Set status of read
            if best_hit != None:
                if '' in best_hit.get_functions():
                    gene.set_status('nofunction')
                    return
                else:
                    gene.set_status('function')
            else:
                gene.set_status('nofunction')
                return
            
            # Filter list of hits by bitscore
            bitscore_lower_cutoff = best_bitscore * (1.0 - bitscore_range_cutoff)
            new_hits = [new_hit for new_hit in new_hit_list.get_hits() if new_hit.get_bitscore() > bitscore_lower_cutoff]
                
            if hit.get_subject_id() not in [new_hit.get_subject_id() for new_hit in new_hits] and hit.get_bitscore() >= best_bitscore:
                    new_hits.append(hit)

            # Collect taxonomy IDs of all hits for LCA inference
            taxonomy_ids = set()
            # If rank-specific AAI cutoffs are not set
            if len(rank_cutoffs) == 0:
                taxonomy_ids = set([ref_data.lookup_protein_tax(h.get_subject_id()) for h in new_hits])

            # If rank-specific AAI cutoffs were calculated for the reference dataset:
            else:
                for h in new_hits:
                    subject_taxon_id = ref_data.lookup_protein_tax(h.get_subject_id())
                    hit_identity = h.get_identity()
                    subject_rank = taxonomy_data.get_taxonomy_rank(subject_taxon_id)
                    while subject_taxon_id != taxonomy_data.ROOT:
                        if subject_rank not in rank_cutoffs:
                            subject_taxon_id, subject_rank = taxonomy_data.get_upper_level_taxon(subject_taxon_id)
                        elif hit_identity < rank_cutoffs[subject_rank]:
                            subject_taxon_id, subject_rank = taxonomy_data.get_upper_level_taxon(subject_taxon_id)
                        else:
                            taxonomy_ids.add(subject_taxon_id)
                            break

            # Make non-redundant list of functions from hits after filtering
            new_functions = {}
            new_functions_counter = Counter()
            new_functions_dict = defaultdict(dict)
            # Find best hit for each function: only one hit with highest bitscore to be reported for each function
            for h in new_hits:
                for f in h.get_functions():
                    new_functions_counter[f] += 1
                    if f in new_functions_dict:
                        if h.get_bitscore() > new_functions_dict[f]['bit_score']:
                            new_functions_dict[f]['bit_score'] = h.get_bitscore()
                            new_functions_dict[f]['hit'] = h
                    else:
                        new_functions_dict[f]['bit_score'] = h.get_bitscore()
                        new_functions_dict[f]['hit'] = h

            # If the most common function in new hits is unknown, set status "nofunction" and return
            if new_functions_counter.most_common(1)[0][0] == '':
                gene.set_status('nofunction')
                return

            # Calculate RPK scores for functions
            for function in new_functions_dict:
                if function == '':
                    continue
                new_functions[function] = get_abundance(1.0, average_coverage, coverage)

            gene.set_functions(new_functions)

            # Set new list of hits
            _hit_list = DiamondHitList(gene.gene_id)
            for f in new_functions_dict:
                if f == '':
                    continue
                good_hit = new_functions_dict[f]['hit']
                good_hit.query_id = gene.gene_id
                good_hit.annotate_hit(ref_data)
                _hit_list.add_hit(good_hit)
            
            gene.set_hit_list(_hit_list)
            # Set read taxonomy ID 
            gene.taxonomy = taxonomy_data.get_lca(taxonomy_ids)


def cleanup_read_id(read_id):
    if read_id.endswith('.1'):
        return read_id[:-2]
    elif read_id.endswith('.2'):
        return read_id[:-2]
    else:
        return read_id
