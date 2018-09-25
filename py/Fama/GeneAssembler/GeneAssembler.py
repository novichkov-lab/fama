import os, csv, operator, shutil
from subprocess import Popen, PIPE, CalledProcessError
from collections import Counter, defaultdict

from Fama.DiamondParser.hit_utils import autovivify,cleanup_protein_id
from Fama.GeneAssembler.Contig import Contig
from Fama.GeneAssembler.Gene import Gene
from Fama.GeneAssembler.GeneAssembly import GeneAssembly
from Fama.DiamondParser.DiamondHitList import DiamondHitList
from Fama.DiamondParser.DiamondHit import DiamondHit
from Fama.OutputUtil.JSONUtil import export_gene_assembly
from Fama.TaxonomyProfile import TaxonomyProfile
from Fama.OutputUtil.KronaXMLWriter import generate_assembly_taxonomy_chart
from Fama.OutputUtil.Report import create_assembly_xlsx,generate_assembly_report
from Fama.ReferenceLibrary.TaxonomyData import TaxonomyData
from Fama.ReferenceLibrary.UniprotData import UniprotData

class GeneAssembler:
    def __init__(self, project):
        self.project = project
        self.assembly = GeneAssembly()
        self.assembly_dir = os.path.join(self.project.options.get_work_dir(),'assembly')
        if not os.path.isdir(self.assembly_dir):
            os.mkdir(self.assembly_dir)
        if not os.path.isdir(os.path.join(self.assembly_dir,'out')):
            os.mkdir(os.path.join(self.assembly_dir,'out'))
        self.uniprot = UniprotData(self.project.config)


    def assemble_contigs(self):
        # Export reads in FASTQ formats
        for sample in sorted(self.project.list_samples()):
            for end in ('pe1','pe2'):
                print ('Loading mapped reads: ', sample, end)
                self.project.load_annotated_reads(sample, end) # Lazy load
                for read_id in self.project.samples[sample][end]:
                    read = self.project.samples[sample][end][read_id]
#                    print(read_id)
                    if read.get_status() == 'function,besthit' or read.get_status() == 'function':
                        #read_functions = project.samples[sample][end][read].get_functions()
                        for function in read.get_functions():
                            if read_id in self.assembly.reads[function]:
                                continue
                            self.assembly.reads[function][read_id] = sample
                            outfile1 = os.path.join(self.assembly_dir,function + '_pe1.fastq')
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
                            with open(outfile2, 'a') as of2:
                                of2.write(read.pe_id + '\n')
                                of2.write(read.pe_sequence + '\n')
                                of2.write(read.pe_line3 + '\n')
                                of2.write(read.pe_quality + '\n')
                                of2.closed
                # Delete reads from memory
                self.project.samples[sample][end] = None
        
        # Run Assembler ('megahit' for Megahit or 'metaSPAdes' for metaSPAdes)
        #run_assembler(sorted(self.assembly.reads.keys()), 'megahit', self.assembly_dir)
        run_assembler(sorted(self.assembly.reads.keys()), 'metaSPAdes', self.assembly_dir)
        self.filter_contigs_by_length()

        # Run Bowtie
        
        run_mapper_indexing(sorted(self.assembly.reads.keys()), self.assembly_dir)
        run_mapper(sorted(self.assembly.reads.keys()), self.assembly_dir)

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
                            line_tokens = line.split(' ')
                            current_id = line_tokens[0][1:]
                            sequence = ''
                        else:
                            sequence += line
                    f.closed
            else:
                print('File ' + contig_file + ' does not exist.')

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
                    f.closed
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
                        print(current_id,function_id,contig_id,gene_id)
                        print('Genes:',self.assembly.contigs[function_id][contig_id].genes.keys())
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
        
        mean_coverage = self.assembly.calculate_mean_coverage()
        
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
                    #print (read_id, hit_start, hit_end, biscore_range_cutoff)
                    #print (_hit_list.print_hits())
                    if gene_id in self.assembly.contigs[function_id][contig_id].genes:
                        coverage = self.assembly.contigs[function_id][contig_id].get_coverage()
                        compare_hits(self.assembly.contigs[function_id][contig_id].genes[gene_id], hit_start, hit_end, _hit_list, biscore_range_cutoff, mean_coverage, coverage) # here should be all the magic
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
                compare_hits(self.assembly.contigs[function_id][contig_id].genes[gene_id], hit_start, hit_end, _hit_list, biscore_range_cutoff, mean_coverage, coverage) # here should be all the magic
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
                
                if gene_id in self.assembly.contigs[function_id][contig_id].genes and taxonomy_id:
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
        
        with open('out.txt', 'w') as of:
            for function in sorted(self.assembly.contigs.keys()):
                for contig in sorted(self.assembly.contigs[function].keys()):
                    of.write(function + '\t' + contig  + '\t' + 
                            str(self.assembly.contigs[function][contig].get_coverage()) + '\t' +
                            str(len(self.assembly.contigs[function][contig].sequence)) + '\t' +
                            ','.join(self.assembly.contigs[function][contig].read_count.keys()) + '\t' +
                            ','.join(str(x) for x in self.assembly.contigs[function][contig].read_count.values()) + '\t' +
                            ','.join(str(x) for x in self.assembly.contigs[function][contig].read_segments.values()) + '\t' +
                             '\n')
            of.closed
        
        # Run Prodigal
        prodigal_outfile = os.path.join(self.assembly_dir,'all_contigs.prodigal.out.faa')
        run_prodigal(prodigal_infile, prodigal_outfile)

        with open (prodigal_outfile, 'r') as f:
            current_id = None
            sequence = ''
            for line in f:
                line = line.rstrip('\n\r')
                if line.startswith('>'):
                    if current_id:
                        function_id,contig_id,gene_id = parse_gene_id(current_id)
                        print (current_id,function_id,contig_id,gene_id)
                        gene = Gene(contig_id = contig_id, gene_id = current_id, sequence = sequence)
                        self.assembly.contigs[function_id][contig_id].add_gene(gene)
                    line_tokens = line.split(' ')
                    current_id = line_tokens[0][1:]
                    sequence = ''
                else:
                    sequence += line
            function_id,contig_id,gene_id = parse_gene_id(current_id)
            gene = Gene(contig_id = contig_id, gene_id = current_id, sequence = sequence)
            self.assembly.contigs[function_id][contig_id].add_gene(gene)
            f.closed

        # Search in reference database
        run_ref_search(self.project)
        
        # Process output of reference DB search
        self.parse_reference_output()

        export_gene_assembly(self.assembly, os.path.join(self.assembly_dir, 'all_contigs_assembly.json'))
        
        ##Import sequence data for selected sequence reads
        print ('Reading FASTQ file')
        self.export_hit_fasta()
        
        # Search in background database
        run_bgr_search(self.project)

        # Process output of reference DB search
        self.parse_background_output()
        
        print('Exporting JSON')
        export_gene_assembly(self.assembly, os.path.join(self.assembly_dir, 'all_contigs_assembly.json'))

        
        # Generate output
        #print('Generating reports')
        #generate_assembly_report(self, taxonomy_data)
        
    def generate_taxonomy_chart(self, taxonomy_data):


        functions_list = set()
        samples_list = sorted(self.project.list_samples())
        genes = autovivify(2) # genes[gene][function][parameter] = parameter_value 
        scores = autovivify(2) # scores[taxonomy ID][function][parameter] = parameter_value 
        
#        functions_of_interest = ['GH5_4','GH10','GH55', 'GH1', 'GH6', 'GH9']

        total_read_count = 0
        for sample in self.project.list_samples():
            total_read_count += self.project.options.get_fastq1_readcount(sample)

        for function in self.assembly.contigs:
#            if function not in functions_of_interest:
#                continue
            functions_list.add(function)
            for contig in self.assembly.contigs[function]:
                for gene_id in self.assembly.contigs[function][contig].genes:
                    gene = self.assembly.contigs[function][contig].genes[gene_id]
                    if gene.get_status() == 'function,besthit' or gene.get_status() == 'function':
                        for hit in gene.hit_list.get_hits():
                            identity = hit.get_identity()
                            taxonomy_id = gene.get_taxonomy_id()
                            if not taxonomy_id:
#                                print ('Taxonomy ID is missing for gene ', gene_id)
                                taxonomy_id = self.project.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                                gene.set_taxonomy_id(taxonomy_id)
#                            else:
#                                print ('Taxonomy ID for gene ', gene_id, ' is ', taxonomy_id)
                            hit_functions = hit.get_functions()
                            for hit_function in hit_functions:
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
                                genes[gene_id][hit_function]['identity'] = '{0:.1f}'.format(identity)
                                genes[gene_id][hit_function]['rpkm'] = '{0:.6f}'.format(self.assembly.contigs[function][contig].get_rpkm(total_read_count) * len(gene.protein_sequence) * 3 / len(self.assembly.contigs[function][contig].sequence))
                                genes[gene_id][hit_function]['count'] = '{0:.0f}'.format(self.assembly.contigs[function][contig].get_read_count() * len(gene.protein_sequence) * 3 / len(self.assembly.contigs[function][contig].sequence))
                                genes[gene_id][hit_function]['coverage'] = '{0:.1f}'.format(self.assembly.contigs[function][contig].get_coverage())
        
        taxonomic_profile = TaxonomyProfile()
        taxonomic_profile.build_assembly_taxonomic_profile(taxonomy_data, scores)
#        taxonomic_profile.print_taxonomy_profile()
        outfile = os.path.join(self.assembly_dir, 'assembly_taxonomic_profile.xml')
        generate_assembly_taxonomy_chart(taxonomic_profile, genes, sorted(functions_list), outfile, score='rpkm')
    
    def generate_output(self):
        taxonomy_data = TaxonomyData(self.project.config)
        taxonomy_data.load_taxdata(self.project.config)

        create_assembly_xlsx(self, taxonomy_data)
        #self.generate_taxonomy_chart(taxonomy_data)

def run_assembler(functions, assembler, output_dir):
    if assembler == 'megahit':
        run_megahit(functions, output_dir)
    elif assembler == 'metaSPAdes':
        run_spades(functions, output_dir)

def run_megahit(functions, output_dir):
    print ('Starting assembly')
    assembler_command = '/home/aekazakov/Soft/Megahit/megahit'
    for function in functions:
        print ('Run assembler for function', function)
        assembler_args = [assembler_command,
                        '-1',
                        os.path.join(output_dir,function + '_pe1.fastq'),
                        '-2',
                        os.path.join(output_dir,function + '_pe2.fastq'),
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


def run_spades(functions, output_dir):
    print ('Starting metaSPAdes')
    assembler_command = '/home/aekazakov/Soft/SPAdes/SPAdes-3.10.1-Linux/bin/metaspades.py'
    tmp_dir = '/home/aekazakov/Soft/SPAdes/SPAdes-3.10.1-Linux/bin/tmp'
    for function in functions:
        print ('Run metaSPAdes for function', function)
        assembler_args = [assembler_command,
                        '--meta',
                        '-t',
                        '12',
                        '-m',
                        '50',
                        '-o',
                        os.path.join(output_dir,function),
                        '--tmp-dir',
                        tmp_dir,
                        '-1',
                        os.path.join(output_dir,function + '_pe1.fastq'),
                        '-2',
                        os.path.join(output_dir,function + '_pe2.fastq')
                        ]

        with Popen(assembler_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
            for line in p.stdout:
                print(line, end='')
        if p.returncode != 0:
            print ('metaSpades finished with error for function ', function)
            #raise CalledProcessError(p.returncode, p.args)
        if os.path.exists(os.path.join(output_dir,function,'contigs.fasta')):
            shutil.copyfile(os.path.join(output_dir,function,'contigs.fasta'),os.path.join(output_dir,function,'final.contigs.fa'))
    print ('metaSPAdes finished')

def run_mapper_indexing(functions, output_dir):
    mapper_command = 'bowtie2-build'
    
    for function in functions:
        if not os.path.exists(os.path.join(output_dir,function, 'final.contigs.filtered.fa')):
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

def run_mapper(functions, output_dir):
    mapper_command = 'bowtie2'
    
    for function in functions:
        if not os.path.exists(os.path.join(output_dir,function, 'final.contigs.filtered.fa')):
            continue
        if (os.path.getsize(os.path.join(output_dir,function, 'final.contigs.filtered.fa')) > 0):
            print ('Run read mapping for function', function)
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

            with Popen(mapper_args, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
                for line in p.stdout:
                    print(line, end='')
            if p.returncode != 0:
                raise CalledProcessError(p.returncode, p.args)

def run_prodigal(infile, outfile):
    print ('Starting Prodigal')
    prodigal_args = ['prodigal',
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
    diamond_args = ['/usr/bin/diamond',
                    'blastp',
                    '--db',
                    project.config.get_reference_diamond_db(project.options.get_collection()),
                    '--query',
                    os.path.join(project.options.get_work_dir(), 'assembly', 'all_contigs.prodigal.out.faa'),
                    '--out',
                    os.path.join(project.options.get_work_dir(), 'assembly', 'all_contigs_' + project.options.get_ref_output_name()),
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
    diamond_args = ['/usr/bin/diamond',
                    'blastp',
                    '--db',
                    project.config.get_background_diamond_db(project.options.get_collection()),
                    '--query',
                    os.path.join(project.options.get_work_dir(), 'assembly', 'all_contigs_'+ project.options.get_ref_hits_fastq_name()),
                    '--out',
                    os.path.join(project.options.get_work_dir(), 'assembly', 'all_contigs_'+ project.options.get_background_output_name()),
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
    uniprot_db_path = '/mnt/data2/Databases/UniRef/uniref100.20171108.dmnd'
    outfile = os.path.join(project.options.get_work_dir(), 'assembly', 'all_contigs_proteins.uniprot.diamondout.txt')
    
    diamond_args = ['/usr/bin/diamond',
                    'blastp',
                    '--db',
                    uniprot_db_path,
                    '--query',
                    os.path.join(project.options.get_work_dir(), 'assembly', 'all_contigs.prodigal.out.faa'),
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



def get_abundance(function_fraction, mean_coverage, coverage):
    if function_fraction > 1.0:
        print('FUNCTION FRACTION TOO BIG!', function_fraction)
    else:
        print('FUNCTION FRACTION ', function_fraction)
    ret_val = coverage * function_fraction/mean_coverage
    return ret_val


def compare_hits(gene, hit_start, hit_end, new_hit_list, bitscore_range_cutoff, mean_coverage, coverage):
    # This functions compares hits assigned to an annotated read with functions
    # from a Diamond hit list
    #
    # hit_start and hit_end parameters are used for identification of hit for
    # comparison, since multiple hits can be associated with a read 
    #
    # This function does not return anything. It sets status of read and 
    # function counter of the read through read methods
    
    for hit in gene.hit_list.get_hits():
        #print(str(hit))
        #print (str(hit.get_query_start()), str(hit_start), str(hit.get_query_end()), str(hit_end))
        #print (type(hit.get_query_start()), type(hit_start), type(hit.get_query_end()), type(hit_end))
        if hit.get_query_start() == hit_start and hit.get_query_end() == hit_end:
            #print ('Start comparison')
            bitscore = hit.get_bitscore()
            bitscore_lower_cutoff = bitscore * (1 - bitscore_range_cutoff)
            bitscore_upper_cutoff = bitscore * (1 + bitscore_range_cutoff)
#            print('Cutoffs:',bitscore_lower_cutoff,bitscore_upper_cutoff)
            # first, make a list of hits with acceptable bitscore values (i.e. within given range):
            new_hits = [new_hit for new_hit in new_hit_list.get_hits() if new_hit.get_bitscore() > bitscore_lower_cutoff]
#            print ('Hits found: ', len(new_hits) or 0)
            if not new_hits:
                print ('case 0')
                print (hit)
                gene.set_status('nofunction')
                # nothing to do here
                
            elif len(new_hits) == 1:
#                print(new_hits[0])
                # if only one hit left, function assignment is very easy
#                print ('case 1: single hit')
#                print (cleanup_protein_id(new_hits[0].get_subject_id()),', ', cleanup_protein_id(hit.get_subject_id()))
                new_functions = {}
                functions = compare_functions(hit, new_hits)
                if cleanup_protein_id(new_hits[0].get_subject_id()) == cleanup_protein_id(hit.get_subject_id()):
                    # this is the same top hit as before
#                    print ('case 1.1')
                    gene.set_status('function,besthit')                    
                    print(functions)
                    total_count = sum(functions.values())#len(new_hits[0].get_functions())
                    for function in new_hits[0].get_functions():
                        new_functions[function] = get_abundance(functions[function]/total_count, mean_coverage, coverage)
                elif '' in functions:
                    # function unknown
#                    print ('case 1.3')
                    gene.set_status('nofunction')
                    total_count = sum(functions.values())
                    return
                else:
#                    print ('case 1.2')
                    gene.set_status('function')
                    total_count = sum(functions.values())
                    for function in functions:
                        new_functions[function] = get_abundance(functions[function]/total_count, mean_coverage, coverage)
                gene.set_functions(new_functions)
                
            else:
 #               print ('case 2: multiple hits')
                # But what if top hit in background DB search is different from the top hit in reference DB search?
                #
                # Basically, several cases are possible:
                # 1. True best hits are not in reference DB, i.e. read function is different.
                #       We must check if a function of top refDB hit is present in list of functions of new_hits list.
                #       If most of proteins are not in the reference database, this read must have no function assigned.
                # 2. There are two close proteins in reference DB, and they switched places in background DB search.
                #       In this case, function of top hits would remain the same. Compare two lists of functions.
                # 3. Hit sequence is nearly equally distant from proteins of interesting function and proteins with other functions.
                #       Compare lists of functions. If most of proteins are not in the reference database, this read must have no function assigned.
                # 4. Top hit in background DB was misannotated. In this case, next hits close to top will have good function.
                #       Compare lists of functions. If most of proteins ARE in the reference database, this read must have right function assigned.
                # 

                functions = compare_functions(hit, new_hits)
                if '' in functions and functions[''] == 0:
#                        print ('case 2.0')
                        gene.set_status('nofunction')
                        return

                if new_hits[0].get_bitscore() > bitscore_upper_cutoff:
                    # we need to refine new_hits list
                    new_bitscore_lower_cutoff = new_hits[0].get_bitscore() * (1 - bitscore_range_cutoff)
                    new_hits = [new_hit for new_hit in new_hits if new_hit.get_bitscore() > new_bitscore_lower_cutoff]
                    new_functions = {}
                    functions = compare_functions(hit, new_hits)
                    if '' in functions and functions[''] == 0: 
#                        print ('case 2.0') # very unlikely
                        gene.set_status('nofunction')
                        return
                    if len(functions) == 1 and '' in functions:
#                        print ('case 2.5')
                        gene.set_status('nofunction')
                        return
                    else:
#                        print ('case 2.3')
                        gene.set_status('function')
                        total_count = sum(functions.values())
                        for function in functions:
                            new_functions[function] = get_abundance(functions[function]/total_count, mean_coverage, coverage)
                    gene.set_functions(new_functions)
                else:
#                    for hit1 in new_hits:
#                        print(hit1)
#                    print(functions)
                    new_functions = {}
                    if len(functions) == 1 and '' in functions:
#                        print ('case 2.4')
                        gene.set_status('nofunction')
                        return
                    elif cleanup_protein_id(new_hits[0].get_subject_id()) == cleanup_protein_id(hit.get_subject_id()):
#                        print ('case 2.1')
                        gene.set_status('function,besthit')
                        total_count = sum(functions.values())#len(new_hits[0].get_functions())
                        for function in functions:
                            if function in new_hits[0].get_functions():
                                new_functions[function] = get_abundance(functions[function]/total_count, mean_coverage, coverage)
                        if not new_functions:
                            for function in functions:
                                new_functions[function] = get_abundance(functions[function]/total_count, mean_coverage, coverage)
#                        print(hit)
#                        print(new_functions)
                    else:
                        # the most interesting: best hit is close to top hit in reference DB search
                        print ('case 2.2')
                        gene.set_status('function')
                        total_count = sum(functions.values())
                        for function in functions:
                            new_functions[function] = get_abundance(functions[function]/total_count, mean_coverage, coverage)
                    gene.set_functions(new_functions)
        else:
#            print('Skipping hit',hit.get_query_id())
            pass
            

def compare_functions(hit, new_hits):
    # This function compares two lists of functions: one list assigned to a single hit
    # and other list of functions assigned to a list of hits. 
    # It returns dictionary of functions and counts for each function
    ret_val = {}
    old_functions = hit.get_functions()
    new_functions_counter = Counter()
    for hit in new_hits:
        for function in hit.get_functions():
            new_functions_counter[function] += 1
    # first, choose minimal count of hit for a function. List of functions may be very long,
    # but we consider only top of the list. The size of the top depends on number of functions
    # assigned to the old hit (typically, one). But if we have more than one top function with equal 
    # count of genes, the original function will not be the top one. So, we should consider
    # all functions with hit counts equal to the count of the top hit.
    minimal_count = 0
    if len(new_functions_counter) > len(old_functions):
        minimal_count = new_functions_counter.most_common(len(old_functions))[-1][1]
    else:
        minimal_count = new_functions_counter.most_common()[-1][1]
    # second, let's truncate new_functions_counter, taking only elements with count equal or above minimal_count
    new_functions = {i[0]:i[1] for i in new_functions_counter.most_common() if i[1] >= minimal_count}
    # if new_functions dict is empty after all, add empty value into ret_val and return
    if not new_functions:
        ret_val[''] = 0
        return ret_val
    else:
        # next, compare keys of new_functions dict with the list of functions of the old hit
        # if most of genes have no function, return only one element
        for old_function in old_functions:
            if old_function in new_functions:
                ret_val[old_function] = new_functions[old_function]

        # if new_functions dict is empty after that (i.e. old functions are not at the top of 
        # new functions list), return count of the top function 
        if not ret_val:
            top_function = max(new_functions.items(), key=operator.itemgetter(1))[0]
            ret_val[top_function] = new_functions[top_function]
        return ret_val

