import os, csv, re
import gzip
from Fama.ProjectUtil.ProgramConfig import ProgramConfig
from Fama.ProjectUtil.ProjectOptions import ProjectOptions
from Fama.ReferenceLibrary.ReferenceData import ReferenceData
from Fama.DiamondParser.DiamondHit import DiamondHit
from Fama.DiamondParser.DiamondHitList import DiamondHitList
from Fama.ReadUtil.AnnotatedRead import AnnotatedRead
from Fama.DiamondParser.hit_utils import cleanup_protein_id,get_rpkm_score,compare_hits,get_paired_end,get_paired_read_id,compare_hits_naive

class DiamondParser:

    def __init__(self, sample, end, config_file=None, project_file=None, config=None, project=None, ref_data=None):
        self.reads = {}
        self.sample = sample
        self.end = end
        self.config = config
        if not self.config:
            self.config = ProgramConfig(config_file)
        self.project = project
        if not self.project:
            self.project = ProjectOptions(project_file)
        collection = self.project.get_collection(self.sample)
        if collection not in self.config.list_collections():
            raise Exception ('Collection ' + collection + ' not found. Available collections are: ' + (',').join(colelctions))
        self.collection = collection
        self.ref_data = ref_data
        if not self.ref_data:
            self.ref_data = ReferenceData(self.config)
            self.ref_data.load_reference_data(self.collection)

    def parse_reference_output(self):
        
        tsvfile = os.path.join(self.project.get_project_dir(self.sample), self.sample + '_' + self.end + '_'+ self.project.get_ref_output_name())
        #add paired-end option
        
        current_sequence_read_id = ''
        _hit_list = DiamondHitList(current_sequence_read_id)
        identity_cutoff = self.config.get_identity_cutoff(self.collection)
        length_cutoff = self.config.get_length_cutoff(self.collection)
        print ('Identity cutoff: ', identity_cutoff, ', Length cutoff: ', length_cutoff)
        
        with open(tsvfile, 'r', newline='') as f:
            tsvin = csv.reader(f, delimiter='\t')
            for row in tsvin:
                hit = DiamondHit()
                (row[0], _ ) = self.parse_fastq_seqid(row[0])
                hit.create_hit(row)
                # filtering by identity and length
                if hit.get_identity() < identity_cutoff:
                    continue # skip this line
                if hit.get_length() < length_cutoff:
                    continue # skip this line

                if hit.get_query_id() != current_sequence_read_id:
                    # filter list for overlapping hits
                    _hit_list.filter_list(self.config.get_overlap_cutoff(self.collection))
                    if _hit_list.get_hits_number() != 0:
                        # annotate_hits
                        _hit_list.annotate_hits(self.ref_data)
                        read = AnnotatedRead(current_sequence_read_id)
                        read.set_hit_list(_hit_list)
                        self.reads[current_sequence_read_id] = read

                    current_sequence_read_id = hit.get_query_id()
                    _hit_list = DiamondHitList(current_sequence_read_id)
                _hit_list.add_hit(hit)
            if _hit_list.get_hits_number() != 0:
                _hit_list.filter_list(self.config.get_overlap_cutoff(self.collection))
                # annotate_hits
                _hit_list.annotate_hits(self.ref_data)
                read = AnnotatedRead(current_sequence_read_id)
                read.set_hit_list(_hit_list)
                self.reads[current_sequence_read_id] = read

    def parse_background_output(self):
        
        if len(self.reads) == 0:
            self.reads = self.import_hit_list()
        
        tsvfile = os.path.join(self.project.get_project_dir(self.sample), self.sample + '_' + self.end + '_'+ self.project.get_background_output_name())
        #add paired-end option
        
        current_query_id = None
        _hit_list = None
        identity_cutoff = self.config.get_identity_cutoff(self.collection)
        length_cutoff = self.config.get_length_cutoff(self.collection)
        biscore_range_cutoff = self.config.get_biscore_range_cutoff(self.collection)
        print ('Identity cutoff: ', identity_cutoff, ', Length cutoff: ', length_cutoff)
        
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
                    _hit_list.annotate_hits(self.ref_data)
                    # compare list of hits from search in background DB with existing hit from search in reference DB
                    (read_id, hit_start, hit_end) = current_query_id.split('|')
                    hit_start= int(hit_start)
                    hit_end = int(hit_end)
                    #print (read_id, hit_start, hit_end, biscore_range_cutoff)
                    #print (_hit_list.print_hits())
                    if read_id in self.reads.keys():
                        #compare_hits(self.reads[read_id], hit_start, hit_end, _hit_list, biscore_range_cutoff, length_cutoff, self.project.get_fastq1_readcount(self.sample)) # here should be all the magic
                        compare_hits_naive(self.reads[read_id], hit_start, hit_end, _hit_list, biscore_range_cutoff, length_cutoff, self.project.get_fastq1_readcount(self.sample)) # here should be all the magic
                    else:
                        print ('Read not found: ', read_id)
#                        raise TypeError
                    current_query_id = hit.get_query_id()
                    _hit_list = DiamondHitList(current_query_id)
                _hit_list.add_hit(hit)
            _hit_list.annotate_hits(self.ref_data)
            (read_id, hit_start, hit_end) = current_query_id.split('|')
            hit_start= int(hit_start)
            hit_end = int(hit_end)
            if read_id in self.reads.keys():
                #compare_hits(self.reads[read_id], hit_start, hit_end, _hit_list, biscore_range_cutoff, length_cutoff, self.project.get_fastq1_readcount(self.sample)) # here should be all the magic
                compare_hits_naive(self.reads[read_id], hit_start, hit_end, _hit_list, biscore_range_cutoff, length_cutoff, self.project.get_fastq1_readcount(self.sample)) # here should be all the magic
            else:
                print ('Read not found: ', read_id)

    
    def import_fastq(self):
        fastq_file = self.project.get_fastq_path(self.sample,self.end)
        line_counter = 0
        current_read = None
        fh = None
        if fastq_file.endswith('.gz'):
            fh = gzip.open(fastq_file, 'rb')
        else:
            fh = open(fastq_file, 'rb')
        if fh:
            for line in fh:
                line_counter += 1
                if line_counter == 5:
                    line_counter = 1
                line = line.decode('utf8').rstrip('\n\r')
                if line_counter == 1:
                    (read_id, end) = self.parse_fastq_seqid(line)
                    
                    if read_id in self.reads:
                        current_read = read_id
                        self.reads[current_read].set_read_id_line(line)
                    else: 
                        current_read = None
                elif line_counter == 2:
                    if current_read:
                        self.reads[current_read].set_sequence(line)
                elif line_counter == 3:
                    if current_read:
                        self.reads[current_read].set_line3(line)
                elif line_counter == 4:
                    if current_read:
                        self.reads[current_read].set_quality(line)
        fh.close()
        
    def import_fasta(self):
        fasta_file = self.project.get_fastq_path(self.sample,self.end)
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
                        self.reads[current_id[1:]].set_read_id_line(current_id)
                        self.reads[current_id[1:]].set_sequence(''.join(sequence))
                    sequence = []
                    seq_id = line[1:]
                    if seq_id in self.reads:
                        current_id = line
                    else: 
                        current_id = None
                else:
                    if current_id:
                        sequence.append(line)
            if current_id:
                self.reads[seq_id].set_read_id_line(current_id)
                self.reads[seq_id].set_sequence(''.join(sequence))
            fh.close()

    def export_read_fastq(self):
        outdir = self.project.get_project_dir(self.sample)
        with open(os.path.join(outdir, self.sample + '_' + self.end + '_' + self.project.get_reads_fastq_name()), 'w') as of:
            for read_id in sorted(self.reads.keys()):
                if self.reads[read_id].get_status() == 'function' or self.reads[read_id].get_status() == 'function,besthit':
                    of.write(self.reads[read_id].get_read_id_line() + '\n')
                    of.write(self.reads[read_id].get_sequence() + '\n') 
                    of.write(self.reads[read_id].get_line3() + '\n') 
                    of.write(self.reads[read_id].get_quality() + '\n') 

    def export_read_fasta(self):
        outdir = self.project.get_project_dir(self.sample)
        with open(os.path.join(outdir, self.sample + '_' + self.end + '_' + self.project.get_reads_fastq_name()), 'w') as of:
            for read_id in sorted(self.reads.keys()):
                if self.reads[read_id].get_status() == 'function' or self.reads[read_id].get_status() == 'function,besthit':
                    of.write(self.reads[read_id].get_read_id_line() + '\n')
                    of.write(self.reads[read_id].get_sequence() + '\n') 

    def export_hit_fastq(self):
        outdir = self.project.get_project_dir(self.sample)
        with open(os.path.join(outdir, self.sample + '_' + self.end + '_' + self.project.get_ref_hits_fastq_name()), 'w') as of:
            for read_id in self.reads.keys():
                for hit in self.reads[read_id].get_hit_list().get_hits():
                    start = hit.get_query_start()
                    end = hit.get_query_end()
                    of.write("@" + self.reads[read_id].get_read_id() + '|' + \
                        str(start) + '|' + str(end) + '\n')
                    if start < end:
                        # hit on + strand
                        start = start - 1
                        end= end
                    else:
                        # hit on - strand
                        t = start
                        start = end - 1
                        end = t
                    try:
                        of.write(self.reads[read_id].get_sequence()[start:end] + '\n') 
                        of.write(self.reads[read_id].get_line3() + '\n') 
                        of.write(self.reads[read_id].get_quality()[start:end] + '\n') 
                    except TypeError:
                        print ('TypeError occurred while exporting ', read_id)

    def export_hit_fasta(self):
        outdir = self.project.get_project_dir(self.sample)
        with open(os.path.join(outdir, self.sample + '_' + self.end + '_' + self.project.get_ref_hits_fastq_name()), 'w') as of:
            for read_id in self.reads.keys():
                for hit in self.reads[read_id].get_hit_list().get_hits():
                    start = hit.get_query_start()
                    end = hit.get_query_end()
                    of.write(">" + self.reads[read_id].get_read_id() + '|' + \
                        str(start) + '|' + str(end) + '\n')
                    if start < end:
                        # hit on + strand
                        start = start - 1
                        end= end
                    else:
                        # hit on - strand
                        t = start
                        start = end - 1
                        end = t
                    try:
                        of.write(self.reads[read_id].get_sequence()[start:end] + '\n') 
                    except TypeError:
                        print ('TypeError occurred while exporting ', read_id)

    def export_hit_list(self):
        outfile = os.path.join(self.project.get_project_dir(self.sample), self.sample + '_' + self.end + '_' + self.project.get_ref_hits_list_name())
        with open(outfile, 'w') as of:
            for read in self.reads.keys():
                for hit in self.reads[read].get_hit_list().get_hits():
                    of.write(str(hit) + '\n')
    
    def get_project(self):
        return self.project

    def get_config(self):
        return self.config

    def get_reads(self):
        return self.reads

    def set_reads(self, reads):
        self.reads = reads
    
    def parse_fastq_seqid(self,line):
        # This function returns read id and read end (if available) for different versions of FASTQ
        if line.startswith('@'):
            line = line[1:]
        if ' ' in line:
            line_tokens = line.split(' ')
            if len(line_tokens) == 2:
                # Casava 1.8+ format
                end = line_tokens[1]
                end = end[0]
                return (line_tokens[0],end)
            elif len(line_tokens) == 3:
                # SRA format?
                if line_tokens[0].endswith('.1') or line_tokens[0].endswith('.2'):
                    # SRA format
                    return (line_tokens[0][:-2], line_tokens[0][-1])
                else:
                    # unknown format
                    return (line_tokens[0], '')
            else:
                # unknown format
                return (line_tokens[0], '')
            # return (line.split('\s')[0], line.split('\s')[1][0])
        elif line.endswith('/1') or line.endswith('/2'):
            # Old Ilumina format
            return (line[:-2], line[-1])
        elif line.endswith('.1') or line.endswith('.2'):
            # Converted SRA
            return (line[:-2], line[-1])
        else:
            return (line, '')
    
    def export_paired_end_reads_fastq(self):
        # This function not only exports FASTQ data, but also imports sequence data from paired-end FASTQ file
        fastq_file = self.project.get_fastq_path(self.sample,get_paired_end(self.end))
        outdir = self.project.get_project_dir(self.sample)
        read_ids = {}
        for read_id in sorted(self.reads.keys()):
            #read_ids[get_paired_read_id(read_id)] = read_id
            read_ids[read_id] = read_id
        line_counter = 0
        with open(os.path.join(outdir, self.sample + '_' + self.end + '_' + self.project.get_pe_reads_fastq_name()), 'w') as of:
            current_read = None
            fh = None
            if fastq_file.endswith('.gz'):
                fh = gzip.open(fastq_file, 'rb')
            else:
                fh = open(fastq_file, 'rb')
            if fh:
                for line in fh:
                    line_counter += 1
                    if line_counter == 5:
                        line_counter = 1
                    line = line.decode('utf8').rstrip('\n\r')
                    if line_counter == 1:
                        (read_id, end) = self.parse_fastq_seqid(line)
                        if read_id in read_ids:
                            current_read = read_id
                            self.reads[current_read].set_pe_id(line)
                            of.write(line + '\n')
                        else: 
                            current_read = None
                    elif current_read:
                        of.write(line + '\n')
                        if line_counter == 2:
                            self.reads[current_read].set_pe_sequence(line)
                        elif line_counter == 3:
                            self.reads[current_read].set_pe_line3(line)
                        elif line_counter == 4:
                            self.reads[current_read].set_pe_quality(line)
                            
                fh.close()
            of.closed
                

    def import_hit_list(self):
        infile = os.path.join(os.path.join(self.project.get_project_dir(self.sample), self.sample + '_' + self.end + '_'+ self.project.get_ref_hits_list_name()))
        ret_val = {}
        _hit_list = None
        current_read_id = None
        
        with open(infile, 'r', newline='') as f:
            tsvin = csv.reader(f, delimiter='\t')
            for row in tsvin:
                if current_read_id is None:
                    # initialize
                    current_read_id = row[0]
                    _hit_list = DiamondHitList(current_read_id)
                elif current_read_id != row[0]:
                    ret_val[current_read_id] = AnnotatedRead(current_read_id)
                    ret_val[current_read_id].set_hit_list(_hit_list)
                    current_read_id = row[0]
                    _hit_list = DiamondHitList(current_read_id)
                hit = DiamondHit()
                hit.import_hit(row)
                _hit_list.add_hit(hit)
            ret_val[current_read_id] = AnnotatedRead(current_read_id)
            ret_val[current_read_id].set_hit_list(_hit_list)
        return ret_val

