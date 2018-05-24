import csv, re
import gzip
from lib.ProjectUtil.ProgramConfig import ProgramConfig
from lib.ProjectUtil.ProjectOptions import ProjectOptions
from lib.ReferenceLibrary.ReferenceData import ReferenceData
from lib.DiamondParser.DiamondHit import DiamondHit
from lib.DiamondParser.DiamondHitList import DiamondHitList
from lib.ReadUtil.AnnotatedRead import AnnotatedRead

class DiamondParser:

    def __init__(self, config_file, project_file, sample, end):
        self.reads = {}
        self.sample = sample
        self.end = end
        self.config = ProgramConfig(config_file)
        collections = self.config.list_collections()
        self.project = ProjectOptions(project_file)
        collection = self.project.get_collection(self.sample)
        if collection not in collections:
            raise Error ('Collection ' + self.collection + ' not found')
        self.collection = collection
        self.ref_data = ReferenceData(self.config)
        self.ref_data.load_reference_data(self.collection)

    def parse_reference_output(self):
        
        tsvfile = self.project.get_project_dir(self.sample) + '/' + self.sample + '_' + self.end + '_'+ self.project.get_ref_output_name()
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
                        read.set_hits(_hit_list)
                        self.reads[current_sequence_read_id] = read

                    current_sequence_read_id = hit.get_query_id()
                    _hit_list = DiamondHitList(current_sequence_read_id)
                _hit_list.add_hit(hit)

    def parse_background_output(self):
        
        if not self.reads:
            self.reads = import_hit_list()
        
        tsvfile = self.project.get_project_dir(self.sample) + '/' + self.sample + '_' + self.end + '_'+ self.project.get_background_output_name()
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
                    self.reads[read_id].compare_hits(hit_start, hit_end, _hit_list, biscore_range_cutoff) # here should be all the magic
                    
                    current_query_id = hit.get_query_id()
                    _hit_list = DiamondHitList(current_query_id)
                _hit_list.add_hit(hit)

    
    def import_fastq(self):
        fastq_file = self.project.get_fastq1_path(self.sample)
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
        
    def export_hit_fastq(self):
        outdir = self.project.get_project_dir(self.sample)
        with open(outdir+'reference_hits.fastq', 'w') as of:
            for read_id in self.reads.keys():
                for hit in self.reads[read_id].get_hits().get_hits():
                    start = hit.get_query_start()
                    end = hit.get_query_end()
                    of.write("@" + self.reads[read_id].get_read_id() + '|' + \
                        str(start) + '|' + str(end) + '\n')
                    if start < end:
                        # hit on + strand
                        start = start - 1
                        end= end - 1
                    else:
                        # hit on - strand
                        t = start
                        start = end - 1
                        end = t - 1
                    try:
                        of.write(self.reads[read_id].get_sequence()[start:end] + '\n') 
                        of.write(self.reads[read_id].get_line3() + '\n') 
                        of.write(self.reads[read_id].get_quality()[start:end] + '\n') 
                    except TypeError:
                        print ('TypeError occurred while exporting ', read_id)

    def export_hit_list(self):
        outfile = self.project.get_project_dir(self.sample) + '/' + self.sample + '_' + self.end + '_' + self.project.get_ref_hits_list_name(self.sample)
        with open(outdir+'reference_hits.txt', 'w') as of:
            for read in self.reads.keys():
                for hit in self.reads[read].get_hits().get_hits():
                    of.write(str(hit) + '\n')
    
    def import_hit_list(self.sample):
        ret_val = {}
        _hit_list = None
        current_read_id = None
        tsvfile = self.project.get_project_dir(self.sample) + '/' + self.sample + '_' + self.end + '_'+ self.project.get_ref_hits_list_name(self.sample)
        with open(tsvfile, 'r', newline='') as f:
            tsvin = csv.reader(f, delimiter='\t')
            for row in tsvin:
                if current_read_id is None:
                    # initialize
                    current_read_id = row[0]
                    _hit_list = DiamondHitList(current_read_id)
                elif current_read_id != row[0]:
                    read = AnnotatedRead(current_sequence_read_id)
                    read.set_hits(_hit_list)
                    ret_val[current_read_id] = read
                    current_read_id = row[0]
                    _hit_list = DiamondHitList(current_read_id)
                hit = DiamondHit()
                hit.import_hit(row)
                _hit_list.add_hit(hit)

    def get_reads(self):
        return self.reads
    
    def parse_fastq_seqid(self,line):
        #to be implemented
        return (line.split('\s')[0][1:], '1')
