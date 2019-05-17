import os
from collections import defaultdict

from Fama.lib_est import get_lib_est

class Sample(object):

    def __init__(self, sample_id = '', sample_name = None, is_paired_end = True, \
                fastq_fwd_path = None, fastq_rev_path = None, fastq_fwd_readcount = 0, \
                fastq_rev_readcount = 0, rpkm_scaling_factor = 0.0, \
                fastq_fwd_basecount = 0, fastq_rev_basecount = 0, \
                rpkg_scaling_factor = 0.0, insert_size = None,\
                work_directory = None, replicate = '0'):
        self.sample_id = sample_id
        self.sample_name = sample_name
        self.is_paired_end = is_paired_end
        self.fastq_fwd_path = fastq_fwd_path
        self.fastq_rev_path = fastq_rev_path
        self.fastq_fwd_readcount = fastq_fwd_readcount
        self.fastq_rev_readcount = fastq_rev_readcount
        self.fastq_fwd_basecount = 0
        self.fastq_rev_basecount = 0
        self.work_directory = work_directory
        self.rpkm_scaling_factor = rpkm_scaling_factor
        self.rpkg_scaling_factor = rpkg_scaling_factor
        self.replicate = replicate
        self.insert_size = insert_size
        self.reads = defaultdict(dict)
        
    def load_sample(self,options):
        self.sample_name = options.parser.get(self.sample_id,'sample_id')
        self.fastq_fwd_path = options.parser.get(self.sample_id,'fastq_pe1')
        self.fastq_rev_path = options.parser.get(self.sample_id,'fastq_pe2', fallback = None)
        if self.fastq_rev_path is None or self.fastq_rev_path == '':
            self.is_paired_end = False
        else:
            self.is_paired_end = True
        self.fastq_fwd_readcount = options.parser.getint(self.sample_id,'fastq_pe1_readcount', fallback=0)
        self.fastq_rev_readcount = options.parser.getint(self.sample_id,'fastq_pe2_readcount', fallback=0)
        self.fastq_fwd_basecount = options.parser.getint(self.sample_id,'fastq_pe1_basecount', fallback=0)
        self.fastq_rev_basecount = options.parser.getint(self.sample_id,'fastq_pe2_basecount', fallback=0)
        self.work_directory = options.parser.get(self.sample_id,'sample_dir')
        if self.fastq_fwd_readcount > 0:
            self.rpkm_scaling_factor = 1000000/self.fastq_fwd_readcount
        if self.is_paired_end:
            self.insert_size = options.parser.getfloat(self.sample_id,'insert_size', fallback=0)
        self.rpkg_scaling_factor = options.parser.getfloat(self.sample_id,'rpkg_scaling', fallback=0.0)
        self.replicate = options.parser.get(self.sample_id,'replicate')
        
    def import_rpkg_scaling_factor(self):
        mc_outfile = os.path.join(self.work_directory, 'microbecensus.out.txt')
        if os.path.exists(mc_outfile):
            with open(mc_outfile, 'r') as f:
                for line in f:
                    if line.startswith('average_genome_size:'):
                        avg = float(line.rstrip('\n\r').split('\t')[-1])
                        if self.fastq_fwd_basecount > 0:
                            self.rpkg_scaling_factor = avg/self.fastq_fwd_basecount
                        elif self.fastq_rev_basecount > 0:
                            self.rpkg_scaling_factor = avg/self.fastq_rev_basecount
                        else:
                            raise ValueError('Total bases count required for RPKG normalization')
                f.closed
        else:
            return
    
    def get_avg_read_length(self,end):
        if end == 'pe1':
            if self.fastq_fwd_readcount == 0:
                return 0.0
            else:
                return self.fastq_fwd_basecount / self.fastq_fwd_readcount
        elif end == 'pe2' and self.fastq_rev_readcount != 0.0:
            return self.fastq_rev_basecount / self.fastq_rev_readcount
        else:
            return 0.0

    def estimate_average_insert_size(self, alignment_length_threshold):
        if not self.is_paired_end:
            return None
        if len(self.reads['pe1']) == 0:
            return None
        if len(self.reads['pe2']) == 0:
            return None
        gene_length_threshold = self.get_avg_read_length('pe1')
        print('gene_length_threshold', gene_length_threshold)
        print('alignment_length_threshold',alignment_length_threshold)
        read_data = defaultdict(dict)
        print ('pe1 reads', str(len(self.reads['pe1'])))
        print ('pe2 reads', str(len(self.reads['pe2'])))
        for read_id,read1 in self.reads['pe1'].items():
            if read1.get_status() != 'function':
                continue
            if read_id not in self.reads['pe2']:
                continue
#            print ('Found read with two mapped ends')
            read2 =self.reads['pe2'][read_id]
            if read2.get_status() != 'function':
                continue
            for hit in read1.get_hit_list().get_hits():
                if hit.get_subject_id() not in [h.get_subject_id() for h in read2.get_hit_list().get_hits()]:
#                    print ('Different target proteins: skipped')
                    continue
                if hit.s_len*3 < gene_length_threshold:
#                    print ('Target protein shorter than threshold: skipped')
                    continue
                if hit.s_end - hit.s_start < alignment_length_threshold:
                    continue
                for hit2 in read2.get_hit_list().get_hits():
                    if hit.get_subject_id() != hit2.get_subject_id():
                        continue
#                    print ('Found read with two hits in one protein')
                    if hit2.s_end - hit2.s_start < alignment_length_threshold:
                        continue
#                    print ('Found read with two hits in one protein longer than alignment cutoff')
                    if (hit.s_end - hit2.s_start) > (hit2.s_end - hit.s_start):
                        # Do not count overhangs
                        insert_size = 3 * (hit.s_end - hit2.s_start)
                        # Count overhangs
                        #insert_size = 3 * (hit.s_end - hit2.s_start) + hit2.q_start - 1 + len(read1.sequence)  - hit.q_end
                    else:
                        # Do not count overhangs
                        insert_size = 3 * (hit2.s_end - hit.s_start)
                        # Count overhangs
                        #insert_size = 3 * (hit2.s_end - hit.s_start) + hit.q_start - 1 + len(read2.sequence)  - hit2.q_end
                    if insert_size > alignment_length_threshold * 3:
                        read_data[read_id]['tlen'] = insert_size
                        read_data[read_id]['rlen'] = (len(read1.sequence) + len(read2.sequence)) / 2
                        read_data[read_id]['ref_len'] = hit.s_len*3
                        read_data[read_id]['ref_name'] = hit.get_subject_id()
                    break
        print(str(len(read_data)), 'fragments found')
        avg_insert_size = get_lib_est(read_data, self.work_directory)
        self.insert_size = avg_insert_size
        return avg_insert_size

