import os
from collections import defaultdict

class Sample(object):

    def __init__(self, sample_id = '', sample_name = None, is_paired_end = True, \
                fastq_fwd_path = None, fastq_rev_path = None, fastq_fwd_readcount = 0, \
                fastq_rev_readcount = 0, rpkm_scaling_factor = 0.0, \
                fastq_fwd_basecount = 0, fastq_rev_basecount = 0, \
                rpkg_scaling_factor = 0.0, fragment_length = None,\
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
        self.fragment_length = fragment_length
        self.reads = defaultdict(dict)
        
    def load_sample(self,options):
        self.sample_name = options.parser.get(self.sample_id,'sample_id')
        self.fastq_fwd_path = options.parser.get(self.sample_id,'fastq_pe1')
        self.fastq_rev_path = options.parser.get(self.sample_id,'fastq_pe2', fallback = None)
        if self.fastq_rev_path is None or self.fastq_rev_path == '':
            self.is_paired_end = False
        else:
            self.is_paired_end = True
        self.fastq_fwd_readcount = options.parser.getint(self.sample_id,'fastq_pe1_readcount')
        self.fastq_rev_readcount = options.parser.getint(self.sample_id,'fastq_pe2_readcount', fallback=0)
        self.fastq_fwd_basecount = options.parser.getint(self.sample_id,'fastq_pe1_basecount', fallback=0)
        self.fastq_rev_basecount = options.parser.getint(self.sample_id,'fastq_pe2_basecount', fallback=0)
        self.work_directory = options.parser.get(self.sample_id,'sample_dir')
        if self.fastq_fwd_readcount > 0:
            self.rpkm_scaling_factor = 1000000/self.fastq_fwd_readcount
        if self.is_paired_end:
            self.fragment_length = options.parser.getint(self.sample_id,'fragment_length')
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
            return self.fastq_fwd_basecount / self.fastq_fwd_readcount
        elif end == 'pe2' and self.fastq_rev_readcount != 0.0:
            return self.fastq_rev_basecount / self.fastq_rev_readcount
        else:
            return 0.0

        
