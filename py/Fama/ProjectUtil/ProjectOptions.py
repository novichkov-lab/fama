import os,configparser

class ProjectOptions:

    def __new__(cls, project_file):
        if not hasattr(cls, 'instance'):
            cls.instance = super(ProjectOptions, cls).__new__(cls)
        cls.instance.parser = configparser.ConfigParser()
        return cls.instance

    def __init__(self, project_file):
        self.project_file = project_file
        self.parser.read(project_file)
    
    # Properties 
    @property
    def project_name(self):
        return self.parser['DEFAULT']['project_name']

    @property
    def work_dir(self):
        return self.parser['DEFAULT']['work_dir']
    
    @property
    def ref_output_name(self):
        return self.parser['DEFAULT']['ref_output_name']

    @property
    def background_output_name(self):
        return self.parser['DEFAULT']['background_output_name']

    @property
    def ref_hits_list_name(self):
        return self.parser['DEFAULT']['ref_hits_list_name']

    @property
    def ref_hits_fastq_name(self):
        return self.parser['DEFAULT']['ref_hits_fastq_name']

    @property
    def reads_fastq_name(self):
        return self.parser['DEFAULT']['reads_fastq_name']
    
    @property
    def pe_reads_fastq_name(self):
        return self.parser['DEFAULT']['pe_reads_fastq_name']

    @property
    def report_name(self):
        return self.parser['DEFAULT']['report_name']

    @property
    def xml_name(self):
        return self.parser['DEFAULT']['xml_name']

    @property
    def html_name(self):
        return self.parser['DEFAULT']['html_name']

    @property
    def reads_json_name(self):
        return self.parser['DEFAULT']['reads_json_name']

    @property
    def assembly_dir(self):
        return os.path.join(self.work_dir(), self.parser['DEFAULT']['assembly_subdir'])
    

    # Methods
    
    def load_project(self, project_file):
        self.parser.read(project_file)

    def get_collection (self, sample=None):
        if sample:
            if self.parser.has_option(sample,'collection'):
                return self.parser[sample]['collection']
            else:
                return self.parser['DEFAULT']['collection']
        else:
            return self.parser['DEFAULT']['collection']

    def get_sample_name(self, sample):
        return self.parser[sample]['sample_id']
    
    def list_samples(self):
        return self.parser.sections()

    def get_fastq1_readcount(self, sample):
        return int(self.parser[sample]['fastq_pe1_readcount'])
    
    def get_fastq2_readcount(self, sample):
        return int(self.parser[sample]['fastq_pe2_readcount'])
    
    def get_project_dir(self, sample):
        if self.parser.has_option(sample,'sample_dir'):
            return self.parser[sample]['sample_dir']
        else:
            return self.parser['DEFAULT']['work_dir']

    def get_output_subdir(self, sample):
        if self.parser.has_option(sample,'output_subdir'):
            return self.parser[sample]['output_subdir']
        else:
            return self.parser['DEFAULT']['output_subdir']
    
    def get_fastq_path(self, sample, end):
        if end == 'pe1':
            if self.parser.has_option(sample,'fastq_pe1'):
                return self.parser[sample]['fastq_pe1']
            else:
                return None
        elif end == 'pe2':
            if self.parser.has_option(sample,'fastq_pe2'):
                return self.parser[sample]['fastq_pe2']
            else:
                return None
    
    def get_coverage_path(self, sample):
        if self.parser.has_option(sample, 'coverage'):
            return self.parser[sample]['coverage']
        elif self.parser.has_option('', 'coverage'):
            return self.parser['DEFAULT']['coverage']
        else:
            return None

    def get_rpkg_scaling_factor(self, sample):
        if self.parser.has_option(sample,'rpkg_scaling'):
            return float(self.parser[sample]['rpkg_scaling'])
        else:
            return None

    def set_sample_data(self,sample):
        # sets parameters of underlying configParser. Call for each sample 
        # before calling save_options()
        self.parser[sample.sample_id]['sample_id'] = sample.sample_name
        self.parser[sample.sample_id]['fastq_pe1'] = sample.fastq_fwd_path
        self.parser[sample.sample_id]['fastq_pe1_readcount'] = str(sample.fastq_fwd_readcount)
        self.parser[sample.sample_id]['fastq_pe1_basecount'] = str(sample.fastq_fwd_basecount)
        if sample.is_paired_end:
            self.parser[sample.sample_id]['fastq_pe2'] = sample.fastq_rev_path
            self.parser[sample.sample_id]['fastq_pe2_readcount'] = str(sample.fastq_rev_readcount)
            self.parser[sample.sample_id]['fastq_pe2_basecount'] = str(sample.fastq_rev_basecount)
        self.parser[sample.sample_id]['sample_dir'] = sample.work_directory
        if sample.rpkg_scaling_factor is not None:
            self.parser[sample.sample_id]['rpkg_scaling'] = str(sample.rpkg_scaling_factor)
        self.parser[sample.sample_id]['replicate'] = sample.replicate
        if sample.insert_size is not None:
            self.parser[sample.sample_id]['insert_size'] = str(sample.insert_size)
        
    def save_options(self):
        i = 0
        while True:
            if os.path.exists(self.project_file + '.old.' + str(i)):
                i += 1
            else:
                backup_file = self.project_file + '.old.' + str(i)
                break
        os.rename(self.project_file, backup_file)
        with open (self.project_file, 'w') as f:
            self.parser.write(f)
            f.closed
        print('Old version of the project file copied to',backup_file)
        
