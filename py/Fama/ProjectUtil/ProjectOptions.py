import os,configparser

class ProjectOptions:

    def __new__(cls, project_file):
        if not hasattr(cls, 'instance'):
            cls.instance = super(ProjectOptions, cls).__new__(cls)
        cls.instance.parser = configparser.ConfigParser()
        return cls.instance

    def __init__(self, project_file):
        self.parser.read(project_file)
    
    def load_project(self, project_file):
        self.parser.read(project_file)

    def get_collection (self, sample=None):
        if sample:
            if 'collection' in self.parser[sample]:
                return self.parser[sample]['collection']
            else:
                return self.parser['DEFAULT']['collection']
        else:
            return self.parser['DEFAULT']['collection']

    def get_name(self):
        return self.parser['DEFAULT']['project_name']

    def get_work_dir(self):
        return self.parser['DEFAULT']['work_dir']
    
    def get_sample_name(self, sample):
        return self.parser[sample]['sample_id']
    
    def list_samples(self):
        return self.parser.sections()

    def get_fastq1_readcount(self, sample):
        return int(self.parser[sample]['fastq_pe1_readcount'])
    
    def get_fastq2_readcount(self, sample):
        return int(self.parser[sample]['fastq_pe2_readcount'])
    
    # Directory names
    
    def get_project_dir(self, sample):
        if self.parser[sample]['sample_dir']:
            return self.parser[sample]['sample_dir']
        else:
            return self.parser['DEFAULT']['work_dir']

    def get_output_subdir(self, sample):
        if self.parser[sample]['output_subdir']:
            return self.parser[sample]['output_subdir']
        else:
            return self.parser['DEFAULT']['output_subdir']
    
    def get_assembly_dir(self):
        return os.path.join(self.get_work_dir(), self.parser['DEFAULT']['assembly_subdir'])
    
    # File paths
    
    def get_fastq_path(self, sample, end):
        if end == 'pe1':
            if self.parser[sample]['fastq_pe1']:
                return self.parser[sample]['fastq_pe1']
            else:
                return None
        elif end == 'pe2':
            if self.parser[sample]['fastq_pe2']:
                return self.parser[sample]['fastq_pe2']
            else:
                return None
    
    def get_coverage_path(self, sample):
        if self.parser[sample]['coverage']:
            return self.parser[sample]['coverage']
        else:
            return None

    # File names
    
    def get_ref_output_name(self, sample = None):
        if not sample:
            return self.parser['DEFAULT']['ref_output_name']
        elif self.parser[sample]['ref_output_name']:
            return self.parser[sample]['ref_output_name']
        else:
            return self.parser['DEFAULT']['ref_output_name']

    def get_background_output_name(self):
        if self.parser[sample]['background_output_name']:
            return self.parser[sample]['background_output_name']
        else:
            return self.parser['DEFAULT']['background_output_name']

    def get_background_output_name(self):
        return self.parser['DEFAULT']['background_output_name']

    def get_ref_hits_list_name(self, sample):
        if self.parser[sample]['ref_hits_list_name']:
            return self.parser[sample]['ref_hits_list_name']
        else:
            return self.parser['DEFAULT']['ref_hits_list_name']

    def get_ref_hits_list_name(self):
        return self.parser['DEFAULT']['ref_hits_list_name']

    def get_ref_hits_fastq_name(self):
        return self.parser['DEFAULT']['ref_hits_fastq_name']

    def get_reads_fastq_name(self):
        return self.parser['DEFAULT']['reads_fastq_name']
    
    def get_pe_reads_fastq_name(self):
        return self.parser['DEFAULT']['pe_reads_fastq_name']

    def get_report_name(self):
        return self.parser['DEFAULT']['report_name']

    def get_xml_name(self):
        return self.parser['DEFAULT']['xml_name']

    def get_html_name(self):
        return self.parser['DEFAULT']['html_name']

    def get_reads_json_name(self):
        return self.parser['DEFAULT']['reads_json_name']

    def get_rpkg_scaling_factor(self, sample):
        if self.parser[sample]['rpkg_scaling']:
            return float(self.parser[sample]['rpkg_scaling'])
        else:
            return None
    
    def set_sample_data(self,sample):
        self.parser[sample.sample_id]['sample_id'] = sample.sample_name
        self.parser[sample.sample_id]['fastq_pe1'] = sample.fastq_fwd_path
        self.parser[sample.sample_id]['fastq_pe2'] = sample.fastq_rev_path
        self.parser[sample.sample_id]['fastq_pe1_readcount'] = str(sample.fastq_fwd_readcount)
        self.parser[sample.sample_id]['fastq_pe2_readcount'] = str(sample.fastq_rev_readcount)
        self.parser[sample.sample_id]['fastq_pe1_basecount'] = str(sample.fastq_fwd_basecount)
        self.parser[sample.sample_id]['fastq_pe2_basecount'] = str(sample.fastq_rev_basecount)
        self.parser[sample.sample_id]['sample_dir'] = sample.work_directory
        self.parser[sample.sample_id]['rpkg_scaling'] = str(sample.rpkg_scaling_factor)
        self.parser[sample.sample_id]['replicate'] = sample.replicate
        
    def save_options(self, project_file):
        with open (project_file, 'w') as f:
            self.parser.write(f)
            f.closed
        
