import os,configparser

class ProjectOptions:

    def __new__(cls, project_file):
        if not hasattr(cls, 'instance'):
            cls.instance = super(ProjectOptions, cls).__new__(cls)
        cls.instance.project = configparser.ConfigParser()
        return cls.instance

    def __init__(self, project_file):
        self.project.read(project_file)
    
    def load_project(self, project_file):
        self.project.read(project_file)

    def get_collection (self, sample=None):
        if sample:
            if 'collection' in self.project[sample]:
                return self.project[sample]['collection']
            else:
                return self.project['DEFAULT']['collection']
        else:
            return self.project['DEFAULT']['collection']

    def get_name(self):
        return self.project['DEFAULT']['project_name']

    def get_work_dir(self):
        return self.project['DEFAULT']['work_dir']
    
    def get_sample_id(self, sample):
        return self.project[sample]['sample_id']
    
    def list_samples(self):
        return self.project.sections()

    def get_fastq1_readcount(self, sample):
        return int(self.project[sample]['fastq_pe1_readcount'])
    
    def get_fastq2_readcount(self, sample):
        return int(self.project[sample]['fastq_pe2_readcount'])
    
    # Directory names
    
    def get_project_dir(self, sample):
        if self.project[sample]['sample_dir']:
            return self.project[sample]['sample_dir']
        else:
            return self.project['DEFAULT']['work_dir']

    def get_output_subdir(self, sample):
        if self.project[sample]['output_subdir']:
            return self.project[sample]['output_subdir']
        else:
            return self.project['DEFAULT']['output_subdir']
    
    def get_assembly_dir(self):
        return os.path.join(self.get_work_dir(), self.project['DEFAULT']['assembly_subdir'])
    
    # File paths
    
    def get_fastq_path(self, sample, end):
        if end == 'pe1':
            if self.project[sample]['fastq_pe1']:
                return self.project[sample]['fastq_pe1']
            else:
                return None
        elif end == 'pe2':
            if self.project[sample]['fastq_pe2']:
                return self.project[sample]['fastq_pe2']
            else:
                return None
    
    def get_coverage_path(self, sample):
        if self.project[sample]['coverage']:
            return self.project[sample]['coverage']
        else:
            return None

    # File names
    
    def get_ref_output_name(self, sample = None):
        if not sample:
            return self.project['DEFAULT']['ref_output_name']
        elif self.project[sample]['ref_output_name']:
            return self.project[sample]['ref_output_name']
        else:
            return self.project['DEFAULT']['ref_output_name']

    def get_background_output_name(self):
        if self.project[sample]['background_output_name']:
            return self.project[sample]['background_output_name']
        else:
            return self.project['DEFAULT']['background_output_name']

    def get_background_output_name(self):
        return self.project['DEFAULT']['background_output_name']

    def get_ref_hits_list_name(self, sample):
        if self.project[sample]['ref_hits_list_name']:
            return self.project[sample]['ref_hits_list_name']
        else:
            return self.project['DEFAULT']['ref_hits_list_name']

    def get_ref_hits_list_name(self):
        return self.project['DEFAULT']['ref_hits_list_name']

    def get_ref_hits_fastq_name(self):
        return self.project['DEFAULT']['ref_hits_fastq_name']

    def get_reads_fastq_name(self):
        return self.project['DEFAULT']['reads_fastq_name']
    
    def get_pe_reads_fastq_name(self):
        return self.project['DEFAULT']['pe_reads_fastq_name']

    def get_report_name(self):
        return self.project['DEFAULT']['report_name']

    def get_xml_name(self):
        return self.project['DEFAULT']['xml_name']

    def get_html_name(self):
        return self.project['DEFAULT']['html_name']

    def get_reads_json_name(self):
        return self.project['DEFAULT']['reads_json_name']

