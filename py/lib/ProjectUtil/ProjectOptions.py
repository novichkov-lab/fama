import configparser

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

    def get_collection (self, sample):
        if self.project[sample]['collection']:
            return self.project[sample]['collection']
        else:
            return self.project['DEFAULT']['collection']

    def get_project_dir(self, sample):
        if self.project[sample]['sample_dir']:
            return self.project[sample]['sample_dir']
        else:
            return self.project['DEFAULT']['work_dir']
        
    def get_fastq1_path(self, sample):
        if self.project[sample]['fastq_pe1']:
            return self.project[sample]['fastq_pe1']
        else:
            return None
    
    def get_fastq2_path(self, sample):
        if self.project[sample]['fastq_pe2']:
            return self.project[sample]['fastq_pe2']
        else:
            return None
    
    def get_fastq1_readcount(self, sample):
        return self.project[sample]['fastq_pe1_readcount']
    
    def get_fastq2_readcount(self, sample):
        return self.project[sample]['fastq_pe2_readcount']
    
    def get_ref_output_name(self, sample):
        if self.project[sample]['ref_output_name']:
            return self.project[sample]['ref_output_name']
        else:
            return self.project['DEFAULT']['ref_output_name']

    def get_ref_output_name(self):
        return self.project['DEFAULT']['ref_output_name']

    def get_background_output_name(self, sample):
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

    def get_ref_hits_fastq1_name(self, sample):
        return self.project['DEFAULT']['fastq_pe1']

    def get_ref_hits_fastq2_name(self, sample):
        return self.project['DEFAULT']['fastq_pe2']

    def get_report_name(self):
        return self.project['DEFAULT']['report_name']
