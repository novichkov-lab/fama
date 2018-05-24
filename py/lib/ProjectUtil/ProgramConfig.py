import configparser

class ProgramConfig:

    def __new__(cls, config_file):
        if not hasattr(cls, 'instance'):
            cls.instance = super(ProgramConfig, cls).__new__(cls)
        cls.instance.config = configparser.ConfigParser()
        return cls.instance

    def __init__(self, config_file):
        self.config.read(config_file)
    
    def load_config(self, config_file):
        self.config.read(config_file)

    def get_functions_file (self, collection):
        if self.config[collection]['functions_file']:
            return self.config[collection]['functions_file']
        else:
            return self.config['DEFAULT']['functions_file']

    def get_proteins_list_file (self, collection):
        if self.config[collection]['proteins_list_file']:
            return self.config[collection]['proteins_list_file']
        else:
            return self.config['DEFAULT']['proteins_list_file']

    def get_identity_cutoff(self, collection):
        if self.config[collection]['identity_cutoff']:
            return float(self.config[collection]['identity_cutoff'])
        else:
            return float(self.config['DEFAULT']['identity_cutoff'])

    def get_length_cutoff(self, collection):
        if self.config[collection]['length_cutoff']:
            return int(self.config[collection]['length_cutoff'])
        else:
            return int(self.config['DEFAULT']['length_cutoff'])

    def get_overlap_cutoff(self, collection):
        if self.config[collection]['hits_overlap_cutoff']:
            return int(self.config[collection]['hits_overlap_cutoff'])
        else:
            return int(self.config['DEFAULT']['hits_overlap_cutoff'])

    def get_biscore_range_cutoff(self, collection):
        if self.config[collection]['biscore_range_cutoff']:
            return int(self.config[collection]['biscore_range_cutoff'])
        else:
            return int(self.config['DEFAULT']['biscore_range_cutoff'])

    def list_collections(self):
        return self.config.sections()
    
    
