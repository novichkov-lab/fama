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
    
    def get_threads(self):
        return self.config['DEFAULT']['threads']

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

    def get_evalue_cutoff(self, collection):
        if self.config[collection]['evalue_cutoff']:
            return float(self.config[collection]['evalue_cutoff'])
        else:
            return float(self.config['DEFAULT']['evalue_cutoff'])

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
            return float(self.config[collection]['biscore_range_cutoff'])
        else:
            return float(self.config['DEFAULT']['biscore_range_cutoff'])
            
    def get_reference_diamond_db(self, collection):
        if self.config[collection]['reference_diamond_db']:
            return self.config[collection]['reference_diamond_db']
        else:
            return self.config['DEFAULT']['reference_diamond_db']

    def get_background_diamond_db(self, collection):
        if self.config[collection]['background_diamond_db']:
            return self.config[collection]['background_diamond_db']
        else:
            return self.config['DEFAULT']['background_diamond_db']

    def get_reference_db_size(self, collection):
        if self.config[collection]['reference_db_size']:
            return int(self.config[collection]['reference_db_size'])
        else:
            return int(self.config['DEFAULT']['reference_db_size'])

    def get_background_db_size(self, collection):
        if self.config[collection]['background_db_size']:
            return int(self.config[collection]['background_db_size'])
        else:
            return int(self.config['DEFAULT']['background_db_size'])

    def get_taxonomy_names_file (self):
            return self.config['DEFAULT']['taxonomy_names_file']

    def get_taxonomy_nodes_file (self):
            return self.config['DEFAULT']['taxonomy_nodes_file']

    def get_taxonomy_merged_file (self):
            return self.config['DEFAULT']['taxonomy_merged_file']

    def list_collections(self):
        return self.config.sections()
    
    def get_uniprot_diamond_db(self):
        return self.config['DEFAULT']['uniprot_db_path']
    
    def get_diamond_path(self):
        return self.config['DEFAULT']['aligner_path']

    def get_microbecensus_datadir(self):
        return self.config['DEFAULT']['microbecensus_data']
        
    def get_krona_path(self):
        return self.config['DEFAULT']['krona_path']

    def get_megahit_path(self):
        return self.config['DEFAULT']['megahit_path']

    def get_metaspades_path(self):
        return self.config['DEFAULT']['metaspades_path']

    def get_bowtie_indexer_path(self):
        return self.config['DEFAULT']['bowtie_indexer_path']

    def get_bowtie_path(self):
        return self.config['DEFAULT']['bowtie_path']

    def get_prodigal_path(self):
        return self.config['DEFAULT']['prodigal_path']
