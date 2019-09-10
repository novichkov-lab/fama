import configparser

class ProgramConfig(object):
    """Configuration parameters for Fama run. Wraps a ConfigParser instance,
    which reads config ini file.
    
    This class is intended to be instantiated only once.
    
    """

    def __new__(cls, config_file):
        if not hasattr(cls, 'instance'):
            cls.instance = super(ProgramConfig, cls).__new__(cls)
        cls.instance.config = configparser.ConfigParser()
        return cls.instance

    def __init__(self, config_file):
        """ 
        Args:
            config_file (str): path to config ini file
        """
        self.config.read(config_file)
    
    def load_config(self, config_file):
        """Loads(reloads) config ini file

        Args:
            config_file (str): path to config ini file
        """
        self.config.read(config_file)
    
    @property
    def threads(self):
        """Number of threads available for multi-threading programs"""
        return self.config['DEFAULT']['threads']

    def get_functions_file (self, collection):
        """Returns path to file with list of functions for a collection.
        If collection has no link to such file, returns default path.
        
        Args:
            collection (str): identifier of collection
            
        """
        if self.config.has_option(collection,'functions_file'):
            return self.config[collection]['functions_file']
        else:
            return self.config['DEFAULT']['functions_file']

    def get_protein_list_file (self, collection):
        """Returns path to file with list of proteins for a collection.
        If collection has no link to such file, returns default path.
        
        Args:
            collection (str): identifier of collection
            
        """
        if self.config.has_option(collection,'proteins_list_file'):
            return self.config[collection]['proteins_list_file']
        else:
            return self.config['DEFAULT']['proteins_list_file']

    def get_identity_cutoff(self, collection):
        """Returns minimal amino acid identity % for DIAMOND hits.
        
        Args:
            collection (str): identifier of collection
            
        """
        if self.config.has_option(collection,'identity_cutoff'):
            return float(self.config[collection]['identity_cutoff'])
        else:
            return float(self.config['DEFAULT']['identity_cutoff'])

    def get_evalue_cutoff(self, collection):
        """Returns e-value threshold for DIAMOND search.
        
        Args:
            collection (str): identifier of collection
            
        """
        if self.config.has_option(collection,'evalue_cutoff'):
            return float(self.config[collection]['evalue_cutoff'])
        else:
            return float(self.config['DEFAULT']['evalue_cutoff'])

    def get_length_cutoff(self, collection):
        """Returns minimal length of amino acid alignment of DIAMOND hits.
        
        Args:
            collection (str): identifier of collection
            
        """
        if self.config.has_option(collection,'length_cutoff'):
            return int(self.config[collection]['length_cutoff'])
        else:
            return int(self.config['DEFAULT']['length_cutoff'])

    def get_overlap_cutoff(self, collection):
        """Returns maximal length of overlapping part for two alignments 
        to be considered the same hit.
        
        Args:
            collection (str): identifier of collection
            
        """
        if self.config.has_option(collection,'hits_overlap_cutoff'):
            return int(self.config[collection]['hits_overlap_cutoff'])
        else:
            return int(self.config['DEFAULT']['hits_overlap_cutoff'])

    def get_biscore_range_cutoff(self, collection):
        """Returns lowest acceptaple relative bitscore (relative to top hit bit-score)
        
        Args:
            collection (str): identifier of collection
            
        """
        if self.config.has_option(collection,'biscore_range_cutoff'):
            return float(self.config[collection]['biscore_range_cutoff'])
        else:
            return float(self.config['DEFAULT']['biscore_range_cutoff'])
            
    def get_ranks_cutoffs(self, collection):
        """Returns dictionary of amino acid identity % values for each taxonomy level
        
        Args:
            collection (str): identifier of collection
            
        """
        ret_val = {}
        if self.config.has_option(collection,'rank_cutoffs'):
            cutoffs = self.config[collection]['rank_cutoffs']
            vals = cutoffs.split(',')
            for k_v in zip(vals[::2],vals[1::2]):
                ret_val[k_v[0]] = float(k_v[1])
            ret_val['norank'] = self.identity_cutoff(collection)
        else:
            pass
        return ret_val
        
    def get_reference_diamond_db(self, collection):
        """Returns path to DIAMOND database for reference search
        
        Args:
            collection (str): identifier of collection
            
        """
        if self.config.has_option(collection,'reference_diamond_db'):
            return self.config[collection]['reference_diamond_db']
        else:
            return self.config['DEFAULT']['reference_diamond_db']

    def get_background_diamond_db(self, collection):
        """Returns path to DIAMOND database for background search
        
        Args:
            collection (str): identifier of collection
            
        """
        if self.config.has_option(collection,'background_diamond_db'):
            return self.config[collection]['background_diamond_db']
        else:
            return self.config['DEFAULT']['background_diamond_db']

    def get_reference_db_size(self, collection):
        """Returns size of DIAMOND database for reference search
        
        Args:
            collection (str): identifier of collection
            
        """
        if self.config.has_option(collection,'reference_db_size'):
            return int(self.config[collection]['reference_db_size'])
        else:
            return int(self.config['DEFAULT']['reference_db_size'])

    def get_background_db_size(self, collection):
        """Returns size of DIAMOND database for background search
        
        Args:
            collection (str): identifier of collection
            
        """
        if self.config.has_option(collection,'background_db_size'):
            return int(self.config[collection]['background_db_size'])
        else:
            return int(self.config['DEFAULT']['background_db_size'])

    @property
    def taxonomy_names_file (self):
        """Path to NCBI taxonomy names.dmp file"""
        return self.config['DEFAULT']['taxonomy_names_file']

    @property
    def taxonomy_nodes_file (self):
        """Path to NCBI taxonomy nodes.dmp file"""
        return self.config['DEFAULT']['taxonomy_nodes_file']

    @property
    def taxonomy_merged_file (self):
        """Path to NCBI taxonomy merged.dmp file"""
        return self.config['DEFAULT']['taxonomy_merged_file']

    @property
    def collections(self):
        """List of collection names"""
        return self.config.sections()
    
    @property
    def uniprot_diamond_db(self):
        """Path to Uniprot DIAMOND database"""
        return self.config['DEFAULT']['uniprot_db_path']
    
    @property
    def diamond_path(self):
        """Path to DIAMOND binary"""
        return self.config['DEFAULT']['aligner_path']

    @property
    def microbecensus_datadir(self):
        """Path to MicrobeCensus data directory"""
        return self.config['DEFAULT']['microbecensus_data']
        
    @property
    def krona_path(self):
        """Path to KronaTools"""
        return self.config['DEFAULT']['krona_path']

    @property
    def megahit_path(self):
        """Path to MEGAHIT binary"""
        return self.config['DEFAULT']['megahit_path']

    @property
    def metaspades_path(self):
        """Path to metaSPAdes binary"""
        return self.config['DEFAULT']['metaspades_path']

    @property
    def bowtie_indexer_path(self):
        """Path to Bowtie indexer binary"""
        return self.config['DEFAULT']['bowtie_indexer_path']

    @property
    def bowtie_path(self):
        """Path to Bowtie aligner binary"""
        return self.config['DEFAULT']['bowtie_path']

    @property
    def prodigal_path(self):
        """Path to Prodigal binary"""
        return self.config['DEFAULT']['prodigal_path']
