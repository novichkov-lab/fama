"""Describes ProgramConfig class"""
import configparser
from lib.utils.utils import singleton


@singleton
class ProgramConfig(object):
    """Configuration parameters for Fama run. Wraps a ConfigParser instance,
    which reads config ini file. This is a singleton class.

    Attributes:
        parser (:obj:ConfigParser): underlying ConfigParser instance
    """

    def __init__(self, config_file):
        """
        Args:
            config_file (str): path to config ini file
        """
        self.parser = configparser.ConfigParser()
        self.parser.read(config_file)

    @property
    def threads(self):
        """Number of threads available for multi-threading programs"""
        return self.parser['DEFAULT']['threads']

    def get_functions_file(self, collection):
        """Returns path to file with list of functions for a collection.
        If collection has no link to such file, returns default path.

        Args:
            collection (str): identifier of collection
        """
        try:
            result = self.parser[collection]['functions_file']
        except KeyError:
            result = self.parser['DEFAULT']['functions_file']
        return result

    def get_protein_list_file(self, collection):
        """Returns path to file with list of proteins for a collection.
        If collection has no link to such file, returns default path.

        Args:
            collection (str): identifier of collection
        """
        try:
            result = self.parser[collection]['proteins_list_file']
        except KeyError:
            result = self.parser['DEFAULT']['proteins_list_file']
        return result

    def get_taxonomy_file(self, collection):
        """Returns path to file with list of taxonomy identifiers for collection.
        If collection has no link to such file, returns default path.

        Args:
            collection (str): identifier of collection
        """
        try:
            result = self.parser[collection]['taxonomy_file']
        except KeyError:
            print('Taxonomy file for collection' + collection + 
                  ' not found. Loading entire NCBI taxonomy instead.')
            result = self.parser['DEFAULT']['taxonomy_file']
        return result

    def get_identity_cutoff(self, collection):
        """Returns minimal amino acid identity % for DIAMOND hits.

        Args:
            collection (str): identifier of collection
        """
        try:
            result = float(self.parser[collection]['identity_cutoff'])
        except KeyError:
            result = float(self.parser['DEFAULT']['identity_cutoff'])
        return result

    def get_evalue_cutoff(self, collection):
        """Returns e-value threshold for DIAMOND search.

        Args:
            collection (str): identifier of collection
        """
        try:
            result = float(self.parser[collection]['evalue_cutoff'])
        except KeyError:
            result = float(self.parser['DEFAULT']['evalue_cutoff'])
        return result

    def get_length_cutoff(self, collection):
        """Returns minimal length of amino acid alignment of DIAMOND hits.

        Args:
            collection (str): identifier of collection
        """
        try:
            result = int(self.parser[collection]['length_cutoff'])
        except KeyError:
            result = int(self.parser['DEFAULT']['length_cutoff'])
        return result

    def get_overlap_cutoff(self, collection):
        """Returns maximal length of overlapping part for two alignments
        to be considered the same hit.

        Args:
            collection (str): identifier of collection
        """
        try:
            result = int(self.parser[collection]['hits_overlap_cutoff'])
        except KeyError:
            result = int(self.parser['DEFAULT']['hits_overlap_cutoff'])
        return result

    def get_biscore_range_cutoff(self, collection):
        """Returns lowest acceptaple relative bitscore (relative to top hit bit-score)

        Args:
            collection (str): identifier of collection
        """
        try:
            result = float(self.parser[collection]['biscore_range_cutoff'])
        except KeyError:
            result = float(self.parser['DEFAULT']['biscore_range_cutoff'])
        return result

    def get_ranks_cutoffs(self, collection):
        """Returns dictionary of amino acid identity % values for each taxonomy level

        Args:
            collection (str): identifier of collection
        """
        ret_val = {}
        if self.parser.has_option(collection, 'rank_cutoffs'):
            cutoffs = self.parser[collection]['rank_cutoffs']
            vals = cutoffs.split(',')
            for k_v in zip(vals[::2], vals[1::2]):
                ret_val[k_v[0]] = float(k_v[1])
            ret_val['norank'] = self.get_identity_cutoff(collection)
        else:
            pass
        return ret_val

    def get_reference_diamond_db(self, collection):
        """Returns path to DIAMOND database for reference search

        Args:
            collection (str): identifier of collection
        """
        try:
            result = self.parser[collection]['reference_diamond_db']
        except KeyError:
            result = self.parser['DEFAULT']['reference_diamond_db']
        return result

    def get_background_diamond_db(self, collection):
        """Returns path to DIAMOND database for background search

        Args:
            collection (str): identifier of collection
        """
        try:
            result = self.parser[collection]['background_diamond_db']
        except KeyError:
            result = self.parser['DEFAULT']['background_diamond_db']
        return result

    def get_reference_db_size(self, collection):
        """Returns size of DIAMOND database for reference search

        Args:
            collection (str): identifier of collection
        """
        try:
            result = int(self.parser[collection]['reference_db_size'])
        except KeyError:
            result = int(self.parser['DEFAULT']['reference_db_size'])
        return result

    def get_background_db_size(self, collection):
        """Returns size of DIAMOND database for background search

        Args:
            collection (str): identifier of collection
        """
        try:
            result = int(self.parser[collection]['background_db_size'])
        except KeyError:
            result = int(self.parser['DEFAULT']['background_db_size'])
        return result

    #~ @property
    #~ def taxonomy_names_file(self):
        #~ """Path to NCBI taxonomy names.dmp file"""
        #~ return self.parser['DEFAULT']['taxonomy_names_file']

    #~ @property
    #~ def taxonomy_nodes_file(self):
        #~ """Path to NCBI taxonomy nodes.dmp file"""
        #~ return self.parser['DEFAULT']['taxonomy_nodes_file']

    #~ @property
    #~ def taxonomy_merged_file(self):
        #~ """Path to NCBI taxonomy merged.dmp file"""
        #~ return self.parser['DEFAULT']['taxonomy_merged_file']

    @property
    def collections(self):
        """List of collection names"""
        return self.parser.sections()

    @property
    def diamond_path(self):
        """Path to DIAMOND binary"""
        return self.parser['DEFAULT']['aligner_path']

    @property
    def microbecensus_datadir(self):
        """Path to MicrobeCensus data directory"""
        return self.parser['DEFAULT']['microbecensus_data']

    @property
    def krona_path(self):
        """Path to KronaTools"""
        return self.parser['DEFAULT']['krona_path']

    @property
    def megahit_path(self):
        """Path to MEGAHIT binary"""
        return self.parser['DEFAULT']['megahit_path']

    @property
    def metaspades_path(self):
        """Path to metaSPAdes binary"""
        return self.parser['DEFAULT']['metaspades_path']

    @property
    def bowtie_indexer_path(self):
        """Path to Bowtie indexer binary"""
        return self.parser['DEFAULT']['bowtie_indexer_path']

    @property
    def bowtie_path(self):
        """Path to Bowtie aligner binary"""
        return self.parser['DEFAULT']['bowtie_path']

    @property
    def prodigal_path(self):
        """Path to Prodigal binary"""
        return self.parser['DEFAULT']['prodigal_path']
