"""Describes ReferenceData class"""
from collections import defaultdict
from lib.utils.utils import singleton
from lib.utils.const import RANKS


@singleton
class ReferenceData(object):
    """ReferenceData object stores functional reference data for a collection.
    This is a singleton class.

    Attributes:
        functions_dict (:obj:'defaultdict'[str,dict[str,str]]): dictionary
            of functions, outer key is function identifier, inner keys are 'name' and 'group name'
        proteins_dict (:obj:'defaultdict'[str,dict[str,str]]): dictionary
            of reference proteins, outer key is protein identifier, inner
            keys are 'taxid' (for taxonomy ID), 'function' (for concatenated
            function IDs) and 'source' (source DB)
    """

    def __init__(self, config, collection):
        """Args:
           config (:obj:'ProgramConfig'): Fama configuration parameters
           collection (str): collection identifier
        """
        self.default_identity_threshold = config.get_identity_cutoff(collection)
        self.default_ranks_thresholds = config.get_ranks_cutoffs(collection)
        self.functions_dict = defaultdict(dict)
        self.proteins_dict = defaultdict(dict)
        self.load_reference_data(config, collection)

    def load_reference_data(self, config, collection):
        """Loads reference data for a given collection"""
        self.initialize_functions_dict(
            config.get_functions_file(collection)
            )
        self.initialize_proteins_dict(
            config.get_protein_list_file(collection)
            )

    def initialize_functions_dict(self, infile):
        """ Reads reference data and populates functions_dict"""
        print('Loading ', infile)
        with open(infile, 'r') as file_handle:
            for line in file_handle:
                if line.startswith('#'):
                    continue  # skip header and comments
                else:
                    line = line.rstrip('\n\r')
                    row = line.split('\t')
                    function = row[0]
                    if len(row) > 2:
                        self.functions_dict[function]['name'] = row[1]
                        self.functions_dict[function]['group'] = row[2]
                    if len(row) > 3:
                        self.functions_dict[function]['function_cutoff'] = float(row[3])
                        self.functions_dict[function]['norank'] = float(row[3])
                        n = 4
                        for rank in reversed(RANKS[1:]):
                            self.functions_dict[function][rank] = float(row[n])
                            n += 1
        print(len(self.functions_dict), ' functions found')

    def initialize_proteins_dict(self, infile):
        """ Reads reference data and populates proteins_dict"""
        print('Loading ', infile)
        with open(infile, 'r') as file_handle:
            for line in file_handle:
                if line.startswith('#'):
                    continue  # skip header and comments
                else:
                    line = line.rstrip('\n\r')
                    line_tokens = line.split('\t')
                    if len(line_tokens) > 3:
                        self.proteins_dict[line_tokens[0]]['taxid'] = line_tokens[1]
                        self.proteins_dict[line_tokens[0]]['function'] = line_tokens[-1]
                        self.proteins_dict[line_tokens[0]]['source'] = line_tokens[-2]
        print(len(self.proteins_dict), ' reference proteins found')

    def lookup_protein_function(self, protein):
        """Returns list of function identifiers assigned to a reference protein"""
        ret_val = []
        try:
            ret_val = self.proteins_dict[protein]['function'].split('|')
        except KeyError:
            # print('Protein', protein, 'not found in reference database. Unknown function reported')
            ret_val.append('')
        return ret_val

    def lookup_protein_tax(self, protein):
        """Returns taxonomy identifier assigned to a reference protein"""
        ret_val = ''
#        print('Lookup tax for ', protein)
        try:
            ret_val = self.proteins_dict[protein]['taxid']
        except KeyError:
            print('Protein not found in reference database', protein, 'Taxonomy ID set to zero')
            ret_val = '0'
        return ret_val

    def lookup_function_name(self, func_id):
        """Returns function name for a given function identifier"""
        ret_val = ''
        if func_id in self.functions_dict:
            ret_val = self.functions_dict[func_id]['name']
        return ret_val

    def lookup_function_group(self, func_id):
        """Returns function group name for a given function identifier"""
        ret_val = ''
        if func_id in self.functions_dict:
            ret_val = self.functions_dict[func_id]['group']
        return ret_val

    def get_functions_in_group(self, group):
        """Returns list of function identifiers having the same function group name"""
        ret_val = []
        try:
            ret_val = [function for function in self.functions_dict if
                       self.functions_dict[function]['group'] == group]
        except KeyError:
            print('Group not found:', group)
            raise
        return ret_val

    def list_functions(self):
        """Returns list of all available functions"""
        return sorted(self.functions_dict.keys())

    def lookup_identity_threshold(self, function=None, rank=None):
        """Returns amino acid identity threshold for a given function and taxonomy rank.
        If taxonomy rank is None, returns amino acid identity threshold for function prediction.
        If function is None and rank is not None, returns amino acid identity threshold
        for the taxonomy rank defined in config.ini.

        If both function and taxonomy rank are None, returns default threshold defined in config.ini

        Args:
            function (str): function identifier
            rank (str): taxonomy rank defined in lib.utils.const.RANKS

        Returns:
            result (float): amino acid identity cutoff

        """
        result = self.default_identity_threshold
        if function == '':
            return result
        try:
            if function is not None and rank is not None:
                result = self.functions_dict[function][rank]
            elif function is not None:
                result = self.functions_dict[function]['function_cutoff']
            elif rank is not None:
                result = self.default_ranks_thresholds[rank]
        except KeyError:
            pass
        return result
