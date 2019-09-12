"""Describes ReferenceData class"""
from collections import defaultdict

def singleton(cls):
    """Implements singleton design pattern"""
    instances = {}
    def getinstance(*args, **kwargs):
        """Creates singleton instance of cls, if not exists, and returns it"""
        if cls not in instances:
            instances[cls] = cls(*args, **kwargs)
        return instances[cls]
    return getinstance

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
           config (:obj:'ProgramConfig') Fama configuration parameters
        """
        self.functions_dict = defaultdict(dict)
        self.proteins_dict = defaultdict(dict)
        self.load_reference_data (config, collection)

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
                    continue #skip header and comments
                else:
                    line = line.rstrip('\n\r')
                    line_tokens = line.split('\t')
                    if len(line_tokens) > 2:
                        self.functions_dict[line_tokens[0]]['name'] = line_tokens[1]
                        self.functions_dict[line_tokens[0]]['group'] = line_tokens[2]
        print(len(self.functions_dict), ' functions found')

    def initialize_proteins_dict(self, infile):
        """ Reads reference data and populates proteins_dict"""
        print('Loading ', infile)
        with open(infile, 'r') as file_handle:
            for line in file_handle:
                if line.startswith('#'):
                    continue #skip header and comments
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
        if protein in self.proteins_dict:
            ret_val = self.proteins_dict[protein]['function'].split('|')
            return ret_val
        ret_val.append('')
        return ret_val

    def lookup_protein_tax(self, protein):
        """Returns taxonomy identifier assigned to a reference protein"""
        ret_val = ''
#        print('Lookup tax for ', protein)
        if protein in self.proteins_dict:
#            print ('found ', self.proteins_dict[protein]['taxid'])
            ret_val = self.proteins_dict[protein]['taxid']
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
        """Returns list of function identifiers haveing the same function group name"""
        return [function for function in self.functions_dict if
                self.functions_dict[function]['group'] == group]

    def list_functions(self):
        """Returns list of all available functions"""
        return sorted(self.functions_dict.keys())
