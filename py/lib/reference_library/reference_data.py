from collections import defaultdict
from lib.project.program_config import ProgramConfig

class ReferenceData(object):
    """ReferenceData object stores functional reference data for a collection.
    This is a singleton.
    
    Attributes:
        functions_dict (:obj:'defaultdict'[str,dict[str,str]]): dictionary 
            of functions, outer key is function identifier, inner keys are 'name' and 'group name'
        proteins_dict (:obj:'defaultdict'[str,dict[str,str]]): dictionary 
            of reference proteins, outer key is protein identifier, inner
            keys are 'taxid' (for taxonomy ID), 'function' (for concatenated 
            function IDs) and 'source' (source DB)
    """

    def __new__(cls, config):
        """Args:
           config (:obj:'ProgramConfig') Fama configuration parameters
        """
        if not hasattr(cls, 'instance'):
            cls.instance = super(ReferenceData, cls).__new__(cls)
        cls.instance.config = config
        return cls.instance

    def __init__(self,config):
        """Args:
           config (:obj:'ProgramConfig') Fama configuration parameters
        """
        self.functions_dict = defaultdict(dict)
        self.proteins_dict = defaultdict(dict)

    def load_reference_data(self,collection):
        """Loads reference data for a given collection"""
        self.functions_dict = self.initialize_functions_dict(self.config.get_functions_file(collection))
        self.proteins_dict = self.initialize_proteins_dict(self.config.get_protein_list_file(collection))

    def initialize_functions_dict(self,infile):
        """ Reads reference data and populates functions_dict"""
        result_dict = defaultdict(dict)
        print ('Loading ', infile)
        with open(infile, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue #skip header and comments
                else:
                    line = line.rstrip('\n\r')
                    line_tokens = line.split('\t')
                    if len(line_tokens) > 2:
                        result_dict[line_tokens[0]]['name'] = line_tokens[1]
                        result_dict[line_tokens[0]]['group'] = line_tokens[2]
        print (len(result_dict), ' functions found')
        return result_dict
        
    def initialize_proteins_dict(self,infile):
        """ Reads reference data and populates proteins_dict"""
        result_dict = defaultdict(dict)
        print ('Loading ', infile)
        with open(infile, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue #skip header and comments
                else:
                    line = line.rstrip('\n\r')
                    line_tokens = line.split('\t')
                    if len(line_tokens) > 3:
                        result_dict[line_tokens[0]]['taxid'] = line_tokens[1]
                        result_dict[line_tokens[0]]['function'] = line_tokens[-1]
                        result_dict[line_tokens[0]]['source'] = line_tokens[-2]
        print (len(result_dict), ' reference proteins found')
        return result_dict
        
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
        return [function for function in self.functions_dict if self.functions_dict[function]['group'] == group]
        
    def list_functions(self):
        """Returns list of all available functions"""
        return sorted(self.functions_dict.keys())
        
