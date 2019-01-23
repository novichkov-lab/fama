from collections import defaultdict
from Fama.ProjectUtil.ProgramConfig import ProgramConfig

class ReferenceData:

    def __new__(cls, options):
        if not hasattr(cls, 'instance'):
            cls.instance = super(ReferenceData, cls).__new__(cls)
        cls.instance.options = options
        return cls.instance

    def __init__(self,options):
        self.functions_dict = defaultdict(dict)
        self.proteins_dict = defaultdict(dict)

    def load_reference_data(self,collection):
        self.functions_dict = self.initialize_functions_dict(self.options.get_functions_file(collection))
        self.proteins_dict = self.initialize_proteins_dict(self.options.get_proteins_list_file(collection))

    def initialize_functions_dict(self,infile):
        _dict = defaultdict(dict)
        print ('Loading ', infile)
        with open(infile, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue #skip header and comments
                else:
                    line = line.rstrip('\n\r')
                    line_tokens = line.split('\t')
                    if len(line_tokens) > 2:
                        _dict[line_tokens[0]]['name'] = line_tokens[1]
                        _dict[line_tokens[0]]['group'] = line_tokens[2]
        print (len(_dict), ' functions found')
        return _dict
        
    def initialize_proteins_dict(self,infile):
        _dict = defaultdict(dict)
        print ('Loading ', infile)
        with open(infile, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue #skip header and comments
                else:
                    line = line.rstrip('\n\r')
                    line_tokens = line.split('\t')
                    if len(line_tokens) > 3:
                        _dict[line_tokens[0]]['taxid'] = line_tokens[1]
                        _dict[line_tokens[0]]['function'] = line_tokens[-1]
                        _dict[line_tokens[0]]['source'] = line_tokens[-2]
        print (len(_dict), ' reference proteins found')
        return _dict
        
    def lookup_protein_function(self, protein):
        ret_val = []
        if protein in self.proteins_dict:
            ret_val = self.proteins_dict[protein]['function'].split('|')
            return ret_val        
        ret_val.append('')
        return ret_val

    def lookup_protein_tax(self, protein):
        ret_val = ''
#        print('Lookup tax for ', protein)
        if protein in self.proteins_dict:
#            print ('found ', self.proteins_dict[protein]['taxid'])
            ret_val = self.proteins_dict[protein]['taxid']
        return ret_val
            
    def lookup_function_name(self, func_id):
        ret_val = ''
        if func_id in self.functions_dict:
            ret_val = self.functions_dict[func_id]['name']
        return ret_val

    def lookup_function_group(self, func_id):
        ret_val = ''
        if func_id in self.functions_dict:
            ret_val = self.functions_dict[func_id]['group']
        return ret_val

    def get_functions_in_group(self, group):
        return [function for function in self.functions_dict if self.functions_dict[function]['group'] == group]
