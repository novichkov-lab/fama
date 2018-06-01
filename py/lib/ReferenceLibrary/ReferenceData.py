from lib.ProjectUtil.ProgramConfig import ProgramConfig

class ReferenceData:

    def __new__(cls, options):
        if not hasattr(cls, 'instance'):
            cls.instance = super(ReferenceData, cls).__new__(cls)
        cls.instance.options = options
        return cls.instance

    def __init__(self,options):
        self.functions_dict = {}
        self.proteins_dict = {}

    def load_reference_data(self,collection):
        self.functions_dict = self.initialize_functions_dict(self.options.get_functions_file(collection))
        self.proteins_dict = self.initialize_proteins_dict(self.options.get_proteins_list_file(collection))

    def initialize_functions_dict(self,infile):
        _dict = {}
        print ('Loading ', infile)
        with open(infile, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue #skip header and comments
                else:
                    line = line.rstrip('\n\r')
                    line_tokens = line.split('\t')
                    if len(line_tokens) > 1:
                        _dict[line_tokens[0]] = line_tokens[1]
        print (len(_dict), ' functions found')
        return _dict
        
    def initialize_proteins_dict(self,infile):
        _dict = {}
        print ('Loading ', infile)
        with open(infile, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue #skip header and comments
                else:
                    line = line.rstrip('\n\r')
                    line_tokens = line.split('\t')
                    if len(line_tokens) > 3:
                        _dict[line_tokens[0]] = line_tokens[3]
        print (len(_dict), ' reference proteins found')
        return _dict
        
    def lookup_protein_function(self, protein):
        ret_val = []
        if protein in self.proteins_dict:
            ret_val = self.proteins_dict[protein].split('|')
            return ret_val
        else: 
            #if protein ID contains starts with function ID
            protein_tokens = protein.split('_')
            protein = protein_tokens[-1]
            if protein in self.proteins_dict:
                ret_val = self.proteins_dict[protein].split('|')
                return ret_val
        
        ret_val.append('')
        return ret_val
            
    def lookup_function_name(self, func_id):
        ret_val = ''
        if func_id in self.functions_dict:
            ret_val = self.functions_dict[func_id]
        return ret_val
