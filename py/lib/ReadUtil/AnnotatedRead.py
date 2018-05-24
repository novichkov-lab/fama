from lib.DiamondParser.DiamondHit import DiamondHit
from lib.DiamondParser.DiamondHitList import DiamondHitList

class AnnotatedRead:
    def __init__(self, read_id):
        self.read_id = read_id
        self.sequence = None
        self.quality = None
        self.line3 = None
        self.hits = None
        self.status = 'Unaccounted'
        self.functions = {}
        
    def get_read_id(self):
        return self.read_id
        
    def set_hits(self,hits):
        self.hits = hits
        
    def get_hits(self):
        return self.hits
        
    def set_sequence(self, seq):
        self.sequence = seq
        
    def set_quality(self, quality):
        self.quality = quality
    
    def set_line3(self,line3):
        self.line3 = line3
        
    def get_sequence(self):
        return self.sequence
        
    def get_quality(self):
        return self.quality
        
    def get_line3(self):
        return self.line3
    
    def show_hits(self):
        self.hits.print_hits()
    
    def compare_hits(self, hit_start, hit_end, new_hit_list, biscore_range_cutoff):
        for hit in self.hits.get_hits():
            new_
            if hit.get_query_start() = hit_start and hit.get_query_end() = hit_end:
                bitscore = hit.get_bitscore()
                bitscore_lower_cutoff = bitscore*(1-biscore_range_cutoff)
                bitscore_upper_cutoff = bitscore*(1+biscore_range_cutoff)
                
                # first, check if best new hit has the same subject as the existing hit
                if new_hit_list.get_hits()[0].get_subject_id() == hit.get_subject_id():
                    #how far we are from next best hit?
                    if new_hit_list.get_hits()[0].get_bitscore() > bitscore_lower_cutoff:
                        # top hits are too close, select the most frequent function
                        self.status = 'True best hit, close to others'
                        new_functions_dict = {}
                        for new_hit in new_hit_list.get_hits():
                            if new_hit.get_bitscore() > bitscore_lower_cutoff:
                                for new_fuction in new_hit.get_functions:
                                    if new_fuction in new_functions_dict:
                                        new_functions_dict[new_fuction] += 1
                                    else:
                                        new_functions_dict[new_fuction] = 1
                        best_new_function_count = 0
                        best_new_function = None
                        for new_fuction in new_functions_dict:
                            if new_functions_dict[new_fuction] > best_new_function_count:
                                best_new_function_count = new_functions_dict[new_fuction]
                                best_new_function = new_fuction
                        if best_new_function: 
                            self.functions[best_new_function] = 1
                    else:
                        self.status = 'True best hit, far better than others'
                        for function in hit.get_functions():
                            self.functions[function] = 1
                    
                # second, check if best new hit is better than existing hit
                elif new_hit_list.get_hits()[0].get_bitscore() > bitscore_upper_cutoff:
                    new_bitscore_lower_cutoff = new_hit_list.get_hits()[0].get_bitscore()*(1-biscore_range_cutoff)
                    #compare functions of old and new best hits
                    if new_hit_list.get_hits()[0].get_functions()[0] == '':
                        # new hit was not in reference DB, select the most frequent function
                        new_functions_dict = {}
                        for new_hit in new_hit_list.get_hits():
                            if new_hit.get_bitscore() > new_bitscore_lower_cutoff:
                                for new_fuction in new_hit.get_functions:
                                    if new_fuction in new_functions_dict:
                                        new_functions_dict[new_fuction] += 1
                                    else:
                                        new_functions_dict[new_fuction] = 1
                        best_new_function_count = 0
                        best_new_function = None
                        for new_fuction in new_functions_dict:
                            if new_functions_dict[new_fuction] > best_new_function_count:
                                best_new_function_count = new_functions_dict[new_fuction]
                                best_new_function = new_fuction
                        if best_new_function: 
                            self.functions[best_new_function] = 1
                    elif len(hit.get_functions()) == len(new_hit_list.get_hits()[0].get_functions()):
                        # new hit has the same number of functions as old hit. Compare it one-by-one
                        flag = True
                        for index,function in enumerate(sorted(hit.get_functions())):
                            if function != sorted(new_hit_list.get_hits()[0].get_functions())[index]:
                                #different functions, stop here
                                self.status = 'Best hit not in reference DB'
                                flag = False
                                break
                        if flag:
                            # different best hits, but identical functions
                            self.status = 'Function from reference DB, different best hits'
                        else:
                            # different best hits, different functions
                            self.status = 'Function from reference DB, differs from function identified by reference DB search'
                    else:
                        # different best hits, different functions
                        self.status = 'Function from reference DB, differs from function identified by reference DB search'
                else:
                    #what if best score in the second search is not so high? 
                    hit_functions = 
                    
                        
                    
                
        
    def __str__(self):
        return str(self.hits) + '\t' + self.function_id
