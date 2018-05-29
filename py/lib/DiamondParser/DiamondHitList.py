from lib.DiamondParser.DiamondHit import DiamondHit

class DiamondHitList:
    def __init__(self, query_id):
        self.query_id = query_id
        self.data = []

    def add_hit(self, hit):
        # this function adds a hit into the list. After adding a hit, the list 
        # may need filtering, since this function performs no checks for overlap. 
        if hit.get_query_id() == self.query_id:
            self.data.append(hit)
        else:
            print ('Diamond hit was not added into the list: different query IDs:' + self.query_id + ','+hit.get_query_id)

    def get_hits(self):
        return self.data

    def set_query_id(self, query_id):
        if self.data:
            for hit in self.data:
                if hit.get_query_id() != query_id:
                    raise ValueError
        self.query_id = query_id
            
    def get_hits_number(self):
        return len(self.data)
    
    def filter_list(self, overlap_cutoff):
        temp_list = []
        
        for hit in self.data:
            if len(temp_list) == 0:
                temp_list.append(hit)
            elif self.is_not_overlapping(hit, temp_list, overlap_cutoff):
                temp_list.append(hit)
            elif self.has_best_score(hit, temp_list, overlap_cutoff):
                temp_list = self.replace_hit(hit, temp_list, overlap_cutoff)
        self.data = temp_list
        
    def is_not_overlapping(self, new_hit, hit_list, overlap_cutoff):
        for hit in hit_list:
            if self.are_two_hits_overlap(new_hit, hit, overlap_cutoff):
                return False
        return True
    
    def are_two_hits_overlap (self, new_hit, hit, overlap_cutoff):
        new_hit_start = new_hit.get_query_start()
        new_hit_end = new_hit.get_query_end()
        start = hit.get_query_start()
        end = hit.get_query_end()
        # check if hits are on different strands
        if new_hit_start < new_hit_end:
            #existing hit on + strand
            if start > end:
                    # hits on different strands: no overlap
                    return False
        elif new_hit_start > new_hit_end:
            #existing hit on - strand
            if start < end:
                    # hits on different strands: no overlap
                    return False
        # check if hits on one strand do not overlap or have overlap below cutoff
        if new_hit_start < new_hit_end:
            # hits on + strand
            if end < (new_hit_start + overlap_cutoff):
                return False
            elif start > (new_hit_end - overlap_cutoff):
                return False
            else: 
                return True #overlap
        if new_hit_start > new_hit_end:
            # hits on - strand
            if start < (new_hit_end + overlap_cutoff):
                return False
            elif end > (new_hit_start - overlap_cutoff):
                return False
            else: 
                return True #overlap
        return False
        
    def has_best_score(self, new_hit, hit_list, overlap_cutoff):
        new_bitscore = new_hit.get_bitscore()
        for hit in hit_list:
            if self.are_two_hits_overlap(new_hit, hit, overlap_cutoff):
                if new_bitscore <= hit.get_bitscore():
                    return False
        return True
            
    def replace_hit(self, new_hit, hit_list, overlap_cutoff):
        # this functions check if _hit_list_ contains Diamond hits overlapping with _new_hit_ that have lower bitscore than _new_hit_ 
        # if it finds such hits, they are removed from list and _new_hit_ is appended to the end of the list
        new_bitscore = new_hit.get_bitscore()
        ret_val = hit_list
        append_flag = False
        for index,hit in reversed(list(enumerate(hit_list))):
            if self.are_two_hits_overlap(new_hit, hit, overlap_cutoff):
                if new_bitscore > hit.get_bitscore():
                    del ret_val[index]
                    append_flag = True
        if append_flag:
            ret_val.append(new_hit)
        return ret_val
    
    def remove_hit(self, hit_to_remove):
        index_remove = None
        for index,hit in (enumerate(self.data)):
            if hit == hit_to_remove:
                index_remove = index
                break
        if index_remove != None:
            del self.data[index]
    
    
    def annotate_hits(self, reference_data):
        for hit in self.data:
            hit.annotate_hit(reference_data)
    
    def print_hits(self):
        for hit in self.data:
            print(hit)
    
    def __str__(self):
        return ('Hit List '+'/n'.join(str(hit) for hit in self.data))
