from Fama.DiamondParser.DiamondHit import DiamondHit

class DiamondHitList:
    def __init__(self, query_id = None):
        self._query_id = query_id
        self._data = []

    def add_hit(self, hit):
        # this function adds a hit to the list. After hit with such coordinates and subject ID already exists in the list, it will not be added. 
        if hit.query_id == self.query_id:
            hit_exists = False
            for h in self.hits:
                if h.subject_id == hit.subject_id and h.q_start == hit.q_start and h.q_end == hit.q_end:
                    hit_exists = True
                    break
            if not hit_exists:
                self.hits.append(hit)
        else:
            print ('Diamond hit was not added to the list: different query IDs:' + self.query_id + ',' + hit.query_id)

    @property
    def hits(self):
        return self._data

    @hits.setter
    def hits(self, var):
        self._data = var

    @property
    def query_id(self):
        return self._query_id
    
    @query_id.setter
    def query_id(self, var):
        if self.hits:
            for hit in self.hits:
                if hit.query_id != var:
                    raise ValueError
        self._query_id = var
            
    @property
    def hits_number(self):
        return len(self.hits)
    
    def filter_list(self, overlap_cutoff):
        temp_list = []
        
        for hit in self.hits:
            if len(temp_list) == 0:
                temp_list.append(hit)
            elif self.is_not_overlapping(hit, temp_list, overlap_cutoff):
                temp_list.append(hit)
            elif self.has_best_score(hit, temp_list, overlap_cutoff):
                temp_list = self.replace_hit(hit, temp_list, overlap_cutoff)
        self.hits = temp_list
        
    def is_not_overlapping(self, new_hit, hit_list, overlap_cutoff):
        # returns True if positions in query sequence of new_hit  overlap 
        # by at least overlap_cutoff with any hit from hit_list, False if not
        for hit in hit_list:
            if self.hits_do_overlap(new_hit, hit, overlap_cutoff):
                return False
        return True
    
    def hits_do_overlap (self, new_hit, hit, overlap_cutoff):
        # returns True if positions in query sequence of two hits overlap 
        # by at least overlap_cutoff, False if not
        new_hit_start = new_hit.q_start
        new_hit_end = new_hit.q_end
        start = hit.q_start
        end = hit.q_end
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
        # returns True if new_hit's bitscore is higher than bitscores of all overlapping hits in hit_list
        for hit in hit_list:
            if self.hits_do_overlap(new_hit, hit, overlap_cutoff):
                if new_hit.bitscore <= hit.bitscore:
                    return False
        return True
            
    def replace_hit(self, new_hit, hit_list, overlap_cutoff):
        # this function check if _hit_list_ contains Diamond hits overlapping with _new_hit_ that have lower bitscore than _new_hit_ 
        # if it finds such hits, they are removed from list and _new_hit_ is appended to the end of the list
        ret_val = hit_list
        append_flag = False
        for index,hit in reversed(list(enumerate(hit_list))):
            if self.are_two_hits_overlap(new_hit, hit, overlap_cutoff):
                if new_hit.bitscore > hit.bitscore:
                    del ret_val[index]
                    append_flag = True
        if append_flag:
            ret_val.append(new_hit)
        return ret_val
    
    def remove_hit(self, hit_to_remove):
        index_remove = None
        for index,hit in (enumerate(self.hits)):
            if hit == hit_to_remove:
                index_remove = index
                break
        if index_remove != None:
            del self.hits[index]
    
    
    def annotate_hits(self, reference_data):
        for hit in self.hits:
            hit.annotate_hit(reference_data)
    
    def print_hits(self):
        for hit in self.hits:
            print(hit)
    
    def __str__(self):
        return ('Hit List '+'\n'.join(str(hit) for hit in self.hits))
