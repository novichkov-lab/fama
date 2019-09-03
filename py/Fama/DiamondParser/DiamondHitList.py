from Fama.DiamondParser.DiamondHit import DiamondHit

class DiamondHitList(object):
    """DiamondHitList stores a set of DiamondHit objects for one
    query sequence (usually, query is a sequence read or a protein)
    
    """
    def __init__(self, query_id = None):
        """ Args:
            query_id (str): query sequence identifier
        
        """
        self._query_id = query_id
        self._data = []

    def add_hit(self, hit):
        """Adds a DiamondHit to the DiamondHitList. 
        
        Note: checks if query_id of DiamondHit is identical to query_id of
        the DiamondHitList. If they don't, new DiamondHit will not be added.
        Then, checks if a DiamondHit with the same subject_id, q_start, q_end 
        already exists in the DiamondHitList. If it does, new DiamondHit
        will not be added.
        
        Args:
            hit (:obj:'DiamondHit'): DiamondHit to be added
            
        """
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
        """:obj:'list' of :obj:'DiamondHit': list of DIAMOND hits"""
        return self._data

    @hits.setter
    def hits(self, var):
        self._data = var

    @property
    def query_id(self):
        """str: query sequence identifier"""
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
        """:obj:'int': number of DiamondHit objects in the list"""
        return len(self.hits)
    
    def filter_list(self, overlap_cutoff):
        """Filters list of DiamondHit objects in the DiamondHitList. 
        Removes hits, which overlap by more than 'overlap_cutoff' base 
        pairs with any hit with higher bit-score. 
        
        Args:
            overlap_cutoff(int): minimal length of a common area between hits to be considered overlapping
        """
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
        """Returns False if query sequence of new_hit overlap 
        by at least 'overlap_cutoff' base pairs with query sequence of
        at least one DiamondHit from hit_list list of DiamondHit objects. 
        Otherwise, returns True
        
        Args:
            new_hit(:obj:'DiamondHit'): a hit to be tested
            hit_list(:obj:'list' of :obj:'DiamondHit'): list of hits to be tested on
            overlap_cutoff(int): minimal length of a common area between hits to be considered overlapping
            
        Returns:
            bool

        """
        for hit in hit_list:
            if self.hits_do_overlap(new_hit, hit, overlap_cutoff):
                return False
        return True
    
    def hits_do_overlap (self, new_hit, hit, overlap_cutoff):
        """Returns True if query sequence of new_hit DiamondHit overlap 
        by at least 'overlap_cutoff' base pairs with query sequence of 
        hit DiamondHit. Otherwise, returns False.
        
        Args:
            new_hit(:obj:'DiamondHit'): a hit to be tested
            hit(:obj:'DiamondHit'): a hit to be tested on
            overlap_cutoff(int): minimal length of a common area between hits to be considered overlapping
            
        Returns:
            bool

        """
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
        """Returns True if bitscore of new_hit is higher than bitscore of
        all overlapping DiamondHit objects from hit_list list. 
        Otherwise, returns False
        
        Args:
            new_hit(:obj:'DiamondHit'): a hit to be tested
            hit_list(:obj:'list' of :obj:'DiamondHit'): list of hits to be tested on
            overlap_cutoff(int): minimal length of a common area between hits to be considered overlapping
            
        Returns:
            bool
        
        """
        for hit in hit_list:
            if self.hits_do_overlap(new_hit, hit, overlap_cutoff):
                if new_hit.bitscore <= hit.bitscore:
                    return False
        return True
            
    def replace_hit(self, new_hit, hit_list, overlap_cutoff):
        """Compares all DiamondHit objects from a list with a 
        new_hit DiamondHit object. Removes all hits, which overlap by more 
        than 'overlap_cutoff' base pairs with new_hit and have lower bit-score. 
        If any hits were removed, new_hit is added to the list.
        
        Args:
            new_hit(:obj:'DiamondHit'): a hit to be tested
            hit_list(:obj:'list' of :obj:'DiamondHit'): list of hits to be tested on
            overlap_cutoff(int): minimal length of a common area between hits to be considered overlapping
        
        Returns:
            :obj:'list' of :obj:'DiamondHit': updated list of DIAMOND hits
        
        """

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
        """Removes a DiamondHit from the DiamondHitList.
        
        Finds identical DiamondHit object in the DiamondHitList and deletes it.
        
        Args:
            hit_to_remove(:obj:'DiamondHit'): a hit to be removed
        
        """
        index_remove = None
        for index,hit in (enumerate(self.hits)):
            if hit == hit_to_remove:
                index_remove = index
                break
        if index_remove != None:
            del self.hits[index]
    
    def annotate_hits(self, reference_data):
        """Assigns function to each DiamondHit in the DiamondHitList"""
        for hit in self.hits:
            hit.annotate_hit(reference_data)
    
    def print_hits(self):
        """Prints string representation of each DiamondHit in the DiamondHitList"""
        for hit in self.hits:
            print(hit)
    
    def __str__(self):
        """Returns string representation of DiamondHitList object"""
        return ('Hit List '+'\n'.join(str(hit) for hit in self.hits))
