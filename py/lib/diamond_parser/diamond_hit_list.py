"""Describes DiamondHitList class"""
import lib.diamond_parser.hit_utils as hit_utils


class DiamondHitList(object):
    """DiamondHitList stores a set of DiamondHit objects for one
    query sequence (usually, query is a sequence read or a protein)

    """
    def __init__(self, query_id=None):
        """ Args:
            query_id (str): query sequence identifier

        """
        self._query_id = query_id
        self.data = []

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
            for existing_hit in self.hits:
                if existing_hit.subject_id == hit.subject_id and (
                        existing_hit.q_start == hit.q_start
                ) and (
                    existing_hit.q_end == hit.q_end
                ):
                    hit_exists = True
                    break
            if not hit_exists:
                self.hits.append(hit)
        else:
            print('Diamond hit was not added to the list: different query IDs:'
                  + self.query_id + ',' + hit.query_id)

    @property
    def hits(self):
        """:obj:'list' of :obj:'DiamondHit': list of DIAMOND hits"""
        return self.data

    @hits.setter
    def hits(self, var):
        self.data = var

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
            overlap_cutoff(int): minimal length of a common area between
                hits to be considered overlapping
        """
        temp_list = []

        for hit in self.hits:
            if not temp_list:
                temp_list.append(hit)
            elif not hit_utils.hit_overlaps_any_hits(hit, temp_list, overlap_cutoff):
                temp_list.append(hit)
            elif hit_utils.has_higher_score(hit, temp_list, overlap_cutoff):
                temp_list = hit_utils.replace_hit(hit, temp_list, overlap_cutoff)
        self.hits = temp_list


    def filter_list_by_identity(self, reference_data):
        """Filters list of DiamondHit objects in the DiamondHitList.
        Removes hits, which have amino acid identity below threshold defined
        in reference data or config file.
        WARNING: before filtering, functions must be assigned to assigned 
        to all hits (i.e. call annotate_hits function first)

        Args:
            reference_data (:obj:ReferenceData): functional reference data
        """
        temp_list = []

        for hit in self.hits:
            max_threshold = 0.0
            for function in hit.functions:
                function_threshold = reference_data.lookup_identity_threshold(function=function)
                if max_threshold < function_threshold:
                    max_threshold = function_threshold
            if max_threshold == 0.0:
                max_threshold = reference_data.lookup_identity_threshold()
            if hit.identity >= max_threshold:
                temp_list.append(hit)
        self.hits = temp_list


    def remove_hit(self, hit_to_remove):
        """Removes a DiamondHit from the DiamondHitList.

        Finds identical DiamondHit object in the DiamondHitList and deletes it.

        Args:
            hit_to_remove(:obj:'DiamondHit'): a hit to be removed

        """
        index_remove = None
        for hit_index, hit in enumerate(self.hits):
            if hit == hit_to_remove:
                index_remove = hit_index
                break
        if index_remove is not None:
            del self.hits[index_remove]

    def remove_hit_by_index(self, index_remove):
        """Removes a DiamondHit from the DiamondHitList by its index.

        Finds DiamondHit object in the DiamondHitList and deletes it.

        Args:
            index_remove(int): index of the hit to be removed

        """
        try:
            del self.hits[index_remove]
        except IndexError:
            print('Unable to remove hit with index', index_remove)

    def remove_all_hits(self):
        """Assigns function to each DiamondHit in the DiamondHitList"""
        self.hits = []

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
        return 'Hit List '+'\n'.join(str(hit) for hit in self.hits)
