"""Describes TaxonomyData class"""
from collections import defaultdict
from lib.utils.const import RANKS, UNKNOWN_TAXONOMY_ID, ROOT_TAXONOMY_ID
from lib.utils.utils import singleton


@singleton
class TaxonomyData(object):
    """TaxonomyData object stores NCBI taxonomic data: taxonomy IDs,
    ranks, taxon names, parent IDs. This is a singleton class.

    Attributes:
        names (defaultdict[str, defaultdict[str,str]]): dictionary
            of taxonomy names, external key is taxonomy identifier,
            internal key is 'name'
        nodes (defaultdict[str, defaultdict[str,str]]): dictionary
            of taxonomy ranks and parent IDs, external key is taxonomy
            identifier, internal keys are 'rank' and 'parent'
    """
    def __init__(self, config):
        """Args:
           config (:obj:'ProgramConfig') Fama configuration parameters
        """
        self.names = defaultdict(dict)
        self.nodes = defaultdict(dict)
        self.load_taxdata(config)

    def load_taxdata(self, config):
        """Loads taxonomic data from NCBI files"""
        names_file = config.taxonomy_names_file
        nodes_file = config.taxonomy_nodes_file
        merged_file = config.taxonomy_merged_file

        # initialize self.names
        print('Loading names file', names_file)
        with open(names_file, 'r') as infile:
            for line in infile:
                line = line.rstrip('\n\r')
                line_tokens = line.split('\t|\t')
                if line_tokens[3] == 'scientific name\t|':
                    self.names[line_tokens[0]]['name'] = line_tokens[1]

        if not self.names:
            raise Exception('Taxonomy names load failed')

        # initialize self.nodes
        print('Loading nodes file', nodes_file)
        with open(nodes_file, 'r') as infile:
            for line in infile:
                line = line.rstrip('\n\r')
                line_tokens = line.split('\t|\t')
                taxid = line_tokens[0]
                parent = line_tokens[1]
                rank = line_tokens[2]
                self.nodes[taxid]['parent'] = parent
                self.nodes[taxid]['rank'] = rank

        # merge
        print('Loading merged file', merged_file)
        with open(merged_file, 'r') as infile:
            for line in infile:
                line = line.rstrip('\n\r')
                line_tokens = line.split('\t')
                old_id = line_tokens[0]
                new_id = line_tokens[2]
                if new_id in self.names:
                    self.names[old_id]['name'] = self.names[new_id]['name']
                    self.nodes[old_id]['parent'] = self.nodes[new_id]['parent']
                    self.nodes[old_id]['rank'] = self.nodes[new_id]['rank']

        # inject 'Unknown' entry
        self.names[UNKNOWN_TAXONOMY_ID]['name'] = 'Unknown'
        self.nodes[UNKNOWN_TAXONOMY_ID]['parent'] = ROOT_TAXONOMY_ID
        self.nodes[UNKNOWN_TAXONOMY_ID]['rank'] = 'norank'
        # make rank of root different from others
        self.nodes[ROOT_TAXONOMY_ID]['rank'] = 'norank'

    def is_exist(self, taxonomy_id):
        """ Checks if taxonomy identifier exists in taxonomy data)

        Args:
            taxonomy_id (str): taxonomy identifier

        Returns:
            result (bool): True if taxonomy_id exists, False otherwise
        """
        result = False
        try:
            if taxonomy_id in self.names:
                result = True
        except KeyError:
            print('Taxonomy identifier %s not found' % taxonomy_id)
        return result

    def get_name(self, taxonomy_id):
        """ Look up taxon name by identifier

        Args:
            taxonomy_id (str): taxonomy identifier

        Returns:
            result (str): taxon name
        """
        result = ''
        try:
            result = self.names[taxonomy_id]['name']
        except KeyError:
            print('Taxonomy identifier %s not found' % taxonomy_id)
        return result

    def get_rank(self, taxonomy_id):
        """ Look up taxon rank by identifier

        Args:
            taxonomy_id (str): taxonomy identifier

        Returns:
            result (str): rank or empty string
        """
        result = 'norank'
        try:
            result = self.nodes[taxonomy_id]['rank']
        except KeyError:
            print('Taxonomy identifier %s not found' % taxonomy_id)
        return result

    def get_parent(self, taxonomy_id):
        """ Look up taxon parent by identifier. Returns direct parent reagrdless
        of its rank. If you need a parent node with rank defined in RANKS, call
        get_upper_level_taxon method instead.

        Args:
            taxonomy_id (str): taxonomy identifier

        Returns:
            result (str): taxonomy dentifier of parent node or empty string
        """
        result = '0'
        try:
            result = self.nodes[taxonomy_id]['parent']
        except KeyError:
            print('Taxonomy identifier %s not found' % taxonomy_id)
        return result

    def get_lca(self, taxonomy_id_list):
        """Returns Lowest Common Ancestor (LCA) taxonomy identifier for a list
        of taxonomy identifiers.
        If any taxon in the list is root, LCA is Unknown.
        If all taxa in the list are Unknowns, LCA is Unknown.
        If any taxon in the list is NOT Unknown, any Unknowns are ignored.

        Args:
            taxonomy_id_list (list of str): taxonomy identifiers

        Returns:
            result (str): LCA taxonomy dentifier
        """
        result = UNKNOWN_TAXONOMY_ID
        if len(taxonomy_id_list) == 1:
            taxonomy_id = taxonomy_id_list.pop()
            if taxonomy_id != '':
                result = taxonomy_id
            return result

        taxonomy_id_list = \
            [taxon_id for taxon_id in taxonomy_id_list if taxon_id != UNKNOWN_TAXONOMY_ID]
        if not taxonomy_id_list:
            # if taxonomy_id_list is empty or contains only Unknowns:
            return result

        taxonomic_lineages = {}
        # Calculate length of the shortest path in taxonomic subtree
        min_depth = 1000
        for taxonomy_id in taxonomy_id_list:
            if not self.is_exist(taxonomy_id):
                continue
            depth = 1
            parent_id = self.get_parent(taxonomy_id)
            lineage = [parent_id, taxonomy_id]
            while self.get_parent(parent_id) != ROOT_TAXONOMY_ID:
                parent_id = self.get_parent(parent_id)
                lineage.insert(0, parent_id)
                depth += 1
            taxonomic_lineages[taxonomy_id] = lineage
            if depth < min_depth:
                min_depth = depth
        # Find the deepest common node for all leaves in taxonomic subtree
        upper_level_taxids = set(UNKNOWN_TAXONOMY_ID)  # top-level LCA is Unknown, not root
        for i in range(0, min_depth+1):
            next_level_taxonomy_ids = set()
            # For each level of taxonomy, find non-redundant list of taxonomy IDs
            for taxonomy_id in taxonomy_id_list:
                if self.is_exist(taxonomy_id):
                    next_level_taxonomy_ids.add(taxonomic_lineages[taxonomy_id][i])
            if len(next_level_taxonomy_ids) > 1:
                # If current level of taxonomy subtree has more than one node,
                # return taxonomy ID of the upper level node. Otherwise,
                # go one level lower
                result = upper_level_taxids.pop()
                break
            upper_level_taxids = next_level_taxonomy_ids
        if len(upper_level_taxids) == 1:
            result = upper_level_taxids.pop()
        return result

    def get_lca2(self, taxonomy_id_list):
        """Returns Lowest Common Ancestor identifier for a list of
        taxonomy identifiers.
        If any taxon in the input list is root, LCA is Unknown.
        If all taxa in the list are Unknowns, LCA is Unknown.
        If any taxon in the list is NOT Unknown, any Unknowns are ignored.
        Note: resulting LCA always has one of ranks defined in RANKS.

        Args:
            taxonomy_id_list (list of str): taxonomy identifiers

        Returns:
            result (str): LCA taxonomy dentifier
        """
        result = UNKNOWN_TAXONOMY_ID

        if len(taxonomy_id_list) == 1:
            taxonomy_id = taxonomy_id_list.pop()
            if taxonomy_id != '':
                result = taxonomy_id
            return result

        taxonomy_id_list = \
            [taxon_id for taxon_id in taxonomy_id_list if taxon_id != UNKNOWN_TAXONOMY_ID]
        if not taxonomy_id_list:
            # if taxonomy_id_list is empty or contains only Unknowns:
            return result

        # Collect taxonomy IDs for each taxonomic rank
        taxonomic_levels = defaultdict(set)
        for taxonomy_id in taxonomy_id_list:
            if self.is_exist(taxonomy_id):
                taxonomic_levels[self.get_rank(taxonomy_id)].add(taxonomy_id)
                parent_id = self.get_upper_level_taxon(taxonomy_id)[0]
                while parent_id != ROOT_TAXONOMY_ID:
                    taxonomic_levels[self.get_rank(parent_id)].add(parent_id)
                    parent_id = self.get_upper_level_taxon(parent_id)[0]
            else:
                print('WARNING: taxonomy ID', taxonomy_id, 'not found in NCBI Taxonomy: skipped')
                continue

        if not taxonomic_levels:
            return result
        # Look for lowest taxonomic rank with one taxonomy ID
        last_good_level = set(UNKNOWN_TAXONOMY_ID)  # top-level LCA is Unknown, not root
        for rank in RANKS[1:]:
            if len(taxonomic_levels[rank]) > 1:
                print(rank, 'is not good!')
                break
            else:
                print(rank, 'is good!')
                last_good_level = taxonomic_levels[rank]
        result = last_good_level.pop()

        return result

    def get_upper_level_taxon(self, taxonomy_id):
        """Finds upper level taxon having rank from RANKS and returns its taxonomy ID and rank

        Args:
            taxonomy_id_list (list of str): taxonomy identifiers

        Returns:
            result (tuple(str, str)): LCA taxonomy dentifier and LCA rank
        """
        print('taxonomy_id', taxonomy_id)
        result = (UNKNOWN_TAXONOMY_ID, self.nodes[UNKNOWN_TAXONOMY_ID]['rank'])
        if not self.is_exist(taxonomy_id):
            return result

        upper_level_id = self.get_parent(taxonomy_id)
        upper_level_rank = self.get_rank(upper_level_id)

        if upper_level_id == UNKNOWN_TAXONOMY_ID:
            pass
        elif taxonomy_id == ROOT_TAXONOMY_ID:
            result = (ROOT_TAXONOMY_ID, self.nodes[ROOT_TAXONOMY_ID]['rank'])
        elif upper_level_rank in RANKS:
            result = (upper_level_id, upper_level_rank)
        else:
            while upper_level_id != ROOT_TAXONOMY_ID:
                upper_level_id = self.get_parent(upper_level_id)
                upper_level_rank = self.get_rank(upper_level_id)
                if upper_level_rank in RANKS:
                    result = (upper_level_id, upper_level_rank)
                    break
#            if result[0] == UNKNOWN_TAXONOMY_ID:
#                result = (ROOT_TAXONOMY_ID, self.nodes[ROOT_TAXONOMY_ID]['rank'])
        return result

    def get_taxonomy_lineage(self, taxonomy_id):
        """Returns taxonomic lineage as concatenated string for a given
        taxonomy identifier. Separator is underscore ('_'), spaces and
        parentheses are replaced with underscores.

        Note: Only taxa with ranks defined in RANKS are reported in the lineage

        Args:
            taxonomy_id (str): taxonomy identifier

        Returns:
            result (str)): taxonomic lineage
        """
        result = ''
        if not self.is_exist(taxonomy_id) or taxonomy_id == ROOT_TAXONOMY_ID:
            return result
        lineage = [self.get_name(taxonomy_id)]
        (parent_id, _) = self.get_upper_level_taxon(taxonomy_id)
        while parent_id != '1':
            lineage.insert(0, self.get_name(parent_id))
            (parent_id, _) = self.get_upper_level_taxon(parent_id)
        result = '_'.join(lineage)
        result = result.replace(' ', '_')
        result = result.replace('(', '_')
        result = result.replace(')', '_')
        return result
