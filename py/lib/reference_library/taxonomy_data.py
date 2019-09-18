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
        """Returns Lowest Common Ancestor identifier for a list of taxonomy identifiers"""
        result = UNKNOWN_TAXONOMY_ID
        if len(taxonomy_id_list) == 1:
            taxonomy_id = taxonomy_id_list.pop()
            if taxonomy_id != '':
                result = taxonomy_id
            return result

        taxonomic_lineages = {}
        # Calculate length of the shortest path in taxonomic subtree
        min_depth = 1000
        for taxonomy_id in taxonomy_id_list:
            depth = 1
            if taxonomy_id in self.nodes:
                parent_id = self.nodes[taxonomy_id]['parent']
            elif taxonomy_id == '':
                continue
            else:
                print('WARNING: taxonomy ID', taxonomy_id, 'not found in NCBI Taxonomy: skipped')
                continue
            lineage = [parent_id, taxonomy_id]
            while self.nodes[parent_id]['parent'] != ROOT_TAXONOMY_ID:
                if self.nodes[parent_id]['parent'] in self.names:
                    parent_id = self.nodes[parent_id]['parent']
                else:
                    parent_id = UNKNOWN_TAXONOMY_ID
                lineage.insert(0, parent_id)
                depth += 1
#            print(lineage)
            taxonomic_lineages[taxonomy_id] = lineage
            if depth < min_depth:
                min_depth = depth
#        print(taxonomic_lineages)
        # Find the deepest common node for all leaves in taxonomic subtree
        upper_level_taxids = set(UNKNOWN_TAXONOMY_ID)
        for i in range(0, min_depth+1):
            id_set = set()
            # For each level of taxonomy, find non-redundant list of taxonomy IDs
            for taxonomy_id in taxonomy_id_list:
                if taxonomy_id in self.nodes:
                    id_set.add(taxonomic_lineages[taxonomy_id][i])
#            print (id_set)
            if len(id_set) > 1:
                # If current level of taxonomy subtree has more than one node,
                # return taxonomy ID of the upper level node. Otherwise,
                # go one level lower
                result = upper_level_taxids.pop()
                break
            else:
                upper_level_taxids = id_set
        if len(upper_level_taxids) == 1:
            result = upper_level_taxids.pop()
        return result

    def get_lca2(self, taxonomy_id_list):
        """Returns Lowest Common Ancestor identifier for a list of
        taxonomy identifiers.
        The LCA always has one of ranks defined in RANKS.
        """
        ret_val = UNKNOWN_TAXONOMY_ID

        if len(taxonomy_id_list) == 1:
            taxonomy_id = taxonomy_id_list.pop()
            if taxonomy_id != '':
                ret_val = taxonomy_id
            return ret_val

        # Collect taxonomy IDs for each taxonomic rank
        taxonomic_levels = defaultdict(set)
        for taxonomy_id in taxonomy_id_list:
            if taxonomy_id == '':
                continue
            if taxonomy_id in self.nodes:
                taxonomic_levels[self.nodes[taxonomy_id]['rank']].add(taxonomy_id)
                parent_id = self.nodes[taxonomy_id]['parent']
                while parent_id != ROOT_TAXONOMY_ID:
                    taxonomic_levels[self.nodes[parent_id]['rank']].add(parent_id)
                    parent_id = self.nodes[parent_id]['parent']
            else:
                print('WARNING: taxonomy ID', taxonomy_id, 'not found in NCBI Taxonomy: skipped')
                continue

        if not taxonomic_levels:
            return ret_val
        # Look for lowest taxonomic rank with one taxonomy ID
        last_good_level = set(ROOT_TAXONOMY_ID)
        for rank in RANKS[1:]:
            if len(taxonomic_levels[rank]) > 1:
                print(rank, 'is not good!')
                break
            else:
                print(rank, 'is good!')
                last_good_level = taxonomic_levels[rank]
        ret_val = last_good_level.pop()

        # Additional check if LCA has acceptable rank
        while ret_val != ROOT_TAXONOMY_ID:
            print('LCA', ret_val)
            if self.nodes[ret_val]['rank'] in RANKS:
                break
            ret_val = self.nodes[ret_val]['parent']

        return ret_val

    def get_taxonomy_profile(self, counts, identity, scores):
        """Calculates taxonomy profile for a number of taxonomy identifiers.
        This function takes three dictionaries, assuming that all of
        them have equal size and identical keys. This function is used only for
        generation of text and PDF reports for individual FASTQ/FASTA files.

        Args:
            counts (dict[str,int]): key is taxonomy identifier, value is count
            identity (dict[str,float]): key is taxonomy identifier, value is amino acid % identity
            scores (dict[str,float]): key is taxonomy identifier, value is score

        Returns:
            counts_per_rank (defaultdict[str, defaultdict[str, float]]):
                external key is rank, internal key is name, value is count
            identity_per_rank (defaultdict[str, defaultdict[str, float]]):
                external key is rank, internal key is name, value is amino acid % identity
            scores_per_rank (defaultdict[str, defaultdict[str, float]]): ():
                external key is rank, internal key is name, value is score
        """
        unknown_label = 'Unknown'
        unknown_rank = 'superkingdom'
        cellular_organisms_taxid = '131567'
        non_cellular_organisms_name = 'Non-cellular'
        non_cellular_organisms_rank = 'superkingdom'

        rpkm_per_rank = defaultdict(lambda: defaultdict(float))
        counts_per_rank = defaultdict(lambda: defaultdict(int))
        identity_per_rank = defaultdict(lambda: defaultdict(float))

        for taxid in counts:
            current_id = taxid
            if taxid == 0:
                label = unknown_label
                rpkm_per_rank[unknown_rank][label] += scores[taxid]
                counts_per_rank[unknown_rank][label] += counts[taxid]
                identity_per_rank[unknown_rank][label] += identity[taxid]
                continue
            is_cellular = False
            not_found = False
            while current_id != ROOT_TAXONOMY_ID:
                if current_id == cellular_organisms_taxid:
                    is_cellular = True
                    break
                if current_id not in self.nodes:
                    print('A) ncbi_code not found in ncbi_nodes: \'' + current_id + '\'')
                    not_found = True
                    break
                current_id = self.nodes[current_id]['parent']

            if not_found:
                continue

            if not is_cellular:
                rpkm_per_rank[non_cellular_organisms_rank][non_cellular_organisms_name] \
                    += scores[taxid]
                counts_per_rank[non_cellular_organisms_rank][non_cellular_organisms_name] \
                    += counts[taxid]
                identity_per_rank[non_cellular_organisms_rank][non_cellular_organisms_name] \
                    += identity[taxid]
                continue

            current_id = taxid
            while current_id != ROOT_TAXONOMY_ID:
                if current_id not in self.nodes:
                    print('B) Got nothing for ncbi_code in ncbi_nodes: ' + current_id)
                    break
                parent = self.nodes[current_id]['parent']
                rank = self.nodes[current_id]['rank']
                if rank in RANKS:
                    name = self.names[current_id]['name']
                    rpkm_per_rank[rank][name] += scores[taxid]
                    counts_per_rank[rank][name] += counts[taxid]
                    identity_per_rank[rank][name] += identity[taxid]
                current_id = parent

        for rank in identity_per_rank:
            for taxon in identity_per_rank[rank]:
                identity_per_rank[rank][taxon] = \
                    identity_per_rank[rank][taxon]/counts_per_rank[rank][taxon]

        return counts_per_rank, identity_per_rank, rpkm_per_rank

    def get_upper_level_taxon(self, taxonomy_id):
        """Finds upper level taxon having rank from RANKS and returns its taxonomy ID"""
        result = UNKNOWN_TAXONOMY_ID, self.nodes[UNKNOWN_TAXONOMY_ID]['rank']
        if taxonomy_id not in self.names:
            return result

        current_id = self.nodes[taxonomy_id]['parent']
        current_rank = self.nodes[current_id]['rank']

        if current_id == UNKNOWN_TAXONOMY_ID:
            pass
        elif taxonomy_id == ROOT_TAXONOMY_ID:
            result = ROOT_TAXONOMY_ID, self.nodes[ROOT_TAXONOMY_ID]['rank']
        elif current_rank in RANKS:
            result = current_id, current_rank
        else:
            while current_id != '1':
                current_id = self.nodes[current_id]['parent']
                current_rank = self.nodes[current_id]['rank']
                if current_rank in RANKS:
                    result = current_id, current_rank
            if result[0] == UNKNOWN_TAXONOMY_ID:
                result = ROOT_TAXONOMY_ID, self.nodes[ROOT_TAXONOMY_ID]['rank']
        return result

    def get_lineage_string(self, taxonomy_id):
        """Returns taxonomic lineage as concatenated string for a given
        taxonomy identifier

        Todo: replace with get_taxonomy_lineage
        """
        ret_val = ''
        if taxonomy_id not in self.nodes:
            return ret_val

        lineage = [self.names[taxonomy_id]['name']]
        parent_id = self.nodes[taxonomy_id]['parent']
        while self.nodes[parent_id]['parent'] != '1':
            if self.nodes[parent_id]['rank'] in RANKS:
                lineage.insert(0, self.names[parent_id]['name'])
            if self.nodes[parent_id]['parent'] in self.names:
                parent_id = self.nodes[parent_id]['parent']
            else:
                parent_id = '0'
        ret_val = '_'.join(lineage)
        ret_val = ret_val.replace(' ', '_')
        ret_val = ret_val.replace('(', '_')
        ret_val = ret_val.replace(')', '_')
        return ret_val

    def get_taxonomy_lineage(self, taxonomy_id):
        """Returns taxonomic lineage as concatenated string for a given
        taxonomy identifier
        """
        ret_val = ''
        if taxonomy_id not in self.nodes:
            return ret_val
        lineage = [self.names[taxonomy_id]['name']]
        parent_id = self.nodes[taxonomy_id]['parent']
        while self.nodes[parent_id]['parent'] != '1':
            if self.nodes[parent_id]['rank'] in RANKS:
                lineage.insert(0, self.names[parent_id]['name'])
            if self.nodes[parent_id]['parent'] in self.names:
                parent_id = self.nodes[parent_id]['parent']
            else:
                parent_id = '0'
        ret_val = '_'.join(lineage)
        ret_val = ret_val.replace(' ', '_')
        ret_val = ret_val.replace('(', '_')
        ret_val = ret_val.replace(')', '_')
        return ret_val
