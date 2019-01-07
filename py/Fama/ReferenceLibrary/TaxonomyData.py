#!/usr/bin/python3
import os
from collections import defaultdict,Counter,OrderedDict

class TaxonomyData:
    
    def __new__(cls, options):
        if not hasattr(cls, 'instance'):
            cls.instance = super(TaxonomyData, cls).__new__(cls)
        cls.instance.options = options
        return cls.instance

    def __init__(self,options):
        self.names = defaultdict(dict)
        self.nodes = defaultdict(dict)

    def load_taxdata(self, options):
        names_file = options.get_taxonomy_names_file()
        nodes_file = options.get_taxonomy_nodes_file()
        merged_file = options.get_taxonomy_merged_file()
        
        #initialize self.names
        print ('Loading names file', names_file)
        with open(names_file, 'r') as f:
            for line in f:
                line = line.rstrip('\n\r')
                line_tokens = line.split('\t|\t')
                if line_tokens[3] == 'scientific name\t|':
                    self.names[line_tokens[0]]['name'] = line_tokens[1]
            f.closed
        
        if not self.names:
            raise Exception('Taxonomy names load failed')
        
        

        #initialize self.nodes
        print ('Loading nodes file', nodes_file)
        with open(nodes_file, 'r') as f:
            for line in f:
                line = line.rstrip('\n\r')
                line_tokens = line.split('\t|\t')
                taxid = line_tokens[0]
                parent = line_tokens[1]
                rank = line_tokens[2]
                self.nodes[taxid]['parent'] = parent
                self.nodes[taxid]['rank'] = rank
            f.closed
            
        #merge 
        print ('Loading merged file', merged_file)
        with open(merged_file, 'r') as f:
            for line in f:
                line = line.rstrip('\n\r')
                line_tokens = line.split('\t')
                old_id = line_tokens[0]
                new_id = line_tokens[2]
                if new_id in self.names:
                    self.names[old_id]['name'] = self.names[new_id]['name']
                    self.nodes[old_id]['parent'] = self.nodes[new_id]['parent']
                    self.nodes[old_id]['rank'] = self.nodes[new_id]['rank']
            f.closed
        
        # inject 'Unknown' entry
        self.names['0']['name'] = 'Unknown'
        self.nodes['0']['parent'] = '1'
        self.nodes['0']['rank'] = 'norank'
        
    def get_lca(self, taxonomy_id_list):
        # This function takes list of NCBI Taxonomy IDs and returns ID
        # of the latest common ancestor node in NCBI Taxonomy
        if len(taxonomy_id_list) == 1:
            return taxonomy_id_list[0]

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
                print('WARNING: taxonomy ID',taxonomy_id,'not found in NCBI Taxonomy: skipped')
                continue
            lineage = [parent_id,taxonomy_id]
            while self.nodes[parent_id]['parent'] != '1':
                if self.nodes[parent_id]['parent'] in self.names:
                    parent_id = self.nodes[parent_id]['parent']
                else:
                    parent_id = '0'
                lineage.insert(0,parent_id)
                depth += 1
            taxonomic_lineages[taxonomy_id] = lineage
            if depth < min_depth:
                min_depth = depth
    #    print(taxonomic_lineages)
        # Find the deepest common node for all leaves in taxonomic subtree
        upper_level_taxids = set('0')
        for  i in range(0,min_depth):
            id_set = set()
            # For each level of taxonomy, find non-redundant list of taxonomy IDs
            for taxonomy_id in taxonomy_id_list:
                if taxonomy_id in self.nodes:
                    id_set.add(taxonomic_lineages[taxonomy_id][i])
    #        print (id_set)
            if len(id_set) > 1:
                # If current level of taxonomy subtree has more than one node,
                # return taxonomy ID of the upper level node. Otherwise, 
                # go one level lower
                return upper_level_taxids.pop()
            else:
                upper_level_taxids = id_set
        if len(upper_level_taxids) == 1:
            return upper_level_taxids.pop()
        return '0'
        
    def get_taxonomy_profile(self,counts,identity,scores):
        unknown_label = 'Unknown'
        unknown_rank = 'superkingdom'
        ranks = ['norank', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        
        cellular_organisms_taxid = '131567';
        non_cellular_organisms_name = 'Non-cellular';
        non_cellular_organisms_rank = 'superkingdom';
        
        rpkm_per_rank = defaultdict(lambda : defaultdict(float))
        counts_per_rank = defaultdict(lambda : defaultdict(int))
        identity_per_rank = defaultdict(lambda : defaultdict(float))
        
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
            while current_id != '1':
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
                rpkm_per_rank[non_cellular_organisms_rank][non_cellular_organisms_name] += scores[taxid]
                counts_per_rank[non_cellular_organisms_rank][non_cellular_organisms_name] += counts[taxid]
                identity_per_rank[non_cellular_organisms_rank][non_cellular_organisms_name] += identity[taxid]
                continue
            
            current_id = taxid
            while current_id != '1':
                if current_id not in self.nodes:
                    print('B) Got nothing for ncbi_code in ncbi_nodes: ' + current_id)
                    break
                parent = self.nodes[current_id]['parent']
                rank = self.nodes[current_id]['rank']
                if rank in ranks:
                    name = self.names[current_id]['name']
                    rpkm_per_rank[rank][name] += scores[taxid]
                    counts_per_rank[rank][name] += counts[taxid]
                    identity_per_rank[rank][name] += identity[taxid]
                current_id = parent
        
        for rank in identity_per_rank:
            for taxon in identity_per_rank[rank]:
                identity_per_rank[rank][taxon] = identity_per_rank[rank][taxon]/counts_per_rank[rank][taxon]
        
        return counts_per_rank, identity_per_rank, rpkm_per_rank



