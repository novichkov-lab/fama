#!/usr/bin/python
import os
from Fama.Tree import Node,Tree
from collections import defaultdict,Counter,OrderedDict

RANKS = ['norank','superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

class TaxonomyProfile:
    
    def __init__(self):
        self.tree = Tree()

    def build_taxonomy_profile(self, taxonomy_data, scores):
        # scores is a dict of dicts of floats, like scores[taxonomy_id][attribute_name] = attribute_value

        unknown_organism_id = '0'
        
        for taxid in sorted(scores.keys()):
            #print('1:', len(scores))
            if taxid == '0':
                #print('2:', len(scores))
                self.tree.data[taxid].add_attribute_recursively(unknown_organism_id,'score',scores[taxid])
            elif taxid in taxonomy_data.names:
                #print('3:', len(scores))
                current_id = taxid
                rank = taxonomy_data.nodes[current_id]['rank']
                while True:
                    if rank in RANKS:
                        break
                    else:
                        #print(rank, ' not in RANKS ', taxid)
                        current_id = taxonomy_data.nodes[current_id]['parent']
                        rank = taxonomy_data.nodes[current_id]['rank']
                        
                if self.tree.is_in_tree(current_id):
                    for attribute_key in scores[taxid]:
                        #print('Call add_attribute_recursively', current_id,attribute_key,scores[taxid][attribute_key] )
                        self.tree.add_attribute_recursively(current_id,attribute_key,scores[taxid][attribute_key])
                else:
                    parent_taxid = taxonomy_data.nodes[current_id]['parent']
                    label = taxonomy_data.names[current_id]['name']
                    node = Node(rank = rank, name = label, taxid = current_id, parent = parent_taxid, children = set())
                    #print ('Call add_node_recursively for ', taxid, current_id, parent_taxid, label)
                    self.tree.add_node_recursively(node, taxonomy_data)
                    for attribute_key in scores[taxid]:
                        #print('Call add_attribute_recursively', current_id,attribute_key,scores[taxid][attribute_key] )
                        self.tree.add_attribute_recursively(current_id,attribute_key,scores[taxid][attribute_key])

            else:
                # Tax ID not in NCBI database
                print('TaxID not found in DB:', taxid)
                node = Node(rank = 'norank', name = taxid, taxid = taxid, parent = unknown_organism_id, children = set())
                self.tree.add_node_recursively(node, taxonomy_data)
                for attribute_key in scores[taxid]:
                    #print('Call add_attribute_recursively',taxid,attribute_key,scores[taxid][attribute_key] )
                    self.tree.add_attribute_recursively(taxid,attribute_key,scores[taxid][attribute_key])


    def get_taxonomy_profile(self,counts,identity,scores):
        unknown_label = 'Unknown'
        unknown_rank = 'superkingdom'
        
        
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
                if rank in RANKS:
                    name = self.names[current_id]['name']
                    rpkm_per_rank[rank][name] += scores[taxid]
                    counts_per_rank[rank][name] += counts[taxid]
                    identity_per_rank[rank][name] += identity[taxid]
                current_id = parent
        
        for rank in identity_per_rank:
            for taxon in identity_per_rank[rank]:
                identity_per_rank[rank][taxon] = identity_per_rank[rank][taxon]/counts_per_rank[rank][taxon]
        
        return counts_per_rank, identity_per_rank, rpkm_per_rank
    
    def print_taxonomy_profile(self):
        root_id = '1'
        offset = 0
        ret_val = self.print_node(root_id, offset)
        return ret_val
        
    def print_node(self, taxid, offset):
        if taxid in self.tree.data:
            ret_val = '\t'*offset + taxid + '\t' + self.tree.data[taxid].rank + '\t' + self.tree.data[taxid].name + '\tParent:' + (self.tree.data[taxid].parent or 'None') + '\tChildren:'+(','.join(self.tree.data[taxid].children) or 'None')
            if self.tree.data[taxid].attributes:
                ret_val += '\tScore:' + format((self.tree.data[taxid].attributes['rpkm']), "0.2f") + '\tIdentity:' + format((self.tree.data[taxid].attributes['identity']/self.tree.data[taxid].attributes['count']), "0.1f") + '%\tRead count:' + format(self.tree.data[taxid].attributes['count'], "0.0f") + '\n'
            else:
                ret_val += '\tScore:N/A\tIdentity:N/A\tRead count:N/A\n'
            offset += 1
            if self.tree.data[taxid].children:
                for child_id in sorted(self.tree.data[taxid].children):
                    ret_val += self.print_node(child_id,offset)
            else:
                print('Node has no children:' + taxid + '\t' + self.tree.data[taxid].rank + '\t' + self.tree.data[taxid].name)
            print(ret_val)
            return ret_val
        else:
            print('Node not found:',taxid)
            return ''
        



