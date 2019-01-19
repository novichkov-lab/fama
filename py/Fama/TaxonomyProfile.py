#!/usr/bin/python
import os
from Fama.Tree import Node,Tree
from collections import defaultdict,Counter,OrderedDict
import pandas as pd


class TaxonomyProfile:
    
    def __init__(self):
        self.tree = Tree()
        self.RANKS = ['norank','superkingdom', 'phylum', 'class', 'order', 'family', 'genus']

    def build_taxonomy_profile(self, taxonomy_data, scores):
        # scores is a dict of dicts of floats, like scores[taxonomy_id][attribute_name] = attribute_value
        # attribute_name can be read_count, identity, rpkm etc.

        unknown_organism_id = '0'
        
        for taxid in sorted(scores.keys()):
            #print('1:', len(scores))
            if taxid == '0':
                #print('2:', len(scores))
                unknown_organisms = Node(rank = 'norank',name = 'Unknown', taxid = '0', parent = '1',children = set())
                self.tree.add_node(unknown_organisms)
                for attribute_key in scores[taxid]:
                    self.tree.add_attribute_recursively(unknown_organism_id,attribute_key,scores[taxid][attribute_key],taxonomy_data)
            elif taxid in taxonomy_data.names:
                #print('3:', len(scores))
                current_id = taxid
                rank = taxonomy_data.nodes[current_id]['rank']
                while True:
                    if rank in self.RANKS:
                        break
                    else:
                        current_id = taxonomy_data.nodes[current_id]['parent']
                        rank = taxonomy_data.nodes[current_id]['rank']
                if self.tree.is_in_tree(current_id):
                    for attribute_key in scores[taxid]:
                        #print('Call add_attribute_recursively', current_id,attribute_key,scores[taxid][attribute_key] )
                        self.tree.add_attribute_recursively(current_id,attribute_key,scores[taxid][attribute_key],taxonomy_data)
                else:
                    parent_taxid = taxonomy_data.nodes[current_id]['parent']
                    label = taxonomy_data.names[current_id]['name']
                    node = Node(rank = rank, name = label, taxid = current_id, parent = parent_taxid, children = set())
                    #print ('Call add_node_recursively for ', taxid, current_id, parent_taxid, label)
                    self.tree.add_node_recursively(node, taxonomy_data)
                    for attribute_key in scores[taxid]:
                        #print('Call add_attribute_recursively', current_id,attribute_key,scores[taxid][attribute_key] )
                        self.tree.add_attribute_recursively(current_id,attribute_key,scores[taxid][attribute_key],taxonomy_data)

            else:
                # Tax ID not in NCBI database
                print('TaxID not found in DB:', taxid)
                node = Node(rank = 'norank', name = taxid, taxid = taxid, parent = unknown_organism_id, children = set())
                self.tree.add_node_recursively(node, taxonomy_data)
                for attribute_key in scores[taxid]:
                    #print('Call add_attribute_recursively',taxid,attribute_key,scores[taxid][attribute_key] )
                    self.tree.add_attribute_recursively(taxid,attribute_key,scores[taxid][attribute_key],taxonomy_data)

    def build_functional_taxonomy_profile(self, taxonomy_data, scores):
        # scores is a dict of dicts of dicts of floats, like scores[taxonomy_id][function_id][attribute_name] = attribute_value
        # attribute_name can be read_count, identity, rpkm etc.

        unknown_organism_id = '0'
        for taxid in sorted(scores.keys()):
            #print('1:', len(scores))
            if taxid == '0':
                #print('2:', len(scores))
                unknown_organisms = Node(rank = 'norank',name = 'Unknown', taxid = unknown_organism_id, parent = '1',children = set())
                self.tree.add_node(unknown_organisms)
                for function in scores[taxid]:
                    self.tree.add_attribute_recursively(taxid,function, scores[taxid][function],taxonomy_data)
            elif taxid in taxonomy_data.names:
                #print('3:', len(scores))
                current_id = taxid
                rank = taxonomy_data.nodes[current_id]['rank']
                while True:
                    if rank in self.RANKS:
                        break
                    else:
                        #print(rank, ' not in RANKS ', taxid)
                        current_id = taxonomy_data.nodes[current_id]['parent']
                        rank = taxonomy_data.nodes[current_id]['rank']
                        
                if self.tree.is_in_tree(current_id):
                    for function in scores[taxid]:
                        self.tree.add_attribute_recursively(current_id,function,scores[taxid][function],taxonomy_data)
                else:
                    parent_taxid = taxonomy_data.nodes[current_id]['parent']
                    label = taxonomy_data.names[current_id]['name']
                    node = Node(rank = rank, name = label, taxid = current_id, parent = parent_taxid, children = set())
                    #print ('Call add_node_recursively for ', taxid, current_id, parent_taxid, label)
                    self.tree.add_node_recursively(node, taxonomy_data)
                    for function in scores[taxid]:
                        #print('Calling add_attribute_recursively', current_id,function,scores[taxid][function] )
                        self.tree.add_attribute_recursively(current_id,function,scores[taxid][function],taxonomy_data)

            else:
                # Tax ID not in NCBI database
                print('TaxID not found in DB:', taxid)
                node = Node(rank = 'norank', name = taxid, taxid = taxid, parent = unknown_organism_id, children = set())
                self.tree.add_node_recursively(node, taxonomy_data)
                for function in scores[taxid]:
                    #print('Call add_attribute_recursively',taxid,attribute_key,scores[taxid][attribute_key] )
                    self.tree.add_attribute_recursively(taxid,function,scores[taxid][function],taxonomy_data)

    def build_assembly_taxonomic_profile(self, taxonomy_data, scores):
        # scores is a dict of dicts of dicts of floats, like scores[taxonomy_id][function_id][attribute_name] = attribute_value
        
        unknown_organism_id = '0'
        for taxid in sorted(scores.keys()):
#            print('1:', len(scores))
            if taxid == '0':
#                print('2:', len(scores))
                unknown_organisms = Node(rank = 'norank',name = 'Unknown', taxid = unknown_organism_id, parent = '1',children = set())
                self.tree.add_node(unknown_organisms)
                for function in scores[taxid]:
                    if 'genes' in scores[taxid][function]:
                        self.tree.add_attribute(taxid,function, {'genes':scores[taxid][function]['genes']}, taxonomy_data)
#                        print('Call add_attribute',taxid,function,scores[taxid][function]['genes'] )
                        attributes = {k:v for k,v in scores[taxid][function].items() if k != 'genes'}
                        self.tree.add_attribute_recursively(taxid,function,attributes,taxonomy_data)
                    else:
                        self.tree.add_attribute_recursively(taxid,function, scores[taxid][function],taxonomy_data)
#                    print('Call add_attribute_recursively',taxid,function,scores[taxid][function] )
            elif taxid in taxonomy_data.names:
#                print('3:', len(scores))
                current_id = taxid
                rank = taxonomy_data.nodes[current_id]['rank']
                while True:
                    if rank in self.RANKS:
                        break
                    elif current_id == '1':
                        break
                    else:
#                        print(rank, ' not in RANKS ', taxid)
                        current_id = taxonomy_data.nodes[current_id]['parent']
                        rank = taxonomy_data.nodes[current_id]['rank']
                        
                if self.tree.is_in_tree(current_id):
                    for function in scores[taxid]:
                        if 'genes' in scores[taxid][function]:
                            self.tree.add_attribute(taxid,function, {'genes':scores[taxid][function]['genes']}, taxonomy_data)
#                            print('Call add_attribute',taxid,function,scores[taxid][function]['genes'] )
                            attributes = {k:v for k,v in scores[taxid][function].items() if k != 'genes'}
                            self.tree.add_attribute_recursively(taxid,function,attributes,taxonomy_data)
                        else:
                            self.tree.add_attribute_recursively(taxid,function, scores[taxid][function],taxonomy_data)
#                        print('Call add_attribute_recursively',taxid,function,scores[taxid][function] )

                else:
                    parent_taxid = taxonomy_data.nodes[current_id]['parent']
                    label = taxonomy_data.names[current_id]['name']
                    node = Node(rank = rank, name = label, taxid = current_id, parent = parent_taxid, children = set())
#                    print ('Call add_node_recursively for ', taxid, current_id, parent_taxid, label)
                    self.tree.add_node_recursively(node, taxonomy_data)
                    for function in scores[taxid]:
                        if 'genes' in scores[taxid][function]:
                            self.tree.add_attribute(taxid,function, {'genes':scores[taxid][function]['genes']}, taxonomy_data)
#                            print('Call add_attribute',taxid,function,scores[taxid][function]['genes'] )
                            attributes = {k:v for k,v in scores[taxid][function].items() if k != 'genes'}
                            self.tree.add_attribute_recursively(taxid,function,attributes,taxonomy_data)
                        else:
                            self.tree.add_attribute_recursively(taxid,function, scores[taxid][function],taxonomy_data)
#                        print('Call add_attribute_recursively',taxid,function,scores[taxid][function] )

            else:
                # Tax ID not in NCBI database
                print('TaxID not found in DB:', taxid)
                node = Node(rank = 'norank', name = taxid, taxid = taxid, parent = unknown_organism_id, children = set())
                self.tree.add_node_recursively(node, taxonomy_data)
                for function in scores[taxid]:
                    if 'genes' in scores[taxid][function]:
                        self.tree.add_attribute(taxid,function, {'genes':scores[taxid][function]['genes']}, taxonomy_data)
#                        print('Call add_attribute',taxid,function,scores[taxid][function]['genes'] )
                        attributes = {k:v for k,v in scores[taxid][function].items() if k != 'genes'}
                        self.tree.add_attribute_recursively(taxid,function,attributes,taxonomy_data)
                    else:
                        self.tree.add_attribute_recursively(taxid,function, scores[taxid][function],taxonomy_data)

#                    print('Call add_attribute_recursively',taxid,function,scores[taxid][function] )
                    #self.tree.add_attribute_recursively(taxid,function,scores[taxid][function],taxonomy_data)

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
                if rank in self.RANKS:
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
                for attribute in self.tree.data[taxid].attributes:
                    ret_val += '\t'+attribute+':'+str(self.tree.data[taxid].attributes[attribute])
                ret_val += '\n'
                #ret_val += '\tScore:' + format((self.tree.data[taxid].attributes['rpkm']), "0.3f") + '\tIdentity:' + format((self.tree.data[taxid].attributes['identity']/self.tree.data[taxid].attributes['count']), "0.1f") + '%\tRead count:' + format(self.tree.data[taxid].attributes['count'], "0.0f") + '\n'
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

    def print_functional_taxonomy_profile(self, score='rpkm'):
        root_id = '1'
        offset = 0
        ret_val = self.print_node_functions(root_id, offset, score)
        return ret_val
        
    def print_node_functions(self, taxid, offset, score='rpkm'):
        if taxid in self.tree.data:
            ret_val = '\t'*offset + taxid + '\t' + self.tree.data[taxid].rank + '\t' + self.tree.data[taxid].name + '\tParent:' + (self.tree.data[taxid].parent or 'None') + '\tChildren:'+(','.join(self.tree.data[taxid].children) or 'None') + '\n'
            if self.tree.data[taxid].attributes:
                for function in self.tree.data[taxid].attributes:
                    try:
                        ret_val += '\t'*(offset + 5) + function + '\tScore:' + format((self.tree.data[taxid].attributes[function][score]), "0.3f") + '\tIdentity:' + format((self.tree.data[taxid].attributes[function]['identity']/self.tree.data[taxid].attributes[function]['hit_count']), "0.1f") + '%\tRead count:' + format(self.tree.data[taxid].attributes[function]['count'], "0.0f") + '\n'
                    except TypeError:
                        print(function, taxid, self.tree.data[taxid].attributes)
            else:
                ret_val += '\t'*(offset + 6) + '\tScore:N/A\tIdentity:N/A\tRead count:N/A\n'
            offset += 1
            if self.tree.data[taxid].children:
                for child_id in sorted(self.tree.data[taxid].children):
                    ret_val += self.print_node_functions(child_id,offset)
            else:
                print('Node has no children:' + taxid + '\t' + self.tree.data[taxid].rank + '\t' + self.tree.data[taxid].name)
            #print(ret_val)
            return ret_val
        else:
            print('Node not found:',taxid)
            return ''

    def print_functional_taxonomy_table(self):
        function_list = set()
        for taxid in self.tree.data:
            for function in self.tree.data[taxid].attributes.keys():
                function_list.add(function) 
        
        ret_val = list()
        ret_val.append('No.\tRank\tName\t' + '\t'.join(x+'\t\t' for x in sorted(function_list)))
        ret_val.append('\t\t' + '\tScore\tIdentity\tRead count'*len(function_list))

        root_id = '1'
        line_number = 1
        ret_val.extend(self.print_node_function_table(root_id, function_list, line_number))
        return '\n'.join(ret_val)

    def print_node_function_table(self, taxid, function_list, line_number):
        ret_val = list()
        if taxid in self.tree.data:
            line = [str(line_number), '\t', self.tree.data[taxid].rank, '\t', self.tree.data[taxid].name]
            if self.tree.data[taxid].attributes:
                for function in function_list:
                    if function in self.tree.data[taxid].attributes:
                        line.extend(['\t', format((self.tree.data[taxid].attributes[function]['rpkm']), "0.3f"), '\t', format((self.tree.data[taxid].attributes[function]['identity']/self.tree.data[taxid].attributes[function]['count']), "0.1f"), '%\t', format(self.tree.data[taxid].attributes[function]['count'], "0.0f")])
                    else:
                        line.append('\t'*3)
            else:
                line.append('\t'*len(function_list)*3)
            ret_val.append(''.join(line))
            line_number += 1
            if self.tree.data[taxid].children:
                for child_id in sorted(self.tree.data[taxid].children):
                    children_lines = self.print_node_function_table(child_id, function_list, line_number)
                    ret_val.extend(children_lines)
                    line_number += len(children_lines)
            else:
                print('Node has no children:', taxid, '\t', self.tree.data[taxid].rank, '\t', self.tree.data[taxid].name)
            #print(ret_val)
            return ret_val
        else:
            print('Node not found:',taxid)
            return ret_val

    def convert_function_taxonomic_profile_into_df(self, score='rpkm'):
        function_list = set()
        for taxid in self.tree.data:
            for function in self.tree.data[taxid].attributes.keys():
                function_list.add(function) 
        root_id = '1'
        line_number = 1
        df = pd.DataFrame(self.convert_node_into_dict(root_id, function_list, line_number, score=score))
        return df.transpose()

    def convert_node_into_dict(self, taxid, function_list, line_number, score='rpkm'):
        ret_val = {}
        line_dict = {}
        if taxid in self.tree.data:
            #line_dict[('', 'No.')] = str(line_number)
            line_dict[('', 'Rank')] = self.tree.data[taxid].rank
            line_dict[('', 'Name')] = self.tree.data[taxid].name
            for function in function_list:
                if function in self.tree.data[taxid].attributes:
                    line_dict[(function, '1.Score')] = format((self.tree.data[taxid].attributes[function][score]), "0.3f")
                    if 'identity' in self.tree.data[taxid].attributes[function]:
                        line_dict[(function, '2.Identity')] = format((self.tree.data[taxid].attributes[function]['identity']/self.tree.data[taxid].attributes[function]['hit_count']), "0.1f") + '%'
                    else:
                        line_dict[(function, '2.Identity')] = ''
                    line_dict[(function, '3.Read count')] = format(self.tree.data[taxid].attributes[function]['count'], "0.0f")
                else:
                    line_dict[(function, '1.Score')] = ''
                    line_dict[(function, '2.Identity')] = ''
                    line_dict[(function, '3.Read count')] = ''
            ret_val[line_number] = line_dict
            line_number += 1
            if self.tree.data[taxid].children:
                for child_id in sorted(self.tree.data[taxid].children):
                    children_lines = self.convert_node_into_dict(child_id, function_list, line_number, score)
                    for child_line_number in children_lines:
                        ret_val[child_line_number] = children_lines[child_line_number]
                    line_number += len(children_lines)
            else:
                #print('Node has no children:', taxid, '\t', self.tree.data[taxid].rank, '\t', self.tree.data[taxid].name)
                pass
            #print(ret_val)
            return ret_val
        else:
            print('Node not found:',taxid)
            return ret_val
