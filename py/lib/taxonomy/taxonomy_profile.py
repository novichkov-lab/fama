"""Describes TaxonomyProfile class"""
from collections import defaultdict
import pandas as pd
from lib.utils.const import RANKS, LOWER_RANKS, UNKNOWN_TAXONOMY_ID, ROOT_TAXONOMY_ID
from lib.utils.utils import autovivify
from lib.taxonomy.tree import Node, Tree


class TaxonomyProfile(object):
    """TaxonomyProfile object performs various operations on phylogenetic tree

    Attributes:
        tree (:obj:Tree): phylogenetic tree
    """
    def __init__(self):
        """Init creates a new Tree instance"""
        self.tree = Tree()

    def make_function_taxonomy_profile(self, taxonomy_data, scores):
        """Builds functional taxonomy profile, with two-level node attributes.

        Args:
            taxonomy_data (:obj:TaxonomyData): taxonomic data instance
            scores (dict[str, dict[str, dict[str, obj]]]): outer key is taxonomy identifier,
                mid-level key is function identifier, inner key is attribute name,
                value is attribute value (usually, float or int).
                Attribute keys can be read_count, identity, rpkm etc.
        """
        for taxid in sorted(scores.keys()):
            if taxid == UNKNOWN_TAXONOMY_ID:
                unknown_organisms = Node(rank='norank',
                                         name='Unknown',
                                         taxid=UNKNOWN_TAXONOMY_ID,
                                         parent=ROOT_TAXONOMY_ID)
                self.tree.add_node(unknown_organisms)
                for function in scores[taxid]:
                    self.tree.add_attribute_recursively(taxid, function,
                                                        scores[taxid][function],
                                                        taxonomy_data)
            elif taxonomy_data.is_exist(taxid):
                current_id = taxid
                rank = taxonomy_data.get_rank(current_id)
                # If taxonomical rank ot taxid is not in RANKS, find parent
                # taxon with acceptable rank
                while True:
                    if rank in RANKS:
                        break
                    else:
                        current_id = taxonomy_data.get_parent(current_id)
                        rank = taxonomy_data.get_rank(current_id)
                # Add new node if current_id is not in the tree
                if not self.tree.is_in_tree(current_id):
                    parent_taxid = taxonomy_data.get_parent(current_id)
                    label = taxonomy_data.get_name(current_id)
                    node = Node(rank=rank, name=label, taxid=current_id, parent=parent_taxid)
                    self.tree.add_node_recursively(node, taxonomy_data)
                for function in scores[taxid]:
                    self.tree.add_attribute_recursively(current_id, function,
                                                        scores[taxid][function],
                                                        taxonomy_data)
            else:
                # Tax ID not in NCBI database
                print('TaxID not found in DB:', taxid)
                node = Node(rank='norank', name=taxid, taxid=taxid, parent=UNKNOWN_TAXONOMY_ID)
                self.tree.add_node_recursively(node, taxonomy_data)
                for function in scores[taxid]:
                    self.tree.add_attribute_recursively(taxid, function,
                                                        scores[taxid][function],
                                                        taxonomy_data)

    def make_assembly_taxonomy_profile(self, taxonomy_data, scores):
        """Builds taxonomic profile of gene assembly

        Args:
            taxonomy_data (:obj:TaxonomyData): taxonomic data instance
            scores (dict[str, dic[str, dict[str, obj]]]): outer key is taxonomy identifier,
                mid-level key is function identifier, inner key is attribute name,
                value is attribute value (usually, float or int).
                Attribute keys can be read_count, identity, rpkm etc.
        """
        for taxid in sorted(scores.keys()):
            if taxid == UNKNOWN_TAXONOMY_ID:
                unknown_organisms = Node(rank='norank', name='Unknown',
                                         taxid=UNKNOWN_TAXONOMY_ID,
                                         parent=ROOT_TAXONOMY_ID)
                self.tree.add_node(unknown_organisms)
                for function in scores[taxid]:
                    if 'genes' in scores[taxid][function]:
                        self.tree.add_attribute(
                            taxid, function, {'genes': scores[taxid][function]['genes']},
                            taxonomy_data
                        )
                        attributes = {
                            k: v for k, v in scores[taxid][function].items() if k != 'genes'
                        }
                        self.tree.add_attribute_recursively(
                            taxid, function, attributes, taxonomy_data
                        )
                    else:
                        self.tree.add_attribute_recursively(
                            taxid, function, scores[taxid][function], taxonomy_data
                        )
            elif taxonomy_data.is_exist(taxid):
                # If taxonomical rank ot taxid is not in RANKS, find parent
                # taxon with acceptable rank
                current_id = taxid
                rank = taxonomy_data.get_rank(current_id)
                while True:
                    if rank in RANKS:
                        break
                    else:
                        current_id = taxonomy_data.get_parent(current_id)
                        rank = taxonomy_data.get_rank(current_id)

                if not self.tree.is_in_tree(current_id):
                    parent_taxid = taxonomy_data.get_parent(current_id)
                    label = taxonomy_data.get_name(current_id)
                    node = Node(rank=rank, name=label, taxid=current_id, parent=parent_taxid)
                    self.tree.add_node_recursively(node, taxonomy_data)

                for function in scores[taxid]:
                    if 'genes' in scores[taxid][function]:
                        self.tree.add_attribute(
                            taxid, function, {'genes': scores[taxid][function]['genes']},
                            taxonomy_data
                        )
                        attributes = {
                            k: v for k, v in scores[taxid][function].items() if k != 'genes'
                        }
                        self.tree.add_attribute_recursively(taxid, function,
                                                            attributes,
                                                            taxonomy_data)
                    else:
                        self.tree.add_attribute_recursively(
                            taxid, function, scores[taxid][function], taxonomy_data
                        )
            else:
                # Tax ID not in NCBI database. Create a child node for Unknown
                print('TaxID not found in DB:', taxid)
                node = Node(rank='norank', name=taxid, taxid=taxid, parent=UNKNOWN_TAXONOMY_ID)
                self.tree.add_node_recursively(node, taxonomy_data)
                for function in scores[taxid]:
                    if 'genes' in scores[taxid][function]:
                        self.tree.add_attribute(
                            taxid, function, {'genes': scores[taxid][function]['genes']},
                            taxonomy_data
                        )
                        attributes = {
                            k: v for k, v in scores[taxid][function].items() if k != 'genes'
                        }
                        self.tree.add_attribute_recursively(
                            taxid, function, attributes, taxonomy_data
                        )
                    else:
                        self.tree.add_attribute_recursively(
                            taxid, function, scores[taxid][function], taxonomy_data
                        )

    #  def get_taxonomy_profile(self, counts, identity, scores, taxonomy_data):
        #  """Returns scores for all levels of taxonomy profile.

        #  Args:
            #  counts (:obj:dict[str, int]): key is taxonomy identifier, value is
                #  raw read count
            #  identity (:obj:dict[str, int]): key is taxonomy identifier, value is
                #  amino acid identity %
            #  scores (:obj:dict[str, int]): key is taxonomy identifier, value is
                #  score

        #  Returns:
            #  counts_per_rank [defaultdict[str, defaultdict[str, int]]]: outer
                #  key is taxonomic rank, inner key is taxon name, value is raw
                #  read count
            #  identity_per_rank [defaultdict[str, defaultdict[str, int]]]: outer
                #  key is taxonomic rank, inner key is taxon name, value is amino
                #  acid identity %
            #  rpkm_per_rank [defaultdict[str, defaultdict[str, int]]]: outer
                #  key is taxonomic rank, inner key is taxon name, value is RPKM score
        #  """
        #  unknown_label = 'Unknown'
        #  unknown_rank = 'norank'

        #  cellular_organisms_taxid = '131567'
        #  non_cellular_organisms_name = 'Non-cellular'
        #  non_cellular_organisms_rank = 'superkingdom'

        #  rpkm_per_rank = defaultdict(lambda: defaultdict(float))
        #  counts_per_rank = defaultdict(lambda: defaultdict(int))
        #  identity_per_rank = defaultdict(lambda: defaultdict(float))

        #  for taxid in counts:
            #  current_id = taxid
            #  if taxid == UNKNOWN_TAXONOMY_ID:
                #  label = unknown_label
                #  rpkm_per_rank[unknown_rank][label] += scores[taxid]
                #  counts_per_rank[unknown_rank][label] += counts[taxid]
                #  identity_per_rank[unknown_rank][label] += identity[taxid]
                #  continue
            #  is_cellular = False
            #  not_found = False
            #  while current_id != ROOT_TAXONOMY_ID:
                #  if current_id == cellular_organisms_taxid:
                    #  is_cellular = True
                    #  break
                #  if current_id not in taxonomy_data.nodes:
                    #  print('A) ncbi_code not found in ncbi_nodes: \'' + current_id + '\'')
                    #  not_found = True
                    #  break
                #  current_id = taxonomy_data.nodes[current_id]['parent']

            #  if not_found:
                #  continue

            #  if not is_cellular:
                #  rpkm_per_rank[non_cellular_organisms_rank][non_cellular_organisms_name] += \
                    #  scores[taxid]
                #  counts_per_rank[non_cellular_organisms_rank][non_cellular_organisms_name] += \
                    #  counts[taxid]
                #  identity_per_rank[non_cellular_organisms_rank][non_cellular_organisms_name] += \
                    #  identity[taxid]
                #  continue

            #  current_id = taxid
            #  while current_id != '1':
                #  if current_id not in taxonomy_data.nodes:
                    #  print('B) Got nothing for ncbi_code in ncbi_nodes: ' + current_id)
                    #  break
                #  parent = taxonomy_data.nodes[current_id]['parent']
                #  rank = taxonomy_data.nodes[current_id]['rank']
                #  if rank in RANKS:
                    #  name = taxonomy_data.names[current_id]['name']
                    #  rpkm_per_rank[rank][name] += scores[taxid]
                    #  counts_per_rank[rank][name] += counts[taxid]
                    #  identity_per_rank[rank][name] += identity[taxid]
                #  current_id = parent

        #  for rank in identity_per_rank:
            #  for taxon in identity_per_rank[rank]:
                #  identity_per_rank[rank][taxon] = identity_per_rank[rank][taxon] / \
                    #  counts_per_rank[rank][taxon]

        #  return counts_per_rank, identity_per_rank, rpkm_per_rank

    def __str__(self):
        """Returns taxonomy profile as text: all nodes starting from root
        (top-down depth-first) as tab-separated fields, with all node attributes.

        Returns:
            ret_val(str): string representation of all nodes in the tree
        """
        root_id = ROOT_TAXONOMY_ID
        offset = 0
        ret_val = self.stringify_node(root_id, offset)
        return ret_val

    def stringify_node(self, taxid, offset):
        """Returns string representation of node (tab-separated), which includes
        taxonomy identifier, rank, taxon name, parent identifier, children
        identifiers, stringified attributes. Recursively called for all children
        of the node.

        Args:
            taxid (str): taxonomy identifier of node
            offset (int): number of tabs in the beginning

        Returns:
            ret_val (str): string representation of node
        """
        ret_val = ''
        if taxid in self.tree.data:
            ret_val = '\t'*offset + taxid + '\t' + self.tree.data[taxid].rank \
                + '\t' + self.tree.data[taxid].name + '\tParent:' \
                + (self.tree.data[taxid].parent or 'None') + '\tChildren:' + \
                (','.join(self.tree.data[taxid].children) or 'None')
            if self.tree.data[taxid].attributes:
                for attribute in self.tree.data[taxid].attributes:
                    ret_val += '\t' + attribute + ':' + \
                        str(self.tree.data[taxid].attributes[attribute])
                ret_val += '\n'
            else:
                ret_val += '\tScore:N/A\tIdentity:N/A\tRead count:N/A\n'
            offset += 1
            if self.tree.data[taxid].children:
                for child_id in sorted(self.tree.data[taxid].children):
                    ret_val += self.stringify_node(child_id, offset)
        else:
            ret_val = ''
        return ret_val

    def convert_profile_into_df(self, metric='rpkm'):
        """Converts functional-taxonomic profile into pandas DataFrame object

        Args:
            metric (str): score metric (default value 'rpkm') to be reported

        Returns:
            result (pd.DataFrame): functional-taxonomic profile as pandas
                DataFrame object
        """
        function_list = set()
        for taxid in self.tree.data:
            for function in self.tree.data[taxid].attributes.keys():
                function_list.add(function)
        root_id = ROOT_TAXONOMY_ID
        line_number = 1
        tax_dict, _ = self.convert_node_into_dict(
            root_id, function_list, line_number, metric=metric
            )
        result = pd.DataFrame(tax_dict)
        result = result.transpose()
        return result

    def convert_node_into_dict(self, taxid, function_list, line_number, metric='rpkm'):
        """Returns node of functional-taxonomic profile for conversion into DataFrame.
        Recursively called for all children of the node.

        Args:
            taxid (str): taxonomy identifier of node
            function_list (list of str): function identifiers to be included
                to the table
            line_number (int): sequential number of node printed
            metric (str): score metric (default value 'rpkm') to be reported

        Returns:
            ret_val (dict[str,dict[tuple(str,str),float]]): outer key is line number,
                inner key is a tuple with function identifier or empty string as
                first element and field name as second element, value is a float.
                Field names are 'Rank', 'Taxon name', '1.Score', '2.Identity',
                '3.Read count'. For each function, only the latter three fields are reported.
            attribute_values (defaultdict[str,dict[str,obj]]): outer key is
                function identifier, inner key may be metric, 'hit_count', 'identity'
                'hit_count', value is float.
        """
        # Collect values of all required attributes for reporting to the upper level
        attribute_values = defaultdict(dict)
        ret_val = defaultdict(dict)
        if taxid not in self.tree.data:
            return ret_val, attribute_values
        for function in function_list:
            for attribute_name in ['count', 'identity', 'hit_count', metric]:
                attribute_values[function][attribute_name] = 0.0
            if function in self.tree.data[taxid].attributes:
                if metric in self.tree.data[taxid].attributes[function]:
                    attribute_values[function][metric] = \
                        self.tree.data[taxid].attributes[function][metric]
                if 'count' in self.tree.data[taxid].attributes[function]:
                    attribute_values[function]['count'] = \
                        self.tree.data[taxid].attributes[function]['count']
                if 'identity' in self.tree.data[taxid].attributes[function]:
                    attribute_values[function]['identity'] = \
                        self.tree.data[taxid].attributes[function]['identity']
                    attribute_values[function]['hit_count'] = \
                        self.tree.data[taxid].attributes[function]['hit_count']

        children_values = autovivify(2, float)

        ret_val[line_number][('', 'Rank')] = self.tree.data[taxid].rank
        ret_val[line_number][('', 'Taxon name')] = self.tree.data[taxid].name
        for function in function_list:
            for field_name in ['1.Score', '2.Identity', '3.Read count']:
                ret_val[line_number][(function, field_name)] = 0.0
            if function in self.tree.data[taxid].attributes:
                ret_val[line_number][(function, '1.Score')] = \
                    self.tree.data[taxid].attributes[function][metric]
                ret_val[line_number][(function, '3.Read count')] = \
                    self.tree.data[taxid].attributes[function]['count']
                if 'identity' in self.tree.data[taxid].attributes[function]:
                    ret_val[line_number][(function, '2.Identity')] = \
                        self.tree.data[taxid].attributes[function]['identity'] \
                        / self.tree.data[taxid].attributes[function]['hit_count']
        line_number += 1
        # If node has children, call convert_node_into_dict recursively
        if self.tree.data[taxid].children:
            for child_id in sorted(self.tree.data[taxid].children):
                child_lines, child_attribute_values = \
                    self.convert_node_into_dict(child_id,
                                                function_list,
                                                line_number,
                                                metric)
                for child_line_number, child_line in child_lines.items():
                    ret_val[child_line_number] = child_line
                line_number += len(child_lines)
                for child_function, child_attrib in child_attribute_values.items():
                    for key, val in child_attrib.items():
                        children_values[child_function][key] += val

            # If read count for at least one function is greater than sum of read
            # counts from all children, some reads map to unidentified
            # taxon. Add a child node for fictional unidentified taxon.
            unidentified_flag = False
            for function in function_list:
                if function in self.tree.data[taxid].attributes and (
                        children_values[function]['count']
                        < self.tree.data[taxid].attributes[function]['count']
                ):
                    unidentified_flag = True
                    break
            # For root node, fictional child name is 'Unclassified'
            # For other node, fictional child name is ' Unclassified <node taxon>'
            # For example, 'Unclassified Proteobacteria'
            if unidentified_flag and self.tree.data[taxid].rank in LOWER_RANKS:
                ret_val[line_number][('', 'Rank')] = LOWER_RANKS[self.tree.data[taxid].rank]
                ret_val[line_number][('', 'Taxon name')] = 'Unclassified ' \
                    + self.tree.data[taxid].name
                if taxid == '1':
                    ret_val[line_number][('', 'Taxon name')] = 'Unclassified'
                # Calculate scores for fictional node
                for function in function_list:
                    for field_name in ['1.Score', '2.Identity', '3.Read count']:
                        ret_val[line_number][(function, field_name)] = 0.0
                    if function in self.tree.data[taxid].attributes and (
                            children_values[function]['count']
                            < self.tree.data[taxid].attributes[function]['count']
                    ):
                        ret_val[line_number][(function, '1.Score')] = \
                            self.tree.data[taxid].attributes[function][metric] \
                            - children_values[function][metric]
                        ret_val[line_number][(function, '3.Read count')] = \
                            self.tree.data[taxid].attributes[function]['count'] \
                            - children_values[function]['count']

                        if 'identity' in self.tree.data[taxid].attributes[function] and (
                                self.tree.data[taxid].attributes[function]['hit_count']
                                > children_values[function]['hit_count']
                        ):
                            ret_val[line_number][(function, '2.Identity')] = (
                                self.tree.data[taxid].attributes[function]['identity']
                                - children_values[function]['identity']
                            ) / (
                                self.tree.data[taxid].attributes[function]['hit_count']
                                - children_values[function]['hit_count']
                            )
                line_number += 1
        return ret_val, attribute_values

    def convert_profile_into_score_df(self, metric='rpkm'):
        """Converts functional-taxonomic profile into pandas DataFrame object

        Args:
            metric (str): score metric (default value 'rpkm') to be reported

        Returns:
            result (pd.DataFrame): functional-taxonomic profile as pandas
                DataFrame object
        """
        function_list = set()
        for taxid in self.tree.data:
            for function in self.tree.data[taxid].attributes.keys():
                function_list.add(function)
        root_id = ROOT_TAXONOMY_ID
        line_number = 1
        tax_dict, _ = self.convert_node_into_values_dict(root_id,
                                                         function_list,
                                                         line_number,
                                                         metric=metric)
        result = pd.DataFrame(tax_dict)
        result = result.transpose()
        return result

    def convert_node_into_values_dict(self, taxid, function_list, line_number, metric='efpkg'):
        """Returns node of functional-taxonomic profile for conversion into DataFrame.
        Recursively called for all children of the node.

        Args:
            taxid (str): taxonomy identifier of node
            function_list (list of str): function identifiers to be included
                to the table
            line_number (int): sequential number of node printed
            metric (str): score metric (default value 'efpkg') to be reported

        Returns:
            attribute_values (defaultdict[str,dict[str,float]]): outer key is
                function identifier, inner key is metric,  value is float.
            ret_val (dict[str,tuple(str,str)]): key is line number, value is a
                tuple with function identifier or empty string as first element
                and field name as second element. Field names are 'Rank',
                'Taxon name', metric.
        """
        # Collect all attributes for reporting to the upper level
        attribute_values = defaultdict(dict)
        for function in function_list:
            if function in self.tree.data[taxid].attributes:
                attribute_values[function][metric] = 0.0
                if metric in self.tree.data[taxid].attributes[function]:
                    attribute_values[function][metric] = \
                        self.tree.data[taxid].attributes[function][metric]

        ret_val = defaultdict(dict)
        children_values = autovivify(2, float)
        if taxid in self.tree.data:
            ret_val[line_number][('', 'Rank')] = self.tree.data[taxid].rank
            ret_val[line_number][('', 'Taxon name')] = self.tree.data[taxid].name
            for function in function_list:
                ret_val[line_number][(function, metric)] = 0.0
                if function in self.tree.data[taxid].attributes:
                    ret_val[line_number][(function, metric)] = \
                        self.tree.data[taxid].attributes[function][metric]
            line_number += 1
            if self.tree.data[taxid].children:
                for child_id in sorted(self.tree.data[taxid].children):
                    child_lines, child_values = \
                        self.convert_node_into_values_dict(child_id,
                                                           function_list,
                                                           line_number,
                                                           metric)
                    for child_line_number, child_line in child_lines.items():
                        ret_val[child_line_number] = child_line
                    line_number += len(child_lines)
                    for datapoint in child_values.keys():
                        for key, val in child_values[datapoint].items():
                            children_values[datapoint][key] += val

                # Add a child node for unidentified child taxon, if needed
                unidentified_flag = False
                for function in function_list:
                    if function in self.tree.data[taxid].attributes and (
                            children_values[function][metric]
                            < self.tree.data[taxid].attributes[function][metric]
                    ):
                        unidentified_flag = True
                        break

                if unidentified_flag and self.tree.data[taxid].rank in LOWER_RANKS:
                    ret_val[line_number][('', 'Rank')] = LOWER_RANKS[self.tree.data[taxid].rank]
                    ret_val[line_number][('', 'Taxon name')] = 'Unclassified ' \
                        + self.tree.data[taxid].name
                    if taxid == '1':
                        ret_val[line_number][('', 'Taxon name')] = 'Unclassified'
                    for function in function_list:
                        ret_val[line_number][(function, metric)] = 0.0
                        if function in self.tree.data[taxid].attributes and (
                                children_values[function][metric]
                                < self.tree.data[taxid].attributes[function][metric]
                        ):
                            ret_val[line_number][(function, metric)] = \
                                self.tree.data[taxid].attributes[function][metric] \
                                - children_values[function][metric]
                    line_number += 1
        else:
            print('Node not found:', taxid)

        return ret_val, attribute_values
