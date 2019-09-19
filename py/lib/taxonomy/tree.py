"""Describes Node and Tree classes"""
import queue
from collections import defaultdict
from lib.utils.const import RANKS, ROOT_TAXONOMY_ID


class Node(object):
    """Node object stores data of one taxonomic tree node

    Attributes:
    parent (str): parent node identifier
    children (set of str): child node identifiers
    rank (str): taxonomic rank (one of ranks defined in RANKS)
    name (str): text label of the node
    taxid (str): taxonomy identifier (usually, NCBI Taxonomy ID)
    attributes (defaultdict[str,dict[str,obj]) or dict[str,obj]: polymorphic
        structure for storing various scores associated with the node. Can
        be one-level or two-level dictionary.
    """
    def __init__(self, rank=None, name=None, taxid=None, parent=None, children=None):
        """Args:
            rank (str): taxonomic rank (one of ranks defined in RANKS)
            name (str): text label of the node
            taxid (str): taxonomy identifier (usually, NCBI Taxonomy ID)
            parent (str): parent node identifier
            children (set of str): child node identifiers
        """
        self.parent = parent
        self._children = set()
        if children:
            self._children = set(children)
        if rank in RANKS:
            self.rank = rank
        else:
            raise Exception('Rank not allowed:' + rank)
        self.name = name
        self.taxid = taxid
        self.attributes = {}

    def add_child(self, child_taxid):
        """Adds a child node identifier to children attribute

        Args:
            child_taxid (str): child node identifier
        """
        self._children.add(child_taxid)

    def set_parent(self, parent_taxid):
        """Sets parent identifier if not None

        Args:
            parent_taxid (str): parent node identifier
        """
        if parent_taxid:
            self.parent = parent_taxid

    def set_taxid(self, taxid):
        """Sets taxonomy identifier if not None

        Args:
            taxid (str): taxonomy identifier
        """
        if not self.taxid:
            self.taxid = taxid

    def set_rank(self, rank):
        """Sets taxonomic rank identifier if it defined in RANKS

        Args:
            rank (str): taxonomic rank

        Returns:
            ret_val (bool): True on success, False on failure
        """
        ret_val = False
        if rank in RANKS:
            self.rank = rank
            ret_val = True
        return ret_val

    @property
    def children(self):
        """Set of children identifiers"""
        return self._children

    def set_attribute(self, key, value):
        """Adds a key-value pair to self.attributes. If the key already
        exists, replaces existing value with the new

        Args:
            key (str): key for node attributes
            value (obj): value for node attributes
        """
        self.attributes[key] = value

    def add_attribute(self, key, value):
        """Adds a key-value pair to self.attributes. If the key already
        exists, adds new value to existing one.

        Args:
            key (str): key for node attributes
            value (obj): value for node attributes. May be dict of arbitrary objects
        """
        if isinstance(value, dict) and not self.attributes:
            self.attributes = defaultdict(dict)
            for innner_key in value:
                if isinstance(value[innner_key], dict):
                    print(key, value[innner_key])
                    raise TypeError('Second-level attribute value cannot be dict')
                self.attributes[key][innner_key] = value[innner_key]
        elif key in self.attributes:
            if isinstance(value, dict):
                for innner_key in value:
                    if isinstance(value[innner_key], dict):
                        print(key, value[innner_key])
                        raise TypeError('Second-level attribute value cannot be dict')
                    if innner_key in self.attributes[key]:
                        self.attributes[key][innner_key] += value[innner_key]
                    else:
                        self.attributes[key][innner_key] = value[innner_key]
            else:
                self.attributes[key] += value
        else:
            if isinstance(value, dict):
                for innner_key in value:
                    if isinstance(value[innner_key], dict):
                        print(key, value[innner_key])
                        raise TypeError('Second-level attribute value cannot be dict')
                    self.attributes[key][innner_key] = value[innner_key]
            else:
                self.attributes[key] = value

    def get_attribute(self, key):
        """Returns attribute value for the given key. If key does not exist,
        returns None.

        """
        result = None
        if key in self.attributes:
            result = self.attributes[key]
        return result

    def is_in_children(self, node_id):
        """Returns True if node_id is in children, False if not"""
        return node_id in self.children


class Tree(object):
    """Tree stores phylogenetic tree. Phylogenetic tree is an oriented graph
    of nodes. Tree has root node with identifier equal to ROOT_TAXONOMY_ID.

    Attributes:
        data (dict[str,:obj:Node]): key is taxonomy identifier, value is Node object
        root (:obj:Node): root node
    """
    def __init__(self):
        self.data = {}
        self.root = Node(rank='norank', name='root', taxid=ROOT_TAXONOMY_ID, children=set())
        self.data[ROOT_TAXONOMY_ID] = self.root

    def add_node(self, node):
        """Adds node to phylogenetic tree, if it does not exist already.

        Args:
            node (:obj:Node): new node
        """
        if node.taxid not in self.data and node.parent in self.data:
            self.data[node.taxid] = node
            self.data[node.parent].add_child(node.taxid)
        elif node.parent not in self.data:
            print('Parent ', node.parent, ' for node ', node.taxid, 'not found in the tree')
        elif not self.data[node.parent].is_in_children(node.taxid):
            self.data[node.parent].add_child(node.taxid)

    def is_in_tree(self, taxid):
        """Checks if taxonomy ID exists in tree. Returns True if yes, False if not

        Args:
            taxid (str): taxonomy identifier
        """
        return taxid in self.data

    def add_node_recursively(self, node, taxonomy_data):
        """Adds node to phylogenetic tree, if it does not exist already.
        If parent nodes do not exists, creates them up to root level.
        Updates attributes of all upper nodes (up to root), adding values
        of all attributes of the new node to them.

        Args:
            node (:obj:Node): new node
            taxonomy_data (:obj:TaxonomyData): taxonomic data instance
        """
        if node.taxid in self.data:
            return  # Nothing to do
        elif node.parent in self.data:
            self.data[node.taxid] = node
            self.data[node.parent].add_child(node.taxid)
        else:
            nodes_stack = queue.LifoQueue()
            current_taxid = node.taxid
            nodes_stack.put((current_taxid, True))
            parent_taxid = '0'
            if taxonomy_data.is_exist(current_taxid):
                parent_taxid = taxonomy_data.get_parent(current_taxid)
            while True:
                if taxonomy_data.is_exist(current_taxid):
                    current_taxid = taxonomy_data.get_parent(current_taxid)
                else:
                    current_taxid = '0'

                if taxonomy_data.get_rank(current_taxid) in RANKS:
                    nodes_stack.put((current_taxid, True))
                else:
                    nodes_stack.put((current_taxid, False))

                if taxonomy_data.get_parent(current_taxid) == ROOT_TAXONOMY_ID:
                    parent_taxid = ROOT_TAXONOMY_ID
                    break
                elif self.is_in_tree(taxonomy_data.get_parent(current_taxid)):
                    parent_taxid = taxonomy_data.get_parent(current_taxid)
                    break

            while not nodes_stack.empty():
                node_id = nodes_stack.get()
                if node_id[1]:
                    if taxonomy_data.is_exist(node_id[0]):
                        node = Node(rank=taxonomy_data.get_rank(node_id[0]),
                                    name=taxonomy_data.get_name(node_id[0]),
                                    taxid=node_id[0],
                                    parent=parent_taxid,
                                    children=set())
                    else:
                        node = Node(rank='norank',
                                    name='taxId ' + node_id[0],
                                    taxid=node_id[0],
                                    parent=parent_taxid,
                                    children=set())

                    self.add_node(node)
                    parent_taxid = node_id[0]

    def add_attribute(self, taxid, key, value, taxonomy_data):
        """Adds value for the given key to a node with 'taxid' identifier. If
        such identifier does not exists on the tree, looks for parent identifier and
        tries to add value to it, and so on (up to root).

        Args:
            taxid (str): taxonomy identifier
            key (str): attributes key
            value (obj): value to be added to attributes. Expected to be int or float
            tax_data (:obj:TaxonomyData): taxonomic data instance
        """
        while True:
            if taxid in self.data or taxid == '1':
                break
            taxid = taxonomy_data.get_parent(taxid)
        self.data[taxid].add_attribute(key, value)

    def add_attribute_recursively(self, taxid, key, value, taxonomy_data):
        """Adds value for the given key to a node with 'taxid' identifier and
        ALL its parent nodes, up to root node. If
        such identifier does not exists on the tree, looks for parent identifier and
        tries to add value to it, and so on (up to root).

        Args:
            taxid (str): taxonomy identifier
            key (str): attributes key
            value (obj): value to be added to attributes. Expected to be int or float
            taxonomy_data (:obj:TaxonomyData): taxonomic data instance
        """
        while True:
            if taxid in self.data or taxid == '1':
                break
            taxid = taxonomy_data.get_parent(taxid)

        parent_taxid = self.data[taxid].parent
        while True:
            self.data[taxid].add_attribute(key, value)
            if taxid == '1':
                break
            taxid = parent_taxid
            parent_taxid = self.data[taxid].parent

    def get_parent(self, node, taxonomy_data):
        """Returns parent node for given node.

        Args:
            node (:obj:Node): node of taxonomic tree
            taxonomy_data (:obj:TaxonomyData): taxonomic data instance
        """
        ret_val = None
        if node.taxid in self.data:
            ret_val = self.data[node.parent]
        elif taxonomy_data.is_exist(node.taxid):
            parent_taxid = taxonomy_data.get_parent(node.taxid)
            if taxonomy_data.is_exist(parent_taxid):
                ret_val = Node(rank=taxonomy_data.get_rank(parent_taxid),
                               name=taxonomy_data.get_name(parent_taxid),
                               taxid=parent_taxid,
                               parent=taxonomy_data.get_parent(parent_taxid))
                ret_val.add_child(node)
        return ret_val
