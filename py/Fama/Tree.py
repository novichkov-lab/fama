import queue
from collections import defaultdict
RANKS = set(['norank','superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])

class Node(object):
    
    def __init__(self, rank = None, name = None, taxid = None, parent = None, children = set()):
        self.parent = parent
        self.children = children
        if not rank or rank in RANKS:
            self.rank = rank
        else:
            raise Exception('Rank not allowed:' + rank)
        self.name = name
        self.taxid = taxid
        self.attributes = {}
    
    def add_child(self,child_taxid):
        self.children.add(child_taxid)
    
    def set_parent(self,parent_taxid):
        if parent_taxid:
            self.parent = parent_taxid
    
    def set_taxid(self,taxid):
        if not self.taxid:
            self.taxid = taxid
    
    def set_rank(self, rank):
        ret_val = False
        if rank in RANKS:
            self.rank = rank
            ret_val = True
        return ret_val

    def get_children(self):
        return self.children
    
    def set_attribute(self, k, v):
        self.attributes[k] = v

    def add_attribute(self, k, v):
        if isinstance(v,dict) and not self.attributes:
            self.attributes = defaultdict(lambda : defaultdict(dict))
        if k in self.attributes:
            if isinstance(v,dict):
                for k2 in v:
                    if k2 in self.attributes[k]:
                        self.attributes[k][k2] += v[k2]
                    else:
                        self.attributes[k][k2] = v[k2]
            else:
                self.attributes[k] += v
        else:
            if isinstance(v,dict):
                for k2 in v:
                    self.attributes[k][k2] = v[k2]
            else:
                self.attributes[k] = v

    def get_attribute(self, k):
        if k in self.attributes:
            return self.attributes[k]
        else:
            return None
    
    def is_in_children(self, node_id):
        return node_id in self.children 
        
class Tree(object):
    
    def __init__(self):
        self.data = {}
        self.root = Node(rank = 'norank',name = 'root', taxid = '1', children = set())
        self.data['1'] = self.root
        #cellular_organisms = Node(rank = 'norank',name = 'cellular organisms', taxid = '131567', parent = '1')
        #self.data['131567'] = cellular_organisms
        
        #unknown_organisms = Node(rank = 'norank',name = 'Unknown', taxid = '0', parent = '1',children = set())
        #self.add_node(unknown_organisms)
        #self.data['0'] = unknown_organisms
        #self.data['1'].add_child('131567')
        #self.data['1'].add_child('0')
                
    def add_node(self, node):
        if node.taxid not in self.data and node.parent in self.data:
            self.data[node.taxid] = node
            self.data[node.parent].add_child(node.taxid)
        elif not node.parent in self.data:
            print ('Parent ', node.parent, ' for node ', node.taxid, 'not found in the tree')
        elif not self.data[node.parent].is_in_children(node.taxid):
            self.data[node.parent].add_child(node.taxid)
    
    def is_in_tree(self,taxid):
        return taxid in self.data
    
    def add_node_recursively(self, node, tax_data):
        if node.taxid in self.data:
            return #Nothing to do
        elif node.parent in self.data:
            self.data[node.taxid] = node
            self.data[node.parent].add_child(node.taxid)
        else:
            nodes_stack = queue.LifoQueue()
            current_taxid = node.taxid
            nodes_stack.put((current_taxid,True))
            parent_taxid = tax_data.nodes[current_taxid]['parent']
            while True:
                current_taxid = tax_data.nodes[current_taxid]['parent']
                if tax_data.nodes[current_taxid]['rank'] in RANKS:
                    nodes_stack.put((current_taxid,True))
                else:
                    #print(tax_data.nodes[current_taxid]['rank'], ' not in RANKS, skip ID ', current_taxid)
                    nodes_stack.put((current_taxid,False))
                
                if tax_data.nodes[current_taxid]['parent'] == '1':
                    parent_taxid = '1'
                    #print('Root reached in recursion from node ', node.taxid)
                    #print(list(nodes_stack.queue))
                    break
                elif self.is_in_tree(tax_data.nodes[current_taxid]['parent']):
                    parent_taxid = tax_data.nodes[current_taxid]['parent']
                    break
                
            
            
            while not nodes_stack.empty():
                node_id = nodes_stack.get()
                if node_id[1]:
                    node = Node(rank = tax_data.nodes[node_id[0]]['rank'], 
                                name = tax_data.names[node_id[0]]['name'], 
                                taxid = node_id[0], 
                                parent = parent_taxid, 
                                children = set())
                    #print ('\tCall add_node for ', node_id[0], tax_data.nodes[node_id[0]]['rank'], parent_taxid, tax_data.names[node_id[0]]['name'])
                    self.add_node(node)
                    parent_taxid = node_id[0]
                    

    def add_attribute_recursively(self, taxid, k, v):
        parent_taxid = self.data[taxid].parent
        while True:
            self.data[taxid].add_attribute(k,v)
            if taxid == '1':
                break
            taxid = parent_taxid
            parent_taxid = self.data[taxid].parent

    def print_tree(self):
        print(self.data['1'].taxid)
        level = 1
                
    def get_parent(self, node):
        if node.taxid in tax_data.nodes and node.taxid in tax_data.names: 
            parent_taxid = tax_data.nodes[taxid]['parent']
            if parent_taxid in tax_data.nodes and parent_taxid in tax_data.names:
                ret_val = Node(rank = tax_data.nodes[parent_taxid]['rank'],
                               name = tax_data.names[parent_taxid]['name'], 
                               taxid = parent_taxid,
                               parent = tax_data.nodes[parent_taxid]['parent'])
                ret_val.add_child(node)
                return ret_val
    
