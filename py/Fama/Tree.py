import queue
RANKS = set(['norank','superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])

class Node(object):
    
    def __init__(self, rank = None, name = None, taxid = None, parent = None, children = []):
        self.parent = parent
        self.children = children
        if not rank or rank in RANKS:
            self.rank = rank
        else:
            raise Exception('Rank not allowed:' + rank)
        self.name = name
        self.taxid = taxid
        self.attributes = {}
    
    def add_child(self,node):
        if node.parent and node.parent == self.parent:
            self.children.append(node)
    
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

    def is_in_children(self, node):
        return any(n for n in self.children if n.taxid == node.taxid) 

class Tree(object):
    
    def __init__(self):
        self.data = {}
        self.root = Node(rank = 'norank',name = 'root', taxid = '1')
        self.data['1'] = self.root
        cellular_organisms = Node(rank = 'norank',name = 'cellular organisms', taxid = '131567', parent = '1')
        self.data['131567'] = cellular_organisms
                
    def add_node(self, node):
        if node.taxid not in self.data and node.parent in self.data:
            self.data[node.taxid] = node
        if not self.data[node.parent].is_in_children(node):
            self.data[node.parent].add_child(node)
    
    def is_in_tree(self,node):
        return node.taxid in self.data
    
    def is_parent_in_tree(self,node):
        return node.parent in self.data
    
    def add_node_recursively(self, node):
        if node.taxid in self.data:
            pass #Nothing to do
        elif node.parent in self.data:
            self.data[node.taxid] = node
        else:
            nodes_stack = queue.LifoQueue()
            current_node = node
            while not is_parent_in_tree(current_node):
                nodes_stack.put(current_node)
                current_node = get_parent(current_node)
            parent_taxid = None
            while not nodes_stack.empty():
                
                current_node = nodes_stack.get()
                if current_node.rank in RANKS:
                    if parent_taxid:
                        current_node.set_parent(parent_taxid)
                    self.add_node(current_node)
                parent_taxid = current_node.taxid
                
    def print_tree(self):
        print(self.data['1'].taxid)
        level = 1
                
def get_parent(node):
    if node.taxid in tax_data.nodes and node.taxid in tax_data.names: 
        parent_taxid = tax_data.nodes[taxid]['parent']
        if parent_taxid in tax_data.nodes and parent_taxid in tax_data.names:
            ret_val = Node(rank = tax_data.nodes[parent_taxid]['rank'],
                           name = tax_data.names[parent_taxid]['name'], 
                           taxid = parent_taxid,
                           parent = tax_data.nodes[parent_taxid]['parent'])
            ret_val.add_child(node)
            return ret_val
    
