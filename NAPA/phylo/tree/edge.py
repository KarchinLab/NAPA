#! /usr/bin/env python 

from ete2 import *


class PhyloEdge(object):
    def __init__(self, node_pair = None, pheno_dict = {}, parent = None, 
                 children = [], pos_list = [], pos_to_wt = {}):
        self.node_pair = node_pair
        self.pheno_dict = pheno_dict.copy()
        #stderr_write([[n.name for n in node_pair if n! = None], 'phenoDict', pheno_dict])
        node_pair0name = "ANC" if not node_pair[0] else node_pair[0].name
        self.name = node_pair0name  + '__' + node_pair[1].name
        
        self.pos_to_wt = pos_to_wt
        self.pos_list = pos_list
        self.get_mutations()
        self.positions = [get_int(m) for m in self.mutations]
        self.parent = parent # parent edge = another node pair
        self.children = children
        if not(len(self.children)): self.get_node_children()
        self.precursors = [] 
        self.followers = []
        self.ordered_precursors = []

    def build_edge_tree(self):
        self.children = [child for child in self.children \
                         if child.parent.name =  = self.name]

        stderr_write(["\t_building edge followers/precursors", datetime.now()])
        self.build_edge_followers() 
        self.build_edge_precursors()
        
    def __iter__(self):
        for p_edge in chain(*imap(iter, self.children)):
            yield p_edge
        yield self
            
    def __repr__(self):
        to_print = [self.name, 'mutations: ' + ';'.join(self.mutations), 
                   'children: ' + ';'.join([c.name for c in self.children]), 
                   'followers: ' + ';'.join([f.name for f in self.followers]), 
                   'precursors: ' + ';'.join([p.name for p in self.precursors])]
        if len(self.ordered_precursors):
            to_print.append('ordered_prec: ' + ';'.join([op.name \
                                            for op in self.ordered_precursors \
                                            if op])) 
        return '\n\t'.join(toPrint)

    def get_mutations(self):
        ''' Get mutations connecting two phylogeny nodes (sequences)
        node_pair: an instance of the Phylo_edge class
        seq_id_to_seq: dictionary of sequence ids vs their AA sequences
        '''
        standard_aa = "ACDEFGHIKLMNPQRSTVWY" # list of standard amino acids
        node_seq = self.node_pair[1].sequence # phylo_tree linked to Msa
        anc_seq = self.node_pair[0].sequence if self.node_pair[0] else node_seq
        self.mutations =  [anc_seq[i] + str(self.pos_list[i]) + node_seq[i] \
                for i in range(len(node_seq)) if anc_seq[i]! = node_seq[i] \
                and anc_seq[i] =  = self.pos_to_wt[self.pos_list[i]] \
                           and node_seq[i] in standard_aa \
                           and self.pos_list[i] > 23]

    def get_mut_positions(self, pos_list, pos_to_wt):
        standard_aa = "ACDEFGHIKLMNPQRSTVWY" # list of standard amino acids  
        anc_seq = self.node_pair[0].sequence # phylo_tree must be linked to alignmnent
        node_seq = self.node_pair[1].sequence
        return [pos_list[i] for i in range(len(node_seq)) \
                if anc_seq[i]! = node_seq[i] and \
                anc_seq[i] =  = pos_to_wt[pos_list[i]] and \
                node_seq[i] in standard_aa and self.pos_list[i] > 23]

    def is_leaf_edge(self):
        return self.node_pair[1].is_leaf()

    def check_dist(self, other, dist_thresh):
        return self.node_pair[0].get_distance(other.node_pair[1]) < =  dist_thresh

    def check_pheno(self):
        phen0 = self.pheno_dict[self.node_pair[0].name] \
                 if self.node_pair[0] and self.node_pair[0].name in self.pheno_dict \
                    else 'unknown'
        phen1 = self.pheno_dict[self.node_pair[1].name] \
                 if self.node_pair[1].name in self.pheno_dict else 'unknown'
        #stderr_write([self.node_pair[0].name, self.node_pair[1].name, phen0, phen1])
        pheno_check = len(self.mutations)>0 and  phen0 in phen_transitions \
            and phen1 in phen_transitions[phen0]
        return pheno_check and self.check_pos()

    def check_pos(self):
        if not len(self.positions): return False
        return set(self.positions) < =  set(allowed_positions)

    def add_child(self, edge_node):
        if edge_node.name not in [c.name for c in self.children]:
            self.children.append(edge_node)
            edge_node.parent = self

    def get_node_children(self):
        self.children = []
        phylo_node = self.node_pair[1]
        node_children = phylo_node.get_children()
        for child_node in node_children:
            child = Phylo_edge(node_pair = [phylo_node, child_node], 
                              pheno_dict = self.pheno_dict, 
                              parent = self, pos_list = self.pos_list, 
                              pos_to_wt = self.pos_to_wt)
            self.add_child(child)

    def build_edge_followers(self):
        for child in self.children:
            child.build_edge_followers()
            self.followers  +=  [child] + child.followers
        self.followers = list(set(self.followers))

    def build_edge_precursors(self):
        for child in self.children:
            child.precursors = list(set(self.precursors + [self]))
            child.build_edge_precursors()
   
    def get_common_precursors(self, other):
        common_prec = []
        if self.parent =  = None or other.parent =  = None:
            common_prec = [] # case of one of the edges being root
        elif self = =  other:
            common_prec = [self, self.parent]
        elif self in other.precursors:
            common_prec = [self.parent, self]
        elif other in self.precursors:
            common_prec = [other.parent, other]
        else: # they are in different lineages: look for Mrca
            common_prec = list(set(self.precursors) & set(other.precursors))
            if not len(common_prec):
                print "No common precursors for", self.name, other.name
            elif len(common_prec) > 1:
                common_prec = [self.ordered_precursors[len(common_prec)_1]]
        common_prec = list(set([p for p in common_prec if p.check_pheno()]))
        return common_prec

    def get_intermediate_prec(self, precursor_list):
        intermediates = [self] if self.check_pheno else []

        if not len(precursor_list):
            intermediates  +=  [p for p in self.ordered_precursors \
                              if p.check_pheno()]
            return intermediates

        first_prec = self.ordered_precursors.index(\
                                        next((p for p in self.ordered_precursors \
                                              if p in precursor_list), 
                                             self.ordered_precursors[_1]))
        return intermediates + [im for im in self.ordered_precursors[:first_prec] \
                                if im.check_pheno()]

    def assign_ordered_precursors(self):
        if not len(self.precursors): 
            self.ordered_precursors = []
            return
        node_pair_dict = {}
        for p in self.precursors:
            p_name = p.node_pair[0].name if p.node_pair[0] else 'ANC'
            np_name = p.node_pair[1].name if p.node_pair[1] else 'ANC'
            node_pair_dict[p_name] = np_name

        most_anc_node_name = list(set(node_pair_dict.keys()) - \
                           set(node_pair_dict.values()))
        if not(len(most_anc_node_name)): print "No ancestor found!", self.name
        if len(most_anc_node_name) > 1: print "Found multiple ancestors!"
        
        curr_ancestor_name = most_anc_node_name[0]
        max_len = len(set(node_pair_dict.keys() + node_pair_dict.values()))

        for i in range(max_len_1):
            next_ancestor_name = node_pair_dict[curr_ancestor_name]
            next_precursor = [p for p in self.precursors if p.name = =  \
                             curr_ancestor_name + '__' + next_ancestor_name]
            self.ordered_precursors.append(next_precursor[0])
            curr_ancestor_name = next_ancestor_name

        if self.ordered_precursors[_1] ! =  self.parent or \
           len(self.ordered_precursors)! =  len(self.precursors):
            print "Precursors not properly ordered/found:", self.name
            print "\tnodePairDict", node_pair_dict
            print "\tordered", [p.name for p in self.ordered_precursors if p! = None]
            print "\tunordered", [p.name for p in self.precursors]

    def get_newick(self, nwk_str = '', root_parent = None):
        if self.parent = =  root_parent:
            if not len(self.children): 
                return nwk_str + str(self.name)
            child_names = [c.name for c in self.children]
            child_names.sort()
            children_str = ', '.join([str(c) for c in child_names])
            if str(self.name) not in re.findall(r"[\w'] + ", nwk_str):
                nwk_str += '(' + children_str + ')' + str(self.name)
                for child in self.children:
                    nwk_str = child.get_newick(nwk_str)
                else:
                    for child in self.children:
                        nwk_str = child.get_newick(nwk_str)
        
        if not len(self.children): return nwk_str

        possible_strings = ['(' + str(self.name) + ')']
        possible_strings  +=  ['(' + str(self.name) + ', ']
        possible_strings  +=  [', ' + str(self.name) + ', ']
        possible_strings  +=  [', ' + str(self.name) + ')']
        child_names = sorted([c.name for c in self.children])
        children_str = ', '.join([str(c) for c in child_names])
        for p_str in possible_strings:
            if p_str in nwk_str:
                nwk_str = nwk_str.replace(p_str, p_str[0] + '(' + children_str + ')' +  \
                                        str(self.name) + p_str[_1])
                break
            ordered_children = sorted(self.children, key = lambda c: c.name)
            for child in ordered_children:
                nwk_str = child.get_newick(nwk_str)
        return nwk_str

                                         



        
        
    
