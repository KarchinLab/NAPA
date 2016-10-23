#! /usr/bin/env python 

from ete2 import PhyloTree
from napa.seq.bioseq import *

from napa.utils.serials import * 
from napa.utils.io import *


class PhyloEdge(object):
    ''' 
    PhyloEdge line graph of PhyloTree from ete2.
    PhyloEdge nodes are egdes in the PhyloTree
    Edges are oriented such that parental node always
    comes before child node.
    Sequences on nodes are BioSeq objects.
    '''

    def __init__(self, parent_node = None, child_node = None,
                 aln = None, wt_seq = None,
                 parent = None, children = []):
        
        # parent node in edge, a PhyloTree Object        
        self.parent_node = parent_node 
        self.parent_node_name = "ANC" if not parent_node \
                                else parent_node.name

        # child node in PhyloTree, representing edge target
        self.child_node = child_node
        self.child_node_name = child_node.name

        # name edge by parent-child names
        self.name = '%s__%s' % (self.parent_node_name, 
                                self.child_node_name)

        # assign BioSeq sequence objects
        # to parent and child 
        self.wt_seq = wt_seq
        self.aln = aln
        self.assign_node_seqs()

        # Get mutations for all sequences in alignment
        self.aln.get_seq_muts(wt_seq = self.wt_seq)        
        
        # Get set of mutations/mut positions 
        # along this PhyloTree edge
        self.get_muts()
        self.mut_pos = [get_int_substring(m) \
                        for m in self.muts]

        # parent edge = another node pair 
        # whose target = phyloEdge's source
        self.parent = parent 

        # child edges = node pairs 
        # whose source = phyloEdge's target
        self.get_node_children(children)
                
        self.precs = [] 
        self.followers = []
        self.ordered_precs = []


    def assign_node_seqs(self):
        if self.parent_node_name in self.aln.seqid_to_seq:
            self.parent_node_seq = \
                self.aln.seqid_to_seq[self.parent_node_name]
        else:
            self.parent_node_seq = self.wt_seq

        if self.child_node_name in self.aln.seqid_to_seq:
            self.child_node_seq = \
                self.aln.seqid_to_seq[self.child_node_name] 
        else:
            self.child_node_seq = self.wt_seq
                

    def add_child(self, phylo_edge):
        if phylo_edge.name not in [c.name for c in self.children]:
            self.children.append(phylo_edge)
            phylo_edge.parent = self


    def get_node_children(self, children):
        if len(children):
            self.children = children
            return
        self.children = []
        child_node_children = self.child_node.get_children()
        for child_node_child in child_node_children:
            self.add_child(\
                PhyloEdge(parent_node = self.child_node,
                          child_node = child_node_child,
                          aln = self.aln, wt_seq = self.wt_seq,
                          parent = self))


    def build_edge_tree(self):
        self.children = [child for child in self.children \
                         if child.parent.name == self.name]
        self.build_edge_followers() 
        self.build_edge_precs()


    def build_edge_followers(self):
        for child in self.children:
            child.build_edge_followers()
            self.followers  +=  [child] + child.followers
        self.followers = list(set(self.followers))


    def build_edge_precs(self):
        for child in self.children:
            child.precs = list(set(self.precs + [self]))
            child.build_edge_precs()

   
    def assign_ordered_precs(self):
        '''
        PhyloEdge's precursor edges in order on tree.
        '''
        if not len(self.precs): 
            self.ordered_precs = []
            return

        parent_to_child_node = {}
        for p in self.precs:
            p_name = p.parent_node.name if p.parent_node else 'ANC'
            np_name = p.child_node.name if p.child_node else 'ANC'
            parent_to_child_node[p_name] = np_name

        most_anc_node_name = \
        list(set(parent_to_child_node.keys()) - \
             set(parent_to_child_node.values()))

        if not(len(most_anc_node_name)): 
            stderr_write(['No ancestor found for phylogeny edge:', 
                          self.name])

        if len(most_anc_node_name) > 1: 
            stderr_write(['Found multiple ancestors',
                          'for phylo edge:',
                          self.name])
        
        curr_anc_name = most_anc_node_name[0]
        max_depth = len(set(parent_to_child_node.keys() + \
                            parent_to_child_node.values()))

        for i in range(max_depth - 1):
            next_anc_name = parent_to_child_node[curr_anc_name]
            next_prec = [p for p in self.precs if p.name == \
                             curr_anc_name+'__'+next_anc_name]
            self.ordered_precs.append(next_prec[0])
            curr_anc_name = next_anc_name

        if self.ordered_precs[-1] != self.parent or \
           len(self.ordered_precs)!= len(self.precs):
            print "_precs not properly ordered/found:", self.name
            print "\tparent_to_child_node", parent_to_child_node
            print "\tordered", \
                [p.name for p in self.ordered_precs if p != None]
            print "\tunordered", [p.name for p in self.precs]


    def get_common_precs(self, other, annot_key, 
                         annot_combinations):
        '''
        Common precursors for a two phylo edges. 
        '''
        common_prec = []
        if self.parent == None or other.parent == None:
            common_prec = [] # case of one of the edges being root
        elif self ==  other:
            common_prec = [self, self.parent]
        elif self in other.precs:
            common_prec = [self.parent, self]
        elif other in self.precs:
            common_prec = [other.parent, other]
        else: # they are in different lineages: look for Mrca
            common_prec = list(set(self.precs) & set(other.precs))
            if not len(common_prec):
                print "No common precs for", self.name, other.name
            elif len(common_prec) > 1:
                common_prec = \
                [self.ordered_precs[len(common_prec)-1]]
        
        common_prec = list(set([p for p in common_prec \
                        if p != None and \
                        p.check_node_seq_annot(annot_key, 
                                            annot_combinations)]))
        
        return common_prec



    def get_intermediate_prec(self, prec_list, annot_key, 
                              annot_combinations):
        '''
        Intermediate edges between two phylo edges.
        '''
        intermediates = [self] if self.check_node_seq_annot(\
                            annot_key, annot_combinations) \
                        else []

        if not len(prec_list):
            intermediates  +=  \
            [p for p in self.ordered_precs \
             if p.check_node_seq_annot(annot_key, 
                                       annot_combinations)]
            return intermediates

        first_prec = self.ordered_precs.index(\
            next((p for p in self.ordered_precs if p in prec_list), 
                 self.ordered_precs[-1]))

        return intermediates + \
        [im for im in self.ordered_precs[:first_prec] \
         if im != None and \
         im.check_node_seq_annot(annot_key, annot_combinations)]
        

    def get_muts(self):
        ''' 
        Get mutations connecting two phylogeny nodes (sequences)
        Mutations from a reference sequence are considered
        Overlap is removed.
        '''
        self.muts = []
        # Get mutations from a reference sequence for parent/child 
        # node. Those should have been pre-calculated in 
        # self.aln from external function creating the 
        # PhyloEdge instance
        
        parent_node_mutations, child_node_mutations = [], []
        if self.parent_node_name in self.aln.seqid_to_mut:
            parent_node_mutations = \
                self.aln.seqid_to_mut[self.parent_node_name]
        if self.child_node_name in self.aln.seqid_to_mut:
            child_node_mutations = \
                self.aln.seqid_to_mut[self.child_node_name]
        self.muts = list(set(child_node_mutations) - \
                         set(parent_node_mutations))    


    def is_leaf_edge(self):
        return self.child_node.is_leaf()


    def check_dist(self, other, dist_thresh):
        if dist_thresh <= 0.:
            return True
        elif self.parent_node:
            return \
                (self.parent_node.get_distance(other.child_node) \
                 <=  dist_thresh)
        else:
            return True

    def check_node_seq_annot(self, annot_key, annot_combinations):
        if not len(annot_combinations):
            return True
        if annot_key not in self.parent_node_seq.seq_annot:
            return False
        if annot_key not in self.child_node_seq.seq_annot:
            return False
        
        parent_node_annot = \
                    self.parent_node_seq.seq_annot[annot_key]
        child_node_annot = \
                    self.child_node_seq.seq_annot[annot_key]

        if parent_node_annot in annot_combinations:
            if child_node_annot in \
               annot_combinations[parent_node_annot]:
                return True
        return False


    def check_mut_pos(self):
        if not len(self.mut_pos): return False
        return set(self.mut_pos) <= self.aln.aln_pos

    def __repr__(self):
        to_print = [self.name, 'mutations: ' + ';'.join(self.muts), 
                   'children: ' + \
                    ';'.join([c.name for c in self.children]), 
                   'followers: ' + \
                    ';'.join([f.name for f in self.followers]), 
                   'precs: ' + ';'.join([p.name \
                                         for p in self.precs])]
        if len(self.ordered_precs):
            to_print.append('ordered_prec: ' + ';'.join([op.name \
                        for op in self.ordered_precs if op])) 
        return '\n\t'.join(to_print)
