#! /usr/bin/env python 

# Notes: consider removing ete2/3 dependence
# ete is GPL + only need the PhyloTree class implementation
#
# https://github.com/etetoolkit/ete/blob/master/\
# ete3/phylo/phylotree.py
#
# need to be able to read Newick (*.nwk) files and build a 
# general (nonbinary) phylo tree with attributes on
# nodes and edges 
# FORMAT  DESCRIPTION
# 1       flexible with internal node names

import copy
from scipy.stats import fisher_exact
from ete2 import PhyloTree

from napa.utils.serials import * 
from napa.utils.io import *

from napa.phylo.edge import *
from napa.phylo.tree import get_mrca_subtree

from joblib import Parallel, delayed
import multiprocessing
import dill

num_cores = multiprocessing.cpu_count()

# Mutation pairs are not allowed to share
# same numerical position
list_mut_pairs = lambda mut_list: \
                 [pair for pair in list_pairs(mut_list) \
                  if get_int_substring(str(pair[0])) != \
                  get_int_substring(str(pair[1]))]

#==============================================================#
class TreeMut(object):
    '''
    Similar to AlnMut class, but counts phylo edges instead of 
    alignment sequences of occurrence.
    '''
    def __init__(self, phylo_edge_pairs = set(), mut_str = None):
        self.mut_str = mut_str
        self.phylo_edge_pairs = set(phylo_edge_pairs)
    
    def add_phylo_edge_pair(self,phylo_edge_pair):
        self.phylo_edge_pairs.add(phylo_edge_pair) 
    
    def add_phylo_edge_pairs(self,phylo_edge_pairs):
        self.phylo_edge_pairs.update(phylo_edge_pairs)

    def __repr__(self):
        return self.mut_str


#==============================================================#
class TreeMutPair(object):
    '''
    Pair of mutations from a single tree.
    '''    
    def __init__(self, phylo_edge_pairs = set(), 
                 mut_pair = None, wt_seq = None,
                 all_edges = set(), 
                 prot_func_transitions = {'':''}):

        # List of pairs of phylo_edge_objects
        if not len(phylo_edge_pairs):
            raise ValueError(\
                ' '.join(['In TreeMutPair:',
                          'Set of Phylo edge pairs', 
                          'must be nonempty!']))
        else:
            self.phylo_edge_pairs = set(phylo_edge_pairs)

        # A unique pair of AA muts
        if mut_pair ==  None or not len(mut_pair):
            raise ValueError('In TreeMutPair: No mutations given!')
        else:
            self.mut_pair = mut_pair
            self.name = '%s_%s' % (mut_pair[0].mut_str, 
                                   mut_pair[1].mut_str)
        if not len(all_edges):
            raise ValueError(' '.join(['In TreeMutPair:',
                                       'No edge pairs having each',
                                       'mutation were given!']))
        else:
            self.all_edges = set(all_edges)
        self.prot_func_transitions = prot_func_transitions

        if wt_seq == None:
            raise ValueError('Wild type sequence not provided ' + \
                             str(self))
        else:
            self.wt_seq = wt_seq

    def __repr__(self):
        return '%s\t%s'%(str(self.mut_pair[0]), 
                         str(self.mut_pair[1])) 

    def add_phylo_edge_pair(self, phylo_edge_pair):
        self.phylo_edge_pairs.add(phylo_edge_pair) 
    
    def add_phylo_edge_pairs(self, phylo_edge_pairs):
        self.phylo_edge_pairs.update(phylo_edge_pairs)


    def assign_lineage_distribution(self):
        common_lineage, mut0_lineage, mut1_lineage = [], [], []
        for pep in self.phylo_edge_pairs:
            pe0, pe1 = pep[0], pep[1]                           
            if pe0 == pe1:
                common_lineage  +=  [pe0.parent] + [pe0] + \
                             pe0.followers 
            else:
                common_prec = pe0.get_common_precs(pe1, 
                            'function', self.prot_func_transitions)
                common_follow = list(set([pe0] + pe0.followers) & \
                                    set([pe1] + pe1.followers))
                common_lineage  +=  common_prec + common_follow

                mut0_prec = pe0.get_intermediate_prec(common_prec, 
                            'function', self.prot_func_transitions)
                mut1_prec = pe1.get_intermediate_prec(common_prec, 
                            'function', self.prot_func_transitions)

                mut0_lineage += \
                    list(set(mut0_prec + pe0.followers) \
                         - set(common_prec + common_follow))
                mut1_lineage += \
                    list(set(mut1_prec + pe1.followers) \
                        - set(common_prec + common_follow))
        
        self.common_lineage = set([p for p in common_lineage \
                                   if p and p.muts])  
        self.mut0_lineage = set([p for p in mut0_lineage \
                            if p and p.muts]) - self.common_lineage
        self.mut1_lineage = set([p for p in mut1_lineage \
                            if p and p.muts]) - self.common_lineage
        assoc_lineage = self.common_lineage 
        assoc_lineage |= self.mut0_lineage | self.mut1_lineage
        self.other_lineage = self.all_edges - assoc_lineage
        


    def print_lineages(self):
        stderr_write(["\tCommon Lineage:", 
                      [p.name for p in self.common_lineage]])
        stderr_write(["\tFirst mutation only lineage:", 
                      [p.name for p in self.mut0_lineage]])
        stderr_write(["\tSecond mutations only lineage:", 
                     [p.name for p in self.mut1_lineage]])
        stderr_write(["\tOther mutations lineage:", 
                     [p.name for p in self.other_lineage]])


    def get_contingency_table(self):
        '''
        Contingency table for Fisher's Exact Test.
        Uses overall size of clades where mutations in pair
        occurred. 
        '''
        self.assign_lineage_distribution()
        self.num11 = len(self.common_lineage) # both mutations
        self.num10 = len(self.mut0_lineage) # just the 1st not 2nd
        self.num01 = len(self.mut1_lineage) # just the 2nd not 1st
        self.num00 = len(self.other_lineage) # neither mutation
        self.contingency_table = [[self.num11, self.num01], 
                                  [self.num10, self.num00]]


    def get_contingency_table_edge_counts(self):
        ''' 
        A simpler contingency table based on counts of 
        co-occurrence of mutations in tree edges 
        rather than clades. 
        Looking at clade size tends to be more informative.
        '''
        self.num11 = len(self.phylo_edge_pairs)
        self.num10 = len(self.mut_pair[0].phylo_edge_pairs - \
                         self.phylo_edge_pairs) 
        self.num01 = len(self.mut_pair[1].phylo_edge_pairs - \
                         self.phylo_edge_pairs)
        self.num00 = len(self.all_phylo_edge_pairs\
                          -self.mut_pair[0].phylo_edge_pairs \
                          -self.mut_pair[1].phylo_edge_pairs)
        self.contingency_table = [[self.num11, self.num01], 
                                  [self.num10, self.num00]]


    def get_mod_jaccard(self):
        '''
        Modified Jaccard score for clade / edge count overlap
        '''
        self.get_contingency_table()

        [[num11, num10], [num01, num00]] = self.contingency_table
        # Penalizes for small co-ocurrence num11
        epsilon = float(num00)/(num11 + num10 + num01)

        self.mod_jaccard = max(0., 
            (self.num11 - epsilon) / float(num11 + num01 + num10))


    def get_fisher_pval(self):
        '''
        Uses scipy.stats to calculate Fisher's 
        exact p-values based 
        on the contingency table. For association, only need the 
        right-tailed p-value: pval_more.
        '''
        self.assign_lineage_distribution()
        self.get_contingency_table()
        self.odds_ratio_less, self.p_val_less = \
            fisher_exact(self.contingency_table, 
                         alternative = 'less')

        self.odds_ratio_more, self.p_val_more = \
            fisher_exact(self.contingency_table, 
                         alternative = 'greater')

#==============================================================#
class TreeMutPairSet(object):

    def __init__(self, aln = None, aln_pos = [], pos_subset = [],
                 wt_seq = None, leaf_fasta_file = '',
                 int_node_fasta_file = '',
                 leaf_seqid_to_prot_func = {'':''}, 
                 int_seqid_to_prot_func = {'':''}, 
                 prot_func_transitions = {}, 
                 sel_prot_func_list = [''], 
                 dist_thresh = -1., 
                 mut_pair_type = 'dir', method = 'mod_jaccard',
                 tree_nwk_file = ''):

        stderr_write(['\nMutation set for tree:\n',
                      tree_nwk_file])

        # Wild-type/ancestral sequence
        if wt_seq == None:
            raise ValueError('Wild-type sequence not provided ' + \
                             'for tree mut pair set based on:\n' \
                             + tree_nwk_file)
        else:
            self.wt_seq = wt_seq

        self.tree_file = tree_nwk_file
        # Add reconstructed internal node sequences to leaf aln
        self.aln = aln.copy()
        self.aln.annot_seqs(annot_key = 'function',
                            seqid_to_annot = \
                            leaf_seqid_to_prot_func)

        self.aln.add_sequences_from_file(int_node_fasta_file)

        # Option to annotate sequences by (reconstructed) function
        # Add reconstructed internal node sequences functions
        self.aln.annot_seqs(annot_key = 'function',
                            seqid_to_annot = \
                            int_seqid_to_prot_func)

        # Names/ids of all tree nodes having selected function
        self.get_sel_func_nodes(annot_key = 'function',
                                annot_list = sel_prot_func_list)


        # Allowed protein function transitions along branches
        self.prot_func_transitions = prot_func_transitions


        # Option to focus on protein region(s)/domain(s) 
        # of interest
        #self.aln.subset_aln_pos(pos_subset = pos_subset)



        # Get the sequence ids of all leaf sequences annotated
        # with the selected function and find most recent 
        # common ancestor
        # This isolates subtree with given function
        self.anc_node, self.sel_func_leaves = \
                    get_mrca_subtree(tree_nwk_file = \
                                     tree_nwk_file,
                                     fasta_file = leaf_fasta_file,
                                     sel_annot_nodes = \
                                     self.sel_func_node_ids)

        # Get muts for seqs corresponding to protein function
        self.get_leaf_func_muts()


        # Set of all phylogenetic tree edges 
        # When looking at shared edges not lineages
        #self.all_phylo_edge_pairs = set()

        # The folling generates a linegraph object of the 
        # whole tree
        # It is a tree rooted in the edge between anc_node 
        # and its ancestor
        stderr_write(['Building Phylo Edge Tree...'])
        self.anc_edge = PhyloEdge(parent_node = self.anc_node.up, 
                                  child_node = self.anc_node,
                                  aln = self.aln, 
                                  wt_seq = self.wt_seq)
         # Nodes of PhyloEdge are edges are in the original tree
        self.anc_edge.build_edge_tree()

        # Get set of all pairs of mutations occurring along tree
        # default directed mutation pairs
        stderr_write(['Obtaining mutation pairs along tree...'])
        self.mut_pair_type = 'undir' \
                             if 'undir' in mut_pair_type.lower() \
                             else 'dir'
        self.get_mut_pairs(dist_thresh, self.mut_pair_type)

        for mp in self.mut_pair_to_obj:
            mut_pair = self.mut_pair_to_obj[mp]
            if 'fisher' in method.lower():
                mut_pair.pval_more = 1. # default setting
                mut_pair.get_fisher_pval()
            else:
                mut_pair.mod_jaccard = 0.
                mut_pair.get_mod_jaccard()


    def __repr__(self):
        out_str = ''
        for mp in self.mut_pair_to_weight:
            out_str += '%s\t%.6f\n'%('\t'.join(mut_pair), 
                            self.mut_pair_to_weight[mut_pair])
        return out_str


    def get_sel_func_nodes(self, annot_key = 'function',
                           annot_list = ['']):
        if annot_list == [''] and len(annot_list_file):
            annot_list = parse_column(annot_list_file)
        
        self.sel_func_node_ids = self.aln.seqids_with_annot(\
                                    annot_key = annot_key,
                                    annot_list = annot_list)
        
    def get_leaf_func_muts(self):
        '''
        Get muts associated with a function that occur 
        in selected tree leaves. 
        This is the subset of all mutations on the 
        tree that is most associated with the function of 
        interest.
        '''
        self.func_muts = []
        for seqid in [leaf.name for leaf in self.sel_func_leaves]:
            self.func_muts += self.wt_seq.get_substitutions(\
                                self.aln.seqid_to_seq[seqid])
        self.func_muts = set(self.func_muts)


    def check_edge_pair(self,start_edge, follow_edge, 
                        dist_thresh):
        if start_edge.is_leaf_edge(): 
            return False
        if not self.check_edge(start_edge, dist_thresh):
            return False
        if not self.check_edge(follow_edge, dist_thresh):
            return False
        if not start_edge.check_dist(follow_edge, dist_thresh):
            return False
        return True

    def check_edge(self, edge, dist_thresh):
        if len(self.prot_func_transitions) and \
           not edge.check_node_seq_annot('function', 
                                self.prot_func_transitions):
            return False

        if len(self.func_muts) and \
           not len(set(edge.muts) & self.func_muts):
            return False

        # Edge is too long -- exceeds distance thresholds
        if not edge.check_dist(edge, dist_thresh):
            return False
            
        return True
        

    def add_phylo_edge_pair(self, phylo_edge_pair, mut_list_pair):
        '''
        Add a phylogeny edge pair on which edges each of the 
        two mutations occurs.
        mut_list_pair: mutation lists for first and second edge
        in phylo_edge_pair
        Also keep track of individual mutation occurrence.
        '''
        mut_set = set(mut_list_pair[0]+ mut_list_pair[1])
        for mut_str in mut_set:
            if mut_str not in self.mut_str_to_obj:
                self.mut_str_to_obj[mut_str] = \
                    TreeMut(phylo_edge_pairs = \
                            set(phylo_edge_pair),
                            mut_str = mut_str)
            else:
                self.mut_str_to_obj[mut_str].\
                    add_phylo_edge_pair(phylo_edge_pair)

        for mp in zip(mut_list_pair[0], mut_list_pair[1]):
            if get_int_substring(mp[0]) == get_int_substring(mp[1]):
                continue
            if mp not in self.mut_pair_to_obj:
                mut_pair = \
                TreeMutPair(\
                wt_seq = self.wt_seq,
                phylo_edge_pairs = set([phylo_edge_pair]),
                mut_pair = [self.mut_str_to_obj[mp[0]],
                            self.mut_str_to_obj[mp[1]]],
                all_edges = self.all_edges)
                self.mut_pair_to_obj[mp] = mut_pair
            else:
                mut_pair = self.mut_pair_to_obj[mp]
                mut_pair.add_phylo_edge_pair(phylo_edge_pair)

    def get_mut_pairs(self, dist_thresh, mut_pair_type):
        ''' 
        Extract pairs of mutations occurring along tree edges
        '''
        self.all_edges = [self.anc_edge] + \
                         [fe for fe in self.anc_edge.followers \
                          if fe]

        stderr_write(['\n' + self.tree_file, '\n',
            'Number of tree edges considered:', 
            len(self.all_edges)])

        self.all_edges = [e for e in self.all_edges \
                          if self.check_edge(e, dist_thresh)]

        stderr_write(['\n' + self.tree_file, '\n',
                      'Number of tree edges passing',
                      'distance threshold:', 
                      len(self.all_edges)])

        self.mut_pair_to_obj, self.mut_str_to_obj = {}, {}

        for start_edge in self.all_edges:
            start_edge_sorted_muts = \
                sort_by_digits(list(set(start_edge.muts) & \
                                     self.func_muts))
            if not len(start_edge_sorted_muts): continue
            
            # Make sure mutation pairs occurring on 
            # same edge included in both directed and 
            # undirected networks
            if len(start_edge_sorted_muts) > 1:
                self.add_phylo_edge_pair(\
                    phylo_edge_pair = (start_edge, start_edge), 
                    mut_list_pair = \
                    zip(*list_mut_pairs(start_edge_sorted_muts)))
            
            if 'undir' in mut_pair_type.lower():
                self.add_undir_pairs(start_edge, 
                                     start_edge_sorted_muts, 
                                     dist_thresh)
            else:
                self.add_dir_pairs(start_edge, 
                                   start_edge_sorted_muts, 
                                   dist_thresh)
        
        stderr_write(['\n' + self.tree_file, '\n',       
                      'Number of mutations associated with',
                      'function:', 
                      len(self.mut_str_to_obj.keys())])
        stderr_write(['\n' + self.tree_file, '\n',  
                      'Number of mutation pairs extracted', 
                      len(self.mut_pair_to_obj.keys())])


    def add_undir_pairs(self, start_edge, start_edge_sorted_muts, 
                        dist_thresh):
        for follow_edge in start_edge.followers:
            follow_edge_sorted_muts = \
                sorted(list((set(follow_edge.muts) & \
                             self.func_muts) - \
                            set(start_edge.muts)))

            if self.check_edge_pair(start_edge, follow_edge, 
                                    dist_thresh):
                # For undirected mutation pairs, 
                # mutations are ordered by 
                # residue position
                mut_pairs = [(sm, fm) \
                    for sm in start_edge_sorted_muts \
                    for fm in follow_edge_sorted_muts \
                    if get_int_substring(sm) < get_int_substring(fm)]

                if len(mut_pairs):
                    self.add_phylo_edge_pair(\
                        phylo_edge_pair = (start_edge, follow_edge), 
                        mut_list_pair = zip(*mut_pairs))


            # Reverse order of mutations in cases where follower
            # mutation position is less than start position and distance 
            # from follower to start is less than threshold
            if self.check_edge_pair(follow_edge, start_edge, 
                                    dist_thresh):

                rev_mut_pairs = [(fm, sm) \
                    for sm in start_edge_sorted_muts \
                    for fm in follow_edge_sorted_muts \
                    if get_int_substring(sm) > get_int_substring(fm)]
                
                if len(rev_mut_pairs):
                    self.add_phylo_edge_pair(\
                        phylo_edge_pair = (follow_edge, start_edge), 
                        mut_list_pair = zip(*rev_mut_pairs))


    def add_dir_pairs(self, start_edge, start_edge_sorted_muts, 
                      dist_thresh):
        for follow_edge in start_edge.followers:
            follow_edge_sorted_muts = \
                sorted(list((set(follow_edge.muts) & \
                             self.func_muts) - \
                            set(start_edge.muts)))

            if self.check_edge_pair(start_edge, follow_edge, 
                                    dist_thresh):
                # For undirected mutation pairs, 
                # mutations are ordered by 
                # residue position
                mut_pairs = [(sm, fm) \
                             for sm in start_edge_sorted_muts \
                             for fm in follow_edge_sorted_muts]
                self.add_phylo_edge_pair(\
                    phylo_edge_pair = (start_edge, follow_edge), 
                    mut_list_pair = [start_edge_sorted_muts,
                                     follow_edge_sorted_muts])

#==============================================================#
class TreeEnsembleMutPairs(object):

    def __init__(self, inp):
        self.__dict__.update(vars(inp))
    
        self.num_trees = len(flatten(self.tree_files.values()))
        self.mut_pair_to_weight = defaultdict(float)
        self.mut_pair_to_trees = defaultdict(int)

        self.prefix_i = \
        [(prefix, i) for prefix in self.tree_files \
         for i in range(len(self.tree_files[prefix]))]
        

    def __repr__(self):
        out_str = ''
        for mp in self.mut_pair_to_weight: 

            if self.tree_support_thresh > 0:
                if self.mut_pair_to_trees[mp] < \
                   self.tree_support_thresh:
                    continue

            if 'jaccard' in self.method.lower():
                if self.mut_pair_to_weight[mp] < self.thresh:
                    continue
            if 'fisher' in self.method.lower():
                if self.mut_pair_to_weight[mp] > self.thresh:
                    continue

            out_str += '%s\t%.6f\n' % ('\t'.join(mp), 
                self.mut_pair_to_weight[mp])

        return out_str

    def assign_mut_pair_weights(self, mut_pair_dicts):
        for dict_pair in mut_pair_dicts:
            mut_pair_to_weight = dict_pair[0]
            mut_pair_to_trees = dict_pair[1]

            self.mut_pair_to_weight = \
            num_dict_update_add(self.mut_pair_to_weight,
                                mut_pair_to_weight)
            
            self.mut_pair_to_trees = \
            num_dict_update_add(self.mut_pair_to_trees,
                                mut_pair_to_trees)
       
    def get_tree_mut_pair_set(self, prefix, i):
        mut_pair_to_weight = defaultdict(float)
        mut_pair_to_trees = defaultdict(int)
        tree_file = self.tree_files[prefix][i]

        int_seq_file = self.int_seq_files[prefix][i]

        int_func_file = self.int_func_files[prefix][i] \
                        if len(self.int_func_files) \
                        else ''

        int_seqid_to_prot_func = \
             parse_keyval_dict(int_func_file)

        mut_pair_set = \
            TreeMutPairSet(\
            aln = self.aln, aln_pos = self.pos_list, 
            pos_subset = self.pos_subset,
            wt_seq = self.wt_seq, 
            leaf_fasta_file = self.aln_fasta_file,
            int_node_fasta_file = int_seq_file,
            leaf_seqid_to_prot_func = \
            self.seqid_to_prot_func, 
            int_seqid_to_prot_func = int_seqid_to_prot_func,
            prot_func_transitions = self.func_transitions,
            sel_prot_func_list = self.sel_prot_func,
            dist_thresh = self.dist_thresh,
            mut_pair_type = self.edge_type, 
            method = self.method, 
            tree_nwk_file = tree_file)
        
        for mp in mut_pair_set.mut_pair_to_obj:
            mut_pair = mut_pair_set.mut_pair_to_obj[mp]

            if 'fisher' in self.method.lower():
                
                if mut_pair.p_val_more <= self.thresh:
                    self.mut_pair_to_weight[mp] += \
                                    1./self.num_trees
                    self.mut_pair_to_trees[mp] += 1
                
            else:
                if mut_pair.mod_jaccard > 0.0:
                    self.mut_pair_to_weight[mp] += \
                    mut_pair.mod_jaccard / float(self.num_trees)
                    self.mut_pair_to_trees[mp] += 1
        

    def add_mut_pairs(self, mut_pair_set):
        if 'fisher' in self.method.lower():
            self.add_mut_pairs_fisher(mut_pair_set)
        else:
            self.add_mut_pairs_mod_jaccard(mut_pair_set)


    def add_mut_pairs_fisher(self, mut_pair_set):
        '''
        to build network, add mutation pairs 
        from a set calculated for a single tree
        all mutations passing the Fisher's exact test 
        with a significant right-tailed p-value are added
        '''
        for mp in mut_pair_set.mut_pair_to_obj:
            mut_pair = \
            mut_pair_set.mut_pair_to_obj[mp]

            if mut_pair.p_val_more <= self.thresh:
                self.mut_pair_to_weight[mp] += \
                        1./self.num_trees
                self.mut_pair_to_trees[mp] += 1
                


    def add_mut_pairs_mod_jaccard(self, mut_pair_set):
        '''
        to build a network, will take the sum of modified 
        Jaccard index and add mutation pairs from a set 
        calculated for a single tree all the mutations 
        with non-zero Jaccard indices.
        The weight will be averaged over tree topologies
        '''
        
        for mp in mut_pair_set.mut_pair_to_obj:
            mut_pair = \
                mut_pair_set.mut_pair_to_obj[mp]

            self.mut_pair_to_weight[mp] += \
                mut_pair.mod_jaccard / float(self.num_trees)
            self.mut_pair_to_trees[mp] += 1
        
    def write_network_to_file(self, file_path):
        '''
        Writes network to file with path file_path
        Columns: source_node  target_node  weight
        '''
        with open(file_path, 'wb') as f:
            f.write(str(self))      
        
               
    def write_table_to_file(self, file_path):
        '''
        Prints a table of all co-occurring mutations.
        Columns: Each mutation in pair,
        weight from phylogeny, 
        number of trees in ensemble passing
        weight threshold.
        '''
        with open(file_path, 'wb') as f:
            f.write('\t'.join(\
                ['Mut1', 'Mut2', 'Weight_' + str(self.method),
                 'Tree_Count']) + '\n')
                    
            for mp in self.mut_pair_to_weight:
                weight = self.mut_pair_to_weight[mp]
                num_trees = self.mut_pair_to_trees[mp]
                f.write('%s\t%.6f\t%d\n' % ('\t'.join(mp), 
                                            weight, num_trees))


    def print_stats(self, aln_pos = [], pos_subset = []):

        stderr_write(['\nPHYLOGENY-BASED NETWORK PROPERTIES:'])

        length = max(len(aln_pos), self.aln.length)

        stderr_write(['Alignment depth (num. rows)',
                      str(self.aln.depth) + \
                      '\nAlignment length (num. columns):', 
                      length])
        
        stderr_write([len(self.aln.aln_pos), 
                      'positions considered for network.'])

        stderr_write(['Number of ensemble trees considered:',
                      self.num_trees])

        stderr_write([len(self.mut_pair_to_weight), 
                    'Mutation pairs with non-zero',
                    'edge weight.\n'])


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def get_tree_mut_pair_dict(tree_ensemble, prefix, i):
       
    mut_pair_to_weight = defaultdict(float)
    mut_pair_to_trees = defaultdict(int)

    te = tree_ensemble
    tree_file = te.tree_files[prefix][i]

    int_seq_file = te.int_seq_files[prefix][i]

    int_func_file = te.int_func_files[prefix][i] \
                    if len(te.int_func_files) \
                    else ''

    int_seqid_to_prot_func = \
         parse_keyval_dict(int_func_file)

    mut_pair_set = \
        TreeMutPairSet(\
        aln = te.aln, aln_pos = te.pos_list, 
        pos_subset = te.pos_subset,
        wt_seq = te.wt_seq, 
        leaf_fasta_file = te.aln_fasta_file,
        int_node_fasta_file = int_seq_file,
        leaf_seqid_to_prot_func = \
        te.seqid_to_prot_func, 
        int_seqid_to_prot_func = int_seqid_to_prot_func,
        prot_func_transitions = te.func_transitions,
        sel_prot_func_list = te.sel_prot_func,
        dist_thresh = te.dist_thresh,
        mut_pair_type = te.edge_type, 
        method = te.method, 
        tree_nwk_file = tree_file)

    for mp in mut_pair_set.mut_pair_to_obj:
        mut_pair = mut_pair_set.mut_pair_to_obj[mp]

        if 'fisher' in te.method.lower():
            if mut_pair.p_val_more <= te.thresh:
                mut_pair_to_weight[mp] += \
                                1./te.num_trees
                mut_pair_to_trees[mp] += 1

        else:
            if mut_pair.mod_jaccard > 0.0:
                mut_pair_to_weight[mp] += \
                mut_pair.mod_jaccard / float(te.num_trees)
                mut_pair_to_trees[mp] += 1

    return (mut_pair_to_weight, mut_pair_to_trees)

def get_mut_pair_dicts(tree_ensemble):

    mut_pair_dicts = Parallel(n_jobs=num_cores)\
    (delayed(get_tree_mut_pair_dict)\
     (tree_ensemble, p[0], p[1]) \
     for p in tree_ensemble.prefix_i)
    return mut_pair_dicts
