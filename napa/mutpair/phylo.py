#! /usr/bin/env python 

# Notes: consider removing ete2/3 dependence
# ete is GPL + only need the PhyloTree class implementation
# https://github.com/etetoolkit/ete/blob/master/ete3/phylo/phylotree.py
# need to be able to read Newick (*.nwk) files and build a 
# general (nonbinary) phylo tree with attributes on nodes and edges 
# FORMAT  DESCRIPTION
# 1       flexible with internal node names


from scipy import stats
from ete2 import PhyloTree

from napa.utils.general import *
from napa.phylo.tree.edge import *
from napa.phylo.tree.tree import *

class PhyloMut(object):
    '''
    Similar to AlnMut class, but counts phylo edges instead of 
    alignmen sequences of occurrence.
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


class PhyloMutPair(object):
    '''
    Pair of muts from tree.
    '''    
    def __init__(self, phylo_edge_pairs = set(), mut_pair = None, 
                 all_edges = set(), prot_func_transitions = {'':''}):

        # List of pairs of phylo_edge_objects
        if not len(phylo_edge_pairs):
            raise ValueError('In PhyloMutPair: ' + \
                             'Set of Phylo edge pairs must be nonempty!')
        else:
            self.phylo_edge_pairs = set(phylo_edge_pairs)

        # A unique pair of AA muts
        if mut_pair ==  None or not len(mut_pair):
            raise ValueError('In PhyloMutPair: No mutations given!')
        else:
            self.mut_pair = mut_pair
            self.name = '%s_%s' % (mut_pair[0].mut_str, mut_pair[1].mut_str)
        if not len(all_edges):
            raise ValueError('In PhyloMutPair: No edge pairs having each ' + \
                             'mutation were given!')
        else:
            self.all_edges = set(all_edges)
        self.prot_func_transitions = prot_func_transitions

    def add_phylo_edge_pair(self, phylo_edge_pair):
        self.phylo_edge_pairs.add(phylo_edge_pair) 
    
    def add_phylo_edge_pairs(self, phylo_edge_pairs):
        self.phylo_edge_pairs.update(phylo_edge_pairs)


    def assign_lineage_distribution(self):
        common_lineage, mut0_lineage, mut1_lineage = [], [], []
        for pep in self.phylo_edge_pairs:
            if pep[0] == pep[1]:
                common_lineage  +=  [pep[0].parent] + [pep[0]] + pep[0].followers 
            else:
                common_prec = pep[0].get_common_precs(pep[1], 'function', 
                                                      self.prot_func_transitions)
                common_follow = list(set([pep[0]] + pep[0].followers) & \
                                    set([pep[1]] + pep[1].followers))
                common_lineage  +=  common_prec + common_follow

                mut0_prec = pep[0].get_intermediate_prec(common_prec, 'function',
                                                         self.prot_func_transitions)
                mut1_prec = pep[1].get_intermediate_prec(common_prec, 'function',
                                                         self.prot_func_transitions)

                mut0_lineage += list(set(mut0_prec + pep[0].followers) - \
                                     set(common_prec + common_follow))
                mut1_lineage += list(set(mut1_prec + pep[1].followers) - \
                                     set(common_prec + common_follow))
        
        self.common_lineage = set([p for p in common_lineage if p.muts])  
        self.mut0_lineage = set([p for p in mut0_lineage if p.muts]) - \
                            self.common_lineage
        self.mut1_lineage = set([p for p in mut1_lineage if p.muts]) - \
                            self.common_lineage
        assoc_lineage = self.common_lineage | self.mut0_lineage | self.mut1_lineage
        self.other_lineage = self.all_edges - assoc_lineage
        #self.print_lineages()


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
        common = len(self.common_lineage)
        only0 = len(self.mut0_lineage) 
        only1 = len(self.mut1_lineage)
        other = len(self.other_lineage)
        self.contingency_table = [[common, only0], [only1, other]]


    def get_contingency_table_edge_counts(self):
        ''' 
        A simpler contingency table based on counts of co-occurrence of 
        muts in tree edges rather than clades. 
        Looking at clade size tends to be more informative.
        '''
        occur01 = len(self.phylo_edge_pairs)
        occur0 = len(self.mut_pair[0].phylo_edge_pairs - self.phylo_edge_pairs) 
        occur1 = len(self.mut_pair[1].phylo_edge_pairs - self.phylo_edge_pairs)
        occur_other = len(self.all_phylo_edge_pairs\
                        -self.mut_pair[0].phylo_edge_pairs \
                        -self.mut_pair[1].phylo_edge_pairs)
        self.contingency_table = [[occur01, occur0], [occur1, occur_other]]
        

    def get_fisher_pval(self):
        '''
        Uses scipy.stats to calculate Fisher's exact p-values based 
        on the contingency table. For association, only need the 
        right-tailed p-value: pval_more.
        '''
        self.assign_lineage_distribution()
        self.get_contingency_table()
        self.odds_ratio_less, self.p_val_less = \
                        stats.fisher_exact(self.contingency_table, 
                                           alternative = 'less')

        self.odds_ratio_more, self.p_val_more = \
                        stats.fisher_exact(self.contingency_table, 
                                           alternative = 'greater')

class PhyloMutPairSet(object):

    def __init__(self, leaf_fasta_file = '', 
                 int_node_fasta_file = '', aln_pos = [],
                 wt_seq_str = '', wt_id = '', pos_subset = [], 
                 leaf_seqid_to_prot_func = {'':''}, 
                 leaf_seqid_to_prot_func_file = '', 
                 int_seqid_to_prot_func = {'':''}, 
                 int_seqid_to_prot_func_file = '', 
                 prot_func_transitions = {},
                 prot_func_transitions_file = '',
                 sel_prot_func_list = [''],
                 sel_prot_func_file = '',
                 dist_thresh = 0.,
                 tree_nwk_file = '',
                 mut_pair_type = 'dir'):
        
        # Get sequences for ALL nodes on tree
        # Sequences from MSA used in phylogeny reconstruction
        # (These end up on tree leaves) AND
        # Reconstructed internal node sequences
        self.aln = BioSeqAln(aln_fasta_file = leaf_fasta_file,
                             aln_pos = aln_pos, 
                             annot_key = 'function', 
                             seqid_to_annot = leaf_seqid_to_prot_func,
                             seqid_to_annot_file = leaf_seqid_to_prot_func_file)
        # Add reconstructed internal node sequences
        self.aln.add_sequences_from_file(int_node_fasta_file)

        # Option to annotate sequences by (reconstructed) function
        # Add reconstructed internal node sequences functions
        self.aln.annot_seqs(annot_key = 'function',
                            seqid_to_annot_file = int_seqid_to_prot_func_file,
                            seqid_to_annot = int_seqid_to_prot_func)

        # Option to focus on protein region(s)/domain(s) of interest
        self.aln.subset_aln_pos(pos_subset = pos_subset)

        # Get mutations from ancestral/wt sequence 
        # for each sequence in alignment
        # option to focus on subset of aln./protein positions
        self.get_wt_seq(wt_id, wt_seq_str)
        self.wt_seq.extract_pos(pos_subset)

        # Python dict representing selected transiotions
        # between reconstructed functions in tree/alignment
        # If empty (default) all functions considered
        # Example, when considering evolution of 
        # extended-spectrum resistance (2be)
        # {'2b':'2be', '2be':'2be', '2br':'2be', '2br':'2ber'}
        self.get_prot_func_transitions(prot_func_transitions, 
                                       prot_func_transitions_file)

        # Names/ids of all tree nodes having selected function
        self.get_sel_func_nodes(annot_key = 'function',
                                annot_list_file = sel_prot_func_file,
                                annot_list = sel_prot_func_list)

        # Get the sequence ids of all leaf sequences annotated
        # with the selected function and find most recent common ancestor
        # This isolates subtree with given function
        self.anc_node, self.sel_func_leaves = \
                    get_mrca_subtree(tree_nwk_file = tree_nwk_file,
                                     fasta_file = leaf_fasta_file,
                                     sel_annot_nodes = self.sel_func_node_ids)


        # Get muts for seqs corresponding to protein function
        self.get_leaf_func_muts()

        # Set of all phylogenetic tree edges
        #self.all_phylo_edge_pairs = set()

        # The folling generates a linegraph object of the whole tree
        # it is a tree rooted in the edge between anc_node and its ancestor
        # nodes of PhyloEdge are edges are in the original tree
        self.anc_edge = PhyloEdge(parent_node = self.anc_node.up, 
                                  child_node = self.anc_node,
                                  aln = self.aln, wt_seq = self.wt_seq)

        self.anc_edge.build_edge_tree()

        # Get set of all pairs of mutations occurring along tree
        # default directed mutation pairs
        self.mut_pair_type = 'undir' if 'undir' in mut_pair_type.lower() \
                             else 'dir'
        self.get_mut_pairs(dist_thresh, self.mut_pair_type)
        self.mut_pairs_fisher_pvalues()

    def get_wt_seq(self, wt_id, wt_seq_str):
        '''
        Obtain WT/ancestral sequence to which all other sequences in 
        the alignment/phylogeny should be compared. Can provide id+sequence,
        or the id of a sequence already in the alignment.
        '''
        if len(wt_id) and len(wt_seq_str) == self.aln.length:
            self.wt_seq = BioSeq(seq_id = wt_id, seq_str = wt_seq_str,
                                 seq_type = 'Protein', 
                                 seq_pos_list = self.aln_pos)
        elif len(wt_id) and len(wt_seq_str) != self.aln.length:
            if wt_id in self.aln.seqid_to_seq:
                #stderr_write(['AlnMutPairSet WARNING:',
                #              'Sequence length does not match alignment,',
                #              'but sequence id in alignment. Replacing',
                #              'with corresp. alignment sequence.'])
                self.wt_seq = self.aln.seqid_to_seq[wt_id]
            else:
                stderr_write(['Invalid  input for WT sequence.',
                              'Please enter valid id already in alignment.',
                               'or a new sequence.'])
                              
        elif not len(wt_id) and len(wt_seq_str) == self.aln.length:
            self.wt_seq = BioSeq(seq_id = 'Wild_type', 
                                 seq_str = wt_seq_str,
                                 seq_type = 'Protein', 
                                 seq_pos_list = self.aln_pos)

        else:
            stderr_write(['Invalid  input for WT sequence.',
                          'Please enter valid id already in alignment.',
                          'or a new sequence.'])


    def get_sel_func_nodes(self, annot_key = 'function',
                           annot_list_file = '',
                           annot_list = ['']):
        if annot_list == [''] and len(annot_list_file):
            annot_list = parse_column(annot_list_file)
        
        self.sel_func_node_ids = self.aln.seqids_with_annot(\
                                    annot_key = annot_key,
                                    annot_list = annot_list)
        
    def get_leaf_func_muts(self):
        '''
        Get muts associated with a function that occur in 
        selected tree leaves.
        This is the subset of all mutations on the tree that is most 
        associated with the function of interest.
        '''
        self.func_muts = []
        for seqid in [leaf.name for leaf in self.sel_func_leaves]:
            self.func_muts += self.wt_seq.get_substitutions(\
                                        self.aln.seqid_to_seq[seqid])
        self.func_muts = set(self.func_muts)

    def get_prot_func_transitions(self, prot_func_transitions = {}, 
                                  prot_func_transitions_file = ''):
        if not len(prot_func_transitions):
            self.prot_func_transitions = \
                        parse_keyval_dlist(prot_func_transitions_file)
        else:
            self.prot_func_transitions = prot_func_transitions

    def check_edge_pair(self,start_edge, follow_edge, dist_thresh):
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
        if not edge.check_node_seq_annot('function', 
                                         self.prot_func_transitions):
            return False

        if not len(set(edge.muts) & self.func_muts):
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
                    PhyloMut(phylo_edge_pairs = set(phylo_edge_pair),
                             mut_str = mut_str)
            else:
                self.mut_str_to_obj[mut_str].add_phylo_edge_pair(phylo_edge_pair)

        for mut_str_pair in zip(mut_list_pair[0], mut_list_pair[1]):
            if mut_str_pair not in self.mut_str_pair_to_obj:
                self.mut_str_pair_to_obj[mut_str_pair] = \
                    PhyloMutPair(phylo_edge_pairs = set([phylo_edge_pair]),
                                 mut_pair = [self.mut_str_to_obj[mut_str_pair[0]],
                                             self.mut_str_to_obj[mut_str_pair[1]]],
                                 all_edges = self.all_edges)
            else:
                self.mut_str_pair_to_obj[mut_str_pair].add_phylo_edge_pair(\
                                                                    phylo_edge_pair)


    def get_mut_pairs(self, dist_thresh, mut_pair_type):
        ''' 
        Extract pairs of mutations occurring along tree edges
        '''
        self.all_edges = [self.anc_edge] + \
                         [fe for fe in self.anc_edge.followers]
        self.all_edges = [e for e in self.all_edges \
                          if self.check_edge(e, dist_thresh)]
        
        stderr_write(['Number of tree edges considered:', 
                      len(self.all_edges)])

        # will hold sets of mutation(s)/ (pairs)
        self.mut_str_pair_to_obj, self.mut_str_to_obj = {}, {}

        for start_edge in self.all_edges:
            start_edge_sorted_muts = sort_by_digits(list(\
                                                    set(start_edge.muts) & \
                                                         self.func_muts))
            if not len(start_edge_sorted_muts): continue
            # Make sure mutation pairs occurring on same edge included
            if len(start_edge_sorted_muts) > 1:
                self.add_phylo_edge_pair(phylo_edge_pair = (start_edge, start_edge), 
                                         mut_list_pair = \
                                         zip(*list_pairs(start_edge_sorted_muts)))
            
            if 'undir' in mut_pair_type.lower():
                self.add_undir_pairs(start_edge, start_edge_sorted_muts, 
                                     dist_thresh)
            else:
                self.add_dir_pairs(start_edge, start_edge_sorted_muts, 
                                   dist_thresh)
                
        stderr_write(['Number of mutations associated with function:', 
                      len(self.mut_str_to_obj.keys())])
        stderr_write(['Number of mutation pairs extracted', 
                      len(self.mut_str_pair_to_obj.keys())])


    def add_undir_pairs(self, start_edge, start_edge_sorted_muts, dist_thresh):
        for follow_edge in start_edge.followers:
            follow_edge_sorted_muts = sorted(list((set(follow_edge.muts) & \
                                         self.func_muts) - set(start_edge.muts)))

            if self.check_edge_pair(start_edge, follow_edge, dist_thresh):
                # For undirected mutation pairs, mutations are ordered by 
                # residue position
                mut_pairs = [(sm, fm) for sm in start_edge_sorted_muts \
                             for fm in follow_edge_sorted_muts \
                             if get_int(sm) <= get_int(fm)]
                self.add_phylo_edge_pair(phylo_edge_pair = (start_edge, follow_edge), 
                                        mut_list_pair = [start_edge_sorted_muts,
                                                         follow_edge_sorted_muts])

            # Reverse order of mutations in cases where follower position
            # is less than start position and distance from follower to start
            # is less than threshold
            if self.check_edge_pair(follow_edge, start_edge, dist_thresh):
                rev_mut_pairs = [(fm, sm) for sm in start_edge_sorted_muts \
                                 for fm in follow_edge_sorted_muts \
                                 if get_int(sm) > get_int(fm)]
                self.add_phylo_edge_pair(phylo_edge_pair = (follow_edge, start_edge), 
                                         mut_list_pair = [follow_edge_sorted_muts,
                                                         start_edge_sorted_muts])


    def add_dir_pairs(self, start_edge, start_edge_sorted_muts, dist_thresh):
        for follow_edge in start_edge.followers:
            follow_edge_sorted_muts = sorted(list((set(follow_edge.muts) & \
                                         self.func_muts) - set(start_edge.muts)))

            if self.check_edge_pair(start_edge, follow_edge, dist_thresh):
                # For undirected mutation pairs, mutations are ordered by 
                # residue position
                mut_pairs = [(sm, fm) for sm in start_edge_sorted_muts \
                             for fm in follow_edge_sorted_muts]
                self.add_phylo_edge_pair(phylo_edge_pair = (start_edge, follow_edge), 
                                         mut_list_pair = [start_edge_sorted_muts,
                                                          follow_edge_sorted_muts])


    def mut_pairs_fisher_pvalues(self):
        '''
        Calculate Fisher's Exact Test p-value, 
        for each PhyloMutPair object
        '''
        for mut_str_pair in self.mut_str_pair_to_obj:
            mut_pair = self.mut_str_pair_to_obj[mut_str_pair]
            mut_pair.get_fisher_pval()


    def __repr__(self):
        out_str = '\t'.join(['_Mut1_', '_Mut2_', 'edgePairList', 
                          'contTable(common_1_2_other)', 
                          "FisherPvalLess", "FisherPvalMore"]) + '\n'
        for mut_str_pair in self.mut_str_pair_to_obj:
            mut_pair = self.mut_str_pair_to_obj[mut_str_pair]
            if len(flatten(mut_pair.contingency_table)) < 4:
                continue
            if mut_pair.contingency_table[0][1] == 0:
                if mut_pair.contingency_table[1][0] == 0:
                    if mut_pair.contingency_table[0][0] < 2:
                        continue
            to_print = [mut_pair.name]
            to_print.append(str(len(mut_pair.phylo_edge_pairs)))
            to_print.append(str(mut_pair.contingency_table))
            to_print.append(str(mut_pair.p_val_less))
            to_print.append(str(mut_pair.p_val_more))
            out_str +=  '\t'.join(to_print) + '\n'
        return out_str

    

class CombinedPhyloMutPairs(object):

    def __init__(self, phylo_mut_pair_set, pval_thresh = 0.05):
        self.pval_thresh = pval_thresh
        self.mut_pair_to_tree_num = defaultdict(int)
        self.add_mut_pairs(phylo_mut_pair_set)
    

    def add_mut_pairs(self, phylo_mut_pair_set):
        '''
        add mutation pairs from a set calculated for a single tree
        all mutations passing the Fisher's exact test with a significant 
        right-tailed p-value are added
        '''
        for mut_str_pair in phylo_mut_pair_set.mut_str_pair_to_obj:
            if phylo_mut_pair_set.mut_str_pair_to_obj[mut_str_pair].p_val_more <= \
               self.pval_thresh:
                self.mut_pair_to_tree_num[mut_str_pair] += 1


    def write_phylo_network(self):
        out_str = ''
        for mut_pair in self.mut_pair_to_tree_num:
            out_str += '%s\t%s\t%d\n'%(mut_pair[0], mut_pair[1], 
                                    self.mut_pair_to_tree_num[mut_pair])
        

    def __repr__(self):
        out_str = ''
        for mut_pair in self.mut_pair_to_tree_num:
            out_str += '%s_%s\t%d\n'%(mut_pair[0], mut_pair[1], 
                                    self.mut_pair_to_tree_num[mut_pair])
        return out_str
        
               
