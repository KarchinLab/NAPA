from ete2 import PhyloTree
form NAPA.phylo.tree.edge import *
from NAPA.phylo.tree.tree import *

class PhyloMut(object):
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
                 all_edges = set()):

        # List of pairs of phylo_edge_objects
        if not len(phylo_edge_pairs):
            raise ValueError('Set of Phylo edge pairs must be nonempty!')
        else:
            self.phylo_edge_pairs = set(phylo_edge_pairs)

        # A unique pair of AA muts
        if mut_pair ==  None:
            print "WARNING: No mut pair given. Set to" + \
                " first pair of muts for first " + \
                " Phylo_edge pair!"
            self.mut_pair = (phylo_edge_pairs[0][0].get_muts()[0], 
                                 phylo_edge_pairs[0][1].get_muts()[1])
        else:
            self.mut_pair = mut_pair

        self.all_edges = all_edges



    def assign_lineage_distribution(self):

        common_lineage, mut0_lineage, mut1_lineage = [], [], []
        for pep in self.phylo_edge_pairs:
            if pep[0] =  = pep[1]:
                common_lineage  +=  [pep[0].parent] + [pep[0]] + pep[0].followers 
            else:
                common_prec = pep[0].get_common_precursors(pep[1])
                common_follow = list(set([pep[0]] + pep[0].followers) & \
                                    set([pep[1]] + pep[1].followers))
                common_lineage  +=  common_prec + common_follow

                mut0_prec = pep[0].get_intermediate_prec(common_prec)
                mut1_prec = pep[1].get_intermediate_prec(common_prec)

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

    
    def print_lineages(self):
        print "\tcommonLin:", [p.name for p in self.common_lineage]
        print "\tmut0_lin:", [p.name for p in self.mut0_lineage]
        print "\tmut1_lin:", [p.name for p in self.mut1_lineage]
        print "\totherLin:", [p.name for p in self.other_lineage]
        print 


    def get_contingency_table(self):

        common = len(self.common_lineage)
        only0 = len(self.mut0_lineage) 
        only1 = len(self.mut1_lineage)
        other = len(self.other_lineage)
        self.contingency_table = [[common, only0], [only1, other]]


    def get_contingency_table_edge_counts(self):
        ''' A simpler contingency table based on counts of co-occurrence of 
        muts in tree edges rather than clades. '''

        occur01 = len(self.phylo_edge_pairs)
        occur0 = len(self.mut_pair[0].phylo_edge_pairs - self.phylo_edge_pairs) 
        occur1 = len(self.mut_pair[1].phylo_edge_pairs - self.phylo_edge_pairs)
        occur_other = len(self.all_phylo_edge_pairs\
                        -self.mut_pair[0].phylo_edge_pairs \
                        -self.mut_pair[1].phylo_edge_pairs)
        self.contingency_table = [[occur01, occur0], [occur1, occur_other]]
        

    def get_fisher_pvals(self):
        self.get_contingency_table()
        self.odds_ratio_less, self.p_val_less = \
                        stats.fisher_exact(self.contingency_table, 
                                           alternative = 'less')

        self.odds_ratio_more, self.p_val_more = \
                        stats.fisher_exact(self.contingency_table, 
                                           alternative = 'greater')



class PhyloMutPairUndir(object):

    def __init__(self, nwk_file, fasta_file, pos_list, pos_to_wt, 
                 pos_subset, tip_to_prot_func, sel_prot_funcs, 
                 int_prot_func_file, dist_thresh):

        # Newick and corresponding FASTA aln. of leaf sequences
        self.nwk_file = nwk_file
        self.fasta_file = fasta_file

        # List of positins in a protein, e.g 3,4,5...296
        self.pos_list = pos_list

        # Wildtype residue at a position
        self.pos_to_wt = pos_to_wt

        # Name (id) of each node in tree with associated function
        # Protein function can be added to each node as an attribute instead
        self.node_name_to_prot_func = tip_to_prot_func.copy()

        # List of selected protein functions in tree/alignment
        # If empty (default) all functions considered
        self.sel_prot_funcs = sel_prot_funcs

        # Names/ids of leaf nodes having selected function
        self.prot_func_leaf_names = [node_name \
                                     for node_name in self.node_name_to_prot_func \
                                     if self.node_name_to_prot_func[node_name] in \
                                     sel_prot_funcs]

        # Set of all phylogenetic tree edges
        self.all_phylo_edge_pairs = set()

        # Get muts along phylogenetic tree edges
        self.get_internal_node_seq()

        # Get muts for seqs corresponding to protein function
        self.get_prot_func_muts()

        # Get reconstructed functions for phylo tree internal nodes and
        # add them to known functions of tip nodes
        self.node_to_prot_func.update(\
                            read_keyval_dict(int_prot_func_file))

        # The folling generates a linegraph object of the whole tree
        # it is a tree rooted in the edge between anc_node and its ancestor
        # nodes of PhyloEdge are edges are in the original tree
        self.anc_edge = PhyloEdge(node_pair = [self.anc_node.up, self.anc_node], 
                                  prot_func_dict = self.node_to_prot_func.copy(), 
                                  parent = None, children = [], 
                                  pos_list = self.pos_list, 
                                  pos_to_wt = self.pos_to_wt)
        self.anc_edge.build_edge_tree()

        # Now get all possible undirected pairs
        self.undir_pairs = []
        self.get_undir_pairs(dist_thresh)


    def get_internal_node_seqs(self):
        ''' 
        Read the tree from newick(nwk) file and link nodes
        to their corresponding fasta seqs (uses ete2). 
        Part could be moved to NAPA.phylo.tree module
        '''
        phylo_tree = load_tree_sequences(self.nwk_file, self.fasta_file)

        # Get list of phylo node objects with selected function
        # It starts with leaf nodes
        self.prot_func_leaves = filter_leaves(phylo_tree, 
                                             self.prot_func_leaf_names)

        # Get common ancestor for leaves with function of interest
        self.anc_node = self.get_mrca(phylo_tree, self.prot_func_nodes)


    def get_prot_func_muts(self):
        '''
        Get muts associated with a function that occur in selected tree leaves.
        Uses ete2.
        '''
        standard_aa = "ACDEFGHIKLMNPQRSTVWY" # list of standard amino acids

        muts = set()
        anc_seq = [self.pos_to_wt[p] for p in self.pos_list if p in self.pos_to_wt]
        for node in self.prot_func_leaves:
            node_seq = node.sequence
            muts.update(set([anc_seq[i] + str(self.pos_list[i]) + node_seq[i] \
                             for i in range(len(node_seq)) \
                             if anc_seq[i]! = node_seq[i] \
                             and node_seq[i] in standard_aa \
                             and self.pos_list[i] in self.pos_subset]))
        self.prot_func_muts = muts


    def get_undir_pairs(self, dist_thresh):
        ''' 
        Extract pairs of mutations occurring along tree edges
        '''
        mut_pair_to_phylo_edge_pairs = defaultdict(list) 
        mut_to_phylo_edge_pairs = defaultdict(list)
        all_edges = [self.anc_edge] + [fe for fe in self.anc_edge.followers]
        stderr_write(["\t\t_number of tree edges:", 
                      len(all_edges)])
        
        stderr_write(["\t\t_assigning ordered precurors for all edges"])
        for edge in all_edges:
           edge.assign_ordered_precursors()
        
       
        for start_edge in self.anc_edge.followers:
            if start_edge.is_leaf_edge(): continue
            if not start_edge.check_prot_func(): continue

            if not start_edge.check_dist(start_edge, dist_thresh): continue
            

            start_edge_sorted_muts = sort_by_digits(list(\
                                                    set(start_edge.muts) & \
                                                         self.prot_func_muts))
            if not len(start_edge_sorted_muts): continue

            #stderr_write([start_edge.name, start_edge.muts])
            if len(start_edge_sorted_muts)>1:
                start_edge_mutation_pairs = list_pairs(start_edge_sorted_muts)
                # get muts within same edge for all anc_edge + followers
                for mut_pair in start_edge_mutation_pairs:
                    mut_pair_to_phylo_edge_pairs[mut_pair].append((start_edge, start_edge))
                    mut_to_phylo_edge_pairs[mut_pair[0]].append((start_edge, start_edge))
                    mut_to_phylo_edge_pairs[mut_pair[1]].append((start_edge, start_edge))
            
            # get mutation pairs for all start, follow edge combinations
            for follow_edge in start_edge.followers:
                #if follow_edge.is_leaf_edge(): continue
                if not follow_edge.check_prot_func(): continue
                if not start_edge.check_dist(follow_edge, dist_thresh): continue
               
                
                follow_edge_sorted_muts = list((set(follow_edge.muts) &\
                                             self.prot_func_muts)-\
                                            set(start_edge.muts))
                
                if not len(follow_edge_sorted_muts): continue
        
                for start_mut in start_edge_sorted_muts:
                    for follow_mut in follow_edge_sorted_muts:
                        if get_int(start_mut) < get_int(follow_mut):
                            mut_pair_to_phylo_edge_pairs[(start_mut, follow_mut)].\
                            append((start_edge, follow_edge))
                            mut_to_phylo_edge_pairs[start_mut].\
                                append((start_edge, follow_edge))
                            mut_to_phylo_edge_pairs[follow_mut].\
                                append((start_edge, follow_edge))

                        elif get_int(start_mut) > get_int(follow_mut):
                            mut_pair_to_phylo_edge_pairs[(follow_mut, start_mut)].\
                                append((follow_edge, start_edge))
                            mut_to_phylo_edge_pairs[start_mut].\
                                append((follow_edge, start_edge))
                            mut_to_phylo_edge_pairs[follow_mut].\
                                append((follow_edge, start_edge))


        stderr_write(["muts and prot_func muts", mut_to_phylo_edge_pairs.keys()])
        all_muts = set(mut_to_phylo_edge_pairs.keys()) & set(self.prot_func_muts)
        stderr_write(["\t\t_number of AA muts:", len(all_muts)])
        
        all_a_amut_edges = set([edge for edge in self.anc_edge.followers \
                             if edge.check_prot_func()])
        stderr_write(["\t\t_number of edges with non-silent mutations for protein function:", 
                      len(all_a_amut_edges)])
 
        self.mut_to_obj = {}
        for mut_pair in mut_pair_to_phylo_edge_pairs:
            if mut_pair[0] not in self.mut_to_obj:
                self.mut_to_obj[mut_pair[0]] = Mutation(mutation_str = mut_pair[0])
            if mut_pair[1] not in self.mut_to_obj:
                self.mut_to_obj[mut_pair[1]] = Mutation(mutation_str = mut_pair[1])
            
            self.mut_to_obj[mut_pair[0]].add_phylo_edge_pairs(set(\
                                                mut_pair_to_phylo_edge_pairs[mut_pair]))
            self.mut_to_obj[mut_pair[1]].add_phylo_edge_pairs(set(\
                                                mut_pair_to_phylo_edge_pairs[mut_pair]))

            undir_mut_pair = Mutation_pair(phylo_edge_pairs = list(set(\
                                        mut_pair_to_phylo_edge_pairs[mut_pair])), 
                                        mutation_pair = (self.mut_to_obj[mut_pair[0]], 
                                                        self.mut_to_obj[mut_pair[1]]), 
                                        all_edges = all_a_amut_edges)
            undir_mut_pair.assign_lineage_distribution()
            self.undir_pairs.append(undir_mut_pair)

        stderr_write(["\t\t_pairwise lineage dsitributions assigned."])
        stderr_write(["\t\t_build contingencies."])
        for mut_pair in self.undir_pairs:
            mut_pair.get_contingency_table()
                        

    def __repr__(self):
        out_str = '\t'.join(['_Mut1_', '_Mut2_', 'edgePairList', 
                          'contTable(common_1_2_other)', 
                          "FisherPvalLess", "FisherPvalMore"]) + '\n'
        
        for undir_pair_idx in range(len(self.undir_pairs)):
            undir_pair = self.undir_pairs[undir_pair_idx]
            
            if len(undir_pair.contingency_table)<_to_:
                continue
            
            if len(flatten(undir_pair.contingency_table))<4:
                continue
            if undir_pair.contingency_table[0][1] =  = 0:
                if undir_pair.contingency_table[1][0] =  = 0:
                    if undir_pair.contingency_table[0][0]< = 2:
                        continue
           
            to_print = [str(undir_pair.mutation_pair[0]), 
                       str(undir_pair.mutation_pair[1])]
            to_print.append(str(len(undir_pair.phylo_edge_pairs)))
            to_print.append(str(undir_pair.contingency_table))
            to_print.append(str(undir_pair.p_val_less))
            to_print.append(str(undir_pair.p_val_more))
            out_str +=  '\t'.join(toPrint) + '\n'

        return out_str

# Should set this up to inherit from PhyloMutPairSetUndir
class PhyloMutPairSetDir(object):

    def __init__(self, nwk_file, fasta_file, pos_list, pos_to_wt, tip_to_pheno, pheno_clades, 
                 dist_thresh):
        self.nwk_file = nwk_file
        self.fasta_file = fasta_file
        self.pos_list = pos_list
        self.pos_to_wt = pos_to_wt
        self.node_to_pheno = tip_to_pheno.copy()        
        self.all_phylo_edge_pairs = set()

        # read the tree from newick(nwk) file and link nodes to corresp. fasta seqs
        self.phylo_tree = Phylo_tree(newick = self.nwk_file, format = 1)
        self.phylo_tree.link_to_alignment(alignment = self.fasta_file, alg_format = 'fasta')
        
        #stderr_write(pheno_clades)
        self.pheno_leaf_nodes = filter(lambda phylo_node: \
                                     phylo_node.name in flatten(pheno_clades.values()), 
                                self.phylo_tree.get_leaves())
        self.pheno_name_nodes  = filter(lambda phylo_node: \
                                     phylo_node.name in pheno_clades["_to_be"], 
                                self.phylo_tree.get_leaves())
        self.get_pheno_muts()
        self.anc_node = self.phylo_tree.get_common_ancestor(self.pheno_name_nodes)
        stderr_write(['\t_mrca 2be', self.anc_node.name, 
                          '\t' + str(datetime.now())])

        int_pheno_file =  nwk_file.replace("int_nwk.tree", "internalStates.txt")
        int_pheno_file =  int_pheno_file.replace("int_sub_pheno.tree", "internalStates.txt")
        with open(int_pheno_file, 'rb') as pf:
            stderr_write(["internal node pheno", int_pheno_file])
            line_recs = [line.strip().split() \
                        for line in pf.read().replace('"', '').split('\n')]
            for recs in line_recs:
                if len(recs)>1:
                    self.node_to_pheno[recs[0]] = recs[1]
        
        stderr_write(["nodes and internal nodes with phenotype", len(self.node_to_pheno)])

        # The folling generates a linegraph object of the whole tree
        # it is a tree rooted in the edge between anc_node and its ancestor
        self.anc_edge = Phylo_edge(node_pair = [self.anc_node.up, self.anc_node], 
                                 pheno_dict = self.node_to_pheno.copy(), 
                                 parent = None, children = [], pos_list = self.pos_list, 
                                 pos_to_wt = self.pos_to_wt)
        self.anc_edge.build_edge_tree()

        # Now get all possible directed pairs
        self.dir_pairs = []
        self.get_dir_pairs(dist_thresh)
        stderr_write(["\t_finished extracting pairs:", self.nwk_file.split('/')[-1], 
                      datetime.now()])


    def get_pheno_muts(self):
        standard_aa = "ACDEFGHIKLMNPQRSTVWY" # list of standard amino acids

        mutations = []
        for node in self.pheno_name_nodes:
            node_seq = node.sequence
            anc_seq = [self.pos_to_wt[p] for p in self.pos_list if p in self.pos_to_wt]
            mutations  +=   [anc_seq[i] + str(self.pos_list[i]) + node_seq[i] \
                        for i in range(len(node_seq)) \
                               if anc_seq[i]! = node_seq[i] and \
                           node_seq[i] in standard_aa and self.pos_list[i]>_to_3]
        self.pheno_mutations = set(mutations)
        stderr_write(["phenotype mutations"] + list(self.pheno_mutations))

    def get_dir_pairs(self, dist_thresh):
        mut_pair_to_phylo_edge_pairs = defaultdict(list) 
        mut_to_phylo_edge_pairs = defaultdict(list)
        all_edges = [self.anc_edge] + self.anc_edge.followers
        for edge in all_edges:
            edge.assign_ordered_precursors()
        stderr_write(["\t\t_number of tree edges:", len(all_edges)])

        for start_edge in self.anc_edge.followers:
            if start_edge.is_leaf_edge(): continue
            if not start_edge.check_pheno(): continue
            if start_edge.check_dist(start_edge, dist_thresh):
                start_edge_sorted_mutations = sort_by_digits(list(\
                                                    set(start_edge.mutations) & \
                                                        self.pheno_mutations))
                if not len(start_edge_sorted_mutations): continue
                if len(start_edge_sorted_mutations)>1:
                    start_edge_mutation_pairs = list_pairs(start_edge_sorted_mutations)
                    start_edge_sorted_mutations.reverse()
                    start_edge_mutation_pairs  +=  list_pairs(start_edge_sorted_mutations)
                    for mut_pair in start_edge_mutation_pairs:
                        if get_int(mut_pair[0]) ==  get_int(mut_pair[1]): continue
                        mut_pair_to_phylo_edge_pairs[mut_pair].append((start_edge, start_edge))
                        mut_to_phylo_edge_pairs[mut_pair[0]].append((start_edge, start_edge))
                        mut_to_phylo_edge_pairs[mut_pair[1]].append((start_edge, start_edge))

            # get mutation pairs for all start, follow edge combinations
            for follow_edge in start_edge.followers:
                #if follow_edge.is_leaf_edge(): continue
                if (not follow_edge.check_pheno()):
                    continue
                if not start_edge.check_dist(follow_edge, dist_thresh):
                    continue

                follow_edge_sorted_mutations = list((set(follow_edge.mutations) &\
                                             self.pheno_mutations)-\
                                            set(start_edge.mutations)) 

                if not len(follow_edge_sorted_mutations): continue
                for start_mut in start_edge_sorted_mutations:
                    for follow_mut in follow_edge_sorted_mutations:
                        if get_int(start_mut) ==  get_int(follow_mut):
                            continue
                        mut_pair_to_phylo_edge_pairs[(start_mut, follow_mut)].\
                            append((start_edge, follow_edge))
                        mut_to_phylo_edge_pairs[start_mut].\
                                append((start_edge, follow_edge))
                        mut_to_phylo_edge_pairs[follow_mut].\
                                append((start_edge, follow_edge))
        
        all_mutations = set(mut_to_phylo_edge_pairs.keys()) & set(self.pheno_mutations)
        stderr_write(["\t\t_number of AA mutations for pheno:", len(all_mutations)])

        #stderr_write(list(set(flatten(mut_pair_to_phylo_edge_pairs.keys()))))
        stderr_write(["\t\t_number of AA mutations for phenotype:", 
                      len(all_mutations)])
        
        all_a_amut_edges = set([edge for edge in all_edges if edge.check_pheno()])
        stderr_write(["\t\t_number of edges with non-silent mutations for phenotype:", 
                      len(all_a_amut_edges)])
 
        self.mut_to_obj = {}
        for mut_pair in mut_pair_to_phylo_edge_pairs:
            if mut_pair[0] not in self.mut_to_obj:
                self.mut_to_obj[mut_pair[0]] = Mutation(mutation_str = mut_pair[0])
            if mut_pair[1] not in self.mut_to_obj:
                self.mut_to_obj[mut_pair[1]] = Mutation(mutation_str = mut_pair[1])
            
            self.mut_to_obj[mut_pair[0]].add_phylo_edge_pairs(set(\
                                                mut_pair_to_phylo_edge_pairs[mut_pair]))
            self.mut_to_obj[mut_pair[1]].add_phylo_edge_pairs(set(\
                                                mut_pair_to_phylo_edge_pairs[mut_pair]))

            mut_pair = Mutation_pair(phylo_edge_pairs = list(set(\
                                        mut_pair_to_phylo_edge_pairs[mut_pair])), 
                                        mutation_pair = (self.mut_to_obj[mut_pair[0]], 
                                                        self.mut_to_obj[mut_pair[1]]), 
                                        all_edges = all_a_amut_edges)
            mut_pair.assign_lineage_distribution()
            self.dir_pairs.append(mut_pair)

        stderr_write(["\t\t_pairwise lineage dsitributions assigned"])
        stderr_write(["\t\t_build contingencies."])
        for mut_pair in self.dir_pairs:
            mut_pair.get_contingency_table()
                
    def __repr__(self):
        out_str = '\t'.join(['_Mut1_', '_Mut2_', 'edgePairNum', 
                          'contTable(common_1_2_other)', 
                          "FisherPvalLess", "FisherPvalMore"]) + '\n'
        
        for dir_pair_idx in range(len(self.dir_pairs)):
            dir_pair = self.dir_pairs[dir_pair_idx]

            if len(dir_pair.contingency_table)<_to_:
                continue
            
            if len(flatten(dir_pair.contingency_table))<4:
                continue
            if dir_pair.contingency_table[0][1] =  = 0:
                if dir_pair.contingency_table[1][0] =  = 0:
                    if dir_pair.contingency_table[0][0]< = 2:
                        continue

            to_print = [str(dir_pair.mutation_pair[0]), str(dir_pair.mutation_pair[1])]
            to_print.append(str(len(dir_pair.phylo_edge_pairs)))
            to_print.append(str(dir_pair.contingency_table))
            to_print.append(str(dir_pair.p_val_less))
            to_print.append(str(dir_pair.p_val_more))
            out_str +=  '\t'.join(toPrint) + '\n'

        return out_str
    
                                         
