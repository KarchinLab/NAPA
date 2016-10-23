from ete2 import PhyloTree


def load_tree_sequences(nwk_file, fasta_file):
    ''' 
    Load a tree with associated sequences on leaves. 
    '''
    tree = PhyloTree(newick = nwk_file, format = 1)
    tree.link_to_alignment(alignment = fasta_file, 
                           alg_format = 'fasta')
    return tree


def filter_leaves(tree, node_name_list):
    if not len(node_name_list):
        # no filtering -- all leaves returned
        return tree.get_leaves() 
    return filter(lambda node: node.name in node_name_list, 
                  tree.get_leaves())
    
def get_mrca(tree, nodes):
    return tree.get_common_ancestor(nodes)

def get_mrca_leaf_names(tree, node_name_list):
    node_subset = filter_leaves(tree, node_name_list)
    return get_mrca(tree, node_subset)

def get_mrca_subtree(tree_nwk_file, fasta_file, 
                     sel_annot_nodes):
    ''' 
    Get subtree whose leaves all have the same annotation
    Typically finds subtree where all sequences evolve same 
    function
    '''
    tree = load_tree_sequences(tree_nwk_file, fasta_file)

    # Get list of tree leaf objects with selected function
    sel_annot_leaves = filter_leaves(tree, sel_annot_nodes)

    # Get common ancestor for leaves with function of interest
    anc_node = get_mrca(tree, sel_annot_leaves)
    return anc_node, sel_annot_leaves
