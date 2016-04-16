# NAPA-dev
Network Analysis of Protein Adaptation (intra-protein residue coevolution network construction and analysis)
--------------------------------------------- 
Creates and analyses networks of coevolving protein mutations

Input: Protein alignments (fasta) or protein phylogenetic tree
       ensemble with a fasta file for each tree

Intermediate input: undirected mutation pairs and coevolution weights
       for each pair directed (only when phylogeny available) pairs of
       mutations with association weights

Output: Network of coevolving pairs for further analysis/visualization
	Central paths in the network shortest path betweenness
	centrality k-path "random walk" path betweenness centrality

--------------------------------------------- 
Alignment analysis pipeline:

1. Mutation pairs from alignment: Obtain linked (undirected) pairs of
   mutations and their links' weights Note: Can skip this step when you
   have your own input of mutation pairs and weights

scripts/get_mut_pairs/from_alignment 
./run_undirPairs_single.sh
./run_summarize_pairs.sh


2. Path centralities: Use python packages networkx, igraph, and
   modified kpath.py scripts to analyze path centralities (python
   packages networkX and iGraph)

scripts/net_analysis_aln_undir 
./run_spBetwCent.sh 
./run_nxPathCent.sh

3. (Network vizualization in CytoScape)
See example in:
test_results/network_visualization

-------------------------------------------- 
Phylogeny analysis pipeline:

0. Phylogeny reconstruction and preprocessing [to come]


1. Mutation pairs from phylogeny 1a. Undirected pairs of mutations
Note: Can skip this step when you have your own input of mutation
pairs and weights 

scripts/get_mut_pairs/from_phylo_undir
./run_undirPairs_single.sh 
./run_summarize_pairs.sh

1b. Directed pairs of mutations 

scripts/get_mut_pairs/from_phylo_undir
./run_dirPairs_single.sh 
./run_summarize_dir_pairs.sh


2. Path centralities

2a. Undirected network 

scripts/net_analysis_phylo_undir
./run_spBetwCent.sh 
./run_nxPathCent.sh

2b. Directed network 
scripts/net_analysis_phylo_dir
./run_spBetwCent.sh 
./run_nxPathCent.sh

3. (Network vizualization in CytoScape)

