# NAPA-dev
Network Analysis of Protein Adaptation (intra-protein residue coevolution network construction and analysis)
--------------------------------------------- 
Creates and analyses networks of coevolving protein mutations

Inputs (2 types): 

A FASTA alignment of (functionally related) homologous protein sequences 

OR

A FASTA alignment of (functionally related) homologous protein sequences
Optional: A set of function categories for sequences in alignment
A Phylogenetic Tree ensemble (trees in newick format)
A set of internal node FASTA sequences for each tree
Optional: A set of internal node functions


Outputs:

A mutation pair network 
  undirected for alignment only
  undirected or directed for phylogeny ensemble input

Network analysis
   network communities
   network node(=mutation) centralities
   network path(=sets of mutations) centralities



--------------------------------------------- 
Sample Extended Workflow 
--------------------------------------------- 
1. Sequence collection (functinal protein orthologs)
2. Sequence alignment methods (e.g. PyCogent)
3. MrBayes (phylogeny) set-up
4. MrBayes run
5. Post-process MrBayes equilibrated ensemble
   -burnin removal
   -thinning
6. Ancestral reconstruction of internal node sequences
7. Ancestral reconstruction of internal node functions
8. Run NAPA network construction and analysis of paths


