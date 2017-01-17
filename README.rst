NAPA - Network Analysis of Protein Adaptation
---------------------------------------------

The NAPA package performs intra-protein residue coevolution network construction and analysis.

The package supports two input types:    

1. A FASTA alignment of (functionally related) homologous protein sequences
 
2. A FASTA alignment of (functionally related) homologous protein sequences AND
   a Phylogenetic Tree ensemble (trees in newick format)

NAPA produces the following outputs:

1. A mutation pair network (undirected for alignment only, directed/undirected for phylogeny ensemble)
2. Network analysis - network communities, network node centralities, network path centralities

Documentation
-------------

Please see the `documentation <https://github.com/KarchinLab/NAPA/wiki>`_ on github for complete usage details.

Installation
------------

Use pip to install the napa package:

.. code-block:: bash

    $ pip install napa

NAPA depends on several other packages and these will be installed automatically by pip if they
are not already installed.  

**Required packages:**

* numpy
* scipy
* joblib
* networkx
* ete2
* pyyaml

Note the pip scipy installation can sometimes have issues.  If so, alternate methods of installing scipy may be preferable.  See scipy `install <http://www.scipy.org/install.html>`_ for details.  Then rerun pip install napa.
