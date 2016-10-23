work_dir=`pwd`
napa_dir=${work_dir%/*}

#----------------------------------------#
# To build the undirected TEM 
# alignment based network
./scripts/build_net.run.sh ${napa_dir} config/aln.config.yaml

# To perform graph-theoretical
# analysis of the TEM directed
# alignment based network
./scripts/net_analysis.run.sh ${napa_dir} config/aln.config.yaml
#----------------------------------------#

#----------------------------------------#
# To build the directed TEM 
# phylogeny based network
./scripts/build_net.run.sh ${napa_dir} config/phylo.dir.config.yaml

# To perform graph-theoretical
# analysis of the TEM directed 
# phylogeny based network
./scripts/net_analysis.run.sh ${napa_dir} config/phylo.dir.config.yaml
#----------------------------------------#

#----------------------------------------#
# To build the undirected TEM 
# phylogeny based network
./scripts/build_net.run.sh ${napa_dir} config/phylo.undir.config.yaml

# To perform graph-theoretical
# analysis of the undirected TEM 
# phylogeny based network
./scripts/net_analysis_run.sh ${napa_dir} config/phylo.undir.config.yaml
#----------------------------------------#
