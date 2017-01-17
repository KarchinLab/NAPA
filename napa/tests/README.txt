work_dir=`pwd`
# napa directory and run file
napa_dir=${work_dir%/*}
napa_run=$napa_dir/run_napa.py

#----------------------------------------#
# Generally run this way
#analysis_command=build #analyze
#config_file=<path/to/my_config.yaml>
#./$napa_run -r $analysis_command -c $config_file

#----------------------------------------#
# To build the undirected TEM 
# alignment based network
./$napa_run -r build -c config/aln.config.yaml


# To perform graph-theoretical
# analysis of the TEM directed
# alignment based network
./$napa_run -r analyze -c config/aln.config.yaml
#----------------------------------------#


#----------------------------------------#
# To build the directed TEM 
# phylogeny based network
./$napa_run -r build -c config/phylo.dir.config.yaml

# To perform graph-theoretical
# analysis of the TEM directed 
# phylogeny based network
./$napa_run -r analyze -c config/phylo.dir.config.yaml
#----------------------------------------#


#----------------------------------------#
# To build the undirected TEM 
# phylogeny based network
./$napa_run -r build -c config/phylo.undir.config.yaml


# To perform graph-theoretical
# analysis of the undirected TEM 
# phylogeny based network
./$napa_run -r analyze -c config/phylo.undir.config.yaml
#----------------------------------------#
