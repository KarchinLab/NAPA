#! /bin/bash
#$ -S /bin/bash
#$ -N cogentAncReconstTest
#$ -o /projects/epistasis/intragene-epistasis/TEM-network_dir/results/phylogeny/ancestors/exprs/log
#$ -e /projects/epistasis/intragene-epistasis/TEM-network_dir/results/phylogeny/ancestors/exprs/log 
#$ -M vbeleva@gmail.com
#$ -m abes
#$ -cwd     

hostname
date
currDir=/projects/epistasis/intragene-epistasis/TEM-network_dir/results/phylogeny/ancestors
input=$currDir/${commandargs}  #Tells it which file to use for individual task commands 
pythoncode=$currDir/runPyCogentAncest.py
python $pythoncode -i $SGE_TASK_ID -a $input
date
