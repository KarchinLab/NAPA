#!/bin/bash

projDir=$1
projDir=$projDir/NAPA
dataDir=$projDir/tests/data
resultsDir=$projDir/tests/results
scriptsDir=$projDir/tests/scripts

pathLen=3

# TEM-phylo-net.full is the full version of the phylo directed net:
# it cannot be generated with current test data, but it can be used 
# for performance scaling analysis

# All 'rel' relative centrality algorithms use recursion, hence are 
# slow and memory intensive: relloc and rellglob split centralities into 
# global and local for faster calculation.

mutPairTypes=( "undir"      "undir"         "dir"           "dir")
netBaseList=( "TEM-aln-net" "TEM-phylo-net" "TEM-phylo-net" "TEM-phylo-net.full")
centTypeList=( "node" "path" "rel" "relloc" "relglob" )

i=3
j=2
mutPairType="${mutPairTypes[$i]}"
netBase="${netBaseList[$i]}"
centType="${centTypeList[$j]}"
netFile=$resultsDir/$netBase.$mutPairType.txt
ls $netFile
centOutFile=$resultsDir/$netBase.$mutPairType.$centType-cent.$pathLen.txt
echo $i $j $mutPairType $centType
echo $centOutFile

python $scriptsDir/net_analysis.py -nf $netFile -mt $mutPairType \
    -pl $pathLen -co $centOutFile -ct $centType

