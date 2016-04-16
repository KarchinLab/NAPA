#./run_nxPathCent.sh 2015-07-13 0.01 3 

projDir=/home/vbeleva/phylogeny/pairs
workDir=$projDir/alignmentPairs
summDir=$workDir/results

#fileBases=( "TEM_undirPairsPhenoClin" "TEM_undirPairsPhenoClinDirevo" \
#    "TEM_undirPairsPhenoClin_Ji" "TEM_undirPairsPhenoClinDirevo_Ji" )
fileBases=( "TEM_undirPairsPhenoClin" "TEM_undirPairsPhenoClin_Ji" )


lengths=( 1 2 3 4 )

for f in "${fileBases[@]}"
do
    for l in "${lengths[@]}"
    do
	netFileName=$summDir/$f.txt
	nohup python MutNetAnalysis_networkx.py -d $netFileName -f 0 -l $l  >& kpathCents.l$l.oe &
    done
done


