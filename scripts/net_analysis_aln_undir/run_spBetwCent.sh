#./run_nxPathCent.sh 2015-07-13 0.01 3 

projDir=/home/vbeleva/phylogeny/pairs
workDir=$projDir/alignmentPairs
summDir=$workDir/results

fileBases=( "TEM_undirPairsPhenoClin_nx" "TEM_undirPairsPhenoClinDirevo_nx" \
    "TEM_undirPairsPhenoClin_Ji_nx" "TEM_undirPairsPhenoClinDirevo_Ji_nx" )


lengths=( 1 2 3 4 5 )

for f in "${fileBases[@]}"
do
    for l in "${lengths[@]}"
    do
	netFileName=$summDir/$f.txt
	nohup python MutNetAnalysis.py -d $netFileName -f 1 -l $l  >& SPpathCents.l$l.oe &
    done
done


