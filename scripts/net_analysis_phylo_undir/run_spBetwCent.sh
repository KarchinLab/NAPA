
projDir=/home/vbeleva/phylogeny/pairs
workDir=$projDir/undirectedPairs
summDir=$workDir/results/summary

fileBase=2016-02-26.undirPhenoPair


lengths=( 1 2 3 4 )
pvalThreshs=( 0.01 )
distThreshs=( 0.05 )


for p in "${pvalThreshs[@]}"
do
    for d in "${distThreshs[@]}"
    do
	netFileName=$summDir/$fileBase$p.dist$d.txt
	echo $netFileName
	for l in "${lengths[@]}"
	do
	    nohup python MutNetAnalysis.py -d $netFileName -f 1 -l $l  >& spCents.p$p.d$d.l$l.oe &
	done
   done
done

