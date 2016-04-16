projDir=/home/vbeleva/phylogeny/pairs
workDir=$projDir/directedPairs
summDir=$workDir/results/summary/directedPairs

fileBase="2016-02-18"


lengths=( 1 2 3 4 )
pvalThreshs=( 0.01 0.05 0.1 )
distThreshs=( 0.05 0.1 )


for p in "${pvalThreshs[@]}"
do
    for d in "${distThreshs[@]}"
    do
	netFileName=$summDir/$fileBase".directPhenoPair_p"$p"_d"$d".txt"
	echo $netFileName
	for l in "${lengths[@]}"
	do
	    nohup python MutNetAnalysis_networkx.py -d $netFileName -f False -l $l  >& kpathCents.p$p.d$d.l$l.oe &
	done
   done
done

