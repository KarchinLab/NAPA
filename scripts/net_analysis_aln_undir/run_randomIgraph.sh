projDir=/home/vbeleva/phylogeny/pairs
workDir=$projDir/undirectedPairs
summDir=$workDir/results/summary

outPath=$summDir/randomNets
mkdir -p $outPath

numNets=$1
numSwaps=$2


distThreshs=( 0.02 0.05 )
pvalThreshs=( 0.01 0.05 )
for d in "${distThreshs[@]}" 
do
for p in "${pvalThreshs[@]}"
do
    netFileName=$summDir/2015-09-11.undirPhenoPair$p.dist$d.txt
    echo "python randomizedIgraph.py -f True -d $netFileName -e $numSwaps -n $numNets"
    nohup python randomizedIgraph.py -f True -d $netFileName -e $numSwaps -n $numNets >& \
	randIgraph.p$p.d$d.n$numNets.s$numSwaps.oe &
done
done


