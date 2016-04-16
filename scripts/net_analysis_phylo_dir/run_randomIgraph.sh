projDir=/home/vbeleva/phylogeny/pairs
workDir=$projDir/directedPairs
summDir=$workDir/results/summary/directedPairs

outPath=$summDir/randomNets
mkdir -p $outPath

pvalThresh=$1
numNets=$2
numSwaps=$3

netFileName=$summDir/2015-07-14.directPhenoPair_$pvalThresh.txt 

nohup python randomizedIgraph.py -d $netFileName -e $numSwaps -n $numNets >& randIgraph.p$pvalThresh.n$numNets.s$numSwaps.oe &



