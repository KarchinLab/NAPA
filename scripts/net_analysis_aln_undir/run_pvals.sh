projDir=/home/vbeleva/phylogeny/pairs
workDir=$projDir/epistasisPairs
summDir=$workDir/results/summary/
randNetDir=$summDir/randomNets

outPath=$summDir/pvals
mkdir -p $outPath

numSwaps=$1
pathLen=$2

nohup python get_randpvals.py -r$randNetDir -n$numSwaps -p$pathLen >& randpvals.$numSwaps.$pathLen.oe &
