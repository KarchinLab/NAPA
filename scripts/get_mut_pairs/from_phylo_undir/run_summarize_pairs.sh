projDir=/home/vbeleva/phylogeny/pairs
workDir=$projDir/undirectedPairs
runDir=$workDir/results

maxRuns=$1
maxTrees=$2
#pvalThresh=$3

outPath=$runDir/summary/
mkdir -p $outPath

distThreshs=( 0.05 )
pvalThreshs=( 0.01 )
for d in "${distThreshs[@]}" 
do
for p in "${pvalThreshs[@]}"
do
    nohup python summarize_undirPairs.py -d $runDir -r $maxRuns -m $maxTrees \
    -t $p -dt $d -o $outPath >& summarizeUndirPairs.$maxRuns.$maxTrees.p$p.d$d.oe &
done
done
