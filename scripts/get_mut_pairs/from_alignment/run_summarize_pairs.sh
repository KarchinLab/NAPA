projDir=/home/vbeleva/phylogeny/pairs
workDir=$projDir/undirectedPairs
runDir=$workDir/results

maxRuns=$1
maxTrees=$2
pvalThresh=$3

outPath=$runDir/summary/
mkdir -p $outPath

distThreshs=( 0.02 0.05 0.1 0.2 0.5 )
pvalThreshs=( 0.01 0.05 0.1 )
for d in "${distThreshs[@]}" 
do
for p in "${pvalThreshs[@]}"
do
    nohup python summarize_undirPairs.py -d $runDir -r $maxRuns -m $maxTrees \
    -t $p -dt $d -o $outPath >& summarizeUndirPairs.$maxRuns.$maxTrees.p$p.d$d.oe &
done
done
