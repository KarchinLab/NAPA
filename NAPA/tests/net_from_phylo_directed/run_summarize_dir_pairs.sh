projDir=/home/vbeleva/phylogeny/pairs
workDir=$projDir/directedPairs
runDir=$workDir/results/summary
globalFasta=$runDir/PSE4-TEM-SHV_noLen.cds.prot.fasta

maxRuns=$1
maxTrees=$2

outPath=$workDir/results/summary/directedPairs/
mkdir -p $outPath

pvalThreshs=( 0.01 0.05 0.1 )
distThreshs=( 0.05 0.1 )

for p in "${pvalThreshs[@]}"
do
for d in "${distThreshs[@]}"
do

    nohup python summarize_directedPairs.py -d $runDir -r $maxRuns -m $maxTrees \
	-t $p -s $d -o $outPath >& \
	summarizeDirPairs.$maxRuns.$maxTrees.d$d.p$p.oe &

done
done
