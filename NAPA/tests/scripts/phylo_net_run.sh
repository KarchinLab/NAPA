projDir=$1
projDir=$projDir/NAPA
dataDir=$projDir/tests/data
resultsDir=$projDir/tests/results
scriptsDir=$projDir/tests/scripts

#Input files
alnFASTA=$dataDir/classA.protein.fasta
#fastaIDs=$dataDir/TEM-clin-ids.txt
posList=$dataDir/TEM-pos2aa.txt
posSubset=$dataDir/TEM-posSubset.txt
protFunc=$dataDir/TEM-functions.txt
funcTransitions=$dataDir/TEM-ext-spec-transitions.txt
selProtFunc=$dataDir/sel-funcs.txt
intFuncFileList=$dataDir/int-node-func-file-list.txt
treeFileList=$dataDir/tree-file-list.txt
intSeqFileList=$dataDir/int-anc-seq-file-list.txt
ls -lh $dataDir/run[1-2]/*internalStates.txt | awk '{print $9}' >  $intFuncFileList
ls -lh $dataDir/run[1-2]/*int_nwk.tree  |  awk '{print $9}'   >  $treeFileList
ls -lh $dataDir/run[1-2]/*anc.prot.fasta|  awk '{print $9}'   >  $intSeqFileList

#Input arguments
treeDistThresh=0.5
pvalThresh=0.01
#mutPairType='undir'
mutPairType='dir'
wtId='TEM_1'

#Output file
netFile=$resultsDir/TEM-phylo-net.$mutPairType.txt

# Construct the phylogeny based network
python $scriptsDir/phylo_mut_pairs.py \
    -a $alnFASTA \
    -pl $posList \
    -d $treeDistThresh \
    -wt $wtId \
    -ps $posSubset \
    -af $protFunc \
    -ft $funcTransitions \
    -sf $selProtFunc \
    -tl $treeFileList \
    -is $intSeqFileList \
    -if $intFuncFileList \
    -mt $mutPairType \
    -pt $pvalThresh \
    -nf $netFile \

