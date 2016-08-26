# Not tested on test set
# Renames internal node names from newick tree files
# and alignment files to be consistent with each other

def consistent_internal_node_names(nwkFile,fastaFile, gfStr,
                                   nwkOutFile, fastaOutFile):

    # change formatting of edge ids output from 
    # cogent ancestral reconstruction
    with open(fastaFile, 'rb') as fF: fString = fF.read()
    fStringOut = fString.replace('>edge.','>edge')

    # use regexp to infer internal node names sequentially (as in PyCogent)
    with open(nwkFile, 'rb') as nF: nFstring = nF.read()
    cnt = itertools.count()
    outNwkStr = re.sub(r'\)\:', lambda x: ')edge{}:'.format(next(cnt)), nFstring)
    outNwkStr = re.sub(r'(.*)edge(.*)\:',r'\1root:',outNwkStr)
    
    with open(fastaOutFile, 'wb') as fOut:
        fOut.write(fStringOut)

    with open(nwkOutFile, 'wb') as nOut:
        nOut.write(outNwkStr)

    return fastaOutFile, nwkOutFile
