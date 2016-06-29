from cogent import LoadTree
import itertools, re 

if __name__ == "__main__":
    nwkFile = "run1.26560000.nwk.tree"
    
    # Pycogent still does not print internal node names!
    #tr = LoadTree(nwkFile)
    #print tr.asciiArt()
    #tr.writeToFile("run1.26560000.internal.nwk.tree",with_distances=True)

    # use regexp to infer internal node names sequentially (as in PyCogent)
    with open(nwkFile, 'rb') as nF: nFstring = nF.read()
    cnt = itertools.count()
    outNwkStr = re.sub(r'\)\:', lambda x: ')edge_{}:'.format(next(cnt)), nFstring)
    outNwkStr = re.sub(r'(.*)edge_(.*)\:',r'\1root:',outNwkStr)
    
    with open(nwkFile.replace('nwk.tree','int_nwk.tree'), 'wb') as nOut:
        nOut.write(outNwkStr)

# treeNoInternal = "((D:0.723274,F:0.567784):0.067192,(B:0.279326,H:0.756049):0.807788);"
# print re.sub(r'\)\:', lambda x: ')edge_{}:'.format(next(cnt)), treeNoInternal) 
