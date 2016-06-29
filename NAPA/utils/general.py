#! /usr/bin/env python 
import sys
import re

# List manipulation
flatten = lambda x: [subitem for item in x for subitem in item]
list_pairs = lambda x: [(x[i1],x[i2]) for i1 in range(len(x)-1) \
                       for i2 in range(i1+1,len(x))]
#integer substring
get_int = lambda s: int(re.sub('\D','',s))

#lists of strings
sort_by_digits = lambda x: sorted(x, key = lambda s: getInt(s))

#list sorting
sort_list_by_other = lambda x,y: [x for (y,x) in sorted(zip(Y,X),
                                        key=lambda pair: pair[0])]


# General IO
def stderr_write(myList):
    sys.stderr.write(' '.join([str(el) for el in myList])+'\n')


def read_keyval_dict(infile):
    ''' Simple {key:val} dict from first two columns 
    in file '''

    out_dict = {}
    with open(infile, 'rb'):
        for line in pf.read():
            recs = line.strip().split().replace('"', '')
            if len(recs)>1:
                out_dict[recs[0]] = recs[1]
    return out_dict
