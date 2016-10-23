#! /usr/bin/env python
import string
import re
from collections import defaultdict
from itertools import islice
'''
Python serial objects manipulation.
'''

#=========================================================#
# List/tuple manipulation
#=========================================================#
flatten = lambda x: [subitem for item in x \
                     for subitem in item]
list_pairs = lambda x: [(x[i1],x[i2]) \
                        for i1 in range(len(x)-1) \
                        for i2 in range(i1+1,len(x))]
sort_list_by_other = \
    lambda x,y: [x for (y,x) in \
                 sorted(zip(y,x), 
                        key=lambda pair: pair[0])]

def window(seq, n=2):
    '''
    Returns a sliding window (of width n) over data 
    from the iterable:
    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ... 
    '''
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result    
    for elem in it:
        result = result[1:] + (elem,)
        yield result

#=========================================================#
# Dictionaries
#=========================================================#
# Update dictionary with numerical values from  another
def num_dict_update_add(d0, d1):
    for k in d1:
        d0[k] += d1[k]
    return d0

#=========================================================#
# String manipulation
#=========================================================#
# Remove punctuation except for underscore
format_str = lambda s: \
    s.translate(string.maketrans("", ""), 
                string.punctuation.replace('_', ''))
#---------------------------------------------------------#
# Integer substring
get_int_substring = lambda s: int(re.sub('\D','',s))

# Get number from str

def get_int_from_str(s, default):
    try:
        return int(s)
    except ValueError:
        return default

def get_float_from_str(s, default):
    try:
        return float(s)
    except ValueError:
        return default

#---------------------------------------------------------#
# Lists of strings (sort by integer substring)
sort_by_digits = \
lambda x: sorted(x, key = lambda s: get_int_substring(s))

#=========================================================#
# Strings or numbers to Boolean 
#=========================================================#
def to_bool(v):
    if isinstance(v, bool):
        return v
    
    if isinstance(v, str):
        if v.lower() in [ '', '0', 'false', 'f', 'n', 'no']:
            return False
        else:
            return True

    if not len(v):
        return False
    return True
