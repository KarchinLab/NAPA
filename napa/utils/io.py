#! /usr/bin/env python
import sys
import glob
from collections import defaultdict

from napa.utils.serials import *

'''
Python file objects manipulation and IO.
'''

#=========================================================#
# Standard input/ output/error
#=========================================================#
def stderr_write(inlist):
    sys.stderr.write(\
        ' '.join([str(v) for v in inlist]) + '\n')

#=========================================================#
# Listing files in a directory
#=========================================================#
def list_files(path_before_wildcard = '', 
               path_after_wildcard = '', 
               sort_bool = True):
    if not sort_bool:
        file_list = \
        list(glob.glob(''.join([path_before_wildcard, '*', 
                             path_after_wildcard])))
    else:
        file_list = \
        sorted(list(glob.glob(''.join([path_before_wildcard, 
                            '*', path_after_wildcard]))))
    return file_list

def list_display_files(file_type = '', prefix = '', 
                       suffix = ''):
    stderr_write(['Listing', str(file_type), 'files',
                  prefix, '...', suffix])
    return list_files(prefix, suffix)


#=========================================================#
# General parsing
#=========================================================#
def file_line_iterator(file_name):
    with open(file_name, 'rb') as f:
        for line in f: 
            yield line.strip()

def file_line_list(file_name):
    with open(file_name, 'rb') as f:
        return f.readlines()

def write_table_str(file_name, header, body):
    with open(file_name, 'w') as f:
        f.write(header + body)

#=========================================================#
# Parsing into python serial objects
#=========================================================#
def parse_column(infile, col_num = 0):
    out_list = []
    with open(infile, 'r') as f:
        for line in f:
            recs = line.replace('"','').strip().split()
            if len(recs) > col_num:
                out_list.append(recs[col_num])
    return out_list

def parse_keyval_dict(infile):
    ''' 
    Simple {key:val} dict from first two columns 
    in file 
    '''
    out_dict = {}
    with open(infile, 'rb') as f:
        for line in f:
            recs = [format_str(s) for s in \
                    line.strip().split()]
            if len(recs) > 1:
                out_dict[recs[0]] = recs[1]
    return out_dict

def parse_keyval_dlist(infile):
    ''' 
    List-valued {key:[val0,val1,...]} dict from 
    first two columns in file 
    '''
    out_dict = defaultdict(list)
    with open(infile, 'rb') as f:
        for line in f:
            recs = [format_str(s) for s in \
                    line.strip().split()]
            if len(recs) > 1:
                out_dict[recs[0]].append(recs[1])
    return out_dict


#=========================================================#
# Writing files from python serial objects
#=========================================================#
def write_list_column(out_list = [],
                 outfile = 'list-column.txt', header = ''):
    with open(outfile, 'w') as f:
        if len(header):
            f.write(header + '\n')
        f.write('\n'.join([str(v) for v in out_list]) + '\n')
            

def write_keyval_dict(out_dict = {},
                      outfile = 'dict-keyval.txt', 
                      header = ''):
    with open(outfile, 'w') as f:
        if len(header):
            f.write(header + '\n')
        f.write('\n'.join(['%s\t%s'%(k,v) \
                for k,v in out_dict.iteritems()]) + '\n')

def write_keyval_dlist(out_dict = {},
                       outfile = 'dict-keyvallist.txt', 
                       header = '', val_sep = ';'):

    with open(outfile, 'w') as f:
        if len(header):
            f.write(header + '\n')
        f.write('\n'.join(\
            ['%s\t%s'%(k, 
                       val_sep.join([str(v) for v in l])) \
                       for k, l in out_dict.iteritems()]))
        f.write('\n')
