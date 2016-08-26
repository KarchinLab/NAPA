import argparse
import os, sys

pwd = os.getcwd()
sys.path.append(pwd[0:pwd.find('/NAPA/')]) # NAPA package path


from napa.utils.general import *
from napa.net.net import *

def parseArgs():
    #Notes: perhaps better a yaml file with options for analysis
    parser = argparse.ArgumentParser( \
    description = 'Characterize mutation network and its nodes/links/clusters.')
    
    parser.add_argument('-nf', dest = "netFile", type=str, required=True,
                        help = 'Path to where weighted mutation'+\
                        ' pairs are stored.')

    parser.add_argument('-mt', dest = "mutNetType", type=str, required=True,
                        help = 'Type of mutation pair network: directed or undirected.')
    
    parser.add_argument('-pl', dest = "pathLength", type=int, default=2,
                        help = 'Length of central paths to consider, in number'+\
                        ' of nodes.')

    parser.add_argument('-co', dest = "centOutFile", type=str, 
                        default='net_centralities.txt',
                        help = 'Output file for centralities.')

    parser.add_argument('-ct', dest = "centralityType", type=str, 
                        default='net_path_centralities.txt',
                        help = 'node, path, or relative centrality.')


    args = parser.parse_args()
    return args

                
def main():
    args = parseArgs()
    mn = MutNet(net_file=args.netFile, net_type = args.mutNetType)
    
    cent_type = args.centralityType
    if cent_type == 'node':
        mn.get_node_centralities()
        with open(args.centOutFile, 'wb') as cf:
            cf.write(mn.str_node_centralities(header = True))
            
    elif cent_type == 'path':
        mn.get_path_between_path_cent(path_node_length = args.pathLength)
        mn.write_path_betw_path_cent(args.centOutFile)

    elif 'rel' in cent_type: 
        if cent_type == 'rel':
            # all relative centralities
            centralities =  mn.get_default_cent_list()

        elif cent_type == 'relloc':
            # local centrality types only: degree, clustering...
            centralities =  [c for c in mn.get_default_cent_list() if 'loc' in c]

        elif cent_type == 'relglob':
            # global centrality types only: closeness, eigenvalue, betweenness...
            centralities =  [c for c in mn.get_default_cent_list() if 'glob' in c]


        header = 'node.or.path' + (args.pathLength - 1) * '\t' + '\t'.join(centralities) + \
                 '\tnum.net.nodes\n'
        out_str = mn.get_rel_cent(cent_list = centralities, path_len = args.pathLength)

        with open(args.centOutFile, 'wb') as cf:
            cf.write(header + out_str)

if __name__ == "__main__": main()
