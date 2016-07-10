import argparse

from NAPA.utils.general import *
from NAPA.net.net import *

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

    parser.add_argument('-pc', dest = "pathCentOutFile", type=str, 
                        default='net_path_centralities.txt',
                        help = 'Output file for betweeness centralities of paths.')

    args = parser.parse_args()
    return args

                
def main():
    args = parseArgs()
    mg = MutNet(net_file=args.netFile, net_type = args.mutNetType)
    mg.get_path_between_path_cent(path_node_length = 3, cent_list = ['shortest_path',
                                                                     'k_path'])
    mg.write_path_betw_path_cent(args.pathCentOutFile)


if __name__ == "__main__": main()
