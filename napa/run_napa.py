#!/usr/bin/env python
'''
To construct network:
python run_napa.py -r build -c myconfig_file

Network analysis:
python run_napa.py -r analyze -c myconfig_file
'''

import argparse
from napa.utils.config import Config

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run', required = True,
                        help = ' '.join(['Command to run:',
                        '<build> or <analyze>',
                        'network']))

    parser.add_argument('-c', '--config_file', 
                        required = True,  
                        help = 'Path to configuration file.')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    config = Config(args.config_file, args.run)
    config.print_input_summary()
    
    if 'build' in args.run:
        if config.net_type == 'aln':
            from napa.analyze.aln_mut_pairs \
                import run_aln_mut_pairs
            run_aln_mut_pairs(config)
        elif 'phylo' in config.net_type: 
            from napa.analyze.phylo_mut_pairs \
                import run_phylo_mut_pairs
            run_phylo_mut_pairs(config)
        # create position-specific network
        elif config.net_type == 'aln-pos':
            from napa.analyze.aln_pos_pairs \
                import run_aln_pos_pairs
            run_aln_pos_pairs(config)
    
    if 'analysis' in args.run: 
        from napa.analyze.net_analysis \
            import run_net_analysis
        run_net_analysis(config)
	      

if __name__ == '__main__':
    main()

	
	      
