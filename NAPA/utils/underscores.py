import re
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest = "infile")
    parser.add_argument('-o', dest = "outfile")
    args = parser.parse_args()
    return args

def underscore_convert(word):
    if len(word) < 3:
        return word
    if not word.startswith("\"") and not word.startswith("\'"):
        #s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
        #outname =  re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()
        word = re.sub(r"([A-Z]+)([A-Z][a-z])", r'\1_\2', word)
        word = re.sub(r"([a-z\d])([A-Z])", r'\1_\2', word)
        word = word.replace("-", "_")
        if not word.startswith('_'):
            word = word.replace('__','_')
        if len(word) >1:
            word = word[0]+''.join(word[1:]).lower()
        if '2' in word[1:-1]:
            word = word.replace('2','_to_')
    return word

def heal_line(line):
    outline = line.replace('=',' = ').replace('= =','==').replace('+ =','+=')
    outline = outline.replace('  =  ', ' = ')
    outline = outline.replace(',',', ').replace(',  ',', ')
    outline = outline.replace('+',' + ').replace('  +  ',' + ').replace('+ =','+=')

    if not outline.strip().startswith('class '):
        outline = ' '.join([underscore_convert(r) for r in outline.split(' ')])

    outline = outline.replace(' none',' None').replace(' false',' False').replace(' true',' True')

    return outline


def main():
    args = parse_arguments()
    with open(args.infile, 'r') as fi:
        with open(args.outfile, 'w') as fo:
            for line in fi:
                fo.write(heal_line(line))

if __name__ == '__main__': main()
