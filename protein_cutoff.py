#!/usr/bin/python3

import sys, os, re, argparse, textwrap
from collections import OrderedDict, defaultdict
from itertools import tee, islice, chain, zip_longest

currentdir = os.getcwd()
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='./%(prog)s [options] -o output_basename -p protein.fasta and/or -gf .gff3 and/or -gb .gbk',
    description = '''    A script to remove proteins from FASTA, GFF3, or GenBank
    format that are shorter than a provided threshold.''',
    
    epilog = """Written by Stephen A. Wyka (2019)""",
    formatter_class = MyFormatter)

parser.add_argument(
    '-o',
    '--out',
    required=True,
    help = 'Basename of output files',
    metavar=''
)
parser.add_argument(
    '-p',
    '--protein',
    help = 'protein file in FASTA format',
    metavar=''
)
parser.add_argument(
    '-gf',
    '--gff3',
    help = 'gff3 file',
    metavar=''
)
parser.add_argument(
    '-c',
    '--cutoff',
    type=int,
    default=50,
    help = 'Cutoff for protein length [default = 50]',
    metavar=''
)
args=parser.parse_args()

def get_next(some_iterable):
    items, nexts = tee(some_iterable, 2)
    nexts = chain(islice(nexts, 1, None), [None])
    return zip_longest(items, nexts)

# Create output folder and define paths to folders

if not os.path.isdir(os.path.join(currentdir, args.out+'_reformated')):
    os.makedirs(os.path.join(currentdir, args.out+'_reformated'))
result_dir = os.path.abspath(os.path.join(currentdir, args.out+'_reformated'))

# Checking arguments and making paths

if args.protein:
    protein_file = os.path.abspath(os.path.join(currentdir, args.protein))
    reformat_fasta = os.path.abspath(os.path.join(result_dir, args.out+'_reformat.fasta'))

if args.gff3:
    gff3_file = os.path.abspath(os.path.join(currentdir, args.gff3))
    reformat_gff = os.path.abspath(os.path.join(result_dir, args.out+'_reformat.gff3'))

# Removing proteins from FASTA file

if args.protein:
    print('Removing proteins shorter than {} aa in your FASTA file'.format(args.cutoff))
    protein_dict=defaultdict(str)
    try:
        with open(protein_file, 'r') as fasta_file:
            name=''
            for line in fasta_file:
                if line.startswith('>'):
                    if ' ' in line:
                        raise
                    name=line[1:-1]
                    continue
                protein_dict[name]+=line.strip()
    except:
        print("Error: The script does not like your FASTA file")
    with open(reformat_fasta, 'w') as re_fasta:
        for key, values in protein_dict.items():
            if len(values) >= args.cutoff:
                re_fasta.write('>'+key + '\n')
                re_fasta.write(textwrap.fill(''.join(protein_dict[key]),width=80) + '\n')

# Removing proteins from GFF3 file

if args.gff3:
    print('Removing proteins shorter than {} aa in your GFF3 file'.format(args.cutoff))
    line_list = []
    group_list = []
    gff_list = []
    with open(gff3_file, 'r') as gff_file:
        for line in gff_file:
            if line.startswith('#'):
                continue
            elif 'region' in line:
                continue
            elif 'gap' in line:
                continue
            elif 'pseudo=true' in line:
                continue
            elif 'tRNA' in line:
                continue
            else:
                gff_list.append(line.strip().split('\t'))
    gff_list = list(filter(None, gff_list))
    for line, next_line in get_next(gff_list):
        if next_line == None:
            line_list.append(line)
            continue
        elif line[2] == 'gene' and next_line[2] != 'gene':
            line_list.append(line)
        elif line[2] != 'gene' and next_line[2] == 'gene':
            line_list.append(line)
            group_list.append(line_list)
            line_list = []
        else:
            line_list.append(line)
    group_list.append(line_list)
    group_list = list(filter(None, group_list))
    cut_list = []
    for gene in group_list:
        lengths = []
        for feature in gene:
            if feature[2] == 'exon':
                lengths.append(int(feature[4]) - int(feature[3]))
        if sum(lengths) >= args.cutoff*3:
            cut_list.append(gene)
    with open(reformat_gff, 'w') as re_gff:
        re_gff.write('##gff-version 3' + '\n')
        for x in cut_list:
            for y in x:
                re_gff.write('\t'.join(y) + '\n')
