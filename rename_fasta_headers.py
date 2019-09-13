#!/usr/bin/python3

import os, sys, re, argparse
from collections import OrderedDict

current_dir = os.getcwd()

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='./%(prog)s [options] -i file.fasta -o output_basename -n prefix',
    description = '''    Renames fasta headers into NCBI locus_tag format (XXXXX_#####). 
    Returns new fasta file and a text file that matches old names to new names.''',
    
    epilog = """Written by Stephen A. Wyka (2019)""",
    formatter_class = MyFormatter)

parser.add_argument(
    '-i',
    '--input',
    required=True,
    help = 'Input FASTA file',
    metavar=''
)
parser.add_argument(
    '-o',
    '--out',
    required=True,
    help = 'Output basename',
    metavar=''
)
parser.add_argument(
    '-n',
    '--name',
    required=True,
    help = 'Prefix for FASTA headers [ex. Clpurp,Leacic]',
    metavar=''
)
args = parser.parse_args()

input_fasta = os.path.abspath(os.path.join(current_dir,args.input))
output_fasta = os.path.abspath(os.path.join(current_dir,args.out+'.fasta'))
output_list = os.path.abspath(os.path.join(current_dir,args.out+'_header_dictionary.txt'))

for file in os.listdir(current_dir):
    fasta_dict = OrderedDict()
    with open(input_fasta, 'r') as in_fasta, open(output_fasta, 'w') as out_fasta, open(output_list, 'w') as header_dict:
        count = 0
        for line in in_fasta:
            if line.startswith('>'):
                count += 1
                header = line.strip('>').strip()
                fasta_dict[header] = args.name+'_'+str(count)
                out_fasta.write('>'+args.name + '_' + str(count) + '\n')
            else:
                out_fasta.write(line)
        for key, value in fasta_dict.items():
            header_dict.write(key + ' = ' + value + '\n')
