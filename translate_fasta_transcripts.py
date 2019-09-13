#!/usr/bin/python3

import os, sys, argparse, textwrap, re
from collections import defaultdict

rundir = os.getcwd()

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='./%(prog)s [options] -i orthogroups -o output_folder',
    description = '''    Translates a fasta file of transcripts (nucleotides) into 
    proteins. This will be a direct translation, make sure initial fasta file contains 
    only CDS of the gene (no introns).''',
    
    epilog = """Written by Stephen A. Wyka (2019)""",
    formatter_class = MyFormatter)

parser.add_argument(
    '-i',
    '--input',
    required = True,
    help = 'Input fasta file to alter',
    metavar=''
)
parser.add_argument(
    '-o',
    '--output',
    required = True,
    help = 'Output protein fasta file',
    metavar=''
)
args=parser.parse_args()

input_fasta = os.path.abspath(args.input)
output_fasta = os.path.abspath(os.path.join(rundir, args.output))

def translate(seq, key, dict): 
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'', 
        'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W', 
        'NNN':'X'
    } 
    protein ="" 
    if len(seq)%3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] 
            if 'N' in codon or 'n' in codon:
                print("N's were found in your sequence, these codons were translated",\
                "as X.")
                codon = 'NNN'
            else:
                protein+= table[codon] 
        dict[key] = protein
    else:
        print('Some transcripts were not divisible by 3, thus were not translated')
    return

fasta_dict = defaultdict(str)
with open(input_fasta, 'r') as in_fasta:
    name = ''
    for line in in_fasta:
        if line.startswith('>'):
                name = line[1:-1]
                continue
        fasta_dict[name]+=line.strip()

new_fasta = {}
for k, v in fasta_dict.items():
    translate(v, k, new_fasta)

if len(new_fasta) == len(fasta_dict):
    None
else:
    print('Something went wrong. The number of genes in new fasta file did not ',\
    'equal the number of genes in the old fasta file')

# try to sort new fasta in best possible way
if re.search(r'(\w+_\d+)', next(iter(new_fasta))) and ' ' not in next(iter(new_fasta)): # locus tag format
    if '-' in next(iter(new_fasta)): # locus_tag with -mRNA potentially
        sort_new = sorted(new_fasta.items(), key=lambda k: int(k[0].split('_')[1].split('-')[0]))
    else:
        sort_new = sorted(new_fasta.items(), key=lambda k: int(k[0].split('_')[1]))
elif re.search(r'(\w+\|\w+\|\w+)', next(iter(new_fasta))): # NCBI/ENA format
    sort_new = sorted(new_fasta.items(), key=lambda k: k[0].split('|')[1])
else:
    sort_new = sorted(new_fasta.items()) # last resort

with open(output_fasta, 'w') as out_fasta:
    for k,v in sort_new:
        out_fasta.write('>' + k + '\n')
        out_fasta.write(textwrap.fill(''.join(v),width=80) + '\n')