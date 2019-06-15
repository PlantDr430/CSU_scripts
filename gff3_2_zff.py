#!/usr/bin/python3

"""
Written by Stephen A. Wyka (2019)

A script to turn a gff3 file into a zff file to be used to train SNAP for gene 
prediction. This script will take a gff3 file along with a genome file and output a 
genome.ann (zff) genome.dna (altered genome file) file to be used for SNAP gene predicton
 training.

Make sure your FASTA headers in your genome file match EXACTLY with your column 1 of your
gff3 file. Extra data in your FASTA headers beside the contig/scaffold ID will upset the 
script.


"""

import os, sys, re, argparse, inspect, shutil, textwrap
from collections import OrderedDict, defaultdict
rundir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(rundir)
sys.path.insert(0, parentdir)

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='python3 %(prog)s [options] -in gff3 -g genome.fasta -o output_basename',
    description = '''    A script to turn a gff3 file into a zff file to be used to 
    train SNAP for gene prediction.''',
    
    epilog = """Written by Stephen A. Wyka (2019)""",
    formatter_class = MyFormatter)

parser.add_argument(
    '-in',
    '--input',
    required=True,
    help = 'GFF3 file',
    metavar=''
)
parser.add_argument(
    '-o',
    '--out',
    required=True,
    help = 'Basename of output',
    metavar=''
)
parser.add_argument(
    '-g',
    '--genome',
    required=True,
    help = 'Genome in FASTA format',
    metavar=''
)
args=parser.parse_args()

# Making directory

if not os.path.isdir(args.out+'_SNAP'):
    os.makedirs(os.path.join(args.out+'_SNAP'))

# Create file paths

zff_dir = os.path.abspath(os.path.join(args.out+'_SNAP'))
snap_zff = os.path.abspath(os.path.join(zff_dir, args.out+'_genome.ann'))
snap_dna = os.path.abspath(os.path.join(zff_dir, args.out+'_genome.dna'))

# Checking arguments

if args.input:
    input_gff = os.path.abspath(os.path.join(rundir, args.input))
else:
    print('Error: Please provide a gff3 file, -in or --input')

if args.genome:
    input_genome = os.path.abspath(os.path.join(rundir, args.genome))
else:
    print('Error: Please provide an associated genome.fasta file, -g or --genome')

exon_gff = os.path.abspath(os.path.join(zff_dir, 'exon.gff'))

def create_exon_gff(input_gff, output_tmp):
    with open(input_gff, 'r') as in_gff, open(exon_gff, 'w') as exon_tmp:
        for line in in_gff:
            if not line.startswith('#'):
                sect = line.split('\t')
                if sect[2] == 'exon':
                        parent = re.search('(Parent=)(.+?)(;)', sect[8]).group(2)
                        exon_tmp.write(line.replace(sect[8], parent) + '\n')

create_exon_gff(input_gff, exon_gff)

if os.stat(exon_gff).st_size == 0:
    print("Error: The script does not like your gff3 file", \
    "Please make sure your gff has 'exon' in column 3 and 'Parent=' in column 9")
    sys.exit(1)

# Turn genome into a dictionary
print('Creating genome dictionary and altering gff3 file')
genome_dict=defaultdict(str)
try:
    with open(input_genome, 'r') as genome:
        name=''
        for line in genome:
            if line.startswith('>'):
                if ' ' in line:
                    raise
                name=line[1:-1]
                continue
            genome_dict[name]+=line.strip()
except:
    print("Error: The script does not like your genome file", \
    "Please make sure there are no spaces in your FASTA header and your headers",\
    "match column 1 of your gff3 file exactly.")


snap_list = []
snap_2nd = []
group_list = []
lines = []
gene_dict = {}
contig_dict = OrderedDict()
with open(exon_gff, 'r') as exons:
    lines = exons.read().splitlines() 
    for line in lines:
        contig = line.split('\t')[0]
        gene = line.split('\t')[8]
        if gene not in gene_dict.keys():
            gene_dict[gene] = []
        if contig not in contig_dict.keys():
            contig_dict[contig] = []
    for line in lines:
        gene = line.split('\t')[8]
        if gene in gene_dict.keys():
            gene_dict[gene].append(line.split('\t'))
for values in gene_dict.values():
    group_list.append(values)
snap_tmp = []
for x in group_list:
    for y in x:
        if '-' in y:
            temp = y.pop(3)
            y.insert(4,temp)
            snap_tmp.append(y)
        else:
            snap_tmp.append(y)
    snap_list.append(snap_tmp)
    snap_tmp = []
snap_list.sort(key = lambda x: (x[0][0], -(int(x[0][3]))), reverse=True)
snap_tmp = []
for x in snap_list:
    if len(x) == 1:
        for y in x:
            y = [x.replace('exon','Esngl') for x in y]
            snap_tmp.append(y)
    if len(x) == 2:
        for y in x:
            if y == x[0]:
                y = [x.replace('exon','Einit') for x in y]
                snap_tmp.append(y)
            else:
                y = [x.replace('exon','Eterm') for x in y]
                snap_tmp.append(y)
    elif len(x) > 2:
        for y in x:
            if y == x[0]:
                y = [x.replace('exon','Einit') for x in y]
                snap_tmp.append(y)
            elif y == x[(len(x) - 1)]:
                y = [x.replace('exon','Eterm') for x in y]
                snap_tmp.append(y)
            else:
                y = [x.replace('exon','Exon') for x in y]
                snap_tmp.append(y)
    snap_2nd.append(snap_tmp)
    snap_tmp = []
for x in snap_2nd:
    for y in x:
        if y[0] in contig_dict.keys():
            contig_dict[y[0]].append([y[2],y[3],y[4],y[8]])

print("Writing output files. Won't be long now.")
with open(snap_zff, 'w') as genome_ann, open(snap_dna, 'w') as genome_dna:
    for cont in contig_dict.keys():
        genome_ann.write('>' + cont + '\n')
        for x in contig_dict[cont]:
            genome_ann.write('\t'.join(x) + '\n')
        if cont in genome_dict.keys():
            genome_dna.write('>' + cont + '\n')
            genome_dna.write(textwrap.fill(''.join(genome_dict[cont]),width=80) + '\n')

os.remove(exon_gff)
