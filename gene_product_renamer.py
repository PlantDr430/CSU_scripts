#!/usr/bin/python3

'''
Written by Stephen A. Wyka (2019)

This script was written for use after a run of the annotate function from Funannotate. I noticed that
many of the proteins were still classified as 'hypothetical proteins' even after functional annotation from 
multiple programs. To help give product names to genes I decided to use the *.proteins.fasta and *.annotations.txt
files from the Funannotate annotate_results/ directory to blast against a set of proteins from a reference annotated 
genome(s) and pull gene product names from genes that had a >80% identify and >90% alignment length to my predicted 
genes. This creates a 3 tab deliminated custom annotation file which can then be used with Funannotate annotate with 
the '-a' flag.

'''

import os, sys, re, argparse

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='./%(prog)s [options] -qu query_prot.fasta -an *.annotations.txt -db database_prot.fasta -o output_dir',
    description = '''    Script to create a 3 tab custom annotation file of gene prdocut names for Funannotate. 
    Database protein FASTA file needs to be taken from NCBI/GenBank or have the following format of FASTA headers

    >Protein_ID product [species]

    For example:
    >CCE35419.1 probable amino acid transport protein GAP1 [Claviceps purpurea 20.1]''',
    epilog = """Written by Stephen A. Wyka (2019)""",
    formatter_class = MyFormatter)

parser.add_argument(
    '-qu',
    '--query',
    required=True,
    help = 'Query file in FASTA format',
    metavar=''
)
parser.add_argument(
    '-an',
    '--annotations',
    required=True,
    help = '*.annotations.txt file from Funannotate annotate_results/ output',
    metavar=''
)
parser.add_argument(
    '-db',
    '--database',
    required=True,
    help = 'FASTA file obtained from NCBI/GenBank to be used for database',
    metavar=''
)
parser.add_argument(
    '-o',
    '--out',
    required=True,
    help = 'Basename of output directory and file',
    metavar=''
)
parser.add_argument(
    '--cpus',
    default=1,
    type=int,
    help = 'Number of cores to use in BLAST [default: 1]',
    metavar=''
)
parser.add_argument(
    '-id',
    '--identity',
    default=0.90,
    type=float,
    help = 'Cutoff value for percent identity [default: 0.80]',
    metavar=''
)
parser.add_argument(
    '-L',
    '--align_length',
    default=0.90,
    type=float,
    help = 'Cutoff value for percent alignment length / query cover [default: 0.75]',
    metavar=''
)
parser.add_argument(
    '--MAKEBLASTDB_PATH',
    help = 'Path to makeblastdb exe if not set to root',
    metavar=''
)
parser.add_argument(
    '--BLASTP_PATH',
    help = 'Path to blastp exe if not set to root',
    metavar=''
)
args=parser.parse_args()

def which_path(file_name):
    for path in os.environ["PATH"].split(os.pathsep):
        full_path = os.path.join(path, file_name)
        if os.path.exists(full_path) and os.access(full_path, os.X_OK):
            return full_path
    return None

if args.query:
    query_file = args.query
else:
    print('Error: Please provide a query file, -q or --query')

# Check to see if query file is protein sequences
try:
    with open(args.query) as que:
            PROTEIN = 'DdEeFfHhIiKkLlMmPpQqRrSsVvWwYy'
            for line in que:
                if not line.startswith('>'):
                    if any(i in PROTEIN for i in line):
                        None
                    else:
                        raise()
                    break
except:
    print('Invalid query file, please provide a FASTA file of protein sequences')
    sys.exit(1)

if args.annotations:
    annotations_file = args.annotations
else:
    print('Error: Please provide an annotation file, -a or --annotations')

if args.database:
    database_name = os.path.splitext(os.path.basename(args.database))[0]
else:
    print('Error: Please provide a database FASTA file')

# Check to see if database file is protein sequences
try:
    with open(args.database) as ref:
        PROTEIN = 'DdEeFfHhIiKkLlMmPpQqRrSsVvWwYy'
        for line in ref:
            if not line.startswith('>'):
                if any(i in PROTEIN for i in line):
                    DB_TYPE = 'prot'
                else:
                    raise()
                break
except:
    print('Invalid database file, please provide a FASTA file of protein sequences')
    sys.exit(1)

# Check FASTA header format of database file

try:
    with open(args.database) as ref:
        for line in ref:
            if line.startswith('>'):
                if not re.search(r'>(\S+) (.+) (\[.+])', line):
                    raise()
                else:
                    None
                break
except:
    print('Invalid FASTA header format of database file, see help message')
    sys.exit(1)

# create output folders and define paths to folders
if not os.path.isdir(args.out):
    os.makedirs(os.path.join(args.out, 'blast_db'))
    os.makedirs(os.path.join(args.out, 'misc_files'))
blast_dir = os.path.join(args.out, 'blast_db')
misc_dir = os.path.join(args.out, 'misc_files')

# fidning paths to required programs
if args.MAKEBLASTDB_PATH:
    MAKEBLASTDB = args.MAKEBLASTDB_PATH
else:
    try:
        if which_path('makeblastdb'):
            MAKEBLASTDB = 'makeblastdb'
    except:
        print('makeblastdb not found, please provide path to executable in command with --MAKEBLASTDB_PATH')

if args.BLASTP_PATH:
    BLASTP = args.BLASTP_PATH
else:
    try:
        if which_path('blastp'):
            BLASTP = 'blastp'
    except:
        print('blastp not found, please provide path to executable in command with --BLASTP_PATH')

### starting run ###
# creating blast database
print('Creating blast database')
database_path = os.path.abspath(os.path.join(blast_dir, database_name))
MakeDB=[MAKEBLASTDB, ' -in ', args.database, ' -dbtype ', DB_TYPE, ' -parse_seqids ', '-out ', database_path]
os.system(''.join(MakeDB))

# run blastp -outfmt 6 qseqid sseqid pident length slen
print('Running blastp')
blast_output = os.path.abspath(os.path.join(misc_dir, 'blast_out.txt'))
run_blast = [BLASTP, ' -db ', database_path, ' -query ', args.query, ' -out ', blast_output, 
    " -evalue 0.000001 -outfmt '6 qseqid sseqid pident length slen' -max_target_seqs 1 -max_hsps 1", ' -num_threads ', 
    str(args.cpus)]
os.system(''.join(run_blast))

# filtering blast results based on percent identify and alignment length
# creating dictionary of query genes and reference genes
filtered_blast = {}
with open(os.path.join(misc_dir, 'blast_out.txt'), 'r') as blast_out:
    for line in blast_out:
        columns = line.split('\t')
        if "|" in columns[1]:
            columns[1] = re.search(r'\|(\S+)\|', columns[1]).group(1)
        else:
            continue
        alignment_length = float((int(columns[3]) / int(columns[4])))
        if float(columns[2]) >= (args.identity * 100) and alignment_length >= args.align_length:
            filtered_blast[columns[0].strip()] = columns[1].strip()

# getting only hypothetical or uncharacterized genes from Funannotate annotation output
hypo_genes = []
with open(args.annotations, 'r') as anno:
    for line in anno:
        columns = line.split('\t')
        if columns[7] == 'hypothetical protein' or columns[7] == 'uncharacterized protein':
            hypo_genes.append(columns[0].strip() + '-T1')

# getting gene product names from reference file
reference_dict = {}
with open(args.database, 'r') as reference:
    for line in reference:
        if line.startswith('>'):
            ref_gene = re.search(r'>(\S+) (.+) (\[.+])', line).group(1)
            ref_product = re.search(r'>(\S+) (.+) (\[.+])', line).group(2)
            if 'hypothetical protein' in ref_product or 'uncharacterized protein' in ref_product or \
            'putative protein' in ref_product:
                continue
            else:
                reference_dict[ref_gene] = ref_product

# creating final list of query genes with potential gene products
filtered_genes = set(list(filtered_blast.keys())) & set(hypo_genes)
reference_genes = reference_dict.keys()
final = {}
for g in filtered_genes:
    if filtered_blast[g] in reference_genes:
        final[g] = reference_dict[filtered_blast[g]]


# creating 3 tab deliminted custom annotation output file
with open(os.path.join(args.out,'{}_custom_annotations.txt'.format(args.out)), 'w') as outfile:
    for gene, product in final.items():
        outfile.write(gene + '\t' + 'product' + '\t' + product + '\n')