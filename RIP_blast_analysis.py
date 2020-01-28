#!/usr/bin/python3

'''
This script uses all v. all blast searches of proteins or transcripts of a geneome 
as a method to identify the presence of RIP in the genome. This method has been used 
in Galagan et al. (2003) and Urquhart et al. (2018). The idea behind the method is that 
RIP requires about 400 bp with a greater than ~80% sequence identity between duplicated 
regions to activate. Therefore, if a genome shows absence of protein or nucleotide sequences 
> 80% that can indicate that RIP is present and actively trying to silince the repeatitive DNA. 
Where as if a genome has sequences with identities of > 80% to other sequences in the genome then 
that could indicate that RIP is not present in the genome.

For further reading please see the two paper that use this method.

Galagan, J.E., Clavo, S. E., Borkovich, K.A., et al. 2003. The genome sequence of the filamentous 
fungus Neurospora crassa. Nature. 422:859-868

Urquhart, A.S., Mondo, S.J., Makela, M.R., et al. 2013. Genomic and genetic insights into a 
cosmopolitan fungus, Paecilomyces variotii (Eurotiales). Frontiers in Microbiology. 
https://doi.org/10.3389/fmicb.2018.03058


As for script usage:

The script can take a single protein or transcript fasta file and perform blast (p or n or x) on 
in an all v. all method. It will then remove self-hit and bin the percent identities by rounding 
to the nearest number. The bins and fraction of genes in each % bin to the number of orthologs 
(best hits) found in the blast search are plotted on the x and y-axis, respectively on a 2D bar chart. 

This script can also accept a batch file for multiple species and will plot them all together in 
a 3D plot. The batch file is a tab-deliminted file and should have a structure such as this: (NOTE: 
the 3rd column (Identifier is not needed) this column will be used to group these species by the same 
color on the plot. If identifiers are not provided then each species will have a separate color.) Also,
the 1st column (Genus species) will be used for the tick labels on the Z-axis of the plot.

Genus species   /full/path/to/fasta/file    Identifier for grouping


Example:
Species_1	/home/directory/Species_1_proteins.fasta	Clade_1
Species_2	/home/directory/Species_2_proteins.fasta	Clade_2
Species_3	/home/directory/Species_3_proteins.fasta	Clade_2
Species_4	/home/directory/Species_4_proteins.fasta	Clade_2
Species_5	/home/directory/Species_5_proteins.fasta	Clade_3
Species_6	/home/directory/Species_6_proteins.fasta	Clade_3

'''

import os, sys, argparse, inspect, subprocess, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
from collections import OrderedDict

rundir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(rundir)
currentdir = os.getcwd()
sys.path.insert(0, parentdir)

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='./%(prog)s [options] [-in fasta_file -s species_name] or [-b batch_file] -o output_directory',
    description = '''    Performs all v. all blast analysis and bins 2nd best hits by identity 
    (i.e. next hit besides itself) against the fraction of genes in bins of total orthologs.''',
    
    epilog = """Written by Stephen A. Wyka (2019)""",
    formatter_class = MyFormatter)

parser.add_argument(
    '-in',
    '--input',
    help = 'Input FASTA file (aa or nt)',
    metavar=''
)
parser.add_argument(
    '-o',
    '--out',
    required=True,
    help = 'Output directory name',
    metavar=''
)
parser.add_argument(
    '-s',
    '--species',
    help = '',
    metavar='Name of species (use quotes if using spaces)'
)
parser.add_argument(
    '-b',
    '--batch',
    help = '',
    metavar='Tab-delimited batch file containing species name and fasta file paths'
)
parser.add_argument(
    '-t',
    '--dbtype',
    choices= ['prot', 'nucl'],
    default='prot',
    help = 'dbtype of your FASTA file [prot|nucl] [default: prot]',
    metavar=''
)
parser.add_argument(
    '-bl',
    '--blast_type',
    choices= ['blastp', 'blastn', 'blastx'],
    default='blastp',
    help = 'Which blast to run [blastp|blastn|blastx] [default: blastp]',
    metavar=''
)
parser.add_argument(
    '-c',
    '--cpus',
    default=2,
    type=int,
    help = 'Number of cores to use [default: 2]',
    metavar=''
)
parser.add_argument(
    '-e',
    '--evalue',
    default=1e-5,
    type=float,
    help = 'Cutoff e_value for blastP[default: 1e-5]',
    metavar=''
)
parser.add_argument(
    '--legend',
    default = 8,
    type=int,
    help = 'Size of legend for overlaid genes [default = 8]',
    metavar=''
)
parser.add_argument(
    '--figsize',
    nargs = '+',
    default = [15, 6],
    type=float,
    help = 'Dimensions of figure in inches [default = 15 5]',
    metavar=''
)
parser.add_argument(
    '--dpi',
    default = 150,
    type=int,
    help = 'DPI of figure [default = 150]',
    metavar=''
)
parser.add_argument(
    '--xlimits',
    nargs = '+',
    default = [15, 101],
    type=float,
    help = 'Min and max limits of the X-axis in log scale [default = 15 101]',
    metavar=''
)
parser.add_argument(
    '--elev',
    default = 50,
    type=int,
    help = 'Elevation angle to view 3D plot [default = 50]',
    metavar=''
)
parser.add_argument(
    '--azim',
    default = -50,
    type=int,
    help = 'Angle of z-axis to view 3D plot [default = -50]',
    metavar=''
)
parser.add_argument(
    '--yfontsize',
    default = 12,
    type=int,
    help = 'Fontsize of y-axis labels [default = 12]',
    metavar=''
)
parser.add_argument(
    '--yrotation',
    default = -20,
    type=int,
    help = 'Rotation of y-axis labels [default = -20]',
    metavar=''
)
args=parser.parse_args()

# Check for dependencies
def which_path(file_name):
    for path in os.environ["PATH"].split(os.pathsep):
        full_path = os.path.join(path, file_name)
        if os.path.exists(full_path) and os.access(full_path, os.X_OK):
            return full_path
    return None

try:
    if which_path('makeblastdb'):
        MAKEBLASTDB = 'makeblastdb'
    else:
        raise
except:
    print('ERROR: makeblastdb not found, please make sure parent directory of', \
    'makeblastdb is located in $PATH')
    sys.exit()

if args.blast_type == 'blastp':
    try:
        if which_path('blastp'):
            BLAST = 'blastp'
        else:
            raise
    except:
        print('ERROR: blastp not found, please make sure parent directory of', \
        'blastp is located in $PATH')
        sys.exit()

if args.blast_type == 'blastn':
    try:
        if which_path('blastn'):
            BLAST = 'blastn'
        else:
            raise
    except:
        print('ERROR: blastn not found, please make sure parent directory of', \
        'blastn is located in $PATH')
        sys.exit()

if args.blast_type == 'blastx':
    try:
        if which_path('blastx'):
            BLAST = 'blastx'
        else:
            raise
    except:
        print('ERROR: blastx not found, please make sure parent directory of', \
        'blastx is located in $PATH')
        sys.exit()

# check arguments are correct
if args.input:
    if args.batch:
        parser.error('please provide only an input fasta file or batch file, not both')
    if not args.species:
        parser.error('-s / --species is reuired if input file is provided')
    if ' ' in args.species:
        name = args.species.replace(' ', '_')
    else:
        name = args.species
    input_fasta = os.path.abspath(args.input)

if args.batch:
    input_tuples = []
    categories = []
    batch_file =  os.path.abspath(os.path.join(currentdir, args.batch))
    with open(batch_file ,'r') as batch_in:
        batch_linearr = [line.strip().split('\t') for line in batch_in if not line.startswith('#')]
        for line in batch_linearr:
            if len(batch_linearr[0]) == 3: # if there are categories
                categories.append(line[2])
                input_tuples.append((line[0].replace(' ', '_'),os.path.join(currentdir,line[1]),line[2]))
            else: # if there are no categories
                input_tuples.append((line[0].replace(' ', '_'),os.path.join(currentdir,line[1])))
    if categories:
        input_tuples = input_tuples
        # input_tuples = sorted(input_tuples, key = lambda x: (x[2], x[0])) # sort to group categories
        categories = list(set(categories)) # get set list of all categories

# create output folders and define pathsep
if not os.path.isdir(args.out):
    os.makedirs(os.path.join(args.out, 'databases'))
    os.makedirs(os.path.join(args.out, 'blast_dir'))
    os.makedirs(os.path.join(args.out, 'logfiles'))

dirs = [os.path.join(args.out, 'databases'), os.path.join(args.out, 'blast_dir'), 
os.path.join(args.out, 'logfiles')]
for d in dirs: # double check all directories exist
    if not os.path.isdir(d):
        os.makedirs(d)
result_dir = os.path.abspath(os.path.join(currentdir, args.out))
databases = os.path.abspath(os.path.join(result_dir, 'databases'))
blast_dir = os.path.abspath(os.path.join(result_dir, 'blast_dir'))
log_dir = os.path.abspath(os.path.join(result_dir, 'logfiles'))

##### Start main run #####

def makeblast_database(n, fasta):
    logfile = os.path.abspath(os.path.join(log_dir, n+'_db.log'))
    database_present = None
    extensions = ['.phr', '.pin', '.pog', '.psd', '.psi', '.psq']
    for ext in extensions:
        database_ext = os.path.abspath(os.path.join(databases, n+ext))
        if not os.path.exists(database_ext):
            database_present=False
            break
    if database_present == False:
        if os.path.exists(logfile):
            os.remove(logfile)
        subprocess.call([MAKEBLASTDB, '-in', fasta, '-dbtype', args.dbtype, '-parse_seqids', 
        '-out', n, '-logfile', logfile], cwd=databases)
        with open(logfile, 'r') as log: # check for errors
            if 'error' in log.read():
                print('ERROR: There was an error when making the blastdb for {}, check the logfile '\
                'in {}'.format(n, logfile))
                sys.exit()
    else:
        pass

def perform_blast_search(n, fasta, final_species):
    logfile = os.path.abspath(os.path.join(log_dir, n+'_blast.log'))
    db_path = os.path.abspath(os.path.join(databases, n))
    blast_out = os.path.abspath(os.path.join(blast_dir, n))
    if not os.path.exists(blast_out):
        with open(logfile, 'w') as log:
            subprocess.call([BLAST, '-db', db_path, '-query', fasta, '-out', n, '-evalue', 
            str(args.evalue), '-outfmt', '6', '-max_target_seqs', '2', '-max_hsps', '1', '-num_threads', 
            str(args.cpus)], cwd=blast_dir, stdout=log, stderr=log)
        with open(logfile, 'r') as log: # check for errors
            if 'error' in log.read():
                print('ERROR: There was an error producing the blast results for {}, check the logfile '\
                'in {}'.format(n, logfile))
                sys.exit()
    else:
        print('Previous blast results for {} found, will use those.'.format(n))
    if os.path.getsize(blast_out) > 0:
        blast_dict[n] = []
        final_species.append(n)
        with open(blast_out, 'r') as b_out:
            out_linearr = [line.strip().split('\t') for line in b_out]
            for line in out_linearr:
                if not line[0] == line[1]: # if query gene matches subject gene, skip
                    blast_dict[n].append(float(line[2])) # append % identify to dictionary

    elif os.path.getsize(blast_out) == 0:
        print('WARNING: The blast results for {} were blank, we are skipping this file'.format(n))
        # Since file is empty, do not append to final_species
    return blast_out

def get_counts_greater_80_percent(n, fasta, file):
    '''This function was added specifically for Wyka et al. Submitted for which 
    this script was written. It is provided here for clarity of data sharing but 
    is turned off by default. Currently, this was a quick and dirty method to get total 
    counts of gene pairs with >= 80% BLASTp identity and determine if those gene 
    pairs are separeted by 0 genes (i.e. directly next to each other) or separated by 
    5 or less genes. It conforms strictly to my file strucutre and is not robust to other 
    data strucutres. However, this section can be made more robust if requested.
    
    Alternatively, you can alter your data to conform to these rules:
    
    1. First column of batch file should be Genus_species_strain
    2. Parent directory of FASTA file for particualr strain should be named as the given "strain" (i.e. LM4)
    3. Wihtin the LM4 parent directory there should be a GFF3 file named LM4.gff3
    4. Gene names are ordered and labeled such that gene_001-mRNA is followed by gene_002-mRNA, which is 
    followed by gene_003-mRNA on the genome and in the GFF3 file.
    '''
    all_80_dict[n] = []
    contig_80_dict[n] = []
    gene_list = []
    count_next_to = 0
    count_5_or_less = 0
    with open(file, 'r') as b_out:
        out_linearr = [line.strip().split('\t') for line in b_out]

        iso = n.rsplit('_', 1)[1]+'.gff3' # get prefix of file name from 'strain' info provided in batch file
        # i.e. Claviceps_purpurea_LM4 = LM4.gff3 which is located in the same folder as LM4_proteins.fasta
        if iso == '20.1.gff3':
            iso = 'Cpurp.gff3'
        elif iso == '1980.gff3':
            iso = 'Clfusi.gff3'
        elif iso == '1481.gff3':
            iso = 'Clpasp.gff3'
        else:
            iso = iso

        gff_dict = OrderedDict() # {Gene_name : contig}
        with open(os.path.join(os.path.dirname(fasta),iso), 'r') as in_gff3:
            for gene in in_gff3:
                if not gene.startswith('#'):
                    col = gene.split('\t')
                    if col[2] == 'gene':
                        gene_name = re.search(r'(ID=)([^;]*)', col[8]).group(2)+'-T1'
                        gff_dict[gene_name] = col[0]

        for line in out_linearr:
            if not line[0] == line[1]: # if query gene matches subject gene, skip
                if float(line[2]) >= 80.00: # if blastp identify > 80%
                    pair = (line[0],line[1])
                    alt_pair = (line[1],line[0])
                    if pair in all_80_dict[n] or alt_pair in all_80_dict[n]:
                        pass
                    else:
                        all_80_dict[n] = all_80_dict[n] + [pair]
                        # append genes to gene list to get unique gene number for later
                        gene_list.append(pair[0])
                        gene_list.append(pair[1])
                #if genes are on the same contig and blastp identify > 80%
                if gff_dict[line[0]] == gff_dict[line[1]] and float(line[2]) >= 80.00:
                    pair = (line[0],line[1])
                    alt_pair = (line[1],line[0])
                    if pair in contig_80_dict[n] or alt_pair in contig_80_dict[n]:
                        pass
                    else:
                        contig_80_dict[n] = contig_80_dict[n] + [pair]
                        
                        # Look to see if genes are next to each other or separeted by 5 or less genes. 
                        # NOTE: this isn't that robust to different datasets, but it works here due 
                        # to the fact that Funannotate adds locus tags in order. For example, gene_001 is 
                        # always before gene_002 which is followed by gene_003. This can be confirmed 
                        # via a GFF file. Therefore, if we know the genes are on the same contig, if 
                        # we simple to 1 - 3 = 2 we know that the gene_001 and gene_003 are separeted 
                        # by 1 gene (i.e. gene_002).
                        # This splicing function is also specific to my data
                        gene_1 = int(pair[0].split('_')[1].split('-')[0])
                        gene_2 = int(pair[1].split('_')[1].split('-')[0])
                        if abs(gene_1 - gene_2) == 1:
                            count_next_to += 1
                        if abs(gene_1 - gene_2) <= 6: # if 6 that means 5 genes separete them
                            count_5_or_less += 1
    
    next_prop = round(((count_next_to/len(all_80_dict[n]))*100),2) # proportion of gene pairs separeted by 0 genes
    five_prop = round(((count_5_or_less/len(all_80_dict[n]))*100),2) # proportion of gene pairs separeted by <=5 genes
    uniq_genes = len(set(gene_list)) # number of unique genes from unique pairs

    gene_80_counts_out = os.path.abspath(os.path.join(result_dir, 'Counts_genes_pairs_with_80_identity.txt'))
    with open(gene_80_counts_out, 'a') as count_out:
        count_out.write(n + ':\n')
        count_out.write('Number of unique gene pairs with >= 80% identity = ' + 
            str(len(all_80_dict[n])) + '\n')
        count_out.write('Number of unique genes found in unique gene pairs = ' + 
            str(uniq_genes) + '\n')
        count_out.write('Number of these pairs that are separated by 0 genes = ' + 
            str(count_next_to) + '(' + str(next_prop)+'%)' + '\n')
        count_out.write('Number of these pairs that are separated by <= 5 genes = ' + 
            str(count_5_or_less) + '(' + str(five_prop)+'%)' + '\n')

    return all_80_dict

def bin_data(n):
    '''Bin % identities by closest rounded number'''
    bin_dict[n] = {i: 0 for i in range(0,101)}
    for percent_identity in blast_dict[n]:
        bin_dict[n][round(percent_identity)] += 1 # bin by rounding and get number in each bin
    for i in range(0, 101):
        bin_dict[n][i] = bin_dict[n][i]/len(blast_dict[n]) # get proportion from genes used

def dict_to_dataframe():
    '''Convert dictionary to dataframe'''
    df = pd.DataFrame() # create blank dataframe
    count = 0
    for name in [name[0] for name in input_tuples]: # loop through each species and create dataframe
        count += 1
        if count == 1: # fill in dataframe with first data
            df = pd.DataFrame(data=bin_dict[name], index=[name]).transpose()
        else: # merge first data frame with subsequent dataframes
            tmp_df = pd.DataFrame(data=bin_dict[name], index=[name]).transpose()
            df = df.merge(tmp_df, right_index=True, left_index=True)
    return df

def make_bar_graph():
    '''Generate figure'''
    if args.input:
        fig, ax = plt.subplots(figsize=(args.figsize[0],args.figsize[1]), dpi=args.dpi)
        x_labels = np.array([i for i in range(0,101)])
        data = np.array(df[name])
        plt.bar(x_labels, data, color='blue', edgecolor='black')
        ax.set_axisbelow(True)
        plt.minorticks_on()
        plt.grid(which='minor', axis='y', color='white', linestyle='--', alpha=0.3)
        ax.yaxis.grid(True, linestyle='-', linewidth='1', which='major', color='white')
        ax.tick_params(axis='x', which='minor', bottom=False)
        ax.set_facecolor('gainsboro')
        plt.xlabel('Percent identity')
        plt.ylabel('Fraction of orthologs')
    if args.batch:
        fig = plt.figure(figsize=(args.figsize[0],args.figsize[1]), dpi=args.dpi)
        ax = fig.add_subplot(111, projection='3d')
        y_ticks = [i for i in range(1,len(input_tuples)+1)]
        data = np.array([df[name[0]] for name in input_tuples])
        if categories:
            if len(categories) > 10: # if > 10 select from rainbow
                colors = plt.cm.tab10([1.*i/len(categories) for i in range(len(categories))])
            else: # if <= 10 use more contrasting colors
                colors = plt.cm.tab10([i for i in range(len(categories))])
                # colors = ['tab:red', 'tab:blue', 'tab:green', 'tab:orange']
            color_nest = zip(categories, colors)
            color_dict = {nest[0] : nest[1] for nest in color_nest} # {categoty : color}
            for n, k, z in zip(input_tuples, y_ticks, data):
                xs = np.array([i for i in range(0,101)]) # x_axis data
                ys = z # y_axis data 
                cs = color_dict[n[2]]
                ax.bar(xs, ys, zs=k, zdir='y', color=cs, alpha=0.8)
            # manually create handles and patches for legend
            patch_list = []
            for category in color_dict.keys():
                color_key = mpatches.Patch(color=color_dict[category], label=category)
                patch_list.append(color_key)
            plt.legend(handles=patch_list,framealpha=1.0, prop={'size': args.legend})
        else:
            if len(input_tuples) > 10:
                colors = plt.cm.gist_rainbow([1.*i/len(input_tuples) for i in range(len(input_tuples))])
            else:
                colors = plt.cm.tab10([i for i in range(len(input_tuples))])
            for c, k, z in zip(colors, y_ticks, data):
                xs = np.array([i for i in range(0,101)])
                ys = z
                cs = [c] * len(ys)
                ax.bar(xs, ys, zs=k, zdir='y', color=cs, alpha=0.8)
        y_labels = [name[0].replace('_',' ') for name in input_tuples]
        ax.set_yticks(y_ticks) # had to reset the ticks to get the labels to align properly
        ax.set_yticklabels(y_labels, fontsize=args.yfontsize, va='center_baseline', ha='left', rotation=args.yrotation)
        ax.set_xlabel('Identity (%)', labelpad = 15, fontsize=12)
        ax.set_zlabel('Fraction of orthologs', fontsize=12)
        ax.view_init(elev=args.elev, azim=args.azim)
    ax.set_xlim((args.xlimits[0],args.xlimits[1]))
    # plt.show()
    plt.savefig(figure_output)
    plt.close()

if __name__ == "__main__":

    # make blast database(s)
    print('Making blast database(s)')
    if args.batch:
        for tup in input_tuples:
            makeblast_database(tup[0], tup[1])
    if args.input:
        makeblast_database(name, input_fasta)

    # perform blast search(es)
    print('Performing blast searches')
    blast_dict = {}
    # all_80_dict = {}
    # contig_80_dict = {}
    final_species = [] # list of species that had blast_results
    if args.batch:
        for tup in input_tuples:
            blast_out = perform_blast_search(tup[0], tup[1], final_species)
            # get_counts_greater_80_percent(tup[0], tup[1], blast_out)
        if not len(input_tuples) == len(final_species): # double check to remove blank data
            input_tuples = [tup for tup in input_tuples if tup[0] in final_species]
    if args.input:
        perform_blast_search(name, input_fasta)

    # print out unique gene pairs with >= 80% identity per strain
    # this feature is turned off by default and was specific to a publication
    # a_80_df = pd.DataFrame.from_dict(data=all_80_dict, orient='index').T
    # df_80_out = os.path.abspath(os.path.join(result_dir, 'All_unique_gene_pairs_with_80_identity.tsv'))
    # a_80_df.to_csv(df_80_out, sep='\t')

    # bin the data
    print('Binning data')
    bin_dict = {}
    if args.batch:
        for tup in input_tuples:
            bin_data(tup[0])
    if args.input:
        bin_data(name)

    # create dataframe from bin_dict
    print('Creating dataframe from binning dictionary')
    if args.batch:
        df = dict_to_dataframe()
    if args.input:
        df = pd.DataFrame(data=bin_dict[name], index=[name]).transpose()
    dataframe_output = os.path.abspath(os.path.join(result_dir, args.out+'_dataframe.csv'))
    df.to_csv(dataframe_output)

    # make figure
    print('Creating figure')
    figure_output = os.path.abspath(os.path.join(result_dir, args.out+'.png'))
    make_bar_graph()