#!/usr/bin/python3

'''
Script to take a gff3 or bed file and create hexbin plots of intergenic lengths (kbp) of 
genes. You can also overlay specific genes onto of the hexbins in a scatter plot through 
input of a file containing a list of genes. make sure the gene names in the list file 
match names in the 4th column of the .bed file or in the gff3 'gene' features in the 9th 
column with ID=. By default this script will plot all first genes of the contig along the 
X-axis, as the 5' flanking regions are unknown, and the last gene of the contig along the 
Y-axis, as the 3' flanking regions are unknown. You can turn this off by adding the 
--disregard_contig_ends, in which case the script will ignore the first and last gene of each 
contig. Note: this script also ingores all genes where that gene is the only gene located on 
a contig.
'''

import os, re, sys, argparse
import matplotlib.pyplot as plt

rundir = os.getcwd()

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='./%(prog)s [options] -i gff3/bed file -o output',
    description = '''    Takes gff3 or bed files and creates hexbin plots of intergenic 
    lengths (kbp) of genes.''',
    
    epilog = """Written by Stephen A. Wyka (2019)""",
    formatter_class = MyFormatter)

parser.add_argument(
    '-i',
    '--input',
    required = True,
    help = 'Input file (gff3 or bed file)',
    metavar=''
)
parser.add_argument(
    '-o',
    '--output',
    required = True,
    help = 'Base name of output figure file (will automatically append .png)',
    metavar=''
)
parser.add_argument(
    '-l',
    '--list',
    help = 'List of genes to overlay ontop of plot (i.e. list effector genes)',
    metavar=''
)
parser.add_argument(
    '-n',
    '--name',
    help = "Name to give to genes overlaid ontop of plot for legend, use quotes if using spaces."\
    "(i.e. 'Effector Proteins')",
    metavar=''
)
parser.add_argument(
    '--marker_size',
    type=float,
    default=3,
    help = 'Size of markers to use for overlaid genes [default = 3]',
    metavar=''
)
parser.add_argument(
    '--marker_style',
    default='^',
    help = "Marker style for overlaid genes (make sure to use quotes: 'o', '<', 's',etc) [default = '^']",
    metavar=''
)
parser.add_argument(
    '--marker_color',
    default = 'red',
    help = 'Color of markers to use for overlaid genes (matplotlib color names or values) [default: red]',
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
    '--disregard_contig_ends',
    action='store_true',
    help = 'Turns off the feature to plot first gene in contigs along Y-axis and last'\
    'gene in contigs along the X-axis. [default: ON]',
)
parser.add_argument(
    '--figsize',
    nargs = '+',
    default = [6.4, 4.8],
    type=float,
    help = 'Dimensions of figure in inches [default = 6.4 4.8]',
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
    '--ylimits',
    nargs = '+',
    default = [0.0025, 200],
    type=float,
    help = 'Min and max limits of the Y-axis in log scale [default = 0.0025 200]',
    metavar=''
)
parser.add_argument(
    '--xlimits',
    nargs = '+',
    default = [0.0025, 200],
    type=float,
    help = 'Min and max limits of the X-axis in log scale [default = 0.0025 200]',
    metavar=''
)
parser.add_argument(
    '--gridsize',
    nargs = '+',
    default = [50, 50],
    type=int,
    help = 'Gridsize for hexes. Larger number is smaller hexes, smaller are larger hexes. [default = 50 50]',
    metavar=''
)
args=parser.parse_args()

def parse_gff3(input_gff3):
    contig_group_dict = {} # {Contig : list of tuples (start, stop, gene_name) pertaining to this contig}
    with open(input_gff3, 'r') as in_gff3:
        for line in in_gff3:
            if not line.startswith('#'):
                col = line.split('\t')
                if col[2] == 'gene':
                    try: # check to make sure 9th columns contains ID= feautres to get gene names
                        gene_name = re.search(r'(ID=)([^;]*)', col[8]).group(2)
                    except:
                        print('ERROR: Problem with gff3 file. Cannot find ID feautre for gene')
                        sys.exit()
                    if col[0] in contig_group_dict.keys():
                        contig_group_dict[col[0]].append((int(col[3]), int(col[4]), gene_name))
                    else:
                        contig_group_dict[col[0]] = [(int(col[3]), int(col[4]), gene_name)]
                else:
                    pass
    return contig_group_dict

def parse_bed(input_bed):
    contig_group_dict = {} # {Contig : list of tuples (start, stop, gene_name) pertaining to this contig}
    with open(input_bed, 'r') as in_bed:
        for line in in_bed:
            col = line.split('\t')
            if col[3]:
                gene_name = col[3] # assume 4th column is gene_name. As it should be
                if col[0] in contig_group_dict.keys():
                    contig_group_dict[col[0]].append((int(col[1]), int(col[2]), gene_name))
                else:
                    contig_group_dict[col[0]] = [(int(col[1]), int(col[2]), gene_name)]
            else:
                if args.list: # if they provided a list we expect gene names to be in the 4th column of .bed file
                    print('ERROR: List was provided to overlay specific genes onto plot, however, ,'\
                    'gene names (4th column) are not present in .bed file')
                    sys.exit()
                if col[0] in contig_group_dict.keys():
                    contig_group_dict[col[0]].append((int(col[1]), int(col[2])))
                else:
                    contig_group_dict[col[0]] = [(int(col[1]), int(col[2]))]
    return contig_group_dict

def get_distances_between_genes():
    dist_dict = {} #{Protein_name:(distance to previous gene (5' flank) (kbp), distance to next gene (3' flank) (kbp))}
    first_gene_count = 0
    last_gene_count = 0
    for contig, genes in contig_dict.items():
        if len(genes) > 1: # throw away genes where there is only 1 gene on the contig (can't get distances)
            for gene in range(0, len(genes)):
                if gene == 0: # if gene is the first gene in the contig
                    first_gene_count += 1
                    if not args.disregard_contig_ends: # set flanking 5' end of 1st gene in contig to x-min
                        dist_dict[genes[gene][2]] = [args.xlimits[0], (genes[gene+1][0] - genes[gene][1])/1000]
                    else: # throw away 1st gene in contigs
                        pass
                elif gene == len(genes) - 1: # if gene is the last gene in the contig
                    last_gene_count += 1
                    if not args.disregard_contig_ends: # set flanking 3' end of last gene in contig to y-min
                        dist_dict[genes[gene][2]] = [(genes[gene][0] - genes[gene-1][1])/1000, args.ylimits[0]]
                    else: # throw away last gene in contigs
                        pass
                else: # get distances from start of gene to end of previous gene and end of gene to start of next gene
                    dist_dict[genes[gene][2]] = [(genes[gene][0]-genes[gene-1][1])/1000,
                        (genes[gene+1][0]-genes[gene][1])/1000]
    if args.disregard_contig_ends: # print a warning message to let user know howm any genes they aren't plotting
        print(('You will not be plotting a total of {} genes from the start and end of each contig')
        .format(first_gene_count + last_gene_count))
    return dist_dict

def create_density_plot():
    prime5_flanking_data = []
    prime3_flanking_data = []
    for distances in distance_dict.values(): # create lists of flanking data for plotting
        prime5_flanking_data.append(distances[0])
        prime3_flanking_data.append(distances[1])
    fig, ax = plt.subplots(figsize=(args.figsize[0],args.figsize[1]), dpi=args.dpi)
    ax.set_ylim((args.ylimits[0],args.ylimits[1]))
    ax.set_xlim((args.xlimits[0],args.xlimits[1]))
    plt.hexbin(prime3_flanking_data, prime5_flanking_data, xscale='log', 
        yscale='log', gridsize=((args.gridsize[0],args.gridsize[0])), mincnt=1) # plot hexbins
    plt.colorbar(shrink=0.25, aspect=5) # add color bar for hexbins
    if args.list: # overlay scatter plot of input genes
        overlay_x_data = []
        overlay_y_data = []
        for gene, distances in distance_dict.items():
            if gene in overlay_list:
                overlay_y_data.append(distances[0])
                overlay_x_data.append(distances[1])
        if not len(overlay_x_data) > 0 and len(overlay_y_data) > 0:
            print('ERROR: None of the genes in your list matched the gene names from your input gff3/bed file')
            sys.exit()
        plt.scatter(overlay_x_data, overlay_y_data, s=args.marker_size, c=args.marker_color,marker=args.marker_style,
            label=args.name)
        ax.legend(framealpha=1.0, prop={'size': args.legend}, markerscale=args.marker_size) # add legend for overlay
    ax.set_facecolor('lightgrey')
    ax.set_axisbelow(True)
    plt.minorticks_on()
    ax.grid(which='major', linestyle='-', linewidth='1', color='white')
    plt.grid(which='minor', axis='both', color='white', linestyle='--', alpha=0.3)
    plt.xlabel("3' Flanking intergenic region (kbp)")
    plt.ylabel("5' Flanking intergenic region (kbp)")
    # plt.show()
    plt.savefig(output_fig)

if __name__ == "__main__":
    input_file = os.path.abspath(args.input)
    output_fig = os.path.abspath(os.path.join(rundir, args.output+'.png'))
    if '.bed' in input_file:
        contig_dict = parse_bed(input_file)
    elif '.gff3' in input_file:
        contig_dict = parse_gff3(input_file)
    else:
        print("ERROR: Couldn't distinquish the input file as a .bed or .gff3 file")
        sys.exit()
    if args.list:
        list_file = os.path.abspath(args.list)
        with open(list_file, 'r') as in_list:
            overlay_list = [gene.strip() for gene in in_list]
    distance_dict = get_distances_between_genes()
    create_density_plot()