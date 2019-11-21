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

By default if a list is passed in to overlay specific genes this script will additionally show 
normal distribution plots for the y- and x- axis to show the distributions of the overlaid genes 
versus all other genes that are not in the overlaid list file. In addition, a statistical test of 
distance means will be computed for each axis and averaged together for an overall p-value which 
will be reported in the legend.
'''

import os, re, sys, argparse, random, inspect
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter

rundir = os.getcwd()

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='./%(prog)s [options] -d directory -o output',
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
    '-s',
    '--stats',
    default='welch',
    choices=['welch', 'mannwhitney'],
    help = 'Statistical test to use for mean comparisons of overlaid genes [welch|mannwhitney] [defaut: welch]',
    metavar=''
)
parser.add_argument(
    '--no_norm',
    action='store_true',
    help = 'Do not compute/report normal distributions and statistics [default: ON]'
)
parser.add_argument(
    '-n',
    '--name',
    default = 'Specific genes',
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
    help = 'Gridsize for hexes. Larger number are smaller hexes, smaller are larger hexes. [default = 50 50]',
    type=int,
    metavar=''
)
parser.add_argument(
    '--title',
    help = "Title to be placed over figure (use quotes if using spaces; i.e 'Genus species'",
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
                        gene_name = re.search(r'(ID=)([^;]*)', col[8]).group(2)+'-T1'
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
    overlap_genes = []
    first_gene_count = 0
    last_gene_count = 0
    overlap_count = 0
    for contig, genes in contig_dict.items():
        if len(genes) > 1: # throw away genes where there is only 1 gene on the contig (can't get distances)
            for gene in range(0, len(genes)):
                if gene == 0: # if gene is the first gene in the contig
                    first_gene_count += 1
                    flank_5 = args.ylimits[0] # set flanking 5' end of 1st gene in contig to y-min
                    flank_3 = (genes[gene+1][0] - genes[gene][1])/1000 # distance of start of next gene to stop of current gene
                    if flank_3 < 0: # check for overlap
                        flank_3 = args.xlimits[0] # if overlap set to x-min
                    if not args.disregard_contig_ends: # keep all results
                        if flank_3 == args.xlimits[0]:
                            overlap_count += 1
                            overlap_genes.append(genes[gene][2])
                        dist_dict[genes[gene][2]] = [flank_5, flank_3]
                    else: # throw away 1st gene in contigs
                        pass
                elif gene == len(genes) - 1: # if gene is the last gene in the contig
                    last_gene_count += 1
                    flank_5 = (genes[gene][0] - genes[gene-1][1])/1000 # distance of stop of previous gene to start of current gene
                    if flank_5 < 0:# check for overlap 
                        flank_5 = args.ylimits[0] # if overlap set to y-min
                    flank_3 = args.xlimits[0] # set flanking 3' end of last gene in contig to x-min
                    if not args.disregard_contig_ends: # keep all results
                        if flank_5 == args.ylimits[0]:
                            overlap_count += 1
                            overlap_genes.append(genes[gene][2])
                        dist_dict[genes[gene][2]] = [flank_5, flank_3]
                    else: # throw away last gene in contigs
                        pass
                else: # get distances from start of gene to end of previous gene and end of gene to start of next gene
                    flank_5 = (genes[gene][0]-genes[gene-1][1])/1000 # distance of stop of previous gene to start of current gene
                    if flank_5 < 0:
                        flank_5 = args.ylimits[0]
                    flank_3 = (genes[gene+1][0]-genes[gene][1])/1000 # distance of start of next gene to stop of current gene
                    if flank_3 < 0:
                        flank_3 = args.xlimits[0]
                    if not args.disregard_contig_ends: # keep all results
                        if flank_3 == args.xlimits[0] or flank_5 == args.ylimits[0]:
                            overlap_count += 1
                            overlap_genes.append(genes[gene][2])
                        dist_dict[genes[gene][2]] = [flank_5, flank_3]
                    if args.disregard_contig_ends: # discard altered genes
                        if flank_3 == args.xlimits[0] or flank_5 == args.ylimits[0]:
                            overlap_count += 1
                            overlap_genes.append(genes[gene][2])
                            pass
                        else:
                            dist_dict[genes[gene][2]] = [flank_5, flank_3]
    filt_dist_dict = {}
    for k, v in dist_dict.items():
        if not v[0] == args.xlimits[0] and v[1] == args.ylimits[0]:
            filt_dist_dict[k] = v
        elif v[0] != args.xlimits[0] and v[1] == args.ylimits[0]:
            filt_dist_dict[k] = v
        elif v[0] == args.xlimits[0] and v[1] != args.ylimits[0]:
            filt_dist_dict[k] = v
        elif v[0] != args.xlimits[0] and v[1] != args.ylimits[0]:
            filt_dist_dict[k] = v
    if overlap_genes: # print a warning message to let user know if overlapped genes
        print(("WARNING: {} genes were found to be overlapping, the '5 or 3' ends have been set to 0. The genes that had "
        "overlapping genes are as followed.").format(overlap_count))
        print(overlap_genes)
    if args.disregard_contig_ends: # print a warning message to let user know howm any genes they aren't plotting
        print(('You will not be plotting a total of {} genes from the start and end of each contig')
        .format(first_gene_count + last_gene_count + overlap_count))
    return filt_dist_dict

def create_arrays():
    prime5_flank_data = []
    prime3_flank_data = []
    for distances in distance_dict.values(): # create lists of flanking data for plotting
        prime5_flank_data.append(distances[0])
        prime3_flank_data.append(distances[1])
    if args.list:
        overlay_y_data = []
        overlay_x_data = []
        prime5_flank_no_overlay = []
        prime3_flank_no_overlay = []
        for gene, distances in distance_dict.items():
            if gene in overlay_list:
                overlay_y_data.append(distances[0])
                overlay_x_data.append(distances[1])
            else:
                prime5_flank_no_overlay.append(distances[0])
                prime3_flank_no_overlay.append(distances[1])
        if not len(overlay_x_data) > 0 and len(overlay_y_data) > 0:
            print('ERROR: None of the genes in your list matched the gene names from your input gff3/bed file')
            sys.exit()
        return prime5_flank_data, prime3_flank_data, overlay_y_data, overlay_x_data, prime5_flank_no_overlay, prime3_flank_no_overlay
    else:
        return prime5_flank_data, prime3_flank_data

def create_hexbin_plot():
    fig, ax = plt.subplots(figsize=(args.figsize[0],args.figsize[1]), dpi=args.dpi)

    ax.set_ylim((args.ylimits[0],args.ylimits[1]))
    ax.set_xlim((args.xlimits[0],args.xlimits[1]))
    plt.hexbin(prime3_flank, prime5_flank, xscale='log', 
        yscale='log', gridsize=((args.gridsize[0],args.gridsize[0])), mincnt=1) # plot hexbins
    cb = plt.colorbar(shrink=0.25, aspect=5) # add color bar for hexbins
    cb.ax.set_title('Counts', fontsize=10)
    ax.set_facecolor('lightgrey')
    ax.set_axisbelow(True)
    plt.minorticks_on()
    ax.grid(which='major', linestyle='-', linewidth='1', color='white')
    plt.grid(which='minor', axis='both', color='white', linestyle='--', alpha=0.3)
    plt.xlabel("3' Flanking intergenic region (kbp)")
    plt.ylabel("5' Flanking intergenic region (kbp)")
    if args.list: # overlay scatter plot of input genes
        scatter = plt.scatter(overlay_x, overlay_y, s=args.marker_size, c=args.marker_color,
            marker=args.marker_style)
        if not args.no_norm:
            overall_p, top_axes = create_normal_distribution_plots(ax)
            if overall_p < 0.001:
                final_p = '< 0.001'
            elif overall_p < 0.01:
                final_p = '< 0.01'
            else:
                final_p = '= %.3f' % round(overall_p,3)
                print(final_p)
            scatter.set_label(args.name+ ' ($\it{p}$ %s)' % final_p)
        else:
            scatter.set_label(args.name)

        ax.legend(framealpha=1.0, prop={'size': args.legend}, markerscale=args.marker_size) # add legend for overlay

    if args.title:
        if len(args.title.split(' ')) == 2:
            header = '$\it{}$ $\it{}$'.format(args.title.split(' ')[0], args.title.split(' ')[1])
        elif len(args.title.split(' ')) == 3:
            header = '$\it{}$ $\it{}$ {}'.format(args.title.rsplit(' ',1)[0].split(' ')[0],
                args.title.rsplit(' ',1)[0].split(' ')[1],args.title.rsplit(' ',1)[1])
        else:
            header = args.title
        if args.no_norm:
            plt.title(header)
        else:
            top_axes.set_title(header)
    else:
        pass
    plt.show()
    # plt.savefig(output_fig)
    # plt.close()

def get_normal_distribution(values):
    values = sorted(np.log(values)) # get natural logs of the values
    mean = np.mean(values)
    std = np.std(values)
    lognorm_fit = stats.norm.pdf(values, loc=mean, scale=std) # get norm pdf
    return lognorm_fit

def create_normal_distribution_plots(axes):
    divider = make_axes_locatable(axes)
    ax_norm_x = divider.append_axes('top', 0.5, pad=0.05, sharex=axes)
    ax_norm_y = divider.append_axes('right', 0.5, pad=0.05, sharey=axes)
    
    ax_norm_x.xaxis.set_tick_params(which='both', labelbottom = False, bottom = False)
    ax_norm_y.yaxis.set_tick_params(which='both', labelleft = False, left = False)
    
    # filter lists to make sure we don't use data points that were positioned along axes (i.e. = 0)
    prime3_filt = [x for x in prime3_flank_no_over if x > args.xlimits[0]]
    prime5_filt = [x for x in prime5_flank_no_over if x > args.ylimits[0]]
    overlay_x_filt = [x for x in overlay_x if x > args.xlimits[0]]
    overlay_y_filt = [x for x in overlay_y if x > args.ylimits[0]]
    overall_p = compute_statistics(overlay_x_filt, prime3_filt, overlay_y_filt, prime5_filt)
    
    prime3_norm = get_normal_distribution(prime3_filt)
    prime5_norm = get_normal_distribution(prime5_filt)
    overlay_x_norm = get_normal_distribution(overlay_x_filt)
    overlay_y_norm = get_normal_distribution(overlay_y_filt)
    
    top_axis_max = max(prime3_norm.tolist() + overlay_x_norm.tolist())*1.1
    ax_norm_x.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax_norm_x.set_ylim(ymin=0, ymax=top_axis_max)
    ax_norm_x.plot(sorted(prime3_filt), prime3_norm, color = 'black')
    ax_norm_x.plot(sorted(overlay_x_filt), overlay_x_norm, color = args.marker_color)
    
    right_axis_max = max(prime5_norm.tolist() + overlay_y_norm.tolist())*1.1
    ax_norm_y.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax_norm_y.set_xlim(xmin=0, xmax=right_axis_max)
    ax_norm_y.plot(prime5_norm, sorted(prime5_filt), color = 'black')
    ax_norm_y.plot(overlay_y_norm, sorted(overlay_y_filt), color = args.marker_color)
    return overall_p, ax_norm_x
    
def compute_statistics(over_x, prime3, over_y, prime5):
    '''
    Computes significane differnces between distance means of overlay genes and 
    all other genes on the 5' and 3' side. These two p-values are averages together 
    to report the final p-value of both sides. Default is the independent 
    unequal variance t-test (Welch's test), but user can change to a Mann-Whitney U-test 
    (two-sided) through arguments.
    '''
    if args.stats == 'welch':
        top_t, top_p = stats.ttest_ind(over_x, prime3, equal_var=False)
        right_t, right_p = stats.ttest_ind(over_y, prime5, equal_var=False)
    elif args.stats == 'mannwhitney':
        top_d, top_p = stats.mannwhitneyu(over_x, prime3, alternative='two-sided')
        right_d, right_p = stats.mannwhitneyu(overY, prime5, alternative='two-sided')
    overall_p = np.mean([top_p, right_p])
    return overall_p

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
            overlay_list = [gene.strip() if '\n' in gene else gene for gene in in_list]
    distance_dict = get_distances_between_genes()

    if args.list:
        prime5_flank, prime3_flank, overlay_y, overlay_x, prime5_flank_no_over, prime3_flank_no_over = create_arrays()
    else:
        prime5_flank, prime3_flank = create_arrays()

    create_hexbin_plot()