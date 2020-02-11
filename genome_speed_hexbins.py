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

import os, re, sys, argparse, inspect
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter
from statsmodels.stats.multitest import multipletests
from random import sample

rundir = os.getcwd()

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='./%(prog)s [options] -i GFF3/Bed file -o output',
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
    '--gridsize',
    nargs = '+',
    default = [50, 50],
    help = 'Gridsize for hexes. Larger number are smaller hexes, smaller are larger \
        hexes. [default = 50 50]',
    type=int,
    metavar=''
)
parser.add_argument(
    '--mincount',
    default = 1,
    help = 'Minimum count of genes required to display a hex bin. Increasing this value \
    can improve plots if the genome has multiple compartments. However, if genomes have \
    one compartment then a better respresentation is provided with the default value. [default = 1]',
    type=int,
    metavar=''
)
parser.add_argument(
    '-s',
    '--stats',
    default='resample',
    choices=['resample','allvall'],
    help = 'Type of statistics to use for statistical comparisons. \
    [resample|allvall] [defaut: resample]',
    metavar=''
)
parser.add_argument(
    '-t',
    '--test',
    default='mannwhitney',
    choices=['mannwhitney','students', 'ks_2samp'],
    help = 'Statistical test to use for mean comparisons of overlaid genes \
        [mannwhitney|students|ks_2samp] [defaut: mannwhitney]',
    metavar=''
)
parser.add_argument(
    '-n_s',
    '--n_size',
    default=100,
    type=int,
    help = 'Gene number sample size of overlaid genes for random resampling [defaut: 100]',
    metavar=''
)
parser.add_argument(
    '-re',
    '--resample',
    default=1000,
    type=int,
    help = 'Number of resamplings to perform [defaut: 1000]',
    metavar=''
)
parser.add_argument(
    '-multi',
    '--multitest',
    default='fdr_bh',
    choices=['bonferroni', 'sidak', 'holm-sidak','holm','simes-hochberg','hommel',\
        'fdr_bh','fdr_by','fdr_tsbh','fdr_tsbky'],
    help = 'Multi-test correction to use [bonferroni|sidak|holm-sidak|holm|'\
        'simes-hochberg|hommel|fdr_bh|fdr_by|fdr_tsbh|fdr_tsbky] [defaut: fdr_bh]',
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
    help = "Name to give to genes overlaid ontop of plot for legend, use quotes \
        if using spaces."\
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
    help = "Marker style for overlaid genes (make sure to use quotes: 'o', '<', \
        's',etc) [default = '^']",
    metavar=''
)
parser.add_argument(
    '--marker_color',
    default = 'red',
    help = 'Color of markers to use for overlaid genes (matplotlib color names or \
        values) [default: red]',
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
                        # gene_name = re.search(r'(ID=)([^;]*)', col[8]).group(2)+'-T1'
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
                    if flank_3 <= 0: # check for overlap
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
                    if flank_5 <= 0:# check for overlap 
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
                    if flank_5 <= 0:
                        flank_5 = args.ylimits[0]
                    flank_3 = (genes[gene+1][0]-genes[gene][1])/1000 # distance of start of next gene to stop of current gene
                    if flank_3 <= 0:
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
    ax.tick_params(axis='both', labelsize=13)
    plt.hexbin(prime3_flank, prime5_flank, xscale='log', 
        yscale='log', gridsize=((args.gridsize[0],args.gridsize[0])), mincnt=args.mincount) # plot hexbins
    cb = plt.colorbar(shrink=0.25, aspect=5) # add color bar for hexbins
    cb.ax.set_title('Counts', fontsize=12)
    ax.set_facecolor('lightgrey')
    ax.set_axisbelow(True)
    plt.minorticks_on()
    ax.grid(which='major', linestyle='-', linewidth='1', color='white')
    plt.grid(which='minor', axis='both', color='white', linestyle='--', alpha=0.3)
    plt.xlabel("3' Flanking intergenic region (kbp)", fontsize=14)
    plt.ylabel("5' Flanking intergenic region (kbp)", fontsize=14)
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
            plt.title(header, fontsize=14)
        else:
            top_axes.set_title(header, fontsize=14)
    else:
        pass
    # plt.show()
    plt.savefig(output_fig)
    plt.close()

def get_normal_distribution(values):
    values = sorted(np.log10(values)) # get log of the values
    mean = np.mean(values)
    std = np.std(values)
    lognorm_fit = stats.norm.pdf(values, loc=mean, scale=std) # get norm pdf
    return lognorm_fit

def create_normal_distribution_plots(axes):
    divider = make_axes_locatable(axes)
    ax_norm_x = divider.append_axes('top', 0.5, pad=0.05, sharex=axes)
    ax_norm_y = divider.append_axes('right', 0.5, pad=0.05, sharey=axes)

    ax_norm_x.xaxis.set_tick_params(which='both', labelbottom = False, bottom = False, labelsize=14)
    ax_norm_y.yaxis.set_tick_params(which='both', labelleft = False, left = False, labelsize=14)

    # filter lists to make sure we don't use data points that were positioned along axes (i.e. = 0).
    # These datapoints are either the 5' or 3' ends of a gene that are at the start or end of a contig
    # so the true intergenic length is unknown so we cannot use that data. These data points can also
    # come from genes that had an overlapping 5' and 3' end. Since they overlapped the intergenic length
    # could be classified as 0, but depending on the quality of the genome assembly and annotation this
    # could heavily bias results. While if overlap is generally low, removing these data points should not
    # affect the overall results as the results will still be computed one thousands of genes.
    # Up for debate.
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
    Determines significance between distributions of overlaid genes and all other genes 
    not in the overliad list on the 5' and 3' side. 
    
    Default = Resample (can be more accurate)
    
    Will subsample (--resample number of times) the overlaid genes and all other genes by 
    --n_size and will compute uncorrected p-values. The uncorrected p-values will be corrected 
    using statsmodels.stats.multitest.multipletests correction [Default: Benjamini/Hochberg]. 
    The corrected p-values for the 5' and 3' sides will be averaged and the two resulting 
    p-values will be averaged together for an overall p-value. Note: If there are not enough 
    overlaid genes to select --n_size times, then the n_size will be automatically reduced to 
    75% of the number of genes available from the overlaid list.
    
    Alternative = allvall (potentially higher risk of type I errors)
    
    Will perform a test comparison of ALL overlaid genes versus ALL other genes (i.e. Mann-Whitney U 
    , Welch's test (i.e. unequal variance students), or Kolmogorov-Smirnov 2 sample test.
    '''
    if args.stats == 'resample':
        if len(over_x) < args.n_size or len(over_y) < args.n_size:
            n_size = int(min([len(over_x), len(over_y)])*0.75)
            print(('WARNING: The number of genes to select for random resampling is {}, '
            'however, there were only {} {} to select from. {} number of genes will be used '
            'to select for random resampling').format(args.n_size, n_size/0.75, args.name, n_size))
        else:
            n_size = args.n_size

        top_ps = [] # uncorrected p-values for top plot (x-axis)
        right_ps = [] # uncorrected p-values for right plot (y-axis)
        for i in range(0, args.resample):
            if args.test == 'mannwhitney':
                top_t, top_p = stats.mannwhitneyu(sample(over_x,n_size), 
                    sample(prime3,args.n_size), alternative='two-sided')
                right_t, right_p = stats.mannwhitneyu(sample(over_y,n_size), 
                    sample(prime5,args.n_size), alternative='two-sided')
            elif args.test == 'students':
                top_t, top_p = stats.ttest_ind(sample(over_x,n_size), 
                    sample(prime3,args.n_size), equal_var=True)
                right_t, right_p = stats.ttest_ind(sample(over_y,n_size), 
                    sample(prime5,args.n_size), equal_var=True)
            elif args.test == 'ks_2samp':
                top_d, top_p = stats.ks_2samp(over_x, prime3)
                right_d, right_p = stats.ks_2samp(over_y, prime5)
            top_ps.append(top_p)
            right_ps.append(right_p)
        t_r, t_p, t_sf, t_bf = multipletests(top_ps, alpha=0.05, method=args.multitest)
        r_r, r_p, r_sf, r_bf = multipletests(right_ps, alpha=0.05, method=args.multitest)
        correct_top = np.mean(t_p)
        correct_right = np.mean(r_p)
        overall_p = np.mean([correct_top,correct_right])
    elif args.stats == 'allvall':
        if args.test == 'mannwhitney':
            top_d, top_p = stats.mannwhitneyu(over_x, prime3, alternative='two-sided')
            right_d, right_p = stats.mannwhitneyu(over_y, prime5, alternative='two-sided')
        elif args.test == 'students':
            top_t, top_p = stats.ttest_ind(over_x, prime3, equal_var=False)
            right_t, right_p = stats.ttest_ind(over_y, prime5, equal_var=False)
        elif args.test == 'ks_2samp':
            top_d, top_p = stats.ks_2samp(over_x, prime3)
            right_d, right_p = stats.ks_2samp(over_y, prime5)
        overall_p = np.mean([top_p, right_p])
    return overall_p

if __name__ == "__main__":
    input_file = os.path.abspath(args.input)
    output_fig = os.path.abspath(os.path.join(rundir, args.output+'.png'))
    if '.bed' in input_file:
        contig_dict = parse_bed(input_file)
    elif '.gff3' in input_file or '.gff' in input_file:
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
