#!/usr/bin/python3

'''
This script takes a Gff3 or bed file and a RepeatMasker.output file as input and
calculates the number of TEs surrounding each gene, Surrounding counts are the number of 
TEs on the 5' flanking side between the given gene and the previous gene and the 3' 
flanking side between the given gene and the next gene. A TE_list file can be passed which 
will separate the TEs by type (i.e. LTR, DNA, LINE) and will provide counts for each time 
surrounding each gene. The final output will be a boxplot of the number of TE's 
surrounding each gene as well as a dataframe.tsv file for reproduction purposes.

This script also calcualtes the distance of each gene to the closest TE. This can also 
be separated with a TE_list. The final output will be a boxplot of the distances of the 
closest TE on the 5' flank and 3' flank of all genes and a dataframe.tsv file as well.

You can also provide multiple files of gene lists (-a / --annotations) and the final results 
will be filtered to provide counts or distances for specific genes. For example, if you want 
to add effector and metabolite genes you can pass two files of gene name lists through and 
the resulting data will be report for effector genes, metabolite genes, and all other genes 
that are not effector or metabolite genes.

This script lastly perform Mannwhiteny-U tests for pairwise comparisons if a TE_list or 
annotations files are passed through.
'''


import os, sys, re, argparse, inspect, warnings
import numpy as np
import pandas as pd
from scipy import stats
from collections import OrderedDict
from itertools import combinations
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

rundir = os.getcwd()

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='./%(prog)s [options] -i Gff3/bed -r RepeatMasker.out file -o output',
    description = '''    Takes a Gff3/bed file and a RepeatMasker.output file 
    to calculate the number of TEs surrounding each gene and the distance of each 
    gene to the closest TE.''',
    
    epilog = """Written by Stephen A. Wyka (2019)""",
    formatter_class = MyFormatter)

parser.add_argument(
    '-i',
    '--input',
    required = True,
    help = 'Input Gff3 file or bed file',
    metavar=''
)
parser.add_argument(
    '-r',
    '--repeatmasker',
    required = True,
    help = 'Repeatmasker output file',
    metavar=''
)
parser.add_argument(
    '-l',
    '--te_list',
    help = 'File containing TE types of interest to separate out for bar graphs (one type per line).',
    metavar=''
)
parser.add_argument(
    '-a',
    '--annotations',
    nargs='*',
    help = 'Files containing list of genes to narrow the search and report TEs only surrounding these genes ,\
    (i.e. effector genes, metabolite genes).',
    metavar=''
)
parser.add_argument(
    '-o',
    '--output',
    required = True,
    help = 'Basename of output files (figure and results text file)',
    metavar=''
)
parser.add_argument(
    '-ma',
    '--multialpha',
    default = 0.05,
    type=float,
    help = 'Alpha cut off for multitest correction [default = 0.05]',
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
    '--title',
    help = "Title to be placed over figure (use quotes if using spaces; i.e 'Genus species'",
    metavar=''
)
args=parser.parse_args()

def parse_gff3(input_gff3):
    contig_group_dict = {} # {Contig : list of tuples (start, stop, gene_name) pertaining to this contig}
    gene_list = []
    with open(input_gff3, 'r') as in_gff3:
        for line in in_gff3:
            if not line.startswith('#'):
                col = line.split('\t')
                if col[2] == 'gene':
                    try: # check to make sure 9th columns contains ID= feautres to get gene names
                        gene_name = re.search(r'(ID=)([^;]*)', col[8]).group(2)
                        gene_list.append(gene_name)
                    except:
                        print('ERROR: Problem with gff3 file. Cannot find ID feautre for gene')
                        sys.exit()
                    if col[0] in contig_group_dict.keys():
                        contig_group_dict[col[0]].append((int(col[3]), int(col[4]), gene_name))
                    else:
                        contig_group_dict[col[0]] = [(int(col[3]), int(col[4]), gene_name)]
                else:
                    pass
    return contig_group_dict, gene_list

def parse_bed(input_bed):
    contig_group_dict = {} # {Contig : list of tuples (start, stop, gene_name) pertaining to this contig}
    gene_list = []
    with open(input_bed, 'r') as in_bed:
        for line in in_bed:
            col = line.split('\t')
            if col[3]:
                gene_name = col[3] # assume 4th column is gene_name. As it should be
                gene_list.append(gene_name)
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
    return contig_group_dict, gene_list

def parse_repeatmasker_out(repeat_infile):
    '''
    Example if list is provided:
    repeatdict = {Contig_# : {'LTR' : [(start,stop),(start,stop)], 'DNA' : [(start,stop)], 'LINE' : [(start,stop)]}}
    
    Example if no list is provided:
    This will exclude all simple and low complexity repeats
    repeatdict = {Contig_# : {'repeat' : [(start,stop),(start,stop)]}}
    '''
    with open(repeat_infile, 'r') as infile:
        file_linearr = [line.strip().split() for line in infile][3:] # skip first 3 lines
        repeat_dict = {feat[4] : {} for feat in file_linearr}
        for feat in file_linearr:
            if args.te_list:
                for te_type in type_list:
                    if te_type in feat[10].lower():
                        if te_type in repeat_dict[feat[4]].keys():
                            repeat_dict[feat[4]][te_type].append((int(feat[5]),int(feat[6])))
                        else:
                            repeat_dict[feat[4]][te_type] = [(int(feat[5]),int(feat[6]))]
                    else:
                        pass
            else:
                if 'simple_repeat' not in feat[10].lower():
                    if 'low_complexity' not in feat[10].lower():
                        if 'repeat' in repeat_dict[feat[4]].keys():
                            repeat_dict[feat[4]]['repeat'].append((int(feat[5]),int(feat[6])))
                        else:
                            repeat_dict[feat[4]]['repeat'] = [(int(feat[5]),int(feat[6]))]
                    else:
                        pass
                else:
                    pass
    return repeat_dict

def look_for_closest_5_te(contig, te_type, ge_start, ge):
    try: # look for closet TE from gene if TE stop is < gene start
        tes = [x for x in repeat_dict[contig][te_type]]
        closest_5_te = min([x for x in repeat_dict[contig][te_type] if x[1] < ge_start], 
        key=lambda x:abs(ge_start-x[1]))
        flank_5_dist = ge_start - closest_5_te[1] # distance = start of gene - TE stop
    except:
        flank_5_dist = 0
    return flank_5_dist

def look_for_closest_3_te(contig, te_type, ge_stop):
    try: # look for closet TE from gene is TE start in > gene stop
        tes = [x for x in repeat_dict[contig][te_type] if x[0] > ge_stop]
        closest_3_te = min([x for x in repeat_dict[contig][te_type] if x[0] > ge_stop], 
        key=lambda x:abs(x[0]-ge_stop))
        flank_3_dist = closest_3_te[0] - ge_stop # distance = start of TE - stop of gene
    except:
        flank_3_dist = 0
    return flank_3_dist

def get_average_distance(flank_5, flank_3):
    '''
    If a dist = 0 we can't be sure of distance, so best guess is to just use other distance
    '''
    if flank_5 == 0 and flank_3 != 0:
        avg_flank_te_dist = flank_3/1000
    elif flank_3 == 0 and flank_5 != 0:
        avg_flank_te_dist = flank_5/1000
    elif flank_3 != 0 and flank_5 != 0:
        avg_flank_te_dist = round((flank_5 + flank_3)/2)/1000
    else:
        avg_flank_te_dist = 0 
    return avg_flank_te_dist

def compare_dictionaries():
    '''
    For each gene get the number of TE fragments (by type) that are between that given gene 
    and the next gene, as well as, the given gene and the previous gene. By summing these two 
    numbers we get a total number of TE fragments that surround a given gene. 
    '''
    gene_te_count_dict = {}
    gene_te_dist_dict = {}
    for contig, genes in contig_dict.items():
        if contig in repeat_dict.keys():
            if len(genes) == 1:
                gene_name = genes[0][2]
                gene_te_count_dict[gene_name] = {}
                gene_te_dist_dict[gene_name] = {}
                for te_type in repeat_dict[contig].keys():
                    # flank_5_count if stops (x[1]) of TE's are < than start of gene (genes[0][0])
                    flank_5_count = len([x for x in repeat_dict[contig][te_type] if x[1] < genes[0][0]])
                    flank_5_dist = look_for_closest_5_te(contig, te_type, genes[0][0], genes[0])
                    # flank_3_count if starts (x[0]) of TE's are > than stop of gene (genes[0][1])
                    flank_3_count = len([x for x in repeat_dict[contig][te_type] if x[0] > genes[0][1]])
                    flank_3_dist = look_for_closest_3_te(contig, te_type, genes[0][1])
                    avg_flank_te_dist = get_average_distance(flank_5_dist, flank_3_dist)
                    sum_flank_te = (flank_5_count + flank_3_count)
                    gene_te_count_dict[gene_name][te_type] = sum_flank_te
                    gene_te_dist_dict[gene_name][te_type] = avg_flank_te_dist
            if len(genes) > 1:
                for g in range(0, len(genes)):
                    gene_name = genes[g][2]
                    gene_te_count_dict[gene_name] = {}
                    gene_te_dist_dict[gene_name] = {}
                    for te_type in repeat_dict[contig].keys():
                        if g == 0: # if gene is the first gene in the contig
                            flank_5_count = len([x for x in repeat_dict[contig][te_type] if x[1] < genes[g][0]])
                            flank_5_dist = look_for_closest_5_te(contig, te_type, genes[g][0], genes[g])
                            # if starts (x[0]) of TE's are > stop of gene but < start of next gene
                            flank_3_count = len([x for x in repeat_dict[contig][te_type] if 
                                genes[g+1][0] > x[0] > genes[g][1]])
                            flank_3_dist = look_for_closest_3_te(contig, te_type, genes[g][1])
                            avg_flank_te_dist = get_average_distance(flank_5_dist, flank_3_dist)
                            sum_flank_te = (flank_5_count + flank_3_count)
                            gene_te_count_dict[gene_name][te_type] = sum_flank_te
                            gene_te_dist_dict[gene_name][te_type] = avg_flank_te_dist
                        elif g == len(genes) - 1:
                            # if stops (x[1]) of TE's are < start of gene but > stop of previous gene
                            flank_5_count = len([x for x in repeat_dict[contig][te_type] if 
                                genes[g-1][1] < x[1] < genes[g][0]])
                            flank_5_dist = look_for_closest_5_te(contig, te_type, genes[g][0], genes[g])
                            flank_3_count = len([x for x in repeat_dict[contig][te_type] if x[0] > genes[g][1]])
                            flank_3_dist = look_for_closest_3_te(contig, te_type, genes[g][1])
                            avg_flank_te_dist = get_average_distance(flank_5_dist, flank_3_dist)
                            sum_flank_te = (flank_5_count + flank_3_count)
                            gene_te_count_dict[gene_name][te_type] = sum_flank_te
                            gene_te_dist_dict[gene_name][te_type] = avg_flank_te_dist
                        else:
                            flank_5_count = len([x for x in repeat_dict[contig][te_type] if 
                                genes[g-1][1] < x[1] < genes[g][0]])
                            flank_5_dist = look_for_closest_5_te(contig, te_type, genes[g][0], genes[g])
                            flank_3_count = len([x for x in repeat_dict[contig][te_type] if 
                                genes[g+1][0] > x[0] > genes[g][1]])
                            flank_3_dist = look_for_closest_3_te(contig, te_type, genes[g][1])
                            avg_flank_te_dist = get_average_distance(flank_5_dist, flank_3_dist)
                            sum_flank_te = (flank_5_count + flank_3_count)/2
                            gene_te_count_dict[gene_name][te_type] = sum_flank_te
                            gene_te_dist_dict[gene_name][te_type] = avg_flank_te_dist
        else:
            pass

    return gene_te_count_dict, gene_te_dist_dict

def create_dataframe_from_dictionary(dictionary):
    for gene, te_data in dictionary.items(): # need to add in NaN's for missing data for all genes
        if args.te_list:
            for te_type in type_list:
                if te_type not in te_data.keys():
                    dictionary[gene][te_type] = np.nan
                else:
                    pass
        else:
            if 'repeat' not in te_data.keys():
                dictionary[gene]['repeat'] = np.nan

    if args.te_list and not args.annotations:
        nest_data_list = [['gene_name'] + type_list]
        for gene, te_data in dictionary.items():
            tmp_list = [gene]
            for t in range(0,len(type_list)):
                tmp_list.append(te_data[type_list[t]])
            nest_data_list.append(tmp_list)
        te_df = pd.DataFrame(nest_data_list[1:],columns=nest_data_list[0])
    elif args.annotations and not args.te_list:
        nest_data_list = [['gene_name'] + ['repeat'] + ['Annotations']]
        for gene, te_data in dictionary.items():
            nest_data_list.append([gene, te_data['repeat'], flipped_anno_dict[gene]])
        te_df = pd.DataFrame(nest_data_list[1:],columns=nest_data_list[0])
    elif args.annotations and args.te_list:
        nest_data_list = [['gene_name'] + type_list + ['Annotations']]
        for gene, te_data in dictionary.items():
            tmp_list = [gene]
            for t in range(0,len(type_list)):
                tmp_list.append(te_data[type_list[t]])
            tmp_list.append(flipped_anno_dict[gene])
            nest_data_list.append(tmp_list)
        te_df = pd.DataFrame(nest_data_list[1:],columns=nest_data_list[0])
    else:
        nest_data_list = [['gene_name'] + ['repeat']]
        for gene, te_data in dictionary.items():
            nest_data_list.append([gene, te_data['repeat']])
        te_df = pd.DataFrame(nest_data_list[1:],columns=nest_data_list[0])
    return te_df

def compute_stats(data, labels, type):
    data_combos = list(combinations(data,2))
    label_combos = list(combinations(labels,2))
    stats_output = os.path.abspath(os.path.join(rundir, args.output+'_'+type+'_statistics.txt'))
    uncor_p_list = []
    with open(stats_output, 'w') as stats_out:
        stats_out.write('dataset_1\tdataset_2\tP-value\n')
        for i in range(0, len(data_combos)):
            w, p = stats.mannwhitneyu(data_combos[i][0], data_combos[i][1], alternative = 'two-sided')
            uncor_p_list.append(p)
            stats_out.write('\t'.join(label_combos[i]) + '\t' + str(p) + '\n')
        r, c_p, sf, bf = multipletests(uncor_p_list, alpha=args.multialpha, method=args.multitest)
        corrected_p = [str(p) for p in c_p]
        stats_out.write('Corrected p-values in order as above:\n' + ','.join(corrected_p) + '\n')
        if len(data) >= 3: # try 3 samples and 4 samples. If error tell user to alter this section
            try:
                s, kw_p = stats.kruskal(data[0],data[1],data[2],data[3], nan_policy='omit')
                stats_out.write('Krustal-Wallis\t' + str(kw_p) + '\n')
            except:
                try:
                    s, kw_p = stats.kruskal(data[0],data[1],data[2], nan_policy='omit')
                    stats_out.write('Krustal-Wallis\t' + str(kw_p) + '\n')
                except:
                    print('We took a guess and tried to run Krustal-Wallis with 3 or 4 annotations and tried to. '\
                    'It appears we were wrong and you have more than four sample, to fix this alter lines 343 in '\
                    'the script to support the number of samples you are trying to run statistics on')

def create_boxplot(df, type, outfile):
    fig, ax = plt.subplots(figsize=(args.figsize[0],args.figsize[1]), dpi=args.dpi)

    if args.te_list and not args.annotations:
        x_data = type_list
        colors = plt.cm.tab10([i for i in range(len(x_data))])
        te_data = [df[te_type].tolist() for te_type in type_list]
        data = np.column_stack(te_data)
        mask = ~np.isnan(data)
        filter_dist = [d[m] for d, m in zip(data.T, mask.T)]
        if type == 'distances':
            compute_stats(filter_dist, x_data, 'distance')
        elif type == 'counts':
            compute_stats(filter_dist, x_data, 'counts')
    elif args.annotations and not args.te_list:
        x_data = list(set(df['Annotations']))
        colors = plt.cm.tab10([i for i in range(len(x_data))])
        te_data = [df.loc[df['Annotations'] == x, 'repeat'].tolist() for x in x_data]
        max_length = max([len(x) for x in te_data])
        appended_te_data = []
        for data in te_data: # need to make all list the same length by adding NaN's
            if len(data) == max_length:
                appended_te_data.append(data)
            else:
                add_len = max_length - len(data)
                tmp_data = np.append(data, np.repeat(np.nan, add_len))
                appended_te_data.append(tmp_data.tolist())
        data = np.column_stack(appended_te_data)
        mask = ~np.isnan(data)
        filter_dist = [d[m] for d, m in zip(data.T, mask.T)]
        if type == 'distances':
            compute_stats(filter_dist, x_data, 'distance')
        elif type == 'counts':
            compute_stats(filter_dist, x_data, 'counts')
    # elif args.annotations and args.te_list:
        # currently do not have support for this for figure creation yet
    else:
        x_data = ['repeat']
        colors = plt.cm.tab10([i for i in range(len(x_data))])
        te_data = df['repeat'].tolist()
        filter_dist = [x for x in te_data if str(x) != 'nan']
    flierprops = dict(marker='.', markerfacecolor='black', markersize=2,linestyle='none')

    ax.set_xticklabels(x_data, rotation=45)
    ax.set_axisbelow(True)
    plt.minorticks_on()
    plt.grid(which='minor', axis='y', color='white', linestyle='--', alpha=0.3)
    ax.yaxis.grid(True, linestyle='-', linewidth='1', which='major', color='white')
    ax.tick_params(axis='x', which='minor', bottom=False)
    ax.set_facecolor('gainsboro')
    bp = plt.boxplot(filter_dist, '', patch_artist=True, flierprops=flierprops, widths=0.75)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    plt.setp(bp['boxes'], linewidth=1)
    plt.setp(bp['whiskers'], color='black', linewidth=1)
    plt.setp(bp['medians'], color='black', linewidth=1)
    plt.setp(bp['caps'], color='black', linewidth=1)
    pos = range(len(bp['boxes']))
    if type == 'distances':
        plt.ylabel('Average distance to TE (kbp)')
    if type == 'counts':
        plt.ylabel('Number of TEs')
    ax.set_aspect(0.07)
    if args.title:
        if len(args.title.split(' ')) == 2:
            header = '$\it{}$ $\it{}$'.format(args.title.split(' ')[0], args.title.split(' ')[1])
        elif len(args.title.split(' ')) == 3:
            header = '$\it{}$ $\it{}$ {}'.format(args.title.rsplit(' ',1)[0].split(' ')[0],
                args.title.rsplit(' ',1)[0].split(' ')[1],args.title.rsplit(' ',1)[1])
        else:
            header = args.title
        plt.title(header)
    else:
        pass
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()
    # plt.show()

if __name__ == "__main__":
    # Get set up
    input_file = os.path.abspath(args.input)
    count_figure = os.path.abspath(os.path.join(rundir, args.output+'_sum_counts.png'))
    distance_figure = os.path.abspath(os.path.join(rundir, args.output+'_kbp_distances.png'))
    count_results = os.path.abspath(os.path.join(rundir, args.output+'_sum_count_results.tsv'))
    distance_results = os.path.abspath(os.path.join(rundir, args.output+'_kbp_distance_results.tsv'))
    if args.te_list:
        list_file = os.path.abspath(args.te_list)
        with open(list_file, 'r') as in_list: # added removal of extra \n for precaution
            label_list = sorted([te_type.strip() for te_type in in_list if te_type != '\n'])
            type_list = [te_type.lower() for te_type in label_list]

    # Parse Gff3 or Bed file
    if '.bed' in input_file:
        contig_dict, gene_list = parse_bed(input_file)
    elif '.gff3' in input_file:
        contig_dict, gene_list = parse_gff3(input_file)
    else:
        print("ERROR: Couldn't distinquish the input file as a .bed or .gff3 file")
        sys.exit()

    # Check that genes from annotation input match gene names found in gff3 or bed
    if args.annotations:
        anno_dict = OrderedDict() # {Annotation Name : [list of genes]}
        filtered_gene_list = gene_list.copy()
        for file in args.annotations:
            apath = os.path.abspath(os.path.join(rundir, file))
            anno_name = os.path.basename(apath).split('.')[0]
            with open(apath, 'r') as in_file:
                # anno_list = [gene.strip().split('-T1')[0] for gene in in_file if gene != '\n']
                anno_list = [gene.strip() for gene in in_file if gene != '\n']
                gene_errors = [anno for anno in anno_list if anno not in gene_list]
                filtered_gene_list = list(set(filtered_gene_list) - set(anno_list))
                anno_dict[anno_name] = anno_list
                if not gene_errors:
                    pass
                if gene_errors: # print out warning of mismatched gene names
                    print('Warning: The following genes in annotation file are not found in gff3/bed file')
                    print(gene_errors)
        anno_dict['Other Genes'] = filtered_gene_list

        flipped_anno_dict = OrderedDict() # Flipped anno_dict to {Gene name: Annotation type}
        for anno, genes in anno_dict.items():
            for gene in genes:
                if gene not in flipped_anno_dict.keys():
                    flipped_anno_dict[gene] = anno
                else: # This will also remove duplicates calls based on order of anno input
                    pass

    # Parse RepeatMasker output file
    repeat_in = os.path.abspath(args.repeatmasker)
    repeat_dict = parse_repeatmasker_out(repeat_in)

    # check for contig matches
    if not [k for k in repeat_dict.keys() if any(x == k for x in contig_dict.keys())]:
        print('ERROR: None of the contig ID in your RepeatMasker file match contig IDs '\
        'in your Gff3/bed file.')
        sys.exit()

    # Compare dictionaries and get count dictionary for dataframe conversion
    gene_te_count_dict, gene_te_dist_dict = compare_dictionaries()
    
    # Convert count dictionary to dataframe
    te_count_df = create_dataframe_from_dictionary(gene_te_count_dict)
    te_dist_df = create_dataframe_from_dictionary(gene_te_dist_dict)
    te_count_df.to_csv(count_results, sep='\t', header=True, index=False)
    te_dist_df.to_csv(distance_results, sep='\t', header=True, index=False)
    
    # Create figure
    create_boxplot(te_count_df, 'counts', count_figure)
    create_boxplot(te_dist_df, 'distances', distance_figure)
