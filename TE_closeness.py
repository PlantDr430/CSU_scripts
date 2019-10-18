#!/usr/bin/python3

'''
This script akes a Gff3 or bed file and a RepeatMasker.output file as input and
calculates the number of TEs surrounding each gene, Surrounding counts are the number of 
TEs on the 5' flanking side between the given gene and the previous gene and the 3' 
flanking side between the given gene and the next gene. A TE_list file can be passed which 
will separate the TEs by type (i.e. LTR, DNA, LINE) and will provide counts for each time 
surrounding each gene. The final output will be a bar graph of average number of TE's 
with standard error bars surrounding a gene as well as a dataframe.tsv file for 
reproduction purposes.

This script also calcualtes the distance of each gene to the closest TE. This can also 
be separated with a TE_list. The final output will be a boxplot of the distances of the 
 closest TE on the 5' flank and 3' flank of all genes and a dataframe.tsv file as well.

You can also provide a gene list (-a / --annotations) and the final results will be filtered 
to provide counts or distances for specific genes. For example, if you want to only see the 
counts of TEs and distances of TE's around effectors, then you would provide a list of effector 
gene names (that match the gene_names in your Gff3/bed file) and the final results will 
only show results for your effector genes.
'''


import os, sys, re, argparse, inspect, warnings
import numpy as np
import pandas as pd
from scipy import stats
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
    help = 'List of genes to narrow the search and report TEs only surrounding these genes ,\
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
                if 'Simple_repeat' not in feat[10]:
                    if 'Low_complexity' not in feat[10].lower():
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
        avg_flank_te_dist = 0 # no TE to 5' or 3' end, will remove 0 later
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
                    sum_flank_te = flank_5_count + flank_3_count
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
                            sum_flank_te = flank_5_count + flank_3_count
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
                            sum_flank_te = flank_5_count + flank_3_count
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
                            sum_flank_te = flank_5_count + flank_3_count
                            gene_te_count_dict[gene_name][te_type] = sum_flank_te
                            gene_te_dist_dict[gene_name][te_type] = avg_flank_te_dist
        else:
            pass
        
    return gene_te_count_dict, gene_te_dist_dict

def create_dataframe_from_count_dict(dictionary):
    for gene, te_counts in dictionary.items():
        if args.te_list:
            for te_type in type_list:
                if te_type not in te_counts.keys():
                    dictionary[gene][te_type] = np.nan
                else:
                    pass
        else:
            if 'repeat' not in te_counts.keys():
                dictionary[gene]['repeat'] = np.nan
    if args.te_list:
        nest_data_list = [['gene_name'] + type_list]
        for gene, te_counts in dictionary.items():
            tmp_list = [gene]
            for t in range(0,len(type_list)):
                tmp_list.append(te_counts[type_list[t]])
            nest_data_list.append(tmp_list)
        te_df = pd.DataFrame(nest_data_list[1:],columns=nest_data_list[0])
    else:
        nest_data_list = [['gene_name'] + ['repeat']]
        for gene, te_counts in dictionary.items():
            nest_data_list.append([gene, te_counts['repeat']])
        te_df = pd.DataFrame(nest_data_list[1:],columns=nest_data_list[0])
    if args.annotations: # filter te_df by input genes from annotation list
        try:
            filter_df = te_df[te_df['gene_name'].isin(anno_list)]
        except:
            print('There was a problem with filtering for gene names provided in your'\
            'annotation list. Please make sure gene names in your list make gene names'\
            'in your Gff3/bed file. For now all genes will be used.')
            filter_df = te_df
    else:
        filter_df = te_df
    return filter_df

def create_bar_graphs(df):
    fig, ax = plt.subplots()
    if args.te_list:
        x_data = type_list
        warnings.filterwarnings("ignore") # sometimes we get types with all NaN's so I added this to stop the warnings
        te_data = [np.nanmean(df[te_type]) for te_type in type_list]
        te_se = [stats.sem(df[te_type], nan_policy='omit') for te_type in type_list]
    else:
        x_data = ['repeat']
        te_data = np.nanmean(df['repeat'])
        te_se = stats.sem(df['repeat'], nan_policy='omit')
    plt.bar(x_data, te_data, yerr=te_se, capsize=12, color=colors, edgecolor='black')
    ax.set_xticklabels(x_data)
    plt.xticks(rotation=45)
    ax.set_axisbelow(True)
    plt.minorticks_on()
    plt.grid(which='minor', axis='y', color='white', linestyle='--', alpha=0.3)
    ax.yaxis.grid(True, linestyle='-', linewidth='1', which='major', color='white')
    ax.tick_params(axis='x', which='minor', bottom=False)
    ax.set_facecolor('gainsboro')
    plt.xlabel('Repeat type')
    plt.ylabel('Average number of flanking repeats')
    plt.tight_layout()
    plt.savefig(count_figure)
    plt.close()
    # plt.show()

def create_boxplot(df):
    fig, ax = plt.subplots()
    if args.te_list:
        x_data = type_list
        te_data = [df[te_type].tolist() for te_type in type_list]
        distances = np.column_stack(te_data)
        mask = ~np.isnan(distances)
        filter_dist = [d[m] for d, m in zip(distances.T, mask.T)]
    else:
        x_data = ['repeat']
        te_data = df['repeat'].tolist()
        filter_dist = [x for x in te_data if str(x) != 'nan']
    ax.set_xticklabels(x_data)
    ax.set_axisbelow(True)
    plt.minorticks_on()
    plt.grid(which='minor', axis='y', color='white', linestyle='--', alpha=0.3)
    ax.yaxis.grid(True, linestyle='-', linewidth='1', which='major', color='white')
    ax.tick_params(axis='x', which='minor', bottom=False)
    ax.set_facecolor('gainsboro')
    bp = plt.boxplot(filter_dist, 0, '', patch_artist=True)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    plt.setp(bp['boxes'], linewidth=2)
    plt.setp(bp['whiskers'], color='black', linewidth=2)
    plt.setp(bp['medians'], color='black', linewidth=2)
    plt.setp(bp['caps'], color='black', linewidth=2)
    pos = range(len(bp['boxes']))
    plt.xlabel('TE types')
    plt.ylabel('Distance (kbp)')
    plt.tight_layout()
    plt.savefig(distance_figure)
    plt.close()
    # plt.show()

if __name__ == "__main__":
    # Get set up
    input_file = os.path.abspath(args.input)
    count_figure = os.path.abspath(os.path.join(rundir, args.output+'_counts.png'))
    distance_figure = os.path.abspath(os.path.join(rundir, args.output+'_distances.png'))
    count_results = os.path.abspath(os.path.join(rundir, args.output+'_count_results.tsv'))
    distance_results = os.path.abspath(os.path.join(rundir, args.output+'_distance_results.tsv'))
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
        anno_file = os.path.abspath(args.annotations)
        with open(anno_file, 'r') as in_list:
            anno_list = [gene.strip().split('-T1')[0] for gene in in_list if gene != '\n']
        gene_errors = [anno for anno in anno_list if anno not in gene_list]
        if not gene_errors:
            pass
        if gene_errors:
            print('Warning: The following genes in annotation file are not found in gff3/bed file')
            print(gene_errors)

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
    print(gene_te_dist_dict)
    sys.exit()
    
    # Convert count dictionary to dataframe
    te_df = create_dataframe_from_count_dict(gene_te_count_dict)
    te_dist_df = create_dataframe_from_count_dict(gene_te_dist_dict)
    te_df.to_csv(count_results, sep='\t', header=True, index=False)
    te_dist_df.to_csv(distance_results, sep='\t', header=True, index=False)
    
    # Create figure
    colors = plt.cm.tab10([i for i in range(len(te_df.columns)-1)])
    create_bar_graphs(te_df)
    create_boxplot(te_dist_df)