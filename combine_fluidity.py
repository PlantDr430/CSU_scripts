#!/usr/bin/python3

import os, sys, argparse, csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy import sqrt, abs
from scipy.stats import norm

rundir = os.getcwd()

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='./%(prog)s [options] -i XXX_pangenome_fluidity.txt XXX_pangenome_fluidity.txt -o output_folder',
    description = '''    Turns multiple outputs from genome_fluidity.py into a figure and pvalue matrix. 
    Make sure to attach a prefix "XXX_pangenome_fluidity.txt if you didn't already use the -p / --prefix 
    flag in the previous script, to be used in the legend''',

    epilog = """Written by Stephen A. Wyka (2019)""",
    formatter_class = MyFormatter)

parser.add_argument(
    '-i',
    '--input',
    nargs = '+',
    required = True,
    help = 'All results.txt files from genome_fluidity.py that you wish to include',
    metavar=''
)
parser.add_argument(
    '-o',
    '--output',
    default = 'Combined_fluidity',
    help = 'Name of output figure (.png will be added automatically)',
    metavar=''
)
parser.add_argument(
    '-f',
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
    '-l',
    '--legend',
    default = 6,
    type=int,
    help = 'Size of legend [default = 6]',
    metavar=''
)
args=parser.parse_args()

x_len = []
data = []
curve_top = []
curve_bottom = []
legend_labels = []
z_test_dict = {} # {Name : [fluidity, genome count, variance]}
def parse_input_files(input_file):
    legend_name = input_file.split('_')[0]
    legend_labels.append(legend_name)
    fpath = os.path.abspath(os.path.join(rundir, input_file))
    df=pd.read_csv(fpath, delimiter='\t', header=0)
    x_labels = np.array(df['Genomes_Sampled'])
    genome_count = x_labels[-1]
    x_len.append(genome_count)
    variance = sqrt(df['Total_Variance'].iloc[len(df['Total_Variance'])-1])
    fluid_all = df['Fluidity'].iloc[len(df['Fluidity'])-1]
    z_test_dict[legend_name] = [fluid_all,genome_count,variance]
    fluid_data = np.array([fluid_all for x in range(0,len(df['Fluidity']))])
    data.append((x_labels, fluid_data))
    curve_top.append(np.array(df['Exponential_top']))
    curve_bottom.append(np.array(df['Exponential_bottom']))

def create_fluidity_figure(output_file):
    fig, ax = plt.subplots(figsize=(args.figsize[0],args.figsize[1]), dpi=args.dpi)
    colors = plt.cm.nipy_spectral([1.*i/len(args.input) for i in range(len(args.input))])
    for i in range(0, len(data)):
        num_genomes = str(data[i][0][-1])
        plt.plot(data[i][0], data[i][1], ls='--', lw=1.5, color=colors[i], alpha=1, label='{}'.format(legend_labels[i]+' ('+num_genomes+')'))
        plt.fill_between(data[i][0], curve_top[i], curve_bottom[i], facecolor=colors[i], alpha=0.25)
    ax.set_facecolor('lightgrey')
    ax.set_axisbelow(True)
    plt.minorticks_on()
    plt.grid(which='minor', axis='y', color='white', linestyle='--', alpha=0.3)
    ax.tick_params(axis='x', which='minor', bottom=False)
    ax.grid(which='major', linestyle='-', linewidth='1', color='white')
    plt.xticks(np.arange(3, max(x_len) + 1, 1.0))
    plt.xlim(3, max(x_len))
    plt.xlabel('Number of genomes sampled')
    plt.ylabel('Fluidity, '+u'\u03C6')
    plt.legend(framealpha=1.0, prop={'size': args.legend})
    plt.tight_layout()
    plt.savefig(output_file)

def two_sample_z(X1, n1, var1, X2, n2, var2):
    '''
    Modified z-test based on documentation in Kislyuk et al. 2011 and personal
    communications with the corresponding author Joshua S Weitz.
    '''
    pooledVAR = sqrt(var1 + var2)
    z = ((X1 - X2) - 0)/pooledVAR
    pval = 2*(norm.sf(abs(z)))
    return pval

def determine_significance_matrix(output_file):
    sorted_dict = sorted(z_test_dict.items(), key=lambda x: x[1], reverse=True)
    result_matrix = []
    for i in range(0, len(sorted_dict)): # loop through isolates
        line_lists = ['-'] * (i+1) # add blanks for easy matrix to datafram as we have unequal lengths
        for m in range(i+1, len(sorted_dict)): # compute pvalues for isolate v. others (no duplicates, no self)
            args = sorted_dict[i][1] + sorted_dict[m][1] # combine dict values to list for key v. key comparison
            pval = two_sample_z(args[0],args[1],args[2],args[3],args[4],args[5]) # get pvalue from two-sample z-test
            short_pval = np.format_float_scientific(pval, precision=2) # shorten pvalue singificant figures
            line_lists.append(short_pval) # create list of all pvalues for isolate i v. others
        result_matrix.append([sorted_dict[i][0]] + line_lists) # append all pvalues for isolate i to matrix
    labels = [sorted_dict[x][0] for x in range(0, len(sorted_dict))] # get headers and index for dataframe
    with open(output_file, 'w') as outfile: # would love to print out a table.png but matplotlib tables aren't good
        report = csv.writer(outfile, delimiter = '\t', lineterminator='\n')
        report.writerow([''] + labels)
        for i in range(0, len(sorted_dict)):
            report.writerow(result_matrix[i])

if __name__ == "__main__":
    output_figure = os.path.abspath(os.path.join(rundir, args.output+'.png'))
    output_pval_table = os.path.abspath(os.path.join(rundir, args.output+'_table.tsv'))
    for files in args.input:
        parse_input_files(files)
    create_fluidity_figure(output_figure)
    determine_significance_matrix(output_pval_table)