#!/usr/bin/python3

'''
This script follows formulas put forth in Kislyuk et al. (2011) to calculate genome 
fluidity of a pangenome dataset. Variance and standard error are estimated as total 
variance by averaging the variance calculated from bootstrapping samples of N genomes 
from the total pool of genomes and the variance calculated by the jackknife estimator of 
N-1 for each sample of N genomes. Results are a text file of fluidity, variance, and 
standard error for all N genome samples and a figure of total genome pool fluidity 
with shaded regions showing total standard error with a power law regression fit.

Notes
1. This will only work if you have at least 5 isolates to make up your pangenome.
2. If you have 5 isolates your graph will probably not look pretty as it's difficult 
to fit with such a low number of samples.
'''

import os, sys, re, argparse, random, itertools, scipy, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations
from collections import OrderedDict
from scipy.optimize import curve_fit, differential_evolution

rundir = os.getcwd()

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='./%(prog)s [options] -i orthogroups -o output_folder',
    description = '''    Performs multiple bootstraps and calculates genome fluidity 
    from a pangenome dataset.''',
    
    epilog = """Written by Stephen A. Wyka (2019)""",
    formatter_class = MyFormatter)

parser.add_argument(
    '-i',
    '--input',
    required = True,
    help = 'Orthogroups file, see format in READ.me',
    metavar=''
)
parser.add_argument(
    '-o',
    '--out',
    required = True,
    help = 'Output folder',
    metavar=''
)
parser.add_argument(
    '-p',
    '--prefix',
    help = 'Prefix to append to the result files (such as Genus, species, etc.)',
    metavar=''
)
args=parser.parse_args()

if not os.path.isdir(args.out):
    os.makedirs(os.path.join(args.out))
result_dir = os.path.abspath(os.path.join(rundir, args.out))

if args.input:
    input_file = os.path.abspath(args.input)
else:
    print('ERROR: No orthogroups file was provided please provide on, -i or --input')
    sys.exit()

if args.prefix:
    fluid_results = os.path.abspath(os.path.join(result_dir, args.prefix+'_pangenome_fluidity.txt'))
    fluid_fig = os.path.abspath(os.path.join(result_dir, args.prefix+'_pangenome_fluidity.png'))
else:
    fluid_results = os.path.abspath(os.path.join(result_dir, 'Pangenome_fluidity.txt'))
    fluid_fig = os.path.abspath(os.path.join(result_dir, 'Pangenome_fluidity.png'))

def create_ortho_dictionary(ortho_file): # create dictionary of gene clusters and isolates per cluster
    print('Creating ortholog dictionary')
    ortho_isolates_dict = OrderedDict() # {Protein Cluster : list of isolates represented in cluster}
    with open(ortho_file, 'r') as infile:
        ortho_list = [item.strip() for item in sorted(infile)]
        for line in ortho_list:
            iso_list = []
            if ':' in line:
                cluster, genes = line.split(':')
            elif '\t' in line:
                cluster, genes = line.split('\t', 1)
            else:
                cluster, genes = line.split(' ', 1)
            for match in re.finditer(r'([^\s]+)', genes):
                isolate = match.group(0).split('_')[0]
                iso_list.append(isolate)
            ortho_isolates_dict[cluster] = list(set(iso_list))
    return ortho_isolates_dict

def create_pair_dictionary(ortho_dictionary, isolate_list):
    print('Creating dictionary of paired ratio values')
    pair_dict = {} # {(Isolate1, Isolate2) : [ratio of sum(unique clusters)/sum(all clusters)]}
    for i in range(0, len(isolate_list)):
        for x in range(0, len(isolate_list)):
            if not isolate_list[i] == isolate_list[x]:
                pair = tuple(sorted([isolate_list[i], isolate_list[x]]))
                if not pair in pair_dict.keys():
                    cogs = {'Shared' : 0, 'Uk' : 0, 'Ul' : 0}
                    for k,v in ortho_dictionary.items():
                        if pair[0] in v and pair[1] in v:
                            cogs['Shared'] += 1
                        elif pair[0] in v and pair[1] not in v:
                            cogs['Uk'] += 1
                        elif pair[0] not in v and pair[1] in v:
                            cogs['Ul'] += 1
                        else:
                            pass # don't need to count a cluster if both isolates are not present
                    unique_pair = cogs['Uk'] + cogs['Ul']
                    all_pair = (cogs['Uk'] + cogs['Shared']) + (cogs['Ul'] + cogs['Shared'])
                    pair_dict[pair] = unique_pair/all_pair
    return pair_dict

def compute_fluidity_all_genomes(pair_dictionary, iso_num):
    fluidity_list = [ratio for ratio in pair_dictionary.values()]
    pangenome_fluidity = (2/(iso_num*(iso_num-1)))*sum(fluidity_list)
    return pangenome_fluidity

def jackknife_genome_subsamples(pair_dictionary, iso_num, iso_list):
    jack_var_dict = {} # {N genomes sampled : jackknife variance}
    for N in range(3, iso_num + 1):
        print('Performing jackknife permutations of {} genomes sampled'.format(N))
        jack_var_dict[N] = []
        N_combos = list(combinations(iso_list, N))
        for sample in N_combos:
            jack_sample = list(combinations(sample, N-1))
            fluidity_i_list = []
            for j_sample in jack_sample:
                pairs = tuple(combinations(j_sample,2))
                j_sample_fluidity = [pair_dictionary[tuple(sorted(p))] for p in pairs]
                fluidity_i = (2/((N-1)*(N-2)))*sum(j_sample_fluidity)
                fluidity_i_list.append(fluidity_i)
            fluidity_i_mean = np.mean(fluidity_i_list)
            fluidity_variance = ((N-1)/N)*sum([(i-fluidity_i_mean)**2 for i in fluidity_i_list])
            jack_var_dict[N].append(fluidity_variance)
    return jack_var_dict

def powerlaw(x, a, b, c):
    return a*(x**b) + c

def sumOfSquaredError(parameterTuple, x_values, y_curve_values):
    warnings.filterwarnings("ignore") # do not print warnings by genetic algorithm
    val = powerlaw(x_values, *parameterTuple)
    return np.sum((y_curve_values - val) ** 2.0)

def generate_Initial_Parameters(x_values, y_curve_values,):
    # min and max used for bounds
    maxX = max(x_values)
    minX = min(x_values)
    maxY = max(y_curve_values)
    minY = min(y_curve_values)
    maxXY = max(maxX, maxY)

    parameterBounds = []
    parameterBounds.append([-maxXY, maxXY]) # seach bounds for a
    parameterBounds.append([-maxXY, maxXY]) # seach bounds for b
    parameterBounds.append([-maxXY, maxXY]) # seach bounds for c

    # "seed" the numpy random number generator for repeatable results
    result = differential_evolution(sumOfSquaredError, parameterBounds, args=(x_values,y_curve_values), seed=3)
    return result.x

def create_fluidity_results(variance_dictionary, pangenome_fluidity, iso_num, figure_output, results_output):
    total_variance = np.array([np.mean(variance_dictionary[i]) for i in range(3, isolate_num + 1)])
    total_stderr = np.array([x**(1/2) for x in total_variance])
    y_fluidity_values = np.array([pangenome_fluidity for i in range(3, isolate_num + 1)])
    x_labels = np.array([i for i in range(3, isolate_num + 1)])
    stderr_bottom = np.array([(pangenome_fluidity - v) for v in total_stderr])
    stderr_top = np.array([(pangenome_fluidity + v) for v in total_stderr])
    geneticParameters_top = generate_Initial_Parameters(x_labels, stderr_top)
    geneticParameters_bottom = generate_Initial_Parameters(x_labels, stderr_bottom)
    popt_t, pcov = curve_fit(powerlaw, x_labels, stderr_top, geneticParameters_top, maxfev=10000)
    popt_b, pcov = curve_fit(powerlaw, x_labels, stderr_bottom, geneticParameters_bottom, maxfev=10000)
    fig, ax = plt.subplots()
    ax.set_axisbelow(True)
    plt.minorticks_on()
    plt.grid(which='minor', axis='y', color='white', linestyle='--', alpha=0.3)
    ax.yaxis.grid(True, linestyle='-', linewidth='1', which='major', color='white')
    ax.xaxis.grid(True, linestyle='-', linewidth='1', which='major', color='white', alpha=0.5)
    ax.tick_params(axis='x', which='minor', bottom=False)
    ax.set_facecolor('gainsboro')
    plt.plot(x_labels, y_fluidity_values, ls='--', lw=1, color='black') # plot y-values of fluidity
    plt.fill_between(x_labels, powerlaw(x_labels, *popt_t), powerlaw(x_labels, *popt_b), facecolor='blue', alpha=0.6)
    plt.xticks(np.arange(x_labels[0], x_labels[len(x_labels)-1]+1, 1.0)) # make sure x interval is 1
    plt.xlim(x_labels[0], x_labels[len(x_labels)-1]) # adjust x limit so it starts with 4 at 0
    max_y = max(stderr_top)
    min_y = min(stderr_bottom)
    plt.ylim((min_y - min_y*0.15), (max_y + max_y*0.15))
    plt.xlabel('Number of genomes sampled')
    plt.ylabel('Fluidity, '+u'\u03C6')
    plt.tight_layout()
    plt.savefig(figure_output)

    with open(results_output, 'w') as results: # print out fluidity results
        results.write('Genomes_Sampled\tFluidity\tTotal_Variance\tTotal_Stderr\tPower_top\tPower_bottom\n')
        r_out = []
        for i in range(0, iso_num-2):
            r_out.append([str(i+3), str(pangenome_fluidity), str(total_variance[i]), str(total_stderr[i]), 
            str(powerlaw(x_labels, *popt_t)[i]), str(powerlaw(x_labels, *popt_b)[i])])
        for line in r_out:
            results.write('\t'.join(line) + '\n')

if __name__ == "__main__":
    ortho_dict = create_ortho_dictionary(input_file)
    isolate_num = max([len(v) for v in ortho_dict.values()])
    isolate_list = list(set(itertools.chain.from_iterable([v for v in ortho_dict.values() if len(v) == isolate_num])))
    pair_dict = create_pair_dictionary(ortho_dict, isolate_list)
    pan_fluidity = compute_fluidity_all_genomes(pair_dict, isolate_num)
    jack_varaince_dict = jackknife_genome_subsamples(pair_dict, isolate_num, isolate_list)
    create_fluidity_results(jack_varaince_dict, pan_fluidity, isolate_num, fluid_fig, fluid_results)
