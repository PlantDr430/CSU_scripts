#!/usr/bin/python3

'''
This script follows formulas put forth in Kislyuk et al. (2011) to calculate genome 
fluidity of a pangenome dataset. Variance and standard error are estimated as total 
variance containing both the variance due to subsampling all possible combinations 
(without replacement) of N genomes from the total pool of genomes and the variance 
due to the limited number of sampled genomes (variance of the pangenome)(Kislyuk et al. 2011). 
However, the script has a default max number of subsamples set to 250,000 for each N genomes. 
This can be altered with the -max_sub / --max_subsamples flag or turned off with the --max_off flag. 
Turning the max_off will force calculations to be done on all possible subsample combinations 
of N genomes. For samples of N genomes that were stopped at the max number of subsamples the subsamples 
are sampled WITH replacement and variance is calculated with a degree of freedom = 1 (i.e. n - 1). 
Results are a text file of fluidity, variance, and standard error for all N genome samples 
and a figure of pangenome fluidity with shaded regions showing total standard error with a 
exponential regression fit.

Notes
1. This will only work if you have at least 5 isolates to make up your pangenome.
2. If you have 5 isolates your graph will probably not look pretty as it's difficult 
to fit with such a low number of samples.
'''

import os, sys, re, argparse, random, itertools, scipy, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from multiprocessing import Pool
from itertools import combinations
from collections import OrderedDict
from collections.abc import Iterable
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
    '-c',
    '--cpus',
    type=int,
    default=1,
    help = 'Number of cores to use for multiprocessing [default: 1]',
    metavar=''
)
parser.add_argument(
    '-max_sub',
    '--max_subsamples',
    type=int,
    default=250000,
    help = 'Max number of subsamples to run on N genomes sampled. [default: 250000]',
    metavar=''
)
parser.add_argument(
    '--max_off',
    action='store_true',
    help = 'Turn off the max subsamples. This will cause the script sample ALL possible combinations'\
    'for N genomes',
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

def create_pair_dictionary(ortho_dictionary):
    print('Creating dictionary of paired ratio values')
    pair_dict = {} # {(Isolate1, Isolate2) : [ratio of sum(unique clusters)/sum(all clusters)]}
    for i in range(0, len(iso_list)):
        for x in range(0, len(iso_list)):
            if not iso_list[i] == iso_list[x]:
                pair = tuple(sorted([iso_list[i], iso_list[x]]))
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

def compute_fluidity_all_genomes():
    '''
    Computes the fluidity and variance for the pangenome in question from the max number 
    of genomes in the pangenome.
    fluidity = 
    fluidity_i =
    fluidity_variance = 
    '''
    N = iso_num
    fluidity_list = [ratio for ratio in pair_dict.values()] # list of ratios 
    pangenome_fluidity = (2/(N*(N-1)))*sum(fluidity_list) # get fluidity from average of all ratios
    jack_samples = list(combinations(iso_list, N - 1)) # get list of all combos of N-1 from max num of genomes
    fluidity_i_list = []
    for sample in jack_samples:
        jack_pairs = tuple(combinations(sample,2)) # get all pairs from current jackknife sample
        jack_sample_fluidity = [pair_dict[tuple(sorted(p))] for p in jack_pairs] # get ratios from pair_dict
        fluidity_i = (2/((N-1)*(N-2)))*sum(jack_sample_fluidity) # calculate fluidity_i 
        fluidity_i_list.append(fluidity_i)
    fluidity_i_mean = np.mean(fluidity_i_list) # calculate fluidity_i_mean from all fluidity_i's
    fluidity_variance = ((N-1)/N)*sum([(i-fluidity_i_mean)**2 for i in fluidity_i_list]) # calculate variance
    return pangenome_fluidity, fluidity_variance

def subsample_multiprocess(combo_list):
    '''
    Takes portions of the full combo_list and runs them on separate threads for faster processing. 
    Calcualtes fluidity for each sample and returns list of fluidities.
    '''
    N = len(combo_list[0]) # get N from number of genomes present
    sample_process_list = []
    for sample in combo_list:
        pairs = tuple(combinations(sample,2))
        pair_fluidity_list = [pair_dict[tuple(sorted(p))] for p in pairs]
        sample_fluidity = (2/(N*(N-1)))*sum(pair_fluidity_list)
        sample_process_list.append(sample_fluidity)
    return sample_process_list

def genome_subsamples_fluidities(perm_list):
    '''
    Compute fluidities from all possible combinations of genomes from 3 to N randomly sampled genomes 
    (N is the max number of gneomes in sample, so only sampled once). Has a cut off of max subsamples 
    at which point variances are calcualted as sample variances (n-1) instead of full population 
    variances.
    '''
    sub_fluid_dict = {} # {N genomes sampled : [list of fluidities from subsamples]}
    for N in range(3, iso_num + 1):
        sub_fluid_dict[N] = []
        N_combos = list(combinations(iso_list, N))
        if args.max_off:
            combos = N_combos
        else:
            if len(N_combos) > args.max_subsamples:
                combos = random.choices(N_combos, k=args.max_subsamples)
                perm_list.append(N)
            else:
                combos = N_combos
        print('Performing fluidity calculations on {} subsample combinations of {} genomes'.format(len(combos),N))
        if not len(N_combos) == 1:
            chunk = round(len(combos)/args.cpus)
            split_combos = [combos[i:i + chunk] for i in range(0, len(combos), chunk)]
            pool = Pool(processes=args.cpus)
            results = pool.imap(subsample_multiprocess, split_combos)
            sub_fluid_dict[N].append(results)
        else:
            last_run = subsample_multiprocess(N_combos)
            sub_fluid_dict[N].append(last_run)
        sub_fluid_dict[N]=list(flatten(sub_fluid_dict[N]))
    return sub_fluid_dict

def flatten(lis):
    for item in lis:
     if isinstance(item, Iterable) and not isinstance(item, str):
         for x in flatten(item):
             yield x
     else:
         yield item

def exponential(x, a, b, c):
    return a * np.exp(b * x) + c

def neg_exponential(x, a, b, c):
    return a * np.exp(-b * x) + c

def sumOfSquaredError(parameterTuple, x_values, y_curve_values, func):
    warnings.filterwarnings("ignore") # do not print warnings by genetic algorithm
    val = func(x_values, *parameterTuple)
    return np.sum((y_curve_values - val) ** 2.0)

def generate_Initial_Parameters(x_values, y_curve_values, func):
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
    result = differential_evolution(sumOfSquaredError, parameterBounds, args=(x_values,y_curve_values, func), seed=3)
    return result.x

def create_fluidity_results(figure_output, results_output):
    total_variance = []
    for i in range(3, iso_num + 1):
        if i in permutation_list:
            total_variance.append(np.var(sub_fluid_dict[i], ddof = 1) + pan_variance)
        else:
            total_variance.append(np.var(sub_fluid_dict[i]) + pan_variance)
    total_variance = np.array(total_variance)
    total_stderr = np.array([x**(1/2) for x in total_variance])
    y_fluidity_values = np.array([pan_fluidity for i in range(3, iso_num + 1)])
    x_labels = np.array([i for i in range(3, iso_num + 1)])
    stderr_bottom = np.array([(pan_fluidity - v) for v in total_stderr])
    stderr_top = np.array([(pan_fluidity + v) for v in total_stderr])
    fig, ax = plt.subplots()
    try: # Still got problems sometimes with fitting curves, this temporary solution seems to be working
        geneticParameters_top = generate_Initial_Parameters(x_labels, stderr_top, exponential)
        geneticParameters_bottom = generate_Initial_Parameters(x_labels, stderr_bottom, exponential)
        popt_t, pcov = curve_fit(exponential, x_labels, stderr_top, geneticParameters_top, maxfev=10000)
        popt_b, pcov = curve_fit(exponential, x_labels, stderr_bottom, geneticParameters_bottom, maxfev=10000)
        if len(set(exponential(x_labels, *popt_t))) > 3 and len(set(exponential(x_labels, *popt_b))) > 3:
            plt.fill_between(x_labels, exponential(x_labels, *popt_t), exponential(x_labels, *popt_b), facecolor='blue', alpha=0.6)
            top_curve = exponential(x_labels, *popt_t)
            bottom_curve = exponential(x_labels, *popt_b)
        if len(set(exponential(x_labels, *popt_t))) <= 3:
            geneticParameters_top = generate_Initial_Parameters(x_labels, stderr_top, neg_exponential)
            popt_t, pcov = curve_fit(neg_exponential, x_labels, stderr_top, geneticParameters_top, maxfev=10000)
            plt.fill_between(x_labels, neg_exponential(x_labels, *popt_t), exponential(x_labels, *popt_b), facecolor='blue', alpha=0.6)
            top_curve = neg_exponential(x_labels, *popt_t)
            bottom_curve = exponential(x_labels, *popt_b)
        else:
            pass
        if len(set(exponential(x_labels, *popt_b))) <= 3:
            geneticParameters_bottom = generate_Initial_Parameters(x_labels, stderr_bottom, neg_exponential)
            popt_b, pcov = curve_fit(neg_exponential, x_labels, stderr_bottom, geneticParameters_bottom, maxfev=10000)
            plt.fill_between(x_labels, exponential(x_labels, *popt_t), neg_exponential(x_labels, *popt_b), facecolor='blue', alpha=0.6)
            top_curve = exponential(x_labels, *popt_t)
            bottom_curve = neg_exponential(x_labels, *popt_b)
        else:
            pass
    except:
        pass
    ax.set_axisbelow(True)
    plt.minorticks_on()
    plt.grid(which='minor', axis='y', color='white', linestyle='--', alpha=0.3)
    ax.yaxis.grid(True, linestyle='-', linewidth='1', which='major', color='white')
    ax.xaxis.grid(True, linestyle='-', linewidth='1', which='major', color='white', alpha=0.5)
    ax.tick_params(axis='x', which='minor', bottom=False)
    ax.set_facecolor('gainsboro')
    plt.plot(x_labels, y_fluidity_values, ls='--', lw=1, color='black') # plot y-values of fluidity
    plt.xticks(np.arange(x_labels[0], x_labels[len(x_labels)-1]+1, 1.0)) # make sure x interval is 1
    plt.xlim(x_labels[0], x_labels[len(x_labels)-1]) # adjust x limit so it starts with 3 at 0
    max_y = max(stderr_top)
    min_y = min(stderr_bottom)
    plt.ylim((min_y - min_y*0.15), (max_y + max_y*0.15))
    plt.xlabel('Number of genomes sampled')
    plt.ylabel('Fluidity, '+u'\u03C6')
    plt.tight_layout()
    # plt.show()
    plt.savefig(figure_output)

    with open(results_output, 'w') as results: # print out fluidity results
        results.write('Genomes_Sampled\tFluidity\tTotal_Variance\tTotal_Stderr\tExponential_top\tExponential_bottom\n')
        r_out = []
        for i in range(0, iso_num-2):
            r_out.append([str(i+3), str(pan_fluidity), str(total_variance[i]), str(total_stderr[i]), 
            str(top_curve[i]), str(bottom_curve[i])])
        for line in r_out:
            results.write('\t'.join(line) + '\n')

if __name__ == "__main__":
    ortho_dict = create_ortho_dictionary(input_file)
    iso_num = max([len(v) for v in ortho_dict.values()])
    iso_list = list(set(itertools.chain.from_iterable([v for v in ortho_dict.values() if len(v) == iso_num])))
    pair_dict = create_pair_dictionary(ortho_dict)
    pan_results = compute_fluidity_all_genomes()
    pan_fluidity = pan_results[0]
    pan_variance = pan_results[1]
    permutation_list = []
    sub_fluid_dict = genome_subsamples_fluidities(permutation_list)
    create_fluidity_results(fluid_fig, fluid_results)
