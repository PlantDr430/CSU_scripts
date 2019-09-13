#!/usr/bin/python3

'''
This script follows formulas put forth in Kislyuk et al. (2011) to calculate genome 
fluidity of a pangenome dataset. Variance and standard deviation are estimated as total 
variance by averaging the variance calculated from bootstrapping samples of N genomes 
from the total pool of genomes and the variance calculated by the jackknife estimator of 
N-1 for each sample of N genomes. Results are a text file of fluidity, variance, and 
standard deviation for all N genome samples and a figure of total genome pool fluidity 
with shaded regions showing standard deviations with a power law regression fit.

Warnings!
1. This will only work if you have at least 5 isolates to make up your pangenome.
2. If you have 5 isolates your graph will probably not look pretty as it's difficult 
to fit with such a low number of samples. You could try reducing the bootstrapping.
3. If you get OptimizeWarning, you can disregard. It is just a warning and curves were 
correctly fit to the data.
'''

import os, sys, re, argparse, random, itertools, scipy
from collections import OrderedDict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution
import warnings

rundir = os.getcwd()

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='./%(prog)s [options] -i orthogroups -o output_folder',
    description = '''    Performs multiple simulations and calculates genome fluidity 
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
    '-b',
    '--bootstraps',
    default=25,
    type=int,
    help = 'Number of bootstraps to run [default: 25]',
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

def get_pairs(sample_list, output_list):
    remove_list = []
    for iso in sample_list:
        remove_list.append(iso)
        sub_list = sample_list.copy()
        for r in remove_list:
            sub_list.remove(r)
        for n in range(0, len(sub_list)):
            random_iso = random.sample(sub_list, k=1)
            sub_list.remove(''.join(random_iso))
            output_list.append([iso,''.join(random_iso)])

def get_shared_single(pairs):
    for k,v in ortho_isolates.items():
            if pairs[0] in v and pairs[1] in v:
                para['Shared'] += 1
            elif pairs[0] in v and pairs[1] not in v:
                para['Uk'] += 1
            elif pairs[0] not in v and pairs[1] in v:
                para['Ul'] += 1
    return

def powerlaw(x, c, m, c0):
    return c*(x**m) + c0

def sumOfSquaredError(parameterTuple):
    warnings.filterwarnings("ignore") # do not print warnings by genetic algorithm
    val = powerlaw(x_labels, *parameterTuple)
    return np.sum((err_top - val) ** 2.0)

def generate_Initial_Parameters(level):
    # min and max used for bounds
    maxX = max(x_labels)
    minX = min(x_labels)
    maxY = max(level)
    minY = min(level)
    maxXY = max(maxX, maxY)

    parameterBounds = []
    parameterBounds.append([-maxXY, maxXY]) # seach bounds for a
    parameterBounds.append([-maxXY, maxXY]) # seach bounds for b
    parameterBounds.append([-maxXY, maxXY]) # seach bounds for c

    # "seed" the numpy random number generator for repeatable results
    result = differential_evolution(sumOfSquaredError, parameterBounds, seed=3)
    return result.x

# create dictionary of gene clusters and isolates per cluster
ortho_isolates = OrderedDict()
with open(input_file, 'r') as infile:
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
        ortho_isolates[cluster] = list(set(iso_list))

# get number of isoaltes in sample and create a list of isolates
iso_num = max([len(v) for v in ortho_isolates.values()])
iso_list = list(set(itertools.chain.from_iterable([v for v in ortho_isolates.values() if len(v) == iso_num])))

matrix_dict = {}
for i in range(0, len(iso_list)):
    for x in range(0, len(iso_list)):
        if not iso_list[i] == iso_list[x]:
            pair = tuple(sorted([iso_list[i], iso_list[x]]))
            if not pair in matrix_dict.keys():
                para = {'Shared' : 0, 'Uk' : 0, 'Ul' : 0}
                get_shared_single(pair)
                unique_pair = para['Uk'] + para['Ul']
                all_pair = (para['Uk'] + para['Shared']) + (para['Ul'] + para['Shared'])
                matrix_dict[pair] = [unique_pair, all_pair]

boot_dict = {}
jack_boot_var = {}
jack_boot_dev = {}
for b in range(0, args.bootstraps): # bootstraps
    genome_sample_dict = {}
    fluidity_dict = {}
    jack_sample_dict = {}
    fluid_i_dict = {}
    jack_variance_dict = {}
    jack_stdev_dict = {}
    for N in range(3, iso_num + 1): # main run
        genome_sample_dict[N] = []
        fluidity_dict[N] = []
        fluid_i_dict[N] = []
        jack_variance_dict[N] = []
        jack_stdev_dict[N] = []
        random_pairs = []

        random_genome_sample = random.sample(iso_list, k=N) # random sample of genomes starting with a pool of 4
        get_pairs(random_genome_sample, random_pairs) # loop through random sample and create pairs of genomes

        for pair in random_pairs: # for pairs loop through gene cluster dictionary and pull out num of shared / single
            pair = tuple(sorted(pair))
            unique_sum = matrix_dict[pair][0]
            all_sum = matrix_dict[pair][1]
            pair_fluidity = (unique_sum/all_sum) # calculate fluidity per pair (Uk + Ul)/(Mk + Ml)
            genome_sample_dict[N].append(pair_fluidity) # append all pair fluidities to dictionary for pool sample
        fluid = ((2/(N*(N-1)))*sum(genome_sample_dict[N])) # determine fluidity based on N genomes
        fluidity_dict[N].append((2/(N*(N-1)))*sum(genome_sample_dict[N]))
        boot_list = boot_dict.get(N, [])
        boot_dict[N] = boot_list + [fluid]

        for i in range(0, len(random_genome_sample)):
            jack_pairs = []
            jack_tmp = []
            jackknife_sample = [n for n in random_genome_sample if n != random_genome_sample[i]]
            get_pairs(jackknife_sample, jack_pairs)
            for pair in jack_pairs: # for pairs loop through gene cluster dictionary and pull out num of shared / single
                pair = tuple(sorted(pair))
                unique_sum = matrix_dict[pair][0]
                all_sum = matrix_dict[pair][1]

                jack_fluidity = (unique_sum/all_sum) # calculate fluidity per pair (Uk + Ul)/(Mk + Ml)
                jack_tmp.append(jack_fluidity)
            fluid_i = (2/((N-1)*(N-2)))*sum(jack_tmp)
            fluid_i_dict[N].append(fluid_i)
        fluid_i_mean = np.mean(fluid_i_dict[N])
        fluid_var = ((N-1)/N)*sum([(i-fluid_i_mean)**2 for i in fluid_i_dict[N]])
        fluid_stdev = fluid_var**(1/2)
        jack_variance_dict[N].append(fluid_var)
        jack_stdev_dict[N].append(fluid_stdev)

    # add values to jack_boot dictionaries for each boot run
    jack_variance = [float(''.join(map(str,v))) for v in jack_variance_dict.values()] # get variance for each N genome pool
    jack_stdev = [float(''.join(map(str,v))) for v in jack_stdev_dict.values()] # get stderr for each N genome pool
    jack_boot_var[b] = jack_variance
    jack_boot_dev[b] = jack_stdev

# transpose jack_boot dictionaries to lists of means
jack_var = np.array([v for v in jack_boot_var.values()]).T.tolist()
final_jack_var = np.array([np.mean(x) for x in jack_var])
jack_dev = np.array([v for v in jack_boot_dev.values()]).T.tolist()
final_jack_dev = np.array([np.mean(x) for x in jack_dev])


boot_var = []
boot_stdev = []
for v in range(3, iso_num + 1):
    boot_var.append(np.var(boot_dict[v]))
    boot_stdev.append(np.std(boot_dict[v]))

# get total variance from bootstrap estimates and jackknife estimates
combined_var = np.array([(boot_var[x] + final_jack_var[x]) for x in range(0, iso_num-2)])
combined_stdev = np.array([(boot_stdev[x] + final_jack_dev[x]) for x in range(0, iso_num-2)])

# create files paths and get fluidty data and x-axis labels
fluid_results = os.path.abspath(os.path.join(result_dir, 'Pangenome_fluidity.txt'))
fluid_fig = os.path.abspath(os.path.join(result_dir, 'Pangenome_fluidity.png'))
best_fluidity = fluidity_dict[iso_num][0] # get fludity calculated from all genomes in dataset
overall_data = np.array([best_fluidity for i in range(3, iso_num + 1)]) # create y-value for all x_labels
x_labels = np.array([i for i in range(3, iso_num + 1)]) # get x-axis label

# get top and bottom values of fludidty +/- stdev
err_bottom = np.array([(best_fluidity - v) for v in combined_stdev])
err_top = np.array([(best_fluidity + v) for v in combined_stdev])

# generate initial parameter values
geneticParameters_top = generate_Initial_Parameters(err_top)
geneticParameters_bottom = generate_Initial_Parameters(err_bottom)

# get best fit model parameters. Used large maxfev to make it robust as fitting is 
# difficult with ~5 genomes, even with generated initial values. Since we are setting 
# initial values I don't see a concern of having a high maxfev.
popt_t, pcov = curve_fit(powerlaw, x_labels, err_top, geneticParameters_top, maxfev=10000)
popt_b, pcov = curve_fit(powerlaw, x_labels, err_bottom, geneticParameters_bottom, maxfev=10000)

fig, ax = plt.subplots()
ax.set_axisbelow(True)
plt.minorticks_on()
plt.grid(which='minor', axis='y', color='white', linestyle='--', alpha=0.3)
ax.yaxis.grid(True, linestyle='-', linewidth='1', which='major', color='white')
ax.xaxis.grid(True, linestyle='-', linewidth='1', which='major', color='white', alpha=0.5)
ax.tick_params(axis='x', which='minor', bottom=False)
ax.set_facecolor('gainsboro')
plt.plot(x_labels, overall_data, ls='--', lw=1, color='black') # plot y-values of fluidity
plt.fill_between(x_labels, powerlaw(x_labels, *popt_t), powerlaw(x_labels, *popt_b), facecolor='blue', alpha=0.6)
plt.xticks(np.arange(x_labels[0], x_labels[len(x_labels)-1]+1, 1.0)) # make sure x interval is 1
plt.xlim(x_labels[0], x_labels[len(x_labels)-1]) # adjust x limit so it starts with 4 at 0
stdev_top = max(err_top)
stdev_bottom = min(err_bottom)
plt.ylim((stdev_bottom - stdev_bottom*0.15), (stdev_top + stdev_top*0.15))
plt.xlabel('Genomes sampled')
plt.ylabel('Fluidity, '+u'\u03C6')
plt.tight_layout()
plt.savefig(fluid_fig)

with open(fluid_results, 'w') as results:
    results.write('Genomes_Sampled\tFluidity\tVariance\tStderr\tPower_top\tPower_bottom\n')
    r_out = []
    for i in range(0, iso_num-2):
        r_out.append([str(i+3), str(fluidity_dict[i+3][0]), str(combined_var[i]), str(combined_stdev[i]),
        str(powerlaw(x_labels, *popt_t)[i]), str(powerlaw(x_labels, *popt_b)[i])])
    for line in r_out:
        results.write('\t'.join(line) + '\n')