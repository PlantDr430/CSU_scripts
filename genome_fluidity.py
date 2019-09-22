#!/usr/bin/python3

'''
This script follows formulas put forth in Kislyuk et al. (2011) to calculate genome 
fluidity of a pangenome dataset. Variance and standard error are estimated as total 
variance by averaging the variance calculated from bootstrapping samples of N genomes 
from the total pool of genomes and the variance calculated by the jackknife estimator of 
N-1 for each sample of N genomes. Results are a text file of fluidity, variance, and 
standard error for all N genome samples and a figure of total genome pool fluidity 
with shaded regions showing total standard error with a power law regression fit.

Warnings!
1. This will only work if you have at least 5 isolates to make up your pangenome.
2. If you have 5 isolates your graph will probably not look pretty as it's difficult 
to fit with such a low number of samples. You could try reducing the bootstrapping.
3. If you get OptimizeWarning, you can disregard. It is just a warning and curves were 
correctly fit to the data.
'''

import os, sys, re, argparse, random, itertools, scipy, warnings
from collections import OrderedDict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution

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
    '-b',
    '--bootstraps',
    default=25,
    type=int,
    help = 'Number of bootstraps to run [default: 25]',
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

def get_shared_unique(pairs):
    for k,v in ortho_isolates.items():
            if pairs[0] in v and pairs[1] in v:
                para['Shared'] += 1
            elif pairs[0] in v and pairs[1] not in v:
                para['Uk'] += 1
            elif pairs[0] not in v and pairs[1] in v:
                para['Ul'] += 1
    return

def powerlaw(x, a, b, c):
    return a*(x**b) + c

def sumOfSquaredError(parameterTuple, x_values, y_curve_values):
    warnings.filterwarnings("ignore") # do not print warnings by genetic algorithm
    val = powerlaw(x_values, *parameterTuple)
    return np.sum((y_curve_values - val) ** 2.0)

def generate_Initial_Parameters(x_values, y_curve_values,):
    # min and max used for bounds
    maxX = max(x_labels)
    minX = min(x_labels)
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

# create dictionary of gene clusters and isolates per cluster
print('Creating dictionaries and matrices')
ortho_isolates = OrderedDict() # {Protein Cluster : list of isolates represented in cluster}
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

# create matrix dictionary of pairs for faster searching
matrix_dict = {} # {(Isolate1, Isolate2) : [unique cluster sum: all cluster sum]}
for i in range(0, len(iso_list)):
    for x in range(0, len(iso_list)):
        if not iso_list[i] == iso_list[x]:
            pair = tuple(sorted([iso_list[i], iso_list[x]]))
            if not pair in matrix_dict.keys():
                para = {'Shared' : 0, 'Uk' : 0, 'Ul' : 0}
                get_shared_unique(pair)
                unique_pair = para['Uk'] + para['Ul']
                all_pair = (para['Uk'] + para['Shared']) + (para['Ul'] + para['Shared'])
                matrix_dict[pair] = [unique_pair, all_pair]

boot_dict = {} # {bootstrap # : nested list of fluidity for [N -> iso_num] for all bootstraps}
jack_boot_var = {} # {bootstrap # : nested list of jackknife variances for [N -> iso_num] for all bootstraps}
jack_boot_err = {}# {bootstrap # : nested list of jackknife standard errors for [N -> iso_num] for all bootstraps}
for b in range(0, args.bootstraps): # bootstraps
    print('Running bootstrap {} out of {}'.format(b,args.bootstraps))
    genome_sample_dict = {} # {N genomes sampled : list of 1st part of fluidity calc for all pairs}
    jack_sample_dict = {} # {N genomes sampled : list of 1st part of fluidity calc for all pairs (N-1)}
    fluid_i_dict = {} # {N genomes sampled : list of fluid(i) calcs for each jackknife pair}
    jack_variance_dict = {} # {N genomes sampled : jackknife variance}
    jack_stderr_dict = {} # {N genomes sampled : jackknife standard error}
    for N in range(3, iso_num + 1): # main run
        genome_sample_dict[N] = []
        fluid_i_dict[N] = []
        jack_variance_dict[N] = []
        jack_stderr_dict[N] = []
        random_pairs = []

        random_genome_sample = random.sample(iso_list, k=N) # random sample of genomes starting with a pool of 3
        get_pairs(random_genome_sample, random_pairs) # loop through random sample and create pairs of genomes

        for pair in random_pairs: # for pairs loop through gene cluster dictionary and pull out num of shared & unique
            pair = tuple(sorted(pair))
            unique_sum = matrix_dict[pair][0]
            all_sum = matrix_dict[pair][1]
            pair_fluidity = (unique_sum/all_sum) # calculate fluidity ratio per pair (Uk + Ul)/(Mk + Ml)
            genome_sample_dict[N].append(pair_fluidity) # append all pair fluidities to dictionary for pool sample
        fluid = ((2/(N*(N-1)))*sum(genome_sample_dict[N])) # determine fluidity based on N genomes
        boot_list = boot_dict.get(N, [])
        boot_dict[N] = boot_list + [fluid]

        for i in range(0, len(random_genome_sample)): # jackknife run for N-1 genomes sampled
            jack_pairs = []
            jack_tmp = []
            jackknife_sample = [n for n in random_genome_sample if n != random_genome_sample[i]]
            get_pairs(jackknife_sample, jack_pairs)
            for pair in jack_pairs: # for pairs loop through gene cluster dictionary and pull out num of shared & unique
                pair = tuple(sorted(pair))
                unique_sum = matrix_dict[pair][0]
                all_sum = matrix_dict[pair][1]

                jack_fluidity = (unique_sum/all_sum) # calculate fluidity ratio per pair (Uk + Ul)/(Mk + Ml)
                jack_tmp.append(jack_fluidity)
            fluid_i = (2/((N-1)*(N-2)))*sum(jack_tmp) # get fluid(i) for for given set of pairs
            fluid_i_dict[N].append(fluid_i)
        fluid_i_mean = np.mean(fluid_i_dict[N]) # calculate the mean of fluid(i) across all sets
        fluid_var = ((N-1)/N)*sum([(i-fluid_i_mean)**2 for i in fluid_i_dict[N]]) # calculate variance from each fluid(i)
        fluid_stderr = fluid_var**(1/2) # calculate estimated standard error from jackknife variance
        jack_variance_dict[N].append(fluid_var)
        jack_stderr_dict[N].append(fluid_stderr)

    # add values to boot dictionaries for each boot run
    jack_variance = [float(''.join(map(str,v))) for v in jack_variance_dict.values()] # get variances for each N genome pool
    jack_stderr = [float(''.join(map(str,v))) for v in jack_stderr_dict.values()] # get stderrs for each N genome pool
    jack_boot_var[b] = jack_variance
    jack_boot_err[b] = jack_stderr

# get bootstrap variance and sterr from bootstraps from main fluidiy runs of N genomes sampled
boot_var = []
boot_err = []
for N in range(3, iso_num + 1):
    boot_est_mean = np.mean(boot_dict[N])
    boot_variance = sum([(x - boot_est_mean)**2 / (args.bootstraps - 1) for x in boot_dict[N]])
    boot_stderr = boot_variance**(1/2)
    boot_var.append(boot_variance)
    boot_err.append(boot_stderr)

# transpose jack_boot dictionaries to lists of pooled variance means
jack_var = np.array([v for v in jack_boot_var.values()]).T.tolist()
final_jack_var = np.array([np.mean(x) for x in jack_var])
jack_err = np.array([v for v in jack_boot_err.values()]).T.tolist()
final_jack_err = np.array([np.mean(x) for x in jack_err])

# get total variance from bootstrap estimates and jackknife estimates
combined_var = np.array([(boot_var[x] + final_jack_var[x]) for x in range(0, iso_num-2)])
combined_stderr = np.array([(boot_err[x] + final_jack_err[x]) for x in range(0, iso_num-2)])

# create files paths and get fluidty data and x-axis labels
if args.prefix:
    fluid_results = os.path.abspath(os.path.join(result_dir, args.prefix+'_pangenome_fluidity.txt'))
    fluid_fig = os.path.abspath(os.path.join(result_dir, args.prefix+'_pangenome_fluidity.png'))
else:
    fluid_results = os.path.abspath(os.path.join(result_dir, 'Pangenome_fluidity.txt'))
    fluid_fig = os.path.abspath(os.path.join(result_dir, 'Pangenome_fluidity.png'))
best_fluidity = boot_dict[iso_num][0] # get fludity calculated from all genomes in dataset
overall_data = np.array([best_fluidity for i in range(3, iso_num + 1)]) # create y-value for all x_labels
x_labels = np.array([i for i in range(3, iso_num + 1)]) # get x-axis label

# get top and bottom values of fluidity +/- stderr
err_bottom = np.array([(best_fluidity - v) for v in combined_stderr])
err_top = np.array([(best_fluidity + v) for v in combined_stderr])

# generate initial parameter values
geneticParameters_top = generate_Initial_Parameters(x_labels, err_top)
geneticParameters_bottom = generate_Initial_Parameters(x_labels, err_bottom)

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
stderr_top = max(err_top)
stderr_bottom = min(err_bottom)
plt.ylim((stderr_bottom - stderr_bottom*0.15), (stderr_top + stderr_top*0.15))
plt.xlabel('Genomes sampled')
plt.ylabel('Fluidity, '+u'\u03C6')
plt.tight_layout()
plt.savefig(fluid_fig)

with open(fluid_results, 'w') as results: # print out fluidity results
    results.write('Genomes_Sampled\tFluidity\tTotal_Variance\tTotal_Stderr\tPower_top\tPower_bottom\n')
    r_out = []
    for i in range(0, iso_num-2):
        r_out.append([str(i+3), str(boot_dict[i+3][0]), str(combined_var[i]), str(combined_stderr[i]), 
        str(powerlaw(x_labels, *popt_t)[i]), str(powerlaw(x_labels, *popt_b)[i])])
    for line in r_out:
        results.write('\t'.join(line) + '\n')