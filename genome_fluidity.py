#!/usr/bin/python3

'''
This script follows formulas put forth in Kislyuk et al. (2011) and performs multiple 
simulations and averages the simulations standard error estiamtes from jackknife 
variance estiamtes to output text results of fluidity, variance, and standard error 
as well as a figure.

If you get OptimizeWarning, you can disregard. It is just a wanring and curves were 
correctly fit to the data.
'''

import os, sys, re, argparse, random, itertools, statistics, math
from collections import OrderedDict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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
    '-s',
    '--simulations',
    default=10,
    type=int,
    help = 'Number of simulations to run [default: 10]',
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
        if pair[0] in v and pair[1] in v:
            para['Shared'] += 1
        elif pair[0] in v and pair[1] not in v:
            para['Uk'] += 1
        elif pair[0] not in v and pair[1] in v:
            para['Ul'] += 1
    return

def powerlaw(x, m, c, c0):
    return c*(x**m) + c0

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
        pair = (iso_list[i], iso_list[x])
        para = {'Shared' : 0, 'Uk' : 0, 'Ul' : 0}
        get_shared_single(pair)
        unique_pair = para['Uk'] + para['Ul']
        all_pair = (para['Uk'] + para['Shared']) + (para['Ul'] + para['Shared'])
        matrix_dict[pair] = [unique_pair, all_pair]

sim_var = {}
sim_err = {}
for b in range(0, args.simulations): # simulations
    genome_sample_dict = {}
    fluidity_dict = {}
    jack_sample_dict = {}
    fluid_i_dict = {}
    variance_dict = {}
    stderr_dict = {}
    for N in range(4, iso_num + 1): # main run
        print(N)
        genome_sample_dict[N] = []
        fluidity_dict[N] = []
        fluid_i_dict[N] = []
        variance_dict[N] = []
        stderr_dict[N] = []
        random_pairs = []
        
        random_genome_sample = random.sample(iso_list, k=N) # random sample of genomes starting with a pool of 4
        get_pairs(random_genome_sample, random_pairs) # loop through random sample and create pairs of genomes
        
        for pair in random_pairs: # for pairs loop through gene cluster dictionary and pull out num of shared / single
            pair = tuple(pair)
            unique_sum = matrix_dict[pair][0]
            all_sum = matrix_dict[pair][1]
            pair_fluidity = (unique_sum/all_sum) # calculate fluidity per pair (Uk + Ul)/(Mk + Ml)
            genome_sample_dict[N].append(pair_fluidity) # append all pair fluidities to dictionary for pool sample
        fluid = ((2/(N*(N-1)))*sum(genome_sample_dict[N])) # determine fluidity based on N genomes
        fluidity_dict[N].append((2/(N*(N-1)))*sum(genome_sample_dict[N]))

        for i in range(0, len(random_genome_sample)):
            jack_pairs = []
            jack_tmp = []
            jackknife_sample = [n for n in random_genome_sample if n != random_genome_sample[i]]
            get_pairs(jackknife_sample, jack_pairs)
            for pair in jack_pairs: # for pairs loop through gene cluster dictionary and pull out num of shared / single
                pair = tuple(pair)
                unique_sum = matrix_dict[pair][0]
                all_sum = matrix_dict[pair][1]

                jack_fluidity = (unique_sum/all_sum) # calculate fluidity per pair (Uk + Ul)/(Mk + Ml)
                jack_tmp.append(jack_fluidity)
            fluid_i = (2/((N-1)*(N-2)))*sum(jack_tmp)
            fluid_i_dict[N].append(fluid_i)
        fluid_i_mean = statistics.mean(fluid_i_dict[N])
        fluid_var = ((N-1)/N)*sum([(i-fluid_i_mean)**2 for i in fluid_i_dict[N]])
        fluid_stderr = fluid_var**(1/2)
        variance_dict[N].append(fluid_var)
        stderr_dict[N].append(fluid_stderr)

    var = [float(''.join(map(str,v))) for v in variance_dict.values()] # get variance for each N genome pool
    stderr = [float(''.join(map(str,v))) for v in stderr_dict.values()] # get stderr for each N genome pool
    sim_var[b] = var # add variances to simulation dictionary
    sim_err[b] = stderr # add stderr to simulation dictionary

fluid_results = os.path.abspath(os.path.join(result_dir, 'Pangenome_fluidity.txt'))
fluid_fig = os.path.abspath(os.path.join(result_dir, 'Pangenome_fluidity.png'))
best_fluidity = fluidity_dict[iso_num][0]
overall_data = np.array([best_fluidity for i in range(4, iso_num + 1)]) # create y-value for all x-labels
x_labels = np.array([i for i in range(4, iso_num + 1)]) # get x-axis label

var_list = [v for v in sim_var.values()] # turn simulation varince dict into list
trans_var = np.array(var_list).T.tolist() # transpose list
final_var = np.array([statistics.mean(x) for x in trans_var]) # get final array

err_list = [v for v in sim_err.values()] # turn simulation varince dict into list
trans_err = np.array(err_list).T.tolist() # transpose list
final_err = np.array([statistics.mean(x) for x in trans_err]) # get final array

err_bottom = np.array([(best_fluidity - v) for v in final_err]) # calculate y-values for stderr_top
err_top = np.array([(best_fluidity + v) for v in final_err]) # calculate y-values for stderr_top

popt_t, pcov_t = curve_fit(powerlaw, x_labels, err_top, p0 = np.asarray([0,1,1]), maxfev=2000) # powerlaw fit_top
popt_b, pcov_b = curve_fit(powerlaw, x_labels, err_bottom, p0 = np.asarray([0,1,1]), maxfev=2000) # powerlaw fit_bottom

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
stdev_top = max(powerlaw(x_labels, *popt_t))
stdev_bottom = min(powerlaw(x_labels, *popt_b))
plt.ylim((stdev_bottom - stdev_bottom*0.25), (stdev_top + stdev_top*0.25))
plt.xlabel('Genomes sampled')
plt.ylabel('Fluidity, '+u'\u03C6')
plt.tight_layout()
plt.savefig(fluid_fig)

with open(fluid_results, 'w') as results:
    results.write('Genomes_Sampled\tFluidity\tVariance\tStderr\n')
    r_out = []
    for i in range(0, iso_num-3):
        r_out.append([str(i+4), str(best_fluidity), str(final_var.tolist()[i]), str(final_err.tolist()[i])])
    for line in r_out:
        results.write('\t'.join(line) + '\n')