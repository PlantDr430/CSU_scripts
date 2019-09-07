#!/usr/bin/python3

'''
This script follows formulas put forth in Kislyuk et al. (2011) and outputs bootstrap 
results and averaged overall results along with the fluidity graph over all genomes in 
the dataset, with     jackknifed variances fit to power law curves.
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
    description = '''    Performs a bootstrapping (resampling) of genome fluidity 
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
    '--bootstrap',
    default=5,
    type=int,
    help = 'Bootsraps [default: 5]',
    metavar=''
)
args=parser.parse_args()

if not os.path.isdir(args.out):
    os.makedirs(os.path.join(args.out, 'bootstraps'))
result_dir = os.path.abspath(os.path.join(rundir, args.out))
boot_dir = os.path.abspath(os.path.join(result_dir, 'bootstraps'))
if not os.path.isdir(boot_dir):
    os.makedirs(os.path.join(args.out, 'bootstraps'))

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
    # para = {'Shared' : 0, 'Uk' : 0, 'Ul' : 0}
    for k,v in ortho_isolates.items():
        if pair[0] in v and pair[1] in v:
            para['Shared'] += 1
        elif pair[0] in v and pair[1] not in v:
            para['Uk'] += 1
        elif pair[0] not in v and pair[1] in v:
            para['Ul'] += 1
    return

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

bootstrap_avg = {}
bootstrap_var = {}
for b in range(0, args.bootstrap): # bootstrap runs
    bootstrap_avg[b] = []
    genome_sample_dict = {}
    fluidity_dict = {}
    jack_sample_dict = {}
    jack_dict = {}
    variance_dict = {}
    for N in range(4, iso_num + 1): # main run
        genome_sample_dict[N] = []
        fluidity_dict[N] = []
        jack_sample_dict[N-1] = []
        jack_dict[N-1] = []
        variance_dict[N] = []
        random_pairs = []
        jack_pairs = []
        random_genome_sample = random.sample(iso_list, k=N) # random sample of genomes starting with a pool of 4
        get_pairs(random_genome_sample, random_pairs) # loop through random sample and create pairs of genomes
        jackknife_sample = random.sample(random_genome_sample, k=N-1) # randomly remove 1 genome from random sample 
        
        get_pairs(jackknife_sample, jack_pairs) # loop through estimator sample and create pairs of genomes

        for pair in random_pairs: # for pairs loop through gene cluster dictionary and pull out num of shared / single
            para = {'Shared' : 0, 'Uk' : 0, 'Ul' : 0}
            get_shared_single(pair)
            pair_fluidity = ((para['Uk'] + para['Ul'])/ # calculate fluidity per pair (Uk + Ul)/(Mk + Ml)
            ((para['Uk'] + para['Shared']) + (para['Ul'] + para['Shared'])))
            genome_sample_dict[N].append(pair_fluidity) # append all pair fluidities to dictionary for pool sample
        fluid = ((2/(N*(N-1)))*sum(genome_sample_dict[N])) # determine fluidity based on N genomes
        fluidity_dict[N].append((2/(N*(N-1)))*sum(genome_sample_dict[N]))

        for pair in jack_pairs: # for pairs loop through gene cluster dictionary and pull out num of shared / single
            para = {'Shared' : 0, 'Uk' : 0, 'Ul' : 0}
            get_shared_single(pair)
            jack_fluidity = ((para['Uk'] + para['Ul'])/ # calculate fluidity per pair (Uk + Ul)/(Mk + Ml)
            ((para['Uk'] + para['Shared']) + (para['Ul'] + para['Shared']))) # Make sure to reduce N by 1 genome
            jack_sample_dict[N-1].append(jack_fluidity) # append all pair fluidities to dictionary for pool sample
        fluid_i = ((2/(((N-1)-1)*((N-1)-2)))*sum(jack_sample_dict[N-1])) # calc fluid(i) for estimator sampled pairs
        fluid_i_mean = (1/(N-1))*fluid_i # calc fluid estimated mean from fluid(i)
        fluid_var = (((N-1)-1)/(N-1))*((fluid_i - fluid_i_mean)**2) # calc fluid variance 
        variance_dict[N].append(fluid_var)

    average = statistics.mean([float(''.join(map(str,v))) for v in fluidity_dict.values()]) # get avg fluidity for all
    var = [float(''.join(map(str,v))) for v in variance_dict.values()] # get variance for each N genome pool
    bootstrap_avg[b].append(average) # add averages ro bootstrap dictionary
    bootstrap_var[b] = var # add variances to bootstrap dictionary
    boot_file = os.path.abspath(os.path.join(boot_dir, 'bootstrap_'+str(b)+'_results.txt'))
    output = []
    with open(boot_file, 'w') as boot_out:
        boot_out.write('Genomes_Sampled\tFluidity\tVariance\n')
        for k, v in fluidity_dict.items():
            output.append([str(k), ''.join(map(str,v)), ''.join(map(str,variance_dict[k]))])
        for line in output:
            boot_out.write('\t'.join(line) + '\n')

fluid_fig = os.path.abspath(os.path.join(result_dir, 'Fluidity_chart.png'))
fig_results = os.path.abspath(os.path.join(result_dir, 'Fluidity.txt'))
overall_avg = statistics.mean([float(''.join(map(str,v))) for v in bootstrap_avg.values()]) # get avg across all boots
overall_data = np.array([overall_avg for i in range(4, iso_num + 1)]) # create y-value for all x-labels
var_list = [v for v in bootstrap_var.values()] # turn bootsrap varince dict into list
trans_var = np.array(var_list).T.tolist() # transpose list
final_var = np.array([statistics.mean(x) for x in trans_var]) # avg bootstrap variance across respective genome pools
x_labels = np.array([i for i in range(4, iso_num + 1)]) # get x-axis label

var_down = np.array([(overall_avg - v) for v in final_var]) # create array for y-values + variance
var_up = np.array([(overall_avg + v) for v in final_var]) # create array for y-values - variance

def powerlaw(x, m, c, c0):
    return c0 + x**m * c

popt_down, pcov_down = curve_fit(powerlaw, x_labels, var_down, p0 = np.asarray([-1,iso_num**5,0])) # power law var_down
popt_up, pcov_up = curve_fit(powerlaw, x_labels, var_up, p0 = np.asarray([-1,iso_num**5,0])) # power law var_up

fig, ax = plt.subplots()
plt.plot(x_labels, overall_data, ls='--', lw=1, color='black') # plot y-values of fluidity
plt.fill_between(x_labels, powerlaw(x_labels, *popt_up), powerlaw(x_labels, *popt_down), facecolor='cyan') # plot power
plt.xticks(np.arange(x_labels[0], x_labels[len(x_labels)-1]+1, 1.0)) # make sure x interval is 1
plt.xlim(x_labels[0], x_labels[len(x_labels)-1]) # adjust x limit so it starts with 4 at 0
plt.ylim((min(overall_data) - min(overall_data)*0.25), (max(overall_data) + max(overall_data)*0.25))
plt.xlabel('Genomes sampled')
plt.ylabel('Fluidity, '+u'\u03C6')
plt.tight_layout()
# plt.show()
plt.savefig(fluid_fig)

with open(fig_results, 'w') as results:
    results.write('Genomes_Sampled\tAvg_Fluidity\tAvg_Variance\tPower_up\tPower_down\n')
    r_out = []
    for i in range(0, iso_num-3):
        r_out.append([str(i+4), str(overall_avg), str(final_var.tolist()[i]), str(powerlaw(x_labels, *popt_up).tolist()[i]), 
        str(powerlaw(x_labels, *popt_down).tolist()[i])])
    for line in r_out:
        results.write('\t'.join(line) + '\n')