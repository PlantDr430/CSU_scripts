#!/usr/bin/python3

'''

NO headers on any of the files!

population_file = List of all gene_cluster ID's in the population (one cluster ID per line)
study_file = List of all gene_cluster ID's in the study (one cluster ID per line)
association_file = List of all gene_cluster ID's in the population within with associated domains:

Association file example (separated by a tab):

OG0006472	IPR023213
OG0006473	IPR003593,IPR003439,IPR011527,IPR036640,IPR017871,IPR027417
OG0006659	IPR004045,IPR004046,IPR010987,IPR036249,IPR036282
OG0006688	IPR021102
OG0006705	IPR020806,IPR023213,IPR010071,IPR036736,IPR000873,IPR001242,IPR009081,IPR006162,IPR020845
OG0006868	IPR002227,IPR008922
OG0006869	IPR002227,IPR008922
OG0006890	IPR002933,IPR017144,IPR011650,IPR017439,IPR036264
OG0006891	IPR020801,IPR020841,IPR020807,IPR020806,IPR020843,IPR014043,IPR014030,IPR013154,IPR032821,IPR001227,IPR013968,IPR009081,IPR014031,IPR036736,IPR013149,IPR016039,IPR006162,IPR018201,IPR036291,IPR016036,IPR011032,IPR016035
OG0006892	IPR029058,IPR005645
OG0006893	IPR021054
OG0006894	IPR021054
OG0006895	IPR018523,IPR015813
OG0006896	IPR010699
OG0006897	IPR004358,IPR001789,IPR003661,IPR003594,IPR005467,IPR036097,IPR011006,IPR036890
OG0006993	IPR004046,IPR004045,IPR010987,IPR036282,IPR036249
OG0006994	IPR011059,IPR013108,IPR032466
OG0007094	
OG0007158	IPR021858
OG0007159	IPR036259
OG0007160	
OG0007161	IPR021765
OG0007164	
OG0007165	IPR027443,IPR026992,IPR005123

'''

import os, sys, re, argparse
import pandas as pd
import numpy as np
from collections.abc import Iterable
from collections import OrderedDict
from statsmodels.stats.multitest import multipletests
import scipy.stats as stats

rundir = os.getcwd()

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='./%(prog)s [options] -p population -s study -a association -o output',
    description = '''    Simple enrichment analysis using fishcer exact test 
    and multi-test corrctions.''',
    
    epilog = """Written by Stephen A. Wyka (2020)""",
    formatter_class = MyFormatter)

parser.add_argument(
    '-p',
    '--population',
    required = True,
    help = 'Population file',
    metavar=''
)
parser.add_argument(
    '-s',
    '--study',
    required = True,
    help = 'Study file',
    metavar=''
)
parser.add_argument(
    '-a',
    '--association',
    required = True,
    help = 'Association file',
    metavar=''
)
parser.add_argument(
    '-o',
    '--out',
    required = True,
    help = 'Basename of output stats',
    metavar=''
)
parser.add_argument(
    '-f_p',
    '--fischer_pvalue',
    default = 0.05,
    type = float,
    help = 'Pvalue cutoff for Fischer exact test (annotation bar charts) [default: 0.05]',
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
args=parser.parse_args()

def flatten(lis):
    for item in lis:
     if isinstance(item, Iterable) and not isinstance(item, str):
         for x in flatten(item):
             yield x
     else:
         yield item

def create_association_dict(filename):
    assoc_df = pd.read_csv(filename, delimiter='\t', header=None)
    assoc_dict = {}
    for index, row in assoc_df.iterrows():
        cluster = row[0]
        functions = row[1]
        if not str(functions) == 'nan':
            funcs = functions.split(',')
            if len(funcs) == 1:
                if funcs[0] in assoc_dict.keys():
                    assoc_dict[funcs[0]].append(cluster)
                else:
                    assoc_dict[funcs[0]] = [cluster]
                assoc_dict[funcs[0]] = list(set(assoc_dict[funcs[0]]))
            else:
                for i in range(0,len(funcs)):
                    if funcs[i] in assoc_dict.keys():
                        assoc_dict[funcs[i]].append(cluster)
                    else:
                        assoc_dict[funcs[i]] = [cluster]
                    assoc_dict[funcs[i]] = list(set(assoc_dict[funcs[i]]))
        else:
            pass
    return assoc_dict

def get_ratios(dict, pop, study):
    ratio_dict = OrderedDict()
    for func, clusters in dict.items():
        pop_count = 0
        study_count = 0
        for clus in clusters:
            if clus in study:
                study_count += 1
                pop_count += 1
            elif clus in pop:
                pop_count += 1
            else:
                pass
        ratio_dict[func] = [[study_count, len(study)],[pop_count, len(pop)]]

    return ratio_dict

def run_fischers_test(r_dict,a_dict):
    tmp_dict = OrderedDict()
    uncorr_p_list = []
    for func, r in r_dict.items():
        o, p = stats.fisher_exact([[r[0][1],r[1][1]],[r[0][0],r[1][0]]])
        if p <= args.fischer_pvalue:
            uncorr_p_list.append(p)
            study_ratio = (r[0][0]/r[0][1])
            pop_ratio = (r[1][0]/r[1][1])
            if study_ratio > pop_ratio:
                tmp_dict[func] = ['e','{}/{}'.format(r[0][0],r[0][1]),
                '{}/{}'.format(r[1][0],r[1][1]), str(p)]
            elif pop_ratio > study_ratio:
                tmp_dict[func] = ['p','{}/{}'.format(r[0][0],r[0][1]),
                '{}/{}'.format(r[1][0],r[1][1]), str(p)]
    if uncorr_p_list:
        corr_p_list = run_multitest_correction(uncorr_p_list)
        results_dict = {}
        count = -1
        for func, results in tmp_dict.items():
            count += 1
            if corr_p_list[count] <= args.multialpha:
                new_line = list(flatten(results + [str(corr_p_list[count])])) + [','.join([x for x in a_dict[func] if x in study])]
                results_dict[func] = new_line
            else:
                pass

    else:
        results_dict = {}
        print('We did not find any signficiant enrichment or purification for {}.'.format(args.out))
    return results_dict

def run_multitest_correction(uncorr_list):
    r, p_b, aS, aB = multipletests(uncorr_list, alpha=args.multialpha, method='bonferroni', is_sorted=False, returnsorted=False)
    r, p_f, aS, aB = multipletests(uncorr_list, alpha=args.multialpha, method='fdr_bh', is_sorted=False, returnsorted=False)
    avg_list = []
    for i in range(0, len(uncorr_list)):
        avg_list.append(np.mean([p_b[i],p_f[i]]))
    return avg_list

def write_out_results(dict, filename):
    outfile = os.path.abspath(os.path.join(rundir, filename))
    with open(outfile, 'w') as output:
        output.write('Domain_ID\tEnrich/Purif\tratio_in_study\tratio_in_pop\tp_uncorrected\tp_multitest\tstudy_gene_families\n')
        for func, results in dict.items():
            output.write(func+'\t'+'\t'.join(results[:-1])+'\t'+results[-1]+'\n')

if __name__ == "__main__":
    output = os.path.abspath(args.out)
    with open(args.study, 'r') as st:
        study = [item.strip() for item in st]
    with open(args.population, 'r') as po:
        pop = [item.strip() for item in po]

    assoc_dict = create_association_dict(args.association)
    ratio_dict = get_ratios(assoc_dict, pop, study)
    results_dict = run_fischers_test(ratio_dict, assoc_dict)
    write_out_results(results_dict, output+'.tsv')

