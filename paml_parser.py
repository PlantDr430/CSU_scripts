#!/usr/bin/python3

'''
Loops through a directory containing multiple paml output files, with a naming scheme 
of clusterID.XXXXXX (for example 2356.paml_out). Returns two text file:

1. A working file of Model output. (see below)

Cluster	Model	np	lnL	Number_sites_>95	Number_sites_>99
551	Model 0	-3309.989574
551	Model 1	-3294.362074
551	Model 2	-3280.01743	3	1
551	Model 3	-3279.439089
551	Model 7	-3294.376692
551	Model 8	-3280.047899	7	2
552	Model 0	-1038.052331
552	Model 1	-1009.671555
552	Model 2	-1008.532431	0	0
552	Model 3	-1008.532431
552	Model 7	-1009.761659
552	Model 8	-1008.532415	0	0

2. Sumarized text file for M3vM0, M2vM1, and M8vM7 showing LRT values, p_values, 
and number of BEB significant sites in the M2 and M8 model. (see below)

Cluster	M3vM0_lrt	M3vM0_p	M2vM1_lrt	M2vM1_p	M2_num_sites_>95	M2_num_sites_>99	M8vM7_lrt	M8vM7_p	M8_num_sites_>95	M8_num_sites_>99
551	61.1	0.0	28.69	0.0	3	1	28.66	0.0	7	2
552	59.04	0.0	2.28	0.32	0	0	2.46	0.293	0	0
553	-0.0	1.0	-0.0	1.0	0	0	-0.0	1.0	0	0
554	0.0	1.0	0.0	1.0	0	0	-0.0	1.0	0	0
555	-0.0	1.0	0.0	1.0	0	0	-0.0	1.0	0	0
556	0.0	1.0	0.0	1.0	0	0	0.01	0.997	0	0
557	36.73	0.0	9.83	0.007	2	0	10.27	0.006	2	0
558	0.02	1.0	0.0	0.999	0	0	-0.0	1.0	0	0
559	140.39	0.0	31.93	0.0	5	4	34.95	0.0	8	5
560	0.0	1.0	0.0	0.999	0	0	-0.0	1.0	0	0
'''

import os, sys, re, argparse
from collections import defaultdict
from scipy import stats
import numpy as np

rundir = os.getcwd()

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='./%(prog)s [options] -d directory -o output_basename',
    description = '''    Loops through a directory of paml output files and extracts 
    LRT test statistics and number of BEB significant sites for M2 and M8 model.''',
    
    epilog = """Written by Stephen A. Wyka (2020)""",
    formatter_class = MyFormatter)

parser.add_argument(
    '-d',
    '--directory',
    required = True,
    help = 'Input directory',
    metavar=''
)
parser.add_argument(
    '-o',
    '--output',
    required = True,
    help = 'Base name of output files',
    metavar=''
)
args=parser.parse_args()

model_list = ['Model 0', 'Model 1', 'Model 2', 'Model 3', 'Model 7', 'Model 8']

if __name__ == "__main__":
    input_dir = os.path.abspath(args.directory)
    codeML_dict = {} # {Cluster_ID : {Model 0 : [lnL], Model 1 : [lnL], Model 8 : [lnL, sites>95, sites>99]}
    files = sorted([f for f in os.listdir(input_dir) if '.paml_out' in f], key=lambda x: int(x.split('.')[0]))
    for f in files:
        cluster = f.split('.')[0]
        codeML_dict[cluster] = {m : [] for m in model_list}
        fpath = os.path.join(input_dir, f)
        model_breaks = defaultdict(list) # dictionary with Model's as keys and lists of lines associated with give model output
        with open(fpath, 'r') as f_in:
            model = ''
            for line in f_in:
                model_start = re.search(r'(Model) (\d)(:)', line)
                if model_start != None:
                    model = 'Model '+model_start.group(2)
                    continue
                model_breaks[model].append([line])
        model_breaks.pop('')

        for model, lines in model_breaks.items():
            codeML_dict[cluster][model] = []
            beb_list = []
            l_c = -1
            for l in lines:
                l_c += 1
                l = l[0]
                if 'lnL' in l:
                    l_data = re.search(r'(lnL)(.+np:) (\d+)(\):)\s(.+)\s([+-])',l)
                    lnL = float(l_data.group(5).strip())
                    codeML_dict[cluster][model].append(lnL)
                elif 'Bayes Empirical Bayes (BEB) analysis' in l:
                    num_sig_95 = 0
                    num_sig_99 = 0
                    for i in range(l_c, len(lines)):
                        sites = re.search(r'(\d+) (\w) (.+) (.+) (\+\-) (.+)',lines[i][0])
                        if sites:
                            if '**' in sites.group(3):
                                num_sig_99 += 1
                                num_sig_95 += 1
                            elif '*' in sites.group(3):
                                num_sig_95 += 1
                    codeML_dict[cluster][model].append(num_sig_95)
                    codeML_dict[cluster][model].append(num_sig_99)
                else:
                    pass

    for cluster, results in codeML_dict.items():
        for model in results:
            if not codeML_dict[cluster][model]:
                codeML_dict[cluster][model] = [float('nan')]

    overall_output = os.path.join(rundir, args.output+'_parsred_paml.txt')
    lrt_output = os.path.join(rundir, args.output+'_lrt_stats.txt')
    lrt_results = [['Cluster\tM3vM0_lrt\tM3vM0_p\tM2vM1_lrt\tM2vM1_p\tM2_num_sites_>95\tM2_num_sites_>99\tM8vM7_lrt\tM8vM7_p\tM8_num_sites_>95\tM8_num_sites_>99']]
    with open(overall_output, 'w') as over_out, open(lrt_output, 'w') as lrt_out:
        over_out.write('Cluster\tModel\tnp\tlnL\tNumber_sites_>95\tNumber_sites_>99\n')
        for cluster, results in codeML_dict.items():
            for model in results:
                res = [str(x) for x in codeML_dict[cluster][model]]
                over_out.write(cluster+'\t'+model+'\t'+'\t'.join(res)+'\n')

            tmp_list = []
            m3_m0 = [codeML_dict[cluster]['Model 3'][0], codeML_dict[cluster]['Model 0'][0]]
            m2_m1 = [codeML_dict[cluster]['Model 2'][0], codeML_dict[cluster]['Model 1'][0]]
            if not np.isnan(codeML_dict[cluster]['Model 2'][0]):
                m2_95 = codeML_dict[cluster]['Model 2'][1]
                m2_99 = codeML_dict[cluster]['Model 2'][2]
            else:
                m2_95 = 'nan'
                m2_99 = 'nan'
            m8_m7 = [codeML_dict[cluster]['Model 8'][0], codeML_dict[cluster]['Model 7'][0]]
            if not np.isnan(codeML_dict[cluster]['Model 8'][0]):
                m8_95 = codeML_dict[cluster]['Model 8'][1]
                m8_99 = codeML_dict[cluster]['Model 8'][2]
            else:
                m8_95 = 'nan'
                m8_99 = 'nan'
            lrt_3v0 = 2*(m3_m0[0] - m3_m0[1])
            p_3v0 = (1- stats.chi2.cdf(lrt_3v0, 4))
            lrt_2v1 = 2*(m2_m1[0] - m2_m1[1])
            p_2v1 = (1- stats.chi2.cdf(2*(m2_m1[0] - m2_m1[1]), 2))
            lrt_8v7 = 2*(m8_m7[0] - m8_m7[1])
            p_8v7 = (1- stats.chi2.cdf(2*(m8_m7[0] - m8_m7[1]), 2))
            tmp_list = [cluster, round(lrt_3v0, 2), round(p_3v0, 3), 
                round(lrt_2v1, 2), round(p_2v1, 3), m2_95, m2_99,
                round(lrt_8v7, 2), round(p_8v7, 3), m8_95, m8_99]
            lrt_results.append([str(x) for x in tmp_list])
        for item in lrt_results:
            lrt_out.write('\t'.join(item) + '\n')