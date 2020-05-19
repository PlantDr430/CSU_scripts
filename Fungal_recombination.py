#!/usr/bin/python3

'''
This script is the working pipeline used in Wyka et al. 2020 which was created based 
on the workflow from Stukenbrock and Dutheil 2018 Genetics, 208(3), 1209-1229. It takes 
a directory of alignment windows created from whole genome alingments and filtered using 
MafFilter and output in LDhat format. This also requires a data file from MafFilter to obtain 
correct nucleotide positions within the genome. These files are then used for running with LDhat 
to create a fine-scale recombination landscape around gene features and genes (specific types) if 
provided. It can also prepare the LDhat results for a subsequent run of LDhot to predict recombination 
hotspots.

LDhot is run in parallel, so it is best to used the non-multithreaded version of LDhot. This script will 
farm out multiple alignment blocks to single cores and actually results in much faster LDhot run (the more 
cores you use). 24 cores took around 2 days for 350x 100kb alignments.
'''


import os, sys, re, argparse, subprocess, warnings, math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from sklearn import datasets, linear_model
from sklearn.linear_model import LinearRegression
from scipy import stats
from collections import Counter, defaultdict, OrderedDict
from collections.abc import Iterable
from itertools import tee, count, islice, chain, zip_longest, combinations, repeat
from matplotlib.patches import Patch
from statsmodels.stats.multitest import multipletests
from matplotlib import ticker as mticker
from multiprocessing import Pool
import statsmodels.api as sm
from matplotlib.lines import Line2D

rundir = os.getcwd()
parentdir = os.path.dirname(rundir)
mpl.rc('font',family='Times New Roman')

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='python3 %(prog)s [options] -in input_directory -o output_basename -d data.csv',
    description = '''    ''',
    
    epilog = """Written by Stephen A. Wyka (2020)""",
    formatter_class = MyFormatter)

parser.add_argument(
    '-in',
    '--input',
    required=True,
    help = 'Directory containing alignment files in FASTA / LDhat format from MafFilter',
    metavar=''
)
parser.add_argument(
    '-o',
    '--out',
    required=True,
    help = 'Basename of output directory',
    metavar=''
)
parser.add_argument(
    '-d',
    '--data',
    required=True,
    help = 'Data csv output file from MafFilter corresponding to the input files. This '\
    'is used to rename alignment files by chromosome and get to correct for start '\
    'position in downstream analysis.',
    metavar=''
)
parser.add_argument(
    '-g',
    '--gff3',
    help = 'GFF3 file for extrating recombination estiamtes of certain regions (i.e. exons, itrons)',
    metavar=''
)
parser.add_argument(
    '-f',
    '--features',
    nargs= '+',
    default = ['all'],
    help = 'Features to extract from GFF3 file [all|Exon|Intron|Upstream|Downstream|Intergenic]. '\
    '[default: all]. For multiple, besides all, provide in a list separated by a space '\
    '(i.e. Exon Intron Upstream)',
    metavar=''
)
parser.add_argument(
    '-bp',
    '--basepairs',
    type=int,
    default=500,
    help = 'Amount of base pairs for searching upstream and downstream regions. [default: 500]',
    metavar=''
)
parser.add_argument(
    '-a',
    '--annotations',
    nargs='*',
    help = 'Files containing list of genes to narrow the search and report recombination only surrounding ,\
    these genes (i.e. effector genes, metabolite genes).',
    metavar=''
)
parser.add_argument(
    '-th',
    '--theta',
    type=float,
    default=0.005,
    help = 'Theta value for likelihood tables [default: 0.005]',
    metavar=''
)
parser.add_argument(
    '-its',
    '--iterations',
    type=int,
    default=10000000,
    help = 'Number of MCMC iterations [default: 10000000]',
    metavar=''
)
parser.add_argument(
    '-samp',
    '--sample',
    type=int,
    default=5000,
    help = 'Number of MCMC iterations between samples [default: 5000]',
    metavar=''
)
parser.add_argument(
    '--bpen',
    type=int,
    default=10,
    help = 'Background block penalty, between 0 - 50 [default: 10]',
    metavar=''
)
parser.add_argument(
    '--burnin',
    type=int,
    default=20,
    help = 'Specify the number of samples to be discarded as burn-in. [default: 20] It '\
    'should be noted that a burn-in of 20 with a sample iteration of 5000 is a burn-in '\
    'of 100000 of the rjMCMC scheme',
    metavar=''
)
parser.add_argument(
    '--missfreq',
    type=float,
    default=0.0,
    help = 'Frequency cut-off for missing data, default is to remove all missing data (i.e. gaps). '\
    'Note: if running LDhot this needs to be 0.0 other wise you risk LDhot encountering errors '\
    '[defuult: 0.0]',
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
parser.add_argument(
    '-multi',
    '--multitest',
    default='fdr_bh',
    choices=['bonferroni', 'sidak', 'holm-sidak','holm','simes-hochberg','hommel',\
        'fdr_bh','fdr_by','fdr_tsbh','fdr_tsbky'],
    help = 'Multi-test correction to use [bonferroni|sidak|holm-sidak|holm|'\
        'simes-hochberg|hommel|fdr_bh|fdr_by|fdr_tsbh|fdr_tsbky] [defaut: fdr_bh]',
    metavar=''
)
parser.add_argument(
    '--LDhot',
    action='store_true',
    help = 'Additionally run LDhot using the recombination maps estimated from LDhat [default: OFF].',
)
parser.add_argument(
    '-c',
    '--cpus',
    type=int,
    default=1,
    help = 'Number of cores to use for multiprocessing with LDhot [default: 1]',
    metavar=''
)
parser.add_argument(
    '-nsim',
    '--numbersims',
    type=int,
    default=1000,
    help = 'Number of simulations to run if running LDhot. [default: 1000]',
    metavar=''
)
args=parser.parse_args()

def which_path(file_name):
    for path in os.environ["PATH"].split(os.pathsep):
        full_path = os.path.join(path, file_name)
        if os.path.exists(full_path) and os.access(full_path, os.X_OK):
            return full_path
    return None

def get_next(some_iterable):
    items, nexts = tee(some_iterable, 2)
    nexts = chain(islice(nexts, 1, None), [None])
    return zip_longest(items, nexts)

def flatten(lis):
    for item in lis:
     if isinstance(item, Iterable) and not isinstance(item, str):
         for x in flatten(item):
             yield x
     else:
         yield item

if not os.path.isdir(args.out):
    os.makedirs(os.path.join(args.out, 'working_directory', 'alignments'))
    os.makedirs(os.path.join(args.out, 'LDhat_results'))
result_dir = os.path.abspath(os.path.join(rundir, args.out))
work_dir = os.path.abspath(os.path.join(result_dir, 'working_directory'))
align_dir = os.path.abspath(os.path.join(work_dir, 'alignments'))
chr_dir = os.path.abspath(os.path.join(result_dir, 'LDhat_results'))
# double check directories are present and made
dirs = [os.path.join(args.out, 'working_directory'),os.path.join(args.out, 'LDhat_results'), 
os.path.join(work_dir, 'alignments')]
for d in dirs:
    if not os.path.isdir(d):
        os.makedirs(d)

make sure LDhat and LDhot (if called) programs are callable from $PATH
try:
    if which_path('complete'):
        COMPLETE = 'complete'
    else:
        raise
except:
    print('complete from LDhat not found, please make sure parent directory of ',\
    'complete is located in $PATH')
    sys.exit()

try:
    if which_path('convert'):
        CONVERT = 'convert'
    else:
        raise
except:
    print('convert from LDhat not found, please make sure parent directory of ',\
    'convert is located in $PATH')
    sys.exit()

try:
    if which_path('interval'):
        INTERVAL = 'interval'
    else:
        raise
except:
    print('interval from LDhat not found, please make sure parent directory of ',\
    'interval is located in $PATH')
    sys.exit()

try:
    if which_path('stat'):
        STAT = 'stat'
    else:
        raise
except:
    print('stat from LDhat not found, please make sure parent directory of ',\
    'stat is located in $PATH')
    sys.exit()

if args.LDhot:
    try:
        if which_path('ldhot'):
            LDHOT = 'ldhot'
        else:
            raise
    except:
        print('ldhot from LDhot not found, please make sure parent directory ',\
        'of ldhot is located in $PATH')
        sys.exit()
    try:
        if which_path('ldhot_summary'):
            LDHOT_SUM = 'ldhot_summary'
        else:
            raise
    except:
        print('ldhot_summary from LDhot not found, please make sure parent directory ',\
        'of ldhot_summary is located in $PATH')
        sys.exit()

def uniquify(seq, suffs = count(1)):
    '''Make all the items unique by adding a suffix (1, 2, etc).

    `seq` is mutable sequence of strings.
    `suffs` is an optional alternative suffix iterable.
    '''
    not_unique = [k for k,v in Counter(seq).items() if v>0] # so we have: ['foo', 'bar']
    # suffix generator dict - e.g., {'foo': <my_gen>, 'bar': <my_gen>}
    suff_gens = dict(zip(not_unique, tee(suffs, len(not_unique))))  
    for idx,s in enumerate(seq):
        try:
            suffix = str(next(suff_gens[s]))
        except:
            # s was unique
            continue
        else:
            seq[idx] += suffix

def create_unqiue_chromosome_list(filename):
    '''Get chromosomes from data.csv and append unique suffixes 
    to generate list for renameing MafFilter alingment files.
    '''
    align_start_dict = {} # {Chromosome-subset : start position from alignment}
    align_length_dict = {} # {Chromosome-subset : length of alignment block}
    df = pd.read_csv(filename, delimiter='\t')
    chr_list = [str(x) for x in df['Chr'].tolist()]
    uniquify(chr_list, (f'-{x!s}' for x in range(1, 100)))
    for i in range(0,len(chr_list)):
        align_start_dict[chr_list[i]] = df['Start'].tolist()[i]
        align_length_dict[chr_list[i]] = df['BlockLength'].tolist()[i]
    return chr_list, align_start_dict, align_length_dict

def check_fasta_LDhat_format(filename):
    '''Create fasta dictionary to check sequence lengths with 
    LDhat header sequence lengths. Also, grab sample size number.
    '''
    fasta_dict = defaultdict(str)
    with open(filename, 'r') as in_fasta:
        next(in_fasta) # skip LDhat header row
        name = ''
        for line in in_fasta:
            if line.startswith('>'):
                name = line[1:-1]
                continue
            fasta_dict[name]+=line.strip()
    with open(filename, 'r') as in_fasta:
        ldhat_line = in_fasta.readlines()[0]
        sample_size, length, ploidy = ldhat_line.split()
    for isolate, sequence in fasta_dict.items():
        if len(sequence) == int(length):
            pass
        else:
            print('ERROR: Sequence lengths in LDhat formatted fasta files do ',\
            'not match sequence lengths of actual sequences')
            sys.exit()
    return sample_size

def create_parallel_dict(len_dict, n):
    '''Create dictionary of ~equal number of alignment blocks (based on alignment length) 
    for running LDhot with multiprocessing. This tries to create chunks of around similar 
    sized alignments so that LDhot should finish roughly around the same time for all chunks.
    '''
    parallel_dict = {k : [] for k in range(n)} # dictionary of chuncks of alignments based on cpu number
    order_list = list(np.arange(n)) + list(np.arange(n)[::-1]) # asecending and descending list of cpus used
    sort_len_dict = {k: v for k, v in sorted(len_dict.items(), key=lambda item: item[1], reverse=True)}
    chunk_size = math.ceil((len(sort_len_dict)/n)/2) # number of alignment blocks per chuink (rounded up)
    full_order_list = list(flatten(list(repeat(order_list, chunk_size))))[:len(sort_len_dict)] # repeat order by chunk_size

    c = -1 # keep count for indexing
    for k in sort_len_dict.keys():
        c += 1
        chunk = full_order_list[c]
        parallel_dict[chunk].append(k)

    return parallel_dict

def check_inputs():
    '''Check for correct input and formats and create some dictionaries and lists 
    for future purposes.
    '''
    if os.path.isdir(args.input):
        input_dir = os.path.abspath(args.input)
    else:
        parser.error('ERROR: Input is not a directory, please make sure to pass a directory ',\
        'to -in or --input')
    if args.data:
        data_csv = os.path.abspath(args.data)
    else:
        parser.error('ERROR: Please provide data.csv file from MafFilter to -d or --data')

    chr_list, align_start_dict, align_length_dict = create_unqiue_chromosome_list(data_csv)
    # create sorted list of alignment fasta files
    align_list = sorted([f for f in os.listdir(input_dir)], key=lambda k: int(k.split('t')[1].split('.')[0]))
    if not len(os.listdir(align_dir)) == len(os.listdir(input_dir)):
        # loop through lists and create symlinks with alignment names corresponding to chromosomes
        for i in range(0, len(chr_list)):
            subprocess.call(['ln', '-s', os.path.abspath(os.path.join(input_dir, align_list[i])), 
            os.path.join(align_dir, 'alingment_'+chr_list[i]+'.fasta')], cwd=rundir)

    elif args.LDhot: # create dictionary for splitting folders for multiprocessing
        parallel_dict = create_parallel_dict(align_length_dict, args.cpus)

    else:
        pass

    sample_size = check_fasta_LDhat_format(os.path.join(align_dir, os.listdir(align_dir)[0]))
    chr_start_dict = {} # {Chromosome number (no subset) : starting positon from alingment}
    for k, v in align_start_dict.items():
        if '-1' in k: # subset one is start of the contig/chromosome
            chr_start_dict[k.split('-')[0]] = v
    return align_start_dict, sample_size, chr_start_dict, parallel_dict

def create_likelihood_table(n):
    '''Run complete function of LDhat to create likelihood table(s)'''
    print('Creating likelihood tables, this can take some time depending on the sample size')
    prefix = 'theta_'+str(args.theta)+'_lk_'
    if not os.path.exists(os.path.join(work_dir, prefix+'new_lk.txt')):
        try:
            with open(os.path.join(work_dir, 'lk_tables.log'), 'w') as logfile:
                subprocess.call([COMPLETE, '-n', n, '-rhomax', '100', '-n_pts', '101', '-theta', 
                str(args.theta), '-prefix', prefix], cwd=work_dir, stdout = logfile, stderr = logfile)
        except:
            print(('There was an error in likelihood table creation, check the logfile located at {}')
            .format(os.path.join(work_dir, 'lk_tables.log')))
    else: # files exists so we don't need to re-run
        pass
    return prefix

def convert_fasta_files():
    '''Run convert function of LDhat'''
    print('Converting FASTA files')
    for files in os.listdir(align_dir):
        file_path = os.path.join(align_dir, files)
        dir_path = os.path.join(chr_dir, files.split('.')[0])
        site_file = os.path.join(dir_path, 'sites.txt')
        loc_file = os.path.join(dir_path, 'locs.txt')
        freq_file = os.path.join(dir_path, 'freqs.txt')
        if not os.path.isdir(dir_path):
            os.makedirs(dir_path)
        if not os.path.exists(site_file) and not os.path.exists(loc_file) and not os.path.exists(freq_file):
            try:
                with open(os.path.join(dir_path, 'convert.log'), 'w') as logfile:
                    subprocess.call([CONVERT, '-seq', file_path, '-2only', 
                    '-missfreqcut', str(args.missfreq)
                    ], cwd=dir_path, stdout=logfile, stderr=logfile)
            except:
                print(('There was an error in converting the fasta files, check the logfile located at {}')
                .format(os.path.join(dir_path, 'convert.log')))
                sys.exit()
        else: # files exists so we don't need to re-run
            pass

def run_interval():
    '''Run interval function of LDhat'''
    print('Running interval program')
    for dirs in os.listdir(chr_dir):
        dir_path = os.path.join(chr_dir, dirs)
        rate_out = os.path.join(dir_path, lk_prefix+'rates.txt')
        bound_out = os.path.join(dir_path, lk_prefix+'bounds.txt')
        table_out = os.path.join(dir_path, lk_prefix+'type_table.txt')
        if not os.path.exists(rate_out) and not os.path.exists(bound_out) and not os.path.exists(table_out):
            loc_file = os.path.join(dir_path, 'locs.txt')
            site_file = os.path.join(dir_path, 'sites.txt')
            lk_file = os.path.join(work_dir, lk_prefix+'new_lk.txt')
            try:
                with open(os.path.join(dir_path, 'interval.log'), 'w') as logfile:
                    subprocess.call([INTERVAL, '-seq', site_file, '-loc', loc_file, '-lk', lk_file, 
                    '-its', str(args.iterations), '-samp', str(args.sample), '-bpen', str(args.bpen), 
                    '-seed', '1234', '-prefix', lk_prefix], cwd=dir_path, stdout=logfile, stderr=logfile)
            except:
                print(('There was an error in running interval, check the logfile located at {}')
                .format(os.path.join(dir_path, 'interval.log')))
                sys.exit()
        else: # files exists so we don't need to re-run
            pass

def summarize_LDhat_results():
    '''Run stat function of LDhat'''
    print('Summarizing results from interval')
    for dirs in os.listdir(chr_dir):
        dir_path = os.path.join(chr_dir, dirs)
        summarized_out = os.path.join(dir_path, lk_prefix+'summarized_rates_res.txt')
        if not os.path.exists(summarized_out):
            rate_results = os.path.join(dir_path, lk_prefix+'rates.txt')
            loc_file = os.path.join(dir_path, 'locs.txt')
            try:
                with open(os.path.join(dir_path, 'stat.log'), 'w') as logfile:
                    subprocess.call([STAT, '-input', rate_results, '-burn', str(args.burnin), 
                    '-loc', loc_file, '-prefix', lk_prefix+'summarized_rates_'
                    ], cwd=dir_path, stdout=logfile, stderr=logfile)
                if not os.path.exists(summarized_out): # additional check as sometimes errors slip through
                    raise
            except:
                print(('There was an error in running stat, check the logfile located at {}')
                .format(os.path.join(dir_path, 'stat.log')))
                sys.exit()
        else: # files exists so we don't need to re-run
            pass

def merge_summarized_LDhat_results():
    '''Merge summarized outputs based on chromosome. Also, add in correct loci column 
    that will allow proper matching to GFF file for gene region extraction. 
    This also filters the summarized files to remove estimates where the width 
    of the confidence interval is >= 2x the mean (Stukenbrock and Dutheil 2018).
    
    C_Loci = Correct Loci corresponding to what positions should be in reference to 
    the GFF3 positions. O_Loci = Original Loci corresponding to results from LDhat. During 
    processing of whole genome alignments I noticed that even after projecting my alignments 
    to a reference genome and merging alignments by chromosome/contig sometimes the starts of 
    the alignments were not "1" but could be 1,515 (for example). While the position 1,515 is 
    saved to the data.csv outputs from MafFilter, when these files are processed by LDhat the 1,515 
    correct start position would be assumed to be 1 in LDhat and the SNPs would be called sequentially 
    so that 1,515 is position 1, 1,516 is position 1, etc. This meant that if we want to extract gene 
    regions from a corresponding GFF3 we needed to correct the Loci positions in the LDhat outputs. This 
    is done by using the data.csv output file from MafFilter corresponding to the alignments used in LDhat 
    that have the correct positions of the start sites.
    
    For example, if alignment_1-1 for conting 1 subset 1 had a MafFilter start position of 1,515 and LDhat 
    start site of 1. Using the dictionary created earlier we have the {Chr-subset : start site} so we would 
    add 1,515 to all loci of the alignment_1-1 LDhat output. For chr's or contigs with multiple subsets. The 
    correct start sites of those would be the LDhat output loci (1) + correct start site of the chr/contig (1,515) 
    + the correct start site of the subset (i.e. 100,015 for example).'''
    print('Merging summarized LDhat results')
    filter_merged_out = os.path.join(work_dir, 'Filtered_merged_summarized_LDhat_results.tsv')
    unfilter_merged_out = os.path.join(work_dir, 'Unfiltered_merged_summarized_LDhat_results.tsv')
    if os.path.exists(filter_merged_out):
        merged_df = pd.read_csv(filter_merged_out, delimiter='\t')
    if os.path.exists(unfilter_merged_out):
        unfil_merged_df = pd.read_csv(unfilter_merged_out, delimiter='\t')
    else:
        merge_dict = OrderedDict()
        unfiltered_merge_dict = OrderedDict()
        chr_name_list = []
        dir_prefix = None
        for dirs in os.listdir(chr_dir): # get list of chromosomes from files
            dir_path = os.path.join(chr_dir, dirs)
            chr_num = dirs.split('_')[1]
            dir_prefix = dirs.split('_')[0]
            chr_name_list.append(chr_num)
        # sort by chromosome number and subsets
        sorted_chr_name_list = sorted(chr_name_list, key=lambda x: (int(x.split('-')[0]),int(x.split('-')[1])))

        for chr in sorted_chr_name_list:
            sum_results_path = os.path.join(chr_dir, dir_prefix+'_'+chr, lk_prefix+'summarized_rates_res.txt')
            chr_num = chr.split('-')[0]
            chr_subset = chr.split('-')[1]
            tmp_list = []
            unfiltered_tmp = []
            with open(sum_results_path, 'r') as sum_r:
                for i in range(2):
                    next(sum_r)
                for line, next_line in get_next(sum_r):
                    loci, mean, median, L95, U95 = line.split()
                    loci = str(int(float(loci)))
                    if next_line == None: # reached the end of the file
                        pass
                    else: # filter estimates by removing estiamtes where the width of the CI is >= 2x the mean
                        next_loci = str(int(float(next_line.split()[0])))
                        if (float(U95) - float(L95)) >= float(mean)*2:
                            if int(chr_subset) == 1: # add column of correct start positions for each loci from MafFilter output
                                cor_pos = str(int(float(loci)) + int(chr_start_dict[chr_num]))
                                cor_next = str(int(float(next_loci)) + int(chr_start_dict[chr_num]))
                            elif int(chr_subset) > 1:
                                cor_pos = str(int(float(loci))+int(chr_start_dict[chr_num])+int(align_start_dict[chr]))
                                cor_next = str(int(float(next_loci))+int(chr_start_dict[chr_num])+int(align_start_dict[chr]))
                            unfiltered_tmp.append([cor_pos, cor_next,loci,next_loci, mean, L95, U95, str(int(cor_next)-int(cor_pos))])

                        else: # same as above, but creating full (unfiltered version)
                            if int(chr_subset) == 1:
                                cor_pos = str(int(float(loci)) + int(chr_start_dict[chr_num]))
                                cor_next = str(int(float(next_loci)) + int(chr_start_dict[chr_num]))
                            elif int(chr_subset) > 1:
                                cor_pos = str(int(float(loci))+int(chr_start_dict[chr_num])+int(align_start_dict[chr]))
                                cor_next = str(int(float(next_loci))+int(chr_start_dict[chr_num])+int(align_start_dict[chr]))
                            tmp_list.append([cor_pos, cor_next,loci,next_loci, mean, L95, U95, str(int(cor_next)-int(cor_pos))])
                            unfiltered_tmp.append([cor_pos, cor_next,loci,next_loci, mean, L95, U95, str(int(cor_next)-int(cor_pos))])
                if int(chr_subset) == 1:
                    merge_dict[chr_num] = tmp_list
                    unfiltered_merge_dict[chr_num] = unfiltered_tmp
                elif int(chr_subset) > 1:
                    merge_dict[chr_num] = merge_dict[chr_num]+ tmp_list
                    unfiltered_merge_dict[chr_num] = unfiltered_merge_dict[chr_num] + unfiltered_tmp

        # write out the filtered merged results
        with open(filter_merged_out, 'w') as output:
            output.write('Chr\tC_Loci_1\tC_Loci_2\tO_Loci_1\tO_Loci_2\tMean_rho\tL95\tU95\tWeights\n')
            for k, v in merge_dict.items():
                for x in v:
                    output.write(k + '\t' + '\t'.join(x) + '\n')
        merged_df = pd.read_csv(filter_merged_out, delimiter='\t')

        with open(unfilter_merged_out, 'w') as output:
            output.write('Chr\tC_Loci_1\tC_Loci_2\tO_Loci_1\tO_Loci_2\tMean_rho\tL95\tU95\tWeights\n')
            for k, v in unfiltered_merge_dict.items():
                for x in v:
                    output.write(k + '\t' + '\t'.join(x) + '\n')
        unfil_merged_df = pd.read_csv(unfilter_merged_out, delimiter='\t')

    return merged_df, unfil_merged_df

def create_full_genomic_range_LDhat(df1, df2):
    print('Getting full genomic range of recombination rates. This will take a decent amount of time.')
    full_bp_out = os.path.join(work_dir, 'All_basepair_filtered_LDhat_results.tsv')
    unfull_bp_out = os.path.join(work_dir, 'All_basepair_unfiltered_LDhat_results.tsv')
    if os.path.exists(full_bp_out):
        full_bp_df = pd.read_csv(full_bp_out, delimiter='\t')
    else:
        rows_list = [['Chr','Loci','Mean_rho','Weights']]
        for index, row in df1.iterrows():
            chr = row['Chr']
            c_loci_1 = int(row['C_Loci_1'])
            c_loci_2 = int(row['C_Loci_2'])
            means = row['Mean_rho']
            weights = row['Weights']
            loci_range = list(range(c_loci_1, c_loci_2))
            for bp in loci_range:
                tmp_list = []
                tmp_list = [chr,bp,means,weights]
                rows_list.append(tmp_list)
        full_bp_df = pd.DataFrame(rows_list[1:],columns=rows_list[0])
        full_bp_df.to_csv(full_bp_out, sep='\t', index=False)

    if os.path.exists(unfull_bp_out):
        unfull_bp_df = pd.read_csv(unfull_bp_out, delimiter='\t')
    else:
        rows_list = [['Chr','Loci','Mean_rho','Weights']]
        for index, row in df2.iterrows():
            chr = row['Chr']
            c_loci_1 = int(row['C_Loci_1'])
            c_loci_2 = int(row['C_Loci_2'])
            means = row['Mean_rho']
            weights = row['Weights']
            loci_range = list(range(c_loci_1, c_loci_2))
            for bp in loci_range:
                tmp_list = []
                tmp_list = [chr,bp,means,weights]
                rows_list.append(tmp_list)
        unfull_bp_df = pd.DataFrame(rows_list[1:],columns=rows_list[0])
        unfull_bp_df.to_csv(full_bp_out, sep='\t', index=False)

    return full_bp_df, unfull_bp_df

def create_gff_dictionary(gff_file):
    '''Turn GFF3 file into a dictionary for downstream analysis.'''
    g_dict = {} # {Chr : {(start,stop) : (Gene, feature)}}
    gene_dict = {} # just a list of genes as dict keys for faster searching
    with open(gff_file, 'r') as gff_in:
        for line in gff_in:
            if not line.startswith('#'):
                col = line.strip().split('\t')
                chr = int(col[0])
                g_start = int(col[3])
                g_stop = int(col[4])
                g_strand = col[6]
                if chr not in g_dict.keys():
                        g_dict[chr] = {}
                else:
                    g_dict[chr] = g_dict[chr]
                if col[2] == 'gene':
                    name = re.search(r'(ID=)([^;]*)', col[8]).group(2)+'-T1' # specific version
                    # name = re.search(r'(ID=)([^;]*)', col[8]).group(2) # robust version
                    if name not in gene_dict.keys():
                        gene_dict[name] = []
                    prime_5 = ((g_start - args.basepairs), g_start) # start of gene minus length of requested basepairs
                    if prime_5[0] < 0: # if the first coordinate is negative, set to 0
                        prime_5 = (0, g_start)
                    prime_3 = (g_stop, (g_stop + args.basepairs)) # stop of gene plus length of requested basepairs
                    if g_strand == '+':
                        g_dict[chr][prime_5[0],prime_5[1]] = (name,'Upstream')
                        g_dict[chr][prime_3[0],prime_3[1]] = (name,'Downstream')
                    elif g_strand == '-': # technically upstream and downstream are switched
                        g_dict[chr][prime_5[0],prime_5[1]] = (name,'Downstream')
                        g_dict[chr][prime_3[0],prime_3[1]] = (name,'Upstream')
                elif col[2] == 'exon':
                    name = re.search(r'(Parent=)([^;]*)', col[8]).group(2)
                    g_dict[chr][g_start,g_stop] = (name, 'Exon')
                elif col[2] == 'intron':
                    name = re.search(r'(Parent=)([^;]*)', col[8]).group(2)
                    g_dict[chr][g_start,g_stop] = (name, 'Intron')
    for chr, coords in g_dict.items():
        c_list = [c for c in coords.keys() if coords[c][1] == 'Downstream' or coords[c][1] == 'Upstream']
        for i in range(0, len(c_list)):
            if i == 0 and c_list[i][0] != 0:
                first_feat_start = c_list[i][0]
                # first intergenic region is 0 to start of down/upstream of first gene
                intergenic = (0, first_feat_start -1) 
                g_dict[chr][intergenic] = ('inter_{}_{}'.format('chr:'+str(chr),i), 'Intergenic')
            elif i == len(c_list)-1:
                last_feat_stop = c_list[i][1]
                # last intergenic region is stop of down/upstream of last gene to end of contig
                # currently, no easy way to get contig end, setting to 1,000,000 after last gene should be sufficient
                # the true end doesn't matter as we only are using this as a search for SNP pairs found in certain regions.
                # i.e. if a SNP pairs occurs in the region from the last gene to the true end of the contig, it will also 
                # occur in the region from the last gene + 1,000,000 bps (unless the true end if > 1,000,000 bp away from the last gene.
                # Also, there won't be any SNPs called in regions that don't exist in the contig.
                intergenic = (last_feat_stop+1, last_feat_stop+1000000)
                g_dict[chr][intergenic] = ('inter_{}_{}'.format('chr:'+str(chr),i), 'Intergenic')
            else:
                if coords[c_list[i]][0] != coords[c_list[i+1]][0]: # filter for between genes (i.e. if gene1_name != gene2_name)
                    gene1_stop = c_list[i][1]
                    gene2_start = c_list[i+1][0]
                    diff = gene2_start - gene1_stop
                    # sometimes upstreams and downstreams overlap. If overlap (neg.) or 0 diff. between coords, no intergenic region
                    if not diff <= 0:
                        intergenic = (gene1_stop+1, gene2_start-1)
                        g_dict[chr][intergenic] = ('inter_{}_{}'.format('chr:'+str(chr),i), 'Intergenic')

    return g_dict, gene_dict

def weighted_avg(group):
    d = group['Mean_rho']
    w = group['Weights']
    w_avg = np.average(d, weights=w)
    return w_avg

def extraction_multiprocessing(df):
    print('Extracting in multiprocess')
    gff = gff_dict
    processing_list = []
    for index, row in df.iterrows():
        chr = row['Chr']
        loci = row['Loci']
        means = row['Mean_rho']
        weights = row['Weights']
        if chr in gff.keys():
            # search to see if SNP pair is located within the bounds of a gene feature (i.e. exon, intron)
            key = [k for k in gff[chr].keys() if k[0] <= loci <= k[1]]
            if key:
                if len(key) > 1: # if there is overlap (i.e. upstream regions overlapping with downstream regions)
                    # we keep both and classify SNP pair for both regions, since it technically is part of both
                    # However, if there is an overlap between exonic/intronic regiond and down/upstream regions
                    # we will discard the down/upstream classification and classify the SNP as part of exonic/intronic
                    genes = [gff[chr][k][0] for k in key]
                    feats = [gff[chr][k][1] for k in key]
                    if not genes[0] == genes[1]:
                        if (feats[0] == 'Exon' or feats[0] == 'Intron') and (feats[1] == 'Downstream' or feats[1] == 'Upstream'):
                            tmp_list = [chr,genes[0],feats[0],loci,means,weights]
                            processing_list.append(tmp_list)
                        elif (feats[0] == 'Downstream' or feats[0] == 'Upstream') and (feats[1] == 'Exon' or feats[1] == 'Intron'):
                            tmp_list = [chr,genes[1],feats[1],loci,means,weights]
                            processing_list.append(tmp_list)
                        else:
                            pass
                    else: # if genes are the same this is just simple overlap of the features regions by 1 bp
                        for i in range(0, len(key)):
                            k = key[i]
                            tmp_list = [chr,gff[chr][k][0],gff[chr][k][1],loci,means,weights]
                            processing_list.append(tmp_list)
                else:
                    k = key[0]
                    tmp_list = [chr,gff[chr][k][0],gff[chr][k][1],loci,means,weights]
                    processing_list.append(tmp_list)
            else: # no gene regions found for SNP pair so return N/A's
                tmp_list = [chr,'N/A','N/A',loci,means,weights]
                processing_list.append(tmp_list)

    return processing_list

def extract_regions_LDhat(df, gff):
    '''Loop through GFF dictionary and calculate weighted mean_rhos for features from full_bp_df 
    and create new dictionary for downstream data anlysis.
    '''
    print('Extracting gene regions. This will take a while')
    classified_region_out = os.path.join(work_dir, 'Classified_regions_LDhat_results.tsv')
    summarized_region_out = os.path.join(result_dir, 'Summarized_regions_LDhat_results.tsv')
    if os.path.exists(summarized_region_out) and os.path.exists(classified_region_out):
        summary_df = pd.read_csv(summarized_region_out, delimiter='\t')
        filtered_class_df = pd.read_csv(classified_region_out, delimiter='\t')
    if os.path.exists(classified_region_out):
        filtered_class_df = pd.read_csv(classified_region_out, delimiter='\t')
        summary_df = filtered_class_df.groupby(['Chr','Gene','Feature']).apply(weighted_avg).reset_index(name='Weighted_rho_means')
        summary_df.to_csv(summarized_region_out, sep='\t', index=False)
    else:
        rows_list = [['Chr','Gene','Feature','Loci','Mean_rho','Weights']]
        df_list = np.array_split(df, args.cpus)
        pool = Pool(processes=args.cpus)
        results = pool.map(extraction_multiprocessing, df_list)
        pool.close()
        pool.join()
        processed_list = [[rows_list.append(y) for y in x] for x in list(results)]

        class_df = pd.DataFrame(rows_list[1:],columns=rows_list[0])
        count_df = class_df.drop(columns=['Loci','Mean_rho']).groupby(['Chr','Gene','Feature']).agg(['nunique']).reset_index(col_level=0)
        count_df.columns = count_df.columns.droplevel(1) # drop nunique column name
        count_dict = count_df.set_index('Chr').to_dict()['Weights'] # create dictionary of chr and unique SNP counts
        keep_list = [k for k in count_dict.keys() if count_dict[k] >= 3] # keep only regions that confidently have 2 different SNP pairs
        filtered_class_df = class_df[class_df['Chr'].isin(keep_list)]
        filtered_class_df.to_csv('NEW_Classified_regions_LDhat_results.tsv', sep='\t', index=False)

        print('Summarizing extracted gene regions')
        summary_df = filtered_class_df.groupby(['Chr','Gene','Feature']).apply(weighted_avg).reset_index(name='Weighted_rho_means')
        summary_df.to_csv(summarized_region_out, sep='\t', index=False)

    return filtered_class_df, summary_df

def fill_annotation_dict(gene_list, anno_files):
    '''Populates anno_dict with annotation type as keys and list of associated 
    genes. Returns a filtered_gene_list which is a list of genes not associated 
    with any of the given annotation genes (i.e. all other genes in the genome).
    '''
    filtered_gene_list = gene_list
    for file in anno_files:
        apath = os.path.abspath(os.path.join(rundir, file))
        anno_name = os.path.basename(apath).split('.')[0]
        with open(apath, 'r') as in_file:
            anno_list = [gene.strip() for gene in in_file if gene != '\n']
            filtered_gene_list = list(set(filtered_gene_list) - set(anno_list))
            anno_dict[anno_name] = anno_list
    return filtered_gene_list

def fill_anno_feat_dict(df, label):
    '''Loops through summarized dataframe and populates anno_feat_dict 
    for figure creation. We get anno_type by searching the rev_anno_dict 
    with the given gene.
    '''
    for index, row in df.iterrows():
        gene = row['Gene']
        feat = row['Feature']
        w_rho = row['Weighted_rho_means']
        if not str(feat) == 'nan':
            if label == 'anno_feat':
                if not 'inter_chr' in gene:
                    anno_type = rev_anno_dict[gene]
                    anno_feat_dict[anno_type][feat].append(w_rho)
            elif label == 'feat':
                feat_dict[feat].append(w_rho)
        else:
            pass

def set_box_colors(bp):
    plt.setp(bp['boxes'], linewidth=0.75, facecolor='white', alpha=0.5)
    plt.setp(bp['whiskers'], color='black', linewidth=0.75)
    plt.setp(bp['medians'], color='black', linewidth=0.75)
    plt.setp(bp['caps'], color='black', linewidth=0.75)

def compute_stats(data, labels, alt_labels, type):
    label_combos = list(combinations(labels,2))
    stats_output = os.path.abspath(os.path.join(result_dir, type+'_statistics.txt'))

    if type == 'Features':
        data_combos = list(combinations(data,2))
        uncor_p_list = []
        for i in range(0, len(data_combos)):
            w, p = stats.mannwhitneyu(data_combos[i][0], data_combos[i][1], alternative = 'two-sided')
            uncor_p_list.append(p)
        r, c_p, sf, bf = multipletests(uncor_p_list, alpha=args.multialpha, method=args.multitest) # fdr multitest
        corrected_p = [str(x) for x in c_p] # turn to strings to printing out
        if len(data) >= 3: # try 3 samples and 4 samples. If error tell user to alter this section
            try:
                s, kw_p = stats.kruskal(data[0],data[1],data[2],data[3],data[4], nan_policy='omit')
            except:
                try:
                    s, kw_p = stats.kruskal(data[0],data[1],data[2],data[3], nan_policy='omit')
                except:
                    try:
                        s, kw_p = stats.kruskal(data[0],data[1],data[2], nan_policy='omit')
                    except:
                        print('We took a guess and tried to run Krustal-Wallis with 3 to 5 samples. '\
                        'It appears we were wrong and you have more than 5 sample, to fix this alter lines 746 in '\
                        'the script to support the number of samples you are trying to run statistics on')
        with open(stats_output, 'w') as stats_out:
            for i in range(0, len(alt_labels)):
                stats_out.write(alt_labels[i]+':\n')
                if kw_p:
                    stats_out.write('Krustal-Wallis\t' + str(kw_p) + '\n')
                stats_out.write('dataset_1\tdataset_2\tP-value\n')
                for j in range(0, len(data_combos)):
                    stats_out.write('\t'.join(label_combos[j]) + '\t' + corrected_p[j] + '\n')

    elif type == 'Annotations_by_features' or type == 'Features_by_annotations':
        uncor_p_nest = []
        krustal_wallis_list = []
        for i in range(0, len(data)):
            dataset = []
            data_combos = list(combinations(data[i],2))
            for j in range(0, len(data_combos)):
                w, p = stats.mannwhitneyu(data_combos[j][0], data_combos[j][1], alternative = 'two-sided')
                dataset.append(p)
            uncor_p_nest.append(dataset)
            if len(data[i]) >= 3: # try 3 samples and 4 samples. If error tell user to alter this section
                try:
                    s, kw_p = stats.kruskal(data[i][0],data[i][1],data[i][2],data[i][3],data[i][4], nan_policy='omit')
                    krustal_wallis_list.append(kw_p)
                except:
                    try:
                        s, kw_p = stats.kruskal(data[i][0],data[i][1],data[i][2],data[i][3], nan_policy='omit')
                        krustal_wallis_list.append(kw_p)
                    except:
                        try:
                            s, kw_p = stats.kruskal(data[i][0],data[i][1],data[i][2], nan_policy='omit')
                            krustal_wallis_list.append(kw_p)
                        except:
                            print('We took a guess and tried to run Krustal-Wallis with 3 to 5 samples. '\
                            'It appears we were wrong and you have more than 5 sample, to fix this alter lines 746 in '\
                            'the script to support the number of samples you are trying to run statistics on')
        data_len = len(uncor_p_nest[0]) # get length of dataset for later
        uncor_p_list = list(flatten(uncor_p_nest)) # flatten list for multicorrection
        r, c_p_m, sf, bf = multipletests(uncor_p_list, alpha=args.multialpha, method=args.multitest) # fdr multitest mannwhitney
        corrected_p = [str(x) for x in c_p_m] # turn to strings to printing out
        corrected_p_nest = [corrected_p[i:i+data_len] for i in range(0, len(corrected_p), data_len)] # re-nest list

        with open(stats_output, 'w') as stats_out:
            for i in range(0, len(alt_labels)):
                if i == 0:
                    stats_out.write(alt_labels[i]+':\n')
                else: 
                    stats_out.write('\n'+alt_labels[i]+':\n')
                if krustal_wallis_list:
                    stats_out.write('Krustal-Wallis\t' + str(krustal_wallis_list[i]) + '\n')
                stats_out.write('dataset_1\tdataset_2\tP-value\n')
                for j in range(0, data_len):
                    stats_out.write('\t'.join(label_combos[j]) + '\t' + corrected_p_nest[i][j] + '\n')

def create_boxplots(d):
    print('Creating plots and running statistics')
    warnings.filterwarnings("ignore")
    if any(isinstance(i,dict) for i in d.values()) == False: # if not a nested dictionary
        figure_output = os.path.join(result_dir, 'Feature_plot.pdf')
        colors = plt.cm.tab10([i for i in range(len(d.keys()))])
        fig, ax = plt.subplots(1, 1, figsize=(6.69,2), dpi=1600)
        x_labels = [x for x in d.keys()]
        data = [[np.log10(x) for x in v] for v in d.values()] # log transform for better KDE
        non_log_data = [v for v in d.values()] # use for statistics
        # total_data = sum([len(v) for v in d.values()])
        # widths = [len(v)/total_data for v in d.values()]
        num_samples = [len(v) for v in data]
        vio = sns.violinplot(data=data, cut=0, palette=colors, linewidth=0.75, inner=None)
        positions = list(range(0,len(x_labels))) # adjust x-ticks for matplot boxplot
        bp = plt.boxplot(data, 0, '', widths=0.25, patch_artist=True, positions=positions, 
            showcaps=False)
        set_box_colors(bp)
        ax.set_axisbelow(True)
        ax.set_xticklabels(x_labels, fontname='Times New Roman', fontsize=12)
        plt.xlabel('Regions', fontname='Times New Roman', fontsize=12)
        plt.minorticks_on()
        plt.grid(which='minor', axis='y', color='white', linestyle='--', alpha=0.3, linewidth='0.75')
        ax.yaxis.grid(True, linestyle='-', linewidth='0.75', which='major', color='white')
        ax.yaxis.grid(True, linestyle='--', which='minor', color='white',  alpha=0.3, linewidth='0.75')
        ax.tick_params(axis='x', which='minor', bottom=False)
        ax.set_facecolor('gainsboro')
        plt.ylabel(r'Population Recombination Rate ($\rho$)', fontname='Times New Roman', fontsize=12)
        ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
        ax.yaxis.set_ticks([np.log10(x) for p in range(-5,1) for x in np.linspace(10**p, 10**(p+1), 10)], minor=True)
        for j in range(0,len(x_labels)): # add sample size data
            ax.text(j, -5.25, num_samples[j], horizontalalignment='center', fontsize=6, fontname='Times New Roman')
        # legend_ele = []
        # for i in range(0, len(x_labels)):
            # legend_ele.append(Patch(facecolor=colors[i],label=x_labels[i], edgecolor='black'))
        # ax.legend(handles=legend_ele, framealpha=1.0)

        compute_stats(non_log_data, x_labels, ['Genome'], 'Features')
        plt.yticks(fontname = 'Times New Roman', fontsize=12)
        plt.savefig(figure_output, bbox_inches = 'tight')
        plt.close()
        # plt.tight_layout()
        # plt.show()

    elif any(isinstance(i,dict) for i in d.values()) == True: # if a nested dictionary
        for k, v in d.items():
            d[k].pop('Intergenic') # remove intergenic keys as we can't match a particualt intergenic region to each gene (since they are between genes)
        alt_list = ['abf', 'fba'] # list for loop ['annotations by feature', 'features by annotations']
        for alt in alt_list:
            log_data = [[[np.log10(x) for x in v] for v in d[x].values()] for x in d.keys()] # log transform for better KDE
            non_log_data = [[v for v in d[x].values()] for x in d.keys()] # use for statistics
            if alt == 'abf':
                num_plots = len(d.keys())
                fig, ax = plt.subplots(1, num_plots, figsize=(6.69,2), dpi=1600, sharey=True, gridspec_kw={'wspace': 0.1})
                figure_output = os.path.join(result_dir, 'Annotation_by_feature_plot.pdf')
                x_ticks = [x for x in d[list(d.keys())[0]].keys()]
                x_labels = [x for x in d.keys()]
                colors = plt.cm.tab10([i for i in range(len(x_ticks))])
                data = log_data
                orig_data = non_log_data
                compute_stats(orig_data, x_ticks, x_labels, 'Annotations_by_features')
                legend_ele = []
                for i in range(0, len(x_ticks)):
                    legend_ele.append(Patch(facecolor=colors[i],label=x_ticks[i], edgecolor='black'))
            if alt == 'fba':
                num_plots = len(d[list(d.keys())[0]].keys())
                fig, ax = plt.subplots(1, num_plots, figsize=(6.69,2), dpi=1600, sharey=True, gridspec_kw={'wspace': 0.1})
                figure_output = os.path.join(result_dir, 'Features_by_annotation_plot.pdf')
                x_labels = [x for x in d[list(d.keys())[0]].keys()]
                x_ticks = [x for x in d.keys()]
                colors = plt.cm.tab10([i for i in range(len(x_ticks))])
                data = np.column_stack(log_data) # stack to group by feature instead of annotations
                orig_data = np.column_stack(non_log_data)
                compute_stats(orig_data, x_ticks, x_labels, 'Features_by_annotations')
                legend_ele = []
                for i in range(0, len(x_ticks)):
                    legend_ele.append(Patch(facecolor=colors[i],label=x_ticks[i], edgecolor='black', linewidth=0.75))
            plt.minorticks_on()
            fig.text(0.04, .5, r'Population Recombination Rate ($\rho$)', ha='center', va='center', rotation='vertical', 
                fontname='Times New Roman', fontsize=10)
            count = 0
            for i in range(0, num_plots):
                count += 1
                # total_data = sum([len(v) for v in data[i]])
                # widths = [len(v)/total_data for v in data[i]]
                num_samples = [len(v) for v in data[i]]
                ax[i].yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
                ax[i].yaxis.set_ticks([np.log10(x) for p in range(-5,1) for x in np.linspace(10**p, 10**(p+1), 10)], minor=True)
                vio = sns.violinplot(data=data[i], cut=0, palette=colors, linewidth=0.75, inner=None, ax=ax[i])
                positions = list(range(0,len(x_ticks))) # adjust x-ticks for matplot boxplot
                bp = ax[i].boxplot(data[i], 0, '', widths=0.25, patch_artist=True, 
                    positions=positions, showcaps=False)
                set_box_colors(bp)
                # ax[i].set_xticklabels(x_ticks, fontname='Times New Roman', fontsize=8, rotation=45, 
                    # horizontalalignment='right', verticalalignment='top')
                ax[i].get_xaxis().set_visible(False)
                ax[i].set_axisbelow(True)
                ax[i].yaxis.grid(True, linestyle='-', linewidth='0.75', which='major', color='white')
                ax[i].yaxis.grid(True, linestyle='--', which='minor', color='white',  alpha=0.3, linewidth='0.75')
                ax[i].tick_params(axis='x', which='minor', bottom=False)
                ax[i].set_facecolor('gainsboro')
                ax[i].set_title(x_labels[i], fontname='Times New Roman', fontsize=7)
                for j in range(0,len(x_ticks)): # add sample size data
                    ax[i].text(j, -5.25, num_samples[j], horizontalalignment='center', fontsize=5, fontname='Times New Roman')
                if count > 1: # remove y-ticks from other subplots
                    ax[i].tick_params(axis='y', which='both', left=False)

            fig.legend(handles=legend_ele, framealpha=1.0, loc='lower center', ncol=len(x_ticks), fontsize='xx-small', columnspacing=0.5)
            plt.yticks(fontname = 'Times New Roman', fontsize=10)
            plt.savefig(figure_output, bbox_inches = 'tight')
            plt.close()

def convert_to_kbp(filename, dir, name):
    '''Converts LDhat output files so that the loci are in terms of kbp instead of bp. 
    LDhot requires loci to be in kbp as it starts by constructing overlapping 3 kb windows 
    across a region and the 'stepping' (--step) 1 kb windows within to look for hotspots. If 
    these loci are not altered to kbp then LDhot still assumes kbp, but in reality you will be 
    creating overlapping windows of 3 bp and "stepping" 1 bp at a time. '''
    new_file = os.path.join(dir, name)
    with open(filename, 'r') as in_file, open(new_file, 'w') as out_file:
        if name == 'locs_kbp.txt':
            for i in range(1):
                header = next(in_file).split()
                kbp_header = float(header[1])/1000
                out_file.write(' '.join([header[0],str(kbp_header), header[2]]) + '\n')
        else:
            for i in range(1): 
                out_file.write(next(in_file)) # write out header
        for line in in_file:
            cols = line.split()
            loci = int(float(cols[0]))
            kbp_loci = str(loci/1000)

            if name == 'sum_rates_res_kbp.txt':
                if float(kbp_loci) < 0: # special print out for 2nd line
                    kbp_loci = '-1.000'
                    out_file.write(kbp_loci+'\t'+'\t'.join(cols[1:])+'\n')
                else:
                    out_file.write(kbp_loci+'\t'+'\t'.join(cols[1:])+'\n')
            elif name == 'locs_kbp.txt':
                out_file.write(kbp_loci+'\n')

    return new_file

def run_LDhot_multiprocess(list_of_directories):
    '''Runs LDhot in parallel. Make sure to use the NON-threaded version 
    of LDhot or this might overload your system.
    
    LDhot is run on unfiltered LDhat results, per Stukenbrock and Dutheil 2018 
    and personal communications with the corresponding author Julien Dutheil. In addition, 
    if running on LDhat results post-processing you would have to alter the locs.file and the 
    sites files to remove the SNPs that were filtered out of the LDhat results. As the summarized 
    res.files from LDhat need to match both the locs.file and sites.files for LDhot.
    '''
    for dir in list_of_directories:
        dir_name = 'alingment_'+dir
        print('Running LDhot on {}'.format(dir_name))
        dir_path = os.path.join(chr_dir, dir_name)
        hot_out = os.path.join(dir_path, lk_prefix+'ldhot.hotspots.txt')
        hot_log = os.path.join(dir_path, lk_prefix+'ldhot.log')
        if not os.path.exists(hot_out) and not os.path.exists(hot_log):
            loc_file = os.path.join(dir_path, 'locs.txt')
            site_file = os.path.join(dir_path, 'sites.txt')
            lk_file = os.path.join(work_dir, lk_prefix+'new_lk.txt')
            res_file = os.path.join(dir_path, lk_prefix+'summarized_rates_res.txt')

            # convert locs.txt and summarized_rates_res.txt so that loci are in kbp instead of bp for LDhot
            kbp_res = convert_to_kbp(res_file, dir_path, 'sum_rates_res_kbp.txt')
            kbp_loc = convert_to_kbp(loc_file, dir_path, 'locs_kbp.txt')

            try:
                with open(os.path.join(dir_path, 'LDhot.log'), 'w') as logfile:
                    subprocess.call([LDHOT, '--seq', site_file, '--loc', kbp_loc, '--lk', lk_file, 
                    '--res', kbp_res, '--nsim', str(args.numbersims), '--out', lk_prefix+'ldhot',
                    '--seed', '1234', '--windist', '10'
                    ], cwd=dir_path, stdout=logfile, stderr=logfile)
            except:
                print(('There was an error in running LDhot, check the logfile located at {}')
                .format(os.path.join(dir_path, 'LDhot.log')))

        else:
            pass

def run_LDhot_sum():
    print('Summarizing LDhot results')
    for dirs in os.listdir(chr_dir):
        dir_path = os.path.join(chr_dir, dirs)
        sum_hot_log = os.path.join(dir_path, lk_prefix+'ldhot_.summary.log')
        sum_hot_out = os.path.join(dir_path, lk_prefix+'ldhot_.hot_summary.txt')
        if not os.path.exists(sum_hot_out) and not os.path.exists(sum_hot_log):
            kbp_res = os.path.join(dir_path, 'sum_rates_res_kbp.txt')
            hot_results = os.path.join(dir_path, lk_prefix+'ldhot.hotspots.txt')
            try:
                with open(os.path.join(dir_path, 'LDhot_summary.log'), 'w') as logfile:
                    subprocess.call([LDHOT_SUM, '--res', kbp_res, '--hot', hot_results, 
                    '--out', lk_prefix+'ldhot_0.001'
                    ], cwd=dir_path, stdout=logfile, stderr=logfile)
            except:
                print(('There was an error in running LDhot_summary, check the logfile located at {}')
                .format(os.path.join(dir_path, 'LDhot_summary.log')))
                sys.exit()
        else:
            pass

def merge_summarized_LDhot_results():
    print('Merging summarized LDhot results')
    merged_ldhot_out = os.path.join(work_dir, 'Full_merged_summarized_LDhot_results.tsv')
    filter_merged_ldhot_out = os.path.join(work_dir, 'Filtered_merged_summarized_LDhot_results.tsv')
    if os.path.exists(filter_merged_ldhot_out):
        filter_merged_ldhot_df = pd.read_csv(filter_merged_ldhot_out, delimiter='\t')
    else:
        merge_dict = OrderedDict()
        chr_name_list = []
        dir_prefix = None
        for dirs in os.listdir(chr_dir): # get list of chromosomes from files
            dir_path = os.path.join(chr_dir, dirs)
            chr_num = dirs.split('_')[1]
            dir_prefix = dirs.split('_')[0]
            chr_name_list.append(chr_num)
        # sort by chromosome number and subsets
        sorted_chr_name_list = sorted(chr_name_list, key=lambda x: (int(x.split('-')[0]),int(x.split('-')[1])))

        ldhot_header = []
        merged_list = []
        filtered_merged_list = []
        for chr in sorted_chr_name_list:
            sum_results_path = os.path.join(chr_dir, dir_prefix+'_'+chr, lk_prefix+'ldhot.hot_summary.txt')
            if os.path.exists(sum_results_path): # LDhot will refuse to run if <100 SNPs, so these files won't exist if that happened
                chr_num = str(chr.split('-')[0])
                chr_subset = str(chr.split('-')[1])
                tmp_list = []
                with open(sum_results_path, 'r') as sum_r:
                    for line in sum_r:
                        if line.startswith('#'):
                            ldhot_header = line.strip()
                        else: # add column of correct start positions for each loci from MafFilter output
                            start, end, p, rho, peak = line.strip().split('\t')
                            start_bp = float(start)*1000 #convert back to bp
                            end_bp = float(end)*1000 #convert back to bp
                            if int(chr_subset) == 1:
                                cor_start = str(int(start_bp) + int(chr_start_dict[chr_num]))
                                cor_end = str(int(end_bp) + int(chr_start_dict[chr_num]))
                                size = str(int(cor_end) - int(cor_start))
                            elif int(chr_subset) > 1:
                                cor_start = str(int(start_bp)+int(chr_start_dict[chr_num])+int(align_start_dict[chr]))
                                cor_end = str(int(end_bp)+int(chr_start_dict[chr_num])+int(align_start_dict[chr]))
                                size = str(int(cor_end) - int(cor_start))
                            merged_list.append([chr_num,cor_start,cor_end,start,end,p,rho,peak,size])
                            if 100 >= float(rho) >= 5 and int(size) <= 20000:
                                filtered_merged_list.append([chr_num,cor_start,cor_end,start,end,p,rho,peak,size])

        with open(merged_ldhot_out, 'w') as m_out, open(filter_merged_ldhot_out, 'w') as f_m_out:
            m_out.write('Chr\tCorrect_start\tCorrect_end\t{}\tHotSpot_size(bp)\n'.format(ldhot_header))
            for item in merged_list:
                m_out.write('\t'.join(item) + '\n')
            f_m_out.write('Chr\tCorrect_start\tCorrect_end\t{}\tHotSpot_size(bp)\n'.format(ldhot_header))
            for item in filtered_merged_list:
                f_m_out.write('\t'.join(item)  + '\n')
        filter_merged_ldhot_df = pd.read_csv(filter_merged_ldhot_out, delimiter='\t')

    return filter_merged_ldhot_df

if __name__ == "__main__":
    align_start_dict, sample_size, chr_start_dict, parallel_dict = check_inputs()
    lk_prefix = create_likelihood_table(sample_size)
    convert_fasta_files()
    run_interval()
    summarize_LDhat_results()
    merged_df, unfil_merged_df, snp_dict = merge_summarized_LDhat_results()
    full_bp_df, unfull_bp_df = create_full_genomic_range_LDhat(merged_df, unfil_merged_df)

    whole_genome_rho = weighted_avg(merged_df)
    print('Whole genomic recombination = ' + str(whole_genome_rho))

    if args.gff3:
        gff_dict, gene_dict = create_gff_dictionary(args.gff3)
        classified_df, summarized_df = extract_regions_LDhat(full_bp_df, gff_dict)
        summarized_region_out = os.path.join(result_dir, 'Summarized_regions_LDhat_results.tsv')
        summarized_df = pd.read_csv(summarized_region_out, delimiter='\t')

        if ''.join(args.features) == 'all':
            feat_list = ['Exon','Intron','Upstream','Downstream','Intergenic']
        else:
            feat_list = args.features

        if args.annotations:
            anno_dict = OrderedDict() # {Annotation Name : [list of genes]}
            gene_list = [g for g in gene_dict.keys()]
            filtered_gene_list = fill_annotation_dict(gene_list, args.annotations)
            anno_dict['Other Genes'] = filtered_gene_list

            rev_anno_dict = {} # Reverse anno_dict to {Gene name: Annotation type}
            for anno, genes in anno_dict.items():
                for gene in genes:
                    if gene not in rev_anno_dict.keys():
                        rev_anno_dict[gene] = anno
                    else: # this will also remove duplicate calls of annotation types based on order of input files
                        pass

            anno_feat_dict = {a : {f : [] for f in feat_list} for a in anno_dict.keys()}
            # Example:
            # anno_feat_dict = {'Effectors' : {'Exon' : [list of rhos], 'Intron' : [list of rhos]}, 'Secreted' :...}
            fill_anno_feat_dict(summarized_df, 'anno_feat')
            create_boxplots(anno_feat_dict)

        else:
            pass
            feat_dict = {f : [] for f in feat_list}
            fill_anno_feat_dict(summarized_df, 'feat')
            create_boxplots(feat_dict)

    else:
        pass

    if args.LDhot: # run LDhot on recombination map estimated from LDhat interval if called
        parallel_list = [v for v in parallel_dict.values()]
        pool = Pool(processes=args.cpus)
        results = pool.imap(run_LDhot_multiprocess, parallel_list)
        pool.close()
        pool.join()
        run_LDhot_sum()

        filter_merged_ldhot_df = merge_summarized_LDhot_results()

    else:
        pass
