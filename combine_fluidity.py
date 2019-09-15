#!/usr/bin/python3

import os, sys, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

rundir = os.getcwd()

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='./%(prog)s [options] -i XXX_pangenome_fluidity.txt XXX_pangenome_fluidity.txt -o output_folder',
    description = '''    Turns multiple outputs from genome_fluidity.py into 
    one figure. Make sure to attach a prefix "XXX_pangenome_fluidity.txt if you 
    didn't already use the -p / --prefix flag in the previous script, to be 
    used in the legend''',
    
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
    default = 8,
    type=int,
    help = 'Size of legend [default = 8]',
    metavar=''
)
args=parser.parse_args()

output_figure = os.path.abspath(os.path.join(rundir, args.output+'.png'))

x_len = []
data = []
curve_top = []
curve_bottom = []
legend_labels = []
for file in args.input:
    legend_labels.append(file.split('_')[0])
    fpath = os.path.abspath(os.path.join(rundir, file))
    df=pd.read_csv(fpath, delimiter='\t', header=0)
    x_labels = np.array(df['Genomes_Sampled'])
    x_len.append(max(x_labels))
    fluid_best = df['Fluidity'].iloc[len(df['Fluidity'])-1]
    fluid_data = np.array([fluid_best for x in range(0,len(df['Fluidity']))])
    data.append((x_labels, fluid_data))
    curve_top.append(np.array(df['Power_top']))
    curve_bottom.append(np.array(df['Power_bottom']))

fig, ax = plt.subplots(figsize=(args.figsize[0],args.figsize[1]), dpi=args.dpi)
colors = plt.cm.nipy_spectral([1.*i/len(args.input) for i in range(len(args.input))])
for i in range(0, len(data)):
    plt.plot(data[i][0], data[i][1], ls='--', lw=1.5, color=colors[i], alpha=1, label='{}'.format(legend_labels[i]))
    plt.fill_between(data[i][0], curve_top[i], curve_bottom[i], facecolor=colors[i], alpha=0.25)
ax.set_facecolor('gainsboro')
ax.set_axisbelow(True)
plt.minorticks_on()
ax.tick_params(axis='x', which='minor', bottom=False)
ax.grid(which='major', linestyle='-', linewidth='1', color='white')
plt.xticks(np.arange(3, max(x_len) + 1, 1.0))
plt.xlim(3, max(x_len))
plt.xlabel('Genomes sampled')
plt.ylabel('Fluidity, '+u'\u03C6')
plt.legend(framealpha=1.0, prop={'size': args.legend})
plt.tight_layout()
plt.savefig(output_figure)
