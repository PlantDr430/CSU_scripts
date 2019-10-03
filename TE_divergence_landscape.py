#!/usr/bin/python3

'''
Not much flexibility with this script as it is more of a quick and dirty script to 
examine the TE landscape. Can add flexibility if requested, however, the script should 
be pretty simple to follow and easy to alter to tailor it for your specific interests. 
'''

import os, sys, argparse
import numpy as np
import matplotlib.pyplot as plt

rundir = os.getcwd()

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='./%(prog)s [options] -d directory -o output',
    description = '''    Takes RepeatMasker .out file and plots a landscape of 
    transposable elements based on percent divergence and transposable element fragment 
    lengths.''',
    
    epilog = """Written by Stephen A. Wyka (2019)""",
    formatter_class = MyFormatter)

parser.add_argument(
    '-i',
    '--input',
    required = True,
    help = 'Repeatmasker output file',
    metavar=''
)
parser.add_argument(
    '-o',
    '--output',
    required = True,
    help = 'Base name of output figure file (will automatically append .png)',
    metavar=''
)parser.add_argument(
    '-xlim',
    '--xlimits',
    nargs = '+',
    default = [-0.5, 50],
    type=float,
    help = 'X-axis data range limits [default: -0.5 50]',
    metavar=''
)
)parser.add_argument(
    '-t',
    '--title',
    help = "Title to be placed above figure, i.e. Species name.(use 'quotes if using spaces')",
    metavar=''
)
args=parser.parse_args()

def parse_repeatmasker_output(input_file):
    # create dictionaries of keys as bins (0 to 100) with values set to 0.
    # this allows us to easily use stacked bar graphs and to quickly add up 
    # sequence lengths to their respective bins
    ltr_dict = {i: 0 for i in range(0,100)}
    dna_dict = {i: 0 for i in range(0,100)}
    line_dict = {i: 0 for i in range(0,100)}
    sine_dict = {i: 0 for i in range(0,100)}
    unclass_dict = {i: 0 for i in range(0,100)}
    with open(input_file, 'r') as infile:
        file_linearr = [line.strip().split() for line in infile][3:] # skip first 3 lines
        for feat in file_linearr:
            div = round(float(feat[1])) # bin %divergence by rounding {keys}
            length = (int(feat[6]) - int(feat[5]) + 1) # gen legnth of fragment, don't forget to add 1
            if 'ltr' in feat[10].lower(): # if the TE is an LTR
                ltr_dict[div] = ltr_dict[div] + length
            elif 'dna' in feat[10].lower(): # if the TE is a DNA
                dna_dict[div] = dna_dict[div] + length
            elif 'line' in feat[10].lower(): # if the TE is a LINE
                line_dict[div] = line_dict[div] + length
            elif 'sine' in feat[10].lower(): # if the TE is a SINE
                sine_dict[div] = sine_dict[div] + length
            elif 'unknown' in feat[10].lower(): # if the TE is unclassified
                unclass_dict[div] = unclass_dict[div] + length
            else:
                pass
    return ltr_dict,dna_dict,line_dict,sine_dict,unclass_dict

def create_stacked_bar():
    fig, ax = plt.subplots()

    # from dictionaries keys are x_values and values are y_values
    # can arrange in any order, just make sure 'bottom=' is correct for placement
    sine_x = np.array([x for x in sine_dict.keys()])
    sine_y = np.array([y for y in sine_dict.values()])
    plt.bar(sine_x, sine_y, edgecolor='black',width=1, color='orange', label='SINE')
    line_x = np.array([x for x in line_dict.keys()])
    line_y = np.array([y for y in line_dict.values()])
    plt.bar(line_x, line_y, edgecolor='black',width=1, color='fuchsia', label='LINE', 
            bottom=sine_y)
    dna_x = np.array([x for x in dna_dict.keys()])
    dna_y = np.array([y for y in dna_dict.values()])
    plt.bar(dna_x, dna_y, edgecolor='black',width=1, color='blue', label='DNA', 
            bottom=sine_y+line_y)
    ltr_x = np.array([x for x in ltr_dict.keys()])
    ltr_y = np.array([y for y in ltr_dict.values()])
    plt.bar(ltr_x, ltr_y, edgecolor='black',width=1, color='red', label='LTR', 
            bottom=sine_y+line_y+dna_y)
    unclass_x = np.array([x for x in unclass_dict.keys()])
    unclass_y = np.array([y for y in unclass_dict.values()])
    plt.bar(unclass_x, unclass_y, edgecolor='black',width=1, color='green', label='Unclassified', 
            bottom=sine_y+line_y+dna_y+ltr_y)

    ax.set_axisbelow(True)
    plt.minorticks_on()
    plt.grid(which='minor', axis='y', color='white', linestyle='--', alpha=0.3)
    ax.yaxis.grid(True, linestyle='-', which='major', color='white')
    ax.xaxis.grid(True, linestyle='-', which='major', color='white')
    ax.tick_params(axis='x', which='minor', bottom=False)
    ax.set_facecolor('gainsboro')
    plt.xlim((args.xlimits[0],args.xlimits[1])) # 
    plt.legend(framealpha=1.0)
    plt.xlabel('Divergence (%)')
    plt.ylabel('Sequence length (bp)')
    if args.species:
        plt.title(args.species)
    else:
        pass # no title
    plt.tight_layout()
    # plt.show()
    plt.savefig(figure_output)
    plt.close()
    

if __name__ == "__main__":
    figure_output = os.path.abspath(os.path.join(rundir, args.output+'.png'))
    result_dicts = parse_repeatmasker_output(os.path.abspath(args.input))
    ltr_dict = result_dicts[0]
    dna_dict = result_dicts[1]
    line_dict = result_dicts[2]
    sine_dict = result_dicts[3]
    unclass_dict = result_dicts[4]
    create_stacked_bar()