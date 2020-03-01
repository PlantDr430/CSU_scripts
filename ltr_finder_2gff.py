#!/usr/bin/python3

import os, sys, re, argparse
from collections import defaultdict

rundir = os.getcwd()
parentdir = os.path.dirname(rundir)

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='./%(prog)s [options] -in LTF_finder_output -o output.gff',
    description = '''Turns LTR_finder output into GFF3 format.''',
    
    epilog = """Written by Stephen A. Wyka (2020)""",
    formatter_class = MyFormatter)

parser.add_argument(
    '-in',
    '--input',
    required=True,
    help = 'LTF_finder output file',
    metavar=''
)
parser.add_argument(
    '-o',
    '--output',
    required=True,
    help = 'Output GFF file',
    metavar=''
)
args=parser.parse_args()

def parse_file_to_dictionary(fpath):
    '''Break up LTR_finder output into a dictionary of found LTR's based on contig 
    and subset (i.e. [1], [2], etc). The contigs and subset act as keys and the lines 
    associated with the LTR are captured as a nested list for parsing.
    '''
    line_breaks = defaultdict(list) # {<contig_ID>-<subset_number> : [Nested list of lines from LTR_output]}
    with open(fpath, 'r') as f_in:
        contig_sub = ''
        for line in f_in:
            match_start = re.match(r'(\[\d+\])\s(.+)\s(.+)', line)
            if match_start != None:
                contig = match_start.group(2)
                sub = ' '.join(re.findall(r'\[(\d+)\]',match_start.group(1)))
                contig_sub = contig+'-'+sub
                continue
            line_breaks[contig_sub].append([line])
    line_breaks.pop('')
    return line_breaks

def dictionary_to_gff(d):
    gff_list = []
    for match, lines in d.items():
        contig = match.split('-')[0]
        subset = match.split('-')[1]
        prog = 'ltr_finder'
        score = ' '.join(re.findall(r'\:\s(\d+)\s\[',lines[1][0]))
        phase = '.'
        strand = re.search(r'(Strand:)(.)',lines[0][0]).group(2)
        # id = contig+'_'+prog+'_model0001'
        parent = contig+'_'+prog+'_model000{}'.format(subset)
        for l in lines:
            l = l[0] # turn list into string
            if l.startswith('Location'):
                sect = re.search(r'(:)\s(\d+)\s-\s(\d+)', l)
                start = sect.group(2)
                stop = sect.group(3)
                attribute = 'ID={}'.format(parent)
                gff_line = [contig, prog, 'LTR_retrotransposon', start, stop, score, strand, phase, attribute]
                gff_list.append(gff_line)
            elif l.startswith("5'-LTR"):
                sect = re.search(r'(:)\s(\d+)\s-\s(\d+)', l)
                start = sect.group(2)
                stop = sect.group(3)
                attribute = 'ID={};Name=Five Prime LTR;Parent={}'.format(parent+'_five_prime_LTR',parent)
                gff_line = [contig, prog, 'five_prime_LTR', start, stop, score, strand, phase, attribute]
                gff_list.append(gff_line)
            elif l.startswith("3'-LTR"):
                sect = re.search(r'(:)\s(\d+)\s-\s(\d+)', l)
                start = sect.group(2)
                stop = sect.group(3)
                attribute = 'ID={};Name=Three Prime LTR;Parent={}'.format(parent+'_three_prime_LTR',parent)
                gff_line = [contig, prog, 'three_prime_LTR', start, stop, score, strand, phase, attribute]
                gff_list.append(gff_line)
            elif l.startswith('TSR'):
                if not 'NOT FOUND' in l: # if target site duplications are found
                    sect = re.search(r'(:)\s(\d+)\s-\s(\d+)\s(,)\s(\d+)\s-\s(\d+)', l)
                    start1 = sect.group(2)
                    stop1 = sect.group(3)
                    start2 = sect.group(5)
                    stop2 = sect.group(6)
                    attribute = 'ID={};Name=Target Site Duplication;Parent={}'.format(parent+'_tsd5',parent)
                    gff_line = [contig, prog, 'target_site_duplication', start1, stop1, score, strand, phase, attribute]
                    gff_list.append(gff_line)
                    attribute = 'ID={};Name=Target Site Duplication;Parent={}'.format(parent+'_tsd3',parent)
                    gff_line = [contig, prog, 'target_site_duplication', start2, stop2, score, strand, phase, attribute]
                    gff_list.append(gff_line)
                else:
                    pass
            elif l.startswith('PPT'):
                sect = re.search(r'(:)\s(\[.+/.+\])\s(\d+)\s-\s(\d+)', l)
                start = sect.group(3)
                stop = sect.group(4)
                attribute = 'ID={};Name=PPT;Parent={}'.format(parent+'_RR_tract',parent)
                gff_line = [contig, prog, 'RR_tract', start, stop, score, strand, phase, attribute]
                gff_list.append(gff_line)
            elif l.startswith('Domain'):
                sect = re.search(r'(:)\s(\d+)\s-\s(\d+)', l)
                start = sect.group(2)
                stop = sect.group(3)
                attribute = 'ID={};Name=Reverse Transcriptase;Parent={}'.format(parent+'_rvt',parent)
                gff_line = [contig, prog, 'transposable_element_gene', start, stop, score, strand, phase, attribute]
                gff_list.append(gff_line)
            else:
                pass
    return gff_list

def write_out_results(opath):
    with open(opath, 'w') as f_out:
        f_out.write('##gff-version 3\n')
        for line in gff_list:
            f_out.write('\t'.join(line) + '\n')

if __name__ == "__main__":
    input_file = os.path.abspath(args.input)
    out_file = os.path.abspath(args.output)
    line_breaks = parse_file_to_dictionary(input_file)
    gff_list = dictionary_to_gff(line_breaks)
    write_out_results(out_file)
