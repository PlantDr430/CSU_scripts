# CSU_scripts
#### Scripts created at CSU as a PhD student for genomic annotations and comparisons.

## Examples of figure output from specified scripts.

### RIP_blast_analysis.py
Blast all v. all method to detect presence of RIP (Galagan et al. 2003, Urquhart et al. 2018).

![RIP_blast](https://github.com/PlantDr430/images/blob/master/rip_blast_examples.png)

### genome_speed_hexbins.py
Takes a Gff3 or Bed file to produce a density plot map of intergenic regions lengths (kbp) like so:

![Hexbins](https://github.com/PlantDr430/images/blob/master/speed_norm_dist_example.png)


### genome_fluidity.py
Takes an orthogroups.txt file from programs such as Orthofinder, OrthoMCL, SiLiX and calculates genome fluidity of a group of at least 5 genomes. Plots genome fluidity and standard error regions fit to an exponential regression curve y = Ae^Bx + C. Also reports a text file to be used with combine_fluidity.py.

![Single_fluidity](https://github.com/PlantDr430/images/blob/master/Pangenome_fluidity.png)

### combine_fluidity.py
Takes multiple text output files from genome_fluidity.py and merges them into one graph. Also report tab-deliminated file of each v. each p-values calculated from a two sample two-sided z-test. Note: This is mock data.

![Combined_fluidity](https://github.com/PlantDr430/images/blob/master/Example_combined_fluidity.png)

### TE_divergence_landscape.py
Takes RepeatMasker .out file and create a TE landscape based on percent divergence and sequence lengths of transposable element fragments. Note: This script does not have a lot of flexibility from the command line as it was mainly used for personal figure creation, however, the script should be simple enough to follow and alter it for your specific interests (i.e. TE classes to show, etc).

![TE_landscape](https://github.com/PlantDr430/images/blob/master/TE_divergence_landscape.png)
