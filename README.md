# CSU_scripts
#### Scripts created at CSU as a PhD student for genomic annotations and comparisons.

## Examples of figure output from specified scripts.

### genome_speed_hexbins.py
Takes a Gff3 or Bed file to produce a density plot map of intergenic regions lengths (kbp) like so:

![Hexbins](https://github.com/PlantDr430/images/blob/master/metabolite_speed.png)


### genome_fluidity.py
Takes an orthogroups.txt file from progrmas such as Orthofinder, OrthoMCL, SiLiX and calcualtes genome fluidity of a group of at least 5 genomes. Plots genome fluidity and standard error regions fit to an exponential regression curve y = Ae^Bx + C. Also reports a text file to be used with combine_fluidity.py.

![Single_fluidity](https://github.com/PlantDr430/images/blob/master/Species_9_pangenome_fluidity.png)

### combine_fluidity.py
Takes multiple text output files from genome_fluidity.py and merges them into one graph. Also report tab-deliminated file of each v. each p-values calculated from a two sample two-sided z-test. Note: This is mock data.

![Combined_fluidity](https://github.com/PlantDr430/images/blob/master/Example_combined_fluidity.png)
