# Combine a phylogenetic tree with a locus coverage heatmap

### This script uses R to combine a presence absence matrix (1/0), and a NEWICK formated phylogenetic tree to generate a heatmap plotted against the tree. 

Typically, I've used this for examining coverage of UCE loci / target enrichment coverage by taxon. It can be used for anything that has a presence/absence matrix and a tree.  

For the example, I've pulled data from the [Phyluce](http://phyluce.readthedocs.io/en/latest/index.html) pipeline, after the [match\_contigs\_to\_probe](http://phyluce.readthedocs.io/en/latest/uce-processing.html#match-contigs-to-probes) script step.

**To get a presence/absence matrix from SQLite database file generated by script:**

```
sqlite3 probe3.matches.sqlite
.mode columns
.mode csv
.headers on
.nullvalue 0
.output matches_uce.csv
SELECT * FROM matches;
.exit
```

## Tree_Heatmap.r
This script will generate a file that looks like this: 
![Alt text](example_plots/Tree_Heatmap_Example.jpg?raw=true)

Inputs are:  
1. 1/0 presence absence matrix.  
2. Newik formatted tree file. 

**Adjust file names in the script as needed.**

Note: R is smart enough to still plot the tree even if the matrix is missing taxa. The missing taxa just wont show up on the plot. 

### Limitations: 
This can only generate a dendrogram style tree, so branch lengths are lost. 

## Missing_Plot.r

This was a first attempt at the problem. It will generate a heatmap, which is useful if you need to order the heatmap by presence etc., but does not include the option to add a tree.

It will generate plots like this: 
![Alt text](example_plots/Missing_Plot_Example.jpg?raw=true)

Input is:
1. 1/0 presence absence matrix. 

**Adjust file name in the script as needed.**