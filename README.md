# GO analysis of GWAS hits

## Requirements:
- .faa file (or .fa files - modify interproscan command for this) of background genes (e.g. pangenome).
- A file with the names of the GWAS hits (these need to be corresponding to the pangenome names used for interproscan). This should be tab separated and the gene names should be found in a column called 'Gene'
- the required packages downloaded

## Stages
1. Run interscanpro on the background file

2. Run consolidate_interproscan.py to reformat the output tsv files (at some point this should probably just be included in the R script)

3. Run plotting_go.R in R


## Choice of background and gene hits
Background should proabably be the pangenome, but it would also work using a normal genome (this should probably the first reference genome used in pyseer).
