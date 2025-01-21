# GO analysis of GWAS hits

## Requirements:
- .faa files (or .fa files - modify interproscan command for this) of background genes (e.g. pangenome) and genome hits.
- the required packages downloaded

## Stages
1. Run interscanpro on the 2 files (hits and background)

2. Run consolidate interproscan to reformat the output tsv files (it may also be possible to do this locally)

3. Run plotting_go.R in R

This could be made more efficient by getting the interproscan output directly processed in R, rather than requiring the additional processing step. However, it's good enough as it is.