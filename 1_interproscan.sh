#!/bin/bash
#SBATCH --job-name=interproscan
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --partition=short
#SBATCH --time=6:00:00
#SBATCH --cluster=arc
#SBATCH --mail-user=trop0670@ox.ac.uk
#SBATCH --mail-type=all

module purge
module load Anaconda3/2023.09-0
source activate /data/biol-micro-genomics/trop0670/env/interproscan

# If using nucleotide sequence change the -t section from p to n
bash /data/biol-micro-genomics/trop0670/interproscan/interproscan-5.72-103.0/interproscan.sh -b pangenome_interpro -cpu 8 -f TSV -goterms -i /data/biol-micro-genomics/trop0670/241015_campy_gwas/gwas/chicken_vs_wb_2/pyseer_outputs/interproscan2/pan_genome_reference.faa -iprlookup -pa -t p
