#!/bin/bash

#SBATCH
#SBATCH --job-name=dadaEnd
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --mem=960G
#SBATCH --exclusive
#SBATCH --partition=lrgmem
#SBATCH --mail-type=END
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --error=merge_taxa.err
#SBATCH --output=merge_taxa.out

module load R/3.5.1

# these will be filled in automatically
BASE_DIR=~/scratch/CB_V3
SCRIPTS_=$BASE_DIR/otu_scripts/utility_scripts
LIBdirDF=$BASE_DIR/otu_data/seq_tabs.csv
OUT_DIR=$BASE_DIR/otu_data/dada2_outputs
TAX_DB=~/scratch/DADA2_Silva_DBs

Rscript $SCRIPTS_/merge_chim_tax.R $TAX_DB $LIBdirDF $OUT_DIR
Rscript $SCRIPTS_/rds_to_tsv.R $OUT_DIR tax_sp_final.rds seqtab_final.rds
