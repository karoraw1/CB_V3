#!/bin/bash

#SBATCH
#SBATCH --job-name=mktreel
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --exclusive
#SBATCH --partition=shared
#SBATCH --mail-type=END
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --error=../Logs/mktree2l.err
#SBATCH --output=../Logs/mktree2l.out

BASE_DIR=/media/login/KeithSSD/CB_V3
TREE_DIR=$BASE_DIR/otu_data/tree_data
SCRIPTS_=$BASE_DIR/otu_scripts/utility_scripts
TABLE_DIR=$BASE_DIR/otu_data/dada2_outputs
source activate otu_caller

# Write query file into tree folder
python $SCRIPTS_/pca_samples.py $TREE_DIR_ create_fasta $TABLE_DIR abundance_table.tsv taxa_table.tsv

# check query to see everything is a 16S gene
tr "[ -%,;\(\):=\.\\\[]\"\']" "_" < $TREE_DIR/query.fasta > $TREE_DIR/query.clean.fasta
cat $TREE_DIR/hug_tol.clean.align.fasta $TREE_DIR/query.clean.fasta > $TREE_DIR/query.hug_tol.clean.fasta
seqmagick mogrify --ungap $TREE_DIR/query.hug_tol.clean.fasta
# may produce error message, but will still succeed

cmsearch --cpu 7 --tblout $TREE_DIR/cm_report.txt --noali -o $TREE_DIR/cm_stdout.txt $TREE_DIR/16S_bacteria.cm $TREE_DIR/query.hug_tol.clean.fasta

COV_MOD_DIR=$TREE_DIR/cov_model_data
mkdir -p $COV_MOD_DIR

python $SCRIPTS_/read_cmsearch_report.py $COV_MOD_DIR cm_report.txt $TREE_DIR

mv $TREE_DIR/cm_report.txt $COV_MOD_DIR
mv $TREE_DIR/cm_stdout.txt $COV_MOD_DIR

python $SCRIPTS_/pca_samples.py $TREE_DIR_ cm_data $TABLE_DIR abundance_table.tsv taxa_table.tsv



