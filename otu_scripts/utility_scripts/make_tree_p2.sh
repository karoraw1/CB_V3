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

LOC_DIR=/home-3/karoraw1@jhu.edu/scratch/ChesBayTransect/data/TREEs

source activate treemaker
cd $LOC_DIR
# check query to see everything is a 16S gene
tr "[ -%,;\(\):=\.\\\[]\"\']" "_" < query.fasta > query.clean.fasta.s
cat hug_tol.clean.align.fasta query.clean.fasta.s > query.hug_tol.clean.fasta.s
seqmagick mogrify --ungap query.hug_tol.clean.fasta.s
cmsearch --cpu 24 --tblout cm_report.txt.s --noali -o cm_stdout.txt.s 16S_bacteria.cm query.hug_tol.clean.fasta.s
