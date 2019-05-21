#!/bin/bash

#SBATCH
#SBATCH --job-name=mktree3
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=800G
#SBATCH --exclusive
#SBATCH --partition=lrgmem
#SBATCH --mail-type=END
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --error=../Logs/mktree32.err
#SBATCH --output=../Logs/mktree32.out

LOC_DIR=/home-3/karoraw1@jhu.edu/scratch/ChesBayTransect/data/TREEs

source activate treemaker
cd $LOC_DIR
QUERY=query_CMSearched

#tr "[ -%,;\(\):=\.\\\[]\"\']" "_" < $QUERY.fasta > $QUERY.clean.fasta
#cat hug_tol.clean.align.fasta $QUERY.clean.fasta > $QUERY.hug_tol.clean.fasta
#seqmagick mogrify --ungap $QUERY.hug_tol.clean.fasta
#cmalign --dna -o $QUERY.hug_tol.clean.align.sto --outformat Pfam 16S_bacteria.cm $QUERY.hug_tol.clean.fasta
#seqmagick convert $QUERY.hug_tol.clean.align.sto $QUERY.hug_tol.clean.align.fasta

pplacer -o $QUERY.hug_tol.clean.align.jplace -p -c hug_tol.refpkg $QUERY.hug_tol.clean.align.fasta
guppy to_csv --point-mass --pp -o $QUERY.hug_tol.clean.align.csv $QUERY.hug_tol.clean.align.jplace
guppy fat --node-numbers --point-mass --pp -o $QUERY.hug_tol.clean.align.phyloxml $QUERY.hug_tol.clean.align.jplace
guppy distmat -o $QUERY.hug_tol.clean.align.dist.tab $QUERY.hug_tol.clean.align.jplace 
