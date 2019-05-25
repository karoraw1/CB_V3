#!/bin/bash

#SBATCH --job-name=mt3g
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=500G
#SBATCH --partition=lrgmem
#SBATCH --exclusive
#SBATCH --mail-type=END 
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --error=mt33p.err
#SBATCH --output=mt33p.out

LOC_DIR=~/scratch/CB_V3/otu_data/tree_data 
NEW_TD=$LOC_DIR/full_tree

mkdir -p $NEW_TD

cd $NEW_TD

QUERY=query_cmsearched

source activate otu_caller

tr "[ -%,;\(\):=\.\\\[]\"\']" "_" < ../$QUERY.fasta > $QUERY.clean.fasta

cat ../hug_tol.clean.align.fasta $QUERY.clean.fasta > $QUERY.hug_tol.clean.fasta

seqmagick mogrify --ungap $QUERY.hug_tol.clean.fasta

cmalign --mxsize 4092 --cpu 24 --dna -o $QUERY.hug_tol.clean.align.sto --sfile cm_scores.txt --outformat Pfam ../16S_bacteria.cm $QUERY.hug_tol.clean.fasta

seqmagick convert $QUERY.hug_tol.clean.align.sto $QUERY.hug_tol.clean.align.fasta

pplacer -o $QUERY.hug_tol.clean.align.jplace -p -c ../hug_tol.refpkg $QUERY.hug_tol.clean.align.fasta

guppy to_csv --point-mass --pp -o $QUERY.hug_tol.clean.align.csv $QUERY.hug_tol.clean.align.jplace

guppy fat --node-numbers --point-mass --pp -o $QUERY.hug_tol.clean.align.phyloxml $QUERY.hug_tol.clean.align.jplace

guppy distmat -o $QUERY.hug_tol.clean.align.dist.tab $QUERY.hug_tol.clean.align.jplace

#raxmlHPC-PTHREADS -T 24 -m GTRGAMMA -s $QUERY.hug_tol.clean.align.fasta -n ref.tre -f d -p 12345 

#raxmlHPC-PTHREADS -T 24 -m GTRGAMMA -f I -t RAxML_bestTree.ref.tre -n root.ref.tre

#raxmlHPC-PTHREADS -T 24 -m GTRGAMMA -f J -p 12345 -t RAxML_rootedTree.root.ref.tre -n conf.root.ref.tre -s $QUERY.hug_tol.clean.align.fasta

#raxmlHPC-PTHREADS -T 24 -f x -p 12345 -t RAxML_rootedTree.root.ref.tre -s $QUERY.hug_tol.clean.align.fasta -m GTRGAMMA -n $QUERY.distmat

