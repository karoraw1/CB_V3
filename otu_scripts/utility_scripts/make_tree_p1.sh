#!/bin/bash

#SBATCH
#SBATCH --job-name=mktree
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --exclusive
#SBATCH --partition=shared
#SBATCH --mail-type=END
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --error=../Logs/mktree.err
#SBATCH --output=../Logs/mktree.out

wget http://www.nature.com/article-assets/npg/nmicrobiol/2016/nmicrobiol201648/extref/nmicrobiol201648-s7.txt
# get this as well http://rfam.xfam.org/family/RF00177/cm and name it 16S_bacteria.cm 

mv nmicrobiol201648-s7.txt ../data/hug_tol.fasta
tr "[ -%,;\(\):=\.\\\[]\"\']" "_" < hug_tol.fasta > hug_tol.clean.fasta
seqmagick mogrify --ungap hug_tol.clean.fasta
cmalign --dna -o hug_tol.clean.align.sto --outformat Pfam 16S_bacteria.cm hug_tol.clean.fasta

source activate otu_caller
cd /home-3/karoraw1@jhu.edu/scratch/ChesBayTransect/data/TREEs
seqmagick convert hug_tol.clean.align.sto hug_tol.clean.align.fasta
seqmagick mogrify --deduplicate-sequences hug_tol.clean.align.fasta
raxmlHPC-PTHREADS -T 20 -m GTRGAMMA -s hug_tol.clean.align.fasta -n ref.tre -f d -p 12345
raxmlHPC-PTHREADS -T 20 -m GTRGAMMA -f I -t RAxML_bestTree.ref.tre -n root.ref.tre
raxmlHPC-PTHREADS -T 20 -m GTRGAMMA -f J -p 12345 -t RAxML_rootedTree.root.ref.tre -n conf.root.ref.tre -s hug_tol.clean.align.fasta
taxit create -l 16S_rRNA -P hug_tol.refpkg --aln-fasta hug_tol.clean.align.fasta --tree-stats RAxML_info.ref.tre --tree-file RAxML_fastTreeSH_Support.conf.root.ref.tre
