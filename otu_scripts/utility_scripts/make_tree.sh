#wget http://www.nature.com/article-assets/npg/nmicrobiol/2016/nmicrobiol201648/extref/nmicrobiol201648-s7.txt
# get this as well http://rfam.xfam.org/family/RF00177/cm and name it 16S_bacteria.cm 

#mv nmicrobiol201648-s7.txt ../data/hug_tol.fasta
#cd ../data
#conda create -n treemaker python=3.6
#conda install -c bioconda infernal seqmagick raxml pplacer taxtastic
source actvate treemaker

#tr "[ -%,;\(\):=\.\\\[]\"\']" "_" < hug_tol.fasta > hug_tol.clean.fasta
#seqmagick mogrify --ungap hug_tol.clean.fasta

cmalign --dna -o hug_tol.clean.align.sto --outformat Pfam 16S_bacteria.cm hug_tol.clean.fasta
seqmagick convert hug_tol.clean.align.sto hug_tol.clean.align.fasta
seqmagick mogrify --deduplicate-sequences hug_tol.clean.align.fasta
raxmlHPC-PTHREADS-AVX2 -T 8 -m GTRGAMMA -s hug_tol.clean.align.fasta -n ref.tre -f d -p 12345
raxmlHPC-PTHREADS-AVX2 -T 2 -m GTRGAMMA -f I -t RAxML_bestTree.ref.tre -n root.ref.tre
raxmlHPC-PTHREADS-AVX2 -T 8 -m GTRGAMMA -f J -p 12345 -t RAxML_rootedTree.root.ref.tre -n conf.root.ref.tre -s hug_tol.clean.align.fasta
taxit create -l 16S_rRNA -P hug_tol.refpkg --aln-fasta hug_tol.clean.align.fasta --tree-stats RAxML_info.ref.tre --tree-file RAxML_fastTreeSH_Support.conf.root.ref.tre

# write out a query file 
tr "[ -%,;\(\):=\.\\\[]\"\']" "_" < query.fasta > query.clean.fasta

## Concatenate the query reads and reference alignment
cat hug_tol.clean.align.fasta query.clean.fasta > query.hug_tol.clean.fasta

## Remove gaps
seqmagick mogrify --ungap query.hug_tol.clean.fasta

## Align
cmalign --dna -o query.hug_tol.clean.align.sto --outformat Pfam 16S_bacteria.cm query.hug_tol.clean.fasta

## Convert to fasta
seqmagick convert query.hug_tol.clean.align.sto query.hug_tol.clean.align.fasta
