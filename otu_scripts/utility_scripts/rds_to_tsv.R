library(dada2)

# command line arguments
args <- commandArgs(TRUE)

# read in location
data_path <- args[1]
# data_path <- '/home-3/karoraw1@jhu.edu/scratch/CB_V3/otu_data/dada2_outputs'
taxa_file <- args[2]
# taxa_file <- 'tax_sp_final.rds'
abund_file <- args[3]
# abund_file = 'seqtab_final.rds'

# read in taxa table
t_path <- file.path(data_path, taxa_file)
taxa <- readRDS(t_path)

# read in abund table
a_path <- file.path(data_path, abund_file)
seqtab <- readRDS(a_path)

#  write out in human readable format
t_out <- file.path(data_path, 'taxa_table.tsv')
write.table(taxa, file=t_out, quote=FALSE, sep='\t', col.names = NA)

# write out in human readable format
a_out <- file.path(data_path, 'abundance_table.tsv')
write.table(seqtab, file=a_out, quote=FALSE, sep='\t', col.names = NA)


