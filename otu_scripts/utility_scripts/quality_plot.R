library(dada2)
args <- commandArgs(TRUE)
base_path <- args[1]
seq_ID <- args[2]
str_patt <- args[3]
lib_state <- args[4]

# library(dada2)
# base_path = "/home-3/karoraw1@jhu.edu/scratch/16S_Libraries"
# seq_ID = "esakows1_132789"
# sub_folder = "Demux" | "Trim"

path=file.path(base_path, seq_ID)
fnFs <- sort(list.files(path, pattern=str_patt, full.names = TRUE))
file_sizes_f = sort(sapply(fnFs, file.size), decreasing = T, na.last = NA)

PNGname = file.path(base_path, seq_ID, paste(lib_state, "R1_and_R2_quals.png", sep="_"))
to_plot=row.names(as.data.frame(file_sizes_f[c(1,2)]))
png(filename=PNGname)
plotQualityProfile(to_plot)
dev.off()

