import os, sys
from math import ceil
from Bio import SeqIO
import random
random.seed(42)

libs_dir, str_flag = sys.argv[-2:]
# libs_dir = '/scratch/groups/sprehei1/Keith_Files/Processed_data_group/esakows1_132789/FASTQ/Demux'
# str_flag = "R1.fastq"

files_ = [os.path.join(libs_dir, f) for f in os.listdir(libs_dir) if str_flag in f]

for fq in files_:
    fqp = SeqIO.parse(fq, "fastq")
    counter = sum((1 for rec in fqp))
    reads_per_lib = int(ceil(5e+5/len(files_)))
    # subsample should be less than counter
    reads_per_lib = min([counter, reads_per_lib])-1
    to_get = set(random.sample(range(counter), reads_per_lib))
    print("{}:{} reads".format(os.path.basename(fq), counter))
    out_file = fq+".sample"
    fqp = SeqIO.parse(fq, "fastq")
    wanted = (rec for idx, rec in enumerate(fqp) if idx in to_get)
    SeqIO.write(wanted, out_file, "fastq")
