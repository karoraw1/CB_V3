# coding=utf-8
import yaml, os, sys
config_file = "config.yml"

with open(config_file, 'r') as stream:
    cfg_dict = yaml.safe_load(stream)

data_dir = cfg_dict['data_directory']
sample_sheet_fn = cfg_dict['sample_sheet']

if not os.path.exists(cfg_dict['batch_dir']):
    os.mkdir(cfg_dict['batch_dir'])

import pandas as pd
import numpy as np

pd.set_option('mode.chained_assignment', None)

ssu_df4 = pd.read_csv(sample_sheet_fn, sep="\t")

for i in ["I", "R1", "R2"]:
    ssu_df4["BulkDataFile_{}".format(i)] = pd.Series(index=ssu_df4.index, data=["None"]*ssu_df4.shape[0])

ssu_df4["MappingFile"] = pd.Series(index=ssu_df4.index, data=["None"]*ssu_df4.shape[0])
ssu_df4["DemuxFile_R1"] = pd.Series(index=ssu_df4.index, data=["None"]*ssu_df4.shape[0])
ssu_df4["DemuxFile_R2"] = pd.Series(index=ssu_df4.index, data=["None"]*ssu_df4.shape[0])

for s_id in ssu_df4.loc[:, 'sequencing ID'].unique():
    row_set = ssu_df4.loc[:, 'sequencing ID'] == s_id
    # sequences dir
    s_id_dir = data_dir+"/"+s_id
    seq_dir = s_id_dir+"/"+"FASTQ"
    # check to see if demuxed
    # TODO: Change Demux detection to automatic
    dmux = ssu_df4.loc[row_set, "Demux_Bool"]
    dmux_q = list(set(dmux.values))
    assert len(dmux_q) == 1
    print("\n{}\n--\n\tAre libs are demuxed: {}".format(s_id, dmux_q[0]))
    if dmux_q[0]:
        demux_dir = seq_dir + "/Demux"
        ssu_df4.loc[row_set, "DemuxFile_R1"] = ssu_df4.loc[row_set, 'SampleID'].apply(lambda x: demux_dir+"/"+x+".R1.fastq")
        ssu_df4.loc[row_set, "DemuxFile_R2"] = ssu_df4.loc[row_set, 'SampleID'].apply(lambda x: demux_dir+"/"+x+".R2.fastq")
        for d_ix in ssu_df4[row_set].index:
            assert os.path.exists(ssu_df4.loc[d_ix, "DemuxFile_R1"])
            assert os.path.exists(ssu_df4.loc[d_ix, "DemuxFile_R2"])
        print("\tAll {} forward and reverse read files exist!".format(row_set.sum()*2))
    
    # this locates the sequence files, unique to our lab
    for k in os.listdir(seq_dir):
        if "Undetermined" in k and k.endswith("fastq"):
            bdf = os.path.join(seq_dir,k)
            if "L001_I1_001" in k:
                ssu_df4.loc[row_set, "BulkDataFile_I"] = bdf
                print("\tBulk index file identified")
            elif "L001_R1_001" in os.path.basename(bdf):
                ssu_df4.loc[row_set, "BulkDataFile_R1"] = bdf
                print("\tBulk R1 file identified")
            elif "L001_R2_001" in os.path.basename(bdf):
                ssu_df4.loc[row_set, "BulkDataFile_R2"] = bdf
                print("\tBulk R2 file identified")
    
    # this creates an index file of barcodes
    subdf = ssu_df4[row_set].loc[:, ["SampleID", '2nd step barcode sequence']].copy()
    subdf.columns = ["SampleName", "RCBarcode"]
    bcode_fn = s_id_dir+"/"+s_id+"_barcodes.txt"
    subdf.to_csv(bcode_fn, sep="\t", index=False)
    print("\tWrote {} barcodes to \n\t\t{}".format(subdf.shape[0], bcode_fn))
    for pip_i_f in ["demux_skeleton.sh", "filter_skeleton.sh", "callOTUs_skeleton.sh"]:
        pip_i = os.path.join("utility_scripts", pip_i_f)
        task_ = pip_i_f.split("_")[0]
        pipe_fn = cfg_dict['batch_dir']+"/"+s_id+"_"+task_+".sh"
        with open(pip_i, 'r') as fh:
            content = fh.read()
        
        n_cont = content.replace("@SID@", s_id)
        with open(pipe_fn , 'w') as ofh:
            ofh.write(n_cont)
        
        print("\tWrote {}  to\n\t\t{}".format(pip_i, pipe_fn))
