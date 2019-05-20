# coding=utf-8
import pandas as pd
import numpy as np
import os, sys

pd.set_option('mode.chained_assignment', None)

# This is the old data sheet with locations of old seq files
def load_sample_sheet(sample_sheet_path):
    ssu_data_fn = sample_sheet_path
    ssu_df = pd.read_csv(ssu_data_fn, sep="\t")
    ssu_df0 = ssu_df[ssu_df.loc[:, 'Sampling date (MMDDYY)'].notnull()]
    ssu_df0.loc[:, "new_date"] = ssu_df0.loc[:, 'Sampling date (MMDDYY)'].apply(lambda x: x[-2:]+x[:2]+x[2:-2])
    nullVal = ssu_df0[ssu_df0.DepthName.isnull()].DepthName.values[0]
    ssu_df0.loc[ssu_df0.DepthName.isin(['PC', 'NC', 'EB', 'SB']), 'DepthName'] = nullVal
    
    # This is a new data sheet with locations of new seq files
    new_ssu_fn = '../data/LibraryPrepProtocols/ThisRunToSampleProcessing.tsv'
    n_ssu_df = pd.read_csv(new_ssu_fn, sep="\t")
    n_ssu_df.loc[n_ssu_df.Station.isnull(), 'Station'] = "NA"
    
    # edit dates and depths of ctrl samples
    ctrl_samps = n_ssu_df.Exp_Ctrl_NA == 'CTRL'
    n_ssu_df.loc[(ctrl_samps) & (n_ssu_df.Date.isnull()) , 'Date'] = 181101.0
    n_ssu_df.loc[(ctrl_samps) & (~n_ssu_df.Station.str.contains("CB")) , 'Station'] = "LAB"
    
    # add these to the sample list
    dateify = lambda x: str(x)[2:-4]+str(x)[-4:-2]+str(x)[:2]
    end_of_sheet = ssu_df0.index.max()
    for c_idx, c_ix in enumerate(n_ssu_df[ctrl_samps].index):
        ssu_ix = end_of_sheet+1+c_idx
        n_ssu_df.IndexMatch[c_ix] = ssu_ix
        c_df = pd.DataFrame(index=[ssu_ix], columns=ssu_df0.columns)
        c_df.loc[ssu_ix, 'Short sample name'] = n_ssu_df.loc[c_ix, 'Sample Name']
        c_df.loc[ssu_ix, 'SampleID'] = n_ssu_df.loc[c_ix, 'Sample Name']
        c_df.loc[ssu_ix, 'DateMMDDYY'] = n_ssu_df.Date.apply(dateify)[c_ix]
        c_df.loc[ssu_ix, 'Sampling date (MMDDYY)'] = '110118'
        c_df.loc[ssu_ix, 'new_date'] = str(int(n_ssu_df.Date[c_ix]))
        c_df.loc[ssu_ix, 'StationName'] = n_ssu_df.Station[c_ix]
        c_df.loc[ssu_ix, 'DepthName'] = n_ssu_df.Depth[c_ix]
        ssu_df0 = ssu_df0.append(c_df, verify_integrity=True)
    
    # create barcode:well lookups
    bNw = set(n_ssu_df.loc[:, ['Barcode Well', 'RCBarcode']].apply(tuple, axis=1).tolist())
    well2bcode = {i.replace("-", ""): j for i, j in bNw}
    bcode2well = {j: i for i, j in well2bcode.items()}
    # remove controls & other peoples data
    n_ssu_df_exp = n_ssu_df[n_ssu_df.Exp_Ctrl_NA.isin(['Exp', 'CTRL'])]
    n_ssu_df_exp.Date = n_ssu_df_exp.Date.apply(int).apply(str)
    n_ssu_df_exp.IndexMatch = n_ssu_df_exp.IndexMatch.apply(int)
    n_ssu_df_exp.Date = n_ssu_df_exp.Date.apply(lambda x : x.replace('170410', '170411'))
    n_ssu_df_exp.loc[n_ssu_df_exp.Depth == 'Surface', 'Depth'] = 'Surf'
    # all of the stations are matched between sheets!
    # april 2017 samples are all listed as 04/11, 
    new_samps = n_ssu_df_exp.loc[:, ['Date', 'Station', 'Depth', 'IndexMatch']].apply(tuple, axis=1)
    
    for ix in new_samps.index:
        #print "Matched new lib #", ix, new_samps[ix]
        date_, stat_, depth_, ix_2 = new_samps[ix]
        assert ssu_df0.new_date[ix_2] == date_
        assert ssu_df0.StationName[ix_2] == stat_
        if depth_ == 'Surf':
            assert 'urface' in ssu_df0.loc[ix_2, 'Short sample name']
        elif depth_ == 'Deep':
            assert 'ottom' in ssu_df0.loc[ix_2, 'Short sample name']
        elif str(depth_) == str(nullVal):
            pass
        else:
           assert int(ssu_df0.DepthName[ix_2]) == int(depth_)
        
        ssu_df0.loc[ix_2, 'sequencing date'] = '01/02/2019'
        ssu_df0.loc[ix_2, 'sequencing ID'] = n_ssu_df_exp.Run[ix]
        ssu_df0.loc[ix_2, '2nd step well'] = n_ssu_df_exp.loc[ix, 'Barcode Well'].replace("-", "")
        ssu_df0.loc[ix_2, '2nd step barcode sequence'] = n_ssu_df_exp.RCBarcode[ix]
    
    # with new files added, we can remove anything without a sequencing date
    ssu_df1 = ssu_df0[ssu_df0.loc[:, 'sequencing date'].notnull()]
    # we can ignore any pointer to old files that needed reseqeuncing 
    ssu_df2 = ssu_df0[ssu_df0.loc[:, 'Resequencing files'].notnull()]
    ssu_df1.loc[ssu_df2.index, 'sequencing ID'] = ssu_df2.loc[ssu_df2.index, 'Resequencing files']
    # rows with "KAW" as a sequencing ID are shotgun metagenomes and should be removed
    ssu_df3 = ssu_df1[ssu_df1.loc[:, 'sequencing ID'] != 'KAW']
    # the barcode well is mis-entered on the sheet, confirmed in mapping file
    ssu_df3.loc[515, '2nd step well'] = 'D10'
    ssu_df3['Demux_Bool'] = pd.Series(index=ssu_df3.index, data=[False]*ssu_df3.shape[0])
    #ssu_df3.loc[ssu_df3.loc[:, 'sequencing ID'] == 'Miseq_data_SarahPreheim_Sept2016', 'Demux_Bool'] = True
    # TODO: make detection here automatic
    ssu_df3.loc[:, 'Demux_Bool'] = True
    
    for idx, ix in enumerate(ssu_df3.index):
        this_bcode = str(ssu_df3.loc[ix, '2nd step barcode sequence'])
        this_well = ssu_df3.loc[ix, '2nd step well']
        
        if this_well not in well2bcode.keys() and this_bcode in bcode2well.keys():
            ssu_df3.loc[ix, '2nd step well'] = bcode2well[this_bcode]
        elif this_well in well2bcode.keys() and this_bcode not in bcode2well.keys():
            ssu_df3.loc[ix, '2nd step barcode sequence'] = well2bcode[this_well]
        elif this_well in well2bcode.keys() and this_bcode in bcode2well.keys():
            pass
        elif this_bcode not in ['MISSING', 'nan']:
            pass
        elif ssu_df3.loc[ix,'Demux_Bool']:
            pass
        else:
            print("{} is an orphan {}, {} w/o a bcode".format(ix, ssu_df3.loc[ix, 'Short sample name'], ssu_df3.loc[ix, 'sequencing ID']))
            # TODO: Change this to a warning
            raise Exception('A library record without a discernable barcode was detected')
    
    ssu_dfx = ssu_df3.copy()
    
    for ix in ssu_dfx.index:
        id_ = ssu_dfx.loc[ ix, 'SampleID']
        date_ = ssu_dfx.loc[ ix, 'DateMMDDYY']
        sdate_ = ssu_dfx.loc[ ix, 'Sampling date (MMDDYY)']
        if date_ not in id_ and id_.startswith("SB"):
            new_sid = id_[:2]+date_+id_[8:]
            #print ix, id_, "->", new_sid
            ssu_dfx.loc[ix, 'SampleID'] = new_sid
    
    ssu_dfx.loc[ssu_dfx.StationName == 'CB54_1', 'StationName'] = 'CB54'
    
    return ssu_dfx
    
