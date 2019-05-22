import matplotlib.pyplot as plt
import matplotlib as mpl
from deicode.optspace import OptSpace
from deicode.preprocessing import rclr
import numpy as np
import os, sys
import pandas as pd
import seaborn as sns
from skbio.stats.ordination import pcoa
from skbio.stats.distance import mantel
from skbio.diversity import beta_diversity
from sklearn.decomposition import PCA
plt.style.use('seaborn-poster')


print("Reading OTU and taxa tables")
taxa_file = sys.argv[-1]
abund_file = sys.argv[-2]
data_path = sys.argv[-3]
task_name = sys.argv[-4]
out_path = sys.argv[-5]

#taxa_file = "taxa_table.tsv"
#abund_file = "abundance_table.tsv"
#data_path = "../../otu_data/dada2_outputs"
#task_name = "create_fasta"
#out_path = "../../otu_data/tree_data"

abund_f = os.path.join(data_path, abund_file)
tax_f = os.path.join(data_path, taxa_file)
taxa_df = pd.read_csv(tax_f, sep="\t")
abund_df = pd.read_csv(abund_f, sep="\t")

print("Renaming rows from the entire sequences to OTU# format")
print("\tStoring sequences in dictionary, accesible by OTU name")
OTU_Seqs = {taxa_df.loc[idx, taxa_df.columns[0]]:idx for idx in taxa_df.index}
OTU_Names = {idx:"OTU{}".format(idx+1) for idx in taxa_df.index }
OTU_name2seq = {OTU_Names[num]:seq for seq, num in OTU_Seqs.items()}
taxa_df.loc[:, taxa_df.columns[0]] = taxa_df.loc[:, taxa_df.columns[0]].apply(lambda x: OTU_Names[OTU_Seqs[x]])
taxa_df = taxa_df.set_index(taxa_df.columns[0])
new_cols = ['Samples']+[OTU_Names[OTU_Seqs[y]] for y in abund_df.columns[1:]]
abund_df.columns = new_cols
abund_df2 = abund_df.set_index('Samples')
# get rid of empty samples
abund_df3 = abund_df2[abund_df2.sum(1)!=0]

# export them for alignment to RFAM covariance model, then stop

if this_step == 'create_fasta':
    heads = sorted(list(abund_df3.columns))
    tails = [OTU_name2seq[i] for i in heads]
    with open(os.path.join(out_path, 'query.fasta'), "w") as wofh:
        wofh.write("".join([">{}\n{}\n".format(i, j) for i, j in zip(heads, tails)]))
elif this_step == 2:
    # read back in data
    with open("../data/TREEs/cmsearch_figs/poor_aligners.txt") as paFH:
        poor_aligners = [i for i in paFH.read().split("\n") if i != ""]
    with open("../data/TREEs/cmsearch_figs/reverse_strand_aligners.txt") as rsFH:
        rev_strand_algn = [i for i in rsFH.read().split("\n") if i != ""]

    # filter out poorly aligned seqs
    taxa_df = taxa_df[~taxa_df.index.isin(poor_aligners)]
    abund_df3 = abund_df3.loc[:, ~abund_df3.columns.isin(poor_aligners)]
    # reverse-complement flipped seqs
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    for r_otu in rev_strand_algn:
        rev_comp = str(Seq(OTU_name2seq[r_otu], generic_dna).reverse_complement())
        OTU_name2seq[r_otu] = rev_comp
    # write out good OTU sequences
    heads = sorted(list(abund_df3.columns))
    tails = [OTU_name2seq[i] for i in heads]
    with open('../data/TREEs/query_cmsearched.fasta', "w") as wofh:
        wofh.write("".join([">{}\n{}\n".format(i, j) for i, j in zip(heads, tails)]))
elif this_step == 3:
    with open("../data/TREEs/cmsearch_figs/poor_aligners.txt") as paFH:
        poor_aligners = [i for i in paFH.read().split("\n") if i != ""]

    with open("../data/TREEs/cmsearch_figs/reverse_strand_aligners.txt") as rsFH:
        rev_strand_algn = [i for i in rsFH.read().split("\n") if i != ""]

    # filter out poorly aligned seqs
    taxa_df = taxa_df[~taxa_df.index.isin(poor_aligners)]
    abund_df3 = abund_df3.loc[:, ~abund_df3.columns.isin(poor_aligners)]
    # reverse-complement flipped seqs
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    for r_otu in rev_strand_algn:
        rev_comp = str(Seq(OTU_name2seq[r_otu], generic_dna).reverse_complement())
        OTU_name2seq[r_otu] = rev_comp

    # read back distmat produced by GUPPY and pplacer info 
    dist_file = "../data/TREEs/query_CMSearched2.hug_tol.clean.align.dist.tab"
    data_file = "../data/TREEs/query_CMSearched2.hug_tol.clean.align.csv"
    
    # ready cluster frame
    known_subseqs = sorted(['OTU32682', 'OTU30', 'OTU2'])
    place_df = pd.read_csv(data_file, usecols=range(1,10))
    place_df['Seq_Len'] = place_df.name.apply(lambda x: len(OTU_name2seq[x]))
    place_df['OTU_Num'] = place_df.name.apply(lambda x: int(x[3:]))
    place_df.sort_values(by=['Seq_Len', 'OTU_Num'], inplace=True)
    place_df.set_index('name', inplace=True)

    # remove unique placements
    edge_counts = place_df.groupby('edge_num').size()
    solo_leaves = edge_counts[edge_counts == 1].index
    centroids = {place_df[place_df.edge_num == _e_].index[0]:[] for _e_ in solo_leaves}
    place_df_ns = place_df[~place_df.edge_num.isin(solo_leaves)]

    # create presence or absence df
    p_a_abund3 = (abund_df3 > 0) * 1.0
    
    # define fxn to calculate % of libraries in which one or the other OTU is observed, but not both. 
    def overlap_score(df):
        temp2 = df[df.sum(1)!=0].sum(1)
        temp3 = (temp2 == 1).sum() / temp2.shape[0]
        return temp3
    
    # function to measure homology %
    def ldist(seq1, seq2):
        ops = Levenshtein.editops(seq1, seq2)
        return 1 - (len(ops) / (len(seq1) + len([o for o in ops if o[0] == 'delete'])))
    
    # enumerate clusters 
    cl_n = place_df_ns.edge_num.unique().shape[0]
    
    # iterate over clusters absorbing matches
    for e_i, _e_ in enumerate(place_df_ns.edge_num.unique()):
        subdf = place_df_ns[place_df_ns.edge_num == _e_]
        processed = set()
        toproc_n = subdf.index.shape[0]
        for _otu_ in subdf.index:
            toproc = set(subdf.index) - processed
            if _otu_ in toproc:
                  # calculate distances within cluster
                  subseq_matches = [(i, ldist(OTU_name2seq[_otu_], OTU_name2seq[i])) for i in toproc if _otu_ != i]
                  # add overlap score and filter by 97.5% homology
                  filt_matches = [i[0] for i in subseq_matches if i[1]+overlap_score(p_a_abund3.loc[:,[_otu_, i[0]]]) > 1.975]
                  processed.update(filt_matches+[_otu_])
                  centroids[_otu_] = filt_matches
        print("{:.2%} complete. {} centroids now from {} recs in cluster {}".format((e_i+1)/cl_n, len(centroids.keys()), toproc_n, _e_))

    fixed_centroids = {}
    for otu_k, otu_cl in centroids.items():
        cluster_members = [otu_k] + otu_cl
        other_dists = sorted([(clm_i, abs(252 - place_df.Seq_Len[clm_i])) for clm_i in cluster_members], key=lambda x: x[1])
        new_center = other_dists[0][0]
        cluster_members.remove(new_center)
        fixed_centroids[new_center] = cluster_members

    abund_df4 = abund_df3.copy()
    for otu_fc, otu_fcl in fixed_centroids.items():
        if len(otu_fcl) != 0:
            abund_df4.loc[:, otu_fc] += abund_df4.loc[:, otu_fcl].sum(1)
            abund_df4.drop(otu_fcl, axis=1, inplace=True)

    # double check length distribution
    sLens, n_sLens = np.unique([len(OTU_name2seq[i]) for i in abund_df4.columns], return_counts=True)
    all_counts = abund_df4.sum().sum()
    wrong_lengths = set()
    for i, j in zip(sLens, n_sLens):
        these_idxs = place_df[place_df.Seq_Len == i].index
        matched_cols = abund_df4.columns[abund_df4.columns.isin(these_idxs)]
        if ((i < 250) or (i > 253)):
            wrong_lengths.update(matched_cols)
        else:
            total_counts = abund_df4.loc[:, matched_cols].sum().sum()
            print("{} bp, {} OTUs, {:.2%} counts".format(i, j, total_counts/all_counts))

    # clean up length distribution
    abund_df5 = abund_df4.loc[:, ~abund_df4.columns.isin(wrong_lengths)]
    abund_df5.to_csv("../data/otu_data_pca/clustered_otu_abundances.tsv", sep="\t")
else:
    # do PCA
    abund_df5 = pd.read_csv("../data/otu_data_pca/clustered_otu_abundances.tsv", sep="\t", index_col=0)
    # write out taxa
    heads = sorted(list(abund_df5.columns))
    tails = [OTU_name2seq[i] for i in heads]
    with open('../data/otu_data_pca/clustered_otus.fa', "w") as wofh:
        wofh.write("".join([">{}\n{}\n".format(i, j) for i, j in zip(heads, tails)]))
        
    # load sample sheet and latest seq run data
    sample_sheet = load_sample_sheet()
    control_indexes = [270, 603, 607, 608, 609, 610, 668, 674] + [sample_sheet.index[55]] + [671, 672]
    control_sids = sample_sheet.loc[control_indexes, 'SampleID'].tolist()
    new_ssu_fn = '../data/LibraryPrepProtocols/ThisRunToSampleProcessing.tsv'
    n_ssu_df = pd.read_csv(new_ssu_fn, sep="\t", index_col=0)

    # Fix weird date
    sample_sheet.loc[sample_sheet['DateMMDDYY'] == 'Mix9', 'DateMMDDYY'] = '100516'

    # make weird samples (not mine or controls, some not in sample sheet) their own group 
    weird_samps = set(abund_df3.index) - set(sample_sheet.SampleID.unique())

    # make a list of "stations" corresponding to otutable 
    stat_unencoded = []
    date_unencoded = []
    seq_run_unencoded = []
    for s_ix in abund_df5.index:
        if s_ix in weird_samps:
            this_stat, this_dat, this_run = ["LAB"], ["LAB"], [n_ssu_df.loc[s_ix, 'Run']]
        elif s_ix in sample_sheet.SampleID.values:
            this_stat = sample_sheet.loc[sample_sheet.SampleID == s_ix, 'StationName'].values
            this_dat = sample_sheet.loc[sample_sheet.SampleID == s_ix, 'DateMMDDYY'].values
            this_run = sample_sheet.loc[sample_sheet.SampleID == s_ix, 'sequencing ID'].values
            assert (this_stat.shape[0] == 1) and (this_dat.shape[0] == 1) and (this_run.shape[0] == 1)
            assert (type(this_run[0]) != type(0.1)) and (type(this_run[0]) != type(n_ssu_df.loc['ANF8_C2', 'Date']))
        else:
            raise ValueError("Illegal index discovered")
        stat_unencoded.append(this_stat[0])
        date_unencoded.append(this_dat[0])
        seq_run_unencoded.append(this_run[0])

    # station name to encoding integer for categorical
    stat_encoding = {stat: idx for idx, stat in enumerate(sample_sheet.StationName.unique())}
    stat_encoded = [stat_encoding[i] for i in stat_unencoded]
    run_encoding = {sqr: idx for idx, sqr in enumerate(set(seq_run_unencoded))}
    run_rev_enc = {idx: sqr for idx, sqr in enumerate(set(seq_run_unencoded))}
    run_encoded = [run_encoding[i] for i in seq_run_unencoded]
    date_encoding = {date_: idx for idx, date_ in enumerate(set(date_unencoded))}
    date_encoded = [date_encoding[i] for i in date_unencoded]

    metadata = np.array([date_encoded, run_encoded, stat_encoded]).T
    meta_df = pd.DataFrame(index=abund_df3.index, columns=['date', 'run', 'station'], data=metadata)

    #negative_idxs = feature_loading[feature_loading.iloc[:, 0] < -.35].index
    #normal_idxs = feature_loading[~feature_loading.index.isin(negative_idxs)].index
    #normal_subset = np.random.choice(normal_idxs, size=negative_idxs.shape)

    control_libs = list(abund_df3.index[abund_df3.index.str.contains("_SJC")])
    control_libs += ['178A_WaterBathControlA', '178B_WaterBathControlB']
    control_libs += list(abund_df3.index[abund_df3.index.str.startswith("ANF")])
    control_libs += list(abund_df3.index[abund_df3.index.str.startswith("BF")])
    control_libs += list(abund_df3.index[abund_df3.index.str.contains("Blank")])
    control_libs += list(abund_df3.index[abund_df3.index.str.contains("Zymo")])
    control_libs += list(abund_df3.index[abund_df3.index.str.contains("Mix9")])
    control_libs += list(abund_df3.index[abund_df3.index.str.contains("EMPTY")])
    control_libs += list(abund_df3.index[abund_df3.index.str.contains("CDSBBR")])

    otus_to_strip_c, otus_to_strip_nc = set(), set()
    for c in control_libs:
        otus_to_strip_nc.update(abund_df3.columns[abund_df3.loc[c, :] > 0])
        otus_to_strip_c.update(abund_df5.columns[abund_df5.loc[c, :] > 0])
        print("{}/{} are to be removed".format(len(otus_to_strip_nc), len(otus_to_strip_c)))

    # remove OTUS in 50 unecessary samples 
    abund_df_c = abund_df5.loc[:, ~abund_df5.columns.isin(otus_to_strip_c)]
    abund_df_c2 = abund_df_c[abund_df_c.sum(1)!=0]
    abund_df_nc = abund_df3.loc[:, ~abund_df3.columns.isin(otus_to_strip_nc)]
    abund_df_nc2 = abund_df_nc[abund_df_nc.sum(1)!=0]

    rename_cols = {i - 1: 'PC' + str(i) for i in range(1, 3)}
    rclr_mat = rclr().fit_transform(abund_df_c2.values)
    rclr_mat2 = rclr().fit_transform(abund_df_nc2.values)
    opt=OptSpace(rank=2).fit(rclr_mat)
    # PCA of new mat
    feature_loading = pd.DataFrame(opt.feature_weights, index=abund_df_c2.columns).rename(columns=rename_cols)
    feature_loading.sort_values('PC1', inplace=True, ascending=True)
    sample_loading = pd.DataFrame(opt.sample_weights, index=abund_df_c2.index).rename(columns=rename_cols)
    # PCA of old mat
    opt2=OptSpace(rank=2).fit(rclr_mat2)
    fl_clean = pd.DataFrame(opt2.feature_weights, index=abund_df_nc2.columns).rename(columns=rename_cols)
    fl_clean.sort_values('PC1', inplace=True, ascending=True)
    sl_clean = pd.DataFrame(opt2.sample_weights, index=abund_df_nc2.index).rename(columns=rename_cols)
    sub_meta_c_st = meta_df.loc[abund_df_c2.index, 'station'].tolist()
    sub_meta_c_date = meta_df.loc[abund_df_c2.index, 'date'].tolist()
    sub_meta_c = meta_df.loc[abund_df_c2.index, 'run'].tolist()
    sub_meta_nc = meta_df.loc[abund_df_nc2.index, 'run'].tolist()

    plt.clf()
    fig, ax_arr = plt.subplots(nrows=2, ncols=2, figsize=(16,16), dpi=250)
    ax0, ax1, ax2, ax3 = ax_arr[0,0], ax_arr[0,1], ax_arr[1,0], ax_arr[1,1]
    ax0.set_title('Before clustering')
    im = ax0.scatter(sl_clean.iloc[:, 0], sl_clean.iloc[:, 1], c=sub_meta_nc, edgecolor='none', alpha=0.5, cmap=plt.cm.get_cmap('Spectral', len(set(sub_meta_nc))))
    fig.colorbar(im, ax=ax0);
    ax1.set_title('After clustering')
    im1 = ax1.scatter(sample_loading.iloc[:, 0], sample_loading.iloc[:, 1], c=sub_meta_c, edgecolor='none', alpha=0.5, cmap=plt.cm.get_cmap('Spectral', len(set(sub_meta_c))))
    fig.colorbar(im1, ax=ax1);
    ax2.set_title('Color by station')
    im2 = ax2.scatter(sample_loading.iloc[:, 0], sample_loading.iloc[:, 1], c=sub_meta_c_st, edgecolor='none', alpha=0.5, cmap=plt.cm.get_cmap('Spectral', len(set(sub_meta_c_st))))
    fig.colorbar(im2, ax=ax2);
    ax3.set_title('Color by date')
    im3 = ax3.scatter(sample_loading.iloc[:, 0], sample_loading.iloc[:, 1], c=sub_meta_c_date, edgecolor='none', alpha=0.5, cmap=plt.cm.get_cmap('Spectral', len(set(sub_meta_c_date))))
    fig.colorbar(im3, ax=ax3);
    for ax in [ax0, ax1, ax2, ax3]:
        ax.set_xlabel(sl_clean.columns[0])
        ax.set_ylabel(sl_clean.columns[1])
    
    plt.savefig("../data/otu_data_pca/PCA_before_after_clustering_plus_other_colors.png", bbox_inches='tight')

    runs_to_keep = sorted([code for code in run_rev_enc.keys() if len(meta_df.loc[meta_df.run == code, 'station'].unique()) > 2])
    zee_mats = [abund_df_c2, abund_df_nc2]
    zee_cols = [meta_df.loc[abund_df_c2.index, :], meta_df.loc[abund_df_nc2.index, :]]
    zee_results = []
    for adf, mdf in zip(zee_mats, zee_cols):
        shared_otus = pd.DataFrame(index=runs_to_keep, columns=runs_to_keep)
        shared_group = set(adf.columns)
        for run_grp1 in shared_otus.columns:
            for run_grp2 in shared_otus.index:
                otus_in_1 = adf.loc[mdf.run == run_grp1, :].sum()
                otus_in_2 = adf.loc[mdf.run == run_grp2, :].sum()
                foldchange = otus_in_1 / (otus_in_2)
                num_shared = ((foldchange < 100) & (foldchange > 0.01)).sum()
                # optional shared across all groups
                shared_here = foldchange.index[(foldchange < 100) & (foldchange > 0.01)]
                shared_group = shared_group.intersection(set(shared_here))
                shared_otus.loc[run_grp2, run_grp1] = num_shared
        zee_results.append(shared_otus.copy())
        print("{} otus shared across all runs".format(len(shared_group)))
        del shared_otus

    plt.clf()
    f, _axes_ = plt.subplots(nrows=1, ncols=2, figsize=(18, 12))
    sns.heatmap(zee_results[0], annot=True, fmt="d", linewidths=.5, ax=_axes_[0], cbar=False)
    sns.heatmap(zee_results[1], annot=True, fmt="d", linewidths=.5, ax=_axes_[1], cbar=False)
    plt.savefig("../data/otu_data_pca/OTU_Overlap_groups_postclustering.png", bbox_inches='tight')
    plt.clf()
    plt.close()















    # make sequencing run and date categories (constrained_layout=True)
    mat_pairs = [(sample_loading, sl_clean), (sample_loading, sl_nn)]
    col_pairs = [(run_encoded, sub_run_encoded), (run_encoded, run_encoded)]
    label_pairs = ['knowncontrols', 'negsampleloadings']

for mat_pair, col_pair, label_ in zip(mat_pairs, col_pairs, label_pairs):
    plt.clf()
    fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, sharex=True, figsize=(16,8), dpi=250)
    im = ax0.scatter(mat_pair[0].iloc[:, 0], mat_pair[0].iloc[:, 1],
                     c=col_pair[0], edgecolor='none', alpha=0.5,
                     cmap=plt.cm.get_cmap('Spectral', len(set(col_pair[0]))))
    fig.colorbar(im, ax=ax0);
    im1 = ax1.scatter(mat_pair[1].iloc[:, 0], mat_pair[1].iloc[:, 1],
                      c=col_pair[1], edgecolor='none', alpha=0.5,
                      cmap=plt.cm.get_cmap('Spectral', len(set(col_pair[1]))))
    fig.colorbar(im1, ax=ax1);
    for ax in [ax0, ax1]:
        ax.set_xlabel(mat_pair[0].columns[0])
        ax.set_ylabel(mat_pair[0].columns[1])
    
    plt.savefig("../data/otu_data_pca/PCA_before_after_{}_removed.png".format(label_), bbox_inches='tight')
    
plt.clf()
plt.close()

# just remove negative sample loadings (all samples retained)
abund_df7 = abund_df3.loc[:, ~abund_df3.columns.isin(negative_idxs)]
abund_df8 = abund_df7[abund_df7.sum(1)!=0]
rclr_mat3 = rclr().fit_transform(abund_df8.values)
opt3=OptSpace(rank=2).fit(rclr_mat3)
fl_nn = pd.DataFrame(opt3.feature_weights, index=abund_df8.columns)
fl_nn = fl_nn.rename(columns=rename_cols)
fl_nn.sort_values('PC1', inplace=True, ascending=True)
sl_nn = pd.DataFrame(opt3.sample_weights, index=abund_df8.index)
sl_nn = sl_nn.rename(columns=rename_cols)

sum_list = [abund_df3.sum(1), abund_df3.loc[:, negative_idxs].sum(1), 
            abund_df3.loc[:, normal_idxs].sum(1)]
pr_ab_list = [(abund_df3 > 0.0).sum(1), (abund_df3.loc[:, negative_idxs] > 0.0).sum(1),
              (abund_df3.loc[:, normal_idxs] > 0.0).sum(1)]

df = pd.concat(sum_list+pr_ab_list, axis=1, verify_integrity=True)
df.columns = ['A. All Counts', "B. Negative Loading Counts", "C Positive Loading Counts"] + \
              ['D All Observed OTUs', "E. Negative Loading Observed OTUs", "F. Postiive Loading Observed OTUs"]
# this will be the sequencing run encoding
df['Run'] = meta_df.run
# this will plot it 
fig1, ax1 = plt.subplots(nrows=2, ncols=3, sharey='row', sharex='col', figsize=(16,8), dpi=250, constrained_layout=True)
ax2 = df.boxplot(by='Run', ax=ax1, return_type='axes')
for n, l in zip(ax2.index, [300000]*3 + [2000]*3):
    ax2[n].set_ylim([0,l])

left = .1
bottom = .03
fig1.text(left, bottom, ", ".join(["{}: {}".format(i, j) for i, j in run_rev_enc.items()]))
plt.savefig('../data/otu_data_pca/otu_counts_by_type_and_run.png')
plt.clf()
plt.close()

abund_df9 = abund_df5.loc[:, ~abund_df5.columns.isin(negative_idxs)]
abund_df10 = abund_df9[abund_df9.sum(1)!=0]










# even sampling 
np.random.seed(42)
from collections import Counter
def even_sampling(composition, depth):
    pile = [j for i in composition.index for j in [i]*composition[i]]
    assert len(pile) == composition.sum()
    sample = np.random.choice(pile, size=(depth,))
    counted = Counter(sample)
    print(len(counted.keys()), "remain of", (composition > 0).sum(), "from {}".format(composition.name))
    new_comp = pd.Series(index=composition.index, data=np.zeros(composition.shape))
    for i, j in counted.items():
        new_comp[i]= j
    return new_comp


#rare_df = pd.DataFrame(index=abund_df5.index, columns=abund_df5.columns)
#
#for samp in abund_df5.index:
#    samp_srs = abund_df5.loc[samp, :].copy()
#    rare_df.loc[samp, :] = even_sampling(samp_srs, 5000)

# how many OTUs are observed in and outside my last two runs

#rare_df2 = rare_df.loc[:, rare_df.sum() > 0]

#shared_otus = pd.DataFrame(index=sub_meta.run.unique(), columns=sub_meta.run.unique())


#shared_otus = shared_otus.rename(columns=run_rev_enc)
#shared_otus = shared_otus.rename(index=run_rev_enc)

# they do cluster into same large groups but the cluster's correlation across samples isn't very high 
heads = sorted(list(negative_idxs))
tails = [OTU_name2seq[i] for i in heads]
with open('../data/otu_data_pca/weird_otus.fa', "w") as wofh:
    wofh.write("".join([">{}\n{}\n".format(i, j) for i, j in zip(heads, tails)]))
sLens, n_sLens = np.unique([len(i) for i in tails], return_counts=True)
print("Sequence Length Counts in weird otus are:\n"+"".join(["\t{}:{}\n".format(i, j) for i, j in zip(sLens, n_sLens)]))



# taxanomic assignment seems to work better on comtaminants ? so they prob are not artefacts
for class_level in taxa_df.columns:
    weird_null = taxa_df.loc[negative_idxs, class_level].isnull().sum() / taxa_df.loc[negative_idxs, class_level].shape[0]
    normal_null = taxa_df.loc[normal_idxs, class_level].isnull().sum() / taxa_df.loc[normal_idxs, class_level].shape[0]
    print("# of null {} classes in contaminations is {:.2%} and in normals is {:.2%}".format(class_level, weird_null, normal_null))

# how many are in blank
for c in ["93_PBS_Blank_Control", '92_Mix93_Control_R2', '92_Mix93_Control_R1', '96_ZymoControl_R1']:
    presence_absence = (abund_df3.loc[c, negative_idxs] > 0).sum() /  abund_df3.loc[c, :].shape[0]
    proportion = abund_df3.loc[c, negative_idxs].sum() / abund_df3.loc[c, :].sum()
    print("Weird indexes make up {:.2%} of the otus and {:.2%} of the counts in {}".format(presence_absence, proportion, c))


# look at bulk values as well such as observed otu #
# one histogram of all libraries
# one histogram of libraries seperated by run 

# this will be replaced by total OTUs, weird OTUS, and random sample of non-weird OTUs


# calculate the average per sample of these otus
outlier_sum_by_sample = abund_df3.loc[:, negative_idxs].sum(axis=1)

# take a rnadom sample of other otus and do a comparison of average values across libraries
# also look at the level of taxanomic assignment among both groups

# look at other cleaner libraries for an indiciatino of what to keep!
# also look at taxa 


# cluster seqs: nothing
# rarify 

# find what samples are closest to the weird controls
# steal tsne plot
rclr_mat4 = rclr().fit_transform(abund_df10.values + 0.00001)
model = PCA(n_components=rclr_mat4.shape[0], svd_solver='randomized', random_state=42)
model.fit(rclr_mat4)
projected = model.fit_transform(rclr_mat4)
pc_df = pd.DataFrame(index=abund_df10.index, data=projected)
sns.set(font_scale=1.5)
plt.clf()
corrplot = sns.clustermap(pc_df, method='ward', figsize=(120,120), cmap="coolwarm", robust=True)
corrplot.savefig('../data/otu_data_pca/pca_of_rclr_ward_clustermap.png', dpi=250)
plt.clf(); plt.close();

bc_dm = beta_diversity("braycurtis", abund_df3.values, list(abund_df3.index))

bc_df = pd.DataFrame(index=abund_df3.index, columns=abund_df3.index, data=bc_dm._data)

cs_df = bc_df.loc[control_sids, :]
for x in cs_df.index:
    print(x)
    print(cs_df.loc[x, :].sort_values().head())

bc_dm = beta_diversity("braycurtis", abund_df5.values, list(abund_df5.index))

bc_pc = pcoa(bc_dm)
fig = bc_pc.plot(sub_meta, 'run', axis_labels=('PC 1', 'PC 2', 'PC 3'), title='Samples colored by run', cmap='jet', s=50)
plt.savefig("../data/otu_data_pca/braycurtis_pcoa_byrun_after.png", bbox_inches='tight')
plt.clf()

# Feature Loadings
feature_loading = pd.DataFrame(opt.feature_weights, index=abund_df3.columns)
feature_loading = feature_loading.rename(columns=rename_cols)
feature_loading.sort_values('PC1', inplace=True, ascending=True)









model = PCA(n_components=8, svd_solver='randomized', random_state=42)
model.fit(df4_scaled)
projected = model.fit_transform(df4_scaled)
mevr = model.explained_variance_ratio_

fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, sharex=True, figsize=(16,8), dpi=250)
color_types = [('depth', c_cat_encoded2), ('latitude', c_cat_encoded)]
fig_n_comps = [(10, 0, 1), (11, 2, 1), (12, 0, 2)]



plt.clf()
# decomp methods


plt.clf()
for ct in color_types:
    for fnc in fig_n_comps:
        cmp1, cmp2 = fnc[1]+1, fnc[2]+1
        plt.figure(fnc[0], figsize=(8,8), dpi=250)
        plt.scatter(projected[:, fnc[1]], projected[:, fnc[2]],
                    c=ct[1], edgecolor='none', alpha=0.5,
                    cmap=plt.cm.get_cmap('Spectral', len(set(ct[1]))))
        plt.xlabel("component {} ({:.1%} variance)".format(cmp1, mevr[fnc[1]]))
        plt.ylabel("component {} ({:.1%} variance)".format(cmp2, mevr[fnc[2]]))
        plt.colorbar();
        plt.savefig("../data/wq_pca_{}_{}_c{}_c{}.png".format(fnc[0], ct[0], cmp1, cmp2), bbox_inches='tight')
        plt.clf()



# cluster map
#rows = compositions and columns = components
"""

