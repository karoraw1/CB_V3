import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def read_cmsearch_file(fn):
    with open(fn, 'r') as chode:
        content = [i for i in chode.read().split("\n") if i!=""]
    colnames = ['#target_name', 'accession', 'query_name', 'accession', 'mdl', 'mdl_from', 'mdl_to',
                'seq_from', 'seq_to', 'strand', 'trunc', 'pass', 'gc', 'bias', 'score', 'E-value',
                'inc', 'description_of_target']
    data = np.array([i.split() for i in content[2:-10]])
    df = pd.DataFrame(index=data[:, 0], columns=colnames[1:], data=data[:, 1:])
    int_cols = ['seq_to', 'seq_from', 'mdl_to', 'mdl_from']
    df['E-value'] = np.log(df['E-value'].astype(np.float64))*-1
    df['score'] = df['score'].astype(np.float64)
    temp = df['E-value']
    perf_matches = (temp == np.inf).sum()
    temp[temp == np.inf] = 0.
    print(" {} perfect matches detected. Max value set at {}.".format(perf_matches,temp.max()))
    truncate_val = temp.max()+10
    df.loc[df['E-value'] == np.inf, 'E-value'] = truncate_val
    for col in int_cols:
        df[col] = df[col].astype(int)
    return df

data_path = sys.argv[-1]
report_file = sys.argv[-2]
out_dir = sys.argv[-3]

the_fn = os.path.join(data_path, report_file);
assert os.path.exists(the_fn)
cm_df = read_cmsearch_file(the_fn)

cm_df["Group"] = (cm_df['E-value']*0.0).apply(int)
cm_df.loc[cm_df.index.str.startswith("OTU"), "Group"] = 1

s1 = cm_df[cm_df.Group == 0]
s2 = cm_df[cm_df.Group == 1]
diff_figs = [('Ref_Seqs', s1, 100), ('OTUs', s2, 500)]
plottable_columns = ['seq_from', 'seq_to', 'mdl_from', 'mdl_to', 'E-value', 'score']
axis_tuples = [(0,0), (0,1), (1,0), (1,1), (2,0), (2,1)]

for _n_, _df_, _bins_ in diff_figs:
    plt.clf()
    fig1, ax1 = plt.subplots(nrows=3, ncols=2, figsize=(16,8), dpi=200, constrained_layout=True);
    for (r_a, c_a), col_x in zip(axis_tuples, plottable_columns):
        df4 = _df_[col_x]
        df4.plot.hist(bins=_bins_, ax=ax1[r_a, c_a])
        ax1[r_a, c_a].set_title(col_x)
    plt.savefig(os.path.join(out_dir, 'cmsearch_data_{}.png'.format(_n_)))

plt.clf();plt.close();

fig2, ax2 = plt.subplots(nrows=1, ncols=1, figsize=(8,8), dpi=200)
df4.plot.hist(bins=500, title='Exp_E-vals', ax=ax2)
ax2.set_ylim([0, 200])
ax2.set_xlim([125, 150])
plt.savefig(os.path.join(out_dir, 'cmsearch_Exp_Evalues.png'))

e_lim = 132.0
to_drop = s2[s2['E-value'] <= e_lim].index
rev_comps = s2[(s2.strand == '-') and (s2['E-value'] > e_lim)].index
print("{} outliers detected, which make up {:.2%} of OTUs".format(to_drop.shape[0], to_drop.shape[0]/cm_df.shape[0]))
print("{} are reverse complimented".format(rev_comps.shape[0]))

for i_list, ls_name in zip([to_drop, rev_comps], ['poor_aligners.txt', 'reverse_strand_aligners.txt']):
    to_write = "\n".join(list(i_list))
    with open(os.path.join(out_dir, ls_name), 'w') as pa_fh:
       pa_fh.write(to_write)
