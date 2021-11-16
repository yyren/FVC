import matplotlib
from seaborn import palettes
matplotlib.use('Agg')
from matplotlib import rcParams
import matplotlib.pyplot as plt
from config import PLOT_PATH, TAB_PATH, DATABASE_PATH
from config import order_caller, order_model, order_strategy, order_mutation
from config import feature_idx2name, feat_group_id, feat_group_name
from utils import dfs_to_sheet, noparam_wrm
import seaborn as sns
import pandas as pd
import numpy as np
import os

# setting plot parameters
snsdeep = sns.color_palette('deep', 10)
snsdark = sns.color_palette('dark', 10)
snsset2 = sns.color_palette('Set2', 10)

sns.set_theme(
        style='ticks',
        context='talk',
        palette='deep')
params = {
    # 'font.family':'serif',
    # 'font.serif':'Times New Roman',
    'font.family': 'sans-serif',
    'font.sans-serif': 'Arial',
    'font.style':'normal',
    'font.weight':'normal', #or 'blod'
    'pdf.fonttype':42
    # 'font.size':'medium',#or large,small
}
rcParams.update(params)

# for leave-one-patient-out setting 
alias_dict = {**order_caller,**order_model,**order_strategy,**order_mutation}

## read metrics table of fnvc features
fnvc_metr_df = pd.melt(
    pd.read_pickle(os.path.join(DATABASE_PATH, 'fnvc_feature_importance.pkl')),
    id_vars=['caller','mutation','patient'], 
    var_name='feature_id', 
    value_name='Importance Score').dropna(axis=0)
fnvc_metr_df['Feature'] = fnvc_metr_df[['caller', 'feature_id']].apply(lambda x: feature_idx2name(x[0],x[1]), axis=1)
fnvc_metr_df['Feature Group'] = fnvc_metr_df['Feature'].apply(lambda x: feat_group_name[feat_group_id[x]])
fnvc_metr_df[['Caller', 'Mutation']] = fnvc_metr_df[['caller', 'mutation']].applymap(lambda x: alias_dict[x])
# dfs_to_sheet([fnvc_metr_df],['feat imp with HG007'], 'featimp_fnvc.xlsx', TAB_PATH)


g = sns.FacetGrid(
    data=fnvc_metr_df.loc[fnvc_metr_df['patient']!='HG007'], 
    row='Mutation', col='Caller', col_order=['GATK', 'Mutect2', 'Varscan2', 'DeepVariant'],
    aspect=1, height=6)
g.map(
    sns.barplot, 'Importance Score', 'Feature', #'Feature Group',
    order=list(feat_group_id.keys()),
    # hue_order=['Sequence context', 'Sequencing experiment', 'Bioinformatics'],
    palette=[snsdeep[i] for i in [*[5]*4,*[8]*8,*[7]*8]],
    errwidth=2, capsize=0.3)
g.add_legend()
g.map(
    sns.stripplot, 'Importance Score', 'Feature', #'Feature Group',
    order=list(feat_group_id.keys()),
    # hue_order=['Sequence context', 'Sequencing experiment', 'Bioinformatics'],
    palette=[snsdark[i] for i in [*[5]*4,*[8]*8,*[7]*8]],
    edgecolor='white', linewidth=1, alpha=0.7,
    dodge=True)
# sns.catplot(
#     data=bar_plot_df, x='Caller', y='MCC', hue='Feature Set', col='Mutation',
#     kind='bar', errwidth=2, capsize=0.1)
plt.savefig(os.path.join(PLOT_PATH, 'featimp_fnvc_bar.pdf'))
# plt.show()

