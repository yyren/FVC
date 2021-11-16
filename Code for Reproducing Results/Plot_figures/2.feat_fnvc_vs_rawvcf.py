import matplotlib
from seaborn import palettes
matplotlib.use('Agg')
from matplotlib import rcParams
import matplotlib.pyplot as plt
from config import PLOT_PATH, TAB_PATH, DATABASE_PATH
from config import order_caller, order_model, order_strategy, order_mutation
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
fnvc_metr_df = pd.read_pickle(os.path.join(DATABASE_PATH, 'reb_imb_metr_df.pkl'))
select_row = ((fnvc_metr_df['software']=='gatk') | (fnvc_metr_df['software']=='mutect2')) & (fnvc_metr_df['balance_strategy']=='original') & (fnvc_metr_df['patient']!='HG007')
fnvc_metr_df = fnvc_metr_df.loc[select_row]
fnvc_metr_df['Feature Set'] = 'FNVC'

## read metrics table of rawvcf features
## NOTICE because of misunderstanding, the `rawfeat_metr_df` is metrics table of rawvcf features actually
rawfeat_metr_df = pd.read_pickle('../plotdata/rawfeat_imb_metr_df.pkl')
rawfeat_metr_df['Feature Set'] = 'without Annot'

## extract metrics table comparing rawvcf and fnvc features
metr_df = pd.concat((fnvc_metr_df, rawfeat_metr_df))
metr_df.insert(metr_df.columns.to_list().index('TNR')+1, 'TPR', 1 - metr_df['FNR'])
metr_df[['Model','Caller','Mutation','Balance Strategy']] = metr_df[['model','software','mutation','balance_strategy']].applymap(lambda x: alias_dict[x])
metr_df_pp = metr_df.groupby(by=['Caller','Model','Balance Strategy','Mutation','patient','Feature Set']).mean()
metr_df_mn = metr_df.groupby(by=['Caller','Model','Balance Strategy','Mutation','Feature Set']).mean()
metr_df_sd = metr_df.groupby(by=['Caller','Model','Balance Strategy','Mutation','Feature Set']).std()
# dfs_to_sheet([metr_df, metr_df_pp, metr_df_mn, metr_df_sd], ['raw', 'per patient', 'mean', 'std'], filename='feat_fnvc_vs_rawvcf.xlsx', root_path=TAB_PATH)

## plot XGBoost
model = 'XGBoost'

target_col = ['MCC','Feature Set','patient','Caller','Model','Mutation']
bar_plot_df = metr_df[target_col].loc[metr_df['Model']==model]
g = sns.FacetGrid(data=bar_plot_df, row='Mutation', aspect=1, height=4)
g.map(
    sns.barplot, 'Caller', 'MCC', 'Feature Set',
    order=['GATK', 'Mutect2'],
    hue_order=['without Annot', 'FNVC'],
    palette=[snsdeep[i] for i in [1,0]],
    errwidth=2, capsize=0.1)
g.add_legend()
g.map(
    sns.stripplot, 'Caller', 'MCC', 'Feature Set',
    order=['GATK', 'Mutect2'],
    hue_order=['without Annot', 'FNVC'],
    palette=[snsdark[i] for i in [1,0]],
    edgecolor='white', linewidth=1, alpha=0.7,
    dodge=True)
# sns.catplot(
#     data=bar_plot_df, x='Caller', y='MCC', hue='Feature Set', col='Mutation',
#     kind='bar', errwidth=2, capsize=0.1)
plt.savefig(os.path.join(PLOT_PATH, 'feat_fnvc_vs_rawvcf_{}.pdf'.format(model)))
# plt.show()

### compare two two feature settings and generate p value
y_metrname = 'MCC'
test_tab_df_list = []
for mutation in set(bar_plot_df['Mutation']):
    for caller in set(bar_plot_df['Caller']):
        x = bar_plot_df.loc[
            (bar_plot_df['Mutation']==mutation)&
            (bar_plot_df['Caller']==caller)&
            (bar_plot_df['Feature Set']=='FNVC')][y_metrname].values
        y = bar_plot_df.loc[
            (bar_plot_df['Mutation']==mutation)&
            (bar_plot_df['Caller']==caller)&
            (bar_plot_df['Feature Set']=='without Annot')][y_metrname].values
        # print(x.shape, y.shape)
        test_tab_df_list.append([mutation, caller, *noparam_wrm(x,y)])
test_tab_df = pd.DataFrame(np.array(test_tab_df_list), columns=['Mutation','Caller','FNVC Ave','FNVC SD','without Annot Ave','without Annot SD','statistic','pvalue'])
# test_tab_df.to_excel(os.path.join(TAB_PATH, 'feat_pvalue_fnvc_vs_rawvcf_{}.xlsx'.format(model)))


#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------

##! [2021/10/16]LOCO is not necessary for feature comparison.

# # for leave-one-chromosome-out setting 
# alias_dict = {**order_caller,**order_model,**order_strategy,**order_mutation}

# ## read metrics table of fnvc features
# fnvc_metr_df = pd.read_pickle(os.path.join(DATABASE_PATH, 'fnvc_imb_30x_chr_metr_df.pkl'))
# select_row = ((fnvc_metr_df['software']=='gatk') | (fnvc_metr_df['software']=='mutect2')) & (fnvc_metr_df['balance_strategy']=='original')
# fnvc_metr_df = fnvc_metr_df.loc[select_row]
# fnvc_metr_df['Feature Set'] = 'FNVC'

# ## read metrics table of rawvcf features
# rawfeat_metr_df = pd.read_pickle('../plotdata/rawvcf_filter_chr_metr_df.pkl')
# select_row = ((rawfeat_metr_df['software']=='gatk') | (rawfeat_metr_df['software']=='mutect2')) & (rawfeat_metr_df['filter']=='FNVC') & (rawfeat_metr_df['depth']=='30X')
# rawfeat_metr_df = rawfeat_metr_df.loc[select_row]
# rawfeat_metr_df['model'] = 'xgbdef'
# rawfeat_metr_df['chr'] = rawfeat_metr_df['patient']
# rawfeat_metr_df['balance_strategy'] = 'original'
# rawfeat_metr_df['Feature Set'] = 'without Annot'
# rawfeat_metr_df.drop(['filter','patient'], inplace=True, axis=1)

# ## extract metrics table comparing rawvcf and fnvc features
# metr_df = pd.concat((fnvc_metr_df, rawfeat_metr_df))
# metr_df.insert(metr_df.columns.to_list().index('TNR')+1, 'TPR', 1 - metr_df['FNR'])
# metr_df[['Model','Caller','Mutation','Balance Strategy']] = metr_df[['model','software','mutation','balance_strategy']].applymap(lambda x: alias_dict[x])
# metr_df_pp = metr_df.groupby(by=['Caller','Model','Balance Strategy','Mutation','chr','Feature Set']).mean()
# metr_df_mn = metr_df.groupby(by=['Caller','Model','Balance Strategy','Mutation','Feature Set']).mean()
# metr_df_sd = metr_df.groupby(by=['Caller','Model','Balance Strategy','Mutation','Feature Set']).std()
# dfs_to_sheet([metr_df, metr_df_pp, metr_df_mn, metr_df_sd], ['raw', 'per chr', 'mean', 'std'], filename='feat_chr_fnvc_vs_rawvcf.xlsx', root_path=TAB_PATH)

# ## plot XGBoost
# model = 'XGBoost'

# target_col = ['MCC','Feature Set','chr','Caller','Model','Mutation']
# bar_plot_df = metr_df[target_col].loc[metr_df['Model']==model]
# sns.catplot(
#     data=bar_plot_df, x='Caller', y='MCC', hue='Feature Set', col='Mutation',
#     kind='bar', errwidth=2, capsize=0.1)
# plt.savefig(os.path.join(PLOT_PATH, 'feat_chr_fnvc_vs_rawvcf_{}.pdf'.format(model)))
# plt.show()

# ### compare two two feature settings and generate p value
# y_metrname = 'MCC'
# test_tab_df_list = []
# for mutation in set(bar_plot_df['Mutation']):
#     for caller in set(bar_plot_df['Caller']):
#         x = bar_plot_df.loc[
#             (bar_plot_df['Mutation']==mutation)&
#             (bar_plot_df['Caller']==caller)&
#             (bar_plot_df['Feature Set']=='FNVC')][y_metrname].values
#         y = bar_plot_df.loc[
#             (bar_plot_df['Mutation']==mutation)&
#             (bar_plot_df['Caller']==caller)&
#             (bar_plot_df['Feature Set']=='without Annot')][y_metrname].values
#         # print(x.shape, y.shape)
#         test_tab_df_list.append([mutation, caller, *noparam_wrm(x,y)])
# test_tab_df = pd.DataFrame(np.array(test_tab_df_list), columns=['Mutation','Caller','FNVC Ave','FNVC SD','without Annot Ave','without Annot SD','statistic','pvalue'])
# test_tab_df.to_excel(os.path.join(TAB_PATH, 'feat_chr_pvalue_fnvc_vs_rawvcf_{}.xlsx'.format(model)))
