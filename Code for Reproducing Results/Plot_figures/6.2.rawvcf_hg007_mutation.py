import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams
import matplotlib.pyplot as plt
from config import PLOT_PATH, TAB_PATH, DATABASE_PATH, normal_dict, dark_dict
from config import order_caller, order_model, order_strategy, order_mutation, order_filter
from utils import dfs_to_sheet, noparam_wrm, get_pal
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
alias_dict = {**order_caller,**order_model,**order_strategy,**order_mutation,**order_filter}

## metr table
metr_df = pd.read_pickle(os.path.join(DATABASE_PATH, 'rawvcf_filter_metr_df.pkl'))
select_row = metr_df['patient']=='HG007'
metr_df = metr_df.loc[select_row]
metr_df.insert(metr_df.columns.to_list().index('TNR')+1, 'TPR', 1 - metr_df['FNR'])
metr_df[['Caller','Mutation','Filter']] = metr_df[['software','mutation','filter']].applymap(lambda x: alias_dict[x])
metr_df_pp = metr_df.groupby(by=['depth','Caller','Mutation','Filter','patient']).mean()
metr_df_mn = metr_df.groupby(by=['depth','Caller','Mutation','Filter']).mean()
metr_df_sd = metr_df.groupby(by=['depth','Caller','Mutation','Filter']).std()
# dfs_to_sheet([metr_df, metr_df_pp, metr_df_mn, metr_df_sd], ['raw', 'per patient', 'mean', 'std'], filename='rawvcf_filter_hg007_mutation.xlsx', root_path=TAB_PATH)


## roc plot
roc_df = pd.read_pickle(os.path.join(DATABASE_PATH, 'rawvcf_filter_roc_df.pkl'))
select_row = roc_df['patient']=='HG007'
roc_df = roc_df.loc[select_row]
roc_df[['Caller','Mutation','Filter']] = roc_df[['software','mutation','filter']].applymap(lambda x: alias_dict[x])

for depth in ['30X', '50X']:
    g = sns.relplot(data=metr_df.loc[metr_df['depth']==depth], row='Mutation', col='Caller', kind='scatter', height=5, aspect=1,
        x='FPR', y='TPR', style='Filter', markers=['x']*len(set(metr_df['Filter'])),
        hue='Filter', hue_order=['FNVC', 'VEF', 'Frequency', 'Garfield', 'VQSR', 'Hard Filter'], 
        palette=[snsdeep[i] for i in [0,9,1,3,4,8]])
    g.set(xlim=(-0.1,1.1), ylim=(-0.1,1.1))
    plt.savefig(os.path.join(PLOT_PATH, 'rawvcf_filter_hg007_mutation_roc_scatter_{}.pdf'.format(depth)))
    # plt.show()

    g = sns.relplot(data=roc_df.loc[roc_df['depth']==depth], row='Mutation', col='Caller', kind='line', height=5, aspect=1,
        x='fpr', y='tpr', 
        hue='Filter', hue_order=['FNVC', 'VEF', 'Frequency', 'Garfield', 'VQSR', 'Hard Filter'], 
        palette=[snsdeep[i] for i in [0,9,1,3,4,8]])
    g.set(xlim=(-0.1,1.1), ylim=(-0.1,1.1))
    plt.savefig(os.path.join(PLOT_PATH, 'rawvcf_filter_hg007_mutation_roc_line_{}.pdf'.format(depth)))
    # plt.show()


## pvalue
test_tab_df_list = []
test_tab_name_list = []
for depth in ['30X', '50X']:
    for y_metrname in ['MCC', 'G1-score']:
        tab_list = []
        for mutation in set(metr_df['Mutation']):
            for caller in set(metr_df['Caller']):
                for filter in set(metr_df['Filter']).difference(['FNVC']):
                    x = metr_df.loc[
                        (metr_df['depth']==depth)&
                        (metr_df['Mutation']==mutation)&
                        (metr_df['Caller']==caller)&
                        (metr_df['Filter']=='FNVC')][y_metrname].values
                    y = metr_df.loc[
                        (metr_df['depth']==depth)&
                        (metr_df['Mutation']==mutation)&
                        (metr_df['Caller']==caller)&
                        (metr_df['Filter']==filter)][y_metrname].values
                    print(x.shape, y.shape)
                    if x.shape[0] == y.shape[0]:
                        tab_list.append([mutation, caller, filter, *noparam_wrm(x,y)])
        test_tab_df_list.append(pd.DataFrame(np.array(tab_list), columns=['Mutation','Caller', 'Filter', 'FNVC Ave','FNVC SD','Filter Ave','Filter SD','statistic','pvalue']))
        test_tab_name_list.append('{}-{}'.format(depth,y_metrname))
# dfs_to_sheet(test_tab_df_list, test_tab_name_list, 'rawvcf_filter_hg007_mutation_pvalue.xlsx', TAB_PATH)

