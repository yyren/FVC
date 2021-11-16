import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams
import matplotlib.pyplot as plt
from config import PLOT_PATH, TAB_PATH, DATABASE_PATH, normal_dict, dark_dict
from config import order_caller, order_model, order_strategy, order_mutation
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
alias_dict = {**order_caller,**order_model,**order_strategy,**order_mutation}

## read metrics table of fnvc features
metr_df = pd.read_pickle(os.path.join(DATABASE_PATH, 'reb_imb_metr_df.pkl'))
select_row = (metr_df['patient']!='HG007') & (metr_df['depth']=='30X')
metr_df = metr_df.loc[select_row]
metr_df.insert(metr_df.columns.to_list().index('TNR')+1, 'TPR', 1 - metr_df['FNR'])
metr_df[['Model','Caller','Mutation','Balance Strategy']] = metr_df[['model','software','mutation','balance_strategy']].applymap(lambda x: alias_dict[x])
metr_df_pp = metr_df.groupby(by=['Caller','Model','Balance Strategy','Mutation','patient']).mean()
metr_df_mn = metr_df.groupby(by=['Caller','Model','Balance Strategy','Mutation']).mean()
metr_df_sd = metr_df.groupby(by=['Caller','Model','Balance Strategy','Mutation']).std()
# dfs_to_sheet([metr_df, metr_df_pp, metr_df_mn, metr_df_sd], ['raw', 'per patient', 'mean', 'std'], filename='reb_original_vs_bal.xlsx', root_path=TAB_PATH)

## plot XGBoost
model = 'XGBoost'

target_col = ['MCC','patient','Caller','Model','Balance Strategy','Mutation']
bar_plot_df = metr_df[target_col].loc[metr_df['Model']==model]
g = sns.FacetGrid(data=bar_plot_df, row='Mutation', aspect=2, height=4)
g.map(
    sns.barplot, 'Caller', 'MCC', 'Balance Strategy',
    order=['GATK', 'Mutect2', 'Varscan2', 'DeepVariant'],
    hue_order=['original', 'SMOTE', 'oversample', 'undersample', 'NearMiss-1'],
    palette=get_pal(normal_dict, ['black','red','orange','purple','blue']),
    errwidth=2, capsize=0.05)
g.add_legend()
g.map(
    sns.stripplot, 'Caller', 'MCC', 'Balance Strategy',
    order=['GATK', 'Mutect2', 'Varscan2', 'DeepVariant'],
    hue_order=['original', 'SMOTE', 'oversample', 'undersample', 'NearMiss-1'],
    palette=get_pal(dark_dict, ['black','red','orange','purple','blue']),
    edgecolor='white', linewidth=0.5, alpha=0.7,
    dodge=True)
# sns.catplot(
#     data=bar_plot_df, x='Caller', y='MCC', hue='Balance Strategy', col='Mutation',
#     order=['GATK', 'Mutect2', 'Varscan2', 'DeepVariant'],
#     hue_order=['original', 'SMOTE', 'oversample', 'undersample', 'NearMiss-1'],
#     palette=get_pal(normal_dict, ['black','red','orange','purple','blue']),
#     kind='bar', errwidth=2, capsize=0.1)
plt.savefig(os.path.join(PLOT_PATH, 'reb_original_vs_bal_{}.pdf'.format(model)))
# plt.show()

### compare two two feature settings and generate p value
y_metrname = 'MCC'
test_tab_df_list = []
for mutation in set(bar_plot_df['Mutation']):
    for caller in set(bar_plot_df['Caller']):
        for strategy in set(bar_plot_df['Balance Strategy']).difference(['original']):
            x = bar_plot_df.loc[
                (bar_plot_df['Mutation']==mutation)&
                (bar_plot_df['Caller']==caller)&
                (bar_plot_df['Balance Strategy']=='original')][y_metrname].values
            y = bar_plot_df.loc[
                (bar_plot_df['Mutation']==mutation)&
                (bar_plot_df['Caller']==caller)&
                (bar_plot_df['Balance Strategy']==strategy)][y_metrname].values
            # print(x.shape, y.shape)
            test_tab_df_list.append([mutation, caller, strategy, *noparam_wrm(x,y)])
test_tab_df = pd.DataFrame(np.array(test_tab_df_list), columns=['Mutation','Caller', 'Balance Strategy', 'original Ave','original SD','strategy Ave','strategy SD','statistic','pvalue'])
# test_tab_df.to_excel(os.path.join(TAB_PATH, 'reb_pvalue_original_vs_bal_{}.xlsx'.format(model)))
