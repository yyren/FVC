import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams
import matplotlib.pyplot as plt
from config import PLOT_PATH, TAB_PATH, DATABASE_PATH, normal_dict, dark_dict
from config import order_caller, order_model, order_strategy, order_coding, order_filter
from utils import dfs_to_sheet, noparam_wrm, get_pal, nanlog10
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
alias_dict = {**order_caller,**order_model,**order_strategy,**order_coding,**order_filter}

## metr table
metr_df = pd.read_pickle(os.path.join(DATABASE_PATH, 'rawfeat_coding_filter_metr_df.pkl'))
select_row = metr_df['patient']!='HG007'
metr_df = metr_df.loc[select_row]
metr_df.insert(metr_df.columns.to_list().index('TNR')+1, 'TPR', 1 - metr_df['FNR'])
metr_df.insert(metr_df.columns.to_list().index('OFO')+1, 'log OFO', [nanlog10(x) for x in metr_df['OFO'].to_numpy()])
metr_df[['Caller','Coding','Filter']] = metr_df[['software','coding','filter']].applymap(lambda x: alias_dict[x])
metr_df_pp = metr_df.groupby(by=['depth','Caller','Coding','Filter','patient']).mean()
metr_df_mn = metr_df.groupby(by=['depth','Caller','Coding','Filter']).mean()
metr_df_sd = metr_df.groupby(by=['depth','Caller','Coding','Filter']).std()
# dfs_to_sheet([metr_df, metr_df_pp, metr_df_mn, metr_df_sd], ['raw', 'per patient', 'mean', 'std'], filename='rawfeat_filter_coding.xlsx', root_path=TAB_PATH)


## bar plot
for depth in ['30X', '50X']:
    for y_metrname in ['MCC', 'log OFO']:
        bar_plot_df = metr_df.loc[metr_df['depth']==depth]
        g = sns.FacetGrid(data=bar_plot_df, row='Coding', aspect=1.6, height=4)
        g.map(
            sns.barplot, 'Caller', y_metrname, 'Filter',
            order=['GATK', 'Varscan2', 'Mutect2', 'DeepVariant'],
            hue_order=['FNVC', 'VEF', 'Frequency', 'Garfield', 'Hard Filter', 'VQSR'],
            palette=[snsdeep[i] for i in [0,9,1,3,4,8]],
            errwidth=1, capsize=0.05)
        g.add_legend()
        g.map(
            sns.stripplot, 'Caller', y_metrname, 'Filter',
            order=['GATK', 'Varscan2', 'Mutect2', 'DeepVariant'],
            hue_order=['FNVC', 'VEF', 'Frequency', 'Garfield', 'Hard Filter', 'VQSR'],
            palette=[snsdark[i] for i in [0,9,1,3,4,8]],
            edgecolor='white', linewidth=0.5, alpha=0.7,
            dodge=True)

        # sns.catplot(
        #     data=metr_df.loc[metr_df['depth']==depth], x='Caller', y=y_metrname,
        #     hue='Filter', row='Coding',
        #     order=['GATK', 'Varscan2', 'Mutect2', 'DeepVariant'],
        #     hue_order=['FNVC','Frequency', 'Garfield', 'VQSR', 'Hard Filter', 'VEF'],
        #     palette=[snsdeep[i] for i in [0,1,3,4,8,9]],
        #     kind='bar', aspect=1.6, height=4, errwidth=1, capsize=0.05)
        plt.savefig(os.path.join(PLOT_PATH, 'rawfeat_filter_coding_bar_{}-{}.pdf'.format(depth, y_metrname.replace(' ', '.'))))
        # plt.show()

## pvalue
test_tab_df_list = []
test_tab_name_list = []
for depth in ['30X', '50X']:
    for y_metrname in ['MCC', 'log OFO']:
        tab_list = []
        for coding in set(metr_df['Coding']):
            for caller in set(metr_df['Caller']):
                for filter in set(metr_df['Filter']).difference(['FNVC']):
                    x = metr_df.loc[
                        (metr_df['depth']==depth)&
                        (metr_df['Coding']==coding)&
                        (metr_df['Caller']==caller)&
                        (metr_df['Filter']=='FNVC')][y_metrname].values
                    y = metr_df.loc[
                        (metr_df['depth']==depth)&
                        (metr_df['Coding']==coding)&
                        (metr_df['Caller']==caller)&
                        (metr_df['Filter']==filter)][y_metrname].values
                    print(x.shape, y.shape)
                    if x.shape[0] == y.shape[0]:
                        tab_list.append([coding, caller, filter, *noparam_wrm(x,y)])
        test_tab_df_list.append(pd.DataFrame(np.array(tab_list), columns=['Coding','Caller', 'Filter', 'FNVC Ave','FNVC SD','Filter Ave','Filter SD','statistic','pvalue']))
        test_tab_name_list.append('{}-{}'.format(depth,y_metrname))
# dfs_to_sheet(test_tab_df_list, test_tab_name_list, 'rawfeat_filter_coding_pvalue.xlsx', TAB_PATH)

