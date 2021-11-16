from copy import copy
import numpy as np
import pandas as pd
from os import path, makedirs
from glob import glob
import pickle
from scipy.stats import ttest_rel



def printx(output):
    print(output, flush = True)


def save_to_pickle(var, filename='temp', root_path='./'):
    '''
    Save variable to pickle file defaultly in `DATABASE_PATH`
    '''
    fn = path.join(root_path, filename)
    if path.exists(path.split(fn)[0]) == False:
        makedirs(path.split(fn)[0])
    printx("saving {} into {}".format(filename, root_path))
    f = open(fn, 'wb')
    pickle.dump(var, f)
    f.close()
    return var


def load_from_pickle(filename='temp', root_path='./'):
    '''
    Load pickle file defaultly from `DATABASE_PATH`
    '''
    printx("loading {} from {}".format(filename, root_path))
    f = open(path.join(root_path, filename), 'rb')
    var = pickle.load(f)
    f.close()
    return var

def dfs_to_sheet(df_list, sheetname_list, filename='temp', root_path='./'):
    '''
    Save a list of csvs into a multi-sheet excel,
        providing list of df and their expected sheet names.
    Return None
    '''
    excel_fn = path.join(root_path, filename)
    excelWriter = pd.ExcelWriter(excel_fn)
    for df,fn in zip(df_list, sheetname_list):
        df.to_excel(excelWriter, sheet_name=fn)
    excelWriter.save()

def get_pal(col_dict, nmidx):
    if isinstance(nmidx, int):
        return list(col_dict.values())[nmidx]
    elif isinstance(nmidx, str):
        return col_dict[nmidx]
    elif isinstance(nmidx, (list, np.ndarray)):
        try:
            return [col_dict[i] for i in nmidx]
        except:
            return [list(col_dict.values())[i] for i in nmidx]


def noparam_wrm(x,y):
    '''
    H0: X is less than Y.
    alternative H: X is greater than Y.
    '''
    return [np.mean(x), np.std(x), np.mean(y), np.std(y), *ttest_rel(x,y,alternative='greater')]

def nanlog10(x):
    if x > 0:
        return np.log10(x)
    else:
        return np.nan