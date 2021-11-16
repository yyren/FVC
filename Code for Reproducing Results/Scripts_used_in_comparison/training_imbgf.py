from config import RESULT_PATH, DATABASE_PATH, DATABASE_TWO_PATH
from utils import load_splitdata, random_sample, save_to_pickle, \
    load_from_pickle, check_loadsave, training_series, evaluating_series, \
        load_30xonly, load_50xonly
from utils import mutypes, pcodes
from os import path, makedirs
from sys import argv
import numpy as np
'''
seed = 0
mt = 'indel'
pc = 'HG001'
software = 'gatk'
gfcode = 'tc'
'''

seed = int(argv[1])
mt = argv[2]
pc = argv[3]
software = argv[4]
gfcode = argv[5]

if gfcode == 'tc':
    root_path = DATABASE_PATH
elif gfcode == 'gf':
    root_path = DATABASE_TWO_PATH

prespath = path.join(RESULT_PATH, software, mt, pc)

X_train, y_train, X_test30x, y_test30x = check_loadsave(
    path.join(prespath, 'data_suitcase_{}_30x.pkl'.format(gfcode)), 
    load_30xonly, {'pcode': pc, 'mutype': mt, 'software': software, 'root_path':root_path})

X_train50x, y_train50x, X_test50x, y_test50x = check_loadsave(
    path.join(prespath, 'data_suitcase_{}_50x.pkl'.format(gfcode)), 
    load_50xonly, {'pcode': pc, 'mutype': mt, 'software': software, 'load_train': False, 'root_path':root_path})

Xs_test = [
    X_test30x, 
    X_test50x
    ]
ys_test = [
    y_test30x, 
    y_test50x
    ]

print('training on 30x data')
clfkit_list, normalizer = check_loadsave(
    path.join(prespath, 'clf', 'clfkits_{}_{}.pkl'.format(seed, gfcode)),
    training_series, {
        'X_train': X_train,
        'y_train': y_train,
        'model_list': ['xgbdef'],
        'model_params': [{}],
        'seed': seed})

# print('testing on 30x')
print('testing on 30x and 50x')
ys_df, metrs_df = check_loadsave(
    path.join(prespath, 'metr', 'metrs_{}_{}.pkl'.format(seed, gfcode)), 
    evaluating_series, {
        'Xs_test': Xs_test, 
        'ys_test': ys_test,
        'clfkit_list': clfkit_list, 
        'normalizer': normalizer})

