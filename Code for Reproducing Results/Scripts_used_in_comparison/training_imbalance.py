from config import RESULT_PATH, DATABASE_PATH
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
ratio_train = 0
ratio_test = 0
software = 'gatk'

'''

seed = int(argv[1])
mt = argv[2]
pc = argv[3]
ratio_train = float(argv[4])
ratio_test = float(argv[5])
software = argv[6]

prespath = path.join(RESULT_PATH, software, mt, pc)

X_train_ori, y_train_ori, X_test30x_ori, y_test30x_ori = check_loadsave(
    path.join(prespath, 'data_suitcase_30x.pkl'), 
    load_30xonly, {'pcode': pc, 'mutype': mt, 'software': software})

# X_train50x_ori, y_train50x_ori, X_test50x_ori, y_test50x_ori = check_loadsave(
#     path.join(prespath, 'data_suitcase_50x.pkl'), 
#     load_50xonly, {'pcode': pc, 'mutype': mt, 'software': software, 'load_train': False})


seed_train = seed_test = seed
if ratio_train <= 0:
    seed_train = 0
if ratio_test <= 0:
    seed_test = 0

# load mask of under sample
msk_train = check_loadsave(
    path.join(prespath, 'msk', 'msk_train_{}_{:.1f}.pkl'.format(seed_train, ratio_train)), 
    random_sample, {'y': y_train_ori, 'ratio': ratio_train, 'seed': seed_train})
msk_test30x = check_loadsave(
    path.join(prespath, 'msk', 'msk_test30x_{}_{:.1f}.pkl'.format(seed_train, ratio_test)), 
    random_sample, {'y': y_test30x_ori, 'ratio': ratio_test, 'seed': seed_train})
# msk_test50x = check_loadsave(
#     path.join(prespath, 'msk', 'msk_test50x_{}_{:.1f}.pkl'.format(seed_train, ratio_test)), 
#     random_sample, {'y': y_test50x_ori, 'ratio': ratio_test, 'seed': seed_train})

# col_sel = [i for i in range(19) if i not in [16]]
col_sel = np.arange(X_train_ori.shape[1]) # complement for LGB
X_train, y_train = X_train_ori[msk_train][:,col_sel], y_train_ori[msk_train]
Xs_test = [
    X_test30x_ori[msk_test30x][:,col_sel], 
    # X_test50x_ori[msk_test50x][:,col_sel]
    ]
ys_test = [
    y_test30x_ori[msk_test30x], 
    # y_test50x_ori[msk_test50x]
    ]

print('training on 30x data')
clfkit_list, normalizer = check_loadsave(
    path.join(prespath, 'clf', 'clfkits_{}_{:.1f}.pkl'.format(seed, ratio_train)),
    training_series, {
        'X_train': X_train,
        'y_train': y_train,
        'model_list': ['logireg', 'lsvm', 'nn', 'rf', 'xgbdef', 'lgbdef'],
        # 'model_list': ['xgbdef'],
        # 'model_list': ['nn', 'lgbdef'],
        'model_params': [{},{},
            {'hidden_layer_sizes':(50,15)}, 
            {},{},{}],
        # 'model_params': [{}],
        'seed': seed})

print('testing on 30x')
# print('testing on 30x and 50x')
ys_df, metrs_df = check_loadsave(
    path.join(prespath, 'metr', 'metrs_{}_30x_{:.1f}-{:.1f}.pkl'.format(seed, ratio_train, ratio_test)), 
    evaluating_series, {
        'Xs_test': Xs_test, 
        'ys_test': ys_test,
        'clfkit_list': clfkit_list, 
        'normalizer': normalizer})

