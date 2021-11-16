from config import RESULT_PATH, DATABASE_PATH
from utils import load_splitdata, balancer_block, save_to_pickle, \
    load_from_pickle, check_loadsave, training_series, evaluating_series, \
        load_30xonly, load_50xonly, random_sample
from utils import mutypes, pcodes
from os import path, makedirs
from sys import argv
import numpy as np

seed = int(argv[1])
mt = argv[2]
pc = argv[3]
balance_strategy = argv[4]
ratio_test = int(argv[5])
software = argv[6]


prespath = path.join(RESULT_PATH, software, mt, pc)

X_train_ori, y_train_ori, X_test30x_ori, y_test30x_ori = check_loadsave(
    path.join(prespath, 'data_suitcase_30x.pkl'), 
    load_30xonly, {'pcode': pc, 'mutype': mt, 'software': software})

# X_train50x_ori, y_train50x_ori, X_test50x_ori, y_test50x_ori = check_loadsave(
#     path.join(prespath, 'data_suitcase_50x.pkl'), 
#     load_50xonly, {'pcode': pc, 'mutype': mt, 'software': software})

msk_test30x = check_loadsave(
    path.join(prespath, 'msk', 'msk_test30x_{}_{:.1f}.pkl'.format(seed, ratio_test)), 
    random_sample, {'y': y_test30x_ori, 'ratio': ratio_test, 'seed': seed})

# no_vqsr = [i for i in range(18) if i not in [16]]
# no_vqsr = np.arange(18) # complement for LGB
col_sel = np.arange(X_train_ori.shape[1])
X_train, y_train = X_train_ori[:,col_sel], y_train_ori

Xs_test = [
    X_test30x_ori[msk_test30x][:,col_sel], 
    # X_test50x_ori[:,col_sel]
    ]
ys_test = [
    y_test30x_ori[msk_test30x], 
    # y_test50x_ori
    ]

# load mask of undersample
print('training on 30x data')
clfkit_list, normalizer = check_loadsave(
    path.join(prespath, 'clf', 'clfkits_{}_{}.pkl'.format(seed, balance_strategy)),
    training_series, {
        'X_train': X_train,
        'y_train': y_train,
        'model_list': ['logireg', 'lsvm', 'nn', 'rf', 'xgbdef', 'lgbdef'],
        # 'model_list': ['nn', 'lgbdef'],
        'model_params': [{},{},
            {'hidden_layer_sizes':(50,15)}, 
            {},{},{}],
        'seed': seed,
        'balance_strategy': balance_strategy})

print('testing on 30x')
# print('testing on 30x and 50x')
ys_df, metrs_df = check_loadsave(
    path.join(prespath, 'metr', 'metrs_{}_30x_{}-{:.1f}.pkl'.format(seed, balance_strategy, ratio_test)), 
    evaluating_series, {
        'Xs_test': Xs_test, 
        'ys_test': ys_test,
        'clfkit_list': clfkit_list, 
        'normalizer': normalizer})

