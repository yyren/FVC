from config import RESULT_PATH, DATABASE_PATH
from utils import load_splitdata, balancer_block, save_to_pickle, \
    load_from_pickle, check_loadsave, training_series, evaluating_series, \
        load_30xonly, load_50xonly
from utils import mutypes, pcodes
from os import path, makedirs
from sys import argv

seed = int(argv[1])
mt = argv[2]
pc = argv[3]
balance_strategy = argv[4]


prespath = path.join(RESULT_PATH, mt, pc)


X_train, y_train, X_test30x, y_test30x = check_loadsave(
    path.join(prespath, 'data_suitcase_30x.pkl'), 
    load_30xonly, {'pcode': pc, 'mutype': mt})

X_train50x, y_train_50x, X_test50x, y_test50x = check_loadsave(
    path.join(prespath, 'data_suitcase_50x.pkl'), 
    load_50xonly, {'pcode': pc, 'mutype': mt})

Xs_test = [
    X_test50x, 
    # X_test50x
    ]
ys_test = [
    y_test50x, 
    # y_test50x
    ]

# load mask of undersample
print('training on 30x data')
clfkit_list, normalizer = check_loadsave(
    path.join(prespath, 'clf', 'clfkits_{}_{}.pkl'.format(seed, balance_strategy)),
    training_series, {
        'X_train': X_train,
        'y_train': y_train,
        'model_list': ['logireg', 'lsvm', 'nn', 'rf', 'xgbdef'],
        'seed': seed,
        'balance_strategy': balance_strategy})

print('testing on 50x only')
ys_df, metrs_df = check_loadsave(
    path.join(prespath, 'metr', 'metrs_{}_50x_{}.pkl'.format(seed, balance_strategy)), 
    evaluating_series, {
        'Xs_test': Xs_test, 
        'ys_test': ys_test,
        'clfkit_list': clfkit_list, 
        'normalizer': normalizer})

