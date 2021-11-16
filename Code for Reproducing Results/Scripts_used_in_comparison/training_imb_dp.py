from config import RESULT_PATH, DATABASE_PATH
from utils import load_splitdata, random_sample, save_to_pickle, \
    load_from_pickle, check_loadsave, training_series, evaluating_series, \
        load_30xonly, readData
from utils import mutypes, pcodes
from os import path, makedirs
from sys import argv

'''
seed = int(1)
mt = mutypes[0]
pc = pcodes[0]
ratio_train = float(1)
ratio_test = float(1)

'''

seed = 0
mt = 'indel'
pc = 'HG007'
ratio_train = 0
ratio_test = 0

prespath = path.join(RESULT_PATH, 'test', mt, pc)
inst_path = '/data1/TumorGroup/DATA/public_database/FVC/Training/annotation_DP/HG007'

X_train, y_train = readData(
    path.join(inst_path, 'Training_gatk_tp_indel_feature.record'),
    path.join(inst_path, 'Training_gatk_fp_indel_feature.record'))
X_test, y_test = readData(
    path.join(inst_path, 'HG007_gatk_tp_indel_DP_feature.record'),
    path.join(inst_path, 'HG007_gatk_fp_indel_DP_feature.record'))

# X_train, y_train = X_train_ori[msk_train], y_train_ori[msk_train]
Xs_test = [X_test]
ys_test = [y_test]


enrol_col = [i for i in range(X_train.shape[1]) if i not in [16, 17]]

print('training on 30x data')
clfkit_list, normalizer = check_loadsave(
    path.join(prespath, 'clf', 'clfkits_{}_nodpvqsr_{:.1f}.pkl'.format(seed, ratio_train)),
    training_series, {
        'X_train': X_train[:,enrol_col],
        'y_train': y_train,
        'model_list': ['logireg', 'lsvm', 'nn', 'rf', 'xgbdef'],
        'seed': seed})

print('testing on 30x only')
ys_df, metrs_df = check_loadsave(
    path.join(prespath, 'metr', 'metrs_{}_nodpvqsr_30x_{:.1f}-{:.1f}.pkl'.format(seed, ratio_train, ratio_test)), 
    evaluating_series, {
        'Xs_test': [_X_test[:,enrol_col] for _X_test in Xs_test], 
        'ys_test': ys_test,
        'clfkit_list': clfkit_list, 
        'normalizer': normalizer})




####################################################################################
# /data1/TumorGroup/DATA/public_database/FVC/Training/annotation_DP/HG007

inst_path = '/export/home/zhouxiaocheng/project/test/snapshot20201204_ryy/main'
seed = 0
mt = 'indel'
pc = 'HG007'
ratio_train = 0
ratio_test = 0

prespath = path.join(RESULT_PATH, 'test', mt, pc)

X_train, y_train, X_test, y_test, X_test50x, y_test50x = load_from_pickle(
    path.join(inst_path, mt, pc, 'data_suitcase.pkl')
)

# X_train, y_train = X_train_ori[msk_train], y_train_ori[msk_train]
Xs_test = [X_test, X_test50x]
ys_test = [y_test, y_test50x]


enrol_col = [i for i in range(X_train.shape[1]) if i not in [1]]

print('training on 30x data')
clfkit_list, normalizer = check_loadsave(
    path.join(prespath, 'clf', 'clfkits_{}_10fnodp_{:.1f}.pkl'.format(seed, ratio_train)),
    training_series, {
        'X_train': X_train[:,enrol_col],
        'y_train': y_train,
        'model_list': ['xgbdef'],
        'seed': seed})

print('testing on 30x only')
ys_df, metrs_df = check_loadsave(
    path.join(prespath, 'metr', 'metrs_{}_10fnodp_{:.1f}-{:.1f}.pkl'.format(seed, ratio_train, ratio_test)), 
    evaluating_series, {
        'Xs_test': [_X_test[:,enrol_col] for _X_test in Xs_test], 
        'ys_test': ys_test,
        'clfkit_list': clfkit_list, 
        'normalizer': normalizer})

