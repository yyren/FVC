from config import RESULT_PATH, DATABASE_PATH, DATABASE_TWO_PATH
from utils import load_from_pickle, check_loadsave, evaluating_series, load_outside
from utils import pcodes
from os import path, makedirs
from sys import argv

'''
seed = 0
mt = 'indel'
opc = 'HCC1395BL_hg38'
software = 'gatk'
gfcode = 'gf'
'''

seed = int(argv[1])
mt = argv[2]
opc = argv[3]
software = argv[4]
gfcode = argv[5]

if gfcode == 'tc':
    root_path = DATABASE_PATH + '/Additional_data'
elif gfcode == 'gf':
    root_path = DATABASE_TWO_PATH + '/Additional_data'


prespath = path.join(RESULT_PATH, software, mt, opc)
X_test, y_test = check_loadsave(
    path.join(prespath, 'data_suitcase_{}.pkl'.format(gfcode)), 
    load_outside, {'pcode': opc, 'mutype': mt, 'software': software, 'root_path':root_path})

Xs_test = [
    X_test
    ]
ys_test = [
    y_test
    ]

for pc in pcodes:
    model_path = path.join(RESULT_PATH, software, mt, pc)
    clfkit_list, normalizer = load_from_pickle(
        path.join(model_path, 'clf', 'clfkits_{}_{}.pkl'.format(seed, gfcode)))
    ys_df, metrs_df = check_loadsave(
        path.join(prespath, 'metr', 'metrs_{}_{}_{}.pkl'.format(seed, gfcode, pc)), 
        evaluating_series, {
            'Xs_test': Xs_test, 
            'ys_test': ys_test,
            'clfkit_list': clfkit_list, 
            'normalizer': normalizer})

