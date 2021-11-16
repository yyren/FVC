from config import DATABASE_PATH, RESULT_PATH
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import LinearSVC, SVC
import xgboost as xgb
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
from sklearn.metrics import confusion_matrix, precision_score, recall_score, \
    roc_auc_score, accuracy_score, balanced_accuracy_score, f1_score, \
        matthews_corrcoef, average_precision_score
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
from imblearn.over_sampling import RandomOverSampler, SMOTE, ADASYN
from imblearn.under_sampling import RandomUnderSampler, NearMiss, TomekLinks, EditedNearestNeighbours

# from imblearn.under_sampling import RandomUnderSampler

from copy import copy
import numpy as np
import pandas as pd
from os import path, makedirs
from glob import glob
import pickle
import time

pcodes = ['HG00{}'.format(i) for i in [1,2,3,4,6,7]]
mutypes = ['snp', 'indel']
xgbspe_params = {'booster':'gbtree',
        'objective': 'binary:logistic',
        'eval_metric': 'auc',
        'max_depth':8,
        'n_estimators':500,
        'lambda':15,
        'subsample':0.75,
        'colsample_bytree':0.75,
        'min_child_weight':1,
        'eta': 0.025,
        'seed':0,
        'verbosity':1,
        'gamma':0.15,
        'nthread':12,
        'learning_rate' : 0.01}

metrname = ['AUC', 'AUPRC', 'ACC', 'BACC', 'MCC', 'F1-score',
    'FDR', 'FNR', 'FPR', 'TNR', 'OFO',
    'Sensitivity(Recall)', 'Specificity', 'Precision',
    'TP', 'FP', 'TN', 'FN']

def printx(output):
    print(output, flush = True)


def imb_ratio(y):
    return (y==1).sum()/(y==0).sum()


def save_to_pickle(var, filename='temp', root_path=DATABASE_PATH):
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


def load_from_pickle(filename='temp', root_path=DATABASE_PATH):
    '''
    Load pickle file defaultly from `DATABASE_PATH`
    '''
    printx("loading {} from {}".format(filename, root_path))
    f = open(path.join(root_path, filename), 'rb')
    var = pickle.load(f)
    f.close()
    return var


def check_loadsave(fpath, func, params):
    """check the existance of file (fpath)
    load the var if exist;
    generate and save the var if non-exist.
    """
    if path.exists(fpath):
        var = load_from_pickle(fpath)
    else:
        var = save_to_pickle(func(**params), fpath)
    return var


def readData(file_path1, file_path2):
    # Read file and combine
    data1 = pd.read_csv(file_path1, header=None, sep=' ', low_memory=False)
    data2 = pd.read_csv(file_path2, header=None, sep=' ', low_memory=False)
    
    combined = pd.concat([data1, data2])
    # Return X, y
    return combined.iloc[:, 6:].values, combined.iloc[:, 0].values


def readData_stage1(file_path1, file_path2, file_path3, file_path4):
    # Read file and combine
    data1 = pd.read_csv(file_path1, header=None, sep=' ', low_memory=False)
    data2 = pd.read_csv(file_path2, header=None, sep=' ', low_memory=False)
    data3 = pd.read_csv(file_path3, header=None, sep=' ', low_memory=False)
    data4 = pd.read_csv(file_path4, header=None, sep=' ', low_memory=False)
    hard_region = pd.concat([data3, data4])
    hard_region.ix[:,0]=2
    combined_temp = pd.concat([data1, data2])
    combined = pd.concat([combined_temp, hard_region])
    # Return X, y
    return combined.iloc[:, 6:].values, combined.iloc[:, 0].values


def mcc_score(y_true, y_pred):
    """Matthews correlation coefficient
    Ref: https://en.wikipedia.org/wiki/Matthews_correlation_coefficient
    """
    cmt = confusion_matrix(y_true, y_pred)
    TP, FP, TN, FN = cmt[1,1], cmt[0,1], cmt[0,0], cmt[1,0]
    return (TP*TN-FP*FN) / np.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))


def spe_score(y_true, y_pred):
    TN = ((y_true == 0) & (y_pred == 0)).sum()
    return TN / (y_true == 0).sum()


def sen_score(y_true, y_pred):
    TP = ((y_true == 1) & (y_pred == 1)).sum()
    return TP / (y_true == 1).sum()


def load_splitdata(pcode, mutype='snp'):
    """load specified mutation type data by patient
    Return X_train, y_train, X_test, y_test, X_50Xtest, y_50Xtest
    """
    X_train, y_train = readData(
        path.join(DATABASE_PATH, '30X', pcode, 'Training_FP_{}_annotaion_all.record'.format(mutype)),
        path.join(DATABASE_PATH, '30X', pcode, 'Training_TP_{}_annotaion_all.record'.format(mutype)))
    X_test, y_test = readData(
        path.join(DATABASE_PATH, '30X', pcode, 'Testing_FP_{}_annotaion_all.record'.format(mutype)),
        path.join(DATABASE_PATH, '30X', pcode, 'Testing_TP_{}_annotaion_all.record'.format(mutype)))
    X_50Xtest, y_50Xtest = readData(
        path.join(DATABASE_PATH, '50X', pcode, 'Testing_FP_{}_annotaion_all.record'.format(mutype)),
        path.join(DATABASE_PATH, '50X', pcode, 'Testing_TP_{}_annotaion_all.record'.format(mutype)))
    return X_train, y_train, X_test, y_test, X_50Xtest, y_50Xtest

def load_30xonly(pcode, mutype='indel', software='gatk', load_train=True, root_path=DATABASE_PATH):
    """
    load specified mutation type (default indel) data by patient
    Return X_train, y_train, X_test, y_test (30x)
    """
    if load_train:
        X_train, y_train = readData(
            path.join(root_path, '30X', pcode, 'Training_{}_fp_{}_feature.record'.format(software, mutype)),
            path.join(root_path, '30X', pcode, 'Training_{}_tp_{}_feature.record'.format(software, mutype)))
    else:
        X_train, y_train = None, None
    X_test, y_test = readData(
        path.join(root_path, '30X', pcode, '{}_{}_fp_{}_feature.record'.format(pcode, software, mutype)),
        path.join(root_path, '30X', pcode, '{}_{}_tp_{}_feature.record'.format(pcode, software, mutype)))
    return X_train, y_train, X_test, y_test


def load_50xonly(pcode, mutype='indel', software='gatk', load_train=True, root_path=DATABASE_PATH):
    """
    load specified mutation type (default indel) data by patient
    Return X_train, y_train, X_test, y_test (50x)
    """
    if load_train:
        X_train, y_train = readData(
            path.join(root_path, '50X', pcode, 'Training_{}_fp_{}_feature.record'.format(software, mutype)),
            path.join(root_path, '50X', pcode, 'Training_{}_tp_{}_feature.record'.format(software, mutype)))
    else:
        X_train, y_train = None, None
    X_test, y_test = readData(
        path.join(root_path, '50X', pcode, '{}_{}_fp_{}_feature.record'.format(pcode, software, mutype)),
        path.join(root_path, '50X', pcode, '{}_{}_tp_{}_feature.record'.format(pcode, software, mutype)))
    return X_train, y_train, X_test, y_test


def load_outside(pcode, mutype='indel', software='gatk', root_path=DATABASE_PATH):
    """
    load specified mutation type (default indel) data by patient
    Return X_train, y_train, X_test, y_test (30x)
    """
    # X_train, y_train = readData(
    #     path.join(DATABASE_PATH, pcode, 'Training_{}_fp_{}_feature.record'.format(software, mutype)),
    #     path.join(DATABASE_PATH, pcode, 'Training_{}_tp_{}_feature.record'.format(software, mutype)))
    X_test, y_test = readData(
        path.join(root_path, pcode, '{}_{}_fp_{}_feature.record'.format(pcode, software, mutype)),
        path.join(root_path, pcode, '{}_{}_tp_{}_feature.record'.format(pcode, software, mutype)))
    return X_test, y_test
# TODO load easy-hard pair
# def load_twostagedata(pcode, mutype='snp'):

#     X_train_stage1, y_train_stage1 = readData_stage1(
#         path.join(DATABASE_PATH, '', train_TP_easy),
#         train_FP_easy, 
#         train_TP_hard, 
#         train_FP_hard)
#     X_test_stage1, y_test_stage1 = readData_stage1(test_TP_easy,test_FP_easy, test_TP_hard, test_FP_hard)

    

def random_sample(y, ratio=1, seed=None):
    """the ratio is true sample (1) size to false sample (0) size:
    undersample the false samples if ratio > original ratio;
    undersample the true samples if ratio < original ratio
    """
    print('########## random sampling ... ##########')
    samp1_idx = np.where(y == 1)[0]
    samp0_idx = np.where(y == 0)[0]
    ori_ratio = imb_ratio(y)
    printx('original ratio: {:.2f}'.format(ori_ratio))
    if ratio <= 0:
        pass
    elif ratio > ori_ratio:
        np.random.seed(seed)
        samp0_idx = np.random.choice(samp0_idx, size=round(len(samp1_idx)/ratio), replace=False)
    elif ratio < ori_ratio:
        np.random.seed(seed)
        samp1_idx = np.random.choice(samp1_idx, size=round(len(samp0_idx)*ratio), replace=False)
    sampler_msk = np.repeat(False, len(y))
    sampler_msk[np.hstack((samp1_idx,samp0_idx))] = True
    printx('select {:.2f}% data'.format(sampler_msk.mean()*100))
    printx('get new ratio: {:.2f}'.format(imb_ratio(y[sampler_msk])))
    return sampler_msk

# TODO add more parameters
balancer_dict = {
    'randomover': RandomOverSampler(),
    'smote': SMOTE(n_jobs=-1),
    'adasyn': ADASYN(n_jobs=-1),
    'randomunder': RandomUnderSampler(),
    'nearmiss1': NearMiss(version=1, n_jobs=-1),
    'nearmiss2': NearMiss(version=2, n_jobs=-1),
    'nearmiss3': NearMiss(version=3, n_jobs=-1),
    'tomelinks': TomekLinks(n_jobs=-1),
    'enn': EditedNearestNeighbours(n_jobs=-1)
}


def balancer_block(X, y, seed=None, technique='smote'):
    printx('########## technique-based balancing ... ##########')
    printx('Technique: {}'.format(technique))
    printx('original ratio: {:.2f}'.format(imb_ratio(y)))
    t0 = time.time()
    if technique:
        balancer = balancer_dict[technique.lower()]
        if technique in ['randomover', 'smote', 'adasyn', 'randomunder']:
            balancer.set_params(random_state=seed)
        X_sampled, y_sampled = balancer.fit_sample(X, y)
    else:
        X_sampled, y_sampled = X, y

    printx('new data percentage: {:.2f}%'.format(len(y_sampled)/len(y)*100))
    printx('get new ratio: {:.2f}'.format(imb_ratio(y_sampled)))
    printx('use {:.0f} s\n'.format(time.time()-t0))
    return X_sampled, y_sampled


def unpack_clfkit(clfkit):
    model = list(clfkit.keys())[0]
    clf = clfkit[model]
    printx('model name: {}\nmodel parameters: \n{}'.format(model, clf.get_params()))
    return model, clf


def trainer_block(X_train, y_train, model='lsvm', seed=None, class_weight=None, hidden_layer_sizes=(9,)):
    """X_train should be scaled
    """
    printx('########### training {} model ... ############'.format(model))
    printx('training set shape: {}'.format(X_train.shape))
    t0 = time.time()
    if model == 'xgbspe': # XGBoost style API
        dtrain = xgb.DMatrix(X_train,label=y_train)
        watchlist = [(dtrain,'train')]
        clf = xgb.train(xgbspe_params, dtrain, num_boost_round=500, evals=watchlist)
    else: # sklearn-like API
        if model == 'logireg':
            clf = LogisticRegression(random_state=seed, n_jobs=-1, class_weight=class_weight)
        elif model in ['linear', 'rbf', 'poly', 'sigmoid']:
            clf = SVC(kernel=model, random_state=seed, class_weight=class_weight)
        elif model == 'rf':
            clf = RandomForestClassifier(random_state=seed, n_jobs=-1, class_weight=class_weight)
        elif model == 'xgbdef':
            clf = XGBClassifier(n_estimators=200, random_state=seed, n_jobs=-1, class_weight=class_weight)
        elif model == 'lgbdef':
            clf = LGBMClassifier(random_state=seed, n_jobs=-1, class_weight=class_weight)
        elif 'nn' in model:
            clf = MLPClassifier(hidden_layer_sizes=hidden_layer_sizes, activation='logistic', random_state=seed)
        else:
            clf = LinearSVC(random_state=seed, class_weight=class_weight)
        clf.fit(X_train, y_train)
        printx('model parameters: \n{}'.format(clf.get_params()))
    clfkit = {model:clf}
    printx('use {:.0f} s\n'.format(time.time()-t0))
    return clfkit


def predictor_block(Xs_test, ys_test, clfkit):
    """Xs_test should be scaled
    clfkit is a single-item dict
    """
    model = list(clfkit.keys())[0]
    clf = clfkit[model]
    printx('########### predicting using {} model ... ############'.format(model))
    t0 = time.time()
    if model == 'xgbspe': # XGBoost style API
        dstest = [xgb.DMatrix(X_test) for X_test in Xs_test]
        ys_score = [clf.predict(dtest) for dtest in dstest]
        ys_pred = [(y_score >= 0.5)*1 for y_score in ys_score]
    else: # sklearn-like API
        ys_pred = [clf.predict(_X_test) for _X_test in Xs_test]
        try:
            ys_score = [clf.predict_proba(X_test)[:,-1] for X_test in Xs_test]
        except:
            ys_score = [clf.decision_function(X_test) for X_test in Xs_test]
    printx('use {:.0f} s\n\n'.format(time.time()-t0))
    return ys_pred, ys_score


def evaluator_block(y_true, y_pred, y_score):
    """notice: this block is item-wise function, not by list, and 
    return a dict of all kinds of metrics.
    """
    cmt = confusion_matrix(y_true, y_pred)
    tp = cmt[1,1]
    fp = cmt[0,1]
    tn = cmt[0,0]
    fn = cmt[1,0]
    metr = {
        'AUC': roc_auc_score(y_true, y_score), 
        'AUPRC': average_precision_score(y_true, y_score),
        'ACC': accuracy_score(y_true, y_pred),
        'BACC': balanced_accuracy_score(y_true, y_pred),
        'MCC': matthews_corrcoef(y_true, y_pred),
        'F1-score': f1_score(y_true, y_pred),
        'G1-score': gbeta_(tp, fp, tn, fn),
        'NPV': tn/(tn+fn),
        'FDR': fp/(fp+tp),
        'FNR': fn/(fn+tp),
        'FPR': fp/(fp+tn),
        'TNR': tn/(fp+tn),
        'OFO': fn/tn,
        'Sensitivity(Recall)': sen_score(y_true, y_pred),
        'Specificity': spe_score(y_true, y_pred),
        'Precision': precision_score(y_true, y_pred), 
        # 'recall': recall_score(y_true, y_pred), 
        'TP': tp,
        'FP': fp,
        'TN': tn,
        'FN': fn}
    for ky in metr.keys():
        printx('{}: {:.5f}'.format(ky, metr[ky]))
    printx('compound matrix:\n{}\n\n'.format(cmt))
    return metr


def training_series(X_train, y_train, model_list=None, seed=None, balance_strategy=None, normalize=True, model_params=None):
    if balance_strategy == 'balanced':
        class_weight = 'balanced'
        balance_technique = None
    elif balance_strategy in balancer_dict.keys():
        class_weight = None
        balance_technique = balance_strategy
    else:
        class_weight = None
        balance_technique = None
    X_train, y_train = balancer_block(X_train, y_train, seed, balance_technique)
    if model_list:
        pass
    else:
        model_list = ['logireg', 'lsvm', 'nn', 'rf', 'xgbdef']
    if model_params == None:
        model_params = [{} for _ in model_list]
    normalizer = StandardScaler()
    if normalize:
        X_train = normalizer.fit_transform(X_train)
    clfkit_list = [trainer_block(X_train, y_train, model, seed, class_weight, **params) for model, params in zip(model_list, model_params)]
    return clfkit_list, normalizer


def evaluating_series(Xs_test, ys_test, clfkit_list, normalizer=None):
    # WARNING filter out the specified XGBoost classifier
    clfkit_list = [clfkit for clfkit in clfkit_list if 'xgbspe' not in clfkit.keys()]
    if normalizer:
        Xs_test = [normalizer.transform(X_test) for X_test in Xs_test]
    else:
        Xs_test = [X_test for X_test in Xs_test]
        
    ys_df = [pd.DataFrame({'y_true':y_test}) for y_test in ys_test]
    metrs_dict = [{} for _ in ys_test]
    for clfkit in clfkit_list:
        ys_pred, ys_score = predictor_block(Xs_test, ys_test, clfkit)
        model, _ = unpack_clfkit(clfkit)
        for i,(y_df,metr_dict, y_test,y_pred,y_score) in enumerate(zip(ys_df,metrs_dict, ys_test,ys_pred,ys_score)):
            printx('------ evaling the {}th test set -------'.format(i))
            y_df[model+'_ypred'] = y_pred
            y_df[model+'_yscore'] = y_score
            metr_dict[model] = evaluator_block(y_test, y_pred, y_score)
    return ys_df, [pd.DataFrame(metr_dict).loc[metrname] for metr_dict in metrs_dict]


cond_dict = {
    'cond_1_var': ['{:.1f}-{:.1f}'.format(tr,te) for tr, te in [(1,1),(1,50),(1,100),(1,200),(1,0)]],
    'cond_var_x': ['{:.1f}-{:.1f}'.format(tr,te) for tr, te in [(1,0),(50,0),(100,0),(200,0),(0,0)]],
    'cond_x_var': ['{:.1f}-{:.1f}'.format(tr,te) for tr, te in [(0,1),(0,50),(0,100),(0,200),(0,0)]],
    'resample': ['0.0-0.0', 'balanced', '1.0-0.0', 'nearmiss1', 'randomover', 'smote']}

depth = ['30x', '50x']

def fbeta_(TP, FP, TN, FN, beta=1):
    return 1 / (1 + np.divide((FP + beta**2*FN), (TP + beta**2*TP)))

def gbeta_(TP, FP, TN, FN, beta=1):
    TP, FP, TN, FN = TN, FN, TP, FP
    return fbeta_(TP, FP, TN, FN, beta)

def useful_metr(metr_df, markers_dict, extra_metr=[]):
    pregb_dict = metr_df.loc[['TP', 'FP', 'TN', 'FN']].to_dict()
    gb_dict = {}
    for ky, vl in pregb_dict.items():
        gb = gbeta_(**vl)
        gb_dict[ky] = gb
    usedmetrs_df = pd.concat((metr_df.loc[['AUC', 'MCC', *extra_metr]], pd.DataFrame(gb_dict, index=['G1-score']))).T
    usedmetrs_df['model'] = usedmetrs_df.index.to_numpy()
    for ky,vl in markers_dict.items():
        usedmetrs_df[ky] = vl
    return usedmetrs_df.loc[usedmetrs_df.model != 'xgbspe']

def sigmoid(X):
   return 1/(1+np.exp(-X))

def get_Hlayer(nnclf, X, layer_idx):
    coefs = nnclf.coefs_
    intercepts = nnclf.intercepts_
    Hi = X
    for i in range(layer_idx):
        Hi = sigmoid(Hi @ coefs[i] + intercepts[i])
    return Hi
