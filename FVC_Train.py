import argparse
from FVC_utils import load_from_pickle, save_to_pickle, readData, printx, time, path
from xgboost import XGBClassifier
import warnings
warnings.filterwarnings('ignore')
# from sklearn.preprocessing import StandardScaler

def training_series(X_train, y_train, normalizer=None, xgb_params={}):
    printx('########### training FVC model ... ############')
    printx('training set shape: {}'.format(X_train.shape))
    if normalizer:
        X_train = normalizer.fit_transform(X_train)
    t0 = time.time()
    clf = XGBClassifier(**xgb_params)
    printx(clf.get_params())
    clf.fit(X_train, y_train)
    # printx('model parameters: \n{}'.format(clf.get_params()))
    clfkit = {'FVC':clf}
    printx('use {:.0f} s\n'.format(time.time()-t0))
    return [clfkit], normalizer


parser = argparse.ArgumentParser(
  description='extract the complex region' )

parser.add_argument(
    '--in_tp', type=str, 
    default="training_tp_tensor.record", 
    help="wait to describe")

parser.add_argument(
    '--in_fp', type=str, 
    default="training_fp_tensor.record", 
    help="wait to describe")

parser.add_argument(
    '--model', type=str, 
    default="no",
    help="the absolute path of pretrained model, default no pretrained model.")

parser.add_argument(
    '--out_model', type=str, 
    default="retrain.model", 
    help="wait to describe")

parser.add_argument(
    '--random_seed', type=int, 
    default=0, 
    help="random state, default 0")

parser.add_argument(
    '--scale_strategy', type=str,
    default='standard', choices=['standard', 'no'],
    help='set normalization strategy, default standard scale. \nIf set no, oringal features will be used. \nIf --model is not no, the parameter will be ignored.')

args = parser.parse_args()

train_TP_file = args.in_tp
train_FP_file = args.in_fp
premodel_file = args.model
out_model = args.out_model
random_seed = args.random_seed
scale_strategy = args.scale_strategy

xgb_params = {
    'n_estimators':200, 
    'n_jobs':-1, 
    'random_state':random_seed,
    # 'tree_method':'approx'
}

if premodel_file != 'no':
    printx("loading the pre-trained model: {} ...".format(premodel_file))
    premodel, normalizer = load_from_pickle(premodel_file, '../pretrain')
    xgb_params['xgbdef'] = list(premodel[0].values())[0]
elif premodel_file == 'no' and scale_strategy == 'standard':
    printx("standardly scale all features ...")
    from sklearn.preprocessing import StandardScaler
    normalizer = StandardScaler()
else:
    normalizer = None

printx('load files, {} and {}'.format(train_TP_file, train_FP_file))
X_train, y_train = readData(train_TP_file, train_FP_file)
outmodel, normalizer = save_to_pickle(training_series(X_train, y_train, normalizer=normalizer, xgb_params=xgb_params),
    out_model)


