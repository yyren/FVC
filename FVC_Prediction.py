import argparse
from FVC_utils import load_from_pickle, printx, time, path
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
def predictor_series(X_test, clfkit, normalizer=None):
    """Xs_test should be scaled
    clfkit is a single-item dict
    """
    if normalizer:
        X_test = normalizer.transform(X_test)
    model = list(clfkit.keys())[0]
    clf = clfkit[model]
    printx('########### predicting using {} model ... ############'.format(model))
    t0 = time.time()
    y_pred = clf.predict(X_test)
    y_score = clf.predict_proba(X_test)[:,-1]
    printx('use {:.0f} s\n\n'.format(time.time()-t0))
    return y_pred, y_score


parser = argparse.ArgumentParser(
  description='extract the complex region' )

parser.add_argument('--in_file', type=str, 
    help="snp/indel feature file")

parser.add_argument('--model', type=str, default="pretrain/gatk.model", 
    help="the model trained on your own data or the pretrain model")

parser.add_argument('--out_file', type=str, default="prediction.txt", 
    help="the filtered mutation (labeled by 0) and the propability that it is true mutation")

args = parser.parse_args()

in_file = args.in_file
model_file = args.model 
result_file = args.out_file

mut_feat = pd.read_csv(in_file, header=None, sep=' ', low_memory=False).values
result = pd.DataFrame(mut_feat[:,1:6], columns=['chr','loc0', 'loc1','ref','target'])
X_test = mut_feat[:,6:]
clfkit, normalizer = load_from_pickle(model_file)
y_pred, y_score = predictor_series(X_test, clfkit[0], normalizer)
result['prediction'] = y_pred
result['probability'] = y_score
result.to_csv(result_file, index=False, sep='\t')

