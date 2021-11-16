PLOT_PATH = '/mnt/d/ProjectRYY/plotrevise/plotpdf'
TAB_PATH = '/mnt/d/ProjectRYY/plotrevise/plottab'
DATABASE_PATH = '/mnt/d/ProjectRYY/plotrevise/plotdata'


dark_dict = {
    'red':'#612d21',
    'blue':'#435a97',
    'purple':'#423a75',
    'green':'#4c723c',
    'orange':'#875b31',
    'gray':'#272b29',
    'black':'#10100e'}
normal_dict = {
    'red':'#ab4842',
    'blue':'#6c94cf',
    'purple':'#574d93',
    'green':'#708d62',
    'orange':'#d39d51',
    'gray':'#989e94',
    'black':'#3d3e43'}
light_dict = {
    'red':'#d0a0ac',
    'blue':'#8db4e0',
    'purple':'#7d88b3',
    'green':'#a7d393',
    'orange':'#d39d51',
    'gray':'#ededed',
    'black':'#babcb7'}
special_dict = {
    'red':'#832522',
    'blue':'#2d4b93',
    'orange':'#c38b3c'}
nature_dict = {
    'red':'#832522',
    'blue':'#435a97',
    'orange':'#c38b3c',
    'pink':'#ca7e8d',
    'green':'#4c723c',
    'black':'#101010'}

patient_list = ['HG00'+str(i) for i in [1,3,4,6]]
chr_list = ['chr{}'.format(i) for i in range(1,23)]
depthlist = ['30X', '50X']
order_caller = {'gatk':'GATK', 'varscan': 'Varscan2', 'mutect2':'Mutect2', 'deepvariant':'DeepVariant'}
order_model = {'logireg':'LR', 'lsvm':'LSVM', 'nn':'MLP', 'lgbdef':'LGBM', 'rf':'RF', 'xgbdef':'XGBoost'}
order_strategy = {'original':'original', 'smote':'SMOTE', 'randomover':'oversample', 'nearmiss1':'NearMiss-1', 'randomunder':'undersample'}
order_mutation = {'indel':'INDEL', 'snp':'SNV'}
order_filter = {'FNVC': 'FNVC', 'Frequency': 'Frequency', 'Gar': 'Garfield', 'VQSR': 'VQSR', 'HF': 'Hard Filter', 'VEF': 'VEF'}
order_frequency = {'H':'High Frequency', 'L':'Low Frequency'}
order_consistency = {'inconsistent':'Hard-to-Detect', 'consistent':'Easy-to-Detect'}
order_coding = {'coding':'Coding', 'non-coding':'Non-Coding'}

gatkdeep_feature = {'MQ':0, 'MQRankSum':1, 'MQ0':2, 'MBQ':3, 'BaseQRankSum':4, 'MFRL':5, 'MPOS':6, 'ReadPosRankSum':7, 'RPA':8, 'RPA2':9, 'LRU':10, 'SOR':11, 'FS':12, 'QD':13, 'ExcessHet':14, 'GQ_MEAN':15, 'AS_SB':16, 'AS_UNIQ_ALT_READ_COUNT':17, 'AF':18, 'CT':19}
varscan_feature = {'MQ':0, 'MQRankSum':1, 'MQ0':2, 'MBQ':3, 'BaseQRankSum':4, 'MFRL':5, 'MPOS':6, 'ReadPosRankSum':7, 'RPA':8, 'RPA2':9, 'LRU':10, 'SOR':11, 'FS':12, 'GQ_MEAN':13, 'AS_SB':14, 'AS_UNIQ_ALT_READ_COUNT':15, 'AF':16, 'CT':17}
mutect2_feature = {'MQ':0, 'MQRankSum':1, 'MQ0':2, 'MBQ':3, 'BaseQRankSum':4, 'MFRL':5, 'MPOS':6, 'ReadPosRankSum':7, 'RPA':8, 'RPA2':9, 'LRU':10, 'SOR':11, 'FS':12, 'AS_SB':13, 'AS_UNIQ_ALT_READ_COUNT':14, 'AF':15, 'CT':16}


#Features Group(1: Sequence context; 2: Sequencing experiment; 3: Bioinformatics)
feat_group_id = {'RPA':1, 'RPA2':1, 'LRU':1, 'CT':1, 'MBQ':2, 'MPOS':2, 'ReadPosRankSum':2, 'SOR':2, 'FS':2, 'BaseQRankSum':2, 'MFRL':2, 'AS_SB':2, 'MQ':3, 'MQ0':3, 'MQRankSum':3, 'AS_UNIQ_ALT_READ_COUNT':3, 'ExcessHet':3, 'GQ_MEAN':3, 'QD':3, 'AF':3}


feat_id_dict = {
    'gatk': {v:k for k,v in gatkdeep_feature.items()},
    'deepvariant': {v:k for k,v in gatkdeep_feature.items()},
    'varscan': {v:k for k,v in varscan_feature.items()},
    'mutect2': {v:k for k,v in mutect2_feature.items()}
}

feat_group_name = {
    1: 'Sequence Context',
    2: 'Sequencing Experiment',
    3: 'Bioinformatics'
}

def feature_idx2name(caller, idx):
    return feat_id_dict[caller][idx]