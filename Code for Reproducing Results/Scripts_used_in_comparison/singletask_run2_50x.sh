# seed: any integer
# mutype: snp, indel
# pcode: HG001, HG002, HG003, HG004, HG006, HG007
# train_ratio: 0, 200, 100, 50, 1 # 0 means the original ratio
# test_ratio: 0, 200, 100, 50, 1 
. ./config.sh

timestamp=`date +%Y%m%d%H%M`
seed=$1
mutype=$2
pcode=$3
# ratio_train=$4
# ratio_test=$5
balance_strategy=$4

workpath=${result_dir}/${mutype}/${pcode}
logpath=${workpath}/log
mkdir -p $logpath

# $pypath ${script_dir}/training_imbalance.py \
# $seed $mutype $pcode $ratio_train $ratio_test \
# 1>${logpath}/${seed}_${ratio_train}-${ratio_test}_${timestamp}.log \
# 2>${logpath}/${seed}_${ratio_train}-${ratio_test}_${timestamp}.err

$pypath ${script_dir}/training_balance_50x.py \
$seed $mutype $pcode $balance_strategy \
1>${logpath}/${seed}_${balance_strategy}_${timestamp}.log \
2>${logpath}/${seed}_${balance_strategy}_${timestamp}.err

