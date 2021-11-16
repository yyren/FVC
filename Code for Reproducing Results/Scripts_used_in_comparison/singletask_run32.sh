# seed: any integer
# mutype: snp, indel
# pcode: HG001, HG002, HG003, HG004, HG006, HG007
# train_ratio: 0, 200, 100, 50, 1 # 0 means the original ratio
# test_ratio: 0, 200, 100, 50, 1 
. ./config.sh

timestamp=`date +%Y%m%d%H%M`
seed=0
mutype=$1
pcode=$2
ratio_train=0
ratio_test=0
software=gatk
# balance_strategy=$4

workpath=${result_dir}/${software}/${mutype}/${pcode}
logpath=${workpath}/log
mkdir -p $logpath

$pypath ${script_dir}/gatk_teston50.py \
$seed $mutype $pcode $ratio_train $ratio_test $software \
1>${logpath}/${seed}_${ratio_train}-${ratio_test}_50x_${timestamp}.log \
2>${logpath}/${seed}_${ratio_train}-${ratio_test}_50x_${timestamp}.err

# $pypath ${script_dir}/training_balance.py \
# $seed $mutype $pcode $balance_strategy \
# 1>${logpath}/${seed}_${balance_strategy}_${timestamp}.log \
# 2>${logpath}/${seed}_${balance_strategy}_${timestamp}.err

