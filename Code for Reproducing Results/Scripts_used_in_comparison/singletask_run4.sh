# seed: any integer
# mutype: snp, indel
# pcode: HCC1395BL_hg38, HCC1395FFPE_hg38, HCC1395_hg38, NA12877
# train_ratio: 0, 200, 100, 50, 1 # 0 means the original ratio
# test_ratio: 0, 200, 100, 50, 1 
. ./config.sh

timestamp=`date +%Y%m%d%H%M`
seed=$1
mutype=$2
pcode=$3
software=$4
gfcode=$5
# balance_strategy=$4

workpath=${result_dir}/${software}/${mutype}/${pcode}
logpath=${workpath}/log
mkdir -p $logpath

$pypath ${script_dir}/testing_imbgf.py \
$seed $mutype $pcode $software $gfcode \
1>${logpath}/${seed}_${gfcode}_${timestamp}.log \
2>${logpath}/${seed}_${gfcode}_${timestamp}.err

# $pypath ${script_dir}/training_balance.py \
# $seed $mutype $pcode $balance_strategy \
# 1>${logpath}/${seed}_${balance_strategy}_${timestamp}.log \
# 2>${logpath}/${seed}_${balance_strategy}_${timestamp}.err

