# seed: any integer
# mutype: snp, indel
# pcode: HG001, HG002, HG003, HG004, HG006, HG007
# train_ratio: 0, 200, 100, 50, 1 # 0 means the original ratio
# test_ratio: 0, 200, 100, 50, 1 
# seed=0
# mutype=snp
# pcode=HG001
# train_ratio=1
# test_ratio=0

. ./config.sh

seed=0
software=gatk
# ratio_train=$1
# ratio_test=$2
balance_strategy=$1
# balanced
# randomover
# smote
# nearmiss1
# tomelinks

for mutype in 'snp' 'indel'; do
    for pcode in 'HG001' 'HG002' 'HG003' 'HG004' 'HG006' 'HG007'; do
        # ${script_dir}/singletask_run.sh $seed $mutype $pcode $ratio_train $ratio_test
        ${script_dir}/singletask_run2.sh $seed $mutype $pcode $balance_strategy $software
    done
done

