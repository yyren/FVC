#!/bin/bash
# Description: filter the variant calls or build pre-trained models
# Author: Yongyong Ren
# Usage:
# filtering: bash supervised_learning_filtering.sh -i input.vcf -o output.vcf -s snv_features.record -l indel_features.record -m snv.model -n indel.model
# training:  bash supervised_learning_filtering.sh -s snv_features.record -l indel_features.record -d /model -p HG007_FVC
# output models: /model/HG007_FVC.snv.model /model/HG007_FVC.indel.model


usage() {
	echo ""
	echo "Usage: $0 <OPTIONS>"
	echo "Required Parameters:"
	echo "-f <function>                     filter or train"
	echo "-i <input_vcf>                    unfiltered VCF file"
	echo "-o <output_vcf>                   output filtered VCF file"
	echo "-p <inrecord_snv>                 input record file with SNV features (/data/snv.record)"
	echo "-q <inrecord_indel>               Output file with INDEL features (/data/indel.record)"
	echo "-m <pre-trained snv model>        pre-trained snv model file (/model/HG007_snv.model)"
	echo "-n <pre-trained indel model>      pre-trained indel model file (/model/HG007_indel.model)"
	echo "-a <tp_snv_record>                input record file with True SNV features (/data/tp_snv.record)"
	echo "-b <fp_snv_record>                input record file with False SNV features (/data/fp_snv.record)"
	echo "-c <tp_indel_record>              input record file with True INDEL features (/data/tp_indel.record)"
	echo "-d <fp_indel_record>              input record file with False INDEL features (/data/fp_indel.record)"
	echo "-j <output_snv_model>             output pre-trained model for snv (/model/HG007.model)"
	echo "-k <output_indel_model>             output pre-trained model for indel (/model/HG007.model)"
	echo ""
	exit 1
}

usage_train() {
	echo ""
	echo "Usage: $0 <OPTIONS>"
	echo "Required Parameters:"
	echo "-f <function>                     filter or train"
	echo "-a <tp_snv_record>                input record file with True SNV features (/data/tp_snv.record)"
	echo "-b <fp_snv_record>                input record file with False SNV features (/data/fp_snv.record)"
	echo "-c <tp_indel_record>              input record file with True INDEL features (/data/tp_indel.record)"
	echo "-d <fp_indel_record>              input record file with False INDEL features (/data/fp_indel.record)"
	echo "-j <output_snv_model>             output pre-trained model for snv (/model/HG007.model)"
	echo "-k <output_indel_model>             output pre-trained model for indel (/model/HG007.model)"
	echo ""
	exit 1
}

usage_filter() {
	echo ""
	echo "Usage: $0 <OPTIONS>"
	echo "Required Parameters:"
	echo "-f <function>                     filter or train"
	echo "-i <input_vcf>                    unfiltered VCF file"
	echo "-o <output_vcf>                   output filtered VCF file"
	echo "-p <inrecord_snv>                 input record file with SNV features (/data/snv.record)"
	echo "-q <inrecord_indel>               Output file with INDEL features (/data/indel.record)"
	echo "-m <pre-trained snv model>        pre-trained snv model file (/model/HG007_snv.model)"
	echo "-n <pre-trained indel model>      pre-trained indel model file (/model/HG007_indel.model)"
	echo ""
	exit 1
}


# Parse input arguments
while getopts "h:f:i:o:p:q:m:n:a:b:c:d:j:k:" opts; do
	case "$opts" in
	h) usage;;
	f) function=$OPTARG;;
	i) in_vcf=$OPTARG;;
	o) out_vcf=$OPTARG;;
	p) inrecord_snv=$OPTARG;;
	q) inrecord_indel=$OPTARG;;
	m) in_snv_model=$OPTARG;;
	n) in_indel_model=$OPTARG;;
	a) tp_snv_record=$OPTARG;;
	b) fp_snv_record=$OPTARG;;
	c) tp_indel_record=$OPTARG;;
	d) fp_indel_record=$OPTARG;;
	j) output_snv_model=$OPTARG;;
	k) output_indel_model=$OPTARG;;
	esac
done


# Check input arguments
if [[ "$function" == "" ]]; then
	echo "cannot find the parameter -f"
	exit 1
elif [[ $function == "train" ]]; then
	if [[ $tp_snv_record == "" || $fp_snv_record == "" || $tp_indel_record == "" || $fp_indel_record == "" || $output_snv_model == "" || $output_indel_model == "" ]]; then
		usage_train
	fi
elif [[ $function == "filter" ]]; then
	if [[ $in_vcf == "" || $out_vcf == "" || $inrecord_snv == "" || $inrecord_indel == "" || $in_snv_model == "" || $in_indel_model == "" ]]; then
		usage_filter
	fi
fi
out_path=${out_vcf%/*}
if [[ ! -z $out_vcf ]]; then
	out_path=${out_vcf%/*}
else
	out_path=${out_vcf%/*}
fi
current_bash_dir=$(cd $(dirname ${BASH_SOURCE:-$0});pwd)

# Training or filtering
if [[ $function = 'filter' ]]; then
	# Predicting
	python $current_bash_dir/FVC_Prediction.py --in_file $inrecord_snv --model $in_snv_model --out_file $out_path/FVC_predict_snv.txt
	python $current_bash_dir/FVC_Prediction.py --in_file $inrecord_indel --model $in_indel_model --out_file $out_path/FVC_predict_indel.txt
	
	# Add filtering information into vcf file
	perl $current_bash_dir/merge_predict_results.pl $out_path/FVC_predict_snv.txt $out_path/FVC_predict_indel.txt $in_vcf $out_vcf
elif [[ $function = 'train' ]]; then
	python $current_bash_dir/FVC_Train.py --in_tp $tp_snv_record --in_fp $fp_snv_record --out_model $output_snv_model
	python $current_bash_dir/FVC_Train.py --in_tp $tp_indel_record --in_fp $fp_indel_record --out_model $output_indel_model
	
fi

