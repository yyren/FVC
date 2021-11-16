#!/bin/bash
# Description: Features construction module (multiple threads) #########
# Usage: singularity exec FVC_image.sif bash features_construction.sh -i input_raw_file.vcf -s out_snv.record -l out_indel.record -r hg19.fa -b HG001.bam -t 40
# Author: Yongyong Ren

usage() {
        echo ""
        echo "Usage: $0 <OPTIONS>"
        echo "Required Parameters:"
        echo "-i <input_vcf>              VCF file with path (/data/input.vcf)"
        echo "-b <bamfile>                bamfile used in variant calling step (HG001.bam)"
        echo "-r <reference>              refseq in fasta format (hg19.fasta)"
        echo "-p <output_snv_features>    Output file with SNV features (/data/snv.record)"
        echo "-q <output_indel_features>  Output file with INDEL features (/data/indel.record)"
        echo "-t <threads>                cpu cores(about 1.5h/10cores  for each WGS VCF file)"
        echo ""
        exit 1
}

while getopts ":i:p:q:r:b:t:h" opts; do
        case "$opts" in
        i) input_vcf=$OPTARG;;
        p) snv_feature=$OPTARG;;
        q) indel_feature=$OPTARG;;
        r) reference=$OPTARG;;
        b) bamfile=$OPTARG;;
        t) threads=$OPTARG;;
        h) usage;;
        esac
done

# Parse input arguments
if [[ "$input_vcf" == "" || "$snv_feature" == "" || "$indel_feature" == "" || "$reference" == "" || "$bamfile" == "" ]] ; then
    usage
fi

if [[ "$threads" == "" ]] ; then
    echo "Warning - without setting threads, use default: 4"
    threads='4'
fi
GATK="/home/bmap/software/GATK/gatk-4.1.9.0/gatk"
samtools_software='/home/bmap/software/samtools/samtools-1.9/samtools'
current_bash_dir=$(cd $(dirname ${BASH_SOURCE:-$0});pwd)

# running with multiple threads
metainfor_lines=`grep '^#' $input_vcf | wc -l |sed 's/\r|\n//g'`
data_lines=`grep -v '^#' $input_vcf | wc -l |sed 's/\r|\n//g'`
block_lines=`awk -v data_lines=$data_lines -v threads=$threads 'BEGIN{printf "%d", (data_lines/threads)}'`

date_log=`date +'%Y%m%d_%H_%M_%S'`
echo "Start features construction: $date_log"

timestamp=`date +'%Y%m%d%H%M%S'`
out_dir=${snv_feature%/*}
temp_folder="$out_dir/temp_${timestamp}"
echo "Step 1/4 Creat temp folder: $temp_folder"
mkdir $temp_folder
cd $temp_folder

# adjust the sampleID according to the bam file
bam_sid=`$samtools_software view -H $bamfile | grep '@RG' |awk 'BEGIN{FS="SM:"}{print $2}'|cut -f 1`
vcf_sid=`grep '#CHROM' $input_vcf | cut -f 10 | sed 's/\n//'`
in_file=$input_vcf
if [[ $bam_sid != $vcf_sid ]]; then
	echo "Warning: The sampleID is not consistent between VCF and Bam"
	echo "Warning: change the SID in vcfID to bamID: $vcf_sid to $bam_sid"
	perl $current_bash_dir/change_header.pl $in_file $temp_folder/input_file.vcf $bam_sid
	in_file="$temp_folder/input_file.vcf"
fi

for((i=0;i<=$[threads - 1];i++))
do
    s_idx=$[i * block_lines + 1]
    e_idx=$[i * block_lines + block_lines]
    if [[ $i = $[threads - 1] ]]; then
        e_idx=$data_lines
    fi
    start_idx=$[metainfor_lines + s_idx]
    end_idx=$[metainfor_lines + e_idx]
    #echo "$start_idx |$end_idx"
    awk -v metainfor_lines=$metainfor_lines -v start_idx=$start_idx -v end_idx=$end_idx '{if((NR<=metainfor_lines) || ((NR>=start_idx) && (NR<=end_idx)))print}' $in_file > $temp_folder/${i}_raw.vcf
done

###### features construction ####
## VariantAnnotator and features construction ########
date_log=`date +'%Y%m%d_%H_%M_%S'`
echo "Step 2/4 Start VariantAnnotator: $date_log"
for((i=0;i<=$[threads - 1];i++))
do
	$GATK --java-options "-Djava.io.tmpdir=./" VariantAnnotator -R $reference -I $bamfile -V $temp_folder/${i}_raw.vcf -O $temp_folder/${i}_raw_add_feature.vcf --enable-all-annotations true > gatk_log_${i}.txt 2>&1 &
done
wait

annotation_data_lines=`grep -v '^#' $temp_folder/0_raw_add_feature.vcf | wc -l |sed 's/\r|\n//g'`
if [[ $annotation_data_lines -lt 1 ]]; then
	echo "Error - VariantAnnotator Failed"
	exit 1
fi

## get features type
gq_mean='no'
qd='no'
variant_caller_type='gatk'
gq_mean=`awk '{if(/GQ\_MEAN\=/){print "yes";exit}}' $temp_folder/0_raw_add_feature.vcf`
qd=`awk '{if(/QD\=/){print "yes";exit}}' $temp_folder/0_raw_add_feature.vcf`
if [[ $gq_mean = 'yes' && $qd = 'yes' ]]; then
	variant_caller_type='gatk'
elif [[ $gq_mean = 'yes' ]]; then
	variant_caller_type='varscan'
else
	variant_caller_type='mutect2'
fi

## features construction
date_log=`date +'%Y%m%d_%H_%M_%S'`
echo "Step 3/4 Features construction: $date_log"
for((i=0;i<=$[threads - 1];i++))
do
	# construct features
	python $current_bash_dir/get_vcf_features.py --in_file $temp_folder/${i}_raw_add_feature.vcf --out_file $temp_folder/${i}_features.record --caller $variant_caller_type --tag 0 &
done
wait


# Merge files
for((i=0;i<=$[threads - 1];i++))
do
	if [[ $i = '0' ]]; then
		cat $temp_folder/${i}_features.record > $temp_folder/merged_features_temp.record
	else
		grep -v '^#' $temp_folder/${i}_features.record >> $temp_folder/merged_features_temp.record
	fi
	rm $temp_folder/${i}_raw_add_feature.vcf
	rm $temp_folder/${i}_raw.vcf
	rm $temp_folder/${i}_features.record
	rm $temp_folder/gatk_log_${i}.txt
done

# Add feature 'CT'
perl $current_bash_dir/separate_snp_indel.pl $temp_folder/merged_features_temp.record $temp_folder/merged_features_snv_temp.record $temp_folder/merged_features_indel_temp.record feature
perl $current_bash_dir/add_region_feature.pl $temp_folder/merged_features_snv_temp.record $temp_folder/merged_features_indel_temp.record $snv_feature $indel_feature

# remove temp files
if [[ -d $temp_folder ]]; then
	rm -r $temp_folder
fi

date_log=`date +'%Y%m%d_%H_%M_%S'`
echo "Step 4/4 Finish: $date_log"
