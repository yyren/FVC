#!/bin/bash
# Description: labeled the variant calls and construct features
# Output: training_tp_snv.record training_fp_snv.record training_tp_indel.record training_fp_indel.record

usage() {
        echo ""
        echo "Usage: $0 <OPTIONS>"
        echo "Required Parameters:"
        echo "-i <input_config_json>      config file with data information (/data/config.json)"
        echo "-t <threads>                cpu cores(about 1.5h/10cores  for each WGS VCF file)"
        echo ""
        exit 1
}

while getopts "h:i:t:" opts; do
	case $opts in
	i) json_file=$OPTARG;;
	t) threads=$OPTARG;;
	h) usage;;
	esac
done

# Parse input arguments
if [[ "$json_file" == "" ]] ; then
    usage
fi

if [[ "$threads" == "" ]] ; then
    echo "Warning - without setting threads, use default: 4"
    threads='4'
fi
if [[ ! -f $json_file ]]; then
	echo "Error: given the absolute path of $json_file"
	exit 1
fi

# get path in singularity image
GATK="/home/bmap/software/GATK/gatk-4.1.9.0/gatk"
rtg_software="/home/bmap/software/RTG/rtg-core-non-commercial-3.10/rtg"
samtools_software='/home/bmap/software/samtools/samtools-1.9/samtools'
current_bash_dir=$(cd $(dirname ${BASH_SOURCE:-$0});pwd)
vcfFiles=`jq '.data | .[] | .unfiltered_VCF' $json_file | sed 's/\"//g'`
reference=`jq '.refSeq' $json_file | sed 's/\"//g'`
outDir_temp=`jq '.outDir' $json_file | sed 's/\"//g'`
outDir=`echo $outDir_temp | sed 's/\/$//'`

if [[ ! -d $outDir ]]; then
	echo "Error - not exist folder: $outDir"
	exit 1
fi

timestamp=`date +'%Y%m%d%H%M%S'`
temp_folder="$outDir_temp/temp_${timestamp}"
mkdir $temp_folder
cd $temp_folder


date_log=`date +'%Y%m%d_%H_%M_%S'`
echo "Start training data construction: $date_log"
#check sdf file
sdf_prefix=${reference%.fa*}
if [[ ! -d ${sdf_prefix}.sdf ]]; then
	echo "Preparing sdf format file: ${sdf_prefix}.sdf"
	$rtg_software format -o ${sdf_prefix}.sdf $reference
fi

idx=0
for input_vcf in ${vcfFiles[*]}
do
	# Add missing features by VariantAnnotator (multiple threads)
	if [[ ! -f $input_vcf ]]; then
		echo "can not find the input_vcf: $input_vcf"
		exit 1
	else
		echo "Processing: $input_vcf"
	fi
	file_name_only=${input_vcf##*/}
	prefix_name=${file_name_only%.*}
	
	idx=$[idx + 1]
	sample_folder=$temp_folder/${prefix_name}_$idx
	mkdir $sample_folder
	bamfile=`jq '.data | .[] | select(.unfiltered_VCF == "'$input_vcf'") | .BAM_file' $json_file | sed 's/\"//g'`
	goldStandard_VCF=`jq '.data | .[] | select(.unfiltered_VCF == "'$input_vcf'") | .goldStandard_VCF' $json_file | sed 's/\"//g'`
	goldStandard_bedFile=`jq '.data | .[] | select(.unfiltered_VCF == "'$input_vcf'") | .goldStandard_bedFile' $json_file | sed 's/\"//g'`
	# running with multiple threads
	metainfor_lines=`grep '^#' $input_vcf | wc -l |sed 's/\r|\n//g'`
	data_lines=`grep -v '^#' $input_vcf | wc -l |sed 's/\r|\n//g'`
	block_lines=`awk -v data_lines=$data_lines -v threads=$threads 'BEGIN{printf "%d", (data_lines/threads)}'`
	echo "$file_name_only|$prefix_name|$bamfile|$goldStandard_VCF|$goldStandard_bedFile|$metainfor_lines|$data_lines|$block_lines"
	
	# adjust the sampleID according to the bam file
	bam_sid=`$samtools_software view -H $bamfile | grep '@RG' |awk 'BEGIN{FS="SM:"}{print $2}'|cut -f 1`
	vcf_sid=`grep '#CHROM' $input_vcf | cut -f 10 | sed 's/\n//'`
	in_file=$input_vcf
	if [[ $bam_sid != $vcf_sid ]]; then
		echo "Warning: The sampleID is not consistent between VCF and Bam"
		echo "Warning: change the SID in vcfID to bamID :$vcf_sid to $bam_sid"
		perl $current_bash_dir/change_header.pl $input_vcf $sample_folder/input_file.vcf $bam_sid
		in_file="$sample_folder/input_file.vcf"
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
		awk -v metainfor_lines=$metainfor_lines -v start_idx=$start_idx -v end_idx=$end_idx '{if((NR<=metainfor_lines) || ((NR>=start_idx) && (NR<=end_idx)))print}' $in_file > $sample_folder/${i}_raw.vcf
	done
	
	echo "Step 1/ VariantAnnotator with multiple threads: $threads"
	## VariantAnnotator and features construction ########
	date_log=`date +'%Y%m%d_%H:%M:%S'`
	echo "Start VariantAnnotator: $date_log"
	for((i=0;i<=$[threads - 1];i++))
	do
		$GATK --java-options "-Djava.io.tmpdir=$sample_folder/" VariantAnnotator -R $reference -I $bamfile -V $sample_folder/${i}_raw.vcf -O $sample_folder/${i}_raw_add_feature.vcf --enable-all-annotations true >$sample_folder/gatk_log_${i}.txt 2>&1 &
	done
	wait
	
	## Merge files
	for((i=0;i<=$[threads - 1];i++))
	do
		if [[ $i = '0' ]]; then
			cat $sample_folder/${i}_raw_add_feature.vcf > $sample_folder/merged_add_features.vcf
		else
			grep -v '^#' $sample_folder/${i}_raw_add_feature.vcf >> $sample_folder/merged_add_features.vcf
		fi
		rm $sample_folder/${i}_raw_add_feature.vcf
		rm $sample_folder/${i}_raw.vcf
	done
	date_log=`date +'%Y%m%d_%H:%M:%S'`
	echo "Finish VariantAnnotator: $date_log"
	
	rtg_data_lines=`grep -v '^#' ${sample_folder}/merged_add_features.vcf | wc -l |sed 's/\r|\n//g'`
	if [[ $rtg_data_lines -lt 1 ]]; then
		echo "Error - VariantAnnotator Failed"
		exit 1
	fi
	
	# filter by bed region and removed the duplicate meta headers
	perl $current_bash_dir/Filter_by_bed.pl $sample_folder/merged_add_features.vcf $goldStandard_bedFile ${sample_folder}/merged_add_features_inbed.vcf ${sample_folder}/merged_add_features_outbed.vcf
	perl $current_bash_dir/Filter_by_bed.pl $goldStandard_VCF $goldStandard_bedFile ${sample_folder}/goldStandard_VCF_inbed.vcf ${sample_folder}/goldStandard_VCF_outbed.vcf
	
	echo "Step 2/4 label the variants"
	# Add labeles (true or false) based on the goldStandard_VCF and the goldStandard_bedFile
	bgzip -c $sample_folder/merged_add_features_inbed.vcf > $sample_folder/merged_add_features_inbed.vcf.gz
	tabix -p vcf $sample_folder/merged_add_features_inbed.vcf.gz
	bgzip -c $sample_folder/goldStandard_VCF_inbed.vcf > $sample_folder/goldStandard_VCF_inbed.vcf.gz
	tabix -p vcf $sample_folder/goldStandard_VCF_inbed.vcf.gz
	
	if [[ -d ${sample_folder}/rtg_out_squash_ploidy ]]; then
		rm -r ${sample_folder}/rtg_out_squash_ploidy
	fi
	$rtg_software vcfeval -c $sample_folder/merged_add_features_inbed.vcf.gz -b $sample_folder/goldStandard_VCF_inbed.vcf.gz -t ${sdf_prefix}.sdf -o ${sample_folder}/rtg_out_squash_ploidy --squash-ploidy --all-records --threads 2
	
	if [[ ! -f ${sample_folder}/rtg_out_squash_ploidy/tp.vcf.gz ]]; then
		echo "Error - RTG vcfeval Failed"
		exit 1
	fi
	gunzip -c ${sample_folder}/rtg_out_squash_ploidy/tp.vcf.gz > ${sample_folder}/tp_add_feature.vcf
	gunzip -c ${sample_folder}/rtg_out_squash_ploidy/fp.vcf.gz > ${sample_folder}/fp_add_feature.vcf
	
	
	# split the labeled vcf file into multiple vcf files
	labeled_files=(tp_add_feature.vcf fp_add_feature.vcf)
	for file in ${labeled_files[*]}
	do
		# running with multiple threads
		metainfor_lines=`grep '^#' ${sample_folder}/$file | wc -l |sed 's/\r|\n//g'`
		data_lines=`grep -v '^#' ${sample_folder}/$file | wc -l |sed 's/\r|\n//g'`
		block_lines=`awk -v data_lines=$data_lines -v threads=$threads 'BEGIN{printf "%d", (data_lines/threads)}'`
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
			awk -v metainfor_lines=$metainfor_lines -v start_idx=$start_idx -v end_idx=$end_idx '{if((NR<=metainfor_lines) || ((NR>=start_idx) && (NR<=end_idx)))print}' ${sample_folder}/$file > $sample_folder/${i}_${file}
		done
	done
	
	date_log=`date +'%Y%m%d_%H:%M:%S'`
	echo "Finish labelling: $date_log"
	echo "Step 3/4 features construction"
	## features construction
	### get features type
	gq_mean='no'
	qd='no'
	variant_caller_type='gatk'
	gq_mean=`awk '{if(/GQ\_MEAN\=/){print "yes";exit}}' $sample_folder/0_fp_add_feature.vcf`
	qd=`awk '{if(/QD\=/){print "yes";exit}}' $sample_folder/0_fp_add_feature.vcf`
	if [[ $gq_mean = 'yes' && $qd = 'yes' ]]; then
		variant_caller_type='gatk'
	elif [[ $gq_mean = 'yes' ]]; then
		variant_caller_type='varscan'
	else
		variant_caller_type='mutect2'
	fi
	
	if [[ ! -f $sample_folder/0_fp_add_feature.vcf ]]; then
		echo "Error - split failed"
		exit 1
	fi
	## extract and construct features
	tags=(tp fp)
	for tag in ${tags[*]}
	do
		for((i=0;i<=$[threads - 1];i++))
		do
			# construct features
			
			python $current_bash_dir/get_vcf_features.py --in_file $sample_folder/${i}_${tag}_add_feature.vcf --out_file $sample_folder/${i}_${tag}.record --caller $variant_caller_type --tag 0 &
		done
		wait
	done
	wait
	
	if [[ -f $sample_folder/tp_merged_features_temp.record ]]; then
		rm $sample_folder/tp_merged_features_temp.record
	fi
	
	if [[ -f $sample_folder/fp_merged_features_temp.record ]]; then
		rm $sample_folder/fp_merged_features_temp.record
	fi
	
	## Merge multiple files into one
	for tag in ${tags[*]}
	do
		for((i=0;i<=$[threads - 1];i++))
		do
			if [[ $i = '0' ]]; then
				cat $sample_folder/${i}_${tag}.record > $sample_folder/${tag}_merged_features_temp.record
			else
				grep -v '^#' $sample_folder/${i}_${tag}.record >> $sample_folder/${tag}_merged_features_temp.record
			fi
			rm $sample_folder/${i}_${tag}_add_feature.vcf
			rm $sample_folder/${i}_${tag}.record
		done
	done
	
	if [[ ! -s $sample_folder/tp_merged_features_temp.record ]]; then
		echo "Error - get_vcf_features.py"
		exit 1
	fi
	## Add feature 'CT'
	for tag in ${tags[*]}
	do
		perl $current_bash_dir/separate_snp_indel.pl $sample_folder/${tag}_merged_features_temp.record $sample_folder/${tag}_merged_features_snv_temp.record $sample_folder/${tag}_merged_features_indel_temp.record feature
		perl $current_bash_dir/add_region_feature.pl $sample_folder/${tag}_merged_features_snv_temp.record $sample_folder/${tag}_merged_features_indel_temp.record $sample_folder/${prefix_name}_${tag}_snv.record $sample_folder/${prefix_name}_${tag}_indel.record
	done
	# remove temp files
	files=(merged_features_temp.record merged_features_snv_temp.record merged_features_indel_temp.record)
	for file in ${files[*]}
	do
		tags=(tp fp)
		for tag in ${tags[*]}
		do
			if [[ -f $sample_folder/${tag}_${file} ]]; then
				rm $sample_folder/${tag}_${file}
			fi
		done
	done
done

# merged the training data from different individuals
traing_files=(Training_tp_snv.record Training_fp_snv.record Training_tp_indel.record Training_fp_indel.record)
for out_file in ${traing_files[*]}
do
	if [[ -f $outDir/$out_file ]]; then
		rm $outDir/$out_file
	fi
done

echo "Step 4/4 merge the training data"
idx=0
for vcf_file in ${vcfFiles[*]}
do
	# Add missing features by VariantAnnotator (multiple threads)
	file_name_only=${vcf_file##*/}
	prefix_name=${file_name_only%.*}
	idx=$[idx + 1]
	sample_folder=$temp_folder/${prefix_name}_$idx
	tags=(tp fp)
	
	for tag in ${tags[*]}
	do
		cat $sample_folder/${prefix_name}_${tag}_snv.record >> $outDir/Training_${tag}_snv.record
		cat $sample_folder/${prefix_name}_${tag}_indel.record >> $outDir/Training_${tag}_indel.record
	done
done

# if [[ -d $temp_folder ]]; then
	# rm -r $temp_folder
# fi

# Finish
date_log=`date +'%Y%m%d_%H:%M:%S'`
echo "Finish constructing data: $date_log"

