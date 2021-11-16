#!/bin/bash

#get arguments
selected_variant_caller=$1 #gatk or mutect2 or varscan2
docker_image_file=$2 # full path
bam_file=$3 #
reference=$4 # eg. hg19.fa
sampleID=$5 #the SM tag in the bam file, eg. HG001
out_folder=$6 # path to the output files
threads=$7 # provided cpu cores

# SHELL_FOLDER=$(cd "$(dirname "$0")";pwd)

#######################################################
#################### Checking environment ################

if ! type singularity >/dev/null 2>&1; then
    echo 'please install: singularity'
    exit
else
    echo 'already install: singularity'
fi

#######################################################
#################### Checking arguments ##################

file_list=($docker_image_file $bam_file $reference)
for file in ${file_list[*]}
do
    if [[ ! -f $docker_image_file ]]; then
        echo "Not exist: $file"
        exit
    else
        echo "Checking: $file exist"
    fi
done

########################################################
################## Define new argument #################

Samtools="singularity exec $docker_image_file samtools"
GATK="singularity exec $docker_image_file /home/bmap/software/GATK4.0.11/gatk-4.0.11.0/gatk"
GATK_latest="singularity exec $docker_image_file /home/bmap/software/GATK4.1.9/gatk-4.1.9.0/gatk "
Varscan="singularity exec $docker_image_file java -jar /home/bmap/software/varscan/VarScan.v2.3.9.jar"
DeepVariant="singularity exec $deepvariant_image /opt/deepvariant/bin/run_deepvariant"
#######################################################
################# Call variants #######################

cd $out_folder
echo "start analysis HC variant calling"

if [[ $selected_variant_caller = 'gatk' ]]; then
    time $GATK HaplotypeCaller \
   -R $reference -I $bam_file -O ${sampleID}_gatk_HC_raw.vcf
    echo "Finish variant calling"
elif [[ $selected_variant_caller = 'mutect2' ]]; then
    time $GATK_latest --java-options "-Xmx30G -Djava.io.tmpdir=./" Mutect2 \
    -R $reference -I $bam_file -O ${sampleID}_mutect2.vcf \
    -tumor $sampleID --enable-all-annotations true --native-pair-hmm-threads 8
    echo "Finish variant calling"
elif [[ $selected_variant_caller = 'varscan2' ]]; then
    $Samtools mpileup -d 100 -B -f $reference $bam_file | $Varscan mpileup2snp --min-var-freq 0.01 --min-coverage 3 --p-value 0.1 --output-vcf 1 > ${sampleID}_varscan_snp.vcf &
    $Samtools mpileup -d 100 -B -f $reference $bam_file | $Varscan mpileup2indel --min-var-freq 0.01 --min-coverage 3 --p-value 0.1 --output-vcf 1 > ${sampleID}_varscan_indel.vcf
    wait
    echo "Finish variant calling"
elif [[ $selected_variant_caller = 'deepvariant' ]]; then
    if [[ ! -d deepvariant_temp ]]; then
        mkdir deepvariant_temp
    else
        rm -r deepvariant_temp
        mkdir deepvariant_temp
    fi
    if [[ ! -d deepvariant_logs ]]; then
        mkdir deepvariant_logs
    fi
    
    $DeepVariant --model_type=WGS \
       --call_variants_extra_args="use_openvino=true" \
       --ref=$reference \
       --reads=$bam_file \
       --output_vcf=${sampleID}_deepvariant.vcf \
       --output_gvcf=${sampleID}_deepvariant.gvcf \
       --num_shards=$threads \
       --logging_dir=$out_folder/deepvariant_logs \
       --intermediate_results_dir $out_folder/deepvariant_temp
else
    echo "Exit with abnormal selected variant caller: $selected_variant_caller"
    echo "Only support: gatk, varscan2, mutect2"
    exit
fi

