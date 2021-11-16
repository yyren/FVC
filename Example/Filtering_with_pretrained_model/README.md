# FNVC filtering with exist pre-trained model
**step1. Features construction**<br>
```bash
singularity exec FNVC_image.sif bash ${absolute_path}/Features_construction.sh \
    -i ${absolute_path}/input.vcf \
    -b ${absolute_path}/input.bam \
    -r ${absolute_path}/hg1kv37.fa \
    -t 4 \
    -p ${absolute_path}/out_snv.record \
    -q ${absolute_path}/out_indel.record
```
The 'input_raw_file.vcf' contains the variants derived from the reads aligned file 'input.bam'.<br>
Please use the 'hg1kv37.fa' in this example

**step2. Filtering**<br>
```bash
singularity exec FNVC_image.sif bash ${absolute_path}/Supervised_learning_filtering.sh \
    -f filter
    -i ${absolute_path}/input.vcf \
    -o ${absolute_path}/input_filtered.vcf \
    -p ${absolute_path}/out_snv.record \
    -q ${absolute_path}/out_indel.record \
    -m ${absolute_path}/gatk_xgbdef_snv.model \
    -n ${absolute_path}/gatk_xgbdef_snv.model
```
pre-trained models: '-m': snv.model; '-n': indel.model<br>