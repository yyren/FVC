# FNVC filtering without pre-trained model

**step1. Constructing training data and features** <br>
```bash
singularity exec FNVC_image.sif bash ${absolute_path}/Data_construction.sh \
    -i ${absolute_path}/training_data_config.json \
    -t 40
```

**step2. Supervised learning** <br>
```bash
singularity exec FNVC_image.sif bash ${absolute_path}/Supervised_learning_filtering.sh \
    -f train \
    -a ${outDir}/Training_tp_snv.record \
    -b ${outDir}/Training_fp_snv.record \
    -c ${outDir}/Training_tp_indel.record \
    -d ${outDir}/Training_fp_indel.record \
    -j ${outDir}/pipeline_adapted_snv.model \
    -k ${outDir}/pipeline_adapted_indel.model
```

**step3. Constructing FNVC features for the testing data** <br>
```bash
singularity exec FNVC_image.sif bash ${absolute_path}/Features_construction.sh \
    -i ${absolute_path}/input.vcf \
    -b ${absolute_path}/input.bam \
    -r ${absolute_path}/hg1kv37.fa \
    -t 40 \
    -p ${absolute_path}/out_snv.record \
    -q ${absolute_path}/out_indel.record
```
-t: cpu cores (about 23min/40cores  for each WGS VCF file), we suggest use as many as you can in this step.<br>
The 'input_raw_file.vcf' contains the variants derived from the reads aligned file 'input.bam'.<br>

**step4. Filtering** <br>
```bash
singularity exec FNVC_image.sif bash ${absolute_path}/Supervised_learning_filtering.sh \
    -f filter
    -i ${absolute_path}/input.vcf \
    -o ${absolute_path}/input_filtered.vcf \
    -p ${absolute_path}/out_snv.record \
    -q ${absolute_path}/out_indel.record \
    -m ${absolute_path}/pipeline_adapted_snv.model \
    -n ${absolute_path}/pipeline_adapted_indel.model
```