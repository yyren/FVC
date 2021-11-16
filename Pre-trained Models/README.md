# Pre-trained models
The pre-trained models used for filtering variants identified by GATK HaplotypeCaller, Mutect2, Varscan2, and DeepVariant were released in this folder.<br>

Training datasets: HG001, HG003, HG004, and HG006. <br>

Algorithm
------------
These models were trained by using XGBoost with the parameters in the 'Code for Reproducing Results/Scripts_used_in_comparison/ML_parameters.json'
XGBoost demonstrated a significant improvement than other machine learning (ML) methods in both SNV classification and Indel classification.
![](https://github.com/yyren/FVC/raw/master/Picture/Comparison_of_different_ML_methods.jpg)<br>


Evaluation
------------

These models were tested on an additional test dataset (HG007).
![](https://github.com/yyren/FVC/raw/master/Picture/HG007_AUC.jpg)<br>

