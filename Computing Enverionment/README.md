# System environment
FVC provides a singularity docker image with all dependent environments. This make the results reproducible and make it easier to be used in different system environment.<br>

Users can download the singularity image from our docker image sofware responsibility<br>


### Docker image for FVC
------------
http://bmap.sjtu.edu.cn/softstorage/details/31

![](../Pictures/FVC_docker_image.png)<br>

### Manually install

##### Python Requirements
* python v3.6.10
* xgboost v1.1.1
* scikit-learn v0.23.0
* pandas v1.0.4
* numpy v1.18.4
* lightgbm v3.1.0
* imbalanced-learn v0.7.0
* re v2.2.1
* argparse v1.1
* os 
* collections 
* math
* itertools
* datetime
* sys

##### Perl Requirements
* perl v5.0
* Getopt::Long
* File::Spec
* FindBin
* other customized models (###.pm) are released in FVC folder

##### Java Requirements
* jdk1.8

##### Others
GATK(version 4.1.9)

### Docker image for variant calling

http://bmap.sjtu.edu.cn/softstorage/details/21

![](Pictures/Variant_calling_docker_image.png)<br>
