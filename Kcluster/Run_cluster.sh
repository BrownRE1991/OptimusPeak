#!/bin/bash

#rm Cluster_outlists/outASRG*.txt
#rm ../Deconvolution/Datasets/outASRG*.txt
rm -r ../Deconvolution/Datasets
mkdir ../Deconvolution/Datasets
#rm out*.txt
#make
SECONDS=0
./Cluster Outlist.txt 1 10 100 70
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
#mv outASRG*.txt Cluster_outlists
#./Run3.sh
#rm ../Deconvolution/Datasets/outASRG*.txt
#cp Cluster_outlists/outASRG*.txt ../Deconvolution/Datasets/
# rm ../DE/Cluster_outlists/outASRG*.txt
# cp Cluster_outlists/outASRG*.txt ../DE/Cluster_outlists
