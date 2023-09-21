#!/bin/bash

rm -r ../Deconvolution/Datasets
mkdir ../Deconvolution/Datasets

#SECONDS=0
./Cluster Outlist.txt 1 10 100 70
#duration=$SECONDS
#echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
