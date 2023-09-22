#!/bin/bash

#if previous iterations output is present, remove it
if [ -f "Threshold/Outlist.bin" ]; then
	rm Threshold/Outlist.bin
fi
w=$1
n=$2

#if not given a specific number of cores, then use nproc
if [ -z "$2" ]
	then
    	echo "No argument supplied"
    	n=$(nproc)
	#if given a number of cores, but that number is less than 1, set number of cores to 6
	else
		if [ $2 -lt 1 ]
		then
			n=6
		else
			n=$2
	fi
fi

rm Thresholding_output.txt

#threshold Spectrum
python3 ThresholderGuas2.py $1 #>> Thresholding_output.txt

#if thresholding was successful, continue
if [ -f "Threshold/Outlist.bin" ]; then
	cp "Threshold/Outlist.bin" "../Kcluster/Outlist.bin"
	cd ../Kcluster
	./Run_cluster.sh
	cd ../Deconvolution
	./Run_Decon.sh $w $n
	cd ../FinalPeaks
	./RunEval.sh $w 0
fi

#go back to start
cd ../Threshold
