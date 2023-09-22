#!/bin/bash
#this script complies the C++ code in the pipeline, turns all scripts to executables and returns to the top to run a test.

op=$(uname)
case $op in 

	'Linux')
	echo "Linux"
	;;
	
	'Darwin')
	echo "Darwin"
	;;
	*)
	echo "Operating System not currently supported. Please email author if you want your operating system added."
	;;
esac

chmod +x OptimusPeak.sh

cd Threshold
chmod +x Run2.sh

cd ../Kcluster
make
chmod +x Run_cluster.sh

cd ../Deconvolution
make
chmod +x Run_Decon.sh
chmod +x RunMod_G.sh

cd ../FinalPeaks
chmod +x RunEval.sh

cd ..

cd Threshold
./Run2.sh $PWD/Spectra/ThresholdingTestSpec049_mth1743_Nhsqc.ft2 6 #~37 mins
cd ..

echo "Now implementing symbolic link for Optimus Peak. This will allow Optimus Peak to be called from anywhere in the directory. Would you like this initailized? (y/n)"
read answer

if [ $answer == "y" ]; then
	sudo cp -s $PWD/OptimusPeak.sh /usr/bin/OptimusPeak
fi

if [ $answer == "yes" ]; then
	sudo cp -s $PWD/OptimusPeak.sh /usr/bin/OptimusPeak
fi

echo "Install and Test complete"

