#!/bin/bash
#this script complies the C++ code in the pipeline, turns all scripts to executables and returns to the top to run a test.

#chmod +x Run_all_Steffen.sh
#chmod +x Run_all_Allan.sh
#chmod +x Run_all.sh
chmod +x OptimusPeak.sh

cd Threshold
chmod +x Run2.sh

##cluster the spectra using ASRG, copy clusters to DE and GoodnessOfFitTests
cd ../Kcluster
make
chmod +x Run_cluster.sh
#
##run deconvolution algorithm on all clusters
cd ../Deconvolution
make
chmod +x Run_Decon.sh
chmod +x RunMod_G.sh
#chmod +x RunMod0.sh
#chmod +x RunMod1.sh
#chmod +x RunMod2.sh
#chmod +x RunMod3.sh
#chmod +x Run_oneCluster.sh
#chmod +x Run_oneProcessor.sh

cd ../FinalPeaks
chmod +x RunEval.sh

cd ..

cd Threshold
./Run2.sh 49 6 #37 mins
#./Run2.sh 20 6 #383 mins
cd ..

sudo cp -s $PWD/OptimusPeak.sh /usr/bin/OptimusPeak

echo "Install and Test complete"
