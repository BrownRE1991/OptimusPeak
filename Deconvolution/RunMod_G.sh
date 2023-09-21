#!/bin/bash

w=1
file=Datasets/outASRG$w.txt
cd Datasets
n=$(ls | wc -l)

cd ..
if [ $n -gt 0 ]; then
	#echo "Datasets/*.txt"
	while [ ! -f $file ]
	do
		((w=w+1))
		file=Datasets/outASRG$w.txt
	done
fi


rm thing$2.txt

while [ -f $file ]
do
	current=$(($w%$1))
	if [ $current -eq $2 ]; then
		if [ -f "Datasets/outASRG$w.txt" ]; then
			if [ -f "genes$2.txt" ]; then
				rm genes$2.txt
			fi
			./K-models Datasets/outASRG$w.txt $2 >> thing$2.txt
			if [ -f "genes$2.txt" ]; then
				#cp -fr genes$2.txt Models_Decon/Models_2ASRG$w.txt
				cp -fr genes$2.txt ../FinalPeaks/Models_Decon/Models_2ASRG$w.txt
			fi
		fi
		echo "model: $w"
 	fi
	((w=w+1))
	file=Datasets/outASRG$w.txt
done
