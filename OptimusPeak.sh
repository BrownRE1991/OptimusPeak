#!/bin/bash

#current=$PWD
#echo $current

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

spectrum=$1
output=$2
n=$3

echo $spectrum
echo $output
echo $n

if test -f $spectrum; then
	cd Threshold
	./Run2.sh $spectrum $n
	cd ..
	cp FinalPeaks/Combined_Results.csv $output
else
	echo "$spectrum not found."
	return 1
fi

echo "Optimus Peak completed successfully."



