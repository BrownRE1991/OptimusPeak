#!/bin/bash
#make

#a=$(lscpu | grep "socket" | cut -d':' -f2)

echo $2
if [ -z "$2" ]
	then
    	echo "No argument supplied"
    	a=$(nproc)
		((a=a / 2 -1))
		b=$(nproc)
		((b= b / 2))
	else
		if [ $2 -eq 0 ]
		then
			a=6
		else
			a=$2
			b=$2
		fi
		((a= a / 2))
		((a=a-1))
		((b= b / 2))
fi


echo $a
echo $b


SECONDS=0
rm -r Models_Decon
mkdir Models_Decon
rm -r ../FinalPeaks/Models_Decon
mkdir ../FinalPeaks/Models_Decon
#rm ../FinalPeaks/Models_Decon/Models_2ASRG*.txt

cd Datasets
n=$(ls | wc -l)
echo $n
cd ..

if [ $n -gt $b ]; then
	for i in `seq 0 $a`; do
    	./RunMod_G.sh $b $i &
    	pids[${i}]=$!
    	#echo "spawned $i"
	done
fi

if [ $n -le $b ]; then
	if [ $n -gt 0 ]; then
		((a=n-1))
		echo $a
		for i in `seq 0 $a`; do
    		./RunMod_G.sh $n $i &
    		pids[${i}]=$!
    		#echo "spawned $i"
		done
	fi
fi

for pid in ${pids[*]}; do
    wait $pid
done

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed. Spectra $1" >> Times.txt

