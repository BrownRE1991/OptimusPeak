#!/bin/bash
#make

#a=$(lscpu | grep "socket" | cut -d':' -f2)

op=$(uname)
case $op in 

	'Linux')
	echo "Linux"
	;;
	
	'Darwin')
	echo "Darwin"
	;;
	*)
	echo "Operating System not currently supported for multiprocessing. Please email author if you want your operating system added."
	;;
esac

#check if number of cores supplied. If no cores supplied, check operating system and use use appropriate call to find the number of cores in the system, and then use half. If the number of core supplied is less than 1, default to 6 cores.
echo $2
if [ -z "$2" ]
	then
    	echo "No argument supplied"
    	case  $op in
    	
    		"Linux")
    		a=$(nproc) #this one pulls number of logical cores. I found people using linux tend to use it as a work computer only, so using all logical cores makes slightly more sense here as a default. You can ask it to use less.
			((a=a / 2 -1))
			b=$(nproc)
			((b= b / 2))
			;;
			'Darwin')
			a=$(sysctl -n hw.physicalcpu) #this one pulls the number of physical cores. 
			((a=a / 2 -1))
			b=$(sysctl -n hw.physicalcpu)
			((b= b / 2))
			;;
			*)
			echo "Operating System not currently supported for multiprocessing. Please email author if you want your operating system added."
			;;
		esac
	else
		if [ $2 -lt 1 ]
		then
			a=6
			b=6
		else
			a=$2
			b=$2
			case  $op in
    	
    		"Linux")
    		z=$(nproc) #this one pulls number of logical cores. I found people using linux tend to use it as a work computer only, so using all logical cores makes slightly more sense here as a default. You can ask it to use less.
    		if [ $2 -lt $z ] 
    		then
    			a=$2
    			b=$2
    		else
    			a=$z
    			b=$z
    		fi
			;;
			'Darwin')
			z=$(sysctl -n hw.physicalcpu) #this one pulls the number of physical cores. 
			if [ $2 -lt $z ] 
			then
    			a=$2
    			b=$2
    		else
    			a=$z
    			b=$z
    		fi
			;;
			*)
				a=6
				b=6
			;;
		esac
		fi
		((a= a / 2))
		((a=a-1))
		((b= b / 2))
fi


echo $a
echo $b


#SECONDS=0 #If you want to time how long it takes put this line and lines 103 and 104 back in. 


#ensure the previous runs files don't mess with the current run by deleting all previous files.
rm -r Models_Decon
mkdir Models_Decon
rm -r ../FinalPeaks/Models_Decon
mkdir ../FinalPeaks/Models_Decon


#set up several processes to run the deconvolution algorithm. Each takes the script number and the total number of cores to be able to tell which files to process. 
cd Datasets
n=$(ls | wc -l) #the number of clusters to be processed
echo $n
cd ..

#if the number of clusters is greater than the number of cores to be used, set up all scripts
if [ $n -gt $b ]; then
	for i in `seq 0 $a`; do
    	./RunMod_G.sh $b $i &
    	pids[${i}]=$!
	done
fi

#if the number of clusters is less than the number of cores to be used, use only enough scripts to run each cluster. You don't always know how many clusters there will be.
if [ $n -le $b ]; then
	if [ $n -gt 0 ]; then
		((a=n-1))
		echo $a
		for i in `seq 0 $a`; do
    		./RunMod_G.sh $n $i &
    		pids[${i}]=$!
		done
	fi
fi

#wait for all processes to complete
for pid in ${pids[*]}; do
    wait $pid
done

#duration=$SECONDS
#echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed. Spectra $1" >> Times.txt

