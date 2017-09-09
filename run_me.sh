#!/usr/bin/bash
echo Hi there! The program expects friendly use, meaning that it won\'t check for the presence of the input files.
echo If you don\'t have a rehashed weighted edge list, enter 1 or enter anything else to pass. 
echo Note: If you are running the script for the first time for a dataset, you should consider option 1
read choice
echo Choice entered = $choice
inv_hash_needed=0
if (($choice == 1))
	then echo Please enter the filename of the unweighted edge list which will be automatically rehashed, and weights would be assigned using Jaccard coefficient
	read filename 
	inv_hash_needed=1
	./pre $filename 
	echo
else
	echo Enter the filename of the rehashed weighted edge list of the form \'rehashed_weighted_NAME\'
	read filename
	filename=${filename:18} #to extract the name of the edge list 
fi

echo Enter the stretch \'t\' of the spanner. t = 3, 5, 7, 9, ...
read t

if ((t % 2 == 0)) 
	then echo t must be odd. Quitting the program.
else 
	((k = (t + 1) / 2 ))
	echo Constructing the community graph corresponding to a $t-spanner of \'$filename\' 
	SECONDS=0
	echo
	./final_cpp rehashed_weighted_$filename $k
	t1=$SECONDS
	echo Construction complete in $t1 secs. 
	echo Running modularity maximization on the community graph
	echo
	SECONDS=0
	./convert -i rehashed_$filename -o graph.bin
	p_filename="_cover.part"
	./louvain graph.bin -p rehashed_weighted_$filename$p_filename -v
	echo 
	echo Modularity maximization took $SECONDS secs. The whole process took $(($SECONDS + $t1)) seconds 
	echo Thanks for using the program! 
fi
