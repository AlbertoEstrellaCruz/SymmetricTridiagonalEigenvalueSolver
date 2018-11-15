#!/bin/bash

iniSize=18000
stepSize=-1000
limitSize=1000
numexp=10

path_matrices=/home/equipo2/maestria/matrices

if [ $# -eq 4 ]
then
	iniSize=$1
	stepSize=$2
	limitSize=$3
	numexp=$4
else
	if [ $# -ne 0 ]
	then
		echo "Usage : mes_time_v2.sh [ <ini> <step> <limit> <numexp> ]"
		exit 1
	fi
fi

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

rm time_table_1t.txt

#make clean
#make

for ((size = $iniSize;  size >= $limitSize; size += $stepSize))
do

	for ((i = 0; i < $numexp; i++))
	do

		/usr/lib64/openmpi/1.4-gcc/bin/mpirun ./cuppen_v4 $size $path_matrices/s_mtab_$size $path_matrices/s_slq_$size 1>>time_table_1t.txt

	done

	echo 1>>time_table_1t.txt

done

