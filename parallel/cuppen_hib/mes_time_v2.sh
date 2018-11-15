#!/bin/bash

iniSize=1000
stepSize=1000
limitSize=18000
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

export OMP_NUM_THREADS=6
export MKL_NUM_THREADS=6

rm time_table.txt

for ((size = $iniSize;  size <= $limitSize; size += $stepSize))
do

	for ((i = 0; i < $numexp; i++))
	do

		/usr/lib64/openmpi/1.4-gcc/bin/mpirun -np 2 ./cuppen_par_v4 $size $path_matrices/s_mtab_$size $path_matrices/h_slq_$size 1>>time_table.txt

	done

	echo 1>>time_table.txt

done

