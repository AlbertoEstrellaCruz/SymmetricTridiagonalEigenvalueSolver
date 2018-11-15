#!/bin/bash

iniSize=1000
stepSize=1000
limitSize=18000

path_matrices=/home/equipo2/maestria/matrices

rm nohup.out time_table.txt

for ((size = $iniSize;  size <= $limitSize; size += $stepSize))
do

	/usr/lib64/openmpi/1.4-gcc/bin/mpirun -np 2 ./cuppen_par_v4 $size $path_matrices/s_mtab_$size $path_matrices/h_slq_$size 1>>time_table.txt

done

