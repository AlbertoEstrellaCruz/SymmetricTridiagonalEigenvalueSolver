#!/bin/bash

iniSize=18000
stepSize=-1000
limitSize=1000

path_matrices=/home/equipo2/maestria/matrices

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

make clean
make

rm time_table.txt

for ((size = $iniSize;  size >= $limitSize; size += $stepSize))
do

	/usr/lib64/openmpi/1.4-gcc/bin/mpirun ./dstedc_v1 $size $path_matrices/s_mtab_$size $path_matrices/dstedc_slq_$size 1>>time_table.txt

done

