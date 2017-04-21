#!/bin/sh
#PBS -l nodes=1:ppn=1
cd $PBS_O_WORKDIR
NPROCS=`wc -l <$PBS_NODEFILE`
../../../two_pot/bin/fast/two_pot two_pot.in.json | tee two_pot.out
 
