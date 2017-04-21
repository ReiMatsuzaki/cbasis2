#!/bin/sh
#PBS -l nodes=1:ppn=1
cd $PBS_O_WORKDIR
NPROCS=`wc -l <$PBS_NODEFILE`
sh run.sh
 
 
