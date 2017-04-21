#!/bin/sh
#PBS -l nodes=1:ppn=1
cd $PBS_O_WORKDIR
NPROCS=`wc -l <$PBS_NODEFILE`
../../../../rhf/bin/fast/rhf rhf.in.json | tee rhf.out.json 
 
 
