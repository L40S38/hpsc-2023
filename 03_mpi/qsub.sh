#!/bin/sh
#$ -cwd
#$ -l f_node=1
#$ -l h_rt=0:01:00
. /etc/profile.d/modules.sh
module load gcc openmpi
module load intel-mpi
make 00_hello