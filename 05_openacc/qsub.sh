#!/bin/sh
#$ -cwd
#$ -l f_node=1
#$ -l h_rt=0:01:00
. /etc/profile.d/modules.sh
module load nvhpc/22.2
make 00_loop
PGI_ACC_TIME=1 ./a.out
