#!/bin/sh
#$ -cwd
#$ -l f_node=1
#$ -l h_rt=0:01:00
. /etc/profile.d/modules.sh
module purge
module load cuda
make 00_hello
./a.out
