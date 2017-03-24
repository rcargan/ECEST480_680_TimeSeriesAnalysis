#!/bin/bash
#$ -j y
#$ -cwd
#$ -M rac89@drexel.edu
#$ -l h_rt=24:00:00
#$ -P rosenclassPrj
#$ -pe shm 64
#$ -l mem_free=12G
#$ -l h_vmem=16G
#$ -q all.q

export R_LIBS=/mnt/HA/groups/rosenclassGrp/r_libs

 . /etc/profile.d/modules.sh
 module load shared
 module load proteus
 module load sge/univa
 module load gcc/4.8.1
 module load bowtie2/2.2.5

 R CMD BATCH /mnt/HA/groups/rosenclassGrp/TimeSeriesGroup/dada_sequencetable_script_rc.R

 exit 0

