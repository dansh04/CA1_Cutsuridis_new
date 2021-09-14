#!/bin/bash -l

#$ -P riseprac

#$ -pe omp 8

#$ -l h_rt=1:00:00

#$ -m ea

#$ -N ca1mem_sim1

#$ -j y
#$ -o log_ca1mem_sim1.qlog

echo "============================================"
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID: $JOB_ID $SGE_TASK_ID"
echo "============================================"

module purge
module load python3/3.6.5
module load openmpi
module load neuron/7.6.7

mpirun -np 8 ./x86_64/special -mpi -c '"sim1"' -c "0.0"  main.py