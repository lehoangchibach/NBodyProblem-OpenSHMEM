#!/bin/bash
#SBATCH -n 1680
#SBATCH -N 35
#SBATCH --ntasks-per-node=48

srun --mpi=pmi2 -o $HOME/comp380/lab3/output/1680.out.%j $HOME/comp380/lab3/lab3 -n 300000 -t 10
