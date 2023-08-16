#!/bin/bash
#SBATCH --partition=hyper
#SBATCH --job-name=720
#SBATCH --nodes=1
#SBATCH --cpus-per-task=112
#SBATCH --ntasks-per-node=1
./main 
