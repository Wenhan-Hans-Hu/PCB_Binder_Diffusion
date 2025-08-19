#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=4gb:ngpus=1
module load miniforge/3
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate /rds/general/user/wh621/home/miniforge3/envs/diffusion
cd $PBS_O_WORKDIR
echo 'Starting at:' $(date +%F) $(date +%T)
 
conda env create -y -f mlfold.yml

echo 'Finishing at:' $(date +%F) $(date +%T)
