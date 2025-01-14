#!/bin/bash

#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=16
#SBATCH --mem=0
#SBATCH --time 48:00:00
#SBATCH -o test_output_%j.txt
#SBATCH -e test_error_%j.txt

# Upload modules
module purge all
module add spack
module add cuda
module load openmpi/3.1.6-cuda-pmi-ucx-slurm-jhklron

# MPI specific exports
export OMPI_MCA_pml=^ucx
export OMPI_MCA_osc=^ucx
export OMPI_MCA_btl_openib_allow_ib=true

# Julia specific enviromental variables
export JULIA_CUDA_MEMORY_POOL=none

# Number of threads in SLURM mode
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}

dir='/home/glwagner/OceananigansPaper/sandbox'

mpiexec -n 4 julia --project $dir/test_interpolate.jl

