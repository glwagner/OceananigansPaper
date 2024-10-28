#!/bin/bash 

#SBATCH --time=48:00:00
#SBATCH --partition=sched_mit_raffaele_gpu
#SBATCH --mem=0
#SBATCH --nodes=1
#SBATCH --gres=gpu:4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH -o omip_simulation_output_%j.txt
#SBATCH -e omip_simulation_error_%j.txt

export JULIA_NUM_PRECOMPILE_TASKS=64
export JULIA_NUM_THREADS=64

export PATH=/orcd/data/raffaele/001/glwagner/Software/julia-1.10.5/bin${PATH:+:${PATH}}
export PATH=/orcd/data/raffaele/001/glwagner/.julia/bin${PATH:+:${PATH}}

mpiexecjl -n 4 julia --project prototype_omip_simulation.jl

