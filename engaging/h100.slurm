#!/bin/bash 

#SBATCH --time=48:00:00
#SBATCH --partition=pi_raffaele
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --gres=gpu:2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=16
#SBATCH -o omip_simulation_output_%j.txt
#SBATCH -e omip_simulation_error_%j.txt

export JULIA_NUM_PRECOMPILE_TASKS=64
export JULIA_NUM_THREADS=64

export PATH=/orcd/data/raffaele/001/glwagner/Software/julia-1.10.5/bin${PATH:+:${PATH}}
export PATH=/orcd/data/raffaele/001/glwagner/.julia/bin${PATH:+:${PATH}}

mpiexecjl -n 2 julia --project prototype_omip_simulation.jl
