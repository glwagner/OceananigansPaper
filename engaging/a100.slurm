#!/bin/bash 

#SBATCH --time=48:00:00
#SBATCH --partition=sched_mit_raffaele_gpu
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --gres=gpu:2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=16
#SBATCH -o omip_simulation_output_%j.txt
#SBATCH -e omip_simulation_error_%j.txt

export JULIA_NUM_PRECOMPILE_TASKS=32
export JULIA_NUM_THREADS=32

export PATH=/orcd/data/raffaele/001/glwagner/Software/julia-1.10.5/bin${PATH:+:${PATH}}
export PATH=/orcd/data/raffaele/001/glwagner/.julia/bin${PATH:+:${PATH}}

#module load mpich
#module load openmpi
#module load cuda

#julia --project -e 'using MPIPreferences; MPIPreferences.use_system_binary()'
julia --project -e 'using MPIPreferences; MPIPreferences.use_jll_binary()'
julia --project -e 'using Pkg; Pkg.build("MPI")'
julia --project -e 'using Pkg; Pkg.build("CUDA")'

# mpiexecjl -n 4 julia --project prototype_omip_simulation.jl
# mpiexecjl -n 2 julia --project test_interpolate.jl

mpiexec -n 2 julia --project test_interpolate.jl

