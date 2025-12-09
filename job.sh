#!/bin/bash
#
#SBATCH --job-name=submission
#SBATCH --output=result.txt
#
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=4
#SBATCH --time=10:00

THREADS=${1:-32}

module use /usr/local.nfs/sgs/modulefiles
module load gcc/10.2
module load openmpi/3.1.6-gcc-10.2
module load vtk/9.0.1
module load cmake/3.18.2

srun -n ${THREADS} ./build/numsim_parallel scenarios/lid_driven_cavity.txt
