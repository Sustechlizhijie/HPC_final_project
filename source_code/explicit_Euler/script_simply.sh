#!/bin/bash
#BSUB -J explicit_heat_problem
#BSUB -q ser
#BSUB -n 1


module purge
module load intel/2018.4
module load mpi/intel/2018.4

mpirun -np 1 ./explicit_heat.out > $LSB_JOBID.log 2>&1

