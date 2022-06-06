#!/bin/bash
#BSUB -J valgrind_explicit_heat_problem
#BSUB -q ser
#BSUB -n 1


module purge
module load intel/2018.4
module load mpi/intel/2018.4
module load valgrind/3.14.0

valgrind mpirun ./explicit_heat.out > $LSB_JOBID-valgrind.log 2>&1
