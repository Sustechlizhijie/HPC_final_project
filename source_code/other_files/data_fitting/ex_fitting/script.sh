#!/bin/bash
#BSUB -J explicit_heat_problem
#BSUB -q ser
#BSUB -n 1


module purge
module load intel/2018.4
module load mpi/intel/2018.4

mpirun -np 1 ./explicit_heat.out -ksp_type gmres \
  -ksp_gmres_restart 30 -ksp_rtol 1.0e-10 \
  -ksp_atol 1.0e-50 -ksp_max_it 1500 \
  -ksp_gmres_modifiedgramschmidt \
  -pc_type asm \
  -sub_ksp_type richardson \
  -sub_pc_type icc -ksp_monitor_short \
  -ksp_converged_reason \
  -ksp_view \
  -log_view > $LSB_JOBID.log 2>&1
