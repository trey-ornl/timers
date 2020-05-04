#!/bin/bash

module load PrgEnv-cray
export MV2_USE_CUDA=1
export MV2_ENABLE_AFFINITY=0
ulimit -c unlimited
rm -f core*
srun -n1 ./a.out


