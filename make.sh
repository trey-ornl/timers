#!/bin/bash

module load PrgEnv-cray
set -x
CC -fopenmp -g -Wall -Werror test.cpp timers.cpp -ferror-limit=1


