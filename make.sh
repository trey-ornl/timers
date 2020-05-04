#!/bin/bash

module load PrgEnv-cray
set -x
CC -g -Wall -Werror test.cpp timers.cpp -ferror-limit=1


