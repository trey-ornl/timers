#pragma once

#include <cstdio>
#include <mpi.h>
#include <string>

void printTimers(FILE *out, MPI_Comm comm);
void startTimer(const std::string &func, const std::string &file, int line);
void stopTimer();

#define START_TIMER startTimer(__func__,__FILE__,__LINE__);
#define STOP_TIMER stopTimer();

