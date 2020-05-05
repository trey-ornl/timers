#pragma once

#ifdef O_TIMERS

#include <cstdio>
#include <mpi.h>
#include <string>

void printTimers(FILE *out, MPI_Comm comm);
void startTimer(const std::string &func, const std::string &file, int line);
void stopTimer();

struct TimerGuard {
  TimerGuard(const std::string &func, const std::string &file, int line)
  { startTimer(func,file,line); }
  ~TimerGuard() { stopTimer(); }
};

#define START_TIMER startTimer(__func__,__FILE__,__LINE__);
#define STOP_TIMER stopTimer();
#define TIME_TO_END_OF_BLOCK const TimerGuard timerGuard_ ## __LINE__ (__func__,__FILE__,__LINE__);

#else

#define START_TIMER
#define STOP_TIMER
#define TIME_TO_END_OF_BLOCK

#define printTimers(X,Y) {}

#endif
