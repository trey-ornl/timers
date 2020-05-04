#ifdef _OPENMP
#include <omp.h>
#endif
#include <unistd.h>
#include "timers.hpp"

int main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);
  //STOP_TIMER;
  START_TIMER;
  sleep(1);
  START_TIMER;
  sleep(2);
#pragma omp parallel
  {
    START_TIMER;
#ifdef _OPENMP
    const int id = omp_get_thread_num();
#else
    const int id = 0;
#endif
    sleep(id);
    //START_TIMER;
    STOP_TIMER;
  }
  //STOP_TIMER;
  STOP_TIMER;
  STOP_TIMER;
  printTimers(stdout,MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
