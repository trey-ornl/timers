#ifdef _OPENMP
#include <omp.h>
#endif
#include <unistd.h>
#include "timers.hpp"

int main(int argc, char **argv)
{
  constexpr int ds = 100000;
  MPI_Init(&argc,&argv);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //STOP_TIMER;
  START_TIMER;
  usleep(ds);
  START_TIMER;
  usleep(ds+ds);
#pragma omp parallel
  {
    START_TIMER;
#ifdef _OPENMP
    const int id = omp_get_thread_num();
#else
    const int id = 0;
#endif
    usleep((id+rank)*ds);
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
