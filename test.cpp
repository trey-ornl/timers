#include <mpi.h>
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
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  //if (rank == size-1) STOP_TIMER;
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
    //if (rank == 2) START_TIMER;
    STOP_TIMER;
  }
  {
    TIME_TO_END_OF_BLOCK;
    for (int i = 0; i < size; i++) {
      START_TIMER;
      usleep(ds);
      STOP_TIMER;
    }
    usleep(ds);
  }
  //STOP_TIMER;
  STOP_TIMER;
  STOP_TIMER;
  printTimers(stdout,MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
