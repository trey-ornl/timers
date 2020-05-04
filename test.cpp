#include <unistd.h>
#include "timers.hpp"

int main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);
  START_TIMER;
  sleep(1);
  STOP_TIMER;
  printTimers(stdout,MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
