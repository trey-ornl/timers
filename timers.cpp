#include "timers.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <map>
#include <stack>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

struct Timer {
  Timer(): calls(0), children(0), total(0) {}
  long calls;
  double children;
  std::stack<double> starts;
  double total;
};

static std::stack<Timer *> stack;
#pragma omp threadprivate(stack)
static std::map<std::string,Timer> timers;
#pragma omp threadprivate(timers)

static bool extraStops = false;

struct Stat {
  Stat(): calls(0), exclusive(0), inclusive(0), on(false), thread(0) {}
  long calls;
  double exclusive;
  double inclusive;
  bool on;
  int thread;
};

struct Line {
  std::string name;
  Stat stat;
  bool operator<(const Line &that) const
  { return stat.exclusive > that.stat.exclusive; }
};

void printTimers(FILE *const out, const MPI_Comm comm)
{
  MPI_Barrier(comm);

  int rank = MPI_PROC_NULL;
  MPI_Comm_rank(comm,&rank);

  std::map<std::string,Stat> stats;
#pragma omp parallel
  {
#ifdef _OPENMP
    const int id = omp_get_thread_num();
#else
    const int id = 0;
#endif
#pragma omp critical(taskStats)
    for (const auto &pair: timers) {
      const std::string &name = pair.first;
      const Timer &t = pair.second;
      Stat &s = stats[name];
      if (!t.starts.empty()) s.on = true;
      const double exclusive = t.total-t.children;
      if (exclusive > s.exclusive) {
        s.calls = t.calls;
        s.exclusive = exclusive;
        s.inclusive = t.total;
        s.thread = id;
      }
    }
  }

  std::vector<Line> lines;
  for (const auto &pair: stats) lines.push_back({pair.first,pair.second});
  std::sort(lines.begin(),lines.end());

  if (rank == 0) {
    fprintf(out,"#TIMER exclusive inclusive calls timer\n");
    if (extraStops) fprintf(out,"#TIMER WARNING: UNMATCHED STOPS\n");
    for (const auto &line: lines) {
      fprintf(out,"#TIMER  %.2f  %.2f  %ld  %s thread %d",line.stat.exclusive,line.stat.inclusive,line.stat.calls,line.name.c_str(),line.stat.thread);
      if (line.stat.on) fprintf(out," WARNING: STILL RUNNING");
      fprintf(out,"\n");
    }
    fflush(out);
  }

  MPI_Barrier(comm);
}

void startTimer(const std::string &func, const std::string &file, const int line)
{
  const std::string name = func+'('+file+':'+std::to_string(line)+')';
  Timer &t = timers[name];
  stack.push(&t);
  t.calls++;
  t.starts.push(MPI_Wtime());
}

void stopTimer()
{
  if (stack.empty()) {
    extraStops = true;
    return;
  }
  Timer &t = *stack.top();
  stack.pop();
  assert(!t.starts.empty());
  const double dt = MPI_Wtime()-t.starts.top();
  t.starts.pop();
  t.total += dt;
  if (stack.empty()) return;
  stack.top()->children += dt;
}


