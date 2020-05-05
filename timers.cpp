#include "timers.hpp"

#ifdef O_TIMERS

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
  Stat(): calls(0), exclusive(0), inclusive(0), on(false), task(0), tasks(0), taskSum(0), thread(0), threads(0), threadSum(0) {}
  long calls;
  double exclusive;
  double inclusive;
  bool on;
  int task;
  int tasks;
  double taskSum;
  int thread;
  int threads;
  double threadSum;
};

struct Line {
  std::string name;
  Stat stat;
  bool operator<(const Line &that) const
  { return stat.exclusive > that.stat.exclusive; }
};

static void getStats(const MPI_Comm comm, const int from, std::map<std::string,Stat> &stats)
{
  int iStops = 0;
  MPI_Recv(&iStops,1,MPI_INT,from,from,comm,MPI_STATUS_IGNORE);
  if (iStops) extraStops = true;
  int nStats = 0;
  MPI_Recv(&nStats,1,MPI_INT,from,from,comm,MPI_STATUS_IGNORE);
  for (int i = 0; i < nStats; i++) {
    int nameSize = 0;
    MPI_Recv(&nameSize,1,MPI_INT,from,from,comm,MPI_STATUS_IGNORE);
    std::string name;
    name.resize(nameSize);
    MPI_Recv(&name.front(),nameSize,MPI_CHAR,from,from,comm,MPI_STATUS_IGNORE);
    Stat &s = stats[name];
    double d[4] = {0,0,0,0};
    MPI_Recv(d,4,MPI_DOUBLE,from,from,comm,MPI_STATUS_IGNORE);
    long l[6] = {0,0,0,0,0};
    MPI_Recv(l,6,MPI_LONG,from,from,comm,MPI_STATUS_IGNORE);
    if (l[1]) s.on = true;
    if (s.exclusive < d[0]) {
      s.exclusive = d[0];
      s.inclusive = d[1];
      s.calls = l[0];
      s.task = l[2];
      s.thread = l[4];
      s.threads = l[5];
      s.threadSum = d[3];
    }
    s.tasks += l[3];
    s.taskSum += d[2];
  }
}

void printTimers(FILE *const out, const MPI_Comm comm)
{
  MPI_Barrier(comm);

  int rank = 0;
  MPI_Comm_rank(comm,&rank);
  int size = 0;
  MPI_Comm_size(comm,&size);

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
        s.task = rank;
        s.tasks = 1;
        s.taskSum = exclusive;
        s.thread = id;
      }
      s.threads++;
      s.threadSum += exclusive;
    }
  }

  if (2*rank+1 < size) getStats(comm,2*rank+1,stats);
  if ((rank > 0) && (2*rank < size)) getStats(comm,2*rank,stats);
  if (rank > 0) {
    const int to = rank/2;
    const int iStops = extraStops;
    MPI_Send(&iStops,1,MPI_INT,to,rank,comm);
    const int nStats = stats.size();
    MPI_Send(&nStats,1,MPI_INT,to,rank,comm);
    for (const auto &pair: stats) {
      const int nameSize = pair.first.size();
      MPI_Send(&nameSize,1,MPI_INT,to,rank,comm);
      MPI_Send(pair.first.data(),nameSize,MPI_CHAR,to,rank,comm);
      const Stat &s = pair.second;
      const double d[4] = {s.exclusive,s.inclusive,s.taskSum,s.threadSum};
      MPI_Send(d,4,MPI_DOUBLE,to,rank,comm);
      const long l[6] = {s.calls,s.on,s.task,s.tasks,s.thread,s.threads};
      MPI_Send(l,6,MPI_LONG,to,rank,comm);
    }
  }

  if (rank == 0) {
    std::vector<Line> lines;
    for (const auto &pair: stats) lines.push_back({pair.first,pair.second});
    std::sort(lines.begin(),lines.end());

    if (extraStops) fprintf(out,"#TIMER *** WARNING: UNMATCHED STOPS\n");
    fprintf(out,"#TIMER exclusive inclusive calls timer\n");
    for (const auto &line: lines) {
      fprintf(out,"#TIMER  %.2f  %.2f  %ld  %s",line.stat.exclusive,line.stat.inclusive,line.stat.calls,line.name.c_str());
      if (size > 1) {
        const double avg = line.stat.taskSum/double(line.stat.tasks);
        const int load = round(100.0*avg/line.stat.exclusive);
        fprintf(out," task %d (%dx%d%%)",line.stat.task,line.stat.tasks,load);
      }
      if (line.stat.threads > 1) {
        const double avg = line.stat.threadSum/double(line.stat.threads);
        const int load = round(100.0*avg/line.stat.exclusive);
        fprintf(out," thread %d (%dx%d%%)",line.stat.thread,line.stat.threads,load);
      }
      if (line.stat.on) fprintf(out," *** WARNING: NOT STOPPED");
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

#endif
