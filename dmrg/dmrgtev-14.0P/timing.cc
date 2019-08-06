//      timing.cc			(C) 2009 Fabio Ortolani	fbo 090716
//	==================================================================
#include "timing.hh"
#include <unistd.h>
#include <time.h>
#include <sys/times.h>
#include <sstream>
#include <iomanip>
//
//============================================================================
static struct tms  actual_cpu;
static struct tm * actual_clock;
static clock_t     uptime;
static clock_t	   tck_cent, tck_sec = 0, tck_min, tck_hour, tck_day;
static long        total_cpu;
static string      timestr;
//
//============================================================================
long timecpu ()
{
  //
  //    Returns the actual cpu time spent by the current process
  //    and child processes in internal clock units.    
  //
  uptime    = times (& actual_cpu);
  total_cpu = actual_cpu.tms_utime  + actual_cpu.tms_stime +
              actual_cpu.tms_cutime + actual_cpu.tms_cstime;
  return total_cpu;
}
//
//____________________________________________________________________________
string timecpustr ()
{
  //
  //    Returns the actual cpu time spent by the current process
  //    and child processes as a human readable string.
  //
  return timecpustr (timecpu());
}
//
//____________________________________________________________________________
string timecpustr (long tck)
{
  //
  //    Returns a human readable string representation of the time 
  //    'tck' given in internal clock units.
  //
  if (! tck_sec) {
    tck_sec     = sysconf (_SC_CLK_TCK);
    tck_cent    = tck_sec  /100;
    tck_min     = tck_sec  * 60;
    tck_hour    = tck_min  * 60;
    tck_day     = tck_hour * 24;
  }
  stringstream s (timestr);
  long days     = tck / tck_day;                tck %= tck_day;
  long hours    = tck / tck_hour;               tck %= tck_hour;
  if (days) {
    s << setw(2) << setfill(' ') << days << "d";
    if (days < 100)  s << setw(2) << setfill('0') << hours << "h";
    return s.str();
  }
  long minutes  = tck / tck_min;                tck %= tck_min;
  if (hours) {
    s << setw(2) << setfill(' ') << hours   << "h" 
      << setw(2) << setfill('0') << minutes << "m" ;
    return s.str();
  }
  long seconds  = tck / tck_sec;                tck %= tck_sec;
  if (minutes) {
    s << setw(2) << setfill(' ') << minutes << "m" 
      << setw(2) << setfill('0') << seconds << "s" ;
    return s.str();
  }
  long cents = tck / tck_cent;
  s << setw(2) << setfill(' ') << seconds << "." 
    << setw(2) << setfill('0') << cents   << "s" ;
  return s.str();
}
//
//____________________________________________________________________________
string timenow ()
{
  //
  //    Returns a string representing the actual watch time of the day
  //
  time_t now = time (NULL);
  actual_clock  = localtime (& now);
  stringstream s (timestr);
  s << setw(2) << setfill('0') << actual_clock -> tm_hour << ":"
    << setw(2) << setfill('0') << actual_clock -> tm_min  << ":"
    << setw(2) << setfill('0') << actual_clock -> tm_sec ;
  return s.str();
}
//       
//____________________________________________________________________________
string timediff (long start, long end)
{
  //
  //    Returns the time difference |end - start| (inputs in seconds) 
  //
  end -= start;
  if (end < 0) end = -end;
  end *=  sysconf(_SC_CLK_TCK);
  return timecpustr (end);
}
//============================================================================
