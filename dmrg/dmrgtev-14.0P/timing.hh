//	timing.hh			(C) 2009 Fabio Ortolani fbo 090716
//	==================================================================
#ifndef TIMING_HH
#define TIMING_HH
#include <string>
using namespace std;
//
//============================================================================
long		timecpu		();				// [timing.cc]
string		timecpustr	();				// [timing.cc]
string		timecpustr	(long);				// [timing.cc]
string		timenow		();				// [timing.cc]
string		timediff	(long, long);			// [timing.cc]
//
//============================================================================
#endif	// TIMING_HH
