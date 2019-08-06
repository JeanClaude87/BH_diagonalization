//	main.cc			     	(C) 2009 Fabio Ortolani fbo 090716
//	==================================================================
#include "version.hh"
#include "storage.hh"
#include "timing.hh"
#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
//
//============================================================================
void	bottom_line ();						  // [main.cc]
void	dmrg_initialize ();					  // [dmrg.cc]
void	dmrg_process ();					  // [dmrg.cc]
void	header_line ();						  // [main.cc]
void	input_collect (int, char * []);				 // [input.cc]
void	input_process ();		       			 // [input.cc]
void	input_show ();						 // [input.cc]
//
//============================================================================
static string	header = NAME"-"VERSION" ("HOSTTYPE") "COPYRIGHT;
static time_t	tstart;
static time_t	tend;
//
//============================================================================
int main (int argc, char * argv [])
{
  srand48(time(0)*getpid());	// initialising random numbers
  tstart = time (0);			// starting time
  atexit (bottom_line);			// print overall statistics on exit
  ios::sync_with_stdio ();		// syncronize c and c++ i/o
  setbuf (stdout, 0);			// unbuffered output
  input_collect (argc, argv);		// collect input
  header_line ();			// output start
  input_show ();			// show collected input lines 
  dmrg_initialize ();			// initialize
  input_process ();			// parse and process input
  dmrg_process  ();			// DMRG algorithm
  return 0;
}
//
//____________________________________________________________________________
void bottom_line ()
{
  //
  //	Prints out an ending line with some statistics
  //
  string 	sout;
  Storage	report;
  //	Collect time statistics ...
  tend = time (0);
  sout = "Cpu " + timecpustr () + ", Run " + timediff (tstart, tend);
  //	Collect storage statistics ...
  sout += ", Storage " + report .usage ();
  //	Unlink swap files
  report .swap_unlink ();
  //	Print bottom line (last output line)
  cout << flush;
  cout << sout << " " << setw (68 - sout .size ()) << left << setfill ('_')
       << "" << setw (0) << " " << timenow () << endl;
}
//
//____________________________________________________________________________
void    header_line ()
{
  //
  //    Prints out a copyright header line
  //
  cout << flush;
  cout <<  left << header << " " << setw(68 - header.size()) << setfill('_') 
       << "" << setfill(' ') << setw(0) << " " << timenow() << endl;
}
//
//============================================================================
