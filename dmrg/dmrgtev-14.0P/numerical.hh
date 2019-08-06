//	numerical.hh			(C) 2009 Fabio Ortolani	fbo 090716
//	==================================================================
#ifndef NUMERICAL_HH
#define NUMERICAL_HH
#include <cstdlib>
//
//============================================================================
void	householder	  (double *, double *, double *, 
			   size_t, size_t);    		     // [numerical.cc]
void	householder	  (double *, double *, double *, 
			   double *, size_t, size_t);        // [numerical.cc]
double 	machine_precision ();				     // [numerical.cc]
void	random_vector	  (double *, size_t);		     // [numerical.cc]
void	reorder_high_low  (double *, double *, double *,
			   double *, size_t, size_t);  	     // [numerical.cc]
void	reorder_low_high  (double *, double *, double *,
			   double *, size_t, size_t);  	     // [numerical.cc]
size_t	tqli		  (double *, double *, double *,
			   double *, long, long);     	     // [numerical.cc] 
//
//============================================================================
#endif // NUMERICAL_HH
