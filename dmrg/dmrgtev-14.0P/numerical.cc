//	numerical.cc			(C) 2009 Fabio Ortolani	fbo 090716
//	==================================================================
//
//	Numerical routines.
//	Main reference: Nunerical Recipes in C++ 2ed. [NR]
//
//============================================================================
#include "numerical.hh"
#include <iostream>
#include <cmath>
using namespace std;
//
//============================================================================
static void donothing () {}	// dirty trick to avoid compiler optimization
//
static double  pythag (const double, const double);	     // [numerical.cc]
//
inline static double square (const double);		     // [numerical.cc]
inline static double sign   (const double, const double);    // [numerical.cc]
//
//============================================================================
static const long max_iterations 	= 60;
//
//============================================================================
void householder (double * a, double * d, double * e, size_t n, size_t na)
{
  //
  //	Transforms symmetric matrix a (n x n) into tridiagonal form
  //	with a sequence of Househoder transformations
  //	
  //	On return arrays d and e contains diagonals and subdiagonal 
  //	elements, a is replaced by global transformation matrix
  //
  size_t i, j, k;
  double eps = machine_precision ();
  if (n < 2) {
    d [0] = a [0]; 
    e [0] = 0.0;
    a [0] = 1.0;
    return;
  }
  for (i = n - 1; i > 0; i--) {
    size_t i_1 = i - 1;
    double h   = 0.0;
    for (j = 0; j < i; ++j) h += square (a [i + j * na]);
    if (h < eps) {
      d [i]   = a [i + i * na];
      e [i_1] = 0.0;
      for (j = 0; j < i; ++j) a [i + j * na] = a [j + i * na] = 0.0;
    }
    else {
      double sigma = sqrt (h);
      if (a [i + i_1 * na] < 0.0) sigma = -sigma;
      d [i]   =  a [i + i * na];
      e [i_1] = -sigma;
      h = h + sigma * a [i + i_1 * na];
      a [i + i_1 * na] += sigma;
      double uau = 0.0;
      for (j = 0; j < i; ++j) {
	double au = 0.0;
	for (k = 0; k < i; ++k) au += a [j + k * na] * a [i + k * na];
	uau += au * a [i + j * na];
	d [j] = au;
      }
      for (j = 0; j < i; ++j) {
	d [j] = d [j] / h - 0.5 * uau * a [i + j * na] / square (h);
	a [j + i * na] = a [i + j * na] / h;
      }
      for (j = 0; j < i; ++j)
	for (k = 0; k < i; ++k) 
	  a [j + k * na] -= (a [i + j * na] * d [k] + d [j] * a [i + k * na]);
    }
  }
  d [0] = a [0];
  a [0] = 1.0;
  for (i = 1; i < n; ++i) {
    for (k = 0; k < i; ++k) {
      double uq = 0.0;
      for (j = 0; j < i; ++j) uq += a [i + j * na] * a [j + k * na];
      for (j = 0; j < i; ++j) a [j + k * na] -= uq * a [j + i * na];
    }
    a [i + i * na] = 1.0;
    for (j = 0; j < i; ++j) a [i + j * na] = a [j + i * na] = 0.0;
  }
} 
//
//____________________________________________________________________________
void householder (double * ar, double * ai, double * d, double * e,
		  size_t n, size_t na)
{
  //
  //	Householder reduction of a hermitian (complex) matrix to tridiagonal
  //	form with a phase transformation to obtain tridiagonal real matrix 
  //	(the resulting offdiagonal elements are real not negative)
  //
  //	On input ar and ai contain real and imaginary matrix elements.
  //
  //	On output the array d and e contain diagonal and offdiagonal 
  //	elements of tridiagonal form.
  //	ar and ai contain the real and imaginary part of corresponding
  //	transformation matrix.
  //	
  //	Asuming that for k > 0 we have a basis |e_j> such that
  //
  //		<e_i, A e_j> = 0 when i < k < j  and j < k < i
  //
  //	we define a vector |u> with components in elements i < k
  //
  //		|u> = sum{j<k} |e_j><e_j,u>
  //
  //	a (complex) phase factor and a new basis:
  //	
  //	|q_i>     =  |e_i>     - 2 /(|u|^2) |u> <u,e_i>,    	  i < k-1
  //	|q_{k-1}> = [|e_{k-1}> - 2 /(|u|^2) |u> <u,e_{k-1}>] eta, |eta| = 1
  //	|q_i>     =  |e_i>,					  i >= k
  //	
  //	definining a unitary transformation Q_k 
  //	
  //		|q_i> = Q_k |e_i>
  //
  //	with the phase factor eta and |u> such that
  //	
  //	<q_i,     A q_k> = 0;	i < k-1
  //	<q_{k-1}, A q_k> = sigma >= 0
  //	
  if (ai == 0) {
    //
    //	For real (symmetric) matrices use real version
    //
    householder (ar, d, e, n, na);
    return;
  }
  //
  double   eps = machine_precision ();
  //
  double * vr  = d;			// temporary (complex) vector |v> 
  double * vi  = e;			// aliased in (d,e)
  //
  size_t i, j, k, km1;
  double etar, etai;
  double sigma = 0;
  for (k = n-1; k > 0; k--) {
    km1   = k-1;
    //
    //	Set resulting diagonal element
    //
    d [k] = ar [k + k *na];
    //
    //	sigma^2 = sum{i<k} |a[k-1,k]|^2
    //	
    double alpha = 0.0;
    for (i = 0; i < k; ++i) 
      alpha += (ar [i + k *na] * ar [i + k *na] + 
		ai [i + k *na] * ai [i + k *na]);
    sigma = sqrt (alpha);
    if (sigma < eps) {
      //
      //  column k is diagonal, no need to change basis
      //
      for (i = 0; i < km1; ++i) 
	ar [i + k *na] = ar [k + i *na] = 
	  ai [i + k *na] = ai [k + i *na] = 0.0;
      ar [k + k *na] = 1.0;
      ai [k + k *na] = 0.0;
      e  [km1] = 0;
      continue;
    }
    //
    //	Phase factor (complex) = (etar,etai) = - a[k-1,k] / |a[k-1,k]|
    //
    etar = -ar [km1 + k *na];
    etai = -ai [km1 + k *na];
    //
    //	eta = |a[k-1,k]|
    //
    double eta  = pythag (etar, etai);	
    if (eta < eps) {
      //
      //  null element, phase factor irrilevant: choose -1
      //
      etar = -1.0;
      etai = 0.0;
      eta  = 0.0;
      ar [km1 + k *na] = 0.0;
      ai [km1 + k *na] = 0.0;
    }
    else {
      etar /= eta;
      etai /= eta;
    }
    //
    //	alpha = |u|^2 / 2 = sigma^2 + sigma |a[k-1,k]|
    //
    alpha = alpha + sigma * eta;
    //	
    //	Store phase factor for next transformations
    //
    ar [k + k *na] = etar;
    ai [k + k *na] = etai;
    //
    //	Store vector |u>:
    //		<e_i,u>     = a[i,k],			i  < k-1
    //		<e_{k-1},u> = a[k-1,k] - sigma eta,	i  = k-1
    //		<e_i,u>	    = 0,			i >= k
    //
    ar [km1 + k *na] -= sigma * etar; 
    ai [km1 + k *na] -= sigma * etai;
    //
    //	Compute:
    //		<e_i, A u> = sum{j<k} <e_i, A e_j> <e_j,u>, i < k
    //		<u, A u>   = sum{i<k} <u,e_i> <e_i, A u>
    //
    double uau = 0.0;
    for (i = 0; i < k; ++i) {
      double aur = 0.0;
      double aui = 0.0;
      for (j = 0; j < k; ++j) {
	aur += (ar [i + j *na] * ar [j + k *na] - 
		ai [i + j *na] * ai [j + k *na]);
	aui += (ar [i + j *na] * ai [j + k *na] + 
		ai [i + j *na] * ar [j + k *na]);
      }
      uau += (ar [i + k *na] * aur + ai [i + k*na] * aui);
      vr [i] = aur;
      vi [i] = aui;
    }
    //
    //	|v> = |A u> / alpha - |u> <u, A u> / (2 alpha^2)
    //
    //	Remember |u> / alpha in row k for transformation
    //
    for (i = 0; i < k; ++i) {
      vr [i] = vr [i] / alpha - 0.5 * uau * ar [i + k *na] / (alpha * alpha);
      vi [i] = vi [i] / alpha - 0.5 * uau * ai [i + k *na] / (alpha * alpha);
      ar [k + i *na] = ar [i + k *na] / alpha;
      ai [k + i *na] = ai [i + k *na] / alpha;
    }
    //	
    //	Transformed submatrix
    //
    //		<q_i, A q_j> = <e_i, {A - |u> <v| - |v> <u|} e_j>
    //
    for (i = 0; i < k; ++i)
      for (j = 0; j < k; ++j) {
	ar [i + j *na] -= (ar [i + k *na] * vr [j] + ai [i + k *na] * vi [j] +
			   vr [i] * ar [j + k *na] + vi [i] * ai [j + k *na]);
	ai [i + j *na] -= (ai [i + k *na] * vr [j] - ar [i + k *na] * vi [j] +
			   vi [i] * ar [j + k *na] - vr [i] * ai [j + k *na]);
      }
    //
    //	Apply phase factor
    //
    for (i = 0; i < k; ++i) {
      double oldr = ar [i + km1 *na];
      double oldi = ai [i + km1 *na];
      ar [i + km1 *na] = (oldr * etar - oldi * etai);
      ai [i + km1 *na] = (oldi * etar + oldr * etai);
      oldr = ar [km1 + i *na];
      oldi = ai [km1 + i *na];
      ar [km1 + i *na] = (oldr * etar + oldi * etai);
      ai [km1 + i *na] = (oldi * etar - oldr * etai);
    }
    //
    //	off diagonal
    //
    e [km1] = sigma;
  }	
  //	
  //	k = 0
  //
  d  [0] = ar [0];
  //
  //	Global transformation
  //
  //		Q = Q_{n-1} * Q_{n-2} * ... * Q_k *... * Q_1
  //
  //		<e_i, Q e_j> = <e_i, q_j> 
  //
  //	The transformation is returned in matrix a
  //
  ar [0] = 1.0;	       
  ai [0] = 0.0;		
  for (k = 1; k < n; ++k) {
    km1 = k - 1;
    //
    //	Apply phase factor to row k-1 
    //	(other elements of row and column k-1 are null from previous step)
    //
    ar [km1 + km1 *na] = ar [k + k *na];
    ai [km1 + km1 *na] = ai [k + k *na];
    //
    //	Accumulate Q_k * Q = Q - |u> <u| Q /alpha 
    //
    for (j = 0; j < k; ++j) {
      //
      //     uq = <u, Q e_j> / alpha
      //
      double uqr = 0.0;
      double uqi = 0.0;
      for (i = 0; i < k; ++i) {
	uqr +=  (ar [k + i *na] * ar [i + j *na] +
		 ai [k + i *na] * ai [i + j *na]);  
	uqi +=  (ar [k + i *na] * ai [i + j *na] -
		 ai [k + i *na] * ar [i + j *na]);
      }
      for (i = 0; i < k; ++i) {
	ar [i + j *na] -= (ar [i + k *na] * uqr - ai [i + k *na] * uqi);
	ai [i + j *na] -= (ai [i + k *na] * uqr + ar [i + k *na] * uqi);
      }
    }
    //	
    //	Next elements are initialized as identity elements
    //
    for (i = 0; i < k; ++i) 
      ar [i + k *na] = ai [i + k *na] = ar [k + i *na] = ai [k + i *na] = 0.0;
    ar [k + k *na] = 1.0;
    ai [k + k *na] = 0.0;
  }
}
//
//____________________________________________________________________________
double machine_precision ()
{
  //
  //	Returns machine precision. 
  //
  static double eps = 1.0;
  if (eps < 1.0) return eps;
  double oneeps	= 2.0;
  double one	= 1.0;
  //
  //	We must perform a dirty trick to avoid compiler optimization
  //	of the while loop
  //
  void (* donothingtrick) () = donothing;
  while (oneeps != one) {
    eps = oneeps - one;
    oneeps  = 0.5 * (one + oneeps);
    donothingtrick ();			// !!! fake the compiler !!!
  }
  return eps;
}
//
//____________________________________________________________________________
static double pythag (const double a, const double b)
{
  //
  //	Computes sqrt (a^2 + b^2) without destructive underflow or overflow
  //	[NR 2.6]
  //
  double absa, absb;
  absa = fabs (a);
  absb = fabs (b);
  if (absa > absb) return absa * sqrt (1.0 + square (absb / absa));
  else return (absb == 0.0 ? 0.0 : absb * sqrt (1.0 + square (absa/absb)));
}
//
//____________________________________________________________________________
void random_vector (double * v, size_t dimension)
{
  //
  //	Fills array v (with length = dimension) with random values 
  //
  if (dimension == 0) return;
  static bool first_time = true;
  size_t i;
  //
  if (first_time) {
    first_time = false;
//    srand48 (43825644L);
  }
  double meanvalue = drand48 () - 0.5;
  double norm      = 0.0;
  for (i = 0; i < dimension; ++i) {
    v [i] = drand48 () + meanvalue;
    norm += v [i] * v [i];
  }
  //
  norm = 1.0 / sqrt (norm);
  for (i = 0; i < dimension; ++i) v [i] *= norm;
}
//
//____________________________________________________________________________
void reorder_high_low (double * a, double * b, double * yr, double * yi, 
		       size_t n, size_t ny)
{
  //
  //	Reorders arrays a, b and columns of matrix y in descending order of a
  //
  for (size_t i = 0; i < n; ++i) {
    double pivot = a [i];
    size_t k = i;
    for (size_t j = i+1; j < n; ++j) 
      if (a [j] > pivot) {
	k = j;
	pivot = a [j];
      }
    if (k > i) {
      a [k] = a [i];
      a [i] = pivot;
      pivot = b [k];
      b [k] = b [i];
      b [i] = pivot;
      if (yr) {
	for (size_t j = 0; j < ny; ++j) {
	  pivot = yr [j + k * ny];
	  yr [j + k * ny] = yr [j + i * ny];
	  yr [j + i * ny] = pivot;
	}
      }
      if (yi) {
	for (size_t j = 0; j < ny; ++j) {
	  pivot = yi [j + k * ny];
	  yi [j + k * ny] = yi [j + i * ny];
	  yi [j + i * ny] = pivot;
	}
      }
    }
  }
}
//
//____________________________________________________________________________
void reorder_low_high (double * a, double * b, double * yr, double * yi, 
		       size_t n, size_t ny)
{
  //
  //	Reorders arrays a, b and columns of matrix y in ascending order of a
  //
  for (size_t i = 0; i < n; ++i) {
    double pivot = a [i];
    size_t k     = i;
    for (size_t j = i+1; j < n; ++j) 
      if (a [j] < pivot) {
	k = j;
	pivot = a [j];
      }
    if (k > i) {
      a [k] = a [i];
      a [i] = pivot;
      pivot = b [k];
      b [k] = b [i];
      b [i] = pivot;
      if (yr) {
	for (size_t j = 0; j < ny; ++j) {
	  pivot = yr [j + k * ny];
	  yr [j + k * ny] = yr [j + i * ny];
	  yr [j + i * ny] = pivot;
	}
      }
      if (yi) {
	for (size_t j = 0; j < ny; ++j) {
	  pivot = yi [j + k * ny];
	  yi [j + k * ny] = yi [j + i * ny];
	  yi [j + i * ny] = pivot;
	}
      }
    }
  }
}
//
//____________________________________________________________________________
inline static double sign (const double a, const double b)
{
  //
  //	Returns |a| times the sign of b
  //	
  return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}
  
//
//____________________________________________________________________________
inline static double square (const double x)
{
  //
  //	Returns x * x
  //
  return x * x;
}
//
//____________________________________________________________________________
size_t tqli (double * d, double * e, double * yr, double *yi, long n, long ny)
{
  //
  //	Tridiagonal simmetric matrix (diagonal d, updiagonal e), using
  //	QL decomposition with implicit shift method. [NR 11.3]
  //	n is the size of arrays d and e. On return the arry d contains
  //	the eigenvalues. 
  //	The iterative transformation are applied to n columns (length ny)
  //	of matrix y.
  //
  //	Return 0 on success, not 0 on error.
  //
  long   i, iter, k, l, m;
  double dd, g, r, c, s, p, f, b;
  //
  if (n < 2)	return 0;
  e [n-1] = 0.0;
  for (l = 0; l < n; ++l) {
    iter = 0;
    do {
      //	
      //	find first null off-diagonal with respect to nearest diagonal
      //
      for (m = l; m < n-1; ++m) {
        dd = fabs (d [m]) + fabs (d [m+1]);
	if ((fabs (e [m]) + dd) == dd) break;
      }
      if (m != l) {
	if (iter++ >= max_iterations) return 1;
	g = (d [l+1] - d [l]) / (2.0 * e [l]);
	r = pythag (g, 1.0);
	g = d [m] - d [l] + e [l] / (g + sign (r, g));
	s = c = 1.0;
	p = 0.0;
	for (i = m-1; i >= l; i--) {
	  f = s * e [i];
	  b = c * e [i];
	  r = pythag (f, g);
	  e [i + 1] = r;
	  if (r == 0.0) {
	    d [i + 1] -= p;
	    e [m]      = 0.0;
	    break;
	  }
	  s = f / r;
	  c = g / r;
	  g = d [i + 1] - p;
	  r = s * (d [i] - g) + 2.0 * c * b;
	  p = s * r;
	  d [i + 1] = g + p;
	  g = c * r - b;
	  for (k = 0; k < ny; ++k) {
	    f = yr [k + (i+1) * ny];
	    yr [k + (i+1) * ny] =  c * f + s * yr [k + i * ny];
	    yr [k +     i * ny] = -s * f + c * yr [k + i * ny];
	  }
	  if (yi) {
	    for (k = 0; k < ny; ++k) {
	      f = yi [k + (i+1) * ny];
	      yi [k + (i+1) * ny] =  c * f + s * yi [k + i * ny];
	      yi [k +     i * ny] = -s * f + c * yi [k + i * ny];
	    }
	  }
	}
	if (r == 0.0 && i >= l) continue;
	d [l] -= p;
	e [l]  = g;
	e [m]  = 0.0;
      }
    } while (m != l); 
  }	
  return 0;
}
//
//============================================================================
