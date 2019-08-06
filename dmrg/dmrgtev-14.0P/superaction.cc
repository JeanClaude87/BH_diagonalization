//	superaction.cc		       	(C) 2011 Fabio Ortolani fbo 110201
//	==================================================================
#include "timing.hh"
#include "numerical.hh"
#include "superaction.hh"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <iomanip>
//
//============================================================================
//
//============================================================================
static const size_t csize = 					// [action.cc]
  (sizeof (complex<double>) + sizeof (double) -1) / sizeof (double);
double initial_residual	= 0.0;
double initial_value   	= 0.0;
//
//============================================================================
Biaction::Biaction ()
  : ba_size 	(0),  
    ba_lftop	(),
    ba_rgtop	(),
    ba_lftab	(),
    ba_rgtab	()
{ 
  //
  //	Default ctor.
  //
}
//
//____________________________________________________________________________
Biaction::Biaction (const Biaction & other)
  : ba_size	(other .ba_size),
    ba_lftop	(other .ba_lftop),
    ba_rgtop	(other .ba_rgtop),
    ba_lftab	(other .ba_lftab),
    ba_rgtab	(other .ba_rgtab)
{ 
  //
  //	Copy ctor.
  //
}
//
//____________________________________________________________________________
Biaction::~Biaction ()
{ 
  //
  //	Default dtor.
  //
}
//
//============================================================================
Bipoli::Bipoli ()
  : bi_coeff 	(0.0,0.0),
    bi_lft	(complex<double> (1.0, 0.0)),
    bi_rgt	(complex<double> (1.0, 0.0))
{ 
  //
  //	Default ctor.
  //
}
//
//____________________________________________________________________________
Bipoli::Bipoli (const Bipoli & other)
  : bi_coeff	(other .bi_coeff),
    bi_lft	(other .bi_lft),
    bi_rgt	(other .bi_rgt)
{ 
  //
  //	Copy ctor.
  //
}
//
//____________________________________________________________________________
Bipoli::~Bipoli ()
{ 
  //
  //	default dtor.
  //
}
//
//============================================================================
Lanczos::Lanczos  ()
  : l_vector	 	(0),
    l_value		(0),
    l_tolerance		(0),
    l_error		(0),
    l_size		(0),
    l_found		(0),
    l_strategy		(0),
    l_steps		(0),
    l_repeats		(0),
    l_zerotolerance	(0.0),
    l_zeronorm		(0.0)
{ 
  //
  //	Default ctor.
  //
}
//
//____________________________________________________________________________
Lanczos::~Lanczos ()
{
  //
  //	Non standard destructor.
  //
  //	Destroy allocated memory
  //
  if (l_vector) delete [] l_vector;
  if (l_value)  delete [] l_value;
}
//
//____________________________________________________________________________
void Lanczos::krylov (double * a, double * b, double * y, 
		      size_t n, size_t m)
{
  //
  //	Diagonalization in Krylov subspace 
  //
  size_t i, j;
  for (j = 0; j < m; j++) {
    for (i = 0; i < m; i++) y [i + j * m] = 0.0;
    y [j + j * m] = 1.0;
  }
  if (n > 1) {
    for (j = 0; j < n; j++) {
      y [j + j * m] = a [j];
      y [n + j * m] = y [j + n * m] = b [j];
    }
    y [n + n * m] = a [n]; 
    householder (y, a, b, n + 1, m);
  }
  if (tqli (a, b, y, 0, m, m)) {
    cout << "Error in Krylov subspace diagonalization!" << endl;
    exit (0);
  }
  reorder_low_high (a, b, y, 0, m, m);
}
//
//____________________________________________________________________________
double * Lanczos::operator [] (size_t index)
{
  //
  //	Array-like access to memory of index-th vector
  //	with allocation if needed.
  //
  if (l_vector [index] .storage () == 0) { 
    l_vector [index] .storage (l_size * sizeof (double));
    ((complex<double> *) (l_vector [index] .storage ())) [0] = 1.0;
  }
  return (double *) l_vector [index] .storage ();
}
//
//____________________________________________________________________________
void Lanczos::parameters (size_t memsize, size_t ransize,
			  size_t wanted, size_t  dimension,
			  size_t strategy, size_t steps, size_t repeats,
			  double zero_value, double zero_norm)
{
  //
  //	Set Lancos parameters
  //
  //	Size (in double units) of Lanczos vectors
  //
  l_size   = memsize;
  l_random = ransize;
  //
  //	Number of wanted (found) lowest eigenvalues
  //
  l_found  = wanted;
  if (wanted > dimension) l_found = dimension;
  //
  //	Convergency strategy:
  //		0:		only lowest eigenvalues
  //		1:		convergency on lowest, keep highest
  //		2:		convergency on lowest and highest
  //
  l_strategy = strategy % 3;
  //
  //	Maximum number of Lanczos algorithm steps
  //
  l_steps = steps;
  if (l_steps < 4 * l_found + 24)  l_steps = 4 * l_found + 24;
  if (l_steps > dimension) l_steps = dimension;
  //
  //	Maximum number of restarts of Lanczos algoritm
  //
  l_repeats = repeats;
  //
  //	Acceptance tolerance and residual of eigenvalues and eigenvectors
  //
  l_zerotolerance = zero_value;
  l_zeronorm      = zero_norm;
  //
  //	Prepare to store Lanczos vectors (Krylov subspace base)
  //	Make virtual space more larger than l_steps vector only 
  //	to accomodate deflated vectors (at most l_steps) and some 
  //	working space.
  //	
  if (l_vector)	    delete [] l_vector;
  l_vector  = new Storage [3 * l_steps + 1];
  if (l_vector == 0) {
    cout << "No room for Lanczos Storage's!" << endl;
    exit (0);
  }
  memset (l_vector, 0, (3 * l_steps + 1) * sizeof (Storage));
  //
  //	Resulting eigenvalues, tolerances, residuals
  //
  if (l_value)	    delete [] l_value;
  l_value         = new double [l_steps * 6];
  if (l_value == 0) {
    cout << "No room for Lancos values!" << endl;
  }
  memset (l_value, 0, 6 * l_steps * sizeof (double));
  l_tolerance 	  = l_value +     2 * l_steps;
  l_error	  = l_tolerance + 2 * l_steps;
}
//
//____________________________________________________________________________
void Lanczos::random (size_t index)
{
  //
  //	Set random values in index-th vector
  //
  random_vector ((*this) [index] + csize, l_random - csize);
}
//
//____________________________________________________________________________
void Lanczos::storage (size_t index, const Storage & other)
{
  //
  //	Set vector [index] from other
  //
  size_t ssize = l_size * sizeof (double);
  if (other .size () != ssize) {
    cout << "Invalid Lanczos vector [" << index << "] initialization!" 
	 << endl;
    exit (0);
  }
  l_vector [index] << other;
}
//
//____________________________________________________________________________
size_t Lanczos::thick (Superaction & superaction, Action & guess)
{
  //
  //	Lanczos diagonalization using thick restart algorithm
  //
  size_t i, j, k, m, mold, lfound, hfound;
  //
  long tm = timecpu ();
  long ta = 0;
  //
  // machine precsion times small factor
  //
  double eps = machine_precision () * 10.0;	
  //
  //	Allocate memory for Lanczos coefficients, and Krylov eigenvectors
  //
  size_t  memsize = l_steps * l_steps + 6 * l_steps;
  Storage memory (memsize * sizeof (double));
  //
  //	The matrix of eigenvectors (in the Lancos vectors base)
  //
  double  * y        = (double *) memory .storage ();
  //	
  //	Lanczos coefficients (then eigenvalues)
  //
  double  * a        = y + l_steps * l_steps; // Ritz values
  //
  //	Lanczos off diagonal coefficients (residual norms)
  //
  double  * b        = a + 2 * l_steps;	// Lanczos offdiagonal coefficients
  //
  //	Service space
  //
  double  * c	     = b + 2 * l_steps;
  //
  //	Final eigenvalues, tolerances, and residuals (errors of eigenvalues)
  //
  double  * aold     = l_value;	 	  // eigenvalues
  double  * delta    = l_tolerance;
  double  * residual = l_error;
  //
  size_t applyes    = 0;
  size_t iterations = 0;
  //
  size_t scalars    = 0;
  size_t locals	    = 0;
  size_t globals    = 0;
  size_t randoms    = 0;
  size_t start	    = 0;
  size_t deflated   = 0;
  size_t deflate    = 0;
  size_t low        = l_found;
  size_t lowextra   = 6;
  size_t high	    = l_found;
  size_t highextra  = 6;
  if (l_strategy == 1) {
    high = 0;
    highextra = l_found + 6;
  }
  if (l_strategy == 0) 	high = highextra = 0;
  size_t lwanted = low;
  size_t hwanted = high; 
  //
  //	Setup for scalar products and complex multiple
  //	
  vector<Abinary> scalar;
  binary (scalar, guess, guess);
  complex<double> z;
  if (guess .iscomplex ()) z = complex<double> (1.0, 1.0);
  else			   z = 1.0;
  Action dummy (guess .range (), z);
  dummy .scalar (z);
  vector<Abinary> multiple;
  binary (multiple, guess, dummy, guess);
  //
  Lanczos & vector = *this;
  double  * v;
  double  * q	  = vector [0];
  double  * w	  = guess .storage ();
  size_t    dim   = guess .dimension ();
  //
  double alpha 	    = 0.0;
  double beta	    = 0.0;
  double beta2	    = 0.0;
  double eta	    = 0.0;
  //
  for (i = csize; i < l_size; i++) q [i] = w [i];  
  //
  size_t kmax = l_steps;
  deflate = deflated = 0;
  mold = m = lfound = hfound = 0;
  //
  while (iterations++ < l_repeats) {
    //
    for (k = start; k < kmax; k++) {
      m = k + 1;
      w = vector [m]; 
      //
      //	Initialization
      //
      if (k > start) {
	q = vector [k-1];
	for (i = csize; i < l_size; i++) w [i] = -beta * q [i];
	release (k-1);
      }
      else {
	for (i = csize; i < l_size; i++) w [i] = 0.0;
	for (j = deflated; j < start; j++) {
	  v = vector [j];
	  for (i = csize; i < l_size; i++) w [i] -= b [j] * v [i];
	  release (j);
	}
      }
      //
      q = vector [k];                   // q = v_k is orthogonal to w = v_{k-1}
      long te = timecpu ();
      biapply (w, superaction, q);	// w  --->  w + superaction * q
      ta += (timecpu () - te);
      applyes++;
      //
      //	Computing diagonal matrix element
      //
      alpha = 0.0;
      for (i = csize; i < l_size; i++) alpha += q [i] * w [i];
      if (applyes == 1) initial_value = alpha;
      scalars++;
      //
      //	Ortogonalization to obtain the residual
      //
      for (i = csize; i < l_size; i++) w [i] -= alpha * q [i];
      eta = alpha * alpha;
      if (k > start) eta += beta * beta;
      else for (j = 0; j < start; j++) eta += b [j] * b [j];
      //
      for (i = csize, beta2 = 0.0; i < l_size; i++) beta2 += w [i] * w [i];
      scalars++;
      //
      //	Are norms near zero?
      //
      if (beta2 < eps) beta2 = 0.0;
      if (beta  < eps) beta  = 0.0;
      if (eta   < eps) eta   = 0.0;
      //
      //	Re-orthogonalization
      //
      if (beta2 > eta) {
	//
	//	Local re-orthogonalization (little correction)
	//
	locals++;
	eta = 0.0;
	for (i = csize; i < l_size; i++) eta += q [i] * w [i];
	alpha -=eta;
	for (i = csize; i < l_size; i++) w [i] -= eta * q [i];
	if (k > start) {
	  v = vector [k-1];
	  z = multiply (v, w, scalar); 
	  scalars++;
	  multiply (w, -z, v, multiple);
	  release (k-1);
	}
	//	
	//	Normalize
	//
	for (i = csize, beta2 = 0.0; i < l_size; i++) beta2 += w [i] * w [i];
	scalars++;
	if (beta2 > eps) {
	  beta = sqrt (beta2);
	  for (i = csize; i < l_size; i++) w [i] /= beta;
	}
	else beta = beta2 = 0.0;
      }
      else if (beta2 > eps * eta) {
	//
	//	Global re-othogonalization (global orthogonality correction)
	//
	globals++;
	eta = 0.0;
	for (i = csize; i < l_size; i++) eta += q [i] * w [i];
	alpha -= eta;
	for (i = csize; i < l_size; i++) w [i] -= eta * q [i];
	for (j = 0; j < k; j++) {
	  v = vector [j];
	  z = multiply (v, w, scalar);
	  scalars++;
	  multiply (w, -z, v, multiple);
	  release (j);
	}
	//	
	//	Normalize
	//
	for (i = csize, beta2 = 0.0; i < l_size; i++) beta2 += w [i] * w [i];
	scalars++;
	if (beta2 > eps) {
	  beta = sqrt (beta2);
	  for (i = csize; i < l_size; i++) w [i] /= beta;
	}
	else beta = beta2 = 0.0;
      }
      else beta = beta2 = 0.0;
      if (beta2 <= eps) {
	//
	//	Krylov space exausted. 
	//	Generate a random vector (orthogonal to actual Krylov 
	//	space) to continue algorithm.
	//	
	randoms++;
	random (m);
	for (j = 0; j < m; j++) {
	  v = vector [j];
	  z = multiply (v, w, scalar);
	  scalars++;
	  multiply (w, -z, v, multiple);
	  if (j != k) release (j);
	}
	//	
	//	Normalize
	//
	for (i = csize, beta2 = 0.0; i < l_size; i++) beta2 += w [i] * w [i];
	scalars++;
	if (beta2 > eps) {
	  beta = sqrt (beta2);
	  for (i = csize; i < l_size; i++) w [i] /= beta;
	}
	else beta2 = 0.0;
	beta = 0.0;
      }
      //
      //	Lanczos coefficients
      //
      a [k] = alpha;
      b [k] = beta;
      if (applyes == 1) initial_residual = beta;
      //
      //	A little cleaning
      //
      if ((k > start) && (b [k-1] < eps * (fabs (a [k-1]) + fabs (a [k])))) 
	b [k-1] = 0.0;
      release (m);
      //
      //	If no residual vector terminate 
      //
      if (beta2 == 0.0) break;
      //
    } // for (k = start; ...
    //
    //  Mark last deflated vectors with a special value for the residual
    //  to recognize them in output.
    //
    for (j = 0; j < deflated; j++) b [j] = 7.77e-77;
    //
    //	Record number of deflated vectors in last iteration
    //
    deflate = deflated;
    //
    //	Diagonalization in the Krylov subspace
    //
    krylov (a + deflated, b + deflated, y, start - deflated, m - deflated);
    //
    // 	Evaluate tolerances and residuals
    //
    for (j = deflated; j < m; j++) {
      b [j] = beta * y [m - deflated - 1 + (j - deflated) * (m - deflated)];
      delta [j] = 0.0;
      if (mold) {
	// 
	//  Search for nearest old eigenvalue
	//
	delta [j] = a [m-1] - aold [deflated];
	for (i = deflated; i < mold; i++) 
	  if (fabs (a [j] - aold [i]) < fabs (delta [j])) 
	    delta [j] = a [j] - aold [i];
      }  
    }
    //
    //	Choose vectors to keep (for restart or deflating) and
    //	count accepted eigenvalues (taking into account previous
    //	deflated vectors)
    //
    lfound = lwanted - low;
    hfound = hwanted - high;
    size_t jlow  = deflated + low;
    size_t jhigh = m - high - 1;
    for (j = deflated; j < m; j++) {
      c [j] = 0.0;
      double vale = fabs (b [j]);
      if (vale < eps) {
	//
	//	Mark converged vector for deflation.
	//
	c [j] = -1;
	if (j < jlow)  {
	  //
	  //  	accept deflated low 
	  low--;
	  lfound++;
	}
	if (j > jhigh) {
	  //
	  //	accept deflated high
	  high--;
	  hfound++;
	}
	continue;
      }
      //
      //  Keep lowest eigenvalues (wanted + 6)
      //
      if (j < jlow + lowextra) c [j] = 1.0;
      if (j < jlow)  {
	//
	//	Check if vector was acceptable
	//
	double vald = fabs (delta [j]);
	if (mold == 0) vald = 1.0;	// first iteration: 	don't accept
	if (m == dim)  vald = 0.0;	// no more iterations:	accept
	//
	if (vald < l_zerotolerance && vale < l_zeronorm) lfound++; 
      }
      //
      //  Keep also highest eigenvalues.
      //
      if (j > jhigh - highextra) c [j] = 1.0;
      //
      if (j > jhigh) {
	//
	//  Check if vector was acceptable
	//
	double vald = fabs (delta [j]);
	if (mold == 0) vald = 1.0;	// first iteration: 	don't accept
	if (m == dim)  vald = 0.0;	// no more iterations:	accept
	//
	if (vald < l_zerotolerance && vale < l_zeronorm) hfound++; 
      }
    }
    //
    //	Set the number of inner vectors to keep
    //
    size_t select = high + low;
    //
    //	Check if we need to restart
    //
    if ((lfound + hfound) == (lwanted + hwanted) || beta2 == 0.0) {
      //
      //  We have found all wanted or possible eigenvalues or all wanted 
      //  eigenvectors are deflated or deflating. 
      //  No need to restart (and to compute extra eigenvectors), unmark
      //  extra vectors. 
      //
      for (j = jlow; j < m; j++) c [j] = 0;
      //
      //   Avoid next search.
      //
      select = 0;
      lfound = lwanted;
    }
    //
    //	select some other vector 
    //
    while (select) {
      double emin = 1.00 + beta;
      size_t jmin = m;
      for (j = deflated; j < m; j++) {
	if (c [j]) continue;
	//
	if (fabs (b [j]) < emin) {
	  emin = fabs (b [j]);
	  jmin = j;
	}
      }
      if (jmin == m) select = 0;
      else {
	c [jmin] = 1.0;
	select--;
      }
    }
    //
    //	Compute selected Ritz vectors
    //
    select = 0;
    for (k = deflated; k < m; k++) {
      if (c [k] == 0.0) continue;
      //
      select++;
      //
      w = vector [m + select];
      for (i = csize; i < l_size; i++) w [i] = 0.0;
      //
      for (j = deflated; j < m; j++) {
	v = vector [j];
	eta = y [j - deflated + (k - deflated) * (m - deflated)];
	for (i = csize; i < l_size; i++) w [i] += eta * v [i];
	release (j);
      }						     
      release (m + select);
    }
    //
    //	Put vectors and values in right places
    //		
    select = 0;
    start  = deflated;
    for (k = deflated; k < m; k++) {
      if (c [k] == 0.0) continue;
      //
      select++;
      //
      //  Find insert position for vectors to keep
      //
      size_t insert = start;
      if (c [k] < 0.0) {
	//
	//  Add to deflated vectors and look for insertion position
	//
	for (insert = deflated; insert; insert--) 
	  if (a [insert-1] <= a [k]) break;
	deflated++;
      }
      //
      //  Put vector in place
      //
      for (j = start; j > insert; j--) storage (j) << storage (j-1);
      storage (insert) << storage (m + select);
      //
      //  Put values in place
      //
      double aa = a [k];
       double bb = b [k];
      double dd = delta [k];
      for (j = k; j > insert; j--) {
	a [j]     = a [j-1];
	b [j]     = b [j-1];
	delta [j] = delta [j-1];
      }
      a [insert]     = aa;
      b [insert]     = bb;
      delta [insert] = dd;
      //
      //	Update start position
      start++;
    }
    //
    //	move residual vector in start position 
    //
    storage (start) << storage (m);
    //
    //	Remove unused vectors to avoid meaningless swap reads
    //
    for (k = start + 1; k < m + select + 1; k++) remove (k);
    //
    //	Update eigenvalues and residuals for next step and in case of exit
    //	from while (tolerances were reordered above)
    //
    mold = m;
    for (j = 0; j < m; j++) {
      aold     [j] = a [j];
      residual [j] = b [j];
    }
    //
    //	Check convergengy
    //
    if ((lfound + hfound) == (lwanted + hwanted) || beta2 == 0.0) break;
    //
    //	Verify to have room for next iteration
    //
    kmax = l_steps + deflated;
    if (kmax > 2 * l_steps) kmax = 2 * l_steps;
    if (kmax - start < 2) break;
    //
  } // while (iterations++ ... 
  //
  //	Reorder eigenvalues and wanted eigenvectors
  //
  //	The deflated vectors and eigenvalues are ordered, but the while 
  //	loop may break with some eigenvectors (from deflated to 
  //	deflated + low) not in position, so we must take care of them.
  //	Note that the true deflated vectors number of last iteration 
  //	is recorded in variable 'deflate', not 'deflated'.
  //
  for (k = deflated; k < m; k++) {
    for (j = k; j; j--) if (aold [j-1] <= a [k]) break;
    if (j < k) {
      double aa = aold     [k];
      double dd = delta    [k];
      double bb = residual [k];
      for (i = k; i > j; i--) {
	aold     [i] = aold     [i-1];
	delta    [i] = delta    [i-1];
	residual [i] = residual [i-1];
      }
      aold     [j] = aa;
      delta    [j] = dd;
      residual [j] = bb;
      if (k < deflated + low) {
	//
	//	Move k far away
	//
	storage (m) << storage (k);
	//
	//	Make room in j
	//
	for (i = k; i > j; i--) storage (i) << storage (i-1);
	storage (i) << storage (m);
      }
    }
  }
  //
  //	Remove unwanted vectors
  //
  for (k = l_found; k < 3 * l_steps + 1; k++) remove (k); 
  //
  tm = timecpu () - tm;
  stringstream outp;
  outp << "Found " << lfound << "/" << lwanted; 
  if (hwanted) outp << "+" << hfound << "/" << hwanted;
  outp << " vals " << applyes << "/" << iterations << " its " 
       << locals << "+" << globals << "+" << randoms << " orts ";
  if (deflate) outp << deflate << " deflated ";
  cout << setw (60) << left << outp .str () 
       << " Cpu " << timecpustr (ta) << "/" << timecpustr (tm) << endl;
  //
  //	Set actual found and last iteration steps on return.
  //
  l_found = lfound;
  l_steps = m;
  return l_found;
}
//
//============================================================================
Superaction::Superaction ()
  : sa_bipoli	(),
    sa_biaction	(),
    sa_inner	()
{ 
  //
  //	Default ctor.
  //
}
//
//____________________________________________________________________________
Superaction::Superaction (const Superaction & other)
  : sa_bipoli	(other .sa_bipoli),
    sa_biaction	(other .sa_biaction),
    sa_inner	(other .sa_inner)
{ 
  //
  //	Copy ctor.
  //
}
//
//____________________________________________________________________________
Superaction::Superaction (const Apoli & poli, 
			  Block & system, Block & universe, bool releasing)  
  : sa_bipoli	(),
    sa_biaction	(),
    sa_inner	()
{ 
  //
  //	Transform the formal polinomial poli into a sum of tensor 
  //	products of operators acting on system and universe Block's
  //
  //	Express polinomial poli as a sum of products of polinomials
  //	of Action's relative to system and universe 
  //	storing both formal expressions and true Action
  //	
  size_t m, n, ml, mr, nl, nr;
  //
  //	Load sa_bipoli with splitted polinomial
  //
  long split = system .sites ();
  for (m = 0; m < poli .size (); m++) {
    Bipoli factor;
    const Amono & mono     = poli [m];
    factor .bi_coeff = mono .am_coeff;
    for (n = 0; n < mono .order (); n++) {
      Afactor af (mono [n]);
      if (af .af_st < split)  	
	factor .bi_lft *= af;
      else {
	af .af_st -= split;
	factor .bi_rgt *= af;
      }
    }
    sa_bipoli .push_back (factor);
  } 
  //
  while (true) {
    ml = mr = nl = nr = 0;
    //
    //	Find highest rank common factor
    //
    for (m = 0 ; m < sa_bipoli .size (); m++) {
      if ((sa_bipoli [m] .bi_lft .order () == 0) &&
	  (sa_bipoli [m] .bi_rgt .order () == 0))	continue;
      size_t cl = 0;
      size_t cr = 0;
      //	
      //  compute number of times lft and right factor appear 
      //  in next terms.
      //
      for (n = m+1; n < sa_bipoli .size (); n++) {
	if (sa_bipoli [n] .bi_lft == sa_bipoli [m] .bi_lft) cl++;
	if (sa_bipoli [n] .bi_rgt == sa_bipoli [m] .bi_rgt) cr++;
      }
      if (cl > nl) {
	ml = m; 	
	nl = cl;
      }
      if (cr > nr) {
	mr = m;
	nr = cr;
      }
    }
    //
    //	if no multiple common factors exit from while
    //
    if (nl == 0 && nr == 0) break;
    //
    //	Collect highest rank common factor 
    //
    if (nl < nr) {
      //
      //  common right factor
      //
      m = mr;
      nl = 0;
    }
    else {	
      //
      //  common left factor
      //
      m = ml;
      nr = 0;
    }
    if (nr) {
      //
      //  collect left factors
      //
      sa_bipoli [m] .bi_lft *= sa_bipoli [m] .bi_coeff;
      sa_bipoli [m] .bi_coeff = 1.0;
      for (n = m+1; n < sa_bipoli .size (); n++)
	if (sa_bipoli [n] .bi_rgt == sa_bipoli [m] .bi_rgt) {
	  sa_bipoli [n] .bi_lft *= sa_bipoli [n] .bi_coeff;
	  sa_bipoli [m] .bi_lft += sa_bipoli [n] .bi_lft;
	  sa_bipoli .erase (sa_bipoli .begin () + n);
	  n--;
	}
      //
      //  normalize
      //
      if (sa_bipoli [m] .bi_lft .size () == 0) 
	sa_bipoli .erase (sa_bipoli .begin () + m);
      else if (sa_bipoli [m] .bi_lft [0] .am_coeff != 1.0) { 
	complex<double> coeff = sa_bipoli [m] .bi_lft [0] .am_coeff;
	sa_bipoli [m] .bi_coeff *= coeff;
	sa_bipoli [m] .bi_lft   /= coeff;
      }
    }
    if (nl) {
      //
      //  collect right factors
      //
      sa_bipoli [m] .bi_rgt *= sa_bipoli [m] .bi_coeff;
      sa_bipoli [m] .bi_coeff = 1.0;
      for (n = m+1; n < sa_bipoli .size (); n++)
	if (sa_bipoli [n] .bi_lft == sa_bipoli [m] .bi_lft) {
	  sa_bipoli [n] .bi_rgt *= sa_bipoli [n] .bi_coeff;
	  sa_bipoli [m] .bi_rgt += sa_bipoli [n] .bi_rgt;
	  sa_bipoli .erase (sa_bipoli .begin () + n);
	  n--;
	}
      // 
      // normalize
      //
      if (sa_bipoli [m] .bi_rgt .size () == 0)
	sa_bipoli .erase (sa_bipoli .begin () + m);
      else if (sa_bipoli [m] .bi_rgt [0] .am_coeff != 1.0) {
	complex<double> coeff = sa_bipoli [m] .bi_rgt [0] .am_coeff;
	sa_bipoli [m] .bi_coeff *= coeff;
	sa_bipoli [m] .bi_rgt   /= coeff;
      }
    }
  }
  for (m = 0; m < size (); m++) {
    if ((sa_bipoli [m] .bi_lft .order ()  > 0) &&
	(sa_bipoli [m] .bi_rgt .order () == 0)) {
      sa_bipoli [m] .bi_lft *= sa_bipoli [m] .bi_coeff;
      sa_bipoli [m] .bi_coeff = 1.0;
    }
    if ((sa_bipoli [m] .bi_rgt .order ()  > 0) &&
	(sa_bipoli [m] .bi_lft .order () == 0)) {
      sa_bipoli [m] .bi_rgt *= sa_bipoli [m] .bi_coeff;
      sa_bipoli [m] .bi_coeff = 1.0;
    }
  }
  //
  //	Special collect for unitary factors	
  //
  n = nr = nl = sa_bipoli .size ();
  for (m = 0; m < sa_bipoli .size (); m++) {
    if ((sa_bipoli [m] .bi_lft .order () == 0) &&
	(sa_bipoli [m] .bi_rgt .order () == 0)) n = m; 
    if ((sa_bipoli [m] .bi_lft .order () >  0) &&
	(sa_bipoli [m] .bi_rgt .order () == 0)) nr = m; 
    if ((sa_bipoli [m] .bi_lft .order () == 0) &&
	(sa_bipoli [m] .bi_rgt .order () >  0)) nl = m; 
  }
  //
  if ((n < sa_bipoli .size ()) && 
      (! ((nl == sa_bipoli .size ()) && (nr == sa_bipoli .size ()))) ) {
    //
    //	Distribute pure scalar between left and right scalar term  
    //
     complex<double> coeff = sa_bipoli [n] .bi_coeff *
      sa_bipoli [n] .bi_lft [0] .am_coeff * 
      sa_bipoli [n] .bi_rgt [0] .am_coeff;
    if ((nl < sa_bipoli .size ()) && (nr < sa_bipoli .size ())) 
      coeff *= 0.5;
    if (nl < sa_bipoli .size ()) {
      //
      //	Add scalar to right true polinomial
      //
      sa_bipoli [nl] .bi_rgt *= sa_bipoli [nl] .bi_coeff;
      sa_bipoli [nl] .bi_coeff = 1.0;
      sa_bipoli [nl] .bi_rgt += coeff;
    }
    if (nr < sa_bipoli .size ()) {
      //
      //	Add scalar to left true polinomial
      //
      sa_bipoli [nr] .bi_lft *= sa_bipoli [nr] .bi_coeff;
      sa_bipoli [nr] .bi_coeff = 1.0;
      sa_bipoli [nr] .bi_lft += coeff;
    }
    sa_bipoli .erase (sa_bipoli .begin () + n);
  }
  /*
  //
  //	Optimizazion of basis
  //
  for (m = 0; m < sa_bipoli .size (); m++) 
    if ((sa_bipoli [m] .bi_lft .order () > 0) &&
	(sa_bipoli [m] .bi_rgt .order () > 0)) {
      system   .optimize (sa_bipoli [m] .bi_lft);
      universe .optimize (sa_bipoli [m] .bi_rgt);
  }
  */
  //
  //	Load operators
  //
  for (n = 0; n < sa_bipoli .size (); n++) {
    sa_biaction .push_back (Biaction ());
    Action & lop = sa_biaction [n] .ba_lftop;
    Action & rop = sa_biaction [n] .ba_rgtop;
    //
    //	Build Action's
    //
    system   .action (sa_bipoli [n] .bi_lft, lop);
    universe .action (sa_bipoli [n] .bi_rgt, rop);
    if (lop .isscalar () && rop .isscalar ()) {
      lop *= sa_bipoli [n] .bi_coeff;
      sa_bipoli [n] .bi_lft   *= sa_bipoli [n] .bi_coeff;
      sa_bipoli [n] .bi_coeff = 1.0;
    }
    else if (lop .isscalar ()) {
      rop *= sa_bipoli [n] .bi_coeff;
      sa_bipoli [n] .bi_rgt   *= sa_bipoli [n] .bi_coeff;
      sa_bipoli [n] .bi_coeff = 1.0;
    }
    else {
      lop *= sa_bipoli [n] .bi_coeff;
      sa_bipoli [n] .bi_lft   *= sa_bipoli [n] .bi_coeff;
      sa_bipoli [n] .bi_coeff = 1.0;
    }
    if (! lop .isscalar ()) lop .normalize ();
    if (! rop .isscalar ()) rop .normalize ();
    rop .transpose ();
    //lop .sparse ();
    //rop .sparse ();
    //
    if (releasing) lop .release ();
    if (releasing) rop .release ();
  }
}
//
//____________________________________________________________________________
Superaction::~Superaction ()
{ 
  //
  //	Default dtor.
  //
}
//
//____________________________________________________________________________
bool Superaction::iscomplex () const
{
  //
  //	Check for some imaginary component
  //
  for (size_t n = 0; n < size (); n++) 
    if (sa_biaction [n] .ba_lftop .iscomplex () ||
	sa_biaction [n] .ba_rgtop .iscomplex ()) return true;
  return false;
}
//
//____________________________________________________________________________
Superaction & Superaction::operator = (const Superaction & other)
{
  //
  //	copy assignment
  //
  if (this != & other) {
    sa_bipoli    = other .sa_bipoli;
    sa_biaction  = other .sa_biaction;
    sa_inner    << other .sa_inner;	// grab inner storage area
  }
  return *this;
}
//
//____________________________________________________________________________
void Superaction::show (const string & s) const
{
  //
  //	Prints a representation of Superaction 
  //
  if (s .size ()) cout << s << endl;
  for (size_t m = 0; m < size (); m++)  {
    if (coefficient (m) == 1.0) ;
    else {
      cout << setw (10) << right << showpos << setprecision (4)  
	 << coefficient (m)  << " * " << setprecision (0);
      cout .unsetf (ios_base::showpos);
    }
    if (m) cout << "+ "; 
    cout << "(" << endl;
    leftpolinomial (m) .show ("", 8);
    sa_biaction [m] .ba_lftop .show ();
    cout << ") * ( transpose " << endl;
    rightpolinomial (m) .show ("", 8);
    sa_biaction [m] .ba_rgtop .show ();
    cout << ") " << endl;
  }
}
//
//____________________________________________________________________________
void Superaction::storage (size_t n)
{
  //
  //	Allocate inner storage area to n doubles
  //
  sa_inner .storage (n * sizeof (double));
}
//
//============================================================================
size_t biaction (Action & result, 
		 Superaction & superaction, const Action & state,
		 bool releasing)
{
  //
  //	Setup superaction and result for fast application to state.
  //
  //	Returns the maximum size of needed intermediate memory
  //	which is allocated inside superaction
  //	
  //	This routine must be called once before application to a 
  //	new state (or when the structure of state is changed).
  //
  size_t maxsize = 0;
  size_t tocopy =  0;
  long   rstat = 0;
  //
  //	Initialize result as a null scalar Action if not defined
  //
  Action apply = Action (state);
  if (result .storage () == 0) 
    result = Action (state .range (), state .domain ());
  if (result .storage () == state .storage ()) tocopy = 1;
  //
  //	Loop on tensor products: compute Action's applications  
  //	binary descriptions and result's structure (merging 
  //	for the sum of all tensor products).
  //
  //	No memory allocation for the result, required allocation size
  //	is only computed.
  //
  for (size_t n = 0; n < superaction .size (); n++) {
    Biaction & ba = superaction [n];
    Action inner, left;
    long stat = apply .statistic ();
    inner << apply;
    ba .ba_rgtab .clear ();
    ba .ba_lftab .clear ();
    if (ba .ba_rgtop .storage ()) {
      inner = Action (apply, ba .ba_rgtop);
      stat = inner .statistic ();
      if (inner .blocks ()) inner .scalar (1.0);
      binary (ba .ba_rgtab, inner, apply, ba .ba_rgtop);
      if (releasing) ba .ba_rgtop .release ();
    }
    if (inner .size () > maxsize) maxsize = inner .size ();
    ba .ba_size = inner .size ();
    if (ba .ba_lftop .storage ()) {
      left = Action (ba .ba_lftop, inner);
      stat = left .statistic ();
      if (left .blocks ()) left .scalar (1.0);
      result .merge (left, tocopy);
      binary (ba .ba_lftab, result, ba .ba_lftop, inner);
      if (releasing) ba .ba_lftop .release ();
    }
    else result. merge (inner);
    //
    //	Update global statistic
    //
    if (n == 0) rstat  = stat;
    if (rstat != stat) rstat = 0;
  }
  //
  result      .scalar    (1.0);
  result      .statistic (rstat);
  superaction .storage   (maxsize);
  return maxsize;
}
//
//____________________________________________________________________________
void biapply (double * mr, Superaction & superaction, double * ms,
	      bool releasing) 
{
  //
  //	Fast apply a superblock operator (assuming that	all is well 
  //	defined) using allocated memory pointers.
  //
  double * mi = superaction .storage ();
  for (size_t n = 0; n < superaction .size (); n++) {
    Biaction & ba = superaction [n];
    double * minner  = ms;
    if (ba .ba_rgtop .storage ()) {
      memset (mi, 0, ba .ba_size * sizeof (double));
      multiply (mi, ms, 0, 
		ba .ba_rgtop .storage (), ba .ba_rgtop .sparsed (), 
		ba .ba_rgtab, true);
      if (releasing) ba .ba_rgtop .release ();
      minner = mi;
    }
    if (ba .ba_lftop .storage ()) {
      multiply (mr, ba .ba_lftop .storage (), ba .ba_lftop .sparsed (),
		minner, 0, ba .ba_lftab);
      if (releasing) ba .ba_lftop .release ();
    }
    else 
      for (size_t m = 0; m < ba .ba_size; m++) mr [m] += minner [m];
  }
}
//
//____________________________________________________________________________
void biapply (Action & result, 
	      Superaction & superaction, const Action & state, bool releasing)
{
  //
  //	Apply a superblock operator superaction represented by a sum
  //	of tensor products to state:
  //	
  //	 	result = superaction * state
  //
  size_t nbi;
  //
  //	Setup superaction and result 
  //
  nbi = biaction (result, superaction, state, releasing);
  result .storage (result .size ());
  if (result .blocks ()) result .scalar (1.0);
  double * ms = state  .storage ();
  double * mr = result .storage ();
  biapply (mr, superaction, ms, releasing);
}
//
//============================================================================

