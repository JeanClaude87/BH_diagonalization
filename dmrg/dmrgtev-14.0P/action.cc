//	action.cc			(C) 2009 Fabio Ortolani	fbo 090716
//	==================================================================
#include "action.hh"
#include "numerical.hh"
#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>
//
//============================================================================
//
//	Initial memory of defined Action's must contain al least a 
//	complex number. csize is the size of a complex number in units 
//	of double size.
//
//	This constant can be accessed by function minasize () from 
//      other program units
//	
static const size_t csize = 
  (sizeof (complex<double>) + sizeof (double) -1) / sizeof (double);
//
//      sparse_num / sparse_den gives the fraction to accept sparsing 
//
static const size_t sparse_num = 5;
static const size_t sparse_den = 8;
extern vector<double>  renyientropy;		   	   // [superaction.cc]
extern vector<double>  renyiarg;		   	   // [superaction.cc]
double * autovalori;
//
void fbmm (const long, const long, const long, const long, const long, double,
	   const double *, const long, const double *, const long,
	   double *, const long);
//============================================================================
Aarray::Aarray  ()
  : aa_action	()
{
  //
  //	Default constructor: empty array
  //
}
//
//____________________________________________________________________________
Aarray::Aarray  (const Aarray & other)
  : aa_action	(other .aa_action)
{ 
  //
  //	Copy constructor: transfer (not copy) of other array into this.
  //
  //	We need to reset pointers of other Aarray (only one pointer
  //	to a given action must be active at any time).
  //	With a dirty trick to fake the compiler (STL methods assume 
  //	other a const Aarray).
  //
  vector<Action *> * va = (vector<Action *> *) & (other .aa_action);
  for (size_t i = 0; i < other .size (); i++) (*va) [i] = 0;
}
//
//____________________________________________________________________________
Aarray::~Aarray ()
{
  //
  //	Default destructor: destroy pointed Action's (created via new
  //	operator)
  //
  for (size_t i = 0; i < size (); i++) 
    if (aa_action [i]) delete aa_action [i];
}
//
//____________________________________________________________________________
const Action & Aarray::operator [] (size_t index) const
{
  return *(aa_action [index]);
}
//
//____________________________________________________________________________
Action & Aarray::operator [] (size_t index)
{	
  //
  //	Array like access with dynamical growing (no check of index range).
  //	Gives access to the pointed content, not to the pointer itself.
  //
  if (size () <= index) aa_action .resize (index + 1, (Action *) 0);
  //
  //	If null pointer create a new content
  //
  if (aa_action [index] == 0) aa_action [index] = new Action;
  if (aa_action [index] == 0) {
    cout << "No room for new Action in array!" << endl;
    exit (0);
  }
  return *(aa_action [index]);
} 
//
//____________________________________________________________________________
Aarray & Aarray::operator = (const Aarray & other) 
{
  //
  //	copy assignment: effective transfer, not copy.
  //
  //	Old content (if any) of this array must be destroyed if not 
  //	used in new array (only one pointer to a given Action must
  //	be active in memory).
  //	
  //	To make a true copy an explicit loop on other array must be 
  //	done with an assignment of content (generating, if needed, 
  //	new pointers with access operator).
  //
  size_t i, j;
  for (i = 0; i < size (); i++) {
    if (aa_action [i]) {
      //
      //  Check if pointer is used in other array
      //
      for (j = 0; j < other .size (); j++) 
	if (aa_action [i] == other .aa_action [j]) break;
      //
      // if not used delete pointer content
      //
      if (j >= other .size ()) 	delete aa_action [i];
    }
    //
    //  To be sure and avoid further deletions of content
    //
    aa_action [i] = 0;
  }
  //
  //	Now we can copy the new array into the old one
  //
  aa_action = other .aa_action;
  //
  //	We need to reset other's pointers wth a dirty trick to fake
  //	the compiler.
  //
  vector<Action *> * va = (vector<Action *> *) & (other .aa_action);
  for (size_t i = 0; i < size (); i++) {
    //aa_action [i] = other .aa_action [i];
    (*va) [i] = 0;
  }
  return *this;
}
//
//============================================================================
Abinary::Abinary ()
  : ab_lftro	(0),
    ab_lftio	(0),
    ab_rgtro	(0),
    ab_rgtio	(0),
    ab_resro	(0),
    ab_resio	(0),
    ab_lftsz	(0),
    ab_rgtsz	(0),
    ab_intsz	(0),
    ab_lftsp	(0),
    ab_rgtsp	(0)
{ 
  //
  //	Default ctor:	all zero.
  //
}
//
//____________________________________________________________________________
Abinary::~Abinary ()
{ 
  //
  //	Default dtor.
  //
}
//
//============================================================================
Ablock::Ablock ()
  : ab_roffset	(0),
    ab_ioffset  (0),
    ab_domain	(0),
    ab_range	(0),
    ab_cindex	(0),
    ab_rindex	(0)
{ 
  //
  //	Default ctor:	all zero
  //
}
//
//____________________________________________________________________________
Ablock::~Ablock ()
{ 
  //
  //	Default dtor.
  //
}
//
//============================================================================
Action::Action  ()
  : a_range	(),
    a_domain	(), 
    a_size	(0),
    a_blocks	(0),
    a_statistic	(0),
    a_storage	(),
    a_block	(),
    a_sparsed	()
{ 
  //
  //	Default ctor. 
  //	Creates an undefined action, all fields are null or empty
  //
}
//
//____________________________________________________________________________
Action::Action  (const Action & other)
  : a_range	(other .a_range),
    a_domain	(other .a_domain), 
    a_size	(other .a_size),
    a_blocks	(other .a_blocks),
    a_statistic (other .a_statistic),
    a_storage	(),
    a_block	(),
    a_sparsed	()
{ 
  //
  //	Copy ctor: Memory and structure are copied from other Action
  //
  if (other .storage ()) {
    size_t osize = other .a_storage .size ();
    a_storage .storage (osize);
    memcpy (storage (), other .storage (), osize);
  }
  else a_storage .storage (0);
  if (other .block ()) {
    size_t osize = other .a_block .size ();
    a_block .storage (osize);
    memcpy (block (), other .block (), osize);
  }
  else a_block .storage (0);
  if (other .sparsed ()) {
    size_t osize = other .a_sparsed .size ();
    a_sparsed .storage (osize);
    memcpy (sparsed (), other .sparsed (), osize);
  }
  else a_sparsed .storage (0);
}
//
//____________________________________________________________________________
Action::Action (const Action & a, const Action & b) 
  : a_range	(a. a_range),
    a_domain	(b .a_domain),
    a_size	(csize),
    a_blocks	(0),
    a_statistic	(a .a_statistic * b .a_statistic),
    a_storage	(csize * sizeof (double)),
    a_block	(),
    a_sparsed	()
{
  //
  //	Special constructor.
  //
  //	Creates an Action and setup block structure for the matrix 
  //	product of two given Action's. (No memory is allocated, 
  //   	required size for non null blocks is computed and stored 
  //	in the a_size field).
  //
  //	The two Action factors must be in normalized form.
  //
  //	If one of the two operands is null or both operands are scalar, 
  //	the result is a null scalar matrix, else it is a non scalar
  //	null matrix (equivalent to a scalar matrix).
  //
  if (! (a .domain () == b .range ())) {
    cout << "Incompatible Action product! "
	 << a .width () << " .ne. " << b .height () << endl;
    a .show ("a");
    b .show ("b");
    exit (0);
  }
  if (! (a .isnormal () && b .isnormal ())) {
    cout << "Product of non normal Action's!" << endl;
    exit (0);
  }
  scalar (0.0);
  //
  //	Scalar * scalar --> scalar; 
  //
  if (a .isscalar () && b .isscalar ()) return;
  //
  //	zero factor	--> zero (scalar)
  //
  if (a .isnull () || b .isnull ()) return;
  //
  //	Setup Ablock structure and compute required memory	
  //
  size_t na, nb, nc, nna, nnb;
  nna = a .blocks ();
  nnb = b .blocks ();
  Ablock * ba = a .block ();
  Ablock * bb = b. block ();
  //
  //	New blocks structure
  //
  vector<Ablock> blk;
  //
  if (a .isscalar ()) {
    //
    //	Scalar * matrix --> matrix
    //
    double zr = a .scalar () .real ();
    double zi = a .scalar () .imag ();
    for (nb = 0; nb < nnb; nb++) {
      Ablock bc = bb [nb];
      size_t mem = height (bc .ab_range) * width (bc .ab_domain);
      bc .ab_roffset = 0;
      bc .ab_ioffset = 0;
      if ((zr && bb [nb] .ab_roffset) || (zi && bb [nb] .ab_ioffset)) {
	//
	//	Real part
	//
	bc .ab_roffset = a_size;
	a_size += mem;
      }
      if ((zr && bb [nb] .ab_ioffset) || (zi && bb [nb] .ab_roffset)) {
	//
	//	Imaginary part
	//
	bc .ab_ioffset = a_size;
	a_size += mem;
      }
      if (bc .ab_roffset || bc .ab_ioffset) blk  .push_back (bc);
    }
    nna = nnb = 0;	// (to be sure)
  }
  //
  if (b .isscalar ()) {
    //
    //	Matrix * scalar	--> matrix
    //
    double zr = b .scalar () .real ();
    double zi = b .scalar () .imag ();
    for (na = 0; na < nna; na++) {
      Ablock bc = ba [na];
      size_t mem = height (bc .ab_range) * width (bc .ab_domain);
      bc .ab_roffset = 0;
      bc .ab_ioffset = 0;
      if ((zr && ba [na] .ab_roffset) || (zi && ba [na] .ab_ioffset)) {
	//
	//	Real part
	//
	bc .ab_roffset = a_size;
	a_size += mem;
      }
      if ((zr && ba [na] .ab_ioffset) || (zi && ba [na] .ab_roffset)) {
	//
	//	Imaginary part
	//
	bc .ab_ioffset = a_size;
	a_size += mem;
      }
      if (bc .ab_roffset || bc .ab_ioffset ) blk  .push_back (bc);
    }
    nna = nnb = 0;	// (to be sure)
  }
  //
  //	Look for possible non null block products
  //
  for (nb = 0; nb < nnb; nb++) {
    if ((bb [nb] .ab_roffset == 0) && (bb [nb] .ab_ioffset == 0)) continue;
    for (na = 0; na < nna; na++) {
      if ((ba [na] .ab_roffset == 0) && (ba [na] .ab_ioffset == 0)) continue;
      if (ba  [na] .ab_domain != bb [nb] .ab_range) continue;
      //
      //  Check block definition in the list
      //
      for (nc = 0; nc < blk .size (); nc++) 
	if ((blk [nc] .ab_range  == ba [na] .ab_range) &&
	    (blk [nc] .ab_domain == bb [nb] .ab_domain)) break;
      //
      //  add new product block
      //
      if (nc == blk .size ()) blk .push_back (Ablock ());
      //
      Ablock & bc = blk [nc];
      bc .ab_domain = bb [nb] .ab_domain;
      bc .ab_range  = ba [na] .ab_range;
      size_t ms = width (bc .ab_domain) * height (bc .ab_range);
      if ((bc .ab_roffset == 0) && 
	  ((ba [na] .ab_roffset && bb [nb] .ab_roffset) ||
	   (ba [na] .ab_ioffset && bb [nb] .ab_ioffset))) {
	bc .ab_roffset = a_size;
	a_size += ms;
      }
      if ((bc .ab_ioffset == 0) && 
	  ((ba [na] .ab_roffset && bb [nb] .ab_ioffset) ||
	   (ba [na] .ab_ioffset && bb [nb] .ab_roffset))) {
	bc .ab_ioffset = a_size;
	a_size += ms;
      }
    }
  }
  //
  //	Fill blocks structure
  //
  blocks (blk .size ());
  block  (blk .size ());
  Ablock * bl = block ();
  for (nc = 0; nc < blk .size (); nc++) bl [nc] = blk [nc];
}
//
//____________________________________________________________________________
Action::Action (const Space & range, const Space & domain)
  : a_range	(range),
    a_domain	(domain),
    a_size	(csize),
    a_blocks	(0),
    a_statistic	(0),
    a_storage	(csize * sizeof (double)),
    a_block	(),
    a_sparsed	()
{ 
  //
  //	Special ctor.
  //
  //	Build a null Action from domain to range Space's
  //
  scalar (0.0);	// null Action
}
//
//____________________________________________________________________________
Action::Action (const Space & space, complex<double> factor)
  : a_range	(space),
    a_domain	(space),
    a_size	(csize),
    a_blocks	(0),
    a_statistic	(1),
    a_storage	(csize * sizeof (double)),
    a_block	(),
    a_sparsed	()
{
  //
  //	Special ctor.
  //
  //	Build a constant times the identity in space
  //
  scalar (factor);
}
//
//____________________________________________________________________________
Action::Action (const Qspace &  range,  const Qspace & domain,
		const Qnumber & select, bool cmplx, const double qa)
  : a_range	(range),
    a_domain	(domain),
    a_size	(csize),
    a_blocks	(0),
    a_statistic (0),
    a_storage	(csize * sizeof (double)),
    a_block	(),
    a_sparsed	()
{
  //
  //	Special ctor.
  //
  //	Build a full structured Action from domain to range spaces
  //	and setup block entries such that the sum of corresponding
  //	domain and range quantum numbers are equal to selected quanum 
  //	number.	No memory is allocated, required size is only computed.
  //	If cmplx is true imaginary offsets are set also.
  // 
  //	This is the structure of an Action representing a state in 
  //	the subspace of tensor product of range and domain with select
  //	quantum numbers.
  //
  vector<Ablock> blk;
  //
  //cout << "select qnumber " << select .str ()<< endl;
  for (size_t db = 1; db < domain .subspaces (); db++) {
    size_t w  = domain .states    (db);
    long   ds = domain .statistic (db);
    //	cout << domain.number (db) .str () << endl;
    for (size_t rb = 1; rb < range .subspaces (); rb++) {
      if (range .rgtlnk (rb) != domain .lftlnk (db)) continue;
      size_t h  = range .states    (rb);
      size_t rs = range .statistic (rb);
      // cout << range .number (rb) .str () << " ";
      Qnumber rd = compose (range .number (rb), domain .number (db), qa);
      //cout << rd .str () << endl;
      if (select == rd) {
	//cout << range .number (rb) .str () << " "
	//   << domain.number (db) .str () << endl;
	Ablock b;
	b .ab_range   = rb;
	b .ab_domain  = db;
	b .ab_roffset = a_size;
	a_size += h * w;
	if (cmplx) {
	  b .ab_ioffset = a_size;
	  a_size += h * w;
	}
	//
	//	Update statistic
	//
	long s = rs * ds;
	if (blk .size () == 0) a_statistic = s;
	if (s != a_statistic) a_statistic = 0;
	blk .push_back (b);
      }
    }
  }
  //
  //	Allocate and fill blocks
  //
  //cout << "selected " << a_size << endl;
  blocks (blk .size ());
  block  (blk .size ());
  Ablock * selected = block ();
  for (size_t ns = 0; ns < a_blocks; ns++) selected [ns] = blk [ns];
}
//
//____________________________________________________________________________
Action::~Action ()
{
  //
  //	Default dtor.
  //
}
//
//____________________________________________________________________________
void Action::add (complex<double> c, const Action & other)
{
  //
  //	Adds Action other times scalar c
  //
  if (this == & other) {
    //
    //	Matrix part is the same
    //
    scalar (scalar () + c);
    return;
  }
  Action adding (other);
  adding *= c;
  *this += adding;
}
//____________________________________________________________________________
void Action::clean    ()
{
  //
  //  Check for effective null blocks putting corresponding 
  //  pointers to zero
  //
  vector<Ablock> blk;
  Ablock * b = block   ();
  double * m = storage ();
  //
  //	Counting not null elements for sparsing allocation
  //
  double eps = machine_precision ();
  size_t update = 0;
  size_t n;
  for (n = 0; n < blocks (); n++) {
    double * mr;
    double * mi;
    size_t h = height (b [n] .ab_range);
    size_t w = width  (b [n] .ab_domain);
    mr = 0;
    if (b [n] .ab_roffset) mr = m + b [n] .ab_roffset;
    mi = 0;
    if (b [n] .ab_ioffset) mi = m + b [n] .ab_ioffset;
    size_t nonzeror = 0;
    size_t nonzeroi = 0;
    for (size_t i = 0; i < h * w; i++) {
      double vr = 0.0;
      double vi = 0.0;
      if (mr) {
	if (fabs (mr [i]) < eps) mr [i] = 0.0;
	else vr = mr [i];
      }
      if (mi) {
	if (fabs (mi [i]) < eps) mi [i] = 0.0;
	else vi = mi [i];
      }
      if (vr) nonzeror++;
      if (vi) nonzeroi++;
    }
    if (nonzeror == 0) b [n] .ab_roffset = 0;
    if (nonzeroi == 0) b [n] .ab_ioffset = 0;
    if (nonzeror + nonzeroi == 0) update = 1;
    else blk .push_back (b[n]);
  }
  if (update) {
    blocks (blk .size ());
    block  (blk .size ());
    b = block ();
    for (n = 0; n < blocks (); n++) b [n] = blk [n];
  }
}
//
//____________________________________________________________________________
void Action::compress (complex<double> * full)
{
  //
  //	Fill block matrix entries from a full matrix. 
  //
  size_t ib, jb, ie, je, i, j, len, h, w, io, jo, rstore, istore, nonzero;
  ie  = a_range  .subspaces ();
  je  = a_domain .subspaces ();
  len = height ();
  vector<Ablock> blk;
  //
  //	Find non null blocks and compute memory size
  //
  a_size = csize;
  for (jb = 1; jb < je; jb++) {
    w  = width (jb);
    jo = a_domain .offset (jb);
    for (ib = 1; ib < ie; ib++) {
      h  = height (ib);
      io = a_range .offset (ib);
      rstore = istore = 0;
      for (j = jo; j < jo+w; j++)
	for (i = io; i < io+h; i++) {
	  complex<double> z = full [i + j*len];
	  if (z .real ()) rstore = 1;
	  if (z .imag ()) istore = 1;
	}
      if (rstore || istore) {
	Ablock bk;
	bk .ab_domain = jb;
	bk .ab_range  = ib;
	bk .ab_roffset = 0;
	bk .ab_ioffset = 0;
	if (rstore) {
	  bk .ab_roffset = a_size;
	  a_size += h * w;
	}
	if (istore) {
	  bk .ab_ioffset = a_size;
	  a_size += h * w;
	}
	blk .push_back (bk);
      }
    }
  }
  //
  //	Allocate memory and fill it
  //	
  storage   (a_size);
  double * a = storage ();
  blocks    (blk .size ());
  block     (blk .size ());
  Ablock * b = block   ();
  a_statistic = 0;
  long st = 0;
  scalar (0.0);
  for (nonzero = 0; nonzero < a_blocks; nonzero++) {
    b [nonzero] = blk [nonzero];
    ib = b [nonzero] .ab_range;
    jb = b [nonzero] .ab_domain;
    w  = width  (jb);
    h  = height (ib);
    jo = a_domain .offset (jb);
    io = a_range  .offset (ib);    
    if ((rstore = b [nonzero] .ab_roffset)) 
      for (j = jo; j < jo+w; j++) 
	for (i = io; i < io+h; i++) 
	  a [rstore++] = full [i + j*len] .real ();
    if ((istore = b [nonzero] .ab_ioffset))
      for (j = jo; j < jo+w; j++) 
	for (i = io; i < io+h; i++) 
	  a [istore++] = full [i + j*len] .imag ();
    //
    //	Compute Fermi-Bose statistic
    //
    st = a_range .statistic (ib) * a_domain .statistic (jb);
    if (nonzero == 0) a_statistic = st;
    else if (st != a_statistic) a_statistic = 0;
  }
  //
  //	if no blocks it is a null matrix
  //
  if (a_blocks) scalar (1.0);
  else		scalar (0.0);
}
//
//____________________________________________________________________________
void Action::dagger ()
{
  //
  //	Hermitian conjugate transformation
  //
  double * m = storage ();
  if (m == 0) {
    cout << "Adjoint of undefined Action!" << endl;
    exit (0);
  }
  //
  //	define new memory area
  //
  Storage  t    (a_size * sizeof (double));
  double * mt = (double *) t  .storage ();
  ((complex<double> *) mt) [0] = conj (((complex<double> *) m) [0]);
  //
  //	All done for a scalar
  //
  if (a_size == csize) return;
  //
  //	define new blocks
  //
  Storage  tb (blocks () * sizeof (Ablock));
  Ablock * bt = (Ablock *) tb .storage ();
  //
  Ablock * b = block ();
  for (size_t nb = 0; nb < blocks (); nb++) {
    bt [nb]     = b [nb];
    size_t roff = b [nb] .ab_roffset;
    size_t ioff = b [nb] .ab_ioffset;
    size_t is   = b [nb] .ab_range;
    size_t js   = b [nb] .ab_domain;
    size_t w    = width  (js);
    size_t h    = height (is);
    if (roff) 
      for (size_t j = 0; j < w; j++) 
	for (size_t i = 0; i < h; i++) 
	  mt [roff + j + i*w] =  m [roff +i +j*h];
    if (ioff) 
      for (size_t j = 0; j < w; j++) 
	for (size_t i = 0; i < h; i++) 
	  mt [ioff + j + i*w] = -m [ioff +i +j*h];
    bt [nb] .ab_range  = js;
    bt [nb] .ab_domain = is;
  }
  //
  //	Grab new memory and block
  //
  storage (t);
  block   (tb);
  //
  //	Remove any sparsing
  //
  sparsed (0);
  Space s = a_domain;
  a_domain = a_range;
  a_range  = s;
}
//
//____________________________________________________________________________
size_t Action::dimension () const
{
  //
  //	Returns the number of non null (complex) matrix entries
  //	pointed by actual blocks
  //
  Ablock * b = block ();
  size_t dim = 0;
  for (size_t n = 0; n < blocks (); n++) 
    if (b [n] .ab_roffset || b [n] .ab_ioffset) 
      dim += height (b [n] .ab_range) * width (b [n] .ab_domain);
  return dim;
}
//
//____________________________________________________________________________
double Action::entropy () const
{
  //
  //	Returns entropy of this Action considered as a super state
  //
  size_t   i, nb, sub, offset, dimension;
  //
  //	compute range reduced density matrix (tracing on domain) 
  //
  Action s = *this;
  Action density = s;
  s .dagger ();
  density *= s;
  density .clean ();
  //
  Ablock * b  = density .block   ();
  double * m  = density .storage ();
  //
  size_t states = height ();
  size_t ss = (2 * states * sizeof (double));
  Storage  mem (ss);
  double * eigen = (double *) mem .storage ();
  autovalori = new double  [states];
  double * offd  = eigen + states;
  for (i = 0; i < 2 * states; i++) eigen [i] = 0.0;
  //
  for (nb = 0; nb < density .blocks (); nb++) {
    //	
    //	A simple check if density is block diagonal
    //
    if (b [nb] .ab_domain != b [nb] .ab_range) {
      cout << "Not block diagonal density matrix!" << endl;
      exit (0);
    }
    sub = b [nb] .ab_domain;
    double * mr = 0;
    double * mi = 0;
    if (b [nb] .ab_roffset) mr = m + b [nb] .ab_roffset;
    if (b [nb] .ab_ioffset) mi = m + b [nb] .ab_ioffset;
    dimension = density .width  (b [nb] .ab_domain);
    offset    = density .domain () .offset (sub);
    //
    //	Diagonalize diagonal block
    //
    householder (mr, mi, eigen + offset, offd + offset, dimension, dimension);
    if (tqli (eigen + offset, offd + offset, mr, mi, dimension, dimension)) {
      cout << "Block::select_states (): tqli error!" << endl;
      exit (0);
    }
  }
  //
  //	Entropy = - Trace {density * Log_2 (density)}
  //
  double ent = 0.0;
  double trace = 0.0; 
  double eps = machine_precision ();
  for (i = 0; i < states; i++) trace += eigen [i];
  for (i = 0; i < states; i++){
    eigen [i] /= trace;
    autovalori[i] =eigen[i];
    if (eigen [i] > eps) ent -= eigen [i] * log (eigen [i]);
  }
  //==============Ordina solo se richiesti =========
  if ((renyiarg.size() > 1) && (renyiarg[1]<0 ))  {
    double temp=0.0;
    for(size_t k=1;k<states;k++) { 
      for(size_t j=0;j<states-k;j++) { 
	if(autovalori[j] < autovalori[j+1]) { 
	  temp=autovalori[j]; 
	  autovalori[j]=autovalori[j+1]; 
	  autovalori[j+1]=temp; 
	} 
      } 
    }
  }  
  ent /= log (2.0);
  renyientropy .clear ();
  renyientropy .push_back (ent);
  //
  //	Renyi entropy (alpha) = 1/(1-alpha) * Trace {density^alpha}
  //
  //	For alpha == 1 Renyi entropy is equal to previous entropy
  //
  //	Requested alpha values are accumulated with input keyword renyi 
  //	in vector renyiarg
  //
  for (size_t j = 1; j < renyiarg .size (); j++) {
    double alpha = renyiarg [j];
    if (abs (alpha - 1.0) <= eps) renyientropy .push_back (ent);
    else {
      double sum = 0.0;
      for (i = 0; i < states; i++) 
	if (eigen [i] > eps) sum += exp (alpha * log (eigen [i]));
      sum = 1.0/(1.0 - alpha) * log (sum) /log (2.0);
      renyientropy .push_back (sum);
    }  
  }
  return ent;
}
//
//____________________________________________________________________________
void Action::expand (complex<double> * full) const
{
  //
  //	Expand Action to a full complex matrix.
  //	No check on memory, take care!
  //	
  Ablock * b  = block   ();
  double * ma = storage ();
  size_t stride = height ();
  for (size_t n = 0; n < blocks (); n++) {
    size_t l = b [n] .ab_range;
    size_t m = b [n] .ab_domain;
    size_t r = b [n] .ab_roffset;
    size_t i = b [n] .ab_ioffset;
    size_t h = height (l);
    size_t w = width  (m);
    l = a_range  .offset (l);
    m = a_domain .offset (m);
    if (r) 
      for (size_t q = 0; q < w; q++) 
	for (size_t p = 0; p < h; p++) 
	  full [l + p + (m + q) * stride] += ma [r + p + q * h];
    if (i) 
      for (size_t q = 0; q < w; q++) 
	for (size_t p = 0; p < h; p++) 
	  full [l + p + (m + q) * stride] += 
	    complex<double> (0, ma [i + p + q * h]);
  }
}
//
//____________________________________________________________________________
bool Action::iscomplex () const
{
  //
  //	Check for some imaginary part
  //
  if (scalar () .imag () != 0.0) return true;
  Ablock * b = block ();
  for (size_t n = 0; n < blocks (); n++) if (b [n] .ab_ioffset) return true;
  return false;
}
//
//____________________________________________________________________________
bool Action::isnormal () const
{
  //
  //	Check for normalized form
  //
  return ((a_size == csize) || (scalar () == 1.0) || (scalar () == 0.0));
}
//
//____________________________________________________________________________
bool Action::isnull () const
{
  //
  //	Check for null matrix
  //
  return (scalar () == 0.0);
}
//
//____________________________________________________________________________
bool Action::isscalar () const
{
  //
  //	Check for scalar
  //
  return (a_size == csize);
}
//
//____________________________________________________________________________
void Action::merge (const Action & other, size_t tocopy)
{
  //
  //	Merge other's block structure into actual block structure
  //	(block entries are modified) computing new memory size
  //	Non null original memory pointers are not modified, new non 
  //	null memory pointers are added referencing additional memory
  //
  //	If tocopy is null no reallocation of memory is done, the 
  //	required size is only computed and stored in a_size field.
  //
  //	If tocopy is non null memory is reallocated and old elements
  //	are copied into new area.
  //
  if (! ((a_domain == other .a_domain) && (a_range == other .a_range))) {
    cout << "Merging incompatible Actions!" << endl;
    exit (0);
  }
  //
  //	Nothing to do for null other matrix or both scalars
  //
  if (other .isnull ()) return;
  if (isscalar () && other .isscalar ()) return;
  //
  if (other .isscalar ()) {
    //
    //	Transform scalar into a true matrix and merge with diagonal
    //  matrix
    //
    Action b (other);
    b .normalize ();
    if (tocopy) tocopy = a_size;
    merge (b, tocopy);
    return;
  }
  //
  if (isscalar () && (! isnull ())) {
    //
    //	Transform this into a true matrix
    //
    normalize ();
  }
  //
  //	Adjust size to be preserved
  //
  if (tocopy) tocopy = a_size;
  //
  size_t mold = a_size;
  size_t bold = a_blocks;
  vector<Ablock> blk;
  Ablock * b = block ();
  size_t n, nm;
  for (n = 0; n < bold; n++) blk .push_back (b [n]);
  b = other .block ();
  for (n = 0; n < other .blocks (); n++) {
    for (nm = 0; nm < blk .size (); nm++) 
      if ((blk [nm] .ab_range  == b [n] .ab_range) &&
	  (blk [nm] .ab_domain == b [n] .ab_domain)) break;
    if (nm == blk .size ()) blk .push_back (Ablock ());
    blk [nm] .ab_range  = b [n] .ab_range;
    blk [nm] .ab_domain	= b [n] .ab_domain;
    size_t mem = width (blk [nm] .ab_domain) * height (blk [nm] .ab_range);
    if ((blk [nm] .ab_roffset == 0) && b [n] .ab_roffset) {
      blk [nm] .ab_roffset = a_size;
      a_size += mem;
    }
    if ((blk [nm] .ab_ioffset == 0) && b [n] .ab_ioffset) {
      blk [nm] .ab_ioffset = a_size;
      a_size += mem;
    }
  }
  //
  //	If memory size is the same also the structure of blocks is
  //	the same.
  //
  if (a_size == mold)	return;
  //
  //	New memory added. Reallocate and fill blocks 
  //	
  blocks (blk .size ());
  block  (blk .size ());
  b = block ();
  for (n = 0; n < blocks (); n++) b [n] = blk [n];
  //
  //	Update global statistic
  //
  if (a_statistic && (a_statistic != other .statistic ())) a_statistic = 0;
  //
  if (tocopy && (tocopy < a_size)) {
    //
    //	Reallocate memory
    //
    Storage mnew  (a_size * sizeof (double));
    double * mn = (double *) mnew .storage ();
    double * mo = (double *)       storage ();
    while (tocopy--) mn [tocopy] = mo [tocopy];
    storage (mnew);
  }
  //
  // deallocate sparsing
  //
  sparsed (0);
}
//
//____________________________________________________________________________
void Action::normalize ()
{
  //
  //	Modify Action's storage area in order to have a unitary
  //	common factor and an effective matrix part
  //
  // 	Special cases: a scalar is transformed into a diagonal 
  //	matrix.
  //
  if (isscalar ()) scalaridentity (scalar ());
  //
  //	If Action is yet normalized, nothing to do
  //
  if (scalar () == 1.0) return;
  //
  //	Effective scaling. New memory and blocks are allocated and 
  //	defined if needed (this can be the case for complex scalar 
  //	factors).
  //	
  double zr = scalar () .real ();
  double zi = scalar () .imag ();
  Storage bnew (a_blocks * sizeof (Ablock));
  Ablock * bn = (Ablock *) bnew .storage ();
  Ablock * b  = block ();
  //	
  //	redefine blocks entries and compute memory size.
  //
  size_t ssize = csize;
  for (size_t n = 0; n < a_blocks; n++) {
    bn [n] .ab_domain = b [n] .ab_domain;
    bn [n] .ab_range  = b [n] .ab_range;
    size_t ms = height (b [n] .ab_range) * width (b [n] .ab_domain);
    if ((zr && b [n] .ab_roffset) || (zi && b [n] .ab_ioffset)) {
      bn [n] .ab_roffset = ssize;
      ssize += ms;
    }
    if ((zr && b [n] .ab_ioffset) || (zi && b [n] .ab_roffset)) {
      bn [n] .ab_ioffset = ssize;
      ssize += ms;
    }
  }
  //
  //	Set new memory
  //
  Storage  mnew (ssize * sizeof (double));
  double * mn = (double *) mnew .storage ();
  double * m  = storage ();
  for (size_t n = 0; n < a_blocks; n++) {
    size_t ms = height (b [n] .ab_range) * width (b [n] .ab_domain);
    size_t oldro = b  [n] .ab_roffset;
    size_t oldio = b  [n] .ab_ioffset;
    size_t newro = bn [n] .ab_roffset;
    size_t newio = bn [n] .ab_ioffset;
    if (zr && oldro)
      for (size_t n = 0; n < ms; n++) mn [newro + n] += zr * m [oldro + n];
    if (zi && oldio)
      for (size_t n = 0; n < ms; n++) mn [newro + n] -= zi * m [oldio + n];
    if (zr && oldio)
      for (size_t n = 0; n < ms; n++) mn [newio + n] += zr * m [oldio + n];
    if (zi && oldro)
      for (size_t n = 0; n < ms; n++) mn [newio + n] += zi * m [oldro + n];
  }
  //
  //	reload storage grabbing new defined ones
  //
  a_size = ssize;
  storage (mnew);
  scalar  (1.0);
  block   (bnew);
  //
  //	deallocate sparsing
  //
  sparsed (0);
}
//
//____________________________________________________________________________
Action & Action::operator = (const Action & other)
{
  //
  //	Copy assignment
  //
  if (this != & other) {
    a_domain 	= other .a_domain;
    a_range  	= other .a_range;
    a_size   	= other .a_size;
    a_blocks 	= other .a_blocks;
    a_statistic = other .a_statistic;
    if (other .storage ()) {
      size_t size = other .a_storage .size ();
      a_storage .storage (size);
      memcpy (storage (), other .storage (), size);
    }
    else a_storage .storage (0);
    if (other .block ()) {
      size_t size = other .a_block .size ();
      a_block .storage (size);
      memcpy (block (), other .block (), size);
    }
    else a_block .storage (0);
    if (other .sparsed ()) {
      size_t size = other .a_sparsed .size ();
      a_sparsed .storage (size);
      memcpy (sparsed (), other .sparsed (), size);
    }
    else a_sparsed .storage (0);
  }
  return *this;
}
//
//____________________________________________________________________________
Action & Action::operator << (const Action & other)
{
  //
  //	Grab assignment.
  //
  //	The structure and memory is grabbed (no allocation) from other Action
  //
  if (this != & other) {
    a_domain 	= other .a_domain;
    a_range  	= other .a_range;
    a_size   	= other .a_size;
    a_blocks 	= other .a_blocks;
    a_statistic = other .a_statistic;
    a_storage  << other .a_storage;
    a_block    << other .a_block;
    a_sparsed  << other .a_sparsed;
  }
  return *this;
}
//
//____________________________________________________________________________
Action & Action::operator += (const Action & other)
{
  //
  //	sum other to actual 
  //
  //
  //	Special cases
  //
  if (isscalar () && other .isscalar ()) {
    //
    //  scalar + scalar = scalar
    //
    scalar (scalar () + other .scalar ());
    if (scalar () != 0.0) a_statistic = 1;
    return *this;
  }
  if (this == & other) {
    //
    //	same Action --> rescale scalar factor
    //
    scalar (scalar () * 2.0);
    return *this;
  }
  if (isnull ()) {
    //
    //	0 + Action = Action
    //
    *this = other;
    return *this;
  }
  if (other .isnull ())
    //
    //	Action + 0 = Action
    //
    return *this;
  if (! ((a_domain == other .a_domain) && (a_range == other .a_range))) {
    cout << "Adding incompatible Actions!" << endl;
    exit (0);
  }
  //
  //	At least one operand is non scalar, result is non scalar.
  //	Grab operand to a temporary operand and normalize operands.	
  //	(the result is normalized). 
  //	
  normalize ();
  Action b;
  b << other;
  b .normalize ();
  //
  //	Check statistic for the result
  //
  if (b .statistic () != a_statistic) a_statistic = 0;
  //
  //	Merge blocks and realloc matrix entries (old a_size entries
  //	are copied into new entries, in the same positions);
  //
  merge (b, a_size);
  //
  //	Add matrix part of right operand
  //
  double * ma =    storage ();
  double * mb = b .storage ();
  Ablock * ba =    block ();
  Ablock * bb = b .block ();
  size_t na, nb;
  for (nb = 0; nb < b .blocks (); nb++) {
    for (na = 0; na < blocks (); na++)
      if ((bb [nb] .ab_range  == ba [na] .ab_range ) &&
	  (bb [nb] .ab_domain == ba [na] .ab_domain)) break;
    size_t mem = width (ba [na] .ab_domain) * height (ba [na] .ab_range);
    size_t bro = bb [nb] .ab_roffset;
    size_t bio = bb [nb] .ab_ioffset;
    size_t aro = ba [na] .ab_roffset;
    size_t aio = ba [na] .ab_ioffset;
    if (bro) for (size_t j = 0; j < mem; j++) ma [aro + j] += mb [bro +j];
    if (bio) for (size_t j = 0; j < mem; j++) ma [aio + j] += mb [bio +j];
  }
  return *this;
}
//
//____________________________________________________________________________
Action & Action::operator *= (complex<double> factor)
{
  //
  //	Multiply by a scalar
  //
  scalar (scalar () * factor);
  return *this;
}
//
//____________________________________________________________________________
Action & Action::operator *= (const Action & other)
{
  //
  //	Right multiplication
  //
  multiply (*this, *this, other);
  return *this;
}
//
//____________________________________________________________________________
void Action::random ()
{
  //
  //	Fills matrix entries with random number and normalize 
  //	
  if (storage ())
    random_vector (storage () + csize, a_size - csize);
}
//
//____________________________________________________________________________
void Action::release ()
{
  //
  //	Releases storage areas
  //
  a_block   .release ();
  a_storage .release ();
  a_sparsed .release ();
}
//
//____________________________________________________________________________
complex<double> Action::scalar () const
{
  //
  //	get scalar common factor
  //
  if (storage ()) return ((complex<double> *) storage ()) [0];
  return complex<double> (0.0, 0.0);
}
//
//____________________________________________________________________________
void Action::scalar (complex<double> z)
{	
  //
  //	Set scalar factor
  //
  if (storage () == 0) {
    cout << "Setting scalar for Undefined Action!" << endl;
    exit (0);
  }
  ((complex<double> *) storage ()) [0] = z;
}
//
//____________________________________________________________________________
void Action::scalaridentity (complex<double> factor)
{
  //
  //	Transform this to an Action representing the identity times 
  //	scalar factor (acting on space a_domain);
  //	Previous structure destroyed and an identity matrix defined.
  //
  if (storage () == 0) {
    cout << "Undefined Action domain for identity definition!" << endl;
    exit (0);
  }
  a_range = a_domain;
  size_t nb = a_domain .subspaces ();
  blocks (nb-1);
  block  (nb-1);
  Ablock * bb = block ();
  a_size = csize;
  for (size_t n = 1; n < nb; n++) {
    size_t mem  = width (n) * height (n);
    Ablock & b = bb [n-1];
    b .ab_domain  = n;
    b .ab_range   = n;
    b .ab_roffset = a_size;
    b .ab_ioffset = 0;
    a_size += mem;
  }
  storage (a_size);
  double * m = storage ();
  for (size_t n = 1; n < nb; n++) {
    size_t off = bb [n-1] .ab_roffset;
    size_t h   = height (n);
    for (size_t j = 0; j < h; j++) m [off + j + j*h] = 1.0;
  }
  scalar (factor);
  a_statistic = 1;
}
//
//____________________________________________________________________________
void Action::show (const string & s) const
{
  //
  //	Prints out the matrix represented by an Action class
  //
  //	This routine can be very slow with heavy memory request
  //	but is used mainly for control and debugging purpose, so
  //	we don't care.
  //
  ios_base::fmtflags reset = cout .flags ();
  Ablock * b = block  ();
  double * m = storage ();
  //
  //	Special cases
  //
  if (a_size == 0) {
    if (s .size ()) cout << s << " " 
			 << height () << "x" << width () << ": ";
    cout << "Undefined" << endl;
    return;
  }
  if (isnull ()) {
    if (s .size ()) cout << s << " " 
			 << height () << "x" << width () << ": ";
    cout << "Null" << endl;
    return;
  }
  complex<double> & factor = ((complex<double> *) m) [0];
  if (isscalar ()) {
    if (s .size ()) cout << s << " " 
			 << height () << "x" << width () << ": ";
    if (factor != 1.0) cout << factor << " ";
    cout << "Identity" << endl;
    return;
  }
  if (s .size ()) {
    cout << s << " " << height () << "x" << width ();
    //
    //	Add statistic information
    //
    if (a_statistic < 0)	cout << " (fermionic)";
    if (a_statistic > 0)	cout << " (bosonic)";
    cout << ": ";
  }
  if (factor != 1.0) cout << "Common factor = " << factor;
  if (s .size () || (factor != 1.0)) cout << endl;
  //
  //		
  //
  Storage table (2 * width () * (height () + 1) * sizeof (size_t));
  size_t * rc = (size_t *) table .storage ();
  size_t * ic = rc + width ();
  size_t * rt = ic + width ();
  size_t * it = rt + width () * height ();
  //
  //	Scan blocks
  //
  for (size_t bl = 0; bl < blocks (); bl++) {
    size_t roff = b [bl] .ab_roffset;
    size_t ioff = b [bl] .ab_ioffset;
    size_t ib   = b [bl] .ab_range;
    size_t jb   = b [bl] .ab_domain;
    size_t h    = a_range  .states (ib);
    size_t hoff = a_range  .offset (ib);
    size_t w 	= a_domain .states (jb);
    size_t woff = a_domain .offset (jb);
    size_t ij  = 0;
    for (size_t j = 0; j < w; j++) {
      size_t ijt = (woff + j) * height () + hoff; 
      for (size_t i = 0; i < h; i++, ij++, ijt++) {
	if (roff && m [roff + ij]) {
	  rt [ijt] = roff + ij;
	  rc [woff + j]++;
	}
	if (ioff && m [ioff + ij]) {
	  it [ijt] = ioff + ij;
	  ic [woff + j]++;
	}
      }
    }
  }
  //
  //	Print non null elements
  //
  long   index [7];
  long   vv    [7];
  for (size_t j = 0; j < 7; j++) index [j] = vv [j] = 0;
  size_t nj = 0;
  for (size_t j = 0; j < width (); j++) {
    if (rc [j]) index [nj++] =  (j+1);
    if (ic [j]) index [nj++] = -(j+1);
    if ((nj == 6) || 
	(nj && (j == width () - 1)) ||
	((nj == 5) && (j < width () - 1) && rc [j+1] && ic [j+1])) {
      //
      //   Print domain indexes
      //
      cout << setw (6) << "";
      for (size_t jj = 0; jj < nj; jj++) {
	if (index [jj] == -index [jj+1]) {
	  cout << setw (23) << index [jj] - 1 << setw (1) << ")"; 
	  jj++;
	}
	else cout << right << setw (11) << abs (index [jj]) - 1 << ")";
      }
      cout << endl;
      //
      //	Print range indexes and non null values
      //
      for (size_t i = 0; i < height (); i++) {
	size_t ni = 0;
	for (size_t jj = 0; jj < nj; jj++) {
	  vv [jj] = 0;
	  long ij = index [jj];
	  if (ij > 0)      ij = rt [i + (ij-1) * height ()];
	  else if (ij < 0) ij = it [i - (ij+1) * height ()];
	  if (ij) {
	    ni++;
	    vv [jj] = ij;
	  }
	}
	if (ni) {
	  //
	  //	Print non null elements in line i
	  //
	  cout << setw (5) << right << i << setw (1) << ")";
	  for (size_t jj = 0; jj < nj; jj++) {
	    if (index [jj] == -index [jj+1]) {
	      //
	      //	Print full complex 
	      //
	      if (vv [jj] + vv [jj+1]) {
		double vr = 0.0;
		if (vv [jj]) vr = m [vv [jj]];
		double vi = 0.0;
		if (vv [jj+1]) vi = m [vv [jj+1]];
		if ((vv [jj] == 0) || 
		    ((fabs (vr) + 0.00000005 < 10.0) &&
		     (fabs (vr) + 0.00000005 > 0.01)))
		  cout << setw (12) << right << fixed << setprecision (7)
		       << vr;
		else cout << setw (12) << right << scientific 
			  << setprecision (3) << vr;
		if ((vv [jj+1] == 0) ||
		    ((fabs (vi) + 0.00000005 < 10.0) &&
		     (fabs (vi) + 0.00000005 > 0.01)))
		  cout << setw (10) << right << fixed << setprecision (7)
		       << showpos << vi << "*i" << noshowpos;
		else cout << setw (10) << right << scientific 
			  << setprecision (3) << showpos << vi 
			  << "*i" << noshowpos;
	      }
	      else cout << setw (24) << "";
	      jj++;
	      continue;
	    }	
	    if (vv [jj] == 0) cout << setw (12) << "";
	    else {
	      double v = m [vv [jj]];
	      if ((fabs (v) + 0.00000005 < 10.0) && 
		  (fabs (v) + 0.00000005 > 0.01)) {
		if (index [jj] > 0) 
		  cout << setw (12) << right << fixed << setprecision (7)
		       << v;
		else cout << setw (10) << showpos << right << fixed 
			  << setprecision (7) << v << "*i"
			  << noshowpos;
	      }
	      else {
		if (index [jj] > 0) 
		  cout << setw (12) << right << scientific 
		       << setprecision (3) << v;
		else cout << setw (10) << showpos << right << scientific 
			  << setprecision (3) << v << "*i" << noshowpos;
	      }      
	    }
	  }
	  cout << endl;
	}
      }
      nj = 0;
    }
  }
  //
  //	reset output options
  //
  cout .flags (reset);
  cout .precision (6);
}
//
//____________________________________________________________________________
void Action::showblocks (const string & s) const
{
  //
  //   prints the structure of blocks 
  //
  size_t n, nb;
  nb = blocks ();
  if (s .size ()) cout << s << " ";
  cout << "blocks () = " <<  nb << endl;
  Ablock * b = block ();
  for (n = 0; n < nb; n++) 
    cout << "block[" << n << "] " << b[n] .ab_domain 
	 << "-->" << b[n] .ab_range 
	 << " (" <<  b[n] .ab_roffset
	 << "," <<   b[n] .ab_ioffset
	 << ") " <<  b[n] .ab_cindex 
	 << " " <<   b[n] .ab_rindex << endl;
}
//
//____________________________________________________________________________
void Action::sparse ()
{
  //
  //	Build sparse indexes list 
  //
  if (blocks () == 0) return;
  //
  //	deallocate sparsing list (just in case)
  //
  sparsed (0);
  //
  Ablock * b = block   ();
  double * m = storage ();
  //
  //	Counting not null elements for sparsing allocation
  //
  double eps = machine_precision ();
  size_t spsize = 0;
  size_t all    = 0;
  size_t cp     = 2;
  for (size_t n = 0; n < blocks (); n++) {
    double * mr;
    double * mi;
    size_t h = height (b [n] .ab_range);
    size_t w = width  (b [n] .ab_domain);
    mr = 0;
    if (b [n] .ab_roffset) mr = m + b [n] .ab_roffset;
    mi = 0;
    if (b [n] .ab_ioffset) mi = m + b [n] .ab_ioffset;
    size_t nonzero = 0.0;
    for (size_t i = 0; i < h * w; i++) {
      double vr = 0.0;
      double vi = 0.0;
      if (mr) {
	if (fabs (mr [i]) < eps) mr [i] = 0.0;
	else vr = mr [i];
      }
      if (mi) {
	if (fabs (mi [i]) < eps) mi [i] = 0.0;
	else vi = mi [i];
      }
      if (vr || vi) nonzero++;
    }
    if (sparse_den * nonzero < h * w * sparse_num) {
      b [n] .ab_cindex = cp;
      cp += (w + 1);
      b [n] .ab_rindex = cp;
      cp += (h + 1);
      spsize += nonzero;
    }
    all    += h * w;
  }
  if (cp == 2) return;
  //
  //	Build indexes list
  //
  sparsed (2 * spsize + cp);
  size_t * p  = sparsed ();
  p [0] = spsize;
  p [1] = all;
  //
  spsize = cp;
  for (size_t n = 0; n < blocks (); n++) {
    double * mr = 0;
    double * mi = 0;
    if (b [n] .ab_roffset) mr = m + b [n] .ab_roffset;
    if (b [n] .ab_ioffset) mi = m + b [n] .ab_ioffset;
    size_t h = height (b [n] .ab_range);
    size_t w = width  (b [n] .ab_domain);
    //
    if (b [n] .ab_cindex) {
      //
      //	column-indexed sparsing indexes
      //
      cp = b [n] .ab_cindex;
      for (size_t j = 0; j < w; j++) {
	p [cp++] = spsize;
	for (size_t i = 0; i < h; i++) {
	  double vr = 0.0;
	  double vi = 0.0;
	  if (mr) vr = mr [i + j * h];
	  if (mi) vi = mi [i + j * h];
	  if (vr || vi) p [spsize++] = i;
	}
      }
      p [cp++] = spsize;
    }
    if (b [n] .ab_rindex) {
      //
      //	row-inedexed sparsing indexes
      //
      cp = b [n] .ab_rindex;
      for (size_t i = 0; i < h; i++) {
	p [cp++] = spsize;
	for (size_t j = 0; j < w; j++) {
	  double vr = 0.0;
	  double vi = 0.0;
	  if (mr) vr = mr [i + j * h];
	  if (mi) vi = mi [i + j * h];
	  if (vr || vi) p [spsize++] = j;
	}
      }
      p [cp++] = spsize;
    }
  }
}
//
//____________________________________________________________________________
void Action::sparsed (size_t size)
{
  //
  //	Allocate sparse indexes array
  //
  a_sparsed .storage (size * sizeof (size_t));
  if ((size == 0) & blocks ()) {
    //
    //  Remove sparse info
    //
    Ablock * b = block ();
    for (size_t n = 0; n < blocks (); n++) 
      b [n] .ab_cindex = b [n] .ab_rindex = 0;
  }
}
//
//____________________________________________________________________________
void Action::transpose ()
{
  //
  //	Transpose transformation (no adjoint)
  //
  double * m = storage ();
  if (m == 0) {
    cout << "Transpose of invalid Action!" << endl;
    exit (0);
  }
  //
  //	All done for a scalar
  //
  if (isscalar ()) return;
  Storage t (a_size * sizeof (double));
  double * mt = (double *) t .storage ();
  ((complex<double> *) mt) [0] = ((complex<double> *) m) [0];
  Storage  tb   (blocks () * sizeof (Ablock));
  Ablock * bt = (Ablock *) tb .storage ();
  Ablock * b = block ();
  for (size_t nb = 0; nb < blocks (); nb++) {
    bt [nb] = b [nb];
    size_t roff = b [nb] .ab_roffset;
    size_t ioff = b [nb] .ab_ioffset;
    size_t is   = b [nb] .ab_range;
    size_t js   = b [nb] .ab_domain;
    size_t w    = width  (js);
    size_t h    = height (is);
    if (roff) 
      for (size_t j = 0; j < w; j++) 
	for (size_t i = 0; i < h; i++) 
	  mt [roff + j + i*w] =  m [roff +i +j*h];
    if (ioff) 
      for (size_t j = 0; j < w; j++) 
	for (size_t i = 0; i < h; i++) 
	  mt [ioff + j + i*w] = m [ioff +i +j*h];
    bt [nb] .ab_range  = js;
    bt [nb] .ab_domain = is;
  }
  storage (t);
  block   (tb);
  Space s = a_domain;
  a_domain = a_range;
  a_range  = s;
  //
  //	deallocate sparsing
  //
  sparsed (0);
}
//
//============================================================================
void binary (vector<Abinary> & list, const Action & res, 
	     const Action & lft, const Action & rgt)
{
  //
  //	Setup list for binary product of Action's:
  //
  //		lft * rgt --> res
  //
  list .clear ();
  if (lft .isnull () || rgt .isnull ()) return;
  if (lft .isscalar () && rgt .isscalar ()) {
    //
    //	scalar * scalar = scalar
    //
    list .push_back (Abinary ());
    return;
  }
  Ablock * blft = lft .block ();
  Ablock * brgt = rgt .block ();
  Ablock * bres = res .block ();
  if (lft .isscalar ()) {
    //
    //	scalar * matrix = matrix
    //
    double zr = lft .scalar () .real ();
    double zi = lft .scalar () .imag ();
    for (size_t nrgt = 0; nrgt < rgt .blocks (); nrgt++) {
      size_t nres;
      for (nres = 0; nres < res .blocks (); nres++)
	if ((brgt [nrgt] .ab_domain == bres [nres] .ab_domain) &&
	    (brgt [nrgt] .ab_range  == bres [nres] .ab_range)) break;
      Abinary p;
      p .ab_intsz = res .width (bres [nres] .ab_domain) * 
	res .height (bres [nres] .ab_range);
      if (zr && brgt [nrgt] .ab_roffset) {
	p .ab_resro = bres [nres] .ab_roffset;
	p .ab_rgtro = brgt [nrgt] .ab_roffset;
      }
      if (zr && brgt [nrgt] .ab_ioffset) {
	p .ab_resio = bres [nres] .ab_ioffset;
	p .ab_rgtio = brgt [nrgt] .ab_ioffset;
      }
      if (zi && brgt [nrgt] .ab_roffset) {
	p .ab_resio = bres [nres] .ab_ioffset;
	p .ab_rgtro = brgt [nrgt] .ab_roffset;
      }
      if (zi && brgt [nrgt] .ab_ioffset) {
	p .ab_resro = bres [nres] .ab_roffset;
	p .ab_rgtio = brgt [nrgt] .ab_ioffset;
      }
      if (p .ab_rgtro || p .ab_rgtio) list .push_back (p);
    }
    return;
  }
  if (rgt .isscalar ()) {
    //
    //	matrix * scalar = matrix
    //
    double zr = rgt .scalar () .real ();
    double zi = rgt .scalar () .imag ();
    for (size_t nlft = 0; nlft < lft .blocks (); nlft++) {
      size_t nres;
      for (nres = 0; nres < res .blocks (); nres++)
	if ((blft [nlft] .ab_domain == bres [nres] .ab_domain) &&
	    (blft [nlft] .ab_range  == bres [nres] .ab_range)) break;
      Abinary p;
      p .ab_intsz = res .width (bres [nres] .ab_domain) * 
	res .height (bres [nres] .ab_range);
      if (zr && blft [nlft] .ab_roffset) {
	p .ab_resro = bres [nres] .ab_roffset;
	p .ab_lftro = blft [nlft] .ab_roffset;
      }
      if (zr && blft [nlft] .ab_ioffset) {
	p .ab_resio = bres [nres] .ab_ioffset;
	p .ab_lftio = blft [nlft] .ab_ioffset;
      }
      if (zi && blft [nlft] .ab_roffset) {
	p .ab_resio = bres [nres] .ab_ioffset;
	p .ab_lftro = blft [nlft] .ab_roffset;
      }
      if (zi && blft [nlft] .ab_ioffset) {
	p .ab_resro = bres [nres] .ab_roffset;
	p .ab_lftio = blft [nlft] .ab_ioffset;
      }
      if (p .ab_lftro || p .ab_lftio) list .push_back (p);
    }
    return;
  }
  //
  //	matrix * matrix = matrix
  //
  for (size_t nrgt = 0; nrgt < rgt .blocks (); nrgt++) {
    if ((brgt [nrgt] .ab_roffset == 0) && (brgt [nrgt] .ab_ioffset == 0)) 
      continue;
    for (size_t nlft = 0; nlft < lft .blocks (); nlft++) {
      if ((blft [nlft] .ab_roffset == 0) && (blft [nlft] .ab_ioffset == 0)) 
	continue;
      if (blft [nlft] .ab_domain != brgt [nrgt] .ab_range) continue;
      size_t nres;
      for (nres = 0; nres < res .blocks (); nres++) 
	if ((bres [nres] .ab_range  == blft [nlft] .ab_range) && 
	    (bres [nres] .ab_domain == brgt [nrgt] .ab_domain)) break;
      Abinary p;
      p .ab_lftsz = lft .height (blft [nlft] .ab_range);
      p .ab_intsz = rgt .height (brgt [nrgt] .ab_range);
      p .ab_rgtsz = rgt .width  (brgt [nrgt] .ab_domain);
       if (lft .range () .statistic (blft [nlft] .ab_range) < 0) 
	p .ab_lftsz = - p .ab_lftsz;
      if (rgt .range () .statistic (brgt [nrgt] .ab_range) < 0) 
	p .ab_intsz = - p .ab_intsz;
      if (rgt .domain () .statistic (brgt [nrgt] .ab_domain) < 0) 
	p .ab_rgtsz = - p .ab_rgtsz;
      //
      //  Set sparse indexes
      //
      if (blft [nlft] .ab_cindex) p .ab_lftsp = blft [nlft] .ab_cindex;
      if (brgt [nrgt] .ab_rindex) p .ab_rgtsp = brgt [nrgt] .ab_rindex;
      //
      if (blft [nlft] .ab_roffset && brgt [nrgt] .ab_roffset) {
	p .ab_resro = bres [nres] .ab_roffset;
	p. ab_lftro = blft [nlft] .ab_roffset;
	p .ab_rgtro = brgt [nrgt] .ab_roffset;
      }
      if (blft [nlft] .ab_roffset && brgt [nrgt] .ab_ioffset) {
	p .ab_resio = bres [nres] .ab_ioffset;
	p. ab_lftro = blft [nlft] .ab_roffset;
	p .ab_rgtio = brgt [nrgt] .ab_ioffset;
      }
      if (blft [nlft] .ab_ioffset && brgt [nrgt] .ab_roffset) {
	p .ab_resio = bres [nres] .ab_ioffset;
	p. ab_lftio = blft [nlft] .ab_ioffset;
	p .ab_rgtro = brgt [nrgt] .ab_roffset;
      }
      if (blft [nlft] .ab_ioffset && brgt [nrgt] .ab_ioffset) {
	p .ab_resro = bres [nres] .ab_roffset;
	p. ab_lftio = blft [nlft] .ab_ioffset;
	p .ab_rgtio = brgt [nrgt] .ab_ioffset;
      }
      if (p .ab_resro || p .ab_resio) list .push_back (p);
    }
  }
}
//
//____________________________________________________________________________
void binary (vector<Abinary> & list, 
	     const Action & lft, const Action & rgt)
{
  //
  //	Setup for scalar product of two actions
  //
  list .clear ();
  Ablock * blft = lft .block ();
  Ablock * brgt = rgt .block ();
  for (size_t nrgt = 0; nrgt < rgt .blocks (); nrgt++) {
    if ((brgt [nrgt] .ab_roffset == 0) && (brgt [nrgt] .ab_ioffset == 0)) 
      continue;
    for (size_t nlft = 0; nlft < lft .blocks (); nlft++) {
      if ((blft [nlft] .ab_roffset == 0) && (blft [nlft] .ab_ioffset == 0)) 
	continue;
      if ((blft [nlft] .ab_range  != brgt [nrgt] .ab_range) ||
	  (blft [nlft] .ab_domain != brgt [nrgt] .ab_domain)) continue;
      Abinary p;
      p .ab_intsz = rgt .width (brgt [nrgt] .ab_domain) * 
	rgt .height (brgt [nrgt] .ab_range);
      p .ab_lftro = blft [nlft] .ab_roffset;
      p .ab_rgtro = brgt [nrgt] .ab_roffset;
      p .ab_lftio = blft [nlft] .ab_ioffset;
      p .ab_rgtio = brgt [nrgt] .ab_ioffset;
      if ((p .ab_lftro || p .ab_lftio) && (p .ab_rgtro || p .ab_rgtio)) 
	list .push_back (p);
    }
  }
}
//
//____________________________________________________________________________
size_t minasize ()
{
  //
  //	Returns the minimum size for a defined Action
  //
  return csize;
}
//
//____________________________________________________________________________
complex<double> multiply (double * a, double * b, vector<Abinary> & list)
{
  //
  //	Computes the scalar product of two Action's memory area
  //	(Schmidt scalar product)
  //
  double scr = 0.0;
  double sci = 0.0;
  for (size_t np = 0; np < list .size (); np++) {
    Abinary & p = list [np];
    size_t lftro = p .ab_lftro;
    size_t lftio = p. ab_lftio;
    size_t rgtro = p .ab_rgtro;
    size_t rgtio = p .ab_rgtio;
    size_t intsz = p .ab_intsz;
    if (lftro && rgtro) 
      for (size_t k = 0; k < intsz; k++) scr += a [lftro + k] * b [rgtro + k];
    if (lftro && rgtio) 
      for (size_t k = 0; k < intsz; k++) sci += a [lftro + k] * b [rgtio + k];
    if (lftio && rgtro) 
      for (size_t k = 0; k < intsz; k++) sci -= a [lftio + k] * b [rgtro + k];
    if (lftio && rgtio) 
      for (size_t k = 0; k < intsz; k++) scr += a [lftio + k] * b [rgtio + k];
  }
  return complex<double> (scr,sci);
}
//
//____________________________________________________________________________
void multiply (double * r, complex<double> z, double * b, 
	       vector<Abinary> & list)
{
  //
  //	Adds to r the product z * b
  //
  double zr = z .real ();
  double zi = z .imag ();
  for (size_t np = 0; np < list .size (); np++) {
    Abinary & p = list [np];
    size_t resro = p .ab_resro;
    size_t resio = p. ab_resio;
    size_t rgtro = p .ab_rgtro;
    size_t rgtio = p .ab_rgtio;
    size_t intsz = p .ab_intsz;
    if (zr && rgtro) 
      for (size_t k = 0; k < intsz; k++) r [resro + k] += zr * b [rgtro + k];
    if (zr && rgtio) 
      for (size_t k = 0; k < intsz; k++) r [resio + k] += zr * b [rgtio + k];
    if (zi && rgtro) 
      for (size_t k = 0; k < intsz; k++) r [resio + k] += zi * b [rgtro + k];
    if (zi && rgtio) 
      for (size_t k = 0; k < intsz; k++) r [resro + k] -= zi * b [rgtio + k];
  }
}
//
//____________________________________________________________________________
void multiply (double * p, double * a, double * b, 
	       vector<Abinary> & list, bool fermi)
{
  //
  //   	Performs a set of multiplication of block matrices according to
  //	list, taking into account Fermi-Bose statistic if fermi is true
  //	
  //	The non-scalar Actions are assumed normalized to unitary scalar
  //	factors
  //
  double zar = ((complex<double> *) a) [0] .real ();
  double zai = ((complex<double> *) a) [0] .imag ();
  double zbr = ((complex<double> *) b) [0] .real ();
  double zbi = ((complex<double> *) b) [0] .imag ();
  for (size_t np = 0; np < list .size (); np++) {
    Abinary & pr     = list [np];
    size_t  lftro  = pr .ab_lftro;
    size_t  lftio  = pr .ab_lftio;
    size_t  rgtro  = pr .ab_rgtro;
    size_t  rgtio  = pr .ab_rgtio;
    size_t  resro  = pr .ab_resro;
    size_t  resio  = pr .ab_resio;
    long    lftsz  = pr .ab_lftsz;
    long    rgtsz  = pr .ab_rgtsz;
    long    intsz  = pr .ab_intsz;
    double  factor = 1.0;
    if (fermi && (lftsz < 0) && (rgtsz * intsz < 0)) factor = -1.0;
    if (lftsz < 0) lftsz = - lftsz;
    if (rgtsz < 0) rgtsz = - rgtsz;
    if (intsz < 0) intsz = - intsz;
    if ((lftro || lftio) && (rgtro || rgtio)) {
      //
      //	Non scalar * non scalar (non scalar result)
      //
      if (lftro && rgtro && resro) 
	for (long k = 0; k < intsz; k++)
	  for (long j = 0; j < rgtsz; j++)
	    for (long i = 0; i < lftsz; i++)
	      p [resro + i + j*lftsz] +=
		a [lftro + i + k*lftsz] * b [rgtro + k + j*intsz] * factor;
      if (lftio && rgtio && resro) 
	for (long k = 0; k < intsz; k++)
	  for (long j = 0; j < rgtsz; j++)
	    for (long i = 0; i < lftsz; i++)
	      p [resro + i + j*lftsz] -=
		a [lftio + i + k*lftsz] * b [rgtio + k + j*intsz] * factor;
      if (lftro && rgtio && resio) 
	for (long k = 0; k < intsz; k++)
	  for (long j = 0; j < rgtsz; j++)
	    for (long i = 0; i < lftsz; i++)
	      p [resio + i + j*lftsz] +=
		a [lftro + i + k*lftsz] * b [rgtio + k + j*intsz] * factor;
      if (lftio && rgtro && resio) 
	for (long k = 0; k < intsz; k++)
	  for (long j = 0; j < rgtsz; j++)
	    for (long i = 0; i < lftsz; i++)
	      p [resio + i + j*lftsz] +=
		a [lftio + i + k*lftsz] * b [rgtro + k + j*intsz] * factor;
    }
    else if (rgtro || rgtio) {
      //	
      //	Scalar * non scalar	(non scalar result)
      //
      if (resro && rgtro) 
	for (long k = 0; k < intsz; k++) p [resro + k] += zar * b [rgtro + k];
      if (resro && rgtio) 
	for (long k = 0; k < intsz; k++) p [resro + k] -= zai * b [rgtio + k];
      if (resio && rgtro) 
	for (long k = 0; k < intsz; k++) p [resio + k] += zai * b [rgtro + k];
      if (resio && rgtio) 
	for (long k = 0; k < intsz; k++) p [resio + k] += zar * b [rgtio + k];
    }
    else if (lftro || lftio) {
      //
      //	Non scalar * scalar	(non scalar result)
      //
      if (resro && lftro) 
	for (long k = 0; k < intsz; k++) p [resro + k] += a [lftro + k] * zbr;
      if (resro && lftio) 
	for (long k = 0; k < intsz; k++) p [resro + k] -= a [lftio + k] * zbi;
      if (resio && lftro) 
	for (long k = 0; k < intsz; k++) p [resio + k] += a [lftro + k] * zbi;
      if (resio && lftio) 
	for (long k = 0; k < intsz; k++) p [resio + k] += a [lftio + k] * zbr;
    }
    else 
      //
      //	Scalar * scalar	(scalar result)
      //
      ((complex<double> *) p) [0] += 
	complex<double> (zar * zbr - zai * zbi, zar * zbi + zai * zbr);
  }
}
//
//____________________________________________________________________________
void multiply (double * p, double * a, size_t * ia, double * b, size_t * ib, 
	       vector<Abinary> & list, bool fermi)
{
  //
  //   	Performs a set of multiplication of block matrices according to
  //	list, taking into account Fermi-Bose statistic if fermi is true
  //	
  //	The non-scalar Actions are assumed normalized to unitary scalar
  //	factors
  //
  double zar = ((complex<double> *) a) [0] .real ();
  double zai = ((complex<double> *) a) [0] .imag ();
  double zbr = ((complex<double> *) b) [0] .real ();
  double zbi = ((complex<double> *) b) [0] .imag ();
  for (size_t np = 0; np < list .size (); np++) {
    Abinary & pr     = list [np];
    size_t  lftro  = pr .ab_lftro;
    size_t  lftio  = pr .ab_lftio;
    size_t  lftsp  = pr .ab_lftsp;
    size_t  rgtro  = pr .ab_rgtro;
    size_t  rgtio  = pr .ab_rgtio;
    size_t  rgtsp  = pr .ab_rgtsp;
    size_t  resro  = pr .ab_resro;
    size_t  resio  = pr .ab_resio;
    long    lftsz  = pr .ab_lftsz;
    long    rgtsz  = pr .ab_rgtsz;
    long    intsz  = pr .ab_intsz;
    double  factor = 1.0;
    if (fermi && (lftsz < 0) && (rgtsz * intsz < 0)) factor = -1.0;
    if (lftsz < 0) lftsz = - lftsz;
    if (rgtsz < 0) rgtsz = - rgtsz;
    if (intsz < 0) intsz = - intsz;
    //
    //	left and right sparsed is not implemented, ignore right sparsed
    //
    if (lftsp && rgtsp) rgtsp = 0;
    //
    if ((lftro || lftio) && (rgtro || rgtio)) {
      //
      //	Non scalar * non scalar (non scalar result)
      //
      if (lftro && rgtro && resro && lftsp)
	//
	//	left sparsed (column-indexed)
	//
	for (long k = 0; k < intsz; k++) {
	  double * bb = b + rgtro + k;
	  for (size_t ii = ia [lftsp + k]; ii < ia [lftsp + k + 1]; ii++) {
	    size_t i = ia [ii];
	    double   aa = a [lftro + i + k *lftsz] * factor;
	    double * pp = p + resro + i;
	    for (long j = 0; j < rgtsz; j++) 
	      pp [j * lftsz] += aa * bb [j * intsz];
	  }
	}
      else if (lftro && rgtro && resro && rgtsp) 
	//
	//	right sparsed (row-indexed)
	//
	for (long k = 0; k < intsz; k++) {
	  double * aa = a + lftro + k * lftsz;
	  for (size_t jj = ib [rgtsp + k]; jj < ib [rgtsp + k + 1]; jj++) {
	    size_t j = ib [jj];
	    double   bb = b [rgtro + k + j *intsz] * factor;
	    double * pp = p + resro + j * lftsz;
	    for (long i = 0; i < lftsz; i++)
	      pp [i] += aa [i] * bb;
	  }
	}
      else if (lftro && rgtro && resro) {
	fbmm (0, 0, lftsz, rgtsz, intsz, factor,
	      a + lftro, lftsz, b + rgtro, intsz, p + resro, lftsz);
	/*
	for (long k = 0; k < intsz; k++) {
	  double * aa = a + lftro + k * lftsz;
	  for (long j = 0; j < rgtsz; j++) {
	    double   bb = b [rgtro + k + j * intsz] * factor;
	    double * pp = p + resro + j * lftsz;
	    for (long i = 0; i < lftsz; i++) pp [i] += aa [i] * bb;
	  }
	}
	*/
      }
      //
      if (lftio && rgtio && resro && lftsp)
	//
	//	left sparsed (column-indexed)
	//
	for (long k = 0; k < intsz; k++) {
	  double * bb = b + rgtio + k;
	  for (size_t ii = ia [lftsp + k]; ii < ia [lftsp + k + 1]; ii++) {
	    size_t i = ia [ii];
	    double  aa = a [lftio + i + k *lftsz] * factor;
	    double * pp = p + resro + i;
	    for (long j = 0; j < rgtsz; j++) 
	      pp [j * lftsz] -= aa * bb [j * intsz];
	  }
	}
      else if (lftio && rgtio && resro && rgtsp) 
	//
	//	right sparsed (row-indexed)
	//
	for (long k = 0; k < intsz; k++) {
	  double * aa = a + lftio + k * lftsz;
	  for (size_t jj = ib [rgtsp + k]; jj < ib [rgtsp + k + 1]; jj++) {
	    size_t   j = ib [jj];
	    double   bb = b [rgtio + k + j *intsz] * factor;
	    double * pp = p + resro + j * lftsz;
	    for (long i = 0; i < lftsz; i++) pp [i] -= aa [i] * bb;
	  }
	}
      else if (lftio && rgtio && resro) {
	fbmm (0, 0, lftsz, rgtsz, intsz, -factor,
	      a + lftio, lftsz, b + rgtio, intsz, p + resro, lftsz);
	/*
	for (long k = 0; k < intsz; k++) {
	  double * aa = a + lftio + k * lftsz;
	  for (long j = 0; j < rgtsz; j++) {
	    double   bb = b [rgtio + k + j * intsz] * factor;
	    double * pp = p + resro + j * lftsz;
	    for (long i = 0; i < lftsz; i++) pp [i] -= aa [i] * bb;
	  }
	}
	*/
      }
      //
      if (lftro && rgtio && resio && lftsp)
	//
	//	left sparsed (column-indexed)
	//
	for (long k = 0; k < intsz; k++) {
	  double * bb = b + rgtio + k;
	  for (size_t ii = ia [lftsp + k]; ii < ia [lftsp + k + 1]; ii++) {
	    size_t   i  = ia [ii];
	    double   aa = a [lftro + i + k * lftsz] * factor;
	    double * pp = p + resio + i;
	    for (long j = 0; j < rgtsz; j++) 
	      pp [j * lftsz] += aa * bb [j * intsz];
	  }
	}
      else if (lftro && rgtio && resio && rgtsp) 
	//
	//	right sparsed (row-indexed)
	//
	for (long k = 0; k < intsz; k++) {
	  double * aa = a + lftro + k * lftsz;
	  for (size_t jj = ib [rgtsp + k]; jj < ib [rgtsp + k + 1]; jj++) {
	    size_t   j  = ib [jj];
	    double   bb = b [rgtio + k + j * intsz] * factor;
	    double * pp = p + resio + j * lftsz;
	    for (long i = 0; i < lftsz; i++) pp [i] += aa [i] * bb;
	  }
	}
      else if (lftro && rgtio && resio) {
	fbmm (0, 0, lftsz, rgtsz, intsz, factor,
	      a + lftro, lftsz, b + rgtio, intsz, p + resio, lftsz);
	/*
	for (long k = 0; k < intsz; k++) {
	  double * aa = a + lftro + k * lftsz;
	  for (long j = 0; j < rgtsz; j++) {
	    double   bb = b [rgtio + k + j * intsz] * factor;
	    double * pp = p + resio + j * lftsz;
	    for (long i = 0; i < lftsz; i++) pp [i] += aa [i] * bb;
	  }
	}
	*/
      }
      //
      if (lftio && rgtro && resio && lftsp)
	//
	//	left sparsed (column-indexed)
	//
	for (long k = 0; k < intsz; k++) {
	  double * bb = b + rgtro + k;
	  for (size_t ii = ia [lftsp + k]; ii < ia [lftsp + k + 1]; ii++) {
	    size_t   i = ia [ii];
	    double   aa = a [lftio + i + k *lftsz] * factor;
	    double * pp = p + resio + i;
	    for (long j = 0; j < rgtsz; j++) 
	      pp [j * lftsz] += aa * bb [j * intsz];
	  }
	}
      else if (lftio && rgtro && resio && rgtsp) 
	//
	//	right sparsed (row-indexed)
	//
	for (long k = 0; k < intsz; k++) {
	  double * aa = a + lftio + k * lftsz;
	  for (size_t jj = ib [rgtsp + k]; jj < ib [rgtsp + k + 1]; jj++) {
	    size_t   j  = ib [jj];
	    double   bb = b [rgtro + k + j *intsz] * factor;
	    double * pp = p + resio + j * lftsz;
	    for (long i = 0; i < lftsz; i++) pp [i] += aa [i] * bb;
	  }
	}
      else if (lftio && rgtro && resio) {
	fbmm (0, 0, lftsz, rgtsz, intsz, factor,
	      a + lftio, lftsz, b + rgtro, intsz, p + resio, lftsz);
	/*
	for (long k = 0; k < intsz; k++) {
	  double * aa = a + lftio + k * lftsz;
	  for (long j = 0; j < rgtsz; j++) {
	    double   bb = b [rgtro + k + j * intsz] * factor;
	    double * pp = p + resio + j * lftsz;
	    for (long i = 0; i < lftsz; i++) pp [i] += aa [i] * bb;
	  }
	}
	*/
      }
    }
    else if (rgtro || rgtio) {
      //	
      //	Scalar * non scalar	(non scalar result)
      //
      if (resro && rgtro) 
	for (long k = 0; k < intsz; k++) p [resro + k] += zar * b [rgtro + k];
      if (resro && rgtio) 
	for (long k = 0; k < intsz; k++) p [resro + k] -= zai * b [rgtio + k];
      if (resio && rgtro) 
	for (long k = 0; k < intsz; k++) p [resio + k] += zai * b [rgtro + k];
      if (resio && rgtio) 
	for (long k = 0; k < intsz; k++) p [resio + k] += zar * b [rgtio + k];
    }
    else if (lftro || lftio) {
      //
      //	Non scalar * scalar	(non scalar result)
      //
      if (resro && lftro) 
	for (long k = 0; k < intsz; k++) p [resro + k] += a [lftro + k] * zbr;
      if (resro && lftio) 
	for (long k = 0; k < intsz; k++) p [resro + k] -= a [lftio + k] * zbi;
      if (resio && lftro) 
	for (long k = 0; k < intsz; k++) p [resio + k] += a [lftro + k] * zbi;
      if (resio && lftio) 
	for (long k = 0; k < intsz; k++) p [resio + k] += a [lftio + k] * zbr;
    }
    else
      //
      //	Scalar left * scalar right	(scalar result)
      //
      ((complex<double> *) p) [0] += 
	complex<double> (zar * zbr - zai * zbi, zar * zbi + zai * zbr);
 }
}
//
//____________________________________________________________________________
void multiply (Action & result, const Action & a, const Action & b, 
	       bool fermi)
{
  //
  //	Compute result = a * b
  //
  vector<Abinary> list;
  // 
  //	use an intermediate Action to take into account the cases 
  //	result == a or result == b
  //
  Action r;	
  if (a .isscalar ()) {
    r = b;
    r *= a .scalar ();
    result << r;
    return;
  }
  if (b .isscalar ()) {
    r = a;
    r *= b .scalar ();
    result << r;
    return;
  }
  r = Action (a, b);
  binary (list, r, a, b);
  r .storage (r .size ());
  if (r .size () > csize) r .scalar  (a .scalar () * b .scalar ());
  else			  r .scalar  (0.0);
  double * ma = a .storage ();
  double * mb = b .storage ();
  double * mr = r .storage ();
  multiply (mr, ma, mb, list, fermi);
  r .statistic (a .statistic () * b .statistic ());
  result << r;
}
//
//____________________________________________________________________________
complex<double> multiply (const Action & bra, const Action & ket)
{
  //
  //	Compute scalar product <bra|ket> (Schmidt scalar product)
  //
  vector<Abinary> scalar;
  //
  //	Prepare Abinary list
  //
  binary (scalar, bra, ket);
  //
  //	Effective scalar product of memory areas
  //
  complex<double> product = 
    multiply (bra .storage (), ket .storage (), scalar);
  product *= conj(bra .scalar ()) * ket .scalar ();
  return product;
}
//
//============================================================================
