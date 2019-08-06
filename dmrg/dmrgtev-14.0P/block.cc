//	block.cc			(C) 2009 Fabio Ortolani	fbo 090716
//	==================================================================
#include "algebra.hh"
#include "block.hh"
#include "numerical.hh"
#include <cstring>
#include <iostream>
#include <sstream>
#include <iomanip>
//
//============================================================================
bool quantum_alternate (size_t);				  // [menu.cc]
//
//============================================================================
extern size_t	idop_none;
extern size_t	idop_base;
extern size_t	idop_hamiltonian;
extern bool	show_states;
extern bool	show_selection;
extern size_t	show_density;
extern double	truncation;
//
//____________________________________________________________________________
Barray	blockempty;
//
//============================================================================
Barray::Barray ()
  : ba_block ()
{ 
  //
  //	Default ctor.
  //
}
//
//____________________________________________________________________________
Barray::~Barray ()
{ 
  //
  //	Default dtor.
  //
}
//
//____________________________________________________________________________
Block & Barray::operator [] (size_t index)
{
  //	
  //	Access index-th block with dynamical grow of vector size.
  //
  if (size () <= index) {
    ba_block .resize (index + 1, Block (*this));
    if (index == 1) 	ba_block [1] = Block (1, *this);
  }
  return ba_block [index];
}
//
//============================================================================
Block::Block (Barray & parent) 
  : b_sites		(0),
    b_states		(0),
    b_parent		(parent), 
    b_lft		(0),
    b_rgt		(0),
    b_tensor		(),
    b_tensororder	(),
    b_tensorsub		(),
    b_siteorder		(),
    b_optimized		(),
    b_quantum		(),
    b_quantumsub	(),
    b_action		()
{ 
  //
  //	Default ctor.
  //
}
//
//____________________________________________________________________________
Block::Block (const Block & other)
  : b_sites		(other .b_sites),
    b_states		(other .b_states),
    b_parent		(other .b_parent),
    b_lft		(other .b_lft),
    b_rgt		(other .b_rgt),
    b_tensor		(other .b_tensor),
    b_tensororder	(other .b_tensororder),
    b_tensorsub		(other .b_tensorsub),
    b_siteorder		(other .b_siteorder),
    b_optimized		(other .b_optimized),
    b_quantum		(other .b_quantum),
    b_quantumsub	(other .b_quantumsub),
    b_action		(other .b_action)
{ 
  //
  //	Copy ctor.
  //
}
//
//____________________________________________________________________________
Block::Block (const Block & other, Barray & parent)
  : b_sites		(other .b_sites),
    b_states		(other .b_states),
    b_parent		(parent),
    b_lft		(other .b_lft),
    b_rgt		(other .b_rgt),
    b_tensor		(other .b_tensor),
    b_tensororder	(other .b_tensororder),
    b_tensorsub		(other .b_tensorsub),
    b_siteorder		(other .b_siteorder),
    b_optimized		(other .b_optimized),
    b_quantum		(other .b_quantum),
    b_quantumsub	(other .b_quantumsub),
    b_action		(other .b_action)
{ 
  //
  //	Copy ctor.
  //
}
//
//____________________________________________________________________________
Block::Block (const long &, Barray & parent)
  : b_sites		(1),
    b_states		(0),
    b_parent		(parent),
    b_lft		(1),
    b_rgt		(0),
    b_tensor		(),
    b_tensororder	(),
    b_tensorsub		(),
    b_siteorder		(1,0),
    b_optimized		(),
    b_quantum		(),
    b_quantumsub	(),
    b_action		()
{ 
  //
  //	Special ctor: creates a one site default block
  //
}
//
//____________________________________________________________________________
Block::Block (const Block & lft, const Block & rgt)
  : b_sites		(lft .b_sites + rgt .b_sites),
    b_states		(0),
    b_parent		(lft .b_parent),
    b_lft		(lft .b_sites),
    b_rgt		(rgt .b_sites),
    b_tensor		(),
    b_tensororder	(),
    b_tensorsub		(),
    b_siteorder		(lft .b_sites + rgt .b_sites),
    b_optimized		(),
    b_quantum		(),
    b_quantumsub	(),
    b_action		()
{
  //
  //	Tensor Block constructor
  //
  //	Builds the structure of subspaces as the tensor product of 
  //	component Block's subspaces with a reordering of states collected
  //	according to quantum numbers. 
  //
  //	It is assumed that Block's are stored in a global vector<Block> array
  //	labelled by the number of sites. The fields b_lft and b_rgt 
  //	remember the construction and are used as indexes for the glonal 
  //	array during the recursive build of operators.
  //
  //	Alternatively we could remember the tensor components via Block
  //	pointers or references, but we must be sure that these pointers 
  //	does not change during the scope of the program. The use of a global 
  //	vector<Block> container (from standard c++ library) with internal
  //	does not assure that the content is stored in fixed positions.
  //
  size_t lb, rb, le, re, w, j;
  //
  //	Inherit lattice sites from components
  //
  for (lb = 0; lb < lft .b_sites; lb++) 
    b_siteorder [lb] = lft .b_siteorder [lb];
  for ( ; lb < b_sites; lb++) 
    b_siteorder [lb] = rgt .b_siteorder [lb - b_lft] + b_lft;
  //
  //	Pogressively add tensor product of states from the left and 
  //	right Block.
  //
  //	The order of insertion is important and must be consistent 
  //	with the building of operators as tensor product of left and 
  //	right Block operators. 
  //	We choose to insert states |i1;i2> with a double loop on i1, i2 with 
  //	left index i1 running more quickly. 
  //
  //	The resulting states are not in this order because they are grouped 
  //	and reordered with respect to the composition of quantum numbers:
  //		q(i1,i2) = q(i1) + q(i2)
  //
  //	The reordering of states is remembered in the array b_tensororder
  //	(built in Block::addstates ()).
  //
  //	Adding states defines also the quantum space structure of the 
  //	lattice block. This is recorded in the Qspace b_quantum and
  //	b_tensor which are initially (here) identical. b_tensor
  //	is never changed, while b_quantum can be changed by a redefinition 
  //	of states. A particular Action represents the correspomdding 
  //	transformation (as an action from b_tensor Space to b_quantum Space).
  //	Initially this matrix would be the identity matrix and its definition
  //	is avoided.
  //
  Qspace select;
  le = lft .subspaces ();
  re = rgt .subspaces ();
  double alternate = 1.0;
  if (b_lft % 2) alternate = -1.0;
  for (rb = 1; rb < re; rb++) {
    w = rgt .b_quantum .states (rb);
    for (j = 0; j < w; j++) 
      for (lb = 1; lb < le; lb++) {
	Qnumber q  = compose (lft .number (lb), rgt .number (rb), alternate);
	long    fb = lft .quantumstat (lb) * rgt .quantumstat (rb);
	size_t  states = lft .b_quantum .states (lb);
	add_states (q, states, fb);
	if (lft .quantum () .rgtlnk (lb) == rgt .quantum () .lftlnk (rb))
	  select .add_states (q, states, fb);
      }
  }
  if (select .states () < b_tensor .states ()) {
    b_states  = select .states ();
    b_quantum = select;
    Action & base = action (idop_base, 0);
    base = Action (b_tensor, b_quantum);
    base .block (subspaces () - 1);
    Ablock * bb = base .block ();
    size_t msize = minasize ();
    for (size_t sub = 1; sub < subspaces (); sub++) {
      size_t order = b_quantum .offset (sub);
      Ablock & b = bb [sub - 1];
      b .ab_domain = sub;
      for (size_t rsub = 1; rsub < b_tensor .subspaces (); rsub++)
	if (b_tensor .number (rsub) - b_quantum .number (sub) == 0)
	  b .ab_range = rsub;
      b .ab_roffset = msize;
      b .ab_ioffset = 0;
      w = b_quantum .states (sub);
      msize += w * w;
      for (j = 0; j < w; j++)	b_quantumsub [order++] = sub;
    }
    base .size (msize);
    base .storage (msize);
    double * m = base .storage ();    
    for (size_t sub = 1; sub < subspaces (); sub++) {
      size_t off = bb [sub - 1] .ab_roffset;
      w	= b_quantum .states (sub);
      for (j = 0; j < w; j++) m [off + j + j * w] = 1.0;
    }
    base .scalar (1.0);
    base .statistic (1);
  }
}
//
//____________________________________________________________________________
Block::~Block ()
{ 
  //
  //	Default dtor.
  //
}
//
//____________________________________________________________________________
const Action & Block::action (size_t id, size_t index) const
{
  return (b_action [id]) [index];
}
//
//____________________________________________________________________________
Action & Block::action (size_t id, size_t index)
{
  //
  //	Return the action with name id and label index
  //	The action is not built
  //
  if (b_action .size () <= id) b_action .resize (id + 1);
  return (b_action [id]) [index];
}
//
//____________________________________________________________________________
Action & Block::action (const Amono & mono)
{
  //
  //	Return an Action representing the monomial mono
  //
  size_t idmono = 0;
  size_t index  = 0;
  if ((mono .order () == 1) && (mono [0] .af_dg == 0)) {
    idmono = mono [0] .af_op;
    index  = mono [0] .af_st;
  }
  else   idmono = name_action (mono .str ());
  return action (mono, action (idmono, index));
}
//
//____________________________________________________________________________
Action & Block::action (const Amono & mono, Action & a)
{
  //
  //	Build (if needed) Action a representing monomial mono
  //	(only Afactor parts of the monomial, coefficient is ignored)
  //
  if (a .storage ()) {
    if (sites () > 1) a .release ();
    return a;
  }
  //
  //	Action must be built
  //
  if (mono .order () == 0) {
    //
    //	No Afactor: Identity	
    //
    a = Action (b_quantum,1.0);
    return a;
  }
  if (b_sites == 1) {
    //
    //	If only one lattice site all Afactor's must be defined
    //	from input.
    //
    //	Initialize monomial to Identity
    //
    a = Action (b_quantum, 1.0);
    for (size_t n = 0; n < mono .order (); n++) {
      Action factor = action (mono [n] .af_op, mono [n] .af_st);
      if (factor .storage () == 0) {
	cout << "Operator " << name_action (mono [n] .af_op) 
	     << " undefined!" << endl;
	exit (0);
      }
      if (mono [n] .af_dg)  factor .dagger ();
      a *= factor;      
    }
    return a;
  }
  if ((mono .order () == 1) && (mono [0] .af_dg)) {
    //
    //	Avoid recursion on adjoint of single Afactor
    //
    Afactor af (mono [0]);
    af .af_dg = 0;
    a = action (Amono (af));
    a .dagger ();
    a .release ();
    return a;
  }
  //
  //	The monomial is built as a tensor product of operators 
  //	acting on subblocks parts
  //
  Amono mlft (1.0);
  Amono mrgt (1.0);
  long  split = b_lft;
  for (size_t n = 0; n < mono .order (); n++) {
    Afactor af (mono [n]); 
    if (af .af_st < split) mlft *= af;
    else {
      af .af_st -= split;
      mrgt *= af;
    }
  }
  Action lft = b_parent [b_lft] .action (mlft);
  Action rgt = b_parent [b_rgt] .action (mrgt);
  tensor (a, lft, rgt);
  if (b_lft > 1) lft .release ();
  if (b_rgt > 1) rgt .release ();
  return a;
}
//
//____________________________________________________________________________
Action & Block::action (const Apoli & poli, Action & a)
{
  //
  //	Build operator representing polinomial poli
  //
  //	The monomials are reordered according to internal ordering of sites.
  //	This is in general effective only for universe component of superblock.
  //
  Apoli polinomial;
  for (size_t m = 0; m < poli .size (); m++) {
    Amono mono (poli [m]);
    for (size_t n = 0; n < mono .order (); n++) {
      long & st = mono .am_factor [n] .af_st;
      st = b_siteorder [st];
    }
    mono .reorder (b_sites);
    polinomial += mono;
  }
  //
  //  Initialize requested Action
  //
  a = Action (b_quantum, b_quantum);
  for (size_t m = 0; m < polinomial .size (); m++) {
    const Amono & mono = polinomial [m];
    complex<double> z = mono .am_coeff;
    //
    //	If mono is the hermitian of a preceding operator 
    //	fetch the preceding operator and operate the hermitian
    //	transform.
    //
    Amono monodagger (1.0);
    Action b;
    for (size_t n = mono .order (); n; n--) {
      //
      //  Accumulate hermitian factors (in inverse order)
      //
      Afactor af = mono [n-1];
      af .af_dg  ^= 1;
      monodagger *= af;
    }
    monodagger .reorder ();
    if (monodagger < mono) {
      b = action (monodagger, b); 
      b .dagger ();
      z *= monodagger .am_coeff;
    }
    else b = action (mono, b);
    a .add (z, b);
  }
  a .release ();
  return a;
}
//
//____________________________________________________________________________
void Block::actionclear (size_t begin, size_t end)
{
  //
  //	Clear Action's with id's in [begin to end)
  //
  if (b_action .size () < end) end = b_action .size ();
  for (size_t n = begin; n < end; n++)  b_action [n] .clear ();
}
//
//____________________________________________________________________________
void Block::add_states	(const Qnumber & q, size_t states, long fb)
{
  //
  //	Add to block a set of states with common quantum numbers q.
  //	They are insertedand ordered according to quantum numbers, so we 
  //	must remember the insertion order.
  //
  size_t order, i;
  size_t oldstates = b_tensor .states ();
  b_states	   = oldstates + states;
  b_tensororder .resize (b_states);
  b_tensorsub   .resize (b_states);
  //
  //	Add states to b_tensor quantum space and get insertion position
  //
  order = b_tensor .add_states (q, states, fb);
  //
  //	Remember insertion ordering to allow a correct building of 
  //	Block actions from actions on Block sublattices.
  //	The ordering is remembered in the vector b_tensororder and 
  //	the corresponding subspace in the vector b_tensorsub:
  //
  //		b_tensororder [entry order] = new order
  //		b_tensorsub   [new order]   = b_tensor subspace index
  //
  for (i=0; i < oldstates; i++) 
    if (b_tensororder [i] >= order) b_tensororder [i] += states;
  for (i = oldstates; i < b_states; i++, order++) b_tensororder [i] = order;
  //
  //	Compute subspace for every state
  //
  for (size_t sub = 1; sub < b_tensor .subspaces (); sub++) {
    order = b_tensor .offset (sub);
    for (i = 0; i < b_tensor .states (sub); i++) b_tensorsub [order++] = sub;
  }
  //
  //	b_tensor is never defined or changed outside here, while b_quantum
  //	can change according to base transformations
  //
  b_quantum    = b_tensor;
  b_quantumsub = b_tensorsub;
}
//
//____________________________________________________________________________
Action &  Block::base ()
{
  // 
  // 	Returns (and defines) inner Action representing Block's basis states
  //
  return base (action (idop_base, 0));
}
//
//____________________________________________________________________________
Action &  Block::base (Action & b)
{
  //
  //	Returns (and builds) Action representing Block's basis states
  //
  if (action (idop_base, 0) .storage ()) {
    //
    //	Basis states were selected from some density matrix, 
    //	grab the corresponding matrix
    //	
    b << action (idop_base, 0);
    return b;
  }
  else {
    //
    //	If base states are not defined spaces b_quantum and b_tensor
    //	are the same and states are implicitly represented by identity 
    //	matrix
    //
    b = Action (b_tensor, 1.0);
    b .scalaridentity (1.0);
    return b;
  }
}
//
//____________________________________________________________________________
Action & Block::basereflected ()
{
  //
  //	Returns (and builds) Action represent
  //
  return basereflected (action (idop_base, 1));
}
//
//____________________________________________________________________________
Action & Block::basereflected (Action & reflected) 
{
  //
  //	Builds Action representing reflected b_tensor states in the 
  //	b_quantum base
  //
  //		R [l p;j] = <p l |j>	(l,p) states in b_tensor, 
  //					j = state in b_quantum
  //
  if (action (idop_base, 1) .storage ()) {
    //
    //	Previously built:	grab previous
    //
    reflected << action (idop_base,1);
    return reflected;
  }
  //
  //	Action must be built
  //
  if (b_sites == 1) {
    //
    //	For 1-site Block we don't have an effective tensor space
    //	The reflection is not effective and the matrix is simply the 
    //	identity matrix 
    //	
    reflected = Action (b_tensor, 1.0);
    reflected .scalaridentity (1.0);
    return reflected;
  }
  //
  //	The matrix is recursively built from left Block component.
  //
  //	For k sites we have (sums are over repeated indexes):
  //
  //	|l_k> = \sum |l_{k-1} p_k> <l_{k-1} p_k |l_k>
  //
  //	<p_k l_{k-1} | l_k> = 
  //	  \sum <l_{k-1}|l_{k_2} p> <p_k l_{k_2}| l_{k-1}'> <l_{k-1}' p| l_k>
  //
  //	Id est:
  //
  //	R_k [l_{k-1} p_k; l_k]  = \sum U_{k-1}^dagger [l_{k-1}; l_{k-2} p] 
  //	         * R_{k-1} [l_{k-2} p_k; m_{k-1}] * U_k [m_{k-1} p; l_k]
  //
  //	With:
  //		U_k [l_{k-1} p; l_k] = <l_{k-1} p | l_k>
  //
  Block  & l  = b_parent [b_lft];
  Block  & p  = b_parent [1];
  size_t lstates = l .states ();
  size_t pstates = p .states ();
  //
  Action ul, rl, ur;
  ul = l .base  ();
  ul .dagger ();
  rl = l .basereflected ();
  ur = base ();
  //
  //	Setup reflected Ablocks
  //
  const Space & lft = ul .range  ();
  const Space & rgt = ur .domain ();
  //
  size_t  nl  = b_tensor  .subspaces ();
  size_t  nr  = b_quantum .subspaces ();
  Storage bstorage (nl * nr * sizeof (Ablock));
  Ablock * b = (Ablock *) bstorage .storage ();
  size_t rsize, nblocks, n;
  rsize = minasize ();
  //
  //	contraction and pseudo contracion
  //
  Action id (lft, 1.0);
  id .scalaridentity (1.0);
  //
  //	rc [q]   [l'; r] = <l' q| r>
  //	lc [q p] [l; l'] = <l| m q> <p m| l'> (summed over m)
  //
  Aarray rc, lc;
  contraction    (id, 0, pstates, ur, rc);
  l .contraction (ul, pstates, pstates, rl, lc);
  //
  for (size_t p = 0; p < pstates; p++) {
    Action r (lft, rgt);
    for (size_t q = 0; q < pstates; q++) {
      lc [q + p * pstates] *= rc [q];
      r  += lc [q + p * pstates];
    }
    //
    // r [i; j] = <p i| j> = R_k [i p; j]
    //
    Ablock * br = r .block ();
    for (n = 0; n < r .blocks (); n++) {
      size_t i  = br [n] .ab_range;
      size_t j  = br [n] .ab_domain;
      size_t ro = br [n] .ab_roffset;
      size_t io = br [n] .ab_ioffset;
      //
      i  = lft .offset      (i);
      i  = b_tensororder    [i + p * lstates];
      i  = b_tensorsub      [i];
      //
      Ablock & bb = b [i + j * nl];
      bb .ab_range  = i; 
      bb .ab_domain = j;
      size_t ms = b_tensor .states (i) * b_quantum .states (j);
      if ((bb .ab_roffset == 0) && ro) {
	bb .ab_roffset = rsize;
	rsize += ms;
      }
      if ((bb .ab_ioffset == 0) && io) {
	bb .ab_ioffset = rsize;
	rsize += ms;
      }
     }
    lc [p] << r;
   }
  //
  for (n = 0, nblocks = 0; n < nl * nr; n++) 
    if (b [n] .ab_roffset || b [n] .ab_ioffset) nblocks++;
  //
  //	Build resulting Action
  //
  reflected = Action (b_tensor, b_quantum);
  reflected .size    (rsize);
  reflected .storage (rsize);
  reflected .blocks  (nblocks);
  reflected .block   (nblocks);
  //
  double * mr = reflected .storage ();
  Ablock * br = reflected .block   ();
  for (n = 0, nblocks = 0; n < nl * nr; n++) 
    if (b [n] .ab_roffset || b [n] .ab_ioffset) br [nblocks++] = b [n];
  //
  for (size_t p = 0; p < pstates; p++) {
    Ablock * ba = lc [p] .block   ();
    double * ma = lc [p] .storage ();
    for (n = 0; n < lc [p] .blocks (); n++) {
      size_t ra = ba [n] .ab_roffset;
      size_t ia = ba [n] .ab_ioffset;
      size_t ls = ba [n] .ab_range;
      size_t rs = ba [n] .ab_domain;
      size_t w  = rgt .states (rs);
      size_t lh = lft .states (ls);
      size_t lo = lft .offset (ls);
      lo  = b_tensororder [lo + p * lstates];
      ls  = b_tensorsub   [lo];
      lo -= b_tensor .offset (ls);
      size_t h    = b_tensor .states (ls);
      Ablock & bb = b [ls + rs * nl];
      size_t roff = bb .ab_roffset + lo;
      size_t ioff = bb .ab_ioffset + lo;
      if (ra) 
	for (size_t j = 0; j < w; j++, roff += h)
	  for (size_t i = 0; i < lh; i++, ra++) mr [roff + i] += ma [ra];
      if (ia) 
	for (size_t j = 0; j < w; j++, ioff += h)
	  for (size_t i = 0; i < lh; i++, ia++) mr [ioff + i] += ma [ia];
    }
  }
  //	
  //	Check statistic (like identity)
  //
  reflected .statistic (1);
  reflected .scalar (1.0);
  return reflected;
}
//
//____________________________________________________________________________
void Block::contraction (const Action & lft, size_t pstates, 
			 size_t qstates, const Action & rgt,
			 Aarray & result, bool fermi) const
{
  //
  //	Operate the partial matrix product contracting left component
  //	subspspace of b_tensor:
  //
  //	result [p q] [i; j] = \sum_k lft [i;k p] * rgt [k q; j] 
  //
  //	taking into account the exchange of k and p if fermi is true
  //	
  //	The result is an Aarray of Action's indexed by 1-site states 
  //	p q.
  //	
  //	pstates is the number of states p, qstates is the number of
  //	q states. If pstates (qstates) is zero or 1, the k p (k q) are 
  //	the states of left subspace of b_tensor
  //
  if (pstates == 0) pstates = 1;
  if (qstates == 0) qstates = 1;
  //
  //	Special cases
  //
  if ((pstates == 1) && (qstates == 1)) {
    //
    //	Contraction is a simple matrix product
    //
    result [0]  = lft;
    result [0] *= rgt;
    return;
  }
  //	
  //	lft's range and rgt's domain are arbitrary
  //
  const Space & arange  = lft .range  ();
    const Space & adomain = rgt .domain ();
  //
  size_t nr = arange  .subspaces ();
  size_t nd = adomain .subspaces ();
  size_t nb = nr * nd;
  //
  if (b_sites == 1) {
    //
    //	No contraction, simple entries product 
    //
    if ((pstates != b_states) || (qstates != b_states)) {
      cout << "Invalid contraction (1) " 
	   << pstates << " != " << b_states 
	   << " || " << b_states << " != " << qstates << endl;
      exit (0);
    }
    //
    //	Matrices are small, no special trick used, expand all and 
    //	then compress results.
    //
    size_t h = arange  .states ();
    size_t w = adomain .states ();
    Storage lfull (h * pstates * sizeof (complex<double>));
    Storage rfull (w * qstates * sizeof (complex<double>));
    complex<double> * mlft = (complex<double> *) lfull .storage ();
    complex<double> * mrgt = (complex<double> *) rfull .storage ();
    lft .expand (mlft);
    rgt .expand (mrgt);
    for (size_t q = 0; q < qstates; q++) 
      for (size_t p = 0; p < pstates; p++) {
	size_t pq = p + q * pstates;
	Storage full (h * w * sizeof (complex<double>));
	complex<double> * m = (complex<double> *) full .storage ();
	for (size_t j = 0; j < w; j++) 
	  for (size_t i = 0; i < h; i++)
	    m [i + j * h] = mlft [i + p * h] * mrgt [q + j * qstates];
	result [pq] = Action  (arange, adomain);
	result [pq] .compress (m);
      }
    return;
  }
  //
  //	Some simple consistency check
  //	
  const Space & sub = b_parent [b_lft] .quantum ();
  const Space & psp = b_parent [1]     .quantum ();
  size_t  ss  = sub .states ();
  //
  if ((pstates * ss != lft .width  ()) ||
      (qstates * ss != rgt .height ())) {
    cout << "Invalid contraction (" << b_sites << "=" << b_lft 
	 << "+" << b_rgt <<") " 
	 << pstates << "x" << ss << " != " << lft .width  () << " || "
	 << rgt .height () << " != " << ss << "x" << qstates << endl;
    exit (0);
  }
  //
  //	Ablock structures
  //
  size_t ns = sub .subspaces ();
  Storage bstorage (nb * sizeof (Ablock));
  Ablock * blk  = (Ablock *) bstorage .storage ();
  Ablock * blft = lft .block ();
  Ablock * brgt = rgt .block ();
  //
  double * mlft = lft .storage ();
  double * mrgt = rgt .storage ();
  //
  for (size_t q = 0; q < qstates; q++) 
    for (size_t p = 0; p < pstates; p++) {
      //
      //   compute statistic of p
      //
      size_t psub  = b_parent [1] .quantumsub (p);
      long   pstat = psp .statistic (psub);
      size_t pq = p + q * pstates;
      size_t n, rsize, nblocks, bs, k, kp, kq, bl, br, ar, ad, h, w;
      //
      // Clear Ablock structure
      //
      for (n = 0; n < nb; n++) blk [n] = Ablock ();
      rsize = minasize ();
      //
      //  loop on contraction subspaces
      //
      for (bs = 1; bs < ns; bs++) {
	kp = kq = bs;
	k = sub .offset (bs);
	if (pstates > 1) kp = b_tensorsub [b_tensororder [k + p * ss]];
	if (qstates > 1) kq = b_tensorsub [b_tensororder [k + q * ss]];
	//
	//  Look for submatrices to contract
	//
	for (br = 0; br < rgt .blocks (); br++) {
	  if (brgt [br] .ab_range != kq) continue;
	  ad = brgt [br] .ab_domain;
	  w  = adomain .states (ad);
	  for (bl = 0; bl < lft .blocks (); bl++) {
	    if (blft [bl] .ab_domain != kp) continue;
	    ar = blft [bl] .ab_range;
	    h  = arange .states (ar);
	    Ablock & b = blk [ar + ad * nr];
	    b .ab_range  = ar;
	    b .ab_domain = ad;
	    if ((b .ab_roffset == 0) &&
		((blft [bl] .ab_roffset && brgt [br] .ab_roffset) ||
		 (blft [bl] .ab_ioffset && brgt [br] .ab_ioffset))) {
	      b .ab_roffset = rsize;
	      rsize += h * w;
	    }
	    if ((b .ab_ioffset == 0) &&
		((blft [bl] .ab_roffset && brgt [br] .ab_ioffset) ||
		 (blft [bl] .ab_ioffset && brgt [br] .ab_roffset))) {
	      b .ab_ioffset = rsize;
	      rsize += h * w;
	    }
	  }
	}
      }
      //
      //	Count blocks and compute statistic
      //
      long pqstat = 1;
      for (n = 0, nblocks = 0; n < nb; n++) 
	if (blk [n] .ab_roffset || blk [n] .ab_ioffset) {
	  long astat = arange .statistic (blk [n] .ab_range) *
	      adomain .statistic (blk [n] .ab_domain);
	  if (nblocks == 0) pqstat = astat;
	  if (astat != pqstat) pqstat = 0;
	  nblocks++;
	}
      //
      //	Define result
      //
      result [pq] = Action (arange, adomain);
      result [pq] .size    (rsize);
      result [pq] .storage (rsize);
      result [pq] .blocks  (nblocks);
      result [pq] .block   (nblocks);
      double * ma  = result [pq] .storage ();
      Ablock * ba  = result [pq] .block   ();
      for (n = 0, nblocks = 0; n < nb; n++) 
	if (blk [n] .ab_roffset || blk [n] .ab_ioffset) 
	  ba [nblocks++] = blk [n];   
      //
      //	Fill memory
      //
      for (bs = 1; bs < ns; bs++) {
	kp = kq = bs;
	k  = sub .offset (bs);
	size_t nk     = sub .states (bs);
	size_t stride = nk;
	//
	size_t poff = 0;
	size_t qoff = 0;
	double factor = 1.0;
	long kstat = sub .statistic (bs);
	if (fermi && (pstat < 0) && (kstat < 0)) factor = -1.0; 
	if (pstates > 1) {
	  poff  = b_tensororder [k + p * ss];
	  kp    = b_tensorsub   [poff];
	  poff -= b_tensor .offset (kp);
	}
	if (qstates > 1) {
	  qoff   = b_tensororder [k + q * ss];
	  kq     = b_tensorsub   [qoff];
	  qoff  -= b_tensor .offset (kq);
	  stride = b_tensor .states (kq);
	}
	//
	//  Look for submatrices to contract
	//
	for (br = 0; br < rgt .blocks (); br++) {
	  if (brgt [br] .ab_range != kq) continue;
	  ad = brgt [br] .ab_domain;
	  w  = adomain .states (ad);
	  size_t rgtro  = brgt [br] .ab_roffset;
	  size_t rgtio  = brgt [br] .ab_ioffset;
	  for (bl = 0; bl < lft .blocks (); bl++) {
	    if (blft [bl] .ab_domain != kp) continue;
	    ar = blft [bl] .ab_range;
	    h  = arange .states (ar);
	    size_t lftro = blft [bl] .ab_roffset;
	    size_t lftio = blft [bl] .ab_ioffset;
	    Ablock & b   = blk [ar + ad * nr];
	    size_t roff  = b .ab_roffset;
	    size_t ioff  = b .ab_ioffset;
	    //
	    //	real * real 
	    //
	    if (roff && lftro && rgtro) 
	      for (size_t l = 0; l < nk; l++) {
		size_t lp = lftro + (poff + l) * h;
		size_t lq = rgtro +  qoff + l;
		for (size_t j = 0; j < w; j++) {
		  double r = mrgt [lq + j * stride] * factor;
		  for (size_t i = 0; i < h; i++) 
		    ma [roff + i + j * h] += mlft [i + lp] * r;
		}
	      }
	    //
	    //	imag * imag 
	    //
	    if (roff && lftio && rgtio) 
	      for (size_t l = 0; l < nk; l++) {
		size_t lp = lftio + (poff + l) * h;
		size_t lq = rgtio +  qoff + l;
		for (size_t j = 0; j < w; j++) {
		  double r = mrgt [lq + j * stride] * factor;
		  for (size_t i = 0; i < h; i++) 
		    ma [roff + i + j * h] -= mlft [i + lp] * r;
		}
	      }
	    //
	    //	real * imag
	    //
	    if (ioff && lftro && rgtio) 
	      for (size_t l = 0; l < nk; l++) {
		size_t lp = lftro + (poff + l) * h;
		size_t lq = rgtio +  qoff + l;
		for (size_t j = 0; j < w; j++) {
		  double r = mrgt [lq + j * stride] * factor;
		  for (size_t i = 0; i < h; i++) 
		    ma [ioff + i + j * h] += mlft [i + lp] * r;
		}
	      }
	    //
	    //	imag * real 
	    //
	    if (ioff && lftio && rgtro) 
	      for (size_t l = 0; l < nk; l++) {
		size_t lp = lftio + (poff + l) * h;
		size_t lq = rgtro +  qoff + l;
		for (size_t j = 0; j < w; j++) {
		  double r = mrgt [lq + j * stride] * factor;
		  for (size_t i = 0; i < h; i++) 
		    ma [ioff + i + j * h] += mlft [i + lp] * r;
		}
	      }
	  }
	}
      }
      if (rsize > minasize ()) result [pq] .scalar (1.0);
      else result [pq] .scalar (0.0);
      //
      //  Set statistic
      //
      result [pq] .statistic (pqstat);
    }
  return;
}
//
//____________________________________________________________________________
Block & Block::operator = (const Block & other) 
{
  //
  //	Assignment operator
  //	Note that operator content of other is removed after transfer
  //	into this block
  //
  new (this) Block (other);
  return *this;
}
//
//____________________________________________________________________________
Block & Block::operator << (const Block & other) 
{
  //
  //	grab operator
  //	The structure of the other block is copied to this block
  //	with grabbing of operator content.
  //
  b_sites   = other .b_sites;
  b_states  = other .b_states;
  //
  // b_parent is retained 
  //
  b_lft         = other .b_lft;
  b_rgt         = other .b_rgt;
  b_tensor      = other .b_tensor;
  b_tensororder = other .b_tensororder;
  b_tensorsub   = other .b_tensorsub;
  b_siteorder   = other. b_siteorder;
  b_optimized   = other .b_optimized;
  b_quantum     = other .b_quantum;
  b_quantumsub  = other .b_quantumsub;
  //
  //	grab operator content (we don't want duplicate actions memory)
  //
  for (size_t opid = 0; opid < other .b_action .size (); opid++) {
    for (size_t index = 0; index < other .b_action [opid] .size (); index++) {
      if (other .action (opid,index) .storage ()) {
	const Action & origin = other .action (opid, index);
	Action & destination = action (opid, index);
	destination << origin;
      }
    }
  }
  //
  return *this;
}
//
//____________________________________________________________________________
void Block::optimize (const Apoli & poli)
{
  // 
  // 	Try to optimize basis states to minimize interaction part
  //	matrix elements of poli operator
  //
  Apoli polinomial;
  bool  left = false, right = false;
  for (size_t m = 0; m < poli .size (); m++) {
    Amono mono (poli [m]);
    for (size_t n = 0; n < mono .order (); n++) {
      long & st = mono .am_factor [n] .af_st;
      st = b_siteorder [st];
      if (st < (long) b_lft) 	left = true;
      else	  		right = true;
    }
    mono .reorder (b_sites);
    polinomial += mono;
  }
  if (left && right) return;
  if (left && (b_lft < b_rgt)) return;
  if (right && (b_rgt < b_lft)) return;
  if (left) {
    if (b_parent [b_lft] .action (idop_base) .storage () == 0) return;
    Action mlft;
    polinomial .show ("Optimize left ", 4);
    b_parent [b_lft] .action   (polinomial, mlft);
    b_parent [b_lft] .optimize (mlft); 
  }
  if (right) {
    if (b_parent [b_rgt] .action (idop_base) .storage () == 0) return;
    Action mrgt;
    polinomial .show ("Optimize right", 4);
    b_parent [b_rgt] .action   (polinomial, mrgt);
    b_parent [b_rgt] .optimize (mrgt);
  }
}
//
//____________________________________________________________________________
void Block::optimize (Action & a)
{
  //
  //	Try to optimize basis states to minimize matrix elements 
  //	of given action
  //
  size_t   nb, ir, jd, ro, io, w, h, nm;
  Ablock * b = a .block   ();
  double * m = a .storage ();
  for (nb = 0; nb < a .blocks (); nb++) {
    ir = b [nb] .ab_range;
    jd = b [nb] .ab_domain;
    ro = b [nb] .ab_roffset;
    io = b [nb] .ab_ioffset;
    h  = b_quantum .states (ir);
    w  = b_quantum .states (jd);
    nm = 0;
    for (size_t n = 0; n < h*w; n++) {
      double mc = 0.0;
      if (ro) mc += m [ro+n] * m [ro+n];
      if (io) mc += m [io+n] * m [io+n];
      if (mc > 1.e-12) nm++;
    }
    cout << "(" << ir << " " << h << " " << b_optimized [ir] << ","
	 << "(" << jd << " " << w << " " << b_optimized [jd] << ") "
	 << nm << "/" << h*w << endl;
    b_optimized [ir] = b_optimized [jd] = 1;
  }
}
//
//____________________________________________________________________________
void Block::reflect ()
{
  //
  //	Globally reflects the order of sites representing the lattice
  //
  vector<size_t> reflected (b_siteorder);
  for (size_t i = 0; i < b_sites; i++)
    reflected [i] = b_siteorder [b_sites - 1 - i];
  b_siteorder = reflected;
  //if (b_sites % 2 == 0)	b_quantum .reflect ();
}
//
//____________________________________________________________________________
void Block::reflectlft ()
{
  //
  //	Locally reflects the order of sites representing left lattice
  //
  vector<size_t> reflected (b_siteorder);
  for (size_t i = 0; i < b_lft; i++)
    reflected [i] = b_siteorder [b_lft -1 - i];
  b_siteorder = reflected;  
}
//
//____________________________________________________________________________
void Block::reflectreset ()
{
  //
  //	Reset to natural order of sites
  //
  vector<size_t> reflected (b_siteorder);
  for (size_t i = 0; i < b_sites; i++)
    reflected [i] = i;
  b_siteorder = reflected;  
}
//
//____________________________________________________________________________
void Block::reflectrgt ()
{
  //
  //	Locally reflects the order of sites representing right lattice
  //
  vector<size_t> reflected (b_siteorder);
  for (size_t i = b_lft; i < b_sites; i++)
    reflected [i] = b_siteorder [b_lft + b_sites - 1 - i];
  b_siteorder = reflected;  
}
//
//____________________________________________________________________________
void Block::release ()
{
  //
  //  release Action memory
  //
  for (size_t n = 0; n < b_action .size (); n++) b_action [n] .release ();
}
//
//______________________________________________________________________SORTING
template <typename T>
void quickSort(T *array, size_t left, size_t right)
{
    size_t l = left;
    size_t r = right - 1;
    size_t size = right - left;

    if (size > 1) {
        T pivot = array[rand() % size + l];

        while (l < r) {
            while (array[r] > pivot && r > l) {
                r--;
            }

            while (array[l] < pivot && l <= r) {
                l++;
            }

            if (l < r) {
                std::swap(array[l], array[r]);
                l++;
            }
        }

        quickSort(array, left, l);
        quickSort(array, r, right);
    }
}
//
//____________________________________________________________________________
void Block::select_states (const Action & density, size_t min, size_t max, 
			   double n_cut)
{
  //
  //	Select states for the block as the eigenvectors of given
  //	density matrix corresponding to highest eigenvalues.
  //
  //	The number of choiced states is between the limits min and max
  //	(if it is possible). The number is chosen trying to minimize
  //	symmetry breakings in the basis.  The consecutive eigenvalues
  //	differences are computed and the cut is chosen at the highest 
  //	jump satisfying the limits (I don't know if this is the optimal 
  //	way but it helps).
  //
  if (b_states <= 50) return;
  size_t   i, j, nb, sub, offset, dimension;
  //size_t   extra =  max + density .blocks ();
  //if (extra >= b_states) extra = (max + b_states)/2; 
  //
  //	Do yo want to show density matrix?
  //
  if (b_sites <= show_density) density .show ("Density matrix");
  Ablock * b  = density .block   ();
  double * m  = density .storage ();
  //
  //	Temporary storage area
  //
  size_t ss = 2 * b_states * sizeof (double);
  size_t sw = 2 * b_tensor .subspaces () * sizeof (long);
  size_t eig_ord = b_states * sizeof (double);
  Storage  mem (ss);
  Storage  mw  (sw);
  Storage  m_ord (eig_ord);
  double * eigen  = (double *) mem .storage ();
  double * offd   = eigen  + b_states;
  //
  //    We add the possibility of choice the maximum error insted the 
  //    max and min number of states. eig_ord is the vector where we 
  //    order the eigenvalues.
  //
  double * eigord = (double *) m_ord .storage ();
  long   * weight = (long *)   mw  .storage ();
  long   * occ    = weight + b_tensor .subspaces ();
  for (i = 0; i < 2 * b_states; ++i) eigen [i] = 0.0;
  //
  for (nb = 0; nb < density .blocks (); nb++) {
    //	
    //	A simple check if density is block diagonal
    //
    // cout  << "BLOCK " << nb << " " << b [nb] .ab_domain << " "
	//   << b [nb] .ab_range << endl;
    if (b [nb] .ab_domain != b [nb] .ab_range) {
      cout << "Not block diagonal density matrix! " 
	   << "block " << nb << " " << b [nb] .ab_domain << " "
	   << b [nb] .ab_range << endl;
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
    householder (mr, mi, eigen + offset, offd + offset, dimension, dimension);
    if (tqli (eigen + offset, offd + offset, mr, mi, dimension, dimension)) {
      cout << "Block::select_states (): tqli error!" << endl;
      exit (0);
    }
    reorder_high_low (eigen + offset, offd + offset, mr, mi, 
		      dimension, dimension);
  }
  //
  //	Compute jumps
  //
  double eps = machine_precision ();
  size_t selected = 0;
  double emax = -2.0;
  double trace = 0.0;
  for (i = 0; i < b_states; ++i)  {
    //
    //	cleanup 
    //
    if (eigen [i] < eps*1.e-6) eigen [i] = 0.0;
    trace += eigen [i];
  }
  for (i = 0, j = b_states; i < b_states; ++i) {
    //
    //	normalize and find highest eigenvalues
    //
    eigen [i] /= trace;
    //	cout << "vecchi:   " << eigen[i] << endl;    
    offd [i] = 0.0;
	eigord[i]=eigen[i];
    if (eigen [i] > emax) {
      emax = eigen [i];
      j = i;
    }
  }

    //	Use offd as a mark for eigenvalues selected
    //
    if (j < b_states) offd [j] = ++selected;
    //
    
  //
  //        Here we order the eigenvalues
  //

/*    for (int i1=b_states-1; i1>=0; i1--) {
        for (int j1=i1-1; j1>=0 ; j1--) {
            if (eigord[j1+1] > eigord[j1]){
               double temp = eigord[j1+1];
               eigord[j1+1]=eigord[j1];
               eigord[j1]=temp;
            }
        }
    }*/

//
//	SORTING
//

    quickSort(eigord, 0, b_states);
/*		for (i = 0; i<b_states; ++i) 
			{
	cout << "autovalore" << "   " << eigord[i] << "   " << i << endl;		
			}*/
   
    //
	//		Set variable number of states  
	//
	
if (n_cut > 0 && n_cut < 1){

    //    Sum of the eigenvalues
    
    double sum = 0;
    int cut_states_number = 0;    
    for ( i = b_states-1; i > 0; i--) {
      sum += eigord[i];   
	  cut_states_number += 1;
//              cout << "Ciao gattone = " << 1-n_cut << "  " << sum << "  " << eigord[i] << "  " << i << endl;   
      if(sum > 1 - n_cut) {
		break;
      }
    }

        
	if (cut_states_number > min && cut_states_number < max){    
        min = cut_states_number;
        max = cut_states_number;
 	}            
    
    if (cut_states_number < min) max = min ;
    if (cut_states_number > max) min = max ;
    
    
//        cout << "n_cut = " << setprecision(9) << n_cut << endl;
    cout << "States -> Cut / Min / Max = " << cut_states_number << "  " << min << "  " << max << endl;
//        cout << "precision = " << setprecision(9) << sum << endl;
};

	//	Print Entanglement SPECTRUM3

 //   for ( i=0; i<15; ++i) {
//   	cout << "autovalore  " << i << "  " << eigord[i] << endl;
 //   }

  double discard = -1.;
  double jump  = -1.0;
  double vcut  = -1.0;
  double vlast = -1.0;
  size_t cut  = b_states;
  size_t minselected = selected;
  size_t maxselected = 0;
  while (selected < b_states) {
    double dmax = -1;
    for (i = 0, j = b_states; i < b_states; ++i) {
      if (offd  [i] > 0.0) continue;
      if (eigen [i] > dmax) {
	dmax = eigen [i];
	j = i;
      }
    }
    double delta = emax - dmax;
    if (delta < 8.0 * eps) delta = 0.0;
    if (delta > 0.0) minselected = selected;
    if ((delta > jump) && (selected >= min)) {
      //
      //  Here we have a possible candidate for a cut on eigenvalues
      //  (We want to select the highest jump between consecutive 
      //  eigenvalues, with a cut position inside the requested 
      //  range or a little above).
      //  
      jump = delta;
      vlast = dmax + 8.0 * eps;
      maxselected = selected;
      discard = dmax;
      //
      //  Forced cut 
      //
      if (min == max) {
	vcut = vlast;
	break;
      }
      //
      //  If the selection is inside the range (borders excluded)
      //  we record the position (updating it for highest jump)
      //
      if ((min < selected) && (selected < max)) vcut = vlast;
      else if ((selected == max) && (delta > eps)) {
	//
	//  The jump is a true high jump at the upper border. Accept
	//  it (and skip other jumps).
	//
	vcut = vlast;
	break;
      }
      else if (selected > max) {
	//
	//  The jump is in a position out of borders. Accept if last
	//  jump was exactly on the borders (this is signaled by the
	//  initialization vcut < 0) and the jump is no too far away
	//  from the borders. Skip any other jump.
	//
	if ((vcut < 0.0) && (selected < b_states))  vcut = vlast;
	break;	
      }
    }
    if (delta > 0.0) emax = dmax;
    if (j < b_states) offd [j] = ++selected;
  }
  if (vcut < 0.0) vcut = vlast;
  if (discard < 0.0) discard = 0.0;
  //
  //	vcut is the cutting value for eigenvalues. 
  //	All eigenvalues above vcut are accepted, other are discarded
  //	Special care must be taken in the case of null eigenvalues
  //	and null vcut (or near). In this case only non null eigenvalues
  //	are chosen and the number can be lower than the requested minimum.
  //
  //	Use offd to mark chosen eigenvalues
  //
  if (min < max) {
    cut = 0;
    for (i = 0; i < b_states; ++i) {
      offd [i] = 0.0;
      if (eigen [i] > vcut) {
	offd [i] = 1.0;
	cut++;
      }
    }
    minselected = cut;
  }
  if (cut < min) {
    //
    //	Simple comparison with vcut gives a number lower than requested 
    //	minimum. Add values in order to reach minimum balancing (if 
    //	possible) among subspaces. We assign a weight to every subspace 
    //	proportional to its emptyness.
    //
    long qHare = (double (b_states - cut)) / ((double) (min - cut));
    long socc;
    long maxw = 0;
    for (sub = 1; sub < b_tensor .subspaces (); ++sub) {
      offset    = b_tensor .offset (sub);
      dimension = b_tensor .states (sub);
      socc = 0;
      for (i = offset; i < offset + dimension; ++i) if (offd [i]) socc++;
      weight [sub] = dimension - socc;
      if (weight [sub] > maxw) maxw = weight [sub];
      occ    [sub] = socc;
    }
    //
    //	Distribution of missing values like seats in a parliament
    //  (Hare method)
    //
    while (cut < min) {
      //
      socc = 0;
      for (sub = 1; sub < b_tensor .subspaces (); ++sub) {
	dimension = b_tensor .states (sub);
	if ((weight [sub] == maxw) && (occ [sub] < (long) dimension)) {
	  occ    [sub] += 1;
	  weight [sub] -= qHare;
	  cut++;
	  if ((min == max) && (min == cut)) break;
	}
	if (weight [sub] > socc)	socc = weight [sub];
      }
      maxw = socc;
      if ((min == max) && (min == cut)) break;
      if (maxw <= 0) break;
    }
    //
    vcut = 0.0;
    for (sub = 1; sub < b_tensor .subspaces (); ++sub) {
      offset    = b_tensor .offset (sub);
      dimension = b_tensor .states (sub);
      vlast = 0.0;
      for (i = offset, socc = 0; socc < occ [sub]; i++, socc++) 
	offd [i] = 1.0;
      if (socc < (long) dimension) vlast = eigen [offset + socc];
      if (vlast > vcut) vcut = vlast;
    }
  }
  //
  //	Build new space with chosen states
  //
  Qspace select;
  vector<Ablock> blk;
  size_t msize = b [0] .ab_roffset;
  size_t nsub = 0;
  double accept  = 0.0;
  if (show_selection) {
    stringstream out1, out2;
    out1 << "(pre-selection "
	 << right << minselected << "-" 
	 << left << maxselected << "/";
    out2 << left << min << "-" << left << max << ")";
    cout << setw (7) << left << "States "
         << setw (26) << right << out1 .str ()
	 << setw (11) << left << out2.str  ()
	 << setw (6) << right << "trace"
	 << setw (10) << setprecision (6) << right << fixed << trace
	 << setw (5) << right << "cut " 
	 << setw (13) << setprecision (6) << right << scientific 
	 << discard << endl;
  }
  cut = 0;
  discard = 0.0;
  for (sub = 1; sub < b_tensor .subspaces (); ++sub) {
    ss = 0;
    offset    = b_tensor .offset (sub);
    dimension = b_tensor .states (sub);
    double upper = 0.0;
    double lower = 0.0;
    for (i = offset; i < offset + dimension; ++i) 
      if (offd [i]) {
	cut++;
	ss++;
	upper += eigen [i];
      } 
      else lower += eigen [i];
    //
    if (show_states) { 
      cout << " Space " << setw (20) << left 
	   << b_tensor .number (sub) .str () 
	   << setw (5) << right << ss << "/" << setw (5) << left  
	   << dimension;
      if (ss == 0) cout << setw (26) << " ";
      else {
	cout << setw (13) << setprecision (6) << right << scientific
	     << eigen [offset];
	if (ss == 1) cout << setw (13) << " ";
	else cout << setw (13) << setprecision (6) << right << scientific
		  << eigen [offset + ss - 1];
      }
      cout << " ";
      if (ss < dimension) 
	cout << setw (13) << setprecision (6) << right << scientific
	     << eigen [offset + ss]; 	
      else cout << setw (13) << " ";
      cout << endl;
    }
    accept += upper;
    discard += lower;
    //
    if (ss) {
      nsub++;
      select .add_states (b_tensor .number (sub), ss, 
			  b_tensor .statistic (sub));
      for (nb = 0; nb < density .blocks (); nb++) 
	if (b [nb] .ab_range == sub) break;
      Ablock bb;      
      if (nb == density .blocks ()) {
	bb .ab_roffset = msize;
	msize += ss * dimension;
      }
      else {
	if (b [nb] .ab_roffset) {
	  bb .ab_roffset = msize;
	  msize += ss * dimension;
	}
	if (b [nb] .ab_ioffset) {
	  bb .ab_ioffset = msize;
	  msize += ss * dimension;
	}
      }
      bb .ab_range  = sub;
      bb .ab_domain = nsub;
      blk .push_back (bb);
    }
  }
  if (show_selection) 
    cout << setw (27) << left << "Total selected:"  
	 << setw (5) << right << cut << "/" << setw (5) << left  
	 << b_states << setw (26) << right 
	 << "Truncation error:" 
	 << setw (14) << setprecision (6) << right << scientific
	 << discard / trace << endl;
  truncation = discard /trace;
  b_states  = cut;
  b_quantum = select;
  //
  //	Compute subspace for every state
  //
  for (sub = 1; sub < subspaces (); ++sub) {
    size_t order  = b_quantum .offset (sub);
    for (i = 0; i < b_quantum .states (sub); ++i) 
      b_quantumsub [order++] = sub;
  }
  //
  //	Define Action representing selected states;
  //
  Action & matrix = action (idop_base, 0);
  matrix = Action (b_tensor, select);
  matrix .size    (msize);
  matrix .storage (msize);
  matrix .blocks  (blk .size ());
  matrix .block   (blk .size ());
  double * mn = matrix .storage ();
  Ablock * mb = matrix .block   ();
  //
  //	Fill matrix memory and blocks
  //
  for (nsub = 0; nsub < blk .size (); nsub++) {
    mb [nsub] = blk [nsub];
    i  = mb [nsub] .ab_range;
    j  = mb [nsub] .ab_domain;
    ss = b_quantum .states (j) * b_tensor .states (i) * sizeof (double);
    for (nb = 0; nb < density .blocks (); nb++) 
      if (b [nb] .ab_range == i) break;
    if (nb == density .blocks ()) {
      memset (mn + mb [nsub] .ab_roffset, 0, ss);
      dimension = b_tensor .states (i);
      for (size_t k = 0; k < b_quantum .states (j) * dimension; 
	   k += (dimension + 1))  mn [mb [nsub] .ab_roffset + k] = 1.0; 
    }
    else {
      if (b [nb] .ab_roffset) 
	memcpy (mn + mb [nsub] .ab_roffset, m + b [nb] .ab_roffset, ss);
      if (b [nb] .ab_ioffset) 
	memcpy (mn + mb [nsub] .ab_ioffset, m + b [nb] .ab_ioffset, ss);
    }
  }
  matrix .scalar (1.0);
  // Action & reflected = action (idop_base, 1);
  // basereflected (reflected);
  //	
  //	Mark basis states as not optimized
  //
  b_optimized .resize (subspaces ());
  for (nsub = 0; nsub < subspaces (); nsub++) b_optimized [nsub] = 0;
  //
  //	Update defined Action's
  //
  Action mdagger = matrix;
  mdagger .dagger ();
  for (size_t id = 0; id < b_action .size (); id++) {
    if (id == idop_base) continue;
    for (size_t index = 0; index < b_action [id] .size (); index++) {
      Action & op = action (id, index);
      if (op .isscalar ()) {
	complex<double> z = op .scalar ();
	op = Action (b_quantum, z);
      }
      else if (action (id,index) .size ()) {
	Action b = mdagger;
	b *= op;
	b *= matrix;
	op << b;
      }
      op .release ();
    }
  }
}
//
//____________________________________________________________________________
void Block::show (const string & s)
{
  //
  //	Print out relevant Block data
  //
  if (s .size ()) cout << s << ": ";
  cout << "sites=" << b_sites << " (" << b_lft << "+" << b_rgt 
       << "), states=" << b_states << ", subspaces=" 
       << subspaces () - 1 << endl;
  //
  //	Subspaces details
  //
  b_quantum .show (" subspace "); 
  //
  //	Defined actions 
  //
  for (size_t id = 0; id < b_action .size (); id++) {
    string name = name_action (id);
    for (size_t index = 0; index < b_action [id] .size (); index++) {
      Action & op = action (id, index);
      if (op .storage ()) {
	stringstream s;
	s << "action(" << id << ") = " << name;
	if (index || (name [name .size () - 1] != ']')) {
	  if (name [name .size () - 1] != '+') s << " ";
	  s << "[" << index << "]";
	}
	op .show (s .str ());
      }
    }      
  }
}
//
//____________________________________________________________________________
Action & Block::symmetry (size_t id) 
{
  //
  //	Returns (and build) symmetry operator id
  //
  return symmetry (id, action (id, 0));
}
//
//____________________________________________________________________________
Action & Block::symmetry (size_t id, Action & s)
{
  //
  //	Builds (if needed) s as the symmetry Action id
  //
  if (action (id, 0) .storage ()) {
    s = action (id, 0);
    if (sites () > 1) s .release ();
    return s;
  }
  //
  //	symmetry is the tensor product of symmetries of Block's 
  //	components.
  //
  if (b_sites == 1) {
    cout << "Undefined symmetry " << name_action (id) << endl;
    exit (0);
  }
  Action & lop = b_parent [b_lft] .symmetry (id);
  Action & rop = b_parent [b_rgt] .symmetry (id);
  tensor (s, lop, rop);
  return s;
}
//
//____________________________________________________________________________
void Block::tensor (Action & t, const Action & lop, const Action & rop)
{
  //
  //	Build Action t as tensor product of two Action's lop and rop
  //	belonging to left and right Block components
  //
  if (lop .isnull () || rop .isnull ()) {
    //
    //	null tensor product
    //
    t = Action (b_quantum, b_quantum);
    return;
  }
  if (lop .isscalar () && rop .isscalar ()) {
    //
    //	Two scalars product
    //
    complex<double> zl = lop .scalar ();
    complex<double> zr = rop .scalar ();
    t = Action (b_quantum, zl * zr);
    return;
  }
  if (lop .isscalar ()) {
    //
    //	lop is a scalar. Build identity action times the scalar
    //  and use it instead of lop
    //
    Action a (lop);
    a .scalaridentity (a .scalar ());
    tensor (t, a, rop);
    return;
  }
  if (rop .isscalar ()) {
    //
    //	rop is a scalar. Build identity action times the scalar    
    //	and use it instead of rop
    //
    Action a (rop);
    a .scalaridentity (a .scalar ());
    tensor (t, lop, a);
    return;
  }
  //  
  size_t lstates = b_parent [b_lft] .states ();
  //
  //	Initialize Action structure according to b_tensor space description
  //
  t = Action (b_tensor, b_tensor);
  size_t tsize = t .size ();
  //
  //	Determine block structure and memory size of tensor product
  //
  vector<Ablock> blk;
  Ablock * lb = lop .block ();
  Ablock * rb = rop .block ();
  size_t nl, nr, nt, il, jl, ir, jr, it, jt, lr, li, rr, ri;
  //
  //	Scan right and letf Action
  //
  for (nr = 0; nr < rop .blocks (); nr++) {
    ir = rb [nr] .ab_range;
    jr = rb [nr] .ab_domain;
    //
    //   right memory offsets
    //
    rr = rb [nr] .ab_roffset;
    ri = rb [nr] .ab_ioffset;
    //
    //	 right index offsets
    //
    ir = rop .range  () .offset (ir);
    jr = rop .domain () .offset (jr);
    //
    for (nl = 0; nl < lop .blocks (); nl++) {
      il = lb [nl] .ab_range;
      jl = lb [nl] .ab_domain;
      //
      //   left memory offsets
      //
      lr = lb [nl] .ab_roffset;	
      li = lb [nl] .ab_ioffset;
      //
      //   left index offsets
      //
      il = lop .range  () .offset (il);
      jl = lop .domain () .offset (jl);
      //
      //  tensor block subspaces
      //
      it = b_tensororder [il + ir * lstates];
      jt = b_tensororder [jl + jr * lstates];
      it = b_tensorsub   [it];
      jt = b_tensorsub   [jt];
      //
      //   check if tensor block was yet initialized
      //
      for (nt = 0; nt < blk .size (); nt++) 
	if ((blk [nt] .ab_domain == jt) && (blk [nt] .ab_range == it)) break;
      if (nt == blk. size ()) blk .push_back (Ablock ());
      //
      //   update tensor block
      //
      Ablock & b = blk [nt];
      b .ab_range   = it;
      b .ab_domain  = jt;
      size_t ms = t .height (it) * t .width (jt);
      if ((b .ab_roffset == 0) && ((lr && rr) || (li && ri))) {
	b .ab_roffset = tsize;
	tsize += ms;
      }
      if ((b .ab_ioffset == 0) && ((lr && ri) || (li && rr))) {
	b .ab_ioffset = tsize;
	tsize += ms;
      }
    }
  }	  
  //
  //	Allocate blocks inside tensor
  //
  nt = blk .size ();
  t .blocks (nt);
  t .block  (nt);
  Ablock * tb = t .block ();
  for (nt = 0; nt < blk .size (); nt++) tb [nt] = blk [nt];
  t .statistic (lop .statistic () * rop .statistic ());
  //
  //	Allocate tensor memory
  //
  t .size    (tsize);
  t .storage (tsize);
  //
  //	Fill memory
  //
  double * ml = lop .storage ();
  double * mr = rop .storage ();
  double * mt = t   .storage ();
  t .scalar (lop .scalar () * rop .scalar ());
  //
  //	Repeat previous loop filling memory
  //
  for (nr = 0; nr < rop .blocks (); nr++) {
    ir = rb [nr] .ab_range;
    jr = rb [nr] .ab_domain;
    long fr = rop .domain () .statistic (jr)  * rop .range () .statistic (ir);
    rr = rb [nr] .ab_roffset;
    ri = rb [nr] .ab_ioffset;
    //
    //	right subspaces offsets and sizes
    //
    size_t rh = rop .height (ir);
    size_t rw = rop .width  (jr);
    ir = rop .range  () .offset (ir);
    jr = rop .domain () .offset (jr);
    //
    for (nl = 0; nl < lop .blocks (); nl++) {
      il = lb [nl] .ab_range;
      jl = lb [nl] .ab_domain;
      long fl = lop .domain () .statistic (jl);
      double fermi = 1.0;
      if ((fl < 0) && (fr < 0)) fermi = -1.0;
      lr = lb [nl] .ab_roffset;
      li = lb [nl] .ab_ioffset;
      //
      //  left subspaces offsets and sizes
      //
      size_t lh = lop .height (il);
      size_t lw = lop .width  (jl);
      il = lop .range  () .offset (il);
      jl = lop .domain () .offset (jl);
      //
      //   tensor subspaces 
      //
      it = b_tensororder [il + ir * lstates];
      jt = b_tensororder [jl + jr * lstates];
      it = b_tensorsub   [it];
      jt = b_tensorsub   [jt];
      //
      //  find tensor block
      //
      for (nt = 0; nt < blk .size (); nt++) 
	if ((blk [nt] .ab_domain == jt) && (blk [nt] .ab_range == it)) break;
      //
      //  tensor block offsets and size
      //
      size_t h = b_tensor .states (it);
      it = b_tensor .offset (it);
      jt = b_tensor .offset (jt);
      size_t roff = blk [nt] .ab_roffset;
      size_t ioff = blk [nt] .ab_ioffset;
      //
      if (lr && rr) 
	// 
	//	real * real
	//
	for (size_t j2 = 0; j2 < rw; j2++) 
	  for (size_t j1 = 0; j1 < lw; j1++) {
	    //
	    //	tensor domain index inside block
	    //
	    size_t j = b_tensororder [j1 + jl + (j2 + jr) * lstates] - jt;
	    for (size_t i2 = 0; i2 < rh; i2++) 
	      for (size_t i1 = 0; i1 < lh; i1++) {
		//
		// tensor range index inside block
		//
		size_t i = b_tensororder [i1 + il + (i2 + ir) * lstates ] - it;
		mt [roff + i + j * h] += fermi *
		  ml [lr + i1 + j1 * lh] * mr [rr + i2 + j2 * rh];
	      }
	  }
      if (li && ri)
	//
	//	imag * imag	
	//
	for (size_t j2 = 0; j2 < rw; j2++) 
	  for (size_t j1 = 0; j1 < lw; j1++) {
	    size_t j = b_tensororder [j1 + jl + (j2 + jr) * lstates] - jt;
	    for (size_t i2 = 0; i2 < rh; i2++) 
	      for (size_t i1 = 0; i1 < lh; i1++) {
		size_t i = b_tensororder [i1 + il + (i2 + ir) * lstates ] - it;
		mt [roff + i + j * h] -= fermi *
		  ml [li + i1 + j1 * lh] * mr [ri + i2 + j2 * rh];
	      }
	  }
      if (li && rr) 
	//
	//	imag * real
	//
	for (size_t j2 = 0; j2 < rw; j2++) 
	  for (size_t j1 = 0; j1 < lw; j1++) {
	    size_t j = b_tensororder [j1 + jl + (j2 + jr) * lstates] - jt;
	    for (size_t i2 = 0; i2 < rh; i2++) 
	      for (size_t i1 = 0; i1 < lh; i1++) {
		size_t i = b_tensororder [i1 + il + (i2 + ir) * lstates ] - it;
		mt [ioff + i + j * h] += fermi *
		  ml [li + i1 + j1 * lh] * mr [rr + i2 + j2 * rh];
	      }
	  }
      if (lr && ri) 
	//
	//	real * imag
	//
	for (size_t j2 = 0; j2 < rw; j2++) 
	  for (size_t j1 = 0; j1 < lw; j1++) {
	    size_t j = b_tensororder [j1 + jl + (j2 + jr) * lstates] - jt;
	    for (size_t i2 = 0; i2 < rh; i2++) 
	      for (size_t i1 = 0; i1 < lh; i1++) {
		size_t i = b_tensororder [i1 + il + (i2 + ir) * lstates ] - it;
		mt [ioff + i + j * h] += fermi *
		  ml [lr + i1 + j1 *  lh] * mr [ri + i2 + j2 * rh];
	      }
	  }
    }
  }
  //
  if (action (idop_base,0) .storage ()) {
    //
    //	b_quantum is different from b_tensor.  A different basis 
    //  is used and action idop_base (from b_quantum to b_tensor)
    //	represent the new base
    //
    Action u = action (idop_base, 0);
    t *= u;
    u .dagger ();
    u *= t;
    t << u;
    action (idop_base, 0) .release ();
  }
  t .release ();
}
//
//============================================================================
