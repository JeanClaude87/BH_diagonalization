//	block.hh			(C) 2009 Fabio Ortolani	fbo 090716
//	==================================================================
#ifndef BLOCK_HH
#define BLOCK_HH
#include "action.hh"
#include "algebra.hh"
using namespace std;
//
//============================================================================
class 	Barray;	// forward reference
extern  Barray 	blockempty;
//
//____________________________________________________________________________
class Block
{
  //
  //	Block class describes the structure of a quantum lattice.
  //	Its members give the description of the corresponding Hilbert space
  //	(and subspaces) and of operators acting on this space (hamiltonian,
  //	global and local site operators).
  //
private:
  //
  size_t		b_sites;       	// n. of lattice stites
  size_t		b_states;	// n. of Hilbert space states
  //
  //	Remember original blocks components and construction
  //
  Barray &		b_parent;	// reference to list of parent blocks
  size_t		b_lft;		// n. of left block sites
  size_t		b_rgt;		// n. of right block sites
  Qspace		b_tensor;	// original tensor structure
  vector<size_t>	b_tensororder;	// generation order of tensor states
  vector<size_t>	b_tensorsub;	// tensor subspace of states
  //
  vector<size_t>	b_siteorder;	// index of lattice sites
  vector<size_t>	b_optimized;	// basis states optimization
  //
  Qspace		b_quantum;	// quantum space structure
  vector<size_t>	b_quantumsub;	// quantum subspace of states 
  //
  vector<Aarray>	b_action;	// operators associated to this 
  //					   block spaces
  //
public:
  //
  Block (Barray & = blockempty);				 // [block.cc]
  Block (const Block &);     					 // [block.cc]
  Block (const Block &, Barray &);		       		 // [block.cc]
  ~Block ();							 // [block.cc]
  //
  //	Special single site Block constructor
  //
  Block (const long &, Barray & = blockempty);         		 // [block.cc]
  //
  //	Every lattice block with more than one site is built composing two
  //	sublattices and the Hilbert space is the tensor product of the two
  //	Hilbert spaces (taking into account of the antisymmetry of fermionic 
  //	states) 
  //
  Block (const Block &, const Block &);				 // [block.cc]
  //
  inline size_t			sites  		() const;	 // [block.hh]
  inline size_t			siteslft       	() const;        // [block.hh]
  inline size_t			sitesrgt       	() const;        // [block.hh]
  inline size_t			states 		() const;	 // [block.hh]
  //
  inline const vector<size_t> & siteorder  	() const;      	 // [block.hh]
  //
  inline const Qspace &		quantum 	() const;      	 // [block.hh]
  inline const vector<size_t> & quantumsub	() const;	 // [block.hh]
  inline size_t			quantumsub	(size_t) const;	 // [block.hh]
  inline long int              	quantumstat 	(size_t) const;  // [block.hh]
  //
  inline const Qnumber & 	number 		(size_t) const;	 // [block.hh]
  inline size_t 		subspaces 	() const;    	 // [block.hh]
  //
  inline const vector<size_t> & tensororder 	() const;	 // [block.hh]
  inline size_t			tensororder	(size_t) const;  // [block.hh]
  inline const vector<size_t> & tensorsub  	() const;	 // [block.hh]
  inline size_t			tensorsub	(size_t) const;	 // [block.hh]
  inline long int		tensorstat	(size_t) const;	 // [block.hh]
  //
  Block &	operator =  (const Block &);			 // [block.cc]
  Block &	operator << (const Block &);	       		 // [block.cc]
  //
  //	Block's Action's
  //
  //	Reference (without construction)
  //
  Action       & action (size_t, size_t = 0);			 // [block.cc]
  const Action & action (size_t, size_t = 0) const;    		 // [block.cc]
  //
  //	Build Action representing a polinomial of operators 
  //
  Action & action   (const Apoli &, Action &);			 // [block.cc]
  //
  // 	Optimization of basis
  //
  void	   optimize (const Apoli &);				 // [block.cc]
  void	   optimize (Action &);					 // [block.cc]
  //
  //	Builds Action representing Block's basis states
  //	(as an effective Action from b_quantum to b_tensor)
  //
  Action & base ();						 // [block.cc]
  Action & base (Action &);				         // [block.cc]
  //
  //	Builds Action representing Block's reflected b_tensor
  //	states projected in b_quantum space
  //
  Action & basereflected ();					 // [block.cc]
  Action & basereflected (Action &);		       		 // [block.cc]
  //
  Action & symmetry (size_t);					 // [block.cc]
  Action & symmetry (size_t, Action &);				 // [block.cc]
  //
  void actionclear   (size_t, size_t);				 // [block.cc]
  void add_states    (const Qnumber &, size_t, long);	    	 // [block.cc]
  //
  void reflect       ();	       				 // [block.cc]
  void reflectlft    ();	       				 // [block.cc]
  void reflectreset  ();       					 // [block.cc]
  void reflectrgt    ();       					 // [block.cc]
  void release       ();                                         // [block.cc]
  //
  void select_states (const Action &, size_t, size_t, double);		 // [block.cc]
  void show          (const string & = "");			 // [block.cc]
  void tensor	     (Action &, const Action &, const Action &); // [block.cc]
  //
  //	Make a partial product of two actions contracting left
  //	b_tensor components and returns an Aarray indexed by 
  //	right components states of b_tensor (the domain of
  //	left operand or the range of the right operand must be
  //	equal to b_tensor Space)
  //
  void contraction (const Action &, size_t, size_t, const Action &, 
		    Aarray &, bool = false) const;     		 // [block.cc]
  //
private:
  //
  //	These methods are private to Block.
  //
  Action & action  (const Amono &, Action &);          		 // [block.cc] 
  Action & action  (const Amono &);				 // [block.cc] 
  //
};
//
//____________________________________________________________________________
class Barray
{
  //
  //	Barray is a container for a dynamic array of Block's accessible 
  //	like an array
  //
private:
  //
  vector<Block>		ba_block;		
  //
public:
  //
  Barray  ();							 // [block.cc]
  ~Barray ();							 // [block.cc]
  //
  inline void   clear ();					 // [block.hh]
  inline size_t	size  () const;					 // [block.hh]
  //
  Block  & operator [] (size_t);   // array-like access		    [block.cc]
  //
};
//
//____________________________________________________________________________
struct	Brule
{
  //
  //	Describes a rule to select DMRG Block's density states
  //
  double    br_time;    // initial time for rule application
  size_t	br_zip;		// initial zip for rule application
  size_t	br_sites;	// initial chain length for rule application
  size_t	br_min;		// minimum n. of states
  size_t	br_max;		// maximum n. of states
  double	br_cut;		// maximum truncation error
  //
};
//
//============================================================================
inline void Barray::clear ()
{
  //
  //	Clear content
  //
  ba_block .clear ();
}
//
//____________________________________________________________________________
inline size_t Barray::size () const
{
  //
  //	size of ba_block
  //
  return ba_block .size ();
}
//
//============================================================================
inline const Qnumber & Block::number (size_t sub) const
{
  //	
  //	Quantum number of sub subspace
  //
  return b_quantum .number (sub);
}
//
//____________________________________________________________________________
inline const Qspace & Block::quantum () const
{
  //
  //	Quantum space
  //
  return b_quantum;
}
//
//____________________________________________________________________________
inline const vector<size_t> & Block::quantumsub () const
{
  // 
  //	The array of subspace for every tensor state
  //
  return b_quantumsub;
}
//
//____________________________________________________________________________
inline size_t Block::quantumsub (size_t state) const
{
  // 
  //	The subspace of tensor state
  //
  return b_quantumsub [state];
}
//
//____________________________________________________________________________
inline long int Block::quantumstat (size_t subspace) const
{
  // 
  //	The statistic of tensor subsspace
  //
  return b_quantum .statistic (subspace);
}
//
//____________________________________________________________________________
inline const vector<size_t> & Block::siteorder () const
{
  // 
  //	The array of sites ordering
  //
  return b_siteorder;
}
//
//____________________________________________________________________________
inline size_t Block::sites () const
{
  //
  //	Returns n. of sites 
  //
  return b_sites;
}
//
//____________________________________________________________________________
inline size_t Block::siteslft () const
{
  //
  //	Returns n. of sites of left component Block taking into account
  // 	possible reflections
  //
  if (b_sites == 0) return 0;
  if (b_siteorder [0] == 0)   return b_lft;	// no reflection
  if (b_siteorder [0] == b_sites - 1) return b_rgt;  // global reflection
  return b_lft;
}
//
//____________________________________________________________________________
inline size_t Block::sitesrgt () const
{
  //
  //	Returns n. of sites of right component Block taking into account
  //	possible reflections
  //
  if (b_sites == 0) return 0;
  if (b_siteorder [0] == 0) return b_rgt; 		// no reflection
  if (b_siteorder [0] == b_sites - 1)   return b_lft;	// global reflection
  return b_rgt;
}
//
//____________________________________________________________________________
inline size_t Block::states () const
{
  //
  //	N. of Block's states
  //
  return b_states;
}
//
//____________________________________________________________________________
inline size_t Block::subspaces () const
{
  //
  //	Returns the number of quantum subspaces
  //
  return b_quantum .subspaces ();
}
//
//____________________________________________________________________________
inline const vector<size_t> & Block::tensororder () const
{
  // 
  //	The array of sites ordering
  //
  return b_tensororder;
}
//
//____________________________________________________________________________
inline size_t Block::tensororder (size_t index) const
{
  //
  //	the insertion order of index state
  //
  return b_tensororder [index];
}
//
//____________________________________________________________________________
inline long int Block::tensorstat (size_t subspace) const
{
  // 
  //	The statistic of tensor subsspace
  //
  return b_tensor .statistic (subspace);
}
//
//____________________________________________________________________________
inline const vector<size_t> & Block::tensorsub () const
{
  // 
  //	The array of subspace for every tensor state
  //
  return b_tensorsub;
}
//
//____________________________________________________________________________
inline size_t Block::tensorsub (size_t state) const
{
  // 
  //	The subspace of tensor state
  //
  return b_tensorsub [state];
}
//
//============================================================================
#endif // BLOCK_HH
