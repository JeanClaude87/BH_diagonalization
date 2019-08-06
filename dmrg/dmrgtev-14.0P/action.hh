//	action.hh			(C) 2009 Fabio Ortolani fbo 090716
//	==================================================================
#ifndef ACTION_HH
#define ACTION_HH
#include "storage.hh"
#include "space.hh"
//
//============================================================================
struct Abinary
{
  //
  //	Abinary structure contains informations suitable to perform a binary
  //	product on two Action's blocks using only non null memory area.
  //
  size_t  	ab_lftro;	// left operand memory real offset
  size_t  	ab_lftio;	// left operand memory imaginary offset
  size_t  	ab_rgtro;	// right operand memory real offset
  size_t  	ab_rgtio;	// right operand memory imaginary offset
  size_t  	ab_resro;	// result memory real offset
  size_t  	ab_resio;	// result memory imaginary offset
  //
  //	sizes for sum operation and product operations
  //
  long	  	ab_lftsz;	// left operand range size
  long	  	ab_rgtsz;	// right operand domain size
  long	  	ab_intsz;	// common left domain and right range size
  //
  //	The sign of sizes is negative if che corresponding states
  //	are fermionic.
  //
  //	Sparsed matrix indexes (if not null the operand is sparsed
  //	and the field give the offset of the corresponding indexes
  //	list, column-indexed for left operand, row-indexed for right
  //	operand).
  //	These fields can be ignored when sparse algorithms are not used.
  //
  size_t	ab_lftsp;	// left operand sparse indexes
  size_t	ab_rgtsp;	// right operand sparse indexes
  //
  Abinary ();							// [action.cc]
  ~Abinary ();							// [action.cc]
  //
};
//
//____________________________________________________________________________
struct Ablock
{
  //
  //	Ablock structure contains informations relative to a block
  //	submatrix.
  //	Null values denote missing elements (real part, imaginary part,
  //	sparsing)
  //
  size_t		ab_roffset;	// real part memory offset
  size_t		ab_ioffset;	// imaginary part memory offset 
  size_t		ab_domain;	// domain subspace index
  size_t		ab_range;	// range  subspace index
  size_t		ab_cindex;	// column-indexed sparse index list
  size_t	       	ab_rindex; 	// row-indexed sparse index list
  //
  Ablock ();							// [action.cc]
  ~Ablock ();							// [action.cc]
  //
  //	Comparison (based on domain and range indicators)
  //
  inline bool operator  < (const Ablock &) const;              	// [action.hh]
  inline bool operator <= (const Ablock &) const;      		// [action.hh]
  inline bool operator == (const Ablock &) const;      		// [action.hh]
  inline bool operator >= (const Ablock &) const;      		// [action.hh]
  inline bool operator  > (const Ablock &) const;      		// [action.hh]
  //
};
//
//____________________________________________________________________________
class Action
{
  //
  //	Action class contains a matrix divided in blocks.
  //
  //	The domain and range spaces of the matrix are partitioned into 
  //	subspaces and every block of the matrix expresses matrix's
  //	action from a subspace of the domain to a subspace of the range.
  //	The true subspaces are labelled by a non null index while 
  //	a null index labels the null subspace (see [space.hh]).
  //
  //	Ideally every subspace of the domain and range is spanned by a set
  //	of states (with the same quantum numbers) and the block structure 
  //	is suitable for operators which transform one subspace into only 
  //	one subspace and most blocks are null (this is the case when 
  //	operators carry well defined quantum properties).
  //
  //	The non null blocks informations are stored into an array of 
  //	Ablock's (see above).
  //
  //	The single blocks can be considered as sparse matrices if the 
  //	content of a_sparsed is not null. We do not implement a true
  //	sparsing of memory area, we build only the indexing part (both
  //	column-indexed and row-indexed part) and leave memory area 
  //	unchanged. This scheme is derived from Numeical Recipes in C++. 
  //	
private:
  //
  //	Domain (row specifications) and range (column specifications)
  //
  Space		a_range;       	// range  subspaces description	    [space.hh]
  Space 	a_domain;	// domain subspaces description	    [space.hh]
  //
  //	Memory and non null blocks sizes. These fields denote only 
  //	needed sizes of storage allocation. Actual storage can be
  //	different.
  //
  size_t	a_size;		// storage allocation (in double units)
  size_t	a_blocks;	// n. of non null blocks
  //
  //	Fermi-Bose Action character:	-1	Fermionic operator
  //					+1	Bosonic   operator
  //					 0	Undefined (treated as bosonic)
  //
  //	This field is used for reordering of formal products of operators. 
  //	It is well defined only when all non null blocks have the same 
  //	statistical property (based on statistical properties of range and 
  //	domain subspaces of the block).
  //
  long		a_statistic;	// Fermi-Bose statistic of operator
  //	
  //	Storage area for matrix entries.
  //
  //	Memory storage is non null if the Action is defined (as an 
  //	effective operator from a domain to a range).
  //
  //	First memory entry is a common multiplicatve complex factor 
  //	(a scalar). This common factor multiply all non null entries. 
  //
  //   	An Action is composed of two multiplicative parts: a scalar and 
  //	a matrix. 
  //	
  //	       	Action = scalar * matrix
  //
  //   	The matrix part structure is described by the Block storage.
  //
  //	If the common factor is the only memory area (a_size == csize 
  //	[action.cc]) the Action is considered a scalar (factor * Identity) 
  //	and the matrix part can be ignored.
  //
  //	An Action is said in normalized form if the scalar factor is 
  //	unitary (or zero as special case) or the matrix part is unitary 
  //	(i.e. Identity: the Action is a scalar). 
  //	We need normalized Action's to perform some binary operations
  //	(mainly for the sum of Action's: the scalar part or the matrices 
  //	must be the same).
  //
  //	An Action is null if the scalar factor is null. It can be 
  //	considered a special case of normalized Action). 
  //	      
  Storage	a_storage;	// memory storage for matrix entries
  //
  //	Blocks storage area is an array of Ablock's (at least
  //	one for each non null submatrix of the full matrix). 
  //	It can be void for a null matrix or a scalar (in such
  //	case the non null blocks are supposed to be Identity
  //	subblocks).
  //
  Storage	a_block;	// matrix blocks infos (Ablock array)
  //
  //	Sparse storage area is an array of size_t's with (if not void)
  //	indexes of non null elements of the matrix. If the matrix is 
  //	considered sparse in every Ablock the ab_cindex and ab_rindex
  //	are not null. The memory storage of matrix entries is not
  //	modified, so the sparse informations can be safely ignored.
  //	
  Storage	a_sparsed;	// sparse matrix infos (size_t array)
  //
public:
  //
  Action  ();							// [action.cc]
  Action  (const Action &);					// [action.cc]
  //
  //	Null Action ctor:	a null Action between two spaces
  //
  Action  (const Space &, const Space &);			// [action.cc]
  Action  (const Qspace &, const Qspace &, const Qnumber &, 
	   bool = false, const double = 1.0);			// [action.cc]
  //
  //	Scalar Action ctor:	scalar * Identity
  //
  Action  (const Space &, complex<double>);			// [action.cc]
  //
  //	Special Action ctor:	Action structure is created for 
  //	the product of two given Actions.
  //	Only Ablock structure is defined. Memory size for the product is 
  //	computed and registred in field a_size, not allocated.
  //	Formally it is a null matrix with a non null block structure. 
  //	(scalar factor is null).
  //
  Action  (const Action &, const Action &);		 	// [action.cc]
  //
  ~Action ();							// [action.cc]
  //
  //	getting fields content and properties.
  //
  inline const Space    & domain    () const;	       		// [action.hh]
  inline const Subspace & domain    (size_t) const;    		// [action.hh]
  inline const Space    & range     () const;   	       	// [action.hh]
  inline const Subspace & range     (size_t) const;    		// [action.hh]
  inline size_t           height    (size_t = 0) const;       	// [action.hh]
  inline size_t   	  width     (size_t = 0) const;       	// [action.hh]
  inline size_t           size      () const;   // required memory [action.hh]
  inline size_t	  	  blocks    () const;   // required blocks [action.hh]
  inline long	          statistic () const;  // global statistic [action.hh]
  inline double *         storage   () const;  // allocated memory [action.hh]
  inline Ablock * 	  block     () const;  // allocated Ablock [action.hh]
  inline size_t *	  sparsed   () const;			// [action.hh]
  //
  complex<double>         scalar    () const;        // get factor [action.cc]
  bool			  iscomplex () const;			// [action.cc]
  bool			  isnormal  () const;			// [action.cc]
  bool 			  isnull    () const;			// [action.cc]
  bool 			  isscalar  () const;	      		// [action.cc]
  double		  entropy   () const;		       	// [action.cc]
  size_t		  dimension () const;			// [action.cc]
  //
  //	setting methods
  //
  inline void	  size	    (size_t);		  		// [action.hh]
  inline void	  blocks    (size_t);				// [action.hh]
  inline void	  statistic (long);	       			// [action.hh]
  inline void 	  storage   (size_t); 	        // allocate memory [action.hh]
  inline void	  storage   (const Storage &);     // grab storage [action.hh]
  inline void	  block     (size_t);	        // allocate Ablock [action.hh]
  inline void	  block     (const Storage &);      // grab Ablock [action.hh]
  void            scalar    (complex<double>);       // set factor [action.cc]
  void 	          sparsed   (size_t); 	        // allocate sparse [action.cc]
  //
  //	modify methods
  //
  void add	      (complex<double>, const Action &);       	// [action.cc]
  void clean          ();                                       // [action.cc]
  void compress       (complex<double> *);             		// [action.cc]
  void dagger	      ();	       				// [action.cc]
  void expand	      (complex<double> *) const;       		// [action.cc]
  void merge	      (const Action &, size_t = 0);		// [action.cc]
  void normalize      ();					// [action.cc]
  void random	      ();					// [action.cc]
  void release        ();		        // release storage [action.cc]
  void scalaridentity (complex<double> = 1.0);       		// [action.cc]
  void sparse	      ();					// [action.cc]
  void transpose      (); 				        // [action.cc]
  //
  //	operations (modify memory and structure)
  //
  Action & operator  = (const Action &); 	// copy assignment [action.cc]
  Action & operator << (const Action &); 	// grab assignment [action.cc]
  Action & operator += (const Action &); 	       // addition [action.cc]
  Action & operator *= (complex<double>);       // scalar multiply [action.cc]
  Action & operator *= (const Action &);               // multiply [action.cc]
  //
  //	utilities
  //
  void show  (const string & = "") const;		       	// [action.cc]
  void showblocks (const string & = "") const;                  // [action.cc]
  //
};
//
//____________________________________________________________________________
class Aarray
{
  //
  //	Aarray is a container for a dynamic array of Action's accessible 
  //	like an ordinary array (with dynamic grow of array size).
  //	
  //	A single element (an Action) can be a large object in memory
  //	and dynamic resizing of array can involve allocation and moving
  //	of its elements. To avoid moving of large Action's the container
  //	is implemented as a dynamic array of pointers to Action's.
  //	Resizing and moving deals only with pointers, and 
  //	pointed Action's are created the first time the element
  //	is accessed and destroyed when the pointer is destroyed.
  //
  //	We don't want copies of the same Action living together 
  //	in memory. 	
  //	As a side effect we have to care when removing, moving,  and 
  //	changing array elements because we have to delete expicitly the 
  //	Action pointed by the pointer. Every Action pointed must be 
  //	deleted only once so we have to be sure that only one active 
  //	pointer to an Action is living in all actual Aarray's.
  //
  //	Note that the pointers are private and not accessible 
  //	from outside class. The array-like access returns the 
  //	content and not the pointer.
  //
private:
  //
  vector<Action *>	aa_action;
  //
public:
  //
  Aarray  ();							// [action.cc]
  Aarray  (const Aarray &);					// [action.cc]
  ~Aarray ();							// [action.cc]
  //
  inline void   clear   ();					// [action.hh]
  inline void   release ();					// [action.hh]	
  inline size_t size    () const;      				// [action.hh]
  //
  const Action & operator [] (size_t) const;	       		// [action.cc]
  Action & operator [] (size_t);	   	   		// [action.cc]
  Aarray & operator =  (const Aarray &);       			// [action.cc]
  //
};
//
//============================================================================
//
//	Action related procedures
//
//____________________________________________________________________________
void binary (vector<Abinary> &, const Action &,
	     const Action &, const Action &);			// [action.cc]
void binary (vector<Abinary> &,
	     const Action &, const Action &);			// [action.cc]
//
size_t minasize ();						// [action.cc]
//
void multiply (Action &, const Action &, const Action &,
	       bool = false); 				   	// [action.cc] 
void multiply (double *, double *, double *, 
	       vector<Abinary> &, bool = false);	   	// [action.cc]
void multiply (double *, double *, size_t *, double *, size_t *,
	       vector<Abinary> &, bool = false);	   	// [action.cc]
	       
//
void multiply (double *, complex<double>, 
	       double *, vector<Abinary> &);			// [action.cc]
//
complex<double> multiply (const Action &, const Action &);	// [action.cc]
complex<double>	multiply (double *, double *, 
			  vector<Abinary> &);			// [action.cc]
//
//============================================================================
//
//	Inline methods
//
//============================================================================
inline void Aarray::clear ()
{
  //
  //	Clear content
  //
  for (size_t i = 0; i < size (); i++)
    if (aa_action [i]) delete aa_action [i];
  //
  //	... and pointers
  //
  aa_action .clear ();
}
//
//____________________________________________________________________________
inline void Aarray::release ()
{
  //
  //	Releases all storage areas of pointed Action's 
  //
  for (size_t i = 0; i < size (); i++) 
    if (aa_action [i]) 	(*(aa_action [i])) .release ();
}
//
//____________________________________________________________________________
inline size_t Aarray::size () const
{
  //
  //	Returns size of Action's array
  //
  return aa_action .size ();
} 
//
//============================================================================
inline bool Ablock::operator  < (const Ablock & other) const
{
  //
  //	Least compare
  //
  return ((ab_domain < other .ab_domain) ||
	  ((ab_domain == other .ab_domain) && (ab_range < other .ab_range)));
}
//
//____________________________________________________________________________
inline bool Ablock::operator <= (const Ablock & other) const
{
  //
  //	Least/equal compare
  //
  return ((ab_domain < other .ab_domain) ||
	  ((ab_domain == other .ab_domain) && (ab_range <= other .ab_range)));
}
//
//____________________________________________________________________________
inline bool Ablock::operator == (const Ablock & other) const
{
  //
  //	Equality compare
  //
  return ((ab_domain == other .ab_domain) && (ab_range == other .ab_range)); 
}
//
//____________________________________________________________________________
inline bool Ablock::operator >= (const Ablock & other) const
{
  //
  //	greter/equal compare
  //
  return ((ab_domain > other .ab_domain) ||
	  ((ab_domain == other .ab_domain) && (ab_range >= other .ab_range)));
}
//
//____________________________________________________________________________
inline bool Ablock::operator > (const Ablock & other) const
{
  //
  //	Greater compare
  //
  return ((ab_domain > other .ab_domain) ||
	  ((ab_domain == other .ab_domain) && (ab_range > other .ab_range)));
}
//
//============================================================================
inline Ablock * Action::block () const
{
  //
  //	Allocated Ablock's
  //
  return (Ablock *) a_block .storage ();
}
//
//____________________________________________________________________________
inline void Action::block (size_t size)
{
  //
  //	Allocate size Ablock's
  //
  a_block .storage (size * sizeof (Ablock));
  a_blocks = size;
}
//
//____________________________________________________________________________
inline void Action::block (const Storage & blocks)
{
  //
  //	Grab Ablock's from blocks
  //
  a_block << blocks;
  a_blocks = a_block .size () / sizeof (Ablock);
}
//
//____________________________________________________________________________
inline size_t Action::blocks () const
{
  //
  //	n. of required blocks
  //
  return a_blocks;
}
//
//____________________________________________________________________________
inline void Action::blocks (size_t size)
{
  //
  //	Set n. of required blocks
  //
  a_blocks = size;
}
//
//____________________________________________________________________________
inline const Space & Action::domain () const
{
  //
  //	Returns domain field
  //	
  return a_domain;
} 
//
//____________________________________________________________________________
inline const Subspace & Action::domain (size_t subspace) const
{
  //
  //	Returns domain subspace
  //
  return a_domain [subspace];
}
//
//____________________________________________________________________________
inline size_t Action::height (size_t subspace) const
{	
  //
  //	Return column height (n. of rows) of subspace 
  //	(full matrix if subspace == 0)
  //
  return a_range .states (subspace);
}
//
//____________________________________________________________________________
inline const Space & Action::range () const
{
  //
  //	Returns range field
  //	
  return a_range;
} 
//
//____________________________________________________________________________
inline const Subspace & Action::range (size_t subspace) const
{
  //
  //	Returns range subspace
  //
  return a_range [subspace];
}
//
//____________________________________________________________________________
inline size_t Action::size () const
{
  //	
  //	Requested memory allocation (in double units)
  //
  return a_size;
}
//
//____________________________________________________________________________
inline void Action::size (size_t size)
{
  //	
  //	set size of required allocation
  //
  a_size = size;
}
//
//____________________________________________________________________________
inline size_t * Action::sparsed () const
{
  //
  //	Returns sparse indexes array
  //
  return (size_t *) (a_sparsed .storage ());
}
//
//____________________________________________________________________________
inline long Action::statistic () const
{
  //
  //	Fermi-Bose statistic
  //	
  return a_statistic;
}
//
//____________________________________________________________________________
inline void Action::statistic (long stat) 
{
  //
  //	set Fermi-Bose statistic
  //	
  a_statistic = stat;
}
//
//____________________________________________________________________________
inline double * Action::storage () const
{
  //
  //	Allocated memory pointer
  //	
  return (double *) (a_storage .storage ());
}
//
//____________________________________________________________________________
inline void Action::storage (size_t size)
{
  //
  //	Allocate memory for size double's
  //
  a_storage .storage (size * sizeof(double));
}
//
//____________________________________________________________________________
inline void Action::storage (const Storage & store)
{
  //
  //	Grab memory from store
  //
  a_storage << store;
}
//
//____________________________________________________________________________
inline size_t Action::width (size_t subspace) const
{	
  //
  //	Return rows width (n. of columns) of subspace
  //	(full matrix width if subspace == 0)
  //
  return a_domain .states (subspace);
}
//
//============================================================================
#endif // ACTION_HH
