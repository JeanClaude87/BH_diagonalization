//	superaction.hh			(C) 2009 Fabio Ortolani	fbo 090707
//	==================================================================
#ifndef SUPERACTION_HH
#define SUPERACTION_HH
#include "block.hh"
//
//============================================================================
struct Biaction 
{
  //
  //	Tensor product of two Action's acting on superblock
  //
  size_t		ba_size;	// size of intermediate memory
  Action       		ba_lftop;	// left Action
  Action		ba_rgtop;	// right Action
  vector<Abinary> 	ba_lftab;	// right binary application
  vector<Abinary> 	ba_rgtab;	// left binary application
  //
  Biaction  ();						   // [superaction.cc]
  Biaction (const Biaction &);				   // [superaction.cc]
  ~Biaction ();						   // [superaction.cc]
  //
};
//
//____________________________________________________________________________
struct Bipoli
{
  //
  //	A formal product of two polinomials of site actions
  //	acting on different lattice blocks
  //
  complex<double> 	bi_coeff;
  Apoli			bi_lft;
  Apoli			bi_rgt;
  //
  Bipoli  ();						   // [superaction.cc]
  Bipoli  (const Bipoli &);				   // [superaction.cc]
  ~Bipoli ();						   // [superaction.cc]
  //
};
//
//____________________________________________________________________________
class Superaction
{
  //
  //	Superaction represents an operator acting on a superblock
  //	tensor product of two Block's (left and right components)
  //
  //	The operator is expressed as a sum of tensor products of
  //	two Action's acting on left and right Block's spaces
  //
  //	We store both formal expression and true Action's tensor
  //	components
  //
private:
  vector<Bipoli>	sa_bipoli;	
  vector<Biaction>	sa_biaction;
  Storage		sa_inner;
  //
public:
  Superaction  ();					   // [superaction.cc]
  Superaction  (const Superaction &);			   // [superaction.cc]
  //
  //	Build a superaction from its formal polinomial
  //	expression
  //
  Superaction  (const Apoli &, Block &, Block &, 
		bool = false);      	  		   // [superaction.cc]
  ~Superaction ();					   // [superaction.cc]
  //
  //	Array-like access to true Action's tensor part
  //
  inline Biaction & operator [] (size_t);      		   // [superaction.hh]
  //
  inline complex<double> coefficient     (size_t) const;   // [superaction.hh]
  inline const Apoli &   leftpolinomial  (size_t) const;   // [superaction.hh]
  inline const Apoli &   rightpolinomial (size_t) const;   // [superaction.hh]
  inline size_t	    	 size            ()       const;   // [superaction.hh]
  inline double *	 storage	 ()	  const;   // [superaction.hh]
  //
  void storage   (size_t);			       	   // [superaction.cc]
  bool iscomplex () const;				   // [superaction.cc]
  void show      (const string & s = "") const;		   // [superaction.cc]
  //
  Superaction & operator = (const Superaction &);	   // [superaction.cc]
  //
};
//
//____________________________________________________________________________
class Lanczos
{
  //
  //	Lanczos class is a container for Lanczos vectors, eigenvalues
  //	and Lanczos algorithm
  //
private:
  //	
  Storage *    	l_vector;	// Lanczos vectors or eigenvectors
  double  *  	l_value;	// Found eigenvalues
  double  *    	l_tolerance;	// Eigenvalue's tolerances
  double  *	l_error;	// Eigenvalue's estimated errors
  //
  //	Parameters
  //
  size_t       	l_size;		  // secular dimension
  size_t        l_random;         // random vector dimension
  size_t       	l_found;	  // n. of found or wanted lowest eigenvalues
  //
  //	Convergency strategy:
  //
  //		0:	only lowest eigenvalues
  //		1:	convercengy on lowest, keep also highest for restart
  //		2:	lowest and highest eigenvalues (very demanding)
  //
  size_t	l_strategy;	  
  size_t       	l_steps;	  // n. of steps
  size_t       	l_repeats;	  // n. of repeatitions 
  //
  double       	l_zerotolerance;  // precision thereshold
  double      	l_zeronorm;	  // zero norm thereshold
  //
public:
  //
  Lanczos  ();						   // [superaction.cc]
  ~Lanczos ();						   // [superaction.cc]
  //
  //	Set methods
  //
  void 	      parameters (size_t, size_t, 
			  size_t, size_t, size_t, 
			  size_t, size_t, double, double); // [superaction.cc]
  void	      random	 (size_t);			   // [superaction.cc]
  inline void release	 (size_t);			   // [superaction.hh]
  inline void remove	 (size_t);			   // [superaction.hh]
  void 	      storage	 (size_t, const Storage &); 	   // [superaction.cc]
  //
  inline size_t   steps	    () const;			   // [superaction.hh]
  inline double * value     () const; 			   // [superaction.hh]
  inline double   value	    (size_t) const;   		   // [superaction.hh]
  inline double   tolerance (size_t) const;   		   // [superaction.hh]
  inline double   error     (size_t) const;		   // [superaction.hh]
  //
  //	Get methods
  //
  double *  	   operator [] (size_t);	       	   // [superaction.cc]
  inline Storage & storage     (size_t);       		   // [superaction.hh]
  //
  //	Diagonalization algorithm
  //
  void	  krylov (double *, double *, double *, 
		  size_t, size_t);			   // [superaction.cc]
  size_t  thick  (Superaction &, Action &);		   // [superaction.cc]
  //
};
//
//============================================================================
size_t biaction (Action &, Superaction &, const Action &,
		 bool = false); 			   // [superaction.cc]
void   biapply  (double *, Superaction &, double *,
		 bool = false);			           // [superaction.cc]
void   biapply  (Action &, Superaction &, const Action &,
		 bool = false); 			   // [superaction.cc]
//
//============================================================================
inline double Lanczos::error (size_t index) const
{
  //	
  //	Returns error
  //   
  return l_error [index];
}
//
//____________________________________________________________________________
inline void Lanczos::release (size_t index)
{
  //
  //	Release index-th vector memory area
  //
  l_vector [index] .release ();
}
//
//____________________________________________________________________________
inline void Lanczos::remove (size_t index)
{
  //
  //	Remove allocation for index-th vector
  //
  l_vector [index] .storage (0);
}
//
//____________________________________________________________________________
inline size_t Lanczos::steps () const
{
  //
  //	Returns steps
  //
  return l_steps;
}
//
//____________________________________________________________________________
inline Storage & Lanczos::storage (size_t index)
{
  //
  //	get index-th Storage
  //
  return l_vector [index];
}
//
//____________________________________________________________________________
inline double Lanczos::tolerance (size_t index) const
{
  //	
  //	Returns tolerance
  //	
  return l_tolerance [index];
}
//
//____________________________________________________________________________
inline double * Lanczos::value () const
{
  //	
  //	Returns eigenvalue pointer
  //	
  return l_value;
}
//
//____________________________________________________________________________
inline double Lanczos::value (size_t index) const
{
  //	
  //	Returns eigenvalue
  //	
  return l_value [index];
}
//
//============================================================================
inline complex<double> Superaction::coefficient (size_t index) const
{
  //
  //	coefficient of index-th tensor
  //
  return sa_bipoli [index] .bi_coeff;
}
//
//____________________________________________________________________________
inline const Apoli & Superaction::leftpolinomial (size_t index) const
{
  //
  //	left polinomial of index-th tensor product
  //
  return sa_bipoli [index] .bi_lft;
}
//
//____________________________________________________________________________
inline Biaction & Superaction::operator [] (size_t index) 
{
  //
  //	Array-like access	
  //
  return sa_biaction [index];
}
//
//____________________________________________________________________________
inline const Apoli & Superaction::rightpolinomial (size_t index) const
{
  //
  //	right polinomial of index-th tensor product
  //
  return sa_bipoli [index] .bi_rgt;
}
//
//____________________________________________________________________________
inline size_t Superaction::size () const
{
  //
  //	n. of tensor products
  //	
  return sa_biaction .size ();
}
//
//____________________________________________________________________________
inline double * Superaction::storage () const
{
  //
  //	Inner memory
  //	
  return ((double *) sa_inner .storage ());
}
//
//============================================================================
#endif // SUPERACTION_HH
