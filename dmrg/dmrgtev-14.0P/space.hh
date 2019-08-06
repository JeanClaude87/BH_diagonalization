//	space.hh			(C) 2009 Fabio Ortolani	fbo 090716
//	==================================================================
#ifndef SPACE_HH
#define SPACE_HH
#include <complex>
#include <vector>
//
using namespace std;
//
//============================================================================
class	Subspace
{
  //
  //	Contains algebraic specifications of a subspace
  //
private:
  //
  size_t ss_offset;	// offset of first state of subspace
  size_t ss_states;	// number of states
  long   ss_statistic;	// fermi-bose statistic of all states in subspace
  //			   0 none, -1 fermi, 1 bose
public:
  //
  inline Subspace ();						  // [space.hh]
  //
private:
  //
  friend class Space;	// Only Space class can access Space_entry
  //
};
//
//____________________________________________________________________________
class Space
{
  //
  //	Container for a collection of subspaces defining a linear space
  //	It is implemented as an array of SubSpace classes. 
  //	Every space contains always the null subspace labelled by a null 
  //	index. It is used to store global informations.
  //
protected:
  //
  vector<Subspace>	s_subspace;	// Only true member
  //
public:
  //
  Space ();							 // [space.cc]
  Space (const Space &);					 // [space.cc]
  ~Space ();							 // [space.cc]
  //
  inline const Subspace & operator [] (size_t) const;		 // [space.hh]
  inline const Subspace & subspace    (size_t) const;		 // [space.hh]
  //
  inline size_t offset	  (size_t)     const;     		 // [space.hh]
  inline size_t states    (size_t = 0) const; 		  	 // [space.hh]
  inline long	statistic (size_t)     const;			 // [space.hh]
  inline size_t	subspaces () 	       const;			 // [space.hh]
  //
  bool operator == (const Space &) const;			 // [space.cc]
  //
protected:
  //
  inline void offsets	();					 // [space.hh]
  inline void states	(size_t, size_t);      			 // [space.hh]
  inline void statistic	(size_t, long);				 // [space.hh]
  //
};
//
//____________________________________________________________________________
class	Qnumber
{
  //
  //	Qnumber class collects quantum numbers for a state or subspace
  //	We admit only two types of quantum properties: additive real numbers 
  //	and multiplicative complex numbers with respect to composition
  //	(tensor product) of states and subspaces.
  //	The sizes of collections of numbers is constant inside the program
  //	(it depends on input and verified in input parsing).
  //
private:
  //
  size_t		   qn_id;      	// identifier (if not null)
  vector<double>	   qn_value;    // additive quantum numbers
  vector<complex<double> > qn_parity; 	// multiplicative properties
  double		   qn_lftlnk;
  double		   qn_rgtlnk;
  //
public:
  //
  Qnumber  ();							 // [space.cc]
  Qnumber  (const Qnumber &);					 // [space.cc]
  Qnumber  (size_t, const vector<double> &,			 //
	    const vector<complex<double> > &,			 // 
	    double, double);					 // [space.cc]
  ~Qnumber ();							 // [space.cc]
  //
  inline size_t id () const;					 // [space.hh]
  inline void	id (size_t);					 // [space.hh]
  string        str () const;  					 // [space.cc]
  void      	reflect ();			       		 // [space.cc]
  //
  bool operator == (const Qnumber &) const;			 // [space.cc]
  //
private:
  //
  friend class Qspace;
  friend Qnumber compose    (const Qnumber &, const Qnumber &,
			     const double);	                 // [space.cc]
  friend Qnumber operator + (const Qnumber &, const Qnumber &);	 // [space.cc]
  friend int	 operator - (const Qnumber &, const Qnumber &);	 // [space.cc]
  //
};
Qnumber compose    (const Qnumber &, const Qnumber &,
			     const double = 1.0);                // [space.cc]
//
//____________________________________________________________________________
class Qspace : public Space
{
  //
  //	Qspace class describes the quantum structure of a linear space
  //	divided into quantum subspaces.
  //	It is implemented as a derived class of Space adding quantum
  //	number informations
  //
private:
  //
  vector<Qnumber>	qs_number;
  //
public:
  //
  Qspace ();						     	 // [space.cc]
  Qspace (const Qspace &);					 // [space.cc]
  ~Qspace ();							 // [space.cc]
  //
  Qspace (const Qspace &, const Qspace &, 
	  const double = 1.0);       	                         // [space.cc]
  //
  inline size_t 	 id	(size_t)          const;       	 // [space.hh]
  inline void		 id 	(size_t, size_t);	       	 // [space.hh]
  inline const Qnumber & number (size_t) 	  const;         // [space.hh]
  inline double 	 lftlnk (size_t)	  const;	 // [space.hh]
  inline double		 rgtlnk (size_t)	  const;	 // [space.hh]
  //
  void show (const string & s = "") const;			 // [space.cc]
  //
  //	Setup
  //
  size_t add_states (const Qnumber &, size_t, long = 0);       	 // [space.cc]
  void   reflect    ();						 // [space.cc]
  //
private:
  //
  friend void quantum_select (Qspace &, const Qspace &);	 // [space.cc]
  //
};
//
//____________________________________________________________________________
class Qtarget : public Qspace
{
  //
  //	Qtarget class describes the quantum structure of the set of
  //	target states (spanning a linear space described by parent class
  //	Qspace) together with associated controls and results.
  //
private:
  //
  size_t		qt_sites;
  vector<double> 	qt_energy;	
  //
public:
  //
  Qtarget ();							 // [space.cc]
  Qtarget (const Qtarget &);					 // [space.cc]
  ~Qtarget ();							 // [space.cc]
  //
  //	Special ctors
  //
  Qtarget (const Qspace &, size_t);				 // [space.cc]
  Qtarget (const Qspace &, const Qspace &, size_t, 
	   const double = 1.0);                                  // [space.cc]
  //
  inline size_t 	  sites  () const;	       		 // [space.hh]
  inline vector<double> & energy ();				 // [space.hh]
  inline double		  energy (size_t) const;	      	 // [space.hh]
  void			  energy (double *, size_t, size_t);	 // [space.cc]
  //
};
//
//============================================================================
inline size_t Qnumber::id () const
{
  //
  //	Returns identifier
  //
  return qn_id;
}
//
//____________________________________________________________________________
inline void Qnumber::id (size_t ident)
{
  //
  //	Set identifier
  //
  qn_id = ident;
}
//
//============================================================================
inline	Subspace::Subspace ()
  : ss_offset		(0),
    ss_states		(0),
    ss_statistic	(0)
{ 
  //
  //	Default ctor.
  //
}
//
//============================================================================
inline size_t Space::offset (size_t index) const
{
  //
  //	Returns offset of first state for subspace index
  //
  return s_subspace [index] .ss_offset;
}
//
//____________________________________________________________________________
inline void Space::offsets ()
{
  //
  //	Updates the offsets of subspaces and counts states.
  //	This method must be called every time the order or the number 
  //	of states of some subspace changes.
  //
  size_t offset = 0;
  for (size_t i = 1; i < s_subspace .size (); i++) {
    s_subspace [i] .ss_offset = offset;
    offset += s_subspace [i] .ss_states;
  }
  s_subspace [0] .ss_states = offset;
}
//
//____________________________________________________________________________
inline const Subspace & Space::operator [] (size_t index) const
{
  //
  //	Array-like behaviour
  //
  return s_subspace [index];
}
//
//____________________________________________________________________________
inline size_t Space::states (size_t index) const
{
  //
  //	Get number of states of subspace index (space dimensions)
  //	or total number of states if index == 0.
  //
  return s_subspace [index] .ss_states;
}
//
//____________________________________________________________________________
inline void Space::states (size_t index, size_t states)
{
  //
  //	Set number of states for subspace index
  //
  s_subspace [index] .ss_states  = states;
  offsets ();
}
//
//____________________________________________________________________________
inline long Space::statistic (size_t index) const
{
  //
  //	Get statistic of states for subspace index
  //
  return s_subspace [index] .ss_statistic;
}
//
//____________________________________________________________________________
inline void Space::statistic (size_t index, long statistic)
{
  //
  //	Set statistic of states for subspace index
  //
  s_subspace [index] .ss_statistic = statistic;
}
//
//____________________________________________________________________________
inline const Subspace & Space::subspace (size_t index) const
{
  //
  //	Get subspace index
  //
  return s_subspace [index];
}
//
//____________________________________________________________________________
inline size_t Space::subspaces () const
{
  //
  //	Get number of subspaces (including null subspace)
  //
  return s_subspace .size ();
}
//
//============================================================================
inline size_t Qspace::id (size_t sub) const
{
  //
  //	identifier of subspace 
  //
  return qs_number [sub] .id ();
}
//
//____________________________________________________________________________
inline void Qspace::id (size_t sub, size_t identifier)
{
  //
  //	set identifier of subspace
  //
  qs_number [sub] .id (identifier);
}
//
//____________________________________________________________________________
inline double Qspace::lftlnk (size_t sub) const
{
  //
  //	Returns left link value
  //
  return qs_number [sub] .qn_lftlnk;
}
//
//____________________________________________________________________________
inline const Qnumber & Qspace::number (size_t sub) const
{
  //
  //	Returns subspace quantum number
  //
  return qs_number [sub];
}
//
//____________________________________________________________________________
inline double Qspace::rgtlnk (size_t sub) const
{
  //
  //	Returns left link value
  //
  return qs_number [sub] .qn_rgtlnk;
}
//
//============================================================================
inline vector<double> & Qtarget::energy ()
{
  //
  //	All energies
  //	
  return qt_energy;
}
//
//____________________________________________________________________________
inline double Qtarget::energy (size_t index) const
{
  //
  //	index-th energy
  //
  return qt_energy [index];
}
//
//____________________________________________________________________________
inline size_t Qtarget::sites () const
{
  //
  //	Returns n. of sites
  //
  return qt_sites;
}
//
//============================================================================
#endif // SPACE_HH
