//	algebra.hh			(C) 2009 Fabio Ortolani fbo 090716
//	==================================================================
#ifndef	ALGEBRA_HH
#define ALGEBRA_HH
#include <complex>
#include <string>
#include <vector>
using namespace std;
//
//============================================================================
struct	Afactor
{
  //
  //	Describes a single Action factor
  //	
  size_t	af_op;	// action name
  size_t	af_dg;	// dagger flag
  long		af_st;	// action index (we admit also negative indexes)
  //
  string  str () const;					       // [algebra.cc]
  //
  inline bool operator == (const Afactor &) const;	       // [algebra.hh]
  inline bool operator != (const Afactor &) const;	       // [algebra.hh]
  inline bool operator >  (const Afactor &) const;	       // [algebra.hh]
  inline bool operator <  (const Afactor &) const;	       // [algebra.hh]
  //
};
//
//____________________________________________________________________________
struct	Amono
{
  //
  //	a monomial of Afactor's
  //
  complex<double>	am_coeff;
  vector<Afactor>	am_factor;
  //
  Amono  ();						       // [algebra.cc]
  Amono  (const Amono &);      				       // [algebra.cc]
  Amono  (complex<double>);    				       // [algebra.cc]
  Amono  (const Afactor &);	       			       // [algebra.cc]
  ~Amono ();						       // [algebra.cc]
  //
  inline complex<double> coeff () const;		       // [algebra.hh]
  inline size_t 	 order () const;		       // [algebra.hh]
  inline size_t 	 size  () const;		       // [algebra.hh]
  //
  inline const Afactor & operator [] (size_t) const;	       // [algebra.hh]
  //
  inline bool operator == (const Amono &) const;	       // [algebra.hh]
  inline bool operator != (const Amono &) const;	       // [algebra.hh]
  inline bool operator >  (const Amono &) const;	       // [algebra.hh]
  inline bool operator <  (const Amono &) const;	       // [algebra.hh]
  //
  inline Amono & operator *= (complex<double>);		       // [algebra.hh]
  inline Amono & operator *= (const Afactor &);		       // [algebra.hh]
  inline Amono & operator *= (const Amono &);		       // [algebra.hh]
  //
  string str 	 () const;     				       // [algebra.cc]
  void	 reorder (size_t = 0);	       			       // [algebra.cc]
  //
};
//
//____________________________________________________________________________
struct Apoli
{
  //
  //	A formal polinomial of Amono
  //
  vector<Amono>	ap_mono;
  //
  Apoli  ();						       // [algebra.cc]
  Apoli  (const Apoli &);				       // [algebra.cc]
  Apoli	 (complex<double>);				       // [algebra.cc]
  Apoli	 (const Amono &);				       // [algebra.cc]
  ~Apoli ();						       // [algebra.cc]
  //
  inline size_t size () const;				       // [algebra.hh]
  //
  inline const Amono & operator [] (size_t) const;     	       // [algebra.hh]
  //
  inline bool operator == (const Apoli &) const;	       // [algebra.hh]
  //
  Apoli & operator += (const Apoli &);			       // [algebra.cc]
  Apoli & operator += (const Amono &);			       // [algebra.cc]
  Apoli & operator += (complex<double>);       		       // [algebra.cc]
  Apoli & operator += (const Afactor &);       		       // [algebra.cc]
  //
  Apoli & operator *= (complex<double>);		       // [algebra.cc]
  Apoli & operator *= (const Afactor &);	       	       // [algebra.cc]
  //
  Apoli & operator /= (complex<double>);		       // [algebra.cc]
  //
  size_t order   () const;				       // [algebra.cc]
  void   show    (const string & = "", size_t = 0) const;      // [algebra.cc]
  void	 reorder (size_t = 0);				       // [algebra.cc]
  //
};
//
//____________________________________________________________________________
struct Aproperty : Apoli
{
  //
  //	A formal property adding property specifications to a formal 
  //	polinomial Apoly
  //
  size_t		ap_sites;
  size_t		ap_id;
  long			ap_index;
  size_t		ap_braid;
  size_t		ap_braindex;
  size_t		ap_ketid;
  size_t		ap_ketindex;
  complex<double>	ap_value;
  //
  Aproperty ();						       // [algebra.cc]
  Aproperty (const Aproperty &);       			       // [algebra.cc]
  ~Aproperty ();					       // [algebra.cc]
  //
};
//
//============================================================================
bool		  boolean_expression (const string &);	       // [algebra.cc]
complex<double>	  complex_expression (const string &); 	       // [algebra.cc]
string		  complex_str 	     (const complex<double> &);// [algebra.cc]
double 		  double_expression  (const string &);         // [algebra.cc]
//
complex<double> * define_entry (const string &);       	       // [algebra.cc]
complex<double> * define_entry (const string &,
				const complex<double> &);      // [algebra.cc]
//
size_t 	       	name_action ();			       	       // [algebra.cc]
size_t 	       	name_action (const string &);  	       	       // [algebra.cc]
const string & 	name_action (size_t);   		       // [algebra.cc]
size_t 	       	name_parity ();			       	       // [algebra.cc]
size_t 	       	name_parity (const string &);  	       	       // [algebra.cc]
const string & 	name_parity (size_t);                 	       // [algebra.cc]
size_t 	       	name_quantum ();	       	       	       // [algebra.cc]
size_t 	       	name_quantum (const string &);         	       // [algebra.cc]
const string & 	name_quantum (size_t);                 	       // [algebra.cc]
//
size_t 	       	  name_define ();      		       	       // [algebra.cc]
size_t 	       	  name_define (const string &, 
			       complex<double> * = 0);         // [algebra.cc]
const string & 	  name_define (size_t);                	       // [algebra.cc]
complex<double> * name_value  (const string &, 
			       complex<double> * = 0);	       // [algebra.cc]
complex<double> * name_value  (size_t);			       // [algebra.cc] 
//
//============================================================================
inline bool Afactor::operator == (const Afactor & other) const
{
  //
  //	Compare Afactor's for equality
  //
  return ((af_op == other .af_op) &&
	  (af_dg == other .af_dg) &&
	  (af_st == other .af_st));
}
//
//____________________________________________________________________________
inline bool Afactor::operator != (const Afactor & other) const
{
  //
  //	Compare Afactor's for inequality
  //
  return ((af_op != other .af_op) ||
	  (af_dg != other .af_dg) ||
	  (af_st != other .af_st));
}
//
//____________________________________________________________________________
inline bool Afactor::operator > (const Afactor & other) const
{
  //
  //	Compare Afactor's for descending order
  //
  if (af_st > other .af_st) return true;
  if (af_st < other .af_st) return false;
  if (af_op > other .af_op) return true;
  if (af_op < other .af_op) return false;
  return (af_dg > other .af_dg);  // not dagger precedes dagger
}
//
//____________________________________________________________________________
inline bool Afactor::operator < (const Afactor & other) const
{
  //
  //	Compare Afactor's for ascending order
  //
  if (af_st < other .af_st) return true;
  if (af_st > other .af_st) return false;
  if (af_op < other .af_op) return true;
  if (af_op > other .af_op) return false;	
  return (af_dg < other .af_dg);  // not dagger precedes dagger
}
//
//============================================================================
inline complex<double> Amono::coeff () const
{
  //
  //	Get monomial coefficient
  //
  return am_coeff;
}
//
//____________________________________________________________________________
inline const Afactor & Amono::operator [] (size_t index) const
{
  //
  //	Get index-th Afactor
  //
  return am_factor [index];
}
//
//____________________________________________________________________________
inline bool Amono::operator == (const Amono & other) const
{
  //
  //	Compare Amono's factors
  //
  if (order () != other .order ()) return false;
  for (size_t m = 0; m < order (); ++m) 
    if (am_factor [m] != other .am_factor [m]) return false;
  return true;
}
//
//____________________________________________________________________________
inline bool Amono::operator != (const Amono & other) const
{
  //
  //	Compare Amono's factors
  //
  if (*this == other) return false;
  return true;
}
//
//____________________________________________________________________________
inline bool Amono::operator >  (const Amono & other) const
{
  //
  //	Compare Amono's factors
  //
  if (order () > other .order ()) return true;
  if (order () < other .order ()) return false;
  for (size_t m = 0; m < order (); ++m) {
    if (am_factor [m] > other .am_factor [m]) return true;
    if (am_factor [m] < other .am_factor [m]) return false;
  }
  return false;
}
//
//____________________________________________________________________________
inline bool Amono::operator <  (const Amono & other) const
{
  //
  //	Compare Amono's factors
  //
  if (order () < other .order ()) return true;
  if (order () > other .order ()) return false;
  for (size_t m = 0; m < order (); ++m) {
    if (am_factor [m] < other .am_factor [m]) return true;
    if (am_factor [m] > other .am_factor [m]) return false;
  }
  return false;
}
//
//____________________________________________________________________________
inline Amono & Amono::operator *= (complex<double> coeff)
{
  //
  //	Multiply by a scalar
  //
  am_coeff *= coeff;
  return *this;
}
//
//____________________________________________________________________________
inline Amono & Amono::operator *= (const Afactor & af)
{
  //
  //	Multiply by a Afactor
  //
  am_factor .push_back (af);
  reorder ();
  return *this;
}
//
//____________________________________________________________________________
inline Amono & Amono::operator *= (const Amono & other)
{
  //
  //	Multiply by a Amono
  //
  am_coeff *= other .am_coeff;
  for (size_t m = 0; m < other .order (); ++m) *this *= other [m];
  return *this;
}
//
//____________________________________________________________________________
inline size_t Amono::order () const
{
  //
  //	Get monomial order (n. of Afactor's)
  //
  return am_factor .size ();
}
//
//____________________________________________________________________________
inline size_t Amono::size () const
{
  //
  //	Get monomial order (n. of Afactor's)
  //
  return am_factor .size ();
}
//
//============================================================================
inline const Amono & Apoli::operator [] (size_t n) const
{
  //
  //	n-th monomial
  //
  return ap_mono [n];
} 
//
//____________________________________________________________________________
inline bool Apoli::operator == (const Apoli & other) const
{
  //
  //	Apoli's comparison
  //
  if (size () != other .size ()) return false;
  for (size_t m = 0; m < size (); ++m) {
    if (ap_mono [m] != other [m]) return false;
    if (ap_mono [m] .am_coeff != other [m] .am_coeff) return false;
  }
  return true;
}
//
//____________________________________________________________________________
inline size_t Apoli::size () const
{
  //
  //	Number of terms
  //
  return ap_mono .size ();
}
//
//============================================================================
#endif // ALGEBRA_HH
