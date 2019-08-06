//	space.cc			(C) 2009 Fabio Ortolani fbo 090716
//	==================================================================
#include "algebra.hh"
#include "space.hh"
#include <cstdlib>
#include <iostream>
#include <iomanip>
//
//============================================================================
bool quantum_alternate (size_t);				  // [menu.cc]
//
//============================================================================
Qnumber::Qnumber ()
  : qn_id	(0),
    qn_value	(),
    qn_parity	(),
    qn_lftlnk	(0),
    qn_rgtlnk	(0)
{ 
  //
  //	Default ctor.
  //
}
//
//____________________________________________________________________________
Qnumber::Qnumber (const Qnumber & other)
  : qn_id	(other .qn_id),
    qn_value	(other .qn_value),
    qn_parity	(other .qn_parity),
    qn_lftlnk	(other .qn_lftlnk),
    qn_rgtlnk	(other .qn_rgtlnk)
{ 
  //
  //	Copy ctor.
  //
}
//
//____________________________________________________________________________
Qnumber::Qnumber (size_t id, const vector<double> & qv,
		  const vector<complex<double> >  & qp,
		  double ql, double qr)
  : qn_id	(id),
    qn_value	(qv),
    qn_parity	(qp),
    qn_lftlnk	(ql),
    qn_rgtlnk 	(qr)
{ 
  //
  //	Initialized ctor.
  //
}
//
//____________________________________________________________________________
Qnumber::~Qnumber ()
{ 
  //
  //	Default ctor.
  //
}
//
//____________________________________________________________________________
bool Qnumber::operator == (const Qnumber & other) const
{
  //
  //	Compare with other
  //
  for (size_t i = 0; i < qn_value .size (); i++)
    if (qn_value [i] != other .qn_value [i]) return false;
  for (size_t i = 0; i < qn_parity .size (); i++)
    if (qn_parity [i] != other .qn_parity [i]) return false;
  if (qn_lftlnk != other .qn_lftlnk)	return false;
  if (qn_rgtlnk != other .qn_rgtlnk)	return false;
  
  return true;
}
//
//____________________________________________________________________________
void Qnumber::reflect ()
{
  //
  //	Reflect (change of sign) alternate quantum numbers
  //
  for (size_t i = 0; i < qn_value .size (); i++) 
    if (quantum_alternate (i)) qn_value [i] = -qn_value [i];
  double tmp = qn_lftlnk;
  qn_lftlnk = qn_rgtlnk;
  qn_rgtlnk = tmp;
}
//
//____________________________________________________________________________
string Qnumber::str () const
{
  //
  //	Returns a string with quantum values
  //
  stringstream s;
  size_t i;
  if (qn_id) s << name_define (qn_id) << " ";
  s << "(" << qn_lftlnk << ";";
  for (i = 0; i < qn_value .size (); i++) { 
    s << qn_value [i];
    if (i < qn_value .size () - 1) s << " ";
  }
  if (qn_parity .size ()) {
    if (qn_value .size ()) s << ";";
    for (i = 0; i < qn_parity .size (); i++) {
      if (imag (qn_parity [i])) {
	if (real (qn_parity [i])) s << real (qn_parity [i]) << "+ ";
	s << "*i" << imag (qn_parity [i]);
      }
      else s << real (qn_parity [i]);
    }
    if (i < qn_parity .size () - 1) s << " ";
  }
  s << ";" << qn_rgtlnk << ")";
  return s .str ();
}
//
//=============================================================================
Space::Space ()
  : s_subspace	(1)
{ 
  //
  //	Default ctor
  //
}
//
//____________________________________________________________________________
Space::Space (const Space & other)
  : s_subspace	(other .s_subspace)
{ 
  //
  //	Copy ctor.
  //
}
//
//____________________________________________________________________________
Space::~Space ()
{ 
  //
  //	Default dtor.
  //
}
//
//____________________________________________________________________________
bool Space::operator == (const Space & other) const
{
  //
  //	Compare two Space's
  //
  if (subspaces () != other .subspaces ()) return false;
  for (size_t i = 0; i < subspaces (); i++) {
    if (s_subspace [i] .ss_offset    != other [i] .ss_offset)    return false;
    if (s_subspace [i] .ss_states    != other [i] .ss_states)    return false;
    if (s_subspace [i] .ss_statistic != other [i] .ss_statistic) return false;
  }
  return true;
}
//
//============================================================================
Qspace::Qspace ()
  : Space     (),
    qs_number (1)
{ 
  //
  //	Default ctor.
  //
}
//
//____________________________________________________________________________
Qspace::Qspace (const Qspace & other)
  : Space     (other),
    qs_number (other .qs_number)
{ 
  //
  //	Copy ctor.
  //
}
//
//____________________________________________________________________________
Qspace::~Qspace ()
{ 
  //
  //	Default dtor.
  //
}
//
//____________________________________________________________________________
Qspace::Qspace (const Qspace & lft, const Qspace & rgt, const double qa)
  : Space (),
    qs_number (1)
{
  //
  //	Create a Qspace from the tensor product of lft and rgt Qspace's
  //	The order of inserted states does not matter
  //
  for (size_t rb = 1; rb < rgt .subspaces (); rb++)
    for (size_t lb = 1; lb < lft .subspaces (); lb++)
      if (lft .rgtlnk (lb) == rgt .lftlnk (rb)) {
	Qnumber  sum     = compose (lft .number (lb), rgt .number (rb), qa);
	size_t states    = lft .states    (lb) * rgt .states    (rb);
	size_t statistic = lft .statistic (lb) * rgt .statistic (rb);
	add_states (sum, states, statistic);
      }  
}
//
//____________________________________________________________________________
size_t Qspace::add_states (const Qnumber & qn, size_t n, long fb)
{
  //
  //	Adds n states with given fermi-bose statistic fb and quantum 
  //	numbers qn, and returns the insertion position of new states
  //
  size_t sub, position=0;
  int    diff = 1;
  for (sub = 1; sub < subspaces (); sub++) {
    diff = qn - qs_number [sub];
    if (diff <= 0) break;
    position += states (sub);
  }
  if (diff) {
    //
    //	Quantum nubers not found. Insert a new empty subspace
    //
    s_subspace .insert (s_subspace .begin () + sub, Subspace ());
    qs_number  .insert (qs_number  .begin () + sub, qn);
    statistic (sub ,fb);
    offsets   ();
  }
  position = offset (sub) + states (sub);
  states (sub, states (sub) + n);
  if (fb && statistic (sub) != fb) {
    cout << "Added states has wrong fermi-bose statistic!" << endl;
    exit (0);
  }
  return position;
}
//
//____________________________________________________________________________
void Qspace::reflect ()
{
  //
  //	Reflect (change of sign) alternate quantum numbers
  //
  for (size_t sub = 1; sub < subspaces (); sub++) 
    qs_number [sub] .reflect ();
}
//
//____________________________________________________________________________
void Qspace::show (const string & s) const
{
  //
  //	Show subspaces details
  //
  for (size_t i = 1; i < subspaces (); i++) {
    if (s .size ()) cout << s;
    cout << setw (38 -s .size()) << left << number (i) .str () << " " 
	 << "offset=" << setw (5) << right << offset (i) << " "
	 << "states=" << setw (5) << right << states (i) << " "
	 << "statistic= " << setw (2) << right << statistic (i) 
	 << setw (0) << endl;
  }
}
//
//============================================================================
Qnumber compose (const Qnumber & a, const Qnumber & b, double qa)
{
  //
  //	Composes two quantum numbers as a linear combination of additive values 
  //    and mutiplying parities.	
  //	The two Qnumber's are assumed compatible as arrays.
  //
  Qnumber c (a);
  c .qn_id = 0;
  for (size_t i = 0; i < c .qn_value .size (); i++) { 
    double cb = 1.0;
    if (quantum_alternate (i)) cb = qa;
    c .qn_value [i] +=  cb * b .qn_value [i];
  }
  for (size_t i = 0; i < c .qn_parity .size (); i++) 
    c .qn_parity [i] *= b .qn_parity [i];
  c .qn_rgtlnk = b .qn_rgtlnk;
  return c;
}
//____________________________________________________________________________
Qnumber operator + (const Qnumber & a, const Qnumber & b)
{
  //
  //	Composes two quantum numbers adding additive values and mutiplying
  //	parities.	
  //	The two Qnumber's are assumed compatible as arrays.
  //
  Qnumber c (a);
  c .qn_id = 0;
  for (size_t i = 0; i < c .qn_value .size (); i++) 
    c .qn_value [i] += b .qn_value [i];
  for (size_t i = 0; i < c .qn_parity .size (); i++) 
    c .qn_parity [i] *= b .qn_parity [i];
  c .qn_lftlnk = a .qn_lftlnk;
  c .qn_rgtlnk = b .qn_rgtlnk;
  return c;
}
//____________________________________________________________________________
int operator - (const Qnumber & a, const Qnumber & b)
{
  //
  //	Compares two Qnumber's a and b.
  //	Returns	0 if identical, -1 if a precedes b, +1 otherwise.
  //
  size_t i, j;
  j = a .qn_value .size ();
  if (j > b .qn_value .size ()) return  1;
  if (j < b .qn_value .size ()) return -1;
  for (i = 0; i < j; i++) {
    if (a .qn_value [i] > b .qn_value [i]) return  1;
    if (a .qn_value [i] < b .qn_value [i]) return -1;
  }
  j = a .qn_parity .size ();
  if (j > b .qn_parity .size ()) return -1;
  if (j < b .qn_parity .size ()) return  1;
  for (i = 0; i < j; i++) {
    if (real (a .qn_parity [i]) > real (b .qn_parity [i])) return  1;
    if (real (a .qn_parity [i]) < real (b .qn_parity [i])) return -1;
    if (imag (a .qn_parity [i]) > imag (b .qn_parity [i])) return  1;
    if (imag (a .qn_parity [i]) < imag (b .qn_parity [i])) return -1;
  }
  if (a .qn_lftlnk > b .qn_lftlnk) return  1;
  if (a .qn_lftlnk < b .qn_lftlnk) return -1;
  if (a .qn_rgtlnk > b .qn_rgtlnk) return  1;
  if (a .qn_rgtlnk < b .qn_rgtlnk) return -1;
  return 0;
}
//
//____________________________________________________________________________
void quantum_select (Qspace & quantum, const Qspace & wanted)
{
  //
  //	Selects in quantum the subspaces with quantum numbers nearest 
  //	to wanted subspaces.
  //	Subspaces in quantum are assumed ordered (with some order).
  //	On return quantum contains only nearest to wanted subspaces.
  //
  size_t sub;
  //
  //	Clear subspaces
  //
  for (sub = 1; sub < quantum .subspaces (); sub++) {
    quantum .id     (sub, 0);
    quantum .states (sub, 0);
  }
  //
  //	Scan wanted subspaces
  //
  for (sub = 1; sub < wanted .subspaces (); sub++) {
    Qnumber qw    = wanted .number (sub);
    size_t idw    = wanted .number (sub) .id ();
    size_t states = wanted .states (sub);
    //
    //	Find the nearest quantum number values (not equal to qw)
    //
    size_t ilow   = 0;
    size_t igreat = quantum .subspaces ();
    for (size_t ib = 1; ib < quantum .subspaces (); ib++) {
      int diff = qw - quantum .number (ib);
      if (diff < 0 && ib < igreat) igreat = ib;
      if (diff > 0 && ib > ilow)   ilow   = ib;
    }
    if (igreat - ilow == 1) {
      if (igreat < quantum .subspaces ()) igreat++;
    }
    else ilow++;
    if (ilow == 0) ilow++;
    for (size_t ib = ilow; ib < igreat; ib++) {
      quantum .id (ib, idw);
      if (quantum .states (ib) < states) quantum .states (ib, states);
    }
  }
  Qspace reduced;
  for (sub = 1; sub < quantum .subspaces (); sub++) 
    if (quantum .states (sub))
      reduced .add_states (quantum .number (sub),
			   quantum .states (sub), 
			   quantum .statistic (sub));
  quantum = reduced;
}
//
//============================================================================
Qtarget::Qtarget ()
  : Qspace 	(),
    qt_sites	(0),
    qt_energy	()
{ 
  //
  //	Default ctor.
  //
}
//
//____________________________________________________________________________
Qtarget::Qtarget (const Qtarget & other)
  : Qspace	(other),
    qt_sites	(other .qt_sites),
    qt_energy	(other .qt_energy)
{ 
  //
  //	Copy ctor.
  //
}
//
//____________________________________________________________________________
Qtarget::~Qtarget ()
{ 
  //
  //	Default dtor.
  //
}
//
//____________________________________________________________________________
Qtarget::Qtarget (const Qspace & qspace, size_t sites)
  : Qspace 	(qspace),
    qt_sites	(sites),
    qt_energy	()
{
  //
  //	Create a Qtarget from qspace for sites-lattice
  //
}
//
//____________________________________________________________________________
Qtarget::Qtarget (const Qspace & lft, const Qspace & rgt, size_t sites, 
		  const double qa)
  : Qspace	(lft, rgt, qa),
    qt_sites	(sites),
    qt_energy	()
{ 
  //
  //	Create a Qtarget from tensor product of lft and rgt Qspace's
  //
}
//
//____________________________________________________________________________
void Qtarget::energy (double * value, size_t sub, size_t add)
{
  //
  //	Set add energies for subspace sub
  //
  states (sub, add);
  for (size_t j = 0; j < add; j++) 
    qt_energy .insert (qt_energy .begin () + offset (sub) + j, value [j]);
}
//
//============================================================================
