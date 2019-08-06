//	algebra.cc			(C) 2009 Fabio Ortolani fbo 090716
//	==================================================================
//
//	Manipulation of algebraic expressions, assignment and resolution 
//	of formal names.
//
//============================================================================
#include "algebra.hh"
#include "block.hh"
#include "numerical.hh"
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <iomanip>
//
//============================================================================
static bool            boolean_value (const string &);	       // [algebra.cc]
static complex<double> complex_value (const string &);	       // [algebra.cc]
static size_t 	       name_entry    (const string &, 
				      vector<string> &);       // [algebra.cc] 
static const string &  name_entry    (size_t, 
				      const vector<string> &); // [algebra.cc] 
//
static size_t end_factor (const string &, size_t); 	       // [algebra.cc] 
static size_t end_par    (const string &, size_t); 	       // [algebra.cc]
static size_t end_square (const string &, size_t); 	       // [algebra.cc]
//
//============================================================================
//
//	List of names for operators, quantum numbers, parity operators,
//	and variable names (with associated pointers to values)
//
static vector<string>			action_names;
static vector<string>			quantum_names;
static vector<string>			parity_names;
static vector<string>			define_names;
static vector<complex<double> *>       	define_values;
static string 				voidstr;
//
//============================================================================
extern Barray 				blocklft;
//
//============================================================================
string Afactor::str () const
{
  //
  //	return a string representing the operator
  //
  stringstream s;
  s << name_action (af_op);
  if (af_dg) s << "+";
  else 	     s << " ";
  s << "[" << af_st << "]";
  return s .str ();
}
//
//============================================================================
Amono::Amono ()
  : am_coeff	(0.0,0.0),
    am_factor	()
{ 
  //
  //	Default ctor.
  //
}
//
//____________________________________________________________________________
Amono::Amono (const Amono & other)
  : am_coeff 	(other .am_coeff),
    am_factor	(other .am_factor)
{ 
  //
  //	Copy ctor.
  //
}
//
//____________________________________________________________________________
Amono::Amono (complex<double> coeff)
  : am_coeff	(coeff),
    am_factor	()
{ 
  //
  //	Scalar monomial
  //
}
//
//____________________________________________________________________________
Amono::Amono (const Afactor & af)
  : am_coeff	(1.0,0.0),
    am_factor	(1, af)
{ 
  //
  //	Order 1 single factor monomial
  //
}
//
//____________________________________________________________________________
Amono::~Amono ()
{ 
  //
  //	Default dtor.
  //
}
//
//____________________________________________________________________________
void Amono::reorder (size_t modulus)
{
  //
  //	Reorder monomial factors taking into account Fermi-Bose statistic
  //	The operators are reordered in ascending order.
  //
  size_t i, j;
  Block & point = blocklft [1];
  //
  //	Fix indexes in range [0, modulus)
  //
  if (modulus)
    for (i = 0; i < order (); ++i) {
      long & site = am_factor [i] .af_st;
      while (site < 0) site += modulus;
      site %= modulus;
    }
  //
  //	Reorder Afactor's with respect to site index checking 
  //	Fermi-Bose statistic from primitive action definitions in one 
  //	site Block.
  //
  for (i = 1; i < order (); ++i) 
    if (am_factor [i-1] .af_st > am_factor [i] .af_st) {
      Afactor tmp = am_factor [i];
      long tstat  = point .action (tmp .af_op, 0) .statistic ();
      for (j = i; j && tmp .af_st < am_factor [j-1] .af_st; j--) {
	long jstat = point .action (am_factor [j-1] .af_op, 0) .statistic ();
	if (tstat < 0 && jstat < 0) am_coeff = -am_coeff;
	am_factor [j] = am_factor [j-1];
      }
      am_factor [j] = tmp;
    }
}
//
//____________________________________________________________________________
string Amono::str () const
{
  //
  //	Returns a string representation of the monomial factors
  //
  string s = "";
  if (order () == 0) s = "[Id]";
  for (size_t m = 0; m < order (); ++m) {
    if (m) s = s + " ";
    s = s + am_factor [m] .str ();
  }
  return s;
}
//
//============================================================================
Apoli::Apoli ()
  : ap_mono ()
{ 
  //
  //	Default ctor.
  //
}
//
//____________________________________________________________________________
Apoli::Apoli (const Apoli & other)
  : ap_mono (other .ap_mono)
{ 
  //
  //	Copy ctor.
  //
}
//
//____________________________________________________________________________
Apoli::Apoli (complex<double> number)
  : ap_mono (1, Amono(number))
{ 
  //
  //	Constant polinomial
  //
}
//
//____________________________________________________________________________
Apoli::Apoli (const Amono & mono)
  : ap_mono (1, mono)
{ 
  //
  //	Single moomial polinomial
  //
}
//
//____________________________________________________________________________
Apoli::~Apoli ()
{ 
  //
  //	Default dtor.
  //
}
//
//____________________________________________________________________________
Apoli & Apoli::operator += (const Amono & mono)
{
  //
  //	Adds a monomial to this polinomial (reordering monomial)
  //
  size_t insert;
  if (abs (mono. am_coeff) == 0.0) return *this;
  for (insert = 0; insert < size (); ++insert) {
    if (ap_mono [insert] == mono) {
      //
      //  Monomial present. Update coefficient.
      //
      ap_mono [insert] .am_coeff += mono .am_coeff;
      //
      // check for null coefficient
      //
      if (abs(ap_mono [insert] .am_coeff) == 0.0) 
	ap_mono .erase (ap_mono .begin () + insert);
      return *this;
    }
    if (ap_mono [insert] > mono) break;
  }
  //
  //	Insert monomial in right place
  //
  ap_mono .insert (ap_mono .begin () + insert, mono);
  return *this;
}
//
//____________________________________________________________________________
Apoli & Apoli::operator += (complex<double>  number)
{
  //
  //	Adds a scalar
  //
  *this += Amono (number);
  return *this;
}
//
//____________________________________________________________________________
Apoli & Apoli::operator += (const Afactor & af)
{
  //
  //	Adds a single factor
  //
  *this += Amono (af);
  return *this;
}
//
//____________________________________________________________________________
Apoli & Apoli::operator += (const Apoli & other)
{
  //
  //	Polinomial sum
  //
  for (size_t m = 0; m < other .size (); ++m) *this += other [m];
  return *this;
}
//
//____________________________________________________________________________
Apoli & Apoli::operator *= (complex<double> number)
{
  //
  //	Multiply by a constant
  //
  for (size_t m = 0; m < size (); ++m) ap_mono [m] *= number;
  return *this;
}
//
//____________________________________________________________________________
Apoli & Apoli::operator *= (const Afactor & af)
{
  //
  //	Multiply by a factor
  //
  for (size_t m = 0; m < size (); ++m) ap_mono [m] *= af;
  return *this;
}
//
//____________________________________________________________________________
Apoli & Apoli::operator /= (complex<double> number)
{
  //
  //	Divide by a constant.
  //
  if (number == 0.0) {
    cout << "Error: division by zero!" << endl;
    exit (0);
  }
  number = 1.0 / number;
  for (size_t m = 0; m < size (); ++m) ap_mono [m] *= number;
  return *this;
}
//
//____________________________________________________________________________
size_t Apoli::order () const
{
  //
  //	Returns the degree of the polinomial
  //
  size_t degree = 0;
  for (size_t m = 0; m < size (); ++m) 
    if (ap_mono [m] .order () > degree) degree = ap_mono [m] .order ();
  return degree;
}
//
//____________________________________________________________________________
void Apoli::reorder (size_t modulus)
{
  //
  //	reorder monomial operators (and whole polinomial)
  //
  Apoli next;
  for (size_t n = 0; n < size (); ++n) {
    ap_mono [n] .reorder (modulus);
    next += ap_mono [n];
  }
  *this = next;
}
//
//____________________________________________________________________________
void Apoli::show (const string & s, size_t offset) const
{
  //
  //	Output a representation of polinomial
  //
  if (s .size ()) cout << s << endl;
  for (size_t m = 0; m < size (); ++m)  {
    if (offset) cout << setw (offset) << "";
    string coeff = complex_str (ap_mono [m] .am_coeff) + " ";
    cout << setw (25) << right << coeff << ap_mono [m] .str () << endl;
  }
}
//
//============================================================================
Aproperty::Aproperty ()
  : Apoli	(),
    ap_sites	(0),
    ap_id	(0),
    ap_index	(0),
    ap_braid	(0),
    ap_braindex	(0),
    ap_ketid	(0),
    ap_ketindex	(0),
    ap_value	(0.0,0.0)
{ 
  //	
  //	Default ctor.
  //
}
//
//____________________________________________________________________________
Aproperty::Aproperty (const Aproperty & other)
  : Apoli	(other),
    ap_sites	(other .ap_sites),
    ap_id	(other .ap_id),
    ap_index	(other .ap_index),
    ap_braid	(other .ap_braid),
    ap_braindex	(other .ap_braindex),
    ap_ketid	(other .ap_ketid),
    ap_ketindex	(other .ap_ketindex),
    ap_value	(other .ap_value)
{ 
  //
  //	Copy ctor.
  //
}
//
//____________________________________________________________________________
Aproperty::~Aproperty ()
{
  //
  //	Default dtor.
  //
}
//
//============================================================================
//
//	Evaluation of formal expressions
//
//____________________________________________________________________________
bool boolean_expression (const string & s)
{
  //
  //	Parses a simple string for a logical expression returning a 
  //	boolean value.
  //	The expression is a sequence of propositions or subexpressions 
  //	(delimited by matched pars) separated by logical operators '||'
  //	(or) && (and). Every proposition can be a number (resolved to
  //	true if not null), the negation of a proposition (with a unary
  //	operator not !), or a comparison between two numbers (with binary 
  //	operators ==, !=, <, <=, >, >=). 
  //
  if (s .size () == 0) {
    cout << "Error: no logical expression" << endl;
    exit (0);
  }
  //
  //	The expression is a logical sum of logical products of factors.
  //
  //	The sum and partial product are initialized to zero (false)
  //
  bool	sum    = false;
  bool prod    = false;
  size_t start = 0;
  size_t end   = 0;
  while (start < s .size ()) {
    //
    //	Check logical operator (if any)
    //
    if ((s [start] == '|') && (s [start+1] == '|')) {
      //
      //  Accumulate in the sum the last product and initialize next
      //  product to true value
      //
      sum = sum || prod;
      prod = true;
      start += 2;
    }
    else if ((s [start] == '&') && (s [start+1] == '&'))
      //
      //  We must do a logical product with previous factors
      //  contained in prod
      //
      start += 2;
    else 
      //
      //  No logical operator starts the factor.  
      //  Initialize the product to true.
      //
      prod = true;
    //
    //	Check if the operation found is not at the beginning of
    //	expression
    //
    if (start && (end == 0)) {
      //
      //  Expression starts with a logical operator
      //
      cout << "Invalid unary logical operator " << string (s, end, start) 
	   << " in expression " << s << endl;
      exit (0);
    }
    //
    if (start >= end) {
      //
      //  Look for end of proposition or subexpression
      //
      end = start;
      while (end < s .size()) {
	//
	//  stop at a logical sum or product
	//
	if ((s [end] == '|') && (s [end+1] == '|')) break;
	if ((s [end] == '&') && (s [end+1] == '&')) break;
	if (s  [end] == ')') {
	  cout << "Unmatched pars in boolean expression " << s << endl;
	  exit (0);
	}
	if (s  [end] == ']') {
	  cout << "Unmatched square pars in boolean expression " << s << endl;
	  exit (0);
	}
	//
	//   Pars denote subexpressios (numerical or logical)
	//
	if (s [end] == '(')      end = end_par    (s, end+1);
	else if (s [end] == '[') end = end_square (s, end+1);
	end++;
      }
    }
    if (start >= end) {
      //
      //  Missing factor
      //
      cout << "Error parsng boolean expression " << s << endl;
      exit (0);
    }
    //
    //  Update product with proposition logical value  
    //
    prod = prod && boolean_value (string (s, start, end - start));
    //
    //	... and continue with next logical operation (if any)
    //
    start = end;
  }
  //
  //	return accumulated value (updated with last product)
  //
  return (sum || prod);
}
//
//____________________________________________________________________________
static bool boolean_value (const string & s)
{
  //
  //	evaluate string s representing a logical proposition 
  //	to a boolean value
  //
  if (s .size () == 0) {
    cout << "Error: null proposition!" << endl;
    exit (0);
  }
  //
  //	Simple values
  //
  if (s == "true") 	return true;
  if (s == "false") 	return false;
  //
  //	Negated expression
  //
  if ((s [0] == '!') && (s [1] != '!')) 
    return ! boolean_expression (string (s, 1, s .size () - 1));
  //
  //	Checks for logical subexpressions (enclosed in pars)
  //
  if ((s [0] == '(') && (end_par (s, 1) == s .size () - 1))
    return boolean_expression (string (s, 1, s .size () - 2));
  //
  //	Look for binary relation (if any)
  //
  int op 	= 0;	// type of relation
  int ops 	= 0;	// number of relations found (only one)
  size_t binary = 0;	// position of relation (size of left operand)
  size_t second = 0;	// start of right operand
  //
  size_t p = 0; 
  while (p < s .size ()) {
    //
    //	Check for comparison operations (only one comparison)
    //
    if (s [p] == '=') {
      op = -1;		// error if only one =
      ops++;	
      binary = p;
      p++;
      if (s [p] == '=') {
	op = 1;		// ==
	p++;
      }
      second = p;
    }
    else if (s [p] == '!') {
      op = -1;		// error if ! something else
      ops++;
      binary = p;
      p++;
      if (s [p] == '=') {
	op = 2;		// !=
	p++;
      }
      second = p;	
    }
    else if (s [p] == '>') {
      op = 3;		// >
      ops++;
      binary = p;
      p++;
      if (s [p] == '=') {
	op = 4;		// >=
	p++;
      }
      second = p;
    }
    else if (s [p] == '<') {
      op = 5;		// <
      ops++;
      binary = p;
      p++;
      if (s [p] == '=') {
	op = 6;		// <=
	p++;
      }
      second = p;
    }
    else if (s [p] == '(') p = end_par    (s, p+1) + 1;
    else if (s [p] == '[') p = end_square (s, p+1) + 1;
    else p++;
  }
  if ((ops > 1) || 
      (ops == 1 && ((binary == 0) || (second == s .size ()))) || 
      (op < 0)) {
    //
    //	We can't verify two consecutive comparisons 
    //	or a missing operand or an invalid operator
    //
    cout << "Invalid logical expression " << s << endl;
    exit (0);
  }
  //
  //	Evaluation
  //
  if (op) {
    //
    //	Comparison of two (real) values:
    //
    double a = double_expression (string (s, 0, binary));
    double b = double_expression (string (s, second, s .size () - second));
    //
    if (op == 1) return a == b;
    if (op == 2) return a != b;
    if (op == 3) return a >  b;
    if (op == 4) return a >= b;
    if (op == 5) return a <  b;
    if (op == 6) return a <= b;
  }
  //
  //	No logical operation found: 
  //
  return (abs (complex_expression (s)) != 0.0);
}
//
//____________________________________________________________________________
complex<double> complex_expression (const string & s)
{
  //
  //	Parses a simple string for an arithmetic expression returning a 
  //	complex value evaluation of the string expression.
  //	The expression is a sequence of numbers, variables and elementary
  //	functions (which are resolved to a complex value), separated by
  //	arithmetic operations ('+', '-', '*', '/'). No blank space is 
  //	allowed, and expression can contain matched pars () enclosing
  //	subexpressions or the argument (as a single subexpression) of 
  //	simple functions.
  //	The expression must be a sum of products of factors.
  //
  complex<double> sum	 = 0.0;
  complex<double> prod	 = 0.0;	
  complex<double> factor = 0.0;	
  size_t start = 0;
  size_t end   = 0;
  if (s .size () == 0) {	// error with null strings
    cout << "Error: no expression! " << endl;
    exit (0);
  }
  int	 operation  = 1;	// multiply 
  while (start < s.size ()) {
    //
    if (s [start] == '+') {
      //
      // unary or binary sum
      //
      sum += prod;	// accumulate previous value (initially zero)
      operation = 1;	// set operation to multiply
      prod	= 1.0;	// initialize product
      start++;
    }
    else if (s [start] == '-') {
      sum += prod;
      operation = 1;
      prod      = -1.0;
      start++;
    }
    else if (s [start] == '*') {
      if (operation) {
	//
	//	missing left operand
	//
	cout << "Error: missing left factor in expression " << s << endl;
	exit (0);
      }
      operation = 1;
      start++;
    }
    else if (s [start] == '/') {
      if (operation) {
	//
	// missing left operand
	//
	cout << "Error: missing dividend in expression " << s << endl;
	exit (0);
      }
      operation = 2;	// set operation to divide
      start++;
    }
    else if (s [start] == '%') {
      if (operation) {
	//
	// missing left operand	
	//
	cout << "Error: missing moduland in expression " << s << endl;
	exit (0);
      }
      operation = 3;	// set operation to module
      start++;
    }
    else {
      //
      //  initialize for a product if no unary operation provided
      //
      prod = 1.0;
      if (operation == 0) {
	//
	// missing binary operation
	//
	cout << "Error: missing binary operator in expression " << s << endl;
	exit (0);
      }
    }
    end = end_factor (s, start);
    if (end == start) {
      cout << "Error: invalid expression " << s << endl;
      exit (0);
    }
    factor = complex_value (string (s, start, end - start));
    if (operation == 1)  prod *= factor;
    else {
      //
      //  divide
      //
      if (factor == 0.0) {
	cout << "Error: division by zero in expression " << s << endl;
	exit (0);
      }
      if (operation == 2) prod /= factor;
      if (operation == 3) {
	if (factor .imag () || prod .imag ()) {
	  cout << "Error: invalid module factor in expression " << s << endl;
	  exit (0);
	}
	long module = factor .real ();
	long number = prod   .real ();
	prod = complex<double>(number % module);
      }
    }
    start = end;
    //
    //	Reset operation to none (next char must be a valid operation)
    //
    operation = 0;	// reset operation to none
  }
  return (sum + prod);
}
//
//____________________________________________________________________________
string	complex_str (const complex<double> & z)
{
  //
  //	Returns a string representation of a complex number
  //	
  stringstream s;
  double vr = z .real ();
  double vi = z .imag ();
  if ((vi == 0.0) || (vr != 0.0)) {
    if ((vr == 0.0) ||
	((fabs (vr) + 0.00000005 < 10.0) &&
	 (fabs (vr) + 0.00000005 > 0.01))) 
      s << setw (0) << setprecision (7) << fixed << showpos << vr;
    else s << setw (0) << setprecision (3) << scientific << showpos << vr;
  }
  if (vi != 0.0) {
    if ((fabs (vi) + 0.00000005 < 10.0) && 
	(fabs (vi) + 0.00000005 > 0.01))
      s << setw  (0) << setprecision (7) << showpos << fixed << vi << "*i";
    else s << setw (0) << setprecision (3) << showpos << scientific
	   << vi << "*i";
  }
  else s << "  ";
  return s .str ();
}
//
//____________________________________________________________________________
static complex<double>	complex_value (const string & s)
{
  //
  //	Resolve a variable name, function name, subexpression, to a
  //	complex value.	
  //
  //	String was checked for matching pars via end_factor, so 
  //	we are sure that an ending par is matched from first par 
  //	in the string),
  //
  complex<double> z;
  double eps = machine_precision () * 4.0;
  if (s [s .size () - 1] == ')') {
    //
    //	The string contains a subexpression or function argument
    //
    size_t par = 0;
    while (s [par] != '(') par++;
    z = complex_expression (string (s, par + 1, s .size () - par - 2));
    if (par) {
      //
      //	The string is a function expression
      //	Valid names are from STL available complex functions
      //
      string name (s, 0, par);
      if      (name  == "abs")   z = abs   (z);
      else if (name  == "arg" )  z = arg   (z);
      else if (name  == "conj")  z = conj  (z);
      else if (name  == "cos")   z = cos   (z);
      else if (name  == "cosh")  z = cosh  (z);
      else if (name  == "exp")   z = exp   (z);
      else if (name  == "log")   z = log   (z);
      else if (name  == "log10") z = log10 (z);
      else if (name  == "sin")   z = sin   (z);
      else if (name  == "sinh")  z = sinh  (z);
      else if (name  == "sqrt")  z = sqrt  (z);
      else if (name  == "tan")   z = tan   (z);
      else if (name  == "tanh")  z = tanh  (z);
      else if (name  == "integer") 	{int y = z.real(); z=y;}
      else if (name  == "random1m1")	{double y = (2*drand48() - 1); z=y; }
      else if (name  == "random01")		{double y = drand48(); z=y; } 
      else if (name  == "theta") {
	if (z .imag () )  z = 0.0 ;
	else if (z .real () > 0.0) z = 1.0;
	else z = 0.0;
      }			    
      else {
	cout << "Unknown function " << s 
	     << " (complex) ... (aborting)" << endl;
	exit (0);
      }
    }
  }
  else if (s [0] == '[') {
    //
    //	Square pars truncate to integer value
    //
    z = complex_expression (string (s, 1, s .size () - 2));
    long n = z .real ();
    long m = z .imag ();
    z = complex<double> (n,m);
  }
  else if (isalpha (s [0]) || (s [0] == '_')) {
    //
    //	string is a variable or array element (if string contains square
    //	pars)
    //
    complex<double> * vz = name_value (s);
    if (vz == 0) {
      cout << "Undefined variable " << s << endl;
      exit (0);
    }
    z = * vz;
  }
  else {
    stringstream val (s);
    double v;
    val >> v;
    z = v;
  }
  if (fabs (z .real ()) < eps) z = complex<double> (0.0, z .imag ());
  if (fabs (z .imag ()) < eps) z = complex<double> (z .real (), 0.0);
  return z;
}
//
//____________________________________________________________________________
double double_expression (const string & s)
{
  //
  //	Returns the real part evaluation of an algebric expression
  //
  return (complex_expression (s)) .real ();
}
//
//____________________________________________________________________________
static size_t end_factor (const string & s, size_t pos)
{
  //
  //	Returns position after the end of the factor starting in 
  //	position pos in string s.
  //
  if (pos >= s .size ()) return pos;
  if (isalpha (s [pos]) || (s [pos] == '(') || (s [pos] == '[')) {
    //
    //	The factor is a variable or function
    //
    while ((pos < s .size ()) &&
	   (isalnum (s [pos]) || (s [pos] == '_'))) pos++;
    if (s [pos] == '(')      pos = end_par    (s, pos+1) + 1;
    else if (s [pos] == '[') pos = end_square (s, pos+1) + 1;
    return pos;
  }
  if ((s [pos] == '0') && ((s [pos+1] == 'X') || (s [pos+1] == 'x'))) {
    //
    //	The factor is hexadecimal value
    //
    pos += 2;
    if (! isxdigit (s [pos])) {
      cout << "Invalid hexadecimal string: " << s << endl;
      exit (0);
    }
    while ((pos < s .size ()) && isxdigit (s [pos])) pos++;
    return pos;
  }
  //
  //	The factor is a numerical value
  //
  int digit = pos - 1;
  int start = pos;
  while ((pos < s .size ()) && isdigit (s [pos])) pos++;
  if (s [pos] == '.') digit = pos++;
  while ((pos < s .size ()) && isdigit (s [pos])) pos++;
  if ((digit > start) || (pos > (unsigned) (digit+1))) {
    //
    //	The factor represents a decimal value
    //
    if ((pos < s .size ()) && ((s [pos] == 'e') || (s [pos] == 'E'))) {
      //
      //  The factor represents a real value with exponent part
      //
      pos++;
      if ((s [pos] == '+') || (s [pos] == '-')) pos++;
      if (pos >= s .size () || (! isdigit (s [pos]))) {
	//
	//  missing exponetial digits
	//
	cout << "Error parsing numeric string " << s << endl;
	exit (0);
      }
      while ((pos < s .size ()) && isdigit (s [pos])) pos++;
    }
    return pos;
  }
  //
  //	the string contains an invalid expression
  //
  if ((s [pos] == ')') || (s [pos] == ']'))
    //
    //	easy error: unmatched par
    //
    cout << "Unmatched pars in expression " << s << endl;
  else cout << "Error parsing expression " << s << endl;
  exit (0);	
  return pos;
}
//
//____________________________________________________________________________
static size_t end_par (const string & s, size_t pos)
{
  //
  //	Returns position (>= pos) of first unbalanced closing par.
  //
  while (pos < s .size ()) {
    if (s [pos] == ')') return pos;
    else if (s [pos] == '(') pos = end_par (s, pos+1);
    else if (s [pos] == '[') pos = end_square (s, pos+1);
    pos++;
  }
  cout << "Error: unbalanced pars in expression " << s << endl;
  exit (0);
  return pos;
}
//
//____________________________________________________________________________
static size_t end_square (const string & s, size_t pos)
{
  //
  //	Returns position (>= pos) of first unbalanced closing square par.
  //
  while (pos < s .size ()) {
    if (s [pos] == ']') return pos;
    else if (s [pos] == '(') pos = end_par    (s, pos+1);
    else if (s [pos] == '[') pos = end_square (s, pos+1);
    pos++;
  }
  cout << "Error: unbalanced square pars in expression " << s << endl;
  exit (0);
  return pos;
}
//
//============================================================================
//
//	Formal variables defintion and resolution
//
//____________________________________________________________________________
complex<double> * define_entry (const string & s)
{
  //
  //	returns a pointer correspoding to name s (allocating space 
  //	if undefined)
  //
  complex<double> * v = name_value (s);
  if (v == 0) v = name_value (s, new complex<double>);
  return v;
}
//
//____________________________________________________________________________
complex<double> * define_entry (const string & s, const complex<double> & z)
{
  //
  //	returns a pointer correspoding to name s (allocating space 
  //	if undefined) and assign value to z
  //
  complex<double> * v = define_entry (s);
  v [0] = z;
  return v;
}
//
//____________________________________________________________________________
size_t name_action ()
{
  //
  //	Return total number of action names
  //
  return action_names .size ();
}
//
//____________________________________________________________________________
size_t name_action (const string & name)
{
  //
  //	Return an index for an action name (defining a new action name
  //	if needed)
  //
  return name_entry (name, action_names);
}
//
//____________________________________________________________________________
const string & name_action (size_t index)
{
  //
  //	Return the string name for action index
  //
  return name_entry (index, action_names);
}
//
//____________________________________________________________________________
size_t name_define ()
{
  //
  //	Return total number of defined names
  //
  return define_names .size ();
}
//
//____________________________________________________________________________
size_t name_define (const string & name, complex<double> * v)
{
  //
  //	Return an index for a defined name (initializing a new definition
  //	if needed).
  //	The define_names list is associated with a list of pointers to 
  //	complex values and we must be sure the two lists are syncronized
  //
  //	If the name contains a square open par (matching a closed square 
  //	at the end of the string) the string indicates an array and the 
  //	content of the squares is resolved to an index value)
  //
  string aname = name;
  //
  //	look for a open square. 
  //
  size_t square = 0;
  while ((square < name .size ()) && (name [square] != '[')) square++;
  if (square < name .size ()) {
    size_t end = end_square (name, square + 1) + 1;
    if ((square == 0) || (end != name .size ())) {
      cout << "Invalid name definition:" << name << endl;
      exit (1);
    }
    if ((square < name .size ()) && (name [square - 1] != ' ')) {
      //
      //  Resolve square par content to an index and redefine name
      //  with added blank space before open square to denote that
      //  index resolution was done.
      //
      stringstream s;
      complex<double> z = 
      complex_expression (string (name, square, end - square)); 
      if (z .real ()) s << setprecision (0) << z .real ();
      if (z .imag ()) s << " " << setprecision (0) << z .imag ();
      aname = string (name, 0, square) + " [" + s .str () + "]";
    }
  }
  //
  //	Get or set anane in the list of names 
  //
  size_t old   = name_define ();
  size_t found = name_entry  (aname, define_names);
  while (old <= found) {
    //
    //	the list of pointers is synchronized with define_names list
    //	adding null pointers
    //
    define_values .push_back (0);
    old++;
  }
  if (v) define_values [found] = v;
  return found;
}
//
//____________________________________________________________________________
const string & name_define (size_t index)
{
  //
  //	Return the string name for definition index
  //
  return name_entry (index, define_names);
}
//
//____________________________________________________________________________
static size_t name_entry (const string & name, vector<string> & namelist)
{
  //
  //	Look for name in the list and return the index of found entry
  //	If the name is not found it is added to the list
  //
  size_t i;
  for (i = 0; i < namelist .size (); ++i) 
    if (namelist [i] == name) return i;
  namelist .push_back (name);	// add name to list
  return i;
}
//
//____________________________________________________________________________
static const string & name_entry (size_t index, 
				  const vector<string> & namelist)
{
  //
  //	Return the string indexed by index in the list or a null string
  // 	if not found
  //
  if (index < namelist .size ()) return namelist [index];
  return voidstr;
}
//
//____________________________________________________________________________
size_t name_parity ()
{
  //
  //	Return total number of action names
  //
  return parity_names .size ();
}
//
//____________________________________________________________________________
size_t name_parity (const string & name)
{
  //
  //	Return an index for an action name (defining a new action name
  //	if needed)
  //
  return name_entry (name, parity_names);
}
//
//____________________________________________________________________________
const string & name_parity (size_t index)
{
  //
  //	Return the string name for action index
  //
  return name_entry (index, parity_names);
}
//
//____________________________________________________________________________
size_t name_quantum ()
{
  //
  //	Return total number of quamtum names
  //
  return quantum_names .size ();
}
//
//____________________________________________________________________________
size_t name_quantum (const string & name)
{
  //
  //	Return an index for an quantum name (defining a new quantum name
  //	if needed)
  //
  return name_entry (name, quantum_names);
}
//
//____________________________________________________________________________
const string & name_quantum (size_t index)
{
  //
  //	Return the string name for quantum index
  //
  return name_entry (index, quantum_names);
}
//
//____________________________________________________________________________
complex<double> * name_value (const string & name, complex<double> * v)
{
  //
  //	Define or set variable name and return a pointer to a complex
  //	associated to name. If v is not null it is associated to name.
  //	If v is null previous associated pointer is returned (marking 
  //	with a null pointer an undefined variable).
  //
  size_t found   = name_define  (name, v);
  return define_values [found]; 
}
//
//____________________________________________________________________________
complex<double> * name_value (size_t index)
{
  //
  //	Return the pointer to a double associated to index-th name 
  //	if defined.
  //
  if (index < define_values .size ()) return define_values [index];
  return 0;
}
//
//============================================================================
