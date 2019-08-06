//	input.cc			(C) 2009 Fabio Ortolani fbo 090716
//	==================================================================
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <iterator>
using namespace std;
//
//============================================================================
typedef vector<string> Slist;
//
//____________________________________________________________________________
void	          input_collect (int, char * []);    		 // [input.cc]
string &    	  input_next (bool = false);	       	       	 // [input.cc]
string &          input_peek ();	       			 // [input.cc]
Slist::iterator & input_position ();				 // [input.cc]
void              input_position (Slist::iterator);	       	 // [input.cc]
void 	    	  input_position (Slist::iterator, 
				  Slist::iterator);		 // [input.cc]
void        	  input_process ();      	       		 // [input.cc]
void	    	  input_show ();            	       		 // [input.cc]
void	 	  menu_parse (); 	     		          // [menu.cc]
static Slist &    input_clean (Slist &);       			 // [input.cc]
static bool 	  input_search (Slist &, Slist::iterator &, 
				string &);			 // [input.cc]
static void 	  input_split ();			       	 // [input.cc]
//
//____________________________________________________________________________
static Slist	       	input_list;
static Slist::iterator 	input_token;
static Slist::iterator 	input_end;
//
//============================================================================
static Slist &	input_clean (Slist & ilist)
{
  //
  //  Removes comments and trailing blanks from the input list of strings.
  //  A comment is delimited by '/*', '*/' (c-style comments) or started
  //  by '#' (script-style comments) or '//' (c++-style comments) which
  //  extends to end of string.
  //
  int c_comment = 0;
  for (Slist::iterator ln = ilist .begin (); ln != ilist .end (); ln++ ) {
    string & s = *ln;		// the line to be cleaned 
    string::size_type p;	// position in line
    bool comment = false;	// c++ or script style comment
    //
    //	scan input line for c++ or script style comments
    //
    for (p = 0; p < s. size ();  p++) {
      if (comment) {
	//
	// inside comment line (c++ or script style comment)
	//
        s[p] = ' ';	// clear to end of line
      }
      else if (s[p] == '/' && s [p+1] == '*') {
	//
	// beginning of c-style comment
	//
        c_comment++;	// mark and count nesting 
        s [p++] = ' ';	// start blanking comment
        s [p]   = ' ';
      }
      else if (c_comment) {
	//
	// inside c-style comment
	//
        if (s[p] == '*' && s[p+1] == '/') {
	  //
	  // closing c-style command
	  //
          c_comment--;		// update nesting
          s[p++] = ' ';		// blank comment end (star char)
        }
        s[p] = ' ' ;		// blank comment
      }
      else if (s [p] == '#' || ( s [p] == '/' && s [p+1] == '/') ) {
	//
	// beginning of c++ or script style command
	//
        s[p] = ' ';		
        comment = true;
      }
    }
    //
    //	Is there something in line?
    //
    p = s .size ();
    while (p && (s .at (p-1) <= ' ')) p--;
    if (p != s. size ()) s .erase (p);	// remove trailing blanks
    if (! p) {
      //
      // delete blank lines	
      //
      Slist::iterator del = ln;
      ln--;		// to continue loop
      ilist .erase (del);
    }
  }
  //
  //	Check for proper closing c-style comments
  //
  if (c_comment) {
    cout << "Error parsing input: unmatched /* ... */" << endl;
    exit (1);
  }
  return ilist;
}
//
//____________________________________________________________________________
void	input_collect (int argc, char * argv[])
{
  //
  //  Collects input from command line, then from standard input and from
  //  any other file specified via 'input' command.
  //
  string s;
  Slist	 ilist;		// temporary list
  Slist::iterator p;
  //
  //  First of all we initialize input from command line arguments
  //
  while (--argc) s = s + (++argv)[0] + ' ';	// read command line
  if (s .size())  ilist.push_back(s);
  ilist = input_clean (ilist);			// remove comments
  input_list .insert (input_list .end(), ilist .begin(), ilist .end());
  //
  //	reset temporary
  //
  ilist .clear();
  //
  //  Adds standard input (if no input command is found among arguments)
  //
  if (! input_search (input_list, p, s="input")) {
    while (getline (cin,s)) ilist .push_back (s);  // read stdin
    ilist = input_clean (ilist);			
    input_list .insert (input_list .end (), ilist .begin (), ilist .end ());
    ilist .clear ();
  }
  //
  //	recursively add input from other files
  //
  while (input_search (input_list, p, s = "input")) {
    //
    //	open file	
    //
    ifstream input (s .c_str ());
    while (getline (input,s)) ilist .push_back (s);	// read input
    input .close();
    ilist = input_clean (ilist);
    //
    //	insert input in place of include command
    //
    p = input_list .erase (p);		
    input_list .insert(p, ilist .begin (), ilist . end());
    ilist .clear();
  }
  //
  //  Search for output redirections (only last output command is active)
  //
  bool newout = false;
  string outfile;
  while (input_search (input_list, p, s = "output"))
  {
    //	
    //	remove output command and record filename
    //
    p = input_list .erase(p);
    newout = true;
    outfile = s;
  }
  if (newout)
  {
    //
    //	I don't know how to redirect output with c++ streams !
    //	So I use standard c routine
    //
    freopen (outfile .c_str(), "w", stdout);
  }
}
//
//____________________________________________________________________________
string & input_next (bool must_exist)
{
  //
  //	move to next token and return it.
  //	If at end of input and must_exist is true a fatal error is issued,
  //	else at end a null string is returned.
  //
  static string nulls;
  if (input_token != input_end) return *input_token++;
  if (must_exist) {
    cout << "Unexpected end of input found!" << endl;
    exit (0);
  }
  return nulls;
}
//
//____________________________________________________________________________
string & input_peek ()
{
  //
  //	return next input token without moving to next token 
  //	if at end of input returns a null string 
  //
  static string nulls;
  if (input_token != input_end) return *(input_token + 1);
  return nulls;
}
//
//____________________________________________________________________________
Slist::iterator & input_position ()
{
  //
  //	return actual input position 
  //
  return input_token;
}
//
//____________________________________________________________________________
void input_position (Slist::iterator pos)
{
  //
  //	set input position for next input 
  //
  input_token = pos;
}
//
//____________________________________________________________________________
void input_position (Slist::iterator pos, Slist::iterator endpos)
{
  //
  //	set input start and end position 
  //
  input_token = pos;
  input_end   = endpos;
}
//
//____________________________________________________________________________
void input_process ()
{
  //
  //	Splits input lines into tokens 
  //
  input_split ();
  //
  //	initialize default input start and end positions
  //
  input_position (input_list .begin (), input_list .end ());
  //
  //	Parse input
  //
  menu_parse ();
}
//
//____________________________________________________________________________
static bool input_search (Slist & ilist,
                          Slist::iterator & ln, string & name)
{
  //
  //	Searches the word name followed by a filename in the input list.
  //	If found adjust the input list in order to have the string
  //	"name filename" as a single element (line), set the iterator ln 
  //	pointing to this  element, set name to filename and returns true.
  //	If no valid command is found, returns false.
  //
  Slist::iterator n;
  string to_find (name);
  for (n = ilist .begin (); n != ilist .end (); ++n) {
    string::size_type p, e;
    string s = *n;
    if (s .find (to_find) == string::npos) continue;
    //
    while ((p = s .find (to_find)) != string::npos)  {
      //
      //  Check if the word name is isolated and not part of other words
      //
      s .erase (p, to_find .size ());
      if (p && (s [p-1] > ' '))	continue;  // something before
      if ((p < s .size ()) && (s [p] > ' '))  
	continue;      	// something after
      //
      // Looking for filename
      //
      while ((p < s .size ()) && (s[p] <= ' ')) p++;
      //
      // filename must follow command on the same line
      //
      if (p == s .size ()) {
	cout << "Missing filename for " << to_find << " command!" << endl;
	exit (1);
      }
      e = p;
      while ((e < s .size ()) && (s [e] > ' ')) e++;
      name = s .substr (p, e-p);          // extracting filename
      //
      //  Skip next blanks
      //
      while ((e < s .size ()) &&  (s [e] <= ' ')) e++;
      //
      //  Adjust list
      //
      Slist ins;
      if (p) ins .push_back (s .substr (0, p));	
      ins .push_back (to_find + " " + s .substr (p, e-p));
      if (e < s.size ()) ins .push_back (s .substr (e, s.size () - e));
      //
      ln = n; 
      ln--;
      n = ilist .erase (n);
      ilist .insert (n, ins .begin (), ins .end ());
      //
      //  Get new list position
      //
      while (ln ->compare (0, to_find .size (), to_find)) ln++;
      return true;
    }
  }
  return false;
}
//
//_____________________________________________________________________________
void input_show ()
{
  //
  //  Prints out the collected input (if any)
  //
  Slist::iterator sin = input_list.begin ();
  if (sin == input_list .end () ) return;
  cout << setw (26) << setfill ('_') << "" 
       << " Input (without comments) " << setw (26) << setfill ('_') << "" 
       << setfill (' ') << setw (0) << endl;
  cout << *sin << endl;
  while (++sin != input_list .end ()) cout << *sin << endl ;
  cout << setw(78) << setfill('_') << "" << setfill (' ') << endl;
  cout << flush;
}
//
//_____________________________________________________________________________
static void input_split ()
{
  //
  //	Splits input list of lines into a list of single string tokens
  //
  string token = "";
  Slist  token_list;
  bool nostring = true;
  bool noop   	= true;
  long parlevel = 0;
  for (Slist::iterator l = input_list .begin (); l != input_list .end (); ++l)	 {
    string & line = *l;
    for (string::size_type p = 0; p < line .size (); p++) {
      if (nostring && (line [p] == '(')) parlevel++;
      if (nostring && (line [p] == '[')) parlevel++;
      if (nostring && (line [p] == '{')) {
	parlevel++;
	continue;
      }
      if (nostring && (line [p] == '}')) { 
	parlevel--;
	continue;
      }
      if (nostring && (line [p] == ']')) parlevel--;
      if (nostring && (line [p] == ')')) parlevel--;
      if (line [p] == '"' ) {
	if (nostring) {
	  nostring = false;
	  parlevel++;
	}
	else {	
	  nostring = true;
	  parlevel--;
	}
      }
      else if ((! nostring) || (line [p] > ' ')) {
	token += line [p];
	noop = true;
	if (nostring && (line [p] == '+'))  noop = false;
	if (nostring && (line [p] == '-'))  noop = false;
	if (nostring && (line [p] == '*'))  noop = false;
	if (nostring && (line [p] == '/'))  noop = false;
	if (nostring && (line [p] == '%'))  noop = false;
	if (nostring && (line [p] == '!'))  noop = false;
	if (nostring && (line [p] == '='))  noop = false;
	if (nostring && (line [p] == '>'))  noop = false;
	if (nostring && (line [p] == '<'))  noop = false;
	if (nostring && (line [p] == '|'))  noop = false;
	if (nostring && (line [p] == '&'))  noop = false;
      }
      else if ((parlevel == 0) && noop && token .size ()) {
	token_list .push_back (token);
	token = "";
      }
    }
    if ((parlevel == 0) && noop && token .size ()) {
      token_list .push_back (token);
      token = "";
    }
  }
  if (parlevel) {
    cout << "Unmatched pars:" << endl << token << endl;
    exit (1);
  }
  if (token .size ()) {
    cout << "Incomplete expression: " << token << endl;
    exit (1);
  }
  input_list = token_list;
}
//
//============================================================================
