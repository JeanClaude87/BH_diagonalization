//	menu.cc				(C) 2011 Fabio Ortolani fbo 110201
//	==================================================================
#include "block.hh"
#include "algebra.hh"
#include <cstdlib>
#include <iostream>
#include <iomanip>
//
//============================================================================
struct	Menu
{
  //
  //	To every keyword associate a function (if not null) and a pointer
  //	as argument for the function
  //
  string	keyword;
  void (*	action) (void *);
  void *	data;
};
//
//____________________________________________________________________________
struct	Loop
{
  //
  //	Control structure for a loop, while or if
  //
  double *			l_index;
  double 			l_start;
  double			l_end;
  double			l_step;
  vector<string>::iterator	l_body;
  vector<string>::iterator	l_close;
};
//
//===========================================================================
void		ham_parse (size_t);				 // [menu.cc]
void            hamt_parse (size_t, double);                     // [menu.cc] 
string & 	input_next (bool = false);		       	// [input.cc]
vector<string>::iterator & input_position ();			// [input.cc]
void 		input_position (vector<string>::iterator);     	// [input.cc]
void 		input_position (vector<string>::iterator,
				vector<string>::iterator);     	// [input.cc]
static void	menu_close (const string &, const string &);	 // [menu.cc]
static Menu * 	menu_find (const string &, Menu *);		 // [menu.cc]
void		menu_parse ();					 // [menu.cc]
//
//	valid menu actions
//
static void	menu_action	(void *);			 // [menu.cc]
static void	menu_boolean   	(void *);      			 // [menu.cc]
static void	menu_complex	(void *);      			 // [menu.cc]
static void	menu_define	(void *);      			 // [menu.cc]
static void	menu_density	(void *);      			 // [menu.cc]
static void	menu_dmrgs	(void *);			 // [menu.cc]
static void	menu_dmrgstime	(void *);			 // [menu.cc]
static void	menu_do		(void *);	       		 // [menu.cc]
static void	menu_done	(void *);      			 // [menu.cc]
static void    	menu_endgroup	(void *);	       	       	 // [menu.cc]
static void    	menu_endwhile	(void *);	       	       	 // [menu.cc]
static void    	menu_evolution	(void *);	       	       	 // [menu.cc]
static void    	menu_false	(void *);      		       	 // [menu.cc]
static void    	menu_hamb      	(void *);	       	       	 // [menu.cc]
static void    	menu_hame      	(void *);	       	       	 // [menu.cc]
static void    	menu_hamtb     	(void *);	       	       	 // [menu.cc]
static void    	menu_hamte     	(void *);	       	       	 // [menu.cc]
static void    	menu_if		(void *);	       	       	 // [menu.cc]
static void    	menu_integer	(void *);	       	       	 // [menu.cc]
static void    	menu_memory	(void *);      		       	 // [menu.cc]
static void    	menu_name	(void *);      		       	 // [menu.cc]
static void    	menu_parse	(void *);		       	 // [menu.cc]
static void	menu_pro_begin 	(void *);      			 // [menu.cc]
static void	menu_pro_ibegin	(void *);      			 // [menu.cc]
static void	menu_pro_bra 	(void *);      			 // [menu.cc]
static void	menu_pro_ket	(void *);      			 // [menu.cc]
static void    	menu_pro_mono	(void *);      		       	 // [menu.cc]
static void    	menu_pro_op	(void *);      		       	 // [menu.cc]
static void	menu_real	(void *);			 // [menu.cc]
static void	menu_renyi	(void *);			 // [menu.cc]
static void    	menu_show_bool	(void *);		       	 // [menu.cc]
static void    	menu_show_endl	(void *);		       	 // [menu.cc]
static void    	menu_show_str	(void *);		       	 // [menu.cc]
static void    	menu_show_val	(void *);		       	 // [menu.cc]
static void    	menu_size	(void *);      		       	 // [menu.cc]
static void	menu_space	(void *);			 // [menu.cc]
static void	menu_space_a	(void *);			 // [menu.cc]
static void	menu_space_b	(void *);			 // [menu.cc]
static void	menu_space_f	(void *);			 // [menu.cc]
static void	menu_space_p	(void *);			 // [menu.cc]
static void	menu_space_q	(void *);			 // [menu.cc]
static void	menu_symmetry	(void *);			 // [menu.cc]
static void	menu_target	(void *);			 // [menu.cc]
static void    	menu_true	(void *);	       	       	 // [menu.cc]
static void     menu_ttarget    (void *);			 // [menu.cc]
static void    	menu_while   	(void *);      		       	 // [menu.cc]
//
static void	space_reset	();				 // [menu.cc]
//
//===========================================================================
extern size_t		length;					 // [dmrg.cc]
extern size_t		sites_step;				 // [dmrg.cc]
extern complex<double>	periodic;      				 // [dmrg.cc]
extern bool		reflect_blocks;				 // [dmrg.cc]
extern bool		reflect_universe;      			 // [dmrg.cc]
extern size_t		zips;					 // [dmrg.cc]
extern size_t		zipinner;      				 // [dmrg.cc]
extern bool		halfsweep;				 // [dmrg.cc]
extern size_t  	       	show_blocks;				 // [dmrg.cc]
extern size_t		show_density;				 // [dmrg.cc]
extern size_t 		show_hamiltonian;			 // [dmrg.cc]
extern size_t	       	show_interaction;			 // [dmrg.cc]
extern bool		show_point;				 // [dmrg.cc]
extern bool		show_memory;				 // [dmrg.cc]
extern size_t		show_system;				 // [dmrg.cc]
extern bool		show_target;				 // [dmrg.cc]
extern bool		show_states;				 // [dmrg.cc]
extern bool		show_selection;				 // [dmrg.cc]
//
extern Barray		blocklft;      				 // [dmrg.cc]
extern Barray		blockrgt;		       		 // [dmrg.cc]
//extern vector<Brule>	brule;					 // [dmrg.cc]
extern vector<Brule>	lftrule;			       	 // [dmrg.cc]
extern vector<Brule>	rgtrule;			       	 // [dmrg.cc]
extern vector<Brule>	lfttimerule;			       	 // [dmrg.cc]
extern vector<Brule>	rgttimerule;			       	 // [dmrg.cc]
extern vector<Qtarget>	target_request;				 // [dmrg.cc]
extern vector<size_t>	symmetry;				 // [dmrg.cc]
//
extern Apoli			hamiltonian;   		   // [superblock.cc]
extern Apoli                    timehamiltonian;           // [superblock.cc]
extern Aproperty                start_action;              // [superblock.cc]
extern size_t                   ttarget_id;		   // [superblock.cc]
extern size_t                   ttarget_index;             // [superblock.cc]
extern double                   timestep;		   // [superblock.cc]
extern size_t                   steps_number;		   // [superblock.cc]
extern size_t			step_divisions;            // [superblock.cc]
extern size_t                   time_zips;                 // [superblock.cc]
//
extern vector<Aproperty>  	correlation;		   // [superblock.cc]
extern vector<double>		renyiarg;		   // [superblock.cc]
extern vector<Aproperty>	densityoperator;	   // [superblock.cc]
//
extern double		zero_energy;			  // [superaction.cc]
extern double		zero_norm;			  // [superaction.cc]
extern size_t		lanczos_iters;			  // [superaction.cc]
extern size_t		lanczos_repeats;       		  // [superaction.cc]
extern size_t		lanczos_strategy;		  // [superaction.cc]
//
extern bool 		evolution_active;		         // [dmrg.cc]
extern double		tensorweight;			   // [superblock.cc]
//
//===========================================================================
static string			keyword;
static int     			loop_nest      	= -1;
static vector<Loop>		loop_list;
static size_t			qstates	       	= 0;
static long			qstatistic     	= 0;
static size_t			qvalues	       	= 0;
static size_t			qparities      	= 0;
static size_t			qid    		= 0;
static size_t			qsites	       	= 0;
static vector<double>		qvalue;
static vector<bool>             qalternate;
static vector<complex<double> > qparity;
static double			qlftlnk;
static double			qrgtlnk;
static size_t			site_states    	= 0;
static size_t			site_actions   	= 0;
static vector<string>		ham_definition;
static vector<string>		hamt_definition;
static vector<string>		densityaction;	
static string                   timestring; 
//
static Aproperty		property;
static Amono			pro_mono;
//
//____________________________________________________________________________
static Menu menu_common [] =
  {
    //
    //	keywords common to all menus
    //
    { "define",			menu_define,	0			},
    { "do",			menu_do,	0			},
    { "done",			menu_done,	0			},
    { "while",			menu_while,	0			},
    { "endwhile",		menu_endwhile,	0			},
    { "if",			menu_if,	0			},
    { "endif",			0,		0			},
    { "endgroup",		menu_endgroup,	0			},
    { "", 0, 0 }
  };
static Menu menu_monomial [] =
  {
    //
    //	Monomial entries for hamiltonian and properies
    //
    { "coeff",			menu_complex,	& pro_mono .am_coeff	},
    { "factor",			menu_complex,	& pro_mono .am_coeff	},
    { "action",			menu_pro_op,	0			},
    { "operator",      		menu_pro_op,  	0			},
    { "", 0, 0 }
  };
static Menu menu_property [] =
  {
    //
    //	Polinomial definition
    //
    { "term",			menu_pro_mono,	0			},
    { "add",			menu_pro_mono,	0			},
    //
    //	Property specs
    //
    { "sites",      		menu_size,  	& property .ap_sites  	},
    { "id",      		menu_name,  	& property .ap_id  	},
    { "name",      		menu_name,  	& property .ap_id  	},
    { "index",			menu_integer,	& property .ap_index	},
    { "bra",      		menu_pro_bra,	0			},
    { "bra_id",      		menu_name,  	& property .ap_braid  	},
    { "bra_index",    		menu_size,  	& property .ap_braindex },
    { "ket",      		menu_pro_ket,	0			},
    { "ket_id",      		menu_name,  	& property .ap_ketid  	},
    { "ket_index",    		menu_size,  	& property .ap_ketindex },    
    { "", 0, 0 }
  };
static Menu menu_lanczos [] =
  {
    //
    //	Lanczos algorithm specifications 
    //
    { "iterations",		menu_size,	& lanczos_iters		},
    { "repeatitions",		menu_size,      & lanczos_repeats	},
    { "strategy",		menu_size,	& lanczos_strategy	},
    { "zero_energy",		menu_real,	& zero_energy		},
    { "zero_norm",		menu_real,	& zero_norm		},
    { "", 0, 0 }
  };
static Menu menu_show [] =
  {
    //
    //	Output control 
    //
    { "nopoint",		menu_false,    	& show_point   		},
    { "point",			menu_true,    	& show_point   		},
    { "notarget",		menu_false,    	& show_target  		},
    { "target",			menu_true,    	& show_target  		},
    { "memory",			menu_true,	& show_memory		},
    { "nomemory",      		menu_false,	& show_memory		},
    { "states",			menu_true,	& show_states		},
    { "nostates",      		menu_false,	& show_states		},
    { "selection",		menu_true,	& show_selection	},
    { "noselection",		menu_false,	& show_selection	},
    { "blocks",			menu_size,	& show_blocks		},
    { "system",			menu_size,	& show_system		},
    { "density",       		menu_size,	& show_density		},
    { "hamiltonian",   		menu_size,	& show_hamiltonian     	},
    { "interaction",   		menu_size,	& show_interaction     	},
    { "boolean",	       	menu_show_bool,	0 },
    { "endl",			menu_show_endl,	0 },
    { "newline",       		menu_show_endl,	0 },
    { "string",			menu_show_str,	0 },
    { "value",			menu_show_val,	0 },
    { "", 0, 0 }
  };
static Menu menu_subspace [] = 
  {
    //	
    //	Subspace and target definition
    //
    //	Statistic
    //
    { "statistic",		menu_integer,	& qstatistic		},
    { "fermi",			menu_space_f, 	0			},
    { "fermionic",     		menu_space_f, 	0			},
    { "fermions",      		menu_space_f, 	0			},
    { "bose",			menu_space_b, 	0			},
    { "bosonic",     		menu_space_b, 	0			},
    { "bosons",      		menu_space_b, 	0			},
    //
    //	N. of states (subspace dimensions)
    //
    { "states",			menu_size,	& qstates		},
    { "dimensions",		menu_size,	& qstates		},
    { "dimension",		menu_size,	& qstates		},
    //
    //	Target active sites
    //
    { "sites",			menu_size,	& qsites		},
    //
    //	Quantum name and value
    //
    { "quantum",		menu_space_q,	0			},
    //
    //	Quantum name and value with alternate sum
    //
    { "alternate",		menu_space_a,	0			},
    { "qalternate",		menu_space_a,	0			},
    { "quantum_alternate",     	menu_space_a,	0			},
    //
    //	Parity name and value
    //
    { "parity",			menu_space_p,	0			},
    { "", 0, 0 }
  };    
//
//____________________________________________________________________________
static Menu menu_time [] =
  {
    //
    //	evolution DMRG algorithm
    //
    { "timestep",               menu_real,  	 & timestep              },
    { "steps_number",           menu_size,  	 & steps_number          },
    { "step_divisions",         menu_size,  	 & step_divisions        },
    { "time_zips",              menu_size,  	 & time_zips             },
    { "dmrg_timestates",	menu_dmrgstime,	 0			},
    { "timeham_begin",          menu_hamtb, 	 0                       },
    { "timeham_end",            menu_hamte, 	 0                       },
    { "target",                 menu_ttarget, 	 0                       },
    { "initial_state",          menu_ttarget, 	 0                       },
    { "initial_action",         menu_pro_ibegin, 0                       },
    { "initial_operator",       menu_pro_ibegin, 0                       },
    { "property",               menu_pro_begin,  0                       },
    { "", 0, 0 }          
  };
//____________________________________________________________________________
static Menu menu_main [] =
  {
    //
    //	Main menu.
    //	Input is parsed startig with keywords of this menu
    //
    //	Lattice specs and structure
    //
    { "sites",			menu_size, 	& length		},
    { "length",			menu_size, 	& length		},
    { "total_sites",	       	menu_size, 	& length		},
    { "sites_step",		menu_size,	& sites_step		},
    { "periodic",		menu_complex,	& periodic		},
    { "boundary",		menu_complex,	& periodic		},
    { "reflect_blocks",		menu_boolean,	& reflect_blocks	},
    { "reflect_universe",      	menu_boolean,	& reflect_universe	},
    //
    //	DMRG algorithm control
    //
    { "iterations",		menu_size,	& zips			},
    { "zips",			menu_size,	& zips			},
    { "zippers",		menu_size,	& zips			},
    { "zipinner",		menu_size,	& zipinner     		},
    { "innerzip",		menu_size,	& zipinner     		},
    { "half_sweep",		menu_boolean,	& halfsweep		},
    { "halfsweep",		menu_boolean,	& halfsweep		},
    { "half_zip",		menu_boolean,	& halfsweep		},
    { "halfzip",		menu_boolean,	& halfsweep		},
    { "block_states",		menu_dmrgs,	0			},
    { "dmrg_states",		menu_dmrgs,	0			},
    { "density_action",		menu_density,	0			},
    { "densityaction",		menu_density,	0			},
    //
    //	Lattice single site subspace definitions
    //
    { "space",			menu_space,	0			},
    { "subspace",      		menu_space,	0			},
    //
    // 	Target definition
    //
    { "target",			menu_target,	0			},
    { "hamweight",		menu_real,	& tensorweight		},
    //
    //	Lattice single site Action definition
    //
    { "action",			menu_action,	0			},
    { "operator",      		menu_action,	0			},
    { "symmetry",      		menu_symmetry,	0			},
    //
    //	Hamiltonian definition
    //
    { "ham_begin",		menu_hamb,     	0			},
    { "ham_end",		menu_hame,     	0		    	},
    //
    //	Property definition
    //
    { "property",		menu_pro_begin,	0			},
    { "correlation",		menu_pro_begin,	0			},
    { "renyi", 			menu_renyi,	0			},
    //
    //  Evolution flag
    //
    { "evolution",              menu_evolution, 0                       },
    { "timestep",               menu_real,      & timestep              },
    { "steps_number",           menu_size,      & steps_number          },
    { "step_divisions",         menu_size,  	& step_divisions        },
    { "time_zips",              menu_size,  	& time_zips             },
    { "dmrg_timestates",	menu_dmrgstime,	0			},
    //
    //	Memory management
    //
    { "memory",	            	menu_memory, 	0			},
    //
    //	Lanczos specs
    //
    { "lanczos",		menu_parse,	menu_lanczos		},
    //
    //	show specs
    //
    { "show",			menu_parse,	menu_show		},
    { "", 0, 0 }
  };
//
//============================================================================
void densityop_parse (size_t actualsys, size_t actualuni)
{
  //
  //	Parses operators to be applied to states for added
  //	reduced density matrices
  // 	
  vector<string>::iterator  start = densityaction .begin ();
  vector<string>::iterator  end   = densityaction .end ();
  densityoperator .clear ();
  if (start == end)	return;
  //
  //	Fix input position
  //
  input_position (start, end);
  keyword = input_next ();
  while (keyword .size ()) {
    if (keyword == "density_operator") {
      keyword = "";
      property = Aproperty ();
      //
      //   Set lattice sites variables
      //
      property .ap_sites = actualsys + actualuni;
      string syssites (input_next (true));
      define_entry (syssites, actualsys);
      string unisites (input_next (true));
      define_entry (unisites, actualuni);
      menu_name (& property .ap_id);
      menu_parse (menu_property);
      property .reorder (actualsys + actualuni);
      densityoperator .push_back (property);
    }
    else {
      cout << "Invalid keyword " << keyword << " parsing density operator."
	   << endl;
      exit (0);
    }
  }
}
//
//____________________________________________________________________________
void ham_parse (size_t actual_sites)
{
  //
  //	Parses lattice hamiltonian definitions according to menu_property
  //
  vector<string>::iterator  start = ham_definition .begin ();
  vector<string>::iterator  end   = ham_definition .end ();
  if (start == end)	return;
  //
  //	Fix input position
  //
  input_position (start, end);
  //
  //	First word is variable name for n. of sites
  //
  string nsites (input_next (true));
  define_entry (nsites, actual_sites);
  //
  //	Clear previous hamiltonian definitions (if any) and initialize 
  //	lattice sites 
  //
  property = Aproperty ();
  property .ap_sites = actual_sites;
  //
  //	Parse next input as a general property but ignore all
  //	property specs, only polinomial expression is used
  //	
  menu_parse (menu_property); 
  if (keyword .size ()) {
    cout << "Parsing hamiltonian error: keyword " << keyword << " unknown!"
	 << endl;
    exit (0);
  }
  hamiltonian = property;
}
//
//____________________________________________________________________________
void hamt_parse (size_t actual_sites, double actual_time)
{
  //
  //	Parses lattice hamiltonian definitions according to menu_property
  //
  vector<string>::iterator  start = hamt_definition .begin ();
  vector<string>::iterator  end   = hamt_definition .end ();
  if (start == end)	return;
  //
  //	Fix input position
  //
  input_position (start, end);
  //
  //	First word is variable name for n. of sites
  //
  string nsites (input_next (true));
  define_entry (nsites, actual_sites);
  define_entry (timestring, actual_time);
  //
  //	Clear previous hamiltonian definitions (if any) and initialize 
  //	lattice sites 
  //
  property = Aproperty ();
  property .ap_sites = actual_sites;
  //
  //	Parse next input as a general property but ignore all
  //	property specs, only polinomial expression is used
  //	
  menu_parse (menu_property); 
  if (keyword .size ()) {
    cout << "Parsing hamiltonian error: keyword " << keyword << " unknown!"
	 << endl;
    exit (0);
  }
  hamiltonian = property;
}
//
//____________________________________________________________________________
bool quantum_alternate (size_t index)
{
  return qalternate [index];
}
//
//____________________________________________________________________________
bool quantum_alternate ()
{
  for (size_t i = 0; i < qalternate .size (); ++i)
    if (qalternate [i]) return true;
  return false;
}
//
//____________________________________________________________________________
static void menu_action (void *)
{
  //
  //	Defines a one-site lattice block operator
  //
  if (site_states == 0) {
    cout << "An operator must be defined after space definition!" << endl;
    exit (0);
  }
  //
  //	Get action name and values (as a full complex matrix)
  //
  string aname (input_next (true));
  size_t states = blocklft [1] .states ();
  const vector<size_t> & order = blocklft [1] .tensororder ();
  complex<double> * m = new complex<double> [states * states];
  //
  //	Get values and reorder them according to internal order
  //
  for (size_t i = 0; i < states; ++i) 
    for (size_t j = 0; j < states; ++j) 
      menu_complex (m + order[i] + order[j] * states);
  //
  //	Build the action
  //
  Action op (blocklft [1] .quantum (), blocklft [1] .quantum ());
  op .compress (m);
  delete [] m;
  if (op .statistic () == 0) 
    cout << "Warning: Action " << aname
	 << " does not have a Fermi-Bose statistic, treated as bosonic!"
	 << endl;
  //
  //	Store op inside point Block
  //
  size_t index = name_action (aname);
  blocklft [1] .action (index, 0) << op;
  blockrgt [1] .action (index, 0) << op;
  site_actions = index + 1;
}
//
//____________________________________________________________________________
static void menu_boolean (void * data)
{
  //
  //	Gets boolean value and assigns to *data
  //
  bool & v = ((bool *) data) [0];
  v = boolean_expression (input_next (true));
}
//
//____________________________________________________________________________
static void menu_close (const string & close, const string & open)
{
  //
  //	Move to a closing control keyword close (matching previous 
  //	keyword open), with possible nestings of opening and closing
  //	control keys.
  //
  string s;
  while ((s = input_next (false)) .size ()) {
    if (s == close) return;
    //	invalid nestings
    if ((s == "done") || (s == "endwhile") || (s == "endif")) break;
    //	valid nestings
    if (s == "do")  	menu_close ("done", "do");
    if (s == "while") 	menu_close ("endwhile", "while");
    if (s == "if")	menu_close ("endif", "if");
  }
  if (s .size ()) cout << "Found wrong closing key " << s;
  else cout << "Not found closing key " << close ;
  cout << " matching opening key " << open << " (aborting) ..." << endl;
  exit (0);
}
//
//____________________________________________________________________________
static void menu_complex (void * data)
{
  //
  //	Assign complex value
  //
  complex<double> & v = ((complex<double> *) data) [0];
  v = complex_expression (input_next (true));
}
//
//____________________________________________________________________________
static void menu_define (void *)
{
  //
  //	Defines or updates a variable from expression
  //
  string dname (input_next (true));
  define_entry (dname, complex_expression (input_next (true)));
}
//
//____________________________________________________________________________
static void menu_density (void *) 
{
  //
  //	Store next input up to keyword "endgroup" as density operators
  //	definition which will be parsed during dmrg steps
  //
  string s = "density_operator";
  densityaction .push_back (s);
  while ((s = input_next (true)) .size ()) {
    if (s == "endgroup") return;
    densityaction .push_back (s);
  } 
}
//
//____________________________________________________________________________
static void menu_dmrgs (void *)
{
  //
  //	Adds a rule to choose dmrg sates
  //
  Brule	rule;
  rule .br_time = -1.0;
  menu_size (& rule .br_zip);
  menu_size (& rule .br_sites);
  menu_size (& rule .br_min);
  menu_size (& rule .br_max);
  menu_real (& rule .br_cut);
  lftrule .push_back (rule);
  rgtrule .push_back (rule);
  //brule   .push_back (rule);
}
//
//____________________________________________________________________________
static void menu_dmrgstime (void *)
{
  //
  //	Adds a rule to choose dmrg sates during time evolution
  //
  Brule	rule;
  menu_real (& rule .br_time);
  menu_size (& rule .br_zip);
  menu_size (& rule .br_sites);
  menu_size (& rule .br_min);
  menu_size (& rule .br_max);
  menu_real (& rule .br_cut);
  lfttimerule .push_back (rule);
  rgttimerule .push_back (rule);
  //brule   .push_back (rule);
}
//
//____________________________________________________________________________
static void menu_do (void *)
{
  //
  //	Initializes a do loop
  //
  if (loop_nest < -1) loop_nest = -1;
  loop_nest++;
  if ((size_t) loop_nest >= loop_list .size ()) 
    loop_list .resize ((size_t) loop_nest + 1);
  string lname (input_next (true));
  loop_list [loop_nest] .l_index = (double *) define_entry (lname);
  loop_list [loop_nest] .l_start = double_expression (input_next (true));
  loop_list [loop_nest] .l_end   = double_expression (input_next (true));
  loop_list [loop_nest] .l_step  = double_expression (input_next (true));
  loop_list [loop_nest] .l_body	 = input_position ();
  //
  //	Look for closing end of loop
  //
  menu_close ("done", "do");
  loop_list [loop_nest] .l_close = input_position ();
  //
  //	Loop init
  //
  * loop_list [loop_nest] .l_index = loop_list [loop_nest] .l_start;
  //
  //	Check for end of loop and set input position
  //
  if (loop_list [loop_nest] .l_start > loop_list [loop_nest] .l_end)
    input_position (loop_list [loop_nest] .l_close);
  else input_position (loop_list [loop_nest] .l_body);
}
//
//____________________________________________________________________________
static void menu_done (void *)
{
  //
  //	Closes a do loop body
  //
  if (loop_nest < 0) {
    cout << "No start point for end of loop!" << endl;
    exit (0);
  }
  * loop_list [loop_nest] .l_index += loop_list [loop_nest] .l_step;
  if (* loop_list [loop_nest] .l_index > loop_list [loop_nest] .l_end)
    loop_nest--;
  else input_position (loop_list [loop_nest] .l_body);
}
//
//____________________________________________________________________________
static void menu_endgroup (void *)
{
  //
  //	Closing a storage list for next parsing
  //
  cout << "Parsing error: No start for closing endgroup." << endl;
  exit (0);
}
//
//____________________________________________________________________________
static void menu_endwhile (void *)
{
  //
  //	Closes a do while loop body
  //
  if (loop_nest < 0) {
    cout << "No start point for end of while!" << endl;
    exit (0);
  }
  //
  //	Set position at control expression
  //
  input_position (loop_list [loop_nest] .l_body);
  //	evaluate control expression
  bool v = boolean_expression (input_next (true));
  if (v) return;
  //
  //	Go to end of while
  //
  input_position (loop_list [loop_nest] .l_close);
  loop_nest--;
}
//
//____________________________________________________________________________
static void menu_evolution (void *)
{
  //
  //     Activate evolution and define evolution parameters
  //
  evolution_active = true;
  timestring = input_next(true);
  menu_parse(menu_time);
}
//____________________________________________________________________________
static void menu_false (void * data)
{
  //
  //	Set true a given boolean
  //
  bool & v = ((bool *) data) [0];
  v = false;
}
//
//____________________________________________________________________________
static Menu * menu_find (const string & key, Menu * menu)
{
  //
  //	Looks for key in menu (and menu_common) and return menu entry
  //	(or null if not found)
  //
  for ( ; menu ->keyword .size (); menu++) 
    if (menu ->keyword == key) return menu;
  for (menu = menu_common; menu ->keyword .size (); menu++) 
    if (menu ->keyword == key) return menu;
  return 0;
}
//
//____________________________________________________________________________
static void menu_hamb (void*)
{
  //
  //	Store next input up to keyword "ham_end" as hamiltonian definitions
  //	which will be parsed during dmrg steps
  //
  string s;
  ham_definition .clear ();
  while ((s = input_next (true)) .size ()) {
    if (s == "ham_end") return;
    ham_definition .push_back (s);
  } 
}
//
//____________________________________________________________________________
static void menu_hamtb (void*)
{
  //
  //	Store next input up to keyword "timeham_end" as hamiltonian 
  //	definitions which will be parsed during dmrg steps
  //
  string s;
  hamt_definition .clear ();
  while ((s = input_next (true)) .size ()) {
    if (s == "timeham_end") return;
    hamt_definition .push_back (s);
  } 
}
//
//____________________________________________________________________________
static void menu_hame (void *)
{
  //
  //	Calling this routine is a parsing error!
  //
  cout << "Parsing error: keyword ham_end does not close ham_begin." << endl;
  exit (0);
}
//
//____________________________________________________________________________
static void menu_hamte (void *)
{
  //
  //	Calling this routine is a parsing error!
  //
  cout << "Parsing error: keyword timeham_end does not close timeham_begin." 
       << endl;
  exit (0);
}
//
//____________________________________________________________________________

static void menu_if (void *)
{
  //
  //	get next expression as if control variable
  //
  bool v = boolean_expression (input_next (true));
  //
  //	if control variable evaluates to true continue execution
  //
  if (v) return;
  //
  //	Else skip to closing endif
  //
  menu_close ("endif", "if");
}
//
//____________________________________________________________________________
static void menu_integer (void * data)
{
  //
  //	Assigns to *data the (long) value from next expression
  //
  long & v = ((long *) data) [0];
  v = (long) double_expression (input_next (true));
}
//
//____________________________________________________________________________
static void menu_memory (void *)
{
  //
  //	set size of allowed memory (in Mb)
  //
  size_t mem = (size_t) double_expression (input_next (true));
  Storage m;
  m .top (mem);
}
//
//____________________________________________________________________________
static void menu_name (void * data)
{
  //
  //	Set data to the identifier corresponding to next string
  //	The identifier is taken from the list of defined names.
  //
  size_t & v = ((size_t *) data) [0];
  v = name_define (input_next (true));
}
//
//____________________________________________________________________________
void menu_parse ()
{
  //
  //	parses input starting from main menu keys
  //
  keyword = "";
  menu_parse (menu_main);
  if (keyword .size()) {
    cout << "Unrecognized input keyword " << keyword << endl;
    exit (0);
  }
}
//
//____________________________________________________________________________
static void menu_parse (void * data)
{
  //
  //	Parses input keywords according to Menu data entries until a valid
  //	keyword is found
  //
  Menu * menu = (Menu *) data;
  while ((keyword = input_next ()) .size ()) {
    while (keyword .size ()) {
      //
      //	Check for forced end of menu
      //
      if (keyword == "end") {
	keyword = "";
	return;
      }
      //
      //	Looks for keyword in given menu (and menu_common, see above)
      //
      Menu * item = menu_find (keyword, menu);
      //
      //	If not found exit from this menu and leave keyword live
      //	for upper menu (if any)
      //
      if (item == 0) return;
      //
      //	Remove the keyword and process menu entry.
      //	Processing can involve the parsing of input according 
      //	to a different menu and on return keyword can contain 
      // 	again a valid menu entry for this menu.
      //
      keyword = "";
      if (item ->action) item ->action (item ->data);
    }
  }    
}
//
//____________________________________________________________________________
static void menu_pro_begin (void *) 
{
  //
  //	Initialize and define a property
  //
  property = Aproperty ();
  //
  //	Next input must be the property name
  //
  menu_name  (& property .ap_id);
  menu_parse (menu_property);
  //
  //	Add defined property to property list
  //
  size_t insert;
  for (insert = 0; insert < correlation .size (); insert++) {
    if (property .ap_sites    < correlation [insert] .ap_sites)	   break;
    if (property .ap_sites    > correlation [insert] .ap_sites)    continue;
    if (property .ap_id       < correlation [insert] .ap_id)       break;
    if (property .ap_id       > correlation [insert] .ap_id)       continue;
    if (property .ap_ketid    < correlation [insert] .ap_ketid)    break;
    if (property .ap_ketid    > correlation [insert] .ap_ketid)    continue;
    if (property .ap_ketindex < correlation [insert] .ap_ketindex) break;
    if (property .ap_ketindex > correlation [insert] .ap_ketindex) continue;
    if (property .ap_braid    < correlation [insert] .ap_braid)    break;
    if (property .ap_braid    > correlation [insert] .ap_braid)    continue;
    if (property .ap_braindex < correlation [insert] .ap_braindex) break;
    if (property .ap_braindex > correlation [insert] .ap_braindex) continue;
    if (property .ap_index    < correlation [insert] .ap_index)    break;
  }
  correlation .insert (correlation .begin () + insert, property);
}
//
//____________________________________________________________________________
static void menu_pro_ibegin (void *) 
{
  //
  //	Initialize and define an initial Superaction for time evolution
  //
  property = Aproperty ();
  //
  //	Next input must be the property name
  //
  menu_name  (& property .ap_id);
  menu_parse (menu_property);
  //
  //	set start_action property
  //
  start_action = property;
}
//
//____________________________________________________________________________
static void menu_pro_bra (void *)
{
  //
  //	Set bra choice for property
  //
  menu_name (& property .ap_braid);
  menu_size (& property .ap_braindex);
}
//
//____________________________________________________________________________
static void menu_pro_ket (void *)
{
  //
  //	Set bra choice for property
  //
  menu_name (& property .ap_ketid);
  menu_size (& property .ap_ketindex);
}
//
//____________________________________________________________________________
static void menu_pro_mono (void *)
{
  //
  //	Adds a monomial to property polinomial 	
  //
  pro_mono = Amono ();
  pro_mono .am_coeff = 1.0;
  menu_parse (menu_monomial);
  //
  //	Add contribute to polinomial expression
  //
  pro_mono .reorder (property .ap_sites);
  property += pro_mono;
}
//
//____________________________________________________________________________
static void menu_pro_op (void *)
{
  //
  //	Adds a factor to property monomial
  //	
  Afactor af;
  string aname (input_next (true));
  //
  size_t old   = name_action ();
  size_t found = name_action (aname);
  if (found == old) {
    cout << "Input error: operator " << aname << " not defined!" << endl;
    exit (0);
  }
  af .af_op = found;
  aname = input_next (true);
  af .af_dg = 0;
  if (aname == "dagger") {
    af .af_dg = 1;
    aname = input_next (true);
  }
  af .af_st = (long) (complex_expression (aname) .real ());
  pro_mono *= af;
}
//
//____________________________________________________________________________
static void menu_real (void * data)
{
  //
  //	Assigns a real value
  //
  double & r = ((double *) data) [0];
  complex<double> v = complex_expression (input_next (true));
  if (v .imag ()) {
    cout << "Real value expected!" << endl;
    exit (0);
  }
  r = v .real ();
}
//
//____________________________________________________________________________
static void menu_renyi (void *)
{
  //
  //	Adds a double value as argument for computation of Renyi entropy
  //
  if (renyiarg .size () == 0) renyiarg .push_back (1.0);
  complex<double> v = complex_expression (input_next (true));
  if (v .imag ()) {
    cout << "Real value expected!" << endl;
    exit (0);
  }
  if (v. real () != 1.0) renyiarg .push_back (v .real ());
}
//
//____________________________________________________________________________
static void menu_show_bool (void *)
{
  //
  //	Immediate output of next input boolean expression
  //
  if (boolean_expression (input_next (true))) cout << "true ";
  else cout << "false ";
}
//
//____________________________________________________________________________
static void menu_show_endl (void *)
{
  //
  //	Immediate output of newline
  //
  cout << endl;
}
//
//____________________________________________________________________________
static void menu_show_str (void *)
{
  //
  //	Immediate output of next input string
  //
  cout << input_next (true) << " ";
}
//
//____________________________________________________________________________
static void menu_show_val (void *)
{
  //
  //	Immediate output of next input value (complex)
  //
  complex<double> v = complex_expression (input_next (true));
  if (v .imag())  cout << v << " ";
  else		  cout << v .real () << " ";
}
//
//____________________________________________________________________________
static void menu_size (void * data)
{
  //
  //	Assigns to *data the (size_t) value from next expression
  //
  size_t & v = ((size_t *) data) [0];
  v = (size_t) double_expression (input_next (true));
}
//
//____________________________________________________________________________
static void menu_space (void *)
{
  //
  //	Processes space keyword entry defining a subspace of a 
  //	single site block
  //
  if (site_actions) {
    cout << "space definition must precede operator definition!" << endl;
    exit (0);
  }
  space_reset ();	// reset space variables
  //
  //	Get the identifier
  //
  string sname (input_next (true));
  qid = name_define (sname);
  //
  //	Use menu_subspace to parse space definitions
  //
  menu_parse (menu_subspace);
  //
  //	Assume bose statistic if no specification was given
  //
  if (qstatistic == 0) qstatistic = 1;
  //
  //	Consistency checks with previous definitions (if any)
  //
  if ((qvalues != name_quantum ()) || (qparities != name_parity ())) {
    cout << "Incompatible quantum specifications! Subspace " << sname << endl;
    exit (0);
  }
  //
  //	Add subspace to one-site lattice block
  //
  blocklft [1] .add_states (Qnumber (qid, qvalue, qparity, qlftlnk, qrgtlnk), 
			    qstates, qstatistic);
  blockrgt [1] .add_states (Qnumber (qid, qvalue, qparity, qlftlnk, qrgtlnk), 
			    qstates, qstatistic);
  site_states = blocklft [1] .states ();
}
//
//____________________________________________________________________________
static void menu_space_a (void *)
{
  //
  //	Adds a quantum name and value to subspace properties
  //
  string qname (input_next (true));
  size_t qindex = name_quantum (qname);
  if (qvalue .size () <= qindex) qvalue .resize (qindex + 1);
  if (qalternate .size () <= qindex)  qalternate .resize (qindex + 1);
  menu_real (& qvalue [qindex]);
  qalternate [qindex] = true;
  qvalues++;
}
//
//____________________________________________________________________________
static void menu_space_b (void *)
{
  //
  //	Set bose statistic for actual subspace
  //	
  qstatistic = 1;
}
//
//____________________________________________________________________________
static void menu_space_f (void *)
{
  //
  //	Set fermi statistic for actual subspace
  //	
  qstatistic = -1;
}
//
//____________________________________________________________________________
static void menu_space_p (void *)
{
  //
  //	Adds a parity name and value to subspace properties
  //
  string pname (input_next (true));
  size_t pindex = name_parity (pname);
  if (qparity .size () <= pindex) qparity .resize (pindex + 1);
  menu_complex (& qparity [pindex]);
  qparities++;
}
//
//____________________________________________________________________________
static void menu_space_q (void *)
{
  //
  //	Adds a quantum name and value to subspace properties
  //
  string qname (input_next (true));
  size_t qindex = name_quantum (qname);
  if (qvalue .size () <= qindex) qvalue .resize (qindex + 1);
  if (qalternate .size () <= qindex)  qalternate .resize (qindex + 1);
  menu_real (& qvalue [qindex]);
  qalternate [qindex] = false;
  qvalues++;
}
//
//____________________________________________________________________________
static void menu_symmetry (void *)
{
  //
  //	Defines a one-site lattice symmetry
  //
  menu_action (0);
  //
  // 	Update symmetries list
  //
  size_t sym = site_actions - 1;
  //
  size_t found;
  for (found = 0; found < symmetry .size (); found++)
    if (symmetry [found] == sym) break;
  if (found == symmetry .size ()) symmetry .push_back (sym);
}
//
//____________________________________________________________________________
static void menu_target (void *)
{
  //
  //	Processes target keyword defining a subspace of target states
  //
  space_reset ();	// resets space definition parameters
  //
  //	Get target identifyer
  //
  string tname (input_next (true));
  qid = name_define (tname);
  //
  //	Use menu_space to parse target definitions
  //
  menu_parse (menu_subspace);
  //
  //	Assume bose statistic if no specification was given
  //
  if (qstatistic == 0) qstatistic = 1;
  //
  //	Consistency check with previous definitions (if any)
  //
  if ((qvalues != name_quantum ()) || (qparities != name_parity ())) {
    cout << "Incompatible quantum specifications! Target " << tname << endl;
    exit (0);
  }
  //
  //	Look if a target for given chain len was yet defined
  //
  size_t found;
  for (found = 0; found < target_request .size (); found++)
    if (target_request [found] .sites () == qsites) break;
  if (found == target_request .size ()) 
    target_request .push_back (Qtarget (Qspace (), qsites));
  target_request [found] .add_states (Qnumber (qid, qvalue, qparity, 
					       qlftlnk, qrgtlnk),
				      qstates, qstatistic);
}//
//____________________________________________________________________________
static void menu_true (void * data)
{
  //
  //	Set true a given boolean
  //
  bool & v = ((bool *) data) [0];
  v = true;
}
//
//____________________________________________________________________________
static void menu_ttarget (void *)
{
  //
  //	Set ket choice for initial state of time evolution
  //
  menu_name (& ttarget_id);
  menu_size (& ttarget_index);
}
//
//____________________________________________________________________________
static void menu_while (void *)
{
  //
  //	Initializes a while loop (using a new Loop control)
  //
  if (loop_nest < -1) loop_nest = -1;
  loop_nest++;
  if ((size_t) loop_nest >= loop_list .size ()) 
    loop_list .resize ((size_t) loop_nest + 1);
  loop_list [loop_nest] .l_index = 0;
  loop_list [loop_nest] .l_start = 0;
  loop_list [loop_nest] .l_end   = 0;
  loop_list [loop_nest] .l_step  = 0;
  //
  //	We mark start and end of while body 
  //	(first item is the control expression)
  //
  loop_list [loop_nest] .l_body	 = input_position ();
  //
  //	Look for closing end of while
  //
  menu_close ("endwhile", "while");
  loop_list [loop_nest] .l_close = input_position ();
  //
  //	We are positioned at the end of while so we can use 
  //	the endwhile procedure.
  //
  menu_endwhile (0);
}
//
//____________________________________________________________________________
void space_reset ()
{
  //
  //	Reset local variables defining a subspace of state 
  //	These variable are used to define a subspace of one-site lattice
  //	block or a superblock target.
  //
  qstates = 0;
  qstatistic = 0;
  qvalues = 0;
  qparities = 0;
  qid = 0;
  qsites = 0;
  for (size_t i = 0; i < qvalue  .size (); ++i) qvalue  [i] = 0.0;
  for (size_t i = 0; i < qparity .size (); ++i) qparity [i] = 0.0;
}
//
//============================================================================
