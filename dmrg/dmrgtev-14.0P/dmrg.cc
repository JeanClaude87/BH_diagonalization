#include "timing.hh"
#include "algebra.hh"
#include "block.hh"
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
using namespace std;
//
//============================================================================
void dmrg_showinfo        (const Block &, const Block &);      	  // [dmrg.cc]
void hamiltonian_energies (Block &, Block &, Qtarget &);    // [superblock.cc]
void properties		  (Block &, Block &);	   	    // [superblock.cc]
void evolution            (Block &, Block &);		    // [superblock.cc]	
bool quantum_alternate	  ();					  // [menu.cc]
void updatebase		  (Block &, Block &, long);	    // [superblock.cc]
//
//============================================================================
size_t			length			= 8;
size_t			sites_step		= 2;
long			sites_diff		= 0;
complex<double>		periodic		= 1.0;
bool			reflect_blocks		= false;
bool			reflect_universe	= false;
size_t			reflection		= 0;
double			truncation		= 0.0;
//
size_t	   	       	show_blocks		= 0;
size_t			show_density		= 0;
size_t 			show_hamiltonian	= 0;
size_t   	       	show_interaction	= 0;
bool			show_point		= false;
bool			show_memory		= false;
size_t			show_system		= 0;
bool			show_target		= false;
bool			show_states		= true;
bool			show_selection		= true;
//
bool                    evolution_active        = false;
//
size_t			sites			= 0;
size_t			old_sites		= 0;
size_t			sites_exact		= 1;
size_t			sites_lft		= 1;
size_t			exact_lft		= 1;
size_t			sites_rgt		= 1;
size_t			exact_rgt		= 1;
size_t			old_sys			= 0;
size_t			old_uni			= 0;
size_t			old_lft			= 0;
size_t			old_rgt			= 0;
//
size_t			zip             = 0;
size_t			zips			= 0;
size_t			zipinner		= 0;
size_t			zipmax			= 0;
long			zipdir			= 0;
size_t			zipodd			= 0;
bool			halfsweep		= false;
//
size_t			min_lft			= 0;
size_t			max_lft			= 0;
size_t			min_rgt			= 0;
size_t			max_rgt			= 0;
double          n_cutlft        = 0;
double          n_cutrgt        = 0;
//
size_t			idop_none		= 0;
size_t			idop_base		= 0;
size_t			idop_hamiltonian       	= 0;
size_t			id_none			= 0;
size_t			id_base			= 0;
size_t			id_hamiltonian		= 0;
//
Barray			blocklft;
Barray			blockrgt;
//vector<Brule>		brule;
vector<Brule>		lftrule;
vector<Brule>		rgtrule;
vector<Qtarget>		target_request;		
Qtarget			target;
vector<size_t>		symmetry;
//
extern Aarray		super_state;		       	   // [superaction.cc]
Storage 		report;
//
//============================================================================
void dmrg_initialize ()
{
  //
  //	Initialization (before input processing)
  //
  //
  //	Initialize global arrays
  //
  blocklft .clear ();
  blockrgt .clear ();
  target_request .clear ();
  //
  //	Define common names
  //
  idop_none	   		= name_action ("No name");
  idop_base 	   	= name_action ("Block base");
  idop_hamiltonian	= name_action ("Hamilton H");
  //
  id_none 	   		= name_define ("No name");
  id_base	   		= name_define ("Block base");
  id_hamiltonian	= name_define ("Hamilton H");
}
//
//____________________________________________________________________________
void dmrg_process ()
{
  //
  //	Main DMRG control routine
  //
  if (blocklft .size () < 2) {
    cout << "No single site block defined!" << endl;
    exit (0);
  }
  if (show_point && (show_blocks == 0))  blocklft [1] .show ("Point block");
  //
  //	DMRG loop
  //
  Block system, universe;
  bool  updatelft = false;
  bool  updatergt = false;
  long  sites_diff = ((long) sites_lft) - ((long) sites_rgt);
  zipdir = 1;  // Start with right sweep for full sweep else left sweep
  if (halfsweep) zipdir = -1;
  bool middle = true;
  while (true) {
    //
    //	zip count is incremented in the middle of the chain
    //
    if (middle) zip++;
    //
    //	Actual sites
    //
    sites = sites_lft + sites_rgt + 2;
    //
    Block & lft	  = blocklft [sites_lft];
    Block & rgt   = blockrgt [sites_rgt];
    //
    system   = Block (lft, blocklft [1]);
    universe = Block (rgt, blockrgt [1]);
    //
    //	Reflection in universe (left blocks are never reflected)
    //
    if (reflect_blocks) {
        if (quantum_alternate ()) {
	cout << "reflect_blocks not supported with alternate quantum numbers"
	    << endl;
	exit (0);
      }
      universe .reflectlft ();
      universe .reflectrgt ();
      reflection |= 2;
    }
    if (reflect_universe) {
      universe .reflect ();
      reflection |= 1;
    }
    //
    //	Show iteration line info
    //
    dmrg_showinfo (system, universe);
    //
    //	Find possible targets
    //
    double qa = 1.0;
    if (system .sites () % 2) qa = -1.0;
    if ((reflection & 1) && (universe .sites () % 2 == 0)) qa = -qa;
    target = Qtarget (system. quantum (), universe .quantum (), sites, qa);
    //
    //	Select targets for actual chain length
    //
    size_t sites_max = 0;	
    //
    // 	special attention for the presence of a general target which has the 
    // 	precedence when sites == length
    //
    size_t wanted = target_request .size ();	
    for (size_t index = 0; index < target_request .size (); index++) {
      size_t s = target_request [index] .sites ();
      if (s == 0 && sites >= length) {
	wanted = index;
	break;
      }
      //
      //   get the rule with highest sites less than actual sites
      //
      if (sites >= s && s >= sites_max) {
	sites_max = s;
	wanted = index;
      }
    }
    if (wanted == target_request .size ()) {
      cout << "No target found (add a valid target specification)!" << endl;
      exit (0);
    }
    //
    //	Select possible targets near wanted targets
    //
    quantum_select (target, target_request [wanted]);
    if (target .states () == 0) {
      cout << "No valid target found! Check your input target specifications."
	   << endl;
      exit (0);
    }
    //
    //	Compute hamiltonian and target states (with relative energies)
    //
    hamiltonian_energies (system, universe, target);
    //
    //	Record actual lattice decomposition
    //
    old_sys 	= system   .sites ();
    old_uni 	= universe .sites (); 
    old_sites 	= sites;
    old_lft	= sites_lft;
    old_rgt 	= sites_rgt;
    //
    //	Find next step site numbers and check which list is to update.
    //  During warmup we sweep if actual sites >= zipinner
    //	With final size use finite size algorithm. 
    //  
    zipmax = 0;
    if ((sites < length) && zipinner && (zipinner <= sites)) zipmax = zips;
    if (sites >= length) zipmax = zips;
    if ((sites_lft <= exact_lft) && (sites_rgt <= exact_rgt)) zipmax = 0;
    //

    //
    if (zip < zipmax) {
    
      //
      //	Finite size algorithm
      //
      
// ////////////////////////////////////////////////////////////////////////////////////////////////// To have properties at each sweep
      
  if ((sites_lft==sites_rgt) &&  (zip==5)) { 
      		 properties (system, universe);
            cout <<endl;
      }

      // //////////////////////////////////////////////////////////////////////////////////////////////////
      
      sites_lft += zipdir;
      sites_rgt -= zipdir;
      if (zipdir > 0) {
	updatelft = true;
	updatergt = false;
      }
      else {
	updatelft = false;
	updatergt = true;
      }
      if (sites_lft <= exact_lft) zipdir =  1;
      if (sites_rgt <= exact_rgt) zipdir = -1;
      middle = (sites_diff == (((long) sites_lft) - ((long) sites_rgt)));
      if (halfsweep && middle) {	zipdir = -1;}

    }
    else {
      //
      //	system and universe have the same length or near.
      //	It is a good place for property evaluation.
      //	

      properties (system, universe);
      //
      //	End of DMRG algorithm
      //
      if (sites >= length) break;
      //
      //	Infinite size algorithm (lattice increment)
      //
      if ((length - sites) == 1) sites_step = 1;
      if (sites_step > 1) {
	sites_lft++;
	sites_rgt++;
	updatelft = true;
	updatergt = true;
      }
      else { 
	//
	// 	sites_step = 1 
	//
	//	Only one of system or universe is growing but we
	//	want to update both block lists in next step, so
	//	we upgrade only the old number of sites of growing 
	//	part.
	//	
	if (sites_lft <= sites_rgt) {
	  sites_lft++;
	  updatelft = true;
	  updatergt = false;
	}
	else {
	  sites_rgt++;
	  updatelft = false;
	  updatergt = true;
	}
      }
      sites_diff = ((long) sites_lft) - ((long) sites_rgt); 
      middle = true;
      zip = 0;
      zipdir = 1;
      if (halfsweep) {
	zipdir = -1;
	updatergt = false;
      }
    }
    //
    //	Set the number of wanted dmrg states according to input rules.
    //	A rule specify the minimum and maximum number of DMRG states 
    //	starting from a given sweep (zip) and number of sites.
    //
    size_t actual_zip = zip;
    if (zip == 0) actual_zip = 1;
    min_lft = max_lft = min_rgt = max_rgt = 0;
    //
    for (size_t rule = 0; rule < lftrule .size (); rule++) {
      size_t z = lftrule [rule] .br_zip;
      size_t s = lftrule [rule] .br_sites;
      if ((actual_zip >= z) && (system .sites () >= s)) {
          if (min_lft < lftrule [rule] .br_min) min_lft = lftrule [rule] .br_min;
          if (max_lft < lftrule [rule] .br_max) max_lft = lftrule [rule] .br_max;
          n_cutlft = lftrule [rule] .br_cut;
      }
    }
    for (size_t rule = 0; rule < rgtrule .size (); rule++) {
      size_t z = rgtrule [rule] .br_zip;
      size_t s = rgtrule [rule] .br_sites;
      if ((actual_zip >= z) && (universe .sites () >= s)) {
          if (min_rgt < rgtrule [rule] .br_min) min_rgt = rgtrule [rule] .br_min;
          if (max_rgt < rgtrule [rule] .br_max) max_rgt = rgtrule [rule] .br_max;
          n_cutrgt = rgtrule [rule] .br_cut;
      }
    }
    //
    //	Update the number of sites without states contraction.
    //	The corresponding Hilbert space is complete and no basis
    //	transformation will be done.	
    //
    if ((system .states () <= max_lft) && 
	(exact_lft == system .sites() - 1)) {
      exact_lft = system .sites ();
    }
    if ((universe .states () <= max_rgt) && 
	(exact_rgt == universe .sites () -1)) {
      exact_rgt = universe .sites ();
    }
    //
    //	Check if we need changement of basis
    //
    bool updatesys = (system   .sites () > exact_lft) && updatelft; 
    bool updateuni = (universe .sites () > exact_rgt) && updatergt;
    long updatetag = 0;
    if (updatesys) updatetag += 1;
    if (updateuni) updatetag -= 1;
    size_t extra_lft = 2 * max_lft - min_lft;
    if (extra_lft >= system .states ()) 
      extra_lft = (max_lft + system .states ())/2;
    size_t extra_rgt = 2 * max_rgt - min_rgt;
    if (extra_rgt >= universe .states ()) 
      extra_rgt = (max_lft + universe .states ())/2;
    if (updatesys || updateuni)	updatebase (system, universe, updatetag);
    //
    //	Upgrade min and max numbers of states with these block sites 
    //	for next sweep or next chain length.
    //
    

    size_t rlft = max_lft - min_lft;
    size_t rrgt = max_rgt - min_rgt;
    
    if ( n_cutlft == 0 ){    
    	if (updatesys && 
		(min_lft < system   .states ()) &&
		(system .states () < max_lft)) {
      	min_lft = system   .states ();
      	max_lft = min_lft + rlft;
      	if (max_lft > extra_lft) max_lft = extra_lft;
      	  Brule rule;
      	  rule .br_zip   = actual_zip;
      	  rule .br_sites = system .sites ();
      	  rule .br_min   = min_lft;
      	  rule .br_max   = max_lft;
     	  rule .br_cut   = n_cutlft;
      	  lftrule .push_back (rule);
    	}
    }	

	if ( n_cutrgt == 0 ){   
    	if (updateuni && 
		(min_rgt < universe .states ()) &&
		(universe .states () < max_rgt)) {
      	min_rgt = universe .states ();
      	max_rgt = min_rgt + rrgt;
      	if (max_rgt > extra_rgt) max_rgt = extra_rgt;
      	  Brule rule;
      	  rule .br_zip   = actual_zip;
      	  rule .br_sites = universe .sites ();
      	  rule .br_min   = min_rgt;
      	  rule .br_max   = max_rgt;
      	  rule .br_cut   = n_cutrgt;
      	  rgtrule .push_back (rule);
    	}
	}
    //
    //	Some output if requested
    //
    if (system   .sites () <= show_system) system   .show ("System   block");
    if (universe .sites () <= show_system) universe .show ("Universe block");
    //
    //	Update block lists
    //
    if (updatelft) {
      if (system .sites () < blocklft .size ()) 
	blocklft [system .sites ()] .actionclear (0, name_action ());
      blocklft [system .sites ()] = system;
      cout << "Updated blocklft [" << system .sites () << "] " 
	   << system .states () << " states. Update = " << updatesys << endl;
    }
    if (halfsweep && middle) {
      for (size_t n = 2; n <= sites_rgt; ++n) { 
	old_sites = 0; 	// to avoid projection
	if (n < blockrgt .size ()) 
	  blockrgt [n] .actionclear (0, name_action ()); 
	blockrgt [n] << blocklft [n];
      }
      cout << "Copyed  blocklft [2..." << sites_rgt << "] to blockrgt." 
	   << endl; 
    }
    else if (updatergt) {
      if (universe .sites () < blockrgt .size ())
	blockrgt [universe .sites ()] .actionclear (0, name_action ());
      blockrgt [universe .sites ()] = universe;
      blockrgt [universe .sites ()] .reflectreset ();
      cout << "Updated blockrgt [" << universe .sites () << "] " 
	   << universe .states () << " states. Update = " << updateuni << endl;
    }
    for (size_t n = 2; n < blocklft .size (); ++n) blocklft [n] .release ();
    for (size_t n = 2; n < blockrgt .size (); ++n) blockrgt [n] .release ();
    } // while (true)
  //
  if (evolution_active) evolution(system, universe);
}
//
//____________________________________________________________________________
void	dmrg_showinfo (const Block & system, const Block & universe)
{
  //
  //	Prints out a line with dmrg iteration info
  //
  stringstream out;
  size_t sl, sr, ul, ur;
  sl = system .siteslft ();
  sr = system .sitesrgt ();
  out << "> Sweep " << zip << ": " << sl << "+" << sr;
  ul = universe .siteslft ();
  ur = universe .sitesrgt ();
  if (ul) out << "+" << ul;
  if (ur) out << "+" << ur;
  out << "=" << sites << " ";
  size_t hs = blocklft [sl] .states () * blocklft [sr] .states ();
  out << blocklft [sl] .states () << "x" << blocklft [sr] .states ();
  if (ul) {
    out << "x" << blockrgt [ul] .states ();
    hs *= blockrgt [ul] .states ();
  }
  if (ur) {
    out << "x" << blockrgt [ur] .states ();
    hs *= blockrgt [ur] .states ();
  }
  out << "=" << hs << " states ";
  if (truncation > 0.0) 
    out << "(" << scientific << setprecision (2) << truncation << ") ";
  out << "<";
  size_t outl, outr;
  string sout = out .str ();
  outl = (69 - sout .size ())/2;
  outr =  69 - sout .size () - outl;
  cout << setw (outl) << setfill ('=') << "" << sout 
       << setw (outr) << setfill ('=') << "" << setw (0) << setfill (' ');
  cout << " " << timenow () << endl;
  out .str ("");
  if (show_memory) 
      cout << "Actual memory usage: " << report .usage ()  << endl; 
  //
  //
  //	Some output if requested
  //
  if (sl == 1 && sl <= show_blocks) blocklft [1]  .show ("Point block"); 
  if (sl > 1  && sl <= show_blocks) blocklft [sl] .show ("Left  block");
  if (ul > 1  && ul <= show_blocks) blockrgt [ul] .show ("Right block");
  if (ur > 1  && ur <= show_blocks) blockrgt [ur] .show ("Right block");
}
//
//============================================================================
