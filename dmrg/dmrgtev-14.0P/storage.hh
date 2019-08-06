//	storage.hh			(C) 2009 Fabio Ortolani fbo 090716
//	==================================================================
#ifndef STORAGE_HH
#define STORAGE_HH
#include <vector>
#include <string>
using namespace std;
//
//============================================================================
class Storage_chunk
{
  //
  //	Storage_chunck class stores informations about a chunk of memory 
  //	used by the program.
  //	This class is accessed only by friend Storage class (see below).
  //	Every chunk resides in the heap if it is active, otherway it can be 
  //	stored in a file and retrieved in the heap when needed. 
  //	A Storage_chunk can be in common with many objects (each object must
  //	refer to the chunk via a Storage element).
  //	For every chunk of memory we store the size (bytes), how many 
  //	objects are using the chunk, its status (if it is active, constant, 
  //	i.e. no more modified) and its memory or file position.
  //	A chunk can be marked as removable from memory only when it is no 
  //	more modified by any object. If a chunk is removable and not active 
  //	for computation it can be stored in a swap file to free memory 
  //	and its position is recorded.
  //	A chunk is deallocated only when no active object use it. 
  //
private:
  //
  Storage_chunk * sc_prev;	// all chunks are linked together in a 
  Storage_chunk * sc_next;	// bidiretional list.
  size_t 	  sc_size;	// size of memory chunk (bytes)
  long   	  sc_links;	// number of references to this chunk
  long		  sc_status;	// status flag
  void * 	  sc_storage;  	// memory pointer (if not null)
  long	 	  sc_file;     	// swap file number
  long	 	  sc_position; 	// swap file position
  //
  inline Storage_chunk ();	// default ctor 		  [storage.hh]
  //
  friend class Storage;		
  //
};
//
//	Storage_chunk	status flags
//
#define SC_UNUSED	0x0001
#define SC_USED		0xfffe
#define SC_SWAPPED	0x0002
//
//____________________________________________________________________________
class Storage_swap
{
  //
  //	Storage_swap class describes a file used to swap memory chunks.
  //	It is accessed only by Storage class (see below).
  //
private:
  //
  string  ss_name;		// file name
  int     ss_descriptor;	// file descriptor (when opened)
  long    ss_position;		// free file position (size of file in bytes)
  //
public:
  //
  inline Storage_swap ();	// default ctor 		  [storage.hh]
  //
private:
  //
  friend class Storage;	
  //
};
//
//____________________________________________________________________________
class Storage
{
  //
  //	Storage class is a controller for administration of memory chunks.
  //	Every chunk of memory is accessed and controlled by one or more 
  //	instances of Storage class which is responsible of allocation and
  //	deallocation or swapping of memory chunks.
  //	Every Storage class accesses only one chunk of memory (or none).
  //
private:
  //
  //	Only one true slot, a private pointer to a Storage_chunk
  //
  Storage_chunk * s_chunk;	// Only one true slot
  //
public:
  //
  inline Storage  (size_t = 0);		// default ctor           [storage.hh]
  //
  //	Copy constructor does not allocate a new chunk of memory.
  //	The chunk is grabbed and its reference counter is incremented.
  //
  inline Storage  (const Storage &);	// copy constructor       [storage.hh]
  inline ~Storage ();			// destructor
  //
  //	No direct assignment is allowed. The only assignment is the
  //	grabbing of memory. To be sure of this we define a grab 
  //	operator << and redefine the assignment as a private operator 
  //	equivalent to grabbing.
  //
  //	The grab operator removes actual memory (if any) from destination
  //	and grabs the requested chunk of memory.
  //
  //	A true copy of memory must be done explicitly by methods using 
  //	Storage classes.
  //
  Storage & operator << (const Storage &);	               // [storage.hh]
  //
  inline size_t size     () const; // actual memory size (bytes)  [storage.hh]
  inline void 	release  () const; // mark chunk as removable	  [storage.hh]
  void *        storage  () const; // actual memory (marked used) [storage.cc]
  void          storage  (size_t); // (de)allocation of memory	  [storage.cc]
  //
  void   swap_garbage () const;				       // [storage.cc] 
  void   swap_unlink  () const;				       // [storage.cc]
  void   top          (size_t) const;          		       // [storage.cc]
  string usage        () const;				       // [storage.cc]
  size_t released     () const; // get released memory amount	  [storage.cc] 
  //
private:
  //
  void * allocate (size_t) const;			       // [storage.cc]
  void	 create   (size_t);			       	       // [storage.cc]
  void   free	  (size_t = 0) const;			       // [storage.cc]
  void   restore  (Storage_chunk *) const;		       // [storage.cc]
  void   save     (Storage_chunk *) const;		       // [storage.cc]
  void   unlink   ();					       // [storage.cc]
  //
  inline Storage & operator  = (const Storage &);	       // [storage.hh]
  //
};
//
//============================================================================
inline Storage_chunk::Storage_chunk ()
  : sc_prev	(0),
    sc_next	(0),
    sc_size	(0),
    sc_links	(0),
    sc_status	(0),
    sc_storage	(0),
    sc_file	(-1),
    sc_position	(0)
{ 
  //
  //	Default ctor.
  //
}
//
//============================================================================
inline Storage_swap::Storage_swap ()
  : ss_name		(""),
    ss_descriptor	(0),
    ss_position		(0)
{ 
  //
  //	Copy ctor.
  //
}
//
//============================================================================
inline Storage::Storage (size_t size)
{
  //
  //	Allocate a new chunk of memory (if size > 0)
  //
  create (size);
}
//
//____________________________________________________________________________
inline Storage::Storage  (const Storage & other)
  : s_chunk (other.s_chunk)
{ 
  //
  //	Grabs other chunk
  //
  if (s_chunk) s_chunk -> sc_links++; 
}
//
//____________________________________________________________________________
inline Storage::~Storage ()
{
  //
  //	Removes connection to actual chunk (removing the chunk if this is
  //	last connection).
  //
  unlink ();
}
//
//____________________________________________________________________________
inline Storage & Storage::operator = (const Storage & other)
{
  //
  //	Define assignment as grabbing.
  //
  *this << other;
  return *this;
}
//
//____________________________________________________________________________
inline void Storage::release () const
{
  //
  //	Mark the chunk as removable 
  //
  if (s_chunk		 &&
      s_chunk -> sc_size &&
      s_chunk -> sc_storage)  s_chunk -> sc_status |= SC_UNUSED;
}
//
//____________________________________________________________________________
inline size_t Storage::size () const
{
  //
  //	Returns size of memory chunk (bytes).
  //
  if (s_chunk) return s_chunk -> sc_size;
  return 0;
}
//
//============================================================================
#endif	// STORAGE_HH
