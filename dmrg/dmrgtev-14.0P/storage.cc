//	storage.cc			(C) 2009 Fabio Ortolani fbo 090716
//	==================================================================
#include "storage.hh"
#include <cstdlib>
#include <fcntl.h>
#include <cstring>
#include <unistd.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
//
//============================================================================
static Storage_chunk *	storage_first	= 0;	// head of chunk list
static Storage_chunk *	storage_last 	= 0;	// tail of chunk list
//
//	Maximum allocated storage:	256 Mb (can be changed from input)
//
static const size_t 	megab 		= 1024L * 1024L;       	// 1 Mb
static size_t		memory_top     	= 256L * megab;
//
//	Memory statistic
//
static size_t memory_use 		= 0;
static size_t memory_max 		= 0;
//
//	Maximum swap file size:	about 16 Gb
//	Maximum number of swap files: 25 (Total swap space: about 400 Gb)
//	The basename of swap files can be changed from input, size and 
//	number of swap files are hardcoded here. 
//
static const long swap_top		= 16L * 1024L * 1024L * 1024L;
static const long swap_files 		= 26;
static string 	  swap_name		= "dmrgcc"; // basename for swap files
//
//	Actual swap file (index) in use
//
static size_t 		    swap_actual	= 1;
static vector<Storage_swap> swap_file	(swap_files);
//
//	Swap file statistic
//
static double swap_total		= 0.0;
static double swap_max			= 0.0;
static double swap_removed		= 0.0;
static double swap_active		= 0.0;
//
//============================================================================
void * Storage::allocate (size_t size) const
{
  //
  //	Allocates size bytes of memory checking memory_top limit
  //
  if (memory_use + size > memory_top) {
    //
    //	Reached top of allowed memory. Try to free memory swapping  
    //	removable chunks to file.
    //
    free (size);
    //
    //	If unable to free memory print an error message and abort
    //
    if (memory_use + size > memory_top) {
      cout << "No room for storage! Used "
	   << usage () << endl
	   << "(requested " << ((double) memory_use /megab) 
	   << "+" << ((double) size /megab) 
	   << " Mb, actual top limit " << ((double) memory_top/megab) 
	   << " Mb)" << endl;
      exit (0);
    }
  }
  //
  //	Actual memory allocation (initialized to zero).
  //
  char * memory = new char [size];
  memset (memory, 0, size);
  //
  //	Update memory statistics
  //
  memory_use += size;
  if (memory_max < memory_use)	memory_max = memory_use;
  return (void *) memory;
}
//
//____________________________________________________________________________
void Storage::create (size_t size)
{
  //
  //	Allocates a new chunk of memory of size bytes if size > 0.
  //	This routine is called only when s_chunk is undefined.
  //
  if (size > 0) {
    s_chunk = new Storage_chunk;
    //
    //	Allocate memory for the chunk and initialize references count
    //
    s_chunk ->sc_storage = allocate (size);
    s_chunk ->sc_size = size;
    s_chunk ->sc_links = 1;
    //
    //	Connect the chunk to the list of chunks.
    //	New allocated chunks are positioned at the start of the list
    //	(they are moved to the end when swapped to a file).
    //
    s_chunk ->sc_prev = 0;
    s_chunk ->sc_next = storage_first;
    if (storage_first) 	storage_first -> sc_prev = s_chunk;
    else		storage_last = s_chunk;
    storage_first = s_chunk;
  }  
  else s_chunk = 0;
}
//
//____________________________________________________________________________
void Storage::free (size_t size) const
{
  //
  //	Looks for a removable storage chunk greater of given size
  //	and save it to swap. 
  //	If size == 0 swap all removable chunks
  //
  Storage_chunk * entry;
  Storage_chunk * next = 0;
  for (entry = storage_first; entry; entry = next) {
    next = entry -> sc_next;
    if (entry ->sc_storage && 
	(entry ->sc_status & SC_UNUSED) &&
	entry ->sc_size >= size) {
      save (entry);
      //
      //	Deallocate memory and update memory statistic
      //
      memory_use -= entry -> sc_size;
      delete [] (char *) entry -> sc_storage;
      entry -> sc_storage = 0;
      if (size) return;
     }
  }
  //
  //	if no chunk of given size was found try to swap all removable chunks
  //
  if (size) free (0);
}
//
//____________________________________________________________________________
Storage & Storage::operator << (const Storage & other)
{
  //
  //	Removes actual chunk of left operand and grab right memory chunk
  //
  if (this != & other) {
    unlink ();
    s_chunk = other .s_chunk;
    if (s_chunk) s_chunk -> sc_links++;
  }
  return *this;
}
//
//____________________________________________________________________________
size_t Storage::released () const
{	
  //
  //	Returns amount of released storage
  //
  size_t rel = 0;
  Storage_chunk * entry;
  Storage_chunk * next = 0;
  for (entry = storage_first; entry; entry = next) {
    next = entry -> sc_next;
    if (entry ->sc_storage && 
	entry ->sc_status & SC_UNUSED ) 
      rel += entry ->sc_size;
  }
  return rel;
}
//
//____________________________________________________________________________
void Storage::restore (Storage_chunk * chunk) const
{
  //
  //	Restore a chunk of memory from swap file
  //
  long size = chunk -> sc_size;
  int 	fd = swap_file [chunk -> sc_file] .ss_descriptor;
  long  position = chunk -> sc_position;
  chunk -> sc_storage = allocate (chunk -> sc_size);
  lseek (fd, position, SEEK_SET);
  if (read (fd, chunk -> sc_storage, size) != size) {
    cout << "Error reading from swap file " 
	 << swap_file [chunk -> sc_file] .ss_name << endl;
    exit (0);
  }
  //
  // mark chunk as used (remove unused flag)
  //
  chunk -> sc_status &= SC_USED;	
}
//
//____________________________________________________________________________
void Storage::save (Storage_chunk * chunk) const
{
  //
  //	Save a chunk of memory on a file 
  //	(if there is too much garbage in the swap files clean them)
  //
  if ((swap_removed > 0.5 * swap_total) ||
      (swap_removed > (1024.0 * megab))) swap_garbage ();
  long  size = chunk -> sc_size;
  if (chunk ->sc_file <= 0) {
    if (size > swap_top) {
      cout << "Chunk size too large for swap!" << endl;
      exit (0);
    }
    //
    //	Look for a swap file and position
    //
    if (swap_file [swap_actual] .ss_position + size > swap_top) {
      swap_actual++;
      if (swap_actual >= swap_file .size ()) {
	cout << "Reached swap size limit! " << endl;
	exit (0);
      }
    }
    //
    //	Record file and position
    //
    chunk -> sc_file = swap_actual;
    chunk -> sc_position = swap_file [swap_actual] .ss_position;
    swap_file [swap_actual] .ss_position += size;
    swap_total += size;
    //
    //	Update swap statistics
    //
    if (swap_max < swap_total) swap_max = swap_total;
    if (swap_active < swap_total - swap_removed) 
      swap_active = swap_total - swap_removed;
    //
    //	Move the chunk to the end of the list of so that swapped 
    //	chunks are ordered and it is easy to clean swap files from
    //	garbage (removed swapped chunks)
    //
    if (chunk ->sc_prev) chunk ->sc_prev ->sc_next = chunk ->sc_next;
    else		       storage_first             = chunk ->sc_next;
    if (chunk ->sc_next) chunk ->sc_next ->sc_prev = chunk ->sc_prev;
    else		       storage_last              = chunk ->sc_prev;
    chunk ->sc_next = 0;
    chunk ->sc_prev = storage_last;
    if (storage_last) storage_last ->sc_next = chunk;
    else	      storage_first = chunk;
    storage_last = chunk;
    //
    //	Check if the file is opened
    //
    if (swap_file [swap_actual] .ss_descriptor <= 0) {
      //
      //   Set file name and open
      //
      stringstream name;
      name << swap_name << "." << setfill ('0') << setw(3) << swap_actual;
      swap_file [swap_actual] .ss_descriptor = 
	open (name .str () .c_str (), O_TRUNC | O_CREAT | O_RDWR, S_IRWXU);
      if (swap_file [swap_actual] .ss_descriptor < 0) {
	cout << "Error opening swap file " << name .str () << endl;
	exit (0);
      }
      swap_file [swap_actual] .ss_name = name .str ();
    }
  }
  //
  //	Write chunk 
  //
  int 	fd = swap_file [chunk ->sc_file] .ss_descriptor;
  long 	position = chunk ->sc_position;
  lseek (fd, position, SEEK_SET);
  if (write (fd, chunk -> sc_storage, size) != size) {
    cout << "Error writing swap file "
	 << swap_file [chunk -> sc_file] .ss_name << endl;
    exit (0);
  }
  chunk ->sc_status |= SC_SWAPPED;
}
//
//____________________________________________________________________________
void * Storage::storage () const
{
  //
  //	Returns a pointer to actual memory (swapping it from file if needed)
  //	and marks memory chunk as active.
  //
  if (s_chunk == 0 || s_chunk -> sc_size == 0) return 0;
  if (s_chunk -> sc_storage == 0) restore (s_chunk);
  //
  // set active (remove unused flag)
  //
  s_chunk -> sc_status &= SC_USED;	
  return s_chunk -> sc_storage;
}
//
//____________________________________________________________________________
void Storage::storage (size_t size)
{
  //
  //	Creates a new chunk of memory (removing actual)
  //
  unlink ();
  create (size);
}
//
//____________________________________________________________________________
void Storage::swap_garbage () const
{
  //
  //	remove from swap files removed chunks
  //
  Storage_chunk * chunk;
  long position = 0;
  swap_actual = 1;
  swap_total  = 0.00;
  //
  //	Swapped chunks (marked by a positive sc_file field)
  //	are ordered in ascending order of file position.
  //
  for (chunk = storage_first; chunk; chunk = chunk -> sc_next) {
    if ((chunk -> sc_file <= 0) || (chunk -> sc_size == 0)) continue;
    long size = chunk -> sc_size;
    if (position + size > swap_top) {
      swap_actual++;
      position = 0;
    }
    if ((chunk -> sc_file > (long) swap_actual) ||
	(chunk -> sc_position > position)) {
      //
      //  Buffered copy from old position to new position
      //
      int fdst = swap_file [swap_actual] 	.ss_descriptor;
      int fsrc = swap_file [chunk -> sc_file] 	.ss_descriptor;
      char buffer [4096];
      long oldp = chunk -> sc_position;
      long newp = position;
      long endp = position + size;
      while (newp < endp) {
	long count = 4096;
	if (count > endp - newp)  count = endp - newp;
	lseek (fsrc, oldp, SEEK_SET);
	if (read (fsrc, buffer, count) != count) {
	  cout << "Error reading in swap_garbage!" << endl;
	  exit (0);
	}
	lseek (fdst, newp, SEEK_SET);
	if (write (fdst, buffer, count) != count) {
	  cout << "Error writing in swap_garbage!" << endl;
	  exit (0);
	}
	oldp += count;
	newp += count;
      }
    }
    chunk ->sc_file = swap_actual;
    chunk ->sc_position = position;
    position += size;
    swap_file [swap_actual] .ss_position = position;
    swap_total += size; 
  }
  swap_removed = 0.0;
}
//
//____________________________________________________________________________
void Storage::swap_unlink () const
{
  //
  //	Close and delete used swap files 
  //	(This is meaningfull only on exiting from the program)
  //
  for (size_t sf = 0; sf < swap_file .size (); sf++) 
    if (swap_file [sf] .ss_descriptor > 0) {
      close  (swap_file [sf] .ss_descriptor);
      ::unlink (swap_file [sf] .ss_name .c_str ());
    }
}
//
//____________________________________________________________________________
void Storage::top (size_t size) const
{
  //
  //	Set the maximum dynamic storage to size megabytes
  //
  memory_top = size * megab;
}
//
//____________________________________________________________________________
void Storage::unlink ()
{
  //
  //	Removes a link to a chunk of memory and the chunk itself if no more 
  //	links to the chunk
  //
  if (s_chunk && --(s_chunk ->sc_links) == 0 ) {
    //
    //	remove chunk from list of chunks
    //
    if (s_chunk ->sc_prev) s_chunk ->sc_prev ->sc_next = s_chunk ->sc_next;
    else		   storage_first = s_chunk ->sc_next;
    if (s_chunk ->sc_next) s_chunk ->sc_next ->sc_prev = s_chunk ->sc_prev;
    else		   storage_last = s_chunk ->sc_prev;
    //
    //	remove allocated memory
    //
    if (s_chunk -> sc_storage) {
      delete [] (char *) s_chunk ->sc_storage;
      memory_use -= s_chunk ->sc_size;
    }
    //
    //	If chunk was saved in swap file the corresponding area is garbage
    //
    if (s_chunk ->sc_status & SC_SWAPPED) swap_removed += s_chunk ->sc_size;
    //
    //	remove Storage_chunk 
    //
    delete s_chunk;
    s_chunk = 0;
  }
}
//
//____________________________________________________________________________
string Storage::usage () const
{
  //
  //	Return a string with a statistic on memory (swap) used
  //	in Mb units as ratios of actual to maximum.
  //	For memory reports resident+released / maximum used
  //	For swap reports the maximum active written / maximum written
  //	The maximum active written gives a rough indication of memory 
  //	needed for a run without swapping (un upper bound to the amount 
  //	to add to maximum used).	
  //
  stringstream usestr;
  double ma = ((double) memory_use)  / ((double) megab) + 0.5;
  double mt = ((double) memory_max)  / ((double) megab) + 0.5;
  double mv = ((double) released ()) / ((double) megab) + 0.5;
  double sa = swap_active / megab + 0.5;
  double st = swap_max    / megab + 0.5;
  usestr << fixed << setprecision(0) 
	 << right << ma - mv << "+" << left << mv << "/" << left << mt 
	 << " (" << right << sa << "/" << left << st << ") Mb";
  return usestr .str ();
}
//
//============================================================================
