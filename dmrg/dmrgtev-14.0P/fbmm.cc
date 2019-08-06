//	bfmm.cc				(C) 2015 Fabio Ortolani		150308
//============================================================================
//
//      NP:     Maximum number of generated threads
//              Hint: set as a multiple of actual physical cores.
//              (hyperthreading is useless)
//
//      NB:     Dimension of inner matrix blocks
//
#define NP              16
#define NB              224
//
//      CLNB:           lenght in bytes of level 1 cache line
//      CLND:           lengh in doubles of level 1 cache line
//                      NB must be a multiple of CLND for better performance
//
//      CSHB:           log2 (CLNB)
//      CSHD:           log2 (CLND)
//
#define CLNB            256
#define CLND            32
#define CSHB            6
#define CSHD            3
//
//      Clined (n):     nearest next mutiple of CLND
//      Calign (p):     nearest next aligned to CLND
//
#define Clined(_n)  ((((_n) + CLND - 1) >> CSHD) << CSHD)
#define Calign(_p) ((double *)(((size_t)(_p) >> CSHB) << CSHB) + (CLND))
//
extern "C" {
#include <pthread.h>
}
#ifdef __SSE3__
extern "C" {
#include <pmmintrin.h>
}
#endif // __SSE3__
#include <string.h>
#include <cstdio>
//____________________________________________________________________________
void 	    kmmtn111	(const long, const long,
			 const long, const long, const long,
			 const double *, const long,
			 const double *, const long,
			 double *, const long);
void 	    kmmtn218	(const long, const long,
			 const long, const long, const long,
			 const double *, const long,
			 const double *, const long,
			 double *, const long);
//____________________________________________________________________________
struct mmarg {
  long			md;
  long			nd;
  const double *	a;
  long			la;
  const double *	b;
  long			lb;
  double *		c;
  long			lc;
};
//
struct tharg {
  long			th;
  long			ta;
  long 			tb;
  long			kd;
  double		alpha;
  mmarg *		ma;
  long			na;
};
//
struct mpool {
  long			tr;
  double		factor;
  const double *	ms;
  double *		ma;
  double *	       	mm;
  struct mpool *	next;
  //
  mpool  ();
  ~mpool ();
};
//
typedef void (* mmproc)	(const long, const long,
			 const long, const long, const long,
			 const double *, const long,
			 const double *, const long,
			 double *, const long);
void mmcopy (const long, double, const double *, const long, double *,
	     const long, const long, const long);
//
//____________________________________________________________________________
long 			Nb	= NB;
long			Np	= NP;
mmproc			kmmtn	= 0;
//____________________________________________________________________________
static mpool *		pool	= 0;
static long		poold	= 0;
static pthread_mutex_t	poolmtx	= PTHREAD_MUTEX_INITIALIZER;
static tharg		ptarg 	[NP];
static pthread_t	thread 	[NP];
//____________________________________________________________________________
mpool::mpool ()
  : tr	(0), factor (0.0), ms (0), ma (0), mm (0), next (0)
{
  ma = new double [Nb * Nb + CLND];
  if (ma == 0) {
    fprintf (stderr, "mpool::mpool(): No room for matrix %ldx%ld\n", Nb, Nb);
    return;
  }
  mm = Calign (ma);
}
//____________________________________________________________________________
mpool::~mpool ()
{
  if (ma) delete [] ma;
  ms = ma = mm = 0;
}
//____________________________________________________________________________
void mpool_clear ()
{
  mpool * p = pool;
  while (p) {
    pool = p -> next;
    delete p;
    p = pool;
  }
  pool = 0;
}
//____________________________________________________________________________
void mpool_release (mpool * p)
{
  p ->tr     = 0;
  p ->factor = 0.0;
  p ->ms     = 0;
}
//____________________________________________________________________________
void mpool_release ()
{
  mpool * p = pool;
  while (p) {
    p ->tr = 0;
    p ->factor = 0.0;
    p ->ms = 0;
    p = p ->next;
  }
}
//____________________________________________________________________________
mpool * mpool_set (long tr, double factor, const double * mm,
		   long lm, long m, long n)
{
  mpool * found = pool;
  while (found) {
    if ((found ->tr == tr) && (found ->factor == factor) &&
	(found ->ms == mm)) break;
    found = found ->next;
  }
  if (found == 0) {
    found = pool;
    while (found) {
      if (found -> ms == 0) break;
      found = found -> next;
    }
    if (found == 0) {
      found = new mpool;
      if (found) {
	found ->next = pool;
	pool = found;
      }
    }
    if (found) {
      found ->tr = tr;
      found ->factor = factor;
      found ->ms = mm;
      //mpool_copy (tr, mm, lm, found ->mm, Nb, m, n);
      mmcopy (tr, factor, mm, lm, found-> mm, Nb, m, n);  
    }
  } 
  return found;
}
//_____________________________________________________________________________
void mmcopy (const long tr, double factor,
	     const double * __restrict__ src, const long ls,
	     double * __restrict__ dst, const long ld,
             const long m, const long n)
{
  long i, j, ii, jj, ir, jr;
  if (tr) {
    for (jj = 0; jj < n; jj += 8) {
      jr = jj + 8;
      if (jr > n) jr = n; 
      for (ii = 0; ii < m; ii += 8) {
        ir = ii + 8;
        if (ir > m) ir = m;
        for (j = jj; j < jr; j++)
          for (i = ii; i < ir; i++)
            dst [i + j * ld ] = factor * src [j + i * ls];
      }
    }
  }
  else {
    for (j = 0; j < n; j++) {
      memcpy (dst + j * ld, src + j * ls, m * sizeof (double));
      if (factor != 1.0)
	for (i = 0; i < m; i++) dst [i + j * ld] *= factor;
    }
  }
}
//____________________________________________________________________________
void * mmwork (void * dum)
{
  tharg * 	targ 	= static_cast<tharg *> (dum);
  long 		th	= targ ->th;
  long 		ta 	= targ ->ta;
  long		tb 	= targ ->tb;
  long  	kd 	= targ ->kd;
  double 	alpha 	= targ ->alpha;
  mmarg *	arg	= targ ->ma;
  long		na 	= targ ->na;
  long nc, m, n, kk, k, ka, kb, la, lb, lc;
  const double * a, * b;
  double * c;
  mpool * ap, * bp, * cp;
  //
  if (ta) 	ta = 0;
  else		ta = 1;
  for (nc = 0; nc < na; nc++) {
    m	= arg [nc] .md;
    n 	= arg [nc] .nd;
    a	= arg [nc] .a;
    la  = arg [nc] .la;
    b 	= arg [nc] .b;
    lb	= arg [nc] .lb;
    c	= arg [nc] .c;
    lc 	= arg [nc] .lc;
    pthread_mutex_lock 	 (& poolmtx);
    cp = mpool_set (0, 1.0, c, lc, m, n);
    pthread_mutex_unlock (& poolmtx);
    for (kk = 0; kk < kd; kk += Nb) {
      k = Nb;
      if ((kk + Nb) > kd) k = kd - kk;
      ka = kk;
      if (ta)  ka = kk * la;
      kb = kk;
      if (tb) kb = kk * lb;
      pthread_mutex_lock 	(& poolmtx);
      ap = mpool_set (ta, alpha, a + ka, la, k, m);
      bp = mpool_set (tb, 1.0,   b + kb, lb, k, n);
      pthread_mutex_unlock 	(& poolmtx);
      if (ap == 0) {
	fprintf (stderr, "Thread %ld at %ld: no room for block a\n", th, nc);
	return (void *) pool;
      }
      if (bp == 0) {
	fprintf (stderr, "Thread %ld at %ld: no room for block b\n", th, nc);
	return (void *) pool;
      }
      kmmtn (ta, tb, m, n, k, ap ->mm, Nb, bp ->mm, Nb, cp ->mm, Nb);
    }
    mmcopy (0, 1.0, cp ->mm, Nb, c, lc, m, n);
  }
  return 0;
}
//____________________________________________________________________________
void fbmm (const long ta, const long tb,
	   const long md, const long nd,  const long kd,
	   double   alpha,
	   const double * __restrict__ a, const long la,
	   const double * __restrict__ b, const long lb,
	   double *       __restrict__ c, const long lc)
{
  if (poold != Nb) {
    poold = Nb;
    mpool_clear ();
  }
  if (kmmtn == 0) kmmtn = kmmtn218;
  long 	mbl  	= (md + Nb -1) / Nb;
  long  nbl  	= (nd + Nb -1) / Nb;
  long  narg 	= mbl * nbl;
  mmarg * arg	= new mmarg [narg];
  long nc, nt, mm, m, nn, n, nth;
  nc = 0;
  for (nn = 0; nn < nd; nn += Nb) {
    n = Nb;
    if (nn + Nb > nd) n = nd - nn; 
    for (mm = 0; mm < md; mm += Nb, nc++) {
      m = Nb;
      if (mm + Nb > md) m = md - mm;
      arg [nc] .md = m;
      arg [nc] .nd = n;
      if (ta)	arg [nc] .a = a + mm * la;
      else	arg [nc] .a = a + mm;
      arg [nc] .la = la;
      if (tb) 	arg [nc] .b = b + nn;
      else	arg [nc] .b = b + nn * lb;
      arg [nc] .lb = lb;
      arg [nc] .c  = c + mm + nn * lc;
      arg [nc] .lc = lc;
    }
  }
  nth = Np;
  if (nc < nth) nth = nc;
  long n0 = nc / nth;
  long k0 = nc % nth;
  mmarg * marg = arg;
  for (nt = 0; nt < nth; nt++) {
    ptarg [nt] .th = nt;
    ptarg [nt] .ta = ta;
    ptarg [nt] .tb = tb;
    ptarg [nt] .kd = kd;
    ptarg [nt] .alpha = alpha;
    ptarg [nt] .ma = marg;
    if (nt < k0) {
      ptarg [nt] .na = n0 + 1;
      marg += (n0 + 1);
    }
    else {
      ptarg [nt] .na = n0;
      marg += n0;
    }
  }
  if (nth == 1) 	mmwork (ptarg);
  else {
    for (nt = 0; nt < nth; nt++)
      pthread_create (thread + nt, 0, mmwork, ptarg + nt);
    for (nt = 0; nt < nth; nt++)
      pthread_join (thread [nt], 0);
  }
  mpool_release ();
  delete [] arg;
}
//____________________________________________________________________________
void kmmtn111   (const long ta, const long tb,
		 const long md,  const long nd, const long kd,
		 const double * __restrict__ a, const long la,
		 const double * __restrict__ b, const long lb,
		 double * __restrict__ c, const long lc)
{
  const double *	enda =  a + la * md;
  const double *	endb =  b + lb * nd;
  long   		sam  =  la - kd;
  long   		sbm  = -kd;
  long   		san  = -la * md;
  long   		sbn  =  lb;
  long   		scn  =  lc - md;
  const double *	pa0  =  a;
  const double *	pb0  =  b;
  double *      	pc0  =  c;
  register long  	k;
  register double 	a0, b0, c00;
  do {  			// n loop
    do {  			// m loop
      c00 = *pc0;
      for (k = kd; k; k--) {  	// k loop
        a0   = *pa0;
        b0   = *pb0;
        c00 +=  a0 * b0;
        pa0++;
        pb0++;
      } 			// k end
      *pc0  = c00;
      pa0  += sam;
      pb0  += sbm;
      pc0++;
    } while (pa0 != enda);      // m end
    pa0 += san;
    pb0 += sbn;
    pc0 += scn;
  } while (pb0 != endb);        // n end
}
#ifdef __SSE3__
//____________________________________________________________________________
void kmmtn218 	(const long ta, const long tb,
		 const long md,  const long nd, const long kd,
		 const double * __restrict__ a, const long la,
		 const double * __restrict__ b, const long lb,
		 double * __restrict__ c, const long lc)
{
  long          	mb      =  (md >> 1) << 1;
  long          	mr      =  (md - mb);
  long          	kb      =  (kd >> 3) << 3;
  long          	kr      =  (kd - kb);
  const double *        enda    =  a + (la * mb);
  const double *        endb    =  b + (lb * nd);
  long          	sam     =  2 * la - kb;
  long          	sbm     = -kb;
  long          	san     = -(la * mb);
  long          	sbn     =  lb;
  long          	scn     =  (lc - mb);
  const double *        pa0     =  a;
  const double *        pa1     =  a + la;
  const double *        pb0     =  b;
  double *              pc0     =  c;
  register long         k;
  register __m128d      a0, a1, b0, c00;
  do {                                                          // n loop
    if (mb) do {                                                // m loop
      c00 = _mm_load_pd (pc0 + 0); 
      for (k = kb; k; k -= 8) {                                 // k loop
        b0   = _mm_load_pd (pb0 + 0);     
        a0   = _mm_load_pd (pa0 + 0);     
        a1   = _mm_load_pd (pa1 + 0);     
        a0   = _mm_mul_pd  (a0, b0);  
        a1   = _mm_mul_pd  (a1, b0);  
        b0   = _mm_hadd_pd (a0, a1);  
        c00  = _mm_add_pd  (b0, c00);
        //
        b0   = _mm_load_pd (pb0 + 2);     
        a0   = _mm_load_pd (pa0 + 2);     
        a1   = _mm_load_pd (pa1 + 2);     
        a0   = _mm_mul_pd  (a0, b0);  
        a1   = _mm_mul_pd  (a1, b0);  
        b0   = _mm_hadd_pd (a0, a1);  
        c00  = _mm_add_pd  (b0, c00);
        //
         b0   = _mm_load_pd (pb0 + 4);     
        a0   = _mm_load_pd (pa0 + 4);     
        a1   = _mm_load_pd (pa1 + 4);     
        a0   = _mm_mul_pd  (a0, b0);  
        a1   = _mm_mul_pd  (a1, b0);  
        b0   = _mm_hadd_pd (a0, a1);  
        c00  = _mm_add_pd  (b0, c00);
        //
        b0   = _mm_load_pd (pb0 + 6);     
        a0   = _mm_load_pd (pa0 + 6);     
        a1   = _mm_load_pd (pa1 + 6);     
        a0   = _mm_mul_pd  (a0, b0);  
        a1   = _mm_mul_pd  (a1, b0);  
        b0   = _mm_hadd_pd (a0, a1);  
        c00  = _mm_add_pd  (b0, c00);
        //
        pa0 += 8;
        pa1 += 8;
        pb0 += 8;
      }                                                         //  k end
      switch (kr) {
      case 7:
        b0   = _mm_load_pd (pb0 + 0);     
        a0   = _mm_load_pd (pa0 + 0);     
        a1   = _mm_load_pd (pa1 + 0);     
        a0   = _mm_mul_pd  (a0, b0);  
        a1   = _mm_mul_pd  (a1, b0);  
        b0   = _mm_hadd_pd (a0, a1);  
        c00  = _mm_add_pd  (b0, c00);
        //
        b0   = _mm_load_pd (pb0 + 2);     
        a0   = _mm_load_pd (pa0 + 2);     
        a1   = _mm_load_pd (pa1 + 2);     
        a0   = _mm_mul_pd  (a0, b0);  
        a1   = _mm_mul_pd  (a1, b0);  
        b0   = _mm_hadd_pd (a0, a1);  
        c00  = _mm_add_pd  (b0, c00);
        //
        b0   = _mm_load_pd (pb0 + 4);     
        a0   = _mm_load_pd (pa0 + 4);     
        a1   = _mm_load_pd (pa1 + 4);     
        a0   = _mm_mul_pd  (a0, b0);  
        a1   = _mm_mul_pd  (a1, b0);  
        b0   = _mm_hadd_pd (a0, a1);  
        c00  = _mm_add_pd  (b0, c00);
        //
        b0   = _mm_loaddup_pd (pb0 + 6);
        a0   = _mm_load_sd    (pa0 + 6);
        a0   = _mm_loadh_pd   (a0, pa1 + 6); 
        b0   = _mm_mul_pd     (a0, b0);  
        c00  = _mm_add_pd     (b0, c00); 
        break;
      case 6:
        b0   = _mm_load_pd (pb0 + 0);     
        a0   = _mm_load_pd (pa0 + 0);     
        a1   = _mm_load_pd (pa1 + 0);     
        a0   = _mm_mul_pd  (a0, b0);  
        a1   = _mm_mul_pd  (a1, b0);  
        b0   = _mm_hadd_pd (a0, a1);  
        c00  = _mm_add_pd  (b0, c00);
        //
        b0   = _mm_load_pd (pb0 + 2);     
        a0   = _mm_load_pd (pa0 + 2);     
        a1   = _mm_load_pd (pa1 + 2);     
        a0   = _mm_mul_pd  (a0, b0);  
        a1   = _mm_mul_pd  (a1, b0);  
        b0   = _mm_hadd_pd (a0, a1);  
        c00  = _mm_add_pd  (b0, c00);
        //
        b0   = _mm_load_pd (pb0 + 4);     
        a0   = _mm_load_pd (pa0 + 4);     
        a1   = _mm_load_pd (pa1 + 4);     
        a0   = _mm_mul_pd  (a0, b0);  
        a1   = _mm_mul_pd  (a1, b0);  
        b0   = _mm_hadd_pd (a0, a1);  
        c00  = _mm_add_pd  (b0, c00);
        break;
      case 5:
        b0   = _mm_load_pd (pb0 + 0);     
        a0   = _mm_load_pd (pa0 + 0);     
        a1   = _mm_load_pd (pa1 + 0);     
        a0   = _mm_mul_pd  (a0, b0);  
        a1   = _mm_mul_pd  (a1, b0);  
        b0   = _mm_hadd_pd (a0, a1);  
        c00  = _mm_add_pd  (b0, c00);
        //
        b0   = _mm_load_pd (pb0 + 2);     
        a0   = _mm_load_pd (pa0 + 2);     
        a1   = _mm_load_pd (pa1 + 2);     
        a0   = _mm_mul_pd  (a0, b0);  
        a1   = _mm_mul_pd  (a1, b0);  
        b0   = _mm_hadd_pd (a0, a1);  
        c00  = _mm_add_pd  (b0, c00);
        //
        b0   = _mm_loaddup_pd (pb0 + 4);
        a0   = _mm_load_sd    (pa0 + 4);
        a0   = _mm_loadh_pd   (a0, pa1 + 4); 
        b0   = _mm_mul_pd     (a0, b0);  
        c00  = _mm_add_pd     (b0, c00); 
        break;
      case 4:
        b0   = _mm_load_pd (pb0 + 0);     
        a0   = _mm_load_pd (pa0 + 0);     
        a1   = _mm_load_pd (pa1 + 0);     
        a0   = _mm_mul_pd  (a0, b0);  
        a1   = _mm_mul_pd  (a1, b0);  
        b0   = _mm_hadd_pd (a0, a1);  
        c00  = _mm_add_pd  (b0, c00);
        //
        b0   = _mm_load_pd (pb0 + 2);     
        a0   = _mm_load_pd (pa0 + 2);     
        a1   = _mm_load_pd (pa1 + 2);     
        a0   = _mm_mul_pd  (a0, b0);  
        a1   = _mm_mul_pd  (a1, b0);  
        b0   = _mm_hadd_pd (a0, a1);  
        c00  = _mm_add_pd  (b0, c00);
        //
        break;
      case 3:
        b0   = _mm_load_pd (pb0 + 0);     
        a0   = _mm_load_pd (pa0 + 0);     
        a1   = _mm_load_pd (pa1 + 0);     
        a0   = _mm_mul_pd  (a0, b0);  
        a1   = _mm_mul_pd  (a1, b0);  
        b0   = _mm_hadd_pd (a0, a1);  
        c00  = _mm_add_pd  (b0, c00);
        //
        b0   = _mm_loaddup_pd (pb0 + 2);
        a0   = _mm_load_sd    (pa0 + 2);
        a0   = _mm_loadh_pd   (a0, pa1 + 2); 
        b0   = _mm_mul_pd     (a0, b0);  
        c00  = _mm_add_pd     (b0, c00); 
        break;
      case 2:
        b0   = _mm_load_pd (pb0 + 0);     
        a0   = _mm_load_pd (pa0 + 0);     
        a1   = _mm_load_pd (pa1 + 0);     
        a0   = _mm_mul_pd  (a0, b0);  
        a1   = _mm_mul_pd  (a1, b0);  
        b0   = _mm_hadd_pd (a0, a1);  
        c00  = _mm_add_pd  (b0, c00);
        //
        break;
      case 1:
        b0   = _mm_loaddup_pd (pb0 + 0);
        a0   = _mm_load_sd    (pa0 + 0);
        a0   = _mm_loadh_pd   (a0, pa1 + 0); 
        b0   = _mm_mul_pd     (a0, b0);  
        c00  = _mm_add_pd     (b0, c00); 
        break;
      case 0:
        break;
      }
      _mm_store_pd (pc0, c00);
      pa0 += sam;
      pa1 += sam;
      pb0 += sbm;
      pc0 += 2;
    } while (pa0 != enda);                                      // m end
    if (mr) {
      for (k = 0; k < kd; k++) *pc0 += pa0 [k] * pb0 [k];
    }
    pa0 += san;
    pa1 += san;
    pb0 += sbn;
    pc0 += scn;
  } while (pb0 != endb);                                        // n end
}
#else // ! __SSE3__
#define kmmtn218 kmmtn111
#endif // __SSE3__
//============================================================================

    
  
      
