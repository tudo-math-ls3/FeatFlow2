/* ========================================================================== */
/* === umf4_f77wrapper_port ================================================= */
/* ========================================================================== */

/* *MK* */

/* This is a modified version of the umf4_f77wrapper.c file to make the library
   more portable. It introduces the following changes:

   - Every wrapper routine is available in both an uppercase and lowercase 
     version. Either one of them with and without trailing underline.
   - To support architectures with 64 bit pointers with 32 bit integers,
     matrices are identified passed by/to Fortran as handles while the
     routines in this file wrap integer-handles to the pointers of the
     matrices.
*/

/* -------------------------------------------------------------------------- */
/* UMFPACK Copyright (c) Timothy A. Davis, CISE,                              */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* FORTRAN interface for the C-callable UMFPACK library (double / int version
 * only and double / UF_long versions only).  This is HIGHLY non-portable.  You
 * will need to modify this depending on how your FORTRAN and C compilers
 * behave.  This has been tested in Linux, Sun Solaris, SGI IRIX, and IBM AIX,
 * with various compilers.  It has not been exhaustively tested on all possible
 * combinations of C and FORTRAN compilers.  The UF_long version works on
 * Solaris, SGI IRIX, and IBM AIX when the UMFPACK library is compiled in
 * 64-bit mode.
 *
 * Only a subset of UMFPACK's capabilities are provided.  Refer to the UMFPACK
 * User Guide for details.
 *
 * For some C and FORTRAN compilers, the FORTRAN compiler appends a single
 * underscore ("_") after each routine name.  C doesn't do this, so the
 * translation is made here.  Other FORTRAN compilers treat underscores
 * differently.  For example, a FORTRAN call to a_b gets translated to a call
 * to a_b__ by g77, and to a_b_ by most other FORTRAN compilers.  Thus, the
 * FORTRAN names here do not use underscores.  The xlf compiler in IBM AIX
 * doesn't add an underscore.
 *
 * The matrix A is passed to UMFPACK in compressed column form, with 0-based
 * indices.  In FORTRAN, for an m-by-n matrix A with nz entries, the row
 * indices of the first column (column 1) are in Ai (Ap (1) + 1 ... Ap (2)),
 * with values in Ax (Ap (1) + 1 ... Ap (2)).  The last column (column n) is
 * in Ai (Ap (n) + 1 ... Ap (n+1)) and Ax (Ap (n) + 1 ... Ap (n+1)).  The row
 * indices in Ai are in the range 0 to m-1.  They must be sorted, with no
 * duplicate entries allowed.  Refer to umfpack_di_triplet_to_col for a more
 * flexible format for the input matrix.  The following defintions apply
 * for each of the routines in this file:
 *
 *	integer m, n, Ap (n+1), Ai (nz), symbolic, numeric, filenum, status
 *	double precision Ax (nz), control (20), info (90), x (n), b (n)
 *
 * UMFPACK's status is returned in either a status argument, or in info (1).
 * It is zero if everything is OK, 1 if the matrix is singular (this is a
 * warning, not an error), and negative if an error occurred.  See umfpack.h
 * for more details on the contents of the control and info arrays, and the
 * value of the sys argument.
 *
 * For the Numeric and Symbolic handles, it's probably safe to assume that a
 * FORTRAN integer is sufficient to store a C pointer.  If that doesn't work,
 * try defining numeric and symbolic as integer arrays of size 2, or as
 * integer*8, in the FORTRAN routine that calls these wrapper routines.
 * The latter is required on Solaris, SGI IRIX, and IBM AIX when UMFPACK is
 * compiled in 64-bit mode.
 *
 * If you want to use 64-bit integers, try compiling this file with the -DDLONG
 * compiler option (via "make fortran64").  First modify your UFconfig.mk
 * file to compile UMFPACK in LP64 mode (see the User Guide for details).
 * Your FORTRAN code should use integer*8.  See umf4hb64.f for an example.
 *
 * Tested with the following compilers:
 *	* Solaris with cc and f77 from Sun WorkShop 6 update 1
 *	    (32-bit and 64-bit modes)
 *	* SGI Irix with MIPSpro cc and f77 compilers version 7.4
 *	    (32-bit and 64-bit modes)
 *	* Linux with GNU gcc and Intel's icc, and GNU g77 and Intel's
 *	    ifc FORTRAN compiler.  See the comments above about g77 and
 *	    underscores.  Only supports 32-bit mode.
 *	* IBM AIX xlc and xlf compilers.
 *	    (32-bit and 64-bit modes)
 *
 * This interface files when using the ABSOFT Fortran compiler on both
 * Linux and Windows.  There is no known fix.
 */

#include "umfpack.h"
#include <ctype.h>
#include <stdio.h>
#ifdef NULL
#undef NULL
#endif
#define NULL 0
#define LEN 200

/* -------------------------------------------------------------------------- */
/* integer type: int or UF_long */
/* -------------------------------------------------------------------------- */

#if defined (DLONG)

#define Int UF_long
#define UMFPACK_defaults	 umfpack_dl_defaults
#define UMFPACK_free_numeric	 umfpack_dl_free_numeric
#define UMFPACK_free_symbolic	 umfpack_dl_free_symbolic
#define UMFPACK_numeric		 umfpack_dl_numeric
#define UMFPACK_report_control	 umfpack_dl_report_control
#define UMFPACK_report_info	 umfpack_dl_report_info
#define UMFPACK_save_numeric	 umfpack_dl_save_numeric
#define UMFPACK_save_symbolic	 umfpack_dl_save_symbolic
#define UMFPACK_load_numeric	 umfpack_dl_load_numeric
#define UMFPACK_load_symbolic	 umfpack_dl_load_symbolic
#define UMFPACK_scale		 umfpack_dl_scale
#define UMFPACK_solve		 umfpack_dl_solve
#define UMFPACK_symbolic	 umfpack_dl_symbolic

#else

#define Int int
#define UMFPACK_defaults	 umfpack_di_defaults
#define UMFPACK_free_numeric	 umfpack_di_free_numeric
#define UMFPACK_free_symbolic	 umfpack_di_free_symbolic
#define UMFPACK_numeric		 umfpack_di_numeric
#define UMFPACK_report_control	 umfpack_di_report_control
#define UMFPACK_report_info	 umfpack_di_report_info
#define UMFPACK_save_numeric	 umfpack_di_save_numeric
#define UMFPACK_save_symbolic	 umfpack_di_save_symbolic
#define UMFPACK_load_numeric	 umfpack_di_load_numeric
#define UMFPACK_load_symbolic	 umfpack_di_load_symbolic
#define UMFPACK_scale		 umfpack_di_scale
#define UMFPACK_solve		 umfpack_di_solve
#define UMFPACK_symbolic	 umfpack_di_symbolic

#endif

/* -------------------------------------------------------------------------- */
/* construct a file name from a file number (not user-callable) */
/* -------------------------------------------------------------------------- */

static void make_filename (Int filenum, char *prefix, char *filename)
{
    char *psrc, *pdst ;
#ifdef DLONG
    sprintf (filename, "%s%ld.umf", prefix, filenum) ;
#else
    sprintf (filename, "%s%d.umf", prefix, filenum) ;
#endif
    /* remove any spaces in the filename */
    pdst = filename ;
    for (psrc = filename ; *psrc ; psrc++)
    {
	if (!isspace (*psrc)) *pdst++ = *psrc ;
    }
    *pdst = '\0' ;
}

/* ========================================================================== */
/* wrapper routines, standard version: lowercase without underscore           */
/* ========================================================================== */

/* first we define a static pointer array which collects all pointers and
   wrap integer-handles to pointers                                           */

#define N_HANDLECOUNT (4096)

void* umf4handles [N_HANDLECOUNT];
int initialized = 0;

/* Some routines to work with the handles */

/* Clear all handles */
void ClearHandles () {
  int i;
  for (i=0; i < N_HANDLECOUNT; i++) {
    umf4handles [i] = NULL;
  }
}

/* Store a pointer, return a handle.
 * Return -1 if there's no more handle available! */
int StorePointer (void* obj) {
  int i;

  /* Check if the Handle-Array has to be initialized */
  if (initialized==0) {
    ClearHandles ();
    initialized = 1;
  }

  /* Search for the next available handle */
  i=0;
  while ((i < N_HANDLECOUNT) && (umf4handles [i] != NULL)) i++;

  if (i==N_HANDLECOUNT) {
    i = -1;                /* no handle available */
  } else {  
    umf4handles [i] = obj; /* store the handle */
  }
  return i;
}

/* To get a pointer from a handle: Look into the array umf4handles.
   To delete a pointer of a handle: Overwrite the entry in umf4handley by NULL. */

/* umf4def: set default control parameters */
/* call umf4def (control) */

void umf4def (double Control [UMFPACK_CONTROL])
{
    UMFPACK_defaults (Control) ;
}

/* -------------------------------------------------------------------------- */
/* umf4pcon: print control parameters */
/* -------------------------------------------------------------------------- */
/* call umf4pcon (control) */

void umf4pcon (double Control [UMFPACK_CONTROL])
{
    fflush (stdout) ;
    UMFPACK_report_control (Control) ;
    fflush (stdout) ;
}

/* -------------------------------------------------------------------------- */
/* umf4sym: pre-ordering and symbolic factorization */
/* -------------------------------------------------------------------------- */
/* call umf4sym (m, n, Ap, Ai, Ax, symbolic, control, info) */
/* The result of this function is the variable SYMBOLIC - a handle to the allocated
   symbolic object. */

void umf4sym (Int *m, Int *n, Int Ap [ ], Int Ai [ ], double Ax [ ], 
    Int *Symbolic,
    double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO])
{
    void* sym = NULL;
    (void) UMFPACK_symbolic (*m, *n, Ap, Ai, Ax, &sym, Control, Info) ;

    if (Info [0] == UMFPACK_OK) {
      *Symbolic = StorePointer (sym);
    } else *Symbolic = -1;
}

/* -------------------------------------------------------------------------- */
/* umf4num: numeric factorization */
/* -------------------------------------------------------------------------- */

/* call umf4num (Ap, Ai, Ax, symbolic, numeric, control, info) */

void umf4num (Int Ap [ ], Int Ai [ ], double Ax [ ],
    Int *Symbolic, Int *Numeric,
    double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO])
{
    void* num = NULL;
    (void) UMFPACK_numeric (Ap, Ai, Ax, umf4handles [*Symbolic], &num, 
                            Control, Info);

    if (Info [0] == UMFPACK_OK) {
      *Numeric = StorePointer (num);
    } else *Numeric = -1;
}

/* -------------------------------------------------------------------------- */
/* umf4solr: solve a linear system with iterative refinement */
/* -------------------------------------------------------------------------- */

/* call umf4solr (sys, Ap, Ai, Ax, x, b, numeric, control, info) */

void umf4solr (Int *sys, Int Ap [ ], Int Ai [ ], double Ax [ ],
    double x [ ], double b [ ], Int *Numeric,
    double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO])
{
    (void) UMFPACK_solve (*sys, Ap, Ai, Ax, x, b, umf4handles [*Numeric], 
                          Control, Info) ;
}

/* -------------------------------------------------------------------------- */
/* umf4sol: solve a linear system without iterative refinement */
/* -------------------------------------------------------------------------- */

/* call umf4sol (sys, x, b, numeric, control, info) */

void umf4sol (Int *sys, double x [ ], double b [ ], Int *Numeric,
    double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO])
{
    Control [UMFPACK_IRSTEP] = 0 ;
    (void) UMFPACK_solve (*sys, (Int *) NULL, (Int *) NULL, (double *) NULL,
	                        x, b, umf4handles [*Numeric], Control, Info) ;
}

/* -------------------------------------------------------------------------- */
/* umf4scal: scale a vector using UMFPACK's scale factors */
/* -------------------------------------------------------------------------- */

/* call umf4scal (x, b, numeric, status) */

void umf4scal (double x [ ], double b [ ], Int *Numeric, Int *status)
{
    *status = UMFPACK_scale (x, b, umf4handles [*Numeric]) ;
}

/* -------------------------------------------------------------------------- */
/* umf4pinf: print info */
/* -------------------------------------------------------------------------- */

/* call umf4pinf (control) */

void umf4pinf (double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO])
{
    fflush (stdout) ;
    UMFPACK_report_info (Control, Info) ;
    fflush (stdout) ;
}

/* -------------------------------------------------------------------------- */
/* umf4fnum: free the Numeric object */
/* -------------------------------------------------------------------------- */

/* call umf4fnum (numeric) */

void umf4fnum (Int *Numeric)
{
    UMFPACK_free_numeric (&umf4handles [*Numeric]) ;
    umf4handles [*Numeric] = NULL;
}

/* -------------------------------------------------------------------------- */
/* umf4fsym: free the Symbolic object */
/* -------------------------------------------------------------------------- */

/* call umf4fsym (symbolic) */

void umf4fsym (Int *Symbolic)
{
    UMFPACK_free_symbolic (&umf4handles [*Symbolic]) ;
    umf4handles [*Symbolic] = NULL;
}

/* -------------------------------------------------------------------------- */
/* umf4snum: save the Numeric object to a file */
/* -------------------------------------------------------------------------- */

/* call umf4snum (numeric, filenum, status) */

void umf4snum (Int *Numeric, Int *filenum, Int *status)
{
    char filename [LEN] ;
    make_filename (*filenum, "n", filename) ;
    *status = UMFPACK_save_numeric (umf4handles [*Numeric], filename) ;
}

/* -------------------------------------------------------------------------- */
/* umf4ssym: save the Symbolic object to a file */
/* -------------------------------------------------------------------------- */

/* call umf4ssym (symbolic, filenum, status) */

void umf4ssym (Int *Symbolic, Int *filenum, Int *status)
{
    char filename [LEN] ;
    make_filename (*filenum, "s", filename) ;
    *status = UMFPACK_save_symbolic (umf4handles [*Symbolic], filename) ;
}

/* -------------------------------------------------------------------------- */
/* umf4lnum: load the Numeric object from a file */
/* -------------------------------------------------------------------------- */

/* call umf4lnum (numeric, filenum, status) */

void umf4lnum (Int *Numeric, Int *filenum, Int *status)
{
    char filename [LEN];
    void* num;
    make_filename (*filenum, "n", filename);
    *status = UMFPACK_load_numeric (&num, filename);
    if (*status == UMFPACK_OK) {
      *Numeric = StorePointer (num);
    } else *Numeric = -1;
}

/* -------------------------------------------------------------------------- */
/* umf4lsym: load the Symbolic object from a file */
/* -------------------------------------------------------------------------- */

/* call umf4lsym (symbolic, filenum, status) */

void umf4lsym (Int *Symbolic, Int *filenum, Int *status)
{
    char filename [LEN] ;
    void* sym;
    make_filename (*filenum, "s", filename) ;
    *status = UMFPACK_load_symbolic (&sym, filename);
    if (*status == UMFPACK_OK) {
      *Symbolic = StorePointer (sym);
    } else *Symbolic = -1;
}

/* -------------------------------------------------------------------------- */
/* uppercase variants */
/* -------------------------------------------------------------------------- */

void UMF4DEF (double Control [UMFPACK_CONTROL]) {
  umf4def (Control);
}

void UMF4PCON (double Control [UMFPACK_CONTROL]) {
  umf4pcon (Control);
}

void UMF4SYM (Int *m, Int *n, Int Ap [ ], Int Ai [ ], double Ax [ ], Int *Symbolic,
              double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]) {
  umf4sym (m, n, Ap, Ai, Ax, Symbolic,Control, Info);
}

void UMF4NUM (Int Ap [ ], Int Ai [ ], double Ax [ ],  Int *Symbolic, Int *Numeric,
              double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]) {
  umf4num (Ap, Ai, Ax, Symbolic, Numeric, Control, Info);
}

void UMF4SOLR (Int *sys, Int Ap [ ], Int Ai [ ], double Ax [ ],
    double x [ ], double b [ ], Int *Numeric,
    double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]) {
  umf4solr (sys, Ap, Ai, Ax, x, b, Numeric, Control, Info);
}

void UMF4SOL (Int *sys, double x [ ], double b [ ], Int *Numeric,
              double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]) {
  umf4sol (sys, x, b, Numeric, Control, Info);
}

void UMF4SCAL (double x [ ], double b [ ], Int *Numeric, Int *status) {
  umf4scal (x, b, Numeric, status);
}

void UMF4PINF (double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]) {
  umf4pinf (Control, Info);
}

void UMF4FNUM (Int *Numeric) {
  umf4fnum (Numeric);
}

void UMF4FSYM (Int *Symbolic) {
  umf4fsym (Symbolic);
}

void UMF4SNUM (Int *Numeric, Int *filenum, Int *status) {
  umf4snum (Numeric, filenum, status);
}

void UMF4SSYM (Int *Symbolic, Int *filenum, Int *status) {
  umf4ssym (Symbolic, filenum, status);
}
void UMF4LNUM (Int *Numeric, Int *filenum, Int *status) {
  umf4lnum (Numeric, filenum, status);
}

void UMF4LSYM (Int *Symbolic, Int *filenum, Int *status) {
  umf4lsym (Symbolic, filenum, status);
}

/* -------------------------------------------------------------------------- */
/* lowercase, with underscore */
/* -------------------------------------------------------------------------- */

void umf4def_ (double Control [UMFPACK_CONTROL]) {
  umf4def (Control);
}

void umf4pcon_ (double Control [UMFPACK_CONTROL]) {
  umf4pcon (Control);
}

void umf4sym_ (Int *m, Int *n, Int Ap [ ], Int Ai [ ], double Ax [ ], Int *Symbolic,
              double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]) {
  umf4sym (m, n, Ap, Ai, Ax, Symbolic,Control, Info);
}

void umf4num_ (Int Ap [ ], Int Ai [ ], double Ax [ ],  Int *Symbolic, Int *Numeric,
              double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]) {
  umf4num (Ap, Ai, Ax, Symbolic, Numeric, Control, Info);
}

void umf4solr_ (Int *sys, Int Ap [ ], Int Ai [ ], double Ax [ ],
    double x [ ], double b [ ], Int *Numeric,
    double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]) {
  umf4solr (sys, Ap, Ai, Ax, x, b, Numeric, Control, Info);
}

void umf4sol_ (Int *sys, double x [ ], double b [ ], Int *Numeric,
              double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]) {
  umf4sol (sys, x, b, Numeric, Control, Info);
}

void umf4scal_ (double x [ ], double b [ ], Int *Numeric, Int *status) {
  umf4scal (x, b, Numeric, status);
}

void umf4pinf_ (double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]) {
  umf4pinf (Control, Info);
}

void umf4fnum_ (Int *Numeric) {
  umf4fnum (Numeric);
}

void umf4fsym_ (Int *Symbolic) {
  umf4fsym (Symbolic);
}

void umf4snum_ (Int *Numeric, Int *filenum, Int *status) {
  umf4snum (Numeric, filenum, status);
}

void umf4ssym_ (Int *Symbolic, Int *filenum, Int *status) {
  umf4ssym (Symbolic, filenum, status);
}
void umf4lnum_ (Int *Numeric, Int *filenum, Int *status) {
  umf4lnum (Numeric, filenum, status);
}

void umf4lsym_ (Int *Symbolic, Int *filenum, Int *status) {
  umf4lsym (Symbolic, filenum, status);
}

/* -------------------------------------------------------------------------- */
/* uppercase, with underscore */
/* -------------------------------------------------------------------------- */

void UMF4DEF_ (double Control [UMFPACK_CONTROL]) {
  umf4def (Control);
}

void UMF4PCON_ (double Control [UMFPACK_CONTROL]) {
  umf4pcon (Control);
}

void UMF4SYM_ (Int *m, Int *n, Int Ap [ ], Int Ai [ ], double Ax [ ], Int *Symbolic,
              double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]) {
  umf4sym (m, n, Ap, Ai, Ax, Symbolic,Control, Info);
}

void UMF4NUM_ (Int Ap [ ], Int Ai [ ], double Ax [ ],  Int *Symbolic, Int *Numeric,
              double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]) {
  umf4num (Ap, Ai, Ax, Symbolic, Numeric, Control, Info);
}

void UMF4SOLR_ (Int *sys, Int Ap [ ], Int Ai [ ], double Ax [ ],
    double x [ ], double b [ ], Int *Numeric,
    double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]) {
  umf4solr (sys, Ap, Ai, Ax, x, b, Numeric, Control, Info);
}

void UMF4SOL_ (Int *sys, double x [ ], double b [ ], Int *Numeric,
              double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]) {
  umf4sol (sys, x, b, Numeric, Control, Info);
}

void UMF4SCAL_ (double x [ ], double b [ ], Int *Numeric, Int *status) {
  umf4scal (x, b, Numeric, status);
}

void UMF4PINF_ (double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]) {
  umf4pinf (Control, Info);
}

void UMF4FNUM_ (Int *Numeric) {
  umf4fnum (Numeric);
}

void UMF4FSYM_ (Int *Symbolic) {
  umf4fsym (Symbolic);
}

void UMF4SNUM_ (Int *Numeric, Int *filenum, Int *status) {
  umf4snum (Numeric, filenum, status);
}

void UMF4SSYM_ (Int *Symbolic, Int *filenum, Int *status) {
  umf4ssym (Symbolic, filenum, status);
}
void UMF4LNUM_ (Int *Numeric, Int *filenum, Int *status) {
  umf4lnum (Numeric, filenum, status);
}

void UMF4LSYM_ (Int *Symbolic, Int *filenum, Int *status) {
  umf4lsym (Symbolic, filenum, status);
}

