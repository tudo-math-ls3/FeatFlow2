#ifndef _COPROC_CORE_H_
#define _COPROC_CORE_H_

/*#############################################################################
 ******************************************************************************
 * <name> coproc_core </name>
 ******************************************************************************
 *
 * <purpose>
 * This header file contains some core declarations for coprocessor support.
 * </purpose>
 *
 *#############################################################################
 */

using namespace std;

/*******************************************************************************
 * The macro FNAME converts a C function name to the form that is
 * required when it is called by a FORTRAN program. On many UNIX
 * systems, subprograms called from FORTRAN have an implied underscore
 * character at the ends of their names. This macro takes care of this
 * operating system quirk.
 *******************************************************************************/
#ifdef VMS
#define FNAME(name)	name

#else

#ifdef __APPLE__
#define FNAME(name)	name

#else

#ifdef __STDC__
#define FNAME(name)	name##_

#else
#define FNAME(name)	name/**/_

#endif
#endif
#endif

/*******************************************************************************
 * The following macros should ensure the interoperability of data types
 *******************************************************************************/

#define __SP      float
#define __DP      double
#ifdef ENABLE_QUADPREC
#define __QP      long double
#else
#define __QP      __DP
#endif
#define __I8      signed char
#define __I16     short
#define __I32     int
#define __I64     long long int
#ifdef ENABLE_LARGEINT
#define __INT     __I64
#else
#define __INT     __I32
#endif
#define __LOGICAL int
#define __CHAR    char

#if defined(_WIN64) || defined(_LP64) || defined(__LP64__)
// LP64 machine, OS X or Linux or Unix or
// LLP64 machine, Windows
#define __SIZET   __I64
#elif defined(_WIN32) || defined(_ILD32)
// 32-bit machine, Windows or Linux or OS X or Unix
#define __SIZET   __I32
#else
// Machine not detected uniquely assume 64-bit
#define __SIZET   __I64
#endif


/*******************************************************************************/

#define __coproc__error__(label)							\
  cerr << "CUDA Error: " << label << "(at " << __LINE__ << " of " __FILE__ << ")" << endl;

/*******************************************************************************/

extern "C"
{
  void coproc_checkErrors(__CHAR *label);

  int coproc_init(int deviceNumber);
  
  int coproc_getSizeOf(int cdatatype,
		       size_t isize);

  int coproc_createStream(cudaStream_t *pStream);
}

#endif
