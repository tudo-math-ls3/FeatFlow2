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

void coproc_checkErrors(char *label);

extern "C"
{
  int coproc_init();
  int FNAME(coproc_init)();
}

#endif
