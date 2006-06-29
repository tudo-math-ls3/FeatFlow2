/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdf.ncsa.uiuc.edu/HDF5/doc/Copyright.html.  If you do not have     *
 * access to either file, you may request a copy from hdfhelp@ncsa.uiuc.edu. *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 *  Header file for error values, etc.
 */
#ifndef _H5Eprivate_H
#define _H5Eprivate_H

#include "H5Epublic.h"

/* Private headers needed by this file */
#include "H5private.h"

#define H5E_NSLOTS	32	/*number of slots in an error stack	     */

/*
 * The list of error messages in the system is kept as an array of
 * error_code/message pairs, one for major error numbers and another for
 * minor error numbers.
 */
typedef struct H5E_major_mesg_t {
    H5E_major_t error_code;
    const char	*str;
} H5E_major_mesg_t;

typedef struct H5E_minor_mesg_t {
    H5E_minor_t error_code;
    const char	*str;
} H5E_minor_mesg_t;

/* An error stack */
typedef struct H5E_t {
    int	nused;			/*num slots currently used in stack  */
    H5E_error_t slot[H5E_NSLOTS];	/*array of error records	     */
    H5E_auto_t auto_func;       /* Function for 'automatic' error reporting */
    void *auto_data;            /* Callback data for 'automatic' error reporting */
} H5E_t;

/*
 * HERROR macro, used to facilitate error reporting between a FUNC_ENTER()
 * and a FUNC_LEAVE() within a function body.  The arguments are the major
 * error number, the minor error number, and a description of the error.
 */
#define HERROR(maj, min, str) H5E_push(maj, min, FUNC, __FILE__, __LINE__, str)

/*
 * HCOMMON_ERROR macro, used by HDONE_ERROR and HGOTO_ERROR
 * (Shouldn't need to be used outside this header file)
 */
#define HCOMMON_ERROR(maj, min, str)  				              \
   HERROR (maj, min, str);						      \
   (void)H5E_dump_api_stack((int)H5_IS_API(FUNC));

/*
 * HDONE_ERROR macro, used to facilitate error reporting between a
 * FUNC_ENTER() and a FUNC_LEAVE() within a function body, but _AFTER_ the
 * "done:" label.  The arguments are
 * the major error number, the minor error number, a return value, and a
 * description of the error.
 * (This macro can also be used to push an error and set the return value
 *      without jumping to any labels)
 */
#define HDONE_ERROR(maj, min, ret_val, str) {				      \
   HCOMMON_ERROR (maj, min, str);					      \
   ret_value = ret_val;                                                       \
}

/*
 * HGOTO_ERROR macro, used to facilitate error reporting between a
 * FUNC_ENTER() and a FUNC_LEAVE() within a function body.  The arguments are
 * the major error number, the minor error number, the return value, and an
 * error string.  The return value is assigned to a variable `ret_value' and
 * control branches to the `done' label.
 */
#define HGOTO_ERROR(maj, min, ret_val, str) {				      \
   HCOMMON_ERROR (maj, min, str);					      \
   HGOTO_DONE (ret_val)						              \
}

/*
 * HGOTO_DONE macro, used to facilitate normal return between a FUNC_ENTER()
 * and a FUNC_LEAVE() within a function body. The argument is the return
 * value which is assigned to the `ret_value' variable.	 Control branches to
 * the `done' label.
 */
#define HGOTO_DONE(ret_val) {ret_value = ret_val; goto done;}

/* Library-private functions defined in H5E package */
H5_DLL herr_t H5E_push (H5E_major_t maj_num, H5E_minor_t min_num,
			 const char *func_name, const char *file_name,
			 unsigned line, const char *desc);
H5_DLL herr_t H5E_clear (void);
H5_DLL herr_t H5E_dump_api_stack (int is_api);

/*
 * Macros handling system error messages as described in C standard.
 * These macros assume errnum is a valid system error code.
 */

/* Retrieve the error code description string and push it onto the error
 * stack.
 */
#define	HSYS_ERROR(errnum){						      \
    HERROR(H5E_INTERNAL, H5E_SYSERRSTR, HDstrerror(errnum));                  \
}
#define	HSYS_DONE_ERROR(majorcode, minorcode, retcode, str){				      \
    HSYS_ERROR(errno);							      \
    HDONE_ERROR(majorcode, minorcode, retcode, str);			      \
}
#define	HSYS_GOTO_ERROR(majorcode, minorcode, retcode, str){				      \
    HSYS_ERROR(errno);							      \
    HGOTO_ERROR(majorcode, minorcode, retcode, str);			      \
}

#ifdef H5_HAVE_PARALLEL
/*
 * MPI error handling macros.
 */

extern	char	H5E_mpi_error_str[MPI_MAX_ERROR_STRING];
extern	int	H5E_mpi_error_str_len;

#define	HMPI_ERROR(mpierr){						      \
    MPI_Error_string(mpierr, H5E_mpi_error_str, &H5E_mpi_error_str_len);      \
    HERROR(H5E_INTERNAL, H5E_MPIERRSTR, H5E_mpi_error_str);                   \
}
#define	HMPI_DONE_ERROR(retcode, str, mpierr){				      \
    HMPI_ERROR(mpierr);							      \
    HDONE_ERROR(H5E_INTERNAL, H5E_MPI, retcode, str);			      \
}
#define	HMPI_GOTO_ERROR(retcode, str, mpierr){				      \
    HMPI_ERROR(mpierr);							      \
    HGOTO_ERROR(H5E_INTERNAL, H5E_MPI, retcode, str);			      \
}
#endif

#endif
