// C++ informative line for the emacs editor: -*- C++ -*-
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

// Added this line for CC to compile at this time.  Will remove it when
// the problem of "Multiple declaration for RcsId" is fixed. BMR - 10/30/00

// This problem is removed.  I could replace all #include "H5Include.h"
// by #include <hdf5.h>, but decide not to. BMR - 3/22/01

#include <hdf5.h>

// Define bool type for platforms that don't support bool yet
#ifdef BOOL_NOTDEFINED
#ifdef false
#undef false
#endif
#ifdef true
#undef true
#endif
typedef int bool;
const bool  false = 0;
const bool  true  = 1;
#endif


