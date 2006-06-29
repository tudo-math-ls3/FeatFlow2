#							-*- shell-script -*-
#
# Copyright by the Board of Trustees of the University of Illinois.
# All rights reserved.
#
# This file is part of HDF5.  The full HDF5 copyright notice, including
# terms governing use, modification, and redistribution, is contained in
# the files COPYING and Copyright.html.  COPYING can be found at the root
# of the source code distribution tree; Copyright.html can be found at the
# root level of an installed copy of the electronic HDF5 document set and
# is linked from the top-level documents page.  It can also be found at
# http://hdf.ncsa.uiuc.edu/HDF5/doc/Copyright.html.  If you do not have
# access to either file, you may request a copy from hdfhelp@ncsa.uiuc.edu.


# This file is part of the HDF5 build script.  It is processed shortly
# after configure starts and defines, among other things, flags for
# the various compile modes.
#
# See BlankForm in this directory for detailed information.

# The default compiler is `cc'
if test "X-" =  "X-$CC"; then
    CC=cc
    CC_BASENAME=cc
fi

# Try GNU compiler flags.
. $srcdir/config/gnu-flags

# Try native DEC compiler
ARCH=${ARCH:='-arch host -tune host'}
. $srcdir/config/dec-flags
hdf5_mpi_complex_derived_datatype_works=${hdf5_mpi_complex_derived_datatype_works='no'}
