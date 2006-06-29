HDF5 library for storing scientific data in platform independent format.

The original software is available at:   http://hdf.ncsa.uiuc.edu/HDF5/

Due to the fact that this library strongly depends on the hardware, a complicated
configure procedure has to be called before the library can be compiled and linked.
Therefore, this library is not integrated into the standard Feat-style build process.
In contrast, the library is compiled and installed into a temporary directory.
Afterwards, all required files are copied to their Feat-style folders and the
temporary installation is removed.

