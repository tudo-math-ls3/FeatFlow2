***********************************************************************
* This file contains the array names of the B-matrix blocks.
*
* As this file defines the filenames for all levels, the file
* CBASICMG.INC has to be included before!
* These information are only provided as output parameters
* for memory allocation routines.
***********************************************************************

      CHARACTER CARRDB*12
      DIMENSION CARRDB(NNLEV)
      DATA CARRDB/'data/VB.1  ','data/VB.2  ','data/VB.3  ',
     *            'data/VB.4  ','data/VB.5  ','data/VB.6  ',
     *            'data/VB.7  ','data/VB.8  ','data/VB.9  '/
      SAVE CARRDB

