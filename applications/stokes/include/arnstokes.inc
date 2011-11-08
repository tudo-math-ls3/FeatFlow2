***********************************************************************
* This file contains the array names of the Stokes matrices.
*
* As this file defines the filenames for all levels, the file
* CBASICMG.INC has to be included before!
* These information are only provided as output parameters
* for memory allocation routines.
***********************************************************************

      CHARACTER CARRST*12
      DIMENSION CARRST(NNLEV)
      DATA CARRST/'data/VST.1 ','data/VST.2 ','data/VST.3 ',
     *            'data/VST.4 ','data/VST.5 ','data/VST.6 ',
     *            'data/VST.7 ','data/VST.8 ','data/VST.9 '/
      SAVE CARRST
