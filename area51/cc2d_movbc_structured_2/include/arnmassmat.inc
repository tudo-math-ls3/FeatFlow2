***********************************************************************
* This file contains the array names of the mass matrices.
*
* As this file defines the filenames for all levels, the file
* CBASICMG.INC has to be included before!
* These information are only provided as output parameters
* for memory allocation routines.
***********************************************************************

      CHARACTER CARRM*12
      DIMENSION CARRM(NNLEV)

      DATA CARRM /'data/VM.1  ','data/VM.2  ','data/VM.3  ',
     *            'data/VM.4  ','data/VM.5  ','data/VM.6  ',
     *            'data/VM.7  ','data/VM.8  ','data/VM.9  '/

      SAVE CARRM
      