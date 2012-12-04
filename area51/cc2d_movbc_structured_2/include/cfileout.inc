************************************************************************
* This file contains the /FILOUT/ COMMON block. This COMMON block
* saves general information about output to files.
************************************************************************

C     The variable MFILE describes a handle to a textfile where to
C     write standard output to.
C     Is opened at program startup and available to all routines!

      INTEGER MFILE
      
C     Name of the protocol file

      CHARACTER CFLPRT*60
      
      COMMON /FILOUT/ MFILE,CFLPRT
 
      SAVE /FILOUT/
