!************************************************************************
!* This file contains the /FILOUT/ COMMON block. This COMMON block
!* saves general information about output to files.
!************************************************************************

!     The variable MFILE describes a handle to a textfile where to
!     write standard output to.
!     Is opened at program startup and available to all routines!

      INTEGER :: MFILE
      
!     Name of the protocol file

      CHARACTER CFLPRT*60
      
      COMMON /FILOUT/ MFILE,CFLPRT
 
      SAVE /FILOUT/
