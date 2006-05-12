************************************************************************
* Read parameters for nonlinear solver of CC2D
*
* This routine reads the parameters for the nonlinear solver used
* in the nonlinear CC2D-iteration from a given file and stores 
* them in the corresponding COMMON blockdiscretization.
*
* In:
*   MDATA  - Unit number to use for reading process
*   MSHOW  - Level of output in the initialization phase
*   CFNAME - Name of the file
************************************************************************

      SUBROUTINE RDNLCC (MDATA,MSHOW,CFNAME)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cnonlinsol_cc2d.inc'
      
C     parameters
      
      CHARACTER CFNAME*(*)
      INTEGER MSHOW,MDATA
      
C     local variables

      INTEGER I,IFMTS
      CHARACTER CSTR*(255)

C     Open DAT file
      
      IFMTS = 1
      CALL OF0 (MDATA,CFNAME,IFMTS)
      
C     Ignore the first 6 lines

      DO I=1,6
        READ (MDATA,*)
      END DO
      
C     Check version number

      READ (MDATA,*) I
      IF (I.NE.100) THEN
        WRITE (*,'(A)') 'Error in reading nonlinear-solver parameters'
        WRITE (*,'(A)') 'Version number incorrect'
        STOP
      END IF
      
      WRITE (CSTR,'(A)') ' Stationary nonlinear solver parameters'
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      WRITE (CSTR,'(A)') '----------------------------------------'
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     Ignore separator
      
      READ (MDATA,*)
      
C     Minimum/Maximum number of iterations
      
      READ (MDATA,*) INLMIN      
      READ (MDATA,*) INLMAX      

      INLMIN = ABS(INLMIN)
      INLMAX = ABS(INLMAX)
      IF (INLMAX.LT.INLMIN) INLMAX=INLMIN
      
      WRITE (CSTR,1000) 
     *      'Min. # of nonlin. iterations  : INLMIN = ', INLMIN
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      WRITE (CSTR,1000) 
     *      'Max. # of nonlin. iterations  : INLMAX = ', INLMAX
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
C     Error bounds, stopping criteria

      READ (MDATA,*) EPSD      
      WRITE (CSTR,1001) 
     *      'Limit for U-defects           : EPSD   = ', EPSD
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) EPSDIV      
      WRITE (CSTR,1001) 
     *      'Limit for DIV-defects         : EPSDIV = ', EPSDIV
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) EPSUR      
      WRITE (CSTR,1001) 
     *      'Limit for U-changes           : EPSUR  = ', EPSUR
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) EPSPR      
      WRITE (CSTR,1001) 
     *      'Limit for P-changes           : EPSPR  = ', EPSPR
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) DMPD      
      WRITE (CSTR,1001) 
     *      'Limit for defect improvement  : DMPD   = ', DMPD
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     Separator
        
      READ (MDATA,*)
      WRITE (MTERM,9001)

C     Parameters for optimal correction

      READ (MDATA,*) OMGMIN      
      READ (MDATA,*) OMGMAX     
      
      IF (OMGINI.LT.OMGMIN) OMGINI=OMGMIN
      IF (OMGINI.GT.OMGMAX) OMGINI=OMGMAX

      WRITE (CSTR,1001) 
     *      'Lower limit for optimal OMEGA : OMGMIN = ', OMGMIN
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      WRITE (CSTR,1001) 
     *      'Upper limit for optimal OMEGA : OMGMAX = ', OMGMAX
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) OMGINI     
      WRITE (CSTR,1001) 
     *      'Start value for optimal OMEGA : OMGINI = ', OMGINI
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     That's it.
      
      WRITE (MTERM,*)
      
      CLOSE (MDATA)
      
1000  FORMAT (A,I4)
1001  FORMAT (A,E16.8)
1002  FORMAT (A,2I4)
1003  FORMAT (A,A)

9000  FORMAT(79('-'))
9001  FORMAT(60('-'))

      END
      