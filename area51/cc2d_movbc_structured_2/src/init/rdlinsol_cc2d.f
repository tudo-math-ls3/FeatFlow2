************************************************************************
* Read parameters for linear solver of CC2D
*
* This routine reads the parameters for the linear solver used
* in the nonlinear CC2D-iteration from a given file and stores 
* them in the corresponding COMMON blockdiscretization.
*
* In:
*   MDATA  - Unit number to use for reading process
*   MSHOW  - Level of output in the initialization phase
*   CFNAME - Name of the file
************************************************************************

      SUBROUTINE RDLSCC (MDATA,MSHOW,CFNAME)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'clinsol_cc2d.inc'
      
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
        WRITE (*,'(A)') 'Error in reading linear-solver parameters'
        WRITE (*,'(A)') 'Version number incorrect'
        STOP
      END IF
      
      WRITE (CSTR,'(A)') ' Linear solver parameters'
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      WRITE (CSTR,'(A)') '--------------------------'
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     Ignore separator
      
      READ (MDATA,*)
      
C     Type of linear solver:

      READ (MDATA,*) ILINSL      
      WRITE (CSTR,1000) 
     *      'Type of linear solver         : ILINSL = ', ILINSL
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
C     Minimum/Maximum number of iterations
      
      READ (MDATA,*) ILMIN
      READ (MDATA,*) ILMAX      

      ILMIN=ABS(ILMIN)      
      ILMAX=ABS(ILMAX)
      IF (ILMAX.LT.ILMIN) ILMAX=ILMIN

      WRITE (CSTR,1000) 
     *      'Min. # of steps in lin. sol.  : ILMIN  = ', ILMIN
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      WRITE (CSTR,1000) 
     *      'Max. # of steps in lin. sol.  : ILMAX  = ', ILMAX
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     Number of preconditioning steps if MG is used as preconditioner

      READ (MDATA,*) IPCMG      
      
      WRITE (CSTR,1000) 
     *      '#precond. steps of MG precond.: IPCMG  = ', IPCMG
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
C     Error bounds, stopping criteria

      READ (MDATA,*) DMPMG      
      WRITE (CSTR,1001) 
     *      'Damping of residuals          : DMPMG  = ', DMPMG
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) EPSMG      
      WRITE (CSTR,1001) 
     *      'Limit of residuals            : EPSMG  = ', EPSMG
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     General damping parameter

      READ (MDATA,*) OMGLIN      
      WRITE (CSTR,1001) 
     *      'Damping parameter for solver  : OMGLIN = ', OMGLIN
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
C     Separator
        
      READ (MDATA,*)
      WRITE (MTERM,9001)

C     Multigrid cycle      
      
      READ (MDATA,*) ICYCLE      
      ICYCLE=ABS(ICYCLE)
      WRITE (CSTR,1000) 
     *      'Type of MG cycle              : ICYCLE = ', ICYCLE
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
C     Minimum/Maximum correction parameter in the coarse grid correction
      
      READ (MDATA,*) AMINMG      
      WRITE (CSTR,1001) 
     *      'Lower limit optimal MG-ALPHA  : AMINMG = ', AMINMG
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) AMAXMG      
      WRITE (CSTR,1001) 
     *      'Upper limit optimal MG-ALPHA  : AMAXMG = ', AMAXMG
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
C     Separator
        
      READ (MDATA,*)
      WRITE (MTERM,9001)

C     Parameters for prolongation/restriction/creation of coarse grid
C     matrices

      READ (MDATA,*) IINT      
      IF ((ABS(IINT).LT.1).OR.(ABS(IINT).GT.4)) IINT=2
      WRITE (CSTR,1000) 
     *      'Type of interpolation         : IINT   = ', IINT
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) IAVPR      
      IF ((IAVPR.LT.0)) IAVPR=3
      WRITE (CSTR,1000) 
     *      'Type of averaging in pr/rest  : IAVPR  = ', IAVPR
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) DPREP      
      WRITE (CSTR,1001) 
     *      'Threshold-par. ext. prol/rest.: DPREP  = ', DPREP
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) IAPR      
      IF ((IAVPR.LT.0)) IAPR=3
      WRITE (CSTR,1000) 
     *      'Type of adaptive prol./rest.  : IAPR   = ', IAPR
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     Separator
        
      READ (MDATA,*)
      WRITE (MTERM,9001)

C     Parameters for smoother

      READ (MDATA,*) ISM      
      WRITE (CSTR,1000) 
     *      'Type of smoother              : ISM    = ', ISM
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) NSM      
      NSM=ABS(NSM)
      WRITE (CSTR,1000) 
     *      'Number of smoothing steps     : NSM    = ', NSM
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) RLXSM      
      WRITE (CSTR,1001) 
     *      'Relaxation for the smoother   : RLXSM  = ', RLXSM
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) NSMFAC      
      IF (NSMFAC.LT.1) NSMFAC=1
      WRITE (CSTR,1000) 
     *      'Factor sm. steps on coarser lv: NSMFAC = ', NSMFAC
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
C     Separator
        
      READ (MDATA,*)
      WRITE (MTERM,9001)

C     Coarse grid solver

      READ (MDATA,*) ISL      
      WRITE (CSTR,1000) 
     *      'Type of coarse grid solver    : ISL    = ', ISL
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) NSL      
      NSL=ABS(NSL)
      WRITE (CSTR,1000) 
     *      'Max. # of steps in coarse-g.s.: NSL    = ', NSL
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) RLXSL      
      WRITE (CSTR,1001) 
     *      'Relaxation for the c.g.-solver: RLXSL  = ', RLXSL
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) DMPSL     
      WRITE (CSTR,1001) 
     *      'Damping of residuals in c.g.-s: DMPSL  = ', DMPSL
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) EPSSL      
      WRITE (CSTR,1001) 
     *      'Limit of residuals in c.g.-s. : EPSSL  = ', EPSSL
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
      