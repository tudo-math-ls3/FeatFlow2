************************************************************************
* Read parameters for time discretization
*
* This routine reads the parameters for the time discretization
* and stores them in the corresponding COMMON blockdiscretization.
*
* In:
*   MDATA  - Unit number to use for reading process
*   MSHOW  - Level of output in the initialization phase
*   CFNAME - Name of the file
************************************************************************

      SUBROUTINE RDTIDC (MDATA,MSHOW,CFNAME)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'ctimediscr.inc'
      
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
        WRITE (*,'(A)') 
     *        'Error in reading time-discretization parameters'
        WRITE (*,'(A)') 'Version number incorrect'
        STOP
      END IF
      
      WRITE (CSTR,'(A)') ' Time discretization parameters'
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      WRITE (CSTR,'(A)') '--------------------------------'
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     Ignore separator
      
      READ (MDATA,*)
      
C     Stationary or nonstationary prolem

      READ (MDATA,*) ISTAT
      IF ((ISTAT.LT.0).OR.(ISTAT.GT.1)) ISTAT=0      
      WRITE (CSTR,1000) 
     *      'Time dependency               : ISTAT  = ', ISTAT
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
C     Separator
        
      READ (MDATA,*)
      WRITE (MTERM,9001)

C     Mesh reconstruction/adaption in time

      READ (MDATA,*) IDMESH
      WRITE (CSTR,'(A,I10)') 
     *      'Mesh reconstruction/adaption:   IDMESH = ', IDMESH
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) IMSCRT
      WRITE (CSTR,'(A,I10)') 
     *      'Mesh adaption criterion:        IMSCRT = ', IMSCRT
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) IMSREP
      WRITE (CSTR,'(A,I10)') 
     *      'Adaption repetitions:           IMSREP = ', IMSREP
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) IMSSTP
      WRITE (CSTR,'(A,I10)') 
     *      'Timesteps between mesh adapt.:  IMSSTP = ', IMSSTP
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) DMSMTM      
      WRITE (CSTR,1001) 
     *      'Max. Time for grid adaption   : DMSMTM = ', DMSMTM
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     Separator
        
      READ (MDATA,*)
      WRITE (MTERM,9001)

C     Stopping criteria on time-dependent solver

      READ (MDATA,*) NITNS
      WRITE (CSTR,'(A,I10)') 
     *      'Max. Number of time steps     : NITNS  = ', NITNS
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) TIMEMX      
      WRITE (CSTR,1001) 
     *      'Max. Time                     : TIMEMX = ', TIMEMX
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) TIMENS
      WRITE (CSTR,1001) 
     *      'Absolute start time           : TIMENS = ', TIMENS
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) TSTEP
      WRITE (CSTR,1001) 
     *      '(Initial) time step           : TSTEP  = ', TSTEP
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) IFRSTP
      WRITE (CSTR,1000) 
     *      'Time stepping scheme          : IFRSTP = ', IFRSTP
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) THETA
      WRITE (CSTR,1001) 
     *      'THETA for one-step scheme     : THETA  = ', THETA
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) IEPSAD
      IF ((IEPSAD.LT.0).OR.(IEPSAD.GT.8)) IEPSAD=1
      WRITE (CSTR,1000) 
     *      'Acceptance criterion          : IEPSAD = ', IEPSAD
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) EPSNS
      WRITE (CSTR,1001) 
     *      'Limit for time derivative     : EPSNS  = ', EPSNS
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     Separator
        
      READ (MDATA,*)
      WRITE (MTERM,9001)

C     Parameters for the adaptive time stepping

      READ (MDATA,*) IADTIM
      IF (ABS(IADTIM).GT.3) IADTIM=0
      WRITE (CSTR,1000) 
     *      'Type of adaptivity            : IADTIM = ', IADTIM
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) IREPIT
      IREPIT = ABS(IREPIT)
      WRITE (CSTR,1000) 
     *      'Max.numbers of repetitions    : IREPIT = ', IREPIT
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) DTMIN
      WRITE (CSTR,1001) 
     *      'Min. Timestep                 : DTMIN  = ', DTMIN
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) DTMAX
      WRITE (CSTR,1001) 
     *      'Max. Timestep                 : DTMAX  = ', DTMAX
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) DTFACT
      DTFACT=ABS(DTFACT)
      IF (DTFACT.LT.1D0) DTFACT=1D0
      WRITE (CSTR,1001) 
     *      'Max. Timestep change          : DTFACT = ', DTFACT
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) EPSADU
      WRITE (CSTR,1001) 
     *      'EPS for not acceptance        : EPSADU = ', EPSADU
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     Separator
        
      READ (MDATA,*)
      WRITE (MTERM,9001)

C     Start procedure

      READ (MDATA,*) IADIN
      IF ((IADIN.LT.0).OR.(IADIN.GT.2)) IADIN=0
      WRITE (CSTR,1000) 
     *      'Type of start procedure       : IADIN  = ', IADIN
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) TIMEIN
      WRITE (CSTR,1001) 
     *      'Time length of start procedure: TIMEIN = ', TIMEIN
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) EPSADI
      WRITE (CSTR,1001) 
     *      'EPS for accep. in start proc. : EPSADI = ', EPSADI
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) EPSADL
      WRITE (CSTR,1001) 
     *      'EPS for accept. after start pr: EPSADL = ', EPSADL
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
      