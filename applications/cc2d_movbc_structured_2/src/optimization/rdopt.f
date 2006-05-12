***********************************************************************
* This file contains a procedure for reading the preferences
* for the optimisation routines from a DAT file.
***********************************************************************

***********************************************************************
* Read optimisation DAT file
*
* Reads in the file CFNAME from hard disc into optimisation
* COMMON blocks. It will open the file CFNAME, read the parameter and
* close the file.
* Initialises the moving-boundary-geometry if a complex 
* example geometry is prescribed in the DAT-file to be used.
*
* In:
*   MDATA  - IO file handle to use.
*   MSHOW  - Level of output in the initialization phase
*   CFNAME - name of the file
*
***********************************************************************

      SUBROUTINE RDOPT (MDATA,MSHOW,CFNAME)
      
      IMPLICIT NONE
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'
      INCLUDE 'ciniopt.inc'
      
      INCLUDE 'sgeometries.inc'
      INCLUDE 'scompgeometries.inc'
      
C     parameters
      
      CHARACTER CFNAME*(*)
      INTEGER MDATA,MSHOW
      
C     local variables

      INTEGER IFMTS,IBPOST
      CHARACTER CSTR*255
      
      WRITE (CSTR,*) 
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      WRITE (CSTR,'(A)') ' Optimization parameters'
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      WRITE (CSTR,'(A)') '-------------------------'
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
C Like above the file "optimization.dat" specifies the type of simulation.
C
C Line  1: "107"            -> Version 1.0.7
C Line  2: -------------------------------------------------------------
C Line  3: "0"/"1"/"2","3"  -> Optimization type.
C                              0=one test (Parameters line 4-8),
C                              1=brute force (Parameters Line 4,13-15)
C                                for horizontal channel
C                              2=brute force, all X-coords and Y-coords
C                                as configured below (line 16-21)
C Line  4: -------------------------------------------------------------
C Line  5: "-1","0","1"  -> inflow profile; 0=bench, 1=vert. channel,
C                                          -1=user defined
C Line  6: inflow velocity Umax
C Line  7: -------------------------------------------------------------
C Line  8: Type of functional to optimize for
C Line  9: Norm to use in the calculation of the functional
C Line 10: Bitfield to configure the optimization algorithm
C Line 11: Fitness-parameter EPS to reach for the norm of the functional
C Line 12: -------------------------------------------------------------
C Line 13: Minimum number of optimization iterations
C Line 14: Maximum number of optimization iterations
C Line 15: Number of optimization iterations to skip.
C Line 16: Initial stepsize for 1D-based optimization algorithms
C Line 17: minimum X-coordinate 
C Line 18: maximum X-coordinate 
C Line 19: initial X-coordinate
C Line 20: X-stepsize
C Line 21: minimum Y-coordinate
C Line 22: maximum Y-coordinate
C Line 23: initial Y-coordinate
C Line 24: Y-stepsize
C Line 25: minimum rotation angle
C Line 26: maximum rotation angle
C Line 27: initial rotation
C Line 28: rotation angle stepsize
C Line 29: -------------------------------------------------------------
C Line 30: "0"/"1"/"2"/"3"   -> 0=drag/lift normal
C                               1=d/l interpolation without recalc.
C                               2=d/l interpolation with recalc.
C                               3=d/l, minimise vol.eff. before
C Line 31: Relaxation parameter for calculation of ALPHA-vector
C          in volume integration
        
      IFMTS = 1
      CALL  OF0 (MDATA,CFNAME,IFMTS)

C     Version number

      CALL GETINT (MDATA,MSHOW,
     *            'Version number of DAT-file (optim.) '//
     *            'IBPOST = ',IBPOST)
      
      IF (IBPOST.NE.109) THEN
        PRINT *,'Wrong version of optimization.dat!'
        CALL EXIT (1)
      END IF

C     separation line

      READ (MDATA,*)
      
      WRITE (CSTR,9001) 
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
C     Type of simulation/optimization

      CALL GETINT (MDATA,MSHOW,
     *            'Type of simulation/optimization:    '//
     *            'IOPTTP = ',IOPTTP)

      CALL GETINT (MDATA,MSHOW,
     *            'Problem to optimize:                '//
     *            'IOPTPR = ',IOPTPR)

C     separation line

      WRITE (CSTR,9001) 
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*)

C     type of functional

      CALL GETINT (MDATA,MSHOW,
     *            'Type of functional for optimization:'//
     *            'IOPTFL = ',IOPTFL)

C     type of functional

      CALL GETINT (MDATA,MSHOW,
     *            'Type of norm for the functional:    '//
     *            'IOPTNM = ',IOPTNM)

C     configuration bitfield for the opt. algorithm

      CALL GETINT (MDATA,MSHOW,
     *            'Configuration of the opt. algorithm:'//
     *            'IOPTCF = ',IOPTCF)

C     Fitness parameter for the optimizer

      CALL GETDBL (MDATA,MSHOW,
     *            'Optimization fitness parameter EPS: '//
     *            'DOPEPS = ',DOPEPS)

C     separation line

      WRITE (CSTR,9001) 
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*)

C     minimum number of iterations

      CALL GETINT (MDATA,MSHOW,
     *            'Min. no. of opt. iterations:        '//
     *            'OITMIN = ',OITMIN)
      IF (OITMIN.LT.0) OITMIN=0

C     maximum number of iterations

      CALL GETINT (MDATA,MSHOW,
     *            'Max. no. of opt. iterations:        '//
     *            'OITMAX = ',OITMAX)

C     number of iterations to skip

      CALL GETINT (MDATA,MSHOW,
     *            'Number of iterations to skip:       '//
     *            'OITSKP = ',OITSKP)
      IF (OITSKP.LT.0) OITSKP=0

C     Initial step-size for 1D-based optimization algorithms

      CALL GETDBL (MDATA,MSHOW,
     *            'Stepsize-par. for 1D-optimization:  '//
     *            'DOPINS = ',DOPINS)

C     read the X- and Y-bounds/step-sizes 

      CALL GETDBL (MDATA,MSHOW,
     *            'X-minimum:                          '//
     *            'OBXMIN = ',OBXMIN)
      CALL GETDBL (MDATA,MSHOW,
     *            'X-maximum:                          '//
     *            'OBXMAX = ',OBXMAX)
      CALL GETDBL (MDATA,MSHOW,
     *            'initial X (if not single-test):     '//
     *            'OBXINI = ',OBXINI)
      CALL GETDBL (MDATA,MSHOW,
     *            'X-stepsize:                         '//
     *            'OBXSTP = ',OBXSTP)
      CALL GETDBL (MDATA,MSHOW,
     *            'Y-minimum:                          '//
     *            'OBYMIN = ',OBYMIN)
      CALL GETDBL (MDATA,MSHOW,
     *            'Y-maximum:                          '//
     *            'OBYMAX = ',OBYMAX)
      CALL GETDBL (MDATA,MSHOW,
     *            'initial Y (if not single-test):     '//
     *            'OBYINI = ',OBYINI)
      CALL GETDBL (MDATA,MSHOW,
     *            'Y-stepsize:                         '//
     *            'OBYSTP = ',OBYSTP)
      CALL GETDBL (MDATA,MSHOW,
     *            'rot. angle minimum:                 '//
     *            'OBRMIN = ',OBRMIN)
      CALL GETDBL (MDATA,MSHOW,
     *            'rot. angle maximum:                 '//
     *            'OBRMAX = ',OBRMAX)
      CALL GETDBL (MDATA,MSHOW,
     *            'initial rot. (if not single-test):  '//
     *            'OBRINI = ',OBRINI)
      CALL GETDBL (MDATA,MSHOW,
     *            'rot. angle stepsize:                '//
     *            'OBRSTP = ',OBRSTP)

C     separation line

      WRITE (CSTR,9001) 
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*)

C type of drag/lift calculation

      CALL GETINT (MDATA,MSHOW,
     *            'Type of drag-/lift-calculation:     '//
     *            'IOPTCC = ',IOPTCC)

C Relaxation of ALPHA-parameter in drag/lift-calculation

      CALL GETDBL (MDATA,MSHOW,
     *            'ALPHA-rlx-par for drag/lift-calc.:  '//
     *            'DVIARX = ',DVIARX)

      CLOSE (MDATA)
        
9000  FORMAT(79('-'))
9001  FORMAT(60('-'))

      END
      
***********************************************************************
* Initialize optimisation structures
*
* Initializes the optimization structures with standard values.
* If IC2PAR=1, the COMMON-block variables for the optimizer are
* transferred into the structures.
*
*
* In:
*   IC2PAR : =0: initialize IPARAM/DPARAMonly with standard values
*            =1: initialize DPARAM/DPARAM with standard values and
*                transfer values of COMMON-blocks into them
*
* Out:
*   IPARAM  : array [1..SZASMI] of integer
*   DPARAM  : array [1..SZASMD] of double
*             Integer and double prec. parameter block that controls the
*             discretization.
************************************************************************

      SUBROUTINE INIOPT (IPARAM,DPARAM,IC2PAR)
      
      IMPLICIT NONE
      
      INCLUDE 'sinioptimization.inc'
      
      INTEGER IPARAM(SZOPTI),IC2PAR
      DOUBLE PRECISION DPARAM(SZOPTD)

      INCLUDE 'ciniopt.inc'
      INCLUDE 'cinigeometry.inc'
            
C     Clear the structures at first:
            
      CALL LCL3(IPARAM,SZOPTI)
      CALL LCL1(DPARAM,SZOPTD)
      
C     Set nonzero entries to standard values:

      DPARAM(ODOPEPS) = 0.001
      IPARAM(OOITMAX) = -1
      
C     For IC2PAR<>0, copy COMMON-block variables to the structure:
      
      IF (IC2PAR.NE.0) THEN

        IPARAM(OIOPTTP) = IOPTTP
        IPARAM(OIOPTPR) = IOPTPR
        
        IPARAM(OIOPTFL) = IOPTFL
        IPARAM(OIOPTNM) = IOPTNM
        IPARAM(OIOPTCF) = IOPTCF
        DPARAM(ODOPEPS) = DOPEPS
        
        IPARAM(OOITMIN) = OITMIN
        IPARAM(OOITMAX) = OITMAX
        IPARAM(OOITSKP) = OITSKP
        DPARAM(ODOPINS) = DOPINS
        DPARAM(OOBXMIN) = OBXMIN
        DPARAM(OOBXMAX) = OBXMAX
        DPARAM(OOBXINI) = OBXINI
        DPARAM(OOBXSTP) = OBXSTP
        DPARAM(OOBYMIN) = OBYMIN
        DPARAM(OOBYMAX) = OBYMAX
        DPARAM(OOBYINI) = OBYINI
        DPARAM(OOBYSTP) = OBYSTP
        DPARAM(OOBRMIN) = OBRMIN
        DPARAM(OOBRMAX) = OBRMAX
        DPARAM(OOBRINI) = OBRINI
        DPARAM(OOBRSTP) = OBRSTP
        
        IPARAM(OIOPTCC) = IOPTCC
        DPARAM(ODVIARX) = DVIARX       
        
      END IF
            
!      IF (IOPTTP.EQ.1) THEN
!
!C       In this brute-force-test we test all vertical positions between
!C       0.06 and 0.35. and all horizontal radii between 0.005 and 0.1.
!C       Initialise the start position/start size
!
!        DCYPOS = 0.06D0
!        DCRADX = 0.0025D0
!      
!      ELSE IF (IOPTTP.EQ.2) THEN
!        
!C       Brute force test with explicit definition of X-/Y-position 
!C       and step-size
!
!        DCXPOS = OBXMIN
!        DCYPOS = OBYMIN
!        DCRPOS = OBRMIN
!        
!      END IF
      
      END