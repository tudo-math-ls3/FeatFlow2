************************************************************************
* Specify input data for inflow/outflow/...
*
* User specified routine. Prescribed data for files coeff.f and bndry.f.
* With this routine the user can specify input data for the flow
* simulation like inflow, outflow, ...
* The framework repeatly calls this routine with different
* ITYP and IBLOC parameters. Depending on these parameters this
* routine has to return the desired infomation about the flow
* the framework needs.
*
* In:
*  ITYP   - Type of information the framework needs.
*           = 1: velocity dirichlet value
*           = 2: velocity x-derivative
*           = 3: velocity y-derivative
*           = 4: exact pressure
*           = 5: Right hand side for momentum equation
*           = 6: Right hand side for continuity equation
*           = 7: Normal stress in a point (X,Y)
*                =Mean pressure value along edge INODE with
*                 midpoint (X,Y)
*  IBLOC  - Current matrix block; corresponds to the solution component
*           that is currently being assembled by the framework. This
*           must be used by this routine to determine the "direction"
*           of the desired information
*           = 0: matrix block of the pressure
*           = 1: matrix block for x-velocity
*           = 2: matrix block for y-velocity
*  X,Y    - Coordinates of the point where the framework needs
*           the information
*  TIMENS - Current time in instationary Navier-Stokes calculation
*  RE     - Reynolds number; from the DAT file
*  IPARAM - array [1..*] of integer 
*  DPARAM - array [1..*] of integer 
*           TIntAssembly/TDoubleAssembly assembly structures; gives
*           additional information about the discretization.
*  IGEOM  - array [1..*] of integer 
*  DGEOM  - array [1..*] of double 
*           Integer- and double-precision parameter blocks with
*           geometry information. Must be passed to fictitious boundary
*           routines if necessary.
*
* "Optional" parameters:
*  INODE  - Number of the node corresponding to (X,Y). Can be 0.
*           If INODE=0, this function must create the necessary 
*              information only by the given (X,Y)-coordinates of 
*              the point
*           If INODE>0, the parameter TRIA must also be defined.
*              In this case INODE is a node number in the triangulation
*              TRIA with coordinates (X,Y). Then this routine can also
*              access additional information about the triangulation
*              (like precomputed values) to compute the result.
*  TRIA   - array [1..SZTRIA] of integer
*           If INODE>0, TRIA defines the triangulation INODE refers to.
*           If INODE=0, this parameter can be a dummy parameter.
*
* These parameter denote the typical used parameters for FDATIN.
* Any information that is more problem dependent can be accessed
* with IPARAM/DPARAM!
*
* Out:
*  Return value = desired information
************************************************************************

      DOUBLE PRECISION FUNCTION FDATIN(ITYP,IBLOC,X,Y,TIMENS,RE,
     *                          INODE,TRIA,IPARAM,DPARAM,IGEOM,DGEOM)

      IMPLICIT NONE
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'stria.inc'
      
      INCLUDE 'sassembly.inc'
      
C     parameters

      DOUBLE PRECISION X,Y,TIMENS,RE,DPARAM(*),DGEOM(*)
      INTEGER ITYP, IBLOC
      INTEGER INODE,TRIA(SZTRIA),IPARAM(*),IGEOM(*)
      
C     external routines
      
      INTEGER ISFBDY,ISFBDM,NFBDYC
      EXTERNAL ISFBDY,ISFBDM,NFBDYC
      
C     local variables

      INTEGER INPR,LIARR,LDARR
      DOUBLE PRECISION DPAR
      
C     The following function implements a parabolic inflow profile 
C     with maximum velocity = UMAX. LEN denotes the length of the 
C     profile (e.g. length of an edge, depending on the 
C     configuration). The parameter value t in [0,LEN] denotes 
C     the position where the profile should be evaluated.

      DOUBLE PRECISION UMAX,T,LEN,PPROF
      PPROF (T,LEN,UMAX) = 4*UMAX*T*(LEN-T)/(LEN*LEN)
      
C     -------------------------------------------------
C     set the standard result of this function to zero:
C     -------------------------------------------------

      FDATIN=0D0

C --------------------------------------------------------------------
C Now analyse the input parameters and return the corresponding result
C for the current point in the geometry. This has to be specified
C by the user depending in the current geometry and simulation!
C
C In short words: Find out, which DIRICHLET information FEATFLOW 
C wants to know and return it!
C
C The type of information is specified by ITYP, the current point
C by X/Y, the current boundary component on IBLOC, the current time
C of the Navier Stokes simulation in TIMENS,...
C --------------------------------------------------------------------

C     Info: INFLOW CONFIGURATION
C     --------------------------
C     The parameters IINFCF and DPUMAX from the DAT file are passed
C     through the IPARAM/DPARAM parameter blocks to here and can be 
C     accessed if desired!

C     =================================================================
C     *** Case 1: Velocity boundary values and/or exact solution
C     =================================================================

      IF (ITYP.EQ.1) THEN

        IF (IPARAM(OIINFCF).LE.-1) THEN

C         add user defined inflow here!
    
          IF (IBLOC.EQ.1) THEN
C           IBLOC=1: X-velocity
          END IF
          IF (IBLOC.EQ.2) THEN
C           IBLOC=2: Y-velocity
          END IF
        
        ELSE IF (IPARAM(OIINFCF).EQ.0) THEN

C         benchmark-like inflow

          IF (IBLOC.EQ.1) THEN
C           parabolic inflow profile from the left - 
C           horizontal channel (bench1)

C            IF (X.EQ.0.0D0)     FDATIN= 4D0*0.3D0/0.1681D0*Y*(0.41D0-Y)
            IF (X.EQ.0.0D0)     FDATIN= PPROF(Y,0.41D0,DPARAM(ODPUMAX)) 

C            IF (X.EQ.0.0D0)     FDATIN= PPROF(Y,1.0D0,DPARAM(ODPUMAX)) 
C            IF (X.EQ.0.0D0)     FDATIN= PPROF(Y,1D0,DPUMAX) 

          END IF
          
c          IF (IBLOC.EQ.2) THEN
c            IF ((Y.EQ.0D0).AND.(X.LE.0.5D0)) 
c     *          FDATIN=(4.0D0*X*(0.5D0-X))
c          END IF
      
C         Set the velocity on all boundary nodes of all moving boundary
C         objects to zero. Handles to precalculated information are to 
C         be found in the user-defined part of the TRIA structure

          LIARR = TRIA(OTRIUD+2)
          LDARR = TRIA(OTRIUD+3)

          IF (NFBDYC(IGEOM,DGEOM).GT.0) THEN
            IF ((INODE.EQ.0).OR.(LIARR.EQ.0).OR.(LDARR.EQ.0)) THEN
              IF (ISFBDY (X,Y,0,IGEOM,DGEOM).GT.0) THEN
                FDATIN=0D0
              END IF
            ELSE
              IF (ISFBDM (INODE, TRIA, LIARR, LDARR,0,
     *                    IGEOM,DGEOM).GT.0) THEN
                FDATIN=0D0
              END IF
            END IF
          END IF

        ELSE IF (IPARAM(OIINFCF).EQ.1) THEN

C         vertical channel
        
          IF (IBLOC.EQ.1) THEN
          END IF

          IF (IBLOC.EQ.2) THEN

C         IBLOC=2: Y-velocity
C         parabolic inflow profile from the bottom -
C         vertical channel (lqsd12m)

C          IF ((X.EQ.0.0D0).OR.(X.EQ.2D0)) FDATIN= 1D0

            IF (Y.EQ.0D0) THEN
             FDATIN = PPROF(X,2D0,DPARAM(ODPUMAX))
C         FDATIN=4.0D0*10.0D0*
C     *      (X/(2.0-0D0))*(1.0D0-X/(2D0-0D0))
            END IF  
          END IF

C         Set the velocity on all boundary nodes of all moving boundary
C         objects to zero:

          IF (NFBDYC(IGEOM,DGEOM).GT.0) THEN
            IF (ISFBDY (X,Y,0,IGEOM,DGEOM).GT.0) THEN
              FDATIN=0D0
            END IF
          END IF

        END IF
      END IF

C     =================================================================
C     *** Case 2: Velocity x-derivative of exact solution
C     =================================================================

      IF (ITYP.EQ.2) THEN

        IF (IBLOC.EQ.1) THEN
          IF (X.EQ.0.0D0) FDATIN=0D0
        ENDIF

        IF (IBLOC.EQ.2) THEN
          IF (X.EQ.0.0D0) FDATIN=0D0
        ENDIF

      ENDIF

C     =================================================================
C     *** Case 3: Velocity y-derivative of exact solution
C     =================================================================

      IF (ITYP.EQ.3) THEN

        IF (IBLOC.EQ.1) THEN
          FDATIN=0D0
        ENDIF

        IF (IBLOC.EQ.2) THEN
          FDATIN=0D0
        ENDIF

      ENDIF

C     =================================================================
C     *** Case 4: Exact pressure solution
C     =================================================================

      IF (ITYP.EQ.4) THEN
        FDATIN=0D0
      ENDIF

C     =================================================================
C     *** Case 5: Right hand side for momentum equation
C     =================================================================

      IF (ITYP.EQ.5) THEN

        IF (IBLOC.EQ.1) THEN
          FDATIN=0D0
        ENDIF

        IF (IBLOC.EQ.2) THEN
          FDATIN=0D0
        ENDIF

      ENDIF

C     =================================================================
C     *** Case 6: Right hand side for continuity equation
C     =================================================================

      IF (ITYP.EQ.6) THEN
        FDATIN=0D0
      ENDIF

C     =================================================================
C     *** Case 7: Normal stress in a point
C                 = Mean pressure values along the edge
C     =================================================================

      IF (ITYP.EQ.7) THEN
        DPAR=X
        INPR=IBLOC

        IF ((DPAR.GT.1D0).AND.(DPAR.LT.2D0).AND.(INPR.EQ.1)) THEN
          FDATIN=0D0
        ENDIF

        IF ((DPAR.GT.3D0).AND.(DPAR.LT.4D0).AND.(INPR.EQ.1)) THEN
          FDATIN=0D0
        ENDIF

      ENDIF

99999 END

************************************************************************
* Specify Neumann boundary parts
*
* User specified routine. Prescribed data for files coeff.f and bndry.f.
* With this routine the user can specify boundary segments as
* Neumann boundary.
* The framework repeatly calls this routine with different
* parameters. Depending on these parameters this
* routine has to return the desired infomation about the flow
* the framework needs.
*
* In:
*  INPART - Type of information the framework needs.
*           = 0: Number of Neumann boundary segments in the simulation
*           > 0: Information about boundary segment INPART
*  TIMENS - Current time in instationary Navier-Stokes calculation
*  RE     - Reynolds number; from the DAT file
*
* "Optional" parameters:
*  INODE  - Number of the node corresponding to (X,Y). Can be 0.
*           If INODE=0, this function must create the necessary 
*              information only by the given (X,Y)-coordinates of 
*              the point
*           If INODE>0, the parameter TRIA is also be defined.
*              In this case INODE is a node number in the triangulation
*              TRIA with coordinates (X,Y). Then this routine can also
*              access additional information about the triangulation
*              (like precomputed values) to compute the result.
*  TRIA   - array [1..SZTRIA] of integer
*           If INODE>0, TRIA defines the triangulation INODE refers to.
*           If INODE=0, this parameter can be a dummy parameter.
*  IPARAM - array [1..*] of integer 
*  DPARAM - array [1..*] of integer 
*           TIntAssembly/TDoubleAssembly assembly structures; gives
*           additional information about the discretization.
*
* These parameter denote the typical used parameters for FDATIN.
* Any information that is more problem dependent can be accessed
* with IPARAM/DPARAM!
*
* The structure of the parameter blocks depends on the type of the
* solver. The basic structure is of type TSolverIParams/TSolverDParams.
* For stationary simulation, the structure is of type TNSDefIParams/
* TNSDefDParams, while for instationary simulation the structure
* is of type TInsteadDefIParams/TInsteadDefDParams. If necessary, the
* routine can figure out the type using the solver identification
* number IPARAM.SLTAG, which is different depending on the type of
* solver used!
*
* Out:
*  If INPART = 0:
*    INPART - Number of boundary segments
*
*  If INPART > 0:
*    INPRN  - Number of the boundary component, DPAR1/DPAR2 refer to.
*    DPAR1,DPAR2
*           - Parameter values of the boundary segments that should be
*             treated as Neumann boudary
************************************************************************

      SUBROUTINE NEUDAT(INPART,INPRN,DPARN1,DPARN2,TIMENS,
     *                  IPARAM,DPARAM)

      IMPLICIT NONE
      
      INCLUDE 'stria.inc'
      
C     parameters

      INTEGER INPART, INPRN, INODE, TRIA(SZTRIA),IPARAM(*)
      DOUBLE PRECISION DPARN1, DPARN2, TIMENS, DPARAM(*)
      INTEGER IINFCF

C --------------------------------------------------------------------
C Now analyse the input parameters and return the corresponding result
C for the current point in the geometry. This has to be specified
C by the user depending in the current geometry and simulation!
C
C In short words: Find out, which NEUMANN information FEATFLOW 
C wants to know and return it!
C --------------------------------------------------------------------

C     =================================================================
C     *** Case 0: Specify number of Neumann-boundary parts
C     =================================================================

      IINFCF = 0

      IF (INPART.EQ.0) THEN

        IF (IINFCF.LE.-1) THEN
        
C         specify user defined number of Neumann boundary segments here

        ELSE IF (IINFCF.EQ.0) THEN
        
C         benchmark-like configuration; Neumann part on the right side

          INPART=1    
          
        ELSE IF (IINFCF.EQ.1) THEN
        
C         vertical channel; Neumann part on the top side

          INPART=1
          
        END IF

C       =================================================================
C       *** Case <>0: Specify Neumann-boundary parts
C       =================================================================

      ELSE IF (INPART.GT.0) THEN

C       Neumann boundary on outflow

        IF (IINFCF.LE.-1) THEN
        
C         specify user defined Neumann boundary segments here

        ELSE IF (IINFCF.EQ.0) THEN
        
C         benchmark-like configuration

          IF (INPART.EQ.1) THEN
            INPRN =1
            DPARN1=1D0
            DPARN2=2D0
          END IF

          IF (INPART.EQ.2) THEN
            INPRN =1
            DPARN1=3D0
            DPARN2=4D0
          END IF
          
        ELSE IF (IINFCF.EQ.1) THEN
        
C         vertical channel (lqsd12m)

          INPRN =1
          DPARN1=96D0
          DPARN2=120D0
          
        END IF

      ENDIF

99999 END

************************************************************************
* Specify input data for body force and pressure calculation
*
* User specified routine. This routine tells the postprocessing routine
* of any solver whether or not and where to calculate body forces
* (drag/lift) on the real boundary.
* The framework repeatly calls this routine with different
* ITYP parameters. Depending on these parameters this routine has to 
* return the desired infomation what the framework should calculate.
*
* The calculated information is then printed to screen and/or written
* to an external file.
*
* In:
*  ITYP   - Type of information the framework needs.
*           = 1: Number of boundary components where body
*                forces should be evaluated
*           = 2: Boundary component and minimum+maximum parameter value
*                on the boundary between where body forces should
*                be evaluated
*           = 3: Number of fictitious boundary components where to
*                evaluate body forces
*           = 4: Number of the ICNT'th fictitious boundary component
*                where to evaluate body forces
*           = 5: Coefficients in the body force integral
*           = 6: Number of vertices in the grid where the velocity 
*                should be evaluated
*           = 7: Number and/or coordinates of the ICNT'th vertex where 
*                the velocity should be evaluated
*           = 8: Number of vertices in the grid where the pressure 
*                should be evaluated
*           = 9: Number and/or coordinates of the ICNT'th vertex where 
*                the pressure should be evaluated
*           =10: Number and/or coordinates of vertices in the grid where 
*                the streamfunction should be evaluated
*           =11: Number and/or coordinates of the ICNT'th vertex where 
*                the streamfunction should be evaluated
*           =12: Number of boundary components where the integral
*                pressure should be evaluated
*           =13: Boundary component and minimum+maximum parameter value
*                on the boundary between where integral pressure should
*                be evaluated
*  ICNT   - Additional information about for an information the
*           framework needs; depending on ITYP. =0 if not used.
*  TIMENS - Current time in instationary Navier-Stokes calculation
*  RE     - Reynolds number; from the DAT file
*  IPARAM - array [1..*] of integer 
*  DPARAM - array [1..*] of integer 
*           TIntAssembly/TDoubleAssembly assembly structures; gives
*           additional information about the discretization.
*  IGEOM  - array [1..*] of integer 
*  DGEOM  - array [1..*] of double 
*           Integer- and double-precision parameter blocks with
*           geometry information. Must be passed to fictitious boundary
*           routines if necessary.
*
*  TRIA   - array [1..SZTRIA] of integer
*           Current triangulation structure.
*
* Any information that is more problem dependent can be accessed
* with IPARAM/DPARAM!
*
* Out:
*  Depending on ITYP:
*   ITYP= 1: IINFO  = Number of boundary components where body
*                     forces should be evaluated
*       = 2: IINFO  = Boundary component 
*            DINFO1 = minimum parameter value
*            DINFO2 = maximum parameter value
*       = 3: IINFO  = Number of fictitious boundary components where
*                     to evaluate body forces. =0: Don't evaluate.
*                     =-1: Evaluate on all and sum body forces together.
*       = 4: IINFO  = Number of the ICNT'th fictitious boundary
*                     component where to evaluate body forces.
*       = 5: DINFO1,
*            DINFO2 = Coefficients in the body force integral
*       = 6: IINFO  = Number of vertices where the velocity should
*                     be evaluated
*       = 7: IINFO  = Number of the ICNT'th vertex where the velocity
*                     should be evaluated or 0, if the velocity should
*                     be computed in the point at coordinates 
*                     (DINFO1,DINFO2)
*       = 8: IINFO  = Number of vertices where the pressure should be
*                     evaluated
*       = 9: IINFO  = Number of the ICNT'th vertex where the pressure
*                     should be evaluated or 0, if the pressure should
*                     be computed in the point at coordinates 
*                     (DINFO1,DINFO2)
*       =10: IINFO  = Number of vertices where the streamfunction should
*                     be evaluated
*       =11: IINFO  = Number of the ICNT'th vertex where the
*                     streamfunction should be evaluated or 0, if the 
*                     streamfunction should be computed in the 
*                     point at coordinates (DINFO1,DINFO2)
*       =12: IINFO  = Number of boundary components where the integral
*                     pressure should be evaluated
*       =13: IINFO  = Boundary component 
*            DINFO1 = minimum parameter value
*            DINFO2 = maximum parameter value
*
* Other values than the one(s) corresponding to ITYP must not be
* changed! I.e. e.g. if ITYP requests data in IINFO, neighter DINFO1
* nor DINFO2 must be changed!
*
* For evaluating point values (IINFO=7,9,11,13), the IINFO-value can
* be set to 0. In this case, the framework expects (DINFO1,DINFO2) to
* be set to the X/Y coordinates of the point where to evaluate.
* The coordinates must be inside of an element in the domain, otherwise
* the evaluation may fail. The point should be in an element on a rather
* coarse level, otherwise the evaluation takes long for searching
* the point inside of the computational grid! 
************************************************************************

      SUBROUTINE FPTSIN(ITYP,ICNT,TIMENS,RE,
     *                  IPARAM,DPARAM,IGEOM,DGEOM,TRIA,
     *                  IINFO,DINFO1,DINFO2)

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'stria.inc'

      INCLUDE 'sinigeometry.inc'
      
      INCLUDE 'sassembly.inc'
      
C     parameters

      DOUBLE PRECISION TIMENS,RE,DPARAM(*)
      INTEGER ITYP, ICNT
      INTEGER TRIA(SZTRIA),IPARAM(*),IGEOM(*)
      INTEGER IINFO
      DOUBLE PRECISION DINFO1,DINFO2,DGEOM(*)
      
C     externals

      DOUBLE PRECISION TMAX
      EXTERNAL TMAX
      
C     local variables

      DOUBLE PRECISION RHO, DIST, UMEAN, DNU
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      
C     What for information does the framework need?
      
      IF (ITYP.EQ.1) THEN

C       Number of boundary components where to evaluate body forces.

        IINFO = 1

      ELSE IF (ITYP.EQ.2) THEN
      
C       Get information about the boundary components where to
C       evaluate body forces:
C     
C       We evaluate on the complete 2nd boundary component, which
C       normally describes the inner circle of a benchmark like
C       configuration:

        IF (ICNT.EQ.1) THEN
          IINFO = 2
          DINFO1 = 0D0
          DINFO2 = TMAX(IINFO)
        END IF

      ELSE IF (ITYP.EQ.3) THEN
      
C       Number of fictitious boundary components where to evaluate 
C       body forces.
C       Standard handling is: evaluate on all and sum together

        IINFO = -1
      
      ELSE IF (ITYP.EQ.4) THEN
      
C       Number of the ICNT'th fictitious boundary component.
C       Is not used with IINFO=-1. We return 0 here.
C       If IINFO > 0, we have to return a list 1,2,3,...

        IINFO = 0
      
      ELSE IF (ITYP.EQ.5) THEN

C       Coefficients in the body force integral.
C       Parameters for lift (DFW=DINFO1) and drag (DAW=DINFO2).
C       The body force integral is defined as:
C       
C         dfw=2 int_s [dpf1 dut/dn n_y - p n_x] ds / dpf2
C         daw=2 int_s [dpf1 dut/dn n_x + p n_y] ds / dpf2
C
C       Benchmark-like configuration
C       ----------------------------
C       RHO    = density of the fluid; normally normed to 1D0
C       UMEAN  = mean velocity of the parabolic profile
C       DIST   = length of the obstacle that is facing the flow;
C                depending on the direction of the flow!
C
C        RHO  =1.0D0
C        DIST =0.1D0
C        UMEAN=0.2D0

        DNU = 1D0/RE

        RHO   = 1.0D0
        DIST  = DGEOM(ODCRDY)*2D0
        UMEAN = DPARAM(ODPUMAX)*2D0/3D0

        DINFO1 = RHO*DNU
        DINFO2 = RHO*DIST*UMEAN**2

C       Benchmark-case, oscillating cylinder

C        DINFO1 = RHO*DNU
C        DINFO2 = RHO*DIST*(2D0*PI*0.25*0.25)**2

      ELSE IF (ITYP.EQ.6) THEN

C       In how many vertices should the velocity being tracked?

        IINFO = 2

      ELSE IF (ITYP.EQ.7) THEN

C       Get the points where to track the velocity. Write the vertex
C       number to IINFO.
C       These node numbers refer to the standard grid "bench1.tri"!

        IF (ICNT.EQ.1) THEN
          IINFO = 77
        ELSE IF (ICNT.EQ.2) THEN
          IINFO = 298
        END IF
        
C       Alternatively, one can use:
C
C        IINFO = 0
C        IF (ICNT.EQ.1) THEN
C          DINFO1 = 0.15
C          DINFO2 = 0.2
C        ELSE IF (ICNT.EQ.2) THEN
C          DINFO1 = 0.25
C          DINFO2 = 0.2
C        END IF
        
      ELSE IF (ITYP.EQ.8) THEN

C       In how many vertices should the pressure being tracked?

        IINFO = 4

      ELSE IF (ITYP.EQ.9) THEN

C       Get the points where to track the pressure. Write the vertex
C       number to IINFO:

        IF (ICNT.EQ.1) THEN
          IINFO = 1 !77
        ELSE IF (ICNT.EQ.2) THEN
          IINFO = 7 !298
        ELSE IF (ICNT.EQ.3) THEN
          IINFO = 72
        ELSE IF (ICNT.EQ.4) THEN
          IINFO = 74
        END IF

C       Alternatively, one can use:
C
C        IINFO = 0
C        IF (ICNT.EQ.1) THEN
C          DINFO1 = 0.15
C          DINFO2 = 0.2
C        ELSE IF (ICNT.EQ.2) THEN
C          DINFO1 = 0.25
C          DINFO2 = 0.2
C        ELSE IF (ICNT.EQ.3) THEN
C          DINFO1 = 0.2
C          DINFO2 = 0.15
C        ELSE IF (ICNT.EQ.4) THEN
C          DINFO1 = 0.2
C          DINFO2 = 0.25
C        END IF

      ELSE IF (ITYP.EQ.10) THEN

C       In how many vertices should the streamfunction being tracked?

        IINFO = 2

      ELSE IF (ITYP.EQ.11) THEN

C       Get the points where to track the streamfunction. Write the
C       vertex number to IINFO:

        IF (ICNT.EQ.1) THEN
          IINFO = 77
        ELSE IF (ICNT.EQ.2) THEN
          IINFO = 298
        END IF

C       Alternatively, one can use:
C
C        IINFO = 0
C        IF (ICNT.EQ.1) THEN
C          DINFO1 = 0.15
C          DINFO2 = 0.2
C        ELSE IF (ICNT.EQ.2) THEN
C          DINFO1 = 0.25
C          DINFO2 = 0.2
C        END IF

      ELSE IF (ITYP.EQ.12) THEN

C       Number of boundary components where to evaluate the pressure
C       integral

        IINFO = 2

      ELSE IF (ITYP.EQ.13) THEN

C       Get information about the boundary components where to
C       evaluate the pressure integral.
C     
C       We evaluate on one hand on the complete 2nd boundary component,
C       on the other hand on the half second boundary component.
C       In benchmark comfiguration, this refers on one hand to the
C       full inner circle and on the other hand to the face of the
C       circle facing the stream.

        IF (ICNT.EQ.1) THEN
          IINFO = 2
          DINFO1 = 0D0
          DINFO2 = TMAX(IINFO)
        ELSE IF (ICNT.EQ.2) THEN
          IINFO = 2
          DINFO1 = 0D0
          DINFO2 = 0.5D0*TMAX(IINFO)
        END IF

      END IF
      
      END
      