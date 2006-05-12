************************************************************************
* This file contains the standard pre-/postprocessing routines for the
* nonstationary solver.
* The user can write its own postprocessing routine and call this
* old one to perform standard postprocessing.
************************************************************************

************************************************************************
* Default Pre- and postprocessing for nonstationary solver
*
* This is the default callback routine for pre- and postprocessing
* in the nonstationary Navier-Stokes solver NONST. It's called at
* different positions inside of the algorithm to do something user-
* defined with what is calculated.
*
* The central parameter in this routine is IFINL. This variable
* gives an identifier for the position where the solver called this
* routine. Then, depending on IFINL, the routine can decide what to do.

* This subroutine performs basic postprocessing. It has nearly
* all information available for that task, the nonstationary
* solver is working with, but should avoid changing any information
* in these arrays. Nevertheless, activating a MG solver for any
* kind of postprocessing is of course fine, e.g., as long as the
* configuration is not changed...
*
* In:
*   NLMIN  : minimum level 
*   NLMAX  : maximum level 
*   TRIAS  : array [1..SZTRIA,1..NLEV] of integer
*            Triangulation structures for all levels.
*   MATDAT : array [1..SZN2MI,1..NNLEV] of integer
*            TNS2DMatrixParams-structures for level NLMIN..NLMAX, 
*            initialized with data
*   VECDAT : array [1..SZN2VI,1..NNLEV] of integer
*            TNS2DVectorParams-structures for level NLMIN..NLMAX. 
*            This structure array must specify the structure of
*            the vectors on each level. 
*
*   IPARAM : array [1..SZISDI] of integer
*   DPARAM : array [1..SZISDI] of double
*            Integer and double prec. parameter blocks that define the
*            behaviour of the nonstationary solver. 
*   ISTPAR : array [1..SZNSDI] of integer
*   DSTPAR : array [1..SZNSDD] of double
*            Integer and double prec. parameter block for the stationary 
*            sub-solver NSDEF2.
*   IMGPAR : array [1..SZ020I+3*SZSLVI] of integer
*   DMGPAR : array [1..SZ020D+3*SZSLVD] of double
*            Integer and double parameter blocks for the multigrid 
*            sub-solver M020. 
*   IASMBL : array [1..SZASMI] of integer
*   DASMBL : array [1..SZASMD] of double
*            Integer and double prec. parameter block that controls the
*            discretization. 
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to fictitious boundary
*            routines. Not used in this routine.
*   IADTS  : array [1..SZADTI] of integer
*   DADTS  : array [1..SZADTD] of double
*            Integer and double prec. parameter block that configures
*            the behaviour of the adaptive time stepping.
*   IPPDAT : array [1..SZISPI] of integer
*   DPPDAT : array [1..SZISPD] of double
*            TInstatPostprocessingIParams and TInstatPostprocessingDParams
*            structures that define the behaviour of the postprocessing
*            routines; contains input and status variables for DFPSTD.
*   ITRACK : array {1..SZTRKI] of integer
*            TTrackingIParams structure that holds the output handles
*            for writing tracked solution values to files.
*
*   NUVP   : Total length of solution vector
*   DUP    : array [1..NUVP] of double
*            Current solution vector
*   DRHS   : array [1..NUVP] of double
*            Right hand side for the next time step
*   DAUX   : array [1..NUVP] of double
*            Auxiliary vector.
*
*   TIMEST - Initial simulation time
*   TIMENS - Current simulation time
*   TSTEP  - Current time step size. Can vary from time-step to
*            time-step because of the Fractional-Step scheme
*            or the adaptive time stepping.
*   ICORST - Number of substep in the time stepping scheme
*   NCORST - Total number of substeps in the time stepping scheme
*   ITNSR  - Count on how often the current time step has been repeated
*            =0: first calculation, =1: first repetition,...
*   STSTP1 - Status of nonlinear solver of prediction step.
*            Might be also 0 if not used.
*   STSTP2 - Status of nonlinear solver of correction step.
*            =0, if DUP is the solution of the predictor step.
*   IFINL  - Identifier that informes the routine about when it was
*            called and for what purpose.
*            =-30: The routine is called at the very beginning.
*                  Time-stepping is initialized and memory for the
*                  vectors is allocated. Yet, RHS is not initialized
*                  and the computational grid is not prepared. Initial
*                  preparations for the first time step can be made.
*            =-20: The routine is called before the time stepping 
*                  starts. Matrices/vectors/... are prepared for the
*                  first time step.
*            =-10: The routine is called before a time step starts.
*                  Matrices/vectors/... are prepared for the next time
*                  step.
*            = -2: The routine is called before the predictor step;
*                  might perform some preprocessing.
*            =  2: DUP is the solution of the predictor step,
*                  which will be used for the calculation of the
*                  adaptive time stepping. The routine might perform
*                  some postprocessing.
*            = -1: The routine is called before a corrector step;
*                  might perform some preprocessing.
*            =  1: DUP is the solution of a corrector step which might 
*                  be rejected if the error analysis decides to repeat
*                  the time step.
*                  The routine might perform some postprocessing.
*            =  0: DUP is a confirmed solution of the corrector
*                  step and will definitely be used for further
*                  computation.
*                  The routine might perform some postprocessing.
*            =  5: DUP is rejected. The simulation time will be reset
*                  to the last accepted solution, calculation will
*                  restart from there. 
*                  The routine might perform some cleanup/reset.
*            = 20: The simulation is finished, DUP contains the 
*                  solution of the last performed time step.
*   ISTATN - Status for the previous run; might be modified by PREPOS.
*
* Out:
*   ISTATN - New status for the next run. Defines what to do in the
*            preparation routine of the next time step. Can be used
*            to tell the solver to skip some calculation in the
*            preparation for the next time step. Bitfield.
*            The solver does not touch ISTATN. At the beginning,
*            all bits of ISTATN are set. If PREPOS modifies ISTATN,
*            it must reset it on its own.
*
*            Bit0: Perform mesh adaption in the preparation of the
*                  next time step as configured by the parameters.
*   IRPSTP - Repeat time step.
*            Can be set to 1 by PREPOS when IFINL=1 or =2.
*            This forces the solver to repeat the time step.
*            The time step will only be repeated if adaptive time
*            stepping is additionally activated.
*
* Remarks: 
* 1.) Multiple calls for the same solution
* ----------------------------------------
*  The pre-/postprocessing routine is often called twice or
*  more often for one and the same solution. On one hand the routine 
*  is called after a new solution is calculated - but at this stage the 
*  solution is not confirmed, as adaptive time stepping might reject it 
*  and force it to be recalculated. So at this stage, the routine is 
*  called with IFINL=1.
*  In a later step when a solution is confirmed, the routine is again
*  called for this solution, this time with IFINL=0.
*  Note that depending on the time stepping scheme, not every solution
*  can be a confirmed one. When using adaptive time stepping, only every
*  3rd step or so (more exactly, every NCORST'st step) a comparison
*  with a predictor step happens, so only every 3rd solution can be
*  a confirmed one!
*
* 2.) Predictor-/Corrector steps
* ------------------------------
*  This routine is called for the solution of the predictor step as well
*  as for standard time steps, which are also called "corrector steps"
*  in this sense.
*  Of course, the predicted solution is intermediate and is only
*  used for the calculation of the time step. So for this solution,
*  PREPOS is only called with IFINL=2.
*
* 3.) Pre-/Postprocessing
* -----------------------
*  The routine has two duties: Pre- and postprocessing.
*  For this purpose, it's called twice - before and after each time
*  step.
*
*  When the routine is called before a time step, there is IFINL<0.
*  In this case the routine can to preparations for the next time
*  step at time TIMENS+TSTEP. The time will be increased to that
*  during the next time step. For the preprocessing before the time
*  stepping, TSTEP will be 0, so TIMENS+TSTEP is still the correct
*  formula for calculating the time of the very first time step.
*
*  When the routine is called after a time step, there is IFINL>=0.
*  In this case the routine can perform any postprocessing for the
*  solution that is computed at time TIMENS.
*
* 4.) Feedback to the solver
* --------------------------
*  The variable ISTATN allowes to give a feedback from the
*  pre-/postprocessing routine to the solver to skip some calculation in
*  the preparation of the next time step. E.g. it allowes to skip
*  the mesh adaption, which is typically not necessary if the geometry
*  does not move. If ISTATN is not modified by PREPOS, the Navier
*  Stokes solver carries out its standard behaviour.
*
*  The vectors identified by VECDAT can be arbitrarily used, since
*  when postprocessing is called, no MG solver is active.
************************************************************************
      
      SUBROUTINE DFPSTD (NLMIN,NLMAX,
     *                   TRIAS,MATDAT,VECDAT,
     *                   IPARAM,DPARAM,ISTPAR,DSTPAR,
     *                   IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *                   IADTS,DADTS,IPPDAT,DPPDAT,ITRACK,
     *                   NUVP,DUP,DRHS,DAUX,
     *                   TIMEST, TIMENS, TSTEP, ICORST, NCORST,
     *                   ITNSR,STSTP1,STSTP2,IFINL,ISTATN,IRPSTP)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'
      
      INCLUDE 'stria.inc'
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'smat2dns.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      
      INCLUDE 'ssolvers.inc'
      INCLUDE 'stiming.inc'
      INCLUDE 'snonsteadpostproc.inc'
      
      INCLUDE 'snonstead.inc'
      
      INCLUDE 'sassembly.inc'
      
      INCLUDE 'stracking.inc'
      
C parameters

      INTEGER NUVP
      INTEGER IPARAM(*),IMGPAR(*),ISTPAR(*),IASMBL(*),IADTS(*),IGEOM(*)
      INTEGER ITRACK(*)
      INTEGER TRIAS(SZTRIA,NNLEV),IPPDAT(SZISPI)
      INTEGER MATDAT(SZN2MI,NNLEV),VECDAT(SZN2VI,NNLEV),NLMIN,NLMAX
      DOUBLE PRECISION DPARAM (*),DMGPAR(*),DSTPAR(*),DASMBL(*),DADTS(*)
      DOUBLE PRECISION DGEOM(*),DPPDAT(SZISPD)
      DOUBLE PRECISION TIMEST, TIMENS, TSTEP
      
      DOUBLE PRECISION DUP(*),DRHS(*),DAUX(*)
      INTEGER ICORST,NCORST,ITNSR,STSTP1,STSTP2
      INTEGER IFINL,ISTATN,IRPSTP
      
C     Coefficient of exact solution

      DOUBLE PRECISION UE,PE,UEX,UEY
      EXTERNAL UE,PE,UEX,UEY
      
C     finite elements

      EXTERNAL E030,E031,EM30,EM31,E010,E011
      
C     externals for GMV output

      EXTERNAL DFGMMC,DFGMMV
      
C     maximal number of tracked solutions
      
      INTEGER MAXTRK
      PARAMETER (MAXTRK=32)
      
C     local variables

      INTEGER I,J,K,II,JJ,MSHOW,NEQU,NEQP
      INTEGER LU,LV,LP,LISO,IEL,LERR
      INTEGER KCORVG,KXNPR
      INTEGER KVERT,KMID,KVBD,KEBD,KVBDP,KBCT,KCORMG,NCELLS,NVERTS
      INTEGER IINFO(MAXTRK)
      DOUBLE PRECISION DINFO1(MAXTRK),DINFO2(MAXTRK),DINFO3(MAXTRK)
      DOUBLE PRECISION D1,D2,COEF1,COEF2,DTMP1,DTMP2
      CHARACTER CFN*60

C     Ok, this routine seems to get a little bit stretchy now, but this
C     lies in the nature of postprocessing, i.e. handling of all the 
C     cases :)

      MSHOW = IPARAM(OMSGTRM)
      
C     Set LU,LV,LP,LISO to 0. These are handles to vertex-based solution
C     vectors which might be calculated during the postprocessing if
C     needed; in this case we have to release them later.

      LU = 0
      LV = 0
      LP = 0
      LISO = 0
      LERR = 0

C     Get the vector size of velocity and pressure part on the finest
C     level; they are frequently used.

      NEQU = VECDAT(ONU,NLMAX)
      NEQP = VECDAT(ONP,NLMAX)

C     =================================================================
C     Calculation of the body forces.
C     Output of tracked solution values.

C     On the beginning of the time stepping, reset the number of
C     solutions we have written out:

      IF (IFINL.EQ.-20) THEN
        IPPDAT(OTRKCNT) = 0
      END IF

C     We calculate the forces and output tracked solution values
C     only on accepted solutions!

      IF (IFINL.EQ.0) THEN 
      
C       Fetch some variables from the geometry for less stuff to 
C       write :) We are working here onb the maximum level.

        KVERT  = L(TRIAS(OLVERT,NLMAX))
        KMID   = L(TRIAS(OLMID,NLMAX))
        KVBD   = L(TRIAS(OLVBD,NLMAX))
        KEBD   = L(TRIAS(OLEBD,NLMAX))
        KCORVG = L(TRIAS(OLCORVG,NLMAX))
        KCORMG = L(TRIAS(OLCORMG,NLMAX))
        KVBDP  = L(TRIAS(OLVBDP,NLMAX))
        KBCT   = L(TRIAS(OLBCT,NLMAX))
        KXNPR  = L(TRIAS(OLXNPR,NLMAX))

C       --------------------------------------------------------------      
C       At first the standard boundary forces evaluation on real 
C       boundary components.
C       Is there a cubature formula for the "real boundary evaluation"
C       given?

        IF (IPPDAT(OIBDFBD).GT.0) THEN
        
C         In how many boundary components should the body force be
C         evaluated? We store that number in J.

          CALL FPTSIN(1,0,TIMENS,DASMBL(ORE),
     *                IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                J,D1,D2)
     
C         We only track a limited number of points

          J = MIN(J,MAXTRK)
     
C         Get the coefficients of the boundary integral into COEF1
C         and COEF2:

          CALL FPTSIN(5,0,TIMENS,DASMBL(ORE),
     *                IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                0,COEF1,COEF2)
     
C         Get the boundary component segments where to evaluate the
C         forces. For now we restrict to a maximum of MAXTRK segments.
C         Evaluate there!

          DO I=1,J
          
C           Request data about the boundary component.
C           Its number is stored in K and the parameter values in D1/D2.

            CALL FPTSIN(2,I,TIMENS,DASMBL(ORE),
     *                  IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                  K,D1,D2)
     
            IINFO(I) = K
     
            IF (IASMBL(OIELEMT).EQ.0) THEN
            
              CALL BDFORX(DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *                TRIAS(1,NLMAX),KWORK(KVERT),
     *                KWORK(KMID),KWORK(KVBD),KWORK(KEBD),KWORK(KBCT),
     *                DWORK(KCORVG),DWORK(KVBDP),E031,.FALSE.,
     *                K,D1,D2,COEF1,COEF2,
     *                1,DINFO1(I),DINFO2(I),DINFO3(I))   
     
            ELSE IF (IASMBL(OIELEMT).EQ.1) THEN
            
              CALL BDFORX(DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *                TRIAS(1,NLMAX),KWORK(KVERT),
     *                KWORK(KMID),KWORK(KVBD),KWORK(KEBD),KWORK(KBCT),
     *                DWORK(KCORVG),DWORK(KVBDP),E030,.FALSE.,
     *                K,D1,D2,COEF1,COEF2,
     *                1,DINFO1(I),DINFO2(I),DINFO3(I))   
     
            ELSE IF (IASMBL(OIELEMT).EQ.2) THEN
            
              CALL BDFORX(DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *                TRIAS(1,NLMAX),KWORK(KVERT),
     *                KWORK(KMID),KWORK(KVBD),KWORK(KEBD),KWORK(KBCT),
     *                DWORK(KCORVG),DWORK(KVBDP),EM31,.TRUE.,
     *                K,D1,D2,COEF1,COEF2,
     *                1,DINFO1(I),DINFO2(I),DINFO3(I))   
     
            ELSE IF (IASMBL(OIELEMT).EQ.3) THEN
            
              CALL BDFORX(DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *                TRIAS(1,NLMAX),KWORK(KVERT),
     *                KWORK(KMID),KWORK(KVBD),KWORK(KEBD),KWORK(KBCT),
     *                DWORK(KCORVG),DWORK(KVBDP),EM30,.TRUE.,
     *                K,D1,D2,COEF1,COEF2,
     *                1,DINFO1(I),DINFO2(I),DINFO3(I))  
     
            END IF 
          
          END DO
          
C         Print the results to screen/file:

          CALL TV2TER (MSHOW,MFILE,J,IINFO,DINFO1,DINFO2,2,
     *      'Body forces real bd., constant pressure, bdc/horiz/vert')
        
C         Write to tracking file, if handles are given.
C         Write horizontal and vertical force.

          IF (ITRACK(OHBFB1).GT.0) THEN
            CALL TV2FIL (ITRACK(OHBFB1),J,IINFO,DINFO1,
     *                   IPPDAT(OTRKCNT).EQ.0,TIMENS,-1,'BDC')
          END IF

          IF (ITRACK(OHBFB2).GT.0) THEN
            CALL TV2FIL (ITRACK(OHBFB2),J,IINFO,DINFO2,
     *                   IPPDAT(OTRKCNT).EQ.0,TIMENS,-1,'BDC')
          END IF
        
        END IF ! IBDFBD <> 0

C       --------------------------------------------------------------      
C       Boundary forces evaluation on fictitious boundary 
C       by volume integration.
C       Is there a cubature formula for the "real boundary evaluation"
C       given?

        IF (IPPDAT(OIBDFVI).GT.0) THEN
        
C         In how many boundary components should the body force be
C         evaluated? We store that number in J.

          CALL FPTSIN(3,0,TIMENS,DASMBL(ORE),
     *                IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                J,D1,D2)
     
C         We only track a limited number of points

          J = MIN(J,MAXTRK)
          
          IF (J.NE.0) THEN
          
C           Get the coefficients of the boundary integral into COEF1
C           and COEF2:

            CALL FPTSIN(5,0,TIMENS,DASMBL(ORE),
     *                  IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                  0,COEF1,COEF2)
     
C           Get the boundary component segments where to evaluate the
C           forces. For now we restrict to a maximum of MAXTRK segments.
C           Evaluate there!
C           If J is -1, we have to evaluate only once - for all fictitious
C           boundary components at the same time!

            DO I=1,MAX(1,J)
            
C             Request data about the boundary component.
C             The number of the fict. bdry component is stored to K.
C             If J=-1, we set K=0 to evaluate everywhere.

              IF (J.EQ.-1) THEN
                K = 0
              ELSE
                CALL FPTSIN(4,I,TIMENS,DASMBL(ORE),
     *                    IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                    K,D1,D2)
              END IF
              
C             Save the boundary component number to IINFO:              
              
              IINFO(I) = K
     
C             Call the evaluation routine, calculate forces on that
C             boundary component
     
              IF (IASMBL(OIELEMT).EQ.0) THEN
              
                CALL BDFVLG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),NEQU,
     *                   DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *                   KWORK(KMID),KWORK(KXNPR),TRIAS(1,NLMAX),
     *                   E031,.FALSE.,DINFO1(I),DINFO2(I),COEF1,COEF2,
     *                   K,0.5D0,IGEOM,DGEOM)
              
              ELSE IF (IASMBL(OIELEMT).EQ.1) THEN
              
                CALL BDFVLG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),NEQU,
     *                   DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *                   KWORK(KMID),KWORK(KXNPR),TRIAS(1,NLMAX),
     *                   E030,.FALSE.,DINFO1(I),DINFO2(I),COEF1,COEF2,
     *                   K,0.5D0,IGEOM,DGEOM)

              ELSE IF (IASMBL(OIELEMT).EQ.2) THEN
              
                CALL BDFVLG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),NEQU,
     *                   DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *                   KWORK(KMID),KWORK(KXNPR),TRIAS(1,NLMAX),
     *                   EM31,.TRUE.,DINFO1(I),DINFO2(I),COEF1,COEF2,
     *                   K,0.5D0,IGEOM,DGEOM)

              ELSE IF (IASMBL(OIELEMT).EQ.3) THEN
              
                CALL BDFVLG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),NEQU,
     *                   DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *                   KWORK(KMID),KWORK(KXNPR),TRIAS(1,NLMAX),
     *                   EM30,.TRUE.,DINFO1(I),DINFO2(I),COEF1,COEF2,
     *                   K,0.5D0,IGEOM,DGEOM)

              END IF 
            
            END DO
            
C           Print the results to screen/file:

            CALL TV2TER (MSHOW,MFILE,MAX(J,1),IINFO,DINFO1,DINFO2,2,
     *        'Body forces fict.bd., constant pressure, bdc/horiz/vert')
        
C           Write to tracking file, if handles are given.
C           Write horizontal and vertical force.

            IF (ITRACK(OHBFV1).GT.0) THEN
              CALL TV2FIL (ITRACK(OHBFV1),MAX(J,1),IINFO,DINFO1,
     *                     IPPDAT(OTRKCNT).EQ.0,TIMENS,-1,'BDC')
            END IF

            IF (ITRACK(OHBFV2).GT.0) THEN
              CALL TV2FIL (ITRACK(OHBFV2),MAX(J,1),IINFO,DINFO2,
     *                     IPPDAT(OTRKCNT).EQ.0,TIMENS,-1,'BDC')
            END IF

          END IF
        
        END IF ! IBDFVL <> 0

C       --------------------------------------------------------------      
C       Boundary forces evaluation on fictitious boundary 
C       by line integration.
C       Is there a cubature formula for the "real boundary evaluation"
C       given?

        IF (IPPDAT(OIBDFLI).GT.0) THEN
        
C         In how many boundary components should the body force be
C         evaluated? We store that number in J.

          CALL FPTSIN(3,0,TIMENS,DASMBL(ORE),
     *                IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                J,D1,D2)
     
C         We only track a limited number of points

          J = MIN(J,MAXTRK)
          
          IF (J.NE.0) THEN
          
C           Get the coefficients of the boundary integral into COEF1
C           and COEF2:

            CALL FPTSIN(5,0,TIMENS,DASMBL(ORE),
     *                  IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                  0,COEF1,COEF2)
     
C           Get the boundary component segments where to evaluate the
C           forces. For now we restrict to a maximum of MAXTRK segments.
C           Evaluate there!
C           If J is -1, we have to evaluate only once - for all fictitious
C           boundary components at the same time!

            DO I=1,MAX(1,J)
            
C             Request data about the boundary component.
C             The number of the fict. bdry component is stored to K.
C             If J=-1, we set K=0 to evaluate everywhere.

              IF (J.EQ.-1) THEN
                K = 0
              ELSE
                CALL FPTSIN(4,I,TIMENS,DASMBL(ORE),
     *                    IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                    K,D1,D2)
              END IF
     
C             Save the boundary component number to IINFO:              
              
              IINFO(I) = K
     
C             Call the evaluation routine, calculate forces on that
C             boundary component.
C             Pass the assembly data structures as INFO-block to
C             these routines. That way all called subroutines
C             can access information about the current assembly
C             status.
     
              IF (IASMBL(OIELEMT).EQ.0) THEN
              
                CALL BDFRIG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *                       DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *                       KWORK(KMID),TRIAS(1,NLMAX),
     *                       E031,.FALSE.,DINFO1(I),DINFO2(I),COEF1,
     *                       COEF2,K,IGEOM,DGEOM)
              
              ELSE IF (IASMBL(OIELEMT).EQ.1) THEN
              
                CALL BDFRIG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *                       DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *                       KWORK(KMID),TRIAS(1,NLMAX),
     *                       E030,.FALSE.,DINFO1(I),DINFO2(I),COEF1,
     *                       COEF2,K,IGEOM,DGEOM)

              ELSE IF (IASMBL(OIELEMT).EQ.2) THEN
              
                CALL BDFRIG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *                       DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *                       KWORK(KMID),TRIAS(1,NLMAX),
     *                       EM31,.TRUE.,DINFO1(I),DINFO2(I),COEF1,
     *                       COEF2,K,IGEOM,DGEOM)

              ELSE IF (IASMBL(OIELEMT).EQ.3) THEN
              
                CALL BDFRIG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *                       DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *                       KWORK(KMID),TRIAS(1,NLMAX),
     *                       EM30,.TRUE.,DINFO1(I),DINFO2(I),COEF1,
     *                       COEF2,K,IGEOM,DGEOM)

              END IF 
            
            END DO
            
C           Print the results to screen/file:

            CALL TV2TER (MSHOW,MFILE,MAX(J,1),IINFO,DINFO1,DINFO2,2,
     *        'Body forces f.b.l.i., constant pressure, bdc/horiz/vert')
        
C           Write to tracking file, if handles are given.
C           Write horizontal and vertical force.

            IF (ITRACK(OHBFL1).GT.0) THEN
              CALL TV2FIL (ITRACK(OHBFL1),MAX(J,1),IINFO,DINFO1,
     *                     IPPDAT(OTRKCNT).EQ.0,TIMENS,-1,'BDC')
            END IF

            IF (ITRACK(OHBFL2).GT.0) THEN
              CALL TV2FIL (ITRACK(OHBFL2),MAX(J,1),IINFO,DINFO2,
     *                     IPPDAT(OTRKCNT).EQ.0,TIMENS,-1,'BDC')
            END IF
        
          END IF
        
        END IF ! IBDFVLI <> 0
        
C       --------------------------------------------------------------      
C       Output of tracked solution values.
C
C       At first: Is the output of tracked solution values active?
C       Check if we have a file handle where to output these values
C       to.
C
C       ----------
C       Velocities

        IF ((ITRACK(OHVEL1).GT.0).OR.(ITRACK(OHVEL2).GT.0)) THEN
        
C         Calculate the interpolated velocity, if it's not calculated

          IF ((LU.EQ.0).OR.(LV.EQ.0)) THEN
          
            CALL XINTUV (DUP(1),DUP(1+NEQU),TRIAS(1,NLMAX),LU,LV)
          
C           Implement boundary conditions, since the interpolation
C           does not correctly handle boundaries:

            KCORVG = L(TRIAS(OLCORVG,NLMAX))
            KXNPR  = L(TRIAS(OLXNPR,NLMAX))
            CALL BDRCOR (DWORK(L(LU)),DWORK(L(LV)),
     *                   TRIAS(1,NLMAX),DWORK(KCORVG),
     *                   KWORK(KXNPR),UE,TIMENS,DASMBL(ORE),
     *                   IASMBL,DASMBL,IGEOM,DGEOM)
          END IF

C         Ask the user defined routine in how many points we should
C         track the velocity; write the number into K.

          CALL FPTSIN(6,0,TIMENS,DASMBL(ORE),
     *              IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *              J,D1,D2)
     
C         We only track a limited number of points

          J = MIN(J,MAXTRK)
          
          IF (J.GT.0) THEN
          
C           Get the point numbers and the values to
C           IINFO/DINFO1/DINFO2:

            DO I=1,J
              CALL FPTSIN(7,I,TIMENS,DASMBL(ORE),
     *                    IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                    IINFO(I),D1,D2)
              IF (IINFO(I).NE.0) THEN
                DINFO1(I) = DWORK(L(LU)+IINFO(I)-1)
                DINFO2(I) = DWORK(L(LV)+IINFO(I)-1)
              ELSE
C               Oops, thats a point anywhere in the domain...
C               Search the element where it is and evaluate there.
                CALL PSRCH5(TRIAS,NLMIN,NLMAX,D1,D2,0,IEL)
                IF (IEL.EQ.0) THEN
                  DINFO1(I) = 0.0
                  DINFO2(I) = 0.0
                ELSE
                  KCORVG = L(TRIAS(OLCORVG,NLMAX))
                  KVERT  = L(TRIAS(OLVERT,NLMAX))
                  KMID   = L(TRIAS(OLMID,NLMAX))
                  CALL SCEVLQ (TRIAS(ONVT,NLMAX),DWORK(L(LU)),
     *                         E011,.FALSE.,D1,D2,IEL,
     *                         TRIAS(1,NLMAX),DWORK(KCORVG),
     *                         KWORK(KVERT),KWORK(KMID),
     *                         DINFO1(I), DTMP1, DTMP2) 
                  CALL SCEVLQ (TRIAS(ONVT,NLMAX),DWORK(L(LV)),
     *                         E011,.FALSE.,D1,D2,IEL,
     *                         TRIAS(1,NLMAX),DWORK(KCORVG),
     *                         KWORK(KVERT),KWORK(KMID),
     *                         DINFO2(I), DTMP1, DTMP2) 
                END IF
              END IF
              
            END DO
            
C           Print the results to screen

            CALL TV2TER (MSHOW,MFILE,J,IINFO,DINFO1,DINFO2,2,
     *        'Tracked velocity point-values P(VELO)')

C           Write the values to the file:

            IF (ITRACK(OHVEL1).GT.0) THEN
              CALL TV2FIL (ITRACK(OHVEL1),J,IINFO,DINFO1,
     *                     IPPDAT(OTRKCNT).EQ.0,TIMENS,-1,'PT')
            END IF

            IF (ITRACK(OHVEL2).GT.0) THEN
              CALL TV2FIL (ITRACK(OHVEL2),J,IINFO,DINFO2,
     *                     IPPDAT(OTRKCNT).EQ.0,TIMENS,-1,'PT')
            END IF
            
          END IF ! J>0
          
        END IF ! HVEL1>0 or HVEL2 > 0

C       ----------
C       Pressure

        IF (ITRACK(OHPRES).GT.0) THEN
        
C         Calculate the interpolated pressure, if it's not calculated

          IF (LP.EQ.0) THEN
            CALL XINTPV (DUP(1+2*NEQU),TRIAS(1,NLMAX),LP)
          END IF

C         Ask the user defined routine in how many points we should
C         track the velocity; write the number into K.

          CALL FPTSIN(8,0,TIMENS,DASMBL(ORE),
     *              IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *              J,D1,D2)
     
C         We only track a limited number of points

          J = MIN(J,MAXTRK)
          
          IF (J.GT.0) THEN
          
C           Get the point numbers and the values to IINFO/DINFO1/DINFO2:

            DO I=1,J
              CALL FPTSIN(9,I,TIMENS,DASMBL(ORE),
     *                    IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                    IINFO(I),D1,D2)
              IF (IINFO(I).NE.0) THEN
                DINFO1(I) = DWORK(L(LP)+IINFO(I)-1)
              ELSE
C               Oops, thats a point anywhere in the domain...
C               Search the element where it is and evaluate there.
                CALL PSRCH5(TRIAS,NLMIN,NLMAX,D1,D2,0,IEL)
                IF (IEL.EQ.0) THEN
                  DINFO1(I) = 0.0
                ELSE
                  KCORVG = L(TRIAS(OLCORVG,NLMAX))
                  KVERT  = L(TRIAS(OLVERT,NLMAX))
                  KMID   = L(TRIAS(OLMID,NLMAX))
                  CALL SCEVLQ (TRIAS(ONVT,NLMAX),DWORK(L(LP)),
     *                         E011,.FALSE.,D1,D2,IEL,
     *                         TRIAS(1,NLMAX),DWORK(KCORVG),
     *                         KWORK(KVERT),KWORK(KMID),
     *                         DINFO1(I), DTMP1, DTMP2) 
                END IF
              END IF
            END DO
            
C           Print the results to screen

            CALL TV2TER (MSHOW,MFILE,J,IINFO,DINFO1,DINFO2,1,
     *        'Tracked pressure point-values P(PRES)')

C           Write the values to the file:

            IF (ITRACK(OHPRES).GT.0) THEN
              CALL TV2FIL (ITRACK(OHPRES),J,IINFO,DINFO1,
     *                     IPPDAT(OTRKCNT).EQ.0,TIMENS,-1,'PT')
            END IF

          END IF ! J>0
          
        END IF ! HPRES > 0

C       --------------
C       Streamfunction

        IF (ITRACK(OHSTRM).GT.0) THEN
        
C         Calculate the interpolated streamfunction, if 
C         it's not calculated

          IF (LISO.EQ.0) THEN
            CALL XU2ISO (DUP(1),DUP(1+1*NEQU),TRIAS(1,NLMAX),LISO)
          END IF

C         Ask the user defined routine in how many points we should
C         track the velocity; write the number into K.

          CALL FPTSIN(10,0,TIMENS,DASMBL(ORE),
     *              IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *              J,D1,D2)
     
C         We only track a limited number of points

          J = MIN(J,MAXTRK)
          
          IF (J.GT.0) THEN
          
C           Get the point numbers and the values to IINFO/DINFO1:

            DO I=1,J
              CALL FPTSIN(11,I,TIMENS,DASMBL(ORE),
     *                    IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                    IINFO(I),D1,D2)
              IF (IINFO(I).NE.0) THEN
                DINFO1(I) = DWORK(L(LISO)+IINFO(I)-1)
              ELSE
C               Oops, thats a point anywhere in the domain...
C               Search the element where it is and evaluate there.
                CALL PSRCH5(TRIAS,NLMIN,NLMAX,D1,D2,0,IEL)
                IF (IEL.EQ.0) THEN
                  DINFO1(I) = 0.0
                ELSE
                  KCORVG = L(TRIAS(OLCORVG,NLMAX))
                  KVERT  = L(TRIAS(OLVERT,NLMAX))
                  KMID   = L(TRIAS(OLMID,NLMAX))
                  CALL SCEVLQ (TRIAS(ONVT,NLMAX),DWORK(L(LISO)),
     *                         E011,.FALSE.,D1,D2,IEL,
     *                         TRIAS(1,NLMAX),DWORK(KCORVG),
     *                         KWORK(KVERT),KWORK(KMID),
     *                         DINFO1(I), DTMP1, DTMP2) 
                END IF
              END IF
            END DO
            
C           Print the results to screen

            CALL TV2TER (MSHOW,MFILE,J,IINFO,DINFO1,DINFO2,1,
     *        'Tracked streamfunction point-values P(FLUX)')

C           Write the values to the file:

            IF (ITRACK(OHSTRM).GT.0) THEN
              CALL TV2FIL (ITRACK(OHSTRM),J,IINFO,DINFO1,
     *                     IPPDAT(OTRKCNT).EQ.0,TIMENS,-1,'PT')
            END IF

          END IF ! J>0
          
        END IF ! OHSTRM > 0

C       -----------------------
C       Mean Integreal Pressure

        IF (ITRACK(OHIPRS).GT.0) THEN
        
C         In how many boundary components should the pressure integral
C         be evaluated? We store that number in J.

          CALL FPTSIN(12,0,TIMENS,DASMBL(ORE),
     *                IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                J,D1,D2)
     
C         We only track a limited number of points

          J = MIN(J,MAXTRK)
          
          IF (J.GT.0) THEN
          
C           Get the boundary component segments where to evaluate the
C           forces. For now we restrict to a maximum of MAXTRK segments.
C           Evaluate there!
C           If J is -1, we have to evaluate only once - for all fictitious
C           boundary components at the same time!

            DO I=1,J
            
C             Request data about the boundary component.
C             The number of the fict. bdry component is stored to K.
C             If J=-1, we set K=0 to evaluate everywhere.

              CALL FPTSIN(13,I,TIMENS,DASMBL(ORE),
     *                  IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                  K,D1,D2)
              
C             Save the boundary component number to IINFO:              
              
              IINFO(I) = K
     
C             Call the evaluation routine, calculate integral on that
C             boundary component into DINFO1. DINFO2 receives the length
C             of the boundary we integrate on.
     
              CALL BDINTX (DUP(1+2*NEQU),TRIAS(1,NLMAX),
     *                 KWORK(KVERT),KWORK(KMID),KWORK(KVBD),KWORK(KEBD),
     *                 KWORK(KBCT),DWORK(KCORVG), DWORK(KVBDP),
     *                 E010,.FALSE.,
     *                 K,D1,D2,1,DINFO1(I),DINFO2(I))
     
C             Divide the integral pressure by the length of the boundary
C             component to get the integral mean pressure

              IF (DINFO2(I).NE.0D0) THEN
                DINFO1(I) = DINFO1(I)/DINFO2(I)
              END IF
              
            END DO
            
C           Print the results to screen

            CALL TV2TER (MSHOW,MFILE,J,IINFO,DINFO1,DINFO2,1,
     *        'Tracked integral mean pressure values I(PRES)')

C           Write to tracking file, if handles are given.
C           Write horizontal and vertical force.

            IF (ITRACK(OHIPRS).GT.0) THEN
              CALL TV2FIL (ITRACK(OHIPRS),J,IINFO,DINFO1,
     *                     IPPDAT(OTRKCNT).EQ.0,TIMENS,-1,'BDC')
            END IF

          END IF
        
        END IF ! HPRES > 0

C       Increase our counter that counts how many solutions we have
C       written out. That way, the headlines are e.g. only written
C       once!

        IPPDAT(OTRKCNT) = IPPDAT(OTRKCNT)+1
      
      END IF

C     =================================================================
C     Pre-/postprocessing regarding writing intermediate solutions to
C     disc.
C
C     Write solution vector to disc, depending on the number of the
C     current time step.

      IF (IPPDAT(OISAV).NE.0) THEN

C       On the beginning of the algorithm, initialize our file counter:

        IF (IFINL.EQ.-20) THEN 
          IPPDAT(OIFILEN) = 0
          RETURN
        END IF

C       The first part writes out intermediate solutions to disc
C       in a step-based manner. This checks the CRITE flag of the
C       status parameter block in IPARAM and writes out the solution
C       into the ./NS/DX.xxxx files when CRITE reaches INSVN:

        IF ((IFINL.GE.0).AND.(IFINL.LT.5)) THEN

C         If INSAV=-1, we write every calculated solution to disc.
C         If INSAV>=0, we write only every INSAV'th accepted solution
C         to disc.

          IF ((IPPDAT(OINSAV).EQ.-1).OR.
     *        ((IPPDAT(OINSAV).GE.0).AND.(IFINL.GE.1))) THEN
          
            II = MAX(IPPDAT(OINSAV),1)
            JJ = MAX(IPPDAT(OINSVN),1)
          
C           Only write out every INSAV'th solution.
C           The sign of ISOL defines whether the solution is written in
C           formatted or unformatted form.
          
            IF ( MOD(IPARAM(OCRITE),II).EQ.0 ) THEN

C             We want to write it to disc on level ISAV, i.e. at
C             level:

              I = MAX(NLMIN,MIN(NLMAX,ABS(IPPDAT(OISAV))))
              
C             But before we can write, we must calculate it.
C             We can use e.g. one of the allocated memory vectors in VECDAT
C             on level I to store the vector, as these are only used when
C             solving linear equations - what we clearly not do here!
C             So we can e.g. store it to DWORK(J):

              J = L(VECDAT(OLTMP,I))
              
C             We have to copy two velocity and one pressure vector to that.
C             Because of the 2-level ordering, the first couple of entries
C             on the finest level belong to the coarser levels! So we
C             can concat together:

              K = VECDAT(ONU,I)
              CALL LCP1(DUP(1),DWORK(J),K)
              CALL LCP1(DUP(1+K),DWORK(J+NEQU),K)
              CALL LCP1(DUP(1+2*K),DWORK(J+2*NEQU),VECDAT(ONP,I))
              
C             Write it.
C             The sign of ISAV decides on whether we write in formatted
C             or unformatted form.
          
              IF (IPPDAT(OLFLSAV).NE.0) THEN
                CALL STPUT (IPPDAT(OLFLSAV), CFN)
              ELSE
                WRITE (CFN,'(A)') 'DUPSAV'
              END IF
              CALL PPWRVC (0,VECDAT(ONEQV,I),DWORK(J),
     *                     1+MOD(IPPDAT(OIFILEN),JJ),
     *                     IPPDAT(OISAV),CFN)

            END IF
            
C           Increase our counter

            IPPDAT(OIFILEN) = IPPDAT(OIFILEN) + 1
        
          END IF ! (INSAV=-1) OR (INSAV>=0)
          
        END IF ! (IFINL>=0) and (IFINL < 5) 

      END IF ! ISAV<>0
      
C     -----------------------------------------------------------------
C     We go on with handling the writing out of solutions,
C     but this time another mothod.
C
C     It's getting a little bit more specific now: We want to write
C     out solution vectors in a time-stepping manner.
C     If the simulation time reaches the next point in time in a
C     distance of DUPSAV, we want to write out the solution.
C
C     The main "problem" with this thing is the predictor-corrector
C     approach. Sometimes, a solution is rejected and the time-stepping
C     switches back to an earlier point in time - and so we have to
C     do, too!
C
C     At first: On the beginning of the algorithm, initialize file
C     counter and simulational time:

      IF (IFINL.EQ.-20) THEN 
        IPPDAT(OIFILN)  = IPPDAT(OIFINIT)
        DPPDAT(OTUPSAV) = TIMEST
      END IF

C     On the beginning of each predictor step, make a backup of the
C     file counter and of the simulation time where we have to write
C     out the next step:

      IF (IFINL.EQ.-2) THEN
        DPPDAT(OTUPSBK) = DPPDAT(OTUPSAV)
        IPPDAT(OIFILB)  = IPPDAT(OIFILN)
      END IF
      
C     If we have a new solution (accepted or not) and the time
C     reaches the next checkpoint, save it to a file.

      IF ((IFINL.EQ.1).AND.(TIMENS.GE.DPPDAT(OTUPSAV))) THEN
      
C       -------------------------------------------------------------
C       Should we write out the plain solution vector as it is?

        IF (IPPDAT(OIUPSAV).NE.0) THEN
        
C         We want to write it to disc on level IUPSAV, i.e. at
C         level:

          I = MAX(NLMIN,MIN(NLMAX,ABS(IPPDAT(OIUPSAV))))
          
C         But before we can write, we must calculate it.
C         We can use e.g. one of the allocated memory vectors in VECDAT
C         on level I to store the vector, as these are only used when
C         solving linear equations - what we clearly not do here!
C         So we can e.g. store it to DWORK(J):

          J = L(VECDAT(OLTMP,I))
          
C         We have to copy two velocity and one pressure vector to that.
C         Because of the 2-level ordering, the first couple of entries
C         on the finest level belong to the coarser levels! So we
C         can concat together:

          K = VECDAT(ONU,I)
          CALL LCP1(DUP(1),DWORK(J),K)
          CALL LCP1(DUP(1+K),DWORK(J+NEQU),K)
          CALL LCP1(DUP(1+2*K),DWORK(J+2*NEQU),VECDAT(ONP,I))
          
C         Write it. The sign of IUPSAV decides on whether the output is
C         in formatted or unformatted form.
          
          IF (IPPDAT(OLFLDUP).NE.0) THEN
            CALL STPUT (IPPDAT(OLFLDUP), CFN)
          ELSE
            WRITE (CFN,'(A)') 'DUP'
          END IF
          CALL PPWRVC (0,VECDAT(ONEQV,I),DWORK(J),IPPDAT(OIFILN),
     *                 IPPDAT(OIUPSAV),CFN)
     
        END IF
     
C       -------------------------------------------------------------
C       Should we write out interpolated velocities?

        IF (IPPDAT(OIFUSAV).NE.0) THEN
          
C         We want to write it to disc on level IFUSAV, i.e. at
C         level:

          I = MAX(NLMIN,MIN(NLMAX,ABS(IPPDAT(OIFUSAV))))

C         Calculate the interpolated velocity, if it's not calculated

          IF ((LU.EQ.0).OR.(LV.EQ.0)) THEN
          
            CALL XINTUV (DUP(1),DUP(1+NEQU),TRIAS(1,NLMAX),LU,LV)
          
C           Implement boundary conditions, since the interpolation
C           does not correctly handle boundaries:

            KCORVG = L(TRIAS(OLCORVG,NLMAX))
            KXNPR  = L(TRIAS(OLXNPR,NLMAX))
            CALL BDRCOR (DWORK(L(LU)),DWORK(L(LV)),
     *                   TRIAS(1,NLMAX),DWORK(KCORVG),
     *                   KWORK(KXNPR),UE,TIMENS,DASMBL(ORE),
     *                   IASMBL,DASMBL,IGEOM,DGEOM)
          END IF

C         Write out the values in the vertices for the X-velocity.
C         The sign of IFUSAV decisdes on whether the output is in
C         formatted or unformatted form.

          CALL STPUT (IPPDAT(OLFLDU), CFN)
          CALL PPWRVC (0,TRIAS(ONVT,I),DWORK(L(LU)),IPPDAT(OIFILN),
     *                 IPPDAT(OIFUSAV),CFN)

C         and the Y-velocity
          
          IF (IPPDAT(OLFLSAV).NE.0) THEN
            CALL STPUT (IPPDAT(OLFLDV), CFN)
          ELSE
            WRITE (CFN,'(A)') 'DV'
          END IF
          CALL PPWRVC (0,TRIAS(ONVT,I),DWORK(L(LV)),IPPDAT(OIFILN),
     *                 IPPDAT(OIFUSAV),CFN)

        END IF

C       -------------------------------------------------------------
C       Should we write out interpolated pressure values?

        IF (IPPDAT(OIFPSAV).NE.0) THEN
          
C         We want to write it to disc on level IFPSAV, i.e. at
C         level:

          I = MAX(NLMIN,MIN(NLMAX,ABS(IPPDAT(OIFPSAV))))

C         Calculate the interpolated pressure, if it's not calculated

          IF (LP.EQ.0) THEN
            CALL XINTPV (DUP(1+1*NEQU),TRIAS(1,NLMAX),LP)
          END IF

C         Write out the values in the vertices for the pressure.
C         Make use of the 2-level ordering to write out only level I.
C         The sign of IFPSAV decides on whether we write in formatted
C         or unformatted form.

          IF (IPPDAT(OLFLDP).NE.0) THEN
            CALL STPUT (IPPDAT(OLFLDP), CFN)
          ELSE
            WRITE (CFN,'(A)') 'DP'
          END IF
          
          CALL PPWRVC (0,TRIAS(ONVT,I),DWORK(L(LP)),IPPDAT(OIFILN),
     *                 IPPDAT(OIFPSAV),CFN)

        END IF
     
C       -------------------------------------------------------------
C       Should we write out the streamfunction?

        IF (IPPDAT(OIFXSAV).GT.0) THEN
          
C         We want to write it to disc on level IFPSAV, i.e. at
C         level:

          I = MAX(NLMIN,MIN(NLMAX,ABS(IPPDAT(OIFXSAV))))

C         Calculate the interpolated streamfunction, if 
C         it's not calculated

          IF (LISO.EQ.0) THEN
            CALL XU2ISO (DUP(1),DUP(1+1*NEQU),TRIAS(1,NLMAX),LISO)
          END IF

C         Write out the values in the vertices for the streamfunction.
C         Make use of the 2-level ordering to write out only level I.
C         The sign of IFXSAV decides on whether we write in formatted
C         or unformatted form.

          IF (IPPDAT(OLFLISO).NE.0) THEN
            CALL STPUT (IPPDAT(OLFLISO), CFN)
          ELSE
            WRITE (CFN,'(A)') 'DISO'
          END IF
          
          CALL PPWRVC (0,TRIAS(ONVT,I),DWORK(L(LP)),IPPDAT(OIFILN),
     *                 IPPDAT(OIFXSAV),CFN)

        END IF
     
C       -------------------------------------------------------------
C       That's it, all solution vectors written out.
C
C       Increase the counter
     
        IPPDAT(OIFILN) = IPPDAT(OIFILN) + 1
        
C       Initialize the next time checkpoint

        DPPDAT(OTUPSAV) = DPPDAT(OTUPSAV) + DPPDAT(ODUPSAV)
        
      END IF ! (IFINL=2) and (TIMENS >= TUPSAV)
      
C     Now it can happen that a solution is rejected. In this case,
C     the algorithm returns to a previous state.
C     We do this, too, by resetting our file-counter and file-time
C     to the previous backup. That way, rejected solutions are 
C     overwritten with new ones if we reach the next checkpoint.

      IF (IFINL.EQ.5) THEN
        DPPDAT(OTUPSAV) = DPPDAT(OTUPSBK) 
        IPPDAT(OIFILN)  = IPPDAT(OIFILB)
      END IF
      
C     =================================================================
C     The writing of the solutions is completed here.
C     We go on now to check if we have to write out GMV files.
C     The way of doing this is very similar to the writing of
C     solution vectors...
C
C     At first check if we should write out GMV files at all:

      IF (IPPDAT(OIGMV).GT.0) THEN
      
C       At first: On the beginning of the algorithm, initialize file
C       counter and simulational time:

        IF (IFINL.EQ.-30) THEN 
          IPPDAT(OIGMVN) = IPPDAT(OIGINIT)
          DPPDAT(OTGMV)  = TIMEST
          
C         Initialise tracers, regularly distributed

          IF (IPPDAT(OITRCES).GT.0) THEN
          
            CALL INRTRC(IPPDAT(OITRCES),IPPDAT(OLTRCS),
     *                  DWORK(L(TRIAS(OLCORVG,NLMAX))),
     *                  KWORK(L(TRIAS(OLVERT,NLMAX))),
     *                  TRIAS(ONEL,NLMAX))
     
          ELSE 
     
C           If tracers are not active (ITRCES<=0), set LTRCS=0
C           (ok, should already be 0 because of the initialisation),
C           so the tracers are not calculated.

            IPPDAT(OLTRCS) = 0
            
          END IF
          
        END IF

C       On the beginning of each predictor step, make a backup of the
C       file counter and of the simulation time where we have to write
C       out the next step:

        IF (IFINL.EQ.-2) THEN
          DPPDAT(OTGMVBK) = DPPDAT(OTGMV)
          IPPDAT(OIGMVB)  = IPPDAT(OIGMVN)
        END IF
      
C       If we have a new solution (accepted or not) and the time
C       reaches the next checkpoint, save it to a file.

        IF ((IFINL.EQ.1).AND.(TIMENS.GE.DPPDAT(OTGMV))) THEN

C         We want to write the GMV file to disc on level IGMV, i.e. at
C         level:

          I = MAX(NLMIN,MIN(NLMAX,IPPDAT(OIGMV)))

C         -------------------------------------------------------------
C         Use handle 69 for writing the GMV file
          
          K = 69
          
C         Open it

          IF (IPPDAT(OLFLGMV).NE.0) THEN
            CALL STPUT (IPPDAT(OLFLGMV), CFN)
          ELSE
            WRITE (CFN,'(A)') 'u'
          END IF
          CALL GMVOF0 (K,IPPDAT(OIGMVN),CFN)
          
          IF (IER.EQ.0) THEN
          
C           Write header and triangulation on level I.
C           Obtain NCELLS and NVERTS.

            CALL GMVHEA (K)
            CALL GMVTRI (K,TRIAS(1,I),1,NCELLS,NVERTS)
        
C           Write materials;
C           we don't give material names for now.

            CALL GMVMAT (K,TRIAS(1,I),0,NCELLS,TRIAS(ONEL,I),
     *                   -16,'',DFGMMC,IGEOM,DGEOM)
            CALL GMVMAT (K,TRIAS(1,I),1,NVERTS,
     *                   TRIAS(ONVT,I)+TRIAS(ONMT,I),
     *                   -16,'',DFGMMV,IGEOM,DGEOM)
        
C           Calculate the velocity field, pressure and the stream 
C           function if they are not calculated already:

            IF ((LU.EQ.0).OR.(LV.EQ.0)) THEN
            
              CALL XINTUV (DUP(1),DUP(1+NEQU),TRIAS(1,NLMAX),LU,LV)
            
C             Implement boundary conditions, since the interpolation
C             does not correctly handle boundaries:

              KCORVG = L(TRIAS(OLCORVG,NLMAX))
              KXNPR  = L(TRIAS(OLXNPR,NLMAX))
              CALL BDRCOR (DWORK(L(LU)),DWORK(L(LV)),
     *                     TRIAS(1,NLMAX),DWORK(KCORVG),
     *                     KWORK(KXNPR),UE,1D0,DASMBL(ORE),
     *                     IASMBL,DASMBL,IGEOM,DGEOM)
            END IF

            IF (LP.EQ.0) THEN
              CALL XINTPV (DUP(1+2*NEQU),TRIAS(1,NLMAX),LP)
            END IF

            IF (LISO.EQ.0) THEN
              CALL XU2ISO (DUP(1),DUP(1+1*NEQU),TRIAS(1,NLMAX),LISO)
            END IF

C           Write all these at the desired level to the GMV file:

            CALL GMVVEL (K,TRIAS(1,I),1,NVERTS,TRIAS(ONVT,I),
     *                   DWORK(L(LU)),DWORK(L(LV)))
            CALL GMVSCA (K,TRIAS(1,I),1,NVERTS,TRIAS(ONVT,I),
     *                   DWORK(L(LP)),'pressure')
            CALL GMVSCA (K,TRIAS(1,I),0,NCELLS,TRIAS(ONEL,I),
     *                   DUP(1+2*NEQU),'pressure')
            CALL GMVSCA (K,TRIAS(1,I),1,NVERTS,TRIAS(ONVT,I),
     *                   DWORK(L(LISO)),'streamfunction')
        
C           Calculate the H1-error to LERR

            IF (LERR.EQ.0) THEN
              CALL ZNEW(TRIAS(ONVT,NLMAX),1,LERR,'DERR  ')
            END IF

            CALL ERPQH1(DWORK(L(LU)),DWORK(L(LV)),TRIAS(ONVT,NLMAX),
     *                  TRIAS(1,NLMAX),DWORK(L(LERR)))
            CALL GMVSCA (K,TRIAS(1,I),1,NVERTS,TRIAS(ONVT,I),
     *                   DWORK(L(LERR)),'H1errorVel')
            CALL ERPCH1(DWORK(L(LP)),TRIAS(ONVT,I),TRIAS(1,NLMAX),
     *                  DWORK(L(LERR)))
            CALL GMVSCA (K,TRIAS(1,I),1,NVERTS,TRIAS(ONVT,I),
     *                   DWORK(L(LERR)),'H1errorP')
        
C           Do we have tracers?

            IF (IPPDAT(OLTRCS).NE.0) THEN
            
C             Move the tracers according to the previously calculated
C             velocity field (at time TIMENS-TSTEP) to calculate
C             their position in the current time step TIMENS.
C             In the very first time step, the tracer velocity is 0.

              CALL TRCMOV(IPPDAT(OLTRCS),TRIAS,NLMIN,NLMAX,TSTEP)
              
C             Recycle tracers outside of the domain

              KCORVG = L(TRIAS(OLCORVG,NLMAX))
              KVERT  = L(TRIAS(OLVERT,NLMAX))
              CALL TRCREC(IPPDAT(OLTRCS),TRIAS(ONEL,NLMAX),
     *                    DWORK(KCORVG),KWORK(KVERT),
     *                    TRIAS(ONVBD,NLMAX),
     *                    KWORK(L(TRIAS(OLEBD,NLMAX))))
            
C             Write out the tracers.
            
              CALL GMVTRS (K,IPPDAT(OLTRCS))
              
C             After writing, calculate the next new tracer 
C             velocity field for the next time step using DU/DV.

              CALL TRCVUP (IPPDAT(OLTRCS),TRIAS,NLMIN,NLMAX,
     *                     TRIAS(ONVT,NLMAX),DWORK(L(LU)),DWORK(L(LV)),
     *                     E011,.FALSE.)

            END IF
          
C           Write geometry:

            CALL FBDGMV (K,2,IGEOM,DGEOM)
            
C           Write the simulation time to the GMV file

            CALL GMVTIM (K,TIMENS)
          
C           Write the footer, finish.

            CALL GMVFOT (K)
            CLOSE (K)
        
          END IF
          
C         -------------------------------------------------------------
C         That's it, GMV file is written out.
C         Increase the counter
     
          IPPDAT(OIGMVN) = IPPDAT(OIGMVN) + 1
          
C         Initialize the next time checkpoint

          DPPDAT(OTGMV) = DPPDAT(OTGMV) + DPPDAT(ODTGMV)
        
        END IF ! ! (IFINL=2) and (TIMENS >= TUPSAV)     

C       Now it can happen that a solution is rejected. In this case,
C       the algorithm returns to a previous state.
C       We do this, too, by resetting our file-counter and file-time
C       to the previous backup. That way, rejected solutions are 
C       overwritten with new ones if we reach the next checkpoint.

        IF (IFINL.EQ.5) THEN
          DPPDAT(OTGMV)  = DPPDAT(OTGMVBK) 
          IPPDAT(OIGMVN) = IPPDAT(OIGMVB)
        END IF
        
C       After the last time step, delete the tracers from the heap
C       if we have some.

        IF (IFINL.EQ.20) THEN 
          IF (IPPDAT(OLTRCS).NE.0) THEN
            CALL ZDISP(0,IPPDAT(OLTRCS),'TRCS  ')
          END IF
        END IF
      
      END IF ! IGMV > 0
      
C     =================================================================
C     Error analysis with analytic function.
C
C     Call the error analysis routine to print current errors
C     to the terminal. This is done only for accepted solutions.

      IF ((IFINL.EQ.0).AND.(IPPDAT(OIERANA).NE.0)) THEN
        CALL ERATRM (MFILE,IPPDAT(OIERANA),VECDAT(1,NLMAX),
     *               DUP,DRHS,DAUX,
     *               NLMIN,NLMAX,TRIAS,TIMENS,
     *               IASMBL,DASMBL,IGEOM,DGEOM,
     *               UE,PE,UEX,UEY)
      END IF
      
C     Release any memory we might have used

      IF (LERR.NE.0) CALL ZDISP(0,LERR,'DERR  ')
      IF (LISO.NE.0) CALL ZDISP(0,LISO,'DISO  ')
      IF (LU.NE.0) CALL ZDISP(0,LU,'DU    ')
      IF (LV.NE.0) CALL ZDISP(0,LV,'DV    ')
      IF (LP.NE.0) CALL ZDISP(0,LP,'DP    ')
      
      END
      
************************************************************************
* Standard time output for nonstationary solver
*
* This routine prints out the timing statistics on the nonstationary
* solver to terminal/file.
*
* In:
*   MSHOW  : Level of output
*   TMINI  : Time needed for preprocessing
*   TMTOT  : Total time of the algorithm including initialization/
*            postprocessing/everything.
*   IPARAM : array [1..SZISDI] of integer
*   DPARAM : array [1..SZISDD] of double
*            Integer and double prec. parameter blocks that define the
*            behaviour of the nonstationary solver. 
*
* Out:
*   If MSHOW>=0, the output is written to the file.
*   If MSHOW>=2, the output is written to the standard terminal
************************************************************************

      SUBROUTINE STATIS (MSHOW,TMINI,TMTOT,IPARAM,DPARAM)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'

      INCLUDE 'ssolvers.inc'
      
      INCLUDE 'stiming.inc'
      INCLUDE 'snsdef.inc'
      INCLUDE 'ststepping.inc'
      INCLUDE 'snonstead.inc'
      
      INTEGER MSHOW,IPARAM(*)
      DOUBLE PRECISION DPARAM(*)
      
      INTEGER IOFS
      DOUBLE PRECISION TTMG,T,TMINI,TMTOT
      CHARACTER CSTR*(255)
      
      IOFS = OTNSTIM-1
      
      IF (MSHOW.GT.0) THEN

        WRITE(CSTR,'(A)')
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A)') 'STATISTICS :'
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

        WRITE(CSTR,'(A,I12)') 'NWORK :      ',NWORK
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,I12)') 'IWORK :      ',IWORK
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,I12)') 'IWMAX :      ',IWMAX
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,1)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C       DSTPAR(TMTOT) in the solver structure shows the total
C       time of the solver + postprocessing. But as we want to
C       show the actual total time for everything, we print
C       out TMTOT from the parameter, not from DSTPAR!

        WRITE(CSTR,'(A,F20.10)') 'total time       : ', TMTOT
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        
        WRITE(CSTR,'(A,F20.10)') 'init. time       : ',TMINI
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        
        T = DPARAM(OTTGRID     ) + DPARAM(IOFS+OTTPOST)+
     *      DPARAM(IOFS+OTTADF ) + DPARAM(IOFS+OTTUPW )+
     *      DPARAM(IOFS+OTTBDR ) + DPARAM(IOFS+OTTLC  )+
     *      DPARAM(IOFS+OTTLSOL)
        
        WRITE(CSTR,'(A,F20.10)') 'appr. time       : ', T
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,F20.10)') '-> grid   time   : ', DPARAM(OTTGRID )
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

        T = DPARAM(IOFS+OTTADF )+DPARAM(IOFS+OTTUPW )+
     *      DPARAM(IOFS+OTTBDR )+DPARAM(IOFS+OTTLC  )
     
        WRITE(CSTR,'(A,F20.10)') '-> lin.   time   : ',    T
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,F20.10)') '   -> mavec time : ',
     *                           DPARAM(IOFS+OTTADF)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,F20.10)') '   -> konv. time : ',
     *                           DPARAM(IOFS+OTTUPW)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,F20.10)') '   -> bdry  time : ',
     *                           DPARAM(IOFS+OTTBDR)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,F20.10)') '   -> LC    time : ',
     *                           DPARAM(IOFS+OTTLC)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

        WRITE(CSTR,'(A,F20.10)') '-> l.sol. time   : ',
     *                            DPARAM(IOFS+OTTLSOL)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,F20.10)') '-> post.  time   : ', 
     *                            DPARAM(IOFS+OTTPOST )
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

        WRITE(CSTR,'(A,I12)') '#timesteps            : ',IPARAM(OTITTS)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,I12)') '#nonlinear iterations : ',IPARAM(OTITNL)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,I12)') '#linear iterations    : ',IPARAM(OTITLI)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,1)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        
        TTMG = DPARAM(IOFS+OTTMG)
        IF (TTMG.EQ.0D0) TTMG=1D0

        WRITE(CSTR,'(A)')
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A)') ' MULTIGRID COMPONENTS [in percent]:'
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,F20.10)') ' smoothing     :', 
     *                 1.D2*DPARAM(IOFS+OTTSMTH)/TTMG
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,F20.10)') ' solver        :', 
     *                 1.D2*DPARAM(IOFS+OTTCGC )/TTMG
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,F20.10)') ' defect calc.  :', 
     *                 1.D2*DPARAM(IOFS+OTTDEF )/TTMG
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,F20.10)') ' prolongation  :', 
     *                 1.D2*DPARAM(IOFS+OTTPROL)/TTMG
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,F20.10)') ' restriction   :', 
     *                 1.D2*DPARAM(IOFS+OTTREST)/TTMG
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,1)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        
      ENDIF

1     FORMAT(79('-'))

      END

************************************************************************
* Initialise postprocessing structure for nonstationary
* postprocessing routines.
*
* This initializes the IPARAM/DPARAM structure arrays, which are
* passed through the solver NONST2 to the postprocessing routines above.
*
* If desired, values of the COMMON-blocks are transferred to these 
* structure arrays. In this case, memory will be allocated for the
* filenames of the output files.
*
* In:
*   IPARAM : array [1..SZNSPI] of integer
*            TInstatPostprocessingIParams with Integer parameters.
*   DPARAM : array [1..SZNSPD] of double
*            TInstatPostprocessingDParams with Double precision 
*             parameters
*   IC2PAR : =0: initialize IPARAM/DPARAM only with standard values
*            =1: initialize IPARAM/DPARAM with standard values and
*                transfer values of COMMON-blocks into them
*
* Out:
*   IPARAM,
*   DPARAM : Initialized structures
*
************************************************************************

      SUBROUTINE ININSP (IPARAM,DPARAM,IC2PAR)
      
      IMPLICIT NONE

      INCLUDE 'dstrings.inc'
      
      INCLUDE 'snonsteadpostproc.inc'
      INCLUDE 'cpostproc.inc'
      
      INTEGER IPARAM(SZISPI),IC2PAR
      DOUBLE PRECISION DPARAM(SZISPD)

C     Initialise the structures completely with zero:

      CALL LCL3(IPARAM,SZISPI)
      CALL LCL1(DPARAM,SZISPD)
      
C     Transfer COMMON block variables?

      IF (IC2PAR.NE.0) THEN

        IPARAM (OISAV)   = ISAV
        IPARAM (OINSAV)  = INSAV
        IPARAM (OINSVN)  = INSAVN
        IPARAM (OIFINIT) = IFINIT
        IPARAM (OIUPSAV) = IUPSAV
        IPARAM (OIFUSAV) = IFUSAV
        IPARAM (OIFPSAV) = IFPSAV
        IPARAM (OIFXSAV) = IFXSAV
        IPARAM (OIGMV)   = IGMV
        IPARAM (OIGINIT) = IGINIT
        IPARAM (OITRCES) = ITRCES
        IPARAM (OIERANA) = IERANA

        DPARAM (ODUPSAV)  = DTFILM
        DPARAM (ODTGMV)   = DTGMV

C       Initialize standard calculation of body forces

        IPARAM (OIBDFBD) = 1
        IPARAM (OIBDFVI) = 4
        IPARAM (OIBDFLI) = 2
      
C       Initialize filenames for output files

        IPARAM (OLFLSAV)  = STNEWC (.TRUE.,CFLSAV)
        IPARAM (OLFLDUP)  = STNEWC (.TRUE.,CFLDUP)
        IPARAM (OLFLDU )  = STNEWC (.TRUE.,CFLDU )
        IPARAM (OLFLDV )  = STNEWC (.TRUE.,CFLDV )
        IPARAM (OLFLDP )  = STNEWC (.TRUE.,CFLDP )
        IPARAM (OLFLISO)  = STNEWC (.TRUE.,CFLISO)
        IPARAM (OLFLGMV)  = STNEWC (.TRUE.,CFLGMV)
      
      END IF
      
      END
      
************************************************************************
* Clean up postprocessing structure for nonstationary
* postprocessing routines.
*
* This routine cleans up the postprocessing structures, thus releasing
* memory that was used for storing filenames of output files.
*
* In:
*   IPARAM : array [1..SZNSPI] of integer
*            TInstatPostprocessingIParams with Integer parameters.
*
* Out:
*   IPARAM : All handles to dynamic information are released.
************************************************************************

      SUBROUTINE DONNSP (IPARAM)
      
      IMPLICIT NONE
      
      INCLUDE 'snonsteadpostproc.inc'
      
      INTEGER IPARAM(SZISPI)

C     Release the string handles for the basic filenames of the
C     output files.

      IF (IPARAM (OLFLSAV).NE.0) CALL STDIS (IPARAM (OLFLSAV))
      IF (IPARAM (OLFLDUP).NE.0) CALL STDIS (IPARAM (OLFLDUP))
      IF (IPARAM (OLFLDU ).NE.0) CALL STDIS (IPARAM (OLFLDU ))
      IF (IPARAM (OLFLDV ).NE.0) CALL STDIS (IPARAM (OLFLDV ))
      IF (IPARAM (OLFLDP ).NE.0) CALL STDIS (IPARAM (OLFLDP ))
      IF (IPARAM (OLFLISO).NE.0) CALL STDIS (IPARAM (OLFLISO))
      IF (IPARAM (OLFLGMV).NE.0) CALL STDIS (IPARAM (OLFLGMV))
      
      END
