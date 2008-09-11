************************************************************************
* Solver for the nonstationary Navier Stokes equations
*
* Extended calling convention.
*
* This routine performs a time-discretization on the nonstationary
* Navier Stokes equations.
*
* In:
*   NLMIN  : minimum level where calculation is allowed
*   NLMAX  : maximum level where calculation is allowed
*   TRIAS  : array [1..SZTRIA,1..NLEV] of integer
*            Triangulation structures for the underlying meshes
*            on all levels. TRIAS(.,ILEV) represents the triangulation
*            structure on level ILEV in its ummodified form, i.e.
*            without any boundary conditions implemented, probably
*            undeformed (if the user wants to start the grid
*            deformation in every time step from the original grids).
*   MATDAT : array [1..SZN2MI,1..NNLEV] of integer
*            TNS2DMatrixParams-structures for level NLMIN..NLMAX, 
*            initialized with data
*   VECDAT : array [1..SZN2VI,1..NNLEV] of integer
*            TNS2DVectorParams-structures for level NLMIN..NLMAX. 
*            This structure array must specify the structure of
*            the vectors on each level. Furthermore it must contain
*            handles to preallocated arrays. These arrays are used
*            as auxiliary arrays for solving a linear system on each 
*            level.
*            VECDAT is used as auxiliary structure when solving
*            the system on the finest level with the Multigrid
*            solver. The content of the vectors in VECDAT are treated
*            as undefined when entering this routine and will be
*            undefined when this routine is left.
*
*   IPARAM : array [1..SZISDI] of integer
*   DPARAM : array [1..SZISDI] of double
*            Integer and double prec. parameter blocks that define the
*            behaviour of the nonstationary solver. 
*   ISTPAR : array [1..SZNSDI] of integer
*   DSTPAR : array [1..SZNSDD] of double
*            Integer and double prec. parameter block for the stationary 
*            sub-solver NSDEF2.
*            The input-variables of the array must be initialized, e.g.
*            with INISTS.
*            The structures will be modified problem-specifically for the
*            solution of the subproblems, but restored onto return
*            of the routine.
*   IMGPAR : array [1..SZ020I+3*SZSLVI] of integer
*   DMGPAR : array [1..SZ020D+3*SZSLVD] of double
*            Integer and double parameter blocks for the multigrid 
*            sub-solver M020. These must be initialized as described
*            in NSDEF2.F. The variables in these blocks are not used
*            in MGSTP2, they are directly passed to NSDEF2.
*            The structures will be modified problem-specifically for the
*            solution of the subproblems, but restored onto return
*            of the routine.
*   IASMBL : array [1..SZASMI] of integer
*   DASMBL : array [1..SZASMD] of double
*            Integer and double prec. parameter block that controls the
*            discretization. This tells all assembly-routines how to 
*            set up the nonlinearity in the system matrices, which 
*            cubature formula to use, etc.
*            The discretization structures are modified during the
*            time stepping according to the situation of the simulation,
*            but restored to their original state upon returning of
*            this routine.
*   IGEOM  : array [1..*] of integer 
*   DGEOM  : array [1..*] of double 
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
*            routines; contains input and status variables for PREPOS.
*   ITRACK : array {1..SZTRKI] of integer
*            TTrackingIParams structure that holds the output handles
*            for writing tracked solution values to files.
*
*   NUVP   : Total nength of solution vector; usually 2*NEQU+NEQP
*   DUP    : array [1..NUVP] of double
*            Start vector for nonlinear iteration, DUP=(u,v,p)
*            = u^n in a time-dependent simulation.
*            Remark: DUP must not coincide with VECDAT[NLMAX].LSOL !
*
*   TIMENS : Initial simulation time; usually = 0D0
*
*   GENRHS : Callback routine that generates a RHS for a time step.
*            SUBROUTINE GENRHS (TRIA,IPARAM,DPARAM,
*                               IASMBL,DASMBL,IGEOM,DGEOM,VECDAT,DRHS)
*            Must not implement any boundary conditions into the RHS.
*
*   PREPOS : Callback routine that performs the pre- and postprocessing
*            after each time step. DFPSTD is the default routine.
*            For a documentation of the parameters, see DFPSTD!
*            SUBROUTINE PREPOS (NLMIN,NLMAX,
*                               TRIAS,MATDAT,VECDAT,
*                               IPARAM,DPARAM,ISTPAR,DSTPAR,
*                               IMGPAR,DMGPAR,IASMBL,DASMBL,
*                               IGEOM,DGEOM,IADTS,IPPDAT,DPPDAT,ITRACK,
*                               NUVP,DUP,DRHS,DAUX,
*                               TIMENS, TSTEP, ICORST, NCORST,
*                               ITNSR,STSTP1,STSTP2,IFINL,ISTATN,IREP)
*
*   PRCNTS : Proceed to next time step.
*            This routine is called onb one hand directly at the 
*            beginning of the simulation and on the other hand
*            directly before the time step increases from 
*            TIMENS to TIMENS+TSTEP. It allowes the
*            caller to make all arranges for the next time step,
*            e.g. modify the geometry, position of objects,...
*            As default routine, DFPNTS can be used!
*            SUBROUTINE PRCNTS (NLMIN,NLMAX,
*                               TRIAS,MATDAT,VECDAT,
*                               IPARAM,DPARAM,ISTPAR,DSTPAR,
*                               IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
*                               IADTS,DADTS,ITRACK,
*                               NUVP,DUP,DRHS,DAUX,
*                               TIMEST, TIMENS, TSTEP)
*
*   GENGRI : Callback routine that performs the grid adaption/generation
*            before each time step. Is only called if mesh generation
*            is activated by IPARAM.IMESH<>0. Can be set to I000 if
*            not used.
*            The routine is called after setting TIMENS to the next
*            simulation time. It can adapt the grid TRIAS to the current
*            geometry settings. If necessary, it can modify the RHS of
*            the previous time step or the current solution vector
*            to represent valid values on the new grid.
*            For a description of the parameters see DFPSTD, which
*            uses nearly an identical list of parameters. 
*            SUBROUTINE GENGRI (NLMIN,NLMAX,TRIAS,
*                               IPARAM,DPARAM,ISTPAR,DSTPAR,
*                               IMGPAR,DMGPAR,IASMBL,DASMBL,
*                               IGEOM,DGEOM
*                               VECDAT,DUP,DRHS,DAUX,
*                               TIMEST, TIMENS, TSTEP)
*
* Out:
*   TIMENS : Final simulation time
*   DUP    : receives the new approximation to the solution vector
*
* For time-dependent simulation (IPARAM.INONST=1):
*   DUP    : receives the final iterate of the last time step
* 
************************************************************************
* Customization - Implementational details
* ----------------------------------------
* Here we give a general overview about some implementational details
* that differ from the standard implementation of structures/algorithms.
* Some of the structures are customized to the special situation
* of CC2D with moving boundaries and geometry support.
*
* 1.) Triangulation structures
* ----------------------------
* The STRIA-structures for all levels are generally kept as described in
* STRIA.INC. The following optional arrays are activated and used
* during the nonlinear iteration:
*   TRIA.DMBDP  - boundary midpoints
*   TRIA.DCORMG - midpoint coordinates
*   TRIA.DAREA  - Area of the elements
*
* The user-defined part of the TRIA structure is maintained as follows:
*   TRIA.TRIUD(1) = NBDMT
*                 = Number of boundary vertices
*   TRIA.TRIUD(2) = LSCNPR
*                 = Handle to array [1..NBDMT] of integer = KSCNPR
*     Shortcut nodal property for all edges. This contains a list of
*     all boundary nodes. The first NVBD entries contain the numbers
*     of all boundary edges. Then all fictitious boundary edges
*     follow, attached to that list. NBDMT is the actual size of
*     this array, as it may vary over time if the fictitious boundary
*     changes. As the maximum number of edges in the geometry is NMT,
*     NBDMT<=NMT. The entries follow the following rule:
*      |Entry|   : Number of the edge on the boundary
*      Entry > 0 : The edge is a Dirichlet edge
*      Entry < 0 : The edge is a Neumann edge
*   TRIA.TRIUD(3) = LIARR
*                   Handle to precalculated integer information about
*                   the current geometry or 0 if not used
*   TRIA.TRIUD(4) = LDARR
*                   Handle to precalculated double information about
*                   the current geometry or 0 if not used
*
************************************************************************

      SUBROUTINE NONST2 (NLMIN,NLMAX,
     *                   TRIAS,MATDAT,VECDAT,
     *                   IPARAM,DPARAM,
     *                   ISTPAR,DSTPAR,
     *                   IMGPAR,DMGPAR,
     *                   IASMBL,DASMBL,IGEOM,DGEOM,
     *                   IADTS,DADTS,IPPDAT,DPPDAT,ITRACK,
     *                   NUVP,DUP,
     *                   TIMENS,
     *                   GENRHS,GENGRI,PREPOS,PRCNTS)

      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'ssolvers.inc'
      
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'
      
      INCLUDE 'm020.inc'

      INCLUDE 'stiming.inc'
      INCLUDE 'snsdef.inc'
      INCLUDE 'ststepping.inc'
      INCLUDE 'snonstead.inc'
      
      INCLUDE 'sassembly.inc'
      
      INCLUDE 'sadtstep.inc'
      
      include 'cmem2.inc'
      
C     externals
C
C     Parametrization of the domain
      DOUBLE PRECISION PARX,PARY,TMAX
      EXTERNAL PARX,PARY,TMAX
      
C     definition of finite elements
      EXTERNAL E030,E031,EM30,EM31,E010
      
C     dirichlet values in the domain

      EXTERNAL UE
      DOUBLE PRECISION UE
      
C parameters

      INTEGER TRIAS(SZTRIA,NNLEV),NUVP,ITRACK(*),IPPDAT(*)
      INTEGER IPARAM (SZISDI),IMGPAR(*),ISTPAR(*),IASMBL(*),IADTS(*)
      INTEGER MATDAT(SZN2MI,NNLEV),VECDAT(SZN2VI,NNLEV),NLMIN,NLMAX
      DOUBLE PRECISION DPARAM (SZISDD),DMGPAR(*),DSTPAR(*),DASMBL(*)
      DOUBLE PRECISION DADTS(*),DPPDAT(*),TIMENS,DTMP,DUP(*)
      INTEGER IGEOM(*)
      DOUBLE PRECISION DGEOM(*)
      
C     Generation of RHS, grid generation, Postprocessing

      EXTERNAL GENRHS, GENGRI, PREPOS, PRCNTS
      
C     local variables

      INTEGER I,J,LTM3,NMG,NNL,STPST1,STPST2,ICONV,NITNS,ITNS,LAUX
      INTEGER STSTP1,STSTP2,ISTATN,LRHS,INEUM,NEQU,NEQP
      DOUBLE PRECISION TIMEST
      
C     variables that control the behaviour and save the state of the
C     adaptive time stepping
      
      INTEGER LBKRHS,LTM0
      INTEGER TSORD,ICORST,NCORST,ITSSRP,IEXPLT,IADTIM,BCKITE
      DOUBLE PRECISION BCKTNS
      
C     DTIMIN is used as auxiliary array for GTMAUX.
C     DTIMBK saves a backup of the timing statistics after each time step,
C     so we can calculate the difference to the last time step.
      
      DOUBLE PRECISION DTIMIN(SZISDD),DTIMBK(SZISDD)
      
C     Current triangulation for a time step

      INTEGER CURTRI(SZTRIA,NNLEV)
      
C     Configuration of current time step in Theta scheme

      DOUBLE PRECISION TIMSTP(SZTSTD),TSTEP
      
      INTEGER NNITR,ITNSR
      INTEGER IRPSTP,IUDREP

      DOUBLE PRECISION DRELTN,CTPROJ,RELU2,RELP2,RELUM
      DOUBLE PRECISION RELPM,RELT

C     Clear timing information

      IPARAM(OTITTS) = 0
      IPARAM(OTITNL) = 0
      IPARAM(OTITLI) = 0

      CALL LCL1(DPARAM(OTNSTIM),SZTIMG)
      DPARAM(OTMTOT) = 0D0
      
      CALL LCL1(DTIMBK,SZISDD)
      
C     Stop the complete time of this algorithm.
C     The initial wall clock time is saved to DTIMIN(OTMTOT).
      
      CALL GTMAUX (DTIMIN,DPARAM,OTMTOT,0)
      
C     Get the vector size of velocity and pressure part on the finest
C     level; the information is frequently used.

      NEQU = VECDAT(ONU,NLMAX)
      NEQP = VECDAT(ONP,NLMAX)

C     Duplicate the original meshes into CURTRI. Depending on how we
C     deform the mesh, we don't have to duplicate everything...

C     For IMESH=0, we can share all information with TRIAS, since
C     nothing will be changed:

      IF (IPARAM(OIDMESH).EQ.0) J = NOT(0)
      
C     For time-dependend boundary conditions (which include time-
C     dependent fictitious boundary objects), make backups of KNPR/KXNPR
C     additionally to the vertex coordinates:

      IF (IPARAM(OIBDR).NE.0) J = NOT(2**7 + 2**15 + 2**0)

C     For IMESH=1-3 don't duplicate structural information (like KVERT), 
C     as this is not changed in the grid deformation process. Other
C     information like coordinates,... we have to back up so we can
C     restore it when necessary.

      IF ((IPARAM(OIDMESH).GE.1).AND.(IPARAM(OIDMESH).LE.3)) 
     *  J = NOT(2**0 + 2**1 + 2**7 + 2**12 + 2**13 + 2**15)
     
C     For IMESH=4, we expect (if ever implemented) that everything might
C     change - so backup everything

      IF (IPARAM(OIDMESH).GE.4) J = 0
      
C     I hope i evaluated the above bitfield correctly :)
      
C     Start the grid duplication for all levels
      
      CALL GTMAUX (DTIMIN,DPARAM,OTTGRID,0)
      
      DO I=NLMAX,NLMIN,-1

C       Copy DCORVG only on the maximum level, since the information
C       is shared between all levels!

        IF (I.EQ.NLMAX) THEN
          CALL TRIDUP(TRIAS(1,I),CURTRI(1,I),0,J)
        ELSE
          CALL TRIDUP(TRIAS(1,I),CURTRI(1,I),0,IOR(J,2**0))
          CURTRI(OLCORVG,I) = CURTRI(OLCORVG,NLMAX)
        END IF
        
      END DO
      
      CALL GTMAUX (DTIMIN,DPARAM,OTTGRID,1)
      
C     For IMESH >= 2, we initialize the discretization structure to
C     use the ALE method for reconstruction of solutions. Then we
C     have to allocate some memory that stores the grid velocity.

      IF (IPARAM(OIDMESH).GE.2) THEN
        IASMBL(OIALE) = 1
        CALL ZNEW (2*TRIAS(ONVT,NLMAX),1,IASMBL(OLGRVEL),'GRVEL ')
      ELSE
        IASMBL(OIALE) = 0
      END IF
      
C     TIMEST will hold a backup of the initial simulational time.
C     This represents the start time of the simulation.
C     TIMENS, the current simulation time, will be increased during 
C     the simulation, every time, a time step is accepted.

      TIMEST = TIMENS
      
C     We save this also to DPARAM in case anyone is interested in that

      DPARAM(OTIMEST) = TIMEST

C     Fetch the initial time step length into TSTEP; will be modified
C     due to adaptive time stepping probably

      TSTEP = DADTS(OTIMSTP)

C     Initialization of the time stepping scheme. 
C     Theoretically we could do this directly inside of the code.
C     But a future enhancement could be to control the time stepping
C     outside of this routine using callback routines. To prepare the
C     support of this, we prepare the configuration of the time
C     stepping here - this can later be replaced by callback routines
C     if necessary without changing the rest of the code...
C
C     Do we have a 1st order scheme or a 2nd order scheme?

      IF ((IPARAM(OFRSTP).EQ.1).OR.(DPARAM(OSSTHETA).NE.1D0)) THEN
        TSORD = 2
      ELSE
        TSORD = 1
      END IF
      
C     Are we using adaptive time stepping at all?

      IF (IADTS(OIADTIM).NE.0) THEN
      
        IADTIM = 1
      
C       The adaptive time stepping consists of one predictor step
C       and (for now) 3 corrector steps (although these are no real
C       corrector steps...)

        NCORST = 3
      
C       Do we have to repeat the time step if necessary?

        IF (ABS(IADTS(OIADTIM)).EQ.3) THEN
          ITSSRP = 1
        ELSE
          ITSSRP = 0
        END IF
        
C       Should we use extrapolation in time?
        
        IEXPLT = IADTS(OIEXTIM)
          
      ELSE
      
        IADTIM = 0
        NCORST = 1
        ITSSRP = 0
        IEXPLT = 0
        
C       If we use Fractional Step, we always have 3 "corrector" steps -
C       regardless of whether we repeat the time step or not.
C       That is because the time-stepping scheme needs 3 substeps for one
C       "large" step of length 3*TSTEP.

        IF (IPARAM(OFRSTP).EQ.1) NCORST = 3
        
      END IF ! (IADTS(OIADTIM).NE.0)
      
C     Additional memory for an auxiliary vector during the calculation

      CALL ZNEW (NUVP,1,LAUX,'DAUX  ')
      
C     Allocate LTMP0 - a temporary vector that holds the solution
C     in the beginning of the current time step.
C     That way we can later either restore a previous time step
C     (if adaptive time stepping with repetition is active or
C     the time step is repeated because the solver broke down, resp.),
C     or we can compare two solutions if the nonstationary simulation
C     got stationary.

      CALL ZNEW(NUVP,-1,LTM0  ,'DMT0  ')
      IF (IER.NE.0) GOTO 99998
      
C     Allocate a vector for the right hand side. In the instationary
C     case, the right hand side is always created by the GENRHS
C     subroutine for the corresponding time step.
C
C     As we are in the beginning of the simulation, we have to generate
C     our first RHS vector f^(0) before we can start the time stepping...
C     But this comes later.

      CALL ZNEW(NUVP,1,LRHS,'DRHS  ')
      IF (IER.NE.0) GOTO 99998
      
C     Now comes the preparation of temporary variables for the adaptive
C     time stepping. If we use adaptive time stepping, we have to
C     remember the solution vector and the RHS vector of the beginning
C     of the time step, so we can switch back after the predictor step.

      IF (IADTIM.NE.0) THEN
      
C       For adaptive time-stepping we need a vector DTM3 that
C       saves the result of the predictor-step.

        CALL ZNEW(NUVP,-1,LTM3  ,'DMT3  ')
        IF (IER.NE.0) GOTO 99998

C       Allocate a vector to save a backup of the RHS vector in
C       case that the right hand side is not homogenuous (IRHS <> 0)
C       or time dependent.

        LBKRHS=0
        IF (IPARAM(OIRHS).GE.1) THEN
          CALL ZNEW(NUVP,1,LBKRHS  ,'DRHS2 ')
          IF (IER.NE.0) GOTO 99998
        ENDIF
        
      ELSE
      
C       Otherwise assign LTM3 and LBKRHS a dummy handle - just
C       to make sure that all routines using these handles
C       don't crash because the arrays don't exist

        LTM3   = LTM0
        LBKRHS = LTM0
        
      END IF

C     To summarize the use of the vectors allocated above:
C
C     DU      - current solution vector u^n, updated to u^(n+1)
C               during NSDEF
C     DTM0    - holds a copy of u^n so that it can be restores;
C               must be restored for recalculation in case of too
C               big time-step and when switching from prediction
C               step to standard step-size.
C     DTM3    - Receives the calculated solution of the predictor
C               step in a time-dependent calculation

C=======================================================================
C     Preparation of nonsteady loop
C=======================================================================

C     NNITR = maximum number of repetitions for IADTIM<>0.

      NNITR = MAX(1,ABS(IADTS(OIREPIT)))

C***********************************************************************

C     Set the number of time steps we have to compute.

      NITNS = IPARAM(ONITMAX)

C     Begin with the time stepping.
C
C     Start with step 0 and use a DO-WHILE until we are at the end.
C     Increase the counter on every beginning of a time step.
C     We not count the time steps directly, as in adaptive time
C     stepping it might happen, that we have to go backward in time
C     to repeat a time step:

      ITNS = 0

C     ICORST counts from 1 to NCORST. On 1 and NCORST, the adaptive
C     time stepping is handled. We initialize ICORST with 0 and add
C     1 on every start of a time step.

      ICORST = 0

C     ITNSR counts the number of repetitions we had for one time step.
C     It counts from 0 to NNITR. If NNITR is reached, the time stepping
C     loop is left with an error.
      
      ITNSR = 0
      
C     ICONV is a flag that can be set to 1 if convergence (to a stationary
C     limit) is detected. The time stepping loop will finish then.
      
      ICONV = 0
      
C     ISTATN is a status flag that allowes the postprocessing routine to
C     take influence to the preparation to the next time step. It's
C     analyzed in every time step by MGSTP to check if some calculation
C     can be skipped.
C     Standard is NOT(0) which lets MGSTP carry out all preparation 
C     tasks.

      ISTATN = NOT(0);
      
C     After the memory for all the things is allocated and the
C     time stepping is initialized, we can now start to prepare
C     the matrix/vector/triangulation data for the first time step.

C     -----------------------------------------------------------------
C     Preparations for the first time step
C     -----------------------------------------------------------------

C     User defined preprocessing before setting up data for the first
C     time step. This allowes e.g. to inform the grid adaption
C     routine about the current simulational time:

      IUDREP = 0
      CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,0)
      CALL PREPOS (NLMIN,NLMAX,
     *           CURTRI,MATDAT,VECDAT,
     *           IPARAM,DPARAM,ISTPAR,DSTPAR,
     *           IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *           IADTS,DADTS,IPPDAT,DPPDAT,ITRACK,
     *           NUVP,DWORK(L(LTM3)),
     *           DWORK(L(LRHS)),DWORK(L(LAUX)),
     *           TIMEST, TIMENS, 0D0, ICORST, NCORST,
     *           ITNSR,0,0,-30,
     *           ISTATN,IUDREP)
      CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,1)
      
C     Call the user defined routine for processing the time
C     stepping to inform the caller that we just started 
C     with the 0'th time step.

      CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,0)
      CALL PRCNTS (NLMIN,NLMAX,
     *             CURTRI,MATDAT,VECDAT,
     *             IPARAM,DPARAM,ISTPAR,DSTPAR,
     *             IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *             NUVP,DWORK(L(LTM3)),
     *             DWORK(L(LRHS)),DWORK(L(LAUX)),
     *             TIMEST, TIMENS, 0D0)
      CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,1)

C     Do we have to adapt the mesh?

      IF (IAND(ISTATN,1).NE.0) THEN
      
C       Only adapt it if geometric features are to be included.

        IF ((IPARAM(OIDMESH).GE.1).AND.(IPARAM(OIDMESH).LE.3)) THEN
        
C         Start the grid deformation.
C         
C         CURTRI is a copy of TRIA in this situation, so we
C         use that for the first adaption.
C         Generate the new grid into CURTRI.

          CALL GTMAUX (DTIMIN,DPARAM,OTTGRID,0)

          CALL GENGRI (NLMIN,NLMAX,CURTRI,CURTRI,0,
     *                 IPARAM,DPARAM,ISTPAR,DSTPAR,
     *                 IMGPAR,DMGPAR,IASMBL,DASMBL,
     *                 IGEOM,DGEOM,
     *                 VECDAT(1,NLMAX),NUVP,DUP,DWORK(L(LRHS)),
     *                 DWORK(L(LAUX)),
     *                 DPARAM(OTIMEST),DASMBL(OTIMENS), TSTEP)
     
          CALL GTMAUX (DTIMIN,DPARAM,OTTGRID,1)
     
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTBDR,1)

C         A modified grid leads directly to modified Laplace/Pressure/
C         Mass matrices! We have to rebuild everything, based on CURTRI!

          CALL GTMAUX (DTIMIN,DPARAM,OTTLC,1)
          CALL GENNSM (NLMIN,NLMAX,CURTRI,IASMBL,DASMBL,7,MATDAT)
          CALL GTMAUX (DTIMIN,DPARAM,OTTLC,1)
     
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTBDR,0)

        END IF

      END IF ! AND(ISTATN,1) <> 0

C     Implement the new boundary conditions into the mesh
C     information for the first time step.
      
      DO I=NLMIN,NLMAX
      
C       Duplicate the prototype of the old grid, reconstruct
C       time-independent boundary conditions:

        CALL TRIDUP (TRIAS(1,I),CURTRI(1,I),1,
     *               NOT(2**7 + 2**15))
      
      END DO ! I
      
C     Calculate boundary conditions, store them into CURTRI.
C
C     This will modify CURTRI.KXNPR according to Dirichlet/
C     Neumann boundary information. As KXNPR we specify KXNPR of
C     CURTRI, overwriting the old CURTRI.KXNPR where necessary.
C
C     Calculate the number of boundary vertices NBDMT and the
C     shortcut nodal property array - both are stored in the
C     user defined block at CURTRI.TRIUD(1) / CURTRI.TRIUD(2).
C
C     Furthermore this routine updates INEUM in the assembly
C     structure whether there are Neumann components at all in 
C     our problem.
C
C     As simulation time we use the time of the next time step -
C     we just switched to that. GNBDGE uses the time identifier
C     in DASMBL to identify the current simulation time.

      CALL GNBDGE (NLMIN,NLMAX,CURTRI,
     *             IASMBL,DASMBL,IGEOM,DGEOM)
    
C     Calculate RHS-vector f^(0) - without implementing any boundary
C     conditions, so it fits to the calling convention of MGSTP.
        
      CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,0)
      CALL GENRHS (CURTRI(1,NLMAX),
     *             IPARAM,DPARAM,IASMBL,DASMBL,IGEOM,DGEOM,
     *             VECDAT(1,NLMAX),DWORK(L(LRHS)))
      CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,1)
      
C     Implement initial Dirichlet boundary conditions into initial
C     solution vector u^0:

      CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTBDR,0)
      CALL BDRSTX (DUP(1),DUP(1+NEQU),
     *             KWORK(L(CURTRI(OLXNPR,NLMAX))),
     *             PARX,PARY,UE,DASMBL(OTIMENS),
     *             DASMBL(ORE),
     *             CURTRI(1,NLMAX),IASMBL,DASMBL,IGEOM,DGEOM)   
      CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTBDR,1)

C=======================================================================
C     Nonsteady loop
C=======================================================================
C
C     User defined preprocessing before the loop:

      CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,0)
      IUDREP = 0
      CALL PREPOS (NLMIN,NLMAX,
     *           CURTRI,MATDAT,VECDAT,
     *           IPARAM,DPARAM,ISTPAR,DSTPAR,
     *           IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *           IADTS,DADTS,IPPDAT,DPPDAT,ITRACK,
     *           NUVP,DWORK(L(LTM3)),
     *           DWORK(L(LRHS)),DWORK(L(LAUX)),
     *           TIMEST, TIMENS, 0D0, ICORST, NCORST,
     *           ITNSR,0,0,-20,
     *           ISTATN,IUDREP)
      CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,1)
      
C     Here it goes...
      
      DO WHILE ((ITNS.LT.NITNS).AND.
     *          (TIMENS.LT.DPARAM(OTIMEMX)).AND.
     *          ((IADTIM.LT.3).OR.(ITNSR.LE.NNITR)).AND.
     *          (ICONV.EQ.0))
      
C       IPARAM.TITTS counts the total number of time iterations -
C       regardless of whether a timestep is repeated or not.

        IPARAM(OTITTS) = IPARAM(OTITTS)+1
        
C       Increase the "correction step counter"

        ICORST = ICORST+1
        IF (ICORST.GT.NCORST) ICORST=1
      
C       Handle the adaptive time stepping. ICORST=1..NCORST acts as a 
C       state-machine that indicates the state of the time stepping.
C       In substep 1 we make a snapshot of the current solution/RHS, so we
C       can restore the solution after the predictor step and if we have
C       to repeat the time step.
      
        IF (ICORST.EQ.1) THEN  
         
C         Reset the status of the nonlinear solver in the prediction step
C         - whether we use it or not.

          STSTP1 = 0
      
C         Snapshot the solution
         
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,0)
         
          CALL LCP1(DUP,DWORK(L(LTM0)),NUVP)

          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,1)
            
C         Remember the time when we had that backup and the
C         iteration counter there. BCKTNS and BCKITE always correspond
C         to the solution vector in LTM0!

          BCKTNS = TIMENS
          BCKITE = ITNS
         
          IF (IADTIM.NE.0) THEN
          
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,0)
         
C           If the RHS is time dependent, snapshot the RHS

            IF (LBKRHS.NE.0) THEN
              CALL LCP1(DWORK(L(LRHS)),DWORK(L(LBKRHS)),NUVP)
            END IF
            
            CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,1)
            
          ENDIF ! (IADTIM.NE.0)
          
        END IF ! (ICORST.EQ.1)

C       Increase the current iteration counter to represent the
C       new time step. This is done after making the backup above,
C       so when repeating a time step and the backup is restored, 
C       ITNS will again represent the correct number of the step
C       after increasing it here.

        ITNS = ITNS + 1

C       We set IRPSTP=0. This indicates at the end of the routine, 
C       if a time step has to be repeated:
C        =0: time step accepted
C        =1: time step has to be repeated because of time error
C            out of bounds
C        =2: time step has to be repeated because of a critical
C            solver broke down in the corrector step
C        =3: time step has to be repeated because of a critical
C            solver broke down in the predictor step

        IRPSTP = 0
        
C       A little bit for nice output :)

        IF (MT.GE.2) WRITE(MTERM,*)
        IF (MT.GE.1) WRITE(MFILE,*)

C       User defined preprocessing at the beginning of the time step

        CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,0)
        IUDREP = 0
        CALL PREPOS (NLMIN,NLMAX,
     *             CURTRI,MATDAT,VECDAT,
     *             IPARAM,DPARAM,ISTPAR,DSTPAR,
     *             IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *             IADTS,DADTS,IPPDAT,DPPDAT,ITRACK,
     *             NUVP,DWORK(L(LTM3)),
     *             DWORK(L(LRHS)),DWORK(L(LAUX)),
     *             TIMEST, TIMENS, TSTEP, ICORST, NCORST,
     *             ITNSR,0,0,-10,
     *             ISTATN,IUDREP)
        CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,1)
      
C       In step ICORST=1, we handle the adaptive time stepping.
C       If this is active, we calculate a predicted solution...

        IF ((IADTIM.NE.0).AND.(ICORST.EQ.1)) THEN

C         ==============================================================
C         Predictor step
C         ==============================================================

C         At this point we know, that any kind of adaptive time-stepping
C         should be performed. At first we will calculate one large step
C         with time-step size 3*TSTEP. This is the predictor-step
C         in the adaptive time-stepping algorithm.
C
C         This calculation branch is only called if a really 
C         time-dependent simulation should be performed. In case of a
C         stationary simulation, IADTIM=0 and so this code is not
C         executed. In that case the solution will be calculated later.
C
C         Initialize the time stepping variables in the DPARAM block.
C         The predictor step always performs a simple, "quick-and-dirty"
C         calculation with Crank Nicolson or explicit Euler and
C         a time step that is 3x as large as the original. 
C         Remember,as this is only a predictor step, it's not so
C         important how accurate it is!
C
C         If the standard time stepping is 1st order, use Euler.
C         If the standard time stepping is 2nd order, use Crank
C         Nicolson.

          IF (TSORD.EQ.1) THEN
          
C           Use standard (1st order) time stepping for the predictor
C           step.
          
            CALL INITSC (TIMSTP,TIMENS,0,DPARAM(OSSTHETA),3D0*TSTEP,1)

          ELSE
          
C           Use Crank Nicolson as 2nd order predictor if 2nd order
C           accuracy is used.

            CALL INITSC (TIMSTP,TIMENS,0,0.5D0,3D0*TSTEP,1)
          END IF

C         Update IPARAM.CRITE to reflect the current step

          IPARAM(OCRITE) = ITNS

C         and DASMBL.TIMENS to reflect the current simulation time

          DASMBL(OTIMENS) = TIMENS

C         ... which will be overwritten in the MGSTP2 routine below,
C         as this goes on one time step.
C
C         Copy the current solution vector to DTM3. When the
C         calculation is successful, DTM3 therefore receives
C         the predicted solution
          
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,0)
          CALL LCP1(DUP,DWORK(L(LTM3)),NUVP)
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,1)

C         User defined preprocessing:

          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,0)
          IUDREP = 0
          CALL PREPOS (NLMIN,NLMAX,
     *               CURTRI,MATDAT,VECDAT,
     *               IPARAM,DPARAM,ISTPAR,DSTPAR,
     *               IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *               IADTS,DADTS,IPPDAT,DPPDAT,ITRACK,
     *               NUVP,DWORK(L(LTM3)),
     *               DWORK(L(LRHS)),DWORK(L(LAUX)),
     *               TIMEST, DASMBL(OTIMENS),TIMSTP(OTSTEP), 
     *               ICORST, NCORST,ITNSR,0,0,-2,
     *               ISTATN,IUDREP)
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,1)

C         Tell the user what we are going to do...

          IF (MT.GE.2) WRITE(MTERM,2)
          IF (MT.GE.2) 
     *      WRITE(MTERM,2001) ITNS,DASMBL(OTIMENS),TIMSTP(OTSTEP)
          IF (MT.GE.2) WRITE(MTERM,2)
          IF (MT.GE.1) WRITE(MFILE,2)
          IF (MT.GE.1) 
     *      WRITE(MFILE,2001) ITNS,DASMBL(OTIMENS),TIMSTP(OTSTEP)
          IF (MT.GE.1) WRITE(MFILE,2)

C         Start the calculation routine to calculate one 
C         macro-time-step with stepsize 3xTSTEP.
C         This will overwrite DTM3 and the RHS vector
C         (if we have instationary RHS).

          CALL MGSTP2(TRIAS,MATDAT,VECDAT,NLMIN,NLMAX,
     *                IPARAM,DPARAM,
     *                ISTPAR,DSTPAR,
     *                IMGPAR,DMGPAR,
     *                IASMBL,DASMBL,IGEOM,DGEOM,
     *                NUVP,DWORK(L(LTM3)),
     *                DWORK(L(LRHS)),DWORK(L(LAUX)),
     *                TIMSTP,CURTRI,
     *                NNL,NMG,ISTATN,STPST1,
     *                GENRHS,GENGRI,PRCNTS) 

C         Direct postprocessing of the time step, output statistical
C         data to terminal

          CALL MGSPOS (VECDAT(1,NLMAX),DWORK(L(LTM3)),
     *                 ITNS,1,DASMBL(OTIMENS))

C         Sum up the number of linear and nonlinear iterations

          IPARAM(OTITLI) = IPARAM(OTITLI)+NMG
          IPARAM(OTITNL) = IPARAM(OTITNL)+NNL

C         Result of the calculation:
C           DU     -> solution vector, updated from u^n to u^(n+1)
C           LRHS   -> updated RHS-vector
C           STPST1 -> indicates success of the solver in the
C                     prediction step
C         STPST1 we need later for the adaptive time stepping...
C
C         Call the postprocessing routine to let it analyze the
C         predictor step if desired:        
        
          IF (MT.GE.3) WRITE (MTERM,*) 
          IF (MT.GE.0) WRITE (MFILE,*) 
          IF (MT.GE.3) WRITE (MTERM,*) 
     *                'Starting postprocessing of predictor step...'
          IF (MT.GE.0) WRITE (MFILE,*) 
     *                'Starting postprocessing of predictor step...'
          IF (MT.GE.3) WRITE (MTERM,*) 
          IF (MT.GE.0) WRITE (MFILE,*) 

          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,0)
          IUDREP = 0
          CALL PREPOS (NLMIN,NLMAX,
     *               CURTRI,MATDAT,VECDAT,
     *               IPARAM,DPARAM,ISTPAR,DSTPAR,
     *               IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *               IADTS,DADTS,IPPDAT,DPPDAT,ITRACK,
     *               NUVP,DWORK(L(LTM3)),
     *               DWORK(L(LRHS)),DWORK(L(LAUX)),
     *               TIMEST, DASMBL(OTIMENS),TIMSTP(OTSTEP), 
     *               ICORST, NCORST, ITNSR,STSTP1,0,2,
     *               ISTATN,IUDREP)
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,1)
          
C         If the linear solver broke down for any reason (Bit 1 of 
C         STATUS <> 0). we directly set the repetition flag to 2 to
C         indicate that the step should be repeated if possible.
C
C         If only the nonlinear solver broke down (or if we don't have
C         repetitions left), let's hope that nonlinear solver of the
C         correction step still works -- so don't care...

          IF ((IAND(STPST1,2).NE.0).AND.(ITNSR.LE.NNITR)) THEN 
            
            IRPSTP = 3

C           Calculate a new stepsize into DTMP, based only 
C           on the result of the solver in the prediction step
C           (-> only STSTP1 has a value).
C           RELT, RELU2, RELUM, RELP2, RELPM serve as dummy variables
C           here.

            CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,0)
            CALL DTSCTR (IADTS, DADTS, 0,
     *                   NEQU, NEQP, NUVP, DWORK(L(LTM3)), DUP, 
     *                   DWORK(L(LAUX)),
     *                   RELT, RELU2, RELUM, RELP2, RELPM,
     *                   TIMEST,TIMENS,TSTEP,STPST1,0,ITNSR,TSORD,
     *                   IPARAM,DPARAM,
     *                   DTMP)
            CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,1)

C           Tell the user the new step size

            IF (MT.GE.2) WRITE(MTERM,3)
            IF (MT.GE.1) WRITE(MFILE,3)
            IF (MT.GE.2) WRITE(MTERM,3001) IRPSTP,ITNS,DTMP,TSTEP
            IF (MT.GE.1) WRITE(MFILE,3001) IRPSTP,ITNS,DTMP,TSTEP
            IF (MT.GE.2) WRITE(MTERM,3)
            IF (MT.GE.2) WRITE(MTERM,*)
            IF (MT.GE.1) WRITE(MFILE,3)
            IF (MT.GE.1) WRITE(MFILE,*)

C           and accept the step size.

            TSTEP = DTMP
            
C           While IRPSTP<>0, no further calculation will be done,
C           the step will be repeated directly.
C
C           We don't have to recover any information about RHS/solution.
C           The repetition handling will restore the vectors for us
C           later.

          ELSE IF (IUDREP.NE.0) THEN
          
C           The postprocessing routine forces us to repeat the step
          
            IRPSTP = 3
          
C           Tell the user that we are going to repeat

            IF (MT.GE.2) WRITE(MTERM,3)
            IF (MT.GE.1) WRITE(MFILE,3)
            IF (MT.GE.2) WRITE(MTERM,3001) IRPSTP,ITNS,DTMP,TSTEP
            IF (MT.GE.1) WRITE(MFILE,3001) IRPSTP,ITNS,DTMP,TSTEP
            IF (MT.GE.2) WRITE(MTERM,3)
            IF (MT.GE.2) WRITE(MTERM,*)
            IF (MT.GE.1) WRITE(MFILE,3)
            IF (MT.GE.1) WRITE(MFILE,*)
          
          ELSE

C           That worked, we now have a predicted solution in DTM3 and
C           the time step is still acceptable.
C
C           Now restore the old state so that the corrector steps will
C           work fine.
C          
C           If the RHS is really overwritten, restore it.

            IF (LBKRHS.NE.0) THEN
              CALL LCP1(DWORK(L(LBKRHS)),DWORK(L(LRHS)),NUVP)
            END IF
            
C           The solution vector must not be restored - we modified DTM3
C           and have not touched the original solutiuon :)

          END IF 

C         Some nice output to the user.
C
C         Stop the time. 
      
          CALL GTMAUX (DTIMIN,DPARAM,OTMTOT,1)
          CALL GTMAUX (DTIMIN,DPARAM,OTMTOT,0)
                
          IF (MT.GE.3) WRITE (MTERM,*) 
     *                'Timing statistics - Total/Last step'
          IF (MT.GE.0) WRITE (MFILE,*) 
     *                'Timing statistics - Total/Last step'
          IF (MT.GE.3) WRITE (MTERM,*) 
     *                '-----------------------------------'
          IF (MT.GE.0) WRITE (MFILE,*) 
     *                '-----------------------------------'
                                                       
          IF (MT.GE.3) WRITE (MTERM,*) 'total time : ',
     *      DPARAM(OTMTOT),'/',
     *      DPARAM(OTMTOT)-DTIMBK(OTMTOT)
          IF (MT.GE.0) WRITE(MFILE,*) 'total time : ', 
     *      DPARAM(OTMTOT),'/',
     *      DPARAM(OTMTOT)-DTIMBK(OTMTOT)
          IF (MT.GE.3) WRITE(MTERM,*) 'mavec time : ', 
     *      DPARAM(OTNSTIM-1+OTTADF),'/',
     *      DPARAM(OTNSTIM-1+OTTADF)-DTIMBK(OTNSTIM-1+OTTADF)
          IF (MT.GE.0) WRITE(MFILE,*) 'mavec time : ', 
     *      DPARAM(OTNSTIM-1+OTTADF),'/',
     *      DPARAM(OTNSTIM-1+OTTADF)-DTIMBK(OTNSTIM-1+OTTADF)
          IF (MT.GE.3) WRITE(MTERM,*) 'konv. time : ', 
     *      DPARAM(OTNSTIM-1+OTTUPW),'/',
     *      DPARAM(OTNSTIM-1+OTTUPW)-DTIMBK(OTNSTIM-1+OTTUPW)
          IF (MT.GE.0) WRITE(MFILE,*) 'konv. time : ', 
     *      DPARAM(OTNSTIM-1+OTTUPW),'/',
     *      DPARAM(OTNSTIM-1+OTTUPW)-DTIMBK(OTNSTIM-1+OTTUPW)
          IF (MT.GE.3) WRITE(MTERM,*) 'bdry  time : ', 
     *      DPARAM(OTNSTIM-1+OTTBDR),'/',
     *      DPARAM(OTNSTIM-1+OTTBDR)-DTIMBK(OTNSTIM-1+OTTBDR)
          IF (MT.GE.0) WRITE(MFILE,*) 'bdry  time : ', 
     *      DPARAM(OTNSTIM-1+OTTBDR),'/',
     *      DPARAM(OTNSTIM-1+OTTBDR)-DTIMBK(OTNSTIM-1+OTTBDR)
          IF (MT.GE.3) WRITE(MTERM,*) 'LC    time : ', 
     *      DPARAM(OTNSTIM-1+OTTLC),'/',
     *      DPARAM(OTNSTIM-1+OTTLC)-DTIMBK(OTNSTIM-1+OTTLC)
          IF (MT.GE.0) WRITE(MFILE,*) 'LC    time : ', 
     *      DPARAM(OTNSTIM-1+OTTLC),'/',
     *      DPARAM(OTNSTIM-1+OTTLC)-DTIMBK(OTNSTIM-1+OTTLC)
          IF (MT.GE.3) WRITE(MTERM,*) 'LinSl time : ', 
     *      DPARAM(OTNSTIM-1+OTTLSOL),'/',
     *      DPARAM(OTNSTIM-1+OTTLSOL)-DTIMBK(OTNSTIM-1+OTTLSOL)
          IF (MT.GE.0) WRITE(MFILE,*) 'LinSl time : ', 
     *      DPARAM(OTNSTIM-1+OTTLSOL),'/',
     *      DPARAM(OTNSTIM-1+OTTLSOL)-DTIMBK(OTNSTIM-1+OTTLSOL)

          IF (MT.GE.3) WRITE (MTERM,*) 
          IF (MT.GE.0) WRITE (MFILE,*) 
          
C         Back up the timing statistics for the next step.
C         Actually we make a backup of the whole parametre block, but
C         it contains the timing block and that's enough :)

          CALL LCP1(DPARAM,DTIMBK,SZISDD)

        END IF ! (IADTIM.NE.0)

C       ================================================================
C       Corrector steps
C       ================================================================

C       Don't do anything if it's already clear, that the time step has
C       to be repeated. Would be useless...

        IF (IRPSTP.EQ.0) THEN

C         Now we return to the calculation of the "real" solution with 
C         the standard time-step size TSTEP. The actual time step size,
C         which my vary from TSTEP according to the time stepping
C         scheme, will be initialized by the time stepping routines,
C         so we don't care here about that.
C         We simply calculate one time step.
C
C         Update IPARAM.CRITE to reflect the current step

          IPARAM(OCRITE) = ITNS

C         and DASMBL.TIMENS to reflect the current simulation time

          DASMBL(OTIMENS) = TIMENS

C         Use the time stepping routines to initialize the time
C         step weights. ICORST counts our current "microstep"
C         where we are in the scheme.

          CALL INITSC (TIMSTP,TIMENS,IPARAM(OFRSTP),DPARAM(OSSTHETA),
     *                 TSTEP,ICORST)
        
C         User defined preprocessing:

          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,0)
          IUDREP = 0
          CALL PREPOS (NLMIN,NLMAX,
     *               CURTRI,MATDAT,VECDAT,
     *               IPARAM,DPARAM,ISTPAR,DSTPAR,
     *               IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *               IADTS,DADTS,IPPDAT,DPPDAT,ITRACK,
     *               NUVP,DWORK(L(LTM3)),
     *               DWORK(L(LRHS)),DWORK(L(LAUX)),
     *               TIMEST, DASMBL(OTIMENS),TIMSTP(OTSTEP), 
     *               ICORST, NCORST,ITNSR,0,0,-1,
     *               ISTATN,IUDREP)
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,1)

C         Tell the user what we are going to do...

          IF (MT.GE.2) WRITE(MTERM,2)
          IF (MT.GE.2) 
     *      WRITE(MTERM,2002) ITNS,DASMBL(OTIMENS),TIMSTP(OTSTEP)
          IF (MT.GE.2) WRITE(MTERM,2)
          IF (MT.GE.1) WRITE(MFILE,2)
          IF (MT.GE.1) 
     *      WRITE(MFILE,2002) ITNS,DASMBL(OTIMENS),TIMSTP(OTSTEP)
          IF (MT.GE.1) WRITE(MFILE,2)

C         Start the calculation routine to calculate one 
C         time step:
          
          CALL MGSTP2(TRIAS,MATDAT,VECDAT,NLMIN,NLMAX,
     *                IPARAM,DPARAM,
     *                ISTPAR,DSTPAR,
     *                IMGPAR,DMGPAR,
     *                IASMBL,DASMBL,IGEOM,DGEOM,
     *                NUVP,DUP,
     *                DWORK(L(LRHS)),DWORK(L(LAUX)),
     *                TIMSTP,CURTRI,
     *                NNL,NMG,ISTATN,STPST2,
     *                GENRHS,GENGRI,PRCNTS)
     
C         Result of the calculation:
C           DUP  -> solution vector, updated from u^n to u^(n+1)
C           RHS  -> updated RHS-vector
C           DASMBL.TIMENS -> new simulation time
C           STPST2 -> indicates whether the solver in the correction
C                     step was successful; we need this for adaptive
C                     time stepping later.
C
C         Direct postprocessing of the time step, output statistical
C         data to terminal

          CALL MGSPOS (VECDAT(1,NLMAX),DUP,
     *                 ITNS,ICORST,DASMBL(OTIMENS))

C         Sum up the number of linear and nonlinear iterations

          IPARAM(OTITLI) = IPARAM(OTITLI)+NMG
          IPARAM(OTITNL) = IPARAM(OTITNL)+NNL
          
C         Fetch the new simulation time from the status variable
C         in the parameter block to TIMENS so that we are also
C         "up to date" here :)

          TIMENS = DASMBL(OTIMENS)

C         Corrector step finished
C
C         If the solver (linear or nonlinear) broke down, repeat the 
C         step if possible.

          IF (STPST2.NE.0) THEN 
          
            IRPSTP = 2
            
C           Calculate the new step sizeinto DTMP. If we are
C           using adaptive time stepping, also use the result of
C           the predictor step to choose a new time step. Otherwise
C           the result is only based on the success of the current
C           solver, and STSTP1=0 as it is initialized that way.
C           RELT, RELU2, RELUM, RELP2, RELPM serve as dummy variables
C           here.

            CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,0)
            CALL DTSCTR(IADTS, DADTS, 0,
     *                 NEQU, NEQP, NUVP, DWORK(L(LTM3)), DUP, 
     *                 DWORK(L(LAUX)),
     *                 RELT, RELU2, RELUM, RELP2, RELPM,
     *                 TIMEST,TIMENS,TSTEP,STPST1,STPST2,ITNSR,TSORD,
     *                 IPARAM,DPARAM,
     *                 DTMP)
            CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,1)
            
C           Tell the user the new step size

            IF (MT.GE.2) WRITE(MTERM,3)
            IF (MT.GE.1) WRITE(MFILE,3)
            IF (MT.GE.2) WRITE(MTERM,3001) IRPSTP,ITNS,DTMP,TSTEP
            IF (MT.GE.1) WRITE(MFILE,3001) IRPSTP,ITNS,DTMP,TSTEP
            IF (MT.GE.2) WRITE(MTERM,3)
            IF (MT.GE.2) WRITE(MTERM,*)
            IF (MT.GE.1) WRITE(MFILE,3)
            IF (MT.GE.1) WRITE(MFILE,*)

C           and accept the step size
            
            TSTEP = DTMP
            
C           Otherwise we can continue with the time step control:
            
          ELSE IF ((IADTIM.NE.0).AND.(ICORST.EQ.NCORST)) THEN
          
C           ==============================================================
C           Time step control
C           ==============================================================
C
C           The time step control takes place if it's active 
C           (of course :) ) and if ICORST reaches NCORST. In this case
C           we have everything together - the predicted solution and
C           the corrected one after NCORST steps.
C
C           Call the time step routine to calculate a new step size
C           into DTMP.
C           This is based on the solution DUP and the predicted solution
C           DTM3. When Crank Nicolson/Fractional step is active we have 
C           2nd order accuracy in time, otherwise 1st order:

            IF (TSORD.EQ.1) THEN
              CTPROJ=2D0
            ELSE
              CTPROJ=8D0
            END IF
            
            CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,0)
            CALL DTSCTR (IADTS, DADTS, 1,
     *                   NEQU, NEQP, NUVP, DWORK(L(LTM3)), DUP, 
     *                   DWORK(L(LAUX)),
     *                   RELT, RELU2, RELUM, RELP2, RELPM,
     *                   TIMEST,TIMENS,TSTEP,STPST1,STPST2,ITNSR,TSORD,
     *                   IPARAM,DPARAM,
     *                   DTMP)
            CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,1)
            
C           Calculate the relation of the previous and new step size

            DRELTN = DTMP/TSTEP

C           When the adaptive time stepping control allowes us to repeat
C           the step and the relation of the new time step to the
C           old one is too large, set IRSTP to 1 to indicate that
C           the time-step has to be repeated.
        
            IF ((ITSSRP.NE.0).AND.(DRELTN.LT.DADTS(OEPSADU))) 
     *        IRPSTP = 1

C           Some output to screen about the current status:

            IF (MT.GE.2) WRITE(MTERM,3)
            IF (MT.GE.1) WRITE(MFILE,3)
            IF (MT.GE.2) WRITE(MTERM,10002) TSTEP,
     *         RELU2/CTPROJ,RELUM/CTPROJ,RELP2/CTPROJ,RELPM/CTPROJ
            IF (MT.GE.1) WRITE(MFILE,10002) TSTEP,
     *         RELU2/CTPROJ,RELUM/CTPROJ,RELP2/CTPROJ,RELPM/CTPROJ

            IF (MT.GE.2) WRITE(MTERM,3001) IRPSTP,ITNS,DTMP,TSTEP
            IF (MT.GE.1) WRITE(MFILE,3001) IRPSTP,ITNS,DTMP,TSTEP
            IF (MT.GE.2) WRITE(MTERM,3)
            IF (MT.GE.2) WRITE(MTERM,*)
            IF (MT.GE.1) WRITE(MFILE,3)
            IF (MT.GE.1) WRITE(MFILE,*)

C           Accept the new step size

            TSTEP  = DTMP

C           Try to improve the solution further by extrapolation.
C           Check if we are allowed to do extrapolation. If yes, call
C           the extrapolation routine to carry out that task.

            IF ((IRPSTP.EQ.0).AND.(IADTS(OIEXTIM).NE.0)) THEN

C             =========================================================
C             Extrapolation step
C             =========================================================
C             Call the extrapolation routine to extrapolate a solution
C             in time. 

              CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,0)

              CALL DTEXPR (IADTS, DADTS,VECDAT(1,NLMAX),
     *                     NUVP, DWORK(L(LTM3)), DUP, 
     *                     DWORK(L(LAUX)),
     *                     TIMEST,TIMENS,TSTEP,STSTP1,STSTP2,
     *                     ITNSR,TSORD,
     *                     IPARAM,DPARAM)
     
              CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,1)
     
C             That's it. Solution completely calculated.
     
            END IF

          END IF ! ((IADTIM.NE.0).AND.(ICORST.EQ.NCORST))
          
        END IF ! (IRPSTP.EQ.0)
        
C       ================================================================
C       Error handling and repetition of the time step
C       ================================================================
C
C       The parameter IREPIT from the DAT-file (where NNITR holds a
C       local copy from) grants up to IREPIT repetitions of the
C       calculation. If IRPSTP<>1, we'd like to repeat the current 
C       time-step - as long as there are repetitions left.

        IF (IRPSTP.NE.0) THEN
        
C         Are we allowed to repeat?
C         We are allowed to repeat if
C          a) we have repetitions left and
C          b) The user has activated any kind of adaptive time stepping and
C          c1) the parameters allow us to repeat a step or
C          c2) the solver for the nonlinear iteration (predictor
C              or corrector step) broke down
C         c2) is an exception to the rule that we don't repeat
C         time steps in adaptive time stepping technique 1!
        
          IF ((ITNSR.LE.NNITR) .AND.
     *        (IADTIM.NE.0)    .AND.
     *        ((ITSSRP.NE.0).OR.(IRPSTP.EQ.2).OR.(IRPSTP.EQ.3))) THEN
          
            CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,0)
          
C           Restore the solution vector from the backup DTM0
              
            CALL LCP1(DWORK(L(LTM0)),DUP,NUVP)
            
C           And restore the RHS if we have an time-dependend RHS

            IF (LBKRHS.NE.0) THEN
              CALL LCP1(DWORK(L(LBKRHS)),DWORK(L(LRHS)),NUVP)
            ENDIF
            
            CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,1)

C           Increase repetition counter

            ITNSR=ITNSR+1

C           Reset the simulation time to where we had the backup
          
            TIMENS = BCKTNS
            ITNS   = BCKITE
            
C           Restore the "correction step counter" - will be incremented
C           on start of the recalculation of the time step.
            
            ICORST = 0
            
C           The loop will then repeat the step...

          ELSE
          
C           No, we mustn't repeat. Then at least, start with a new
C           predictor step by resetting the "correction step counter"
C           to 0.

            ICORST = 0

C           Set IRPSTP back to 0, don't repeat the step.
C           This will later be an indicator for the postprocessing
C           routine whether the step is accepted or not.
        
            IRPSTP = 0
            
          END IF

        END IF

C       Otherwise, IRPSTP is =0, so the step is not repeated.
C       This is important, as it will later be an indicator for the 
C       postprocessing routine whether the step is accepted or not.
        
C       ================================================================
C       Error analysis and postprocessing
C       ================================================================
       
C       Don't perform the error analysis if it's clear that we have to
C       repeat the step:
       
        IF (IRPSTP.EQ.0) THEN
       
          IF (MT.GE.3) WRITE (MTERM,*) 
          IF (MT.GE.0) WRITE (MFILE,*) 
          IF (MT.GE.3) WRITE (MTERM,*) 
     *                'Starting postprocessing...'
          IF (MT.GE.0) WRITE (MFILE,*) 
     *                'Starting postprocessing...'
          IF (MT.GE.3) WRITE (MTERM,*) 
          IF (MT.GE.0) WRITE (MFILE,*) 

C         Call the postprocessing routine. This will calculate errors,
C         print them on screen,write GMV files,...
C         I bet this parameter list is the longest in the whole program :)
C         No wonder - postprocessing might need information about
C         everything...
C         The variable ISTATN defines the status for the next run and must
C         be set in PREPOS.

          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,0)
          
          IUDREP = 0
          CALL PREPOS (NLMIN,NLMAX,
     *                 CURTRI,MATDAT,VECDAT,
     *                 IPARAM,DPARAM,ISTPAR,DSTPAR,
     *                 IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *                 IADTS,DADTS,IPPDAT,DPPDAT,ITRACK,
     *                 NUVP,DUP,DWORK(L(LRHS)),DWORK(L(LAUX)),
     *                 TIMEST, DASMBL(OTIMENS),TIMSTP(OTSTEP), 
     *                 ICORST, NCORST, ITNSR,STSTP1,STSTP2,1,
     *                 ISTATN,IUDREP)
     
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,1)
     
          IF (IUDREP.NE.0) THEN

C           Oops, the caller forces us to repeat the step!

            IRPSTP = 2

C           Are we allowed to repeat?
C           We are allowed to repeat if
C            a) we have repetitions left and
C            b) The user has activated any kind of adaptive time stepping and
C            c) the parameters allow us to repeat a step
        
            IF ((ITNSR.LE.NNITR) .AND.
     *          (IADTIM.NE.0)    .AND.
     *          ((ITSSRP.NE.0).OR.(IRPSTP.EQ.2))) THEN
            
              CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,0)
            
C             Restore the solution vector from the backup DTM0
                
              CALL LCP1(DWORK(L(LTM0)),DUP,NUVP)
              
C             And restore the RHS if we have an time-dependend RHS

              IF (LBKRHS.NE.0) THEN
                CALL LCP1(DWORK(L(LBKRHS)),DWORK(L(LRHS)),NUVP)
              ENDIF
              
              CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,1)

C             Increase repetition counter

              ITNSR=ITNSR+1

C             Reset the simulation time to where we had the backup
            
              TIMENS = BCKTNS
              ITNS   = BCKITE
              
C             Restore the "correction step counter" - will be incremented
C             on start of the recalculation of the time step.
              
              ICORST = 0
              
C             The loop will then repeat the step...

            END IF
            
          END IF ! IUDREP <> 0
     
          IF (IRPSTP.EQ.0) THEN
     
C           Additionally check, if the solution is accepted.
C           If that's the case, call the postprocessing routine again
C           with an indicator whether the solution is accepted
C           or not. A solution counts as accepted if
C            a) We are in the last substep of the time stepping scheme and
C            b) The step should not be repeated
C           ICORST=0 happens if the step is forcefully accepted because
C           we don't have repetitions left.

            IF ((ICORST.EQ.0).OR.(ICORST.EQ.NCORST)) THEN
            
              CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,0)
              IUDREP = 0
              CALL PREPOS (NLMIN,NLMAX,
     *                   CURTRI,MATDAT,VECDAT,
     *                   IPARAM,DPARAM,ISTPAR,DSTPAR,
     *                   IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *                   IADTS,DADTS,IPPDAT,DPPDAT,ITRACK,
     *                   NUVP,DUP,DWORK(L(LRHS)),
     *                   DWORK(L(LAUX)),
     *                   TIMEST, DASMBL(OTIMENS),TIMSTP(OTSTEP), 
     *                   ICORST, NCORST,ITNSR,STSTP1,STSTP2,0,
     *                   ISTATN,IUDREP)
              CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,1)

C             Furthermore, as the solution is accepted,            
C             set the repetition counter back to 0. We now have again
C             NNITR repetitions left for the next time step.

              ITNSR = 0

            END IF
            
          END IF ! IRPSTP = 0
          
C         Stationary solution convergence criterion.
C
C         In the next step we compute the solution of the previous time
C         step with our current solution to check if the time dependent
C         simulation got stationary. If yes, we can finish.
C
C         Calculate the time derivative into RELT. Use DTM0 (the solution
C         of the beginning of the time step) and DUP for that purpose.
C         Use the information about the simulation time to calculate
C         their distance in time.

          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,0)
          CALL DTCHKS (IADTS, DADTS, 
     *                 NEQU, NEQP, NUVP, DWORK(L(LTM0)), DUP, 
     *                 DWORK(L(LAUX)),
     *                 TIMENS-BCKTNS, IPARAM, DPARAM,
     *                 RELT, RELU2, RELUM, RELP2, RELPM)
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,1)

C         Print out the results to screen

          IF (MT.GE.2) WRITE(MTERM,2)
          IF (MT.GE.2) WRITE(MTERM,20001) ITNS,ICORST,
     *                                    TIMENS,RELU2,RELP2,RELT
          IF (MT.GE.2) WRITE(MTERM,2)
          IF (MT.GE.2) WRITE(MTERM,*)

          IF (MT.GE.1) WRITE(MFILE,2)
          IF (MT.GE.1) WRITE(MFILE,20001) ITNS,ICORST,
     *                                    TIMENS,RELU2,RELP2,RELT
          IF (MT.GE.1) WRITE(MFILE,2)
          IF (MT.GE.1) WRITE(MFILE,*)

C         Stop the calculation if the time-derivative is too low
C         (-> simulation got stationary), or if the maximum simulation
C         time is reached.

          IF (RELT.LE.DADTS(OEPSNS)) ICONV = 1

        ELSE ! (IRPSTP.EQ.0)
        
C         Solution is rejected. Inform the pre/postprocessing
C         routine.

          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,0)
          IUDREP = 0
          CALL PREPOS (NLMIN,NLMAX,
     *                 CURTRI,MATDAT,VECDAT,
     *                 IPARAM,DPARAM,ISTPAR,DSTPAR,
     *                 IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *                 IADTS,DADTS,IPPDAT,DPPDAT,ITRACK,
     *                 NUVP,DUP,DWORK(L(LRHS)),DWORK(L(LAUX)),
     *                 TIMEST, DASMBL(OTIMENS),TIMSTP(OTSTEP),
     *                 ICORST, NCORST, ITNSR,STSTP1,STSTP2,5,
     *                 ISTATN,IUDREP)
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,1)

        END IF ! (IRPSTP.EQ.0)

C       When the predictor step is not the reason for a repetition,
C       print out statistical data.
C       If it is the reason for a repetition, no corrector step
C       is performed and thus we would print the statistics twice...

        IF (IRPSTP.NE.3) THEN

C         Finally some nice output to the user.
C
C         Stop the time. 
      
          CALL GTMAUX (DTIMIN,DPARAM,OTMTOT,1)
          CALL GTMAUX (DTIMIN,DPARAM,OTMTOT,0)                
          
          IF (MT.GE.3) WRITE (MTERM,*) 
     *                'Timing statistics - Total/Last step'
          IF (MT.GE.0) WRITE (MFILE,*) 
     *                'Timing statistics - Total/Last step'
          IF (MT.GE.3) WRITE (MTERM,*) 
     *                '-----------------------------------'
          IF (MT.GE.0) WRITE (MFILE,*) 
     *                '-----------------------------------'
                                                       
          IF (MT.GE.3) WRITE(MTERM,*) 'total time : ', 
     *      DPARAM(OTMTOT),'/',
     *      DPARAM(OTMTOT)-DTIMBK(OTMTOT)
          IF (MT.GE.0) WRITE(MFILE,*) 'total time : ', 
     *      DPARAM(OTMTOT),'/',
     *      DPARAM(OTMTOT)-DTIMBK(OTMTOT)
          IF (MT.GE.3) WRITE(MTERM,*) 'mavec time : ', 
     *      DPARAM(OTNSTIM-1+OTTADF),'/',
     *      DPARAM(OTNSTIM-1+OTTADF)-DTIMBK(OTNSTIM-1+OTTADF)
          IF (MT.GE.0) WRITE(MFILE,*) 'mavec time : ', 
     *      DPARAM(OTNSTIM-1+OTTADF),'/',
     *      DPARAM(OTNSTIM-1+OTTADF)-DTIMBK(OTNSTIM-1+OTTADF)
          IF (MT.GE.3) WRITE(MTERM,*) 'konv. time : ', 
     *      DPARAM(OTNSTIM-1+OTTUPW),'/',
     *      DPARAM(OTNSTIM-1+OTTUPW)-DTIMBK(OTNSTIM-1+OTTUPW)
          IF (MT.GE.0) WRITE(MFILE,*) 'konv. time : ', 
     *      DPARAM(OTNSTIM-1+OTTUPW),'/',
     *      DPARAM(OTNSTIM-1+OTTUPW)-DTIMBK(OTNSTIM-1+OTTUPW)
          IF (MT.GE.3) WRITE(MTERM,*) 'bdry  time : ', 
     *      DPARAM(OTNSTIM-1+OTTBDR),'/',
     *      DPARAM(OTNSTIM-1+OTTBDR)-DTIMBK(OTNSTIM-1+OTTBDR)
          IF (MT.GE.0) WRITE(MFILE,*) 'bdry  time : ', 
     *      DPARAM(OTNSTIM-1+OTTBDR),'/',
     *      DPARAM(OTNSTIM-1+OTTBDR)-DTIMBK(OTNSTIM-1+OTTBDR)
          IF (MT.GE.3) WRITE(MTERM,*) 'LC    time : ', 
     *      DPARAM(OTNSTIM-1+OTTLC),'/',
     *      DPARAM(OTNSTIM-1+OTTLC)-DTIMBK(OTNSTIM-1+OTTLC)
          IF (MT.GE.0) WRITE(MFILE,*) 'LC    time : ', 
     *      DPARAM(OTNSTIM-1+OTTLC),'/',
     *      DPARAM(OTNSTIM-1+OTTLC)-DTIMBK(OTNSTIM-1+OTTLC)
          IF (MT.GE.3) WRITE(MTERM,*) 'LinSl time : ', 
     *      DPARAM(OTNSTIM-1+OTTLSOL),'/',
     *      DPARAM(OTNSTIM-1+OTTLSOL)-DTIMBK(OTNSTIM-1+OTTLSOL)
          IF (MT.GE.0) WRITE(MFILE,*) 'LinSl time : ', 
     *      DPARAM(OTNSTIM-1+OTTLSOL),'/',
     *      DPARAM(OTNSTIM-1+OTTLSOL)-DTIMBK(OTNSTIM-1+OTTLSOL)
     
          IF (MT.GE.3) WRITE (MTERM,*) 
          IF (MT.GE.0) WRITE (MFILE,*) 

C         Back up the timing statistics for the next step.
C         Actually we make a backup of the whole parametrr block, but
C         it contains the timing block and that's enough :)

          CALL LCP1(DPARAM,DTIMBK,SZISDD)

        END IF ! (IRPSTP.NE.2)
     
C       That's it with the time-dependent, adaptive time-stepping
C       simulation loop. Proceed with the next timestep.

      END DO ! WHILE

C     At this point the time-dependent simulation is completed.
C     Only the time-dependent simulation will reach this point,
C     as the stationary simulation jumps out of the caluculation
C     above!
C     Reset ITNS to NITNS as at the end of the simulation,
C     the DO-loop brings ITNS to NITNS+1...

      ITNS=NITNS
      
C     User defined postprocessing at the end

      CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,0)
      IUDREP = 0
      CALL PREPOS (NLMIN,NLMAX,
     *           CURTRI,MATDAT,VECDAT,
     *           IPARAM,DPARAM,ISTPAR,DSTPAR,
     *           IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *           IADTS,DADTS,IPPDAT,DPPDAT,ITRACK,
     *           NUVP,DWORK(L(LTM3)),
     *           DWORK(L(LRHS)),DWORK(L(LAUX)),
     *           TIMEST, DASMBL(OTIMENS),0D0, ICORST, NCORST,
     *           ITNSR,0,0,20,
     *           ISTATN,IUDREP)
      CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,1)
      
C     Release boundary information

      CALL DNBDGE (NLMIN,NLMAX,
     *             CURTRI,IASMBL,DASMBL)
     
C     Delete copies of the triangulation.
C     Node that this will only delete information that is maintained
C     by the CURTRI structures! Information that is shared between
C     CURTRI and TRIAS will not be deleted!

      DO I=NLMAX,NLMIN,-1
        CALL TRIDEL (CURTRI(1,I))
      END DO
      
C     Release grid velocity vector if ALE method was used:

      IF (IASMBL(OIALE).NE.0) THEN
        CALL ZDISP(0,IASMBL(OLGRVEL),'GRVEL ')
        IASMBL(OIALE) = 0
      END IF
      
C     Release the memory used.

      CALL ZDISP(0,LTM0  ,'DMT0  ')
      IF (IER.NE.0) GOTO 99998

C     Don't release memory referred to by LBKRHS and LTM3 if adaptive
C     time stepping is not active - because in this case these
C     handles are dummy handles! (see above)

      IF (IADTIM.NE.0) THEN
        CALL ZDISP(0,LBKRHS  ,'DRHS  ')
        IF (IER.NE.0) GOTO 99998
        CALL ZDISP(0,LTM3  ,'DMT3  ')
        IF (IER.NE.0) GOTO 99998
      ENDIF

C     Release auxiliary vector

      CALL ZDISP(0,LAUX,'DAUX  ')

C     Release main RHS vector

      CALL ZDISP(0,LRHS,'DRHS  ')
      
999   GOTO 99999

C     ==================================================================
C     Error case
C     ==================================================================

99998 WRITE(MTERM,*) 'IER', IER
      WRITE(MTERM,*) 'IN SUBROUTINE ',SUB

99999 CONTINUE

C     Stop the final time.
      
      CALL GTMAUX (DTIMIN,DPARAM,OTMTOT,1)
      
C     That's it - nonstationary solver finished.

   1  FORMAT(79('-'))
   2  FORMAT(79('='))
   3  FORMAT(79('@'))
1000  FORMAT (6E12.5)
1001  FORMAT(' IT DIV-L2',3X,'RHOMGP')
1003  FORMAT(I3,2(D9.2))
2001  FORMAT('Predictor-Step at Time-Step ',I4,', Time =',D10.3,
     *       ', Stepsize: DT1 = ',D10.3)
2002  FORMAT('Time-Step ',I4,' at Time = ',D10.3,
     *       ' with Stepsize : DT3 = ',D10.3)
3001  FORMAT('Timestepping by ',I2,'(',I4,'), New Stepsize = ',
     *       D10.3,', Old Stepsize = ',D10.3)
10002 FORMAT ('OLD DT=',D9.2,
     *        1X,'U(L2)=',D9.2,1X,'U(MX)=',D9.2,
     *        1X,'P(L2)=',D9.2,1X,'P(MX)=',D9.2)
10003 FORMAT ('OLD DT=',D9.2,
     *        1X,'FWREL=',D9.2,1X,'AWREL=',D9.2,
     *        1X,'FW3  =',D9.2,1X,'FA3  =',D9.2)
20001 FORMAT ('#',I4,1X,'(',I4,')',1X,'TIME=',D10.3,1X,'RELU(L2)=',
     *        D9.2,1X,'RELP(L2)=',D9.2,1X,'REL=',D9.2)

      END

************************************************************************
* Default Time step control and error analysis
*
* This routine calculates, based on a predictor and a corrector step
* solution, the errors in U and P and tries to find a new suitable
* time step size.
*
* In:
*   IADTS  : array [1..SZADTI] of integer
*   DADTS  : array [1..SZADTD] of double
*            Integer and double prec. parameter block that configures
*            the behaviour of the adaptive time stepping.
*   ITANL  - Type of time error control with the solution vectors
*            =0: Don't perform any time analysis
*                This is used to indicate that one of the solvers
*                broke down, so time analysis is not possible.
*                Calculate a new step size without respecting
*                time information.
*            =1: Use time analysis, calculate RELT, RELU2,...
*   NEQU   - Length of velocity component in solution vector
*   NEQP   - Length of pressure component in solution vector
*   NUVP   - Total length of solution vector; usually 2*NEQU+NEQP
*   DUPPRE - array [1..NUVP] of double
*            Solution vector after predictor step
*   DUP    - array [1..NUVP] of double
*            Real Solution vector
*   DAUX   - array [1..NUVP] of double
*            Auxiliary vector
*   TIMEST - Initial simulation time
*   TIMENS - Current simulation time
*   TSTEP  - Current (theoretical) time step size.
*            Can vary from the actual time step size, depending on
*            the time stepping scheme.
*   STPRED - Error indicator of the solver in the prediction step.
*            Bitfield. Can be set to 0 if not used. 
*            Bit0: =1: nonlinear solver broke down
*            Bit1: =1: linear solver broke down
*   STCORR - Error indicator of the solver in the correction step.
*            Bitfield. Can be set to 0 if not used.
*            Bit0: =1: nonlinear solver broke down
*            Bit1: =1: linear solver broke down
*   IREPET - =0: DUP was calculated for the first time
*            >0: DUP was calculated while repeating the time step
*   TAPORD - Order of the time approximation
*            =1: time approximation is of 1st order(Euler)
*            =2: time approximation is of 2nd order
*                (Crank Nicolson, Fractional Step)
*
*   IPARAM : array [1..SZISDI] of integer
*   DPARAM : array [1..SZISDI] of double
*            Integer and double prec. parameter blocks that define the
*            behaviour of the nonstationary solver. 
*
* Out:
*   TSTEPN - New time step size 
*
* If no solver broke down, or if at most the nonlinear solver in the
* prediction step broke down:
*
*   RELT   - Time derivative, depending on RELU2,RELUM,RELP2,RELPM
*
*   RELU2  - Relative L2 error in U
*   RELUM  - Maximum error in U
*   RELP2  - Relative L2 error in P
*   RELPM  - Maximum error in P
*
* If not computed, the appropriate result is set to 0D0.
************************************************************************

      SUBROUTINE DTSCTR (IADTS, DADTS, ITANL,
     *                   NEQU, NEQP, NUVP, DUPPRE, DUP, DAUX,
     *                   RELT, RELU2, RELUM, RELP2, RELPM,
     *                   TIMEST,TIMENS,TSTEP,STPRED,STCORR,IREPET,
     *                   TAPORD,IPARAM,DPARAM,
     *                   TSTEPN)
     
      IMPLICIT NONE
      
      INCLUDE 'ststepping.inc'
      INCLUDE 'sadtstep.inc'

C     parameters      

      INTEGER ITANL, NEQU, NEQP, TAPORD,IREPET, NUVP
      DOUBLE PRECISION DUPPRE(*), DUP(*), DAUX(*),TIMENS,TIMEST
      DOUBLE PRECISION RELT,RELU2, RELUM, RELP2, RELPM ,TSTEP,TSTEPN
      INTEGER STPRED,STCORR
      INTEGER IPARAM(*),IADTS(*)
      DOUBLE PRECISION DPARAM(*),DADTS(*)
      
C     local variables
      
      INTEGER IADTIM,IEPSAD,IUP,INDMAX
      DOUBLE PRECISION RELU20,RELP20,RELUM0,RELPM0,DSXN
      DOUBLE PRECISION DTFACT,DTMIN,DTMAX,EPSAD,HTFACT
      
C     The implementation follows p. 160ff, Turek's CFD-book.
      
      IADTIM = IADTS(OIADTIM)
      DTFACT = DADTS(ODTFACT)
      DTMIN = DADTS(ODTMIN)
      DTMAX = DADTS(ODTMAX)
      IEPSAD = IADTS(OIEPSAD)
      
C     As standard, take the old stepsize as the new
      
      TSTEPN = TSTEP

C     Check if we can use time analysis...
      
      IF (ITANL.EQ.0) THEN
 
C       Ehm, we really don't have much information in this case.
C       We can just "guess" a new time step - if we are at all allowed
C       to change it!

        IF (IADTIM.NE.0) THEN

          IF (IOR(IAND(STPRED,2),IAND(STCORR,1)).NE.0) THEN
          
C           If one the linear solvers of the prediction step or
C           the nonlinear solver of the correction step broke down,
C           broke down, reduce the time step by sqrt(DTFACT)
          
            IF (IADTIM.NE.0) TSTEPN = TSTEP / SQRT(DTFACT)
          
          ELSE IF (IAND(STCORR,2).NE.0) THEN

C           If one the linear solvers of the correction step broke down,
C           broke down, reduce the time step by DTFACT
          
            IF (IADTIM.NE.0) TSTEPN = TSTEP / DTFACT
          
          END IF
          
        END IF !  (IADTIM.NE.0)

C       Otherwise we leave the time step as it is...
      
      ELSE
      
C       Use time analysis to calculate a step size.

        RELU2=0D0
        RELUM=0D0
        RELP2=0D0
        RELPM=0D0

C       DUP contains the result of the calculation with the standard
C       time-steps, DUPPRE the result of the predictor-step.
C       Subtract both solutions and save the result in DAUX

        DO IUP=1,NUVP
          DAUX(IUP) = DUP(IUP)-DUPPRE(IUP)
        END DO

C       Calculate some residual information. This will help us later
C       to calculate the time-step.

        IF ((IEPSAD.GE.1).OR.(IEPSAD.EQ.-1)) THEN
          CALL LL21(DUP,2*NEQU,DSXN)
          RELU20=DSXN/(SQRT(DBLE(2*NEQU)))
          IF (RELU20.LE.1D0) RELU20=1D0
          CALL LL21(DAUX,2*NEQU,DSXN)
          RELU2=DSXN/(SQRT(DBLE(2*NEQU))*RELU20)
        ENDIF

        IF ((IEPSAD.GE.1).OR.(IEPSAD.EQ.-2)) THEN
          CALL LLI1(DUP,2*NEQU,DSXN,INDMAX)
          RELUM0=DSXN
          IF (RELUM0.LE.1D0) RELUM0=1D0
          CALL LLI1(DAUX,2*NEQU,DSXN,INDMAX)
          RELUM=DSXN/RELUM0
        ENDIF

        IF ((IEPSAD.GE.1).OR.(IEPSAD.EQ.-3)) THEN
          CALL LL21(DUP(1+2*NEQU),NEQP,DSXN)
          RELP20=DSXN/(SQRT(DBLE(NEQP)))
          IF (RELP20.LE.1D0) RELP20=1D0
          CALL LL21(DAUX(1+2*NEQU),NEQP,DSXN)
          RELP2=DSXN/(SQRT(DBLE(NEQP))*RELP20)
        ENDIF

        IF ((IEPSAD.GE.1).OR.(IEPSAD.EQ.-4)) THEN
          CALL LLI1(DUP(1+2*NEQU),NEQP,DSXN,INDMAX)
          RELPM0=DSXN
          IF (RELPM0.LE.1D0) RELPM0=1D0
          CALL LLI1(DAUX(1+2*NEQU),NEQP,DSXN,INDMAX)
          RELPM=DSXN/RELPM0
        ENDIF
        
C       Depending on IEPSAD, set up RELT=J(.), our error functional:
        
        IF (ABS(IEPSAD).EQ.1)  RELT=RELU2
        IF (ABS(IEPSAD).EQ.2)  RELT=RELUM
        IF (ABS(IEPSAD).EQ.3)  RELT=RELP2
        IF (ABS(IEPSAD).EQ.4)  RELT=RELPM
        IF (     IEPSAD.EQ.5)  RELT=MAX(RELU2,RELP2)
        IF (     IEPSAD.EQ.6)  RELT=MAX(RELUM,RELPM)
        IF (     IEPSAD.EQ.7)  RELT=MAX(RELU2,RELUM,RELP2,RELPM)
        IF (     IEPSAD.EQ.8)  RELT=MIN(RELU2,RELUM,RELP2,RELPM)
        
C       Now, should we compute a new step size?

        IF (IADTIM.NE.0) THEN
        
C         Using the values in TIMExx, EPSADx and IADIN, adjust (i.e.
C         calculate) a new stopping criterion EPSAD.
            
          CALL CRITAD(DADTS(OTIMEIN),TIMENS,TIMEST,
     *                DADTS(OEPSADI),DADTS(OEPSADL),EPSAD,
     *                IADTS(OIADIN))

          IF (TAPORD.LE.1) THEN

C           Time approximation of 1st order; used for all one-step
C           schemes. The error can be represented as:
C             J(v_k) - J(v) = k e(v) + O(k^2)

            TSTEPN=TSTEP*2D0*EPSAD/RELT
            
          ELSE
          
C           Time approximation of 1st order; used for Fractional
C           step Theta scheme. The error can be represented as:
C             J(v_k) - J(v) = k^2 e(v) + O(k^4)

            TSTEPN=TSTEP*SQRT(8D0*EPSAD/RELT)
            
          ENDIF

C         The nonlinear solver in the prediction step might break
C         down, while the other solvers work.
C         When the nonlinear solver broke down and we are in a time
C         stepping mode that allowes repetition of the time step,
C         bound the time step: TSTEPN must result in the interval
C         [ TSTEP/SQRT(DTFACT) .. TSTEP ]:

          IF ((ABS(IADTIM).GT.1).AND.(IAND(STPRED,1).NE.0)) THEN
            TSTEPN = MAX( TSTEP/SQRT(DTFACT), MIN(TSTEPN,TSTEP) )
          ENDIF

C         TSTEPN must not be smaller than TSTEP/DTFACT

          TSTEPN = MAX( TSTEPN, TSTEP/DTFACT)

C         Use another upper bound if the solution was calculated
C         by repetition of a time step:

          IF (IREPET.GT.0) THEN
            HTFACT=DTFACT**(1D0/DBLE(IREPET+1))
            TSTEPN = MIN(TSTEPN,TSTEP*HTFACT) 
          ELSE       
            TSTEPN = MIN(TSTEPN,TSTEP*DTFACT) 
          ENDIF       
          
        END IF ! (IADTIM.NE.0)
        
      END IF ! (IAND(STPRED,2).NE.0)
      
C     Bound the time step to the interval DTMIN..DTMAX - for sure
          
      TSTEPN = MIN(DTMAX, MAX(TSTEPN, DTMIN))

      END
      
************************************************************************
* Default Time extrapolation
*
* This routine extrapolates the solution of a predictor step with
* the solution of a corrector step to generate a new (better) solution.
*
* In:
*   IADTS  - array [1..SZADTI] of integer
*   DADTS  - array [1..SZADTD] of double
*            Integer and double prec. parameter block that configures
*            the behaviour of the adaptive time stepping.
*   VECDAT - array [1..SZN2VI] of integer 
*            TNS2DVectorParams-structure, 
*            Defines the form of the RHS vector.
*   NUVP   - Length of the vectors
*   DUPPRE - array [1..NUVP] of double
*            Solution vector after predictor step
*   DUP    - array [1..NUVP] of double
*            Real Solution vector
*   DAUX   - array [1..NUVP] of double
*            Auxiliary vector
*   TIMEST - Initial simulation time
*   TIMENS - Current simulation time
*   TSTEP  - Current (theoretical) time step size.
*            Can vary from the actual time step size, depending on
*            the time stepping scheme.
*   STPRED - Error indicator of the solver in the prediction step.
*            Bitfield. Can be set to 0 if not used. 
*            Bit0: =1: nonlinear solver broke down
*            Bit1: =1: linear solver broke down
*   STCORR - Error indicator of the solver in the correction step.
*            Bitfield. Can be set to 0 if not used.
*            Bit0: =1: nonlinear solver broke down
*            Bit1: =1: linear solver broke down
*   IREPET - =0: DUP was calculated for the first time
*            >0: DUP was calculated while repeating the time step
*   TAPORD - Order of the time approximation
*            =1: time approximation is of 1st order(Euler)
*            =2: time approximation is of 2nd order
*                (Crank Nicolson, Fractional Step)
*
* Out:
*   DUP    - New solution vector
************************************************************************

      SUBROUTINE DTEXPR (IADTS, DADTS,
     *                   VECDAT, NUVP, DUPPRE, DUP, DAUX,
     *                   TIMEST,TIMENS,TSTEP,STPRED,STCORR,
     *                   IREPET,TAPORD,
     *                   IPARAM,DPARAM)
      
      IMPLICIT NONE

      INCLUDE 'sadtstep.inc'

      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'

C     parameters
      
      INTEGER TAPORD,STPRED,STCORR,IADTS(*),IREPET,NUVP
      DOUBLE PRECISION DUPPRE(*),DUP(*),TIMENS,TIMEST,TSTEP,DADTS(*)
      INTEGER IPARAM(*),VECDAT(SZN2VI)
      DOUBLE PRECISION DPARAM(*),DAUX(*)

C     local variables

      INTEGER IADTIM,IEXTIM
      DOUBLE PRECISION EX1,EX2

      IADTIM = IADTS(OIADTIM)
      IEXTIM = IADTS(OIEXTIM)
      
      IF (TAPORD.EQ.1) THEN

C       1st order extrapolation
      
        EX1= 4D0/3D0
        EX2=-1D0/3D0
        
      ELSE
      
C       2nd order extrapolation
      
        EX1= 9D0/8D0
        EX2=-1D0/8D0
        
      END IF

      IF ((IEXTIM.NE.0).AND.(IADTIM.NE.0).AND.(STPRED.EQ.0)) THEN
      
C       Calculate the new solution vector by weighted combination
C       of the standard-time-step solution DU and the predictor-step
C       in DUPPRE. 

        CALL LLC1(DUPPRE,DUP,NUVP,EX2,EX1)
      
      ENDIF
      
      END
      
************************************************************************
* Default Check for stationary solution
*
* This routine compares the solution of the current time step with
* the previous solution to calculate the time derivative.
* This is used to figure out i fthe solution got stationary.
*
* In:
*   IADTS  : array [1..SZADTI] of integer
*   DADTS  : array [1..SZADTD] of double
*            Integer and double prec. parameter block that configures
*            the behaviour of the adaptive time stepping.
*   NEQU   - Length of velocity component in solution vector
*   NEQP   - Length of pressure component in solution vector
*   DUPOLD - array [1..*] of double
*            Solution vector of beginning of the time step
*   DUP    - array [1..*] of double
*            Solution vector after the time step
*   DAUX   - array [1..*] of double
*            Auxiliary vector
*   TSTEP  - Length of the time step.
*
*   IPARAM : array [1..SZISDI] of integer
*   DPARAM : array [1..SZISDI] of double
*            Integer and double prec. parameter blocks that define the
*            behaviour of the nonstationary solver. 
*
* Out:
*   RELT   - Time derivative, depending on RELU2,RELUM,RELP2,RELPM
*
*   RELU2  - Relative L2 error in U
*   RELUM  - Maximum error in U
*   RELP2  - Relative L2 error in P
*   RELPM  - Maximum error in P
************************************************************************

      SUBROUTINE DTCHKS (IADTS, DADTS, 
     *                   NEQU, NEQP, NUVP, DUPOLD, DUP, DAUX,
     *                   TSTEP, IPARAM, DPARAM,
     *                   RELT, RELU2, RELUM, RELP2, RELPM)
     
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      
      INCLUDE 'ststepping.inc'
      INCLUDE 'sadtstep.inc'

C     parameters      

      INTEGER NEQU, NEQP,NUVP
      DOUBLE PRECISION DUPOLD(*), DUP(*), DAUX(*)
      DOUBLE PRECISION RELT, RELU2, RELUM, RELP2, RELPM, TSTEP
      INTEGER IPARAM(*),IADTS(*)
      DOUBLE PRECISION DPARAM(*),DADTS(*)
      
C     local variables
      
      INTEGER IEPSAD,INDMAX
      DOUBLE PRECISION DSXN
      
      IEPSAD = IADTS(OIEPSAD)

C     Calculate the difference vector in DAUX

      CALL LCP1(DUPOLD,DAUX,NUVP)
      CALL LLC1(DUP,DAUX,NUVP,-1D0,1D0)

C     Calculate some residual information for output

      RELU2=0D0
      RELP2=0D0
      RELUM=0D0
      RELPM=0D0

      CALL LL21(DAUX,2*NEQU,DSXN)
      RELU2=DSXN/(SQRT(DBLE(2*NEQU))*TSTEP)
      CALL LL21(DAUX(1+2*NEQU),NEQP,DSXN)
      RELP2=DSXN/(SQRT(DBLE(NEQP))*TSTEP)

      IF ((ABS(IEPSAD).EQ.2).OR.(IEPSAD.GE.5)) THEN
        CALL LLI1(DAUX,2*NEQU,DSXN,INDMAX)
        RELUM=DSXN/TSTEP
      ENDIF

      IF ((ABS(IEPSAD).EQ.4).OR.(IEPSAD.GE.5)) THEN
        CALL LLI1(DAUX(1+2*NEQU),NEQP,DSXN,INDMAX)
        RELPM=DSXN/TSTEP
      ENDIF
      
      IF (ABS(IEPSAD).EQ.1)  RELT=RELU2
      IF (ABS(IEPSAD).EQ.2)  RELT=RELUM
      IF (ABS(IEPSAD).EQ.3)  RELT=RELP2
      IF (ABS(IEPSAD).EQ.4)  RELT=RELPM
      IF (    IEPSAD .EQ.5)  RELT=MAX(RELU2,RELP2)
      IF (    IEPSAD .EQ.6)  RELT=MAX(RELUM,RELPM)
      IF (    IEPSAD .EQ.7)  RELT=MAX(RELU2,RELUM,RELP2,RELPM)
      IF (    IEPSAD .EQ.8)  RELT=MIN(RELU2,RELUM,RELP2,RELPM)

      END
      