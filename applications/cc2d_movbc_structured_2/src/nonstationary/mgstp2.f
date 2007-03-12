************************************************************************
* The routines in this file are responsible for performing
* single macrosteps for the instationary Navier Stokes equations.
*
* The stationary solver NSDEF2 is used as sub-solver.
************************************************************************
      
      SUBROUTINE MGSTP2(TRIAS,MATDAT,VECDAT,NLMIN,NLMAX,
     *                  IPARAM,DPARAM,
     *                  ISTPAR,DSTPAR,
     *                  IMGPAR,DMGPAR,
     *                  IASMBL,DASMBL,IGEOM,DGEOM,
     *                  NUVP,DUP,DRHS,DAUX,
     *                  TIMSTP,CURTRI,
     *                  NNL,NMG,ISTATN,STATUS,
     *                  GENRHS,GENGRI,PRCNTS) 
      
************************************************************************
* Solver for one time step of the coupled method for
* solving the (in)stationary incompressible Navier Stokes.
*
* Performs 1 time step of coupled method.
*
* Extended calling convention.
*
* This routine calculates the solution u^(n+1) for the next time step
* based on the solution u^n in the current time step. It builds up
* right-hand-side vectors according to the time stepping and then
* calls the stationary solver to calculate the solution.
*
* The routine is configured by the information in the double/integer
* parameter blocks IPARAM/DPARAM. Parameters of the nonstationary
* solver as well as parameters for the nonlinear solver are needed
* to adapt to the problem.
* IMGPAR/DMGPAR that defines the behaviour of the multigrid
* solution engine for the linear subproblems. Apart from that,
* a parameter block TRIAS must define the underlying triangulations
* on level NLMIN..NLMAX, MATDAT must provide information about
* the system matrices and VECDAT must provide space for the
* solution of the linear systems with multigrid.
*
* In:
*   TRIAS  : array [1..SZTRIA,1..NLEV] of integer
*            Triangulation structures for the underlying meshes
*            on all levels. TRIAS(.,ILEV) represents the triangulation
*            structure on level ILEV in its ummodified form, i.e.
*            without any boundary conditions implemented, probably
*            undeformed (if the user wants to start the grid
*            deformation in every time step from the original grids).
*   MATDAT : array [1..SZN2MI,1..NNLEV] of integer
*            TNS2DMatrixParams-structures for level NLMIN..NLMAX, 
*            initialised with data
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
*   NLMIN  : minimum level in MATDAT/VECDAT that is filled with data
*   NLMAX  : maximum level in MATDAT/VECDAT that is filled with data
*
*   IPARAM : array [1..SZISDI] of integer
*   DPARAM : array [1..SZISDD] of double
*            Integer and double prec. parameter blocks that define the
*            behaviour of the nonstationary solver. 
*   ISTPAR : array [1..SZNSDI] of integer
*   DSTPAR : array [1..SZNSDD] of double
*            Integer and double prec. parameter block for the stationary 
*            sub-solver NSDEF2.
*            The input-variables of the array must be initialised, e.g.
*            with INISTS.
*            The structures will be modified problem-specifically for the
*            solution of the subproblems, but restored onto return
*            of the routine.
*   IMGPAR : array [1..SZ020I+3*SZSLVI] of integer
*   DMGPAR : array [1..SZ020D+3*SZSLVD] of double
*            Integer and double parameter blocks for the multigrid 
*            sub-solver M020. These must be initialised as described
*            in NSDEF2.F. The variables in these blocks are not used
*            in MGSTP2, they are directly passed to NSDEF2.
*   IASMBL : array [1..SZASMI] of integer
*   DASMBL : array [1..SZASMD] of double
*            Integer and double prec. parameter block that controls the
*            discretisation. This tells all assembly-routines how to 
*            set up the nonlinearity in the system matrices, which 
*            cubature formula to use, etc.
*            The discretisation structures are modified during the
*            time stepping according to the situation of the simulation,
*            but restored to their original state upon returning of
*            this routine.
*   IGEOM  : array [1..*] of integer 
*   DGEOM  : array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to fictitious boundary
*            routines. Not used in this routine.
*
*   NUVP   : total length of solution vector; usually 2*NEQU+NEQP
*   DUP    : array [1..NUVP] of double
*            Start vector for nonlinear iteration, DUP=(u,v,p)
*            = u^n in a time-dependent simulation.
*            Remark: DUP must not coincide with VECDAT[NLMAX].LSOL !
*   DRHS   : array [1..NUVP] of double
*            Right hand side for the nonlinear equation.
*            Remark: DRHS must not coincide with VECDAT[NLMAX].LRHS !
*            When the RHS is time-dependent (IRHS=2) or the Neumann
*            boundary conditions are time-dependent (IPARAM.IBDR<>0), the
*            boundary conditions must not be implemented into the RHS.
*            Otherwise the boundary conditions must be implemented
*            by the caller.
*   DAUX   : array [1..NUVP] of double
*            Auxiliary vector.
*            Remark: DAUX must not coincide with VECDAT[NLMAX].LRHS !
*   TIMSTP : array [1..SZTSTD] of double
*            Time-step structure that defines current simulation time and
*            the weights that should be used for the RHS/matrix to 
*            calculate u^n+1 (for nonstationary simulation) or u^1 
*            (for stationary simulation), respectively.
*
*   CURTRI : array [1..SZTRIA,1..NNLEV] of integer
*            Triangulation structures for the underlying meshes
*            on all levels. CURTRI(.,ILEV) represents the triangulation
*            structure on level ILEV in the current time step.
*            Remarks:
*            1.) CURTRI must point to the IPARAM.CURTRI variable in the
*             parameter structure! This is for easier and faster access.
*            2.) CURTRI is assumed to be a duplicate from the TRIAS
*             structure set, created by TRIDUP.
*            3.) For time-dependent boundary conditions or moving meshes,
*             these structures will be modified to include boundary
*             conditions. For moving meshes, the new mesh is generated
*             into these structures.
*            4.) The boundary conditions saved in the triangulation
*             structures must be up-to-date for the current simulation
*             time. This can be done by a call to UPDBDX for all 
*             triangulations in CURTRI before calling this routine.
*
*   ISTATN : Status flag that tells MGSTP which tasks must be performed
*            in the preparation and which can be skipped. Bitfield.
*            Bit0: call mesh adaption in the preparation phase
*
*   GENRHS : Callback routine that generates a RHS for a time step.
*            SUBROUTINE GNRHSV (TRIA,IPARAM,DPARAM,
*                               IASMBL,DASMBL,IGEOM,DGEOM,VECDAT,DRHS)
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
*                               IGEOM,DGEOM,
*                               VECDAT,DUP,DRHS,DAUX,
*                               TIMEST, TIMENS, TSTEP)
*
*   PRCNTS : Proceed to next time step.
*            This routine is called directly before the time step
*            increases from TIMENS to TIMENS+TSTEP. It allowes the
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
* Out:
*   STATUS : = 0: Time step successfully carried out
*            = 1: Nonlinear solver broke down
*            = 2: Linear solver broke down
*   NNL    : The number of iterations of the nonlinear solver
*   NMG    : The number of iterations of the linear solver.
*   DUP    : receives the new approximation to the solution vector
*
*   IPARAM.CURTIM:
*            The current simulation time is updated according to the
*            time step setting in TIMSTP.
*   DPARAM.TMTOT,
*   DPARAM.TNSTIM,
*   DPARAM.TTGRID:
*            Will be modified according to the time the time step needs.
*            The necessary time is added to the structure/variable.
*
* For time-dependent simulation (IPARAM.INONST=1):
*   DUP    : receives the new iterate u^(n+1)
*   DRHS   : Right hand side f^n+1 in the next time step without
*            boundary conditions being implemented.
*
* Depending on the configuration of the problem, the following
* structures are changed:
* 1.) IBND >= 2 - time-dependent boundary conditions or moving 
*                 fictitious boundary components
*     => The boundary condition structures in CURTRI are updated
*        according to the time of the next time step.
*        TRIAS will serve as prototype for the creation of CURTRI.
* 2.) IMESH = 1,2 - moving mesh in every time step
*     => The boundary condition structures in CURTRI are updated
*        according to the time corresponding to u^n+1.
*        TRIAS will serve as prototype for the creation of CURTRI.
*     => The mesh coordinates CURTRI.DCORVG/DCORMG are changed
*        according to the time corresponding to u^n+1.
*        TRIAS will serve as prototype for the creation of CURTRI.
*        Whether or not the coordinates TRIAS.DCORVG/DCORMG serve
*        as initial state for the grid deformation depends on
*        whether CURTRI was created from TRIAS by duplicating
*        the arrays DCORVG/DCORMG or simply by copying the handles
*        (Remember that CURTRI is assumed to be created from
*         TRIAS by TRIDUP).
*
*        Remark that if the DCORVG/DCORMG arrays are shared
*        between TRIAS and CURTRI, the original vertex coordinates
*        in TRIAS are changed during the grid deformation!
* 3.) IALE > 0 - ALE method
*     => The grid velocity field in the discretisation structure
*        is updated.
************************************************************************

      IMPLICIT NONE
      
C main COMMON blocks
      
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
      
C     parameters
      
      INTEGER IPARAM (SZISDI),NNL,NMG,ISTATN,NUVP
      DOUBLE PRECISION DPARAM (SZISDD),DMGPAR(*),DSTPAR(*),DASMBL(*)
      DOUBLE PRECISION DUP(*),DGEOM(*)
      DOUBLE PRECISION DRHS(*),DAUX(*)
      INTEGER TRIAS(SZTRIA,NNLEV),CURTRI(SZTRIA,NNLEV)
      INTEGER MATDAT(SZN2MI,NNLEV),VECDAT(SZN2VI,NNLEV),NLMIN,NLMAX
      INTEGER IMGPAR(*),ISTPAR(*),IASMBL(*),IGEOM(*)
      DOUBLE PRECISION TIMSTP(SZTSTD)

C     generation of RHS and grid

      EXTERNAL GENRHS, GENGRI, PRCNTS

C     externals

C     Parametrisation of the domain

      DOUBLE PRECISION PARX,PARY,TMAX
      EXTERNAL PARX,PARY,TMAX
      
C     Coefficient of stiffness matrix, right hand side, exact solution

      DOUBLE PRECISION UE,UEX,UEY
      EXTERNAL UE,UEX,UEY
      
C     definition of finite elements

      EXTERNAL E030,E031,EM31,EM30

C     local variables

      INTEGER STATUS,I,IREP,MSREP,NEQU,NEQP
C     INTEGER KRHS

      INTEGER ILEV,KGRVEL
      DOUBLE PRECISION TOSTEP
      INTEGER KM1,KST1,KA1,KCOLA,KLDA,NA,KVERT,KMID,KCORVG
      LOGICAL BGRADAP

C     Backup of the discretisation structures
      
      INTEGER IASMBK (SZASMI)
      DOUBLE PRECISION DASMBK(SZASMD)

C     Backup of the configuration of the stationary solver
      
      INTEGER ISTBCK(SZNSDI)
      DOUBLE PRECISION DSTBCK(SZNSDD)
      
C     Possible backup of solution vector, RHS, simulation time

      INTEGER LUBK,LRHSBK,J,SOLTRI(SZTRIA,NNLEV)
      DOUBLE PRECISION TIMEBK
      
C     Allocate a solver structure for statistical timing information:

      DOUBLE PRECISION DTIMIN(SZISDD)

C=======================================================================

C     Stop the time we need for the time step:

      CALL GTMAUX (DTIMIN,DPARAM,OTMTOT,0)
      
C     Get the vector size of velocity and pressure part on the finest
C     level; they are frequently used.

      NEQU = VECDAT(ONU,NLMAX)
      NEQP = VECDAT(ONP,NLMAX)

C     If the ALE method is used, store the address of the grid velocity
C     vector to KGRVEL. Otherwise store 1 into KGRVEL, which is used
C     as a dummy in the call to the upwinding routines to prevent
C     specifying a DWORK(.) address out of bounds.

      IF (IASMBL(OIALE).NE.0) THEN
        KGRVEL = L(IASMBL(OLGRVEL))
      ELSE
        KGRVEL = 1
      END IF

C     BGRADAP receives a flag whether or not grid adaption takes place
C     in this time step.

      BGRADAP = ((IPARAM(OIDMESH).EQ.2)   .OR.
     *           (IPARAM(OIDMESH).EQ.3))                     .AND.
     *          ((IPARAM(OIMSSTP).EQ.0)  .OR.
     *           (MOD(IPARAM(OCRITE),IPARAM(OIMSSTP)).EQ.0)).AND.
     *          ((DPARAM(ODMSMTM).LT.0D0).OR.
     *           (DASMBL(OTIMENS).LE.DPARAM(ODMSMTM)))

C     Probably the user wants to adapt the grid every IMSSTP timesteps.
C     If we are in such a step where we must adapt the grid, we have to
C     make a backup of the current grid. This backup must serve as basis
C     for the FE function in the evaluation of the monitor function, while
C     the real computational grid is deformed (and therefore invalid at
C     that moment).
      
      IF ((IPARAM(OIMSSTP).GE.1).AND.BGRADAP) THEN
 
C       In the duplication of TRIAS, we basically backup everything:
        
        J = 0
        
C       For time-dependend boundary conditions (which include time-
C       dependent fictitious boundary objects), make backups of KNPR/KXNPR
C       additionally to the vertex coordinates:

        IF (IPARAM(OIBDR).NE.0) J = NOT(2**0 + 2**15)

C       For IMESH=1-3 don't duplicate structural information (like KVERT), 
C       as this is not changed in the grid deformation process. Other
C       information like coordinates,... we have to back up so we can
C       restore it when necessary.

        IF ((IPARAM(OIDMESH).GE.1).AND.(IPARAM(OIDMESH).LE.3)) 
     *    J = NOT(2**0 + 2**1 + 2**12 + 2**13 + 2**15)
     
      ELSE
      
C       In the duplication of TRIAS, we can share all information 
C       with TRIAS, since nothing will be changed.

        J = NOT(0)
        
      END IF

C     If the user wants to repeat every time step with a grid that
C     is adapted depending on an error criterion, we have to be able
C     to return to the beginning of this time step. In that case,
C     make backups of the current simulation time, the RHS vector
C     and the solution vector.

      IF ((IPARAM(OIMSREP).GE.1).AND.
     *    ((IPARAM(OIMSSTP).GE.1).AND.BGRADAP)) THEN
      
        CALL ZNEW(NUVP,-1,LUBK  ,'DUBK  ')
        IF (IER.NE.0) RETURN
        CALL LCP1(DUP,DWORK(L(LUBK)),NUVP)
        
C       The RHS has only to be duplicated if it's time dependent...

        LRHSBK = 0
        IF (IPARAM(OIRHS).GE.2) THEN
          CALL ZNEW(NUVP,-1,LRHSBK,'DRHSBK')
          IF (IER.NE.0) RETURN
          CALL LCP1(DRHS,DWORK(L(LRHSBK)),NUVP)
        END IF
        
        TIMEBK = DASMBL(OTIMENS)
 
      ELSE

C       Set the handles to 0 to indicate that we don't have a backup...

        LUBK   = 0
        LRHSBK = 0
        TIMEBK = 0D0
        
      END IF
      
C     Start the grid duplication for all levels.
C     Duplicate (based on the variable J which was set above) the TRIAS
C     structures to SOLTRI. When mesh-readaption in a time step
C     is activated, SOLTRI then defines the triangulations where the
C     current solution vector / RHS is defines -- which allowes error
C     calculation and mesh adaption based on that.
      
      CALL GTMAUX (DTIMIN,DPARAM,OTTGRID,0)
      
      DO I=NLMAX,NLMIN,-1

C       Copy DCORVG only on the maximum level, since the information
C       is shared between all levels!

        IF (I.EQ.NLMAX) THEN
          CALL TRIDUP(CURTRI(1,I),SOLTRI(1,I),0,J)
        ELSE
          CALL TRIDUP(CURTRI(1,I),SOLTRI(1,I),0,IOR(J,2**0))
          SOLTRI(OLCORVG,I) = SOLTRI(OLCORVG,NLMAX)
        END IF
        
      END DO
      
      CALL GTMAUX (DTIMIN,DPARAM,OTTGRID,1)
      
C     To calculate the solution for the new time step, we make
C     IMSREP retries - except for if there is no grid adaption in time...

      MSREP = MAX(0,IPARAM(OIMSREP))
      IF (IPARAM(OIDMESH).LE.1) MSREP = 0
      IF ((IPARAM(OIMSSTP).EQ.0).OR.
     *    (MOD(IPARAM(OCRITE),IPARAM(OIMSSTP)).NE.0)) MSREP = 0

      DO IREP = 0,MSREP
      
C       In nonstationary simulations we have to prepare
C       a lot of stuff: New system matrices have to be generated,
C       new RHS have to be generated,...

C***********************************************************************
C       The Nav.St.-Equation (cf. p. 38 (43), Turek's book)
C
C           u_t + N(u)u + grad(p) = f,   div(u)=0
C
C       is discretised in time
C
C           [ I + Theta k N(u)]u + k grad(p) = g,   div(u)=0
C
C       and space
C
C           [ M + Theta k N(u_h)]u_h + k B p_h = g_h,   B^T u_h = 0
C
C       resulting in the abstract system
C
C           S(u_h)u_h  +  k B p_h  =  g_h,   B^T u_h = 0
C
C       with
C
C           S(u) = [alpha*M + Theta_1*nu*k*L + Theta_2*k*K(u)] u
C
C       Here, alpha=0 for stationary and alpha=1 for
C       instationary problems, L=KST1=Laplace-Matrix.
C       The time-step constant k=TSTEP in front of the pressure 
C       B p_h in the abstract system above can be realised by 
C       scaling the pressure in every time-step.
C***********************************************************************

C***********************************************************************
C       What do we have to prepare before calling the nonlinear
C       stationary solver? Well, basically we have to do only
C       one thing: Prepare the Right-Hand-Side for the new time step!
C
C       But this decomposes in a couple of subtasks, which all happen
C       at different points in the simulational time. We have to
C       - adapt the mesh if necessary
C       - include boundary conditions
C       - plug together the RHS of the previous time step with a
C         new one, if the RHS is time dependent.
C       We structure the assembling here according to the simulation
C       time. At first we do everything that must be done with the data
C       of the old time step, then we increase our time step and
C       do the rest.
C***********************************************************************
C       We give a rough overview about what we have to assemble and
C       which notation we use. The assembling of the right-hand-side
C       vector for the nonlinear iteration is based on the vector f^n
C       (currently in DRHS), f^n+1 (which has to be created)
C       and the nonlinear term depending on the current solution 
C       (see page 165/166 (151/152) in Turek's CFD-book):  
C
C       [I - theta_1 k N(u^n)]u^n + Theta_2 k f^(n+1) + Theta_3 k f^n
C
C       with k=Delta(t) being the step size of the substep. 
C
C       We assemble the RHS in the DAUX array.
C
C       Classical 1-step scheme
C       -----------------------
C       In the classical 1-step scheme, the RHS has the form 
C       (cf. p.152):  
C
C       [I - (1-Theta) k N(u^n)]u^n + Theta k f^(n+1) + (1-Theta) k f^n
C          ^^^^^^^^^^^^^              ^^^^^^^           ^^^^^^^^^^^
C            TRMWGH                   TR1WGH             TR2WGH
C
C       In case we have steady (in-)homogenuous RHS, we have
C       f^n = f^(n+1) = DRHS and therefore the latter both terms
C       reduce to:
C
C         Theta k f^(n+1) + (1-Theta) k f^n
C       = Theta k f^(n+1) + (1-Theta) k f^(n+1)
C       = k f^(n+1)
C       =    k      DRHS
C         ^^^^^^^
C         TRSWGH
C
C       Fractional step Theta-scheme
C       ----------------------------
C       This is a little bit more complicated due to the different
C       substeps. As RHS we (again) sum together:
C
C       [I + TRMWGH N(u^n)]u^n + TR1WGH f^(n+1) + TR2WGH f^n
C
C       But this time the step-variables have different values
C       depending on the step. With TR2WGH=-TRMWGH we have:
C
C           Step:       1                2                  3
C       Variable:
C          TRMWGH   -beta*Theta*K  alpha*Theta'*K      -beta*Theta*K
C          TR1WGH   alpha*Theta*K  -beta*Theta'*K      alpha*Theta*K
C
C       or in formulas, respectively:
C
C           Step:       1                2                  3
C       Variable:
C          TRMWGH   -beta*Theta*K  alpha*(1-2Theta)*K   -beta*Theta*K
C          TR1WGH   alpha*Theta*K  -beta*(1-2Theta)*K   alpha*Theta*K
C
C       with:   alpha = FALPHA = (1-2Theta)/(1-Theta)
C               beta  = FBETA  = Theta/(1-Theta)
C               Theta = 1-2/sqrt(2)
C               Theta'= (1-2Theta)
C
C         WARNING: The documentation on page 165 (151) in Turek's book 
C          on the Theta Scheme is different than implemented here!!!
C          This version of the Theta Scheme has never been published,
C          but was derived independently of Glowinski by Rannacher/
C          Turek! It makes more sense than the Glowinski-scheme
C          as the coefficients before f^n+1 and f^n add up to 1 and
C          the RHS f^(n+1-Theta) is used neither in the second nor
C          the third step!
C***********************************************************************

C       The RHS handling is easiest if the RHS is stationary.
C
C       If we have homogeneous RHS, simply force the RHS to 0.
        
        IF (IPARAM(OIRHS).EQ.0) THEN
        
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,0)
          CALL LCL1(DAUX,NUVP)
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,1)
          
        ENDIF

C       If we have steady inhomogenuous RHS, copy the calculated
C       RHS-vector from DRHS to DF1, weighted by the coefficients
C       of the Theta-scheme.    
C       The RHS is normally build in the INIT1 with the help
C       of the INDATxD.F-file.

        IF (IPARAM(OIRHS).EQ.1) THEN

C         The RHS reduces as stated above. We only have to multiply by
C         
C                    TRSWGH  =  TR1WGH + TR2WGH  =  k

          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,0)
          CALL LLC1(DRHS,DAUX,NUVP,TIMSTP(OTRSWGH),0D0)
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,1)
          
        ENDIF

C       Implement pressure drop conditions of the current simulation
C       time into the RHS if the RHS is stationary. Note that the
C       pressure drop is generally time dependent, so we can not
C       simply add the pressure drop conditions once with weight TRSWGH.

        IF ( ((IPARAM(OIRHS).EQ.0).OR.(IPARAM(OIRHS).EQ.1)) .AND.
     *       (IAND(IPARAM(OIBDR),1).NE.0) ) THEN
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTBDR,0)
          CALL PDSETX (DAUX(1),DAUX(1+NEQU),
     *                 DASMBL(OTIMENS),TIMSTP(OTR2WGH),DASMBL(ORE),
     *                 CURTRI(1,NLMAX),IASMBL,DASMBL,
     *                 IGEOM,DGEOM)
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTBDR,1)
        END IF

C       If we have nonsteady inhomogenuous RHS, we have to calculate 
C       it - at different points in time. Start the assembling with
C       f^n in the current time step:

        IF (IPARAM(OIRHS).EQ.2) THEN
        
C         We use the DAUX vector for assembling the RHS. 
C
C         We go "from back to front" when adding the terms
C         to the RHS in the Theta scheme (cf. p. 164 (150), Turek's 
C         book). Take the current RHS f^n, weight it according 
C         to the Theta-scheme with (1-Theta)*k and write it into DAUX
C         as start point of the assembling.

          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,0)
          CALL LLC1(DRHS,DAUX,NUVP,TIMSTP(OTR2WGH),0D0)
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,0)

C         Implement pressure-drop boundary conditions into
C         RHS-vectors DF=(DF1,DF2,DP).
C         Use the current simulational time, which corresponds to f^n.
          
          IF (IPARAM(OIBDR).GE.1) THEN
            CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTBDR,0)
            CALL PDSETX (DAUX(1),DAUX(1+NEQU),
     *                   TIMSTP(OTIMENS),TIMSTP(OTR2WGH),DASMBL(ORE),
     *                   CURTRI(1,NLMAX),IASMBL,DASMBL,
     *                   IGEOM,DGEOM)
            CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTBDR,1)
          ENDIF ! (IPARAM(OIBDR).GE.1)
          
        END IF ! (IPARAM(OIRHS).EQ.2)

C-----------------------------------------------------------------------

C       Now that we completed everything of the current time step,
C       we proceed to the next one. This involves...
C        - increasing the simulation time
C        - handle the pressure
C        - adapt the mesh if necessary
C        - calculate new boundary conditions
C        - Assemble the rest of the right hand side
C
C       At first tell the caller that we proceed to the next time step...

        CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,0)
        CALL PRCNTS (NLMIN,NLMAX,
     *               CURTRI,MATDAT,VECDAT,
     *               IPARAM,DPARAM,ISTPAR,DSTPAR,
     *               IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *               NUVP,DUP,DRHS,DAUX,
     *               DASMBL(OTIMEST), DASMBL(OTIMENS), TIMSTP(OTSTEP))
        CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTPOST,1)

C       ...and then increase the simulation time in the assembly
C       structure.

        DASMBL(OTIMENS) = DASMBL(OTIMENS) + TIMSTP(OTSTEP)

C-----------------------------------------------------------------------

C       Scale the pressure by the step-length.
C       That's all for the handling of the pressure - it's independent
C       of the current simulation time.

        CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,0)
        CALL LLC1(DUP(1+2*NEQU),DUP(1+2*NEQU),NEQP,0D0,TIMSTP(OTSTEP))
        CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,1)

C-----------------------------------------------------------------------
C       Mesh adaption part with ALE method.
C
C       Check if the ALE-method is active. In this case, make a backup of
C       the corner coordinates of the grid on the finest level into the
C       grid velocity vector in the IASMBL structure.
C       This will later be used to calculate the grid velocity.

        IF (IASMBL(OIALE).NE.0) THEN
          CALL GTMAUX (DTIMIN,DPARAM,OTTGRID,0)
          CALL LCP1(DWORK(L(CURTRI(OLCORVG,NLMAX))),
     *              DWORK(KGRVEL),
     *              CURTRI(ONVT,NLMAX))
          CALL GTMAUX (DTIMIN,DPARAM,OTTGRID,1)
        END IF
          
C       Now: do we have to adapt the mesh for the new time?

        IF (IAND(ISTATN,1).NE.0) THEN

          IF (((IPARAM(OIDMESH).EQ.2).OR.(IPARAM(OIDMESH).EQ.3)).AND.
     *        ((IPARAM(OIMSSTP).EQ.0).OR.
     *         (MOD(IPARAM(OCRITE),IPARAM(OIMSSTP)).EQ.0)) ) THEN
     
C           Start the grid deformation.
C           
C           Should we reset our grid to the original?

            CALL GTMAUX (DTIMIN,DPARAM,OTTGRID,0)

            IF (IPARAM(OIDMESH).EQ.2) THEN
            
C             Go through all levels and call the grid-restore routine
C             to restore the original mesh.
C             Restore everything we have duplicated.

              DO I=NLMIN,NLMAX
                CALL TRIRST (CURTRI(1,I),TRIAS(1,I))
              END DO
            
            END IF
            
C           ELSE: stay at the current grid and modify that one.
C
C           Generate the new grid into CURTRI.
C           In the first adaption process, never use adaption by
C           error control -- only in step 2,3,4,...

            IF ((IREP.EQ.0).AND.(MSREP.GT.0)) THEN
              CALL GENGRI (NLMIN,NLMAX,CURTRI,SOLTRI,0,
     *                   IPARAM,DPARAM,ISTPAR,DSTPAR,
     *                   IMGPAR,DMGPAR,IASMBL,DASMBL,
     *                   IGEOM,DGEOM,
     *                   VECDAT(1,NLMAX),NUVP,DUP,DRHS,DAUX,
     *                   DPARAM(OTIMEST),DASMBL(OTIMENS), 
     *                   TIMSTP(OTSTEP))
            ELSE
              CALL GENGRI (NLMIN,NLMAX,CURTRI,SOLTRI,IPARAM(OIMSCRT),
     *                   IPARAM,DPARAM,ISTPAR,DSTPAR,
     *                   IMGPAR,DMGPAR,IASMBL,DASMBL,
     *                   IGEOM,DGEOM,
     *                   VECDAT(1,NLMAX),NUVP,DUP,DRHS,DAUX,
     *                   DPARAM(OTIMEST),DASMBL(OTIMENS), 
     *                   TIMSTP(OTSTEP))
            END IF
     
            CALL GTMAUX (DTIMIN,DPARAM,OTTGRID,1)
     
C           A modified grid leads directly to modified Laplace/Pressure/
C           Mass matrices! We have to rebuild everything, based on CURTRI!

            CALL GTMAUX (DTIMIN,DPARAM,OTTLC,1)
            CALL GENNSM (NLMIN,NLMAX,CURTRI,IASMBL,DASMBL,7,MATDAT)
            CALL GTMAUX (DTIMIN,DPARAM,OTTLC,1)
     
          END IF

        END IF ! AND(ISTATN,1) <> 0

C       In case the ALE method is used, calculate the grid velocity
C       based on the previous coordinates and the new ones.
C       The LGRVEL handle points to the old coordinates, CURTRI to
C       the new ones, so (new coordinates - old coordinates)/steplength
C       is a 1st order approximation to the grid velocity.

        IF (IASMBL(OIALE).NE.0) THEN
          CALL LLC1(DWORK(L(CURTRI(OLCORVG,NLMAX))),
     *              DWORK(KGRVEL),CURTRI(ONVT,NLMAX),
     *              1D0/TIMSTP(OTSTEP),-1D0/TIMSTP(OTSTEP))
        END IF

C       Check if our boundary conditions or the mesh moves...

        IF ((IAND(IPARAM(OIBDR),2).NE.0).OR.(IPARAM(OIDMESH).GE.2)) THEN
        
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTBDR,0)
        
C         Either our mesh has changed due to grid deformation, or
C         we have time-dependent Dirichlet/Neumann boundary 
C         conditions! (IBDR>=2 from the DAT-file).
C         Implement the new boundary conditions into the mesh
C         information.
        
          DO ILEV=NLMIN,NLMAX
          
C           Duplicate the prototype of the old grid, reconstruct
C           time-independent boundary conditions:

            CALL TRIDUP (TRIAS(1,ILEV),CURTRI(1,ILEV),1,
     *                   NOT(2**7+2**15))
          
          END DO ! ILEV
     
C         Calculate boundary conditions, store them into CURTRI.
C
C         This will modify CURTRI.KXNPR according to Dirichlet/
C         Neumann boundary information. As KXNPR we specify KXNPR of
C         CURTRI, overwriting the old CURTRI.KXNPR where necessary.
C
C         Calculate the number of boundary vertices NBDMT and the
C         shortcut nodal property array - both are stored in the
C         user defined block at CURTRI.TRIUD(1) / CURTRI.TRIUD(2).
C
C         Furthermore this routine updates INEUM in the assembly
C         structure whether there are Neumann components at all in 
C         our problem.
C
C         As simulation time we use the time of the next time step -
C         we just switched to that.

          CALL GNBDGE (NLMIN,NLMAX,CURTRI,
     *                 IASMBL,DASMBL,IGEOM,DGEOM)
     
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTBDR,1)
          
        ENDIF ! (IBDR.GE.2)
     
C***********************************************************************

C       Implement pressure drop conditions of the current simulation
C       time into the RHS if the RHS is stationary. Note that the
C       pressure drop is generally time dependent, so we cannot
C       simply add the pressure drop conditions once with weight TRSWGH.

        IF ( ((IPARAM(OIRHS).EQ.0).OR.(IPARAM(OIRHS).EQ.1)) .AND.
     *       (IAND(IPARAM(OIBDR),1).NE.0) ) THEN
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTBDR,0)
          CALL PDSETX (DAUX(1),DAUX(1+NEQU),
     *                 DASMBL(OTIMENS),TIMSTP(OTR1WGH),DASMBL(ORE),
     *                 CURTRI(1,NLMAX),IASMBL,DASMBL,
     *                 IGEOM,DGEOM)
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTBDR,1)
        END IF

C       When the RHS is time-dependent, we have to continue 
C       assembling it. The previous assembling only handled everything
C       related to f^n - now we must handle f^n+1:

        IF (IPARAM(OIRHS).EQ.2) THEN
        
C         Assemble the RHS-vector f^(n+1) in RHS overwriting the old one.
          
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,0)
          
          CALL GENRHS (CURTRI(1,NLMAX),
     *                 IPARAM,DPARAM,IASMBL,DASMBL,IGEOM,DGEOM,
     *                 VECDAT(1,NLMAX),DRHS)
     
C         Add the newly assembled RHS-vectors to the RHS-vectors from
C         above, weighted according to the Theta-scheme with Theta*k.

          CALL LLC1 (DRHS,DAUX,2*NEQU,TIMSTP(OTR2WGH),1D0)

C         We don't implement any boundary conditions into DRHS so DRHS is
C         the proper RHS vector for the next time step.

          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,1)

C         Implement pressure drop conditions for the new simulation 
C         time into the RHS vectors in DAUX. This will simply add the
C         missing information for the new time step.

          IF (IPARAM(OIBDR).GE.1) THEN
            CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTBDR,0)
            CALL PDSETX (DAUX(1),DAUX(1+NEQU),
     *                   DASMBL(OTIMENS),TIMSTP(OTR1WGH),DASMBL(ORE),
     *                   CURTRI(1,NLMAX),IASMBL,DASMBL,
     *                   IGEOM,DGEOM )
            CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTBDR,1)
          ENDIF ! (IBDR.GE.1)

C         Implement boundary conditions into that RHS vector.
C         We'll do that later again when the complete RHS is assembled;
C         but it has to be done here, too, to prevent the later
C         XADF6 to calculate wrong values.

          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTBDR,0)
    
          CALL BDRSTX (DAUX(1),DAUX(1+NEQU),
     *               KWORK(L(CURTRI(OLXNPR,NLMAX))),
     *               PARX,PARY,UE,DASMBL(OTIMENS),
     *               DASMBL(ORE),
     *               CURTRI(1,NLMAX),IASMBL,DASMBL,IGEOM,DGEOM) 
    
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTBDR,1)

        ENDIF ! (IRHS.EQ.2)

C======================================================================
C
C       The f^n and f^(n+1)-parts are now assembled together into
C       DAUX. What's still missing is the nonlinear block with the
C       data of the "old" solution DU=u^n. This we have to completely
C       reassemble...
C                     
C       Fetch some information about the matrices from the matrix block

        KM1    = L(MATDAT(OLM,NLMAX))
        KST1   = L(MATDAT(OLST,NLMAX))
        KA1    = L(MATDAT(OLA1,NLMAX))
        KCOLA  = L(MATDAT(OLCLA1,NLMAX))
        KLDA   = L(MATDAT(OLLDA1,NLMAX))
        NA     = MATDAT(ONA1,NLMAX)
        
        KVERT  = L(CURTRI(OLVERT,NLMAX))
        KMID   = L(CURTRI(OLMID,NLMAX))
        KCORVG = L(CURTRI(OLCORVG,NLMAX))
                     
C       Set up linear part of RHS-vector of Theta-scheme as well as
C       linear part of nonlinear system matrix. Use the current
C       vector in KFx as source and overwrite it.
C       -->   KA1 = [  M*I  + TRMWGH * (-nu * Laplace(.))  ]
C       -->   DF  = [ M*u^n + TRMWGH*(-nu*Laplace(u^n)) ]  +  KFx
C                                                            ^^^^^
C                  Remember:  KFx = TR1WGH*f^(n+1) + TR2WGH*RHS
C                                 = TRSWGH*RHS  (in steady case)
C                             TR2WGH = (1-Theta)*k

        CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTADF,0)
        CALL XMADF6(CURTRI(1,NLMAX),MATDAT(1,NLMAX),DAUX,DUP,
     *              TIMSTP(OTRMWGH),IASMBL(OIPRECA),NLMAX)
        CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTADF,1)
        
C       Now treat the convective part.
C
C       Modify the current RHS-vector KF such that
C           OTRMSTP * u grad(u)
C       is subtracted - but because we want to *add* this term,
C
C           KFx = KFx + TRMWGH * u grad(u)
C                     ^^^^^^^^^^^^^^^^^^^^^^
C               = [ M*I + TRMWGH*N(u) ]u + TR1WGH*f^(n+1) + TR2WGH*f^n
C
C       we have to switch the sign of THSTEP!
C       The matrix is not to be modified (we don't need it), so we
C       call the assembling routine with IDEF=2.
        
        IF (TIMSTP(OTRMWGH).NE.0D0) THEN

          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTUPW,0)
          
          IF (IASMBL(OIUPW).EQ.1) THEN
C            CALL GUPWD (DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
C       *           1D0,0D0,DUP(1),DUP(1+NEQU),
C       *           DAUX(1),DAUX(1+NEQU),DWORK(KA1),KWORK(KCOLA),
C       *           KWORK(KLDA),KWORK(KVERT),KWORK(KMID),
C       *           DWORK(KCORVG),CURTRI(ONEL,NLMAX),CURTRI(ONVT,NLMAX),2,
C       *           DASMBL(OUPSAM),DASMBL(ORE),-TIMSTP(OTRMWGH))
            CALL GUPWDX (DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
     *             1D0,0D0,DUP(1),DUP(1+NEQU),
     *             DAUX(1),DAUX(1+NEQU),DWORK(KA1),KWORK(KCOLA),
     *             KWORK(KLDA),CURTRI(1,NLMAX),
     *             KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),2,
     *             DASMBL(OUPSAM),DASMBL(ORE),-TIMSTP(OTRMWGH),
     *             IASMBL(OIALE),DWORK(KGRVEL))
          ELSE
          
            IF (IASMBL(OIELEMT).EQ.0) THEN
C              CALL SUPGPX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
C       *                1D0,0D0,DUP(1),DUP(1+NEQU),
C       *                DAUX(1),DAUX(1+NEQU),
C       *                DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
C       *                KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
C       *                CURTRI(1,NLMAX),E031,IASMBL(OICUBN),
C       *                (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
C       *                IASMBL(OISTOK), DASMBL(OUPSAM), DASMBL(ORE),
C       *                2,-1D0,-TIMSTP(OTRMWGH))
              CALL SUPAPX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
     *                  DAUX(1),DAUX(1+NEQU),
     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
     *                  CURTRI(1,NLMAX),E031,IASMBL(OICUBN),
     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
     *                  IASMBL(OISTOK), DASMBL(OUPSAM), DASMBL(ORE),
     *                  2,-1D0,-TIMSTP(OTRMWGH),
     *                  IASMBL(OIALE),DWORK(KGRVEL))
            END IF

            IF (IASMBL(OIELEMT).EQ.1) THEN
C              CALL SUPGPX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
C       *                1D0,0D0,DUP(1),DUP(1+NEQU),
C       *                DAUX(1),DAUX(1+NEQU),
C       *                DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
C       *                KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
C       *                CURTRI(1,NLMAX),E030,IASMBL(OICUBN),
C       *                (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
C       *                IASMBL(OISTOK), DASMBL(OUPSAM), DASMBL(ORE),
C       *                2,-1D0,-TIMSTP(OTRMWGH))
              CALL SUPAPX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
     *                  DAUX(1),DAUX(1+NEQU),
     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
     *                  CURTRI(1,NLMAX),E030,IASMBL(OICUBN),
     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
     *                  IASMBL(OISTOK), DASMBL(OUPSAM), DASMBL(ORE),
     *                  2,-1D0,-TIMSTP(OTRMWGH),
     *                  IASMBL(OIALE),DWORK(KGRVEL))
            END IF

            IF (IASMBL(OIELEMT).EQ.2) THEN
C              CALL SUPGNX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
C       *                1D0,0D0,DUP(1),DUP(1+NEQU),
C       *                DAUX(1),DAUX(1+NEQU),
C       *                DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
C       *                KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
C       *                CURTRI(1,NLMAX),EM31,IASMBL(OICUBN),
C       *                (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
C       *                IASMBL(OISTOK), DASMBL(OUPSAM), DASMBL(ORE),
C       *                2,-1D0,-TIMSTP(OTRMWGH))
              CALL SUPANX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
     *                  DAUX(1),DAUX(1+NEQU),
     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
     *                  CURTRI(1,NLMAX),EM31,IASMBL(OICUBN),
     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
     *                  IASMBL(OISTOK), DASMBL(OUPSAM), DASMBL(ORE),
     *                  2,-1D0,-TIMSTP(OTRMWGH),
     *                  IASMBL(OIALE),DWORK(KGRVEL))
            END IF

            IF (IASMBL(OIELEMT).EQ.3) THEN
C              CALL SUPGNX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
C       *                1D0,0D0,DUP(1),DUP(1+NEQU),
C       *                DAUX(1),DAUX(1+NEQU),
C       *                DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
C       *                KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
C       *                CURTRI(1,NLMAX),EM30,IASMBL(OICUBN),
C       *                (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
C       *                IASMBL(OISTOK), DASMBL(OUPSAM), DASMBL(ORE),
C       *                2,-1D0,-TIMSTP(OTRMWGH))
              CALL SUPANX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
     *                  DAUX(1),DAUX(1+NEQU),
     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
     *                  CURTRI(1,NLMAX),EM30,IASMBL(OICUBN),
     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
     *                  IASMBL(OISTOK), DASMBL(OUPSAM), DASMBL(ORE),
     *                  2,-1D0,-TIMSTP(OTRMWGH),
     *                  IASMBL(OIALE),DWORK(KGRVEL))
            END IF

          ENDIF ! (IUPW.EQ.1)
        
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTUPW,1)
          
        ELSE

C         THSTEP=0D0 - This is the case for backward Euler only.
C         (THSTEP = c(1-THETA) = 0  <->  THETA=1)

          IF ((IASMBL(OIPRECA).EQ.4).AND.(IASMBL(OIMASS).EQ.1)) THEN

C           The handling is slightly different because we have
C           to set THSTEP=1 during the upwinding in this case,
C           and we must use DCMASS=0.

            CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTUPW,0)

C           Use a step length of 1D0 in call to the creation of the
C           convective part:
            
            TOSTEP=1D0
            
            IF (IASMBL(OIELEMT).EQ.0) 
     *        CALL SUPGPX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
     *                  DAUX(1),DAUX(1+NEQU),
     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
     *                  CURTRI(1,NLMAX),E031,IASMBL(OICUBN),
     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
     *                  IASMBL(OISTOK), DASMBL(OUPSAM), DASMBL(ORE),
     *                  2,0D0,TOSTEP)

            IF (IASMBL(OIELEMT).EQ.1) 
     *        CALL SUPGPX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
     *                  DAUX(1),DAUX(1+NEQU),
     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
     *                  CURTRI(1,NLMAX),E030,IASMBL(OICUBN),
     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
     *                  IASMBL(OISTOK), DASMBL(OUPSAM), DASMBL(ORE),
     *                  2,0D0,TOSTEP)

            IF (IASMBL(OIELEMT).EQ.2) 
     *        CALL SUPGNX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
     *                  DAUX(1),DAUX(1+NEQU),
     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
     *                  CURTRI(1,NLMAX),E031,IASMBL(OICUBN),
     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
     *                  IASMBL(OISTOK), DASMBL(OUPSAM), DASMBL(ORE),
     *                  2,0D0,TOSTEP)

            IF (IASMBL(OIELEMT).EQ.3) 
     *        CALL SUPGNX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
     *                  DAUX(1),DAUX(1+NEQU),
     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
     *                  CURTRI(1,NLMAX),E030,IASMBL(OICUBN),
     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
     *                  IASMBL(OISTOK), DASMBL(OUPSAM), DASMBL(ORE),
     *                  2,0D0,TOSTEP)

          ENDIF ! ((IPRECA.EQ.4).AND.(IMASS.EQ.1))
          
          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTUPW,1)
          
        ENDIF ! (THSTEP.NE.0D0)

C       Incorporate Dirichlet boundary conditions into current
C       solution and RHS vector. This finishes the assembling of
C       the right hand side...

        CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTBDR,0)
        
        CALL BDRSTX (DUP(1),DUP(1+NEQU),
     *               KWORK(L(CURTRI(OLXNPR,NLMAX))),
     *               PARX,PARY,UE,DASMBL(OTIMENS),
     *               DASMBL(ORE),
     *               CURTRI(1,NLMAX),IASMBL,DASMBL,IGEOM,DGEOM)   
            
        CALL BDRSTX (DAUX(1),DAUX(1+NEQU),
     *               KWORK(L(CURTRI(OLXNPR,NLMAX))),
     *               PARX,PARY,UE,DASMBL(OTIMENS),
     *               DASMBL(ORE),
     *               CURTRI(1,NLMAX),IASMBL,DASMBL,IGEOM,DGEOM) 
              
        CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTBDR,1)

************************************************************************
C *** fixed point defect correction for stationary NS equations
C***********************************************************************

C       Prepare the nonlinear solver.
C
C       Make a backup of the stationary solver structures so we
C       can restore it on exit. For solving stationary problems,
C       we modify the structures temporarily.

        CALL LCP3(ISTPAR,ISTBCK,SZNSDI)
        CALL LCP1(DSTPAR,DSTBCK,SZNSDD)

C       Make a backup of the assembly structures, so that we can restore 
C       the old state when we are finished.

        CALL LCP3 (IASMBL,IASMBK,SZASMI)
        CALL LCP1 (DASMBL,DASMBK,SZASMD)

C       In a nonstationary simulation, the mass matrix has to be added
C       to the system matrix following the discretisation rule
C                  [ M  +  THWEIG * N(u) ] u + ... = ...
C       This is because of the time discretisation du/dt. So we have to
C       impose IALPHA=1 in the discretisation structure to enforce
C       NSDEF not to forget this mass matrix.

        IASMBL(OIALPHA) = 1
        
C       Set the THSTEP-parameter in the parameter block of the
C       nonlinear solver to TMWGH from our time stepping scheme.
C       This is the weight for the N(.) matrix on the left hand side
C       of the equation! As example to identify it, we consider the
C       equation in the first step of the Fractional-Step Theta-
C       scheme (cf. p. 167 (153), Turek's book):
C
C          [ I + alpha k Theta A_(n+theta) ] u_(n+Theta)  =  RHS
C                ^^^^^^^^^^^^^ ^^^^^^^^^^^
C                   TMSTEP     Created in NSDEF  

        DASMBL(OTHWEIG) = TIMSTP(OTMWGH)

C       Call NSDEF2 and perform fixed point defect correction
C       to calculate the next iterate of this substep.
C       This will update the solution vector in DU1/DU2/DP.
C       As right hand side for this problem, we use DAUX;
C       this is either the original right hand side (for stationary
C       simulations/simulations with steady RHS) or the newly
C       assembled one (for time dependend RHS):
        
        CALL NSDEF2 (NLMIN,NLMAX,CURTRI,
     *               MATDAT,VECDAT,
     *               ISTPAR,DSTPAR,
     *               IMGPAR,DMGPAR,
     *               IASMBL,DASMBL,
     *               NUVP,DUP,DAUX)
        
C       Add the timing information from NSDEF to our current timing
C       statistics

        DO I=0,SZTIMG-1
          DPARAM(OTNSTIM+I) = DPARAM(OTNSTIM+I) + DSTPAR(OTNLTIM+I)
        END DO
        
C       Remember number of nonlinear steps and iterations of the linear
C       solver:

        NNL = ISTPAR(OITE)
        NMG = ISTPAR(ONLIN)
        
        IPARAM(OTITNL) = IPARAM(OTITNL) + NNL
        IPARAM(OTITLI) = IPARAM(OTITLI) + NMG
        
        IF (IER.NE.0) GOTO 99999

C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================

C       Reinterpret the status of the nonlinear solver as result for
C       the timestep:

        STATUS = ISTPAR(OSTATUS)
        
C       If everything was fine, do some post-correction of the pressure
C       before finish
        
        IF (STATUS.EQ.0) THEN

          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,0)

C***********************************************************************

C         Scale back the pressure by the time-step.
C         Remember: p was scaled down by the time-step above,
C         so we have now to scale by the inverse of the time-step
C         to scale it back.

          CALL LLC1(DUP(1+2*NEQU),DUP(1+2*NEQU),NEQP,0D0,
     *              1D0/TIMSTP(OTSTEP))

C***********************************************************************

          CALL GTMAUX (DTIMIN,DPARAM,OTNSTIM-1+OTTLC,1)
          
        END IF

C       Time step completed - but probably we have to repeat it...
C
C       Every time we want to recalculate the current time step,
C       we have to return to our backup...

        IF (IREP.LT.MSREP) THEN
        
          CALL LCP1(DWORK(L(LUBK)),DUP,NUVP)
          
          IF (LRHSBK.NE.0) THEN
            CALL LCP1(DWORK(L(LRHSBK)),DRHS,NUVP)
          END IF
          
          DASMBL(OTIMENS) = TIMEBK
          
        END IF

      END DO ! IREP
      
99999 CONTINUE

C     Release memory that was allocated for backups of grids:

      DO I=NLMAX,NLMIN,-1
        CALL TRIDEL (SOLTRI(1,I))
      END DO

C     Restore the original discretisation setup
      
      CALL LCP3(IASMBK,IASMBL,SZASMI)
      CALL LCP1(DASMBK,DASMBL,SZASMD)
      
C     Restore the stationary solver configuration
      
      CALL LCP3(ISTBCK,ISTPAR,SZNSDI)
      CALL LCP1(DSTBCK,DSTPAR,SZNSDD)
      
C     Stop the time we need for the time step:

      CALL GTMAUX (DTIMIN,DPARAM,OTMTOT,1)
        
C     Finally, release all temporary arrays we probably 
C     allocated...

      IF (LRHSBK.NE.0) CALL ZDISP(0,LRHSBK,'DRHSBK')
      IF (LUBK.NE.0)   CALL ZDISP(0,LUBK,'DUBK  ')

      END

************************************************************************
* Postprocessing of time step
*
* This routine calculates some statistical data based on the
* result of a call to MGSTP2 and writes them to screen.
*
* In:
*   VECDAT : array [1..SZN2VI] of integer 
*            TNS2DVectorParams-structure,
*            defines the form of the RHS vector.
*   DUP    - array [1..*] of double
*            Solution vector
*   ITNS   - Number of current time step
*   ICORST - Number of current substep in Theta scheme
*   TIMENS - Current simulation time
*
************************************************************************

      SUBROUTINE MGSPOS (VECDAT,DUP,ITNS,ICORST,TIMENS)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'
      
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      
      INTEGER ITNS,ICORST,VECDAT(SZN2VI)
      DOUBLE PRECISION DUP(*),TIMENS
      
      DOUBLE PRECISION DSXN,RELU2,RELP2
      
C     local variables

      INTEGER NEQU,NEQP

C     Get the vector size of velocity and pressure part:

      NEQU = VECDAT(ONU)
      NEQP = VECDAT(ONP)

C     Calculate norm of U and P:      
      
      CALL LL21(DUP,2*NEQU,DSXN)
      RELU2=DSXN/SQRT(DBLE(NEQU))
      CALL LL21(DUP(1+2*NEQU),NEQP,DSXN)
      RELP2=DSXN/SQRT(DBLE(NEQP))
      
      IF (MT.GE.2) WRITE(MTERM,3)
      IF (MT.GE.2) WRITE(MTERM,10002) ITNS,ICORST,TIMENS,RELU2,RELP2
      IF (MT.GE.2) WRITE(MTERM,3)
      IF (MT.GE.2) WRITE(MTERM,*)

      IF (MT.GE.1) WRITE(MFILE,3)
      IF (MT.GE.1) WRITE(MFILE,10002) ITNS,ICORST,TIMENS,RELU2,RELP2
      IF (MT.GE.1) WRITE(MFILE,3)
      IF (MT.GE.1) WRITE(MFILE,*)

   1  FORMAT(79('-'))
   3  FORMAT(79('+'))
1000  FORMAT (6E12.5)
1001  FORMAT(' IT DIV-L2',3X,'RHOMGP')
1003  FORMAT(I4,2(D9.2))
10002 FORMAT ('#',I4,'(',I4,')',1X,'TIME=',D10.3,2X,'NORM(U)=',
     *        D14.7,2X,'NORM(P)=',D14.7)

      END
      
