************************************************************************
* The routines in this file capsule the solver for the stationary
* incomporessible Navier Stokes equations.
*
* For solving linear sub-problems, NSDEF calls a linear sub-solver
* in NSDEFLINSOL, which typically makes use of the multigrid solver 
* M020. The corresponding callback routines of M020 are defined
* in NSDEFMGROUT.F.
************************************************************************

      SUBROUTINE  NSDEF2 (NLMIN,NLMAX,TRIAS,
     *                    MATDAT,VECDAT,
     *                    IPARAM,DPARAM,
     *                    IMGPAR,DMGPAR,
     *                    IASMBL,DASMBL,
     *                    NUVP,DUP,DRHS)  

************************************************************************
* Solver for the stationary incompressible Navier Stokes.
* Fixed point defect correction method, nonlinear iteration.
* Multigrid as solver for linear auxiliary Oseen problems.
* Nonlinear parameter optimization for the correction from the linear
* solver step.
*
* Extended calling convention.
*
* This routine calculates a solution DU = (u,v,p)
* for the  stationary incompressible Navier Stokes by nonlinear
* iteration. DU is used as a start vector of the iteration and will be
* overwritten by the new iterate.
*
* The routine is configured by the information in the double/integer
* parameter blocks IPARAM/DPARAM as described in SNSDEF.INC.
*
* In:
*   NLMIN  : minimum level in MATDAT/VECDAT that is filled with data
*   NLMAX  : maximum level in MATDAT/VECDAT that is filled with data
*   TRIAS  : array [1..SZTRIA,1..NNLEV] of integer
*            Triangulation structures for the underlying meshes
*            on all levels. TRIAS(.,ILEV) represents the triangulation
*            structure on level ILEV.
*   MATDAT : array {1..SZN2MI,1..NNLEV] of integer
*            TNS2DMatrixParams-structures for level NLMIN..NLMAX, 
*            initialized with data. MATDAT(1..NLMAX).LA1 must hold
*            handles to a temporary space, which is used to create
*            linear approximations to the nonlinear system matrices
*            during the iteration.
*   VECDAT : array {1..SZN2VI,1..NNLEV] of integer
*            TNS2DVectorParams-structures for level NLMIN..NLMAX. 
*            This structure array must specify the structure of
*            the vectors on each level. Furthermore it must contain
*            handles to preallocated arrays. These arrays are used
*            as auxiliary arrays for solving a linear system on each 
*            level.
*            VECDAT is used as auxiliary structure when solving
*            the system on the finest level with the Multigrid
*            solver.
*
*   IPARAM : array [1..SZNSDI] of integer
*            Integer parameter block with input and output parameters.
*            The input-variables of the array must be initialized, e.g.
*            with INISTS.
*   DPARAM : array [1..SZNSDD] of double
*            Double precision parameter block input and output
*            parameters.
*            The input-variables of the array must be initialized, e.g.
*            with INISTS.
*   IMGPAR : array [1..SZ020I+3*SZSLVI] of integer
*   DMGPAR : array [1..SZ020D+3*SZSLVD] of double
*            Integer and double parameter blocks for the multigrid 
*            sub-solver M020. 
*            The proble-specific variables in this structure will 
*            initialized by NSDEF,
*                        KOFFx, KNEQ, NLMIN, NLMAX
*            and reset when the routine is finished.
*            The general solver parameters, e.g. smoother, coarse 
*            grid solver,... as well as the number of smoothing steps 
*            on each level must be initialized by the caller, e.g. 
*            with INMGST.
*            NSDEF will frequently call the multigrid solver with 
*            this structure to calculate solutions for auxiliary 
*            problems.
*            The data in these structures must correspond to
*            the rest of the parameters of this routine,
*            namely VECDAT
*            
*            xMGPAR contains two different substructures:
*            1.) The first part of the structure [1..SZ020x] 
*                represents the structure of the main multigrid 
*                solver and defines its behaviour
*            2.) The second part of the structure defines the coarse
*                grid solver. Up to now, only one coarse grid solver
*                is available (VANCA). Nevertheless, the structure has
*                to be completely initialized like for every solver,
*                as it defines the behaviour of the coarse grid solver
*                (number of iterations,...)
*   IASMBL : array [1..SZASMI] of integer
*   DASMBL : array [1..SZASMD] of double
*            Integer and double prec. parameter block that controls the
*            discretization. This tells all assembly-routines how to 
*            set up the nonlinearity in the system matrices, which 
*            cubature formula to use, etc.
*
*   NUVP   : total length of solution vector; usually 2*NEQU+NEQP
*   DUP    : array [1..NUVP] of double
*            Start vector for nonlinear iteration.
*            = u^n in a time-dependent simulation.
*            Remark: DUP must not coincide with VECDAT[NLMAX].LSOL !
*   DRHS   : array [1..NUVP] of double
*            Right hand side for the nonlinear equation
*            Remark: DRHS must not coincide with VECDAT[NLMAX].LRHS !
*  
* Out:
*   IPARAM,
*   DPARAM : The "output" part is changed accordingly to the outcome
*            of the algorithm
*   DUP    : receives the new approximation to the solution vector
*
*   The memory for the system matrices, referenced by the handles 
*   MATDAT(1..NLMAX).LA1, is used as temporary space for the
*   calculation of the linear representations of the nonlinear system
*   matrices on every level in every iteration.
*   The content of these arrays is undefined when this routine finishes.
************************************************************************

      IMPLICIT NONE
      
C main COMMON blocks
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'ssolvers.inc'
      
      INCLUDE 'sassembly.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'
      
      INCLUDE 'stiming.inc'
      INCLUDE 'snsdef.inc'
      
C parameters
      
      INTEGER IPARAM (SZNSDI)
      DOUBLE PRECISION DPARAM (SZNSDD),DMGPAR(*),DASMBL(SZASMD)
      DOUBLE PRECISION DUP(*)
      DOUBLE PRECISION DRHS(*)
      INTEGER TRIAS(SZTRIA,NNLEV)
      INTEGER MATDAT(SZN2MI,NNLEV),VECDAT(SZN2VI,NNLEV),NLMIN,NLMAX
      INTEGER IMGPAR(*),IASMBL(SZASMI)
      
C local variables

      INTEGER NITER,STATUS,NEQU,NEQP
      DOUBLE PRECISION RHOLI

C     Allocate a 2nd double parameter structure - this one is
C     only used for measuring timing information!

      DOUBLE PRECISION DTIMIN(SZNSDD),DXOTIM(5)
      INTEGER KVERT,KMID,KCORVG,KSCNPR,NBDMT,KD,KAUX,KX,NUVP
      INTEGER KA1, KCOLA, KLDA, IELT, NA, ILEV, MSHOW
      
C     Some temporary variables to transfer information between levels:

      INTEGER KA1X, KCOLAX, KLDAX, KVERTX, KMIDX, KCORVX, KADJX
      INTEGER KU1X, KU2X, KD1X, KD2X, NELX, NVTX, NEQAX, NEQBX
      INTEGER KSCNPX, NBDMTX, NAX, NEQX

      INTEGER KA1F, KCOLAF, KLDAF, KVERTF, KMIDF, KCORVF, KADJF
      INTEGER KU1F, KU2F, KD1F, KD2F, NELF, NVTF, NEQAF, NEQBF
      INTEGER NAF,  NEQF

      DOUBLE PRECISION UPSAM,RE,THWEIG,TTMG

      LOGICAL BNLEND

      INTEGER INL,I1,KGRVEL,KGRVLX
      DOUBLE PRECISION RESU, RESDIV, DELP
      DOUBLE PRECISION RHO,DELU
      LOGICAL BMG
      DOUBLE PRECISION RESOLD,RES0,RES,EPSRES,OMEGA
      
      CHARACTER CFN*60
      
C     Timing for linear solver

      DOUBLE PRECISION TIMGLS (9)

C     definition of finite elements

      EXTERNAL E030,E031,EM30,EM31

C=======================================================================
C     Initialization
C=======================================================================

C     Clear timing and status information
      
      CALL LCL1(DPARAM(OTNLTIM),SZTIMG)
      
      IPARAM(ONLIN) = 0
      
      IPARAM(OSTATUS) = 0

C     The TIMGLS-array calculates the total time for the different
C     solver components...

      CALL LCL1(TIMGLS,9)

C     Stop the complete time of this algorithm, save into
C     the timing structure.
      
      CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTNL,0)
      
C     Save some array variables to local variables for easier access,
C     especially information about our current system matrix.

      UPSAM  = DASMBL(OUPSAM)
      RE     = DASMBL(ORE)
      THWEIG = DASMBL(OTHWEIG)
      
      IF (THWEIG.EQ.0D0) THEN
        WRITE (MTERM,'(A)') 'NSDEF: Error; THETA=0.0 is not allowed!'
        RETURN
      END IF
      
      IELT   = IASMBL(OIELEMT)

      KCORVG = L(TRIAS(OLCORVG,NLMAX))
      KMID   = L(TRIAS(OLMID,NLMAX))
      KVERT  = L(TRIAS(OLVERT,NLMAX))
      KCOLA  = L(MATDAT(OLCLA1,NLMAX))
      KLDA   = L(MATDAT(OLLDA1,NLMAX))
      NA     = MATDAT(ONA1,NLMAX)
      
C     Get the vector size of velocity and pressure part:

      NEQU = VECDAT(ONU,NLMAX)
      NEQP = VECDAT(ONP,NLMAX)

C     If the ALE method is used, store the address of the grid velocoity
C     vector to KGRVEL. Otherwise store 1 into KGRVEL, which is used
C     as a dummy in the call to the upwinding routines to prevent
C     specifying a DWORK(.) address out of bounds.

      IF (IASMBL(OIALE).NE.0) THEN
        KGRVEL = L(IASMBL(OLGRVEL))
      ELSE
        KGRVEL = 1
      END IF

C     Use the shortcut nodal property array to decide which points are
C     on the boundary, which are Neumann boundary,...

      KSCNPR = L(TRIAS(OTRIUD+1,NLMAX))
      NBDMT  = TRIAS(OTRIUD,NLMAX)

C     The system matrix LA1 on the finest level is a temporary memory
C     block that is frequently overwritten during the nonlinear
C     iteration. As the system matrix is nonlinear, anyway there exists
C     no definite (linear) representation...
C     The matrix can be accessed in DWORK at:

      KA1    = L(MATDAT(OLA1,NLMAX))
 
C     Message level

      MSHOW = IPARAM (OMSGTRM)
      
C     Save the position of the RHS-vector on the finest level into KD.
C     The defect vector will be build up there, then the MG solver
C     will be used to precondition that vector, and finally the
C     solution of the MG solver will be used to update our current
C     solution DU.

      KD = L(VECDAT(OLRHS,NLMAX))
C     Also save the position of the auxiliary vector and of the
C     solution vector on finest level to KAUX and KX:
      
      KAUX = L(VECDAT(OLTMP,NLMAX))
      KX = L(VECDAT(OLSOL,NLMAX))
      
C=======================================================================
C     First generation of the nonlinear block A on level NLMAX
C     and generation of the defect vectors.
C=======================================================================

C     BNLEND will be set to true if the nonlinear solver converged
C     successfully

      BNLEND = .FALSE.

C     Copy the current RHS to the vector DD, which is used to 
C     hold the defect during the nonlinear iteration.

      CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTLC,0)
      CALL LCP1 (DRHS,DWORK(KD),NUVP)
      CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTLC,1)
      
C     Now calculate the linear part of the system matrix as well as
C     the defect vector for velocity/pressure arising from the linear
C     parts:
C
C      -->   KA1 = [  ALPHA*M*I  +  THWEIG * (-nu * Laplace(.))  ]
C      -->   (KD1) = (KD1)                                               ( KU1 )
C            (KD2) = (KD2) - [ ALPHA*M*u^n + THWEIG*(-nu*Laplace(u^n)) ] ( KU2 )
C            (KDP) = (KDP)                                               ( KP  )
C
C     with ALPHA=ISTAT. The matrix KA1 will later also serve as a
C     preconditioner in the Oseen equation. The nonlinear part
C     is subtracted later...

      CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTADF,0)
      CALL LCL1(DWORK(KA1+2750),250)
      CALL XMADF4(TRIAS(1,NLMAX),MATDAT(1,NLMAX),DUP,DWORK(KD),
     *            IPARAM(ONITMAX)-1,IASMBL(OIPRECA),
     *            IASMBL(OINEUM),THWEIG,IASMBL(OIALPHA),NLMAX)
      CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTADF,1)

C     Check if we have Stokes flow or Navier Stokes flow.
C
C     When calculating Navier Stokes,
C
C         du/dt + nu*Laplace(u) + u*grad(u) + grad(p) = f
C
C     we have to treat the nonlinearity u*grad(u). In contrast when
C     treating Stokes flow ("slow flow"),
C
C         du/dt + nu*Laplace(u) + grad(p) = f
C
C     the nonlinear term can be neglected and so we have nothing to do
C     there to treat it.
C
C     Only exception: Stokes Flow with ORCSMT=1=build always everything
C     without relying on precalculated data in the KST1-matrices.
      
      CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTUPW,0)
      
      IF ((IASMBL(OISTOK).NE.1).OR.(IASMBL(OIPRECA).NE.4)) THEN

C       Add the convective part u*grad(u) to the system matrix KA1 -
C       the linear Oseen matrix that will later serve as a
C       preconditioner.
C
C       Up to now our system matrix has consisted only of linear
C       terms:
C
C           KA1 = [  M*I  +  THWEIG * (-nu * Laplace(.))  ]
C
C       but this is not enough... the nonlinear part is still missing.
C       So we have to add the following term to KA1:
C
C           THWEIG * u grad(.)
C
C       what will finally result in the system matrix
C
C           KA1 = [  M*I  +  THWEIG * (-nu * Laplace(.))  ] + THWEIG * u grad(.)
C               = [  M*I  +  THWEIG * N(u) ]
C
C       Add the corresponding nonlinearity to the defect vectors, 
C       what is still missing.
C
C       This overall procedure then results in a matrix for the
C       "linearizes Navier Stokes" equation, i.e. the Oseen equation.
C       This matrix is later used as a preconditioner for the defect
C       vector (cf. p. 168, Turek's book).
C
C       Additionally to setting up the nonlinear part, also update
C       the defect vector; subtract (Nonliner Part)*(Solution) to
C       incorporate the nonlinear defect according to
C
C      -->   (KD1) = (KD1)                 ( KU1 )
C            (KD2) = (KD2) - [ u grad(.) ] ( KU2 )
C            (KDP) = (KDP)                 ( KP  )
        
        IF (IASMBL(OIUPW).EQ.1) THEN
          CALL GUPWDX (DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
     *               1D0,0D0,DUP(1),DUP(1+NEQU),
     *               DWORK(KD),DWORK(KD+NEQU),
     *               DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *               TRIAS(1,NLMAX),KWORK(KVERT),KWORK(KMID),
     *               DWORK(KCORVG),1,
     *               UPSAM,RE,THWEIG,
     *               IASMBL(OIALE),DWORK(KGRVEL))
        ELSE
          IF (IELT.EQ.0) THEN 
C            CALL SUPGPX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
C     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
C     *                  DWORK(KD),DWORK(KD+NEQU),
C     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
C     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
C     *                  TRIAS(1,NLMAX),E031,IASMBL(OICUBN),
C     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS), 
C     *                  IASMBL(OISTOK), UPSAM, RE,  
C     *                  1,1D0,THWEIG)
            CALL SUPAPX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
     *                  DWORK(KD),DWORK(KD+NEQU),
     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
     *                  TRIAS(1,NLMAX),E031,IASMBL(OICUBN),
     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
     *                  IASMBL(OISTOK),UPSAM, RE,
     *                  1,1D0,THWEIG,
     *                  IASMBL(OIALE),DWORK(KGRVEL))
          END IF
     
          IF (IELT.EQ.1) THEN 
C            CALL SUPGPX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
C     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
C     *                  DWORK(KD),DWORK(KD+NEQU),
C     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
C     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
C     *                  TRIAS(1,NLMAX),E030,IASMBL(OICUBN),
C     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
C     *                  IASMBL(OISTOK),UPSAM, RE,  
C     *                  1,1D0,THWEIG)
            CALL SUPAPX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
     *                  DWORK(KD),DWORK(KD+NEQU),
     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
     *                  TRIAS(1,NLMAX),E030,IASMBL(OICUBN),
     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
     *                  IASMBL(OISTOK),UPSAM, RE,
     *                  1,1D0,THWEIG,
     *                  IASMBL(OIALE),DWORK(KGRVEL))
          END IF
          
          IF (IELT.EQ.2) THEN
C            CALL SUPGNX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
C     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
C     *                  DWORK(KD),DWORK(KD+NEQU),
C     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
C     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
C     *                  TRIAS(1,NLMAX),EM31,IASMBL(OICUBN),
C     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
C     *                  IASMBL(OISTOK),UPSAM, RE,
C     *                  1,1D0,THWEIG)
            CALL SUPANX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
     *                  DWORK(KD),DWORK(KD+NEQU),
     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
     *                  TRIAS(1,NLMAX),EM31,IASMBL(OICUBN),
     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
     *                  IASMBL(OISTOK),UPSAM, RE,
     *                  1,1D0,THWEIG,
     *                  IASMBL(OIALE),DWORK(KGRVEL))
          END IF
          
          IF (IELT.EQ.3) THEN
C            CALL SUPGNX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
C     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
C     *                  DWORK(KD),DWORK(KD+NEQU),
C     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
C     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
C     *                  TRIAS(1,NLMAX),EM30,IASMBL(OICUBN),
C     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
C     *                  IASMBL(OISTOK),UPSAM, RE,
C     *                  1,1D0,THWEIG)
            CALL SUPANX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
     *                  DWORK(KD),DWORK(KD+NEQU),
     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
     *                  TRIAS(1,NLMAX),EM30,IASMBL(OICUBN),
     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
     *                  IASMBL(OISTOK),UPSAM, RE,
     *                  1,1D0,THWEIG,
     *                  IASMBL(OIALE),DWORK(KGRVEL))
          END IF

        ENDIF ! (IUPW.EQ.1)

C       That's it for the nonlinearity. Now we come to the real
C       iteration.

      ENDIF ! ((ISTOK.NE.1).OR.((ISTOK.EQ.1).AND.(IPRECA.EQ.4)))
        
      CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTUPW,1)
   
C=======================================================================
C     Treat Dirichlet-nodes:
C=======================================================================

C     Replace all rows in the system matrix corresponding to Dirichlet
C     nodes by unit vectors:
   
      CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTBDR,0)
      CALL BDRYA (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *            KWORK(KSCNPR),NBDMT)

C     Also set the entries in the defect vector corresponding to
C     Dirichlet nodes to 0 - because as solution and RHS are equal
C     there and the corresponding line in the matrix is a unit vector,
C     the defect is 0 there.

      CALL BDRY0 (DWORK(KD),DWORK(KD+NEQU),KWORK(KSCNPR),NBDMT)
      CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTBDR,0)

C=======================================================================
C     Calculation of initial defects
C=======================================================================

      CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTLC,0)
      
C     Calculate the norm of the initial residuum. The residuum itself
C     was calculated in the assembling-routines above (XMADF4, ...)
C     The results of this routine are saved into RESU,RESDIV.
      
      CALL  RESDFK(DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *             DWORK(KD),DWORK(KD+NEQU),DWORK(KD+2*NEQU),
     *             DRHS(1),DRHS(1+NEQU),DRHS(1+2*NEQU),NEQU,
     *             NEQP,RESU,RESDIV)
     
      RESOLD=SQRT(RESU*RESU+RESDIV*RESDIV)
      RES0  =RESOLD
      RES   =RESOLD
      EPSRES=DPARAM(ODMPD  )*RESOLD
      
C     Some nice output :)

      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)

      IF (MSHOW.GE.2) WRITE(MTERM,1001)
      IF (MSHOW.GE.1) WRITE(MFILE,1001)
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)
      INL=0
      IF (MSHOW.GE.2) WRITE(MTERM,1002)  INL,RESU,RESDIV,RESOLD
      IF (MSHOW.GE.1) WRITE(MFILE,1002)  INL,RESU,RESDIV,RESOLD
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)

C     Initialize the OMEGA-parameter, our damping parameter for the
C     defect correction

      OMEGA = DPARAM(OOMGINI)

      CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTLC,1)

C=======================================================================
C *** Loop of nonlinear iteration
C=======================================================================

C     Now we start with the actual nonlinear iteration to solve the
C     problem. Out problem is to find a solution u that fulfills
C     (cf p. 163 in Turek's book):
C
C         T(u)u = f
C
C     for an appropriate operator T(.) and a solution u. In our context,
C     T(.) is the discrete Navier-Stokes operator and u=(u,v,p) is the
C     velocity/pressure solution vector.
C
C     The general approach for the nonlinear here is the defect
C     correction:
C
C           u^(l+1)  =  u^l  -  OMEGA * C * ( T(u^l)u^l - f )
C
C     with l=1..INLMAX for a start vector u^1 and C being an 
C     appropriate preconditioner. (In fact, as preconditioner
C     we choose C=T(u^l)^{-1} here!) The parameter OMEGA
C     can be prescribed by the user by setting OMGMIN=OMGMAX=anything,
C     or can be calculated automatically when OMGMIN < OMGMAX with
C     the "adaptive fixed point defect correction" method.
C
C     Perform a maximum of INLMAX nonlinear iterations:
      
      DO INL=1,IPARAM(ONITMAX)
      
        IPARAM(OCRITE) = INL
      
C=======================================================================
C *** Generate the A blocks for all coarse levels
C=======================================================================

        IF (NLMAX.GT.NLMIN)  THEN

C         Copy the current solution vector to the initial solution
C         vector of the MG solver - it denotes the initial state of
C         the MG iteration.
C         Below, this vector is restricted to the lower levels
C         and used to build the nonlinear system matrices on the lower
C         levels!

          CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTLC,0)
          KU1F = L(VECDAT(OLSOL,NLMAX))
          CALL LCP1 (DUP,DWORK(KU1F),NUVP)
          CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTLC,1)
       
C         Loop through all levels, looking from the lower level
C         to the upper level always.
        
          DO ILEV=NLMAX-1,NLMIN,-1
          
            CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTLC,0)
            
C           Initialize local variables with information about the
C           current (coarser) level for easier access (so that the
C           calls are not so long in the following ;) )

            KCORVX = L(TRIAS(OLCORVG,ILEV))
            KMIDX  = L(TRIAS(OLMID,ILEV))
            KVERTX = L(TRIAS(OLVERT,ILEV))
            KADJX  = L(TRIAS(OLADJ,ILEV))
            KA1X   = L(MATDAT(OLA1,ILEV))
            KCOLAX = L(MATDAT(OLCLA1,ILEV))
            KLDAX  = L(MATDAT(OLLDA1,ILEV))
            NAX    = MATDAT(ONA1,ILEV)
            NEQAX  = MATDAT(ONEQA,ILEV)
            NEQBX  = MATDAT(ONEQB,ILEV)
            NEQX   = VECDAT(ONEQV,ILEV)
            KU1X   = L(VECDAT(OLSOL,ILEV))
            KU2X   = KU1X+NEQAX
            KD1X   = L(VECDAT(OLTMP,ILEV))
            KD2X   = KD1X+NEQAX
            NELX   = TRIAS(ONEL,ILEV)
            NVTX   = TRIAS(ONVT,ILEV)

C           and some other local variables with information about the
C           finer one:

            I1=ILEV+1

            KCORVF = L(TRIAS(OLCORVG,I1))
            KMIDF  = L(TRIAS(OLMID,I1))
            KVERTF = L(TRIAS(OLVERT,I1))
            KADJF  = L(TRIAS(OLADJ,I1))
            KA1F   = L(MATDAT(OLA1,I1))
            KCOLAF = L(MATDAT(OLCLA1,I1))
            KLDAF  = L(MATDAT(OLLDA1,I1))
            NAF    = MATDAT(ONA1,I1)
            NEQAF  = MATDAT(ONEQA,I1)
            NEQBF  = MATDAT(ONEQB,I1)
            NEQF   = VECDAT(ONEQV,I1)
            KU1F   = L(VECDAT(OLSOL,I1))
            KU2F   = KU1F+NEQAF
            KD1F   = L(VECDAT(OLTMP,I1))
            KD2F   = KD1F+NEQAF
            NELF   = TRIAS(ONEL,I1)
            NVTF   = TRIAS(ONVT,I1)
            
C           The grid velocity vector starts at the same position
C           as the one on the finest level, as the sorting always
C           ensures that the first set of vertex coordinates on the
C           finest level corresponds to the vertex coordinates on
C           the lower levels.
            
            KGRVLX = KGRVEL

C           Restrict the solution from the finer level to the coarser
C           one. Above we initialized the solution vector on the finest
C           level, which is now transported to the lower ones. So we have
C           a representation of the solution on each level.

C            CALL  RESTRU (DWORK(KU1X),DWORK(KU2X),
C     *                    DWORK(KU1F),DWORK(KU2F),
C     *                    KWORK(KVERTX),KWORK(KMIDX),KWORK(KADJX),
C     *                    NEQAX,NEQBX,NVTX,
C     *                    KWORK(KVERTF),KWORK(KMIDF),KWORK(KADJF),
C     *                    NEQAF,NEQBF,NVTF)

            CALL  RESTRU (DWORK(KU1X),DWORK(KU2X),
     *                    DWORK(KU1F),DWORK(KU2F),
     *                    KWORK(KVERTX),KWORK(KMIDX),KWORK(KADJX),
     *                    NEQAX,NELX,NVTX,
     *                    KWORK(KVERTF),KWORK(KMIDF),KWORK(KADJF),
     *                    NEQAF,NELF,NVTF)

            CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTLC,1)

C           Construct the linear part of the nonlinear matrix on the lower
C           level. Use the restricted solution from above to construct
C           the nonlinearity.
C           There is no defect vector to be calculated, we use a 
C           routine here that only builds the nonlinearity of the 
C           KA1-matrix.

            CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTADF,0)
            CALL XMADF5(MATDAT(1,ILEV),IASMBL(OIPRECA),THWEIG,
     *                  IASMBL(OIALPHA),ILEV)
            CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTADF,1)

C           In the case that we have to calculate Navier Stokes (ISTOK=0)
C           or that we have a Stokes calculation with IPRECA=4,
C           add the convective part u*grad(u) to the system matrix KA1.
C
C           Up to now our system matrix has consisted only of linear
C           terms:
C
C               KA1 = [  M*I  +  THWEIG * (-nu * Laplace(.))  ]
C
C           but this is not enough... the nonlinear part is still missing.
C           So we have to add the following term to KA1:
C
C               THWEIG * u grad(.)
C
C           what will finally result in the system matrix
C
C               KA1 = [  M*I  +  THWEIG * (-nu * Laplace(.))  ] + THWEIG * u grad(.)
C                   = [  M*I  +  THWEIG * N(u) ]

            CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTUPW,0)
            
            IF ((IASMBL(OISTOK).NE.1).OR.(IASMBL(OIPRECA).NE.0)) THEN
              IF (IASMBL(OIUPW).EQ.1) THEN
               CALL GUPWDX (DWORK(KU1X),DWORK(KU2X),
     *                     DWORK(KU1X),DWORK(KU2X),
     *                     1D0,0D0,DWORK(KU1X),DWORK(KU2X),
     *                     DWORK(KD1X),DWORK(KD2X),
     *                     DWORK(KA1X),KWORK(KCOLAX),KWORK(KLDAX),
     *                     TRIAS(1,ILEV),
     *                     KWORK(KVERTX),KWORK(KMIDX),
     *                     DWORK(KCORVX),0,
     *                     UPSAM,RE,THWEIG,
     *                     IASMBL(OIALE),DWORK(KGRVLX))
              ELSE
              
               IF (IELT.EQ.0) THEN
C                 CALL SUPGPX(NEQAX,DWORK(KU1X),DWORK(KU2X),
C     *                  DWORK(KU1X),DWORK(KU2X),1D0,0D0,
C     *                  DWORK(KU1X),DWORK(KU2X),
C     *                  DWORK(KD1X),DWORK(KD2X),
C     *                  DWORK(KA1X),NAX,KWORK(KCOLAX),KWORK(KLDAX),
C     *                  KWORK(KVERTX),KWORK(KMIDX),DWORK(KCORVX),
C     *                  TRIAS(1,ILEV),E031,IASMBL(OICUBN),
C     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS), 
C     *                  IASMBL(OISTOK), DASMBL(OUPSAM), RE,  
C     *                  0,1D0,THWEIG)
                 CALL SUPAPX(NEQAX,DWORK(KU1X),DWORK(KU2X),
     *                  DWORK(KU1X),DWORK(KU2X),1D0,0D0,
     *                  DWORK(KU1X),DWORK(KU2X),
     *                  DWORK(KD1X),DWORK(KD2X),
     *                  DWORK(KA1X),NAX,KWORK(KCOLAX),KWORK(KLDAX),
     *                  KWORK(KVERTX),KWORK(KMIDX),DWORK(KCORVX),
     *                  TRIAS(1,ILEV),E031,IASMBL(OICUBN),
     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS), 
     *                  IASMBL(OISTOK), DASMBL(OUPSAM), RE,  
     *                  0,1D0,THWEIG,
     *                  IASMBL(OIALE),DWORK(KGRVLX))
               END IF

               IF (IELT.EQ.1) THEN
C                 CALL SUPGPX(NEQAX,DWORK(KU1X),DWORK(KU2X),
C     *                  DWORK(KU1X),DWORK(KU2X),1D0,0D0,
C     *                  DWORK(KU1X),DWORK(KU2X),
C     *                  DWORK(KD1X),DWORK(KD2X),
C     *                  DWORK(KA1X),NAX,KWORK(KCOLAX),KWORK(KLDAX),
C     *                  KWORK(KVERTX),KWORK(KMIDX),DWORK(KCORVX),
C     *                  TRIAS(1,ILEV),E030,IASMBL(OICUBN),
C     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS), 
C     *                  IASMBL(OISTOK), DASMBL(OUPSAM), RE,  
C     *                  0,1D0,THWEIG)
                 CALL SUPAPX(NEQAX,DWORK(KU1X),DWORK(KU2X),
     *                  DWORK(KU1X),DWORK(KU2X),1D0,0D0,
     *                  DWORK(KU1X),DWORK(KU2X),
     *                  DWORK(KD1X),DWORK(KD2X),
     *                  DWORK(KA1X),NAX,KWORK(KCOLAX),KWORK(KLDAX),
     *                  KWORK(KVERTX),KWORK(KMIDX),DWORK(KCORVX),
     *                  TRIAS(1,ILEV),E030,IASMBL(OICUBN),
     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS), 
     *                  IASMBL(OISTOK), DASMBL(OUPSAM), RE,  
     *                  0,1D0,THWEIG,
     *                  IASMBL(OIALE),DWORK(KGRVLX))
               END IF

               IF (IELT.EQ.2) THEN 
C                 CALL SUPGNX(NEQAX,DWORK(KU1X),DWORK(KU2X),
C     *                  DWORK(KU1X),DWORK(KU2X),1D0,0D0,
C     *                  DWORK(KU1X),DWORK(KU2X),
C     *                  DWORK(KD1X),DWORK(KD2X),
C     *                  DWORK(KA1X),NAX,KWORK(KCOLAX),KWORK(KLDAX),
C     *                  KWORK(KVERTX),KWORK(KMIDX),DWORK(KCORVX),
C     *                  TRIAS(1,ILEV),EM31,IASMBL(OICUBN),
C     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS), 
C     *                  IASMBL(OISTOK), DASMBL(OUPSAM), RE,  
C     *                  0,1D0,THWEIG)
                 CALL SUPANX(NEQAX,DWORK(KU1X),DWORK(KU2X),
     *                  DWORK(KU1X),DWORK(KU2X),1D0,0D0,
     *                  DWORK(KU1X),DWORK(KU2X),
     *                  DWORK(KD1X),DWORK(KD2X),
     *                  DWORK(KA1X),NAX,KWORK(KCOLAX),KWORK(KLDAX),
     *                  KWORK(KVERTX),KWORK(KMIDX),DWORK(KCORVX),
     *                  TRIAS(1,ILEV),EM31,IASMBL(OICUBN),
     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS), 
     *                  IASMBL(OISTOK), DASMBL(OUPSAM), RE,  
     *                  0,1D0,THWEIG,
     *                  IASMBL(OIALE),DWORK(KGRVLX))
               END IF

               IF (IELT.EQ.3) THEN 
C                 CALL SUPGNX(NEQAX,DWORK(KU1X),DWORK(KU2X),
C     *                  DWORK(KU1X),DWORK(KU2X),1D0,0D0,
C     *                  DWORK(KU1X),DWORK(KU2X),
C     *                  DWORK(KD1X),DWORK(KD2X),
C     *                  DWORK(KA1X),NAX,KWORK(KCOLAX),KWORK(KLDAX),
C     *                  KWORK(KVERTX),KWORK(KMIDX),DWORK(KCORVX),
C     *                  TRIAS(1,ILEV),EM30,IASMBL(OICUBN),
C     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS), 
C     *                  IASMBL(OISTOK), DASMBL(OUPSAM), RE,  
C     *                  0,1D0,THWEIG)
                 CALL SUPANX(NEQAX,DWORK(KU1X),DWORK(KU2X),
     *                  DWORK(KU1X),DWORK(KU2X),1D0,0D0,
     *                  DWORK(KU1X),DWORK(KU2X),
     *                  DWORK(KD1X),DWORK(KD2X),
     *                  DWORK(KA1X),NAX,KWORK(KCOLAX),KWORK(KLDAX),
     *                  KWORK(KVERTX),KWORK(KMIDX),DWORK(KCORVX),
     *                  TRIAS(1,ILEV),EM30,IASMBL(OICUBN),
     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS), 
     *                  IASMBL(OISTOK), DASMBL(OUPSAM), RE,  
     *                  0,1D0,THWEIG,
     *                  IASMBL(OIALE),DWORK(KGRVLX))
               END IF

              ENDIF
            ENDIF
            CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTUPW,1)

C           Call matrix restriction routine to rebuild the matrix
C           entrys belonging to anisotropic cells. This procedure is 
C           applied to all matrices except for the finest level

            IF (ILEV.NE.NLMAX) THEN
              CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTUPW,0)
              
              CALL IMPMRX (TRIAS(1,ILEV),MATDAT(1,ILEV),
     *                     TRIAS(1,ILEV+1),MATDAT(1,ILEV),
     *                     IASMBL(OIAPRM),DASMBL(ODMTEP))

              CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTUPW,1)
            END IF

C           Because the coarse grid matrix is constructed with the 
C           help of the fine grid matrix, the implementation of the 
C           boundary conditions into the matrix can't be done here! 
C           This is done later in a second step...  

          END DO

C         Implement Dirichlet boundary conditions into the matrices on
C         all levels. Replace the rows corresponding to Dirichlet nodes
C         by unit vectors.

          DO ILEV=NLMAX-1,NLMIN,-1
          
            CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTBDR,0)
            
            KA1X   = L(MATDAT(OLA1,ILEV))
            KCOLAX = L(MATDAT(OLCLA1,ILEV))
            KLDAX  = L(MATDAT(OLLDA1,ILEV))
            KSCNPX = L(TRIAS(OTRIUD+1,ILEV))
            NBDMTX = TRIAS(OTRIUD,ILEV)
            
            CALL BDRYA (DWORK(KA1X),KWORK(KCOLAX),KWORK(KLDAX),
     *                  KWORK(KSCNPX),NBDMTX)

            CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTBDR,1)
            
          END DO

        ENDIF

C=======================================================================
C     Ok, small summary - what have we calculated up to now and where?
C     (cf. p. 164 in Turek's book)
C
C     DWORK(KD)    ->  currently contains d^l = T(u^l)-f,
C                      was calculated with XMADF4 at the beginning.
C     KF1,KF2,KFP  ->  Right hand side of the equation, was
C                      provided by the caller
C     KST1,KB1,KB2 ->  Matrix for the Oseen equation = linearized
C                      Navier Stokes
C                      (denoted as T~(u^l) in the book).
C                      This is now used as a preconditioner for
C                      the defect(=update) vector.
C                      The matrices are available through the MATDAT()-
C                      structure.
C
C     Remember, the basic nonlinear iteration, given by
C
C       u^(l+1)  =  u^l  -  OMEGA * T~(u^l)^(-1) * ( T(u^l)u^l - f )
C
C     - with the preconditioner C = T~(u^l)^{-1} - can be decomposed 
C     into three steps:
C
C      1.) Calculate the defect (done above):
C             d^l = T(u^l) u^l - f
C
C      2.) Solve an auxiliary problem (Oseen equation):
C             T~(u^l) y^l = d^l
C
C      3.) Update the solution:
C             u^(l+1) = u^l - Omega * y^l
C
C     So we are now at step 2 - we must solve the Oseen equation
C     to apply the preconditioner T~(u^l)^{-1} to the defect vector.
C
C     Remark: The precontidioner T~(u^l) was set up the same way
C      as the system matrix T(u^l) that was used for the defect
C      correction. Theoretically (cf. p. 168) the preconditioner
C      can be
C
C       [ ALPHA*M + nu*L + K~(u^l)                          B1 ] ^ {-1}
C   C = [                          ALPHA*M + nu*L + K~(u^l) B2 ]
C       [        B1^T                      B2^T             0  ]
C
C      with K~ not necessarily build up the same way as K in the system
C      matrix
C
C       [ ALPHA*M + nu*L + K(u^l)                          B1 ] ^ {-1}
C   T = [                         ALPHA*M + nu*L + K~(u^l) B2 ]
C       [        B1^T                     B2^T             0  ]
C
C       [ KA1       B1 ]
C     = [      KA1  B1 ]
C       [ B1^T B2^T 0  ]
C
C      For T, K must be as exact as possible, so theat the approximate
C      solution converges to the correct solution - but not for C.
C      Here we theoretically could use another scheme for setting up
C      the convective part.
C      Nevertheless, K~ := K is the natural choice - and as we already
C      built it for the defect correction, we can also use it here,
C      also to save some time.
C=======================================================================

C=======================================================================
C *** Perform MG for Oseen equations
C=======================================================================

C       At this point, the linearized nonlinear system matrix
C       is set up, so we have an stationary Oseen-equation here:
C
C       [ ALPHA*M + THETA*(-nu Laplace(.)) + u^n*grad(.) ] y = (KD1,KD2,KDP)
C                                            ^^^^^^^^^^^   ^
C                                          Linearization   to calc.
C
C       This is just a linear system of the form Au=b, which we can
C       solve with our multigrid solver - so call that.
C
C       Note that solving the linear system will destroy our current
C       solution. Fortunately we normally make a backup of that 
C       (in the if-clause (LU1OLD.NE.0) above), so we can restore it
C       later for the defect correction.
C
C       Set the initial solution vector completely to 0. This expecially
C       sets all the entries corresponding to Dirichlet components to 0.
C       These entries must not be changed from 0 to anything else
C       during the solution process, as the "solution" vector here is
C       a defect vector which will later be added to the real solution!

        CALL LCL1 (DWORK(KX),NUVP)
        
C       Remember: Our defect vector DWORK(KD) is prepared as RHS for the
C       linear Oseen equation. Now we can call the muligrid solver to
C       precondition that vector.
C
C       Solve the Oseen equation to obtain ( Y1, Y2, YP )^T.
C       This is saved in the solution vectors KU1,KU2,KP on the
C       finest level.
C       This will return a solution in VECDAT[NLMAX].SOL
C       as well as NITER, STATUS and RHOLI.
C     
C       NSLINS must transmit the parameter blocks IPARAM/DPARAM to any
C       callback routines to allow them accessing our current
C       status, matrix configuration and so on.
C
C       The callback routine will add timing information to TIMGLS
C       about how much time is needed by which solver component.
C       At the end of NSDEF, TIMGLS then contains the total time
C       for all the components.

        CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTLSOL,0)
        
        CALL NSLINS (NLMIN,NLMAX,TRIAS,
     *               MATDAT,VECDAT,
     *               IPARAM,DPARAM,
     *               IMGPAR,DMGPAR,
     *               IASMBL,DASMBL,
     *               INL,NITER,STATUS,RHOLI,TIMGLS)
        
        CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTLSOL,1)
        
C       This overwrote KU1/KU2/KUP with the update vector,
C       arising from solving Oseen with the defect as a RHS.
C       The vector Y=(KU1,KU2,KP) has to be added to our solution 
C       vector. The current solution vector is still saved in the
C       array of LU1OLD/LU2OLD/LPOLD.

C       Did the MG-solver broke down? How many iterations did we need?

        BMG = STATUS.NE.0
        IPARAM(ONLIN) = IPARAM(ONLIN) + NITER
               
C       That's it, we have our new iterate u^new inside of the
C       nonlinear iteration. If INLMAX=1, this is even our
C       solution vector u^(n+1).

C=======================================================================
C    End of multigrid iteration
C=======================================================================

C=======================================================================
C *** Calculate the optimal correction
C=======================================================================

C       The multigrid solver returns on the finest level the
C       preconditioned correction vector Y=(YU,YV,YP). This is
C       accessible in DWORK at position KOFFX(NLMAX)+1.
C
C       We now want to use that vector to correct our iterate
C       DUP=(DU,DV,DP) using the formula
C
C                      DUP = DUP - OMEGA*Y
C
C       The following routine will calculate and return an optimal
C       OMEGA and correct our solution vector for the next nonlinear
C       iteration. Update the system matrix KA1 on the finest level
C       for the next nonlinear iteration.
C
C       Use the defect vector of the MG structure on the maximum
C       level as temporary space - this is completely save as the state
C       is undefined on exit of multigrid, so the vector is unused.
C       Note that DWORK(KD) and the RHS of the MG iteration on the max.
C       level coincide, so we must use DWORK(KAUX) as temporary vector.

        CALL XOPTCX (TRIAS(1,NLMAX), MATDAT(1,NLMAX), 
     *               NUVP, DWORK(KX), DUP, DRHS, 
     *               DWORK(KD), DWORK(KAUX), 
     *               THWEIG, IASMBL(OIALPHA), NLMAX,
     *               IASMBL,DASMBL,
     *               DPARAM(OOMGMIN), DPARAM(OOMGMAX), 
     *               OMEGA, DELU, DELP, DXOTIM)
     
        DPARAM(OTNLTIM-1+OTTLC)  = DPARAM(OTNLTIM-1+OTTLC)  + DXOTIM(1)
        DPARAM(OTNLTIM-1+OTTADF) = DPARAM(OTNLTIM-1+OTTADF) + DXOTIM(2)
        DPARAM(OTNLTIM-1+OTTUPW) = DPARAM(OTNLTIM-1+OTTUPW) + DXOTIM(3)
        DPARAM(OTNLTIM-1+OTTBDR) = DPARAM(OTNLTIM-1+OTTBDR) + DXOTIM(4)
        DPARAM(OTNLTIM-1+OTTCOR) = DPARAM(OTNLTIM-1+OTTCOR) + DXOTIM(5)
     
C       The LD1,LD2,LDP-vectors are now written with undefined data.
C       Now we have our new iterate (KU1,KU2,KP)=u^(l+1) - hopefully...

C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================

        IF (BMG.AND.(ABS(OMEGA).GE.1D-2).AND.
     *      (IASMBL(OIALPHA).EQ.1)) THEN
          WRITE (MTERM,'(A)') 'Calculation canceled, MG-solver broke '//
     *                        'down'
          WRITE (MFILE,'(A)') 'Calculation canceled, MG-solver broke '//
     *                        'down'
          IPARAM(OSTATUS) = 2
        END IF

C=======================================================================
C *** Calculation of defects
C=======================================================================

C       Copy the RHS of the system to the defect vector KD for the
C       computation. After the call to XOPTCX above, the content of
C       KD is no more needed...
        
        CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTLC,0)
        
        CALL LCP1 (DRHS,DWORK(KD),NUVP)
        
        CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTLC,1)

C       Build the linear part of the system matrix into KST1.
C       Calculate the linear part of the defect vector into LD1:

        CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTADF,0)
        CALL XMADF4(TRIAS(1,NLMAX),MATDAT(1,NLMAX),DUP,DWORK(KD),
     *              IPARAM(ONITMAX)-1,IASMBL(OIPRECA),
     *              IASMBL(OINEUM),THWEIG,IASMBL(OIALPHA),NLMAX)
        CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTADF,1)

C       Incorporate the nonlinearity into the system matrix:
C
C           KA1 = [  ALPHA*M*I  +  THWEIG * (-nu * Laplace(.))  ] + THWEIG * u grad(.)
C               = [  ALPHA*M*I  +  THWEIG * N(u) ]
C
C       (with ALPHA=ISTAT)
C
C       Incorporate the nonlinearity into the defect vector.

        CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTUPW,0)
        
        IF ((IASMBL(OISTOK).NE.1).OR.(IASMBL(OIPRECA).EQ.4)) THEN
          IF (IASMBL(OIUPW).EQ.1) THEN
            CALL GUPWDX (DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
     *                  DWORK(KD),DWORK(KD+NEQU),
     *                  DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *                  TRIAS(1,NLMAX),KWORK(KVERT),KWORK(KMID),
     *                  DWORK(KCORVG),1,
     *                  UPSAM,RE,THWEIG,
     *                  IASMBL(OIALE),DWORK(KGRVEL))
          ELSE
            IF (IELT.EQ.0) THEN 
C              CALL SUPGPX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
C     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
C     *                  DWORK(KD),DWORK(KD+NEQU),
C     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
C     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
C     *                  TRIAS(1,NLMAX),E031,IASMBL(OICUBN),
C     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS), 
C     *                  IASMBL(OISTOK), DASMBL(OUPSAM), RE,  
C     *                  1,1D0,THWEIG)
              CALL SUPAPX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
     *                  DWORK(KD),DWORK(KD+NEQU),
     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
     *                  TRIAS(1,NLMAX),E031,IASMBL(OICUBN),
     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS), 
     *                  IASMBL(OISTOK), DASMBL(OUPSAM), RE,  
     *                  1,1D0,THWEIG,
     *                  IASMBL(OIALE),DWORK(KGRVEL))
            END IF

            IF (IELT.EQ.1) THEN
C              CALL SUPGPX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
C     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
C     *                  DWORK(KD),DWORK(KD+NEQU),
C     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
C     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
C     *                  TRIAS(1,NLMAX),E030,IASMBL(OICUBN),
C     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS), 
C     *                  IASMBL(OISTOK), DASMBL(OUPSAM), RE,  
C     *                  1,1D0,THWEIG)
              CALL SUPAPX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
     *                  DWORK(KD),DWORK(KD+NEQU),
     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
     *                  TRIAS(1,NLMAX),E030,IASMBL(OICUBN),
     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS), 
     *                  IASMBL(OISTOK), DASMBL(OUPSAM), RE,  
     *                  1,1D0,THWEIG,
     *                  IASMBL(OIALE),DWORK(KGRVEL))
            END IF

            IF (IELT.EQ.2) THEN
C              CALL SUPGNX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
C     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
C     *                  DWORK(KD),DWORK(KD+NEQU),
C     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
C     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
C     *                  TRIAS(1,NLMAX),E031,IASMBL(OICUBN),
C     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS), 
C     *                  IASMBL(OISTOK), DASMBL(OUPSAM), RE,  
C     *                  1,1D0,THWEIG)
              CALL SUPANX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
     *                  DWORK(KD),DWORK(KD+NEQU),
     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
     *                  TRIAS(1,NLMAX),EM31,IASMBL(OICUBN),
     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS), 
     *                  IASMBL(OISTOK), DASMBL(OUPSAM), RE,  
     *                  1,1D0,THWEIG,
     *                  IASMBL(OIALE),DWORK(KGRVEL))
            END IF

            IF (IELT.EQ.3) THEN
C              CALL SUPGNX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
C     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
C     *                  DWORK(KD),DWORK(KD+NEQU),
C     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
C     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
C     *                  TRIAS(1,NLMAX),EM30,IASMBL(OICUBN),
C     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS), 
C     *                  IASMBL(OISTOK), DASMBL(OUPSAM), RE,  
C     *                  1,1D0,THWEIG)
              CALL SUPANX(NEQU,DUP(1),DUP(1+NEQU),DUP(1),DUP(1+NEQU),
     *                  1D0,0D0,DUP(1),DUP(1+NEQU),
     *                  DWORK(KD),DWORK(KD+NEQU),
     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),
     *                  TRIAS(1,NLMAX),EM30,IASMBL(OICUBN),
     *                  (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS), 
     *                  IASMBL(OISTOK), DASMBL(OUPSAM), RE,  
     *                  1,1D0,THWEIG,
     *                  IASMBL(OIALE),DWORK(KGRVEL))
            END IF
            
          ENDIF
          
        ENDIF
        
        CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTUPW,1)

C       Incorporate Dirichlet conditions into the matrix and the
C       defect vector.

        CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTBDR,0)
     
        CALL BDRYA (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *              KWORK(KSCNPR),NBDMT)
        CALL BDRY0 (DWORK(KD),DWORK(KD+NEQU),KWORK(KSCNPR),NBDMT)
     
        CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTBDR,1)

C=======================================================================
C *** Calculation of defects and norms
C=======================================================================

C       Use RESDFK to compute the norms RESU and RESDIV from
C       the solution, the RHS and the residuum:
C
C                     || (D1,D2) ||_E 
C       RESU = -----------------------------
C              max ( ||F1||_E , ||F2||_E )
C
C                   || P ||_E
C       RESDIV = ----------------
C                || (U1,U2) ||_E
C
C       RES = || (RESU,RESDIV) ||_E


        CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTLC,0)
        
        CALL RESDFK(DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *              DWORK(KD),DWORK(KD+NEQU),DWORK(KD+2*NEQU),
     *              DRHS(1),DRHS(1+NEQU),
     *              DRHS(1+2*NEQU),NEQU,NEQP,RESU,RESDIV)

        RES = SQRT(RESU*RESU+RESDIV*RESDIV)

C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================

        IF ( RES/RESOLD.GT.1D2)
     *    IPARAM(OSTATUS) = IOR(IPARAM(OSTATUS),1)

        RESOLD=RES
        RHO   =(RES/RES0)**(1D0/DBLE(INL))
        
        CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTLC,1)
        
C=======================================================================
C *** Control of terminating the nonlinear iteration
C=======================================================================

C       Check if the nonlinear iteration can prematurely terminate.
C       In that case sei BNLENT=true - a later IF-case will leave
C       the iteration in this case.

        IF((DELU.LE.DPARAM(OEPSUR)).AND.(DELP.LE.DPARAM(OEPSPR))   .AND.
     *     (RESU.LE.DPARAM(OEPSD)) .AND.(RESDIV.LE.DPARAM(OEPSDIV)).AND.
     *     (RES.LE.EPSRES)         .AND.(INL.GE.IPARAM(ONITMIN))) 
     *     BNLEND=.TRUE.

        IF (MSHOW.GE.2) 
     *    WRITE(MTERM,1003) INL,DELU,DELP,RESU,RESDIV,RES,RHO,
     *                      OMEGA,RHOLI
        IF (MSHOW.GE.1) 
     *    WRITE(MFILE,1003) INL,DELU,DELP,RESU,RESDIV,RES,RHO,
     *                      OMEGA,RHOLI
C
C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================

C       Cancel the iteration here if one of the solvers above
C       broke down.

        IF (IPARAM(OSTATUS).NE.0) GOTO 221

C=======================================================================
C *** Autosave
C Save the solution to disc every IAUSAV iterations
C=======================================================================

        IF (IPARAM(OIAUSAV).NE.0) THEN
          IF (MOD(INL,IPARAM(OIAUSAV)).EQ.0) THEN      
            CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTPOST,0)
            
            CALL STPUT (IPARAM(OLFLAUT), CFN)
            CALL PPWRVC (0,NUVP,DUP,INL,0,CFN)
            
            CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTPOST,1)
          ENDIF
        ENDIF

C=======================================================================
C *** Return if BNLEND=true
C=======================================================================

C       Leave the nonlinear iteration if 
C       - BNLEND=true (solution converged properly or
C       - OMEGA is too small

        IF (BNLEND.OR.(ABS(OMEGA).LT.1D-1)) GOTO 221

C       That's it with the nonlinear iteration

      END DO ! INL

C=======================================================================
C *** End of the nonlinear loop
C=======================================================================

221   CONTINUE

C     In case that the maximum number of nonlinear steps were done,
C     INL = INLMAX+1 because of the DO-loop - reset INL=INLMAX
C     for output purposes.

      IF (INL.GT.IPARAM(ONITMAX)) THEN
        IPARAM(OITE) = IPARAM(ONITMAX)
      ELSE
        IPARAM(OITE) = INL
      END IF

C     A little bit of nice output :)

      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)

C     Call the linear solver routine a final time to do cleanup 
C     if necessary - this actually does not solve anything

      CALL NSLINS (NLMIN,NLMAX,TRIAS,
     *             MATDAT,VECDAT,
     *             IPARAM,DPARAM,
     *             IMGPAR,DMGPAR,
     *             IASMBL,DASMBL,
     *             -1,NITER,STATUS,RHOLI,TIMGLS)
      
C     Add the timing information of our local timing structure
C     to the timing information in the parameter block:

      DPARAM(OTNLTIM-1+OTTMG  ) = DPARAM(OTNLTIM-1+OTTMG  ) + TIMGLS(2)
      DPARAM(OTNLTIM-1+OTTPROL) = DPARAM(OTNLTIM-1+OTTPROL) + TIMGLS(3) 
      DPARAM(OTNLTIM-1+OTTREST) = DPARAM(OTNLTIM-1+OTTREST) + TIMGLS(4) 
      DPARAM(OTNLTIM-1+OTTDEF ) = DPARAM(OTNLTIM-1+OTTDEF ) + TIMGLS(5) 
      DPARAM(OTNLTIM-1+OTTSMTH) = DPARAM(OTNLTIM-1+OTTSMTH) + TIMGLS(6) 
      DPARAM(OTNLTIM-1+OTTCGSL) = DPARAM(OTNLTIM-1+OTTCGSL) + TIMGLS(7) 
      DPARAM(OTNLTIM-1+OTTBC  ) = DPARAM(OTNLTIM-1+OTTBC  ) + TIMGLS(8) 
      DPARAM(OTNLTIM-1+OTTCGC ) = DPARAM(OTNLTIM-1+OTTCGC ) + TIMGLS(9) 
      
C     Finally stop the time of this algorithm, write the result to TTNL
C     in the timing structure.

      CALL GTMAUX (DTIMIN,DPARAM,OTNLTIM-1+OTTNL,1)
      
C     TTNL represents the complete time of the algorithm. The solver
C     structure on the other hand has the general field TMTOT
C     representing the total time. So copy TTNL to TMTOT for
C     compatibility to the structure!

      DPARAM(OTMTOT) = DPARAM(OTNLTIM-1+OTTNL)

C     Print out some statistical information, finish.

      IF (MSHOW.GE.4) THEN
      
        TTMG = TIMGLS(2)
        IF (TTMG.EQ.0D0) TTMG=1D0
      
        WRITE(MTERM,*) ' Iterations nonlinear / linear :',
     *                 IPARAM(OITE),'/',IPARAM(ONLIN)
        WRITE(MTERM,*) ' Time linear solver            :',
     *                 TIMGLS(1)
        WRITE(MTERM,*) ' Time multigrid                :',
     *                 TIMGLS(2)
        WRITE(MTERM,*) ' smoothing     :',  
     *                 1.D2*TIMGLS(6)/TTMG
        WRITE(MTERM,*) ' solver        :',  
     *                 1.D2*TIMGLS(9)/TTMG
        WRITE(MTERM,*) ' defect calc.  :',  
     *                 1.D2*TIMGLS(5)/TTMG
        WRITE(MTERM,*) ' prolongation  :',  
     *                 1.D2*TIMGLS(3)/TTMG
        WRITE(MTERM,*) ' restriction   :',  
     *                 1.D2*TIMGLS(4)/TTMG
        WRITE(MTERM,1)
      ENDIF

      IF (MSHOW.GE.3) THEN            

        TTMG = TIMGLS(2)
        IF (TTMG.EQ.0D0) TTMG=1D0

        WRITE(MFILE,*) ' Iterations nonlinear / linear :',
     *                 IPARAM(OITE),'/',IPARAM(ONLIN)
        WRITE(MFILE,*) ' Time linear solver:           :',
     *                 TIMGLS(1)
        WRITE(MFILE,*) ' Time multigrid:               :',
     *                 TIMGLS(2)
        WRITE(MFILE,*) ' smoothing     :',  
     *                 1.D2*TIMGLS(6)/TTMG
        WRITE(MFILE,*) ' solver        :',  
     *                 1.D2*TIMGLS(9)/TTMG
        WRITE(MFILE,*) ' defect calc.  :',  
     *                 1.D2*TIMGLS(5)/TTMG
        WRITE(MFILE,*) ' prolongation  :',  
     *                 1.D2*TIMGLS(3)/TTMG
        WRITE(MFILE,*) ' restriction   :',  
     *                 1.D2*TIMGLS(4)/TTMG
        WRITE(MFILE,1)
      ENDIF
      
   1  FORMAT(79('-'))
1001  FORMAT(' IT RELU',5X,'RELP',5X,'DEF-U',4X,'DEF-DIV',2X,
     *       'DEF-TOT',2X,'RHONL ',3X,'OMEGNL',3X,'RHOMG')
1002  FORMAT(I3,18X,3(D9.2))
1003  FORMAT(I3,8(D9.2))

99999 END
