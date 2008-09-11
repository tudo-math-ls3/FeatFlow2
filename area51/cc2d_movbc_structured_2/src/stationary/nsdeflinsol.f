************************************************************************
* This file contains routines for solving the linear sub-problems
* that arise in the nonlinear iteration loop.
*
* The routines here capsule the complete linear solver. If the solver
* has to be exchanged by another one, it only has to be done here
* (probably with a redefinition of the Yxxx-routines in NSDEFMGROUT.F
*  and with a change of the initialization routines - but that's all!)
************************************************************************

************************************************************************
* Solve linear problem for NSDEF
*
* This routine solves a linear system Ax=b, which is prepared by
* NSDEF. The problem is solved on level NLMAX. VECDAT[NLMAX].RHS
* defines the right hand side of the system. The routine will return
* the solution of the system in VECDAT[NLMAX].SOL.
*
* The routine is called multiple times in the nonlinear loop. The
* variable ISOL counts how often the routine is called, thus it's
* possible to do e.g. some preprocessing for the whole solver in 
* the first call. When the nonlinear loop is finished, this routine
* is again called with ISOL=-1 to allow some postprocessing and/or
* cleanup of the solution process.
*
* In:
*   NLMIN  : minimum level in MATDAT/VECDAT that is filled with data
*   NLMAX  : maximum level in MATDAT/VECDAT that is filled with data
*   TRIAS  : array [1..SZTRIA,1..NNLEV] of integer
*            Triangulation structures for the underlying meshes
*            on all levels. TRIAS(.,ILEV) represents the triangulation
*            structure on level ILEV.
*   MATDAT : array {1..SZN2MI,1..NNLEV] of integer
*            TNS2DMatrixParams-structures for level NLMIN..NLMAX. 
*            The matrices on all levels are initialized with data
*            by NSDEF.
*   VECDAT : array {1..SZN2VI,1..NNLEV] of integer
*            TNS2DVectorParams-structures for level NLMIN..NLMAX. 
*            VECDAT[OLRHS,NLMAX] contains the RHS on the finest
*            level of the system.
*   ISTPAR : array [1..SZNSDI] of integer
*   DSTPAR : array [1..SZNSDD] of double
*            Integer/Double precision parameter block of NSDEF.
*   IMGPAR : array [1..SZ020I+3*SZSLVI] of integer
*   DMGPAR : array [1..SZ020D+3*SZSLVD] of double
*            Integer and double parameter blocks for the multigrid 
*            sub-solver M020. 
*            The problem-specific variables in this structure will 
*            be initialized,
*                        KOFFx, KNEQ, NLMIN, NLMAX
*            and the routine is allowed to change the solver structure
*            to its needs for solving the system.
*            The structures must be reset to original values when 
*            NSDEF is finished (indicated by a call with ISOL=-1).
*
*            The general solver parameters, e.g. smoother, coarse 
*            grid solver,... as well as the number of smoothing steps 
*            on each level must be initialized by the caller, e.g. 
*            with INMGST.
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
*            discretization. 
*
*   ISOL   - Identifier that counts how often this routine is called.
*            =1,2,3,... : This is the ISOL'th call to this routine
*            =-1        : Final call to this routine for postprocessing.
*                         No linear system is actually to solve.
*   VECDAT[NLMAX].RHS   : Right hand side of the system on finest level
*
* Out: (only if ISOL>=0 !)
*   VECDAT[NLMAX].SOL - solution of the system
*   NITER             - number of iterations of the solver
*   STATUS            - result of the solver
*                       =0: solution calculated successfully
*                       >0: error
*   RHO               - convergence rate of the solver
*   TIMING            - array [1..9] of double
*                       Timing information. Every entry in TIMING
*                       represents a special part of the solver,
*                       which time is measured by NSDEF independently.
*                       This routine has to add the time for the
*                       different components to the entries. The entries
*                       are defined as follows:
*                        TIMING(1) - Total time
*                        TIMING(2) - Total time for MG solver
*                        TIMING(3) - Time for prolongation
*                        TIMING(4) - Time for restriction
*                        TIMING(5) - Time for defect calculation
*                        TIMING(6) - Time for smoothing operations
*                        TIMING(7) - Time for coarse-grid solver
*                        TIMING(8) - Time for boundary conditions
*                        TIMING(9) - Time for coarse-grid corrections
*
* Remarks:
* a) The VECDAT structure is used as auxiliary structure during the
*  solution process. The content of the structure is undefined
*  upon leaving this routine (except for the solution vector in
*  VECDAT[NLMAX].SOL and the not-touched RHS-vector in
*  VECDAT[NLMAX].RHS).
* b) During the solution process, callback routines of the solver might
*  access information in the Status-part of the nonlinear solver.
*  NSDEF has prepared some information about the problem there to
*  allow the callback routines to fulfill their duty (e.g. all 
*  triangulation/matrix information can be found there as well as
*  whether there are Neumann boundaries in the problem,...)
* c) The initial solution vector VECDAT[NLMAX].SOL is filled with 0
*  by the caller. For all Dirichlet components, the entries in
*  this vector must not be changed!
************************************************************************

      SUBROUTINE NSLINS (NLMIN,NLMAX,TRIAS,
     *                   MATDAT,VECDAT,
     *                   ISTPAR,DSTPAR,
     *                   IMGPAR,DMGPAR,
     *                   IASMBL,DASMBL,
     *                   ISOL,NITER,STATUS,RHO,TIMING)  

      IMPLICIT NONE

      INCLUDE 'cout.inc'
      INCLUDE 'cmem.inc'

      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'ssolvers.inc'
      
      INCLUDE 'sassembly.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'
      
      INCLUDE 'm020.inc'

      INCLUDE 'snsdeflinsol.inc'

C     parameters
      
      DOUBLE PRECISION DSTPAR (*),DMGPAR(*),DASMBL(SZASMD)
      INTEGER TRIAS(SZTRIA,NNLEV)
      INTEGER MATDAT(SZN2MI,NNLEV),VECDAT(SZN2VI,NNLEV),NLMIN,NLMAX
      INTEGER ISTPAR (*),IMGPAR(*),IASMBL(SZASMI)
      INTEGER ISOL,STATUS,NITER
      DOUBLE PRECISION RHO,TIMING(9)
      
C     local variables

      INTEGER I1,NSLISL(SZNLSI),NEQ,LTMP,LTMSOL,LTMRHS
      DOUBLE PRECISION DUMMY
      
C     local variables for UMFPACK solver

      DOUBLE PRECISION CNTRU4(20),INFOU4(90)
      INTEGER KU4SYM,KU4NUM,KSYS
      DOUBLE PRECISION T

C     Backup for MG structures:

      INTEGER IMGBCK (SZ020I+3*SZSLVI)
      DOUBLE PRECISION DMGBCK (SZ020D+3*SZSLVD)
      
C     Temporary arrays for vector handles for the preparation of M020

      INTEGER LOFFX(NNLEV),LOFFD(NNLEV),LOFFB(NNLEV)
      
C     Multigrid components

      EXTERNAL YAX2,YPROL2,YREST2,YSM2,YEX2,YDBC2,YSTEP2,I000,YMG0C
      
C     Leave the routine when called with ISOL=-1.
C     There is nothing to solve and nothing to clean up in that case.

      IF (ISOL.LT.0) RETURN
      
      IF (IMGPAR(OSLTAG).EQ.13) THEN
      
C       Use a direct UMFPACK solver.

        CALL GTMAUX (T,TIMING,1,0)

C       Factorize the matrix on the maximum level:

        CALL PRUMF4 (MATDAT(1,NLMAX),TRIAS(1,NLMAX),3,
     *               CNTRU4,KU4SYM,KU4NUM)

C       Solve the system. Note that UMFPACK expects the matrix in
C       CSR format, which is transposed to our matrix format 7 --
C       So we solve the transposed system:
        
        KSYS = 1
        CALL UMF4SOL(KSYS,DWORK(L(VECDAT(OLSOL,NLMAX))),
     *                    DWORK(L(VECDAT(OLRHS,NLMAX))),
     *               KU4NUM,CNTRU4,INFOU4)
  
C       Release UMFPACK4 objects from heap again:
  
        CALL UMF4FNUM(KU4NUM)
        CALL UMF4FSYM(KU4SYM)
  
C       Initialize return status:
  
        STATUS = 0
        NITER  = 1
        RHO    = 0.0D0
        
        CALL GTMAUX (T,TIMING,1,1)
  
      ELSE IF (IMGPAR(OSLTAG).EQ.6) THEN
      
C       Use the BiCGStab-solver with Multigrod preconditioning.
C      
C       Make a backup of the multigrid structures, so that we can
C       restore it onto return. Notice that the structure has the size
C       of two standard-structure plus a M020 structure
C       (-> BiCGStab-solver, Coarse grid solver, M020 precond.).

        CALL LCP3(IMGPAR,IMGBCK,SZ020I+3*SZSLVI)
        CALL LCP1(DMGPAR,DMGBCK,SZ020D+3*SZSLVD)
        
C       The multigrid structure we have to initialize is to be found
C       at position SZSLVx of the IMGPAR/DMGPAR array since it's
C       the second structure in the block! We have to be aware
C       of that when initializing the parameters.
C
C       Initialization of the offset arrays KOFFX,KOFFB,KOFFD and KNEQ  
C       for the multigrid solver.
C       The offset for the solution and RHS on the finest level will
C       be modified later, as it's used in the preconditioning and
C       must not point to the real solution/RHS.

        DO I1=NLMIN,NLMAX
          LOFFX(I1)=VECDAT(OLSOL,I1)
          LOFFB(I1)=VECDAT(OLRHS,I1)
          LOFFD(I1)=VECDAT(OLTMP,I1)
          IMGPAR(SZSLVI+OKNEQ+I1-1) = VECDAT(ONEQV,I1)
        END DO
        
C       Allocate memory for two additional arrays on the finest 
C       level for the MG structure.

        NEQ = VECDAT(ONEQV,NLMAX)

        CALL ZNEW (NEQ,-1,LTMSOL,'DTMSOL')
        CALL ZNEW (NEQ,-1,LTMRHS,'DTMRHS')

C       And modify the MG structure on the finest level with 
C       the handles. MG will use these vectors for preconditioning,
C       while the real solution/RHS vectors (which are given in
C       the VECDAT structure on the finest level)
C       are handled by BiCGStab.

        LOFFX(NLMAX) = LTMSOL
        LOFFB(NLMAX) = LTMRHS
        
C       Inform the MG solver on which levels we want to solve:

        IMGPAR(SZSLVI+ONLMIN) = NLMIN
        IMGPAR(SZSLVI+ONLMAX) = NLMAX
        
C       Prepare the user defined parameter structure to match the
C       problem specific data.
C
C       Make a copy of the MATDAT structure into parameter block. 
C       The linear solver callback routines can use that to access the 
C       system matrices.
        
        CALL LCP3(MATDAT,NSLISL(OIMATS),NNLEV*SZN2MI)
        
C       Also copy information about our triangulations into the
C       parameter block, as matrix-vector multiplication routines 
C       need that.
        
        CALL LCP3(TRIAS,NSLISL(OITRIAS),NNLEV*SZTRIA)
        
C       Inform our callback routines for the linear solver whether 
C       there are Neumann boundary components to handle; they'll
C       access our solver structure.

        NSLISL(OICNEUM) = IASMBL(OINEUM)
        
C       What's up with the coarse grid solver? Do we have to use
C       UMFPACK4?
C       The coarse grid solver structure is at position 3 in the list:

        IF (IMGPAR(SZSLVI+SZ020I+OSLTAG).EQ.13) THEN
        
C         Ok, we should use UMFPACK4. We symbolically and numerically 
C         factorize theglobal matrix on the coarse-grid solver level. 
C         We create the CONTROL-structure of UMFPACK in the user-defined 
C         block of the solver structure of the coarse grid solver. This
C         allows us later to use it in the callback routines.
C         The handle of the symbolically factorized matrix is also
C         stored into the user-defined area of the coarse grid solver:

          CALL PRUMF4 (MATDAT(1,NLMIN),TRIAS(1,NLMIN),3,
     *                 DMGPAR(SZSLVI+SZ020D+SZSLVD+ODUDSLA),
     *                 KU4SYM,KU4NUM)
          
          IMGPAR(SZSLVI+SZ020I+OIUDSLA)   = KU4SYM
          IMGPAR(SZSLVI+SZ020I+OIUDSLA+1) = KU4NUM
        
        END IF
        
C       Call the preparation routine to initialize the rest of the
C       multigrid structure so we can call M020.
C       Don't initialize the pre-/postsmoothing steps; this had to
C       be done by the caller.
        
        CALL PRM020 (IMGPAR(1+SZSLVI), DMGPAR(1+SZSLVD), NLMIN,NLMAX,
     *               -1, -1, LOFFX, LOFFB, LOFFD)
        
C       Now M020 is ready for action...
C
C       The starting vector of the iteration can be filled with 0,
C       but we can leave out this step. By definition, the starting
C       vector is already filled with 0 by the caller!
C       CALL LCL1 (DWORK(L(LOFFX(NLMAX))),VECDAT(ONEQV,NLMAX))

C       Allocate temporary space for BiCGStab solver.

        CALL ZNEW (5*NEQ,-1,LTMP,'BCGSTM')
        
C       Clear the first 16 entries in the user defined area
C       of the double-prec. data block of the BiCGStab solver.
C       The timing information of the MG solver will be
C       created there.

        CALL LCL1(DMGPAR(ODUDSLA+1),16)

C       Call BiCGStab to solve the system. We don't have a user
C       defined double precision parameter block, so we use a dummy
C       variable for that.

        DUMMY = 0
        
        CALL II01X(IMGPAR,DMGPAR, NSLISL,DUMMY, NEQ,
     *             DWORK(L(VECDAT(OLSOL,NLMAX))),
     *             DWORK(L(VECDAT(OLRHS,NLMAX))),
     *             DWORK(L(LTMP)),DWORK(L(LTMP)+NEQ),
     *             DWORK(L(LTMP)+2*NEQ),DWORK(L(LTMP)+3*NEQ),
     *             DWORK(L(LTMP)+4*NEQ),
     *             YAX2,YMG0C,I000)
        
C       Release auxiliary memory for BiCGStab.

        CALL ZDISP (0,LTMP,'BCGSTM')
        CALL ZDISP (0,LTMRHS,'DTMRHS')
        CALL ZDISP (0,LTMSOL,'DTMSOL')
        
C       Was UMFPACK4 used as coarse grid solver?

        IF (IMGPAR(SZSLVI+SZ020I+OSLTAG).EQ.13) THEN
        
C         Yes, it was used. Release the numeric and 
C         symbolic factorization.
        
          CALL UMF4FNUM(KU4NUM)
          CALL UMF4FSYM(KU4SYM)
        
        END IF
     
C       Return status, number of iterations and convergence rate
        
        STATUS = IMGPAR(OSTATUS)
        NITER  = IMGPAR(OITE)
        RHO    = DMGPAR(ORHO)
        
C       Add the timing information of the solver to the timing
C       information parameter block.
C       The time for the solver is to be found in 
C       DMGPAR(OTMTOT ). The overall time for MG and its components
C       is created in the user defined area of the BiCGStab solver.

        TIMING(1) = TIMING(1) + DMGPAR(OTMTOT )
        TIMING(2) = TIMING(2) + DMGPAR(ODUDSLA+1)
        TIMING(3) = TIMING(3) + DMGPAR(ODUDSLA+2)
        TIMING(4) = TIMING(4) + DMGPAR(ODUDSLA+3)
        TIMING(5) = TIMING(5) + DMGPAR(ODUDSLA+4)
        TIMING(6) = TIMING(6) + DMGPAR(ODUDSLA+5)
        TIMING(7) = TIMING(7) + DMGPAR(ODUDSLA+6)
        TIMING(8) = TIMING(8) + DMGPAR(ODUDSLA+7)
        TIMING(9) = TIMING(9) + DMGPAR(ODUDSLA+8)
                                
C       Restore original multigrid configuration

        CALL LCP3(IMGBCK,IMGPAR,SZ020I+3*SZSLVI)
        CALL LCP1(DMGBCK,DMGPAR,SZ020D+3*SZSLVD)
        
      ELSE IF (IMGPAR(OSLTAG).EQ.11) THEN
      
C       Use the multigrid solver.
C      
C       Make a backup of the multigrid structures, so that we can
C       restore it onto return.

        CALL LCP3(IMGPAR,IMGBCK,SZ020I+3*SZSLVI)
        CALL LCP1(DMGPAR,DMGBCK,SZ020D+3*SZSLVD)

C       Initialization of the offset arrays KOFFX,KOFFB,KOFFD and KNEQ  
C       for the multigrid solver

        DO I1=NLMIN,NLMAX
          LOFFX(I1)=VECDAT(OLSOL,I1)
          LOFFB(I1)=VECDAT(OLRHS,I1)
          LOFFD(I1)=VECDAT(OLTMP,I1)
          IMGPAR(OKNEQ+I1-1) = VECDAT(ONEQV,I1)
        END DO
        
C       Inform the MG solver on which levels we want to solve:

        IMGPAR(ONLMIN) = NLMIN
        IMGPAR(ONLMAX) = NLMAX
        
C       Prepare the user defined parameter structure to match the
C       problem specific data.
C
C       Make a copy of the MATDAT structure into parameter block. 
C       The linear solver callback routines can use that to access the 
C       system matrices.
        
        CALL LCP3(MATDAT,NSLISL(OIMATS),NNLEV*SZN2MI)
        
C       Also copy information about our triangulations into the
C       parameter block, as matrix-vector multiplication routines 
C       need that.
        
        CALL LCP3(TRIAS,NSLISL(OITRIAS),NNLEV*SZTRIA)
        
C       Inform our callback routines for the linear solver whether 
C       there are Neumann boundary components to handle; they'll
C       access our solver structure.

        NSLISL(OICNEUM) = IASMBL(OINEUM)
        
C       What's up with the coarse grid solver? Do we have to use
C       UMFPACK4?
C       The coarse grid solver structure is at position 2 in the list:

        IF (IMGPAR(SZ020I+OSLTAG).EQ.13) THEN
        
C         Ok, we should use UMFPACK4. We symbolically and numerically 
C         factorize theglobal matrix on the coarse-grid solver level. 
C         We create the CONTROL-structure of UMFPACK in the user-defined 
C         block of the solver structure of the coarse grid solver. This
C         allows us later to use it in the callback routines.
C         The handle of the symbolically factorized matrix is also
C         stored into the user-defined area of the coarse grid solver:

          CALL PRUMF4 (MATDAT(1,NLMIN),TRIAS(1,NLMIN),3,
     *                 DMGPAR(SZ020D+SZSLVD+ODUDSLA),
     *                 KU4SYM,KU4NUM)
          
          IMGPAR(SZ020I+OIUDSLA)   = KU4SYM
          IMGPAR(SZ020I+OIUDSLA+1) = KU4NUM
        
        END IF
        
C       Call the preparation routine to initialize the rest of the
C       multigrid structure so we can call M020.
C       Don't initialize the pre-/postsmoothing steps; this had to
C       be done by the caller.
        
        CALL PRM020 (IMGPAR, DMGPAR, NLMIN,NLMAX,
     *               -1, -1, LOFFX, LOFFB, LOFFD)
        
C       Now M020 is ready for action...
C
C       To make the starting vector of the MG iteration not represent
C       a totally random vector, we fill it with 0:

        CALL LCL1 (DWORK(L(LOFFX(NLMAX))),VECDAT(ONEQV,NLMAX))

C       Call multigrid to solve the system. We don't have a user
C       defined double precision parameter block, so we use a dummy
C       variable for that.

        DUMMY = 0
        
        CALL M020 (IMGPAR,DMGPAR, NSLISL,DUMMY, 
     *             DWORK(1),DWORK(1),DWORK(1),
     *             YAX2,YPROL2,YREST2,YSM2,YSM2,
     *             YEX2,YEX2,YDBC2,YSTEP2,I000)
     
C       Was UMFPACK4 used as coarse grid solver?

        IF (IMGPAR(SZ020I+OSLTAG).EQ.13) THEN
        
C         Yes, it was used. Release the numeric and 
C         symbolic factorization.
        
          CALL UMF4FNUM(KU4NUM)
          CALL UMF4FSYM(KU4SYM)
        
        END IF
     
C       Return status, number of iterations and convergence rate
        
        STATUS = IMGPAR(OSTATUS)
        NITER  = IMGPAR(OITE)
        RHO    = DMGPAR(ORHO)
        
C       Add the timing information of the solver to the timing
C       information parameter block.
C       As we directly use MG as solver, TIMING(1) and TIMING(2)
C       are identical.

        TIMING(1) = TIMING(1) + DMGPAR(OTMTOT )
        TIMING(2) = TIMING(2) + DMGPAR(OTMTOT )
        TIMING(3) = TIMING(3) + DMGPAR(OTMPROL)
        TIMING(4) = TIMING(4) + DMGPAR(OTMREST)
        TIMING(5) = TIMING(5) + DMGPAR(OTMDEF )
        TIMING(6) = TIMING(6) + DMGPAR(OTMSMTH)
        TIMING(7) = TIMING(7) + DMGPAR(OTMCGSL)
        TIMING(8) = TIMING(8) + DMGPAR(OTMBC  )
        TIMING(9) = TIMING(9) + DMGPAR(OTMCGC )
        
C       Restore original multigrid configuration

        CALL LCP3(IMGBCK,IMGPAR,SZ020I+3*SZSLVI)
        CALL LCP1(DMGBCK,DMGPAR,SZ020D+3*SZSLVD)
        
      ELSE
        
        WRITE (MTERM,'(A,I6)') 'FATAL: Unknown linear solver:',
     *                         IMGPAR(OSLTAG)
        STOP
        
      END IF

      END
      
************************************************************************
* Prepare UMFPACK4 solver
*
* This routine prepares the call to UMFPACK as linear solver - for
* the main problem, as coarse grid solver or what else...
* It accepts the matrix structure of one level and calculates the
* symbolical factorization and - if desired - the numerical one, too.
*
* In:
*   IMAT    : TNS2DMatrixParams structure for the matrices on the
*             current level
*   TRIA    : array [1..SZTRIA] of integer
*             Triangulation structure of the mesh corresponding to IMAT
*   IMETH   : decides on what to do:
*             =1: Build the main matrix and factorize symbolically.
*             =2: Build the main matrix and factorize numerically.
*                 KU4SYM must be a handle to a previously symbolically
*                 factorized matrix with the same structure!
*             =3: Build the main matrix and factorize both, 
*                 symbolically and numerically.
*   KU4SYM  : If IMETH=2, KU4SYM must be a handle to a previously 
*             symbolically factorized matrix. Otherwise, KU4SYM can
*             be undefined on call.
*   CNTRU4  : array [1..20] of double
*             If IMETH=2, existing CONTROL structure for UMFPACK.
*             Otherwise CNTRU4 can be undefined on call.
*   INFOU4  : array [1..90] of double
*             If IMETH=2, existing INFO structure for UMFPACK
*             Otherwise INFOU4 can be undefined on call.
* 
* Out:
*   CNTRU4  : array [1..20] of double
*             CONTROL structure for UMFPACK
*   KU4SYM  : Handle of the symbolically factorized matrix;
*             Unchanged, if only numerical factorization should
*             be performed
*   KU4NUM  : Handle of the numerically factorized matrix;
*             =0, if no numerical factorization should be performed
************************************************************************

      SUBROUTINE PRUMF4 (IMAT,TRIA,IMETH,
     *                   CNTRU4,KU4SYM,KU4NUM)
     
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'smat2dns.inc'
      INCLUDE 'stria.inc'
     
C     parameters
     
      DOUBLE PRECISION CNTRU4(20)
      INTEGER KU4SYM,KU4NUM,IMAT(*),IMETH
      INTEGER TRIA(SZTRIA)
      
C     local variables

      DOUBLE PRECISION INFOU4(90)
      INTEGER LA,LCOL,LLD,NEQ,LTMP,LB1,LB2
      INTEGER KLA(3,3),KCOLA(3,3),KLDA(3,3),KNEQ(3,3)
      INTEGER LLA(3,3),LCOLA(3,3),LLDA(3,3)
      
C     Prepare factorization in general.
C     First clear the structure arrays that define the global matrix:
      
      CALL LCL3(KCOLA,9)
      CALL LCL3(KLDA,9)
      CALL LCL3(KLA,9)
      CALL LCL3(KNEQ,9)
      
C     Initialize the index-arrays in KLxxx with the positions
C     of the sub-matrices in KWORK/DWORK.
C     Don't initialize the handle of the entries of B1/B2 as we
C     have to build modified B-matrices later!
      
      KLA(1,1) = L(IMAT(OLA1))
      KLA(2,2) = L(IMAT(OLA1))
      
      KCOLA(1,1) = L(IMAT(OLCLA1))
      KCOLA(2,2) = L(IMAT(OLCLA1))
      KCOLA(1,3) = L(IMAT(OLCLB1))
      KCOLA(2,3) = L(IMAT(OLCLB2))

      KLDA(1,1) = L(IMAT(OLLDA1))
      KLDA(2,2) = L(IMAT(OLLDA1))
      KLDA(1,3) = L(IMAT(OLLDB1))
      KLDA(2,3) = L(IMAT(OLLDB2))
      
C     What's missing is the transposed B-matrices.
C     Allocate temporary arrays for the creation
C     of B1^T and B2^T:
      
      CALL ZNEW(IMAT(ONEQA),-3,LTMP,'KTMP  ')
      CALL ZNEW(IMAT(ONEQA)+1,-3,LLDA(3,1),'KTMP  ')
      CALL ZNEW(IMAT(ONEQA)+1,-3,LLDA(3,2),'KTMP  ')
      CALL ZNEW(IMAT(ONB1),-1,LLA(3,1),'DTMP  ')
      CALL ZNEW(IMAT(ONB2),-1,LLA(3,2),'DTMP  ')
      CALL ZNEW(IMAT(ONB1),-3,LCOLA(3,1),'KTMP  ')
      CALL ZNEW(IMAT(ONB2),-3,LCOLA(3,2),'KTMP  ')
      
C     Initialize the missing KLxxx-entries with information,
C     where the transposed matrices will be saved to:
      
      KLDA(3,1) = L(LLDA(3,1))
      KLDA(3,2) = L(LLDA(3,2))
      KLA(3,1) = L(LLA(3,1))
      KLA(3,2) = L(LLA(3,2))
      KCOLA(3,1) = L(LCOLA(3,1))
      KCOLA(3,2) = L(LCOLA(3,2))
      KNEQ(1,1) = IMAT(ONEQA)
      KNEQ(2,2) = IMAT(ONEQA)
      KNEQ(1,3) = IMAT(ONEQA)
      KNEQ(2,3) = IMAT(ONEQA)      
      KNEQ(3,1) = IMAT(ONEQB)
      KNEQ(3,2) = IMAT(ONEQB)      
      
C     Transpose B1:
      
      CALL TRM79 (IMAT(ONEQA),
     *            IMAT(ONEQB),
     *            DWORK(L(IMAT(OLB1))),
     *            KWORK(L(IMAT(OLCLB1))),
     *            KWORK(L(IMAT(OLLDB1))),
     *            KWORK(L(LTMP)),
     *            DWORK(KLA(3,1)),
     *            KWORK(KCOLA(3,1)),
     *            KWORK(KLDA(3,1)))
     
C     Transpose B2:
     
      CALL TRM79 (IMAT(ONEQA),
     *            IMAT(ONEQB),
     *            DWORK(L(IMAT(OLB2))),
     *            KWORK(L(IMAT(OLCLB2))),
     *            KWORK(L(IMAT(OLLDB2))),
     *            KWORK(L(LTMP)),
     *            DWORK(KLA(3,2)),
     *            KWORK(KCOLA(3,2)),
     *            KWORK(KLDA(3,2)))
     
C     Build modified B-matrices.
C     Take the original B1/B2 matrices and duplicate them:
     
      LB1 = 0
      LB2 = 0
      CALL ZCPY(IMAT(OLB1),'B1    ',LB1,'B1    ')
      CALL ZCPY(IMAT(OLB2),'B2    ',LB2,'B2    ')
      
C     Initialize the handles
      
      KLA(1,3) = L(LB1)
      KLA(2,3) = L(LB2)
      
C     Assemble all matrices together to a global one.
C     This will give LA, LCOL, LLD and NEQ:
      
      CALL MTXA79 (3,3,DWORK,KWORK,KWORK,
     *             KLA,KCOLA,KLDA,KNEQ,
     *             .TRUE.,LA,LCOL,LLD,NEQ)
     
C     Modify global matrix. Implement Dirichlet boundary conditions.
C     Overwrite all lines of the velocity components that correspond
C     to Dirichlet boundary segments by unit vectors - the velocity
C     A-matrices as well as the B-part in each line:

      CALL BDRYGA (IMAT(ONEQA),
     *             DWORK(L(LA)),KWORK(L(LCOL)),KWORK(L(LLD)),
     *             KWORK(L(TRIA(OTRIUD+1))),TRIA(OTRIUD))

C     We can release all the temporary matrices now:

      CALL ZDISP(0,LB2,'B2    ')
      CALL ZDISP(0,LB1,'B1    ')
      CALL ZDISP(0,LTMP,'KTMP  ')
      CALL ZDISP(0,LLDA(3,1),'KTMP  ')
      CALL ZDISP(0,LLDA(3,2),'KTMP  ')
      CALL ZDISP(0,LLA(3,1),'DTMP  ')
      CALL ZDISP(0,LLA(3,2),'DTMP  ')
      CALL ZDISP(0,LCOLA(3,1),'KTMP  ')
      CALL ZDISP(0,LCOLA(3,2),'KTMP  ')
      
C     For the factorization, UMFPACK needs all arrays 0-based;
C     so convert the KCOL/KLD array from our matrix from
C     1-based to 0-based:
      
      CALL M7IDSH (KWORK(L(LCOL)),KWORK(L(LLD)),NEQ)

      IF (IMETH.EQ.1) THEN
      
C       Perform symbolic factorization. 
C       Initialize the control structure:
        
        CALL UMF4DEF(CNTRU4)
        
C       factorize symbolically:
        
        CALL UMF4SYM(NEQ,NEQ,KWORK(L(LLD)),KWORK(L(LCOL)),
     *               DWORK(L(LA)),KU4SYM,CNTRU4,INFOU4)
     
        IF (INFOU4(1).LT.0) THEN
          WRITE (*,'(A,F7.0)') 'UMFPACK-Error: ',INFOU4(1)
          RETURN
        END IF

      END IF
        
      IF (IMETH.EQ.2) THEN

C       Compute the numeric factorization

        CALL UMF4NUM(KWORK(L(LLD)),KWORK(L(LCOL)),DWORK(L(LA)),
     *               KU4SYM,KU4NUM,CNTRU4,INFOU4)

        IF (INFOU4(1).LT.0) THEN
          WRITE (*,'(A,F7.0)') 'UMFPACK-Error: ',INFOU4(1)
          RETURN
        END IF
        
      END IF

      IF (IMETH.EQ.3) THEN
      
C       Initialize the control structure:
        
        CALL UMF4DEF(CNTRU4)
        
C       factorize symbolically:
        
        CALL UMF4SYM(NEQ,NEQ,KWORK(L(LLD)),KWORK(L(LCOL)),
     *               DWORK(L(LA)),KU4SYM,CNTRU4,INFOU4)
     
        IF (INFOU4(1).LT.0) THEN
          WRITE (*,'(A,F7.0)') 'UMFPACK-Error: ',INFOU4(1)
          RETURN
        END IF
        
C       Compute the numeric factorization

        CALL UMF4NUM(KWORK(L(LLD)),KWORK(L(LCOL)),DWORK(L(LA)),
     *               KU4SYM,KU4NUM,CNTRU4,INFOU4)

        IF (INFOU4(1).LT.0) THEN
          WRITE (*,'(A,F7.0)') 'UMFPACK-Error: ',INFOU4(1)
          RETURN
        END IF

      END IF

C     Release the global matrix again
      
      CALL ZDISP(0,LA,'DA   ')
      CALL ZDISP(0,LCOL,'KCOL ')
      CALL ZDISP(0,LLD,'KLD   ')

      END
