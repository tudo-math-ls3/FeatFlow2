************************************************************************
*   M020 - extended multigrid solver
************************************************************************
* The M020 routine (based on M011) is an extended version of the
* multigrid solver that is designed to work
* - without COMMON blocks
* - using less parameters in the call
* - with more influence of the user to the solver by parameters
* To use this solver, the user has to provide a set of callback-routines
* for prolongation, restriction, matrix multiplication, solution of
* linear systems on any level,...
*
* The solver itself does not use any COMMON blocks and makes no use
* of any memory management routines. M020 is designed as a Black-Box
* solver. The intelligence of the solver comes from the callback-
* routines, which implement the real work. The callback routines are
* the ones that know about the structure of the matrices, possible
* FEM solution vectors,...
*
* The behaviour of the solver is defined by a parameter block of
* integer and double variables. This parameter block is also passed
* by the solver to all callback-routines. The first part of the
* parameter block is reserved for the parameters of the solver.
* The caller can add more variables to the end of this parameter
* block with information for the callback routines, e.g. the handles
* to matrices, parameters for the solver components,...
*
* Apart of these both multigrid parameter blocks, a user defined
* integer and double precision data block is passed to the
* callback routines as well. These can be used by the caller to
* define the behaviour of "its" callback routines, as it 
* allows the callback routines to directly access information
* from the caller.
*
* All in all the following arguments have to be passed to M020.
* Here, SZMGRx >= SZM020x are theoretical constants, defined by the
* user when creating the structures IPARAM/DPARAM, as the user
* can define these arrays larger, containing additional information 
* for coarse grid solver(s), smoothers,...
*
* In:
*   IPARAM  - array [1..SZMGRI] of integer
*             Integer parameter structure for M020
*   DPARAM  - array [1..SZMGRI] of integer
*             Double precision parameter structure for M020
*   IDATA   - array [1..*] of integer
*             Start address of a user defined array, which is 
*             passed to all callback routines.
*   DDATA   - array [1..*] of double
*             Start address of a user defined array, which is 
*             passed to all callback routines.
*   DX      - array [1..*] of double precision
*             Starting address of solution vectors
*   DB      - array [1..*] of double precision
*             Starting address of right-hand-side vectors
*   DD      - array [1..*] of double precision
*             Starting address of temporary vectors
*   DAX     - SUBROUTINE (DX,DAX,NEQ,A1,A2,IPARAM,DPARAM,IDATA,DDATA)
*             Callback-routine, Matrix-vector multiplication.
*             Returns:   DAX := A1*A*DX+A2*DAX  
*   DPROL   - SUBROUTINE (DX,DFINE,IPARAM,DPARAM,IDATA,DDATA)
*             Callback-routine, Prolongation.
*             Calculates DFINE := Prolongation(DX) of solution on
*             corser level ILEV-1 to current (finer) level ILEV
*   DREST   - SUBROUTINE (DFINE,DX,IPARAM,DPARAM,IDATA,DDATA)
*             Callback-routine, Restriction.
*             Calculates DX := Restriction(DFINE) of solution on finer
*             level ILEV+1 to current (coarser) level ILEV
*   DPRSM   - SUBROUTINE (DX,DB,DD,NEQ,NSM,IPARAM,DPARAM,IDATA,DDATA)
*             Callback-routine, Pre-smoothing.
*             Performs NSM smoothing steps to DX. DD can be used
*             as auxiliary vector.
*   DPOSM   - SUBROUTINE (DX,DB,DD,NEQ,NSM,IPARAM,DPARAM,IDATA,DDATA)
*             Callback-routine, Post-smoothing.
*             Performs NSM smoothing steps to DX. DD can be used
*             as auxiliary vector.
*   DBC     - SUBROUTINE (DX,NEQ,IPARAM,DPARAM,IDATA,DDATA)
*             Callback-routine, Boundary-implementation
*             Implements boundary conditions into solution vector
*   DSTEP   - SUBROUTINE (DX,DD,DB,DSTEPP,IPARAM,DPARAM,IDATA,DDATA)
*             Callback-routine, Step-Size control.
*             Calculates DSTEPP = step-length parameter.
*                               = 1, if no step length control is used
*   DEX     - SUBROUTINE (DX,DB,DD,NEQ,RHO,ITE,IPARAM,DPARAM,IDATA,DDATA)
*             Callback-routine, coarse grid solver if more than
*             one level is used. The solver should solve up to
*             the accuracy the MG-solver should solve to (can be
*             taken from DPARAM), but is not forced to do that.
*             Returns: DX = A^-1 * DB
*                      ITE = number of iterations, if iterative
*                            solver is used; 1 otherwise.
*                      RHO = convergence rate
*   DEXS    - SUBROUTINE (DX,DB,DD,NEQ,RHO,ITE,IPARAM,DPARAM,IDATA,DDATA)
*             Callback-routine, Coarse grid solver if there's only one 
*             level; i.e. fallback-routine to standard one-level solver.
*             The solver should solve up to the accuracy the MG-solver 
*             should solve to (can be taken from DPARAM).
*             Returns: DX = A^-1 * DB
*                      ITE = number of iterations, if iterative
*                            solver is used; 1 otherwise.
*                      RHO = convergence rate
*   DFILT   - SUBROUTINE (DX,NEQ,IALGP,IPARAM,DPARAM,IDATA,DDATA)
*             Callback-routine, adaptive filtering.
*             Performs filtering on different positions of the
*             algorithm. Is only used if IPARAM(OFILT)<>0.
* 
* The behaviour of the algorithm is determined by the "Input-Variables"
* part in the IPARAM/DPARAM parameter blocks.
*
* The routine uses a couple of vectors for the computation: DX, DB and
* DD. The variables in the parameter list of M020 are pointers
* to a large array containing all these vectors. The starting
* indices of the vectors on the different levels are given by the
* KOFFx offset array in the integer parameters. E.g. the space for
* the solution vectors for level 1,2,3,... are expected to be at
* DX(1+KOFFX(1)), DX(1+KOFFX(2)), DX(1+KOFFX(2)),...
* Concerning the content of the vectors themselves, the caller has to
* initialise the following:
* - The initial start vector, which is expected at DX(1+KOFFX(NLMAX)).
* - The right hand side of the system, which is expected in 
*   DB(1+KOFFB(NLMAX)).
*
* The output of the algorithm is returned by:
* - the "Output-Variables" block in IPARAM/DPARAM
* - the DX-vector on the finest level, which is a new approximation
*   to the solution vector; can be found at DX(1+KOFFX(NLMAX))
*
* The algorithm uses DD (on levels NLMIN..NLMAX) and DX (on level 
* NLMIN..NLMAX-1) as temporary vectors during the computation. 
* The parameter blocks IPARAM/DPARAM contain further temporary 
* variables which are only valid during the computation, and which 
* inform the callback-routine about the current status of the 
* algorithm (current level,...). All these variables are initialized 
* internally. The user has only to specify the input parameters and 
* gets the result in the output variables and the solution vector.
*
* The multigrid algorithm stops if 
* - the iteration is divergent (both criteria introduced by
*   DIVREL and DIVABS are fulfilled)
* - the iteration has been convergent (both criteria introduced by
*   EPSREL and EPSABS are fulfilled) and the minimum number of
*   iterations is reached
* - the maximum number of iterations is reached
* If NODEF is given <> 0, the residuals are not checked (except for the
* maximum norm of the initial residuum) and there are always
* NITMAX iterations performed.
************************************************************************

***********************************************************************
* Description of multigrid filtering by DFILT:
* 
* In:
*  DX     - Solution vector; array [1..NEQ] of double
*  NEQ    - length of solution vector
*  IALGP  - Position in the algorithm
*           0=undefined - enforce the filtering
*           1=on start of the algorithm on finest level; 
*             DX is solution vector
*           2=on start of the MG sweep on the current level;
*             DX is the first solution vector before the
*             sweep starts (0D0-array normally, escept for the
*             finest level, there it's the result of the previous
*             call with IALGP=1) 
*           3=before smoothing; DX is solution vector 
*           4=after smoothing; DX is solution vector 
*           5=before restriction; DX is defect vector
*           6=after restriction; DX is defect vector
*           7=before coarse grid solver;
*             DX is RHS vector on coarsest level
*             (the result of IALGP=6 !)
*           8=before coarse grid solver;
*             DX is start vector on coarsest level (normally filled
*             with 0)
*           9=after coarse grid solver; 
*             DX is calculated solution vector on coarsest level
*          10=before prolongation;
*             DX is update-vector on coarser level
*          11=after prolongation; DX is prolongated update vector
*          12=after coarse grid correction, before post smoothing;
*             DX is solution vector
*          13=after post-smoothing; DX is solution vector
*          14=after one step of MG; DX is the solution vector on the
*             finest level.
*          15=after one step of MG; DX is the last calculated residuum
*   IPARAM- array [1..SZMGRI] of integer
*           Integer parameter structure of multigrid solver
*   DPARAM- array [1..SZMGRI] of integer
*           Double precision parameter structure of multigrid solver
*   IDATA - array [1..*] of integer
*           User defined integer array
*   DDATA - array [1..*] of double
*           User defined double array
*
* Out (to be returned by DFILT):
*  DX     - updated vector.
***********************************************************************

      SUBROUTINE  M020 (IPARAM,DPARAM, IDATA,DDATA, DX,DB,DD,
     *                  DAX,DPROL,DREST,DPRSM,DPOSM,
     *                  DEX,DEXS,DBC,DSTEP,DFILT)

      IMPLICIT NONE

      INCLUDE 'cout.inc'

      INCLUDE 'cbasicmg.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'm020.inc'
      
      DOUBLE PRECISION DX(*),DB(*),DD(*)
      INTEGER IPARAM(SZ020I),IDATA(*)
      DOUBLE PRECISION DPARAM(SZ020D),DDATA(*)
      EXTERNAL DAX,DPROL,DREST,DPRSM,DPOSM,DEX,DEXS,DBC,DSTEP,DFILT
      
C     local variables:

      INTEGER IASRLN, ITIM, NLMIN, NLMAX, INRM, ICYCLE, IFILT, NIT0
      INTEGER KOFFX(NNLEV), KOFFB(NNLEV), KOFFD(NNLEV), KNEQ(NNLEV)
      DOUBLE PRECISION TIMIN(SZ020D),DMAX,DEFOLD,R,DSTEPP
      INTEGER I,ILEV,NODEFC,ITE,MTDV,J

C     The queue saves the current residual and the two previous
C     and the two previous residuals

      DOUBLE PRECISION RESQUE(32)

C     Length of the queue of last residuals for the computation of
C     the asymptotic convergence rate

      IASRLN = MAX(MIN(32,IPARAM (OIASRLN)),1)
      
C     Put some variables into local variables for faster access.
C     This is typically optimized away by the compiler...

      ITIM = IPARAM   (OTIM)
      INRM = IPARAM   (OINRM)
      ICYCLE = IPARAM (OICYCLE)
      IFILT = IPARAM  (OFILT)
      NODEFC = IPARAM (ONODEFC)
      CALL LCP3(IPARAM(OKOFFX),KOFFX,NNLEV)
      CALL LCP3(IPARAM(OKOFFB),KOFFB,NNLEV)
      CALL LCP3(IPARAM(OKOFFD),KOFFD,NNLEV)
      CALL LCP3(IPARAM(OKNEQ),KNEQ,NNLEV)

C     Iteration when the residuum is printed:

      MTDV = MAX(1,IPARAM(OMTRMRS))
      
C     minimum number of MG-steps:

      NIT0 = MAX(IPARAM(ONITMIN),0)
      
C     If timing information should be computed, initialise the variables
C     in the parameter block and measure the start time.
C     Timing information is always measured by stopping a start- and an
C     end-timestep. We use a temporary array TIMIN for measuring the start
C     timings, while after measuring the end-timings, the difference
C     of start- and end-time is added to the final timing information
C     in the output block.
      
      IF (ITIM.NE.0) THEN
        DPARAM(OTMTOT ) = 0D0
        DPARAM(OTMMG  ) = 0D0
        DPARAM(OTMPROL) = 0D0
        DPARAM(OTMREST) = 0D0
        DPARAM(OTMDEF ) = 0D0
        DPARAM(OTMSMTH) = 0D0
        DPARAM(OTMCGSL) = 0D0
        DPARAM(OTMFILT) = 0D0
        DPARAM(OTMBC  ) = 0D0                              
        DPARAM(OTMCGC ) = 0D0
        CALL LCL1(TIMIN,SZ020D)
        CALL GTMAUX (TIMIN,DPARAM,OTMTOT,0)
      ENDIF

C     initialise output parameters to standard values

      IPARAM(OSTATUS) = 0
      IPARAM(OITE)    = 0
      DPARAM(ODEFINI) = 0D0
      DPARAM(ODEFFIN) = 0D0
      DPARAM(ORHO   ) = 0D0
      DPARAM(ORHOASM) = 0D0
      
C     Clear queue with old residuals
      
      CALL LCL1(RESQUE,IASRLN)
      
C     Check the cycle

      IF (ICYCLE.LT.0) THEN
      
C       Wrong parameters

        IPARAM(OSTATUS) = 1
        IF (IPARAM(OMSGTRM).GE.1) THEN
          WRITE (MTERM,'(A)') 'M020: Invalid cycle'
        END IF
        GOTO 99999
        
      END IF
      
      NLMIN = IPARAM(ONLMIN)
      NLMAX = IPARAM(ONLMAX)
      
      IF (NLMAX.LT.NLMIN.OR.NLMIN.LE.0.OR.NLMAX.GT.NNLEV) THEN
        IPARAM(OSTATUS) = 1
        IF (IPARAM(OMSGTRM).GE.1) THEN
          WRITE (MTERM,'(A)') 'M020: invalid NLMIN/NLMAX'
        END IF
        GOTO 99999
      ENDIF

C     The NEQ entry in the solver structure is initialized to the
C     number of equations on the finest level. It will stay so till the
C     end of the routine, as this information is actually not used by
C     this type of solver.

      IPARAM(ONEQ) = KNEQ(NLMAX)

C     Set the starting address of the auxiliary array for the callback
C     routines to DD on finest level. This vector is used in every
C     MG sweep only once: To calculate the defect, smooth it and
C     restrict it to the coarser level. Afterwards it's no more used.
C     So a couple of callback routines can use it for intermediate
C     calculations, which saves some memory!
C     We set KCBAUX<>0 only if it's save to use it. At the moment,
C     it's not -- set it to 0.

      IPARAM(OKCBAUX) = 0
      
C     We start on the maximum level. Set this here, because the callback
C     routines might need it.
      
      IPARAM (OILEV) = NLMAX      
      
C     Test, if the RHS-vector is 0; in this case the solution is
C     also 0 -> special case.
C     This is done by checking the maximum norm of the vector.

      IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMDEF,0)

      CALL LLI1(DB(1+KOFFB(NLMAX)),KNEQ(NLMAX),DMAX,I)

      IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMDEF,1)

      IF (DMAX.LE.DPARAM(OVECZER)) THEN
      
C       final defect is 0, as initialised in the output variable above

        CALL LCL1(DX(1+KOFFX(NLMAX)),KNEQ(NLMAX))
        GOTO 1000
        
      ENDIF
      
C     The same way we check the initial defect for being zero.

C     Test if there's only one level. In this case directly activate
C     that coarse grid solver, which is responsible if there's only
C     one level.

      IF (NLMIN.EQ.NLMAX) THEN
        IF (IPARAM(OMSGTRM).GT.1) THEN
           WRITE (MTERM,'(A)') 'M020: Only one level;'//
     *                         ' switching back to standard solver.'
        END IF
        
C       Calculate initial defect

        IF (NODEFC.EQ.0) THEN

          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMDEF,0)

          CALL LCP1(DB(1+KOFFB(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX))
          CALL DAX (DX(1+KOFFX(NLMAX)),DD(1+KOFFD(NLMAX)),
     *              KNEQ(NLMAX),-1D0,1D0,IPARAM,DPARAM,IDATA,DDATA)
          CALL LL21(DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),DPARAM(ODEFINI))
     
C         Scaling for the vector (1111...) to have norm 1 (weighted l2-norm)

          IF (INRM.GT.0) 
     *      DPARAM(ODEFFIN) = DPARAM(ODEFINI) / SQRT (DBLE(KNEQ(NLMAX)))

          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMDEF,1)
  
        END IF
        
C       activate DAXS coarse grid solver

        IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMCGSL,0)
        IPARAM(OITE) = NIT0
        CALL DEXS(DX(1+KOFFX(NLMAX)),DB(1+KOFFB(NLMAX)),
     *            DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),
     *            DPARAM(ORHO),IPARAM(OITE),IPARAM,DPARAM,IDATA,DDATA)
        IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMCGSL,1)
        
C       Calculate final defect

        IF (NODEFC.EQ.0) THEN

          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMDEF,0)

          CALL LCP1(DB(1+KOFFB(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX))
          CALL DAX (DX(1+KOFFX(NLMAX)),DD(1+KOFFD(NLMAX)),
     *              KNEQ(NLMAX),-1D0,1D0,IPARAM,DPARAM,IDATA,DDATA)
          CALL LL21(DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),DPARAM(ODEFFIN))
          
C         Scaling for the vector (1111...) to have norm 1 (weighted l2-norm)

          IF (INRM.GT.0) 
     *      DPARAM(ODEFFIN) = DPARAM(ODEFFIN) / SQRT (DBLE(KNEQ(NLMAX)))

          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMDEF,1)

C         In this case we can't get the asymptotic residual, otherwise 
C         we would have to change the code of the coarse grid solver 
C         itself.

          DPARAM(ORHOASM) = DPARAM(ORHO)

        END IF
       
C       That's it, single-grid solver = complete solver, finished
       
        GOTO 1000
       
      ENDIF

C     The KIT/KIT0-variables in the IPARAM block count the current sweep
C     on each level, thus realising the F/V/W-cycle without recursion.
C
C     KIT0 is initialised with the number of MG-sweeps on every level.
C     KIT counts backwards the number of MG-sweeps left on each level.
C     An F-cycle is initialised the same way as the W-cycle, but later
C     in decreasing the KIT()-entries it's handled differently.

      IPARAM(OKIT0+NLMAX-1) = 1
      
      DO I = NLMIN+1,NLMAX-1
        IF (ICYCLE.EQ.0) THEN
          IPARAM(OKIT0+I-1) = 2
        ELSE
          IPARAM(OKIT0+I-1) = ICYCLE
        ENDIF
      END DO

C     On start of the algorithm perform the first filtering:

      IF (IFILT.NE.0) THEN
        IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,0)
        CALL DFILT(DX(1+KOFFX(NLMAX)),KNEQ(NLMAX),1,IPARAM,DPARAM,
     *             IDATA,DDATA)
        IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,1)
      END IF

      IF (NODEFC.EQ.0) THEN
      
C       After the filtering calculate the initial defect:

        IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMDEF,0)

        CALL LCP1(DB(1+KOFFB(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX))
        CALL DAX(DX(1+KOFFX(NLMAX)),DD(1+KOFFD(NLMAX)),
     *           KNEQ(NLMAX),-1D0,1D0,IPARAM,DPARAM,IDATA,DDATA)
        CALL LL21(DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),DPARAM(ODEFINI))
        
C       Scaling for the vector (1111...) to have norm 1 (weighted l2-norm)

        IF (INRM.GT.0) 
     *    DPARAM(ODEFINI) = DPARAM(ODEFINI) / SQRT (DBLE(KNEQ(NLMAX)))

        IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMDEF,1)

C       DEFOLD always saves the defect from the last step. The current
C       defect is directly calculated into the output block of DPARAM.
C       In the first iteration the "previous" defect is initialised
C       by the initial defect.
      
        DEFOLD = DPARAM(ODEFINI)
        DPARAM(ODEFFIN) = DPARAM(ODEFINI)
       
C       Print out the initial residuum

        IF (IPARAM(OMSGTRM).GE.2) THEN
          WRITE (MTERM,'(A,I7,A,D25.16)') 
     *     'M020: Iteration ',0,',  !!RES!! = ',DPARAM(ODEFINI)
        END IF

C       Check if out initial defect is zero. This may happen if 
C       the filtering routine filters "everything out"!
C       In that case we can directly stop our computation.

        IF ( (DPARAM(ODEFINI).LT.DPARAM(OVECZER)) ) THEN
     
C         final defect is 0, as initialised in the output variable above

          CALL LCL1(DX(1+KOFFX(NLMAX)),KNEQ(NLMAX))
          GOTO 1000
          
        END IF
        
C       Initialize the queue of the last residuals with the 
C       initial defect

        DO I=1,IASRLN
          RESQUE(I) = DPARAM(ODEFINI)
        END DO  

      END IF

C     Start multigrid iteration; perform at most IPARAM(ONITMAX) iterations.

      DO ITE = 1, IPARAM(ONITMAX)
      
        IPARAM(OCRITE) = ITE
        
C       Initialize level counts for all levels.
C       Transfer the KIT0-array to the KIT-array; it will be decreased 
C       consecutively, until 0 is reached.

        DO I = NLMIN,NLMAX
          IPARAM(OKIT+I-1) = IPARAM(OKIT0+I-1)
        END DO

C       Start on the maximum level

        IPARAM (OILEV) = NLMAX
        ILEV = NLMAX
      
C       Perform the filtering for the current solution before the MG sweep.

        IF (IFILT.NE.0) THEN
          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,0)
          CALL DFILT(DX(1+KOFFX(NLMAX)),KNEQ(NLMAX),2,IPARAM,DPARAM,
     *               IDATA,DDATA)
          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,1)
        END IF
      
C       On the maximum level we already built out defect vector. If we are
C       on a lower level than NLMAX, perform smoothing+restriction down to the
C       lowest level NLMIN.
C
C       Crippled WHILE-DO-loop to allow decreasing of ILEV in the inner
C       of the loop, not at the end.
      
110     IF (ILEV.NE.NLMIN) THEN

C       First perform filtering of the solution before smoothing.

          IF (IFILT.NE.0) THEN
            IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,0)
            CALL DFILT(DX(1+KOFFX(ILEV)),KNEQ(ILEV),3,IPARAM,DPARAM,
     *                 IDATA,DDATA)
            IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,1)
          END IF

C         Perform the pre-smoothing with the current solution vector

          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMSMTH,0)
          IF (IPARAM(OKPRSM+ILEV-1).GT.0) 
     *      CALL DPRSM(DX(1+KOFFX(ILEV)),DB(1+KOFFB(ILEV)),
     *              DD(1+KOFFD(ILEV)),KNEQ(ILEV),IPARAM(OKPRSM+ILEV-1),
     *              IPARAM,DPARAM,IDATA,DDATA)
          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMSMTH,1)

C         Perform filtering

          IF (IFILT.NE.0) THEN
            IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,0)
            CALL DFILT(DX(1+KOFFX(ILEV)),KNEQ(ILEV),4,IPARAM,DPARAM,
     *                 IDATA,DDATA)
            IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,1)
          END IF

C         Build the defect vector

          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMDEF,0)
          CALL LCP1(DB(1+KOFFB(ILEV)),DD(1+KOFFD(ILEV)),KNEQ(ILEV))
          CALL DAX(DX(1+KOFFX(ILEV)),DD(1+KOFFD(ILEV)),
     *             KNEQ(ILEV),-1D0,1D0,IPARAM,DPARAM,IDATA,DDATA)
          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMDEF,1)

C         Filter the defect vector

          IF (IFILT.NE.0) THEN
            IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,0)
            CALL DFILT(DD(1+KOFFD(ILEV)),KNEQ(ILEV),5,IPARAM,DPARAM,
     *                 IDATA,DDATA)
            IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,1)
          END IF

C         Go down one level

          ILEV = ILEV-1
          IPARAM (OILEV) = ILEV

C         Restriction of the defect. The restricted defect is placed
C         in DB as the right hand side of the lower level.

          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMREST,0)
          CALL DREST(DD(1+KOFFD(ILEV+1)),DB(1+KOFFB(ILEV)),
     *               IPARAM,DPARAM,IDATA,DDATA)
          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMREST,1)
       
C         Filter the restricted defect vector

          IF (IFILT.NE.0) THEN
            IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,0)
            CALL DFILT(DB(1+KOFFB(ILEV)),KNEQ(ILEV),6,IPARAM,DPARAM,
     *                 IDATA,DDATA)
            IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,1)
          END IF

C         From now on (i.e. after the first restriction)
C         it's save to use "DD(NLMAX)" as auxiliary vector:

          IPARAM(OKCBAUX) = KOFFD(NLMAX)

C         Choose zero as initial vector on lower level. Implement boundary
C         conditions into the just calculated right hand side.

          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMBC,0)
          CALL LCL1(DX(1+KOFFX(ILEV)),KNEQ(ILEV))
          CALL DBC(DB(1+KOFFB(ILEV)),KNEQ(ILEV),IPARAM,DPARAM,
     *             IDATA,DDATA)
          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMBC,1)
        
C         Perform the filtering on the start solution

          IF (IFILT.NE.0) THEN
            IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,0)
            CALL DFILT(DX(1+KOFFX(ILEV)),KNEQ(ILEV),2,IPARAM,DPARAM,
     *                 IDATA,DDATA)
            IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,1)
          END IF

C         If we are not on the lowest level, repeat the smoothing of 
C         the solution/restriction of the new defect:

          GOTO 110

        END IF

C       The previous IF/GOTO sweep ensures that we are on the lowest level now.
C       In DD there is the defect on the lowest level, DX is filled with zero
C       (plus eventually filtering).

C       Do some probably filtering for coarse grid solution and RHS vector

        IF (IFILT.NE.0) THEN
          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,0)
          CALL DFILT(DB(1+KOFFB(ILEV)),KNEQ(ILEV),7,IPARAM,DPARAM,
     *               IDATA,DDATA)
          CALL DFILT(DX(1+KOFFX(ILEV)),KNEQ(ILEV),8,IPARAM,DPARAM,
     *               IDATA,DDATA)
          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,1)
        END IF

C       Solve the system on lowest level:

        IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMCGSL,0)
        I = NIT0
        CALL DEX(DX(1+KOFFX(NLMIN)),DB(1+KOFFB(NLMIN)),
     *           DD(1+KOFFD(NLMIN)),KNEQ(NLMIN),R,I,IPARAM,DPARAM,
     *           IDATA,DDATA)
        IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMCGSL,1)

C       Filter the solution on lowest level

        IF (IFILT.NE.0) THEN
          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,0)
          CALL DFILT(DX(1+KOFFX(ILEV)),KNEQ(ILEV),9,IPARAM,DPARAM,
     *               IDATA,DDATA)
          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,1)
        END IF

C       Prolongate the solution vector, perform the coarse grid
C       correction and realise the MG-cycles until we have reached
C       the fine grid again:

130     IF (ILEV.NE.NLMAX) THEN
    
C         go to the next higher level
    
          ILEV = ILEV+1
          IPARAM (OILEV)=ILEV
          
C         When we reach NLMAX, KCBAUX must no more be used:

          IF (ILEV.EQ.NLMAX) THEN
            IPARAM(OKCBAUX) = 0
          END IF

C         First perform filtering to the non-prolongated update vector

          IF (IFILT.NE.0) THEN
            IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,0)
            CALL DFILT(DX(1+KOFFX(ILEV)),KNEQ(ILEV),10,IPARAM,DPARAM,
     *                 IDATA,DDATA)
            IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,1)
          END IF

C         Prolongate the update vector; DPROL  returns  DD:=PROL(DX)

          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMPROL,0)
          CALL DPROL(DX(1+KOFFX(ILEV-1)),DD(1+KOFFD(ILEV)),
     *               IPARAM,DPARAM,IDATA,DDATA)
          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMPROL,1)

C         implement boundary conditions into the prolongated vector
      
          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMBC,0)
          CALL DBC(DD(1+KOFFD(ILEV)),KNEQ(ILEV),IPARAM,DPARAM,
     *             IDATA,DDATA)
          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMBC,1)

C         Perform filtering of the prolongated update vector

          IF (IFILT.NE.0) THEN
            IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,0)
            CALL DFILT(DD(1+KOFFD(ILEV)),KNEQ(ILEV),11,IPARAM,DPARAM,
     *                 IDATA,DDATA)
            IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,1)
          END IF

C         Calculate the step length parameter for the coarse grid correction

          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMCGC,0)

          CALL DSTEP(DX(1+KOFFX(ILEV)),DD(1+KOFFD(ILEV)),
     *               DB(1+KOFFB(ILEV)),KNEQ(ILEV),DSTEPP,IPARAM,DPARAM,
     *               IDATA,DDATA)
     
C         Perform the coarse grid correction by adding the coarse grid 
C         solution (with the calculated step-length parameter) to
C         the current solution
     
          CALL LLC1(DD(1+KOFFD(ILEV)),DX(1+KOFFX(ILEV)),KNEQ(ILEV),
     *              DSTEPP,1D0)

          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMCGC,1)

C         Perform filtering of the updated solution vector before
C         post-smoothing

          IF (IFILT.NE.0) THEN
            IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,0)
            CALL DFILT(DX(1+KOFFX(ILEV)),KNEQ(ILEV),12,IPARAM,DPARAM,
     *                 IDATA,DDATA)
            IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,1)
          END IF

C         Perform the post-smoothing with the current solution vector

          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMSMTH,0)
          IF (IPARAM(OKPOSM+ILEV-1).GT.0) 
     *      CALL DPOSM(DX(1+KOFFX(ILEV)),DB(1+KOFFB(ILEV)),
     *               DD(1+KOFFD(ILEV)),KNEQ(ILEV),IPARAM(OKPOSM+ILEV-1),
     *               IPARAM,DPARAM,IDATA,DDATA)
          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMSMTH,1)

C         Filter the currect solution vector after post-smoothing

          IF (IFILT.NE.0) THEN
            IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,0)
            CALL DFILT(DX(1+KOFFX(ILEV)),KNEQ(ILEV),13,IPARAM,DPARAM,
     *                 IDATA,DDATA)
            IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,1)
          END IF

C         Update the iteration counter(s) for realising the MG-cycle(s).
C         Then either jump to 130 to perform the next prolongation or
C         jump to 110 to do perform a next MG sweep on the current
C         level.
C
C         Here ICYCLE defines how the KIT()-entry is updated.
C         For a W-cycle it's resetted to 2 if the sweep is fulfilled
C         on the current level, for F-cycle it's set to 1 to not
C         perform more that 1 cycle on the current level anymore.
        
          J = OKIT
          IPARAM(OKIT+ILEV-1) = IPARAM(OKIT+ILEV-1)-1
          IF (IPARAM(OKIT+ILEV-1).LE.0) THEN
            IF (ICYCLE.EQ.0) THEN
              IPARAM(OKIT+ILEV-1) = 1
            ELSE
              IPARAM(OKIT+ILEV-1) = IPARAM(OKIT0+ILEV-1)
            END IF
            GOTO 130
          ELSE
            GOTO 110
          END IF

        END IF

C       We have (hopefully) successfully performed one MG-sweep, starting
C       and ending on the finest level. As we are now on the finest level
C       again, we can update our defect vector to test the current
C       residuum...

C       But first perform some possible filtering with the current
C       solution:

        IF (IFILT.NE.0) THEN
          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,0)
          CALL DFILT(DX(1+KOFFX(ILEV)),KNEQ(ILEV),14,IPARAM,DPARAM,
     *               IDATA,DDATA)
          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,1)
        END IF

        IF (NODEFC.EQ.0) THEN

C         Calculate the residuum and its norm; the result can be
C         found in DD:

          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMDEF,0)
          CALL LCP1(DB(1+KOFFB(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX))
          CALL DAX(DX(1+KOFFX(NLMAX)),DD(1+KOFFD(NLMAX)),
     *             KNEQ(NLMAX),-1D0,1D0,IPARAM,DPARAM,IDATA,DDATA)

C         Filter the residuum before calculating the norm!

          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,0)
          CALL DFILT(DD(1+KOFFD(ILEV)),KNEQ(ILEV),15,IPARAM,DPARAM,
     *               IDATA,DDATA)
          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMFILT,1)

          CALL LL21(DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),DPARAM(ODEFFIN))
C         Scaling for the vector (1111...) to have norm 1 (weighted l2-norm)
          IF (INRM.GT.0) 
     *      DPARAM(ODEFFIN) = DPARAM(ODEFFIN) / SQRT (DBLE(KNEQ(NLMAX)))

          IF (ITIM.GT.1) CALL GTMAUX (TIMIN,DPARAM,OTMDEF,1)
       
C         Shift the queue with the last residuals and add the new
C         residual to it

          DO I=1,IASRLN-1
            RESQUE(I) = RESQUE(I+1)
          END DO  
          RESQUE(IASRLN) = DPARAM(ODEFFIN)
         
C         Check if the iteration is diverging - by checking the absolute
C         as well as the relative residuum.
C         Use "not" instead of "ge" because this handles NaN/Infinity-cases
C         more robust!

          IF ( .NOT. ((DPARAM(ODEFFIN).LT.DPARAM(ODIVABS)).AND.
     *         (DPARAM(ODEFFIN).LT.DPARAM(ODEFINI)*DPARAM(ODIVREL))) ) 
     *    THEN
            IPARAM(OSTATUS) = 2
            GOTO 1000
          END IF

C         Ok, we are still alive. Check if we are even convergent

          IF (ITE.GE.NIT0) THEN

C           jump out of the loop if we reached the convergence criterion;
C           in standard programming languages i would use a BREAK here :-)
          
            IF ( (DPARAM(OEPSABS).EQ.0D0) .AND.
     *           (DPARAM(ODEFFIN).LT.DPARAM(ODEFINI)*DPARAM(OEPSREL) ) )
     *        GOTO 1000

            IF ( (DPARAM(ODEFFIN).LT.DPARAM(OEPSABS)).AND.
     *           (DPARAM(OEPSREL).EQ.0D0) ) 
     *        GOTO 1000
          
            IF ( (DPARAM(ODEFFIN).LT.DPARAM(OEPSABS)).AND.
     *           (DPARAM(ODEFFIN).LT.DPARAM(ODEFINI)*DPARAM(OEPSREL)) ) 
     *        GOTO 1000
            
          END IF

          IF ((IPARAM(OMSGTRM).GE.2).AND.(MOD(ITE,MTDV).EQ.0)) THEN
            WRITE (MTERM,'(A,I7,A,D25.16)') 
     *         'M020: Iteration ',ITE,',  !!RES!! = ',DPARAM(ODEFFIN)
          END IF

        END IF

C       No, not yet convergent - we have to perform the next sweep.
C       Save the current defect as "old" defect

        DEFOLD = DPARAM(ODEFFIN)

      END DO ! ITE

C     Ok, the multigrid sweep has finished - either successfully (i.e.
C     convergent), or because the maximum number of iterations has
C     been reached.
C     Remark that at this point, ITE=NITMAX+1 -- but in CRITE the
C     correct number of iterations is noted.

1000  CONTINUE

C     Finish - either with an error or converged.
C
C     Calculation of statistical data is only done if we did not use
C     the single-gris solver:

      IF (NLMIN.NE.NLMAX) THEN

C       Print the last residuum, if we finished before reaching the
C       maximum number of iterations. The DO loop has the property
C       that ITE=NITMAX+1 if it runs through completely!

        IF ((IPARAM(OMSGTRM).GE.2).AND.
     *      (ITE.GE.1).AND.(ITE.LE.IPARAM(ONITMAX))) THEN
          WRITE (MTERM,'(A,I7,A,D25.16)') 
     *          'M020: Iteration ',ITE,
     *          ',  !!RES!! = ',DPARAM(ODEFFIN)
        END IF

C       We now gather and print some statistical data, before we close
C       this algorithm. From now on we use the real number of iterations,
C       which is counted in CRITE:

        IPARAM(OITE)    = IPARAM(OCRITE)

        DPARAM(ORHO)    = 0D0
        DPARAM(ORHOASM) = 0D0
        
        IF (NODEFC.EQ.0) THEN
        
C         Don't calculate anything if the final residuum is out of 
C         bounds - would result in NaN's,...
        
          IF (DPARAM(ODEFFIN).LT.1D99) THEN
        
C           Calculate asymptotic convergence rate
        
            IF (RESQUE(1).GE.1D-70) THEN
              I = MIN(IPARAM(OITE),IASRLN-1)
              DPARAM(ORHOASM) = 
     *              (DPARAM(ODEFFIN)/RESQUE(1))**(1D0/DBLE(I))
            END IF

C           If the initial defect was zero, the solver immediately
C           exits - and so the final residuum is zero and we performed
C           no steps; so the resulting multigrid convergence rate stays zero.
C           In the other case the multigrid convergence rate computes as
C           (final defect/initial defect) ** 1/nit :

            IF (DPARAM(ODEFINI).GT.DPARAM(OVECZER)) THEN
              DPARAM(ORHO) = (DPARAM(ODEFFIN) / DPARAM(ODEFINI)) ** 
     *                       (1D0/DBLE(IPARAM(OITE)))
            END IF
            
C           If the convergence rate is really > 0, we treat the iteration
C           as diverging - the error was getting larger than the original
C           one!

            IF (DPARAM(ORHO).GT.1D0) THEN
              IPARAM(OSTATUS) = 2
            END IF
            
          END IF
          
        END IF
        
      END IF
        
C     Print statistical data
        
      IF (NODEFC.EQ.0) THEN
      
C       Don't calculate anything if the final residuum is out of 
C       bounds - would result in NaN's,...
      
        IF (DPARAM(ODEFFIN).LT.1D99) THEN

          IF (IPARAM(OMSGTRM).GE.2) THEN
            WRITE (MTERM,'(A)') ''
            IF (NLMIN.NE.NLMAX) THEN
              WRITE (MTERM,'(A)') 'Multigrid statistics:'
            ELSE
              WRITE (MTERM,'(A)') 'Single grid solver statistics:'
            END IF
            WRITE (MTERM,'(A)') ''
            WRITE (MTERM,'(A,I5)')     'Iterations              : ',
     *            IPARAM(OITE)
            WRITE (MTERM,'(A,D24.12)') '!!INITIAL RES!!         : ',
     *            DPARAM(ODEFINI)
            WRITE (MTERM,'(A,D24.12)') '!!RES!!                 : ',
     *            DPARAM(ODEFFIN)
            IF (DPARAM(ODEFINI).GT.DPARAM(OVECZER)) THEN     
              WRITE (MTERM,'(A,D24.12)') '!!RES!!/!!INITIAL RES!! : ',
     *              DPARAM(ODEFFIN) / DPARAM(ODEFINI)
            ELSE
              WRITE (MTERM,'(A,D24.12)') '!!RES!!/!!INITIAL RES!! : ',
     *              0D0
            END IF
            WRITE (MTERM,'(A)') ''
            WRITE (MTERM,'(A,D24.12)') 'Rate of convergence     : ',
     *            DPARAM(ORHO)
            WRITE (MTERM,'(A)') ''
          END IF

          IF (IPARAM(OMSGTRM).EQ.1) THEN
            WRITE (MTERM,'(A,I5,A,D24.12)') 
     *            'Multigrid: Iterations/Rate of convergence: ',
     *            IPARAM(OITE),' /',DPARAM(ORHO)
          END IF
        ELSE
        
C         DEF=Infinity; RHO=Infinity, set to 1

          DPARAM(ORHO) = 1D0
          DPARAM(ORHOASM) = 1D0
          
        END IF  
        
      END IF

99999 CONTINUE

C     Temporary array no more available

      IPARAM(OKCBAUX) = 0

C     Gather some timing information, finish

      IF (ITIM.GT.0) THEN
        CALL GTMAUX (TIMIN,DPARAM,OTMTOT,1)
        DPARAM(OTMMG  ) = DPARAM(OTMTOT ) - DPARAM(OTMPROL) 
     *                  - DPARAM(OTMREST) - DPARAM(OTMDEF )
     *                  - DPARAM(OTMSMTH) - DPARAM(OTMCGSL)
     *                  - DPARAM(OTMFILT) - DPARAM(OTMBC  )
     *                  - DPARAM(OTMCGC )
      END IF

      END

************************************************************************
* M020 initialization
*
* The following routine can be used to initialise the IPARAM/DPARAM
* array structures with default values for the computation. After
* calling this routine, the user has to do the following initialisations
* before calling M020:
* - initialise NLMIN and NLMAX
* - initialise KOFFX, KOFFB, KOFFD, KNEQ, KPRSM, KPOSM
* - initialise RHS vector and start vector
* - initialise any user-defined variables attached to the structures
*   IPARAM/DPARAM
*
* In:
*   -
* Out:
*   IPARAM  - array [1..SZMGRI] of integer
*             Integer parameter structure
*   DPARAM  - array [1..SZMGRI] of integer
*             Double precision parameter structure
************************************************************************

      SUBROUTINE INM020 (IPARAM, DPARAM)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'm020.inc'

      INTEGER IPARAM(SZ020I)
      DOUBLE PRECISION DPARAM(SZ020D)
      
C     Clear the structures

      CALL LCL3 (IPARAM,SZ020I)
      CALL LCL1 (DPARAM,SZ020D)
      
C     Initialize standard-parameters:

      CALL INGSLV (IPARAM,DPARAM)
      
C     Set the non-zero standard values
      
      IPARAM (OICYCLE) = 0
      DPARAM (OSTPMIN) = 1D0
      DPARAM (OSTPMAX) = 1D0
      
C     Standard solver tag is 11 for MG

      IPARAM (OSLTAG) = 11

      END
      
************************************************************************
* Prepare M020 vectors
*
* This is another preparation routine for the IPARAM/DPARAM array
* structures. It can be called after INM020 to initialise
* the KOFFX/KOFFB/KOFFD/KPRSM/KPOSM subarrays with starting addresses
* of arrays in DWORK and the number of smoothing steps, resp.
*
* PRM020 accepts a couple of handles generated by the dynamic
* memory management. It initialises KOFFX/KOFFB/KOFFD according to
* the starting addresses of the corresponding arrays in DWORK.
* The later call to M020 has then to be made using
* DX=DB=DD=DWORK(1):
*
*  CALL M020 (IPARAM,DPARAM,DWORK(1),DWORK(1),DWORK(1),...)
*
* After calling this routine, the caller must make sure that there's
* no memory deallocation between this function and the cann to M020!
* Otherwise the starting addresses of the solution/RHS/aux. arrays
* may not be valid anymore!
*
* IPARAM/DPARAM should be filled with 0 before calling this routine,
* as here only the nonzero parameters are set.
*
* Warning: The routine does not initialize KNEQ!
*          This must still be done by the caller!
*
* In:
*   NLMIN  : minimum level, where the coarse grid solver should solve;
*            >= 1.
*   NLMAX  : maximum level, where the solution should be computed; 
*            <= NNLEV!
*   NPRSM   - number of pre-smoothing steps on each level.
*             =-1: don't initialise
*   NPOSM   - number of post-smoothing steps on each level.
*             =-1: don't initialise
*   LOFFX   - array [1..NNLEV] of integer
*             Array of handles to the DX-vectors
*   LOFFB   - array [1..NNLEV] of integer
*             Array of handles to the DB-vectors
*   LOFFD   - array [1..NNLEV] of integer
*             Array of handles to the DD-vectors
*
* Out:
*   The IPARAM/DPARAM structure will be modified in the following
*   variables:
*
*   KOFFX   - array [1..NNLEV] of integer
*             Starting offsets of the DX-vectors relative to DWORK(1)
*   KOFFB   - array [1..NNLEV] of integer
*             Starting offsets of the DX-vectors relative to DWORK(1)
*   KOFFD   - array [1..NNLEV] of integer
*             Starting offsets of the DX-vectors relative to DWORK(1)
*
*   If NPRSM <> -1:
*   KPRSM   - array [1..NNLEV] of integer
*             Number of pre-smoothing steps on each level
*
*   If NPOSM <> -1:
*   KPOSM   - array [1..NNLEV] of integer
*             Number of post-smoothing steps on each level
************************************************************************

      SUBROUTINE PRM020 (IPARAM, DPARAM,NLMIN,NLMAX, 
     *                   NPRSM, NPOSM, LOFFX, LOFFB, LOFFD)

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'm020.inc'

      INTEGER IPARAM(SZ020I)
      DOUBLE PRECISION DPARAM(SZ020D)
      INTEGER LOFFX(NNLEV),LOFFB(NNLEV),LOFFD(NNLEV),NPRSM,NPOSM
      INTEGER NLMIN,NLMAX

      INTEGER I
      
      IF (NPRSM.GT.0) THEN
        DO I=0,NNLEV-1
          IPARAM (OKPRSM+I) = NPRSM
        END DO
      END IF
      
      IF (NPOSM.GT.0) THEN
        DO I=0,NNLEV-1
          IPARAM (OKPOSM+I) = NPOSM
        END DO
      END IF
      
      DO I=1,NNLEV
        IF ((I.GE.NLMIN).AND.(I.LE.NLMAX)) THEN
          IPARAM (OKOFFX+I-1) = L(LOFFX(I))-1
          IPARAM (OKOFFB+I-1) = L(LOFFB(I))-1
          IPARAM (OKOFFD+I-1) = L(LOFFD(I))-1
        ELSE
          IPARAM (OKOFFX+I-1) = 0
          IPARAM (OKOFFB+I-1) = 0
          IPARAM (OKOFFD+I-1) = 0
        END IF
      END DO
      
      END
      