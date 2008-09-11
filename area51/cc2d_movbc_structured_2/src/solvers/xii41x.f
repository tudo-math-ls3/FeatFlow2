************************************************************************
* II01X
*
* BiCGStab solver with extended parameter interface
*
* This is a modified version of the BiCGStab solver, based on XII41N.
* This version uses a modified calling interface. The solver is desiged
* as a Black-Box solver. The user has to provide some callback routines,
* which implement the real intelligence.
*
* The implementation follows the original paper introducing the
* BiCGStab:
*   van der Vorst, H.A.; BiCGStab: A Fast and Smoothly Converging
*   Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems;
*   SIAM J. Sci. Stat. Comput. 1992, Vol. 13, No. 2, pp. 631-644
*
* The solver is configured with an integer and a double precision
* parameter block. As the standard BiCGStab solver makes only use
* of standard parameters in the general solver parameter blocks,
* there is no extension necessary for the standard case.
*
* In the case that there is additional information needed for the
* callback routines (e.g. matrix handles), the caller should append this
* information to the end of the general parameter block. The callback
* routines have then to be adapted to access this information properly.
*
* To use the solver, the caller has to do the following:
* - initialize the solver structure, e.g. with INX41X
* - append necessary information to the end of the structure if necessary
* - call the solver
* - find the results in the output variables and the output
*   parameters
* 
* Apart from the input-parameters in the parameter blocks, the routine
* uses the parameters:
*
* In:
*   IPARAM - array [1..SZMGRI] of integer
*            Integer parameter structure
*   DPARAM - array [1..SZMGRI] of integer
*            Double precision parameter structure
*   IDATA  - array [1..*] of integer
*            User defined integer array; passed to callback routines
*   DDATA  - array [1..*] of double
*            User defined double array; passed to callback routines
*   NEQ    - size of vectors
*   DX     - array [1..NEQ] of double precision
*            Start vector
*   DB     - array [1..NEQ] of double precision
*            Right hand side
*   DR,
*   DR0,
*   DP,
*   DPA,
*   DSA    - array [1..NEQ] of double precision
*            Temporary vectors
*   DCG0C  - SUBROUTINE (DG,NEQ,IPARAM,DPARAM,IDATA,DDATA)
*            Extended preconditioning routine; only used if IPREC-
*            variable in parameter block is set to 2.
*            Results  DG := C**(-1) * DG  for the precondioning    
*            matrix C.
*   YMVMUL - SUBROUTINE YMVMUL (DX,DAX,NEQ,A1,A2,IPARAM,DPARAM,
*                               IDATA,DDATA)
*            Performes matrix vector multiplication
*            DAX = A1*DX + A2*DAX
*   DFILT  - SUBROUTINE (DX,NEQ,IALGP,IPARAM,DPARAM,IDATA,DDATA)
*            Callback-routine, adaptive filtering.
*            Performs filtering on different positions of the
*            algorithm. Is only used if IPARAM(OFILT)<>0.
*
* The following variables from the parameter blocks have a specialized
* meaning:
*   IPREC  - =0  : no preconditioning.
*            <> 0: Call DCG0C for preconditioning
*            If standard callback routines are used for
*            preconditioning, the PRECI/PRECD parameter block must
*            be initialized according to the preconditioner
*   EPSREL : relative stopping criterion. 
*   EPSABS : absolute stopping criterion. 
*            If EPSREL<1D99, the relative stopping criterion will be 
*            used, otherwise the absolute one.
*
* Out:
*   DX     - Updated solution vector
************************************************************************

***********************************************************************
* Description of filtering by DFILT:
* 
* In:
*  DX     - Solution vector; array [1..NEQ] of double
*  NEQ    - length of solution vector
*  IALGP  - Position in the algorithm
*           0=undefined - enforce the filtering
*  IPARAM - array [1..SZMGRI] of integer
*           Integer parameter structure of the solver
*  DPARAM - array [1..SZMGRI] of integer
*           Double precision parameter structure of the solver
*  IDATA  - array [1..*] of integer
*           User defined integer array
*  DDATA  - array [1..*] of double
*           User defined double array
*
* Out (to be returned by DFILT):
*  DX     - updated vector.
***********************************************************************

      SUBROUTINE II01X(IPARAM, DPARAM, IDATA,DDATA,
     *                 NEQ,DX,DB,DR,DR0,DP,DPA,DSA,
     *                 YMVMUL,DCG0C,DFILT)

      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cmem.inc'
      
      INCLUDE 'ssolvers.inc'
      
      INTEGER IPARAM(SZSLVI),IDATA(*)
      DOUBLE PRECISION DPARAM(SZSLVD),DDATA(*)
      
      INTEGER NEQ
      DOUBLE PRECISION DX(NEQ),DB(NEQ)
      DOUBLE PRECISION DR(NEQ),DR0(NEQ),DP(NEQ),DPA(NEQ),DSA(NEQ)
      EXTERNAL YMVMUL,DCG0C,DFILT
      
C     local variables

      INTEGER IASRLN,MTDV,I,NIT0,ITE
      LOGICAL BPREC,BFILT
      DOUBLE PRECISION RHO0,DALPHA,OMEGA0,RES,DEFINI,RHO1,DBETA
      DOUBLE PRECISION OMEGA1, OMEGA2, FR

C     The queue saves the current residual and the two previous
C     and the two previous residuals

      DOUBLE PRECISION RESQUE(32)

      IPARAM(OSTATUS) = 0

      IF (NEQ.EQ.0) THEN
      
C       no linear system !?!

        IPARAM(OSTATUS) = 2
        RETURN
        
      END IF

C     Length of the queue of last residuals for the computation of
C     the asymptotic convergence rate

      IASRLN = MAX(0,MIN(32,IPARAM(OIASRLN)))

C     Minimum number of iterations
 
      NIT0 = MAX(IPARAM(ONITMIN),0)
      
C     Use preconditioning? Filtering?

      BPREC = IPARAM(OPCTAG).NE.0
      BFILT = IPARAM(OFILT).NE.0
      
C     Iteration when the residuum is printed:

      MTDV = MAX(1,IPARAM(OMTRMRS))

C Initialize used vectors with zero
      
      CALL LCL1(DP,NEQ)
      CALL LCL1(DPA,NEQ)

C *** Initialization

      RHO0  =1D0
      DALPHA=1D0
      OMEGA0=1D0

C Matrix-vector multiplication with preconditioning

      CALL LCP1(DB,DR,NEQ)
      CALL YMVMUL(DX,DR,NEQ,-1D0,1D0,IPARAM,DPARAM,IDATA,DDATA)
      IF (BFILT) CALL DFILT(DR,NEQ,0,IPARAM,DPARAM,IDATA,DDATA)
      IF (BPREC) CALL DCG0C(DR,NEQ,IPARAM,DPARAM,IDATA,DDATA)
      CALL LL21(DR,NEQ,RES)

C Scaling for the vektor (1111...) to have norm = 1
      
      IF (IPARAM(OINRM).EQ.0) THEN
        RES = RES / SQRT (DBLE(NEQ))
        IF (.NOT.((RES.GE.1D-99).AND.(RES.LE.1D99))) RES = 0D0
      END IF

      DPARAM(ODEFINI) = RES

C Initialize starting residuum
      
      DEFINI = RES

C initialize the queue of the last residuals with RES

      DO I=1,IASRLN
        RESQUE(I)=RES
      END DO  

C Check if out initial defect is zero. This may happen if the filtering
C routine filters "everything out"!
C In that case we can directly stop our computation.

      IF ( (DPARAM(ODEFINI).LT.DPARAM(OVECZER)) ) THEN
     
C final defect is 0, as initialised in the output variable above
        CALL LCL1(DX,NEQ)
        ITE = 0
        DPARAM(ODEFFIN) = RES
        GOTO 200
          
      END IF

      IF (IPARAM(OMSGTRM).GE.2) THEN
        WRITE (MTERM,'(A,I7,A,D25.16)') 
     *    'II01X: Iteration ',0,',  !!RES!! = ',DPARAM(ODEFINI)
      END IF

      CALL LCP1(DR,DR0,NEQ)

C Perform at most NIT loops to get a new vector

      DO ITE = 1,IPARAM(ONITMAX)
      
        IPARAM(OCRITE) = ITE

        CALL LSP1(DR0,DR,NEQ,RHO1)

        IF (RHO0*OMEGA0.EQ.0D0) THEN
C Should not happen
          IF (IPARAM(OMSGTRM).GE.2) THEN
            WRITE (MTERM,'(A)') 
     *'II01X: Iteration prematurely stopped! Correction vector is zero!'
          END IF

C         Some tuning for the output, then cancel.

          IPARAM(OSTATUS) = 2
          IPARAM(OITE) = ITE-1
          GOTO 220
          
        END IF

        DBETA=(RHO1*DALPHA)/(RHO0*OMEGA0)
        RHO0 =RHO1

        CALL LLC1(DR ,DP,NEQ,1D0,DBETA)
        CALL LLC1(DPA,DP,NEQ,-DBETA*OMEGA0,1D0)

        CALL YMVMUL(DP,DPA,NEQ,1D0,0D0,IPARAM,DPARAM,IDATA,DDATA)
        IF (BFILT) CALL DFILT(DPA,NEQ,0,IPARAM,DPARAM,IDATA,DDATA)
        IF (BPREC) CALL DCG0C(DPA,NEQ,IPARAM,DPARAM,IDATA,DDATA)

        CALL LSP1(DR0,DPA,NEQ,DALPHA)
        
        IF (DALPHA.EQ.0D0) THEN
C We are below machine exactness - we can't do anything more...
C May happen with very small problems with very few unknowns!
          IF (IPARAM(OMSGTRM).GE.2) THEN
            WRITE (MTERM,'(A)') 'II01X: Convergence failed!'
            IPARAM(OSTATUS) = 2
            GOTO 200
          END IF
        END IF
        
        DALPHA=RHO1/DALPHA

        CALL LLC1(DPA,DR,NEQ,-DALPHA,1D0)

        CALL YMVMUL(DR,DSA,NEQ,1D0,0D0,IPARAM,DPARAM,IDATA,DDATA)
        IF (BFILT) CALL DFILT(DSA,NEQ,0,IPARAM,DPARAM,IDATA,DDATA)
        IF (BPREC) CALL DCG0C(DSA,NEQ,IPARAM,DPARAM,IDATA,DDATA)

        CALL LSP1(DSA,DR ,NEQ,OMEGA1)
        CALL LSP1(DSA,DSA,NEQ,OMEGA2)
        IF (OMEGA1.EQ.0D0) THEN
          OMEGA0 = 0D0
        ELSE
          IF (OMEGA2.EQ.0D0) THEN
            IF (IPARAM(OMSGTRM).GE.2) THEN
              WRITE (MTERM,'(A)') 'II01X: Convergence failed!'
              IPARAM(OSTATUS) = 2
              GOTO 200
            END IF
          END IF
          OMEGA0=OMEGA1/OMEGA2
        END IF

        CALL LLC1(DP ,DX ,NEQ,DALPHA,1D0)
        CALL LLC1(DR ,DX ,NEQ,OMEGA0,1D0)

        CALL LLC1(DSA,DR,NEQ,-OMEGA0,1D0)

        CALL LL21(DR,NEQ,FR)

C Scaling for the vektor (1111...) to have norm = 1
        
        IF (IPARAM(OINRM).EQ.0) THEN
          FR = FR / SQRT (DBLE(NEQ))
          IF (.NOT.((FR.GE.1D-99).AND.(FR.LE.1D99))) FR = 0D0
        END IF
     
C shift the queue with the last residuals and add the new
C residual to it

        DO I=1,IASRLN-1
          RESQUE(I)=RESQUE(I+1)
        END DO  
        RESQUE(IASRLN) = FR

        DPARAM(ODEFFIN) = FR
     
C At least perform NIT0 iterations

        IF (ITE.GE.NIT0) THEN
        
C         Both stopping criteria given? Stop if both are fulfilled.

          IF ((DPARAM(OEPSREL).NE.0D0) .AND.
     *        (DPARAM(OEPSABS).NE.0D0)) THEN

            IF (( FR.LE.RES*DPARAM(OEPSREL) ) .AND.
     *          ( FR.LE.DPARAM(OEPSABS) )) GOTO 200

          ELSE IF (DPARAM(OEPSREL).NE.0D0) THEN
                    
C           Use only relative stopping criterion
                    
            IF ( FR.LE.RES*DPARAM(OEPSREL) ) GOTO 200
            
          ELSE
          
C           Use only absolute stopping criterion

            IF ( FR.LE.DPARAM(OEPSABS) ) GOTO 200        
            
          END IF
          
        END IF

C print out the current residuum

        IF ((IPARAM(OMSGTRM).GE.2).AND.(MOD(ITE,MTDV).EQ.0)) THEN
          WRITE (MTERM,'(A,I7,A,D25.16)') 
     *        'II01X: Iteration ',ITE,',  !!RES!! = ',DPARAM(ODEFFIN)
        END IF

      END DO

C Set ITE to NIT to prevent printing of "NIT+1" of the loop was
C completed

      ITE = IPARAM(ONITMAX)

200   CONTINUE

C Finish - either with an error or if converged.
C Print the last residuum.

      IF ((IPARAM(OMSGTRM).GE.2).AND.
     *    (ITE.GE.1).AND.(ITE.LT.IPARAM(ONITMAX))) THEN
        WRITE (MTERM,'(A,I7,A,D25.16)') 
     *        'II01X: Iteration ',ITE,',  !!RES!! = ',DPARAM(ODEFFIN)
      END IF

      IPARAM(OITE) = ITE
      
220   CONTINUE

C Don't calculate anything if the final residuum is out of bounds -
C would result in NaN's,...
      
      IF (DPARAM(ODEFFIN).LT.1D99) THEN
      
C Calculate asymptotic convergence rate
      
        IF (RESQUE(1).GE.1D-70) THEN
          I = MAX(1,MIN(IPARAM(OITE),IASRLN-1))
          DPARAM(ORHOASM) = (DPARAM(ODEFFIN)/RESQUE(1))**(1D0/DBLE(I))
        END IF

C If the initial defect was zero, the solver immediately
C exits - and so the final residuum is zero and we performed
C no steps; so the resulting multigrid convergence rate stays zero.
C In the other case the multigrid convergence rate computes as
C (final defect/initial defect) ** 1/nit :

        IF (DPARAM(ODEFINI).GT.DPARAM(OVECZER)) THEN
          DPARAM(ORHO) = (DPARAM(ODEFFIN) / DPARAM(ODEFINI)) ** 
     *                     (1D0/DBLE(IPARAM(OITE)))
        END IF
        
        IF (IPARAM(OMSGTRM).GE.2) THEN
          WRITE (MTERM,'(A)') ''
          WRITE (MTERM,'(A)') 'BiCGStab statistics:'
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

        END IF

        IF (IPARAM(OMSGTRM).EQ.1) THEN
          WRITE (MTERM,'(A,I5,A,D24.12)') 
     *          'BiCGStab: Iterations/Rate of convergence: ',
     *          IPARAM(OITE),' /',DPARAM(ORHO)
        END IF
        
      ELSE
C DEF=Infinity; RHO=Infinity, set to 1
        DPARAM(ORHO) = 1D0
        DPARAM(ORHOASM) = 1D0
      END IF  

99999 END

************************************************************************
* II21X - BiCGStab smoother with extended parameter interface
*
* This routine implements smoothing operations with BiCGStab,
* as used e.g. in a multigrid algorithm. The calling convention
* supports passing IPARAM/DPARAM to the matrix/vector and
* preconditioning routines. IPARAM/DPARAM is expected to be
* initialized with standard parameters as usual for a solver.
* In contrast to a solver, the stopping criteria parameters
* are not used.
*
* In:
*   IPARAM - array [1..SZMGRI] of integer
*            Integer parameter structure
*   DPARAM - array [1..SZMGRI] of integer
*            Double precision parameter structure
*   IDATA  - array [1..*] of integer
*            User defined integer array; passed to callback routines
*   DDATA  - array [1..*] of double
*            User defined double array; passed to callback routines
*   NEQ    - size of vectors
*   DX     - array [1..NEQ] of double precision
*            Start vector
*   DB     - array [1..NEQ] of double precision
*            Right hand side
*   DR,
*   DR0,
*   DP,
*   DPA,
*   DSA    - array [1..NEQ] of double precision
*            Temporary vectors
*   DCG0C  - SUBROUTINE (DG,NEQ,IPARAM,DPARAM,IDATA,DDATA)                     
*            Extended preconditioning routine; only used if IPREC-
*            variable in parameter block is set to 2.
*            Results  DG := C**(-1) * DG  for the precondioning    
*            matrix C.
*   YMVMUL - SUBROUTINE YMVMUL (DX,DAX,NEQ,A1,A2,IPARAM,DPARAM,
*                               IDATA,DDATA)
*            Performes matrix vector multiplication
*            DAX = A1*DX + A2*DAX
*   DFILT  - SUBROUTINE (DX,NEQ,IALGP,IPARAM,DPARAM)
*            Callback-routine, adaptive filtering.
*            Performs filtering on different positions of the
*            algorithm. Is only used if IPARAM(OFILT)<>0.
*   NIT    - Number of smoothing steps to perform
*
* The following variables from the parameter blocks have a specialized
* meaning:
*   IPREC  - =0  : no preconditioning.
*            <> 0: Call DCG0C for preconditioning
*            If standard callback routines are used for
*            preconditioning, the PRECI/PRECD parameter block must
*            be initialized according to the preconditioner
*   IFILT  : Whether to use filtering with DFILT
*   STATUS : Not used
*   ITE    : Returned numer of performed smoothing steps.
*            <> NIT if an error has happened.
*   NITMIN : Not used
*   NITMAX : Not used
*   INRM   : Not used
*   ITIM   : Not used
*   IASRLN : Not used
*   MSGTRM : Not used
*   MFLTRM : Not used
*   MTRMRS : Not used
*   CRITE  : Not used
*
*   DEFINI : Not used
*   DEFFIN : Not used
*   RHO    : Not used
*   RHOASM : Not used
*   TMTOT  : Not used
*   TMFILT : Not used
*   EPSREL : Not used
*   EPSABS : Not used
*   DIVREL : Not used
*   DIVABS : Not used
*   VECZER : Not used
*
* Out:
*   DX     - Smoothed vector
************************************************************************

      SUBROUTINE II21X(IPARAM, DPARAM, IDATA,DDATA,
     *                 NEQ,DX,DB,DR,DR0,DP,DPA,DSA,
     *                 YMVMUL,DCG0C,DFILT, NIT)

      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      
      INCLUDE 'ssolvers.inc'
      
      INTEGER IPARAM(SZSLVI),IDATA(*)
      DOUBLE PRECISION DPARAM(SZSLVD),DDATA(*)
      
      INTEGER NEQ,NIT
      DOUBLE PRECISION DX(NEQ),DB(NEQ)
      DOUBLE PRECISION DR(NEQ),DR0(NEQ),DP(NEQ),DPA(NEQ),DSA(NEQ)
      DOUBLE PRECISION RHO1, DBETA, OMEGA1, OMEGA2
      EXTERNAL YMVMUL,DCG0C
      
C local variables

      INTEGER ITE
      LOGICAL BPREC,BFILT
      DOUBLE PRECISION RHO0,DALPHA,OMEGA0

C Use preconditioning? Filtering?

      BPREC = IPARAM(OPCTAG).NE.0
      BFILT = IPARAM(OFILT).NE.0

C Initialize used vectors with zero
      
      CALL LCL1(DP,NEQ)
      CALL LCL1(DPA,NEQ)

C *** Initialization
      RHO0  =1D0
      DALPHA=1D0
      OMEGA0=1D0

      CALL LCP1(DB,DR,NEQ)
      CALL YMVMUL(DX,DR,NEQ,-1D0,1D0,IPARAM,DPARAM,IDATA,DDATA)
      IF (BFILT) CALL DFILT(DR,NEQ,0,IPARAM,DPARAM,IDATA,DDATA)
      IF (BPREC) CALL DCG0C(DR,NEQ,IPARAM,DPARAM,IDATA,DDATA)

      CALL LCP1(DR,DR0,NEQ)

C *** Iterative correction
      DO ITE=1,NIT

        IPARAM(OCRITE) = ITE
        IPARAM(OITE) = ITE

        CALL LSP1(DR0,DR,NEQ,RHO1)
C Cancel here if the norm is nearly zero
        IF (RHO0*OMEGA0.LT.1D-99) GOTO 99999

        DBETA=(RHO1*DALPHA)/(RHO0*OMEGA0)
        RHO0 =RHO1

        CALL LLC1(DR ,DP,NEQ,1D0,DBETA)
        CALL LLC1(DPA,DP,NEQ,-DBETA*OMEGA0,1D0)

        CALL YMVMUL(DP,DPA,NEQ,1D0,0D0,IPARAM,DPARAM,IDATA,DDATA)
        IF (BFILT) CALL DFILT(DPA,NEQ,0,IPARAM,DPARAM,IDATA,DDATA)
        IF (BPREC) CALL DCG0C(DPA,NEQ,IPARAM,DPARAM,IDATA,DDATA)

        CALL LSP1(DR0,DPA,NEQ,DALPHA)
        
C Cancel here if the norm is nearly zero
        IF (DALPHA.LT.1D-99) GOTO 99999
        
        DALPHA=RHO1/DALPHA

        CALL LLC1(DPA,DR,NEQ,-DALPHA,1D0)

        CALL YMVMUL(DR,DSA,NEQ,1D0,0D0,IPARAM,DPARAM,IDATA,DDATA)
        IF (BFILT) CALL DFILT(DSA,NEQ,0,IPARAM,DPARAM,IDATA,DDATA)
        IF (BPREC) CALL DCG0C(DSA,NEQ,IPARAM,DPARAM,IDATA,DDATA)

        CALL LSP1(DSA,DR ,NEQ,OMEGA1)
        CALL LSP1(DSA,DSA,NEQ,OMEGA2)
        
C Cancel here if the norm is nearly zero
        IF (OMEGA2.LT.1D-99) GOTO 99999

        OMEGA0=OMEGA1/OMEGA2

        CALL LLC1(DP ,DX ,NEQ,DALPHA,1D0)
        CALL LLC1(DR ,DX ,NEQ,OMEGA0,1D0)

        CALL LLC1(DSA,DR,NEQ,-OMEGA0,1D0)
        
      END DO

99999 END
