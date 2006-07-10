************************************************************************
* Computes the new OMEGA-parameter for the defect correction.
* Updates the system matrix on the finest level for the next iteration.
* Updates the KUx-arrays with the defect correction approach:
*
*   u^(l+1) = u^l + OMEGA * Y
*
* with Y=(KU1,KU2,KP) the solution from the Oseen equation and
* u^l=(KU1OLD,KU2OLD,KPOLD) the old solution.
* Calculates the norms of the changes into DELxx.
*
*    Input:
*      - vectors U,UOLD,F
*      - OMEGA for the new calculation of the nonlinear block A
*      - matrices and pointers A,...,ST,NA
*      - number of equations NU
*    Output:
*      - updated vector U
*      - optimal relaxation parameter OMEGA
*      - maximum relative changes  DELU 
*      - note that D is changed to  D=K*V with correction V=U-UOLD 
*
* Extended calling convention.
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Triangulation structure of current level
*   MATPAR - array [1..SZN2VI] of integer
*            TNS2DMatrixParams structure containing all information about
*            system-/mass-/Laplace matrices on the current level.
*            Must correspond to TRIA
*   NUVP   - length of vectors
*   DY     - array [1..NUVP] of double
*            Update vector Y for the new time step, calculated
*            by solvíng the Oseen equation
*   DU     - array [1..NUVP] of double
*            Current solution that should be updated.
*   DRHS   - array [1..NUVP] of double
*            Right hand side of the nonlinear iteration
*   DD     - array [1..NUVP] of double
*            Defect vector = RHS of the Oseen equation
*
*   DAUX   - array [1..NUVP] of double
*            temporary vector used for some computations
*
*   OMGMIN - minimal OMEGA
*   OMGMAX - maximal OMEGA
*
*   INEUM  - whether there are Neumann boundary components in the
*            problem or not
*   THSTEP - double. Theta-scheme identifier.
*            =0: Forward Euler, =1: Backward Euler, =0.5: Crank Nicolson
*            For Fractional-Step theta scheme this is set to different
*            values, depending on the current step in the scheme.
*            For stationary simulations this parameter must be set
*            to 1 to include the full Laplacian matrix into the
*            system matrix:
*                    [alpha*M + THETA*KST1 ] 
*   INONST - Whether the simulation is stationary or not.
*            As this routine generally builds
*                    [alpha*M + Theta_1*nu*k*L ] u
*            (the linear part of S(u)), this parameter decides whether
*            ALPHA=0 or ALPHA=1, INONST=ALPHA, so whether the mass 
*            matrix is added or not.   
*   ILEV   - Current level. Only used if IASMBL.IPRECA=3 to read the
*            correct matrix into the memory - otherwise ignored.
*   IASMBL : array [1..SZASMI] of integer
*   DASMBL : array [1..SZASMD] of double
*            Integer and double prec. parameter block that controls the
*            discretization. This tells all assembly-routines how to 
*            set up the nonlinearity in the system matrices, which 
*            cubature formula to use, etc.
*
* Out:
*   DU     - Updated iteration vector
*   OMEGA  - new relaxation parameter OMEGA for defect correction
*   DELU,
*   DELP,
*   DELT   - Norm of new solution vectors / relative change of
*            solution vectors - see below at the end of the file.
*
*   DAUX   - is overwritten with some temporary data (cf. COMMON block)
*   DD     - Is overwritten with some temporary data
*
*   KA1 (System matrix)
*          - is rebuild at the point u^l-omegaold*Y
*
*   DTIM   - array [1..5] of double
*            Time necessary for building matrices,...
*            DTIM(1) = TTLC
*            DTIM(2) = TTADF
*            DTIM(3) = TTUPW
*            DTIM(4) = TTBDR
*            DTIM(5) = TTCOR
************************************************************************

      SUBROUTINE XOPTCX (TRIA, MATPAR, NUVP, DY, DU, DRHS, DD, DAUX, 
     *                   THSTEP, INONST, ILEV,
     *                   IASMBL,DASMBL,
     *                   OMGMIN, OMGMAX, 
     *                   OMEGA, DELU, DELP, DTIM)

      IMPLICIT NONE
      
C common blocks

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'stria.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'
      
      INCLUDE 'sassembly.inc'
      
C parameters

      INTEGER TRIA(SZTRIA), MATPAR(SZN2MI), NUVP , IUPW
      INTEGER INONST, ILEV, KVERT, KMID, KCORVG, IELEMT
      DOUBLE PRECISION DY(NUVP), DU(NUVP), DRHS(NUVP), DTIM(5)
      DOUBLE PRECISION DD(NUVP), DAUX(NUVP)
      DOUBLE PRECISION OMGMIN, OMGMAX, OMEGA, THSTEP, DELU,DELP
      
      INTEGER IASMBL(SZASMI)
      DOUBLE PRECISION DASMBL(SZASMD)
      
      INTEGER KA1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,KST1,KM1,NA,NU,NP,
     *        KSCNPR,NBDMT

C local variables

      DOUBLE PRECISION AA1, AA2
      DOUBLE PRECISION SKV1, SKV2, DELU1, DELU2
      DOUBLE PRECISION DELT1, DELT2, DELT, DTIM2(5)
      INTEGER INDU1, INDU2, INDP, INDT1, INDT2, INDT, NVT, NEL, INEUM
      INTEGER KGRVEL

C     externals:

C     definition of finite elements

      EXTERNAL E030,E031,EM31,EM30

C     Stop the time we need for the calculation of the optimal
C     correction.

      CALL GTMAUX(DTIM2,DTIM,5,0)

C     At first transfer the starting addresses of the matrices
C     to local variables for faster access.
C
C     We assume here, that we only have one system and one pressure
C     matrix structure.

      KM1    = L(MATPAR(OLM))
      
      KST1   = L(MATPAR(OLST))

      KA1    = L(MATPAR(OLA1))
      KCOLA  = L(MATPAR(OLCLA1))
      KLDA   = L(MATPAR(OLLDA1))
      NA     = MATPAR(ONA1)

      KB1    = L(MATPAR(OLB1))
      KB2    = L(MATPAR(OLB2))
      KCOLB  = L(MATPAR(OLCLB1))
      KLDB   = L(MATPAR(OLLDB1))
      
      NU     = MATPAR(ONEQA)
      NP     = MATPAR(ONEQB)
      
      KVERT  = L(TRIA(OLVERT))
      KMID   = L(TRIA(OLMID))
      KCORVG = L(TRIA(OLCORVG))
      NEL    = TRIA(ONEL)
      NVT    = TRIA(ONVT)
      
      INEUM  = IASMBL(OINEUM)
      IELEMT = IASMBL(OIELEMT)
      
      IUPW   = IASMBL(OIUPW)

C     Use the shortcut nodal property array to decide which points are
C     on the boundary, which are Neumann boundary,...

      KSCNPR = L(TRIA(OTRIUD+1))
      NBDMT  = TRIA(OTRIUD)

C     If the ALE method is used, store the address of the grid velocoity
C     vector to KGRVEL. Otherwise store 1 into KGRVEL, which is used
C     as a dummy in the call to the upwinding routines to prevent
C     specifying a DWORK(.) address out of bounds.

      IF (IASMBL(OIALE).NE.0) THEN
        KGRVEL = L(IASMBL(OLGRVEL))
      ELSE
        KGRVEL = 1
      END IF
      
C     Clear the timing information array
      
      DTIM(1) = 0D0
      DTIM(2) = 0D0
      DTIM(3) = 0D0
      DTIM(4) = 0D0
      DTIM(5) = 0D0

C-----------------------------------------------------------------------

      IF (OMGMIN.EQ.OMGMAX) THEN
      
        IF (OMGMIN.LT.0D0) THEN

C         If OMGMIN=OMGMAX<0, we have nothing to do - no 
C         relative changes are calculated. Cancel here.

          OMEGA=ABS(OMGMIN)
          DELU=0D0
          DELP =0D0
          GOTO 99999
        ELSE

C         For OMGMIN=OMGMAX>0, there's no choice for OMEGA - 
C         prescribe it directly and skip its calculation.

          OMEGA=OMGMIN
          GOTO 999
          
        ENDIF
        
      ENDIF

C***********************************************************************

C     The first big part now is to calculate a new OMEGA parameter
C     with OMGMIN < OMEGA < OMGMAX.
C
C     The defect correction for a problem like T(u)u=f has the form
C
C           u^(l+1)  =  u^l  -  OMEGA * C * ( T(u^l)u^l - f )
C
C     with an appropriate preconditioner C (see below).
C     In our case, this iteration system can be written as:
C
C     (KU1)    (KU1)               ( [ KST1        KB1] (KU1)   (KF1) )
C     (KU2) := (KU2) - OMEGA * C * ( [       KST1  KB2] (KU2) - (KF2) )
C     (KP)     (KP)                ( [ KB1^T KB2^T  0 ] (KP)    (KFP) )
C
C                                  |----------------------------------|
C                                           = (KD1,KD2,KDP)^T
C                              |--------------------------------------|
C                                          = -Y = -(KU1,KU2,KP)
C
C     with KST1=KST1(KU1,KU2,KP) and Y being the solution from
C     the Oseen equation with 
C
C                      [ KST1        KB1 ]
C        C = T(u^l) =  [       KST1  KB2 ]
C                      [ KB1^T KB2^T  0  ]
C
C     The parameter OMEGA is calculated as the result of the 1D
C     minimization problem:
C
C       OMEGA = min_omega || T(u^l-omega*Y)*(u^l-omega*Y) - f ||_E
C
C               < T(u^l-omegaold*Y)Y , f - T(u^l-omegaold*Y)u^l >
C            ~= -------------------------------------------------
C                  < T(u^l-omegaold*Y)Y , T(u^l-omegaold*Y)Y >
C
C     when choosing omegaold=previous omega, which is a good choice
C     as one can see by linearization (see p. 187 (170), Turek's book).
C
C     Here, ||.||_E denotes the the Euclidian norm to the Euclidian 
C     scalar product <.,.>.

C=======================================================================
C *** Calculate on AUX1/2 the linearization point: AUX = UOLD+OMEGA*Y
C=======================================================================

C     We calculate DAUX=UOLD+OMEGA*Y instead of UOLD-OMEGA*Y, because
C     the "-"-sign is already incorporated in Y -- Y is calculated
C     as residuum "b-Ax" rather than as a defect "Ax-b"...

      AA1=1.0D0
      AA2=OMEGA
      CALL LCP1 (DY,DAUX,NUVP)

      CALL LLC1 (DU,DAUX,NUVP,AA1,AA2)
      
C=======================================================================
C *** First term of scalar product in the nominator
C
C     Calculate the new nonlinear block A at the
C     point AUX = u^l-omegaold*Y
C=======================================================================

C     Construct the linear part of the nonlinear matrix on the current
C     level. As there is no defect vector to be calculated (we
C     obtained it already by restriction), we use a routine here
C     that only builds the KA1-matrix.

      CALL GTMAUX (DTIM2,DTIM,2,0)
      CALL XMADF5(MATPAR,IASMBL(OIPRECA),THSTEP,INONST,ILEV)
      CALL GTMAUX (DTIM2,DTIM,2,1)

C     Up to now our system matrix has consisted only of linear
C     terms:
C
C         KA1 = [  M*I  +  THSTEP * (-nu * Laplace(.))  ]
C
C     but this is not enough... the nonlinear part is still missing.
C     So we have to add the following term to KA1:
C
C         THSTEP * u grad(.)
C
C     what will finally result in the system matrix
C
C         KA1 = [  M*I  +  THSTEP * (-nu * Laplace(.))  ] + THSTEP * u grad(.)
C             = [  M*I  +  THSTEP * N(u) ]
C
C     where u=KAUX in this case. Don't perform any correction of defect
C     vectors (IDEF=0), as we only want to have the matrix modified.

      CALL GTMAUX (DTIM2,DTIM,3,0)
      IF (IUPW.EQ.1) THEN
       CALL GUPWDX (DAUX(1),DAUX(1+NU),DAUX(1),DAUX(1+NU),
     *             1D0,0D0,DAUX(1),DAUX(1+NU),
     *             DAUX(1),DAUX(1+NU),DWORK(KA1),KWORK(KCOLA),
     *             KWORK(KLDA),TRIA,KWORK(KVERT),KWORK(KMID),
     *             DWORK(KCORVG),0,
     *             DASMBL(OUPSAM), DASMBL(ORE), THSTEP,
     *             IASMBL(OIALE),DWORK(KGRVEL))
      ELSE IF (IELEMT.EQ.0) THEN
C        CALL SUPGPX(NU,DAUX(1),DAUX(1+NU),DAUX(1),DAUX(1+NU),
C     *             1D0,0D0,DAUX(1),DAUX(1+NU),
C     *             DAUX(1),DAUX(1+NU),DWORK(KA1),NA,KWORK(KCOLA),
C     *             KWORK(KLDA),KWORK(KVERT),KWORK(KMID),
C     *             DWORK(KCORVG),TRIA,E031,IASMBL(OICUBN),
C     *             (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
C     *             IASMBL(OISTOK), DASMBL(OUPSAM), DASMBL(ORE),
C     *             0,1D0,THSTEP)
        CALL SUPAPX(NU,DAUX(1),DAUX(1+NU),DAUX(1),DAUX(1+NU),
     *             1D0,0D0,DAUX(1),DAUX(1+NU),
     *             DAUX(1),DAUX(1+NU),DWORK(KA1),NA,KWORK(KCOLA),
     *             KWORK(KLDA),KWORK(KVERT),KWORK(KMID),
     *             DWORK(KCORVG),TRIA,E031,IASMBL(OICUBN),
     *             (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
     *             IASMBL(OISTOK), DASMBL(OUPSAM), DASMBL(ORE),
     *             0,1D0,THSTEP,
     *             IASMBL(OIALE),DWORK(KGRVEL))
       ELSE IF (IELEMT.EQ.1) THEN
C        CALL SUPGPX(NU,DAUX(1),DAUX(1+NU),DAUX(1),DAUX(1+NU),
C     *             1D0,0D0,DAUX(1),DAUX(1+NU),
C     *             DAUX(1),DAUX(1+NU),DWORK(KA1),NA,KWORK(KCOLA),
C     *             KWORK(KLDA),KWORK(KVERT),KWORK(KMID),
C     *             DWORK(KCORVG),TRIA,E030,IASMBL(OICUBN),
C     *             (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
C     *             IASMBL(OISTOK), DASMBL(OUPSAM), DASMBL(ORE),
C     *             0,1D0,THSTEP)
        CALL SUPAPX(NU,DAUX(1),DAUX(1+NU),DAUX(1),DAUX(1+NU),
     *             1D0,0D0,DAUX(1),DAUX(1+NU),
     *             DAUX(1),DAUX(1+NU),DWORK(KA1),NA,KWORK(KCOLA),
     *             KWORK(KLDA),KWORK(KVERT),KWORK(KMID),
     *             DWORK(KCORVG),TRIA,E030,IASMBL(OICUBN),
     *             (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
     *             IASMBL(OISTOK), DASMBL(OUPSAM), DASMBL(ORE),
     *             0,1D0,THSTEP,
     *             IASMBL(OIALE),DWORK(KGRVEL))
      ELSE IF (IELEMT.EQ.2) THEN
C        CALL SUPGNX(NU,DAUX(1),DAUX(1+NU),DAUX(1),DAUX(1+NU),
C     *             1D0,0D0,DAUX(1),DAUX(1+NU),
C     *             DAUX(1),DAUX(1+NU),DWORK(KA1),NA,KWORK(KCOLA),
C     *             KWORK(KLDA),KWORK(KVERT),KWORK(KMID),
C     *             DWORK(KCORVG),TRIA,EM31,IASMBL(OICUBN),
C     *             (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
C     *             IASMBL(OISTOK), DASMBL(OUPSAM), DASMBL(ORE),
C     *             0,1D0,THSTEP)
        CALL SUPANX(NU,DAUX(1),DAUX(1+NU),DAUX(1),DAUX(1+NU),
     *             1D0,0D0,DAUX(1),DAUX(1+NU),
     *             DAUX(1),DAUX(1+NU),DWORK(KA1),NA,KWORK(KCOLA),
     *             KWORK(KLDA),KWORK(KVERT),KWORK(KMID),
     *             DWORK(KCORVG),TRIA,EM31,IASMBL(OICUBN),
     *             (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
     *             IASMBL(OISTOK), DASMBL(OUPSAM), DASMBL(ORE),
     *             0,1D0,THSTEP,
     *             IASMBL(OIALE),DWORK(KGRVEL))
      ELSE IF (IELEMT.EQ.3) THEN
C        CALL SUPGNX(NU,DAUX(1),DAUX(1+NU),DAUX(1),DAUX(1+NU),
C     *             1D0,0D0,DAUX(1),DAUX(1+NU),
C     *             DAUX(1),DAUX(1+NU),DWORK(KA1),NA,KWORK(KCOLA),
C     *             KWORK(KLDA),KWORK(KVERT),KWORK(KMID),
C     *             DWORK(KCORVG),TRIA,EM30,IASMBL(OICUBN),
C     *             (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
C     *             IASMBL(OISTOK), DASMBL(OUPSAM), DASMBL(ORE),
C     *             0,1D0,THSTEP)
        CALL SUPANX(NU,DAUX(1),DAUX(1+NU),DAUX(1),DAUX(1+NU),
     *             1D0,0D0,DAUX(1),DAUX(1+NU),
     *             DAUX(1),DAUX(1+NU),DWORK(KA1),NA,KWORK(KCOLA),
     *             KWORK(KLDA),KWORK(KVERT),KWORK(KMID),
     *             DWORK(KCORVG),TRIA,EM30,IASMBL(OICUBN),
     *             (IASMBL(OIPRECA).EQ.4),IASMBL(OIMASS),
     *             IASMBL(OISTOK), DASMBL(OUPSAM), DASMBL(ORE),
     *             0,1D0,THSTEP,
     *             IASMBL(OIALE),DWORK(KGRVEL))
      ENDIF
      CALL GTMAUX (DTIM2,DTIM,3,1)

C     Incorporate Dirichlet-boundary into the matrix. Replace all
C     rows corresponding to Dirichlet nodes by unit vectors.

      CALL GTMAUX (DTIM2,DTIM,4,0)
      CALL BDRYA (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *            KWORK(KSCNPR),NBDMT)
      CALL GTMAUX (DTIM2,DTIM,4,1)

C=======================================================================
C *** Second term of the scalar product in the nominator
C     Calculate the defect  D = F-T*UOLD,
C     store it in KD1,KD2,KDP, overwriting the old defect
C=======================================================================

      CALL LCP1 (DRHS,DD,NUVP)

      CALL GTMAUX (DTIM2,DTIM,2,0)
      
      CALL MTMUL2(DU,DD,-1D0,1D0,MATPAR)
     
C     Force the Dirichlet boundary components of the velocity
C     vector to zero:

      CALL BDRY0 (DD(1),DD(1+NU),KWORK(KSCNPR),NBDMT)

C=======================================================================
C     All terms in the fraction:
C     Calculate the value  AUX = T*Y
C=======================================================================

      CALL MTMUL2(DY,DAUX,1D0,0D0,MATPAR)
     
C     Force the Dirichlet boundary components of the velocity
C     defect vector to zero:

      CALL BDRY0 (DAUX(1),DAUX(1+NU),KWORK(KSCNPR),NBDMT)

      CALL GTMAUX (DTIM2,DTIM,2,1)

C=======================================================================
C *** Calculation of the fraction terms.
C     Calculate nominator:    SKV1:= (T*Y,D)   = (AUX,D)
C     Calculate denominator:  SKV2:= (T*Y,T*Y) = (AUX,AUX)
C=======================================================================

      CALL LSP1 (DAUX,DD,NUVP,SKV1)
      CALL LSP1 (DAUX,DAUX,NUVP,SKV2)

      IF (SKV2.LT. 1.0D-40) THEN
        WRITE(MTERM,*) 'ERROR in OPTCOR: SKVKV is nearly zero'
        STOP
      ENDIF

C     Ok, we have the nominator and the denominator. Divide them
C     by each other to calculate the new OMEGA.

      OMEGA=SKV1/SKV2
      
C     And make sure it's in the allowed range:

      IF (OMEGA.LT.ABS(OMGMIN)) OMEGA=ABS(OMGMIN)
      IF (OMEGA.GT.ABS(OMGMAX)) OMEGA=ABS(OMGMAX)

C     That's it, we have our new Omega.

999   CONTINUE

C=======================================================================
C *** Update of the solution vector.
C     Calculate  Unew  :=  u^l + Omega*Y  =  UOLD+OMEGA*Y 
C=======================================================================

C     Calculate maximum changes   DELU
C
C     Calculate the maximum norm of Y=(KU1, KU2, KP) and save the
C     it individually in DELU1, DELU2, DELP.
C     (INDxxx is the index of the component giving rise to the maximum)

      CALL  LLI1 (DY(1),NU,DELU1,INDU1)
      CALL  LLI1 (DY(1+NU),NU,DELU2,INDU2)
      CALL  LLI1 (DY(1+2*NU) ,NP,DELP ,INDP )

C=======================================================================
C *** Defect correction.
C     Update the solution   U := UOLD + OMEGA*Y
C=======================================================================

C     The following weights are exchanged in contrast to XOPTCN
C     because Y is now in DY, not in DU.

      AA1= OMEGA
      AA2= 1.D0
      CALL  LLC1 (DY,DU,NUVP,AA1,AA2)

C     This gives us our new solution DU, overwriting our
C     old solution vector - this one is now history...

C=======================================================================
C *** relative maximum changes   
C=======================================================================

C     This simply calculates some postprocessing values of the relative
C     change in the solution.
C
C     Maximum norm of solution vector:
C
C       DELT := || (KU1,KU2) ||_max
C
C                   || P ||_max
C       DELP := -------------------
C               || (KU1,KU2) ||_max
C
C
C     Relative change of solution vector:
C
C                || (Y1,Y2) ||_max
C       DELU := -------------------
C               || (KU1,KU2) ||_max

      CALL  LLI1 (DU(1),NU,DELT1,INDT1)
      CALL  LLI1 (DU(1+NU),NU,DELT2,INDT2)
      DELT=MAX(DELT1,DELT2)
      IF (ABS(DELT).LT.1D-8) DELT=1D0
      DELU=MAX(DELU1,DELU2)/DELT

      CALL  LLI1 (DU(1+2*NU),NP,DELT,INDT)
      IF (ABS(DELT).LT.1D-8) DELT=1D0
      DELP=DELP/DELT

C     Finally stop the time for the calculation of the correction

      CALL GTMAUX (DTIM2,DTIM,5,1)
      
C     The calculation of the LC time is a little bit tricky as
C     LC is usually done very fast. We calculate it by subtracting
C     the time of all non-LC operations from the total time.
      
      DTIM(1) = DTIM(5) - DTIM(2)-DTIM(3)-DTIM(4)

99999 END
