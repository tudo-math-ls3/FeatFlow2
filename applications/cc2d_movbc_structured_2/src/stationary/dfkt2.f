************************************************************************
* Computes the norms RESU, RESDIV
*
* This routine computes the "U-residual" RESU and the "divergence
* residual" RESDIV.
*
* In:
*   U1,
*   U2,
*   P      - Velocity/pressure vector of the solution of 
*            the linearised system    A * (U1,U2,P) = (F1,F2,FP)
*   D1,
*   D2,
*   DP     - Defect vector   (D1,D2,DP) = (F1,F2,FP) - A*(U1,U2,P)
*   F1,
*   F2,
*   FP     - RHS vector of the system   A * (U1,U2,P) = (F1,F2,FP)
*   NU     - length of the velocity vectors
*   NP     - length of the pressure vectors
* 
* Out:
*   RESU   - Norm of the velocity residual:
*                               || (D1,D2) ||_E
*               RESU   = -----------------------------
*                        max ( ||F1||_E , ||F2||_E )
*
*   RESDIV - Norm of the pressure residual
*                           || P ||_E
*               RESDIV = ----------------
*                        || (U1,U2) ||_E
************************************************************************

      SUBROUTINE RESDFK(U1,U2,P,D1,D2,DP,F1,F2,FP,NU,NP,RESU,RESDIV)
      
      IMPLICIT NONE
      
C parameters      

      INTEGER NU,NP
      DOUBLE PRECISION U1(*),U2(*),P(*),D1(*),D2(*),DP(*),F1(*),F2(*),
     *                 FP(*)
      DOUBLE PRECISION RESU, RESDIV

C local variables
      
      DOUBLE PRECISION RESF, RESF1, RESF2, RESU1, RESU2
      DOUBLE PRECISION DNORMU, DNRMU1, DNRMU2

C-----------------------------------------------------------------------
C     Compute the relative l2-norms  RESU,RESDIV
C-----------------------------------------------------------------------

C     RESF := max ( ||F1||_E , ||F2||_E )

      CALL LL21 (F1,NU,RESF1)
      CALL LL21 (F2,NU,RESF2)
      RESF=MAX(RESF1,RESF2)
      IF (ABS(RESF).LT.1D-8) RESF=1D0

C                   || (D1,D2) ||_E
C     RESU = -----------------------------
C            max ( ||F1||_E , ||F2||_E )

      CALL LL21 (D1,NU,RESU1)
      CALL LL21 (D2,NU,RESU2)
      RESU=SQRT(RESU1*RESU1+RESU2*RESU2)/RESF

C     DNORMU = || (U1,U2) ||_l2 

      CALL LL21 (U1,NU,DNRMU1)
      CALL LL21 (U2,NU,DNRMU2)
      DNORMU=SQRT(DNRMU1*DNRMU1+DNRMU2*DNRMU2)
      IF (ABS(DNORMU).LT.1D-8) DNORMU=1D0

C                 || P ||_E
C     RESDIV = ----------------
C              || (U1,U2) ||_E

      CALL LL21 (DP,NP,RESDIV)
      RESDIV=RESDIV/DNORMU
C
      END

************************************************************************
* Defect-RHS and matrix construction routine 4
*
* Extended calling convention
*
* Build the linear part of the defect vector for the nonlinear 
* iteration. Combinate Stokes- and Mass-matrix to build the
* linear part of the nonlinear system matrix.
*
* This routine is called in NSDEF2 to assemble mass-/system
* matrices for the nonlinear iteration. It has the same duty as
* XMADF2 but uses a different calling convention. 
*
* The shortcut nodal property array is used to decide on which
* DOF's are boundary and which not!
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Triangulation structure on current level
*   MATPAR - array [1..SZN2VI] of integer
*            TNS2DMatrixParams structure containing all information about
*            system-/mass-/Laplace matrices on the current level
*   DU     - array [1..2*NEQA+NEQB] of double
*            This is the the X-,Y- and pressure-component of the current 
*            solution vector u^n of the equation.
*            Must correspond to MATPAR.
*   DD     - array [1..2*NEQA+NEQB] of double
*            Vector that contains the current RHS-vector of the linear 
*            system for velocity and pressure. This will be overwritten
*            by the resulting defect vector, if IDEF=1.
*            Must correspond to MATPAR.
*   THSTEP - double. Theta-scheme identifier.
*            =0: Forward Euler, =1: Backward Euler, =0.5: Crank Nicolson
*            For Fractional-Step theta scheme this is set to different
*            values, depending on the current step in the scheme.
*            For stationary simulations this parameter must be set
*            to 1 to include the full Laplacian matrix into the
*            system matrix:
*                    [alpha*M + THETA*KST1 ] 
*   IALPHA - Whether the simulation is stationary or not.
*            As this routine generally builds
*                    [alpha*M + Theta_1*nu*k*L ] u
*            (the linear part of S(u)), this parameter decides whether
*            ALPHA=0 or ALPHA=1, so whether the mass 
*            matrix is added or not.   
*   IDEF   - Whether to update the defect vector according to
*             DEF := DEF - [ ALPHA*M*I - THSTEP*(-nu*Laplace(u^n)) ] u 
*            = 0: don't update defect vector
*            <>0: update defect vector
*   IPRECA - Defines how to set up the system matrix.
*            = 0,
*            = 1: Copy a precalculated Stokes matrix KST to KA
*                 and add mass matrix if desired
*            = 2,
*            = 3: Read precalculated Stokes matrix for level ILEV from
*                 disc and add mass matrix if desired
*            = 4: Clear the System matrix. In a stationary simulation
*                 (IALPHA=1) with lumped mass matrix (IMASS=0), 
*                 add the mass matrix to the system matrix.
*   INEUM  - Whether there are Neumann boundary components in our
*            problem or not.
*   ILEV   - Current level. Only used if IPRECA=3 to read the correct
*            matrix into the memory - otherwise ignored.
*
* Out:
*   - The system matrix (identified by MATPAR(ILEV).LA1) is filled with 
*     a linear approximation to the nonlinear system matrix
*   - DD   - Will receive the defect vector, calculated
*            with the RHS, the solution and the linear part of the
*            nonlinear system matrix:
*                  DD = DD - [alpha*M*u-THSTEP*{-nu*Laplace(u)}]*u
************************************************************************

      SUBROUTINE XMADF4(TRIA,MATPAR,DU,DD,IDEF,IPRECA,
     *                  INEUM,THSTEP,IALPHA,ILEV)
     
************************************************************************

      IMPLICIT NONE

C main COMMON blocks

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'stria.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'

C parameters

      INTEGER TRIA(SZTRIA),MATPAR(SZN2MI)
      INTEGER INEUM,IALPHA,IDEF,IPRECA,ILEV
      DOUBLE PRECISION DU(*),DD(*),THSTEP

C local variables

      INTEGER KM1,KST1,KA1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB
      INTEGER NA,KSCNPR,NBDMT,NU
      INTEGER INU, INUA

C constants for filenames

      CHARACTER CFILST*12,CFILM*12,CFILE*12
      DIMENSION CFILST(NNLEV),CFILM(NNLEV)
      DATA CFILST/'ns/ST1     ','ns/ST2     ','ns/ST3     ',
     *            'ns/ST4     ','ns/ST5     ','ns/ST6     ',
     *            'ns/ST7     ','ns/ST8     ','ns/ST9     '/
      DATA CFILM /'ns/MA1     ','ns/MA2     ','ns/MA3     ',
     *            'ns/MA4     ','ns/MA5     ','ns/MA6     ',
     *            'ns/MA7     ','ns/MA8     ','ns/MA9     '/
      SAVE CFILST,CFILM,CFILE

C     At first transfer the starting addresses of the matrices
C     to local variables for faster access.
C
C     We assume here, that we only have one system and one pressure
C     matrix structure.

      KM1    = L(MATPAR(OLM))
      KA1    = L(MATPAR(OLA1))
      KCOLA  = L(MATPAR(OLCLA1))
      KLDA   = L(MATPAR(OLLDA1))
      KB1    = L(MATPAR(OLB1))
      KB2    = L(MATPAR(OLB2))
      KCOLB  = L(MATPAR(OLCLB1))
      KLDB   = L(MATPAR(OLLDB1))
      NA     = MATPAR(ONA1)

C     Use the shortcut nodal property array to decide which points are
C     on the boundary, which are Neumann boundary,...
      
      NBDMT   = TRIA(OTRIUD)
      KSCNPR  = L(TRIA(OTRIUD+1))
      
C     We known the number of velocity equations by the size of A

      NU = MATPAR(ONEQA)

C     Standard matrix assembling branch:

      IF ((IPRECA.EQ.0).OR.(IPRECA.EQ.1)) THEN

C       We have a Stokes matrix

        KST1   = L(MATPAR(OLST))

C       Decide whether or not to include the mass matrix (for
C       instationary simulation) into the system matrix:

        IF (IALPHA.EQ.1) THEN

C         We have to build the linear part of the nonlinear system 
C         matrix with the help of the Stokes- and mass-matrices:
C
C         KA1 = [  M*I  +  THSTEP * N(u)  ]
C             = [  M*I  +  THSTEP * (-nu * Laplace(.))  +  u * grad(.)  ]
C                  ^^^              ^^^^^^^^^^^^^^^^^     ^^^^^^^^^^^
C                 ->KM1                ->KST1             -> ignored here
C
C         Check if we have a real or a lumped mass matrix.

          IF (MATPAR(OIMALMP).EQ.1) THEN

C           We have a real mass matrix in the structure of the 
C           system matrix.
C           Now check the Theta-parameter. If it's 0, we can skip the 
C           linear combination with the Stokes-matrix. Otherwise
C           simply add the matrices:
C
C                          KM1 + THSTEP*KST1

            IF (THSTEP.NE.0D0) THEN
              CALL LLC1(DWORK(KST1),DWORK(KA1),NA,THSTEP,0D0)
              CALL LLC1(DWORK(KM1 ),DWORK(KA1),NA,1D0,1D0)
            ELSE
              CALL LCP1(DWORK(KM1),DWORK(KA1),NA)
            ENDIF
            
          ELSE
          
C           Using lumped mass matrices we have only to tackle
C           the diagonal when adding it to the system matrix.
C           Again add the matrices:
C
C                          KM1 + THSTEP*KST1

            IF (THSTEP.NE.0D0) THEN
              CALL LLC1(DWORK(KST1),DWORK(KA1),NA,THSTEP,0D0)
              DO INU=1,NU
                INUA=KWORK(KLDA+INU-1)-1
                DWORK(KA1+INUA)=DWORK(KA1+INUA)+DWORK(KM1+INU-1)
              END DO
            ELSE
              CALL LCL1(DWORK(KA1),NA)
              DO INU=1,NU
                INUA=KWORK(KLDA+INU-1)-1
                DWORK(KA1+INUA)=DWORK(KM1+INU-1)
              END DO
            ENDIF
          ENDIF
        ELSE
        
C         The case of a stationary simulation is much easier. 
C         Remember, we have to solve a system of the form
C
C             S(u_h)u_h  +  k B p_h  =  g_h,   B^T u_h = 0
C
C         where
C
C             S(u) = [alpha*M + Theta_1*nu*k*L + Theta_2*k*K(u)] u
C         
C         is the matrix whose linear parts we are building here.
C         In the stationary approach, there is alpha=0, so adding
C         the mass matrix can be dropped completely!
C
C         Therefore we only have to build
C
C         KA1 = [  THSTEP * (-nu * Laplace(.))  +  u * grad(.)  ]
C                           ^^^^^^^^^^^^^^^^^     ^^^^^^^^^^^
C                              ->KST1             -> ignored here
C
C         what can simply be done by copying KST1 to KA1 with the 
C         correct scaling factor.
        
          CALL LLC1(DWORK(KST1),DWORK(KA1),NA,THSTEP,0D0)
          
        ENDIF

C       Multiply the solution vector with the built matrix, and 
C       subtract it from the RHS-vector (the current setting of KDx)
C       to build the defect vector.  
C
C         DEF := DEF - [ ALPHA*M*I - THSTEP*(-nu*Laplace(u^n)) ] u

        IF (IDEF.NE.0) THEN
          CALL  LAX17 (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DU(1),DD(1),-1D0,1D0)
          CALL  LAX17 (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DU(1+NU),DD(1+NU),-1D0,1D0)
        ENDIF

      ENDIF

C=======================================================================

C     Check if the matrices are so big that they have to be read in
C     from a file. The duty of this branch is the same as above.

      IF ((IPRECA.EQ.2).OR.(IPRECA.EQ.3)) THEN

C       Read data about level ILEV into memory.

        IF (IALPHA.EQ.1) THEN
          IF (THSTEP.NE.0D0) THEN
            CALL  OF0 (59,CFILST(ILEV),0)
            CFILE='STMAT '
            CALL  ORA1 (DWORK(KA1),CFILE,59,0)
            REWIND(59)
            CLOSE (59)
            IF (IER.NE.0) RETURN
            CALL LLC1(DWORK(KA1),DWORK(KA1),NA,THSTEP,0D0)
          ELSE
            CALL LCL1(DWORK(KA1),NA)
          ENDIF

          IF (MATPAR(OIMALMP).EQ.1) THEN
            CALL  OF0 (59,CFILM(ILEV),0)
            CFILE='MASMAT'
            CALL  ORALC1 (DWORK(KA1),1D0,CFILE,59,0)
            REWIND(59)
            CLOSE (59)
            IF (IER.NE.0) RETURN
          ELSE
            IF (THSTEP.NE.0D0) THEN
              DO INU=1,NU
                INUA=KWORK(KLDA+INU-1)-1
                DWORK(KA1+INUA)=DWORK(KA1+INUA)+DWORK(KM1+INU-1)
              END DO
            ELSE
              DO INU=1,NU
                INUA=KWORK(KLDA+INU-1)-1
                DWORK(KA1+INUA)=DWORK(KM1+INU-1)
              END DO
            ENDIF
          ENDIF
        ELSE
          CALL  OF0 (59,CFILST(ILEV),0)
          CFILE='STMAT '
          CALL  ORA1 (DWORK(KA1),CFILE,59,0)
          REWIND(59)
          CLOSE (59)
          IF (IER.NE.0) RETURN
          CALL LLC1(DWORK(KA1),DWORK(KA1),NA,THSTEP,0D0)
        ENDIF

        IF (IDEF.NE.0) THEN
          CALL  LAX17 (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DU(1),DD(1),-1D0,1D0)
          CALL  LAX17 (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DU(1+NU),DD(1+NU),-1D0,1D0)
        ENDIF

      ENDIF

C=======================================================================

C     In this branch we don't rely on any precalculated information
C     in KST1. If we have a lumped mass matrix, we can build a part of
C     the RHS-vector and of KA1 - otherwise simply clear KA1 as is
C     will be calculated later by the caller.

      IF (IPRECA.EQ.4) THEN

        IF (IALPHA.EQ.1) THEN
        
          IF (MATPAR(OIMALMP).EQ.0) THEN
          
            IF (IDEF.NE.0) THEN
            
C             Treat both components separately instead of 
C             simultaneously in one loop. That is a little faster
C             as the compiler can unroll it...
            
              DO INU=1,NU
                DD(INU)= DD(INU)-DWORK(KM1+INU-1)*DU(INU)
              END DO
              
              DO INU=1,NU
                DD(NU+INU)= DD(NU+INU)-DWORK(KM1+INU-1)*DU(NU+INU)
              END DO
              
            ENDIF ! (IDEF.NE.0)

            CALL LCL1(DWORK(KA1),NA)

            DO INU=1,NU
              INUA=KWORK(KLDA+INU-1)-1
              DWORK(KA1+INUA)=DWORK(KM1+INU-1)
            END DO
            
          ELSE
          
            CALL LCL1(DWORK(KA1),NA)
            
          ENDIF ! (IMASS.EQ.0)
          
        ELSE
        
          CALL LCL1(DWORK(KA1),NA)
          
        ENDIF ! (IALPHA.EQ.1)

      ENDIF ! (IPRECA.EQ.4)

C     That's it for the velocity part of the system matrix
C
C=======================================================================
C
C     Now turn over to handle the pressure part. Up to now we
C     constructed the defect vector KD1/KD2 only using the velocity
C     part - but of course that's not enough. Also the pressure-part
C     has to be included into the defect calculation.
C
C     So now update the defect vector by
C
C         KD1 = KD1 - KB1*KP
C
C     to get the correct defect for the velocity. Remember, the
C     whole thing we are calculating here is:
C
C        (def1)   (d1)   [ KST1        KB1 ] ( u1 ) 
C        (def2) = (d2) - [       KST1  KB2 ] ( u2 ) 
C        (defp)   (dp)   [ KB1^T KB2^T  0  ] ( p  ) 
C
C     Standard handling: IPRECB<>2

C      IF (IPRECB.NE.2) THEN
      
C       Subtract KB1*p from KD1 to get the correct defect for U
        
        CALL  LAX19 (DWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,
     *               DU(1+2*NU),DD(1),-1D0,1D0)
     
C       Subtract KB2*p from KD2 to get the correct defect for V
        
        CALL  LAX19 (DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,
     *               DU(1+2*NU),DD(1+NU),-1D0,1D0)
        

C      ELSE
C
C        Special treatment: IPRECB=2=elementwise application.
C        Might not be implemented completely...
C      
C        CALL BMUL1 (KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
C     *              DWORK(L(LCORVG)),DU(1+2*NU),DD(1),DD(1+NU)
C     *              NEL,NVT,NMT,-1D0,1D0)
C     
C      ENDIF

C=======================================================================

C     Implement Dirichlet-values into the defect.
C     On all Dirichlet-nodes, the defect must be exactly 0.

      IF (IDEF.NE.0) CALL BDRY0 (DD(1),DD(1+NU),KWORK(KSCNPR),NBDMT)

C=======================================================================

C     Finally calculate the defect vector of the pressure by subtracting
C     the velocity components (multiplied with KB) from KDP.
C
C       KDP = KDP  -  KB1^T KU1  -  KB2^T KU2

C      IF (IPRECB.NE.2) THEN

       CALL  LTX19 (DWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,
     *              MATPAR(ONEQB),DU(1),DD(1+2*NU),-1D0,1D0)
       CALL  LTX19 (DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,
     *              MATPAR(ONEQB),DU(1+NU),DD(1+2*NU),-1D0,1D0)

C      ELSE

C       Special treatment: IPRECB=2=elementwise application.
C       Might not be implemented completely...

C       WRITE(6,*) 'BTMUL ???'
C       CALL BTMUL1(KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
C     *             DWORK(L(LCORVG)),DU(1),DU(1+NU),DD(1+2*NU),
C     *             NEL,NVT,NMT,-1D0)
C      ENDIF

99998 END

************************************************************************
* Matrix construction routine 6
*
* Extended calling convention
*
* Combinate Stokes- and Mass-matrix to build the
* linear part of the nonlinear system matrix.
*
* This routine is called in NGSTP2 to assemble mass-/system
* matrices for the nonlinear iteration on all lower levels than NLMAX.
*
* It's basically the same routine as XMADF3, except for that it uses
* a different calling convention.
*
* In:
*   MATPAR - array [1..SZN2VI] of integer
*            TNS2DMatrixParams structure containing all information about
*            system-/mass-/Laplace matrices on the current level
*   IPRECA - Defines how to set up the system matrix.
*            = 0,
*            = 1: Copy a precalculated Stokes matrix KST to KA
*                 and add mass matrix if desired
*            = 2,
*            = 3: Read precalculated Stokes matrix for level ILEV from
*                 disc and add mass matrix if desired
*            = 4: Clear the System matrix. In a stationary simulation
*                 (IALPHA=1) with lumped mass matrix (IMASS=0), 
*                 add the mass matrix to the system matrix.
*   THSTEP - double. Theta-scheme identifier.
*            =0: Forward Euler, =1: Backward Euler, =0.5: Crank Nicolson
*            For Fractional-Step theta scheme this is set to different
*            values, depending on the current step in the scheme.
*            For stationary simulations this parameter must be set
*            to 1 to include the full Laplacian matrix into the
*            system matrix:
*                    [alpha*M + THETA*KST1 ] 
*   IALPHA - Whether the simulation is stationary or not.
*            As this routine generally builds
*                    [alpha*M + Theta_1*nu*k*L ] u
*            (the linear part of S(u)), this parameter decides whether
*            ALPHA=0 or ALPHA=1, so whether the mass 
*            matrix is added or not.   
*   ILEV   - Current level. Only used if IPRECA=3 to read the correct
*            matrix into the memory - otherwise ignored.
*
* Out:
*   - The system matrix (identified by MATPAR(ILEV).LA1) is filled with 
*     a linear approximation to the nonlinear system matrix
************************************************************************

      SUBROUTINE XMADF5(MATPAR,IPRECA,THSTEP,IALPHA,ILEV)
      
      IMPLICIT NONE

C main COMMON blocks

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'stria.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'

C parameters

      INTEGER MATPAR(SZN2MI)
      INTEGER IALPHA,IPRECA,ILEV
      DOUBLE PRECISION THSTEP

C local variables

      INTEGER KM1,KST1,KA1,KCOLA,KLDA,NA,NU
      INTEGER INU, INUA

C constants

      CHARACTER CFILST*12,CFILM*12,CFILE*12
      DIMENSION CFILST(NNLEV),CFILM(NNLEV)
      DATA CFILST/'ns/ST1     ','ns/ST2     ','ns/ST3     ',
     *            'ns/ST4     ','ns/ST5     ','ns/ST6     ',
     *            'ns/ST7     ','ns/ST8     ','ns/ST9     '/
      DATA CFILM /'ns/MA1     ','ns/MA2     ','ns/MA3     ',
     *            'ns/MA4     ','ns/MA5     ','ns/MA6     ',
     *            'ns/MA7     ','ns/MA8     ','ns/MA9     '/
      SAVE CFILST, CFILM, CFILE

C     At first transfer the starting addresses of the matrices
C     to local variables for faster access.
C
C     We assume here, that we only have one system and one pressure
C     matrix structure.

      KM1    = L(MATPAR(OLM))
      KA1    = L(MATPAR(OLA1))
      KCOLA  = L(MATPAR(OLCLA1))
      KLDA   = L(MATPAR(OLLDA1))
      NA     = MATPAR(ONA1)
      NU     = MATPAR(ONEQA)
      
C     Standard matrix assembling branch:

      IF ((IPRECA.EQ.0).OR.(IPRECA.EQ.1)) THEN

C       We have a Stokes matrix

        KST1   = L(MATPAR(OLST))

C       Decide whether or not to include the mass matrix (for
C       instationary simulation) into the system matrix:

        IF (IALPHA.EQ.1) THEN

C         We have to build the linear part of the nonlinear system 
C         matrix with the help of the Stokes- and mass-matrices:
C
C         KA1 = [  M*I  +  THSTEP * N(u)  ]
C             = [  M*I  +  THSTEP * (-nu * Laplace(.))  +  u * grad(.)  ]
C                  ^^^              ^^^^^^^^^^^^^^^^^     ^^^^^^^^^^^
C                 ->KM1                ->KST1             -> ignored here
C
C         Check if we have a real or a lumped mass matrix.

          IF (MATPAR(OIMALMP).EQ.1) THEN

C           We have a real mass matrix in the structure of the 
C           system matrix.
C           Now check the Theta-parameter. If it's 0, we can skip the 
C           linear combination with the Stokes-matrix. Otherwise
C           simply add the matrices:
C
C                          KM1 + THSTEP*KST1

            IF (THSTEP.NE.0D0) THEN
              CALL LLC1(DWORK(KST1),DWORK(KA1),NA,THSTEP,0D0)
              CALL LLC1(DWORK(KM1),DWORK(KA1),NA,1D0,1D0)
            ELSE
              CALL LCP1(DWORK(KM1),DWORK(KA1),NA)
            ENDIF
          ELSE

C           Using lumped mass matrices we have only to tackle
C           the diagonal when adding it to the system matrix.
C           Again add the matrices:
C
C                          KM1 + THSTEP*KST1

            IF (THSTEP.NE.0D0) THEN
              CALL LLC1(DWORK(KST1),DWORK(KA1),NA,THSTEP,0D0)
              DO INU=1,NU
                INUA=KWORK(KLDA+INU-1)-1
                DWORK(KA1+INUA)=DWORK(KA1+INUA)+DWORK(KM1+INU-1)
              END DO
            ELSE
              CALL LCL1(DWORK(KA1),NA)
              DO INU=1,NU
                INUA=KWORK(KLDA+INU-1)-1
                DWORK(KA1+INUA)=DWORK(KM1+INU-1)
              END DO
            ENDIF
          ENDIF
        ELSE
        
C         The case of a stationary simulation is much easier. 
C         Remember, we have to solve a system of the form
C
C             S(u_h)u_h  +  k B p_h  =  g_h,   B^T u_h = 0
C
C         where
C
C             S(u) = [alpha*M + Theta_1*nu*k*L + Theta_2*k*K(u)] u
C         
C         is the matrix whose linear parts we are building here.
C         In the stationary approach, there is alpha=0, so adding
C         the mass matrix can be dropped completely!
C
C         Therefore we only have to build
C
C         KA1 = [  THSTEP * (-nu * Laplace(.))  +  u * grad(.)  ]
C                           ^^^^^^^^^^^^^^^^^     ^^^^^^^^^^^
C                              ->KST1             -> ignored here
C
C         what can simply be done by copying KST1 to KA1 with the 
C         correct scaling factor.

          CALL LLC1(DWORK(KST1),DWORK(KA1),NA,THSTEP,0D0)
          
        ENDIF

      ENDIF

C=======================================================================

C     Check if the matrices are so big that they have to be read in
C     from a file. The duty of this branch is the same as above.

      IF ((IPRECA.EQ.2).OR.(IPRECA.EQ.3)) THEN

        IF (THSTEP.NE.0D0) THEN
          CALL  OF0 (59,CFILST(ILEV),0)
          CFILE='STMAT '
          CALL  ORA1 (DWORK(KA1),CFILE,59,0)
          REWIND(59)
          CLOSE (59)
          IF (IER.NE.0) RETURN
          CALL LLC1(DWORK(KA1),DWORK(KA1),NA,THSTEP,0D0)
        ELSE
          CALL LCL1(DWORK(KA1),NA)
        ENDIF

        IF (MATPAR(OIMALMP).EQ.1) THEN
          CALL  OF0 (59,CFILM(ILEV),0)
          CFILE='MASMAT'
          CALL  ORALC1 (DWORK(KA1),1D0,CFILE,59,0)
          REWIND(59)
          CLOSE (59)
          IF (IER.NE.0) RETURN
        ELSE
          IF (THSTEP.NE.0D0) THEN
            DO INU=1,NU
              INUA=KWORK(KLDA+INU-1)-1
              DWORK(KA1+INUA)=DWORK(KA1+INUA)+DWORK(KM1+INU-1)
            END DO
          ELSE
            DO INU=1,NU
              INUA=KWORK(KLDA+INU-1)-1
              DWORK(KA1+INUA)=DWORK(KM1+INU-1)
            END DO
          ENDIF
        ENDIF

      ENDIF
C=======================================================================

C     In this branch we don't rely on any precalculated information
C     in KST1. If we have a lumped mass matrix, we can build a part of
C     the RHS-vector and of KA1 - otherwise simply clear KA1 as is
C     will be calculated later by the caller.

      IF (IPRECA.EQ.4) THEN

        IF (MATPAR(OIMALMP).EQ.0) THEN
          CALL LCL1(DWORK(KA1),NA)
          DO INU=1,NU
            INUA=KWORK(KLDA+INU-1)-1
            DWORK(KA1+INUA)=DWORK(KM1+INU-1)
          END DO
        ELSE
          CALL LCL1(DWORK(KA1),NA)
        ENDIF

      ENDIF

C=======================================================================

99998 END
