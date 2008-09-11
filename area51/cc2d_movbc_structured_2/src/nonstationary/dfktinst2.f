************************************************************************
* RHS and matrix construction routine 6
*
* Recombinate the linear part of the RHS for the Theta-scheme.
* for one macro-step. Combinate Stokes- and Mass-matrix to build the
* linear part of the nonlinear system matrix.
*
* This routine is called in MGSTP2 to assemble mass-/system
* matrices for the nonlinear iteration. It has the same duty as
* XMADF1 but uses a different calling convention. 
*
* In:
*   DST1   - array [1..?] of double.
*            This "Stokes"-matrix represents the linear / Laplacian
*            part of the nonlinear matrix
*              N(u) = -nu * Laplace (.)  +  u * grad(.)
*            This precalculated matrix is used as a bases for N that
*            way that other terms are added to it to build
*            the nonlinear matrix.
*   CM1    - array [1..?] of double.
*            This is the precalculated mass-matrix of the system
*   DA1    - array [1..NA] of double.
*            Will be filled with data about the linear part of the
*            nonlinear system matrix for the Theta-scheme:
*                     [ M*I + THSTEP * (-nu*Laplace) ] 
*   KCOLA,
*   KLDA   - Structure of the system matrix
*   DF     - array [1..2*NU] of double
*            X- any Y-component of the RHS of the equation
*   DU     - array [1..2*NU] of double
*            X- and Y-component of the current solution vector u^n 
*            of the equation
*   THSTEP - double. Theta-scheme identifier.
*            =0: Forward Euler, =1: Backward Euler, =0.5: Crank Nicolson
*            For Fractional-Step theta scheme this is set to different
*            values, depending on the current step in the scheme.
*   IPRECA - Determines whether the linear part of the system matrix/RHS
*            is constructed with precalculated data.
*            If IPRECA=4, the Stokes-part -nu*Laplace(u) is not added
*             to KA1. If there is a lumped mass matrix, that part of KA1
*             is build and a part of the RHS-vector is constructed:
*               RHS = M*u + KF,
*               KA1 = M*I
*             If the mass matrix is not lumped, KA1 is filled with 0 and
*             the RHS-vector KFx is not touched.
*            If IPRECA=2,3, KM1/KST1 is read from disc before
*             reconstruction.   
*   ILEV   - Current level. Only used if IPRECA=3 to read the correct
*            matrix into the memory - otherwise ignored.
*
* Out:
*   DA1 - array for the nonlinear system matrix, 
*         filled with data about the linear part
*   DF  - By linear combination with the system matrix, these vectors
*         are updated now to represent the "linear part" of the right 
*         hand side of the Theta-scheme:
*               DF = [M*u-THSTEP*{-nu*Laplace(u)}] + DF
************************************************************************

      SUBROUTINE XMADF6(TRIA,MATPAR,DF,DU,
     *                  THSTEP,IPRECA,ILEV)

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

      INTEGER TRIA(SZTRIA),MATPAR(*)
      DOUBLE PRECISION DF(*),DU(*)
      DOUBLE PRECISION THSTEP
      INTEGER IPRECA,ILEV

C constants

      CHARACTER CFILST*12,CFILM*12,CFILE*12
      DIMENSION CFILST(NNLEV),CFILM(NNLEV)
      DATA CFILST/'ns/ST1     ','ns/ST2     ','ns/ST3     ',
     *            'ns/ST4     ','ns/ST5     ','ns/ST6     ',
     *            'ns/ST7     ','ns/ST8     ','ns/ST9     '/
      DATA CFILM /'ns/MA1     ','ns/MA2     ','ns/MA3     ',
     *            'ns/MA4     ','ns/MA5     ','ns/MA6     ',
     *            'ns/MA7     ','ns/MA8     ','ns/MA9     '/
      SAVE CFILST,CFILM,CFILE

C local variables

      INTEGER KM1,KST1,KA1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB
      INTEGER NA,KMBD,NMBD,NU
      INTEGER INU, INUA

C     At first transfer the starting addresses of the matrices
C     to local variables for faster access.
C
C     We assume here, that we only have one system and one pressure
C     matrix structure.

      KM1    = L(MATPAR(OLM))
      KA1    = L(MATPAR(OLA1))
      KST1   = L(MATPAR(OLST))
      KCOLA  = L(MATPAR(OLCLA1))
      KLDA   = L(MATPAR(OLLDA1))
      KB1    = L(MATPAR(OLB1))
      KB2    = L(MATPAR(OLB2))
      KCOLB  = L(MATPAR(OLCLB1))
      KLDB   = L(MATPAR(OLLDB1))
      NA     = MATPAR(ONA1)
      
      KMBD   = L(TRIA(OLMBD))
      NMBD   = TRIA(ONVBD)

C     We known the number of velocity equations by the size of A

      NU = MATPAR(ONEQA)

C     Standard matrix assembling branch:

      IF ((IPRECA.EQ.0).OR.(IPRECA.EQ.1)) THEN
      
C         We have to build the linear part of the nonlinear system 
C         matrix with the help of the Stokes- and mass-matrices:
C
C       KA1 = [  M*I  +  THSTEP * N(u)  ]
C           = [  M*I  +  THSTEP * (-nu * Laplace(.))  +  u * grad(.)  ]
C                ^^^              ^^^^^^^^^^^^^^^^^     ^^^^^^^^^^^
C               ->KM1                ->KST1             -> ignored here
C
C       Check if we have a real or a lumped mass matrix.

        IF (MATPAR(OIMALMP).EQ.1) THEN
        
C         We have a real mass matrix in the structure of the 
C         system matrix.
C         Now check the Theta-parameter. If it's 0, we can skip the 
C         linear combination with the Stokes-matrix. Otherwise
C         simply add the matrices:
C
C                        KM1 + THSTEP*KST1
        
          IF (THSTEP.NE.0D0) THEN
            CALL LLC1(DWORK(KST1),DWORK(KA1),NA,THSTEP,0D0)
            CALL LLC1(DWORK(KM1 ),DWORK(KA1),NA,1D0,1D0)
          ELSE
            CALL LCP1(DWORK(KM1),DWORK(KA1),NA)
          ENDIF
          
        ELSE
        
C         Using lumped mass matrices we have only to tackle
C         the diagonal when adding it to the system matrix.
C         Again add the matrices:
C
C                        KM1 + THSTEP*KST1
        
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

C       Multiply the current solution- and RHS-vector to build
C       the "linear" part of the RHS-vector of the Theta-scheme:
C
C         RHS := [ M*I - THSTEP*(-nu*Laplace(u^n)) ] u^n  +  RHS

        CALL LAX17 (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *              DU(1),DF(1),1D0,1D0)
        CALL LAX17 (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *              DU(1+NU),DF(1+NU),1D0,1D0)

      ENDIF

C=======================================================================

C     Check if the matrices are so big that they have to be read in
C     from a file. The duty of this branch is the same as above.

      IF ((IPRECA.EQ.2).OR.(IPRECA.EQ.3)) THEN

C       Read data about level ILEV into memory.

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

        CALL  LAX17 (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *               DU(1),DF(1),1D0,1D0)
        CALL  LAX17 (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *               DU(1+NU),DF(1+NU),1D0,1D0)

      ENDIF

C=======================================================================

C     In this branch we don't rely on any precalculated information
C     in KST1. If we have a lumped mass matrix, we can build a part of
C     the RHS-vector and of KA1 - otherwise simply clear KA1 as is
C     will be calculated later by the caller.

      IF (IPRECA.EQ.4) THEN

        IF (MATPAR(OIMALMP).EQ.0) THEN

C         Treat both components separately instead of 
C         simultaneously in one loop. That is a little faster
C         as the compiler can unroll it...

          DO INU=1,NU
            DF(INU)   = DF(INU)    + DWORK(KM1+INU-1)*DU(INU)
          END DO

          DO INU=1,NU
            DF(NU+INU)= DF(NU+INU) + DWORK(KM1+INU-1)*DU(NU+INU)
          END DO

          CALL LCL1(DWORK(KA1),NA)

          DO INU=1,NU
            INUA            = KWORK(KLDA+INU-1)-1
            DWORK(KA1+INUA) = DWORK(KM1+INU-1)
          END DO
          
        ELSE
        
          CALL LCL1(DWORK(KA1),NA)
          
        ENDIF

      ENDIF

99998 END
