************************************************************************
* Invert a 4x4 system
*
* Outdated routine - LAPACK is faster :)
*
* In:
*   A    - array [1..4,1..4] of double
*          4x4-matrix; if IPAR=1: factorized matrix.
*   F    - array [1..4] of double
*          RHS of the system Ax=F
*   IPAR - =0: Factorize matrix A, don't calculate x
*          =1: Calculate x using the factorized matrix A
*
* Out:
*   If IPAR=0:
*      A    - Factorized matrix
*   If IPAR=1:
*      X    - array [1..4] of double
*             Solution x to Ax=F 
************************************************************************

      SUBROUTINE INVERT(A,F,X,IPAR)
      
      IMPLICIT NONE
      
C parameters
      
      INTEGER NDIM
      PARAMETER (NDIM=4)

      DOUBLE PRECISION A(NDIM,NDIM),B(NDIM,NDIM),F(NDIM),X(NDIM)
      INTEGER IPAR

C local variables

      INTEGER MERKX(NDIM),MERKY(NDIM)
      INTEGER IA,IB,IFEHL
      
C     IPAR=0: factorize A
      
      IF (IPAR.EQ.0) THEN
        CALL AUSTAU(NDIM,NDIM,A,B,MERKX,MERKY,IFEHL)
        DO IA=1,NDIM
          DO IB=1,NDIM
            A(IA,IB)=B(IA,IB)
          END DO
        END DO
      ENDIF

C     IPAR=1: Solve Ax=F

      IF (IPAR.EQ.1) THEN
        DO IA=1,NDIM
          X(IA)=0D0
          DO IB=1,NDIM
            X(IA)=X(IA)+A(IA,IB)*F(IB)
          END DO
        END DO
      ENDIF

      END

************************************************************************
      SUBROUTINE AUSTAU(NDIM,N,A,B,MERKX,MERKY,IFEHL)
      
      IMPLICIT NONE

C parameters

      INTEGER NDIM, N      
      DOUBLE PRECISION A(NDIM,N),B(NDIM,N)
      INTEGER MERKX(N),MERKY(N)
      
C local variables

      INTEGER I,J,K,L,IFEHL,IX,IY,INDX,INDY,M
      DOUBLE PRECISION HILF,PIVOT

      IFEHL=1
      DO I=1,N
        MERKX(I)=0
        MERKY(I)=0
        DO L=1,N
          B(I,L)=A(I,L)
        END DO
      END DO

      DO 400 I=1,N
      PIVOT=0D0
      DO 200 IX=1,N
      IF (MERKX(IX).NE.0) GOTO 200
      DO 180 IY=1,N
      IF (MERKY(IY).NE.0) GOTO 180
      IF (ABS(B(IX,IY)).LE.ABS(PIVOT)) GOTO 180
      PIVOT=B(IX,IY)
      INDX=IX
      INDY=IY
180   CONTINUE
200   CONTINUE
C
      IF (ABS(PIVOT).LE.0.0) GOTO 770
      MERKX(INDX)=INDY
      MERKY(INDY)=INDX
      B(INDX,INDY)=1D0/PIVOT
      DO 300 L=1,N
      IF (L.EQ.INDX) GOTO 300
      DO 280 M=1,N
      IF (M.EQ.INDY) GOTO 280
      B(L,M)=B(L,M)-B(L,INDY)*B(INDX,M)/PIVOT
280   CONTINUE
300   CONTINUE
C
      DO 390 IX=1,N
      IF (IX.NE.INDX) B(IX,INDY)=B(IX,INDY)/PIVOT
390   CONTINUE
      DO 400 IY=1,N
      IF (IY.NE.INDY) B(INDX,IY)=-B(INDX,IY)/PIVOT
400   CONTINUE
C
C
      DO 500 I=2,N
      IX=I-1
      IF (MERKX(IX).EQ.IX) GOTO 500
      DO 450 J=1,N
      IY=J
      IF (MERKX(IY).EQ.IX) GOTO 460
450   CONTINUE
460   DO 490 K=1,N
      HILF=B(IX,K)
      B(IX,K)=B(IY,K)
490   B(IY,K)=HILF
      MERKX(IY)=MERKX(IX)
500   MERKX(IX)=IX
      DO 600 I=2,N
      IX=I-1
      IF (MERKY(IX).EQ.IX) GOTO 600
      DO 550 J=1,N
      IY=J
      IF (MERKY(IY).EQ.IX) GOTO 560
550   CONTINUE
560   DO 590 K=1,N
      HILF=B(K,IX)
      B(K,IX)=B(K,IY)
590   B(K,IY)=HILF
      MERKY(IY)=MERKY(IX)
600   MERKY(IX)=IX
C
      IFEHL=0
770   RETURN
      END
