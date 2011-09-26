************************************************************************
      SUBROUTINE IF037 (DA,VC,KCOL,KLD,DX,DB,DD,NEQ,NIT,ITE,EPS,OMEGA,
     *                  RHO)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),VC(*),KCOL(*),KLD(*),DX(*),DB(*),DD(*)
C
C
C
      CALL LCP1 (DB,DD,NEQ)
      CALL LAX17(DA,KCOL,KLD,NEQ,DX,DD,-1D0,1D0)
      CALL LL21 (DD,NEQ,FD)
      IF (ABS(FD).LT.1D-12) THEN
       ITE=1
       RHO=0D0
       RETURN
      ENDIF       
C
      DO 10 ITE=1,NIT
C
      CALL LCP1 (DB,DD,NEQ)
      CALL LAX17(DA,KCOL,KLD,NEQ,DX,DD,-1D0,1D0)
      CALL LL21 (DD,NEQ,RES)
      IF (ABS(RES).LT.EPS*FD) THEN
       RHO=(RES/FD)**(1D0/DBLE(ITE-1))
       RETURN       
      ENDIF
C
      CALL IF137(VC,KCOL,KLD,DD,NEQ)
C
      CALL LLC1 (DD,DX,NEQ,OMEGA,1D0)
C
10    CONTINUE
C
      RHO=(RES/FD)**(1D0/DBLE(NIT))
C
C
      END
