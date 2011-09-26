************************************************************************
      SUBROUTINE IF227 (VA,VC,KCOL,KLD,DX,DB,DD,NEQ,NIT,OMEGA)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),VC(*),KCOL(*),KLD(*),DX(*),DB(*),DD(*)
C
      DO 10 ITE=1,NIT
      CALL LCP1 (DB,DD,NEQ)
      CALL LAX37(VA,KCOL,KLD,NEQ,DX,DD,-1D0,1D0)
      CALL IF137(VC,KCOL,KLD,DD,NEQ)
      CALL LLC1 (DD,DX,NEQ,OMEGA,1D0)
10    CONTINUE
C
      END
