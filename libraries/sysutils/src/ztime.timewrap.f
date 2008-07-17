      SUBROUTINE ZTIME(ERG)
C
C
      EXTERNAL TIME
      DOUBLE PRECISION ERG
      INTEGER R,C,T
      INTEGER TIME
      
      T=TWRAP(C,R)
C           
      ERG=DBLE(R)/DBLE(C)
      
      END


