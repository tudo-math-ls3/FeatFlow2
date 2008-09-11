      SUBROUTINE MP010 (DP1,DP2,KADJ1,KADJ2,NEL1,NEL2)

C *** Prolongation of pressure

      IMPLICIT NONE

      INCLUDE 'cbasictria.inc'
      
C parameters
      
      DOUBLE PRECISION DP1(*),DP2(*)
      INTEGER KADJ1(NNVE,*),KADJ2(NNVE,*)

      INTEGER NEL1, NEL2
      
C local variables

      INTEGER IEL1, IELH1, IELH2, IELH3, IELH4
      DOUBLE PRECISION DPH

C-----------------------------------------------------------------------
C
      DO 10 IEL1=1,NEL1
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)

      DPH=DP1(IEL1)
      DP2(IELH1)=DPH
      DP2(IELH2)=DPH
      DP2(IELH3)=DPH
      DP2(IELH4)=DPH
C
10    CONTINUE
C
C
      END

************************************************************************
      SUBROUTINE MR010 (DP1,DP2,KADJ1,KADJ2,NEL1,NEL2)

C *** Restriction of pressure

      IMPLICIT NONE

      INCLUDE 'cbasictria.inc'
      
C parameters
      
      DOUBLE PRECISION DP1(*),DP2(*)
      INTEGER KADJ1(NNVE,*),KADJ2(NNVE,*)

      INTEGER NEL1, NEL2
      
C local variables

      INTEGER IEL1, IELH1, IELH2, IELH3, IELH4

C-----------------------------------------------------------------------
C
      DO 10 IEL1=1,NEL1
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)

      DP1(IEL1)= DP2(IELH1)+DP2(IELH2)+DP2(IELH3)+DP2(IELH4)

10    CONTINUE

      END
