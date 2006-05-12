      SUBROUTINE MR311(DF2,DF1,KVERT2,KVERT1,KADJ2,KADJ1,
     *                 NVT2,NVT1,NEL2,NEL1)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION DF1(1),DF2(1),KVERT2(NNVE,*),KVERT1(NNVE,*),
     *          KADJ2(NNAE,*),KADJ1(NNAE,*)
C
C
C
      CALL LCP1(DF2,DF1,NVT1)
      DO 10 IEL2=1,NEL2
C
      I1=KVERT2(1,IEL2)
      I3=KVERT2(3,IEL2)
      I6=KVERT2(6,IEL2)
      I7=KVERT2(7,IEL2)
      I8=KVERT2(8,IEL2)
C
      IF (KADJ2(1,IEL2).EQ.0) THEN
       A3=0.25D0
      ELSE
       A3=0.125D0
      ENDIF
C
      IF (KADJ2(2,IEL2).EQ.0) THEN
       A6=0.25D0
      ELSE
       A6=0.125D0
      ENDIF
C
      IF (KADJ2(5,IEL2).EQ.0) THEN
       A8=0.25D0
      ELSE
       A8=0.125D0
      ENDIF
C
      DF1(I1)=DF1(I1)+A3*DF2(I3)+A6*DF2(I6)+0.125D0*DF2(I7)+A8*DF2(I8)
C
10    CONTINUE
C
      IVT=NVT1
C
      DO 20 IEL1=1,NEL1
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(3,IELH5)
      IELH7=KADJ2(3,IELH6)
      IELH8=KADJ2(3,IELH7)
C
      I1 =KVERT2(2,IELH1)
      I2 =KVERT2(2,IELH2)
      I3 =KVERT2(2,IELH3)
      I4 =KVERT2(2,IELH4)
      I5 =KVERT2(5,IELH1)
      I6 =KVERT2(5,IELH2)
      I7 =KVERT2(5,IELH3)
      I8 =KVERT2(5,IELH4)
      I9 =KVERT2(2,IELH5)
      I10=KVERT2(2,IELH6)
      I11=KVERT2(2,IELH7)
      I12=KVERT2(2,IELH8)
C
      J1=KVERT1(1,IEL1)
      J2=KVERT1(2,IEL1)
      J3=KVERT1(3,IEL1)
      J4=KVERT1(4,IEL1)
      J5=KVERT1(5,IEL1)
      J6=KVERT1(6,IEL1)
      J7=KVERT1(7,IEL1)
      J8=KVERT1(8,IEL1)
C
      IF (I1.GT.IVT) THEN
       DF1(J1)=DF1(J1)+0.5D0*DF2(I1)
       DF1(J2)=DF1(J2)+0.5D0*DF2(I1)
       IVT=IVT+1
      ENDIF
C
      IF (I2.GT.IVT) THEN
       DF1(J2)=DF1(J2)+0.5D0*DF2(I2)
       DF1(J3)=DF1(J3)+0.5D0*DF2(I2)
       IVT=IVT+1
      ENDIF
C
      IF (I3.GT.IVT) THEN
       DF1(J3)=DF1(J3)+0.5D0*DF2(I3)
       DF1(J4)=DF1(J4)+0.5D0*DF2(I3)
       IVT=IVT+1
      ENDIF
C
      IF (I4.GT.IVT) THEN
       DF1(J1)=DF1(J1)+0.5D0*DF2(I4)
       DF1(J4)=DF1(J4)+0.5D0*DF2(I4)
       IVT=IVT+1
      ENDIF
C
      IF (I5.GT.IVT) THEN
       DF1(J1)=DF1(J1)+0.5D0*DF2(I5)
       DF1(J5)=DF1(J5)+0.5D0*DF2(I5)
       IVT=IVT+1
      ENDIF
C
      IF (I6.GT.IVT) THEN
       DF1(J2)=DF1(J2)+0.5D0*DF2(I6)
       DF1(J6)=DF1(J6)+0.5D0*DF2(I6)
       IVT=IVT+1
      ENDIF
C
      IF (I7.GT.IVT) THEN
       DF1(J3)=DF1(J3)+0.5D0*DF2(I7)
       DF1(J7)=DF1(J7)+0.5D0*DF2(I7)
       IVT=IVT+1
      ENDIF
C
      IF (I8.GT.IVT) THEN
       DF1(J4)=DF1(J4)+0.5D0*DF2(I8)
       DF1(J8)=DF1(J8)+0.5D0*DF2(I8)
       IVT=IVT+1
      ENDIF
C
      IF (I9.GT.IVT) THEN
       DF1(J5)=DF1(J5)+0.5D0*DF2(I9)
       DF1(J6)=DF1(J6)+0.5D0*DF2(I9)
       IVT=IVT+1
      ENDIF
C
      IF (I10.GT.IVT) THEN
       DF1(J6)=DF1(J6)+0.5D0*DF2(I10)
       DF1(J7)=DF1(J7)+0.5D0*DF2(I10)
       IVT=IVT+1
      ENDIF
C
      IF (I11.GT.IVT) THEN
       DF1(J7)=DF1(J7)+0.5D0*DF2(I11)
       DF1(J8)=DF1(J8)+0.5D0*DF2(I11)
       IVT=IVT+1
      ENDIF
C
      IF (I12.GT.IVT) THEN
       DF1(J5)=DF1(J5)+0.5D0*DF2(I12)
       DF1(J8)=DF1(J8)+0.5D0*DF2(I12)
       IVT=IVT+1
      ENDIF
C
20    CONTINUE
C
C
C
      END
