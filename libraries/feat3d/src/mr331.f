      SUBROUTINE MR331(DF2,DF1,KAREA2,KAREA1,KADJ2,KADJ1,NAT2,NAT1,
     *                 NEL2,NEL1)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNAE=6)
      PARAMETER (A1= 11D0/24D0,A2= 7D0/48D0,A3=-5D0/48D0,A4=-1D0/24D0)
      PARAMETER (A5= 11D0/24D0,A6= 1D0/12D0,A7=-1D0/24D0)
      DIMENSION DF1(1),DF2(1),KAREA2(NNAE,*),KAREA1(NNAE,*),
     *          KADJ1(NNAE,*),KADJ2(NNAE,*)
C
C
C
      CALL LCL1(DF1,NAT1)
C
      DO 10 IEL1=1,NEL1
C
      IA1=KAREA1(1,IEL1)
      IA2=KAREA1(2,IEL1)
      IA3=KAREA1(3,IEL1)
      IA4=KAREA1(4,IEL1)
      IA5=KAREA1(5,IEL1)
      IA6=KAREA1(6,IEL1)
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
      DY1= DF2(KAREA2(1,IELH1))
      DY2= DF2(KAREA2(1,IELH2))
      DY3= DF2(KAREA2(1,IELH3))
      DY4= DF2(KAREA2(1,IELH4))
      DY5= DF2(KAREA2(2,IELH1))
      DY6= DF2(KAREA2(5,IELH2))
      DY7= DF2(KAREA2(5,IELH6))
      DY8= DF2(KAREA2(2,IELH5))
      DY9= DF2(KAREA2(2,IELH2))
      DY10=DF2(KAREA2(5,IELH3))
      DY11=DF2(KAREA2(5,IELH7))
      DY12=DF2(KAREA2(2,IELH6))
      DY13=DF2(KAREA2(2,IELH3))
      DY14=DF2(KAREA2(5,IELH4))
      DY15=DF2(KAREA2(5,IELH8))
      DY16=DF2(KAREA2(2,IELH7))
      DY17=DF2(KAREA2(2,IELH4))
      DY18=DF2(KAREA2(5,IELH1))
      DY19=DF2(KAREA2(5,IELH5))
      DY20=DF2(KAREA2(2,IELH8))
      DY21=DF2(KAREA2(1,IELH5))
      DY22=DF2(KAREA2(1,IELH6))
      DY23=DF2(KAREA2(1,IELH7))
      DY24=DF2(KAREA2(1,IELH8))
      DY25=DF2(KAREA2(3,IELH1))
      DY26=DF2(KAREA2(3,IELH2))
      DY27=DF2(KAREA2(3,IELH3))
      DY28=DF2(KAREA2(3,IELH4))
      DY29=DF2(KAREA2(6,IELH1))
      DY30=DF2(KAREA2(6,IELH2))
      DY31=DF2(KAREA2(6,IELH3))
      DY32=DF2(KAREA2(6,IELH4))
      DY33=DF2(KAREA2(3,IELH5))
      DY34=DF2(KAREA2(3,IELH6))
      DY35=DF2(KAREA2(3,IELH7))
      DY36=DF2(KAREA2(3,IELH8))
C
      IF (KADJ1(1,IEL1).NE.0) THEN
       DF1(IA1)= DF1(IA1)
     *          +A1*(DY1+DY2+DY3+DY4)
     *          +A2*(DY5+DY6+DY9+DY10+DY13+DY14+DY17+DY18)
     *          +A3*(DY7+DY8+DY11+DY12+DY15+DY16+DY19+DY20)
     *          +A4*(DY21+DY22+DY23+DY24)
     *          +A5*(DY25+DY26+DY27+DY28)
     *          +A6*(DY29+DY30+DY31+DY32)
     *          +A7*(DY33+DY34+DY35+DY36)     
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN
       DF1(IA2)= DF1(IA2)
     *          +A1*(DY5+DY6+DY7+DY8)
     *          +A2*(DY1+DY2+DY9+DY12+DY21+DY22+DY18+DY19)
     *          +A3*(DY3+DY4+DY10+DY11+DY23+DY24+DY17+DY20)
     *          +A4*(DY13+DY14+DY15+DY16)
     *          +A5*(DY25+DY29+DY30+DY33)
     *          +A6*(DY26+DY28+DY34+DY36)
     *          +A7*(DY27+DY31+DY32+DY35)     
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN
       DF1(IA3)= DF1(IA3)
     *          +A1*(DY9+DY10+DY11+DY12)
     *          +A2*(DY2+DY3+DY6+DY7+DY22+DY23+DY13+DY16)
     *          +A3*(DY1+DY4+DY5+DY8+DY21+DY24+DY14+DY15)
     *          +A4*(DY17+DY18+DY19+DY20)
     *          +A5*(DY26+DY30+DY31+DY34)
     *          +A6*(DY25+DY27+DY33+DY35)
     *          +A7*(DY28+DY29+DY32+DY36)     
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN
       DF1(IA4)= DF1(IA4)
     *          +A1*(DY13+DY14+DY15+DY16)
     *          +A2*(DY3+DY4+DY10+DY11+DY23+DY24+DY17+DY20)
     *          +A3*(DY1+DY2+DY9+DY12+DY21+DY22+DY18+DY19)
     *          +A4*(DY5+DY6+DY7+DY8)
     *          +A5*(DY27+DY31+DY32+DY35)  
     *          +A6*(DY26+DY28+DY34+DY36)
     *          +A7*(DY25+DY29+DY30+DY33)
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN
       DF1(IA5)= DF1(IA5)
     *          +A1*(DY17+DY18+DY19+DY20)
     *          +A2*(DY1+DY4+DY5+DY8+DY21+DY24+DY14+DY15)
     *          +A3*(DY2+DY3+DY6+DY7+DY22+DY23+DY13+DY16)
     *          +A4*(DY9+DY10+DY11+DY12)
     *          +A5*(DY28+DY29+DY32+DY36)     
     *          +A6*(DY25+DY27+DY33+DY35)
     *          +A7*(DY26+DY30+DY31+DY34)
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN
       DF1(IA6)= DF1(IA6)
     *          +A1*(DY21+DY22+DY23+DY24)
     *          +A2*(DY7+DY8+DY11+DY12+DY15+DY16+DY19+DY20)
     *          +A3*(DY5+DY6+DY9+DY10+DY13+DY14+DY17+DY18)
     *          +A4*(DY1+DY2+DY3+DY4)
     *          +A5*(DY33+DY34+DY35+DY36)     
     *          +A6*(DY29+DY30+DY31+DY32)
     *          +A7*(DY25+DY26+DY27+DY28)
      ENDIF
C
C
10    CONTINUE
C
C
      END
