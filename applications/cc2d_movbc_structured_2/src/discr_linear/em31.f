      SUBROUTINE EM31(XI1,XI2,IPAR)
      
      IMPLICIT NONE

C include decessary COMMON blocks

      INCLUDE 'cbasictria.inc'

      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      INCLUDE 'ccub.inc'

C parameters
      
      DOUBLE PRECISION DXM(4),DYM(4),A(4,4),F(4)
      DOUBLE PRECISION CKH(4),CK(4,4)
      DOUBLE PRECISION XI1,XI2
      INTEGER IPAR

C local variables

      INTEGER IVE,IA,IK1,IK2,IK

C     Constants.
C     The EM30 element does only work for elements with NVE=4!

      INTEGER NVE
      PARAMETER (NVE=4)

C implicit function definitions

      DOUBLE PRECISION F1, F2, F3, F4
      DOUBLE PRECISION X,Y,CA1,CA2,CA3,CB1,CB2,CB3,CC3

      F1(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)=1D0
      F2(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)=CA1*X  +CB1*Y
      F3(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)=CA2*X  +CB2*Y
      F4(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)=CA3*X*X+CB3*X*Y+CC3*Y*Y
C
C
      SUB='EM31'
      IF (ICHECK.GE.998) CALL OTRC('EM31  ','01/08/94')
C
C
C *** Dummy call
      IF (IPAR.EQ.-1) THEN
       IER=0
       IPAR=31
       GOTO 99999
      ENDIF
C
C
      IF (IPAR.EQ.-2) THEN
       DO 20 IVE=1,NVE
       DXM(IVE)=0.5D0*(DX(IVE)+DX(MOD(IVE,4)+1))
       DYM(IVE)=0.5D0*(DY(IVE)+DY(MOD(IVE,4)+1))
20     CONTINUE
C
       CA1=(DXM(2)-DXM(4))/SQRT((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
       CB1=(DYM(2)-DYM(4))/SQRT((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
       CA2=(DXM(3)-DXM(1))/SQRT((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)
       CB2=(DYM(3)-DYM(1))/SQRT((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)
       CA3=CA1**2-CA2**2
       CB3=2D0*(CA1*CB1-CA2*CB2)
       CC3=CB1**2-CB2**2
C
       DO 22 IA=1,4
       A(IA,1)=F1(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)       
       A(IA,2)=F2(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)       
       A(IA,3)=F3(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)       
       A(IA,4)=F4(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
22     CONTINUE
C
       CALL INVERT(A,F,CKH,0)       
C
       DO 24 IK1=1,4
       DO 24 IK2=1,4
24     CK(IK1,IK2)=A(IK2,IK1)
C
       DO 26 IK=1,4
       COB(IK,1)=CK(IK,4)*CA3
       COB(IK,2)=CK(IK,4)*CC3
       COB(IK,3)=CK(IK,4)*CB3
       COB(IK,4)=CK(IK,2)*CA1+CK(IK,3)*CA2
       COB(IK,5)=CK(IK,2)*CB1+CK(IK,3)*CB2
       COB(IK,6)=CK(IK,1)
26     CONTINUE
      ENDIF
C
C

      IER=-1
C *** No second order derivatives available
C *** Used for second order problems only
      IF (BDER(4).OR.BDER(5).OR.BDER(6)) GOTO 99999
C
      IER=0
C
      IF (.NOT.BDER(1)) GOTO 101
C
C *** Function values
      DBAS(1,1)= COB(1,1)*XI1**2+COB(1,2)*XI2**2+COB(1,3)*XI1*XI2
     *          +COB(1,4)*XI1   +COB(1,5)*XI2   +COB(1,6)
      DBAS(2,1)= COB(2,1)*XI1**2+COB(2,2)*XI2**2+COB(2,3)*XI1*XI2
     *          +COB(2,4)*XI1   +COB(2,5)*XI2   +COB(2,6)
      DBAS(3,1)= COB(3,1)*XI1**2+COB(3,2)*XI2**2+COB(3,3)*XI1*XI2
     *          +COB(3,4)*XI1   +COB(3,5)*XI2   +COB(3,6)
      DBAS(4,1)= COB(4,1)*XI1**2+COB(4,2)*XI2**2+COB(4,3)*XI1*XI2
     *          +COB(4,4)*XI1   +COB(4,5)*XI2   +COB(4,6)
101   IF (.NOT.(BDER(2).OR.BDER(3))) GOTO 99999
C
C *** First order derivatives
      IF (.NOT.BDER(2)) GOTO 102
      DBAS(1,2)= 2D0*COB(1,1)*XI1+COB(1,3)*XI2+COB(1,4)
      DBAS(2,2)= 2D0*COB(2,1)*XI1+COB(2,3)*XI2+COB(2,4)
      DBAS(3,2)= 2D0*COB(3,1)*XI1+COB(3,3)*XI2+COB(3,4)
      DBAS(4,2)= 2D0*COB(4,1)*XI1+COB(4,3)*XI2+COB(4,4)
C
102   IF (.NOT.BDER(3)) GOTO 99999
      DBAS(1,3)= 2D0*COB(1,2)*XI2+COB(1,3)*XI1+COB(1,5)
      DBAS(2,3)= 2D0*COB(2,2)*XI2+COB(2,3)*XI1+COB(2,5)
      DBAS(3,3)= 2D0*COB(3,2)*XI2+COB(3,3)*XI1+COB(3,5)
      DBAS(4,3)= 2D0*COB(4,2)*XI2+COB(4,3)*XI1+COB(4,5)
C
C
99999 END
