************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* LTX3n                                                                *
*                                                                      *
* Purpose  Matrix vector product with transposed matrix                *
*          Matrix stored in technique  n  (see Reference Manual)       *
*          Single/double precision version                             *
*                                                                      *
* Subroutines/functions called   LSC1, LCL1, LVM3                      *
*                                                                      *
* Version from  11/15/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VA       R*4    Matrix stored in technique  n                        *
* NEQ      I*4    Number of equations (length of DX, DTX)              *
*                 For storage technique 9: NEQ1 (length of DX) and     *
*                 NEQ2 (length of DTX) are needed                      *
* DX       R*8    Vector                                               *
* A1,A2    R*8    DTX := A1*VAT*DX + A2*DTX                            *
* KLD      I*4    Pointer vectors corresponding to the                 *
* KCOL     I*4    storage technique                                    *
* KOP      I*4                                                         *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DTX      R*8    Resulting vector                                     *
*                                                                      *
************************************************************************
C
      SUBROUTINE LTX33(VA,KDIA,KDIAS,NDIA,NEQ,DX,DTX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KDIA(*),KDIAS(*),DX(*),DTX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LTX33 ','11/19/90')
C
      IF (A1.NE.0D0) THEN
       IF (A2.EQ.0D0) THEN
	  CALL LVM3(DX,VA,DTX,NEQ,1D0,0D0)
       ELSE
        A=A2/A1
        IF (A.NE.1D0) CALL LSC1(DTX,NEQ,A)
	  CALL LVM3(DX,VA,DTX,NEQ,1D0,1D0)
       ENDIF
       DO 1 IDIA=2,NDIA
       J1=-KDIA(IDIA)
       IF (J1.LT.0) THEN
        I1=1-J1
	  NEQ1=NEQ+J1
       ELSE
        I1=1
	  NEQ1=NEQ-J1
       ENDIF
       J0=KDIAS(IDIA)-I1
       CALL LVM3(DX(I1+J1),VA(I1+J0),DTX(I1),NEQ1,1D0,1D0)
1      CONTINUE
       IF (A1.NE.1D0) CALL LSC1(DTX,NEQ,A1)
      ELSE
       CALL LSC1(DTX,NEQ,A2)
      ENDIF
C
      END
C
C
C
      SUBROUTINE LTX37(VA,KCOL,KLD,NEQ,DX,DTX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),DX(*),DTX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LTX37 ','01/02/89')
C
      IF (A1.NE.0D0) THEN
       IF (A2.EQ.0D0) THEN
        DO 1 IROW=1,NEQ
1       DTX(IROW)=DX(IROW)*DBLE(VA(KLD(IROW)))
       ELSE
        A=A2/A1
        IF (A.NE.1D0) CALL LSC1(DTX,NEQ,A)
        DO 2 IROW=1,NEQ
2       DTX(IROW)=DX(IROW)*DBLE(VA(KLD(IROW)))+DTX(IROW)
       ENDIF
       DO 4 IROW=1,NEQ
       DO 3 ICOL=KLD(IROW)+1,KLD(IROW+1)-1
       JCOL=KCOL(ICOL)
3      DTX(JCOL)=DTX(JCOL)+DBLE(VA(ICOL))*DX(IROW)
4      CONTINUE
       IF (A1.NE.1D0) CALL LSC1(DTX,NEQ,A1)
      ELSE
       CALL LSC1(DTX,NEQ,A2)
      ENDIF
C
      END
C
C
C
      SUBROUTINE LTX39(VA,KCOL,KLD,NEQ1,NEQ2,DX,DTX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),DX(*),DTX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LTX39 ','11/15/89')
C
      IF (A1.NE.0D0) THEN
       IF (A2.EQ.0D0) THEN
        CALL LCL1(DTX,NEQ2)
       ELSE
        A=A2/A1
        IF (A.NE.1D0) CALL LSC1(DTX,NEQ2,A)
       ENDIF
       DO 4 IROW=1,NEQ1
       DO 3 ICOL=KLD(IROW),KLD(IROW+1)-1
       JCOL=KCOL(ICOL)
3      DTX(JCOL)=DTX(JCOL)+DBLE(VA(ICOL))*DX(IROW)
4      CONTINUE
       IF (A1.NE.1D0) CALL LSC1(DTX,NEQ2,A1)
      ELSE
       CALL LSC1(DTX,NEQ2,A2)
      ENDIF
C
      END
C
C
C
      SUBROUTINE LTX3A(VA,KCOL,KLD,KOP,NEQ,DX,DTX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),KOP(*),DX(*),DTX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LTX3A ','01/02/89')
C
      IF (A1.NE.0D0) THEN
       IF (A2.EQ.0D0) THEN
        CALL LCL1(DTX,NEQ)
       ELSE
        A=A2/A1
        IF (A.NE.1D0) CALL LSC1(DTX,NEQ,A)
       ENDIF
       DO 3 IROW=1,NEQ
       JOP=KOP(IROW)
       DO 4 ICOL=KLD(JOP),KLD(JOP+1)-1
       JCOL=IROW+KCOL(ICOL)
4      DTX(JCOL)=DTX(JCOL)+DBLE(VA(ICOL))*DX(IROW)
3      CONTINUE
       IF (A1.NE.1D0) CALL LSC1(DTX,NEQ,A1)
      ELSE
       CALL LSC1(DTX,NEQ,A2)
      ENDIF
C
      END
