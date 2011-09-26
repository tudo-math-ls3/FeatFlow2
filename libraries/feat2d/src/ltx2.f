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
* LTX2n                                                                *
*                                                                      *
* Purpose  Matrix vector product with transposed matrix                *
*          Matrix stored in technique  n  (see Reference Manual)       *
*          Single/Single precision version                             *
*                                                                      *
* Subroutines/functions called   LSC2, LCL2, LVM2                      *
*                                                                      *
* Version from  01/15/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VA       R*4    Matrix stored in technique  n                        *
* NEQ      I*4    Number of equations (length of VX, VTX)              *
*                 For storage technique 9: NEQ1 (length of VX) and     *
*                 NEQ2 (length of VTX) are needed                      *
* VX       R*4    Vector                                               *
* A1,A2    R*8    VTX := A1*VAT*VX + A2*VTX                            *
* KLD      I*4    Pointer vectors corresponding to the                 *
* KCOL     I*4    storage technique                                    *
* KOP      I*4                                                         *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* VTX      R*4    Resulting vector                                     *
*                                                                      *
************************************************************************
C
      SUBROUTINE LTX23(VA,KDIA,KDIAS,NDIA,NEQ,VX,VTX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KDIA(*),KDIAS(*),VX(*),VTX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LTX23 ','11/19/90')
C
      IF (A1.NE.0D0) THEN
       IF (A2.EQ.0D0) THEN
	  CALL LVM2(VX,VA,VTX,NEQ,1D0,0D0)
       ELSE
        A=A2/A1
        IF (A.NE.1D0) CALL LSC2(VTX,NEQ,A)
	  CALL LVM2(VX,VA,VTX,NEQ,1D0,1D0)
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
       CALL LVM2(VX(I1+J1),VA(I1+J0),VTX(I1),NEQ1,1D0,1D0)
1      CONTINUE
       IF (A1.NE.1D0) CALL LSC2(VTX,NEQ,A1)
      ELSE
       CALL LSC2(VTX,NEQ,A2)
      ENDIF
C
      END
C
C
C
      SUBROUTINE LTX27(VA,KCOL,KLD,NEQ,VX,VTX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),VX(*),VTX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LTX27 ','01/02/89')
C
      IF (A1.NE.0D0) THEN
       IF (A2.EQ.0D0) THEN
        DO 1 IROW=1,NEQ
1       VTX(IROW)=DBLE(VX(IROW))*DBLE(VA(KLD(IROW)))
       ELSE
        A=A2/A1
        IF (A.NE.1D0) CALL LSC2(VTX,NEQ,A)
        DO 2 IROW=1,NEQ
2       VTX(IROW)=DBLE(VX(IROW))*DBLE(VA(KLD(IROW)))+DBLE(VTX(IROW))
       ENDIF
       DO 4 IROW=1,NEQ
       DO 3 ICOL=KLD(IROW)+1,KLD(IROW+1)-1
       JCOL=KCOL(ICOL)
3      VTX(JCOL)=DBLE(VTX(JCOL))+DBLE(VA(ICOL))*DBLE(VX(IROW))
4      CONTINUE
       IF (A1.NE.1D0) CALL LSC2(VTX,NEQ,A1)
      ELSE
       CALL LSC2(VTX,NEQ,A2)
      ENDIF
C
      END
C
C
C
      SUBROUTINE LTX29(VA,KCOL,KLD,NEQ1,NEQ2,VX,VTX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),VX(*),VTX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LTX29 ','11/15/89')
C
      IF (A1.NE.0D0) THEN
       IF (A2.EQ.0D0) THEN
        CALL LCL2(VTX,NEQ2)
       ELSE
        A=A2/A1
        IF (A.NE.1D0) CALL LSC2(VTX,NEQ2,A)
       ENDIF
       DO 4 IROW=1,NEQ1
       DO 3 ICOL=KLD(IROW),KLD(IROW+1)-1
       JCOL=KCOL(ICOL)
3      VTX(JCOL)=DBLE(VTX(JCOL))+DBLE(VA(ICOL))*DBLE(VX(IROW))
4      CONTINUE
       IF (A1.NE.1D0) CALL LSC2(VTX,NEQ2,A1)
      ELSE
       CALL LSC2(VTX,NEQ2,A2)
      ENDIF
C
      END
C
C
C
      SUBROUTINE LTX2A(VA,KCOL,KLD,KOP,NEQ,VX,VTX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),KOP(*),VX(*),VTX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LTX2A ','01/02/89')
C
      IF (A1.NE.0D0) THEN
       IF (A2.EQ.0D0) THEN
        CALL LCL2(VTX,NEQ)
       ELSE
        A=A2/A1
        IF (A.NE.1D0) CALL LSC2(VTX,NEQ,A)
       ENDIF
       DO 3 IROW=1,NEQ
       JOP=KOP(IROW)
       DO 4 ICOL=KLD(JOP),KLD(JOP+1)-1
       JCOL=IROW+KCOL(ICOL)
4      VTX(JCOL)=DBLE(VTX(JCOL))+
     *           DBLE(VA(ICOL))*DBLE(VX(IROW))
3      CONTINUE
       IF (A1.NE.1D0) CALL LSC2(VTX,NEQ,A1)
      ELSE
       CALL LSC2(VTX,NEQ,A2)
      ENDIF
C
      END
