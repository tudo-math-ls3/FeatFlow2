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
* LAX1n                                                                *
*                                                                      *
* Purpose  Matrix vector product                                       *
*          Matrix stored in technique  n  (see Reference Manual)       *
*          Double/double precision version                             *
*                                                                      *
* Subroutines/functions called   LSC1, LCL1, LVM1                      *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DA       R*8    Matrix stored in technique  n                        *
* NEQ      I*4    Number of equations (length of DX, DAX)              *
* DX       R*8    Vector                                               *
* A1,A2    R*8    DAX := A1*DA*DX + A2*DAX                             *
* KLD      I*4    Pointer vectors corresponding to the                 *
* KCOL     I*4    storage technique                                    *
* KOP      I*4                                                         *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DAX      R*8    Resulting vector                                     *
*                                                                      *
************************************************************************
C
      SUBROUTINE LAX13(DA,KDIA,KDIAS,NDIA,NEQ,DX,DAX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KDIA(*),KDIAS(*),DX(*),DAX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LAX13 ','11/19/90')
C
      IF (A1.NE.0D0) THEN
       IF (A2.EQ.0D0) THEN
        CALL LVM1(DX,DA,DAX,NEQ,1D0,0D0)
       ELSE
        A=A2/A1
        IF (A.NE.1D0) CALL LSC1(DAX,NEQ,A)
        CALL LVM1(DX,DA,DAX,NEQ,1D0,1D0)
       ENDIF
       DO 1 IDIA=2,NDIA
       J1=KDIA(IDIA)
       IF (J1.GT.0) THEN
        I1=1
	  NEQ1=NEQ-J1
       ELSE
        I1=1-J1
	  NEQ1=NEQ+J1
       ENDIF
       J0=KDIAS(IDIA)-I1
       CALL LVM1(DX(I1+J1),DA(I1+J0),DAX(I1),NEQ1,1D0,1D0)
1      CONTINUE
       IF (A1.NE.1D0) CALL LSC1(DAX,NEQ,A1)
      ELSE
       CALL LSC1(DAX,NEQ,A2)
      ENDIF
C
      END
C
C
C
      SUBROUTINE LAX14(DA,KDIA,KDIAS,NDIA,NEQ,DX,DAX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KDIA(*),KDIAS(*),DX(*),DAX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LAX14 ','11/19/90')
C
      IF (A1.NE.0D0) THEN
       IF (A2.EQ.0D0) THEN
        CALL LVM1(DX,DA,DAX,NEQ,1D0,0D0)
       ELSE
        A=A2/A1
        IF (A.NE.1D0) CALL LSC1(DAX,NEQ,A)
        CALL LVM1(DX,DA,DAX,NEQ,1D0,1D0)
       ENDIF
       DO 1 IDIA=2,NDIA
       J0=KDIAS(IDIA)-1
       J1=KDIA(IDIA)
	 NEQ1=NEQ-J1
       CALL LVM1(DX(1+J1),DA(1+J0),DAX,NEQ1,1D0,1D0)
       CALL LVM1(DX(1+J1),DA(1+J0),DAX(1+J1),NEQ1,1D0,1D0)
1      CONTINUE
       IF (A1.NE.1D0) CALL LSC1(DAX,NEQ,A1)
      ELSE
       CALL LSC1(DAX,NEQ,A2)
      ENDIF
C
      END
C
C
C
      SUBROUTINE LAX17(DA,KCOL,KLD,NEQ,DX,DAX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),DX(*),DAX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LAX17 ','01/02/89')
C
      IF (A1.NE.0D0) THEN
       IF (A2.EQ.0D0) THEN
        DO 1 IROW=1,NEQ
1       DAX(IROW)=DX(IROW)*DA(KLD(IROW))
       ELSE
        A=A2/A1
        IF (A.NE.1D0) CALL LSC1(DAX,NEQ,A)
        DO 2 IROW=1,NEQ
2       DAX(IROW)=DX(IROW)*DA(KLD(IROW))+DAX(IROW)
       ENDIF
       DO 6 IROW=1,NEQ
       DO 5 ICOL=KLD(IROW)+1,KLD(IROW+1)-1
5      DAX(IROW)=DAX(IROW)+DA(ICOL)*DX(KCOL(ICOL))
6      CONTINUE
       IF (A1.NE.1D0) CALL LSC1(DAX,NEQ,A1)
      ELSE
       CALL LSC1(DAX,NEQ,A2)
      ENDIF
C
      END
C
C
C
      SUBROUTINE LAX18(DA,KCOL,KLD,NEQ,DX,DAX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),DX(*),DAX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LAX18 ','01/02/89')
C
      IF (A1.NE.0D0) THEN
       IF (A2.EQ.0D0) THEN
        DO 1 IROW=1,NEQ
1       DAX(IROW)=DA(KLD(IROW))*DX(IROW)
       ELSE
        A=A2/A1
        IF (A.NE.1D0) CALL LSC1(DAX,NEQ,A)
        DO 2 IROW=1,NEQ
2       DAX(IROW)=DX(IROW)*DA(KLD(IROW))+DAX(IROW)
       ENDIF
       DO 5 IROW=1,NEQ-1
       FAK=0D0
       DO 6 ICOL=KLD(IROW)+1,KLD(IROW+1)-1
       JCOL=KCOL(ICOL)
       FAK=FAK+DA(ICOL)*DX(JCOL)
       DAX(JCOL)=DAX(JCOL)+DA(ICOL)*DX(IROW)
6      CONTINUE
5      DAX(IROW)=DAX(IROW)+FAK
       IF (A1.NE.1D0) CALL LSC1(DAX,NEQ,A1)
      ELSE
       CALL LSC1(DAX,NEQ,A2)
      ENDIF
C
      END
C
C
C
      SUBROUTINE LAX19(DA,KCOL,KLD,NEQ,DX,DAX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),DX(*),DAX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LAX19 ','01/02/89')
C
      IF (A1.NE.0D0) THEN
       IF (A2.EQ.0D0) THEN
        DO 1 IROW=1,NEQ
        ICOL=KCOL(KLD(IROW))
1       DAX(IROW)=DX(ICOL)*DA(KLD(IROW))
       ELSE
        A=A2/A1
        IF (A.NE.1D0) CALL LSC1(DAX,NEQ,A)
        DO 2 IROW=1,NEQ
        ICOL=KCOL(KLD(IROW))
2       DAX(IROW)=DX(ICOL)*DA(KLD(IROW))+DAX(IROW)
       ENDIF
       DO 6 IROW=1,NEQ
       DO 5 ICOL=KLD(IROW)+1,KLD(IROW+1)-1
5      DAX(IROW)=DAX(IROW)+DA(ICOL)*DX(KCOL(ICOL))
6      CONTINUE
       IF (A1.NE.1D0) CALL LSC1(DAX,NEQ,A1)
      ELSE
       CALL LSC1(DAX,NEQ,A2)
      ENDIF
C
      END
C
C
C
      SUBROUTINE LAX1A(DA,KCOL,KLD,KOP,NEQ,DX,DAX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),KOP(*),DX(*),DAX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LAX1A ','01/02/89')
C
      IF (A1.NE.0D0) THEN
       IF (A2.EQ.0D0) THEN
        CALL LCL1(DAX,NEQ)
       ELSE
        A=A2/A1
        IF (A.NE.1D0) CALL LSC1(DAX,NEQ,A)
       ENDIF
       DO 4 IROW=1,NEQ
       JOP=KOP(IROW)
       DO 5 ICOL=KLD(JOP),KLD(JOP+1)-1
5      DAX(IROW)=DAX(IROW)+DA(ICOL)*DX(KCOL(ICOL)+IROW)
4      CONTINUE
       IF (A1.NE.1D0) CALL LSC1(DAX,NEQ,A1)
      ELSE
       CALL LSC1(DAX,NEQ,A2)
      ENDIF
C
      END
