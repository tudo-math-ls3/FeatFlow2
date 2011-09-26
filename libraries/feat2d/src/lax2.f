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
* LAX2n                                                                *
*                                                                      *
* Purpose  Matrix vector product                                       *
*          Matrix stored in technique  n  (see Reference Manual)       *
*          Single/Single precision version                             *
*                                                                      *
* Subroutines/functions called   LSC2, LCL2, LVM2                      *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VA       R*4    Matrix stored in technique  n                        *
* NEQ      I*4    Number of equations (length of VX, VAX)              *
* VX       R*4    Vector                                               *
* A1,A2    R*8    VAX := A1*VA*VX + A2*VAX                             *
* KLD      I*4    Pointer vectors corresponding to the                 *
* KCOL     I*4    storage technique                                    *
* KOP      I*4                                                         *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* VAX      R*4    Resulting vector                                     *
*                                                                      *
************************************************************************
C
      SUBROUTINE LAX23(VA,KDIA,KDIAS,NDIA,NEQ,VX,VAX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KDIA(*),KDIAS(*),VX(*),VAX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LAX23 ','11/19/90')
C
      IF (A1.NE.0D0) THEN
       IF (A2.EQ.0D0) THEN
        CALL LVM2(VX,VA,VAX,NEQ,1D0,0D0)
       ELSE
        A=A2/A1
        IF (A.NE.1D0) CALL LSC2(VAX,NEQ,A)
        CALL LVM2(VX,VA,VAX,NEQ,1D0,1D0)
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
       CALL LVM2(VX(I1+J1),VA(I1+J0),VAX(I1),NEQ1,1D0,1D0)
1      CONTINUE
       IF (A1.NE.1D0) CALL LSC2(VAX,NEQ,A1)
      ELSE
       CALL LSC2(VAX,NEQ,A2)
      ENDIF
C
      END
C
C
C
      SUBROUTINE LAX24(VA,KDIA,KDIAS,NDIA,NEQ,VX,VAX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KDIA(*),KDIAS(*),VX(*),VAX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LAX24 ','11/19/90')
C
      IF (A1.NE.0D0) THEN
       IF (A2.EQ.0D0) THEN
        CALL LVM2(VX,VA,VAX,NEQ,1D0,0D0)
       ELSE
        A=A2/A1
        IF (A.NE.1D0) CALL LSC2(VAX,NEQ,A)
        CALL LVM2(VX,VA,VAX,NEQ,1D0,1D0)
       ENDIF
       DO 1 IDIA=2,NDIA
       J0=KDIAS(IDIA)-1
       J1= KDIA(IDIA)
	 NEQ1=NEQ-J1
       CALL LVM2(VX(1+J1),VA(1+J0),VAX,NEQ1,1D0,1D0)
       CALL LVM2(VX(1+J1),VA(1+J0),VAX(1+J1),NEQ1,1D0,1D0)
1      CONTINUE
       IF (A1.NE.1D0) CALL LSC2(VAX,NEQ,A1)
      ELSE
       CALL LSC2(VAX,NEQ,A2)
      ENDIF
C
      END
C
C
C
      SUBROUTINE LAX27(VA,KCOL,KLD,NEQ,VX,VAX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),VX(*),VAX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LAX27 ','01/02/89')
C
      IF (A1.NE.0D0) THEN
       IF (A2.EQ.0D0) THEN
        DO 1 IROW=1,NEQ
1       VAX(IROW)=DBLE(VX(IROW))*DBLE(VA(KLD(IROW)))
       ELSE
        A=A2/A1
        IF (A.NE.1D0) CALL LSC2(VAX,NEQ,A)
        DO 2 IROW=1,NEQ
2       VAX(IROW)=DBLE(VX(IROW))*DBLE(VA(KLD(IROW)))+DBLE(VAX(IROW))
       ENDIF
       DO 6 IROW=1,NEQ
       DO 5 ICOL=KLD(IROW)+1,KLD(IROW+1)-1
5      VAX(IROW)=DBLE(VAX(IROW))+DBLE(VA(ICOL))*DBLE(VX(KCOL(ICOL)))
6      CONTINUE
       IF (A1.NE.1D0) CALL LSC2(VAX,NEQ,A1)
      ELSE
       CALL LSC2(VAX,NEQ,A2)
      ENDIF
C
      END
C
C
C
      SUBROUTINE LAX28(VA,KCOL,KLD,NEQ,VX,VAX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),VX(*),VAX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LAX28 ','01/02/89')
C
      IF (A1.NE.0D0) THEN
       IF (A2.EQ.0D0) THEN
        DO 1 IROW=1,NEQ
1       VAX(IROW)=DBLE(VA(KLD(IROW)))*DBLE(VX(IROW))
       ELSE
        A=A2/A1
        IF (A.NE.1D0) CALL LSC2(VAX,NEQ,A)
        DO 2 IROW=1,NEQ
2       VAX(IROW)=DBLE(VX(IROW))*DBLE(VA(KLD(IROW)))+DBLE(VAX(IROW))
       ENDIF
       DO 5 IROW=1,NEQ-1
       FAK=0D0
       DO 6 ICOL=KLD(IROW)+1,KLD(IROW+1)-1
       JCOL=KCOL(ICOL)
       FAK=FAK+DBLE(VA(ICOL))*DBLE(VX(JCOL))
       VAX(JCOL)=DBLE(VAX(JCOL))+DBLE(VA(ICOL))*DBLE(VX(IROW))
6      CONTINUE
5      VAX(IROW)=DBLE(VAX(IROW))+FAK
       IF (A1.NE.1D0) CALL LSC2(VAX,NEQ,A1)
      ELSE
       CALL LSC2(VAX,NEQ,A2)
      ENDIF
C
      END
C
C
C
      SUBROUTINE LAX29(VA,KCOL,KLD,NEQ,VX,VAX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),VX(*),VAX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LAX29 ','01/02/89')
C
      IF (A1.NE.0D0) THEN
       IF (A2.EQ.0D0) THEN
        DO 1 IROW=1,NEQ
        ICOL=KCOL(KLD(IROW))
1       VAX(IROW)=DBLE(VX(ICOL))*DBLE(VA(KLD(IROW)))
       ELSE
        A=A2/A1
        IF (A.NE.1D0) CALL LSC2(VAX,NEQ,A)
        DO 2 IROW=1,NEQ
        ICOL=KCOL(KLD(IROW))
2       VAX(IROW)=DBLE(VX(ICOL))*DBLE(VA(KLD(IROW)))+DBLE(VAX(IROW))
       ENDIF
       DO 6 IROW=1,NEQ
       DO 5 ICOL=KLD(IROW)+1,KLD(IROW+1)-1
5      VAX(IROW)=DBLE(VAX(IROW))+DBLE(VA(ICOL))*DBLE(VX(KCOL(ICOL)))
6      CONTINUE
       IF (A1.NE.1D0) CALL LSC2(VAX,NEQ,A1)
      ELSE
       CALL LSC2(VAX,NEQ,A2)
      ENDIF
C
      END
C
C
C
      SUBROUTINE LAX2A(VA,KCOL,KLD,KOP,NEQ,VX,VAX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),KOP(*),VX(*),VAX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LAX2A ','01/02/89')
C
      IF (A1.NE.0D0) THEN
       IF (A2.EQ.0D0) THEN
        CALL LCL2(VAX,NEQ)
       ELSE
        A=A2/A1
        IF (A.NE.1D0) CALL LSC2(VAX,NEQ,A)
       ENDIF
       DO 4 IROW=1,NEQ
       JOP=KOP(IROW)
       DO 5 ICOL=KLD(JOP),KLD(JOP+1)-1
5      VAX(IROW)=DBLE(VAX(IROW))
     *           +DBLE(VA(ICOL))*DBLE(VX(KCOL(ICOL)+IROW))
4      CONTINUE
       IF (A1.NE.1D0) CALL LSC2(VAX,NEQ,A1)
      ELSE
       CALL LSC2(VAX,NEQ,A2)
      ENDIF
C
      END
