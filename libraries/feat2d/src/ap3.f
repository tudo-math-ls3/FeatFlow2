************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.0)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* XAP3                                                                 *
*                                                                      *
* Purpose  Call AP3                                                    *
*                                                                      *
* Subroutines/functions called  AP3, EA00, ZNEW, ZDISP                 *
*                                                                      *
* Version from  11/19/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* ELE      SUBR   EXTERNAL Subroutine - evaluation of basis functions  *
*                 for dummy call only                                  *
* For the description of the remaining input parameter see AP3         *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* LDIA     I*4    Number of KDIA                                       *
* LDIAS    I*4    Number of KDIAS                                      *
*                                                                      *
************************************************************************
*                                                                      *
* AP3                                                                  *
*                                                                      *
* Purpose  Calculation of the pointer vectors                          *
*          for a matrix corresponding to a given element type          *
*                                                                      *
* Subroutines/functions called  NDFL, NDFGL, EA00                      *
*                                                                      *
* Version from  11/19/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NEQ      I*4    Number of equations                                  *
* ELE      SUBR   EXTERNAL Subroutine - evaluation of basis functions -*
*                 for dummy call only                                  *
* ISYMM    I*4    1  matrix symmetric - storage technique 4            *
* KVERT    I*4    Arrays describing the triangulation                  *
* KMID     I*4                                                         *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KDIA     I*4    Pointer vector containing column indices             *
* KDIAS    I*4    Pointer vector to diagonal elements                  *
* NDIA     I*4    Number of diagonal rows                              *
* NEQ      I*4    Number of equations                                  *
* NA       I*4    Number of elements in KDIA                           *
* IER      I*4    Error indicator                                      *
*                 -118  Not enough space for KCOL                      *
*                       Error occured on element IEL                   *
*                                                                      *
************************************************************************
C
      SUBROUTINE XAP3(LDIA,LDIAS,NDIA,NA,NEQ,ELE,ISYMM)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL ELE
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/
C
      SUB='XAP3'
      IF (ICHECK.GE.997) CALL OTRC('XAP3  ','11/19/90')
      IER=0

C     Put the COMMON-block NVE to a local NVE. Helps to fix a compiler
C     bug with Intel Fortran when compiling with 64-bit integers,
C     128 Bit doubles!

      NVE1 = NVE

C
C *** Determine total number of degrees of freedom
      CALL EA00(ELE,ELE,NVE1,IELTYP)
      NEQ=NDFG(IELTYP)
      IF (IER.NE.0) GOTO 99999
C
      IF (ISYMM.EQ.1) THEN
       NEQ0=NEQ
      ELSE
       NEQ0=2*NEQ-1
      ENDIF
C
      CALL ZNEW(NEQ0  ,-3,LDIA ,'KDIA  ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NEQ0+1, 3,LDIAS,'KDIAS ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NEQ0  , 3,LDIAH,'KDIAH ')
      IF (IER.NE.0) GOTO 99999
C
      CALL AP3(KWORK(L(LDIA)),KWORK(L(LDIAS)),NDIA,NA,NEQ,ELE,ISYMM,
     *         KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LDIAH)))
      IF (IER.NE.0) GOTO 99999
C
      CALL ZDISP(0     ,LDIAH,'KDIAH ')
      CALL ZDISP(NDIA+1,LDIAS,'KDIAS ')
      CALL ZDISP(NDIA  ,LDIA ,'KDIA  ')
C
99999 END
C
C
C
      SUBROUTINE AP3(KDIA,KDIAS,NDIA,NA,NEQ,ELE,ISYMM,KVERT,KMID,KDIAH)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=21,NNVE=4)
      DIMENSION KDIAH(*),KDIAS(*),KDIA(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KDFG(NNBAS),KDFL(NNBAS)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      EXTERNAL ELE
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='AP3'
      IF (ICHECK.GE.997) CALL OTRC('AP3   ','11/19/90')
C
      IER=0
      JDFL=1
      BSYMM=ISYMM.EQ.1

C     Put the COMMON-block NVE to a local NVE. Helps to fix a compiler
C     bug with Intel Fortran when compiling with 64-bit integers,
C     128 Bit doubles!

      NVE1 = NVE

C
C *** Set element number by dummy call ***
      CALL EA00(ELE,ELE,NVE1,IELTYP)
C
C *** Determine number of degrees of freedom per element ***
      IDFL=NDFL(IELTYP)
      IF (IER.NE.0) GOTO 99999
C
      IF (BSYMM) THEN
       IOFSET=1
      ELSE
       IOFSET=NEQ
      ENDIF
C
C *** Determine differences between all node numbers
C *** Each difference corresponds to a diagonal
      DO 100 IEL=1,NEL
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.NE.0) GOTO 99999
C
      DO 110 IDFL1=1,IDFL
      IDFG1=KDFG(IDFL1)
      IF (BSYMM) JDFL=IDFL1
C
      DO 120 IDFL2=JDFL,IDFL
      IDFG2=KDFG(IDFL2)
      IDFG=IDFG2-IDFG1+IOFSET
      KDIAH(IDFG)=1
120   CONTINUE
C
110   CONTINUE
C
100   CONTINUE
C
C *** Store the resulting differences in KDIA (main diagonal
C *** first, then the other diagonals from left to right)
      NDIA=1
      KDIA(NDIA)=0
      IF (BSYMM) THEN
        DO 130 IEQ=2,NEQ
          IF (KDIAH(IEQ) .NE. 0) THEN
            NDIA=NDIA+1
            KDIA(NDIA)=IEQ-1
          ENDIF
130     CONTINUE
      ELSE
        DO 140 IEQ=1,NEQ-1
          IF (KDIAH(IEQ) .NE. 0) THEN
            NDIA=NDIA+1
            KDIA(NDIA)=IEQ-NEQ
          ENDIF
140     CONTINUE
        DO 150 IEQ=NEQ+1,2*NEQ-1
          IF (KDIAH(IEQ) .NE. 0) THEN
            NDIA=NDIA+1
            KDIA(NDIA)=IEQ-NEQ
          ENDIF
150     CONTINUE
      ENDIF
C
      KDIAS(1)=1
      NA=NEQ
      DO 160 IDIA=2,NDIA
      KDIAS(IDIA)=NA+1
      NA=NA+NEQ-ABS(KDIA(IDIA))
160   CONTINUE
      KDIAS(NDIA+1)=NA+1
C
99999 END
