************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT3D  (Release 1.1)               *
*                                                                      *
* Authors: J. Harig, S. Turek                                          *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* XAP3                                                                 *
*                                                                      *
* Purpose  Call AP3                                                    *
*          Allocate KDIAS and KDIA on DWORK                            *
*                                                                      *
* Subroutines/functions called  AP3, ZNEW, ZDISP                       *
*                                                                      *
* Version from  03/07/90                                               *
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
* Purpose  Calculation of the pointer vectors KDIAS and KDIA           *
*          for a matrix corresponding to a given element type          *
*          Storage technique 3/4                                       *
*                                                                      *
* Subroutines/functions called  NDFL, NDFGL                            *
*                                                                      *
* Version from  03/07/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NEQ      I*4    Number of equations                                  *
* ELE      SUBR   EXTERNAL Subroutine - evaluation of basis functions -*
*                 for dummy call only                                  *
* ISYMM    I*4    =1: Matrix symmetric - storage technique 4           *
* KVERT    I*4    Arrays describing the triangulation                  *
* KEDGE    I*4                                                         *
* KAREA    I*4                                                         *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KDIA     I*4    Pointer vector containing column indices             *
* KDIAS    I*4    Pointer vector to diagonal elements                  *
* NDIA     I*4    Number of diagonal rows                              *
* NA       I*4    Number of elements in KDIA                           *
* IER      I*4    Error indicator                                      *
*                 -118  Not enough space for KDIA                      *
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
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL ELE
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/
C
      SUB='XAP3'
      IF (ICHECK.GE.997) CALL OTRC('XAP3  ','03/07/90')
C
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
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
     *         KWORK(L(LVERT)),KWORK(L(LEDGE)),KWORK(L(LAREA)),
     *         KWORK(L(LDIAH)))
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
      SUBROUTINE AP3(KDIA,KDIAS,NDIA,NA,NEQ,ELE,ISYMM,
     *               KVERT,KEDGE,KAREA,KDIAH)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=27,NNVE=8,NNEE=12,NNAE=6)
      DIMENSION KDIA(*),KDIAS(*),KDIAH(*)
      DIMENSION KVERT(NNVE,*),KEDGE(NNEE,*),KAREA(NNAE,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      EXTERNAL ELE
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='AP3'
      IF (ICHECK.GE.997) CALL OTRC('AP3   ','03/07/90')
C
      IER=0
      JDFL=1
      BSYMM=ISYMM.EQ.1
C
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
C
      IDFL=NDFL(IELTYP)
      IF (IER.NE.0) GOTO 99999
C
      IF (BSYMM) THEN
       IOFSET=1
      ELSE
       IOFSET=NEQ
      ENDIF
C
      DO 100 IEL=1,NEL
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.NE.0) GOTO 99999
C
      DO 110 IDFL1=1,IDFL
      IDFG1=KDFG(IDFL1)
      IF (BSYMM) JDFL=IDFL1
      DO 120 IDFL2=JDFL,IDFL
      IDFG2=KDFG(IDFL2)
      IDFG=IDFG2-IDFG1+IOFSET
      KDIAH(IDFG)=1
120   CONTINUE
110   CONTINUE
C
100   CONTINUE
C
      NDIA=1
      KDIA(NDIA)=0
      IF (BSYMM) THEN
        DO 130 IEQ=2,NEQ
          IF (KDIAH(IEQ).NE.0) THEN
            NDIA=NDIA+1
            KDIA(NDIA)=IEQ-1
          ENDIF
130     CONTINUE
      ELSE
        DO 140 IEQ=1,NEQ-1
          IF (KDIAH(IEQ).NE.0) THEN
            NDIA=NDIA+1
            KDIA(NDIA)=IEQ-NEQ
          ENDIF
140     CONTINUE
        DO 150 IEQ=NEQ+1,2*NEQ-1
          IF (KDIAH(IEQ).NE.0) THEN
            NDIA=NDIA+1
            KDIA(NDIA)=IEQ-NEQ
          ENDIF
150     CONTINUE
      ENDIF
C
      KDIAS(1)=1
      NA=NEQ
C$DIR SCALAR
      DO 160 IDIA=2,NDIA
      KDIAS(IDIA)=NA+1
      NA=NA+NEQ-ABS(KDIA(IDIA))
160   CONTINUE
      KDIAS(NDIA+1)=NA+1
C
99999 END
