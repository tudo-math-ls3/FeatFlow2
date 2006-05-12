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
* XAP7                                                                 *
*                                                                      *
* Purpose  Call AP7                                                    *
*          Allocate KLD and KCOL on DWORK                              *
*                                                                      *
* Subroutines/functions called  AP7, ZNEW, ZDISP                       *
*                                                                      *
* Version from  12/20/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* ELE      SUBR   EXTERNAL Subroutine - evaluation of basis functions  *
*                 for dummy call only                                  *
* For the description of the remaining input parameter see AP7         *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* LCOL     I*4    Number of KCOL                                       *
* LLD      I*4    Number of KLD                                        *
* NA       I*4    Length of KCOL (and VA (DA))                         *
* NEQ      I*4    Number of equations                                  *
*                                                                      *
************************************************************************
*                                                                      *
* AP7                                                                  *
*                                                                      *
* Purpose  Calculation of the pointer vectors KLD and KCOL             *
*          for a matrix corresponding to a given element type          *
*          Storage technique 7/8                                       *
*                                                                      *
* Subroutines/functions called  NDFL, NDFGL                            *
*                                                                      *
* Version from  12/20/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NA       I*4    Maximal length of KCOL                               *
* NEQ      I*4    Number of equations                                  *
* ELE      SUBR   EXTERNAL Subroutine - evaluation of basis functions -*
*                 for dummy call only                                  *
* ISYMM    I*4    =1: Matrix symmetric - storage technique 8           *
* NEROW    I*4    Assumed maximum number of nonzero entries in a row   *
* KVERT    I*4    Arrays describing the triangulation                  *
* KEDGE    I*4                                                         *
* KAREA    I*4                                                         *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KCOL     I*4    Pointer vector containing column indices             *
* KLD      I*4    Pointer vector to diagonal elements                  *
* NA       I*4    Number of elements in KCOL                           *
* IER      I*4    Error indicator                                      *
*                 -118  Not enough space for KCOL                      *
*                       Error occured on element IEL                   *
*                                                                      *
************************************************************************
C
      SUBROUTINE XAP7(LCOL,LLD,NA,NEQ,ELE,ISYMM,NEROW)
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
      SUB='XAP7'
      IF (ICHECK.GE.997) CALL OTRC('XAP7  ','12/20/89')
C
C *** Determine total number of degrees of freedom
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      NEQ=NDFG(IELTYP)
      IF (IER.NE.0) GOTO 99999
C
C *** Allocate KLD on DWORK ***
      CALL ZNEW(NEQ+1,-3,LLD,'KLD   ')
      IF (IER.NE.0) GOTO 99999
C
C *** Reserve free part of DWORK for KCOL
      NA=0
      CALL ZNEW(NA,-3,LCOL,'KCOL  ')
      IF (IER.NE.0) GOTO 99999
C
C *** Do not change input parameter
      IEROW=NEROW
C *** NA contains number of free elements of type I*4 ***
      IF (NA.LT.NEQ*IEROW) IEROW=NA/NEQ
      CALL AP7(KWORK(L(LCOL)),KWORK(L(LLD)),NA,NEQ,ELE,ISYMM,
     *         IEROW,KWORK(L(LVERT)),KWORK(L(LEDGE)),KWORK(L(LAREA)))
      IF (IER.NE.0) GOTO 99999
C
C *** Release space on DWORK not used for KCOL ***
      CALL ZDISP(NA,LCOL,'KCOL  ')
C
99999 END
C
C
C
      SUBROUTINE AP7(KCOL,KLD,NA,NEQ,ELE,ISYMM,NEROW,KVERT,KEDGE,KAREA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=27,NNVE=8,NNEE=12,NNAE=6)
      DIMENSION KCOL(*),KLD(*)
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
      SUB='AP7'
      IF (ICHECK.GE.997) CALL OTRC('AP7   ','03/01/90')
C
      IER=0
      BSYMM=ISYMM.EQ.1
C *** Do not change input parameter
      IEROW=NEROW
      IF (IEROW.LE.0) IEROW=27
C *** Standard number of entries for trilinear elements
C
      NA1=IEROW*NEQ
C     IFREE: free number of entries in KCOL
C            for the moment corresponding to the choice of NEROW
      IFREE=(IEROW-1)*NEQ
C *** Initialization of KCOL and KLD ***
      DO 10 IA1=1,NA1
10    KCOL(IA1)=0
      DO 20 IEQ=1,NEQ
      KLD(IEQ)=(IEQ-1)*IEROW+1
20    KCOL((IEQ-1)*IEROW+1)=IEQ
      KLD(NEQ+1)=NA1+1
C
C *** Set element number by dummy call ***
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
C *** Determine number of degrees of freedom per element ***
      IDFL=NDFL(IELTYP)
      IF (IER.NE.0) GOTO 99999
C
      JDOFE1=1
      DO 100 IEL=1,NEL
C
C *** KDFG returns the global degrees of freedom in increasing order
      CALL NDFGL(IEL,0,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.NE.0) GOTO 99999
C
      DO 110 JDOFE=1,IDFL
      IROW=KDFG(JDOFE)
      IF (BSYMM) THEN
       IF (IROW.EQ.NEQ) GOTO 110
       JDOFE1=JDOFE
      ENDIF
      IPOS2=KLD(IROW+1)-1
      JPOSP=KLD(IROW)+1
C
      DO 120 JDOFP=JDOFE1,IDFL
      IF (JDOFE.EQ.JDOFP) GOTO 120
      JCOL=KDFG(JDOFP)
C *** JCOL must be inserted
      IF (JPOSP.GT.IPOS2) THEN
       JINSP=JPOSP
       GOTO 140
      ENDIF
C
      DO 130 JINSP=JPOSP,IPOS2
      JCOL0=KCOL(JINSP)
      IF (JCOL0.EQ.0) THEN
C *** Insert as last entry in row IROW
       JPOSP=JINSP+1
       KCOL(JINSP)=JCOL
       GOTO 190
      ENDIF
      IF (JCOL-JCOL0) 132,120,130
130   CONTINUE
C
132   IF (KCOL(IPOS2).EQ.0) THEN
C *** Insert JCOL at position JINSP
      JPOSP=JINSP+1
      DO 133 JP=JPOSP,IPOS2
133   IF (KCOL(JP).EQ.0) GOTO 134
134   DO 135 MM=JP,JPOSP,-1
C *** Shift in row IROW needed
135   KCOL(MM)=KCOL(MM-1)
      KCOL(JINSP)=JCOL
      GOTO 190
      ENDIF
C
140   IF (IFREE.EQ.0) THEN
C ***  DWORK exhausted
       IF (NA1.EQ.NA) THEN
        WRITE (CPARAM,'(I15)') IEL
        CALL WERR(-118,'AP7   ')
        GOTO 99999
       ENDIF
C ***  Increase length of KCOL taking another NEQ elements
       NA2=MIN(NA1+NEQ,NA)
       KLD(NEQ+1)=NA2+1
       IFREE=NA2-NA1
C ***  Initialization of the new elements
       DO 141 I=NA1+1,NA2
141    KCOL(I)=0
       NA1=NA2
       IF (IROW.EQ.NEQ) IPOS2=NA2
       WRITE (CPARAM,'(I15)') NA2
       CALL OMSG(20,'AP7   ')
      ENDIF
C
C *** Global shift needed - look for free position
      IF (IROW.EQ.NEQ.AND.KCOL(NA1).NE.0) GOTO 150
      DO 142 JP=IROW,NEQ
      I2=KLD(JP+1)-1
142   IF (KCOL(I2).EQ.0) GOTO 143
C *** No free position to the right
      GOTO 150
143   DO 144 MM=IROW+1,JP
144   KLD(MM)=KLD(MM)+1
      IPOS2=IPOS2+1
      DO 145 JP=I2-1,JINSP,-1
145   KCOL(JP+1)=KCOL(JP)
      KCOL(JINSP)=JCOL
      JPOSP=JINSP+1
      GOTO 190
C
C *** Look for free position to the left
C *** It must be found since IFREE > 0
150   DO 151 JP=IROW,2,-1
      I2=KLD(JP)-1
151   IF (KCOL(I2).LE.0) GOTO 152
152   DO 153 MM=JP,IROW
153   KLD(MM)=KLD(MM)-1
      DO 154 JP=I2,JINSP-2
154   KCOL(JP)=KCOL(JP+1)
      KCOL(JINSP-1)=JCOL
      JPOSP=JINSP
C
190   IFREE=IFREE-1
120   CONTINUE
110   CONTINUE
100   CONTINUE
C
C *** Compress KCOL ***
C
      NA=0
      DO 200 IEQ=1,NEQ
      IPOS1=KLD(IEQ)
      DO 201 IPOS2=KLD(IEQ+1)-1,IPOS1,-1
201   IF (KCOL(IPOS2).NE.0) GOTO 202
202   IF (IEQ.EQ.1) GOTO 204
      DO 203 JP=IPOS1,IPOS2
203   KCOL(JP-IPOS1+NA+1)=KCOL(JP)
204   KLD(IEQ)=NA+1
200   NA=NA+IPOS2-IPOS1+1
      KLD(NEQ+1)=NA+1
C
99999 END
