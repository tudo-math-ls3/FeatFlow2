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
* XAP9                                                                 *
*                                                                      *
* Purpose  Call AP9                                                    *
*          Allocate KLD and KCOL on DWORK                              *
*                                                                      *
* Subroutines/functions called  AP9, EA00, ZNEW, ZDISP, ZFREE         *
*                                                                      *
* Version from  12/20/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* ELE1     SUBR   EXTERNAL Subroutine - evaluation of basis functions  *
* ELE2     SUBR   for dummy call only                                  *
* For the description of the remaining input parameter see AP9         *
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
* AP9                                                                  *
*                                                                      *
* Purpose  Calculation of the pointer vectors KLD and KCOL             *
*          for a matrix corresponding to given element types           *
*          Storage technique 9                                         *
*                                                                      *
* Subroutines/functions called  NDFL, NDFGL, EA00                      *
*                                                                      *
* Version from  12/20/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* KCOL1    I*4    Auxiliary array                                      *
* KIND     I*4    Auxiliary array                                      *
* NA       I*4    Maximal length of KCOL                               *
* NEQ      I*4    Number of equations                                  *
* ELE1     SUBR   EXTERNAL Subroutine - evaluation of basis functions  *
* ELE2     SUBR   for dummy call only                                  *
* KVERT    I*4    Arrays describing the triangulation                  *
* KMID     I*4                                                         *
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
      SUBROUTINE XAP9(LCOL,LLD,NA,NEQ,ELE1,ELE2)
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
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL ELE1,ELE2
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/
C
      SUB='XAP9'
      IF (ICHECK.GE.997) CALL OTRC('XAP9  ','12/20/89')

C     Put the COMMON-block NVE to a local NVE. Helps to fix a compiler
C     bug with Intel Fortran when compiling with 64-bit integers,
C     128 Bit doubles!

      NVE1 = NVE

C
C *** Determine total number of degrees of freedom
      CALL EA00(ELE1,ELE1,NVE1,ITYPE1)
      NEQ=NDFG(ITYPE1)
      IF (IER.NE.0) GOTO 99999
C
C *** Allocate KLD on DWORK ***
      CALL ZNEW(NEQ+1,-3,LLD,'KLD   ')
      IF (IER.NE.0) GOTO 99999
C
C *** Determine free part of DWORK
      IWMAX0=IWMAX
      CALL ZFREE(3,IFREE)
      IF (IER.NE.0) GOTO 99999
C
      NA=IFREE/3-3
      CALL ZNEW(NA,-3,LCOL,'KCOL  ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NA,-3,LCOL1,'KCOL1 ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NA,-3,LIND,'KIND  ')
      IF (IER.NE.0) GOTO 99999
C
C *** NA contains number of free elements of type I*4 ***
      CALL AP9(KWORK(L(LCOL)),KWORK(L(LCOL1)),
     *         KWORK(L(LIND)),KWORK(L(LLD)),
     *         NA,NEQ,ELE1,ELE2,KWORK(L(LVERT)),KWORK(L(LMID)))
      IF (IER.NE.0) GOTO 99999
C
C *** Release space on DWORK not used for KCOL ***
      CALL ZDISP(NA,LIND,'KIND  ')
      CALL ZDISP(NA,LCOL1,'KCOL1 ')
      CALL ZDISP(NA,LCOL,'KCOL  ')
C
      IWMAX=MAX(IWORK,IWMAX0)
C
      CALL ZDISP(0,LIND,'KIND  ')
      CALL ZDISP(0,LCOL1,'KCOL1 ')
99999 END
C
C
C
      SUBROUTINE AP9(KCOL,KCOL1,KIND,KLD,NA,NEQ,ELE1,ELE2,KVERT,KMID)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=21,NNDER=6,NNVE=4)
      DIMENSION KLD(*),KCOL(*),KCOL1(*),KIND(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),ITYP1(2)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /COAUX2/ DBAS1(NNBAS,NNDER,3),KDFG1(NNBAS,3),
     *                KDFL1(NNBAS,3),IDFL1(3)
      EXTERNAL ELE1,ELE2
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/COAUX2/
C
      SUB='AP9'
      IF (ICHECK.GE.997) CALL OTRC('AP9   ','12/20/89')
C
      IER=0
C
      IF (NA.LT.NEQ) THEN
C *** DWORK exhausted
       WRITE (CPARAM,'(I15)') 1
       CALL WERR(-118,'AP9   ')
       GOTO 99999
      ENDIF
      IFREE=NA-NEQ
      NA=NEQ
C
C *** Initialization of KIND and KCOL1
      DO 10 IEQ=1,NEQ
      KIND(IEQ)=0
10    KCOL1(IEQ)=0
C
C *** Set element numbers by dummy call ***
      CALL EA00(ELE1,ELE1,NVE,ITYP1(1))
      CALL EA00(ELE2,ELE2,NVE,ITYP1(2))
      IF (ITYP1(1).EQ.ITYP1(2)) THEN
       NELE=1
      ELSE
       NELE=2
      ENDIF
C
C *** Determine number of degrees of freedom per element ***
      DO 21 I=1,NELE
      IDFL1(I)=NDFL(ITYP1(I))
      IF (IER.NE.0) GOTO 99999
21    CONTINUE
C
      DO 100 IEL=1,NEL
C
C *** KDFG returns the global degrees of freedom in increasing order
      DO 101 I=1,NELE
      CALL NDFGL(IEL,0,ITYP1(I),KVERT,KMID,KDFG1(1,I),KDFL1(1,I))
      IF (IER.NE.0) GOTO 99999
101   CONTINUE
C
      DO 110 JDOFE=1,IDFL1(1)
      IROW=KDFG1(JDOFE,1)
C
      DO 120 JDOFP=1,IDFL1(NELE)
      JCOL=KDFG1(JDOFP,NELE)
C *** JCOL is the global number of d.o.f. to be inserted into row IROW
C
C *** Look whether entry (IROW,JCOL) is already provided
C
      IF (KCOL1(IROW).EQ.0) THEN
C *** Case 1: No entry in row IROW
       KCOL1(IROW)=JCOL
C
      ELSE
C
       IPOS=IROW
121    IF (KIND(IPOS).EQ.0) THEN
C *** final entry in row IROW up to now
C
        IF (IFREE.LE.0) THEN
C *** DWORK exhausted
         WRITE (CPARAM,'(I15)') IEL
         CALL WERR(-118,'AP9   ')
         GOTO 99999
        ENDIF
C
C *** insert new entry
        NA=NA+1
        IFREE=IFREE-1
        KCOL1(NA)=JCOL
        KIND(IPOS)=NA
        KIND(NA)=0
C
       ELSE
C
        IPOS=KIND(IPOS)
        IF (KCOL1(IPOS).NE.JCOL) GOTO 121
C
       ENDIF
      ENDIF
C
120   CONTINUE
110   CONTINUE
100   CONTINUE
C
C
C *** Collect entries on KCOL1 separately for each row
C
      NA=0
      DO 200 IEQ=1,NEQ
      NA=NA+1
      KLD(IEQ)=NA
      KCOL(NA)=KCOL1(IEQ)
      IPOS=IEQ
C
201   IF (KIND(IPOS).NE.0) THEN
       NA=NA+1
       IPOS=KIND(IPOS)
       KCOL(NA)=KCOL1(IPOS)
       GOTO 201
      ENDIF
C
200   CONTINUE
      KLD(NEQ+1)=NA+1
C
C
C *** Sort entries on KCOL separately for each row
C
      DO 300 IEQ=1,NEQ
C
301   BSORT=.TRUE.
      DO 302 ICOL=KLD(IEQ),KLD(IEQ+1)-2
      IF (KCOL(ICOL).GT.KCOL(ICOL+1)) THEN
       IHELP=KCOL(ICOL)
       KCOL(ICOL)=KCOL(ICOL+1)
       KCOL(ICOL+1)=IHELP
       BSORT=.FALSE.
      ENDIF
302   CONTINUE
      IF (.NOT.BSORT) GOTO 301
C
300   CONTINUE
C
99999 END
