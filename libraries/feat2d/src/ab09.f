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
* XAB09                                                                *
*                                                                      *
* Purpose  Allocation of the blocks in vector DA , KOFF and COECON     *
*          Call of AB09                                                *
*                                                                      *
* Subroutines/functions called  AB09, ZNEW, ZDISP, ZLEN, ZTYPE, ZCLEAR *
*                                                                      *
* Version from  12/04/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LA       I*4    Vector of numbers of the arrays                      *
*                 New arrays are allocated on DWORK for LA(IBLOC)=0    *
* LCOL     I*4    Numbers of the pointer arrays KCOL and KLD           *
* LLD      I*4                                                         *
* ICLEAR   I*4    ICLEAR=1  all matrices are cleared                   *
* ARR      C*6    Names of blocks (for error messages only)            *
* BSNGL    L*4    NBLOC logical values                                 *
*                 .TRUE. means that corresponding array is converted   *
*                 to single precision after completion                 *
* For the description of the remaining parameters see AB09             *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* LA       I*4    Vector of numbers of arrays on DA                    *
* IER      I*4    Error indicator                                      *
*                 -114 Type of at least one array is not double prec.  *
*                 -115 Length of at least one array is < NA            *
*                                                                      *
************************************************************************
*                                                                      *
* AB09                                                                 *
*                                                                      *
* Purpose  Calculation of NBLOC matrices corresponding to              *
*          bilinear forms  a(v,u)                                      *
*          Quadrilateral elements - bilinear transformation            *
*          Mixed trial and test functions                              *
*          Storage technique 9                                         *
*          Each matrix makes use of the same pointer vectors           *
*          Double precision version                                    *
*                                                                      *
* Subroutines/functions called NDFL, NDFGL, CB2Q                       *
*                                                                      *
* Version from  12/04/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* KCOL     I*4    Pointer vectors for each matrix on DA corresponding  *
* KLD      I*4    to storage technique 9 (calculated by AP9)           *
* NA       I*4    Number of elements per matrix                        *
* NBLOC    I*4    Number of matrices stored on DA                      *
* KOFF     I*4    Matrix IBLOC starts at position KOFF(IBLOC)+1 on DA  *
* KVERT    I*4    Arrays describing the triangulation                  *
* KMID     I*4                                                         *
* DCORVG   R*8                                                         *
* ELE1     SUBR   EXTERNAL SUBROUTINE - values of basis functions      *
*                 corresponding to first argument of a(.,.)            *
* ELE2     SUBR   EXTERNAL SUBROUTINE - values of basis functions      *
*                 corresponding to second argument of a(.,.)           *
* ELE3     SUBR   EXTERNAL SUBROUTINE - values of basis functions      *
*                 corresponding to nonlinearity in COEFF               *
* COEFF    R*8    EXTERNAL FUNCTION - coefficients of the bilinear form*
* BCON     LOG    BCON(IBLOC)=.TRUE. means constant coefficients in    *
*                 matrix IBLOC                                         *
* KAB      I*4    Pairs of multiindices occuring in the bilinear forms *
*                 specified separately for each matrix                 *
* KABN     I*4    Number of additive terms in each bilinear form       *
* ICUB     I*4    Number of cubature formula in CB2Q                   *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DA       R*8    Calculated matrices                                  *
* IER      I*4    Error indicator                                      *
*                 -116 Wrong value in array KAB                        *
*                                                                      *
************************************************************************
C
      SUBROUTINE XAB09(LA,LCOL,LLD,NA,NBLOC,ICLEAR,
     *                 ELE1,ELE2,ELE3,COEFF,BCON,KAB,KABN,
     *                 ICUB,ARR,BSNGL)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6
C
      PARAMETER (NNARR=299,NNAB=21,NNDER=6)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION LA(*),KAB(2,NNAB,*),KABN(*),BCON(*),ARR(*),BSNGL(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL COEFF,ELE1,ELE2,ELE3
      SAVE /TRIAA/,/OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='XAB09'
      IF (ICHECK.GE.997) CALL OTRC('XAB09 ','12/04/90')
      IER=0
C
      DO 1 IBLOC=1,NBLOC
      IF (LA(IBLOC).EQ.0) THEN
       CALL ZNEW(NA,1,LA(IBLOC),ARR(IBLOC))
       IF (IER.NE.0) GOTO 99999
      ELSE
C ***  Check input parameter
       CALL ZTYPE(LA(IBLOC),ITYPE)
       IF (ITYPE.NE.1) THEN
        WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
        CALL WERR(-114,'XAB09 ')
        GOTO 99999
       ENDIF
       CALL ZLEN(LA(IBLOC),ILEN)
       IF (ILEN.LT.NA) THEN
        WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
        CALL WERR(-115,'XAB09 ')
        GOTO 99999
       ENDIF
       IF (ICLEAR.EQ.1) CALL ZCLEAR(LA(IBLOC),ARR(IBLOC))
      ENDIF
1     CONTINUE
C
      CALL ZNEW(NBLOC,-3,LOFF,'KOFF  ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NBLOC*NNDER**2,1,LOECON,'COECON')
      IF (IER.NE.0) GOTO 99999
C
      DO 2 IBLOC=1,NBLOC
      KWORK(L(LOFF)+IBLOC-1)=L(LA(IBLOC))-1
2     CONTINUE
C
      CALL AB09(DWORK(1),KWORK(L(LCOL)),KWORK(L(LLD)),NA,NBLOC,
     *          KWORK(L(LOFF)),KWORK(L(LVERT)),KWORK(L(LMID)),
     *          DWORK(L(LCORVG)),ELE1,ELE2,
     *          ELE3,COEFF,BCON,DWORK(L(LOECON)),KAB,KABN,ICUB)
      IF (IER.NE.0) GOTO 99999
C
      CALL ZDISP(0,LOECON,'COECON')
      CALL ZDISP(0,LOFF,'KOFF  ')
C
      DO 3 IBLOC=NBLOC,1,-1
      IF (BSNGL(IBLOC)) THEN
       CALL ZCTYPE(2,LA(IBLOC),ARR(IBLOC))
       IF (IER.NE.0) GOTO 99999
      ENDIF
3     CONTINUE
C
99999 END
C
C
C
      SUBROUTINE AB09(DA,KCOL,KLD,NA,NBLOC,KOFF,KVERT,KMID,
     *                DCORVG,ELE1,ELE2,ELE3,COEFF,
     *                BCON,COECON,KAB,KABN,ICUB)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21)
      DIMENSION KCOL(*),KVERT(NNVE,*),KMID(NNVE,*),KOFF(*)
      DIMENSION KLD(*),DA(*),DCORVG(2,*)
      DIMENSION BCON(*),KAB(2,NNAB,*),KABN(*),COECON(NNDER,NNDER,*)
      DIMENSION BDER1(NNDER,3),ITYP1(3),KENTRY(NNBAS,NNBAS),BELE(3)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX2/ DBAS1(NNBAS,NNDER,3),KDFG1(NNBAS,3),
     *                KDFL1(NNBAS,3),IDFL1(3)
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/ELEM/,/TRIAD/,/CUB/,/COAUX2/
C
      SUB='AB09'
      IF (ICHECK.GE.997) CALL OTRC('AB09  ','12/04/90')
C
C *** Preparation - evaluation of parameters
      IER=0
C *** Which derivatives of basis functions are needed?
      DO 1 IDER=1,NNDER
      DO 1 J=1,2
1     BDER1(IDER,J)=.FALSE.
      DO 2 IBLOC=1,NBLOC
      DO 3 IABN=1,KABN(IBLOC)
      DO 3 J=1,2
      IAB=KAB(J,IABN,IBLOC)
      IF (IAB.LE.0.OR.IAB.GT.NNDER) THEN
       WRITE (CPARAM,'(I15)') IBLOC
       CALL WERR(-116,'AB09  ')
       GOTO 99999
      ENDIF
3     BDER1(IAB,J)=.TRUE.
2     CONTINUE
C *** Dummy call of ELEn sets number of element
      DO 4 I=1,3
      ITYP1(I)=-1
      GOTO (41,42,43) I
41    CALL ELE1(0D0,0D0,ITYP1(1))
      GOTO 44
42    CALL ELE2(0D0,0D0,ITYP1(2))
      GOTO 44
43    CALL ELE3(0D0,0D0,ITYP1(3))
44    IF (IER.NE.0) GOTO 99999
      IDFL1(I)=NDFL(ITYP1(I))
4     CONTINUE
C
      BELE(1)=.TRUE.
      BELE(2)=ITYP1(1).NE.ITYP1(2)
      IF (BELE(2)) THEN
       NELE=2
      ELSE
       NELE=1
      ENDIF
      BELE(3)=ITYP1(3).NE.ITYP1(1).AND.ITYP1(3).NE.ITYP1(2)
C *** Evaluation of cubature points and weights
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C
      BFIRST=.TRUE.
C *** Dummy call of COEFF for nonlinear problems
C *** COEFF must set BDER1(IDER,.)=.TRUE. if derivative IDER is needed
      AUX=COEFF(0D0,0D0,-1,-1,0,BFIRST)
C
      BCON0=.TRUE.
      DO 5 IBLOC=1,NBLOC
      IF (BCON(IBLOC)) THEN
       DO 6 I=1,KABN(IBLOC)
       IA=KAB(1,I,IBLOC)
       IB=KAB(2,I,IBLOC)
6      COECON(IA,IB,IBLOC)=COEFF(0D0,0D0,IA,IB,IBLOC,BFIRST)
      ELSE
       BCON0=.FALSE.
      ENDIF
5     CONTINUE
************************************************************************
C *** Calculation of the matrix - storage technique 9
************************************************************************
C *** Dummy call - ELEn may save arithmetic operations
      ICUBP=ICUB
      DO 10 I=1,NELE
      IF (.NOT.BELE(I)) GOTO 10
      DO 11 IDER=1,NNDER
11    BDER(IDER)=BDER1(IDER,I)
      GOTO (12,13,14) I
12    CALL ELE1(0D0,0D0,-2)
      GOTO 10
13    CALL ELE2(0D0,0D0,-2)
      GOTO 10
14    CALL ELE3(0D0,0D0,-2)
10    CONTINUE
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
C
      DO 101 I=1,3
      IF (.NOT.BELE(I)) GOTO 101
      CALL NDFGL(IEL,1,ITYP1(I),KVERT,KMID,KDFG1(1,I),KDFL1(1,I))
101   CONTINUE
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL1(1)
      JCOL0=KLD(KDFG1(JDOFE,1))
      DO 110 IDOFE=1,IDFL1(NELE)
      IDFG=KDFG1(IDOFE,NELE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOL(JCOL).EQ.IDFG) GOTO 111
112   CONTINUE
111   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices
      DO 120 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
120   CONTINUE
C
      DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
      DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
      DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
      DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
C
C *** Loop over all cubature points
      DO 200 ICUBP = 1, NCUBP
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
C *** Jacobian of the bilinear mapping onto the reference element
      DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
      DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
      DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
      DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
      DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
      OM=DOMEGA(ICUBP)*DETJ
C *** ELEn needs the information ICUBP because of preceeding
C *** dummy call using IPAR = -2
      DO 20 I=1,3
      IF (.NOT.BELE(I)) GOTO 20
      DO 21 IDER=1,NNDER
21    BDER(IDER)=BDER1(IDER,I)
      GOTO (22,23,24) I
22    CALL ELE1(XI1,XI2,-3)
      GOTO 25
23    CALL ELE2(XI1,XI2,-3)
      GOTO 25
24    CALL ELE3(XI1,XI2,-3)
25    DO 27 IDER=1,NNDER
      IF (.NOT.BDER1(IDER,I)) GOTO 27
      DO 26 IBAS=1,IDFL1(I)
26    DBAS1(IBAS,IDER,I)=DBAS(IBAS,IDER)
27    CONTINUE
20    CONTINUE
C
      IF (.NOT.BCON0) THEN
       XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *   +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
       YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1+0.5D0*
     *          (DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2
      ENDIF
C
C *** Summing up over all pairs of multiindices
      BFIRST=.TRUE.
      DO 300 IBLOC=1,NBLOC
C
      DO 301 IALBET=1,KABN(IBLOC)
      IA=KAB(1,IALBET,IBLOC)
      IB=KAB(2,IALBET,IBLOC)
      IF (.NOT.BCON(IBLOC)) THEN
       AUX=COEFF(XX,YY,IA,IB,IBLOC,BFIRST)*OM
      ELSE
       AUX=COECON(IA,IB,IBLOC)*OM
      ENDIF
C
      DO 310 JDOFE=1,IDFL1(1)
      DB=DBAS1(KDFL1(JDOFE,1),IA,1)
      DO 320 IDOFE=1,IDFL1(NELE)
      JCOLB=KENTRY(JDOFE,IDOFE)+KOFF(IBLOC)
      DA(JCOLB)=DA(JCOLB)+DB*DBAS1(KDFL1(IDOFE,NELE),IB,NELE)*AUX
320   CONTINUE
310   CONTINUE
301   CONTINUE
C
      BFIRST=.FALSE.
300   CONTINUE
C
200   CONTINUE
100   CONTINUE
C
99999 END
