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
* XAA03                                                                *
*                                                                      *
* Purpose  Allocation of the blocks in vector DA , KOFF and COECON     *
*          Call of AA03                                                *
*                                                                      *
* Subroutines/functions called  AA03, ZNEW, ZDISP, ZLEN, ZTYPE, ZCLEAR *
*                                                                      *
* Version from  12/04/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LA       I*4    Vector of numbers of the arrays                      *
*                 New arrays are allocated on DWORK for LA(IBLOC)=0    *
* LDIA     I*4    Numbers of the pointer arrays KDIA and KDIAS         *
* LDIAS    I*4    calculated by AP3                                    *
* ICLEAR   I*4    ICLEAR=1  all matrices are cleared                   *
* ARR      C*6    Names of blocks (for error messages only)            *
* BSNGL    L*4    NBLOC logical values                                 *
*                 .TRUE. means that corresponding array is converted   *
*                 to single precision after completion                 *
* For the description of the remaining parameters see AA03             *
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
* AA03                                                                 *
*                                                                      *
* Purpose  Calculation of NBLOC matrices corresponding to              *
*          bilinear forms  a(v,u)                                      *
*          Triangular elements - affine linear transformation          *
*          Same trial and test functions                               *
*          Storage technique 3                                         *
*          Each matrix makes use of the same pointer vectors           *
*          Double precision version                                    *
*                                                                      *
* Subroutines/functions called NDFL, NDFGL, CB2T                       *
*                                                                      *
* Version from  12/04/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* KDIA     I*4    Pointer vectors for each matrix on DA corresponding  *
* KDIAS    I*4    to storage technique 3 (calculated by AP3)           *
* NA       I*4    Number of elements per matrix                        *
* NDIA     I*4    Number of diagonal rows                              *
* NEQ      I*4    Number of equations per matrix                       *
* NBLOC    I*4    Number of matrices stored on DA                      *
* KOFF     I*4    Matrix IBLOC starts at position KOFF(IBLOC)+1 on DA  *
* KVERT    I*4    Arrays describing the triangulation                  *
* KMID     I*4                                                         *
* DCORVG   R*8                                                         *
* ELE      SUBR   EXTERNAL SUBROUTINE - values of basis functions      *
* COEFF    R*8    EXTERNAL FUNCTION - coefficients of the bilinear form*
* BCON     LOG    BCON(IBLOC)=.TRUE. means constant coefficients in    *
*                 matrix IBLOC                                         *
* KAB      I*4    Pairs of multiindices occuring in the bilinear forms *
*                 specified separately for each matrix                 *
* KABN     I*4    Number of additive terms in each bilinear form       *
* ICUB     I*4    Number of cubature formula in CB2T                   *
* ISYMM    I*4    0  storage technique 3                               *
*                 1  storage technique 4                               *
*                 2  storage technique 3 but for symmetric matrices    *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DA       R*8    Calculated matrices                                  *
* IER      I*4    Error indicator                                      *
*                 -116 Wrong value in array KAB                        *
*                                                                      *
************************************************************************
C
      SUBROUTINE XAA03(LA,LDIA,LDIAS,NDIA,NA,NEQ,NBLOC,ICLEAR,ELE,
     *                 COEFF,BCON,KAB,KABN,ICUB,ISYMM,ARR,BSNGL)
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
      EXTERNAL COEFF,ELE
      SAVE /TRIAA/,/OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='XAA03'
      IF (ICHECK.GE.997) CALL OTRC('XAA03 ','12/04/90')
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
        CALL WERR(-114,'XAA03 ')
        GOTO 99999
       ENDIF
       CALL ZLEN(LA(IBLOC),ILEN)
       IF (ILEN.LT.NA) THEN
        WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
        CALL WERR(-115,'XAA03 ')
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
      CALL AA03(DWORK(1),KWORK(L(LDIA)),KWORK(L(LDIAS)),NDIA,NA,NEQ,
     *          NBLOC,KWORK(L(LOFF)),KWORK(L(LVERT)),KWORK(L(LMID)),
     *          DWORK(L(LCORVG)),ELE,COEFF,
     *          BCON,DWORK(L(LOECON)),KAB,KABN,ICUB,ISYMM)
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
      SUBROUTINE AA03(DA,KDIA,KDIAS,NDIA,NA,NEQ,NBLOC,KOFF,KVERT,KMID,
     *                DCORVG,ELE,COEFF,BCON,COECON,
     *                KAB,KABN,ICUB,ISYMM)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21)
      DIMENSION KDIA(*),KVERT(NNVE,*),KMID(NNVE,*),KOFF(*)
      DIMENSION KDIAS(*),KDFG(NNBAS),KDFL(NNBAS),DA(*),DCORVG(2,*)
      DIMENSION BCON(*),KAB(2,NNAB,*),KABN(*),COECON(NNDER,NNDER,*)
      DIMENSION KENTRY(NNBAS,NNBAS)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/ELEM/,/TRIAD/,/CUB/,/COAUX1/
	DATA IOFF/0/
C
      SUB='AA03'
      IF (ICHECK.GE.997) CALL OTRC('AA03  ','12/04/90')
C
C *** Preparation - evaluation of parameters
      IER=0
      BSYMM=ISYMM.GE.1
      JDFG=1
C *** Which derivatives of basis functions are needed?
      DO 1 I = 1, NNDER
1     BDER(I)=.FALSE.
      DO 2 IBLOC=1,NBLOC
      DO 3 I=1,KABN(IBLOC)
      DO 3 J=1,2
      I1=KAB(J,I,IBLOC)
      IF (I1.LE.0.OR.I1.GT.NNDER) THEN
       WRITE (CPARAM,'(I15)') IBLOC
       CALL WERR(-116,'AA03  ')
       GOTO 99999
      ENDIF
3     BDER(I1)=.TRUE.
2     CONTINUE
C *** Dummy call of ELE sets number of element
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2T(ICUB)
      IF (IER.NE.0) GOTO 99999
C
      BFIRST=.TRUE.
C *** Dummy call of COEFF for nonlinear problems
C *** COEFF must set BDER(IDER)=.TRUE. if derivative IDER is needed
      AUX=COEFF(0D0,0D0,-1,-1,0,BFIRST)
      BCON0=.TRUE.
C
      DO 4 IBLOC=1,NBLOC
      IF (BCON(IBLOC)) THEN
       DO 5 I=1,KABN(IBLOC)
       IA=KAB(1,I,IBLOC)
       IB=KAB(2,I,IBLOC)
5      COECON(IA,IB,IBLOC)=COEFF(0D0,0D0,IA,IB,IBLOC,BFIRST)
      ELSE
       BCON0=.FALSE.
      ENDIF
4     CONTINUE
************************************************************************
C *** Calculation of the matrix - storage technique 3 or 4
************************************************************************
C *** Dummy call - ELE may save arithmetic operations
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
C
C *** Loop over all elements
      JDOFE1=1
      DO 100 IEL=1,NEL
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine entry positions in matrix
      DO 110 I=1,IDFL
      IDFG1=KDFG(I)
      IF (BSYMM) JDFG=I
      DO 111 J=JDFG,IDFL
      IDFG2=KDFG(J)
      IDFG=IDFG2-IDFG1
C
C *** Determine offset from diagonal 
      DO 112 IDIA=1,NDIA
      IF (IDFG .EQ. KDIA(IDIA)) THEN
         IOFF=KDIAS(IDIA)
         GOTO 113
      ENDIF
112   CONTINUE
113   IF (IDFG .GE. 0) THEN
       IGLOB=IOFF+IDFG1-1
      ELSE
       IGLOB=IOFF+IDFG2-1
      ENDIF
      KENTRY(I,J)=IGLOB
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices
      DO 120 IVE = 1, NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
120   CONTINUE
C
C *** Jacobian of the affine linear mapping onto the reference element
      DJAC(1,1)=DX(2)-DX(1)
      DJAC(1,2)=DX(3)-DX(1)
      DJAC(2,1)=DY(2)-DY(1)
      DJAC(2,2)=DY(3)-DY(1)
      DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
C
C *** Loop over all cubature points
      DO 200 ICUBP=1,NCUBP
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
      OM=DOMEGA(ICUBP)*DETJ
C *** ELE needs the information ICUBP because of preceeding
C *** dummy call using IPAR = -2
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
      IF (.NOT.BCON0) THEN
       XX=DX(1)+XI2*DJAC(1,1)+XI3*DJAC(1,2)
       YY=DY(1)+XI2*DJAC(2,1)+XI3*DJAC(2,2)
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
      DO 310 JDOFE=1,IDFL
      DB=DBAS(KDFL(JDOFE),IA)
      IF (BSYMM) JDOFE1=JDOFE
      DO 310 IDOFE=JDOFE1,IDFL
      JCOLB=KENTRY(JDOFE,IDOFE)+KOFF(IBLOC)
      DA(JCOLB)=DA(JCOLB)+DB*DBAS(KDFL(IDOFE),IB)*AUX
310   CONTINUE
301   CONTINUE
C
      BFIRST=.FALSE.
300   CONTINUE
C
200   CONTINUE
100   CONTINUE
C
      IF (ISYMM.GT.1) THEN
       DO 400 IBLOC=1,NBLOC
C
       DO 401 IDIA=2,(NDIA+1)/2
       IDIA1=NDIA+2-IDIA
       DO 402 IDG=1,KDIAS(IDIA+1)-KDIAS(IDIA)
       IDIAB=KDIAS(IDIA)-1+IDG
       IDIA1B=KDIAS(IDIA1)-1+IDG
C ***  Elements IDIA and IDIA1 in block IBLOC
       DA(IDIAB+KOFF(IBLOC))=DA(IDIA1B+KOFF(IBLOC))
402    CONTINUE
401    CONTINUE
C
400    CONTINUE
      ENDIF
C
99999 END
