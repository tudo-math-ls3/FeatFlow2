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
* XVA1                                                                 *
*                                                                      *
* Purpose  Allocation of the blocks in vector DB on DWORK              *
*          Call of VA1                                                 *
*                                                                      *
* Subroutines/functions called  SIVB, VA1, ZNEW, ZDISP, ZLEN,          *
*                               ZTYPE, ZCLEAR                          *
*                                                                      *
* Version from  12/04/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LB       I*4    Vector of numbers of the arrays                      *
*                 New arrays are allocated on DWORK for LB(IBLOC)=0    *
* NEQ      I*4    Length of each vector                                *
* ICLEAR   I*4    ICLEAR=1  all vectors are cleared                    *
* TMAX     R*8    EXTERNAL FUNCTION - parametrization of the domain    *
* DPAR1    R*8    Description of boundary part GAMMA1                  *
* DPAR2    R*8                                                         *
* ARR      C*6    Names of blocks (for error messages only)            *
* BSNGL    L*4    NBLOC logical values                                 *
*                 .TRUE. means that corresponding array is converted   *
*                 to single precision after completion                 *
* For the description of the remaining parameters see VA1              *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* LB       I*4    Vector of numbers of arrays on DB                    *
* IER      I*4    Error indicator                                      *
*                 -114 Type of at least one array is not double prec.  *
*                 -115 Length of at least one array is < NEQ           *
*                                                                      *
************************************************************************
*                                                                      *
* VA1                                                                  *
*                                                                      *
* Purpose  Calculation of NBLOC vectors corresponding to               *
*          linear forms  l(u) (boundary integrals over gamma1)         *
*          Triangular elements - affine linear transformation          *
*          The vectors are stored subsequently on the vector DB        *
*          Double precision version                                    *
*                                                                      *
* Subroutines/functions called NDFL, NDFGL, CB1                        *
*                                                                      *
* Version from  12/04/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NBLOC    I*4    Number of vectors stored on DB                       *
* KOFF     I*4    Vector IBLOC starts at position KOFF(IBLOC)+1 on DB  *
* KVERT    I*4                                                         *
* KMID     I*4                                                         *
* DCORVG   R*8    Arrays describing the triangulation                  *
* KBCT     I*4                                                         *
* IBCT     I*4                                                         *
* KVBD     I*4    Vector of vertices on the boundary                   *
* KEBD     I*4    Vector of corresponding element numbers              *
* IVBD1    I*4    Minimum and maximum index on KVBD -                  *
* IVBD2    I*4    calculated by SIVB                                   *
* ELE      SUBR   EXTERNAL SUBROUTINE - values of basis functions      *
* COEFF    R*8    EXTERNAL FUNCTION - coefficients of the linear forms *
* BCON     LOG    BCON(IBLOC)=.TRUE. means constant coefficients in    *
*                 vector IBLOC                                         *
* KB       I*4    Multiindices occuring in the linear forms            *
*                 specified separately for each vector                 *
* KBN      I*4    Number of additive terms in each linear form         *
* ICUB     I*4    Number of cubature formula in CB1                    *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DB       R*8    Calculated vectors                                   *
* IER      I*4    Error indicator                                      *
*                 -117 Wrong value in array KB                         *
*                                                                      *
************************************************************************
C
      SUBROUTINE XVA1(LB,NEQ,NBLOC,ICLEAR,TMAX,DPAR1,DPAR2,
     *                IBCT,ELE,COEFF,BCON,KB,KBN,ICUB,ARR,BSNGL)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6
C
      PARAMETER (NNARR=299,NNAB=21,NNDER=6)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION LB(*),KB(NNAB,*),KBN(*),BCON(*),ARR(*),BSNGL(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL TMAX,COEFF,ELE
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAA/
C
      SUB='XVA1'
      IF (ICHECK.GE.997) CALL OTRC('XVA1  ','12/04/90')
      IER=0
C
      DO 1 IBLOC=1,NBLOC
      IF (LB(IBLOC).EQ.0) THEN
       CALL ZNEW(NEQ,1,LB(IBLOC),ARR(IBLOC))
       IF (IER.NE.0) GOTO 99999
      ELSE
C ***  Check input parameter
       CALL ZTYPE(LB(IBLOC),ITYPE)
       IF (ITYPE.NE.1) THEN
        WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
        CALL WERR(-114,'XVA1  ')
        GOTO 99999
       ENDIF
       CALL ZLEN(LB(IBLOC),ILEN)
       IF (ILEN.LT.NEQ) THEN
        WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
        CALL WERR(-115,'XVA1  ')
        GOTO 99999
       ENDIF
       IF (ICLEAR.EQ.1) CALL ZCLEAR(LB(IBLOC),ARR(IBLOC))
      ENDIF
1     CONTINUE
C
      CALL ZNEW(NBLOC,-3,LOFF,'KOFF  ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NBLOC*NNDER,1,LOECON,'COECON')
      IF (IER.NE.0) GOTO 99999
C
      DO 2 IBLOC=1,NBLOC
      KWORK(L(LOFF)+IBLOC-1)=L(LB(IBLOC))-1
2     CONTINUE
C
C *** Determine first and last index on gamma1
      CALL SIVB(KWORK(L(LVBD)),IBCT,KWORK(L(LBCT)),TMAX,
     *          DWORK(L(LCORVG)),DPAR1,DPAR2,IVBD1,IVBD2)
      IF (IER.NE.0) GOTO 99999
C
      CALL VA1(DWORK(1),NBLOC,KWORK(L(LOFF)),KWORK(L(LVERT)),
     *         KWORK(L(LMID)),DWORK(L(LCORVG)),
     *         KWORK(L(LBCT)),IBCT,KWORK(L(LVBD)),KWORK(L(LEBD)),
     *         IVBD1,IVBD2,ELE,COEFF,BCON,DWORK(L(LOECON)),
     *         KB,KBN,ICUB)
      IF (IER.NE.0) GOTO 99999
C
      CALL ZDISP(0,LOECON,'COECON')
      CALL ZDISP(0,LOFF,'KOFF  ')
C
      DO 3 IBLOC=NBLOC,1,-1
      IF (BSNGL(IBLOC)) THEN
       CALL ZCTYPE(2,LB(IBLOC),ARR(IBLOC))
       IF (IER.NE.0) GOTO 99999
      ENDIF
3     CONTINUE
C
99999 END
C
C
C
      SUBROUTINE VA1(DB,NBLOC,KOFF,KVERT,KMID,DCORVG,
     *               KBCT,IBCT,KVBD,KEBD,IVBD1,IVBD2,
     *               ELE,COEFF,BCON,COECON,KB,KBN,ICUB)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KOFF(*)
      DIMENSION KVBD(*),KEBD(*),KBCT(*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS),DB(*),DCORVG(2,*)
      DIMENSION BCON(*),KB(NNAB,*),KBN(*),COECON(NNDER,*)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/ELEM/,/CUB/,/COAUX1/
	DATA NEBD/0/
C
      SUB='VA1'
      IF (ICHECK.GE.997) CALL OTRC('VA1   ','12/04/90')
C
C *** Preparation - evaluation of parameters
      IER=0
      Q2=1.0D0/SQRT(2.0D0)
C
C *** Which derivatives of basis functions are needed?
      DO 1 IDER=1,NNDER
1     BDER(IDER)=.FALSE.
      DO 2 IBLOC=1,NBLOC
      DO 3 IBN=1,KBN(IBLOC)
      IB=KB(IBN,IBLOC)
      IF (IB.LE.0.OR.IB.GT.NNDER) THEN
       WRITE (CPARAM,'(I15)') IBLOC
       CALL WERR(-117,'VA1   ')
       GOTO 99999
      ENDIF
3     BDER(IB)=.TRUE.
2     CONTINUE
C
C *** Dummy call of ELE sets number of element
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IF (IER.NE.0) GOTO 99999
      IDFL=NDFL(IELTYP)
      IF (IER.LT.0) GOTO 99999
      CALL CB1(ICUB)
      IF (IER.NE.0) GOTO 99999
      BFIRST=.TRUE.
C *** Dummy call of COEFF for nonlinear problems
C *** COEFF must set BDER(IDER)=.TRUE. if derivative IDER is needed
      AUX=COEFF(0D0,0D0,-1,0,BFIRST)
C
      BCON0=.TRUE.
      DO 4 IBLOC=1,NBLOC
      IF (BCON(IBLOC)) THEN
       DO 5 IBN=1,KBN(IBLOC)
       IB=KB(IBN,IBLOC)
5      COECON(IB,IBLOC)=COEFF(0D0,0D0,IB,IBLOC,BFIRST)
      ELSE
       BCON0=.FALSE.
      ENDIF
4     CONTINUE
************************************************************************
C *** Calculation of the linear form
************************************************************************
C *** Dummy call - ELE may save arithmetic operations
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
C
C *** Loop over all nodes on the boundary gamma1
      DO 100 IVBD=IVBD1,IVBD2
      IEL=KEBD(IVBD)
      JVBD1=KVBD(IVBD)
      JVBD2=KVBD(IVBD+1)
      IF (IVBD2.EQ.(KBCT(IBCT+1)-1)) JVBD2=KBCT(IBCT)
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Evaluation of coordinates of the vertices
      DO 110 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
C *** Determine local number NEBD of the edge
      IF (JP.EQ.JVBD1) NEBD=IVE
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
110   CONTINUE
C
      DS=SQRT((DX(MOD(NEBD,NVE)+1)-DX(NEBD))**2+
     *        (DY(MOD(NEBD,NVE)+1)-DY(NEBD))**2)
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
      GOTO (210,220,230),NEBD
210   XI1= 0.5D0*(DXI(ICUBP,1)+1D0)
      XI2=-0.5D0*(DXI(ICUBP,1)-1D0)
      XI3= 0D0
      OM=0.5D0*DS*DOMEGA(ICUBP)
      GOTO 240
220   XI1= 0D0
      XI2= 0.5D0*(DXI(ICUBP,1)+1D0)
      XI3=-0.5D0*(DXI(ICUBP,1)-1D0)
      OM=0.5D0*Q2*DS*DOMEGA(ICUBP)
      GOTO 240
230   XI1=-0.5D0*(DXI(ICUBP,1)-1D0)
      XI2= 0D0
      XI3= 0.5D0*(DXI(ICUBP,1)+1D0)
      OM=0.5D0*DS*DOMEGA(ICUBP)
240   CONTINUE
C
C *** ELE needs the information ICUBP because of preceeding
C *** dummy call using IPAR = -2
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
      IF (.NOT.BCON0) THEN
       XX=DX(1)+XI2*DJAC(1,1)+XI3*DJAC(1,2)
       YY=DY(1)+XI2*DJAC(2,1)+XI3*DJAC(2,2)
      ENDIF
C
C *** Summing up over all multiindices
      BFIRST=.TRUE.
      DO 300 IBLOC=1,NBLOC
      DO 301 IBN=1,KBN(IBLOC)
      IB=KB(IBN,IBLOC)
      IF (.NOT.BCON(IBLOC)) THEN
       AUX=COEFF(XX,YY,IB,IBLOC,BFIRST)*OM
      ELSE
       AUX=COECON(IB,IBLOC)*OM
      ENDIF
C
      DO 310 JDOFE=1,IDFL
      IGLOB=KDFG(JDOFE)+KOFF(IBLOC)
      DB(IGLOB)=DB(IGLOB)+DBAS(KDFL(JDOFE),IB)*AUX
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
