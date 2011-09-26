************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT3D  (Release 1.1)               *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Turek                                 *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* XVB0                                                                 *
*                                                                      *
* Purpose  Allocation of the blocks in vector DB on DWORK              *
*          Call of VB0                                                 *
*                                                                      *
* Subroutines/functions called  VB0, ZNEW, ZDISP, ZLEN, ZTYPE, ZCLEAR  *
*                                                                      *
* Version from  10/01/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LB       I*4    Vector of numbers of the arrays                      *
*                 New arrays are allocated on DWORK for LB(IBLOC)=0    *
* NEQ      I*4    Length of each vector                                *
* ICLEAR   I*4    ICLEAR=1  all vectors are cleared                    *
* BSNGL    LOG                                                         *
* ARR      C*6    Names of blocks (for error messages only)            *
* For the description of the remaining parameters see VB0              *
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
* VB0                                                                  *
*                                                                      *
* Purpose  Calculation of NBLOC vectors corresponding to               *
*          linear forms  l(u)                                          *
*          Quadrilateral elements - trilinear transformation           *
*                                                                      *
* Subroutines/functions called NDFL, NDFGL, CB3H                       *
*                                                                      *
* Version from  10/01/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NBLOC    I*4    Number of vectors stored on DB                       *
* KOFF     I*4    Vector IBLOC starts at position KOFF(IBLOC)+1 on DB  *
* KVERT    I*4                                                         *
* KEDGE    I*4                                                         *
* KAREA    I*4                                                         *
* DCORVG   R*8    Arrays describing the triangulation                  *
* ELE      SUBR   EXTERNAL SUBROUTINE - values of basis functions      *
* COEFF    R*8    EXTERNAL FUNCTION - coefficients of the linear forms *
* BCON     LOG    BCON(IBLOC)=.TRUE. means constant coefficients in    *
*                 vector IBLOC                                         *
* KB       I*4    Multiindices occuring in the linear forms            *
*                 specified separately for each vector                 *
* KBN      I*4    Number of additive terms in each linear form         *
* ICUB     I*4    Number of cubature formula in CB3H                   *
* ILINT    I*4    0  full trilinear transformation                     *
*                 1  linear transformation only                        *
*                 2  axiparallel grid                                  *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DB       R*8    Calculated vectors                                   *
* IER      I*4    Error indicator                                      *
*                 -117 Wrong value in array KB                         *
*                                                                      *
************************************************************************
C
      SUBROUTINE XVB0(LB,NEQ,NBLOC,ICLEAR,ELE,COEFF,BCON,KB,KBN,
     *                ICUB,ILINT,BSNGL,ARR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6
      PARAMETER (NNARR=299,NNAB=21,NNDER=10)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION LB(*),KB(NNAB,*),KBN(*),BCON(*),ARR(*),BSNGL(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL COEFF,ELE
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAA/
C
      SUB='XVB0'
      IF (ICHECK.GE.997) CALL OTRC('XVB0  ','10/01/89')
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
        CALL WERR(-114,'XVB0  ')
        GOTO 99999
       ENDIF
       CALL ZLEN(LB(IBLOC),ILEN)
       IF (ILEN.LT.NEQ) THEN
        WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
        CALL WERR(-115,'XVB0  ')
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
C$DIR SCALAR
      DO 2 IBLOC=1,NBLOC
      KWORK(L(LOFF)+IBLOC-1)=L(LB(IBLOC))-1
2     CONTINUE
C
      CALL VB0(DWORK(1),NBLOC,KWORK(L(LOFF)),KWORK(L(LVERT)),
     *         KWORK(L(LEDGE)),KWORK(L(LAREA)),DWORK(L(LCORVG)),
     *         ELE,COEFF,BCON,DWORK(L(LOECON)),KB,KBN,ICUB,ILINT)
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
      SUBROUTINE VB0(DB,NBLOC,KOFF,KVERT,KEDGE,KAREA,DCORVG,
     *               ELE,COEFF,BCON,COECON,KB,KBN,ICUB,ILINT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNAB=21,NNDIM=3)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
      DIMENSION KVERT(NNVE,*),KEDGE(NNEE,*),KAREA(NNAE,*),KOFF(*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS),DB(*),DCORVG(3,*)
      DIMENSION BCON(*),KB(NNAB,*),KBN(*),COECON(NNDER,*)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,
     *                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),
     *                IEL,NDIM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/ELEM/,/CUB/,/COAUX1/
C
      SUB='VB0'
      IF (ICHECK.GE.997) CALL OTRC('VB0   ','10/01/89')
C
C *** Preparation - evaluation of parameters
      IER=0
C
C *** Which derivatives of basis functions are needed?
      DO 1 IDER=1,NNDER
1     BDER(IDER)=.FALSE.
      DO 2 IBLOC=1,NBLOC
      DO 3 IBN=1,KBN(IBLOC)
      IB=KB(IBN,IBLOC)
      IF (IB.LE.0.OR.IB.GT.NNDER) THEN
       WRITE (CPARAM,'(I15)') IBLOC
       CALL WERR(-117,'VB0   ')
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
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
      BFIRST=.TRUE.
C *** Dummy call of COEFF for nonlinear problems
C *** COEFF must set BDER(IDER)=.TRUE. if derivative IDER is needed
      AUX=COEFF(0D0,0D0,0D0,-1,0,BFIRST)
C
      BCON0=.TRUE.
      DO 4 IBLOC=1,NBLOC
      IF (BCON(IBLOC)) THEN
       DO 5 IBN=1,KBN(IBLOC)
       IB=KB(IBN,IBLOC)
5      COECON(IB,IBLOC)=COEFF(0D0,0D0,0D0,IB,IBLOC,BFIRST)
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
C *** Set zero elements of Jacobian for axiparallel grid
      IF (ILINT.EQ.2) THEN
       DJAC(1,3)=0D0
       DJAC(2,3)=0D0
       DJAC(3,1)=0D0
       DJAC(3,2)=0D0
      ENDIF
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Evaluation of coordinates of the vertices
      DO 110 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
      DZ(IVE)=DCORVG(3,JP)
110   CONTINUE
C
      IF (ILINT.EQ.2) THEN
       DJ11=(DX(2)+DX(4))*Q2
       DJ12=(DY(2)+DY(4))*Q2
       DJ13=(DZ(1)+DZ(5))*Q2
       DJAC(1,1)=(-DX(1)+DX(2))*Q2
       DJAC(2,1)=(-DY(1)+DY(2))*Q2
       DJAC(1,2)=(-DX(1)+DX(4))*Q2
       DJAC(2,2)=(-DY(1)+DY(4))*Q2
       DJAC(3,3)=(-DZ(1)+DZ(5))*Q2
       DETJ=DJAC(3,3)*(DJAC(1,1)*DJAC(2,2)-DJAC(2,1)*DJAC(1,2))
      ELSE IF (ILINT.EQ.1) THEN
       DJ11=(DX(1)+DX(2)+DX(3)+DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
       DJ12=(DY(1)+DY(2)+DY(3)+DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
       DJ13=(DZ(1)+DZ(2)+DZ(3)+DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
       DJAC(1,1)=(-DX(1)+DX(2)+DX(3)-DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
       DJAC(2,1)=(-DY(1)+DY(2)+DY(3)-DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
       DJAC(3,1)=(-DZ(1)+DZ(2)+DZ(3)-DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
       DJAC(1,2)=(-DX(1)-DX(2)+DX(3)+DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
       DJAC(2,2)=(-DY(1)-DY(2)+DY(3)+DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
       DJAC(3,2)=(-DZ(1)-DZ(2)+DZ(3)+DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
       DJAC(1,3)=(-DX(1)-DX(2)-DX(3)-DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
       DJAC(2,3)=(-DY(1)-DY(2)-DY(3)-DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
       DJAC(3,3)=(-DZ(1)-DZ(2)-DZ(3)-DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
       DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      ELSE
       DJ11=( DX(1)+DX(2)+DX(3)+DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
       DJ12=( DY(1)+DY(2)+DY(3)+DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
       DJ13=( DZ(1)+DZ(2)+DZ(3)+DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
       DJ21=(-DX(1)+DX(2)+DX(3)-DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
       DJ22=(-DY(1)+DY(2)+DY(3)-DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
       DJ23=(-DZ(1)+DZ(2)+DZ(3)-DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
       DJ31=(-DX(1)-DX(2)+DX(3)+DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
       DJ32=(-DY(1)-DY(2)+DY(3)+DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
       DJ33=(-DZ(1)-DZ(2)+DZ(3)+DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
       DJ41=(-DX(1)-DX(2)-DX(3)-DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
       DJ42=(-DY(1)-DY(2)-DY(3)-DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
       DJ43=(-DZ(1)-DZ(2)-DZ(3)-DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
       DJ51=( DX(1)-DX(2)+DX(3)-DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
       DJ52=( DY(1)-DY(2)+DY(3)-DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
       DJ53=( DZ(1)-DZ(2)+DZ(3)-DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
       DJ61=( DX(1)-DX(2)-DX(3)+DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
       DJ62=( DY(1)-DY(2)-DY(3)+DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
       DJ63=( DZ(1)-DZ(2)-DZ(3)+DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
       DJ71=( DX(1)+DX(2)-DX(3)-DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
       DJ72=( DY(1)+DY(2)-DY(3)-DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
       DJ73=( DZ(1)+DZ(2)-DZ(3)-DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
       DJ81=(-DX(1)+DX(2)-DX(3)+DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
       DJ82=(-DY(1)+DY(2)-DY(3)+DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
       DJ83=(-DZ(1)+DZ(2)-DZ(3)+DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
      ENDIF
C
C *** Loop over all cubature points
      DO 200 ICUBP=1,NCUBP
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C *** Jacobian of the trilinear mapping onto the reference element
      IF (ILINT.EQ.0) THEN
       DJAC(1,1)=DJ21+DJ51*XI2+DJ61*XI3+DJ81*XI2*XI3
       DJAC(1,2)=DJ31+DJ51*XI1+DJ71*XI3+DJ81*XI1*XI3
       DJAC(1,3)=DJ41+DJ61*XI1+DJ71*XI2+DJ81*XI1*XI2
       DJAC(2,1)=DJ22+DJ52*XI2+DJ62*XI3+DJ82*XI2*XI3
       DJAC(2,2)=DJ32+DJ52*XI1+DJ72*XI3+DJ82*XI1*XI3
       DJAC(2,3)=DJ42+DJ62*XI1+DJ72*XI2+DJ82*XI1*XI2
       DJAC(3,1)=DJ23+DJ53*XI2+DJ63*XI3+DJ83*XI2*XI3
       DJAC(3,2)=DJ33+DJ53*XI1+DJ73*XI3+DJ83*XI1*XI3
       DJAC(3,3)=DJ43+DJ63*XI1+DJ73*XI2+DJ83*XI1*XI2
       DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      ENDIF
C      DETJ=ABS(DETJ)
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C *** ELE needs the information ICUBP because of preceeding
C *** dummy call using IPAR = -2
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
      IF (.NOT.BCON0) THEN
       IF (ILINT.EQ.2) THEN
        XX=DJ11+DJAC(1,1)*XI1+DJAC(1,2)*XI2
        YY=DJ12+DJAC(2,1)*XI1+DJAC(2,2)*XI2
        ZZ=DJ13+DJAC(3,3)*XI3
       ELSE IF (ILINT.EQ.1) THEN
        XX=DJ11+DJAC(1,1)*XI1+DJAC(1,2)*XI2+DJAC(1,3)*XI3
        YY=DJ12+DJAC(2,1)*XI1+DJAC(2,2)*XI2+DJAC(2,3)*XI3
        ZZ=DJ13+DJAC(3,1)*XI1+DJAC(3,2)*XI2+DJAC(3,3)*XI3
       ELSE
        XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
        YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
        ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2
       ENDIF
      ENDIF
C
C *** Summing up over all multiindices
      BFIRST=.TRUE.
      DO 300 IBLOC=1,NBLOC
C
      DO 301 IBN=1,KBN(IBLOC)
      IB=KB(IBN,IBLOC)
      IF (.NOT.BCON(IBLOC)) THEN
       AUX=COEFF(XX,YY,ZZ,IB,IBLOC,BFIRST)*OM
      ELSE
       AUX=COECON(IB,IBLOC)*OM
      ENDIF
C
      DO 302 JDOFE=1,IDFL
      IGLOB=KDFG(JDOFE)+KOFF(IBLOC)
C$DIR SCALAR
      DO 303 IDIM=1,NDIM
303   DB(IGLOB)=DB(IGLOB)+DBAS(IDIM,KDFL(JDOFE),IB)*AUX
302   CONTINUE
301   CONTINUE
C
      BFIRST=.FALSE.
300   CONTINUE
C
200   CONTINUE
100   CONTINUE
C
99999 END
