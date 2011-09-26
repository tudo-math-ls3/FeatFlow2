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
* XAB07                                                                *
*                                                                      *
* Purpose  Allocation of the blocks in vector DA , KOFF and COECON     *
*          Call of AB07                                                *
*                                                                      *
* Subroutines/functions called  AB07, ZNEW, ZDISP, ZLEN, ZTYPE, ZCLEAR *
*                                                                      *
* Version from  12/18/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LA       I*4    Vector of numbers of the arrays                      *
*                 New arrays are allocated on DWORK for LA(IBLOC)=0    *
* LCOL     I*4    Numbers of the pointer arrays KCOL and KLD           *
* LLD      I*4    calculated by AP7                                    *
* ICLEAR   I*4    ICLEAR=1  all matrices are cleared                   *
* BSNGL    LOG    =.TRUE.  conversion to single precision              *
* ARR      C*6    Names of blocks (for error messages only)            *
* For the description of the remaining parameters see AB07             *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* LA       I*4    Vector of numbers of arrays on DA(VA)                *
* IER      I*4    Error indicator                                      *
*                 -114 Type of at least one array is not double prec.  *
*                 -115 Length of at least one array is < NA            *
*                                                                      *
************************************************************************
*                                                                      *
* AB07                                                                 *
*                                                                      *
* Purpose  Calculation of NBLOC matrices corresponding to              *
*          bilinear forms  a(v,u)                                      *
*          Quadrilateral elements - trilinear transformation           *
*          Same trial and test functions                               *
*          Storage technique 7/8                                       *
*          Each matrix makes use of the same pointer vectors           *
*                                                                      *
* Subroutines/functions called NDFL, NDFGL, CB3H                       *
*                                                                      *
* Version from  12/18/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* KCOL     I*4    Pointer vectors for each matrix on DA corresponding  *
* KLD      I*4    to storage technique 7/8 (calculated by AP7)         *
* NA       I*4    Number of elements per matrix                        *
* NEQ      I*4    Number of equations per matrix                       *
* NBLOC    I*4    Number of matrices stored on DA                      *
* KOFF     I*4    Matrix IBLOC starts at position KOFF(IBLOC)+1 on DA  *
* KVERT    I*4    Arrays describing the triangulation                  *
* KEDGE    I*4                                                         *
* KAREA    I*4                                                         *
* DCORVG   R*8                                                         *
* ELE      SUBR   EXTERNAL SUBROUTINE - values of basis functions      *
* COEFF    R*8    EXTERNAL FUNCTION - coefficients of the bilinear form*
* BCON     LOG    BCON(IBLOC)=.TRUE. means constant coefficients in    *
*                 matrix IBLOC                                         *
* KAB      I*4    Pairs of multiindices occuring in the bilinear forms *
*                 specified separately for each matrix                 *
* KABN     I*4    Number of additive terms in each bilinear form       *
* ICUB     I*4    Number of cubature formula in CB3H                   *
* ISYMM    I*4    0  storage technique 7                               *
*                 1  storage technique 8                               *
*                 2  storage technique 7 but for symmetric matrices    *
* ILINT    I*4    0  full trilinear transformation                     *
*                 1  linear transformation only                        *
*                 2  axiparallel grid                                  *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DA       R*8    Calculated matrices                                  *
* IER      I*4    Error indicator                                      *
*                 -116 Wrong value in array KAB                        *
*                                                                      *
************************************************************************
C
      SUBROUTINE XAB07(LA,LCOL,LLD,NA,NEQ,NBLOC,ICLEAR,ELE,
     *                 COEFF,BCON,KAB,KABN,ICUB,ISYMM,ILINT,BSNGL,ARR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6
      PARAMETER (NNARR=299,NNAB=21,NNDER=10)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION LA(*),KAB(2,NNAB,*),KABN(*),BCON(*),ARR(*),BSNGL(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL COEFF,ELE
      SAVE /TRIAA/,/OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='XAB07'
      IF (ICHECK.GE.997) CALL OTRC('XAB07 ','12/18/89')
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
        CALL WERR(-114,'XAB07 ')
        GOTO 99999
       ENDIF
       CALL ZLEN(LA(IBLOC),ILEN)
       IF (ILEN.LT.NA) THEN
        WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
        CALL WERR(-115,'XAB07 ')
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
C$DIR SCALAR
      DO 3 IBLOC=1,NBLOC
      KWORK(L(LOFF)+IBLOC-1)=L(LA(IBLOC))-1
3     CONTINUE
C
      CALL AB07(DWORK(1),KWORK(L(LCOL)),KWORK(L(LLD)),NA,NEQ,
     *          NBLOC,KWORK(L(LOFF)),KWORK(L(LVERT)),KWORK(L(LEDGE)),
     *          KWORK(L(LAREA)),DWORK(L(LCORVG)),ELE,COEFF,BCON,
     *          DWORK(L(LOECON)),KAB,KABN,ICUB,ISYMM,ILINT)
      IF (IER.NE.0) GOTO 99999
C
      CALL ZDISP(0,LOECON,'COECON')
      CALL ZDISP(0,LOFF,'KOFF  ')
C
      DO 4 IBLOC=NBLOC,1,-1
      IF (BSNGL(IBLOC)) THEN
       CALL ZCTYPE(2,LA(IBLOC),ARR(IBLOC))
       IF (IER.NE.0) GOTO 99999
      ENDIF
4     CONTINUE
C
99999 END
C
C
C
      SUBROUTINE AB07(DA,KCOL,KLD,NA,NEQ,NBLOC,KOFF,KVERT,KEDGE,KAREA,
     *                DCORVG,ELE,COEFF,BCON,COECON,KAB,KABN,ICUB,ISYMM,
     *                ILINT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNAB=21,NNDIM=3)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
      DIMENSION KCOL(*),KVERT(NNVE,*),KEDGE(NNEE,*),KAREA(NNAE,*)
      DIMENSION KOFF(*),KLD(*),KDFG(NNBAS),KDFL(NNBAS),DA(*),DCORVG(3,*)
      DIMENSION BCON(*),KAB(2,NNAB,*),KABN(*),COECON(NNDER,NNDER,*)
      DIMENSION DB(NNDIM),KENTRY(NNBAS,NNBAS)
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
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/ELEM/,/TRIAD/,/CUB/,/COAUX1/
C
      SUB='AB07'
      IF (ICHECK.GE.997) CALL OTRC('AB07  ','12/18/89')
C
C *** Preparation - evaluation of parameters
      IER=0
      NDIM=1
      BSYMM=ISYMM.GE.1
C *** Which derivatives of basis functions are needed?
      DO 1 I=1,NNDER
1     BDER(I)=.FALSE.
      DO 2 IBLOC=1,NBLOC
      DO 3 I=1,KABN(IBLOC)
      DO 3 J=1,2
      I1=KAB(J,I,IBLOC)
      IF (I1.LE.0.OR.I1.GT.NNDER) THEN
       WRITE (CPARAM,'(I15)') IBLOC
       CALL WERR(-116,'AB07  ')
       GOTO 99999
      ENDIF
3     BDER(I1)=.TRUE.
2     CONTINUE
C *** Dummy call of ELE sets number of element
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
C
      BFIRST=.TRUE.
C *** Dummy call of COEFF for nonlinear problems
C *** COEFF must set BDER(IDER)=.TRUE. if derivative IDER is needed
      AUX=COEFF(0D0,0D0,0D0,-1,-1,0,BFIRST)
      BCON0=.TRUE.
C
      DO 4 IBLOC=1,NBLOC
      IF (BCON(IBLOC)) THEN
       DO 5 I=1,KABN(IBLOC)
       IA=KAB(1,I,IBLOC)
       IB=KAB(2,I,IBLOC)
5      COECON(IA,IB,IBLOC)=COEFF(0D0,0D0,0D0,IA,IB,IBLOC,BFIRST)
      ELSE
       BCON0=.FALSE.
      ENDIF
4     CONTINUE
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
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
      JDOFE1=1
      DO 100 IEL=1,NEL
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLD(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
      IF (BSYMM) JDOFE1=JDOFE
      JCOL0=ILD
      DO 111 IDOFE=JDOFE1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOL(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices
      DO 120 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
      DZ(IVE)=DCORVG(3,JP)
120   CONTINUE
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
C       IF (ICHECK.GT.0) THEN
C       WRITE(*,'(A,I3,A,8I4)') 'Element ',IEL,' Nodes ',(KVE(I),I=1,8)
C       WRITE(*,'(A,8D10.3)') 'X',(DX(I),I=1,8)
C       WRITE(*,'(A,8D10.3)') 'Y',(DY(I),I=1,8)
C       WRITE(*,'(A,8D10.3)') 'Z',(DZ(I),I=1,8)
C       WRITE(*,'(A,3D12.4)') 'D0', DJ11,DJ12,DJ13
C       WRITE(*,'(A,3D12.4)') 'D1',(DJAC(1,I),I=1,3)
C       WRITE(*,'(A,3D12.4)') 'D2',(DJAC(2,I),I=1,3)
C       WRITE(*,'(A,3D12.4)') 'D3',(DJAC(3,I),I=1,3)
C       ENDIF
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
C
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
C *** Summing up over all pairs of multiindices
      BFIRST=.TRUE.
      DO 300 IBLOC=1,NBLOC
C
      DO 301 IALBET=1,KABN(IBLOC)
      IA=KAB(1,IALBET,IBLOC)
      IB=KAB(2,IALBET,IBLOC)
      IF (.NOT.BCON(IBLOC)) THEN
       AUX=COEFF(XX,YY,ZZ,IA,IB,IBLOC,BFIRST)*OM
      ELSE
       AUX=COECON(IA,IB,IBLOC)*OM
      ENDIF
C
      DO 302 JDOFE=1,IDFL
C$DIR SCALAR
      DO 303 IDIM=1,NDIM
303   DB(IDIM)=DBAS(IDIM,KDFL(JDOFE),IA)
      IF (BSYMM) JDOFE1=JDOFE
      DO 302 IDOFE=JDOFE1,IDFL
      JCOLB=KENTRY(JDOFE,IDOFE)+KOFF(IBLOC)
C$DIR SCALAR
      DO 304 IDIM=1,NDIM
304   DA(JCOLB)=DA(JCOLB)+DB(IDIM)*DBAS(IDIM,KDFL(IDOFE),IB)*AUX
302   CONTINUE
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
       DO 401 IROW=1,NEQ
       DO 402 ICOL=KLD(IROW)+1,KLD(IROW+1)-1
       J1=KCOL(ICOL)
       IF (J1.GE.IROW) GOTO 401
       DO 403 ICOL1=KLD(J1+1)-1,KLD(J1),-1
       IF (KCOL(ICOL1).EQ.IROW) GOTO 404
403    CONTINUE
       GOTO 402
404    ICOLB=ICOL+KOFF(IBLOC)
       ICOL1B=ICOL1+KOFF(IBLOC)
C ***  Elements ICOL and ICOL1 in block IBLOC
       DA(ICOLB)=DA(ICOL1B)
402    CONTINUE
401    CONTINUE
C
400    CONTINUE
      ENDIF
C
99999 END
