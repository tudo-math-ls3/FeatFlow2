************************************************************************
* WARNING:
*   The discretisation routines in this file have been replaced
*   by the cleaned-up routines in AP79X.F and AB79X.F.
*   The routines here are not used anymore, only kept for reference.
************************************************************************

************************************************************************
* XAP7X                                                                *
*                                                                      *
* Purpose  Call AP7                                                    *
*          Allocate KLD and KCOL on DWORK                              *
*                                                                      *
* Subroutines/functions called  AP7, EA00, ZNEW, ZDISP, ZFREE          *
*                                                                      *
* Version from  12/11/89                                               *
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

      SUBROUTINE XAP7X(LCOL,LLD,NA,NEQ,TRIA,ELE,ISYMM)

      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'stria.inc'
      
C     parameters

      INTEGER LCOL,LLD,NA,NEQ,ISYMM,TRIA(SZTRIA)
      
      EXTERNAL ELE
      
C     local variables
      
      INTEGER IELTYP,IWMAX0,IFREE,LCOL1,LIND
      INTEGER KVERT,KMID
      
      INTEGER NDFGX
      EXTERNAL NDFGX

      SUB='XAP7'
      IF (ICHECK.GE.997) CALL OTRC('XAP7  ','12/11/89')

C *** Determine total number of degrees of freedom
      CALL EA00(ELE,ELE,TRIA(ONVE),IELTYP)
      NEQ=NDFGX(IELTYP,TRIA)
      IF (IER.NE.0) GOTO 99999

C *** Allocate KLD on DWORK ***
      CALL ZNEW(NEQ+1,-3,LLD,'KLD   ')
      IF (IER.NE.0) GOTO 99999

C *** Determine free space on DWORK 
      IWMAX0=IWMAX
      CALL ZFREE(3,IFREE)
      IF (IER.NE.0) GOTO 99999

      NA=IFREE/3-3
      CALL ZNEW(NA,-3,LCOL,'KCOL  ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NA,-3,LCOL1,'KCOL1 ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NA,-3,LIND,'KIND  ')
      IF (IER.NE.0) GOTO 99999

      KVERT  = L(TRIA(OLVERT))
      KMID   = L(TRIA(OLMID))

      CALL AP7X(KWORK(L(LCOL)),KWORK(L(LCOL1)),
     *         KWORK(L(LIND)),KWORK(L(LLD)),
     *         NA,NEQ,TRIA,ELE,ISYMM, KWORK(KVERT),KWORK(KMID))
      IF (IER.NE.0) GOTO 99999

C *** Release space on DWORK not used for KCOL ***
      CALL ZDISP(NA,LIND,'KIND  ')
      CALL ZDISP(NA,LCOL1,'KCOL1 ')
      CALL ZDISP(NA,LCOL,'KCOL  ')

      IWMAX=MAX(IWORK,IWMAX0)
C
      CALL ZDISP(0,LIND,'KIND  ')
      CALL ZDISP(0,LCOL1,'KCOL1 ')
99999 END
      
      
************************************************************************
*                                                                      *
* AP7                                                                  *
*                                                                      *
* Purpose  Calculation of the pointer vectors KLD and KCOL             *
*          for a matrix corresponding to a given element type          *
*          Storage technique 7/8                                       *
*                                                                      *
* Subroutines/functions called  NDFL, NDFGL, EA00                      *
*                                                                      *
* Version from  12/11/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* KCOL1    I*4    Auxiliary array                                      *
* KIND     I*4    Auxiliary array                                      *
* NA       I*4    Maximal length of KCOL                               *
* NEQ      I*4    Number of equations                                  *
* ELE      SUBR   EXTERNAL Subroutine - evaluation of basis functions- *
*                 for dummy call only                                  *
* ISYMM    I*4    >=1 matrix symmetric - storage technique 8           *
* KVERT    I*4    Arrays describing the triangulation                  *
* KMID     I*4                                                         *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KCOL     I*4    Pointer vector containing column indices             *
* KLD      I*4    Pointer vector to diagonal elements                  *
* NA       I*4    Effective number of elements in KCOL                 *
* IER      I*4    Error indicator                                      *
*                 -118  Not enough space for KCOL                      *
*                       Error occured on element IEL                   *
*                                                                      *
************************************************************************
      
      SUBROUTINE AP7X(KCOL,KCOL1,KIND,KLD,NA,NEQ,TRIA,ELE,ISYMM,
     *                KVERT,KMID)
      
      IMPLICIT NONE

      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasicelem.inc'
      
      INCLUDE 'stria.inc'
      
C     parameters 

      INTEGER KCOL(*),KCOL1(*),KIND(*),KLD(*),NA,NEQ,ISYMM
      INTEGER KVERT(*),KMID(*)
      
      EXTERNAL ELE
      
C     local variables

      INTEGER IFREE,IEQ,IEL,IDOFE,JDOFE,IELTYP,IDOFE1,IROW
      INTEGER IPOS,ICOL,JCOL,IHELP,TRIA(SZTRIA)
      LOGICAL BSYMM,BSORT
      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL
      
      INTEGER NDFL
      EXTERNAL NDFL

      SUB='AP7'
      IF (ICHECK.GE.997) CALL OTRC('AP7   ','12/11/89')

      IER=0

      IF (NA.LT.NEQ) THEN
C *** DWORK exhausted
       WRITE (CPARAM,'(I15)') 1
       CALL WERR(-118,'AP7   ')
       GOTO 99999
      ENDIF
      IFREE=NA-NEQ
      NA=NEQ
C *** Initialization of KIND and KCOL1
      DO 10 IEQ=1,NEQ
      KIND(IEQ)=0
10    KCOL1(IEQ)=IEQ
C
C *** Set element number by dummy call ***
      CALL EA00(ELE,ELE,TRIA(ONVE),IELTYP)
C *** Determine number of degrees of freedom per element ***
      IDFL=NDFL(IELTYP)
      IF (IER.NE.0) GOTO 99999
C
      BSYMM=ISYMM.EQ.1
      IDOFE1=1
C *** Loop over elements
      DO 100 IEL=1,TRIA(ONEL)
C
C *** KDFG returns the global degrees of freedom in increasing order
      CALL NDFGLX(TRIA,IEL,0,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.NE.0) GOTO 99999
C
C *** Loop over local number of degrees of freedom
      DO 110 IDOFE=1,IDFL
      IROW=KDFG(IDOFE)
      IF (BSYMM) THEN
       IF (IROW.EQ.NEQ) GOTO 110
       IDOFE1=IDOFE
      ENDIF
C
C *** Loop over off-diagonal elements
C *** only upper triangular part in symmetric case
      DO 120 JDOFE=IDOFE1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 120
      JCOL=KDFG(JDOFE)
C *** JCOL is the global number of d.o.f. to be inserted into row IROW
C
C *** Look whether entry (IROW,JCOL) is already provided
      IPOS=IROW
121   IF (KIND(IPOS).EQ.0) THEN
C *** final entry in row IROW up to now
C
      IF (IFREE.LE.0) THEN
C *** DWORK exhausted
        WRITE (CPARAM,'(I15)') IEL
        CALL WERR(-118,'AP7   ')
        GOTO 99999
       ENDIF
C
C *** insert new entry
       NA=NA+1
       IFREE=IFREE-1
       KCOL1(NA)=JCOL
       KIND(IPOS)=NA
       KIND(NA)=0
       GOTO 120
C
      ELSE
C
       IPOS=KIND(IPOS)
       IF (KCOL1(IPOS).NE.JCOL) GOTO 121
C
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
      KCOL(NA)=IEQ
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
C *** Sort off-diagonal entries on KCOL separately for each row
C
      DO 300 IEQ=1,NEQ
C
301   BSORT=.TRUE.
      DO 302 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-2
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

************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek, M. Koester         *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
* XABM7X                                                               *
*                                                                      *
* Purpose  Allocation of the blocks in vector DA , KOFF and COECON     *
*          Call of ABM7                                                *
*                                                                      *
* Version from  10/16/05                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LA       I*4    Vector of numbers of the arrays                      *
*                 New arrays are allocated on DWORK for LA(IBLOC)=0    *
* LCOL     I*4    Numbers of the pointer arrays KCOL and KLD           *
* LLD      I*4    calculated by AP7                                    *
* ICLEAR   I*4    ICLEAR=1  all matrices are cleared                   *
* ARR      C*6    Names of blocks (for error messages only)            *
* BNONPR   L*4    =1 if the element is nonparametric                   *
* TRIA     S      structure of the triangulation                       *
* For the description of the remaining parameters see AB07             *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* LA       I*4    Vector of numbers of arrays on DA                    *
* IER      I*4    Error indicator                                      *
*                 -114 Type of at least one array is not double prec.  *
*                 -115 Length of at least one array is < NA            *
*                                                                      *
************************************************************************

      SUBROUTINE XAB7X(LA,LCOL,LLD,NA,NEQ,NBLOC,ICLEAR,TRIA,ELE,BNONPR,
     *                 COEFF,BCON,KAB,KABN,ICUB,ISYMM,ARR,IPARAM,DPARAM)

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasicelem.inc'
      
      INCLUDE 'stria.inc'

C parameters

      INTEGER LA(*),KAB(2,NNAB,*),KABN(*)
      LOGICAL BCON(*)

      CHARACTER ARR*12
      DIMENSION ARR(*)
      
      INTEGER LCOL, LLD, NA, NEQ, NBLOC, ICLEAR, ICUB, ISYMM
      INTEGER TRIA(SZTRIA)

      INTEGER IPARAM(*)
      DOUBLE PRECISION DPARAM(*)

      DOUBLE PRECISION COEFF
      EXTERNAL COEFF,ELE
      
      LOGICAL BNONPR

C local variables
      
      INTEGER IBLOC, ITYPE, ILEN, LOFF, LOECON
      INTEGER KVERT,KMID,KCORVG

      SUB='XABM7'
      IF (ICHECK.GE.997) CALL OTRC('XABM7 ','12/08/94')
      IER=0

      DO IBLOC=1,NBLOC
        IF (LA(IBLOC).EQ.0) THEN
          CALL ZNEW(NA,1,LA(IBLOC),ARR(IBLOC))
          IF (IER.NE.0) GOTO 99999
        ELSE
C ***  Check input parameters
          CALL ZTYPE(LA(IBLOC),ITYPE)
          IF (ITYPE.NE.1) THEN
            WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
            CALL WERR(-114,'XABM7 ')
            GOTO 99999
          ENDIF
          CALL ZLEN(LA(IBLOC),ILEN)
          IF (ILEN.LT.NA) THEN
            WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
            CALL WERR(-115,'XABM7 ')
            GOTO 99999
          ENDIF
          IF (ICLEAR.EQ.1) CALL ZCLEAR(LA(IBLOC),ARR(IBLOC))
        ENDIF
      END DO

      CALL ZNEW(NBLOC,-3,LOFF,'KOFF  ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NBLOC*NNDER**2,1,LOECON,'COECON')
      IF (IER.NE.0) GOTO 99999

      DO IBLOC=1,NBLOC
        KWORK(L(LOFF)+IBLOC-1)=L(LA(IBLOC))-1
      END DO
      
      KVERT  = L(TRIA(OLVERT))
      KMID   = L(TRIA(OLMID))
      KCORVG = L(TRIA(OLCORVG))
      
      IF (.NOT.BNONPR) THEN
      
        CALL AB7X(DWORK(1),KWORK(L(LCOL)),KWORK(L(LLD)),NA,NEQ,
     *          NBLOC,KWORK(L(LOFF)),KWORK(KVERT),KWORK(KMID),
     *          DWORK(KCORVG),TRIA,ELE,COEFF,
     *          BCON,DWORK(L(LOECON)),KAB,KABN,ICUB,ISYMM,
     *          IPARAM,DPARAM)

      ELSE 

        CALL AB7NX(DWORK(1),KWORK(L(LCOL)),KWORK(L(LLD)),NA,NEQ,
     *          NBLOC,KWORK(L(LOFF)),KWORK(KVERT),KWORK(KMID),
     *          DWORK(KCORVG),TRIA,ELE,COEFF,
     *          BCON,DWORK(L(LOECON)),KAB,KABN,ICUB,ISYMM,
     *          IPARAM,DPARAM)
     
      END IF

      IF (IER.NE.0) GOTO 99999

      CALL ZDISP(0,LOECON,'COECON')
      CALL ZDISP(0,LOFF,'KOFF  ')

99999 END
      
************************************************************************
* ABM7X                                                                *
*                                                                      *
* Purpose  Calculation of NBLOC matrices corresponding to              *
*          bilinear forms  a(v,u)                                      *
*          Nonparametric quadrilateral elements                        *
*          Same trial and test functions                               *
*          Storage technique 7/8                                       *
*          Each matrix makes use of the same pointer vectors           *
*          Double precision version                                    *
*                                                                      *
* Subroutines/functions called NDFL, NDFGL, CB2Q                       *
*                                                                      *
* Version from  10/16/05                                               *
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
* KMID     I*4                                                         *
* DCORVG   R*8                                                         *
* ELE      SUBR   EXTERNAL SUBROUTINE - values of basis functions      *
* COEFF    R*8    EXTERNAL FUNCTION - coefficients of the bilinear form*
* BCON     LOG    BCON(IBLOC)=.TRUE. means constant coefficients in    *
*                 matrix IBLOC                                         *
* KAB      I*4    Pairs of multiindices occuring in the bilinear forms *
*                 specified separately for each matrix                 *
* KABN     I*4    Number of additive terms in each bilinear form       *
* ICUB     I*4    Number of cubature formula in CB2Q                   *
* ISYMM    I*4    0  storage technique 7                               *
*                 1  storage technique 8                               *
*                 2  storage technique 7 but for symmetric matrices    *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DA       R*8    Calculated matrices                                  *
* IER      I*4    Error indicator                                      *
*                 -116 Wrong value in array KAB                        *
*                                                                      *
************************************************************************

C     Parametric version

      SUBROUTINE AB7X(DA,KCOL,KLD,NA,NEQ,NBLOC,KOFF,KVERT,KMID,
     *                DCORVG,TRIA,ELE,COEFF,BCON,COECON,
     *                KAB,KABN,ICUB,ISYMM,IPARAM,DPARAM)

      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      INCLUDE 'ccub.inc'
      
C parameters
      
      INTEGER KCOL(*),KVERT(NNVE,*),KMID(NNVE,*),KOFF(*)
      INTEGER KLD(*),KAB(2,NNAB,*),KABN(*),TRIA(SZTRIA)
      INTEGER KENTRY(NNBAS,NNBAS)
      DOUBLE PRECISION DA(*),DCORVG(2,*),COECON(NNDER,NNDER,*)
      LOGICAL BCON(*)
      INTEGER NEQ,NA,ICUB,ISYMM,NBLOC

      INTEGER IPARAM(*)
      DOUBLE PRECISION DPARAM(*)

      DOUBLE PRECISION COEFF
      EXTERNAL ELE,COEFF

C externals

      INTEGER NDFL
      EXTERNAL NDFL

C local variables

      LOGICAL BSYMM,BFIRST,BCON0
      INTEGER I,J,IBLOC,I1,IELTYP,IA,IB,ILD,ICOL,JCOL
      INTEGER IDOFE,JDOFE,JDOFE1,JCOL0,IDFG
      INTEGER IVE,JP,IALBET,JCOLB,IROW,J1,ICOL1,ICOLB,ICOL1B
      DOUBLE PRECISION AUX,DJ1,DJ2,DJ3,DJ4,XI1,XI2,XX,YY,OM,DB
      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL

      SUB='ABM7'
      IF (ICHECK.GE.997) CALL OTRC('ABM7  ','12/08/94')
C
C *** Preparation - evaluation of parameters
      IER=0
      BSYMM=ISYMM.GE.1
C *** Which derivatives of basis functions are needed?
      DO 1 I = 1,NNDER
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
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C
      BFIRST=.TRUE.
C *** Dummy call of COEFF for nonlinear problems
C *** COEFF must set BDER(IDER)=.TRUE. if derivative IDER is needed
      AUX=COEFF(0D0,0D0,-1,-1,0,BFIRST,TRIA,IPARAM,DPARAM)
      BCON0=.TRUE.
C
      DO 4 IBLOC=1,NBLOC
      IF (BCON(IBLOC)) THEN
       DO 5 I=1,KABN(IBLOC)
       IA=KAB(1,I,IBLOC)
       IB=KAB(2,I,IBLOC)
5      COECON(IA,IB,IBLOC)=COEFF(0D0,0D0,IA,IB,IBLOC,BFIRST,TRIA,
     *                           IPARAM,DPARAM)
      ELSE
       BCON0=.FALSE.
      ENDIF
4     CONTINUE
C********************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
C********************************************************************
C *** Dummy call - ELE may save arithmetic operations
      CALL ELE(0D0,0D0,-2)
      IF (IER.LT.0) GOTO 99999
C
C *** Loop over all elements
      JDOFE1=1
      DO 100 IEL=1,TRIA(ONEL)
      CALL NDFGLX(TRIA,IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
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
      DO 120 IVE = 1, TRIA(ONVE)
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
      DO 200 ICUBP = 1,NCUBP
C
C *** Cubature points on the reference element
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
C
C *** Jacobian of the bilinear mapping onto the reference element
      DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
      DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
      DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
      DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
      DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
C
C *** Cubature points on the actual element + weights
      XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *  +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
      YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1+0.5D0*
     *         (DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2
      OM=DOMEGA(ICUBP)*DETJ
C
      CALL ELE(XI1,XI2,-3)
      IF (IER.LT.0) GOTO 99999
C
C *** Summing up over all pairs of multiindices
      BFIRST=.TRUE.
      DO 300 IBLOC=1,NBLOC
C
      DO 301 IALBET=1,KABN(IBLOC)
      IA=KAB(1,IALBET,IBLOC)
      IB=KAB(2,IALBET,IBLOC)
      IF (.NOT.BCON(IBLOC)) THEN
       AUX=COEFF(XX,YY,IA,IB,IBLOC,BFIRST,TRIA,IPARAM,DPARAM)*OM
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

C     Nonparametric version

      SUBROUTINE AB7NX(DA,KCOL,KLD,NA,NEQ,NBLOC,KOFF,KVERT,KMID,
     *                DCORVG,TRIA,ELE,COEFF,BCON,COECON,
     *                KAB,KABN,ICUB,ISYMM,IPARAM,DPARAM)

      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      INCLUDE 'ccub.inc'
      
C parameters
      
      INTEGER KCOL(*),KVERT(NNVE,*),KMID(NNVE,*),KOFF(*)
      INTEGER KLD(*),KAB(2,NNAB,*),KABN(*),TRIA(SZTRIA)
      INTEGER KENTRY(NNBAS,NNBAS)
      DOUBLE PRECISION DA(*),DCORVG(2,*),COECON(NNDER,NNDER,*)
      LOGICAL BCON(*)
      INTEGER NEQ,NA,ICUB,ISYMM,NBLOC

      INTEGER IPARAM(*)
      DOUBLE PRECISION DPARAM(*)

      DOUBLE PRECISION COEFF
      EXTERNAL ELE,COEFF

C externals

      INTEGER NDFL
      EXTERNAL NDFL

C local variables

      LOGICAL BSYMM,BFIRST,BCON0
      INTEGER I,J,IBLOC,I1,IELTYP,IA,IB,ILD,ICOL,JCOL
      INTEGER IDOFE,JDOFE,JDOFE1,JCOL0,IDFG
      INTEGER IVE,JP,IALBET,JCOLB,IROW,J1,ICOL1,ICOLB,ICOL1B
      DOUBLE PRECISION AUX,DJ1,DJ2,DJ3,DJ4,XI1,XI2,XX,YY,OM,DB
      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL

      SUB='ABM7'
      IF (ICHECK.GE.997) CALL OTRC('ABM7  ','12/08/94')
C
C *** Preparation - evaluation of parameters
      IER=0
      BSYMM=ISYMM.GE.1
C *** Which derivatives of basis functions are needed?
      DO 1 I = 1,NNDER
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
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C
      BFIRST=.TRUE.
C *** Dummy call of COEFF for nonlinear problems
C *** COEFF must set BDER(IDER)=.TRUE. if derivative IDER is needed
      AUX=COEFF(0D0,0D0,-1,-1,0,BFIRST,TRIA,IPARAM,DPARAM)
      BCON0=.TRUE.
C
      DO 4 IBLOC=1,NBLOC
      IF (BCON(IBLOC)) THEN
       DO 5 I=1,KABN(IBLOC)
       IA=KAB(1,I,IBLOC)
       IB=KAB(2,I,IBLOC)
5      COECON(IA,IB,IBLOC)=COEFF(0D0,0D0,IA,IB,IBLOC,BFIRST,TRIA,
     *                           IPARAM,DPARAM)
      ELSE
       BCON0=.FALSE.
      ENDIF
4     CONTINUE
C********************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
C********************************************************************
C
C *** Loop over all elements
      JDOFE1=1
      DO 100 IEL=1,TRIA(ONEL)
      CALL NDFGLX(TRIA,IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
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
      DO 120 IVE = 1, TRIA(ONVE)
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
C *** Dummy call - ELE may save arithmetic operations
      CALL ELE(0D0,0D0,-2)
      IF (IER.LT.0) GOTO 99999
C
C *** Loop over all cubature points
      DO 200 ICUBP = 1,NCUBP
C
C *** Cubature points on the reference element
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
C
C *** Jacobian of the bilinear mapping onto the reference element
      DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
      DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
      DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
      DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
      DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
C
C *** Cubature points on the actual element + weights
      XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *  +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
      YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1+0.5D0*
     *         (DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2
      OM=DOMEGA(ICUBP)*DETJ
C
      CALL ELE(XX,YY,-3)
      IF (IER.LT.0) GOTO 99999
C
C *** Summing up over all pairs of multiindices
      BFIRST=.TRUE.
      DO 300 IBLOC=1,NBLOC
C
      DO 301 IALBET=1,KABN(IBLOC)
      IA=KAB(1,IALBET,IBLOC)
      IB=KAB(2,IALBET,IBLOC)
      IF (.NOT.BCON(IBLOC)) THEN
       AUX=COEFF(XX,YY,IA,IB,IBLOC,BFIRST,TRIA,IPARAM,DPARAM)*OM
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

************************************************************************
* XAP9                                                                 *
*                                                                      *
* Purpose  Call AP9                                                    *
*          Allocate KLD and KCOL on DWORK                              *
*                                                                      *
* Subroutines/functions called  AP9, EA00, ZNEW, ZDISP, ZFREE          *
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

      SUBROUTINE XAP9X(LCOL,LLD,NA,NEQ,TRIA,ELE1,ELE2)

      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'stria.inc'
      
C     parameters

      INTEGER LCOL,LLD,NA,NEQ,TRIA(SZTRIA)
      
      EXTERNAL ELE1,ELE2
      
C     local variables
      
      INTEGER IWMAX0,IFREE,LCOL1,LIND,ITYPE1
      INTEGER KVERT,KMID
      
      INTEGER NDFGX
      EXTERNAL NDFGX

      SUB='XAP9'
      IF (ICHECK.GE.997) CALL OTRC('XAP9  ','12/20/89')
C
C *** Determine total number of degrees of freedom
      CALL EA00(ELE1,ELE1,TRIA(ONVE),ITYPE1)
      NEQ=NDFGX(ITYPE1,TRIA)
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

      KVERT  = L(TRIA(OLVERT))
      KMID   = L(TRIA(OLMID))

      CALL AP9X(KWORK(L(LCOL)),KWORK(L(LCOL1)),
     *         KWORK(L(LIND)),KWORK(L(LLD)),
     *         NA,NEQ,TRIA,ELE1,ELE2,KWORK(KVERT),KWORK(KMID))
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


************************************************************************
*                                                                      *
* AP9X                                                                 *
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

      SUBROUTINE AP9X(KCOL,KCOL1,KIND,KLD,NA,NEQ,TRIA,ELE1,ELE2,
     *                KVERT,KMID)
      
      IMPLICIT NONE

      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'stria.inc'
      
C     parameters 

      INTEGER KCOL(*),KCOL1(*),KIND(*),KLD(*),NA,NEQ
      INTEGER KVERT(*),KMID(*),TRIA(SZTRIA)
      
      EXTERNAL ELE1,ELE2
      
C     local variables

      INTEGER IFREE,IEQ,IEL,JDOFE,IROW
      INTEGER IPOS,ICOL,JCOL,IHELP,ITYP1(3),NELE,I,JDOFP
      LOGICAL BSORT
      INTEGER KDFG1(NNBAS,3),KDFL1(NNBAS,3),IDFL1(3)
      
      INTEGER NDFL
      EXTERNAL NDFL

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
      CALL EA00(ELE1,ELE1,TRIA(ONVE),ITYP1(1))
      CALL EA00(ELE2,ELE2,TRIA(ONVE),ITYP1(2))
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
      DO 100 IEL=1,TRIA(ONEL)
C
C *** KDFG returns the global degrees of freedom in increasing order
      DO 101 I=1,NELE
      CALL NDFGLX(TRIA,IEL,0,ITYP1(I),KVERT,KMID,KDFG1(1,I),KDFL1(1,I))
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

      SUBROUTINE XAB09X(LA,LCOL,LLD,NA,NBLOC,ICLEAR,TRIA,
     *                 ELE1,ELE2,ELE3,COEFF,BCON,KAB,KABN,
     *                 ICUB,ARR,IPARAM,DPARAM)

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasicelem.inc'
      
      INCLUDE 'stria.inc'

C parameters

      INTEGER LA(*),KAB(2,NNAB,*),KABN(*)
      LOGICAL BCON(*)

      CHARACTER ARR*12
      DIMENSION ARR(*)
      
      INTEGER LCOL, LLD, NA, NBLOC, ICLEAR, ICUB
      INTEGER TRIA(SZTRIA)
      INTEGER IPARAM(*)
      DOUBLE PRECISION DPARAM(*)

      DOUBLE PRECISION COEFF
      EXTERNAL COEFF,ELE1,ELE2,ELE3

C local variables
      
      INTEGER IBLOC, ITYPE, ILEN, LOFF, LOECON
      INTEGER KVERT,KMID,KCORVG

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

      KVERT  = L(TRIA(OLVERT))
      KMID   = L(TRIA(OLMID))
      KCORVG = L(TRIA(OLCORVG))

      CALL AB09X(DWORK(1),KWORK(L(LCOL)),KWORK(L(LLD)),NA,NBLOC,
     *          KWORK(L(LOFF)),KWORK(KVERT),KWORK(KMID),
     *          DWORK(KCORVG),TRIA,ELE1,ELE2,
     *          ELE3,COEFF,BCON,DWORK(L(LOECON)),KAB,KABN,ICUB,
     *          IPARAM,DPARAM)
      IF (IER.NE.0) GOTO 99999
C
      CALL ZDISP(0,LOECON,'COECON')
      CALL ZDISP(0,LOFF,'KOFF  ')

99999 END

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

      SUBROUTINE AB09X(DA,KCOL,KLD,NA,NBLOC,KOFF,KVERT,KMID,
     *                DCORVG,TRIA,ELE1,ELE2,ELE3,COEFF,
     *                BCON,COECON,KAB,KABN,ICUB,IPARAM,DPARAM)

      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      INCLUDE 'ccub.inc'
      
C parameters
      
      INTEGER KCOL(*),KVERT(NNVE,*),KMID(NNVE,*),KOFF(*)
      INTEGER KLD(*),KAB(2,NNAB,*),KABN(*)
      INTEGER KENTRY(NNBAS,NNBAS),TRIA(SZTRIA)
      DOUBLE PRECISION DA(*),DCORVG(2,*),COECON(NNDER,NNDER,*)
      LOGICAL BCON(*)
      INTEGER NA,ICUB,NBLOC

      INTEGER IPARAM(*)
      DOUBLE PRECISION DPARAM(*)

      DOUBLE PRECISION COEFF
      EXTERNAL ELE1,ELE2,ELE3,COEFF

C externals

      INTEGER NDFL
      EXTERNAL NDFL

C local variables

      LOGICAL BFIRST,BCON0
      INTEGER I,J,IBLOC,IA,IB,JCOL,NELE,IBAS
      INTEGER IDOFE,JDOFE,JCOL0,IDFG,IDER,ITYP1(3)
      INTEGER IVE,JP,IALBET,JCOLB,IABN,IAB
      DOUBLE PRECISION AUX,DJ1,DJ2,DJ3,DJ4,XI1,XI2,XX,YY,OM,DB
      DOUBLE PRECISION DBAS1(NNBAS,NNDER,3)
      LOGICAL BDER1(NNDER,3),BELE(3)
      INTEGER KDFG1(NNBAS,3),KDFL1(NNBAS,3),IDFL1(3)

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
6      COECON(IA,IB,IBLOC)=COEFF(0D0,0D0,IA,IB,IBLOC,BFIRST,TRIA,
     *                           IPARAM,DPARAM)
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
      DO 100 IEL=1,TRIA(ONEL)
C
      DO 101 I=1,3
      IF (.NOT.BELE(I)) GOTO 101
      CALL NDFGLX(TRIA,IEL,1,ITYP1(I),KVERT,KMID,KDFG1(1,I),KDFL1(1,I))
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
      DO 120 IVE=1,TRIA(ONVE)
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
       AUX=COEFF(XX,YY,IA,IB,IBLOC,BFIRST,TRIA,IPARAM,DPARAM)*OM
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

