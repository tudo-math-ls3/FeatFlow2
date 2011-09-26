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
* XSA0                                                                 *
*                                                                      *
* Purpose  Adjust dimensions of DCORVG, KVERT, KADJ and KNPR           *
*          Call of SA0                                                 *
*                                                                      *
* Subroutines/functions called  SA0, ZLEN, ZDISP, ZNEW, ZCPY           *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* IDISP    I*4    =1 Release free space on the arrays                  *
*                    after determination of the new subdivision        *
* For the description of the remaining parameters see SA0              *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
*                                                                      *
************************************************************************
*                                                                      *
* SA0                                                                  *
*                                                                      *
* Purpose  Regular refinement of a given triangular subdivision        *
*          of a two-dimensional domain                                 *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DCORVG   R*8    Cartesian coordinates of interior vertices           *
*                 Parameter of vertices on the boundary                *
* KVERT    I*4    Numbers of vertices in each element, counterclockwise*
* KADJ     I*4    Number of neighbouring element, 0 for boundary edges *
*                                                                      *
*                 Convention for numbers of vertices and neighbours    *
*                                                                      *
*                                     P3                               *
*                                      *                               *
*                                 N3  . . N2                           *
*                                    .   .                             *
*                                   *.....*                            *
*                                 P1   N1  P2                          *
*                                                                      *
* KNPR     I*4    0    for interior vertices                           *
*                 IBCT for nodal points on boundary component IBCT     *
* KMM      I*4    KMM(1,IBCT) and KMM(2,IBCT) contain the numbers      *
*                 of vertices with minimum and maximum parameter       *
*                                                                      *
* NNEL     I*4    Maximum dimension provided for KVERT                 *
* NNVT     I*4    Maximum dimension provided for DCORVG, KNPR          *
*                                                                      *
* NFINE    I*4    Number of regular subdivisions                       *
* S2DI     SUBR   S2DI(X1,Y1,X2,Y1,XT,YT) gives the coordinates        *
*                 of the new vertex on the line (X1,Y1)-(X2,Y2)        *
*                 S2DI=S2DI0 takes the midpoint                        *
* S2DB     SUBR   S2DB(X1,X2,XT,BMM,IBCT,PARX,PARY,TMAX(IBCT))         *
*                 gives the parameter of the new boundary vertex       *
*                 on the boundary segment joining the vertices with    *
*                 parameter X1 and X2                                  *
*                 BMM=.TRUE. means that X1,X2 are the minimum and      *
*                 maximum parameter on boundary component IBCT         *
* PARX     FUNC                                                        *
* PARY     FUNC   Double precision functions describing the domain     *
* TMAX     FUNC                                                        *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DCORVG   R*8    As above                                             *
* KVERT    I*4    As above                                             *
* KADJ     I*4    As above                                             *
* KNPR     I*4    As above                                             *
* KMM      I*4    As above                                             *
*                                                                      *
* NNEL     I*4    Dimension needed for KVERT, KADJ                     *
* NNVT     I*4    Dimension needed for DCORVG, KNPR                    *
*                                                                      *
* NEL      I*4    Number of elements of final triangulation            *
* NVT      I*4    Number if vertices of final triangulation            *
*                                                                      *
* IER      I*4     0  No error                                         *
*                  1  NNEL or NNVT too small                           *
*                     Requirements                                     *
*                     NNEL >=NEL(NFINE)                                *
*                     NNVT >=NVT(NFINE)                                *
*                                                                      *
************************************************************************
C
      SUBROUTINE XSA0(NFINE,IDISP,S2DI,S2DB,PARX,PARY,TMAX)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=4,NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL S2DI,S2DB,PARX,PARY,TMAX
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/
C
      SUB='XSA0'
      IF (ICHECK.GE.997) CALL OTRC('XSA0  ','01/02/89')
C
      CALL ZLEN(LVERT,ILEN)
      NNEL=ILEN/NNVE
      CALL ZLEN(LADJ,ILEN)
      NNEL1=ILEN/NNVE
      NNEL=MIN(NNEL,NNEL1)
      CALL ZLEN(LCORVG,ILEN)
      NNVT=ILEN/2
      NNEL1=NNEL
      NNVT1=NNVT
C
      CALL SA0(DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *         KWORK(L(LADJ)),KWORK(L(LNPR)),KWORK(L(LMM)),
     *         NFINE,NNEL,NNVT,S2DI,S2DB,PARX,PARY,TMAX)
C
      IF (IER) 99999,1,2
C
1     IF (IDISP.EQ.1) THEN
       CALL ZDISP(NNVE*NNEL,LVERT,'KVERT ')
       CALL ZDISP(NNVE*NNEL,LADJ,'KADJ  ')
       CALL ZDISP(2*NNVT,LCORVG,'DCORVG')
       CALL ZDISP(NNVT,LNPR,'KNPR  ')
      ENDIF
      GOTO 99999
C
2     IF (NNEL.GT.NNEL1) THEN
       CALL ZNEW(NNVE*NNEL,3,LV1,'KVERT ')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LVERT,'KVOLD ',LV1,'KVERT ')
       CALL ZDISP(0,LVERT,'KVOLD ')
       LVERT=LV1
C
       CALL ZNEW(NNVE*NNEL,3,LA1,'KADJ  ')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LADJ,'KAOLD ',LA1,'KADJ  ')
       CALL ZDISP(0,LADJ,'KAOLD ')
       LADJ=LA1
      ENDIF
C
      IF (NNVT.GT.NNVT1) THEN
       CALL ZNEW(2*NNVT,1,LCV1,'DCORVG')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LCORVG,'DCVOLD',LCV1,'DCORVG')
       CALL ZDISP(0,LCORVG,'DCVOLD')
       LCORVG=LCV1
      ENDIF
C
      IF (NNVT.GT.NNVT1) THEN
      CALL ZNEW(NNVT,3,LNPR1,'KNPR  ')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LNPR,'KNPOLD',LNPR1,'KNPR  ')
       CALL ZDISP(0,LNPR,'KNPOLD')
       LNPR=LNPR1
      ENDIF
C
      CALL SA0(DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *         KWORK(L(LADJ)),KWORK(L(LNPR)),KWORK(L(LMM)),
     *         NFINE,NNEL,NNVT,S2DI,S2DB,PARX,PARY,TMAX)
C
99999 END
C
C
C
      SUBROUTINE SA0(DCORVG,KVERT,KADJ,KNPR,KMM,
     *               NFINE,NNEL,NNVT,S2DI,S2DB,PARX,PARY,TMAX)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNVE=4)
      DIMENSION DCORVG(2,*)
      DIMENSION KVERT(NNVE,*),KADJ(NNVE,*),KNPR(*),KMM(2,*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      EXTERNAL PARX,PARY
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
      DATA JEL1/0/,IVT/0/
C
      SUB='SA0'
      IF (ICHECK.GE.997) CALL OTRC('SA0   ','01/02/89')
C
      NNEL1=NNEL
      NNVT1=NNVT
      NNEL=NEL
      NNVT=NVT
      NNMT=NNEL+NNVT+NBCT-2
C
      IF (NFINE.GT.0) THEN
       DO 1 IFINE=1,NFINE
       NNEL=4*NNEL
       NNVT=NNVT+NNMT
1      NNMT=NNVT+NNEL+NBCT-2
      ENDIF
C
      IF (NNEL.GT.NNEL1.OR.NNVT.GT.NNVT1) THEN
       IER=1
       GOTO 99999
      ENDIF
C
      IER=0
C
      DO 100 IFINE=1,NFINE
      NVT1=NVT
      NEL1=NEL
      DO 110 IEL=1,NEL1
C
C Calculate coordinates of parameter of new vertices
C
      DO 120 IVE1=1,NVE
C
      IF (KADJ(IVE1,IEL).LT.IEL.AND.KADJ(IVE1,IEL).NE.0) GOTO 120
      NVT1=NVT1+1
      IVT1=KVERT(IVE1,IEL)
      IBCT1=KNPR(IVT1)
      IVE2=MOD(IVE1,NVE)+1
      IVT2=KVERT(IVE2,IEL)
      IBCT2=KNPR(IVT2)
C
      IF (KADJ(IVE1,IEL).NE.0) THEN
       IF (IBCT1.EQ.0) THEN
        X1=DCORVG(1,IVT1)
        Y1=DCORVG(2,IVT1)
       ELSE
        X1=PARX(DCORVG(1,IVT1),IBCT1)
        Y1=PARY(DCORVG(1,IVT1),IBCT1)
       ENDIF
       IF (IBCT2.EQ.0) THEN
        X2=DCORVG(1,IVT2)
        Y2=DCORVG(2,IVT2)
       ELSE
        X2=PARX(DCORVG(1,IVT2),IBCT2)
        Y2=PARY(DCORVG(1,IVT2),IBCT2)
       ENDIF
C
       CALL S2DI(X1,Y1,X2,Y2,DCORVG(1,NVT1),DCORVG(2,NVT1))
       KNPR(NVT1)=0
C
      ELSE
C
       BMM=(IVT1.EQ.KMM(1,IBCT1).AND.IVT2.EQ.KMM(2,IBCT1)).OR.
     *     (IVT2.EQ.KMM(1,IBCT1).AND.IVT1.EQ.KMM(2,IBCT1))
       TTMAX=TMAX(IBCT1)
       CALL S2DB(DCORVG(1,IVT1),DCORVG(1,IVT2),DCORVG(1,NVT1),
     *           BMM,IBCT1,PARX,PARY,TTMAX)
       KNPR(NVT1)=IBCT1
C
       IF (BMM) THEN
        BMM=.FALSE.
        IF (DCORVG(1,NVT1).GT.DCORVG(1,KMM(2,IBCT1))) THEN
         KMM(2,IBCT1)=NVT1
        ELSE
         KMM(1,IBCT1)=NVT1
        ENDIF
       ENDIF
C
      ENDIF
120   CONTINUE
C
C *** Actualization of KVERT
C
      DO 130 IVE1=1,NVE
130   KVERT(1,NEL+IVE1)=KVERT(IVE1,IEL)
C
      DO 140 IVE1=1,NVE
      IVE2=MOD(IVE1,NVE)+1
      IEL1=KADJ(IVE1,IEL)
C
      IF (IEL1.GE.IEL.OR.IEL1.EQ.0) THEN
       NVT=NVT+1
       KVERT(2,NEL+IVE1)=NVT
       KVERT(3,NEL+IVE2)=NVT
       KVERT(IVE1,IEL)=NVT
C
      ELSE
C
C *** Find out the number of the vertex (already exists)
C
       DO 141 IVE=1,NVE
       JEL=KADJ(IVE,IEL1)
       DO 141 JVE=1,NVE
       IF (KADJ(JVE,JEL).EQ.IEL) THEN
        IF (JVE.EQ.NVE) THEN
         IVT=KVERT(NVE,JEL)
        ELSE
         IVT=KVERT(2,JEL)
        ENDIF
        GOTO 145
       ENDIF
141    CONTINUE
C
145    KVERT(2,NEL+IVE1)=IVT
       KVERT(3,NEL+IVE2)=IVT
       KVERT(IVE1,IEL)=IVT
      ENDIF
C
140   CONTINUE
C
C *** Actualization of KADJ
C
      DO 150 IVE1=1,NVE
      IVE2=MOD(IVE1,NVE)+1
      JEL=KADJ(IVE1,IEL)
      IF (JEL.EQ.0.OR.JEL.GE.IEL) THEN
       KADJ(1,NEL+IVE1)=JEL
       KADJ(NVE,NEL+IVE2)=JEL
      ELSE
       DO 151 JVE=1,NVE
       JEL1=KADJ(JVE,JEL)
       IF (KADJ(1,JEL1).EQ.IEL) GOTO 152
151    CONTINUE
152    JEL=JEL1+1
       IF (JVE.EQ.2) JEL=JEL1-2
       KADJ(1,NEL+IVE1)=JEL
       KADJ(NVE,NEL+IVE2)=JEL1
       KADJ(1,JEL1)=NEL+IVE2
       KADJ(NVE,JEL)=NEL+IVE1
      ENDIF
150   CONTINUE
C
      NEL=NEL+3
      KADJ(1,IEL)=NEL-1
      KADJ(2,IEL)=NEL
      KADJ(3,IEL)=NEL-2
      KADJ(2,NEL-1)=IEL
      KADJ(2,NEL)  =IEL
      KADJ(2,NEL-2)=IEL
C
110   CONTINUE
100   CONTINUE
C
99999 END
