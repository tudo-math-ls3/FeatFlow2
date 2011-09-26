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
* XSC0                                                                 *
*                                                                      *
* Purpose  Adjust dimensions of DCORVG, KVERT, KADJ and KNPR           *
*          Call of SC0                                                 *
*                                                                      *
* Subroutines/functions called  SC0, ZLEN, ZDISP, ZNEW, ZCPY           *
*                                                                      *
* Version from  27/01/07 (M.Möller)                                    *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* IDISP    I*4    =1 Release free space on the arrays                  *
*                    after determination of the new subdivision        *
* For the description of the remaining parameters see SC0              *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
*                                                                      *
************************************************************************
*                                                                      *
* SC0                                                                  *
*                                                                      *
* Purpose  Regular refinement of a given mixed triangular and          *
*          quadrilateral subdivision of a two-dimensional domain       *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  27/01/07 (M.Möller)                                    *
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
*                 For triangular elements:                             *
*                                                                      *
*                                     P3                               *
*                                      *                               *
*                                 N3  . . N2                           *
*                                    .   .                             *
*                                   *.....*                            *
*                                 P1   N1  P2                          *
*                                                                      *
*                                                                      *
*                 For quadrilateral elements:                          *
*                                                                      *
*                             P4     N3      P3                        *
*                                * - - - - *                           *
*                                !         !                           *
*                             N4 !         ! N2                        *
*                                !         !                           *
*                                * - - - - *                           *
*                             P1     N1      P2                        *
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
      SUBROUTINE XSC0(NFINE,IDISP,S2DI,S2DB,PARX,PARY,TMAX)
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
      SUB='XSC0'
      IF (ICHECK.GE.997) CALL OTRC('XSC0  ','27/01/07')
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
      CALL SC0(DWORK(L(LCORVG)),KWORK(L(LVERT)),
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
      CALL SC0(DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *         KWORK(L(LADJ)),KWORK(L(LNPR)),KWORK(L(LMM)),
     *         NFINE,NNEL,NNVT,S2DI,S2DB,PARX,PARY,TMAX)
C
99999 END
C
C
C
      SUBROUTINE SC0(DCORVG,KVERT,KADJ,KNPR,KMM,
     *               NFINE,NNEL,NNVT,S2DI,S2DB,PARX,PARY,TMAX)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=4)
      DIMENSION DCORVG(2,*),KEL(NNVE),KEL1(NNVE),X(NNVE),Y(NNVE)
      DIMENSION KVERT(NNVE,*),KADJ(NNVE,*),KNPR(*),KMM(2,*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      EXTERNAL PARX,PARY
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
      DATA IVT/0/
C
      SUB='SC0'
      IF (ICHECK.GE.997) CALL OTRC('SC0   ','01/02/89')
C
      NNEL1=NNEL
      NNVT1=NNVT
      NNEL=NEL
      NNVT=NVT
      NNMT=NNEL+NNVT+NBCT-2
C
C *** Determine number of quadrilaterals
C
      NQUAD=0
      DO 10 IEL=1,NEL
      IF (KVERT(4,IEL).NE.0) THEN
       NQUAD=NQUAD+1
      ENDIF
10    CONTINUE
      IF (NFINE.GT.0) THEN
       DO 1 IFINE=1,NFINE
       NNEL=4*NNEL
       NNVT=NNVT+NNMT+NQUAD
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
C *** Determine type of element: triangle/quadrilateral
C
      IF (KVERT(4,IEL).EQ.0) THEN
       NVE=3
       IOFF=2
      ELSE
       NVE=4
       IOFF=1
      ENDIF
C
C *** Calculate coordinates of parameter of new vertices
C
      DO 120 IVE1=1,NVE
C
C *** Neighboring elements which have been refined previously, i.e.,
C *** those which have a smaller number, and boundary components
C *** can be neglected.
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
C
C *** Current edge is interior edge
C
       IF (IBCT1.EQ.0) THEN
        X1=DCORVG(1,IVT1)
        Y1=DCORVG(2,IVT1)
       ELSE
        X1=PARX(DCORVG(1,IVT1),IBCT1)
        Y1=PARY(DCORVG(1,IVT1),IBCT1)
       ENDIF
C
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
C *** Current edge is boundary edge
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
C *** Now that the vertex has been added, the vertex surrounding 
C *** element list KVERT must be updated.
C
      KEL(1)=IEL
      DO 130 I=1,3
130   KEL(I+1)=NEL+I
C
C *** Set first vertex for each subelement
C
      DO 131 I=IOFF,4
131   KVERT(1,KEL(I))=KVERT(I-IOFF+1,IEL)
C
C *** Loop over all edges connecting element IEL with its neighbors IEL1
C
      DO 140 IVE1=1,NVE
      IVE2=MOD(IVE1,NVE)+1
      IEL1=KADJ(IVE1,IEL)
C
C *** Check if adjacent element IEL1 needs to be refined
C
      IF (IEL1.GE.IEL.OR.IEL1.EQ.0) THEN
C
C *** If yes, then add new vertex
C
       NVT=NVT+1
       IVT=NVT

      ELSE
C
C *** Otherwise, determine type of adjacent element and find
C *** the vertex number of the already existing midpoint
C
       IF (KVERT(4,IEL1).EQ.0) THEN
C
C *** Neighboring element is triangle
C
        DO 141 IVE=1,3
        JEL=KADJ(IVE,IEL1)
        IF(KADJ(1,JEL).EQ.IEL) THEN
         IVT=KVERT(2,JEL)
         GOTO 145
        ENDIF
141     CONTINUE
       ELSE
C
C *** Neighboring element is quadrilateral
C
        KEL1(1)=IEL1
        DO 142 I=1,3
142     KEL1(I+1)=KADJ(2,KEL1(I))
C
        DO 143 IVE=1,4
        IF(KADJ(4,KEL1(IVE)).EQ.IEL) THEN
         IVT=KVERT(4,KEL1(IVE))
         GOTO 145
        ENDIF
143     CONTINUE
       ENDIF
      ENDIF
C
C *** Update KVERT for element IEL
C
145   IF(NVE.EQ.3) THEN
       KVERT(2,NEL+IVE1)=IVT
       KVERT(3,NEL+IVE2)=IVT
       KVERT(IVE1,IEL)  =IVT
      ELSE
       KVERT(2,KEL(IVE1))=IVT
       KVERT(4,KEL(IVE2))=IVT
      ENDIF      
C
140   CONTINUE
C
C Now that the vertex surrounding element list KVERT has been
C updated the list of adjacent elements KADJ needs to be modified.
C Loop over all edges connecting element IEL with is neighbor IEL1
C
      DO 150 IVE1=1,NVE
      IVE2=MOD(IVE1,NVE)+1
      IEL1=KADJ(IVE1,IEL)
C
C *** Check if adjacent element IEL1 needs to be refined
C
      IF (IEL1.GE.IEL.OR.IEL1.EQ.0) THEN
C
C *** If yes, then add new elements
       IF (NVE.EQ.3) THEN
C
C *** Element IEL is triangle
C
        KADJ(1,NEL+IVE1)=IEL1
        KADJ(3,NEL+IVE2)=IEL1
       ELSE
C
C *** Element IEL is quadrilateral
C
        KADJ(1,KEL(IVE1))=IEL1
        KADJ(4,KEL(IVE2))=IEL1
       ENDIF
      ELSE
C
C *** Otherwise, determine type of adjacent element and find
C *** the number of the already existing element neighbor
       IF(KVERT(4,IEL1).EQ.0) THEN
C
C *** Neighboring element IEL1 is triangle
C
        NVE1=3
        DO 151 IVE=1,3
151     KEL1(IVE)=KADJ(IVE,IEL1)

        DO 152 IVE=1,3
        IF(KADJ(3,KEL1(IVE)).EQ.IEL) THEN
         IF(IVE.EQ.1) THEN
          IVE3=3
         ELSE
          IVE3=IVE-1
         ENDIF
         GOTO 155
        ENDIF
152     CONTINUE
       ELSE
C
C *** Neighboring element IEL1 is quadrilateral
C
        NVE1=4
C
        KEL1(1)=IEL1
        DO 153 I=1,3
153     KEL1(I+1)=KADJ(2,KEL1(I))
         
        DO 154 IVE=1,4
        IF(KADJ(4,KEL1(IVE)).EQ.IEL) THEN
         IF(IVE.EQ.1) THEN
          IVE3=4
         ELSE
          IVE3=IVE-1
         ENDIF
         GOTO 155
        ENDIF
154     CONTINUE
       ENDIF
C
C *** Update KADJ
C
155    KADJ(NVE1,KEL1(IVE))      =KEL(IVE1+4-NVE)
       KADJ(1,   KEL(IVE1+4-NVE))=KEL1(IVE)
       KADJ(1,   KEL1(IVE3))     =KEL(IVE2+4-NVE)
       KADJ(NVE,KEL(IVE2+4-NVE)) =KEL1(IVE3)
C
      ENDIF  
150   CONTINUE
C
C *** Update number of elements
C
      NEL=NEL+3
C
C *** Update internal adjacency list
C
      IF(NVE.EQ.3) THEN
       KADJ(1,IEL)  =NEL-1
       KADJ(2,IEL)  =NEL
       KADJ(3,IEL)  =NEL-2
       KADJ(2,NEL-1)=IEL
       KADJ(2,NEL)  =IEL
       KADJ(2,NEL-2)=IEL
      ELSE
       DO 156 I=1,4
       I1=MOD(I,4)+1
       KADJ(2,KEL(I)) =KEL(I1)
       KADJ(3,KEL(I1))=KEL(I)
156    CONTINUE
      ENDIF
C
110   CONTINUE
C
C *** Calculation of coordinates and numbers of midpoints
C *** but only for quadrilateral elements
C
      DO 160 IEL=1,NEL1
C
C *** Check if element IEL is quadrilateral
C
      IF(KVERT(4,IEL).EQ.0) GOTO 160
C
      KEL(1)=IEL
      DO 161 I=1,3
161   KEL(I+1)=KADJ(2,KEL(I))
C
      DO 162 IVE=1,4
      IVT=KVERT(2,KEL(IVE))
      IBCT=KNPR(IVT)
      IF (IBCT.EQ.0) THEN
       X(IVE)=DCORVG(1,IVT)
       Y(IVE)=DCORVG(2,IVT)
      ELSE
       X(IVE)=PARX(DCORVG(1,IVT),IBCT)
       Y(IVE)=PARY(DCORVG(1,IVT),IBCT)
      ENDIF
162   CONTINUE
C
      NVT=NVT+1
      X21=X(2)-X(1)
      X31=X(3)-X(1)
      X24=X(2)-X(4)
      Y21=Y(2)-Y(1)
      Y31=Y(3)-Y(1)
      Y24=Y(2)-Y(4)
      ALPHA=(X21*Y24-Y21*X24)/(X31*Y24-Y31*X24)
      DCORVG(1,NVT)=X(1)+ALPHA*X31
      DCORVG(2,NVT)=Y(1)+ALPHA*Y31
      KNPR(NVT)=0
C
      DO 163 I=1,4
163   KVERT(3,KEL(I))=NVT
160   CONTINUE
C
100   CONTINUE
C
99999 END
