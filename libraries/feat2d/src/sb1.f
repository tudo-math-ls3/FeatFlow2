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
* XSB1                                                                 *
*                                                                      *
* Purpose  Adjust dimensions of DCORVG, KVERT, KADJ and KNPR           *
*          Call of SB1                                                 *
*                                                                      *
* Subroutines/functions called  SB1, ZLEN, ZDISP, ZNEW, ZCPY           *
*                                                                      *
* Version from  04/12/91                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* IDISP    I*4    =1 Release free space on the arrays                  *
*                    after determination of the new subdivision        *
* For the description of the remaining parameters see SB0              *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
*                                                                      *
************************************************************************
*                                                                      *
* SB1                                                                  *
*                                                                      *
* Purpose  Regular refinement of a given quadrilateral subdivision     *
*          of a two-dimensional domain                                 *
*          Each macroelement subdivided into NFINE**2 quadrilaterals   *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  04/12/91                                               *
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
* NFINE    I*4    NFINE**2 elements in each macroquadrilateral         *
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
* Conventions                                                          *
* -----------                                                          *
* Numbers of vertices - 1. Macrovertices                               *
*                       2. Interior vertices of macroedges             *
*                       3. Interior vertices of macroelements          *
* Numbers of elements - (I-1)*NFINE**2+1,...,I*NFINE**2 for I-th       *
*                       macroelement                                   *
*                                                                      *
************************************************************************
C
      SUBROUTINE XSB1(NFINE,IDISP,PARX,PARY,TMAX)
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
      COMMON /MACROD/ NMAEL,NMAVT,NMAEDG,NMAVE,NMAVEL,NMABCT,NMAVBD
      COMMON /MACROA/ LMACVG,LMACMG,LMAVT,LMAMID,LMAADJ,LMAVEL,LMAMEL,
     *                LMANPR,LMAMM,LMAVBD,LMAEBD,LMABCT,LMAVBP,LMAMBP,
     *                LMAVE
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL PARX,PARY,TMAX
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/,/MACROA/,/MACROD/
C
      SUB='XSB1'
      IF (ICHECK.GE.997) CALL OTRC('XSB1  ','04/12/91')
C
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
      CALL ZNEW(NNVE*(NFINE+1),-3,LVAUX,'KVAUX ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(2*(NFINE+1),-3,LLU,'KLU   ')
      IF (IER.NE.0) GOTO 99999
C
      CALL SB1(KWORK(L(LMAVT)),KWORK(L(LMAADJ)),KWORK(L(LMAMID)),
     *         KWORK(L(LMAVE)),DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *         KWORK(L(LADJ)),KWORK(L(LNPR)),KWORK(L(LMM)),
     *         NFINE,NNEL,NNVT,PARX,PARY,TMAX,
     *         KWORK(L(LVAUX)),KWORK(L(LLU)))
C
      IF (IER) 99999,1,2
C
1     CALL ZDISP(0,LLU,'KLU   ')
      IF (IER.NE.0) GOTO 99999
      CALL ZDISP(0,LVAUX,'KVAUX ')
      IF (IER.NE.0) GOTO 99999
      IF (IDISP.EQ.1) THEN
       CALL ZDISP(NNVE*NNEL,LVERT,'KVERT ')
       CALL ZDISP(NNVE*NNEL,LADJ,'KADJ  ')
       CALL ZDISP(2*NNVT,LCORVG,'DCORVG')
       CALL ZDISP(NNVT,LNPR,'KNPR  ')
      ENDIF
      GOTO 1000
C
2     CALL ZDISP(0,LLU,'KLU   ')
      IF (IER.NE.0) GOTO 99999
      CALL ZDISP(0,LVAUX,'KVAUX ')
      IF (IER.NE.0) GOTO 99999
C
      IF (NNEL.GT.NNEL1) THEN
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
      CALL ZNEW(NNVE*(NFINE+1),-3,LVAUX,'KVAUX ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(2*(NFINE+1),-3,LLU,'KLU   ')
      IF (IER.NE.0) GOTO 99999
C
C
      CALL SB1(KWORK(L(LMAVT)),KWORK(L(LMAADJ)),KWORK(L(LMAMID)),
     *         KWORK(L(LMAVE)),DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *         KWORK(L(LADJ)),KWORK(L(LNPR)),KWORK(L(LMM)),
     *         NFINE,NNEL,NNVT,PARX,PARY,TMAX,
     *         KWORK(L(LVAUX)),KWORK(L(LLU)))
C
      CALL ZDISP(0,LLU,'KLU   ')
      IF (IER.NE.0) GOTO 99999
      CALL ZDISP(0,LVAUX,'KVAUX ')
      IF (IER.NE.0) GOTO 99999
      CALL ZDISP(0,LADJ,'KADJ  ')
      IF (IER.NE.0) GOTO 99999
C
C *** Determine KADJ for the new mesh
C
1000  CALL XS2A
C
99999 END
C
C
C
      SUBROUTINE SB1(KMAVT,KMAADJ,KMAMID,KMAVE,DCORVG,KVERT,KADJ,
     *               KNPR,KMM,NFINE,NNEL,NNVT,PARX,PARY,TMAX,
     *               KVAUX,KLU)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNVE=4)
      DIMENSION KMAVT(NNVE,*),KMAADJ(NNVE,*),KMAVE(2,*),DCORVG(2,*)
      DIMENSION KMAMID(NNVE,*)
      DIMENSION KVERT(NNVE,*),KADJ(NNVE,*),KNPR(*),KMM(2,*)
      DIMENSION XX(NNVE),YY(NNVE),KV(NNVE),KVAUX(NNVE,*),KLU(2,*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /MACROD/ NMAEL,NMAVT,NMAEDG,NMAVE,NMAVEL,NMABCT,NMAVBD
      EXTERNAL PARX,PARY
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/MACROD/
      DATA XT1/0D0/,XT2/0D0/
C
      SUB='SB1'
      IF (ICHECK.GE.997) CALL OTRC('SB1   ','04/12/91')
C
C
C *** Step 0 - Determine size of arrays needed for final mesh
      NNEL1=NNEL
      NNVT1=NNVT
      NNEL=NMAEL*NFINE**2
      NNVT=NMAVT+NMAEDG*(NFINE-1)+NMAEL*(NFINE-1)**2
C
C *** Return to XSB1 to adjust array dimensions
      IF (NNEL.GT.NNEL1.OR.NNVT.GT.NNVT1) THEN
       IER=1
       GOTO 99999
      ENDIF
C
      IER=0
C
C
C *** Step 1 (Loop 100)
C *** Calculate numbers and coordinates of new vertices on macroedges
C
      NVT=NMAVT
      DO 100 IMAEL=1,NMAEL
      DO 100 IVE=1,NVE
      IMAADJ=KMAADJ(IVE,IMAEL)
      IF (IMAADJ.NE.0.AND.IMAADJ.LT.IMAEL) GOTO 100
      IMAEDG=KMAMID(IVE,IMAEL)
C
C *** Determine numbers of endpoints of macroedge
      IVT1=KMAVE(1,IMAEDG)
      IVT2=KMAVE(2,IMAEDG)
C *** Find coordinates and/or parameters of endpoints
      IBCT1=KNPR(IVT1)
      IF (IBCT1.LE.0) THEN
       X1=DCORVG(1,IVT1)
       Y1=DCORVG(2,IVT1)
      ELSE
       XT1=DCORVG(1,IVT1)
       X1=PARX(XT1,IBCT1)
       Y1=PARY(XT1,IBCT1)
      ENDIF
      IBCT2=KNPR(IVT2)
      IF (IBCT2.LE.0) THEN
       X2=DCORVG(1,IVT2)
       Y2=DCORVG(2,IVT2)
      ELSE
       XT2=DCORVG(1,IVT2)
       X2=PARX(XT2,IBCT2)
       Y2=PARY(XT2,IBCT2)
      ENDIF
C
C *** Case 1 - Boundary edge
C
      IF (IMAADJ.EQ.0) THEN
       BMM=(IVT1.EQ.KMM(1,IBCT1).AND.IVT2.EQ.KMM(2,IBCT1)).OR.
     *     (IVT2.EQ.KMM(1,IBCT1).AND.IVT1.EQ.KMM(2,IBCT1))
       TTMAX=TMAX(IBCT1)
       IF (BMM) THEN
        XS1=MAX(XT1,XT2)
        H=(MIN(XT1,XT2)-XS1+TTMAX)/NFINE
       ELSE
        XS1=MIN(XT1,XT2)
        H=(MAX(XT1,XT2)-XS1)/NFINE
       ENDIF
C
       DO 101 IV=1,NFINE-1
       NVT=NVT+1
       DCORVG(2,NVT)=0D0
       KNPR(NVT)=IBCT1
       DCORVG(1,NVT)=XS1+H
C *** Update values of KMM, if necessary
       IF (BMM) THEN
        IF (DCORVG(1,NVT).GE.TTMAX) THEN
         DCORVG(1,NVT)=DCORVG(1,NVT)-TTMAX
         KMM(1,IBCT1)=NVT
         BMM=.FALSE.
        ELSE
         KMM(2,IBCT1)=NVT
        ENDIF
       ENDIF
       XS1=DCORVG(1,NVT)
101    CONTINUE
C
C *** Case 2 - Interior edge
C
      ELSE
C
       HX=(X2-X1)/NFINE
       HY=(Y2-Y1)/NFINE
       DO 102 IV=1,NFINE-1
       NVT=NVT+1
       X1=X1+HX
       Y1=Y1+HY
       DCORVG(1,NVT)=X1
       DCORVG(2,NVT)=Y1
       KNPR(NVT)=0
102   CONTINUE
      ENDIF
100   CONTINUE
C
C
C *** Step 2 (Loop 200)
C *** Calculate numbers and coordinates of new vertices in macroelements
C
      DO 200 IMAEL=1,NMAEL
C
C *** Store numbers and coordinates of macrovertices on auxiliary arrays
      DO 201 IVE=1,NVE
      KV(IVE)=KMAVT(IVE,IMAEL)
      IBCT=KNPR(KV(IVE))
      IF (IBCT.LE.0) THEN
       XX(IVE)=DCORVG(1,KV(IVE))
       YY(IVE)=DCORVG(2,KV(IVE))
      ELSE
       XX(IVE)=PARX(DCORVG(1,KV(IVE)),IBCT)
       YY(IVE)=PARY(DCORVG(1,KV(IVE)),IBCT)
      ENDIF
201   CONTINUE
C
C *** Determine increments in XI2 on left and right edge
      HX1=(XX(4)-XX(1))/NFINE
      HY1=(YY(4)-YY(1))/NFINE
      HX2=(XX(3)-XX(2))/NFINE
      HY2=(YY(3)-YY(2))/NFINE
C *** Loop over (NFINE-1)**2 interior vertices
      DO 202 IV1=1,NFINE-1
C *** Left and right lower edge on next XI2-level
      XX(1)=XX(1)+HX1
      YY(1)=YY(1)+HY1
      XX(2)=XX(2)+HX2
      YY(2)=YY(2)+HY2
C *** Stepsize in XI1 direction
      HX=(XX(2)-XX(1))/NFINE
      HY=(YY(2)-YY(1))/NFINE
      X1=XX(1)
      Y1=YY(1)
      DO 202 IV2=1,NFINE-1
      NVT=NVT+1
C *** Proceeding in XI1-direction
      X1=X1+HX
      Y1=Y1+HY
      DCORVG(1,NVT)=X1
      DCORVG(2,NVT)=Y1
      KNPR(NVT)=0
202   CONTINUE
C
200   CONTINUE
C
C
C
C *** Step 3 (Loop 300)
C *** Calculate vertices of new elements
C
      NEL=0
C
      DO 300 IMAEL=1,NMAEL
C
C *** NFINE**2 elements must be determined in each macroelement
C
C *** Step 3.1 - Store vertices on corresponding macroedges on
C ***            auxiliary array KVAUX
C
      DO 301 IVE=1,NVE
C *** Numbers of macroedges
      IMAEDG=KMAMID(IVE,IMAEL)
C *** Macro-vertices of element on macroedge
      KVAUX(IVE,1)=KMAVT(IVE,IMAEL)
      KVAUX(IVE,NFINE+1)=KMAVT(MOD(IVE,NVE)+1,IMAEL)
C
C *** Vertices on macroedge IMAEDG have the numbers
C *** NMAVT+(IMAEDG-1)*(NFINE-1)+1,...,...+NFINE-1
      IV1=NMAVT+(IMAEDG-1)*(NFINE-1)
      IF (KMAADJ(IVE,IMAEL).EQ.0) THEN
C
C *** Boundary macro-edge
C *** New vertices on macroedge are ordered with respect
C *** to the parameter (see loop 100)
       DO 310 IV=1,NFINE-1
310    KVAUX(IVE,IV+1)=IV1+IV
C
      ELSE
C *** Interior macroedge
C *** New vertices are ordered from smaller to larger
C *** number of macrovertex
C
       IF (KMAVT(IVE,IMAEL).EQ.KMAVE(1,IMAEDG)) THEN
        DO 311 IV=1,NFINE-1
311     KVAUX(IVE,IV+1)=IV1+IV
       ELSE
        DO 312 IV=1,NFINE-1
312     KVAUX(IVE,NFINE-IV+1)=IV1+IV
       ENDIF
      ENDIF
C
301   CONTINUE
C
C
C *** Outer loop over levels in XI2-direction
C *** Inner loop over levels in XI1-direction
C
C *** First, extract the numbers of vertices for the lower
C *** and upper XI2-level
C
C *** Interior vertices have the numbers
C *** NMAVT+NMAEDG*(NFINE-1)+(IMAEL-1)*(NFINE**2+1,...
C
      IV1=NMAVT+NMAEDG*(NFINE-1)+(IMAEL-1)*(NFINE-1)**2
C
      DO 320 ILEV=1,NFINE
C *** Number of vertices on upper level
C *** Find numbers on lower level (copy numbers from upper level)
      IF (ILEV.EQ.1) THEN
       DO 321 IV=1,NFINE+1
321    KLU(1,IV)=KVAUX(1,IV)
      ELSE
       DO 322 IV=1,NFINE+1
322    KLU(1,IV)=KLU(2,IV)
      ENDIF
C
C *** Find numbers on upper level
      IF (ILEV.EQ.NFINE) THEN
       DO 323 IV=1,NFINE+1
323    KLU(2,NFINE-IV+2)=KVAUX(3,IV)
      ELSE
       KLU(2,1)=KVAUX(4,NFINE-ILEV+1)
       KLU(2,NFINE+1)=KVAUX(2,ILEV+1)
       DO 324 IV=2,NFINE
       IV1=IV1+1
324    KLU(2,IV)=IV1
      ENDIF
C
C
C *** Determine KVERT for elements on level ILEV
C
      DO 325 ILU=1,NFINE
      NEL=NEL+1
      KVERT(1,NEL)=KLU(1,ILU)
      KVERT(2,NEL)=KLU(1,ILU+1)
      KVERT(3,NEL)=KLU(2,ILU+1)
      KVERT(4,NEL)=KLU(2,ILU)
325   CONTINUE
C
320   CONTINUE
C
300   CONTINUE
      NVT=NNVT
C
99999 END
