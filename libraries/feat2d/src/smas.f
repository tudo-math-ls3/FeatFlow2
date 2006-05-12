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
* XSMAS                                                                *
*                                                                      *
* Purpose  Save information about macrotriangulation                   *
*                                                                      *
* Subroutines/functions called  SMAS, ZCPY                             *
*                                                                      *
* Version from  04/12/91                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----   ----                                                         *
*                                                                      *
* OUTPUT  TYPE                                                         *
* ------  ----                                                         *
*                                                                      *
************************************************************************
C
      SUBROUTINE XSMAS
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNARR=299,NNVE=4)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /MACROD/ NMAEL,NMAVT,NMAEDG,NMAVE,NMAVEL,NMABCT,NMAVBD
      COMMON /MACROA/ LMACVG,LMACMG,LMAVT,LMAMID,LMAADJ,LMAVEL,LMAMEL,
     *                LMANPR,LMAMM,LMAVBD,LMAEBD,LMABCT,LMAVBP,LMAMBP,
     *                LMAVE
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/,/MACROD/,/MACROA/
C
      SUB='XSMAS'
      IF (ICHECK.GE.997) CALL OTRC('XSMAS ','04/12/91')
      IER=0
C
      NMAEL=NEL
      NMAVT=NVT
      NMAEDG=NVT+NEL+NBCT-2
      NMAVE=NVE
	NMABCT=NBCT
	NMAVBD=NVBD
C
C *** DMACVG  - Coordinates of macrovertices
      CALL ZCPY(LCORVG,'DCORVG',LMACVG,'DMACVG')
      IF (IER.NE.0) GOTO 99999
C *** DMACVG  - Coordinates of macrovertices
      IF (LCORMG.NE.0) THEN
       CALL ZCPY(LCORMG,'DCORMG',LMACMG,'DMACMG')
       IF (IER.NE.0) GOTO 99999
	ENDIF
C *** KMAVT  - Vertices of macroelements
      CALL ZCPY(LVERT,'KVERT ',LMAVT,'KMAVT ')
      IF (IER.NE.0) GOTO 99999
C *** KMAADJ - Neighbours of macroelements
      CALL ZCPY(LADJ,'KADJ  ',LMAADJ,'KMAADJ')
      IF (IER.NE.0) GOTO 99999
C *** KMANPR - Nodal properties of macrovertices
      CALL ZCPY(LNPR,'KNPR  ',LMANPR,'KMANPR')
      IF (IER.NE.0) GOTO 99999
C *** KMAMM - Macrovertices with minimum and maximum parameter
      CALL ZCPY(LMM,'KMM   ',LMAMM,'KMAMM ')
      IF (IER.NE.0) GOTO 99999
C *** KMAVBD - Macrovertices on the boundary
      IF (LVBD.GT.0) THEN
       CALL ZCPY(LVBD,'KVBD  ',LMAVBD,'KMAVBD')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C *** KMAEBD - Macroelement at the boundary
      IF (LEBD.GT.0) THEN
       CALL ZCPY(LEBD,'KEBD  ',LMAEBD,'KMAEBD')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C *** KMABCT - Pointer field
      IF (LBCT.GT.0) THEN
       CALL ZCPY(LBCT,'KBCT  ',LMABCT,'KMABCT')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C *** DMAVBP, DMAMBP - Parameter values
      IF (LVBDP.GT.0) THEN
       CALL ZCPY(LVBDP,'DVBDP ',LMAVBP,'DMAVBP')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LMBDP.GT.0) THEN
       CALL ZCPY(LMBDP,'DMBDP ',LMAMBP,'DMAMBP')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
C *** KMAMID - Numbers of macroedges for each element
      CALL ZNEW(NMAEL*NNVE,3,LMAMID,'KMAMID')
      IF (IER.NE.0) GOTO 99999
C *** KMAVE  - Numbers of endpoints for each macroedge
C ***          Smaller number first
      CALL ZNEW(2*NMAEDG,-3,LMAVE,'KMAVE ')
      IF (IER.NE.0) GOTO 99999
C
      CALL SMAS(KWORK(L(LMAVT)),KWORK(L(LMAADJ)),
     *          KWORK(L(LMAMID)),KWORK(L(LMAVE)))
C
      CALL XSMA2V
C
C
99999 END
C
C
C
      SUBROUTINE SMAS(KMAVT,KMAADJ,KMAMID,KMAVE)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=4)
      DIMENSION KMAVT(NNVE,*),KMAADJ(NNVE,*),KMAMID(NNVE,*),KMAVE(2,*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /MACROD/ NMAEL,NMAVT,NMAEDG,NMAVE,NMAVEL,NMABCT,NMAVBD
      SAVE /ERRCTL/,/CHAR/,/MACROD/,/TRIAD/
	DATA JVE2/0/
C
      SUB='SMAS'
      IF (ICHECK.GE.997) CALL OTRC('SMAS  ','04/12/91')
      IER=0
C
      NMAEDG=0
C
      DO 100 IEL=1,NMAEL
      DO 110 IVE1=1,NVE
C *** Number of following vertex
      IVE2=MOD(IVE1,NVE)+1
C *** Number of neighbouring element
      JEL=KMAADJ(IVE1,IEL)
C *** Numbers of endpoints of macroedge
      IVT1=KMAVT(IVE1,IEL)
      IVT2=KMAVT(IVE2,IEL)
C
      IF (JEL.EQ.0.OR.JEL.GE.IEL) THEN
C *** New macroedge
C *** Update number
       NMAEDG=NMAEDG+1
C *** Store numbers of endpoints of macroedge NMAEDG
       KMAVE(1,NMAEDG)=MIN(IVT1,IVT2)
       KMAVE(2,NMAEDG)=MAX(IVT1,IVT2)
C *** Store number of macroedge for element IEL
       KMAMID(IVE1,IEL)=NMAEDG
C
      ELSE
C
C *** Macroedge already available
C *** Find matching macroedge of neighbouring element
C *** Recall that vertices occur in opposite order
       DO 120 JVE1=1,NVE
       JVE2=JVE1-1
       IF (JVE2.EQ.0) JVE2=NVE
       IF (KMAVT(JVE1,JEL).EQ.KMAVT(IVE1,IEL).AND.
     *     KMAVT(JVE2,JEL).EQ.KMAVT(IVE2,IEL)) GOTO 121
120    CONTINUE
121    KMAMID(IVE1,IEL)=KMAMID(JVE2,JEL)
C
      ENDIF
C
110   CONTINUE
100   CONTINUE
C
      END
