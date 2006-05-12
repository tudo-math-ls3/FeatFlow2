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
* XSB1X                                                                *
*                                                                      *
* Purpose  Call of the following subdivision routines                  *
*                                                                      *
* Subroutines/functions called  XSAMS, XS2C, XSB1, XS2M, XS2V, ZDISP   *
*                               XSVEB                                  *
*                                                                      *
* Version from  04/12/91                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* ICHECK   I*4    =0 Skip check of subdivision                         *
* NFINE    I*4    Each macroelement subdivided into NFINE**2 triangles *
* IMID     I*4    >=1 Determine numbers of midpoints                   *
*                 >=2 Determine coordinates of midpoints               *
* IADJ     I*4    =0 Release array KADJ on return                      *
*                    after determination of the new subdivision        *
* IVEL     I*4    =1 Determine numbers of elements meeting at each vert*
* IDISP    I*4    =1 Release free space on all arrays after completion *
* IBDP     I*4    =-1 Adjust only NVBD and KBCT, (CALL XSVEB(-1,TMAX)) *
*                    after completion                                  *
*                 >=0 Store numbers of boundary vertices on KVBD       *                       *                    (CALL XSVEB(0,TMAX))                              *
*                 >=1 Calculate sorted arrays KVBD and KEBD            *
*                    (CALL XSVEB(1,TMAX))                              *
*                 >=2 Store parameter values for boundary vertices     *
*                     (and midpoints) on DVBDP and DMBDP               *
* ISMAS    I*4    =0  use current macro-subdivision                    *
*                 >=1 accept current subdivision as macro-subdivision  *
*                 = 2 save current subdivision on file CFILE           *
*                 =-1 use previously saved subdivision from scratch    *
*                     file as macro-subdivision                        *
* MFILE    I*4                                                         *
* CFILE    C*(*)   Parameters for XOWS and XORS                        *
* IFMT     I*4                                                         *
*                                                                      *
*                                                                      *
* For the description of the remaining parameters see SB1              *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DCORVG   R*8   Cartesian coordinates of vertices                     *                                                                      *                                                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE XSB1X(NFINE,IMID,IADJ,IVEL,IDISP,IBDP,
     *                 ISMAS,MFILE,CFILE,IFMT,PARX,PARY,TMAX)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER CFILE*(*)
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /MACROD/ NMAEL,NMAVT,NMAEDG,NMAVE,NMAVEL,NMABCT,NMAVBD
      COMMON /MACROA/ LMACVG,LMACMG,LMAVT,LMAMID,LMAADJ,LMAVEL,LMAMEL,
     *                LMANPR,LMAMM,LMAVBD,LMAEBD,LMABCT,LMAVBP,LMAMBP,
     *                LMAVE
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL PARX,PARY,TMAX,S2DI0,S2DB0
      SAVE /ERRCTL/,/TRIAA/,/OUTPUT/,/MACROA/,/MACROD/,/TRIAD/
C
      IF (ICHECK.GE.997) CALL OTRC('XSB1X ','04/12/91')
C
      IF (ISMAS.EQ.0) THEN
C *** use current macro-decomposition
C
       IF (LMAVBD.EQ.0.OR.LMAVBP.EQ.0.OR.NMAVBD.EQ.0) THEN
        CALL WERR(-164,'XOWA  ')
        GOTO 99999
       ENDIF
       CALL XSCL
       IF (IER.NE.0) GOTO 99999
       NEL=NMAEL
       NVT=NMAVT
       NMT=NMAEDG
       NVE=NMAVE
       NVEL=0
       NBCT=NMABCT
       NVBD=NMAVBD
C
       CALL ZCPY(LMACVG,'DMACVG',LCORVG,'DCORVG')
       IF (IER.NE.0) GOTO 99999
       IF (LMACMG.GT.0) THEN
        CALL ZCPY(LMACMG,'DMACMG',LCORMG,'DCORMG')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       CALL ZCPY(LMAVT,'KMAVT ',LVERT,'KVERT ')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LMAADJ,'KMAADJ ',LADJ,'KADJ  ')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LMANPR,'KMANPR',LNPR,'KNPR  ')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LMAMM,'KMAMM ',LMM,'KMM   ')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LMAVBD,'KMAVBD',LVBD,'KVBD  ')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LMAVBP,'DMAVBP',LVBDP,'DVBDP ')
       IF (IER.NE.0) GOTO 99999
C
C
      ELSE IF (ISMAS.GE.1) THEN
C *** accept current subdivision as macro-decomposition
       IF (LVBD.EQ.0.OR.LVBDP.EQ.0.OR.NVBD.EQ.0) THEN
        CALL WERR(-164,'XOWA  ')
        GOTO 99999
       ENDIF
       IF (ISMAS.GT.1) THEN
C *** save current mesh
        CALL XOWS(MFILE,CFILE,IFMT)
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       CALL XSMAS
C
C
      ELSE
C *** read macro mesh from file
C *** clean up possible prior decomposition and macro information
       CALL XSCL
       IF (IER.NE.0) GOTO 99999
       CALL XSMACL
       IF (IER.NE.0) GOTO 99999
       CALL XORS(MFILE,CFILE,IFMT)
       IF (LVBD.EQ.0.OR.LVBDP.EQ.0.OR.NVBD.EQ.0) THEN
        CALL WERR(-164,'XOWA  ')
        GOTO 99999
       ENDIF
C
       CALL XSMAS
C
      ENDIF
C
C
C
      CALL SBD02(KWORK(L(LVBD)),DWORK(L(LVBDP)),DWORK(L(LMBDP)),
     *           DWORK(L(LCORVG)))
C
      IF (NMT.GT.0) THEN
       CALL ZDISP(0,LMBDP,'DMBDP ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      CALL ZDISP(0,LVBDP,'DVBDP ')
      IF (IER.NE.0) GOTO 99999
      CALL ZDISP(0,LVBD,'KVBD  ')
      IF (IER.NE.0) GOTO 99999
      IF (LEBD.NE.0) CALL ZDISP(0,LEBD,'KEBD  ')
      IF (IER.NE.0) GOTO 99999
      NVBD=0
C
      IF (ICHECK.GT.0) THEN
       CALL XS2C(PARX,PARY,TMAX)
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      CALL XSB1(NFINE,IDISP,PARX,PARY,TMAX)
      IF (IER.NE.0) GOTO 99999
C
      IF (IVEL.GT.0) THEN
       CALL XS2V
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      IF (IMID.GT.0) THEN
       CALL XS2M(IMID,IADJ,IDISP,S2DI0,S2DB0,PARX,PARY,TMAX)
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      IF (IADJ.EQ.0.AND.LADJ.GT.0) CALL ZDISP(0,LADJ,'KADJ  ')
C
      IBDP1=IBDP
      IF (IBDP1.EQ.2) IBDP1=1
      CALL XSVEB(IBDP1,TMAX)
C
      IF (IBDP.EQ.2) THEN
       CALL ZNEW(NVBD,-1,LVBDP,'DVBDP ')
       IF (IER.NE.0) GOTO 99999
       IF (NMT.GT.0) THEN
        CALL ZNEW(NVBD,-1,LMBDP,'DMBDP ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       CALL SBD03(DWORK(L(LCORVG)),KWORK(L(LVBD)),
     *            DWORK(L(LVBDP)),DWORK(L(LMBDP)))
      ENDIF
C
      CALL SBD04(DWORK(L(LCORVG)),KWORK(L(LNPR)),PARX,PARY)
C      
99999 END
