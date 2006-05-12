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
* XMSB1                                                                *
*                                                                      *
* Purpose  Generate sequence of meshes for multigrid applications      *
*          by successive calls of XSB1X  (block ordering)              *
*                                                                      *
* Subroutines/functions called  XSB0X, WERR, ZCPY                      *
*                                                                      *
* Version from  04/12/91                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NFINE0   I*4    degree of refinemenet of basic mesh                  *
* NLEV     I*4    desired number of l                                  *
* IMID     I*4                                                         *
* IADJ     I*4                                                         *
* IVEL     I*4                                                         *
* IDISP    I*4                                                         *
* IBDP     I*4     See parameter list of XSB0X                         *
* PARX     FNCT                                                        *
* PARY     FNCT                                                        *
* TMAX     FNCT                                                        *
* Coarse grid on /TRIAD and /TRIAA/                                    *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DCORVG   R*8   Cartesian coordinates of vertices                     *                                                                      *                                                                      *
************************************************************************
C
      SUBROUTINE XMSB1(NFINE0,IMID,IADJ,IVEL,IDISP,IBDP,
     *                 MFILE,CFILE,IFMT,PARX,PARY,TMAX)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299,NNLEV=9)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER*(*) CFILE
C
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EXTERNAL PARX,PARY,TMAX
      SAVE /ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/,/MGTRD/,/MGTRA/,/MGPAR/
C
      SUB='XMSB1 '
      IF (ICHECK.GE.997) CALL OTRC('XMSB1 ','04/12/91')
C
      IF (NLEV.GT.NNLEV) THEN
       CALL WERR(-180,'XMSB1 ')
       GOTO 99999
      ENDIF
C
C
      IADJ0=1
      IBDP0=2
      NFINE1=NFINE0
C
      CALL OF0(MFILE,CFILE,IFMT)
      IF (IER.NE.0) GOTO 99999
C
      DO 10 ILEV=1,NLEV
      IF (ILEV.EQ.NLEV) THEN
       IADJ0=IADJ
       IBDP0=IBDP
      ENDIF
      REWIND (MFILE)
      CALL XSB1X(NFINE1,IMID,IADJ0,IVEL,IDISP,IBDP0,
     *           -1,MFILE,CFILE,IFMT,PARX,PARY,TMAX)
      NFINE1=2*NFINE1
C
C
C *** Save dimensions for all levels
C
      KNEL(ILEV) =NEL
      KNVT(ILEV) =NVT
      KNMT(ILEV) =NMT
      KNVEL(ILEV)=NVEL
      KNVBD(ILEV)=NVBD     
C
C
C *** Save arrays for all levels
C
C *** DCORVG need not necessarily be saved because of two-level ordering
C
      IF (ILEV.LT.NLEV) THEN
C
       CALL ZCPY(LCORVG,'DCORVG',KLCVG(ILEV),'DCVG0 ')
       IF (IER.NE.0) GOTO 99999
       IF (LCORMG.NE.0) THEN
        CALL ZCPY(LCORMG,'DCORMG',KLCMG(ILEV),'DCMG0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       CALL ZCPY(LVERT,'KVERT ',KLVERT(ILEV),'KVERT0')
       IF (IER.NE.0) GOTO 99999
       IF (LMID.NE.0) THEN
        CALL ZCPY(LMID,'KMID  ',KLMID(ILEV),'KMID0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (IADJ.EQ.1) THEN
        CALL ZCPY(LADJ,'KADJ  ',KLADJ(ILEV),'KADJ0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (LVEL.NE.0) THEN
        CALL ZCPY(LVEL,'KVEL  ',KLVEL(ILEV),'KVEL0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (LMEL.NE.0) THEN
        CALL ZCPY(LMEL,'KMEL  ',KLMEL(ILEV),'KMEL0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       CALL ZCPY(LNPR,'KNPR  ',KLNPR(ILEV),'KNPR0 ')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LMM,'KMM   ',KLMM(ILEV),'KMM0  ')
       IF (IER.NE.0) GOTO 99999
       IF (IBDP.GE.0) THEN
        CALL ZCPY(LVBD,'KVBD  ',KLVBD(ILEV),'KVBD0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (IBDP.GE.1) THEN
       CALL ZCPY(LEBD,'KEBD  ',KLEBD(ILEV),'KEBD0 ')
       IF (IER.NE.0) GOTO 99999
       ENDIF
       CALL ZCPY(LBCT,'KLBCT ',KLBCT(ILEV),'KLBCT0')
       IF (IER.NE.0) GOTO 99999
       IF (IBDP.GE.2) THEN
        CALL ZCPY(LVBDP,'DVBDP ',KLVBDP(ILEV),'DVBDP0')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (IBDP.GE.2.AND.LMBDP.NE.0) THEN
        CALL ZCPY(LMBDP,'DMBDP ',KLMBDP(ILEV),'DMBDP0')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
      ELSE
C
       KLCVG(ILEV) =LCORVG
       KLCMG(ILEV) =LCORMG
       KLVERT(ILEV)=LVERT
       KLMID(ILEV) =LMID
       KLADJ(ILEV) =LADJ
       KLVEL(ILEV) =LVEL
       KLMEL(ILEV) =LMEL
       KLNPR(ILEV) =LNPR
       KLMM(ILEV)  =LMM
       KLVBD(ILEV) =LVBD
       KLEBD(ILEV) =LEBD
       KLBCT(ILEV) =LBCT
       KLVBDP(ILEV)=LVBDP
       KLMBDP(ILEV)=LMBDP
C
      ENDIF
C
10    CONTINUE
C
CC      DO 20 ILEV=1,NLEV
CC20    KLCVG(ILEV)=LCORVG
C
99999 END
