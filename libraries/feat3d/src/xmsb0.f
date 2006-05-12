************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT3D                              *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
* based on the 2-D routine, modificated by P.Schreiber                 *
* XMSB0                                                                *
*                                                                      *
* Purpose  Generate sequence of meshes for multigrid applications      *
*          by successive calls of XSB0X  (two-level ordering)          *
*                                                                      *
* Subroutines/functions called  XSB0X, WERR, ZCPY                      *
*                                                                      *
* Version from  03/02/93                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NLEV     I*4    desired number of l                                  *
* ISCAD    I*4    =1 Determine array KADJ from coarse grid             *
* ISE      I*4    =1 Determine number of midpoints                     *
* ISA      I*4    =0 Release array KADJ on return                      *
*                    after determination of the new subdivisions       *
* ISEEL    I*4    =1 Determine numbers of elements meeting at each     *
*                    edge                                              *
* ISAEL    I*4    =1 Determine numbers of elements meeting at each     *
*                    area                                              *
* ISVEL    I*4    =1 Determine numbers of elements meeting at each     *
*                    vertex                                            *
* IDISP    I*4    =1 Release free space on all arrays after completion *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DCORVG   R*8   Cartesian coordinates of vertices                     *                                                                      *                                                                      *
************************************************************************
C
      SUBROUTINE XMSB0(ISCAD,ISE,ISA,ISVEL,ISEEL,ISAEL,
     *                 ISVED,ISAED,ISVAR,ISEAR,ISEVE,ISAVE,
     *                 ISVBD,ISEBD,ISABD,IDISP,PARX,PARY,PARZ,
     *                 SEDB,SADB)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299,NNLEV=9)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EXTERNAL SEDB,SADB,PARX,PARY,PARZ
      SAVE /ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/,/MGTRD/,/MGTRA/,/MGPAR/
C
      SUB='XMSB0 '
      IF (ICHECK.GE.997) CALL OTRC('XMSB0 ','04/12/91')
C
      IF (NLEV.GT.NNLEV) THEN
       CALL WERR(-180,'XMSB0 ')
       GOTO 99999
      ENDIF
C
C
      DO 10 ILEV=1,NLEV
      IF (ILEV.EQ.1) THEN
      CALL XSB0X(0,1,MAX(1,ISE),ISA,ISVEL,ISEEL,ISAEL,
     *           ISVED,ISAED,ISVAR,ISEAR,ISEVE,ISAVE,
     *           ISVBD,ISEBD,ISABD,IDISP,PARX,PARY,PARZ,
     *           SEDB,SADB)
      ELSE
      CALL XSB0X(1,1,MAX(1,ISE),ISA,ISVEL,ISEEL,ISAEL,
     *           ISVED,ISAED,ISVAR,ISEAR,ISEVE,ISAVE,
     *           ISVBD,ISEBD,ISABD,IDISP,PARX,PARY,PARZ,
     *           SEDB,SADB)      
      ENDIF
C
C *** Save dimensions for all levels
C
      KNEL(ILEV) =NEL
      KNVT(ILEV) =NVT
      KNET(ILEV) =NET
      KNAT(ILEV) =NAT
      KNVE(ILEV) =NVE
      KNEE(ILEV) =NEE
      KNAE(ILEV) =NAE
      KNVEL(ILEV) =NVEL
      KNEEL(ILEV) =NEEL
      KNVED(ILEV) =NVED
      KNVAR(ILEV) =NVAR
      KNEAR(ILEV) =NEAR
      KNBCT(ILEV) =NBCT
      KNVBD(ILEV) =NVBD
      KNEBD(ILEV) =NEBD
      KNABD(ILEV) =NABD 
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
       IF (ISA.NE.0) THEN
        CALL ZCPY(LCORAG,'DCORAG',KLCAG(ILEV),'DCAG0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       CALL ZCPY(LVERT,'KVERT ',KLVERT(ILEV),'KVERT0')
       IF (IER.NE.0) GOTO 99999
       IF (LEDGE.EQ.1) THEN
        CALL ZCPY(LEDGE,'KEDGE  ',KLEDGE(ILEV),'KEDG0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (ISA.NE.0) THEN
        CALL ZCPY(LAREA,'KAREA  ',KLAREA(ILEV),'KARE0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (ISA.NE.0) THEN
        CALL ZCPY(LADJ,'KADJ  ',KLADJ(ILEV),'KADJ0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (LVEL.NE.0) THEN
        CALL ZCPY(LVEL,'KVEL  ',KLVEL(ILEV),'KVEL0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (LEEL.NE.0) THEN
        CALL ZCPY(LEEL,'KEEL  ',KLEEL(ILEV),'KEEL0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (LAEL.NE.0) THEN
        CALL ZCPY(LAEL,'KAEL  ',KLAEL(ILEV),'KAEL0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (LVED.NE.0) THEN
        CALL ZCPY(LVED,'KVED  ',KLVED(ILEV),'KVED0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (LAED.NE.0) THEN
        CALL ZCPY(LAED,'KAED  ',KLAED(ILEV),'KAED0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (LVAR.NE.0) THEN
        CALL ZCPY(LVAR,'KVAR  ',KLVAR(ILEV),'KVAR0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (LEAR.NE.0) THEN
        CALL ZCPY(LEAR,'KEAR  ',KLEAR(ILEV),'KEAR0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (LEVE.NE.0) THEN
        CALL ZCPY(LEVE,'KEVE  ',KLEVE(ILEV),'KEVE0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (LAVE.NE.0) THEN
        CALL ZCPY(LAVE,'KAVE  ',KLAVE(ILEV),'KAVE0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       CALL ZCPY(LNPR,'KNPR  ',KLNPR(ILEV),'KNPR0 ')
       IF (IER.NE.0) GOTO 99999
      IF (LBCT.NE.0) THEN
       CALL ZCPY(LBCT,'KLBCT ',KLBCT(ILEV),'KBCT0 ')
       IF (IER.NE.0) GOTO 99999
       ENDIF
      IF (LVBD.NE.0) THEN
       CALL ZCPY(LVBD,'KVBD  ',KLVBD(ILEV),'KVBD0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      IF (LEBD.GE.1) THEN
       CALL ZCPY(LEBD,'KEBD  ',KLEBD(ILEV),'KEBD0 ')
       IF (IER.NE.0) GOTO 99999
       ENDIF
      IF (LABD.GE.1) THEN
       CALL ZCPY(LABD,'KABD  ',KLABD(ILEV),'KABD0 ')
       IF (IER.NE.0) GOTO 99999
      ENDIF   
C
      ELSE
C
       KLCVG(ILEV) =LCORVG
       KLCMG(ILEV) =LCORMG
       KLCAG(ILEV) =LCORAG
       KLVERT(ILEV)=LVERT
       KLEDGE(ILEV)=LEDGE
       KLAREA(ILEV)=LAREA
       KLADJ(ILEV)=LADJ
       KLVEL(ILEV)=LVEL
       KLEEL(ILEV)=LEEL
       KLAEL(ILEV)=LAEL
       KLVED(ILEV)=LVED
       KLAED(ILEV)=LAED
       KLVAR(ILEV)=LVAR
       KLEAR(ILEV)=LEAR
       KLEVE(ILEV)=LEVE
       KLAVE(ILEV)=LAVE
       KLNPR(ILEV)=LNPR
       KLBCT(ILEV)=LBCT
       KLVBD(ILEV)=LVBD
       KLEBD(ILEV)=LEBD
       KLABD(ILEV)=LABD
C     
      ENDIF
C
10    CONTINUE
C
99999 END
