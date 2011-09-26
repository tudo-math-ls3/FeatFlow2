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
* XSB0                                                                 *
*                                                                      *
* Purpose  Adjust dimensions of DCORVG, KVERT, KEDGE, KAREA, KADJ      *
*          and KNPR                                                    *
*                                                                      *
* Subroutines/functions called  SB0, SBVEL, SBE, SBEEL, SBA, SBAEL,    *
*                               ZLEN, ZDISP, ZNEW, ZCPY                *
*                                                                      *
* Version from  12/11/89                                               *
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
* SB0                                                                  *
*                                                                      *
* Purpose  Regular refinement of a given hexahedral subdivision        *
*          of a three-dimensional domain                               *
*                                                                      *
* Subroutines/functions called   SBV, SBVEL, SBE                       *
*                                                                      *
* Version from  12/11/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DCORVG   R*8    Cartesian coordinates of interior vertices           *
*                 Parameter of vertices on the boundary                *
* KVERT    I*4    Numbers of vertices in each element, counterclockwise*
* KEDGE    I*4    Numbers of edges in each element, counterclockwise   *
* KAREA    I*4    Numbers of areas in each element, counterclockwise   *
* KADJ     I*4    Number of neighbouring element, 0 for boundary edges *
* KNPR     I*4    0    for interior vertices                           *
*                 IBCT for nodal points on boundary component IBCT     *
* NNEL     I*4    Maximum dimension provided for KVERT                 *
* NNVT     I*4    Maximum dimension provided for DCORVG, KNPR          *
*                                                                      *
* NFINE    I*4    Number of regular subdivisions                       *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DCORVG   R*8    As above                                             *
* KVERT    I*4    As above                                             *
* KEDGE    I*4    As above                                             *
* KAREA    I*4    As above                                             *
* KADJ     I*4    As above                                             *
* KNPR     I*4    As above                                             *
* KMM      I*4    As above                                             *
*                                                                      *
* NNEL     I*4    Dimension needed for KVERT, KADJ                     *
* NNVT     I*4    Dimension needed for DCORVG, KNPR                    *
*                                                                      *
* NEL      I*4    Number of elements of final triangulation            *
* NVT      I*4    Number if vertices of final triangulation            *
* NET      I*4    Number if edges of final triangulation               *
* NAT      I*4    Number if areas of final triangulation               *
*                                                                      *
* IER      I*4     0  No error                                         *
*                  1  NNEL or NNVT too small                           *
*                     Requirements                                     *
*                     NNEL >=NEL(NFINE)                                *
*                     NNVT >=NVT(NFINE)                                *
*                                                                      *
************************************************************************
C
      SUBROUTINE XSB0(NFINE,ISE,ISA,ISVEL,ISEEL,ISAEL,ISVED,ISAED,
     *                ISVAR,ISEAR,ISEVE,ISAVE,ISVBD,ISEBD,ISABD,IDISP,
     *                PARX,PARY,PARZ,SEDB,SADB)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNVE=8,NNAE=6,NNEE=12,NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAH/  LCORH
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL PARX,PARY,PARZ,SEDB,SADB
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/,/TRIAH/
C
      SUB='XSB0'
      IF (ICHECK.GE.997) CALL OTRC('XSB0  ','12/01/93')
C
      CALL ZLEN(LVERT,ILEN)
      NNEL=ILEN/NNVE
      CALL ZLEN(LADJ,ILEN)
      NNEL1=ILEN/NNAE
      CALL ZLEN(LEDGE,ILEN)
      NNEL2=ILEN/NNEE
      NNEL=MIN0(NNEL,NNEL1,NNEL2)
      CALL ZLEN(LCORVG,ILEN)
      NNVT=ILEN/3
      NNEL1=NNEL
      NNVT1=NNVT
C
      IF (LCORH.NE.0) THEN
       CALL ZLEN(LCORVG,ILEN1)
       CALL ZLEN(LCORH ,ILEN2)
       IF (ILEN1.NE.ILEN2) GOTO 99999
       CALL ZDISP(0,LCORVG,'DCVORG')
       IF (IER.NE.0) GOTO 99999
       LCORVG=LCORH
       LCORH =0
      ENDIF
C
      NVEL=0
      CALL ZLEN(LVEL,ILEN)
      IF (ILEN.GT.0) CALL ZDISP(0,LVEL,'KVEL ')
      CALL ZNEW(NVEL,-3,LVEL,'KVEL  ')
      IF (IER.NE.0) GOTO 99999
C
      CALL SB0(DWORK(L(LCORVG)),DWORK(L(LCORMG)),KWORK(L(LVERT)),
     *         KWORK(L(LADJ)),KWORK(L(LEDGE)),KWORK(L(LNPR)),
     *         KWORK(L(LVEL)),NFINE,NNEL,NNVT,PARX,PARY,PARZ,
     *         SEDB,SADB)
C
      IF (IER) 99999,1,2
C
1     IF ((ISVEL.GT.0).OR.(ISE.GT.0)) THEN
       CALL SBVEL(KWORK(L(LVERT)),KWORK(L(LVEL)),0)
       CALL SBVEL(KWORK(L(LVERT)),KWORK(L(LVEL)),1)
       CALL ZDISP(NVEL*NVT,LVEL,'KVEL  ')
       IF (IER.NE.0) GOTO 99999
      ELSE
       CALL ZDISP(0,LVEL,'KVEL  ')
       IF (IER.NE.0) GOTO 99999
       LVEL=0
      ENDIF
C
      CALL ZNEW(3*NVT,-1,LCORH,'DCORH ')
      IF (IER.NE.0) GOTO 99999
      CALL ZCPY (LCORVG,'DCORVG',LCORH,'DCORH ')
      IF (IER.NE.0) GOTO 99999
      CALL SVC(DWORK(L(LCORVG)),KWORK(L(LNPR)),NVT,PARX,PARY,PARZ)
      IF (IER.NE.0) GOTO 99999
C
      GOTO 99995
C
C
2     CALL ZDISP(0,LVEL,'KVEL  ')
C
      IF (NNEL.GT.NNEL1) THEN
       CALL ZNEW(NNVE*NNEL,3,LV1,'KVERT ')
       IF (IER.NE.0) GOTO 99999
       IF (LVERT.NE.0) THEN
        CALL ZCPY (LVERT,'KVOLD ',LV1,'KVERT ')
        CALL ZDISP(0,LVERT,'KVOLD ')
       ENDIF
       LVERT=LV1
C
       CALL ZNEW(NNAE*NNEL,3,LA1,'KADJ  ')
       IF (IER.NE.0) GOTO 99999
       IF (LADJ.NE.0) THEN
        CALL ZCPY (LADJ,'KAOLD ',LA1,'KADJ  ')
        CALL ZDISP(0,LADJ,'KAOLD ')
       ENDIF
       LADJ=LA1
C
       CALL ZNEW(NNEE*NNEL,3,LE1,'KEDGE ')
       IF (IER.NE.0) GOTO 99999
       IF (LEDGE.NE.0) THEN
        CALL ZCPY (LEDGE,'KEOLD ',LE1,'KEDGE ')
        CALL ZDISP(0,LEDGE,'KEOLD ')
       ENDIF
       LEDGE=LE1
      ENDIF
C
      IF (NNVT.GT.NNVT1) THEN
       CALL ZNEW(3*NNVT,1,LCV1,'DCORVG')
       IF (IER.NE.0) GOTO 99999
       IF (LCORVG.NE.0) THEN
        CALL ZCPY (LCORVG,'DCVOLD',LCV1,'DCORVG')
        CALL ZDISP(0,LCORVG,'DCVOLD')
       ENDIF
       LCORVG=LCV1
C
       CALL ZNEW(NNVT,3,LNPR1,'KNPR  ')
       IF (IER.NE.0) GOTO 99999
       IF (LNPR.NE.0) THEN
        CALL ZCPY (LNPR,'KNPOLD',LNPR1,'KNPR  ')
        CALL ZDISP(0,LNPR,'KNPOLD')
       ENDIF
       LNPR=LNPR1
       LNPR1=0
      ENDIF
C
      NVEL=0
      CALL ZLEN(LVEL,ILEN)
      IF (ILEN.GT.0) CALL ZDISP(0,LVEL,'KVEL ')
      CALL ZNEW(NVEL,-3,LVEL,'KVEL  ')
      IF (IER.NE.0) GOTO 99999
C
      CALL SB0(DWORK(L(LCORVG)),DWORK(L(LCORMG)),KWORK(L(LVERT)),
     *         KWORK(L(LADJ)),KWORK(L(LEDGE)),KWORK(L(LNPR)),
     *         KWORK(L(LVEL)),NFINE,NNEL,NNVT,PARX,PARY,PARZ,
     *         SEDB,SADB)
C
C
      IF ((ISVEL.GT.0).OR.(ISE.GT.0)) THEN
       CALL SBVEL(KWORK(L(LVERT)),KWORK(L(LVEL)),0)
       CALL SBVEL(KWORK(L(LVERT)),KWORK(L(LVEL)),1)
       CALL ZDISP(NVEL*NVT,LVEL,'KVEL  ')
       IF (IER.NE.0) GOTO 99999
      ELSE
       CALL ZDISP(0,LVEL,'KVEL  ')
       IF (IER.NE.0) GOTO 99999
       LVEL=0
      ENDIF
C
      CALL ZNEW(3*NVT,-1,LCORH,'DCORH ')
      IF (IER.NE.0) GOTO 99999
      CALL ZCPY (LCORVG,'DCORVG',LCORH,'DCORH ')
      IF (IER.NE.0) GOTO 99999
      CALL SVC(DWORK(L(LCORVG)),KWORK(L(LNPR)),NVT,PARX,PARY,PARZ)
      IF (IER.NE.0) GOTO 99999
C
C
99995 IF (ISE.GT.0) THEN
       IF (ISE.GE.2) THEN
        LNPR1=0
        CALL ZNEW(NVT+8*NET,3,LNPR1,'KNPR1 ')
        IF (IER.NE.0) GOTO 99999
        CALL ZCPY (LNPR,'KNPOLD',LNPR1,'KNPR  ')
        IF (IER.NE.0) GOTO 99999
        CALL ZDISP(0,LNPR,'KNPOLD')
        IF (IER.NE.0) GOTO 99999
        LNPR=LNPR1
        LNPR1=0
       ENDIF
C        
       IF (ISE.GE.3) THEN
        IF (LCORMG.NE.0) THEN
         CALL ZDISP(0,LCORMG,'DCORMG')
         IF (IER.NE.0) GOTO 99999
        ENDIF
        CALL ZNEW(3*8*NET,1,LCORMG,'DCOREG')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C        
       CALL SBE(KWORK(L(LVERT)),KWORK(L(LNPR)),KWORK(L(LVEL)),
     *          DWORK(L(LCORH)),KWORK(L(LEDGE)),DWORK(L(LCORMG)),ISE,
     *          PARX,PARY,PARZ,SEDB)
       IF (IER.NE.0) GOTO 99999
C
       IF (ISE.GE.3) THEN
        CALL ZDISP(3*NET,LCORMG,'DCOREG')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C        
       IF (ISE.GE.2) THEN
        CALL ZDISP(NVT+NET,LNPR,'KNPR  ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C        
      ENDIF
C
C
      IF ((ISVEL.EQ.0).AND.(ISE.GT.0)) THEN
       CALL ZDISP(0,LVEL,'KVEL  ')
       IF (IER.NE.0) GOTO 99999
       LVEL=0
      ENDIF
C
C
      IF ((ISEEL.GT.0).AND.(ISE.GT.0)) THEN
       NEEL=0
       CALL ZLEN(LEEL,ILEN)
       IF (ILEN.GT.0) CALL ZDISP(0,LEEL,'KEEL ')
       CALL ZNEW(NEEL,-3,LEEL,'KEEL  ')
       IF (IER.NE.0) GOTO 99999
       CALL SBEEL(KWORK(L(LEDGE)),KWORK(L(LEEL)),0)
       CALL SBEEL(KWORK(L(LEDGE)),KWORK(L(LEEL)),1)
       CALL ZDISP(NEEL*NET,LEEL,'KEEL  ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
C
      IF (ISA.GT.0) THEN
       CALL ZLEN(LAREA,ILEN)
       IF (ILEN.GT.0) CALL ZDISP(0,LAREA,'KAREA')
       CALL ZNEW(NNAE*NEL,3,LAREA,'KAREA ')
       IF (IER.NE.0) GOTO 99999
C
       IF (ISA.GE.2) THEN
        LNPR1=0
        IF (ISE.GE.2) THEN
         INDA=NVT+NET
        ELSE
         INDA=NVT
        ENDIF
        CALL ZNEW(INDA+8*NAT,3,LNPR1,'KNPR1 ')
        IF (IER.NE.0) GOTO 99999
        CALL ZCPY (LNPR,'KNPOLD',LNPR1,'KNPR  ')
        IF (IER.NE.0) GOTO 99999
        CALL ZDISP(0,LNPR,'KNPOLD')
        IF (IER.NE.0) GOTO 99999
        LNPR=LNPR1
        LNPR1=0
       ENDIF
C        
       IF (ISA.GE.3) THEN
        IF (LCORAG.NE.0) THEN
         CALL ZDISP(0,LCORAG,'DCORAG')
         IF (IER.NE.0) GOTO 99999
        ENDIF
        CALL ZNEW(3*8*NAT,1,LCORAG,'DCORAG')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C        
       CALL SBA(KWORK(L(LVERT)),KWORK(L(LNPR)),KWORK(L(LADJ)),
     *          DWORK(L(LCORH)),KWORK(L(LAREA)),DWORK(L(LCORAG)),
     *          ISA,INDA,PARX,PARY,PARZ,SADB)
       IF (IER.NE.0) GOTO 99999
C
       IF (ISA.GE.3) THEN
        CALL ZDISP(3*NAT,LCORAG,'DCORAG')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C        
       IF (ISA.GE.2) THEN
        CALL ZDISP(INDA+NAT,LNPR,'KNPR  ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C        
      ENDIF
C
C
      IF ((ISAEL.GT.0).AND.(ISA.GT.0)) THEN
       CALL ZLEN(LAEL,ILEN)
       IF (ILEN.GT.0) CALL ZDISP(0,LAEL,'KAEL ')
       CALL ZNEW(2*NAT,3,LAEL,'KAEL  ')
       IF (IER.NE.0) GOTO 99999
       CALL SBAEL(KWORK(L(LAREA)),KWORK(L(LAEL)))
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
C
      IF (LCORH.NE.0) THEN
       CALL ZDISP(0,LCORH,'DCVOH ')
       IF (IER.NE.0) GOTO 99999
       LCORH =0
      ENDIF
C
C
C
99999 END
C
C
C
      SUBROUTINE SB0(DCORVG,DCOREG,KVERT,KADJ,KEDGE,KNPR,KVEL,NFINE,
     *               NNEL,NNVT,PARX,PARY,PARZ,SEDB,SADB)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=8,NNAE=6,NNEE=12)
      DIMENSION DCORVG(3,*),DCOREG(3,*)
      DIMENSION KVERT(NNVE,*),KADJ(NNAE,*),KEDGE(NNEE,*),KNPR(*)
      DIMENSION KVEL(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      EXTERNAL PARX,PARY,PARZ,SEDB,SADB
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='SB0'
      IF (ICHECK.GE.997) CALL OTRC('SB0   ','12/11/89')
C
      NNEL1=NNEL
      NNVT1=NNVT
      NNEL=NEL
      NNVT=NVT
      NNET=NET
      NNAT=NAT
C
      IF (NFINE.GT.0) THEN
       DO 1 IFINE=1,NFINE
       NNVT=NNVT+NNET+NNAT+NNEL
       NNET=2*NNET+6*NNEL+4*NNAT
       NNAT=4*NNAT+12*NNEL
1      NNEL=8*NNEL
      ENDIF
C
      IF (NNEL.GT.NNEL1.OR.NNVT.GT.NNVT1) THEN
       IER=1
       GOTO 99999
      ENDIF
C
      IER=0
C
      DO 100 IFINE=1,NFINE-1
      CALL SBV(DCORVG,KVERT,KEDGE,KADJ,KNPR,PARX,PARY,PARZ,SEDB,SADB)
      CALL SBVEL(KVERT,KVEL,0)
      CALL SBVEL(KVERT,KVEL,1)
      CALL SBE(KVERT,KNPR,KVEL,DCORVG,KEDGE,DCOREG,1,PARX,PARY,PARZ,
     *         SEDB)
100   CONTINUE
C
      IF (NFINE.GT.0) CALL SBV(DCORVG,KVERT,KEDGE,KADJ,KNPR,
     *                         PARX,PARY,PARZ,SEDB,SADB)
C
C
C
99999 END
