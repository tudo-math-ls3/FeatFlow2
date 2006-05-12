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
* GORSM                                                                *
*                                                                      *
* Purpose  Read subdivision in MOVIE.BYU format from unit MFILE        *
*          Data are assumed to be written by GOWSM                     *
*                                                                      *
* Subroutines/functions called  OF0, XORA                              *
*                                                                      *
* Version from  12/11/89                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----   ----                                                         *
* MFILE   I*4    Unit number of external file                          *
* CFILE   C*15   Filename                                              *
*                                                                      *
* OUTPUT  TYPE                                                         *
* ------  ----                                                         *
* MFILE   I*4    Unit number used                                      *
* CFILE   C*15   Filename used                                         *
* IER     I*4    Error indicator                                       *
*                -110  MFILE  exceeds maximum range 16,...,80          *
*                -110  Read error on unit MFILE                        *
*                -112  File CFILE could not be opened                  *
* The remaining output quantities are contained in                     *
* COMMON blocks /TRIAA/ and /TRIAD/                                    *
*                                                                      *
************************************************************************
C
      SUBROUTINE GORSM(MFILE,CFILE)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CFILE*15
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/
C
      SUB='GORSM'
      IF (ICHECK.GE.997) CALL OTRC('GORSM ','12/11/89')
      IER=0
C
C *** Open I/O-file
C
      IFMT=1
      CALL OF0(MFILE,CFILE,IFMT)
      IF (IER.NE.0) GOTO 99999
C
      READ (MFILE,'(4I8)',ERR=99997,END=99997) NP,NVT,NEL,NEDGE
      READ (MFILE,'(2I8)',ERR=99997,END=99997) NPTH,NPT
C
      NVE=8
      NEE=12
      NAE=6
      NBCT=1
C
      CALL ZNEW(3*NVT,-1,LCORVG,'DCORVG')
      IF (IER.NE.0) GOTO 99997
      CALL ZNEW(3*NVT,-2,LCORMO,'VCORMO')
      IF (IER.NE.0) GOTO 99997
      READ (MFILE,'(1P6E12.5)') (VWORK(L(LCORMO)+I),I=0,3*NVT-1)
      DO 1 I=1,3*NVT
1     DWORK(L(LCORVG)+I-1)=DBLE(VWORK(L(LCORMO)+I-1))
      CALL ZDISP(0,LCORMO,'VCORMO')
      IF (IER.NE.0) GOTO 99997
C
      CALL ZNEW(NEDGE,-3,LVERT,'KVERT ')
      IF (IER.NE.0) GOTO 99997
      READ (MFILE,'(10I8)') (KWORK(L(LVERT)+I),I=0,NEDGE-1)
      DO 2 I=1,NEDGE
2     KWORK(L(LVERT)+I-1)=IABS(KWORK(L(LVERT)+I-1))
C
      GOTO 99999
C
99997 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'GORSM ')
C
99999 END
C
C
C
      SUBROUTINE VMOV(KVERT,KVERMO)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KVERMO(*),KVERT(8,*)
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      SAVE /TRIAD/
C
      DO 1 IEL=1,NEL
      KVERT(1,IEL)=IABS(KVERMO((IEL-1)*NAE*4+1))
      KVERT(2,IEL)=IABS(KVERMO((IEL-1)*NAE*4+2))
      KVERT(3,IEL)=IABS(KVERMO((IEL-1)*NAE*4+3))
      KVERT(4,IEL)=IABS(KVERMO((IEL-1)*NAE*4+4))
      KVERT(5,IEL)=IABS(KVERMO((IEL-1)*NAE*4+8))
      KVERT(6,IEL)=IABS(KVERMO((IEL-1)*NAE*4+7))
      KVERT(7,IEL)=IABS(KVERMO((IEL-1)*NAE*4+11))
      KVERT(8,IEL)=IABS(KVERMO((IEL-1)*NAE*4+15))
1     CONTINUE
C
      END
