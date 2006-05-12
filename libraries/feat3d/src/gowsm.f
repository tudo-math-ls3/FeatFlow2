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
* GOWSM                                                                *
*                                                                      *
* Purpose  Write subdivision in MOVIE.BYU format onto unit MFILE       *
*                                                                      *
* Subroutines/functions called  XOWA                                   *
*                                                                      *
* Version from  11/21/90                                               *
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
*                -111  MFILE  exceeds maximum range 16,...,80          *
*                -111  Write error on unit  MFILE                      *
*                -112  File  CFILE  could not be opened                *
*                                                                      *
************************************************************************
C
      SUBROUTINE GOWSM(MFILE,CFILE)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CFILE*15
C
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
      SUB='GOWSM'
      IF (ICHECK.GE.997) CALL OTRC('GOWSM ','11/21/90')
      IER=0
C
C *** Open I/O-file
C
      IFMT=1
      CALL OF0(MFILE,CFILE,IFMT)
      IF (IER.NE.0) GOTO 99999
C
      NP=1
      NJ=NVT
      NPT=NEL*NAE
      NEDGE=NEL*NAE*4
      NPTH=1
C
      WRITE (MFILE,'(4I8)',ERR=99997) NP,NJ,NPT,NEDGE
      WRITE (MFILE,'(2I8)',ERR=99997) NPTH,NPT
C
      WRITE (MFILE,'(1P6E12.5)') (REAL(DWORK(L(LCORVG)+I)),I=0,3*NVT-1)
      IF (IER.NE.0) GOTO 99997
C
      CALL ZNEW(NEDGE,-3,LMOVV,'KMOVV ')
      IF (IER.NE.0) GOTO 99997
      CALL MOVV(KWORK(L(LMOVV)),KWORK(L(LVERT)))
      WRITE (MFILE,'(10I8)') (KWORK(L(LMOVV)+I),I=0,NEDGE-1)
      CALL ZDISP(0,LMOVV,'KMOVV ')
      IF (IER.NE.0) GOTO 99997
      GOTO 99999
C
99997 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-111,'GOWSM ')
C
99999 END
C
C
C
      SUBROUTINE MOVV(KVERMO,KVERT)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KVERMO(*),KVERT(8,*)
      DIMENSION KIAD(4,6)
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      SAVE /TRIAD/
      DATA KIAD/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
C
      DO 1 IEL=1,NEL
      DO 2 IAR=1,NAE
      DO 3 IVE=1,3
3     KVERMO((IEL-1)*NAE*4+(IAR-1)*4+IVE)= KVERT(KIAD(IVE,IAR),IEL)
2     KVERMO((IEL-1)*NAE*4+IAR*4)=        -KVERT(KIAD(  4,IAR),IEL)
1     CONTINUE
C
      END
