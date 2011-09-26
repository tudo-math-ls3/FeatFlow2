************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
* based on the 2-D routine, modificated by P.Schreiber                 *
* XMOWS                                                                *
*                                                                      *
* Purpose  Write multiple triangulations (multigrid version)           *
*                                                                      *
* Subroutines/functions called  XOWS                                   *
*                                                                      *
* Version from  08/25/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NLEV     I*4    Vector of numbers of the arrays                      *
*                 New arrays are allocated on DWORK for LA(IBLOC)=0    *
* MFILE    I*4    Unit number used for output                          *
* CCFILE   C*(*)  Array of filenames used for output                   *
* Meshes on COMMON /MGTRD/ and /MGTRA/                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Error indicator                                      *
*                 Set by OF0 or OWA                                    *
*                                                                      *
************************************************************************
C
      SUBROUTINE XMOWS(MFILE,CCFILE,IFMT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER*(*) CCFILE(*)
C
      PARAMETER (NNLEV=9)
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDG(NNLEV),KLARE(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /TRIAA/,/TRIAD/,/MGTRD/,/MGTRA/,/MGPAR/,/ERRCTL/,/CHAR/
C
      SUB='XMOWS'
      IF (ICHECK.GE.997) CALL OTRC('XMOWS ','08/25/90')
      IER=0
C
      DO 10 ILEV=1,NLEV
C
      NEL =KNEL(ILEV)    
      NVT =KNVT(ILEV)
      NET =KNET(ILEV)
      NAT =KNAT(ILEV)
      NVE =KNVE(ILEV)
      NEE =KNEE(ILEV)
      NAE =KNAE(ILEV)
      NVEL =KNVEL(ILEV)
      NEEL =KNEEL(ILEV)
      NVED =KNVED(ILEV)
      NVAR =KNVAR(ILEV)
      NEAR =KNEAR(ILEV)
      NBCT =KNBCT(ILEV)
      NVBD =KNVBD(ILEV)
      NEBD =KNEBD(ILEV)
      NABD =KNABD(ILEV)
C
      LCORVG =KLCVG(ILEV) 
      LCORMG =KLCMG(ILEV) 
      LCORAG =KLCAG(ILEV) 
      LVERT  =KLVERT(ILEV)
      LEDGE  =KLEDG(ILEV) 
      LAREA  =KLARE(ILEV) 
      LADJ   =KLADJ(ILEV) 
      LVEL   =KLVEL(ILEV) 
      LEEL   =KLEEL(ILEV)
      LAEL   =KLAEL(ILEV)   
      LVED   =KLVED(ILEV) 
      LAED   =KLAED(ILEV) 
      LVAR   =KLVAR(ILEV) 
      LEAR   =KLEAR(ILEV) 
      LEVE   =KLEVE(ILEV) 
      LAVE   =KLAVE(ILEV) 
      LNPR   =KLNPR(ILEV) 
      LBCT   =KLBCT(ILEV) 
      LVBD   =KLVBD(ILEV) 
      LEBD   =KLEBD(ILEV) 
      LABD   =KLABD(ILEV) 
C
      CALL XOWS(MFILE,CCFILE(ILEV),IFMT)
      IF (IER.NE.0) GOTO 99999
      CLOSE(MFILE)
C     
10    CONTINUE   
C     
99999 END
      
      
