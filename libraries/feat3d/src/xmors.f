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
* XMORS                                                                *
*                                                                      *
* Purpose  Write multiple triangulations (multigrid version)           *
*                                                                      *
* Subroutines/functions called  XORS                                   *
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
      SUBROUTINE XMORS(MFILE,CCFILE,IFMT)
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
      SAVE /TRIAA/,/TRIAD/,/MGTRD/,/MGTRA/,/MGPAR/,/ERRCTL/,/CHAR/
C
      SUB='XMORS'
      IF (ICHECK.GE.997) CALL OTRC('XMORS ','08/25/90')
      IER=0
C
      CALL XMSCL
      IF (IER.NE.0) GOTO 99999
      CALL XSCL
      IF (IER.NE.0) GOTO 99999
C      
      DO 10 ILEV=1,NLEV
C
      NEL =0
      NVT =0
      NET =0
      NAT =0
      NVE =0
      NEE =0
      NAE =0
      NVEL =0
      NEEL =0
      NVED =0
      NVAR =0
      NEAR =0
      NBCT =0
      NVBD =0
      NEBD =0
      NABD =0
C
      LCORVG=0
      LCORMG=0
      LCORAG=0
      LVERT =0
      LEDGE =0
      LAREA =0
      LADJ =0
      LVEL =0
      LEEL =0
      LAEL =0
      LVED =0
      LAED =0
      LVAR =0
      LEAR =0
      LEVE =0
      LAVE =0
      LNPR =0
      LBCT =0
      LVBD =0
      LEBD =0
      LABD =0
C
      CALL XORS(MFILE,CCFILE(ILEV),IFMT)
      IF (IER.NE.0) GOTO 99999
      CLOSE(MFILE)
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
      KLCVG(ILEV) =LCORVG
      KLCMG(ILEV) =LCORMG
      KLCAG(ILEV) =LCORAG
      KLVERT(ILEV)=LVERT
      KLEDGE(ILEV) =LEDGE
      KLAREA(ILEV) =LAREA
      KLADJ(ILEV) =LADJ
      KLVEL(ILEV) =LVEL
      KLEEL(ILEV) =LEEL
      KLAEL(ILEV) =LAEL  
      KLVED(ILEV) =LVED
      KLAED(ILEV) =LAED 
      KLVAR(ILEV) =LVAR
      KLEAR(ILEV) =LEAR 
      KLEVE(ILEV) =LEVE
      KLAVE(ILEV) =LAVE
      KLNPR(ILEV) =LNPR
      KLBCT(ILEV) =LBCT
      KLVBD(ILEV) =LVBD
      KLEBD(ILEV) =LEBD
      KLABD(ILEV) =LABD 
C
10    CONTINUE   
C     
99999 END
      
      
