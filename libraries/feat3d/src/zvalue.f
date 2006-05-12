      BLOCK DATA ZVALUE
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299,NNLEV=9)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
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
      COMMON /TABLE/  KTYPE(NNARR),KLEN(NNARR),KLEN8(NNARR),IFLAG
C
      DATA M/2/,MT/2/,MKEYB/5/,MTERM/6/,IER/0/,ICHECK/1/
      DATA MERR/11/,MPROT/12/,MSYS/13/,MTRC/14/,IRECL8/512/
      DATA SUB/'MAIN  '/
      DATA FMT/'(3D24.16)','(5E14.7)','(6I12)'/
      DATA CPARAM/'                                                   '/
      DATA NEL/0/,NVT/0/,NET/0/,NAT/0/,NVE/0/,NEE/0/,NAE/0/,NVEL/0/,
     *     NEEL/0/,NVED/0/,NVAR/0/,NEAR/0/,NBCT/0/,NVBD/0/,NEBD/0/,
     *     NABD/0/
      DATA LCORVG/0/,LCORMG/0/,LCORAG/0/,LVERT/0/,LEDGE/0/,LAREA/0/,
     *     LADJ/0/,LVEL/0/,LEEL/0/,LAEL/0/,LVED/0/,LAED/0/,LVAR/0/,
     *     LEAR/0/,LEVE/0/,LAVE/0/,LNPR/0/,LBCT/0/,LVBD/0/,LEBD/0/,
     *     LABD/0/
      DATA KNEL/NNLEV*0/,KNVT/NNLEV*0/,KNET/NNLEV*0/,
     *     KNAT/NNLEV*0/,KNVE/NNLEV*0/,KNEE/NNLEV*0/,    
     *     KNAE/NNLEV*0/,KNVEL/NNLEV*0/,KNEEL/NNLEV*0/, 
     *     KNVED/NNLEV*0/,KNVAR/NNLEV*0/,KNEAR/NNLEV*0/,
     *     KNBCT/NNLEV*0/,KNVBD/NNLEV*0/,KNEBD/NNLEV*0/, 
     *     KNABD/NNLEV*0/          
      DATA KLCVG/NNLEV*0/,KLCMG/NNLEV*0/,KLCAG/NNLEV*0/,
     *     KLVERT/NNLEV*0/,KLEDGE/NNLEV*0/,KLAREA/NNLEV*0/,
     *     KLADJ/NNLEV*0/,KLVEL/NNLEV*0/,KLEEL/NNLEV*0/,
     *     KLAEL/NNLEV*0/,KLVED/NNLEV*0/,KLAED/NNLEV*0/,
     *     KLVAR/NNLEV*0/,KLEAR/NNLEV*0/,KLEVE/NNLEV*0/,
     *     KLAVE/NNLEV*0/,KLNPR/NNLEV*0/,KLBCT/NNLEV*0/,
     *     KLVBD/NNLEV*0/,KLEBD/NNLEV*0/,KLABD/NNLEV*0/
      DATA KTYPE/NNARR*0/,KLEN/NNARR*0/,KLEN8/NNARR*0/,IFLAG/0/
C
      END
