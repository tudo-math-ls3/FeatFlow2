************************************************************************
* This file contains the COMMON blocks describing the (depending on
* the situation) current triangulation. In multigrid routines
* the routine SETLEV initializes these variables depending on the
* desired level.
************************************************************************

      INTEGER NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      INTEGER LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *        LVBD,LEBD,LBCT,LVBDP,LMBDP

      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP

      SAVE /TRIAD/, /TRIAA/
      