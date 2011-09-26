************************************************************************
* This file contains COMMON blocks that store information about the
* parametrization. It's used by PARQ2D.F.
************************************************************************

      INTEGER NBCT,NNCOMP,NNNPAR,LNCOMP,LICPTR,LITYP,LNSPLN,LNPAR
      INTEGER LIPPTR,LXPAR,LYPAR

      COMMON /TDATA/ NBCT,NNCOMP,NNNPAR,
     *               LNCOMP,LICPTR,LITYP,LNSPLN,LNPAR,
     *               LIPPTR,LXPAR,LYPAR

      SAVE /TDATA/