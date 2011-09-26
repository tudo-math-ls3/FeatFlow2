!************************************************************************
!* This file contains the COMMON blocks that are used by elements to
!* return calculated values.
!*
!* It's an extension to files CBASICELEM.INC and CBASICTRIA.INC, which
!* have to be included prior to this file!
!************************************************************************

      DOUBLE PRECISION DX(NNVE),DY(NNVE),DJAC(2,2), DBAS(NNBAS,NNDER)
      DOUBLE PRECISION COB,DETJ
      INTEGER KVE(NNVE), IEL
      LOGICAL BDER(NNDER)

      COMMON /ELEM/   DX,DY,DJAC,DETJ,DBAS,BDER,KVE,IEL
      COMMON /COFBAS/ COB(NNBAS,NNCOF)

      SAVE /ELEM/, /COFBAS/