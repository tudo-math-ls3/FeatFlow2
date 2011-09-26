************************************************************************
* This file contains common block with variables that are used during
* numerical integration.
************************************************************************

      INTEGER NNCUBP
      PARAMETER (NNCUBP=36)
      
      DOUBLE PRECISION DXI, DOMEGA
      INTEGER NCUBP, ICUBP
      
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP

      SAVE /CUB/
      