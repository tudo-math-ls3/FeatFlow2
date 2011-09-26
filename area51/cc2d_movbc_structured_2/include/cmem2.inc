************************************************************************
* This file is an extension to the file CMEM.FOR. It contains
* additional COMMON blocks that are internally usen in the
* pseudodynamical memory management.
*
* As this file is an extension to CMEM.FOR, the CMEM.FOR file must
* be included prior to this file!
************************************************************************

      INTEGER KTYPE, KLEN, KLEN8, IFLAG
      
      COMMON /TABLE/  KTYPE(NNARR),KLEN(NNARR),KLEN8(NNARR),IFLAG
      
      SAVE /TABLE/
