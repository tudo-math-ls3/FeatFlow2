!************************************************************************
!* This file contains COMMON blocks responsible for the output to the
!* screen / to a file.
!************************************************************************

      INTEGER M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8

      CHARACTER SUB*6,FMT*15,CPARAM

      COMMON /CHAR/   SUB,FMT(3),CPARAM(120)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8

      SAVE /CHAR/, /OUTPUT/
      