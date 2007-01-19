************************************************************************
*   this subroutine generates random values between 0.0 and 1.0 using
*   an integer seed
*   it is based on the imsl routine ggubs.
************************************************************************

      DOUBLE PRECISION FUNCTION USRAN(IR)
      
      IMPLICIT NONE
      
      INTEGER IR
      
      DOUBLE PRECISION DA,DB,DC
      
      PARAMETER(DA=16807.D0,DB=2147483647.D0,DC=2147483648.D0)
      
      IR=ABS(MOD(DA*IR,DB)+0.5D0)
      USRAN=DFLOAT(IR)/DC
      
      END