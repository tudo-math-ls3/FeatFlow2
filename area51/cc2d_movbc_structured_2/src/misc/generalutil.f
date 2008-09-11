************************************************************************
* General collection of routines which might help somewhere...
************************************************************************

************************************************************************
* Linear rescale
*
* Scales a coordinate x linearly from the interval [a,b] to the
* interval [c,d].
*
* In:
*   x     - coordinate to be rescaled
*   [a,b] - source interval
*   [c,d] - destination interval
*
* Out:
*   y     - Rescaled coordinate
************************************************************************

      SUBROUTINE LRSCLE (X,A,B,C,D,Y)
      
      IMPLICIT NONE
      
      DOUBLE PRECISION X,A,B,C,D,Y
      
      DOUBLE PRECISION D1,D2,D3
      
C     Calculate the coefficients of the transformation 
C         D1*A+D2 = C, D1*B+D2 = D.
C     Use them to calculate Y=D1*X+D2.

      IF (A.EQ.B) THEN
        Y = C
        RETURN
      END IF

      D3 = 1D0/(A-B)
      D1 = (C-D)*D3
      D2 = (-B*C+A*D)*D3
      
      Y = D1*X+D2
      
      END
