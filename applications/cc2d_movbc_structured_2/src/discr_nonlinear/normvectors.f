************************************************************************
* Transform Q0 pressure vector to L2_0
*
* This routine transforms the vector P into the space L2_0.
* For this purpose, the vector AREA with the areas of all elements
* is used. AREA(NEL+1) is assumed to be the sum of all AREA(IEL).
*
* In:
*   NEL    - Number of elements in the geometry
*   P      - array [1..NEL] of double
*            Pressure vector in the space P0, one pressure value
*            per element
*   AREA   - array [1..NEL+1] of double
*            Area of all elements; AREA[NEL+1] is the sum of the
*            areas of all elements
*
* Out:
*   P      - the modified vector
************************************************************************

      SUBROUTINE TOL20A (P,AREA,NEL)

      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'

C     parameters

      INTEGER INEUM, NEL
      DOUBLE PRECISION AREA(*)
      DOUBLE PRECISION P(*)
      
C     local variables

      DOUBLE PRECISION PINT,C
      INTEGER IEL

C     Build the integral
C       int_Omega p dx
C     This is approximated by
C       PINT = SUM_Elements P(Element)*Volume(Element)

      PINT=0D0
      DO IEL=1,NEL
        PINT=PINT+P(IEL)*AREA(IEL)
      END DO

C     Divide PINT by the volume of the domain; this gives the integral
C     mean value of the pressure:

      C = PINT/AREA(NEL+1)

C     Subtract the integral mean value C of the pressure from all
C     pressure components. Afterwards, P has integral mean value = 0.

      DO IEL=1,NEL
        P(IEL)=P(IEL)-C
      END DO

      END
