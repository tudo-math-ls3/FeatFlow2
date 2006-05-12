************************************************************************
* This file contains the parametrisation routines for the
* computational domain:
*   PARX - Obtain the X-coordinate of a point from its parameter value
*   PARX - Obtain the Y-coordinate of a point from its parameter value
*   TMAX - Obtain the maximum parameter value on a boundary component
*   TNBC - Obtain number of real boundary components
*
* The routines here are basically wrapper routines choosing the
* "real" parametrisation routines. The governing variable is the
* variable IMESH of the DAT file:
* If IMESH=1, these routines will use the standard OMEGA implementation 
*             routines for defining the computational domain by an 
*             external file.
* If IMESH=0, these routines will use the user defined FEAT
*             implementation FPARX, FPARY, FTMAX that is provided
*             in the file PARQ2D.F.
* Compatibility remark: If importing an "old" PARQ2D.F into the
* new implementation, the "old" PARX/PARY/TMAX/TNBC-routines have simply
* to be renamed to FPARX/FPARY/FTMAX/FTNBC in order to work.
*
* In contrast to the original parq2d.f from earlier FEATFLOW versions,
* an additional routine TNBC was added that returns the number of
* boundary components NBCT. This makes the parametrization independent
* of any information in the triangulation structures, as NBC was
* preciously stored in the /TRADx/ COMMON blocks!
************************************************************************

************************************************************************
* Get X-coordinate of parameter value.
*
* This routine returns the X-coordinate of a parameter value.
*
* In:
*   T      - parameter value
*   IBCT   - number of boundary component
*
* Return:
*   X-coordinate of parameter value T on boundary component IBCT.
************************************************************************

      DOUBLE PRECISION FUNCTION PARX(T,IBCT)

      IMPLICIT NONE
      
      INCLUDE 'cparametrization.inc'
      
C parameters

      DOUBLE PRECISION T
      INTEGER IBCT
      
C externals

      DOUBLE PRECISION FPARX,OPARX
      EXTERNAL FPARX,OPARX
      
C call the "correct" PARX-implementation

      IF (IMESH.EQ.0) THEN
        PARX = FPARX(T,IBCT)
      ELSE IF (IMESH.EQ.1) THEN
        PARX = OPARX (T,IBCT)
      END IF
      
      END


************************************************************************
* Get Y-coordinate of parameter value.
*
* This routine returns the Y-coordinate of a parameter value.
*
* In:
*   T      - parameter value
*   IBCT   - number of boundary component
*
* Return:
*   Y-coordinate of parameter value T on boundary component IBCT.
************************************************************************

      DOUBLE PRECISION FUNCTION PARY(T,IBCT)

      IMPLICIT NONE
      
      INCLUDE 'cparametrization.inc'

C parameters

      DOUBLE PRECISION T
      INTEGER IBCT
      
C externals

      DOUBLE PRECISION FPARY,OPARY
      EXTERNAL FPARY,OPARY

C call the "correct" PARY-implementation

      IF (IMESH.EQ.0) THEN
        PARY = FPARY(T,IBCT)
      ELSE IF (IMESH.EQ.1) THEN
        PARY = OPARY (T,IBCT)
      END IF

      END


************************************************************************
* Get maximum parameter value
*
* This routine returns the maximum parameter value on a boundary
* component
*
* In:
*   IBCT   - number of boundary component
*
* Return:
*   Maximum parameter value on that boundary component.
************************************************************************

      DOUBLE PRECISION FUNCTION TMAX(IBCT)

      IMPLICIT NONE
      
      INCLUDE 'cparametrization.inc'

C parameters

      INTEGER IBCT
      
C externals

      DOUBLE PRECISION FTMAX,OTMAX
      EXTERNAL FTMAX,OTMAX

C call the "correct" TMAX-implementation

      IF (IMESH.EQ.0) THEN
        TMAX = FTMAX(IBCT)
      ELSE IF (IMESH.EQ.1) THEN
        TMAX = OTMAX (IBCT)
      END IF

      END

************************************************************************
* Get number of real boundary components
*
* This routine returns the number of real boundary components in the
* parametrization.
*
* In:
*   -
*
* Return:
*   Number of boundary components
************************************************************************

      INTEGER FUNCTION TNBC ()

      IMPLICIT NONE
      
      INCLUDE 'cparametrization.inc'

C parameters

      INTEGER IBCT
      
C externals

      INTEGER FTNBC,OTNBC
      EXTERNAL FTNBC,OTNBC

C call the "correct" TNBC-implementation

      IF (IMESH.EQ.0) THEN
        TNBC = FTNBC()
      ELSE IF (IMESH.EQ.1) THEN
        TNBC = OTNBC()
      END IF

      END
