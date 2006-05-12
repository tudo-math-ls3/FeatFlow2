***********************************************************************
* This file provides filtering routines for pure Neumann problems.
***********************************************************************

***********************************************************************
* Calculate the integral mean value of the given vector
***********************************************************************

      DOUBLE PRECISION FUNCTION DCLIMV (DX,NEQ)
      
      IMPLICIT NONE
      
C parameters
      
      INTEGER NEQ
      DOUBLE PRECISION DX(NEQ)
      
C local variables

      DOUBLE PRECISION SM
      INTEGER I
      
C Sum up the coefficients, divide by NEQ to get the mean value.
C For piecewise constant functions this is exactly the integral
C mean value (because of midpoint rule). For non-constant
C functions this is an approximation to the mean value of
C 2nd order.

      SM=0D0
      DO I=1,NEQ
        SM=SM+DX(I)
      END DO
      
      SM=SM/DBLE(NEQ)

      DCLIMV = SM

      END

***********************************************************************
* The following routine changes the given coefficient vector to have
* integral mean value=0.
*
* In:
*  DX     - Solution vector; array [1..NEQ] of double
*  NEQ    - length of solution vector
*
* Out:
*  DX     - updated solution vector.
*
* As this function interpretes the entries in the solution vector
* as function values in nodes, this technique only ensures
* integral-mean-value=0 with conforming finite element spaces!
***********************************************************************

      SUBROUTINE IMVZER (DX,NEQ)
      
      IMPLICIT NONE
      
C parameters
      
      INTEGER NEQ
      DOUBLE PRECISION DX(NEQ)
      
C local variables

      DOUBLE PRECISION SM
      INTEGER I
      
      DOUBLE PRECISION DCLIMV
      EXTERNAL DCLIMV
      
C Sum up the coefficients, divide by NEQ to get the mean value.
C For piecewise constant functions this is exactly the integral
C mean value (because of midpoint rule). For non-constant
C functions this is an approximation to the mean value of
C 2nd order.

      SM=DCLIMV(DX,NEQ)
      
C Subtract this from the solution vector to ensure integral mean
C value = 0 (at least approximatively...)

      DO I=1,NEQ
        DX(I)=DX(I)-SM
      END DO

C      SM=DCLIMV(DX,NEQ)
C Debug: Recalculate the integral mean value      
C      SM=0D0
C      DO I=1,NEQ
C        SM=SM+DX(NEQ)
C      END DO
C      SM=SM/DBLE(NEQ)

      END
      
***********************************************************************
* Neumann filtering routine for multigrid.
* 
* This routine is called in multigrid on different positions of the
* algorithm to perform filtering.
* This version performs filtering to integral-mean-value=0
*
* In:
*  DX     - Solution vector; array [1..NEQ] of double
*  NEQ    - length of solution vector
*  IALGP  - Position in the algorithm
*           0=undefined - enforce the filtering
*           1=on start of the algorithm on finest level; 
*             DX is solution vector
*           2=on start of the MG sweep on the current level;
*             DX is the first solution vector before the
*             sweep starts (0D0-array normally, escept for the
*             finest level, there it's the result of the previous
*             call with IALGP=1) 
*           3=before presmoothing; DX is solution vector 
*           4=after presmoothing; DX is solution vector 
*           5=before restriction; DX is defect vector
*           6=after restriction; DX is defect vector
*           7=before coarse grid solver;
*             DX is RHS vector on coarsest level
*             (the result of IALGP=6 !)
*           8=before coarse grid solver;
*             DX is start vector on coarsest level (normally filled
*             with 0)
*           9=after coarse grid solver; 
*             DX is calculated solution vector on coarsest level
*          10=before prolongation;
*             DX is update-vector on coarser level
*          11=after prolongation; DX is prolongated update vector
*          12=after coarse grid correction, before post smoothing;
*             DX is solution vector
*          13=after post-smoothing; DX is solution vector
*          14=after one step of MG; DX is the solution vector on the
*             finest level.
*  IPARAM - array [1..SZMGRI] of integer
*           Integer parameter structure of the solver
*  DPARAM - array [1..SZMGRI] of integer
*           Double precision parameter structure of the solver
*  IDATA  - array [1..*] of integer
*           User defined integer array
*  DDATA  - array [1..*] of double
*           User defined double array
*
* Out:
*  DX     - updated vector.
***********************************************************************
      
      SUBROUTINE YNMFLT (DX,NEQ,IALGP,IPARAM,DPARAM,IDATA,DDATA)
      
      IMPLICIT NONE

C parameters
      
      INTEGER NEQ,IALGP
      DOUBLE PRECISION DX(NEQ)
      
      DOUBLE PRECISION DCLIMV
      EXTERNAL DCLIMV
      
      INTEGER IPARAM(*),IDATA(*)
      DOUBLE PRECISION DPARAM(*),DDATA(*)

C Perform the filtering after the coarse grid correction and after the 
C post-smoothing. This is the point where the new solution is generated,
C which is only defined up to an additive constant. Here we ensure
C integral-mean-value=0.

C Caution: The solution must ensure Du/Dn=0, so it's not possible to
C prescribe an arbitrary RHS! But: A RHS to a solution with Du/Dn=0
C automatically provides integral mean value=0...

      IF ((IALGP.EQ.0).OR.(IALGP.EQ.13)) THEN
        CALL IMVZER (DX,NEQ)
      END IF

C      PRINT *,'MG-Pos: ',IALGP,', IMV=',DCLIMV(DX,NEQ)
      
99999 END

***********************************************************************
* Neumann filtering routine for multigrid M020.
* 
* This routine is called in multigrid on different positions of the
* algorithm to perform filtering.
* This version performs filtering to integral-mean-value=0
*
* In:
*  DX     - Solution vector; array [1..NEQ] of double
*  NEQ    - length of solution vector
*  IALGP  - Position in the algorithm
*           0=undefined - enforce the filtering
*           1=on start of the algorithm on finest level; 
*             DX is solution vector
*           2=on start of the MG sweep on the current level;
*             DX is the first solution vector before the
*             sweep starts (0D0-array normally, escept for the
*             finest level, there it's the result of the previous
*             call with IALGP=1) 
*           3=before presmoothing; DX is solution vector 
*           4=after presmoothing; DX is solution vector 
*           5=before restriction; DX is defect vector
*           6=after restriction; DX is defect vector
*           7=before coarse grid solver;
*             DX is RHS vector on coarsest level
*             (the result of IALGP=6 !)
*           8=before coarse grid solver;
*             DX is start vector on coarsest level (normally filled
*             with 0)
*           9=after coarse grid solver; 
*             DX is calculated solution vector on coarsest level
*          10=before prolongation;
*             DX is update-vector on coarser level
*          11=after prolongation; DX is prolongated update vector
*          12=after coarse grid correction, before post smoothing;
*             DX is solution vector
*          13=after post-smoothing; DX is solution vector
*          14=after one step of MG; DX is the solution vector on the
*             finest level.
*   IPARAM- array [1..SZMGRI] of integer
*           Integer parameter structure of multigrid solver
*   DPARAM- array [1..SZMGRI] of integer
*           Double precision parameter structure of multigrid solver
*   IDATA - array [1..*] of integer
*           User defined integer array
*   DDATA - array [1..*] of double
*           User defined double array
*
* Out:
*  DX     - updated vector.
***********************************************************************
      
      SUBROUTINE YNMFL2 (DX,NEQ,IALGP,IPARAM,DPARAM,IDATA,DDATA)
      
      IMPLICIT NONE

C parameters
      
      INTEGER NEQ,IALGP
      DOUBLE PRECISION DX(NEQ)

      INTEGER IPARAM(*)
      DOUBLE PRECISION DPARAM(*)
      INTEGER IDATA(*)
      DOUBLE PRECISION DDATA(*)
      
      DOUBLE PRECISION DCLIMV
      EXTERNAL DCLIMV

C Perform the filtering after the coarse grid correction and after the 
C post-smoothing. This is the point where the new solution is generated,
C which is only defined up to an additive constant. Here we ensure
C integral-mean-value=0.

C Caution: The solution must ensure Du/Dn=0, so it's not possible to
C prescribe an arbitrary RHS! But: A RHS to a solution with Du/Dn=0
C automatically provides integral mean value=0...

      IF ((IALGP.EQ.0).OR.(IALGP.EQ.13)) THEN
        CALL IMVZER (DX,NEQ)
      END IF

C      PRINT *,'MG-Pos: ',IALGP,', IMV=',DCLIMV(DX,NEQ)
      
99999 END
