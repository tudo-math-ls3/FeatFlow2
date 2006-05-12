************************************************************************
* This file contains routines for the maintainance of the time
* stepping scheme. It provides routines to fill the time stepping
* structures with values for Euler, Crank-Nicolson and Fractional
* step.
* The underlying structure is TTimeStepDParams in the file
* STSTEPPING.INC.
************************************************************************

************************************************************************
* Initialize time stepping
*
* This routine initializes the time stepping scheme. Based on the
* parameters about the type of the simulation, the time stepping
* structure will be filled with data about the weights in that time
* step.
*
* In:
*   TIMENS - Current simulation time of the time step.
*            Standard = 0D0.
*   IFRSTP - Initialize for fractional step.
*            = 0: Don't use fractional step, use standard One-step
*                 scheme
*            = 1: Initialize for Fractional step
*   THETA  - time stepping identifier.
*            =0: Forward Euler
*            =1: Backward Euler
*            =0.5: Crank Nicolson
*   TSTEP  - (Theoretical) Time step length of the simulation.
*            Might vary from the real time step size by the type of the
*            scheme. The actual time step size will be stored into
*            TSCHEM by this routine.
*   ISSTP  - Number of substep, the weights should be initialized for.
*            Ignored for stationary simulation.
*            Depending on the time-stepping technique used, the caller
*            has to use this parameter to select the correct weights.
*            Which values are allowed for ISSTP is depending on the
*            time-stepping technique:
*             IFRSTP=0: ISSTP=1..1 - One-Step scheme; there doesn't 
*                       exist any substeps.
*             IFRSTP=1: ISSTP=1..3 - Fractional Step scheme; the
*                       weights are 3-periodic.
*
* Out:
*   TSCHEM - array [1..SZTSTD,1..NSSTPS] of double
*            TTimeStepDParams structure, filled with data for every
*            substep.
************************************************************************

      SUBROUTINE INITSC (TSCHEM,TIMENS,IFRSTP,THETA,TSTEP,
     *                   ISSTP)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'ststepping.inc'
      
      INTEGER ISSTP
      INTEGER IFRSTP
      DOUBLE PRECISION THETA,TSTEP,TIMENS,TSCHEM(SZTSTD)
      
      DOUBLE PRECISION THETA1, THETP1, ALPHA, BETA
      
C     Clear the structure

      CALL LCL1 (TSCHEM,SZTSTD)
      
C     Do we have a stationary simulation?
C
C     --> definitely not, so we comment the if-clause out.
C         Stationary simulations are handled by a separate solver
C         and not as a special case of the nonstationary solver!

C      IF (INONST.EQ.0) THEN
      
C       No, that's a stationary one. Initialize the time stepping
C       in a simple manner for performing one step.
C
C       For stationary simulations we perform only "1 step" - some
C       kind of pseudo time step with time step size 1D0.
      
C        THETA1 = 1D0
C      
C        TSCHEM(OTHETA) = THETA1
C      
C        TSCHEM(OTHSTPS) = 1D0
C        TSCHEM(OTSTEP)  = 1D0
C        
C       Initialize weights for matrices and RHS:

C        TSCHEM(OTMWGH)  = TSCHEM(OTSTEP) *  THETA1
C        TSCHEM(OTRMWGH) = TSCHEM(OTSTEP) * (THETA1-1D0)
C        TSCHEM(OTR1WGH) = TSCHEM(OTSTEP) *  THETA1
C        TSCHEM(OTR2WGH) = -TSCHEM(OTR1WGH)
C        TSCHEM(OTRSWGH) = TSCHEM(OTSTEP)

C      ELSE 
      
C       Ok, the nonstationary case is more complex.

        TSCHEM(OCURTIM) = TIMENS
      
C       Do we have standard time stepping or Fractional Step?
      
        IF (IFRSTP.EQ.0) THEN

C         Standard time stepping. Here the parameters are a little bit
C         easier to initialize. Loop through the substeps to initialize:

          TSCHEM(OTHETA)  = THETA
          TSCHEM(OTHSTEP) = TSTEP * THETA

C         Initialize the substep length:

          TSCHEM(OTSTEP)  = TSTEP
      
C         Initialize weights for matrices and RHS:

          TSCHEM(OTMWGH)  = TSTEP *  THETA
          TSCHEM(OTRMWGH) = TSTEP * (THETA-1D0)
          TSCHEM(OTR1WGH) = TSTEP *  THETA
          TSCHEM(OTR2WGH) = -TSCHEM(OTRMWGH)
          TSCHEM(OTRSWGH) = TSTEP
            
        ELSE IF (IFRSTP.EQ.1) THEN
        
C         Fractional Step. Now it gets interesting...
C         We can't use the data provided by the caller, since that has
C         no influence to the Fractional Step scheme - Fractional step
C         uses predefined values!

          IF ((ISSTP.LT.1).OR.(ISSTP.GT.3)) THEN
            WRITE (MTERM,'(A)') 'Error: Irregular value for substep!'
            STOP
          END IF

C         The FS Theta-Scheme uses by theory 4 parameters: 
C
C           Theta   = 1 - sqrt(2) / 2
C           Theta'  = 1 - 2 * Theta
C           alpha   = ( 1 - 2 * Theta ) / ( 1 - Theta )
C           beta    = 1 - alpha
C
C         The parameter THETA in the DAT-file is ignored and replaced
C         by a hard-coded setting.
        
          THETA1 = 1D0-SQRT(0.5D0)
          THETP1 = 1D0-2D0*THETA1
          ALPHA  = THETP1 / (1D0-THETA1)
          BETA   = THETA1 / (1D0-THETA1)
      
          TSCHEM(OTHETA)  = THETA1
          TSCHEM(OTHETP)  = THETP1
          TSCHEM(OALPHA)  = ALPHA
          TSCHEM(OBETA)   = BETA

C         For fractional step the handling of the time step size TSTEP
C         is slightly different.
C         There we are orienting on the length of the macrostep of
C         step length 3xTSTEP and break up that into three different
C         substeps at different points in time, not corresponding to
C         TSTEP. Depending on the number of the current substep, we
C         have two settings for the weights and time length:

          IF (ISSTP.NE.2) THEN
          
C           1st and 3rd substep:

            TSCHEM(OTSTEP)  =  3D0*TSTEP * THETA1

            TSCHEM(OTMWGH)  =  3D0*TSTEP * ALPHA * THETA1
            TSCHEM(OTRMWGH) = -3D0*TSTEP * BETA  * THETA1
            TSCHEM(OTR1WGH) =  3D0*TSTEP * ALPHA * THETA1
            TSCHEM(OTR2WGH) =  3D0*TSTEP * BETA  * THETA1
            TSCHEM(OTRSWGH) =  3D0*TSTEP * THETA1

          ELSE
          
C           2nd substep

            TSCHEM(OTSTEP)  =  3D0*TSTEP * THETP1

            TSCHEM(OTMWGH)  =  3D0*TSTEP * ALPHA * THETA1
            TSCHEM(OTRMWGH) = -3D0*TSTEP * ALPHA * THETP1
            TSCHEM(OTR1WGH) =  3D0*TSTEP * BETA  * THETP1
            TSCHEM(OTR2WGH) =  3D0*TSTEP * ALPHA * THETP1
            TSCHEM(OTRSWGH) =  3D0*TSTEP * THETP1
            
          END IF ! (ISSTP.NE.2)
          
        ELSE

          WRITE (MTERM,'(A)') 'Error: Unknown time stepping technique!'
          STOP

        END IF ! (IFRSTP.EQ.1)
      
C      END IF ! (INONST.EQ.0)
      
      END
      