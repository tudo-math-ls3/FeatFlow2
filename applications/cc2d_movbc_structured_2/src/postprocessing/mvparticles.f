************************************************************************
* This file contains auxiliary routines for the movement of 
* fictitious boundary particles by the flow.
************************************************************************

************************************************************************
* Rotate a particle
*
* Thus subroutine rotates a particle.
* For this purpose, an ODE is solved by a simple explicit Euler
* approach.
*
* In:
*   TORQUE - Torque force acting on the particle. Can be calculated
*            e.g. by BDFVLT.
*   INERT  - Inertia of the particle. Particle dependent. This is
*            a constant assigned to each particle
*   ROT    - Current rotation angle of the particle
*   TSTEP  - Time step size
*
* Out:
*   ROTNEW - New rotation angle of the particle.
************************************************************************

      SUBROUTINE MVPROT (TORQUE,INERT,ROT,TSTEP,ROTNEW)
      
      IMPLICIT NONE
      
      DOUBLE PRECISION TORQUE,INERT,ROT,TSTEP,ROTNEW
      
C     To rotate the particle, we solve a somple ODE:
C
C       T = I * dw/dt, w(0)=ROT
C
C     with w=new rotation angle. We do this by a simple explicit Euler.

      IF (INERT.NE.0D0) THEN
        ROTNEW = ROT + TSTEP*TORQUE/INERT
      ELSE
        ROTNEW = ROT
      END IF
      
      END
      