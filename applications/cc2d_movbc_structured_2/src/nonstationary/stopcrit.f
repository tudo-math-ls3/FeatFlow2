************************************************************************
* Calculate stopping criterion for adaptive time stepping
*
* Based on the configuration, CRITAD will compute a bound EPSAD
* that can be used in the adaptive time stepping as stopping criterion.
*
* In:
*   TIMEIN  - Length of startup phase
*   TIMENS  - Current simulation time
*   TIMEST  - Initial simulation time
*   EPSADI  - error limit during startup phase
*   EPSADL  - standard error limit after startup phase
*   IADIN   - Identifier for the error control in the startup phase
*             =0: use EPSADI during startup phase
*             =1: use linear blending from EPSADI to EPSADL during
*                 startup phase
*             =2: use logarithmic blending from EPSADI to EPSADL during
*                 startup phase
*
* Out:
*   EPSAD   - Stopping criterion for adaptive time step control.
*             During the startup phase (time T in the range 
*             TIMEST..TIMEST+TIMEIN), the stopping criterion
*             will be calculated according to IADIN. After the
*             startup phase, EPSADL will be used.
************************************************************************

      SUBROUTINE CRITAD(TIMEIN,TIMENS,TIMEST,EPSADI,EPSADL,EPSAD,IADIN)

      IMPLICIT NONE
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

C parameters

      DOUBLE PRECISION TIMEIN,TIMENS,TIMEST,EPSADI,EPSADL,EPSAD
      INTEGER IADIN

C local variables

      DOUBLE PRECISION TDIFF
             
C     Standard result: EPSADL
             
      EPSAD=EPSADL

C     Do we have a startup phase at all? Would be better, otherwise
C     Navier-Stokes might stagnate with decreasing time steps at the
C     beginning...

      IF (TIMEIN.GT.0) THEN
      
C       Calculate the "current position" of the simulation time in
C       the interval TIMEST..TIMEST+TIMEIN
      
        TDIFF=TIMENS-TIMEST

C       Standard error control in the startup phase

        IF (IADIN.EQ.0) THEN
        
C         Use EPSADI during startup phase and EPSADL afterwards:
        
          IF (TDIFF.LE.TIMEIN) THEN
            EPSAD=EPSADI
          ELSE
            EPSAD=EPSADL
          ENDIF
          
        ENDIF

C       Linear blending as control in the startup phase.
C       Blend linearly between EPSADI and EPSADL from the initial
C       simulation time TIMNEST to the end of the startup phase
C       TIMEST+TIMEIN.
C       After the startup phase, use EPSADL.

        IF (IADIN.EQ.1) THEN
          IF (TDIFF.LE.TIMEIN) THEN
            EPSAD=EPSADI+TDIFF/TIMEIN*(EPSADL-EPSADI)
          ELSE
            EPSAD=EPSADL
          ENDIF
        ENDIF

C       Logarithmic blending as control in the startup phase.
C       Blend logarithmically between EPSADI and EPSADL from the initial
C       simulation time TIMNEST to the end of the startup phase
C       TIMEST+TIMEIN.
C       After the startup phase, use EPSADL.

        IF (IADIN.EQ.2) THEN
          IF (TDIFF.LE.TIMEIN) THEN
            EPSAD=EPSADI**(1D0-TDIFF/TIMEIN)*EPSADL**(TDIFF/TIMEIN)
          ELSE
            EPSAD=EPSADL
          ENDIF
        ENDIF

      ENDIF

      END
