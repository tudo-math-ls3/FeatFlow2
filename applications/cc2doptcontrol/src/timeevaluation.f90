!##############################################################################
!# ****************************************************************************
!# <name> timeevaluation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides functions to evaluate a space-time function given
!# by a space-time vector in time.
!#
!# The module contains the following routines:
!#
!# 1.) fetevl_evaluate
!#     -> Evaluate a space-time FE function at a given time. 
!# </purpose>
!##############################################################################

MODULE timeevaluation

  USE fsystem
  USE genoutput
  USE externalstorage
  USE spatialdiscretisation
  USE linearsystemscalar
  USE linearsystemblock
  USE collection
  USE vectorio
  
  USE timediscretisation
  USE spacetimevectors

  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE tmevl_evaluate(rspaceTimeVector,dtime,rvector)
  
!<description>
  ! Evaluates a space-time function at the time dtime. The result is
  ! written to rvector, which is a vector of a function in space.
!</description>

!<input>
  ! Space-time vector to be evaluated
  TYPE(t_spaceTimeVector), INTENT(IN) :: rspaceTimeVector
  
  ! Time where to evaluate
  REAL(DP), INTENT(IN) :: dtime
!</input>

!<inputoutput>
  ! A block vector in space that receives the result.
  ! If the block vector is empty, a new vector is created.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    REAL(DP) :: dabstime,dntimesteps,dreltime
    INTEGER :: itimestep,ntimesteps
    
    ! Rescale the time to the interval 0..ntimeintervals
    ntimesteps = rspaceTimeVector%p_rtimeDiscretisation%nintervals
    dntimesteps = REAL(ntimesteps,DP)
    CALL mprim_linearRescale(dtime,&
      rspaceTimeVector%p_rtimeDiscretisation%dtimeInit,&
      rspaceTimeVector%p_rtimeDiscretisation%dtimeMax,&
      0.0_DP,dntimesteps,dabstime)

    ! Is the destination vector empty? If yes, create a new one.
    IF (rvector%NEQ .EQ. 0) THEN
      CALL lsysbl_createVecBlockByDiscr (rspaceTimeVector%p_rblockDiscretisation,&
          rvector,.FALSE.)
    END IF
    
    ! Clear the destination vector
    CALL lsysbl_clearVector (rvector)
    
    ! Now we have to take a look on the time discretisation.
    ! Depending on that, we choose the best suitable evaluation
    ! method.
    IF (rspaceTimeVector%p_rtimeDiscretisation%ctype .EQ. TDISCR_THETA) THEN
      
      ! Get the time step which is closest to the time stamp.
      itimestep = INT(dabstime + 0.5_DP)
      
      IF (dabstime .EQ. REAL(itimestep,DP)) THEN
        ! Nice coincidence, we have exactly timestep itimestep. Ok, then we 
        ! can call the routine to get that timestep; this saves us some
        ! time as the interpolation can be omitted.
        CALL sptivec_getTimestepData (rspaceTimeVector, 1+itimestep, rvector)
        RETURN
      END IF
      
      ! Get the 'relative' evaluation time; this in the interval -1..1
      dreltime = dabstime-REAL(itimestep,DP)
      
      IF (rspaceTimeVector%NEQtime .EQ. 2) THEN
        ! Special case: only one timestep!
        ! Interpolate linearly.
        CALL interpolateLinear (dreltime,0,1,rspaceTimeVector,rvector)
      ELSE
        ! Is this the first or the last timestep?
        IF (itimestep .EQ. 0) THEN
          ! First timestep. Interpolate between timesteps 0,1 and 2, evaluate 
          ! near timestep 0.
          CALL interpolateQuadratic (dreltime,0,1,2,rspaceTimeVector,rvector)
        ELSE IF (itimestep .EQ. ntimesteps) THEN
          ! Last timestep. Interpolate between timesteps n-2,n-1 and n, evaluate 
          ! near timestep n.
          CALL interpolateQuadratic (dreltime,&
            ntimesteps-2,ntimesteps-1,ntimesteps,rspaceTimeVector,rvector)
        ELSE
          ! Somewhere in the inner. Get the number of the previous and next timestep
          ! and interpolate there.
          CALL interpolateQuadratic (dreltime,&
            itimestep-1,itimestep,itimestep+1,rspaceTimeVector,rvector)
        END IF
      END IF
      
    ELSE IF (rspaceTimeVector%p_rtimeDiscretisation%ctype .EQ. TDISCR_DG0) THEN
    
      ! Get the time step midpoint which is closest to the time stamp.
      itimestep = INT(dabstime)
      
      IF ((dabstime-0.5_DP) .EQ. REAL(itimestep,DP)) THEN
        ! Nice coincidence, we have exactly timestep itimestep. Ok, then we 
        ! can call the routine to get that timestep; this saves us some
        ! time as the interpolation can be omitted.
        CALL sptivec_getTimestepData (rspaceTimeVector, 1+itimestep, rvector)
        RETURN
      END IF
        
      ! Get the 'relative' evaluation time; this in the interval -1..1
      dreltime = dabstime-0.5_DP-REAL(itimestep,DP)
      
      IF (rspaceTimeVector%NEQtime .EQ. 2) THEN
        ! Special case: only one timestep!
        ! Get the one and only solution
        CALL sptivec_getTimestepData (rspaceTimeVector, 1+0, rvector)
        
      ELSE IF (rspaceTimeVector%NEQtime .EQ. 3) THEN
        ! Special case: only two timestep!
        ! Interpolate linearly.
        CALL interpolateLinear (dreltime,0,1,rspaceTimeVector,rvector)
        
      ELSE
        ! Is this the first or the last timestep?
        IF (itimestep .EQ. 0) THEN
          ! First timestep. Interpolate between timesteps 0,1 and 2, evaluate 
          ! near timestep 0.
          CALL interpolateQuadratic (dreltime,0,1,2,rspaceTimeVector,rvector)
        ELSE IF (itimestep .EQ. ntimesteps-1) THEN
          ! Last timestep. Interpolate between timesteps n-2,n-1 and n, evaluate 
          ! near timestep n.
          CALL interpolateQuadratic (dreltime,&
            ntimesteps-1-2,ntimesteps-1-1,ntimesteps-1,rspaceTimeVector,rvector)
        ELSE
          ! Somewhere in the inner. Get the number of the previous and next timestep
          ! and interpolate there.
          CALL interpolateQuadratic (dreltime,&
            itimestep-1,itimestep,itimestep+1,rspaceTimeVector,rvector)
        END IF
      END IF
              
    ELSE
      CALL output_line ('Unsupported time discretisation.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fetevl_evaluate')    
      CALL sys_halt()
    END IF
    
  CONTAINS
  
    SUBROUTINE interpolateLinear (dt,istep1,istep2,rspaceTimeVector,rvector)
    
    ! Calculates the linear interpolation of timestep istep1 and timestep
    ! istep2 of the space time vector rspaceTimeVector. The result is
    ! written to rvector.
    
    ! Interpolation weight; range -1..1 with -1~istep1 and 1~istep2.
    REAL(DP), INTENT(IN) :: dt
    
    ! 'Left' time step corresponding to dt=0.
    INTEGER, INTENT(IN) :: istep1

    ! 'Right' time step corresponding to dt=1.
    INTEGER, INTENT(IN) :: istep2
    
    ! Space time vector containing the data.
    TYPE(t_spaceTimeVector), INTENT(IN) :: rspaceTimeVector
    
    ! Block vector receiving the result
    TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
    
      ! local variables
      TYPE(t_vectorBlock) :: rvec
      
      ! Get a temporary vector
      CALL lsysbl_createVecBlockIndirect (rvector,rvec,.FALSE.)
      
      ! Get the two timesteps
      CALL sptivec_getTimestepData (rspaceTimeVector, 1+istep1, rvec)
      CALL sptivec_getTimestepData (rspaceTimeVector, 1+istep2, rvector)
      
      ! Interpolate
      CALL lsysbl_vectorLinearComb (rvec,rvector,&
          0.5_DP*(dt+1.0_DP),0.5_DP*(1.0_DP-dt))
      
      ! Release unused data
      CALL lsysbl_releaseVector (rvec)
      
    END SUBROUTINE
    
    ! ---------------------------------------------------------------
    
    SUBROUTINE interpolateQuadratic (dt,istep1,istep2,istep3,rspaceTimeVector,rvector)
    
    ! Calculates the linear interpolation of timestep istep1, istep2 and istep3
    ! of the space time vector rspaceTimeVector. The result is
    ! written to rvector.
    
    ! Interpolation weight; range -1..1 with 0~istep1 and 0~istep2 and 1~istep3
    REAL(DP), INTENT(IN) :: dt
    
    ! 'Left' time step corresponding to dt=-1.
    INTEGER, INTENT(IN) :: istep1

    ! 'Middle' time step corresponding to dt=0.
    INTEGER, INTENT(IN) :: istep2

    ! 'Right' time step corresponding to dt=1.
    INTEGER, INTENT(IN) :: istep3
    
    ! Space time vector containing the data.
    TYPE(t_spaceTimeVector), INTENT(IN) :: rspaceTimeVector
    
    ! Block vector receiving the result
    TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
    
      ! local variables
      INTEGER(PREC_DOFIDX) :: i
      TYPE(t_vectorBlock) :: rvec1,rvec2
      REAL(DP), DIMENSION(:), POINTER :: p_Ddata1,p_Ddata2,p_Ddata3
      
      ! Get some temporary vectors
      CALL lsysbl_createVecBlockIndirect (rvector,rvec1,.FALSE.)
      CALL lsysbl_createVecBlockIndirect (rvector,rvec2,.FALSE.)
      
      ! Get the two timesteps
      CALL sptivec_getTimestepData (rspaceTimeVector, 1+istep1, rvec1)
      CALL sptivec_getTimestepData (rspaceTimeVector, 1+istep2, rvec2)
      CALL sptivec_getTimestepData (rspaceTimeVector, 1+istep3, rvector)

      ! Get the data arrays
      CALL lsysbl_getbase_double (rvec1,p_Ddata1)
      CALL lsysbl_getbase_double (rvec2,p_Ddata2)
      CALL lsysbl_getbase_double (rvector,p_Ddata3)
      
      ! Interpolate
      DO i=1,SIZE(p_Ddata3)
        CALL mprim_quadraticInterpolation (dreltime,&
            p_Ddata1(i),p_Ddata2(i),p_Ddata3(i),p_Ddata3(i))
      END DO
      
      ! Release unused data
      CALL lsysbl_releaseVector (rvec2)
      CALL lsysbl_releaseVector (rvec1)
      
    END SUBROUTINE

  END SUBROUTINE

END MODULE
