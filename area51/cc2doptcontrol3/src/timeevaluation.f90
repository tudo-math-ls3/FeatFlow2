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

module timeevaluation

  use fsystem
  use genoutput
  use externalstorage
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use collection
  use vectorio
  
  use timediscretisation
  use spacetimevectors
  use mprimitives

  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine tmevl_evaluate(rspaceTimeVector,dtime,rvector)
  
!<description>
  ! Evaluates a space-time function at the time dtime. The result is
  ! written to rvector, which is a vector of a function in space.
!</description>

!<input>
  ! Space-time vector to be evaluated
  type(t_spaceTimeVector), intent(IN) :: rspaceTimeVector
  
  ! Time where to evaluate
  real(DP), intent(IN) :: dtime
!</input>

!<inputoutput>
  ! A block vector in space that receives the result.
  ! If the block vector is empty, a new vector is created.
  type(t_vectorBlock), intent(INOUT) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP) :: dabstime,dntimesteps,dreltime
    integer :: itimestep,ntimesteps
    
    ! Rescale the time to the interval 0..ntimeintervals
    ntimesteps = rspaceTimeVector%p_rtimeDiscr%nintervals
    dntimesteps = real(ntimesteps,DP)
    call mprim_linearRescale(dtime,&
        rspaceTimeVector%p_rtimeDiscr%dtimeInit,&
        rspaceTimeVector%p_rtimeDiscr%dtimeMax,&
        0.0_DP,dntimesteps,dabstime)
    dabstime = min(max(0.0_DP,dabstime),1.0_DP)

    ! Is the destination vector empty? If yes, create a new one.
    if (rvector%NEQ .eq. 0) then
      call lsysbl_createVecBlockByDiscr (rspaceTimeVector%p_rspaceDiscr,&
          rvector,.false.)
    end if
    
    ! Clear the destination vector
    call lsysbl_clearVector (rvector)
    
    ! Now we have to take a look on the time discretisation.
    ! Depending on that, we choose the best suitable evaluation
    ! method.
    if (rspaceTimeVector%p_rtimeDiscr%ctype .eq. TDISCR_ONESTEPTHETA) then
      
      ! Get the time step which is closest to the time stamp.
      itimestep = int(dabstime + 0.5_DP)
      
      if (dabstime .eq. real(itimestep,DP)) then
        ! Nice coincidence, we have exactly timestep itimestep. Ok, then we
        ! can call the routine to get that timestep; this saves us some
        ! time as the interpolation can be omitted.
        call sptivec_getTimestepData (rspaceTimeVector, 1+itimestep, rvector)
        return
      end if
      
      if (rspaceTimeVector%NEQtime .eq. 2) then
        ! Special case: only one timestep!
        !
        ! Get the "relative" evaluation time; this in the interval 0..1
        dreltime = dabstime/dntimesteps
        
        ! Interpolate linearly.
        call interpolateLinear (dreltime,0,1,rspaceTimeVector,rvector)
      else
        ! Get the "relative" evaluation time; this in the interval -1..1
        dreltime = dabstime-real(itimestep,DP)
        
        ! Is this the first or the last timestep?
        if (itimestep .eq. 0) then
          ! First timestep. Interpolate between timesteps 0,1 and 2, evaluate
          ! near timestep 0.
          call interpolateQuadratic (dreltime-1.0_DP,0,1,2,rspaceTimeVector,rvector)
        else if (itimestep .eq. ntimesteps) then
          ! Last timestep. Interpolate between timesteps n-2,n-1 and n, evaluate
          ! near timestep n.
          call interpolateQuadratic (dreltime+1.0_DP,&
            ntimesteps-2,ntimesteps-1,ntimesteps,rspaceTimeVector,rvector)
        else
          ! Somewhere in the inner. Get the number of the previous and next timestep
          ! and interpolate there.
          call interpolateQuadratic (dreltime,&
            itimestep-1,itimestep,itimestep+1,rspaceTimeVector,rvector)
        end if
      end if
      
    else if (rspaceTimeVector%p_rtimeDiscr%ctype .eq. TDISCR_DG0) then
    
      ! Get the time step midpoint which is closest to the time stamp.
      itimestep = int(dabstime)
      
      if ((dabstime-0.5_DP) .eq. real(itimestep,DP)) then
        ! Nice coincidence, we have exactly timestep itimestep. Ok, then we
        ! can call the routine to get that timestep; this saves us some
        ! time as the interpolation can be omitted.
        call sptivec_getTimestepData (rspaceTimeVector, 1+itimestep, rvector)
        return
      end if
        
      ! Get the "relative" evaluation time; this in the interval -1..1
      dreltime = dabstime-0.5_DP-real(itimestep,DP)
      
      if (rspaceTimeVector%NEQtime .eq. 2) then
        ! Special case: only one timestep!
        ! Get the one and only solution
        call sptivec_getTimestepData (rspaceTimeVector, 1+0, rvector)
        
      else if (rspaceTimeVector%NEQtime .eq. 3) then
        ! Special case: only two timestep!
        ! Interpolate linearly.
        call interpolateLinear (dreltime,0,1,rspaceTimeVector,rvector)
        
      else
        ! Is this the first or the last timestep?
        if (itimestep .eq. 0) then
          ! First timestep. Interpolate between timesteps 0,1 and 2, evaluate
          ! near timestep 0.
          call interpolateQuadratic (dreltime-1.0_DP,0,1,2,rspaceTimeVector,rvector)
        else if (itimestep .eq. ntimesteps-1) then
          ! Last timestep. Interpolate between timesteps n-2,n-1 and n, evaluate
          ! near timestep n.
          call interpolateQuadratic (dreltime+1.0_DP,&
            ntimesteps-1-2,ntimesteps-1-1,ntimesteps-1,rspaceTimeVector,rvector)
        else
          ! Somewhere in the inner. Get the number of the previous and next timestep
          ! and interpolate there.
          call interpolateQuadratic (dreltime,&
            itimestep-1,itimestep,itimestep+1,rspaceTimeVector,rvector)
        end if
      end if
              
    else
      call output_line ("Unsupported time discretisation.", &
                        OU_CLASS_ERROR,OU_MODE_STD,"fetevl_evaluate")
      call sys_halt()
    end if
    
  contains
  
    subroutine interpolateLinear (dt,istep1,istep2,rspaceTimeVector,rvector)
    
    ! Calculates the linear interpolation of timestep istep1 and timestep
    ! istep2 of the space time vector rspaceTimeVector. The result is
    ! written to rvector.
    
    ! Interpolation weight; range 0..1 with 0~istep1 and 1~istep2.
    real(DP), intent(IN) :: dt
    
    ! "Left" time step corresponding to dt=0.
    integer, intent(IN) :: istep1

    ! "Right" time step corresponding to dt=1.
    integer, intent(IN) :: istep2
    
    ! Space time vector containing the data.
    type(t_spaceTimeVector), intent(IN) :: rspaceTimeVector
    
    ! Block vector receiving the result
    type(t_vectorBlock), intent(INOUT) :: rvector
    
      ! local variables
      type(t_vectorBlock) :: rvec
      
      ! Get a temporary vector
      call lsysbl_createVecBlockIndirect (rvector,rvec,.false.)
      
      ! Get the two timesteps
      call sptivec_getTimestepData (rspaceTimeVector, 1+istep1, rvec)
      call sptivec_getTimestepData (rspaceTimeVector, 1+istep2, rvector)
      
      ! Interpolate
      call lsysbl_vectorLinearComb (rvec,rvector,(1.0_DP-dt),dt)
      
      ! Release unused data
      call lsysbl_releaseVector (rvec)
      
    end subroutine
    
    ! ---------------------------------------------------------------
    
    subroutine interpolateQuadratic (dt,istep1,istep2,istep3,rspaceTimeVector,rvector)
    
    ! Calculates the linear interpolation of timestep istep1, istep2 and istep3
    ! of the space time vector rspaceTimeVector. The result is
    ! written to rvector.
    
    ! Interpolation weight; range -1..1 with 0~istep1 and 0~istep2 and 1~istep3
    real(DP), intent(IN) :: dt
    
    ! "Left" time step corresponding to dt=-1.
    integer, intent(IN) :: istep1

    ! "Middle" time step corresponding to dt=0.
    integer, intent(IN) :: istep2

    ! "Right" time step corresponding to dt=1.
    integer, intent(IN) :: istep3
    
    ! Space time vector containing the data.
    type(t_spaceTimeVector), intent(IN) :: rspaceTimeVector
    
    ! Block vector receiving the result
    type(t_vectorBlock), intent(INOUT) :: rvector
    
      ! local variables
      integer :: i
      type(t_vectorBlock) :: rvec1,rvec2
      real(DP), dimension(:), pointer :: p_Ddata1,p_Ddata2,p_Ddata3
      
      ! Get some temporary vectors
      call lsysbl_createVecBlockIndirect (rvector,rvec1,.false.)
      call lsysbl_createVecBlockIndirect (rvector,rvec2,.false.)
      
      ! Get the two timesteps
      call sptivec_getTimestepData (rspaceTimeVector, 1+istep1, rvec1)
      call sptivec_getTimestepData (rspaceTimeVector, 1+istep2, rvec2)
      call sptivec_getTimestepData (rspaceTimeVector, 1+istep3, rvector)

      ! Get the data arrays
      call lsysbl_getbase_double (rvec1,p_Ddata1)
      call lsysbl_getbase_double (rvec2,p_Ddata2)
      call lsysbl_getbase_double (rvector,p_Ddata3)
      
      ! Interpolate
      do i=1,size(p_Ddata3)
        call mprim_quadraticInterpolation (dt,&
            p_Ddata1(i),p_Ddata2(i),p_Ddata3(i),p_Ddata3(i))
      end do
      
      ! Release unused data
      call lsysbl_releaseVector (rvec2)
      call lsysbl_releaseVector (rvec1)
      
    end subroutine

  end subroutine

end module
