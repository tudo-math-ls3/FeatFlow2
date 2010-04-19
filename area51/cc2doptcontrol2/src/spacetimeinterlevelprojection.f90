!##############################################################################
!# ****************************************************************************
!# <name> spacetimeinterevelprojection </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises the interlevel projection (prolongation and 
!# restriction) of space-time coupled solution vectors. The t_sptiProjection
!# structure configures the interlevel projection. For the actual projection,
!# the following routines can be found in this module:
!#
!# 1.) sptipr_initProjection
!#     -> Initialises a space-time interlevel projection structure
!#
!# 2.) sptipr_doneProjection
!#     -> Cleans up a space-time interlevel projection structure
!#
!# 3.) sptipr_performProlongation
!#     -> Prolongation of a space-time coupled solution vector to a higher 
!#        level
!#
!# 4.) sptipr_performRestriction
!#     -> Retriction of a space-time coupled defect vector to a lower level
!#
!# 5.) sptipr_performInterpolation
!#     -> Interpolation of a space-time coupled solution vector to a lower
!#        level
!#
!# </purpose>
!##############################################################################

module spacetimeinterlevelprojection

  use fsystem
  use genoutput
  use storage
  use linearsystemscalar
  use linearsystemblock
  
  use spacetimevectors
  use multilevelprojection

  use spatialdiscretisation
  use spacetimehierarchy
  use multilevelprojection
  
  use timediscretisation

  implicit none
  
  private
  
  public :: t_sptiProjHierarchy
  public :: sptipr_initProjection
  public :: sptipr_doneProjection
  public :: sptipr_performProlongation
  public :: sptipr_performRestriction
  public :: sptipr_performInterpolation
  
!<types>

!<typeblock>
  
  ! A hierarchy of space-time projection definitions for a hierarchy
  ! of space-time levels.
  type t_sptiProjHierarchy
    ! Order of projection in time.
    ! =1: first order (linear), =2: second order (quadratic)
    integer :: itimeOrder = 0

    ! A space-time hierarchy that describes the discretisation in space and time
    ! for all levels.
    type(t_spaceTimeHierarchy), pointer :: p_rspaceTimeHierarchy

    ! Pointer to a hierarchy of interlevel projection structures in space.
    type(t_interlevelProjectionHier), pointer :: p_rprojHierarchySpace    
  end type

!</typeblock>

!</types>


contains

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_initProjection (rprojHier,rspaceTimeHierarchy,&
      rprojHierarchySpace,iorderTimeProlRest)
  
!<description>
  ! 
!</description>

!<input>
  ! A space-time hierarchy that describes the discretisation in space and time
  ! for all levels.
  type(t_spaceTimeHierarchy), intent(in), target :: rspaceTimeHierarchy
  
  ! Interlevel projection structure for space levels.
  type(t_interlevelProjectionHier), intent(in), target :: rprojHierarchySpace

  ! Order of the prolongation/restriction in time.
  ! =-1: Automatic
  ! =1: 1st order (bilinear)
  ! =2: 2nd order
  integer, intent(in), optional :: iorderTimeProlRest
  
!</input>

!<output>
  ! Space-time projection structure.
  type(t_sptiProjHierarchy), intent(OUT) :: rprojHier
!</output>

!</subroutine>

    ! Set itimeOrder <> 0 -> structure initialised.
    rprojHier%itimeOrder = -1
    if (present(iorderTimeProlRest)) rprojHier%itimeOrder = iorderTimeProlRest

    ! Remember the discretisation and projection hierarchy in space.  
    rprojHier%p_rspaceTimeHierarchy => rspaceTimeHierarchy
    rprojHier%p_rprojHierarchySpace => rprojHierarchySpace
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_doneProjection (rprojHier)
  
!<description>
  ! Releases memory allocated in sptipr_initProjection and cleans up rprojHier.
!</description>

!<inputoutput>
  ! A space/time interlevel projection structure that configures the 
  ! prolongation/restriction in space/time. 
  ! The structure is cleaned up.
  type(t_sptiProjHierarchy), intent(INOUT) :: rprojHier
!</output>

!</inputoutput>

    ! Clean up the spatial projection structure
    nullify(rprojHier%p_rprojHierarchySpace)
    nullify(rprojHier%p_rspaceTimeHierarchy)
    
    ! Set itimeOrder=0 -> structure not initialised anymore.
    rprojHier%itimeOrder = 0

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_performProlongation (rprojHier,ilevelfine,rcoarseVector, &
      rfineVector,rtempVecCoarse,rtempVecFine)
  
!<description>
  ! Performs a prolongation for a given space/time vector (i.e. a projection
  ! in the primal space where the solution lives). The vector
  ! rcoarseVector on a coarser space/time mesh is projected to the vector
  ! rfineVector on a finer space/time mesh. 
  ! rprojHier configures how the transfer is performed.
  ! This projection structure rprojHier must correspond to the space/time
  ! discretisation of rcoarseVector and rfineVector.
!</description>

!<input>
  ! A space/time interlevel projection structure that configures the 
  ! prolongation/restriction in space/time.
  type(t_sptiProjHierarchy), intent(IN) :: rprojHier
  
  ! Index of the fine mesh.
  integer, intent(in) :: ilevelfine

  ! Coarse grid vector
  type(t_spacetimeVector), intent(INOUT) :: rcoarseVector
!</input>

!<inputoutput>
  ! Temporary space-vector, specifying the discretisation and vector shape
  ! on the coarse grid.
  type(t_vectorBlock), intent(INOUT) :: rtempVecCoarse

  ! Temporary space-vector, specifying the discretisation and vector shape
  ! on the fine grid.
  type(t_vectorBlock), intent(INOUT) :: rtempVecFine
!</inputoutput>

!<output>
  ! Fine grid vector
  type(t_spacetimeVector), intent(INOUT) :: rfineVector
!</output>

!</subroutine>

    ! local variables
    integer :: istep,ispacelevelfine,ispacelevelcoarse,iorder
    type(t_vectorBlock) :: rx1,rx2,rx3
    type(t_vectorScalar) :: rtempVecFineScalar,rtempVecFineScalar2

    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Dx1,p_Dx3
    real(DP), dimension(:), pointer :: p_DtempVecCoarse,p_DtempVecFine
    real(DP), dimension(:), pointer :: p_DtempVecFineSca

    ! Get the space-time coarse and fine discretisations.
    call sth_getLevel (rprojHier%p_rspaceTimeHierarchy,ilevelfine-1,&
      ispaceLevel=ispacelevelcoarse)
    call sth_getLevel (rprojHier%p_rspaceTimeHierarchy,ilevelfine,&
      ispaceLevel=ispacelevelfine)

    ! Select the prolongation/restriction to use
    iorder = rprojHier%itimeOrder
    if (iorder .eq. -1) then
      ! automatic mode. Select the order based on the time stepping scheme.
      call tdiscr_getOrder(rcoarseVector%p_rtimeDiscr,iorder)
    end if
    
    select case (iorder)
    case (0)
      ! Constant prolongation
      !
      ! We need two more temp vectors:
      !
      ! One for the current timestep
      call lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.FALSE.)
      
      ! And one for the next timestep
      call lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.FALSE.)
      
      ! We need a scalar representation of the temp vector
      call lsysbl_createScalarFromVec (rtempVecFine,rtempVecFineScalar)
      
      ! DEBUG!!!
      call lsysbl_getbase_double (rx1,p_Dx1)
      call lsysbl_getbase_double (rx3,p_Dx3)
      call lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
      call lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
      call lsyssc_getbase_double (rtempVecFineScalar,p_DtempVecFineSca)

      ! Load timestep 0 into the temp vector and interpolate to the current level.
      ! Put the result to rx3.
      if (ispacelevelcoarse .ne. ispacelevelfine) then
        ! Space + time
        call sptivec_getTimestepData (rcoarseVector, 0, rtempVecCoarse)
        call mlprj_performProlongation (&
            rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
            rtempVecCoarse,rx3,rtempVecFineScalar)
      else
        ! Only time
        call sptivec_getTimestepData (rcoarseVector, 0, rx3)
      end if
      
      ! Save that vector as initial vector on the fine grid.
      call sptivec_setTimestepData (rfineVector, 0, rx3)
      
      if (rcoarseVector%NEQtime .NE. rfineVector%NEQtime) then
        
        ! Prolongation in time.
        !
        ! Loop through the time steps.
        do istep = 1,rcoarseVector%NEQtime
      
          ! rx3 was the vector from the last timestep. Shift it to rx1 and load
          ! the vector for the new timestep in rx3
          call lsysbl_copyVector (rx3,rx1)

          if (ispacelevelcoarse .ne. ispacelevelfine) then
            ! Space + time
            call sptivec_getTimestepData (rcoarseVector, istep, rtempVecCoarse)
            call mlprj_performProlongation (&
                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
                rtempVecCoarse,rx3,rtempVecFineScalar)
          else
            ! Only time
            call sptivec_getTimestepData (rcoarseVector, istep, rx3)
          end if
          
          ! Save that vector as new vector on the fine grid.
          call sptivec_setTimestepData (rfineVector, 2*istep, rx3)
          
          ! In the temp vector, create the prolongation. Use the constant prolongation.
          ! Primal vector of rx1 -> fine vector <- Dual vector of rx3.
          ! That's the prolongated vector.
          call lsysbl_copyVector (rx1,rtempVecFine)
          call lsyssc_copyVector (rx3%RvectorBlock(4),rtempVecFine%RvectorBlock(4))
          call lsyssc_copyVector (rx3%RvectorBlock(5),rtempVecFine%RvectorBlock(5))
          call lsyssc_copyVector (rx3%RvectorBlock(6),rtempVecFine%RvectorBlock(6))

          ! Save that vector as new vector on the fine grid.
          call sptivec_setTimestepData (rfineVector, 2*istep-1, rtempVecFine)
        
        end do
        
      else
      
        if (ispacelevelcoarse .ne. ispacelevelfine) then
        
          ! Prolongation only in space. 
      
          ! Loop through the time steps.
          do istep = 1,rcoarseVector%NEQtime
      
            ! Space prolongation
            call sptivec_getTimestepData (rcoarseVector, istep, rtempVecCoarse)
            call mlprj_performProlongation (&
                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
                rtempVecCoarse,rx3,rtempVecFineScalar)

            ! Save that vector as new vector on the fine grid.
            call sptivec_setTimestepData (rfineVector, istep, rx3)

          end do

        end if
      
      end if
      
      ! Release the temp vectors
      call lsyssc_releaseVector (rtempVecFineScalar)
      call lsysbl_releaseVector (rx3)
      call lsysbl_releaseVector (rx1)

    case (1,2)
      ! Linear prolongation:
      !
      ! We need two more temp vectors:
      !
      ! One for the current timestep
      call lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.false.)
      
      ! And one for the next timestep
      call lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.false.)
      
      ! We need a scalar representation of the temp vector
      call lsysbl_createScalarFromVec (rtempVecFine,rtempVecFineScalar)
      
      ! DEBUG!!!
      call lsysbl_getbase_double (rx1,p_Dx1)
      call lsysbl_getbase_double (rx3,p_Dx3)
      call lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
      call lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
      call lsyssc_getbase_double (rtempVecFineScalar,p_DtempVecFineSca)

      ! Load timestep 0 into the temp vector and interpolate to the current level.
      ! Put the result to rx3.
      if (ispacelevelcoarse .ne. ispacelevelfine) then
        ! Space + time
        call sptivec_getTimestepData (rcoarseVector, 1+0, rtempVecCoarse)
        call mlprj_performProlongation (&
            rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
            rtempVecCoarse,rx3,rtempVecFineScalar)
      else
        ! Only time
        call sptivec_getTimestepData (rcoarseVector, 1+0, rx3)
      end if
      
      ! Save that vector as initial vector on the fine grid.
      call sptivec_setTimestepData (rfineVector, 1+0, rx3)
      
      if (rcoarseVector%NEQtime .ne. rfineVector%NEQtime) then
        
        ! Prolongation in time.
        !
        ! Loop through the time steps.
        do istep = 1,rcoarseVector%NEQtime-1
      
          ! rx3 was the vector from the last timestep. Shift it to rx1 and load
          ! the vector for the new timestep in rx3
          call lsysbl_copyVector (rx3,rx1)

          if (ispacelevelcoarse .ne. ispacelevelfine) then
            ! Space + time
            call sptivec_getTimestepData (rcoarseVector, 1+istep, rtempVecCoarse)
            call mlprj_performProlongation (&
                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
                rtempVecCoarse,rx3,rtempVecFineScalar)
          else
            ! Only time
            call sptivec_getTimestepData (rcoarseVector, 1+istep, rx3)
          end if
          
          ! Save that vector as new vector on the fine grid.
          call sptivec_setTimestepData (rfineVector, 1+2*istep, rx3)
          
          ! In the temp vector, create the interpolation between rx1 and rx3.
          ! THat's the prolongated vector.
          call lsysbl_copyVector (rx1,rtempVecFine)
          call lsysbl_vectorLinearComb (rx3,rtempVecFine,0.5_DP,0.5_DP)

          ! Save that vector as new vector on the fine grid.
          call sptivec_setTimestepData (rfineVector, 1+2*istep-1, rtempVecFine)
        
        end do
        
      else
      
        if (ispacelevelcoarse .ne. ispacelevelfine) then
        
          ! Prolongation only in space. 
      
          ! Loop through the time steps.
          do istep = 1,rcoarseVector%NEQtime-1
      
            ! Space prolongation
            call sptivec_getTimestepData (rcoarseVector, 1+istep, rtempVecCoarse)
            call mlprj_performProlongation (&
                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
                rtempVecCoarse,rx3,rtempVecFineScalar)

            ! Save that vector as new vector on the fine grid.
            call sptivec_setTimestepData (rfineVector, 1+istep, rx3)

          end do

        end if
      
      end if
      
      ! Release the temp vectors
      call lsyssc_releaseVector (rtempVecFineScalar)
      call lsysbl_releaseVector (rx3)
      call lsysbl_releaseVector (rx1)
      
      ! DEBUG!!!
      !call cc_postprocSpaceTimeGMV (rproblem,rdiscrFine,rfineVector,'gmv/fine.gmv')
      !call cc_postprocSpaceTimeGMV (rproblem,rdiscrCoarse,rcoarseVector,'gmv/coarse.gmv')

    case (3)
      ! Quadratic prolongation.
      !
      ! We need vectors from three timesteps.
      call lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.false.)
      call lsysbl_createVecBlockIndirect (rtempVecFine,rx2,.false.)
      call lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.false.)
      
      ! We need a scalar representation of the temp vector
      call lsysbl_createScalarFromVec (rtempVecFine,rtempVecFineScalar)
      call lsysbl_createScalarFromVec (rtempVecFine,rtempVecFineScalar2)
      
      ! DEBUG!!!
      call lsysbl_getbase_double (rx1,p_Dx1)
      call lsysbl_getbase_double (rx3,p_Dx3)
      call lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
      call lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
      call lsyssc_getbase_double (rtempVecFineScalar,p_DtempVecFineSca)

      ! Load timestep 0 into the temp vector and interpolate to the current level.
      ! Put the result to rx3.
      if (ispacelevelcoarse .ne. ispacelevelfine) then
        ! Space + time
        call sptivec_getTimestepData (rcoarseVector, 1+0, rtempVecCoarse)
        call mlprj_performProlongation (&
            rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
            rtempVecCoarse,rx3,rtempVecFineScalar)
      else
        ! Only time
        call sptivec_getTimestepData (rcoarseVector, 1+0, rx3)
      end if
      
      ! Save vector rx1 as initial vector on the fine grid.
      call sptivec_setTimestepData (rfineVector, 1+0, rx3)
      
      if (rcoarseVector%NEQtime .ne. rfineVector%NEQtime) then
        
        ! Quadratic Prolongation in time.
        !
        !                v------- -0.125 -----------------+
        !     +------------------ -0.125 ------v
        !                v-- 0.75 --+-- 0.75 --v
        !     +- 0.375 --v                     v-- 0.375 -+
        !
        !     X----------o----------X----------o----------X
        !     rx1                  rx2                  rx3
        !
        ! Loop through the time steps. We loop in sets of size 2 and
        ! apply quadratic interpolation in each set.
        do istep = 1,rcoarseVector%NEQtime-2,2
      
          ! rx3 was the vector from the last timestep. Shift it to rx1 and load
          ! the vector for the new timestep patch into rx3, the data from the
          ! middle of the timestep patch to rx2.
          call lsysbl_copyVector (rx3,rx1)

          if (ispacelevelcoarse .ne. ispacelevelfine) then
            ! Space + time
            call sptivec_getTimestepData (rcoarseVector, 1+istep+1, rtempVecCoarse)
            call mlprj_performProlongation (&
                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
                rtempVecCoarse,rx3,rtempVecFineScalar)

            call sptivec_getTimestepData (rcoarseVector, 1+istep, rtempVecCoarse)
            call mlprj_performProlongation (&
                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
                rtempVecCoarse,rx2,rtempVecFineScalar)
          else
            ! Only time
            call sptivec_getTimestepData (rcoarseVector, 1+istep+1, rx3)
            call sptivec_getTimestepData (rcoarseVector, 1+istep, rx2)
          end if
          
          ! Save the vector as new vectors on the fine grid.
          call sptivec_setTimestepData (rfineVector, 1+2*istep, rx2)
          call sptivec_setTimestepData (rfineVector, 1+2*istep+2, rx3)
          
          ! In the temp vector, create the interpolation between rx1,rx2 and rx3.
          ! That's the prolongated vector.
          ! Save that vector as new vector on the fine grid.
          call lsysbl_copyVector (rx2,rtempVecFine)
          call lsysbl_vectorLinearComb (rx1,rtempVecFine,0.375_DP,0.75_DP)
          call lsysbl_vectorLinearComb (rx3,rtempVecFine,-0.125_DP,1.0_DP)
          call sptivec_setTimestepData (rfineVector, 1+2*istep-1, rtempVecFine)

          call lsysbl_copyVector (rx2,rtempVecFine)
          call lsysbl_vectorLinearComb (rx1,rtempVecFine,-0.125_DP,0.75_DP)
          call lsysbl_vectorLinearComb (rx3,rtempVecFine,0.375_DP,1.0_DP)
          call sptivec_setTimestepData (rfineVector, 1+2*istep+1, rtempVecFine)
        
        end do
        
        if (mod(rcoarseVector%NEQtime,2) .eq. 0) then
        
          ! Number of time solutions is equal, that means we have an odd number
          ! of timesteps -- and so we have to process the last timestep which
          ! could not be reached by the above loop as it's only half a patch.
        
          ! We proceed only half a patch and create a new patch from the
          ! time solutions NEQTIME, NEQTIME-1 and NEQTIME-2.
          ! The solution at NEQTIME-2 is currently present in rx2,
          ! the solution at NEQTIME-1 in rx3.
          call lsysbl_copyVector (rx2,rx1)
          call lsysbl_copyVector (rx3,rx2)

          if (ispacelevelcoarse .ne. ispacelevelfine) then
            ! Space + time
            call sptivec_getTimestepData (rcoarseVector, rcoarseVector%NEQtime, rtempVecCoarse)
            call mlprj_performProlongation (&
                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
                rtempVecCoarse,rx3,rtempVecFineScalar)
          else
            ! Only time
            call sptivec_getTimestepData (rcoarseVector, rcoarseVector%NEQtime, rx3)
          end if
          
          ! Save the last vector
          call sptivec_setTimestepData (rfineVector, rfineVector%NEQtime, rx3)

          ! Save the quadratic interpolation.          
          call lsysbl_copyVector (rx2,rtempVecFine)
          call lsysbl_vectorLinearComb (rx1,rtempVecFine,-0.125_DP,0.75_DP)
          call lsysbl_vectorLinearComb (rx3,rtempVecFine,0.375_DP,1.0_DP)
          call sptivec_setTimestepData (rfineVector, rfineVector%NEQtime-1, rtempVecFine)
        
        end if

      else
      
        if (ispacelevelcoarse .ne. ispacelevelfine) then
        
          ! Prolongation only in space. 
      
          ! Loop through the time steps.
          do istep = 1,rcoarseVector%NEQtime-1
      
            ! Space prolongation
            call sptivec_getTimestepData (rcoarseVector, 1+istep, rtempVecCoarse)
            call mlprj_performProlongation (&
                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
                rtempVecCoarse,rx3,rtempVecFineScalar)

            ! Save that vector as new vector on the fine grid.
            call sptivec_setTimestepData (rfineVector, 1+istep, rx3)

          end do

        end if
      
      end if
      
      ! Release the temp vectors
      call lsyssc_releaseVector (rtempVecFineScalar2)
      call lsyssc_releaseVector (rtempVecFineScalar)
      call lsysbl_releaseVector (rx3)
      call lsysbl_releaseVector (rx2)
      call lsysbl_releaseVector (rx1)

    end select

  end subroutine
     
  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_performRestriction (rprojHier,ilevelfine,rcoarseVector, &
      rfineVector,rtempVecCoarse,rtempVecFine)
  
!<description>
  ! Performs a restriction for a given space/time vector (i.e. a projection
  ! in the dual space where the RHS vector lives).The vector
  ! rfineVector on a finer grid is projected to the vector
  ! rcoarseVector on a coarser space/time mesh. 
  ! rprojHier configures how the transfer is performed.
  ! This projection structure rprojHier must correspond to the space/time
  ! discretisation of rcoarseVector and rfineVector.
!</description>

!<input>
  ! A space/time interlevel projection structure that configures the 
  ! prolongation/restriction in space/time.
  type(t_sptiProjHierarchy), intent(IN) :: rprojHier

  ! Index of the fine mesh.
  integer, intent(in) :: ilevelfine

  ! Fine grid vector
  type(t_spacetimeVector), intent(INOUT) :: rfineVector
!</input>

!<inputoutput>
  ! Temporary space-vector, specifying the discretisation and vector shape
  ! on the coarse grid.
  type(t_vectorBlock), intent(INOUT) :: rtempVecCoarse

  ! Temporary space-vector, specifying the discretisation and vector shape
  ! on the fine grid.
  type(t_vectorBlock), intent(INOUT) :: rtempVecFine
!</inputoutput>

!<output>
  ! Coarse grid vector
  type(t_spacetimeVector), intent(INOUT) :: rcoarseVector
!</output>

!</subroutine>

    ! local variables
    integer :: istep,ispacelevelfine,ispacelevelcoarse,iorder
    type(t_vectorBlock) :: rx1,rx2,rx3,rx4,rx5
    type(t_vectorScalar) :: rx1scalar
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Dx1,p_Dx3
    real(DP), dimension(:), pointer :: p_DtempVecCoarse,p_DtempVecFine

    ! Get the space-time coarse and fine discretisations.
    call sth_getLevel (rprojHier%p_rspaceTimeHierarchy,ilevelfine-1,&
      ispaceLevel=ispacelevelcoarse)
    call sth_getLevel (rprojHier%p_rspaceTimeHierarchy,ilevelfine,&
      ispaceLevel=ispacelevelfine)

    ! Select the prolongation/restriction to use
    iorder = rprojHier%itimeOrder
    if (iorder .eq. -1) then
      ! automatic mode. Select the order based on the time stepping scheme.
      call tdiscr_getOrder(rcoarseVector%p_rtimeDiscr,iorder)
    end if

    select case (iorder)
    case (0)
      ! Constant restriction:
      !
      ! Prolongation 'distributes' information from the coarse grid nodes
      ! to the fine grid nodes as follows:
      !
      ! Timestep:    n-1         n          n+1        (fine grid)
      !               x <--1/2-- X --1/2--> x
      !                         ^ \
      !                        /   \
      !                       +--1--+
      !
      ! Restriction is the 'adjoint' of the prolongation, so the corresponding
      ! matrix is the transposed. Written nodewise this means, that restriction
      ! has to 'collect' the distributed information with the same weights:
      !
      ! Timestep:    n-1         n          n+1        (fine grid)
      !               x --1/2--> X <--1/2-- x
      !                         ^ \
      !                        /   \
      !                       +--1--+
      !
      ! But because this is a restriction of a Finite-Difference RHS vector,
      ! the RHS must be divided by h/(h/2)=2 !

      ! One for the current timestep
      call lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.false.)
      
      ! And one for the next timestep
      call lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.false.)
      
      ! Create a scalar representation of rx1 for auxiliary reasons
      call lsysbl_createScalarFromVec (rx1,rx1Scalar)
     
      ! DEBUG!!!
      call lsysbl_getbase_double (rx1,p_Dx1)
      call lsysbl_getbase_double (rx3,p_Dx3)
      call lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
      call lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
         
      if (rcoarseVector%NEQtime .NE. rfineVector%NEQtime) then       
      
        ! Load timestep 0,1 (fine grid) into the temp vectors.
        call sptivec_getTimestepData (rfineVector, 0, rtempVecFine)
        call sptivec_getTimestepData (rfineVector, 1, rx3)
        
        ! Perform restriction for the first time step.
        !                          X <--1/2-- x
        !                         ^ \
        !                        /   \
        !                       +--1--+
        call lsyssc_vectorLinearComb (rx3%RvectorBlock(4),rtempVecFine%RvectorBlock(4),&
            0.5_DP/2.0_DP,1.0_DP/2.0_DP)
        call lsyssc_vectorLinearComb (rx3%RvectorBlock(5),rtempVecFine%RvectorBlock(5),&
            0.5_DP/2.0_DP,1.0_DP/2.0_DP)
        call lsyssc_vectorLinearComb (rx3%RvectorBlock(6),rtempVecFine%RvectorBlock(6),&
            0.5_DP/2.0_DP,1.0_DP/2.0_DP)
        
        ! Probably restrict the vector in space to the lower level.
        ! Save the result.
        if (ispacelevelcoarse .ne. ispacelevelfine) then
          ! Space + time
          call mlprj_performRestriction (&
                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
                rtempVecCoarse,rtempVecFine,rx1Scalar)
          call sptivec_setTimestepData (rcoarseVector, 0, rtempVecCoarse)
        else
          ! Only time
          call sptivec_setTimestepData (rcoarseVector, 0, rtempVecFine)
        end if
        
        ! Loop through the time steps -- on the coarse grid!
        do istep = 1,rcoarseVector%NEQtime-1
        
          ! rx3 was the 'right inner fine grid node x' in the above picture.
          ! Copy it to rx1, so it gets the new 'left inner fine grid node x'
          ! and load the new 'right inner fine grid node x' to rx1
          ! as well as the 'inner fine grid node X' to rtempVecFine.
          call lsysbl_copyVector (rx3,rx1)
          
          call sptivec_getTimestepData (rfineVector,2*istep, rtempVecFine)
          call sptivec_getTimestepData (rfineVector,2*istep+1, rx3)
          
          ! In rtempVecFine, create the restriction of the fine grid vectors
          ! rx1 and rx3 according to
          !
          ! Timestep:    n-1          n          n+1        (fine grid)
          !              rx3 --1/2--> X <--1/2-- rx3
          !             (old)        ^ \
          !             =rx1        /   \
          !                        +--1--+
          
          call lsyssc_vectorLinearComb (rx1%RvectorBlock(1),rtempVecFine%RvectorBlock(1),&
              0.5_DP,1.0_DP/2.0_DP)
          call lsyssc_vectorLinearComb (rx1%RvectorBlock(2),rtempVecFine%RvectorBlock(2),&
              0.5_DP,1.0_DP/2.0_DP)
          call lsyssc_vectorLinearComb (rx1%RvectorBlock(3),rtempVecFine%RvectorBlock(3),&
              0.5_DP,1.0_DP/2.0_DP)

          call lsyssc_vectorLinearComb (rx3%RvectorBlock(4),rtempVecFine%RvectorBlock(4),&
              0.5_DP,1.0_DP/2.0_DP)
          call lsyssc_vectorLinearComb (rx3%RvectorBlock(5),rtempVecFine%RvectorBlock(5),&
              0.5_DP,1.0_DP/2.0_DP)
          call lsyssc_vectorLinearComb (rx3%RvectorBlock(6),rtempVecFine%RvectorBlock(6),&
              0.5_DP,1.0_DP/2.0_DP)

          ! Probably restrict the vector in space to the lower level.
          ! Save the result.
          if (ispacelevelcoarse .ne. ispacelevelfine) then
            ! Space + time
            call mlprj_performRestriction (&
                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
                rtempVecCoarse,rtempVecFine,rx1Scalar)
            call sptivec_setTimestepData (rcoarseVector, istep, rtempVecCoarse)
          else
            ! Only time
            call sptivec_setTimestepData (rcoarseVector, istep, rtempVecFine)
          end if
        
        end do

        ! Last time step.

        ! rx3 is now the new 'left inner fine grid node x'.
        ! Load the 'inner fine grid node X' to rtempVecFine.
        call sptivec_getTimestepData (rfineVector,2*rcoarseVector%NEQtime, rtempVecFine)
        
        ! In rtempVecFine, create the restriction of the fine grid vectors
        ! rx1 and rx3 according to
        !
        ! Timestep:    n-1          n       (fine grid)
        !              rx3 --1/2--> X             
        !                          ^ \
        !                         /   \
        !                        +--1--+
        
        call lsyssc_vectorLinearComb (rx3%RvectorBlock(1),rtempVecFine%RvectorBlock(1),&
            0.5_DP,1.0_DP/2.0_DP)
        call lsyssc_vectorLinearComb (rx3%RvectorBlock(2),rtempVecFine%RvectorBlock(2),&
            0.5_DP,1.0_DP/2.0_DP)
        call lsyssc_vectorLinearComb (rx3%RvectorBlock(3),rtempVecFine%RvectorBlock(3),&
            0.5_DP,1.0_DP/2.0_DP)
        
        ! Probably restrict the vector in space to the lower level.
        ! Save the result.
        if (ispacelevelcoarse .ne. ispacelevelfine) then
          ! Space + time
          call mlprj_performRestriction (&
              rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
              rtempVecCoarse, rtempVecFine,rx1Scalar)
          call sptivec_setTimestepData (rcoarseVector, &
              rcoarseVector%NEQtime, rtempVecCoarse)
        else
          ! Only time
          call sptivec_setTimestepData (rcoarseVector, &
              rcoarseVector%NEQtime, rtempVecFine)
        end if
        
      else
      
        if (ispacelevelcoarse .ne. ispacelevelfine) then
        
          ! Restriction only in space.
          !
          ! Loop through the time steps
          do istep = 0,rcoarseVector%NEQtime
          
            ! Load the data
            call sptivec_getTimestepData (rfineVector,istep, rtempVecFine)

            ! Process the restriction and save the data to the coarse vector.
            call mlprj_performRestriction (&
                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
                rtempVecCoarse, rtempVecFine,rx1Scalar)
                
            call sptivec_setTimestepData (rcoarseVector, istep, rtempVecCoarse)
          
          end do
        
        end if
      
      end if

      ! Release the temp vectors
      call lsyssc_releaseVector (rx1Scalar)
      call lsysbl_releaseVector (rx3)
      call lsysbl_releaseVector (rx1)

    case (1,2)
      ! Linear restriction:
      !
      ! Prolongation 'distributes' information from the coarse grid nodes
      ! to the fine grid nodes as follows:
      !
      ! Timestep:    n-1         n          n+1        (fine grid)
      !               x <--1/2-- X --1/2--> x
      !                         ^ \
      !                        /   \
      !                       +--1--+
      !
      ! Restriction is the 'adjoint' of the prolongation, so the corresponding
      ! matrix is the transposed. Written nodewise this means, that restriction
      ! has to 'collect' the distributed information with the same weights:
      !
      ! Timestep:    n-1         n          n+1        (fine grid)
      !               x --1/2--> X <--1/2-- x
      !                         ^ \
      !                        /   \
      !                       +--1--+
      !
      ! But because this is a restriction of a Finite-Difference RHS vector,
      ! the RHS must be divided by h/(h/2)=2 !

      ! We need two more temp vectors:
      !
      ! One for the current timestep
      call lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.false.)
      
      ! And one for the next timestep
      call lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.false.)
      
      ! Create a scalar representation of rx1 for auxiliary reasons
      call lsysbl_createScalarFromVec (rx1,rx1Scalar)
     
      ! DEBUG!!!
      call lsysbl_getbase_double (rx1,p_Dx1)
      call lsysbl_getbase_double (rx3,p_Dx3)
      call lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
      call lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
         
      if (rcoarseVector%NEQtime .ne. rfineVector%NEQtime) then       
      
        ! Load timestep 0,1 (fine grid) into the temp vectors.
        call sptivec_getTimestepData (rfineVector, 1+0, rtempVecFine)
        call sptivec_getTimestepData (rfineVector, 1+1, rx3)
        
        ! Perform restriction for the first time step.
        !                          X <--1/2-- x
        !                         ^ \
        !                        /   \
        !                       +--1--+
        call lsysbl_vectorLinearComb (rx3,rtempVecFine,0.5_DP/2.0_DP,1.0_DP/2.0_DP)
        
        ! Probably restrict the vector in space to the lower level.
        ! Save the result.
        if (ispacelevelcoarse .ne. ispacelevelfine) then
          ! Space + time
          call mlprj_performRestriction (&
              rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
              rtempVecCoarse,rtempVecFine,rx1Scalar)
          call sptivec_setTimestepData (rcoarseVector, 1+0, rtempVecCoarse)
        else
          ! Only time
          call sptivec_setTimestepData (rcoarseVector, 1+0, rtempVecFine)
        end if
        
        ! Loop through the time steps -- on the coarse grid!
        do istep = 1,rcoarseVector%NEQtime-1-1
        
          ! rx3 was the 'right inner fine grid node x' in the above picture.
          ! Copy it to rx1, so it gets the new 'left inner fine grid node x'
          ! and load the new 'right inner fine grid node x' to rx1
          ! as well as the 'inner fine grid node X' to rtempVecFine.
          call lsysbl_copyVector (rx3,rx1)
          
          call sptivec_getTimestepData (rfineVector,1+2*istep, rtempVecFine)
          call sptivec_getTimestepData (rfineVector,1+2*istep+1, rx3)
          
          ! In rtempVecFine, create the restriction of the fine grid vectors
          ! rx1 and rx3 according to
          !
          ! Timestep:    n-1          n          n+1        (fine grid)
          !              rx3 --1/2--> X <--1/2-- rx3
          !             (old)        ^ \
          !                         /   \
          !                        +--1--+
          
          call lsysbl_vectorLinearComb (rx1,rtempVecFine,0.5_DP/2.0_DP,1.0_DP/2.0_DP)
          call lsysbl_vectorLinearComb (rx3,rtempVecFine,0.5_DP/2.0_DP,1.0_DP)
          
          ! Probably restrict the vector in space to the lower level.
          ! Save the result.
          if (ispacelevelcoarse .ne. ispacelevelfine) then
            ! Space + time
            call mlprj_performRestriction (&
                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
                rtempVecCoarse,rtempVecFine,rx1Scalar)
            call sptivec_setTimestepData (rcoarseVector, 1+istep, rtempVecCoarse)
          else
            ! Only time
            call sptivec_setTimestepData (rcoarseVector, 1+istep, rtempVecFine)
          end if
        
        end do

        ! Last time step.

        ! rx3 is now the new 'left inner fine grid node x'.
        ! Load the 'inner fine grid node X' to rtempVecFine.
        call sptivec_getTimestepData (rfineVector,1+2*(rcoarseVector%NEQtime-1), &
          rtempVecFine)
        
        ! In rtempVecFine, create the restriction of the fine grid vectors
        ! rx1 and rx3 according to
        !
        ! Timestep:    n-1          n       (fine grid)
        !              rx3 --1/2--> X             
        !                          ^ \
        !                         /   \
        !                        +--1--+
        
        call lsysbl_vectorLinearComb (rx3,rtempVecFine,0.5_DP/2.0_DP,1.0_DP/2.0_DP)
        
        ! Probably restrict the vector in space to the lower level.
        ! Save the result.
        if (ispacelevelcoarse .ne. ispacelevelfine) then
          ! Space + time
          call mlprj_performRestriction (&
              rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
              rtempVecCoarse, rtempVecFine,rx1Scalar)
          call sptivec_setTimestepData (rcoarseVector, &
              rcoarseVector%NEQtime, rtempVecCoarse)
        else
          ! Only time
          call sptivec_setTimestepData (rcoarseVector, &
              rcoarseVector%NEQtime, rtempVecFine)
        end if
        
      else
      
        if (ispacelevelcoarse .ne. ispacelevelfine) then
        
          ! Restriction only in space.
          !
          ! Loop through the time steps
          do istep = 0,rcoarseVector%NEQtime-1
          
            ! Load the data
            call sptivec_getTimestepData (rfineVector,1+istep, rtempVecFine)

            ! Process the restriction and save the data to the coarse vector.
            call mlprj_performRestriction (&
                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
                rtempVecCoarse, rtempVecFine,rx1Scalar)
                
            call sptivec_setTimestepData (rcoarseVector,1+istep, rtempVecCoarse)
          
          end do
        
        end if
      
      end if

      ! Release the temp vectors
      call lsyssc_releaseVector (rx1Scalar)
      call lsysbl_releaseVector (rx3)
      call lsysbl_releaseVector (rx1)

    case (3)
      ! Quadratic restriction:
      !
      ! Prolongation 'distributes' information from the coarse grid nodes
      ! to the fine grid nodes as follows:
      !
      !                v------- -0.125 -----------------+
      !     +------------------ -0.125 ------v
      !                v-- 0.75 --+-- 0.75 --v
      !     +- 0.375 --v                     v-- 0.375 -+
      !
      !     X----------o----------X----------o----------X
      ! Step1                   Step2                   Step3
      !
      ! Restriction is the 'adjoint' of the prolongation, so the corresponding
      ! matrix is the transposed. Written nodewise this means, that restriction
      ! has to 'collect' the distributed information with the same weights:
      !
      !                +------- -0.125 -----------------v
      !     v------------------ -0.125 ------+
      !                +-- 0.75 --v-- 0.75 --+
      !     v- 0.375 --+                     +-- 0.375 -v
      !
      !     X----------o----------X----------o----------X
      ! Step1                   Step2                   Step3
      !
      ! But because this is a restriction of a Finite-Difference RHS vector,
      ! the RHS must be divided by h/(h/2)=2 !

      ! We need some temporary vectors:
      call lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.false.)
      call lsysbl_createVecBlockIndirect (rtempVecFine,rx2,.false.)
      call lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.false.)
      call lsysbl_createVecBlockIndirect (rtempVecFine,rx4,.false.)
      call lsysbl_createVecBlockIndirect (rtempVecFine,rx5,.false.)
      
      ! Create a scalar representation of rx1 for auxiliary reasons
      call lsysbl_createScalarFromVec (rx1,rx1Scalar)
     
      ! DEBUG!!!
      call lsysbl_getbase_double (rx1,p_Dx1)
      call lsysbl_getbase_double (rx3,p_Dx3)
      call lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
      call lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
         
         
      ! DEBUG:
!      call lsysbl_clearVector (rtempVecFine,1.0_DP)
!      call sptivec_setTimestepData (rfineVector, 1+0, rtempVecFine)
!      call sptivec_setTimestepData (rfineVector, 1+1, rtempVecFine)
!      call sptivec_setTimestepData (rfineVector, 1+2, rtempVecFine)
!      call sptivec_setTimestepData (rfineVector, 1+3, rtempVecFine)
!      call sptivec_setTimestepData (rfineVector, 1+4, rtempVecFine)
!      call sptivec_setTimestepData (rfineVector, 1+5, rtempVecFine)
!      call sptivec_setTimestepData (rfineVector, 1+6, rtempVecFine)
!      call lsysbl_clearVector (rtempVecCoarse,-3.0_DP)
!      call sptivec_setTimestepData (rcoarseVector, 1+5, rtempVecCoarse)
!      call lsysbl_clearVector (rtempVecCoarse,1.0_DP)
!      call sptivec_setTimestepData (rcoarseVector, 1+1, rtempVecCoarse)
!      call sptivec_setTimestepData (rcoarseVector, 1+3, rtempVecCoarse)
!      call lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
      
      if (rcoarseVector%NEQtime .ne. rfineVector%NEQtime) then       
      
        ! Load timestep 0
        call sptivec_getTimestepData (rfineVector, 1+0, rx5)
        
        ! Prepare the temp vector
        call lsysbl_clearVector (rtempVecFine)
        
        ! Loop through the timesteps in chunks of size 2
        do istep = 0,rcoarseVector%NEQtime-3,2
        
          ! Get timestep i from the previous loop
          call lsysbl_copyVector (rx5,rx1)
        
          ! Load timestep i+1..i+4 of the fine mesh
          call sptivec_getTimestepData (rfineVector, 1+2*istep+1, rx2)
          call sptivec_getTimestepData (rfineVector, 1+2*istep+2, rx3)
          call sptivec_getTimestepData (rfineVector, 1+2*istep+3, rx4)
          call sptivec_getTimestepData (rfineVector, 1+2*istep+4, rx5)
          
          ! Do the restriction in the temp vector. The vector is already
          ! prepared from the last loop.
          if (istep .eq. 0) then
            ! In the very first step, put the (weighted) initial condition into
            ! the first coarse timestep vector. In later timesteps, rtempVecFine contains
            ! already the contribution from the last (left) interval and we only have
            ! to update it with the remaining contributions (from the right).
            call lsysbl_vectorLinearComb  (rx1,rtempVecFine,1.0_DP/2.0_DP,0.0_DP)
          end if
          call lsysbl_vectorLinearComb  (rx2,rtempVecFine,0.375_DP/2.0_DP,1.0_DP)
          call lsysbl_vectorLinearComb  (rx4,rtempVecFine,-0.125_DP/2.0_DP,1.0_DP)
          
          ! Probably restrict the vector in space to the lower level.
          ! Save the result.
          if (ispacelevelcoarse .ne. ispacelevelfine) then
            ! Space + time
            call mlprj_performRestriction (&
                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
                rtempVecCoarse,rtempVecFine,rx1Scalar)
            call sptivec_setTimestepData (rcoarseVector, 1+istep, rtempVecCoarse)
          else
            ! Only time
            call sptivec_setTimestepData (rcoarseVector, 1+istep, rtempVecFine)
          end if
          
          ! Second subvector in the patch. Prepare the temp vector
          call lsysbl_copyVector (rx3,rtempVecFine)
          call lsysbl_vectorLinearComb  (rx2,rtempVecFine,0.75_DP/2.0_DP,1.0_DP/2.0_DP)
          call lsysbl_vectorLinearComb  (rx4,rtempVecFine,0.75_DP/2.0_DP,1.0_DP)
          if (ispacelevelcoarse .ne. ispacelevelfine) then
            ! Space + time
            call mlprj_performRestriction (&
                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
                rtempVecCoarse,rtempVecFine,rx1Scalar)
            call sptivec_setTimestepData (rcoarseVector, 1+istep+1, rtempVecCoarse)
          else
            ! Only time
            call sptivec_setTimestepData (rcoarseVector, 1+istep+1, rtempVecFine)
          end if
          
          ! The 3rd subvector is not calculated here, but in the beginning of
          ! the next loop. We have to prepare the temp vector for that.
          call lsysbl_copyVector (rx5,rtempVecFine)
          call lsysbl_vectorLinearComb  (rx2,rtempVecFine,-0.125_DP/2.0_DP,1.0_DP/2.0_DP)
          call lsysbl_vectorLinearComb  (rx4,rtempVecFine,0.375_DP/2.0_DP,1.0_DP)
        
        end do
        
        if (mod(rcoarseVector%NEQtime,2) .ne. 0) then
          ! Equal number of timesteps / odd number of time-DOF's.
          ! Save the last vector, that's it.
          if (ispacelevelcoarse .ne. ispacelevelfine) then
            ! Space + time
            call mlprj_performRestriction (&
                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
                rtempVecCoarse,rtempVecFine,rx1Scalar)
            call sptivec_setTimestepData (rcoarseVector, rcoarseVector%NEQtime, rtempVecCoarse)
          else
            ! Only time
            call sptivec_setTimestepData (rcoarseVector, rcoarseVector%NEQtime, rtempVecFine)
          end if
        else
          ! Ok, the number of timesteps is even, that means there is only half a
          ! patch at the end. We have to shift the solutions by half a patch
          ! instead of a full patch and load the missing data.
          call lsysbl_copyVector (rx3,rx1)
          call lsysbl_copyVector (rx4,rx2)
          call lsysbl_copyVector (rx5,rx3)
          call sptivec_getTimestepData (rfineVector, rfineVector%NEQtime-1, rx4)
          call sptivec_getTimestepData (rfineVector, rfineVector%NEQtime, rx5)
          
          ! The temp vector rtempVecFine now corresponds to a partially prepared
          ! vector in the 'middle' of the patch like rx3. 
          call lsysbl_vectorLinearComb  (rx4,rtempVecFine,0.75_DP/2.0_DP,1.0_DP)
          if (ispacelevelcoarse .ne. ispacelevelfine) then
            ! Space + time
            call mlprj_performRestriction (&
                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
                rtempVecCoarse,rtempVecFine,rx1Scalar)
            call sptivec_setTimestepData (rcoarseVector, rcoarseVector%NEQtime-1, rtempVecCoarse)
          else
            ! Only time
            call sptivec_setTimestepData (rcoarseVector, rcoarseVector%NEQtime-1, rtempVecFine)
          end if
          
          ! Assemble the last vector in time.
          call lsysbl_copyVector (rx5,rtempVecFine)
          call lsysbl_vectorLinearComb (rx2,rtempVecFine,-0.125_DP/2.0_DP,1.0_DP/2.0_DP)
          call lsysbl_vectorLinearComb (rx4,rtempVecFine,0.375_DP/2.0_DP,1.0_DP)
        
          ! Save the last vector, that's it.
          if (ispacelevelcoarse .ne. ispacelevelfine) then
            ! Space + time
            call mlprj_performRestriction (&
                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
                rtempVecCoarse,rtempVecFine,rx1Scalar)
            call sptivec_setTimestepData (rcoarseVector, rcoarseVector%NEQtime, rtempVecCoarse)
          else
            ! Only time
            call sptivec_setTimestepData (rcoarseVector, rcoarseVector%NEQtime, rtempVecFine)
          end if

        end if
        
      end if

      ! DEBUG!!!
!      call lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
!      do istep = 1,rcoarseVector%NEQtime
!        call sptivec_getTimestepData (rcoarseVector, istep, rtempVecCoarse)
!      end do
      
      ! Release the temp vectors
      call lsyssc_releaseVector (rx1Scalar)
      call lsysbl_releaseVector (rx5)
      call lsysbl_releaseVector (rx4)
      call lsysbl_releaseVector (rx3)
      call lsysbl_releaseVector (rx2)
      call lsysbl_releaseVector (rx1)

    end select
    
!    call sptivec_saveToFileSequence(rcoarseVector,"(""rrhscoarse.txt."",I5.5)",.true.)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_performInterpolation (rprojHier,ilevelfine,rcoarseVector, &
      rfineVector,rtempVecCoarse,rtempVecFine)
  
!<description>
  ! Performs an interpolation for a given space/time vector to a lower level.
  ! The solution vector rfineVector on a finer grid is projected to the vector
  ! rcoarseVector on a coarser space/time mesh. 
  ! rprojHier configures how the transfer is performed.
  ! This projection structure rprojHier must correspond to the space/time
  ! discretisation of rcoarseVector and rfineVector.
!</description>

!<input>
  ! A space/time interlevel projection structure that configures the 
  ! prolongation/restriction in space/time.
  type(t_sptiProjHierarchy), intent(IN) :: rprojHier

  ! Index of the fine mesh.
  integer, intent(in) :: ilevelfine

  ! Fine grid vector
  type(t_spacetimeVector), intent(INOUT) :: rfineVector
!</input>

!<inputoutput>
  ! Temporary space-vector, specifying the discretisation and vector shape
  ! on the coarse grid.
  type(t_vectorBlock), intent(INOUT) :: rtempVecCoarse

  ! Temporary space-vector, specifying the discretisation and vector shape
  ! on the fine grid.
  type(t_vectorBlock), intent(INOUT) :: rtempVecFine
!</inputoutput>

!<output>
  ! Coarse grid vector
  type(t_spacetimeVector), intent(INOUT) :: rcoarseVector
!</output>

!</subroutine>

    ! local variables
    integer :: istep,ispacelevelcoarse,ispacelevelfine
    type(t_vectorBlock) :: rx1,rx3
    type(t_vectorScalar) :: rx1scalar
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Dx1,p_Dx3
    real(DP), dimension(:), pointer :: p_DtempVecCoarse,p_DtempVecFine

    ! Get the space-time coarse and fine discretisations.
    call sth_getLevel (rprojHier%p_rspaceTimeHierarchy,ilevelfine-1,&
      ispaceLevel=ispacelevelcoarse)
    call sth_getLevel (rprojHier%p_rspaceTimeHierarchy,ilevelfine,&
      ispaceLevel=ispacelevelfine)

    ! Prolongation 'distributes' information from the coarse grid nodes
    ! to the fine grid nodes as follows:
    !
    ! Timestep:    n-1         n          n+1        (fine grid)
    !               x <--1/2-- X --1/2--> x
    !                         ^ \
    !                        /   \
    !                       +--1--+
    !
    ! Interpolation is now taking the mean of adjacent nodes.
    ! Written nodewise this means:
    !
    ! Timestep:    n-1         n          n+1        (fine grid)
    !               x --1/4--> X <--1/4-- x
    !                         ^ \
    !                        /   \
    !                       +-1/2-+

    ! We need two more temp vectors:
    !
    ! One for the current timestep
    call lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.false.)
    
    ! And one for the next timestep
    call lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.false.)
    
    ! Create a scalar representation of rx1 for auxiliary reasons
    call lsysbl_createScalarFromVec (rx1,rx1Scalar)
   
    ! DEBUG!!!
    call lsysbl_getbase_double (rx1,p_Dx1)
    call lsysbl_getbase_double (rx3,p_Dx3)
    call lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
    call lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
       
    if (rcoarseVector%NEQtime .ne. rfineVector%NEQtime) then       
      ! Load timestep 0,1 (fine grid) into the temp vectors.
      call sptivec_getTimestepData (rfineVector, 1+0, rtempVecFine)
      call sptivec_getTimestepData (rfineVector, 1+1, rx3)
      
      ! Perform restriction for the first time step.
      !                          X <--1/4-- x
      !                         ^ \
      !                        /   \
      !                       +-3/4-+
      call lsysbl_vectorLinearComb (rx3,rtempVecFine,1._DP/4._DP,3._DP/4._DP)
      
      ! Probably restrict the vector in space to the lower level.
      ! Save the result.
      if (ispacelevelcoarse .ne. ispacelevelfine) then
        ! Space + time
        call mlprj_performInterpolation (&
            rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
            rtempVecCoarse,rtempVecFine,rx1Scalar)
        call sptivec_setTimestepData (rcoarseVector, 1+0, rtempVecCoarse)
      else
        ! Only time
        call sptivec_setTimestepData (rcoarseVector, 1+0, rtempVecFine)
      end if
      
      ! Loop through the time steps -- on the coarse grid!
      do istep = 1,rcoarseVector%NEQtime-1-1
      
        ! rx3 was the 'right inner fine grid node x' in the above picture.
        ! Copy it to rx1, so it gets the new 'left inner fine grid node x'
        ! and load the new 'right inner fine grid node x' to rx1
        ! as well as the 'inner fine grid node X' to rtempVecFine.
        call lsysbl_copyVector (rx3,rx1)
        
        call sptivec_getTimestepData (rfineVector,1+2*istep, rtempVecFine)
        call sptivec_getTimestepData (rfineVector,1+2*istep+1, rx3)
        
        ! In rtempVecFine, create the restriction of the fine grid vectors
        ! rx1 and rx3 according to
        !
        ! Timestep:    n-1          n          n+1        (fine grid)
        !              rx3 --1/4--> X <--1/4-- rx3
        !             (old)        ^ \
        !                         /   \
        !                        +-1/2-+
        
        call lsysbl_vectorLinearComb (rx1,rtempVecFine,0.25_DP,0.5_DP)
        call lsysbl_vectorLinearComb (rx3,rtempVecFine,0.25_DP,1.0_DP)
        
        ! Probably restrict the vector in space to the lower level.
        ! Save the result.
        if (ispacelevelcoarse .ne. ispacelevelfine) then
          ! Space + time
          call mlprj_performInterpolation (&
              rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
              rtempVecCoarse,rtempVecFine,rx1Scalar)
          call sptivec_setTimestepData (rcoarseVector, 1+istep, rtempVecCoarse)
        else
          ! Only time
          call sptivec_setTimestepData (rcoarseVector, 1+istep, rtempVecFine)
        end if
      
      end do

      ! Last time step.

      ! rx3 is now the new 'left inner fine grid node x'.
      ! Load the 'inner fine grid node X' to rtempVecFine.
      call sptivec_getTimestepData (rfineVector,1+2*(rcoarseVector%NEQtime-1), &
          rtempVecFine)
      
      ! In rtempVecFine, create the restriction of the fine grid vectors
      ! rx1 and rx3 according to
      !
      ! Timestep:    n-1          n       (fine grid)
      !              rx3 --1/4--> X             
      !                          ^ \
      !                         /   \
      !                        +-3/4-+
      
      call lsysbl_vectorLinearComb (rx3,rtempVecFine,1._DP/4._DP,3._DP/4._DP)
      
      ! Probably restrict the vector in space to the lower level.
      ! Save the result.
      if (ispacelevelcoarse .ne. ispacelevelfine) then
        ! Space + time
        call mlprj_performInterpolation (&
            rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
            rtempVecCoarse, rtempVecFine,rx1Scalar)
        call sptivec_setTimestepData (rcoarseVector, &
            rcoarseVector%NEQtime, rtempVecCoarse)
      else
        ! Only time
        call sptivec_setTimestepData (rcoarseVector, &
            rcoarseVector%NEQtime, rtempVecFine)
      end if

    else
    
      if (ispacelevelcoarse .ne. ispacelevelfine) then
      
        ! Interpolation only in space.
      
        !
        ! Loop through the time steps
        do istep = 0,rcoarseVector%NEQtime-1
        
          ! Load the data
          call sptivec_getTimestepData (rfineVector,1+istep, rtempVecFine)

          ! Process the restriction and save the data to the coarse vector.
          call mlprj_performInterpolation (&
              rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
              rtempVecCoarse, rtempVecFine,rx1Scalar)
              
          call sptivec_setTimestepData (rcoarseVector,1+istep, rtempVecCoarse)
        
        end do

      end if
    
    end if

    ! Release the temp vectors
    call lsyssc_releaseVector (rx1Scalar)
    call lsysbl_releaseVector (rx3)
    call lsysbl_releaseVector (rx1)

  end subroutine

end module
