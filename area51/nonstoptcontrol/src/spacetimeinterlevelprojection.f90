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
  
  use matrixio
  
  use physics

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

    ! Type of projection in time.
    ! =0: constant prolongation/restriction
    ! =1: linear prolongation/restriction. primal+dual solutions located
    !     at the endpoints of the time interval.
    ! =2: linear prolonggation/restriction, primal+dual solutions located
    !     according to the 1-step theta scheme. No time negotiation, so
    !     the restriction uses a constant interpolation of the primal
    !     space to synchronise in time with the dual space.
    ! =3: linear prolongation/restriction, primal+dual solutions located
    !     according to the 1-step theta scheme. Full linear restriction
    !     with interpolation in time.
    ! =4: Piecewise cubic interpolation in the inner. Piecewise quadratic
    !     extrapolation at the time endpoints. Solutions located at the
    !     time endpoints. (not completely implemented)
    ! =5: Piecewise linear with larger stencil. Solutions located at the
    !     time endpoints. (not completely implemented)
    ! =6: linear prolongation/restriction, primal+dual solutions located
    !     according to the 1-step theta scheme. No time negotiation, so
    !     the restriction uses a constant interpolation of the primal
    !     space to synchronise in time with the dual space. Simple
    !     restriction approach, i.e. R_p=P_p^T and R_d=P_d^T without
    !     respecting that the matrices must be exchanged.
    ! =7: linear prolongation/constant restriction. primal+dual solutions 
    !     located at the endpoints of the time interval.
    ! =8: linear prolongation/constant restriction. primal+dual solutions 
    !     located according to the timestep scheme.
    integer :: ctimeProjection = -1
    
    ! Underlying physics
    type(t_physics), pointer :: p_rphysics
    
    ! Pointer to the time coarse mesh.
    type(t_timeDiscretisation), pointer :: p_rtimeCoarseDiscr

    ! A space-time hierarchy that describes the discretisation in space and time
    ! for all levels.
    type(t_spaceTimeHierarchy), pointer :: p_rspaceTimeHierarchy

    ! Pointer to a hierarchy of interlevel projection structures in space.
    type(t_interlevelProjectionHier), pointer :: p_rprojHierarchySpace    
    
    ! An array of prolongation matrices.
    ! rprojectionMat(i) defines the weights for the projection
    ! between level i and level i+1.
    type(t_matrixScalar), dimension(:), pointer :: p_RprolongationMatPrimal
    type(t_matrixScalar), dimension(:), pointer :: p_RprolongationMatDual

    ! An array of restriction matrices.
    ! rprojectionMat(i) defines the weights for the projection
    ! between level i+1 and level i.
    type(t_matrixScalar), dimension(:), pointer :: p_RrestrictionMatPrimal
    type(t_matrixScalar), dimension(:), pointer :: p_RrestrictionMatDual
    
    ! An array of interpolation matrices.
    ! rprojectionMat(i) defines the weights for the projection
    ! between level i+1 and level i.
    type(t_matrixScalar), dimension(:), pointer :: p_RinterpolationMatPrimal
    type(t_matrixScalar), dimension(:), pointer :: p_RinterpolationMatDual

    ! An array of embedding matrices that embed the primal in the dual space
    ! and vice versa. Used to convert solutions between the spaces on every level.
    ! The matrices are saved transposed as they are only used that way.
    type(t_matrixScalar), dimension(:), pointer :: p_RembedMatPrimalInDualT
    type(t_matrixScalar), dimension(:), pointer :: p_RembedMatDualInPrimalT

  end type

!</typeblock>

!</types>


contains

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_initProjection (rprojHier,rspaceTimeHierarchy,&
      rprojHierarchySpace,rphysics,ctimeProjection)
  
!<description>
  ! 
!</description>

!<input>
  ! A space-time hierarchy that describes the discretisation in space and time
  ! for all levels.
  type(t_spaceTimeHierarchy), intent(in), target :: rspaceTimeHierarchy
  
  ! Interlevel projection structure for space levels.
  type(t_interlevelProjectionHier), intent(in), target :: rprojHierarchySpace

  ! Underlying physics
  type(t_physics), intent(in), target :: rphysics

  ! Type of projection in time.
  ! =-1: automatic.
  ! =0: constant prolongation/restriction
  ! =1: linear prolongation/restriction. primal+dual solutions located
  !     at the endpoints of the time interval.
  ! =2: linear prolonggation/restriction, primal+dual solutions located
  !     according to the 1-step theta scheme. No time negotiation, so
  !     the restriction uses a constant interpolation of the primal
  !     space to synchronise in time with the dual space.
  ! =3: linear prolonggation/restriction, primal+dual solutions located
  !     according to the 1-step theta scheme. Full linear restriction
  !     with interpolation in time.
  ! =4: Piecewise cubic interpolation in the inner. Piecewise quadratic
  !     extrapolation at the time endpoints. Solutions located at the
  !     time endpoints.
  ! =5: Piecewise linear with larger stencil. Solutions located at the
  !     time endpoints. (not completely implemented)
  ! =6: linear prolongation/restriction, primal+dual solutions located
  !     according to the 1-step theta scheme. No time negotiation, so
  !     the restriction uses a constant interpolation of the primal
  !     space to synchronise in time with the dual space. Simple
  !     restriction approach, i.e. R_p=P_p^T and R_d=P_d^T without
  !     respecting that the matrices must be exchanged.
  ! =7: linear prolongation/constant restriction. primal+dual solutions 
  !     located at the endpoints of the time interval.
  ! =8: linear prolongation/constant restriction. primal+dual solutions 
  !     located according to the timestep scheme.
  integer, intent(in), optional :: ctimeProjection
  
!</input>

!<output>
  ! Space-time projection structure.
  type(t_sptiProjHierarchy), intent(OUT) :: rprojHier
!</output>

!</subroutine>
    integer :: i,itimeOrder,ithetaschemetype

    ! Remember the physics; necessary so we know how and what to project
    rprojHier%p_rphysics => rphysics

    ! Remember the discretisation and projection hierarchy in space.  
    rprojHier%p_rspaceTimeHierarchy => rspaceTimeHierarchy
    rprojHier%p_rprojHierarchySpace => rprojHierarchySpace
    rprojHier%p_rtimeCoarseDiscr => rspaceTimeHierarchy%p_rtimeHierarchy%p_RtimeLevels(1)
    
    ithetaschemetype = rprojHier%p_rtimeCoarseDiscr%itag
        
    ! Set ctimeProjection
    rprojHier%ctimeProjection = -1
    if (present(ctimeProjection)) rprojHier%ctimeProjection = ctimeProjection

    if (rprojHier%ctimeProjection .eq. -1) then
      ! Automatic mode. Select the order based on the time stepping scheme.
      call tdiscr_getOrder(rprojHier%p_rtimeCoarseDiscr,itimeOrder)
      select case (itimeOrder)
      case (0)
        rprojHier%ctimeProjection = 0
      case (1)
        rprojHier%ctimeProjection = 1
      case (2)
        rprojHier%ctimeProjection = 3
      case default
        rprojHier%ctimeProjection = 1
      end select
      
      ! If our special 1-step scheme is activated, reduce the order to 1 in order
      ! to activate the corresponding prol/rest.
      if (ithetaschemetype .eq. 1) then
        rprojHier%ctimeProjection = 2
      end if
      
    end if
    
    ! Create prolongation and restriction matrices.
    allocate (rprojHier%p_RprolongationMatPrimal(max(1,rspaceTimeHierarchy%nlevels-1)))
    allocate (rprojHier%p_RprolongationMatDual(max(1,rspaceTimeHierarchy%nlevels-1)))
    allocate (rprojHier%p_RrestrictionMatPrimal(max(1,rspaceTimeHierarchy%nlevels-1)))
    allocate (rprojHier%p_RrestrictionMatDual(max(1,rspaceTimeHierarchy%nlevels-1)))
    allocate (rprojHier%p_RinterpolationMatPrimal(max(1,rspaceTimeHierarchy%nlevels-1)))
    allocate (rprojHier%p_RinterpolationMatDual(max(1,rspaceTimeHierarchy%nlevels-1)))
    
    nullify(rprojHier%p_RembedMatPrimalInDualT)
    nullify(rprojHier%p_RembedMatDualInPrimalT)
    select case (rprojHier%ctimeProjection)
    case (3)
      allocate (rprojHier%p_RembedMatPrimalInDualT(max(1,rspaceTimeHierarchy%nlevels)))
      allocate (rprojHier%p_RembedMatDualInPrimalT(max(1,rspaceTimeHierarchy%nlevels)))
    end select
    
    do i=1,rspaceTimeHierarchy%nlevels
      ! The restriction matrices are given as their transpose of the prolongation
      ! matrices.
      !
      ! WARNING!!!
      ! The primal restriction matrix is the transpose of the dual prolongation matrix!
      ! The dual restriction matrix is the transpose of the primal prolongation matrix!
      ! This is because the primal RHS is located at the timesteps of the dual
      ! solution (between the primal timesteps) and vice versa!!!
      select case (rprojHier%ctimeProjection)
      case (0:2,4,6)
      
        if (i .lt. rspaceTimeHierarchy%nlevels) then
          call sptipr_getProlMatrixPrimal(rspaceTimeHierarchy,i,rprojHier%ctimeProjection,&
              rprojHier%p_RprolongationMatPrimal(i))
              
          !call matio_writeMatrixHR (rprojHier%p_RprolongationMatPrimal(i), "pmat",&
          !    .true., 0, "matrixp."//trim(sys_siL(i,10)), "(E20.10)")
          
          call sptipr_getProlMatrixDual(rspaceTimeHierarchy,i,rprojHier%ctimeProjection,&
              rprojHier%p_RprolongationMatDual(i))

          !call matio_writeMatrixHR (rprojHier%p_RprolongationMatDual(i), "dmat",&
          !    .true., 0, "matrixd."//trim(sys_siL(i,10)), "(E20.10)")
              
          if (rprojHier%ctimeProjection .ne. 6) then
            ! Restriction matrices are obtained by transposing the prolongation
            ! matrices and exchanging the primal/dual matrices.
            call lsyssc_transposeMatrix (rprojHier%p_RprolongationMatPrimal(i),&
                rprojHier%p_RrestrictionMatDual(i),LSYSSC_TR_ALL)
                
            call lsyssc_transposeMatrix (rprojHier%p_RprolongationMatDual(i),&
                rprojHier%p_RrestrictionMatPrimal(i),LSYSSC_TR_ALL)
          else
            ! Restriction matrices are obtained by transposing the prolongation
            ! matrices and NOT exchanging the primal/dual matrices.
            !
            ! This is actually the wrong approach, but needed for tests.
            call lsyssc_transposeMatrix (rprojHier%p_RprolongationMatPrimal(i),&
                rprojHier%p_RrestrictionMatPrimal(i),LSYSSC_TR_ALL)
                
            call lsyssc_transposeMatrix (rprojHier%p_RprolongationMatDual(i),&
                rprojHier%p_RrestrictionMatDual(i),LSYSSC_TR_ALL)
          end if          
        end if

      case (7)
      
        if (i .lt. rspaceTimeHierarchy%nlevels) then
        
          ! Get constant restriction matrices.
          ! This is the interpolation matrix divided by 2.
          call sptipr_getInterpMatrixPrimal(rspaceTimeHierarchy,i,0,&
              rprojHier%p_RrestrictionMatPrimal(i))
          rprojHier%p_RrestrictionMatPrimal(i)%dscaleFactor = 2.0_DP
              
          call sptipr_getInterpMatrixPrimal(rspaceTimeHierarchy,i,0,&
              rprojHier%p_RrestrictionMatDual(i))
          rprojHier%p_RrestrictionMatDual(i)%dscaleFactor = 2.0_DP

          ! Get linear prolongation matrices.
          call sptipr_getProlMatrixPrimal(rspaceTimeHierarchy,i,1,&
              rprojHier%p_RprolongationMatPrimal(i))
              
          call sptipr_getProlMatrixDual(rspaceTimeHierarchy,i,1,&
              rprojHier%p_RprolongationMatDual(i))

        end if

      case (8)
      
        if (i .lt. rspaceTimeHierarchy%nlevels) then
        
          ! Get constant restriction matrices.
          ! This is the interpolation matrix divided by 2.
          call sptipr_getInterpMatrixPrimal(rspaceTimeHierarchy,i,0,&
              rprojHier%p_RrestrictionMatPrimal(i))
          rprojHier%p_RrestrictionMatPrimal(i)%dscaleFactor = 2.0_DP
              
          call sptipr_getInterpMatrixPrimal(rspaceTimeHierarchy,i,0,&
              rprojHier%p_RrestrictionMatDual(i))
          rprojHier%p_RrestrictionMatDual(i)%dscaleFactor = 2.0_DP

          ! Get linear prolongation matrices.
          call sptipr_getProlMatrixPrimal(rspaceTimeHierarchy,i,2,&
              rprojHier%p_RprolongationMatPrimal(i))
              
          call sptipr_getProlMatrixDual(rspaceTimeHierarchy,i,2,&
              rprojHier%p_RprolongationMatDual(i))

        end if

      case (3)
        ! Same as in the other case, but this needs additional embedding matrices.
        
        if (i .lt. rspaceTimeHierarchy%nlevels) then
          call sptipr_getProlMatrixPrimal(rspaceTimeHierarchy,i,rprojHier%ctimeProjection,&
              rprojHier%p_RprolongationMatPrimal(i))
              
          !call matio_writeMatrixHR (rprojHier%p_RprolongationMatPrimal(i), "pmat",&
          !    .true., 0, "matrixp."//trim(sys_siL(i,10)), "(E20.10)")
          
          call sptipr_getProlMatrixDual(rspaceTimeHierarchy,i,rprojHier%ctimeProjection,&
              rprojHier%p_RprolongationMatDual(i))

          !call matio_writeMatrixHR (rprojHier%p_RprolongationMatDual(i), "dmat",&
          !    .true., 0, "matrixd."//trim(sys_siL(i,10)), "(E20.10)")
              
          call lsyssc_transposeMatrix (rprojHier%p_RprolongationMatPrimal(i),&
              rprojHier%p_RrestrictionMatDual(i),LSYSSC_TR_ALL)
              
          call lsyssc_transposeMatrix (rprojHier%p_RprolongationMatDual(i),&
              rprojHier%p_RrestrictionMatPrimal(i),LSYSSC_TR_ALL)

        end if
            
        ! Get an additional embedding matrix
        call sptipr_getEmbedMatPrimalToDual (rspaceTimeHierarchy,i,rprojHier%ctimeProjection,&
            rprojHier%p_RembedMatPrimalInDualT(i),0)

        call sptipr_getEmbedMatDualToPrimal (rspaceTimeHierarchy,i,rprojHier%ctimeProjection,&
            rprojHier%p_RembedMatDualInPrimalT(i),0)

        call lsyssc_transposeMatrixInSitu(rprojHier%p_RembedMatPrimalInDualT(i))
        call lsyssc_transposeMatrixInSitu(rprojHier%p_RembedMatDualInPrimalT(i))

        ! The reverse is its transpose.
!        call lsyssc_transposeMatrix (rprojHier%p_RembedMatPrimalInDualT(i),&
!            rprojHier%p_RembedMatDualInPrimalT(i),LSYSSC_TR_ALL)
            
      case (5)
        ! NOT COMPLETELY IMPLEMENTED.
        ! Uses some different settings for the time interpolation matrices...

        if (i .lt. rspaceTimeHierarchy%nlevels) then
          call sptipr_getProlMatrixPrimal(rspaceTimeHierarchy,i,1,&
              rprojHier%p_RprolongationMatPrimal(i))
              
          call sptipr_getProlMatrixDual(rspaceTimeHierarchy,i,3,&
              rprojHier%p_RprolongationMatDual(i))

          call sptipr_getProlMatrixPrimal(rspaceTimeHierarchy,i,2,&
              rprojHier%p_RrestrictionMatDual(i))
              
          call sptipr_getProlMatrixDual(rspaceTimeHierarchy,i,4,&
              rprojHier%p_RrestrictionMatPrimal(i))
              
          call lsyssc_transposeMatrixInSitu(rprojHier%p_RrestrictionMatPrimal(i))
          call lsyssc_transposeMatrixInSitu(rprojHier%p_RrestrictionMatDual(i))
        end if
      end select
          
      if (i .lt. rspaceTimeHierarchy%nlevels) then

        ! The restriction matrices have to be divided by 2 as they are
        ! finite difference restrictions, not finite element restrictions!
        call lsyssc_scaleMatrix (rprojHier%p_RrestrictionMatPrimal(i),0.5_DP)
        call lsyssc_scaleMatrix (rprojHier%p_RrestrictionMatDual(i),0.5_DP)
        
  !      call lsyssc_getbase_double (rprojHier%p_RrestrictionMatPrimal(i),p_Da)
  !      call lsyssc_getbase_Kcol (rprojHier%p_RrestrictionMatPrimal(i),p_Kcol)
  !      call lsyssc_getbase_Kld (rprojHier%p_RrestrictionMatPrimal(i),p_Kld)
  !      p_Da(1) = 2.0_DP*p_Da(1)
  !      p_Da(rprojHier%p_RrestrictionMatPrimal(i)%NA) = 2.0_DP*p_Da(rprojHier%p_RrestrictionMatPrimal(i)%NA)

        ! Finally, calculate the interpolation matrices.
        call sptipr_getInterpMatrixPrimal(rspaceTimeHierarchy,i,rprojHier%ctimeProjection,&
            rprojHier%p_RinterpolationMatPrimal(i))

        !call matio_writeMatrixHR (rprojHier%p_RinterpolationMatPrimal(i), "pmat",&
        !    .true., 0, "imatrixp."//trim(sys_siL(i,10)), "(E20.10)")

        call sptipr_getInterpMatrixDual(rspaceTimeHierarchy,i,rprojHier%ctimeProjection,&
            rprojHier%p_RinterpolationMatDual(i))

        !call matio_writeMatrixHR (rprojHier%p_RinterpolationMatDual(i), "pmat",&
        !    .true., 0, "imatrixd."//trim(sys_siL(i,10)), "(E20.10)")
        
      end if
    end do

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

    integer :: i

    ! Release memory
    do i=1,rprojHier%p_rspaceTimeHierarchy%nlevels-1
      call lsyssc_releaseMatrix(rprojHier%p_RprolongationMatPrimal(i))
      call lsyssc_releaseMatrix(rprojHier%p_RprolongationMatDual(i))
      call lsyssc_releaseMatrix(rprojHier%p_RrestrictionMatPrimal(i))
      call lsyssc_releaseMatrix(rprojHier%p_RrestrictionMatDual(i))
      call lsyssc_releaseMatrix(rprojHier%p_RinterpolationMatPrimal(i))
      call lsyssc_releaseMatrix(rprojHier%p_RinterpolationMatDual(i))
    end do

    do i=1,rprojHier%p_rspaceTimeHierarchy%nlevels
      if (associated(rprojHier%p_RembedMatDualInPrimalT)) then
        call lsyssc_releaseMatrix(rprojHier%p_RembedMatDualInPrimalT(i))
      end if
      if (associated(rprojHier%p_RembedMatPrimalInDualT)) then
        call lsyssc_releaseMatrix(rprojHier%p_RembedMatPrimalInDualT(i))
      end if
    end do

    if (associated(rprojHier%p_RembedMatPrimalInDualT)) then
      deallocate (rprojHier%p_RembedMatPrimalInDualT)
      deallocate (rprojHier%p_RembedMatDualInPrimalT)
    end if

    deallocate (rprojHier%p_RprolongationMatPrimal)
    deallocate (rprojHier%p_RprolongationMatDual)
    deallocate (rprojHier%p_RrestrictionMatPrimal)
    deallocate (rprojHier%p_RrestrictionMatDual)
    deallocate (rprojHier%p_RinterpolationMatPrimal)
    deallocate (rprojHier%p_RinterpolationMatDual)
    
    ! Clean up the spatial projection structure
    nullify(rprojHier%p_rprojHierarchySpace)
    nullify(rprojHier%p_rspaceTimeHierarchy)
    nullify(rprojHier%p_rtimeCoarseDiscr)
    
    ! Set itimeOrder=0 -> structure not initialised anymore.
    rprojHier%ctimeProjection = -1

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_getProlMatrixPrimal (rspaceTimeHierarchy,ilevel,ctimeProjection,rprolMatrix)
  
!<description>
  ! Creates a prolongation matrix for the primal space between time-level 
  ! ilevel and ilevel+1 of the given space-time hierarchy.
!</description>

!<input>
  ! Underlying space-time hierarchy.
  type(t_spacetimeHierarchy), intent(in) :: rspaceTimeHierarchy
  
  ! Id of the coarse level
  integer, intent(in) :: ilevel
  
  ! Type of prolongation.
  integer, intent(in) :: ctimeProjection
!</input>

!<output>
  ! Matrix with the weights how to combine the vectors on the coarse mesh
  ! to get vectors on the fine mesh. The columns in the matrix correspond
  ! to the spatial vectors in all the timesteps.
  type(t_matrixScalar), intent(out) :: rprolMatrix
!</output>

!</subroutine>

    ! local variables
    type(t_timeDiscretisation), pointer :: p_rtimeDiscrCoarse
    type(t_timeDiscretisation), pointer :: p_rtimeDiscrFine
    integer :: ndofCoarse,ndofFine
    integer, dimension(:), pointer :: p_Kld, p_Kcol
    real(DP), dimension(:), pointer :: p_Da
    integer :: irow,icol
    integer, dimension(10) :: Kcol
    real(DP), dimension(10) :: Da

    ! Get the time levels of the two levels where we have to interpolate
    ! inbetween. The number of timesteps gives us the size of the matrix.
    call sth_getLevel (rspaceTimeHierarchy,ilevel,p_rtimeDiscr=p_rtimeDiscrCoarse)
    call sth_getLevel (rspaceTimeHierarchy,ilevel+1,p_rtimeDiscr=p_rtimeDiscrFine)

    ! At first, create an empty matrix.
    ! #rows = #time dofs of the fine level.
    ! #columns = #time dofs of the coarse level.
    ndofCoarse = tdiscr_igetNDofGlob(p_rtimeDiscrCoarse)
    ndofFine = tdiscr_igetNDofGlob(p_rtimeDiscrFine)
    
    if (ndofCoarse .eq. ndofFine) then
  
      ! Coarse level = fine level.
      ! Use the identity matrix!
      call lsyssc_createDiagMatrixStruc (rprolMatrix,ndofCoarse,LSYSSC_MATRIX9)
      call lsyssc_initialiseIdentityMatrix (rprolMatrix)
    
    else

      call lsyssc_createEmptyMatrixStub (rprolMatrix,LSYSSC_MATRIX9,ndofFine,ndofCoarse)
        
      ! Now depending on the order, create the matrix.
      select case (ctimeProjection)
      
      case (0)
        ! Constant prolongation. Just shift the value from
        ! the coarse mesh to the fine mesh.
      
        rprolMatrix%NA = ndofFine
        call storage_new ('sptipr_getProlMatrixPrimal', 'KLD', &
            ndofFine+1, ST_INT, rprolMatrix%h_KLD,ST_NEWBLOCK_ZERO)
        call storage_new ('sptipr_getProlMatrixPrimal', 'KCOL', &
            rprolMatrix%NA, ST_INT, rprolMatrix%h_KCOL,ST_NEWBLOCK_ZERO)
        call storage_new ('sptipr_getProlMatrixPrimal', 'DA', &
            rprolMatrix%NA, ST_DOUBLE, rprolMatrix%h_Da,ST_NEWBLOCK_ZERO)

        call lsyssc_getbase_double (rprolMatrix,p_Da)
        call lsyssc_getbase_Kcol (rprolMatrix,p_Kcol)
        call lsyssc_getbase_Kld (rprolMatrix,p_Kld)
        
        ! Fill KLD.
        do irow = 1,ndofFine+1
          p_Kld(irow) = irow
        end do
        
        ! Fill KCOL and DA
        p_Da(1) = 1.0_DP
        p_Kcol(1) = 1
        
        do icol = 1,ndofCoarse-1
          p_Kcol(2+2*(icol-1)) = icol
          p_Da(2+2*(icol-1)) = 1.0_DP

          p_Kcol(3+2*(icol-1)) = icol+1
          p_Da(3+2*(icol-1)) = 1.0_DP
        end do

      case (1,2,3,6,7,8)
        ! Simplest case.
        ! For simplicity, we take the lineaer interpolation, which is also
        ! 2nd order. In contrast to case 2, this matrix has a reduced
        ! matrix stencil and works only for implicit Euler!
        ! 
        ! Timestep:     1                     2        (coarse grid)
        !               x                     x          ...
        !                 --1/2--> + <--1/2--   --1/2--> ...
        !               |          |          |
        !               | 1        |          | 1
        !               V          V          V
        !
        !               x          x          x        (fine grid)
        !               1          2          3
        !
        ! The corresponding matrix looks like this:
        !
        !   1   
        !   1/2 1/2 
        !       1   
        !       1/2 1/2
        !           1   
        !           1/2 1/2
        !               1   
        !   ...
        !
        ! and is organised in 2x2 blocks
        
        rprolMatrix%NA = 3*ndofCoarse-2
        call storage_new ('sptipr_getProlMatrixPrimal', 'KLD', &
            ndofFine+1, ST_INT, rprolMatrix%h_KLD,ST_NEWBLOCK_ZERO)
        call storage_new ('sptipr_getProlMatrixPrimal', 'KCOL', &
            rprolMatrix%NA, ST_INT, rprolMatrix%h_KCOL,ST_NEWBLOCK_ZERO)
        call storage_new ('sptipr_getProlMatrixPrimal', 'DA', &
            rprolMatrix%NA, ST_DOUBLE, rprolMatrix%h_Da,ST_NEWBLOCK_ZERO)

        call lsyssc_getbase_double (rprolMatrix,p_Da)
        call lsyssc_getbase_Kcol (rprolMatrix,p_Kcol)
        call lsyssc_getbase_Kld (rprolMatrix,p_Kld)
        
        ! Fill KLD.
        p_Kld(1) = 1
        do irow = 2,ndofCoarse
          p_Kld(2*irow-2) = 2+(irow-2)*3
          p_Kld(2*irow-1) = 4+(irow-2)*3
        end do
        p_Kld(ndoffine+1) = rprolMatrix%NA+1
        
        ! Fill KCOL and DA
        p_Da(1) = 1.0_DP
        p_Kcol(1) = 1
        
        do icol = 1,ndofCoarse-1
          p_Kcol(-1+3*icol+0) = icol
          p_Kcol(-1+3*icol+1) = icol+1
          p_Kcol(-1+3*icol+2) = icol+1
          
          p_Da(-1+3*icol+0) = 0.5_DP
          p_Da(-1+3*icol+1) = 0.5_DP
          p_Da(-1+3*icol+2) = 1.0_DP
        end do

      case (4)
        ! Piecewise quadratic interpolation at the beginning/end.
        ! Piecewise cubic interpolation in the inner.

        rprolMatrix%NA = 1 + 3 + max(0,ndofcoarse-3)*5 + 1 + 3 + 1
        call storage_new ('sptipr_getProlMatrixPrimal', 'KLD', &
            ndofFine+1, ST_INT, rprolMatrix%h_KLD,ST_NEWBLOCK_ZERO)
        call storage_new ('sptipr_getProlMatrixPrimal', 'KCOL', &
            rprolMatrix%NA, ST_INT, rprolMatrix%h_KCOL,ST_NEWBLOCK_ZERO)
        call storage_new ('sptipr_getProlMatrixPrimal', 'DA', &
            rprolMatrix%NA, ST_DOUBLE, rprolMatrix%h_Da,ST_NEWBLOCK_ZERO)

        call lsyssc_getbase_double (rprolMatrix,p_Da)
        call lsyssc_getbase_Kcol (rprolMatrix,p_Kcol)
        call lsyssc_getbase_Kld (rprolMatrix,p_Kld)
        
        ! Fill KCOL and DA
        p_Kld(1) = 1
        p_Kcol(1) = 1
        p_Da(1) = 1.0_DP
        
        p_Kld(2) = 2
        p_Kcol(2) = 1
        p_Kcol(3) = 2
        p_Kcol(4) = 3
        p_Da(2) = 3.0_DP/8.0_DP
        p_Da(3) = 3.0_DP/4.0_DP
        p_Da(4) = -1.0_DP/8.0_DP
        
        do icol = 2,ndofCoarse-2
          p_Kld(3+(icol-2)*2) = 5+(icol-2)*5+0
        
          p_Kcol(5+(icol-2)*5+0) = icol
          p_Da(5+(icol-2)*5+0) = 1.0_DP
          
          p_Kld(3+(icol-2)*2+1) = 5+(icol-2)*5+1

          p_Kcol(5+(icol-2)*5+1) = icol-1
          p_Kcol(5+(icol-2)*5+2) = icol
          p_Kcol(5+(icol-2)*5+3) = icol+1
          p_Kcol(5+(icol-2)*5+4) = icol+2
          
          p_Da(5+(icol-2)*5+1) = -1.0_DP/16.0_DP
          p_Da(5+(icol-2)*5+2) = 9.0_DP/16.0_DP
          p_Da(5+(icol-2)*5+3) = 9.0_DP/16.0_DP
          p_Da(5+(icol-2)*5+4) = -1.0_DP/16.0_DP
        end do

        p_Kld(ndoffine-2) = rprolMatrix%NA-4
        p_Kcol(rprolMatrix%NA-4) = ndofcoarse-1
        p_Da(rprolMatrix%NA-4) = 1.0
        
        p_Kld(ndoffine-1) = rprolMatrix%NA-3
        p_Kcol(rprolMatrix%NA-3) = ndofcoarse-2
        p_Kcol(rprolMatrix%NA-2) = ndofcoarse-1
        p_Kcol(rprolMatrix%NA-1) = ndofcoarse
        p_Da(rprolMatrix%NA-3) = -1.0_DP/8.0_DP
        p_Da(rprolMatrix%NA-2) = 3.0_DP/4.0_DP
        p_Da(rprolMatrix%NA-1) = 3.0_DP/8.0_DP

        p_Kld(ndoffine) = rprolMatrix%NA
        p_Kcol(rprolMatrix%NA) = ndofcoarse
        p_Da(rprolMatrix%NA) = 1.0_DP

        p_Kld(ndoffine+1) = rprolMatrix%NA+1
        
      case (5)
        ! Piecewise linear with larger stencil and interpolation to the
        ! points in time.
        call lsyssc_createEmptyMatrix9 (rprolMatrix,ndofFine,(ndofFine-5)*3+1+2+3+2+2,ndofCoarse)
        
        Da(1) = 1.0_DP
        Kcol(1) = 1
        call lsyssc_setRowMatrix9 (rprolMatrix,1,1,Kcol,Da)
        
        Da(1:2) = (/0.5_DP,0.5_DP/)
        Kcol(1:2) = (/1,2/)
        call lsyssc_setRowMatrix9 (rprolMatrix,2,2,Kcol,Da)
        
        Da(1:3) = (/0.25_DP,0.375_DP,0.375_DP/)
        Kcol(1:3) = (/1,2,3/)
        call lsyssc_setRowMatrix9 (rprolMatrix,3,3,Kcol,Da)
        
        do irow = 3,ndofcoarse-1
          Kcol(1:3) = (/irow-1,irow,irow+1/)
          
          Da(1:3) = (/0.375_DP,0.5_DP,0.125_DP/)
          call lsyssc_setRowMatrix9 (rprolMatrix,4+2*(irow-3),3,Kcol,Da)
          
          Da(1:3) = (/0.125_DP,0.5_DP,0.375_DP/)
          call lsyssc_setRowMatrix9 (rprolMatrix,4+2*(irow-3)+1,3,Kcol,Da)
        end do

        Kcol(1:2) = (/ndofcoarse-1,ndofcoarse/)
        
        Da(1:2) = (/0.25_DP,0.75_DP/)
        call lsyssc_setRowMatrix9 (rprolMatrix,ndofFine-1,2,Kcol,Da)

        Da(1:2) = (/-0.25_DP,1.25_DP/)
        call lsyssc_setRowMatrix9 (rprolMatrix,ndofFine,2,Kcol,Da)

  !    case (2)
  !      ! Slightly harder case.
  !      ! Take linear interpolation but enlarge the blocks a bit
  !      ! such that the block size is aligned with that of the dual
  !      ! prolongation.
  !      !
  !      ! The matrix looks like this:
  !      !
  !      !   1
  !      !   1/2 1/2 0
  !      !       1   0
  !      !       1/2 1/2  
  !      !           1   0
  !      !           1/2 1/2  
  !      !               1   0
  !      !   ...
  !      !               1   0
  !      !               1/2 1/2
  !      !               0   1  
  !      !
  !      ! and is organised in 2x2 blocks with a special case at the 
  !      ! beginning and at the end.
  !      
  !      rprolMatrix%NA = 4+(ndofFine-2)*2
  !      call storage_new ('sptipr_getProlMatrixPrimal', 'KLD', &
  !          ndofFine+1, ST_INT, rprolMatrix%h_KLD,ST_NEWBLOCK_ZERO)
  !      call storage_new ('sptipr_getProlMatrixPrimal', 'KCOL', &
  !          rprolMatrix%NA, ST_INT, rprolMatrix%h_KCOL,ST_NEWBLOCK_ZERO)
  !      call storage_new ('sptipr_getProlMatrixPrimal', 'DA', &
  !          rprolMatrix%NA, ST_DOUBLE, rprolMatrix%h_Da,ST_NEWBLOCK_ZERO)
  !
  !      call lsyssc_getbase_double (rprolMatrix,p_Da)
  !      call lsyssc_getbase_Kcol (rprolMatrix,p_Kcol)
  !      call lsyssc_getbase_Kld (rprolMatrix,p_Kld)
  !      
  !      ! Fill KLD.
  !      p_Kld(1) = 1
  !      p_Kld(2) = 2
  !      do irow = 3,ndofFine+1
  !        p_Kld(irow) = 5+(irow-3)*2
  !      end do
  !      
  !      ! Set up the first two rows.
  !      p_Da(1) = 1.0_DP
  !      p_Da(2) = 0.5_DP
  !      p_Da(3) = 0.5_DP
  !      
  !      p_Kcol(1) = 1
  !      p_Kcol(2) = 1
  !      p_Kcol(3) = 2
  !      p_Kcol(4) = 3
  !      
  !      ! Then the remaining rows except for the end.
  !      do icol = 2,ndofCoarse-1
  !        p_Kcol(5+4*(icol-2)+0) = icol
  !        p_Kcol(5+4*(icol-2)+1) = icol+1
  !        p_Kcol(5+4*(icol-2)+2) = icol
  !        p_Kcol(5+4*(icol-2)+3) = icol+1
  !
  !        p_Da(5+4*(icol-2)+0) = 1.0_DP
  !        p_Da(5+4*(icol-2)+2) = 0.5_DP
  !        p_Da(5+4*(icol-2)+3) = 0.5_DP
  !      end do
  !
  !      ! The last row.
  !      p_Kcol(5+4*(ndofCoarse-2)+0) = ndofCoarse-1
  !      p_Kcol(5+4*(ndofCoarse-2)+1) = ndofCoarse
  !
  !      p_Da(5+4*(ndofCoarse-2)+1) = 1.0_DP
        
      end select
      
    end if
            
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_getProlMatrixDual (rspaceTimeHierarchy,ilevel,ctimeProjection,rprolMatrix)
  
!<description>
  ! Creates a prolongation matrix for the dual space between time-level 
  ! ilevel and ilevel+1 of the given space-time hierarchy.
!</description>

!<input>
  ! Underlying space-time hierarchy.
  type(t_spacetimeHierarchy), intent(in) :: rspaceTimeHierarchy
  
  ! Id of the coarse level
  integer, intent(in) :: ilevel
  
  ! Type of the prolongation.
  integer, intent(in) :: ctimeProjection
!</input>

!<output>
  ! Matrix with the weights how to combine the vectors on the coarse mesh
  ! to get vectors on the fine mesh. The columns in the matrix correspond
  ! to the spatial vectors in all the timesteps.
  type(t_matrixScalar), intent(out) :: rprolMatrix
!</output>

!</subroutine>

    ! local variables
    type(t_timeDiscretisation), pointer :: p_rtimeDiscrCoarse
    type(t_timeDiscretisation), pointer :: p_rtimeDiscrFine
    integer :: ndofCoarse,ndofFine
    integer, dimension(:), pointer :: p_Kld, p_Kcol
    real(DP), dimension(:), pointer :: p_Da
    integer :: irow,icol
    real(DP) :: dtheta

    ! Get the time levels of the two levels where we have to interpolate
    ! inbetween. The number of timesteps gives us the size of the matrix.
    call sth_getLevel (rspaceTimeHierarchy,ilevel,p_rtimeDiscr=p_rtimeDiscrCoarse)
    call sth_getLevel (rspaceTimeHierarchy,ilevel+1,p_rtimeDiscr=p_rtimeDiscrFine)

    ! At first, create an empty matrix.
    ! #rows = #time dofs of the fine level.
    ! #columns = #time dofs of the coarse level.
    ndofCoarse = tdiscr_igetNDofGlob(p_rtimeDiscrCoarse)
    ndofFine = tdiscr_igetNDofGlob(p_rtimeDiscrFine)

    if (ndofCoarse .eq. ndofFine) then
  
      ! Coarse level = fine level.
      ! Use the identity matrix!
      call lsyssc_createDiagMatrixStruc (rprolMatrix,ndofCoarse,LSYSSC_MATRIX9)
      call lsyssc_initialiseIdentityMatrix (rprolMatrix)
    
    else

      call lsyssc_createEmptyMatrixStub (rprolMatrix,LSYSSC_MATRIX9,ndofFine,ndofCoarse)
          
      ! Now depending on the order, create the matrix.
      select case (ctimeProjection)

      case (2,3,6,7)
        ! Slightly harder case.
        ! Because the dual solution is aligned between the primal,
        ! the prolongation stencil looks different.
        ! There is still a linear interpolation involved,
        ! but it has to be evaluated at different points.
        !
        ! The restriction is the adjoint of the prolongation.
        ! For Crank-Nicolson, we have:
        !
        !     [1]         [-----2-----]           [-----3-----] |         [-----4-----]
        !     X-----------o-----------X-----------o-----------X-----------o-----------X
        !     y0                      y1                      y2                      y3
        !     xi0                     xi1                     xi2                     xi3
        !     l0          l1                      l2                      l3
        !     p0          p1                      p2                      p3          
        !
        ! The prolongation will be the weighted mean for y and xi:
        !
        !     y0                      y1                      y2                      y3
        !     xi0                     xi1                     xi2                     xi3
        !     X  --1/2--> + <--1/2--- X ---1/2--> + <--1/2--- X ---1/2--> + <--1/2--- X
        !     |           |           |           |           |           |           |
        !     V           V           V           V           V           V           V

        !     X-----------X-----------X-----------X-----------X-----------X-----------X
        !     y0         y1           y2          y3          y4          y5          y6
        !     xi0        xi1          xi2         xi3         xi4         xi5         xi6
        !
        ! The situation is more complicated for the primal pressure and dual velocity.
        ! We use linear interpolation in the intervals and linear extrapolation at the
        ! endpoints of the interval.
        !
        !
        !
        !     l0          l1                      l2                      l3
        !     p0          p1                      p2                      p3          
        !     X-----------o-----------X-----------o-----------X-----------o-----------X
        !                  -3/4> <---- 1/4 ------- -3/4> <----- 1/4 ------ 
        !
        !                  ----- 1/4 ------> <3/4- ----- 1/4 ------> <3/4-
        !            <5/4--                                               --5/4> 
        !            <--------- -1/4 ------------- ---- -1/4 ------------------> 
        !
        !     X-----o-----X-----o-----X-----o-----X-----o-----X-----o-----X-----o-----X
        !     l0    l1          l2          l3          l4          l5          l6
        !     p0    p1          p2          p3          p4          p5          p6
        !     [1]   [--2--]     [--3--]     [--4--]     [--5--]     [--6--]     [--7--]
        !
        ! So at the end, we have:
        !
        !     [1]   [--2--]     [--3--]     [--4--]     [--5--]     [--6--]     [--7--]
        !     X-----o-----X-----o-----X-----o-----X-----o-----X-----o-----X-----o-----X
        !     y0          y1          y2          y3          y4          y5          y6
        !     xi0         xi1         xi2         xi3         xi4         xi5         xi6
        !     l0    l1          l2          l3          l4          l5          l6
        !     p0    p1          p2          p3          p4          p5          p6
        !
        ! The matrix therefore looks like this (for Crank-Nicolson):
        !
        !   1
        !   0   5/4 -1/4
        !       3/4 1/4
        !       1/4 3/4  
        !           3/4 1/4
        !           1/4 -3/4 
        !               3/4 1/4
        !   ...
        !               1/4 3/4
        !              -1/4 5/4
        
        dtheta = rspaceTimeHierarchy%p_rtimeHierarchy%p_RtimeLevels(1)%dtheta
        
        rprolMatrix%NA = 4+(ndofFine-2)*2
        call storage_new ('sptipr_getProlMatrixPrimal', 'KLD', &
            ndofFine+1, ST_INT, rprolMatrix%h_KLD,ST_NEWBLOCK_ZERO)
        call storage_new ('sptipr_getProlMatrixPrimal', 'KCOL', &
            rprolMatrix%NA, ST_INT, rprolMatrix%h_KCOL,ST_NEWBLOCK_ZERO)
        call storage_new ('sptipr_getProlMatrixPrimal', 'DA', &
            rprolMatrix%NA, ST_DOUBLE, rprolMatrix%h_Da,ST_NEWBLOCK_ZERO)

        call lsyssc_getbase_double (rprolMatrix,p_Da)
        call lsyssc_getbase_Kcol (rprolMatrix,p_Kcol)
        call lsyssc_getbase_Kld (rprolMatrix,p_Kld)
        
        ! Fill KLD.
        p_Kld(1) = 1
        p_Kld(2) = 2
        do irow = 3,ndofFine+1
          p_Kld(irow) = 5+(irow-3)*2
        end do
        
        ! Set up the first two rows.
        p_Da(1) = 1.0_DP
        p_Da(2) = 0.0_DP
        p_Da(3) = 1.0_DP+0.5_DP*dtheta  ! 1.25_DP
        p_Da(4) = -0.5_DP*dtheta        ! -0.25_DP

        p_Kcol(1) = 1
        p_Kcol(2) = 1
        p_Kcol(3) = 2
        p_Kcol(4) = 3
        
        ! Then the remaining rows except for the end.
        do icol = 2,ndofCoarse-1
          p_Kcol(5+4*(icol-2)+0) = icol
          p_Kcol(5+4*(icol-2)+1) = icol+1
          p_Kcol(5+4*(icol-2)+2) = icol
          p_Kcol(5+4*(icol-2)+3) = icol+1

          p_Da(5+4*(icol-2)+0) = 0.5_DP + 0.5_DP*dtheta ! 0.75_DP
          p_Da(5+4*(icol-2)+1) = 0.5_DP - 0.5_DP*dtheta ! 0.25_DP
          p_Da(5+4*(icol-2)+2) = 0.5_DP*dtheta          ! 0.25_DP
          p_Da(5+4*(icol-2)+3) = 1.0_DP - 0.5_DP*dtheta ! 0.75_DP
        end do

        ! The last row.
        p_Kcol(5+4*(ndofCoarse-2)+0) = ndofCoarse-1
        p_Kcol(5+4*(ndofCoarse-2)+1) = ndofCoarse

        p_Da(5+4*(ndofCoarse-2)) = -0.5_DP+0.5_DP*dtheta  ! -0.25_DP
        p_Da(5+4*(ndofCoarse-2)+1) = 1.5_DP-0.5_DP*dtheta ! 1.25_DP
        
      case (5)
        call sptipr_getProlMatrixPrimal (rspaceTimeHierarchy,ilevel,5,rprolMatrix)

      case default
        call sptipr_getProlMatrixPrimal (rspaceTimeHierarchy,ilevel,ctimeProjection,rprolMatrix)

      end select
      
    end if
            
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_getInterpMatrixPrimal (rspaceTimeHierarchy,ilevel,&
      ctimeProjection,rprolMatrix)
  
!<description>
  ! Creates a interpolation matrix for the primal space between time-level 
  ! ilevel and ilevel+1 of the given space-time hierarchy.
!</description>

!<input>
  ! Underlying space-time hierarchy.
  type(t_spacetimeHierarchy), intent(in) :: rspaceTimeHierarchy
  
  ! Id of the coarse level
  integer, intent(in) :: ilevel
  
  ! Type of the projection
  integer, intent(in) :: ctimeProjection
!</input>

!<output>
  ! Matrix with the weights how to combine the vectors on the coarse mesh
  ! to get vectors on the fine mesh. The columns in the matrix correspond
  ! to the spatial vectors in all the timesteps.
  type(t_matrixScalar), intent(out) :: rprolMatrix
!</output>

!</subroutine>

    ! local variables
    type(t_timeDiscretisation), pointer :: p_rtimeDiscrCoarse
    type(t_timeDiscretisation), pointer :: p_rtimeDiscrFine
    integer :: ndofCoarse,ndofFine
    integer, dimension(:), pointer :: p_Kld, p_Kcol
    real(DP), dimension(:), pointer :: p_Da
    integer :: irow,icol

    ! Get the time levels of the two levels where we have to interpolate
    ! inbetween. The number of timesteps gives us the size of the matrix.
    call sth_getLevel (rspaceTimeHierarchy,ilevel,p_rtimeDiscr=p_rtimeDiscrCoarse)
    call sth_getLevel (rspaceTimeHierarchy,ilevel+1,p_rtimeDiscr=p_rtimeDiscrFine)

    ! At first, create an empty matrix.
    ! #rows = #time dofs of the fine level.
    ! #columns = #time dofs of the coarse level.
    ndofCoarse = tdiscr_igetNDofGlob(p_rtimeDiscrCoarse)
    ndofFine = tdiscr_igetNDofGlob(p_rtimeDiscrFine)

    if (ndofCoarse .eq. ndofFine) then
  
      ! Coarse level = fine level.
      ! Use the identity matrix!
      call lsyssc_createDiagMatrixStruc (rprolMatrix,ndofCoarse,LSYSSC_MATRIX9)
      call lsyssc_initialiseIdentityMatrix (rprolMatrix)
    
    else
      call lsyssc_createEmptyMatrixStub (rprolMatrix,LSYSSC_MATRIX9,ndofCoarse,ndofFine)
          
      ! Now depending on the order, create the matrix.
      select case (ctimeProjection)
      case (0,1,2,3,4,6,7,8)
        ! The interpolation is just taking the values in the points in time.
        !
        ! Timestep: 
        !               1          2          3      
        !               x          x          x        (fine grid)
        !
        !               |                     |
        !               | 1                   | 1
        !               V                     V
        !
        !               x                     x          ...
        !               1                     2        (coarse grid)
        !
        !
        ! The matrix looks like this:
        !
        !   1   
        !   .   .   1
        !           .   .   1
        !   ...
        !
        ! 
        
        rprolMatrix%NA = ndofCoarse
        call storage_new ('sptipr_getProlMatrixPrimal', 'KLD', &
            ndofCoarse+1, ST_INT, rprolMatrix%h_KLD,ST_NEWBLOCK_ZERO)
        call storage_new ('sptipr_getProlMatrixPrimal', 'KCOL', &
            rprolMatrix%NA, ST_INT, rprolMatrix%h_KCOL,ST_NEWBLOCK_ZERO)
        call storage_new ('sptipr_getProlMatrixPrimal', 'DA', &
            rprolMatrix%NA, ST_DOUBLE, rprolMatrix%h_Da,ST_NEWBLOCK_ZERO)

        call lsyssc_getbase_double (rprolMatrix,p_Da)
        call lsyssc_getbase_Kcol (rprolMatrix,p_Kcol)
        call lsyssc_getbase_Kld (rprolMatrix,p_Kld)
        
        ! Fill KLD.
        p_Kld(1) = 1
        do irow = 1,ndofCoarse+1
          p_Kld(irow) = irow
        end do
        
        ! Fill KCOL and DA
        do icol = 1,ndofCoarse
          p_Kcol(icol) = 1+2*(icol-1)
          p_Da(icol) = 1.0_DP
        end do

      end select
    
    end if
            
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_getInterpMatrixDual (rspaceTimeHierarchy,ilevel,ctimeProjection,rprolMatrix)
  
!<description>
  ! Creates a prolongation matrix for the dual space between time-level 
  ! ilevel and ilevel+1 of the given space-time hierarchy.
!</description>

!<input>
  ! Underlying space-time hierarchy.
  type(t_spacetimeHierarchy), intent(in) :: rspaceTimeHierarchy
  
  ! Id of the coarse level
  integer, intent(in) :: ilevel
  
  ! Type of the projection.
  integer, intent(in) :: ctimeProjection
!</input>

!<output>
  ! Matrix with the weights how to combine the vectors on the coarse mesh
  ! to get vectors on the fine mesh. The columns in the matrix correspond
  ! to the spatial vectors in all the timesteps.
  type(t_matrixScalar), intent(out) :: rprolMatrix
!</output>

!</subroutine>

    ! local variables
    type(t_timeDiscretisation), pointer :: p_rtimeDiscrCoarse
    type(t_timeDiscretisation), pointer :: p_rtimeDiscrFine
    integer :: ndofCoarse,ndofFine
    integer, dimension(:), pointer :: p_Kld, p_Kcol
    real(DP), dimension(:), pointer :: p_Da
    integer :: irow,icol
    real(DP) :: dtheta

    ! Get the time levels of the two levels where we have to interpolate
    ! inbetween. The number of timesteps gives us the size of the matrix.
    call sth_getLevel (rspaceTimeHierarchy,ilevel,p_rtimeDiscr=p_rtimeDiscrCoarse)
    call sth_getLevel (rspaceTimeHierarchy,ilevel+1,p_rtimeDiscr=p_rtimeDiscrFine)

    ! At first, create an empty matrix.
    ! #rows = #time dofs of the fine level.
    ! #columns = #time dofs of the coarse level.
    ndofCoarse = tdiscr_igetNDofGlob(p_rtimeDiscrCoarse)
    ndofFine = tdiscr_igetNDofGlob(p_rtimeDiscrFine)

    if (ndofCoarse .eq. ndofFine) then
  
      ! Coarse level = fine level.
      ! Use the identity matrix!
      call lsyssc_createDiagMatrixStruc (rprolMatrix,ndofCoarse,LSYSSC_MATRIX9)
      call lsyssc_initialiseIdentityMatrix (rprolMatrix)
    
    else
      call lsyssc_createEmptyMatrixStub (rprolMatrix,LSYSSC_MATRIX9,ndofCoarse,ndofFine)
          
      ! Now depending on the order, create the matrix.
      select case (ctimeProjection)

      case (2,3,6)
        ! Slightly harder case.
        ! Because the dual solution is aligned between the primal,
        ! the interpolation stencil looks different.
        !
        ! For Crank-Nicolson, we have:
        !
        !     [1]   [--2--]     [--3--]     [--4--]     [--5--]     [--6--]     [--7--]
        !     X-----o-----X-----o-----X-----o-----X-----o-----X-----o-----X-----o-----X
        !     y0          y1          y2          y3          y4          y5          y6
        !     xi0         xi1         xi2         xi3         xi4         xi5         xi6
        !     l0    l1          l2          l3          l4          l5          l6
        !     p0    p1          p2          p3          p4          p5          p6
        !
        ! The primal solution for y and xi uses the standard interpolation:
        !
        !     X-----------X-----------X-----------X-----------X-----------X-----------X
        !     y0         y1           y2          y3          y4          y5          y6
        !     xi0        xi1          xi2         xi3         xi4         xi5         xi6
        !     |                       |                       |                       |
        !     V                       V                       V                       V
        !     X-----------------------X-----------------------X-----------------------X
        !     y0                      y1                      y2                      y3
        !     xi0                     xi1                     xi2                     xi3
        !
        ! The situation is more complicated for the primal pressure and dual velocity.
        ! Here, we average the solution of two neighboring points in time to get a value
        ! in the middle:
        !
        !     [1]   [--2--]     [--3--]     [--4--]     [--5--]     [--6--]     [--7--]
        !     l0    l1          l2          l3          l4          l5          l6
        !     p0    p1          p2          p3          p4          p5          p6
        !     X-----o-----X-----o-----X-----o-----X-----o-----X-----o-----X-----o-----X
        !           |-1/2> <1/2-|           |-1/2> <1/2-|           |-1/2> <1/2-|
        !                 |                       |                       |
        !                 V                       V                       V
        !     X-----------o-----------X-----------o-----------X-----------o-----------X
        !     l0          l1                      l2                      l3
        !     p0          p1                      p2                      p3          
        !
        ! So at the end, we have:
        !
        !     [1]         [-----2-----]           [-----3-----] |         [-----4-----]
        !     X-----------o-----------X-----------o-----------X-----------o-----------X
        !     y0                      y1                      y2                      y3
        !     xi0                     xi1                     xi2                     xi3
        !     l0          l1                      l2                      l3
        !     p0          p1                      p2                      p3          
        !
        !
        ! The matrix looks like this (for Crank-Nicolson):
        !
        !   1
        !       1/2 1/2
        !               1/2 1/2
        !   ...
        ! The idea is that we take a linear interpolation between
        ! the time nodes to get the value in the time midpoint.
        
        dtheta = rspaceTimeHierarchy%p_rtimeHierarchy%p_RtimeLevels(1)%dtheta
        
        rprolMatrix%NA = ndofFine
        call storage_new ('sptipr_getProlMatrixPrimal', 'KLD', &
            ndofCoarse+1, ST_INT, rprolMatrix%h_KLD,ST_NEWBLOCK_ZERO)
        call storage_new ('sptipr_getProlMatrixPrimal', 'KCOL', &
            rprolMatrix%NA, ST_INT, rprolMatrix%h_KCOL,ST_NEWBLOCK_ZERO)
        call storage_new ('sptipr_getProlMatrixPrimal', 'DA', &
            rprolMatrix%NA, ST_DOUBLE, rprolMatrix%h_Da,ST_NEWBLOCK_ZERO)

        call lsyssc_getbase_double (rprolMatrix,p_Da)
        call lsyssc_getbase_Kcol (rprolMatrix,p_Kcol)
        call lsyssc_getbase_Kld (rprolMatrix,p_Kld)
        
        ! Fill KLD.
        p_Kld(1) = 1
        do irow = 2,ndofCoarse+1
          p_Kld(irow) = 2+(irow-2)*2
        end do
        
        ! Set up the first row.
        p_Da(1) = 1.0_DP
        p_Kcol(1) = 1
        
        ! Then the remaining rows.
        do icol = 2,ndofCoarse
          p_Kcol(2+2*(icol-2)) = 2+2*(icol-2)
          p_Da(2+2*(icol-2)) = dtheta

          p_Kcol(2+2*(icol-2)+1) = 2+2*(icol-2)+1
          p_Da(2+2*(icol-2)+1) = 1.0_DP-dtheta
        end do

      case default
        ! Default case, standard handling. Take the matrix from the primal.
        call sptipr_getInterpMatrixPrimal (rspaceTimeHierarchy,ilevel,&
            ctimeProjection,rprolMatrix)

      end select
      
    end if
            
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_getEmbedMatPrimalToDual (&
      rspaceTimeHierarchy,ilevel,ctimeProjection,rmatrix,itype)
  
!<description>
  ! Creates an embedding matrix that embeds the primal solution in the dual space.
!</description>

!<input>
  ! Underlying space-time hierarchy.
  type(t_spacetimeHierarchy), intent(in) :: rspaceTimeHierarchy
  
  ! Id of the coarse level
  integer, intent(in) :: ilevel
  
  ! Type of the prolongation.
  integer, intent(in) :: ctimeProjection
  
  ! Type of the embedding matrix. Depends on the projection type.
  ! ctimeProjection = 3:
  !   itype = 0: standard linear embedding matrix.
  !   itype = 1: quadratic embedding matrix.
  integer, intent(in) :: itype
!</input>

!<output>
  ! Embedding matrix
  type(t_matrixScalar), intent(out) :: rmatrix
!</output>

!</subroutine>

    ! local variables
    type(t_timeDiscretisation), pointer :: p_rtimeDiscrFine
    integer :: ndofFine
    integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Kdiagonal
    !real(DP), dimension(:), pointer :: p_Da
    !integer :: irow
    integer :: icol
    real(DP) :: dtheta
    
    integer, dimension(3) :: Kcol
    real(DP), dimension(3) :: Da

    ! Get the time levels of the two levels where we have to interpolate
    ! inbetween. The number of timesteps gives us the size of the matrix.
    call sth_getLevel (rspaceTimeHierarchy,ilevel,p_rtimeDiscr=p_rtimeDiscrFine)

    ! At first, create an empty matrix.
    ! #rows = #time dofs of the fine level.
    ! #columns = #time dofs of the coarse level.
    ndofFine = tdiscr_igetNDofGlob(p_rtimeDiscrFine)
        
    ! Now depending on the order, create the matrix.
    select case (ctimeProjection)

    case (3)
      select case (itype)
      case (0)
        ! The embedding matrix for the primal space is just taking the average of
        ! the points in time in the interval plus a constant interpolation at the beginning.
        !
        !                             y0                      y1                      y2
        !                            xi0                     xi1                     xi2
        !                             X-----------|-----------X-----------|-----------X
        !                  <-1--|           |-1/2> <1/2-|           |-1/2> <1/2-|
        !                 |                       |                       |
        !                 V                       V                       V
        !                 X-----------|-----------X-----------|-----------X-----------|
        !               y-1/2                   y1/2                     y3/2
        !               xi-1/2                  xi1/2                    yi3/2        
        !
        ! The matrix looks like this (for Crank-Nicolson):
        !
        !   1
        !       1/2 1/2
        !               1/2 1/2
        !   ...
        ! The idea is that we take a linear interpolation between
        ! the time nodes to get the value in the time midpoint.

        dtheta = rspaceTimeHierarchy%p_rtimeHierarchy%p_RtimeLevels(1)%dtheta
        
        call lsyssc_createEmptyMatrix9 (rmatrix,ndofFine,2*(ndofFine-1)+1)
        
        Kcol(1) = 1
        Da(1) = 1.0_DP
        call lsyssc_setRowMatrix9 (rmatrix,1,1,Kcol,Da)
        
        Da(1:2) = (/1.0_DP-dtheta,dtheta/)
        do icol = 2,ndoffine
          Kcol(1:2) = (/icol-1,icol/)
          call lsyssc_setRowMatrix9 (rmatrix,icol,2,Kcol,Da)
        end do
        ! Calculate Kdiagonal
        call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
        call lsyssc_getbase_Kld (rmatrix,p_Kld)
        call lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)
        call lsyssc_rebuildKdiagonal (p_Kcol, p_Kld, p_Kdiagonal, ndofFine)
        
      case (1)
        ! Quadratic embedding matrix. Only for CN at the moment!!!
        call lsyssc_createEmptyMatrix9 (rmatrix,ndofFine,3*(ndofFine-1)+1)
        
        Kcol(1) = 1
        Da(1) = 1.0_DP
        call lsyssc_setRowMatrix9 (rmatrix,1,1,Kcol,Da)
        
        Kcol(1:3) = (/1,2,3/)
        Da(1:3) = (/0.375_DP,0.75_DP,-0.125_DP/)
        call lsyssc_setRowMatrix9 (rmatrix,2,3,Kcol,Da)

        Da(1:3) = (/-0.125_DP,0.75_DP,0.375_DP/)
        do icol = 3,ndoffine
          Kcol(1:3) = (/icol-2,icol-1,icol/)
          call lsyssc_setRowMatrix9 (rmatrix,icol,3,Kcol,Da)
        end do
        
        ! Calculate Kdiagonal
        call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
        call lsyssc_getbase_Kld (rmatrix,p_Kld)
        call lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)
        call lsyssc_rebuildKdiagonal (p_Kcol, p_Kld, p_Kdiagonal, ndofFine)
                    
      end select

    case default
      call lsyssc_releaseMatrix (rmatrix)

    end select
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_getEmbedMatDualToPrimal (&
      rspaceTimeHierarchy,ilevel,ctimeProjection,rmatrix,itype)
  
!<description>
  ! Creates an embedding matrix that embeds the dual solution in the primal space.
!</description>

!<input>
  ! Underlying space-time hierarchy.
  type(t_spacetimeHierarchy), intent(in) :: rspaceTimeHierarchy
  
  ! Id of the coarse level
  integer, intent(in) :: ilevel
  
  ! Type of the prolongation.
  integer, intent(in) :: ctimeProjection
  
  ! Type of the embedding matrix. Depends on the projection type.
  ! ctimeProjection = 3:
  !   itype = 0: standard linear embedding matrix.
  !   itype = 1: quadratic embedding matrix.
  integer, intent(in) :: itype
!</input>

!<output>
  ! Embedding matrix
  type(t_matrixScalar), intent(out) :: rmatrix
!</output>

!</subroutine>

    ! local variables
    type(t_timeDiscretisation), pointer :: p_rtimeDiscrFine
    integer :: ndofFine
    integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Kdiagonal
    !integer :: irow
    integer :: icol
    real(DP) :: dtheta
    
    integer, dimension(3) :: Kcol
    real(DP), dimension(3) :: Da

    ! Get the time levels of the two levels where we have to interpolate
    ! inbetween. The number of timesteps gives us the size of the matrix.
    call sth_getLevel (rspaceTimeHierarchy,ilevel,p_rtimeDiscr=p_rtimeDiscrFine)

    ! At first, create an empty matrix.
    ! #rows = #time dofs of the fine level.
    ! #columns = #time dofs of the coarse level.
    ndofFine = tdiscr_igetNDofGlob(p_rtimeDiscrFine)
        
    ! Now depending on the order, create the matrix.
    select case (ctimeProjection)

    case (3)
      select case (itype)
      case (0)
        dtheta = rspaceTimeHierarchy%p_rtimeHierarchy%p_RtimeLevels(1)%dtheta
        
        call lsyssc_createEmptyMatrix9 (rmatrix,ndofFine,2*(ndofFine-1))
        
        ! Constant in the first row
        Kcol(1) = 1
        Da(1) = 1.0_DP
        call lsyssc_setRowMatrix9 (rmatrix,1,1,Kcol,Da)
        
        ! Mean value in all intermediate rows
        Da(1:2) = (/dtheta,1.0_DP-dtheta/)
        do icol = 2,ndoffine-1
          Kcol(1:2) = (/icol,icol+1/)
          call lsyssc_setRowMatrix9 (rmatrix,icol,2,Kcol,Da)
        end do

        ! Constant in the last row
        Kcol(1) = 1
        Da(1) = 1.0_DP
        call lsyssc_setRowMatrix9 (rmatrix,ndoffine,1,Kcol,Da)
        
        ! Calculate Kdiagonal
        call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
        call lsyssc_getbase_Kld (rmatrix,p_Kld)
        call lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)
        call lsyssc_rebuildKdiagonal (p_Kcol, p_Kld, p_Kdiagonal, ndofFine)
        
      case (1)

        call output_line ("Not implemented.")
        call sys_halt()                    
      end select

    case default
      call lsyssc_releaseMatrix (rmatrix)

    end select
    
  end subroutine

  ! ***************************************************************************

    subroutine getSpatialVector (rprojHier,rcoarseVector,iindex,rx,&
        rtempVecCoarse,rtempVecFineScalar,&
        ispacelevelcoarse,ispacelevelfine)
    
    ! Extracts the spatial subvector iindex from rcoarseVector and puts it
    ! into rx. If necessary, the vector is prolongated in space.
    
    ! A space/time interlevel projection structure that configures the 
    ! prolongation/restriction in space/time.
    type(t_sptiProjHierarchy), intent(IN) :: rprojHier

    ! Space-time source vector
    type(t_spaceTimeVector), intent(in) :: rcoarseVector
    
    ! Index of the subvector
    integer, intent(in) :: iindex
    
    ! Space vector output.
    type(t_vectorBlock), intent(inout) :: rx
    
    ! Temp vectors
    type(t_vectorBlock), intent(inout) :: rtempVecCoarse
    type(t_vectorScalar), intent(inout) :: rtempVecFineScalar
    
    ! Level of the coarse and (destination) fine vector.
    integer, intent(in) :: ispacelevelcoarse
    integer, intent(in) :: ispacelevelfine
  
      ! Load timestep iindex into the temp vector and interpolate to the current level.
      ! Put the result to rx3.
      if (ispacelevelcoarse .ne. ispacelevelfine) then
        ! Space + time
        call sptivec_getTimestepData (rcoarseVector, iindex, rtempVecCoarse)
        call mlprj_performProlongation (&
            rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
            rtempVecCoarse,rx,rtempVecFineScalar)
      else
        ! Only time
        call sptivec_getTimestepData (rcoarseVector, iindex, rx)
      end if
      
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

    ! Local variables
    type(t_spaceTimeVectorAccess) :: raccessPool
    integer :: ispacelevelfine,ispacelevelcoarse,itimelevelfine,itimelevelcoarse
    type(t_vectorScalar) :: rtempVecFineScalar
    type(t_vectorBlock), pointer :: p_rx
    integer :: irow, icol
    type(t_matrixScalar), pointer :: p_rprolMatrixPrim,p_rprolMatrixDual
    real(DP), dimension(:), pointer :: p_DaPrim,p_DaDual,p_DdataCoarse,p_DdataFine
    integer, dimension(:), pointer :: p_KcolPrim, p_KldPrim, p_KcolDual, p_KldDual
    real(DP) :: dscalePrim, dscaleDual

    ! DEBUG!!!
    !do istep=1,rcoarseVector%NEQtime
    !  call sptivec_setSubvectorConstant(rcoarseVector,istep,real(istep,DP))
    !end do      

    ! Get the space-time coarse and fine discretisations.
    call sth_getLevel (rprojHier%p_rspaceTimeHierarchy,ilevelfine-1,&
      ispaceLevel=ispacelevelcoarse,itimeLevel=itimelevelcoarse)
    call sth_getLevel (rprojHier%p_rspaceTimeHierarchy,ilevelfine,&
      ispaceLevel=ispacelevelfine,itimeLevel=itimelevelfine)
    
    ! We need a scalar representation of the temp vector
    call lsysbl_createScalarFromVec (rtempVecFine,rtempVecFineScalar)
    
    ! Allocate temp memory
    select case (rprojHier%ctimeProjection)
    case (0)
      call sptivec_createAccessPool (rfineVector%p_rspaceDiscr,raccessPool,3)
    case (1)
      call sptivec_createAccessPool (rfineVector%p_rspaceDiscr,raccessPool,3)
    case (2,3,6,7,8)
      call sptivec_createAccessPool (rfineVector%p_rspaceDiscr,raccessPool,3)
    case (4)
      call sptivec_createAccessPool (rfineVector%p_rspaceDiscr,raccessPool,4)
    case default
      call sptivec_createAccessPool (rfineVector%p_rspaceDiscr,raccessPool,10)
    end select
    
    ! Prolongation means, we multiply with the prolongation matrix.
    ! y and xi have to be multiplied with the primal prol. matrix,
    ! lambda and p with the dual. Note that for implicit Euler,
    ! the primal and dual prolongation matrix is the same, while
    ! for CN, y and xi are at the same points in time as well as
    ! lambda and p: p is between y and thus at the same points
    ! in time as lambda, not as y!
    !
    ! Get the matrix.
    p_rprolMatrixPrim => rprojHier%p_RprolongationMatPrimal(ilevelfine-1)
    p_rprolMatrixDual => rprojHier%p_RprolongationMatDual(ilevelfine-1)
    
    call lsyssc_getbase_double (p_rprolMatrixPrim,p_DaPrim)
    call lsyssc_getbase_double (p_rprolMatrixDual,p_DaDual)

    call lsyssc_getbase_Kcol (p_rprolMatrixPrim,p_KcolPrim)
    call lsyssc_getbase_Kld (p_rprolMatrixPrim,p_KldPrim)

    call lsyssc_getbase_Kcol (p_rprolMatrixDual,p_KcolDual)
    call lsyssc_getbase_Kld (p_rprolMatrixDual,p_KldDual)
    
    dscalePrim = p_rprolMatrixPrim%dscaleFactor
    dscaleDual = p_rprolMatrixDual%dscaleFactor
    
    ! Apply the multiplication.
    ! The rows in the matrix correspond to the time fine mesh, the columns
    ! to the time coarse mesh.
    do irow = 1,p_rprolMatrixPrim%NEQ
    
      ! Clear the destination
      call lsysbl_clearVector (rtempVecFine)
      
      ! DEBUG!!!
      call lsysbl_getbase_double (rtempVecFine,p_DdataFine)
      
      select case (rprojHier%p_rphysics%cequation)
      case (0,2)
        ! Heat equation
      
        ! Primal space: y,xi
        do icol = p_KldPrim(irow),p_KldPrim(irow+1)-1
        
          ! Try to get the vector from the vector pool. Saves time.
          call sptivec_getVectorFromPool(raccessPool,p_KcolPrim(icol),p_rx)
          if (.not. associated(p_rx)) then
            ! No, we have to fetch/calculate the vector.
            !
            ! Get a buffer where to save it.
            call sptivec_getFreeBufferFromPool (raccessPool,p_KcolPrim(icol),p_rx)
            
            ! Read the source vector. The column of the matrix specifies
            ! the timestep.
            call getSpatialVector (rprojHier,rcoarseVector,p_KcolPrim(icol),p_rx,&
                rtempVecCoarse,rtempVecFineScalar,ispacelevelcoarse,ispacelevelfine)
          end if
          
          ! DEBUG!!!
          call lsysbl_getbase_double (p_rx,p_DdataCoarse)
          
          ! Now, rx is the time vector at timestep icol. Weighted multiplication
          ! into rtempVecFine for y and xi.
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(1),rtempVecFine%RvectorBlock(1),dscalePrim*p_DaPrim(icol),1.0_DP)
              
        end do
        
        ! Dual space: lambda/p
        do icol = p_KldDual(irow),p_KldDual(irow+1)-1
        
          ! Try to get the vector from the vector pool. Saves time.
          call sptivec_getVectorFromPool(raccessPool,p_KcolDual(icol),p_rx)
          if (.not. associated(p_rx)) then
            ! No, we have to fetch/calculate the vector.
            !
            ! Get a buffer where to save it.
            call sptivec_getFreeBufferFromPool (raccessPool,p_KcolDual(icol),p_rx)
            
            ! Read the source vector. The column of the matrix specifies
            ! the timestep.
            call getSpatialVector (rprojHier,rcoarseVector,p_KcolDual(icol),p_rx,&
                rtempVecCoarse,rtempVecFineScalar,ispacelevelcoarse,ispacelevelfine)
          end if
          
          ! DEBUG!!!
          call lsysbl_getbase_double (p_rx,p_DdataCoarse)
          
          ! Now, rx is the time vector at timestep icol. Weighted multiplication
          ! into rtempVecFine for y and xi.
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(2),rtempVecFine%RvectorBlock(2),dscaleDual*p_DaDual(icol),1.0_DP)
              
        end do
        
      case (1)
        ! Stokes equation
      
        ! Primal space: y,xi
        do icol = p_KldPrim(irow),p_KldPrim(irow+1)-1
        
          ! Try to get the vector from the vector pool. Saves time.
          call sptivec_getVectorFromPool(raccessPool,p_KcolPrim(icol),p_rx)
          if (.not. associated(p_rx)) then
            ! No, we have to fetch/calculate the vector.
            !
            ! Get a buffer where to save it.
            call sptivec_getFreeBufferFromPool (raccessPool,p_KcolPrim(icol),p_rx)
            
            ! Read the source vector. The column of the matrix specifies
            ! the timestep.
            call getSpatialVector (rprojHier,rcoarseVector,p_KcolPrim(icol),p_rx,&
                rtempVecCoarse,rtempVecFineScalar,ispacelevelcoarse,ispacelevelfine)
          end if
          
          ! DEBUG!!!
          call lsysbl_getbase_double (p_rx,p_DdataCoarse)
          
          ! Now, rx is the time vector at timestep icol. Weighted multiplication
          ! into rtempVecFine for y and xi.
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(1),rtempVecFine%RvectorBlock(1),dscalePrim*p_DaPrim(icol),1.0_DP)
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(2),rtempVecFine%RvectorBlock(2),dscalePrim*p_DaPrim(icol),1.0_DP)
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(6),rtempVecFine%RvectorBlock(6),dscalePrim*p_DaPrim(icol),1.0_DP)
              
        end do
        
        ! Dual space: lambda/p
        do icol = p_KldDual(irow),p_KldDual(irow+1)-1
        
          ! Try to get the vector from the vector pool. Saves time.
          call sptivec_getVectorFromPool(raccessPool,p_KcolDual(icol),p_rx)
          if (.not. associated(p_rx)) then
            ! No, we have to fetch/calculate the vector.
            !
            ! Get a buffer where to save it.
            call sptivec_getFreeBufferFromPool (raccessPool,p_KcolDual(icol),p_rx)
            
            ! Read the source vector. The column of the matrix specifies
            ! the timestep.
            call getSpatialVector (rprojHier,rcoarseVector,p_KcolDual(icol),p_rx,&
                rtempVecCoarse,rtempVecFineScalar,ispacelevelcoarse,ispacelevelfine)
          end if
          
          ! DEBUG!!!
          call lsysbl_getbase_double (p_rx,p_DdataCoarse)
          
          ! Now, rx is the time vector at timestep icol. Weighted multiplication
          ! into rtempVecFine for y and xi.
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(4),rtempVecFine%RvectorBlock(4),dscaleDual*p_DaDual(icol),1.0_DP)
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(5),rtempVecFine%RvectorBlock(5),dscaleDual*p_DaDual(icol),1.0_DP)
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(3),rtempVecFine%RvectorBlock(3),dscaleDual*p_DaDual(icol),1.0_DP)
              
        end do

      case default
      
        call output_line ("Equation not supported.")
        call sys_halt()

      end select

      ! Vector finished.
      call sptivec_setTimestepData (rfineVector, irow, rtempVecFine)

    end do
    
    ! Release the buffer.
    call sptivec_releaseAccessPool(raccessPool)
    
    call lsyssc_releaseVector (rtempVecFineScalar)
    
  end subroutine
     
  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_performEmbedding (rprojHier,ilevel,bdualToPrimal,rsrcVector, &
      rdestVector,rtempVector)
  
!<description>
  ! Embeds a solution vector in another space.
  ! Primal vectors are interpolated from the time-primal space to the 
  ! time-dual space and dual vectors are interpolated from the time-dual space
  ! to the time primal space.
  ! If bdualToPrimal=.true., the embedding is the other way around.
!</description>

!<input>
  ! A space/time interlevel projection structure that configures the 
  ! prolongation/restriction in space/time.
  type(t_sptiProjHierarchy), intent(IN) :: rprojHier

  ! Level corresponding to the solution vectors.
  integer, intent(in) :: ilevel

  ! Source vector to be projected.
  type(t_spacetimeVector), intent(INOUT) :: rsrcVector
  
  ! Type of embedding.
  ! =true: For primal vectors: interpolate from primal to dual time space.
  !        For dual vectors: interpolate from dual to primal time space.
  ! =false: For primal vectors: interpolate from dual to primal time space.
  !         For dual vectors: interpolate from primal to dual time space.
  logical, intent(in) :: bdualToPrimal
!</input>

!<inputoutput>
  ! Temporary space-vector, specifying the discretisation and vector shape.
  type(t_vectorBlock), intent(INOUT) :: rtempVector
!</inputoutput>

!<output>
  ! Destination vector that receives the projected solution.
  type(t_spacetimeVector), intent(INOUT) :: rdestVector
!</output>

!</subroutine>

    ! We have to multiply with the embedding matrices.
    ! The primal vectors have to be multiplied with E^T, the dual with E^-T.
    ! If bdualToPrimal=true, the primal have to be multiplied with E^T and
    ! the dual with E.

    if (.not. bdualToPrimal) then
      call sptipr_timeMatVec (&
          rprojHier%p_RembedMatPrimalInDualT(ilevel),&
          rprojHier%p_RembedMatDualInPrimalT(ilevel),3,&
          rprojHier%p_rphysics,rsrcVector,rdestVector,1.0_DP,0.0_DP,rtempVector)
    else
      call sptipr_timeMatVec (&
          rprojHier%p_RembedMatDualInPrimalT(ilevel),&
          rprojHier%p_RembedMatPrimalInDualT(ilevel),3,&
          rprojHier%p_rphysics,rsrcVector,rdestVector,1.0_DP,0.0_DP,rtempVector)
    end if

!    call sptivec_copyVector (rsrcVector,rdestVector)
!
!    if (.not. bdualToPrimal) then
!      call sptipr_timeMatVec (&
!          rprojHier%p_RembedMatDualInPrimalT(ilevel),&
!          rprojHier%p_RembedMatDualInPrimalT(ilevel),1,&
!          rprojHier%p_rphysics,rsrcVector,rdestVector,1.0_DP,0.0_DP,rtempVector)
!      call sptipr_timeBackwardSolve (&
!          rprojHier%p_RembedMatDualInPrimalT(ilevel),&
!          rprojHier%p_RembedMatDualInPrimalT(ilevel),2,&
!          rprojHier%p_rphysics,rdestVector,rtempVector)
!    else
!      call sptipr_timeMatVec (&
!          rprojHier%p_RembedMatDualInPrimalT(ilevel),&
!          rprojHier%p_RembedMatDualInPrimalT(ilevel),2,&
!          rprojHier%p_rphysics,rsrcVector,rdestVector,1.0_DP,0.0_DP,rtempVector)
!      call sptipr_timeBackwardSolve (&
!          rprojHier%p_RembedMatDualInPrimalT(ilevel),&
!          rprojHier%p_RembedMatDualInPrimalT(ilevel),1,&
!          rprojHier%p_rphysics,rdestVector,rtempVector)
!    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_timeMatVec (rmatrixPrim,rmatrixDual,cprimdual,rphysics,rx,ry,dcx,dcy,rtempVector)
  
!<description>
  ! Performs a matrix-vector multiplication in time with a time-weighting
  ! matrix. rmatrix defines a format-9 matrix containing time-weights.
  ! The routine computes ry = dcx*A*rx + dcy*ry in time.
  ! rmatrixPrimal is used for the primal space, rmatrixDual for the dual space.
!</description>

!<input>
  ! Matrix with weights in time for the primal space.
  type(t_matrixScalar), intent(in) :: rmatrixPrim

  ! Matrix with weights in time for the dual space.
  type(t_matrixScalar), intent(in) :: rmatrixDual
  
  ! Underlying physics of the vectors
  type(t_physics), intent(in) :: rphysics

  ! Determins whether multiply the primal and/or dual part.
  ! =1: Multiply the primal part.
  ! =2: Multiply the dual part.
  ! =3: Multiply the both, primal and dual part.
  integer, intent(in) :: cprimdual

  ! Source vector rx
  type(t_spacetimeVector), intent(INOUT) :: rx
  
  ! Weight for rx
  real(DP), intent(in) :: dcx
  
  ! Weight for ry.
  real(DP), intent(in) :: dcy
!</input>

!<inputoutput>
  ! Temporary vector.
  type(t_vectorBlock), intent(INOUT) :: rtempVector
!</inputoutput>

!<output>
  ! rhs and destination vector.
  type(t_spacetimeVector), intent(INOUT) :: ry
!</output>

!</subroutine>

    ! Local variables
    type(t_spaceTimeVectorAccess) :: raccessPool
    type(t_vectorBlock), pointer :: p_rx
    integer :: irow, icol
    real(DP), dimension(:), pointer :: p_DaPrim,p_DaDual,p_DdataDest,p_DdataSrc
    integer, dimension(:), pointer :: p_KcolPrim, p_KldPrim, p_KcolDual, p_KldDual
    real(DP) :: dscalePrim, dscaleDual

    ! DEBUG!!!
    !call sptivec_setConstant(rfineVector,1.0_DP)

    ! We have to multiply with the embedding matrices.
    ! The primal vectors have to be multiplied with E, the dual with E^T.
    ! If bdualToPrimal=true, the primal have to be multiplied with E^T and
    ! the dual with E.

    ! Get the matrix data arrays.
    call lsyssc_getbase_double (rmatrixPrim,p_DaPrim)
    call lsyssc_getbase_double (rmatrixDual,p_DaDual)

    call lsyssc_getbase_Kcol (rmatrixPrim,p_KcolPrim)
    call lsyssc_getbase_Kld (rmatrixPrim,p_KldPrim)

    call lsyssc_getbase_Kcol (rmatrixDual,p_KcolDual)
    call lsyssc_getbase_Kld (rmatrixDual,p_KldDual)

    ! Allocate temp memory
    call sptivec_createAccessPool (rx,raccessPool,&
        2*max(p_KldPrim(3)-p_KldPrim(2),p_KldDual(3)-p_KldDual(3)))
    
    dscalePrim = rmatrixPrim%dscaleFactor
    dscaleDual = rmatrixDual%dscaleFactor
    
    ! Scale the destination.
    call sptivec_scaleVector (ry,dcy)
    
    ! DEBUG!!!
    call lsysbl_getbase_double (rtempVector,p_DdataDest)
    
    ! Apply the multiplication.
    ! The rows in the matrix correspond to the time fine mesh, the columns
    ! to the time coarse mesh.
    do irow = 1,rmatrixPrim%NEQ
    
      ! Get the destination
      call sptivec_getTimestepData(ry,irow,rtempVector)
      
      select case (rphysics%cequation)
      case (0,2)
        ! Heat equation

        if (iand(cprimdual,1) .ne. 0) then
          ! Primal space: y,xi
          do icol = p_KldPrim(irow),p_KldPrim(irow+1)-1
          
            ! Get the fine grid vector using the vector pool as buffer. Saves time.
            call sptivec_getVectorFromPool(raccessPool,p_KcolPrim(icol),p_rx)
            
            call lsysbl_getbase_double (p_rx,p_DdataSrc)
            
            ! Now, rx is the time vector at timestep icol. Weighted multiplication
            ! into rtempVecFine for y and xi.
            call lsyssc_vectorLinearComb(&
                p_rx%RvectorBlock(1),rtempVector%RvectorBlock(1),&
                dscalePrim*p_DaPrim(icol)*dcx,1.0_DP)
                
          end do
        end if

        if (iand(cprimdual,2) .ne. 0) then        
          ! Dual space: lambda/p
          do icol = p_KldDual(irow),p_KldDual(irow+1)-1
          
            ! Try to get the vector from the vector pool. Saves time.
            call sptivec_getVectorFromPool(raccessPool,p_KcolDual(icol),p_rx)
            
            call lsysbl_getbase_double (p_rx,p_DdataSrc)
            
            ! Now, rx is the time vector at timestep icol. Weighted multiplication
            ! into rtempVecFine for y and xi.
            call lsyssc_vectorLinearComb(&
                p_rx%RvectorBlock(2),rtempVector%RvectorBlock(2),&
                dscaleDual*p_DaDual(icol)*dcx,1.0_DP)
                
          end do
        end if
        
      case (1)
        ! Stokes equation

        if (iand(cprimdual,1) .ne. 0) then
          ! Primal space: y,xi
          do icol = p_KldPrim(irow),p_KldPrim(irow+1)-1
          
            ! Get the fine grid vector using the vector pool as buffer. Saves time.
            call sptivec_getVectorFromPool(raccessPool,p_KcolPrim(icol),p_rx)
            
            call lsysbl_getbase_double (p_rx,p_DdataDest)
            
            ! Now, rx is the time vector at timestep icol. Weighted multiplication
            ! into rtempVecFine for y and xi.
            call lsyssc_vectorLinearComb(&
                p_rx%RvectorBlock(1),rtempVector%RvectorBlock(1),&
                dscalePrim*p_DaPrim(icol)*dcx,1.0_DP)
            call lsyssc_vectorLinearComb(&
                p_rx%RvectorBlock(2),rtempVector%RvectorBlock(2),&
                dscalePrim*p_DaPrim(icol)*dcx,1.0_DP)
            call lsyssc_vectorLinearComb(&
                p_rx%RvectorBlock(6),rtempVector%RvectorBlock(6),&
                dscalePrim*p_DaPrim(icol)* dcx,1.0_DP)
                
          end do
        end if
        
        if (iand(cprimdual,2) .ne. 0) then
          ! Dual space: lambda/p
          do icol = p_KldDual(irow),p_KldDual(irow+1)-1
          
            ! Try to get the vector from the vector pool. Saves time.
            call sptivec_getVectorFromPool(raccessPool,p_KcolDual(icol),p_rx)
            
            call lsysbl_getbase_double (p_rx,p_DdataDest)
            
            ! Now, rx is the time vector at timestep icol. Weighted multiplication
            ! into rtempVecFine for y and xi.
            call lsyssc_vectorLinearComb(&
                p_rx%RvectorBlock(4),rtempVector%RvectorBlock(4),&
                dscaleDual*p_DaDual(icol)*dcx,1.0_DP)
            call lsyssc_vectorLinearComb(&
                p_rx%RvectorBlock(5),rtempVector%RvectorBlock(5),&
                dscaleDual*p_DaDual(icol)*dcx,1.0_DP)
            call lsyssc_vectorLinearComb(&
                p_rx%RvectorBlock(3),rtempVector%RvectorBlock(3),&
                dscaleDual*p_DaDual(icol)*dcx,1.0_DP)
                
          end do
        end if

      case default
      
        call output_line ("Equation not supported.")
        call sys_halt()
      end select
      
      ! Vector finished.
      call sptivec_setTimestepData (ry, irow, rtempVector)

    end do
    
    ! Release the buffer.
    call sptivec_releaseAccessPool(raccessPool)
    
    ! DEBUG!!!
    !call sptivec_saveToFileSequence (rcoarseVector,"(""./ns/coarse.txt."",I5.5)",.true.,&
    !    rtempVecCoarse)
    !call sys_halt()
        
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_timeBackwardSolve (rmatrixPrim,rmatrixDual,cprimdual,rphysics,ry,rtempVector)
  
!<description>
  ! Performs a backward solve to get $ry(new) = A^{-1} ry$.
  ! Only elements on the upper diagonal of the matrix are respected.
!</description>

!<input>
  ! Matrix with weights in time for the primal space.
  type(t_matrixScalar), intent(in) :: rmatrixPrim

  ! Matrix with weights in time for the dual space.
  type(t_matrixScalar), intent(in) :: rmatrixDual
  
  ! Underlying physics of the vectors
  type(t_physics), intent(in) :: rphysics
  
  ! Determins whether to solve for the primal and/or dual part.
  ! =1: Solve for primal part.
  ! =2: Solve for dual part.
  ! =3: Solve for both, primal and dual part.
  integer, intent(in) :: cprimdual
!</input>

!<inputoutput>
  ! Temporary vector.
  type(t_vectorBlock), intent(INOUT) :: rtempVector
!</inputoutput>

!<output>
  ! Source and destination vector.
  type(t_spacetimeVector), intent(INOUT) :: ry
!</output>

!</subroutine>

    ! Local variables
    type(t_spaceTimeVectorAccess) :: raccessPool
    type(t_vectorBlock), pointer :: p_rx
    integer :: irow, icol
    real(DP), dimension(:), pointer :: p_DaPrim,p_DaDual,p_DdataDest,p_DdataSrc
    integer, dimension(:), pointer :: p_KcolPrim, p_KldPrim, p_KcolDual, p_KldDual
    integer, dimension(:), pointer :: p_KdiagonalPrim, p_KdiagonalDual
    real(DP) :: dscalePrim, dscaleDual

    ! DEBUG!!!
    !call sptivec_setConstant(rfineVector,1.0_DP)

    ! We have to multiply with the embedding matrices.
    ! The primal vectors have to be multiplied with E, the dual with E^T.
    ! If bdualToPrimal=true, the primal have to be multiplied with E^T and
    ! the dual with E.

    ! Get the matrix data arrays.
    call lsyssc_getbase_double (rmatrixPrim,p_DaPrim)
    call lsyssc_getbase_double (rmatrixDual,p_DaDual)

    call lsyssc_getbase_Kcol (rmatrixPrim,p_KcolPrim)
    call lsyssc_getbase_Kld (rmatrixPrim,p_KldPrim)
    call lsyssc_getbase_Kdiagonal (rmatrixPrim,p_KdiagonalPrim)

    call lsyssc_getbase_Kcol (rmatrixDual,p_KcolDual)
    call lsyssc_getbase_Kld (rmatrixDual,p_KldDual)
    call lsyssc_getbase_Kdiagonal (rmatrixDual,p_KdiagonalDual)

    ! Allocate temp memory
    call sptivec_createAccessPool (ry,raccessPool,&
        2*max(p_KldPrim(3)-p_KldPrim(2),p_KldDual(3)-p_KldDual(3)))
    
    dscalePrim = rmatrixPrim%dscaleFactor
    dscaleDual = rmatrixDual%dscaleFactor
    
    ! DEBUG!!!
    call lsysbl_getbase_double (rtempVector,p_DdataDest)
    
    ! Apply the multiplication.
    ! The rows in the matrix correspond to the time fine mesh, the columns
    ! to the time coarse mesh.
    do irow = rmatrixPrim%NEQ,1,-1
    
      ! Get the source/destination
      call sptivec_getTimestepData(ry,irow,rtempVector)
      
      select case (rphysics%cequation)
      case (0,2)
        ! Heat equation

        if (iand(cprimdual,1) .ne. 0) then
          ! Primal space: y,xi
          do icol = p_KdiagonalPrim(irow)+1,p_KldPrim(irow+1)-1
          
            ! Get the fine grid vector using the vector pool as buffer. Saves time.
            call sptivec_getVectorFromPool(raccessPool,p_KcolPrim(icol),p_rx)
            
            call lsysbl_getbase_double (p_rx,p_DdataSrc)
            
            ! Subtract
            call lsyssc_vectorLinearComb(&
                p_rx%RvectorBlock(1),rtempVector%RvectorBlock(1),&
                -dscalePrim*p_DaPrim(icol),1.0_DP)
                
          end do
          
          ! Divide by the diagonal.
          icol = p_KdiagonalPrim(irow)
          call lsyssc_scaleVector(rtempVector%RvectorBlock(1),1.0_DP/(dscalePrim*p_DaPrim(icol)))
        end if
        
        if (iand(cprimdual,2) .ne. 0) then
          ! Dual space: lambda/p
          do icol = p_KdiagonalDual(irow)+1,p_KldDual(irow+1)-1
          
            ! Get the fine grid vector using the vector pool as buffer. Saves time.
            call sptivec_getVectorFromPool(raccessPool,p_KcolDual(icol),p_rx)
            
            call lsysbl_getbase_double (p_rx,p_DdataSrc)
            
            ! Subtract
            call lsyssc_vectorLinearComb(&
                p_rx%RvectorBlock(2),rtempVector%RvectorBlock(2),&
                -dscaleDual*p_DaDual(icol),1.0_DP)
                
          end do
          
          ! Divide by the diagonal.
          icol = p_KdiagonalDual(irow)
          call lsyssc_scaleVector(rtempVector%RvectorBlock(2),1.0_DP/(dscaleDual*p_DaDual(icol)))
        end if
        
      case (1)
        ! Stokes equation

        if (iand(cprimdual,1) .ne. 0) then
          ! Primal space: y,xi
          do icol = p_KdiagonalPrim(irow)+1,p_KldPrim(irow+1)-1
          
            ! Get the fine grid vector using the vector pool as buffer. Saves time.
            call sptivec_getVectorFromPool(raccessPool,p_KcolPrim(icol),p_rx)
            
            call lsysbl_getbase_double (p_rx,p_DdataSrc)
            
            ! Subtract
            call lsyssc_vectorLinearComb(&
                p_rx%RvectorBlock(1),rtempVector%RvectorBlock(1),&
                -dscalePrim*p_DaPrim(icol),1.0_DP)
            call lsyssc_vectorLinearComb(&
                p_rx%RvectorBlock(2),rtempVector%RvectorBlock(2),&
                -dscalePrim*p_DaPrim(icol),1.0_DP)
            call lsyssc_vectorLinearComb(&
                p_rx%RvectorBlock(6),rtempVector%RvectorBlock(6),&
                -dscalePrim*p_DaPrim(icol),1.0_DP)
                
          end do
          
          ! Divide by the diagonal.
          icol = p_KdiagonalPrim(irow)
          call lsyssc_scaleVector(rtempVector%RvectorBlock(1),1.0_DP/(dscalePrim*p_DaPrim(icol)))
          call lsyssc_scaleVector(rtempVector%RvectorBlock(2),1.0_DP/(dscalePrim*p_DaPrim(icol)))
          call lsyssc_scaleVector(rtempVector%RvectorBlock(6),1.0_DP/(dscalePrim*p_DaPrim(icol)))
        end if
        
        if (iand(cprimdual,2) .ne. 0) then

          ! Dual space: lambda/p
          do icol = p_KdiagonalDual(irow)+1,p_KldDual(irow+1)-1
          
            ! Get the fine grid vector using the vector pool as buffer. Saves time.
            call sptivec_getVectorFromPool(raccessPool,p_KcolDual(icol),p_rx)
            
            call lsysbl_getbase_double (p_rx,p_DdataSrc)
            
            ! Subtract
            call lsyssc_vectorLinearComb(&
                p_rx%RvectorBlock(4),rtempVector%RvectorBlock(4),&
                -dscaleDual*p_DaDual(icol),1.0_DP)
            call lsyssc_vectorLinearComb(&
                p_rx%RvectorBlock(5),rtempVector%RvectorBlock(5),&
                -dscaleDual*p_DaDual(icol),1.0_DP)
            call lsyssc_vectorLinearComb(&
                p_rx%RvectorBlock(3),rtempVector%RvectorBlock(3),&
                -dscaleDual*p_DaDual(icol),1.0_DP)
                
          end do
          
          ! Divide by the diagonal.
          icol = p_KdiagonalDual(irow)
          call lsyssc_scaleVector(rtempVector%RvectorBlock(4),1.0_DP/(dscaleDual*p_DaDual(icol)))
          call lsyssc_scaleVector(rtempVector%RvectorBlock(5),1.0_DP/(dscaleDual*p_DaDual(icol)))
          call lsyssc_scaleVector(rtempVector%RvectorBlock(3),1.0_DP/(dscaleDual*p_DaDual(icol)))
        end if

      case default
      
        call output_line ("Equation not supported.")
        call sys_halt()
      end select
      
      ! Vector finished.
      call sptivec_setTimestepData (ry, irow, rtempVector)
      call sptivec_invalidateVecInPool (raccessPool,irow)

    end do
    
    ! Release the buffer.
    call sptivec_releaseAccessPool(raccessPool)
    
    ! DEBUG!!!
    !call sptivec_saveToFileSequence (rcoarseVector,"(""./ns/coarse.txt."",I5.5)",.true.,&
    !    rtempVecCoarse)
    !call sys_halt()
        
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_performRestriction (rprojHier,ilevelfine,rcoarseVector, &
      rfineVector,rtempVecCoarse,rtempVecFine,&
      rtempVecSpaceTimeCoarse,rtempVecSpaceTimeFine)
  
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

  ! Temporary space-time-vector, specifying the discretisation and vector shape
  ! on the coarse grid.
  type(t_spacetimeVector), intent(INOUT) :: rtempVecSpaceTimeCoarse

  ! Temporary space--timevector, specifying the discretisation and vector shape
  ! on the fine grid.
  type(t_spacetimeVector), intent(INOUT) :: rtempVecSpaceTimeFine
!</inputoutput>

!<output>
  ! Coarse grid vector
  type(t_spacetimeVector), intent(INOUT) :: rcoarseVector
!</output>

!</subroutine>

    ! Local variables
    type(t_spaceTimeVectorAccess) :: raccessPool
    integer :: ispacelevelfine,ispacelevelcoarse,itimelevelfine,itimelevelcoarse
    type(t_vectorScalar) :: rtempVecFineScalar
    type(t_vectorBlock), pointer :: p_rx
    integer :: irow, icol
    type(t_matrixScalar), pointer :: p_rrestMatrixPrim,p_rrestMatrixDual
    real(DP), dimension(:), pointer :: p_DaPrim,p_DaDual,p_DdataCoarse,p_DdataFine
    integer, dimension(:), pointer :: p_KcolPrim, p_KldPrim, p_KcolDual, p_KldDual
    real(DP) :: dscalePrim, dscaleDual

    ! DEBUG!!!
    !call sptivec_setConstant(rfineVector,1.0_DP)

    ! Get the space-time coarse and fine discretisations.
    call sth_getLevel (rprojHier%p_rspaceTimeHierarchy,ilevelfine-1,&
      ispaceLevel=ispacelevelcoarse,itimeLevel=itimelevelcoarse)
    call sth_getLevel (rprojHier%p_rspaceTimeHierarchy,ilevelfine,&
      ispaceLevel=ispacelevelfine,itimeLevel=itimelevelfine)
    
    ! We need a scalar representation of the temp vector
    call lsysbl_createScalarFromVec (rtempVecFine,rtempVecFineScalar)
    
    ! Allocate temp memory
    select case (rprojHier%ctimeProjection)
    case (0)
      call sptivec_createAccessPool (rfineVector,raccessPool,3)
    case (1)
      call sptivec_createAccessPool (rfineVector,raccessPool,3)
    case (2,3,6,7,8)
      call sptivec_createAccessPool (rfineVector,raccessPool,7)
    case (4)
      call sptivec_createAccessPool (rfineVector,raccessPool,10)
    case default
      call sptivec_createAccessPool (rfineVector,raccessPool,10)
    end select
    
    ! Prolongation means, we multiply with the prolongation matrix.
    ! y and xi have to be multiplied with the primal prol. matrix,
    ! lambda and p with the dual. Note that for implicit Euler,
    ! the primal and dual prolongation matrix is the same, while
    ! for CN, y and xi are at the same points in time as well as
    ! lambda and p: p is between y and thus at the same points
    ! in time as lambda, not as y!
    !
    ! Get the matrix.
    p_rrestMatrixPrim => rprojHier%p_rrestrictionMatPrimal(ilevelfine-1)
    p_rrestMatrixDual => rprojHier%p_rrestrictionMatDual(ilevelfine-1)
    
    call lsyssc_getbase_double (p_rrestMatrixPrim,p_DaPrim)
    call lsyssc_getbase_double (p_rrestMatrixDual,p_DaDual)

    call lsyssc_getbase_Kcol (p_rrestMatrixPrim,p_KcolPrim)
    call lsyssc_getbase_Kld (p_rrestMatrixPrim,p_KldPrim)

    call lsyssc_getbase_Kcol (p_rrestMatrixDual,p_KcolDual)
    call lsyssc_getbase_Kld (p_rrestMatrixDual,p_KldDual)
    
    dscalePrim = p_rrestMatrixPrim%dscaleFactor
    dscaleDual = p_rrestMatrixDual%dscaleFactor
    
    ! DEBUG!!!
    call lsysbl_getbase_double (rtempVecFine,p_DdataFine)
    
    select case (rprojHier%ctimeProjection)
    case (3)
      ! Apply the embedding since the space of the vector has to be corrected.
      call sptivec_copyVector (rfineVector,rtempVecSpaceTimeFine)
      call sptipr_performEmbedding (rprojHier,ilevelFine,.false.,rtempVecSpaceTimeFine, &
          rfineVector,rtempVecFine)
    end select
    
    ! Apply the multiplication.
    ! The rows in the matrix correspond to the time fine mesh, the columns
    ! to the time coarse mesh.
    do irow = 1,p_rrestMatrixPrim%NEQ
    
      ! Clear the destination
      call lsysbl_clearVector (rtempVecFine)
      
      select case (rprojHier%p_rphysics%cequation)
      case (0,2)
        ! Heat equation

        ! Primal space: y,xi
        do icol = p_KldPrim(irow),p_KldPrim(irow+1)-1
        
          ! Get the fine grid vector using the vector pool as buffer. Saves time.
          call sptivec_getVectorFromPool(raccessPool,p_KcolPrim(icol),p_rx)
          
          call lsysbl_getbase_double (p_rx,p_DdataCoarse)
          
          ! Now, rx is the time vector at timestep icol. Weighted multiplication
          ! into rtempVecFine for y and xi.
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(1),rtempVecFine%RvectorBlock(1),dscalePrim*p_DaPrim(icol),1.0_DP)
              
        end do
        
        ! Dual space: lambda/p
        do icol = p_KldDual(irow),p_KldDual(irow+1)-1
        
          ! Try to get the vector from the vector pool. Saves time.
          call sptivec_getVectorFromPool(raccessPool,p_KcolDual(icol),p_rx)
          
          call lsysbl_getbase_double (p_rx,p_DdataCoarse)
          
          ! Now, rx is the time vector at timestep icol. Weighted multiplication
          ! into rtempVecFine for y and xi.
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(2),rtempVecFine%RvectorBlock(2),dscaleDual*p_DaDual(icol),1.0_DP)
              
        end do
        
      case (1)
        ! Stokes equation

        ! Primal space: y,xi
        do icol = p_KldPrim(irow),p_KldPrim(irow+1)-1
        
          ! Get the fine grid vector using the vector pool as buffer. Saves time.
          call sptivec_getVectorFromPool(raccessPool,p_KcolPrim(icol),p_rx)
          
          call lsysbl_getbase_double (p_rx,p_DdataCoarse)
          
          ! Now, rx is the time vector at timestep icol. Weighted multiplication
          ! into rtempVecFine for y and xi.
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(1),rtempVecFine%RvectorBlock(1),dscalePrim*p_DaPrim(icol),1.0_DP)
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(2),rtempVecFine%RvectorBlock(2),dscalePrim*p_DaPrim(icol),1.0_DP)
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(6),rtempVecFine%RvectorBlock(6),dscalePrim*p_DaPrim(icol),1.0_DP)
              
        end do
        
        ! Dual space: lambda/p
        do icol = p_KldDual(irow),p_KldDual(irow+1)-1
        
          ! Try to get the vector from the vector pool. Saves time.
          call sptivec_getVectorFromPool(raccessPool,p_KcolDual(icol),p_rx)
          
          call lsysbl_getbase_double (p_rx,p_DdataCoarse)
          
          ! Now, rx is the time vector at timestep icol. Weighted multiplication
          ! into rtempVecFine for y and xi.
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(4),rtempVecFine%RvectorBlock(4),dscaleDual*p_DaDual(icol),1.0_DP)
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(5),rtempVecFine%RvectorBlock(5),dscaleDual*p_DaDual(icol),1.0_DP)
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(3),rtempVecFine%RvectorBlock(3),dscaleDual*p_DaDual(icol),1.0_DP)
              
        end do

      case default
      
        call output_line ("Equation not supported.")
        call sys_halt()
      end select
      
      ! Vector finished.
      call setSpatialVector (rprojHier,rtempVecFine,rcoarseVector,irow,&
          rtempVecCoarse,rtempVecFineScalar,ispacelevelcoarse,ispacelevelfine)

    end do

    select case (rprojHier%ctimeProjection)
    case (3)
      ! Apply the reverse embedding.
      call sptivec_copyVector (rcoarseVector,rtempVecSpaceTimeCoarse)
      call sptipr_performEmbedding (rprojHier,ilevelFine-1,.true.,rtempVecSpaceTimeCoarse, &
          rcoarseVector,rtempVecCoarse)
    end select
    
    ! Release the buffer.
    call sptivec_releaseAccessPool(raccessPool)
    
    call lsyssc_releaseVector (rtempVecFineScalar)

    ! DEBUG!!!
    !call sptivec_saveToFileSequence (rcoarseVector,"(""./ns/coarse.txt."",I5.5)",.true.,&
    !    rtempVecCoarse)
    !call sys_halt()
        
  contains

    subroutine setSpatialVector (rprojHier,rx,rcoarseVector,iindex,&
        rtempVecCoarse,rtempVecFineScalar,ispacelevelcoarse,ispacelevelfine)
    
    ! Saves the spatial subvector iindex from rx to rfineVector.
    ! If necessary, the vector is restricted in space.
    
    ! A space/time interlevel projection structure that configures the 
    ! prolongation/restriction in space/time.
    type(t_sptiProjHierarchy), intent(in) :: rprojHier

    ! Space-time destination vector
    type(t_spaceTimeVector), intent(inout) :: rcoarseVector
    
    ! Index of the subvector
    integer, intent(in) :: iindex
    
    ! Space vector output.
    type(t_vectorBlock), intent(inout) :: rx
    
    ! Temp vectors
    type(t_vectorBlock), intent(inout) :: rtempVecCoarse
    type(t_vectorScalar), intent(inout) :: rtempVecFineScalar
    
    ! Level of the coarse and (destination) fine vector.
    integer, intent(in) :: ispacelevelcoarse
    integer, intent(in) :: ispacelevelfine
  
      if (ispacelevelcoarse .ne. ispacelevelfine) then
        ! Space + time
        call mlprj_performRestriction (&
            rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
            rtempVecCoarse,rx,rtempVecFineScalar)
        call sptivec_setTimestepData (rcoarseVector, iindex, rtempVecCoarse)
      else
        ! Only time
        call sptivec_setTimestepData (rcoarseVector, iindex, rx)
      end if

    end subroutine
    
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

    ! Local variables
    type(t_spaceTimeVectorAccess) :: raccessPool
    integer :: ispacelevelfine,ispacelevelcoarse,itimelevelfine,itimelevelcoarse
    type(t_vectorScalar) :: rtempVecFineScalar
    type(t_vectorBlock), pointer :: p_rx
    integer :: irow, icol
    type(t_matrixScalar), pointer :: p_rinterpMatrixPrim,p_rinterpMatrixDual
    real(DP), dimension(:), pointer :: p_DaPrim,p_DaDual,p_DdataCoarse,p_DdataFine
    integer, dimension(:), pointer :: p_KcolPrim, p_KldPrim, p_KcolDual, p_KldDual
    real(DP) :: dscalePrim, dscaleDual

    ! DEBUG!!!
    !call sptivec_setConstant(rfineVector,1.0_DP)

    ! Get the space-time coarse and fine discretisations.
    call sth_getLevel (rprojHier%p_rspaceTimeHierarchy,ilevelfine-1,&
      ispaceLevel=ispacelevelcoarse,itimeLevel=itimelevelcoarse)
    call sth_getLevel (rprojHier%p_rspaceTimeHierarchy,ilevelfine,&
      ispaceLevel=ispacelevelfine,itimeLevel=itimelevelfine)
    
    ! We need a scalar representation of the temp vector
    call lsysbl_createScalarFromVec (rtempVecFine,rtempVecFineScalar)
    
    ! Allocate temp memory
    select case (rprojHier%ctimeProjection)
    case (0)
      call sptivec_createAccessPool (rfineVector,raccessPool,3)
    case (1)
      call sptivec_createAccessPool (rfineVector,raccessPool,3)
    case (2,3,6,7,8)
      call sptivec_createAccessPool (rfineVector,raccessPool,7)
    case (4)
      call sptivec_createAccessPool (rfineVector,raccessPool,10)
    case default
      call sptivec_createAccessPool (rfineVector,raccessPool,10)
    end select
    
    ! Prolongation means, we multiply with the prolongation matrix.
    ! y and xi have to be multiplied with the primal prol. matrix,
    ! lambda and p with the dual. Note that for implicit Euler,
    ! the primal and dual prolongation matrix is the same, while
    ! for CN, y and xi are at the same points in time as well as
    ! lambda and p: p is between y and thus at the same points
    ! in time as lambda, not as y!
    !
    ! Get the matrix.
    p_rinterpMatrixPrim => rprojHier%p_RinterpolationMatPrimal(ilevelfine-1)
    p_rinterpMatrixDual => rprojHier%p_RinterpolationMatDual(ilevelfine-1)
    
    call lsyssc_getbase_double (p_rinterpMatrixPrim,p_DaPrim)
    call lsyssc_getbase_double (p_rinterpMatrixDual,p_DaDual)

    call lsyssc_getbase_Kcol (p_rinterpMatrixPrim,p_KcolPrim)
    call lsyssc_getbase_Kld (p_rinterpMatrixPrim,p_KldPrim)

    call lsyssc_getbase_Kcol (p_rinterpMatrixDual,p_KcolDual)
    call lsyssc_getbase_Kld (p_rinterpMatrixDual,p_KldDual)
    
    dscalePrim = p_rinterpMatrixPrim%dscaleFactor
    dscaleDual = p_rinterpMatrixDual%dscaleFactor
    
    ! DEBUG!!!
    call lsysbl_getbase_double (rtempVecFine,p_DdataFine)
    
    ! Apply the multiplication.
    ! The rows in the matrix correspond to the time fine mesh, the columns
    ! to the time coarse mesh.
    do irow = 1,p_rinterpMatrixPrim%NEQ
    
      ! Clear the destination
      call lsysbl_clearVector (rtempVecFine)
      
      select case (rprojHier%p_rphysics%cequation)
      case (0,2)
        ! Heat equation

        ! Primal space: y,xi
        do icol = p_KldPrim(irow),p_KldPrim(irow+1)-1
        
          ! Get the fine grid vector using the vector pool as buffer. Saves time.
          call sptivec_getVectorFromPool(raccessPool,p_KcolPrim(icol),p_rx)
          
          call lsysbl_getbase_double (p_rx,p_DdataCoarse)
          
          ! Now, rx is the time vector at timestep icol. Weighted multiplication
          ! into rtempVecFine for y and xi.
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(1),rtempVecFine%RvectorBlock(1),dscalePrim*p_DaPrim(icol),1.0_DP)
              
        end do
        
        ! Dual space: lambda/p
        do icol = p_KldDual(irow),p_KldDual(irow+1)-1
        
          ! Try to get the vector from the vector pool. Saves time.
          call sptivec_getVectorFromPool(raccessPool,p_KcolDual(icol),p_rx)
          
          call lsysbl_getbase_double (p_rx,p_DdataCoarse)
          
          ! Now, rx is the time vector at timestep icol. Weighted multiplication
          ! into rtempVecFine for y and xi.
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(2),rtempVecFine%RvectorBlock(2),dscaleDual*p_DaDual(icol),1.0_DP)
              
        end do
        
      case (1)
        ! Stokes equation

        ! Primal space: y,xi
        do icol = p_KldPrim(irow),p_KldPrim(irow+1)-1
        
          ! Get the fine grid vector using the vector pool as buffer. Saves time.
          call sptivec_getVectorFromPool(raccessPool,p_KcolPrim(icol),p_rx)
          
          call lsysbl_getbase_double (p_rx,p_DdataCoarse)
          
          ! Now, rx is the time vector at timestep icol. Weighted multiplication
          ! into rtempVecFine for y and xi.
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(1),rtempVecFine%RvectorBlock(1),dscalePrim*p_DaPrim(icol),1.0_DP)
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(2),rtempVecFine%RvectorBlock(2),dscalePrim*p_DaPrim(icol),1.0_DP)
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(6),rtempVecFine%RvectorBlock(6),dscalePrim*p_DaPrim(icol),1.0_DP)
              
        end do
        
        ! Dual space: lambda/p
        do icol = p_KldDual(irow),p_KldDual(irow+1)-1
        
          ! Try to get the vector from the vector pool. Saves time.
          call sptivec_getVectorFromPool(raccessPool,p_KcolDual(icol),p_rx)
          
          call lsysbl_getbase_double (p_rx,p_DdataCoarse)
          
          ! Now, rx is the time vector at timestep icol. Weighted multiplication
          ! into rtempVecFine for y and xi.
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(4),rtempVecFine%RvectorBlock(4),dscaleDual*p_DaDual(icol),1.0_DP)
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(5),rtempVecFine%RvectorBlock(5),dscaleDual*p_DaDual(icol),1.0_DP)
          call lsyssc_vectorLinearComb(&
              p_rx%RvectorBlock(3),rtempVecFine%RvectorBlock(3),dscaleDual*p_DaDual(icol),1.0_DP)
              
        end do

      case default
      
        call output_line ("Equation not supported.")
        call sys_halt()
      end select
      
      ! Vector finished.
      call setSpatialVector (rprojHier,rtempVecFine,rcoarseVector,irow,&
          rtempVecCoarse,rtempVecFineScalar,ispacelevelcoarse,ispacelevelfine)

    end do
    
    ! Release the buffer.
    call sptivec_releaseAccessPool(raccessPool)
    
    call lsyssc_releaseVector (rtempVecFineScalar)

    ! DEBUG!!!
    !call sptivec_saveToFileSequence (rcoarseVector,"(""coarse.txt."",I5.5)",.true.,&
    !    rtempVecFine)

  contains

    subroutine setSpatialVector (rprojHier,rx,rcoarseVector,iindex,&
        rtempVecCoarse,rtempVecFineScalar,ispacelevelcoarse,ispacelevelfine)
    
    ! Saves the spatial subvector iindex from rx to rfineVector.
    ! If necessary, the vector is restricted in space.
    
    ! A space/time interlevel projection structure that configures the 
    ! prolongation/restriction in space/time.
    type(t_sptiProjHierarchy), intent(in) :: rprojHier

    ! Space-time destination vector
    type(t_spaceTimeVector), intent(inout) :: rcoarseVector
    
    ! Index of the subvector
    integer, intent(in) :: iindex
    
    ! Space vector output.
    type(t_vectorBlock), intent(inout) :: rx
    
    ! Temp vectors
    type(t_vectorBlock), intent(inout) :: rtempVecCoarse
    type(t_vectorScalar), intent(inout) :: rtempVecFineScalar
    
    ! Level of the coarse and (destination) fine vector.
    integer, intent(in) :: ispacelevelcoarse
    integer, intent(in) :: ispacelevelfine
  
      if (ispacelevelcoarse .ne. ispacelevelfine) then
        ! Space + time
        call mlprj_performInterpolation (&
            rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
            rtempVecCoarse,rx,rtempVecFineScalar)
        call sptivec_setTimestepData (rcoarseVector, iindex, rtempVecCoarse)
      else
        ! Only time
        call sptivec_setTimestepData (rcoarseVector, iindex, rx)
      end if

    end subroutine

  end subroutine

end module
