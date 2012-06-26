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
!# 6.) sptipr_initProjectionBlock
!#     -> Initialises a space-time projection block structure
!#
!# 7.) sptipr_doneProjectionBlock
!#     -> Releases a space-time projection block structure
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
  
  use constantsdiscretisation
  use structuresdiscretisation
  use structuresoptcontrol

  implicit none
  
  private
  
  public :: t_sptiProjHierarchy
  public :: t_sptiProjHierarchyBlock
  public :: sptipr_initProjection
  public :: sptipr_doneProjection
  public :: sptipr_performProlongation
  public :: sptipr_performRestriction
  public :: sptipr_performInterpolation
  public :: sptipr_initProjectionBlock
  public :: sptipr_doneProjectionBlock
  
!<types>

!<typeblock>
  
  ! A hierarchy of space-time projection definitions for a hierarchy
  ! of space-time levels.
  type t_sptiProjHierarchy

    ! Type of projection in time.
    ! =-1: undefined
    ! =0: constant prolongation/restriction
    ! =1: linear prolongation/restriction. Central differences. 
    !     Solutions located at the endpoints of the time interval.
    ! =2: linear prolongation/restriction. Central differences.
    !     Solutions located in time according to the theta scheme.
    !     Full linear restriction with interpolation in time.
    integer :: ctimeProjection = -1
    
    ! List of components in the vector to which this type
    ! of projection is to be applied.
    ! If not associated, the projection is applied to all components.
    integer, dimension(:), pointer :: p_Icomponents => null()
    
    ! Underlying physics
    type(t_settings_physics), pointer :: p_rphysics => null()
    
    ! Pointer to the time coarse mesh.
    type(t_timeDiscretisation), pointer :: p_rtimeCoarseDiscr => null()

    ! A space-time hierarchy that describes the discretisation in space and time
    ! for all levels.
    type(t_spaceTimeHierarchy), pointer :: p_rspaceTimeHierarchy => null()

    ! Pointer to a hierarchy of interlevel projection structures in space.
    type(t_interlevelProjectionHier), pointer :: p_rprojHierarchySpace => null()
    
    ! An array of prolongation matrices.
    ! rprojectionMat(i) defines the weights for the projection
    ! between level i and level i+1.
    type(t_matrixScalar), dimension(:), pointer :: p_RprolongationMat => null()

    ! An array of restriction matrices.
    ! rprojectionMat(i) defines the weights for the projection
    ! between level i+1 and level i.
    type(t_matrixScalar), dimension(:), pointer :: p_RrestrictionMat => null()
    
    ! An array of interpolation matrices.
    ! rprojectionMat(i) defines the weights for the projection
    ! between level i+1 and level i.
    type(t_matrixScalar), dimension(:), pointer :: p_RinterpolationMat => null()

  end type

!</typeblock>

!<typeblock>

  ! A set of space-time projection hierarchy that allows to specify a
  ! separate prolongation/restriction in space and time for every component
  ! of a solution
  type t_sptiProjHierarchyBlock
  
    ! A set of projection hierarchies. Every structure defines a projection
    ! for s set of solution components.
    type(t_sptiProjHierarchy), dimension(:), pointer :: p_RprojHier => null()
  
    ! Number of entries in the list.
    integer :: ncount = 0
  
  end type
  
!</typeblock>

!</types>


contains

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_initProjectionBlock (rprojHierBlock,rspaceTimeHierarchy,&
      rprojHierarchySpace,rphysics,roptcontrol,cspace,ctimeProjection)
  
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
  type(t_settings_physics), intent(in), target :: rphysics

  ! Optimal control parameters
  type(t_settings_optcontrol), intent(in) :: roptcontrol

  ! Type of space, this projection is set up for.
  ! =CCSPACE_PRIMAL  : Primal space, forward in time
  ! =CCSPACE_DUAL    : Dual space, backward in time
  ! =CCSPACE_CONTROL : Control space
  integer, intent(in) :: cspace

  ! Type of projection in time.
  ! =-1: automatic.
  ! =0: constant prolongation/restriction
  ! =1: linear prolongation/restriction. Central differences. 
  !     Solutions located at the endpoints in time.
  ! =2: linear prolongation/restriction. Central differences. 
  !     Dual Solutions located in time according to the theta scheme.
  integer, intent(in), optional :: ctimeProjection
!</input>

!<output>
  ! Space-time projection block structure.
  type(t_sptiProjHierarchyBlock), intent(out) :: rprojHierBlock
!</output>

!</subroutine>

    ! local variables
    integer :: i,ncount
    
    ! Number of components depends on the type of control,...
    ncount = 0

    select case (cspace)
    ! ===============================================
    ! Primal space
    ! ===============================================
    case (CCSPACE_PRIMAL)
      ncount = 1
    
    ! ===============================================
    ! Dual space
    ! ===============================================
    case (CCSPACE_DUAL)
      ncount = 1

    ! ===============================================
    ! Control space
    ! ===============================================
    case (CCSPACE_CONTROL)
      if (roptcontrol%dalphaC .ge. 0.0_DP) ncount = ncount + 1
    end select
    
    rprojHierBlock%ncount = ncount
    allocate(rprojHierBlock%p_RprojHier(ncount))

    ! Initialise all substructures
    do i=1,rprojHierBlock%ncount
      call sptipr_initProjection (rprojHierBlock%p_RprojHier(i),rspaceTimeHierarchy,&
          rprojHierarchySpace,rphysics,roptcontrol,cspace,ctimeProjection,i)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_doneProjectionBlock (rprojHierBlock)
  
!<description>
  ! Cleans up a block of projection structures.
!</description>

!<output>
  ! Space-time projection block structure to be cleaned up
  type(t_sptiProjHierarchyBlock), intent(inout) :: rprojHierBlock
!</output>

!</subroutine>

    ! local variables
    integer :: i
    
    ! Deallocate all substructures
    do i=rprojHierBlock%ncount,1,-1
      call sptipr_doneProjection (rprojHierBlock%p_RprojHier(i))
    end do

    ! Deallocate the structure itself
    deallocate(rprojHierBlock%p_RprojHier)
    rprojHierBlock%ncount = 0

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_initProjection (rprojHier,rspaceTimeHierarchy,&
      rprojHierarchySpace,rphysics,roptcontrol,cspace,ctimeProjection,isubspace)
  
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
  type(t_settings_physics), intent(in), target :: rphysics
  
  ! Optimal control parameters
  type(t_settings_optcontrol), intent(in) :: roptcontrol

  ! Type of space, this projection is set up for.
  ! =CCSPACE_PRIMAL  : Primal space, forward in time
  ! =CCSPACE_DUAL    : Dual space, backward in time
  ! =CCSPACE_CONTROL : Control space
  integer, intent(in) :: cspace

  ! Type of projection in time.
  ! =-1: automatic.
  ! =0: constant prolongation/restriction
  ! =1: linear prolongation/restriction. Central differences. 
  !     Solutions located at the endpoints in time.
  ! =1: linear prolongation/restriction. Central differences. 
  !     Dual Solutions located in time according to the theta scheme.
  integer, intent(in) :: ctimeProjection
  
  ! Part of the solution which should be projected with
  ! the above projection. This is usually =1 except for some controls, e.g.,
  ! where one part has to be projected in a different way than another.
  integer, intent(in) :: isubspace
!</input>

!<output>
  ! Space-time projection structure.
  type(t_sptiProjHierarchy), intent(out) :: rprojHier
!</output>

!</subroutine>
    integer :: i,itimeOrder,cthetaschemetype
    type(t_matrixScalar) :: rprolmat2

    ! Remember the physics; necessary so we know how and what to project
    rprojHier%p_rphysics => rphysics

    ! Remember the discretisation and projection hierarchy in space.
    rprojHier%p_rspaceTimeHierarchy => rspaceTimeHierarchy
    rprojHier%p_rprojHierarchySpace => rprojHierarchySpace
    rprojHier%p_rtimeCoarseDiscr => rspaceTimeHierarchy%p_rtimeHierarchy%p_RtimeLevels(1)
    
    cthetaschemetype = rprojHier%p_rtimeCoarseDiscr%itag
    
    ! The component part is =0 usually
    nullify(rprojHier%p_Icomponents)

    ! Set ctimeProjection
    rprojHier%ctimeProjection = ctimeProjection

    if (rprojHier%ctimeProjection .eq. -1) then
      ! Automatic mode. Select the order based on the time stepping scheme.
      call tdiscr_getOrder(rprojHier%p_rtimeCoarseDiscr,itimeOrder)
      select case (itimeOrder)
      case (0)
        rprojHier%ctimeProjection = 0
      case (1)
        rprojHier%ctimeProjection = 1
      case (2)
        rprojHier%ctimeProjection = 2
      case default
        rprojHier%ctimeProjection = 1
      end select
      
      ! If our special 1-step scheme is activated, reduce the order to 1 in order
      ! to activate the corresponding prol/rest.
      if (cthetaschemetype .eq. 1) then
        rprojHier%ctimeProjection = 2
      end if
      
    end if
    
    ! Create prolongation and restriction matrices.
    allocate (rprojHier%p_RprolongationMat(max(1,rspaceTimeHierarchy%nlevels-1)))
    allocate (rprojHier%p_RrestrictionMat(max(1,rspaceTimeHierarchy%nlevels-1)))
    allocate (rprojHier%p_RinterpolationMat(max(1,rspaceTimeHierarchy%nlevels-1)))
    
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
      case (0,1,2)
      
        if (i .lt. rspaceTimeHierarchy%nlevels) then
        
          select case (cspace)
          
          ! Create the prolongation matrix in p_RprolongationMat.
          ! Create the transposed restriction matrix in rprolmat2 and transpose.
          ! Create the interpolation matrix in p_RinterpolationMat.
          
          ! ===============================================
          ! Primal space
          ! ===============================================
          case (CCSPACE_PRIMAL)
            call sptipr_getProlMatrixPrimal(rspaceTimeHierarchy,i,rprojHier%ctimeProjection,&
                rprojHier%p_RprolongationMat(i))
                
            call sptipr_getProlMatrixDual(rspaceTimeHierarchy,i,rprojHier%ctimeProjection,&
                rprolmat2)

            call sptipr_getInterpMatrixPrimal(rspaceTimeHierarchy,i,rprojHier%ctimeProjection,&
                rprojHier%p_RinterpolationMat(i))

            !call matio_writeMatrixHR (rprojHier%p_RprolongationMatPrimal(i), "pmat",&
            !    .true., 0, "matrixp."//trim(sys_siL(i,10)), "(E20.10)")
          
            !call matio_writeMatrixHR (rprojHier%p_RinterpolationMat(i), "pmat",&
            !    .true., 0, "imatrixp."//trim(sys_siL(i,10)), "(E20.10)")
            
          ! ===============================================
          ! Dual space
          ! ===============================================
          case (CCSPACE_DUAL)
            call sptipr_getProlMatrixDual(rspaceTimeHierarchy,i,rprojHier%ctimeProjection,&
                rprojHier%p_RprolongationMat(i))

            call sptipr_getProlMatrixPrimal(rspaceTimeHierarchy,i,rprojHier%ctimeProjection,&
                rprolmat2)

            call sptipr_getInterpMatrixDual(rspaceTimeHierarchy,i,rprojHier%ctimeProjection,&
                rprojHier%p_RinterpolationMat(i))

            !call matio_writeMatrixHR (rprojHier%p_RprolongationMatDual(i), "dmat",&
            !    .true., 0, "matrixd."//trim(sys_siL(i,10)), "(E20.10)")

            !call matio_writeMatrixHR (rprojHier%p_RinterpolationMat(i), "pmat",&
            !    .true., 0, "imatrixp."//trim(sys_siL(i,10)), "(E20.10)")

          ! ===============================================
          ! Control space
          ! ===============================================
          case (CCSPACE_CONTROL)

            ! If there are multiple controls, we have to allocate rprojHier%p_Icomponents
            ! and put there which components have to be projected like what.
            ! isubspace defines the number of the control currently being defined,
            ! so,. e.g., 1=distributed control, 2=boundary control,...

            select case (rphysics%cequation)
            
            ! -------------------------------------------------------------
            ! Stokes/Navier Stokes.
            ! -------------------------------------------------------------
            case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
              if ((roptcontrol%dalphaC .ge. 0.0_DP) .and. (isubspace .eq. 1)) then
                call sptipr_getProlMatrixDual(rspaceTimeHierarchy,i,rprojHier%ctimeProjection,&
                    rprojHier%p_RprolongationMat(i))

                call sptipr_getProlMatrixPrimal(rspaceTimeHierarchy,i,rprojHier%ctimeProjection,&
                    rprolmat2)

                call sptipr_getInterpMatrixDual(rspaceTimeHierarchy,i,rprojHier%ctimeProjection,&
                    rprojHier%p_RinterpolationMat(i))
              end if

            ! -------------------------------------------------------------
            ! Heat equation
            ! -------------------------------------------------------------
            case (CCEQ_HEAT2D)
              if ((roptcontrol%dalphaC .ge. 0.0_DP) .and. (isubspace .eq. 1)) then
                call sptipr_getProlMatrixDual(rspaceTimeHierarchy,i,rprojHier%ctimeProjection,&
                    rprojHier%p_RprolongationMat(i))

                call sptipr_getProlMatrixPrimal(rspaceTimeHierarchy,i,rprojHier%ctimeProjection,&
                    rprolmat2)

                call sptipr_getInterpMatrixDual(rspaceTimeHierarchy,i,rprojHier%ctimeProjection,&
                    rprojHier%p_RinterpolationMat(i))
              end if
                  
            end select
            
            !call matio_writeMatrixHR (rprojHier%p_RprolongationMatDual(i), "dmat",&
            !    .true., 0, "matrixc."//trim(sys_siL(i,10)), "(E20.10)")

            !call matio_writeMatrixHR (rprojHier%p_RinterpolationMat(i), "pmat",&
            !    .true., 0, "imatrixp."//trim(sys_siL(i,10)), "(E20.10)")
            
          end select
          
          ! Restriction matrices are obtained by transposing the prolongation
          ! matrices and exchanging the primal/dual matrices.
          call lsyssc_transposeMatrix (rprolmat2,&
              rprojHier%p_RrestrictionMat(i),LSYSSC_TR_ALL)
              
          call lsyssc_releaseMatrix (rprolmat2)

          if (rprojHier%p_RrestrictionMat(i)%neq .ne. rprojHier%p_RrestrictionMat(i)%ncols) then
            ! The restriction matrices have to be divided by 2 as they are
            ! finite difference restrictions, not finite element restrictions!
            ! Exception: In case of the identity matrix (pure restriction in space),
            ! no division must be applied.
            call lsyssc_scaleMatrix (rprojHier%p_RrestrictionMat(i),0.5_DP)
          end if
          
        end if
        
      case default
        call output_line("Undefined projection.",&
            OU_CLASS_ERROR,OU_MODE_STD,"sptipr_initProjection")
        call sys_halt()

      end select
          
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
      call lsyssc_releaseMatrix(rprojHier%p_RprolongationMat(i))
      call lsyssc_releaseMatrix(rprojHier%p_RrestrictionMat(i))
      call lsyssc_releaseMatrix(rprojHier%p_RinterpolationMat(i))
    end do

    deallocate (rprojHier%p_RprolongationMat)
    deallocate (rprojHier%p_RrestrictionMat)
    deallocate (rprojHier%p_RinterpolationMat)
    
    ! Clean up the spatial projection structure
    nullify(rprojHier%p_rprojHierarchySpace)
    nullify(rprojHier%p_rspaceTimeHierarchy)
    nullify(rprojHier%p_rtimeCoarseDiscr)
    
    if (associated(rprojHier%p_Icomponents)) deallocate(rprojHier%p_Icomponents)
    
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

    ! Get the time levels of the two levels where we have to interpolate
    ! inbetween. The number of timesteps gives us the size of the matrix.
    call sth_getLevel (rspaceTimeHierarchy,ilevel,p_rtimeDiscr=p_rtimeDiscrCoarse)
    call sth_getLevel (rspaceTimeHierarchy,ilevel+1,p_rtimeDiscr=p_rtimeDiscrFine)

    ! At first, create an empty matrix.
    ! #rows = #time dofs of the fine level.
    ! #columns = #time dofs of the coarse level.
    ndofCoarse = tdiscr_igetNDofGlob(p_rtimeDiscrCoarse)
    ndofFine = tdiscr_igetNDofGlob(p_rtimeDiscrFine)
    call lsyssc_createEmptyMatrixStub (rprolMatrix,LSYSSC_MATRIX9,ndofFine,ndofCoarse)
    
    if (ndofCoarse .eq. ndofFine) then
      ! No time prolongation, identity matrix.
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
      
      ! Fill KLD, KCOL and DA
      do icol = 1,ndofFine
        p_Da(icol) = 1.0_DP
        p_Kcol(icol) = icol
        p_Kld(icol) = icol
      end do
      
      p_Kld(ndofFine+1) = ndofFine+1
      
      return
    end if
    
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

    case (1,2)
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
      
    end select
            
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

    call lsyssc_createEmptyMatrixStub (rprolMatrix,LSYSSC_MATRIX9,ndofFine,ndofCoarse)

    if (ndofCoarse .eq. ndofFine) then
      ! No time prolongation, identity matrix.
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
      
      ! Fill KLD, KCOL and DA
      do icol = 1,ndofFine
        p_Da(icol) = 1.0_DP
        p_Kcol(icol) = icol
        p_Kld(icol) = icol
      end do
      
      p_Kld(ndofFine+1) = ndofFine+1
      
      return
    end if
            
    ! Now depending on the order, create the matrix.
    select case (ctimeProjection)

    case (2)
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
      
      if (ndofCoarse .lt. 3) then
        call output_line("For this projection, there must be at least 3 timesteps available" &
            //" on the time coarse mesh!",&
            OU_CLASS_ERROR,OU_MODE_STD,"sptipr_getProlMatrixDual")
        call sys_halt()
      end if

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
      
    case default
      call sptipr_getProlMatrixPrimal (rspaceTimeHierarchy,ilevel,ctimeProjection,rprolMatrix)

    end select
            
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
    call lsyssc_createEmptyMatrixStub (rprolMatrix,LSYSSC_MATRIX9,ndofCoarse,ndofFine)

    if (ndofCoarse .eq. ndofFine) then
      ! No time prolongation, identity matrix.
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
      
      ! Fill KLD, KCOL and DA
      do icol = 1,ndofFine
        p_Da(icol) = 1.0_DP
        p_Kcol(icol) = icol
        p_Kld(icol) = icol
      end do
      
      p_Kld(ndofFine+1) = ndofFine+1
      
      return
    end if
            
    ! Now depending on the order, create the matrix.
    select case (ctimeProjection)
    case (0,1,2)
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
    call lsyssc_createEmptyMatrixStub (rprolMatrix,LSYSSC_MATRIX9,ndofCoarse,ndofFine)

    if (ndofCoarse .eq. ndofFine) then
      ! No time prolongation, identity matrix.
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
      
      ! Fill KLD, KCOL and DA
      do icol = 1,ndofFine
        p_Da(icol) = 1.0_DP
        p_Kcol(icol) = icol
        p_Kld(icol) = icol
      end do
      
      p_Kld(ndofFine+1) = ndofFine+1
      
      return
    end if
            
    ! Now depending on the order, create the matrix.
    select case (ctimeProjection)

    case (2)
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
            
  end subroutine

  ! ***************************************************************************

    subroutine getSpatialVector (rprojHier,rcoarseVector,iindex,rx,&
        rtempVecFineScalar,ispacelevelcoarse,ispacelevelfine)
    
    ! Extracts the spatial subvector iindex from rcoarseVector and puts it
    ! into rx. If necessary, the vector is prolongated in space.
    
    ! A space/time interlevel projection structure that configures the
    ! prolongation/restriction in space/time.
    type(t_sptiProjHierarchy), intent(IN) :: rprojHier

    ! Space-time source vector
    type(t_spaceTimeVectorAccess), intent(inout) :: rcoarseVector
    
    ! Index of the subvector
    integer, intent(in) :: iindex
    
    ! Space vector output.
    type(t_vectorBlock), intent(inout) :: rx
    
    ! Temp vectors
    type(t_vectorScalar), intent(inout) :: rtempVecFineScalar
    
    ! Level of the coarse and (destination) fine vector.
    integer, intent(in) :: ispacelevelcoarse
    integer, intent(in) :: ispacelevelfine

      ! local variables
      type(t_vectorBlock), pointer :: p_rtempVecCoarse
      
      ! Get the coarse grid vector  
      call sptivec_getVectorFromPool (rcoarseVector, iindex, p_rtempVecCoarse)

      ! Interpolate to the current level. Put the result to rx.
      if (ispacelevelcoarse .ne. ispacelevelfine) then
        ! Space + time
        call mlprj_performProlongation (&
            rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
            p_rtempVecCoarse,rx,rtempVecFineScalar)
      else
        ! Only time
        call lsysbl_copyVector (p_rtempVecCoarse,rx)
      end if
      
    end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_performProlongation (rprojHierBlock,ilevelfine,rcoarseVector, &
      rfineVector)
  
!<description>
  ! Performs a prolongation for a given space/time vector (i.e. a projection
  ! in the primal space where the solution lives). The vector
  ! rcoarseVector on a coarser space/time mesh is projected to the vector
  ! rfineVector on a finer space/time mesh.
  ! rprojHierBlock%p_RprojHier configures how the transfer is performed.
  ! This projection structure rprojHierBlock%p_RprojHier must correspond to the space/time
  ! discretisation of rcoarseVector and rfineVector.
!</description>

!<input>
  ! Array of space/time interlevel projection structures that configure the
  ! prolongation/restriction in space/time.
  type(t_sptiProjHierarchyBlock), intent(in) :: rprojHierBlock
  
  ! Index of the fine mesh.
  integer, intent(in) :: ilevelfine

  ! Coarse grid vector
  type(t_spacetimeVectorAccess), intent(inout) :: rcoarseVector
!</input>

!<output>
  ! Fine grid vector
  type(t_spacetimeVectorAccess), intent(inout) :: rfineVector
!</output>

!</subroutine>

    ! Local variables
    type(t_spaceTimeVectorAccess) :: raccessPool
    integer :: ispacelevelfine,ispacelevelcoarse,itimelevelfine,itimelevelcoarse
    type(t_vectorScalar) :: rtempVecFineScalar
    type(t_vectorBlock), pointer :: p_rx,p_rdestVector,p_rtempVecFine
    integer :: irow, icol, neq, icomp, imat, icompidx, i, neqfine
    type(t_matrixScalar), pointer :: p_rtimeMatrix
    real(DP), dimension(:), pointer :: p_Da,p_DdataCoarse,p_DdataFine
    integer, dimension(:), pointer :: p_Kcol, p_Kld
    real(DP) :: dscale

    ! DEBUG!!!
    !do istep=1,rcoarseVector%NEQtime
    !  call sptivec_setSubvectorConstant(rcoarseVector,istep,real(istep,DP))
    !end do

    ! Get the space-time coarse and fine discretisations.
    call sth_getLevel (rprojHierBlock%p_RprojHier(1)%p_rspaceTimeHierarchy,ilevelfine-1,&
      ispaceLevel=ispacelevelcoarse,itimeLevel=itimelevelcoarse)
    call sth_getLevel (rprojHierBlock%p_RprojHier(1)%p_rspaceTimeHierarchy,ilevelfine,&
      ispaceLevel=ispacelevelfine,itimeLevel=itimelevelfine)
    
    ! Allocate temp memory
    i = 0
    do icompidx = 1,rprojHierBlock%ncount
      select case (rprojHierBlock%p_RprojHier(icompidx)%ctimeProjection)
      case (0,1,2)
        i = max(i,3)
      case default
        i = max(i,10)
      end select
    end do
    
    ! The access pool receives the prolongated solution and is
    ! not bounded to and space-time vector!
    call sptivec_createAccessPool (rfineVector%p_rspaceDiscr,raccessPool,i+1)
    
    ! We need one temp vector. Take it from the pool and and associate
    ! the index neqfine -- this will for sure not appear in the loop.
    ! Lock the vector so it is not overwritten.
    neqfine = rfineVector%p_rspaceTimeVector%NEQtime
    call sptivec_getFreeBufferFromPool (raccessPool,neqfine+1,p_rtempVecFine)
    call sptivec_lockVecInPool (raccessPool,neqfine+1)
    
    ! We need a scalar representation of the temp vector
    call lsysbl_createScalarFromVec (p_rtempVecFine,rtempVecFineScalar)
    
    ! Prolongation means, we multiply with the prolongation matrix.
    ! y and xi have to be multiplied with the primal prol. matrix,
    ! lambda and p with the dual. Note that for implicit Euler,
    ! the primal and dual prolongation matrix is the same, while
    ! for CN, y and xi are at the same points in time as well as
    ! lambda and p: p is between y and thus at the same points
    ! in time as lambda, not as y!
    
    neq = rprojHierBlock%p_RprojHier(1)%p_RprolongationMat(ilevelfine-1)%neq
    
    ! Apply the multiplication.
    ! The rows in the matrix correspond to the time fine mesh, the columns
    ! to the time coarse mesh.
    do irow = 1,neq
    
      ! Clear the destination
      call sptivec_getFreeBufferFromPool (rfineVector,irow,p_rdestVector)
      call lsysbl_clearVector (p_rdestVector)
      
      ! DEBUG!!!
      call lsysbl_getbase_double (p_rdestVector,p_DdataFine)
      
      ! Loop over the matrices configuring the prolongation for
      ! all the components
      do imat = 1,rprojHierBlock%ncount
      
        ! Get the matrix
        p_rtimeMatrix => rprojHierBlock%p_RprojHier(imat)%p_RprolongationMat(ilevelfine-1)
        call lsyssc_getbase_double (p_rtimeMatrix,p_Da)
        call lsyssc_getbase_Kcol (p_rtimeMatrix,p_Kcol)
        call lsyssc_getbase_Kld (p_rtimeMatrix,p_Kld)
        dscale = p_rtimeMatrix%dscaleFactor

        do icol = p_Kld(irow),p_Kld(irow+1)-1
        
          ! Try to get the vector from the vector pool. Saves time.
          call sptivec_getVectorFromPool(raccessPool,p_Kcol(icol),p_rx)
          if (.not. associated(p_rx)) then
            ! No, we have to fetch/calculate the vector.
            !
            ! Get a buffer where to save it.
            call sptivec_getFreeBufferFromPool (raccessPool,p_Kcol(icol),p_rx)
            
            ! Read the source vector. The column of the matrix specifies
            ! the timestep.
            call getSpatialVector (rprojHierBlock%p_RprojHier(imat),&
                rcoarseVector,p_Kcol(icol),p_rx,&
                rtempVecFineScalar,ispacelevelcoarse,ispacelevelfine)
          end if
          
          ! DEBUG!!!
          call lsysbl_getbase_double (p_rx,p_DdataCoarse)
          
          ! Now, rx is the time vector at timestep icol. Weighted multiplication
          ! into rtempVecFine for y and xi.
          ! Multiply either all components or the specified subcomponents.
          if (.not. associated(rprojHierBlock%p_RprojHier(imat)%p_Icomponents)) then

            call lsysbl_vectorLinearComb(p_rx,p_rdestVector,dscale*p_Da(icol),1.0_DP)

          else

            do icompidx = 1,ubound(rprojHierBlock%p_RprojHier(imat)%p_Icomponents,1)
              icomp = rprojHierBlock%p_RprojHier(imat)%p_Icomponents(icompidx)
              call lsyssc_vectorLinearComb(&
                  p_rx%Rvectorblock(icomp),p_rdestVector%RvectorBlock(icomp),&
                  dscale*p_Da(icol),1.0_DP)
            end do
            
          end if
          
        end do
        
      end do
      
      ! Vector finished.
      call sptivec_commitVecInPool (rfineVector,irow)

    end do
    
    ! Release the buffer.
    call sptivec_releaseAccessPool(raccessPool)
    
    call lsyssc_releaseVector (rtempVecFineScalar)
    
  end subroutine
     
  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_performRestriction (rprojHierBlock,ilevelfine,rcoarseVector, &
      rfineVector)
  
!<description>
  ! Performs a restriction for a given space/time vector (i.e. a projection
  ! in the dual space where the RHS vector lives). The vector
  ! rfineVector on a finer grid is projected to the vector
  ! rcoarseVector on a coarser space/time mesh.
  ! rprojHierBlock%p_RprojHier configures how the transfer is performed.
  ! This projection structure rprojHierBlock%p_RprojHier must correspond to the space/time
  ! discretisation of rcoarseVector and rfineVector.
!</description>

!<input>
  ! Array of space/time interlevel projection structures that configure the
  ! prolongation/restriction in space/time.
  type(t_sptiProjHierarchyBlock), intent(in) :: rprojHierBlock

  ! Index of the fine mesh.
  integer, intent(in) :: ilevelfine

  ! Fine grid vector
  type(t_spacetimeVectorAccess), intent(inout) :: rfineVector
!</input>

!<output>
  ! Coarse grid vector
  type(t_spacetimeVectorAccess), intent(INOUT) :: rcoarseVector
!</output>

!</subroutine>

    ! Local variables
    type(t_spaceTimeVectorAccess) :: raccessPool
    integer :: ispacelevelfine,ispacelevelcoarse,itimelevelfine,itimelevelcoarse
    type(t_vectorScalar) :: rtempVecFineScalar
    type(t_vectorBlock), pointer :: p_rx,p_rdestVector, p_rtempVecFine
    integer :: irow, icol, neq, icomp, imat, icompidx, i, neqfine
    type(t_matrixScalar), pointer :: p_rtimeMatrix
    real(DP), dimension(:), pointer :: p_Da,p_DdataCoarse,p_DdataFine
    integer, dimension(:), pointer :: p_Kcol, p_Kld
    real(DP) :: dscale

    ! DEBUG!!!
    !call sptivec_setConstant(rfineVector,1.0_DP)

    ! Get the space-time coarse and fine discretisations.
    call sth_getLevel (rprojHierBlock%p_RprojHier(1)%p_rspaceTimeHierarchy,ilevelfine-1,&
      ispaceLevel=ispacelevelcoarse,itimeLevel=itimelevelcoarse)
    call sth_getLevel (rprojHierBlock%p_RprojHier(1)%p_rspaceTimeHierarchy,ilevelfine,&
      ispaceLevel=ispacelevelfine,itimeLevel=itimelevelfine)
    
    ! Allocate temp memory
    i = 0
    do icompidx = 1,rprojHierBlock%ncount
      select case (rprojHierBlock%p_RprojHier(icompidx)%ctimeProjection)
      case (0,1,2)
        i = max(i,3)
      case default
        i = max(i,10)
      end select
    end do
    
    call sptivec_createAccessPool (rfineVector%p_rspaceDiscr,raccessPool,i+1)
    
    ! We need one temp vector. Take it from the pool and and associate
    ! the index neqfine -- this will for sure not appear in the loop.
    ! Lock the vector so it is not overwritten.
    neqfine = rfineVector%p_rspaceTimeVector%NEQtime
    call sptivec_getFreeBufferFromPool (raccessPool,neqfine+1,p_rtempVecFine)
    call sptivec_lockVecInPool (raccessPool,neqfine+1)
    
    ! We need a scalar representation of the temp vector
    call lsysbl_createScalarFromVec (p_rtempVecFine,rtempVecFineScalar)
    
    ! Prolongation means, we multiply with the prolongation matrix.
    ! y and xi have to be multiplied with the primal prol. matrix,
    ! lambda and p with the dual. Note that for implicit Euler,
    ! the primal and dual prolongation matrix is the same, while
    ! for CN, y and xi are at the same points in time as well as
    ! lambda and p: p is between y and thus at the same points
    ! in time as lambda, not as y!

    neq = rprojHierBlock%p_RprojHier(1)%p_RrestrictionMat(ilevelfine-1)%neq
    
    ! Apply the multiplication.
    ! The rows in the matrix correspond to the time fine mesh, the columns
    ! to the time coarse mesh.
    do irow = 1,neq

      ! Clear the destination
      call sptivec_getFreeBufferFromPool (rfineVector,irow,p_rdestVector)
      call lsysbl_clearVector (p_rdestVector)
      
      ! DEBUG!!!
      call lsysbl_getbase_double (p_rdestVector,p_DdataFine)
      
      ! Loop over the matrices configuring the prolongation for
      ! all the components
      do imat = 1,rprojHierBlock%ncount
      
        ! Get the matrix
        p_rtimeMatrix => rprojHierBlock%p_RprojHier(imat)%p_RrestrictionMat(ilevelfine-1)
        call lsyssc_getbase_double (p_rtimeMatrix,p_Da)
        call lsyssc_getbase_Kcol (p_rtimeMatrix,p_Kcol)
        call lsyssc_getbase_Kld (p_rtimeMatrix,p_Kld)
        dscale = p_rtimeMatrix%dscaleFactor

        ! Primal space: y,xi
        do icol = p_Kld(irow),p_Kld(irow+1)-1
        
          ! Get the fine grid vector using the vector pool as buffer. Saves time.
          call sptivec_getVectorFromPool(rfineVector,p_Kcol(icol),p_rx)
          
          call lsysbl_getbase_double (p_rx,p_DdataCoarse)
          
          ! Now, rx is the time vector at timestep icol. Weighted multiplication
          ! into rtempVecFine for y and xi.

          ! Multiply either all components or the specified subcomponents.
          if (.not. associated(rprojHierBlock%p_RprojHier(imat)%p_Icomponents)) then

            call lsysbl_vectorLinearComb(p_rx,p_rdestVector,dscale*p_Da(icol),1.0_DP)

          else

            do icompidx = 1,ubound(rprojHierBlock%p_RprojHier(imat)%p_Icomponents,1)
              icomp = rprojHierBlock%p_RprojHier(imat)%p_Icomponents(icompidx)
              call lsyssc_vectorLinearComb(&
                  p_rx%Rvectorblock(icomp),p_rdestVector%RvectorBlock(icomp),&
                  dscale*p_Da(icol),1.0_DP)
            end do
            
          end if
              
        end do
        
      end do
        
      ! Vector finished.
      ! Pass rprojHierBlock%p_RprojHier(1); it defines the spatial projection.
      call setSpatialVector (rprojHierBlock%p_RprojHier(1),p_rdestVector,rcoarseVector,irow,&
          rtempVecFineScalar,ispacelevelcoarse,ispacelevelfine)

    end do

    ! Release the buffer.
    call sptivec_releaseAccessPool(raccessPool)
    
    call lsyssc_releaseVector (rtempVecFineScalar)

    ! DEBUG!!!
    !call sptivec_saveToFileSequence (rcoarseVector,"(""./ns/coarse.txt."",I5.5)",.true.,&
    !    rtempVecCoarse)
    !call sys_halt()
        
  contains

    subroutine setSpatialVector (rprojHier,rx,rcoarseVector,iindex,&
        rtempVecFineScalar,ispacelevelcoarse,ispacelevelfine)
    
    ! Saves the spatial subvector iindex from rx to rfineVector.
    ! If necessary, the vector is restricted in space.
    
    ! A space/time interlevel projection structure that configures the
    ! prolongation/restriction in space/time.
    type(t_sptiProjHierarchy), intent(in) :: rprojHier

    ! Space-time destination vector
    type(t_spaceTimeVectorAccess), intent(inout) :: rcoarseVector
    
    ! Index of the subvector
    integer, intent(in) :: iindex
    
    ! Space vector output.
    type(t_vectorBlock), intent(inout) :: rx
    
    ! Temp vectors
    type(t_vectorScalar), intent(inout) :: rtempVecFineScalar
    
    ! Level of the coarse and (destination) fine vector.
    integer, intent(in) :: ispacelevelcoarse
    integer, intent(in) :: ispacelevelfine

      ! local variables
      type(t_vectorBlock), pointer :: p_rtempVecCoarse
  
      ! Get the coarse grid vector
      call sptivec_getVectorFromPool (rcoarseVector, iindex, p_rtempVecCoarse)
  
      if (ispacelevelcoarse .ne. ispacelevelfine) then
        ! Space + time
        call mlprj_performRestriction (&
            rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
            p_rtempVecCoarse,rx,rtempVecFineScalar)
      else
        ! Only time
        call lsysbl_copyVector (rx,p_rtempVecCoarse)
      end if
      
      ! Save the vector
      call sptivec_commitVecInPool (rcoarseVector, iindex)

    end subroutine
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_performInterpolation (rprojHierBlock,ilevelfine,rcoarseVector, &
      rfineVector)
  
!<description>
  ! Performs an interpolation for a given space/time vector to a lower level.
  ! The solution vector rfineVector on a finer grid is projected to the vector
  ! rcoarseVector on a coarser space/time mesh.
  ! rprojHierBlock%p_RprojHier configures how the transfer is performed.
  ! This projection structure rprojHierBlock%p_RprojHier must correspond to the space/time
  ! discretisation of rcoarseVector and rfineVector.
!</description>

!<input>
  ! Array of space/time interlevel projection structures that configure the
  ! prolongation/restriction in space/time.
  type(t_sptiProjHierarchyBlock), intent(in) :: rprojHierBlock

  ! Index of the fine mesh.
  integer, intent(in) :: ilevelfine

  ! Fine grid vector
  type(t_spacetimeVectorAccess), intent(inout) :: rfineVector
!</input>

!<output>
  ! Coarse grid vector
  type(t_spacetimeVectorAccess), intent(inout) :: rcoarseVector
!</output>

!</subroutine>

    ! Local variables
    type(t_spaceTimeVectorAccess) :: raccessPool
    integer :: ispacelevelfine,ispacelevelcoarse,itimelevelfine,itimelevelcoarse
    type(t_vectorScalar) :: rtempVecFineScalar
    type(t_vectorBlock), pointer :: p_rx,p_rdestVector,p_rtempVecFine
    integer :: irow, icol, neq, icomp, imat, icompidx, i, neqfine
    type(t_matrixScalar), pointer :: p_rtimeMatrix
    real(DP), dimension(:), pointer :: p_Da,p_DdataCoarse,p_DdataFine
    integer, dimension(:), pointer :: p_Kcol, p_Kld
    real(DP) :: dscale

    ! DEBUG!!!
    !call sptivec_setConstant(rfineVector,1.0_DP)

    ! Get the space-time coarse and fine discretisations.
    call sth_getLevel (rprojHierBlock%p_RprojHier(1)%p_rspaceTimeHierarchy,ilevelfine-1,&
      ispaceLevel=ispacelevelcoarse,itimeLevel=itimelevelcoarse)
    call sth_getLevel (rprojHierBlock%p_RprojHier(1)%p_rspaceTimeHierarchy,ilevelfine,&
      ispaceLevel=ispacelevelfine,itimeLevel=itimelevelfine)
    
    ! Allocate temp memory
    i = 0
    do icompidx = 1,rprojHierBlock%ncount
      select case (rprojHierBlock%p_RprojHier(icompidx)%ctimeProjection)
      case (0,1,2)
        i = max(i,3)
      case default
        i = max(i,10)
      end select
    end do
    
    call sptivec_createAccessPool (rfineVector%p_rspaceDiscr,raccessPool,i+1)
    
    ! We need one temp vector. Take it from the pool and and associate
    ! the index neqfine -- this will for sure not appear in the loop.
    ! Lock the vector so it is not overwritten.
    neqfine = rfineVector%p_rspaceTimeVector%NEQtime
    call sptivec_getFreeBufferFromPool (raccessPool,neqfine+1,p_rtempVecFine)
    call sptivec_lockVecInPool (raccessPool,neqfine+1)
    
    ! We need a scalar representation of the temp vector
    call lsysbl_createScalarFromVec (p_rtempVecFine,rtempVecFineScalar)

    ! Prolongation means, we multiply with the prolongation matrix.
    ! y and xi have to be multiplied with the primal prol. matrix,
    ! lambda and p with the dual. Note that for implicit Euler,
    ! the primal and dual prolongation matrix is the same, while
    ! for CN, y and xi are at the same points in time as well as
    ! lambda and p: p is between y and thus at the same points
    ! in time as lambda, not as y!

    neq = rprojHierBlock%p_RprojHier(1)%p_RinterpolationMat(ilevelfine-1)%neq
    
    ! Apply the multiplication.
    ! The rows in the matrix correspond to the time fine mesh, the columns
    ! to the time coarse mesh.
    do irow = 1,neq

      ! Clear the destination
      call sptivec_getFreeBufferFromPool (rfineVector,irow,p_rdestVector)
      call lsysbl_clearVector (p_rdestVector)
      
      ! DEBUG!!!
      call lsysbl_getbase_double (p_rdestVector,p_DdataFine)
      
      ! Loop over the matrices configuring the prolongation for
      ! all the components
      do imat = 1,rprojHierBlock%ncount
      
        ! Get the matrix
        p_rtimeMatrix => rprojHierBlock%p_RprojHier(imat)%p_RinterpolationMat(ilevelfine-1)
        call lsyssc_getbase_double (p_rtimeMatrix,p_Da)
        call lsyssc_getbase_Kcol (p_rtimeMatrix,p_Kcol)
        call lsyssc_getbase_Kld (p_rtimeMatrix,p_Kld)
        dscale = p_rtimeMatrix%dscaleFactor

        ! Primal space: y,xi
        do icol = p_Kld(irow),p_Kld(irow+1)-1
        
          ! Get the fine grid vector using the vector pool as buffer. Saves time.
          call sptivec_getVectorFromPool(rfineVector,p_Kcol(icol),p_rx)
          
          call lsysbl_getbase_double (p_rx,p_DdataCoarse)
          
          ! Now, rx is the time vector at timestep icol. Weighted multiplication
          ! into rtempVecFine for y and xi.

          ! Multiply either all components or the specified subcomponents.
          if (.not. associated(rprojHierBlock%p_RprojHier(imat)%p_Icomponents)) then

            call lsysbl_vectorLinearComb(p_rx,p_rdestVector,dscale*p_Da(icol),1.0_DP)

          else

            do icompidx = 1,ubound(rprojHierBlock%p_RprojHier(imat)%p_Icomponents,1)
              icomp = rprojHierBlock%p_RprojHier(imat)%p_Icomponents(icompidx)
              call lsyssc_vectorLinearComb(&
                  p_rx%Rvectorblock(icomp),p_rdestVector%RvectorBlock(icomp),&
                  dscale*p_Da(icol),1.0_DP)
            end do
            
          end if
              
        end do
        
      end do
      
      ! Vector finished.
      ! Pass rprojHierBlock%p_RprojHier(1); it defines the spatial projection.
      call setSpatialVector (rprojHierBlock%p_RprojHier(1),p_rdestVector,rcoarseVector,irow,&
          rtempVecFineScalar,ispacelevelcoarse,ispacelevelfine)

    end do
    
    ! Release the buffer.
    call sptivec_releaseAccessPool(raccessPool)
    
    call lsyssc_releaseVector (rtempVecFineScalar)

    ! DEBUG!!!
    !call sptivec_saveToFileSequence (rcoarseVector,"(""coarse.txt."",I5.5)",.true.,&
    !    rtempVecFine)

  contains

    subroutine setSpatialVector (rprojHier,rx,rcoarseVector,iindex,&
        rtempVecFineScalar,ispacelevelcoarse,ispacelevelfine)
    
    ! Saves the spatial subvector iindex from rx to rfineVector.
    ! If necessary, the vector is restricted in space.
    
    ! A space/time interlevel projection structure that configures the
    ! prolongation/restriction in space/time.
    type(t_sptiProjHierarchy), intent(in) :: rprojHier

    ! Space-time destination vector
    type(t_spaceTimeVectorAccess), intent(inout) :: rcoarseVector
    
    ! Index of the subvector
    integer, intent(in) :: iindex
    
    ! Space vector output.
    type(t_vectorBlock), intent(inout) :: rx
    
    ! Temp vectors
    type(t_vectorScalar), intent(inout) :: rtempVecFineScalar
    
    ! Level of the coarse and (destination) fine vector.
    integer, intent(in) :: ispacelevelcoarse
    integer, intent(in) :: ispacelevelfine
  
      ! local variables
      type(t_vectorBlock), pointer :: p_rtempVecCoarse

      ! Get the coarse grid vector
      call sptivec_getVectorFromPool (rcoarseVector, iindex, p_rtempVecCoarse)

      if (ispacelevelcoarse .ne. ispacelevelfine) then
        ! Space + time
        call mlprj_performInterpolation (&
            rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
            p_rtempVecCoarse,rx,rtempVecFineScalar)
      else
        ! Only time
        call lsysbl_copyVector (rx,p_rtempVecCoarse)
      end if

      ! Save the vector
      call sptivec_commitVecInPool (rcoarseVector, iindex)

    end subroutine

  end subroutine

end module
