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
    ! =1: first order implicit Euler, 
    ! =2: Full 1-step Theta scheme up to second order (Crank Nicolson)
    integer :: itimeOrder = 0
    
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
    integer :: i
    real(DP), dimension(:), pointer :: p_Da
    integer, dimension(:), pointer :: p_Kcol, p_Kld

    ! Remember the discretisation and projection hierarchy in space.  
    rprojHier%p_rspaceTimeHierarchy => rspaceTimeHierarchy
    rprojHier%p_rprojHierarchySpace => rprojHierarchySpace
    rprojHier%p_rtimeCoarseDiscr => rspaceTimeHierarchy%p_rtimeHierarchy%p_RtimeLevels(1)
    
    ! Set itimeOrder <> 0 -> structure initialised.
    rprojHier%itimeOrder = -1
    if (present(iorderTimeProlRest)) rprojHier%itimeOrder = iorderTimeProlRest

    if (rprojHier%itimeOrder .eq. -1) then
      ! Automatic mode. Select the order based on the time stepping scheme.
      call tdiscr_getOrder(rprojHier%p_rtimeCoarseDiscr,rprojHier%itimeOrder)
      
      ! If our special 1-step scheme is activated, reduce iorder to 1 in order
      ! to activate the corresponding prol/rest.
      if (rprojHier%p_rtimeCoarseDiscr%itag .eq. 1) then
        rprojHier%itimeOrder = 1
      end if
      
    end if
    
    ! Create prolongation and restriction matrices.
    allocate (rprojHier%p_RprolongationMatPrimal(max(1,rspaceTimeHierarchy%nlevels-1)))
    allocate (rprojHier%p_RprolongationMatDual(max(1,rspaceTimeHierarchy%nlevels-1)))
    allocate (rprojHier%p_RrestrictionMatPrimal(max(1,rspaceTimeHierarchy%nlevels-1)))
    allocate (rprojHier%p_RrestrictionMatDual(max(1,rspaceTimeHierarchy%nlevels-1)))
    allocate (rprojHier%p_RinterpolationMatPrimal(max(1,rspaceTimeHierarchy%nlevels-1)))
    allocate (rprojHier%p_RinterpolationMatDual(max(1,rspaceTimeHierarchy%nlevels-1)))
    
    do i=1,rspaceTimeHierarchy%nlevels-1
      call sptipr_getProlMatrixPrimal(rspaceTimeHierarchy,i,rprojHier%itimeOrder,&
          rprojHier%p_RprolongationMatPrimal(i))
          
      !call matio_writeMatrixHR (rprojHier%p_RprolongationMatPrimal(i), "pmat",&
      !    .true., 0, "matrixp."//trim(sys_siL(i,10)), "(E20.10)")
      
      call sptipr_getProlMatrixDual(rspaceTimeHierarchy,i,rprojHier%itimeOrder,&
          rprojHier%p_RprolongationMatDual(i))

      !call matio_writeMatrixHR (rprojHier%p_RprolongationMatDual(i), "dmat",&
      !    .true., 0, "matrixd."//trim(sys_siL(i,10)), "(E20.10)")
          
      ! The restriction matrices are given as their transpose...
      !
      ! WARNING!!!
      ! The primal restriction matrix is the transpose of the dual prolongation matrix!
      ! The dual restriction matrix is the transpose of the primal prolongation matrix!
      ! This is because the primal RHS is located at the timesteps of the dual
      ! solution (between the primal timesteps) and vice versa!!!
      call lsyssc_transposeMatrix (rprojHier%p_RprolongationMatPrimal(i),&
          rprojHier%p_RrestrictionMatDual(i),LSYSSC_TR_ALL)
          
      call lsyssc_transposeMatrix (rprojHier%p_RprolongationMatDual(i),&
          rprojHier%p_RrestrictionMatPrimal(i),LSYSSC_TR_ALL)
          
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
      call sptipr_getInterpMatrixPrimal(rspaceTimeHierarchy,i,rprojHier%itimeOrder,&
          rprojHier%p_RinterpolationMatPrimal(i))

      !call matio_writeMatrixHR (rprojHier%p_RinterpolationMatPrimal(i), "pmat",&
      !    .true., 0, "imatrixp."//trim(sys_siL(i,10)), "(E20.10)")

      call sptipr_getInterpMatrixDual(rspaceTimeHierarchy,i,rprojHier%itimeOrder,&
          rprojHier%p_RinterpolationMatDual(i))

      !call matio_writeMatrixHR (rprojHier%p_RinterpolationMatDual(i), "pmat",&
      !    .true., 0, "imatrixd."//trim(sys_siL(i,10)), "(E20.10)")
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
    rprojHier%itimeOrder = 0

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_getProlMatrixPrimal (rspaceTimeHierarchy,ilevel,iorder,rprolMatrix)
  
!<description>
  ! Creates a prolongation matrix for the primal space between time-level 
  ! ilevel and ilevel+1 of the given space-time hierarchy.
!</description>

!<input>
  ! Underlying space-time hierarchy.
  type(t_spacetimeHierarchy), intent(in) :: rspaceTimeHierarchy
  
  ! Id of the coarse level
  integer, intent(in) :: ilevel
  
  ! Order of the prolongation.
  integer, intent(in) :: iorder
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
    integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Kdiagonal
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
        
    ! Now depending on the order, create the matrix.
    select case (iorder)
    
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

    case (3)
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
            
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_getProlMatrixDual (rspaceTimeHierarchy,ilevel,iorder,rprolMatrix)
  
!<description>
  ! Creates a prolongation matrix for the dual space between time-level 
  ! ilevel and ilevel+1 of the given space-time hierarchy.
!</description>

!<input>
  ! Underlying space-time hierarchy.
  type(t_spacetimeHierarchy), intent(in) :: rspaceTimeHierarchy
  
  ! Id of the coarse level
  integer, intent(in) :: ilevel
  
  ! Order of the prolongation.
  integer, intent(in) :: iorder
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
        
    ! Now depending on the order, create the matrix.
    select case (iorder)

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
      p_Da(2) = dtheta-0.5_DP          ! 0.0_DP
      p_Da(3) = 2.0_DP-1.5_DP*dtheta   ! 1.25_DP
      p_Da(4) = -0.5_DP+0.5_DP*dtheta  ! -0.25_DP
      
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

        p_Da(5+4*(icol-2)+0) = 1.0_DP - 0.5_DP*dtheta ! 0.75_DP
        p_Da(5+4*(icol-2)+1) = 0.5_DP*dtheta          ! 0.25_DP
        p_Da(5+4*(icol-2)+2) = 0.5_DP*dtheta          ! 0.25_DP
        p_Da(5+4*(icol-2)+3) = 1.0_DP - 0.5_DP*dtheta ! 0.75_DP
      end do

      ! The last row.
      p_Kcol(5+4*(ndofCoarse-2)+0) = ndofCoarse-1
      p_Kcol(5+4*(ndofCoarse-2)+1) = ndofCoarse

      p_Da(5+4*(ndofCoarse-2)) = -0.5_DP+0.5_DP*dtheta  ! -0.25_DP
      p_Da(5+4*(ndofCoarse-2)+1) = 2.0_DP-1.5_DP*dtheta ! 1.25_DP
      
    case default
      call sptipr_getProlMatrixPrimal (rspaceTimeHierarchy,ilevel,iorder,rprolMatrix)

    end select
            
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_getInterpMatrixPrimal (rspaceTimeHierarchy,ilevel,iorder,rprolMatrix)
  
!<description>
  ! Creates a interpolation matrix for the primal space between time-level 
  ! ilevel and ilevel+1 of the given space-time hierarchy.
!</description>

!<input>
  ! Underlying space-time hierarchy.
  type(t_spacetimeHierarchy), intent(in) :: rspaceTimeHierarchy
  
  ! Id of the coarse level
  integer, intent(in) :: ilevel
  
  ! Order of the prolongation.
  integer, intent(in) :: iorder
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
    integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Kdiagonal
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
        
    ! Now depending on the order, create the matrix.
    select case (iorder)
    case (0,1,2,3)
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

  subroutine sptipr_getInterpMatrixDual (rspaceTimeHierarchy,ilevel,iorder,rprolMatrix)
  
!<description>
  ! Creates a prolongation matrix for the dual space between time-level 
  ! ilevel and ilevel+1 of the given space-time hierarchy.
!</description>

!<input>
  ! Underlying space-time hierarchy.
  type(t_spacetimeHierarchy), intent(in) :: rspaceTimeHierarchy
  
  ! Id of the coarse level
  integer, intent(in) :: ilevel
  
  ! Order of the prolongation.
  integer, intent(in) :: iorder
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
        
    ! Now depending on the order, create the matrix.
    select case (iorder)

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
      call sptipr_getInterpMatrixPrimal (rspaceTimeHierarchy,ilevel,iorder,rprolMatrix)

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
    integer :: istep,ispacelevelfine,ispacelevelcoarse,itimelevelfine,itimelevelcoarse
    type(t_vectorScalar) :: rtempVecFineScalar,rtempVecFineScalar2
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
    select case (rprojHier%itimeOrder)
    case (0)
      call sptivec_createAccessPool (rfineVector%p_rspaceDiscr,raccessPool,3)
    case (1)
      call sptivec_createAccessPool (rfineVector%p_rspaceDiscr,raccessPool,3)
    case (2)
      call sptivec_createAccessPool (rfineVector%p_rspaceDiscr,raccessPool,3)
    case (3)
      call sptivec_createAccessPool (rfineVector%p_rspaceDiscr,raccessPool,4)
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
            p_rx%RvectorBlock(3),rtempVecFine%RvectorBlock(3),dscaleDual*p_DaDual(icol),1.0_DP)
            
        call lsyssc_vectorLinearComb(&
            p_rx%RvectorBlock(4),rtempVecFine%RvectorBlock(4),dscaleDual*p_DaDual(icol),1.0_DP)

        call lsyssc_vectorLinearComb(&
            p_rx%RvectorBlock(5),rtempVecFine%RvectorBlock(5),dscaleDual*p_DaDual(icol),1.0_DP)
        
      end do

      ! Vector finished.
      call sptivec_setTimestepData (rfineVector, irow, rtempVecFine)

    end do
    
    ! Release the buffer.
    call sptivec_releaseAccessPool(raccessPool)
    
    call lsyssc_releaseVector (rtempVecFineScalar)

    ! DEBUG!!!
    !call sptivec_saveToFileSequence (rfineVector,"(""fine.txt."",I5.5)",.true.,&
    !    rtempVecFine)

!    ! local variables
!    integer :: istep,ispacelevelfine,ispacelevelcoarse,iorder
!    type(t_vectorBlock) :: rx1,rx2,rx3
!    type(t_vectorScalar) :: rtempVecFineScalar,rtempVecFineScalar2
!
!    ! DEBUG!!!
!    real(DP), dimension(:), pointer :: p_Dx1,p_Dx3
!    real(DP), dimension(:), pointer :: p_DtempVecCoarse,p_DtempVecFine
!    real(DP), dimension(:), pointer :: p_DtempVecFineSca
!
!    ! Get the space-time coarse and fine discretisations.
!    call sth_getLevel (rprojHier%p_rspaceTimeHierarchy,ilevelfine-1,&
!      ispaceLevel=ispacelevelcoarse)
!    call sth_getLevel (rprojHier%p_rspaceTimeHierarchy,ilevelfine,&
!      ispaceLevel=ispacelevelfine)
!
!    ! Select the prolongation/restriction to use
!    iorder = rprojHier%itimeOrder
!    if (iorder .eq. -1) then
!      ! automatic mode. Select the order based on the time stepping scheme.
!      call tdiscr_getOrder(rcoarseVector%p_rtimeDiscr,iorder)
!    end if
!    
!    select case (iorder)
!    case (0)
!      ! Constant prolongation
!      !
!      ! We need two more temp vectors:
!      !
!      ! One for the current timestep
!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.FALSE.)
!      
!      ! And one for the next timestep
!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.FALSE.)
!      
!      ! We need a scalar representation of the temp vector
!      call lsysbl_createScalarFromVec (rtempVecFine,rtempVecFineScalar)
!      
!      ! DEBUG!!!
!      call lsysbl_getbase_double (rx1,p_Dx1)
!      call lsysbl_getbase_double (rx3,p_Dx3)
!      call lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
!      call lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
!      call lsyssc_getbase_double (rtempVecFineScalar,p_DtempVecFineSca)
!
!      ! Load timestep 0 into the temp vector and interpolate to the current level.
!      ! Put the result to rx3.
!      if (ispacelevelcoarse .ne. ispacelevelfine) then
!        ! Space + time
!        call sptivec_getTimestepData (rcoarseVector, 0, rtempVecCoarse)
!        call mlprj_performProlongation (&
!            rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!            rtempVecCoarse,rx3,rtempVecFineScalar)
!      else
!        ! Only time
!        call sptivec_getTimestepData (rcoarseVector, 0, rx3)
!      end if
!      
!      ! Save that vector as initial vector on the fine grid.
!      call sptivec_setTimestepData (rfineVector, 0, rx3)
!      
!      if (rcoarseVector%NEQtime .NE. rfineVector%NEQtime) then
!        
!        ! Prolongation in time.
!        !
!        ! Loop through the time steps.
!        do istep = 1,rcoarseVector%NEQtime
!      
!          ! rx3 was the vector from the last timestep. Shift it to rx1 and load
!          ! the vector for the new timestep in rx3
!          call lsysbl_copyVector (rx3,rx1)
!
!          if (ispacelevelcoarse .ne. ispacelevelfine) then
!            ! Space + time
!            call sptivec_getTimestepData (rcoarseVector, istep, rtempVecCoarse)
!            call mlprj_performProlongation (&
!                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!                rtempVecCoarse,rx3,rtempVecFineScalar)
!          else
!            ! Only time
!            call sptivec_getTimestepData (rcoarseVector, istep, rx3)
!          end if
!          
!          ! Save that vector as new vector on the fine grid.
!          call sptivec_setTimestepData (rfineVector, 2*istep, rx3)
!          
!          ! In the temp vector, create the prolongation. Use the constant prolongation.
!          ! Primal vector of rx1 -> fine vector <- Dual vector of rx3.
!          ! That's the prolongated vector.
!          call lsysbl_copyVector (rx1,rtempVecFine)
!          call lsyssc_copyVector (rx3%RvectorBlock(4),rtempVecFine%RvectorBlock(4))
!          call lsyssc_copyVector (rx3%RvectorBlock(5),rtempVecFine%RvectorBlock(5))
!          call lsyssc_copyVector (rx3%RvectorBlock(6),rtempVecFine%RvectorBlock(6))
!
!          ! Save that vector as new vector on the fine grid.
!          call sptivec_setTimestepData (rfineVector, 2*istep-1, rtempVecFine)
!        
!        end do
!        
!      else
!      
!        if (ispacelevelcoarse .ne. ispacelevelfine) then
!        
!          ! Prolongation only in space. 
!      
!          ! Loop through the time steps.
!          do istep = 1,rcoarseVector%NEQtime
!      
!            ! Space prolongation
!            call sptivec_getTimestepData (rcoarseVector, istep, rtempVecCoarse)
!            call mlprj_performProlongation (&
!                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!                rtempVecCoarse,rx3,rtempVecFineScalar)
!
!            ! Save that vector as new vector on the fine grid.
!            call sptivec_setTimestepData (rfineVector, istep, rx3)
!
!          end do
!
!        end if
!      
!      end if
!      
!      ! Release the temp vectors
!      call lsyssc_releaseVector (rtempVecFineScalar)
!      call lsysbl_releaseVector (rx3)
!      call lsysbl_releaseVector (rx1)
!
!    case (1)
!      ! Linear prolongation:
!      !
!      ! We need two more temp vectors:
!      !
!      ! One for the current timestep
!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.false.)
!      
!      ! And one for the next timestep
!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.false.)
!      
!      ! We need a scalar representation of the temp vector
!      call lsysbl_createScalarFromVec (rtempVecFine,rtempVecFineScalar)
!      
!      ! DEBUG!!!
!      call lsysbl_getbase_double (rx1,p_Dx1)
!      call lsysbl_getbase_double (rx3,p_Dx3)
!      call lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
!      call lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
!      call lsyssc_getbase_double (rtempVecFineScalar,p_DtempVecFineSca)
!
!      if (rcoarseVector%NEQtime .ne. rfineVector%NEQtime) then
!        
!        ! Prolongation in time.
!        !
!        ! Load timestep 0 into the temp vector and interpolate to the current level.
!        ! Put the result to rx3.
!        if (ispacelevelcoarse .ne. ispacelevelfine) then
!          ! Space + time
!          call sptivec_getTimestepData (rcoarseVector, 1+0, rtempVecCoarse)
!          call mlprj_performProlongation (&
!              rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!              rtempVecCoarse,rx3,rtempVecFineScalar)
!        else
!          ! Only time
!          call sptivec_getTimestepData (rcoarseVector, 1+0, rx3)
!        end if
!        
!        ! Save that vector as initial vector on the fine grid.
!        call sptivec_setTimestepData (rfineVector, 1+0, rx3)
!      
!        ! Loop through the time steps.
!        do istep = 1,rcoarseVector%NEQtime-1
!      
!          ! rx3 was the vector from the last timestep. Shift it to rx1 and load
!          ! the vector for the new timestep in rx3
!          call lsysbl_copyVector (rx3,rx1)
!
!          if (ispacelevelcoarse .ne. ispacelevelfine) then
!            ! Space + time
!            call sptivec_getTimestepData (rcoarseVector, 1+istep, rtempVecCoarse)
!            call mlprj_performProlongation (&
!                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!                rtempVecCoarse,rx3,rtempVecFineScalar)
!          else
!            ! Only time
!            call sptivec_getTimestepData (rcoarseVector, 1+istep, rx3)
!          end if
!          
!          ! Save that vector as new vector on the fine grid.
!          call sptivec_setTimestepData (rfineVector, 1+2*istep, rx3)
!          
!          ! In the temp vector, create the interpolation between rx1 and rx3.
!          ! THat's the prolongated vector.
!          call lsysbl_copyVector (rx1,rtempVecFine)
!          call lsysbl_vectorLinearComb (rx3,rtempVecFine,0.5_DP,0.5_DP)
!
!          ! Save that vector as new vector on the fine grid.
!          call sptivec_setTimestepData (rfineVector, 1+2*istep-1, rtempVecFine)
!        
!        end do
!        
!      else
!      
!        if (ispacelevelcoarse .ne. ispacelevelfine) then
!        
!          ! Prolongation only in space. 
!      
!          ! Loop through the time steps.
!          do istep = 1,rcoarseVector%NEQtime
!      
!            ! Space prolongation
!            call sptivec_getTimestepData (rcoarseVector, istep, rtempVecCoarse)
!            call mlprj_performProlongation (&
!                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!                rtempVecCoarse,rx3,rtempVecFineScalar)
!
!            ! Save that vector as new vector on the fine grid.
!            call sptivec_setTimestepData (rfineVector, istep, rx3)
!
!          end do
!          
!        else
!        
!          ! Copy the vector, it has the same size.
!          call sptivec_copyVector (rcoarseVector,rfineVector)
!
!        end if
!      
!      end if
!      
!      ! Release the temp vectors
!      call lsyssc_releaseVector (rtempVecFineScalar)
!      call lsysbl_releaseVector (rx3)
!      call lsysbl_releaseVector (rx1)
!      
!      ! DEBUG!!!
!      !call cc_postprocSpaceTimeGMV (rproblem,rdiscrFine,rfineVector,'gmv/fine.gmv')
!      !call cc_postprocSpaceTimeGMV (rproblem,rdiscrCoarse,rcoarseVector,'gmv/coarse.gmv')
!
!    case (2)
!      ! Linear prolongation respecting the Theta scheme.
!      !
!      ! This is a slightly more involved prolongation routine as the points in time
!      ! of the primal pressure and dual velocity might be off the timestepping.
!      ! Therefore, we have to interpolate / extrapolate to calculate some of the
!      ! vectors!
!      !
!      ! We need two more temp vectors:
!      !
!      ! One for the current timestep
!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.false.)
!      
!      ! And one for the next timestep
!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.false.)
!      
!      ! We need a scalar representation of the temp vector
!      call lsysbl_createScalarFromVec (rtempVecFine,rtempVecFineScalar)
!      
!      ! DEBUG!!!
!      call lsysbl_getbase_double (rx1,p_Dx1)
!      call lsysbl_getbase_double (rx3,p_Dx3)
!      call lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
!      call lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
!      call lsyssc_getbase_double (rtempVecFineScalar,p_DtempVecFineSca)
!
!      ! The fine grid solution corresponding to timestep 1 is created
!      ! by extrapolation and interpolation. The background is that
!      ! the dual velocity and primal pressure are located at different
!      ! points in time than the primal velocity and dual pressure.
!      ! For Crank-Nicolson, we have:
!      !
!      !
!      !     [1]         [-----2-----]           [-----3-----] |         [-----4-----]
!      !     X-----------o-----------X-----------o-----------X-----------o-----------X
!      !     y0                      y1                      y2                      y3
!      !     xi0                     xi1                     xi2                     xi3
!      !     l0          l1                      l2                      l3
!      !     p0          p1                      p2                      p3          
!      !
!      !     rx1         |-----rx3-----|         |-------------|         |-------------|
!      !
!      ! The prolongation will interpolate:
!      !
!      !     X-----------------------X-----------------------X-----------------------X
!      !     y0                      y1                      y2                      y3
!      !     xi0                     xi1                     xi2                     xi3
!      !        --1/2-->   <--1/2--     --1/2-->   <--1/2--     --1/2-->   <--1/2--
!      !     |                       |                       |                       |
!      !     v                       v                       v                       v
!      !     X-----------X-----------X-----------X-----------X-----------X-----------X
!      !     y0         y1           y2          y3          y4          y5          y6
!      !     xi0        xi1          xi2         xi3         xi4         xi5         xi6
!      !
!      ! The situation is more complicated for the primal pressure and dual velocity.
!      ! We use linear interpolation in the intervals and linear extrapolation at the
!      ! endpoints of the interval.
!      !
!      !     [1]         [-----2-----]           [-----3-----] |         [-----4-----]
!      !     l0          l1                      l2                      l3            
!      !     p0          p1                      p2                      p3             
!      !     X-----------o-----------X-----------o-----------X-----------o-----------X
!      !
!      !                 |-3/4-><------1/4-------|-3/4-><------1/4-------|
!      !                 |------1/4-------><-3/4-|------1/4-------><-3/4-|
!      !           <-5/4-|                                               |-5/4->
!      !           <---------- -1/4 -------------|---- -1/4 ------------------->
!      !
!      !     X-----o-----X-----o-----X-----o-----X-----o-----X-----o-----X-----o-----X
!      !     l0    l1          l2          l3          l4          l5          l6
!      !     p0    p1          p2          p3          p4          p5          p6
!      !
!      ! So at the end, we have:
!      !
!      !     [1]   [--2--]     [--3--]     [--4--]     [--5--]     [--6--]     [--7--]
!      !     X-----o-----X-----o-----X-----o-----X-----o-----X-----o-----X-----o-----X
!      !     y0          y1          y2          y3          y4          y5          y6
!      !     xi0         xi1         xi2         xi3         xi4         xi5         xi6
!      !     l0    l1          l2          l3          l4          l5          l6
!      !     p0    p1          p2          p3          p4          p5          p6
!      !
!      
!      if (rcoarseVector%NEQtime .ne. rfineVector%NEQtime) then
!        
!        ! Prolongation in time.
!
!        ! Load timestep 0 and 1 into the temp vector and interpolate to the current level.
!        ! Put the result to rx1 and rx3.
!        call getSpatialVector (rprojHier,rcoarseVector,1,rx1,&
!            rtempVecCoarse,rtempVecFineScalar,&
!            ispacelevelcoarse,ispacelevelfine)
!
!        call getSpatialVector (rprojHier,rcoarseVector,2,rx3,&
!            rtempVecCoarse,rtempVecFineScalar,&
!            ispacelevelcoarse,ispacelevelfine)
!        
!        ! Save rx1 as initial vector on the fine grid.
!        call sptivec_setTimestepData (rfineVector, 1, rx1)
!        
!        ! Prolongate y and xi.
!        call lsyssc_vectorLinearComb (rx1%RvectorBlock(1),rx3%RvectorBlock(1),&
!            0.5_DP,0.5_DP,rtempVecFine%RvectorBlock(1))
!        call lsyssc_vectorLinearComb (rx1%RvectorBlock(2),rx3%RvectorBlock(2),&
!            0.5_DP,0.5_DP,rtempVecFine%RvectorBlock(2))
!        call lsyssc_vectorLinearComb (rx1%RvectorBlock(6),rx3%RvectorBlock(6),&
!            0.5_DP,0.5_DP,rtempVecFine%RvectorBlock(6))
!        
!        ! Shift rx3 to rx1 and load the next solution.
!        call lsysbl_copyVector (rx3,rx1)
!        call getSpatialVector (rprojHier,rcoarseVector,3,rx3,&
!            rtempVecCoarse,rtempVecFineScalar,&
!            ispacelevelcoarse,ispacelevelfine)
!        
!        ! Create the first interpolated solution and save.
!        call lsyssc_vectorLinearComb (rx1%RvectorBlock(3),rx3%RvectorBlock(3),&
!            1.25_DP,-0.25_DP,rtempVecFine%RvectorBlock(3))
!        call lsyssc_vectorLinearComb (rx1%RvectorBlock(4),rx3%RvectorBlock(4),&
!            1.25_DP,-0.25_DP,rtempVecFine%RvectorBlock(4))
!        call lsyssc_vectorLinearComb (rx1%RvectorBlock(5),rx3%RvectorBlock(5),&
!            1.25_DP,-0.25_DP,rtempVecFine%RvectorBlock(5))
!            
!        call sptivec_setTimestepData (rfineVector, 2, rtempVecFine)
!
!        ! Loop through the time steps.
!        do istep = 3,rcoarseVector%NEQtime
!      
!          if (istep .gt. 3) then
!            ! rx3 was the vector from the last timestep. Shift it to rx1 and load
!            ! the vector for the new timestep in rx3.
!            call lsysbl_copyVector (rx3,rx1)
!
!            call getSpatialVector (rprojHier,rcoarseVector,1+istep,rx3,&
!                rtempVecCoarse,rtempVecFineScalar,&
!                ispacelevelcoarse,ispacelevelfine)
!          end if
!          
!          ! Create the y and xi. These can be copied from rx1.
!          call lsyssc_copyVector (rx1%RvectorBlock(1),rtempVecFine%RvectorBlock(1))
!          call lsyssc_copyVector (rx1%RvectorBlock(2),rtempVecFine%RvectorBlock(2))
!          call lsyssc_copyVector (rx1%RvectorBlock(6),rtempVecFine%RvectorBlock(6))
!          
!          ! Create the corresponding lambda and p by interpolation.
!          call lsyssc_vectorLinearComb (rx1%RvectorBlock(3),rtempVecFine%RvectorBlock(3),&
!              0.75_DP,0.25_DP)
!          call lsyssc_vectorLinearComb (rx1%RvectorBlock(4),rtempVecFine%RvectorBlock(4),&
!              0.75_DP,0.25_DP)
!          call lsyssc_vectorLinearComb (rx1%RvectorBlock(5),rtempVecFine%RvectorBlock(5),&
!              0.75_DP,0.25_DP)
!          
!          ! Save that vector as new vector on the fine grid.
!          call sptivec_setTimestepData (rfineVector, 2*(istep-1)-1, rtempVecFine)
!          
!          ! Create the next y and xi by interpolation.
!          call lsyssc_vectorLinearComb (rx1%RvectorBlock(1),rx3%RvectorBlock(1),&
!              0.5_DP,0.5_DP,rtempVecFine%RvectorBlock(1))
!          call lsyssc_vectorLinearComb (rx1%RvectorBlock(2),rx3%RvectorBlock(2),&
!              0.5_DP,0.5_DP,rtempVecFine%RvectorBlock(2))
!          call lsyssc_vectorLinearComb (rx1%RvectorBlock(6),rx3%RvectorBlock(6),&
!              0.5_DP,0.5_DP,rtempVecFine%RvectorBlock(6))
!          
!          ! And the corresponding lambda and p.
!          call lsyssc_vectorLinearComb (rx1%RvectorBlock(3),rx3%RvectorBlock(3),&
!              0.25_DP,0.75_DP,rtempVecFine%RvectorBlock(3))                               
!          call lsyssc_vectorLinearComb (rx1%RvectorBlock(4),rx3%RvectorBlock(4),&
!              0.25_DP,0.75_DP,rtempVecFine%RvectorBlock(4))                               
!          call lsyssc_vectorLinearComb (rx1%RvectorBlock(5),rx3%RvectorBlock(5),&
!              0.25_DP,0.75_DP,rtempVecFine%RvectorBlock(5))
!              
!          ! Save it.
!          call sptivec_setTimestepData (rfineVector, 2*(istep-1), rtempVecFine)
!          
!        end do
!        
!        ! The vector at the end is missing. To create it, we take y and xi from
!        ! the end and extrapolate lambda and p.
!        call lsyssc_copyVector (rx3%RvectorBlock(1),rtempVecFine%RvectorBlock(1))
!        call lsyssc_copyVector (rx3%RvectorBlock(2),rtempVecFine%RvectorBlock(2))
!        call lsyssc_copyVector (rx3%RvectorBlock(6),rtempVecFine%RvectorBlock(6))
!        
!        call lsyssc_vectorLinearComb (rx1%RvectorBlock(3),rx3%RvectorBlock(3),&
!            -0.25_DP,1.25_DP,rtempVecFine%RvectorBlock(3))
!        call lsyssc_vectorLinearComb (rx1%RvectorBlock(4),rx3%RvectorBlock(4),&
!            -0.25_DP,1.25_DP,rtempVecFine%RvectorBlock(4))
!        call lsyssc_vectorLinearComb (rx1%RvectorBlock(5),rx3%RvectorBlock(5),&
!            -0.25_DP,1.25_DP,rtempVecFine%RvectorBlock(5))
!            
!        call sptivec_setTimestepData (rfineVector, 2*rcoarseVector%NEQtime-1, rtempVecFine)
!        
!      else
!      
!        if (ispacelevelcoarse .ne. ispacelevelfine) then
!        
!          ! Prolongation only in space. 
!      
!          ! Loop through the time steps.
!          do istep = 1,rcoarseVector%NEQtime
!      
!            ! Space prolongation
!            call sptivec_getTimestepData (rcoarseVector, istep, rtempVecCoarse)
!            call mlprj_performProlongation (&
!                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!                rtempVecCoarse,rx3,rtempVecFineScalar)
!
!            ! Save that vector as new vector on the fine grid.
!            call sptivec_setTimestepData (rfineVector, istep, rx3)
!
!          end do
!          
!        else
!        
!          ! Copy the vector, it has the same size.
!          call sptivec_copyVector (rcoarseVector,rfineVector)
!
!        end if
!      
!      end if
!      
!      ! Release the temp vectors
!      call lsyssc_releaseVector (rtempVecFineScalar)
!      call lsysbl_releaseVector (rx3)
!      call lsysbl_releaseVector (rx1)
!      
!      ! DEBUG!!!
!      !call cc_postprocSpaceTimeGMV (rproblem,rdiscrFine,rfineVector,'gmv/fine.gmv')
!      !call cc_postprocSpaceTimeGMV (rproblem,rdiscrCoarse,rcoarseVector,'gmv/coarse.gmv')
!
!    case (3)
!      ! Quadratic prolongation.
!      !
!      ! We need vectors from three timesteps.
!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.false.)
!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx2,.false.)
!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.false.)
!      
!      ! We need a scalar representation of the temp vector
!      call lsysbl_createScalarFromVec (rtempVecFine,rtempVecFineScalar)
!      call lsysbl_createScalarFromVec (rtempVecFine,rtempVecFineScalar2)
!      
!      ! DEBUG!!!
!      call lsysbl_getbase_double (rx1,p_Dx1)
!      call lsysbl_getbase_double (rx3,p_Dx3)
!      call lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
!      call lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
!      call lsyssc_getbase_double (rtempVecFineScalar,p_DtempVecFineSca)
!
!      ! Load timestep 0 into the temp vector and interpolate to the current level.
!      ! Put the result to rx3.
!      if (ispacelevelcoarse .ne. ispacelevelfine) then
!        ! Space + time
!        call sptivec_getTimestepData (rcoarseVector, 1+0, rtempVecCoarse)
!        call mlprj_performProlongation (&
!            rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!            rtempVecCoarse,rx3,rtempVecFineScalar)
!      else
!        ! Only time
!        call sptivec_getTimestepData (rcoarseVector, 1+0, rx3)
!      end if
!      
!      ! Save vector rx1 as initial vector on the fine grid.
!      call sptivec_setTimestepData (rfineVector, 1+0, rx3)
!      
!      if (rcoarseVector%NEQtime .ne. rfineVector%NEQtime) then
!        
!        ! Quadratic Prolongation in time.
!        !
!        !                v------- -0.125 -----------------+
!        !     +------------------ -0.125 ------v
!        !                v-- 0.75 --+-- 0.75 --v
!        !     +- 0.375 --v                     v-- 0.375 -+
!        !
!        !     X----------o----------X----------o----------X
!        !     rx1                  rx2                  rx3
!        !
!        ! Loop through the time steps. We loop in sets of size 2 and
!        ! apply quadratic interpolation in each set.
!        do istep = 1,rcoarseVector%NEQtime-2,2
!      
!          ! rx3 was the vector from the last timestep. Shift it to rx1 and load
!          ! the vector for the new timestep patch into rx3, the data from the
!          ! middle of the timestep patch to rx2.
!          call lsysbl_copyVector (rx3,rx1)
!
!          if (ispacelevelcoarse .ne. ispacelevelfine) then
!            ! Space + time
!            call sptivec_getTimestepData (rcoarseVector, 1+istep+1, rtempVecCoarse)
!            call mlprj_performProlongation (&
!                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!                rtempVecCoarse,rx3,rtempVecFineScalar)
!
!            call sptivec_getTimestepData (rcoarseVector, 1+istep, rtempVecCoarse)
!            call mlprj_performProlongation (&
!                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!                rtempVecCoarse,rx2,rtempVecFineScalar)
!          else
!            ! Only time
!            call sptivec_getTimestepData (rcoarseVector, 1+istep+1, rx3)
!            call sptivec_getTimestepData (rcoarseVector, 1+istep, rx2)
!          end if
!          
!          ! Save the vector as new vectors on the fine grid.
!          call sptivec_setTimestepData (rfineVector, 1+2*istep, rx2)
!          call sptivec_setTimestepData (rfineVector, 1+2*istep+2, rx3)
!          
!          ! In the temp vector, create the interpolation between rx1,rx2 and rx3.
!          ! That's the prolongated vector.
!          ! Save that vector as new vector on the fine grid.
!          call lsysbl_copyVector (rx2,rtempVecFine)
!          call lsysbl_vectorLinearComb (rx1,rtempVecFine,0.375_DP,0.75_DP)
!          call lsysbl_vectorLinearComb (rx3,rtempVecFine,-0.125_DP,1.0_DP)
!          call sptivec_setTimestepData (rfineVector, 1+2*istep-1, rtempVecFine)
!
!          call lsysbl_copyVector (rx2,rtempVecFine)
!          call lsysbl_vectorLinearComb (rx1,rtempVecFine,-0.125_DP,0.75_DP)
!          call lsysbl_vectorLinearComb (rx3,rtempVecFine,0.375_DP,1.0_DP)
!          call sptivec_setTimestepData (rfineVector, 1+2*istep+1, rtempVecFine)
!        
!        end do
!        
!        if (mod(rcoarseVector%NEQtime,2) .eq. 0) then
!        
!          ! Number of time solutions is equal, that means we have an odd number
!          ! of timesteps -- and so we have to process the last timestep which
!          ! could not be reached by the above loop as it's only half a patch.
!        
!          ! We proceed only half a patch and create a new patch from the
!          ! time solutions NEQTIME, NEQTIME-1 and NEQTIME-2.
!          ! The solution at NEQTIME-2 is currently present in rx2,
!          ! the solution at NEQTIME-1 in rx3.
!          call lsysbl_copyVector (rx2,rx1)
!          call lsysbl_copyVector (rx3,rx2)
!
!          if (ispacelevelcoarse .ne. ispacelevelfine) then
!            ! Space + time
!            call sptivec_getTimestepData (rcoarseVector, rcoarseVector%NEQtime, rtempVecCoarse)
!            call mlprj_performProlongation (&
!                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!                rtempVecCoarse,rx3,rtempVecFineScalar)
!          else
!            ! Only time
!            call sptivec_getTimestepData (rcoarseVector, rcoarseVector%NEQtime, rx3)
!          end if
!          
!          ! Save the last vector
!          call sptivec_setTimestepData (rfineVector, rfineVector%NEQtime, rx3)
!
!          ! Save the quadratic interpolation.          
!          call lsysbl_copyVector (rx2,rtempVecFine)
!          call lsysbl_vectorLinearComb (rx1,rtempVecFine,-0.125_DP,0.75_DP)
!          call lsysbl_vectorLinearComb (rx3,rtempVecFine,0.375_DP,1.0_DP)
!          call sptivec_setTimestepData (rfineVector, rfineVector%NEQtime-1, rtempVecFine)
!        
!        end if
!
!      else
!      
!        if (ispacelevelcoarse .ne. ispacelevelfine) then
!        
!          ! Prolongation only in space. 
!      
!          ! Loop through the time steps.
!          do istep = 1,rcoarseVector%NEQtime-1
!      
!            ! Space prolongation
!            call sptivec_getTimestepData (rcoarseVector, 1+istep, rtempVecCoarse)
!            call mlprj_performProlongation (&
!                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!                rtempVecCoarse,rx3,rtempVecFineScalar)
!
!            ! Save that vector as new vector on the fine grid.
!            call sptivec_setTimestepData (rfineVector, 1+istep, rx3)
!
!          end do
!
!        end if
!      
!      end if
!      
!      ! Release the temp vectors
!      call lsyssc_releaseVector (rtempVecFineScalar2)
!      call lsyssc_releaseVector (rtempVecFineScalar)
!      call lsysbl_releaseVector (rx3)
!      call lsysbl_releaseVector (rx2)
!      call lsysbl_releaseVector (rx1)
!
!    end select
!
!    call sptivec_saveToFileSequence (rfineVector,"(""fineold.txt."",I5.5)",.true.,&
!        rtempVecFine)

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

    ! Local variables
    type(t_spaceTimeVectorAccess) :: raccessPool
    integer :: istep,ispacelevelfine,ispacelevelcoarse,itimelevelfine,itimelevelcoarse
    type(t_vectorScalar) :: rtempVecFineScalar,rtempVecFineScalar2
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
    select case (rprojHier%itimeOrder)
    case (0)
      call sptivec_createAccessPool (rfineVector,raccessPool,3)
    case (1)
      call sptivec_createAccessPool (rfineVector,raccessPool,3)
    case (2)
      call sptivec_createAccessPool (rfineVector,raccessPool,7)
    case (3)
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
    
    ! Apply the multiplication.
    ! The rows in the matrix correspond to the time fine mesh, the columns
    ! to the time coarse mesh.
    do irow = 1,p_rrestMatrixPrim%NEQ
    
      ! Clear the destination
      call lsysbl_clearVector (rtempVecFine)
      
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
            p_rx%RvectorBlock(3),rtempVecFine%RvectorBlock(3),dscaleDual*p_DaDual(icol),1.0_DP)
            
        call lsyssc_vectorLinearComb(&
            p_rx%RvectorBlock(4),rtempVecFine%RvectorBlock(4),dscaleDual*p_DaDual(icol),1.0_DP)

        call lsyssc_vectorLinearComb(&
            p_rx%RvectorBlock(5),rtempVecFine%RvectorBlock(5),dscaleDual*p_DaDual(icol),1.0_DP)
        
      end do
      
      ! Vector finished.
      call setSpatialVector (rprojHier,rtempVecFine,rcoarseVector,irow,&
          rtempVecCoarse,rtempVecFineScalar,ispacelevelcoarse,ispacelevelfine)

    end do
    
    ! Release the buffer.
    call sptivec_releaseAccessPool(raccessPool)
    
    call lsyssc_releaseVector (rtempVecFineScalar)

    ! DEBUG!!!
!    call sptivec_saveToFileSequence (rcoarseVector,"(""coarse.txt."",I5.5)",.true.,&
!        rtempVecFine)
!    call sys_halt()
        
! ********************************************************************************************

!    ! local variables
!    integer :: istep,ispacelevelfine,ispacelevelcoarse,iorder
!    type(t_vectorBlock) :: rx1,rx2,rx3,rx4,rx5
!    type(t_vectorScalar) :: rx1scalar
!    type(t_vectorBlock) :: rx1,rx2,rx3,rx4,rx5
!    type(t_vectorScalar) :: rx1scalar
!    real(DP), dimension(:), pointer :: p_Dx1,p_Dx3
!    real(DP), dimension(:), pointer :: p_DtempVecCoarse,p_DtempVecFine
!    
!    ! DEBUG!!!
!    real(DP), dimension(:), pointer :: p_Dx1,p_Dx3
!    real(DP), dimension(:), pointer :: p_DtempVecCoarse,p_DtempVecFine
!
!    ! Get the space-time coarse and fine discretisations.
!    call sth_getLevel (rprojHier%p_rspaceTimeHierarchy,ilevelfine-1,&
!      ispaceLevel=ispacelevelcoarse)
!    call sth_getLevel (rprojHier%p_rspaceTimeHierarchy,ilevelfine,&
!      ispaceLevel=ispacelevelfine)
!
!!    ! Select the prolongation/restriction to use
!!    iorder = rprojHier%itimeOrder
!!    if (iorder .eq. -1) then
!!      ! automatic mode. Select the order based on the time stepping scheme.
!!      call tdiscr_getOrder(rcoarseVector%p_rtimeDiscr,iorder)
!!    end if
!
!    select case (rprojHier%itimeOrder)
!    case (0)
!      ! Constant restriction:
!      !
!      ! Prolongation 'distributes' information from the coarse grid nodes
!      ! to the fine grid nodes as follows:
!      !
!      ! Timestep:    n-1         n          n+1        (fine grid)
!      !               x <--1/2-- X --1/2--> x
!      !                         ^ \
!      !                        /   \
!      !                       +--1--+
!      !
!      ! Restriction is the 'adjoint' of the prolongation, so the corresponding
!      ! matrix is the transposed. Written nodewise this means, that restriction
!      ! has to 'collect' the distributed information with the same weights:
!      !
!      ! Timestep:    n-1         n          n+1        (fine grid)
!      !               x --1/2--> X <--1/2-- x
!      !                         ^ \
!      !                        /   \
!      !                       +--1--+
!      !
!      ! But because this is a restriction of a Finite-Difference RHS vector,
!      ! the RHS must be divided by h/(h/2)=2 !
!
!      ! One for the current timestep
!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.false.)
!      
!      ! And one for the next timestep
!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.false.)
!      
!      ! Create a scalar representation of rx1 for auxiliary reasons
!      call lsysbl_createScalarFromVec (rx1,rx1Scalar)
!     
!      ! DEBUG!!!
!      call lsysbl_getbase_double (rx1,p_Dx1)
!      call lsysbl_getbase_double (rx3,p_Dx3)
!      call lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
!      call lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
!         
!      if (rcoarseVector%NEQtime .NE. rfineVector%NEQtime) then       
!      
!        ! Load timestep 0,1 (fine grid) into the temp vectors.
!        call sptivec_getTimestepData (rfineVector, 0, rtempVecFine)
!        call sptivec_getTimestepData (rfineVector, 1, rx3)
!        
!        ! Perform restriction for the first time step.
!        !                          X <--1/2-- x
!        !                         ^ \
!        !                        /   \
!        !                       +--1--+
!        call lsyssc_vectorLinearComb (rx3%RvectorBlock(4),rtempVecFine%RvectorBlock(4),&
!            0.5_DP/2.0_DP,1.0_DP/2.0_DP)
!        call lsyssc_vectorLinearComb (rx3%RvectorBlock(5),rtempVecFine%RvectorBlock(5),&
!            0.5_DP/2.0_DP,1.0_DP/2.0_DP)
!        call lsyssc_vectorLinearComb (rx3%RvectorBlock(6),rtempVecFine%RvectorBlock(6),&
!            0.5_DP/2.0_DP,1.0_DP/2.0_DP)
!        
!        ! Probably restrict the vector in space to the lower level.
!        ! Save the result.
!        if (ispacelevelcoarse .ne. ispacelevelfine) then
!          ! Space + time
!          call mlprj_performRestriction (&
!                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!                rtempVecCoarse,rtempVecFine,rx1Scalar)
!          call sptivec_setTimestepData (rcoarseVector, 0, rtempVecCoarse)
!        else
!          ! Only time
!          call sptivec_setTimestepData (rcoarseVector, 0, rtempVecFine)
!        end if
!        
!        ! Loop through the time steps -- on the coarse grid!
!        do istep = 1,rcoarseVector%NEQtime-1
!        
!          ! rx3 was the 'right inner fine grid node x' in the above picture.
!          ! Copy it to rx1, so it gets the new 'left inner fine grid node x'
!          ! and load the new 'right inner fine grid node x' to rx1
!          ! as well as the 'inner fine grid node X' to rtempVecFine.
!          call lsysbl_copyVector (rx3,rx1)
!          
!          call sptivec_getTimestepData (rfineVector,2*istep, rtempVecFine)
!          call sptivec_getTimestepData (rfineVector,2*istep+1, rx3)
!          
!          ! In rtempVecFine, create the restriction of the fine grid vectors
!          ! rx1 and rx3 according to
!          !
!          ! Timestep:    n-1          n          n+1        (fine grid)
!          !              rx3 --1/2--> X <--1/2-- rx3
!          !             (old)        ^ \
!          !             =rx1        /   \
!          !                        +--1--+
!          
!          call lsyssc_vectorLinearComb (rx1%RvectorBlock(1),rtempVecFine%RvectorBlock(1),&
!              0.5_DP,1.0_DP/2.0_DP)
!          call lsyssc_vectorLinearComb (rx1%RvectorBlock(2),rtempVecFine%RvectorBlock(2),&
!              0.5_DP,1.0_DP/2.0_DP)
!          call lsyssc_vectorLinearComb (rx1%RvectorBlock(3),rtempVecFine%RvectorBlock(3),&
!              0.5_DP,1.0_DP/2.0_DP)
!
!          call lsyssc_vectorLinearComb (rx3%RvectorBlock(4),rtempVecFine%RvectorBlock(4),&
!              0.5_DP,1.0_DP/2.0_DP)
!          call lsyssc_vectorLinearComb (rx3%RvectorBlock(5),rtempVecFine%RvectorBlock(5),&
!              0.5_DP,1.0_DP/2.0_DP)
!          call lsyssc_vectorLinearComb (rx3%RvectorBlock(6),rtempVecFine%RvectorBlock(6),&
!              0.5_DP,1.0_DP/2.0_DP)
!
!          ! Probably restrict the vector in space to the lower level.
!          ! Save the result.
!          if (ispacelevelcoarse .ne. ispacelevelfine) then
!            ! Space + time
!            call mlprj_performRestriction (&
!                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!                rtempVecCoarse,rtempVecFine,rx1Scalar)
!            call sptivec_setTimestepData (rcoarseVector, istep, rtempVecCoarse)
!          else
!            ! Only time
!            call sptivec_setTimestepData (rcoarseVector, istep, rtempVecFine)
!          end if
!        
!        end do
!
!        ! Last time step.
!
!        ! rx3 is now the new 'left inner fine grid node x'.
!        ! Load the 'inner fine grid node X' to rtempVecFine.
!        call sptivec_getTimestepData (rfineVector,2*rcoarseVector%NEQtime, rtempVecFine)
!        
!        ! In rtempVecFine, create the restriction of the fine grid vectors
!        ! rx1 and rx3 according to
!        !
!        ! Timestep:    n-1          n       (fine grid)
!        !              rx3 --1/2--> X             
!        !                          ^ \
!        !                         /   \
!        !                        +--1--+
!        
!        call lsyssc_vectorLinearComb (rx3%RvectorBlock(1),rtempVecFine%RvectorBlock(1),&
!            0.5_DP,1.0_DP/2.0_DP)
!        call lsyssc_vectorLinearComb (rx3%RvectorBlock(2),rtempVecFine%RvectorBlock(2),&
!            0.5_DP,1.0_DP/2.0_DP)
!        call lsyssc_vectorLinearComb (rx3%RvectorBlock(3),rtempVecFine%RvectorBlock(3),&
!            0.5_DP,1.0_DP/2.0_DP)
!        
!        ! Probably restrict the vector in space to the lower level.
!        ! Save the result.
!        if (ispacelevelcoarse .ne. ispacelevelfine) then
!          ! Space + time
!          call mlprj_performRestriction (&
!              rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!              rtempVecCoarse, rtempVecFine,rx1Scalar)
!          call sptivec_setTimestepData (rcoarseVector, &
!              rcoarseVector%NEQtime, rtempVecCoarse)
!        else
!          ! Only time
!          call sptivec_setTimestepData (rcoarseVector, &
!              rcoarseVector%NEQtime, rtempVecFine)
!        end if
!        
!      else
!      
!        if (ispacelevelcoarse .ne. ispacelevelfine) then
!        
!          ! Restriction only in space.
!          !
!          ! Loop through the time steps
!          do istep = 0,rcoarseVector%NEQtime
!          
!            ! Load the data
!            call sptivec_getTimestepData (rfineVector,istep, rtempVecFine)
!
!            ! Process the restriction and save the data to the coarse vector.
!            call mlprj_performRestriction (&
!                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!                rtempVecCoarse, rtempVecFine,rx1Scalar)
!                
!            call sptivec_setTimestepData (rcoarseVector, istep, rtempVecCoarse)
!          
!          end do
!        
!        end if
!      
!      end if
!
!      ! Release the temp vectors
!      call lsyssc_releaseVector (rx1Scalar)
!      call lsysbl_releaseVector (rx3)
!      call lsysbl_releaseVector (rx1)
!
!    case (1)
!      ! Linear restriction:
!      !
!      ! Prolongation 'distributes' information from the coarse grid nodes
!      ! to the fine grid nodes as follows:
!      !
!      ! Timestep:    n-1         n          n+1        (fine grid)
!      !               x <--1/2-- X --1/2--> x
!      !                         ^ \
!      !                        /   \
!      !                       +--1--+
!      !
!      ! Restriction is the 'adjoint' of the prolongation, so the corresponding
!      ! matrix is the transposed. Written nodewise this means, that restriction
!      ! has to 'collect' the distributed information with the same weights:
!      !
!      ! Timestep:    n-1         n          n+1        (fine grid)
!      !               x --1/2--> X <--1/2-- x
!      !                         ^ \
!      !                        /   \
!      !                       +--1--+
!      !
!      ! But because this is a restriction of a Finite-Difference RHS vector,
!      ! the RHS must be divided by h/(h/2)=2 !
!
!      ! We need two more temp vectors:
!      !
!      ! One for the current timestep
!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.false.)
!      
!      ! And one for the next timestep
!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.false.)
!      
!      ! Create a scalar representation of rx1 for auxiliary reasons
!      call lsysbl_createScalarFromVec (rx1,rx1Scalar)
!     
!      ! DEBUG!!!
!      call lsysbl_getbase_double (rx1,p_Dx1)
!      call lsysbl_getbase_double (rx3,p_Dx3)
!      call lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
!      call lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
!         
!      if (rcoarseVector%NEQtime .ne. rfineVector%NEQtime) then       
!      
!        ! Load timestep 0,1 (fine grid) into the temp vectors.
!        call sptivec_getTimestepData (rfineVector, 1+0, rtempVecFine)
!        call sptivec_getTimestepData (rfineVector, 1+1, rx3)
!        
!        ! Perform restriction for the first time step.
!        !                          X <--1/2-- x
!        !                         ^ \
!        !                        /   \
!        !                       +--1--+
!        call lsysbl_vectorLinearComb (rx3,rtempVecFine,0.5_DP/2.0_DP,1.0_DP/2.0_DP)
!        
!        ! Probably restrict the vector in space to the lower level.
!        ! Save the result.
!        if (ispacelevelcoarse .ne. ispacelevelfine) then
!          ! Space + time
!          call mlprj_performRestriction (&
!              rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!              rtempVecCoarse,rtempVecFine,rx1Scalar)
!          call sptivec_setTimestepData (rcoarseVector, 1+0, rtempVecCoarse)
!        else
!          ! Only time
!          call sptivec_setTimestepData (rcoarseVector, 1+0, rtempVecFine)
!        end if
!        
!        ! Loop through the time steps -- on the coarse grid!
!        do istep = 1,rcoarseVector%NEQtime-1-1
!        
!          ! rx3 was the 'right inner fine grid node x' in the above picture.
!          ! Copy it to rx1, so it gets the new 'left inner fine grid node x'
!          ! and load the new 'right inner fine grid node x' to rx1
!          ! as well as the 'inner fine grid node X' to rtempVecFine.
!          call lsysbl_copyVector (rx3,rx1)
!          
!          call sptivec_getTimestepData (rfineVector,1+2*istep, rtempVecFine)
!          call sptivec_getTimestepData (rfineVector,1+2*istep+1, rx3)
!          
!          ! In rtempVecFine, create the restriction of the fine grid vectors
!          ! rx1 and rx3 according to
!          !
!          ! Timestep:    n-1          n          n+1        (fine grid)
!          !              rx3 --1/2--> X <--1/2-- rx3
!          !             (old)        ^ \
!          !                         /   \
!          !                        +--1--+
!          
!          call lsysbl_vectorLinearComb (rx1,rtempVecFine,0.5_DP/2.0_DP,1.0_DP/2.0_DP)
!          call lsysbl_vectorLinearComb (rx3,rtempVecFine,0.5_DP/2.0_DP,1.0_DP)
!          
!          ! Probably restrict the vector in space to the lower level.
!          ! Save the result.
!          if (ispacelevelcoarse .ne. ispacelevelfine) then
!            ! Space + time
!            call mlprj_performRestriction (&
!                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!                rtempVecCoarse,rtempVecFine,rx1Scalar)
!            call sptivec_setTimestepData (rcoarseVector, 1+istep, rtempVecCoarse)
!          else
!            ! Only time
!            call sptivec_setTimestepData (rcoarseVector, 1+istep, rtempVecFine)
!          end if
!        
!        end do
!
!        ! Last time step.
!
!        ! rx3 is now the new 'left inner fine grid node x'.
!        ! Load the 'inner fine grid node X' to rtempVecFine.
!        call sptivec_getTimestepData (rfineVector,1+2*(rcoarseVector%NEQtime-1), &
!          rtempVecFine)
!        
!        ! In rtempVecFine, create the restriction of the fine grid vectors
!        ! rx1 and rx3 according to
!        !
!        ! Timestep:    n-1          n       (fine grid)
!        !              rx3 --1/2--> X             
!        !                          ^ \
!        !                         /   \
!        !                        +--1--+
!        
!        call lsysbl_vectorLinearComb (rx3,rtempVecFine,0.5_DP/2.0_DP,1.0_DP/2.0_DP)
!        
!        ! Probably restrict the vector in space to the lower level.
!        ! Save the result.
!        if (ispacelevelcoarse .ne. ispacelevelfine) then
!          ! Space + time
!          call mlprj_performRestriction (&
!              rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!              rtempVecCoarse, rtempVecFine,rx1Scalar)
!          call sptivec_setTimestepData (rcoarseVector, &
!              rcoarseVector%NEQtime, rtempVecCoarse)
!        else
!          ! Only time
!          call sptivec_setTimestepData (rcoarseVector, &
!              rcoarseVector%NEQtime, rtempVecFine)
!        end if
!        
!      else
!      
!        if (ispacelevelcoarse .ne. ispacelevelfine) then
!        
!          ! Restriction only in space.
!          !
!          ! Loop through the time steps
!          do istep = 0,rcoarseVector%NEQtime-1
!          
!            ! Load the data
!            call sptivec_getTimestepData (rfineVector,1+istep, rtempVecFine)
!
!            ! Process the restriction and save the data to the coarse vector.
!            call mlprj_performRestriction (&
!                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!                rtempVecCoarse, rtempVecFine,rx1Scalar)
!                
!            call sptivec_setTimestepData (rcoarseVector,1+istep, rtempVecCoarse)
!          
!          end do
!        
!        end if
!      
!      end if
!
!      ! Release the temp vectors
!      call lsyssc_releaseVector (rx1Scalar)
!      call lsysbl_releaseVector (rx3)
!      call lsysbl_releaseVector (rx1)
!
!    case (2)
!!      ! Linear restriction with support for the general Theta scheme.
!!      !
!!      ! We need some temporary vectors:
!!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.false.)
!!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx2,.false.)
!!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.false.)
!!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx4,.false.)
!!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx5,.false.)
!!
!!      ! One for the current timestep
!!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.false.)
!!      
!!      ! And one for the next timestep
!!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.false.)
!!      
!!      ! Create a scalar representation of rx1 for auxiliary reasons
!!      call lsysbl_createScalarFromVec (rx1,rx1Scalar)
!!     
!!      ! DEBUG!!!
!!      call lsysbl_getbase_double (rx1,p_Dx1)
!!      call lsysbl_getbase_double (rx3,p_Dx3)
!!      call lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
!!      call lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
!!         
!!      if (rcoarseVector%NEQtime .ne. rfineVector%NEQtime) then       
!!      
!!        ! The restriction is the adjoint of the prolongation.
!!        ! For Crank-Nicolson, we have:
!!        !
!!        !     [1]   [--2--]     [--3--]     [--4--]     [--5--]     [--6--]     [--7--]
!!        !     X-----o-----X-----o-----X-----o-----X-----o-----X-----o-----X-----o-----X
!!        !     y0          y1          y2          y3          y4          y5          y6
!!        !     xi0         xi1         xi2         xi3         xi4         xi5         xi6
!!        !     l0    l1          l2          l3          l4          l5          l6
!!        !     p0    p1          p2          p3          p4          p5          p6
!!        !
!!        ! The restriction will to the weighten mean for y and xi:
!!        !
!!        !     X-----------X-----------X-----------X-----------X-----------X-----------X
!!        !     y0         y1           y2          y3          y4          y5          y6
!!        !     xi0        xi1          xi2         xi3         xi4         xi5         xi6
!!        !        <-1/2---   ---1/2-->   <--1/2---   ---1/2-->   <--1/2---   ---1/2-->
!!        !     |                       |                       |                       |
!!        !     v                       v                       v                       v
!!        !     X-----------------------X-----------------------X-----------------------X
!!        !     y0                      y1                      y2                      y3
!!        !     xi0                     xi1                     xi2                     xi3
!!        !
!!        ! The situation is more complicated for the primal pressure and dual velocity.
!!        ! We use linear interpolation in the intervals and linear extrapolation at the
!!        ! endpoints of the interval.
!!        !
!!        !     [1]   [--2--]     [--3--]     [--4--]     [--5--]     [--6--]     [--7--]
!!        !     l0    l1          l2          l3          l4          l5          l6
!!        !     p0    p1          p2          p3          p4          p5          p6
!!        !     X-----o-----X-----o-----X-----o-----X-----o-----X-----o-----X-----o-----X
!!        !
!!        !                 <-3/4-|------1/4-------><-3/4-|-------1/4------->
!!        !                 <------1/4--------|-3/4><------1/4--------|-3/4->
!!        !           |-5/4->                                               <-5/4-|
!!        !           |---------- -1/4 ------------><---- -1/4 -------------------|
!!        !
!!        !     X-----------o-----------X-----------o-----------X-----------o-----------X
!!        !     l0          l1                      l2                      l3
!!        !     p0          p1                      p2                      p3          
!!        !
!!        ! So at the end, we have:
!!        !
!!        !     [1]         [-----2-----]           [-----3-----] |         [-----4-----]
!!        !     X-----------o-----------X-----------o-----------X-----------o-----------X
!!        !     y0                      y1                      y2                      y3
!!        !     xi0                     xi1                     xi2                     xi3
!!        !     l0          l1                      l2                      l3
!!        !     p0          p1                      p2                      p3          
!!
!!        ! Load the vector 1..5. We need them for the restriction at the start.
!!        ! The initial condition is taken as it is.
!!        call sptivec_getTimestepData (rfineVector, 1, rx1)
!!        call sptivec_getTimestepData (rfineVector, 2, rx2)
!!        call sptivec_getTimestepData (rfineVector, 3, rx3)
!!        call sptivec_getTimestepData (rfineVector, 4, rx4)
!!        call sptivec_getTimestepData (rfineVector, 5, rx5)
!!
!!        ! Restrict the initial condition.
!!        ! There is a contribtion for y and xi, not for lambda and xi!
!!        call lsysbl_copyVector (rx1,rtempVecFine)
!!        call lsyssc_vectorLinearComb  (
!!            rx2%RvectorBlock(1),rtempVecFine%RvectorBlock(1),0.5_DP,1.0_DP)
!!        call lsyssc_vectorLinearComb  (
!!            rx2%RvectorBlock(2),rtempVecFine%RvectorBlock(2),0.5_DP,1.0_DP)
!!        call lsyssc_vectorLinearComb  (
!!            rx2%RvectorBlock(6),rtempVecFine%RvectorBlock(6),0.5_DP,1.0_DP)
!!        
!!        call setSpatialVector (rprojHier,rx1,rcoarseVector,1,&
!!            rtempVecCoarse,rtempVecFineScalar,ispacelevelcoarse,ispacelevelfine)
!!
!!        ! Now restrict the first two intervals.
!!        !
!!        ! First y/xi.
!!        call lsyssc_copyVector (rx3%RvectorBlock(1),rtempVecFine%RvectorBlock(1))
!!        call lsyssc_vectorLinearComb  (
!!            rx2%RvectorBlock(1),rtempVecFine%RvectorBlock(1),0.5_DP,1.0_DP)
!!        call lsyssc_vectorLinearComb  (
!!            rx4%RvectorBlock(1),rtempVecFine%RvectorBlock(1),0.5_DP,1.0_DP)
!!        
!!        call lsyssc_copyVector (rx3%RvectorBlock(2),rtempVecFine%RvectorBlock(2))
!!        call lsyssc_vectorLinearComb  (
!!            rx2%RvectorBlock(2),rtempVecFine%RvectorBlock(2),0.5_DP,1.0_DP)
!!        call lsyssc_vectorLinearComb  (
!!            rx4%RvectorBlock(2),rtempVecFine%RvectorBlock(2),0.5_DP,1.0_DP)
!!
!!        call lsyssc_copyVector (rx3%RvectorBlock(6),rtempVecFine%RvectorBlock(6))
!!        call lsyssc_vectorLinearComb  (
!!            rx2%RvectorBlock(6),rtempVecFine%RvectorBlock(6),0.5_DP,1.0_DP)
!!        call lsyssc_vectorLinearComb  (
!!            rx4%RvectorBlock(6),rtempVecFine%RvectorBlock(6),0.5_DP,1.0_DP)
!!        
!!        ! lambda/p is more involved.
!!        call lsyssc_vectorLinearComb  (
!!            rx2%RvectorBlock(3),rtempVecFine%RvectorBlock(3),1.25_DP,0.0_DP)
!!        call lsyssc_vectorLinearComb  (
!!            rx3%RvectorBlock(3),rtempVecFine%RvectorBlock(3),0.75_DP,1.0_DP)
!!        call lsyssc_vectorLinearComb  (
!!            rx4%RvectorBlock(3),rtempVecFine%RvectorBlock(3),0.25_DP,1.0_DP)
!!
!!        call lsyssc_vectorLinearComb  (
!!            rx2%RvectorBlock(4),rtempVecFine%RvectorBlock(4),1.25_DP,0.0_DP)
!!        call lsyssc_vectorLinearComb  (
!!            rx3%RvectorBlock(4),rtempVecFine%RvectorBlock(4),0.75_DP,1.0_DP)
!!        call lsyssc_vectorLinearComb  (
!!            rx4%RvectorBlock(4),rtempVecFine%RvectorBlock(4),0.25_DP,1.0_DP)
!!
!!        call lsyssc_vectorLinearComb  (
!!            rx2%RvectorBlock(5),rtempVecFine%RvectorBlock(5),1.25_DP,0.0_DP)
!!        call lsyssc_vectorLinearComb  (
!!            rx3%RvectorBlock(5),rtempVecFine%RvectorBlock(5),0.75_DP,1.0_DP)
!!        call lsyssc_vectorLinearComb  (
!!            rx4%RvectorBlock(5),rtempVecFine%RvectorBlock(5),0.25_DP,1.0_DP)
!!            
!!        ! That is the second vector, save it.
!!        call setSpatialVector (rprojHier,rtempVecFine,rcoarseVector,2,&
!!            rtempVecCoarse,rtempVecFineScalar,ispacelevelcoarse,ispacelevelfine)
!!        
!!        
!!        
!!      else
!!      
!!        if (ispacelevelcoarse .ne. ispacelevelfine) then
!!        
!!          ! Restriction only in space.
!!          !
!!          ! Loop through the time steps
!!          do istep = 0,rcoarseVector%NEQtime-1
!!          
!!            ! Load the data
!!            call sptivec_getTimestepData (rfineVector,1+istep, rtempVecFine)
!!
!!            ! Process the restriction and save the data to the coarse vector.
!!            call mlprj_performRestriction (&
!!                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!!                rtempVecCoarse, rtempVecFine,rx1Scalar)
!!                
!!            call sptivec_setTimestepData (rcoarseVector,1+istep, rtempVecCoarse)
!!          
!!          end do
!!        
!!        end if
!!      
!!      end if
!!
!!      ! Release the temp vectors
!!      call lsyssc_releaseVector (rx1Scalar)
!!      call lsysbl_releaseVector (rx5)
!!      call lsysbl_releaseVector (rx4)
!!      call lsysbl_releaseVector (rx3)
!!      call lsysbl_releaseVector (rx2)
!!      call lsysbl_releaseVector (rx1)
!
!    case (3)
!      ! Quadratic restriction:
!      !
!      ! Prolongation 'distributes' information from the coarse grid nodes
!      ! to the fine grid nodes as follows:
!      !
!      !                v------- -0.125 -----------------+
!      !     +------------------ -0.125 ------v
!      !                v-- 0.75 --+-- 0.75 --v
!      !     +- 0.375 --v                     v-- 0.375 -+
!      !
!      !     X----------o----------X----------o----------X
!      ! Step1                   Step2                   Step3
!      !
!      ! Restriction is the 'adjoint' of the prolongation, so the corresponding
!      ! matrix is the transposed. Written nodewise this means, that restriction
!      ! has to 'collect' the distributed information with the same weights:
!      !
!      !                +------- -0.125 -----------------v
!      !     v------------------ -0.125 ------+
!      !                +-- 0.75 --v-- 0.75 --+
!      !     v- 0.375 --+                     +-- 0.375 -v
!      !
!      !     X----------o----------X----------o----------X
!      ! Step1                   Step2                   Step3
!      !
!      ! But because this is a restriction of a Finite-Difference RHS vector,
!      ! the RHS must be divided by h/(h/2)=2 !
!
!      ! We need some temporary vectors:
!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.false.)
!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx2,.false.)
!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.false.)
!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx4,.false.)
!      call lsysbl_createVecBlockIndirect (rtempVecFine,rx5,.false.)
!      
!      ! Create a scalar representation of rx1 for auxiliary reasons
!      call lsysbl_createScalarFromVec (rx1,rx1Scalar)
!     
!      ! DEBUG!!!
!      call lsysbl_getbase_double (rx1,p_Dx1)
!      call lsysbl_getbase_double (rx3,p_Dx3)
!      call lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
!      call lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
!         
!         
!      ! DEBUG:
!!      call lsysbl_clearVector (rtempVecFine,1.0_DP)
!!      call sptivec_setTimestepData (rfineVector, 1+0, rtempVecFine)
!!      call sptivec_setTimestepData (rfineVector, 1+1, rtempVecFine)
!!      call sptivec_setTimestepData (rfineVector, 1+2, rtempVecFine)
!!      call sptivec_setTimestepData (rfineVector, 1+3, rtempVecFine)
!!      call sptivec_setTimestepData (rfineVector, 1+4, rtempVecFine)
!!      call sptivec_setTimestepData (rfineVector, 1+5, rtempVecFine)
!!      call sptivec_setTimestepData (rfineVector, 1+6, rtempVecFine)
!!      call lsysbl_clearVector (rtempVecCoarse,-3.0_DP)
!!      call sptivec_setTimestepData (rcoarseVector, 1+5, rtempVecCoarse)
!!      call lsysbl_clearVector (rtempVecCoarse,1.0_DP)
!!      call sptivec_setTimestepData (rcoarseVector, 1+1, rtempVecCoarse)
!!      call sptivec_setTimestepData (rcoarseVector, 1+3, rtempVecCoarse)
!!      call lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
!      
!      if (rcoarseVector%NEQtime .ne. rfineVector%NEQtime) then       
!      
!        ! Load timestep 0
!        call sptivec_getTimestepData (rfineVector, 1+0, rx5)
!        
!        ! Prepare the temp vector
!        call lsysbl_clearVector (rtempVecFine)
!        
!        ! Loop through the timesteps in chunks of size 2
!        do istep = 0,rcoarseVector%NEQtime-3,2
!        
!          ! Get timestep i from the previous loop
!          call lsysbl_copyVector (rx5,rx1)
!        
!          ! Load timestep i+1..i+4 of the fine mesh
!          call sptivec_getTimestepData (rfineVector, 1+2*istep+1, rx2)
!          call sptivec_getTimestepData (rfineVector, 1+2*istep+2, rx3)
!          call sptivec_getTimestepData (rfineVector, 1+2*istep+3, rx4)
!          call sptivec_getTimestepData (rfineVector, 1+2*istep+4, rx5)
!          
!          ! Do the restriction in the temp vector. The vector is already
!          ! prepared from the last loop.
!          if (istep .eq. 0) then
!            ! In the very first step, put the (weighted) initial condition into
!            ! the first coarse timestep vector. In later timesteps, rtempVecFine contains
!            ! already the contribution from the last (left) interval and we only have
!            ! to update it with the remaining contributions (from the right).
!            call lsysbl_vectorLinearComb  (rx1,rtempVecFine,1.0_DP/2.0_DP,0.0_DP)
!          end if
!          call lsysbl_vectorLinearComb  (rx2,rtempVecFine,0.375_DP/2.0_DP,1.0_DP)
!          call lsysbl_vectorLinearComb  (rx4,rtempVecFine,-0.125_DP/2.0_DP,1.0_DP)
!          
!          ! Probably restrict the vector in space to the lower level.
!          ! Save the result.
!          if (ispacelevelcoarse .ne. ispacelevelfine) then
!            ! Space + time
!            call mlprj_performRestriction (&
!                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!                rtempVecCoarse,rtempVecFine,rx1Scalar)
!            call sptivec_setTimestepData (rcoarseVector, 1+istep, rtempVecCoarse)
!          else
!            ! Only time
!            call sptivec_setTimestepData (rcoarseVector, 1+istep, rtempVecFine)
!          end if
!          
!          ! Second subvector in the patch. Prepare the temp vector
!          call lsysbl_copyVector (rx3,rtempVecFine)
!          call lsysbl_vectorLinearComb  (rx2,rtempVecFine,0.75_DP/2.0_DP,1.0_DP/2.0_DP)
!          call lsysbl_vectorLinearComb  (rx4,rtempVecFine,0.75_DP/2.0_DP,1.0_DP)
!          if (ispacelevelcoarse .ne. ispacelevelfine) then
!            ! Space + time
!            call mlprj_performRestriction (&
!                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!                rtempVecCoarse,rtempVecFine,rx1Scalar)
!            call sptivec_setTimestepData (rcoarseVector, 1+istep+1, rtempVecCoarse)
!          else
!            ! Only time
!            call sptivec_setTimestepData (rcoarseVector, 1+istep+1, rtempVecFine)
!          end if
!          
!          ! The 3rd subvector is not calculated here, but in the beginning of
!          ! the next loop. We have to prepare the temp vector for that.
!          call lsysbl_copyVector (rx5,rtempVecFine)
!          call lsysbl_vectorLinearComb  (rx2,rtempVecFine,-0.125_DP/2.0_DP,1.0_DP/2.0_DP)
!          call lsysbl_vectorLinearComb  (rx4,rtempVecFine,0.375_DP/2.0_DP,1.0_DP)
!        
!        end do
!        
!        if (mod(rcoarseVector%NEQtime,2) .ne. 0) then
!          ! Equal number of timesteps / odd number of time-DOF's.
!          ! Save the last vector, that's it.
!          if (ispacelevelcoarse .ne. ispacelevelfine) then
!            ! Space + time
!            call mlprj_performRestriction (&
!                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!                rtempVecCoarse,rtempVecFine,rx1Scalar)
!            call sptivec_setTimestepData (rcoarseVector, rcoarseVector%NEQtime, rtempVecCoarse)
!          else
!            ! Only time
!            call sptivec_setTimestepData (rcoarseVector, rcoarseVector%NEQtime, rtempVecFine)
!          end if
!        else
!          ! Ok, the number of timesteps is even, that means there is only half a
!          ! patch at the end. We have to shift the solutions by half a patch
!          ! instead of a full patch and load the missing data.
!          call lsysbl_copyVector (rx3,rx1)
!          call lsysbl_copyVector (rx4,rx2)
!          call lsysbl_copyVector (rx5,rx3)
!          call sptivec_getTimestepData (rfineVector, rfineVector%NEQtime-1, rx4)
!          call sptivec_getTimestepData (rfineVector, rfineVector%NEQtime, rx5)
!          
!          ! The temp vector rtempVecFine now corresponds to a partially prepared
!          ! vector in the 'middle' of the patch like rx3. 
!          call lsysbl_vectorLinearComb  (rx4,rtempVecFine,0.75_DP/2.0_DP,1.0_DP)
!          if (ispacelevelcoarse .ne. ispacelevelfine) then
!            ! Space + time
!            call mlprj_performRestriction (&
!                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!                rtempVecCoarse,rtempVecFine,rx1Scalar)
!            call sptivec_setTimestepData (rcoarseVector, rcoarseVector%NEQtime-1, rtempVecCoarse)
!          else
!            ! Only time
!            call sptivec_setTimestepData (rcoarseVector, rcoarseVector%NEQtime-1, rtempVecFine)
!          end if
!          
!          ! Assemble the last vector in time.
!          call lsysbl_copyVector (rx5,rtempVecFine)
!          call lsysbl_vectorLinearComb (rx2,rtempVecFine,-0.125_DP/2.0_DP,1.0_DP/2.0_DP)
!          call lsysbl_vectorLinearComb (rx4,rtempVecFine,0.375_DP/2.0_DP,1.0_DP)
!        
!          ! Save the last vector, that's it.
!          if (ispacelevelcoarse .ne. ispacelevelfine) then
!            ! Space + time
!            call mlprj_performRestriction (&
!                rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!                rtempVecCoarse,rtempVecFine,rx1Scalar)
!            call sptivec_setTimestepData (rcoarseVector, rcoarseVector%NEQtime, rtempVecCoarse)
!          else
!            ! Only time
!            call sptivec_setTimestepData (rcoarseVector, rcoarseVector%NEQtime, rtempVecFine)
!          end if
!
!        end if
!        
!      end if
!
!      ! DEBUG!!!
!!      call lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
!!      do istep = 1,rcoarseVector%NEQtime
!!        call sptivec_getTimestepData (rcoarseVector, istep, rtempVecCoarse)
!!      end do
!      
!      ! Release the temp vectors
!      call lsyssc_releaseVector (rx1Scalar)
!      call lsysbl_releaseVector (rx5)
!      call lsysbl_releaseVector (rx4)
!      call lsysbl_releaseVector (rx3)
!      call lsysbl_releaseVector (rx2)
!      call lsysbl_releaseVector (rx1)
!
!    end select
!    
!    !call sptivec_saveToFileSequence (rcoarseVector,"(""coarseold.txt."",I5.5)",.true.,&
!    !    rtempVecFine)

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
    integer :: istep,ispacelevelfine,ispacelevelcoarse,itimelevelfine,itimelevelcoarse
    type(t_vectorScalar) :: rtempVecFineScalar,rtempVecFineScalar2
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
    select case (rprojHier%itimeOrder)
    case (0)
      call sptivec_createAccessPool (rfineVector,raccessPool,3)
    case (1)
      call sptivec_createAccessPool (rfineVector,raccessPool,3)
    case (2)
      call sptivec_createAccessPool (rfineVector,raccessPool,7)
    case (3)
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
            p_rx%RvectorBlock(3),rtempVecFine%RvectorBlock(3),dscaleDual*p_DaDual(icol),1.0_DP)
            
        call lsyssc_vectorLinearComb(&
            p_rx%RvectorBlock(4),rtempVecFine%RvectorBlock(4),dscaleDual*p_DaDual(icol),1.0_DP)

        call lsyssc_vectorLinearComb(&
            p_rx%RvectorBlock(5),rtempVecFine%RvectorBlock(5),dscaleDual*p_DaDual(icol),1.0_DP)
        
      end do
      
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

!    ! local variables
!    integer :: istep,ispacelevelcoarse,ispacelevelfine
!    type(t_vectorBlock) :: rx1,rx3
!    type(t_vectorScalar) :: rx1scalar
!    
!    ! DEBUG!!!
!    real(DP), dimension(:), pointer :: p_Dx1,p_Dx3
!    real(DP), dimension(:), pointer :: p_DtempVecCoarse,p_DtempVecFine
!
!    ! Get the space-time coarse and fine discretisations.
!    call sth_getLevel (rprojHier%p_rspaceTimeHierarchy,ilevelfine-1,&
!      ispaceLevel=ispacelevelcoarse)
!    call sth_getLevel (rprojHier%p_rspaceTimeHierarchy,ilevelfine,&
!      ispaceLevel=ispacelevelfine)
!
!    ! Prolongation 'distributes' information from the coarse grid nodes
!    ! to the fine grid nodes as follows:
!    !
!    ! Timestep:    n-1         n          n+1        (fine grid)
!    !               x <--1/2-- X --1/2--> x
!    !                         ^ \
!    !                        /   \
!    !                       +--1--+
!    !
!    ! Interpolation is now taking the mean of adjacent nodes.
!    ! Written nodewise this means:
!    !
!    ! Timestep:    n-1         n          n+1        (fine grid)
!    !               x --1/4--> X <--1/4-- x
!    !                         ^ \
!    !                        /   \
!    !                       +-1/2-+
!
!    ! We need two more temp vectors:
!    !
!    ! One for the current timestep
!    call lsysbl_createVecBlockIndirect (rtempVecFine,rx1,.false.)
!    
!    ! And one for the next timestep
!    call lsysbl_createVecBlockIndirect (rtempVecFine,rx3,.false.)
!    
!    ! Create a scalar representation of rx1 for auxiliary reasons
!    call lsysbl_createScalarFromVec (rx1,rx1Scalar)
!   
!    ! DEBUG!!!
!    call lsysbl_getbase_double (rx1,p_Dx1)
!    call lsysbl_getbase_double (rx3,p_Dx3)
!    call lsysbl_getbase_double (rtempVecCoarse,p_DtempVecCoarse)
!    call lsysbl_getbase_double (rtempVecFine,p_DtempVecFine)
!       
!    if (rcoarseVector%NEQtime .ne. rfineVector%NEQtime) then       
!      ! Load timestep 0,1 (fine grid) into the temp vectors.
!      call sptivec_getTimestepData (rfineVector, 1+0, rtempVecFine)
!      call sptivec_getTimestepData (rfineVector, 1+1, rx3)
!      
!      ! Perform restriction for the first time step.
!      !                          X <--1/4-- x
!      !                         ^ \
!      !                        /   \
!      !                       +-3/4-+
!      call lsysbl_vectorLinearComb (rx3,rtempVecFine,1._DP/4._DP,3._DP/4._DP)
!      
!      ! Probably restrict the vector in space to the lower level.
!      ! Save the result.
!      if (ispacelevelcoarse .ne. ispacelevelfine) then
!        ! Space + time
!        call mlprj_performInterpolation (&
!            rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!            rtempVecCoarse,rtempVecFine,rx1Scalar)
!        call sptivec_setTimestepData (rcoarseVector, 1+0, rtempVecCoarse)
!      else
!        ! Only time
!        call sptivec_setTimestepData (rcoarseVector, 1+0, rtempVecFine)
!      end if
!      
!      ! Loop through the time steps -- on the coarse grid!
!      do istep = 1,rcoarseVector%NEQtime-1-1
!      
!        ! rx3 was the 'right inner fine grid node x' in the above picture.
!        ! Copy it to rx1, so it gets the new 'left inner fine grid node x'
!        ! and load the new 'right inner fine grid node x' to rx1
!        ! as well as the 'inner fine grid node X' to rtempVecFine.
!        call lsysbl_copyVector (rx3,rx1)
!        
!        call sptivec_getTimestepData (rfineVector,1+2*istep, rtempVecFine)
!        call sptivec_getTimestepData (rfineVector,1+2*istep+1, rx3)
!        
!        ! In rtempVecFine, create the restriction of the fine grid vectors
!        ! rx1 and rx3 according to
!        !
!        ! Timestep:    n-1          n          n+1        (fine grid)
!        !              rx3 --1/4--> X <--1/4-- rx3
!        !             (old)        ^ \
!        !                         /   \
!        !                        +-1/2-+
!        
!        call lsysbl_vectorLinearComb (rx1,rtempVecFine,0.25_DP,0.5_DP)
!        call lsysbl_vectorLinearComb (rx3,rtempVecFine,0.25_DP,1.0_DP)
!        
!        ! Probably restrict the vector in space to the lower level.
!        ! Save the result.
!        if (ispacelevelcoarse .ne. ispacelevelfine) then
!          ! Space + time
!          call mlprj_performInterpolation (&
!              rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!              rtempVecCoarse,rtempVecFine,rx1Scalar)
!          call sptivec_setTimestepData (rcoarseVector, 1+istep, rtempVecCoarse)
!        else
!          ! Only time
!          call sptivec_setTimestepData (rcoarseVector, 1+istep, rtempVecFine)
!        end if
!      
!      end do
!
!      ! Last time step.
!
!      ! rx3 is now the new 'left inner fine grid node x'.
!      ! Load the 'inner fine grid node X' to rtempVecFine.
!      call sptivec_getTimestepData (rfineVector,1+2*(rcoarseVector%NEQtime-1), &
!          rtempVecFine)
!      
!      ! In rtempVecFine, create the restriction of the fine grid vectors
!      ! rx1 and rx3 according to
!      !
!      ! Timestep:    n-1          n       (fine grid)
!      !              rx3 --1/4--> X             
!      !                          ^ \
!      !                         /   \
!      !                        +-3/4-+
!      
!      call lsysbl_vectorLinearComb (rx3,rtempVecFine,1._DP/4._DP,3._DP/4._DP)
!      
!      ! Probably restrict the vector in space to the lower level.
!      ! Save the result.
!      if (ispacelevelcoarse .ne. ispacelevelfine) then
!        ! Space + time
!        call mlprj_performInterpolation (&
!            rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!            rtempVecCoarse, rtempVecFine,rx1Scalar)
!        call sptivec_setTimestepData (rcoarseVector, &
!            rcoarseVector%NEQtime, rtempVecCoarse)
!      else
!        ! Only time
!        call sptivec_setTimestepData (rcoarseVector, &
!            rcoarseVector%NEQtime, rtempVecFine)
!      end if
!
!    else
!    
!      if (ispacelevelcoarse .ne. ispacelevelfine) then
!      
!        ! Interpolation only in space.
!      
!        !
!        ! Loop through the time steps
!        do istep = 0,rcoarseVector%NEQtime-1
!        
!          ! Load the data
!          call sptivec_getTimestepData (rfineVector,1+istep, rtempVecFine)
!
!          ! Process the restriction and save the data to the coarse vector.
!          call mlprj_performInterpolation (&
!              rprojHier%p_rprojHierarchySpace%p_Rprojection(ispacelevelfine),&
!              rtempVecCoarse, rtempVecFine,rx1Scalar)
!              
!          call sptivec_setTimestepData (rcoarseVector,1+istep, rtempVecCoarse)
!        
!        end do
!
!      end if
!    
!    end if
!
!    ! Release the temp vectors
!    call lsyssc_releaseVector (rx1Scalar)
!    call lsysbl_releaseVector (rx3)
!    call lsysbl_releaseVector (rx1)

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
