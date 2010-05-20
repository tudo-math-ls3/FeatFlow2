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

    ! Order of projection in time.
    ! =1: first order implicit Euler, 
    ! =2: Full 1-step Theta scheme up to second order (Crank Nicolson)
    integer :: itimeOrder = 0
    
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

  end type

!</typeblock>

!</types>


contains

  ! ***************************************************************************

!<subroutine>

  subroutine sptipr_initProjection (rprojHier,rspaceTimeHierarchy,&
      rprojHierarchySpace,rphysics,iorderTimeProlRest)
  
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

    ! Remember the physics; necessary so we know how and what to project
    rprojHier%p_rphysics => rphysics

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
      if (rprojHier%itimeOrder .ne. 4) then
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
      else
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
      
    case (4)
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
      
    case (4)
      call sptipr_getProlMatrixPrimal (rspaceTimeHierarchy,ilevel,4,rprolMatrix)

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
    select case (rprojHier%itimeOrder)
    case (0)
      call sptivec_createAccessPool (rfineVector%p_rspaceDiscr,raccessPool,3)
    case (1)
      call sptivec_createAccessPool (rfineVector%p_rspaceDiscr,raccessPool,3)
    case (2)
      call sptivec_createAccessPool (rfineVector%p_rspaceDiscr,raccessPool,3)
    case (3)
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
    select case (rprojHier%itimeOrder)
    case (0)
      call sptivec_createAccessPool (rfineVector,raccessPool,3)
    case (1)
      call sptivec_createAccessPool (rfineVector,raccessPool,3)
    case (2)
      call sptivec_createAccessPool (rfineVector,raccessPool,7)
    case (3)
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
    
    ! Release the buffer.
    call sptivec_releaseAccessPool(raccessPool)
    
    call lsyssc_releaseVector (rtempVecFineScalar)

    ! DEBUG!!!
    call sptivec_saveToFileSequence (rcoarseVector,"(""./ns/coarse.txt."",I5.5)",.true.,&
        rtempVecCoarse)
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
    select case (rprojHier%itimeOrder)
    case (0)
      call sptivec_createAccessPool (rfineVector,raccessPool,3)
    case (1)
      call sptivec_createAccessPool (rfineVector,raccessPool,3)
    case (2)
      call sptivec_createAccessPool (rfineVector,raccessPool,7)
    case (3)
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
