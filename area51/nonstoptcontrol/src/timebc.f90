!##############################################################################
!# ****************************************************************************
!# <name> timebc </name>
!# ****************************************************************************
!#
!# <purpose>
!# Realises the initial condition.
!# </purpose>
!##############################################################################

module timebc

  use fsystem
  use storage
  use genoutput
  use linearalgebra
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use collection
  use vectorio
  use derivatives
  use timediscretisation
  
  use collection
  
  use analyticprojection
  
  use spacetimehierarchy
  use fespacehierarchy
  use spacetimematrices

  use physics
  
  use callback
  
  implicit none

contains

  ! ***************************************************************************

  subroutine spop_applyInitCondSolRhs (rmatrix,rsolution, rrhs)

  ! Incorporates the initial condition into the solution and right hand side.
  
  ! Space-time matrix of the discrete problem.
  type(t_spaceTimeMatrix), intent(in) :: rmatrix
  
  ! Solution vector.
  type(t_spaceTimeVector), intent(inout) :: rsolution
  
  ! OPTIONAL: RHS vector.
  type(t_spaceTimeVector), intent(inout), optional :: rrhs
  
    ! local varibales
    type(t_vectorBlock) :: rvector,rvector2
    type(t_collection) :: rcollection
    type(t_matrixBlock) :: rsubmatrix
    integer :: ispaceLevel

    ! Create a temp vector, clear it.
    call lsysbl_createVectorBlock(rsolution%p_rspaceDiscr,rvector,.true.)
    call sptivec_getTimestepData (rsolution, 1, rvector)
    
    ! Initialise the collection with information about our problem.
    rcollection%DquickAccess (1) = rmatrix%p_rphysics%dtimemin ! dtime
    
    rcollection%IquickAccess (2) = rmatrix%p_rphysics%cequation
    rcollection%IquickAccess (3) = rmatrix%p_rphysics%creferenceProblem
    rcollection%DquickAccess (2) = rmatrix%p_rphysics%doptControlAlpha
    rcollection%DquickAccess (3) = rmatrix%p_rphysics%dpar
    rcollection%DquickAccess (4) = rmatrix%p_rphysics%dtimemin
    rcollection%DquickAccess (5) = rmatrix%p_rphysics%dtimemax
    
    ! Project the analytical solution to this vector.
    select case (rmatrix%p_rphysics%cequation)
    case (0,2)
      ! 1D/2D Heat equation
      call lsyssc_clearVector (rvector%RvectorBlock(1))
      
      rcollection%IquickAccess (1) = 1 ! icomponent
      call anprj_discrDirect (rvector%RvectorBlock(1),&
          ferrFunction, rcollection)
      
    case (1)
      ! 2D Stokes equation
      call lsyssc_clearVector (rvector%RvectorBlock(1))
      call lsyssc_clearVector (rvector%RvectorBlock(2))
      
      rcollection%IquickAccess (1) = 1 ! icomponent
      call anprj_discrDirect (rvector%RvectorBlock(1),&
          ferrFunction, rcollection)

      rcollection%IquickAccess (1) = 2
      call anprj_discrDirect (rvector%RvectorBlock(2),&
          ferrFunction, rcollection)
      
    end select
    
    ! Save it.
    call sptivec_setTimestepData (rsolution, 1, rvector)
    
    ! If a RHS is present, incorporate the initial condition into the RHS.
    if (present(rrhs)) then
      call lsysbl_createVectorBlock(rsolution%p_rspaceDiscr,rvector2,.true.)
      
      ! Multiply the solution with the matrix in the 1st timestep to
      ! get the correct RHS.
      call sth_getLevel(rmatrix%p_rspaceTimeHierarchy,rmatrix%ilevel,&
          ispaceLevel=ispaceLevel)
      call stmat_allocSubmatrix (rmatrix%cmatrixType,rmatrix%p_rphysics,&
        rmatrix%p_rmatVecTempl(ispaceLevel),rsubmatrix)
      call stmat_getSubmatrix (rmatrix, ispaceLevel, 1, 1, rsubmatrix)
      
      call lsysbl_blockMatVec (rsubmatrix,rvector,rvector2,1.0_DP,0.0_DP)

      call lsysbl_releaseMatrix (rsubmatrix)
      
      ! Incorporate the primal RHS into the given RHS.
      call sptivec_getTimestepData (rrhs, 1, rvector)
      
      select case (rmatrix%p_rphysics%cequation)
      case (0,2)
        ! 1D/2D Heat equation
        call lsyssc_copyVector (rvector2%RvectorBlock(1),rvector%RvectorBlock(1))
        
      case (1)
        ! 2D Stokes equation
        call lsyssc_copyVector (rvector2%RvectorBlock(1),rvector%RvectorBlock(1))
        call lsyssc_copyVector (rvector2%RvectorBlock(2),rvector%RvectorBlock(2))
        
      end select
      
      call sptivec_setTimestepData (rrhs, 1, rvector)
      
      call lsysbl_releaseVector (rvector2)
    end if
    
    call lsysbl_releaseVector (rvector)

  end subroutine

  ! ***************************************************************************

  subroutine spop_applyInitCondDef (rphysics,rdefect)

  ! Incorporates the initial condition into the defect.

  ! Underlying physics
  type(t_physics), intent(in) :: rphysics
  
  ! Defect vector.
  type(t_spaceTimeVector), intent(inout) :: rdefect

    ! local varibales
    type(t_vectorBlock) :: rvector

    ! The defect for the initial condition is just zero in the primal equation.
    call lsysbl_createVectorBlock(rdefect%p_rspaceDiscr,rvector,.false.)
    call sptivec_getTimestepData (rdefect, 1, rvector)
    
    select case (rphysics%cequation)
    case (0,2)
      ! 1D/2D Heat equation
      call lsyssc_clearVector (rvector%RvectorBlock(1))
    case (1)
      ! 2D Stokes equation
      call lsyssc_clearVector (rvector%RvectorBlock(1))
      call lsyssc_clearVector (rvector%RvectorBlock(2))
    end select
    
    call sptivec_setTimestepData (rdefect, 1, rvector)
    
    call lsysbl_releaseVector (rvector)

  end subroutine

end module
