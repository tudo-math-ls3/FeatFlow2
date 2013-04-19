!##############################################################################
!# ****************************************************************************
!# <name> stokesdbg2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# </purpose>
!##############################################################################

module stokesdbg2d

use fsystem
use storage
use genoutput
use collection, only: t_collection

use stokesdbg_core
use stokesdbg_driver2d

implicit none

private

public :: stokesdbg_run2d
  
contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine stokesdbg_run2d(rparam, idriver)
  type(t_parlist), intent(inout) :: rparam
  integer, intent(in) :: idriver

!</subroutine>

  type(t_problem) :: rproblem
  type(t_system) :: rsystem
  type(t_collection) :: rcollect
  type(t_linearForm) :: rform
  real(DP), dimension(2), target :: DerrorUL2, DerrorUH1
  real(DP), dimension(1), target :: DerrorPL2
  type(t_errorScVec) :: rerrorU, rerrorP
  integer :: ilvl, imesh, iinfo

    ! allocate levels
    call stdbg_initProblem(rproblem, rparam)
    
    ! store driver id
    rproblem%idriver = idriver

    ! store problem parameters
    rcollect%DquickAccess(1) = rproblem%dnu
    rcollect%DquickAccess(2) = rproblem%dalpha
    rcollect%DquickAccess(3) = rproblem%dbeta
    rcollect%DquickAccess(4) = rproblem%dgamma
    rcollect%IquickAccess(1) = idriver
    rcollect%IquickAccess(2) = 0
    rcollect%IquickAccess(3) = rproblem%ilevelMin
    rcollect%IquickAccess(4) = rproblem%ilevelMax


    ! fetch additional info
    call parlst_getvalue_int(rparam, '', 'VORTEX_VELO', rcollect%IquickAccess(6), 0)
    call parlst_getvalue_int(rparam, '', 'VORTEX_PRES', rcollect%IquickAccess(7), 0)

    ! read in triangulations
    call parlst_getvalue_int(rparam, '', 'MESH', imesh, 0)
    select case(imesh)
    case (0)
      call stdbg_initTriangulation2D(rproblem, rparam, 'QUAD')
    case (1)
      call stdbg_initTriangulation2D(rproblem, rparam, 'QUADX4R')
    case (2)
      call stdbg_initTriangulation2D(rproblem, rparam, 'QUADIRR2')
    end select

    ! initialise discretisations
    call stdbg_initDiscretisation(rproblem, rparam)
    
    ! initialise multilevel projections
    call stdbg_initProjections(rproblem)
    
    ! initialise matrices
    call stdbg_initMatrices(rproblem)
    
    ! initialise BCs
    call stdrv_initBoundaryConditions(rproblem, rparam, idriver)
    
    ! loop over all levels
    do ilvl = rproblem%ilevelMin, rproblem%ilevelMax
      
      call output_line("Processing Level " // trim(sys_sil(ilvl, 4)) // "...")

      ! store current level in collection
      rcollect%IquickAccess(2) = ilvl

      ! initialise system
      call stdbg_initSystem(rproblem, rsystem, ilvl)

      ! Assemble the RHS vector
      rform%itermCount = 1
      rform%Idescriptors(1) = DER_FUNC
      rcollect%IquickAccess(5) = 1
      call linf_buildVectorScalar (rform,.true., rsystem%rvecRhs%RvectorBlock(1), &
          rproblem%Rlevels(ilvl)%rcubInfo, stdrv_funcRhs, rcollect)
      rcollect%IquickAccess(5) = 2
      call linf_buildVectorScalar (rform,.true., rsystem%rvecRhs%RvectorBlock(2), &
          rproblem%Rlevels(ilvl)%rcubInfo, stdrv_funcRhs, rcollect)

      ! Filter solution and rhs vectors
      call vecfil_discreteBCsol(rsystem%rvecSol)
      call vecfil_discreteBCrhs(rsystem%rvecRhs)

      ! initialise filter chain
      call stdrv_initFilterChain(rproblem, rsystem, rparam)

      ! initialise multigrid solver
      call stdbg_initMultigrid(rproblem, rsystem, rparam)
    
      ! solve system
      call stdbg_solve(rproblem, rsystem)
      
      ! compute errors
      ! Set up the error structure for velocity
      rerrorU%p_RvecCoeff => rsystem%rvecSol%RvectorBlock(1:2)
      rerrorU%p_DerrorL2 => DerrorUL2
      rerrorU%p_DerrorH1 => DerrorUH1
    
      ! Set up the error structure for pressure
      rerrorP%p_RvecCoeff => rsystem%rvecSol%RvectorBlock(3:3)
      rerrorP%p_DerrorL2 => DerrorPL2
    
      ! Calculate errors of velocity and pressure against analytic solutions.
      call pperr_scalarVec(rerrorU, stdrv_funcVelocity2D, rcollect);
      call pperr_scalarVec(rerrorP, stdrv_funcPressure2D, rcollect);
    
      ! Calculate errors of velocity field
      rproblem%p_Dstat(1,rsystem%ilevel) = sqrt(DerrorUL2(1)**2 + DerrorUL2(2)**2)
      rproblem%p_Dstat(2,rsystem%ilevel) = sqrt(DerrorUH1(1)**2 + DerrorUH1(2)**2)
      
      ! Store error of pressure
      rproblem%p_Dstat(4,rsystem%ilevel) = derrorPL2(1)
    
      ! Calculate divergence of velocity field
      call stdbg_aux_calcDiv2D(rproblem%p_Dstat(3,rsystem%ilevel), &
        rsystem%rvecSol%RvectorBlock(1:2), rproblem%Rlevels(ilvl)%rcubInfo)

      ! write VTK if desired
      call stdbg_writeVTK(rproblem, rsystem, rparam)

      ! release system
      call stdbg_doneSystem(rsystem)
    
    end do ! next level
    
    ! print statistics
    call stdbg_printMeshStatistics2D(rproblem)
    call stdbg_printMatrixStatistics(rproblem)
    call stdbg_printErrorStatistics(rproblem)
  
    ! release problem
    call stdbg_doneProblem(rproblem)

  end subroutine

end module
