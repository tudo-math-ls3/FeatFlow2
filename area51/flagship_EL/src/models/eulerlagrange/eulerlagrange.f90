!##############################################################################
!# ****************************************************************************
!# <name> euler_lagrange </name>
!# ****************************************************************************


module euler_lagrange

  use afcstabilisation
  use bilinearformevaluation
  use boundary
  use boundaryfilter
  use collection
  use derivatives
  use element
  use eulerlagrange_basic
  use eulerlagrange_callback
  use eulerlagrange_callback1d
  use eulerlagrange_callback2d
  use eulerlagrange_callback3d
  use flagship_basic
  use fparser
  use fsystem
  use genoutput
  use graph
  use groupfemsystem
  use hadaptaux
  use hadaptivity
  use linearformevaluation
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use pprocerror
  use pprocgradients
  use pprocindicator
  use pprocsolution
  use problem
  use solveraux
  use spatialdiscretisation
  use statistics
  use stdoperators
  use storage
  use thermodynamics
  use timestep
  use timestepaux
  use ucd

  implicit none

  private
  public :: eulerlagrange_init
  public :: eulerlagrange_timestep

  

contains

subroutine eulerlagrange_init(rparlist,p_rproblemLevel,rsolution,rtimestep)

    type(t_parlist), intent(inout) :: rparlist

    ! problem structure
    type(t_problem), intent(inout) :: rproblem

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! solver struchture
    type(t_solver), intent(inout), target :: rsolver

    ! primal solution vector
    type(t_vectorBlock), intent(inout), target :: rsolution

    ! collection structure
    type(t_collection), intent(inout) :: rcollection    

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel
    
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    
    real(dp) :: nvt, ndim, nel
    
    
    
   nvt  = p_rproblemLevel%rtriangulation%nvt
   ndim = p_rproblemLevel%rtriangulation%ndim
   nel  = p_rproblemLevel%rtriangulation%nel
    
    
    
    storage_getbase_double(p_rproblemLevel%rtriangulation%h_DvertexCoords,p_DvertexCoords)



end subroutine eulerlagrange_init

subroutine eulerlagrange_timestep(rparlist,p_rproblemLevel,rsolution,rtimestep)

    type(t_parlist), intent(inout) :: rparlist

    ! problem structure
    type(t_problem), intent(inout) :: rproblem

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! solver struchture
    type(t_solver), intent(inout), target :: rsolver

    ! primal solution vector
    type(t_vectorBlock), intent(inout), target :: rsolution

    ! collection structure
    type(t_collection), intent(inout) :: rcollection    

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

write *,* 'Blub'
stop

end subroutine eulerlagrange_timestep



end module euler_lagrange
