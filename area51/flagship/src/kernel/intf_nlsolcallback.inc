! This file describes the interfaces of the callback routines of the nonlinear
! defect-correction loop solver.
! ***************************************************************************  
  
  interface

  ! *************************************************************************** 

!<subroutine>

  subroutine fcb_setBoundary(rproblemLevel, rtimestep, rsolver,&
                             ru, rres, ru0, istatus)

  use linearsystemblock
  use problem
  use solver

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(IN)         :: rtimestep

    ! solver structure
    type(t_solver), intent(IN)           :: rsolver

    ! initial solution vector
    type(t_vectorBlock), intent(IN)      :: ru0
!</input>

!<inputoutput>
    ! multigrid level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! solution block vector
    type(t_vectorBlock), intent(INOUT)   :: ru

    ! nonlinear residual block vector
    type(t_vectorBlock), intent(INOUT)   :: rres

    ! status of the callback function
    integer, intent(INOUT), optional     :: istatus
!</inputoutput>
!</subroutine

  end subroutine fcb_setBoundary

  ! *************************************************************************** 

!<subroutine>

  subroutine fcb_calcResidual(rproblemLevel, rtimestep, rsolver,&
                              ru, ru0, rrhs, rres, ite, imode, istatus)

  use linearsystemblock
  use problem
  use solver

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(IN) :: rtimestep

    ! current solution block vector
    type(t_vectorBlock), intent(IN) :: ru

    ! initial solution block vector
    type(t_vectorBlock), intent(IN) :: ru0

    ! number of nonlinear iteration
    integer, intent(IN) :: ite

    ! mode: (1) primal or (2) dual problem
    integer, intent(IN) :: imode
!</input>

!<inputoutput>
    ! multigrid level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel
    
    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver

    ! constant right-hand side block vector
    type(t_vectorBlock), intent(INOUT) :: rrhs

    ! nonlinear residual block vector
    type(t_vectorBlock), intent(INOUT) :: rres

    ! status of the callback function
    integer, intent(INOUT), optional :: istatus
!</inputoutput>
!</subroutine>

  end subroutine fcb_calcResidual

  ! *************************************************************************** 

!<subroutine>

  subroutine fcb_calcRHS(rproblemLevel, rtimestep, rsolver,&
                         ru, ru0, rrhs, istep, imode, istatus)

  use linearsystemblock
  use problem
  use solver

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(IN) :: rtimestep

    ! current solution block vector
    type(t_vectorBlock), intent(IN) :: ru

    ! initial solution block vector
    type(t_vectorBlock), intent(IN) :: ru0

    ! number of explicit step
    integer, intent(IN) :: istep

    ! mode: (1) primal or (2) dual problem
    integer, intent(IN) :: imode
!</input>

!<inputoutput>
    ! multigrid level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel
    
    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver

    ! constant right-hand side block vector
    type(t_vectorBlock), intent(INOUT) :: rrhs

    ! status of the callback function
    integer, intent(INOUT), optional :: istatus
!</inputoutput>
!</subroutine>

  end subroutine fcb_calcRHS

  ! ***************************************************************************  

!<subroutine>

  subroutine fcb_calcJacobian(rproblemLevel, rtimestep, rsolver, ru, ru0,&
                              imode, bfailure, istatus)

  use linearsystemblock
  use problem
  use solver

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(IN) :: rtimestep

    ! current solution block vector
    type(t_vectorBlock), intent(IN) :: ru

    ! initial solution block vector
    type(t_vectorBlock), intent(IN) :: ru0

    ! mode: (1) primal or (2) dual problem
    integer, intent(IN) :: imode

    ! OPTIONAL: Newton subiteration failed, 
    !           return to defect correction
    logical, intent(IN), optional :: bfailure
!</input>

!<inputoutput>
    ! multigrid level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel
    
    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver

    ! status of the callback function
    integer, intent(INOUT), optional :: istatus
!</inputoutput>
!</subroutine>

  end subroutine fcb_calcJacobian

  ! ***************************************************************************
  
  subroutine fcb_applyJacobian(rproblemLevel, rx, ry, cx, cy, istatus)

  use linearsystemblock
  use problem

!<input>
    ! multigrid level structure
    type(t_problemLevel), intent(IN)  :: rproblemLevel

    ! vector to which the Jacobian should be applied to
    type(t_vectorBlock), intent(IN)    :: rx

    ! factor by which vector rx should be scaled
    real(DP), intent(IN)               :: cx

    ! factor by which vector ry should be scaled
    real(DP), intent(IN)               :: cy
!</input>

!<inputoutput>
    ! vector to which the result should be added
    type(t_vectorBlock), intent(INOUT) :: ry

    ! status of the callback function
    integer, intent(INOUT), optional   :: istatus
!</inputoutput>
!</subroutine>

  end subroutine fcb_applyJacobian

  end interface
