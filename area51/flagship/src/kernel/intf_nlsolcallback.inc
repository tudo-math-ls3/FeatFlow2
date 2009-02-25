! This file describes the interfaces of the callback routines of the nonlinear
! defect-correction loop solver.
! ***************************************************************************  
  
  interface

  ! *************************************************************************** 

!<subroutine>

  subroutine fcb_setBoundary(rproblemLevel, rtimestep, rsolver,&
                             ru, rres, ru0, rcollection)

  use collection
  use linearsystemblock
  use problem
  use solveraux
  use timestepaux

!<input>
    ! time-stepping structure
    type(t_timestep), intent(IN) :: rtimestep

    ! solver structure
    type(t_solver), intent(IN) :: rsolver

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: ru0
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! solution vector
    type(t_vectorBlock), intent(INOUT) :: ru

    ! nonlinear residual vector
    type(t_vectorBlock), intent(INOUT) :: rres

    ! collection structure
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine

  end subroutine fcb_setBoundary

  ! *************************************************************************** 

!<subroutine>

  subroutine fcb_calcResidual(rproblemLevel, rtimestep, rsolver,&
                              ru, ru0, rrhs, rres, ite, rcollection)

  use collection
  use linearsystemblock
  use problem
  use solveraux
  use timestepaux

!<input>
    ! time-stepping structure
    type(t_timestep), intent(IN) :: rtimestep

    ! solution vector
    type(t_vectorBlock), intent(IN) :: ru

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: ru0

    ! number of nonlinear iteration
    integer, intent(IN) :: ite
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel
    
    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver

    ! right-hand side vector
    type(t_vectorBlock), intent(INOUT) :: rrhs

    ! residual vector
    type(t_vectorBlock), intent(INOUT) :: rres

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

  end subroutine fcb_calcResidual

  ! *************************************************************************** 

!<subroutine>

  subroutine fcb_calcRHS(rproblemLevel, rtimestep, rsolver,&
                         ru, ru0, rrhs, istep, rcollection)

  use collection
  use linearsystemblock
  use problem
  use solveraux
  use timestepaux

!<input>
    ! time-stepping structure
    type(t_timestep), intent(IN) :: rtimestep

    ! solution vector
    type(t_vectorBlock), intent(IN) :: ru

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: ru0

    ! number of explicit step
    integer, intent(IN) :: istep
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel
    
    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver

    ! right-hand side vector
    type(t_vectorBlock), intent(INOUT) :: rrhs

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

  end subroutine fcb_calcRHS

  ! ***************************************************************************  

!<subroutine>

  subroutine fcb_calcJacobian(rproblemLevel, rtimestep, rsolver,&
                              ru, ru0, bfailure, rcollection)
  use collection
  use linearsystemblock
  use problem
  use solveraux
  use timestepaux

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(IN) :: rtimestep

    ! current solution block vector
    type(t_vectorBlock), intent(IN) :: ru

    ! initial solution block vector
    type(t_vectorBlock), intent(IN) :: ru0

    ! Newton subiteration failed, return to defect correction
    logical, intent(IN) :: bfailure
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel
    
    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

  end subroutine fcb_calcJacobian

  ! ***************************************************************************
  
  subroutine fcb_applyJacobian(rproblemLevel, rx, ry, cx, cy, rcollection)

  use collection
  use fsystem
  use linearsystemblock
  use problem

!<input>
    ! problem level structure
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

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

  end subroutine fcb_applyJacobian

  end interface
