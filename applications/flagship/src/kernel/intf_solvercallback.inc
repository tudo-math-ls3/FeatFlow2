! -*- mode: f90; -*-

!##############################################################################
!# ****************************************************************************
!# <name> intf_solvercallback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This file describes the interfaces of the callback routines of the
!# nonlinear solver.
!# </purpose>
!##############################################################################

  interface

    ! ***************************************************************************

!<subroutine>

    subroutine fcb_nlsolverCallback(rproblemLevel, rtimestep, rsolver,&
        rsolution, rsolution0, rrhs, rres, istep,&
        ioperationSpec, rcollection, istatus, rsource)

      use collection
      use fsystem
      use linearsystemblock
      use problem
      use solverbase
      use timestepbase

!<input>
      ! old solution vector from previous time step
      type(t_vectorBlock), intent(IN) :: rsolution0

      ! number of solver step
      integer, intent(IN) :: istep

      ! specifier for operations
      integer(I32), intent(IN) :: ioperationSpec

      ! OPTIONAL: given source vector
      type(t_vectorBlock), intent(IN), optional :: rsource
!</input>

!<inputoutput>
      ! problem level structure
      type(t_problemLevel), intent(INOUT) :: rproblemLevel

      ! time-stepping structure
      type(t_timestep), intent(INOUT) :: rtimestep

      ! solver structure
      type(t_solver), intent(INOUT) :: rsolver

      ! solution vector
      type(t_vectorBlock), intent(INOUT) :: rsolution

      ! right-hand side vector
      type(t_vectorBlock), intent(INOUT) :: rrhs

      ! residual vector
      type(t_vectorBlock), intent(INOUT) :: rres

      ! collection structure
      type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>

!<output>
      ! status flag
      integer, intent(OUT) :: istatus
!</output>
!</subroutine>

    end subroutine fcb_nlsolverCallback

  end interface
