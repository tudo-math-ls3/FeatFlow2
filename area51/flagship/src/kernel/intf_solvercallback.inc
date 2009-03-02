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
                                    rsolution, rsolutionInitial,&
                                    rvector, istep, ioperationSpec,&
                                    rcollection, istatus)

      use collection
      use linearsystemblock
      use problem
      use solveraux
      use timestepaux

!<input>
      ! initial solution vector
      type(t_vectorBlock), intent(IN) :: rsolutionInitial

      ! number of solver step
      integer, intent(IN) :: istep

      ! specifier for operations
      integer(I32), intent(IN) :: ioperationSpec
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

      ! vector to calculate
      type(t_vectorBlock), intent(INOUT) :: rvector

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
