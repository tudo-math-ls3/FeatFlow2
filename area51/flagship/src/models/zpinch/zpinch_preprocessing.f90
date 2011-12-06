!##############################################################################
!# ****************************************************************************
!# <Name> zpinch_preprocessing </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all preprocessing routines which are required
!# to solve the time-dependent magnetohydrodynamic equations in the
!# one-, two- or three-dimensional domain $\Omega$.
!#
!# The following routines are available:
!#
!# 1.) zpinch_initSolvers
!#     -> Initialises the solve structures from the parameter list.
!#
!# 2.) zpinch_initProblemDescriptor
!#     -> Initialises the abstract problem descriptor based on the
!#        parameter settings given by the parameter list.
!#
!# </purpose>
!##############################################################################

module zpinch_preprocessing

#include "hydro.h"

  use collection
  use fparser
  use genoutput
  use hydro_preprocessing
  use paramlist
  use problem
  use solveraux
  use timestep
  use timestepaux
  use transport_preprocessing

  implicit none

  private

  public :: zpinch_initSolvers
  public :: zpinch_initProblemDescriptor

contains

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_initSolvers(rparlist, ssectionName, rtimestep, rsolver)

!<description>
    ! This subroutine initialises the time-stepping structure and
    ! the top-level solver structure from the parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section names in parameter list
    character(LEN=*), intent(in) :: ssectionName
!</input>

!<output>
    ! time-stepping structure
    type(t_timestep), intent(out) :: rtimestep

    ! solver struchture
    type(t_solver), intent(out) :: rsolver
!</output>
!</subroutine>

    ! section name for the top-level solver
    character(LEN=SYS_STRLEN) :: ssolverName

    ! section name for time-stepping scheme
    character(LEN=SYS_STRLEN) :: stimestepName


    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'timestep', stimestepName)

    ! Initialise time-stepping
    call tstep_createTimestep(rparlist, stimestepName, rtimestep, 2)


    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'solver', ssolverName)

    ! Initialise solver structures
    call solver_createSolver(rparlist, ssolverName, rsolver)
    call solver_adjustHierarchy(rsolver)
    call solver_updateStructure(rsolver)

  end subroutine zpinch_initSolvers

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_initProblemDescriptor(rparlist, ssectionName,&
      nlmin, nlmax, rproblemDescriptor)

!<description>
    ! This subroutine initialises the abstract problem descriptor
    ! using the parameters settings defined in the parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! minimum/maximum problem level
    integer, intent(in) :: nlmin, nlmax
!</input>

!<output>
    ! problem descriptor
    type(t_problemDescriptor), intent(out) :: rproblemDescriptor
!</output>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: ssectionNameHydro
    character(LEN=SYS_STRLEN) :: ssectionNameTransport

    ! abstract problem descriptor
    type(t_problemDescriptor) :: rproblemDescriptorHydro
    type(t_problemDescriptor) :: rproblemDescriptorTransport

    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'subapplication', ssectionNameHydro, isubstring=1)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'subapplication', ssectionNameTransport, isubstring=2)

    call hydro_initProblemDescriptor(rparlist, ssectionNameHydro,&
        nlmin, nlmax, rproblemDescriptorHydro)
    call transp_initProblemDescriptor(rparlist, ssectionNameTransport,&
        nlmin, nlmax, rproblemDescriptorTransport)
    
    rproblemDescriptor = problem_combineDescriptors(rproblemDescriptorHydro,&
                                                    rproblemDescriptorTransport)

  end subroutine zpinch_initProblemDescriptor
end module zpinch_preprocessing
