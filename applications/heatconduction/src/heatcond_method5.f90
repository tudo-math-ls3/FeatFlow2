!##############################################################################
!# ****************************************************************************
!# <name> heatcond_method5 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstation program how to solve a simple 
!# heat conduction problem with constant coefficients 
!# on a simple domain.
!#
!# The routine splits up the tasks of reading the domain, creating 
!# triangulations, discretisation, solving, postprocessing and creanup into
!# different subroutines. The communication between these subroutines
!# is done using an application-specific structure saving problem data
!# as well as a collection structure for the communication with callback
!# routines.
!#
!# The different tasks of setting up the triangulation, discretisation, solver,
!# etc. are separated into different sub-modules as follows:
!#
!# - heatcond_basic
!#   -> Defines a basic problem structure that contains all parameters for
!#      the simulation
!#
!# - heatcond_callback
!#   -> Contains user-defined routines for setting up the RHS or for
!#      configuring values on the boundary.
!#
!# - heatcond_partridiscr
!#   -> Contains routines to read in the parametrisation, triangulation and
!#      to set up the discretisation (element, cubature rules,...=
!#
!# - heatcond_matvec
!#   -> Contains routines to initialise matrices and vectors.
!#
!# - heatcond_bounadrycondition
!#   -> Contains routines to set up analytical and discrete boudary conditions.
!#
!# - heatcond_solver
!#   -> Configures the linear solver for the problem.
!#
!# - heatcond_timeloop
!#   -> Contains the actual time stepping algorithm and nonstationary solver.
!# </purpose>
!##############################################################################

module heatcond_method5

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use sortstrategy
  use coarsegridcorrection
  use ucd
  use timestepping
  use genoutput
  
  use collection
  use paramlist
    
  use heatcond_callback
  
  use heatcond_basic
  use heatcond_matvec
  use heatcond_boundarycondition
  use heatcond_partridiscr
  use heatcond_solver
  use heatcond_timeloop
  
  implicit none
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine heatcond5
  
!<description>
  ! This is a 'separated' heatcond solver for solving a nonstationary heat 
  ! conduction problem. The different tasks of the problem are separated into
  ! subroutines. The problem uses a problem-specific structure for the 
  ! communication: All subroutines add their generated information to the
  ! structure, so that the other subroutines can work with them.
  ! (This is somehow a cleaner implementation than using a collection!).
  ! For the communication to callback routines of black-box subroutines
  ! (matrix-assembly), a collection is used.
  !
  ! The following tasks are performed by the subroutines:
  !
  ! 1.) Read in parametrisation
  ! 2.) Read in triangulation
  ! 3.) Set up RHS
  ! 4.) Set up matrix
  ! 5.) Create solver structure
  ! 6.) Solve the problem
  ! 7.) Write solution to GMV file
  ! 8.) Release all variables, finish
!</description>

!</subroutine>

    ! A paramlist structure with parameters from the dat file
    type(t_parlist) :: rparams

    ! A problem structure for our problem
    type(t_problem), target :: rproblem
    
    ! An initial RHS vector and a solution vector
    type(t_vectorBlock) :: rrhs,rvector
    
    integer :: i
    
    ! Initialise the parameter list
    call parlst_init(rparams)
    
    ! Initialise the collection.
    call collct_init (rproblem%rcollection)
    call collct_setvalue_parlst (rproblem%rcollection, 'PARAMS', rparams, .true.)
    
    ! Read in the parameters from the DAT file and initialise the basic
    ! structures with these.
    call hc5_initparameters (rparams,rproblem)

    ! Add space for level information in the collection
    do i=1,rproblem%ilvmax
      call collct_addlevel_all (rproblem%rcollection)
    end do
    
    ! Allocate memory for the level information
    allocate (rproblem%RlevelInfo(rproblem%ilvmax))
    
    ! So now the different steps - one after the other.
    !
    ! Initialisation
    call hc5_initParamTriang (rproblem%ilvmin,rproblem%ilvmax,rproblem)
    call hc5_initDiscretisation (rproblem)    
    call hc5_initMatVec (rproblem,rparams)    

    ! Use the auxiliary RHS vector on the finest level to create an
    ! initial RHS and solution vector, which we pass later to the timeloop.
    call lsysbl_createVecBlockIndirect (rproblem%rrhs,rrhs,.false.)
    call lsysbl_createVecBlockIndirect (rproblem%rrhs,rvector,.true.)

    ! Calculate the initial RHS, don't incorporate any BC's.
    call hc5_calcRHS (rproblem,rrhs)  
    
    ! Discretise the boundary conditions
    call hc5_initDiscreteBC (rproblem)
    
    ! Implement them into the initial solution vector, as we have a zero
    ! vector as initial solution.
    rvector%p_rdiscreteBC => rproblem%rrhs%p_rdiscreteBC
    call vecfil_discreteBCsol (rvector)
    
    ! Initialise the solver
    call hc5_initSolver (rproblem)
    
    ! Call the timeloop to solve the problem
    call hc5_timeloop (rproblem,rvector,rrhs)
    
    ! Release the solver, we donÄt need it anymore
    call hc5_doneSolver (rproblem)
    
    ! Cleanup
    call hc5_doneMatVec (rproblem)
    call hc5_doneBC (rproblem)
    call hc5_doneDiscretisation (rproblem)
    call hc5_doneParamTriang (rproblem)
    
    ! Release memory for level information
    deallocate (rproblem%RlevelInfo)
    
    ! Release parameter list
    call collct_deletevalue (rproblem%rcollection,'PARAMS')
    call parlst_done (rparams)
    
    ! Release RHS and solution vector
    call lsysbl_releaseVector (rvector)
    call lsysbl_releaseVector (rrhs)

    ! Print some statistical data about the collection - anything forgotten?
    call output_lbrk ()
    call output_line ('Remaining collection statistics:')
    call output_line ('--------------------------------')
    call output_lbrk ()
    call collct_printStatistics (rproblem%rcollection)
    
    ! Finally release the collection.
    call collct_done (rproblem%rcollection)
    
  end subroutine

end module
