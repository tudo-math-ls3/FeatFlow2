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

MODULE heatcond_method5

  USE fsystem
  USE storage
  USE linearsolver
  USE boundary
  USE bilinearformevaluation
  USE linearformevaluation
  USE cubature
  USE matrixfilters
  USE vectorfilters
  USE bcassembly
  USE triangulation
  USE spatialdiscretisation
  USE sortstrategy
  USE coarsegridcorrection
  USE ucd
  USE timestepping
  USE genoutput
  
  USE collection
  USE paramlist
    
  USE heatcond_callback
  
  USE heatcond_basic
  USE heatcond_matvec
  USE heatcond_boundarycondition
  USE heatcond_partridiscr
  USE heatcond_solver
  USE heatcond_timeloop
  
  IMPLICIT NONE
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE heatcond5
  
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
    TYPE(t_parlist) :: rparams

    ! A problem structure for our problem
    TYPE(t_problem), TARGET :: rproblem
    
    ! An initial RHS vector and a solution vector
    TYPE(t_vectorBlock) :: rrhs,rvector
    
    INTEGER :: i
    
    ! Initialise the parameter list
    CALL parlst_init(rparams)
    
    ! Initialise the collection.
    CALL collct_init (rproblem%rcollection)
    CALL collct_setvalue_parlst (rproblem%rcollection, 'PARAMS', rparams, .TRUE.)
    
    ! Read in the parameters from the DAT file and initialise the basic
    ! structures with these.
    CALL hc5_initparameters (rparams,rproblem)

    ! Add space for level information in the collection
    DO i=1,rproblem%ilvmax
      CALL collct_addlevel_all (rproblem%rcollection)
    END DO
    
    ! Allocate memory for the level information
    ALLOCATE (rproblem%RlevelInfo(rproblem%ilvmax))
    
    ! So now the different steps - one after the other.
    !
    ! Initialisation
    CALL hc5_initParamTriang (rproblem%ilvmin,rproblem%ilvmax,rproblem)
    CALL hc5_initDiscretisation (rproblem)    
    CALL hc5_initMatVec (rproblem,rparams)    

    ! Use the auxiliary RHS vector on the finest level to create an
    ! initial RHS and solution vector, which we pass later to the timeloop.
    CALL lsysbl_createVecBlockIndirect (rproblem%rrhs,rrhs,.FALSE.)
    CALL lsysbl_createVecBlockIndirect (rproblem%rrhs,rvector,.TRUE.)

    ! Calculate the initial RHS, don't incorporate any BC's.
    CALL hc5_calcRHS (rproblem,rrhs)  
    
    ! Discretise the boundary conditions
    CALL hc5_initDiscreteBC (rproblem)
    
    ! Implement them into the initial solution vector, as we have a zero
    ! vector as initial solution.
    rvector%p_rdiscreteBC => rproblem%rrhs%p_rdiscreteBC
    CALL vecfil_discreteBCsol (rvector)
    
    ! Initialise the solver
    CALL hc5_initSolver (rproblem)
    
    ! Call the timeloop to solve the problem
    CALL hc5_timeloop (rproblem,rvector,rrhs)
    
    ! Release the solver, we donÄt need it anymore
    CALL hc5_doneSolver (rproblem)
    
    ! Cleanup
    CALL hc5_doneMatVec (rproblem)
    CALL hc5_doneBC (rproblem)
    CALL hc5_doneDiscretisation (rproblem)
    CALL hc5_doneParamTriang (rproblem)
    
    ! Release memory for level information
    DEALLOCATE (rproblem%RlevelInfo)
    
    ! Release parameter list
    CALL collct_deletevalue (rproblem%rcollection,'PARAMS')
    CALL parlst_done (rparams)
    
    ! Release RHS and solution vector
    CALL lsysbl_releaseVector (rvector)
    CALL lsysbl_releaseVector (rrhs)

    ! Print some statistical data about the collection - anything forgotten?
    CALL output_lbrk ()
    CALL output_line ('Remaining collection statistics:')
    CALL output_line ('--------------------------------')
    CALL output_lbrk ()
    CALL collct_printStatistics (rproblem%rcollection)
    
    ! Finally release the collection.
    CALL collct_done (rproblem%rcollection)
    
  END SUBROUTINE

END MODULE
