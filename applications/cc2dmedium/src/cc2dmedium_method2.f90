!##############################################################################
!# ****************************************************************************
!# <name> cc2dmedium_method2 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstation program how to solve a stationary
!# Navier-Stokes problem 
!#
!#              $$- \nu Laplace(u) + u*grad(u) + \Nabla p = f $$
!#              $$ \Nable \cdot p = 0$$
!#
!# on a 2D domain for a 2D function $u=(u_1,u_2)$ and a pressure $p$.
!#
!# The routine splits up the tasks of reading the domain, creating 
!# triangulations, discretisation, solving, postprocessing and creanup into
!# different subroutines. The communication between these subroutines
!# is done using an application-specific structure saving problem data
!# as well as a collection structure for the communication with callback
!# routines.
!#
!# For the nonlinearity, the nonlinear solver is invoked. The
!# defect that is setted up there is preconditioned by a linear Multigrid
!# solver with a simple-VANCA smoother/preconditioner for
!# 2D saddle point problems, Jacobi-Type. As coarse grid solver,
!# UMFPACK is used.
!# </purpose>
!##############################################################################

MODULE cc2dmedium_method2

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
  USE coarsegridcorrection
  USE spdiscprojection
  USE nonlinearsolver
  USE paramlist
  
  USE collection
  USE convection
    
  USE cc2dmediumm2basic
  USE cc2dmediumm2init
  USE cc2dmediumm2boundary
  USE cc2dmediumm2discretisation
  USE cc2dmediumm2postprocessing
  USE cc2dmediumm2stationary
  
  IMPLICIT NONE
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc2dmedium2
  
!<description>
  ! This is a 'separated' Navier-Stokes solver for solving a Navier-Stokes
  ! problem. The different tasks of the problem are separated into
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

    ! A problem structure for our problem
    TYPE(t_problem), POINTER :: p_rproblem
    
    ! A structure for the solution vector and the RHS vector of the problem.
    TYPE(t_vectorBlock) :: rvector,rrhs
    
    INTEGER :: i
    
    ! Ok, let's start. 
    !
    ! Allocate memory fo rthe problem structure -- it's rather large!
    ALLOCATE (p_rproblem)
    
    ! Initialise the collection
    CALL collct_init (p_rproblem%rcollection)
    DO i=1,NNLEV
      CALL collct_addlevel_all (p_rproblem%rcollection)
    END DO
    
    ! Initialise the parameter list object. This creates an empty parameter list.
    CALL parlst_init (p_rproblem%rparamList)
    
    ! Add the parameter list to the collection so that the parameters
    ! from the DAT/INI files are available everywhere where we have the 
    ! collection.
    CALL collct_setvalue_parlst(p_rproblem%rcollection,'INI',&
                                p_rproblem%rparamList,.TRUE.)

    ! Read parameters from the INI/DAT files into the parameter list. 
    ! Each 'readfromfile' command adds the parameter of the specified file 
    ! to the parameter list.
    CALL parlst_readfromfile (p_rproblem%rparamList, './data/discretisation.dat')
    CALL parlst_readfromfile (p_rproblem%rparamList, './data/linsol_cc2d.dat')
    CALL parlst_readfromfile (p_rproblem%rparamList, './data/nonlinsol_cc2d.dat')
    CALL parlst_readfromfile (p_rproblem%rparamList, './data/output.dat')
    CALL parlst_readfromfile (p_rproblem%rparamList, './data/paramtriang.dat')
    CALL parlst_readfromfile (p_rproblem%rparamList, './data/bdconditions.dat')
    CALL parlst_readfromfile (p_rproblem%rparamList, './data/timediscr.dat')
    
    ! Ok, parameters are read in.
    ! Get the output levels during the initialisation phase and during the program.
    CALL c2d2_initOutput (p_rproblem)
    
    ! Evaluate these parameters and initialise global data in the problem
    ! structure for global access.
    CALL c2d2_initParameters (p_rproblem)
    
    ! So now the different steps - one after the other.
    !
    ! Initialisation
    CALL c2d2_initParamTriang (p_rproblem)
    CALL c2d2_initDiscretisation (p_rproblem)    
    CALL c2d2_allocMatVec (p_rproblem,rvector,rrhs)    
    CALL c2d2_initAnalyticBC (p_rproblem)   

    ! On all levels, generate the static matrices and the basic
    ! system matrix.
    CALL c2d2_generateBasicMatrices (p_rproblem)

    ! Now choose the algorithm. Stationary or time-dependent simulation?
    IF (p_rproblem%itimedependence .EQ. 0) THEN
    
      ! Stationary simulation
      !
      ! Generate the RHS vector.
      CALL c2d2_generateBasicRHS (p_rproblem,rrhs)
      
      ! Generate discrete boundary conditions
      CALL c2d2_initDiscreteBC (p_rproblem,rvector,rrhs)

      ! Implementation of boundary conditions
      CALL c2d2_implementBC (p_rproblem,rvector,rrhs,.TRUE.,.TRUE.,.TRUE.)
    
      ! Solve the problem
      CALL c2d2_solve (p_rproblem,rvector,rrhs)
    
      ! Postprocessing
      CALL c2d2_postprocessingStationary (p_rproblem,rvector)
      
    ELSE
    
      ! Time dependent simulation with explicit time stepping.
      !
      ! Generate the RHS vector for the first time step.
      CALL c2d2_generateBasicRHS (p_rproblem,rrhs)
      
      ! Initialise the boundary conditions, but 
      ! don't implement any boundary conditions as the nonstationary solver
      ! doesn't like this.
      CALL c2d2_initDiscreteBC (p_rproblem,rvector,rrhs)
      
      ! Call the nonstationary solver to solve the problem.
      CALL c2d2_solveNonstationary (p_rproblem,rvector,rrhs)
      
    END IF
    
    ! Cleanup
    CALL c2d2_doneMatVec (p_rproblem,rvector,rrhs)
    CALL c2d2_doneBC (p_rproblem)
    CALL c2d2_doneDiscretisation (p_rproblem)
    CALL c2d2_doneParamTriang (p_rproblem)
    
    ! Release parameters from the DAT/INI files from the problem structure.
    CALL c2d2_doneParameters (p_rproblem)

    ! Release the parameter list
    CALL collct_deleteValue (p_rproblem%rcollection,'INI')
    CALL parlst_done (p_rproblem%rparamList)
    
    ! Print some statistical data about the collection - anything forgotten?
    CALL output_lbrk ()
    CALL output_line ('Remaining collection statistics:')
    CALL output_line ('--------------------------------')
    CALL output_lbrk ()
    CALL collct_printStatistics (p_rproblem%rcollection)
    
    ! Finally release the collection and the problem structure.
    CALL collct_done (p_rproblem%rcollection)
    
    DEALLOCATE(p_rproblem)
    
  END SUBROUTINE

END MODULE
