!##############################################################################
!# ****************************************************************************
!# <name> cc2dmini_method2 </name>
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
!# different subroutines, which are separated from the main application;
!# they can be found in the files cc2dminim2basic.f90, cc2dminim2boundary,...
!# The communication between these subroutines is done using an 
!# application-specific structure saving problem data
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

MODULE cc2dmini_method2

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
    
  USE cc2dminim2basic
  USE cc2dminim2init
  USE cc2dminim2boundary
  USE cc2dminim2discretisation
  USE cc2dminim2postprocessing
  USE cc2dminim2stationary
  
  IMPLICIT NONE
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc2dmini2
  
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
    
    INTEGER :: i
    
    ! Ok, let's start. 
    !
    ! Allocate the problem structure on the heap -- it's rather large.
    ALLOCATE(p_rproblem)

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
    
    ! Ok, parameters are read in.
    ! Evaluate these parameters and initialise global data in the problem
    ! structure for global access.
    CALL c2d2_initParameters (p_rproblem)
    
    ! So now the different steps - one after the other.
    !
    ! Initialisation
    CALL c2d2_initParamTriang (p_rproblem)
    CALL c2d2_initDiscretisation (p_rproblem)    
    CALL c2d2_initMatVec (p_rproblem)    
    CALL c2d2_initDiscreteBC (p_rproblem)
    
    ! Implementation of boundary conditions
    CALL c2d2_implementBC (p_rproblem)
    
    ! Solve the problem
    CALL c2d2_solve (p_rproblem)
    
    ! Postprocessing
    CALL c2d2_postprocessing (p_rproblem)
    
    ! Cleanup
    CALL c2d2_doneMatVec (p_rproblem)
    CALL c2d2_doneBC (p_rproblem)
    CALL c2d2_doneDiscretisation (p_rproblem)
    CALL c2d2_doneParamTriang (p_rproblem)
    
    ! Release parameters from the DAT/INI files from the problem structure.
    CALL c2d2_doneParameters (p_rproblem)

    ! Release the parameter list
    CALL collct_deleteValue (p_rproblem%rcollection,'INI')
    CALL parlst_done (p_rproblem%rparamList)
    
    ! Print some statistical data about the collection - anything forgotten?
    PRINT *
    PRINT *,'Remaining collection statistics:'
    PRINT *,'--------------------------------'
    PRINT *
    CALL collct_printStatistics (p_rproblem%rcollection)
    
    ! Finally release the collection and the problem structure.
    CALL collct_done (p_rproblem%rcollection)
    
    DEALLOCATE(p_rproblem)
    
  END SUBROUTINE

END MODULE
