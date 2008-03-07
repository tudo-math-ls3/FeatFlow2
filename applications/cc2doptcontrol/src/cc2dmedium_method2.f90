!##############################################################################
!# ****************************************************************************
!# <name> cc2dmedium_method2 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module solves an optimal control problem for the stationary and
!# nonstationary Navier-Stokes optimal control problem 
!#
!#  $$ min J(y,u) = 1/2||y-z||_{L^2} + \gamma/2||y(T)-z(T)||_{L^2} + \alphga/2||u||^2 $$
!#
!#  $$- \nu Delta(y) + y*\Nabla(y) + \Nabla p = f $$
!#  $$ \Nabla \cdot y = 0$$
!#  $$- \nu Delta(\lambda) - y*\Nabla(\lambda) + \lambda\Nabla y + \Nabla \xi = y-z $$
!#  $$ \Nabla \cdot \lambda = 0$$
!#              
!#
!# on a 2D domain for a 2D function $y=(y_1,y_2)$, a pressure $p$,
!# a dual velocity $\lambda$ and a dual pressure $\xi$. $u$ is the control
!# and $z$ a desired flow field.
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
  
  USE externalstorage
  
  IMPLICIT NONE
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc2dmedium2_getLogFiles (slogfile,serrorfile)
  
!<description>
  ! Temporarily reads the output DAT file to get the names of the output
  ! files.
!</description>

!<output>
  ! Name of the message log file.
  CHARACTER(LEN=*), INTENT(OUT) :: slogfile
  
  ! Name of the error log file.
  CHARACTER(LEN=*), INTENT(OUT) :: serrorfile
!</output>

!</subroutine>

    TYPE(t_parlist) :: rparlist
    CHARACTER(LEN=SYS_STRLEN) :: sstring

    ! Init parameter list that accepts parameters for output files
    CALL parlst_init (rparlist)

    ! Read parameters that configure the output
    CALL parlst_readfromfile (rparlist, './data/output.dat')
    
    ! Now the real initialisation of the output including log file stuff!
    CALL parlst_getvalue_string (rparlist,'GENERALOUTPUT',&
                                'smsgLog',sstring,'')
    READ(sstring,*) slogfile

    CALL parlst_getvalue_string (rparlist,'GENERALOUTPUT',&
                                'serrorLog',sstring,'')
    READ(sstring,*) serrorfile
    
    ! That temporary parameter list is not needed anymore.
    CALL parlst_done (rparlist)
    
    END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc2dmedium2_getDAT (rparamList)
  
!<description>
  ! Reads in all DAT files into the parameter list rparlist
!</description>

!<inputoutput>
  ! The parameter list where the values of the DAT files should be stored.
  ! The structure must have been initialised, the parameters are just added
  ! to the list.
  TYPE(t_parlist), INTENT(INOUT) :: rparamList
!</inputoutput>

!</subroutine>

    ! Each 'readfromfile' command adds the parameter of the specified file 
    ! to the parameter list.
    CALL parlst_readfromfile (rparamList, './data/discretisation.dat')
    CALL parlst_readfromfile (rparamList, './data/linsol_cc2d.dat')
    ! CALL parlst_readfromfile (rparamList, './data/nonlinsol_cc2d.dat')
    CALL parlst_readfromfile (rparamList, './data/output.dat')
    CALL parlst_readfromfile (rparamList, './data/paramtriang.dat')
    CALL parlst_readfromfile (rparamList, './data/bdconditions.dat')
    CALL parlst_readfromfile (rparamList, './data/timediscr.dat')
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc2dmedium2optc
  
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
    ! Initialise the external storage management.
    
    CALL exstor_init (999,100)
    !CALL exstor_attachDirectory('./ff2storage')
    
    ! Allocate memory for the problem; it's rather large.
    ALLOCATE (p_rproblem)
    
    ! Initialise the collection
    CALL collct_init (p_rproblem%rcollection)

    ! Initialise the parameter list object. This creates an empty parameter list.
    CALL parlst_init (p_rproblem%rparamList)
    
    ! Read parameters from the INI/DAT files into the parameter list. 
    CALL cc2dmedium2_getDAT (p_rproblem%rparamList)
    
    ! Ok, parameters are read in.
    ! Get the output levels during the initialisation phase and during the program.
    CALL cc_initOutput (p_rproblem)
    OU_LINE_LENGTH = 132
    
    ! Evaluate these parameters and initialise global data in the problem
    ! structure for global access.
    CALL cc_initParameters (p_rproblem)
    
    DO i=1,p_rproblem%NLMAX
      CALL collct_addlevel_all (p_rproblem%rcollection)
    END DO
    
    ! So now the different steps - one after the other.
    !
    ! Initialisation
    CALL cc_initParamTriang (p_rproblem)
    CALL cc_initDiscretisation (p_rproblem)    
    CALL cc_allocMatVec (p_rproblem,rvector,rrhs)   
    
    ! Print information about the discretisation
    CALL output_line ('Discretisation statistics:')
    CALL output_line ('--------------------------')
    DO i=p_rproblem%NLMIN,p_rproblem%NLMAX
      CALL output_lbrk ()
      CALL output_line ('Level '//sys_siL(i,10))
      CALL output_line ('---------')
      CALL dof_infoDiscrBlock (p_rproblem%RlevelInfo(i)%rdiscretisation,.FALSE.)
    END DO
     
    ! On all levels, generate the static matrices used as templates
    ! for the system matrix (Laplace, B, Mass,...)
    CALL cc_generateBasicMatrices (p_rproblem)

    ! Create the solution vector -- zero or read from file.
    CALL cc_initInitialSolution (p_rproblem,rvector)
    
    ! Now choose the algorithm. Stationary or time-dependent simulation?
    IF (p_rproblem%itimedependence .EQ. 0) THEN
    
      ! Stationary simulation

      ! Read the (stationary) target flow.
      CALL cc_initTargetFlow (p_rproblem)

      ! Generate the RHS vector.
      CALL cc_generateBasicRHS (p_rproblem,rrhs)
      
      ! Generate discrete boundary conditions
      CALL cc_initDiscreteBC (p_rproblem,rvector,rrhs)

      ! Implementation of boundary conditions
      CALL cc_implementBC (p_rproblem,rvector=rvector,rrhs=rrhs)
    
      ! Solve the problem
      CALL cc_solve (p_rproblem,rvector,rrhs)
    
      ! Postprocessing
      CALL cc_postprocessingStationary (p_rproblem,rvector)
      
      ! Release the target flow
      CALL cc_doneTargetFlow (p_rproblem)
      
    ELSE
    
      ! Time dependent simulation with explicit time stepping.
      
      ! Initialise the boundary conditions for the 0th time step, but 
      ! don't implement any boundary conditions as the nonstationary solver
      ! doesn't like this.
      CALL cc_initDiscreteBC (p_rproblem,rvector,rrhs)
      
      ! Don't read the target flow, this is done in 
      ! cc_solveNonstationaryDirect!

      ! Call the nonstationary solver to solve the problem.
      CALL cc_solveNonstationaryDirect (p_rproblem)
      
    END IF
    
    ! (Probably) write final solution vector
    CALL cc_writeSolution (p_rproblem,rvector)
    
    ! Cleanup
    CALL cc_doneMatVec (p_rproblem,rvector,rrhs)
    CALL cc_doneBC (p_rproblem)
    CALL cc_doneDiscretisation (p_rproblem)
    CALL cc_doneParamTriang (p_rproblem)
    
    ! Release parameters from the DAT/INI files from the problem structure.
    CALL cc_doneParameters (p_rproblem)

    ! Release the parameter list
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
    
    ! Information about external storage usage
    CALL output_lbrk ()
    CALL exstor_info ()
    
    ! Clean up the external storage management
    CALL exstor_done ()
    
  END SUBROUTINE

END MODULE
