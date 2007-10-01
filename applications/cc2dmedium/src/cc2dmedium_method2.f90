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
    CALL parlst_readfromfile (rparamList, './data/nonlinsol_cc2d.dat')
    CALL parlst_readfromfile (rparamList, './data/output.dat')
    CALL parlst_readfromfile (rparamList, './data/paramtriang.dat')
    CALL parlst_readfromfile (rparamList, './data/bdconditions.dat')
    CALL parlst_readfromfile (rparamList, './data/timediscr.dat')
    CALL parlst_readfromfile (rparamList, './data/postprocessing.dat')
  
  END SUBROUTINE

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
    
    ! A structure for the postprocessing.
    TYPE(t_c2d2postprocessing) :: rpostprocessing
    
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
    CALL cc2dmedium2_getDAT (p_rproblem%rparamList)
    
    ! Ok, parameters are read in.
    ! Get the output levels during the initialisation phase and during the program.
    CALL c2d2_initOutput (p_rproblem)
    
    ! Print the configuration to the terminal
    IF (p_rproblem%MSHOW_Initialisation .GE. 2) THEN
      CALL output_line ('Parameters:')
      CALL parlst_info (p_rproblem%rparamList)
    END IF
    
    ! Evaluate these parameters and initialise global data in the problem
    ! structure for global access.
    CALL c2d2_initParameters (p_rproblem)
    
    ! So now the different steps - one after the other.
    !
    ! Initialisation
    !
    ! Parametrisation & Triangulation
    IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
      CALL output_separator (OU_SEP_MINUS)
      CALL output_line('Initialising parametrisation / triangulation...')
    END IF
    
    CALL c2d2_initParamTriang (p_rproblem)
    
    ! Print mesh information
    IF (p_rproblem%MSHOW_Initialisation .GE. 2) THEN
      CALL output_lbrk ()
      CALL output_line ('Mesh statistics:')
      CALL output_lbrk ()
      DO i=p_rproblem%NLMIN,p_rproblem%NLMAX
        CALL tria_infoStatistics (p_rproblem%RlevelInfo(i)%rtriangulation,&
            i .EQ. p_rproblem%NLMIN,i)
      END DO
    END IF
    
    ! Discretisation
    IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
      CALL output_separator (OU_SEP_MINUS)
      CALL output_line('Initialising discretisation...')
    END IF
    CALL c2d2_initDiscretisation (p_rproblem)    

    IF (p_rproblem%MSHOW_Initialisation .GE. 2) THEN
      CALL output_lbrk ()
      CALL output_line ('Discretisation statistics:')
      DO i=p_rproblem%NLMIN,p_rproblem%NLMAX
        CALL output_lbrk ()
        CALL output_line ('Level '//sys_siL(i,5))
        CALL dof_infoDiscrBlock (p_rproblem%RlevelInfo(i)%p_rdiscretisation,.FALSE.)
      END DO
    END IF
    
    ! And all the other stuff...
    IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
      CALL output_separator (OU_SEP_MINUS)
      CALL output_line('Initialising postprocessing...')
    END IF
    CALL c2d2_initPostprocessing (p_rproblem,rpostprocessing)
    
    IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
      CALL output_separator (OU_SEP_MINUS)
      CALL output_line('Initialising matrices/vectors...')
    END IF
    CALL c2d2_allocMatVec (p_rproblem,rvector,rrhs)    

    IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
      CALL output_separator (OU_SEP_MINUS)
      CALL output_line('Initialising analytic boundary conditions...')
    END IF
    CALL c2d2_initAnalyticBC (p_rproblem)   

    ! On all levels, generate the static matrices used as templates
    ! for the system matrix (Laplace, B, Mass,...)
    IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
      CALL output_separator (OU_SEP_MINUS)
      CALL output_line('Generating basic matrices...')
    END IF
    CALL c2d2_generateBasicMatrices (p_rproblem)

    ! Create the solution vector -- zero or read from file.
    IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
      CALL output_separator (OU_SEP_MINUS)
      CALL output_line('Initialising initial solution vector...')
    END IF
    CALL c2d2_initInitialSolution (p_rproblem,rvector)

    ! Now choose the algorithm. Stationary or time-dependent simulation?
    IF (p_rproblem%itimedependence .EQ. 0) THEN
    
      ! Stationary simulation
      !
      ! Generate the RHS vector.
      IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
        CALL output_separator (OU_SEP_MINUS)
        CALL output_line('Generating RHS vector...')
      END IF
      CALL c2d2_generateBasicRHS (p_rproblem,rrhs)
      
      ! Generate discrete boundary conditions
      IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
        CALL output_separator (OU_SEP_MINUS)
        CALL output_line('Generating discrete boundary conditions...')
      END IF
      CALL c2d2_initDiscreteBC (p_rproblem,rvector,rrhs)

      ! Implementation of boundary conditions
      CALL c2d2_implementBC (p_rproblem,rvector,rrhs,.TRUE.,.TRUE.)
    
      ! Solve the problem
      IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
        CALL output_separator (OU_SEP_MINUS)
        CALL output_line('Invoking stationary solver...')
        CALL output_separator (OU_SEP_MINUS)
      END IF
      CALL c2d2_solve (p_rproblem,rvector,rrhs)
    
      ! Postprocessing
      CALL c2d2_postprocessingStationary (p_rproblem,rvector,rpostprocessing)
      
    ELSE
    
      ! Time dependent simulation with explicit time stepping.
      !
      ! Generate the RHS vector for the first time step.
      IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
        CALL output_separator (OU_SEP_MINUS)
        CALL output_line('Generating RHS vector...')
      END IF
      CALL c2d2_generateBasicRHS (p_rproblem,rrhs)
      
      ! Initialise the boundary conditions, but 
      ! don't implement any boundary conditions as the nonstationary solver
      ! doesn't like this.
      IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
        CALL output_separator (OU_SEP_MINUS)
        CALL output_line('Generating discrete boundary conditionsof first time step...')
      END IF
      CALL c2d2_initDiscreteBC (p_rproblem,rvector,rrhs)
      
      ! Call the nonstationary solver to solve the problem.
      IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
        CALL output_separator (OU_SEP_MINUS)
        CALL output_line('Invoking nonstationary solver...')
        CALL output_separator (OU_SEP_MINUS)
      END IF
      CALL c2d2_solveNonstationary (p_rproblem,rvector,rrhs,rpostprocessing)
      
    END IF
    
    ! (Probably) write final solution vector
    CALL c2d2_writeSolution (p_rproblem,rvector)
    
    ! Cleanup
    CALL c2d2_doneMatVec (p_rproblem,rvector,rrhs)
    CALL c2d2_doneBC (p_rproblem)
    CALL c2d2_doneDiscretisation (p_rproblem)
    CALL c2d2_donepostprocessing (rpostprocessing)    
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
