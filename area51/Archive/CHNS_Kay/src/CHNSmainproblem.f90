!##############################################################################
!# ****************************************************************************
!# <name> ccmainproblem </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a program how to solve a nonstationary
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
!# routines. For the nonlinearity, the nonlinear solver is invoked.
!# </purpose>
!##############################################################################

module CHNSmainproblem

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
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use statistics
  use dofmapping
  
  use collection
  use convection
    
  use ccbasic
  use ccinitgeneralparameters
  use ccinitparamtriang
  use ccgeneraldiscretisation
  use ccpostprocessing
  use ccnonstationary
  use ccboundarycondition
!~~~~~~~~~~~~~~~~~use CahnHilliard module~~~~~~~~~~~~~~~~~~~~~
  use CahnHilliard
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  implicit none
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine CHNS2dmain
  
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
    type(t_problem), pointer :: p_rproblem
    
    ! A structure for the solution vector and the RHS vector of the problem.
    type(t_vectorBlock) :: rvector,rrhs
    
    ! A structure for the postprocessing.
    type(t_c2d2postprocessing) :: rpostprocessing

!~~~~~~~~~~~~~~~~~~~~~~~CH probelm~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! A paramlist structure with parameters from the dat file
    type(t_parlist) :: rCHparams

    ! A problem structure saving problem-dependent information.
    type(t_CHproblem) :: rCHproblem
!
!    ! A structure for the solution vector and the RHS vector of the problem.
    type(t_vectorBlock) :: rCHvector,rCHrhs
!
!    ! A structure for the postprocessing.
    type(t_c2d2postprocessing) :: rCHpostprocessing
!~~~~~~~~~~~~~~~~~~~~~~~End of CH problem parameters~~~~~~~~~~~~~~~~~
!
    
    ! Timer objects for stopping time
    type(t_timer) :: rtimerTotal
    type(t_timer) :: rtimerGridGeneration
    type(t_timer) :: rtimerMatrixGeneration
    type(t_timer) :: rtimerSolver
    
    integer :: i

!  added~~~~for initialize the solution of CH problem
    real(DP), dimension(:,:), pointer ::  p_DvertexCoords
    real(DP), dimension(:), pointer ::  p_vectordata
    real(DP), dimension(:), pointer ::  p_data
!    integer :: mcai
    
    ! Ok, let's start.

    ! Initialise the timers by zero:
    call stat_clearTimer(rtimerTotal)
    call stat_clearTimer(rtimerGridGeneration)
    call stat_clearTimer(rtimerMatrixGeneration)
    call stat_clearTimer(rtimerSolver)

    ! Start the timer
    call stat_startTimer(rtimerTotal)
    
    ! Allocate memory fo rthe problem structure -- it's rather large!
    allocate (p_rproblem)
    
	! Initialise the collection
    call collct_init (p_rproblem%rcollection)
    
    ! Initialise the parameter list object. This creates an empty parameter list.
    call parlst_init (p_rproblem%rparamList)
    
    ! Read parameters from the INI/DAT files into the parameter list.
    call cc2d_getDAT (p_rproblem%rparamList)
    
    ! Ok, parameters are read in.
    ! Get the output levels during the initialisation phase and during the program.
    call cc_initOutput (p_rproblem)

    ! Print the configuration to the terminal
    if (p_rproblem%MSHOW_Initialisation .ge. 2) then
      call output_line ('Parameters:')
      call output_lbrk ()
      call parlst_info (p_rproblem%rparamList)
    end if
    
    ! Evaluate these parameters and initialise global data in the problem
    ! structure for global access.
    call cc_initParameters (p_rproblem)
    
    ! So now the different steps - one after the other.
    !
    ! Initialisation
    !
    ! Parametrisation & Triangulation
    if (p_rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line('Initialising parametrisation / triangulation...')
    end if
    
    call stat_startTimer(rtimerGridGeneration)
    
    call cc_initParamTriang (p_rproblem)
    
    call stat_stopTimer(rtimerGridGeneration)
    call output_lbrk ()
    call output_line ("Time for mesh generation: "//&
      trim(sys_sdL(rtimerGridGeneration%delapsedReal,10)))
      
    p_rproblem%rstatistics%dtimeGridGeneration = &
      p_rproblem%rstatistics%dtimeGridGeneration + rtimerGridGeneration%delapsedReal
    
    ! Print mesh information
    if (p_rproblem%MSHOW_Initialisation .ge. 2) then
      call output_lbrk ()
      call output_line ('Mesh statistics:')
      call output_lbrk ()
      do i=p_rproblem%NLMIN,p_rproblem%NLMAX
        call tria_infoStatistics (p_rproblem%RlevelInfo(i)%rtriangulation,&
            i .eq. p_rproblem%NLMIN,i)
      end do
    end if

    ! Discretisation
    if (p_rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line('Initialising discretisation...')
    end if
    call cc_initDiscretisation (p_rproblem)

    if (p_rproblem%MSHOW_Initialisation .ge. 2) then
      call output_lbrk ()
      call output_line ('Discretisation statistics:')
      do i=p_rproblem%NLMIN,p_rproblem%NLMAX
        call output_lbrk ()
        call output_line ('Level '//sys_siL(i,5))
        call dof_infoDiscrBlock (p_rproblem%RlevelInfo(i)%rdiscretisation,.false.)
      end do
    end if

    ! And all the other stuff...
    if (p_rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line('Initialising postprocessing...')
    end if
    call cc_initPostprocessing (p_rproblem,rpostprocessing)
    
    if (p_rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line('Initialising matrices/vectors...')
    end if
    
    call stat_startTimer(rtimerMatrixGeneration)
    
    call cc_allocMatVec (p_rproblem,rvector,rrhs)
    
    call stat_stopTimer(rtimerMatrixGeneration)
    call output_lbrk ()
    call output_line ("Time for matrix initialisation: "//&
      trim(sys_sdL(rtimerGridGeneration%delapsedReal,10)))

    ! On all levels, generate the static matrices used as templates
    ! for the system matrix (Laplace, B, Mass,...)
    if (p_rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line('Generating basic matrices...')
    end if
    call cc_generateBasicMatrices (p_rproblem)

    ! Create the solution vector -- zero or read from file.
    if (p_rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line('Initialising initial solution vector...')
    end if
    call cc_initInitialSolution (p_rproblem,rvector)

    write(*,*)'End of NS problem initialization, we have no problem in NS'


!~~~~~CH prolem, read in data and initilzation~~~~~~~~~~~~~~~~~~~~~~
    ! Initialise the parameter list
    call parlst_init(rCHparams)
    call parlst_init(rCHproblem%rparamList)

    ! read in coefficients
    call  parlst_readfromfile(rCHparams, 'data/CahnHilliard.dat')

    ! Initialise the collection.
    call collct_init (rCHproblem%rcollection)
    call collct_setvalue_parlst (rCHproblem%rcollection, 'PARAMS', rCHparams, .TRUE.)
    
    ! Read in the parameters from the DAT file and initialise the basic
    ! structures with these.
	! including time discretisation and spatial discretisation
    call CH_initparameters (rCHparams,rCHproblem)

    ! Add space for level information in the collection
    do i=1,rCHproblem%NLMAX
      call collct_addlevel (rCHproblem%rcollection)
    end do

    ! Allocate memory for the level information
    ALLOCATE (rCHproblem%RlevelInfo(rCHproblem%NLMAX))
    
    ! So now the different steps - one after the other.
    !
    ! Initialisation
    call CH_initParamTriang (rCHproblem%NLMIN,rCHproblem%NLMAX,rCHproblem)

    ! Take care of init discretisation
    call CH_initDiscretisation (rCHproblem)

     !MCai, do we need the following code?
     call CH_initCollectForAssembly (rCHproblem,rCHproblem%rcollection)

    ! Allocate memory first, and create rCHproblem%rhs in terms of discretisation
    call CH_allocMatVec (rCHproblem, rCHparams)

    do i=rCHproblem%NLMIN,rCHproblem%NLMAX
        call CH_generateStaticMatrices (&
            rCHproblem,rCHproblem%RlevelInfo(i))
    end do

! Q: why it does not work? Debug.
    ! first check that rCHproblem%RlevelInfo(NLMAX)%rmatrixLaplace and
	! rCHproblem%RlevelInfo(NLMAX)%rmatrixMass have been assembled.
    ! Debug
!     call lsysbl_scaleVector(rCHproblem%rrhs, 0.0_DP)
!     print *, 'make reuss '
!     call lsyssc_scalarMatVec (rCHproblem%RlevelInfo(rCHproblem%NLMAX)%rmatrixMass, &
! 	 rCHproblem%rrhs%RvectorBlock(1),  rCHproblem%rrhs%RvectorBlock(1), 1.0_DP, 1.0_DP)
!     PRINT *, 'OK?'
!     call lsyssc_scalarMatVec (rCHproblem%RlevelInfo(rCHproblem%NLMAX)%rmatrixLaplace, &
! 	 rCHproblem%rrhs%RvectorBlock(1),  rCHproblem%rrhs%RvectorBlock(1), 1.0_DP, 1.0_DP)
!
!     stop

    ! MCai,
    ! use the auxiliary RHS vector on the finest level to create an
    ! initial RHS and solution vector, which we pass later to the timeloop.
    !MCai, please have a clear picture that rCHrhs exists from here, it is used to
    ! store system rrhs. Also rCHvector: for solution vector.
    call lsysbl_createVecBlockIndirect (rCHproblem%rrhs, rCHrhs, .FALSE.)
    call lsysbl_createVecBlockIndirect (rCHproblem%rrhs, rCHvector,.TRUE.)

    ! MCai, Initialise the initial solution of CH problem.
    call CH_initSolution(rCHproblem, rCHvector)

!MCai, do we need to Initialise postprocessing, How to do this?
! we first check that initial solution is correct.
!    call CH_postprocessing (rCHproblem,rCHvector,0,0.0_DP)
!    stop

    ! Calculate the initial RHS, don't incorporate any BC's.
    call CH_calcRHS (rCHproblem,rCHrhs)

    ! Discretise the boundary conditions
    call CH_initDiscreteBC (rCHproblem)
  
    ! Implement them into the initial solution vector, as we have a zero
    ! vector as initial solution.
    rCHvector%p_rdiscreteBC => rCHproblem%rrhs%p_rdiscreteBC
    call vecfil_discreteBCsol (rCHvector)

!~~~~~~~~~~~~~~~~~~~~~~~~End of Initialization of CH problem~~~~~~~~~~~~~~


!~~~Time loop, we only have time dependent problem~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Now choose the algorithm. Stationary or time-dependent simulation?
    if (p_rproblem%itimedependence .ne. 0) then
      ! Time dependent simulation with explicit time stepping.
      !
      ! Generate the RHS vector for the first time step.
      if (p_rproblem%MSHOW_Initialisation .ge. 1) then
        call output_separator (OU_SEP_MINUS)
        call output_line('Generating RHS vector...')
      end if

!
!Mcai, we should use the following, because at the very beginning, we need rCHvector
!therefore, we need to have the initial solution of rCHvector.
      call cc_generateBasicRHS (p_rproblem, rrhs, &
	           rvector, rCHproblem, rCHvector, rCHrhs)
      
!MCai, debug here, show that the derivatives of rCHvector is 0
!      call CH_calcRHS (rCHproblem, rCHvector, rCHrhs,&
!                          p_rproblem, rvector, rrhs)
!      stop
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! Initialise the boundary conditions, but
      ! don't implement any boundary conditions as the nonstationary solver
      ! doesn't like this.
      if (p_rproblem%MSHOW_Initialisation .ge. 1) then
        call output_separator (OU_SEP_MINUS)
        call output_line('Generating discrete boundary conditions of first time step...')
      end if
      call cc_initDiscreteBC (p_rproblem,rvector,rrhs)
      
      ! Call the nonstationary solver to solve the problem.
      if (p_rproblem%MSHOW_Initialisation .ge. 1) then
        call output_separator (OU_SEP_MINUS)
        call output_line('Invoking nonstationary solver...')
        call output_separator (OU_SEP_MINUS)
      end if
      
      call stat_startTimer(rtimerSolver)

!~~~~~~~We need to modify the code here, solve both NS and Allen-Cahn equation~~
!~~~We may need more input parameters for rCHproblem, rCHvector, rCHrhs,
!~~~~rCHpostprocessing etc.
!      call CHNS_solveNonstationary (p_rproblem,rvector,rrhs,rpostprocessing)

!~~~~~~~~~~Before, we go to next step, stop, make sure all ... are correct~~~
  
! Debug, make sure that there is no problem before time loop.
! Stop
       call CHNS_solveNonstationary (p_rproblem,rvector,rrhs,rpostprocessing, &
                               rCHproblem, rCHvector, rCHrhs, rCHpostprocessing)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call stat_stopTimer(rtimerSolver)
      
    end if
    
!***********************MCai, we need to release memory here**********
    ! release memory from CH problem
    ! Cleanup
    call CH_doneMatVec (rCHproblem)
    call CH_doneBC (rCHproblem)
    call CH_doneDiscretisation (rCHproblem)
    call CH_doneParamTriang (rCHproblem)

    ! Release memory for level information
    DEALLOCATE (rCHproblem%RlevelInfo)
    
    ! Release parameter list
    call collct_deletevalue (rCHproblem%rcollection,'PARAMS')
    call parlst_done (rCHparams)
    call parlst_done(rCHproblem%rparamList)
    
    ! Release RHS and solution vector
    call lsysbl_releaseVector (rCHvector)
    call lsysbl_releaseVector (rCHrhs)
    
    ! Print some statistical data about the collection - anything forgotten?
    call output_lbrk ()
    call output_line ('Remaining collection statistics:')
    call output_line ('--------------------------------')
    call output_lbrk ()
    call collct_printStatistics (rCHproblem%rcollection)
    
    ! Finally release the collection.
    call collct_done (rCHproblem%rcollection)

!~~~~~~~~~~Release memory of NS problem~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Gather statistics
    p_rproblem%rstatistics%dtimeSolver = &
      p_rproblem%rstatistics%dtimeSolver + rtimerSolver%delapsedReal
    
    ! (Probably) write final solution vector
    call cc_writeSolution (p_rproblem,rvector)
    
    ! Cleanup
    call cc_doneMatVec (p_rproblem,rvector,rrhs)
    call cc_doneBC (p_rproblem)
    call cc_doneDiscretisation (p_rproblem)
    call cc_donepostprocessing (rpostprocessing)
    call cc_doneParamTriang (p_rproblem)
    
    ! Release parameters from the DAT/INI files from the problem structure.
    call cc_doneParameters (p_rproblem)

    ! Release the parameter list
    call parlst_done (p_rproblem%rparamList)
    
    ! Print some statistical data about the collection - anything forgotten?
    call output_lbrk ()
    call output_line ('Remaining collection statistics:')
    call output_line ('--------------------------------')
    call output_lbrk ()
    call collct_printStatistics (p_rproblem%rcollection)
    
    ! Finally release the collection and the problem structure.
    call collct_done (p_rproblem%rcollection)
    
    ! Stop the timer
    call stat_stopTimer(rtimerTotal)

    ! Gather statistics
    p_rproblem%rstatistics%dtimeTotal = &
      p_rproblem%rstatistics%dtimeTotal + rtimerTotal%delapsedReal
    
    ! Print the time for the total computation
    call output_lbrk ()
    call output_line ("Total time:                             "//&
        trim(sys_sdL(p_rproblem%rstatistics%dtimeTotal,10)))

    call output_line ("Time for initial mesh generation:       "//&
        trim(sys_sdL(rtimerGridGeneration%delapsedReal,10)))
        
    call output_line ("Time for initial matrix assembly:       "//&
      trim(sys_sdL(rtimerMatrixGeneration%delapsedReal,10)))
      
    call output_line ("Total time for grid generation:         "//&
      trim(sys_sdL(p_rproblem%rstatistics%dtimeGridGeneration,10)))
      
    call output_line ("Total time for complete solver:         "//&
      trim(sys_sdL(p_rproblem%rstatistics%dtimeSolver,10)))
      
    call output_line ("Total time for nonlinear solver:        "//&
      trim(sys_sdL(p_rproblem%rstatistics%dtimeNonlinearSolver,10)))
      
    call output_line ("Total time for defect calculation:      "//&
      trim(sys_sdL(p_rproblem%rstatistics%dtimeDefectCalculation,10)))

    call output_line ("Total time for optimal damping:         "//&
      trim(sys_sdL(p_rproblem%rstatistics%dtimeOptimalCorrection,10)))
      
    call output_line ("Total time for matrix assembly:         "//&
      trim(sys_sdL(p_rproblem%rstatistics%dtimeMatrixAssembly,10)))
      
    call output_line ("Total time for linear solver:           "//&
      trim(sys_sdL(p_rproblem%rstatistics%dtimeLinearSolver,10)))
      
    call output_line ("Total time for factorisation:           "//&
      trim(sys_sdL(p_rproblem%rstatistics%dtimeLinearSolverFactorisation,10)))

    call output_line ("Total time for postprocessing:          "//&
      trim(sys_sdL(p_rproblem%rstatistics%dtimePostprocessing,10)))
      
    call output_line ("Total #iterations nonlinear solver:     "//&
      trim(sys_siL(p_rproblem%rstatistics%nnonlinearIterations,10)))
      
    call output_line ("Total #iterations linear solver:        "//&
      trim(sys_siL(p_rproblem%rstatistics%nlinearIterations,10)))

    call output_line ("Total number of calculated timesteps:   "//&
      trim(sys_siL(p_rproblem%rstatistics%ntimesteps,10)))
    ! That's it.
    deallocate(p_rproblem)

!~~~~~~~~~~~~~~~~~~~~~~End of postprocessing of NS Equation~~~~~~~~~~~~~~~~~

  end subroutine

end module
