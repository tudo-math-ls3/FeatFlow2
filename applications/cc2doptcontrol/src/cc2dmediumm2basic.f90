!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2basic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic problem definitions for the cc2dmini_method2
!# solver. The basic structure and content of the different structures
!# are described here.
!# </purpose>
!##############################################################################

MODULE cc2dmediumm2basic

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
  USE timestepping
  USE timediscretisation
  USE spacetimevectors
  
  USE collection
  
  IMPLICIT NONE
  
!<types>

!<typeblock description="Type block defining all information about one level">

  TYPE t_problem_lvl
  
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation) :: rtriangulation

    ! An object specifying the block discretisation
    ! (size of subvectors in the solution vector, trial/test functions,...)
    TYPE(t_blockDiscretisation) :: rdiscretisation

    ! An object specifying the block discretisation structure only for the
    ! primal space.
    TYPE(t_blockDiscretisation) :: rdiscretisationPrimal

    ! A template FEM matrix that defines the structure of Laplace/Stokes/...
    ! matrices. The matrix contains only a stucture, no content.
    TYPE(t_matrixScalar) :: rmatrixTemplateFEM

    ! A template FEM matrix that defines the structure of gradient
    ! matrices (B1/B2) matrices. The matrix contains only a stucture, no content.
    TYPE(t_matrixScalar) :: rmatrixTemplateGradient
    
    ! Stokes matrix for that specific level (=nu*Laplace)
    TYPE(t_matrixScalar) :: rmatrixStokes
    
    ! An identity matrix in the size of the pressure
    TYPE(t_matrixScalar) :: rmatrixIdentityPressure
    
    ! B1-matrix for that specific level. 
    TYPE(t_matrixScalar) :: rmatrixB1

    ! B2-matrix for that specific level. 
    TYPE(t_matrixScalar) :: rmatrixB2

    ! Three temp vectors for the full system.
    TYPE(t_vectorBlock) :: rtempVector1
    TYPE(t_vectorBlock) :: rtempVector2
    TYPE(t_vectorBlock) :: rtempVector3
    
    ! Reference ro rtempVector(1:3), which corresponds to the primal solution.
    TYPE(t_vectorBlock) :: rtempVectorPrimal

    ! Reference ro rtempVector(4:6), which corresponds to the dual solution.
    TYPE(t_vectorBlock) :: rtempVectorDual

    ! A variable describing the discrete boundary conditions for the velocity.
    ! Points to NULL() until the BC's are discretised for the first time.
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC => NULL()
  
    ! A structure for discrete fictitious boundary conditions
    ! Points to NULL() until the BC's are discretised for the first time.
    TYPE(t_discreteFBC), POINTER :: p_rdiscreteFBC => NULL()
    
    ! Mass matrix
    TYPE(t_matrixScalar) :: rmatrixMass

    ! Nonstationary simulation: A scalar discretisation structure that 
    ! specifies how to generate the mass matrix.
    TYPE(t_spatialDiscretisation) :: rdiscretisationMass
    
    ! This flag signales whether there are Neumann boundary components
    ! visible on the boundary of this level or not. If there are no
    ! Neumann boundary components visible, the equation gets indefinite
    ! for the pressure.
    LOGICAL :: bhasNeumannBoundary
    
  END TYPE
  
!</typeblock>


!<typeblock>
  
  ! Application-specific type block for the nonstationary Nav.St. problem
  TYPE t_problem_explTimeStepping
  
    ! Number of current time step; changes during the simulation.
    INTEGER :: itimeStep           = 0
    
    ! Current simulation time; changes during the simulation.
    REAL(DP) :: dtime              = 0.0_DP
  
    ! Number of time steps; former NITNS
    INTEGER :: niterations         = 0
    
    ! Absolute start time of the simulation
    REAL(DP) :: dtimeInit          = 0.0_DP     
    
    ! Maximum time of the simulation
    REAL(DP) :: dtimeMax           = 0.0_DP
    
    ! Time-stepping scheme (IFRSTP);
    ! 0=one step scheme, 1=fractional step, 2=dG(0)
    INTEGER :: ctimeStepScheme     = 0
    
    ! Parameter for one step scheme (THETA) if ctimeStepScheme=0;
    ! =0:Forward Euler(instable), =1: Backward Euler, =0.5: Crank-Nicolson
    REAL(DP) :: dtimeStepTheta     = 1.0_DP
    
  END TYPE

!</typeblock>


!<typeblock>

  ! This block saves parameters of the optimal control problem
  TYPE t_problem_optcontrol
  
    ! $\alpha$ parameter of the optimal control functional
    REAL(DP) :: dalphaC = 1.0_DP
    
    ! $\gamma$ parameter of the nonstationary optimal control functional
    REAL(DP) :: dgammaC = 0.0_DP
  
    ! Type of target flow.
    ! =0: analytically given.
    ! =1: stationary solution, read from file.
    ! =2: nonstationary solution, read from a sequence of files.
    ! =3: stationary solution, read from file. Arbitrary level.
    ! =4: nonstationary solution, read from a sequence of files. Arbitrary level.
    INTEGER :: itypeTargetFlow = 0
    
    ! Formulation of the Space-time problem.
    ! =0: usual formulation as specified in the DFG applicance
    ! =1: Formulation for the generation of reference results from papers
    ! The two formulations differ in a "-"-sign in front of the dual velocity.
    INTEGER :: ispaceTimeFormulation = 0
  
    ! Whether to treat the convection explicitly or implicitly.
    ! =0: Treat the convection implicitely.
    ! =1: Treat the convection explicitely. (This is the formulation of the paper
    !     of Baerwolff and Hinze!)
    INTEGER :: iconvectionExplicit = 0

    ! Type of implementation of the terminal condition.
    ! =0: implement terminal condition in a weak sense by filtering.
    ! =1: implement terminal condition in a strong sense by modifying the matrix.
    INTEGER :: itypeTerminalCondition = 0
  
    ! Name of the file with the target flow
    CHARACTER(SYS_STRLEN) :: stargetFlow = ''
    
    ! Refinement level of the target flow.
    INTEGER :: ilevelTargetFlow = 0
    
    ! Name of the TRI file with the mesh corresponding to target flow.
    ! ='': Use the same mesh as for the computation of the solution
    CHARACTER(SYS_STRLEN) :: smeshTargetFlow = ''

    ! Underlying mesh for the target flow if itypeTargetFlow<>0.
    TYPE(t_triangulation), POINTER :: p_rtriangulationTargetFlow => NULL()
    
    ! Discretisation structure specifying the discretisation of the
    ! target flow.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscrTargetFlow => NULL()
    
    ! Solution vector containing a stationary target flow.
    ! Only used if itypeTargetFlow <> 0.
    TYPE(t_vectorBlock) :: rtargetFlow
  
    ! Time discretisation for a nonstationary target flow.
    ! Only used if itypeTargetFlow = 2 or 4
    TYPE(t_timeDiscretisation) :: rtargetTimeDiscr
    
    ! Solution vector containing a nonstationary target flow.
    ! Only used if itypeTargetFlow = 2 or 4
    TYPE(t_spaceTimeVector) :: rtargetFlowNonstat
    
    ! Target flow granilarity for the filenames of a nonstationary target flow.
    ! Specifies how the number in the filename is increased.
    INTEGER :: itargetFlowDelta = 1
    
    ! Number of files that describe the target flow.
    ! -1=As many files as there are solution vectors on the finest time level.
    INTEGER :: itargetFlowTimesteps = -1
  
    ! Type of constraints to apply to the control u.
    ! =0: No constraints.
    ! =1: Constant constraints on u active: 
    !     dumin1 <= u_1 <= dumax1, dumin2 <= u_2 <= dumax2
    INTEGER :: ccontrolContraints = 0

    ! Constraints on u_1
    real(DP) :: dumin1 = -1.0E10
    real(DP) :: dumax1 = 1.0E10

    ! Constraints in u_2
    real(DP) :: dumin2 = -1.0E10
    real(DP) :: dumax2 = 1.0E10
  
  END TYPE

!</typeblock>


!<typeblock>
  
  ! Structure with some pointers to internal variables of a special algorithm.
  ! Used to provide data to the parser so the user can advice the program to do
  ! something with this data.
  ! Only valid during the execution of the command parser.
  TYPE t_problem_oneshot
  
    ! Current global iteration
    INTEGER :: iglobIter
  
    ! Norm of the defect vector
    REAL(DP) :: ddefNorm
    
    ! Pointer to the current space-time solution vector
    TYPE(t_spaceTimeVector), POINTER :: p_rx

    ! Pointer to the current space-time RHS vector
    TYPE(t_spaceTimeVector), POINTER :: p_rb
    
  END TYPE

!</typeblock>


!<typeblock>

  ! Application-specific type block for the Nav.St. problem
  TYPE t_problem
  
    ! Output level during the initialisation phase.
    INTEGER :: MSHOW_Initialisation
  
    ! Output level of the application.
    INTEGER :: MT_OutputLevel
  
    ! Minimum refinement level; = Level i in RlevelInfo
    INTEGER :: NLMIN
    
    ! Maximum refinement level
    INTEGER :: NLMAX
    
    ! Viscosity parameter nu = 1/Re
    REAL(DP) :: dnu
    
    ! Type of problem.
    ! =0: Navier-Stokes.
    ! =1: Stokes.
    INTEGER :: iequation
    
    ! Type of subproblem of the main problem. Depending on iequationType.
    ! If iequationType=0 or =1:
    ! =0: (Navier-)Stokes with gradient tensor
    ! =1: (Navier-)Stokes with deformation tensor
    INTEGER :: isubEquation
        
    ! An object for saving the domain:
    TYPE(t_boundary) :: rboundary
    
    ! A variable describing the analytic boundary conditions.    
    TYPE(t_boundaryConditions), POINTER   :: p_rboundaryConditions

    ! A variable describing the analytic boundary conditions for the primal
    ! equation.
    TYPE(t_boundaryConditions), POINTER   :: p_rboundaryConditionsPrimal

    ! A variable describing the analytic boundary conditions for the dual
    ! equation.
    TYPE(t_boundaryConditions), POINTER   :: p_rboundaryConditionsDual

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER           :: p_rsolverNode

    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation. There is currently
    ! only one level supported, identified by NLMAX!
    TYPE(t_problem_lvl), DIMENSION(:), POINTER :: RlevelInfo
    
    ! Type of simulation.
    ! =0: stationary simulation.
    ! =1: time-dependent simulation with explicit time stepping configured 
    !     by rtimedependence
    INTEGER                               :: itimedependence
    
    ! A parameter block for everything that controls the time dependence.
    ! Only valid if itimedependence=1!
    TYPE(t_problem_explTimeStepping)      :: rtimedependence
    
    ! A parameter block for everything that controls the optimal control
    ! problem.
    TYPE(t_problem_optcontrol)            :: roptcontrol
    
    ! A collection object that saves structural data and some 
    ! problem-dependent information which is e.g. passed to 
    ! callback routines.
    TYPE(t_collection)                    :: rcollection
    
    ! A param list that saves all parameters from the DAT/INI file(s).
    TYPE(t_parlist)                       :: rparamList

    ! A t_problem_oneshot structure that is only valid during the
    ! execution of the command parser.
    TYPE(t_problem_oneshot)               :: rdataOneshot

  END TYPE

!</typeblock>

!</types>

!<globals>

  ! Directory containing the data files.
  CHARACTER(LEN=SYS_STRLEN), SAVE :: DIR_DATA = "./data";

!</globals>

!******************************************************************************
! Documentation
!******************************************************************************

!******************************************************************************
! The problem structure t_problem
!******************************************************************************
!
! The problem structure t_problem collects all information of importance
! for the main CC2D solver module. It contains:
!  - Analytic definition of the boundary
!  - Analytic boundary conditions
!  - Information from the INI/DAT files
!  - Main solution and RHS vector
!  - Information for preconditioning in nonlinear loop
!  - For every level in the discretisation:
!    - Triangulation
!    - Block discretisation
!    - Matrices and
!    - Temporary vectors
! This problem structure is available only in the top-level routines of
! the CC2D module; it's not available in any callback routines.
! Parameters for callback routines are passed through the collection
! structure rcollection, which is part of the problem structure.
!
!
!******************************************************************************
! The collection structure t_problem%rcollection
!******************************************************************************
!
! The collection structure collects all information that is necessary to be
! passed to callback routines. This e.g. allows to pass matrices/vectors/
! constants from the main problem.
! This structure contains the following informaion, which is added by the
! initialisation routines to the collection:
!
! Global, level independent data:
!
! Name                  | Description
! ----------------------+------------------------------------------------------
! NU                    | Reciprocal 1/RE of parameter RE from the DAT file
! IBOUNDARY             ! =0: stationary Dirichlet boundary conditions
!                       ! =1: nonstationary boundary conditions, possibly with
!                       !     pressure drop
!                       ! =2: nonstationary boundary conditions, possibly with
!                       !     prssure drop and/or Neumann boundary parts
! UPSAM                 | Stabilisation parameter
!
! On every level between NLMIN and NLMAX:
!
! Expressions for boundary conditions are saved in special 
! sections in the collection:
!
! Section [BDEXPRESSIONS]:
! ------------------------
! Saves the expressions that are to be evaluated on the boundary.
!
! During the assembly of RHS vectors and boundary conditions, the following
! additional information is valid in the collection:
!
!   collection%IquickAccess(1)   = 0: stationary, 
!                                  1: nonstationary with explicit time stepping
!   collection%DquickAccess(1)   = current simulation time
!   collection%DquickAccess(2)   = minimum simulation time
!   collection%DquickAccess(3)   = maximum simulation time
  
END MODULE
