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
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation

    ! An object specifying the block discretisation structure only for the
    ! primal space.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisationPrimal

    ! An object specifying the block discretisation structure only for the
    ! dual space.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisationDual

    ! A template FEM matrix that defines the structure of Laplace/Stokes/...
    ! matrices. The matrix contains only a stucture, no content.
    TYPE(t_matrixScalar) :: rmatrixTemplateFEM

    ! A template FEM matrix that defines the structure of gradient
    ! matrices (B1/B2) matrices. The matrix contains only a stucture, no content.
    TYPE(t_matrixScalar) :: rmatrixTemplateGradient
    
    ! A template matrix for the system matrix for that specific level.
    ! Provides memory for intermediate calculations of the system matrix.
    ! This matrix is used e.g. as template matrix for a Multigrid preconditioner
    ! during the nonlinear iteration; in this case, the 'preconditioner matrices'
    ! on all levels share some information with this to prevent frequent 
    ! reallocation of memory. On the other hand, the matrix might have to be
    ! evaluated for some reason (e.g. the evaluation of damping parameters)
    ! which can be done with this variable to avoid memory allocation.
    !
    ! The system matrix at the same time defines on the one hand the shape
    ! of the global system (number of DOF's, submatrices for gradient and/or
    ! deformation tensor). On the other hand, the boundary conditions are associated
    ! to this matrix to allow the implementation of boundary conditions into
    ! a globally assembled system matrix.
    !
    ! Note that the system matrix does not have to be assembled for calculating
    ! the defect! Routines to assemble the system matrix or the defect
    ! can be found in the module cc2dmediumm2matvecassembly.
    TYPE(t_matrixBlock) :: rpreallocatedSystemMatrix
    
    ! This is a reference to rpreallocatedSystemMatrix(1:3,1:3), which realises
    ! the pure primal system.
    TYPE(t_matrixBlock) :: rpreallocatedSystemMatrixPrimal
    
    ! This is a reference to rpreallocatedSystemMatrix(4:6,4:6), which realises
    ! the pure dual system.
    TYPE(t_matrixBlock) :: rpreallocatedSystemMatrixDual
    
    ! Stokes matrix for that specific level (=nu*Laplace)
    TYPE(t_matrixScalar) :: rmatrixStokes
    
    ! An identity matrix in the size of the pressure
    TYPE(t_matrixScalar) :: rmatrixIdentityPressure
    
    ! B1-matrix for that specific level. 
    TYPE(t_matrixScalar) :: rmatrixB1

    ! B2-matrix for that specific level. 
    TYPE(t_matrixScalar) :: rmatrixB2

    ! Temp vector for the full system.
    TYPE(t_vectorBlock) :: rtempVector
    
    ! Reference ro rtempVector(1:3), which corresponds to the primal solution.
    TYPE(t_vectorBlock) :: rtempVectorPrimal

    ! Reference ro rtempVector(4:6), which corresponds to the dual solution.
    TYPE(t_vectorBlock) :: rtempVectorDual

    ! A variable describing the discrete boundary conditions for the velocity
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC
  
    ! A structure for discrete fictitious boundary conditions
    TYPE(t_discreteFBC), POINTER :: p_rdiscreteFBC
    
    ! A variable describing the discrete boundary conditions for the primal velocity
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBCprimal

    ! A variable describing the discrete boundary conditions for the dual velocity
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBCdual
  
    ! A structure for discrete fictitious boundary conditions (primal solution)
    TYPE(t_discreteFBC), POINTER :: p_rdiscreteFBCprimal

    ! A structure for discrete fictitious boundary conditions (dual solution)
    TYPE(t_discreteFBC), POINTER :: p_rdiscreteFBCdual
    
    ! Mass matrix
    TYPE(t_matrixScalar) :: rmatrixMass

    ! Nonstationary simulation: A scalar discretisation structure that 
    ! specifies how to generate the mass matrix.
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisationMass
    
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
  
    ! Type of implementation of the terminal condition.
    ! =0: implement terminal condition in a weak sense by filtering.
    ! =1: implement terminal condition in a strong sense by modifying the matrix.
    INTEGER :: itypeTerminalCondition = 0
  
    ! Name of the file with the target flow
    CHARACTER(SYS_STRLEN) :: stargetFlow = ''
    
    ! Refinement level of the target flow.
    INTEGER :: ilevelTargetFlow = 0
    
    ! Underlying mesh for the target flow if itypeTargetFlow<>0.
    TYPE(t_triangulation), POINTER :: p_rtriangulation => NULL()
    
    ! Discretisation structure specifying the discretisation of the
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
    TYPE(t_boundary), POINTER :: p_rboundary
    
    ! Flag indicating if the X- and Y-velocity is decoupled (i.e. yield different 
    ! matrices). This is the case e.g. for no-slip boundary conditions
    ! where then implementation of the BC's into the first velocity
    ! matrix must not affect the 2nd velocity matrix!
    ! This value must be initialised before the matrices are set
    ! up and must not be changed afterwards.
    LOGICAL :: bdecoupledXY

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
! INI                   | t_paramlist object
!                       | Contains parameters from the DAT/INI files
!                       |
! NU                    | Reciprocal 1/RE of parameter RE from the DAT file
! NLMIN                 | Minimum level of the discretisation
! NLMAX                 | Maximum level of the discretisation
! IUPWIND               | Type of stabilisation. =0:streamline diff, =1:upwind,
!                       | =2:jump stabil
! IBOUNDARY             ! =0: stationary Dirichlet boundary conditions
!                       ! =1: nonstationary boundary conditions, possibly with
!                       !     pressure drop
!                       ! =2: nonstationary boundary conditions, possibly with
!                       !     prssure drop and/or Neumann boundary parts
! UPSAM                 | Stabilisation parameter
!
! On every level between NLMIN and NLMAX:
!
! Name                  | Description
! ----------------------+------------------------------------------------------
! STOKES                | Stokes matrix (=nu*Laplace) if nu=constant
! SYSTEMMAT             | Nonlinear system matrix
! RTEMPVEC              | Temporary vector, compatible to matrix
!
!
! Global, level independent data, available during the nonlinear iteration:
!
! Name                  | Description
! ----------------------+------------------------------------------------------
! RHS                   | t_vectorBlock object 
!                       | Current RHS vector on maximum level
!                       |
! SOLUTION              | t_vectorBlock object 
!                       | Current solution vector on maximum level
!                       |
! ILVPROJECTION         | t_interlevelProjectionBlock structure.
!                       | Configures prolongation/restriction for multigrid
!                       | solver component.
!                       |
! RTEMPSCALAR           | t_vectorScalar object 
!                       | Temporary vector
!                       |
! RTEMP2SCALAR          | t_vectorScalar object 
!                       | Temporary vector
!                       |
! LINSOLVER             | t_linsol object
!                       | Configures the linear solver for preconditioning
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
