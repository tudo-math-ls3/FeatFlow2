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

module cc2dmediumm2basic

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
  use timestepping
  use linearsystemscalar
  use linearsystemblock
  use discretebc
  use discretefbc

  use collection
  
  use adaptivetimestep
    
  implicit none
  
  ! Maximum allowed level in this application; must be =9 for 
  ! FEAT 1.x compatibility (still)!
  integer, parameter :: NNLEV = 9
  
!<types>

!<typeblock description="Type block defining all information about one level">

  type t_problem_lvl
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the block discretisation
    ! (size of subvectors in the solution vector, trial/test functions,...)
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! A template FEM matrix that defines the structure of Laplace/Stokes/...
    ! matrices. The matrix contains only a stucture, no content.
    type(t_matrixScalar) :: rmatrixTemplateFEM

    ! A template FEM matrix that defines the structure of gradient
    ! matrices (B1/B2) matrices. The matrix contains only a stucture, no content.
    type(t_matrixScalar) :: rmatrixTemplateGradient
    
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
    type(t_matrixBlock) :: rpreallocatedSystemMatrix
    
    ! Stokes matrix for that specific level (=nu*Laplace)
    type(t_matrixScalar) :: rmatrixStokes
    
    ! B1-matrix for that specific level. 
    type(t_matrixScalar) :: rmatrixB1

    ! B2-matrix for that specific level. 
    type(t_matrixScalar) :: rmatrixB2

    ! Temporary vector in the size of the RHS/solution vector on that level.
    type(t_vectorBlock) :: rtempVector

    ! A variable describing the discrete boundary conditions fo the velocity
    type(t_discreteBC), pointer :: p_rdiscreteBC
  
    ! A structure for discrete fictitious boundary conditions
    type(t_discreteFBC), pointer :: p_rdiscreteFBC
    
    ! Nonstationary simulation: Mass matrix
    type(t_matrixScalar) :: rmatrixMass

    ! Nonstationary simulation: A scalar discretisation structure that 
    ! specifies how to generate the mass matrix.
    type(t_spatialDiscretisation), pointer :: p_rdiscretisationMass

    ! This flag signales whether there are Neumann boundary components
    ! visible on the boundary of this level or not. If there are no
    ! Neumann boundary components visible, the equation gets indefinite
    ! for the pressure.
    logical :: bhasNeumannBoundary
    
  end type
  
!</typeblock>


!<typeblock description="Application-specific type block for the nonstationary Nav.St. problem">

  type t_problem_explTimeStepping
  
    ! Number of current time step; changes during the simulation.
    integer :: itimeStep           = 0
    
    ! Current simulation time; changes during the simulation.
    real(DP) :: dtime              = 0.0_DP
  
    ! Maximum number of time steps; former NITNS
    integer :: niterations         = 0
    
    ! Absolute start time of the simulation
    real(DP) :: dtimeInit          = 0.0_DP     
    
    ! Time step size; former TSTEP
    real(DP) :: dtimeStep          = 0.0_DP       
    
    ! Maximum time of the simulation
    real(DP) :: dtimeMax           = 0.0_DP
  
    ! Lower limit for the time derivative to be treated as zero. Former EPSNS.
    ! Simulation stops if time derivative drops below this value.
    real(DP) :: dminTimeDerivative = 0.00001_DP
    
    ! Configuration block for the adaptive time stepping.
    type(t_adaptimeTimeStepping) :: radaptiveTimeStepping
    
  end type

!</typeblock>


!<typeblock description="Application-specific type block for the Nav.St. problem">

  type t_problem
  
    ! Output level during the initialisation phase.
    integer :: MSHOW_Initialisation
  
    ! Output level of the application.
    integer :: MT_OutputLevel
  
    ! Minimum refinement level; = Level i in RlevelInfo
    integer :: NLMIN
    
    ! Maximum refinement level
    integer :: NLMAX
    
    ! Viscosity parameter nu = 1/Re
    real(DP) :: dnu
    
    ! Type of problem.
    ! =0: Stokes.
    ! =1: Navier-Stokes.
    integer :: iequation
    
    ! Type of subproblem of the main problem. Depending on iequationType.
    ! If iequationType=0 or =1:
    ! =0: (Navier-)Stokes with gradient tensor
    ! =1: (Navier-)Stokes with deformation tensor
    integer :: isubEquation
        
    ! An object for saving the domain:
    type(t_boundary), pointer :: p_rboundary

    ! Flag indicating if the X- and Y-velocity is decoupled (i.e. yield different 
    ! matrices). This is the case e.g. for no-slip boundary conditions
    ! where then implementation of the BC's into the first velocity
    ! matrix must not affect the 2nd velocity matrix!
    ! This value must be initialised before the matrices are set
    ! up and must not be changed afterwards.
    logical :: bdecoupledXY

    ! iboundary parameter from the DAT file. Specifies whether to update
    ! boundary conditions in nonstationary simulations or not.
    ! =0: stationary Dirichlet boundary conditions
    ! =1: nonstationary boundary conditions, possibly with
    !     pressure drop
    ! =2: nonstationary boundary conditions, possibly with
    !     prssure drop and/or Neumann boundary parts
    integer :: iboundary

    ! A solver node that accepts parameters for the linear solver    
    type(t_linsolNode), pointer           :: p_rsolverNode

    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation. There is currently
    ! only one level supported, identified by NLMAX!
    type(t_problem_lvl), dimension(NNLEV) :: RlevelInfo
    
    ! Type of simulation.
    ! =0: stationary simulation.
    ! =1: time-dependent simulation with explicit time stepping configured 
    !     by rtimedependence
    integer                               :: itimedependence
    
    ! A parameter block for everything that controls the time dependence.
    ! Only valid if itimedependence=1!
    type(t_problem_explTimeStepping)      :: rtimedependence
    
    ! A collection object that saves structural data and some 
    ! problem-dependent information which is e.g. passed to 
    ! callback routines.
    type(t_collection)                    :: rcollection
    
    ! A param list that saves all parameters from the DAT/INI file(s).
    type(t_parlist)                       :: rparamList

  end type

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
! This structure contains the following information, which is added by the
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
! UPSAM                 | Stabilisation parameter
!
! Before the nonlinear iteration, some parameters are added to the collection.
! These are available during the nonlinear iteration and released afterwards.
! The routine c2d2_saveNonlinearLoop saves all these parameters while
! they can be restored from the collection by c2d2_restoreNonlinearLoop.
! c2d2_removeNonlinearLoop finally deletes these entries from the collection
! again.
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

end module
