!##############################################################################
!# ****************************************************************************
!# <name> spacematvecassembly </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the very basic matrix assembly routines for the
!# core equation. It's independent of any nonlinear iteration and provides
!# just one functionality: Assemble a matrix or a vector based on a given
!# set of parameters.
!#
!# The discretised core equation reads at the moment:
!#
!#  $$           A_1 y   +  \eta_1 B p   +          R lambda      = f_1 $$
!#  $$ \tau_{11} B^T y   +  \kappa_{11} I p                       = f_2 $$
!#
!#  $$ (R_{22}+R2_{22}) y +  A_{22} \lambda  + \eta_{22}   B \xi           = f_3 $$
!#  $$               \tau_{22} B^T \lambda  + \kappa_{22} I \xi           = f_4 $$
!#
!# with
!#
!#   $$ A_1 = \iota_{11} I  +  \alpha_{11} M  +  \theta_{11} L  +
!#            \gamma_{11} N(y) + dnewton_{11} N* + gammaT_{11} N^t(y) + dnewtonT_{11} N*^t(y)(y)$$
!#   $$ A_2 = \iota_{22} I  +  \alpha_{22} M  +  \theta_{22} L  +
!#            \gamma_{22} N(y) + dnewton_2 N*(y) + \gammaT_{22} N^t(y) + dnewtonT_{22} N*^t(y)$$
!#   $$ R_1 = \alpha_{12} (P)M(\lambda) $$
!#   $$ R_2 = \alpha_{21} M  +
!#            \gamma_{21} N(\lambda) + dnewton_{21} N*(\lambda) +
!#            gammaT_{21} N^t(\lambda) + dnewtonT_{21} N*^t(\lambda) $$
!#   $$ R2_2 =  dnewton_{21} N*(\lambda) + gammaT_{21} N^t(\lambda) $$
!#
!# and
!#
!#   $I$     = identity matrix,
!#   $M$     = mass matrix,
!#   $L$     = Stokes matrix ($\nu$*Laplace),
!#   $N(y)$  = Nonlinearity including stabilisation, depending on the
!#             primal velocity, i.e.
!#                   $$ (y\Delta)\cdot $$
!#   $N*(y)$ = Newton matrix, depending on the primal velocity, i.e.
!#                  $$ (\Delta y)\cdot $$
!#   $N^t(y)$ = Transposed nonlinearity including stabilisation, depending on the
!#              primal velocity, i.e.
!#                   $$ (y\Delta)^t\cdot $$
!#   $N*^t(y)$ = Transposed Newton matrix, depending on the primal velocity, i.e.
!#                  $$ (\Delta y)^t\cdot $$
!#
!#   $\iota_{ij}$  = 0/1     - switches the identity matrix on/off,
!#   $\alpha_{ij}$ = 0/1     - switches the mass matrices on/off;
!#                             =0 for stationary problem,
!#   $\theta_{íj}$           - weight for the Laplace matrix,
!#   $\gamma_{ij}$ = 0/1     - Switches the nonlinearity on/off;
!#                             =0 for Stokes system,
!#   $dnewton_{ij} \in R$    - Switches the Newton matrix on/off.
!#   $\eta_{ij}$   = 0/1     - Switches the 'B'-term on/off,
!#   $\tau_{ij}$   = 0/1     - Switches the 'B^T'-term on/off,
!#   $\kappa_{ij}$ = 0/1     - Switches of the identity matrix I for the pressure
!#                             in the continuity equation
!#
!# Note that the nonlinear part is always dependent on the primal velocity --
!# for the primal equation as well as for the dual one!
!#
!# To assemble such a matrix, the application has to follow two steps:
!#
!# 1.) Create a structure of type t_nonlinearSpatialMatrix and set the
!#     parameters in this structure according to the matrix which should
!#     be assembled.
!#
!# 2.) Call the matrix/vector assembly routine with this structure as input.
!#
!# The matrix assembly routine will return an initialised matrix structure
!# of a system matrix based on the provided parameters.
!#
!# The module contains the following routines:
!#
!# 1.) smva_assembleMatrix
!#     -> Assembles a matrix based on a set of input parameters.
!#
!# 2.) smva_assembleDefect
!#     -> Set up a defect vector d:=b-A(x)x
!#
!# 3.) smva_initNonlinMatrix
!#     -> Initialises a nonlinear-matrix structure with basic parameters.
!#
!# 4.) smva_initNonlinearData
!#     -> Initialises a nonlinear-data structure.
!#
!# 5.) smva_clearMatrix
!#     -> Clears all weights in a matrix.
!#
!# 6.) smva_disableSubmatrix
!#     -> Disables a submatrix by setting the corresponding weights to zero.
!#
!# 7.) smva_clearMatrix
!#     -> Sets all weights in a matrix to zero
!#
!# 8.) smva_prepareViscoAssembly
!#     -> Prepare a collection for the use in ffunctionViscoModel
!#
!# 9.) ffunctionViscoModel
!#     -> Auxiliary function that defines the nonconstant viscosity
!# </purpose>
!##############################################################################

module spacematvecassembly

  use fsystem
  use storage
  use genoutput
  use basicgeometry
  use boundary
  use cubature
  use elementpreprocessing
  use element
  use derivatives
  use matrixfilters
  use vectorfilters
  use bcassemblybase
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use scalarpde
  use linearsystemscalar
  use linearsystemblock
  use linearsolver
  use transformation
  use bilinearformevaluation
  use linearformevaluation
  use linearsolverautoinitialise
  use matrixrestriction
  use trilinearformevaluation
  use matrixmodification
  use matrixio
  use feevaluation
  use convection
  use collection
  use mprimitives
  use dofmapping
  use domainintegration
  use spatialdiscretisation
  use timediscretisation
  use discretebc
  use bcassembly
  use matrixfilters
  
  use constantsoptc
  use assemblytemplates
  use assemblytemplatesoptc
  use structuresoptc
  
  use structuresoptflow
  use user_callback

  use spatialbcdef

  use optcontrolconvection
  use newtonderivative
    
  implicit none
  
  private
  
!<constants>

!<constantblock description="Identifiers for the 'coperation' input parameter of the matrix assembly routine">

  ! Allocate memory if necessary.
  integer(I32), parameter, public :: CCMASM_ALLOCMEM              = 1
  
  ! Compute all matrix entries.
  integer(I32), parameter, public :: CCMASM_COMPUTE               = 2
  
  ! Allocate memory and compute matrix entries.
  integer(I32), parameter, public :: CCMASM_ALLOCANDCOMPUTE       = 3
  
  ! Bypass memory allocation for matrices.
  integer(I32), parameter, public :: CMASM_QUICKREFERENCES        = 4
  
  ! Don't clear the old matrix
  integer(i32), parameter, public :: CMASM_NOCLEAR                = 8
  
!</constantblock>

!<constantblock description="Identifiers for the IUPWIND parameter that specifies how to set up the nonlinearity or stabilisation.">

  ! Streamline diffusion; configured by dupsam
  integer, parameter, public :: CCMASM_STAB_STREAMLINEDIFF    = 0

  ! 1st-order upwind; configured by dupsam
  integer, parameter, public :: CCMASM_STAB_UPWIND            = 1
  
  ! Edge-oriented stabilisation; configured by dupsam as 'gamma'
  integer, parameter, public :: CCMASM_STAB_EDGEORIENTED      = 2

  ! Streamline diffusion; configured by dupsam, new implementation
  integer, parameter, public :: CCMASM_STAB_STREAMLINEDIFF2   = 3

  ! Edge-oriented stabilisation; configured by dupsam as 'gamma', new implementation
  integer, parameter, public :: CCMASM_STAB_EDGEORIENTED2     = 4

  ! Edge-oriented stabilisation; configured by dupsam as 'gamma', new implementation.
  ! Precomputed matrix.
  integer, parameter, public :: CCMASM_STAB_EDGEORIENTED3     = 5
!</constantblock>

!<constantblock description="Matrix type ID's specifying the general matrix class to set up.">

  ! Standard matrix.
  integer, parameter, public :: CCMASM_MTP_AUTOMATIC         = 0
  
  ! Standard matrix with decoupled velocity blocks
  integer, parameter, public :: CCMASM_MTP_DECOUPLED         = 1
  
  ! Extended 'full-tensor' matrix with submatrices A11, A12, A21, A22, all independent from
  ! each other.
  integer, parameter, public :: CCMASM_MTP_FULLTENSOR        = 2

!</constantblock>

!<constantblock description="Matrix types"

  ! Optimal control matrix (6 equations)
  integer, parameter, public :: MATT_OPTCONTROL = 0
  
  ! Linearised optimal control matrix (6 equations, Newton type)
  integer, parameter, public :: MATT_LINOPTCONTROL = 1

  ! Primal equation matrix (3 equations)
  integer, parameter, public :: MATT_PRIMAL = 2
  
  ! Linearised primal equation matrix (3 equations, Newton type)
  integer, parameter, public :: MATT_LINPRIMAL = 3

  ! Dual equation matrix (3 equations)
  integer, parameter, public :: MATT_DUAL = 4
  
  ! Linearised dual equation matrix (3 equations, Newton type)
  integer, parameter, public :: MATT_LINDUAL = 5


!</constantblock>

!</constants>

  public :: t_nonlinearSpatialMatrix
  public :: t_spatialMatrixNonlinearData
  public :: t_matrixAssemblyFlags
  public :: smva_initNonlinMatrix
  public :: smva_assembleMatrix
  public :: smva_assembleDefect
  public :: smva_getDiscrData
  public :: smva_initNonlinearData
  public :: smva_addBdEOJvector
  public :: smva_addBdEOJOperator
  public :: smva_disableSubmatrix
  public :: smva_clearMatrix
  public :: smva_prepareViscoAssembly
  public :: ffunctionViscoModel

!<types>

!<typeblock>
  ! Structure that encapsules all nonlinearities that may appear in the nonlinear
  ! spatial matrix.
  type t_spatialMatrixNonlinearData
  
    ! Specifies where to evaluate the nonlinearity and must contain the data
    ! for the 'previous' timestep. If there is no previous timestep (e.g.
    ! like in the 0th timestep), the vector can be undefined.
    type(t_vectorBlock), pointer :: p_rvector1 => null()

    ! Specifies where to evaluate the nonlinearity and must contain the data
    ! for the 'current' timestep.
    type(t_vectorBlock), pointer :: p_rvector2 => null()

    ! Specifies where to evaluate the nonlinearity and must contain the data
    ! for the 'next' timestep. If there is no next timestep (e.g.
    ! like in the last timestep), the vector can be undefined.
    type(t_vectorBlock), pointer :: p_rvector3 => null()
    
    ! Point in time associated to this matrix.
    real(DP) :: dassociatedTimePrimalVel = 0.0_DP
    real(DP) :: dassociatedTimeDualVel = 0.0_DP

    ! Definition of the Neumann boundary conditions.
    ! WARNING: IMPLEMENTATION INCOMPLETE!!!!
    ! ACTUALLY, THERE HAS TO BE ONE NEUMANN BOUDNARY STRUCTURE FOR EVERY NONLINEAR
    ! EVALUATION VECTOR FROM ABOVE. THE ASSEMBLY ROUTINE FOR THE BOUNDARY INTEGRAL
    ! IN THE DUAL EQUATION MUST ACTUALLY BE ABLE TO INCORPORATE THE NEUMANN
    ! BOUNDARY DEFINITION ACCORDING TO THE TIME OF THE NONLINEARITY; REPRESENTED
    ! BY THE ABOVE VECTOR(S). AS LONG AS THE BC'S ARE CONSTANT IN TIME, THE NEUMENN
    ! BC'S ARE THE SAME FOR ALL TIMESTEPS AND THIS ONE AND ONLY STRUCTURE IS ENOGH!!!
    type(t_boundaryRegionList), pointer :: p_rneumannBoundary => null()

    ! Definition of the Dirichlet control boundary conditions.
    type(t_boundaryRegionList), pointer :: p_rdirichletBCCBoundary => null()

  end type
  
!</typeblock>

!<typeblock>

  ! All discretisation and level related data structures that are necessary
  ! to evaluate the matrix or create a defect.
  type t_spatialMatrixDiscrData
  
    ! The physics of the problem
    type(t_settings_physics) :: rphysicsPrimal
    
    ! Discretisation settings of the problem (element, cubature rule,...)
    type(t_settings_discr) :: rsettingsSpaceDiscr

    ! Stabilisation parameters for the primal and dual system.
    type(t_settings_stabil) :: rstabilPrimal
    type(t_settings_stabil) :: rstabilDual
    
    ! Structure defining constraints of the problem.
    type(t_optcConstraintsSpace) :: rconstraints

    ! Pointer to the observation area or null(), if the complete
    ! area is to be observed.
    real(DP), dimension(:), pointer :: p_DobservationArea => null()
    
    ! Discretisation of the level, primal space.
    type(t_blockDiscretisation), pointer :: p_rdiscrPrimal => null()

    ! Discretisation of the level, primal/dual space.
    type(t_blockDiscretisation), pointer :: p_rdiscrPrimalDual => null()
    
    ! Assembly template data on that level.
    type(t_staticSpaceAsmTemplates), pointer :: p_rstaticAsmTemplates => null()

    ! Assembly template data on that level for the optimal control problem.
    type(t_staticSpaceAsmTemplatesOptC), pointer :: p_rstaticAsmTemplatesOptC => null()

    ! Pointer to program wide debug flags.
    type(t_optcDebugFlags), pointer :: p_rdebugFlags => null()

  end type

!</typeblock>

  public :: t_spatialMatrixDiscrData

!<typeblock>

!<!--
! The following structure encapsules all coefficients in front of the terms in each
! submatrix of the global matrix. Here, we don't address each submatrix in the
! global matrix but provide coefficients always for a complete block of all dimensions
! for the primal as well as for the dual term.
!
! Example: The Dmass, DBmat and DBTmat coefficients are used in the following way:
!
!    ( Dmass(1,1)*M   Dmass(1,1)*M   DBmat(1)*B1   Dmass(1,2)*M   Dmass(1,2)*M               )
!    ( Dmass(1,1)*M   Dmass(1,1)*M   DBmat(1)*B2   Dmass(1,2)*M   Dmass(1,2)*M               )
!    ( DBTmat(1)*B1^t DBTmat(1)*B2^t                                                         )
!    ( Dmass(2,1)*M   Dmass(2,1)*M                 Dmass(2,2)*M   Dmass(2,2)*M   DBmat(2)*B1 )
!    ( Dmass(2,1)*M   Dmass(2,1)*M                 Dmass(2,2)*M   Dmass(2,2)*M   DBmat(2)*B2 )
!    (                                             DBTmat(2)*B1^t DBTmat(2)*B2^t             )
!
! The use of the other parameters is similar.
! -->

  ! This structure describes a nonlinear matrix in space.
  ! The nonlinearity usually stems from a dependence of the Navier-Stokes operator
  ! on a solution (Oseen system) and stabilisation operators.
  type t_nonlinearSpatialMatrix
  
    ! IOTA-parameters that switch the identity on/off.
    real(DP), dimension(2,2) :: DidentityY = 0.0_DP
  
    ! ALPHA-parameters that switch the mass matrix on/off.
    real(DP), dimension(2,2) :: Dmass = 0.0_DP
  
    ! THETA-parameters that switch the Stokes matrix on/off
    real(DP), dimension(2,2) :: Dstokes = 0.0_DP
    
    ! GAMMA-parameters that switch the nonlinearity (y*grad(.)) on/off
    real(DP), dimension(2,2) :: Dygrad = 0.0_DP

    ! GAMMAADJ-parameters that switch the adjoint nonlinearity (.,y*grad(.)) on/off
    real(DP), dimension(2,2) :: DygradAdj = 0.0_DP
  
    ! NEWTON-parameters that switch the Newton term ((.)*grad(y)) on/off
    real(DP), dimension(2,2) :: Dgrady = 0.0_DP

    ! NEWTON-parameters that switch a 2nd Newton term ((.)*grad(y)) on/off
    real(DP), dimension(2,2) :: Dgrady2 = 0.0_DP

    ! NEWTON-parameters that switch the adjoint Newton term (.,(.*grad)(y)) on/off
    real(DP), dimension(2,2) :: DgradyAdj = 0.0_DP

    ! NEWTON-parameters that switch a 2nd adjoint Newton term (.,(.*grad)(y)) on/off
    real(DP), dimension(2,2) :: DgradyAdj2 = 0.0_DP

    ! GAMMAT-parameters that switch the transposed nonlinearity (y*grad(.)^T) on/off
    real(DP), dimension(2,2) :: DygradT = 0.0_DP

    ! GAMMATADJ-parameters that switch the adjoint transposed nonlinearity (.,y*grad(.)^T) on/off
    real(DP), dimension(2,2) :: DygradTAdj = 0.0_DP

    ! GAMMAT-parameters that switch a 2nd transposed nonlinearity (y*grad(.)^T) on/off
    real(DP), dimension(2,2) :: DygradT2 = 0.0_DP

    ! GAMMATADJ-parameters that switch a 2nd transposed nonlinearity (.,y*grad(.)^T) on/off
    real(DP), dimension(2,2) :: DygradTAdj2 = 0.0_DP
  
    ! NEWTONT-parameters that switch the transposed Newton term ((.)*grad(y)^T) on/off
    real(DP), dimension(2,2) :: DgradyT = 0.0_DP
    
    ! ETA-parameters that switch the B-terms on/off.
    real(DP), dimension(2,2) :: DBmat = 0.0_DP
    
    ! TAU-parameters that switch the B^T-terms on/off
    real(DP), dimension(2,2) :: DBTmat = 0.0_DP
    
    ! KAPPA-parameters that switch the I matrix in the continuity equation
    ! on/off.
    real(DP), dimension(2,2) :: DidentityP = 0.0_DP
    
    ! Weight in front of Y in the Dirichlet boundary control term 1.
    real(DP), dimension(2,2) :: DdirichletBCCY = 0.0_DP
    
    ! Weight in front of Dlambda in the Dirichlet boundary control term 2.
    real(DP), dimension(2,2) :: DdirichletBCCLambda = 0.0_DP
    
    ! Weight in front of xi in the Dirichlet boundary control term 3.
    real(DP), dimension(2,2) :: DdirichletBCCXi = 0.0_DP
    
    ! When evaluating nonlinear terms, the evaluation routine accepts
    ! in timestep i three solution vectors -- one corresponding to
    ! the matrix Aii-1, one corresponding to matrix Aii and one corresponding
    ! to matrix Aii+1. The following id's specify which primal solution is
    ! used for which matrix. A value of 1 correspponds to solution i-1,
    ! a value of 2 to solution i and a value of 3 to solution i+1.
    ! iprimalSol=1 means e.g. that A=A(xi-1). Similarly,
    ! iprimalSol=3 would mean that A=A(xi+1).
    ! =0: undefined.
    integer :: iprimalSol  = 2

    ! The following three variables specify the id of the dual solution to be
    ! evaluated. 1=lambda_{i-1}, 2=lambda_i, 3=lambda_{i+1}.
    ! =0: undefined.
    integer :: idualSol  = 2

    ! The following three variables specify the id of a 2nd dual solution to be
    ! evaluated. 1=lambda_{i-1}, 2=lambda_i, 3=lambda_{i+1}.
    ! =0: undefined.
    ! This is used for the discretisation of the 2nd Newton/GammT-term which
    ! appears on the main diagonal of the system for the dual equation
    ! if Crank-Nicolson is active.
    integer :: idualSol2  = 2

    ! Regularisation parameter ALPHA which controls the transformation
    ! from the dual variable $\lambda$ to the control $u$:
    ! $u=-1/\alpha \lambda$.
    real(DP) :: dalphaC = 0.0_DP
    
    ! Regularisation parameter ALPHA which controls the transformation
    ! from the dual variable $\lambda$ to the control $u$ on the boundary:
    ! $u=-1/\beta \lambda$.
    real(DP) :: dbetaC = 0.0_DP
    
    ! Penalty parameter for the dirichlet boundary control
    real(DP) :: ddirichletBCPenalty = 0.0_DP
    
    ! Type of this matrix. One of the MATT_xxxx constants.
    integer :: cmatrixType = 0

    ! Discretisation related data.
    type(t_spatialMatrixDiscrData) :: rdiscrData
    
    ! Pointer to global data.
    type(t_globalData), pointer :: p_rglobalData => null()

    ! Pointer to a structure specifying all nonlinearities.
    ! May point to NULL if there is no nonlinearity.
    type(t_spatialMatrixNonlinearData), pointer :: p_rnonlinearity => null()

  end type

!</typeblock>

!<typeblock>
  ! Additional status flags that configure how to assemble a matrix.
  type t_matrixAssemblyFlags
  
    ! MATRIX RESTRICTION: Parameter to activate matrix restriction.
    ! Can be used to generate parts of the matrices on coarse grids where the
    ! aspect ratio of the cells is large. Only applicable for $\tilde Q_1$
    ! discretisations.
    ! Standard = 0 = deactivate matrix restriction
    integer :: iadaptiveMatrices = 0
    
    ! MATRIX RESTRICTION: Threshold parameter for adaptive matrix generation
    ! of coarse grid matrices (here: maximum aspect ratio).
    ! Only applicable if iadaptiveMatrices <> 0.
    ! Standard = 20.0
    real(DP) :: dadmatthreshold = 20.0_DP
    
    ! ASSEMBLY SPECIALS: If set to TRUE, the D-matrices are assembled as
    ! virtually transposed B-matrices. With the standard setting FALSE, the
    ! D-matrices are assembled as standard divergence matrices.
    ! WARNING: If matrices are to be assembled as virtually transposed
    ! matrices, they must not be overwritten, as they share their data with
    ! other matrices.
    ! If the matrix is assembled in the standard way without virtual
    ! transposition, separate memory is requested for the matrices and so
    ! they are allowed to be overwritten if necessary. (Neessary for UMFPACK
    ! on a pure Dirichlet problem.)
    logical :: bvirtualTransposedD = .false.

  end type
!</typeblock>

!</types>

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine smva_getDiscrData (rsettings, ilevel, rdiscrData,rphysics)
  
!<description>
  ! Fetches all discretisation data from the main program structure that is
  ! necessary to set up nonlinear matrices.
!</description>

!<input>
  ! The structure of the main solver
  type(t_settings_optflow), intent(in), target :: rsettings
  
  ! Level where the discretisation takes place.
  integer, intent(in) :: ilevel

  ! OPTIONAL: Alternative physics definition to use.
  ! If not present, the standard global physics settings are used.
  type(t_settings_physics), intent(in), optional :: rphysics
!</input>

!<output>
  ! Discretisation related data, can be used to generate nonlinear matrices.
  type(t_spatialMatrixDiscrData), intent(out) :: rdiscrData
!</output>

!</subroutine>

    ! Get level-independent data.
    rdiscrData%rphysicsPrimal = rsettings%rphysicsPrimal
    if (present(rphysics)) then
      rdiscrData%rphysicsPrimal = rphysics
    end if
    rdiscrData%rsettingsSpaceDiscr = rsettings%rsettingsSpaceDiscr
    rdiscrData%rstabilPrimal = rsettings%rstabilPrimal
    rdiscrData%rstabilDual = rsettings%rstabilDual
    rdiscrData%rconstraints = rdiscrData%rconstraints

    rdiscrData%p_rdiscrPrimal => &
        rsettings%rfeHierPrimal%p_rfeSpaces(ilevel)%p_rdiscretisation

    rdiscrData%p_rdiscrPrimalDual => &
        rsettings%rfeHierPrimalDual%p_rfeSpaces(ilevel)%p_rdiscretisation
    
    rdiscrData%p_rstaticAsmTemplates => &
        rsettings%rspaceAsmHierarchy%p_RasmTemplList(ilevel)

    rdiscrData%p_rstaticAsmTemplatesOptC => &
        rsettings%rspaceAsmHierarchyOptC%p_RasmTemplList(ilevel)

    rdiscrData%p_rdebugFlags => rsettings%rdebugFlags
    
    rdiscrData%p_DobservationArea => rsettings%rsettingsOptControl%p_DobservationArea

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine smva_initNonlinMatrix (rnonlinearSpatialMatrix,rglobalData,rdiscrData,rnonlinearity)

!<description>
  ! Initialises the rnonlinearCCMatrix structure.
!</description>

!<input>
  ! Global program data
  type(t_globalData), target :: rglobalData

  ! Discretisation related data with information about the level, the
  ! matrix acts on.
  type(t_spatialMatrixDiscrData) :: rdiscrData
  
  ! Structure specifying the nonlinearity evaluation points of the matrix.
  type(t_spatialMatrixNonlinearData), target :: rnonlinearity
!</input>

!<inputoutput>
  ! Nonlinear matrix structure.
  ! Basic parameters in this structure are filled with data.
  type(t_nonlinearSpatialMatrix), intent(inout) :: rnonlinearSpatialMatrix
!</inputoutput>
              
!</subroutine>
    
    ! Initialise the matrix assembly structure
    ! with basic global information.
    rnonlinearSpatialMatrix%rdiscrData = rdiscrData

    ! Change the sign of dupsam2 for a consistent stabilisation.
    ! Reason: The stablisation is added to the dual operator by the SD/
    ! EOJ stabilisation in the following way:
    !
    !    ... - (u grad lamda + dupsam2*stabilisation) + ... = rhs
    !
    ! We want to *add* the stabilisation, so we have to introduce a "-" sign
    ! in dupsam2 to get
    !
    !    ... - (u grad lamda) - (-dupsam2*stabilisation) + ... = rhs
    ! <=>
    !    ... - (u grad lamda) + dupsam2*stabilisation + ... = rhs
    
    rnonlinearSpatialMatrix%rdiscrData%rstabilDual%dupsam = &
        -rnonlinearSpatialMatrix%rdiscrData%rstabilDual%dupsam * &
        mprim_signum(rnonlinearSpatialMatrix%rdiscrData%p_rdebugFlags%dweightDualConvection)
        
    ! Remember the evaluation point of the nonlinearity.
    rnonlinearSpatialMatrix%p_rnonlinearity => rnonlinearity
    
    ! Global program data
    rnonlinearSpatialMatrix%p_rglobalData => rglobalData
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine smva_nonlinMatrixSetNonlin (rnonlinearSpatialMatrix,rnonlinearity)

!<description>
  ! Defines the nonlinearity of a nonlinear matrix. Allows to change the
  ! nonlinerity in rnonlinearSpatialMatrix to rnonlinearity.
!</description>

!<input>
  ! Structure specifying the nonlinearity evaluation points of the matrix.
  type(t_spatialMatrixNonlinearData), target :: rnonlinearity
!</input>

!<inputoutput>
  ! Nonlinear matrix structure.
  ! Basic parameters in this structure are filled with data.
  type(t_nonlinearSpatialMatrix), intent(inout) :: rnonlinearSpatialMatrix
!</inputoutput>
              
!</subroutine>
    
    ! Remember the evaluation point of the nonlinearity.
    rnonlinearSpatialMatrix%p_rnonlinearity => rnonlinearity
    
  end subroutine

! ***************************************************************************
  !<subroutine>

  subroutine coeff_ProjMass (rdiscretisationTrial,rdiscretisationTest,rform, &
      nelements,npointsPerElement,Dpoints, IdofsTrial,IdofsTest,rdomainIntSubset, &
      Dcoefficients,rcollection)
    
  !<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form
    ! that assembles the projective mass matrix.
    !
    ! The coefficients is c=c(lambda_1) or =c(lambda_2) with
    ! c=1/alpha if a < -1/alpha lambda_i < b and c=0 otherwise. This is the derivative
    ! of the projection operator "-P[a,b](-1/alpha lambda_i)" on the left hand
    ! side of the equation.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(IN)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                        :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in trial space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    ! local variables
    type(t_vectorBlock), pointer :: p_rvector
    type(t_vectorScalar), pointer :: p_rsubvector
    real(dp), dimension(:,:), allocatable :: Dfunc
    integer(I32) :: celement
    real(DP) :: da, db, dalpha, dp1
    integer :: ipt, iel
    
    ! Get the bounds and the multiplier from the collection
    dalpha = rcollection%DquickAccess(3)
    dp1 = rcollection%DquickAccess(4)
    
    ! Scale the bounds by -alpha as we analyse lambda
    ! and not u. Change the role of min/max because of the "-"
    ! sign!
    da = -rcollection%DquickAccess(2)*dalpha
    db = -rcollection%DquickAccess(1)*dalpha
    
    ! Get a pointer to the FE solution from the collection.
    ! The routine below wrote a pointer to the vector T to the
    ! first quick-access vector pointer in the collection.
    p_rvector => rcollection%p_rvectorQuickAccess1

    ! Do we have to analyse lambda_1 or lambda_2?
    if (rcollection%IquickAccess(1) .eq. 1) then
      p_rsubvector => p_Rvector%RvectorBlock(4)
    else
      p_rsubvector => p_Rvector%RvectorBlock(5)
    end if
  
    ! Allocate memory for the function values in the cubature points:
    allocate(Dfunc(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))
    
    ! Calculate the function value of the solution vector in all
    ! our cubature points:
    !
    ! Figure out the element type, then call the
    ! evaluation routine for a prepared element set.
    ! This works only if the trial space of the matrix coincides
    ! with the FE space of the vector T we evaluate!
    
    celement = rdomainIntSubset%celement
    
    call fevl_evaluate_sim (p_rsubvector, &
        rdomainIntSubset%p_revalElementSet, &
        celement, rdomainIntSubset%p_IdofsTrial, DER_FUNC, Dfunc)
    
    ! Now check the function values lambda.
    ! If b < -1/alpha lambda < a, return 1/alpha.
    ! Otherwise, return 0.
    do iel = 1,ubound(Dcoefficients,3)
      do ipt = 1,ubound(Dcoefficients,2)
        ! Check if the dual variable is in the bounds for the control.
        if ((Dfunc(ipt,iel) .gt. da) .and. (Dfunc(ipt,iel) .lt. db)) then
          Dcoefficients(1,ipt,iel) = dp1
        else
          Dcoefficients(1,ipt,iel) = 0.0_DP
        end if
      end do
    end do
    
    ! Release memory
    deallocate(Dfunc)

  end subroutine

! ***************************************************************************
  !<subroutine>

  subroutine coeff_ProjMassCollect (rdiscretisationTrial,rdiscretisationTest,rform, &
      nelements,npointsPerElement,Dpoints, IdofsTrial,IdofsTest,rdomainIntSubset, &
      Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form
    ! that assembles the projective mass matrix.
    !
    ! The coefficients is c=c(lambda_1) or =c(lambda_2) with
    ! c=1/alpha if a < -1/alpha lambda_i < b and c=0 otherwise. This is the derivative
    ! of the projection operator "-P[a,b](-1/alpha lambda_i)" on the left hand
    ! side of the equation.
    !
    ! Collects the elements of the active set for later re-assembly.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(IN)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                        :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in trial space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    ! local variables
    type(t_vectorBlock), pointer :: p_rvector
    type(t_vectorScalar), pointer :: p_rsubvector
    real(dp), dimension(:,:), allocatable :: Dfunc
    integer(I32) :: celement
    real(DP) :: da, db, dalphaC, dp1
    integer :: ipt, iel
    integer :: nptsInactive
    integer, dimension(:), pointer :: p_IelementList
    
    ! Get the bounds and the multiplier from the collection
    dalphaC = rcollection%DquickAccess(3)
    dp1 = rcollection%DquickAccess(4)
    
    ! Forget it if alpha<=0. Not possible.
    if (dalphaC .le. 0.0_DP) return
    
    ! Scale the bounds by -alpha as we analyse lambda
    ! and not u. Change the role of min/max because of the "-"
    ! sign!
    da = -rcollection%DquickAccess(2)*dalphaC
    db = -rcollection%DquickAccess(1)*dalphaC
    
    ! Get a pointer to the FE solution from the collection.
    ! The routine below wrote a pointer to the vector T to the
    ! first quick-access vector pointer in the collection.
    p_rvector => rcollection%p_rvectorQuickAccess1

    ! Do we have to analyse lambda_1 or lambda_2?
    if (rcollection%IquickAccess(1) .eq. 1) then
      p_rsubvector => p_Rvector%RvectorBlock(4)
    else
      p_rsubvector => p_Rvector%RvectorBlock(5)
    end if
  
    ! Allocate memory for the function values in the cubature points:
    allocate(Dfunc(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))
    
    ! Calculate the function value of the solution vector in all
    ! our cubature points:
    !
    ! Figure out the element type, then call the
    ! evaluation routine for a prepared element set.
    ! This works only if the trial space of the matrix coincides
    ! with the FE space of the vector T we evaluate!
    
    celement = rdomainIntSubset%celement
    
    call fevl_evaluate_sim (p_rsubvector, &
        rdomainIntSubset%p_revalElementSet, &
        celement, rdomainIntSubset%p_IdofsTrial, DER_FUNC, Dfunc)
    
    call storage_getbase_int (rcollection%IquickAccess(2),p_IelementList)
    
    ! Now check the function values lambda.
    ! If b < -1/alpha lambda < a, return 1/alpha.
    ! Otherwise, return 0.
    !
    ! If we detect an element, which is on the border of the active set,
    ! return 0 everywhere and collect the element to the element list,
    ! so the assembly routine can assemble there on a 2nd pass later.
    do iel = 1,ubound(Dcoefficients,3)
      
      ! Count the inactive points.
      nptsInactive = 0
      do ipt = 1,ubound(Dcoefficients,2)
        ! Check if the dual variable is in the bounds for the control.
        if ((Dfunc(ipt,iel) .gt. da) .and. (Dfunc(ipt,iel) .lt. db)) then
          nptsInactive = nptsInactive + 1
        end if
      end do
      
      ! All points inactive? Ok.
      ! Partially active? Remember the element.
      ! Completely active? Return 0 everywhere.
      if (nptsInactive .eq.  ubound(Dcoefficients,2)) then
        do ipt = 1,ubound(Dcoefficients,2)
          Dcoefficients(1,ipt,iel) = dp1
        end do
      else
        do ipt = 1,ubound(Dcoefficients,2)
          Dcoefficients(1,ipt,iel) = 0.0_DP
        end do
        if (nptsInactive .gt. 0) then
          rcollection%IquickAccess(3) = rcollection%IquickAccess(3) + 1
          p_IelementList(rcollection%IquickAccess(3)) = &
              rdomainIntSubset%p_Ielements(iel)
        end if
      end if
    end do
    
    ! Release memory
    deallocate(Dfunc)

  end subroutine

! ***************************************************************************
  !<subroutine>

  subroutine coeff_boxIdentity (rdiscretisationTrial,rdiscretisationTest,rform, &
      nelements,npointsPerElement,Dpoints, IdofsTrial,IdofsTest,rdomainIntSubset, &
      Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form
    ! that assembles the projective mass matrix.
    !
    ! Assembles the box-constrained identity, i.e., returns
    ! =1 for all points in a defined box coming from the collection.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(IN)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                        :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in trial space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    integer :: iel,ipt
    real(DP) :: dx1, dy1, dx2, dy2, dweight
    
    ! Get the box
    dx1 = rcollection%DquickAccess(1)
    dy1 = rcollection%DquickAccess(2)
    dx2 = rcollection%DquickAccess(3)
    dy2 = rcollection%DquickAccess(4)
    dweight = rcollection%DquickAccess(5)

    ! Loop through the elements and cubature points.
    ! Set everything to zero which is out of our bounds.
    do iel = 1,nelements
      do ipt = 1,npointsPerElement
        if ( (Dpoints(1,ipt,iel) .lt. dx1) .or. (Dpoints(1,ipt,iel) .gt. dx2) .or. &
              (Dpoints(2,ipt,iel) .lt. dy1) .or. (Dpoints(2,ipt,iel) .gt. dy2) ) then
          Dcoefficients (1,ipt,iel) = 0.0_DP
        else
          Dcoefficients (1,ipt,iel) = dweight
        end if
      end do
    end do

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine assembleBoxIdentity (Dbox,dweight,rmatrix,rcubatureInfo)

!<description>
  ! Assembles the characteristic function of a box.
!</description>  
  
!<input>
  ! The box. Format: (x1, y1, x2, y2)
  real(DP), dimension(:), intent(in) :: Dbox
  
  ! Weight for the operator
  real(DP), intent(in) :: dweight
  
  ! Cubature information structure
  type(t_scalarCubatureInfo), intent(in) :: rcubatureInfo
!</input>

!<inputoutput>
  ! Matrix which receives the operator.
  type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>
  
!</subroutine>

    ! local variables
    type(t_bilinearForm) :: rform
    type(t_collection) :: rcollection
    
    ! Prepare a bilinear form.
    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_FUNC

    ! In this case, we have nonconstant coefficients.
    rform%ballCoeffConstant = .false.
    rform%BconstantCoeff(:) = .false.

    ! Prepare a collection structure with the box.
    rcollection%DquickAccess(1) = Dbox(1)
    rcollection%DquickAccess(2) = Dbox(2)
    rcollection%DquickAccess(3) = Dbox(3)
    rcollection%DquickAccess(4) = Dbox(4)
    rcollection%DquickAccess(5) = dweight

    ! Now we can build the matrix entries.
    call bilf_buildMatrixScalar (rform,.true.,rmatrix,&
        rcubatureInfo,coeff_boxIdentity,rcollection)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine smva_assembleMatrix (coperation,cmatrixType,rflags,&
      rnonlinearSpatialMatrix,rmatrix,rfineMatrix)

!<description>
  ! This routine assembles a matrix space. The caller must initialise
  ! the rnonlinearSpatialMatrix according to how the matrix should look like.
  ! The 'coperation' parameter tells the routine what to do.
  ! The destination matrix rmatrix is then set up or updated.
  !
  ! The parameters rvector and rfineMatrix are optional. rvector must be
  ! specified, if the nonlinearity is activated (parameter $\gamma\not=0$ in
  ! rnonlinearSpatialMatrix). This vector specifies the 'solution' where the nonlinearity
  ! $u\nabla u$ is evaluated.
  ! rfineMatrix allows to specify a matrix of a 'one level refined mesh'. This
  ! is usually used when setting up preconditioners over multiple levels.
  ! Specifying such a matrix allows the routine (depending on the discretisation)
  ! to include some special stabilisation terms into the matrix rmatrix.
  !
  ! The routine will not include any boundary conditions in the matrix.
!</description>

!<input>

  ! One of the CCMASM_xxxx-constants. This parameter specifies 'what to do'.
  ! Using the CCMASM_ALLOCMEM constant, the routine will allocate memory in
  ! rmatrix for all the matrix blocks but will not compute any entries.
  ! Using the CCMASM_COMPUTE constant, the routine will compute the matrix
  ! entries, while assuming that the memory was allocated previously.
  ! Using the CCMASM_ALLOCANDCOMPUTE constant, the routine will do both,
  ! allocate memory and compute the entries.
  !
  ! If any of the CCMASM_ALLOCxxxx constants is specified here and rmatrix is
  ! already initialised, rmatrix is released and completely rebuild in memory!
  ! Therefore, CCMASM_COMPUTE should be used to update an existing matrix.
  !
  ! The constant CMASM_QUICKREFERENCES may be specified additional to one of
  ! the other constants (e.g. as 'CCMASM_ALLOCANDCOMPUTE+CMASM_QUICKREFERENCES').
  ! If specified, the routine tries to avoid memory allocation. This means e.g.
  ! that references to the original gradient (B-)matrices from rnonlinearSpatialMatrix
  ! are written to rmatrix; matrix entries are not copied!
  ! (This can be used e.g. for setting up a matrix for building a defect
  !  vector without copying matrix data.)
  ! In this case, the caller MUST NOT CHANGE rmatrix in any way, otherwise
  ! the original (template) matrices would be changed!
  integer(I32), intent(IN) :: coperation

  ! Type of matrix that should be set up in rmatrix. One of the CCMASM_MTP_xxxx
  ! constants.
  ! Usually, CCMASM_MTP_AUTOMATIC is used here. This will automatically determine
  ! the 'smallest possible' matrix structure fitting the needs of the input
  ! parameters. By specifying another matrix type, the caller can explicitly
  ! take influence on the general matrix structure.
  !
  ! If the matrix already exists and only its entries are to be computed,
  ! CCMASM_MTP_AUTOMATIC should be specified here.
  integer, intent(IN) :: cmatrixType
  
  ! Additional matrix assembly flags that configure how to generate the matrix.
  type(t_matrixAssemblyFlags), intent(in) :: rflags

  ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
  ! about how to set up the matrix.
  !
  ! Note that if coperation=CCMASM_ALLOCxxxx is specified, p_rmatrixTemplateXXXX
  ! must be initialised as well as p_rdiscretisation!
  ! The new matrix is created based p_rmatrixTemplateXXXX as well as p_rdiscretisation.
  ! Memory is automatically allocated if it's missing.
  type(t_nonlinearSpatialMatrix), intent(in) :: rnonlinearSpatialMatrix

  ! OPTIONAL: This parameter allows to specify a 'fine grid matrix'. This is
  ! usually done when assembling matrices on multiple levels. If specified, the
  ! routine will (if possible) try to include a level-dependent stabilisation
  ! term into the matrix (-> e.g. constant matrix restriction for nonparametric
  ! Rannacher-Turek element if cells are too anisotropic).
  type(t_matrixBlock), intent(IN), optional :: rfineMatrix
  
!</input>

!<inputoutput>

  ! The destination matrix which should be set up.
  ! If not initialised, a new matrix is created (as if CCMASM_ALLOCxxxx
  ! was specified).
  ! If initialised, the existing matrix is updated or recreated, depending on
  ! coperation.
  type(t_matrixBlock), intent(inout) :: rmatrix
  
!</inputoutput>
  
!</subroutine>

    ! local variables
    logical :: ballocate
    type(t_matrixBlock) :: rtempMatrix
    type(t_vectorBlock) :: rvectorPrimal, rvectorDual, rvectorDual2
    type(t_settings_stabil) :: rstabilisation
    integer :: i,j
    type(t_optcoperator) :: roptcoperator
    type(t_vectorBlock), pointer :: p_rprimalSol, p_rdualSol
    type(t_blockDiscretisation) :: rvelDiscr
    real(dp), dimension(:), pointer :: p_Ddata
    real(DP) :: dweightConvection, dweightDualConvection, dweightNaturalBdcDual
    real(DP) :: dweightDualNewtonT
    type(t_scalarCubatureInfo) :: rcubatureInfo
    
    logical, parameter :: bnewmethod = .false.
    
    ! Debug weights for the convection
    dweightConvection = &
        rnonlinearSpatialMatrix%rdiscrData%p_rdebugFlags%dweightConvection
    dweightDualConvection = &
        rnonlinearSpatialMatrix%rdiscrData%p_rdebugFlags%dweightDualConvection
    dweightNaturalBdcDual = &
        rnonlinearSpatialMatrix%rdiscrData%p_rdebugFlags%dweightNaturalBdcDual
    dweightDualNewtonT = &
        rnonlinearSpatialMatrix%rdiscrData%p_rdebugFlags%dweightDualNewtonT

    ballocate = .false.
    if ((rmatrix%NEQ .le. 0) .or. &
        iand(coperation,CCMASM_ALLOCMEM) .ne. 0) then
      ballocate = .true.
    end if
    
    ! What should we do? Allocate memory?
    if (ballocate) then
    
      ! Release the matrix if present.
      call lsysbl_releaseMatrix (rmatrix)
    
      ! Allocate the system matrix
      call lsysbl_createMatBlockByDiscr(&
          rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimalDual,rmatrix)
      do j=1,2
        do i=1,2
          call lsysbl_createEmptyMatrix (rtempMatrix,NDIM2D+1,NDIM2D+1)
          
          call allocSubmatrix (rnonlinearSpatialMatrix,cmatrixType,rflags,rtempMatrix,&
              rnonlinearSpatialMatrix%DidentityY(i,j),rnonlinearSpatialMatrix%Dmass(i,j),&
              rnonlinearSpatialMatrix%Dstokes(i,j),&
              rnonlinearSpatialMatrix%Dygrad(i,j),rnonlinearSpatialMatrix%DygradAdj(i,j),&
              rnonlinearSpatialMatrix%DygradT(i,j),rnonlinearSpatialMatrix%DygradTAdj(i,j),&
              rnonlinearSpatialMatrix%Dgrady(i,j),rnonlinearSpatialMatrix%DgradyT(i,j),&
              rnonlinearSpatialMatrix%DgradyAdj(i,j),&
              rnonlinearSpatialMatrix%DBmat(i,j),rnonlinearSpatialMatrix%DBTMat(i,j),&
              rnonlinearSpatialMatrix%DidentityP(i,j),&
              rnonlinearSpatialMatrix%DdirichletBCCY(i,j),&
              rnonlinearSpatialMatrix%DdirichletBCCLambda(i,j),&
              rnonlinearSpatialMatrix%DdirichletBCCXi(i,j))
          call lsysbl_updateMatStrucInfo(rtempMatrix)
          call lsysbl_moveToSubmatrix (rtempMatrix,rmatrix,(i-1)*3+1,(j-1)*3+1)
          call lsysbl_releaseMatrix (rtempMatrix)
        end do
      end do
      
      call lsysbl_allocEmptyMatrix (rmatrix,LSYSSC_SETM_UNDEFINED,.true.)

    end if
   
    if (iand(coperation,CCMASM_COMPUTE) .ne. 0) then
    
      ! Create a cubature information structure that defines how to
      ! set up the nonlinear operator.
      call spdiscr_createDefCubStructure (&
          rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimalDual%RspatialDiscr(1),&
          rcubatureInfo,rnonlinearSpatialMatrix%rdiscrData%rsettingsSpaceDiscr%icubStokes)

!#ifndef TESTCODE
!      call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,3),p_Ddata)
!#endif

      if (.not. bnewmethod) then
        ! The system matrix in the whole looks like:
        !
        !    ( A11  A12  B1  R1           )
        !    ( A21  A22  B2       R1      )
        !    ( B1^T B2^T I                )
        !    ( R2            A44  A45  B1 )
        !    (      R2       A54  A55  B2 )
        !    (               B1^T B2^T I  )
        !
        ! With some multiplication factors in front of the matrices.
        !
        ! At first switch off all submatrices. Those submatrices that we
        ! assign some data are explicitely switched on again.
        rmatrix%RmatrixBlock(:,:)%dscaleFactor = 0.0_DP
        
        ! Derive a discretisation structure for either primal or
        ! dual variables. Create a temporary matrix based on that.
        call spdiscr_deriveBlockDiscr (&
            rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimalDual, rvelDiscr, 1, 3)
        call lsysbl_createMatBlockByDiscr (rvelDiscr,rtempMatrix)

        ! Get the primal velocity
        select case (rnonlinearSpatialMatrix%iprimalSol)
        case (1)
          call lsysbl_deriveSubvector(&
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1,rvectorPrimal, 1,2,.true.)
        case (2)
          call lsysbl_deriveSubvector(&
            rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2,rvectorPrimal, 1,2,.true.)
        case (3)
          call lsysbl_deriveSubvector(&
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3,rvectorPrimal, 1,2,.true.)
        end select

        ! Get the dual velocity.
        ! For pure forward/backward equations, this may be missing if not used.
        ! In this case, the temp vector stays uninitialised.
        select case (rnonlinearSpatialMatrix%idualSol)
        case (1)
          if (rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1%nblocks .gt. 3) then
            call lsysbl_deriveSubvector(&
                rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1,rvectorDual, 4,5,.true.)
          end if
        case (2)
          if (rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2%nblocks .gt. 3) then
            call lsysbl_deriveSubvector(&
                rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2,rvectorDual, 4,5,.true.)
          end if
        case (3)
          if (rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3%nblocks .gt. 3) then
            call lsysbl_deriveSubvector(&
                rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3,rvectorDual, 4,5,.true.)
          end if
        end select
        
        ! Get the 2nd dual velocity for reactive terms.
        ! For pure forward/backward equations, this may be missing if not used.
        ! In this case, the temp vector stays uninitialised.
        select case (rnonlinearSpatialMatrix%idualSol2)
        case (1)
          if (rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1%nblocks .gt. 3) then
            call lsysbl_deriveSubvector(&
                rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1,rvectorDual2, 4,5,.true.)
          end if
        case (2)
          if (rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2%nblocks .gt. 3) then
            call lsysbl_deriveSubvector(&
                rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2,rvectorDual2, 4,5,.true.)
          end if
        case (3)
          if (rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3%nblocks .gt. 3) then
            call lsysbl_deriveSubvector(&
                rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3,rvectorDual2, 4,5,.true.)
          end if
        end select
        
        ! Primal equation
        ! ---------------
        ! In the first step, assemble the linear parts
        ! (Laplace, Mass, B/B^T) on the main diagonal of the primal equation.
        call lsysbl_extractSubmatrix (rmatrix,rtempMatrix,1,3)
        call assembleLinearSubmatrix (rnonlinearSpatialMatrix,cmatrixType,rflags,&
            rtempMatrix,rnonlinearSpatialMatrix%DidentityY(1,1),rnonlinearSpatialMatrix%Dmass(1,1),&
            rnonlinearSpatialMatrix%Dstokes(1,1),&
            rnonlinearSpatialMatrix%DBmat(1,1),rnonlinearSpatialMatrix%DBTmat(1,1))


        ! Assemble the nonlinearity u*grad(.) or the Newton nonlinearity
        ! u*grad(.)+grad(u)*(.) to the velocity.
        
        call assembleNonlinearity (&
            rnonlinearSpatialMatrix,rtempMatrix,&
            rvectorPrimal,rvectorPrimal,&
            rnonlinearSpatialMatrix%Dstokes(1,1),&
            dweightConvection*rnonlinearSpatialMatrix%Dygrad(1,1),&
            dweightConvection*rnonlinearSpatialMatrix%DygradAdj(1,1),&
            dweightConvection*rnonlinearSpatialMatrix%DygradT(1,1),&
            dweightConvection*rnonlinearSpatialMatrix%DygradTAdj(1,1),&
            dweightConvection*rnonlinearSpatialMatrix%Dgrady(1,1),&
            dweightConvection*rnonlinearSpatialMatrix%DgradyT(1,1),&
            dweightConvection*rnonlinearSpatialMatrix%DgradyAdj(1,1),&
            rcubatureInfo,rnonlinearSpatialMatrix%rdiscrData%rstabilPrimal,&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplatesOptC%rmatrixEOJ1)

        ! Reintegrate the computed matrix
        call lsysbl_moveToSubmatrix (rtempMatrix,rmatrix,1,1)

        ! Include the mass matrix blocks
        !
        !    (  .    .    .    R1   .    . )
        !    (  .    .    .    .    R1   . )
        !    (  .    .    .    .    .    . )
        
        call lsysbl_extractSubmatrix (rmatrix,rtempMatrix,1,3,4,6)
        select case (rnonlinearSpatialMatrix%idualSol)
        case (1)
          call assembleProjectedMassBlocks (rnonlinearSpatialMatrix,rtempMatrix, &
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1, &
              rnonlinearSpatialMatrix%Dmass(1,2))
        case (2)
          call assembleProjectedMassBlocks (rnonlinearSpatialMatrix,rtempMatrix, &
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2, &
              rnonlinearSpatialMatrix%Dmass(1,2))
        case (3)
          call assembleProjectedMassBlocks (rnonlinearSpatialMatrix,rtempMatrix, &
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3,&
              rnonlinearSpatialMatrix%Dmass(1,2))
        end select

        ! Reintegrate the computed matrix
        call lsysbl_moveToSubmatrix (rtempMatrix,rmatrix,1,4)
        
        ! Dirichlet boudary control operator.
        call smva_assembleDirichletBCC (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal,&
            rmatrix%RmatrixBlock(1,1),rmatrix%RmatrixBlock(2,2),&
            rmatrix%RmatrixBlock(1,4),rmatrix%RmatrixBlock(2,5),&
            rmatrix%RmatrixBlock(1,6),rmatrix%RmatrixBlock(2,6),&
            rnonlinearSpatialMatrix%DdirichletBCCY(1,1),&
            rnonlinearSpatialMatrix%DdirichletBCCLambda(1,2),&
            rnonlinearSpatialMatrix%DdirichletBCCXi(1,2),&
            rnonlinearSpatialMatrix%ddirichletBCPenalty,&
            rnonlinearSpatialMatrix%p_rnonlinearity%p_rdirichletBCCBoundary)
        
        ! Dual equation
        ! -------------
        ! In the first step, assemble the linear parts
        ! (Laplace, Mass, B/B^T) on the main diagonal of the primal equation.
        call lsysbl_extractSubmatrix (rmatrix,rtempMatrix,4,6)
        call assembleLinearSubmatrix (rnonlinearSpatialMatrix,cmatrixType,rflags,rtempMatrix,&
            rnonlinearSpatialMatrix%DidentityY(2,2),rnonlinearSpatialMatrix%Dmass(2,2),&
            rnonlinearSpatialMatrix%Dstokes(2,2),&
            rnonlinearSpatialMatrix%DBmat(2,2),rnonlinearSpatialMatrix%DBTmat(2,2))

        ! Assemble the nonlinearity u*grad(.) or the Newton nonlinearity
        ! u*grad(.)+grad(u)*(.) to the velocity.
        call assembleNonlinearity (&
            rnonlinearSpatialMatrix,rtempMatrix,&
            rvectorPrimal,rvectorPrimal,&
            rnonlinearSpatialMatrix%Dstokes(2,2),&
            dweightDualConvection*rnonlinearSpatialMatrix%Dygrad(2,2),&
            dweightDualConvection*rnonlinearSpatialMatrix%DygradAdj(2,2),&
            dweightConvection*rnonlinearSpatialMatrix%DygradT(2,2),&
            dweightConvection*rnonlinearSpatialMatrix%DygradTAdj(2,2),&
            dweightConvection*rnonlinearSpatialMatrix%Dgrady(2,2),&
            dweightDualNewtonT*rnonlinearSpatialMatrix%DgradyT(2,2),&
            dweightConvection*rnonlinearSpatialMatrix%DgradyAdj(2,2),&
            rcubatureInfo,rnonlinearSpatialMatrix%rdiscrData%rstabilDual,&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplatesOptC%rmatrixEOJ2)

        ! Reintegrate the computed matrix
        call lsysbl_moveToSubmatrix (rtempMatrix,rmatrix,4,4)

        ! Include the mass matrix blocks
        !
        !    (  R2   .    .    .    .    .  )
        !    (  .    R2   .    .    .    .  )
        !    (  .    .    .    .    .    .  )
        !
        ! Remember that they are nonlinear if we have a preconditioner with
        ! Newton activated!
        
        call lsysbl_extractSubmatrix (rmatrix,rtempMatrix,4,6,1,3)
        
        ! Assemble linear parts.
        if (.not. associated(rnonlinearSpatialMatrix%rdiscrData%p_DobservationArea)) then
          ! No specific observation area
          call assembleLinearSubmatrix (rnonlinearSpatialMatrix,cmatrixType,rflags,rtempMatrix,&
              rnonlinearSpatialMatrix%DidentityY(2,1),rnonlinearSpatialMatrix%Dmass(2,1),&
              rnonlinearSpatialMatrix%Dstokes(2,1),&
              rnonlinearSpatialMatrix%DBmat(2,1),rnonlinearSpatialMatrix%DBTmat(2,1))
        else
          ! A special observation area. Two steps. In the second, assemble the
          ! nonlinear mass matrix.
          call assembleLinearSubmatrix (rnonlinearSpatialMatrix,cmatrixType,rflags,rtempMatrix,&
              rnonlinearSpatialMatrix%DidentityY(2,1),0.0_DP,&
              rnonlinearSpatialMatrix%Dstokes(2,1),&
              rnonlinearSpatialMatrix%DBmat(2,1),rnonlinearSpatialMatrix%DBTmat(2,1))
              
          call assembleBoxIdentity (rnonlinearSpatialMatrix%rdiscrData%p_DobservationArea,&
              rnonlinearSpatialMatrix%Dmass(2,1),rtempMatrix%RmatrixBlock(1,1),rcubatureInfo)
          call assembleBoxIdentity (rnonlinearSpatialMatrix%rdiscrData%p_DobservationArea,&
              rnonlinearSpatialMatrix%Dmass(2,1),rtempMatrix%RmatrixBlock(2,2),rcubatureInfo)
    
          ! Switch on the mass matrices
          rtempMatrix%RmatrixBlock(1,1)%dscaleFactor = 1.0_DP
          rtempMatrix%RmatrixBlock(2,2)%dscaleFactor = 1.0_DP
        end if

        ! Co stabilisation in the convective parts here.
        ! rstabilisation = t_convecStabilisation(&
        !    rnonlinearSpatialMatrix%iupwind2,rnonlinearSpatialMatrix%dupsam2)
        rstabilisation = t_settings_stabil(3,0.0_DP,1,1,1)

        call assembleNonlinearity (&
            rnonlinearSpatialMatrix,rtempMatrix,&
            rvectorPrimal,rvectorDual,&
            rnonlinearSpatialMatrix%Dstokes(2,1),&
            dweightConvection*rnonlinearSpatialMatrix%Dygrad(2,1),&
            dweightConvection*rnonlinearSpatialMatrix%DygradAdj(2,1),&
            dweightConvection*rnonlinearSpatialMatrix%DygradT(2,1),&
            dweightConvection*rnonlinearSpatialMatrix%DygradTAdj(2,1),&
            dweightConvection*rnonlinearSpatialMatrix%Dgrady(2,1),&
            dweightConvection*rnonlinearSpatialMatrix%DgradyT(2,1),&
            dweightConvection*rnonlinearSpatialMatrix%DgradyAdj(2,1),&
            rcubatureInfo,rstabilisation)
            
        ! There is probably a 2nd reactive term stemming from the next time step.
        ! Assemble it.
        
        call assembleNonlinearity (&
            rnonlinearSpatialMatrix,rtempMatrix,&
            rvectorPrimal,rvectorDual2,&
            0.0_DP,0.0_DP,0.0_DP,&
            dweightConvection*rnonlinearSpatialMatrix%DygradT2(2,1),&
            dweightConvection*rnonlinearSpatialMatrix%DygradTAdj2(2,1),&
            dweightConvection*rnonlinearSpatialMatrix%Dgrady2(2,1),&
            0.0_DP,dweightConvection*rnonlinearSpatialMatrix%DgradyAdj2(2,1),&
            rcubatureInfo,rstabilisation)

        ! Reintegrate the computed matrix
        call lsysbl_moveToSubmatrix (rtempMatrix,rmatrix,4,1)

        ! Switch the I-matrix in the continuity equation on/off.
        ! The matrix always exists -- but is usually filled with zeroes
        ! or switched off.
        rmatrix%RmatrixBlock(3,3)%dscaleFactor = rnonlinearSpatialMatrix%DidentityP(1,1)
        if (rnonlinearSpatialMatrix%DidentityP(1,1) .ne. 0.0_DP) then
          call lsyssc_initialiseIdentityMatrix (rmatrix%RmatrixBlock(3,3))
        end if
        
        rmatrix%RmatrixBlock(6,6)%dscaleFactor = rnonlinearSpatialMatrix%DidentityP(2,2)
        if (rnonlinearSpatialMatrix%DidentityP(2,2) .ne. 0.0_DP) then
          call lsyssc_initialiseIdentityMatrix (rmatrix%RmatrixBlock(6,6))
        end if

        ! Matrix restriction
        ! ---------------------------------------------------

        ! For the construction of matrices on lower levels, call the matrix
        ! restriction. In case we have a uniform discretisation with Q1~,
        ! iadaptivematrix may be <> 0 and so this will rebuild some matrix entries
        ! by a Galerkin approach using constant prolongation/restriction.
        ! This helps to stabilise the solver if there are elements in the
        ! mesh with high aspect ratio.
        if (present(rfineMatrix)) then
        
          ! Primal system:
          call mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(1,1), &
              rmatrix%RmatrixBlock(1,1), &
              rflags%iadaptiveMatrices, &
              rflags%dadmatthreshold)
              
          if (.not. lsyssc_isMatrixContentShared(&
              rfineMatrix%RmatrixBlock(1,1),rfineMatrix%RmatrixBlock(2,2))) then
            call mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(2,2), &
                rmatrix%RmatrixBlock(1,1), &
                rflags%iadaptiveMatrices, &
                rflags%dadmatthreshold)
          end if
            
          ! Dual system
          call mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(4,4), &
              rmatrix%RmatrixBlock(4,4), &
              rflags%iadaptiveMatrices, &
              rflags%dadmatthreshold)
              
          if (.not. lsyssc_isMatrixContentShared(&
            rmatrix%RmatrixBlock(4,4),rmatrix%RmatrixBlock(5,5))) then
            call mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(5,5), &
                rmatrix%RmatrixBlock(5,5), &
                rflags%iadaptiveMatrices, &
                rflags%dadmatthreshold)
          end if
          
        end if

        ! Release memory
        call lsysbl_releaseMatrix(rtempMatrix)
        call spdiscr_releaseBlockDiscr(rvelDiscr)

        !call matio_writeBlockMatrixHR (rmatrix, 'matrix',&
        !    .true., 0, 'matrix1.txt', '(E12.5)', 1E-10_DP)
        
      else

        ! New implementation: Use the 'direct' assembly routine.
        
        rmatrix%RmatrixBlock(1:2,1:2)%dscaleFactor = 1.0_DP
        rmatrix%RmatrixBlock(4:5,4:5)%dscaleFactor = 1.0_DP
        rmatrix%RmatrixBlock(1:2,4:5)%dscaleFactor = 1.0_DP
        rmatrix%RmatrixBlock(4:5,1:2)%dscaleFactor = 1.0_DP

        ! Clear the old matrix
        if (iand(coperation,CMASM_NOCLEAR) .eq. 0) then
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,1))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,2))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,1))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,2))
          
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(4,4))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(4,5))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(5,4))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(5,5))
          
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,4))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,4))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,5))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,5))
          
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(4,1))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(4,2))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(5,1))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(5,2))
          
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(3,3))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(6,6))
        end if
        
        ! The B1/B2 matrices exist up to now only in rnonlinearSpatialMatrix.
        ! Put a copy of them into the block matrix.
        !
        ! Note that we share the structure of B1/B2 with those B1/B2 of the
        ! block matrix, while we create copies of the entries. The B-blocks
        ! are already prepared and memory for the entries is already allocated;
        ! so we only have to copy the entries.
        !
        ! Note that idubContent = LSYSSC_DUP_COPY will automatically allocate
        ! memory if necessary.
        if (rnonlinearSpatialMatrix%DBmat(1,1) .ne. 0.0_DP) then
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
              rmatrix%RmatrixBlock(1,3),LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)

          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
              rmatrix%RmatrixBlock(2,3),LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
        end if
        
        if (rnonlinearSpatialMatrix%DBTmat(1,1) .ne. 0.0_DP) then
          ! Furthermore, put B1^T and B2^T to the block matrix.
          ! These matrices are always 'shared'.
          if (rflags%bvirtualTransposedD) then
            call lsyssc_transposeMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
                rmatrix%RmatrixBlock(3,1),LSYSSC_TR_VIRTUAL)
            call lsyssc_transposeMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
                rmatrix%RmatrixBlock(3,2),LSYSSC_TR_VIRTUAL)
          else
            call lsyssc_duplicateMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD1, &
                rmatrix%RmatrixBlock(3,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

            call lsyssc_duplicateMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD2, &
                rmatrix%RmatrixBlock(3,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
          end if
        end if
                                        
        if (rnonlinearSpatialMatrix%DBTmat(2,2) .ne. 0.0_DP) then
          if (rflags%bvirtualTransposedD) then
            call lsyssc_transposeMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
                rmatrix%RmatrixBlock(6,4),LSYSSC_TR_VIRTUAL)
            call lsyssc_transposeMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
                rmatrix%RmatrixBlock(6,5),LSYSSC_TR_VIRTUAL)
          else
            call lsyssc_duplicateMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD1, &
                rmatrix%RmatrixBlock(6,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

            call lsyssc_duplicateMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD2, &
                rmatrix%RmatrixBlock(6,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
          end if
        end if

        rmatrix%RmatrixBlock(1,3)%dscaleFactor = rnonlinearSpatialMatrix%DBmat(1,1)
        rmatrix%RmatrixBlock(2,3)%dscaleFactor = rnonlinearSpatialMatrix%DBmat(1,1)

        rmatrix%RmatrixBlock(4,6)%dscaleFactor = rnonlinearSpatialMatrix%DBmat(2,2)
        rmatrix%RmatrixBlock(5,6)%dscaleFactor = rnonlinearSpatialMatrix%DBmat(2,2)
        
        rmatrix%RmatrixBlock(3,1)%dscaleFactor = rnonlinearSpatialMatrix%DBTmat(1,1)
        rmatrix%RmatrixBlock(3,2)%dscaleFactor = rnonlinearSpatialMatrix%DBTmat(1,1)

        rmatrix%RmatrixBlock(6,4)%dscaleFactor = rnonlinearSpatialMatrix%DBTmat(2,2)
        rmatrix%RmatrixBlock(6,5)%dscaleFactor = rnonlinearSpatialMatrix%DBTmat(2,2)

        ! ---------------------------------------------------
        ! Now a slightly more advanced task for which we use a separate
        ! routine and some submatrices/vectors: The nonlinearity.

        ! Initialise the operator structure for what we need.
        roptcoperator%dupsamPrimal = rnonlinearSpatialMatrix%rdiscrData%rstabilPrimal%dupsam
        roptcoperator%dupsamDual = rnonlinearSpatialMatrix%rdiscrData%rstabilDual%dupsam
        
        ! Timestep-weights
        roptcoperator%dprimalAlpha = rnonlinearSpatialMatrix%Dmass(1,1)
        roptcoperator%ddualAlpha   = rnonlinearSpatialMatrix%Dmass(2,2)

        ! Stokes operator
        roptcoperator%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst
        roptcoperator%dprimalBeta = rnonlinearSpatialMatrix%Dstokes(1,1)
        roptcoperator%ddualBeta   = rnonlinearSpatialMatrix%Dstokes(2,2)
        
        ! Nonlinearity
        if (rnonlinearSpatialMatrix%Dygrad(1,1) .ne. 0.0_DP) then
          roptcoperator%dprimalDelta = dweightConvection * rnonlinearSpatialMatrix%Dygrad(1,1)
          roptcoperator%ddualDelta   = dweightDualConvection * rnonlinearSpatialMatrix%Dygrad(2,2)
          roptcoperator%ddualNewtonTrans = dweightDualNewtonT * rnonlinearSpatialMatrix%DgradyT(2,2)
          
          ! Whether or not Newton is active has no influence to the
          ! defect, so the following lines are commented out.
          ! if (rparams%bnewton) then
          roptcoperator%dprimalNewton    = dweightConvection * rnonlinearSpatialMatrix%Dgrady(1,1)
          roptcoperator%ddualRDeltaTrans = dweightConvection * rnonlinearSpatialMatrix%DygradT(2,1)
          roptcoperator%ddualRNewton     = dweightConvection * rnonlinearSpatialMatrix%Dgrady(2,1)
          ! end if
          
        end if
        
        ! Coupling matrices
        !if (rparams%bdualcoupledtoprimal) then
          roptcoperator%ddualRAlpha = rnonlinearSpatialMatrix%Dmass(2,1)
        !end if

        !if (rparams%bcontrolactive) then
          roptcoperator%dcontrolWeight = &
              -rnonlinearSpatialMatrix%Dmass(1,2)*rnonlinearSpatialMatrix%dalphaC
          roptcoperator%dcontrolMultiplier = -1.0_DP/rnonlinearSpatialMatrix%dalphaC
        !end if
        
        if (rnonlinearSpatialMatrix%rdiscrData%rconstraints%ccontrolConstraints .ne. 0) then

          roptcoperator%ccontrolProjection = &
              rnonlinearSpatialMatrix%rdiscrData%rconstraints%ccontrolConstraints

          roptcoperator%cconstraintsType = &
              rnonlinearSpatialMatrix%rdiscrData%rconstraints%cconstraintsType
              
          roptcoperator%dmin1 = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1
          roptcoperator%dmax1 = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1
          roptcoperator%dmin2 = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2
          roptcoperator%dmax2 = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2

          roptcoperator%p_rumin1 => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rumin1
          roptcoperator%p_rumax1 => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rumax1
          roptcoperator%p_rumin2 => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rumin2
          roptcoperator%p_rumax2 => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rumax2

          roptcoperator%p_rvectorumin => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin
          roptcoperator%p_rvectorumax => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax
        end if
        
        select case (rnonlinearSpatialMatrix%iprimalSol)
        case (1)
          p_rprimalSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1
        case (2)
          p_rprimalSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2
        case (3)
          p_rprimalSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3
        end select
        
        select case (rnonlinearSpatialMatrix%idualSol)
        case (1)
          p_rdualSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1
        case (2)
          p_rdualSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2
        case (3)
          p_rdualSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3
        end select
        
        ! Calculate the velocity-dependent part of the system matrix.
        call conv_strdiffOptC2dgetMatrix (rmatrix,roptcoperator,1.0_DP,&
            p_rprimalSol,p_rdualSol,rcubatureInfo=rcubatureInfo)
            
        !call matio_writeBlockMatrixHR (rmatrix, 'matrix',&
        !    .true., 0, 'matrix2.txt', '(E12.5)', 1E-10_DP)
      
      end if
      
      ! Release cubature related stuff
      call spdiscr_releaseCubStructure(rcubatureInfo)
      
    end if
    
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,3),p_Ddata)
    
    if (rvectorPrimal%NEQ .ne. 0) &
      call lsysbl_releaseVector (rvectorPrimal)

    if (rvectorDual%NEQ .ne. 0) &
      call lsysbl_releaseVector (rvectorDual)

    if (rvectorDual2%NEQ .ne. 0) &
      call lsysbl_releaseVector (rvectorDual2)
    
    
  contains
  
    ! -----------------------------------------------------
    
    subroutine allocSubmatrix (rnonlinearSpatialMatrix,cmatrixType,rflags,rsubmatrix,&
        diota,dalpha,dtheta,dygrad,dygradAdj,dygradT,dygradTAdj,dgrady,dgradyT,dgradyAdj,&
        deta,dtau,dkappa,ddirichletBCCY,ddirichletBCCLambda,ddirichletBCCXi)
        
    ! Allocates memory for a submatrix of the system matrix.
    ! The submatrix can then be inserted into a larger system matrix.
    
    ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix.
    type(t_nonlinearSpatialMatrix), intent(IN), target :: rnonlinearSpatialMatrix

    ! Type of matrix that should be set up in rmatrix. One of the CCMASM_MTP_xxxx
    ! constants.
    integer, intent(IN) :: cmatrixType

    ! Additional matrix assembly flags that configure how to generate the matrix.
    type(t_matrixAssemblyFlags), intent(in) :: rflags

    ! A block matrix that receives the basic system matrix.
    type(t_matrixBlock), intent(INOUT) :: rsubmatrix
    
    ! Coefficients in front of all operators in the system matrix.
    ! Depending on which operators are 'active', memory is allocated.
    real(DP), intent(in) :: diota,dalpha,dtheta,dygrad,dygradAdj,dygradT,dygradTAdj,&
        dgrady,dgradyT,dgradyAdj,deta,dtau,ddirichletBCCY,ddirichletBCCLambda,ddirichletBCCXi
        
    ! Coefficient in front of a pressure matrix at position (3,3).
    ! Switches that pressure-matrix on/off.
    real(DP), intent(in) :: dkappa
       
      ! local variables
      logical :: bdecoupled,bfulltensor
       
      ! Determine the shape of the matrix
      bdecoupled = (cmatrixType .eq. CCMASM_MTP_DECOUPLED) .or. &
                   (ddirichletBCCLambda .ne. 0.0_DP)
      bfulltensor = cmatrixType .eq. CCMASM_MTP_FULLTENSOR
      
      if (cmatrixType .eq. CCMASM_MTP_AUTOMATIC) then
        ! Should we assemble Newton? If yes, we have a full-tensor matrix.
        bfulltensor = (dgrady .ne. 0.0_DP) .or. bfulltensor
      end if
    
      ! Let's consider the global system in detail. The standard matrix It has
      ! roughly the following shape:
      !
      !    ( A11  .    B1  )
      !    ( .    A22  B2  )
      !    ( B1^T B2^T C   )
      !
      ! All matices may have multiplication factors in their front.
      !
      ! The structure of the matrices A11 and A22 of the global system matrix
      ! is governed by the template FEM matrix.
      ! Initialise them with the same structure, i.e. A11, A22 share (!) their
      ! structure (not the entries) with that of the template matrix.
      !
      ! For this purpose, use the "duplicate matrix" routine.
      ! The structure of the matrix is shared with the template FEM matrix.
      ! For the content, a new empty array is allocated which will later receive
      ! the entries.
      call lsyssc_duplicateMatrix (&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEM,&
          rsubmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
          
      if (.not. bfulltensor) then
           
        ! If X- and Y-velocity is to be treated in a 'coupled' way, the matrix
        ! A22 is identical to A11! So mirror A11 to A22 sharing the
        ! structure and the content.
        call lsyssc_duplicateMatrix (rsubmatrix%RmatrixBlock(1,1),&
                    rsubmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
                    
      else
      
        ! Otherwise, create another copy of the template matrix.
        call lsyssc_duplicateMatrix (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEM,&
            rsubmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                    
      end if
      
      ! A 'full tensor matrix' consists also of blocks A12 and A21.
      if (bfulltensor) then

        ! We have a matrix in the following shape:
        !
        !    ( A11  A12  B1  )
        !    ( A21  A22  B2  )
        !    ( B1^T B2^T C   )
        !
        ! Create A12 and A21.
      
        if (rsubmatrix%RmatrixBlock(1,2)%cmatrixFormat &
            .eq. LSYSSC_MATRIXUNDEFINED) then
            
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEM, &
              rsubmatrix%RmatrixBlock(1,2), LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
            
          ! Allocate memory for the entries; don't initialise the memory.
          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
          ! zero.
          ! call lsyssc_allocEmptyMatrix (&
          !     rsubmatrix%RmatrixBlock(1,2),LSYSSC_SETM_UNDEFINED)
            
        end if

        if (rsubmatrix%RmatrixBlock(2,1)%cmatrixFormat &
            .eq. LSYSSC_MATRIXUNDEFINED) then
            
          ! Create a new matrix A21 in memory. create a new matrix
          ! using the template FEM matrix...
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEM, &
              rsubmatrix%RmatrixBlock(2,1), LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
            
          ! Allocate memory for the entries; don't initialise the memory.
          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
          ! zero.
          ! call lsyssc_allocEmptyMatrix (&
          !     rsubmatrix%RmatrixBlock(2,1),LSYSSC_SETM_UNDEFINED)
            
        end if
        
      end if

      ! The B1/B2 matrices exist up to now only in our local problem structure.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2 with those B1/B2 of the
      ! block matrix, while we create empty space for the entries.
      ! Later, the B-matrices are copied into here and modified for boundary
      ! conditions.
      !
      ! If cdirichletBCCXi>0, we may have boundary control and need that
      ! matrix pair as well.
      if (deta*deta + ddirichletBCCXi*ddirichletBCCXi .ne. 0.0_DP) then
        call lsyssc_duplicateMatrix (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
            rsubmatrix%RmatrixBlock(1,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

        call lsyssc_duplicateMatrix (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
            rsubmatrix%RmatrixBlock(2,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      end if
        
      ! In the same manner, insert an identiy matrix for the pressure
      ! to the system matrix. The matrix is zero by default but may
      ! partially or fully be initialised to an identity matrix depending
      ! on the situation. (It's mostly used for direct solvers/UMFPACK)
      if (dkappa .ne. 0.0_DP) then
        call lsyssc_duplicateMatrix (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEMPressure, &
            rsubmatrix%RmatrixBlock(3,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      end if

      ! Now, prepare B1^T and B2^T. Are these matrices given?
      !
      ! If yes, take the given matrices. If not, create them by
      ! 'virtually transposing' B1 and B2 (i.e. they share the same
      ! data as B1 and B2 but hate the 'transpose'-flag set).

      if (dtau .ne. 0.0_DP) then
        if (rflags%bvirtualTransposedD) then
          call lsyssc_transposeMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
              rsubmatrix%RmatrixBlock(3,1),LSYSSC_TR_VIRTUAL)

          call lsyssc_transposeMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
              rsubmatrix%RmatrixBlock(3,2),LSYSSC_TR_VIRTUAL)
        else
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD1, &
              rsubmatrix%RmatrixBlock(3,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD2, &
              rsubmatrix%RmatrixBlock(3,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        end if
      end if

    end subroutine

    ! -----------------------------------------------------
    
    subroutine assembleLinearSubmatrix (rnonlinearSpatialMatrix,cmatrixType,rflags,&
        rmatrix,diota,dalpha,dtheta,deta,dtau)
        
    ! Assembles a matrix containing only linear terms. The matrix is a
    ! submatrix of a larger system matrix and can then be inserted
    ! into a larger system matrix.
    
    ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix.
    type(t_nonlinearSpatialMatrix), intent(IN) :: rnonlinearSpatialMatrix

    ! Type of matrix that should be set up in rmatrix. One of the CCMASM_MTP_xxxx
    ! constants.
    integer, intent(IN) :: cmatrixType

    ! Additional matrix assembly flags that configure how to generate the matrix.
    type(t_matrixAssemblyFlags), intent(in) :: rflags

    ! A block matrix that receives the matrix.
    type(t_matrixBlock), intent(INOUT) :: rmatrix
    
    ! Coefficient in front of an identity matrix.
    real(DP), intent(in) :: diota
    
    ! Coefficient in front of a mass matrix.
    real(DP), intent(in) :: dalpha
    
    ! Coefficient in front of the Laplace matrix.
    real(DP), intent(in) :: dtheta
    
    ! Coefficient in front of the B-matrices
    real(DP), intent(in) :: deta
    
    ! Coefficient in front of the B^T-matrices
    real(DP), intent(in) :: dtau
        
    ! local variables
    logical :: bshared
    
      ! Is A11=A22 physically?
      bshared = lsyssc_isMatrixContentShared(&
                    rmatrix%RmatrixBlock(1,1),rmatrix%RmatrixBlock(2,2))
                    
      if ((diota .eq. 0.0_DP) .and. (dalpha .eq. 0.0_DP) .and.&
          (dtheta .eq. 0.0_DP)) then
          
        ! Switch off that matrices
        rmatrix%RmatrixBlock(1,1)%dscaleFactor = 0.0_DP
        rmatrix%RmatrixBlock(2,2)%dscaleFactor = 0.0_DP
          
      else
                    
        rmatrix%RmatrixBlock(1,1)%dscaleFactor = 1.0_DP
        rmatrix%RmatrixBlock(2,2)%dscaleFactor = 1.0_DP
      
        ! Let's consider the global system in detail:
        !
        !    ( A11  A12  B1  ) = ( A11  A12  A13 )
        !    ( A21  A22  B2  )   ( A21  A22  A23 )
        !    ( B1^T B2^T 0   )   ( A31  A32  A33 )
        !
        ! Here, we set up the linear parts in Aij and Bi.
        
        ! ---------------------------------------------------
        ! If diota <> 0, initialise the matrix with the identity,
        ! otherwise with zero.
        if (diota .eq. 0.0_DP) then
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,1))
          
          if (.not. bshared) then
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,2))
          end if
          
        else
          
          call lsyssc_initialiseIdentityMatrix (&
              rmatrix%RmatrixBlock(1,1))
          call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(1,1),diota)
          
          if (.not. bshared) then
            call lsyssc_initialiseIdentityMatrix (&
                rmatrix%RmatrixBlock(2,2))
            call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(2,2),diota)
          end if
          
        end if
      
        ! ---------------------------------------------------
        ! Plug in the mass matrix?
        if (dalpha .ne. 0.0_DP) then
         
          call lsyssc_matrixLinearComb (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rmatrix%RmatrixBlock(1,1),&
              dalpha,1.0_DP,&
              .false.,.false.,.true.,.true.,rmatrix%RmatrixBlock(1,1))
              
          if (.not. bshared) then

            call lsyssc_matrixLinearComb (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
                rmatrix%RmatrixBlock(2,2),&
                dalpha,1.0_DP,&
                .false.,.false.,.true.,.true.,rmatrix%RmatrixBlock(2,2))
          end if
          
        end if
        
        ! ---------------------------------------------------
        ! Plug in the Stokes matrix?
        if (dtheta .ne. 0.0_DP) then
          ! Plug in the Stokes matrix in case of the gradient tensor.
          ! In case of the deformation tensor or nonconstant viscosity,
          ! that is done during the assembly of the nonlinearity.
          if ((rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%isubequation .eq. 0) .and. &
              (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%cviscoModel .eq. 0)) then
        
            call lsyssc_matrixLinearComb (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixLaplace,&
                rmatrix%RmatrixBlock(1,1),&
                rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst * dtheta,1.0_DP,&
                .false.,.false.,.true.,.true.,rmatrix%RmatrixBlock(1,1))
                
            if (.not. bshared) then
              call lsyssc_matrixLinearComb (&
                  rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixLaplace,&
                  rmatrix%RmatrixBlock(2,2),&
                  rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst * dtheta,1.0_DP,&
                  .false.,.false.,.true.,.true.,rmatrix%RmatrixBlock(2,2))
            end if
          end if
          
        end if
        
      end if
      
      ! We exclude the velocity submatrices here, so our system looks like:
      !
      !    (           B1 ) = (           A13 )
      !    (           B2 )   (           A23 )
      !    ( B1^T B2^T    )   ( A31  A32      )

      ! The B1/B2 matrices exist up to now only in rnonlinearSpatialMatrix.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2 with those B1/B2 of the
      ! block matrix, while we create copies of the entries. The B-blocks
      ! are already prepared and memory for the entries is already allocated;
      ! so we only have to copy the entries.
      !
      ! Note that idubContent = LSYSSC_DUP_COPY will automatically allocate
      ! memory if necessary.
      if (deta .ne. 0.0_DP) then
        call lsyssc_duplicateMatrix (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
            rmatrix%RmatrixBlock(1,3),LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)

        call lsyssc_duplicateMatrix (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
            rmatrix%RmatrixBlock(2,3),LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
      end if
      
      if (dtau .ne. 0.0_DP) then
        ! Furthermore, put B1^T and B2^T to the block matrix.
        if (rflags%bvirtualTransposedD) then
          ! Virtually transposed matrices are always 'shared'.
          call lsyssc_transposeMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
              rmatrix%RmatrixBlock(3,1),LSYSSC_TR_VIRTUAL)
          call lsyssc_transposeMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
              rmatrix%RmatrixBlock(3,2),LSYSSC_TR_VIRTUAL)
        else
          ! Real matrices are not shared. This allows them to be overwritten
          ! if necessary (e.g. if an UMFPACK solver is used or similar).
          ! In a situation, where the matrix needs to be overwritten,
          ! the caller must therefore always request a not virtually transposed
          ! matix!
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD1, &
              rmatrix%RmatrixBlock(3,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPYOVERWRITE)

          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD2, &
              rmatrix%RmatrixBlock(3,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPYOVERWRITE)
        end if
      end if

      rmatrix%RmatrixBlock(1,3)%dscaleFactor = deta
      rmatrix%RmatrixBlock(2,3)%dscaleFactor = deta
      
      rmatrix%RmatrixBlock(3,1)%dscaleFactor = dtau
      rmatrix%RmatrixBlock(3,2)%dscaleFactor = dtau

    end subroutine
    
    ! -----------------------------------------------------

    subroutine assembleNonlinearity (&
        rnonlinearSpatialMatrix,rmatrix,rvectorPrimalVel,rvectorNonlinearity,&
        dtheta,dygrad,dygradAdj,dygradT,dygradTadj,dgrady,dgradyT,dgradyAdj,&
        rcubatureInfo,rstabilisation,rmatrixEOJ)
        
    ! Assembles the convection matrix in the block matrix rmatrix at position (1,1):
    !
    ! rmatrix := dygrad*N(rvector) + dygradT*N^t(rvector) +
    !            dgrady*N*(rvector) + dgradyT*N*^t(rvector)
    !
    ! Even if no nonlinearity is present, the routine can be used to
    ! add stabilisation into the matrix.
    
    ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix.
    type(t_nonlinearSpatialMatrix), intent(in) :: rnonlinearSpatialMatrix
    
    ! Block matrix where the 2x2-velocity submatrix should be assembled
    type(t_matrixBlock), intent(inout) :: rmatrix
    
    ! Primal velocity vector. Used for nonlinear viscosity.
    type(t_vectorBlock), intent(in) :: rvectorPrimalVel
    
    ! Velocity vector for the nonlinearity. Must be specified if
    ! GAMMA <> 0; can be omitted if GAMMA=0.
    type(t_vectorBlock), intent(in) :: rvectorNonlinearity
    
    ! Weight for the stabilisation and Stokes operator in the Stokes case.
    real(DP), intent(in) :: dtheta
    
    ! Weight for the nonlinear term u\grad(.)
    real(DP), intent(in) :: dygrad

    ! Weight for the adjoint nonlinear term (.,u\grad(.))
    real(DP), intent(in) :: dygradAdj

    ! Weight for the nonlinear term u(\grad(.))^t
    real(DP), intent(in) :: dygradT

    ! Weight for the adjoint nonlinear term (.,u(\grad(.))^t)
    real(DP), intent(in) :: dygradTAdj

    ! Weight for the nonlinear term (\grad(.))u
    real(DP), intent(in) :: dgrady

    ! Weight for the nonlinear term (\grad(.))^t u
    real(DP), intent(in) :: dgradyT

    ! Weight for the adjoint of the nonlinear term (\grad(.))u
    real(DP), intent(in) :: dgradyAdj
    
    ! Cubature info structure that defines how to apply cubature.
    type(t_scalarCubatureInfo), intent(in) :: rcubatureInfo
    
    ! Stabilisation parameters.
    type(t_settings_stabil), intent(in) :: rstabilisation
    
    ! OPTIONAL: EOJ matrix in case the precomputed EOJ stabilisation is active.
    type(t_matrixScalar), intent(in), optional :: rmatrixEOJ
    
    ! local variables
    logical :: bshared
    type(t_convUpwind) :: rupwind
    type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
    type(t_convStreamDiff2) :: rstreamlineDiffusion2
    type(t_convStreamDiff2) :: rstreamlineDiffusion3
    type(t_jumpStabilisation) :: rjumpStabil
    type(t_convUpwind) :: rupwindStabil
    real(dp), dimension(:), pointer :: p_Ddata1,p_Ddata2,p_Ddata3
    integer, dimension(:), pointer :: p_Idofs
    integer :: h_Idofs
    type(t_matrixBlock) :: rmatrixtemp
    type(t_collection), target :: rcollection, ruserCollection
    
      ! Is A11=A22 physically?
      bshared = lsyssc_isMatrixContentShared(&
                    rmatrix%RmatrixBlock(1,1),&
                    rmatrix%RmatrixBlock(2,2))
                    
      if (rstabilisation%cconvectionOnBoundaryMatrix .eq. 0) then
        ! Copy the matrix so we can restore rows if necessary
        call lsysbl_duplicateMatrix (rmatrix,rmatrixtemp,&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      end if
                    
      call lsysbl_getbase_double (rvectorNonlinearity,p_Ddata1)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,1),p_Ddata2)
             
      if ((dygrad .ne. 0.0_DP) .or. (dygradT .ne. 0.0_DP) .or. &
          (dygradAdj .ne. 0.0_DP) .or. (dygradTAdj .ne. 0.0_DP) .or. &
          (dgrady .ne. 0.0_DP) .or. (dgradyT .ne. 0.0_DP) .or. (dgradyAdj .ne. 0.0_DP) .or. &
          (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%isubequation .ne. 0) .or. &
          (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%cviscoModel .ne. 0)) then
                  
        ! Prepare a user defined collection structure for the
        ! callback routines
        call collct_init (ruserCollection)
        call user_initCollectForAssembly (rnonlinearSpatialMatrix%p_rglobalData,&
            rnonlinearSpatialMatrix%p_rnonlinearity%dassociatedTimePrimalVel,ruserCollection)
                      
        ! Switch on the offdiagonal matrices if necessary
        if ((dygradT .ne. 0.0_DP) .or. (dygradTAdj .ne. 0.0_DP) .or. &
            (dgrady .ne. 0.0_DP) .or. (dgradyT .ne. 0.0_DP)) then
          ! Switch on missing matrices. Initialise with zero when they
          ! are created the first time.
          if (rmatrix%RmatrixBlock(1,2)%dscaleFactor .eq. 0.0_DP) then
            rmatrix%RmatrixBlock(1,2)%dscaleFactor = 1.0_DP
            rmatrix%RmatrixBlock(2,1)%dscaleFactor = 1.0_DP
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,2))
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,1))
          end if
        else
          ! If the submatrices A12 and A21 exist, fill them with zero.
          ! If they don't exist, we don't have to do anything.
          if (lsysbl_isSubmatrixPresent (rmatrix,1,2)) then
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,2))
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,1))
          end if
        end if
                      
        select case (rstabilisation%cupwind)
        case (CCMASM_STAB_STREAMLINEDIFF)
        
          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%cviscoModel .ne. 0) then
            ! Not supported by this SD method.
            call output_line (&
                "This assembly method does not support nonconstant viscosity!", &
                OU_CLASS_ERROR,OU_MODE_STD,"cc_assembleMatrix")
            call sys_halt()
          end if

          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst
          
          ! Set stabilisation parameter
          rstreamlineDiffusion%dupsam = rstabilisation%dupsam
          
          ! Matrix weights
          rstreamlineDiffusion%ddelta            = dygrad
          rstreamlineDiffusion%ddeltaTransposed  = dygradT
          rstreamlineDiffusion%dnewton           = dgrady
          rstreamlineDiffusion%dnewtonTransposed = dgradyT
          
          if (rstabilisation%dupsam .ne. 0) then
            call output_line ("assembleNonlinearity: Warning. Please use the")
            call output_line ("alternative SD method for setting up the stabilised operator!!!")
          end if

          if ((dygradAdj .ne. 0.0_DP) .or. (dygradTAdj .ne. 0.0_DP)) then
            call output_line ("Operator not supported by this SD method!", &
                OU_CLASS_ERROR,OU_MODE_STD,"cc_assembleMatrix")
            call sys_halt()
          end if
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvectorNonlinearity!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamlineDiffusionBlk2d (&
                              rvectorNonlinearity, rvectorNonlinearity, &
                              1.0_DP, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODMATRIX, &
                              rmatrix,rcubatureInfo=rcubatureInfo)
                              
        case (CCMASM_STAB_STREAMLINEDIFF2)
        
          ! Set up the SD structure.
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst
          
          ! Set stabilisation parameter
          rstreamlineDiffusion2%dupsam = rstabilisation%dupsam
          
          ! Matrix weights
          rstreamlineDiffusion2%ddelta  = dygrad
          rstreamlineDiffusion2%ddeltaT = dygradT
          rstreamlineDiffusion2%ddeltaAdj  = dygradAdj
          rstreamlineDiffusion2%ddeltaTAdj = dygradTAdj
          rstreamlineDiffusion2%dnewton = dgrady
          rstreamlineDiffusion2%dnewtonT = dgradyT
          rstreamlineDiffusion2%dnewtonAdj = dgradyAdj

          ! Probably, we have nonconstant viscosity.
          ! In that case, init a collection structure for a callback
          ! routine that specifies the viscosity
          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%cviscoModel .ne. 0) then
          
            rstreamlineDiffusion2%bconstnu = .false.
            
            ! Assemble at least the Stokes matrix.
            ! For a more complicated formulation, dbeta/dbetaT
            ! may be changed later.
            rstreamlineDiffusion2%dbeta = dtheta
            
            ! Prepare the collection. The "next" collection points to the user defined
            ! collection.
            rcollection%p_rnextCollection => ruserCollection
            call smva_prepareViscoAssembly (&
                rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal,&
                rcollection,rvectorPrimalVel)
            
          end if

          ! Assemble the deformation tensor?
          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%isubequation .ne. 0) then
            rstreamlineDiffusion2%dbeta = 0.5_DP * dtheta
            rstreamlineDiffusion2%dbetaT = 0.5_DP * dtheta
          end if

          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvectorNonlinearity!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamDiff2Blk2dMat (&
              rstreamlineDiffusion2,rmatrix,rvectorNonlinearity,&
              ffunctionViscoModel,rcollection,rcubatureInfo)

        case (CCMASM_STAB_UPWIND)

          ! Set up the SD structure for the creation of the defect.
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst
          
          ! Set stabilisation parameter
          rstreamlineDiffusion2%dupsam = 0.0_DP
          
          ! Matrix weights
          rstreamlineDiffusion2%ddelta  = 0.0_DP
          rstreamlineDiffusion2%ddeltaT = dygradT
          rstreamlineDiffusion2%ddeltaAdj  = dygradAdj
          rstreamlineDiffusion2%ddeltaTAdj = dygradTAdj
          rstreamlineDiffusion2%dnewton = dgrady
          rstreamlineDiffusion2%dnewtonT = dgradyT
          rstreamlineDiffusion2%dnewtonAdj = dgradyAdj

          ! Probably, we have nonconstant viscosity.
          ! In that case, init a collection structure for a callback
          ! routine that specifies the viscosity
          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%cviscoModel .ne. 0) then
          
            rstreamlineDiffusion2%bconstnu = .false.
            
            ! Assemble at least the Stokes matrix.
            ! For a more complicated formulation, dbeta/dbetaT
            ! may be changed later.
            rstreamlineDiffusion2%dbeta = dtheta
            
            ! Prepare the collection. The "next" collection points to the user defined
            ! collection.
            rcollection%p_rnextCollection => ruserCollection
            call smva_prepareViscoAssembly (&
                rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal,&
                rcollection,rvectorPrimalVel)
            
          end if
                    
          ! Assemble the deformation tensor?
          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%isubequation .ne. 0) then
            rstreamlineDiffusion2%dbeta = 0.5_DP * dtheta
            rstreamlineDiffusion2%dbetaT = 0.5_DP * dtheta
          end if
          
          ! Call the SD method to calculate the nonlinearity for everything except
          ! for the convection. As velocity vector, specify rvectorNonlinearity!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamDiff2Blk2dMat (&
              rstreamlineDiffusion2,rmatrix,rvectorNonlinearity,&
              ffunctionViscoModel,rcollection,rcubatureInfo)

          ! Prepare the upwind structure for the assembly of the convection.
          ! Note: Stabilisation weight is actually wrong, but it is not possible
          ! to specify in rupwindStabil%dupsam whether the stabilisation
          ! is added or subtracted!
          rupwindStabil%dupsam = abs(rstabilisation%dupsam)
          rupwindStabil%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst
          rupwindStabil%dtheta = dygrad
          
          ! Apply the upwind operator.
          call conv_upwind2d (rvectorNonlinearity, rvectorNonlinearity, 1.0_DP, 0.0_DP,&
              rupwindStabil, CONV_MODMATRIX, rmatrix%RmatrixBlock(1,1))

          if (.not. bshared) then
            call conv_upwind2d (rvectorNonlinearity, rvectorNonlinearity, 1.0_DP, 0.0_DP,&
                rupwindStabil, CONV_MODMATRIX, rmatrix%RmatrixBlock(2,2))
          end if

          ! Prepare the upwind structure for the assembly of the convection.
          !rstreamlineDiffusion3%dupsam = rstabilisation%dupsam
          !rstreamlineDiffusion3%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
          !rstreamlineDiffusion3%dtheta = dygrad
          !rstreamlineDiffusion3%ddelta = 1.0_DP
          
          ! Apply the upwind operator.
          !call conv_streamDiff2Blk2dMat (rstreamlineDiffusion3,rmatrix,rvectorNonlinearity)
          
          !call output_line ('Upwind not supported.', &
          !                  OU_CLASS_ERROR,OU_MODE_STD,'assembleNonlinearity')
          !call sys_halt()

        case (CCMASM_STAB_EDGEORIENTED)
          ! Jump stabilisation.

          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%cviscoModel .ne. 0) then
            ! Not supported by this SD method.
            call output_line (&
                "This assembly method does not support nonconstant viscosity!", &
                OU_CLASS_ERROR,OU_MODE_STD,"cc_assembleMatrix")
            call sys_halt()
          end if

          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion%ddelta            = dygrad
          rstreamlineDiffusion%ddeltaTransposed  = dygradT
          rstreamlineDiffusion%dnewton           = dgrady
          rstreamlineDiffusion%dnewtonTransposed = dgradyT
          
          if ((dygradAdj .ne. 0.0_DP) .or. (dygradTAdj .ne. 0.0_DP)) then
            call output_line ("Operator not supported by this SD method!", &
                OU_CLASS_ERROR,OU_MODE_STD,"cc_assembleMatrix")
            call sys_halt()
          end if
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvectorNonlinearity!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,1), p_Ddata1)
          call conv_streamlineDiffusionBlk2d (&
                              rvectorNonlinearity, rvectorNonlinearity, &
                              1.0_DP, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODMATRIX, &
                              rmatrix,rcubatureInfo=rcubatureInfo)
                              
          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rstreamlineDiffusion%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgamma = abs(rstabilisation%dupsam)
          
          ! Matrix weight. Compensate for any "-" sign in dygrad!
          rjumpStabil%dtheta = dygrad * mprim_signum(rstabilisation%dupsam)

          ! Call the jump stabilisation technique to stabilise that stuff.
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
                              rjumpStabil, CONV_MODMATRIX, &
                              rmatrix%RmatrixBlock(1,1))
          if (rstabilisation%ceojStabilOnBoundary .eq. 0) then
            call smva_addBdEOJOperator (rjumpStabil,-1.0_DP,rmatrix%RmatrixBlock(1,1))
          end if

          if (.not. bshared) then
            call conv_jumpStabilisation2d (&
                                rjumpStabil, CONV_MODMATRIX, &
                                rmatrix%RmatrixBlock(2,2))
            if (rstabilisation%ceojStabilOnBoundary .eq. 0) then
              call smva_addBdEOJOperator (rjumpStabil,-1.0_DP,rmatrix%RmatrixBlock(2,2))
            end if
          end if

        case (CCMASM_STAB_EDGEORIENTED2)
          ! Jump stabilisation.
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.

          ! Set up the SD structure.
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion2%dupsam = 0.0_DP
          
          ! Matrix weights
          rstreamlineDiffusion2%ddelta  = dygrad
          rstreamlineDiffusion2%ddeltaT = dygradT
          rstreamlineDiffusion2%ddeltaAdj  = dygradAdj
          rstreamlineDiffusion2%ddeltaTAdj = dygradTAdj
          rstreamlineDiffusion2%dnewton = dgrady
          rstreamlineDiffusion2%dnewtonT = dgradyT
          rstreamlineDiffusion2%dnewtonAdj = dgradyAdj

          ! Probably, we have nonconstant viscosity.
          ! In that case, init a collection structure for a callback
          ! routine that specifies the viscosity
          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%cviscoModel .ne. 0) then
          
            rstreamlineDiffusion2%bconstnu = .false.
            
            ! Assemble at least the Stokes matrix.
            ! For a more complicated formulation, dbeta/dbetaT
            ! may be changed later.
            rstreamlineDiffusion2%dbeta = dtheta
            
            ! Prepare the collection. The "next" collection points to the user defined
            ! collection.
            rcollection%p_rnextCollection => ruserCollection
            call smva_prepareViscoAssembly (&
                rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal,&
                rcollection,rvectorPrimalVel)
            
          end if
                    
          ! Assemble the deformation tensor?
          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%isubequation .ne. 0) then
            rstreamlineDiffusion2%dbeta = 0.5_DP * dtheta
            rstreamlineDiffusion2%dbetaT = 0.5_DP * dtheta
          end if
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvectorNonlinearity!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamDiff2Blk2dMat (&
              rstreamlineDiffusion2,rmatrix,rvectorNonlinearity,&
              ffunctionViscoModel,rcollection,rcubatureInfo)

          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rstreamlineDiffusion%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgamma = abs(rstabilisation%dupsam)
          
          ! Matrix weight. Compensate for any "-" sign in dygrad!
          rjumpStabil%dtheta = dygrad * mprim_signum(rstabilisation%dupsam)

          ! Call the jump stabilisation technique to stabilise that stuff.
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
                              rjumpStabil, CONV_MODMATRIX, &
                              rmatrix%RmatrixBlock(1,1))
          if (rstabilisation%ceojStabilOnBoundary .eq. 0) then
            call smva_addBdEOJOperator (rjumpStabil,-1.0_DP,rmatrix%RmatrixBlock(1,1))
          end if

          if (.not. bshared) then
            call conv_jumpStabilisation2d (&
                                rjumpStabil, CONV_MODMATRIX, &
                                rmatrix%RmatrixBlock(2,2))
            if (rstabilisation%ceojStabilOnBoundary .eq. 0) then
              call smva_addBdEOJOperator (rjumpStabil,-1.0_DP,rmatrix%RmatrixBlock(2,2))
            end if
          end if

        case (CCMASM_STAB_EDGEORIENTED3)
          ! Jump stabilisation.
          
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.

          ! Set up the SD structure.
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion2%dupsam = 0.0_DP
          
          ! Matrix weights
          rstreamlineDiffusion2%ddelta  = dygrad
          rstreamlineDiffusion2%ddeltaT = dygradT
          rstreamlineDiffusion2%ddeltaAdj  = dygradAdj
          rstreamlineDiffusion2%ddeltaTAdj = dygradTAdj
          rstreamlineDiffusion2%dnewton = dgrady
          rstreamlineDiffusion2%dnewtonT = dgradyT
          rstreamlineDiffusion2%dnewtonAdj = dgradyAdj

          ! Probably, we have nonconstant viscosity.
          ! In that case, init a collection structure for a callback
          ! routine that specifies the viscosity
          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%cviscoModel .ne. 0) then
          
            rstreamlineDiffusion2%bconstnu = .false.
            
            ! Assemble at least the Stokes matrix.
            ! For a more complicated formulation, dbeta/dbetaT
            ! may be changed later.
            rstreamlineDiffusion2%dbeta = dtheta
            
            ! Prepare the collection. The "next" collection points to the user defined
            ! collection.
            rcollection%p_rnextCollection => ruserCollection
            call smva_prepareViscoAssembly (&
                rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal,&
                rcollection,rvectorPrimalVel)
            
          end if
                    
          ! Assemble the deformation tensor?
          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%isubequation .ne. 0) then
            rstreamlineDiffusion2%dbeta = 0.5_DP * dtheta
            rstreamlineDiffusion2%dbetaT = 0.5_DP * dtheta
          end if
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvectorNonlinearity!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,1), p_Ddata1)
          call lsyssc_getbase_double (rmatrixEOJ, p_Ddata2)
          call conv_streamDiff2Blk2dMat (&
              rstreamlineDiffusion2,rmatrix,rvectorNonlinearity,&
              ffunctionViscoModel,rcollection,rcubatureInfo)

          ! We use the precomputed EOJ matrix and sum it up to the
          ! existing matrix.
          ! Matrix weight. Compensate for any "-" sign in dygrad!
          call lsyssc_matrixLinearComb(rmatrixEOJ,rmatrix%RmatrixBlock(1,1),&
              dygrad*mprim_signum(rstabilisation%dupsam),1.0_DP,&
              .false.,.false.,.true.,.true.,rmatrix%RmatrixBlock(1,1))

          if (.not. bshared) then
            ! Also for the Y-velocity.
            call lsyssc_matrixLinearComb(rmatrixEOJ,rmatrix%RmatrixBlock(2,2),&
                dygrad*mprim_signum(rstabilisation%dupsam),1.0_DP,&
                .false.,.false.,.true.,.true.,rmatrix%RmatrixBlock(2,2))
          end if

        case default
          print *,'Don''t know how to set up nonlinearity!?!'
          stop
        
        end select

        ! Release the collection stuff
        call user_doneCollectForAssembly (rnonlinearSpatialMatrix%p_rglobalData,ruserCollection)
        call collct_done (ruserCollection)
        
      else
      
        ! That's the Stokes-case. Jump stabilisation is possible...
      
        if ((rstabilisation%dupsam .ne. 0.0_DP)  .and. (dtheta .ne. 0.0_DP)) then

          ! Prepare a user defined collection structure for the
          ! callback routines
          call collct_init (ruserCollection)
          call user_initCollectForAssembly (rnonlinearSpatialMatrix%p_rglobalData,&
              rnonlinearSpatialMatrix%p_rnonlinearity%dassociatedTimePrimalVel,ruserCollection)

          ! Probably, we have nonconstant viscosity.
          ! In that case, init a collection structure for a callback
          ! routine that specifies the viscosity
          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%cviscoModel .ne. 0) then
          
            ! Set up the SD structure for the creation of the defect.
            rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst
            
            ! Set stabilisation parameter
            rstreamlineDiffusion2%dupsam = 0.0_DP
            
            ! Matrix weights
            rstreamlineDiffusion2%ddelta  = 0.0_DP
            rstreamlineDiffusion2%ddeltaT = dygradT
            rstreamlineDiffusion2%ddeltaAdj  = dygradAdj
            rstreamlineDiffusion2%ddeltaTAdj = dygradTAdj
            rstreamlineDiffusion2%dnewton = dgrady
            rstreamlineDiffusion2%dnewtonT = dgradyT
            rstreamlineDiffusion2%dnewtonAdj = dgradyAdj
  
            rstreamlineDiffusion2%bconstnu = .false.
            
            ! Assemble at least the Stokes matrix.
            ! For a more complicated formulation, dbeta/dbetaT
            ! may be changed later.
            rstreamlineDiffusion2%dbeta = dtheta
            
            ! Prepare the collection. The "next" collection points to the user defined
            ! collection.
            rcollection%p_rnextCollection => ruserCollection
            call smva_prepareViscoAssembly (&
                rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal,&
                rcollection,rvectorPrimalVel)
            
            ! Assemble the deformation tensor?
            if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%isubequation .ne. 0) then
              rstreamlineDiffusion2%dbeta = 0.5_DP * dtheta
              rstreamlineDiffusion2%dbetaT = 0.5_DP * dtheta
            end if
          
            ! Call the SD method to calculate the nonlinearity for everything except
            ! for the convection. As velocity vector, specify rvectorNonlinearity!
            ! Therefore, the primal velcity is always used for assembling
            ! that thing!
            call conv_streamDiff2Blk2dMat (&
                rstreamlineDiffusion2,rmatrix,rvectorNonlinearity,&
                ffunctionViscoModel,rcollection,rcubatureInfo)

          end if

          select case (rstabilisation%cupwind)
          case (CCMASM_STAB_EDGEORIENTED,CCMASM_STAB_EDGEORIENTED2)
            
            ! Set up the jump stabilisation structure.
            ! There's not much to do, only initialise the viscosity...
            rjumpStabil%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst
            
            ! Set stabilisation parameter
            rjumpStabil%dgamma = abs(rstabilisation%dupsam)
            
            ! Matrix weight. Take the weight of the Stokes operator.
            rjumpStabil%dtheta = dtheta

            ! Call the jump stabilisation technique to stabilise that stuff.
            ! We can assemble the jump part any time as it's independent of any
            ! convective parts...
            call conv_jumpStabilisation2d (&
                                rjumpStabil, CONV_MODMATRIX, &
                                rmatrix%RmatrixBlock(1,1))
            if (rstabilisation%ceojStabilOnBoundary .eq. 0) then
              call smva_addBdEOJOperator (rjumpStabil,-1.0_DP,rmatrix%RmatrixBlock(1,1))
            end if

            if (.not. bshared) then
              call conv_jumpStabilisation2d (&
                                  rjumpStabil, CONV_MODMATRIX, &
                                  rmatrix%RmatrixBlock(2,2))
              if (rstabilisation%ceojStabilOnBoundary .eq. 0) then
                call smva_addBdEOJOperator (rjumpStabil,-1.0_DP,rmatrix%RmatrixBlock(2,2))
              end if
            end if
            
          case (CCMASM_STAB_EDGEORIENTED3)
            
            ! We use the precomputed EOJ matrix and sum it up to the
            ! existing matrix.
            call lsyssc_matrixLinearComb(rmatrixEOJ,rmatrix%RmatrixBlock(1,1),&
                dtheta,1.0_DP,&
                .false.,.false.,.true.,.true.,rmatrix%RmatrixBlock(1,1))

            if (.not. bshared) then
              ! Also for the Y-velocity.
              call lsyssc_matrixLinearComb(rmatrixEOJ,rmatrix%RmatrixBlock(2,2),&
                  dtheta,1.0_DP,&
                  .false.,.false.,.true.,.true.,rmatrix%RmatrixBlock(2,2))
            end if
            
          case default
            ! No stabilisation
          
          end select

          ! Release the collection stuff
          call user_doneCollectForAssembly (rnonlinearSpatialMatrix%p_rglobalData,ruserCollection)
          call collct_done (ruserCollection)
        end if

      end if
      
      if (rstabilisation%cconvectionOnBoundaryMatrix .eq. 0) then
        ! Restore the matrix values on the boundary.
        h_Idofs = ST_NOHANDLE
        call bcasm_getDOFsOnBoundary (&
            rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimal%RspatialDiscr(1), h_Idofs)
        if (h_Idofs .ne. ST_NOHANDLE) then
          call storage_getbase_int (h_Idofs,p_Idofs)
          
          call lsyssc_replaceRowsMatrix9 (&
              rmatrixtemp%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,1),p_Idofs)

          call lsyssc_replaceRowsMatrix9 (&
              rmatrixtemp%RmatrixBlock(1,2),&
              rmatrix%RmatrixBlock(1,2),p_Idofs)

          call lsyssc_replaceRowsMatrix9 (&
              rmatrixtemp%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,1),p_Idofs)

          call lsyssc_replaceRowsMatrix9 (&
              rmatrixtemp%RmatrixBlock(2,2),&
              rmatrix%RmatrixBlock(2,2),p_Idofs)
              
          call storage_free (h_Idofs)
        end if
        call lsysbl_releaseMatrix (rmatrixtemp)
      end if
      
    end subroutine
      
    ! -----------------------------------------------------

    subroutine assembleProjectedMassBlocks (rnonlinearSpatialMatrix,rmatrix, &
        rvector, dweight)
        
    ! Assembles a 2x2 block matrix with
    ! probably nonlinear mass matrices on the diagonal.
    ! The mass matrices are projected according to where the velocity
    ! vector -1/alphaC*rvector%RvectorBlock(4/5) (from rnonlinearSpatialMatrix)
    ! violates the bounds dmin/dmax.
    
    ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix.
    type(t_nonlinearSpatialMatrix), intent(IN) :: rnonlinearSpatialMatrix
    
    ! Block matrix where the 2x2-velocity submatrix should be assembled.
    ! The matrix is cleared on call to this routine!
    type(t_matrixBlock), intent(INOUT) :: rmatrix
    
    ! Vector that specifies where to evaluate nonlinear terms
    type(t_vectorBlock), intent(IN), target :: rvector
    
    ! Weight for the (projected) mass matrices when adding them to the
    ! system.
    real(DP), intent(in) :: dweight
    
      ! local variables
      type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
      integer, dimension(:), pointer :: p_IelementList
      integer :: ccontrolConstraints
      type(t_scalarCubatureInfo) :: rcubatureInfo, rcubatureInfoAdapt

      ! Assemble A14/A25?
      if (dweight .ne. 0.0_DP) then
          
        ! Switch on these matrices
        rmatrix%RmatrixBlock(1,1)%dscaleFactor = 1.0_DP
        rmatrix%RmatrixBlock(2,2)%dscaleFactor = 1.0_DP
        
        ! Determine how to assemble...
        ccontrolConstraints = rnonlinearSpatialMatrix%rdiscrData%rconstraints%ccontrolConstraints
          
        ! Calculate the usual mass matrix if conrol constraints are deactivated
        ! or if Newton is not active.
        if ((ccontrolConstraints .eq. 0) .or. &
            (rnonlinearSpatialMatrix%cmatrixType .eq. 0)) then
      
          ! Copy the entries of the mass matrix. Share the structure.
          ! We must not share the entries as these might be changed by the caller
          ! e.g. due to boundary conditions!
          
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
              
          ! Scale the entries by dweight.
          if (dweight .ne. 1.0_DP) then
            call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(1,1),dweight)
            call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(2,2),dweight)
          end if
          
        else if (rnonlinearSpatialMatrix%cmatrixType .eq. 1) then
          
          ! Constraints implemented by filtered mass matrices.
          
          select case (rnonlinearSpatialMatrix%rdiscrData%rconstraints%ccontrolConstraints)
          
          case (1)
          
            ! Calculate the operator. This automatically creates new matrices
            ! in the diagonal blocks of rmatrix, since the diagonal blocks
            ! are empty at the moment.
            !
            ! The weight "-rnonlinearSpatialMatrix%dalphaC" is included in
            ! dweight, so it has to be canceled out.

            select case (rnonlinearSpatialMatrix%rdiscrData%rconstraints%cconstraintsType)
            case (0)
              ! Constant bounds
            
              call nwder_minMaxProjByMass (&
                  rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
                  -rnonlinearSpatialMatrix%dalphaC*dweight,rmatrix%RmatrixBlock(1,1),&
                  .true.,&
                  -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvector%RvectorBlock(4),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1)

              call nwder_minMaxProjByMass (&
                  rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
                  -rnonlinearSpatialMatrix%dalphaC*dweight,rmatrix%RmatrixBlock(2,2),&
                  .true.,&
                  -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvector%RvectorBlock(5),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2)

            case (1)
              ! Variable bounds

              call nwder_minMaxProjByMass (&
                  rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
                  -rnonlinearSpatialMatrix%dalphaC*dweight,rmatrix%RmatrixBlock(1,1),&
                  .true.,&
                  -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvector%RvectorBlock(4),&
                  1.0_DP,1.0_DP,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(1),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(1))

              call nwder_minMaxProjByMass (&
                  rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
                  -rnonlinearSpatialMatrix%dalphaC*dweight,rmatrix%RmatrixBlock(2,2),&
                  .true.,&
                  -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvector%RvectorBlock(5),&
                  1.0_DP,1.0_DP,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(2),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(2))

            case default
              ! Not implemented.
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
          
          case (2)
          
            ! Exact reassembly of the mass matrices.
          
            ! In A11/A22 we have to create a 'projective mass matrix'.
            ! This is the derivative of a projection operator
            ! P[a,b](f)=a if f<a, =b if f>b, =f otherwise.
            ! For a<f<b, this is the mass matrix. Everywhere else, this is =0.
            ! We assemble this matrix just as a standard mass matrix with noconstant
            ! coefficients. Whereever u = -1/alpha * lambda is out of bounds,
            ! we return 0 as coefficient, otherwise 1.
            !
            ! For alpha=0, we have bang-bang-control. In this case,
            ! the Newton derivative is the zero operator, so nothing
            ! has to be assembled here.
            !
            ! Alpha < 0 switches the distributed control off, so also in this
            ! case, nothing has to be assembled.
            if (rnonlinearSpatialMatrix%dalphaC .gt. 0.0_DP) then

              ! Create a cubature info structure that defines the cubature
              ! rule to use.            
              call spdiscr_createDefCubStructure (&
                  rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimalDual%RspatialDiscr(1),&
                  rcubatureInfo,rnonlinearSpatialMatrix%rdiscrData%rsettingsSpaceDiscr%icubMass)

              call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,1))
              call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,2))

              ! Calculate the Newton derivative.
              
              select case (rnonlinearSpatialMatrix%rdiscrData%rconstraints%cconstraintsType)
              case (0)
                ! Cónstant bounds

                ! The weight contains -1/alpha, so we have to cancel it out here.
                call nwder_minMaxProjByCubature (-rnonlinearSpatialMatrix%dalphaC*dweight,&
                    rmatrix%RmatrixBlock(1,1),rcubatureInfo,&
                    -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvector%RvectorBlock(4),&
                    rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1,&
                    rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1)

                call nwder_minMaxProjByCubature (-rnonlinearSpatialMatrix%dalphaC*dweight,&
                    rmatrix%RmatrixBlock(2,2),rcubatureInfo,&
                    -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvector%RvectorBlock(5),&
                    rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2,&
                    rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2)
                    
              case (1)
                ! Variable bounds
                
                ! The weight contains -1/alpha, so we have to cancel it out here.
                call nwder_minMaxProjByCubature (-rnonlinearSpatialMatrix%dalphaC*dweight,&
                    rmatrix%RmatrixBlock(1,1),rcubatureInfo,&
                    -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvector%RvectorBlock(4),&
                    1.0_DP,1.0_DP,&
                    rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(1),&
                    rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(1))

                call nwder_minMaxProjByCubature (-rnonlinearSpatialMatrix%dalphaC*dweight,&
                    rmatrix%RmatrixBlock(2,2),rcubatureInfo,&
                    -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvector%RvectorBlock(5),&
                    1.0_DP,1.0_DP,&
                    rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(2),&
                    rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(2))
              end select
            
              ! Release the cubature structure
              call spdiscr_releaseCubStructure(rcubatureInfo)
            
            end if
            
          case (3)
          
            call lsyssc_duplicateMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
                rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          
            call lsyssc_duplicateMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
                rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
                
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,1))
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,2))

            select case (rnonlinearSpatialMatrix%rdiscrData%rconstraints%cconstraintsType)
            case (0)
              ! Constant bounds

              call nwder_minMaxProjByApproxDer (-rnonlinearSpatialMatrix%dalphaC*dweight,&
                  rmatrix%RmatrixBlock(1,1),0.001_DP,&
                  -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvector%RvectorBlock(4),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1)

              call nwder_minMaxProjByApproxDer (-rnonlinearSpatialMatrix%dalphaC*dweight,&
                  rmatrix%RmatrixBlock(2,2),0.001_DP,&
                  -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvector%RvectorBlock(5),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2)

            case (1)
              ! Variable bounds

              call nwder_minMaxProjByApproxDer (-rnonlinearSpatialMatrix%dalphaC*dweight,&
                  rmatrix%RmatrixBlock(1,1),0.001_DP,&
                  -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvector%RvectorBlock(4),&
                  1.0_DP,1.0_DP,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(1),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(1))

              call nwder_minMaxProjByApproxDer (-rnonlinearSpatialMatrix%dalphaC*dweight,&
                  rmatrix%RmatrixBlock(2,2),0.001_DP,&
                  -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvector%RvectorBlock(5),&
                  1.0_DP,1.0_DP,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(2),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(2))
            
            end select

          case (4)
          
            ! Exact reassembly of the mass matrices with adaptive integration.
            
            ! In A11/A22 we have to create a 'projective mass matrix'.
            ! This is the derivative of a projection operator
            ! P[a,b](f)=a if f<a, =b if f>b, =f otherwise.
            ! For a<f<b, this is the mass matrix. Everywhere else, this is =0.
            ! We assemble this matrix just as a standard mass matrix with noconstant
            ! coefficients. Whereever u = -1/alpha * lambda is out of bounds,
            ! we return 0 as coefficient, otherwise 1.
            !
            ! For alpha=0, we have bang-bang-control. In this case,
            ! the Newton derivative is the zero operator, so nothing
            ! has to be assembled here.
            !
            ! Alpha < 0 switches the distributed control off, so also in this
            ! case, nothing has to be assembled.
            if (rnonlinearSpatialMatrix%dalphaC .gt. 0.0_DP) then

              ! Create a cubature info structure that defines the cubature
              ! rule to use.            
              call spdiscr_createDefCubStructure (&
                  rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimalDual%RspatialDiscr(1),&
                  rcubatureInfo,rnonlinearSpatialMatrix%rdiscrData%rsettingsSpaceDiscr%icubMass)
                  
              ! And another one for the adaptive cubature rule on the border
              ! of the active set.
              call spdiscr_createDefCubStructure (&
                  rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimalDual%RspatialDiscr(1),&
                  rcubatureInfoAdapt,cub_getSummedCubType(&
                      rnonlinearSpatialMatrix%rdiscrData%rsettingsSpaceDiscr%icubMass,4))

              call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,1))
              call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,2))

              ! Calculate the Newton derivative.
              
              select case (rnonlinearSpatialMatrix%rdiscrData%rconstraints%cconstraintsType)
              case (0)
                ! Cónstant bounds

                ! The weight contains -1/alpha, so we have to cancel it out here.
                call nwder_minMaxProjByAdaptCub (-rnonlinearSpatialMatrix%dalphaC*dweight,&
                    rmatrix%RmatrixBlock(1,1),rcubatureInfo,rcubatureInfoAdapt,&
                    -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvector%RvectorBlock(4),&
                    rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1,&
                    rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1)

                call nwder_minMaxProjByAdaptCub (-rnonlinearSpatialMatrix%dalphaC*dweight,&
                    rmatrix%RmatrixBlock(2,2),rcubatureInfo,rcubatureInfoAdapt,&
                    -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvector%RvectorBlock(5),&
                    rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2,&
                    rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2)
                    
              case (1)
                ! Variable bounds
                
                ! The weight contains -1/alpha, so we have to cancel it out here.
                call nwder_minMaxProjByAdaptCub (-rnonlinearSpatialMatrix%dalphaC*dweight,&
                    rmatrix%RmatrixBlock(1,1),rcubatureInfo,rcubatureInfoAdapt,&
                    -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvector%RvectorBlock(4),&
                    1.0_DP,1.0_DP,&
                    rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(1),&
                    rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(1))

                call nwder_minMaxProjByAdaptCub (-rnonlinearSpatialMatrix%dalphaC*dweight,&
                    rmatrix%RmatrixBlock(2,2),rcubatureInfo,rcubatureInfoAdapt,&
                    -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvector%RvectorBlock(5),&
                    1.0_DP,1.0_DP,&
                    rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(2),&
                    rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(2))
              end select
            
              ! Release the cubature structure
              call spdiscr_releaseCubStructure(rcubatureInfoAdapt)
              call spdiscr_releaseCubStructure(rcubatureInfo)
            
            end if
            
          case default
          
            call output_line ('Unsupported ccontrolConstraints flag.', &
                              OU_CLASS_ERROR,OU_MODE_STD,'assembleProjectedMassBlocks')
            call sys_halt()
            
          end select

        else
        
          call output_line ('Unsupported cmatrixType.', &
                            OU_CLASS_ERROR,OU_MODE_STD,'assembleProjectedMassBlocks')
          call sys_halt()
          
        end if
        
      else
      
        ! Switch off these matrices
        rmatrix%RmatrixBlock(1,1)%dscaleFactor = 0.0_DP
        rmatrix%RmatrixBlock(2,2)%dscaleFactor = 0.0_DP
        
      end if
            
    end subroutine

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine smva_assembleDefect (rnonlinearSpatialMatrix,rx,rd,cx)

!<description>
  ! This routine assembles the nonlinear defect
  !      rd := rd - cx A(ry) rx
  ! with the system matrix A(.) defined by the configuration in rnonlinearSpatialMatrix.
  ! The caller must initialise the rnonlinearSpatialMatrix according to how the
  ! matrix should look like.
  !
  ! The parameters rvectorI are optional. If specified, these parameter defines where to
  ! evaluate the nonlinearity (if the system matrix $A$ contains a nonlinearity).
  ! The multiplication factor cx is optional as well. If not specified,
  ! cx=1.0 is assumed.
  !
  ! Note that for evaluating the nonlinearity / Newton matrix in
  ! the dual equation, the primal velocity is used!
  ! This is due to the fact that the 'nonlinearity' in the dual equation
  ! is defined as
  !       $$ y \grad(\lambda)  +   \lambda \grad(y)  $$
  ! so both parts are evaluated at the point y of the primal velocity
  ! before they are multiplied by the dual velocity $\lambda$!
  !
  ! The routine will not include any boundary conditions in the defect.
!</description>

  ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
  ! about how to set up the matrix.
  !
  ! The caller must provide either p_rmatrixTemplateXXXX in this structure
  ! or set the p_rmatrixTemplateXXXX as well as p_rdiscretisation to
  ! appropriate values. This is necessary for exploiting then structure
  ! of the matrix.
  type(t_nonlinearSpatialMatrix), intent(IN) :: rnonlinearSpatialMatrix

  ! This vector specifies the 'x' that is multiplied to the matrix.
  type(t_vectorBlock), intent(IN), target :: rx

  ! OPTIONAL: Multiplication factor in front of the term 'A(ry) rx'.
  ! If not specified, cx=1.0 is assumed.
  real(DP), intent(IN), optional :: cx

!</input>

!<inputoutput>
  ! Destination vector. cx*A(ry)*rx is subtracted from this vector.
  type(t_vectorBlock), intent(INOUT) :: rd
!</inputoutput>
  
!</subroutine>

    ! local variables
    real(DP) :: dcx
    type(t_matrixBlock) :: rtempmatrix
    type(t_matrixScalar) :: rtempMassMat
    type(t_vectorBlock) :: rtempVectorX,rtempVectorB
    type(t_vectorBlock) :: rvectorPrimal,rvectorDual,rvectorDual2
    type(t_vectorBlock), pointer :: p_rprimalSol, p_rdualSol
    type(t_settings_stabil) :: rstabilisation
    type(t_optcoperator) :: roptcoperator
    type(t_blockDiscretisation) :: rvelDiscr
    real(DP) :: dweightConvection,dweightDualConvection,dweightNaturalBdcDual
    real(DP) :: dweightDualNewtonT
    type(t_scalarCubatureInfo) :: rcubatureInfo
    integer :: ndofs
    integer, dimension(:), pointer :: p_Idofs
    
    logical, parameter :: bnewmethod = .false.
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Dd,p_Dx
    
    call lsysbl_getbase_double (rd,p_Dd)
    call lsysbl_getbase_double (rx,p_Dx)
    
    ! Debug weights for the convection
    dweightConvection = &
        rnonlinearSpatialMatrix%rdiscrData%p_rdebugFlags%dweightConvection
    dweightDualConvection = &
        rnonlinearSpatialMatrix%rdiscrData%p_rdebugFlags%dweightDualConvection
    dweightNaturalBdcDual = &
        rnonlinearSpatialMatrix%rdiscrData%p_rdebugFlags%dweightNaturalBdcDual
    dweightDualNewtonT = &
        rnonlinearSpatialMatrix%rdiscrData%p_rdebugFlags%dweightDualNewtonT
    
    dcx = 1.0_DP
    if (present(cx)) dcx = cx
    if (dcx .eq. 0.0_DP) return
    
    if (.not. bnewmethod) then
    
      ! The system matrix looks like:
      !
      !    ( A11  A12  B1  M~           )
      !    ( A21  A22  B2       M~      )
      !    ( B1^T B2^T                  )
      !    ( M             A44  A45  B1 )
      !    (      M        A54  A55  B2 )
      !    (               B1^T B2^T    )
      !
      ! We have now to manually perform a matrix-vector multiplication with
      ! all submatrices in that matrix. At first, we start with the
      ! linear part, that's the easiest.
      !
      ! ---------------------------------------------------
      ! 1.) Identity.
      !    ( I                          )
      !    (      I                     )
      !    (                            )
      !    (               I            )
      !    (                    I       )
      !    (                            )

      if (rnonlinearSpatialMatrix%DidentityY(1,1) .ne. 0.0_DP) then
        call lsyssc_vectorLinearComb (&
            rx%RvectorBlock(1), rd%RvectorBlock(1), &
            -rnonlinearSpatialMatrix%DidentityY(1,1)*dcx, 1.0_DP)

        call lsyssc_vectorLinearComb (&
            rx%RvectorBlock(2), rd%RvectorBlock(2), &
            -rnonlinearSpatialMatrix%DidentityY(1,1)*dcx, 1.0_DP)
      end if

      if (rnonlinearSpatialMatrix%DidentityY(2,2) .ne. 0.0_DP) then
        call lsyssc_vectorLinearComb (&
            rx%RvectorBlock(4), rd%RvectorBlock(4), &
            -rnonlinearSpatialMatrix%DidentityY(2,2)*dcx, 1.0_DP)

        call lsyssc_vectorLinearComb (&
            rx%RvectorBlock(5), rd%RvectorBlock(5), &
            -rnonlinearSpatialMatrix%DidentityY(2,2)*dcx, 1.0_DP)
      end if

      ! ---------------------------------------------------
      ! 2.) Mass matrices
      !    ( M                          )
      !    (      M                     )
      !    (                            )
      !    ( M             M            )
      !    (      M             M       )
      !    (                            )
      if (rnonlinearSpatialMatrix%Dmass(1,1) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity, &
            rx%RvectorBlock(1), rd%RvectorBlock(1), &
            -rnonlinearSpatialMatrix%Dmass(1,1)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity, &
            rx%RvectorBlock(2), rd%RvectorBlock(2), &
            -rnonlinearSpatialMatrix%Dmass(1,1)*dcx, 1.0_DP)
      end if

      if (rnonlinearSpatialMatrix%Dmass(2,2) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity, &
            rx%RvectorBlock(4), rd%RvectorBlock(4), &
            -rnonlinearSpatialMatrix%Dmass(2,2)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity, &
            rx%RvectorBlock(5), rd%RvectorBlock(5), &
            -rnonlinearSpatialMatrix%Dmass(2,2)*dcx, 1.0_DP)
      end if

      if (rnonlinearSpatialMatrix%Dmass(2,1) .ne. 0.0_DP) then
        if (.not. associated(rnonlinearSpatialMatrix%rdiscrData%p_DobservationArea)) then
          
          ! No specific observation area.
          ! Observe the complete domain.
          call lsyssc_scalarMatVec (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity, &
              rx%RvectorBlock(1), rd%RvectorBlock(4), &
              -rnonlinearSpatialMatrix%Dmass(2,1)*dcx, 1.0_DP)

          call lsyssc_scalarMatVec (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity, &
              rx%RvectorBlock(2), rd%RvectorBlock(5), &
              -rnonlinearSpatialMatrix%Dmass(2,1)*dcx, 1.0_DP)
              
        else
          ! A box as observation area. Only assemble the defect in
          ! this box. For that purpose, assemble a mass matrix corresponding
          ! to the characteristic function of the box.
          call lsyssc_duplicateMatrix (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
            rtempMassMat,LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          call lsyssc_clearMatrix (rtempMassMat)

          call spdiscr_createDefCubStructure (&
              rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimalDual%RspatialDiscr(1),&
              rcubatureInfo, rnonlinearSpatialMatrix%rdiscrData%rsettingsSpaceDiscr%icubMass)
          
          call assembleBoxIdentity (rnonlinearSpatialMatrix%rdiscrData%p_DobservationArea,&
              1.0_DP,rtempMassMat,rcubatureInfo)

          call spdiscr_releaseCubStructure (rcubatureInfo)
          
          ! Now apply the defect with this matrix.
          call lsyssc_scalarMatVec (rtempMassMat, &
              rx%RvectorBlock(1), rd%RvectorBlock(4), &
              -rnonlinearSpatialMatrix%Dmass(2,1)*dcx, 1.0_DP)

          call lsyssc_scalarMatVec (rtempMassMat, &
              rx%RvectorBlock(2), rd%RvectorBlock(5), &
              -rnonlinearSpatialMatrix%Dmass(2,1)*dcx, 1.0_DP)

          ! Release memory
          call lsyssc_releaseMatrix (rtempMassMat)
        end if
      end if
      
      ! Don't do anything with Dmass(1,2) -- the mass matrices here
      ! are probably nonlinear! We assemble them later!
      
      ! ---------------------------------------------------
      ! 3.) Stokes matrices
      !    ( L                          )
      !    (      L                     )
      !    (                            )
      !    (               L            )
      !    (                    L       )
      !    (                            )
      if (rnonlinearSpatialMatrix%Dstokes(1,1) .ne. 0.0_DP) then
        ! In case of the deformation tensor or nonconstant viscosity,
        ! the defect assembly is done during the assembly of the nonlinearity.
        if ((rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%isubequation .eq. 0) .and. &
            (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%cviscoModel .eq. 0)) then
          call lsyssc_scalarMatVec (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixLaplace, &
              rx%RvectorBlock(1), rd%RvectorBlock(1), &
              -rnonlinearSpatialMatrix%Dstokes(1,1)*dcx*&
                rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst,&
              1.0_DP)

          call lsyssc_scalarMatVec (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixLaplace, &
              rx%RvectorBlock(2), rd%RvectorBlock(2), &
              -rnonlinearSpatialMatrix%Dstokes(1,1)*dcx*&
                rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst,&
              1.0_DP)
        end if
      end if
            
      if (rnonlinearSpatialMatrix%Dstokes(2,2) .ne. 0.0_DP) then
        ! In case of the deformation tensor or nonconstant viscosity,
        ! the defect assembly is done during the assembly of the nonlinearity.
        if ((rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%isubequation .eq. 0) .and. &
            (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%cviscoModel .eq. 0)) then
          call lsyssc_scalarMatVec (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixLaplace, &
              rx%RvectorBlock(4), rd%RvectorBlock(4), &
              -rnonlinearSpatialMatrix%Dstokes(2,2)*dcx*&
                rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst,&
              1.0_DP)

          call lsyssc_scalarMatVec (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixLaplace, &
              rx%RvectorBlock(5), rd%RvectorBlock(5), &
              -rnonlinearSpatialMatrix%Dstokes(2,2)*dcx*&
                rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst,&
              1.0_DP)
        end if
      end if
      
      ! ---------------------------------------------------
      ! 3.) B-matrices
      !    (           B1               )
      !    (           B2               )
      !    (                            )
      !    (                        B1  )
      !    (                        B2  )
      !    (                            )
      
      if (rnonlinearSpatialMatrix%DBmat(1,1) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
            rx%RvectorBlock(3), rd%RvectorBlock(1), &
            -rnonlinearSpatialMatrix%DBmat(1,1)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
            rx%RvectorBlock(3), rd%RvectorBlock(2), &
            -rnonlinearSpatialMatrix%DBmat(1,1)*dcx, 1.0_DP)
      end if
      
      if (rnonlinearSpatialMatrix%DBmat(2,2) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
            rx%RvectorBlock(6), rd%RvectorBlock(4), &
            -rnonlinearSpatialMatrix%DBmat(2,2)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
            rx%RvectorBlock(6), rd%RvectorBlock(5), &
            -rnonlinearSpatialMatrix%DBmat(2,2)*dcx, 1.0_DP)
      end if
      
      ! ---------------------------------------------------
      ! 4.) B^T-matrices
      !    (                            )
      !    (                            )
      !    ( B1^T B2^T                  )
      !    (                            )
      !    (                            )
      !    (              B1^T B2^T     )
      
      if (rnonlinearSpatialMatrix%DBTmat(1,1) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD1, &
            rx%RvectorBlock(1), rd%RvectorBlock(3), &
            -rnonlinearSpatialMatrix%DBTmat(1,1)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD2, &
            rx%RvectorBlock(2), rd%RvectorBlock(3), &
            -rnonlinearSpatialMatrix%DBTmat(1,1)*dcx, 1.0_DP)
      end if
      
      if (rnonlinearSpatialMatrix%DBTmat(2,2) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD1, &
            rx%RvectorBlock(4), rd%RvectorBlock(6), &
            -rnonlinearSpatialMatrix%DBTmat(2,2)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD2, &
            rx%RvectorBlock(5), rd%RvectorBlock(6), &
            -rnonlinearSpatialMatrix%DBTmat(2,2)*dcx, 1.0_DP)
      end if
      
      ! ---------------------------------------------------
      ! Now a slightly more advanced task for which we use a separate
      ! routine and some submatrices/vectors: The nonlinearity.
      !
      ! Get a reference to the correct velocity vector we have to use for the
      ! nonlinearity.
      select case (rnonlinearSpatialMatrix%iprimalSol)
      case (1)
        call lsysbl_deriveSubvector(&
            rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1,rvectorPrimal, 1,2,.true.)
      case (2)
        call lsysbl_deriveSubvector(&
          rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2,rvectorPrimal, 1,2,.true.)
      case (3)
        call lsysbl_deriveSubvector(&
            rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3,rvectorPrimal, 1,2,.true.)
      end select

      select case (rnonlinearSpatialMatrix%idualSol)
      case (1)
        call lsysbl_deriveSubvector(&
            rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1,rvectorDual, 4,5,.true.)
      case (2)
        call lsysbl_deriveSubvector(&
            rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2,rvectorDual, 4,5,.true.)
      case (3)
        call lsysbl_deriveSubvector(&
            rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3,rvectorDual, 4,5,.true.)
      end select

      select case (rnonlinearSpatialMatrix%idualSol2)
      case (1)
        call lsysbl_deriveSubvector(&
            rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1,rvectorDual2, 4,5,.true.)
      case (2)
        call lsysbl_deriveSubvector(&
            rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2,rvectorDual2, 4,5,.true.)
      case (3)
        call lsysbl_deriveSubvector(&
            rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3,rvectorDual2, 4,5,.true.)
      end select

      ! Create a block discretisation by deriving it from the 'full' matrix.
      ! This will serve as a local discretisation structure for all
      ! velocity modifications.
      call spdiscr_deriveBlockDiscr (rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimalDual, &
          rvelDiscr, 1, 2)

      ! Create a 2x2 matrix based on the structure of the FE space.
      ! The matrix does not need any entries, we only need the structure.
      call lsysbl_createMatBlockByDiscr (rvelDiscr,rtempMatrix)
      call lsyssc_duplicateMatrix (&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEM,&
          rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_duplicateMatrix (&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEMOffdiag,&
          rtempMatrix%RmatrixBlock(1,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_duplicateMatrix (&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEMOffdiag,&
          rtempMatrix%RmatrixBlock(2,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_duplicateMatrix (&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEM,&
          rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

      ! 1.) Primal equation, y*grad(.), probably + grad(.)*y
      !
      !    ( N(y) N(y)                  )
      !    ( N(y) N(y)                  )
      !    (                            )
      !    (                            )
      !    (                            )
      !    (                            )
      
      ! For the nonlinearity, apply cubature formula icubStokes
      call spdiscr_createDefCubStructure (&
          rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimalDual%RspatialDiscr(1),&
          rcubatureInfo, rnonlinearSpatialMatrix%rdiscrData%rsettingsSpaceDiscr%icubStokes)
      
      call lsysbl_deriveSubvector(rx,rtempVectorX,1,2,.true.)
      call lsysbl_deriveSubvector(rd,rtempVectorB,1,2,.true.)

      call assembleNonlinearityDefect (&
          rnonlinearSpatialMatrix,rtempMatrix,&
          rvectorPrimal,rvectorPrimal,rtempVectorX,rtempVectorB,&
          rnonlinearSpatialMatrix%Dstokes(1,1),&
          dweightConvection*rnonlinearSpatialMatrix%Dygrad(1,1),&
          dweightConvection*rnonlinearSpatialMatrix%DygradAdj(1,1),&
          dweightConvection*rnonlinearSpatialMatrix%DygradT(1,1),&
          dweightConvection*rnonlinearSpatialMatrix%DygradTAdj(1,1),&
          dweightConvection*rnonlinearSpatialMatrix%Dgrady(1,1),&
          dweightConvection*rnonlinearSpatialMatrix%DgradyT(1,1),&
          dweightConvection*rnonlinearSpatialMatrix%DgradyAdj(1,1),&
          rcubatureInfo,rnonlinearSpatialMatrix%rdiscrData%rstabilPrimal,dcx,&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplatesOptC%rmatrixEOJ1)
      
      call lsysbl_releaseVector (rtempVectorX)
      call lsysbl_releaseVector (rtempVectorB)
      
      ! 2.) Dual equation, y*grad(.) + grad(.)^T*y
      !
      !    (                            )
      !    (                            )
      !    (                            )
      !    (              N*(y) N*(y)   )
      !    (              N*(y) N*(y)   )
      !    (                            )
      
      call lsysbl_deriveSubvector(rx,rtempVectorX,4,5,.true.)
      call lsysbl_deriveSubvector(rd,rtempVectorB,4,5,.true.)

      call assembleNonlinearityDefect (&
          rnonlinearSpatialMatrix,rtempMatrix,&
          rvectorPrimal,rvectorPrimal,rtempVectorX,rtempVectorB,&
          rnonlinearSpatialMatrix%Dstokes(2,2),&
          dweightDualConvection*rnonlinearSpatialMatrix%Dygrad(2,2),&
          dweightDualConvection*rnonlinearSpatialMatrix%DygradAdj(2,2),&
          dweightConvection*rnonlinearSpatialMatrix%DygradT(2,2),&
          dweightConvection*rnonlinearSpatialMatrix%DygradTAdj(2,2),&
          dweightConvection*rnonlinearSpatialMatrix%Dgrady(2,2),&
          dweightDualNewtonT*rnonlinearSpatialMatrix%DgradyT(2,2),&
          dweightConvection*rnonlinearSpatialMatrix%DgradyAdj(2,2),&
          rcubatureInfo,rnonlinearSpatialMatrix%rdiscrData%rstabilDual,dcx,&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplatesOptC%rmatrixEOJ2)
      
      call lsysbl_releaseVector (rtempVectorX)
      call lsysbl_releaseVector (rtempVectorB)
      
      ! 3.) Dual equation, Newton term lambda*grad(.) + grad(.)^T*lambda
      !
      !    (                            )
      !    (                            )
      !    (                            )
      !    (   N*(l) N*(l)              )
      !    (   N*(l) N*(l)              )
      !    (                            )
      
      ! No stabilisation here
      ! rstabilisation = t_convecStabilisation(&
      !     rnonlinearSpatialMatrix%iupwind2,rnonlinearSpatialMatrix%dupsam2)
      rstabilisation = t_settings_stabil(3,0.0_DP,1,1,1)
      
      call lsysbl_deriveSubvector(rx,rtempVectorX,1,2,.true.)
      call lsysbl_deriveSubvector(rd,rtempVectorB,4,5,.true.)

      call assembleNonlinearityDefect (&
          rnonlinearSpatialMatrix,rtempMatrix,&
          rvectorPrimal,rvectorDual,rtempVectorX,rtempVectorB,&
          rnonlinearSpatialMatrix%Dstokes(2,1),&
          dweightConvection*rnonlinearSpatialMatrix%Dygrad(2,1),&
          dweightConvection*rnonlinearSpatialMatrix%DygradAdj(2,1),&
          dweightConvection*rnonlinearSpatialMatrix%DygradT(2,1),&
          dweightConvection*rnonlinearSpatialMatrix%DygradTAdj(2,1),&
          dweightConvection*rnonlinearSpatialMatrix%Dgrady(2,1),&
          dweightConvection*rnonlinearSpatialMatrix%DgradyT(2,1),&
          dweightConvection*rnonlinearSpatialMatrix%DgradyAdj(2,1),&
          rcubatureInfo,rstabilisation,dcx)
          
      ! There is probably a 2nd reactive term involved stemming from
      ! the next timestep when Crank-Nicolson is used.

      rstabilisation = t_settings_stabil(3,0.0_DP,1,1,1)
      
      call assembleNonlinearityDefect (&
          rnonlinearSpatialMatrix,rtempMatrix,&
          rvectorPrimal,rvectorDual2,rtempVectorX,rtempVectorB,0.0_DP,0.0_DP,0.0_DP,&
          dweightConvection*rnonlinearSpatialMatrix%DygradT2(2,1),&
          dweightConvection*rnonlinearSpatialMatrix%DygradTAdj2(2,1),&
          dweightConvection*rnonlinearSpatialMatrix%Dgrady2(2,1),&
          0.0_DP,dweightConvection*rnonlinearSpatialMatrix%DgradyAdj2(2,1),&
          rcubatureInfo,rstabilisation,dcx)
      
      call lsysbl_releaseVector (rtempVectorX)
      call lsysbl_releaseVector (rtempVectorB)
      
      ! Release cubature-related stuff
      call spdiscr_releaseCubStructure (rcubatureInfo)
      
      ! 4.) (Projected) mass matrix
      !
      !    (              PM(l)         )
      !    (                    PM(l)   )
      !    (                            )
      !    (                            )
      !    (                            )
      !    (                            )
      
      ! What's the type of the current matrix? Is this a Newton-matrix?
      if (rnonlinearSpatialMatrix%cmatrixType .eq. 0) then

        ! No, this is a standard matrix. That means, we just have to project
        ! the control u.
        if (rnonlinearSpatialMatrix%Dmass(1,2) .ne. 0.0_DP) then

          ! Apply the projection operator to (-1/alpha lambda) and sum up the 
          ! resulting vector to the current defect.
          select case (rnonlinearSpatialMatrix%rdiscrData%rconstraints%ccontrolConstraints)
          
          case (0)
            ! No constraints.
            !
            ! Multiply with the mass matrix to include it to the defect.
            ! Note that the multiplication factor is -(-cx) = cx because
            ! it's put on the RHS of the system for creating the defect.
            ! d = b - cx A x = b - ... + \nu Laplace(y) - y\grad(y) - grad(p) + P(-1/alpha lambda)
            call lsyssc_scalarMatVec (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity, &
                rx%RvectorBlock(4), rd%RvectorBlock(1), &
                -rnonlinearSpatialMatrix%Dmass(1,2)*dcx, 1.0_DP)
            call lsyssc_scalarMatVec (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity, &
                rx%RvectorBlock(5), rd%RvectorBlock(2), &
                -rnonlinearSpatialMatrix%Dmass(1,2)*dcx, 1.0_DP)
          
          case (2)
            ! Use the dual solution as right hand side and assemble a temporary vector
            ! like a right hand side. The result can be added to the defect.
            
            ! Apply cubature formula icubF
            call spdiscr_createDefCubStructure (&
                rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimalDual%RspatialDiscr(1),&
                rcubatureInfo, rnonlinearSpatialMatrix%rdiscrData%rsettingsSpaceDiscr%icubF)
            
            ! Apply: rd(primal) = rd(primal) + dcx * P(-1/alpha * lambda)
            
            select case (rnonlinearSpatialMatrix%rdiscrData%rconstraints%cconstraintsType)
            case (0)
              ! Constant bounds
              call nwder_rhsMinMaxProjByCubature (dcx,rd%RvectorBlock(1),&
                  rcubatureInfo,-rnonlinearSpatialMatrix%Dmass(1,2),rx%RvectorBlock(4),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1)

              call nwder_rhsMinMaxProjByCubature (dcx,rd%RvectorBlock(2),&
                  rcubatureInfo,-rnonlinearSpatialMatrix%Dmass(1,2),rx%RvectorBlock(5),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2)
                  
            case (1)
              ! Variable bounds
              call nwder_rhsMinMaxProjByCubature (dcx,rd%RvectorBlock(1),&
                  rcubatureInfo,-rnonlinearSpatialMatrix%Dmass(1,2),rx%RvectorBlock(4),&
                  1.0_DP,1.0_DP,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(1),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(1))

              call nwder_rhsMinMaxProjByCubature (dcx,rd%RvectorBlock(2),&
                  rcubatureInfo,-rnonlinearSpatialMatrix%Dmass(1,2),rx%RvectorBlock(5),&
                  1.0_DP,1.0_DP,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(2),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(2))
              
            case default
              ! Not implemented.
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
            
            call spdiscr_releaseCubStructure (rcubatureInfo)
            
          case (4)
            ! Not implemented.
            !
            ! Must be done similar to (2), i.e., a RHS must be calculated.
            ! However, on the cells crossing the active set, the RHS must
            ! be calculated with a different (summed) cubature formula.
            call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
            call sys_halt()
        
          case default
          
            ! Temporary vector.
            call lsysbl_deriveSubvector(rx,rtempVectorX,4,4,.false.)
          
            ! Just project the DOF's.
            
            select case (rnonlinearSpatialMatrix%rdiscrData%rconstraints%cconstraintsType)
            case (0)
              ! Constant bounds
              call nwder_rhsMinMaxProjByMass (dcx,rd%RvectorBlock(1),&
                  rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
                  rtempVectorX%RvectorBlock(1),&
                  -rnonlinearSpatialMatrix%Dmass(1,2),rx%RvectorBlock(4),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1)

              call nwder_rhsMinMaxProjByMass (dcx,rd%RvectorBlock(2),&
                  rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
                  rtempVectorX%RvectorBlock(1),&
                  -rnonlinearSpatialMatrix%Dmass(1,2),rx%RvectorBlock(5),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2)
                  
            case (1)
              ! Variable bounds
              call nwder_rhsMinMaxProjByMass (dcx,rd%RvectorBlock(1),&
                  rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
                  rtempVectorX%RvectorBlock(1),&
                  -rnonlinearSpatialMatrix%Dmass(1,2),rx%RvectorBlock(4),&
                  1.0_DP,1.0_DP,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(1),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(1))

              call nwder_rhsMinMaxProjByMass (dcx,rd%RvectorBlock(2),&
                  rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
                  rtempVectorX%RvectorBlock(1),&
                  -rnonlinearSpatialMatrix%Dmass(1,2),rx%RvectorBlock(5),&
                  1.0_DP,1.0_DP,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(2),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(2))
              
            case default
              ! Not implemented.
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
            
            ! Release temp memory
            call lsysbl_releaseVector (rtempVectorX)
            
          end select
               
        end if

      else

        if (rnonlinearSpatialMatrix%Dmass(1,2) .ne. 0.0_DP) then
          ! Yes, that's a Newton matrix. That means, we have to multiply the
          ! vector with the derivative of the projection operator:
          ! b-(-P[a,b]'(-1/alpha lambda)).
          ! For that purpose, we have to assemble special mass matrices:
          select case (rnonlinearSpatialMatrix%idualSol)
          case (1)
            call assemblePrimalUConstrMassDefect (&
                rnonlinearSpatialMatrix,rx,&
                rd,dcx*rnonlinearSpatialMatrix%Dmass(1,2),&
                rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1)
          case (2)
            call assemblePrimalUConstrMassDefect (&
                rnonlinearSpatialMatrix,rx,&
                rd,dcx*rnonlinearSpatialMatrix%Dmass(1,2),&
                rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2)
          case (3)
            call assemblePrimalUConstrMassDefect (&
                rnonlinearSpatialMatrix,rx,&
                rd,dcx*rnonlinearSpatialMatrix%Dmass(1,2),&
                rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3)
          end select
        end if
        
      end if
      
      ! Release the temp vectors/matrices, that's it.
      call spdiscr_releaseBlockDiscr(rvelDiscr)
      call lsysbl_releaseVector (rvectorDual2)
      call lsysbl_releaseVector (rvectorDual)
      call lsysbl_releaseVector (rvectorPrimal)
      call lsysbl_releaseMatrix (rtempMatrix)

      ! 5.) Special processing: boundary control.
      !
      !    ( QM       v1D       v2M     )
      !    (     QM       v1D   v2M     )
      !    (                            )
      !    (                            )
      !    (                            )
      !    (                            )
      !
      ! with v1. v2 some constants. These matrices are only
      ! applied on the Dirichlet control boundary. The contributions
      ! replace the original lines in the matrices.
      ! Therefore, the corresponding entries in the defect vector
      ! have to be restored and replaced by the new contributions.
      
      ! Create a temporary matrix for the boundary control operator.
      call lsysbl_createMatBlockByDiscr (&
          rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimalDual,rtempMatrix)
      
      call lsyssc_duplicateMatrix (&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEM,&
          rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEM,&
          rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEM,&
          rtempMatrix%RmatrixBlock(1,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEM,&
          rtempMatrix%RmatrixBlock(2,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateGradient,&
          rtempMatrix%RmatrixBlock(1,6),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateGradient,&
          rtempMatrix%RmatrixBlock(2,6),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          
      call lsysbl_clearMatrix (rtempMatrix)
      
      ! Assemble the submatrices
      call smva_assembleDirichletBCC (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal,&
          rtempMatrix%RmatrixBlock(1,1),rtempMatrix%RmatrixBlock(2,2),&
          rtempMatrix%RmatrixBlock(1,4),rtempMatrix%RmatrixBlock(2,5),&
          rtempMatrix%RmatrixBlock(1,6),rtempMatrix%RmatrixBlock(2,6),&
          rnonlinearSpatialMatrix%DdirichletBCCY(1,1),&
          rnonlinearSpatialMatrix%DdirichletBCCLambda(1,2),&
          rnonlinearSpatialMatrix%DdirichletBCCXi(1,2),&
          rnonlinearSpatialMatrix%ddirichletBCPenalty,&
          rnonlinearSpatialMatrix%p_rnonlinearity%p_rdirichletBCCBoundary)
            
        ! Substract
        call lsysbl_blockMatVec (rtempmatrix, rx, rd, -dcx, 1.0_DP)

        ! Release the temp matrix again        
        call lsysbl_releaseMatrix (rtempMatrix)

    else
    
      ! ---------------------------------------------------
      ! 3.) B-matrices
      !    (           B1               )
      !    (           B2               )
      !    (                            )
      !    (                        B1  )
      !    (                        B2  )
      !    (                            )
      
      if (rnonlinearSpatialMatrix%DBmat(1,1) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
            rx%RvectorBlock(3), rd%RvectorBlock(1), &
            -rnonlinearSpatialMatrix%DBmat(1,1)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
            rx%RvectorBlock(3), rd%RvectorBlock(2), &
            -rnonlinearSpatialMatrix%DBmat(1,1)*dcx, 1.0_DP)
      end if
      
      if (rnonlinearSpatialMatrix%DBmat(2,2) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
            rx%RvectorBlock(6), rd%RvectorBlock(4), &
            -rnonlinearSpatialMatrix%DBmat(2,2)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
            rx%RvectorBlock(6), rd%RvectorBlock(5), &
            -rnonlinearSpatialMatrix%DBmat(2,2)*dcx, 1.0_DP)
      end if
      
      ! ---------------------------------------------------
      ! 4.) B^T-matrices
      !    (                            )
      !    (                            )
      !    ( B1^T B2^T                  )
      !    (                            )
      !    (                            )
      !    (              B1^T B2^T     )
      
      if (rnonlinearSpatialMatrix%DBTmat(1,1) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD1, &
            rx%RvectorBlock(1), rd%RvectorBlock(3), &
            -rnonlinearSpatialMatrix%DBTmat(1,1)*dcx, 1.0_DP)
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD2, &
            rx%RvectorBlock(2), rd%RvectorBlock(3), &
            -rnonlinearSpatialMatrix%DBTmat(1,1)*dcx, 1.0_DP)
      end if
      
      if (rnonlinearSpatialMatrix%DBTmat(2,2) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD1, &
            rx%RvectorBlock(4), rd%RvectorBlock(6), &
            -rnonlinearSpatialMatrix%DBTmat(2,2)*dcx, 1.0_DP)
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD2, &
            rx%RvectorBlock(5), rd%RvectorBlock(6), &
            -rnonlinearSpatialMatrix%DBTmat(2,2)*dcx, 1.0_DP)
      end if
      
      ! ---------------------------------------------------
      ! Now a slightly more advanced task for which we use a separate
      ! routine and some submatrices/vectors: The nonlinearity.

      ! Initialise the operator structure for what we need.
      roptcoperator%dupsamPrimal = rnonlinearSpatialMatrix%rdiscrData%rstabilPrimal%dupsam
      roptcoperator%dupsamDual = rnonlinearSpatialMatrix%rdiscrData%rstabilDual%dupsam
      
      ! Timestep-weights
      roptcoperator%dprimalAlpha = rnonlinearSpatialMatrix%Dmass(1,1)
      roptcoperator%ddualAlpha   = rnonlinearSpatialMatrix%Dmass(2,2)

      ! Stokes operator
      roptcoperator%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst
      roptcoperator%dprimalBeta = rnonlinearSpatialMatrix%Dstokes(1,1)
      roptcoperator%ddualBeta   = rnonlinearSpatialMatrix%Dstokes(2,2)
      
      ! Nonlinearity
      if (rnonlinearSpatialMatrix%Dygrad(1,1) .ne. 0.0_DP) then
        roptcoperator%dprimalDelta = dweightConvection*rnonlinearSpatialMatrix%Dygrad(1,1)
        roptcoperator%ddualDelta   = dweightDualConvection*rnonlinearSpatialMatrix%Dygrad(2,2)
        roptcoperator%ddualNewtonTrans = dweightDualNewtonT*rnonlinearSpatialMatrix%DgradyT(2,2)
        
        ! Newton implies additional operators.
        if (rnonlinearSpatialMatrix%Dgrady(1,1) .ne. 0.0_DP) then
          roptcoperator%dprimalNewton    = dweightConvection*1.0_DP
          roptcoperator%ddualRDeltaTrans = dweightConvection*1.0_DP
          roptcoperator%ddualRNewton     = dweightConvection*(-1.0_DP)
        end if
        
      end if
      
      ! Coupling matrices
      !if (rparams%bdualcoupledtoprimal) then
        roptcoperator%ddualRAlpha = rnonlinearSpatialMatrix%Dmass(2,1)
      !end if

      !if (rparams%bcontrolactive) then
        roptcoperator%dcontrolWeight = -rnonlinearSpatialMatrix%Dmass(1,2)*rnonlinearSpatialMatrix%dalphaC
        roptcoperator%dcontrolMultiplier = -1.0_DP/rnonlinearSpatialMatrix%dalphaC
      !end if
      
      if (rnonlinearSpatialMatrix%rdiscrData%rconstraints%ccontrolConstraints .ne. 0) then
        roptcoperator%ccontrolProjection = &
            rnonlinearSpatialMatrix%rdiscrData%rconstraints%ccontrolConstraints

        roptcoperator%cconstraintsType = &
            rnonlinearSpatialMatrix%rdiscrData%rconstraints%cconstraintsType

        roptcoperator%dmin1 = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1
        roptcoperator%dmax1 = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1
        roptcoperator%dmin2 = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2
        roptcoperator%dmax2 = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2

        roptcoperator%p_rumin1 => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rumin1
        roptcoperator%p_rumax1 => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rumax1
        roptcoperator%p_rumin2 => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rumin2
        roptcoperator%p_rumax2 => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rumax2

        roptcoperator%p_rvectorumin => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin
        roptcoperator%p_rvectorumax => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax
      end if
      
      select case (rnonlinearSpatialMatrix%iprimalSol)
      case (1)
        p_rprimalSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1
      case (2)
        p_rprimalSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2
      case (3)
        p_rprimalSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3
      end select
      
      select case (rnonlinearSpatialMatrix%idualSol)
      case (1)
        p_rdualSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1
      case (2)
        p_rdualSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2
      case (3)
        p_rdualSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3
      end select
      
      print *,"This was not updated to alpha <=0"
      call sys_halt()
      
      ! Calculate the velocity-dependent part of the system matrix.
      call conv_strdiffOptC2dgetDefect (&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
          roptcoperator,p_rprimalSol,p_rdualSol,dcx,rx,rd,rcubatureInfo=rcubatureInfo)
    
    end if
    
  contains

    ! -----------------------------------------------------

    subroutine assembleNonlinearityDefect (&
        rnonlinearSpatialMatrix,rmatrix,rvectorPrimalVel,rvectorNonlinearity,&
        rx,rb,dtheta,dygrad,dygradAdj,dygradT,dygradTadj,dgrady,dgradyT,dgradyAdj,&
        rcubatureInfo,rstabilisation,dcx,rmatrixEOJ)
        
    ! Assembles the convection matrix in the block matrix rmatrix at position (1,1):
    !
    ! rmatrix := dygrad*N(rvector) + dygradT*N^t(rvector) +
    !            dgrady*N*(rvector) + dgradyT*N*^t(rvector)
    !
    ! Even if no nonlinearity is present, the routine can be used to
    ! add stabilisation into the matrix.
    
    ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix.
    type(t_nonlinearSpatialMatrix), intent(IN) :: rnonlinearSpatialMatrix
    
    ! 2X2 block matrix that specifies the structure of the velocity FE space.
    type(t_matrixBlock), intent(INOUT) :: rmatrix
    
    ! Primal velocity vector. Used for nonlinear viscosity.
    type(t_vectorBlock), intent(in) :: rvectorPrimalVel
    
    ! Velocity vector for the nonlinearity. Must be specified if
    ! GAMMA <> 0; can be omitted if GAMMA=0.
    type(t_vectorBlock), intent(in) :: rvectorNonlinearity
    
    ! The current solution vector for the velocity (x- and y-velocity)
    type(t_vectorBlock), intent(in) :: rx

    ! The RHS vector; a defect will be created in this vector.
    type(t_vectorBlock), intent(inout) :: rb
    
    ! Weight for the stabilisation and Stokes operator in the Stokes case.
    real(DP), intent(in) :: dtheta

    ! Weight for the nonlinear term u\grad(.)
    real(DP), intent(in) :: dygrad

    ! Weight for the adjoint nonlinear term (.,u\grad(.))
    real(DP), intent(in) :: dygradAdj

    ! Weight for the nonlinear term u(\grad(.))^t
    real(DP), intent(in) :: dygradT

    ! Weight for the adjoint nonlinear term (.,u(\grad(.))^t)
    real(DP), intent(in) :: dygradTAdj

    ! Weight for the nonlinear term (\grad(.))u
    real(DP), intent(in) :: dgrady

    ! Weight for the nonlinear term (\grad(.))^t u
    real(DP), intent(in) :: dgradyT
    
    ! Weight for the adjoint nonlinear term (\grad(.))u
    real(DP), intent(in) :: dgradyAdj

    ! Weight for the operator when multiplying: d = b - dcx * A x. Standard = 1.0_DP
    real(DP), intent(in) :: dcx
    
    ! Cubature information structure that defines how to do the cubature.
    type(t_scalarCubatureInfo), intent(in) :: rcubatureInfo
    
    ! Stabilisation parameters
    type(t_settings_stabil), intent(in) :: rstabilisation

    ! OPTIONAL: EOJ matrix in case the precomputed EOJ stabilisation is active.
    type(t_matrixScalar), intent(in), optional :: rmatrixEOJ
    
    ! local variables
    logical :: bshared
    type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
    type(t_convStreamDiff2) :: rstreamlineDiffusion2
    type(t_convStreamDiff2) :: rstreamlineDiffusion3
    type(t_jumpStabilisation) :: rjumpStabil
    type(t_convupwind) :: rupwindStabil
    type(t_vectorBlock) :: rbtemp
    integer :: h_Idofs
    integer, dimension(:), pointer :: p_Idofs
    type(t_collection), target :: rcollection, ruserCollection
    
    ! DEBUG!!!
    real(dp), dimension(:), pointer :: p_Ddata1,p_Ddata2
    call lsysbl_getbase_double (rx,p_Ddata1)
    call lsysbl_getbase_double (rb,p_Ddata2)
    
      ! Is A11=A22 physically?
      bshared = lsyssc_isMatrixContentShared(&
                    rmatrix%RmatrixBlock(1,1),&
                    rmatrix%RmatrixBlock(2,2))

      if (rstabilisation%cconvectionOnBoundaryDefect .eq. 0) then
        ! Copy the RHS vector so we can restore the values
        call lsysbl_copyVector (rb,rbtemp)
      end if
                    
      if ((dygrad .ne. 0.0_DP) .or. (dygradT .ne. 0.0_DP) .or. &
          (dgrady .ne. 0.0_DP) .or. (dgradyT .ne. 0.0_DP) .or. (dgradyAdj .ne. 0.0_DP) .or. &
          (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%isubequation .ne. 0) .or. &
          (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%cviscoModel .ne. 0)) then
          
        ! Prepare a user defined collection structure for the
        ! callback routines
        call collct_init (ruserCollection)
        call user_initCollectForAssembly (rnonlinearSpatialMatrix%p_rglobalData,&
            rnonlinearSpatialMatrix%p_rnonlinearity%dassociatedTimePrimalVel,ruserCollection)

        ! Type of assembly routine?                      
        select case (rstabilisation%cupwind)
        case (CCMASM_STAB_STREAMLINEDIFF)
        
          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%cviscoModel .ne. 0) then
            ! Not supported by this SD method.
            call output_line (&
                "This assembly method does not support nonconstant viscosity!", &
                OU_CLASS_ERROR,OU_MODE_STD,"cc_assembleMatrix")
            call sys_halt()
          end if

          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst
          
          rstreamlineDiffusion%dtheta = dcx
                    
          ! Set stabilisation parameter
          rstreamlineDiffusion%dupsam = rstabilisation%dupsam
          
          ! Matrix weights
          rstreamlineDiffusion%ddelta            = dygrad
          rstreamlineDiffusion%ddeltaTransposed  = dygradT
          rstreamlineDiffusion%dnewton           = dgrady
          rstreamlineDiffusion%dnewtonTransposed = dgradyT
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvectorNonlinearity!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamlineDiffusionBlk2d (&
                              rvectorNonlinearity, rvectorNonlinearity, &
                              1.0_DP, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODDEFECT, &
                              rmatrix,rx,rb,rcubatureInfo=rcubatureInfo)
                              
          if ((dygradAdj .ne. 0.0_DP) .or. (dygradTAdj .ne. 0.0_DP)) then
            call output_line ("Operator not supported by this SD method!", &
                OU_CLASS_ERROR,OU_MODE_STD,"cc_assembleMatrix")
            call sys_halt()
          end if
                              
        case (CCMASM_STAB_STREAMLINEDIFF2)
        
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst
          
          rstreamlineDiffusion2%dtheta = dcx
                    
          ! Set stabilisation parameter
          rstreamlineDiffusion2%dupsam = rstabilisation%dupsam
          
          ! Matrix weights
          rstreamlineDiffusion2%ddelta  = dygrad
          rstreamlineDiffusion2%ddeltaT = dygradT
          rstreamlineDiffusion2%ddeltaAdj  = dygradAdj
          rstreamlineDiffusion2%ddeltaTAdj = dygradTAdj
          rstreamlineDiffusion2%dnewton = dgrady
          rstreamlineDiffusion2%dnewtonT = dgradyT
          rstreamlineDiffusion2%dnewtonAdj = dgradyAdj
          
          ! Probably, we have nonconstant viscosity.
          ! In that case, init a collection structure for a callback
          ! routine that specifies the viscosity
          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%cviscoModel .ne. 0) then
          
            rstreamlineDiffusion2%bconstnu = .false.
            
            ! Assemble at least the Stokes matrix.
            ! For a more complicated formulation, dbeta/dbetaT
            ! may be changed later.
            rstreamlineDiffusion2%dbeta = dtheta
            
            ! Prepare the collection. The "next" collection points to the user defined
            ! collection.
            rcollection%p_rnextCollection => ruserCollection
            call smva_prepareViscoAssembly (&
                rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal,&
                rcollection,rvectorPrimalVel)
            
          end if
          
          ! Assemble the deformation tensor?
          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%isubequation .ne. 0) then
            rstreamlineDiffusion2%dbeta = 0.5_DP * dtheta
            rstreamlineDiffusion2%dbetaT = 0.5_DP * dtheta
          end if
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvectorNonlinearity!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamDiff2Blk2dDef (rstreamlineDiffusion2,rmatrix,&
              rx,rb,rvectorNonlinearity,ffunctionViscoModel,rcollection,rcubatureInfo)
              
        case (CCMASM_STAB_UPWIND)
        
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst
          
          rstreamlineDiffusion2%dtheta = dcx
                    
          ! Set stabilisation parameter
          rstreamlineDiffusion2%dupsam = 0.0_DP
          
          ! Matrix weights
          rstreamlineDiffusion2%ddelta  = 0.0_DP
          rstreamlineDiffusion2%ddeltaAdj  = dygradAdj
          rstreamlineDiffusion2%ddeltaT = dygradT
          rstreamlineDiffusion2%ddeltaTadj = dygradTadj
          rstreamlineDiffusion2%dnewton = dgrady
          rstreamlineDiffusion2%dnewtonT = dgradyT
          rstreamlineDiffusion2%dnewtonAdj = dgradyAdj
          
          ! Probably, we have nonconstant viscosity.
          ! In that case, init a collection structure for a callback
          ! routine that specifies the viscosity
          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%cviscoModel .ne. 0) then
          
            rstreamlineDiffusion2%bconstnu = .false.
            
            ! Assemble at least the Stokes matrix.
            ! For a more complicated formulation, dbeta/dbetaT
            ! may be changed later.
            rstreamlineDiffusion2%dbeta = dtheta
            
            ! Prepare the collection. The "next" collection points to the user defined
            ! collection.
            rcollection%p_rnextCollection => ruserCollection
            call smva_prepareViscoAssembly (&
                rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal,&
                rcollection,rvectorPrimalVel)
            
          end if

          ! Assemble the deformation tensor?
          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%isubequation .ne. 0) then
            rstreamlineDiffusion2%dbeta = 0.5_DP * dtheta
            rstreamlineDiffusion2%dbetaT = 0.5_DP * dtheta
          end if
          
          ! Call the SD method to calculate the nonlinearity except for the
          ! convection part. As velocity vector, specify rvectorNonlinearity!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamDiff2Blk2dDef (rstreamlineDiffusion2,rmatrix,&
              rx,rb,rvectorNonlinearity,ffunctionViscoModel,rcollection,rcubatureInfo)

          ! Prepare the upwind structure for the assembly of the convection.
          ! Note: Stabilisation weight is actually wrong, but it is not possible
          ! to specify in rupwindStabil%dupsam whether the stabilisation
          ! is added or subtracted!
          rupwindStabil%dupsam = abs(rstabilisation%dupsam)
          rupwindStabil%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst
          rupwindStabil%dtheta = dcx*dygrad
          
          ! Set stabilisation parameter
          rjumpStabil%dgamma = abs(rstabilisation%dupsam)
          
          ! Apply the upwind operator.
          call conv_upwind2d (rvectorNonlinearity, rvectorNonlinearity, 1.0_DP, 0.0_DP,&
              rupwindStabil, CONV_MODDEFECT, rmatrix%RmatrixBlock(1,1), rx, rb)

          ! Prepare the upwind structure for the assembly of the convection.
          !rstreamlineDiffusion3%dupsam = rstabilisation%dupsam
          !rstreamlineDiffusion3%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
          !rstreamlineDiffusion3%dtheta = dygrad*dcx
          !rstreamlineDiffusion3%ddelta = 1.0_DP
          
          ! Apply the upwind operator.
          !call conv_streamDiff2Blk2dDef (rstreamlineDiffusion3,rmatrix,&
          !    rx,rb,rvectorNonlinearity)

        case (CCMASM_STAB_EDGEORIENTED)
          ! Jump stabilisation.

          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%cviscoModel .ne. 0) then
            ! Not supported by this SD method.
            call output_line (&
                "This assembly method does not support nonconstant viscosity!", &
                OU_CLASS_ERROR,OU_MODE_STD,"cc_assembleMatrix")
            call sys_halt()
          end if

          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst
          
          rstreamlineDiffusion%dtheta = dcx
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion%ddelta            = dygrad
          rstreamlineDiffusion%ddeltaTransposed  = dygradT
          rstreamlineDiffusion%dnewton           = dgrady
          rstreamlineDiffusion%dnewtonTransposed = dgradyT
          
          if ((dygradAdj .ne. 0.0_DP) .or. (dygradTAdj .ne. 0.0_DP)) then
            call output_line ("Operator not supported by this SD method!", &
                OU_CLASS_ERROR,OU_MODE_STD,"cc_assembleMatrix")
            call sys_halt()
          end if
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvectorNonlinearity!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamlineDiffusionBlk2d (&
                              rvectorNonlinearity, rvectorNonlinearity, &
                              1.0_DP, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODDEFECT, &
                              rmatrix,rx,rb,rcubatureInfo=rcubatureInfo)
                              
          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rstreamlineDiffusion%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgamma = abs(rstabilisation%dupsam)
          
          ! Matrix weight. Compensate for any "-" sign in dygrad!
          rjumpStabil%dtheta = dcx*dygrad*mprim_signum(rstabilisation%dupsam)

          ! Call the jump stabilisation technique to stabilise that stuff.
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
                              rjumpStabil, CONV_MODDEFECT, &
                              rmatrix%RmatrixBlock(1,1),rx,rb)
          if (rstabilisation%ceojStabilOnBoundary .eq. 0) then
            call smva_addBdEOJvector (rjumpStabil,-1.0_DP,rmatrix%RmatrixBlock(1,1),rx,rb)
          end if

!          if (.not. bshared) then
!            call conv_jumpStabilisation2d (&
!                                rjumpStabil, CONV_MODDEFECT, &
!                                rmatrix%RmatrixBlock(2,2),rx,rb)
!          end if

        case (CCMASM_STAB_EDGEORIENTED2)
          ! Jump stabilisation.
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst
          
          rstreamlineDiffusion2%dtheta = dcx
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion2%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion2%ddelta   = dygrad
          rstreamlineDiffusion2%ddeltaT  = dygradT
          rstreamlineDiffusion2%ddeltaAdj   = dygradAdj
          rstreamlineDiffusion2%ddeltaTAdj  = dygradTAdj
          rstreamlineDiffusion2%dnewton  = dgrady
          rstreamlineDiffusion2%dnewtonT = dgradyT
          rstreamlineDiffusion2%dnewtonAdj = dgradyAdj
          
          ! Probably, we have nonconstant viscosity.
          ! In that case, init a collection structure for a callback
          ! routine that specifies the viscosity
          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%cviscoModel .ne. 0) then
          
            rstreamlineDiffusion2%bconstnu = .false.
            
            ! Assemble at least the Stokes matrix.
            ! For a more complicated formulation, dbeta/dbetaT
            ! may be changed later.
            rstreamlineDiffusion2%dbeta = dtheta
            
            ! Prepare the collection. The "next" collection points to the user defined
            ! collection.
            rcollection%p_rnextCollection => ruserCollection
            call smva_prepareViscoAssembly (&
                rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal,&
                rcollection,rvectorPrimalVel)
            
          end if

          ! Assemble the deformation tensor?
          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%isubequation .ne. 0) then
            rstreamlineDiffusion2%dbeta = 0.5_DP * dtheta
            rstreamlineDiffusion2%dbetaT = 0.5_DP * dtheta
          end if
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvectorNonlinearity!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamDiff2Blk2dDef (rstreamlineDiffusion2,rmatrix,&
              rx,rb,rvectorNonlinearity,ffunctionViscoModel,rcollection,rcubatureInfo)
                              
          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rstreamlineDiffusion%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgamma = abs(rstabilisation%dupsam)
          
          ! Matrix weight
          rjumpStabil%dtheta = dcx*dygrad*mprim_signum(rstabilisation%dupsam)

          ! Call the jump stabilisation technique to stabilise that stuff.
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          !
          ! The routine processes the X as well as the Y velocity!
          call conv_jumpStabilisation2d (&
                              rjumpStabil, CONV_MODDEFECT, &
                              rmatrix%RmatrixBlock(1,1),rx,rb)
          if (rstabilisation%ceojStabilOnBoundary .eq. 0) then
            call smva_addBdEOJvector (rjumpStabil,-1.0_DP,rmatrix%RmatrixBlock(1,1),rx,rb)
          end if

        case (CCMASM_STAB_EDGEORIENTED3)
        
          ! Jump stabilisation, precomputed matrix.
          !
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst
          
          rstreamlineDiffusion2%dtheta = dcx
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion2%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion2%ddelta   = dygrad
          rstreamlineDiffusion2%ddeltaT  = dygradT
          rstreamlineDiffusion2%ddeltaAdj   = dygradAdj
          rstreamlineDiffusion2%ddeltaTAdj  = dygradTAdj
          rstreamlineDiffusion2%dnewton  = dgrady
          rstreamlineDiffusion2%dnewtonT = dgradyT
          rstreamlineDiffusion2%dnewtonAdj = dgradyAdj
          
          ! Probably, we have nonconstant viscosity.
          ! In that case, init a collection structure for a callback
          ! routine that specifies the viscosity
          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%cviscoModel .ne. 0) then
          
            rstreamlineDiffusion2%bconstnu = .false.
            
            ! Assemble at least the Stokes matrix.
            ! For a more complicated formulation, dbeta/dbetaT
            ! may be changed later.
            rstreamlineDiffusion2%dbeta = dtheta
            
            ! Prepare the collection. The "next" collection points to the user defined
            ! collection.
            rcollection%p_rnextCollection => ruserCollection
            call smva_prepareViscoAssembly (&
                rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal,&
                rcollection,rvectorPrimalVel)
            
          end if

          ! Assemble the deformation tensor?
          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%isubequation .ne. 0) then
            rstreamlineDiffusion2%dbeta = 0.5_DP * dtheta
            rstreamlineDiffusion2%dbetaT = 0.5_DP * dtheta
          end if
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvectorNonlinearity!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamDiff2Blk2dDef (rstreamlineDiffusion2,rmatrix,&
              rx,rb,rvectorNonlinearity,ffunctionViscoModel,rcollection,rcubatureInfo)
                              
          ! Just do some mat-vec's to compute the defect.
          !
          ! The stabilisation parameter is already incorporated into
          ! the matrix -- only its sign is not incorporated, so we do that here.
          call lsyssc_scalarMatVec(rmatrixEOJ,&
              rx%RvectorBlock(1),rb%RvectorBlock(1),&
              -dcx*dygrad*mprim_signum(rstabilisation%dupsam),1.0_DP)

          call lsyssc_scalarMatVec(rmatrixEOJ,&
              rx%RvectorBlock(2),rb%RvectorBlock(2),&
              -dcx*dygrad*mprim_signum(rstabilisation%dupsam),1.0_DP)

        case default
          print *,'Don''t know how to set up nonlinearity!?!'
          stop
        
        end select

        ! Release the collection stuff
        call user_doneCollectForAssembly (rnonlinearSpatialMatrix%p_rglobalData,ruserCollection)
        call collct_done (ruserCollection)

      else
      
        ! That's the Stokes-case. Jump stabilisation is possible...
        
        if ((rstabilisation%dupsam .ne. 0.0_DP) .and. (dtheta .ne. 0.0_DP)) then
          
          ! Prepare a user defined collection structure for the
          ! callback routines
          call collct_init (ruserCollection)
          call user_initCollectForAssembly (rnonlinearSpatialMatrix%p_rglobalData,&
              rnonlinearSpatialMatrix%p_rnonlinearity%dassociatedTimePrimalVel,ruserCollection)

          ! Probably, we have nonconstant viscosity.
          ! In that case, init a collection structure for a callback
          ! routine that specifies the viscosity
          if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%cviscoModel .ne. 0) then
          
            ! Set up the SD structure for the creation of the defect.
            rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst
            
            ! Set stabilisation parameter
            rstreamlineDiffusion2%dupsam = 0.0_DP
            
            ! Matrix weights
            rstreamlineDiffusion2%ddelta  = 0.0_DP
            rstreamlineDiffusion2%ddeltaT = dygradT
            rstreamlineDiffusion2%ddeltaAdj   = dygradAdj
            rstreamlineDiffusion2%ddeltaTAdj  = dygradTAdj
            rstreamlineDiffusion2%dnewton = dgrady
            rstreamlineDiffusion2%dnewtonT = dgradyT
            rstreamlineDiffusion2%dnewtonAdj = dgradyAdj
  
            rstreamlineDiffusion2%bconstnu = .false.
            
            ! Assemble at least the Stokes matrix.
            ! For a more complicated formulation, dbeta/dbetaT
            ! may be changed later.
            rstreamlineDiffusion2%dbeta = dtheta
            
            ! Prepare the collection. The "next" collection points to the user defined
            ! collection.
            rcollection%p_rnextCollection => ruserCollection
            call smva_prepareViscoAssembly (&
                rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal,&
                rcollection,rvectorPrimalVel)
            
            ! Assemble the deformation tensor?
            if (rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%isubequation .ne. 0) then
              rstreamlineDiffusion2%dbeta = 0.5_DP * dtheta
              rstreamlineDiffusion2%dbetaT = 0.5_DP * dtheta
            end if
          
            ! Call the SD method to calculate the nonlinearity for everything except
            ! for the convection. As velocity vector, specify rvectorNonlinearity!
            ! Therefore, the primal velcity is always used for assembling
            ! that thing!
            call conv_streamDiff2Blk2dMat (&
                rstreamlineDiffusion2,rmatrix,rvectorNonlinearity,&
                ffunctionViscoModel,rcollection,rcubatureInfo)

          end if
                    
          select case (rstabilisation%cupwind)
          case (CCMASM_STAB_EDGEORIENTED,CCMASM_STAB_EDGEORIENTED2)
            
            ! Set up the jump stabilisation structure.
            ! There's not much to do, only initialise the viscosity...
            rjumpStabil%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnuConst
            
            ! Set stabilisation parameter
            rjumpStabil%dgamma = abs(rstabilisation%dupsam)
            
            ! Matrix weight. Take the weight of the Stokes operator.
            rjumpStabil%dtheta = dcx*dtheta

            ! Call the jump stabilisation technique to stabilise that stuff.
            ! We can assemble the jump part any time as it's independent of any
            ! convective parts...
            call conv_jumpStabilisation2d (&
                                rjumpStabil, CONV_MODDEFECT, &
                                rmatrix%RmatrixBlock(1,1),rx,rb)
            if (rstabilisation%ceojStabilOnBoundary .eq. 0) then
              call smva_addBdEOJvector (rjumpStabil,-1.0_DP,rmatrix%RmatrixBlock(1,1),rx,rb)
            end if

!            if (.not. bshared) then
!              call conv_jumpStabilisation2d (&
!                                  rjumpStabil, CONV_MODDEFECT, &
!                                  rmatrix%RmatrixBlock(2,2),rx,rb)
!            end if
            
          case (CCMASM_STAB_EDGEORIENTED3)
          
            ! Jump stabilisation, precomputed matrix.
            !
            ! Just do some mat-vec's to compute the defect.
            ! The stabilisation parameter is already incorporated into
            ! the matrix -- only its sign is not incorporated, so we do that here.
            call lsyssc_scalarMatVec(rmatrixEOJ,&
                rx%RvectorBlock(1),rb%RvectorBlock(1),&
                -dcx*dtheta,1.0_DP)

            call lsyssc_scalarMatVec(rmatrixEOJ,&
                rx%RvectorBlock(2),rb%RvectorBlock(2),&
                -dcx*dtheta,1.0_DP)

          case default
            ! No stabilisation
          
          end select

          ! Release the collection stuff
          call user_doneCollectForAssembly (rnonlinearSpatialMatrix%p_rglobalData,ruserCollection)
          call collct_done (ruserCollection)

        end if
        
      end if
      
      if (rstabilisation%cconvectionOnBoundaryDefect .eq. 0) then
        ! Restore the RHS values on the boundary.
        h_Idofs = ST_NOHANDLE
        call bcasm_getDOFsOnBoundary (&
            rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimal%RspatialDiscr(1), h_Idofs)
        if (h_Idofs .ne. ST_NOHANDLE) then
          call storage_getbase_int (h_Idofs,p_Idofs)
          call lsyssc_replaceVectorEntries (rbtemp%RvectorBlock(1),rb%RvectorBlock(1),p_Idofs)
          call lsyssc_replaceVectorEntries (rbtemp%RvectorBlock(2),rb%RvectorBlock(2),p_Idofs)
          call storage_free (h_Idofs)
        end if
        call lsysbl_releaseVector (rbtemp)
      end if
      
    end subroutine

    subroutine assemblePrimalUConstrMassDefect (&
        rnonlinearSpatialMatrix,rvector,rdefect,dcx,rvelocityVector) !,cprojectionType)
        
    ! Assembles the defect arising from the projective coupling mass
    ! matrices in the primal equation which comes from constraints on u.
    ! rdefect must have been initialised with the right hand side vector.
    !
    ! Let the mass matrix M~ be the usual mass matrix where u is ok
    ! and the 0-operator where u is out of bounds.
    ! Then, we assemble
    !
    !       rdefect = r(primal)defect - dcx (dmu1 M~ r(dual)vector)

    ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
    ! about how the mass matrix part is weighted (dr11, dr12).
    type(t_nonlinearSpatialMatrix), intent(IN) :: rnonlinearSpatialMatrix

    ! Solution vector.
    type(t_vectorBlock), intent(IN) :: rvector
    
    ! On entry: RHS vector.
    ! Is overwritten by the defect vector in the velocity subsystem.
    type(t_vectorBlock), intent(INOUT) :: rdefect
    
    ! Multiplication factor for the whole operator A*rvector
    real(DP), intent(IN) :: dcx
    
    ! Velocity vector field that should be used for the assembly of the
    ! nonlinearity. Block 1 and 2 in that block vector are used as velocity
    ! field and shall contain the dual velocity field.
    type(t_vectorBlock), intent(IN), target :: rvelocityVector

    ! Type of projection.
    ! =0: Simple DOF-based projection, rows in the mass matrices are
    !     replaced by zero (faster)
    ! =1: Reassembly of the mass matrices, projection based on cubature points.
    !integer, intent(in) :: cprojectionType

      ! local variables
      type(t_matrixBlock) :: rtempmatrix
      type(t_vectorBlock) :: rtempvectorEval,rtempVectorDef
      integer, dimension(:), pointer :: p_IelementList
      integer :: ccontrolConstraints
      type(t_scalarCubatureInfo) :: rcubatureInfo, rcubatureInfoAdapt
      
      ! If we have a reactive coupling mass matrix, it gets interesting...
      if (dcx .ne. 0.0_DP) then

        call lsysbl_createEmptyMatrix (rtempMatrix,2)
        
        ccontrolConstraints = rnonlinearSpatialMatrix%rdiscrData%rconstraints%ccontrolConstraints

        ! The ccontrolConstraints decides on whether we use the 'quick-setup'
        ! method for the mass matrices or rebuild them again.
        ! ccontrolConstraints=0 uses standard mass matrices.
        select case (ccontrolConstraints)
        case (0)
        
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
           
          call lsysbl_updateMatStrucInfo (rtempMatrix)
             
        case (1)
        
          select case (rnonlinearSpatialMatrix%rdiscrData%rconstraints%cconstraintsType)
          case (0)
            ! Constant bounds
          
            call nwder_minMaxProjByMass (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
                -rnonlinearSpatialMatrix%dalphaC,rtempMatrix%RmatrixBlock(1,1),&
                .true.,&
                -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvelocityVector%RvectorBlock(4),&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1,&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1)

            call nwder_minMaxProjByMass (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
                -rnonlinearSpatialMatrix%dalphaC,rtempMatrix%RmatrixBlock(2,2),&
                .true.,&
                -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvelocityVector%RvectorBlock(5),&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2,&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2)

          case (1)
            ! Variable bounds

            call nwder_minMaxProjByMass (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
                -rnonlinearSpatialMatrix%dalphaC,rtempMatrix%RmatrixBlock(1,1),&
                .true.,&
                -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvelocityVector%RvectorBlock(4),&
                1.0_DP,1.0_DP,&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(1),&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(1))

            call nwder_minMaxProjByMass (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
                -rnonlinearSpatialMatrix%dalphaC,rtempMatrix%RmatrixBlock(2,2),&
                .true.,&
                -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvelocityVector%RvectorBlock(5),&
                1.0_DP,1.0_DP,&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(2),&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(2))

          case default
            ! Not implemented.
            call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
            call sys_halt()
          end select

          call lsysbl_updateMatStrucInfo (rtempMatrix)
             
        case (2)
        
          ! Create an empty matrix with the structure we need. Share the structure
          ! of the mass matrix. Entries are not necessary for the assembly
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
              
          call lsysbl_updateMatStrucInfo (rtempMatrix)

          ! Clear the matrices in advance.
          call lsyssc_clearMatrix (rtempMatrix%RmatrixBlock(1,1))
          call lsyssc_clearMatrix (rtempMatrix%RmatrixBlock(2,2))

          ! For alpha=0, we have bang-bang-control. In this case,
          ! the Newton derivative is the zero operator, so nothing
          ! has to be assembled here.
          !
          ! Alpha < 0 switches the distributed control off, so also in this
          ! case, nothing has to be assembled.
          if (rnonlinearSpatialMatrix%dalphaC .gt. 0.0_DP) then

            ! Create a cubature info structure that defines the cubature
            ! rule to use.            
            call spdiscr_createDefCubStructure (&
                rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimalDual%RspatialDiscr(1),&
                rcubatureInfo,rnonlinearSpatialMatrix%rdiscrData%rsettingsSpaceDiscr%icubMass)

            ! Calculate the Newton derivative.
            !
            ! The operator weight is just 1.0 here; the matrix is scaled later.
            ! However, cancel out the -1/alpha from the function weight, which
            ! is included into the operator as part of the Newton derivative.

            select case (rnonlinearSpatialMatrix%rdiscrData%rconstraints%cconstraintsType)
            case (0)
              ! Constant bounds

              call nwder_minMaxProjByCubature (-rnonlinearSpatialMatrix%dalphaC,&
                  rtempmatrix%RmatrixBlock(1,1),rcubatureInfo,&
                  -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvelocityVector%RvectorBlock(4),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1)

              call nwder_minMaxProjByCubature (-rnonlinearSpatialMatrix%dalphaC,&
                  rtempmatrix%RmatrixBlock(2,2),rcubatureInfo,&
                  -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvelocityVector%RvectorBlock(5),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2)

            case (1)
              ! Variable bounds

              call nwder_minMaxProjByCubature (-rnonlinearSpatialMatrix%dalphaC,&
                  rtempmatrix%RmatrixBlock(1,1),rcubatureInfo,&
                  -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvelocityVector%RvectorBlock(4),&
                  1.0_DP,1.0_DP,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(1),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(1))

              call nwder_minMaxProjByCubature (-rnonlinearSpatialMatrix%dalphaC,&
                  rtempmatrix%RmatrixBlock(2,2),rcubatureInfo,&
                  -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvelocityVector%RvectorBlock(5),&
                  1.0_DP,1.0_DP,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(2),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(2))
            
            end select
              
            ! Release the cubature structure
            call spdiscr_releaseCubStructure(rcubatureInfo)
          
          end if
          
        case (3)
        
          ! Create a matrix with the structure we need. Share the structure
          ! of the mass matrix. Entries are not necessary for the assembly
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

          call lsyssc_clearMatrix (rtempMatrix%RmatrixBlock(1,1))
          call lsyssc_clearMatrix (rtempMatrix%RmatrixBlock(2,2))

          call lsysbl_updateMatStrucInfo (rtempMatrix)

          select case (rnonlinearSpatialMatrix%rdiscrData%rconstraints%cconstraintsType)
          case (0)
            ! Constant bounds

            call nwder_minMaxProjByApproxDer (-rnonlinearSpatialMatrix%dalphaC,&
                rtempmatrix%RmatrixBlock(1,1),0.001_DP,&
                -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvelocityVector%RvectorBlock(4),&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1,&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1)

            call nwder_minMaxProjByApproxDer (-rnonlinearSpatialMatrix%dalphaC,&
                rtempmatrix%RmatrixBlock(2,2),0.001_DP,&
                -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvelocityVector%RvectorBlock(5),&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2,&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2)

          case (1)
            ! Variable bounds

            call nwder_minMaxProjByApproxDer (-rnonlinearSpatialMatrix%dalphaC,&
                rtempmatrix%RmatrixBlock(1,1),0.001_DP,&
                -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvelocityVector%RvectorBlock(4),&
                1.0_DP,1.0_DP,&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(1),&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(1))

            call nwder_minMaxProjByApproxDer (-rnonlinearSpatialMatrix%dalphaC,&
                rtempmatrix%RmatrixBlock(2,2),0.001_DP,&
                -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvelocityVector%RvectorBlock(5),&
                1.0_DP,1.0_DP,&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(2),&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(2))
          
          end select

          call lsysbl_updateMatStrucInfo (rtempMatrix)
            
        case (4)
        
          ! Create an empty matrix with the structure we need. Share the structure
          ! of the mass matrix. Entries are not necessary for the assembly
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
              
          call lsysbl_updateMatStrucInfo (rtempMatrix)

          ! Clear the matrices in advance.
          call lsyssc_clearMatrix (rtempMatrix%RmatrixBlock(1,1))
          call lsyssc_clearMatrix (rtempMatrix%RmatrixBlock(2,2))

          ! For alpha=0, we have bang-bang-control. In this case,
          ! the Newton derivative is the zero operator, so nothing
          ! has to be assembled here.
          !
          ! Alpha < 0 switches the distributed control off, so also in this
          ! case, nothing has to be assembled.
          if (rnonlinearSpatialMatrix%dalphaC .gt. 0.0_DP) then

            ! Create a cubature info structure that defines the cubature
            ! rule to use.            
            call spdiscr_createDefCubStructure (&
                rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimalDual%RspatialDiscr(1),&
                rcubatureInfo,rnonlinearSpatialMatrix%rdiscrData%rsettingsSpaceDiscr%icubMass)

            ! And another one for the adaptive cubature rule on the border
            ! of the active set.
            call spdiscr_createDefCubStructure (&
                rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimalDual%RspatialDiscr(1),&
                rcubatureInfoAdapt,cub_getSummedCubType(&
                    rnonlinearSpatialMatrix%rdiscrData%rsettingsSpaceDiscr%icubMass,4))

            ! Calculate the Newton derivative.
            !
            ! The operator weight is just 1.0 here; the matrix is scaled later.
            ! However, cancel out the -1/alpha from the function weight, which
            ! is included into the operator as part of the Newton derivative.

            select case (rnonlinearSpatialMatrix%rdiscrData%rconstraints%cconstraintsType)
            case (0)
              ! Constant bounds

              call nwder_minMaxProjByAdaptCub (-rnonlinearSpatialMatrix%dalphaC,&
                  rtempmatrix%RmatrixBlock(1,1),rcubatureInfo,rcubatureInfoAdapt,&
                  -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvelocityVector%RvectorBlock(4),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1)

              call nwder_minMaxProjByAdaptCub (-rnonlinearSpatialMatrix%dalphaC,&
                  rtempmatrix%RmatrixBlock(2,2),rcubatureInfo,rcubatureInfoAdapt,&
                  -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvelocityVector%RvectorBlock(5),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2)

            case (1)
              ! Variable bounds

              call nwder_minMaxProjByAdaptCub (-rnonlinearSpatialMatrix%dalphaC,&
                  rtempmatrix%RmatrixBlock(1,1),rcubatureInfo,rcubatureInfoAdapt,&
                  -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvelocityVector%RvectorBlock(4),&
                  1.0_DP,1.0_DP,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(1),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(1))

              call nwder_minMaxProjByAdaptCub (-rnonlinearSpatialMatrix%dalphaC,&
                  rtempmatrix%RmatrixBlock(2,2),rcubatureInfo,rcubatureInfoAdapt,&
                  -1.0_DP/rnonlinearSpatialMatrix%dalphaC,rvelocityVector%RvectorBlock(5),&
                  1.0_DP,1.0_DP,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(2),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(2))
            
            end select
              
            ! Release the cubature structure
            call spdiscr_releaseCubStructure(rcubatureInfoAdapt)
            call spdiscr_releaseCubStructure(rcubatureInfo)
          
          end if

        case default
        
          ! Cancel
          call lsysbl_releaseMatrix (rtempMatrix)
          return
          
        end select

        ! Create a temporary block vector that points to the dual velocity.
        ! This has to be evaluated during the assembly.
        call lsysbl_deriveSubvector (rvector,rtempvectorEval,4,5,.true.)
        
        ! Create a temporary block vector for the dual defect.
        ! Matrix*primal velocity is subtracted from this.
        call lsysbl_deriveSubvector (rdefect,rtempvectorDef,1,2,.true.)

        ! Create the defect
        call lsysbl_blockMatVec (rtempmatrix, rtempvectorEval, rtempvectorDef, -dcx, 1.0_DP)
        
        ! Release memory
        call lsysbl_releaseVector (rtempvectorDef)
        call lsysbl_releaseVector (rtempvectorEval)
        call lsysbl_releaseMatrix (rtempMatrix)
      
      end if
      
    end subroutine

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_prjControl (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! For distributed control.
    ! Coefficients in the bilinear form of the projection operator.
    ! Constant constraints.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT) :: Dcoefficients
  !</output>
    
  !</subroutine>
  
      ! local variables
    type(t_vectorBlock), pointer :: p_rvector
    real(dp), dimension(:,:), allocatable :: Dfunc
    integer(I32) :: celement
    real(DP) :: da, db
    integer :: ipt, iel
    
    ! Get the bounds
    da = rcollection%DquickAccess(1)
    db = rcollection%DquickAccess(2)
    
    ! Get a pointer to the FE solution from the collection.
    ! The routine below wrote a pointer to the vector to the
    ! first quick-access vector pointer in the collection.
    p_rvector => rcollection%p_rvectorQuickAccess1

    ! Allocate memory for the function values in the cubature points:
    allocate(Dfunc(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))
    
    ! Calculate the function value of the solution vector in all
    ! our cubature points:
    !
    ! Figure out the element type, then call the
    ! evaluation routine for a prepared element set.
    ! This works only if the trial space of the matrix coincides
    ! with the FE space of the vector T we evaluate!
    
    celement = rdomainIntSubset%celement
    
    call fevl_evaluate_sim (p_rvector%RvectorBlock(1), &
        rdomainIntSubset%p_revalElementSet, &
        celement, rdomainIntSubset%p_IdofsTrial, DER_FUNC, Dfunc)
    
    ! Now check the function values lambda.
    ! Return the projected control in every cubature point.
    do iel = 1,ubound(Dcoefficients,3)
      do ipt = 1,ubound(Dcoefficients,2)
        Dcoefficients(1,ipt,iel) = min(max(Dfunc(ipt,iel),da),db)
      end do
    end do
    
    ! Release memory
    deallocate(Dfunc)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine coeff_prjControlVec (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! For distributed control.
    ! Coefficients in the bilinear form of the projection operator.
    ! Variable constraints, given by a FEM function.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT) :: Dcoefficients
  !</output>
    
  !</subroutine>
  
      ! local variables
    type(t_vectorBlock), pointer :: p_rvector,p_rvectorUmin,p_rvectorUmax
    real(dp), dimension(:,:), allocatable :: Dfunc
    real(dp), dimension(:), allocatable :: Dumin,Dumax
    integer, dimension(:), allocatable :: Ielements
    integer(I32) :: celement
    integer :: ipt, iel
    
    ! Get a pointer to the FE solution from the collection.
    ! The routine below wrote a pointer to the vector to the
    ! first quick-access vector pointer in the collection.
    p_rvector => rcollection%p_rvectorQuickAccess1
    p_rvectorUmin => rcollection%p_rvectorQuickAccess2
    p_rvectorUmax => rcollection%p_rvectorQuickAccess3

    ! Allocate memory for the function values in the cubature points.
    ! Function value, minimum and maximum.
    allocate(Dfunc(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))
    
    ! Allocate temp memory for element hints and min/max values for u.
    allocate(Dumin(ubound(Dcoefficients,3)))
    allocate(Dumax(ubound(Dcoefficients,3)))
    allocate(Ielements(ubound(Dcoefficients,2)))
    
    ! Calculate the function value of the solution vector in all
    ! our cubature points:
    !
    ! Figure out the element type, then call the
    ! evaluation routine for a prepared element set.
    ! This works only if the trial space of the matrix coincides
    ! with the FE space of the vector we evaluate!
    
    celement = rdomainIntSubset%celement
    
    call fevl_evaluate_sim (p_rvector%RvectorBlock(1), &
        rdomainIntSubset%p_revalElementSet, &
        celement, rdomainIntSubset%p_IdofsTrial, DER_FUNC, Dfunc(:,:))
    
    ! Now check the function values lambda.
    ! Return the projected control in every cubature point.
    do iel = 1,ubound(Dcoefficients,3)
    
      ! Evaluate min and max value of the control in the cubature points
      ! on the current element.
      Ielements(:) = rdomainIntSubset%p_Ielements(iel)

      call fevl_evaluate (DER_FUNC, Dumin, p_rvectorUmin%RvectorBlock(1), &
          Dpoints(:,:,iel), IelementsHint=Ielements)

      call fevl_evaluate (DER_FUNC, Dumax, p_rvectorUmin%RvectorBlock(2), &
          Dpoints(:,:,iel), IelementsHint=Ielements)

      ! Calculate the projection in the cubature points
      do ipt = 1,ubound(Dcoefficients,2)
        Dcoefficients(1,ipt,iel) = min(max(Dfunc(ipt,iel),Dumin(ipt)),Dumax(ipt))
      end do
    end do
    
    ! Release memory
    deallocate(Dumin)
    deallocate(Dumax)
    deallocate(Ielements)
    
    deallocate(Dfunc)

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_initNonlinearData (rnonlinearData,rvector1,rvector2,rvector3,&
      dassociatedTimePrimalVel,dassociatedTimeDualVel,&
      rneumannBoundary,rdirichletBCCBoundary)

!<description>
  ! Initialises a nonlinear-data structure that defines the nonlinearity
  ! in a nonlinear matrix.
!</description>

!<input>
    ! Specifies where to evaluate the nonlinearity and must contain the data
    ! for the 'previous' timestep. If there is no previous timestep (e.g.
    ! like in the 0th timestep), the vector can be undefined.
    type(t_vectorBlock), intent(in), target :: rvector1

    ! Specifies where to evaluate the nonlinearity and must contain the data
    ! for the 'current' timestep.
    type(t_vectorBlock), intent(in), target :: rvector2

    ! Specifies where to evaluate the nonlinearity and must contain the data
    ! for the 'next' timestep. If there is no next timestep (e.g.
    ! like in the last timestep), the vector can be undefined.
    type(t_vectorBlock), intent(in), target :: rvector3

    ! Point in time associated to this matrix.
    real(DP), intent(in) :: dassociatedTimePrimalVel
    real(DP), intent(in) :: dassociatedTimeDualVel
    
    ! Specifies the Neumann boundary.
    type(t_boundaryRegionList), intent(in), target :: rneumannBoundary
    
    ! Specifies the Dirichlet control boundary.
    type(t_boundaryRegionList), intent(in), target :: rdirichletBCCBoundary
    
!</input>

!<inputoutput>
  ! Nonlinear-data structure to be created.
  type(t_spatialMatrixNonlinearData), intent(out) :: rnonlinearData
!</inputoutput>

!</subroutine>
 
    ! Set the pointers.
    rnonlinearData%p_rvector1 => rvector1
    rnonlinearData%p_rvector2 => rvector2
    rnonlinearData%p_rvector3 => rvector3
    rnonlinearData%p_rneumannBoundary => rneumannBoundary
    rnonlinearData%p_rdirichletBCCBoundary => rdirichletBCCBoundary

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_addBdEOJOperator (rjumpStabil,dweight,rmatrix)

!<description>
  ! Adds the jump stabilisation on all boundary edges to a matrix.
!</description>

!<input>
    ! Jump stabilisation structure that defines the jump stabilisation operator.
    type(t_jumpStabilisation), intent(in) :: rjumpStabil
    
    ! Weight of the operator. Default = 1.0
    real(DP), intent(in) :: dweight
!</input>

!<inputoutput>
    ! Matrix to be modified.
    type(t_matrixScalar), intent(inout), target :: rmatrix
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_jumpStabilisation) :: rjumpStabilLocal
    type(t_boundaryRegion) :: rboundaryRegion
    integer, dimension(:), allocatable :: p_Iedges, p_Itemp
    integer :: iedge,ibc,nedgecount
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_boundary), pointer :: p_rboundary
    
    p_rtriangulation => rmatrix%p_rspatialDiscrTrial%p_rtriangulation
    p_rboundary => rmatrix%p_rspatialDiscrTrial%p_rboundary
    
    ! prepare a negative operator.
    rjumpStabilLocal = rjumpStabil
    
    rjumpStabilLocal%dtheta = rjumpStabilLocal%dtheta*dweight

    allocate (p_Iedges(p_rtriangulation%NMT))
    allocate (p_Itemp(p_rtriangulation%NMT))

    ! Loop over the edges and boundary components
    do ibc = 1,boundary_igetNBoundComp(p_rboundary)
      do iedge = 1,boundary_igetNsegments(p_rboundary, ibc)

        ! Get a region identifying that boundary edge
        call boundary_createRegion(p_rboundary,ibc,iedge,rboundaryRegion)

        ! Get triangulation edges on that boundary edge
        call bcasm_getEdgesInBdRegion (p_rtriangulation,rboundaryRegion, &
            nedgecount, p_Iedges, p_Itemp)

        ! Subtract the jump
        call conv_JumpStabilisation2d(rjumpStabilLocal,CONV_MODMATRIX,rmatrix, &
            InodeList=p_Iedges(1:nedgecount))

      end do
    end do

    deallocate (p_Iedges,p_Itemp)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_addBdEOJvector (rjumpStabil,dweight,rmatrix,rx,rb)

!<description>
  ! Adds the jump stabilisation on all boundary edges.
!</description>

!<input>
    ! Jump stabilisation structure that defines the jump stabilisation operator.
    type(t_jumpStabilisation), intent(in) :: rjumpStabil
    
    ! Weight of the operator. Default = 1.0.
    real(DP), intent(in) :: dweight

    ! Solution vector.
    type(t_vectorBlock), intent(in) :: rx
    
    ! Matrix that specifies the sparsity pattern of the operator.
    type(t_matrixScalar), intent(inout) :: rmatrix
!</input>

!<inputoutput>
    ! Defect vector to be modified.
    type(t_vectorBlock), intent(inout), target :: rb
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_jumpStabilisation) :: rjumpStabilLocal
    type(t_boundaryRegion) :: rboundaryRegion
    integer, dimension(:), allocatable :: p_Iedges, p_Itemp
    integer :: iedge,ibc,nedgecount
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_boundary), pointer :: p_rboundary
    
    p_rtriangulation => rmatrix%p_rspatialDiscrTrial%p_rtriangulation
    p_rboundary => rmatrix%p_rspatialDiscrTrial%p_rboundary
    
    ! prepare a negative operator.
    rjumpStabilLocal = rjumpStabil
    
    rjumpStabilLocal%dtheta = rjumpStabilLocal%dtheta*dweight

    allocate (p_Iedges(p_rtriangulation%NMT))
    allocate (p_Itemp(p_rtriangulation%NMT))

    ! Loop over the edges and boundary components
    do ibc = 1,boundary_igetNBoundComp(p_rboundary)
      do iedge = 1,boundary_igetNsegments(p_rboundary, ibc)

        ! Get a region identifying that boundary edge
        call boundary_createRegion(p_rboundary,ibc,iedge,rboundaryRegion)

        ! Get triangulation edges on that boundary edge
        call bcasm_getEdgesInBdRegion (p_rtriangulation,rboundaryRegion, &
            nedgecount, p_Iedges, p_Itemp)

        ! Add the jump
        call conv_jumpStabilisation2d (&
                              rjumpStabilLocal, CONV_MODDEFECT, &
                              rmatrix,rx,rb,&
                              InodeList=p_Iedges(1:nedgecount))

      end do
    end do

    deallocate (p_Iedges,p_Itemp)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_replaceRowsMatrix9 (rmatrixSource,rmatrixDest,Irows)
  
!<description>
  ! Replaces all rows Irows in rmatrixDest by those in rmatrixSource.
  ! Source and destination matrix must have the same structure.
!</description>

!<input>
  ! Source matrix. Must be format 9!
  type(t_matrixScalar), intent(in) :: rmatrixSource

  ! Rows to replace
  integer, dimension(:), intent(in) :: Irows
!</input>

!<inputoutput>
  ! Destination matrix
  type(t_matrixScalar), intent(inout) :: rmatrixDest
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,irow
    integer, dimension(:), pointer :: p_Kld1, p_Kld2
    real(DP), dimension(:), pointer :: p_Da1, p_Da2
    
    ! Get pointers to the matrix data.
    call lsyssc_getbase_Kld (rmatrixSource,p_Kld1)
    call lsyssc_getbase_Kld (rmatrixDest,p_Kld2)
    call lsyssc_getbase_double (rmatrixSource,p_Da1)
    call lsyssc_getbase_double (rmatrixDest,p_Da2)
    
    ! Loop over all rows
    do irow = 1,size(Irows)
      if ((p_Kld1(irow) .ne. p_Kld2(irow)) .or. (p_Kld1(irow+1) .ne. p_Kld2(irow+1))) then
        call output_line("Source and destination matrix incompatible!",&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_replaceRowsMatrix9')
        call sys_halt()
      end if
      
      ! Copy the row.
      do i=p_Kld1(Irows(irow)),p_Kld1(Irows(irow)+1)-1
        p_Da2(i) = p_Da1(i)
      end do
    
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_replaceVectorEntries (rvectorSource,rvectorDest,Irows)
  
!<description>
  ! Replaces all entries Irows in rvectorDest by those in rvectorSource.
!</description>

!<input>
  ! Source vector. Must be format 9!
  type(t_vectorScalar), intent(in) :: rvectorSource

  ! Rows to replace
  integer, dimension(:), intent(in) :: Irows
!</input>

!<inputoutput>
  ! Destination vector
  type(t_vectorScalar), intent(inout) :: rvectorDest
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
    real(DP), dimension(:), pointer :: p_Da1, p_Da2
    
    ! Get pointers to the matrix data.
    call lsyssc_getbase_double (rvectorSource,p_Da1)
    call lsyssc_getbase_double (rvectorDest,p_Da2)
    
    ! Loop over all rows
    do i = 1,size(Irows)
      p_Da2(Irows(i)) = p_Da1(Irows(i))
    end do
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lsyssc_vectorLinearCombIndexed (rx,ry,cx,cy,Iidx)

!<description>
  ! Applies a linear combination of parts of two vectors:
  ! ry = cx * rx  +  cy * ry for all entries listed in Iidx
!</description>

  ! First source vector
  type(t_vectorScalar), intent(in) :: rx
  
  ! Scaling factor for Dx
  real(DP), intent(in) :: cx

  ! Scaling factor for Dy
  real(DP), intent(in) :: cy
  
  ! List of entries in rx/ry to process.
  integer, dimension(:), intent(in) :: Iidx
!</input>

!<inputoutput>
  ! Second source vector; also receives the result if rdest is not specified.
  type(t_vectorScalar), intent(inout), target :: ry
!</inputoutput>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dsource, p_Ddest
    integer :: i

    ! Get the pointers and copy the whole data array.
    call lsyssc_getbase_double(rx,p_Dsource)
    call lsyssc_getbase_double(ry,p_Ddest)
    
    do i=1,size(Iidx)
      p_Ddest(Iidx(i)) = cx * p_Dsource(Iidx(i)) + cy * p_Ddest(Iidx(i))
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fcoeff_dirichletbcc1 (rdiscretisationTrial,&
                  rdiscretisationTest, rform, nelements, npointsPerElement,&
                  Dpoints, ibct, DpointPar, IdofsTrial, IdofsTest,&
                  rdomainIntSubset, Dcoefficients, rcollection)
  
!<description>
  ! Calculates the operator
  !    $$ int_gamma c n grad (phi_j) phi_i $$
  ! Used for Dirichlet boundary control.
!</description>
  
  type(t_spatialDiscretisation), intent(in) :: rdiscretisationTrial
  type(t_spatialDiscretisation), intent(in) :: rdiscretisationTest
  type(t_bilinearForm), intent(in) :: rform
  integer, intent(in) :: nelements
  integer, intent(in) :: npointsPerElement
  real(DP), dimension(:,:,:), intent(in) :: Dpoints
  integer, intent(in) :: ibct
  real(DP), dimension(:,:), intent(in) :: DpointPar
  integer, dimension(:,:), intent(in) :: IdofsTrial
  integer, dimension(:,:), intent(in) :: IdofsTest
  type(t_domainIntSubset), intent(in) :: rdomainIntSubset

  type(t_collection), intent(inout), optional :: rcollection

  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients

!</subroutine>

    ! local variables
    integer :: ipt,iel
    real(DP), dimension(npointsPerElement,2) :: DtempVal
    type(t_boundary), pointer :: p_rboundary
    real(DP) :: dc
    
    p_rboundary => rdiscretisationTrial%p_rboundary
    dc = rcollection%DquickAccess(1)
    
    ! Evaluate the FEM functions in all the points.
    do iel=1,nelements
    
      ! Get the normal vector using the boundary object.
      call boundary_getNormalVec2D_mult(&
          p_rboundary, ibct, DpointPar(:,iel), &
          DtempVal(:,1),DtempVal(:,2),cparType=BDR_PAR_LENGTH)

      do ipt = 1,npointsPerElement
        ! c * n
        Dcoefficients(1,ipt,iel) = dc * DtempVal(ipt,1)
        Dcoefficients(2,ipt,iel) = dc * DtempVal(ipt,2)
      end do
      
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fcoeff_dirichletbcc2 (rdiscretisationTrial,&
                  rdiscretisationTest, rform, nelements, npointsPerElement,&
                  Dpoints, ibct, DpointPar, IdofsTrial, IdofsTest,&
                  rdomainIntSubset, Dcoefficients, rcollection)
  
!<description>
  ! Calculates the operator
  !    $$ int_gamma c n grad (phi_j) phi_i $$
  ! Used for Dirichlet boundary control.
!</description>
  
  type(t_spatialDiscretisation), intent(in) :: rdiscretisationTrial
  type(t_spatialDiscretisation), intent(in) :: rdiscretisationTest
  type(t_bilinearForm), intent(in) :: rform
  integer, intent(in) :: nelements
  integer, intent(in) :: npointsPerElement
  real(DP), dimension(:,:,:), intent(in) :: Dpoints
  integer, intent(in) :: ibct
  real(DP), dimension(:,:), intent(in) :: DpointPar
  integer, dimension(:,:), intent(in) :: IdofsTrial
  integer, dimension(:,:), intent(in) :: IdofsTest
  type(t_domainIntSubset), intent(in) :: rdomainIntSubset

  type(t_collection), intent(inout), optional :: rcollection

  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients

!</subroutine>

    ! local variables
    integer :: ipt,iel, idim
    real(DP), dimension(npointsPerElement,2) :: DtempVal
    type(t_boundary), pointer :: p_rboundary
    real(DP) :: dc
    
    p_rboundary => rdiscretisationTrial%p_rboundary
    dc = rcollection%DquickAccess(1)
    
    ! Dimension. 1=X, 2=Y-direction.
    idim = rcollection%IquickAccess(1)
    
    ! Evaluate the FEM functions in all the points.
    do iel=1,nelements
    
      ! Get the normal vector using the boundary object.
      call boundary_getNormalVec2D_mult(&
          p_rboundary, ibct, DpointPar(:,iel), &
          DtempVal(:,1),DtempVal(:,2), cparType=BDR_PAR_LENGTH)

      do ipt = 1,npointsPerElement
        ! c * n
        Dcoefficients(1,ipt,iel) = dc * DtempVal(ipt,idim)
      end do
      
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fcoeff_dirichletbcclf1 (rdiscretisation, rform, &
                  nelements, npointsPerElement, Dpoints, ibct, DpointPar, &
                  IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use collection
    use domainintegration
    use fsystem
    use scalarpde
    use spatialdiscretisation
    use triangulation
    
!<description>
  ! Calculates the operator
  !    $$ int_gamma c n grad (phi_j) phi_i $$
  ! Used for Dirichlet boundary control.
!</description>
    
  !<input>
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    type(t_linearForm), intent(in) :: rform
    integer, intent(in) :: nelements
    integer, intent(in) :: npointsPerElement
    real(DP), dimension(:,:,:), intent(in) :: Dpoints
    integer, intent(in) :: ibct
    real(DP), dimension(:,:), intent(in) :: DpointPar
    integer, dimension(:,:), intent(in) :: IdofsTest
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
  !</input>

  !<inputoutput>
    type(t_collection), intent(inout), optional :: rcollection
  !</inputoutput>
  
  !<output>
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    ! local variables
    integer :: ipt,iel
    real(DP), dimension(npointsPerElement,2) :: DtempVal
    type(t_boundary), pointer :: p_rboundary
    real(DP) :: dc
    
    p_rboundary => rdiscretisation%p_rboundary
    dc = rcollection%DquickAccess(1)
    
    ! Evaluate the FEM functions in all the points.
    do iel=1,nelements
    
      ! Get the normal vector using the boundary object.
      call boundary_getNormalVec2D_mult(&
          p_rboundary, ibct, DpointPar(:,iel), &
          DtempVal(:,1),DtempVal(:,2),cparType=BDR_PAR_LENGTH)

      do ipt = 1,npointsPerElement
        ! c * n
        Dcoefficients(1,ipt,iel) = dc * DtempVal(ipt,1)
        Dcoefficients(2,ipt,iel) = dc * DtempVal(ipt,2)
      end do
      
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fcoeff_dirichletbcclf2 (rdiscretisation, rform, &
                  nelements, npointsPerElement, Dpoints, ibct, DpointPar, &
                  IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use collection
    use domainintegration
    use fsystem
    use scalarpde
    use spatialdiscretisation
    use triangulation
    
!<description>
  ! Calculates the operator
  !    $$ int_gamma c n grad (phi_j) phi_i $$
  ! Used for Dirichlet boundary control.
!</description>
    
  !<input>
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    type(t_linearForm), intent(in) :: rform
    integer, intent(in) :: nelements
    integer, intent(in) :: npointsPerElement
    real(DP), dimension(:,:,:), intent(in) :: Dpoints
    integer, intent(in) :: ibct
    real(DP), dimension(:,:), intent(in) :: DpointPar
    integer, dimension(:,:), intent(in) :: IdofsTest
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
  !</input>

  !<inputoutput>
    type(t_collection), intent(inout), optional :: rcollection
  !</inputoutput>
  
  !<output>
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
  !</output>

!<subroutine>

    ! local variables
    integer :: ipt,iel, idim
    real(DP), dimension(npointsPerElement,2) :: DtempVal
    type(t_boundary), pointer :: p_rboundary
    real(DP) :: dc
    
    p_rboundary => rdiscretisation%p_rboundary
    dc = rcollection%DquickAccess(1)
    
    ! Dimension. 1=X, 2=Y-direction.
    idim = rcollection%IquickAccess(1)
    
    ! Evaluate the FEM functions in all the points.
    do iel=1,nelements
    
      ! Get the normal vector using the boundary object.
      call boundary_getNormalVec2D_mult(&
          p_rboundary, ibct, DpointPar(:,iel), &
          DtempVal(:,1),DtempVal(:,2), cparType=BDR_PAR_LENGTH)

      do ipt = 1,npointsPerElement
        ! c * n
        Dcoefficients(1,ipt,iel) = dc * DtempVal(ipt,idim)
      end do
      
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fcoeff_dirichletbcc_primal (rdiscretisationTrial,&
                  rdiscretisationTest, rform, nelements, npointsPerElement,&
                  Dpoints, ibct, DpointPar, IdofsTrial, IdofsTest,&
                  rdomainIntSubset, Dcoefficients, rcollection)
  
!<description>
  ! Calculates the operator
  !    $$ int_gamma c n grad (phi_j) phi_i $$
  ! Used for Dirichlet boundary control.
!</description>
  
  type(t_spatialDiscretisation), intent(in) :: rdiscretisationTrial
  type(t_spatialDiscretisation), intent(in) :: rdiscretisationTest
  type(t_bilinearForm), intent(in) :: rform
  integer, intent(in) :: nelements
  integer, intent(in) :: npointsPerElement
  real(DP), dimension(:,:,:), intent(in) :: Dpoints
  integer, intent(in) :: ibct
  real(DP), dimension(:,:), intent(in) :: DpointPar
  integer, dimension(:,:), intent(in) :: IdofsTrial
  integer, dimension(:,:), intent(in) :: IdofsTest
  type(t_domainIntSubset), intent(in) :: rdomainIntSubset

  type(t_collection), intent(inout), optional :: rcollection

  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients

!</subroutine>

    ! local variables
    integer :: i,j, iedge, iedgeidx, ivt1, ivt2
    real(DP) :: dedgelength
    integer, dimension(:), pointer :: p_IedgesAtBoundary
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(npointsPerElement,2) :: Dnormal
    
    real(DP) :: dc, dnu, dpenalty
    
    ! Get parameters from the collection
    dc = rcollection%DquickAccess(1)
    dnu = rcollection%DquickAccess(2)
    dpenalty = rcollection%DquickAccess(3)
    
    ! Access to the triangulation
    call storage_getbase_int (&
        rdiscretisationTrial%p_rtriangulation%h_IedgesAtBoundary,p_IedgesAtBoundary)
    call storage_getbase_int2d (&
        rdiscretisationTrial%p_rtriangulation%h_IverticesAtEdge,p_IverticesAtEdge)
    call storage_getbase_double2d (&
        rdiscretisationTrial%p_rtriangulation%h_DvertexCoords,p_DvertexCoords)
    
    do i=1,nelements
      do j=1,npointsPerElement
        ! For every parameter value, fetch the index of the edge in the edge list
        call tria_searchBoundaryEdgePar2D (ibct,DpointPar(j,i),&
            rdiscretisationTrial%p_rtriangulation,rdiscretisationTrial%p_rboundary,&
            iedgeidx,BDR_PAR_LENGTH)
        
        ! Get the actual edge number and adjacent vertices
        iedge = p_IedgesAtBoundary(iedgeidx)
        ivt1 = p_IverticesAtEdge(1,iedge)
        ivt2 = p_IverticesAtEdge(2,iedge)
        
        ! Calculate the length of the edge. this is our local h for the penalty term
        dedgelength = &
            sqrt ((p_DvertexCoords(1,ivt2)-p_DvertexCoords(1,ivt1))**2 + &
                  (p_DvertexCoords(2,ivt2)-p_DvertexCoords(2,ivt1))**2)
                  
        ! 1st coefficient: Penalty term
        Dcoefficients(1,j,i) = dc * dpenalty*dnu/dedgelength
      end do
      
      ! Get the normal vector using the boundary object.
      call boundary_getNormalVec2D_mult(&
          rdiscretisationTrial%p_rboundary, ibct, DpointPar(:,i), &
          Dnormal(:,1),Dnormal(:,2),cparType=BDR_PAR_LENGTH)

      ! 2nd/3rd coefficient: normal derivative in the test function
      do j = 1,npointsPerElement
        ! c * n
        Dcoefficients(2,j,i) = dc * (-dnu) * Dnormal(j,1)
        Dcoefficients(3,j,i) = dc * (-dnu) * Dnormal(j,2)
      end do
      
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fcoeff_dirichletbcc_dual (rdiscretisationTrial,&
                  rdiscretisationTest, rform, nelements, npointsPerElement,&
                  Dpoints, ibct, DpointPar, IdofsTrial, IdofsTest,&
                  rdomainIntSubset, Dcoefficients, rcollection)
  
!<description>
  ! Calculates the operator
  !    $$ int_gamma c n grad (phi_j) phi_i $$
  ! Used for Dirichlet boundary control.
!</description>
  
  type(t_spatialDiscretisation), intent(in) :: rdiscretisationTrial
  type(t_spatialDiscretisation), intent(in) :: rdiscretisationTest
  type(t_bilinearForm), intent(in) :: rform
  integer, intent(in) :: nelements
  integer, intent(in) :: npointsPerElement
  real(DP), dimension(:,:,:), intent(in) :: Dpoints
  integer, intent(in) :: ibct
  real(DP), dimension(:,:), intent(in) :: DpointPar
  integer, dimension(:,:), intent(in) :: IdofsTrial
  integer, dimension(:,:), intent(in) :: IdofsTest
  type(t_domainIntSubset), intent(in) :: rdomainIntSubset

  type(t_collection), intent(inout), optional :: rcollection

  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients

!</subroutine>

    ! local variables
    integer :: i,j, iedge, iedgeidx, ivt1, ivt2
    real(DP) :: dedgelength
    integer, dimension(:), pointer :: p_IedgesAtBoundary
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(npointsPerElement,2) :: Dnormal
    
    real(DP) :: dc, dnu, dpenalty
    
    ! Get parameters from the collection
    dc = rcollection%DquickAccess(1)
    dnu = rcollection%DquickAccess(2)
    dpenalty = rcollection%DquickAccess(3)
    
    ! Access to the triangulation
    call storage_getbase_int (&
        rdiscretisationTrial%p_rtriangulation%h_IedgesAtBoundary,p_IedgesAtBoundary)
    call storage_getbase_int2d (&
        rdiscretisationTrial%p_rtriangulation%h_IverticesAtEdge,p_IverticesAtEdge)
    call storage_getbase_double2d (&
        rdiscretisationTrial%p_rtriangulation%h_DvertexCoords,p_DvertexCoords)
    
    do i=1,nelements
    
      ! Get the normal vector using the boundary object.
      call boundary_getNormalVec2D_mult(&
          rdiscretisationTrial%p_rboundary, ibct, DpointPar(:,i), &
          Dnormal(:,1),Dnormal(:,2),cparType=BDR_PAR_LENGTH)

      do j=1,npointsPerElement
        ! For every parameter value, fetch the index of the edge in the edge list
        call tria_searchBoundaryEdgePar2D (ibct,DpointPar(j,i),&
            rdiscretisationTrial%p_rtriangulation,rdiscretisationTrial%p_rboundary,&
            iedgeidx,BDR_PAR_LENGTH)
        
        ! Get the actual edge number and adjacent vertices
        iedge = p_IedgesAtBoundary(iedgeidx)
        ivt1 = p_IverticesAtEdge(1,iedge)
        ivt2 = p_IverticesAtEdge(2,iedge)
        
        ! Calculate the length of the edge. this is our local h for the penalty term
        dedgelength = &
            sqrt ((p_DvertexCoords(1,ivt2)-p_DvertexCoords(1,ivt1))**2 + &
                  (p_DvertexCoords(2,ivt2)-p_DvertexCoords(2,ivt1))**2)
                  
        ! 1st coefficient: Penalty term
        Dcoefficients(1,j,i) = dc*Dnormal(j,1) * dpenalty*dnu/dedgelength
        Dcoefficients(4,j,i) = dc*Dnormal(j,2) * dpenalty*dnu/dedgelength
      end do
      
      ! 2nd/3rd coefficient: normal derivative in the test function
      do j = 1,npointsPerElement
        ! c * n
        Dcoefficients(2,j,i) = dc * Dnormal(j,1) * (-dnu) * Dnormal(j,1)
        Dcoefficients(3,j,i) = dc * Dnormal(j,1) * (-dnu) * Dnormal(j,2)
        Dcoefficients(5,j,i) = dc * Dnormal(j,2) * (-dnu) * Dnormal(j,1)
        Dcoefficients(6,j,i) = dc * Dnormal(j,2) * (-dnu) * Dnormal(j,2)
      end do
      
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fcoeff_dirichletbcc_xi (rdiscretisationTrial,&
                  rdiscretisationTest, rform, nelements, npointsPerElement,&
                  Dpoints, ibct, DpointPar, IdofsTrial, IdofsTest,&
                  rdomainIntSubset, Dcoefficients, rcollection)
  
!<description>
  ! Calculates the operator
  !    $$ int_gamma c n grad (phi_j) phi_i $$
  ! Used for Dirichlet boundary control.
!</description>
  
  type(t_spatialDiscretisation), intent(in) :: rdiscretisationTrial
  type(t_spatialDiscretisation), intent(in) :: rdiscretisationTest
  type(t_bilinearForm), intent(in) :: rform
  integer, intent(in) :: nelements
  integer, intent(in) :: npointsPerElement
  real(DP), dimension(:,:,:), intent(in) :: Dpoints
  integer, intent(in) :: ibct
  real(DP), dimension(:,:), intent(in) :: DpointPar
  integer, dimension(:,:), intent(in) :: IdofsTrial
  integer, dimension(:,:), intent(in) :: IdofsTest
  type(t_domainIntSubset), intent(in) :: rdomainIntSubset

  type(t_collection), intent(inout), optional :: rcollection

  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients

!</subroutine>

    ! local variables
    integer :: i,j, iedge, iedgeidx, ivt1, ivt2, idim
    real(DP) :: dedgelength
    integer, dimension(:), pointer :: p_IedgesAtBoundary
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(npointsPerElement,2) :: Dnormal
    
    real(DP) :: dc, dnu, dpenalty
    
    ! Get parameters from the collection
    dc = rcollection%DquickAccess(1)
    dnu = rcollection%DquickAccess(2)
    dpenalty = rcollection%DquickAccess(3)

    ! Dimension. 1=X, 2=Y-direction.
    idim = rcollection%IquickAccess(1)
    
    ! Access to the triangulation
    call storage_getbase_int (&
        rdiscretisationTrial%p_rtriangulation%h_IedgesAtBoundary,p_IedgesAtBoundary)
    call storage_getbase_int2d (&
        rdiscretisationTrial%p_rtriangulation%h_IverticesAtEdge,p_IverticesAtEdge)
    call storage_getbase_double2d (&
        rdiscretisationTrial%p_rtriangulation%h_DvertexCoords,p_DvertexCoords)
    
    do i=1,nelements

      ! Get the normal vector using the boundary object.
      call boundary_getNormalVec2D_mult(&
          rdiscretisationTrial%p_rboundary, ibct, DpointPar(:,i), &
          Dnormal(:,1),Dnormal(:,2),cparType=BDR_PAR_LENGTH)

      do j=1,npointsPerElement
        ! For every parameter value, fetch the index of the edge in the edge list
        call tria_searchBoundaryEdgePar2D (ibct,DpointPar(j,i),&
            rdiscretisationTrial%p_rtriangulation,rdiscretisationTrial%p_rboundary,&
            iedgeidx,BDR_PAR_LENGTH)
        
        ! Get the actual edge number and adjacent vertices
        iedge = p_IedgesAtBoundary(iedgeidx)
        ivt1 = p_IverticesAtEdge(1,iedge)
        ivt2 = p_IverticesAtEdge(2,iedge)
        
        ! Calculate the length of the edge. this is our local h for the penalty term
        dedgelength = &
            sqrt ((p_DvertexCoords(1,ivt2)-p_DvertexCoords(1,ivt1))**2 + &
                  (p_DvertexCoords(2,ivt2)-p_DvertexCoords(2,ivt1))**2)
                  
        ! 1st coefficient: Penalty term
        Dcoefficients(1,j,i) = dc*Dnormal(j,idim) * dpenalty*dnu/dedgelength
      
        ! 2nd/3rd coefficient: normal derivative in the test function
        
        ! c * n
        Dcoefficients(2,j,i) = dc*Dnormal(j,idim) * (-dnu) * Dnormal(j,1)
        Dcoefficients(3,j,i) = dc*Dnormal(j,idim) * (-dnu) * Dnormal(j,2)
      end do
      
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine smva_assembleDirichletBCC (rphysics,rmatrixY1,rmatrixY2,&
      rmatrixLambda1,rmatrixLambda2,rmatrixGrad1,rmatrixGrad2,dc1,dc2,dc3,&
      dpenalty,rdirichletBCC)

!<description>
  ! This routine assembles the term
  !    $$ int_gamma dc1 y + dc1 n grad(lambda) + dc2 n xi $$
  ! on all boundary edges.
  ! This term is crucial for the Dirichlet bonudary control.
!</description>

!<input>
  ! Underlying physics
  type(t_settings_physics) :: rphysics

  ! Multiplier for the mass matrix
  real(DP), intent(in) :: dc1

  ! Multiplier for the lambda-matrix
  real(DP), intent(in) :: dc2

  ! Multiplier for the gradient matrix
  real(DP), intent(in) :: dc3
  
  ! Penalty parameter for the Dirichlet boundary conditions
  real(DP), intent(in) :: dpenalty
  
  ! Defines the Boundary segments where to assemble
  type(t_boundaryRegionList), intent(in) :: rdirichletBCC
!</input>

!<inputoutput>
  ! Primal velocity destination matrices
  type(t_matrixScalar), intent(inout), target :: rmatrixY1,rmatrixY2

  ! Dual velocity destination matrices
  type(t_matrixScalar), intent(inout), target :: rmatrixLambda1,rmatrixLambda2

  ! Gradient destination matrices (pressure)
  type(t_matrixScalar), intent(inout), target :: rmatrixGrad1,rmatrixGrad2
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer :: i
    type(t_bilinearForm) :: rform1,rform2,rform3
    type(t_bilinearForm) :: rformweak1,rformweak2,rformweak3
    type(t_linearForm) :: rlinform1,rlinform2,rlinform3
    type(t_bdRegionEntry), pointer :: p_rdirichletBCCBd
    type(t_collection) :: rcollection
    type(t_matrixScalar) :: rmatrixY1temp,rmatrixY2temp
    type(t_matrixScalar) :: rmatrixLambda1temp,rmatrixLambda2temp
    type(t_matrixScalar) :: rmatrixGrad1temp,rmatrixGrad2temp
    integer, dimension(:), allocatable :: Idofs
    integer :: ndofs
    
    ! Cancel if there is nothing to assemble.
    if ((dc1 .eq. 0.0_DP) .and. (dc2 .eq. 0.0_DP) .and. (dc3 .eq. 0.0_DP)) then
      return
    end if
    
    ! If any of the matrices is deactivated by its scaliing factor,
    ! activate it and clear it in advance.
    if (rmatrixY1%dscaleFactor .eq. 0.0_DP) then
      rmatrixY1%dscaleFactor = 1.0_DP
      call lsyssc_clearMatrix (rmatrixY1)
    end if

    if (rmatrixY2%dscaleFactor .eq. 0.0_DP) then
      rmatrixY2%dscaleFactor = 1.0_DP
      call lsyssc_clearMatrix (rmatrixY2)
    end if

    if (rmatrixLambda1%dscaleFactor .eq. 0.0_DP) then
      rmatrixLambda1%dscaleFactor = 1.0_DP
      call lsyssc_clearMatrix (rmatrixLambda1)
    end if

    if (rmatrixLambda2%dscaleFactor .eq. 0.0_DP) then
      rmatrixLambda2%dscaleFactor = 1.0_DP
      call lsyssc_clearMatrix (rmatrixLambda2)
    end if

    if (rmatrixGrad1%dscaleFactor .eq. 0.0_DP) then
      rmatrixGrad1%dscaleFactor = 1.0_DP
      call lsyssc_clearMatrix (rmatrixGrad1)
    end if

    if (rmatrixGrad2%dscaleFactor .eq. 0.0_DP) then
      rmatrixGrad2%dscaleFactor = 1.0_DP
      call lsyssc_clearMatrix (rmatrixGrad2)
    end if
    
!    ! Create temporary copies.
!    call lsyssc_duplicateMatrix (rmatrixY1,rmatrixY1temp,LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!    call lsyssc_duplicateMatrix (rmatrixY2,rmatrixY2temp,LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!    call lsyssc_duplicateMatrix (rmatrixLambda1,rmatrixLambda1temp,LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!    call lsyssc_duplicateMatrix (rmatrixLambda2,rmatrixLambda2temp,LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!    call lsyssc_duplicateMatrix (rmatrixGrad1,rmatrixGrad1temp,LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!    call lsyssc_duplicateMatrix (rmatrixGrad2,rmatrixGrad2temp,LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
    
    ! Assemble the matrices
    
    ! Overwrite all rows by zero which do not correspond to
    ! DOF's on the boundary (HACK!!!)
    
    ! Prepare a bilinear form for y
    rform1%itermCount = 1
    rform1%Idescriptors(1,1) = DER_FUNC2D
    rform1%Idescriptors(2,1) = DER_FUNC2D
    rform1%ballCoeffConstant = .true.
    rform1%BconstantCoeff(1) = .true.
    rform1%Dcoefficients(1:rform1%itermCount)  = dc1 * dpenalty

    ! Prepare a bilinear form for lambda
    rform2%itermCount = 2
    rform2%Idescriptors(1,1) = DER_DERIV2D_X
    rform2%Idescriptors(2,1) = DER_FUNC2D
    rform2%Idescriptors(1,2) = DER_DERIV2D_Y
    rform2%Idescriptors(2,2) = DER_FUNC2D
    rform2%ballCoeffConstant = .false.
    rform2%BconstantCoeff(1:rform2%itermCount) = .false.
    rform2%Dcoefficients(1:rform2%itermCount)  = 1.0_DP
    
    ! Prepare a bilinear form for p
    rform3%itermCount = 1
    rform3%Idescriptors(1,1) = DER_FUNC2D
    rform3%Idescriptors(2,1) = DER_FUNC2D
    rform3%ballCoeffConstant = .false.
    rform3%BconstantCoeff(1:rform3%itermCount) = .false.
    rform3%Dcoefficients(1:rform3%itermCount)  = 1.0_DP

    ! Prepare a linear form for y
    rlinform1%itermCount = 1
    rlinform1%Idescriptors(1) = DER_FUNC2D
    rlinform1%Dcoefficients(1:rlinform1%itermCount)  = dc1*dpenalty

    ! Prepare a linear form for lambda
    rlinform2%itermCount = 2
    rlinform2%Idescriptors(1) = DER_DERIV2D_X
    rlinform2%Idescriptors(2) = DER_DERIV2D_Y
    rlinform2%Dcoefficients(1:rlinform2%itermCount)  = 1.0_DP
   
    ! Prepare a linear form for p
    rlinform3%itermCount = 1
    rlinform3%Idescriptors(1) = DER_FUNC2D
    rlinform3%Dcoefficients(1:rform3%itermCount)  = 1.0_DP

    ! Prepare a bilinear form for y
    rformweak1%itermCount = 3
    rformweak1%Idescriptors(1,1) = DER_FUNC2D
    rformweak1%Idescriptors(2,1) = DER_FUNC2D
    rformweak1%Idescriptors(1,2) = DER_FUNC2D
    rformweak1%Idescriptors(2,2) = DER_DERIV2D_X
    rformweak1%Idescriptors(1,3) = DER_FUNC2D
    rformweak1%Idescriptors(2,3) = DER_DERIV2D_Y
    rformweak1%ballCoeffConstant = .false.
    rformweak1%BconstantCoeff(1:rformweak1%itermCount) = .false.
    rformweak1%Dcoefficients(1:rformweak1%itermCount)  = 1.0_DP

    ! Prepare a bilinear form for lambda
    rformweak2%itermCount = 6
    rformweak2%Idescriptors(1,1) = DER_DERIV2D_X
    rformweak2%Idescriptors(2,1) = DER_FUNC2D
    
    rformweak2%Idescriptors(1,2) = DER_DERIV2D_X
    rformweak2%Idescriptors(2,2) = DER_DERIV2D_X
    
    rformweak2%Idescriptors(1,3) = DER_DERIV2D_X
    rformweak2%Idescriptors(2,3) = DER_DERIV2D_Y
    
    rformweak2%Idescriptors(1,4) = DER_DERIV2D_Y
    rformweak2%Idescriptors(2,4) = DER_FUNC2D
    
    rformweak2%Idescriptors(1,5) = DER_DERIV2D_Y
    rformweak2%Idescriptors(2,5) = DER_DERIV2D_X
    
    rformweak2%Idescriptors(1,6) = DER_DERIV2D_Y
    rformweak2%Idescriptors(2,6) = DER_DERIV2D_Y
    
    rformweak2%ballCoeffConstant = .false.
    rformweak2%BconstantCoeff(1:rformweak2%itermCount) = .false.
    rformweak2%Dcoefficients(1:rformweak2%itermCount)  = 1.0_DP

    ! Prepare a bilinear form for p
    rformweak3%itermCount = 3
    rformweak3%Idescriptors(1,1) = DER_FUNC2D
    rformweak3%Idescriptors(2,1) = DER_FUNC2D
    rformweak3%Idescriptors(1,2) = DER_FUNC2D
    rformweak3%Idescriptors(2,2) = DER_DERIV2D_X
    rformweak3%Idescriptors(1,3) = DER_FUNC2D
    rformweak3%Idescriptors(2,3) = DER_DERIV2D_Y
    rformweak3%ballCoeffConstant = .false.
    rformweak3%BconstantCoeff(1) = .false.
    rformweak3%Dcoefficients(1:rformweak3%itermCount)  = 1.0_DP

    ! Assemble the operator on all boundary components in rneumannBoundary.
    p_rdirichletBCCBd => rdirichletBCC%p_rprimalBdHead
    
    do i = 1,rdirichletBCC%nregionsPrimal

!      ! Get the DOF's in that segment
!      call bcasm_getDOFsInBDRegion (rmatrixY1Temp%p_rspatialDiscrTrial,&
!          p_rdirichletBCCBd%rboundaryRegion,ndofs=ndofs)
!      
!      ! Allocate memory if not done already.
!      if (.not. allocated(Idofs)) allocate (Idofs(ndofs))
!      if (size (Idofs) .lt. ndofs) then
!        deallocate (Idofs)
!        allocate (Idofs(ndofs))
!      end if
!
!      ! Calculate the DOF`s which are affected.          
!      call bcasm_getDOFsInBDRegion (rmatrixY1Temp%p_rspatialDiscrTrial,&
!          p_rdirichletBCCBd%rboundaryRegion,IdofsArray=Idofs)

!#ifdef DEBUG
!      ! Clear the matrices
!      call lsyssc_clearMatrix (rmatrixY1temp)
!      call lsyssc_clearMatrix (rmatrixY2temp)
!      call lsyssc_clearMatrix (rmatrixLambda1temp)
!      call lsyssc_clearMatrix (rmatrixLambda2temp)
!      call lsyssc_clearMatrix (rmatrixGrad1temp)
!      call lsyssc_clearMatrix (rmatrixGrad2temp)
!#else
!      ! Clear the rows in the matrix which are later added to the global matrix.
!      call mmod_replaceLinesByZero (rmatrixY1Temp,Idofs(1:ndofs))
!      call mmod_replaceLinesByZero (rmatrixY2Temp,Idofs(1:ndofs))
!      call mmod_replaceLinesByZero (rmatrixLambda1Temp,Idofs(1:ndofs))
!      call mmod_replaceLinesByZero (rmatrixLambda2Temp,Idofs(1:ndofs))
!      call mmod_replaceLinesByZero (rmatrixGrad1Temp,Idofs(1:ndofs))
!      call mmod_replaceLinesByZero (rmatrixGrad2Temp,Idofs(1:ndofs))
!#endif

      ! ------------------------------------------------------------
      ! Direct weak implementation of the BCC, no temporary matrices
      ! ------------------------------------------------------------
      
      ! Set up the matrix. Primal velocity.
      
      call bilf_buildMatrixScalarBdr2D (rform1, CUB_G4_1D, .false., &
          rmatrixY1,rboundaryRegion=p_rdirichletBCCBd%rboundaryRegion)
      
      if (.not. lsyssc_isMatrixContentShared(rmatrixY1,rmatrixY2)) then
        ! Second matrix if not identical to the first one.
        call bilf_buildMatrixScalarBdr2D (rform1, CUB_G4_1D, .false., &
            rmatrixY2,rboundaryRegion=p_rdirichletBCCBd%rboundaryRegion)
      end if
      
      ! Dual velocity.

      rcollection%DquickAccess(1) = dc2 * dpenalty

      call bilf_buildMatrixScalarBdr2D (rform2, CUB_G4_1D, .false., &
          rmatrixLambda1,fcoeff_dirichletbcc1,rboundaryRegion=p_rdirichletBCCBd%rboundaryRegion,&
          rcollection=rcollection)

      if (.not. lsyssc_isMatrixContentShared(rmatrixLambda1,rmatrixLambda2)) then
        ! Second matrix if not identical to the first one.
        call bilf_buildMatrixScalarBdr2D (rform2, CUB_G4_1D, .false., &
            rmatrixLambda2,fcoeff_dirichletbcc1,rboundaryRegion=p_rdirichletBCCBd%rboundaryRegion,&
            rcollection=rcollection)
      end if

      ! Pressure
      rcollection%DquickAccess(1) = dc3 * dpenalty

      ! X-derivative
      rcollection%IquickAccess(1) = 1

      call bilf_buildMatrixScalarBdr2D (rform3, CUB_G4_1D, .false., &
          rmatrixGrad1,fcoeff_dirichletbcc2,p_rdirichletBCCBd%rboundaryRegion,rcollection)

      ! Y-derivative
      rcollection%IquickAccess(1) = 2

      call bilf_buildMatrixScalarBdr2D (rform3, CUB_G4_1D, .false., &
          rmatrixGrad2,fcoeff_dirichletbcc2,p_rdirichletBCCBd%rboundaryRegion,rcollection)

!      ! ------------------------------
!      ! Weak implementation of the BCC
!      ! ------------------------------
!      
!      ! Set up the matrix. Primal velocity.
!      
!      call bilf_buildMatrixScalarBdr2D (rform1, CUB_G4_1D, .true., &
!          rmatrixY1Temp,rboundaryRegion=p_rdirichletBCCBd%rboundaryRegion)
!      
!      if (.not. lsyssc_isMatrixContentShared(rmatrixY1,rmatrixY2)) then
!        ! Second matrix if not identical to the first one.
!        call bilf_buildMatrixScalarBdr2D (rform1, CUB_G4_1D, .true., &
!            rmatrixY2Temp,rboundaryRegion=p_rdirichletBCCBd%rboundaryRegion)
!      end if
!      
!      ! Dual velocity.
!
!      rcollection%DquickAccess(1) = dc2 * dpenalty
!
!      call bilf_buildMatrixScalarBdr2D (rform2, CUB_G4_1D, .true., &
!          rmatrixLambda1Temp,fcoeff_dirichletbcc1,rboundaryRegion=p_rdirichletBCCBd%rboundaryRegion,&
!          rcollection=rcollection)
!
!      if (.not. lsyssc_isMatrixContentShared(rmatrixLambda1,rmatrixLambda2)) then
!        ! Second matrix if not identical to the first one.
!        call bilf_buildMatrixScalarBdr2D (rform2, CUB_G4_1D, .true., &
!            rmatrixLambda2Temp,fcoeff_dirichletbcc1,rboundaryRegion=p_rdirichletBCCBd%rboundaryRegion,&
!            rcollection=rcollection)
!      end if
!
!      ! Pressure
!      rcollection%DquickAccess(1) = dc3 * dpenalty
!
!      ! X-derivative
!      rcollection%IquickAccess(1) = 1
!
!      call bilf_buildMatrixScalarBdr2D (rform3, CUB_G4_1D, .true., &
!          rmatrixGrad1Temp,fcoeff_dirichletbcc2,p_rdirichletBCCBd%rboundaryRegion,rcollection)
!
!      ! Y-derivative
!      rcollection%IquickAccess(1) = 2
!
!      call bilf_buildMatrixScalarBdr2D (rform3, CUB_G4_1D, .true., &
!          rmatrixGrad2Temp,fcoeff_dirichletbcc2,p_rdirichletBCCBd%rboundaryRegion,rcollection)


!      ! --------------------------------
!      ! Strong implementation of the BCC
!      ! --------------------------------
!
!      ! Clear the rows in the matrix which are affected.
!!      call mmod_replaceLinesByZero (rmatrixY1Temp,Idofs(1:ndofs))
!!      call mmod_replaceLinesByZero (rmatrixY2Temp,Idofs(1:ndofs))
!!      call mmod_replaceLinesByZero (rmatrixLambda1Temp,Idofs(1:ndofs))
!!      call mmod_replaceLinesByZero (rmatrixLambda2Temp,Idofs(1:ndofs))
!!      call mmod_replaceLinesByZero (rmatrixGrad1Temp,Idofs(1:ndofs))
!!      call mmod_replaceLinesByZero (rmatrixGrad2Temp,Idofs(1:ndofs))
!
!      ! Set up the matrix. Primal velocity.
!
!      call linf_getBoundaryOperatorMatrix (&
!          rlinform1,rmatrixY1Temp,.false.,p_rdirichletBCCBd%rboundaryRegion)
!
!      if (.not. lsyssc_isMatrixContentShared(rmatrixY1,rmatrixY2)) then
!        ! Second matrix if not identical to the first one.
!        call linf_getBoundaryOperatorMatrix (&
!            rlinform1,rmatrixY2Temp,.false.,p_rdirichletBCCBd%rboundaryRegion)
!      end if
!
!      ! Dual velocity.
!
!      rcollection%DquickAccess(1) = dc2*dpenalty
!
!      call linf_getBoundaryOperatorMatrix (&
!          rlinform2,rmatrixLambda1Temp,.false.,p_rdirichletBCCBd%rboundaryRegion,&
!          fcoeff_dirichletbcclf1,rcollection)
!
!      if (.not. lsyssc_isMatrixContentShared(rmatrixLambda1,rmatrixLambda2)) then
!        ! Second matrix if not identical to the first one.
!        call linf_getBoundaryOperatorMatrix (&
!            rlinform2,rmatrixLambda2Temp,.false.,p_rdirichletBCCBd%rboundaryRegion,&
!            fcoeff_dirichletbcclf1,rcollection)
!      end if
!
!      ! Pressure
!      rcollection%DquickAccess(1) = dc3*dpenalty
!
!      ! X-derivative
!      rcollection%IquickAccess(1) = 1
!
!      call linf_getBoundaryOperatorMatrix (&
!          rlinform3,rmatrixGrad1Temp,.false.,p_rdirichletBCCBd%rboundaryRegion,&
!          fcoeff_dirichletbcclf2,rcollection)
!
!      ! Y-derivative
!      rcollection%IquickAccess(1) = 2
!
!      call linf_getBoundaryOperatorMatrix (&
!          rlinform3,rmatrixGrad2Temp,.false.,p_rdirichletBCCBd%rboundaryRegion,&
!          fcoeff_dirichletbcclf2,rcollection)

!      ! ------------------------------------------
!      ! Alternative Weak implementation of the BCC
!      ! ------------------------------------------
!      
!      ! Set up the matrix. Primal velocity.
!      
!      rcollection%DquickAccess(1) = dc1
!      rcollection%DquickAccess(2) = rphysics%dnu
!      rcollection%DquickAccess(3) = dpenalty
!      
!      call bilf_buildMatrixScalarBdr2D (rformweak1, CUB_G4_1D, .false., &
!          rmatrixY1Temp,fcoeff_dirichletbcc_primal,&
!          rboundaryRegion=p_rdirichletBCCBd%rboundaryRegion,&
!          rcollection=rcollection)
!
!      if (.not. lsyssc_isMatrixContentShared(rmatrixY1,rmatrixY2)) then
!        call bilf_buildMatrixScalarBdr2D (rformweak1, CUB_G4_1D, .false., &
!            rmatrixY2Temp,fcoeff_dirichletbcc_primal,&
!            rboundaryRegion=p_rdirichletBCCBd%rboundaryRegion,&
!            rcollection=rcollection)
!      end if
!
!      ! Set up the matrix. Dual velocity.
!      
!      rcollection%DquickAccess(1) = dc2
!      rcollection%DquickAccess(2) = rphysics%dnu
!      
!      call bilf_buildMatrixScalarBdr2D (rformweak2, CUB_G4_1D, .false., &
!          rmatrixLambda1Temp,fcoeff_dirichletbcc_dual,&
!          rboundaryRegion=p_rdirichletBCCBd%rboundaryRegion,&
!          rcollection=rcollection)
!
!      if (.not. lsyssc_isMatrixContentShared(rmatrixLambda1,rmatrixLambda2)) then
!        call bilf_buildMatrixScalarBdr2D (rformweak2, CUB_G4_1D, .false., &
!            rmatrixLambda2Temp,fcoeff_dirichletbcc_dual,&
!            rboundaryRegion=p_rdirichletBCCBd%rboundaryRegion,&
!            rcollection=rcollection)
!      end if
!
!      ! Pressure
!      rcollection%DquickAccess(1) = dc3
!
!      ! X-derivative
!      rcollection%IquickAccess(1) = 1
!
!      call bilf_buildMatrixScalarBdr2D (rformweak3, CUB_G4_1D, .false., &
!          rmatrixGrad1Temp,fcoeff_dirichletbcc_xi,p_rdirichletBCCBd%rboundaryRegion,rcollection)
!
!      ! Y-derivative
!      rcollection%IquickAccess(1) = 2
!
!      call bilf_buildMatrixScalarBdr2D (rformweak3, CUB_G4_1D, .false., &
!          rmatrixGrad2Temp,fcoeff_dirichletbcc_xi,p_rdirichletBCCBd%rboundaryRegion,rcollection)

      ! ----------------------------------------
      ! Imponsing of the control to the matrices
      ! ----------------------------------------

!      ! Sum up the matrices
!
!      ! Only the entries on the boundary, replacing the old entries
!      call lsyssc_matrixLinearCombIndexed (&
!          rmatrixY1temp,rmatrixY1,1.0_DP,0.0_DP,Idofs(1:ndofs))
!      if (.not. lsyssc_isMatrixContentShared(rmatrixY1,rmatrixY2)) then
!        call lsyssc_matrixLinearCombIndexed (&
!            rmatrixY2temp,rmatrixY2,1.0_DP,0.0_DP,Idofs(1:ndofs))
!      end if
!      
!      call lsyssc_matrixLinearCombIndexed (&
!          rmatrixLambda1Temp,rmatrixLambda1,1.0_DP,0.0_DP,Idofs(1:ndofs))
!      
!      if (.not. lsyssc_isMatrixContentShared(rmatrixLambda1,rmatrixLambda2)) then
!        call lsyssc_matrixLinearCombIndexed (&
!            rmatrixLambda2temp,rmatrixLambda2,1.0_DP,0.0_DP,Idofs(1:ndofs))
!      end if
!      
!      call lsyssc_matrixLinearCombIndexed (&
!          rmatrixGrad1temp,rmatrixGrad1,1.0_DP,0.0_DP,Idofs(1:ndofs))
!      call lsyssc_matrixLinearCombIndexed (&
!          rmatrixGrad2temp,rmatrixGrad2,1.0_DP,0.0_DP,Idofs(1:ndofs))
          
!      ! Only the entries on the boundary, summing up to the old entries
!      call lsyssc_matrixLinearCombIndexed (&
!          rmatrixY1temp,rmatrixY1,1.0_DP,1.0_DP,Idofs(1:ndofs))
!      if (.not. lsyssc_isMatrixContentShared(rmatrixY1,rmatrixY2)) then
!        call lsyssc_matrixLinearCombIndexed (&
!            rmatrixY2temp,rmatrixY2,1.0_DP,1.0_DP,Idofs(1:ndofs))
!      end if
!      
!      call lsyssc_matrixLinearCombIndexed (&
!          rmatrixLambda1Temp,rmatrixLambda1,1.0_DP,1.0_DP,Idofs(1:ndofs))
!      if (.not. lsyssc_isMatrixContentShared(rmatrixLambda1,rmatrixLambda2)) then
!        call lsyssc_matrixLinearCombIndexed (&
!            rmatrixLambda2temp,rmatrixLambda2,1.0_DP,1.0_DP,Idofs(1:ndofs))
!      end if
!      
!      call lsyssc_matrixLinearCombIndexed (&
!          rmatrixGrad1temp,rmatrixGrad1,1.0_DP,1.0_DP,Idofs(1:ndofs))
!      call lsyssc_matrixLinearCombIndexed (&
!          rmatrixGrad2temp,rmatrixGrad2,1.0_DP,1.0_DP,Idofs(1:ndofs))

!      ! All entries from the integral, summing up
!      call lsyssc_matrixLinearCombIndexed (&
!          rmatrixY1temp,rmatrixY1,1.0_DP,1.0_DP)
!      if (.not. lsyssc_isMatrixContentShared(rmatrixY1,rmatrixY2)) then
!        call lsyssc_matrixLinearCombIndexed (&
!            rmatrixY2temp,rmatrixY2,1.0_DP,1.0_DP)
!      end if
!      
!      call lsyssc_matrixLinearCombIndexed (&
!          rmatrixLambda1Temp,rmatrixLambda1,1.0_DP,1.0_DP)
!      if (.not. lsyssc_isMatrixContentShared(rmatrixLambda1,rmatrixLambda2)) then
!        call lsyssc_matrixLinearCombIndexed (&
!            rmatrixLambda2temp,rmatrixLambda2,1.0_DP,1.0_DP)
!      end if
!      
!      call lsyssc_matrixLinearCombIndexed (&
!          rmatrixGrad1temp,rmatrixGrad1,1.0_DP,1.0_DP)
!      call lsyssc_matrixLinearCombIndexed (&
!          rmatrixGrad2temp,rmatrixGrad2,1.0_DP,1.0_DP)


      ! Next segment
      p_rdirichletBCCBd => p_rdirichletBCCBd%p_nextBdRegion
    end do

!    ! Release memory
!    if (allocated(Idofs)) deallocate (Idofs)

!    ! Release the temporary copies.
!    call lsyssc_releaseMatrix (rmatrixY1temp)
!    call lsyssc_releaseMatrix (rmatrixY2temp)
!    call lsyssc_releaseMatrix (rmatrixLambda1temp)
!    call lsyssc_releaseMatrix (rmatrixLambda2temp)
!    call lsyssc_releaseMatrix (rmatrixGrad1temp)
!    call lsyssc_releaseMatrix (rmatrixGrad2temp)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lsyssc_matrixLinearCombIndexed (rmatrixA,rmatrixB,ca,cb,Irows)

    !<description>
    ! Performs a linear combination:
    !   rmatrixB = ca*rmatrixA + cB*rmatrixB
    !
    ! Only those rows indexed by Irows are processed. 
    ! Both matrices must have the same matrix format and structure.
    ! </description>

!<input>
    ! Source matrix
    type(t_matrixScalar), intent(in) :: rmatrixA

    ! scaling factors
    real(DP), intent(in) :: ca,cb

    ! OPTIONAL: List of rows to be processed.
    ! If not specified, all rows are processed.
    integer, dimension(:), intent(in), optional :: Irows
!</input>

!<inputoutput>
    ! OPTIONAL: Destination matrix.
    type(t_matrixScalar), intent(inout), target, optional :: rmatrixB
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata1, p_Ddata2
    integer, dimension(:), pointer :: p_Kld
    integer :: i,j

    ! Check if all matrices are compatible
    if ((rmatrixA%NEQ   .ne. rmatrixB%NEQ) .or.&
        (rmatrixA%NCOLS .ne. rmatrixB%NCOLS)) then
      call output_line("Number of rows/columns is not compatible!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_matrixLinearCombIndexed")
      call sys_halt()
    end if

    ! Matrices must have the same formatcmatrixFormat
    if (rmatrixA%cmatrixFormat   .ne. rmatrixB%cmatrixFormat) then
      call output_line("Matrices have different format!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_matrixLinearCombIndexed")
      call sys_halt()
    end if

    ! Check if all matrices have the same sorting
    if (rmatrixA%isortStrategy .ne. rmatrixB%isortStrategy) then
      call output_line("Sorting strategies are incompatible!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_matrixLinearCombIndexed")
      call sys_halt()
    end if

    ! Check if destination matrix has compatible data type
    if ((rmatrixA%cdatatype .ne. ST_DOUBLE) .or.&
        (rmatrixB%cdatatype .ne. ST_DOUBLE)) then
      call output_line("Data type of destination matrix not supported!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_matrixLinearCombIndexed")
      call sys_halt()
    end if
    
    select case (rmatrixA%cmatrixFormat)
    case (LSYSSC_MATRIX9)

      select case (rmatrixA%cmatrixFormat)
      case (LSYSSC_MATRIX9)
      
        ! Get the data pointers
        call lsyssc_getbase_double (rmatrixA,p_Ddata1)
        call lsyssc_getbase_double (rmatrixB,p_Ddata2)
        call lsyssc_getbase_Kld (rmatrixA,p_Kld)
        
        ! Linear combination of the rows
        if (present(Irows)) then
          ! Subset of the rows
          if ((ca .eq. 0.0_DP) .and. (cb .eq. 0.0_DP)) then
            do i=1,size(Irows)
              do j=p_Kld(Irows(i)), p_Kld(Irows(i+1)-1)
                p_Ddata2(j) = 0.0_DP
              end do
            end do
          else if (ca .eq. 0.0_DP) then
            do i=1,size(Irows)
              do j=p_Kld(Irows(i)), p_Kld(Irows(i)+1)-1
                p_Ddata2(j) = cb * p_Ddata2(j)
              end do
            end do
          else if (cb .eq. 0.0_DP) then
            do i=1,size(Irows)
              do j=p_Kld(Irows(i)), p_Kld(Irows(i)+1)-1
                p_Ddata2(j) = ca * p_Ddata1(j)
              end do
            end do
          else
            do i=1,size(Irows)
              do j=p_Kld(Irows(i)), p_Kld(Irows(i)+1)-1
                p_Ddata2(j) = ca * p_Ddata1(j) + cb * p_Ddata2(j)
              end do
            end do
          end if
        
        else
        
          ! All rows
          if ((ca .eq. 0.0_DP) .and. (cb .eq. 0.0_DP)) then
            do j=1,rmatrixA%NA
              p_Ddata2(j) = 0.0_DP
            end do
          else if (ca .eq. 0.0_DP) then
            do j=1,rmatrixA%NA
              p_Ddata2(j) = cb * p_Ddata2(j)
            end do
          else if (cb .eq. 0.0_DP) then
            do j=1,rmatrixA%NA
              p_Ddata2(j) = ca * p_Ddata1(j)
            end do
          else
            do j=1,rmatrixA%NA
              p_Ddata2(j) = ca * p_Ddata1(j) + cb * p_Ddata2(j)
            end do
          end if
        end if        
      
      case default
        call output_line(&
            "Combination of source and destination matrix format not supported!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_matrixLinearCombIndexed")
        call sys_halt()
      end select

    case default
      call output_line(&
          "Combination of source and destination matrix format not supported!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_matrixLinearCombIndexed")
      call sys_halt()
    end select      

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine smva_clearMatrix (rnonlinearSpatialMatrix)

!<description>
  ! Resets all operators in a nonlinear matrix to zero.
!</description>

!<inputoutput>
  ! Matrix to clear.
  type(t_nonlinearSpatialMatrix), intent(inout) :: rnonlinearSpatialMatrix
!</inputoutput>
  
!</subroutine>

    rnonlinearSpatialMatrix%DidentityY(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%Dmass(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%Dstokes(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%Dygrad(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%DygradAdj(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%Dgrady(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%DgradyAdj(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%DygradT(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%DygradTAdj(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%Dgrady2(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%DgradyAdj2(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%DygradT2(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%DygradTAdj2(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%DgradyT(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%DBmat(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%DBTmat(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%DidentityP(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%DdirichletBCCY(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%DdirichletBCCLambda(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%DdirichletBCCXi(:,:) = 0.0_DP

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_disableSubmatrix (rnonlinearSpatialMatrix,irow,icolumn)

!<description>
  ! Disables a subbklock in the nonlinear matrix rnonlinearSpatialMatrix.
  ! All weights of the correspopnding subblock are set to 0.
!</description>

!<input>
  ! The row/column of the submatrix to be disabled.
  integer :: irow,icolumn
!</input>

!<inputoutput>
  ! A t_nonlinearSpatialMatrix structure that defines the shape of the core
  ! equation. The weights that specify the submatrices of a small 6x6
  ! block matrix system are initialised depending on the position
  ! specified by isubstep and nsubsteps.
  !
  ! The structure must have been initialised with smva_initNonlinMatrix!
  type(t_nonlinearSpatialMatrix), intent(inout) :: rnonlinearSpatialMatrix
!</inputoutput>

!</subroutine>

    ! Clear the coefficients
    rnonlinearSpatialMatrix%DidentityY(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Dmass(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Dstokes(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Dygrad(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%DygradAdj(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Dgrady(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%DygradT(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%DgradyAdj(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%DygradTAdj(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Dgrady2(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%DygradT2(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%DgradyAdj2(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%DygradTAdj2(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%DgradyT(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%DBmat(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%DBTmat(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%DidentityP(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%DdirichletBCCY(irow,icolumn)= 0.0_DP
    rnonlinearSpatialMatrix%DdirichletBCCLambda(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%DdirichletBCCXi(irow,icolumn) = 0.0_DP

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine elem_getNodeFunctional (celement,idof,npoints,DrefCoords,Dweights)

!<description>
  ! This routine calculates for an element and a local degree of freedom
  ! the so-called `node-functional`. A `node functional` is a set of points
  ! and corresponding weights which allow to calculate a degree of freedom
  ! by evaluation in special points.
  !
  ! In DrefCoords, the routine returns the coordinates of npoints points on
  ! the reference element, in combination with npoints weights in Dweights.
  ! Multiplying the function values at the coordinates DrefCoords with
  ! the weights in Dweights and summing them up gives the value of the
  ! degree of freedom idof.
  !
  ! Example: For Q1, the coordinates of the idof`th corner are returned,
  ! the weights are 1.0.
  ! For the integral mean value based element Q1tilde, the coordinates of
  ! the 2-pt Gauss Formula on the idof`th edge are returned, the weights
  ! are the cubature weights for the 2-point Gauss Formula on that line.
!</description>

!<input>
  ! Element identifier
  integer(I32), intent(in) :: celement
  
  ! Number of the degree of freedom which `node functional` should be calculated
  integer, intent(in) :: idof
  
!</input>

!<output>
  ! Number of points that represent the `node functional`
  integer, intent(out) :: npoints
  
  ! OPTIONAL: Coordinates of the points representing the `node functional`.
  ! The array must be large rnough.
  ! By calling this routine without DrefCoords and Dweights being specified,
  ! the number of points representing the `node functional` can be obtained.
  real(DP), dimension(:,:), intent(out), optional :: DrefCoords
  
  ! OPTIONAL: Weights representing the node functional.
  ! Multiplying these weights with function values at the points in DrefCoords
  ! gives the value of the degree of freedom idof.
  ! Must be specified if DrefCoords is specified.
  real(DP), dimension(:), intent(out), optional :: Dweights
!</output>

!</subroutine>

    ! Position of cubature points for 2-point Gauss formula on an edge.
    ! Used for Q2T.
    real(DP), parameter :: Q2G1 = -0.577350269189626_DP !-SQRT(1.0_DP/3.0_DP)
    real(DP), parameter :: Q2G2 =  0.577350269189626_DP ! SQRT(1.0_DP/3.0_DP)

    select case (elem_getPrimaryElement(celement))
    
    ! -= 1D Line Elements =-
    case (EL_P0_1D,EL_DG_T0_1D)
    
      ! Point evaluation
      npoints = 1
      if (present(DrefCoords)) then
        DrefCoords(1,1) =  0.0_DP
        Dweights(1) = 1.0_DP
      end if
      
    case (EL_P1_1D,EL_DG_T1_1D)

      ! Point evaluation
      npoints = 2
      if (present(DrefCoords)) then
        select case (idof)
        case (1)
          DrefCoords(1,1) = -1.0_DP
        case (2)
          DrefCoords(1,1) =  1.0_DP
        end select
        Dweights(1) = 1.0_DP
      end if
      
    case (EL_P2_1D,EL_DG_T2_1D)

      ! Point evaluation
      npoints = 3
      if (present(DrefCoords)) then
        select case (idof)
        case (1)
          DrefCoords(1,1) = -1.0_DP
        case (2)
          DrefCoords(1,1) =  1.0_DP
        case (3)
          DrefCoords(1,1) =  0.0_DP
        end select
        Dweights(1) = 1.0_DP
      end if
      
    case (EL_S31_1D)

      ! Point evaluation
      npoints = 4
      if (present(DrefCoords)) then
        select case (idof)
        case (1)
          DrefCoords(1,1) = -1.0_DP
        case (2)
          DrefCoords(1,1) =  1.0_DP
        case (3)
          DrefCoords(1,1) = -1.0_DP
        case (4)
          DrefCoords(1,1) =  1.0_DP
        end select
        Dweights(1) = 1.0_DP
      end if
      
    ! -= 2D Triangle Elements =-
    case (EL_P0)
      ! Point evaluation
      npoints = 1
      if (present(DrefCoords)) then
        DrefCoords(1,1) = 1.0_DP/3.0_DP
        DrefCoords(2,1) = 1.0_DP/3.0_DP
        DrefCoords(3,1) = 1.0_DP/3.0_DP
        Dweights(1) = 1.0_DP
      end if

    case (EL_P1)

      ! Point evaluation
      npoints = 1
      if (present(DrefCoords)) then
        select case (idof)
        case (1)
          DrefCoords(1,1) = 1.0_DP
          DrefCoords(2,1) = 0.0_DP
          DrefCoords(3,1) = 0.0_DP
        case (2)
          DrefCoords(1,1) = 0.0_DP
          DrefCoords(2,1) = 1.0_DP
          DrefCoords(3,1) = 0.0_DP
        case (3)
          DrefCoords(1,1) = 0.0_DP
          DrefCoords(2,1) = 0.0_DP
          DrefCoords(3,1) = 1.0_DP
        end select

        Dweights(1) = 1.0_DP
      end if
      
    case (EL_P2)

      ! Point evaluation
      npoints = 1
      if (present(DrefCoords)) then
        select case (idof)
        case (1)
          DrefCoords(1,1) = 1.0_DP
          DrefCoords(2,1) = 0.0_DP
          DrefCoords(3,1) = 0.0_DP
        case (2)
          DrefCoords(1,2) = 0.0_DP
          DrefCoords(2,1) = 1.0_DP
          DrefCoords(3,1) = 0.0_DP
        case (3)
          DrefCoords(1,1) = 0.0_DP
          DrefCoords(2,1) = 0.0_DP
          DrefCoords(3,1) = 1.0_DP
        case (4)
          DrefCoords(1,1) = 0.5_DP
          DrefCoords(2,1) = 0.5_DP
          DrefCoords(3,1) = 0.0_DP
        case (5)
          DrefCoords(1,1) = 0.0_DP
          DrefCoords(2,1) = 0.5_DP
          DrefCoords(3,1) = 0.5_DP
        case (6)
          DrefCoords(1,1) = 0.5_DP
          DrefCoords(2,1) = 0.0_DP
          DrefCoords(3,1) = 0.5_DP
        end select
        Dweights(1) = 1.0_DP

      end if
      
    case (EL_P3)
      ! Point evaluation
      npoints = 1
      if (present(DrefCoords)) then
        select case (idof)
        case (1)
          DrefCoords(1,1) = 1.0_DP
          DrefCoords(2,1) = 0.0_DP
          DrefCoords(3,1) = 0.0_DP
        case (2)
          DrefCoords(1,1) = 0.0_DP
          DrefCoords(2,1) = 1.0_DP
          DrefCoords(3,1) = 0.0_DP
        case (3)
          DrefCoords(1,1) = 0.0_DP
          DrefCoords(2,1) = 0.0_DP
          DrefCoords(3,1) = 1.0_DP
        case (4)
          DrefCoords(1,1) = 0.5_DP
          DrefCoords(2,1) = 0.5_DP
          DrefCoords(3,1) = 0.0_DP
        case (5)
          DrefCoords(1,1) = 0.0_DP
          DrefCoords(2,1) = 0.5_DP
          DrefCoords(3,1) = 0.5_DP
        case (6)
          DrefCoords(1,1) = 0.5_DP
          DrefCoords(2,1) = 0.0_DP
          DrefCoords(3,1) = 0.5_DP
        case (7)
          DrefCoords(1,1) = 0.5_DP
          DrefCoords(2,1) = 0.5_DP
          DrefCoords(3,1) = 0.0_DP
        case (8)
          DrefCoords(1,1) = 0.0_DP
          DrefCoords(2,1) = 0.5_DP
          DrefCoords(3,1) = 0.5_DP
        case (9)
          DrefCoords(1,1) = 0.5_DP
          DrefCoords(2,1) = 0.0_DP
          DrefCoords(3,1) = 0.5_DP
        end select
        Dweights(1) = 1.0_DP
      end if
      
    case (EL_P1T)
      ! Point evaluation
      npoints = 1
      if (present(DrefCoords)) then
        select case (idof)
        case (1)
          DrefCoords(1,1) = 0.5_DP
          DrefCoords(2,1) = 0.5_DP
          DrefCoords(3,1) = 0.0_DP
        case (2)
          DrefCoords(1,1) = 0.0_DP
          DrefCoords(2,1) = 0.5_DP
          DrefCoords(3,1) = 0.5_DP
        case (3)
          DrefCoords(1,1) = 0.5_DP
          DrefCoords(2,1) = 0.0_DP
          DrefCoords(3,1) = 0.5_DP
        end select
        Dweights(1) = 1.0_DP
      end if

    ! -= 2D Quadrilateral Elements =-
    case (EL_Q0)
      ! Point evaluation
      npoints = 1
      if (present(DrefCoords)) then
        DrefCoords(1,1) = 0.0_DP
        DrefCoords(2,1) = 0.0_DP
        Dweights(1) = 1.0_DP
      end if
      
    case (EL_Q1,EL_DG_Q1_2D)
      ! Point evaluation
      npoints = 1
      if (present(DrefCoords)) then
        select case (idof)
        case (1)
          DrefCoords(1,1) = -1.0_DP
          DrefCoords(2,1) = -1.0_DP
        case (2)
          DrefCoords(1,1) =  1.0_DP
          DrefCoords(2,1) = -1.0_DP
        case (3)
          DrefCoords(1,1) =  1.0_DP
          DrefCoords(2,1) =  1.0_DP
        case (4)
          DrefCoords(1,1) = -1.0_DP
          DrefCoords(2,1) =  1.0_DP
        end select
        Dweights(1) = 1.0_DP
      end if

    case (EL_Q2,EL_DG_Q2_2D)
      ! Point evaluation
      npoints = 1
      if (present(DrefCoords)) then
        select case (idof)
        case (1)
          DrefCoords(1,1) = -1.0_DP
          DrefCoords(2,1) = -1.0_DP
        case (2)
          DrefCoords(1,1) =  1.0_DP
          DrefCoords(2,1) = -1.0_DP
        case (3)
          DrefCoords(1,1) =  1.0_DP
          DrefCoords(2,1) =  1.0_DP
        case (4)
          DrefCoords(1,1) = -1.0_DP
          DrefCoords(2,1) =  1.0_DP
        case (5)
          DrefCoords(1,1) =  0.0_DP
          DrefCoords(2,1) = -1.0_DP
        case (6)
          DrefCoords(1,1) =  1.0_DP
          DrefCoords(2,1) =  0.0_DP
        case (7)
          DrefCoords(1,1) =  0.0_DP
          DrefCoords(2,1) =  1.0_DP
        case (8)
          DrefCoords(1,1) = -1.0_DP
          DrefCoords(2,1) =  0.0_DP
        case (9)
          DrefCoords(1,1) =  0.0_DP
          DrefCoords(2,1) =  0.0_DP
        end select
        
        Dweights(1) = 1.0_DP
      end if
        
    case (EL_QP1)
      DrefCoords      = 0.0_DP
      
      ! Not yet implemented
      call sys_halt()      
      
    case (EL_Q1T)
      if (iand(celement,int(2**16,I32)) .eq. 0) then
        ! Conformal element.
        !
        ! Point evaluation
        npoints = 1
        if (present(DrefCoords)) then
          select case (idof)
          case (1)
            DrefCoords(1,1) =  0.0_DP
            DrefCoords(2,1) = -1.0_DP
          case (2)
            DrefCoords(1,1) =  1.0_DP
            DrefCoords(2,1) =  0.0_DP
          case (3)
            DrefCoords(1,1) =  0.0_DP
            DrefCoords(2,1) =  1.0_DP
          case (4)
            DrefCoords(1,1) = -1.0_DP
            DrefCoords(2,1) =  0.0_DP
          end select
          
          Dweights(1) = 1.0_DP
        end if

      else 
      
        ! Nonconformal element.
        !
        ! Use 2-point Gauss formula which calculated the integral mean value.
        npoints = 2
        if (present(DrefCoords)) then
          select case (idof)
          case (1)
            DrefCoords(1,1) =  Q2G1
            DrefCoords(2,1) = -1.0_DP
            DrefCoords(1,2) =  Q2G2
            DrefCoords(2,2) = -1.0_DP
          case (2)
            DrefCoords(1,1) =  1.0_DP
            DrefCoords(2,1) =  Q2G1
            DrefCoords(1,2) =  1.0_DP
            DrefCoords(2,2) =  Q2G2
          case (3)
            DrefCoords(1,1) =  Q2G2
            DrefCoords(2,1) =  1.0_DP
            DrefCoords(1,2) =  Q2G1
            DrefCoords(2,2) =  1.0_DP
          case (4)
            DrefCoords(1,1) = -1.0_DP
            DrefCoords(2,1) =  Q2G2
            DrefCoords(1,2) = -1.0_DP
            DrefCoords(2,2) =  Q2G1
          end select
          
          Dweights(1) = 0.5_DP
          Dweights(2) = 0.5_DP
        end if
      end if
      
    case (EL_Q1TB)
      if (iand(celement,int(2**16,I32)) .eq. 0) then
        ! Conformal element.
        !
        ! Point evaluation
        npoints = 1
        if (present(DrefCoords)) then
          select case (idof)
          case (1)
            DrefCoords(1,1) =  0.0_DP
            DrefCoords(2,1) = -1.0_DP
          case (2)
            DrefCoords(1,1) =  1.0_DP
            DrefCoords(2,1) =  0.0_DP
          case (3)
            DrefCoords(1,1) =  0.0_DP
            DrefCoords(2,1) =  1.0_DP
          case (4)
            DrefCoords(1,1) = -1.0_DP
            DrefCoords(2,1) =  0.0_DP
          case (5)
            DrefCoords(1,1) =  0.0_DP
            DrefCoords(2,1) =  0.0_DP
          end select
          
          Dweights(1) = 1.0_DP
        end if

      else 
      
        ! Nonconformal element.
        !
        ! Use 2-point Gauss formula which calculated the integral mean value.
        npoints = 2
        if (present(DrefCoords)) then
          select case (idof)
          case (1)
            DrefCoords(1,1) =  Q2G1
            DrefCoords(2,1) = -1.0_DP
            DrefCoords(1,2) =  Q2G2
            DrefCoords(2,2) = -1.0_DP

            Dweights(1) = 0.5_DP
            Dweights(2) = 0.5_DP
          case (2)
            DrefCoords(1,1) =  1.0_DP
            DrefCoords(2,1) =  Q2G1
            DrefCoords(1,2) =  1.0_DP
            DrefCoords(2,2) =  Q2G2

            Dweights(1) = 0.5_DP
            Dweights(2) = 0.5_DP
          case (3)
            DrefCoords(1,1) =  Q2G2
            DrefCoords(2,1) =  1.0_DP
            DrefCoords(1,2) =  Q2G1
            DrefCoords(2,2) =  1.0_DP

            Dweights(1) = 0.5_DP
            Dweights(2) = 0.5_DP
          case (4)
            DrefCoords(1,1) = -1.0_DP
            DrefCoords(2,1) =  Q2G2
            DrefCoords(1,2) = -1.0_DP
            DrefCoords(2,2) =  Q2G1

            Dweights(1) = 0.5_DP
            Dweights(2) = 0.5_DP
          case (5)
            ! Centre: Point value
            npoints = 1
            DrefCoords(1,1) =  0.0_DP
            DrefCoords(2,1) =  0.0_DP

            Dweights(1) = 1.0_DP
          end select
          
        end if
      end if
      
    ! -= 3D Tetrahedron Elements =-
    case (EL_P0_3D)
      ! Point evaluation
      npoints = 1
      if (present(DrefCoords)) then
        DrefCoords(1,1) = 1.0_DP/3.0_DP
        DrefCoords(2,1) = 1.0_DP/3.0_DP
        DrefCoords(3,1) = 1.0_DP/3.0_DP
        
        DrefCoords(4,1) = 1.0_DP/3.0_DP
        Dweights(1) = 1.0_DP
      end if
      
    case (EL_P1_3D)
      ! Point evaluation
      npoints = 1
      if (present(DrefCoords)) then
        select case (idof)
        case (1)
          DrefCoords(1,1) = 1.0_DP
          DrefCoords(2,1) = 0.0_DP
          DrefCoords(3,1) = 0.0_DP
          DrefCoords(4,1) = 0.0_DP
        case (2)
          DrefCoords(1,1) = 0.0_DP
          DrefCoords(2,1) = 1.0_DP
          DrefCoords(3,1) = 0.0_DP
          DrefCoords(4,1) = 0.0_DP
        case (3)
          DrefCoords(1,1) = 0.0_DP
          DrefCoords(2,1) = 0.0_DP
          DrefCoords(3,1) = 1.0_DP
          DrefCoords(4,1) = 0.0_DP
        case (4)
          DrefCoords(1,1) = 0.0_DP
          DrefCoords(2,1) = 0.0_DP
          DrefCoords(3,1) = 0.0_DP
          DrefCoords(4,1) = 1.0_DP
        end select
        
        Dweights(1) = 1.0_DP
      end if
          
    ! -= 3D Hexahedron Elements =-
    case (EL_Q0_3D)
      ! Point evaluation
      npoints = 1
      if (present(DrefCoords)) then
        DrefCoords(1,1) = 0.0_DP
        DrefCoords(2,1) = 0.0_DP
        DrefCoords(3,1) = 0.0_DP
        Dweights(1) = 1.0_DP
      end if

    case (EL_Q1_3D)
      ! Point evaluation
      npoints = 1
      if (present(DrefCoords)) then
        select case (idof)
        case (1)
          DrefCoords(1,1) = -1.0_DP
          DrefCoords(2,1) = -1.0_DP
          DrefCoords(3,1) = -1.0_DP
        case (2)
          DrefCoords(1,2) =  1.0_DP
          DrefCoords(2,2) = -1.0_DP
          DrefCoords(3,2) = -1.0_DP
        case (3)
          DrefCoords(1,3) =  1.0_DP
          DrefCoords(2,3) =  1.0_DP
          DrefCoords(3,3) = -1.0_DP
        case (4)
          DrefCoords(1,4) = -1.0_DP
          DrefCoords(2,4) =  1.0_DP
          DrefCoords(3,4) = -1.0_DP
        case (5)
          DrefCoords(1,5) = -1.0_DP
          DrefCoords(2,5) = -1.0_DP
          DrefCoords(3,5) =  1.0_DP
        case (6)
          DrefCoords(1,6) =  1.0_DP
          DrefCoords(2,6) = -1.0_DP
          DrefCoords(3,6) =  1.0_DP
        case (7)
          DrefCoords(1,7) =  1.0_DP
          DrefCoords(2,7) =  1.0_DP
          DrefCoords(3,7) =  1.0_DP
        case (8)
          DrefCoords(1,8) = -1.0_DP
          DrefCoords(2,8) =  1.0_DP
          DrefCoords(3,8) =  1.0_DP
        end select
        
        Dweights(1) = 1.0_DP
      end if

    ! -= 3D Pyramid Elements =-
    
    ! -= 3D Prism Elements =-
      
    case default
      call output_line ("Unsupported element.", &
                        OU_CLASS_ERROR,OU_MODE_STD,"elem_getNodeFunctional")
        call sys_halt()
    end select

  end subroutine


  ! ***************************************************************************
  
!<subroutine>

  subroutine linf_getBoundaryOperatorMatrix (rlinform,rmatrix,bclear,rboundaryRegion,&
      fcoeff_buildVectorScBdr2D_sim,rcollection)

!<description>
  ! For a given set of degrees of freedom, this routine calculates the matrix
  ! of an operator rlinform if being applied to a finite element vector:
  ! Multiplying the matrix rmatrix with a finite element vector x realises
  ! the application of the operator rlinform to the underlying finite element
  ! function in strong form in the boundary refion rboundaryRegion.
!</description>

!<input>
  ! The operator to be realised in matrix form
  type(t_linearForm), intent(in) :: rlinform

  ! Whether or not to clear the matrix in advance.
  logical, intent(in) :: bclear
  
  ! Boundary region where the operator is applied.
  type(t_boundaryRegion), intent(in) :: rboundaryRegion
  
  ! OPTIONAL: Coefficient function
  include '../../../kernel/DOFMaintenance/intf_coefficientVectorScBdr2D.inc'
  optional :: fcoeff_buildVectorScBdr2D_sim
  
  ! OPTIONAL: Collection passed to the coefficient function
  type(t_collection), intent(inout), optional :: rcollection
!</input>

!<inputoutput>
  ! Matrix to be modified.
  type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

    ! lcoal variables
    integer, dimension(:), allocatable :: IdofsAtBoundary
    integer :: ndofs, idofidx
    
    real(DP), dimension(:), pointer :: p_Ddata
    integer, dimension(:), pointer :: p_Kld, p_Kcol
    
    integer(I32) :: celementTrial,celementTest
    integer :: iderType
    real(DP), dimension(2,2,1) :: DrefCoords,DrealCoords
    real(DP), dimension(2) :: Dweights
    real(DP), dimension(2,2,1) :: Dcoefficients
    integer, dimension(32,1), target :: IdofGlobTrial,IdofGlobTest
    real(DP), dimension(:,:), allocatable :: DpointPar
    integer :: ndofLocTrial,ndofLocTest,idoflocTrial,idofLocTest,npoints,nve,ipoint
    
    integer :: nel,ielidx,iel,i,i1,ia
    integer, dimension(:), allocatable, target :: Ielements, Iedges

    ! Transformation
    integer(I32) :: ctrafoType
    real(DP), dimension(TRAFO_MAXDIMREFCOORD) :: DparPoint
    real(DP) :: dpar1, dpar2
    
    ! Triangulation
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    integer :: nvt, iedge
    
    ! Values of basis functions and DOF`s
    real(DP), dimension(EL_MAXNBAS,EL_MAXNDER) :: Dbas
    logical, dimension(EL_MAXNDER) :: Bder

    ! Evaluation structure and tag
    type(t_evalElement) :: revalElement
    integer(I32) :: cevaluationTag
    type(t_domainIntSubset) :: rintSubset

      ! Probably clear the matrix
      if (bclear) then
        call lsyssc_clearMatrix (rmatrix)
      end if

      ! Get information about the matrix.
      
      ! Matrix format 9.
      call lsyssc_getbase_double (rmatrix,p_Ddata)
      call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      call lsyssc_getbase_Kld (rmatrix,p_Kld)
      
      ! Information from the triangulation
      call storage_getbase_int2d (&
          rmatrix%p_rspatialDiscrTrial%p_rtriangulation%h_IedgesAtElement,&
          p_IedgesAtElement)
      call storage_getbase_double (&
          rmatrix%p_rspatialDiscrTrial%p_rtriangulation%h_DvertexParameterValue,&
          p_DvertexParameterValue)
      nvt = rmatrix%p_rspatialDiscrTrial%p_rtriangulation%NVT
      
      ! Which derivatives of basis functions are needed?
      ! Check the descriptors of the linear form and set BDERxxxx
      ! according to these.
      Bder(:) = .false.
      
      ! Loop through the additive terms
      do i = 1,rlinform%itermCount
        ! The desriptor Idescriptors gives directly the derivative
        ! which is to be computed! Build templates for BDER.
        ! We do not compute the actual BDER here, as there might be some special
        ! processing if trial/test functions are identical!
        !
        ! At first build the descriptors for the trial functions
        i1 = rlinform%Idescriptors(i)
        
        if ((i1 .le.0) .or. (i1 .gt. DER_MAXNDER)) then
          call output_line ("Invalid descriptor!",&
              OU_CLASS_ERROR,OU_MODE_STD,"linf_getBoundaryOperatorMatrix")
          call sys_halt()
        endif
        
        Bder(i1)=.true.
      end do
      
      ! Allocate memory if not done already.
      ! Calculate the DOF`s which are affected.
      call bcasm_getDOFsInBDRegion (rmatrix%p_rspatialDiscrTest,&
          rboundaryRegion,ndofs=ndofs)

      allocate (IdofsAtBoundary(ndofs))
      call bcasm_getDOFsInBDRegion (rmatrix%p_rspatialDiscrTest,&
          rboundaryRegion,ndofs=ndofs,IdofsArray=IdofsAtBoundary)

      ! Get a list of cells in the boundary region
      call bcasm_getElementsInBdRegion (&
          rmatrix%p_rspatialDiscrTrial%p_rtriangulation,rboundaryRegion,nel)
      allocate (Ielements(nel))
      allocate (Iedges(nel))
      call bcasm_getElementsInBdRegion (&
          rmatrix%p_rspatialDiscrTrial%p_rtriangulation,rboundaryRegion,nel,&
          Ielements,IedgeLocal=Iedges)
      
      do i = 1,nel
        ! Calculate the actual edge numbers
        iedge = p_IedgesAtElement(Iedges(i),Ielements(i))
        
        ! Find the position of the edge on the boundary.
        call tria_searchBoundaryEdge(iedge, &
            rmatrix%p_rspatialDiscrTrial%p_rtriangulation, Iedges(i))
      end do
      
      ! Get information about the FEM space
      celementTrial = rmatrix%p_rspatialDiscrTrial%RelementDistr(1)%celement
      celementTest = rmatrix%p_rspatialDiscrTest%RelementDistr(1)%celement
      ndofLocTrial = elem_igetNDofLoc(celementTrial)
      ndofLocTest = elem_igetNDofLoc(celementTest)
      nve = elem_igetNVE(celementTrial)
      ctrafoType = elem_igetTrafoType(celementTrial)
      cevaluationTag = elem_getEvaluationTag(celementTrial)
     
      ! The weights are =1 by default, only changed if
      ! the coefficients are nonconstant.
      Dcoefficients(:,:,:) = 1.0_DP
     
      ! Simultaneously loop over the DOF`s and the cells on the boundary
      idofidx = 1
      ielidx = 1
      do
        ! Stop if there is no call or DOF available anymore
        if ((idofidx .gt. ndofs) .or. (ielidx .gt. nel)) exit
        
        ! Get the DOF`s of the test element. These are the columns in
        ! the matrix to be modified.
        call dof_locGlobMapping(rmatrix%p_rspatialDiscrTrial, &
            Ielements(ielidx), IdofGlobTrial(:,1))
        
        ! Get the DOF`s of the test element. These are the rows in
        ! the matrix to be modified.
        call dof_locGlobMapping(rmatrix%p_rspatialDiscrTest, &
            Ielements(ielidx), IdofGlobTest(:,1))
        
        ! Figure out the local DOF of the current global DOF
        ! on the boundary.
        do idoflocTest = 1,ndoflocTest
          if (IdofGlobTest(idoflocTest,1) .eq. IdofsAtBoundary(idofidx)) exit
        end do
        
        ! Not found -> next element, repeat
        if (idoflocTest .gt. ndoflocTest) then
          ielidx = ielidx + 1
          cycle
        end if
        
        ! Otherwise, the local DOF is found.
        !
        ! Calculate the node functional of that point.
        ! The test space specifies the way to evaluate.
        call elem_getNodeFunctional (celementTest,idoflocTest,npoints,DrefCoords(:,:,1),Dweights)
        
        ! If a callback function is given, calculate the weights in the points.
        if (present(fcoeff_buildVectorScBdr2D_sim)) then
          
          ! Allocate memory for the parameter values.
          if (.not. allocated (DpointPar)) then
            allocate (DpointPar(npoints,1))
          end if
          
          if (ubound(DpointPar,1) .ne. npoints) then
            deallocate(DpointPar)
            allocate (DpointPar(npoints,1))
          end if
          
          ! Approximate the parameter values of the points.
          ! Take the start parameter value of the edge of the current element
          ! and add the position of the point to get an approximate parameter value.
          
          ! Start/end parameter value of the edge
          dpar1 = p_DvertexParameterValue(Iedges(ielidx))
          if (Iedges(ielidx) .lt. size (p_DvertexParameterValue)) then
            dpar2 = p_DvertexParameterValue(Iedges(ielidx)+1)
          else
            dpar2 = boundary_dgetMaxParVal(&
                rmatrix%p_rspatialDiscrTrial%p_rboundary,rboundaryRegion%iboundCompIdx)
          end if
          
          ! Quad element
          do ipoint=1,npoints
          
            ! Determine the edge
            if (DrefCoords(2,ipoint,1) .eq. -1.0_DP) then
              ! Bottom edge
              call mprim_linearRescale(DrefCoords(1,ipoint,1),&
                  -1.0_DP,1.0_DP,dpar1,dpar2,DpointPar(ipoint,1))
            else if (DrefCoords(1,ipoint,1) .eq. 1.0_DP) then
              ! Right edge
              call mprim_linearRescale(DrefCoords(2,ipoint,1),&
                  -1.0_DP,1.0_DP,dpar1,dpar2,DpointPar(ipoint,1))
            else if (DrefCoords(2,ipoint,1) .eq. 1.0_DP) then
              ! Top edge
              call mprim_linearRescale(DrefCoords(1,ipoint,1),&
                  1.0_DP,-1.0_DP,dpar1,dpar2,DpointPar(ipoint,1))
            else  if (DrefCoords(1,ipoint,1) .eq. -1.0_DP) then
              ! Left edge
              call mprim_linearRescale(DrefCoords(2,ipoint,1),&
                  1.0_DP,-1.0_DP,dpar1,dpar2,DpointPar(ipoint,1))
            else
              call output_line ("Evaluation point not on the boundary!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"linf_getBoundaryOperatorMatrix")
              call sys_halt()
            end if
            
          end do
        
          ! Convert to length parametrization
          call boundary_convertParameterList(&
              rmatrix%p_rspatialDiscrTrial%p_rboundary, rboundaryRegion%iboundCompIdx, &
              DpointPar,DpointPar,BDR_PAR_01,BDR_PAR_LENGTH)
                  
          ! Prepare the call of the callback routine
          rintSubset%ielementStartIdx      =  ielidx
          rintSubset%p_Ielements           => Ielements
          rintSubset%p_IdofsTrial          => IdofGlobTrial
          rintSubset%celement = celementTrial
          call fcoeff_buildVectorScBdr2D_sim (rmatrix%p_rspatialDiscrTrial, rlinform, &
                  1, npoints, DrefCoords, rboundaryRegion%iboundCompIdx, DpointPar,&
                  IdofGlobTest(:,:), rintSubset, Dcoefficients, rcollection)
                  
        end if
        
        ! Multiplying a FE function in the npoints points DrefCoords, weighted
        ! by Dweights gives the value of the DOF. This summation has to be emulated by
        ! the matrix.
        
        ! We need the values of the basis functions in the above coordinates 
        ! DrefCoords. Multiplied by Dweights, these values define the content
        ! of the matrix.
          
        ! Loop over the points
        do ipoint = 1,npoints
        
          ! Get the element shape information in that point
          call elprep_prepareForEvaluation (revalElement, &
              cevaluationTag, &
              rmatrix%p_rspatialDiscrTrial%p_rtriangulation, Ielements(ielidx), &
              ctrafoType, DrefCoords(:,ipoint,1))
          
          ! Call the element to calculate the values of the basis functions
          ! in the point.
          call elem_generic2 (celementTrial, revalElement, Bder, Dbas)
          
          ! Loop over the terms.
          do ia = 1,rlinform%itermCount
          
            ! Current derivative to be applied to the basis function
            iderType = rlinform%Idescriptors(ia)
          
            ! Sum up Dbas and Dweights to the matrix entries.
            ! The row is given by IdofsAtBoundary(idofidx) which is the global DOF to calculate.
            ! The columns are given by IdofGlob and must be found in Kcol.
            idofloop: do idoflocTrial = 1,ndoflocTrial
            
              ! Find the position in the matrix to be modified.
              do i = p_Kld(IdofsAtBoundary(idofidx)), p_Kld(IdofsAtBoundary(idofidx)+1)-1
              
                if (p_Kcol(i) .eq. IdofGlobTrial(idoflocTrial,1)) then
                  ! Position in the matrix to be modified is i.
                  ! Sum up the weights.
                  p_Ddata(i) = p_Ddata(i) + Dcoefficients(ia,ipoint,1) * rlinform%Dcoefficients(ia) * &
                      Dweights(ipoint) * Dbas(idoflocTrial,iderType)
                  
                  ! Next DOF
                  cycle idofloop
                end if
              
              end do
              
            end do idofloop 
          
          end do
          
        end do
        
        ! All points processed, matrix entries calculating the DOF IdofsAtBoundary(idofidx)
        ! calculated. Proceed with the next DOF (=next row in the matrix)
        idofidx = idofidx + 1
        
      end do
      
      ! Release memory
      deallocate (Ielements)
      deallocate (IdofsAtBoundary)

      if (allocated (DpointPar)) then
        deallocate (DpointPar)
      end if

  end subroutine

! *****************************************************************

!<subroutine>

  subroutine smva_prepareViscoAssembly (rphysics,rcollection,rvelocityvector)
  
  use basicgeometry
  use triangulation
  use scalarpde
  use domainintegration
  use spatialdiscretisation
  use collection
  
!<description>
  ! Based on the the input parameters, this routine prepares a collection
  ! structure rcollection for being used with ffunctionViscoModel.
  ! The structure realises the following setting:
  !
  ! IquickAccess(1) = cviscoModel
  ! IquickAccess(2) = type of the tensor
  ! DquickAccess(1) = nu
  ! DquickAccess(2) = dviscoexponent
  ! DquickAccess(3) = dviscoEps
  ! DquickAccess(4) = dviscoYield
  ! p_rvectorQuickAccess1 => evaluation velocity vector
  !
!</description>
  
!<input>
  ! Physics of the problem
  type(t_settings_physics), intent(in), target :: rphysics
  
  ! OPTIONAL: The current velocity/pressure vector. May not be present.
  type(t_vectorBlock), intent(in), target, optional :: rvelocityvector
!</input>

!<inputoutput>
  ! Collection structure to be prepared.
  type(t_collection), intent(inout) :: rcollection
!</inputoutput>
  
!</subroutine>

    ! IquickAccess(1) = cviscoModel
    rcollection%IquickAccess(1) = rphysics%cviscoModel
    
    ! IquickAccess(2) = type of the tensor
    rcollection%IquickAccess(2) = rphysics%isubequation
    
    ! DquickAccess(1) = nu
    rcollection%DquickAccess(1) = rphysics%dnuConst
    
    ! DquickAccess(2) = dviscoexponent
    ! DquickAccess(3) = dviscoEps
    ! DquickAccess(4) = dviscoYield
    rcollection%DquickAccess(2) = rphysics%dviscoexponent
    rcollection%DquickAccess(3) = rphysics%dviscoEps
    rcollection%DquickAccess(4) = rphysics%dviscoYield
    
    ! The first quick access array specifies the evaluation point
    ! of the velocity -- if it exists.
    nullify(rcollection%p_rvectorQuickAccess1)
    if (present(rvelocityvector)) &
        rcollection%p_rvectorQuickAccess1 => rvelocityvector

  end subroutine
    
! *****************************************************************

!<subroutine>

  subroutine ffunctionViscoModel (cterm,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset, &
                Dcoefficients,rcollection)
  
  use basicgeometry
  use triangulation
  use scalarpde
  use domainintegration
  use spatialdiscretisation
  use collection
  
!<description>
  ! This subroutine is called during the calculation of the SD operator. It has to
  ! compute the coefficients in front of the terms.
  !
  ! The routine accepts a set of elements and a set of points on these
  ! elements (cubature points) in in real coordinates.
  ! According to the terms in the linear form, the routine has to compute
  ! simultaneously for all these points.
  !
  ! The following data must be passed to this routine in the collection in order
  ! to work correctly:
  !
  ! IquickAccess(1) = cviscoModel
  ! IquickAccess(2) = type of the tensor
  ! DquickAccess(1) = nu
  ! DquickAccess(2) = dviscoexponent
  ! DquickAccess(3) = dviscoEps
  ! p_rvectorQuickAccess1 => evaluation velocity vector
  ! p_rnextCollection => user defined collection structure
  !
!</description>
  
!<input>
  ! Term which is to be computed.
  ! =0: Calculate the $\nu$ values in front of the Laplace.
  ! =1: Calculate the $\alpha$ values in front of the Mass matrix.
  integer, intent(in) :: cterm

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in)  :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional      :: rcollection
  
!</input>

!<output>
    ! This array has to receive the values of the coefficients
    ! in all the points specified in Dpoints.
    ! cterm specifies what to evaluate.
    real(DP), dimension(:,:), intent(out) :: Dcoefficients
!</output>
  
!</subroutine>

    ! local variables
    integer :: cviscoModel,i,j,isubEquation
    real(DP) :: dnu,dviscoexponent,dviscoEps,dviscoYield
    type(t_vectorBlock), pointer :: p_rvector
    integer, dimension(2) :: Ibounds
    real(DP), dimension(:,:,:), allocatable :: Ddata
    
    if (.not. present(rcollection)) then
      Dcoefficients(:,:) = 0.0_DP
      return
    end if
    
    ! Get the parameters from the collection,
    ! as specified in the call to the SD assembly routine.
    cviscoModel = rcollection%IquickAccess(1)
    isubEquation = rcollection%IquickAccess(2)
    dnu = rcollection%DquickAccess(1)
    dviscoexponent = rcollection%DquickAccess(2)
    dviscoEps = rcollection%DquickAccess(3)
    dviscoYield = rcollection%DquickAccess(4)
    p_rvector => rcollection%p_rvectorQuickAccess1
    
    ! p_rvector may point to NULL if there is no nonlinearity
    
    ! If we have a nonlinear viscosity, calculate
    !    z := D(u):D(u) = ||D(u)||^2
    ! The isubEquation defines the shape of the tensor, which may me
    !    D(u) = grad(u)
    ! or D(u) = 1/2 ( grad(u) + grad(u)^T )
    select case (cviscoModel)
    case (0)
      ! Constant viscosity. This is actually not used as the routine is
      ! not called if nu is constant.
      Dcoefficients(:,:) = dnu
      
    case (1,2,3)
    
      ! Allocate memory do calculate D(u) in all points
      Ibounds = ubound(Dcoefficients)
      allocate(Ddata(Ibounds(1),Ibounds(2),5))
      
      ! Evaluate D(u).
      call fevl_evaluate_sim (p_rvector%RvectorBlock(1), &
                                rdomainIntSubset, DER_DERIV_X, Ddata(:,:,2))
      call fevl_evaluate_sim (p_rvector%RvectorBlock(1), &
                                rdomainIntSubset, DER_DERIV_Y, Ddata(:,:,3))
      call fevl_evaluate_sim (p_rvector%RvectorBlock(2), &
                                rdomainIntSubset, DER_DERIV_X, Ddata(:,:,4))
      call fevl_evaluate_sim (p_rvector%RvectorBlock(2), &
                                rdomainIntSubset, DER_DERIV_Y, Ddata(:,:,5))
                         
      ! Calculate ||D(u)||^2 to Ddata(:,:,1):
      select case (isubequation)
      case (0)
        ! D(u) = grad(u)
        do i=1,Ibounds(2)
          do j=1,Ibounds(1)
            Ddata(j,i,1) = Ddata(j,i,2)**2 + Ddata(j,i,3)**2 + &
                           Ddata(j,i,4)**2 + Ddata(j,i,5)**2
          end do
        end do
        
      case (1)
        ! D(u) = 1/2 ( grad(u) + grad(u)^T )
        do i=1,Ibounds(2)
          do j=1,Ibounds(1)
            Ddata(j,i,1) = (Ddata(j,i,2)**2 + &
                            0.5_DP * (Ddata(j,i,3) + Ddata(j,i,4))**2 + &
                            Ddata(j,i,5)**2)
          end do
        end do
        
      case default
      
        Ddata(:,:,:) = 0.0_DP
        
      end select
      
      ! Calculate the viscosity.
      ! (WARNING: This inner select-case builds upon the outer select case above.
      !  Be careful that all cviscoModel cases here also appear above!!!)
      select case (cviscoModel)
      case (1)
        ! Power law:
        !   nu = nu_0 * z^(dviscoexponent/2 - 1),
        !   nu_0 = 1/RE,
        !   z = ||D(u)||^2+dviscoEps
        Dcoefficients(:,:) = dnu * (Ddata(:,:,1)+dviscoEps)**(0.5_DP*dviscoexponent-1.0_DP)

      case (2)
        ! Bingham fluid:
        !   nu = nu_0 + sqrt(2)/2 * dviscoyield / sqrt(|D(u)||^2+dviscoEps^2),
        !   nu_0 = 1/RE
        Dcoefficients(:,:) = dnu + 0.5_DP * sqrt(2.0_DP) * &
            dviscoyield / sqrt(Ddata(:,:,1) + dviscoEps**2)
        
      case (3)
        ! General viscoplastic fluid:
        !   nu = nu_0 +  + sqrt(2)/2 * dviscoyield * z^(dviscoexponent/2 - 1),
        !   nu_0 = 1/RE,
        !   z = ||D(u)||^2 + dviscoEps^2
        Dcoefficients(:,:) = dnu + 0.5_DP * sqrt(2.0_DP) * &
            dviscoyield * ( Ddata(:,:,1) + dviscoEps**2 )**( 0.5_DP * dviscoexponent - 1.0_DP)

      end select
      
      ! Deallocate needed memory.
      deallocate(Ddata)

    case default
    
      ! Viscosity specified by the callback function getNonconstantViscosity.
      ! Call it and pass the user-defined collection.
      call user_getNonconstantViscosity (cterm,rdiscretisation, &
          nelements,npointsPerElement,Dpoints, &
          IdofsTest,rdomainIntSubset, &
          Dcoefficients,p_rvector,rcollection%p_rnextCollection)
    
    end select

  end subroutine
  
end module
