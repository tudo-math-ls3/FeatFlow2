!##############################################################################
!# ****************************************************************************
!# <name> cc2dmatrixassembly </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the very basic matrix assembly routines for the
!# core equation. It is independent of any nonlinear iteration and provides
!# just one functionality: Assemble a matrix or a vector based on a given
!# set of parameters.
!#
!# The discretised core equation reads at the moment:
!#
!#  $$        A_1 y   +  \eta B p             = f_1 $$
!#  $$   \tau B^T y   +  \mu P p  + \kappa K  = f_2 $$
!#  $$                               \jota C  = f_3 $$
!#
!# with
!#
!#   $$ A_1 = \alpha M  +  \theta L  +  \gamma N(y) + \text{newton }N*(y)$$
!#
!# and
!#
!#   $M$     = mass matrix,
!#   $L$     = Stokes matrix ($\nu$*Laplace),
!#   $N(y)$  = Nonlinearity $y\delta(\cdot)$ includung stabilisation,
!#   $N*(y)$ = Adjoint term $\cdot\delta(y)$ of the nonlinearity,
!#             used for the Newton matrix
!#
!#   $\alpha$ - weight in front of the mass matrix;
!#                =0 for stationary problem,
!#   $\theta$ - weight for the Laplace matrix,
!#   $\gamma$ - weight in front of the nonlinearity;
!#                =0 for Stokes system,
!#   $\eta$   - Switches the 'B'-term on/off,
!#   $\tau$   - Switches the 'B^T'-term on/off,
!#   $\mu$    - Switches the 'P'-term on/off,
!#   $\kappa$ - Switches the 'K'-term on/off,
!#   $\jota$  - Switches the 'C'-term on/off,
!#   newton   - Weight for the Newton term
!#
!# This equation can be written as a nonlinear system $A(y)(y,p) = (f1,f2)$
!# with a nonlinear matrix $A(\cdot)$. The structure t_nonlinearCCmatrix
!# contains a description of this matrix, With this description, it is possible
!# to do matrix vector multiplication or to 'evaluate' the matrix at a
!# given 'point' $y$ to get the 'linearised' matrix $A(y)$.
!#
!# The module contains the following routines:
!#
!# 1.) cc_assembleMatrix
!#     -> Assembles a matrix based on a set of input parameters, i.e.
!#        evaluates the system matrix $A(\cdot)$ at a point $y$ to create
!#        the linear matrix $A(y)$.
!#
!# 2.) cc_nonlinearMatMul
!#     -> Performs a matrix vector multiplication $d:=A(y)x+d$.
!#
!# 3.) ccmva_prepareViscoAssembly
!#     -> Prepare a collection for the use in ffunctionViscoModel
!#
!# 4.) ffunctionViscoModel
!#     -> Auxiliary function that defines the nonconstant viscosity
!# </purpose>
!##############################################################################

module ccmatvecassembly

  use fsystem
  use storage
  use genoutput
  use basicgeometry
  use boundary
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use derivatives
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use matrixrestriction
  use linearsystemscalar
  use bilinearformevaluation
  use linearformevaluation
  use trilinearformevaluation
  use matrixio
  use statistics
  use linearsystemblock
  use linearsolver
  use linearsolverautoinitialise
  use collection
  use feevaluation
  use convection
  use stdoperators
  
  use ccbasic
  use cccallback
  
  use pprocnavierstokes
  
  implicit none
  
!<constants>

!<constantblock description="Identifiers for the 'coperation' input parameter of the matrix assembly routine">

  ! Allocate memory if necessary.
  integer(I32), parameter :: CCMASM_ALLOCMEM              = 1
  
  ! Compute all matrix entries.
  integer(I32), parameter :: CCMASM_COMPUTE               = 2
  
  ! Allocate memory and compute matrix entries.
  integer(I32), parameter :: CCMASM_ALLOCANDCOMPUTE       = 3
  
  ! Bypass memory allocation for matrices.
  integer(I32), parameter :: CMASM_QUICKREFERENCES        = 4
  
!</constantblock>

!<constantblock description="Identifiers for the IUPWIND parameter that specifies how to set up the nonlinearity or stabilisation.">

  ! Streamline diffusion; configured by dupsam
  integer, parameter :: CCMASM_STAB_STREAMLINEDIFF    = 0

  ! 1st-order upwind; configured by dupsam
  integer, parameter :: CCMASM_STAB_UPWIND            = 1
  
  ! Edge-oriented stabilisation; configured by dupsam as 'gamma'
  integer, parameter :: CCMASM_STAB_EDGEORIENTED      = 2

  ! Fast Edge-oriented stabilisation; configured by dupsam as 'gamma'. Preconputed matrix.
  integer, parameter :: CCMASM_STAB_FASTEDGEORIENTED  = 3

  ! Streamline diffusion; configured by dupsam. New implementation
  integer, parameter :: CCMASM_STAB_STREAMLINEDIFF2   = 4

  ! Edge-oriented stabilisation; configured by dupsam as 'gamma'.
  ! Nonlinear terms set up with the element independend SD routine.
  integer, parameter :: CCMASM_STAB_EDGEORIENTED2     = 5
!</constantblock>

!<constantblock description="Matrix type ID's specifying the general matrix class to set up.">

  ! Standard matrix.
  integer, parameter :: CCMASM_MTP_AUTOMATIC         = 0
  
  ! Standard matrix with decoupled velocity blocks
  integer, parameter :: CCMASM_MTP_DECOUPLED         = 1
  
  ! Extended 'full-tensor' matrix with submatrices A11, A12, A21, A22, all independent from
  ! each other.
  integer, parameter :: CCMASM_MTP_FULLTENSOR        = 2

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! This routine describes the nonlinear system matrix. The system matrix
  ! does actually not exist in memory -- since it is nonlinear! Therefore,
  ! this structure contains all parameters and settings which are necessary
  ! do apply(!) the matrix to a vector or to evaluate it.
  ! ('Evaluate a nonlinear matrix' means: Using a given FE-function $y$,
  ! assemble the linear matrix A(y).)
  type t_nonlinearCCMatrix
  
    ! ALPHA-parameter that controls the weight of the mass matrix in the
    ! core equation. =0.0 for stationary simulations.
    real(DP) :: dalpha = 0.0_DP
    
    ! THETA-parameter that controls the weight of the Stokes matrix
    ! in the core equation. =1.0 for stationary simulations.
    real(DP) :: dtheta = 0.0_DP
    
    ! GAMMA-parameter that controls the weight in front of the
    ! nonlinearity N(u). =1.0 for Navier-Stokes, =0.0 for Stokes equation.
    real(DP) :: dgamma = 0.0_DP

    ! ETA-parameter that switch the B-term on/off in the matrix.
    real(DP) :: deta = 0.0_DP
    
    ! TAU-parameter that switch the B^T-term on/off in the matrix.
    real(DP) :: dtau = 0.0_DP

    ! MU-parameter that switch the P-term on/off in the matrix.
    real(DP) :: dmu = 0.0_DP

    ! JOTA-parameter that switch the C-term on/off in the matrix.
    real(DP) :: djota = 0.0_DP

    ! KAPPA-parameter that switch the K-term on/off in the matrix.
    real(DP) :: dkappa = 0.0_DP

    ! Weight for the Newton matrix N*(u).
    ! = 0.0 deactivates the Newton part.
    real(DP) :: dnewton = 0.0_DP

    ! A structure referring to the physics of the problem.
    type(t_problem_physics), pointer :: p_rphysics => null()

    ! A pointer referring to the stabilisation to use
    type(t_problem_stabilisation), pointer :: p_rstabilisation => null()

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
    logical :: bvirtualTransposedD = .false.
    
    ! An object specifying the block discretisation
    ! (size of subvectors in the solution vector, trial/test functions,...).
    type(t_blockDiscretisation), pointer :: p_rdiscretisation => null()

    ! Pointer to template information on this level.
    type(t_asmTemplates), pointer :: p_rasmTempl => null()

    ! Pointer to dynamic information on this level (e.g. boundary conditions).
    ! This information usually changes with every timestep.
    type(t_dynamicLevelInfo), pointer :: p_rdynamicInfo => null()

  end type

!</typeblock>


!</types>

contains
  
! *****************************************************************

!<subroutine>

  subroutine ccmva_prepareViscoAssembly (rproblem,rphysics,rcollection,rvelocityvector)
  
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
  ! Global problem structure.
  type(t_problem), intent(in), target :: rproblem
  
  ! Physics of the problem
  type(t_problem_physics), intent(in), target :: rphysics
  
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
    rcollection%DquickAccess(1) = rphysics%dnu
    
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
      call getNonconstantViscosity (cterm,rdiscretisation, &
          nelements,npointsPerElement,Dpoints, &
          IdofsTest,rdomainIntSubset, &
          Dcoefficients,p_rvector,rcollection%p_rnextCollection)
    
    end select

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc_assembleMatrix (coperation,cmatrixType,rmatrix,&
      rnonlinearCCMatrix,rproblem,rvector,rfineMatrix)

!<description>
  ! This routine assembles a global matrix, i.e. it evaluates the nonlinear
  ! system matrix $A(\cdot)$ at a point $y$. The rnonlinearCCMatrix
  ! must contain a description of the matrix.
  ! The 'coperation' parameter tells the routine what to do.
  ! The destination matrix rmatrix, which receives the evaluated matrix
  ! A(rvector), is then set up or updated.
  !
  ! The parameters rvector and rfineMatrix are optional. rvector must be
  ! specified, if the nonlinearity is activated (parameter $\gamma\not=0$ in
  ! rnonlinearCCMatrix). This vector specifies the 'solution' where the
  ! nonlinearity $u\nabla u$ is evaluated.
  ! rfineMatrix allows to specify a matrix of a 'one level refined mesh'. This
  ! is usually used when setting up preconditioners over multiple levels.
  ! Specifying such a matrix allows the routine (depending on the discretisation)
  ! to include some special stabilisation terms into the matrix rmatrix.
  !
  ! The routine does not include any boundary conditions in the matrix.
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
  ! that references to the original gradient (B-)matrices from rnonlinearCCMatrix
  ! are written to rmatrix; matrix entries are not copied!
  ! (This can be used e.g. for setting up a matrix for building a defect
  !  vector without copying matrix data.)
  ! In this case, the caller MUST NOT CHANGE rmatrix in any way, otherwise
  ! the original (template) matrices would be changed!
  integer(I32), intent(in) :: coperation

  ! Type of matrix that should be set up in rmatrix. One of the CCMASM_MTP_xxxx
  ! constants.
  ! Usually, CCMASM_MTP_AUTOMATIC is used here. This will automatically determine
  ! the 'smallest possible' matrix structure fitting the needs of the input
  ! parameters. By specifying another matrix type, the caller can explicitly
  ! take influence on the general matrix structure.
  !
  ! If the matrix already exists and only its entries are to be computed,
  ! CCMASM_MTP_AUTOMATIC should be specified here.
  integer, intent(in) :: cmatrixType

  ! A t_nonlinearCCMatrix structure providing all necessary 'source' information
  ! about how to set up the matrix.
  !
  ! Note that if coperation=CCMASM_ALLOCxxxx is specified, p_rmatrixTemplateXXXX
  ! must be initialised as well as p_rdiscretisation!
  ! The new matrix is created based p_rmatrixTemplateXXXX as well as
  ! p_rdiscretisation. Memory is allocated automatically if it is missing.
  type(t_nonlinearCCMatrix), intent(in) :: rnonlinearCCMatrix

  ! Standard problem structure that defines all underlying parameters.
  type(t_problem), intent(inout), target :: rproblem

  ! OPTIONAL: If a nonlinearity is to be set up, this vector must be specified.
  ! It specifies where to evaluate the nonlinearity.
  type(t_vectorBlock), intent(in), optional :: rvector

  ! OPTIONAL: This parameter allows to specify a 'fine grid matrix'. This is
  ! usually done when assembling matrices on multiple levels. If specified, the
  ! routine will (if possible) try to include a level-dependent stabilisation
  ! term into the matrix (-> e.g. constant matrix restriction for nonparametric
  ! Rannacher-Turek element if cells are too anisotropic).
  type(t_matrixBlock), intent(in), optional :: rfineMatrix
  
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
    
    ballocate = .false.
    if ((rmatrix%NEQ .le. 0) .or. &
        iand(coperation,CCMASM_ALLOCMEM) .ne. 0) then
      ballocate = .true.
    end if
    
    ! What should we do? Allocate memory?
    if (ballocate) then
    
      ! Release the matrix if present.
      call lsysbl_releaseMatrix (rmatrix)
    
      ! Create a complete new matrix.
      call allocMatrix (cmatrixType,rnonlinearCCMatrix,rmatrix,rproblem)
    end if
   
    if (iand(coperation,CCMASM_COMPUTE) .ne. 0) then

      ! The system matrix looks like:
      !
      !    ( A11  A12  B1     )
      !    ( A21  A22  B2     )
      !    ( D1   D2   P   K  )
      !    (               C  )
      !
      ! Assemble the velocity submatrices
      !
      !    ( A11  A12   .    )
      !    ( A21  A22   .    )
      !    (  .    .    .    )
      !    (  .    .    .    )
      
      call assembleVelocityBlocks (&
          rnonlinearCCMatrix,rmatrix,rproblem,rvector,1.0_DP)
      
      ! Assemble the gradient submatrices
      !
      !    ( .    .    B1  )
      !    ( .    .    B2  )
      !    ( D1   D2   .   )
      
      call assembleGradientMatrices (rnonlinearCCMatrix,rmatrix,&
        iand(coperation,CMASM_QUICKREFERENCES) .ne. 0,rproblem)

      ! Assemble the pressure submatrix (if it exists)
      !
      !    ( .    .    .  . )
      !    ( .    .    .  . )
      !    ( .    .    P  . )
      !    ( .    .    .  . )
      
      call assemblePressureMatrix (rnonlinearCCMatrix,rmatrix,rproblem)

      ! Assemble the K submatrix (if it exists)
      !
      !    ( .    .    .   . )
      !    ( .    .    .   . )
      !    ( .    .    .   K )
      !    ( .    .    .   . )
      
      call assembleKappaMatrix (rnonlinearCCMatrix,rmatrix,rproblem,rvector)

      ! Assemble the C submatrix (if it exists)
      !
      !    ( .    .    .   . )
      !    ( .    .    .   . )
      !    ( .    .    .   . )
      !    ( .    .    .   C )
      
      call assembleConcentrationMatrix (rnonlinearCCMatrix,rmatrix,rproblem,rvector)

      ! 2.) Initialise the weights for the B-matrices
      !
      !    (  .    .   B1  )
      !    (  .    .   B2  )
      !    ( D1   D2    .  )
      
      rmatrix%RmatrixBlock(1,3)%dscaleFactor = rnonlinearCCMatrix%deta
      rmatrix%RmatrixBlock(2,3)%dscaleFactor = rnonlinearCCMatrix%deta
      
      rmatrix%RmatrixBlock(3,1)%dscaleFactor = rnonlinearCCMatrix%dtau
      rmatrix%RmatrixBlock(3,2)%dscaleFactor = rnonlinearCCMatrix%dtau
      
      ! Matrix restriction
      ! ---------------------------------------------------
      !
      ! For the construction of matrices on lower levels, call the matrix
      ! restriction. In case we have a uniform discretisation with Q1~,
      ! iadaptivematrix is <> 0 and so this will rebuild some matrix entries
      ! by a Galerkin approach using constant prolongation/restriction.
      ! This helps to stabilise the solver if there are elements in the
      ! mesh with high aspect ratio.
      if (present(rfineMatrix)) then
        call mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(1,1), &
            rmatrix%RmatrixBlock(1,1), &
            rnonlinearCCMatrix%iadaptiveMatrices, rnonlinearCCMatrix%dadmatthreshold)
            
        if (.not. lsyssc_isMatrixContentShared(rmatrix%RmatrixBlock(2,2))) then
          call mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(2,2), &
              rmatrix%RmatrixBlock(2,2), &
              rnonlinearCCMatrix%iadaptiveMatrices, rnonlinearCCMatrix%dadmatthreshold)
        end if
      end if
      
      !call matio_writeBlockMatrixHR (rmatrix, 'matrix',&
      !    .true., 0, 'matrix.txt', '(E20.10)')
      
    end if
    
  contains
  
    ! -----------------------------------------------------
  
    subroutine allocMatrix (cmatrixType,rnonlinearCCMatrix,rmatrix,rproblem)
    
    ! Allocates memory for the system matrix. rnonlinearCCMatrix provides information
    ! about the submatrices that are 'plugged into' rmatrix.
    ! Therefore, before this routine is called, rnonlinearCCMatrix must have been set up.

    ! Type of matrix that should be set up in rmatrix. One of the CCMASM_MTP_xxxx
    ! constants.
    integer, intent(in) :: cmatrixType

    ! A t_nonlinearCCMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix.
    type(t_nonlinearCCMatrix), intent(in), target :: rnonlinearCCMatrix

    ! A block matrix that receives the basic system matrix.
    type(t_matrixBlock), intent(inout) :: rmatrix
    
    ! Standard problem structure that defines all underlying parameters.
    type(t_problem), intent(inout), target :: rproblem

      ! local variables
      logical :: bdecoupled, bfulltensor

      ! A pointer to the system matrix and the RHS/solution vectors.
      type(t_matrixScalar), pointer :: p_rmatrixTemplateFEM

      ! A pointer to the discretisation structure with the data.
      type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
      ! Determine the shape of the matrix
      bdecoupled = cmatrixType .eq. CCMASM_MTP_DECOUPLED
      bfulltensor = cmatrixType .eq. CCMASM_MTP_FULLTENSOR
      
      if (cmatrixType .eq. CCMASM_MTP_AUTOMATIC) then
        ! Should we assemble Newton or the deformation tensor?
        ! If yes, we have a full-tensor matrix.
        bfulltensor = (rnonlinearCCMatrix%dnewton .ne. 0.0_DP) .or. &
                      (rnonlinearCCMatrix%p_rphysics%isubequation .ne. 0)
      end if
    
      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rnonlinearCCMatrix%p_rdiscretisation
      
      ! Get a pointer to the template FEM matrix. If that does not exist,
      ! take the Stokes matrix as template.
      p_rmatrixTemplateFEM => rnonlinearCCMatrix%p_rasmTempl%rmatrixTemplateFEM
      if (.not. associated(p_rmatrixTemplateFEM)) &
        p_rmatrixTemplateFEM => rnonlinearCCMatrix%p_rasmTempl%rmatrixStokes
      if (.not. associated(p_rmatrixTemplateFEM)) then
        call output_line ('Cannot set up A matrices in system matrix!', &
            OU_CLASS_ERROR,OU_MODE_STD,'allocMatrix')
        call sys_halt()
      end if

      ! Initialise the block matrix with default values based on
      ! the discretisation.
      if (associated(p_rdiscretisation)) then
        call lsysbl_createMatBlockByDiscr (p_rdiscretisation,rmatrix)
      else
        ! No discretisation structure; create the matrix directly as 3x3 matrix.
        call lsysbl_createEmptyMatrix (rmatrix,NDIM2D+1)
      end if
        
      ! Let us consider the global system in detail. The standard matrix It has
      ! roughly the following shape:
      !
      !    ( A11       B1    ) = ( A11  A12  A13  A14 )
      !    (      A22  B2    )   ( A21  A22  A23  A24 )
      !    ( D1   D2   P   K )   ( A31  A32  A33  A34 )
      !    (               C )   ( A41  A42  A43  A44 )
      !
      ! All matrices may have multiplication factors in their front.
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
      call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                  rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          
      if (.not. bdecoupled .and. .not. bfulltensor) then
           
        ! If X- and Y-velocity is to be treated in a 'coupled' way, the matrix
        ! A22 is identical to A11! So mirror A11 to A22 sharing the
        ! structure and the content.
        call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(1,1),&
                    rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
                    
      else
      
        ! Otherwise, create another copy of the template matrix.
        call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                    rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
                    
      end if
      
      ! Manually change the discretisation structure of the Y-velocity
      ! matrix to the Y-discretisation structure.
      ! Ok, we use the same discretisation structure for both, X- and Y-velocity,
      ! so this is not really necessary - we do this for sure...
      call lsyssc_assignDiscrDirectMat (rmatrix%RmatrixBlock(2,2),&
          p_rdiscretisation%RspatialDiscr(2))

      ! A 'full tensor matrix' consists also of blocks A12 and A21.
      if (bfulltensor) then

        ! We have a matrix in the following shape:
        !
        !    ( A11  A12  B1  . )
        !    ( A21  A22  B2  . )
        !    ( D1   D2   .   . )
        !    ( .    .    .   . )
        !
        ! Create A12 and A21.
      
        if (rmatrix%RmatrixBlock(1,2)%cmatrixFormat &
            .eq. LSYSSC_MATRIXUNDEFINED) then
            
          call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM, &
            rmatrix%RmatrixBlock(1,2), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
            
          ! Allocate memory for the entries; do not initialise the memory.
          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
          ! zero.
          ! CALL lsyssc_allocEmptyMatrix (&
          !     p_rmatrixPreconditioner%RmatrixBlock(1,2),LSYSSC_SETM_UNDEFINED)
            
        end if

        if (rmatrix%RmatrixBlock(2,1)%cmatrixFormat &
            .eq. LSYSSC_MATRIXUNDEFINED) then
            
          ! Create a new matrix A21 in memory. create a new matrix
          ! using the template FEM matrix...
          call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM, &
            rmatrix%RmatrixBlock(2,1), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
            
          ! Allocate memory for the entries; do not initialise the memory.
          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
          ! zero.
          ! CALL lsyssc_allocEmptyMatrix (&
          !     p_rmatrixPreconditioner%RmatrixBlock(2,1),LSYSSC_SETM_UNDEFINED)
            
        end if
        
      end if

      ! The B1/B2 matrices exist up to now only in our local problem structure.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2 with those B1/B2 of the
      ! block matrix, while we create empty space for the entries.
      ! Later, the B-matrices are copied into here and modified for boundary
      ! conditions.
      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixB1, &
          rmatrix%RmatrixBlock(1,3),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixB2, &
          rmatrix%RmatrixBlock(2,3),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        
      ! Now, prepare D1 and D2.
      ! The flag bvirtualTransposedD in the structure decides on whether
      ! these matrices are created by the default Di matrices or by
      ! virtually transposing the B-matrices. (This may be need by some VANKA
      ! variants in the preconditioner e.g.)
      if (rnonlinearCCMatrix%bvirtualTransposedD) then
        call lsyssc_transposeMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixD1T,&
            rmatrix%RmatrixBlock(3,1),LSYSSC_TR_VIRTUAL)
            
        call lsyssc_transposeMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixD2T,&
            rmatrix%RmatrixBlock(3,2),LSYSSC_TR_VIRTUAL)
      else
        call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixD1, &
            rmatrix%RmatrixBlock(3,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

        call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixD2, &
            rmatrix%RmatrixBlock(3,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      end if

      ! Include a matrix for the pressure
      !
      !    ( .    .    .    .  )
      !    ( .    .    .    .  )
      !    ( .    .    C    .  )
      !    ( .    .    .    .  )
      !
      ! which may be used for stabilisation or other features.
      ! This submatrix will be deactived by setting the scaling factor
      ! to 0.
      call lsyssc_duplicateMatrix (&
          rnonlinearCCMatrix%p_rasmTempl%rmatrixTemplateFEMPressure,&
          rmatrix%RmatrixBlock(3,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          
      call lsyssc_clearMatrix (rmatrix%RmatrixBlock(3,3))

      ! Include a K-matrix
      !
      !    ( .    .    .    .  )
      !    ( .    .    .    .  )
      !    ( .    .    .    K  )
      !    ( .    .    .    .  )
      !
      ! which may be used for stabilisation or other features.
      ! This submatrix will be deactived by setting the scaling factor
      ! to 0.
      call lsyssc_duplicateMatrix (&
          rnonlinearCCMatrix%p_rasmTempl%rmatrixTemplateK,&
          rmatrix%RmatrixBlock(3,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          
      call lsyssc_clearMatrix (rmatrix%RmatrixBlock(3,4))

      ! Include a concentration matrix
      !
      !    ( .    .    .    .  )
      !    ( .    .    .    .  )
      !    ( .    .    .    .  )
      !    ( .    .    .    C  )
      !
      ! which may be used for stabilisation or other features.
      ! This submatrix will be deactived by setting the scaling factor
      ! to 0.
      call lsyssc_duplicateMatrix (&
          rnonlinearCCMatrix%p_rasmTempl%rmatrixTemplateFEMConcentration,&
          rmatrix%RmatrixBlock(4,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          
      call lsyssc_clearMatrix (rmatrix%RmatrixBlock(4,4))

      ! That is it, all submatrices are set up.
        
    end subroutine
    
    ! -----------------------------------------------------

    subroutine assembleVelocityBlocks (rnonlinearCCMatrix,rmatrix,rproblem,&
        rvelocityvector,dvectorWeight)
        
    ! Assembles the velocity matrix in the block matrix rmatrix at position (1,1):
    !
    ! rmatrix := dalpha*M + dtheta*Laplace + dgamma*N(p_rvector) +
    !            dnewton*N*(p_rvector)
    
    ! A t_nonlinearCCMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix.
    type(t_nonlinearCCMatrix), intent(in) :: rnonlinearCCMatrix
    
    ! Block matrix where the 2x2-velocity submatrix should be assembled
    type(t_matrixBlock), intent(inout) :: rmatrix
    
    ! Standard problem structure that defines all underlying parameters.
    type(t_problem), intent(inout), target :: rproblem

    ! Velocity vector for the nonlinearity. Must be specified if
    ! GAMMA <> 0; can be omitted if GAMMA=0.
    type(t_vectorBlock), target, optional :: rvelocityvector
    
    ! Weight for the velocity vector; standard = -1.
    real(DP), intent(in), optional :: dvectorWeight
    
    ! local variables
    logical :: bshared
    integer :: iupwind
    type(t_convUpwind) :: rupwind
    type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
    type(t_convStreamDiff2) :: rstreamlineDiffusion2
    type(t_jumpStabilisation) :: rjumpStabil
    real(DP) :: dvecWeight
    type(t_collection) :: rcollection
    integer, dimension(:), pointer :: p_IedgesDirichletBC
    
      ! Standard value for dvectorWeight is = -1.
      dvecWeight = -1.0_DP
      if (present(dvectorWeight)) dvecWeight = dvectorWeight
    
      ! Is A11=A22 physically?
      bshared = lsyssc_isMatrixContentShared(&
                    rmatrix%RmatrixBlock(1,1),&
                    rmatrix%RmatrixBlock(2,2))
                    
      ! Allocate memory if necessary. Normally this should not be necessary...
      ! A11:
      if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(1,1))) then
        call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(1,1),LSYSSC_SETM_UNDEFINED)
      end if
    
      ! A22:
      if (.not. bshared) then
        if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(2,2))) then
          call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(2,2),LSYSSC_SETM_UNDEFINED)
        end if
      end if

      ! A12/ A21:
      if (lsysbl_isSubmatrixPresent (rmatrix,1,2)) then
        if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(1,2))) then
          call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(1,2),LSYSSC_SETM_UNDEFINED)
        end if
        if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(2,1))) then
          call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(2,1),LSYSSC_SETM_UNDEFINED)
        end if
      end if
    
      ! ---------------------------------------------------
      ! Plug in the mass matrix?
      if (rnonlinearCCMatrix%dalpha .ne. 0.0_DP) then
       
        ! Allocate memory if necessary. Normally this should not be necessary...
        if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(1,1))) then
          call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(1,1),LSYSSC_SETM_UNDEFINED)
        end if
      
        call lsyssc_matrixLinearComb (&
            rnonlinearCCMatrix%p_rasmTempl%rmatrixMass,rnonlinearCCMatrix%dalpha,&
            rmatrix%RmatrixBlock(1,1),0.0_DP,&
            rmatrix%RmatrixBlock(1,1),&
            .false.,.false.,.true.,.true.)
            
        if (.not. bshared) then

          ! Allocate memory if necessary. Normally this should not be necessary...
          if (.not. lsyssc_hasMatrixContent (rmatrix%RmatrixBlock(2,2))) then
            call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(2,2),LSYSSC_SETM_UNDEFINED)
          end if

          call lsyssc_matrixLinearComb (&
              rnonlinearCCMatrix%p_rasmTempl%rmatrixMass,rnonlinearCCMatrix%dalpha,&
              rmatrix%RmatrixBlock(2,2),0.0_DP,&
              rmatrix%RmatrixBlock(2,2),&
              .false.,.false.,.true.,.true.)
        end if
        
      else
      
        ! Otherwise, initialise the basic matrix with 0.0
        call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,1))
        
        if (.not. bshared) then
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,2))
        end if
        
      end if
      
      ! If the submatrices A12 and A21 exist, fill them with zero.
      ! If they do not exist, we do not have to do anything.
      if ((rnonlinearCCMatrix%dnewton .ne. 0.0_DP) .or. &
          (rnonlinearCCMatrix%p_rphysics%isubequation .ne. 0)) then
        rmatrix%RmatrixBlock(1,2)%dscaleFactor = 1.0_DP
        rmatrix%RmatrixBlock(2,1)%dscaleFactor = 1.0_DP
      else
        rmatrix%RmatrixBlock(1,2)%dscaleFactor = 0.0_DP
        rmatrix%RmatrixBlock(2,1)%dscaleFactor = 0.0_DP
      end if
      
      if (lsysbl_isSubmatrixPresent (rmatrix,1,2)) then
        call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,2))
        call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,1))
      end if
        
      ! ---------------------------------------------------
      ! Plug in the Stokes matrix?
      if (rnonlinearCCMatrix%dtheta .ne. 0.0_DP) then
        ! Plug in the Stokes matrix in case of the gradient tensor.
        ! In case of the deformation tensor ir nonconstant viscosity,
        ! that is done during the assembly of the nonlinearity.
        if ((rnonlinearCCMatrix%p_rphysics%isubequation .eq. 0) .and. &
            (rnonlinearCCMatrix%p_rphysics%cviscoModel .eq. 0)) then
          call lsyssc_matrixLinearComb (&
              rnonlinearCCMatrix%p_rasmTempl%rmatrixStokes,rnonlinearCCMatrix%dtheta,&
              rmatrix%RmatrixBlock(1,1),1.0_DP,&
              rmatrix%RmatrixBlock(1,1),&
              .false.,.false.,.true.,.true.)
              
          if (.not. bshared) then
            call lsyssc_matrixLinearComb (&
                rnonlinearCCMatrix%p_rasmTempl%rmatrixStokes,rnonlinearCCMatrix%dtheta,&
                rmatrix%RmatrixBlock(2,2),1.0_DP,&
                rmatrix%RmatrixBlock(2,2),&
                .false.,.false.,.true.,.true.)
          end if
        end if
      end if
      
      ! ---------------------------------------------------
      ! That was easy -- the adventure begins now... The nonlinearity!
      if ((rnonlinearCCMatrix%dgamma .ne. 0.0_DP) .or. &
          (rnonlinearCCMatrix%p_rphysics%isubequation .ne. 0) .or. &
          (rnonlinearCCMatrix%p_rphysics%cviscoModel .ne. 0)) then
      
        if (.not. present(rvelocityvector)) then
          call output_line ('Velocity vector not present!', &
                             OU_CLASS_ERROR,OU_MODE_STD,'cc_assembleMatrix')
          stop
        end if
      
        select case (rnonlinearCCMatrix%p_rstabilisation%iupwind)
        case (CCMASM_STAB_STREAMLINEDIFF)
          ! Streamline diffusion.

          ! Set up the SD structure for the creation of the defect.
          ! There is not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearCCMatrix%p_rphysics%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion%dupsam = rnonlinearCCMatrix%p_rstabilisation%dupsam
          
          ! Set calculation method for local H
          rstreamlineDiffusion%clocalH = rnonlinearCCMatrix%p_rstabilisation%clocalH
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion%ddelta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newton part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = rnonlinearCCMatrix%dnewton
          
          ! Call the SD method to calculate the nonlinearity.
          call conv_streamlineDiffusionBlk2d (&
                              rvelocityvector, rvelocityvector, &
                              dvecWeight, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODMATRIX, &
                              rmatrix)
                              
        case (CCMASM_STAB_STREAMLINEDIFF2)

          ! Set up the SD structure for the creation of the defect.
          ! There is not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearCCMatrix%p_rphysics%dnu

          ! Probably, we have nonconstant viscosity.
          ! In that case, init a collection structure for a callback
          ! routine that specifies the viscosity
          if (rnonlinearCCMatrix%p_rphysics%cviscoModel .ne. 0) then
            rstreamlineDiffusion2%bconstnu = .false.
            
            ! Assemble at least the Stokes matrix.
            ! For a more complicated formulation, dbeta/dbetaT
            ! may be changed later.
            rstreamlineDiffusion2%dbeta = 1.0_DP
            
            ! Prepare the collection. The "next" collection points to the user defined
            ! collection.
            rcollection%p_rnextCollection => rproblem%rcollection
            call ccmva_prepareViscoAssembly (rproblem,rnonlinearCCMatrix%p_rphysics,&
                rcollection,rvelocityVector)
            
          end if
          
          ! Set stabilisation parameter
          rstreamlineDiffusion2%dupsam = rnonlinearCCMatrix%p_rstabilisation%dupsam
          
          ! Set calculation method for local H
          rstreamlineDiffusion2%clocalH = rnonlinearCCMatrix%p_rstabilisation%clocalH
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion2%ddelta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newton part; =0 deactivates Newton.
          rstreamlineDiffusion2%dnewton = rnonlinearCCMatrix%dnewton
          
          ! Assemble the deformation tensor?
          if (rnonlinearCCMatrix%p_rphysics%isubequation .ne. 0) then
            rstreamlineDiffusion2%dbeta = 0.5_DP
            rstreamlineDiffusion2%dbetaT = 0.5_DP
          end if
          
          ! Initialise the user defined collection for the assembly.
          call cc_initCollectForAssembly (rproblem,rproblem%rcollection)

          ! Call the SD method to calculate the nonlinearity.
          call conv_streamDiff2Blk2dMat (rstreamlineDiffusion2,rmatrix,rvelocityvector,&
              ffunctionViscoModel,rcollection)
           
          
              
          ! That is it.
          call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

        case (CCMASM_STAB_UPWIND)
          ! Set up the upwind structure for the creation of the defect.
          ! There is not much to do, only initialise the viscosity...
          rupwind%dnu = rnonlinearCCMatrix%p_rphysics%dnu
          
          ! Set stabilisation parameter
          rupwind%dupsam = rnonlinearCCMatrix%p_rstabilisation%dupsam

          ! Matrix weight
          rupwind%dtheta = rnonlinearCCMatrix%dgamma
          
          if (rnonlinearCCMatrix%dnewton .ne. 0.0_DP) then
            call output_line ('Warning: Upwind does not support assembly '&
                //'of the Newton matrix!',OU_CLASS_TRACE1)
          end if
          
          ! Call the upwind method to calculate the nonlinear matrix.
          call conv_upwind2d (rvelocityvector, rvelocityvector, &
                              dvecWeight, 0.0_DP,&
                              rupwind, CONV_MODMATRIX, &
                              rmatrix%RmatrixBlock(1,1))
                              
          if (.not. bshared) then
            ! Modify also the matrix block (2,2)
            call conv_upwind2d (rvelocityvector, rvelocityvector, &
                                dvecWeight, 0.0_DP,&
                                rupwind, CONV_MODMATRIX, &
                                rmatrix%RmatrixBlock(2,2))
          end if

        case (CCMASM_STAB_EDGEORIENTED)
          ! Jump stabilisation.

          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.

          ! There is not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearCCMatrix%p_rphysics%dnu
          
          ! Set stabilisation parameter to zero for central difference
          rstreamlineDiffusion%dupsam = 0.0_DP
          
          ! Set calculation method for local H
          rstreamlineDiffusion%clocalH = rnonlinearCCMatrix%p_rstabilisation%clocalH
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion%ddelta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newton part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = rnonlinearCCMatrix%dnewton
          
          ! Call the SD method to calculate the nonlinearity.
          call conv_streamlineDiffusionBlk2d (&
                              rvelocityvector, rvelocityvector, &
                              dvecWeight, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODMATRIX, &
                              rmatrix)
                              
          ! Set up the jump stabilisation structure.
          ! There is not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rnonlinearCCMatrix%p_rphysics%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgamma = rnonlinearCCMatrix%p_rstabilisation%dupsam
          rjumpStabil%dgammastar = rnonlinearCCMatrix%p_rstabilisation%dupsamstar
          rjumpStabil%deojEdgeExp = rnonlinearCCMatrix%p_rstabilisation%deojEdgeExp
          
          ! Matrix weight
          rjumpStabil%dtheta = rnonlinearCCMatrix%dtheta

          ! Cubature formula
          rjumpStabil%ccubType = rnonlinearCCMatrix%p_rstabilisation%ccubEOJ

          ! Call the jump stabilisation technique to stabilise that stuff.
          ! We can assemble the jump part any time as it is independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
              rjumpStabil, CONV_MODMATRIX, rmatrix%RmatrixBlock(1,1),&
              rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil)

          if (.not. bshared) then
            call conv_jumpStabilisation2d (&
                rjumpStabil, CONV_MODMATRIX,rmatrix%RmatrixBlock(2,2),&
              rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil)
          end if
          
          ! Subtract the EOJ matrix for the Dirichlet boundary conditions.
          if (rnonlinearCCMatrix%p_rdynamicInfo%nedgesDirichletBC .ne. 0) then
            rjumpStabil%dnu = rnonlinearCCMatrix%p_rphysics%dnu
            rjumpStabil%dgamma = rnonlinearCCMatrix%p_rstabilisation%dupsam
            rjumpStabil%dgammastar = rnonlinearCCMatrix%p_rstabilisation%dupsamstar
            rjumpStabil%deojEdgeExp = rnonlinearCCMatrix%p_rstabilisation%deojEdgeExp
            rjumpStabil%dtheta = -rnonlinearCCMatrix%dtheta
            rjumpStabil%ccubType = rnonlinearCCMatrix%p_rstabilisation%ccubEOJ
            call storage_getbase_int(rnonlinearCCMatrix%p_rdynamicInfo%hedgesDirichletBC,&
                p_IedgesDirichletBC)
            call conv_jumpStabilisation2d (&
                rjumpStabil, CONV_MODMATRIX, rmatrix%RmatrixBlock(1,1),&
                rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil,&
                InodeList=p_IedgesDirichletBC)

            if (.not. bshared) then
              call conv_jumpStabilisation2d (&
                  rjumpStabil, CONV_MODMATRIX,rmatrix%RmatrixBlock(2,2),&
                rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil,&
                InodeList=p_IedgesDirichletBC)
            end if
          end if

        case (CCMASM_STAB_FASTEDGEORIENTED)
          ! Jump stabilisation with precomputed matrix.

          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.

          ! There is not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearCCMatrix%p_rphysics%dnu
          
          ! Set stabilisation parameter to zero for central difference
          rstreamlineDiffusion%dupsam = 0.0_DP
          
          ! Set calculation method for local H
          rstreamlineDiffusion%clocalH = rnonlinearCCMatrix%p_rstabilisation%clocalH
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion%ddelta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newton part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = rnonlinearCCMatrix%dnewton
          
          ! Call the SD method to calculate the nonlinearity.
          call conv_streamlineDiffusionBlk2d (&
                              rvelocityvector, rvelocityvector, &
                              dvecWeight, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODMATRIX, &
                              rmatrix)
        
          ! Sum up the precomputed edge stabilisation matrix.
          call lsyssc_matrixLinearComb (&
              rnonlinearCCMatrix%p_rasmTempl%rmatrixStabil,rnonlinearCCMatrix%dtheta,&
              rmatrix%RmatrixBlock(1,1),1.0_DP,&
              rmatrix%RmatrixBlock(1,1),&
              .false.,.false.,.true.,.true.)
          
          if (.not. bshared) then
            call lsyssc_matrixLinearComb (&
                rnonlinearCCMatrix%p_rasmTempl%rmatrixStabil,rnonlinearCCMatrix%dtheta,&
                rmatrix%RmatrixBlock(2,2),1.0_DP,&
                rmatrix%RmatrixBlock(2,2),&
                .false.,.false.,.true.,.true.)
          end if
          
          ! Subtract the EOJ matrix for the Dirichlet boundary conditions.
          if (rnonlinearCCMatrix%p_rdynamicInfo%nedgesDirichletBC .ne. 0) then
            rjumpStabil%dnu = rnonlinearCCMatrix%p_rphysics%dnu
            rjumpStabil%dgamma = rnonlinearCCMatrix%p_rstabilisation%dupsam
            rjumpStabil%dgammastar = rnonlinearCCMatrix%p_rstabilisation%dupsamstar
            rjumpStabil%deojEdgeExp = rnonlinearCCMatrix%p_rstabilisation%deojEdgeExp
            rjumpStabil%dtheta = -rnonlinearCCMatrix%dtheta
            rjumpStabil%ccubType = rnonlinearCCMatrix%p_rstabilisation%ccubEOJ
            call storage_getbase_int(rnonlinearCCMatrix%p_rdynamicInfo%hedgesDirichletBC,&
                p_IedgesDirichletBC)
            call conv_jumpStabilisation2d (&
                rjumpStabil, CONV_MODMATRIX, rmatrix%RmatrixBlock(1,1),&
                rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil,&
                InodeList=p_IedgesDirichletBC)

            if (.not. bshared) then
              call conv_jumpStabilisation2d (&
                  rjumpStabil, CONV_MODMATRIX,rmatrix%RmatrixBlock(2,2),&
                rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil,&
                InodeList=p_IedgesDirichletBC)
            end if
          end if

        case (CCMASM_STAB_EDGEORIENTED2)
          ! Jump stabilisation.

          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! Set up the SD structure for the creation of the defect.
          ! There is not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearCCMatrix%p_rphysics%dnu

          ! Probably, we have nonconstant viscosity.
          ! In that case, init a collection structure for a callback
          ! routine that specifies the viscosity
          if (rnonlinearCCMatrix%p_rphysics%cviscoModel .ne. 0) then
            rstreamlineDiffusion2%bconstnu = .false.
            
            ! Assemble at least the Stokes matrix.
            ! For a more complicated formulation, dbeta/dbetaT
            ! may be changed later.
            rstreamlineDiffusion2%dbeta = 1.0_DP

            ! Prepare the collection. The "next" collection points to the user defined
            ! collection.
            rcollection%p_rnextCollection => rproblem%rcollection
            call ccmva_prepareViscoAssembly (rproblem,rnonlinearCCMatrix%p_rphysics,&
                rcollection,rvelocityVector)

          end if
          
          ! Set UPSAM=0 to deactivate the stabilisation
          rstreamlineDiffusion2%dupsam = 0.0_DP
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion2%ddelta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newton part; =0 deactivates Newton.
          rstreamlineDiffusion2%dnewton = rnonlinearCCMatrix%dnewton
          
          ! Assemble the deformation tensor?
          if (rnonlinearCCMatrix%p_rphysics%isubequation .ne. 0) then
            rstreamlineDiffusion2%dbeta = 0.5_DP
            rstreamlineDiffusion2%dbetaT = 0.5_DP
          end if

          ! Initialise the user defined collection for the assembly.
          call cc_initCollectForAssembly (rproblem,rproblem%rcollection)

          ! Call the SD method to calculate the nonlinearity.
          call conv_streamDiff2Blk2dMat (rstreamlineDiffusion2,rmatrix,rvelocityvector,&
              ffunctionViscoModel,rcollection)

          ! Set up the jump stabilisation structure.
          ! There is not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rnonlinearCCMatrix%p_rphysics%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgamma = rnonlinearCCMatrix%p_rstabilisation%dupsam
          rjumpStabil%dgammastar = rnonlinearCCMatrix%p_rstabilisation%dupsamstar
          rjumpStabil%deojEdgeExp = rnonlinearCCMatrix%p_rstabilisation%deojEdgeExp
          
          ! Matrix weight
          rjumpStabil%dtheta = rnonlinearCCMatrix%dtheta

          ! Cubature formula
          rjumpStabil%ccubType = rnonlinearCCMatrix%p_rstabilisation%ccubEOJ

          ! Call the jump stabilisation technique to stabilise that stuff.
          ! We can assemble the jump part any time as it is independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
              rjumpStabil, CONV_MODMATRIX, rmatrix%RmatrixBlock(1,1),&
              rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil)

          if (.not. bshared) then
            call conv_jumpStabilisation2d (&
                rjumpStabil, CONV_MODMATRIX,rmatrix%RmatrixBlock(2,2),&
              rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil)
          end if

          ! Subtract the EOJ matrix for the Dirichlet boundary conditions.
          if (rnonlinearCCMatrix%p_rdynamicInfo%nedgesDirichletBC .ne. 0) then
            rjumpStabil%dnu = rnonlinearCCMatrix%p_rphysics%dnu
            rjumpStabil%dgamma = rnonlinearCCMatrix%p_rstabilisation%dupsam
            rjumpStabil%dgammastar = rnonlinearCCMatrix%p_rstabilisation%dupsamstar
            rjumpStabil%deojEdgeExp = rnonlinearCCMatrix%p_rstabilisation%deojEdgeExp
            rjumpStabil%dtheta = -rnonlinearCCMatrix%dtheta
            rjumpStabil%ccubType = rnonlinearCCMatrix%p_rstabilisation%ccubEOJ
            call storage_getbase_int(rnonlinearCCMatrix%p_rdynamicInfo%hedgesDirichletBC,&
                p_IedgesDirichletBC)
            call conv_jumpStabilisation2d (&
                rjumpStabil, CONV_MODMATRIX, rmatrix%RmatrixBlock(1,1),&
                rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil,&
                InodeList=p_IedgesDirichletBC)

            if (.not. bshared) then
              call conv_jumpStabilisation2d (&
                  rjumpStabil, CONV_MODMATRIX,rmatrix%RmatrixBlock(2,2),&
                rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil,&
                InodeList=p_IedgesDirichletBC)
            end if
          end if

          ! That is it.
          call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

        case default
          call output_line ('Don''t know how to set up nonlinearity!?!', &
              OU_CLASS_ERROR,OU_MODE_STD,'assembleVelocityBlocks')
          call sys_halt()
        
        end select

      else
      
        ! That is the Stokes-case. Jump stabilisation is possible...
      
        select case (rnonlinearCCMatrix%p_rstabilisation%iupwind)
        case (CCMASM_STAB_EDGEORIENTED,CCMASM_STAB_EDGEORIENTED2)
          ! Jump stabilisation.
        
          ! Set up the jump stabilisation structure.
          ! There is not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rnonlinearCCMatrix%p_rphysics%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgamma = rnonlinearCCMatrix%p_rstabilisation%dupsam
          rjumpStabil%dgammastar = rnonlinearCCMatrix%p_rstabilisation%dupsamstar
          rjumpStabil%deojEdgeExp = rnonlinearCCMatrix%p_rstabilisation%deojEdgeExp
          
          ! Matrix weight
          rjumpStabil%dtheta = rnonlinearCCMatrix%dtheta

          ! Cubature formula
          rjumpStabil%ccubType = rnonlinearCCMatrix%p_rstabilisation%ccubEOJ

          ! Call the jump stabilisation technique to stabilise that stuff.
          ! We can assemble the jump part any time as it is independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
              rjumpStabil, CONV_MODMATRIX,rmatrix%RmatrixBlock(1,1),&
              rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil)

          if (.not. bshared) then
            call conv_jumpStabilisation2d (&
                rjumpStabil,CONV_MODMATRIX,rmatrix%RmatrixBlock(2,2),&
                rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil)
          end if

          ! Subtract the EOJ matrix for the Dirichlet boundary conditions.
          if (rnonlinearCCMatrix%p_rdynamicInfo%nedgesDirichletBC .ne. 0) then
            rjumpStabil%dnu = rnonlinearCCMatrix%p_rphysics%dnu
            rjumpStabil%dgamma = rnonlinearCCMatrix%p_rstabilisation%dupsam
            rjumpStabil%dgammastar = rnonlinearCCMatrix%p_rstabilisation%dupsamstar
            rjumpStabil%deojEdgeExp = rnonlinearCCMatrix%p_rstabilisation%deojEdgeExp
            rjumpStabil%dtheta = -rnonlinearCCMatrix%dtheta
            rjumpStabil%ccubType = rnonlinearCCMatrix%p_rstabilisation%ccubEOJ
            call storage_getbase_int(rnonlinearCCMatrix%p_rdynamicInfo%hedgesDirichletBC,&
                p_IedgesDirichletBC)
            call conv_jumpStabilisation2d (&
                rjumpStabil, CONV_MODMATRIX, rmatrix%RmatrixBlock(1,1),&
                rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil,&
                InodeList=p_IedgesDirichletBC)

            if (.not. bshared) then
              call conv_jumpStabilisation2d (&
                  rjumpStabil, CONV_MODMATRIX,rmatrix%RmatrixBlock(2,2),&
                rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil,&
                InodeList=p_IedgesDirichletBC)
            end if
          end if

        case (CCMASM_STAB_FASTEDGEORIENTED)
          ! Fast Jump stabilisation. Precomputed matrix.
        
          ! Sum up the precomputed edge stabilisation matrix.
          call lsyssc_matrixLinearComb (&
              rnonlinearCCMatrix%p_rasmTempl%rmatrixStabil,rnonlinearCCMatrix%dtheta,&
              rmatrix%RmatrixBlock(1,1),1.0_DP,&
              rmatrix%RmatrixBlock(1,1),&
              .false.,.false.,.true.,.true.)
          
          if (.not. bshared) then
            call lsyssc_matrixLinearComb (&
                rnonlinearCCMatrix%p_rasmTempl%rmatrixStabil,rnonlinearCCMatrix%dtheta,&
                rmatrix%RmatrixBlock(2,2),1.0_DP,&
                rmatrix%RmatrixBlock(2,2),&
                .false.,.false.,.true.,.true.)
          end if

          ! Subtract the EOJ matrix for the Dirichlet boundary conditions.
          if (rnonlinearCCMatrix%p_rdynamicInfo%nedgesDirichletBC .ne. 0) then
            rjumpStabil%dnu = rnonlinearCCMatrix%p_rphysics%dnu
            rjumpStabil%dgamma = rnonlinearCCMatrix%p_rstabilisation%dupsam
            rjumpStabil%dgammastar = rnonlinearCCMatrix%p_rstabilisation%dupsamstar
            rjumpStabil%deojEdgeExp = rnonlinearCCMatrix%p_rstabilisation%deojEdgeExp
            rjumpStabil%dtheta = -rnonlinearCCMatrix%dtheta
            rjumpStabil%ccubType = rnonlinearCCMatrix%p_rstabilisation%ccubEOJ
            call storage_getbase_int(rnonlinearCCMatrix%p_rdynamicInfo%hedgesDirichletBC,&
                p_IedgesDirichletBC)
            call conv_jumpStabilisation2d (&
                rjumpStabil, CONV_MODMATRIX, rmatrix%RmatrixBlock(1,1),&
                rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil,&
                InodeList=p_IedgesDirichletBC)

            if (.not. bshared) then
              call conv_jumpStabilisation2d (&
                  rjumpStabil, CONV_MODMATRIX,rmatrix%RmatrixBlock(2,2),&
                rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil,&
                InodeList=p_IedgesDirichletBC)
            end if
          end if

        case default
          ! No stabilisation
        
        end select
      
      end if ! gamma <> 0
      
    end subroutine
      
    ! -----------------------------------------------------
    
    subroutine assembleGradientMatrices (rnonlinearCCMatrix,rmatrix,bsharedMatrix,rproblem)
    
    ! Initialises the gradient/divergence matrices with entries from
    ! the rnonlinearCCMatrix structure.
    !
    ! The routine copies references from the submatrices tormatrix,
    ! but it does not initialise any matrix weights / scaling factors.
    !
    ! If bsharedMatrix=TRUE, the matrix is created using references to the
    ! matrix building blocks in rlevelInfo, thus sharing all information
    ! with those matrices in rnonlinearCCMatrix. In this case, the caller must
    ! not change the matrix entries, because this would change the
    ! original 'template' matrices!
    ! (This can be used e.g. for setting up a matrix for building a defect
    !  vector without copying matrix data.)
    ! If bsharedMatrix=TRUE on the other hand, the matrix entries of the
    ! original template (B-) matrices are copied in memory,
    ! so the new matrix is allowed to be changed!

    ! A t_nonlinearCCMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix.
    type(t_nonlinearCCMatrix), intent(in) :: rnonlinearCCMatrix

    ! Block matrix where the B-matrices should be set up
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! Whether or not the matrix entries of the source gradient-matrices
    ! should be copied in memory.
    !
    ! If set to FALSE, the routine tries to initialise rmatrix
    ! only with references to the original matrices, thus the caller must not
    ! change the entries. Nevertheless, if rmatrix is the owner of one of the
    ! submatrices, the routine will always copy the matrix entries,
    ! as otherwise memory would have to be deallocated!
    !
    ! If set to TRUE, the entries of the source matrices in rnonlinearCCMatrix are
    ! copied, so the caller can change rmatrix afterwards (e.g. to implement
    ! boundary conditions).
    logical, intent(in) :: bsharedMatrix

    ! Standard problem structure that defines all underlying parameters.
    type(t_problem), intent(inout), target :: rproblem

      ! local variables
      integer :: idubStructure,idubContent
      
      ! Initialise a copy flag that tells the duplicateMatrix-routine whether to
      ! copy the entries or to create references.
      if (bsharedMatrix) then
      
        idubContent = LSYSSC_DUP_SHARE
        
        ! Normally we share entries -- except for if the submatrices belong to
        ! rmatrix! To avoid memory deallocation in this case, we copy
        ! the entries.
        if ((.not. lsyssc_isMatrixContentShared(Rmatrix%RmatrixBlock(1,3))) .or.&
            (.not. lsyssc_isMatrixContentShared(Rmatrix%RmatrixBlock(2,3)))) then
          idubContent = LSYSSC_DUP_COPY
        end if
        
      else
      
        idubContent = LSYSSC_DUP_COPY
        
      end if
      
      idubStructure = LSYSSC_DUP_SHARE
      
      ! Let us consider the global system in detail:
      !
      !    ( A11  A12  B1  ) = ( A11  A12  A13 )
      !    ( A21  A22  B2  )   ( A21  A22  A23 )
      !    ( D1   D2   0   )   ( A31  A32  A33 )
      !
      ! We exclude the velocity submatrices here, so our system looks like:
      !
      !    (           B1 ) = (           A13 )
      !    (           B2 )   (           A23 )
      !    ( D1   D2      )   ( A31  A32  A33 )

      ! The B1/B2 matrices exist up to now only in rnonlinearCCMatrix.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2 with those B1/B2 of the
      ! block matrix, while we create copies of the entries. The B-blocks
      ! are already prepared and memory for the entries is already allocated;
      ! so we only have to copy the entries.
      !
      ! Note that idubContent = LSYSSC_DUP_COPY will automatically allocate
      ! memory if necessary.
      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixB1, &
                                    rmatrix%RmatrixBlock(1,3),&
                                    idubStructure,idubContent)

      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixB2, &
                                    rmatrix%RmatrixBlock(2,3),&
                                    idubStructure,idubContent)
      
      ! Now, prepare D1 and D2. These matrices always share
      ! their data with the 'template' matrices as the data in these
      ! matrices is usually not overwritten by boundary conditions.
      ! Check the flag bvirtualTransposedD; this decides on whether Di are
      ! created as virtually transposed B-matrices or by taking the D-matrices.
      if (rnonlinearCCMatrix%bvirtualTransposedD) then
        call lsyssc_transposeMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixD1T,&
            rmatrix%RmatrixBlock(3,1),LSYSSC_TR_VIRTUAL)
            
        call lsyssc_transposeMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixD2T,&
            rmatrix%RmatrixBlock(3,2),LSYSSC_TR_VIRTUAL)
      else
        call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixD1, &
                                      rmatrix%RmatrixBlock(3,1),&
                                      LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

        call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixD2, &
                                      rmatrix%RmatrixBlock(3,2),&
                                      LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      end if
                                    
    end subroutine

    ! -----------------------------------------------------
    
    subroutine assemblePressureMatrix (rnonlinearCCMatrix,rmatrix,rproblem)
    
    ! Initialises the pressure matrix with entries from
    ! the rnonlinearCCMatrix structure.

    ! A t_nonlinearCCMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix.
    type(t_nonlinearCCMatrix), intent(in) :: rnonlinearCCMatrix

    ! Block matrix where the C-matrix (3,3) should be set up
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! Standard problem structure that defines all underlying parameters.
    type(t_problem), intent(inout), target :: rproblem
      !local variables:
      type(t_bilinearForm) :: rform
      type(t_spatialDiscretisation) :: rdiscretisation
      type(t_collection) :: rcollection
      
      real(DP) :: dweight1
      
      ! Get the weight of the mass matrix from the physics-structure
      dweight1 = rproblem%rphysics%dmuWeight
      
      ! For the moment, there cannot be found much in P.
      ! If the matrix exists (scale factor <> 0), we clear the
      ! content, otherwise we ignore it.
      if (rmatrix%RmatrixBlock(3,3)%dscaleFactor .ne. 0.0_DP) then
      
        ! rnonlinearCCMatrix%dmu defines the weight for the timestepping!

        ! TO BE MODIFIED!

!        call lsyssc_matrixLinearComb (&
!              rnonlinearCCMatrix%p_rasmTempl%rmatrixMassPressure,&
!              rnonlinearCCMatrix%dmu*dweight1*rmatrix%RmatrixBlock(3,3)%dscaleFactor,&
!              rmatrix%RmatrixBlock(3,3),1.0_DP,&
!              rmatrix%RmatrixBlock(3,3),&
!              .false.,.false.,.true.,.true.)          
        rform%itermCount = 1
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_FUNC
        
        rform%ballCoeffConstant = .false.
        rform%BconstantCoeff = .false.
        rform%Dcoefficients(1)  = 1.0
        
        ! Matrix weight
        rcollection%DquickAccess(1) = rnonlinearCCMatrix%dmu*dweight1
        
        call bilf_buildMatrixScalar (rform, .true., rmatrix%RmatrixBlock(3,3),&
                coeff_Domain_features,rcollection)  
        
      else
        
        call output_line ("No mass matrix?!?", &
            OU_CLASS_ERROR,OU_MODE_STD,'assemblePressureMatrix')
        call sys_halt()
        
      end if
                                    
    end subroutine

    ! -----------------------------------------------------
    
    subroutine assembleKappaMatrix (rnonlinearCCMatrix,rmatrix,rproblem,rvector)
    
    ! Initialises the pressure matrix with entries from
    ! the rnonlinearCCMatrix structure.

    ! A t_nonlinearCCMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix.
    type(t_nonlinearCCMatrix), intent(in) :: rnonlinearCCMatrix

    ! Block matrix where the C-matrix (3,3) should be set up
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! Standard problem structure that defines all underlying parameters.
    type(t_problem), intent(inout), target :: rproblem

    ! Solution vector
    type(t_vectorBlock), intent(in), optional :: rvector
    
    !local variables:
      type(t_bilinearForm) :: rform
      type(t_spatialDiscretisation) :: rdiscretisation
      type(t_collection) :: rcollection

      real(DP) :: dweight2
      
      ! Get the weight of the mass matrix from the physics-structure
      dweight2 = rproblem%rphysics%dkappaWeight

      ! For the moment, there cannot be found much in K.
      ! If the matrix exists (scale factor <> 0), we clear the
      ! content, otherwise we ignore it.
      if (rmatrix%RmatrixBlock(3,4)%dscaleFactor .ne. 0.0_DP) then
      
        call lsyssc_clearMatrix (rmatrix%RmatrixBlock(3,4))

        ! rnonlinearCCMatrix%dkappa defines the weight for the timestepping!
      
        ! TO BE MODIFIED!
!        call stdop_assembleSimpleMatrix (rmatrix%RmatrixBlock(3,4),&
!            DER_FUNC,DER_FUNC,&
!            rnonlinearCCMatrix%dkappa*dweight2*rmatrix%RmatrixBlock(3,4)%dscaleFactor,&
!            bclear=.false.)
        
        rform%itermCount = 1
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_FUNC
        
        rform%ballCoeffConstant = .false.
        rform%BconstantCoeff = .false.
        rform%Dcoefficients(1)  = 1.0
        
        rcollection%DquickAccess(1) = &
                rnonlinearCCMatrix%dkappa*dweight2*rmatrix%RmatrixBlock(3,4)%dscaleFactor
        call bilf_buildMatrixScalar (rform, .true., rmatrix%RmatrixBlock(3,4),&
                coeff_Domain_features,rcollection)     
        
      else
        
        call output_line ("No mass matrix?!?", &
            OU_CLASS_ERROR,OU_MODE_STD,'assemblePressureMatrix')
        call sys_halt()
        
      end if
                                    
    end subroutine

    ! -----------------------------------------------------
    
    subroutine assembleConcentrationMatrix (rnonlinearCCMatrix,rmatrix,rproblem,rvector)
    
    ! Initialises the pressure matrix with entries from
    ! the rnonlinearCCMatrix structure.

    ! A t_nonlinearCCMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix.
    type(t_nonlinearCCMatrix), intent(in) :: rnonlinearCCMatrix

    ! Block matrix where the C-matrix (3,3) should be set up
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! Standard problem structure that defines all underlying parameters.
    type(t_problem), intent(inout), target :: rproblem

    ! Solution vector
    type(t_vectorBlock), intent(in), target, optional :: rvector
    
      ! local variables
      type(t_convStreamDiff2) :: rstreamlineDiffusion2
      type(t_matrixBlock) :: rmatrixTempBlock
      type(t_collection) :: rcollection
      type(t_bilinearForm) :: rform

      ! For the moment, there cannot be found much in C.
      ! If the matrix exists (scale factor <> 0), we clear the
      ! content, otherwise we ignore it.
      if (rmatrix%RmatrixBlock(4,4)%dscaleFactor .ne. 0.0_DP) then
      
        call lsyssc_clearMatrix (rmatrix%RmatrixBlock(4,4))

        ! Mass for time-dependent problems!
        call lsyssc_matrixLinearComb (&
              rnonlinearCCMatrix%p_rasmTempl%rmatrixMassConcentration,&
              rnonlinearCCMatrix%dalpha, &
              !+ rnonlinearCCMatrix%djota*rproblem%rphysics%dreactionWeight,&
              rmatrix%RmatrixBlock(4,4),1.0_DP,&
              rmatrix%RmatrixBlock(4,4),&
              .false.,.false.,.true.,.true.)   
        !For the position dependend domain
        rform%itermCount = 1
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_FUNC
        
        rform%ballCoeffConstant = .false.
        rform%BconstantCoeff = .false.
        rform%Dcoefficients(1)  = 1.0
        
        rcollection%DquickAccess(1) = &
                rnonlinearCCMatrix%djota*rproblem%rphysics%dreactionWeight
                
        call bilf_buildMatrixScalar (rform, .FALSE., rmatrix%RmatrixBlock(4,4),&
                coeff_Domain_features,rcollection)  
                             
        ! Assemble convection+Laplace
        call lsyssc_matrixLinearComb (&
              rnonlinearCCMatrix%p_rasmTempl%rmatrixLaplaceConcentration,&
              rnonlinearCCMatrix%djota*rproblem%rphysics%ddiffusionWeight,&
              rmatrix%RmatrixBlock(4,4),1.0_DP,&
              rmatrix%RmatrixBlock(4,4),&
              .false.,.false.,.true.,.true.)

        ! Fixed convection
        call stdop_assembleSimpleMatrix (rmatrix%RmatrixBlock(4,4),&
            DER_DERIV_X,DER_FUNC,&
            rnonlinearCCMatrix%djota*rproblem%rphysics%dconvectionBeta1,&
            bclear=.false.)

        call stdop_assembleSimpleMatrix (rmatrix%RmatrixBlock(4,4),&
            DER_DERIV_Y,DER_FUNC,&
            rnonlinearCCMatrix%djota*rproblem%rphysics%dconvectionBeta2,&
            bclear=.false.)

        ! ---------------------------------------------------
        ! That was easy -- the adventure begins now... The nonlinearity!
        if (rproblem%rphysics%dconvectionWeight .ne. 0.0_DP) then
        
          if (.not. present(rvector)) then
            call output_line ('Velocity vector not present!', &
                              OU_CLASS_ERROR,OU_MODE_STD,'assembleConcentrationMatrix')
            stop
          end if
        
          ! ---------------------------------------------------
          ! That was easy -- the adventure begins now... The nonlinearity!
          if (rproblem%rphysics%dconvectionWeight .ne. 0.0_DP) then
          
            rform%itermCount = 2
            rform%Idescriptors(1,1) = DER_DERIV_X
            rform%Idescriptors(2,1) = DER_FUNC
            rform%Idescriptors(1,2) = DER_DERIV_Y
            rform%Idescriptors(2,2) = DER_FUNC

            rform%ballCoeffConstant = .false.
            rform%BconstantCoeff = .false.
            rform%Dcoefficients(1)  = 1.0
            rform%Dcoefficients(2)  = 1.0
            
            rcollection%DquickAccess(2) = &
                rnonlinearCCMatrix%djota * rproblem%rphysics%dconvectionWeight
            rcollection%p_rvectorQuickAccess1 => rvector
            
            call bilf_buildMatrixScalar2 (rform, .false., rmatrix%RmatrixBlock(4,4),&
                fcoeff_convection,rcollection)

          end if
          
        end if
        
      end if
                                    
    end subroutine


  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_nonlinearMatMul (rnonlinearCCMatrix,rx,rd,dcx,dcd,rproblem,ry)

!<description>
  ! This routine performs a matrix-vector multiplication with a nonlinear
  ! matrix:
  !      rd := cx A(ry) rx + cd rd
  ! with the system matrix A(.) defined by the configuration in rnonlinearCCMatrix.
  ! The caller must initialise the rnonlinearCCMatrix according to how the
  ! matrix should look like.
  !
  ! The parameter ry is optional. If specified, this parameter defines where to
  ! evaluate the nonlinearity (if the system matrix $A$ contains a nonlinearity).
  ! If not specified, ry=rx is assumed.
  !
  ! The routine will not include any boundary conditions in the defect.
!</description>

  ! A t_nonlinearCCMatrix structure providing all necessary 'source' information
  ! about how to set up the matrix.
  !
  ! The caller must provide either p_rmatrixTemplateXXXX in this structure
  ! or set the p_rmatrixTemplateXXXX as well as p_rdiscretisation to
  ! appropriate values. This is necessary for exploiting then structure
  ! of the matrix.
  type(t_nonlinearCCMatrix), intent(in) :: rnonlinearCCMatrix

  ! This vector specifies the 'x' that is multiplied to the matrix.
  type(t_vectorBlock), intent(in), target :: rx

  ! Multiplication factor in front of the term 'A(ry) rx'.
  real(DP), intent(in) :: dcx

  ! Multiplication factor in front of the term 'rd'.
  real(DP), intent(in) :: dcd

  ! Standard problem structure that defines all underlying parameters.
  type(t_problem), intent(inout), target :: rproblem

  ! OPTIONAL: Point where to evaluate the nonlinearity. If not specified,
  ! ry=rx is assumed.
  type(t_vectorBlock), intent(in), target, optional :: ry

!</input>

!<inputoutput>
  ! Destination vector. cx*A(ry)*rx is subtracted from this vector.
  type(t_vectorBlock), intent(inout) :: rd
!</inputoutput>
  
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_ry
    type(t_matrixBlock) :: rmatrix
    
    ! DEBUG!!!
    real(dp), dimension(:), pointer :: p_DdataX,p_DdataD
    
    call lsysbl_getbase_double (rx,p_DdataX)
    call lsysbl_getbase_double (rd,p_DdataD)
    
    p_ry => rx
    if (present(ry)) p_ry => ry

    ! Probably weight the input vector.
    if (dcd .ne. 1.0_DP) then
      call lsysbl_scaleVector (rd,dcd)
    end if
    
    ! The system matrix looks like:
    !
    !    ( A11  A12  B1  .  )
    !    ( A21  A22  B2  .  )
    !    ( D1   D2   P   K  )
    !    ( .    .    .   C  )
    !
    ! Create a temporary matrix that covers this structure.
    call lsysbl_createMatBlockByDiscr (rnonlinearCCMatrix%p_rdiscretisation,rmatrix)
    
    ! Put references to the Stokes- and B-matrices to Aij. assembleVelocityDefect
    ! needs this template matrix to provide the structure for the stabilisation
    ! routines! The B-matrices are needed later.
    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixStokes,&
        rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixStokes,&
        rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    
    if (rnonlinearCCMatrix%dnewton .ne. 0.0_DP) then
      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixStokes,&
          rmatrix%RmatrixBlock(1,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixStokes,&
          rmatrix%RmatrixBlock(2,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    end if
    
    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixB1,&
        rmatrix%RmatrixBlock(1,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixB2,&
        rmatrix%RmatrixBlock(2,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixD1,&
        rmatrix%RmatrixBlock(3,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixD2,&
        rmatrix%RmatrixBlock(3,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixTemplateFEMPressure,&
        rmatrix%RmatrixBlock(3,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixTemplateK,&
        rmatrix%RmatrixBlock(3,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    call lsyssc_duplicateMatrix (rnonlinearCCMatrix%p_rasmTempl%rmatrixTemplateFEMConcentration,&
        rmatrix%RmatrixBlock(4,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

    ! In the first step, we assemble the defect that arises in the velocity
    ! components. This is characterised by the following submatrix:
    !
    !    ( A11  A12  .    . )
    !    ( A21  A22  .    . )
    !    ( .    .    .    . )
    !    ( .    .    .    . )
    !
    ! assembleVelocityDefect handles exactly these submatrices.

    call assembleVelocityDefect (rnonlinearCCMatrix,rmatrix,rx,rd,p_ry,-dcx,rproblem)
    
    ! Go on with the pressure defect:
    !
    !    ( .    .    .    . )
    !    ( .    .    .    . )
    !    ( .    .    P    K )
    !    ( .    .    .    . )
    !
    ! assembleVelocityDefect handles exactly these submatrices.

    call assemblePressureDefect (rnonlinearCCMatrix,rmatrix,rx,rd,p_ry,-dcx,rproblem)
    
    ! And the concentration defect
    !
    !    ( .    .    .    . )
    !    ( .    .    .    . )
    !    ( .    .    .    . )
    !    ( .    .    .    C )
    !
    ! assembleVelocityDefect handles exactly these submatrices.

    call assembleConcentrationDefect (rnonlinearCCMatrix,rmatrix,rx,rd,p_ry,-dcx,rproblem)
    
    ! Now, we treat all the remaining blocks. Let us see what is missing:
    !
    !    ( .    .    B1   . )
    !    ( .    .    B2   . )
    !    ( D1   D2   .    . )
    !    ( .    .    .    . )

    ! To build the appropriate defect, we first remove the velocity blocks:
    
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(1,1))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(1,2))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(2,1))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(2,2))
    
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(3,3))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(3,4))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(4,4))

    ! Initialise the weights for the B/B^T matrices
    rmatrix%RmatrixBlock(1,3)%dscaleFactor = rnonlinearCCMatrix%deta
    rmatrix%RmatrixBlock(2,3)%dscaleFactor = rnonlinearCCMatrix%deta
    
    rmatrix%RmatrixBlock(3,1)%dscaleFactor = rnonlinearCCMatrix%dtau
    rmatrix%RmatrixBlock(3,2)%dscaleFactor = rnonlinearCCMatrix%dtau

    ! ------------------------------------------------
    ! Build the defect by matrix-vector multiplication
    !
    ! Note that no time step or whatever is included here; everything
    ! is initialised with the multiplication factors in the submatrices
    ! from above!
    call lsysbl_blockMatVec (rmatrix, rx, rd, dcx, 1.0_DP)
    
    ! Release the temporary matrix, we do not need it anymore.
    call lsysbl_releaseMatrix (rmatrix)

  contains

    subroutine assembleVelocityDefect (rnonlinearCCMatrix,&
        rmatrix,rvector,rdefect,rvelocityVector,dvectorWeight,rproblem)
        
    ! Assembles the velocity defect in the block matrix rmatrix at position
    ! itop..itop+1 in the velocity vector. rdefect must have been initialised
    ! with the right hand side vector.
    !
    ! With a matrix 'A' of the theoretical form
    !
    !       A := dalpha*M + dtheta*Laplace + dgamma*N(p_rvector) +
    !            dnewton*N*(p_rvector)
    !
    ! and c=dvectorWeight, the routine will construct
    !
    !       rdefect = rdefect - c A rvector
    
    ! A t_nonlinearCCMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix.
    type(t_nonlinearCCMatrix), intent(in) :: rnonlinearCCMatrix

    ! Reference to the system matrix. Only the structure of the matrix
    ! is used to reconstruct the structure of the discretisation.
    ! The content of the matrix is not changed or used.
    type(t_matrixBlock), intent(inout) :: rmatrix
    
    ! Solution vector.
    type(t_vectorBlock), intent(in) :: rvector
    
    ! On entry: RHS vector.
    ! Is overwritten by the defect vector in the velocity subsystem.
    type(t_vectorBlock), intent(inout) :: rdefect
    
    ! Weight for the velocity vector; usually = 1.0
    real(DP), intent(in) :: dvectorWeight
    
    ! Velocity vector field that should be used for the assembly of the
    ! nonlinearity. The first two blocks in that block vector are
    ! used as velocity field.
    type(t_vectorBlock), intent(in), target :: rvelocityVector

    ! Standard problem structure that defines all underlying parameters.
    type(t_problem), intent(inout), target :: rproblem

    ! local variables
    logical :: bshared
    type(t_convUpwind) :: rupwind
    type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
    type(t_convStreamDiff2) :: rstreamlineDiffusion2
    type(T_jumpStabilisation) :: rjumpStabil
    type(t_collection) :: rcollection
    integer, dimension(:), pointer :: p_IedgesDirichletBC
    
    ! DEBUG!!!
    real(dp), dimension(:), pointer :: p_DdataX,p_DdataD
    
    call lsysbl_getbase_double (rvector,p_DdataX)
    call lsysbl_getbase_double (rdefect,p_DdataD)

      ! Is A11=A22 physically?
      bshared = lsyssc_isMatrixContentShared(&
                    rmatrix%RmatrixBlock(1,1),&
                    rmatrix%RmatrixBlock(2,2)) .or.&
                (.not. lsyssc_hasMatrixContent(rmatrix%RmatrixBlock(1,1)) .and.&
                 .not. lsyssc_hasMatrixContent(rmatrix%RmatrixBlock(2,2)))

      ! ---------------------------------------------------
      ! Subtract the mass matrix stuff?
      if (rnonlinearCCMatrix%dalpha .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (rnonlinearCCMatrix%p_rasmTempl%rmatrixMass, &
            rvector%RvectorBlock(1), rdefect%RvectorBlock(1), &
            -rnonlinearCCMatrix%dalpha, 1.0_DP)

        call lsyssc_scalarMatVec (rnonlinearCCMatrix%p_rasmTempl%rmatrixMass, &
            rvector%RvectorBlock(2), rdefect%RvectorBlock(2), &
            -rnonlinearCCMatrix%dalpha, 1.0_DP)
      end if
      
      ! ---------------------------------------------------
      ! Subtract the Stokes matrix stuff?
      if (rnonlinearCCMatrix%dtheta .ne. 0.0_DP) then
      
        ! In case of the gradient tensor, we can directly substract
        ! the Stokes matrix.
        if ((rnonlinearCCMatrix%p_rphysics%isubequation .eq. 0) .and. &
            (rnonlinearCCMatrix%p_rphysics%cviscoModel .eq. 0)) then
      
          call lsyssc_scalarMatVec (rnonlinearCCMatrix%p_rasmTempl%rmatrixStokes, &
              rvector%RvectorBlock(1), rdefect%RvectorBlock(1), &
              -rnonlinearCCMatrix%dtheta, 1.0_DP)

          call lsyssc_scalarMatVec (rnonlinearCCMatrix%p_rasmTempl%rmatrixStokes, &
              rvector%RvectorBlock(2), rdefect%RvectorBlock(2), &
              -rnonlinearCCMatrix%dtheta, 1.0_DP)
              
        end if
        
        ! Otherwise, the operator is substracted during the assembly branch
        ! of the nonlinearity.
        
      end if
      
      ! ---------------------------------------------------
      ! That was easy -- the adventure begins now... The nonlinearity!
      if ((rnonlinearCCMatrix%dgamma .ne. 0.0_DP) .or. &
          (rnonlinearCCMatrix%p_rphysics%isubequation .ne. 0) .or. &
          (rnonlinearCCMatrix%p_rphysics%cviscoModel .ne. 0)) then
      
        ! Type of stablilisation?
        select case (rnonlinearCCMatrix%p_rstabilisation%iupwind)
        case (CCMASM_STAB_STREAMLINEDIFF)
          ! Set up the SD structure for the creation of the defect.
          ! There is not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearCCMatrix%p_rphysics%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion%dupsam = rnonlinearCCMatrix%p_rstabilisation%dupsam
          
          ! Set calculation method for local H
          rstreamlineDiffusion%clocalH = rnonlinearCCMatrix%p_rstabilisation%clocalH
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion%ddelta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newtop part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = rnonlinearCCMatrix%dnewton
          
          ! Call the SD method to calculate the defect of the nonlinearity.
          ! As rrhsTemp shares its entries with rdefect, the result is
          ! directly written to rdefect!
          ! As velocity field, we specify rvelocityVector here. The first two
          ! subvectors are used as velocity field.
          
          call conv_streamlineDiffusionBlk2d (&
                              rvelocityVector, rvelocityVector, &
                              dvectorWeight, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODDEFECT, &
                              rmatrix,rsolution=rvector,rdefect=rdefect)
                              
        case (CCMASM_STAB_STREAMLINEDIFF2)
                  
          ! Set up the SD structure for the creation of the defect.
          ! There is not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearCCMatrix%p_rphysics%dnu
          
          ! Probably, we have nonconstant viscosity.
          ! In that case, init a collection structure for a callback
          ! routine that specifies the viscosity
          if (rnonlinearCCMatrix%p_rphysics%cviscoModel .ne. 0) then
            rstreamlineDiffusion2%bconstnu = .false.
            
            ! Assemble at least the Stokes matrix.
            ! For a more complicated formulation, dbeta/dbetaT
            ! may be changed later.
            rstreamlineDiffusion2%dbeta = 1.0_DP

            ! Prepare the collection. The "next" collection points to the user defined
            ! collection.
            rcollection%p_rnextCollection => rproblem%rcollection
            call ccmva_prepareViscoAssembly (rproblem,rnonlinearCCMatrix%p_rphysics,&
                rcollection,rvelocityVector)

            ! The "next" collection points to the user defined collection.
            rcollection%p_rnextCollection => rproblem%rcollection

          end if
          
          ! Set stabilisation parameter
          rstreamlineDiffusion2%dupsam = rnonlinearCCMatrix%p_rstabilisation%dupsam
          
          ! Set calculation method for local H
          rstreamlineDiffusion2%clocalH = rnonlinearCCMatrix%p_rstabilisation%clocalH
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion2%ddelta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newtop part; =0 deactivates Newton.
          rstreamlineDiffusion2%dnewton = rnonlinearCCMatrix%dnewton
          
          ! Assemble the deformation tensor?
          if (rnonlinearCCMatrix%p_rphysics%isubequation .ne. 0) then
            rstreamlineDiffusion2%dbeta = 0.5_DP
            rstreamlineDiffusion2%dbetaT = 0.5_DP
          end if

          ! Initialise the user defined collection for the assembly.
          call cc_initCollectForAssembly (rproblem,rproblem%rcollection)

          ! Call the SD method to assemble.
          call conv_streamDiff2Blk2dDef (rstreamlineDiffusion2,rmatrix,&
              rvector,rdefect,rvelocityVector,ffunctionViscoModel,rcollection)
                              
          ! That is it.
          call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

        case (CCMASM_STAB_UPWIND)
          ! Set up the upwind structure for the creation of the defect.
          ! There is not much to do, only initialise the viscosity...
          rupwind%dnu = rnonlinearCCMatrix%p_rphysics%dnu
          
          ! Set stabilisation parameter
          rupwind%dupsam = rnonlinearCCMatrix%p_rstabilisation%dupsam

          ! Matrix weight
          rupwind%dtheta = rnonlinearCCMatrix%dgamma
          
          ! Call the upwind method to calculate the nonlinear defect.
          call conv_upwind2d (rvelocityvector, rvelocityvector, &
                              dvectorWeight, 0.0_DP,&
                              rupwind, CONV_MODDEFECT, &
                              rmatrix%RmatrixBlock(1,1),rvector,rdefect)
                              
          if (.not. bshared) then
            call output_line ('Upwind does not support independent A11/A22!', &
                OU_CLASS_ERROR,OU_MODE_STD,'assembleVelocityDefect')
            call sys_halt()
          end if

        case (CCMASM_STAB_EDGEORIENTED)
          ! Jump stabilisation.
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.

          ! There is not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearCCMatrix%p_rphysics%dnu
          
          ! Set stabilisation parameter to zero for central difference
          rstreamlineDiffusion%dupsam = 0.0_DP
          
          ! Set calculation method for local H
          rstreamlineDiffusion%clocalH = rnonlinearCCMatrix%p_rstabilisation%clocalH
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion%ddelta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newtop part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = rnonlinearCCMatrix%dnewton
          
          if (rnonlinearCCMatrix%dnewton .eq. 0.0_DP) then

            ! Deactivate the matrices A12 and A21 by setting the multiplicators
            ! to 0.0. Whatever the content is (if there is content at all),
            ! these matrices are ignored then by the kernel.
            
            rmatrix%RmatrixBlock(1,2)%dscaleFactor = 0.0_DP
            rmatrix%RmatrixBlock(2,1)%dscaleFactor = 0.0_DP
            
          else

            ! Clear A12/A21 that receives parts of the Newton matrix
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,2))
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,1))
          
            ! Activate the submatrices A12 and A21 if they are not.
            rmatrix%RmatrixBlock(1,2)%dscaleFactor = 1.0_DP
            rmatrix%RmatrixBlock(2,1)%dscaleFactor = 1.0_DP
           
          end if
         
          ! Call the SD method to calculate the nonlinearity.
          call conv_streamlineDiffusionBlk2d (&
                              rvelocityVector, rvelocityVector, &
                              dvectorWeight, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODDEFECT, &
                              rmatrix,rsolution=rvector,rdefect=rdefect)
        
          ! Set up the jump stabilisation structure.
          ! There is not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rnonlinearCCMatrix%p_rphysics%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgamma = rnonlinearCCMatrix%p_rstabilisation%dupsam
          rjumpStabil%dgammastar = rnonlinearCCMatrix%p_rstabilisation%dupsamstar
          rjumpStabil%deojEdgeExp = rnonlinearCCMatrix%p_rstabilisation%deojEdgeExp
          
          ! Matrix weight
          rjumpStabil%dtheta = rnonlinearCCMatrix%dtheta

          ! Cubature formula
          rjumpStabil%ccubType = rnonlinearCCMatrix%p_rstabilisation%ccubEOJ

          ! Call the jump stabilisation technique to stabilise that stuff.
          ! We can assemble the jump part any time as it is independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
              rjumpStabil, CONV_MODDEFECT,rmatrix%RmatrixBlock(1,1),&
              rsolution=rvector,rdefect=rdefect,&
              rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil)

          if (.not. bshared) then
            call output_line (&
                'Edge oriented stabilisation does not support independent A11/A22!', &
                OU_CLASS_ERROR,OU_MODE_STD,'assembleVelocityDefect')
            call sys_halt()
          end if
          
          ! Subtract the EOJ matrix for the Dirichlet boundary conditions.
          if (rnonlinearCCMatrix%p_rdynamicInfo%nedgesDirichletBC .ne. 0) then
            rjumpStabil%dnu = rnonlinearCCMatrix%p_rphysics%dnu
            rjumpStabil%dgamma = rnonlinearCCMatrix%p_rstabilisation%dupsam
            rjumpStabil%dtheta = -rnonlinearCCMatrix%dtheta
            rjumpStabil%ccubType = rnonlinearCCMatrix%p_rstabilisation%ccubEOJ
            call storage_getbase_int(rnonlinearCCMatrix%p_rdynamicInfo%hedgesDirichletBC,&
                p_IedgesDirichletBC)
            call conv_jumpStabilisation2d (&
                rjumpStabil, CONV_MODDEFECT,rmatrix%RmatrixBlock(1,1),&
                rsolution=rvector,rdefect=rdefect,&
                rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil,&
                InodeList=p_IedgesDirichletBC)
          end if

        case (CCMASM_STAB_EDGEORIENTED2)
          ! Jump stabilisation.
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! Set up the SD structure for the creation of the defect.
          ! There is not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearCCMatrix%p_rphysics%dnu
          
          ! Probably, we have nonconstant viscosity.
          ! In that case, init a collection structure for a callback
          ! routine that specifies the viscosity
          if (rnonlinearCCMatrix%p_rphysics%cviscoModel .ne. 0) then
            rstreamlineDiffusion2%bconstnu = .false.
            
            ! Assemble at least the Stokes matrix.
            ! For a more complicated formulation, dbeta/dbetaT
            ! may be changed later.
            rstreamlineDiffusion2%dbeta = 1.0_DP

            ! Prepare the collection. The "next" collection points to the user defined
            ! collection.
            rcollection%p_rnextCollection => rproblem%rcollection
            call ccmva_prepareViscoAssembly (rproblem,rnonlinearCCMatrix%p_rphysics,&
                rcollection,rvelocityvector)

          end if
          
          ! Set stabilisation parameter to zero for central difference
          rstreamlineDiffusion2%dupsam = 0.0_DP
          
          ! Set calculation method for local H
          rstreamlineDiffusion2%clocalH = rnonlinearCCMatrix%p_rstabilisation%clocalH
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion2%ddelta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newtop part; =0 deactivates Newton.
          rstreamlineDiffusion2%dnewton = rnonlinearCCMatrix%dnewton
          
          ! Assemble the deformation tensor?
          if (rnonlinearCCMatrix%p_rphysics%isubequation .ne. 0) then
            rstreamlineDiffusion2%dbeta = 0.5_DP
            rstreamlineDiffusion2%dbetaT = 0.5_DP
          end if

          if ((rnonlinearCCMatrix%dnewton .eq. 0.0_DP) .and. &
              (rstreamlineDiffusion2%dbetaT .eq. 0.0_DP)) then

            ! Deactivate the matrices A12 and A21 by setting the multiplicators
            ! to 0.0. Whatever the content is (if there is content at all),
            ! these matrices are ignored then by the kernel.
            
            rmatrix%RmatrixBlock(1,2)%dscaleFactor = 0.0_DP
            rmatrix%RmatrixBlock(2,1)%dscaleFactor = 0.0_DP
            
          else

            ! Clear A12/A21 that receives parts of the Newton matrix
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,2))
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,1))
          
            ! Activate the submatrices A12 and A21 if they are not.
            rmatrix%RmatrixBlock(1,2)%dscaleFactor = 1.0_DP
            rmatrix%RmatrixBlock(2,1)%dscaleFactor = 1.0_DP
           
          end if
         
          ! Initialise the user defined collection for the assembly.
          call cc_initCollectForAssembly (rproblem,rproblem%rcollection)

          ! Call the SD method to assemble.
          call conv_streamDiff2Blk2dDef (rstreamlineDiffusion2,rmatrix,&
              rvector,rdefect,rvelocityVector,ffunctionViscoModel,rcollection)
                              
          ! Set up the jump stabilisation structure.
          ! There is not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rnonlinearCCMatrix%p_rphysics%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgamma = rnonlinearCCMatrix%p_rstabilisation%dupsam
          rjumpStabil%dgammastar = rnonlinearCCMatrix%p_rstabilisation%dupsamstar
          rjumpStabil%deojEdgeExp = rnonlinearCCMatrix%p_rstabilisation%deojEdgeExp
          
          ! Matrix weight
          rjumpStabil%dtheta = rnonlinearCCMatrix%dtheta

          ! Cubature formula
          rjumpStabil%ccubType = rnonlinearCCMatrix%p_rstabilisation%ccubEOJ

          ! Call the jump stabilisation technique to stabilise that stuff.
          ! We can assemble the jump part any time as it is independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
              rjumpStabil, CONV_MODDEFECT,rmatrix%RmatrixBlock(1,1),&
              rsolution=rvector,rdefect=rdefect,&
              rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil)

          if (.not. bshared) then
            call output_line (&
                'Edge oriented stabilisation does not support independent A11/A22!', &
                OU_CLASS_ERROR,OU_MODE_STD,'assembleVelocityDefect')
            call sys_halt()
          end if

          ! Subtract the EOJ matrix for the Dirichlet boundary conditions.
          if (rnonlinearCCMatrix%p_rdynamicInfo%nedgesDirichletBC .ne. 0) then
            rjumpStabil%dnu = rnonlinearCCMatrix%p_rphysics%dnu
            rjumpStabil%dgamma = rnonlinearCCMatrix%p_rstabilisation%dupsam
            rjumpStabil%dgammastar = rnonlinearCCMatrix%p_rstabilisation%dupsamstar
            rjumpStabil%deojEdgeExp = rnonlinearCCMatrix%p_rstabilisation%deojEdgeExp
            rjumpStabil%dtheta = -rnonlinearCCMatrix%dtheta
            rjumpStabil%ccubType = rnonlinearCCMatrix%p_rstabilisation%ccubEOJ
            call storage_getbase_int(rnonlinearCCMatrix%p_rdynamicInfo%hedgesDirichletBC,&
                p_IedgesDirichletBC)
            call conv_jumpStabilisation2d (&
                rjumpStabil, CONV_MODDEFECT,rmatrix%RmatrixBlock(1,1),&
                rsolution=rvector,rdefect=rdefect,&
                rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil,&
                InodeList=p_IedgesDirichletBC)
          end if

          ! That is it.
          call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
          
        case (CCMASM_STAB_FASTEDGEORIENTED)
          ! Fast Jump stabilisation. Precomputed matrix.
          
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.

          ! There is not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearCCMatrix%p_rphysics%dnu
          
          ! Set stabilisation parameter to zero for central difference
          rstreamlineDiffusion%dupsam = 0.0_DP
          
          ! Set calculation method for local H
          rstreamlineDiffusion%clocalH = rnonlinearCCMatrix%p_rstabilisation%clocalH
          
          ! Matrix weight for the nonlinearity
          rstreamlineDiffusion%ddelta = rnonlinearCCMatrix%dgamma
          
          ! Weight for the Newton part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = rnonlinearCCMatrix%dnewton
          
          if (rnonlinearCCMatrix%dnewton .eq. 0.0_DP) then

            ! Deactivate the matrices A12 and A21 by setting the multiplicators
            ! to 0.0. Whatever the content is (if there is content at all),
            ! these matrices are ignored then by the kernel.
            
            rmatrix%RmatrixBlock(1,2)%dscaleFactor = 0.0_DP
            rmatrix%RmatrixBlock(2,1)%dscaleFactor = 0.0_DP
            
          else

            ! Clear A12/A21 that receives parts of the Newton matrix
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,2))
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,1))
          
            ! Activate the submatrices A12 and A21 if they are not.
            rmatrix%RmatrixBlock(1,2)%dscaleFactor = 1.0_DP
            rmatrix%RmatrixBlock(2,1)%dscaleFactor = 1.0_DP
           
          end if
         
          ! Call the SD method to calculate the nonlinearity.
          call conv_streamlineDiffusionBlk2d (&
                              rvelocityvector, rvelocityvector, &
                              dvectorWeight, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODDEFECT, &
                              rmatrix,rsolution=rvector,rdefect=rdefect)
        
          ! Subtract the stabilisation matrix stuff.
          call lsyssc_scalarMatVec (rnonlinearCCMatrix%p_rasmTempl%rmatrixStabil, &
              rvector%RvectorBlock(1), rdefect%RvectorBlock(1), &
              -rnonlinearCCMatrix%dtheta, 1.0_DP)

          call lsyssc_scalarMatVec (rnonlinearCCMatrix%p_rasmTempl%rmatrixStabil, &
              rvector%RvectorBlock(2), rdefect%RvectorBlock(2), &
              -rnonlinearCCMatrix%dtheta, 1.0_DP)

          ! Subtract the EOJ matrix for the Dirichlet boundary conditions.
          if (rnonlinearCCMatrix%p_rdynamicInfo%nedgesDirichletBC .ne. 0) then
            rjumpStabil%dnu = rnonlinearCCMatrix%p_rphysics%dnu
            rjumpStabil%dgamma = rnonlinearCCMatrix%p_rstabilisation%dupsam
            rjumpStabil%dgammastar = rnonlinearCCMatrix%p_rstabilisation%dupsamstar
            rjumpStabil%deojEdgeExp = rnonlinearCCMatrix%p_rstabilisation%deojEdgeExp
            rjumpStabil%dtheta = -rnonlinearCCMatrix%dtheta
            rjumpStabil%ccubType = rnonlinearCCMatrix%p_rstabilisation%ccubEOJ
            call storage_getbase_int(rnonlinearCCMatrix%p_rdynamicInfo%hedgesDirichletBC,&
                p_IedgesDirichletBC)
            call conv_jumpStabilisation2d (&
                rjumpStabil, CONV_MODDEFECT,rmatrix%RmatrixBlock(1,1),&
                rsolution=rvector,rdefect=rdefect,&
                rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil,&
                InodeList=p_IedgesDirichletBC)
          end if

        case default
          call output_line ('Don''t know how to set up nonlinearity!?!', &
              OU_CLASS_ERROR,OU_MODE_STD,'assembleVelocityDefect')
          call sys_halt()
        
        end select
      
      else
      
        ! That is the Stokes-case. Jump stabilisation is possible...
        !
        ! Type of stablilisation?
        select case (rnonlinearCCMatrix%p_rstabilisation%iupwind)
        case (CCMASM_STAB_EDGEORIENTED)
          ! Jump stabilisation.
        
          ! Set up the jump stabilisation structure.
          ! There is not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rnonlinearCCMatrix%p_rphysics%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgamma = rnonlinearCCMatrix%p_rstabilisation%dupsam
          rjumpStabil%dgammastar = rnonlinearCCMatrix%p_rstabilisation%dupsamstar
          rjumpStabil%deojEdgeExp = rnonlinearCCMatrix%p_rstabilisation%deojEdgeExp
          
          ! Matrix weight
          rjumpStabil%dtheta = rnonlinearCCMatrix%dtheta

          ! Cubature formula
          rjumpStabil%ccubType = rnonlinearCCMatrix%p_rstabilisation%ccubEOJ

          ! Call the jump stabilisation technique to stabilise that stuff.
          ! We can assemble the jump part any time as it is independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
              rjumpStabil, CONV_MODDEFECT,rmatrix%RmatrixBlock(1,1),&
              rsolution=rvector,rdefect=rdefect,&
              rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil)

          ! Subtract the EOJ matrix for the Dirichlet boundary conditions.
          if (rnonlinearCCMatrix%p_rdynamicInfo%nedgesDirichletBC .ne. 0) then
            rjumpStabil%dnu = rnonlinearCCMatrix%p_rphysics%dnu
            rjumpStabil%dgamma = rnonlinearCCMatrix%p_rstabilisation%dupsam
            rjumpStabil%dgammastar = rnonlinearCCMatrix%p_rstabilisation%dupsamstar
            rjumpStabil%deojEdgeExp = rnonlinearCCMatrix%p_rstabilisation%deojEdgeExp
            rjumpStabil%dtheta = -rnonlinearCCMatrix%dtheta
            rjumpStabil%ccubType = rnonlinearCCMatrix%p_rstabilisation%ccubEOJ
            call storage_getbase_int(rnonlinearCCMatrix%p_rdynamicInfo%hedgesDirichletBC,&
                p_IedgesDirichletBC)
            call conv_jumpStabilisation2d (&
                rjumpStabil, CONV_MODDEFECT,rmatrix%RmatrixBlock(1,1),&
                rsolution=rvector,rdefect=rdefect,&
                rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil,&
                InodeList=p_IedgesDirichletBC)
          end if

        case (CCMASM_STAB_FASTEDGEORIENTED)
          ! Fast Jump stabilisation. Precomputed matrix.
          
          ! Subtract the stabilisation matrix stuff.
          call lsyssc_scalarMatVec (rnonlinearCCMatrix%p_rasmTempl%rmatrixStabil, &
              rvector%RvectorBlock(1), rdefect%RvectorBlock(1), &
              -rnonlinearCCMatrix%dtheta, 1.0_DP)

          call lsyssc_scalarMatVec (rnonlinearCCMatrix%p_rasmTempl%rmatrixStabil, &
              rvector%RvectorBlock(2), rdefect%RvectorBlock(2), &
              -rnonlinearCCMatrix%dtheta, 1.0_DP)

          ! Subtract the EOJ matrix for the Dirichlet boundary conditions.
          if (rnonlinearCCMatrix%p_rdynamicInfo%nedgesDirichletBC .ne. 0) then
            rjumpStabil%dnu = rnonlinearCCMatrix%p_rphysics%dnu
            rjumpStabil%dgamma = rnonlinearCCMatrix%p_rstabilisation%dupsam
            rjumpStabil%dgammastar = rnonlinearCCMatrix%p_rstabilisation%dupsamstar
            rjumpStabil%deojEdgeExp = rnonlinearCCMatrix%p_rstabilisation%deojEdgeExp
            rjumpStabil%dtheta = -rnonlinearCCMatrix%dtheta
            rjumpStabil%ccubType = rnonlinearCCMatrix%p_rstabilisation%ccubEOJ
            call storage_getbase_int(rnonlinearCCMatrix%p_rdynamicInfo%hedgesDirichletBC,&
                p_IedgesDirichletBC)
            call conv_jumpStabilisation2d (&
                rjumpStabil, CONV_MODDEFECT,rmatrix%RmatrixBlock(1,1),&
                rsolution=rvector,rdefect=rdefect,&
                rdiscretisation=rnonlinearCCMatrix%p_rasmTempl%rdiscretisationStabil,&
                InodeList=p_IedgesDirichletBC)
          end if

        case default
          ! No stabilisation
        
        end select
      
      end if ! gamma <> 0
    
    end subroutine
    
    ! -------------------------------------------------------------------------

    subroutine assemblePressureDefect (rnonlinearCCMatrix,&
        rmatrix,rvector,rdefect,rvelocityVector,dvectorWeight,rproblem)
        
    ! Assembles the pressure defect in the block matrix rmatrix.
    ! rdefect must have been initialised with the right hand side vector.
    !
    ! With a matrix 'A' of the theoretical form
    !
    !       A := mu*P + kappa*K
    !
    ! and c=dvectorWeight, the routine will construct
    !
    !       rdefect = rdefect - c A rvector
    
    ! A t_nonlinearCCMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix.
    type(t_nonlinearCCMatrix), intent(in) :: rnonlinearCCMatrix

    ! Reference to the system matrix. Only the structure of the matrix
    ! is used to reconstruct the structure of the discretisation.
    ! The content of the matrix is not changed or used.
    type(t_matrixBlock), intent(inout) :: rmatrix
    
    ! Solution vector.
    type(t_vectorBlock), intent(in) :: rvector
    
    ! On entry: RHS vector.
    ! Is overwritten by the defect vector in the velocity subsystem.
    type(t_vectorBlock), intent(inout) :: rdefect
    
    ! Weight for the velocity vector; usually = 1.0
    real(DP), intent(in) :: dvectorWeight
     
    ! Velocity vector field that should be used for the assembly of the
    ! nonlinearity. The first two blocks in that block vector are
    ! used as velocity field.
    type(t_vectorBlock), intent(in), target :: rvelocityVector

    ! Standard problem structure that defines all underlying parameters.
    type(t_problem), intent(inout), target :: rproblem

      type(t_matrixScalar) :: rmatrixTemp
      type(t_bilinearForm) :: rform
      type(t_collection) :: rcollection

      real(DP) :: dweight1, dweight2
      
      dweight1 = rproblem%rphysics%dmuWeight
      dweight2 = rproblem%rphysics%dkappaWeight

      ! TO BE MODIFIED!

      ! Generate temp matrix for pressure part
      call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(3,3),&
          rmatrixTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          
      ! Pressure mass matrix
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_FUNC

      rform%ballCoeffConstant = .false.
      rform%BconstantCoeff = .false.
      rform%Dcoefficients(1)  = 1.0
      
      rcollection%DquickAccess(1) = 1.0_DP

      call bilf_buildMatrixScalar (rform, .true., rmatrixTemp,&
              coeff_Domain_features,rcollection)  
      
      call lsyssc_scalarMatVec (rmatrixTemp, &
          rvector%RvectorBlock(3), rdefect%RvectorBlock(3), &
          -rnonlinearCCMatrix%dmu*dweight1, 1.0_DP)
          
      call lsyssc_releaseMatrix (rmatrixTemp)

      ! TO BE MODIFIED!
  
      ! Generate temp matrix for K-part
      call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(3,4),&
          rmatrixTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          
      call lsyssc_clearMatrix (rmatrixTemp)
          
      ! K-matrix
      
       ! Kappa mass matrix
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_FUNC

      rform%ballCoeffConstant = .false.
      rform%BconstantCoeff = .false.
      rform%Dcoefficients(1)  = 1.0

      rcollection%DquickAccess(1) = 1.0_DP

      call bilf_buildMatrixScalar (rform, .true., rmatrixTemp,&
              coeff_Domain_features,rcollection)
!      call stdop_assembleSimpleMatrix (rmatrixTemp,&
!          DER_FUNC,DER_FUNC,rnonlinearCCMatrix%dkappa*dweight2,bclear=.false.)

      call lsyssc_scalarMatVec (rmatrixTemp, &
          rvector%RvectorBlock(4), rdefect%RvectorBlock(3), &
          -rnonlinearCCMatrix%dkappa*dweight2, 1.0_DP)
          
      call lsyssc_releaseMatrix (rmatrixTemp)

    end subroutine

    ! -------------------------------------------------------------------------

    subroutine assembleConcentrationDefect (rnonlinearCCMatrix,&
        rmatrix,rvector,rdefect,rvelocityVector,dvectorWeight,rproblem)
        
    ! Assembles the pressure defect in the block matrix rmatrix.
    ! rdefect must have been initialised with the right hand side vector.
    !
    ! With a matrix 'A' of the theoretical form
    !
    !       A := jota*C
    !
    ! and c=dvectorWeight, the routine will construct
    !
    !       rdefect = rdefect - c A rvector
    
    ! A t_nonlinearCCMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix.
    type(t_nonlinearCCMatrix), intent(in) :: rnonlinearCCMatrix

    ! Reference to the system matrix. Only the structure of the matrix
    ! is used to reconstruct the structure of the discretisation.
    ! The content of the matrix is not changed or used.
    type(t_matrixBlock), intent(inout) :: rmatrix
    
    ! Solution vector.
    type(t_vectorBlock), intent(in) :: rvector
    
    ! On entry: RHS vector.
    ! Is overwritten by the defect vector in the velocity subsystem.
    type(t_vectorBlock), intent(inout) :: rdefect
    
    ! Weight for the velocity vector; usually = 1.0
    real(DP), intent(in) :: dvectorWeight
    
    ! Velocity vector field that should be used for the assembly of the
    ! nonlinearity. The first two blocks in that block vector are
    ! used as velocity field.
    type(t_vectorBlock), intent(in), target :: rvelocityVector

    ! Standard problem structure that defines all underlying parameters.
    type(t_problem), intent(inout), target :: rproblem

      ! local variables
      type(t_matrixScalar) :: rmatrixTemp
      type(t_matrixBlock) :: rmatrixTempBlock
      type(t_vectorBlock) :: rvectorTempBlock,rdefectTempBlock
      type(t_collection) :: rcollection
      type(t_bilinearForm) :: rform
      real(DP), dimension(:), pointer :: p_Ddata

      ! TO BE MODIFIED!
  
      ! Generate Laplace+convection operator
      call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(4,4),&
          rmatrixTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      call lsyssc_getbase_double (rmatrixTemp,p_Ddata)
      
      ! Mass for time-dependent problems!
      call lsyssc_matrixLinearComb (&
            rnonlinearCCMatrix%p_rasmTempl%rmatrixMassConcentration,&
            rnonlinearCCMatrix%dalpha, &
            !+ rnonlinearCCMatrix%djota*rproblem%rphysics%dreactionWeight,&
            rmatrixTemp,0.0_DP,&
            rmatrixTemp,&
            .false.,.false.,.true.,.true.)
      
      !For the position dependend domain
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_FUNC
      
      rform%ballCoeffConstant = .false.
      rform%BconstantCoeff = .false.
      rform%Dcoefficients(1)  = 1.0
      
      rcollection%DquickAccess(1) = rnonlinearCCMatrix%djota*rproblem%rphysics%dreactionWeight
              
      call bilf_buildMatrixScalar (rform, .FALSE., rmatrixTemp,&
              coeff_Domain_features,rcollection)      
              
      ! Generate Laplace+convection operator
      call lsyssc_matrixLinearComb (&
            rnonlinearCCMatrix%p_rasmTempl%rmatrixLaplaceConcentration,&
            rnonlinearCCMatrix%djota*rproblem%rphysics%ddiffusionWeight,&
            rmatrixTemp,1.0_DP,&
            rmatrixTemp,&
            .false.,.false.,.true.,.true.)
      
      ! Constant convection
      call stdop_assembleSimpleMatrix (rmatrixTemp,&
          DER_DERIV_X,DER_FUNC,&
          rnonlinearCCMatrix%djota * rproblem%rphysics%dconvectionBeta1,bclear=.false.)
      call stdop_assembleSimpleMatrix (rmatrixTemp,&
          DER_DERIV_Y,DER_FUNC,&
          rnonlinearCCMatrix%djota * rproblem%rphysics%dconvectionBeta2,bclear=.false.)

      ! ---------------------------------------------------
      ! That was easy -- the adventure begins now... The nonlinearity!
      if (rproblem%rphysics%dconvectionWeight .ne. 0.0_DP) then
      
        rform%itermCount = 2
        rform%Idescriptors(1,1) = DER_DERIV_X
        rform%Idescriptors(2,1) = DER_FUNC
        rform%Idescriptors(1,2) = DER_DERIV_Y
        rform%Idescriptors(2,2) = DER_FUNC

        rform%ballCoeffConstant = .false.
        rform%BconstantCoeff = .false.
        rform%Dcoefficients(1)  = 1.0
        rform%Dcoefficients(2)  = 1.0
        
        rcollection%DquickAccess(1) = &
            rnonlinearCCMatrix%djota * rproblem%rphysics%dconvectionWeight
        rcollection%p_rvectorQuickAccess1 => rvelocityVector
        
        call bilf_buildMatrixScalar2 (rform, .false., rmatrixTemp,&
            fcoeff_convection,rcollection)

      end if

      call lsyssc_scalarMatVec (rmatrixTemp, &
          rvector%RvectorBlock(4), rdefect%RvectorBlock(4), &
          -1.0_DP, 1.0_DP)
          
      call lsyssc_releaseMatrix (rmatrixTemp)

    end subroutine
    
  end subroutine
  
! ***************************************************************************
  !<subroutine>

  subroutine fcoeff_convection (rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the bilinear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints
    
    ! An array accepting the DOF`s on all elements in the trial space.
    ! DIMENSION(#local DOF`s in trial space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTrial
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
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
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rvector
    real(DP) :: dconvectionWeight
    real(DP), dimension(:,:,:), allocatable :: p_Dvalues
    integer :: i,j
  
    dconvectionWeight = rcollection%DquickAccess(1)
    p_rvector => rcollection%p_RvectorQuickAccess1
    allocate (p_Dvalues(ubound(Dcoefficients,2),ubound(Dcoefficients,3),2))

    call fevl_evaluate_sim (DER_FUNC, p_Dvalues(:,:,1), &
        p_rvector%Rvectorblock(1), Dpoints, &
        rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

    call fevl_evaluate_sim (DER_FUNC, p_Dvalues(:,:,2), &
        p_rvector%Rvectorblock(2), Dpoints, &
        rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
    do j=1,ubound(Dcoefficients,3)
      do i=1,ubound(Dcoefficients,2)
        Dcoefficients(1,i,j) = dconvectionWeight * p_Dvalues(i,j,1)
        Dcoefficients(2,i,j) = dconvectionWeight * p_Dvalues(i,j,2)
      end do
    end do
    
    deallocate (p_Dvalues)

  end subroutine
  

end module
