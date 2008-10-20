!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2matrixassembly </name>
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
!#  $$        A_1 y   +  \eta_1 B p   +          R lambda      = f_1 $$
!#  $$ \tau_1 B^T y   +  \kappa_1 I p                                      = f_2 $$
!#
!#  $$        R_2 y   +  A_2 \lambda  + \eta_2   B \xi           = f_3 $$
!#  $$            \tau_2 B^T \lambda  + \kappa_2 I \xi           = f_4 $$
!#
!# with
!#
!#   $$ A_1 = \iota_1 I  +  \alpha_1 M  +  \theta_1 L  +  \gamma_1 N(y) + dnewton_1 N*(y)$$
!#   $$ A_2 = \iota_2 I  +  \alpha_2 M  +  \theta_2 L  +  \gamma_2 N(y) + dnewton_2 N*(y)$$
!#   $$ R_1 =               \mu_1 (P)M  +            dr_{11} N(\lambda) + dr_{12}   N*(\lambda) $$
!#   $$ R_2 =               \mu_2    M  +            dr_{21} N(\lambda) + dr_{22}   N*(\lambda) $$
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
!#   
!#   $\iota_i$  = 0/1     - switches the identity matrix on/off,
!#   $\alpha_i$ = 0/1     - switches the mass matrix on/off;
!#                          =0 for stationary problem,
!#   $\theta_í$           - weight for the Laplace matrix,
!#   $\gamma_i$ = 0/1     - Switches the nonlinearity on/off;
!#                          =0 for Stokes system,
!#   $dnewton_i \in R$    - Switches the Newton matrix on/off.
!#   $\eta_i$   = 0/1     - Switches the 'B'-term on/off,
!#   $\tau_i$   = 0/1     - Switches the 'B^T'-term on/off,
!#   $\mu_i$              - Weight for the 'coupling' mass matrix.
!#   $\kappa_i$ = 0/1     - Switches of the identity matrix I for the pressure
!#                          in the continuity equation
!#   $dr_ij \in R$        - Switches the 'reactive nonlinearity/coupling mass 
!#                          matrix' on/off
!#   $dp$                 - Switches the mass-type matrix of the derivative of the
!#                          projection on/off
!#
!# Note that the nonlinear part is always dependent on the primal velocity --
!# for the primal equation as well as for the dual one!
!#
!# To assemble such a matrix, the application has to follow two steps:
!# 
!# 1.) Create a structure of type t_ccmatrixComponents and set the
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
!# 1.) cc_assembleMatrix
!#     -> Assembles a matrix based on a set of input parameters.
!#
!# 2.) cc_assembleDefect
!#     -> Set up a defect vector d:=b-A(x)x
!#
!# 3.) cc_projectControlTimestep
!#     Projects a vector into a range of numbers.
!#
!# </purpose>
!##############################################################################

module cc2dmediumm2matvecassembly

  use fsystem
  use storage
  use linearsystemblock
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
  use linearsolverautoinitialise
  use matrixrestriction
  use trilinearformevaluation
  use matrixio
  
  use convection
    
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

  ! This structure provides a set of input parametes for the matrix assembly
  ! routine.
  type t_ccmatrixComponents
  
    ! IOTA-parameters that switch the identity in the primal/dual equation on/off.
    real(DP) :: diota1 = 0.0_DP
    real(DP) :: diota2 = 0.0_DP
    
    ! KAPPA-parameters that switch the I matrix in the continuity equation
    ! on/off.
    real(DP) :: dkappa1 = 0.0_DP
    real(DP) :: dkappa2 = 0.0_DP

    ! ALPHA-parameters that control the weight of the mass matrix in the
    ! core equation. =0.0 for stationary simulations.
    real(DP) :: dalpha1 = 0.0_DP
    real(DP) :: dalpha2 = 0.0_DP
    
    ! THETA-parameters that control the weight of the Stokes matrix
    ! in the core equation. =1.0 for stationary simulations.
    real(DP) :: dtheta1 = 0.0_DP
    real(DP) :: dtheta2 = 0.0_DP
    
    ! GAMMA-parameters that control the weight in front of the
    ! nonlinearity. =1.0 for Navier-Stokes, =0.0 for Stokes equation.
    real(DP) :: dgamma1 = 0.0_DP
    real(DP) :: dgamma2 = 0.0_DP
    
    ! DNEWTON-parameters that control the weight in font of the
    ! newton matrix (adjoint of the nonlinearity). =0.0 to dectivate.
    real(DP) :: dnewton1 = 0.0_DP
    real(DP) :: dnewton2 = 0.0_DP
    
    ! ETA-parameters that switch the B-terms on/off.
    real(DP) :: deta1 = 0.0_DP
    real(DP) :: deta2 = 0.0_DP
    
    ! TAU-parameters that switch the B^T-terms on/off
    real(DP) :: dtau1 = 0.0_DP
    real(DP) :: dtau2 = 0.0_DP
    
    ! MU-parameters that weight the coupling mass matrices.
    real(DP) :: dmu1 = 0.0_DP
    real(DP) :: dmu2 = 0.0_DP

    ! R-parameter that weight the reactive coupling mass matrix
    ! in the primal equation.
    real(DP) :: dr11 = 0.0_DP
    real(DP) :: dr12 = 0.0_DP
    
    ! R-parameter that weight the reactive coupling mass matrix
    ! in the dual equation.
    real(DP) :: dr21 = 0.0_DP
    real(DP) :: dr22 = 0.0_DP
    
    ! Regularisation parameter ALPHA which controls the transformation
    ! from the dual variable $\lambda$ to the control $u$:
    ! $u=-1/\alpha \lambda$.
    real(DP) :: dalphaC = 0.0_DP
    
    ! Type of this matrix.
    ! =0: Standard matrix.
    ! =1: Newton matrix.
    integer :: cmatrixType = 0

    ! Type of constraints to apply to the control u.
    ! =0: No constraints.
    ! =1: Constant constraints on u active: dumin1 <= u_1 <= dumax1, dumin2 <= u_2 <= dumax2.
    ! Necessary for applying the correct projection during matrix-vector
    ! multiplication.
    integer :: ccontrolConstraints = 0
    
    ! Minimum and maximum value of u in case projection of u is
    ! active (i.e. dp <> 0). Min/max values for u1.
    real(DP) :: dumin1 = 0.0_DP
    real(DP) :: dumax1 = 0.0_DP
    
    ! Minimum and maximum value of u in case projection of u is
    ! active (i.e. dp <> 0). Min/max values for u2.
    real(DP) :: dumin2 = 0.0_DP
    real(DP) :: dumax2 = 0.0_DP

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

    ! STABILISATION: Parameter that defines how to set up the nonlinearity and 
    ! whether to use some kind of stabilisation. One of the CCMASM_STAB_xxxx 
    ! constants. Standard is CCMASM_STAB_STREAMLINEDIFF.
    integer :: iupwind1 = CCMASM_STAB_STREAMLINEDIFF
    integer :: iupwind2 = CCMASM_STAB_STREAMLINEDIFF
    
    ! STABILISATION: Viscosity parameter. Used for stabilisation schemes when 
    ! a nonlinearity is set up.
    real(DP) :: dnu = 0.0_DP
    
    ! STABILISATION: Stabilisation parameter for streamline diffusion, upwind and 
    ! edge-oriented stabilisation. If iupwind=CCMASM_STAB_STREAMLINEDIFF, a value of 
    ! 0.0 deactivates any stabilisation.
    real(DP) :: dupsam1 = 0.0_DP
    real(DP) :: dupsam2 = 0.0_DP
    
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
    
    ! An object specifying the block discretisation
    ! (size of subvectors in the solution vector, trial/test functions,...).
    ! Only used during matrix creation.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation => null()

    ! Pointer to a template FEM matrix that defines the structure of 
    ! Laplace/Stokes/... matrices. Only used during matrix creation.
    type(t_matrixScalar), pointer :: p_rmatrixTemplateFEM => null()

    ! A template FEM matrix that defines the structure of gradient
    ! matrices (B1/B2) matrices. Only used during matrix creation.
    type(t_matrixScalar), pointer :: p_rmatrixTemplateGradient => null()

    ! Pointer to Stokes matrix (=nu*Laplace). 
    type(t_matrixScalar), pointer :: p_rmatrixStokes => null()

    ! Pointer to a B1-matrix.
    type(t_matrixScalar), pointer :: p_rmatrixB1 => null()

    ! Pointer to a B2-matrix.
    type(t_matrixScalar), pointer :: p_rmatrixB2 => null()

    ! Pointer to a B1^T-matrix.
    ! This pointer may point to NULL(). In this case, B1^T is created
    ! by 'virtually transposing' the B1 matrix.
    !
    ! Note: This information is automatically created when the preconditioner
    ! is initialised! The main application does not have to initialise it!
    type(t_matrixScalar), pointer :: p_rmatrixB1T => null()

    ! Pointer to a B2-matrix.
    ! This pointer may point to NULL(). In this case, B2^T is created
    ! by 'virtually transposing' the B2 matrix.
    !
    ! Note: This information is automatically created when the preconditioner
    ! is initialised! The main application does not have to initialise it!
    type(t_matrixScalar), pointer :: p_rmatrixB2T => null()

    ! Pointer to a Mass matrix.
    ! May point to NULL() during matrix creation.
    type(t_matrixScalar), pointer :: p_rmatrixMass => null()
    
    ! Pointer to an identity matrix for the pressure.
    type(t_matrixScalar), pointer :: p_rmatrixIdentityPressure => null()

  end type

!</typeblock>

!</types>

contains
  
! ***************************************************************************
  !<subroutine>

  subroutine coeff_ProjMass (rdiscretisationTrial,rdiscretisationTest,rform, &
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
    integer(PREC_ELEMENTIDX), intent(IN)                        :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in trial space,nelements)
    integer(PREC_DOFIDX), dimension(:,:), intent(IN) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    integer(PREC_DOFIDX), dimension(:,:), intent(IN) :: IdofsTest
    
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
    
    celement = rdiscretisationTrial%RelementDistr(&
        rdomainIntSubset%ielementDistribution)%celement
    
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

  ! -----------------------------------------------------

  subroutine massmatfilter (rmatrix, rvector, dalphaC, dmin, dmax)
      
  ! Filters a mass matrix. The lines in the matrix rmatrix corresponding
  ! to all entries in the (control-)vector violating the constraints
  ! of the problem.
  
  ! Matrix to be filtered
  type(t_matrixScalar), intent(inout) :: rmatrix

  ! Vector containing a dial solution lambda. Whereever -1/alpha*lambda
  ! violates the control constraints given by rmatrixComponents, the corresponding
  ! lines are set to 0.
  type(t_vectorScalar), intent(in) :: rvector
  
  ! ALPHA regularisation parameter from the space-time matrix
  real(dp), intent(in) :: dalphaC
  
  ! minimum bound for the control
  real(dp), intent(in) :: dmin

  ! maximum bound for the control
  real(dp), intent(in) :: dmax
  
    ! local variables
    real(dp), dimension(:), pointer :: p_Ddata
    integer(i32), dimension(:), allocatable :: p_Idofs
    integer :: i,nviolate
    real(dp) :: du
    
    ! Get the vector data
    call lsyssc_getbase_double (rvector,p_Ddata)
    
    ! Figure out the DOF's violating the constraints
    allocate(p_Idofs(rvector%NEQ))
    
    nviolate = 0
    do i=1,rvector%NEQ
      du = -p_Ddata(i)/dalphaC
      if ((du .le. dmin) .or. (du .ge. dmax)) then
        nviolate = nviolate + 1
        p_Idofs(nviolate) = i
      end if
    end do
    
    if (nviolate .gt. 0) then
      ! Filter the matrix
      call mmod_replaceLinesByZero (rmatrix,p_Idofs(1:nviolate))
    end if
    
    deallocate(p_Idofs)

  end subroutine    

  ! ***************************************************************************

!<subroutine>

  subroutine cc_assembleMatrix (coperation,cmatrixType,rmatrix,rmatrixComponents,&
      rvector1,rvector2,rvector3,rfineMatrix)

!<description>
  ! This routine assembles a global matrix. The caller must initialise the 
  ! rmatrixComponents according to how the matrix should look like.
  ! The 'coperation' parameter tells the routine what to do.
  ! The destination matrix rmatrix is then set up or updated.
  !
  ! The parameters rvector and rfineMatrix are optional. rvector must be
  ! specified, if the nonlinearity is activated (parameter $\gamma\not=0$ in 
  ! rmatrixComponents). This vector specifies the 'solution' where the nonlinearity 
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
  ! that references to the original gradient (B-)matrices from rmatrixComponents
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

  ! A t_ccmatrixComponents structure providing all necessary 'source' information
  ! about how to set up the matrix. 
  !
  ! Note that if coperation=CCMASM_ALLOCxxxx is specified, p_rmatrixTemplateXXXX
  ! must be initialised as well as p_rdiscretisation!
  ! The new matrix is created based p_rmatrixTemplateXXXX as well as p_rdiscretisation.
  ! Memory is automatically allocated if it's missing.
  type(t_ccmatrixComponents), intent(IN) :: rmatrixComponents

  ! OPTIONAL: If a nonlinearity is to be set up, this vector must be specified.
  ! It specifies where to evaluate the nonlinearity and must contain the data
  ! for the 'previous' timestep. If there is no previous timestep (e.g.
  ! like in the 0th timestep), the vector can be undefined.
  type(t_vectorBlock), intent(IN), optional :: rvector1

  ! OPTIONAL: If a nonlinearity is to be set up, this vector must be specified.
  ! It specifies where to evaluate the nonlinearity and must contain the data
  ! for the 'current' timestep. 
  type(t_vectorBlock), intent(IN), optional :: rvector2

  ! OPTIONAL: If a nonlinearity is to be set up, this vector must be specified.
  ! It specifies where to evaluate the nonlinearity and must contain the data
  ! for the 'next' timestep. If there is no next timestep (e.g.
  ! like in the last timestep), the vector can be undefined.
  type(t_vectorBlock), intent(IN), optional :: rvector3

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
  type(t_matrixBlock), intent(INOUT) :: rmatrix
  
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
    
      ! Allocate the system matrix
      call allocMatrix (cmatrixType,rmatrixComponents,rmatrix)
    end if   
   
    if (iand(coperation,CCMASM_COMPUTE) .ne. 0) then

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
      
      ! Primal equation
      ! ---------------
      ! Assemble rows 1..3 of the block matrix:
      !
      !    ( A11  A12  B1  R1           ) 
      !    ( A21  A22  B2       R2      ) 
      !    ( B1^T B2^T I                )
      
      ! In case dnewton1=0, switch off the A12/A21 matrices as we don't
      ! assemble them.
      if (rmatrixComponents%dnewton1 .ne. 0.0_DP) then
        rmatrix%RmatrixBlock(1,2)%dscaleFactor = 1.0_DP
        rmatrix%RmatrixBlock(2,1)%dscaleFactor = 1.0_DP
      else
        rmatrix%RmatrixBlock(1,2)%dscaleFactor = 0.0_DP
        rmatrix%RmatrixBlock(2,1)%dscaleFactor = 0.0_DP
      end if

      ! 1.) Assemble the velocity submatrix
      !
      !    ( A11  A12   .    .    .    . ) 
      !    ( A21  A22   .    .    .      ) 
      !    (  .    .    .    .    .    . )

      select case (rmatrixComponents%iprimalSol)
      case (1)
        call assembleVelocityBlocks (.false.,&
            rmatrixComponents,rmatrix,rvector1,1.0_DP)
      case (2)
        call assembleVelocityBlocks (.false.,&
            rmatrixComponents,rmatrix,rvector2,1.0_DP)
      case (3)
        call assembleVelocityBlocks (.false.,&
            rmatrixComponents,rmatrix,rvector3,1.0_DP)
      end select
      
      ! Include the mass matrix blocks
      !
      !    (  .    .    .    R1   .    . ) 
      !    (  .    .    .    .    R1   . ) 
      !    (  .    .    .    .    .    . )
      
      select case (rmatrixComponents%idualSol)
      case (1)
        call assembleMassBlocks (rmatrixComponents,rmatrix,rvector1)
      case (2)
        call assembleMassBlocks (rmatrixComponents,rmatrix,rvector2)
      case (3)
        call assembleMassBlocks (rmatrixComponents,rmatrix,rvector3)
      end select
      
      ! 3.) Initialise the weights for the idenity- and B-matrices
      !
      !    (  .    .   B1    .    .    . ) 
      !    (  .    .   B2    .    .    . ) 
      !    ( B1^T B2^T  I    .    .    . )
      
      call assembleGradientMatrices (.false.,rmatrixComponents,rmatrix,&
        iand(coperation,CMASM_QUICKREFERENCES) .ne. 0)

      ! 2.) Initialise the weights for the B-matrices
      !
      !    (  .    .   B1  ) 
      !    (  .    .   B2  ) 
      !    ( B1^T B2^T  I  )
      
      ! Initialise the weights for the B/B^T matrices
      rmatrix%RmatrixBlock(1,3)%dscaleFactor = rmatrixComponents%deta1
      rmatrix%RmatrixBlock(2,3)%dscaleFactor = rmatrixComponents%deta1
      
      rmatrix%RmatrixBlock(3,1)%dscaleFactor = rmatrixComponents%dtau1
      rmatrix%RmatrixBlock(3,2)%dscaleFactor = rmatrixComponents%dtau1
      
      ! Switch the I-matrix in the continuity equation on/off.
      ! The matrix always exists -- but is usually filled with zeroes
      ! or switched off.
      rmatrix%RmatrixBlock(3,3)%dscaleFactor = rmatrixComponents%dkappa1
      if (rmatrixComponents%dkappa1 .ne. 0.0_DP) then
        call lsyssc_initialiseIdentityMatrix (rmatrix%RmatrixBlock(3,3))
      end if
    
      ! Dual equation
      ! ---------------
      ! In case dgamma2=0, switch off the A45/A54 matrices as we haven't
      ! assembled them.
      if (rmatrixComponents%dgamma2 .ne. 0.0_DP) then
        rmatrix%RmatrixBlock(4,5)%dscaleFactor = 1.0_DP
        rmatrix%RmatrixBlock(5,4)%dscaleFactor = 1.0_DP
      else
        rmatrix%RmatrixBlock(4,5)%dscaleFactor = 0.0_DP
        rmatrix%RmatrixBlock(5,4)%dscaleFactor = 0.0_DP
      end if
      
      ! Assemble rows 1..3 of the block matrix:
      !
      !    (  R2   .    .  A44  A45  B1 ) 
      !    (  .    R2   .  A54  A55  B2 ) 
      !    (  .    .    .  B1^T B2^T I  ) 
      !
      ! 1.) Assemble the velocity submatrix
      !
      !    (  .    .    .   A44  A45   . ) 
      !    (  .    .    .   A54  A55   . ) 
      !    (  .    .    .    .    .    . )
      !
      ! Note that the Newton part (if assembled) is set up with the
      ! primal velocity! The weight is taken from the GAMMA parameter
      ! and switches of the assembly of Newton. Actually, it's not
      ! a 'Newton' matrix but the adjoint of the nonlinearity
      ! which is assembled here...
      !
      ! Specify the correct primal solution for evaluating the nonlinear
      ! dual matrices.
      select case (rmatrixComponents%iprimalSol)
      case (1)
        call assembleVelocityBlocks (.true.,&
            rmatrixComponents,rmatrix,rvector1,1.0_DP)
      case (2)
        call assembleVelocityBlocks (.true.,&
            rmatrixComponents,rmatrix,rvector2,1.0_DP)
      case (3)
        call assembleVelocityBlocks (.true.,&
            rmatrixComponents,rmatrix,rvector3,1.0_DP)
      end select

      ! 2.) Include the mass matrix blocks
      !
      !    (  R2   .    .    .    .    .  ) 
      !    (  .    R2   .    .    .    .  ) 
      !    (  .    .    .    .    .    .  )
      
      ! Specify the correct dual solution for a possible evaluationp of
      ! nonlinearity in R.
      select case (rmatrixComponents%idualSol)
      case (1)
        call assembleDualMassBlocks (rmatrixComponents,rmatrix,rvector1)
        !CALL assembleMassBlocks (4,1,rmatrixComponents,rmatrix)
      case (2)
        call assembleDualMassBlocks (rmatrixComponents,rmatrix,rvector2)
      case (3)
        call assembleDualMassBlocks (rmatrixComponents,rmatrix,rvector3)
      end select
      
      ! 3.) Initialise the weights for the idenity- and B-matrices
      !
      !    ( .    .    .    .    .   B1 ) 
      !    ( .    .    .    .    .   B2 ) 
      !    ( .    .    .   B1^T B2^T  I )
      
      ! The B/B^T-matrices themself share their structure and data with
      ! those of the primal equation. This is possible as long as the
      ! type of boundary conditions are the same, as then implementing the
      ! BC's into the B/B^T-matrices of the primal equation will affect
      ! those of the dual equation as well and vice versa
      ! (in total the BC's are then implemented 2x into the same matrices
      ! which needs a little bit more time but not much).
      !
      ! Furthermore, the VANKA preconditioner can only handle the
      ! situation where the B/B^T matrices of the dual equation share their
      ! data with that of the primal one...
      call lsyssc_duplicateMatrix ( &
          rmatrix%Rmatrixblock(1,3),rmatrix%Rmatrixblock(4,6), &
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      call lsyssc_duplicateMatrix ( &
          rmatrix%Rmatrixblock(2,3),rmatrix%Rmatrixblock(5,6), &
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      call lsyssc_duplicateMatrix ( &
          rmatrix%Rmatrixblock(3,1),rmatrix%Rmatrixblock(6,4), &
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      call lsyssc_duplicateMatrix ( &
          rmatrix%Rmatrixblock(3,2),rmatrix%Rmatrixblock(6,5), &
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
  
      ! Initialise the weights for the B/B^T matrices
      rmatrix%RmatrixBlock(4,6)%dscaleFactor = rmatrixComponents%deta2
      rmatrix%RmatrixBlock(5,6)%dscaleFactor = rmatrixComponents%deta2
      
      rmatrix%RmatrixBlock(6,4)%dscaleFactor = rmatrixComponents%dtau2
      rmatrix%RmatrixBlock(6,5)%dscaleFactor = rmatrixComponents%dtau2
      
      ! Initialise the weights of the mass matrices
      !rmatrix%RmatrixBlock(4,1)%dscaleFactor = rmatrixComponents%dmu2
      !rmatrix%RmatrixBlock(5,2)%dscaleFactor = rmatrixComponents%dmu2
      
      ! Switch the I-matrix in the dual continuity equation on/off.
      ! The matrix always exists -- but is usually filled with zeroes
      ! or switched off.
      rmatrix%RmatrixBlock(6,6)%dscaleFactor = rmatrixComponents%dkappa2
      if (rmatrixComponents%dkappa2 .ne. 0.0_DP) then
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
            rmatrixComponents%iadaptiveMatrices, &
            rmatrixComponents%dadmatthreshold)
            
        if (.not. lsyssc_isMatrixContentShared(&
            rfineMatrix%RmatrixBlock(1,1),rfineMatrix%RmatrixBlock(2,2))) then
          call mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(2,2), &
              rmatrix%RmatrixBlock(1,1), &
              rmatrixComponents%iadaptiveMatrices, &
              rmatrixComponents%dadmatthreshold)
        end if
          
        ! Dual system
        call mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(4,4), &
            rmatrix%RmatrixBlock(4,4), &
            rmatrixComponents%iadaptiveMatrices, &
            rmatrixComponents%dadmatthreshold)
            
        if (.not. lsyssc_isMatrixContentShared(&
          rmatrix%RmatrixBlock(4,4),rmatrix%RmatrixBlock(5,5))) then
          call mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(5,5), &
              rmatrix%RmatrixBlock(5,5), &
              rmatrixComponents%iadaptiveMatrices, &
              rmatrixComponents%dadmatthreshold)
        end if
        
      end if

    end if
    
  contains
  
    ! -----------------------------------------------------
  
    subroutine allocMatrix (cmatrixType,rmatrixComponents,rmatrix)
    
    ! Allocates memory for the system matrix. rmatrixComponents provides information
    ! about the submatrices that are 'plugged into' rmatrix.
    ! Therefore, before this routine is called, rmatrixComponents must have been set up.

    ! Type of matrix that should be set up in rmatrix. One of the CCMASM_MTP_xxxx
    ! constants.
    integer, intent(IN) :: cmatrixType

    ! A t_ccmatrixComponents structure providing all necessary 'source' information
    ! about how to set up the matrix. 
    type(t_ccmatrixComponents), intent(IN), target :: rmatrixComponents

    ! A block matrix that receives the basic system matrix.
    type(t_matrixBlock), intent(INOUT) :: rmatrix
    
      ! local variables
      logical :: bdecoupled, bfulltensor

      ! A pointer to the system matrix and the RHS/solution vectors.
      type(t_matrixScalar), pointer :: p_rmatrixTemplateFEM,p_rmatrixTemplateGradient

      ! A pointer to the discretisation structure with the data.
      type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rmatrixComponents%p_rdiscretisation
      
      ! Get a pointer to the template FEM matrix. If that doesn't exist,
      ! take the Stokes matrix as template.
      p_rmatrixTemplateFEM => rmatrixComponents%p_rmatrixTemplateFEM
      if (.not. associated(p_rmatrixTemplateFEM)) &
        p_rmatrixTemplateFEM => rmatrixComponents%p_rmatrixStokes
      if (.not. associated(p_rmatrixTemplateFEM)) then
        print *,'allocMatrix: Cannot set up A matrices in system matrix!'
        stop
      end if

      ! In the global system, there are two gradient matrices B1 and B2.
      ! Get a pointer to the template structure for these.
      ! If there is no pointer, try to get use a pointer to one of these
      ! matrices directly.
      p_rmatrixTemplateGradient => rmatrixComponents%p_rmatrixTemplateGradient
      if (.not. associated(p_rmatrixTemplateGradient)) &
        p_rmatrixTemplateGradient => rmatrixComponents%p_rmatrixB1
      if (.not. associated(p_rmatrixTemplateGradient)) then
        print *,'allocMatrix: Cannot set up B matrix in system matrix!'
        stop
      end if

      ! Initialise the block matrix with default values based on
      ! the discretisation.
      if (associated(p_rdiscretisation)) then
        call lsysbl_createMatBlockByDiscr (p_rdiscretisation,rmatrix)    
      else
        ! No discretisation structure; create the matrix directly as 3x3 matrix.
        call lsysbl_createEmptyMatrix (rmatrix,NDIM2D+1)
      end if

      ! -------------------------------------------------------------
      ! Primal equation
      ! -------------------------------------------------------------
        
      ! Determine the shape of the matrix
      bdecoupled = cmatrixType .eq. CCMASM_MTP_DECOUPLED
      bfulltensor = cmatrixType .eq. CCMASM_MTP_FULLTENSOR
      
      if (cmatrixType .eq. CCMASM_MTP_AUTOMATIC) then
        ! Should we assemble Newton? If yes, we have a full-tensor matrix.
        bfulltensor = rmatrixComponents%dnewton1 .ne. 0.0_DP
      end if
    
      ! Let's consider the global system in detail. The standard matrix It has 
      ! roughly the following shape:
      !
      !    ( A11  .    B1  M   .   .   ) = ( A11  .    A13 A14 .   .   )
      !    ( .    A22  B2  .   M   .   )   ( .    A22  A23 .   A25 .   )
      !    ( B1^T B2^T .   .   .   .   )   ( A31  A32  .   .   .   .   )
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
      call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                  rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          
      if (.not. bfulltensor) then     
           
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
      call lsyssc_assignDiscretDirectMat (rmatrix%RmatrixBlock(2,2),&
          p_rdiscretisation%RspatialDiscr(2))
                                       
      ! A 'full tensor matrix' consists also of blocks A12 and A21.
      if (bfulltensor) then

        ! We have a matrix in the following shape:
        !
        !    ( A11  A12  B1  R1  .   .   ) = ( A11  A12  A13 A14 .   .   )
        !    ( A21  A22  B2  .   R1  .   )   ( A21  A22  A23 .   A25 .   )
        !    ( B1^T B2^T .   .   .   .   )   ( A31  A32  .   .   .   .   )
        !
        ! Create A12 and A21.
      
        if (rmatrix%RmatrixBlock(1,2)%cmatrixFormat &
            .eq. LSYSSC_MATRIXUNDEFINED) then
            
          call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM, &
            rmatrix%RmatrixBlock(1,2), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
            
          ! Allocate memory for the entries; don't initialise the memory.
          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
          ! zero.
          ! CALL lsyssc_allocEmptyMatrix (&
          !     rmatrix%RmatrixBlock(1,2),LSYSSC_SETM_UNDEFINED)
            
        end if

        if (rmatrix%RmatrixBlock(2,1)%cmatrixFormat &
            .eq. LSYSSC_MATRIXUNDEFINED) then
            
          ! Create a new matrix A21 in memory. create a new matrix
          ! using the template FEM matrix...
          call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM, &
            rmatrix%RmatrixBlock(2,1), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
            
          ! Allocate memory for the entries; don't initialise the memory.
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
      call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixB1, &
                                    rmatrix%RmatrixBlock(1,3),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixB2, &
                                    rmatrix%RmatrixBlock(2,3),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        
      ! In the same manner, insert an identiy matrix for the pressure
      ! to the system matrix. The matrix is zero by default but may
      ! partially or fully be initialised to an identity matrix depending
      ! on the situation. (It's mostly used for direct solvers/UMFPACK)
      call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixIdentityPressure, &
                                    rmatrix%RmatrixBlock(3,3),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      call lsyssc_clearMatrix (rmatrix%RmatrixBlock(3,3))

      ! Now, prepare B1^T and B2^T. Are these matrices given?
      !
      ! If yes, take the given matrices. If not, create them by
      ! 'virtually transposing' B1 and B2 (i.e. they share the same
      ! data as B1 and B2 but hate the 'transpose'-flag set).
      
      if (associated(rmatrixComponents%p_rmatrixB1T)) then
        call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixB1T, &
                                      rmatrix%RmatrixBlock(3,1),&
                                      LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      else
        call lsyssc_transposeMatrix (rmatrixComponents%p_rmatrixB1, &
                                      rmatrix%RmatrixBlock(3,1),&
                                      LSYSSC_TR_VIRTUAL)
      end if

      if (associated(rmatrixComponents%p_rmatrixB2T)) then
        call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixB2T, &
                                      rmatrix%RmatrixBlock(3,2),&
                                      LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      else
        call lsyssc_transposeMatrix (rmatrixComponents%p_rmatrixB2, &
                                      rmatrix%RmatrixBlock(3,2),&
                                      LSYSSC_TR_VIRTUAL)
      end if

      ! Insert free space for the mass matrices to the matrix.
      !
      ! Note that we share the structure of M, while we create empty space 
      ! for the entries. 
      ! Later, the M-matrices are copied into here and modified for boundary
      ! conditions.
      call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass, &
                                    rmatrix%RmatrixBlock(1,4),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass, &
                                    rmatrix%RmatrixBlock(2,5),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
                                    
      ! -------------------------------------------------------------
      ! Dual equation
      ! -------------------------------------------------------------
        
      ! Determine the shape of the matrix
      bdecoupled = cmatrixType .eq. CCMASM_MTP_DECOUPLED
      bfulltensor = cmatrixType .eq. CCMASM_MTP_FULLTENSOR
      
      if (cmatrixType .eq. CCMASM_MTP_AUTOMATIC) then
        ! Should we assemble Newton? If yes, we have a full-tensor matrix.
        bfulltensor = rmatrixComponents%dnewton2 .ne. 0.0_DP
      end if
    
      ! Let's consider the global system in detail. The standard matrix It has 
      ! roughly the following shape:
      !
      !    ( R2   .    .   A44  .    B1  ) = ( A41  .    .   A44 .   A46 )
      !    ( .    R2   .   .    A55  B2  )   ( .    A52  .   .   A55 A56 )
      !    ( .    .    .   B1^T B2^T .   )   ( .    .    .   A64 A65 .   )
      !
      ! All matices may have multiplication factors in their front.
      !
      ! The structure of the matrices A44 and A55 of the global system matrix 
      ! is governed by the template FEM matrix.
      ! Initialise them with the same structure, i.e. A44, A55 share (!) their
      ! structure (not the entries) with that of the template matrix.
      !
      ! For this purpose, use the "duplicate matrix" routine.
      ! The structure of the matrix is shared with the template FEM matrix.
      ! For the content, a new empty array is allocated which will later receive
      ! the entries.
      call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                  rmatrix%RmatrixBlock(4,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          
      if (.not. bdecoupled .and. .not. bfulltensor) then     
           
        ! If X- and Y-velocity is to be treated in a 'coupled' way, the matrix 
        ! A22 is identical to A44! So mirror A44 to A55 sharing the
        ! structure and the content.
        call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(4,4),&
                    rmatrix%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
                    
      else
      
        ! Otherwise, create another copy of the template matrix.
        call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                    rmatrix%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
                    
      end if
      
      ! Manually change the discretisation structure of the Y-velocity 
      ! matrix to the Y-discretisation structure.
      ! Ok, we use the same discretisation structure for both, X- and Y-velocity,
      ! so this is not really necessary - we do this for sure...
      call lsyssc_assignDiscretDirectMat (rmatrix%RmatrixBlock(5,5),&
          p_rdiscretisation%RspatialDiscr(5))
                                          
      ! A 'full tensor matrix' consists also of blocks A12 and A21.
      if (bfulltensor) then

        ! We have a submatrix in the following shape:
        !
        !    ( R2   .    .   A44  A45  B1  ) = ( A41  .    .   A44 A45 A46 )
        !    ( .    R2   .   A54  A55  B2  )   ( .    A52  .   A54 A55 A56 )
        !    ( .    .    .   B1^T B2^T .   )   ( .    .    .   A64 A65 .   )
        !
        ! Create A45 and A54.
      
        if (rmatrix%RmatrixBlock(4,5)%cmatrixFormat &
            .eq. LSYSSC_MATRIXUNDEFINED) then
            
          call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM, &
            rmatrix%RmatrixBlock(4,5), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
            
          ! Allocate memory for the entries; don't initialise the memory.
          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
          ! zero.
          ! CALL lsyssc_allocEmptyMatrix (&
          !     rmatrix%RmatrixBlock(4,5),LSYSSC_SETM_UNDEFINED)
            
        end if

        if (rmatrix%RmatrixBlock(5,4)%cmatrixFormat &
            .eq. LSYSSC_MATRIXUNDEFINED) then
            
          ! Create a new matrix A54 in memory. create a new matrix
          ! using the template FEM matrix...
          call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM, &
            rmatrix%RmatrixBlock(5,4), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
            
          ! Allocate memory for the entries; don't initialise the memory.
          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
          ! zero.
          ! CALL lsyssc_allocEmptyMatrix (&
          !     p_rmatrixPreconditioner%RmatrixBlock(5,4),LSYSSC_SETM_UNDEFINED)
            
        end if
        
      end if

      ! The B1/B2 matrices exist up to now only in our local problem structure.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2 with those B1/B2 of the
      ! block matrix, while we create empty space for the entries. 
      ! Later, the B-matrices are copied into here and modified for boundary
      ! conditions.
      call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixB1, &
                                    rmatrix%RmatrixBlock(4,6),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixB2, &
                                    rmatrix%RmatrixBlock(5,6),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        
      ! Now, prepare B1^T and B2^T. Are these matrices given?
      !
      ! If yes, take the given matrices. If not, create them by
      ! 'virtually transposing' B1 and B2 (i.e. they share the same
      ! data as B1 and B2 but hate the 'transpose'-flag set).
      
      if (associated(rmatrixComponents%p_rmatrixB1T)) then
        call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixB1T, &
                                      rmatrix%RmatrixBlock(6,4),&
                                      LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      else
        call lsyssc_transposeMatrix (rmatrixComponents%p_rmatrixB1, &
                                      rmatrix%RmatrixBlock(6,4),&
                                      LSYSSC_TR_VIRTUAL)
      end if

      if (associated(rmatrixComponents%p_rmatrixB2T)) then
        call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixB2T, &
                                      rmatrix%RmatrixBlock(6,5),&
                                      LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      else
        call lsyssc_transposeMatrix (rmatrixComponents%p_rmatrixB2, &
                                      rmatrix%RmatrixBlock(6,5),&
                                      LSYSSC_TR_VIRTUAL)
      end if

      ! In the same manner, insert an identiy matrix for the pressure
      ! to the system matrix; as the enties aren't changed, we can share
      ! the entries with the original one.
      call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixIdentityPressure, &
                                    rmatrix%RmatrixBlock(6,6),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      call lsyssc_clearMatrix (rmatrix%RmatrixBlock(6,6))
      
      ! Insert free space for the mass matrices to the matrix.
      !
      ! Note that we share the structure of M, while we create empty space 
      ! for the entries. 
      ! Later, the M-matrices are copied into here and modified for boundary
      ! conditions.
      call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass, &
                                    rmatrix%RmatrixBlock(4,1),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass, &
                                    rmatrix%RmatrixBlock(5,2),&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      ! That's it, all submatrices are basically set up.
      !
      ! Update the structural information of the block matrix, as we manually
      ! changed the submatrices:
      call lsysbl_updateMatStrucInfo (rmatrix)
        
    end subroutine
    
    ! -----------------------------------------------------

    subroutine assembleVelocityBlocks (bdualEquation,&
        rmatrixComponents,rmatrix,rvector,dvectorWeight)
        
    ! Assembles the velocity matrix in the block matrix rmatrix at position (1,1):
    !
    ! rmatrix := dalpha*M + dtheta*Laplace + dgamma*N(p_rvector) +
    !            dnewton*N*(p_rvector)
    
    ! Whether to set up the primal or the dual equation.
    ! FALSE=primal, TRUE=dual equation.
    logical, intent(IN) :: bdualEquation
    
    ! A t_ccmatrixComponents structure providing all necessary 'source' information
    ! about how to set up the matrix. 
    type(t_ccmatrixComponents), intent(IN) :: rmatrixComponents
    
    ! Block matrix where the 2x2-velocity submatrix should be assembled
    type(t_matrixBlock), intent(INOUT) :: rmatrix
    
    ! Velocity vector for the nonlinearity. Must be specified if
    ! GAMMA <> 0; can be omitted if GAMMA=0.
    type(t_vectorBlock), optional :: rvector
    
    ! Weight for the velocity vector; standard = 1.
    real(DP), intent(IN), optional :: dvectorWeight
    
    ! local variables
    logical :: bshared
    integer :: iupwind
    real(DP) :: dupsam
    type(t_convUpwind) :: rupwind
    type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
    type(t_jumpStabilisation) :: rjumpStabil
    type(t_matrixBlock) :: rtempMatrix
    real(DP) :: dvecWeight
    integer :: imatOffset
    real(DP) :: diota, dalpha, dtheta, dnewton, dgamma
    
      if (.not. bdualEquation) then
        ! Set the weights used here according to the primal equation.
        ! Set imatOffset=0 so the submatrix at position 1,1 is tackled.
        imatOffset = 0
        diota = rmatrixComponents%diota1
        dalpha = rmatrixComponents%dalpha1
        dtheta = rmatrixComponents%dtheta1
        dgamma = rmatrixComponents%dgamma1
        dnewton = rmatrixComponents%dnewton1
        iupwind = rmatrixComponents%iupwind1
        dupsam = rmatrixComponents%dupsam1
      else
        ! Set the weights used here according to the primal equation.
        ! Set imatOffset=3 so the submatrix at position 4,4 is tackled.
        imatOffset = 3
        diota = rmatrixComponents%diota2
        dalpha = rmatrixComponents%dalpha2
        dtheta = rmatrixComponents%dtheta2
        dgamma = rmatrixComponents%dgamma2
        dnewton = rmatrixComponents%dnewton2
        iupwind = rmatrixComponents%iupwind2
        dupsam = rmatrixComponents%dupsam2
      end if

      ! Standard value for dvectorWeight is = -1.
      dvecWeight = -1.0_DP
      if (present(dvectorWeight)) dvecWeight = dvectorWeight
    
      ! Is A11=A22 physically?
      bshared = lsyssc_isMatrixContentShared(&
                    rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1),&
                    rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2))
                    
      ! Allocate memory if necessary. Normally this should not be necessary...
      ! A11:
      if (.not. lsyssc_hasMatrixContent (&
          rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1))) then
        call lsyssc_allocEmptyMatrix (&
            rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1),LSYSSC_SETM_UNDEFINED)
      end if
    
      ! A22:
      if (.not. bshared) then
        if (.not. lsyssc_hasMatrixContent (&
            rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2))) then
          call lsyssc_allocEmptyMatrix (&
              rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2),LSYSSC_SETM_UNDEFINED)
        end if
      end if

      ! A12/ A21:
      if (lsysbl_isSubmatrixPresent (rmatrix,imatOffset+1,imatOffset+2)) then
        if (.not. lsyssc_hasMatrixContent (&
            rmatrix%RmatrixBlock(imatOffset+1,imatOffset+2))) then
          call lsyssc_allocEmptyMatrix (&
              rmatrix%RmatrixBlock(imatOffset+1,imatOffset+2),LSYSSC_SETM_UNDEFINED)
        end if
        if (.not. lsyssc_hasMatrixContent (&
            rmatrix%RmatrixBlock(imatOffset+2,imatOffset+1))) then
          call lsyssc_allocEmptyMatrix (&
              rmatrix%RmatrixBlock(imatOffset+2,imatOffset+1),LSYSSC_SETM_UNDEFINED)
        end if
      end if
    
      ! ---------------------------------------------------
      ! If diota <> 0, initialise the matrix with the identity,
      ! otherwise with zero.
      if (diota .eq. 0.0_DP) then
        call lsyssc_clearMatrix (rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1))
        
        if (.not. bshared) then
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2))
        end if
        
      else
        
        call lsyssc_initialiseIdentityMatrix (&
            rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1))
        call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1),diota)
        
        if (.not. bshared) then
          call lsyssc_initialiseIdentityMatrix (&
              rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2))
          call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2),diota)
        end if
        
      end if
    
      ! ---------------------------------------------------
      ! Plug in the mass matrix?
      if (dalpha .ne. 0.0_DP) then
       
        call lsyssc_matrixLinearComb (&
            rmatrixComponents%p_rmatrixMass  ,dalpha,&
            rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1),1.0_DP,&
            rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1),&
            .false.,.false.,.true.,.true.)
            
        if (.not. bshared) then

          call lsyssc_matrixLinearComb (&
              rmatrixComponents%p_rmatrixMass     ,dalpha,&
              rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2),1.0_DP,&
              rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2),&
              .false.,.false.,.true.,.true.)
        end if
        
      end if
      
      ! ---------------------------------------------------
      ! Plug in the Stokes matrix?
      if (dtheta .ne. 0.0_DP) then
        call lsyssc_matrixLinearComb (&
            rmatrixComponents%p_rmatrixStokes     ,dtheta,&
            rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1),1.0_DP,&
            rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1),&
            .false.,.false.,.true.,.true.)
            
        if (.not. bshared) then
          call lsyssc_matrixLinearComb (&
              rmatrixComponents%p_rmatrixStokes   ,dtheta,&
              rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2),1.0_DP,&
              rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2),&
              .false.,.false.,.true.,.true.)
        end if
      end if
      
      ! ---------------------------------------------------
      ! That was easy -- the adventure begins now... The nonlinearity!
      if (dgamma .ne. 0.0_DP) then
      
        if (.not. present(rvector)) then
          call output_line ('Velocity vector not present!', &
                             OU_CLASS_ERROR,OU_MODE_STD,'cc_assembleMatrix')
          call sys_halt()
        end if
      
        select case (iupwind)
        case (CCMASM_STAB_STREAMLINEDIFF)
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rmatrixComponents%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion%dupsam = dupsam
          
          ! Matrix weight
          rstreamlineDiffusion%ddelta = dgamma
          
          ! Weight for the Newton part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = dnewton
          
          if (dnewton .eq. 0.0_DP) then
          
            ! If the submatrices A12 and A21 exist, fill them with zero.
            ! If they don't exist, we don't have to do anything.
            if (lsysbl_isSubmatrixPresent (rmatrix,imatOffset+1,imatOffset+2)) then
              call lsyssc_clearMatrix (rmatrix%RmatrixBlock(imatOffset+1,imatOffset+2))
              call lsyssc_clearMatrix (rmatrix%RmatrixBlock(imatOffset+2,imatOffset+1))
            end if
            
          else

            ! Clear A12/A21 that may receive parts of the Newton matrix
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(imatOffset+1,imatOffset+2))
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(imatOffset+2,imatOffset+1))
          
          end if

          ! Create a temporary block matrix only contining the velocity submatrices
          ! we want to change. Share structure and entries such that changing
          ! the temporary matrix will also change the original matrix.
          call lsysbl_deriveSubmatrix (rmatrix,rtempMatrix,&
                                       LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,&
                                       imatOffset+1,imatOffset+2)
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvector!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamlineDiffusionBlk2d (&
                              rvector, rvector, &
                              dvecWeight, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODMATRIX, &
                              rtempMatrix)
                              
          ! Release the temp matrix.
          call lsysbl_releaseMatrix (rtempMatrix)

        case (CCMASM_STAB_UPWIND)
          ! Set up the upwind structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rupwind%dnu = rmatrixComponents%dnu
          
          ! Set stabilisation parameter
          rupwind%dupsam = dupsam

          ! Matrix weight
          rupwind%dtheta = dgamma
          
          if (dnewton .ne. 0.0_DP) then
            call output_line ('Warning: Upwind does not support assembly '&
                //'of the Newton matrix!',OU_CLASS_TRACE1)
          end if
          
          ! Call the upwind method to calculate the nonlinear matrix.
          ! As velocity vector, specify rvector!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_upwind2d (rvector, rvector, &
                              dvecWeight, 0.0_DP,&
                              rupwind, CONV_MODMATRIX, &
                              rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1)) 
                              
          if (.not. bshared) then
            ! Modify also the matrix block (2,2)
            call conv_upwind2d (rvector, rvector, &
                                dvecWeight, 0.0_DP,&
                                rupwind, CONV_MODMATRIX, &
                                rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2)) 
          end if     

        case (CCMASM_STAB_EDGEORIENTED)
          ! Jump stabilisation.
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rmatrixComponents%dnu
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion%ddelta = dgamma
          
          ! Weight for the Newtop part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = dnewton
          
          if (dnewton .eq. 0.0_DP) then

            ! If the submatrices A12 and A21 exist, fill them with zero.
            ! If they don't exist, we don't have to do anything.
            if (lsysbl_isSubmatrixPresent (rmatrix,imatOffset+1,imatOffset+2)) then
              call lsyssc_clearMatrix (rmatrix%RmatrixBlock(imatOffset+1,imatOffset+2))
              call lsyssc_clearMatrix (rmatrix%RmatrixBlock(imatOffset+2,imatOffset+1))
            end if
            
          else

            ! Clear A12/A21 that receives parts of the Newton matrix
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(imatOffset+1,imatOffset+2))
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(imatOffset+2,imatOffset+1))
          
            ! Activate the submatrices A12 and A21 if they aren't.
            rmatrix%RmatrixBlock(imatOffset+1,imatOffset+2)%dscaleFactor = 1.0_DP
            rmatrix%RmatrixBlock(imatOffset+2,imatOffset+1)%dscaleFactor = 1.0_DP
           
          end if
         
          ! Create a temporary block matrix only contining the velocity submatrices
          ! we want to change. Share structure and entries such that changing
          ! the temporary matrix will also change the original matrix.
          call lsysbl_deriveSubmatrix (rmatrix,rtempMatrix,&
                                       LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,&
                                       imatOffset+1,imatOffset+2)
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvector!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamlineDiffusionBlk2d (&
                              rvector, rvector, &
                              dvecWeight, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODMATRIX, &
                              rtempMatrix)
                             
          ! Release the temp matrix.
          call lsysbl_releaseMatrix (rtempMatrix)          
        
          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rstreamlineDiffusion%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgammastar = dupsam
          rjumpStabil%dgamma = rjumpStabil%dgammastar
          
          ! Matrix weight
          rjumpStabil%dtheta = dgamma

          ! Call the jump stabilisation technique to stabilise that stuff.   
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
                              rjumpStabil, CONV_MODMATRIX, &
                              rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1))   

          if (.not. bshared) then
            call conv_jumpStabilisation2d (&
                                rjumpStabil, CONV_MODMATRIX, &
                                rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2))   
          end if

        case DEFAULT
          print *,'Don''t know how to set up nonlinearity!?!'
          stop
        
        end select

      else
      
        ! That's the Stokes-case. Jump stabilisation is possible...
      
        select case (iupwind)
        case (CCMASM_STAB_EDGEORIENTED)
          
          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rmatrixComponents%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgammastar = dupsam
          rjumpStabil%dgamma = rjumpStabil%dgammastar
          
          ! Matrix weight
          rjumpStabil%dtheta = dgamma

          ! Call the jump stabilisation technique to stabilise that stuff.   
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
                              rjumpStabil, CONV_MODMATRIX, &
                              rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1))   

          if (.not. bshared) then
            call conv_jumpStabilisation2d (&
                                rjumpStabil, CONV_MODMATRIX, &
                                rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2))   
          end if

        case DEFAULT
          ! No stabilisation
        
        end select
      
      end if ! gamma <> 0

    end subroutine  
      
    ! -----------------------------------------------------

    subroutine assembleMassBlocks (rmatrixComponents,rmatrix, rvector)
        
    ! Assembles a 2x2 block matrix with mass matrices and
    ! probably nonlinear submatrices on the diagonal.
    
    ! A t_ccmatrixComponents structure providing all necessary 'source' information
    ! about how to set up the matrix. 
    type(t_ccmatrixComponents), intent(IN) :: rmatrixComponents
    
    ! Block matrix where the 2x2-velocity submatrix should be assembled
    type(t_matrixBlock), intent(INOUT) :: rmatrix
    
    ! Vector that specifies where to evaluate nonlinear terms
    type(t_vectorBlock), intent(IN), target :: rvector
    
      ! local variables
      type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
      type(t_matrixBlock) :: rtempmatrix
      type(t_vectorBlock) :: rtempvector
      type(t_bilinearForm) :: rform
      type(t_collection) :: rcollection

      ! Assemble A14/A25? 
      if (rmatrixComponents%dmu1 .ne. 0.0_DP) then
          
        ! Calculate the usual mass matrix if conrol constraints are deactivated
        ! or if Newton is not active.
        if ((rmatrixComponents%ccontrolConstraints .eq. 0) .or. &
            (rmatrixComponents%cmatrixType .eq. 0)) then
      
          ! Copy the entries of the mass matrix. Share the structure.
          ! We must not share the entries as these might be changed by the caller
          ! e.g. due to boundary conditions!
          
          call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
              rmatrix%RmatrixBlock(1,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        
          call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
              rmatrix%RmatrixBlock(2,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
              
          ! Scale the entries by dmu2.
          if (rmatrixComponents%dmu1 .ne. 1.0_DP) then
            call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(1,4),rmatrixComponents%dmu1)
            call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(2,5),rmatrixComponents%dmu1)
          end if
          
        else if ((rmatrixComponents%ccontrolConstraints .eq. 1) .and. &
                 (rmatrixComponents%cmatrixType .eq. 1)) then
          
!          ! In A14/A25 we have to create a 'projective mass matrix'.
!          ! This is the derivative of a projection operator
!          ! P[a,b](f)=a if f<a, =b if f>b, =f otherwise.
!          ! For a<f<b, this is the mass matrix. Everywhere else, this is =0.
!          ! We assemble this matrix just as a standard mass matrix with noconstant
!          ! coefficients. Whereever u = -1/alpha * lambda is out of bounds,
!          ! we return 0 as coefficient, otherwise 1.
!        
!          rform%itermCount = 1
!          rform%Idescriptors(1,1) = DER_FUNC
!          rform%Idescriptors(2,1) = DER_FUNC
!
!          ! In this case, we have nonconstant coefficients.
!          rform%ballCoeffConstant = .FALSE.
!          rform%BconstantCoeff(:) = .FALSE.
!
!          ! Prepare a collection structure to be passed to the callback
!          ! routine. We attach the vector T in the quick-access variables
!          ! so the callback routine can access it.
!          ! The bounds and the alpha value are passed in the
!          ! quickaccess-arrays.
!          call collct_init(rcollection)
!          rcollection%p_rvectorQuickAccess1 => rvector
!
!          ! Coefficient is dmu1=1/alpha or 0, depending on lambda
!          rcollection%DquickAccess(3)  = rmatrixComponents%dalphaC
!          rcollection%DquickAccess(4)  = rmatrixComponents%dmu1
!          
!          ! At first, set up A14, depending on lambda_1.
!          rcollection%IquickAccess(1) = 1
!          rcollection%DquickAccess(1) = rmatrixComponents%dumin1
!          rcollection%DquickAccess(2) = rmatrixComponents%dumax1
!          
!          ! Now we can build the matrix entries.
!          ! We specify the callback function coeff_Laplace for the coefficients.
!          ! As long as we use constant coefficients, this routine is not used.
!          ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
!          ! the framework will call the callback routine to get analytical
!          ! data.
!          ! The collection is passed as additional parameter. That's the way
!          ! how we get the vector to the callback routine.
!          call bilf_buildMatrixScalar (rform,.TRUE.,rmatrix%RmatrixBlock(1,4),&
!              coeff_ProjMass,rcollection)
!
!          ! Now, set up A25, depending on lambda_2.
!          rcollection%IquickAccess(1)  = 2
!          rcollection%DquickAccess(1) = rmatrixComponents%dumin2
!          rcollection%DquickAccess(2) = rmatrixComponents%dumax2
!
!          call bilf_buildMatrixScalar (rform,.TRUE.,rmatrix%RmatrixBlock(2,5),&
!              coeff_ProjMass,rcollection)
!          
!          ! Now we can forget about the collection again.
!          call collct_done (rcollection)

          ! Copy the entries of the mass matrix. Share the structure.
          ! We must not share the entries as these might be changed by the caller
          ! e.g. due to boundary conditions!
          
          call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
              rmatrix%RmatrixBlock(1,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        
          call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
              rmatrix%RmatrixBlock(2,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
              
          ! Scale the entries by dmu2.
          if (rmatrixComponents%dmu1 .ne. 1.0_DP) then
            call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(1,4),rmatrixComponents%dmu1)
            call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(2,5),rmatrixComponents%dmu1)
          end if

          ! Filter the matrix. All the rows corresponding to DOF's that violate
          ! the bounds must be set to zero.
          call massmatfilter (rmatrix%RmatrixBlock(1,4),rvector%RvectorBlock(4),&
              rmatrixComponents%dalphaC,rmatrixComponents%dumin1,rmatrixComponents%dumax1)
          call massmatfilter (rmatrix%RmatrixBlock(2,5),rvector%RvectorBlock(5),&
              rmatrixComponents%dalphaC,rmatrixComponents%dumin2,rmatrixComponents%dumax2)
        end if
        
        rmatrix%RmatrixBlock(1,4)%dscaleFactor = 1.0_DP
        rmatrix%RmatrixBlock(2,5)%dscaleFactor = 1.0_DP

        if (rmatrixComponents%dr12 .ne. 0.0_DP) then
          ! There is some data in A24/A15, so create empty space there
          ! in case it's missing.
          if (.not. lsysbl_isSubmatrixPresent(rmatrix,2,4)) then
            call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
                rmatrix%RmatrixBlock(2,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
                
            call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
                rmatrix%RmatrixBlock(1,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          end if

          ! Clear the offdiagonal matrices, switch them on
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,5))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,4))

          rmatrix%RmatrixBlock(1,5)%dscaleFactor = 1.0_DP
          rmatrix%RmatrixBlock(2,4)%dscaleFactor = 1.0_DP
        
        else

          rmatrix%RmatrixBlock(1,5)%dscaleFactor = 0.0_DP
          rmatrix%RmatrixBlock(2,4)%dscaleFactor = 0.0_DP
          
        end if
        
      else if ((rmatrixComponents%dr11 .ne. 0.0_DP) .or. &
               (rmatrixComponents%dr12 .ne. 0.0_DP)) then
        
        ! There is some data in A14/A25, so create empty space there
        ! in case it's missing.
        if (.not. lsysbl_isSubmatrixPresent(rmatrix,1,4)) then
          call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
              rmatrix%RmatrixBlock(1,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
              
          call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
              rmatrix%RmatrixBlock(2,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        end if
        
        ! Clear the diagonal matrices, switch them on
        call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,4))
        call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,5))

        rmatrix%RmatrixBlock(1,4)%dscaleFactor = 1.0_DP
        rmatrix%RmatrixBlock(2,5)%dscaleFactor = 1.0_DP
        
        if (rmatrixComponents%dr12 .ne. 0.0_DP) then
          ! There is some data in A42/A51, so create empty space there
          ! in case it's missing.
          if (.not. lsysbl_isSubmatrixPresent(rmatrix,2,4)) then
            call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
                rmatrix%RmatrixBlock(2,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
                
            call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
                rmatrix%RmatrixBlock(1,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          end if

          ! Clear the offdiagonal matrices, switch them on
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,5))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,4))

          rmatrix%RmatrixBlock(1,5)%dscaleFactor = 1.0_DP
          rmatrix%RmatrixBlock(2,4)%dscaleFactor = 1.0_DP
        
        else

          rmatrix%RmatrixBlock(1,5)%dscaleFactor = 0.0_DP
          rmatrix%RmatrixBlock(2,4)%dscaleFactor = 0.0_DP
         
        end if

      else
      
        ! Deactivate this submatrix
        rmatrix%RmatrixBlock(1,4)%dscaleFactor = 0.0_DP
        rmatrix%RmatrixBlock(1,5)%dscaleFactor = 0.0_DP
        rmatrix%RmatrixBlock(2,4)%dscaleFactor = 0.0_DP
        rmatrix%RmatrixBlock(2,5)%dscaleFactor = 0.0_DP

      end if
            
      ! If we have a reactive coupling mass matrix, it gets interesting...
      if ((rmatrixComponents%dr11 .ne. 0.0_DP) .or. &
          (rmatrixComponents%dr12 .ne. 0.0_DP)) then
      
        ! The reactive part is: "dr2 * . * grad(lambda)".
        ! This is exactly the 'Newton' part assembled by the streamline diffusion
        ! method if we use the dual velocity as velocity field!
        ! So prepare to call streamline diffusion.

        ! Viscosity; ok, actually not used.
        rstreamlineDiffusion%dnu = rmatrixComponents%dnu
        
        ! Set stabilisation parameter
        rstreamlineDiffusion%dupsam = rmatrixComponents%dupsam2
        
        ! Weight dr21 of the convective part.
        rstreamlineDiffusion%ddelta = rmatrixComponents%dr11
        
        ! Weight for the Newton part; here, this is the dr22 weight.
        rstreamlineDiffusion%dnewton = rmatrixComponents%dr12
        
        ! Create a temporary block matrix only contining the velocity submatrices
        ! we want to change. Share structure and entries such that changing
        ! the temporary matrix will also change the original matrix.
        call lsysbl_deriveSubmatrix (rmatrix,rtempMatrix,&
                                      LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,1,2,4,5)
                                      
        ! Create a temporary block vector that points to the dual velocity.
        call lsysbl_deriveSubvector (rvector,rtempvector,4,5,.true.)
        
        ! Call the SD method to calculate the nonlinearity.
        ! As velocity vector, specify rtempvector, which points to
        ! the dual velocity.
        call conv_streamlineDiffusionBlk2d (&
                            rtempvector, rtempvector, &
                            1.0_DP, 0.0_DP,&
                            rstreamlineDiffusion, CONV_MODMATRIX, &
                            rtempMatrix)

        ! Release the temp matrix and temp vector.
        call lsysbl_releaseVector (rtempVector)
        call lsysbl_releaseMatrix (rtempMatrix)
      
      end if
    
    end subroutine
    
    ! -----------------------------------------------------

    subroutine assembleDualMassBlocks (rmatrixComponents,rmatrix,rvector)
        
    ! Assembles a 2x2 block matrix with mass matrices on the diagonal
    ! in the dual equation. These matrices consist of a standard mass
    ! matrix and/or a reactive coupling mass matrix depending on the
    ! dual velocity.
    
    ! A t_ccmatrixComponents structure providing all necessary 'source' information
    ! about how to set up the matrix. 
    type(t_ccmatrixComponents), intent(IN) :: rmatrixComponents
    
    ! Block matrix where the 2x2-velocity submatrix should be assembled
    type(t_matrixBlock), intent(INOUT) :: rmatrix
    
    ! Vector that specifies where to evaluate nonlinear terms
    type(t_vectorBlock), intent(IN) :: rvector
    
      ! local variables
      type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
      type(t_matrixBlock) :: rtempmatrix
      type(t_vectorBlock) :: rtempvector
    
      if (rmatrixComponents%dmu2 .ne. 0.0_DP) then
    
        ! Copy the entries of the mass matrix. Share the structure.
        ! We must not share the entries as these might be changed by the caller
        ! e.g. due to boundary conditions!
        
        call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
            rmatrix%RmatrixBlock(4,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      
        call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
            rmatrix%RmatrixBlock(5,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
            
        ! Scale the entries by dmu2.
        if (rmatrixComponents%dmu2 .ne. 1.0_DP) then
          call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(4,1),rmatrixComponents%dmu2)
          call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(5,2),rmatrixComponents%dmu2)
        end if
        
        rmatrix%RmatrixBlock(4,1)%dscaleFactor = 1.0_DP
        rmatrix%RmatrixBlock(5,2)%dscaleFactor = 1.0_DP

        if (rmatrixComponents%dr22 .ne. 0.0_DP) then
          ! There is some data in A42/A51, so create empty space there
          ! in case it's missing.
          if (.not. lsysbl_isSubmatrixPresent(rmatrix,4,2)) then
            call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
                rmatrix%RmatrixBlock(4,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
                
            call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
                rmatrix%RmatrixBlock(5,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          end if

          ! Clear the offdiagonal matrices, switch them on
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(5,1))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(4,2))

          rmatrix%RmatrixBlock(5,1)%dscaleFactor = 1.0_DP
          rmatrix%RmatrixBlock(4,2)%dscaleFactor = 1.0_DP
        
        else

          rmatrix%RmatrixBlock(5,1)%dscaleFactor = 0.0_DP
          rmatrix%RmatrixBlock(4,2)%dscaleFactor = 0.0_DP
          
        end if
        
      else if ((rmatrixComponents%dr21 .ne. 0.0_DP) .or. &
               (rmatrixComponents%dr22 .ne. 0.0_DP)) then
        
        ! There is some data in A41/A52, so create empty space there
        ! in case it's missing.
        if (.not. lsysbl_isSubmatrixPresent(rmatrix,4,1)) then
          call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
              rmatrix%RmatrixBlock(4,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
              
          call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
              rmatrix%RmatrixBlock(5,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        end if
        
        ! Clear the diagonal matrices, switch them on
        call lsyssc_clearMatrix (rmatrix%RmatrixBlock(4,1))
        call lsyssc_clearMatrix (rmatrix%RmatrixBlock(5,2))

        rmatrix%RmatrixBlock(4,1)%dscaleFactor = 1.0_DP
        rmatrix%RmatrixBlock(5,2)%dscaleFactor = 1.0_DP
        
        if (rmatrixComponents%dr22 .ne. 0.0_DP) then
          ! There is some data in A42/A51, so create empty space there
          ! in case it's missing.
          if (.not. lsysbl_isSubmatrixPresent(rmatrix,4,2)) then
            call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
                rmatrix%RmatrixBlock(4,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
                
            call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
                rmatrix%RmatrixBlock(5,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          end if

          ! Clear the offdiagonal matrices, switch them on
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(5,1))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(4,2))

          rmatrix%RmatrixBlock(5,1)%dscaleFactor = 1.0_DP
          rmatrix%RmatrixBlock(4,2)%dscaleFactor = 1.0_DP
        
        else

          rmatrix%RmatrixBlock(5,1)%dscaleFactor = 0.0_DP
          rmatrix%RmatrixBlock(4,2)%dscaleFactor = 0.0_DP
         
        end if

      else
      
        ! Deactivate this submatrix
        rmatrix%RmatrixBlock(4,1)%dscaleFactor = 0.0_DP
        rmatrix%RmatrixBlock(5,1)%dscaleFactor = 0.0_DP
        rmatrix%RmatrixBlock(4,2)%dscaleFactor = 0.0_DP
        rmatrix%RmatrixBlock(5,2)%dscaleFactor = 0.0_DP

      end if
            
      ! If we have a reactive coupling mass matrix, it gets interesting...
      if ((rmatrixComponents%dr21 .ne. 0.0_DP) .or. &
          (rmatrixComponents%dr22 .ne. 0.0_DP)) then
      
        ! The reactive part is: "dr2 * . * grad(lambda)".
        ! This is exactly the 'Newton' part assembled by the streamline diffusion
        ! method if we use the dual velocity as velocity field!
        ! So prepare to call streamline diffusion.

        ! Viscosity; ok, actually not used.
        rstreamlineDiffusion%dnu = rmatrixComponents%dnu
        
        ! Set stabilisation parameter
        rstreamlineDiffusion%dupsam = rmatrixComponents%dupsam2
        
        ! Weight dr21 of the convective part.
        rstreamlineDiffusion%ddelta = rmatrixComponents%dr21
        
        ! Weight for the Newton part; here, this is the dr22 weight.
        rstreamlineDiffusion%dnewton = rmatrixComponents%dr22
        
        ! Create a temporary block matrix only contining the velocity submatrices
        ! we want to change. Share structure and entries such that changing
        ! the temporary matrix will also change the original matrix.
        call lsysbl_deriveSubmatrix (rmatrix,rtempMatrix,&
                                      LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,4,5,1,2)
                                      
        ! Create a temporary block vector that points to the dual velocity.
        call lsysbl_deriveSubvector (rvector,rtempvector,4,5,.true.)
        
        ! Call the SD method to calculate the nonlinearity.
        ! As velocity vector, specify rtempvector, which points to
        ! the dual velocity.
        call conv_streamlineDiffusionBlk2d (&
                            rtempvector, rtempvector, &
                            1.0_DP, 0.0_DP,&
                            rstreamlineDiffusion, CONV_MODMATRIX, &
                            rtempMatrix)

        ! Release the temp matrix and temp vector.
        call lsysbl_releaseVector (rtempVector)
        call lsysbl_releaseMatrix (rtempMatrix)
      
      end if
      
    end subroutine
    
    ! -----------------------------------------------------
    
    subroutine assembleGradientMatrices (bdualEquation,&
        rmatrixComponents,rmatrix,bsharedMatrix)
    
    ! Initialises the gradient/divergence matrices with entries from
    ! the rmatrixComponents structure.
    !
    ! The routine copies references from the submatrices tormatrix,
    ! but it does not initialise any matrix weights / scaling factors.
    !
    ! If bsharedMatrix=TRUE, the matrix is created using references to the
    ! matrix building blocks in rlevelInfo, thus sharing all information
    ! with those matrices in rmatrixComponents. In this case, the caller must
    ! not change the matrix entries, because this would change the
    ! original 'template' matrices!
    ! (This can be used e.g. for setting up a matrix for building a defect
    !  vector without copying matrix data.)
    ! If bsharedMatrix=TRUE on the other hand, the matrix entries of the
    ! original template (B-) matrices are copied in memory,
    ! so the new matrix is allowed to be changed!

    ! Whether to set up the proimal or the dual equation.
    ! FALSE=primal, TRUE=dual equation.
    logical, intent(IN) :: bdualEquation

    ! A t_ccmatrixComponents structure providing all necessary 'source' information
    ! about how to set up the matrix. 
    type(t_ccmatrixComponents), intent(IN) :: rmatrixComponents

    ! Block matrix where the B-matrices should be set up
    type(t_matrixBlock), intent(INOUT) :: rmatrix

    ! Whether or not the matrix entries of the source gradient-matrices 
    ! should be copied in memory. 
    !
    ! If set to TRUE, the routine tries to initialise rmatrix
    ! only with references to the original matrices, thus the caller must not
    ! change the entries. Nevertheless, if rmatrix is the owner of one of the
    ! submatrices, the routine will always copy the matrix entries,
    ! as otherwise memory would have to be deallocated!
    !
    ! If set to FALSE, the entries of the source matrices in rmatrixComponents are
    ! copied, so the caller can change rmatrix afterwards (e.g. to implement
    ! boundary conditions).
    logical, intent(IN) :: bsharedMatrix

      ! local variables
      integer :: idubStructure,idubContent,imatOffset
      
      if (.not. bdualEquation) then
        ! Set imatOffset=0 so the submatrix at position 1,1 is tackled.
        imatOffset = 0
      else
        ! Set the weights used here according to the primal equation.
        imatOffset = 3
      end if

      ! Initialise a copy flag that tells the duplicateMatrix-routine whether to
      ! copy the entries or to create references.
      if (bsharedMatrix) then
      
        idubContent = LSYSSC_DUP_SHARE
        
        ! Normally we share entries -- except for if the submatrices belong to
        ! rmatrix! To avoid memory deallocation in this case, we copy
        ! the entries.
        if ((.not. lsyssc_isMatrixContentShared(&
                Rmatrix%RmatrixBlock(imatOffset+1,imatOffset+3))) .or.&
            (.not. lsyssc_isMatrixContentShared(&
                Rmatrix%RmatrixBlock(imatOffset+2,imatOffset+3)))) then
          idubContent = LSYSSC_DUP_COPY
        end if
        
      else
      
        idubContent = LSYSSC_DUP_COPY
        
      end if
      
      idubStructure = LSYSSC_DUP_SHARE
      
      ! Let's consider the global system in detail:
      !
      !    ( A11  A12  B1  ) = ( A11  A12  A13 )
      !    ( A21  A22  B2  )   ( A21  A22  A23 )
      !    ( B1^T B2^T 0   )   ( A31  A32  A33 )
      !
      ! We exclude the velocity submatrices here, so our system looks like:
      !
      !    (           B1 ) = (           A13 )
      !    (           B2 )   (           A23 )
      !    ( B1^T B2^T    )   ( A31  A32      )

      ! The B1/B2 matrices exist up to now only in rmatrixComponents.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2 with those B1/B2 of the
      ! block matrix, while we create copies of the entries. The B-blocks
      ! are already prepared and memory for the entries is already allocated;
      ! so we only have to copy the entries.
      !
      ! Note that idubContent = LSYSSC_DUP_COPY will automatically allocate
      ! memory if necessary.
      call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixB1, &
                                    rmatrix%RmatrixBlock(imatOffset+1,imatOffset+3),&
                                    idubStructure,idubContent)

      call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixB2, &
                                    rmatrix%RmatrixBlock(imatOffset+2,imatOffset+3),&
                                    idubStructure,idubContent)
      
      ! Furthermore, put B1^T and B2^T to the block matrix.
      ! These matrices are always 'shared'.
      call lsyssc_transposeMatrix (rmatrixComponents%p_rmatrixB1, &
                                    rmatrix%RmatrixBlock(imatOffset+3,imatOffset+1),&
                                    LSYSSC_TR_VIRTUAL)

      call lsyssc_transposeMatrix (rmatrixComponents%p_rmatrixB2, &
                                    rmatrix%RmatrixBlock(imatOffset+3,imatOffset+2),&
                                    LSYSSC_TR_VIRTUAL)

    end subroutine

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_assembleDefect (rmatrixComponents,rx,rd,cx,rvector1,rvector2,rvector3)

!<description>
  ! This routine assembles the nonlinear defect
  !      rd := rd - cx A(ry) rx
  ! with the system matrix A(.) defined by the configuration in rmatrixComponents.
  ! The caller must initialise the rmatrixComponents according to how the 
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

  ! A t_ccmatrixComponents structure providing all necessary 'source' information
  ! about how to set up the matrix. 
  !
  ! The caller must provide either p_rmatrixTemplateXXXX in this structure
  ! or set the p_rmatrixTemplateXXXX as well as p_rdiscretisation to
  ! appropriate values. This is necessary for exploiting then structure
  ! of the matrix.
  type(t_ccmatrixComponents), intent(IN) :: rmatrixComponents

  ! This vector specifies the 'x' that is multiplied to the matrix.
  type(t_vectorBlock), intent(IN), target :: rx

  ! OPTIONAL: Multiplication factor in front of the term 'A(ry) rx'.
  ! If not specified, cx=1.0 is assumed.
  real(DP), intent(IN), optional :: cx

  ! OPTIONAL: If a nonlinearity is to be set up, this vector must be specified.
  ! It specifies where to evaluate the nonlinearity and must contain the data
  ! for the 'previous' timestep. If there is no previous timestep (e.g.
  ! like in the 0th timestep), the vector can be undefined.
  type(t_vectorBlock), intent(IN), optional :: rvector1

  ! OPTIONAL: If a nonlinearity is to be set up, this vector must be specified.
  ! It specifies where to evaluate the nonlinearity and must contain the data
  ! for the 'current' timestep. 
  type(t_vectorBlock), intent(IN), optional :: rvector2

  ! OPTIONAL: If a nonlinearity is to be set up, this vector must be specified.
  ! It specifies where to evaluate the nonlinearity and must contain the data
  ! for the 'next' timestep. If there is no next timestep (e.g.
  ! like in the last timestep), the vector can be undefined.
  type(t_vectorBlock), intent(IN), optional :: rvector3

!</input>

!<inputoutput>
  ! Destination vector. cx*A(ry)*rx is subtracted from this vector.
  type(t_vectorBlock), intent(INOUT) :: rd
!</inputoutput>
  
!</subroutine>

    ! local variables
    real(DP) :: dcx
    type(t_matrixBlock) :: rmatrix
    type(t_vectorScalar) :: rtempVector
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Dd
    
    call lsysbl_getbase_double (rd,p_Dd)
    
    dcx = 1.0_DP
    if (present(cx)) dcx = cx
    
    ! The system matrix looks like:
    !          
    !    ( A11  A12  B1  M            ) 
    !    ( A21  A22  B2       M       ) 
    !    ( B1^T B2^T                  ) 
    !    ( M             A44  A45  B1 ) 
    !    (      M        A54  A55  B2 ) 
    !    (               B1^T B2^T    ) 
    ! 
    ! Create a temporary matrix that covers this structure.
    call lsysbl_createEmptyMatrix (rmatrix,2*(NDIM2D+1))
    
    ! Put references to the Stokes- and B-matrices to Aij. assembleVelocityDefect 
    ! needs this template matrix to provide the structure for the stabilisation
    ! routines! The B-matrices are needed later.
    ! -----
    ! Primal equation
    call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixStokes,&
        rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixStokes,&
        rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    
    if (rmatrixComponents%dnewton1 .ne. 0.0_DP) then
      call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixStokes,&
          rmatrix%RmatrixBlock(1,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixStokes,&
          rmatrix%RmatrixBlock(2,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    end if
    
    call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixB1,&
        rmatrix%RmatrixBlock(1,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixB2,&
        rmatrix%RmatrixBlock(2,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

    call lsyssc_transposeMatrix (rmatrixComponents%p_rmatrixB1, &
                                  rmatrix%RmatrixBlock(3,1),&
                                  LSYSSC_TR_VIRTUAL)

    call lsyssc_transposeMatrix (rmatrixComponents%p_rmatrixB2, &
                                  rmatrix%RmatrixBlock(3,2),&
                                  LSYSSC_TR_VIRTUAL)

    call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
        rmatrix%RmatrixBlock(1,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
        rmatrix%RmatrixBlock(2,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

    ! -----
    ! Dual equation
    call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixStokes,&
        rmatrix%RmatrixBlock(4,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixStokes,&
        rmatrix%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    
    if (rmatrixComponents%dnewton2 .ne. 0.0_DP) then
      call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixStokes,&
          rmatrix%RmatrixBlock(4,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixStokes,&
          rmatrix%RmatrixBlock(5,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    end if
    
    call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixB1,&
        rmatrix%RmatrixBlock(4,6),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixB2,&
        rmatrix%RmatrixBlock(5,6),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

    call lsyssc_transposeMatrix (rmatrixComponents%p_rmatrixB1, &
                                  rmatrix%RmatrixBlock(6,4),&
                                  LSYSSC_TR_VIRTUAL)

    call lsyssc_transposeMatrix (rmatrixComponents%p_rmatrixB2, &
                                  rmatrix%RmatrixBlock(6,5),&
                                  LSYSSC_TR_VIRTUAL)

    call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
        rmatrix%RmatrixBlock(4,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
        rmatrix%RmatrixBlock(5,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

    ! Update the structural information of the block matrix, as we manually
    ! changed the submatrices:
    call lsysbl_updateMatStrucInfo (rmatrix)
    
    ! In the first step, we assemble the defect that arises in the velocity 
    ! components. This is characterised by the following submatrix:
    !
    !    ( A11  A12                   ) 
    !    ( A21  A22                   ) 
    !    (                            ) 
    !    (               A44  A45     ) 
    !    (               A54  A55     ) 
    !    (                            ) 
    ! 
    ! assembleVelocityDefect handles exactly these submatices.
    ! We call the routine twice -- once for the primal and once for the dual
    ! equation. In both cases, we specify the primal velocity p_ry
    ! as velocity field (!).

    select case (rmatrixComponents%iprimalSol)
    case (1)
      call assembleVelocityDefect (&
          .false.,rmatrixComponents,rmatrix,rx,rd,dcx,rvector1,1.0_DP)
      call assembleVelocityDefect (&
          .true.,rmatrixComponents,rmatrix,rx,rd,dcx,rvector1,1.0_DP)
    case (2)
      call assembleVelocityDefect (&
          .false.,rmatrixComponents,rmatrix,rx,rd,dcx,rvector2,1.0_DP)
      call assembleVelocityDefect (&
          .true.,rmatrixComponents,rmatrix,rx,rd,dcx,rvector2,1.0_DP)
    case (3)
      call assembleVelocityDefect (&
          .false.,rmatrixComponents,rmatrix,rx,rd,dcx,rvector3,1.0_DP)
      call assembleVelocityDefect (&
          .true.,rmatrixComponents,rmatrix,rx,rd,dcx,rvector3,1.0_DP)
    end select
    
    ! Now, we treat all the remaining blocks. Let's see what is missing:
    !
    !    ( .    .    B1  M             ) 
    !    ( .    .    B2       M        ) 
    !    ( B1^T B2^T .                 ) 
    !    ( M             .    .    B1  ) 
    !    (      M        .    .    B2  ) 
    !    (               B1^T B2^T .   ) 

    ! To build the appropriate defect, we firat remove the velocity blocks:
    
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(1,1))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(1,2))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(2,1))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(2,2))
    
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(4,4))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(4,5))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(5,4))
    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(5,5))

    ! Initialise the weights for the B/B^T matrices
    rmatrix%RmatrixBlock(1,3)%dscaleFactor = rmatrixComponents%deta1
    rmatrix%RmatrixBlock(2,3)%dscaleFactor = rmatrixComponents%deta1
    
    rmatrix%RmatrixBlock(3,1)%dscaleFactor = rmatrixComponents%dtau1
    rmatrix%RmatrixBlock(3,2)%dscaleFactor = rmatrixComponents%dtau1

    rmatrix%RmatrixBlock(4,6)%dscaleFactor = rmatrixComponents%deta2
    rmatrix%RmatrixBlock(5,6)%dscaleFactor = rmatrixComponents%deta2
    
    rmatrix%RmatrixBlock(6,4)%dscaleFactor = rmatrixComponents%dtau2
    rmatrix%RmatrixBlock(6,5)%dscaleFactor = rmatrixComponents%dtau2
    
    ! Initialise the weights for the mass matrices
    if (rmatrixComponents%ccontrolConstraints .eq. 0) then
      rmatrix%RmatrixBlock(1,4)%dscaleFactor = rmatrixComponents%dmu1
      rmatrix%RmatrixBlock(2,5)%dscaleFactor = rmatrixComponents%dmu1
    else
      rmatrix%RmatrixBlock(1,4)%dscaleFactor = 0.0_DP
      rmatrix%RmatrixBlock(2,5)%dscaleFactor = 0.0_DP
    end if

    rmatrix%RmatrixBlock(4,1)%dscaleFactor = rmatrixComponents%dmu2
    rmatrix%RmatrixBlock(5,2)%dscaleFactor = rmatrixComponents%dmu2

    ! ------------------------------------------------
    ! Build the defect by matrix-vector multiplication
    !
    ! Note that no time step or whatever is included here; everything
    ! is initialised with the multiplication factors in the submatrices
    ! from above!
    call lsysbl_blockMatVec (rmatrix, rx, rd, -dcx, 1.0_DP)
    
    ! The coupling of lambda to the primal equation is a little bit tricky.
    ! Ok, if there are no constraints active, it's easy -- that case was handled
    ! in the MV above...
    if ((rmatrixComponents%ccontrolConstraints .eq. 1) .and. &
        (rmatrixComponents%dmu1 .ne. 0.0_DP)) then
    
      ! But now it get's interesting. When control constraints are active,
      ! we have the system
      !
      !    ( .    .    B1  PM            ) 
      !    ( .    .    B2       PM       ) 
      !    ( B1^T B2^T .                 ) 
      !    ( M             .    .    B1  ) 
      !    (      M        .    .    B2  ) 
      !    (               B1^T B2^T .   ) 
      !
      ! Where P is a projection operator onto the allowed space. The dual variable
      ! \lambda must not be changed, but before adding M\lambda to the RHS,
      ! we have to apply a projection!
      !
      ! That's a little bit ugly because for this we have to carry out the matrix
      ! vector multiplication 'by hand'. 
      

      ! What's the type of the current matrix? Is this a Newton-matrix?
      if (rmatrixComponents%cmatrixType .eq. 0) then

        ! No, this is a standard matrix. That means, we just have to project
        ! the control u and multiply it with the mass matrix.

!      Projection of the solution. Test code. Not used as the inequality
!      defining the coupling must be used in integral form -- for what
!      the implementation below that is correct.

!      ! At first, create a temporary vector that
!      ! receives the projected solution.
!      call lsyssc_duplicateVector (rx%RvectorBlock(4),rtempVector,&
!        LSYSSC_DUP_COPY,LSYSSC_DUP_EMPTY)
!      
!      ! Multiply with the mass matrix, correctly scaled.
!      call lsyssc_scalarMatVec (rmatrixComponents%p_rmatrixMass, rx%RvectorBlock(4), &
!        rtempVector, rmatrixComponents%dmu1, 0.0_DP)
!      
!      ! Project that to the allowed range.
!      call cc_projectControlTimestep (rtempVector,&
!          -rmatrixComponents%dumax1,-rmatrixComponents%dumin1)
!          
!      ! And then finally, carry our the defect calculation for y_1.
!      call lsyssc_vectorLinearComb (rtempVector,rd%RvectorBlock(1),-dcx,1.0_DP)
!      
!      ! The same stuff has to be done for y_1 / lambda_2:
!      call lsyssc_scalarMatVec (rmatrixComponents%p_rmatrixMass, rx%RvectorBlock(5), &
!        rtempVector, rmatrixComponents%dmu1, 0.0_DP)
!      call cc_projectControlTimestep (rtempVector,&
!          -rmatrixComponents%dumax2,-rmatrixComponents%dumin2)
!      call lsyssc_vectorLinearComb (rtempVector,rd%RvectorBlock(2),-dcx,1.0_DP)

        ! Copy our solution vector \lambda_1. Scale it by -1/alpha.
        call lsyssc_duplicateVector (rx%RvectorBlock(4),rtempVector,&
            LSYSSC_DUP_COPYOVERWRITE,LSYSSC_DUP_COPYOVERWRITE)
            
        call lsyssc_scaleVector (rtempVector,-rmatrixComponents%dmu1)
        
        ! Project that to the allowed range to create u_1.
        call cc_projectControlTimestep (rtempVector,&
            rmatrixComponents%dumin1,rmatrixComponents%dumax1)

        ! Now carry out MV and include it to the defect.
        ! Here, we apply an MV where we include again dmu1 into the coefficient.
        call lsyssc_scalarMatVec (rmatrixComponents%p_rmatrixMass, rtempVector, &
            rd%RvectorBlock(1), dcx, 1.0_DP)

        ! The same stuff has to be done for y_2 / lambda_2:
        call lsyssc_duplicateVector (rx%RvectorBlock(5),rtempVector,&
            LSYSSC_DUP_COPYOVERWRITE,LSYSSC_DUP_COPYOVERWRITE)
        call lsyssc_scaleVector (rtempVector,-rmatrixComponents%dmu1)
        call cc_projectControlTimestep (rtempVector,&
            rmatrixComponents%dumin2,rmatrixComponents%dumax2)
        call lsyssc_scalarMatVec (rmatrixComponents%p_rmatrixMass, rtempVector, &
          rd%RvectorBlock(2), dcx, 1.0_DP)

        ! Release the temp vector, that's it.
        call lsyssc_releaseVector (rtempVector)

      else

        ! Yes, that's a Newton matrix. That means, we have to multiply the
        ! vector with the derivative of the projection operator:
        ! b-P[a,b]'(lambda).
        ! For that purpose, we have to assemble special mass matrices:
      
        select case (rmatrixComponents%idualSol)
        case (1)
          call assemblePrimalUConstrMassDefect (rmatrixComponents,rx,&
              rd,dcx,rvector1)
        case (2)
          call assemblePrimalUConstrMassDefect (rmatrixComponents,rx,&
              rd,dcx,rvector2)
        case (3)
          call assemblePrimalUConstrMassDefect (rmatrixComponents,rx,&
              rd,dcx,rvector3)
        end select
        
      end if

    end if
    
    ! If we have a reactive coupling mass matrix, it gets interesting...
    !
    ! Primal equation
    if ((rmatrixComponents%dr11 .ne. 0.0_DP) .or. &
        (rmatrixComponents%dr12 .ne. 0.0_DP)) then
    
      ! Assemble the defect of the reactive coupling mass matrices
      select case (rmatrixComponents%idualSol)
      case (1)
        call assemblePrimalReactMassDefect (rmatrixComponents,rx,&
            rd,dcx,rvector1,1.0_DP)
      case (2)
        call assemblePrimalReactMassDefect (rmatrixComponents,rx,&
            rd,dcx,rvector2,1.0_DP)
      case (3)
        call assemblePrimalReactMassDefect (rmatrixComponents,rx,&
            rd,dcx,rvector3,1.0_DP)
      end select
    
    end if

    ! Dual equation
    if ((rmatrixComponents%dr21 .ne. 0.0_DP) .or. &
        (rmatrixComponents%dr22 .ne. 0.0_DP)) then
    
      ! Assemble the defect of the reactive coupling mass matrices
      select case (rmatrixComponents%idualSol)
      case (1)
        call assembleDualReactMassDefect (rmatrixComponents,rx,&
            rd,dcx,rvector1,1.0_DP)
      case (2)
        call assembleDualReactMassDefect (rmatrixComponents,rx,&
            rd,dcx,rvector2,1.0_DP)
      case (3)
        call assembleDualReactMassDefect (rmatrixComponents,rx,&
            rd,dcx,rvector3,1.0_DP)
      end select
    
    end if
    
    ! Release the temporary matrix, we don't need it anymore.
    call lsysbl_releaseMatrix (rmatrix)

  contains

    subroutine assembleVelocityDefect (bdualEquation,rmatrixComponents,&
        rmatrix,rvector,rdefect,dcx,rvelocityVector,dvectorWeight)
        
    ! Assembles the velocity defect in the block matrix rmatrix at position
    ! itop..itop+1 in the velocity vector. rdefect must have been initialised
    ! with the right hand side vector.
    !
    ! With a matrix 'A' of the theoretical form
    !
    !       A := dalpha*M + dtheta*Laplace + dgamma*N(p_rvector) +
    !            dnewton*N*(p_rvector)
    !
    ! the routine will construct
    !
    !       rdefect = rdefect - dcx * (dtheta A rvector)
    
    ! Whether to set up the primal or the dual equation.
    ! FALSE=primal, TRUE=dual equation.
    logical, intent(IN) :: bdualEquation
        
    ! A t_ccmatrixComponents structure providing all necessary 'source' information
    ! about how to set up the matrix. 
    type(t_ccmatrixComponents), intent(IN) :: rmatrixComponents

    ! Reference to the system matrix. Only the structure of the matrix
    ! is used to reconstruct the structure of the discretisation.
    ! The content of the matrix is not changed or used.
    type(t_matrixBlock), intent(INOUT) :: rmatrix
    
    ! Solution vector.
    type(t_vectorBlock), intent(IN) :: rvector
    
    ! On entry: RHS vector.
    ! Is overwritten by the defect vector in the velocity subsystem.
    type(t_vectorBlock), intent(INOUT) :: rdefect
    
    ! Multiplication factor for the whole operator A*rvector
    real(DP), intent(IN) :: dcx
    
    ! Weight for the velocity vector rvelocityVector; usually = 1.0
    real(DP), intent(IN) :: dvectorWeight
    
    ! Velocity vector field that should be used for the assembly of the
    ! nonlinearity. The first two blocks in that block vector are
    ! used as velocity field.
    type(t_vectorBlock), intent(IN) :: rvelocityVector

    ! local variables
    logical :: bshared
    integer :: iupwind,dupsam
    type(t_convUpwind) :: rupwind
    type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
    type(T_jumpStabilisation) :: rjumpStabil
    type(t_vectorBlock) :: rtempVector,rtempDefect
    type(t_matrixBlock) :: rtempMatrix
    integer :: imatOffset
    real(DP) :: dalpha, dtheta, dnewton, diota, dgamma
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Dd
    
    call lsysbl_getbase_double (rdefect,p_Dd)
    
      if (.not. bdualEquation) then
        ! Set the weights used here according to the primal equation.
        ! Set imatOffset=0 so the submatrix at position 1,1 is tackled.
        imatOffset = 0
        dalpha = rmatrixComponents%dalpha1
        dtheta = rmatrixComponents%dtheta1 
        diota  = rmatrixComponents%diota1
        dgamma = rmatrixComponents%dgamma1 
        dnewton = rmatrixComponents%dnewton1
        iupwind = rmatrixComponents%iupwind1
        dupsam = rmatrixComponents%dupsam1
      else
        ! Set the weights used here according to the primal equation.
        ! Set imatOffset=3 so the submatrix at position 4,4 is tackled.
        imatOffset = 3
        dalpha = rmatrixComponents%dalpha2
        dtheta = rmatrixComponents%dtheta2 
        diota  = rmatrixComponents%diota2
        dgamma = rmatrixComponents%dgamma2
        dnewton = rmatrixComponents%dnewton2
        iupwind = rmatrixComponents%iupwind2
        dupsam = rmatrixComponents%dupsam2
      end if
    
      ! Derive a temporary vector that contains only those velocity
      ! subvectors that might affect the matrix.
      call lsysbl_deriveSubvector(rvector,rtempVector, &
          imatOffset+1,imatOffset+2,.true.)
      call lsysbl_deriveSubvector(rdefect,rtempDefect, &
          imatOffset+1,imatOffset+2,.true.)
      
      ! Create a temporary block matrix only contining the velocity submatrices
      ! we want to change. Share structure and entries such that changing
      ! the temporary matrix will also change the original matrix.
      call lsysbl_deriveSubmatrix (rmatrix,rtempMatrix,&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,&
                                    imatOffset+1,imatOffset+2)

      ! Is A11=A22 physically?
      bshared = lsyssc_isMatrixContentShared(&
                    rtempmatrix%RmatrixBlock(1,1),&
                    rtempmatrix%RmatrixBlock(2,2))

      ! ---------------------------------------------------
      ! Subtract Ix ?
      if (diota*dcx .ne. 0.0_DP) then
        call lsyssc_vectorLinearComb (&
            rvector%RvectorBlock(imatOffset+1), &
            rdefect%RvectorBlock(imatOffset+1), &
            -diota*dcx, 1.0_DP)

        call lsyssc_vectorLinearComb (&
            rvector%RvectorBlock(imatOffset+2), &
            rdefect%RvectorBlock(imatOffset+2), &
            -diota*dcx, 1.0_DP)
      end if

      ! ---------------------------------------------------
      ! Subtract the mass matrix stuff?
      if (dalpha*dcx .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (rmatrixComponents%p_rmatrixMass, &
            rvector%RvectorBlock(imatOffset+1), &
            rdefect%RvectorBlock(imatOffset+1), &
            -dalpha*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (rmatrixComponents%p_rmatrixMass, &
            rvector%RvectorBlock(imatOffset+2), &
            rdefect%RvectorBlock(imatOffset+2), &
            -dalpha*dcx, 1.0_DP)
      end if
      
      ! ---------------------------------------------------
      ! Subtract the Stokes matrix stuff?
      if (dtheta*dcx .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (rmatrixComponents%p_rmatrixStokes, &
            rvector%RvectorBlock(imatOffset+1), &
            rdefect%RvectorBlock(imatOffset+1), &
            -dtheta*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (rmatrixComponents%p_rmatrixStokes, &
            rvector%RvectorBlock(imatOffset+2), &
            rdefect%RvectorBlock(imatOffset+2), &
            -dtheta*dcx, 1.0_DP)
      end if
      
      ! ---------------------------------------------------
      ! That was easy -- the adventure begins now... The nonlinearity!
      if (dgamma*dcx .ne. 0.0_DP) then
      
        select case (iupwind)
        case (0)
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rmatrixComponents%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion%dupsam = dupsam
          
          ! Matrix weight
          rstreamlineDiffusion%ddelta = dgamma*dcx
          
          ! Weight for the Newtop part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = dnewton*dcx
          
          ! Call the SD method to calculate the defect of the nonlinearity.
          ! As rrhsTemp shares its entries with rdefect, the result is
          ! directly written to rdefect!
          ! As velocity field, we specify rvelocityVector here -- the primal
          ! velocity!. The first two subvectors are used as velocity field.
          
          call conv_streamlineDiffusionBlk2d (&
                              rvelocityVector, rvelocityVector, &
                              dvectorWeight, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODDEFECT, &
                              rtempMatrix,rsolution=rtempVector,rdefect=rtempDefect)
                              
        case (1)
          ! Set up the upwind structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rupwind%dnu = rmatrixComponents%dnu
          
          ! Set stabilisation parameter
          rupwind%dupsam = dupsam

          ! Matrix weight
          rupwind%dtheta = dgamma*dcx
          
          ! Call the upwind method to calculate the nonlinear defect.
          call conv_upwind2d (rtempVector, rtempVector, &
                              dvectorWeight, 0.0_DP,&
                              rupwind, CONV_MODDEFECT, &
                              rtempMatrix%RmatrixBlock(1,1),&
                              rtempVector,rtempDefect) 
                              
          if (.not. bshared) then
            print *,'Upwind does not support independent A11/A22!'
            stop
          end if     

        case (2)
          ! Jump stabilisation.
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rmatrixComponents%dnu
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion%ddelta = dgamma*dcx
          
          ! Weight for the Newtop part; =0 deactivates Newton.
          rstreamlineDiffusion%dnewton = dnewton*dcx
          
          if (dnewton .eq. 0.0_DP) then

            ! Deactivate the matrices A12 and A21 by setting the multiplicators
            ! to 0.0. Whatever the content is (if there's content at all),
            ! these matrices are ignored then by the kernel.
            
            rtempMatrix%RmatrixBlock(1,2)%dscaleFactor = 0.0_DP
            rtempMatrix%RmatrixBlock(2,1)%dscaleFactor = 0.0_DP
            
          else

            ! Clear A12/A21 that receives parts of the Newton matrix
            call lsyssc_clearMatrix (rtempMatrix%RmatrixBlock(1,2))
            call lsyssc_clearMatrix (rtempMatrix%RmatrixBlock(2,1))
          
            ! Activate the submatrices A12 and A21 if they aren't.
            rtempMatrix%RmatrixBlock(1,2)%dscaleFactor = 1.0_DP
            rtempMatrix%RmatrixBlock(2,1)%dscaleFactor = 1.0_DP
           
          end if
         
          ! Call the SD method to calculate the nonlinearity.
          call conv_streamlineDiffusionBlk2d (&
                              rtempVector, rtempVector, &
                              dvectorWeight, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODDEFECT, &
                              rtempMatrix,rsolution=rtempVector,rdefect=rtempDefect)  
                              
          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rstreamlineDiffusion%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgammastar = dupsam
          rjumpStabil%dgamma = rjumpStabil%dgammastar
          
          ! Matrix weight
          rjumpStabil%dtheta = dgamma * dcx

          ! Call the jump stabilisation technique to stabilise that stuff.   
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
                              rjumpStabil, CONV_MODDEFECT, &
                              rtempMatrix%RmatrixBlock(1,1),&
                              rsolution=rtempVector,rdefect=rtempDefect)   

          if (.not. bshared) then
            print *,'Edge oriented stabilisation does not support independent A11/A22!'
            stop
          end if

        case DEFAULT
          print *,'Don''t know how to set up nonlinearity!?!'
          stop
        
        end select
      
      else
      
        ! That's the Stokes-case. Jump stabilisation is possible...
      
        select case (iupwind)
        case (2)
          ! Jump stabilisation.

          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rmatrixComponents%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgammastar = dupsam
          rjumpStabil%dgamma = rjumpStabil%dgammastar
          
          ! Matrix weight
          rjumpStabil%dtheta = dgamma * dcx

          ! Call the jump stabilisation technique to stabilise that stuff.   
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
                              rjumpStabil, CONV_MODDEFECT, &
                              rtempMatrix%RmatrixBlock(1,1),&
                              rsolution=rtempVector,rdefect=rtempDefect)   

          if (.not. bshared) then
            call conv_jumpStabilisation2d (&
                                rjumpStabil, CONV_MODDEFECT, &
                                rtempMatrix%RmatrixBlock(2,2),&
                                rsolution=rtempVector,rdefect=rtempDefect)   
          end if

        case DEFAULT
          ! No stabilisation
        
        end select
      
      end if ! gamma <> 0
      
      ! Release the temp matrix
      call lsysbl_releaseMatrix (rtempMatrix)
                          
      ! Release the temp vector if allocated.
      ! Derive a temporary vector that contains only those velocity
      ! subvectors that might affect the matrix.
      call lsysbl_releaseVector (rtempVector)
      call lsysbl_releaseVector (rtempDefect)
    
    end subroutine

    subroutine assemblePrimalReactMassDefect (rmatrixComponents,&
        rvector,rdefect,dcx,rvelocityVector,dvectorWeight)
        
    ! Assembles the defect arising from the reactive coupling mass
    ! matrices in the primal equation. rdefect must have been initialised 
    ! with the right hand side vector.
    !
    ! This special matrix is added in case that Newton is active. It has the
    ! form
    !        $$ M~ = dr_{11} N(\lambda) + dr_{12}   N*(\lambda) $$
    !
    ! the routine will construct
    !
    !       rdefect = r(primal)defect - dcx * (dtheta M~ r(dual)vector)
    
    ! A t_ccmatrixComponents structure providing all necessary 'source' information
    ! about how the mass matrix part is weighted (dr11, dr12).
    type(t_ccmatrixComponents), intent(IN) :: rmatrixComponents

    ! Solution vector. Provides the dual equation for the assembly.
    type(t_vectorBlock), intent(IN) :: rvector
    
    ! On entry: RHS vector.
    ! Is overwritten by the defect vector in the velocity subsystem.
    type(t_vectorBlock), intent(INOUT) :: rdefect
    
    ! Multiplication factor for the whole operator A*rvector
    real(DP), intent(IN) :: dcx
    
    ! Weight for the velocity vector rvelocityVector; usually = 1.0
    real(DP), intent(IN) :: dvectorWeight
    
    ! Velocity vector field that should be used for the assembly of the
    ! nonlinearity. Block 4 and 5 in that block vector are used as velocity 
    ! field.
    type(t_vectorBlock), intent(IN) :: rvelocityVector

      ! local variables
      type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
      type(t_matrixBlock) :: rtempmatrix
      type(t_vectorBlock) :: rtempvectorEval,rtempVectorDef
      
      ! If we have a reactive coupling mass matrix, it gets interesting...
      if ((rmatrixComponents%dr11 .ne. 0.0_DP) .or. &
          (rmatrixComponents%dr12 .ne. 0.0_DP)) then

        call lsysbl_createEmptyMatrix (rtempMatrix,2)

        ! Create a matrix with the structure we need. Share the structure
        ! of the mass matrix. Entries are not necessary for the assembly      
        call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
            rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

        call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
            rtempMatrix%RmatrixBlock(1,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

        call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
            rtempMatrix%RmatrixBlock(2,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

        call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
            rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
            
        ! Switch the matrices on.
        rtempMatrix%RmatrixBlock(1:2,1:2)%dscaleFactor = 1.0_DP
        
        call lsysbl_updateMatStrucInfo (rtempMatrix)
      
        ! The reactive part is: "dr1 * . * grad(lambda)".
        ! This is exactly the 'Newton' part assembled by the streamline diffusion
        ! method if we use the dual velocity as velocity field!
        ! So prepare to call streamline diffusion.

        ! Viscosity; ok, actually not used.
        rstreamlineDiffusion%dnu = rmatrixComponents%dnu
        
        ! Set stabilisation parameter
        rstreamlineDiffusion%dupsam = rmatrixComponents%dupsam2
        
        ! Weight dr21 of the convective part.
        rstreamlineDiffusion%ddelta = rmatrixComponents%dr11*dcx
        
        ! Weight for the Newton part; here, this is the dr22 weight.
        rstreamlineDiffusion%dnewton = rmatrixComponents%dr12*dcx
        
        ! Create a temporary block vector that points to the dual velocity.
        ! This has to be evaluated during the assembly.
        call lsysbl_deriveSubvector (rvelocityVector,rtempvectorEval,4,5,.true.)
        
        ! Create a temporary block vector for the dual defect.
        ! Matrix*primal velocity is subtracted from this.
        call lsysbl_deriveSubvector (rdefect,rtempvectorDef,1,2,.true.)
        
        ! Call the SD method to calculate the nonlinearity.
        ! As velocity vector, specify rtempvector, which points to
        ! the dual velocity.
        call conv_streamlineDiffusionBlk2d (&
                            rtempvectorEval, rtempvectorEval, &
                            dvectorWeight, 0.0_DP,&
                            rstreamlineDiffusion, CONV_MODDEFECT, &
                            rtempMatrix,rsolution=rvector,rdefect=rtempvectorDef)

        ! Release the temp matrix and temp vector.
        call lsysbl_releaseVector (rtempVectorEval)
        call lsysbl_releaseVector (rtempVectorDef)
        call lsysbl_releaseMatrix (rtempMatrix)
      
      end if
      
    end subroutine
      
    subroutine assembleDualReactMassDefect (rmatrixComponents,&
        rvector,rdefect,dcx,rvelocityVector,dvectorWeight)
        
    ! Assembles the defect arising from the reactive coupling mass
    ! matrices. rdefect must have been initialised with the right hand side 
    ! vector.
    !
    ! This special matrix is added in case that Newton is active. It has the
    ! form
    !        $$ M~ = dr_{21} N(\lambda) + dr_{22}   N*(\lambda) $$
    !
    ! the routine will construct
    !
    !       rdefect = r(dual)defect - dcx * (dtheta M~ r(primal)vector)
    
    ! A t_ccmatrixComponents structure providing all necessary 'source' information
    ! about how the mass matrix part is weighted (dr21, dr22).
    type(t_ccmatrixComponents), intent(IN) :: rmatrixComponents

    ! Solution vector. Provides the dual equation for the assembly.
    type(t_vectorBlock), intent(IN) :: rvector
    
    ! On entry: RHS vector.
    ! Is overwritten by the defect vector in the velocity subsystem.
    type(t_vectorBlock), intent(INOUT) :: rdefect
    
    ! Multiplication factor for the whole operator A*rvector
    real(DP), intent(IN) :: dcx
    
    ! Weight for the velocity vector rvelocityVector; usually = 1.0
    real(DP), intent(IN) :: dvectorWeight
    
    ! Velocity vector field that should be used for the assembly of the
    ! nonlinearity. The first two blocks in that block vector are
    ! used as velocity field.
    type(t_vectorBlock), intent(IN) :: rvelocityVector

      ! local variables
      type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
      type(t_matrixBlock) :: rtempmatrix
      type(t_vectorBlock) :: rtempvectorEval,rtempVectorDef
      
      ! If we have a reactive coupling mass matrix, it gets interesting...
      if ((rmatrixComponents%dr21 .ne. 0.0_DP) .or. &
          (rmatrixComponents%dr22 .ne. 0.0_DP)) then

        call lsysbl_createEmptyMatrix (rtempMatrix,2)

        ! Create a matrix with the structure we need. Share the structure
        ! of the mass matrix. Entries are not necessary for the assembly      
        call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
            rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

        call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
            rtempMatrix%RmatrixBlock(1,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

        call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
            rtempMatrix%RmatrixBlock(2,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

        call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
            rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
            
        ! Switch the matrices on.
        rtempMatrix%RmatrixBlock(1:2,1:2)%dscaleFactor = 1.0_DP
        
        call lsysbl_updateMatStrucInfo (rtempMatrix)
      
        ! The reactive part is: "dr2 * . * grad(lambda)".
        ! This is exactly the 'Newton' part assembled by the streamline diffusion
        ! method if we use the dual velocity as velocity field!
        ! So prepare to call streamline diffusion.

        ! Viscosity; ok, actually not used.
        rstreamlineDiffusion%dnu = rmatrixComponents%dnu
        
        ! Set stabilisation parameter
        rstreamlineDiffusion%dupsam = rmatrixComponents%dupsam2
        
        ! Weight dr21 of the convective part.
        rstreamlineDiffusion%ddelta = rmatrixComponents%dr21*dcx
        
        ! Weight for the Newton part; here, this is the dr22 weight.
        rstreamlineDiffusion%dnewton = rmatrixComponents%dr22*dcx
        
        ! Create a temporary block vector that points to the dual velocity.
        ! This has to be evaluated during the assembly.
        call lsysbl_deriveSubvector (rvelocityVector,rtempvectorEval,4,5,.true.)
        
        ! Create a temporary block vector for the dual defect.
        ! Matrix*primal velocity is subtracted from this.
        call lsysbl_deriveSubvector (rdefect,rtempvectorDef,4,5,.true.)
        
        ! Call the SD method to calculate the nonlinearity.
        ! As velocity vector, specify rtempvector, which points to
        ! the dual velocity.
        call conv_streamlineDiffusionBlk2d (&
                            rtempvectorEval, rtempvectorEval, &
                            dvectorWeight, 0.0_DP,&
                            rstreamlineDiffusion, CONV_MODDEFECT, &
                            rtempMatrix,rsolution=rvector,rdefect=rtempvectorDef)

        ! Release the temp matrix and temp vector.
        call lsysbl_releaseVector (rtempVectorEval)
        call lsysbl_releaseVector (rtempVectorDef)
        call lsysbl_releaseMatrix (rtempMatrix)
      
      end if
      
    end subroutine

    subroutine assemblePrimalUConstrMassDefect (rmatrixComponents,&
        rvector,rdefect,dcx,rvelocityVector)
        
    ! Assembles the defect arising from the projective coupling mass
    ! matrices in the primal equation which comes from constraints on u. 
    ! rdefect must have been initialised with the right hand side vector.
    !
    ! Let the mass matrix M~ be the usual mass matrix where u is ok
    ! and the 0-operator where u is out of bounds.
    ! Then, we assemble
    !
    !       rdefect = r(primal)defect - dcx (dmu1 M~ r(dual)vector)
    !
    ! A t_ccmatrixComponents structure providing all necessary 'source' information
    ! about how the mass matrix part is weighted (dr11, dr12).
    type(t_ccmatrixComponents), intent(IN) :: rmatrixComponents

    ! Solution vector. 
    type(t_vectorBlock), intent(IN) :: rvector
    
    ! On entry: RHS vector.
    ! Is overwritten by the defect vector in the velocity subsystem.
    type(t_vectorBlock), intent(INOUT) :: rdefect
    
    ! Multiplication factor for the whole operator A*rvector
    real(DP), intent(IN) :: dcx
    
    ! Velocity vector field that should be used for the assembly of the
    ! nonlinearity. Block 4 and 5 in that block vector are used as velocity 
    ! field.
    type(t_vectorBlock), intent(IN), target :: rvelocityVector

      ! local variables
      type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
      type(t_matrixBlock) :: rtempmatrix
      type(t_vectorBlock) :: rtempvectorEval,rtempVectorDef
      type(t_collection) :: rcollection
      type(t_bilinearForm) :: rform
      
      ! If we have a reactive coupling mass matrix, it gets interesting...
      if (rmatrixComponents%dmu1 .ne. 0.0_DP) then

        call lsysbl_createEmptyMatrix (rtempMatrix,2)

!        ! Create a matrix with the structure we need. Share the structure
!        ! of the mass matrix. Entries are not necessary for the assembly      
!        CALL lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
!            rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!        CALL lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
!            rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!        CALL lsysbl_updateMatStrucInfo (rtempMatrix)
!
!        ! Assemble the modified mass matrices.
!        
!        rform%itermCount = 1
!        rform%Idescriptors(1,1) = DER_FUNC
!        rform%Idescriptors(2,1) = DER_FUNC
!
!        ! In this case, we have nonconstant coefficients.
!        rform%ballCoeffConstant = .FALSE.
!        rform%BconstantCoeff(:) = .FALSE.
!
!        ! Prepare a collection structure to be passed to the callback
!        ! routine. We attach the vector T in the quick-access variables
!        ! so the callback routine can access it.
!        ! The bounds and the alpha value are passed in the
!        ! quickaccess-arrays.
!        call collct_init(rcollection)
!        rcollection%p_rvectorQuickAccess1 => rvelocityVector
!
!        ! Coefficient is dmu1=1/alpha or 0, depending on lambda
!        rcollection%DquickAccess(3)  = rmatrixComponents%dalphaC
!        rcollection%DquickAccess(4)  = rmatrixComponents%dmu1
!        
!        ! At first, set up A14, depending on lambda_1.
!        rcollection%IquickAccess(1) = 1
!        rcollection%DquickAccess(1) = rmatrixComponents%dumin1
!        rcollection%DquickAccess(2) = rmatrixComponents%dumax1
!        
!        ! Now we can build the matrix entries.
!        ! We specify the callback function coeff_Laplace for the coefficients.
!        ! As long as we use constant coefficients, this routine is not used.
!        ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
!        ! the framework will call the callback routine to get analytical
!        ! data.
!        ! The collection is passed as additional parameter. That's the way
!        ! how we get the vector to the callback routine.
!        call bilf_buildMatrixScalar (rform,.TRUE.,rtempmatrix%RmatrixBlock(1,1),&
!            coeff_ProjMass,rcollection)
!
!        ! Now, set up A22, depending on lambda_2.
!        rcollection%IquickAccess(1) = 2
!        rcollection%DquickAccess(1) = rmatrixComponents%dumin2
!        rcollection%DquickAccess(2) = rmatrixComponents%dumax2
!
!        call bilf_buildMatrixScalar (rform,.TRUE.,rtempmatrix%RmatrixBlock(2,2),&
!            coeff_ProjMass,rcollection)
!        
!        ! Now we can forget about the collection again.
!        call collct_done (rcollection)

        call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
            rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      
        call lsyssc_duplicateMatrix (rmatrixComponents%p_rmatrixMass,&
            rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
         
        call lsysbl_updateMatStrucInfo (rtempMatrix)
           
        ! Scale the entries by dmu1.
        if (rmatrixComponents%dmu1 .ne. 1.0_DP) then
          call lsyssc_scaleMatrix (rtempMatrix%RmatrixBlock(1,1),rmatrixComponents%dmu1)
          call lsyssc_scaleMatrix (rtempMatrix%RmatrixBlock(2,2),rmatrixComponents%dmu1)
        end if

        ! Filter the matrix. All the rows corresponding to DOF's that violate
        ! the bounds must be set to zero.
        call massmatfilter (rtempMatrix%RmatrixBlock(1,1),rvelocityVector%RvectorBlock(4),&
            rmatrixComponents%dalphaC,rmatrixComponents%dumin1,rmatrixComponents%dumax1)
        call massmatfilter (rtempMatrix%RmatrixBlock(2,2),rvelocityVector%RvectorBlock(5),&
            rmatrixComponents%dalphaC,rmatrixComponents%dumin2,rmatrixComponents%dumax2)
            

        !CALL lsysbl_deriveSubvector (rvector,rtempvectorEval,4,5,.false.)
        !CALL lsysbl_scaleVector (rtempvectorEval,-rmatrixComponents%dmu1)
        
        !Projection deactivated. For a defect, there is usually no projection!
        !call cc_projectControlTimestep (rtempvectorEval%RvectorBlock(1),&
        !    rmatrixComponents%dumin1,rmatrixComponents%dumax1)
        !call cc_projectControlTimestep (rtempvectorEval%RvectorBlock(2),&
        !    rmatrixComponents%dumin2,rmatrixComponents%dumax2)

        
        ! Create a temporary block vector that points to the dual velocity.
        ! This has to be evaluated during the assembly.
        call lsysbl_deriveSubvector (rvector,rtempvectorEval,4,5,.true.)
        
        ! Create a temporary block vector for the dual defect.
        ! Matrix*primal velocity is subtracted from this.
        call lsysbl_deriveSubvector (rdefect,rtempvectorDef,1,2,.true.)

        ! Create the defect
        ! Negative dcx if mu1 is directly incorporated to the matrix.
        ! Positive dcx if mu1 is incorporated to the vector rtempvectorEval.
        call lsysbl_blockMatVec (rtempmatrix, rtempvectorEval, rtempvectorDef, -dcx, 1.0_DP)
        !CALL lsysbl_blockMatVec (rtempmatrix, rtempvectorEval, rtempvectorDef, dcx, 1.0_DP)
        
        ! Release memory
        call lsysbl_releaseVector (rtempvectorDef)
        call lsysbl_releaseVector (rtempvectorEval)
        call lsysbl_releaseMatrix (rtempMatrix)
      
      end if
      
    end subroutine

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine cc_projectControlTimestep (rdualSolution,dumin,dumax)

!<description>
  ! Projects a dual solution vector u in such a way, that
  ! dumin <= u <= dumax holds.
!</description>

!<input>
  ! Minimum value for u
  real(DP), intent(in) :: dumin

  ! Maximum value for u
  real(DP), intent(in) :: dumax
!</input>

!<inputoutput>
  ! Vector to be restricted
  type(t_vectorScalar), intent(inout) :: rdualSolution
!</inputoutput>

!</subroutine>
 
    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: i
    real(dp) :: dactmin,dactmax
    
    ! Get the vector array
    call lsyssc_getbase_double (rdualSolution,p_Ddata)
    
    ! Restrict the vector
    do i=1,rdualSolution%NEQ
      p_Ddata(i) = min(max(p_Ddata(i),dumin),dumax)
    end do

  end subroutine   

end module
