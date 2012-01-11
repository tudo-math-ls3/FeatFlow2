!##############################################################################
!# ****************************************************************************
!# <name> optcontrolconvection </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# Contains extended assembly routines for the matrix in the optimal control
!# of distributed fluid flow.
!#
!# 1.) conv_strdiffOptC2dinitasm
!#     -> Initialises the assembly
!#
!# 2.) conv_strdiffOptC2ddoneasm
!#     -> Cleans up after an assembly
!#
!# 3.) conv_strdiffOptC2dgetMatrix
!#     -> Calculate the matrix of the operator.
!#
!# 4.) conv_strdiffOptC2dgetDerMatrix
!#     -> Calculate the matrix of the derivative of the operator.
!#
!# 5.) conv_strdiffOptC2dgetDefect
!#     -> Calculate the defect.
!#
!# </purpose>
!##############################################################################

module optcontrolconvection

  use fsystem
  use storage
  use genoutput
  use linearalgebra
  use derivatives
  use element
  use cubature
  use basicgeometry
  use dofmapping
  use triangulation
  use transformation
  use elementpreprocessing
  use spatialdiscretisation
  use bilinearformevaluation
  use linearsystemscalar
  use linearsystemblock
  use analyticsolution
  
  implicit none
  
  private
  
  public :: t_optcoperator
  public :: t_optcassemblyinfo
  public :: conv_strdiffOptC2dinitasm
  public :: conv_strdiffOptC2ddoneasm
  public :: conv_strdiffOptC2dgetMatrix
  public :: conv_strdiffOptC2dgetDerMatrix
  public :: conv_strdiffOptC2dgetDefect

  ! Defines the terms in the Optimal-Control operator to compute.

  type t_optcoperator
  
    ! Stabilisation parameter for the primal equation.
    ! Note: A value of 0.0_DP sets up the convection part without
    ! any stabilisation (central-difference like discretisation).
    real(DP) :: dupsamPrimal = 0.0_DP

    ! Stabilisation parameter for the dual equation.
    ! Note: A value of 0.0_DP sets up the convection part without
    ! any stabilisation (central-difference like discretisation).
    real(DP) :: dupsamDual = 0.0_DP
    
    ! Whether the viscosity is constant.
    logical :: bconstViscosity = .true.
    
    ! Viscosity parameter $\nu = 1/Re$ if viscosity is constant
    real(dp) :: dnu = 0.0_DP
    
    ! Coefficient in front of the Mass matrix in the primal equation
    real(dp) :: dprimalAlpha = 0.0_DP
  
    ! Coefficient in front of the Stokes operator nu*Laplace in the primal equation
    real(dp) :: dprimalBeta = 0.0_DP
    
    ! Coefficient in front of the convection operator \grad(.) y in the
    ! primal equation.
    real(dp) :: dprimalDelta = 0.0_DP

    ! Coefficient in front of the transposed convection operator \grad(.)^t y
    ! in the primal equation.
    real(dp) :: dprimalDeltaTrans = 0.0_DP
  
    ! Coefficient in front of the Newton operator \grad(y)(.) in the
    ! primal equation.
    real(dp) :: dprimalNewton = 0.0_DP
  
    ! Coefficient in front of the Newton operator \grad(y)^t(.) in the
    ! primal equation.
    real(dp) :: dprimalNewtonTrans = 0.0_DP

    ! Coefficient in front of the Mass matrix in the dual equation
    real(dp) :: ddualAlpha = 0.0_DP
  
    ! Coefficient in front of the Stokes operator nu*Laplace in the dual equation
    real(dp) :: ddualBeta = 0.0_DP
    
    ! Coefficient in front of the convection operator \grad(.) y in the
    ! dual equation.
    real(dp) :: ddualDelta = 0.0_DP

    ! Coefficient in front of the transposed convection operator \grad(.)^t y
    ! in the dual equation.
    real(dp) :: ddualDeltaTrans = 0.0_DP
  
    ! Coefficient in front of the Newton operator \grad(y)(.) in the
    ! dual equation.
    real(dp) :: ddualNewton = 0.0_DP
  
    ! Coefficient in front of the Newton operator \grad(y)^t(.) in the
    ! dual equation.
    real(dp) :: ddualNewtonTrans = 0.0_DP
    
    ! Coefficient in front of the transposed reactive convection operator
    ! \lambda\grad(.)^t in the dual equation
    real(dp) :: ddualRDeltaTrans = 0.0_DP
    
    ! Coefficient in front of the reactive Newton operator
    ! \grad(\lambda)(.) in the dual equation
    real(dp) :: ddualRNewton = 0.0_DP
    
    ! Coefficients in front of the reactive mass matrix that couples the dual
    ! to the primal equation
    real(dp) :: ddualRAlpha = 0.0_DP
    
    ! Multiplier to convert the dual velocity to the control.
    ! This is usually = -1/alpha for alpha coming from the optimal control problem.
    real(dp) :: dcontrolMultiplier = 0.0_DP
    
    ! Weight in front of the control u (= -1/alpha lambda).
    real(dp) :: dcontrolWeight = 0.0_DP
    
    ! Type of projection to use when converting the dual velocity to the
    ! control u.
    ! =0: No projection
    ! =1: DOF-based projection
    ! =2: Cubature-point based projection
    integer :: ccontrolProjection = 0
    
    ! Type of definition of the constraints if ccontrolConstraints <> 0.
    ! =0: constants specified in dumin1/2, dumax1/2.
    ! =1: analytical functions defined in p_rumin1/2, p_rumax1/2.
    integer :: cconstraintsType = 0

    ! Minimum/Maximum bound for the X-control.
    real(dp) :: dmin1 = -1E99_DP
    real(dp) :: dmax1 = 1E99_DP

    ! Analytical constraints for u_1.
    type(t_anSolution), pointer :: p_rumin1 => null()
    type(t_anSolution), pointer :: p_rumax1 => null()

    ! Minimum/Maximum bound for the Y-control.
    real(dp) :: dmin2 = -1E99_DP
    real(dp) :: dmax2 = 1E99_DP

    ! Analytical constraints for u_2
    type(t_anSolution), pointer :: p_rumin2 => null()
    type(t_anSolution), pointer :: p_rumax2 => null()
    
    ! Discrete constraints for u_1 and u_2.
    type(t_vectorBlock), pointer :: p_rvectorumin => null()
    type(t_vectorBlock), pointer :: p_rvectorumax => null()

    ! Calculation of local H.
    ! In 3D, there are 2 different methods to calculate the local H which
    ! is needed for the assembly of the convective part.
    ! =0: Use the cube-root of the volume of the hexahedron as local H
    ! =1: Use the length of the way that a particle travels through
    !     the hexahedron in direction of the flow
    integer :: clocalH = 1

  end type

  ! Assembly structure that saves information about the current
  ! status during the assembly of the local/global matrices.

  type t_optcassemblyinfo
  
    ! Block discretisation structure that defines the discretisation
    type(t_spatialDiscretisation), pointer :: p_rdiscretisation
    
    ! Id of the element distribution to be processed.
    integer :: ielemDistribution
    
    ! Pointer to the current element distribution
    type(t_elementDistribution), pointer :: p_relementDistribution

    ! Derivative specifier for evaluating elements
    logical, dimension(EL_MAXNDER) :: Bder
    
    ! An allocateable array accepting the DOF's of a set of elements.
    integer, dimension(:,:), pointer :: Idofs
    
    ! Allocateable arrays for the values of the basis functions -
    ! for test and trial spaces.
    real(DP), dimension(:,:,:,:), pointer :: Dbas

    ! An array receiving the coordinates of cubature points on
    ! the reference element for all elements in a set.
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef

    ! For every cubature point on the reference element,
    ! the corresponding cubature weight
    real(DP), dimension(CUB_MAXCUBP) :: Domega

    ! Type of transformation from the reference to the real element
    integer(I32) :: ctrafoType

    ! An element evaluation set for evaluating elements.
    type(t_evalElementSet) :: revalElementSet
    
    ! Specifier to signal when cubature points are initialised.
    logical :: bcubPtsInitialised
    
    ! Number of local DOF's
    integer :: indof
    
    ! Number of cubature points on the reference element
    integer :: ncubp
    
    ! Maximum number of simultaneously processed elements per element set
    integer :: nelementsPerBlock
    
    ! Local matrices, used during the assembly.
    ! Values and positions of values in the global matrix.
    integer, dimension(:,:,:), pointer :: Kentry
    
    ! Additional contributions for the submatrices A11, A12, A21, A22 stemming from Newton.
    integer, dimension(:,:,:), pointer :: Kentry12
    
    ! Maximum velocity magnitude
    real(DP) :: dumax
    
    ! An array with local DELTA's, each DELTA for one element, primal equation
    real(DP), dimension(:), pointer :: DlocalDeltaPrimal

    ! An array with local DELTA's, each DELTA for one element, dual equation
    real(DP), dimension(:), pointer :: DlocalDeltaDual
    
    ! A pointer to an element-number list
    integer, dimension(:), pointer :: p_IelementList

    ! Pointer to the primal/dual velocity field in the cubature points.
    real(DP), dimension(:,:,:), pointer :: Dpvel,Ddvel
    
    ! Pointer to the velocity X- and Y-derivative of the primal velocity
    ! in the cubature points
    real(DP), dimension(:,:,:), pointer :: DpvelXderiv,DpvelYderiv
    real(DP), dimension(:,:,:), pointer :: DdvelXderiv,DdvelYderiv
    
    ! Pointer to the primal/dual velocity DOF's on the elements where to
    ! evaluate the matrices.
    real(DP), dimension(:,:,:), pointer :: DpvelDofs,DdvelDofs

    ! Pointer to the primal/dual velocity DOF's on the elements which are
    ! multiplied to the local matrices to get the derivative of the
    ! operator.
    real(DP), dimension(:,:,:), pointer :: DpvelDofsAlt,DdvelDofsAlt
    
  end type

contains

!! ***************************************************************************
!
!  !<subroutine>
!
!  subroutine computeLocalMatricesDiag (Dbas,Domega,Ddetj,ndof,ncubp,NEL,&
!      Dvelocity,DvelocityXderiv,DvelocityYderiv,DlocalDelta,&
!      DentryA11,DentryA22,DentryA44,DentryA55,&
!      computeForm,dweight,dnu)
!
!  !<description>
!    ! Computes a local matrix to be incorporated into the global matrix.
!    ! Variant for handling the diagonal matrices.
!  !</description>
!
!  !<input>
!
!    ! Basis functions in all cubature points
!    real(DP), dimension(:,:,:,:), intent(in)  :: Dbas
!
!    ! Cubature weights
!    real(DP), dimension(:), intent(in)        :: Domega
!
!    ! Jacobian determinants in all cubature points
!    real(DP), dimension(:,:), intent(in)      :: Ddetj
!
!    ! Number of local DOF's in each element
!    integer, intent(in)                       :: ndof
!
!    ! Number of cubature points on each element
!    integer, intent(in)                       :: ncubp
!
!    ! Number of elements
!    integer, intent(in)                       :: NEL
!
!    ! Velocity vector and its derivatives in the cubature points
!    real(DP), dimension(:,:,:), intent(in)    :: Dvelocity
!    real(DP), dimension(:,:,:), intent(in)    :: DvelocityXderiv
!    real(DP), dimension(:,:,:), intent(in)    :: DvelocityYderiv
!
!    ! Local delta of the SD method
!    real(DP), dimension(:), intent(in)        :: DlocalDelta
!
!    ! Weight for A11-22, A44-55, A41-52, A14-25
!    real(DP), dimension(2,2), intent(in) :: Dweight
!
!    ! Viscosity parameter
!    real(DP), intent(in) :: dnu
!
!    interface
!      subroutine computeForm (dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY,&
!          du1,du1x,du1y,du2,du2x,du2y,dnu,dlocalDelta,Dweight,da11,da22,da44,da55)
!
!      use fsystem
!
!      !<description>
!        ! Compute the form in a cubature point.
!      !</description>
!
!      !<input>
!        ! Basis function phi_i and its derivatives.
!        real(DP), intent(in) :: dbasI,dbasIX,dbasIY
!
!        ! Basis function phi_j and its derivatives.
!        real(DP), intent(in) :: dbasJ,dbasJX,dbasJY
!
!        ! X-velocity of the given solution vector and its derivatives.
!        real(DP), intent(in) :: du1,du1x,du1y
!
!        ! Y-velocity of the given solution vector and its derivatives.
!        real(DP), intent(in) :: du2,du2x,du2y
!
!        ! Viscosity parameter
!        real(DP), intent(in) :: dnu
!
!        ! Local delta
!        real(DP), intent(in) :: dlocalDelta
!
!        ! Weight for A11-22, A44-55, A41-52, A14-25
!        real(DP), dimension(2,2), intent(in) :: Dweight
!      !</input>
!
!      !<output>
!        ! Values of the form
!        real(DP), intent(out) :: da11,da22,da44,da55
!      !</output>
!
!      end subroutine
!    end interface
!
!    external :: computeForm
!
!  !</input>
!
!  !<output>
!
!    ! Entries in the matrix
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA11
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA22
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA44
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA55
!
!  !</output>
!
!  !</subroutine>
!
!    ! local variables
!    integer :: iel, icubp, idofe, jdofe
!    real(DP) :: OM,du1,du2,du1locx,du1locy,du2locx,du2locY
!    real(DP) :: HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3
!    real(DP) :: AH11,AH22,AH44,AH55
!
!    AH11 = 0.0_DP
!    AH22 = 0.0_DP
!    AH44 = 0.0_DP
!    AH55 = 0.0_DP
!
!    ! Loop over the elements in the current set.
!    do iel=1,NEL
!
!      ! Loop over all cubature points on the current element
!      do icubp = 1, ncubp
!
!        ! Calculate the current weighting factor in the cubature formula
!        ! in that cubature point.
!        !
!        ! Normally, we have to take the absolut value of the determinant
!        ! of the mapping here!
!        ! In 2D, the determinant is always positive, whereas in 3D,
!        ! the determinant might be negative -- that's normal!
!        ! But because this routine only works in 2D, we can skip
!        ! the ABS here!
!
!        OM = Domega(icubp)*Ddetj(icubp,iel)
!
!        ! Current velocity in this cubature point:
!        du1 = Dvelocity (1,icubp,iel)
!        du2 = Dvelocity (2,icubp,iel)
!        du1locx = DvelocityXderiv (1,icubp,iel)
!        du1locy = DvelocityXderiv (2,icubp,iel)
!        du2locx = DvelocityYderiv (1,icubp,iel)
!        du2locy = DvelocityYderiv (2,icubp,iel)
!
!        ! Outer loop over the DOF's i=1..indof on our current element,
!        ! which corresponds to the basis functions Phi_i:
!
!        do idofe=1,ndof
!
!          ! Fetch the contributions of the (test) basis functions Phi_i
!          ! (our "O")  for function value and first derivatives for the
!          ! current DOF into HBASIy:
!
!          HBASI1 = Dbas(idofe,1,icubp,iel)
!          HBASI2 = Dbas(idofe,2,icubp,iel)
!          HBASI3 = Dbas(idofe,3,icubp,iel)
!
!          ! Inner loop over the DOF's j=1..indof, which corresponds to
!          ! the basis function Phi_j:
!
!          do jdofe=1,ndof
!
!            ! Fetch the contributions of the (trial) basis function Phi_j
!            ! (out "X") for function value and first derivatives for the
!            ! current DOF into HBASJy:
!
!            HBASJ1 = Dbas(jdofe,1,icubp,iel)
!            HBASJ2 = Dbas(jdofe,2,icubp,iel)
!            HBASJ3 = Dbas(jdofe,3,icubp,iel)
!
!            ! Finally calculate the contribution to the system
!            ! matrices A11, A12, A21 and A22.
!            call computeForm (HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3,&
!                du1,du2,du1locx,du1locy,du2locx,du2locy,dnu,DlocalDelta(iel),&
!                Dweight,AH11,AH22,AH44,AH55)
!
!            ! Weighten the calculated value AHxy by the cubature
!            ! weight OM and add it to the local matrices. After the
!            ! loop over all DOF's is finished, each entry contains
!            ! the calculated integral.
!
!            DentryA11(jdofe,idofe,iel) = DentryA11(jdofe,idofe,iel)+OM*AH11
!            DentryA22(jdofe,idofe,iel) = DentryA22(jdofe,idofe,iel)+OM*AH22
!            DentryA44(jdofe,idofe,iel) = DentryA44(jdofe,idofe,iel)+OM*AH44
!            DentryA55(jdofe,idofe,iel) = DentryA55(jdofe,idofe,iel)+OM*AH55
!
!          end do ! idofe
!
!        end do ! jdofe
!
!      end do ! icubp
!
!    end do ! iel
!
!  end subroutine
!
!  ! ***************************************************************************
!
!  !<subroutine>
!
!  subroutine computeLocalMatricesFullDiag (Dbas,Domega,Ddetj,ndof,ncubp,NEL,&
!      Dvelocity,DvelocityXderiv,DvelocityYderiv,DlocalDelta,&
!      DentryA11,DentryA22,DentryA44,DentryA55,&
!      DentryA12,DentryA21,DentryA45,DentryA54,&
!      computeForm,dweight,dnu)
!
!  !<description>
!    ! Computes a local matrix to be incorporated into the global matrix.
!    ! Variant for handling the diagonal matrix blocks.
!  !</description>
!
!  !<input>
!
!    ! Basis functions in all cubature points
!    real(DP), dimension(:,:,:,:), intent(in)  :: Dbas
!
!    ! Cubature weights
!    real(DP), dimension(:), intent(in)        :: Domega
!
!    ! Jacobian determinants in all cubature points
!    real(DP), dimension(:,:), intent(in)      :: Ddetj
!
!    ! Number of local DOF's in each element
!    integer, intent(in)                       :: ndof
!
!    ! Number of cubature points on each element
!    integer, intent(in)                       :: ncubp
!
!    ! Number of elements
!    integer, intent(in)                       :: NEL
!
!    ! Velocity vector and its derivatives in the cubature points
!    real(DP), dimension(:,:,:), intent(in)    :: Dvelocity
!    real(DP), dimension(:,:,:), intent(in)    :: DvelocityXderiv
!    real(DP), dimension(:,:,:), intent(in)    :: DvelocityYderiv
!
!    ! Local delta of the SD method
!    real(DP), dimension(:), intent(in)        :: DlocalDelta
!
!    ! Weight for A11-22, A44-55, A41-52, A14-25
!    real(DP), dimension(2,2), intent(in) :: Dweight
!
!    ! Viscosity parameter
!    real(DP), intent(in) :: dnu
!
!    interface
!      subroutine computeForm (dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY,&
!          du1,du1x,du1y,du2,du2x,du2y,dnu,dlocalDelta,Dweight,&
!          da11,da22,da44,da55,da12,da21,da45,da54)
!
!      use fsystem
!
!      !<description>
!        ! Compute the form in a cubature point.
!      !</description>
!
!      !<input>
!        ! Basis function phi_i and its derivatives.
!        real(DP), intent(in) :: dbasI,dbasIX,dbasIY
!
!        ! Basis function phi_j and its derivatives.
!        real(DP), intent(in) :: dbasJ,dbasJX,dbasJY
!
!        ! X-velocity of the given solution vector and its derivatives.
!        real(DP), intent(in) :: du1,du1x,du1y
!
!        ! Y-velocity of the given solution vector and its derivatives.
!        real(DP), intent(in) :: du2,du2x,du2y
!
!        ! Viscosity parameter
!        real(DP), intent(in) :: dnu
!
!        ! Local delta
!        real(DP), intent(in) :: dlocalDelta
!
!        ! Weight for A11-22, A44-55, A41-52, A14-25
!        real(DP), dimension(2,2), intent(in) :: Dweight
!      !</input>
!
!      !<output>
!        ! Values of the form
!        real(DP), intent(out) :: da11,da22,da44,da55,da12,da21,da45,da54
!      !</output>
!
!      end subroutine
!    end interface
!
!    external :: computeForm
!
!  !</input>
!
!  !<output>
!
!    ! Entries in the matrix
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA11
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA22
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA44
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA55
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA12
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA21
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA45
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA54
!
!  !</output>
!
!  !</subroutine>
!
!    ! local variables
!    integer :: iel, icubp, idofe, jdofe
!    real(DP) :: OM,du1,du2,du1locx,du1locy,du2locx,du2locY
!    real(DP) :: HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3
!    real(DP) :: AH11,AH22,AH44,AH55,AH12,AH21,AH45,AH54
!
!    AH11 = 0.0_DP
!    AH22 = 0.0_DP
!    AH44 = 0.0_DP
!    AH55 = 0.0_DP
!    AH12 = 0.0_DP
!    AH21 = 0.0_DP
!    AH45 = 0.0_DP
!    AH54 = 0.0_DP
!
!    ! Loop over the elements in the current set.
!    do iel=1,NEL
!
!      ! Loop over all cubature points on the current element
!      do icubp = 1, ncubp
!
!        ! Calculate the current weighting factor in the cubature formula
!        ! in that cubature point.
!        !
!        ! Normally, we have to take the absolut value of the determinant
!        ! of the mapping here!
!        ! In 2D, the determinant is always positive, whereas in 3D,
!        ! the determinant might be negative -- that's normal!
!        ! But because this routine only works in 2D, we can skip
!        ! the ABS here!
!
!        OM = Domega(icubp)*Ddetj(icubp,iel)
!
!        ! Current velocity in this cubature point:
!        du1 = Dvelocity (1,icubp,iel)
!        du2 = Dvelocity (2,icubp,iel)
!        du1locx = DvelocityXderiv (1,icubp,iel)
!        du1locy = DvelocityXderiv (2,icubp,iel)
!        du2locx = DvelocityYderiv (1,icubp,iel)
!        du2locy = DvelocityYderiv (2,icubp,iel)
!
!        ! Outer loop over the DOF's i=1..indof on our current element,
!        ! which corresponds to the basis functions Phi_i:
!
!        do idofe=1,ndof
!
!          ! Fetch the contributions of the (test) basis functions Phi_i
!          ! (our "O")  for function value and first derivatives for the
!          ! current DOF into HBASIy:
!
!          HBASI1 = Dbas(idofe,1,icubp,iel)
!          HBASI2 = Dbas(idofe,2,icubp,iel)
!          HBASI3 = Dbas(idofe,3,icubp,iel)
!
!          ! Inner loop over the DOF's j=1..indof, which corresponds to
!          ! the basis function Phi_j:
!
!          do jdofe=1,ndof
!
!            ! Fetch the contributions of the (trial) basis function Phi_j
!            ! (out "X") for function value and first derivatives for the
!            ! current DOF into HBASJy:
!
!            HBASJ1 = Dbas(jdofe,1,icubp,iel)
!            HBASJ2 = Dbas(jdofe,2,icubp,iel)
!            HBASJ3 = Dbas(jdofe,3,icubp,iel)
!
!            ! Finally calculate the contribution to the system matrices.
!            call computeForm (HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3,&
!                du1,du2,du1locx,du1locy,du2locx,du2locy,dnu,DlocalDelta(iel),&
!                Dweight,AH11,AH22,AH44,AH55,AH12,AH21,AH45,AH54)
!
!            ! Weighten the calculated value AHxy by the cubature
!            ! weight OM and add it to the local matrices. After the
!            ! loop over all DOF's is finished, each entry contains
!            ! the calculated integral.
!
!            DentryA11(jdofe,idofe,iel) = DentryA11(jdofe,idofe,iel)+OM*AH11
!            DentryA22(jdofe,idofe,iel) = DentryA22(jdofe,idofe,iel)+OM*AH22
!            DentryA44(jdofe,idofe,iel) = DentryA44(jdofe,idofe,iel)+OM*AH44
!            DentryA55(jdofe,idofe,iel) = DentryA55(jdofe,idofe,iel)+OM*AH55
!
!            DentryA12(jdofe,idofe,iel) = DentryA12(jdofe,idofe,iel)+OM*AH12
!            DentryA21(jdofe,idofe,iel) = DentryA21(jdofe,idofe,iel)+OM*AH21
!            DentryA45(jdofe,idofe,iel) = DentryA45(jdofe,idofe,iel)+OM*AH45
!            DentryA54(jdofe,idofe,iel) = DentryA54(jdofe,idofe,iel)+OM*AH54
!
!          end do ! idofe
!
!        end do ! jdofe
!
!      end do ! icubp
!
!    end do ! iel
!
!  end subroutine
!
!  ! ***************************************************************************
!
!  !<subroutine>
!
!  subroutine computeLocalMatricesFull (Dbas,Domega,Ddetj,ndof,ncubp,NEL,&
!      Dvelocity,DvelocityXderiv,DvelocityYderiv,DlocalDelta,&
!      DentryA11,DentryA22,DentryA44,DentryA55,&
!      DentryA12,DentryA21,DentryA45,DentryA54,&
!      DentryA41,DentryA52,DentryA42,DentryA51,&
!      DentryA14,DentryA25,DentryA24,DentryA15,&
!      computeForm,dweight,dnu)
!
!  !<description>
!    ! Computes a local matrix to be incorporated into the global matrix.
!    ! Variant for handling the full matrix.
!  !</description>
!
!  !<input>
!
!    ! Basis functions in all cubature points
!    real(DP), dimension(:,:,:,:), intent(in)  :: Dbas
!
!    ! Cubature weights
!    real(DP), dimension(:), intent(in)        :: Domega
!
!    ! Jacobian determinants in all cubature points
!    real(DP), dimension(:,:), intent(in)      :: Ddetj
!
!    ! Number of local DOF's in each element
!    integer, intent(in)                       :: ndof
!
!    ! Number of cubature points on each element
!    integer, intent(in)                       :: ncubp
!
!    ! Number of elements
!    integer, intent(in)                       :: NEL
!
!    ! Velocity vector and its derivatives in the cubature points
!    real(DP), dimension(:,:,:), intent(in)    :: Dvelocity
!    real(DP), dimension(:,:,:), intent(in)    :: DvelocityXderiv
!    real(DP), dimension(:,:,:), intent(in)    :: DvelocityYderiv
!
!    ! Local delta of the SD method
!    real(DP), dimension(:), intent(in)        :: DlocalDelta
!
!    ! Weight for A11-22, A44-55, A41-52, A14-25
!    real(DP), dimension(2,2), intent(in) :: Dweight
!
!    ! Viscosity parameter
!    real(DP), intent(in) :: dnu
!
!    interface
!      subroutine computeForm (dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY,&
!          du1,du1x,du1y,du2,du2x,du2y,dnu,dlocalDelta,Dweight,&
!          da11,da22,da44,da55,da12,da21,da45,da54,&
!          da41,da52,da42,da51,da14,da25,da24,da15)
!
!      use fsystem
!
!      !<description>
!        ! Compute the form in a cubature point.
!        ! Variant for handling the full matrix.
!      !</description>
!
!      !<input>
!        ! Basis function phi_i and its derivatives.
!        real(DP), intent(in) :: dbasI,dbasIX,dbasIY
!
!        ! Basis function phi_j and its derivatives.
!        real(DP), intent(in) :: dbasJ,dbasJX,dbasJY
!
!        ! X-velocity of the given solution vector and its derivatives.
!        real(DP), intent(in) :: du1,du1x,du1y
!
!        ! Y-velocity of the given solution vector and its derivatives.
!        real(DP), intent(in) :: du2,du2x,du2y
!
!        ! Viscosity parameter
!        real(DP), intent(in) :: dnu
!
!        ! Local delta
!        real(DP), intent(in) :: dlocalDelta
!
!        ! Weight for A11-22, A44-55, A41-52, A14-25
!        real(DP), dimension(2,2), intent(in) :: Dweight
!      !</input>
!
!      !<output>
!        ! Values of the form
!        real(DP), intent(out) :: da11,da22,da44,da55,da12,da21,da45,da54
!        real(DP), intent(out) :: da41,da52,da42,da51,da14,da25,da24,da15
!      !</output>
!
!      end subroutine
!    end interface
!
!    external :: computeForm
!
!  !</input>
!
!  !<output>
!
!    ! Entries in the matrix
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA11
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA22
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA44
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA55
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA12
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA21
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA45
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA54
!
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA41
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA52
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA42
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA51
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA14
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA25
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA24
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA15
!
!  !</output>
!
!  !<subroutine>
!
!    ! local variables
!    integer :: iel, icubp, idofe, jdofe
!    real(DP) :: OM,du1,du2,du1locx,du1locy,du2locx,du2locY
!    real(DP) :: HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3
!    real(DP) :: AH11,AH22,AH44,AH55,AH12,AH21,AH45,AH54
!    real(DP) :: AH41,AH52,AH42,AH51,AH14,AH25,AH24,AH15
!
!    AH11 = 0.0_DP
!    AH22 = 0.0_DP
!    AH44 = 0.0_DP
!    AH55 = 0.0_DP
!    AH12 = 0.0_DP
!    AH21 = 0.0_DP
!    AH45 = 0.0_DP
!    AH54 = 0.0_DP
!
!    AH41 = 0.0_DP
!    AH52 = 0.0_DP
!    AH42 = 0.0_DP
!    AH51 = 0.0_DP
!    AH14 = 0.0_DP
!    AH25 = 0.0_DP
!    AH24 = 0.0_DP
!    AH15 = 0.0_DP
!
!    ! Loop over the elements in the current set.
!    do iel=1,NEL
!
!      ! Loop over all cubature points on the current element
!      do icubp = 1, ncubp
!
!        ! Calculate the current weighting factor in the cubature formula
!        ! in that cubature point.
!        !
!        ! Normally, we have to take the absolut value of the determinant
!        ! of the mapping here!
!        ! In 2D, the determinant is always positive, whereas in 3D,
!        ! the determinant might be negative -- that's normal!
!        ! But because this routine only works in 2D, we can skip
!        ! the ABS here!
!
!        OM = Domega(icubp)*Ddetj(icubp,iel)
!
!        ! Current velocity in this cubature point:
!        du1 = Dvelocity (1,icubp,iel)
!        du2 = Dvelocity (2,icubp,iel)
!        du1locx = DvelocityXderiv (1,icubp,iel)
!        du1locy = DvelocityXderiv (2,icubp,iel)
!        du2locx = DvelocityYderiv (1,icubp,iel)
!        du2locy = DvelocityYderiv (2,icubp,iel)
!
!        ! Outer loop over the DOF's i=1..indof on our current element,
!        ! which corresponds to the basis functions Phi_i:
!
!        do idofe=1,ndof
!
!          ! Fetch the contributions of the (test) basis functions Phi_i
!          ! (our "O")  for function value and first derivatives for the
!          ! current DOF into HBASIy:
!
!          HBASI1 = Dbas(idofe,1,icubp,iel)
!          HBASI2 = Dbas(idofe,2,icubp,iel)
!          HBASI3 = Dbas(idofe,3,icubp,iel)
!
!          ! Inner loop over the DOF's j=1..indof, which corresponds to
!          ! the basis function Phi_j:
!
!          do jdofe=1,ndof
!
!            ! Fetch the contributions of the (trial) basis function Phi_j
!            ! (out "X") for function value and first derivatives for the
!            ! current DOF into HBASJy:
!
!            HBASJ1 = Dbas(jdofe,1,icubp,iel)
!            HBASJ2 = Dbas(jdofe,2,icubp,iel)
!            HBASJ3 = Dbas(jdofe,3,icubp,iel)
!
!            ! Finally calculate the contribution to the system matrices.
!            call computeForm (HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3,&
!                du1,du2,du1locx,du1locy,du2locx,du2locy,dnu,DlocalDelta(iel),&
!                Dweight,AH11,AH22,AH44,AH55,AH12,AH21,AH45,AH54,&
!                AH41,AH52,AH42,AH51,AH14,AH25,AH24,AH15)
!
!            ! Weighten the calculated value AHxy by the cubature
!            ! weight OM and add it to the local matrices. After the
!            ! loop over all DOF's is finished, each entry contains
!            ! the calculated integral.
!
!            DentryA11(jdofe,idofe,iel) = DentryA11(jdofe,idofe,iel)+OM*AH11
!            DentryA22(jdofe,idofe,iel) = DentryA22(jdofe,idofe,iel)+OM*AH22
!            DentryA44(jdofe,idofe,iel) = DentryA44(jdofe,idofe,iel)+OM*AH44
!            DentryA55(jdofe,idofe,iel) = DentryA55(jdofe,idofe,iel)+OM*AH55
!
!            DentryA12(jdofe,idofe,iel) = DentryA12(jdofe,idofe,iel)+OM*AH12
!            DentryA21(jdofe,idofe,iel) = DentryA21(jdofe,idofe,iel)+OM*AH21
!            DentryA45(jdofe,idofe,iel) = DentryA45(jdofe,idofe,iel)+OM*AH45
!            DentryA54(jdofe,idofe,iel) = DentryA54(jdofe,idofe,iel)+OM*AH54
!
!            DentryA14(jdofe,idofe,iel) = DentryA14(jdofe,idofe,iel)+OM*AH14
!            DentryA25(jdofe,idofe,iel) = DentryA25(jdofe,idofe,iel)+OM*AH25
!            DentryA24(jdofe,idofe,iel) = DentryA24(jdofe,idofe,iel)+OM*AH24
!            DentryA15(jdofe,idofe,iel) = DentryA15(jdofe,idofe,iel)+OM*AH15
!
!            DentryA41(jdofe,idofe,iel) = DentryA41(jdofe,idofe,iel)+OM*AH41
!            DentryA52(jdofe,idofe,iel) = DentryA52(jdofe,idofe,iel)+OM*AH52
!            DentryA51(jdofe,idofe,iel) = DentryA51(jdofe,idofe,iel)+OM*AH51
!            DentryA42(jdofe,idofe,iel) = DentryA42(jdofe,idofe,iel)+OM*AH42
!
!          end do ! idofe
!
!        end do ! jdofe
!
!      end do ! icubp
!
!    end do ! iel
!
!  end subroutine
!
!  subroutine computeFormMass (dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY,&
!      du1,du1x,du1y,du2,du2x,du2y,dnu,dlocalDelta,Dweight,da11,da22,da44,da55)
!
!  !<description>
!    ! Compute the form of the mass matrix.
!  !</description>
!
!  !<input>
!    ! Basis function phi_i and its derivatives.
!    real(DP), intent(in) :: dbasI,dbasIX,dbasIY
!
!    ! Basis function phi_j and its derivatives.
!    real(DP), intent(in) :: dbasJ,dbasJX,dbasJY
!
!    ! X-velocity of the given solution vector and its derivatives.
!    real(DP), intent(in) :: du1,du1x,du1y
!
!    ! Y-velocity of the given solution vector and its derivatives.
!    real(DP), intent(in) :: du2,du2x,du2y
!
!    ! Viscosity parameter
!    real(DP), intent(in) :: dnu
!
!    ! Local delta
!    real(DP), intent(in) :: dlocalDelta
!
!    ! Weight for A11-22, A44-55, A41-52, A14-25
!    real(DP), dimension(2,2), intent(in) :: Dweight
!  !</input>
!
!  !<output>
!    ! Values of the form
!    real(DP), intent(out) :: da11,da22,da44,da55
!  !</output>
!
!    real(dp) :: dtemp
!
!    ! dalpha*HBASI1*HBASJ1
!    dtemp = dbasI*dbasJ
!
!    da11 = Dweight(1,1)*dtemp
!    da22 = Dweight(1,1)*dtemp
!    da44 = Dweight(2,2)*dtemp
!    da55 = Dweight(2,2)*dtemp
!
!  end subroutine
!
!  subroutine computeFormStokes (dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY,&
!      du1,du1x,du1y,du2,du2x,du2y,dnu,dlocalDelta,Dweight,da11,da22,da44,da55)
!
!  !<description>
!    ! Compute the form of the Stokes matrix.
!  !</description>
!
!  !<input>
!    ! Basis function phi_i and its derivatives.
!    real(DP), intent(in) :: dbasI,dbasIX,dbasIY
!
!    ! Basis function phi_j and its derivatives.
!    real(DP), intent(in) :: dbasJ,dbasJX,dbasJY
!
!    ! X-velocity of the given solution vector and its derivatives.
!    real(DP), intent(in) :: du1,du1x,du1y
!
!    ! Y-velocity of the given solution vector and its derivatives.
!    real(DP), intent(in) :: du2,du2x,du2y
!
!    ! Viscosity parameter
!    real(DP), intent(in) :: dnu
!
!    ! Local delta
!    real(DP), intent(in) :: dlocalDelta
!
!    ! Weight for A11-22, A44-55, A41-52, A14-25
!    real(DP), dimension(2,2), intent(in) :: Dweight
!  !</input>
!
!  !<output>
!    ! Values of the form
!    real(DP), intent(out) :: da11,da22,da44,da55
!  !</output>
!
!    real(dp) :: dtemp
!
!    ! dny*(grad(phi_j,grad(phi_i))
!    dtemp = dnu*(dbasIX*dbasIX+dbasIY*dbasIY)
!
!    da11 = Dweight(1,1)*dtemp
!    da22 = Dweight(1,1)*dtemp
!    da44 = Dweight(2,2)*dtemp
!    da55 = Dweight(2,2)*dtemp
!
!  end subroutine
!
!  subroutine computeFormConvection (dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY,&
!      du1,du1x,du1y,du2,du2x,du2y,dnu,dlocalDelta,Dweight,da11,da22,da44,da55)
!
!  !<description>
!    ! Compute the form of the Convection matrix.
!  !</description>
!
!  !<input>
!    ! Basis function phi_i and its derivatives.
!    real(DP), intent(in) :: dbasI,dbasIX,dbasIY
!
!    ! Basis function phi_j and its derivatives.
!    real(DP), intent(in) :: dbasJ,dbasJX,dbasJY
!
!    ! X-velocity of the given solution vector and its derivatives.
!    real(DP), intent(in) :: du1,du1x,du1y
!
!    ! Y-velocity of the given solution vector and its derivatives.
!    real(DP), intent(in) :: du2,du2x,du2y
!
!    ! Viscosity parameter
!    real(DP), intent(in) :: dnu
!
!    ! Local delta
!    real(DP), intent(in) :: dlocalDelta
!
!    ! Weight for A11-22, A44-55, A41-52, A14-25
!    real(DP), dimension(2,2), intent(in) :: Dweight
!  !</input>
!
!  !<output>
!    ! Values of the form
!    real(DP), intent(out) :: da11,da22,da44,da55
!  !</output>
!
!    real(dp) :: HSUMI,HSUMJ,dtemp
!
!    ! Delta*(U*grad(Phi_j), U*grad(Phi_i)) + (U*grad(Phi_j),Phi_i)
!
!    HSUMI = dbasIX*du1 + dbasIY*du2
!    HSUMJ = dbasJX*du1 + dbasJY*du2
!
!    dtemp = HSUMJ*(dlocalDelta*HSUMI+dbasI)
!
!    da11 = Dweight(1,1)*dtemp
!    da22 = Dweight(1,1)*dtemp
!    da44 = Dweight(2,2)*dtemp
!    da55 = Dweight(2,2)*dtemp
!
!  end subroutine
!
!  subroutine computeFormNewton (dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY,&
!      du1,du1x,du1y,du2,du2x,du2y,dnu,dlocalDelta,Dweight,&
!      da11,da22,da44,da55,da12,da21,da45,da54)
!
!  !<description>
!    ! Compute the form in a cubature point.
!  !</description>
!
!  !<input>
!    ! Basis function phi_i and its derivatives.
!    real(DP), intent(in) :: dbasI,dbasIX,dbasIY
!
!    ! Basis function phi_j and its derivatives.
!    real(DP), intent(in) :: dbasJ,dbasJX,dbasJY
!
!    ! X-velocity of the given solution vector and its derivatives.
!    real(DP), intent(in) :: du1,du1x,du1y
!
!    ! Y-velocity of the given solution vector and its derivatives.
!    real(DP), intent(in) :: du2,du2x,du2y
!
!    ! Viscosity parameter
!    real(DP), intent(in) :: dnu
!
!    ! Local delta
!    real(DP), intent(in) :: dlocalDelta
!
!    ! Weight for A11-22, A44-55, A41-52, A14-25
!    real(DP), dimension(2,2), intent(in) :: Dweight
!  !</input>
!
!  !<output>
!    ! Values of the form
!    real(DP), intent(out) :: da11,da22,da44,da55,da12,da21,da45,da54
!  !</output>
!
!
!    real(dp) :: dtemp
!
!    dtemp = dbasI*dbasJ
!
!    da11 = Dweight(1,1) * du1x * dtemp
!    da12 = Dweight(1,1) * du1y * dtemp
!    da21 = Dweight(1,1) * du2x * dtemp
!    da22 = Dweight(1,1) * du2y * dtemp
!
!    da44 = Dweight(2,2) * du1x * dtemp
!    da45 = Dweight(2,2) * du1y * dtemp
!    da54 = Dweight(2,2) * du2x * dtemp
!    da55 = Dweight(2,2) * du2y * dtemp
!
!  end subroutine
!
!  subroutine computeFormConvectionTransposed (dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY,&
!      du1,du1x,du1y,du2,du2x,du2y,dnu,dlocalDelta,Dweight,&
!      da11,da22,da44,da55,da12,da21,da45,da54)
!
!  !<description>
!    ! Compute the form in a cubature point.
!  !</description>
!
!  !<input>
!    ! Basis function phi_i and its derivatives.
!    real(DP), intent(in) :: dbasI,dbasIX,dbasIY
!
!    ! Basis function phi_j and its derivatives.
!    real(DP), intent(in) :: dbasJ,dbasJX,dbasJY
!
!    ! X-velocity of the given solution vector and its derivatives.
!    real(DP), intent(in) :: du1,du1x,du1y
!
!    ! Y-velocity of the given solution vector and its derivatives.
!    real(DP), intent(in) :: du2,du2x,du2y
!
!    ! Viscosity parameter
!    real(DP), intent(in) :: dnu
!
!    ! Local delta
!    real(DP), intent(in) :: dlocalDelta
!
!    ! Weight for A11-22, A44-55, A41-52, A14-25
!    real(DP), dimension(2,2), intent(in) :: Dweight
!  !</input>
!
!  !<output>
!    ! Values of the form
!    real(DP), intent(out) :: da11,da22,da44,da55,da12,da21,da45,da54
!  !</output>
!
!
!    real(dp) :: dtemp1,dtemp2
!
!    dtemp1 = dbasJX*dbasI
!    dtemp2 = dbasJY*dbasI
!
!    da11 = Dweight(1,1) * du1 * dtemp1
!    da12 = Dweight(1,1) * du2 * dtemp1
!    da21 = Dweight(1,1) * du1 * dtemp2
!    da22 = Dweight(1,1) * du2 * dtemp2
!
!    da44 = Dweight(2,2) * du1 * dtemp1
!    da45 = Dweight(2,2) * du2 * dtemp1
!    da54 = Dweight(2,2) * du1 * dtemp2
!    da55 = Dweight(2,2) * du2 * dtemp2
!
!  end subroutine
!
!  subroutine computeFormNewtonTransposed (dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY,&
!      du1,du1x,du1y,du2,du2x,du2y,dnu,dlocalDelta,Dweight,&
!      da11,da22,da44,da55,da12,da21,da45,da54)
!
!  !<description>
!    ! Compute the form of the transposed Newton operator.
!  !</description>
!
!  !<input>
!    ! Basis function phi_i and its derivatives.
!    real(DP), intent(in) :: dbasI,dbasIX,dbasIY
!
!    ! Basis function phi_j and its derivatives.
!    real(DP), intent(in) :: dbasJ,dbasJX,dbasJY
!
!    ! X-velocity of the given solution vector and its derivatives.
!    real(DP), intent(in) :: du1,du1x,du1y
!
!    ! Y-velocity of the given solution vector and its derivatives.
!    real(DP), intent(in) :: du2,du2x,du2y
!
!    ! Viscosity parameter
!    real(DP), intent(in) :: dnu
!
!    ! Local delta
!    real(DP), intent(in) :: dlocalDelta
!
!    ! Weight for A11-22, A44-55, A41-52, A14-25
!    real(DP), dimension(2,2), intent(in) :: Dweight
!  !</input>
!
!  !<output>
!    ! Values of the form
!    real(DP), intent(out) :: da11,da22,da44,da55,da12,da21,da45,da54
!  !</output>
!
!
!    real(dp) :: dtemp
!
!    dtemp = dbasJ*dbasI
!
!    da11 = Dweight(1,1) * du1x * dtemp
!    da12 = Dweight(1,1) * du2x * dtemp
!    da21 = Dweight(1,1) * du1y * dtemp
!    da22 = Dweight(1,1) * du2y * dtemp
!
!    da44 = Dweight(2,2) * du1x * dtemp
!    da45 = Dweight(2,2) * du2x * dtemp
!    da54 = Dweight(2,2) * du1y * dtemp
!    da55 = Dweight(2,2) * du2y * dtemp
!
!  end subroutine
!
!!  es fehlt:
!!  * Umsetzung des diskreten Newtons
!
!  ! ***************************************************************************
!
!  !<subroutine>
!
!  subroutine computeLocalMatricesFullExt (Dbas,Domega,Ddetj,ndof,ncubp,NEL,&
!      Dpvel,DpvelXderiv,DpvelYderiv,Ddvel,DdvelXderiv,DdvelYderiv,&
!      DentryA11,DentryA22,DentryA44,DentryA55,&
!      DentryA12,DentryA21,DentryA45,DentryA54,&
!      DentryA41,DentryA52,DentryA42,DentryA51,&
!      DentryA14,DentryA25,DentryA24,DentryA15,&
!      DlocalDelta,roptcoperator)
!
!  !<description>
!    ! Computes a local matrix to be incorporated into the global matrix.
!    ! Variant for handling the full matrix.
!  !</description>
!
!  !<input>
!
!    ! Basis functions in all cubature points
!    real(DP), dimension(:,:,:,:), intent(in)  :: Dbas
!
!    ! Cubature weights
!    real(DP), dimension(:), intent(in)        :: Domega
!
!    ! Jacobian determinants in all cubature points
!    real(DP), dimension(:,:), intent(in)      :: Ddetj
!
!    ! Number of local DOF's in each element
!    integer, intent(in)                       :: ndof
!
!    ! Number of cubature points on each element
!    integer, intent(in)                       :: ncubp
!
!    ! Number of elements
!    integer, intent(in)                       :: NEL
!
!    ! Primal velocity vector and its derivatives in the cubature points
!    real(DP), dimension(:,:,:), intent(in)    :: Dpvel
!    real(DP), dimension(:,:,:), intent(in)    :: DpvelXderiv
!    real(DP), dimension(:,:,:), intent(in)    :: DpvelYderiv
!
!    real(DP), dimension(:,:,:), intent(in)    :: Ddvel
!    real(DP), dimension(:,:,:), intent(in)    :: DdvelXderiv
!    real(DP), dimension(:,:,:), intent(in)    :: DdvelYderiv
!
!    ! Local delta of the SD method
!    real(DP), dimension(:), intent(in)        :: DlocalDelta
!
!    ! Configuration of the operator
!    type(t_optcoperator), intent(in)          :: roptcoperator
!
!  !</input>
!
!  !<output>
!
!    ! Entries in the matrix
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA11
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA22
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA44
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA55
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA12
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA21
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA45
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA54
!
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA41
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA52
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA42
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA51
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA14
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA25
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA24
!    real(DP), dimension(:,:,:), intent(inout) :: DentryA15
!
!  !</output>
!
!  !<subroutine>
!
!    ! local variables
!    integer :: iel, icubp, idofe, jdofe
!    real(DP) :: OM
!    real(DP) :: du1p,du2p,du1locxp,du1locyp,du2locxp,du2locyp
!    real(DP) :: du1d,du2d,du1locxd,du1locyd,du2locxd,du2locyd
!    real(DP) :: HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3
!    real(DP) :: AH11,AH22,AH44,AH55,AH12,AH21,AH45,AH54
!    real(DP) :: AH41,AH52,AH42,AH51,AH14,AH25,AH24,AH15
!
!    AH11 = 0.0_DP
!    AH22 = 0.0_DP
!    AH44 = 0.0_DP
!    AH55 = 0.0_DP
!    AH12 = 0.0_DP
!    AH21 = 0.0_DP
!    AH45 = 0.0_DP
!    AH54 = 0.0_DP
!
!    AH41 = 0.0_DP
!    AH52 = 0.0_DP
!    AH42 = 0.0_DP
!    AH51 = 0.0_DP
!    AH14 = 0.0_DP
!    AH25 = 0.0_DP
!    AH24 = 0.0_DP
!    AH15 = 0.0_DP
!
!    ! Loop over the elements in the current set.
!    do iel=1,NEL
!
!      ! Loop over all cubature points on the current element
!      do icubp = 1, ncubp
!
!        ! Calculate the current weighting factor in the cubature formula
!        ! in that cubature point.
!        !
!        ! Normally, we have to take the absolut value of the determinant
!        ! of the mapping here!
!        ! In 2D, the determinant is always positive, whereas in 3D,
!        ! the determinant might be negative -- that's normal!
!        ! But because this routine only works in 2D, we can skip
!        ! the ABS here!
!
!        OM = Domega(icubp)*Ddetj(icubp,iel)
!
!        ! Current velocity in this cubature point:
!        du1p = Dpvel (1,icubp,iel)
!        du2p = Dpvel (2,icubp,iel)
!        du1locxp = DpvelXderiv (1,icubp,iel)
!        du1locyp = DpvelXderiv (2,icubp,iel)
!        du2locxp = DpvelYderiv (1,icubp,iel)
!        du2locyp = DpvelYderiv (2,icubp,iel)
!
!        du1d = Ddvel (1,icubp,iel)
!        du2d = Ddvel (2,icubp,iel)
!        du1locxd = DdvelXderiv (1,icubp,iel)
!        du1locyd = DdvelXderiv (2,icubp,iel)
!        du2locxd = DdvelYderiv (1,icubp,iel)
!        du2locyd = DdvelYderiv (2,icubp,iel)
!
!        ! Outer loop over the DOF's i=1..indof on our current element,
!        ! which corresponds to the basis functions Phi_i:
!
!        do idofe=1,ndof
!
!          ! Fetch the contributions of the (test) basis functions Phi_i
!          ! (our "O")  for function value and first derivatives for the
!          ! current DOF into HBASIy:
!
!          HBASI1 = Dbas(idofe,1,icubp,iel)
!          HBASI2 = Dbas(idofe,2,icubp,iel)
!          HBASI3 = Dbas(idofe,3,icubp,iel)
!
!          ! Inner loop over the DOF's j=1..indof, which corresponds to
!          ! the basis function Phi_j:
!
!          do jdofe=1,ndof
!
!            ! Fetch the contributions of the (trial) basis function Phi_j
!            ! (out "X") for function value and first derivatives for the
!            ! current DOF into HBASJy:
!
!            HBASJ1 = Dbas(jdofe,1,icubp,iel)
!            HBASJ2 = Dbas(jdofe,2,icubp,iel)
!            HBASJ3 = Dbas(jdofe,3,icubp,iel)
!
!            ! Finally calculate the contribution to the system matrices.
!            !call computeForm (HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3,&
!            !    du1,du2,du1locx,du1locy,du2locx,du2locy,dnu,DlocalDelta(iel),&
!            !    Dweight,AH11,AH22,AH44,AH55,AH12,AH21,AH45,AH54,&
!            !    AH41,AH52,AH42,AH51,AH14,AH25,AH24,AH15)
!
!            ! Weighten the calculated value AHxy by the cubature
!            ! weight OM and add it to the local matrices. After the
!            ! loop over all DOF's is finished, each entry contains
!            ! the calculated integral.
!
!            DentryA11(jdofe,idofe,iel) = DentryA11(jdofe,idofe,iel)+OM*AH11
!            DentryA22(jdofe,idofe,iel) = DentryA22(jdofe,idofe,iel)+OM*AH22
!            DentryA44(jdofe,idofe,iel) = DentryA44(jdofe,idofe,iel)+OM*AH44
!            DentryA55(jdofe,idofe,iel) = DentryA55(jdofe,idofe,iel)+OM*AH55
!
!            DentryA12(jdofe,idofe,iel) = DentryA12(jdofe,idofe,iel)+OM*AH12
!            DentryA21(jdofe,idofe,iel) = DentryA21(jdofe,idofe,iel)+OM*AH21
!            DentryA45(jdofe,idofe,iel) = DentryA45(jdofe,idofe,iel)+OM*AH45
!            DentryA54(jdofe,idofe,iel) = DentryA54(jdofe,idofe,iel)+OM*AH54
!
!            DentryA14(jdofe,idofe,iel) = DentryA14(jdofe,idofe,iel)+OM*AH14
!            DentryA25(jdofe,idofe,iel) = DentryA25(jdofe,idofe,iel)+OM*AH25
!            DentryA24(jdofe,idofe,iel) = DentryA24(jdofe,idofe,iel)+OM*AH24
!            DentryA15(jdofe,idofe,iel) = DentryA15(jdofe,idofe,iel)+OM*AH15
!
!            DentryA41(jdofe,idofe,iel) = DentryA41(jdofe,idofe,iel)+OM*AH41
!            DentryA52(jdofe,idofe,iel) = DentryA52(jdofe,idofe,iel)+OM*AH52
!            DentryA51(jdofe,idofe,iel) = DentryA51(jdofe,idofe,iel)+OM*AH51
!            DentryA42(jdofe,idofe,iel) = DentryA42(jdofe,idofe,iel)+OM*AH42
!
!          end do ! idofe
!
!        end do ! jdofe
!
!      end do ! icubp
!
!    end do ! iel
!
!  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine computeLocalOptCMatrices (Dbas,Domega,Ddetj,ndof,ncubp,&
      Ielements,DpvelDofs,DdvelDofs,&
      Dpvel,DpvelXderiv,DpvelYderiv,Ddvel,DdvelXderiv,DdvelYderiv,&
      DentryA11,DentryA22,DentryA44,DentryA55,&
      DentryA12,DentryA21,DentryA45,DentryA54,&
      DentryA41,DentryA52,DentryA42,DentryA51,&
      DentryA14,DentryA25,DentryA24,DentryA15,&
      DlocalDeltaPrimal,DlocalDeltaDual,roptcoperator,roptcassemblyinfo)
      
  !<description>
    ! Computes a local matrix of the optimal control operator
    ! which is to be incorporated into the global matrix.
  !</description>
      
  !<input>
  
    ! Basis functions in all cubature points
    real(DP), dimension(:,:,:,:), intent(in)  :: Dbas
    
    ! Cubature weights
    real(DP), dimension(:), intent(in)        :: Domega
    
    ! Jacobian determinants in all cubature points
    real(DP), dimension(:,:), intent(in)      :: Ddetj
    
    ! Number of local DOF's in each element
    integer, intent(in)                       :: ndof
    
    ! Number of cubature points on each element
    integer, intent(in)                       :: ncubp
    
    ! List of currently handled elements
    integer, dimension(:), intent(in)         :: Ielements
    
    ! DOF's in the primal/dual velocity
    real(DP), dimension(:,:,:), intent(in)    :: DpvelDofs
    real(DP), dimension(:,:,:), intent(in)    :: DdvelDofs

    ! Temporary arrays for the primal and dual velocity and their
    ! derivatives in the cubature points.
    real(DP), dimension(:,:,:), intent(inout)    :: Dpvel
    real(DP), dimension(:,:,:), intent(inout)    :: DpvelXderiv
    real(DP), dimension(:,:,:), intent(inout)    :: DpvelYderiv

    real(DP), dimension(:,:,:), intent(inout)    :: Ddvel
    real(DP), dimension(:,:,:), intent(inout)    :: DdvelXderiv
    real(DP), dimension(:,:,:), intent(inout)    :: DdvelYderiv
    
    ! Temporary space for local delta of the SD method, primal equation
    real(DP), dimension(:), intent(inout)        :: DlocalDeltaPrimal

    ! Temporary space for local delta of the SD method, dual equation
    real(DP), dimension(:), intent(inout)        :: DlocalDeltaDual
    
    ! Configuration of the operator
    type(t_optcoperator), intent(in)          :: roptcoperator
    
    ! Assembly information structure
    type(t_optcassemblyinfo), intent(in)      :: roptcassemblyinfo
    
  !</input>
    
  !<output>
    
    ! Entries in the matrix.
    ! All local matrices are saved transposed, i.e. we have
    ! DentryIJ(column,row,element). This has to be taken care of
    ! during the assembly and was introduced to gain higher speed
    ! as it avoids extensive jumps in the memory access.
    real(DP), dimension(:,:,:), intent(inout) :: DentryA11
    real(DP), dimension(:,:,:), intent(inout) :: DentryA22
    real(DP), dimension(:,:,:), intent(inout) :: DentryA44
    real(DP), dimension(:,:,:), intent(inout) :: DentryA55
    
    real(DP), dimension(:,:,:), intent(inout) :: DentryA12
    real(DP), dimension(:,:,:), intent(inout) :: DentryA21
    real(DP), dimension(:,:,:), intent(inout) :: DentryA45
    real(DP), dimension(:,:,:), intent(inout) :: DentryA54
    
    real(DP), dimension(:,:,:), intent(inout) :: DentryA41
    real(DP), dimension(:,:,:), intent(inout) :: DentryA52
    real(DP), dimension(:,:,:), intent(inout) :: DentryA42
    real(DP), dimension(:,:,:), intent(inout) :: DentryA51

    real(DP), dimension(:,:,:), intent(inout) :: DentryA14
    real(DP), dimension(:,:,:), intent(inout) :: DentryA25
    real(DP), dimension(:,:,:), intent(inout) :: DentryA24
    real(DP), dimension(:,:,:), intent(inout) :: DentryA15
 
  !</output>
  
  !<subroutine>
    
    ! local variables
    integer :: iel, icubp, idofe, jdofe, NEL
    real(DP) :: OM,dtemp,dsumI,dsumJ
    real(DP) :: du1p,du2p,du1d,du2d, db,dbx,dby
    real(DP) :: du1locxp,du1locyp,du2locxp,du2locyp
    real(DP) :: du1locxd,du1locyd,du2locxd,du2locyd
    real(DP) :: dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY
    real(DP) :: AH11,AH22,AH44,AH55,AH12,AH21,AH45,AH54
    real(DP) :: AH41,AH52,AH42,AH51,AH14,AH25,AH24,AH15
    real(dp) :: dtemp1,dtemp2
    
    NEL = size(Ielements)
    
    ! ===========================================
    ! Solution evaluation phase
    !
    ! Evaluate the solutions for the nonlinearity
    ! ===========================================

    ! Calculate local DELTA's for streamline diffusion method.
    ! (cf. p. 121 in Turek's CFD book).
    ! For every element, we need a local DELTA.
    ! Every local delta is weighted by the global "ddelta".
    ! If ddelta=0, we don't do anything as this disables the
    ! nonlinear term.
    ! If UPSAM=0.0, we have a central-difference like discretisation, which
    ! is one can see as the local stabilisation weight Delta is also = 0.0.
    ! In this case, we even switch of the calculation of the local Delta,
    ! as it is always =0.0, so we save a little bit time.
    if (roptcoperator%dupsamPrimal .ne. 0.0_DP) then
      call getLocalDeltaQuad (DpvelDofs,roptcassemblyinfo%p_rdiscretisation%p_rtriangulation,&
          Ielements,roptcassemblyinfo%dumax,roptcoperator%dupsamPrimal,&
          roptcoperator%dnu,DlocalDeltaPrimal)
    end if

    if (roptcoperator%dupsamDual .ne. 0.0_DP) then
      call getLocalDeltaQuad (DpvelDofs,roptcassemblyinfo%p_rdiscretisation%p_rtriangulation,&
          Ielements,roptcassemblyinfo%dumax,roptcoperator%dupsamDual,&
          roptcoperator%dnu,DlocalDeltaDual)
    end if
    
    ! Loop over all elements in the current set
    do iel=1,NEL
    
      ! Loop over all cubature points on the current element
      do icubp = 1, ncubp
      
        ! Primal/dual velocity
        du1p = 0.0_DP
        du2p = 0.0_DP
        du1d = 0.0_DP
        du2d = 0.0_DP
        
        ! X/Y-derivative of the primal/dual velocity
        du1locxp = 0.0_DP
        du1locyp = 0.0_DP
        du2locxp = 0.0_DP
        du2locyp = 0.0_DP

        du1locxd = 0.0_DP
        du1locyd = 0.0_DP
        du2locxd = 0.0_DP
        du2locyd = 0.0_DP
      
        ! Perform a loop through the trial DOF's.
        do jdofe=1,ndof

          ! Get the value of the (test) basis function
          ! phi_i (our "O") in the cubature point:
          
          db =  Dbas(jdofe,1,icubp,iel)
          dbx = Dbas(jdofe,DER_DERIV2D_X,icubp,iel)
          dby = Dbas(jdofe,DER_DERIV2D_Y,icubp,iel)
          
          ! Sum up to the value in the cubature point
          
          du1p = du1p + DpvelDofs(jdofe,iel,1)*db
          du2p = du2p + DpvelDofs(jdofe,iel,2)*db
          
          du1d = du1d + DdvelDofs(jdofe,iel,1)*db
          du2d = du2d + DdvelDofs(jdofe,iel,2)*db
          
          du1locxp = du1locxp + DpvelDofs(jdofe,iel,1)*dbx
          du1locyp = du1locyp + DpvelDofs(jdofe,iel,1)*dby
          du2locxp = du2locxp + DpvelDofs(jdofe,iel,2)*dbx
          du2locyp = du2locyp + DpvelDofs(jdofe,iel,2)*dby
                                                    
          du1locxd = du1locxd + DdvelDofs(jdofe,iel,1)*dbx
          du1locyd = du1locyd + DdvelDofs(jdofe,iel,1)*dby
          du2locxd = du2locxd + DdvelDofs(jdofe,iel,2)*dbx
          du2locyd = du2locyd + DdvelDofs(jdofe,iel,2)*dby

        end do ! jdofe
        
        ! Save the computed velocities and derivatives
        Dpvel(1,icubp,iel) = du1p
        Dpvel(2,icubp,iel) = du2p

        Ddvel(1,icubp,iel) = du1d
        Ddvel(2,icubp,iel) = du2d

        DpvelXderiv(1,icubp,iel) = du1locxp
        DpvelXderiv(2,icubp,iel) = du1locyp
        DpvelYderiv(1,icubp,iel) = du2locxp
        DpvelYderiv(2,icubp,iel) = du2locyp
      
        DdvelXderiv(1,icubp,iel) = du1locxd
        DdvelXderiv(2,icubp,iel) = du1locyd
        DdvelYderiv(1,icubp,iel) = du2locxd
        DdvelYderiv(2,icubp,iel) = du2locyd

      end do ! icubp
      
    end do ! iel
    
    ! ===========================================
    ! Matrix evaluation phase
    !
    ! Evaluate the matrices using the solutions
    ! from above.
    ! ===========================================
      
    AH11 = 0.0_DP
    AH22 = 0.0_DP
    AH44 = 0.0_DP
    AH55 = 0.0_DP
    
    AH12 = 0.0_DP
    AH21 = 0.0_DP
    AH45 = 0.0_DP
    AH54 = 0.0_DP
    
    AH41 = 0.0_DP
    AH52 = 0.0_DP
    AH42 = 0.0_DP
    AH51 = 0.0_DP
    
    AH14 = 0.0_DP
    AH25 = 0.0_DP
    AH24 = 0.0_DP
    AH15 = 0.0_DP
    
    ! Compute mass matrix part?
    if ((roptcoperator%dprimalAlpha .ne. 0.0_DP) .or. &
        (roptcoperator%ddualAlpha .ne. 0.0_DP) .or. &
        (roptcoperator%ddualRAlpha .ne. 0.0_DP)) then
    
      ! Loop over the elements in the current set.
      do iel=1,NEL

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! Calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Normally, we have to take the absolut value of the determinant
          ! of the mapping here!
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that's normal!
          ! But because this routine only works in 2D, we can skip
          ! the ABS here!

          OM = Domega(icubp)*Ddetj(icubp,iel)

          ! Outer loop over the DOF's i=1..indof on our current element,
          ! which corresponds to the basis functions Phi_i:

          do idofe=1,ndof
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! (our "O")  for function value and first derivatives for the
            ! current DOF into HBASIy:
          
            dbasI = Dbas(idofe,1,icubp,iel)
            
            ! Inner loop over the DOF's j=1..indof, which corresponds to
            ! the basis function Phi_j:

            do jdofe=1,ndof
              
              ! Fetch the contributions of the (trial) basis function Phi_j
              ! (out "X") for function value and first derivatives for the
              ! current DOF into HBASJy:
            
              dbasJ = Dbas(jdofe,1,icubp,iel)

              ! Contribution to the mass matrix -- weighted by the cubature weight
              dtemp = OM * dbasI * dbasJ

              ! Weighten the calculated value AHxy by the cubature
              ! weight OM and add it to the local matrices. After the
              ! loop over all DOF's is finished, each entry contains
              ! the calculated integral.

              DentryA11(jdofe,idofe,iel) = DentryA11(jdofe,idofe,iel) + &
                                           roptcoperator%dprimalAlpha * dtemp
                                           
              DentryA22(jdofe,idofe,iel) = DentryA22(jdofe,idofe,iel) + &
                                           roptcoperator%dprimalAlpha * dtemp
                                           
              DentryA44(jdofe,idofe,iel) = DentryA44(jdofe,idofe,iel) + &
                                           roptcoperator%ddualAlpha * dtemp
                                           
              DentryA55(jdofe,idofe,iel) = DentryA55(jdofe,idofe,iel) + &
                                           roptcoperator%ddualAlpha * dtemp

              DentryA41(jdofe,idofe,iel) = DentryA41(jdofe,idofe,iel) + &
                                           roptcoperator%ddualRAlpha * dtemp
                                           
              DentryA52(jdofe,idofe,iel) = DentryA52(jdofe,idofe,iel) + &
                                           roptcoperator%ddualRAlpha * dtemp
              
            end do ! idofe
            
          end do ! jdofe

        end do ! icubp
      
      end do ! iel
    end if
    
    ! Laplace matrix
    if ((roptcoperator%dprimalBeta .ne. 0.0_DP) .or. &
        (roptcoperator%ddualBeta .ne. 0.0_DP)) then
    
      ! Loop over the elements in the current set.
      do iel=1,NEL

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! Calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Normally, we have to take the absolut value of the determinant
          ! of the mapping here!
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that's normal!
          ! But because this routine only works in 2D, we can skip
          ! the ABS here!

          OM = Domega(icubp)*Ddetj(icubp,iel)

          ! Outer loop over the DOF's i=1..indof on our current element,
          ! which corresponds to the basis functions Phi_i:

          do idofe=1,ndof
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! (our "O")  for function value and first derivatives for the
            ! current DOF into HBASIy:
          
            dbasIX = Dbas(idofe,2,icubp,iel)
            dbasIY = Dbas(idofe,3,icubp,iel)
            
            ! Inner loop over the DOF's j=1..indof, which corresponds to
            ! the basis function Phi_j:

            do jdofe=1,ndof
              
              ! Fetch the contributions of the (trial) basis function Phi_j
              ! (out "X") for function value and first derivatives for the
              ! current DOF into HBASJy:
            
              dbasJX = Dbas(jdofe,2,icubp,iel)
              dbasJY = Dbas(jdofe,3,icubp,iel)

              ! dny*(grad(phi_j,grad(phi_i))
              dtemp = roptcoperator%dnu * (dbasIX*dbasJX + dbasIY*dbasJY)

              ! Weighten the calculated value AHxy by the cubature
              ! weight OM and add it to the local matrices. After the
              ! loop over all DOF's is finished, each entry contains
              ! the calculated integral.

              DentryA11(jdofe,idofe,iel) = DentryA11(jdofe,idofe,iel) + &
                  OM * roptcoperator%dprimalBeta * dtemp
                  
              DentryA22(jdofe,idofe,iel) = DentryA22(jdofe,idofe,iel) + &
                  OM * roptcoperator%dprimalBeta * dtemp
                  
              DentryA44(jdofe,idofe,iel) = DentryA44(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualBeta * dtemp
                  
              DentryA55(jdofe,idofe,iel) = DentryA55(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualBeta * dtemp
              
            end do ! idofe
            
          end do ! jdofe

        end do ! icubp
      
      end do ! iel
    end if

    ! Convection operator grad(.)*y.
    ! This operator implies the transposed Newton \grad(y)^t*(.) in the
    ! dual equation.
    if ((roptcoperator%dprimalDelta .ne. 0.0_DP) .or. &
        (roptcoperator%ddualDelta .ne. 0.0_DP)) then
    
      ! Loop over the elements in the current set.
      do iel=1,NEL

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! Calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Normally, we have to take the absolut value of the determinant
          ! of the mapping here!
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that's normal!
          ! But because this routine only works in 2D, we can skip
          ! the ABS here!

          OM = Domega(icubp)*Ddetj(icubp,iel)
          
          du1p = Dpvel(1,icubp,iel)
          du2p = Dpvel(2,icubp,iel)
          du1locxp = DpvelXderiv(1,icubp,iel)
          du1locyp = DpvelXderiv(2,icubp,iel)
          du2locxp = DpvelYderiv(1,icubp,iel)
          du2locyp = DpvelYderiv(2,icubp,iel)

          ! Outer loop over the DOF's i=1..indof on our current element,
          ! which corresponds to the basis functions Phi_i:

          do idofe=1,ndof
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! (our "O")  for function value and first derivatives for the
            ! current DOF into HBASIy:
          
            dbasI = Dbas(idofe,1,icubp,iel)
            dbasIX = Dbas(idofe,2,icubp,iel)
            dbasIY = Dbas(idofe,3,icubp,iel)
            
            dsumI = dbasIX*du1p + dbasIY*du2p
            
            ! Inner loop over the DOF's j=1..indof, which corresponds to
            ! the basis function Phi_j:

            do jdofe=1,ndof
              
              ! Fetch the contributions of the (trial) basis function Phi_j
              ! (out "X") for function value and first derivatives for the
              ! current DOF into HBASJy:
            
              dbasJ = Dbas(jdofe,1,icubp,iel)
              dbasJX = Dbas(jdofe,2,icubp,iel)
              dbasJY = Dbas(jdofe,3,icubp,iel)

              dsumJ = dbasJX*du1p + dbasJY*du2p

              ! Transposed Newton in the dual equation.
              dtemp = dbasJ*dbasI
            
              ! Delta*(U*grad(Phi_j), U*grad(Phi_i)) + (U*grad(Phi_j),Phi_i)
              dtemp1 = dsumJ * (DlocalDeltaPrimal(iel)*dsumI + dbasI)
              
              ! For the dual equation, the sign of the stabilisation term
              ! must be changed.
              ! Probably because
              !    ( -c grad(phi_j), -c grad(phi_i) ) = + ( c grad(phi_j), c grad(phi_i) )
              ! see the paper of Heinkenschloss!
              !
              ! HINT: Change reverted to '+'. The "-"-sign is now introduced
              ! by a negative dupsam parameter to compensate the minus
              ! sign of the operator.
              ! dtemp2 = dsumJ * (-DlocalDeltaDual(iel)*dsumI + dbasI)
              
              dtemp2 = dsumJ * (DlocalDeltaDual(iel)*dsumI + dbasI)

              ! Weighten the calculated value AHxy by the cubature
              ! weight OM and add it to the local matrices. After the
              ! loop over all DOF's is finished, each entry contains
              ! the calculated integral.

              DentryA11(jdofe,idofe,iel) = DentryA11(jdofe,idofe,iel) + &
                  OM * roptcoperator%dprimalDelta * dtemp1
                  
              DentryA22(jdofe,idofe,iel) = DentryA22(jdofe,idofe,iel) + &
                  OM * roptcoperator%dprimalDelta * dtemp1
                  
              DentryA44(jdofe,idofe,iel) = DentryA44(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualDelta * dtemp2
                  
              DentryA55(jdofe,idofe,iel) = DentryA55(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualDelta * dtemp2


              DentryA44(jdofe,idofe,iel) = DentryA44(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualNewtonTrans * du1locxp * dtemp
                  
              DentryA45(jdofe,idofe,iel) = DentryA45(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualNewtonTrans * du2locxp * dtemp
                  
              DentryA54(jdofe,idofe,iel) = DentryA54(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualNewtonTrans * du1locyp * dtemp
                  
              DentryA55(jdofe,idofe,iel) = DentryA55(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualNewtonTrans * du2locyp * dtemp
              
            end do ! idofe
            
          end do ! jdofe

        end do ! icubp
      
      end do ! iel
    end if
    
    ! Newton operator grad(y)*(.).
    ! This operator implies the reactive operator (.)*grad(\lambda)+(\grad(.)^t)\lambda
    ! in the dual equation.
    if (roptcoperator%dprimalNewton .ne. 0.0_DP) then
    
      ! Loop over the elements in the current set.
      do iel=1,NEL

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! Calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Normally, we have to take the absolut value of the determinant
          ! of the mapping here!
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that's normal!
          ! But because this routine only works in 2D, we can skip
          ! the ABS here!

          OM = Domega(icubp)*Ddetj(icubp,iel)
          
          du1locxp = DpvelXderiv(1,icubp,iel)
          du1locyp = DpvelXderiv(2,icubp,iel)
          du2locxp = DpvelYderiv(1,icubp,iel)
          du2locyp = DpvelYderiv(2,icubp,iel)
          
          du1locxd = DdvelXderiv(1,icubp,iel)
          du1locyd = DdvelXderiv(2,icubp,iel)
          du2locxd = DdvelYderiv(1,icubp,iel)
          du2locyd = DdvelYderiv(2,icubp,iel)
          
          du1d = Ddvel(1,icubp,iel)
          du2d = Ddvel(2,icubp,iel)

          ! Outer loop over the DOF's i=1..indof on our current element,
          ! which corresponds to the basis functions Phi_i:

          do idofe=1,ndof
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! (our "O")  for function value and first derivatives for the
            ! current DOF into HBASIy:
          
            dbasI = Dbas(idofe,1,icubp,iel)
            dbasIX = Dbas(idofe,2,icubp,iel)
            dbasIY = Dbas(idofe,3,icubp,iel)
            
            ! Inner loop over the DOF's j=1..indof, which corresponds to
            ! the basis function Phi_j:

            do jdofe=1,ndof
              
              ! Fetch the contributions of the (trial) basis function Phi_j
              ! (out "X") for function value and first derivatives for the
              ! current DOF into HBASJy:
            
              dbasJ = Dbas(jdofe,1,icubp,iel)
              dbasJX = Dbas(jdofe,2,icubp,iel)
              dbasJY = Dbas(jdofe,3,icubp,iel)

              ! Delta*(U*grad(Phi_j), U*grad(Phi_i)) + (U*grad(Phi_j),Phi_i)
              dtemp = dbasJ*dbasI
              dtemp1 = dbasJX*dbasI
              dtemp2 = dbasJY*dbasI

              ! Weighten the calculated value AHxy by the cubature
              ! weight OM and add it to the local matrices. After the
              ! loop over all DOF's is finished, each entry contains
              ! the calculated integral.

              ! Newton operator in the primal
              DentryA11(jdofe,idofe,iel) = DentryA11(jdofe,idofe,iel) + &
                  OM * roptcoperator%dprimalNewton * du1locxp * dtemp
                  
              DentryA12(jdofe,idofe,iel) = DentryA12(jdofe,idofe,iel) + &
                  OM * roptcoperator%dprimalNewton * du1locyp * dtemp
                  
              DentryA21(jdofe,idofe,iel) = DentryA21(jdofe,idofe,iel) + &
                  OM * roptcoperator%dprimalNewton * du2locxp * dtemp
                  
              DentryA22(jdofe,idofe,iel) = DentryA22(jdofe,idofe,iel) + &
                  OM * roptcoperator%dprimalNewton * du2locyp * dtemp


              ! In the dual equation, we have to set up the reactive
              ! operator.
              !
              ! Newton part \grad(lambda)*(.)
              DentryA41(jdofe,idofe,iel) = DentryA41(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualRNewton * du1locxd * dtemp
                  
              DentryA42(jdofe,idofe,iel) = DentryA42(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualRNewton * du1locyd * dtemp
                  
              DentryA51(jdofe,idofe,iel) = DentryA51(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualRNewton * du2locxd * dtemp
                  
              DentryA52(jdofe,idofe,iel) = DentryA52(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualRNewton * du2locyd * dtemp

              ! Transposed Convective part grad(lambda)^t*(.)
              DentryA41(jdofe,idofe,iel) = DentryA41(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualRDeltaTrans * du1d * dtemp1
                  
              DentryA42(jdofe,idofe,iel) = DentryA42(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualRDeltaTrans * du2d * dtemp1
                  
              DentryA51(jdofe,idofe,iel) = DentryA51(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualRDeltaTrans * du1d * dtemp2
                  
              DentryA52(jdofe,idofe,iel) = DentryA52(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualRDeltaTrans * du2d * dtemp2

            end do ! idofe
            
          end do ! jdofe

        end do ! icubp
      
      end do ! iel
    end if
    
    ! Compute the control mass matrix part?
    if (roptcoperator%dcontrolWeight .ne. 0.0_DP) then
    
      ! Here, we have two cases. In case this is a standard operator,
      ! we just set up a weighted mass matrix that converts the
      ! dual velocity to the control.
      
      if (roptcoperator%dprimalNewton .eq. 0.0_DP) then
    
        ! Loop over the elements in the current set.
        do iel=1,NEL

          ! Loop over all cubature points on the current element
          do icubp = 1, ncubp

            ! Calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Normally, we have to take the absolut value of the determinant
            ! of the mapping here!
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that's normal!
            ! But because this routine only works in 2D, we can skip
            ! the ABS here!

            OM = Domega(icubp)*Ddetj(icubp,iel)
            
            ! Outer loop over the DOF's i=1..indof on our current element,
            ! which corresponds to the basis functions Phi_i:

            do idofe=1,ndof
            
              ! Fetch the contributions of the (test) basis functions Phi_i
              ! (our "O")  for function value and first derivatives for the
              ! current DOF into HBASIy:
            
              dbasI = Dbas(idofe,1,icubp,iel)
              
              ! Inner loop over the DOF's j=1..indof, which corresponds to
              ! the basis function Phi_j:

              do jdofe=1,ndof
                
                ! Fetch the contributions of the (trial) basis function Phi_j
                ! (out "X") for function value and first derivatives for the
                ! current DOF into HBASJy:
              
                dbasJ = Dbas(jdofe,1,icubp,iel)

                ! Contribution to the mass matrix
                dtemp = dbasI * dbasJ

                ! Weighten the calculated value AHxy by the cubature
                ! weight OM and add it to the local matrices. After the
                ! loop over all DOF's is finished, each entry contains
                ! the calculated integral.

                DentryA14(jdofe,idofe,iel) = DentryA14(jdofe,idofe,iel) + &
                     OM * roptcoperator%dcontrolWeight * roptcoperator%dcontrolMultiplier * dtemp
                                             
                DentryA25(jdofe,idofe,iel) = DentryA25(jdofe,idofe,iel) + &
                     OM * roptcoperator%dcontrolWeight * roptcoperator%dcontrolMultiplier * dtemp
                
              end do ! idofe
              
            end do ! jdofe

          end do ! icubp
        
        end do ! iel
      
      else
      
        ! In the other case, we have to set up a projected mass matrix.
        ! The question is, which projection to use.
        
        select case (roptcoperator%ccontrolProjection)
        case (0)
          ! No projection. Same loop as above.
          !
          ! Loop over the elements in the current set.
          do iel=1,NEL

            ! Loop over all cubature points on the current element
            do icubp = 1, ncubp

              ! Calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Normally, we have to take the absolut value of the determinant
              ! of the mapping here!
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that's normal!
              ! But because this routine only works in 2D, we can skip
              ! the ABS here!

              OM = Domega(icubp)*Ddetj(icubp,iel)

              ! Outer loop over the DOF's i=1..indof on our current element,
              ! which corresponds to the basis functions Phi_i:

              do idofe=1,ndof
              
                ! Fetch the contributions of the (test) basis functions Phi_i
                ! (our "O")  for function value and first derivatives for the
                ! current DOF into HBASIy:
              
                dbasI = Dbas(idofe,1,icubp,iel)
                
                ! Inner loop over the DOF's j=1..indof, which corresponds to
                ! the basis function Phi_j:

                do jdofe=1,ndof
                  
                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! (out "X") for function value and first derivatives for the
                  ! current DOF into HBASJy:
                
                  dbasJ = Dbas(jdofe,1,icubp,iel)

                  ! Contribution to the mass matrix
                  dtemp = roptcoperator%dcontrolWeight * roptcoperator%dcontrolMultiplier * &
                      dbasI * dbasJ

                  ! Weighten the calculated value AHxy by the cubature
                  ! weight OM and add it to the local matrices. After the
                  ! loop over all DOF's is finished, each entry contains
                  ! the calculated integral.

                  DentryA14(jdofe,idofe,iel) = DentryA14(jdofe,idofe,iel) + OM * dtemp
                  DentryA25(jdofe,idofe,iel) = DentryA25(jdofe,idofe,iel) + OM * dtemp
                  
                end do ! idofe
                
              end do ! jdofe

            end do ! icubp
          
          end do ! iel
         
        case (1)
          ! DOF-based projection. At first, set up the standard matrix as above.
          !
          ! Loop over the elements in the current set.
          do iel=1,NEL

            ! Loop over all cubature points on the current element
            do icubp = 1, ncubp

              ! Calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Normally, we have to take the absolut value of the determinant
              ! of the mapping here!
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that's normal!
              ! But because this routine only works in 2D, we can skip
              ! the ABS here!

              OM = Domega(icubp)*Ddetj(icubp,iel)

              ! Outer loop over the DOF's i=1..indof on our current element,
              ! which corresponds to the basis functions Phi_i:

              do idofe=1,ndof
              
                ! Fetch the contributions of the (test) basis functions Phi_i
                ! (our "O")  for function value and first derivatives for the
                ! current DOF into HBASIy:
              
                dbasI = Dbas(idofe,1,icubp,iel)
                
                ! Check if the DOF violates the bounds. If that's the case,
                ! it must be ignored. At first, create the control from the dual velocity
                dtemp1 = DdvelDofs(idofe,iel,1) * roptcoperator%dcontrolMultiplier
                dtemp2 = DdvelDofs(idofe,iel,2) * roptcoperator%dcontrolMultiplier
                
                ! Inner loop over the DOF's j=1..indof, which corresponds to
                ! the basis function Phi_j:
                if ((dtemp1 .ge. roptcoperator%dmin1) .and.  &
                    (dtemp1 .le. roptcoperator%dmax1)) then
                  do jdofe=1,ndof
                    
                    ! Fetch the contributions of the (trial) basis function Phi_j
                    ! (out "X") for function value and first derivatives for the
                    ! current DOF into HBASJy:
                  
                    dbasJ = Dbas(jdofe,1,icubp,iel)

                    ! Contribution to the mass matrix
                    dtemp = roptcoperator%dcontrolWeight * roptcoperator%dcontrolMultiplier * &
                        dbasI * dbasJ

                    ! Weighten the calculated value AHxy by the cubature
                    ! weight OM and add it to the local matrices. After the
                    ! loop over all DOF's is finished, each entry contains
                    ! the calculated integral.

                    DentryA14(jdofe,idofe,iel) = DentryA14(jdofe,idofe,iel) + OM * dtemp
                    
                  end do ! idofe
                end if

                if ((dtemp2 .ge. roptcoperator%dmin2) .and.  &
                    (dtemp2 .le. roptcoperator%dmax2)) then
                  do jdofe=1,ndof
                    
                    ! Fetch the contributions of the (trial) basis function Phi_j
                    ! (out "X") for function value and first derivatives for the
                    ! current DOF into HBASJy:
                  
                    dbasJ = Dbas(jdofe,1,icubp,iel)

                    ! Contribution to the mass matrix
                    dtemp = roptcoperator%dcontrolWeight * roptcoperator%dcontrolMultiplier * &
                        dbasI * dbasJ

                    ! Weighten the calculated value AHxy by the cubature
                    ! weight OM and add it to the local matrices. After the
                    ! loop over all DOF's is finished, each entry contains
                    ! the calculated integral.

                    DentryA25(jdofe,idofe,iel) = DentryA25(jdofe,idofe,iel) + OM * dtemp
                    
                  end do ! idofe
                end if
                
              end do ! jdofe

            end do ! icubp
          
          end do ! iel
          
        case (2)
          ! Cubature-point based projection.
          !
          ! Loop over the elements in the current set.
          do iel=1,NEL

            ! Loop over all cubature points on the current element
            do icubp = 1, ncubp

              ! Calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Normally, we have to take the absolut value of the determinant
              ! of the mapping here!
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that's normal!
              ! But because this routine only works in 2D, we can skip
              ! the ABS here!

              OM = Domega(icubp)*Ddetj(icubp,iel)

              du1d = roptcoperator%dcontrolMultiplier * Ddvel(1,icubp,iel)
              du2d = roptcoperator%dcontrolMultiplier * Ddvel(2,icubp,iel)

              ! Outer loop over the DOF's i=1..indof on our current element,
              ! which corresponds to the basis functions Phi_i:

              do idofe=1,ndof
              
                ! Fetch the contributions of the (test) basis functions Phi_i
                ! (our "O")  for function value and first derivatives for the
                ! current DOF into HBASIy:
              
                dbasI = Dbas(idofe,1,icubp,iel)
                
                ! Inner loop over the DOF's j=1..indof, which corresponds to
                ! the basis function Phi_j:

                do jdofe=1,ndof
                  
                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! (out "X") for function value and first derivatives for the
                  ! current DOF into HBASJy:
                
                  dbasJ = Dbas(jdofe,1,icubp,iel)

                  ! Contribution to the mass matrix
                  dtemp = roptcoperator%dcontrolWeight * roptcoperator%dcontrolMultiplier * &
                      dbasI * dbasJ
                      
                  ! Calculate the control in the cubature point
                  dtemp1 = dtemp
                  dtemp2 = dtemp
                  
                  ! Compute dtemp1/dtemp2. Set it to dtemp if the control is
                  ! in the bounds, otherwise to 0.
                  if ((du1d .lt. roptcoperator%dmin1) .or. (du1d .gt. roptcoperator%dmax1)) then
                    dtemp1 = 0.0_DP
                  end if

                  if ((du2d .lt. roptcoperator%dmin2) .or. (du2d .gt. roptcoperator%dmax2)) then
                    dtemp2 = 0.0_DP
                  end if

                  ! Weighten the calculated value AHxy by the cubature
                  ! weight OM and add it to the local matrices. After the
                  ! loop over all DOF's is finished, each entry contains
                  ! the calculated integral.

                  DentryA14(jdofe,idofe,iel) = DentryA14(jdofe,idofe,iel) + OM * dtemp1
                  DentryA25(jdofe,idofe,iel) = DentryA25(jdofe,idofe,iel) + OM * dtemp2
                  
                end do ! idofe
                
              end do ! jdofe

            end do ! icubp
          
          end do ! iel
         
        end select
      
      end if
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine conv_strdiffOptC2dinitasm (rdiscretisation,ielemDistribution,&
      roptcassemblyinfo,balternativeVelocity)
  
!<description>
  ! Initialise the roptcassemblyinfo structure for the assembly of
  ! local matrices.
!</description>

!<input>
  ! Block discretisation structure that defines the discretisation
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisation
  
  ! Id of the element distribution to be processed.
  integer, intent(in) :: ielemDistribution
  
  ! OPTIONAL: Create arrays for an additional velocity, independent of the
  ! evaluation point of the matrices. If not specified, FALSE is assumed.
  logical, intent(in), optional :: balternativeVelocity
!</input>
  
!<output>
  ! Structure collecting information about the assembly.
  type(t_optcassemblyinfo), intent(out) :: roptcassemblyinfo
!</output>

!</subroutine>

    ! local variables
    type(t_elementDistribution), pointer :: p_relementDistribution
    integer :: i,k
    
    ! Cubature point coordinates on the reference element
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi

    ! Get the element distribution. Here, we can find information about
    ! the cubature formula etc...
    p_relementDistribution => rdiscretisation%RelementDistr(ielemDistribution)
    
    ! Remember the discretisation and the element distribution
    roptcassemblyinfo%p_rdiscretisation => rdiscretisation
    roptcassemblyinfo%p_relementDistribution => p_relementDistribution
    roptcassemblyinfo%ielemDistribution = ielemDistribution
    
    ! Get the number of local DOF's for trial/test functions.
    ! We assume trial and test functions to be the same.
    roptcassemblyinfo%indof = elem_igetNDofLoc(p_relementDistribution%celement)

    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    roptcassemblyinfo%nelementsPerBlock = min(BILF_NELEMSIM,rdiscretisation%p_rtriangulation%NEL)
    
    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    roptcassemblyinfo%ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)
    
    ! Allocate some memory to hold the cubature points on the reference element
    allocate(roptcassemblyinfo%p_DcubPtsRef(&
        trafo_igetReferenceDimension(roptcassemblyinfo%ctrafoType),CUB_MAXCUBP))

    ! Initialise the cubature formula,
    ! Get cubature weights and point coordinates on the reference element
    call cub_getCubPoints(roptcassemblyinfo%p_relementDistribution%ccubTypeBilForm, &
        roptcassemblyinfo%ncubp, Dxi, roptcassemblyinfo%Domega)
    
    ! Reformat the cubature points; they are in the wrong shape!
    do i=1,roptcassemblyinfo%ncubp
      do k=1,ubound(roptcassemblyinfo%p_DcubPtsRef,1)
        roptcassemblyinfo%p_DcubPtsRef(k,i) = Dxi(i,k)
      end do
    end do
    
    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  allocate(Dbas(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly,
    ! which reduces the speed by 50%!
    
    allocate(roptcassemblyinfo%Dbas(&
        roptcassemblyinfo%indof,elem_getMaxDerivative(p_relementDistribution%celement), &
        roptcassemblyinfo%ncubp,roptcassemblyinfo%nelementsPerBlock))
        
    ! Allocate memory for the DOF's of all the elements.
    allocate(roptcassemblyinfo%Idofs(roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    
    ! Allocate memory for array with local DELTA's
    allocate(roptcassemblyinfo%DlocalDeltaPrimal(roptcassemblyinfo%nelementsPerBlock))
    allocate(roptcassemblyinfo%DlocalDeltaDual(roptcassemblyinfo%nelementsPerBlock))

    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*BILF_NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*BILF_NELEMSIM is normally much smaller!
    !
    ! Kentry (:,:,:) defines the positions of the local matrices
    ! in the submatrices A11 and A22.
    ! KentryA12 (:,:,:) defines the positions of the local matrices
    ! in the submatrices A12 and A21.
    allocate(roptcassemblyinfo%Kentry(&
        roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    
    allocate(roptcassemblyinfo%Kentry12(&
        roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    
    ! Allocate memory for the velocites in the cubature points.
    allocate(roptcassemblyinfo%Dpvel(NDIM2D,roptcassemblyinfo%ncubp,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(roptcassemblyinfo%Ddvel(NDIM2D,roptcassemblyinfo%ncubp,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(roptcassemblyinfo%DpvelXderiv(NDIM2D,roptcassemblyinfo%ncubp,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(roptcassemblyinfo%DpvelYderiv(NDIM2D,roptcassemblyinfo%ncubp,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(roptcassemblyinfo%DdvelXderiv(NDIM2D,roptcassemblyinfo%ncubp,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(roptcassemblyinfo%DdvelYderiv(NDIM2D,roptcassemblyinfo%ncubp,&
        roptcassemblyinfo%nelementsPerBlock))
    
    ! Allocate memory for the DOF's on the elements.
    allocate (roptcassemblyinfo%DpvelDofs(&
        roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock,NDIM2D))
    allocate (roptcassemblyinfo%DdvelDofs(&
        roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock,NDIM2D))
       
    ! The alternative velocity coincides with the standard velocity.
    roptcassemblyinfo%DpvelDofsAlt => roptcassemblyinfo%DpvelDofs
    roptcassemblyinfo%DdvelDofsAlt => roptcassemblyinfo%DdvelDofs

    if (present(balternativeVelocity)) then
      if (balternativeVelocity) then
        ! Create an additional evaluation point for the velocity.
        allocate (roptcassemblyinfo%DpvelDofsAlt(&
            roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock,NDIM2D))
        allocate (roptcassemblyinfo%DdvelDofsAlt(&
            roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock,NDIM2D))
      end if
    end if
    
    ! Indicate that cubature points must still be initialised in the element set.
    roptcassemblyinfo%bcubPtsInitialised = .false.
    
    ! Initialise the derivative flags
    roptcassemblyinfo%Bder = .false.
    roptcassemblyinfo%Bder(DER_FUNC) = .true.
    roptcassemblyinfo%Bder(DER_DERIV2D_X) = .true.
    roptcassemblyinfo%Bder(DER_DERIV2D_Y) = .true.
    
    ! Initialisation of the element set.
    call elprep_init(roptcassemblyinfo%revalElementSet)

    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    call storage_getbase_int (p_relementDistribution%h_IelementList, &
                              roptcassemblyinfo%p_IelementList)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine conv_strdiffOptC2dinitelemset (rvelmatrix,rvelmatrixoffdiag,&
      rvelocityVectorPrimal,rvelocityVectorDual,&
      roptcoperator,roptcassemblyinfo,istartElement,iendElement,&
      rvelocityVectorPrimalAlt,rvelocityVectorDualAlt)
  
!<description>
  ! Initialise the roptcassemblyinfo structure for the assembly of
  ! local matrices for the elements istartElement..istartElement+ielCountMax-1 (of the
  ! current element set).
!</description>

!<input>
  ! Template FE matrix specifying the connectivity of the primal/dual
  ! velocity space for diagonal blocks.
  type(t_matrixScalar), intent(in) :: rvelmatrix

  ! Template FE matrix specifying the connectivity of the primal/dual
  ! velocity space for offdiagonal blocks.
  type(t_matrixScalar), intent(in) :: rvelmatrixoffdiag
  
  ! Velocity vector, primal equation (Block 1/2) for the evaluation of the matrices
  type(t_vectorBlock), intent(in) :: rvelocityVectorPrimal

  ! Velocity vector, dual equation (Blcok 4/5) for the evaluation of the matrices
  type(t_vectorBlock), intent(in) :: rvelocityVectorDual

  ! Structure defining the operator to set up.
  type(t_optcoperator), intent(in) :: roptcoperator

  ! Index of the start element of the current element set in roptcassemblyinfo.
  integer, intent(in) :: istartElement
  
  ! Index of the last element.
  integer, intent(in) :: iendElement

  ! OPTIONAL: Alternative Velocity vector, primal equation (Block 1/2).
  type(t_vectorBlock), intent(in), optional :: rvelocityVectorPrimalAlt

  ! OPTIONAL: Alternative Velocity vector, dual equation (Blcok 4/5).
  type(t_vectorBlock), intent(in), optional :: rvelocityVectorDualAlt

!</input>

!<inputoutput>
  ! Structure collecting information about the assembly.
  type(t_optcassemblyinfo), intent(inout) :: roptcassemblyinfo
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel, idofe, jdofe, jcol0, jcol, jdfg
    integer, dimension(:), pointer :: p_Kcol
    integer, dimension(:), pointer :: p_Kld
    integer, dimension(:), pointer :: p_Kcol12
    integer, dimension(:), pointer :: p_Kld12

    ! Primal and dual velocity
    real(DP), dimension(:), pointer :: p_DpvelocityX, p_DpvelocityY
    real(DP), dimension(:), pointer :: p_DdvelocityX, p_DdvelocityY
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer(I32) :: cevaluationTag
    
    ! Get pointers to the matrix structure(s).
    call lsyssc_getbase_Kcol (rvelmatrix,p_Kcol)
    call lsyssc_getbase_Kld (rvelmatrix,p_Kld)
    
    call lsyssc_getbase_Kcol (rvelmatrixoffdiag,p_Kcol12)
    call lsyssc_getbase_Kld (rvelmatrixoffdiag,p_Kld12)
   
    ! The outstanding feature with finite elements is: A basis
    ! function for a DOF on one element has common support only
    ! with the DOF's on the same element! E.g. for Q1:
    !
    !        #. . .#. . .#. . .#
    !        .     .     .     .
    !        .  *  .  *  .  *  .
    !        #-----O-----O. . .#
    !        |     |     |     .
    !        |     | iel |  *  .
    !        #-----X-----O. . .#
    !        |     |     |     .
    !        |     |     |  *  .
    !        #-----#-----#. . .#
    !
    ! --> On element iel, the basis function at "X" only interacts
    !     with the basis functions in "O". Elements in the
    !     neighbourhood ("*") have no support, therefore we only have
    !     to collect all "O" DOF's.
    !
    ! Calculate the global DOF's into IdofsTrial / IdofsTest.
    !
    ! More exactly, we call dof_locGlobMapping_mult to calculate all the
    ! global DOF's of our BILF_NELEMSIM elements simultaneously.
    call dof_locGlobMapping_mult(roptcassemblyinfo%p_rdiscretisation, &
        roptcassemblyinfo%p_IelementList(istartElement:iendElement), roptcassemblyinfo%Idofs)
                                
    ! For the assembly of the global matrix, we use a "local"
    ! approach. At first we build a "local" system matrix according
    ! to the current element. This contains all additive
    ! contributions of element iel, which are later added at the
    ! right positions to the elements in the global system matrix.
    !
    ! We have indofTrial trial DOF's per element and
    ! indofTest test DOF's per element. Therefore there are
    ! indofTrial*indofTest tupel of basis-/testfunctions (phi_i,psi_j)
    ! "active" (i.e. have common support) on our current element, each
    ! giving an additive contribution to the system matrix.
    !
    ! We build a quadratic indofTrial*indofTest local matrix:
    ! Kentry(1..indofTrial,1..indofTest) receives the position
    !   in the global system matrix, where the corresponding value
    !   has to be added to.
    ! (The corresponding contrbutions can be saved separately,
    !  but we directly add them to the global matrix in this
    !  approach.)
    !
    ! We build local matrices for all our elements
    ! in the set simultaneously.
    ! Loop through elements in the set and for each element,
    ! loop through the local matrices to initialise them:
    do iel=1,iendelement-istartElement+1
    
      ! For building the local matrices, we have first to
      ! loop through the test functions (the "O"'s), as these
      ! define the rows in the matrix.
      do idofe=1,roptcassemblyinfo%indof
      
        ! Row idofe of the local matrix corresponds
        ! to row=global DOF KDFG(idofe) in the global matrix.
        ! This is one of the the "O"'s in the above picture.
        ! Get the starting position of the corresponding row
        ! to jcol0:

        jcol0=p_KLD(roptcassemblyinfo%Idofs(idofe,iel))
        
        ! Now we loop through the other DOF's on the current element
        ! (the "O"'s).
        ! All these have common support with our current basis function
        ! and will therefore give an additive value to the global
        ! matrix.
        
        do jdofe=1,roptcassemblyinfo%indof
          
          ! Get the global DOF of the "X" which interacts with
          ! our "O".
          
          jdfg=roptcassemblyinfo%Idofs(jdofe,iel)
          
          ! Starting in jcol0 (which points to the beginning of
          ! the line initially), loop through the elements in
          ! the row to find the position of column IDFG.
          ! Jump out of the do loop if we find the column.
          
          do jcol=jcol0,rvelmatrix%NA
            if (p_KCOL(jcol) .eq. jdfg) exit
          end do
          
          ! Because columns in the global matrix are sorted
          ! ascendingly (except for the diagonal element),
          ! the next search can start after the column we just found.
          
          ! jcol0=jcol+1
          
          ! Save the position of the matrix entry into the local
          ! matrix.
          ! Note that a column in Kentry corresponds to a row in
          ! the real matrix. We aligned Kentry/DENTRY this way to get
          ! higher speed of the assembly routine, since this leads
          ! to better data locality.
          
          roptcassemblyinfo%Kentry(jdofe,idofe,iel)=jcol
          
        end do ! idofe
        
      end do ! jdofe
      
    end do ! iel
    
    ! If the Newton part is to be calculated, we also need the matrix positions
    ! in A12 and A21. We can skip this part if the column structure is
    ! exactly the same!
    if (associated(p_Kcol,p_Kcol12)) then
    
      roptcassemblyinfo%Kentry12(:,:,:) = roptcassemblyinfo%Kentry(:,:,:)
      
    else

      do iel=1,iendelement-istartElement+1
      
        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"'s), as these
        ! define the rows in the matrix.
        do idofe=1,roptcassemblyinfo%indof
        
          ! Row idofe of the local matrix corresponds
          ! to row=global DOF KDFG(idofe) in the global matrix.
          ! This is one of the the "O"'s in the above picture.
          ! Get the starting position of the corresponding row
          ! to jcol0:

          jcol0=p_KLD12(roptcassemblyinfo%Idofs(idofe,iel))
          
          ! Now we loop through the other DOF's on the current element
          ! (the "O"'s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.
          
          do jdofe=1,roptcassemblyinfo%indof
            
            ! Get the global DOF of the "X" which interacts with
            ! our "O".
            
            jdfg=roptcassemblyinfo%Idofs(jdofe,iel)
            
            ! Starting in jcol0 (which points to the beginning of
            ! the line initially), loop through the elements in
            ! the row to find the position of column IDFG.
            ! Jump out of the do loop if we find the column.
            
            do jcol=jcol0,rvelmatrixoffdiag%NA
              if (p_KCOL12(jcol) .eq. jdfg) exit
            end do

            ! Because columns in the global matrix are sorted
            ! ascendingly (except for the diagonal element),
            ! the next search can start after the column we just found.
            
            ! jcol0=jcol+1
            
            ! Save the position of the matrix entry into the local
            ! matrix.
            ! Note that a column in Kentry corresponds to a row in
            ! the real matrix. We aligned Kentry/DENTRY this way to get
            ! higher speed of the assembly routine, since this leads
            ! to better data locality.
            
            roptcassemblyinfo%Kentry12(jdofe,idofe,iel)=jcol
            
          end do ! idofe
          
        end do ! jdofe
        
      end do ! iel
      
    end if
    
    ! Ok, we found the positions of the local matrix entries
    ! that we have to change.
    ! To calculate the matrix contributions, we have to evaluate
    ! the elements to give us the values of the basis functions
    ! in all the DOF's in all the elements in our set.

    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag.
    cevaluationTag = elem_getEvaluationTag(roptcassemblyinfo%p_relementDistribution%celement)
                    
    ! In the first loop, calculate the coordinates on the reference element.
    ! In all later loops, use the precalculated information.
    !
    ! Note: Why not using
    !   if (ielset .EQ. 1) then
    ! here, but this strange concept with the boolean variable?
    ! Because the if-command does not work with OpenMP! bcubPtsInitialised
    ! is a local variable and will therefore ensure that every thread
    ! is initialising its local set of cubature points!
    if (.not. roptcassemblyinfo%bcubPtsInitialised) then
      roptcassemblyinfo%bcubPtsInitialised = .true.
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)
    else
      cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
    end if
    
    ! Calculate all information that is necessary to evaluate the finite element
    ! on all cells of our subset. This includes the coordinates of the points
    ! on the cells.
    call elprep_prepareSetForEvaluation (roptcassemblyinfo%revalElementSet,&
        cevaluationTag, roptcassemblyinfo%p_rdiscretisation%p_rtriangulation, &
        roptcassemblyinfo%p_IelementList(istartElement:iendelement), &
        roptcassemblyinfo%ctrafoType, roptcassemblyinfo%p_DcubPtsRef(:,1:roptcassemblyinfo%ncubp))

    ! Calculate the values of the basis functions.
    ! Pass p_DcubPts as point coordinates, which point either to the
    ! coordinates on the reference element (the same for all elements)
    ! or on the real element - depending on whether this is a
    ! parametric or nonparametric element.
    call elem_generic_sim2 (roptcassemblyinfo%p_relementDistribution%celement, &
        roptcassemblyinfo%revalElementSet, roptcassemblyinfo%Bder, roptcassemblyinfo%Dbas)
          
    ! Calculate the primal and dual velocity and its derivatives
    ! in all the cubature points on all the elements.
    ! The calculation routine might need them.
    
    ! Get pointers to the velocity to be evaluated.
    call lsyssc_getbase_double (rvelocityVectorPrimal%RvectorBlock(1),p_DpvelocityX)
    call lsyssc_getbase_double (rvelocityVectorPrimal%RvectorBlock(2),p_DpvelocityY)
    call lsyssc_getbase_double (rvelocityVectorDual%RvectorBlock(4),p_DdvelocityX)
    call lsyssc_getbase_double (rvelocityVectorDual%RvectorBlock(5),p_DdvelocityY)
    
    ! Extract the velocity DOF's omn the elements.
    do iel=1,iendelement-istartElement+1
      do idofe = 1,roptcassemblyinfo%indof
        jdofe = roptcassemblyinfo%Idofs(idofe,iel)
        roptcassemblyinfo%DpvelDofs(idofe,iel,1) = p_DpvelocityX(jdofe)
        roptcassemblyinfo%DpvelDofs(idofe,iel,2) = p_DpvelocityY(jdofe)
        roptcassemblyinfo%DdvelDofs(idofe,iel,1) = p_DdvelocityX(jdofe)
        roptcassemblyinfo%DdvelDofs(idofe,iel,2) = p_DdvelocityY(jdofe)
      end do
    end do

    if (present(rvelocityVectorPrimalAlt)) then
      ! Also fetch the alternative velocity.
      call lsyssc_getbase_double (rvelocityVectorPrimalAlt%RvectorBlock(1),p_DpvelocityX)
      call lsyssc_getbase_double (rvelocityVectorPrimalAlt%RvectorBlock(2),p_DpvelocityY)
      call lsyssc_getbase_double (rvelocityVectorDualAlt%RvectorBlock(4),p_DdvelocityX)
      call lsyssc_getbase_double (rvelocityVectorDualAlt%RvectorBlock(5),p_DdvelocityY)

      ! Extract the velocity DOF's omn the elements.
      do iel=1,iendelement-istartElement+1
        do idofe = 1,roptcassemblyinfo%indof
          jdofe = roptcassemblyinfo%Idofs(idofe,iel)
          roptcassemblyinfo%DpvelDofsAlt(idofe,iel,1) = p_DpvelocityX(jdofe)
          roptcassemblyinfo%DpvelDofsAlt(idofe,iel,2) = p_DpvelocityY(jdofe)
          roptcassemblyinfo%DdvelDofsAlt(idofe,iel,1) = p_DdvelocityX(jdofe)
          roptcassemblyinfo%DdvelDofsAlt(idofe,iel,2) = p_DdvelocityY(jdofe)
        end do
      end do
    end if
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine conv_strdiffOptC2ddoneasm (roptcassemblyinfo)
  
!<description>
  ! Clean up the roptcassemblyinfo structure.
!</description>

!<inputoutput>
  ! Structure collecting information about the assembly.
  type(t_optcassemblyinfo), intent(inout) :: roptcassemblyinfo
!</inputoutput>

!</subroutine>

    deallocate(roptcassemblyinfo%Dbas)
    deallocate(roptcassemblyinfo%Idofs)
    deallocate(roptcassemblyinfo%DlocalDeltaPrimal)
    deallocate(roptcassemblyinfo%DlocalDeltaDual)
    deallocate(roptcassemblyinfo%Kentry)
    deallocate(roptcassemblyinfo%Kentry12)
    
    ! Allocate memory for the velocites in the cubature points.
    deallocate(roptcassemblyinfo%Dpvel)
    deallocate(roptcassemblyinfo%Ddvel)
    deallocate(roptcassemblyinfo%DpvelXderiv)
    deallocate(roptcassemblyinfo%DpvelYderiv)
    deallocate(roptcassemblyinfo%DdvelXderiv)
    deallocate(roptcassemblyinfo%DdvelYderiv)
    
    if (associated(roptcassemblyinfo%DpvelDofsAlt,roptcassemblyinfo%DpvelDofs)) then
      deallocate(roptcassemblyinfo%DpvelDofsAlt)
      deallocate(roptcassemblyinfo%DdvelDofsAlt)
    else
      nullify(roptcassemblyinfo%DpvelDofsAlt)
      nullify(roptcassemblyinfo%DdvelDofsAlt)
    end if
    
    deallocate(roptcassemblyinfo%DpvelDofs)
    deallocate(roptcassemblyinfo%DdvelDofs)
    
    deallocate(roptcassemblyinfo%p_DcubPtsRef)

    ! Clean up of the element set.
    call elprep_releaseElementSet(roptcassemblyinfo%revalElementSet)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  subroutine conv_strdiffOptC2dgetMatrix (rmatrix,roptcoperator,dweight,&
      rvelocityVectorPrimal,rvelocityVectorDual,rcubatureInfo)
!<description>
  ! Calculate the matrix of the nonlinear operator:
  !   dweight*A(rvelocityVector)
!</description>

!<input>

  ! Structure defining the operator to set up.
  type(t_optcoperator), intent(in) :: roptcoperator
  
  ! Weight of the operator. Standard = 1.0
  real(dp), intent(in) :: dweight
      
  ! Velocity vector for the nonlinearity.
  ! The first blocks 1/2 in this vector define the evaluation
  ! point (primal velocity).
  type(t_vectorBlock), intent(in) :: rvelocityVectorPrimal

  ! Velocity vector for the nonlinearity.
  ! The first blocks 4/5 in this vector define the evaluation
  ! point (dual velocity).
  type(t_vectorBlock), intent(in) :: rvelocityVectorDual
  
  ! OPTIONAL: A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), target :: rcubatureInfo

!</input>

!<inputoutput>
  ! The system matrix. The submatrices for the velocity must be in block
  ! A11, A12, A21 and A22 and must be in matrix format 7 or 9.
  ! A11 and A22 must have the same structure. A12 and A21 must have
  ! the same structure.
  type(t_matrixBlock), intent(inout), target :: rmatrix
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: ielset,ielmax
  
  ! Matrix structure arrays
  real(DP), dimension(:), pointer :: p_Da11,p_Da22,p_Da44,p_Da55
  real(DP), dimension(:), pointer :: p_Da12,p_Da21,p_Da45,p_Da54
  real(DP), dimension(:), pointer :: p_Da41,p_Da51,p_Da42,p_Da52
  real(DP), dimension(:), pointer :: p_Da14,p_Da15,p_Da24,p_Da25
  
  ! Local matrices
  real(DP), dimension(:,:,:), allocatable :: DentryA11
  real(DP), dimension(:,:,:), allocatable :: DentryA12
  real(DP), dimension(:,:,:), allocatable :: DentryA21
  real(DP), dimension(:,:,:), allocatable :: DentryA22
  
  real(DP), dimension(:,:,:), allocatable :: DentryA44
  real(DP), dimension(:,:,:), allocatable :: DentryA45
  real(DP), dimension(:,:,:), allocatable :: DentryA54
  real(DP), dimension(:,:,:), allocatable :: DentryA55
  
  real(DP), dimension(:,:,:), allocatable :: DentryA41
  real(DP), dimension(:,:,:), allocatable :: DentryA52
  real(DP), dimension(:,:,:), allocatable :: DentryA42
  real(DP), dimension(:,:,:), allocatable :: DentryA51
  
  real(DP), dimension(:,:,:), allocatable :: DentryA14
  real(DP), dimension(:,:,:), allocatable :: DentryA25
  real(DP), dimension(:,:,:), allocatable :: DentryA24
  real(DP), dimension(:,:,:), allocatable :: DentryA15
  
  ! Temp vector
  type(t_vectorBlock) :: rvectorBlock
  
  ! Assembly status structure
  type(t_optcassemblyinfo) :: roptcassemblyinfo
  
  ! Cubature information structure
  type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

    if (roptcoperator%dnu .eq. 0.0_DP) then
      print *,'SD: NU=0 not allowed! Set dbeta=0 to prevent Stokes operator'// &
              ' from being build!'
      call sys_halt()
    end if
    
    ! Get pointers to the matrix content
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,1),p_Da11)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,2),p_Da22)

    nullify(p_Da12,p_Da21)
    if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,2),p_Da12)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,1),p_Da21)
    end if

    call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,4),p_Da44)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,5),p_Da55)

    nullify(p_Da45,p_Da54)
    if (lsysbl_isSubmatrixPresent(rmatrix,4,5)) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,5),p_Da45)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,4),p_Da54)
    end if

    call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,1),p_Da41)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,2),p_Da52)
    
    nullify(p_Da42,p_Da51)
    if (lsysbl_isSubmatrixPresent(rmatrix,4,2)) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,2),p_Da42)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,1),p_Da51)
    end if

    call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,4),p_Da14)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,5),p_Da25)

    nullify(p_Da24,p_Da15)
    if (lsysbl_isSubmatrixPresent(rmatrix,2,4)) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,4),p_Da24)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,5),p_Da15)
    end if
    
    ! Initialise the asembly of the local matrices.
    call conv_strdiffOptC2dinitasm (rvelocityVectorPrimal%p_rblockDiscr%RspatialDiscr(1),&
        1,roptcassemblyinfo)

    ! Calculate the maximum norm of the actual velocity field
    ! Round up the norm to 1D-8 if it's too small...
    call lsysbl_deriveSubvector(rvelocityVectorPrimal,rvectorBlock,1,2,.true.)
    call lsysbl_getVectorMagnitude (rvectorBlock,dumax=roptcassemblyinfo%dumax)
    call lsysbl_releaseVector (rvectorBlock)
  
    if (roptcassemblyinfo%dumax .lt. 1E-8_DP) roptcassemblyinfo%dumax=1E-8_DP

    ! Allocate memory for the local matrices
    allocate(DentryA11(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA12(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA21(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA22(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))

    allocate(DentryA44(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA45(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA54(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA55(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))

    allocate(DentryA41(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA52(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA42(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA51(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))

    allocate(DentryA14(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA25(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA24(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA15(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
  
!    ! Do we have an assembly structure?
!    ! If we do not have it, create a cubature info structure that
!    ! defines how to do the assembly.
!    if (.not. present(rcubatureInfo)) then
!      call spdiscr_createDefCubStructure(p_rdiscretisation,&
!          rtempCubatureInfo,CUB_GEN_DEPR_EVAL)
!      p_rcubatureInfo => rtempCubatureInfo
!    else
      p_rcubatureInfo => rcubatureInfo
!    end if

    ! Initialise the array with the local Delta values for the stabilisation
    call lalg_clearVectorDble (roptcassemblyinfo%DlocalDeltaPrimal)
    call lalg_clearVectorDble (roptcassemblyinfo%DlocalDeltaDual)

    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so BILF_NELEMSIM local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !%OMP do SCHEDULE(dynamic,1)
    do ielset = 1, size(roptcassemblyinfo%p_IelementList), BILF_NELEMSIM

      ! We always handle BILF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
      ! elements simultaneously.
      
      ielmax = min(size(roptcassemblyinfo%p_IelementList),ielset-1+BILF_NELEMSIM)
    
      ! Initialise the element set, compute the basis functions in the
      ! cubature points.
      if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
        call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
            rmatrix%RmatrixBlock(1,2),rvelocityVectorPrimal,rvelocityVectorDual,&
            roptcoperator,roptcassemblyinfo,ielset,ielmax)
      else
        call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
            rmatrix%RmatrixBlock(1,1),rvelocityVectorPrimal,rvelocityVectorDual,&
            roptcoperator,roptcassemblyinfo,ielset,ielmax)
      end if

      ! Clear the local matrices. If the Newton part is to be calculated,
      ! we must clear everything, otherwise only Dentry.
      DentryA11 = 0.0_DP
      DentryA12 = 0.0_DP
      DentryA21 = 0.0_DP
      DentryA22 = 0.0_DP
      
      DentryA44 = 0.0_DP
      DentryA45 = 0.0_DP
      DentryA54 = 0.0_DP
      DentryA55 = 0.0_DP
      
      DentryA41 = 0.0_DP
      DentryA52 = 0.0_DP
      DentryA42 = 0.0_DP
      DentryA51 = 0.0_DP
      
      DentryA14 = 0.0_DP
      DentryA25 = 0.0_DP
      DentryA24 = 0.0_DP
      DentryA15 = 0.0_DP
      
      ! Now calculate the local matrices on all the elements.
      call computeLocalOptCMatrices (roptcassemblyinfo%Dbas,&
          roptcassemblyinfo%Domega,roptcassemblyinfo%revalElementSet%p_Ddetj,&
          roptcassemblyinfo%indof,roptcassemblyinfo%ncubp,&
          roptcassemblyinfo%p_IelementList(ielset:ielmax),&
          roptcassemblyinfo%DpvelDofs,roptcassemblyinfo%DdvelDofs,&
          roptcassemblyinfo%Dpvel,roptcassemblyinfo%DpvelXderiv,&
          roptcassemblyinfo%DpvelYderiv,roptcassemblyinfo%Ddvel,&
          roptcassemblyinfo%DdvelXderiv,roptcassemblyinfo%DdvelYderiv,&
          DentryA11,DentryA22,DentryA44,DentryA55,&
          DentryA12,DentryA21,DentryA45,DentryA54,&
          DentryA41,DentryA52,DentryA42,DentryA51,&
          DentryA14,DentryA25,DentryA24,DentryA15,&
          roptcassemblyinfo%DlocalDeltaPrimal,roptcassemblyinfo%DlocalDeltaDual,&
          roptcoperator,roptcassemblyinfo)
      
      ! Incorporate the computed local matrices into the global matrix.
      !%OMP CRITICAL
      if (associated(p_Da11)) &
        call incorporateLocalMatrix (DentryA11,p_Da11,roptcassemblyinfo%Kentry,dweight)

      if (associated(p_Da12)) &
        call incorporateLocalMatrix (DentryA12,p_Da12,roptcassemblyinfo%Kentry12,dweight)

      if (associated(p_Da21)) &
        call incorporateLocalMatrix (DentryA21,p_Da21,roptcassemblyinfo%Kentry12,dweight)

      if (associated(p_Da22) .and. .not. associated (p_Da11,p_Da22)) &
        call incorporateLocalMatrix (DentryA22,p_Da22,roptcassemblyinfo%Kentry,dweight)


      if (associated(p_Da44)) &
        call incorporateLocalMatrix (DentryA44,p_Da44,roptcassemblyinfo%Kentry,dweight)

      if (associated(p_Da45)) &
        call incorporateLocalMatrix (DentryA45,p_Da45,roptcassemblyinfo%Kentry12,dweight)

      if (associated(p_Da54)) &
        call incorporateLocalMatrix (DentryA54,p_Da54,roptcassemblyinfo%Kentry12,dweight)

      if (associated(p_Da55)) &
        call incorporateLocalMatrix (DentryA55,p_Da55,roptcassemblyinfo%Kentry,dweight)
      
      
      if (associated(p_Da41)) &
        call incorporateLocalMatrix (DentryA41,p_Da41,roptcassemblyinfo%Kentry,dweight)

      if (associated(p_Da42)) &
        call incorporateLocalMatrix (DentryA42,p_Da42,roptcassemblyinfo%Kentry12,dweight)

      if (associated(p_Da51)) &
        call incorporateLocalMatrix (DentryA51,p_Da51,roptcassemblyinfo%Kentry12,dweight)

      if (associated(p_Da52)) &
        call incorporateLocalMatrix (DentryA52,p_Da52,roptcassemblyinfo%Kentry,dweight)
      
      
      if (associated(p_Da14)) &
        call incorporateLocalMatrix (DentryA14,p_Da14,roptcassemblyinfo%Kentry,dweight)
                                                    
      if (associated(p_Da24)) &
        call incorporateLocalMatrix (DentryA24,p_Da24,roptcassemblyinfo%Kentry12,dweight)
                                                    
      if (associated(p_Da15)) &
        call incorporateLocalMatrix (DentryA15,p_Da15,roptcassemblyinfo%Kentry12,dweight)
                                                    
      if (associated(p_Da25)) &
        call incorporateLocalMatrix (DentryA25,p_Da25,roptcassemblyinfo%Kentry,dweight)
      

    end do ! ielset
    !%OMP end do
    
!    ! Release the assembly structure if necessary.
!    if (.not. present(rcubatureInfo)) then
!      call spdiscr_releaseCubStructure(rtempCubatureInfo)
!    end if
    
    ! Release memory
    deallocate(DentryA11)
    deallocate(DentryA22)
    if (allocated(DentryA12)) deallocate(DentryA12)
    if (allocated(DentryA21)) deallocate(DentryA21)
               
    deallocate(DentryA44)
    deallocate(DentryA55)
    if (allocated(DentryA45)) deallocate(DentryA45)
    if (allocated(DentryA54)) deallocate(DentryA54)
               
    deallocate(DentryA41)
    deallocate(DentryA52)
    if (allocated(DentryA42)) deallocate(DentryA42)
    if (allocated(DentryA51)) deallocate(DentryA51)
               
    deallocate(DentryA14)
    deallocate(DentryA25)
    if (allocated(DentryA24)) deallocate(DentryA24)
    if (allocated(DentryA15)) deallocate(DentryA15)
               
  contains
    
    subroutine incorporateLocalMatrix (DaLocal,Da,Kentry,dweight)
    
    ! Incorporate a local matrix into a global one.
    
    ! The local matrix to incorporate.
    real(dp), dimension(:,:,:), intent(in) :: DaLocal
    
    ! The global matrix
    real(dp), dimension(:), intent(inout) :: Da
    
    ! Positions of the local matrix in the global one.
    integer, dimension(:,:,:), intent(in) :: Kentry
    
    ! Weight of the local matrix
    real(dp) :: dweight

      ! local variables
      integer :: iel,idofe,jdofe
    
      do iel=1,ubound(DaLocal,3)
        do idofe=1,ubound(DaLocal,2)
          do jdofe=1,ubound(DaLocal,1)
            Da(Kentry(jdofe,idofe,iel)) = Da(Kentry(jdofe,idofe,iel)) + &
                dweight * DaLocal(jdofe,idofe,iel)
          end do
        end do
      end do
    
    end subroutine
    
  end subroutine
     
  ! ***************************************************************************

!<subroutine>
  subroutine conv_strdiffOptC2dgetDerMatrix (rmatrix,roptcoperator,dweight,&
      rvelocityVectorPrimal,rvelocityVectorDual,dh,&
      rvelocityVectorPrimalAlt,rvelocityVectorDualAlt,rcubatureInfo)
!<description>
  ! Calculate the derivative matrix of the nonlinear operator:
  !   dweight*A(rvelocityVector)
!</description>

!<remarks>
  ! rvelocityVectorPrimalAlt / rvelocityVectorDualAlt allows to specify an
  ! alternative velocity in the following sense:
  ! If NOT specified, the routine assumes that an operator
  !   R(x) = A(x)x
  ! is given and calculates the derivative with respect to x. On the other hand,
  ! if y:= rvelocityVectorPrimalAlt/rvelocityVectorDualAlt is specified,
  ! the routine assumes that the operator is given in the form
  !   R(x) = R(x,y) = A(x)y
  ! and calculates the derivative only with respect to x.
  ! This can be used to calculate the derivative matrix of subblocks in
  ! a block matrix, where the block A(x) changes while the corresponding
  ! evaluation point y stays fixed.
!</remarks>

!<input>

  ! Structure defining the operator to set up.
  type(t_optcoperator), intent(in) :: roptcoperator
  
  ! Weight of the operator. Standard = 1.0
  real(dp), intent(in) :: dweight
      
  ! Velocity vector for the nonlinearity.
  ! The first blocks 1/2 in this vector define the evaluation
  ! point (primal velocity).
  type(t_vectorBlock), intent(in) :: rvelocityVectorPrimal

  ! Velocity vector for the nonlinearity.
  ! The first blocks 4/5 in this vector define the evaluation
  ! point (dual velocity).
  type(t_vectorBlock), intent(in) :: rvelocityVectorDual
  
  ! Step length for the discrete derivative.
  real(dp), intent(in) :: dh

  ! OPTIONAL: A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), target :: rcubatureInfo

  ! Alternative Velocity vector for the nonlinearity to be multiplied to the
  ! matrix. The first blocks 1/2 in this vector define the evaluation
  ! point (primal velocity).
  ! If not present, this coincides with rvelocityVectorPrimal
  type(t_vectorBlock), intent(in), optional :: rvelocityVectorPrimalAlt

  ! Alternative Velocity vector for the nonlinearity to be multiplied to the
  ! matrix. The first blocks 4/5 in this vector define the evaluation
  ! point (dual velocity).
  ! If not present, this coincides with rvelocityVectorDual.
  type(t_vectorBlock), intent(in), optional :: rvelocityVectorDualAlt
  
!</input>

!<inputoutput>
  ! The system matrix. The submatrices for the velocity must be in block
  ! A11, A12, A21 and A22 and must be in matrix format 7 or 9.
  ! A11 and A22 must have the same structure. A12 and A21 must have
  ! the same structure.
  type(t_matrixBlock), intent(inout), target :: rmatrix
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: ielset,ielmax,idof,iel
  real(dp) :: dweight2
  
  ! Matrix structure arrays
  real(DP), dimension(:), pointer :: p_Da11,p_Da22,p_Da44,p_Da55
  real(DP), dimension(:), pointer :: p_Da12,p_Da21,p_Da45,p_Da54
  real(DP), dimension(:), pointer :: p_Da41,p_Da51,p_Da42,p_Da52
  real(DP), dimension(:), pointer :: p_Da14,p_Da15,p_Da24,p_Da25
  
  ! Local matrices
  real(DP), dimension(:,:,:), allocatable :: DentryA11
  real(DP), dimension(:,:,:), allocatable :: DentryA12
  real(DP), dimension(:,:,:), allocatable :: DentryA21
  real(DP), dimension(:,:,:), allocatable :: DentryA22
                                                      
  real(DP), dimension(:,:,:), allocatable :: DentryA44
  real(DP), dimension(:,:,:), allocatable :: DentryA45
  real(DP), dimension(:,:,:), allocatable :: DentryA54
  real(DP), dimension(:,:,:), allocatable :: DentryA55
                                                      
  real(DP), dimension(:,:,:), allocatable :: DentryA41
  real(DP), dimension(:,:,:), allocatable :: DentryA52
  real(DP), dimension(:,:,:), allocatable :: DentryA42
  real(DP), dimension(:,:,:), allocatable :: DentryA51
                                                      
  real(DP), dimension(:,:,:), allocatable :: DentryA14
  real(DP), dimension(:,:,:), allocatable :: DentryA25
  real(DP), dimension(:,:,:), allocatable :: DentryA24
  real(DP), dimension(:,:,:), allocatable :: DentryA15
  real(DP), dimension(:), allocatable :: Dtemp
  
  ! Temp vector
  type(t_vectorBlock) :: rvectorBlock
  
  ! Assembly status structure
  type(t_optcassemblyinfo) :: roptcassemblyinfo

  ! Cubature information structure
  type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

    if (roptcoperator%dnu .eq. 0.0_DP) then
      print *,'SD: NU=0 not allowed! Set dbeta=0 to prevent Stokes operator'// &
              ' from being build!'
      call sys_halt()
    end if
    
    ! Get pointers to the matrix content
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,1),p_Da11)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,2),p_Da22)

    nullify(p_Da12,p_Da21)
    if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,2),p_Da12)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,1),p_Da21)
    end if

    call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,4),p_Da44)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,5),p_Da55)

    nullify(p_Da45,p_Da54)
    if (lsysbl_isSubmatrixPresent(rmatrix,4,5)) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,5),p_Da45)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,4),p_Da54)
    end if

    call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,1),p_Da41)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,2),p_Da52)
    
    nullify(p_Da42,p_Da51)
    if (lsysbl_isSubmatrixPresent(rmatrix,4,2)) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,2),p_Da42)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,1),p_Da51)
    end if

    call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,4),p_Da14)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,5),p_Da25)

    nullify(p_Da24,p_Da15)
    if (lsysbl_isSubmatrixPresent(rmatrix,2,4)) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,4),p_Da24)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,5),p_Da15)
    end if
    
    ! Initialise the asembly of the local matrices.
    call conv_strdiffOptC2dinitasm (rvelocityVectorPrimal%p_rblockDiscr%RspatialDiscr(1),&
        1,roptcassemblyinfo,present(rvelocityVectorPrimalAlt))

    ! Calculate the maximum norm of the actual velocity field
    ! Round up the norm to 1D-8 if it's too small...
    call lsysbl_deriveSubvector(rvelocityVectorPrimal,rvectorBlock,1,2,.true.)
    call lsysbl_getVectorMagnitude (rvectorBlock,dumax=roptcassemblyinfo%dumax)
    call lsysbl_releaseVector (rvectorBlock)
  
    if (roptcassemblyinfo%dumax .lt. 1E-8_DP) roptcassemblyinfo%dumax=1E-8_DP

    ! Allocate memory for the local matrices
    allocate(DentryA11(roptcassemblyinfo%indof,roptcassemblyinfo%indof,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA12(roptcassemblyinfo%indof,roptcassemblyinfo%indof,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA21(roptcassemblyinfo%indof,roptcassemblyinfo%indof,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA22(roptcassemblyinfo%indof,roptcassemblyinfo%indof,&
        roptcassemblyinfo%nelementsPerBlock))

    allocate(DentryA44(roptcassemblyinfo%indof,roptcassemblyinfo%indof,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA45(roptcassemblyinfo%indof,roptcassemblyinfo%indof,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA54(roptcassemblyinfo%indof,roptcassemblyinfo%indof,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA55(roptcassemblyinfo%indof,roptcassemblyinfo%indof,&
        roptcassemblyinfo%nelementsPerBlock))

    allocate(DentryA41(roptcassemblyinfo%indof,roptcassemblyinfo%indof,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA52(roptcassemblyinfo%indof,roptcassemblyinfo%indof,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA42(roptcassemblyinfo%indof,roptcassemblyinfo%indof,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA51(roptcassemblyinfo%indof,roptcassemblyinfo%indof,&
        roptcassemblyinfo%nelementsPerBlock))

    allocate(DentryA14(roptcassemblyinfo%indof,roptcassemblyinfo%indof,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA25(roptcassemblyinfo%indof,roptcassemblyinfo%indof,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA24(roptcassemblyinfo%indof,roptcassemblyinfo%indof,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA15(roptcassemblyinfo%indof,roptcassemblyinfo%indof,&
        roptcassemblyinfo%nelementsPerBlock))
        
    ! Temp array
    allocate(Dtemp(roptcassemblyinfo%indof))

!    ! Do we have an assembly structure?
!    ! If we do not have it, create a cubature info structure that
!    ! defines how to do the assembly.
!    if (.not. present(rcubatureInfo)) then
!      call spdiscr_createDefCubStructure(p_rdiscretisation,&
!          rtempCubatureInfo,CUB_GEN_DEPR_EVAL)
!      p_rcubatureInfo => rtempCubatureInfo
!    else
      p_rcubatureInfo => rcubatureInfo
!    end if

    ! Initialise the array with the local Delta values for the stabilisation
    call lalg_clearVectorDble (roptcassemblyinfo%DlocalDeltaPrimal)
    call lalg_clearVectorDble (roptcassemblyinfo%DlocalDeltaDual)

    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so BILF_NELEMSIM local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !%OMP do SCHEDULE(dynamic,1)
    do ielset = 1, size(roptcassemblyinfo%p_IelementList), BILF_NELEMSIM

      ! We always handle BILF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
      ! elements simultaneously.
      
      ielmax = min(size(roptcassemblyinfo%p_IelementList),ielset-1+BILF_NELEMSIM)
    
      ! To assemble the discrete Newton matrix, we have to
      ! * Assemble the matrix when modifying each DOF on an element
      ! * Multiplication with the corresponding solution vector gives one
      !   column of the Newton matrix
      !
      ! The global matrix structure of the main matrix looks like:
      !
      !   (A11     A41    )
      !   (    A22     A52)
      !   (A41     A44 A45)
      !   (    A52 A54 A55)
      !
      ! The corresponding Newton matrix has the shape
      !
      !   (A11 A12 A41    )
      !   (A21 A22     A52)
      !   (A41 A42 A44 A45)
      !   (A51 A52 A54 A55)
      !
      ! as A41/A52 are linear.
      !
      ! We calculate the matrix column-wise and element based:
      ! Let Bij := 1/2h * (Aij(u+h*e_k)) and
      !     Cij := 1/2h * (Aij(u-h*e_k)), then we have with
      ! (pdofp1,pdofp2,ddofp1,ddofp2) = u+h*e_k
      ! (pdofn1,pdofn2,ddofn1,ddofn2) = u-h*e_k
      !
      !   (B11     B41    ) (pdofp1) - (C11     C41    ) (pdofn1)  ->  (    X          )
      !   (    B22     B52) (pdofp2)   (    C22     C52) (pdofn2)      (    X          )
      !   (B41     B44 B45) (ddofp1)   (C41     C44 C45) (ddofn1)      (    X          )
      !   (    B52 B54 B55) (ddofp2)   (    C52 C54 C55) (ddofn2)      (    X          )
      !                                                                     ^k
    
      ! Loop through the DOF's. All DOF's must be once increased and once decreased
      ! by h.
      dweight2 = dweight*0.5_DP/dh
      do idof = 1,roptcassemblyinfo%indof
    
        ! Initialise the element set, compute the basis functions in the
        ! cubature points.
        if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,2),rvelocityVectorPrimal,rvelocityVectorDual,&
              roptcoperator,roptcassemblyinfo,ielset,ielmax,&
            rvelocityVectorPrimalAlt,rvelocityVectorDualAlt)
        else
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,1),rvelocityVectorPrimal,rvelocityVectorDual,&
              roptcoperator,roptcassemblyinfo,ielset,ielmax,&
            rvelocityVectorPrimalAlt,rvelocityVectorDualAlt)
        end if
        
        ! Increase the idof'th entry in the velocity evaluation point by h.
        do iel = 1,ielmax
          roptcassemblyinfo%DpvelDofs(idof,iel,1) = roptcassemblyinfo%DpvelDofs(idof,iel,1) + dh
        end do

        ! Clear the local matrices. If the Newton part is to be calculated,
        ! we must clear everything, otherwise only Dentry.
        DentryA11 = 0.0_DP
        DentryA12 = 0.0_DP
        DentryA21 = 0.0_DP
        DentryA22 = 0.0_DP
        
        DentryA44 = 0.0_DP
        DentryA45 = 0.0_DP
        DentryA54 = 0.0_DP
        DentryA55 = 0.0_DP
        
        DentryA41 = 0.0_DP
        DentryA52 = 0.0_DP
        DentryA42 = 0.0_DP
        DentryA51 = 0.0_DP
        
        DentryA14 = 0.0_DP
        DentryA25 = 0.0_DP
        DentryA24 = 0.0_DP
        DentryA15 = 0.0_DP
        
        ! Now calculate the local matrices on all the elements.
        call computeLocalOptCMatrices (roptcassemblyinfo%Dbas,&
            roptcassemblyinfo%Domega,roptcassemblyinfo%revalElementSet%p_Ddetj,&
            roptcassemblyinfo%indof,roptcassemblyinfo%ncubp,&
            roptcassemblyinfo%p_IelementList(ielset:ielmax),&
            roptcassemblyinfo%DpvelDofs,roptcassemblyinfo%DdvelDofs,&
            roptcassemblyinfo%Dpvel,roptcassemblyinfo%DpvelXderiv,&
            roptcassemblyinfo%DpvelYderiv,roptcassemblyinfo%Ddvel,&
            roptcassemblyinfo%DdvelXderiv,roptcassemblyinfo%DdvelYderiv,&
            DentryA11,DentryA22,DentryA44,DentryA55,&
            DentryA12,DentryA21,DentryA45,DentryA54,&
            DentryA41,DentryA52,DentryA42,DentryA51,&
            DentryA14,DentryA25,DentryA24,DentryA15,&
            roptcassemblyinfo%DlocalDeltaPrimal,roptcassemblyinfo%DlocalDeltaDual,&
            roptcoperator,roptcassemblyinfo)
        
        ! Incorporate the computed local matrices into the global matrix.
        if (associated(p_Da11)) then
          call incorporateLocalMatVecCol (DentryA11,p_Da11,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,1,idof,dweight2,Dtemp)

          if (roptcoperator%ccontrolProjection .eq. 0) then
            call incorporateLocalMatVecCol (DentryA14,p_Da11,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)
          else
            select case (roptcoperator%cconstraintsType)
            case (0)
              call incorporateLocalMatVecCol (DentryA14,p_Da11,roptcassemblyinfo%Kentry,&
                                              roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp,&
                                              roptcoperator%dmax1/roptcoperator%dcontrolMultiplier,&
                                              roptcoperator%dmin1/roptcoperator%dcontrolMultiplier)
            case default
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
          end if
        end if

        if (associated(p_Da21)) then
          call incorporateLocalMatVecCol (DentryA22,p_Da21,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,2,idof,dweight2,Dtemp)

          if (roptcoperator%ccontrolProjection .eq. 0) then
            call incorporateLocalMatVecCol (DentryA25,p_Da21,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
          else
            select case (roptcoperator%cconstraintsType)
            case (0)
              call incorporateLocalMatVecCol (DentryA25,p_Da21,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp,&
                                            roptcoperator%dmax2/roptcoperator%dcontrolMultiplier,&
                                            roptcoperator%dmin2/roptcoperator%dcontrolMultiplier)
            case default
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
          end if
        end if


        if (associated(p_Da41)) then
          call incorporateLocalMatVecCol (DentryA41,p_Da41,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA44,p_Da41,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA45,p_Da41,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
        end if


        if (associated(p_Da51)) then
          call incorporateLocalMatVecCol (DentryA52,p_Da51,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,2,idof,dweight2,Dtemp)
                                                      
          call incorporateLocalMatVecCol (DentryA54,p_Da51,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA55,p_Da51,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
        end if
            
      end do
      
      dweight2 = -dweight*0.5_DP/dh
      do idof = 1,roptcassemblyinfo%indof
    
        ! Initialise the element set, compute the basis functions in the
        ! cubature points.
        if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,2),rvelocityVectorPrimal,rvelocityVectorDual,&
              roptcoperator,roptcassemblyinfo,ielset,ielmax,&
            rvelocityVectorPrimalAlt,rvelocityVectorDualAlt)
        else
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,1),rvelocityVectorPrimal,rvelocityVectorDual,&
              roptcoperator,roptcassemblyinfo,ielset,ielmax,&
            rvelocityVectorPrimalAlt,rvelocityVectorDualAlt)
        end if
        
        ! Increase the idof'th entry in the velocity evaluation point by h.
        do iel = 1,ielmax
          roptcassemblyinfo%DpvelDofs(idof,iel,1) = roptcassemblyinfo%DpvelDofs(idof,iel,1) - dh
        end do

        ! Clear the local matrices. If the Newton part is to be calculated,
        ! we must clear everything, otherwise only Dentry.
        DentryA11 = 0.0_DP
        DentryA12 = 0.0_DP
        DentryA21 = 0.0_DP
        DentryA22 = 0.0_DP
        
        DentryA44 = 0.0_DP
        DentryA45 = 0.0_DP
        DentryA54 = 0.0_DP
        DentryA55 = 0.0_DP
        
        DentryA41 = 0.0_DP
        DentryA52 = 0.0_DP
        DentryA42 = 0.0_DP
        DentryA51 = 0.0_DP
        
        DentryA14 = 0.0_DP
        DentryA25 = 0.0_DP
        DentryA24 = 0.0_DP
        DentryA15 = 0.0_DP
        
        ! Now calculate the local matrices on all the elements.
        call computeLocalOptCMatrices (roptcassemblyinfo%Dbas,&
            roptcassemblyinfo%Domega,roptcassemblyinfo%revalElementSet%p_Ddetj,&
            roptcassemblyinfo%indof,roptcassemblyinfo%ncubp,&
            roptcassemblyinfo%p_IelementList(ielset:ielmax),&
            roptcassemblyinfo%DpvelDofs,roptcassemblyinfo%DdvelDofs,&
            roptcassemblyinfo%Dpvel,roptcassemblyinfo%DpvelXderiv,&
            roptcassemblyinfo%DpvelYderiv,roptcassemblyinfo%Ddvel,&
            roptcassemblyinfo%DdvelXderiv,roptcassemblyinfo%DdvelYderiv,&
            DentryA11,DentryA22,DentryA44,DentryA55,&
            DentryA12,DentryA21,DentryA45,DentryA54,&
            DentryA41,DentryA52,DentryA42,DentryA51,&
            DentryA14,DentryA25,DentryA24,DentryA15,&
            roptcassemblyinfo%DlocalDeltaPrimal,roptcassemblyinfo%DlocalDeltaDual,&
            roptcoperator,roptcassemblyinfo)
        
        ! Incorporate the computed local matrices into the global matrix.
        if (associated(p_Da11)) then
          call incorporateLocalMatVecCol (DentryA11,p_Da11,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,1,idof,dweight2,Dtemp)

          if (roptcoperator%ccontrolProjection .eq. 0) then
            call incorporateLocalMatVecCol (DentryA14,p_Da11,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)
          else
            select case (roptcoperator%cconstraintsType)
            case (0)
              call incorporateLocalMatVecCol (DentryA14,p_Da11,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp,&
                                            roptcoperator%dmax1/roptcoperator%dcontrolMultiplier,&
                                            roptcoperator%dmin1/roptcoperator%dcontrolMultiplier)
            case default
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
          end if
        end if

        if (associated(p_Da21)) then
          call incorporateLocalMatVecCol (DentryA22,p_Da21,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,2,idof,dweight2,Dtemp)

          if (roptcoperator%ccontrolProjection .eq. 0) then
            call incorporateLocalMatVecCol (DentryA25,p_Da21,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
          else
            select case (roptcoperator%cconstraintsType)
            case (0)
              call incorporateLocalMatVecCol (DentryA25,p_Da21,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp,&
                                            roptcoperator%dmax2/roptcoperator%dcontrolMultiplier,&
                                            roptcoperator%dmin2/roptcoperator%dcontrolMultiplier)
            case default
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
          end if
        end if


        if (associated(p_Da41)) then
          call incorporateLocalMatVecCol (DentryA41,p_Da41,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA44,p_Da41,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA45,p_Da41,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
        end if


        if (associated(p_Da51)) then
          call incorporateLocalMatVecCol (DentryA52,p_Da51,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,2,idof,dweight2,Dtemp)
                                                      
          call incorporateLocalMatVecCol (DentryA54,p_Da51,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA55,p_Da51,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
        end if

      end do
      
      ! Loop through the DOF's. All DOF's must be once increased and once decreased
      ! by h.
      dweight2 = dweight*0.5_DP/dh
      do idof = 1,roptcassemblyinfo%indof
    
        ! Initialise the element set, compute the basis functions in the
        ! cubature points.
        if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,2),rvelocityVectorPrimal,rvelocityVectorDual,&
              roptcoperator,roptcassemblyinfo,ielset,ielmax,&
            rvelocityVectorPrimalAlt,rvelocityVectorDualAlt)
        else
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,1),rvelocityVectorPrimal,rvelocityVectorDual,&
              roptcoperator,roptcassemblyinfo,ielset,ielmax,&
            rvelocityVectorPrimalAlt,rvelocityVectorDualAlt)
        end if
        
        ! Increase the idof'th entry in the velocity evaluation point by h.
        do iel = 1,ielmax
          roptcassemblyinfo%DpvelDofs(idof,iel,2) = roptcassemblyinfo%DpvelDofs(idof,iel,2) + dh
        end do

        ! Clear the local matrices. If the Newton part is to be calculated,
        ! we must clear everything, otherwise only Dentry.
        DentryA11 = 0.0_DP
        DentryA12 = 0.0_DP
        DentryA21 = 0.0_DP
        DentryA22 = 0.0_DP
        
        DentryA44 = 0.0_DP
        DentryA45 = 0.0_DP
        DentryA54 = 0.0_DP
        DentryA55 = 0.0_DP
        
        DentryA41 = 0.0_DP
        DentryA52 = 0.0_DP
        DentryA42 = 0.0_DP
        DentryA51 = 0.0_DP
        
        DentryA14 = 0.0_DP
        DentryA25 = 0.0_DP
        DentryA24 = 0.0_DP
        DentryA15 = 0.0_DP
        
        ! Now calculate the local matrices on all the elements.
        call computeLocalOptCMatrices (roptcassemblyinfo%Dbas,&
            roptcassemblyinfo%Domega,roptcassemblyinfo%revalElementSet%p_Ddetj,&
            roptcassemblyinfo%indof,roptcassemblyinfo%ncubp,&
            roptcassemblyinfo%p_IelementList(ielset:ielmax),&
            roptcassemblyinfo%DpvelDofs,roptcassemblyinfo%DdvelDofs,&
            roptcassemblyinfo%Dpvel,roptcassemblyinfo%DpvelXderiv,&
            roptcassemblyinfo%DpvelYderiv,roptcassemblyinfo%Ddvel,&
            roptcassemblyinfo%DdvelXderiv,roptcassemblyinfo%DdvelYderiv,&
            DentryA11,DentryA22,DentryA44,DentryA55,&
            DentryA12,DentryA21,DentryA45,DentryA54,&
            DentryA41,DentryA52,DentryA42,DentryA51,&
            DentryA14,DentryA25,DentryA24,DentryA15,&
            roptcassemblyinfo%DlocalDeltaPrimal,roptcassemblyinfo%DlocalDeltaDual,&
            roptcoperator,roptcassemblyinfo)
        
        ! Incorporate the computed local matrices into the global matrix.
        if (associated(p_Da12)) then
          call incorporateLocalMatVecCol (DentryA11,p_Da12,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,1,idof,dweight2,Dtemp)

          if (roptcoperator%ccontrolProjection .eq. 0) then
            call incorporateLocalMatVecCol (DentryA14,p_Da12,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)
          else
            select case (roptcoperator%cconstraintsType)
            case (0)
              call incorporateLocalMatVecCol (DentryA14,p_Da12,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp,&
                                            roptcoperator%dmax1/roptcoperator%dcontrolMultiplier,&
                                            roptcoperator%dmin1/roptcoperator%dcontrolMultiplier)
            case default
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
          end if
        end if

        if (associated(p_Da22)) then
          call incorporateLocalMatVecCol (DentryA22,p_Da22,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,2,idof,dweight2,Dtemp)

          if (roptcoperator%ccontrolProjection .eq. 0) then
            call incorporateLocalMatVecCol (DentryA25,p_Da22,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
          else
            select case (roptcoperator%cconstraintsType)
            case (0)
              call incorporateLocalMatVecCol (DentryA25,p_Da22,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp,&
                                            roptcoperator%dmax2/roptcoperator%dcontrolMultiplier,&
                                            roptcoperator%dmin2/roptcoperator%dcontrolMultiplier)
            case default
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
          end if
        end if


        if (associated(p_Da42)) then
          call incorporateLocalMatVecCol (DentryA41,p_Da42,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA44,p_Da42,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA45,p_Da42,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
        end if


        if (associated(p_Da52)) then
          call incorporateLocalMatVecCol (DentryA52,p_Da52,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,2,idof,dweight2,Dtemp)
                                                      
          call incorporateLocalMatVecCol (DentryA54,p_Da52,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA55,p_Da52,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
        end if
            
      end do

      dweight2 = -dweight*0.5_DP/dh
      do idof = 1,roptcassemblyinfo%indof
    
        ! Initialise the element set, compute the basis functions in the
        ! cubature points.
        if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,2),rvelocityVectorPrimal,rvelocityVectorDual,&
              roptcoperator,roptcassemblyinfo,ielset,ielmax,&
            rvelocityVectorPrimalAlt,rvelocityVectorDualAlt)
        else
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,1),rvelocityVectorPrimal,rvelocityVectorDual,&
              roptcoperator,roptcassemblyinfo,ielset,ielmax,&
            rvelocityVectorPrimalAlt,rvelocityVectorDualAlt)
        end if
        
        ! Increase the idof'th entry in the velocity evaluation point by h.
        do iel = 1,ielmax
          roptcassemblyinfo%DpvelDofs(idof,iel,2) = roptcassemblyinfo%DpvelDofs(idof,iel,2) - dh
        end do

        ! Clear the local matrices. If the Newton part is to be calculated,
        ! we must clear everything, otherwise only Dentry.
        DentryA11 = 0.0_DP
        DentryA12 = 0.0_DP
        DentryA21 = 0.0_DP
        DentryA22 = 0.0_DP
        
        DentryA44 = 0.0_DP
        DentryA45 = 0.0_DP
        DentryA54 = 0.0_DP
        DentryA55 = 0.0_DP
        
        DentryA41 = 0.0_DP
        DentryA52 = 0.0_DP
        DentryA42 = 0.0_DP
        DentryA51 = 0.0_DP
        
        DentryA14 = 0.0_DP
        DentryA25 = 0.0_DP
        DentryA24 = 0.0_DP
        DentryA15 = 0.0_DP
        
        ! Now calculate the local matrices on all the elements.
        call computeLocalOptCMatrices (roptcassemblyinfo%Dbas,&
            roptcassemblyinfo%Domega,roptcassemblyinfo%revalElementSet%p_Ddetj,&
            roptcassemblyinfo%indof,roptcassemblyinfo%ncubp,&
            roptcassemblyinfo%p_IelementList(ielset:ielmax),&
            roptcassemblyinfo%DpvelDofs,roptcassemblyinfo%DdvelDofs,&
            roptcassemblyinfo%Dpvel,roptcassemblyinfo%DpvelXderiv,&
            roptcassemblyinfo%DpvelYderiv,roptcassemblyinfo%Ddvel,&
            roptcassemblyinfo%DdvelXderiv,roptcassemblyinfo%DdvelYderiv,&
            DentryA11,DentryA22,DentryA44,DentryA55,&
            DentryA12,DentryA21,DentryA45,DentryA54,&
            DentryA41,DentryA52,DentryA42,DentryA51,&
            DentryA14,DentryA25,DentryA24,DentryA15,&
            roptcassemblyinfo%DlocalDeltaPrimal,roptcassemblyinfo%DlocalDeltaDual,&
            roptcoperator,roptcassemblyinfo)
        
        ! Incorporate the computed local matrices into the global matrix.
        if (associated(p_Da12)) then
          call incorporateLocalMatVecCol (DentryA11,p_Da12,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,1,idof,dweight2,Dtemp)

          if (roptcoperator%ccontrolProjection .eq. 0) then
            call incorporateLocalMatVecCol (DentryA14,p_Da12,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)
          else
            select case (roptcoperator%cconstraintsType)
            case (0)
              call incorporateLocalMatVecCol (DentryA14,p_Da12,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp,&
                                            roptcoperator%dmax1/roptcoperator%dcontrolMultiplier,&
                                            roptcoperator%dmin1/roptcoperator%dcontrolMultiplier)
            case default
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
          end if
        end if

        if (associated(p_Da22)) then
          call incorporateLocalMatVecCol (DentryA22,p_Da22,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,2,idof,dweight2,Dtemp)

          if (roptcoperator%ccontrolProjection .eq. 0) then
            call incorporateLocalMatVecCol (DentryA25,p_Da22,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
          else
            select case (roptcoperator%cconstraintsType)
            case (0)
              call incorporateLocalMatVecCol (DentryA25,p_Da22,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp,&
                                            roptcoperator%dmax2/roptcoperator%dcontrolMultiplier,&
                                            roptcoperator%dmin2/roptcoperator%dcontrolMultiplier)
            case default
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
          end if
        end if


        if (associated(p_Da42)) then
          call incorporateLocalMatVecCol (DentryA41,p_Da42,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA44,p_Da42,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA45,p_Da42,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
        end if


        if (associated(p_Da52)) then
          call incorporateLocalMatVecCol (DentryA52,p_Da52,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,2,idof,dweight2,Dtemp)
                                                      
          call incorporateLocalMatVecCol (DentryA54,p_Da52,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA55,p_Da52,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
        end if
            
      end do
     
      ! Loop through the DOF's. All DOF's must be once increased and once decreased
      ! by h.
      dweight2 = dweight*0.5_DP/dh
      do idof = 1,roptcassemblyinfo%indof
    
        ! Initialise the element set, compute the basis functions in the
        ! cubature points.
        if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,2),rvelocityVectorPrimal,rvelocityVectorDual,&
              roptcoperator,roptcassemblyinfo,ielset,ielmax,&
            rvelocityVectorPrimalAlt,rvelocityVectorDualAlt)
        else
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,1),rvelocityVectorPrimal,rvelocityVectorDual,&
              roptcoperator,roptcassemblyinfo,ielset,ielmax,&
            rvelocityVectorPrimalAlt,rvelocityVectorDualAlt)
        end if
        
        ! Increase the idof'th entry in the velocity evaluation point by h.
        do iel = 1,ielmax
          roptcassemblyinfo%DdvelDofs(idof,iel,1) = roptcassemblyinfo%DdvelDofs(idof,iel,1) + dh
        end do

        ! Clear the local matrices. If the Newton part is to be calculated,
        ! we must clear everything, otherwise only Dentry.
        DentryA11 = 0.0_DP
        DentryA12 = 0.0_DP
        DentryA21 = 0.0_DP
        DentryA22 = 0.0_DP
        
        DentryA44 = 0.0_DP
        DentryA45 = 0.0_DP
        DentryA54 = 0.0_DP
        DentryA55 = 0.0_DP
        
        DentryA41 = 0.0_DP
        DentryA52 = 0.0_DP
        DentryA42 = 0.0_DP
        DentryA51 = 0.0_DP
        
        DentryA14 = 0.0_DP
        DentryA25 = 0.0_DP
        DentryA24 = 0.0_DP
        DentryA15 = 0.0_DP
        
        ! Now calculate the local matrices on all the elements.
        call computeLocalOptCMatrices (roptcassemblyinfo%Dbas,&
            roptcassemblyinfo%Domega,roptcassemblyinfo%revalElementSet%p_Ddetj,&
            roptcassemblyinfo%indof,roptcassemblyinfo%ncubp,&
            roptcassemblyinfo%p_IelementList(ielset:ielmax),&
            roptcassemblyinfo%DpvelDofs,roptcassemblyinfo%DdvelDofs,&
            roptcassemblyinfo%Dpvel,roptcassemblyinfo%DpvelXderiv,&
            roptcassemblyinfo%DpvelYderiv,roptcassemblyinfo%Ddvel,&
            roptcassemblyinfo%DdvelXderiv,roptcassemblyinfo%DdvelYderiv,&
            DentryA11,DentryA22,DentryA44,DentryA55,&
            DentryA12,DentryA21,DentryA45,DentryA54,&
            DentryA41,DentryA52,DentryA42,DentryA51,&
            DentryA14,DentryA25,DentryA24,DentryA15,&
            roptcassemblyinfo%DlocalDeltaPrimal,roptcassemblyinfo%DlocalDeltaDual,&
            roptcoperator,roptcassemblyinfo)
        
        ! Incorporate the computed local matrices into the global matrix.
        if (associated(p_Da14)) then
          call incorporateLocalMatVecCol (DentryA11,p_Da14,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,1,idof,dweight2,Dtemp)

          if (roptcoperator%ccontrolProjection .eq. 0) then
            call incorporateLocalMatVecCol (DentryA14,p_Da14,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)
          else
            select case (roptcoperator%cconstraintsType)
            case (0)
              call incorporateLocalMatVecCol (DentryA14,p_Da14,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp,&
                                            roptcoperator%dmax1/roptcoperator%dcontrolMultiplier,&
                                            roptcoperator%dmin1/roptcoperator%dcontrolMultiplier)
            case default
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
          end if
        end if

        if (associated(p_Da24)) then
          call incorporateLocalMatVecCol (DentryA22,p_Da24,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,2,idof,dweight2,Dtemp)

          if (roptcoperator%ccontrolProjection .eq. 0) then
            call incorporateLocalMatVecCol (DentryA25,p_Da24,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
          else
            select case (roptcoperator%cconstraintsType)
            case (0)
              call incorporateLocalMatVecCol (DentryA25,p_Da24,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp,&
                                            roptcoperator%dmax2/roptcoperator%dcontrolMultiplier,&
                                            roptcoperator%dmin2/roptcoperator%dcontrolMultiplier)
            case default
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
          end if
        end if


        if (associated(p_Da44)) then
          call incorporateLocalMatVecCol (DentryA41,p_Da44,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA44,p_Da44,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA45,p_Da44,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
        end if


        if (associated(p_Da54)) then
          call incorporateLocalMatVecCol (DentryA52,p_Da54,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,2,idof,dweight2,Dtemp)
                                                      
          call incorporateLocalMatVecCol (DentryA54,p_Da54,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA55,p_Da54,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
        end if
            
      end do
      
      dweight2 = -dweight*0.5_DP/dh
      do idof = 1,roptcassemblyinfo%indof
    
        ! Initialise the element set, compute the basis functions in the
        ! cubature points.
        if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,2),rvelocityVectorPrimal,rvelocityVectorDual,&
              roptcoperator,roptcassemblyinfo,ielset,ielmax,&
            rvelocityVectorPrimalAlt,rvelocityVectorDualAlt)
        else
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,1),rvelocityVectorPrimal,rvelocityVectorDual,&
              roptcoperator,roptcassemblyinfo,ielset,ielmax,&
            rvelocityVectorPrimalAlt,rvelocityVectorDualAlt)
        end if
        
        ! Increase the idof'th entry in the velocity evaluation point by h.
        do iel = 1,ielmax
          roptcassemblyinfo%DdvelDofs(idof,iel,1) = roptcassemblyinfo%DdvelDofs(idof,iel,1) - dh
        end do

        ! Clear the local matrices. If the Newton part is to be calculated,
        ! we must clear everything, otherwise only Dentry.
        DentryA11 = 0.0_DP
        DentryA12 = 0.0_DP
        DentryA21 = 0.0_DP
        DentryA22 = 0.0_DP
        
        DentryA44 = 0.0_DP
        DentryA45 = 0.0_DP
        DentryA54 = 0.0_DP
        DentryA55 = 0.0_DP
        
        DentryA41 = 0.0_DP
        DentryA52 = 0.0_DP
        DentryA42 = 0.0_DP
        DentryA51 = 0.0_DP
        
        DentryA14 = 0.0_DP
        DentryA25 = 0.0_DP
        DentryA24 = 0.0_DP
        DentryA15 = 0.0_DP
        
        ! Now calculate the local matrices on all the elements.
        call computeLocalOptCMatrices (roptcassemblyinfo%Dbas,&
            roptcassemblyinfo%Domega,roptcassemblyinfo%revalElementSet%p_Ddetj,&
            roptcassemblyinfo%indof,roptcassemblyinfo%ncubp,&
            roptcassemblyinfo%p_IelementList(ielset:ielmax),&
            roptcassemblyinfo%DpvelDofs,roptcassemblyinfo%DdvelDofs,&
            roptcassemblyinfo%Dpvel,roptcassemblyinfo%DpvelXderiv,&
            roptcassemblyinfo%DpvelYderiv,roptcassemblyinfo%Ddvel,&
            roptcassemblyinfo%DdvelXderiv,roptcassemblyinfo%DdvelYderiv,&
            DentryA11,DentryA22,DentryA44,DentryA55,&
            DentryA12,DentryA21,DentryA45,DentryA54,&
            DentryA41,DentryA52,DentryA42,DentryA51,&
            DentryA14,DentryA25,DentryA24,DentryA15,&
            roptcassemblyinfo%DlocalDeltaPrimal,roptcassemblyinfo%DlocalDeltaDual,&
            roptcoperator,roptcassemblyinfo)
        
        ! Incorporate the computed local matrices into the global matrix.
        if (associated(p_Da14)) then
          call incorporateLocalMatVecCol (DentryA11,p_Da14,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,1,idof,dweight2,Dtemp)

          if (roptcoperator%ccontrolProjection .eq. 0) then
            call incorporateLocalMatVecCol (DentryA14,p_Da14,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)
          else
            select case (roptcoperator%cconstraintsType)
            case (0)
              call incorporateLocalMatVecCol (DentryA14,p_Da14,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp,&
                                            roptcoperator%dmax1/roptcoperator%dcontrolMultiplier,&
                                            roptcoperator%dmin1/roptcoperator%dcontrolMultiplier)
            case default
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
          end if
        end if

        if (associated(p_Da24)) then
          call incorporateLocalMatVecCol (DentryA22,p_Da24,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,2,idof,dweight2,Dtemp)

          if (roptcoperator%ccontrolProjection .eq. 0) then
            call incorporateLocalMatVecCol (DentryA25,p_Da24,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
          else
            select case (roptcoperator%cconstraintsType)
            case (0)
              call incorporateLocalMatVecCol (DentryA25,p_Da24,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp,&
                                            roptcoperator%dmax2/roptcoperator%dcontrolMultiplier,&
                                            roptcoperator%dmin2/roptcoperator%dcontrolMultiplier)
            case default
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
          end if
        end if


        if (associated(p_Da44)) then
          call incorporateLocalMatVecCol (DentryA41,p_Da44,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA44,p_Da44,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA45,p_Da44,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
        end if


        if (associated(p_Da54)) then
          call incorporateLocalMatVecCol (DentryA52,p_Da54,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,2,idof,dweight2,Dtemp)
                                                      
          call incorporateLocalMatVecCol (DentryA54,p_Da54,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA55,p_Da54,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
        end if
            
      end do
      

      ! Loop through the DOF's. All DOF's must be once increased and once decreased
      ! by h.
      dweight2 = dweight*0.5_DP/dh
      do idof = 1,roptcassemblyinfo%indof
    
        ! Initialise the element set, compute the basis functions in the
        ! cubature points.
        if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,2),rvelocityVectorPrimal,rvelocityVectorDual,&
              roptcoperator,roptcassemblyinfo,ielset,ielmax,&
            rvelocityVectorPrimalAlt,rvelocityVectorDualAlt)
        else
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,1),rvelocityVectorPrimal,rvelocityVectorDual,&
              roptcoperator,roptcassemblyinfo,ielset,ielmax,&
            rvelocityVectorPrimalAlt,rvelocityVectorDualAlt)
        end if
        
        ! Increase the idof'th entry in the velocity evaluation point by h.
        do iel = 1,ielmax
          roptcassemblyinfo%DdvelDofs(idof,iel,2) = roptcassemblyinfo%DdvelDofs(idof,iel,2) + dh
        end do

        ! Clear the local matrices. If the Newton part is to be calculated,
        ! we must clear everything, otherwise only Dentry.
        DentryA11 = 0.0_DP
        DentryA12 = 0.0_DP
        DentryA21 = 0.0_DP
        DentryA22 = 0.0_DP
        
        DentryA44 = 0.0_DP
        DentryA45 = 0.0_DP
        DentryA54 = 0.0_DP
        DentryA55 = 0.0_DP
        
        DentryA41 = 0.0_DP
        DentryA52 = 0.0_DP
        DentryA42 = 0.0_DP
        DentryA51 = 0.0_DP
        
        DentryA14 = 0.0_DP
        DentryA25 = 0.0_DP
        DentryA24 = 0.0_DP
        DentryA15 = 0.0_DP
        
        ! Now calculate the local matrices on all the elements.
        call computeLocalOptCMatrices (roptcassemblyinfo%Dbas,&
            roptcassemblyinfo%Domega,roptcassemblyinfo%revalElementSet%p_Ddetj,&
            roptcassemblyinfo%indof,roptcassemblyinfo%ncubp,&
            roptcassemblyinfo%p_IelementList(ielset:ielmax),&
            roptcassemblyinfo%DpvelDofs,roptcassemblyinfo%DdvelDofs,&
            roptcassemblyinfo%Dpvel,roptcassemblyinfo%DpvelXderiv,&
            roptcassemblyinfo%DpvelYderiv,roptcassemblyinfo%Ddvel,&
            roptcassemblyinfo%DdvelXderiv,roptcassemblyinfo%DdvelYderiv,&
            DentryA11,DentryA22,DentryA44,DentryA55,&
            DentryA12,DentryA21,DentryA45,DentryA54,&
            DentryA41,DentryA52,DentryA42,DentryA51,&
            DentryA14,DentryA25,DentryA24,DentryA15,&
            roptcassemblyinfo%DlocalDeltaPrimal,roptcassemblyinfo%DlocalDeltaDual,&
            roptcoperator,roptcassemblyinfo)
        
        ! Incorporate the computed local matrices into the global matrix.
        if (associated(p_Da15)) then
          call incorporateLocalMatVecCol (DentryA11,p_Da15,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,1,idof,dweight2,Dtemp)

          if (roptcoperator%ccontrolProjection .eq. 0) then
            call incorporateLocalMatVecCol (DentryA14,p_Da15,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)
          else
            select case (roptcoperator%cconstraintsType)
            case (0)
              call incorporateLocalMatVecCol (DentryA14,p_Da15,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp,&
                                            roptcoperator%dmax1/roptcoperator%dcontrolMultiplier,&
                                            roptcoperator%dmin1/roptcoperator%dcontrolMultiplier)
            case default
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
          end if
        end if

        if (associated(p_Da25)) then
          call incorporateLocalMatVecCol (DentryA22,p_Da25,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,2,idof,dweight2,Dtemp)

          if (roptcoperator%ccontrolProjection .eq. 0) then
            call incorporateLocalMatVecCol (DentryA25,p_Da25,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
          else
            select case (roptcoperator%cconstraintsType)
            case (0)
              call incorporateLocalMatVecCol (DentryA25,p_Da25,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp,&
                                            roptcoperator%dmax2/roptcoperator%dcontrolMultiplier,&
                                            roptcoperator%dmin2/roptcoperator%dcontrolMultiplier)
            case default
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
          end if
        end if


        if (associated(p_Da45)) then
          call incorporateLocalMatVecCol (DentryA41,p_Da45,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA44,p_Da45,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA45,p_Da45,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
        end if


        if (associated(p_Da55)) then
          call incorporateLocalMatVecCol (DentryA52,p_Da55,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,2,idof,dweight2,Dtemp)
                                                      
          call incorporateLocalMatVecCol (DentryA54,p_Da55,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA55,p_Da55,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
        end if
            
      end do

      dweight2 = -dweight*0.5_DP/dh
      do idof = 1,roptcassemblyinfo%indof
    
        ! Initialise the element set, compute the basis functions in the
        ! cubature points.
        if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,2),rvelocityVectorPrimal,rvelocityVectorDual,&
              roptcoperator,roptcassemblyinfo,ielset,ielmax,&
            rvelocityVectorPrimalAlt,rvelocityVectorDualAlt)
        else
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,1),rvelocityVectorPrimal,rvelocityVectorDual,&
              roptcoperator,roptcassemblyinfo,ielset,ielmax,&
            rvelocityVectorPrimalAlt,rvelocityVectorDualAlt)
        end if
        
        ! Increase the idof'th entry in the velocity evaluation point by h.
        do iel = 1,ielmax
          roptcassemblyinfo%DdvelDofs(idof,iel,2) = roptcassemblyinfo%DdvelDofs(idof,iel,2) - dh
        end do

        ! Clear the local matrices. If the Newton part is to be calculated,
        ! we must clear everything, otherwise only Dentry.
        DentryA11 = 0.0_DP
        DentryA12 = 0.0_DP
        DentryA21 = 0.0_DP
        DentryA22 = 0.0_DP
        
        DentryA44 = 0.0_DP
        DentryA45 = 0.0_DP
        DentryA54 = 0.0_DP
        DentryA55 = 0.0_DP
        
        DentryA41 = 0.0_DP
        DentryA52 = 0.0_DP
        DentryA42 = 0.0_DP
        DentryA51 = 0.0_DP
        
        DentryA14 = 0.0_DP
        DentryA25 = 0.0_DP
        DentryA24 = 0.0_DP
        DentryA15 = 0.0_DP
        
        ! Now calculate the local matrices on all the elements.
        call computeLocalOptCMatrices (roptcassemblyinfo%Dbas,&
            roptcassemblyinfo%Domega,roptcassemblyinfo%revalElementSet%p_Ddetj,&
            roptcassemblyinfo%indof,roptcassemblyinfo%ncubp,&
            roptcassemblyinfo%p_IelementList(ielset:ielmax),&
            roptcassemblyinfo%DpvelDofs,roptcassemblyinfo%DdvelDofs,&
            roptcassemblyinfo%Dpvel,roptcassemblyinfo%DpvelXderiv,&
            roptcassemblyinfo%DpvelYderiv,roptcassemblyinfo%Ddvel,&
            roptcassemblyinfo%DdvelXderiv,roptcassemblyinfo%DdvelYderiv,&
            DentryA11,DentryA22,DentryA44,DentryA55,&
            DentryA12,DentryA21,DentryA45,DentryA54,&
            DentryA41,DentryA52,DentryA42,DentryA51,&
            DentryA14,DentryA25,DentryA24,DentryA15,&
            roptcassemblyinfo%DlocalDeltaPrimal,roptcassemblyinfo%DlocalDeltaDual,&
            roptcoperator,roptcassemblyinfo)
        
        ! Incorporate the computed local matrices into the global matrix.
        if (associated(p_Da15)) then
          call incorporateLocalMatVecCol (DentryA11,p_Da15,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,1,idof,dweight2,Dtemp)

          if (roptcoperator%ccontrolProjection .eq. 0) then
            call incorporateLocalMatVecCol (DentryA14,p_Da15,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)
          else
            select case (roptcoperator%cconstraintsType)
            case (0)
              ! Note that dcontrolMultiplier is negative, so we have to switch
              ! the min and the max!!!
              call incorporateLocalMatVecCol (DentryA14,p_Da15,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp,&
                                            roptcoperator%dmax1/roptcoperator%dcontrolMultiplier,&
                                            roptcoperator%dmin1/roptcoperator%dcontrolMultiplier)
            case default
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
          end if
        end if

        if (associated(p_Da25)) then
          call incorporateLocalMatVecCol (DentryA22,p_Da25,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,2,idof,dweight2,Dtemp)

          if (roptcoperator%ccontrolProjection .eq. 0) then
            call incorporateLocalMatVecCol (DentryA25,p_Da25,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
          else
            select case (roptcoperator%cconstraintsType)
            case (0)
              call incorporateLocalMatVecCol (DentryA25,p_Da25,roptcassemblyinfo%Kentry,&
                                            roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp,&
                                            roptcoperator%dmax2/roptcoperator%dcontrolMultiplier,&
                                            roptcoperator%dmin2/roptcoperator%dcontrolMultiplier)
            case default
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
          end if
        end if


        if (associated(p_Da45)) then
          call incorporateLocalMatVecCol (DentryA41,p_Da45,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA44,p_Da45,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA45,p_Da45,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
        end if


        if (associated(p_Da55)) then
          call incorporateLocalMatVecCol (DentryA52,p_Da55,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DpvelDofsAlt,2,idof,dweight2,Dtemp)
                                                      
          call incorporateLocalMatVecCol (DentryA54,p_Da55,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,1,idof,dweight2,Dtemp)

          call incorporateLocalMatVecCol (DentryA55,p_Da55,roptcassemblyinfo%Kentry,&
                                          roptcassemblyinfo%DdvelDofsAlt,2,idof,dweight2,Dtemp)
        end if
            
      end do

    end do ! ielset
    !%OMP end do
    
!    ! Release the assembly structure if necessary.
!    if (.not. present(rcubatureInfo)) then
!      call spdiscr_releaseCubStructure(rtempCubatureInfo)
!    end if
   
    deallocate(Dtemp)
    
    ! Release memory
    deallocate(DentryA11)
    deallocate(DentryA22)
    if (allocated(DentryA12)) deallocate(DentryA12)
    if (allocated(DentryA21)) deallocate(DentryA21)
               
    deallocate(DentryA44)
    deallocate(DentryA55)
    if (allocated(DentryA45)) deallocate(DentryA45)
    if (allocated(DentryA54)) deallocate(DentryA54)
               
    deallocate(DentryA41)
    deallocate(DentryA52)
    if (allocated(DentryA42)) deallocate(DentryA42)
    if (allocated(DentryA51)) deallocate(DentryA51)
               
    deallocate(DentryA14)
    deallocate(DentryA25)
    if (allocated(DentryA24)) deallocate(DentryA24)
    if (allocated(DentryA15)) deallocate(DentryA15)
               
  contains
  
    subroutine incorporateLocalMatVecCol (DaLocal,Da,Kentry,Dvelocity,ivelcomp,idof,dweight,Dtemp,&
        dmin,dmax)
    
    ! Incorporate a local matrix into a global one.
    
    ! The local matrix to incorporate.
    real(dp), dimension(:,:,:), intent(in) :: DaLocal
    
    ! The global matrix
    real(dp), dimension(:), intent(inout) :: Da
    
    ! Positions of the local matrix in the global one.
    integer, dimension(:,:,:), intent(in) :: Kentry
    
    ! Velocity vector on the elements
    real(dp), dimension(:,:,:), intent(inout) :: Dvelocity
    
    ! Current velocity component
    integer, intent(in) :: ivelcomp
    
    ! Current DOF. Specifies the target column where to write data to.
    integer, intent(in) :: idof
    
    ! Weight of the local matrix
    real(dp) :: dweight
    
    ! Temporary vector
    real(dp), dimension(:), intent(inout) :: Dtemp

    ! Min/Max value for the velocity. If present, the velocity will be
    ! projected to this range.
    real(dp), intent(in), optional :: dmin,dmax

      ! local variables
      integer :: iel,idofe,jdofe
    
      if (.not. present(dmin)) then
      
        do iel=1,ubound(DaLocal,3)
          ! Apply a local matrix-vector multiplication to get the column
          ! which is to be incorporated into the matrix.
          Dtemp(:) = 0.0_DP
          do idofe = 1,ubound(DaLocal,2)
            do jdofe = 1,ubound(DaLocal,1)
              Dtemp(idofe) = Dtemp(idofe) + DaLocal(jdofe,idofe,iel)*Dvelocity(jdofe,iel,ivelcomp)
            end do
          end do
        
          ! Incorporate the vector to the global matrix.
          do idofe=1,ubound(DaLocal,2)
            Da(Kentry(idof,idofe,iel)) = Da(Kentry(idof,idofe,iel)) + &
                dweight * Dtemp(idofe)
          end do
        end do
        
      else
      
        do iel=1,ubound(DaLocal,3)
          ! Apply a local matrix-vector multiplication to get the column
          ! which is to be incorporated into the matrix.
          Dtemp(:) = 0.0_DP
          do idofe = 1,ubound(DaLocal,2)
            do jdofe = 1,ubound(DaLocal,1)
              Dtemp(idofe) = Dtemp(idofe) + &
                  DaLocal(jdofe,idofe,iel)*min(max(Dvelocity(jdofe,iel,ivelcomp),dmin),dmax)
            end do
          end do
        
          ! Incorporate the vector to the global matrix.
          do idofe=1,ubound(DaLocal,2)
            Da(Kentry(idof,idofe,iel)) = Da(Kentry(idof,idofe,iel)) + &
                dweight * Dtemp(idofe)
          end do
        end do
      
      end if
    
    end subroutine
    
  end subroutine
     
  ! ***************************************************************************

!<subroutine>
  subroutine conv_strdiffOptC2dgetDefect (rvelMatrix,roptcoperator,&
      rvelocityVectorPrimal,rvelocityVectorDual,dweight,rx,rd,rcubatureInfo)
      
!<description>
  ! Calculate the defect of the nonlinear operator:
  !   rd = rd - dweight*A(rvelocityVector)rx
!</description>

!<input>
  ! Template FE matrix specifying the connectivity in the FE space
  ! of the primal/dual velocity.
  type(t_matrixScalar), intent(in) :: rvelMatrix

  ! Structure defining the operator to set up.
  type(t_optcoperator), intent(in) :: roptcoperator
  
  ! Velocity vector for the nonlinearity.
  ! The first blocks 1/2 in this vector define the evaluation
  ! point (primal velocity).
  type(t_vectorBlock), intent(in) :: rvelocityVectorPrimal

  ! Velocity vector for the nonlinearity.
  ! The first blocks 4/5 in this vector define the evaluation
  ! point (dual velocity).
  type(t_vectorBlock), intent(in) :: rvelocityVectorDual
  
  ! Solution vector, to be multiplied with the matrix.
  type(t_vectorBlock), intent(in) :: rx
  
  ! Weight for the solution
  real(DP), intent(in) :: dweight

  ! OPTIONAL: A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), target :: rcubatureInfo
!</input>
  
!<inputoutput>
  ! On entry: RHS. On exit: Defect.
  type(t_vectorBlock), intent(inout) :: rd
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: ielset,ielmax
  
  ! Local matrices
  real(DP), dimension(:,:,:), allocatable :: DentryA11
  real(DP), dimension(:,:,:), allocatable :: DentryA12
  real(DP), dimension(:,:,:), allocatable :: DentryA21
  real(DP), dimension(:,:,:), allocatable :: DentryA22
  
  real(DP), dimension(:,:,:), allocatable :: DentryA44
  real(DP), dimension(:,:,:), allocatable :: DentryA45
  real(DP), dimension(:,:,:), allocatable :: DentryA54
  real(DP), dimension(:,:,:), allocatable :: DentryA55
  
  real(DP), dimension(:,:,:), allocatable :: DentryA41
  real(DP), dimension(:,:,:), allocatable :: DentryA52
  real(DP), dimension(:,:,:), allocatable :: DentryA42
  real(DP), dimension(:,:,:), allocatable :: DentryA51
  
  real(DP), dimension(:,:,:), allocatable :: DentryA14
  real(DP), dimension(:,:,:), allocatable :: DentryA25
  real(DP), dimension(:,:,:), allocatable :: DentryA24
  real(DP), dimension(:,:,:), allocatable :: DentryA15
  
  ! Temp vector
  type(t_vectorBlock) :: rvectorBlock
  
  ! Pointer to the subvectors
  real(dp), dimension(:), pointer :: p_Dx1,p_Dx2,p_Dx4,p_Dx5
  real(dp), dimension(:), pointer :: p_Dd1,p_Dd2,p_Dd4,p_Dd5
  
  ! Assembly status structure
  type(t_optcassemblyinfo) :: roptcassemblyinfo

  ! Cubature formula and cubature info structure
  type(t_scalarCubatureInfo), target :: rtempCubatureInfo
  type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo
  integer(I32) :: ccubature

    if (roptcoperator%dnu .eq. 0.0_DP) then
      print *,'SD: NU=0 not allowed! Set dbeta=0 to prevent Stokes operator'// &
              ' from being build!'
      call sys_halt()
    end if
    
    ! Initialise the asembly of the local matrices.
    call conv_strdiffOptC2dinitasm (rvelocityVectorPrimal%p_rblockDiscr%RspatialDiscr(1),&
        1,roptcassemblyinfo)

    ! Calculate the maximum norm of the actual velocity field
    ! Round up the norm to 1D-8 if it's too small...
    call lsysbl_deriveSubvector(rvelocityVectorPrimal,rvectorBlock,1,2,.true.)
    call lsysbl_getVectorMagnitude (rvectorBlock,dumax=roptcassemblyinfo%dumax)
    call lsysbl_releaseVector (rvectorBlock)
  
    if (roptcassemblyinfo%dumax .lt. 1E-8_DP) roptcassemblyinfo%dumax=1E-8_DP

    ! Allocate memory for the local matrices
    allocate(DentryA11(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA12(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA21(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA22(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))

    allocate(DentryA44(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA45(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA54(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA55(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))

    allocate(DentryA41(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA52(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA42(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA51(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))

    allocate(DentryA14(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA25(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA24(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA15(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
  
    ! Initialise the array with the local Delta values for the stabilisation
    call lalg_clearVectorDble (roptcassemblyinfo%DlocalDeltaPrimal)
    call lalg_clearVectorDble (roptcassemblyinfo%DlocalDeltaDual)
    
    ! Get pointers to the subvectors
    call lsyssc_getbase_double (rx%RvectorBlock(1),p_Dx1)
    call lsyssc_getbase_double (rx%RvectorBlock(2),p_Dx2)
    call lsyssc_getbase_double (rx%RvectorBlock(4),p_Dx4)
    call lsyssc_getbase_double (rx%RvectorBlock(5),p_Dx5)

    call lsyssc_getbase_double (rd%RvectorBlock(1),p_Dd1)
    call lsyssc_getbase_double (rd%RvectorBlock(2),p_Dd2)
    call lsyssc_getbase_double (rd%RvectorBlock(4),p_Dd4)
    call lsyssc_getbase_double (rd%RvectorBlock(5),p_Dd5)

!    ! Do we have an assembly structure?
!    ! If we do not have it, create a cubature info structure that
!    ! defines how to do the assembly.
!    if (.not. present(rcubatureInfo)) then
!      call spdiscr_createDefCubStructure(p_rdiscretisation,&
!          rtempCubatureInfo,CUB_GEN_DEPR_EVAL)
!      p_rcubatureInfo => rtempCubatureInfo
!    else
      p_rcubatureInfo => rcubatureInfo
!    end if

    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so BILF_NELEMSIM local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !%OMP do SCHEDULE(dynamic,1)
    do ielset = 1, size(roptcassemblyinfo%p_IelementList), BILF_NELEMSIM

      ! We always handle BILF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
      ! elements simultaneously.
      
      ielmax = min(size(roptcassemblyinfo%p_IelementList),ielset-1+BILF_NELEMSIM)
    
      ! Initialise the element set, compute the basis functions in the
      ! cubature points.
      call conv_strdiffOptC2dinitelemset (rvelMatrix,rvelMatrix,&
          rvelocityVectorPrimal,rvelocityVectorDual,roptcoperator,&
          roptcassemblyinfo,ielset,ielmax)

      ! Clear the local matrices. If the Newton part is to be calculated,
      ! we must clear everything, otherwise only Dentry.
      DentryA11 = 0.0_DP
      DentryA12 = 0.0_DP
      DentryA21 = 0.0_DP
      DentryA22 = 0.0_DP
      
      DentryA44 = 0.0_DP
      DentryA45 = 0.0_DP
      DentryA54 = 0.0_DP
      DentryA55 = 0.0_DP
      
      DentryA41 = 0.0_DP
      DentryA52 = 0.0_DP
      DentryA42 = 0.0_DP
      DentryA51 = 0.0_DP
      
      DentryA14 = 0.0_DP
      DentryA25 = 0.0_DP
      DentryA24 = 0.0_DP
      DentryA15 = 0.0_DP
      
      ! Now calculate the local matrices on all the elements.
      call computeLocalOptCMatrices (roptcassemblyinfo%Dbas,&
          roptcassemblyinfo%Domega,roptcassemblyinfo%revalElementSet%p_Ddetj,&
          roptcassemblyinfo%indof,roptcassemblyinfo%ncubp,&
          roptcassemblyinfo%p_IelementList(ielset:ielmax),&
          roptcassemblyinfo%DpvelDofs,roptcassemblyinfo%DdvelDofs,&
          roptcassemblyinfo%Dpvel,roptcassemblyinfo%DpvelXderiv,&
          roptcassemblyinfo%DpvelYderiv,roptcassemblyinfo%Ddvel,&
          roptcassemblyinfo%DdvelXderiv,roptcassemblyinfo%DdvelYderiv,&
          DentryA11,DentryA22,DentryA44,DentryA55,&
          DentryA12,DentryA21,DentryA45,DentryA54,&
          DentryA41,DentryA52,DentryA42,DentryA51,&
          DentryA14,DentryA25,DentryA24,DentryA15,&
          roptcassemblyinfo%DlocalDeltaPrimal,roptcassemblyinfo%DlocalDeltaDual,&
          roptcoperator,roptcassemblyinfo)
      
      ! Perform a local matrix-vector multiplication with the local matrices.
      ! Incorporate the computed local matrices into the global matrix.
      call calcDefect (DentryA11,roptcassemblyinfo%Idofs,p_Dx1,p_Dd1,dweight)
      call calcDefect (DentryA12,roptcassemblyinfo%Idofs,p_Dx2,p_Dd1,dweight)
      call calcDefect (DentryA21,roptcassemblyinfo%Idofs,p_Dx1,p_Dd2,dweight)
      call calcDefect (DentryA22,roptcassemblyinfo%Idofs,p_Dx2,p_Dd2,dweight)
      
      call calcDefect (DentryA41,roptcassemblyinfo%Idofs,p_Dx1,p_Dd4,dweight)
      call calcDefect (DentryA42,roptcassemblyinfo%Idofs,p_Dx2,p_Dd4,dweight)
      call calcDefect (DentryA51,roptcassemblyinfo%Idofs,p_Dx1,p_Dd5,dweight)
      call calcDefect (DentryA52,roptcassemblyinfo%Idofs,p_Dx2,p_Dd5,dweight)

      call calcDefect (DentryA44,roptcassemblyinfo%Idofs,p_Dx4,p_Dd4,dweight)
      call calcDefect (DentryA45,roptcassemblyinfo%Idofs,p_Dx5,p_Dd4,dweight)
      call calcDefect (DentryA54,roptcassemblyinfo%Idofs,p_Dx4,p_Dd5,dweight)
      call calcDefect (DentryA55,roptcassemblyinfo%Idofs,p_Dx5,p_Dd5,dweight)

      if (roptcoperator%ccontrolProjection .eq. 0) then
        call calcDefect (DentryA14,roptcassemblyinfo%Idofs,p_Dx4,p_Dd1,dweight)
        call calcDefect (DentryA24,roptcassemblyinfo%Idofs,p_Dx4,p_Dd2,dweight)
        call calcDefect (DentryA15,roptcassemblyinfo%Idofs,p_Dx5,p_Dd1,dweight)
        call calcDefect (DentryA25,roptcassemblyinfo%Idofs,p_Dx5,p_Dd2,dweight)
      else
        call calcDefect (DentryA14,roptcassemblyinfo%Idofs,p_Dx4,p_Dd1,dweight,&
            roptcoperator%dmax1/roptcoperator%dcontrolMultiplier,&
            roptcoperator%dmin1/roptcoperator%dcontrolMultiplier)
        call calcDefect (DentryA24,roptcassemblyinfo%Idofs,p_Dx4,p_Dd2,dweight,&
            roptcoperator%dmax1/roptcoperator%dcontrolMultiplier,&
            roptcoperator%dmin1/roptcoperator%dcontrolMultiplier)
        call calcDefect (DentryA15,roptcassemblyinfo%Idofs,p_Dx5,p_Dd1,dweight,&
            roptcoperator%dmax2/roptcoperator%dcontrolMultiplier,&
            roptcoperator%dmin2/roptcoperator%dcontrolMultiplier)
        call calcDefect (DentryA25,roptcassemblyinfo%Idofs,p_Dx5,p_Dd2,dweight,&
            roptcoperator%dmax2/roptcoperator%dcontrolMultiplier,&
            roptcoperator%dmin2/roptcoperator%dcontrolMultiplier)
      end if
      
    end do ! ielset
    !%OMP end do
    
!    ! Release the assembly structure if necessary.
!    if (.not. present(rcubatureInfo)) then
      call spdiscr_releaseCubStructure(rtempCubatureInfo)
!    end if
    
    ! Release memory
    deallocate(DentryA11)
    deallocate(DentryA22)
    deallocate(DentryA12)
    deallocate(DentryA21)
               
    deallocate(DentryA44)
    deallocate(DentryA55)
    deallocate(DentryA45)
    deallocate(DentryA54)
               
    deallocate(DentryA41)
    deallocate(DentryA52)
    deallocate(DentryA42)
    deallocate(DentryA51)
               
    deallocate(DentryA14)
    deallocate(DentryA25)
    deallocate(DentryA24)
    deallocate(DentryA15)
               
  contains
    
    subroutine calcDefect (DaLocal,KcolLocal,Dx,Dd,dweight,dmin,dmax)
    
    ! Calculate the local defect:
    !   Dd = Dd - DaLocal*Dx
    
    ! The local matrix to incorporate.
    real(dp), dimension(:,:,:), intent(in) :: DaLocal
    
    ! Global column numbers of the local matrix columns on each element.
    integer, dimension(:,:), intent(in) :: KcolLocal
    
    ! Solution vector
    real(dp), dimension(:), intent(in) :: Dx
    
    ! RHS vector; receives the local defect.
    real(dp), dimension(:), intent(inout) :: Dd
    
    ! Weight for the solution vectot
    real(dp), intent(in) :: dweight
    
    ! Min/Max value for the velocity. If present, the velocity will be
    ! projected to this range.
    real(dp), intent(in), optional :: dmin,dmax
    
      ! local variables
      integer :: iel,idofe,jdofe
      real(dp) :: dval
    
      if (.not. present(dmin)) then
        do iel=1,ubound(DaLocal,3)
          do idofe=1,ubound(DaLocal,1)
            dval = 0.0_DP
            ! DaLocal*Dx. Note that DaLocal is saved transposed for quicker access,
            ! so jdofe in the first index describes the column!
            do jdofe=1,ubound(DaLocal,2)
              dval = dval + DaLocal(jdofe,idofe,iel) * Dx(KcolLocal(jdofe,iel))
            end do
            Dd(KcolLocal(idofe,iel)) = Dd(KcolLocal(idofe,iel)) - dweight*dval
          end do
        end do
      else
        do iel=1,ubound(DaLocal,3)
          do idofe=1,ubound(DaLocal,1)
            dval = 0.0_DP
            ! DaLocal*Dx. Note that DaLocal is saved transposed for quicker access,
            ! so jdofe in the first index describes the column!
            do jdofe=1,ubound(DaLocal,2)
              dval = dval + DaLocal(jdofe,idofe,iel) * &
                  min(max(Dx(KcolLocal(jdofe,iel)),dmin),dmax)
            end do
            Dd(KcolLocal(idofe,iel)) = Dd(KcolLocal(idofe,iel)) - dweight*dval
          end do
        end do
      end if
    
    end subroutine
    
  end subroutine
     
  ! ----------------------------------------------------------------------

  subroutine getLocalDeltaQuad (Du,rtriangulation,Ielements,&
      dumax,dupsam,dnu,DlocalDelta)

  ! This routine calculates a local ddelta=DELTA_T for a finite element
  ! T=IEL. This can be used by the streamline diffusion stabilisation
  ! technique as a multiplier of the (local) bilinear form.
  
  ! Velocity DOF's on all elements.
  ! Du(#dofs per element,#elements,ndim2d)
  real(DP), dimension(:,:,:), intent(IN) :: Du
  
  ! Underlying triangulation
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! List of elements to process
  integer, dimension(:), intent(in) :: Ielements
  
  ! Maximum norm of velocity in the domain:
  ! duMaxR = ||u||_Omega
  real(DP), intent(IN) :: dumax
  
  ! Viscosity parameter
  real(DP), intent(IN) :: dnu
  
  ! user defined parameter for configuring the streamline diffusion.
  ! < 0: Simple calculation of ddelta, using
  !      ddelta = |UPSAM| * h_T.
  ! > 0: usually UPSAM = 0.1 .. 2; Samarskji-like calculation of ddelta using:
  !      ddelta = UPSAM * h_t/||u||_T * 2*Re_T/(1+Re_T)
  real(DP), intent(IN) :: dupsam
  
  ! local local delta on all elements in the list
  real(DP), dimension(:), intent(out) :: DlocalDelta

    ! local variables
    real(DP) :: dlocalH,du1,du2,dunorm,dreLocal,dnuR,dumaxR
    integer :: iel
    integer :: idof
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    call storage_getbase_int2d(rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    call storage_getbase_double2d(rtriangulation%h_DvertexCoords,p_DvertexCoords)

    ! We later need the reciprocals of dumax and dre
    dnuR = 1.0_DP/dnu
    dumaxR = 1.0_DP/dumax

    ! Loop over all elements
    do iel = 1,size(Ielements)
    
      ! Calculate the mean velocity on the current element
      ! and using this the norm of the local velocity.
      du1=0.0_dp
      du2=0.0_dp
      do idof=1,ubound(Du,1)
        du1=du1+Du(idof,iel,1)
        du2=du2+Du(idof,iel,2)
      end do
    
      dunorm = sqrt(du1**2+du2**2) / dble(ubound(Du,1))
    
      ! If the norm of the velocity is small, we choose ddelta = 0,
      ! which results in central difference in the streamline diffusion
      ! matrix assembling:

      if (dunorm .le. 1D-8) then
      
        DlocalDelta(iel) = 0.0_DP

      else

        ! u_T defines the "slope" of the velocity through
        ! the element T. At next, calculate the local mesh width
        ! dlocalH = h = h_T on our element T=IEL:

        call getLocalMeshWidthQuad (dlocalH,dunorm, du1, du2, iel, &
            p_IverticesAtElement,p_DvertexCoords)

        ! Calculate ddelta... (cf. p. 121 in Turek's CFD book)

        if (dupsam .lt. 0.0_DP) then

          ! For UPSAM<0, we use simple calculation of ddelta:
        
          DlocalDelta(iel) = abs(dupsam)*dlocalH
          
        else
        
          ! For UPSAM >= 0, we use standard Samarskji-like calculation
          ! of ddelta. At first calculate the local Reynolds number
          ! RELOC = Re_T = ||u||_T * h_T / NU
          
          dreLocal = dunorm*dlocalH*dnuR
          
          ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)
          
          DlocalDelta(iel) = dupsam * dlocalH*duMaxR * 2.0_DP*(dreLocal/(1.0_DP+dreLocal))
          
        end if ! (UPSAM.LT.0.0)
        
      end if ! (dunorm.LE.1D-8)

    end do

  end subroutine

  ! ----------------------------------------------------------------------

  pure subroutine getLocalMeshWidthQuad (dlocalH, dunorm,  XBETA1, &
                      XBETA2, JEL,Kvert,Dcorvg)
  
  ! Determine the local mesh width for an element JEL of a
  ! triangulation.
  
  ! Element where the local h should be calculated
  integer, intent(IN)               :: JEL
  
  integer, dimension(TRIA_MAXNVE2D,*), intent(IN) :: Kvert
  real(DP), dimension(NDIM2D,*), intent(IN)          :: Dcorvg
  
  ! norm ||u||_T = mean velocity through element T=JEL
  real(DP), intent(IN)  :: dunorm
  
  ! mean velocity u_T = (xbeta1,xbeta2) through element T=JEL
  real(DP), intent(IN)  :: XBETA1, XBETA2
  
  ! local mesh width
  real(DP), intent(OUT) :: dlocalH
  
  ! local variables
  real(DP) :: dlambda
  integer :: NECK1,NECK2,NECK3,NECK4
  real(DP) :: X1,Y1,X2,Y2,X3,Y3,X4,Y4
  real(DP) :: dalphaMax, dalpha

    ! Fetch the numbers of the four corners of element JEL

    neck1=Kvert(1,JEL)
    neck2=Kvert(2,JEL)
    neck3=Kvert(3,JEL)
    neck4=Kvert(4,JEL)

    ! Fetch the coordinates of these corners

    x1=Dcorvg(1,neck1)
    y1=Dcorvg(2,neck1)
    x2=Dcorvg(1,neck2)
    y2=Dcorvg(2,neck2)
    x3=Dcorvg(1,neck3)
    y3=Dcorvg(2,neck3)
    x4=Dcorvg(1,neck4)
    y4=Dcorvg(2,neck4)

    ! Scale: (deactivated)

    !  dsp=max(xbeta1,xbeta2)

    !  xbeta1=xbeta1
    !  xbeta2=xbeta2

    dalphaMax=0.0_DP
    
    ! In the next step, we calculate the 'maximum possible mesh with
    ! in direction of the flow'; this is the maximum possible length
    ! that a particle can cross in the current element.
    ! The picture in mind is the following:
    !
    !          G3
    !   +-------------X-------+
    !   |            /        |
    !   |           /         |
    !   |          /          |
    !   |         /           |
    !   |        /            |
    ! G4|       /             | G2
    !   |      ^ (beta1,beta2)|
    !   |     /               |
    !   |    /                |
    !   |   /                 |
    !   |  /                  |
    !   | /                   |
    !   |/                    |
    !   O---------------------+
    !            G1
    !
    ! The vector (beta1,beta2) gives the direction of the flow.
    ! A particle starting in point O and moves at most up to point X.
    ! The length of the line (O,X) is the local mesh with h.
    !
    ! Loop through the four corners of element JEL and check
    ! of a line with slope BETA=(xbeta1,xbeta2) starting in this
    ! corner really intersects with one of the edges of the element.
    ! Remark that we only have to check the two opposite edges
    ! to the current corner!

    ! -----------------------------------------------------------------
    ! Check the first corner:

    call intersectLines2D(X1,Y1,dalpha,XBETA1,XBETA2, &
                X3,Y3,dlambda,X2,Y2)
    dalphaMax=max(dalpha,dalphaMax)

    call intersectLines2D(X1,Y1,dalpha,XBETA1,XBETA2, &
                X3,Y3,dlambda,X4,Y4)
    dalphaMax=max(dalpha,dalphaMax)
    
    ! -----------------------------------------------------------------
    ! The second one...
    
    call intersectLines2D(X2,Y2,dalpha,XBETA1,XBETA2, &
                X4,Y4,dlambda,X1,Y1)
    dalphaMax=max(dalpha,dalphaMax)

    call intersectLines2D(X2,Y2,dalpha,XBETA1,XBETA2, &
                X4,Y4,dlambda,X3,Y3)
    dalphaMax=max(dalpha,dalphaMax)
    
    ! -----------------------------------------------------------------
    ! The third one...
    
    call intersectLines2D(X3,Y3,dalpha,XBETA1,XBETA2, &
                X1,Y1,dlambda,X2,Y2)
    dalphaMax=max(dalpha,dalphaMax)

    call intersectLines2D(X3,Y3,dalpha,XBETA1,XBETA2, &
                X1,Y1,dlambda,X4,Y4)
    dalphaMax=max(dalpha,dalphaMax)
    
    ! -----------------------------------------------------------------
    ! And the fourth=last one...
    
    call intersectLines2D(X4,Y4,dalpha,XBETA1,XBETA2, &
                X2,Y2,dlambda,X1,Y1)
    dalphaMax=max(dalpha,dalphaMax)

    call intersectLines2D(X4,Y4,dalpha,XBETA1,XBETA2, &
                X2,Y2,dlambda,X3,Y3)
    dalphaMax=max(dalpha,dalphaMax)

    ! -----------------------------------------------------------------
    ! finally determine the local h=h_T
    !
    ! dalphaMax is the maximum alpha, normalised as 'parameter value',
    ! i.e. dalphaMax=1.0 corresponds to a vector 1.0*(dbeta1,dbeta2).
    ! We multiply with dunorm=|(dbeta1,dbeta2)| to get the actual length
    ! of the vector which can be placed inside of the element.
    !
    ! Furthermore, we multiply with an additional weight 4. (why ?!?)

    dlocalH=dalphaMax*4.0_DP*dunorm

  end subroutine
  
  ! ----------------------------------------------------------------------

  pure subroutine intersectLines2D (XO,YO,dalpha,BETA1,BETA2, &
                      XA,YA,dlambda,XB,YB)

  ! Intersect two lines in R^2

  ! Origin of line 1
  real(DP), intent(IN) :: XO,YO
  
  ! Direction of line 1
  real(DP), intent(IN) :: BETA1,BETA2
  
  ! One point on the second line
  real(DP), intent(IN) :: XA,YA
  
  ! Another point on the second line
  real(DP), intent(IN) :: XB,YB
  
  ! Parameter value of the intersection point on line 1.
  ! =0.0, if there is no intersection point
  real(DP), intent(OUT) :: dalpha
  
  real(DP), intent(OUT) :: dlambda
  
  ! local variables
  double precision :: dsp

    ! Scalar product of the line (xa,ya)->(xb,yb) with the
    ! counterclockwise normal n1 of (beta1,beta2)
    dsp=BETA2*(XB-XA)-BETA1*(YB-YA)
    
    if (dsp.eq.0.0_DP) then
    
      ! beta and the vector are parallel
      dalpha=0.0_DP
      
    else

      ! Scalar product of (beta1,beta2) with the (inner) normal vector n2
      ! of the line (xo,yo)->(xa,ya).
      dlambda=(BETA1*(YA-YO)-BETA2*(XA-XO))/dsp

      !                    (xb,yb)
      !   +-----------------+
      !   |                 |
      !   |                 |
      !   ^ n2              |
      !   !                 |
      !   !  (beta1,beta2)  |    (beta1,beta2)
      !   !    ^            |    ^
      !   !   /  ^__ n1     |   /
      !   !  /      \__     |  /
      !   ! /          \__  | /
      !   !/              \_|/
      !   +-----------------+
      ! (xo,yo)            (xa,ya)
      !
      ! (What is this? Documentation incomplete. Has someone a good
      ! reference?)

      ! is the intersection point inside of the element?
      if ((dlambda.ge.-1E-1_DP).and.(dlambda.le.1.11E0_DP)) then
        if (BETA1 .ne. 0.0_DP) then
          dalpha=((XA-XO)+dlambda*(XB-XA))/BETA1
        else
          if (BETA2 .ne. 0.0_DP) then
            dalpha=((YA-YO)+dlambda*(YB-YA))/BETA2
          else
            dalpha=0.0_DP
          end if
        end if
      else
        dalpha=0.0_DP
      end if
      
    end if

  end subroutine

end module
