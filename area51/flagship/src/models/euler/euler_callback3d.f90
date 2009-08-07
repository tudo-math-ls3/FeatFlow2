!##############################################################################
!# ****************************************************************************
!# <name> euler_callback3d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the compressible Euler/Navier-Stokes equations in 3D.
!#
!# The following callback functions are available:
!#
!# 1.) euler_calcFluxGalerkin3d
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#
!# 2.) euler_calcFluxGalerkinNoBdr3d
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#        without assembling the symmetric boundary contribution
!#
!# 3.) euler_calcFluxScalarDiss3d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting scalar artificial viscosities
!#
!# 4.) euler_calcFluxDSplitScalarDiss3d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting scalar artificial viscosities based on
!#        dimensional splitting approach
!#
!# 5.) euler_calcFluxTensorDiss3d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting tensorial artificial viscosities
!#
!# 6.) euler_calcFluxDSplitTensorDiss3d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting tensorial artificial viscosities based on
!#        dimensional splitting approach
!#
!# 7.) euler_calcFluxRusanov3d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting the Rusanov artificial diffusion
!#
!# 8.) euler_calcFluxDSplitRusanov3d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting the Rusanov artificial diffusion
!#
!# 9.) euler_calcMatrixDiagonalDiag3d
!#     -> Computes local matrix for diagonal entry
!#
!# 10.) euler_calcMatrixDiagonal3d
!#      -> Computes local matrix for diagonal entry
!#
!# 11.) euler_calcMatrixGalerkinDiag3d
!#      -> Computes local matrices for standard Galerkin scheme
!#
!# 12.) euler_calcMatrixGalerkin3d
!#      -> Computes local matrices for standard Galerkin scheme
!#
!# 13.) euler_calcMatrixScalarDissDiag3d
!#      -> Computes local matrices for low-order discretization
!#         adopting scalar artificial viscosities
!#
!# 14.) euler_calcMatrixScalarDiss3d
!#      -> Computes local matrices for low-order discretization
!#         adopting scalar artificial viscosities
!#
!# 15.) euler_calcMatrixTensorDissDiag3d
!#      -> Computes local matrices for low-order discretization
!#         adopting tensorial artificial viscosities
!#
!# 16.) euler_calcMatrixTensorDiss3d
!#      -> Computes local matrices for low-order discretization
!#         adopting tensorial artificial viscosities
!#
!# 17.) euler_calcMatrixRusanovDiag3d
!#      -> Computes local matrices for low-order discretization
!#         adopting the Rusanov artificial viscosities
!#
!# 18.) euler_calcMatrixRusanov3d
!#      -> Computes local matrices for low-order discretization
!#         adopting the Rusanov flux artificial viscosities
!#
!# 19.) euler_calcCharacteristics3d
!#      -> Computes characteristic variables in 3D
!#
!# 20.) euler_calcBoundaryvalues3d
!#      -> Computes the boundary values for a given node
!#
!# 21.) euler_hadaptCallbackScalar3d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 3D, whereby the vector is stored in interleave format
!#
!# 22.) euler_hadaptCallbackBlock3d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 3D, whereby the vector is stored in block format
!#
!# </purpose>
!##############################################################################

module euler_callback3d

  use boundaryfilter
  use collection
  use euler_basic
  use flagship_callback
  use fsystem
  use genoutput
  use graph
  use groupfemsystem
  use hadaptaux
  use linearsystemblock
  use linearsystemscalar
  use problem
  use solveraux
  use storage
  use thermodynamics

  implicit none

  private
  public :: euler_calcFluxGalerkin3d
  public :: euler_calcFluxGalerkinNoBdr3d
  public :: euler_calcFluxScalarDiss3d
  public :: euler_calcFluxDSplitScalarDiss3d
  public :: euler_calcFluxTensorDiss3d
  public :: euler_calcFluxDSplitTensorDiss3d
  public :: euler_calcFluxRusanov3d
  public :: euler_calcFluxDSplitRusanov3d
  public :: euler_calcMatrixDiagonalDiag3d
  public :: euler_calcMatrixDiagonal3d
  public :: euler_calcMatrixGalerkinDiag3d
  public :: euler_calcMatrixGalerkin3d
  public :: euler_calcMatrixScalarDissDiag3d
  public :: euler_calcMatrixScalarDiss3d
  public :: euler_calcMatrixTensorDissDiag3d
  public :: euler_calcMatrixTensorDiss3d
  public :: euler_calcMatrixRusanovDiag3d
  public :: euler_calcMatrixRusanov3d
  public :: euler_calcCharacteristics3d
  public :: euler_calcBoundaryvalues3d
  public :: euler_hadaptCallbackScalar3d
  public :: euler_hadaptCallbackBlock3d
  
contains
  
  !*****************************************************************************
  
!<subroutine>

  pure subroutine euler_calcFluxGalerkin3d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the standard Galerkin 
    ! discretization in 3D.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling coefficient
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! inviscid fluxes
    real(DP), dimension(:), intent(out) :: F_ij, F_ji
!</output>
!</subroutine>

  end subroutine euler_calcFluxGalerkin3d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxGalerkinNoBdr3d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the TVD
    ! discretization in 3D. The symmetric boundary contributions
    ! are neglected and incorporated in the antidiffusive flux.
    ! Hence, this is simply the standard Galerkin flux for the
    ! skew-symmetric internal contributions.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling coefficient
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! inviscid fluxes
    real(DP), dimension(:), intent(out) :: F_ij, F_ji
!</output>
!</subroutine>

  end subroutine euler_calcFluxGalerkinNoBdr3d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxScalarDiss3d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 3D using scalar dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling coefficient
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! inviscid fluxes
    real(DP), dimension(:), intent(out) :: F_ij, F_ji
!</output>
!</subroutine>

  end subroutine euler_calcFluxScalarDiss3d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxDSplitScalarDiss3d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 3D using scalar dissipation,
    ! whereby dimensional splitting is employed.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling coefficient
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! inviscid fluxes
    real(DP), dimension(:), intent(out) :: F_ij, F_ji
!</output>
!</subroutine>

  end subroutine euler_calcFluxDSplitScalarDiss3d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxTensorDiss3d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 3D using tensorial dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling coefficient
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! inviscid fluxes
    real(DP), dimension(:), intent(out) :: F_ij, F_ji
!</output>
!</subroutine>

  end subroutine euler_calcFluxTensorDiss3d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxDSplitTensorDiss3d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 3D using tensorial dissipation,
    ! whereby dimensional splitting is employed.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling coefficient
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! inviscid fluxes
    real(DP), dimension(:), intent(out) :: F_ij, F_ji
!</output>
!</subroutine>

  end subroutine euler_calcFluxDSplitTensorDiss3d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxRusanov3d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 3D using the Rusanov dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling coefficient
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! inviscid fluxes
    real(DP), dimension(:), intent(out) :: F_ij, F_ji
!</output>
!</subroutine>

  end subroutine euler_calcFluxRusanov3d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxDSplitRusanov3d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 3D using the Rusanov dissipation,
    ! whereby dimensional splitting is employed.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling coefficient
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! inviscid fluxes
    real(DP), dimension(:), intent(out) :: F_ij, F_ji
!</output>
!</subroutine>

  end subroutine euler_calcFluxDSplitRusanov3d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixDiagonalDiag3d(U_i, C_ii, i, dscale, K_ii)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! for the diagonal block of the global operator in 3D
!</description>

!<input>
    ! local solution at node I
    real(DP), dimension(:), intent(in) :: U_i

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ii

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! row number
    integer, intent(in) :: i
!</input>

!<output>
    ! local Jacobian matrix
    real(DP), dimension(:), intent(out) :: K_ii
!</output>
!</subroutine>

  end subroutine euler_calcMatrixDiagonalDiag3d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixDiagonal3d(U_i, C_ii, i, dscale, K_ii)

!<description>
    ! This subroutine computes the Galerkin matrices for the diagonal
    ! block of the global operator in 3D
!</description>

!<input>
    ! local solution at node I
    real(DP), dimension(:), intent(in) :: U_i

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ii

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! row number
    integer, intent(in) :: i
!</input>

!<output>
    ! local Jacobian matrix
    real(DP), dimension(:), intent(out) :: K_ii
!</output>
!</subroutine>

  end subroutine euler_calcMatrixDiagonal3d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixGalerkinDiag3d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices in 3D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(out) :: K_ij,K_ji,D_ij
!</output>
!</subroutine>

  end subroutine euler_calcMatrixGalerkinDiag3d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixGalerkin3d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
     ! This subroutine computes the Galerkin matrices in 3D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(out) :: K_ij,K_ji,D_ij
!</output>
!</subroutine>

  end subroutine euler_calcMatrixGalerkin3d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixScalarDissDiag3d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies scalar artificial viscosities in 3D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij,C_ji
    
    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(out) :: K_ij,K_ji,D_ij
!</output>
!</subroutine>

  end subroutine euler_calcMatrixScalarDissDiag3d

!*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixScalarDiss3d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies scalar artificial viscosities in 3D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(out) :: K_ij,K_ji,D_ij
!</output>
!</subroutine>

  end subroutine euler_calcMatrixScalarDiss3d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixTensorDissDiag3d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

    !<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies tensorial artificial viscosities in 3D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(out) :: K_ij,K_ji,D_ij
!</output>
!</subroutine>

  end subroutine euler_calcMatrixTensorDissDiag3d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixTensorDiss3d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies tensorial artificial viscosities in 3D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(out) :: K_ij,K_ji,D_ij
!</output>
!</subroutine>

  end subroutine euler_calcMatrixTensorDiss3d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixRusanovDiag3d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 3D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(out) :: K_ij,K_ji,D_ij
!</output>
!</subroutine>

  end subroutine euler_calcMatrixRusanovDiag3d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixRusanov3d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 3D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(out) :: K_ij,K_ji,D_ij
!</output>
!</subroutine>

  end subroutine euler_calcMatrixRusanov3d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcCharacteristics3d(&
      U_i, U_j, Dweight, W_ij, Lbd_ij, R_ij, L_ij)

!<description>
    ! This subroutine computes the characteristic variables in 3D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! weighting vector
    real(DP), dimension(:), intent(in) :: Dweight
!</input>

!<output>
    ! vector of characteristic variables
    real(DP), dimension(:), intent(out), optional :: W_ij

    ! OPTIONAL: diagonal matrix of eigenvalues
    real(DP), dimension(:), intent(out), optional :: Lbd_ij

    ! OPTIONAL: transformation matrix into conservative variables
    real(DP), dimension(:), intent(out), optional :: R_ij

    ! OPTIONAL: transformation matrix into characteristic variables
    real(DP), dimension(:), intent(out), optional :: L_ij
!</output>
!</subroutine>

  end subroutine euler_calcCharacteristics3d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcBoundaryvalues3d(DbdrNormal, DpointNormal,&
      DbdrValue, ibdrCondType, Du, Du0, istatus)

!<description>
    ! This subroutine computes the boundary values for a given node in 3D
!</description>

!<input>
    ! normal vector at the boundary
    real(DP), dimension(:), intent(in) :: DbdrNormal

    ! normal vector at the point on the boundary
    real(DP), dimension(:), intent(in) :: DpointNormal

    ! evaluated boundary values
    real(DP), dimension(:), intent(in) :: DbdrValue

    ! initial solution from the previous time step
    real(DP), dimension(:), intent(in) :: Du0

    ! type of boundary condition
    integer, intent(in) :: ibdrCondType
!</input>

!<inputoutput>
    ! computed boundary values
    real(DP), dimension(:), intent(inout) :: Du

    ! OPTIONAL: status of the callback function
    integer, intent(inout), optional :: istatus
!</inputoutput>
!</subroutine>

  end subroutine euler_calcBoundaryvalues3d

  !*****************************************************************************

!<subroutine>

  subroutine euler_hadaptCallbackScalar3d(rcollection, iOperation,&
      Ivertices, Ielements)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 3D. The solution vector is assumed
    ! to be store in scalar interleave format.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(in) :: iOperation

    ! Array of vertices involved in the adaptivity step
    integer, dimension(:), intent(in) :: Ivertices

    ! Array of elements involved in the adaptivity step
    integer, dimension(:), intent(in) :: Ielements
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer, save :: rsolution
    real(DP), dimension(:), pointer, save :: p_Dsolution
    integer :: ivar


    ! What operation should be performed?
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! Retrieve solution vector from colletion
      rsolution => rcollection%p_rvectorQuickAccess1

      ! Check if solution is stored in interleave format
      if (rsolution%nblocks .ne. 1) then
        call output_line('Vector is not in interleave format!',&
            OU_CLASS_WARNING,OU_MODE_STD,'euler_hadaptCallbackScalar3d')
        call sys_halt()
      end if

      ! Set pointer
      call lsysbl_getbase_double(rsolution, p_Dsolution)
      
      ! Call the general callback function
      call flagship_hadaptCallback3d(rcollection, iOperation,&
          Ivertices, Ielements)
      
      
    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback3d(rcollection, iOperation,&
          Ivertices, Ielements)
      

    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. NVAR3D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolution, NVAR3D*Ivertices(1),&
            .false., .true.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR3D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolution, NVAR3D*Ivertices(1),&
            .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      do ivar = 1, NVAR3D
        p_Dsolution((Ivertices(1)-1)*NVAR3D+ivar) = &
            0.5_DP*(p_Dsolution((Ivertices(2)-1)*NVAR3D+ivar)+&
                    p_Dsolution((Ivertices(3)-1)*NVAR3D+ivar))
      end do

      ! Call the general callback function
      call flagship_hadaptCallback3d(rcollection, iOperation,&
          Ivertices, Ielements)


    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR3D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolution, NVAR3D*Ivertices(1),&
            .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      do ivar = 1, NVAR3D
        p_Dsolution((Ivertices(1)-1)*NVAR3D+ivar) = &
            0.25_DP*(p_Dsolution((Ivertices(2)-1)*NVAR3D+ivar)+&
                     p_Dsolution((Ivertices(3)-1)*NVAR3D+ivar)+&
                     p_Dsolution((Ivertices(4)-1)*NVAR3D+ivar)+&
                     p_Dsolution((Ivertices(5)-1)*NVAR3D+ivar))
      end do

      ! Call the general callback function
      call flagship_hadaptCallback3d(rcollection, iOperation,&
          Ivertices, Ielements)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (Ivertices(2) .ne. 0) then
        do ivar = 1, NVAR3D
          p_Dsolution((Ivertices(1)-1)*NVAR3D+ivar) = &
              p_Dsolution((Ivertices(2)-1)*NVAR3D+ivar)
        end do
      else
        do ivar = 1, NVAR3D
          p_Dsolution((Ivertices(1)-1)*NVAR3D+ivar) = 0.0_DP
        end do
      end if
      
      ! Call the general callback function
      call flagship_hadaptCallback3d(rcollection, iOperation,&
          Ivertices, Ielements)


    case DEFAULT
      ! Call the general callback function
      call flagship_hadaptCallback3d(rcollection, iOperation,&
          Ivertices, Ielements)

    end select
    
  end subroutine euler_hadaptCallbackScalar3d

  !*****************************************************************************

!<subroutine>

  subroutine euler_hadaptCallbackBlock3d(rcollection, iOperation,&
      Ivertices, Ielements)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 3D. The solution vector is assumed
    ! to be store in block format.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(in) :: iOperation

    ! Array of vertices involved in the adaptivity step
    integer, dimension(:), intent(in) :: Ivertices

    ! Array of elements involved in the adaptivity step
    integer, dimension(:), intent(in) :: Ielements
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer, save :: rsolution
    real(DP), dimension(:), pointer, save :: p_Dsolution
    integer :: ivar,neq


    ! What operation should be performed?
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! Retrieve solution vector from colletion
      rsolution => rcollection%p_rvectorQuickAccess1

      ! Check if solution is stored in interleave format
      if (rsolution%nblocks .ne. NVAR3D) then
        call output_line('Vector is not in block format!',&
            OU_CLASS_WARNING,OU_MODE_STD,'euler_hadaptCallbackBlock3d')
        call sys_halt()
      end if

      ! Set pointer
      call lsysbl_getbase_double(rsolution, p_Dsolution)
      
      ! Call the general callback function
      call flagship_hadaptCallback3d(rcollection, iOperation,&
          Ivertices, Ielements)
      
      
    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback3d(rcollection, iOperation,&
          Ivertices, Ielements)
      

    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. NVAR3D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolution, NVAR3D*Ivertices(1),&
            .false., .true.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR3D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolution, NVAR3D*Ivertices(1),&
            .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      neq = rsolution%NEQ/NVAR3D
      do ivar = 1, NVAR3D
        p_Dsolution((ivar-1)*neq+Ivertices(1)) = &
            0.5_DP*(p_Dsolution((ivar-1)*neq+Ivertices(2))+&
                    p_Dsolution((ivar-1)*neq+Ivertices(3)) )
      end do

      ! Call the general callback function
      call flagship_hadaptCallback3d(rcollection, iOperation,&
          Ivertices, Ielements)

      
    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR3D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolution, NVAR3D*Ivertices(1),&
            .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      neq = rsolution%NEQ/NVAR3D
      do ivar = 1, NVAR3D
        p_Dsolution((ivar-1)*neq+Ivertices(1)) =&
            0.25_DP*(p_Dsolution((ivar-1)*neq+Ivertices(2))+&
                     p_Dsolution((ivar-1)*neq+Ivertices(3))+&
                     p_Dsolution((ivar-1)*neq+Ivertices(4))+&
                     p_Dsolution((ivar-1)*neq+Ivertices(5)) )
      end do

      ! Call the general callback function
      call flagship_hadaptCallback3d(rcollection, iOperation,&
          Ivertices, Ielements)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (Ivertices(2) .ne. 0) then
        neq = rsolution%NEQ/NVAR3D
        do ivar = 1, NVAR3D
          p_Dsolution((ivar-1)*neq+Ivertices(1)) = &
              p_Dsolution((ivar-1)*neq+Ivertices(2))
        end do
      else
        neq = rsolution%NEQ/NVAR3D
        do ivar = 1, NVAR3D
          p_Dsolution((ivar-1)*neq+Ivertices(1)) = 0.0_DP
        end do
      end if
      
      ! Call the general callback function
      call flagship_hadaptCallback3d(rcollection, iOperation,&
          Ivertices, Ielements)


    case DEFAULT
      ! Call the general callback function
      call flagship_hadaptCallback3d(rcollection, iOperation,&
          Ivertices, Ielements)

    end select
    
  end subroutine euler_hadaptCallbackBlock3d
  
end module euler_callback3d
