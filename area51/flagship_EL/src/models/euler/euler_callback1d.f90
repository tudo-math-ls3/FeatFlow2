!##############################################################################
!# ****************************************************************************
!# <name> euler_callback1d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the compressible Euler/Navier-Stokes equations in 1D.
!#
!# The following callback functions are available:
!#
!# 1.) euler_calcFluxGalerkin1d
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#
!# 2.) euler_calcFluxGalerkinNoBdr1d
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#        without assembling the symmetric boundary contribution
!#
!# 3.) euler_calcFluxScalarDiss1d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting scalar artificial viscosities
!#
!# 4.) euler_calcFluxTensorDiss1d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting tensorial artificial viscosities
!#
!# 5.) euler_calcFluxRusanov1d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting the Rusanov artificial diffusion
!#
!# 6.) euler_calcMatrixDiagonalDiag1d
!#     -> Computes local matrix for diagonal entry
!#
!# 7.) euler_calcMatrixDiagonal1d
!#     -> Computes local matrix for diagonal entry
!#
!# 8.) euler_calcMatrixGalerkinDiag1d
!#     -> Computes local matrices for standard Galerkin scheme
!#
!# 9.) euler_calcMatrixGalerkin1d
!#     -> Computes local matrices for standard Galerkin scheme
!#
!# 10.) euler_calcMatrixScalarDissDiag1d
!#      -> Computes local matrices for low-order discretization
!#         adopting scalar artificial viscosities
!#
!# 11.) euler_calcMatrixScalarDiss1d
!#      -> Computes local matrices for low-order discretization
!#         adopting scalar artificial viscosities
!#
!# 12.) euler_calcMatrixTensorDissDiag1d
!#      -> Computes local matrices for low-order discretization
!#         adopting tensorial artificial viscosities
!#
!# 13.) euler_calcMatrixTensorDiss1d
!#      -> Computes local matrices for low-order discretization
!#         adopting tensorial artificial viscosities
!#
!# 14.) euler_calcMatrixRusanovDiag1d
!#      -> Computes local matrices for low-order discretization
!#         adopting the Rusanov artificial viscosities
!#
!# 15.) euler_calcMatrixRusanov1d
!#      -> Computes local matrices for low-order discretization
!#         adopting the Rusanov flux artificial viscosities
!#
!# 16.) euler_calcCharacteristics1d
!#      -> Computes characteristic variables in 1D
!#
!# 17.) euler_calcBoundaryvalues1d
!#      -> Computes the boundary values for a given node
!#
!# 18.) euler_hadaptCallbackScalar1d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 1D, whereby the vector is stored in interleave format
!#
!# 19.) euler_hadaptCallbackBlock1d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 1D, whereby the vector is stored in block format
!#
!# </purpose>
!##############################################################################

module euler_callback1d

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
  public :: euler_calcFluxGalerkin1d
  public :: euler_calcFluxGalerkinNoBdr1d
  public :: euler_calcFluxScalarDiss1d
  public :: euler_calcFluxTensorDiss1d
  public :: euler_calcFluxRusanov1d
  public :: euler_calcMatrixDiagonalDiag1d
  public :: euler_calcMatrixDiagonal1d
  public :: euler_calcMatrixGalerkinDiag1d
  public :: euler_calcMatrixGalerkin1d
  public :: euler_calcMatrixScalarDissDiag1d
  public :: euler_calcMatrixScalarDiss1d
  public :: euler_calcMatrixTensorDissDiag1d
  public :: euler_calcMatrixTensorDiss1d
  public :: euler_calcMatrixRusanovDiag1d
  public :: euler_calcMatrixRusanov1d
  public :: euler_calcCharacteristics1d
  public :: euler_calcBoundaryvalues1d
  public :: euler_hadaptCallbackScalar1d
  public :: euler_hadaptCallbackBlock1d

contains
  
  !*****************************************************************************
  
!<subroutine>

  pure subroutine euler_calcFluxGalerkin1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the standard
    ! Galerkin discretization in 1D.
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
    
  end subroutine euler_calcFluxGalerkin1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxGalerkinNoBdr1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the TVD
    ! discretization in 1D. The symmetric boundary contributions
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

  end subroutine euler_calcFluxGalerkinNoBdr1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxScalarDiss1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 1D using scalar dissipation.
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

  end subroutine euler_calcFluxScalarDiss1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxTensorDiss1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 1D using tensorial dissipation.
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

  end subroutine euler_calcFluxTensorDiss1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxRusanov1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 1D using the Rusanov dissipation.
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

  end subroutine euler_calcFluxRusanov1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixDiagonalDiag1d(U_i, C_ii, i, dscale, K_ii)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! for the diagonal block of the global operator in 1D
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

  end subroutine euler_calcMatrixDiagonalDiag1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixDiagonal1d(U_i, C_ii, i, dscale, K_ii)

!<description>
    ! This subroutine computes the Galerkin matrices for the diagonal
    ! block of the global operator in 1D
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

  end subroutine euler_calcMatrixDiagonal1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixGalerkinDiag1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices in 1D
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

  end subroutine euler_calcMatrixGalerkinDiag1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixGalerkin1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the Galerkin matrices in 1D
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

  end subroutine euler_calcMatrixGalerkin1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixScalarDissDiag1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies scalar artificial viscosities in 1D
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

    
  end subroutine euler_calcMatrixScalarDissDiag1d

!*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixScalarDiss1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies scalar artificial viscosities in 1D
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

  end subroutine euler_calcMatrixScalarDiss1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixTensorDissDiag1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies tensorial artificial viscosities in 1D
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

  end subroutine euler_calcMatrixTensorDissDiag1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixTensorDiss1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies tensorial artificial viscosities in 1D
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

  end subroutine euler_calcMatrixTensorDiss1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixRusanovDiag1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 1D
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

  end subroutine euler_calcMatrixRusanovDiag1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixRusanov1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 1D
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

  end subroutine euler_calcMatrixRusanov1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcCharacteristics1d(&
      U_i, U_j, Dweight, W_ij, Lbd_ij, R_ij, L_ij)

!<description>
    ! This subroutine computes the characteristic variables in 1D
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

  end subroutine euler_calcCharacteristics1d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcBoundaryvalues1d(DbdrNormal, DpointNormal,&
      DbdrValue, ibdrCondType, Du, Du0, istatus)

!<description>
    ! This subroutine computes the boundary values for a given node in 1D
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

  end subroutine euler_calcBoundaryvalues1d

  !*****************************************************************************

!<subroutine>

  subroutine euler_hadaptCallbackScalar1d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 1D. The solution vector is assumed
    ! to be store in scalar interleave format.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(in) :: iOperation
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
            OU_CLASS_WARNING,OU_MODE_STD,'euler_hadaptCallbackScalar1d')
        call sys_halt()
      end if

      ! Set pointer
      call lsysbl_getbase_double(rsolution, p_Dsolution)
      
      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)
      
      
    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)
      

    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. NVAR1D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR1D*rcollection%IquickAccess(1), .false., .true.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR1D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR1D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      do ivar = 1, NVAR1D
        p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR1D+ivar) = &
            0.5_DP*(p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR1D+ivar)+&
                    p_Dsolution((rcollection%IquickAccess(3)-1)*NVAR1D+ivar))
      end do

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (rcollection%IquickAccess(2) .ne. 0) then
        do ivar = 1, NVAR1D
          p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR1D+ivar) = &
              p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR1D+ivar)
        end do
      else
        do ivar = 1, NVAR1D
          p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR1D+ivar) = 0.0_DP
        end do
      end if
      
      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case default
      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)

    end select
    
  end subroutine euler_hadaptCallbackScalar1d

  !*****************************************************************************

!<subroutine>

  subroutine euler_hadaptCallbackBlock1d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 1D. The solution vector is assumed
    ! to be store in block format.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(in) :: iOperation
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
      if (rsolution%nblocks .ne. NVAR1D) then
        call output_line('Vector is not in block format!',&
            OU_CLASS_WARNING,OU_MODE_STD,'euler_hadaptCallbackBlock1d')
        call sys_halt()
      end if

      ! Set pointer
      call lsysbl_getbase_double(rsolution, p_Dsolution)
      
      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)
      
      
    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)
      

    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. NVAR1D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR1D*rcollection%IquickAccess(1), .false., .true.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR1D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR1D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      neq = rsolution%NEQ/NVAR1D
      do ivar = 1, NVAR1D
        p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = &
            0.5_DP*(p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))+&
                    p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(3)) )
      end do

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (rcollection%IquickAccess(2) .ne. 0) then
        neq = rsolution%NEQ/NVAR1D
        do ivar = 1, NVAR1D
          p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = &
              p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))
        end do
      else
        neq = rsolution%NEQ/NVAR1D
        do ivar = 1, NVAR1D
          p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = 0.0_DP
        end do
      end if
      
      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case default
      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)

    end select
    
  end subroutine euler_hadaptCallbackBlock1d
 
end module euler_callback1d
