!##############################################################################
!# ****************************************************************************
!# <name> stdoperators </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains some routines to assemble basic matrices that are
!# frequently used in most Finite Element based programs.
!#
!# The following routines can be found in this module:
!#
!# 1.) stdop_assembleLaplaceMatrix
!#     -> Assembles a standard Laplace matrix
!#
!# 2.) stdop_assembleSimpleMatrix
!#     -> Assembles simple standard matrices like a mass matrix.
!#
!# </purpose>
!##############################################################################

module stdoperators

!$use omp_lib
  use fsystem
  use genoutput
  use linearsystemscalar
  use scalarpde
  use derivatives
  use spatialdiscretisation
  use bilinearformevaluation
  use perfconfig

  implicit none

  private

  public :: stdop_assembleLaplaceMatrix
  public :: stdop_assembleSimpleMatrix

contains

  ! ***************************************************************************

!<subroutine>

  subroutine stdop_assembleLaplaceMatrix (rmatrix,bclear,dalpha, &
      rcubatureInfo, rperfconfig)

!<description>
  ! This routine assembles a Laplace matrix into rmatrix.
!</description>

!<input>
  ! OPTIONAL: If set to TRUE (standard), the content of rmatrix is set to 0.0
  ! before assembling the matrix.
  logical, intent(in), optional :: bclear

  ! OPTIONAL: Constant coefficient in front of the matrix, which is multiplied
  ! to all entries. If not specified, 1.0 is assumed.
  real(DP), intent(in), optional :: dalpha

  ! OPTIONAL: A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), optional, target :: rcubatureInfo

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! Matrix structure where to save the entries of the Laplace matrix to.
  ! The structure of the matrix (KCOL, KLD,...) as well as the discretisation
  ! structure must already be given.
  ! If the array for the matrix content does not exist, a new array is created.
  ! If the array exist, the new entries of the Laplace operator overwrite
  ! the old entries (if bclear=true) or are added to the old entries
  ! (if bclear=false).
  type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP) :: dalpha1
    logical :: bclear1
    type(t_spatialDiscretisation), pointer :: p_rdiscr => null()
    type(t_scalarCubatureInfo), target :: rtempCubatureInfo
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform

    bclear1 = .true.
    dalpha1 = 1.0_DP
    if (present(bclear)) bclear1=bclear
    if (present(dalpha)) dalpha1=dalpha

    ! Now try to get a pointer to the spatial discretisation - either from the trial
    ! or the test space. We just need to know whether we are in 1D, 2D or 3D.
    p_rdiscr => rmatrix%p_rspatialDiscrTrial

    if(.not. associated(p_rdiscr)) then
      call output_line('Matrix does not have a discretisation!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'stdop_assembleLaplaceMatrix')
      call sys_halt()
    end if

    ! Now we can set up the bilinearform.
    select case(p_rdiscr%ndimension)
    case (1)
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_DERIV1D_X
      rform%Idescriptors(2,1) = DER_DERIV1D_X

    case (2)
      rform%itermCount = 2
      rform%Idescriptors(1,1) = DER_DERIV2D_X
      rform%Idescriptors(2,1) = DER_DERIV2D_X
      rform%Idescriptors(1,2) = DER_DERIV2D_Y
      rform%Idescriptors(2,2) = DER_DERIV2D_Y

    case (3)
      rform%itermCount = 3
      rform%Idescriptors(1,1) = DER_DERIV3D_X
      rform%Idescriptors(2,1) = DER_DERIV3D_X
      rform%Idescriptors(1,2) = DER_DERIV3D_Y
      rform%Idescriptors(2,2) = DER_DERIV3D_Y
      rform%Idescriptors(1,3) = DER_DERIV3D_Z
      rform%Idescriptors(2,3) = DER_DERIV3D_Z

    end select

    ! In the standard case, we have constant coefficients:
    rform%ballCoeffConstant = .true.
    rform%Dcoefficients(1:rform%itermCount)  = dalpha1

    ! If we do not have it, create a cubature info structure that
    ! defines how to do the assembly.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_createDefCubStructure(rmatrix%p_rspatialDiscrTrial,&
          rtempCubatureInfo,CUB_GEN_DEPR_BILFORM)
      p_rcubatureInfo => rtempCubatureInfo
    else
      p_rcubatureInfo => rcubatureInfo
    end if

    ! Now we can build the matrix entries.
    call bilf_buildMatrixScalar (rform,bclear1,rmatrix,&
        p_rcubatureInfo,rperfconfig=rperfconfig)

    ! Release the assembly structure if necessary.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_releaseCubStructure(rtempCubatureInfo)
    end if

  end subroutine stdop_assembleLaplaceMatrix

  ! ***************************************************************************

!<subroutine>

  subroutine stdop_assembleSimpleMatrix (rmatrix,cderivTrial,cderivTest,dalpha,&
      bclear,rcubatureInfo,rperfconfig)

!<description>
  ! This routine assembles a simple Finite element matrix into rmatrix.
  ! "Simple" means here, that a matrix entry consists of exactly one
  ! trial/test function pair and a constant coefficient:
  !   <tex> $$ a_ij = int_{\Omega} \alpha * \psi^j_k * \phi^i_l $$ </tex>
  ! With $k,l$ being derivative quantifiers.
  ! This allows to assemble mass matrices or similar matrices
  ! directly. (A mass matrix for example can be assembled using
  ! cderivTrial = cderivTest = DER_FUNC.)
!</description>

!<input>
  ! A derivative quantifier for the trial function. One of the DER_XXXX constants,
  ! e.g. DER_DERIV_X for X- or DER_DERIV_Y for Y-derivative.
  integer, intent(in) :: cderivTrial

  ! A derivative quantifier for the test function. One of the DER_XXXX constants,
  ! e.g. DER_DERIV_X for X- or DER_DERIV_Y for Y-derivative.
  integer, intent(in) :: cderivTest

  ! Constant coefficient in front of the matrix, which is multiplied
  ! to all entries. If not specified, a value of 1.0 is assumed.
  real(DP), intent(in), optional :: dalpha

  ! OPTIONAL: If set to TRUE (standard), the content of rmatrix is set to 0.0
  ! before assembling the matrix.
  logical, intent(in), optional :: bclear

  ! OPTIONAL: A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), optional, target :: rcubatureInfo

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! Matrix structure where to save the entries of the Laplace matrix to.
  ! The structure of the matrix (KCOL, KLD,...) as well as the discretisation
  ! structure must already be given.
  ! If the array for the matrix content does not exist, a new array is created.
  ! If the array exist, the new entries of the Laplace operator overwrite
  ! the old entries (if bclear=true) or are added to the old entries
  ! (if bclear=false).
  type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP) :: dalpha1
    logical :: bclear1
    type(t_scalarCubatureInfo), target :: rtempCubatureInfo
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform

    bclear1 = .true.
    dalpha1 = 1.0_DP
    if (present(bclear)) bclear1=bclear
    if (present(dalpha)) dalpha1=dalpha

    ! Set up bilinear form for generating the mass matrix
    rform%itermCount = 1
    rform%Idescriptors(1,1) = cderivTrial
    rform%Idescriptors(2,1) = cderivTest

    ! We have constant coefficients:
    rform%ballCoeffConstant = .true.
    rform%Dcoefficients(1)  = dalpha1

    ! If we do not have it, create a cubature info structure that
    ! defines how to do the assembly.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_createDefCubStructure(rmatrix%p_rspatialDiscrTrial,&
          rtempCubatureInfo,CUB_GEN_DEPR_BILFORM)
      p_rcubatureInfo => rtempCubatureInfo
    else
      p_rcubatureInfo => rcubatureInfo
    end if

    ! Now we can build the matrix entries.
    call bilf_buildMatrixScalar (rform,bclear1,rmatrix,&
        p_rcubatureInfo,rperfconfig=rperfconfig)

    ! Release the assembly structure if necessary.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_releaseCubStructure(rtempCubatureInfo)
    end if

  end subroutine stdop_assembleSimpleMatrix

end module stdoperators
