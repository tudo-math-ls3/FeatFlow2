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

  use fsystem
  use genoutput
  use linearsystemscalar
  use bilinearformevaluation
  
  implicit none
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine stdop_assembleLaplaceMatrix (rmatrix,bclear,dalpha)
  
!<description>
  ! This routine assembles a Laplace matrix into rmatrix.
!</description>

!<input>
  ! OPTIONAL: If set to TRUE (standard), the content of rmatrix is set to 0.0
  ! before assembling the matrix.
  logical, intent(IN), optional :: bclear
  
  ! OPTIONAL: Constant coefficient in front of the matrix, which is multiplied
  ! to all entries. If not specified, 1.0 is assumed.
  real(DP), intent(IN), optional :: dalpha
!</input>

!<inputoutput>
  ! Matrix structure where to save the entries of the Laplace matrix to.
  ! The structure of the matrix (KCOL, KLD,...) as well as the discretisation
  ! structure must already be given.
  ! If the array for the matrix content does not exist, a new array is created.
  ! If the array exist, the new entries of the Laplace operator overwrite
  ! the old entries (if bclear=true) or are added to the old entries 
  ! (if bclear=false).
  type(t_matrixScalar), intent(INOUT) :: rmatrix
!</inputoutput>
  
!</subroutine>

    ! local variables
    real(DP) :: dalpha1
    logical :: bclear1
    type(t_spatialDiscretisation), pointer :: p_rdiscr => NULL()
    
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
    rform%BconstantCoeff = .true.
    rform%Dcoefficients(1)  = dalpha1

    ! Now we can build the matrix entries.
    call bilf_buildMatrixScalar (rform,bclear1,rmatrix)

  end subroutine stdop_assembleLaplaceMatrix

  ! ***************************************************************************

!<subroutine>

  subroutine stdop_assembleSimpleMatrix (rmatrix,cderivTrial,cderivTest,dalpha,&
      bclear)
  
!<description>
  ! This routine assembles a simple Finite element matrix into rmatrix.
  ! "Simple" means here, that a matrix entry consists of exactly one
  ! trial/test function pair and a constant coefficient:
  ! $$ a_ij = int_{\Omega} \alpha * \psi^j_k * \phi^i_l $$
  ! With $k,l$ being derivative quantifiers.
  ! This allowes to assemble mass matrices or similar matrices
  ! directly. (A mass matrix for example can be assembled using
  ! cderivTrial = cderivTest = DER_FUNC.)
!</description>

!<input>
  ! A derivative quantifier for the trial function. One of the DER_XXXX constants,
  ! e.g. DER_DERIV_X for X- or DER_DERIV_Y for Y-derivative.
  integer, intent(IN) :: cderivTrial
  
  ! A derivative quantifier for the test function. One of the DER_XXXX constants,
  ! e.g. DER_DERIV_X for X- or DER_DERIV_Y for Y-derivative.
  integer, intent(IN) :: cderivTest
  
  ! Constant coefficient in front of the matrix, which is multiplied
  ! to all entries. If not specified, a value of 1.0 is assumed.
  real(DP), intent(IN), optional :: dalpha
  
  ! OPTIONAL: If set to TRUE (standard), the content of rmatrix is set to 0.0
  ! before assembling the matrix.
  logical, intent(IN), optional :: bclear
!</input>

!<inputoutput>
  ! Matrix structure where to save the entries of the Laplace matrix to.
  ! The structure of the matrix (KCOL, KLD,...) as well as the discretisation
  ! structure must already be given.
  ! If the array for the matrix content does not exist, a new array is created.
  ! If the array exist, the new entries of the Laplace operator overwrite
  ! the old entries (if bclear=true) or are added to the old entries 
  ! (if bclear=false).
  type(t_matrixScalar), intent(INOUT) :: rmatrix
!</inputoutput>
  
!</subroutine>

    ! local variables
    real(DP) :: dalpha1
    logical :: bclear1
    
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
    rform%BconstantCoeff = .true.
    rform%Dcoefficients(1)  = dalpha1

    ! Now we can build the matrix entries.
    call bilf_buildMatrixScalar (rform,bclear1,rmatrix)

  end subroutine stdop_assembleSimpleMatrix

end module stdoperators
