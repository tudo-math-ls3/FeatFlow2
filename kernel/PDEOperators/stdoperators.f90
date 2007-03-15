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

MODULE stdoperators

  USE fsystem
  USE linearsystemscalar
  USE bilinearformevaluation
  
  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE stdop_assembleLaplaceMatrix (rmatrix,bclear,dalpha)
  
!<description>
  ! This routine assembles a Laplace matrix into rmatrix.
!</description>

!<input>
  ! OPTIONAL: If set to TRUE (standard), the content of rmatrix is set to 0.0
  ! before assembling the matrix.
  LOGICAL, INTENT(IN), OPTIONAL :: bclear
  
  ! OPTIONAL: Constant coefficient in front of the matrix, which is multiplied
  ! to all entries. If not specified, 1.0 is assumed.
  REAL(DP), INTENT(IN), OPTIONAL :: dalpha
!</input>

!<inputoutput>
  ! Matrix structure where to save the entries of the Laplace matrix to.
  ! The structure of the matrix (KCOL, KLD,...) as well as the discretisation
  ! structure must already be given.
  ! If the array for the matrix content does not exist, a new array is created.
  ! If the array exist, the new entries of the Laplace operator overwrite
  ! the old entries (if bclear=true) or are added to the old entries 
  ! (if bclear=false).
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
!</inputoutput>
  
!</subroutine>

    ! local variables
    REAL(DP) :: dalpha1
    LOGICAL :: bclear1
    
    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_bilinearForm) :: rform
    
    bclear1 = .TRUE.
    dalpha1 = 1.0_DP
    IF (PRESENT(bclear)) bclear1=bclear
    IF (PRESENT(dalpha)) dalpha1=dalpha

    ! For assembling of the entries, we need a bilinear form, 
    ! which first has to be set up manually.
    ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
    ! scalar system matrix in 2D.
    
    rform%itermCount = 2
    rform%Idescriptors(1,1) = DER_DERIV_X
    rform%Idescriptors(2,1) = DER_DERIV_X
    rform%Idescriptors(1,2) = DER_DERIV_Y
    rform%Idescriptors(2,2) = DER_DERIV_Y

    ! In the standard case, we have constant coefficients:
    rform%ballCoeffConstant = .TRUE.
    rform%BconstantCoeff = .TRUE.
    rform%Dcoefficients(1)  = dalpha1
    rform%Dcoefficients(2)  = dalpha1

    ! Now we can build the matrix entries.
    CALL bilf_buildMatrixScalar (rform,bclear1,rmatrix)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE stdop_assembleSimpleMatrix (rmatrix,cderivTrial,cderivTest,dalpha,&
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
  INTEGER, INTENT(IN) :: cderivTrial
  
  ! A derivative quantifier for the test function. One of the DER_XXXX constants,
  ! e.g. DER_DERIV_X for X- or DER_DERIV_Y for Y-derivative.
  INTEGER, INTENT(IN) :: cderivTest
  
  ! Constant coefficient in front of the matrix, which is multiplied
  ! to all entries. If not specified, a value of 1.0 is assumed.
  REAL(DP), INTENT(IN), OPTIONAL :: dalpha
  
  ! OPTIONAL: If set to TRUE (standard), the content of rmatrix is set to 0.0
  ! before assembling the matrix.
  LOGICAL, INTENT(IN), OPTIONAL :: bclear
!</input>

!<inputoutput>
  ! Matrix structure where to save the entries of the Laplace matrix to.
  ! The structure of the matrix (KCOL, KLD,...) as well as the discretisation
  ! structure must already be given.
  ! If the array for the matrix content does not exist, a new array is created.
  ! If the array exist, the new entries of the Laplace operator overwrite
  ! the old entries (if bclear=true) or are added to the old entries 
  ! (if bclear=false).
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
!</inputoutput>
  
!</subroutine>

    ! local variables
    REAL(DP) :: dalpha1
    LOGICAL :: bclear1
    
    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_bilinearForm) :: rform
    
    bclear1 = .TRUE.
    dalpha1 = 1.0_DP
    IF (PRESENT(bclear)) bclear1=bclear
    IF (PRESENT(dalpha)) dalpha1=dalpha

    ! Set up bilinear form for generating the mass matrix
    rform%itermCount = 1
    rform%Idescriptors(1,1) = cderivTrial
    rform%Idescriptors(2,1) = cderivTest

    ! We have constant coefficients:
    rform%ballCoeffConstant = .TRUE.
    rform%BconstantCoeff = .TRUE.
    rform%Dcoefficients(1)  = dalpha1

    ! Now we can build the matrix entries.
    CALL bilf_buildMatrixScalar (rform,bclear1,rmatrix)

  END SUBROUTINE

END MODULE
