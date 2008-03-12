!##############################################################################
!# ****************************************************************************
!# <name> groupfemsystem </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic routines for applying the group-finite
!# element formulation to systems of conservation laws.
!# The technique was proposed by C.A.J. Fletcher in:
!#
!#     C.A.J. Fletcher, "The group finite element formulation"
!#     Computer Methods in Applied Mechanics and Engineering (ISSN 0045-7825),
!#     vol. 37, April 1983, p. 225-244.
!#
!# The group finite element formulation uses the same basis functions for the
!# unknown solution and the fluxes. This allows for an efficient matrix assemble,
!# whereby the constant coefficient matrices can be assembled once and for all
!# at the beginning of the simulation and each time the grid is modified.
!#
!# Moreover, this module allows to modifying discrete operators by means of
!# the algebraic flux correction (AFC) methodology proposed by Kuzmin, Moeller
!# and Turek in a series of publications. As a starting point for systems of
!# conservation laws, the reader is referred to the book chapter
!#
!#     D. Kuzmin and M. Moeller, "Algebraic flux correction II. Compsressible Euler
!#     Equations", In: D. Kuzmin et al. (eds), Flux-Corrected Transport: Principles, 
!#     Algorithms, and Applications, Springer, 2005, 207-250.
!#
!#
!# A more detailed description of the algorithms is given in the comments of the
!# subroutine implementing the corresponding discretisation schemes. All methods
!# are based on the stabilisation structure t_afcstab which is defined in the 
!# underlying module "afcstabilisation". The initialisation as a system stabilisation
!# structure is done by the routine gfsys_initStabilisation. 
!#
!# There are three types of routines. The gfsys_buildDivOperator routines can be
!# used to assemble the discrete divergence operators resulting from the standard
!# Galerkin finite element discretisation plus some discretely defined artificial
!# viscosities. This technique represents a generalisation of the discrete upwinding
!# approach which has been used to construct upwind finite element scheme for
!# scalar conservation laws (see module groupfemscalar for details).
!#
!# The second type of routines is given by gfsc_buildResidualXXX. They can be
!# used to update/initialise the residual term applying some sort of algebraic
!# flux correction. Importantly, the family of AFC schemes gives rise to nonlinear
!# algebraic equations that need to be solved iteratively. Thus, it is usefull to
!# build the compensating antidiffusion into the residual vector rather than the
!# right hand side. However, it is still possible to give a negative scaling factor.
!#
!# The third type of routines is used to assemble the Jacobian matrix for Newton's
!# method. Here, the exact Jacobian matrix is approximated by means of second-order
!# divided differences whereby the "perturbation parameter" is specified by the
!# user. You should be aware of the fact, that in general the employed  flux limiters
!# are not differentiable globally so that the construction of the Jacobian matrix
!# is somehow "delicate". Even though the routines will produce some matrix without
!# warnings, this matrix may be singular and/or ill-conditioned.
!#
!# The following routines are available:
!#
!# 1.) gfsys_initStabilisation = gfsys_initStabilisationScalar /
!#                               gfsys_initStabilisationBlock
!#     -> initialize the stabilisation structure
!#
!# 2.) gfsys_isMatrixCompatible = gfsys_isMatrixCompatibleScalar /
!#                                gfsys_isMatrixCompatibleBlock
!#     -> checks wether a matrix and a stabilisation structure are compatible
!#
!# 3.) gfsys_isVectorCompatible = gfsys_isVectorCompatibleScalar /
!#                                gfsys_isVectorCompatibleBlock
!#     -> checks wether a vector and a stabilisation structure are compatible
!#
!# 4.) gfsys_buildDivOperator = gfsys_buildDivOperatorScalar /
!#                              gfsys_buildDivOperatorBlock
!#     -> assemble the global operator that results from the discretisation
!#        of the divergence term $div(F)$ by means of the Galerkin method 
!#        and some sort of artificial dissipation (if required)
!#
!# 5.) gfsys_buildResidual = gfsys_buildResScalar /
!#                           gfsys_buildResBlock
!#     -> assemble the residual vector
!#
!# 6.) gfsys_buildResidualFCT = gfsys_buildResScalarFCT /
!#                              gfsys_buildResBlockFCT
!#     -> assemble the residual for FEM-FCT stabilisation
!#
!# 7.) gfsys_buildResidualTVD = gfsys_buildResScalarTVD /
!#                              gfsys_buildResBlockTVD
!#     -> assemble the residual for FEM-TVD stabilisation
!#
!# 8.) gfsys_buildDivJacobian = gfsys_buildDivJacobianScalar /
!#                              gfsys_buildDivJacobianBlock
!#     -> assemble the Jacobian matrix
!#
!# The following internal routines are available:
!#
!# 1.) gfsys_getbase_double
!#     -> assign the pointers of an array to a given block matrix
!#
!# 2.) gfsys_getbase_single
!#     -> assign the pointers of an array to a given block matrix
!#
!# 3.) gfsys_getbase_int
!#     -> assign the pointers of an array to a given block matrix
!#
!# </purpose>
!##############################################################################

MODULE groupfemsystem

  USE afcstabilisation
  USE fsystem
  USE genoutput
  USE linearsystemblock
  USE linearsystemscalar
  USE storage

  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: gfsys_initStabilisation
  PUBLIC :: gfsys_isMatrixCompatible
  PUBLIC :: gfsys_isVectorCompatible
  PUBLIC :: gfsys_buildDivOperator
  PUBLIC :: gfsys_buildResidual
  PUBLIC :: gfsys_buildResidualFCT
  PUBLIC :: gfsys_buildResidualTVD
  PUBLIC :: gfsys_buildDivJacobian

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<constants>
!<constantblock description="Global constants for directional splitting">

  ! unit vector in X-direction in 1D
  REAL(DP), DIMENSION(NDIM1D) :: XDir1D = (/ 1._DP /)

  ! unit vector in X-direction in 2D
  REAL(DP), DIMENSION(NDIM2D) :: XDir2D = (/ 1._DP, 0._DP /)

  ! unit vector in Y-direction in 2D
  REAL(DP), DIMENSION(NDIM2D) :: YDir2D = (/ 0._DP, 1._DP /)

  ! unit vector in X-direction in 3D
  REAL(DP), DIMENSION(NDIM3D) :: XDir3D = (/ 1._DP, 0._DP, 0._DP /)

  ! unit vector in Y-direction in 3D
  REAL(DP), DIMENSION(NDIM3D) :: YDir3D = (/ 0._DP, 1._DP, 0._DP /)

  ! unit vector in Z-direction in 3D
  REAL(DP), DIMENSION(NDIM3D) :: ZDir3D = (/ 0._DP, 0._DP, 1._DP /)

!</constantblock>
!</constants>

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<types>
!<typeblock>

  ! An auxiliary array that is used to store the content of all
  ! scalar submatrices simultaneously. Note that this derived type
  ! is PRIVATE and cannot be accessed from outside of this module.
  PRIVATE :: t_array
  TYPE t_array

    ! Pointer to the double-valued matrix data
    REAL(DP), DIMENSION(:), POINTER :: Da

    ! Pointer to the single-valued matrix data
    REAL(SP), DIMENSION(:), POINTER :: Fa

    ! Pointer to the integer-valued matrix data
    INTEGER, DIMENSION(:), POINTER  :: Ia
  END TYPE t_array
!</typeblock>
!</types>

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

  INTERFACE gfsys_initStabilisation
    MODULE PROCEDURE gfsys_initStabilisationScalar
    MODULE PROCEDURE gfsys_initStabilisationBlock
  END INTERFACE

  INTERFACE gfsys_isMatrixCompatible
    MODULE PROCEDURE gfsys_isMatrixCompatibleScalar
    MODULE PROCEDURE gfsys_isMatrixCompatibleBlock
  END INTERFACE

  INTERFACE gfsys_isVectorCompatible
    MODULE PROCEDURE gfsys_isVectorCompatibleScalar
    MODULE PROCEDURE gfsys_isVectorCompatibleBlock
  END INTERFACE

  INTERFACE gfsys_buildDivOperator
     MODULE PROCEDURE gfsys_buildDivOperatorScalar
     MODULE PROCEDURE gfsys_buildDivOperatorBlock
  END INTERFACE

  INTERFACE gfsys_buildResidual
    MODULE PROCEDURE gfsys_buildResScalar
    MODULE PROCEDURE gfsys_buildResBlock
  END INTERFACE

  INTERFACE gfsys_buildResidualFCT
    MODULE PROCEDURE gfsys_buildResScalarFCT
    MODULE PROCEDURE gfsys_buildResBlockFCT
  END INTERFACE

  INTERFACE gfsys_buildResidualTVD
    MODULE PROCEDURE gfsys_buildResScalarTVD
    MODULE PROCEDURE gfsys_buildResBlockTVD
  END INTERFACE

  INTERFACE gfsys_buildDivJacobian
    MODULE PROCEDURE gfsys_buildDivJacobianScalar
    MODULE PROCEDURE gfsys_buildDivJacobianBlock
  END INTERFACE

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

CONTAINS

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsys_initStabilisationBlock(rmatrixBlockTemplate, rafcstab)

!<description>
    ! This subroutine initialises the discrete stabilisation structure
    ! for use as a scalar stabilisation. The template matrix is used
    ! to determine the number of equations and the number of edges.
    !
    ! Note that the matrix is required as block matrix. If this matrix
    ! contains only one block, then the scalar counterpart of this
    ! subroutine is called with the corresponding scalar submatrix.
!</description>

!<input>
    ! template block matrix
    TYPE(t_matrixBlock), INTENT(IN) :: rmatrixBlockTemplate
!</input>

!<inputoutput>
    ! discrete operator structure
    TYPE(t_afcstab), INTENT(INOUT)  :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER :: i

    ! Check if block matrix has only one block
    IF (rmatrixBlockTemplate%ndiagblocks .EQ. 1) THEN
      CALL gfsys_initStabilisationScalar(&
          rmatrixBlockTemplate%RmatrixBlock(1,1), rafcstab)
      RETURN
    END IF

    ! Check if matrix exhibits group structure
    IF (rmatrixBlockTemplate%imatrixSpec .NE.&
        LSYSBS_MSPEC_GROUPMATRIX) THEN
      CALL output_line('Block matrix must have group structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_initStabilisationBlock')
      CALL sys_halt()
    END IF
    
    ! Set atomic data from first block
    rafcstab%NVAR  = rmatrixBlockTemplate%ndiagblocks
    rafcstab%NEQ   = rmatrixBlockTemplate%RmatrixBlock(1,1)%NEQ
    rafcstab%NEDGE = INT(0.5*(rmatrixBlockTemplate%RmatrixBlock(1,1)%NA-&
                     rmatrixBlockTemplate%RmatrixBlock(1,1)%NEQ), I32)

    ! What kind of stabilisation are we?
    SELECT CASE(rafcstab%ctypeAFCstabilisation)
      
    CASE (AFCSTAB_GALERKIN, AFCSTAB_UPWIND)
      ! do nothing

    CASE (AFCSTAB_FEMTVD)

      ! Nodal vectors
      ALLOCATE(rafcstab%RnodalVectors(6))
      DO i = 1, 6
        CALL lsyssc_createVector(rafcstab%RnodalVectors(i),&
            rafcstab%NEQ, rafcstab%NVAR, .FALSE., ST_DOUBLE)
      END DO
      
    CASE (AFCSTAB_FEMFCT)

      ! Nodal vectors
      ALLOCATE(rafcstab%RnodalVectors(7))
      DO i = 1, 7
        CALL lsyssc_createVector(rafcstab%RnodalVectors(i),&
            rafcstab%NEQ, rafcstab%NVAR, .FALSE., ST_DOUBLE)
      END DO

      ! Edgewise vectors
      ALLOCATE(rafcstab%RedgeVectors(4))
      DO i = 1, 4
        CALL lsyssc_createVector(rafcstab%RedgeVectors(i),&
            rafcstab%NEDGE, rafcstab%NVAR, .FALSE., ST_DOUBLE)
      END DO

    CASE DEFAULT
      CALL output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_initStabilisationBlock')
      CALL sys_halt()
    END SELECT

    ! Set specifier
    rafcstab%iSpec = AFCSTAB_INITIALISED
  END SUBROUTINE gfsys_initStabilisationBlock

  ! *****************************************************************************

  SUBROUTINE gfsys_initStabilisationScalar(rmatrixTemplate, rafcstab)

!<description>
    ! This subroutine initialises the discrete stabilisation structure
    ! for use as a scalar stabilisation. The template matrix is used
    ! to determine the number of equations and the number of edges.
    !
    ! Note that the matrix is required as scalar matrix. It can be
    ! stored in interleave format.
!</description>

!<input>
    ! template matrix
    TYPE(t_matrixScalar), INTENT(IN) :: rmatrixTemplate
!</input>

!<inputoutput>
    ! discrete operator structure
    TYPE(t_afcstab), INTENT(INOUT)   :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER :: i
  
    ! Set atomic data
    rafcstab%NVAR  = rmatrixTemplate%NVAR
    rafcstab%NEQ   = rmatrixTemplate%NEQ
    rafcstab%NEDGE = INT(0.5*(rmatrixTemplate%NA-rmatrixTemplate%NEQ), I32)

    ! What kind of stabilisation are we?
    SELECT CASE(rafcstab%ctypeAFCstabilisation)
      
    CASE (AFCSTAB_GALERKIN, AFCSTAB_UPWIND)
      ! do nothing
      
    CASE (AFCSTAB_FEMTVD)

      ! Nodal vectors
      ALLOCATE(rafcstab%RnodalVectors(6))
      DO i = 1, 6
        CALL lsyssc_createVector(rafcstab%RnodalVectors(i),&
            rafcstab%NEQ, rafcstab%NVAR, .FALSE., ST_DOUBLE)
      END DO

    CASE (AFCSTAB_FEMFCT)

      ! Nodal vectors
      ALLOCATE(rafcstab%RnodalVectors(7))
      DO i = 1, 7
        CALL lsyssc_createVector(rafcstab%RnodalVectors(i),&
            rafcstab%NEQ, rafcstab%NVAR, .FALSE., ST_DOUBLE)
      END DO

      ! Edgewise vectors
      ALLOCATE(rafcstab%RedgeVectors(4))
      DO i = 1, 4
        CALL lsyssc_createVector(rafcstab%RedgeVectors(i),&
            rafcstab%NEDGE, rafcstab%NVAR, .FALSE., ST_DOUBLE)
      END DO

    CASE DEFAULT
      CALL output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_initStabilisationScalar')
      CALL sys_halt()
    END SELECT
    
    ! Set specifier
    rafcstab%iSpec = AFCSTAB_INITIALISED
  END SUBROUTINE gfsys_initStabilisationScalar

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE gfsys_isMatrixCompatibleBlock(rafcstab, rmatrixBlock, bcompatible)

!<description>
    ! This subroutine checks whether a block matrix and a discrete 
    ! stabilisation structure are compatible to each other, 
    ! i.e. if they share the same structure, size and so on.
    !
    ! If the matrix has only one block, then the scalar counterpart of this
    ! subroutine is called with the corresponding scalar submatrix.
    ! Otherwise, the matrix is required to possess group structure.
!</description>

!<input>
    ! block matrix
    TYPE(t_matrixBlock), INTENT(IN) :: rmatrixBlock

    ! stabilisation structure
    TYPE(t_afcstab), INTENT(IN)     :: rafcstab
!</input>

!<output>
    ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
    ! whether matrix and stabilisation are compatible or not.
    ! If not given, an error will inform the user if the matrix/operator are
    ! not compatible and the program will halt.
    LOGICAL, INTENT(OUT), OPTIONAL  :: bcompatible
!</output>
!</subroutine>

    ! Check if matrix has only one block
    IF (rmatrixBlock%ndiagblocks .EQ. 1) THEN
      CALL gfsys_isMatrixCompatibleScalar(rafcstab,&
          rmatrixBlock%RmatrixBlock(1,1), bcompatible)
      RETURN
    END IF

    ! Check if matrix exhibits group structure
    IF (rmatrixBlock%imatrixSpec .NE. LSYSBS_MSPEC_GROUPMATRIX) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        CALL output_line('Block matrix must have group structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_isMatrixCompatibleBlock')
        CALL sys_halt()
      END IF
    END IF

    ! Matrix/operator must have the same size
    IF (rafcstab%NVAR  .NE. rmatrixBlock%ndiagblocks .OR.&
        rafcstab%NEQ   .NE. rmatrixBlock%RmatrixBlock(1,1)%NEQ  .OR.&
        rafcstab%NEDGE .NE. INT(0.5*(rmatrixBlock%RmatrixBlock(1,1)%NA-&
        rmatrixBlock%RmatrixBlock(1,1)%NEQ),I32)) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        CALL output_line('Matrix/Operator not compatible, different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfssy_isMatrixCompatibleBlock')
        CALL sys_halt()
      END IF
    END IF
  END SUBROUTINE gfsys_isMatrixCompatibleBlock

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE gfsys_isMatrixCompatibleScalar(rafcstab, rmatrix, bcompatible)

!<description>
    ! This subroutine checks whether a scalar matrix and a discrete 
    ! stabilisation structure are compatible to each other, 
    ! i.e. if they share the same structure, size and so on.
    !
    ! Note that the matrix is required as scalar matrix. It can be
    ! stored in interleave format.
!</description>

!<input>
    ! scalar matrix
    TYPE(t_matrixScalar), INTENT(IN) :: rmatrix

    ! stabilisation structure
    TYPE(t_afcstab), INTENT(IN)      :: rafcstab
!</input>

!<output>
    ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
    ! whether matrix and stabilisation are compatible or not.
    ! If not given, an error will inform the user if the matrix/operator are
    ! not compatible and the program will halt.
    LOGICAL, INTENT(OUT), OPTIONAL   :: bcompatible
!</output>
!</subroutine>

    ! Matrix/operator must have the same size
    IF (rafcstab%NEQ   .NE. rmatrix%NEQ  .OR.&
        rafcstab%NVAR  .NE. rmatrix%NVAR .OR.&
        rafcstab%NEDGE .NE. INT(0.5*(rmatrix%NA-rmatrix%NEQ),I32)) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        CALL output_line('Matrix/Operator not compatible, different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfssy_isMatrixCompatibleScalar')
        CALL sys_halt()
      END IF
    END IF
  END SUBROUTINE gfsys_isMatrixCompatibleScalar

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsys_isVectorCompatibleBlock(rafcstab, rvectorBlock, bcompatible)

!<description>
    ! This subroutine checks whether a block vector and a stabilisation
    ! structure are compatible to each other, i.e., share the same
    ! structure, size and so on.
    !
    ! If the vectors has only one block, then the scalar counterpart of
    ! this subroutine is called with the corresponding scalar subvector.
!</description>

!<input>
    ! block vector
    TYPE(t_vectorBlock), INTENT(IN) :: rvectorBlock

    ! stabilisation structure
    TYPE(t_afcstab), INTENT(IN)     :: rafcstab
!</input>

!<output>
    ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
    ! whether matrix and stabilisation are compatible or not.
    ! If not given, an error will inform the user if the matrix/operator are
    ! not compatible and the program will halt.
    LOGICAL, INTENT(OUT), OPTIONAL  :: bcompatible
!</output>
!</subroutine>

    ! Check if block vectors has just one block
    IF (rvectorBlock%nblocks .EQ. 1) THEN
      CALL gfsys_isVectorCompatibleScalar(rafcstab,&
          rvectorBlock%RvectorBlock(1), bcompatible)
      RETURN
    END IF

    ! Vector/operator must have the same size
    IF (rafcstab%NEQ*rafcstab%NVAR .NE. rvectorBlock%NEQ) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        CALL output_line('Vector/Operator not compatible, different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_isVectorCompatibleBlock')
        CALL sys_halt()
      END IF
    END IF
  END SUBROUTINE gfsys_isVectorCompatibleBlock

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsys_isVectorCompatibleScalar(rafcstab, rvector, bcompatible)

!<description>
    ! This subroutine checks whether a vector and a stabilisation
    ! structure are compatible to each other, i.e., share the same
    ! structure, size and so on.
    !
    ! Note that the vector is required as scalar vector. It can be
    ! stored in interleave format.
!</description>

!<input>
    ! scalar vector
    TYPE(t_vectorScalar), INTENT(IN) :: rvector

    ! stabilisation structure
    TYPE(t_afcstab), INTENT(IN)      :: rafcstab
!</input>

!<output>
    ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
    ! whether matrix and stabilisation are compatible or not.
    ! If not given, an error will inform the user if the matrix/operator are
    ! not compatible and the program will halt.
    LOGICAL, INTENT(OUT), OPTIONAL   :: bcompatible
!</output>
!</subroutine>

    ! Matrix/operator must have the same size
    IF (rafcstab%NEQ   .NE. rvector%NEQ .OR.&
        rafcstab%NVAR  .NE. rvector%NVAR) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        CALL output_line('Vector/Operator not compatible, different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_isVectorCompatibleScalar')
        CALL sys_halt()
      END IF
    END IF
  END SUBROUTINE gfsys_isVectorCompatibleScalar

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsys_buildDivOperatorBlock(RmatrixC, ru,&
      fcb_calcMatrix, dscale, bclear, rmatrixL)

!<description>
    ! This subroutine assembles the discrete transport operator $K$ so that
    ! $$(Ku)_i=-\sum_{j\ne i}\bf c}_{ij}\cdot({\bf F}_j-{\bf F}_i),$$
    ! whereby the flux difference can be represented as
    ! $${\bf F}_j-{\bf F}_i=\hat{\bf A}_{ij}(u_j-u_i)$$
    ! and the matrix $\hat{\bf A}$ corresponds to the Jacobian tensor
    ! evaluated for some special set of averaged variables.
    !
    ! Note that this routine is designed for block matrices/vectors. 
    ! If there is only one block, then the corresponding scalar routine 
    ! is called. Otherwise, the global operator is treated as block matrix.
    ! This block matrix has to be in group structure, that is, the structure
    ! of subblock(1,1) will serve as template for all other submatrices.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: RmatrixC

    ! solution vector
    TYPE(t_vectorBlock), INTENT(IN)                :: ru

    ! scaling factor
    REAL(DP), INTENT(IN)                           :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    LOGICAL, INTENT(IN)                            :: bclear

    ! callback functions to compute local Roe matrix and
    ! local dissipation matrix
    INCLUDE 'intf_gfsyscallback.inc'
!</input>   

!<inputoutput>
    ! global transport operator
    TYPE(t_matrixBlock), INTENT(INOUT)             :: rmatrixL
!</inputoutput>
!</subroutine>
    
    ! local variables
    TYPE(t_array), DIMENSION(ru%nblocks,ru%nblocks)  :: rarray
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER      :: p_Kld,p_Ksep
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER      :: p_Kdiagonal
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER      :: p_Kcol
    REAL(DP), DIMENSION(:), POINTER                  :: p_Cx,p_Cy,p_Cz,p_u
    INTEGER :: h_Ksep
    INTEGER :: idim,ndim,iblock,jblock
    LOGICAL :: bisFullMatrix

    ! Check if block vector contains only one block and if
    ! global operator is stored in interleave format.
    IF (ru%nblocks .EQ. 1   .AND.&
        rmatrixL%ndiagblocks .EQ. 1) THEN
      CALL gfsys_buildDivOperatorScalar(RmatrixC, ru%RvectorBlock(1),&
          fcb_calcMatrix, dscale, bclear, rmatrixL%RmatrixBlock(1,1))
      RETURN       
    END IF

    ! Check if matrix/vector is compatible
    CALL lsysbl_isMatrixCompatible(ru, rmatrixL)

    ! Check if all coefficient matrices are comptible
    ndim = SIZE(RmatrixC,1)
    DO idim = 1, ndim
      CALL lsyssc_isMatrixCompatible(RmatrixC(idim), rmatrixL%RmatrixBlock(1,1))
    END DO

    ! Check if block matrix exhibits group structure
    IF (rmatrixL%imatrixSpec .NE. LSYSBS_MSPEC_GROUPMATRIX) THEN
      CALL output_line('Block matrix must have group structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
      CALL sys_halt()
    END IF
    
    ! Clear matrix?
    IF (bclear) CALL lsysbl_clearMatrix(rmatrixL)
    
    ! Set pointers
    CALL lsyssc_getbase_Kld   (RmatrixC(1),p_Kld)
    CALL lsyssc_getbase_Kcol  (RmatrixC(1),p_Kcol)
    CALL lsysbl_getbase_double(ru, p_u)
    CALL gfsys_getbase_double (rmatrixL, rarray, bisFullMatrix)

    ! How many dimensions do we have?
    SELECT CASE(ndim)
    CASE (NDIM1D)
      CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)

    CASE (NDIM2D)
      CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
      CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)

    CASE (NDIM3D)
      CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
      CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
      CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)

    CASE DEFAULT
      CALL output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
      CALL sys_halt()
    END SELECT
    
    ! What kind of matrix are we?
    SELECT CASE(rmatrixL%RmatrixBlock(1,1)%cmatrixFormat)
    CASE(LSYSSC_MATRIX7)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      CALL storage_copy(RmatrixC(1)%h_Kld, h_Ksep)
      CALL storage_getbase_int(h_Ksep, p_Ksep, RmatrixC(1)%NEQ+1)
      
      ! What type of matrix are we?
      IF (bisFullMatrix) THEN
        
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doOperatorMat7_1D(p_Kld, p_Kcol, p_Ksep,&
              rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, rarray)

        CASE (NDIM2D)
          CALL doOperatorMat7_2D(p_Kld, p_Kcol, p_Ksep,&
              rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, rarray)

        CASE (NDIM3D)
          CALL doOperatorMat7_3D(p_Kld, p_Kcol, p_Ksep,&
              rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, rarray)

        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
          CALL sys_halt()
        END SELECT
        
      ELSE

        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doOperatorMat7Diag_1D(p_Kld, p_Kcol, p_Ksep,&
              rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, rarray)

        CASE (NDIM2D)
          CALL doOperatorMat7Diag_2D(p_Kld, p_Kcol, p_Ksep,&
              rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, rarray)

        CASE (NDIM3D)
          CALL doOperatorMat7Diag_3D(p_Kld, p_Kcol, p_Ksep,&
              rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, rarray)

        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
          CALL sys_halt()
        END SELECT

      END IF
      
      ! Release diagonal separator
      CALL storage_free(h_Ksep)
      
      
    CASE(LSYSSC_MATRIX9)

      ! Set pointers
      CALL lsyssc_getbase_Kdiagonal(rmatrixC(1), p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      CALL storage_copy(RmatrixC(1)%h_Kld, h_Ksep)
      CALL storage_getbase_int(h_Ksep, p_Ksep, RmatrixC(1)%NEQ+1)

      ! What type of matrix are we?
      IF (bisFullMatrix) THEN
        
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doOperatorMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, rarray)

        CASE (NDIM2D)
          CALL doOperatorMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, rarray)

        CASE (NDIM3D)
          CALL doOperatorMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, rarray)

        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
          CALL sys_halt()
        END SELECT

      ELSE

        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doOperatorMat9Diag_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, rarray)

        CASE (NDIM2D)
          CALL doOperatorMat9Diag_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, rarray)

        CASE (NDIM3D)
          CALL doOperatorMat9Diag_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, rarray)

        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
          CALL sys_halt()
        END SELECT

      END IF
      
      ! Release diagonal separator
      CALL storage_free(h_Ksep)
      
    CASE DEFAULT
      CALL output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
      CALL sys_halt()
    END SELECT
    
    ! Do we have to scale the matrix?
    IF (dscale .NE. 1.0_DP) THEN
      IF (bisFullMatrix) THEN
        DO iblock = 1, rmatrixL%ndiagblocks
          DO jblock = 1, rmatrixL%ndiagblocks
            CALL lsyssc_scaleMatrix(rmatrixL%Rmatrixblock(iblock,jblock), dscale)
          END DO
        END DO
      ELSE
        DO iblock = 1, rmatrixL%ndiagblocks
          CALL lsyssc_scaleMatrix(rmatrixL%Rmatrixblock(iblock,iblock), dscale)
        END DO
      END IF
    END IF
        

  CONTAINS

    ! Here, the working routines follow
    
    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 1D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doOperatorMat7Diag_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,S_ij,u_i,u_j
      REAL(DP), DIMENSION(NDIM1D) :: C_ij,C_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j
      INTEGER                     :: ivar
      
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute matrices
          CALL fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii)+A_ij(ivar)+S_ij(ivar)
            K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij)-A_ij(ivar)-S_ij(ivar)
            K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji)+A_ij(ivar)-S_ij(ivar)
            K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj)-A_ij(ivar)+S_ij(ivar)
          END DO
        END DO
      END DO
    END SUBROUTINE doOperatorMat7Diag_1D

    
    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 2D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doOperatorMat7Diag_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, Cy, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,S_ij,u_i,u_j
      REAL(DP), DIMENSION(NDIM2D) :: C_ij,C_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j
      INTEGER                     :: ivar
      
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute matrices
          CALL fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii)+A_ij(ivar)+S_ij(ivar)
            K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij)-A_ij(ivar)-S_ij(ivar)
            K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji)+A_ij(ivar)-S_ij(ivar)
            K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj)-A_ij(ivar)+S_ij(ivar)
          END DO
        END DO
      END DO
    END SUBROUTINE doOperatorMat7Diag_2D


    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 3D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doOperatorMat7Diag_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, Cy, Cz, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,S_ij,u_i,u_j
      REAL(DP), DIMENSION(NDIM3D) :: C_ij,C_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j
      INTEGER                     :: ivar
      
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute matrices
          CALL fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii)+A_ij(ivar)+S_ij(ivar)
            K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij)-A_ij(ivar)-S_ij(ivar)
            K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji)+A_ij(ivar)-S_ij(ivar)
            K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj)-A_ij(ivar)+S_ij(ivar)
          END DO
        END DO
      END DO
    END SUBROUTINE doOperatorMat7Diag_3D


    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 1D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doOperatorMat9Diag_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,S_ij,u_i,u_j
      REAL(DP), DIMENSION(NDIM1D) :: C_ij,C_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j
      INTEGER                     :: ivar

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute matrices
          CALL fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)

          ! Assemble the global operator
          DO ivar = 1, NVAR
            K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii)+A_ij(ivar)+S_ij(ivar)
            K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij)-A_ij(ivar)-S_ij(ivar)
            K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji)+A_ij(ivar)-S_ij(ivar)
            K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj)-A_ij(ivar)+S_ij(ivar)
          END DO
        END DO
      END DO
    END SUBROUTINE doOperatorMat9Diag_1D


    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 2D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doOperatorMat9Diag_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,S_ij,u_i,u_j
      REAL(DP), DIMENSION(NDIM2D) :: C_ij,C_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j
      INTEGER                     :: ivar

      ! Loop over all rows
      DO i = 1, NEQ

        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute matrices
          CALL fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii)+A_ij(ivar)+S_ij(ivar)
            K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij)-A_ij(ivar)-S_ij(ivar)
            K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji)+A_ij(ivar)-S_ij(ivar)
            K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj)-A_ij(ivar)+S_ij(ivar)
          END DO
        END DO
      END DO
    END SUBROUTINE doOperatorMat9Diag_2D


    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 2D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doOperatorMat9Diag_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, Cz, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,S_ij,u_i,u_j
      REAL(DP), DIMENSION(NDIM3D) :: C_ij,C_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j
      INTEGER                     :: ivar

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute matrices
          CALL fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii)+A_ij(ivar)+S_ij(ivar)
            K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij)-A_ij(ivar)-S_ij(ivar)
            K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji)+A_ij(ivar)-S_ij(ivar)
            K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj)-A_ij(ivar)+S_ij(ivar)
          END DO
        END DO
      END DO
    END SUBROUTINE doOperatorMat9Diag_3D


    !**************************************************************
    ! Assemble divergence operator K in 1D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doOperatorMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,S_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)           :: i,j
      INTEGER                        :: ivar,jvar,idx
      
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute matrices
          CALL fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii)+A_ij(idx)+S_ij(idx)
              K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij)-A_ij(idx)-S_ij(idx)
              K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji)+A_ij(idx)-S_ij(idx)
              K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj)-A_ij(idx)+S_ij(idx)
            END DO
          END DO
        END DO
      END DO
    END SUBROUTINE doOperatorMat7_1D

    
    !**************************************************************
    ! Assemble divergence operator K in 2D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doOperatorMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, Cy, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,S_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)           :: i,j
      INTEGER                        :: ivar,jvar,idx
      
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute matrices
          CALL fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii)+A_ij(idx)+S_ij(idx)
              K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij)-A_ij(idx)-S_ij(idx)
              K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji)+A_ij(idx)-S_ij(idx)
              K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj)-A_ij(idx)+S_ij(idx)
            END DO
          END DO
        END DO
      END DO
    END SUBROUTINE doOperatorMat7_2D


    !**************************************************************
    ! Assemble divergence operator K in 3D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doOperatorMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, Cy, Cz, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,S_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)           :: i,j
      INTEGER                        :: ivar,jvar,idx
      
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute matrices
          CALL fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii)+A_ij(idx)+S_ij(idx)
              K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij)-A_ij(idx)-S_ij(idx)
              K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji)+A_ij(idx)-S_ij(idx)
              K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj)-A_ij(idx)+S_ij(idx)
            END DO
          END DO
        END DO
      END DO
    END SUBROUTINE doOperatorMat7_3D

    
    !**************************************************************
    ! Assemble divergence operator K in 1D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doOperatorMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,S_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)           :: i,j
      INTEGER                        :: ivar,jvar,idx

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute matrices
          CALL fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii)+A_ij(idx)+S_ij(idx)
              K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij)-A_ij(idx)-S_ij(idx)
              K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji)+A_ij(idx)-S_ij(idx)
              K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj)-A_ij(idx)+S_ij(idx)
            END DO
          END DO
        END DO
      END DO
    END SUBROUTINE doOperatorMat9_1D


    !**************************************************************
    ! Assemble divergence operator K in 2D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doOperatorMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,S_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)           :: i,j
      INTEGER                        :: ivar,jvar,idx

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute matrices
          CALL fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii)+A_ij(idx)+S_ij(idx)
              K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij)-A_ij(idx)-S_ij(idx)
              K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji)+A_ij(idx)-S_ij(idx)
              K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj)-A_ij(idx)+S_ij(idx)
            END DO
          END DO
        END DO
      END DO
    END SUBROUTINE doOperatorMat9_2D


    !**************************************************************
    ! Assemble divergence operator K in 3D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doOperatorMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, Cz, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,S_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)           :: i,j
      INTEGER                        :: ivar,jvar,idx

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute matrices
          CALL fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii)+A_ij(idx)+S_ij(idx)
              K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij)-A_ij(idx)-S_ij(idx)
              K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji)+A_ij(idx)-S_ij(idx)
              K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj)-A_ij(idx)+S_ij(idx)
            END DO
          END DO
        END DO
      END DO
    END SUBROUTINE doOperatorMat9_3D
  END SUBROUTINE gfsys_buildDivOperatorBlock

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE gfsys_buildDivOperatorScalar(RmatrixC, ru,&
      fcb_calcMatrix, dscale, bclear, rmatrixL)

!<description>
    ! This subroutine assembles the discrete transport operator $K$ so that
    ! $$(Ku)_i=-\sum_{j\ne i}\bf c}_{ij}\cdot({\bf F}_j-{\bf F}_i),$$
    ! whereby the flux difference can be represented as
    ! $${\bf F}_j-{\bf F}_i=\hat{\bf A}_{ij}(u_j-u_i)$$
    ! and the matrix $\hat{\bf A}$ corresponds to the Jacobian tensor
    ! evaluated for some special set of averaged variables.
    !
    ! Note that this routine requires scalar matrices stored in the interleave
    ! format. 
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: RmatrixC
    
    ! scalar solution vector
    TYPE(t_vectorScalar), INTENT(IN)               :: ru

    ! scaling factor
    REAL(DP), INTENT(IN)                           :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    LOGICAL, INTENT(IN)                            :: bclear

    ! callback functions to compute local Roe matrix and
    ! local dissipation matrix
    INCLUDE 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! scalar transport operator
    TYPE(t_matrixScalar), INTENT(INOUT)            :: rmatrixL
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld,p_Ksep
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kdiagonal
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
    REAL(DP), DIMENSION(:), POINTER             :: p_Cx,p_Cy,p_Cz,p_L,p_u
    INTEGER :: h_Ksep
    INTEGER :: idim,ndim

    ! Check if all matrices/vectors are compatible
    CALL lsyssc_isMatrixCompatible(ru, rmatrixL, .FALSE.)

    ndim = SIZE(RmatrixC,1)
    DO idim = 1, ndim
      CALL lsyssc_isMatrixCompatible(rmatrixC(idim), rmatrixL)
    END DO
        
    ! Clear matrix?
    IF (bclear) CALL lsyssc_clearMatrix(rmatrixL)
    
    ! Set pointers
    CALL lsyssc_getbase_Kld   (RmatrixC(1), p_Kld)
    CALL lsyssc_getbase_Kcol  (RmatrixC(1), p_Kcol)
    CALL lsyssc_getbase_double(rmatrixL, p_L)
    CALL lsyssc_getbase_double(ru, p_u)

    ! How many dimensions do we have?
    SELECT CASE(ndim)
    CASE (NDIM1D)
      CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)

    CASE (NDIM2D)
      CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)
      CALL lsyssc_getbase_double(rmatrixC(2), p_Cy)

    CASE (NDIM3D)
      CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)
      CALL lsyssc_getbase_double(rmatrixC(2), p_Cy)
      CALL lsyssc_getbase_double(rmatrixC(3), p_Cz)

    CASE DEFAULT
      CALL output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
      CALL sys_halt()
    END SELECT
    
    ! What kind of matrix are we?
    SELECT CASE(rmatrixL%cmatrixFormat)
    CASE(LSYSSC_MATRIX7INTL)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      CALL storage_copy(RmatrixC(1)%h_Kld, h_Ksep)
      CALL storage_getbase_int(h_Ksep, p_Ksep, RmatrixC(1)%NEQ+1)
      
      ! What type of matrix are we?
      SELECT CASE(rmatrixL%cinterleavematrixFormat)
        
      CASE (LSYSSC_MATRIX1)
        
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doOperatorMat7_1D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, p_L)

        CASE (NDIM2D)
          CALL doOperatorMat7_2D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, p_L)

        CASE (NDIM3D)
          CALL doOperatorMat7_3D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, p_L)
          
        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
          CALL sys_halt()
        END SELECT
        
      CASE (LSYSSC_MATRIXD)
        
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doOperatorMat7Diag_1D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, p_L)

        CASE (NDIM2D)
          CALL doOperatorMat7Diag_2D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, p_L)

        CASE (NDIM3D)
          CALL doOperatorMat7Diag_3D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, p_L)

        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
          CALL sys_halt()
        END SELECT
        
      CASE DEFAULT
        CALL output_line('Unsupported interleave matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
        CALL sys_halt()
      END SELECT
      
      ! Release diagonal separator
      CALL storage_free(h_Ksep)
      
      
    CASE (LSYSSC_MATRIX9INTL)
      
      ! Set pointers
      CALL lsyssc_getbase_Kdiagonal(rmatrixL, p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      CALL storage_copy(RmatrixC(1)%h_Kld, h_Ksep)
      CALL storage_getbase_int(h_Ksep, p_Ksep, RmatrixC(1)%NEQ+1)
      
      ! What type of matrix are we?
      SELECT CASE(rmatrixL%cinterleavematrixFormat)
        
      CASE (LSYSSC_MATRIX1)
        
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doOperatorMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, p_L)

        CASE (NDIM2D)
          CALL doOperatorMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, p_L)

        CASE (NDIM3D)
          CALL doOperatorMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, p_L)

        CASE DEFAULT
          CALL output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
          CALL sys_halt() 
        END SELECT
        
      CASE (LSYSSC_MATRIXD)
        
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doOperatorMat9Diag_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, p_L)

        CASE (NDIM2D)
          CALL doOperatorMat9Diag_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, p_L)

        CASE (NDIM3D)
          CALL doOperatorMat9Diag_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, p_L)
          
        CASE DEFAULT
          CALL output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
          CALL sys_halt()
        END SELECT
        
      CASE DEFAULT
        CALL output_line('Unsupported interleave matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
        CALL sys_halt()
      END SELECT
      
      ! Release diagonal separator
      CALL storage_free(h_Ksep)
      
      
    CASE DEFAULT
      CALL output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
      CALL sys_halt()
    END SELECT
    

    ! Do we have to scale the matrix?
    IF (dscale .NE. 1._DP) THEN
      CALL lsyssc_scaleMatrix(rmatrixL, dscale)
    END IF

  CONTAINS
    
    ! Here, the working routines follow
    
    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 1D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doOperatorMat7Diag_1D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,S_ij
      REAL(DP), DIMENSION(NDIM1D) :: C_ij,C_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Compute matrices
          CALL fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+S_ij
          K(:,ij) = K(:,ij)-A_ij-S_ij
          K(:,ji) = K(:,ji)+A_ij-S_ij
          K(:,jj) = K(:,jj)-A_ij+S_ij
        END DO
      END DO
    END SUBROUTINE doOperatorMat7Diag_1D


    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 2D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doOperatorMat7Diag_2D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, Cy, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,S_ij
      REAL(DP), DIMENSION(NDIM2D) :: C_ij,C_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute matrices
          CALL fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+S_ij
          K(:,ij) = K(:,ij)-A_ij-S_ij
          K(:,ji) = K(:,ji)+A_ij-S_ij
          K(:,jj) = K(:,jj)-A_ij+S_ij
        END DO
      END DO
    END SUBROUTINE doOperatorMat7Diag_2D


    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 3D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doOperatorMat7Diag_3D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, Cy, Cz, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,S_ij
      REAL(DP), DIMENSION(NDIM3D) :: C_ij,C_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)
          
          ! Compute matrices
          CALL fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+S_ij
          K(:,ij) = K(:,ij)-A_ij-S_ij
          K(:,ji) = K(:,ji)+A_ij-S_ij
          K(:,jj) = K(:,jj)-A_ij+S_ij
        END DO
      END DO
    END SUBROUTINE doOperatorMat7Diag_3D

    
    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 1D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doOperatorMat9Diag_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, Cx, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,S_ij
      REAL(DP), DIMENSION(NDIM1D) :: C_ij,C_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute matrices
          CALL fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+S_ij
          K(:,ij) = K(:,ij)-A_ij-S_ij
          K(:,ji) = K(:,ji)+A_ij-S_ij
          K(:,jj) = K(:,jj)-A_ij+S_ij
        END DO
      END DO
    END SUBROUTINE doOperatorMat9Diag_1D


    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 2D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doOperatorMat9Diag_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, Cx, Cy, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,S_ij
      REAL(DP), DIMENSION(NDIM2D) :: C_ij,C_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute matrices
          CALL fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+S_ij
          K(:,ij) = K(:,ij)-A_ij-S_ij
          K(:,ji) = K(:,ji)+A_ij-S_ij
          K(:,jj) = K(:,jj)-A_ij+S_ij
        END DO
      END DO
    END SUBROUTINE doOperatorMat9Diag_2D


    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 3D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doOperatorMat9Diag_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, Cx, Cy, Cz, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,S_ij
      REAL(DP), DIMENSION(NDIM3D) :: C_ij,C_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)
          
          ! Compute matrices
          CALL fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+S_ij
          K(:,ij) = K(:,ij)-A_ij-S_ij
          K(:,ji) = K(:,ji)+A_ij-S_ij
          K(:,jj) = K(:,jj)-A_ij+S_ij
        END DO
      END DO
    END SUBROUTINE doOperatorMat9Diag_3D

    
    !**************************************************************
    ! Assemble divergence operator K in 1D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doOperatorMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,S_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)           :: i,j
      
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute matrices
          CALL fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+S_ij
          K(:,ij) = K(:,ij)-A_ij-S_ij
          K(:,ji) = K(:,ji)+A_ij-S_ij
          K(:,jj) = K(:,jj)-A_ij+S_ij
        END DO
      END DO
    END SUBROUTINE doOperatorMat7_1D

    
    !**************************************************************
    ! Assemble divergence operator K in 2D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doOperatorMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, Cy, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,S_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)           :: i,j
      
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute matrices
          CALL fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+S_ij
          K(:,ij) = K(:,ij)-A_ij-S_ij
          K(:,ji) = K(:,ji)+A_ij-S_ij
          K(:,jj) = K(:,jj)-A_ij+S_ij
        END DO
      END DO
    END SUBROUTINE doOperatorMat7_2D


    !**************************************************************
    ! Assemble divergence operator K in 3D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doOperatorMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, Cy, Cz, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,S_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)           :: i,j
      
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)
          
          ! Compute matrices
          CALL fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+S_ij
          K(:,ij) = K(:,ij)-A_ij-S_ij
          K(:,ji) = K(:,ji)+A_ij-S_ij
          K(:,jj) = K(:,jj)-A_ij+S_ij
        END DO
      END DO
    END SUBROUTINE doOperatorMat7_3D


    !**************************************************************
    ! Assemble divergence operator K in 1D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doOperatorMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, Cx, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,S_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)           :: i,j
      
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute matrices
          CALL fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+S_ij
          K(:,ij) = K(:,ij)-A_ij-S_ij
          K(:,ji) = K(:,ji)+A_ij-S_ij
          K(:,jj) = K(:,jj)-A_ij+S_ij
        END DO
      END DO
    END SUBROUTINE doOperatorMat9_1D


    !**************************************************************
    ! Assemble divergence operator K in 2D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doOperatorMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, Cx, Cy, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,S_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)           :: i,j
      
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute matrices
          CALL fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)

          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+S_ij
          K(:,ij) = K(:,ij)-A_ij-S_ij
          K(:,ji) = K(:,ji)+A_ij-S_ij
          K(:,jj) = K(:,jj)-A_ij+S_ij
        END DO
      END DO
    END SUBROUTINE doOperatorMat9_2D


    !**************************************************************
    ! Assemble divergence operator K in 3D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doOperatorMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, Cx, Cy, Cz, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,S_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)           :: i,j
      
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)
          
          ! Compute matrices
          CALL fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)

          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+S_ij
          K(:,ij) = K(:,ij)-A_ij-S_ij
          K(:,ji) = K(:,ji)+A_ij-S_ij
          K(:,jj) = K(:,jj)-A_ij+S_ij
        END DO
      END DO
    END SUBROUTINE doOperatorMat9_3D
  END SUBROUTINE gfsys_buildDivOperatorScalar

  ! *****************************************************************************

!<subroutine>
  
  SUBROUTINE gfsys_buildResBlock(RmatrixC, ru,&
      fcb_calcFlux, dscale, bclear, rres)

!<description>
    ! This subroutine assembles the residual vector for block vectors.
    ! If the vector contains only one block, then the scalar counterpart
    ! of this routine is called with the scalar subvector.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: RmatrixC

    ! solution vector
    TYPE(t_vectorBlock), INTENT(IN)                :: ru

    ! scaling factor
    REAL(DP), INTENT(IN)                           :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    LOGICAL, INTENT(IN)                            :: bclear

    ! callback functions to compute local matrices
    INCLUDE 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! residual vector
    TYPE(t_vectorBlock), INTENT(INOUT)             :: rres
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kcol
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kld,p_Ksep,p_Kdiagonal
    REAL(DP), DIMENSION(:), POINTER             :: p_Cx,p_Cy,p_Cz,p_u,p_res
    INTEGER :: h_Ksep

    
    ! Check if block vectors contain only one block.
    IF ((ru%nblocks .EQ. 1) .AND. (rres%nblocks .EQ. 1) ) THEN
      CALL gfsys_buildResScalar(RmatrixC, ru%RvectorBlock(1),&
          fcb_calcFlux, dscale, bclear, rres%RvectorBlock(1))
      RETURN       
    END IF

    ! Check if vectors are compatible
    CALL lsysbl_isVectorCompatible(ru, rres)

    ! Clear vector?
    IF (bclear) CALL lsysbl_clearVector(rres)

    ! Set pointers
    CALL lsyssc_getbase_Kld   (RmatrixC(1), p_Kld)
    CALL lsyssc_getbase_Kcol  (RmatrixC(1), p_Kcol)
    CALL lsysbl_getbase_double(ru,   p_u)
    CALL lsysbl_getbase_double(rres, p_res)

    ! Create diagonal Ksep=Kld
    h_Ksep = ST_NOHANDLE
    CALL storage_copy(RmatrixC(1)%h_Kld, h_Ksep)
    CALL storage_getbase_int(h_Ksep, p_Ksep, RmatrixC(1)%NEQ+1)
    
    
    ! What kind of matrix are we?
    SELECT CASE(RmatrixC(1)%cmatrixFormat)
    CASE(LSYSSC_MATRIX7)
      
      ! How many dimensions do we have?
      SELECT CASE(SIZE(RmatrixC,1))
      CASE (NDIM1D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL doResidualMat7_1D(p_Kld, p_Kcol, p_Ksep,&
            RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, dscale, p_res)
        
      CASE (NDIM2D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
         CALL doResidualMat7_2D(p_Kld, p_Kcol, p_Ksep,&
            RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, dscale, p_res)
        
      CASE (NDIM3D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
        CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
        CALL doResidualMat7_3D(p_Kld, p_Kcol, p_Ksep,&
            RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, dscale, p_res)
        
      CASE DEFAULT
        CALL output_line('Unsupported spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlock')
        CALL sys_halt()
      END SELECT
    

    CASE (LSYSSC_MATRIX9)

      ! Set pointer
      CALL lsyssc_getbase_Kdiagonal(RmatrixC(1), p_Kdiagonal)

      ! How many dimensions do we have?
      SELECT CASE(SIZE(RmatrixC,1))
      CASE (NDIM1D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL doResidualMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
            RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, dscale, p_res)
        
      CASE (NDIM2D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
         CALL doResidualMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
            RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, dscale, p_res)
        
      CASE (NDIM3D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
        CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
        CALL doResidualMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
            RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, dscale, p_res)
        
      CASE DEFAULT
        CALL output_line('Unsupported spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlock')
        CALL sys_halt()
      END SELECT
      

    CASE DEFAULT
      CALL output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlock')
      CALL sys_halt()
    END SELECT
    
    ! Release Ksep
    CALL storage_free(h_Ksep)
    
  CONTAINS
    
    ! Here, the working routines follow
    
    !**************************************************************
    ! Assemble residual vector in 1D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doResidualMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, u, dscale, res)
      
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(INOUT)      :: res
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: u_i,u_j,F_ij,F_ji
      REAL(DP), DIMENSION(NDIM1D) :: C_ij,C_ji
      INTEGER(PREC_MATIDX)        :: ij,ji
      INTEGER(PREC_VECIDX)        :: i,j
   
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute the fluxes
          CALL fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
        END DO
      END DO
    END SUBROUTINE doResidualMat7_1D


    !**************************************************************
    ! Assemble residual vector in 2D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doResidualMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, Cy, u, dscale, res)
      
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(INOUT)      :: res
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j,F_ij,F_ji
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j
   
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute the fluxes
          CALL fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
        END DO
      END DO
    END SUBROUTINE doResidualMat7_2D

    
    !**************************************************************
    ! Assemble residual vector in 3D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doResidualMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, Cy, Cz, u, dscale, res)
      
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(INOUT)      :: res
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j,F_ij,F_ji
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j
   
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute the fluxes
          CALL fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)

          ! Assemble residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
        END DO
      END DO
    END SUBROUTINE doResidualMat7_3D

    
    !**************************************************************
    ! Assemble residual vector in 1D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doResidualMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, u, dscale, res)
      
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(INOUT)      :: res
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j,F_ij,F_ji
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j
   
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute the fluxes
          CALL fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
        END DO
      END DO
    END SUBROUTINE doResidualMat9_1D


    !**************************************************************
    ! Assemble residual vector in 2D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doResidualMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, u, dscale, res)
      
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(INOUT)      :: res
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j,F_ij,F_ji
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j
   
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute the fluxes
          CALL fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
        END DO
      END DO
    END SUBROUTINE doResidualMat9_2D


    !**************************************************************
    ! Assemble residual vector in 3D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doResidualMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, Cz, u, dscale, res)
      
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(INOUT)      :: res
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j,F_ij,F_ji
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j
   
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)


          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute the fluxes
          CALL fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)

          ! Assemble residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
        END DO
      END DO
    END SUBROUTINE doResidualMat9_3D
  END SUBROUTINE gfsys_buildResBlock
  
  ! *****************************************************************************

!<subroutine>
  
  SUBROUTINE gfsys_buildResScalar(RmatrixC, ru,&
      fcb_calcFlux, dscale, bclear, rres)

!<description>
    ! This subroutine assembles the residual vector. Note that the vectors are
    ! required as scalar vectors which are stored in the interleave format.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: RmatrixC

    ! solution vector
    TYPE(t_vectorScalar), INTENT(IN)               :: ru

    ! scaling factor
    REAL(DP), INTENT(IN)                           :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    LOGICAL, INTENT(IN)                            :: bclear

    ! callback functions to compute local matrices
    INCLUDE 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! residual vector
    TYPE(t_vectorScalar), INTENT(INOUT)            :: rres
!</inputoutput>
!</subroutine>
  
    ! local variables
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kcol
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kld,p_Ksep,p_Kdiagonal
    REAL(DP), DIMENSION(:), POINTER             :: p_Cx,p_Cy,p_Cz,p_u,p_res
    INTEGER :: h_Ksep

    ! Check if vectors are compatible
    CALL lsyssc_isVectorCompatible(ru, rres)
    
    ! Clear vector?
    IF (bclear) CALL lsyssc_clearVector(rres)

    ! Set pointers
    CALL lsyssc_getbase_Kld   (RmatrixC(1), p_Kld)
    CALL lsyssc_getbase_Kcol  (RmatrixC(1), p_Kcol)
    CALL lsyssc_getbase_double(ru,   p_u)
    CALL lsyssc_getbase_double(rres, p_res)

    ! Create diagonal Ksep=Kld
    h_Ksep = ST_NOHANDLE
    CALL storage_copy(RmatrixC(1)%h_Kld, h_Ksep)
    CALL storage_getbase_int(h_Ksep, p_Ksep, RmatrixC(1)%NEQ+1)

    
    ! What kind of matrix are we?
    SELECT CASE(RmatrixC(1)%cmatrixFormat)
    CASE(LSYSSC_MATRIX7)
      
      ! How many dimensions do we have?
      SELECT CASE(SIZE(RmatrixC,1))
      CASE (NDIM1D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL doResidualMat7_1D(p_Kld, p_Kcol, p_Ksep,&
            RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, dscale, p_res)
        
      CASE (NDIM2D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
         CALL doResidualMat7_2D(p_Kld, p_Kcol, p_Ksep,&
            RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, dscale, p_res)
        
      CASE (NDIM3D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
        CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
        CALL doResidualMat7_3D(p_Kld, p_Kcol, p_Ksep,&
            RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, dscale, p_res)
        
      CASE DEFAULT
        CALL output_line('Unsupported spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalar')
        CALL sys_halt()
      END SELECT
    

    CASE (LSYSSC_MATRIX9)

      ! Set pointer
      CALL lsyssc_getbase_Kdiagonal(RmatrixC(1), p_Kdiagonal)

      ! How many dimensions do we have?
      SELECT CASE(SIZE(RmatrixC,1))
      CASE (NDIM1D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL doResidualMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
            RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, dscale, p_res)
        
      CASE (NDIM2D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
         CALL doResidualMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
            RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, dscale, p_res)
        
      CASE (NDIM3D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
        CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
        CALL doResidualMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
            RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, dscale, p_res)
        
      CASE DEFAULT
        CALL output_line('Unsupported spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalar')
        CALL sys_halt()
      END SELECT
      

    CASE DEFAULT
      CALL output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalar')
      CALL sys_halt()
    END SELECT

    ! Release Ksep
    CALL storage_free(h_Ksep)

  CONTAINS

    ! Here, the working routines follow
    
    !**************************************************************
    ! Assemble residual vector in 1D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doResidualMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, u, dscale, res)
      
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: res
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)      :: F_ij,F_ji
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j
   
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Compute the fluxes
          CALL fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
        END DO
      END DO
    END SUBROUTINE doResidualMat7_1D


    !**************************************************************
    ! Assemble residual vector in 2D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doResidualMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, Cy, u, dscale, res)
      
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: res
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)      :: F_ij,F_ji
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j
   
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute the fluxes
          CALL fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
        END DO
      END DO
    END SUBROUTINE doResidualMat7_2D

    
    !**************************************************************
    ! Assemble residual vector in 3D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doResidualMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, Cy, Cz, u, dscale, res)
      
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: res
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)      :: F_ij,F_ji
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j
   
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)
          
          ! Compute the fluxes
          CALL fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
        END DO
      END DO
    END SUBROUTINE doResidualMat7_3D

    
    !**************************************************************
    ! Assemble residual vector in 1D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doResidualMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, u, dscale, res)
      
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: res
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)      :: F_ij,F_ji
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j
   
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute the fluxes
          CALL fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)

          ! Assemble residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
        END DO
      END DO
    END SUBROUTINE doResidualMat9_1D


    !**************************************************************
    ! Assemble residual vector in 2D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doResidualMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, u, dscale, res)
      
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,NEQ), INTENT(IN)           :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,NEQ), INTENT(INOUT)        :: res
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)      :: F_ij,F_ji
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j
   
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Compute the fluxes
          CALL fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)

          ! Assemble residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
        END DO
      END DO
    END SUBROUTINE doResidualMat9_2D


    !**************************************************************
    ! Assemble residual vector in 3D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doResidualMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, Cz, u, dscale, res)
      
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: res
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)      :: F_ij,F_ji
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j
   
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute averaged coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Compute the fluxes
          CALL fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
        END DO
      END DO
    END SUBROUTINE doResidualMat9_3D
  END SUBROUTINE gfsys_buildResScalar

  ! *****************************************************************************
  
!<subroutine>
  
  SUBROUTINE gfsys_buildResBlockFCT(RmatrixC, rmatrixML, ru,&
      fcb_calcFlux, fcb_calcRawFlux, fcb_calcCharacteristics, rafcstab,&
      dscale, theta, tstep, bclear, rres, ruPredict, rmatrixMC)

!<description>
    ! This subroutine assembles the residual vector for FEM-FCT schemes.
    ! If the vectors contain only one block, then the scalar counterpart
    ! of this routine is called with the scalar subvectors.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: RmatrixC

    ! lumped mass matrix
    TYPE(t_matrixScalar), INTENT(IN)               :: rmatrixML

    ! solution vector
    TYPE(t_vectorBlock), INTENT(IN)                :: ru

    ! predicted solution vector (low-order)
    TYPE(t_vectorBlock), INTENT(IN), OPTIONAL      :: ruPredict

    ! consistent mass matrix
    TYPE(t_matrixScalar), INTENT(IN), OPTIONAL     :: rmatrixMC

    ! scaling factor
    REAL(DP), INTENT(IN)                           :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    LOGICAL, INTENT(IN)                            :: bclear

    ! implicitness parameter
    REAL(DP), INTENT(IN)                           :: theta

    ! time step size
    REAL(DP), INTENT(IN)                           :: tstep

    ! callback functions to compute local matrices
    INCLUDE 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)                 :: rafcstab

    ! residual vector
    TYPE(t_vectorBlock), INTENT(INOUT)             :: rres
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain only one block.
    IF ((ru%nblocks .EQ. 1) .AND. (rres%nblocks .EQ. 1) ) THEN
      IF (PRESENT(ruPredict)) THEN
        CALL gfsys_buildResScalarFCT(RmatrixC, rmatrixML,&
            ru%RvectorBlock(1), fcb_calcFlux, fcb_calcRawFlux, fcb_calcCharacteristics,&
            rafcstab, dscale, theta, tstep, bclear, rres%RvectorBlock(1),&
            ruPredict%RvectorBlock(1), rmatrixMC)
      ELSE
        CALL gfsys_buildResScalarFCT(RmatrixC, rmatrixML,&
            ru%RvectorBlock(1), fcb_calcFlux, fcb_calcRawFlux, fcb_calcCharacteristics,&
            rafcstab, dscale, theta, tstep, bclear, rres%RvectorBlock(1),&
            rmatrixMC=rmatrixMC)
      END IF
      RETURN       
    END IF
    
    ! Check if vectors are compatible
    CALL lsysbl_isVectorCompatible(ru, rres)
    CALL gfsys_isVectorCompatible(rafcstab, ru)

    ! Clear vector?
    IF (bclear) CALL lsysbl_clearVector(rres)

  END SUBROUTINE gfsys_buildResBlockFCT

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE gfsys_buildResScalarFCT(RmatrixC, rmatrixML, ru,&
      fcb_calcFlux, fcb_calcRawFlux, fcb_calcCharacteristics, rafcstab,&
      dscale, theta, tstep, bclear, rres, ruPredict, rmatrixMC)
    
!<description>
    ! This subroutine assembles the residual vector for FEM-FCT schemes.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: RmatrixC

    ! lumped mass matrix
    TYPE(t_matrixScalar), INTENT(IN)               :: rmatrixML

    ! solution vector
    TYPE(t_vectorScalar), INTENT(IN)               :: ru

    ! predicted solution vector (low-order)
    TYPE(t_vectorScalar), INTENT(IN), OPTIONAL     :: ruPredict

    ! consistent mass matrix
    TYPE(t_matrixScalar), INTENT(IN), OPTIONAL     :: rmatrixMC

    ! scaling parameter
    REAL(DP), INTENT(IN)                           :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    LOGICAL, INTENT(IN)                            :: bclear

    ! implicitness parameter
    REAL(DP), INTENT(IN)                           :: theta

    ! time step size
    REAL(DP), INTENT(IN)                           :: tstep

    ! callback functions to compute local matrices
    INCLUDE 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)                 :: rafcstab
    
    ! residual vector
    TYPE(t_vectorScalar), INTENT(INOUT)            :: rres
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kcol
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kld,p_Ksep,p_Kdiagonal
    REAL(DP), DIMENSION(:), POINTER :: p_Cx,p_Cy,p_Cz,p_MC,p_ML,p_u,p_uPredict,p_res
    REAL(DP), DIMENSION(:), POINTER :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm,p_flux,p_flux0,p_fluxY,p_fluxY0
    INTEGER :: h_Ksep
    LOGICAL :: bmass

    ! Check if vectors are compatible
    CALL lsyssc_isVectorCompatible(ru, rres)
    CALL gfsys_isVectorCompatible (rafcstab, ru)

    ! Check if consistent mass matrix is available
    bmass = PRESENT(rmatrixMC)
    IF (bmass) THEN
      CALL lsyssc_getbase_double(rmatrixMC, p_MC)
    END IF

    ! Check if predictor is available, then store it in
    ! stabilisation structure for subsequent usage
    IF (PRESENT(ruPredict)) THEN
      CALL lsyssc_isVectorCompatible(ru, ruPredict)
      CALL lsyssc_copyVector (ruPredict, rafcstab%RnodalVectors(7))
    END IF

    ! Clear vector?
    IF (bclear) CALL lsyssc_clearVector(rres)

    ! Set pointers
    CALL lsyssc_getbase_Kld   (RmatrixC(1), p_Kld)
    CALL lsyssc_getbase_Kcol  (RmatrixC(1), p_Kcol)
    CALL lsyssc_getbase_double(rmatrixML,   p_ML)
    CALL lsyssc_getbase_double(ru,   p_u)
    CALL lsyssc_getbase_double(rres, p_res)


    ! Create diagonal Ksep=Kld
    h_Ksep = ST_NOHANDLE
    CALL storage_copy(RmatrixC(1)%h_Kld, h_Ksep)
    CALL storage_getbase_int(h_Ksep, p_Ksep, RmatrixC(1)%NEQ+1)

    ! What kind of stabilisation should be applied?
    SELECT CASE(rafcstab%ctypeAFCstabilisation)
      
    CASE (AFCSTAB_FEMFCT)

      ! Check if stabilisation is prepeared
      IF (IAND(rafcstab%iSpec, AFCSTAB_INITIALISED) .EQ. 0) THEN
        CALL output_line('Stabilisation has not been initialised',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarFCT')
        CALL sys_halt()
      END IF
      
      ! Set pointers
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(7), p_uPredict)
      CALL lsyssc_getbase_double(rafcstab%RedgeVectors(1) , p_flux)
      CALL lsyssc_getbase_double(rafcstab%RedgeVectors(2) , p_flux0)
      
      ! Set specifiers for Ps, Qs and Rs
      rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_LIMITER)
     
      ! What kind of matrix are we?
      SELECT CASE(RmatrixC(1)%cmatrixFormat)

      CASE(LSYSSC_MATRIX7)

        ! How many dimensions do we have?
        SELECT CASE(SIZE(RmatrixC,1))
        CASE (NDIM1D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)

        CASE (NDIM2D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)

          IF (PRESENT(ruPredict)) THEN
            
            CALL doInitFCTMat7_2D(p_Kld, p_Kcol, p_Ksep,&
                rafcstab%NEQ, rafcstab%NEDGE, rafcstab%NVAR,&
                p_Cx, p_Cy, p_ML, p_u, p_uPredict, dscale, theta, tstep,&
                p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_flux, p_flux0, p_res, p_MC)
            
          ELSE
            
            CALL doLimitFCTMat7_2D(p_Kld, p_Kcol, p_Ksep,&
                rafcstab%NEQ, rafcstab%NEDGE, rafcstab%NVAR,&
                p_Cx, p_Cy, p_u, p_uPredict, dscale, theta, tstep,&
                p_flux, p_flux0, p_res, p_MC)
            
          END IF

        CASE (NDIM3D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)

        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarFCT')
          CALL sys_halt()
        END SELECT
        
      CASE (LSYSSC_MATRIX9)
        
        ! Set pointer
        CALL lsyssc_getbase_Kdiagonal(RmatrixC(1), p_Kdiagonal)

        ! How many dimensions do we have?
        SELECT CASE(SIZE(RmatrixC,1))
        CASE (NDIM1D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)

        CASE (NDIM2D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          
      CASE (NDIM3D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)

        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarFCT')
          CALL sys_halt()
        END SELECT

      CASE DEFAULT
        CALL output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarFCT')
        CALL sys_halt()
      END SELECT
      
    CASE DEFAULT
      CALL output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarFCT')
      CALL sys_halt()
    END SELECT
    
    ! Release Ksep
    CALL storage_free(h_Ksep)
    
  CONTAINS

    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 2D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doInitFCTMat7_2D(Kld, Kcol, Ksep,&
        NEQ, NEDGE, NVAR, Cx, Cy, ML, u, uPredict,&
        dscale, theta, tstep, pp,pm, qp, qm, rp, rm,&
        flux, flux0, res, MC)
      
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER(PREC_MATIDX), INTENT(IN)                  :: NEDGE
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,ML
      REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL      :: MC
      REAL(DP), DIMENSION(NVAR,NEQ), INTENT(IN)         :: u,uPredict
      REAL(DP), INTENT(IN)                              :: dscale,theta,tstep
      REAL(DP), DIMENSION(NVAR,NEQ), INTENT(OUT)        :: pp,pm,qp,qm,rp,rm
      REAL(DP), DIMENSION(NVAR,NEDGE), INTENT(OUT)      :: flux,flux0
      REAL(DP), DIMENSION(NVAR,NEQ), INTENT(INOUT)      :: res

      REAL(DP), DIMENSION(NVAR*NVAR) :: T_ij
      REAL(DP), DIMENSION(NVAR)      :: F_ij,F_ji,G_ij,W_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji,iedge
      INTEGER(PREC_VECIDX)           :: i,j,ieq
      INTEGER                        :: ivar

      ! Clear nodal vectors
      CALL lalg_clearVectorDble2D(pp)
      CALL lalg_clearVectorDble2D(pm)
      CALL lalg_clearVectorDble2D(qp)
      CALL lalg_clearVectorDble2D(qm)

      ! Initialise edge counter
      iedge = 0

      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Increase edge counter
          iedge = iedge+1

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute the fluxes
          CALL fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
                   
          ! Assemble high-order residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji


          ! Compute the explicit antidiffusive fluxes and store it
          CALL fcb_calcRawFlux(u(:,i), u(:,j), C_ij, C_ji, tstep, F_ij)
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ij)-Cx(ji))
          C_ij(2) = 0.5_DP*(Cy(ij)-Cy(ji))

          ! Compute characteristic variables based on the explicit predictor
          CALL fcb_calcCharacteristics(uPredict(:,i), uPredict(:,j), C_ij, W_ij, L_ij=T_ij)

          W_ij=uPredict(:,j)-uPredict(:,i)

          ! Update the upper/lower bounds
          qp(:,i) = MAX(qp(:,i),W_ij); qp(:,j) = MAX(qp(:,j),-W_ij)
          qm(:,i) = MIN(qm(:,i),W_ij); qm(:,j) = MIN(qm(:,j),-W_ij)

          ! Transform fluxes into characteristic variables
!          CALL DGEMV('n', NVAR, NVAR, 1._DP, T_ij, NVAR, F_ij, 1, 0._DP, G_ij, 1)
          G_ij=F_ij
          flux0(:,iedge) = 0._DP
          flux (:,iedge) = G_ij

          ! Update the sums of positive/negative fluxes
          pp(:,i) = pp(:,i)+MAX(0._DP,G_ij); pp(:,j) = pp(:,j)+MAX(0._DP,-G_ij)
          pm(:,i) = pm(:,i)+MIN(0._DP,G_ij); pm(:,j) = pm(:,j)+MIN(0._DP,-G_ij)
        END DO
      END DO

      ! Adopt the explicit part (if required)
      IF (theta < 1.0_DP) THEN
        CALL lalg_vectorLinearCombDble2D(flux, flux0, 1.0_DP-theta, 1.0_DP)
      END IF

      ! Multiply upper/lower boundary by lumped mass matrix
      DO ieq = 1, NEQ
        qp(:,ieq) =  ML(ieq)*qp(:,ieq)
        qm(:,ieq) =  ML(ieq)*qm(:,ieq)
      END DO
      
      ! Apply the nodal limiter
      rp = afcstab_limit( pp, qp, 0._DP)
      rm = afcstab_limit(-pm,-qm, 0._DP)

      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1

          ! Get explicit raw antidiffusive flux
          G_ij=flux(:,iedge)

          ! Perform flux limiting
          DO ivar = 1, NVAR
            IF (G_ij(ivar) > 0._DP) THEN
              G_ij(ivar) = MIN(rp(ivar,i), rm(ivar,j))*G_ij(ivar)
            ELSE
              G_ij(ivar) = MIN(rm(ivar,i), rp(ivar,j))*G_ij(ivar)
            END IF
          END DO

          flux(:,iedge) = G_ij

          iedge = iedge-1
        END DO
      END DO
    END SUBROUTINE doInitFCTMat7_2D


    SUBROUTINE doLimitFCTMat7_2D(Kld, Kcol, Ksep,&
        NEQ, NEDGE, NVAR, Cx, Cy, u, uPredict, dscale, theta, tstep,&
        flux, flux0, res, MC)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER(PREC_MATIDX), INTENT(IN)                  :: NEDGE
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL      :: MC
      REAL(DP), DIMENSION(NVAR,NEQ), INTENT(IN)         :: u,uPredict
      REAL(DP), INTENT(IN)                              :: dscale,theta,tstep
      REAL(DP), DIMENSION(NVAR,NEDGE), INTENT(IN)       :: flux,flux0
      REAL(DP), DIMENSION(NVAR,NEQ), INTENT(INOUT)      :: res

      REAL(DP), DIMENSION(NVAR*NVAR) :: S_ij,T_ij
      REAL(DP), DIMENSION(NVAR)      :: F_ij,F_ji,G_ij,W_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji,iedge
      INTEGER(PREC_VECIDX)           :: i,j,ieq
      INTEGER                        :: ivar
      
      ! Initialise edge counter
      iedge = 0

      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Increase edge counter
          iedge = iedge+1

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute the fluxes
          CALL fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
                   
          ! Assemble high-order residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji

          
          ! Compute the implicit antidiffusive flux
          CALL fcb_calcRawFlux(u(:,i), u(:,j), C_ij, C_ji, tstep*theta, F_ij)
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ij)-Cx(ji))
          C_ij(2) = 0.5_DP*(Cy(ij)-Cy(ji))

          ! Compute characteristic variables based on the explicit predictor
!          CALL fcb_calcCharacteristics(uPredict(:,i), uPredict(:,j), C_ij, W_ij, R_ij=S_ij, L_ij=T_ij)

          ! Transform fluxes into characteristic variables
!          CALL DGEMV('n', NVAR, NVAR, 1._DP, T_ij, NVAR, F_ij, 1, 0._DP, G_ij, 1)

          G_ij=F_ij

          ! Assemble raw antidiffusive flux from explicit and implicit parts
          G_ij = flux0(:,iedge)+G_ij

          DO ivar = 1, NVAR
            IF (G_ij(ivar) > 0._DP) THEN
              G_ij(ivar) = MIN(G_ij(ivar), MAX(flux(ivar,iedge), 0._DP))
            ELSE
              G_ij(ivar) = MAX(G_ij(ivar), MIN(flux(ivar,iedge), 0._DP))
            END IF
          END DO

 !         CALL DGEMV('n', NVAR, NVAR, 1._DP, S_ij, NVAR, G_ij, 1, 0._DP, F_ij, 1)

          F_ij=G_ij

          ! Update the residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij

        END DO
      END DO
    END SUBROUTINE doLimitFCTMat7_2D
    
  END SUBROUTINE gfsys_buildResScalarFCT

  ! *****************************************************************************
  
!<subroutine>
  
  SUBROUTINE gfsys_buildResBlockTVD(RmatrixC, ru, fcb_calcFlux,&
      fcb_calcCharacteristics, rafcstab, dscale, bclear, rres)

!<description>
    ! This subroutine assembles the residual vector for FEM-TVD schemes.
    ! If the vectors contain only one block, then the scalar counterpart
    ! of this routine is called with the scalar subvectors.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: RmatrixC

    ! solution vector
    TYPE(t_vectorBlock), INTENT(IN)                :: ru

    ! scaling factor
    REAL(DP), INTENT(IN)                           :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    LOGICAL, INTENT(IN)                            :: bclear

    ! callback functions to compute local matrices
    INCLUDE 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)                 :: rafcstab

    ! residual vector
    TYPE(t_vectorBlock), INTENT(INOUT)             :: rres
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kcol
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kld,p_Ksep,p_Kdiagonal
    REAL(DP), DIMENSION(:), POINTER :: p_Cx,p_Cy,p_Cz,p_u,p_res
    REAL(DP), DIMENSION(:), POINTER :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    INTEGER :: h_Ksep

    
    ! Check if block vectors contain only one block.
    IF ((ru%nblocks .EQ. 1) .AND. (rres%nblocks .EQ. 1) ) THEN
      CALL gfsys_buildResScalarTVD(RmatrixC, ru%RvectorBlock(1),&
          fcb_calcFlux, fcb_calcCharacteristics,&
          rafcstab, dscale, bclear, rres%RvectorBlock(1))
      RETURN       
    END IF

    ! Check if vectors are compatible
    CALL lsysbl_isVectorCompatible(ru, rres)
    CALL gfsys_isVectorCompatible(rafcstab, ru)

    ! Clear vector?
    IF (bclear) CALL lsysbl_clearVector(rres)
    
    ! Set pointers
    CALL lsyssc_getbase_Kld   (RmatrixC(1), p_Kld)
    CALL lsyssc_getbase_Kcol  (RmatrixC(1), p_Kcol)
    CALL lsysbl_getbase_double(ru,   p_u)
    CALL lsysbl_getbase_double(rres, p_res)

    ! Create diagonal Ksep=Kld
    h_Ksep = ST_NOHANDLE
    CALL storage_copy(RmatrixC(1)%h_Kld, h_Ksep)
    CALL storage_getbase_int(h_Ksep, p_Ksep, RmatrixC(1)%NEQ+1)
       

    ! What kind of stabilisation should be applied?
    SELECT CASE(rafcstab%ctypeAFCstabilisation)
      
    CASE (AFCSTAB_FEMTVD)
      
      ! Check if stabilisation is prepeared
      IF (IAND(rafcstab%iSpec, AFCSTAB_INITIALISED) .EQ. 0) THEN
        CALL output_line('Stabilisation has not been initialised',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlockTVD')
        CALL sys_halt()
      END IF
      
      ! Set pointers
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
      
      ! Set specifiers for Ps, Qs and Rs
      rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_LIMITER)

       ! What kind of matrix are we?
      SELECT CASE(RmatrixC(1)%cmatrixFormat)
      CASE(LSYSSC_MATRIX7)
        
        ! How many dimensions do we have?
        SELECT CASE(SIZE(RmatrixC,1))
        CASE (NDIM1D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL doLimitTVDMat7_1D(p_Kld, p_Kcol, p_Ksep, RmatrixC(1)%NEQ, ru%nblocks,&
              p_Cx, p_u, dscale, p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
          
        CASE (NDIM2D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL doLimitTVDMat7_2D(p_Kld, p_Kcol, p_Ksep, RmatrixC(1)%NEQ, ru%nblocks,&
              p_Cx, p_Cy, p_u, dscale, p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
          
        CASE (NDIM3D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
          CALL doLimitTVDMat7_3D(p_Kld, p_Kcol, p_Ksep, RmatrixC(1)%NEQ, ru%nblocks,&
              p_Cx, p_Cy, p_Cz, p_u, dscale, p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
          
        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlockTVD')
          CALL sys_halt()
        END SELECT
        
        
      CASE (LSYSSC_MATRIX9)
        
        ! Set pointer
        CALL lsyssc_getbase_Kdiagonal(RmatrixC(1), p_Kdiagonal)

        ! How many dimensions do we have?
        SELECT CASE(SIZE(RmatrixC,1))
        CASE (NDIM1D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL doLimitTVDMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, dscale,&
              p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
          
        CASE (NDIM2D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL doLimitTVDMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, dscale,&
              p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
          
        CASE (NDIM3D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
          CALL doLimitTVDMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, dscale,&
              p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
          
        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlockTVD')
          CALL sys_halt()
        END SELECT
        
      CASE DEFAULT
        CALL output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlockTVD')
        CALL sys_halt()
      END SELECT
      
    CASE DEFAULT
      CALL output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlockTVD')
      CALL sys_halt()
    END SELECT
    
    ! Release Ksep
    CALL storage_free(h_Ksep)

  CONTAINS

    ! Here, the working routines follow
    
    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 1D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doLimitTVDMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, u, dscale, pp, pm, qp, qm, rp, rm, res)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: pp,pm,qp,qm,rp,rm
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(INOUT)      :: res

      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,R_ij,L_ij
      REAL(DP), DIMENSION(NVAR)      :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j,iloc,jloc
      INTEGER                        :: ivar
      

      ! Clear P's and Q's
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0
      
      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute the fluxes
          CALL fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)
           
          ! Assemble high-order residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u_i, u_j, XDir1D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR

            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO
          
        END DO
      END DO

      ! Compute nodal correction factors for X-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)
      
      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u_i, u_j, XDir1D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij   =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij
        END DO
      END DO
    END SUBROUTINE doLimitTVDMat7_1D


    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 2D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doLimitTVDMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, Cy, u, dscale, pp, pm, qp, qm, rp, rm, res)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: pp,pm,qp,qm,rp,rm
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(INOUT)      :: res

      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,R_ij,L_ij
      REAL(DP), DIMENSION(NVAR)      :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j,iloc,jloc
      INTEGER                        :: ivar
      
      ! Clear P's and Q's
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0

      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute the fluxes
          CALL fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)
                   
          ! Assemble high-order residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u_i, u_j, XDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR

            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO
          
        END DO
      END DO
      
      ! Compute nodal correction factors for X-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)
      
      ! Clear P's and Q's for Y-direction
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0

      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u_i, u_j, XDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij   =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij

          ! Compute characteristic fluxes in Y-direction
          CALL fcb_calcCharacteristics(u_i, u_j, YDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ij)-Cy(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR

            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO
          
        END DO
      END DO

      ! Compute nodal correction factors for Y-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)
      
      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute characteristic fluxes in Y-direction
          CALL fcb_calcCharacteristics(u_i, u_j, YDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ij)-Cy(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij  =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar =1, NVAR

            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij
        END DO
      END DO
    END SUBROUTINE doLimitTVDMat7_2D

    
    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 3D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doLimitTVDMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, Cy, Cz, u, dscale, pp, pm, qp, qm, rp, rm, res)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: pp,pm,qp,qm,rp,rm
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(INOUT)      :: res

      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,R_ij,L_ij
      REAL(DP), DIMENSION(NVAR)      :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j,iloc,jloc
      INTEGER                        :: ivar
      
      ! Clear P's and Q's
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0
      
      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)
          
          ! Compute the fluxes
          CALL fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble high-order residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u_i, u_j, XDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR

            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO
          
        END DO
      END DO

      ! Compute nodal correction factors for X-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)

      ! Clear P's and Q's for Y-direction
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0

      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u_i, u_j, XDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij   =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij

          ! Compute characteristic fluxes in Y-direction
          CALL fcb_calcCharacteristics(u_i, u_j, YDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ij)-Cy(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR

            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO
          
        END DO
      END DO

      ! Compute nodal correction factors for Y-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)

      ! Clear P's and Q's for Z-direction
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0

      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute characteristic fluxes in Y-direction
          CALL fcb_calcCharacteristics(u_i, u_j, YDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ij)-Cy(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij  =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar =1, NVAR

            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij

          ! Compute characteristic fluxes in Z-direction
          CALL fcb_calcCharacteristics(u_i, u_j, ZDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cz(ij)-Cz(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cz(ij)+Cz(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij

          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR
            
            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO

        END DO
      END DO

      ! Compute nodal correction factors for Z-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)

      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute characteristic fluxes in Z-direction
          CALL fcb_calcCharacteristics(u_i, u_j, ZDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cz(ij)-Cz(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cz(ij)-Cz(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij  =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar =1, NVAR

            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij
        END DO
      END DO
    END SUBROUTINE doLimitTVDMat7_3D
    
    
    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 1D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doLimitTVDMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, u, dscale, pp, pm, qp, qm, rp, rm, res)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: pp,pm,qp,qm,rp,rm
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(INOUT)      :: res

      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,R_ij,L_ij
      REAL(DP), DIMENSION(NVAR)      :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j,iloc,jloc
      INTEGER                        :: ivar
      

      ! Clear P's and Q's
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0
      
      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute the fluxes
          CALL fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble high-order residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u_i, u_j, XDir1D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR

            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO
          
        END DO
      END DO

      ! Compute nodal correction factors for X-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)
      
      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); Ksep(j) = Ksep(j)-1; ji = Ksep(j)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u_i, u_j, XDir1D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij   =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij
        END DO
      END DO
    END SUBROUTINE doLimitTVDMat9_1D


    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 2D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doLimitTVDMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, u, dscale, pp, pm, qp, qm, rp, rm, res)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: pp,pm,qp,qm,rp,rm
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(INOUT)      :: res

      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,R_ij,L_ij
      REAL(DP), DIMENSION(NVAR)      :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j,iloc,jloc
      INTEGER                        :: ivar

      ! Clear P's and Q's
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0
      
      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute the fluxes
          CALL fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble high-order residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u_i, u_j, XDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR

            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO
          
        END DO
      END DO
      
      ! Compute nodal correction factors for X-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)
      
      ! Clear P's and Q's for Y-direction
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0

      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); Ksep(j) = Ksep(j)-1; ji = Ksep(j)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u_i, u_j, XDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij   =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij

          ! Compute characteristic fluxes in Y-direction
          CALL fcb_calcCharacteristics(u_i, u_j, YDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ij)-Cy(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR

            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO
          
        END DO
      END DO

      ! Compute nodal correction factors for Y-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)
      
      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute characteristic fluxes in Y-direction
          CALL fcb_calcCharacteristics(u_i, u_j, YDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ij)-Cy(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij  =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar =1, NVAR

            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij
        END DO
      END DO
    END SUBROUTINE doLimitTVDMat9_2D

    
    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 3D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doLimitTVDMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, Cz, u, dscale, pp, pm, qp, qm, rp, rm, res)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: pp,pm,qp,qm,rp,rm
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(INOUT)      :: res

      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,R_ij,L_ij
      REAL(DP), DIMENSION(NVAR)      :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j,iloc,jloc
      INTEGER                        :: ivar
      
      ! Clear P's and Q's
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0
      
      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)
          
          ! Compute the fluxes
          CALL fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)
                 
          ! Assemble high-order residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u_i, u_j, XDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR

            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO
          
        END DO
      END DO

      ! Compute nodal correction factors for X-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)

      ! Clear P's and Q's for Y-direction
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0

      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); Ksep(j) = Ksep(j)-1; ji = Ksep(j)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
                    
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u_i, u_j, XDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij   =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij

          ! Compute characteristic fluxes in Y-direction
          CALL fcb_calcCharacteristics(u_i, u_j, YDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ij)-Cy(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR

            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO
          
        END DO
      END DO

      ! Compute nodal correction factors for Y-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)

      ! Clear P's and Q's for Z-direction
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0

      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute characteristic fluxes in Y-direction
          CALL fcb_calcCharacteristics(u_i, u_j, YDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ij)-Cy(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij  =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar =1, NVAR

            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij

          ! Compute characteristic fluxes in Z-direction
          CALL fcb_calcCharacteristics(u_i, u_j, ZDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cz(ij)-Cz(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cz(ij)+Cz(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij

          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR
            
            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO

        END DO
      END DO

      ! Compute nodal correction factors for Z-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)

      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); Ksep(j) = Ksep(j)-1; ji = Ksep(j)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute characteristic fluxes in Z-direction
          CALL fcb_calcCharacteristics(u_i, u_j, ZDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cz(ij)-Cz(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cz(ij)+Cz(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij  =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar =1, NVAR

            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij
        END DO
      END DO
    END SUBROUTINE doLimitTVDMat9_3D
  END SUBROUTINE gfsys_buildResBlockTVD

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE gfsys_buildResScalarTVD(RmatrixC, ru, fcb_calcFlux,&
      fcb_calcCharacteristics, rafcstab, dscale, bclear, rres)
    
!<description>
    ! This subroutine assembles the residual vector for FEM-TVD schemes
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: RmatrixC

    ! solution vector
    TYPE(t_vectorScalar), INTENT(IN)               :: ru

    ! scaling factor
    REAL(DP), INTENT(IN)                           :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    LOGICAL, INTENT(IN)                            :: bclear

    ! callback functions to compute local matrices
    INCLUDE 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)                 :: rafcstab
    
    ! residual vector
    TYPE(t_vectorScalar), INTENT(INOUT)            :: rres
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kcol
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kld,p_Ksep,p_Kdiagonal
    REAL(DP), DIMENSION(:), POINTER :: p_Cx,p_Cy,p_Cz,p_u,p_res
    REAL(DP), DIMENSION(:), POINTER :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    INTEGER :: h_Ksep


    ! Check if vectors are compatible
    CALL lsyssc_isVectorCompatible(ru, rres)
    CALL gfsys_isVectorCompatible(rafcstab, ru)

    ! Clear vector?
    IF (bclear) CALL lsyssc_clearVector(rres)
    
    ! Set pointers
    CALL lsyssc_getbase_Kld   (RmatrixC(1), p_Kld)
    CALL lsyssc_getbase_Kcol  (RmatrixC(1), p_Kcol)
    CALL lsyssc_getbase_double(ru,   p_u)
    CALL lsyssc_getbase_double(rres, p_res)

    ! Create diagonal Ksep=Kld
    h_Ksep = ST_NOHANDLE
    CALL storage_copy(RmatrixC(1)%h_Kld, h_Ksep)
    CALL storage_getbase_int(h_Ksep, p_Ksep, RmatrixC(1)%NEQ+1)
       

    ! What kind of stabilisation should be applied?
    SELECT CASE(rafcstab%ctypeAFCstabilisation)
      
    CASE (AFCSTAB_FEMTVD)
      
      ! Check if stabilisation is prepeared
      IF (IAND(rafcstab%iSpec, AFCSTAB_INITIALISED) .EQ. 0) THEN
        CALL output_line('Stabilisation has not been initialised',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarTVD')
        CALL sys_halt()
      END IF
      
      ! Set pointers
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
      
      ! Set specifiers for Ps, Qs and Rs
      rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_LIMITER)
      
      ! What kind of matrix are we?
      SELECT CASE(RmatrixC(1)%cmatrixFormat)
      CASE(LSYSSC_MATRIX7)
        
        ! How many dimensions do we have?
        SELECT CASE(SIZE(RmatrixC,1))
        CASE (NDIM1D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL doLimitTVDMat7_1D(p_Kld, p_Kcol, p_Ksep, RmatrixC(1)%NEQ, ru%NVAR,&
              p_Cx, p_u, dscale, p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
          
        CASE (NDIM2D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL doLimitTVDMat7_2D(p_Kld, p_Kcol, p_Ksep, RmatrixC(1)%NEQ, ru%NVAR,&
              p_Cx, p_Cy, p_u, dscale, p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
          
        CASE (NDIM3D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
          CALL doLimitTVDMat7_3D(p_Kld, p_Kcol, p_Ksep, RmatrixC(1)%NEQ, ru%NVAR,&
              p_Cx, p_Cy, p_Cz, p_u, dscale, p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
          
        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarTVD')
          CALL sys_halt()
        END SELECT
        
        
      CASE (LSYSSC_MATRIX9)
        
        ! Set pointer
        CALL lsyssc_getbase_Kdiagonal(RmatrixC(1), p_Kdiagonal)

        ! How many dimensions do we have?
        SELECT CASE(SIZE(RmatrixC,1))
        CASE (NDIM1D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL doLimitTVDMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, dscale,&
              p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
          
        CASE (NDIM2D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)

          CALL doLimitTVDMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, dscale,&
              p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
          
        CASE (NDIM3D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
          CALL doLimitTVDMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, dscale,&
              p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
          
        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarTVD')
          CALL sys_halt()
        END SELECT
        
      CASE DEFAULT
        CALL output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarTVD')
        CALL sys_halt()
      END SELECT
      
    CASE DEFAULT
      CALL output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarTVD')
      CALL sys_halt()
    END SELECT
    
    ! Release Ksep
    CALL storage_free(h_Ksep)
    
  CONTAINS

    ! Here, the working routines follow
    
    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 1D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doLimitTVDMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, u, dscale, pp, pm, qp, qm, rp, rm, res)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: pp,pm,qp,qm,rp,rm
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: res

      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,R_ij,L_ij
      REAL(DP), DIMENSION(NVAR)      :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j,iloc,jloc
      INTEGER                        :: ivar
      

      ! Clear P's and Q's
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0
      
      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute the fluxes
          CALL fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble high-order residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), XDir1D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR

            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO
          
        END DO
      END DO

      ! Compute nodal correction factors for X-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)
      
      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1
          
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), XDir1D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij   =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij
        END DO
      END DO
    END SUBROUTINE doLimitTVDMat7_1D


    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 2D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doLimitTVDMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, Cy, u, dscale, pp, pm, qp, qm, rp, rm, res)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: pp,pm,qp,qm,rp,rm
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: res

      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,R_ij,L_ij
      REAL(DP), DIMENSION(NVAR)      :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j,iloc,jloc
      INTEGER                        :: ivar
      

      ! Clear P's and Q's
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0
      
      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Compute the fluxes
          CALL fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble high-order residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), XDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR

            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO
          
        END DO
      END DO
      
      ! Compute nodal correction factors for X-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)
      
      ! Clear P's and Q's for Y-direction
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0

      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1
                    
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), XDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij   =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij

          ! Compute characteristic fluxes in Y-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), YDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ij)-Cy(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR

            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO
          
        END DO
      END DO

      ! Compute nodal correction factors for Y-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)
      
      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute characteristic fluxes in Y-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), YDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ij)-Cy(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij  =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar =1, NVAR

            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij
        END DO
      END DO
    END SUBROUTINE doLimitTVDMat7_2D


    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 3D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doLimitTVDMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, Cy, Cz, u, dscale, pp, pm, qp, qm, rp, rm, res)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: pp,pm,qp,qm,rp,rm
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: res

      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,R_ij,L_ij
      REAL(DP), DIMENSION(NVAR)      :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j,iloc,jloc
      INTEGER                        :: ivar
      

      ! Clear P's and Q's
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0
      
      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
              
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Compute the fluxes
          CALL fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
                    
          ! Assemble high-order residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), XDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR

            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO
          
        END DO
      END DO

      ! Compute nodal correction factors for X-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)

      ! Clear P's and Q's for Y-direction
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0

      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1
          
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), XDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij   =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij

          ! Compute characteristic fluxes in Y-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), YDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ij)-Cy(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR

            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO
          
        END DO
      END DO

      ! Compute nodal correction factors for Y-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)

      ! Clear P's and Q's for Z-direction
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0

      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute characteristic fluxes in Y-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), YDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ij)-Cy(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij  =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar =1, NVAR

            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij

          ! Compute characteristic fluxes in Z-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), ZDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cz(ij)-Cz(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cz(ij)+Cz(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij

          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR
            
            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO

        END DO
      END DO

      ! Compute nodal correction factors for Z-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)

      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1

          ! Compute characteristic fluxes in Z-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), ZDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cz(ij)-Cz(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cz(ij)+Cz(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij  =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar =1, NVAR

            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij
        END DO
      END DO
    END SUBROUTINE doLimitTVDMat7_3D


    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 1D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doLimitTVDMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, u, dscale, pp, pm, qp, qm, rp, rm, res)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: pp,pm,qp,qm,rp,rm
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: res

      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,R_ij,L_ij
      REAL(DP), DIMENSION(NVAR)      :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j,iloc,jloc
      INTEGER                        :: ivar
      

      ! Clear P's and Q's
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0
      
      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
              
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Compute the fluxes
          CALL fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble high-order residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), XDir1D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR

            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO
          
        END DO
      END DO

      ! Compute nodal correction factors for X-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)
      
      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); Ksep(j) = Ksep(j)-1; ji = Ksep(j)
          
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), XDir1D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij   =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij
        END DO
      END DO
    END SUBROUTINE doLimitTVDMat9_1D

    
    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 2D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doLimitTVDMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, u, dscale, pp, pm, qp, qm, rp, rm, res)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: pp,pm,qp,qm,rp,rm
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: res

      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,R_ij,L_ij
      REAL(DP), DIMENSION(NVAR)      :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j,iloc,jloc
      INTEGER                        :: ivar

      ! Clear P's and Q's
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0
      
      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Compute the fluxes
          CALL fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble high-order residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), XDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR

            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO
          
        END DO
      END DO
      
      ! Compute nodal correction factors for X-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)

      ! Clear P's and Q's for Y-direction
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0

      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); Ksep(j) = Ksep(j)-1; ji = Ksep(j)
          
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), XDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij   =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar = 1, NVAR

            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij

          ! Compute characteristic fluxes in Y-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), YDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ij)-Cy(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR

            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO
          
        END DO
      END DO

      ! Compute nodal correction factors for Y-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)

      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute characteristic fluxes in Y-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), YDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ij)-Cy(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij  =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar =1, NVAR

            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij
        END DO
      END DO
    END SUBROUTINE doLimitTVDMat9_2D


    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 3D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doLimitTVDMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, Cz, u, dscale, pp, pm, qp, qm, rp, rm, res)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: pp,pm,qp,qm,rp,rm
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: res

      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,R_ij,L_ij
      REAL(DP), DIMENSION(NVAR)      :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij,C_ji
      INTEGER(PREC_MATIDX)           :: ij,ji
      INTEGER(PREC_VECIDX)           :: i,j,iloc,jloc
      INTEGER                        :: ivar
      

      ! Clear P's and Q's
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0
      
      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Compute the fluxes
          CALL fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble high-order residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), XDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR

            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO
          
        END DO
      END DO

      ! Compute nodal correction factors for X-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)

      ! Clear P's and Q's for Y-direction
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0

      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); Ksep(j) = Ksep(j)-1; ji = Ksep(j)
          
          ! Compute characteristic fluxes in X-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), XDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ij)-Cx(ji))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij   =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij

          ! Compute characteristic fluxes in Y-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), YDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ij)-Cy(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR

            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO
          
        END DO
      END DO

      ! Compute nodal correction factors for Y-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)

      ! Clear P's and Q's for Z-direction
      pp(:,1:NEQ)=0; pm(:,1:NEQ)=0
      qp(:,1:NEQ)=0; qm(:,1:NEQ)=0

      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute characteristic fluxes in Y-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), YDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ij)-Cy(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij  =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar =1, NVAR

            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij

          ! Compute characteristic fluxes in Z-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), ZDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cz(ij)-Cz(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cz(ij)+Cz(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij

          ! Update sums of downstream/upstream edge contributions
          DO ivar = 1, NVAR
            
            ! Set node orientation
            IF (ka_ij(ivar) > 0) THEN
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            ELSE
              iloc = i; jloc = j
            END IF
            
            ! Assemble P's and Q's
            IF (F_ij(ivar) > 0) THEN
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            ELSE
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            END IF
          END DO

        END DO
      END DO

      ! Compute nodal correction factors for Z-direction
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP, 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP, 1._DP)

      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); Ksep(j) = Ksep(j)-1; ji = Ksep(j)

          ! Compute characteristic fluxes in Z-direction
          CALL fcb_calcCharacteristics(u(:,i), u(:,j), ZDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cz(ij)-Cz(ji))*Lbd_ij
          ks_ij =  0.5_DP*(Cz(ij)+Cz(ji))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-ks_ij, 2._DP*ABS(ka_ij)))*W_ij
          W_ij  =  ABS(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          DO ivar =1, NVAR

            ! Limit characteristic fluxes
            IF (ka_ij(ivar) < 0) THEN
              IF (F_ij(ivar) > 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              END IF
            ELSE
              IF (F_ij(ivar) < 0) THEN
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              ELSE
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              END IF
            END IF
          END DO
          
          ! Transform back into conservative variables
          CALL DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0._DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij
        END DO
      END DO
    END SUBROUTINE doLimitTVDMat9_3D
  END SUBROUTINE gfsys_buildResScalarTVD

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsys_buildDivJacobianBlock(RmatrixC, ru,&
      fcb_calcMatrix, hstep, bclear, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix.
    !
    ! Note that this routine is designed for block matrices/vectors. 
    ! If there is only one block, then the corresponding scalar routine 
    ! is called. Otherwise, the global operator is treated as block matrix.
    ! This block matrix has to be in group structure, that is, the structure
    ! of subblock(1,1) will serve as template for all other submatrices.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: RmatrixC

    ! solution vector
    TYPE(t_vectorBlock), INTENT(IN)                :: ru
    
    ! perturbation parameter
    REAL(DP), INTENT(IN)                           :: hstep

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    LOGICAL, INTENT(IN)                            :: bclear

    ! callback functions to compute velocity
    INCLUDE 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! Jacobian matrix
    TYPE(t_matrixBlock), INTENT(INOUT)             :: rmatrixJ
!</inputoutput>
!</subroutine>

    ! Check if block vector contains only one block and if
    ! global operator is stored in interleave format.
    IF ((ru%nblocks           .EQ. 1) .AND.&
        (rmatrixJ%ndiagblocks .EQ. 1)) THEN
      CALL gfsys_buildDivJacobianScalar(RmatrixC, ru%RvectorBlock(1),&
          fcb_calcMatrix, hstep, bclear, rmatrixJ%RmatrixBlock(1,1))
      RETURN       
    END IF

    PRINT *, "Jacobian matrix for systems is not yet implemented!"
    STOP
  END SUBROUTINE gfsys_buildDivJacobianBlock

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsys_buildDivJacobianScalar(RmatrixC, ru,&
      fcb_calcMatrix, hstep, bclear, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix.
    !
    ! Note that this routine requires scalar matrices stored in the interleave
    ! format.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: RmatrixC
    
    ! solution vector
    TYPE(t_vectorScalar), INTENT(IN)               :: ru
    
    ! perturbation parameter
    REAL(DP), INTENT(IN)                           :: hstep

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    LOGICAL, INTENT(IN)                            :: bclear

    ! callback functions to compute velocity
    INCLUDE 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! Jacobian matrix
    TYPE(t_matrixScalar), INTENT(INOUT)            :: rmatrixJ
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld,p_Ksep
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kdiagonal
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
    REAL(DP), DIMENSION(:), POINTER             :: p_Cx,p_Cy,p_Cz,p_J,p_u
    INTEGER :: h_Ksep
    INTEGER :: idim,ndim

    ! Check if all matrices/vectors are compatible
    CALL lsyssc_isMatrixCompatible(ru, rmatrixJ, .FALSE.)

    ndim = SIZE(RmatrixC,1)
    DO idim = 1, ndim
      CALL lsyssc_isMatrixCompatible(rmatrixC(idim), rmatrixJ)
    END DO
        
    ! Clear matrix?
    IF (bclear) CALL lsyssc_clearMatrix(rmatrixJ)
    
    ! Set pointers
    CALL lsyssc_getbase_Kld   (RmatrixC(1), p_Kld)
    CALL lsyssc_getbase_Kcol  (RmatrixC(1), p_Kcol)
    CALL lsyssc_getbase_double(rmatrixJ, p_J)
    CALL lsyssc_getbase_double(ru, p_u)

    ! How many dimensions do we have?
    SELECT CASE(ndim)
    CASE (NDIM1D)
      CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)

    CASE (NDIM2D)
      CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)
      CALL lsyssc_getbase_double(rmatrixC(2), p_Cy)

    CASE (NDIM3D)
      CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)
      CALL lsyssc_getbase_double(rmatrixC(2), p_Cy)
      CALL lsyssc_getbase_double(rmatrixC(3), p_Cz)

    CASE DEFAULT
      CALL output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivJacobianScalar')
      CALL sys_halt()
    END SELECT
    
    ! What kind of matrix are we?
    SELECT CASE(rmatrixJ%cmatrixFormat)
    CASE(LSYSSC_MATRIX7INTL)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      CALL storage_copy(RmatrixC(1)%h_Kld, h_Ksep)
      CALL storage_getbase_int(h_Ksep, p_Ksep, RmatrixC(1)%NEQ+1)
      
      ! What type of matrix are we?
      SELECT CASE(rmatrixJ%cinterleavematrixFormat)
        
      CASE (LSYSSC_MATRIX1)
        
        SELECT CASE(ndim)
!!$        CASE (NDIM1D)
!!$          CALL doJacobianMat7_1D(p_Kld, p_Kcol, p_Ksep,&
!!$              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, p_J)
!!$
!!$        CASE (NDIM2D)
!!$          CALL doJacobianMat7_2D(p_Kld, p_Kcol, p_Ksep,&
!!$              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, p_J)
!!$
!!$        CASE (NDIM3D)
!!$          CALL doJacobianMat7_3D(p_Kld, p_Kcol, p_Ksep,&
!!$              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, p_J)
          
        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivJacobianScalar')
          CALL sys_halt()
        END SELECT
        
      CASE (LSYSSC_MATRIXD)
        
        SELECT CASE(ndim)
!!$        CASE (NDIM1D)
!!$          CALL doJacobianMat7Diag_1D(p_Kld, p_Kcol, p_Ksep,&
!!$              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, p_J)
!!$
!!$        CASE (NDIM2D)
!!$          CALL doJacobianMat7Diag_2D(p_Kld, p_Kcol, p_Ksep,&
!!$              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, p_J)
!!$
!!$        CASE (NDIM3D)
!!$          CALL doJacobianMat7Diag_3D(p_Kld, p_Kcol, p_Ksep,&
!!$              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, p_J)

        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivJacobianScalar')
          CALL sys_halt()
        END SELECT
        
      CASE DEFAULT
        CALL output_line('Unsupported interleave matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivJacobianScalar')
        CALL sys_halt()
      END SELECT
      
      ! Release diagonal separator
      CALL storage_free(h_Ksep)
      
      
    CASE (LSYSSC_MATRIX9INTL)
      
      ! Set pointers
      CALL lsyssc_getbase_Kdiagonal(rmatrixJ, p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      CALL storage_copy(RmatrixC(1)%h_Kld, h_Ksep)
      CALL storage_getbase_int(h_Ksep, p_Ksep, RmatrixC(1)%NEQ+1)
      
      ! What type of matrix are we?
      SELECT CASE(rmatrixJ%cinterleavematrixFormat)
        
      CASE (LSYSSC_MATRIX1)
        
        SELECT CASE(ndim)
!!$        CASE (NDIM1D)
!!$          CALL doJacobianMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
!!$              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, hstep, p_J)
!!$
!!$        CASE (NDIM2D)
!!$          CALL doJacobianMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
!!$              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, hstep, p_J)
!!$
!!$        CASE (NDIM3D)
!!$          CALL doJacobianMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
!!$              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, hstep, p_J)

        CASE DEFAULT
          CALL output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivJacobianScalar')
          CALL sys_halt() 
        END SELECT
        
      CASE (LSYSSC_MATRIXD)
        
        SELECT CASE(ndim)
!!$        CASE (NDIM1D)
!!$          CALL doJacobianMat9Diag_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
!!$              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, hstep, p_J)
!!$
!!$        CASE (NDIM2D)
!!$          CALL doJacobianMat9Diag_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
!!$              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, hstep, p_J)
!!$
!!$        CASE (NDIM3D)
!!$          CALL doJacobianMat9Diag_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
!!$              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, hstep, p_J)
          
        CASE DEFAULT
          CALL output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivJacobianScalar')
          CALL sys_halt()
        END SELECT
        
      CASE DEFAULT
        CALL output_line('Unsupported interleave matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivJacobianScalar')
        CALL sys_halt()
      END SELECT
      
      ! Release diagonal separator
      CALL storage_free(h_Ksep)
      
      
    CASE DEFAULT
      CALL output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivJacobianScalar')
      CALL sys_halt()
    END SELECT

  END SUBROUTINE gfsys_buildDivJacobianScalar

  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsys_getbase_double(rmatrix, rarray, bisFullMatrix)

!<description>
    ! This subroutine assigns the pointers of the array to the scalar
    ! submatrices of the block matrix.
    ! If the optional parameter bisFullMatrix is given, then this routine
    ! returns bisFullMatrix = .TRUE. if all blocks of rmatrix are associated.
    ! Otherwise, bisFullMatrix = .FALSE. is returned if only the diagonal
    ! blocks of the block matrix are associated.
!</description>

!<input>
    ! The block matrix
    TYPE(t_matrixBlock), INTENT(IN)            :: rmatrix
!</input>

!<output>
    ! The array
    TYPE(t_array), DIMENSION(:,:), INTENT(OUT) :: rarray

    ! OPTIONAL: indicator for full block matrix
    LOGICAL, INTENT(OUT), OPTIONAL             :: bisFullMatrix
!</output>
!</subroutine>
    
    ! local variables
    INTEGER :: iblock,jblock
    LOGICAL :: bisFull
    
    ! Check if array is compatible
    IF (rmatrix%ndiagblocks .NE. SIZE(rarray,1) .OR.&
        rmatrix%ndiagblocks .NE. SIZE(rarray,2)) THEN
      CALL output_line('Block matrix and array are not compatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_getbase_double')
      CALL sys_halt()
    END IF

    ! Assign pointers
    bisFull = .TRUE.
    DO iblock = 1, rmatrix%ndiagblocks
      DO jblock = 1, rmatrix%ndiagblocks
        IF (lsyssc_isExplicitMatrix1D(rmatrix%RmatrixBlock(iblock,jblock))) THEN
          CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(iblock,jblock),&
                                     rarray(iblock,jblock)%Da)
        ELSE
          NULLIFY(rarray(iblock,jblock)%Da)
          bisFull = .FALSE.
        END IF
      END DO
    END DO

    IF (PRESENT(bisFullMatrix)) bisFullMatrix = bisFull
  END SUBROUTINE gfsys_getbase_double

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsys_getbase_single(rmatrix, rarray, bisFullMatrix)

!<description>
    ! This subroutine assigns the pointers of the array to the scalar
    ! submatrices of the block matrix.
    ! If the optional parameter bisFullMatrix is given, then this routine
    ! returns bisFullMatrix = .TRUE. if all blocks of rmatrix are associated.
    ! Otherwise, bisFullMatrix = .FALSE. is returned if only the diagonal
    ! blocks of the block matrix are associated.
!</description>

!<input>
    ! The block matrix
    TYPE(t_matrixBlock), INTENT(IN)            :: rmatrix
!</input>

!<output>
    ! The array
    TYPE(t_array), DIMENSION(:,:), INTENT(OUT) :: rarray

    ! OPTIONAL: indicator for full block matrix
    LOGICAL, INTENT(OUT), OPTIONAL             :: bisFullMatrix
!</output>
!</subroutine>
    
    ! local variables
    INTEGER :: iblock,jblock
    LOGICAL :: bisFull
    
    ! Check if array is compatible
    IF (rmatrix%ndiagblocks .NE. SIZE(rarray,1) .OR.&
        rmatrix%ndiagblocks .NE. SIZE(rarray,2)) THEN
      CALL output_line('Block matrix and array are not compatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_getbase_single')
      CALL sys_halt()
    END IF

    ! Assign pointers
    bisFull = .TRUE.
    DO iblock = 1, rmatrix%ndiagblocks
      DO jblock = 1, rmatrix%ndiagblocks
        IF (lsyssc_isExplicitMatrix1D(rmatrix%RmatrixBlock(iblock,jblock))) THEN
          CALL lsyssc_getbase_single(rmatrix%RmatrixBlock(iblock,jblock),&
                                     rarray(iblock,jblock)%Fa)
        ELSE
          NULLIFY(rarray(iblock,jblock)%Fa)
          bisFull = .FALSE.
        END IF
      END DO
    END DO

    IF (PRESENT(bisFullMatrix)) bisFullMatrix = bisFull
  END SUBROUTINE gfsys_getbase_single

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsys_getbase_int(rmatrix, rarray, bisFullMatrix)

!<description>
    ! This subroutine assigns the pointers of the array to the scalar
    ! submatrices of the block matrix.
    ! If the optional parameter bisFullMatrix is given, then this routine
    ! returns bisFullMatrix = .TRUE. if all blocks of rmatrix are associated.
    ! Otherwise, bisFullMatrix = .FALSE. is returned if only the diagonal
    ! blocks of the block matrix are associated.
!</description>

!<input>
    ! The block matrix
    TYPE(t_matrixBlock), INTENT(IN)            :: rmatrix
!</input>

!<output>
    ! The array
    TYPE(t_array), DIMENSION(:,:), INTENT(OUT) :: rarray

    ! OPTIONAL: indicator for full block matrix
    LOGICAL, INTENT(OUT), OPTIONAL             :: bisFullMatrix
!</output>
!</subroutine>
    
    ! local variables
    INTEGER :: iblock,jblock
    LOGICAL :: bisFull
    
    ! Check if array is compatible
    IF (rmatrix%ndiagblocks .NE. SIZE(rarray,1) .OR.&
        rmatrix%ndiagblocks .NE. SIZE(rarray,2)) THEN
      CALL output_line('Block matrix and array are not compatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_getbase_int')
      CALL sys_halt()
    END IF

    ! Assign pointers
    bisFull = .TRUE.
    DO iblock = 1, rmatrix%ndiagblocks
      DO jblock = 1, rmatrix%ndiagblocks
        IF (lsyssc_isExplicitMatrix1D(rmatrix%RmatrixBlock(iblock,jblock))) THEN
          CALL lsyssc_getbase_int(rmatrix%RmatrixBlock(iblock,jblock),&
                                  rarray(iblock,jblock)%Ia)
        ELSE
          NULLIFY(rarray(iblock,jblock)%Da)
          bisFull = .FALSE.
        END IF
      END DO
    END DO

    IF (PRESENT(bisFullMatrix)) bisFullMatrix = bisFull
  END SUBROUTINE gfsys_getbase_int
END MODULE groupfemsystem
