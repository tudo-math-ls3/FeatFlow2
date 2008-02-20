!##############################################################################
!# ****************************************************************************
!# <name> groupfemsystem </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic routines for applying the group-finite
!# element formulation to systems of conservation laws.
!# The group finite element formulation introduced by Fletcher uses the
!# same basis functions for the unknown solution and the fluxes.
!# Moreover, this module allows to modifying discrete operators by means of
!# the algebraic flux correction (AFC) methodology by Kuzmin, Möller and Turek.
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
!# 5.) gfsys_buildResidualGalerkin = gfsys_buildResScalarGalerkin /
!#                                   gfsys_buildResBlockGalerkin
!#     -> assemble the residual for standard Galerkin scheme
!#
!# 6.) gfsys_buildResidualUpwind = gfsys_buildResScalarUpwind /
!#                                 gfsys_buildResBlockUpwind
!#     -> assemble the residual for discrete upwinding scheme
!#
!# 7.) gfsys_buildResidualTVD = gfsys_buildResScalarTVD /
!#                              gfsys_buildResBlockTVD
!#     -> assemble the residual for FEM-TVD stabilisation
!#
!# 6.) gfsys_buildDivJacobian = gfsys_buildDivJacobianScalar /
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
  PUBLIC :: gfsys_buildResidualGalerkin
  PUBLIC :: gfsys_buildResidualUpwind
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

  INTERFACE gfsys_buildResidualGalerkin
    MODULE PROCEDURE gfsys_buildResScalarGalerkin
    MODULE PROCEDURE gfsys_buildResBlockGalerkin
  END INTERFACE

  INTERFACE gfsys_buildResidualUpwind
    MODULE PROCEDURE gfsys_buildResScalarUpwind
    MODULE PROCEDURE gfsys_buildResBlockUpwind
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

!<constants>
!<constantblock description="Global format flags for convective operator">

  ! Galerkin transport operator (only diagonal part)
  INTEGER, PARAMETER, PUBLIC :: GFSYS_GALERKINDIAG = 0

  ! Galerkin transport operator
  INTEGER, PARAMETER, PUBLIC :: GFSYS_GALERKIN     = 1

  ! Upwind by scalar dissipation (only diagonal part)
  INTEGER, PARAMETER, PUBLIC :: GFSYS_SCALARDIAG   = 2

  ! Upwind by scalar dissipation
  INTEGER, PARAMETER, PUBLIC :: GFSYS_SCALAR       = 3

  ! Upwind by tensorial dissipation (only diagonal part)
  INTEGER, PARAMETER, PUBLIC :: GFSYS_TENSORDIAG   = 4

  ! Upwind by tensorial dissipation
  INTEGER, PARAMETER, PUBLIC :: GFSYS_TENSOR       = 5

!</constantblock>
!</constants>

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

    CASE (AFCSTAB_FEMFCT, AFCSTAB_FEMGP, AFCSTAB_FEMTVD)

      ! Nodal vectors
      ALLOCATE(rafcstab%RnodalVectors(6))
      DO i = 1, 6
        CALL lsyssc_createVector(rafcstab%RnodalVectors(i),&
            rafcstab%NEQ, .FALSE., ST_DOUBLE)
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

      
    CASE (AFCSTAB_FEMFCT, AFCSTAB_FEMGP, AFCSTAB_FEMTVD)

      ! Nodal vectors
      ALLOCATE(rafcstab%RnodalVectors(6))
      DO i = 1, 6
        CALL lsyssc_createVector(rafcstab%RnodalVectors(i),&
            rafcstab%NEQ*rafcstab%NVAR, .FALSE., ST_DOUBLE)
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
    LOGICAL, INTENT(OUT), OPTIONAL :: bcompatible
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
    IF (rafcstab%NEQ   .NE. rmatrixBlock%RmatrixBlock(1,1)%NEQ  .OR.&
        rafcstab%NVAR  .NE. rmatrixBlock%RmatrixBlock(1,1)%NVAR .OR.&
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
    LOGICAL, INTENT(OUT), OPTIONAL :: bcompatible
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
    LOGICAL, INTENT(OUT), OPTIONAL :: bcompatible
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
    LOGICAL, INTENT(OUT), OPTIONAL :: bcompatible
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
      fcb_getRoeMatrix, fcb_getDissipation, iassembly, bclear, rmatrixL)

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

    ! type of assembly
    !   iassembly = 0 : diagonal + Galerkin
    !   iassembly = 1 :            Galerkin
    !   iassembly = 2 : diagonal + scalar dissipation
    !   iassembly = 3 :            scalar dissipation
    !   iassembly = 4 : diagonal + tensorial dissipation
    !   iassembly = 5 :            tensorial dissipation
    INTEGER, INTENT(IN)                            :: iassembly

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
    INTEGER :: idim,ndim,ivar,jvar
    LOGICAL :: bisFullMatrix

    ! Check if block vector contains only one block and if
    ! global operator is stored in interleave format.
    IF (ru%nblocks .EQ. 1   .AND.&
        rmatrixL%ndiagblocks .EQ. 1) THEN
      CALL gfsys_buildDivOperatorScalar(RmatrixC, ru%RvectorBlock(1),&
          fcb_getRoeMatrix, fcb_getDissipation, iassembly,&
          bclear, rmatrixL%RmatrixBlock(1,1))
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
      
      ! What type of assembly should be performed
      SELECT CASE(iassembly)
        
      CASE (GFSYS_GALERKINDIAG) ! diagonal + Galerkin
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doGalerkinMat7Diag_1D(p_Kld, p_Kcol, p_Ksep,&
              rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, rarray)
        CASE (NDIM2D)
          CALL doGalerkinMat7Diag_2D(p_Kld, p_Kcol, p_Ksep,&
              rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, rarray)
        CASE (NDIM3D)
          CALL doGalerkinMat7Diag_3D(p_Kld, p_Kcol, p_Ksep,&
              rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, rarray)
        END SELECT
        

       CASE (GFSYS_GALERKIN)     ! Galerkin
         IF (.NOT.bisFullMatrix) THEN
           CALL output_line('Block matrix has some empty blocks!',&
               OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
           CALL sys_halt()
         END IF
         SELECT CASE(ndim)
         CASE (NDIM1D)
           CALL doGalerkinMat7_1D(p_Kld, p_Kcol, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, rarray)
         CASE (NDIM2D)
           CALL doGalerkinMat7_2D(p_Kld, p_Kcol, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, rarray)
         CASE (NDIM3D)
           CALL doGalerkinMat7_3D(p_Kld, p_Kcol, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, rarray)
         END SELECT
          

       CASE (GFSYS_SCALARDIAG)   ! diagonal + scalar dissipation
         SELECT CASE(ndim)
         CASE (NDIM1D)
           CALL doScalarDissMat7Diag_1D(p_Kld, p_Kcol, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, rarray)
         CASE (NDIM2D)
           CALL doScalarDissMat7Diag_2D(p_Kld, p_Kcol, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, rarray)
         CASE (NDIM3D)
           CALL doScalarDissMat7Diag_3D(p_Kld, p_Kcol, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, rarray)
         END SELECT
         
         
       CASE (GFSYS_SCALAR)       ! scalar dissipation
         IF (.NOT.bisFullMatrix) THEN
           CALL output_line('Block matrix has some empty blocks!',&
               OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
           CALL sys_halt()
         END IF
         SELECT CASE(ndim)
         CASE (NDIM1D)
           CALL doScalarDissMat7_1D(p_Kld, p_Kcol, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, rarray)
         CASE (NDIM2D)
           CALL doScalarDissMat7_2D(p_Kld, p_Kcol, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, rarray)
         CASE (NDIM3D)
           CALL doScalarDissMat7_3D(p_Kld, p_Kcol, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, rarray)
         END SELECT
           
          
       CASE (GFSYS_TENSORDIAG)   ! diagonal + tensorial dissipation
         SELECT CASE(ndim)
         CASE (NDIM1D)
           CALL doTensorDissMat7Diag_1D(p_Kld, p_Kcol, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, rarray)
         CASE (NDIM2D)
           CALL doTensorDissMat7Diag_2D(p_Kld, p_Kcol, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, rarray)
         CASE (NDIM3D)
           CALL doTensorDissMat7Diag_3D(p_Kld, p_Kcol, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, rarray)
         END SELECT
         
         
       CASE (GFSYS_TENSOR)       ! tensorial dissipation
         IF (.NOT.bisFullMatrix) THEN
           CALL output_line('Block matrix has some empty blocks!',&
               OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
           CALL sys_halt()
         END IF
         SELECT CASE(ndim)
         CASE (NDIM1D)
           CALL doTensorDissMat7_1D(p_Kld, p_Kcol, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, rarray)
         CASE (NDIM2D)
           CALL doTensorDissMat7_2D(p_Kld, p_Kcol, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, rarray)
         CASE (NDIM3D)
           CALL doTensorDissMat7_3D(p_Kld, p_Kcol, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, rarray)
         END SELECT

          
       CASE DEFAULT
          CALL output_line('Invalid type of assembly!',&
               OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
          CALL sys_halt()
      END SELECT
      
      ! Release diagonal separator
      CALL storage_free(h_Ksep)
      
      
    CASE(LSYSSC_MATRIX9)

      ! Set pointers
      CALL lsyssc_getbase_Kdiagonal(rmatrixC(1), p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      CALL storage_copy(RmatrixC(1)%h_Kld, h_Ksep)
      CALL storage_getbase_int(h_Ksep, p_Ksep, RmatrixC(1)%NEQ+1)

      ! What type of assembly should be performed
      SELECT CASE(iassembly)
        
      CASE (GFSYS_GALERKINDIAG) ! diagonal + Galerkin
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doGalerkinMat9Diag_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, rarray)
        CASE (NDIM2D)
          CALL doGalerkinMat9Diag_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, rarray)
        CASE (NDIM3D)
          CALL doGalerkinMat9Diag_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, rarray)
        END SELECT
        

       CASE (GFSYS_GALERKIN)     ! Galerkin
         IF (.NOT.bisFullMatrix) THEN
           CALL output_line('Block matrix has some empty blocks!',&
               OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
           CALL sys_halt()
         END IF
         SELECT CASE(ndim)
         CASE (NDIM1D)
           CALL doGalerkinMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, rarray)
         CASE (NDIM2D)
           CALL doGalerkinMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, rarray)
         CASE (NDIM3D)
           CALL doGalerkinMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, rarray)
         END SELECT
          

       CASE (GFSYS_SCALARDIAG)   ! diagonal + scalar dissipation
         SELECT CASE(ndim)
         CASE (NDIM1D)
           CALL doScalarDissMat9Diag_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, rarray)
         CASE (NDIM2D)
           CALL doScalarDissMat9Diag_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, rarray)
         CASE (NDIM3D)
           CALL doScalarDissMat9Diag_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, rarray)
         END SELECT
         
         
       CASE (GFSYS_SCALAR)       ! scalar dissipation
         IF (.NOT.bisFullMatrix) THEN
           CALL output_line('Block matrix has some empty blocks!',&
               OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
           CALL sys_halt()
         END IF
         SELECT CASE(ndim)
         CASE (NDIM1D)
           CALL doScalarDissMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, rarray)
         CASE (NDIM2D)
           CALL doScalarDissMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, rarray)
         CASE (NDIM3D)
           CALL doScalarDissMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, rarray)
         END SELECT
           
          
       CASE (GFSYS_TENSORDIAG)   ! diagonal + tensorial dissipation
         SELECT CASE(ndim)
         CASE (NDIM1D)
           CALL doTensorDissMat9Diag_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, rarray)
         CASE (NDIM2D)
           CALL doTensorDissMat9Diag_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, rarray)
         CASE (NDIM3D)
           CALL doTensorDissMat9Diag_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, rarray)
         END SELECT
         
         
       CASE (GFSYS_TENSOR)       ! tensorial dissipation
         IF (.NOT.bisFullMatrix) THEN
           CALL output_line('Block matrix has some empty blocks!',&
               OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
           CALL sys_halt()
         END IF
         SELECT CASE(ndim)
         CASE (NDIM1D)
           CALL doTensorDissMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, rarray)
         CASE (NDIM2D)
           CALL doTensorDissMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, rarray)
         CASE (NDIM3D)
           CALL doTensorDissMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
               rmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, rarray)
         END SELECT

          
       CASE DEFAULT
          CALL output_line('Invalid type of assembly!',&
               OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
          CALL sys_halt()
      END SELECT
      
      ! Release diagonal separator
      CALL storage_free(h_Ksep)

    CASE DEFAULT
      CALL output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
      CALL sys_halt()
    END SELECT
    
  CONTAINS

    ! Here, the working routines follow
    
    !**************************************************************
    ! Assemble block-diagonal high-order Galerkin operator K in 1D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doGalerkinMat7Diag_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR)   :: A_ij,B_ij,u_i,u_j
      REAL(DP), DIMENSION(NDIM1D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1)=C_ij(1)+Cx(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, B_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii)+A_ij(ivar)+B_ij(ivar)
            K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij)-A_ij(ivar)-B_ij(ivar)
            K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji)+A_ij(ivar)-B_ij(ivar)
            K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj)-A_ij(ivar)+B_ij(ivar)
          END DO
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7Diag_1D

    
    !**************************************************************
    ! Assemble block-diagonal high-order Galerkin operator K in 2D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doGalerkinMat7Diag_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR)   :: A_ij,B_ij,u_i,u_j
      REAL(DP), DIMENSION(NDIM2D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1)=C_ij(1)+Cx(ij)
          C_ij(2)=C_ij(2)+Cy(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, B_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii)+A_ij(ivar)+B_ij(ivar)
            K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij)-A_ij(ivar)-B_ij(ivar)
            K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji)+A_ij(ivar)-B_ij(ivar)
            K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj)-A_ij(ivar)+B_ij(ivar)
          END DO
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7Diag_2D


    !**************************************************************
    ! Assemble block-diagonal high-order Galerkin operator K in 3D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doGalerkinMat7Diag_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR)   :: A_ij,B_ij,u_i,u_j
      REAL(DP), DIMENSION(NDIM3D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1)=C_ij(1)+Cx(ij)
          C_ij(2)=C_ij(2)+Cy(ij)
          C_ij(3)=C_ij(3)+Cz(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, B_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii)+A_ij(ivar)+B_ij(ivar)
            K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij)-A_ij(ivar)-B_ij(ivar)
            K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji)+A_ij(ivar)-B_ij(ivar)
            K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj)-A_ij(ivar)+B_ij(ivar)
          END DO
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7Diag_3D


    !**************************************************************
    ! Assemble block-diagonal high-order Galerkin operator K in 1D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doGalerkinMat9Diag_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR)   :: A_ij,B_ij,u_i,u_j
      REAL(DP), DIMENSION(NDIM1D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1) = C_ij(1)+Cx(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, B_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii)+A_ij(ivar)+B_ij(ivar)
            K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij)-A_ij(ivar)-B_ij(ivar)
            K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji)+A_ij(ivar)-B_ij(ivar)
            K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj)-A_ij(ivar)+B_ij(ivar)
          END DO
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9Diag_1D


    !**************************************************************
    ! Assemble block-diagonal high-order Galerkin operator K in 2D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doGalerkinMat9Diag_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR)   :: A_ij,B_ij,u_i,u_j
      REAL(DP), DIMENSION(NDIM2D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1) = C_ij(1)+Cx(ij)
          C_ij(2) = C_ij(2)+Cy(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, B_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii)+A_ij(ivar)+B_ij(ivar)
            K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij)-A_ij(ivar)-B_ij(ivar)
            K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji)+A_ij(ivar)-B_ij(ivar)
            K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj)-A_ij(ivar)+B_ij(ivar)
          END DO
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9Diag_2D


    !**************************************************************
    ! Assemble block-diagonal high-order Galerkin operator K in 2D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doGalerkinMat9Diag_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR)   :: A_ij,B_ij,u_i,u_j
      REAL(DP), DIMENSION(NDIM3D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1) = C_ij(1)+Cx(ij)
          C_ij(2) = C_ij(2)+Cy(ij)
          C_ij(3) = C_ij(3)+Cz(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, B_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii)+A_ij(ivar)+B_ij(ivar)
            K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij)-A_ij(ivar)-B_ij(ivar)
            K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji)+A_ij(ivar)-B_ij(ivar)
            K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj)-A_ij(ivar)+B_ij(ivar)
          END DO
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9Diag_3D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 1D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doGalerkinMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute averaged coefficients for boundary
          C_ij(1) = C_ij(1)+Cx(ij)
          
          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, B_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii)+A_ij(idx)+B_ij(idx)
              K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij)-A_ij(idx)-B_ij(idx)
              K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji)+A_ij(idx)-B_ij(idx)
              K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj)-A_ij(idx)+B_ij(idx)
            END DO
          END DO
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7_1D

    
    !**************************************************************
    ! Assemble high-order Galerkin operator K in 2D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doGalerkinMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute averaged coefficients for boundary
          C_ij(1) = C_ij(1)+Cx(ij)
          C_ij(2) = C_ij(2)+Cy(ij)
          
          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, B_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii)+A_ij(idx)+B_ij(idx)
              K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij)-A_ij(idx)-B_ij(idx)
              K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji)+A_ij(idx)-B_ij(idx)
              K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj)-A_ij(idx)+B_ij(idx)
            END DO
          END DO
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7_2D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 3D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doGalerkinMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute averaged coefficients for boundary
          C_ij(1) = C_ij(1)+Cx(ij)
          C_ij(2) = C_ij(2)+Cy(ij)
          C_ij(3) = C_ij(3)+Cz(ij)
          
          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, B_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii)+A_ij(idx)+B_ij(idx)
              K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij)-A_ij(idx)-B_ij(idx)
              K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji)+A_ij(idx)-B_ij(idx)
              K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj)-A_ij(idx)+B_ij(idx)
            END DO
          END DO
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7_3D

    
    !**************************************************************
    ! Assemble high-order Galerkin operator K in 1D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doGalerkinMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          
          ! Compute averaged coefficients for boundary
          C_ij(1) = C_ij(1)+Cx(ij)
          
          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, B_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii)+A_ij(idx)+B_ij(idx)
              K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij)-A_ij(idx)-B_ij(idx)
              K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji)+A_ij(idx)-B_ij(idx)
              K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj)-A_ij(idx)+B_ij(idx)
            END DO
          END DO
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9_1D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 2D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doGalerkinMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          
          ! Compute averaged coefficients for boundary
          C_ij(1) = C_ij(1)+Cx(ij)
          C_ij(2) = C_ij(2)+Cy(ij)
          
          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, B_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii)+A_ij(idx)+B_ij(idx)
              K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij)-A_ij(idx)-B_ij(idx)
              K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji)+A_ij(idx)-B_ij(idx)
              K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj)-A_ij(idx)+B_ij(idx)
            END DO
          END DO
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9_2D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 3D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doGalerkinMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          
          ! Compute averaged coefficients for boundary
          C_ij(1) = C_ij(1)+Cx(ij)
          C_ij(2) = C_ij(2)+Cy(ij)
          C_ij(3) = C_ij(3)+Cz(ij)
          
          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, B_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii)+A_ij(idx)+B_ij(idx)
              K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij)-A_ij(idx)-B_ij(idx)
              K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji)+A_ij(idx)-B_ij(idx)
              K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj)-A_ij(idx)+B_ij(idx)
            END DO
          END DO
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9_3D

    
    !**************************************************************
    ! Assemble block-diagonal low-order operator L 
    ! with scalar dissipation in 1D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doScalarDissMat7Diag_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,u_i,u_j
      REAL(DP), DIMENSION(1)      :: d_ij
      REAL(DP), DIMENSION(NDIM1D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          
          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, d_ij)

          ! Assemble the global operator
          DO ivar = 1, NVAR
            L(ivar,ivar)%DA(ii) = L(ivar,ivar)%DA(ii)+A_ij(ivar)-d_ij(1)
            L(ivar,ivar)%DA(ij) = L(ivar,ivar)%DA(ij)-A_ij(ivar)+d_ij(1)
            L(ivar,ivar)%DA(ji) = L(ivar,ivar)%DA(ji)+A_ij(ivar)+d_ij(1)
            L(ivar,ivar)%DA(jj) = L(ivar,ivar)%DA(jj)-A_ij(ivar)-d_ij(1)
          END DO
        END DO
      END DO
    END SUBROUTINE doScalarDissMat7Diag_1D


    !**************************************************************
    ! Assemble block-diagonal low-order operator L 
    ! with scalar dissipation in 2D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doScalarDissMat7Diag_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, Cy, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,u_i,u_j
      REAL(DP), DIMENSION(1)      :: d_ij
      REAL(DP), DIMENSION(NDIM2D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          
          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, d_ij)

          ! Assemble the global operator
          DO ivar = 1, NVAR
            L(ivar,ivar)%DA(ii) = L(ivar,ivar)%DA(ii)+A_ij(ivar)-d_ij(1)
            L(ivar,ivar)%DA(ij) = L(ivar,ivar)%DA(ij)-A_ij(ivar)+d_ij(1)
            L(ivar,ivar)%DA(ji) = L(ivar,ivar)%DA(ji)+A_ij(ivar)+d_ij(1)
            L(ivar,ivar)%DA(jj) = L(ivar,ivar)%DA(jj)-A_ij(ivar)-d_ij(1)
          END DO
        END DO
      END DO
    END SUBROUTINE doScalarDissMat7Diag_2D

    
    !**************************************************************
    ! Assemble block-diagonal low-order operator L 
    ! with scalar dissipation in 3D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doScalarDissMat7Diag_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, Cy, Cz, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,u_i,u_j
      REAL(DP), DIMENSION(1)      :: d_ij
      REAL(DP), DIMENSION(NDIM3D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          
          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, d_ij)

          ! Assemble the global operator
          DO ivar = 1, NVAR
            L(ivar,ivar)%DA(ii) = L(ivar,ivar)%DA(ii)+A_ij(ivar)-d_ij(1)
            L(ivar,ivar)%DA(ij) = L(ivar,ivar)%DA(ij)-A_ij(ivar)+d_ij(1)
            L(ivar,ivar)%DA(ji) = L(ivar,ivar)%DA(ji)+A_ij(ivar)+d_ij(1)
            L(ivar,ivar)%DA(jj) = L(ivar,ivar)%DA(jj)-A_ij(ivar)-d_ij(1)
          END DO
        END DO
      END DO
    END SUBROUTINE doScalarDissMat7Diag_3D


    !**************************************************************
    ! Assemble block-diagonal low-order operator L
    ! with scalar dissipation in 1D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doScalarDissMat9Diag_1D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,u_i,u_j
      REAL(DP), DIMENSION(1)      :: d_ij
      REAL(DP), DIMENSION(NDIM1D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          
          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, d_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            L(ivar,ivar)%DA(ii) = L(ivar,ivar)%DA(ii)+A_ij(ivar)-d_ij(1)
            L(ivar,ivar)%DA(ij) = L(ivar,ivar)%DA(ij)-A_ij(ivar)+d_ij(1)
            L(ivar,ivar)%DA(ji) = L(ivar,ivar)%DA(ji)+A_ij(ivar)+d_ij(1)
            L(ivar,ivar)%DA(jj) = L(ivar,ivar)%DA(jj)-A_ij(ivar)-d_ij(1)
          END DO
        END DO
      END DO
    END SUBROUTINE doScalarDissMat9Diag_1D


    !**************************************************************
    ! Assemble block-diagonal low-order operator L
    ! with scalar dissipation in 2D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doScalarDissMat9Diag_2D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, Cy, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,u_i,u_j
      REAL(DP), DIMENSION(1)      :: d_ij
      REAL(DP), DIMENSION(NDIM2D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          
          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, d_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            L(ivar,ivar)%DA(ii) = L(ivar,ivar)%DA(ii)+A_ij(ivar)-d_ij(1)
            L(ivar,ivar)%DA(ij) = L(ivar,ivar)%DA(ij)-A_ij(ivar)+d_ij(1)
            L(ivar,ivar)%DA(ji) = L(ivar,ivar)%DA(ji)+A_ij(ivar)+d_ij(1)
            L(ivar,ivar)%DA(jj) = L(ivar,ivar)%DA(jj)-A_ij(ivar)-d_ij(1)
          END DO
        END DO
      END DO
    END SUBROUTINE doScalarDissMat9Diag_2D


    !**************************************************************
    ! Assemble block-diagonal low-order operator L
    ! with scalar dissipation in 3D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doScalarDissMat9Diag_3D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, Cy, Cz, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,u_i,u_j
      REAL(DP), DIMENSION(1)      :: d_ij
      REAL(DP), DIMENSION(NDIM3D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          
          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, d_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            L(ivar,ivar)%DA(ii) = L(ivar,ivar)%DA(ii)+A_ij(ivar)-d_ij(1)
            L(ivar,ivar)%DA(ij) = L(ivar,ivar)%DA(ij)-A_ij(ivar)+d_ij(1)
            L(ivar,ivar)%DA(ji) = L(ivar,ivar)%DA(ji)+A_ij(ivar)+d_ij(1)
            L(ivar,ivar)%DA(jj) = L(ivar,ivar)%DA(jj)-A_ij(ivar)-d_ij(1)
          END DO
        END DO
      END DO
    END SUBROUTINE doScalarDissMat9Diag_3D

    
    !**************************************************************
    ! Assemble low-order operator L with scalar dissipation in 1D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doScalarDissMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, d_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              L(jvar,ivar)%DA(ii) = L(jvar,ivar)%DA(ii)+A_ij(idx)
              L(jvar,ivar)%DA(ij) = L(jvar,ivar)%DA(ij)-A_ij(idx)
              L(jvar,ivar)%DA(ji) = L(jvar,ivar)%DA(ji)+A_ij(idx)
              L(jvar,ivar)%DA(jj) = L(jvar,ivar)%DA(jj)-A_ij(idx)
            END DO
          END DO
          
          ! Apply the scalar dissipation I*d_ij
          DO ivar = 1, NVAR
            L(ivar,ivar)%DA(ii) = L(ivar,ivar)%DA(ii)-d_ij(1)
            L(ivar,ivar)%DA(ij) = L(ivar,ivar)%DA(ij)+d_ij(1)
            L(ivar,ivar)%DA(ji) = L(ivar,ivar)%DA(ji)+d_ij(1)
            L(ivar,ivar)%DA(jj) = L(ivar,ivar)%DA(jj)-d_ij(1)
          END DO
        END DO
      END DO
    END SUBROUTINE doScalarDissMat7_1D


    !**************************************************************
    ! Assemble low-order operator L with scalar dissipation in 2D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doScalarDissMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, Cy, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, d_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              L(jvar,ivar)%DA(ii) = L(jvar,ivar)%DA(ii)+A_ij(idx)
              L(jvar,ivar)%DA(ij) = L(jvar,ivar)%DA(ij)-A_ij(idx)
              L(jvar,ivar)%DA(ji) = L(jvar,ivar)%DA(ji)+A_ij(idx)
              L(jvar,ivar)%DA(jj) = L(jvar,ivar)%DA(jj)-A_ij(idx)
            END DO
          END DO
          
          ! Apply the scalar dissipation I*d_ij
          DO ivar = 1, NVAR
            L(ivar,ivar)%DA(ii) = L(ivar,ivar)%DA(ii)-d_ij(1)
            L(ivar,ivar)%DA(ij) = L(ivar,ivar)%DA(ij)+d_ij(1)
            L(ivar,ivar)%DA(ji) = L(ivar,ivar)%DA(ji)+d_ij(1)
            L(ivar,ivar)%DA(jj) = L(ivar,ivar)%DA(jj)-d_ij(1)
          END DO
        END DO
      END DO
    END SUBROUTINE doScalarDissMat7_2D


    !**************************************************************
    ! Assemble low-order operator L with scalar dissipation in 3D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doScalarDissMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, Cy, Cz, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, d_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              L(jvar,ivar)%DA(ii) = L(jvar,ivar)%DA(ii)+A_ij(idx)
              L(jvar,ivar)%DA(ij) = L(jvar,ivar)%DA(ij)-A_ij(idx)
              L(jvar,ivar)%DA(ji) = L(jvar,ivar)%DA(ji)+A_ij(idx)
              L(jvar,ivar)%DA(jj) = L(jvar,ivar)%DA(jj)-A_ij(idx)
            END DO
          END DO
          
          ! Apply the scalar dissipation I*d_ij
          DO ivar = 1, NVAR
            L(ivar,ivar)%DA(ii) = L(ivar,ivar)%DA(ii)-d_ij(1)
            L(ivar,ivar)%DA(ij) = L(ivar,ivar)%DA(ij)+d_ij(1)
            L(ivar,ivar)%DA(ji) = L(ivar,ivar)%DA(ji)+d_ij(1)
            L(ivar,ivar)%DA(jj) = L(ivar,ivar)%DA(jj)-d_ij(1)
          END DO
        END DO
      END DO
    END SUBROUTINE doScalarDissMat7_3D

    
    !**************************************************************
    ! Assemble low-order operator L with scalar dissipation in 1D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doScalarDissMat9_1D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, d_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              L(jvar,ivar)%DA(ii) = L(jvar,ivar)%DA(ii)+A_ij(idx)
              L(jvar,ivar)%DA(ij) = L(jvar,ivar)%DA(ij)-A_ij(idx)
              L(jvar,ivar)%DA(ji) = L(jvar,ivar)%DA(ji)+A_ij(idx)
              L(jvar,ivar)%DA(jj) = L(jvar,ivar)%DA(jj)-A_ij(idx)
            END DO
          END DO
          
          ! Apply the scalar dissipation I*d_ij
          DO ivar = 1, NVAR
            L(ivar,ivar)%DA(ii) = L(ivar,ivar)%DA(ii)-d_ij(1)
            L(ivar,ivar)%DA(ij) = L(ivar,ivar)%DA(ij)+d_ij(1)
            L(ivar,ivar)%DA(ji) = L(ivar,ivar)%DA(ji)+d_ij(1)
            L(ivar,ivar)%DA(jj) = L(ivar,ivar)%DA(jj)-d_ij(1)
          END DO
        END DO
      END DO
    END SUBROUTINE doScalarDissMat9_1D


    !**************************************************************
    ! Assemble low-order operator L with scalar dissipation in 2D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doScalarDissMat9_2D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, Cy, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, d_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              L(jvar,ivar)%DA(ii) = L(jvar,ivar)%DA(ii)+A_ij(idx)
              L(jvar,ivar)%DA(ij) = L(jvar,ivar)%DA(ij)-A_ij(idx)
              L(jvar,ivar)%DA(ji) = L(jvar,ivar)%DA(ji)+A_ij(idx)
              L(jvar,ivar)%DA(jj) = L(jvar,ivar)%DA(jj)-A_ij(idx)
            END DO
          END DO
          
          ! Apply the scalar dissipation I*d_ij
          DO ivar = 1, NVAR
            L(ivar,ivar)%DA(ii) = L(ivar,ivar)%DA(ii)-d_ij(1)
            L(ivar,ivar)%DA(ij) = L(ivar,ivar)%DA(ij)+d_ij(1)
            L(ivar,ivar)%DA(ji) = L(ivar,ivar)%DA(ji)+d_ij(1)
            L(ivar,ivar)%DA(jj) = L(ivar,ivar)%DA(jj)-d_ij(1)
          END DO
        END DO
      END DO
    END SUBROUTINE doScalarDissMat9_2D


    !**************************************************************
    ! Assemble low-order operator L with scalar dissipation in 3D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doScalarDissMat9_3D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, Cy, Cz, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, d_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              L(jvar,ivar)%DA(ii) = L(jvar,ivar)%DA(ii)+A_ij(idx)
              L(jvar,ivar)%DA(ij) = L(jvar,ivar)%DA(ij)-A_ij(idx)
              L(jvar,ivar)%DA(ji) = L(jvar,ivar)%DA(ji)+A_ij(idx)
              L(jvar,ivar)%DA(jj) = L(jvar,ivar)%DA(jj)-A_ij(idx)
            END DO
          END DO
          
          ! Apply the scalar dissipation I*d_ij
          DO ivar = 1, NVAR
            L(ivar,ivar)%DA(ii) = L(ivar,ivar)%DA(ii)-d_ij(1)
            L(ivar,ivar)%DA(ij) = L(ivar,ivar)%DA(ij)+d_ij(1)
            L(ivar,ivar)%DA(ji) = L(ivar,ivar)%DA(ji)+d_ij(1)
            L(ivar,ivar)%DA(jj) = L(ivar,ivar)%DA(jj)-d_ij(1)
          END DO
        END DO
      END DO
    END SUBROUTINE doScalarDissMat9_3D


    !**************************************************************
    ! Assemble block-diagonal low-order operator L 
    ! with tensorial dissipation in 1D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doTensorDissMat7Diag_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij
      REAL(DP), DIMENSION(NVAR)   :: D_ij
      REAL(DP), DIMENSION(NVAR)   :: u_i,u_j
      REAL(DP), DIMENSION(NDIM1D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, D_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            L(ivar,ivar)%DA(ii) = L(ivar,ivar)%DA(ii)+A_ij(ivar)-D_ij(ivar)
            L(ivar,ivar)%DA(ij) = L(ivar,ivar)%DA(ij)-A_ij(ivar)+D_ij(ivar)
            L(ivar,ivar)%DA(ji) = L(ivar,ivar)%DA(ji)+A_ij(ivar)+D_ij(ivar)
            L(ivar,ivar)%DA(jj) = L(ivar,ivar)%DA(jj)-A_ij(ivar)-D_ij(ivar)
          END DO
        END DO
      END DO
    END SUBROUTINE doTensorDissMat7Diag_1D


    !**************************************************************
    ! Assemble block-diagonal low-order operator L 
    ! with tensorial dissipation in 2D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doTensorDissMat7Diag_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, Cy, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij
      REAL(DP), DIMENSION(NVAR)   :: D_ij
      REAL(DP), DIMENSION(NVAR)   :: u_i,u_j
      REAL(DP), DIMENSION(NDIM2D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, D_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            L(ivar,ivar)%DA(ii) = L(ivar,ivar)%DA(ii)+A_ij(ivar)-D_ij(ivar)
            L(ivar,ivar)%DA(ij) = L(ivar,ivar)%DA(ij)-A_ij(ivar)+D_ij(ivar)
            L(ivar,ivar)%DA(ji) = L(ivar,ivar)%DA(ji)+A_ij(ivar)+D_ij(ivar)
            L(ivar,ivar)%DA(jj) = L(ivar,ivar)%DA(jj)-A_ij(ivar)-D_ij(ivar)
          END DO
        END DO
      END DO
    END SUBROUTINE doTensorDissMat7Diag_2D


    !**************************************************************
    ! Assemble block-diagonal low-order operator L 
    ! with tensorial dissipation in 3D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doTensorDissMat7Diag_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, Cy, Cz, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij
      REAL(DP), DIMENSION(NVAR)   :: D_ij
      REAL(DP), DIMENSION(NVAR)   :: u_i,u_j
      REAL(DP), DIMENSION(NDIM3D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, D_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            L(ivar,ivar)%DA(ii) = L(ivar,ivar)%DA(ii)+A_ij(ivar)-D_ij(ivar)
            L(ivar,ivar)%DA(ij) = L(ivar,ivar)%DA(ij)-A_ij(ivar)+D_ij(ivar)
            L(ivar,ivar)%DA(ji) = L(ivar,ivar)%DA(ji)+A_ij(ivar)+D_ij(ivar)
            L(ivar,ivar)%DA(jj) = L(ivar,ivar)%DA(jj)-A_ij(ivar)-D_ij(ivar)
          END DO
        END DO
      END DO
    END SUBROUTINE doTensorDissMat7Diag_3D


    !**************************************************************
    ! Assemble block-diagonal low-order operator L 
    ! with tensorial dissipation in 1D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doTensorDissMat9Diag_1D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij
      REAL(DP), DIMENSION(NVAR)   :: D_ij
      REAL(DP), DIMENSION(NVAR)   :: u_i,u_j
      REAL(DP), DIMENSION(NDIM1D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, D_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            L(ivar,ivar)%DA(ii) = L(ivar,ivar)%DA(ii)+A_ij(ivar)-D_ij(ivar)
            L(ivar,ivar)%DA(ij) = L(ivar,ivar)%DA(ij)-A_ij(ivar)+D_ij(ivar)
            L(ivar,ivar)%DA(ji) = L(ivar,ivar)%DA(ji)+A_ij(ivar)+D_ij(ivar)
            L(ivar,ivar)%DA(jj) = L(ivar,ivar)%DA(jj)-A_ij(ivar)-D_ij(ivar)
          END DO
        END DO
      END DO
    END SUBROUTINE doTensorDissMat9Diag_1D


    !**************************************************************
    ! Assemble block-diagonal low-order operator L 
    ! with tensorial dissipation in 2D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doTensorDissMat9Diag_2D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, Cy, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij
      REAL(DP), DIMENSION(NVAR)   :: D_ij
      REAL(DP), DIMENSION(NVAR)   :: u_i,u_j
      REAL(DP), DIMENSION(NDIM2D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, D_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            L(ivar,ivar)%DA(ii) = L(ivar,ivar)%DA(ii)+A_ij(ivar)-D_ij(ivar)
            L(ivar,ivar)%DA(ij) = L(ivar,ivar)%DA(ij)-A_ij(ivar)+D_ij(ivar)
            L(ivar,ivar)%DA(ji) = L(ivar,ivar)%DA(ji)+A_ij(ivar)+D_ij(ivar)
            L(ivar,ivar)%DA(jj) = L(ivar,ivar)%DA(jj)-A_ij(ivar)-D_ij(ivar)
          END DO
        END DO
      END DO
    END SUBROUTINE doTensorDissMat9Diag_2D


     !**************************************************************
    ! Assemble block-diagonal low-order operator L 
    ! with tensorial dissipation in 3D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doTensorDissMat9Diag_3D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, Cy, Cz, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij
      REAL(DP), DIMENSION(NVAR)   :: D_ij
      REAL(DP), DIMENSION(NVAR)   :: u_i,u_j
      REAL(DP), DIMENSION(NDIM3D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, D_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            L(ivar,ivar)%DA(ii) = L(ivar,ivar)%DA(ii)+A_ij(ivar)-D_ij(ivar)
            L(ivar,ivar)%DA(ij) = L(ivar,ivar)%DA(ij)-A_ij(ivar)+D_ij(ivar)
            L(ivar,ivar)%DA(ji) = L(ivar,ivar)%DA(ji)+A_ij(ivar)+D_ij(ivar)
            L(ivar,ivar)%DA(jj) = L(ivar,ivar)%DA(jj)-A_ij(ivar)-D_ij(ivar)
          END DO
        END DO
      END DO
    END SUBROUTINE doTensorDissMat9Diag_3D


    !**************************************************************
    ! Assemble low-order operator L with tensorial dissipation in 1D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doTensorDissMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          
          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, D_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              L(jvar,ivar)%DA(ii) = L(jvar,ivar)%DA(ii)+A_ij(idx)-D_ij(idx)
              L(jvar,ivar)%DA(ij) = L(jvar,ivar)%DA(ij)-A_ij(idx)+D_ij(idx)
              L(jvar,ivar)%DA(ji) = L(jvar,ivar)%DA(ji)+A_ij(idx)+D_ij(idx)
              L(jvar,ivar)%DA(jj) = L(jvar,ivar)%DA(jj)-A_ij(idx)-D_ij(idx)
            END DO
          END DO
        END DO
      END DO
    END SUBROUTINE doTensorDissMat7_1D


    !**************************************************************
    ! Assemble low-order operator L with tensorial dissipation in 2D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doTensorDissMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, Cy, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          
          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, D_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              L(jvar,ivar)%DA(ii) = L(jvar,ivar)%DA(ii)+A_ij(idx)-D_ij(idx)
              L(jvar,ivar)%DA(ij) = L(jvar,ivar)%DA(ij)-A_ij(idx)+D_ij(idx)
              L(jvar,ivar)%DA(ji) = L(jvar,ivar)%DA(ji)+A_ij(idx)+D_ij(idx)
              L(jvar,ivar)%DA(jj) = L(jvar,ivar)%DA(jj)-A_ij(idx)-D_ij(idx)
            END DO
          END DO
        END DO
      END DO
    END SUBROUTINE doTensorDissMat7_2D


    !**************************************************************
    ! Assemble low-order operator L with tensorial dissipation in 3D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doTensorDissMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
        Cx, Cy, Cz, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          
          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, D_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              L(jvar,ivar)%DA(ii) = L(jvar,ivar)%DA(ii)+A_ij(idx)-D_ij(idx)
              L(jvar,ivar)%DA(ij) = L(jvar,ivar)%DA(ij)-A_ij(idx)+D_ij(idx)
              L(jvar,ivar)%DA(ji) = L(jvar,ivar)%DA(ji)+A_ij(idx)+D_ij(idx)
              L(jvar,ivar)%DA(jj) = L(jvar,ivar)%DA(jj)-A_ij(idx)-D_ij(idx)
            END DO
          END DO
        END DO
      END DO
    END SUBROUTINE doTensorDissMat7_3D


    !**************************************************************
    ! Assemble low-order operator L with tensorial dissipation in 1D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doTensorDissMat9_1D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, D_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              L(jvar,ivar)%DA(ii) = L(jvar,ivar)%DA(ii)+A_ij(idx)-D_ij(idx)
              L(jvar,ivar)%DA(ij) = L(jvar,ivar)%DA(ij)-A_ij(idx)+D_ij(idx)
              L(jvar,ivar)%DA(ji) = L(jvar,ivar)%DA(ji)+A_ij(idx)+D_ij(idx)
              L(jvar,ivar)%DA(jj) = L(jvar,ivar)%DA(jj)-A_ij(idx)-D_ij(idx)
            END DO
          END DO
        END DO
      END DO
    END SUBROUTINE doTensorDissMat9_1D


    !**************************************************************
    ! Assemble low-order operator L with tensorial dissipation in 2D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doTensorDissMat9_2D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, Cy, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)      :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, D_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              L(jvar,ivar)%DA(ii) = L(jvar,ivar)%DA(ii)+A_ij(idx)-D_ij(idx)
              L(jvar,ivar)%DA(ij) = L(jvar,ivar)%DA(ij)-A_ij(idx)+D_ij(idx)
              L(jvar,ivar)%DA(ji) = L(jvar,ivar)%DA(ji)+A_ij(idx)+D_ij(idx)
              L(jvar,ivar)%DA(jj) = L(jvar,ivar)%DA(jj)-A_ij(idx)-D_ij(idx)
            END DO
          END DO
        END DO
      END DO
    END SUBROUTINE doTensorDissMat9_2D


    !**************************************************************
    ! Assemble low-order operator L with tensorial dissipation in 3D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doTensorDissMat9_3D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, Cy, Cz, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      TYPE(t_array), DIMENSION(:,:), INTENT(INOUT)       :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, D_ij)
          
          ! Assemble the global operator
          DO ivar = 1, NVAR
            DO jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              L(jvar,ivar)%DA(ii) = L(jvar,ivar)%DA(ii)+A_ij(idx)-D_ij(idx)
              L(jvar,ivar)%DA(ij) = L(jvar,ivar)%DA(ij)-A_ij(idx)+D_ij(idx)
              L(jvar,ivar)%DA(ji) = L(jvar,ivar)%DA(ji)+A_ij(idx)+D_ij(idx)
              L(jvar,ivar)%DA(jj) = L(jvar,ivar)%DA(jj)-A_ij(idx)-D_ij(idx)
            END DO
          END DO
        END DO
      END DO
    END SUBROUTINE doTensorDissMat9_3D
  END SUBROUTINE gfsys_buildDivOperatorBlock

  ! *****************************************************************************

!<subroutine>

  SUBROUTINE gfsys_buildDivOperatorScalar(RmatrixC, ru,&
      fcb_getRoeMatrix, fcb_getDissipation, iassembly, bclear, rmatrixL)

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

    ! type of assembly
    !   iassembly = 0 : diagonal + Galerkin
    !   iassembly = 1 :            Galerkin
    !   iassembly = 2 : diagonal + scalar dissipation
    !   iassembly = 3 :            scalar dissipation
    !   iassembly = 4 : diagonal + tensorial dissipation
    !   iassembly = 5 :            tensorial dissipation
    INTEGER, INTENT(IN)                            :: iassembly

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
      
      ! What type of assembly should be performed
      SELECT CASE(iassembly)
        
      CASE (GFSYS_GALERKINDIAG) ! diagonal + Galerkin
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doGalerkinMat7Diag_1D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, p_L)
        CASE (NDIM2D)
          CALL doGalerkinMat7Diag_2D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, p_L)
        CASE (NDIM3D)
          CALL doGalerkinMat7Diag_3D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, p_L)
        END SELECT
        
        
      CASE (GFSYS_GALERKIN)     ! Galerkin
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doGalerkinMat7_1D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, p_L)
        CASE (NDIM2D)
          CALL doGalerkinMat7_2D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, p_L)
        CASE (NDIM3D)
          CALL doGalerkinMat7_3D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, p_L)
        END SELECT
        
        
      CASE (GFSYS_SCALARDIAG)   ! diagonal + scalar dissipation
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doScalarDissMat7Diag_1D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, p_L)
        CASE (NDIM2D)
          CALL doScalarDissMat7Diag_2D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, p_L)
        CASE (NDIM3D)
          CALL doScalarDissMat7Diag_3D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, p_L)
        END SELECT
        

      CASE (GFSYS_SCALAR)       ! scalar dissipation
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doScalarDissMat7_1D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, p_L)
        CASE (NDIM2D)
          CALL doScalarDissMat7_2D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, p_L)
        CASE (NDIM3D)
          CALL doScalarDissMat7_3D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, p_L)
        END SELECT
        
        
      CASE (GFSYS_TENSORDIAG)   ! diagonal + tensorial dissipation
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doTensorDissMat7Diag_1D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, p_L)
        CASE (NDIM2D)
          CALL doTensorDissMat7Diag_2D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, p_L)
        CASE (NDIM3D)
          CALL doTensorDissMat7Diag_3D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, p_L)
        END SELECT
        

      CASE (GFSYS_TENSOR)       ! tensorial dissipation
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doTensorDissMat7_1D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, p_L)
        CASE (NDIM2D)
          CALL doTensorDissMat7_2D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, p_L)
        CASE (NDIM3D)
          CALL doTensorDissMat7_3D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, p_L)
        END SELECT


      CASE DEFAULT
        CALL output_line('Invalid type of assembly!',&
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
      
      ! What type of assembly should be performed
      SELECT CASE(iassembly)
        
      CASE (GFSYS_GALERKINDIAG) ! diagonal + Galerkin
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doGalerkinMat9Diag_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, p_L)
        CASE (NDIM2D)
          CALL doGalerkinMat9Diag_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, p_L)
        CASE (NDIM3D)
          CALL doGalerkinMat9Diag_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, p_L)
        END SELECT
        
        
      CASE (GFSYS_GALERKIN)     ! Galerkin
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doGalerkinMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, p_L)
        CASE (NDIM2D)
          CALL doGalerkinMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, p_L)
        CASE (NDIM3D)
          CALL doGalerkinMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, p_L)
        END SELECT
        

      CASE (GFSYS_SCALARDIAG)   ! diagonal + scalar dissipation
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doScalarDissMat9Diag_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, p_L)
        CASE (NDIM2D)
          CALL doScalarDissMat9Diag_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, p_L)
        CASE (NDIM3D)
          CALL doScalarDissMat9Diag_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, p_L)
        END SELECT
        
        
      CASE (GFSYS_SCALAR)       ! scalar dissipation
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doScalarDissMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, p_L)
        CASE (NDIM2D)
          CALL doScalarDissMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, p_L)
        CASE (NDIM3D)
          CALL doScalarDissMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, p_L)
        END SELECT

        
      CASE (GFSYS_TENSORDIAG)   ! diagonal + tensorial dissipation
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doTensorDissMat9Diag_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, p_L)
        CASE (NDIM2D)
          CALL doTensorDissMat9Diag_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, p_L)
        CASE (NDIM3D)
          CALL doTensorDissMat9Diag_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, p_L)
        END SELECT
        
        
      CASE (GFSYS_TENSOR)       ! tensorial dissipation
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL doTensorDissMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, p_L)
        CASE (NDIM2D)
          CALL doTensorDissMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, p_L)
        CASE (NDIM3D)
          CALL doTensorDissMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, p_L)
        END SELECT

        
      CASE DEFAULT
        CALL output_line('Invalid type of assembly!',&
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
    
  CONTAINS
    
    ! Here, the working routines follow
    
    !**************************************************************
    ! Assemble block-diagonal high-order Galerkin operator K in 1D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doGalerkinMat7Diag_1D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,B_ij
      REAL(DP), DIMENSION(NDIM1D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1)=C_ij(1)+Cx(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, B_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+B_ij
          K(:,ij) = K(:,ij)-A_ij-B_ij
          K(:,ji) = K(:,ji)+A_ij-B_ij
          K(:,jj) = K(:,jj)-A_ij+B_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7Diag_1D


    !**************************************************************
    ! Assemble block-diagonal high-order Galerkin operator K in 2D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doGalerkinMat7Diag_2D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, Cy, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,B_ij
      REAL(DP), DIMENSION(NDIM2D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1)=C_ij(1)+Cx(ij)
          C_ij(2)=C_ij(2)+Cy(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, B_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+B_ij
          K(:,ij) = K(:,ij)-A_ij-B_ij
          K(:,ji) = K(:,ji)+A_ij-B_ij
          K(:,jj) = K(:,jj)-A_ij+B_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7Diag_2D


    !**************************************************************
    ! Assemble block-diagonal high-order Galerkin operator K in 3D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doGalerkinMat7Diag_3D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, Cy, Cz, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij,B_ij
      REAL(DP), DIMENSION(NDIM3D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1)=C_ij(1)+Cx(ij)
          C_ij(2)=C_ij(2)+Cy(ij)
          C_ij(3)=C_ij(3)+Cz(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, B_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+B_ij
          K(:,ij) = K(:,ij)-A_ij-B_ij
          K(:,ji) = K(:,ji)+A_ij-B_ij
          K(:,jj) = K(:,jj)-A_ij+B_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7Diag_3D

    
    !**************************************************************
    ! Assemble block-diagonal high-order Galerkin operator K in 1D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doGalerkinMat9Diag_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, Cx, u, K)
      
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
      REAL(DP), DIMENSION(NVAR)   :: A_ij,B_ij
      REAL(DP), DIMENSION(NDIM1D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1) = C_ij(1)+Cx(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, B_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+B_ij
          K(:,ij) = K(:,ij)-A_ij-B_ij
          K(:,ji) = K(:,ji)+A_ij-B_ij
          K(:,jj) = K(:,jj)-A_ij+B_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9Diag_1D


    !**************************************************************
    ! Assemble block-diagonal high-order Galerkin operator K in 2D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doGalerkinMat9Diag_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, Cx, Cy, u, K)
      
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
      REAL(DP), DIMENSION(NVAR)   :: A_ij,B_ij
      REAL(DP), DIMENSION(NDIM2D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1) = C_ij(1)+Cx(ij)
          C_ij(2) = C_ij(2)+Cy(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, B_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+B_ij
          K(:,ij) = K(:,ij)-A_ij-B_ij
          K(:,ji) = K(:,ji)+A_ij-B_ij
          K(:,jj) = K(:,jj)-A_ij+B_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9Diag_2D


    !**************************************************************
    ! Assemble block-diagonal high-order Galerkin operator K in 3D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doGalerkinMat9Diag_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, Cx, Cy, Cz, u, K)
      
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
      REAL(DP), DIMENSION(NVAR)   :: A_ij,B_ij
      REAL(DP), DIMENSION(NDIM3D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1) = C_ij(1)+Cx(ij)
          C_ij(2) = C_ij(2)+Cy(ij)
          C_ij(3) = C_ij(3)+Cz(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, B_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+B_ij
          K(:,ij) = K(:,ij)-A_ij-B_ij
          K(:,ji) = K(:,ji)+A_ij-B_ij
          K(:,jj) = K(:,jj)-A_ij+B_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9Diag_3D

    
    !**************************************************************
    ! Assemble high-order Galerkin operator K in 1D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doGalerkinMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute averaged coefficients for boundary
          C_ij(1) = C_ij(1)+Cx(ij)
          
          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, B_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+B_ij
          K(:,ij) = K(:,ij)-A_ij-B_ij
          K(:,ji) = K(:,ji)+A_ij-B_ij
          K(:,jj) = K(:,jj)-A_ij+B_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7_1D

    
    !**************************************************************
    ! Assemble high-order Galerkin operator K in 2D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doGalerkinMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, Cy, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute averaged coefficients for boundary
          C_ij(1) = C_ij(1)+Cx(ij)
          C_ij(2) = C_ij(2)+Cy(ij)
          
          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, B_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+B_ij
          K(:,ij) = K(:,ij)-A_ij-B_ij
          K(:,ji) = K(:,ji)+A_ij-B_ij
          K(:,jj) = K(:,jj)-A_ij+B_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7_2D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 3D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doGalerkinMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, Cy, Cz, u, K)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: K
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute averaged coefficients for boundary
          C_ij(1) = C_ij(1)+Cx(ij)
          C_ij(2) = C_ij(2)+Cy(ij)
          C_ij(3) = C_ij(3)+Cz(ij)
          
          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, B_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+B_ij
          K(:,ij) = K(:,ij)-A_ij-B_ij
          K(:,ji) = K(:,ji)+A_ij-B_ij
          K(:,jj) = K(:,jj)-A_ij+B_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7_3D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 1D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doGalerkinMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, Cx, u, K)
      
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          
          ! Compute averaged coefficients for boundary
          C_ij(1) = C_ij(1)+Cx(ij)
          
          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, B_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+B_ij
          K(:,ij) = K(:,ij)-A_ij-B_ij
          K(:,ji) = K(:,ji)+A_ij-B_ij
          K(:,jj) = K(:,jj)-A_ij+B_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9_1D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 2D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doGalerkinMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, Cx, Cy, u, K)
      
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          
          ! Compute averaged coefficients for boundary
          C_ij(1) = C_ij(1)+Cx(ij)
          C_ij(2) = C_ij(2)+Cy(ij)
          
          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, B_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+B_ij
          K(:,ij) = K(:,ij)-A_ij-B_ij
          K(:,ji) = K(:,ji)+A_ij-B_ij
          K(:,jj) = K(:,jj)-A_ij+B_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9_2D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 3D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doGalerkinMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, Cx, Cy, Cz, u, K)
      
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          
          ! Compute averaged coefficients for boundary
          C_ij(1) = C_ij(1)+Cx(ij)
          C_ij(2) = C_ij(2)+Cy(ij)
          C_ij(3) = C_ij(3)+Cz(ij)
          
          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, B_ij)
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii)+A_ij+B_ij
          K(:,ij) = K(:,ij)-A_ij-B_ij
          K(:,ji) = K(:,ji)+A_ij-B_ij
          K(:,jj) = K(:,jj)-A_ij+B_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9_3D


    !**************************************************************
    ! Assemble block-diagonal low-order operator L 
    ! with scalar dissipation in 1D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doScalarDissMat7Diag_1D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij
      REAL(DP), DIMENSION(1)      :: d_ij
      REAL(DP), DIMENSION(NDIM1D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          
          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, d_ij)

          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij-d_ij(1)
          L(:,ij) = L(:,ij)-A_ij+d_ij(1)
          L(:,ji) = L(:,ji)+A_ij+d_ij(1)
          L(:,jj) = L(:,jj)-A_ij-d_ij(1)
        END DO
      END DO
    END SUBROUTINE doScalarDissMat7Diag_1D


    !**************************************************************
    ! Assemble block-diagonal low-order operator L 
    ! with scalar dissipation in 2D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doScalarDissMat7Diag_2D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, Cy, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij
      REAL(DP), DIMENSION(1)      :: d_ij
      REAL(DP), DIMENSION(NDIM2D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          
          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, d_ij)

          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij-d_ij(1)
          L(:,ij) = L(:,ij)-A_ij+d_ij(1)
          L(:,ji) = L(:,ji)+A_ij+d_ij(1)
          L(:,jj) = L(:,jj)-A_ij-d_ij(1)
        END DO
      END DO
    END SUBROUTINE doScalarDissMat7Diag_2D


    !**************************************************************
    ! Assemble block-diagonal low-order operator L 
    ! with scalar dissipation in 3D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doScalarDissMat7Diag_3D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, Cy, Cz, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij
      REAL(DP), DIMENSION(1)      :: d_ij
      REAL(DP), DIMENSION(NDIM3D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          
          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, d_ij)

          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij-d_ij(1)
          L(:,ij) = L(:,ij)-A_ij+d_ij(1)
          L(:,ji) = L(:,ji)+A_ij+d_ij(1)
          L(:,jj) = L(:,jj)-A_ij-d_ij(1)
        END DO
      END DO
    END SUBROUTINE doScalarDissMat7Diag_3D
    

    !**************************************************************
    ! Assemble block-diagonal low-order operator L
    ! with scalar dissipation in 1D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doScalarDissMat9Diag_1D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij
      REAL(DP), DIMENSION(1)      :: d_ij
      REAL(DP), DIMENSION(NDIM1D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          
          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, d_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij-d_ij(1)
          L(:,ij) = L(:,ij)-A_ij+d_ij(1)
          L(:,ji) = L(:,ji)+A_ij+d_ij(1)
          L(:,jj) = L(:,jj)-A_ij-d_ij(1)
        END DO
      END DO
    END SUBROUTINE doScalarDissMat9Diag_1D


    !**************************************************************
    ! Assemble block-diagonal low-order operator L
    ! with scalar dissipation in 2D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doScalarDissMat9Diag_2D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, Cy, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij
      REAL(DP), DIMENSION(1)      :: d_ij
      REAL(DP), DIMENSION(NDIM2D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          
          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, d_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij-d_ij(1)
          L(:,ij) = L(:,ij)-A_ij+d_ij(1)
          L(:,ji) = L(:,ji)+A_ij+d_ij(1)
          L(:,jj) = L(:,jj)-A_ij-d_ij(1)
        END DO
      END DO
    END SUBROUTINE doScalarDissMat9Diag_2D

    
    !**************************************************************
    ! Assemble block-diagonal low-order operator L
    ! with scalar dissipation in 3D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doScalarDissMat9Diag_3D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, Cy, Cz, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij
      REAL(DP), DIMENSION(1)      :: d_ij
      REAL(DP), DIMENSION(NDIM3D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          
          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, d_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij-d_ij(1)
          L(:,ij) = L(:,ij)-A_ij+d_ij(1)
          L(:,ji) = L(:,ji)+A_ij+d_ij(1)
          L(:,jj) = L(:,jj)-A_ij-d_ij(1)
        END DO
      END DO
    END SUBROUTINE doScalarDissMat9Diag_3D


    !**************************************************************
    ! Assemble low-order operator L with scalar dissipation in 1D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doScalarDissMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
      INTEGER, DIMENSION(NVAR)       :: idiag 
      INTEGER(PREC_MATIDX)           :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)           :: i,j
      INTEGER                        :: ivar

      ! Initialize diagonal entries
      idiag = (/ (NVAR*(ivar-1)+ivar, ivar = 1, NVAR) /)

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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, d_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij
          L(:,ij) = L(:,ij)-A_ij
          L(:,ji) = L(:,ji)+A_ij
          L(:,jj) = L(:,jj)-A_ij

          ! Apply the scalar dissipation I*d_ij 
          L(idiag,ii) = L(idiag,ii)-d_ij(1)
          L(idiag,ij) = L(idiag,ij)+d_ij(1)
          L(idiag,ji) = L(idiag,ji)+d_ij(1)
          L(idiag,jj) = L(idiag,jj)-d_ij(1)
        END DO
      END DO
    END SUBROUTINE doScalarDissMat7_1D


    !**************************************************************
    ! Assemble low-order operator L with scalar dissipation in 2D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doScalarDissMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, Cy, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
      INTEGER, DIMENSION(NVAR)       :: idiag 
      INTEGER(PREC_MATIDX)           :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)           :: i,j
      INTEGER                        :: ivar

      ! Initialize diagonal entries
      idiag = (/ (NVAR*(ivar-1)+ivar, ivar = 1, NVAR) /)

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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, d_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij
          L(:,ij) = L(:,ij)-A_ij
          L(:,ji) = L(:,ji)+A_ij
          L(:,jj) = L(:,jj)-A_ij

          ! Apply the scalar dissipation I*d_ij 
          L(idiag,ii) = L(idiag,ii)-d_ij(1)
          L(idiag,ij) = L(idiag,ij)+d_ij(1)
          L(idiag,ji) = L(idiag,ji)+d_ij(1)
          L(idiag,jj) = L(idiag,jj)-d_ij(1)
        END DO
      END DO
    END SUBROUTINE doScalarDissMat7_2D
    

    !**************************************************************
    ! Assemble low-order operator L with scalar dissipation in 3D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doScalarDissMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, Cy, Cz, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
      INTEGER, DIMENSION(NVAR)       :: idiag 
      INTEGER(PREC_MATIDX)           :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)           :: i,j
      INTEGER                        :: ivar

      ! Initialize diagonal entries
      idiag = (/ (NVAR*(ivar-1)+ivar, ivar = 1, NVAR) /)

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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, d_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij
          L(:,ij) = L(:,ij)-A_ij
          L(:,ji) = L(:,ji)+A_ij
          L(:,jj) = L(:,jj)-A_ij

          ! Apply the scalar dissipation I*d_ij 
          L(idiag,ii) = L(idiag,ii)-d_ij(1)
          L(idiag,ij) = L(idiag,ij)+d_ij(1)
          L(idiag,ji) = L(idiag,ji)+d_ij(1)
          L(idiag,jj) = L(idiag,jj)-d_ij(1)
        END DO
      END DO
    END SUBROUTINE doScalarDissMat7_3D


    !**************************************************************
    ! Assemble low-order operator L with scalar dissipation in 1D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doScalarDissMat9_1D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
      INTEGER, DIMENSION(NVAR)       :: idiag 
      INTEGER(PREC_MATIDX)           :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)           :: i,j
      INTEGER                        :: ivar

      ! Initialize diagonal entries
      idiag = (/ (NVAR*(ivar-1)+ivar, ivar = 1, NVAR) /)

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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, d_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij
          L(:,ij) = L(:,ij)-A_ij
          L(:,ji) = L(:,ji)+A_ij
          L(:,jj) = L(:,jj)-A_ij

          ! Apply the scalar dissipation I*d_ij 
          L(idiag,ii) = L(idiag,ii)-d_ij(1)
          L(idiag,ij) = L(idiag,ij)+d_ij(1)
          L(idiag,ji) = L(idiag,ji)+d_ij(1)
          L(idiag,jj) = L(idiag,jj)-d_ij(1)
        END DO
      END DO
    END SUBROUTINE doScalarDissMat9_1D


    !**************************************************************
    ! Assemble low-order operator L with scalar dissipation in 2D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doScalarDissMat9_2D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, Cy, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
      INTEGER, DIMENSION(NVAR)       :: idiag 
      INTEGER(PREC_MATIDX)           :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)           :: i,j
      INTEGER                        :: ivar

      ! Initialize diagonal entries
      idiag = (/ (NVAR*(ivar-1)+ivar, ivar = 1, NVAR) /)

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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, d_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij
          L(:,ij) = L(:,ij)-A_ij
          L(:,ji) = L(:,ji)+A_ij
          L(:,jj) = L(:,jj)-A_ij

          ! Apply the scalar dissipation I*d_ij 
          L(idiag,ii) = L(idiag,ii)-d_ij(1)
          L(idiag,ij) = L(idiag,ij)+d_ij(1)
          L(idiag,ji) = L(idiag,ji)+d_ij(1)
          L(idiag,jj) = L(idiag,jj)-d_ij(1)
        END DO
      END DO
    END SUBROUTINE doScalarDissMat9_2D


    !**************************************************************
    ! Assemble low-order operator L with scalar dissipation in 3D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doScalarDissMat9_3D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, Cy, Cz, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
      INTEGER, DIMENSION(NVAR)       :: idiag 
      INTEGER(PREC_MATIDX)           :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)           :: i,j
      INTEGER                        :: ivar

      ! Initialize diagonal entries
      idiag = (/ (NVAR*(ivar-1)+ivar, ivar = 1, NVAR) /)

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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, d_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij
          L(:,ij) = L(:,ij)-A_ij
          L(:,ji) = L(:,ji)+A_ij
          L(:,jj) = L(:,jj)-A_ij

          ! Apply the scalar dissipation I*d_ij 
          L(idiag,ii) = L(idiag,ii)-d_ij(1)
          L(idiag,ij) = L(idiag,ij)+d_ij(1)
          L(idiag,ji) = L(idiag,ji)+d_ij(1)
          L(idiag,jj) = L(idiag,jj)-d_ij(1)
        END DO
      END DO
    END SUBROUTINE doScalarDissMat9_3D

    
    !**************************************************************
    ! Assemble block-diagonal low-order operator L 
    ! with tensorial dissipation in 1D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doTensorDissMat7Diag_1D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij
      REAL(DP), DIMENSION(NVAR)   :: D_ij
      REAL(DP), DIMENSION(NDIM1D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, D_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij-D_ij
          L(:,ij) = L(:,ij)-A_ij+D_ij
          L(:,ji) = L(:,ji)+A_ij+D_ij
          L(:,jj) = L(:,jj)-A_ij-D_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat7Diag_1D


    !**************************************************************
    ! Assemble block-diagonal low-order operator L 
    ! with tensorial dissipation in 2D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doTensorDissMat7Diag_2D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, Cy, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij
      REAL(DP), DIMENSION(NVAR)   :: D_ij
      REAL(DP), DIMENSION(NDIM2D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, D_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij-D_ij
          L(:,ij) = L(:,ij)-A_ij+D_ij
          L(:,ji) = L(:,ji)+A_ij+D_ij
          L(:,jj) = L(:,jj)-A_ij-D_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat7Diag_2D


    !**************************************************************
    ! Assemble block-diagonal low-order operator L 
    ! with tensorial dissipation in 3D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doTensorDissMat7Diag_3D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, Cy, Cz, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij
      REAL(DP), DIMENSION(NVAR)   :: D_ij
      REAL(DP), DIMENSION(NDIM3D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, D_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij-D_ij
          L(:,ij) = L(:,ij)-A_ij+D_ij
          L(:,ji) = L(:,ji)+A_ij+D_ij
          L(:,jj) = L(:,jj)-A_ij-D_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat7Diag_3D


    !**************************************************************
    ! Assemble block-diagonal low-order operator L 
    ! with tensorial dissipation in 1D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doTensorDissMat9Diag_1D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij
      REAL(DP), DIMENSION(NVAR)   :: D_ij
      REAL(DP), DIMENSION(NDIM1D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, D_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij-D_ij
          L(:,ij) = L(:,ij)-A_ij+D_ij
          L(:,ji) = L(:,ji)+A_ij+D_ij
          L(:,jj) = L(:,jj)-A_ij-D_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat9Diag_1D


    !**************************************************************
    ! Assemble block-diagonal low-order operator L 
    ! with tensorial dissipation in 2D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doTensorDissMat9Diag_2D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, Cy, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij
      REAL(DP), DIMENSION(NVAR)   :: D_ij
      REAL(DP), DIMENSION(NDIM2D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, D_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij-D_ij
          L(:,ij) = L(:,ij)-A_ij+D_ij
          L(:,ji) = L(:,ji)+A_ij+D_ij
          L(:,jj) = L(:,jj)-A_ij-D_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat9Diag_2D


    !**************************************************************
    ! Assemble block-diagonal low-order operator L 
    ! with tensorial dissipation in 3D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doTensorDissMat9Diag_3D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, Cy, Cz, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR)   :: A_ij
      REAL(DP), DIMENSION(NVAR)   :: D_ij
      REAL(DP), DIMENSION(NDIM3D) :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, D_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij-D_ij
          L(:,ij) = L(:,ij)-A_ij+D_ij
          L(:,ji) = L(:,ji)+A_ij+D_ij
          L(:,jj) = L(:,jj)-A_ij-D_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat9Diag_3D


    !**************************************************************
    ! Assemble low-order operator L with tensorial dissipation in 1D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doTensorDissMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          
          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, D_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij-D_ij
          L(:,ij) = L(:,ij)-A_ij+D_ij
          L(:,ji) = L(:,ji)+A_ij+D_ij
          L(:,jj) = L(:,jj)-A_ij-D_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat7_1D


    !**************************************************************
    ! Assemble low-order operator L with tensorial dissipation in 2D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doTensorDissMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, Cy, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          
          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, D_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij-D_ij
          L(:,ij) = L(:,ij)-A_ij+D_ij
          L(:,ji) = L(:,ji)+A_ij+D_ij
          L(:,jj) = L(:,jj)-A_ij-D_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat7_2D


    !**************************************************************
    ! Assemble low-order operator L with tensorial dissipation in 3D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doTensorDissMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, Cy, Cz, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          
          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, D_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij-D_ij
          L(:,ij) = L(:,ij)-A_ij+D_ij
          L(:,ji) = L(:,ji)+A_ij+D_ij
          L(:,jj) = L(:,jj)-A_ij-D_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat7_3D


    !**************************************************************
    ! Assemble low-order operator L with tensorial dissipation in 1D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doTensorDissMat9_1D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, D_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij-D_ij
          L(:,ij) = L(:,ij)-A_ij+D_ij
          L(:,ji) = L(:,ji)+A_ij+D_ij
          L(:,jj) = L(:,jj)-A_ij-D_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat9_1D


    !**************************************************************
    ! Assemble low-order operator L with tensorial dissipation in 2D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doTensorDissMat9_2D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, Cy, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, D_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij-D_ij
          L(:,ij) = L(:,ij)-A_ij+D_ij
          L(:,ji) = L(:,ji)+A_ij+D_ij
          L(:,jj) = L(:,jj)-A_ij-D_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat9_2D

    
    !**************************************************************
    ! Assemble low-order operator L with tensorial dissipation in 3D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doTensorDissMat9_3D(Kld, Kcol, Kdiagonal,&
        Ksep, NEQ, NVAR, Cx, Cy, Cz, u, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), DIMENSION(NVAR*NVAR,*), INTENT(INOUT)   :: L
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, D_ij)
          
          ! Assemble the global operator
          L(:,ii) = L(:,ii)+A_ij-D_ij
          L(:,ij) = L(:,ij)-A_ij+D_ij
          L(:,ji) = L(:,ji)+A_ij+D_ij
          L(:,jj) = L(:,jj)-A_ij-D_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat9_3D
  END SUBROUTINE gfsys_buildDivOperatorScalar

  ! *****************************************************************************

!<subroutine>
  
  SUBROUTINE gfsys_buildResBlockGalerkin(RmatrixC, ru,&
      fcb_getRoeMatrix, dscale, rres)

!<description>
    ! This subroutine assembles the residual vector for standard
    ! Galerkin scheme. If the vectors contain only one block,
    ! then the scalar counterpart of this routine is called with
    ! the scalar subvectors.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: RmatrixC

    ! solution vector
    TYPE(t_vectorBlock), INTENT(IN)                :: ru

    ! scaling factor
    REAL(DP), INTENT(IN)                           :: dscale

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
    REAL(DP), DIMENSION(:), POINTER :: p_Cx,p_Cy,p_Cz,p_u,p_res
    INTEGER :: h_Ksep

    
    ! Check if block vectors contain only one block.
    IF ((ru%nblocks .EQ. 1) .AND. (rres%nblocks .EQ. 1) ) THEN
      CALL gfsys_buildResScalarGalerkin(RmatrixC, ru%RvectorBlock(1),&
          fcb_getRoeMatrix, dscale, rres%RvectorBlock(1))
      RETURN       
    END IF

    ! Check if vectors are compatible
    CALL lsysbl_isVectorCompatible(ru, rres)

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
        CALL doGalerkinMat7_1D(p_Kld, p_Kcol, p_Ksep,&
            RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, dscale, p_res)
        
      CASE (NDIM2D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
         CALL doGalerkinMat7_2D(p_Kld, p_Kcol, p_Ksep,&
            RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, dscale, p_res)
        
      CASE (NDIM3D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
        CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
        CALL doGalerkinMat7_3D(p_Kld, p_Kcol, p_Ksep,&
            RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, dscale, p_res)
        
      CASE DEFAULT
        CALL output_line('Unsupported spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlockGalerkin')
        CALL sys_halt()
      END SELECT
    

    CASE (LSYSSC_MATRIX9)

      ! Set pointer
      CALL lsyssc_getbase_Kdiagonal(RmatrixC(1), p_Kdiagonal)

      ! How many dimensions do we have?
      SELECT CASE(SIZE(RmatrixC,1))
      CASE (NDIM1D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL doGalerkinMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
            RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, dscale, p_res)
        
      CASE (NDIM2D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
         CALL doGalerkinMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
            RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, dscale, p_res)
        
      CASE (NDIM3D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
        CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
        CALL doGalerkinMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
            RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, dscale, p_res)
        
      CASE DEFAULT
        CALL output_line('Unsupported spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlockGalerkin')
        CALL sys_halt()
      END SELECT
      

    CASE DEFAULT
      CALL output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlockGalerkin')
      CALL sys_halt()
    END SELECT
    
    ! Release Ksep
    CALL storage_free(h_Ksep)
    
  CONTAINS
    
    ! Here, the working routines follow
    
    !**************************************************************
    ! Assemble residual for high-order Galerkin operator in 1D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doGalerkinMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rB_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute solution difference
          diff = u_j-u_i

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1) = C_ij(1)+Cx(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, B_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, B_ij, NVAR, diff, 1, 0._DP, rB_ij, 1)

          ! Assemble residual vector
          res(i,:) = res(i,:)+rA_ij+rB_ij
          res(j,:) = res(j,:)+rA_ij-rB_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7_1D


    !**************************************************************
    ! Assemble residual for high-order Galerkin operator in 2D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doGalerkinMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rB_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute solution difference
          diff = u_j-u_i

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1) = C_ij(1)+Cx(ij)
          C_ij(2) = C_ij(2)+Cy(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, B_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, B_ij, NVAR, diff, 1, 0._DP, rB_ij, 1)

          ! Assemble residual vector
          res(i,:) = res(i,:)+rA_ij+rB_ij
          res(j,:) = res(j,:)+rA_ij-rB_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7_2D

    
    !**************************************************************
    ! Assemble residual for high-order Galerkin operator in 3D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doGalerkinMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rB_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute solution difference
          diff = u_j-u_i

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1) = C_ij(1)+Cx(ij)
          C_ij(2) = C_ij(2)+Cy(ij)
          C_ij(3) = C_ij(3)+Cz(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, B_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, B_ij, NVAR, diff, 1, 0._DP, rB_ij, 1)

          ! Assemble residual vector
          res(i,:) = res(i,:)+rA_ij+rB_ij
          res(j,:) = res(j,:)+rA_ij-rB_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7_3D

    
    !**************************************************************
    ! Assemble residual for high-order Galerkin operator in 1D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doGalerkinMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rB_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute solution difference
          diff = u_j-u_i

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1) = C_ij(1)+Cx(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, B_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, B_ij, NVAR, diff, 1, 0._DP, rB_ij, 1)

          ! Assemble residual vector
          res(i,:) = res(i,:)+rA_ij+rB_ij
          res(j,:) = res(j,:)+rA_ij-rB_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9_1D


    !**************************************************************
    ! Assemble residual for high-order Galerkin operator in 2D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doGalerkinMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rB_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute solution difference
          diff = u_j-u_i

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1) = C_ij(1)+Cx(ij)
          C_ij(2) = C_ij(2)+Cy(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, B_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, B_ij, NVAR, diff, 1, 0._DP, rB_ij, 1)

          ! Assemble residual vector
          res(i,:) = res(i,:)+rA_ij+rB_ij
          res(j,:) = res(j,:)+rA_ij-rB_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9_2D


    !**************************************************************
    ! Assemble residual for high-order Galerkin operator in 3D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doGalerkinMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rB_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute solution difference
          diff = u_j-u_i

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1) = C_ij(1)+Cx(ij)
          C_ij(2) = C_ij(2)+Cy(ij)
          C_ij(3) = C_ij(3)+Cz(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, B_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, B_ij, NVAR, diff, 1, 0._DP, rB_ij, 1)

          ! Assemble residual vector
          res(i,:) = res(i,:)+rA_ij+rB_ij
          res(j,:) = res(j,:)+rA_ij-rB_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9_3D
  END SUBROUTINE gfsys_buildResBlockGalerkin
  
  ! *****************************************************************************

!<subroutine>
  
  SUBROUTINE gfsys_buildResScalarGalerkin(RmatrixC, ru,&
      fcb_getRoeMatrix, dscale, rres)

!<description>
    ! This subroutine assembles the residual vector for standard
    ! Galerkin scheme. Note that the vectors are required as scalar
    ! vectors which are stored in the interleave format.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: RmatrixC

    ! solution vector
    TYPE(t_vectorScalar), INTENT(IN)               :: ru

    ! scaling factor
    REAL(DP), INTENT(IN)                           :: dscale

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
    REAL(DP), DIMENSION(:), POINTER :: p_Cx,p_Cy,p_Cz,p_u,p_res
    INTEGER :: h_Ksep

    ! Check if vectors are compatible
    CALL lsyssc_isVectorCompatible(ru, rres)

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
        CALL doGalerkinMat7_1D(p_Kld, p_Kcol, p_Ksep,&
            RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, dscale, p_res)
        
      CASE (NDIM2D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
         CALL doGalerkinMat7_2D(p_Kld, p_Kcol, p_Ksep,&
            RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, dscale, p_res)
        
      CASE (NDIM3D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
        CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
        CALL doGalerkinMat7_3D(p_Kld, p_Kcol, p_Ksep,&
            RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, dscale, p_res)
        
      CASE DEFAULT
        CALL output_line('Unsupported spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarGalerkin')
        CALL sys_halt()
      END SELECT
    

    CASE (LSYSSC_MATRIX9)

      ! Set pointer
      CALL lsyssc_getbase_Kdiagonal(RmatrixC(1), p_Kdiagonal)

      ! How many dimensions do we have?
      SELECT CASE(SIZE(RmatrixC,1))
      CASE (NDIM1D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL doGalerkinMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
            RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, dscale, p_res)
        
      CASE (NDIM2D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
         CALL doGalerkinMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
            RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, dscale, p_res)
        
      CASE (NDIM3D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
        CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
        CALL doGalerkinMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
            RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, dscale, p_res)
        
      CASE DEFAULT
        CALL output_line('Unsupported spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarGalerkin')
        CALL sys_halt()
      END SELECT
      

    CASE DEFAULT
      CALL output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarGalerkin')
      CALL sys_halt()
    END SELECT

    ! Release Ksep
    CALL storage_free(h_Ksep)

  CONTAINS

    ! Here, the working routines follow
    
    !**************************************************************
    ! Assemble residual for high-order Galerkin operator in 1D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doGalerkinMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rB_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute solution difference
          diff = u(:,j)-u(:,i)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1) = C_ij(1)+Cx(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, B_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, B_ij, NVAR, diff, 1, 0._DP, rB_ij, 1)

          ! Assemble residual vector
          res(:,i) = res(:,i)+rA_ij+rB_ij
          res(:,j) = res(:,j)+rA_ij-rB_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7_1D


    !**************************************************************
    ! Assemble residual for high-order Galerkin operator in 2D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doGalerkinMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rB_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute solution difference
          diff = u(:,j)-u(:,i)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1) = C_ij(1)+Cx(ij)
          C_ij(2) = C_ij(2)+Cy(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, B_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, B_ij, NVAR, diff, 1, 0._DP, rB_ij, 1)

          ! Assemble residual vector
          res(:,i) = res(:,i)+rA_ij+rB_ij
          res(:,j) = res(:,j)+rA_ij-rB_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7_2D

    
    !**************************************************************
    ! Assemble residual for high-order Galerkin operator in 3D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doGalerkinMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rB_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))
          
          ! Compute solution difference
          diff = u(:,j)-u(:,i)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1) = C_ij(1)+Cx(ij)
          C_ij(2) = C_ij(2)+Cy(ij)
          C_ij(3) = C_ij(3)+Cz(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, B_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, B_ij, NVAR, diff, 1, 0._DP, rB_ij, 1)

          ! Assemble residual vector
          res(:,i) = res(:,i)+rA_ij+rB_ij
          res(:,j) = res(:,j)+rA_ij-rB_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7_3D

    
    !**************************************************************
    ! Assemble residual for high-order Galerkin operator in 1D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doGalerkinMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rB_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute solution difference
          diff = u(:,j)-u(:,i)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1) = C_ij(1)+Cx(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, B_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, B_ij, NVAR, diff, 1, 0._DP, rB_ij, 1)

          ! Assemble residual vector
          res(:,i) = res(:,i)+rA_ij+rB_ij
          res(:,j) = res(:,j)+rA_ij-rB_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9_1D


    !**************************************************************
    ! Assemble residual for high-order Galerkin operator in 2D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doGalerkinMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, u, dscale, res)
      
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: res
      
      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rB_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute solution difference
          diff = u(:,j)-u(:,i)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1) = C_ij(1)+Cx(ij)
          C_ij(2) = C_ij(2)+Cy(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, B_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, B_ij, NVAR, diff, 1, 0._DP, rB_ij, 1)

          ! Assemble residual vector
          res(:,i) = res(:,i)+rA_ij+rB_ij
          res(:,j) = res(:,j)+rA_ij-rB_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9_2D


    !**************************************************************
    ! Assemble residual for high-order Galerkin operator in 3D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doGalerkinMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,B_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rB_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))
          
          ! Compute solution difference
          diff = u(:,j)-u(:,i)

          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute averaged coefficients for boundary contribution
          C_ij(1) = C_ij(1)+Cx(ij)
          C_ij(2) = C_ij(2)+Cy(ij)
          C_ij(3) = C_ij(3)+Cz(ij)

          ! Compute local Roe matrix for boundary contribution
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, B_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, B_ij, NVAR, diff, 1, 0._DP, rB_ij, 1)

          ! Assemble residual vector
          res(:,i) = res(:,i)+rA_ij+rB_ij
          res(:,j) = res(:,j)+rA_ij-rB_ij
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9_3D
  END SUBROUTINE gfsys_buildResScalarGalerkin

  ! *****************************************************************************

!<subroutine>
  
  SUBROUTINE gfsys_buildResBlockUpwind(RmatrixC, ru,&
      fcb_getRoeMatrix, fcb_getDissipation, idissipation, dscale, rres)

!<description>
    ! This subroutine assembles the residual vector for discrete
    ! upwinding schemes. If the vectors contain only one block,
    ! then the scalar counterpart of this routine is called with
    ! the scalar subvectors.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: RmatrixC

    ! solution vector
    TYPE(t_vectorBlock), INTENT(IN)                :: ru

    ! type of dissipation
    !   AFCSTAB_SCALARDISSIPATION: scalar dissipation
    !   AFCSTAB_TENSORDISSIPATION: tensorial dissipation
    INTEGER, INTENT(IN)                            :: idissipation
    
    ! scaling factor
    REAL(DP), INTENT(IN)                           :: dscale

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
    REAL(DP), DIMENSION(:), POINTER :: p_Cx,p_Cy,p_Cz,p_u,p_res
    INTEGER :: h_Ksep


    ! Check if block vectors contain only one block.
    IF ((ru%nblocks .EQ. 1) .AND. (rres%nblocks .EQ. 1) ) THEN
      CALL gfsys_buildResScalarUpwind(RmatrixC, ru%RvectorBlock(1),&
          fcb_getRoeMatrix, fcb_getDissipation, idissipation, dscale, rres%RvectorBlock(1))
      RETURN       
    END IF
    
    ! Check if vectors are compatible
    CALL lsysbl_isVectorCompatible(ru, rres)

    ! Set pointers
    CALL lsyssc_getbase_Kld   (RmatrixC(1), p_Kld)
    CALL lsyssc_getbase_Kcol  (RmatrixC(1), p_Kcol)
    CALL lsysbl_getbase_double(ru,   p_u)
    CALL lsysbl_getbase_double(rres, p_res)

    ! Create diagonal Ksep=Kld
    h_Ksep = ST_NOHANDLE
    CALL storage_copy(RmatrixC(1)%h_Kld, h_Ksep)
    CALL storage_getbase_int(h_Ksep, p_Ksep, RmatrixC(1)%NEQ+1)
    

    ! What kind of dissipation should be applied?
    SELECT CASE(idissipation)
    CASE (AFCSTAB_SCALARDISSIPATION)
      
      ! What kind of matrix are we?
      SELECT CASE(RmatrixC(1)%cmatrixFormat)
      CASE(LSYSSC_MATRIX7)
        
        ! How many dimensions do we have?
        SELECT CASE(SIZE(RmatrixC,1))
        CASE (NDIM1D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL doScalarDissMat7_1D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, dscale, p_res)
          
        CASE (NDIM2D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL doScalarDissMat7_2D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, dscale, p_res)
          
        CASE (NDIM3D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
          CALL doScalarDissMat7_3D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, dscale, p_res)
          
        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlockUpwind')
          CALL sys_halt()
        END SELECT
        
        
      CASE (LSYSSC_MATRIX9)
        
        ! Set pointer
        CALL lsyssc_getbase_Kdiagonal(RmatrixC(1), p_Kdiagonal)
        
        ! How many dimensions do we have?
        SELECT CASE(SIZE(RmatrixC,1))
        CASE (NDIM1D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL doScalarDissMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, dscale, p_res)
          
        CASE (NDIM2D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL doScalarDissMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, dscale, p_res)
          
        CASE (NDIM3D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
          CALL doScalarDissMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, dscale, p_res)
          
        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlockUpwind')
          CALL sys_halt()
        END SELECT
      
        
      CASE DEFAULT
        CALL output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarUpwind')
        CALL sys_halt()
      END SELECT


    CASE (AFCSTAB_TENSORDISSIPATION)
      
      ! What kind of matrix are we?
      SELECT CASE(RmatrixC(1)%cmatrixFormat)
      CASE(LSYSSC_MATRIX7)
        
        ! How many dimensions do we have?
        SELECT CASE(SIZE(RmatrixC,1))
        CASE (NDIM1D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL doTensorDissMat7_1D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, dscale, p_res)
          
        CASE (NDIM2D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL doTensorDissMat7_2D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, dscale, p_res)
          
        CASE (NDIM3D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
          CALL doTensorDissMat7_3D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, dscale, p_res)
          
        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlockUpwind')
          CALL sys_halt()
        END SELECT
        
        
      CASE (LSYSSC_MATRIX9)
        
        ! Set pointer
        CALL lsyssc_getbase_Kdiagonal(RmatrixC(1), p_Kdiagonal)
        
        ! How many dimensions do we have?
        SELECT CASE(SIZE(RmatrixC,1))
        CASE (NDIM1D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL doTensorDissMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_u, dscale, p_res)
          
        CASE (NDIM2D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL doTensorDissMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_u, dscale, p_res)
          
        CASE (NDIM3D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
          CALL doTensorDissMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, dscale, p_res)
          
        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlockUpwind')
          CALL sys_halt()
        END SELECT
      
        
      CASE DEFAULT
        CALL output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlockUpwind')
        CALL sys_halt()
      END SELECT


    CASE DEFAULT
      CALL output_line('Unsupported type of dissipation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlockUpwind')
      CALL sys_halt()
    END SELECT
    
    ! Release Ksep
    CALL storage_free(h_Ksep)
    
  CONTAINS

    ! Here, the working routines follow
    
    !**************************************************************
    ! Assemble residual for low-order operator with scalar dissipation in 1D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doScalarDissMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute solution difference
          diff = u_j-u_i

          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, d_ij)
          rD_ij = dscale*d_ij(1)*diff
          
          ! Assemble residual vector
          res(i,:) = res(i,:)+rA_ij+rD_ij
          res(j,:) = res(j,:)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doScalarDissMat7_1D


    !**************************************************************
    ! Assemble residual for low-order operator with scalar dissipation in 2D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doScalarDissMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute solution difference
          diff = u_j-u_i

          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, d_ij)
          rD_ij = dscale*d_ij(1)*diff
          
          ! Assemble residual vector
          res(i,:) = res(i,:)+rA_ij+rD_ij
          res(j,:) = res(j,:)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doScalarDissMat7_2D


    !**************************************************************
    ! Assemble residual for low-order operator with scalar dissipation in 3D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doScalarDissMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute solution difference
          diff = u_j-u_i

          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, d_ij)
          rD_ij = dscale*d_ij(1)*diff
          
          ! Assemble residual vector
          res(i,:) = res(i,:)+rA_ij+rD_ij
          res(j,:) = res(j,:)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doScalarDissMat7_3D


    !**************************************************************
    ! Assemble residual for low-order operator with scalar dissipation in 1D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doScalarDissMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, u, dscale, res)
      
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(INOUT)      :: res

      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute solution difference
          diff = u_j-u_i
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, d_ij)
          rD_ij = dscale*d_ij(1)*diff
          
          ! Assemble residual vector
          res(i,:) = res(i,:)+rA_ij+rD_ij
          res(j,:) = res(j,:)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doScalarDissMat9_1D


    !**************************************************************
    ! Assemble residual for low-order operator with scalar dissipation in 2D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doScalarDissMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, u, dscale, res)
      
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(INOUT)      :: res

      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute solution difference
          diff = u_j-u_i
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, d_ij)
          rD_ij = dscale*d_ij(1)*diff
          
          ! Assemble residual vector
          res(i,:) = res(i,:)+rA_ij+rD_ij
          res(j,:) = res(j,:)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doScalarDissMat9_2D


    !**************************************************************
    ! Assemble residual for low-order operator with scalar dissipation in 3D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doScalarDissMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, Cz, u, dscale, res)
      
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(INOUT)      :: res

      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute solution difference
          diff = u_j-u_i
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, d_ij)
          rD_ij = dscale*d_ij(1)*diff
          
          ! Assemble residual vector
          res(i,:) = res(i,:)+rA_ij+rD_ij
          res(j,:) = res(j,:)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doScalarDissMat9_3D

    
    !**************************************************************
    ! Assemble residual for low-order operator
    ! with tensorial dissipation in 1D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doTensorDissMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute solution difference
          diff = u_j-u_i
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Compute tensorial dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, D_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, D_ij, NVAR, diff, 1, 0._DP, rD_ij, 1)

          ! Assemble residual vector
          res(i,:) = res(i,:)+rA_ij+rD_ij
          res(j,:) = res(j,:)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat7_1D


    !**************************************************************
    ! Assemble residual for low-order operator
    ! with tensorial dissipation in 2D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doTensorDissMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute solution difference
          diff = u_j-u_i
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Compute tensorial dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, D_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, D_ij, NVAR, diff, 1, 0._DP, rD_ij, 1)

          ! Assemble residual vector
          res(i,:) = res(i,:)+rA_ij+rD_ij
          res(j,:) = res(j,:)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat7_2D

    
    !**************************************************************
    ! Assemble residual for low-order operator
    ! with tensorial dissipation in 3D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doTensorDissMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute solution difference
          diff = u_j-u_i
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Compute tensorial dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, D_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, D_ij, NVAR, diff, 1, 0._DP, rD_ij, 1)

          ! Assemble residual vector
          res(i,:) = res(i,:)+rA_ij+rD_ij
          res(j,:) = res(j,:)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat7_3D

    
    !**************************************************************
    ! Assemble residual for low-order operator
    ! with tensorial dissipation in 1D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doTensorDissMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, u, dscale, res)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(INOUT)      :: res

      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute solution difference
          diff = u_j-u_i
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Compute tensorial dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, D_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, D_ij, NVAR, diff, 1, 0._DP, rD_ij, 1)

          ! Assemble residual vector
          res(i,:) = res(i,:)+rA_ij+rD_ij
          res(j,:) = res(j,:)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat9_1D


    !**************************************************************
    ! Assemble residual for low-order operator
    ! with tensorial dissipation in 2D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doTensorDissMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, u, dscale, res)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(INOUT)      :: res

      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute solution difference
          diff = u_j-u_i
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Compute tensorial dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, D_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, D_ij, NVAR, diff, 1, 0._DP, rD_ij, 1)

          ! Assemble residual vector
          res(i,:) = res(i,:)+rA_ij+rD_ij
          res(j,:) = res(j,:)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat9_2D


    !**************************************************************
    ! Assemble residual for low-order operator
    ! with tensorial dissipation in 3D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doTensorDissMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, Cz, u, dscale, res)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(IN)         :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NEQ,NVAR), INTENT(INOUT)      :: res

      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute solution difference
          diff = u_j-u_i
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Compute tensorial dissipation
          CALL fcb_getDissipation(u_i, u_j, C_ij, D_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, D_ij, NVAR, diff, 1, 0._DP, rD_ij, 1)

          ! Assemble residual vector
          res(i,:) = res(i,:)+rA_ij+rD_ij
          res(j,:) = res(j,:)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat9_3D
  END SUBROUTINE gfsys_buildResBlockUpwind

  ! *****************************************************************************

!<subroutine>
  
  SUBROUTINE gfsys_buildResScalarUpwind(RmatrixC, ru,&
      fcb_getRoeMatrix, fcb_getDissipation, idissipation, dscale, rres)

!<description>
    ! This subroutine assembles the residual vector for discrete
    ! upwinding schemes. Note that the vectors are required as scalar
    ! vectors which are stored in the interleave format.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: RmatrixC

    ! solution vector
    TYPE(t_vectorScalar), INTENT(IN)               :: ru

    ! type of dissipation
    !   AFCSTAB_SCALARDISSIPATION: scalar dissipation
    !   AFCSTAB_TENSORDISSIPATION: tensorial dissipation
    INTEGER, INTENT(IN)                            :: idissipation
    
    ! scaling factor
    REAL(DP), INTENT(IN)                           :: dscale

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
    REAL(DP), DIMENSION(:), POINTER :: p_Cx,p_Cy,p_Cz,p_u,p_res
    INTEGER :: h_Ksep

    ! Check if vectors are compatible
    CALL lsyssc_isVectorCompatible(ru,rres)

    ! Set pointers
    CALL lsyssc_getbase_Kld   (RmatrixC(1), p_Kld)
    CALL lsyssc_getbase_Kcol  (RmatrixC(1), p_Kcol)
    CALL lsyssc_getbase_double(ru,   p_u)
    CALL lsyssc_getbase_double(rres, p_res)

    ! Create diagonal Ksep=Kld
    h_Ksep = ST_NOHANDLE
    CALL storage_copy(RmatrixC(1)%h_Kld, h_Ksep)
    CALL storage_getbase_int(h_Ksep, p_Ksep, RmatrixC(1)%NEQ+1)


    ! What kind of dissipation should be applied?
    SELECT CASE(idissipation)
    CASE (AFCSTAB_SCALARDISSIPATION)
      
      ! What kind of matrix are we?
      SELECT CASE(RmatrixC(1)%cmatrixFormat)
      CASE(LSYSSC_MATRIX7)
        
        ! How many dimensions do we have?
        SELECT CASE(SIZE(RmatrixC,1))
        CASE (NDIM1D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL doScalarDissMat7_1D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, dscale, p_res)
          
        CASE (NDIM2D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL doScalarDissMat7_2D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, dscale, p_res)
          
        CASE (NDIM3D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
          CALL doScalarDissMat7_3D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, dscale, p_res)
          
        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarUpwind')
          CALL sys_halt()
        END SELECT
        
        
      CASE (LSYSSC_MATRIX9)
        
        ! Set pointer
        CALL lsyssc_getbase_Kdiagonal(RmatrixC(1), p_Kdiagonal)
        
        ! How many dimensions do we have?
        SELECT CASE(SIZE(RmatrixC,1))
        CASE (NDIM1D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL doScalarDissMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, dscale, p_res)
          
        CASE (NDIM2D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL doScalarDissMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, dscale, p_res)
          
        CASE (NDIM3D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
          CALL doScalarDissMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, dscale, p_res)
          
        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarUpwind')
          CALL sys_halt()
        END SELECT
      
        
      CASE DEFAULT
        CALL output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarUpwind')
        CALL sys_halt()
      END SELECT


    CASE (AFCSTAB_TENSORDISSIPATION)
      
      ! What kind of matrix are we?
      SELECT CASE(RmatrixC(1)%cmatrixFormat)
      CASE(LSYSSC_MATRIX7)
        
        ! How many dimensions do we have?
        SELECT CASE(SIZE(RmatrixC,1))
        CASE (NDIM1D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL doTensorDissMat7_1D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, dscale, p_res)
          
        CASE (NDIM2D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL doTensorDissMat7_2D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, dscale, p_res)
          
        CASE (NDIM3D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
          CALL doTensorDissMat7_3D(p_Kld, p_Kcol, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, dscale, p_res)
          
        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarUpwind')
          CALL sys_halt()
        END SELECT
        
        
      CASE (LSYSSC_MATRIX9)
        
        ! Set pointer
        CALL lsyssc_getbase_Kdiagonal(RmatrixC(1), p_Kdiagonal)
        
        ! How many dimensions do we have?
        SELECT CASE(SIZE(RmatrixC,1))
        CASE (NDIM1D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL doTensorDissMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_u, dscale, p_res)
          
        CASE (NDIM2D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL doTensorDissMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_u, dscale, p_res)
          
        CASE (NDIM3D)
          CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
          CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
          CALL doTensorDissMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              RmatrixC(1)%NEQ, ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, dscale, p_res)
          
        CASE DEFAULT
          CALL output_line('Unsupported spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarUpwind')
          CALL sys_halt()
        END SELECT
      
        
      CASE DEFAULT
        CALL output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarUpwind')
        CALL sys_halt()
      END SELECT

      
    CASE DEFAULT
      CALL output_line('Unsupported type of dissipation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarUpwind')
      CALL sys_halt()
    END SELECT
    
    ! Release Ksep
    CALL storage_free(h_Ksep)

  CONTAINS

    ! Here, the working routines follow
    
    !**************************************************************
    ! Assemble residual for low-order operator with scalar dissipation in 1D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doScalarDissMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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

          ! Compute solution difference
          diff = u(:,j)-u(:,i)

          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, d_ij)
          rD_ij = dscale*d_ij(1)*diff
          
          ! Assemble residual vector
          res(:,i) = res(:,i)+rA_ij+rD_ij
          res(:,j) = res(:,j)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doScalarDissMat7_1D


    !**************************************************************
    ! Assemble residual for low-order operator with scalar dissipation in 2D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doScalarDissMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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

          ! Compute solution difference
          diff = u(:,j)-u(:,i)

          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, d_ij)
          rD_ij = dscale*d_ij(1)*diff
          
          ! Assemble residual vector
          res(:,i) = res(:,i)+rA_ij+rD_ij
          res(:,j) = res(:,j)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doScalarDissMat7_2D


    !**************************************************************
    ! Assemble residual for low-order operator with scalar dissipation in 3D
    ! All matrices are stored in matrix format 7

    SUBROUTINE doScalarDissMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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

          ! Compute solution difference
          diff = u(:,j)-u(:,i)

          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, d_ij)
          rD_ij = dscale*d_ij(1)*diff
          
          ! Assemble residual vector
          res(:,i) = res(:,i)+rA_ij+rD_ij
          res(:,j) = res(:,j)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doScalarDissMat7_3D
    
    
    !**************************************************************
    ! Assemble residual for low-order operator with scalar dissipation in 1D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doScalarDissMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, u, dscale, res)
      
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: res

      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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

          ! Compute solution difference
          diff = u(:,j)-u(:,i)

          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, d_ij)
          rD_ij = dscale*d_ij(1)*diff
          
          ! Assemble residual vector
          res(:,i) = res(:,i)+rA_ij+rD_ij
          res(:,j) = res(:,j)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doScalarDissMat9_1D


    !**************************************************************
    ! Assemble residual for low-order operator with scalar dissipation in 2D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doScalarDissMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, u, dscale, res)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: res

      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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

          ! Compute solution difference
          diff = u(:,j)-u(:,i)

          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, d_ij)
          rD_ij = dscale*d_ij(1)*diff
          
          ! Assemble residual vector
          res(:,i) = res(:,i)+rA_ij+rD_ij
          res(:,j) = res(:,j)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doScalarDissMat9_2D


    !**************************************************************
    ! Assemble residual for low-order operator with scalar dissipation in 3D
    ! All matrices are stored in matrix format 9

    SUBROUTINE doScalarDissMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, Cz, u, dscale, res)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: res

      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(1)         :: d_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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

          ! Compute solution difference
          diff = u(:,j)-u(:,i)

          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)

          ! Compute scalar dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, d_ij)
          rD_ij = dscale*d_ij(1)*diff
          
          ! Assemble residual vector
          res(:,i) = res(:,i)+rA_ij+rD_ij
          res(:,j) = res(:,j)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doScalarDissMat9_3D
    
    
    !**************************************************************
    ! Assemble residual for low-order operator
    ! with tensorial dissipation in 1D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doTensorDissMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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

          ! Compute solution difference
          diff = u(:,j)-u(:,i)
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Compute tensorial dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, D_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, D_ij, NVAR, diff, 1, 0._DP, rD_ij, 1)

          ! Assemble residual vector
          res(:,i) = res(:,i)+rA_ij+rD_ij
          res(:,j) = res(:,j)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat7_1D


    !**************************************************************
    ! Assemble residual for low-order operator
    ! with tensorial dissipation in 2D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doTensorDissMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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

          ! Compute solution difference
          diff = u(:,j)-u(:,i)
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Compute tensorial dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, D_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, D_ij, NVAR, diff, 1, 0._DP, rD_ij, 1)

          ! Assemble residual vector
          res(:,i) = res(:,i)+rA_ij+rD_ij
          res(:,j) = res(:,j)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat7_2D

    
    !**************************************************************
    ! Assemble residual for low-order operator
    ! with tensorial dissipation in 3D
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doTensorDissMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
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
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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

          ! Compute solution difference
          diff = u(:,j)-u(:,i)
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Compute tensorial dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, D_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, D_ij, NVAR, diff, 1, 0._DP, rD_ij, 1)

          ! Assemble residual vector
          res(:,i) = res(:,i)+rA_ij+rD_ij
          res(:,j) = res(:,j)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat7_3D


    !**************************************************************
    ! Assemble residual for low-order operator
    ! with tensorial dissipation in 1D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doTensorDissMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, u, dscale, res)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: res

      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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

          ! Compute solution difference
          diff = u(:,j)-u(:,i)
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Compute tensorial dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, D_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, D_ij, NVAR, diff, 1, 0._DP, rD_ij, 1)

          ! Assemble residual vector
          res(:,i) = res(:,i)+rA_ij+rD_ij
          res(:,j) = res(:,j)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat9_1D

    
    !**************************************************************
    ! Assemble residual for low-order operator
    ! with tensorial dissipation in 2D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doTensorDissMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, u, dscale, res)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: res

      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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

          ! Compute solution difference
          diff = u(:,j)-u(:,i)
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Compute tensorial dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, D_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, D_ij, NVAR, diff, 1, 0._DP, rD_ij, 1)

          ! Assemble residual vector
          res(:,i) = res(:,i)+rA_ij+rD_ij
          res(:,j) = res(:,j)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat9_2D

    
    !**************************************************************
    ! Assemble residual for low-order operator
    ! with tensorial dissipation in 3D
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doTensorDissMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        Cx, Cy, Cz, u, dscale, res)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER, INTENT(IN)                               :: NVAR
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(NVAR,*), INTENT(IN)           :: u
      REAL(DP), INTENT(IN)                              :: dscale
      REAL(DP), DIMENSION(NVAR,*), INTENT(INOUT)        :: res

      ! local variables
      REAL(DP), DIMENSION(NVAR*NVAR) :: A_ij,D_ij
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,rD_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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

          ! Compute solution difference
          diff = u(:,j)-u(:,i)
          
          ! Compute averaged coefficients
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          C_ij(3) = 0.5_DP*(Cz(ji)-Cz(ij))
          
          ! Compute local Roe matrix
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Compute tensorial dissipation
          CALL fcb_getDissipation(u(:,i), u(:,j), C_ij, D_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, D_ij, NVAR, diff, 1, 0._DP, rD_ij, 1)

          ! Assemble residual vector
          res(:,i) = res(:,i)+rA_ij+rD_ij
          res(:,j) = res(:,j)+rA_ij-rD_ij
        END DO
      END DO
    END SUBROUTINE doTensorDissMat9_3D
  END SUBROUTINE gfsys_buildResScalarUpwind

  ! *****************************************************************************
  
  !<subroutine>
  
  SUBROUTINE gfsys_buildResBlockTVD(RmatrixC, ru,&
      fcb_getRoeMatrix, fcb_getDissipation, fcb_getCharacteristics,&
      rafcstab, dscale, rres)

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
    REAL(DP), DIMENSION(:), POINTER :: p_Cx,p_Cy,p_Cz,p_L,p_u,p_res
    REAL(DP), DIMENSION(:), POINTER :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    INTEGER :: h_Ksep

    
    ! Check if block vectors contain only one block.
    IF ((ru%nblocks .EQ. 1) .AND. (rres%nblocks .EQ. 1) ) THEN
      CALL gfsys_buildResScalarTVD(RmatrixC, ru%RvectorBlock(1),&
          fcb_getRoeMatrix, fcb_getDissipation, fcb_getCharacteristics,&
          rafcstab, dscale, rres%RvectorBlock(1))
      RETURN       
    END IF

    ! Check if vectors are compatible
    CALL lsysbl_isVectorCompatible(ru, rres)
    CALL gfsys_isVectorCompatible(rafcstab, ru)
    
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
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,F_ij,W_ij,Lbd_ij,ka_ij,kb_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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
          
          ! Compute solution difference
          diff = u_j-u_i
          
          ! Compute averaged coefficients a_ij:=(c_ji-c_ij)/2
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute local Roe matrix and multiply it with the solution difference
          !
          !   rA_ij:=A_ij*(U_j-U_i),
          !
          ! where A_ij denotes the cumulated Roe matrix 
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Assemble high-order residual vector. Since A_ij is skew-symmetric, the 
          ! contribution of A_ji to node j is identical to that of A_ij to node i
          res(i,:) = res(i,:)+rA_ij
          res(j,:) = res(j,:)+rA_ij
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u_i, u_j, XDir1D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(1)*Lbd_ij
          kb_ij = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
          
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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)
      
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

          ! Compute averaged coefficient
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u_i, u_j, XDir1D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  C_ij(1)*Lbd_ij
          kb_ij  = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,F_ij,W_ij,Lbd_ij,ka_ij,kb_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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

          ! Compute solution difference
          diff = u_j-u_i
          
          ! Compute averaged coefficients a_ij:=(c_ji-c_ij)/2
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix and multiply it with the solution difference
          !
          !   rA_ij:=A_ij*(U_j-U_i),
          !
          ! where A_ij denotes the cumulated Roe matrix 
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Assemble high-order residual vector. Since A_ij is skew-symmetric, the 
          ! contribution of A_ji to node j is identical to that of A_ij to node i
          res(i,:) = res(i,:)+rA_ij
          res(j,:) = res(j,:)+rA_ij
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u_i, u_j, XDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(1)*Lbd_ij
          kb_ij = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
          
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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)
      
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

          ! Compute averaged coefficient
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u_i, u_j, XDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  C_ij(1)*Lbd_ij
          kb_ij  = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
          CALL fcb_getCharacteristics(u_i, u_j, YDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(2)*Lbd_ij
          kb_ij = (C_ij(2)+Cy(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
          
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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)
      
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
          
          ! Compute averaged coefficients
          C_ij(2)=0.5_DP*(Cy(ji)-Cy(ij))

          ! Compute characteristic fluxes in Y-direction
          CALL fcb_getCharacteristics(u_i, u_j, YDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(2)*Lbd_ij
          kb_ij = (C_ij(2)+Cy(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,F_ij,W_ij,Lbd_ij,ka_ij,kb_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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

          ! Compute solution difference
          diff = u_j-u_i
          
          ! Compute averaged coefficients a_ij:=(c_ji-c_ij)/2
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix and multiply it with the solution difference
          !
          !   rA_ij:=A_ij*(U_j-U_i),
          !
          ! where A_ij denotes the cumulated Roe matrix 
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Assemble high-order residual vector. Since A_ij is skew-symmetric, the 
          ! contribution of A_ji to node j is identical to that of A_ij to node i
          res(i,:) = res(i,:)+rA_ij
          res(j,:) = res(j,:)+rA_ij
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u_i, u_j, XDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(1)*Lbd_ij
          kb_ij = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
          
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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)

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

          ! Compute averaged coefficient
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u_i, u_j, XDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  C_ij(1)*Lbd_ij
          kb_ij  = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
          CALL fcb_getCharacteristics(u_i, u_j, YDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(2)*Lbd_ij
          kb_ij = (C_ij(2)+Cy(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
          
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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)

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
          
          ! Compute averaged coefficients
          C_ij(2)=0.5_DP*(Cy(ji)-Cy(ij))

          ! Compute characteristic fluxes in Y-direction
          CALL fcb_getCharacteristics(u_i, u_j, YDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(2)*Lbd_ij
          kb_ij = (C_ij(2)+Cy(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
          CALL fcb_getCharacteristics(u_i, u_j, ZDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(3)*Lbd_ij
          kb_ij = (C_ij(3)+Cz(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij

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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)

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

          ! Compute averaged coefficients
          C_ij(3)=0.5_DP*(Cz(ji)-Cz(ij))

          ! Compute characteristic fluxes in Z-direction
          CALL fcb_getCharacteristics(u_i, u_j, ZDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(3)*Lbd_ij
          kb_ij = (C_ij(3)+Cz(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,F_ij,W_ij,Lbd_ij,ka_ij,kb_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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

          ! Compute solution difference
          diff = u_j-u_i
          
          ! Compute averaged coefficients a_ij:=(c_ji-c_ij)/2
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute local Roe matrix and multiply it with the solution difference
          !
          !   rA_ij:=A_ij*(U_j-U_i),
          !
          ! where A_ij denotes the cumulated Roe matrix 
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Assemble high-order residual vector. Since A_ij is skew-symmetric, the 
          ! contribution of A_ji to node j is identical to that of A_ij to node i
          res(i,:) = res(i,:)+rA_ij
          res(j,:) = res(j,:)+rA_ij
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u_i, u_j, XDir1D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(1)*Lbd_ij
          kb_ij = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
          
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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)
      
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
          
          ! Compute averaged coefficient
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u_i, u_j, XDir1D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  C_ij(1)*Lbd_ij
          kb_ij  = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,F_ij,W_ij,Lbd_ij,ka_ij,kb_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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

          ! Compute solution difference
          diff = u_j-u_i
          
          ! Compute averaged coefficients a_ij:=(c_ji-c_ij)/2
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix and multiply it with the solution difference
          !
          !   rA_ij:=A_ij*(U_j-U_i),
          !
          ! where A_ij denotes the cumulated Roe matrix 
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Assemble high-order residual vector. Since A_ij is skew-symmetric, the 
          ! contribution of A_ji to node j is identical to that of A_ij to node i
          res(i,:) = res(i,:)+rA_ij
          res(j,:) = res(j,:)+rA_ij
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u_i, u_j, XDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(1)*Lbd_ij
          kb_ij = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
          
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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)
      
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

          ! Compute averaged coefficient
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u_i, u_j, XDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  C_ij(1)*Lbd_ij
          kb_ij  = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
          CALL fcb_getCharacteristics(u_i, u_j, YDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(2)*Lbd_ij
          kb_ij = (C_ij(2)+Cy(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
          
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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)
      
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
          
          ! Compute averaged coefficients
          C_ij(2)=0.5_DP*(Cy(ji)-Cy(ij))

          ! Compute characteristic fluxes in Y-direction
          CALL fcb_getCharacteristics(u_i, u_j, YDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(2)*Lbd_ij
          kb_ij = (C_ij(2)+Cy(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,F_ij,W_ij,Lbd_ij,ka_ij,kb_ij
      REAL(DP), DIMENSION(NVAR)      :: u_i,u_j
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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
          
          ! Compute solution difference
          diff = u_j-u_i
          
          ! Compute averaged coefficients a_ij:=(c_ji-c_ij)/2
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix and multiply it with the solution difference
          !
          !   rA_ij:=A_ij*(U_j-U_i),
          !
          ! where A_ij denotes the cumulated Roe matrix 
          CALL fcb_getRoeMatrix(u_i, u_j, C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Assemble high-order residual vector. Since A_ij is skew-symmetric, the 
          ! contribution of A_ji to node j is identical to that of A_ij to node i
          res(i,:) = res(i,:)+rA_ij
          res(j,:) = res(j,:)+rA_ij
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u_i, u_j, XDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(1)*Lbd_ij
          kb_ij = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
          
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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)

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
          
          ! Compute averaged coefficient
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u_i, u_j, XDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  C_ij(1)*Lbd_ij
          kb_ij  = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
          CALL fcb_getCharacteristics(u_i, u_j, YDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(2)*Lbd_ij
          kb_ij = (C_ij(2)+Cy(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
          
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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)

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

          ! Compute averaged coefficients
          C_ij(2)=0.5_DP*(Cy(ji)-Cy(ij))

          ! Compute characteristic fluxes in Y-direction
          CALL fcb_getCharacteristics(u_i, u_j, YDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(2)*Lbd_ij
          kb_ij = (C_ij(2)+Cy(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
          CALL fcb_getCharacteristics(u_i, u_j, ZDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(3)*Lbd_ij
          kb_ij = (C_ij(3)+Cz(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij

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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)

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

          ! Compute averaged coefficients
          C_ij(3)=0.5_DP*(Cz(ji)-Cz(ij))

          ! Compute characteristic fluxes in Z-direction
          CALL fcb_getCharacteristics(u_i, u_j, ZDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(3)*Lbd_ij
          kb_ij = (C_ij(3)+Cz(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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

  SUBROUTINE gfsys_buildResScalarTVD(RmatrixC, ru,&
      fcb_getRoeMatrix, fcb_getDissipation, fcb_getCharacteristics,&
      rafcstab, dscale, rres)
    
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
    REAL(DP), DIMENSION(:), POINTER :: p_Cx,p_Cy,p_Cz,p_L,p_u,p_res
    REAL(DP), DIMENSION(:), POINTER :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    INTEGER :: h_Ksep


    ! Check if vectors are compatible
    CALL lsyssc_isVectorCompatible(ru, rres)
    CALL gfsys_isVectorCompatible(rafcstab, ru)
    
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
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,F_ij,W_ij,Lbd_ij,ka_ij,kb_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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
                   
          ! Compute solution difference
          diff = u(:,j)-u(:,i)
          
          ! Compute averaged coefficients a_ij:=(c_ji-c_ij)/2
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute local Roe matrix and multiply it with the solution difference
          !
          !   rA_ij:=A_ij*(U_j-U_i),
          !
          ! where A_ij denotes the cumulated Roe matrix 
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Assemble high-order residual vector. Since A_ij is skew-symmetric, the 
          ! contribution of A_ji to node j is identical to that of A_ij to node i
          res(:,i) = res(:,i)+rA_ij
          res(:,j) = res(:,j)+rA_ij
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u(:,i), u(:,j), XDir1D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(1)*Lbd_ij
          kb_ij = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
          
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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)
      
      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1
          
          ! Compute averaged coefficient
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u(:,i), u(:,j), XDir1D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  C_ij(1)*Lbd_ij
          kb_ij  = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,F_ij,W_ij,Lbd_ij,ka_ij,kb_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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
                   
          ! Compute solution difference
          diff = u(:,j)-u(:,i)
          
          ! Compute averaged coefficients a_ij:=(c_ji-c_ij)/2
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix and multiply it with the solution difference
          !
          !   rA_ij:=A_ij*(U_j-U_i),
          !
          ! where A_ij denotes the cumulated Roe matrix 
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Assemble high-order residual vector. Since A_ij is skew-symmetric, the 
          ! contribution of A_ji to node j is identical to that of A_ij to node i
          res(:,i) = res(:,i)+rA_ij
          res(:,j) = res(:,j)+rA_ij
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u(:,i), u(:,j), XDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(1)*Lbd_ij
          kb_ij = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
          
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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)
      
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
          
          ! Compute averaged coefficient
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u(:,i), u(:,j), XDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  C_ij(1)*Lbd_ij
          kb_ij  = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
          CALL fcb_getCharacteristics(u(:,i), u(:,j), YDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(2)*Lbd_ij
          kb_ij = (C_ij(2)+Cy(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
          
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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)
      
      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute averaged coefficients
          C_ij(2)=0.5_DP*(Cy(ji)-Cy(ij))

          ! Compute characteristic fluxes in Y-direction
          CALL fcb_getCharacteristics(u(:,i), u(:,j), YDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(2)*Lbd_ij
          kb_ij = (C_ij(2)+Cy(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,F_ij,W_ij,Lbd_ij,ka_ij,kb_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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
                   
          ! Compute solution difference
          diff = u(:,j)-u(:,i)
          
          ! Compute averaged coefficients a_ij:=(c_ji-c_ij)/2
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix and multiply it with the solution difference
          !
          !   rA_ij:=A_ij*(U_j-U_i),
          !
          ! where A_ij denotes the cumulated Roe matrix 
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Assemble high-order residual vector. Since A_ij is skew-symmetric, the 
          ! contribution of A_ji to node j is identical to that of A_ij to node i
          res(:,i) = res(:,i)+rA_ij
          res(:,j) = res(:,j)+rA_ij
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u(:,i), u(:,j), XDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(1)*Lbd_ij
          kb_ij = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
          
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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)

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
          
          ! Compute averaged coefficient
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u(:,i), u(:,j), XDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  C_ij(1)*Lbd_ij
          kb_ij  = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
          CALL fcb_getCharacteristics(u(:,i), u(:,j), YDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(2)*Lbd_ij
          kb_ij = (C_ij(2)+Cy(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
          
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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)

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
          
          ! Compute averaged coefficients
          C_ij(2)=0.5_DP*(Cy(ji)-Cy(ij))

          ! Compute characteristic fluxes in Y-direction
          CALL fcb_getCharacteristics(u(:,i), u(:,j), YDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(2)*Lbd_ij
          kb_ij = (C_ij(2)+Cy(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
          CALL fcb_getCharacteristics(u(:,i), u(:,j), ZDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(3)*Lbd_ij
          kb_ij = (C_ij(3)+Cz(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij

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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)

      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1

          ! Compute averaged coefficients
          C_ij(3)=0.5_DP*(Cz(ji)-Cz(ij))

          ! Compute characteristic fluxes in Z-direction
          CALL fcb_getCharacteristics(u(:,i), u(:,j), ZDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(3)*Lbd_ij
          kb_ij = (C_ij(3)+Cz(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,F_ij,W_ij,Lbd_ij,ka_ij,kb_ij
      REAL(DP), DIMENSION(NDIM1D)    :: C_ij
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
                   
          ! Compute solution difference
          diff = u(:,j)-u(:,i)
          
          ! Compute averaged coefficients a_ij:=(c_ji-c_ij)/2
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute local Roe matrix and multiply it with the solution difference
          !
          !   rA_ij:=A_ij*(U_j-U_i),
          !
          ! where A_ij denotes the cumulated Roe matrix 
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Assemble high-order residual vector. Since A_ij is skew-symmetric, the 
          ! contribution of A_ji to node j is identical to that of A_ij to node i
          res(:,i) = res(:,i)+rA_ij
          res(:,j) = res(:,j)+rA_ij
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u(:,i), u(:,j), XDir1D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(1)*Lbd_ij
          kb_ij = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
          
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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)
      
      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); Ksep(j) = Ksep(j)-1; ji = Ksep(j)
          
          ! Compute averaged coefficient
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u(:,i), u(:,j), XDir1D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  C_ij(1)*Lbd_ij
          kb_ij  = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,F_ij,W_ij,Lbd_ij,ka_ij,kb_ij
      REAL(DP), DIMENSION(NDIM2D)    :: C_ij
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

          ! Compute solution difference
          diff = u(:,j)-u(:,i)
          
          ! Compute averaged coefficients a_ij:=(c_ji-c_ij)/2
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix and multiply it with the solution difference
          !
          !   rA_ij:=A_ij*(U_j-U_i),
          !
          ! where A_ij denotes the cumulated Roe matrix 
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Assemble high-order residual vector. Since A_ij is skew-symmetric, the 
          ! contribution of A_ji to node j is identical to that of A_ij to node i
          res(:,i) = res(:,i)+rA_ij
          res(:,j) = res(:,j)+rA_ij
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u(:,i), u(:,j), XDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(1)*Lbd_ij
          kb_ij = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
          
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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)
      
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
          
          ! Compute averaged coefficient
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u(:,i), u(:,j), XDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  C_ij(1)*Lbd_ij
          kb_ij  = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
          CALL fcb_getCharacteristics(u(:,i), u(:,j), YDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(2)*Lbd_ij
          kb_ij = (C_ij(2)+Cy(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
          
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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)
      
      ! Loop over all rows (forward)
      DO i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute averaged coefficients
          C_ij(2)=0.5_DP*(Cy(ji)-Cy(ij))

          ! Compute characteristic fluxes in Y-direction
          CALL fcb_getCharacteristics(u(:,i), u(:,j), YDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(2)*Lbd_ij
          kb_ij = (C_ij(2)+Cy(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
      REAL(DP), DIMENSION(NVAR)      :: diff,rA_ij,F_ij,W_ij,Lbd_ij,ka_ij,kb_ij
      REAL(DP), DIMENSION(NDIM3D)    :: C_ij
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
                   
          ! Compute solution difference
          diff = u(:,j)-u(:,i)
          
          ! Compute averaged coefficients a_ij:=(c_ji-c_ij)/2
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute local Roe matrix and multiply it with the solution difference
          !
          !   rA_ij:=A_ij*(U_j-U_i),
          !
          ! where A_ij denotes the cumulated Roe matrix 
          CALL fcb_getRoeMatrix(u(:,i), u(:,j), C_ij, A_ij)
          CALL DGEMV('n', NVAR, NVAR, dscale, A_ij, NVAR, diff, 1, 0._DP, rA_ij, 1)
          
          ! Assemble high-order residual vector. Since A_ij is skew-symmetric, the 
          ! contribution of A_ji to node j is identical to that of A_ij to node i
          res(:,i) = res(:,i)+rA_ij
          res(:,j) = res(:,j)+rA_ij
                   
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u(:,i), u(:,j), XDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(1)*Lbd_ij
          kb_ij = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
          
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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)

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
          
          ! Compute averaged coefficient
          C_ij(1) = 0.5_DP*(Cx(ji)-Cx(ij))
          C_ij(2) = 0.5_DP*(Cy(ji)-Cy(ij))
          
          ! Compute characteristic fluxes in X-direction
          CALL fcb_getCharacteristics(u(:,i), u(:,j), XDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  C_ij(1)*Lbd_ij
          kb_ij  = (C_ij(1)+Cx(ij))*Lbd_ij
          F_ij   = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
          CALL fcb_getCharacteristics(u(:,i), u(:,j), YDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(2)*Lbd_ij
          kb_ij = (C_ij(2)+Cy(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
          
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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)

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
          
          ! Compute averaged coefficients
          C_ij(2)=0.5_DP*(Cy(ji)-Cy(ij))

          ! Compute characteristic fluxes in Y-direction
          CALL fcb_getCharacteristics(u(:,i), u(:,j), YDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(2)*Lbd_ij
          kb_ij = (C_ij(2)+Cy(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
          CALL fcb_getCharacteristics(u(:,i), u(:,j), ZDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(3)*Lbd_ij
          kb_ij = (C_ij(3)+Cz(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij

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
      rp(:,1:NEQ) = afcstab_limit( pp(:,1:NEQ), qp(:,1:NEQ), 1._DP)
      rm(:,1:NEQ) = afcstab_limit(-pm(:,1:NEQ),-qm(:,1:NEQ), 1._DP)

      ! Loop over all rows (backward)
      DO i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        DO ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); Ksep(j) = Ksep(j)-1; ji = Ksep(j)

          ! Compute averaged coefficients
          C_ij(3)=0.5_DP*(Cz(ji)-Cz(ij))

          ! Compute characteristic fluxes in Z-direction
          CALL fcb_getCharacteristics(u(:,i), u(:,j), ZDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  C_ij(3)*Lbd_ij
          kb_ij = (C_ij(3)+Cz(ij))*Lbd_ij
          F_ij  = -MAX(0._DP, MIN(ABS(ka_ij)-kb_ij, 2._DP*ABS(ka_ij)))*W_ij
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
      fcb_getRoeMatrix, fcb_getDissipation, iassembly, hstep, bclear, rmatrixJ)

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
    
    ! type of assembly
    !   iassembly = 0 : diagonal + Galerkin
    !   iassembly = 1 :            Galerkin
    !   iassembly = 2 : diagonal + scalar dissipation
    !   iassembly = 3 :            scalar dissipation
    !   iassembly = 4 : diagonal + tensorial dissipation
    !   iassembly = 5 :            tensorial dissipation
    INTEGER, INTENT(IN)                            :: iassembly

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
    IF (ru%nblocks .EQ. 1    .AND.&
        rmatrixJ%ndiagblocks .EQ. 1) THEN
      CALL gfsys_buildDivJacobianScalar(RmatrixC, ru%RvectorBlock(1),&
          fcb_getRoeMatrix, fcb_getDissipation, iassembly, hstep, bclear,&
          rmatrixJ%RmatrixBlock(1,1))
      RETURN       
    END IF

    PRINT *, "Jacobian matrix for systems is not yet implemented!"
    STOP
  END SUBROUTINE gfsys_buildDivJacobianBlock

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsys_buildDivJacobianScalar(RmatrixC, ru,&
      fcb_getRoeMatrix, fcb_getDissipation, iassembly, hstep, bclear, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: RmatrixC
    
    ! solution vector
    TYPE(t_vectorScalar), INTENT(IN)               :: ru
    
    ! type of assembly
    !   iassembly = 0 : diagonal + Galerkin
    !   iassembly = 1 :            Galerkin
    !   iassembly = 2 : diagonal + scalar dissipation
    !   iassembly = 3 :            scalar dissipation
    !   iassembly = 4 : diagonal + tensorial dissipation
    !   iassembly = 5 :            tensorial dissipation
    INTEGER, INTENT(IN)                            :: iassembly

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


    PRINT *, "Jacobian matrix for systems is not yet implemented!"
    STOP
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
