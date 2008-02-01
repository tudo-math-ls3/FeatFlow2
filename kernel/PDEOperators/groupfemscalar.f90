!##############################################################################
!# ****************************************************************************
!# <name> groupfemscalar </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic routines for applying the group-finite
!# element formulation to scalar problems, i.e. conservation laws.
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
!# the algebraic flux correction (AFC) methodology proposed by Kuzmin, Möller
!# and Turek in a series of publications. As a starting point for scalar
!# conservation laws, the reader is referred to the book chapter
!#
!#     D. Kuzmin and M. Möller, "Algebraic flux correction I. Scalar conservation
!#     laws", In: D. Kuzmin et al. (eds), Flux-Corrected Transport: Principles, 
!#     Algorithms, and Applications, Springer, 2005, 155-206.
!#
!# A more detailed description of the algorithms is given in the comments of the
!# subroutine implementing the corresponding discretisation schemes. All methods
!# are based on the stabilisation structure t_afcstab which is defined in the 
!# underlying module "afcstabilisation". The initialisation as a scalar stabilisation
!# structure is done by the routine gfsc_initStabilisation. 
!# 
!# There are three types of routines. The gfsc_buildXXXOperator routines can be
!# used to assemble the discrete convection or diffusion operators resulting
!# from the standard Galerkin finite element discretisation plus some discretely
!# defined artificial diffusion. For convective terms, this technique is termed
!# "discrete upwinding", whereas for physical diffusion operators this approach
!# ensures that the discrete maximum principle holds. Consequently, the term
!# DMP is adopted.
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
!# 1.) gfsc_initStabilisation
!#     -> initialize the stabilisation structure
!#
!# 2.) gfsc_isMatrixCompatible
!#     -> checks wether a matrix and a stabilisation structure are compatible
!#
!# 3.) gfsc_isVectorCompatible
!#     -> checks wether a vector and a stabilisation structure are compatible
!#
!# 4.) gfsc_buildConvectionOperator = gfsc_buildConvOperatorScalar /
!#                                    gfsc_buildConvOperatorBlock
!#     -> assemble the convective part of the transport operator for 
!#        a scalar convection-diffusion-reaction equation
!#
!# 5.) gfsc_buildDiffusionOperator
!#     -> assemble the diffusive part of the transport operator for
!#        scalar convection-diffusion-reaction equation 
!#
!# 6.) gfsc_buildResidualFCT = gfsc_buildResScalarFCT /
!#                             gfsc_buildResBlockFCT
!#     -> assemble the residual vector for AFC stabilisation of FCT type
!#
!# 7.) gfsc_buildResidualTVD = gfsc_buildResScalarTVD /
!#                             gfsc_buildResBlockTVD
!#     -> assemble the residual vector for AFC stabilisation of TVD type
!#
!# 8.) gfsc_buildResidualSymm = gfsc_buildResScalarSymm /
!#                                    gfsc_buildResBlockSymm
!#     -> assemble the residual vector for stabilisation by means of
!#        symmetric flux limiting for diffusion operators
!#
!# 9.) gfsc_buildConvectionJacobian = gfsc_buildConvJacobianScalar /
!#                                    gfsc_buildConvJacobianBlock
!#     -> assemble the Jacobian matrix for the convective part of
!#        the transport operator for a scalar convection equation
!#
!# 10.) gfsc_buildStabLinearJacobian_FCT = gfsc_buildStabJacLinearScalar_FCT /
!#                                         gfsc_buildStabJacLinearBlock_FCT
!#      -> assemble the Jacobian matrix for the stabilisation part of FCT type;
!#         the velocity is assumed to be linear
!#
!# 11.) gfsc_buildStabLinearJacobian_GPTVD = gfsc_buildStabJacLinearScalar_GPTVD /
!#                                           gfsc_buildStabJacLinearBlock_GPTVD
!#      -> assemble the Jacobian matrix for the stabilisation part of TVD type
!#         and/or for the general purpose limiter; the velocity is assumed to be linear
!#
!# 12.) gfsc_buildStabJacobian_FCT = gfsc_buildStabJacobianScalar_FCT /
!#                                   gfsc_buildStabJacobianBlock_FCT
!#      -> assemble the Jacobian matrix for the stabilisation part of FCT type;
!#         the velocity can be arbitrary
!#
!# 13.) gfsc_buildStabJacobian_GPTVD = gfsc_buildStabJacobianScalar_GPTVD /
!#                                     gfsc_buildStabJacobianBlock_GPTVD
!#      -> assemble the Jacobian matrix for the stabilisation part of TVD type
!#         and/or for the general purpose limiter; the velocity can be arbitrary
!#
!# 14.) gfsc_buildStabJacobian_Symmetric = gfsc_buildStabJacobianScalar_Symm /
!#                                         gfsc_buildStabJacobianBlock_Symm
!#      -> assemble the Jacobian matrix for the stabilisation part of symmetric
!#         flux limiting for diffusion operators
!#
!# </purpose>
!##############################################################################

MODULE groupfemscalar

  USE afcstabilisation
  USE fsystem
  USE genoutput
  USE linearalgebra
  USE linearsystemblock
  USE linearsystemscalar
  USE storage

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: gfsc_initStabilisation
  PUBLIC :: gfsc_isMatrixCompatible
  PUBLIC :: gfsc_isVectorCompatible
  PUBLIC :: gfsc_buildConvectionOperator
  PUBLIC :: gfsc_buildDiffusionOperator
  PUBLIC :: gfsc_buildResidualFCT
  PUBLIC :: gfsc_buildResidualTVD
  PUBLIC :: gfsc_buildResidualSymm
  PUBLIC :: gfsc_buildConvectionJacobian
  PUBLIC :: gfsc_buildStabLinearJacobian_FCT
  PUBLIC :: gfsc_buildStabLinearJacobian_GPTVD
  PUBLIC :: gfsc_buildStabJacobian_FCT
  PUBLIC :: gfsc_buildStabJacobian_GPTVD
  PUBLIC :: gfsc_buildStabJacobian_Symmetric

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

  INTERFACE gfsc_buildConvectionOperator
    MODULE PROCEDURE gfsc_buildConvOperatorScalar
    MODULE PROCEDURE gfsc_buildConvOperatorBlock
  END INTERFACE

  INTERFACE gfsc_buildResidualFCT
    MODULE PROCEDURE gfsc_buildResScalarFCT
    MODULE PROCEDURE gfsc_buildResBlockFCT
  END INTERFACE

  INTERFACE gfsc_buildResidualTVD
    MODULE PROCEDURE gfsc_buildResScalarTVD
    MODULE PROCEDURE gfsc_buildResBlockTVD
  END INTERFACE

  INTERFACE gfsc_buildResidualSymm
    MODULE PROCEDURE gfsc_buildResScalarSymm
    MODULE PROCEDURE gfsc_buildResBlockSymm
  END INTERFACE

  INTERFACE gfsc_buildConvectionJacobian
    MODULE PROCEDURE gfsc_buildConvJacobianScalar
    MODULE PROCEDURE gfsc_buildConvJacobianBlock
  END INTERFACE

  INTERFACE gfsc_buildStabLinearJacobian_FCT
    MODULE PROCEDURE gfsc_buildStabJacLinearScalar_FCT
    MODULE PROCEDURE gfsc_buildStabJacLinearBlock_FCT
  END INTERFACE

  INTERFACE gfsc_buildStabLinearJacobian_GPTVD
    MODULE PROCEDURE gfsc_buildStabJacLinearScalar_GPTVD
    MODULE PROCEDURE gfsc_buildStabJacLinearBlock_GPTVD
  END INTERFACE

  INTERFACE gfsc_buildStabJacobian_FCT
    MODULE PROCEDURE gfsc_buildStabJacobianScalar_FCT
    MODULE PROCEDURE gfsc_buildStabJacobianBlock_FCT
  END INTERFACE

  INTERFACE gfsc_buildStabJacobian_GPTVD
    MODULE PROCEDURE gfsc_buildStabJacobianScalar_GPTVD
    MODULE PROCEDURE gfsc_buildStabJacobianBlock_GPTVD
  END INTERFACE

  INTERFACE gfsc_buildStabJacobian_Symmetric
    MODULE PROCEDURE gfsc_buildStabJacobianScalar_Symm
    MODULE PROCEDURE gfsc_buildStabJacobianBlock_Symm
  END INTERFACE
  
  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

CONTAINS

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_initStabilisation(rmatrixTemplate, rafcstab)

!<description>
    ! This subroutine initialises the discrete stabilisation structure
    ! for use as a scalar stabilisation. The template matrix is used
    ! to determine the number of equations and the number of edges.
    ! Note that there exists an edge between two nodes i and j of their
    ! basis functions have overlapping support. For linear finite elements
    ! the number of physical edges coincides with the number of edges used
    ! in the stabilisation structure. In general the number of edges equals
    ! (NA-NEQ)/2, that is, it is the number of nonzero matrix entries not
    ! counting the diagonal divided by two.
!</description>

!<input>
    ! The template matrix
    TYPE(t_matrixScalar), INTENT(IN) :: rmatrixTemplate
!</input>

!<inputoutput>
    ! The stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)   :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER :: i
    INTEGER(I32), DIMENSION(2) :: Isize
    
    ! Set atomic data
    rafcstab%NVAR  = rmatrixTemplate%NVAR
    rafcstab%NEQ   = rmatrixTemplate%NEQ
    rafcstab%NEDGE = INT(0.5*(rmatrixTemplate%NA-rmatrixTemplate%NEQ),I32)

    ! What kind of stabilisation are we?
    SELECT CASE(rafcstab%ctypeAFCstabilisation)
      
    CASE (AFCSTAB_GALERKIN, AFCSTAB_UPWIND, AFCSTAB_DMP)
      ! do nothing

      
    CASE (AFCSTAB_FEMFCT, AFCSTAB_FEMGP, AFCSTAB_FEMTVD)
      
      ! Handle for IsuperdiagonalEdgesIdx
      IF (rafcstab%h_IsuperdiagonalEdgesIdx .NE. ST_NOHANDLE)&
          CALL storage_free(rafcstab%h_IsuperdiagonalEdgesIdx)
      CALL storage_new('gfsc_initStabilisation', 'IsuperdiagonalEdgesIdx',&
          rafcstab%NEQ+1, ST_INT, rafcstab%h_IsuperdiagonalEdgesIdx, ST_NEWBLOCK_NOINIT)
      
      ! Handle for IverticesAtEdge
      Isize = (/4, rafcstab%NEDGE/)
      IF (rafcstab%h_IverticesAtEdge .NE. ST_NOHANDLE)&
          CALL storage_free(rafcstab%h_IverticesAtEdge)
      CALL storage_new('gfsc_initStabilisation', 'IverticesAtEdge',&
          Isize, ST_INT, rafcstab%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)

      ! Handle for DcoefficientsAtEdge
      Isize = (/3, rafcstab%NEDGE/)
      IF (rafcstab%h_DcoefficientsAtEdge .NE. ST_NOHANDLE)&
          CALL storage_free(rafcstab%h_DcoefficientsAtEdge)
      CALL storage_new('gfsc_initStabilisation', 'DcoefficientsAtEdge',&
          Isize, ST_DOUBLE, rafcstab%h_DcoefficientsAtEdge, ST_NEWBLOCK_NOINIT)

      ! Nodal vectors
      ALLOCATE(rafcstab%RnodalVectors(6))
      DO i = 1, 6
        CALL lsyssc_createVector(rafcstab%RnodalVectors(i),&
            rafcstab%NEQ, .FALSE., ST_DOUBLE)
      END DO

      ! Edgewise vectors
      ALLOCATE(rafcstab%RedgeVectors(2))
      DO i = 1, 2
        CALL lsyssc_createVector(rafcstab%RedgeVectors(i),&
            rafcstab%NEDGE, .FALSE., ST_DOUBLE)
      END DO


    CASE (AFCSTAB_SYMMETRIC)

      ! Handle for IsuperdiagonalEdgesIdx
      IF (rafcstab%h_IsuperdiagonalEdgesIdx .NE. ST_NOHANDLE)&
          CALL storage_free(rafcstab%h_IsuperdiagonalEdgesIdx)
      CALL storage_new('gfsc_initStabilisation', 'IsuperdiagonalEdgesIdx',&
          rafcstab%NEQ+1, ST_INT, rafcstab%h_IsuperdiagonalEdgesIdx, ST_NEWBLOCK_NOINIT)
      
      ! Handle for IverticesAtEdge
      Isize = (/2, rafcstab%NEDGE/)
      IF (rafcstab%h_IverticesAtEdge .NE. ST_NOHANDLE)&
          CALL storage_free(rafcstab%h_IverticesAtEdge)
      CALL storage_new('gfsc_initStabilisation', 'IverticesAtEdge',&
          Isize, ST_INT, rafcstab%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)

      ! Handle for DcoefficientsAtEdge
      Isize = (/2, rafcstab%NEDGE/)
      IF (rafcstab%h_DcoefficientsAtEdge .NE. ST_NOHANDLE)&
          CALL storage_free(rafcstab%h_DcoefficientsAtEdge)
      CALL storage_new('gfsc_initStabilisation', 'DcoefficientsAtEdge',&
          Isize, ST_DOUBLE, rafcstab%h_DcoefficientsAtEdge, ST_NEWBLOCK_NOINIT)

      ! Nodal vectors
      ALLOCATE(rafcstab%RnodalVectors(6))
      DO i = 1, 6
        CALL lsyssc_createVector(rafcstab%RnodalVectors(i),&
            rafcstab%NEQ, .FALSE., ST_DOUBLE)
      END DO
      
      ! Edgewise vectors
      ALLOCATE(rafcstab%RedgeVectors(1))
      DO i = 1, 1
        CALL lsyssc_createVector(rafcstab%RedgeVectors(i),&
            rafcstab%NEDGE, .FALSE., ST_DOUBLE)
      END DO

    CASE DEFAULT
      CALL output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_initStabilisation')
      CALL sys_halt()
    END SELECT

    ! Set specifier
    rafcstab%iSpec = AFCSTAB_INITIALISED
  END SUBROUTINE gfsc_initStabilisation

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_isMatrixCompatible(rafcstab, rmatrix, bcompatible)

!<description>
    ! This subroutine checks if a scalar matrix and a discrete 
    ! stabilisation structure are compatible to each other, 
    ! i.e. if they share the same structure, size and so on.
!</description>

!<input>
    ! The scalar matrix
    TYPE(t_matrixScalar), INTENT(IN) :: rmatrix

    ! The stabilisation structure
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
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_isMatrixCompatible')
        CALL sys_halt()
      END IF
    END IF
  END SUBROUTINE gfsc_isMatrixCompatible

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_isVectorCompatible(rafcstab, rvector, bcompatible)

!<description>
    ! This subroutine checks if a vector and a stabilisation
    ! structure are compatible to each other, i.e., share the
    ! same structure, size and so on.
!</description>

!<input>
    ! The scalar vector
    TYPE(t_vectorScalar), INTENT(IN) :: rvector

    ! Teh stabilisation structure
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
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_isVectorCompatible')
        CALL sys_halt()
      END IF
    END IF
  END SUBROUTINE gfsc_isVectorCompatible

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_buildConvOperatorBlock(RmatrixC, ru,&
      fcb_getVelocity, bZeroRowsum, bStabilize, bclear, rmatrixL, rafcstab)
    
!<description>
    ! This subroutine assembles the discrete transport operator which results
    ! from the group finite element formulation of the continuous problem
    !
    !     $$ \nabla\cdot({\bf v}u)$$
    !
    ! Note that this routine serves as a wrapper for block vectors. If there
    ! is only one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
!</description>

!<input>
    ! The array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: RmatrixC

    ! The solution vector
    ! Note that this vector is only required for nonlinear
    ! problems which require the evaluation of the velocity
    TYPE(t_vectorBlock), INTENT(IN)                :: ru
    
    ! Switch for row-sum property
    ! TRUE  : assume zero row-sum, that is, the diegonal coefficient
    !         is computed as the sum of negative off-diagonal coefficients
    ! FALSE : compute the diagonal coefficient by standard finite elements
    LOGICAL, INTENT(IN)                            :: bZeroRowsum

    ! Switch for stabilisation
    ! TRUE  : perform stabilisation
    ! FALSE : perform no stabilisation
    LOGICAL, INTENT(IN)                            :: bStabilize

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    LOGICAL, INTENT(IN)                            :: bclear

    ! callback functions to compute velocity
    INCLUDE 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! The transport operator
    TYPE(t_matrixScalar), INTENT(INOUT)            :: rmatrixL
    
    ! OPTIONAL: the stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT), OPTIONAL       :: rafcstab
!</inputoutput>
!</subroutine>

    ! Check if block vector contains exactly one block
    IF (ru%nblocks .NE. 1) THEN

      CALL output_line('Solution vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvOperatorBlock')
      CALL sys_halt()

    ELSE
      
      CALL gfsc_buildConvOperatorScalar(RmatrixC, ru%RvectorBlock(1),&
          fcb_getVelocity, bZeroRowsum, bStabilize, bclear, rmatrixL, rafcstab)

    END IF
  END SUBROUTINE gfsc_buildConvOperatorBlock
  
  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_buildConvOperatorScalar(RmatrixC, ru, fcb_getVelocity,&
      bZeroRowsum, bStabilize, bclear, rmatrixL, rafcstab)

!<description>
    ! This subroutine assembles the discrete transport operator which results
    ! from the group finite element formulation of the continuous problem
    !
    !     $$ \nabla\cdot({\bf v}u)$$
    !
    ! This routine can be used to apply the following discretisations:
    !
    ! (1) the standard Galerkin finite element method
    !     which will be referred to as high-order approximation
    !
    ! (2) the discrete upwinding operator which results from the conservative
    !     elimination of negative off-diagonal entries from the Galerkin
    !     operator. This technique is for instance described in the reference:
    !
    !     D. Kuzmin and M. Möller, "Algebraic flux correction I. Scalar
    !     conservation laws", In: D. Kuzmin et al. (eds), Flux-Corrected
    !     Transport: Principles,  Algorithms, and Applications,
    !     Springer, 2005, 155-206.
    !
    ! (3) In addition to (2), discrete upwinding is performed and auxiliary
    !     data required for flux limiting is generated in the optional 
    !     stabilisation structure rafcstab. For symmetric flux limiters,
    !     the artificil diffusio coefficient d_ij and the entries of the
    !     low-order operator l_ij=k_ij+d_ij are stored as well as the node
    !     numbers i and j and the matrix positions ij/ji of the edge (i,j).
    !     For upwind-biased flux limiters the same data are stored but the
    !     following orientation convention is applied. Node i is located 
    !     "upwind" and corresponds to the column number whose entry has been
    !     eliminated.
!</description>

!<input>
    ! The array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: rmatrixC

    ! The solution vector
    ! Note that this vector is only required for nonlinear
    ! problems which require the evaluation of the velocity
    TYPE(t_vectorScalar), INTENT(IN)               :: ru

    ! Switch for row-sum property
    ! TRUE  : assume zero row-sum, that is, the diegonal coefficient
    !         is computed as the sum of negative off-diagonal coefficients
    ! FALSE : compute the diagonal coefficient by standard finite elements
    LOGICAL, INTENT(IN)                            :: bZeroRowsum

    ! Switch for stabilisation
    ! TRUE  : perform stabilisation
    ! FALSE : perform no stabilisation
    LOGICAL, INTENT(IN)                            :: bStabilize

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    LOGICAL, INTENT(IN)                            :: bclear

    ! callback functions to compute velocity
    INCLUDE 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! The transport operator
    TYPE(t_matrixScalar), INTENT(INOUT)            :: rmatrixL
    
    ! OPTIONAL: the stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT), OPTIONAL       :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_MATIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER   :: p_IsuperdiagonalEdgesIdx
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_Kld,p_Ksep
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_Kdiagonal
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER   :: p_Kcol
    REAL(DP), DIMENSION(:,:), POINTER             :: p_DcoefficientsAtEdge
    REAL(DP), DIMENSION(:), POINTER               :: p_Cx,p_Cy,p_Cz,p_L,p_u
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

    ! Set pointer
    CALL lsyssc_getbase_Kld   (rmatrixL, p_Kld)
    CALL lsyssc_getbase_Kcol  (rmatrixL, p_Kcol)
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
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvOperatorScalar')
      CALL sys_halt()
    END SELECT
    

    ! What kind of matrix are we?
    SELECT CASE(rmatrixL%cmatrixFormat)
    CASE(LSYSSC_MATRIX7)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      CALL storage_copy(rmatrixL%h_Kld, h_Ksep)
      CALL storage_getbase_int(h_Ksep, p_Ksep, rmatrixL%NEQ+1)

      ! Do we have a stabilisation structure?
      IF (PRESENT(rafcstab)) THEN

        ! Check if stabilisation has been initialised
        IF (IAND(rafcstab%iSpec, AFCSTAB_INITIALISED) .EQ. 0) THEN
          CALL output_line('Stabilisation has not been initialised',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvOperatorScalar')
          CALL sys_halt()
        END IF
        
        ! Check if matrix/vector and stabilisation
        ! structure are compatible to each other
        CALL gfsc_isMatrixCompatible(rafcstab, rmatrixL)
        CALL gfsc_isVectorCompatible(rafcstab, ru)
        
        ! What kind of stabilisation are we?
        SELECT CASE (rafcstab%ctypeAFCstabilisation)
        CASE (AFCSTAB_UPWIND)
          ! Perform discrete upwinding without generating edge structure

          IF (bZeroRowsum) THEN
            
            SELECT CASE(ndim)
            CASE (NDIM1D)
              CALL doUpwindZRSMat7_1D(p_Kld, p_Kcol, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_u, p_L)
            CASE (NDIM2D)
              CALL doUpwindZRSMat7_2D(p_Kld, p_Kcol, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_u, p_L)
            CASE (NDIM3D)
              CALL doUpwindZRSMat7_3D(p_Kld, p_Kcol, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_L)
            END SELECT

          ELSE

            SELECT CASE(ndim)
            CASE (NDIM1D)
              CALL doUpwindMat7_1D(p_Kld, p_Kcol, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_u, p_L)
            CASE (NDIM2D)
              CALL doUpwindMat7_2D(p_Kld, p_Kcol, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_u, p_L)
            CASE (NDIM3D)
              CALL doUpwindMat7_3D(p_Kld, p_Kcol, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_L)
            END SELECT

          END IF
          

        CASE (AFCSTAB_FEMFCT)
          ! Set additional pointers
          CALL afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
          CALL afcstab_getbase_IverticesAtEdge(rafcstab,  p_IverticesAtEdge)
          CALL afcstab_getbase_DcoeffsAtEdge(rafcstab,    p_DcoefficientsAtEdge)
          
          ! Adopt no orientation convention and generate edge structure

          IF (bZeroRowsum) THEN
            
            SELECT CASE(ndim)
            CASE (NDIM1D)
              CALL doUpwindZRSMat7_AFC_1D(p_Kld, p_Kcol, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            CASE (NDIM2D)
              CALL doUpwindZRSMat7_AFC_2D(p_Kld, p_Kcol, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            CASE (NDIM3D)
              CALL doUpwindZRSMat7_AFC_3D(p_Kld, p_Kcol, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            END SELECT

          ELSE

            SELECT CASE(ndim)
            CASE (NDIM1D)
              CALL doUpwindMat7_AFC_1D(p_Kld, p_Kcol, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            CASE (NDIM2D)
              CALL doUpwindMat7_AFC_2D(p_Kld, p_Kcol, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            CASE (NDIM3D)
              CALL doUpwindMat7_AFC_3D(p_Kld, p_Kcol, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            END SELECT

          END IF

          
          ! Set state of stabilisation
          rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)
          rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_EDGEVALUES)

          
        CASE (AFCSTAB_FEMTVD, AFCSTAB_FEMGP)
          ! Set additional pointers
          CALL afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
          CALL afcstab_getbase_IverticesAtEdge(rafcstab,  p_IverticesAtEdge)
          CALL afcstab_getbase_DcoeffsAtEdge(rafcstab,    p_DcoefficientsAtEdge)

          ! Adopt orientation convention IJ, such that L_ij < L_ji and
          ! generate edge structure for the flux limiter

          IF (bZeroRowsum) THEN
            
            SELECT CASE(ndim)
            CASE (NDIM1D)
              CALL doUpwindZRSMat7_ordAFC_1D(p_Kld, p_Kcol, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            CASE (NDIM2D)
              CALL doUpwindZRSMat7_ordAFC_2D(p_Kld, p_Kcol, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            CASE (NDIM3D)
              CALL doUpwindZRSMat7_ordAFC_3D(p_Kld, p_Kcol, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            END SELECT

          ELSE

            SELECT CASE(ndim)
            CASE (NDIM1D)
              CALL doUpwindMat7_ordAFC_1D(p_Kld, p_Kcol, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            CASE (NDIM2D)
              CALL doUpwindMat7_ordAFC_2D(p_Kld, p_Kcol, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            CASE (NDIM3D)
              CALL doUpwindMat7_ordAFC_3D(p_Kld, p_Kcol, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            END SELECT

          END IF


          ! Set state of stabilisation
          rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)
          rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_EDGEVALUES)
          rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_EDGEORIENTATION)

          
        CASE DEFAULT
          CALL output_line('Invalid type of AFC stabilisation!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvOperatorScalar')
          CALL sys_halt()
        END SELECT
      
      ELSEIF (bStabilize) THEN
        
        ! Perform discrete upwinding without generating edge structure
        IF (bZeroRowsum) THEN
          
          SELECT CASE(ndim)
          CASE (NDIM1D)
            CALL doUpwindZRSMat7_1D(p_Kld, p_Kcol, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_u, p_L)
          CASE (NDIM2D)
            CALL doUpwindZRSMat7_2D(p_Kld, p_Kcol, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_Cy, p_u, p_L)
          CASE (NDIM3D)
            CALL doUpwindZRSMat7_3D(p_Kld, p_Kcol, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_L)
          END SELECT
          
        ELSE

          SELECT CASE(ndim)
          CASE (NDIM1D)
            CALL doUpwindMat7_1D(p_Kld, p_Kcol, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_u, p_L)
          CASE (NDIM2D)
            CALL doUpwindMat7_2D(p_Kld, p_Kcol, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_Cy, p_u, p_L)
          CASE (NDIM3D)
            CALL doUpwindMat7_3D(p_Kld, p_Kcol, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_L)
          END SELECT

        END IF

      ELSE
        
        ! Apply standard Galerkin discretisation
        IF (bZeroRowsum) THEN
          
          SELECT CASE(ndim)
          CASE (NDIM1D)
            CALL doGalerkinZRSMat7_1D(p_Kld, p_Kcol, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_u, p_L)
          CASE (NDIM2D)
            CALL doGalerkinZRSMat7_2D(p_Kld, p_Kcol, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_Cy, p_u, p_L)
          CASE (NDIM3D)
            CALL doGalerkinZRSMat7_3D(p_Kld, p_Kcol, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_L)
          END SELECT

        ELSE

          SELECT CASE(ndim)
          CASE (NDIM1D)
            CALL doGalerkinMat7_1D(p_Kld, p_Kcol, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_u, p_L)
          CASE (NDIM2D)
            CALL doGalerkinMat7_2D(p_Kld, p_Kcol, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_Cy, p_u, p_L)
          CASE (NDIM3D)
            CALL doGalerkinMat7_3D(p_Kld, p_Kcol, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_L)
          END SELECT

        END IF
      END IF

      ! Release diagonal separator
      CALL storage_free(h_Ksep)


    CASE(LSYSSC_MATRIX9)

      ! Set pointers
      CALL lsyssc_getbase_Kdiagonal(rmatrixL, p_Kdiagonal)

      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      CALL storage_copy(rmatrixL%h_Kld, h_Ksep)
      CALL storage_getbase_int(h_Ksep, p_Ksep, rmatrixL%NEQ+1)

      ! Do we have a stabilisation structure?
      IF (PRESENT(rafcstab)) THEN
        
        ! Check if stabilisation has been initialised
        IF (IAND(rafcstab%iSpec, AFCSTAB_INITIALISED) .EQ. 0) THEN
          CALL output_line('Stabilisation has not been initialised',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvOperatorScalar')
          CALL sys_halt()
        END IF
        
        ! Check if matrix/vector and stabilisation
        ! structure are compatible to each other
        CALL gfsc_isMatrixCompatible(rafcstab, rmatrixL)
        CALL gfsc_isVectorCompatible(rafcstab, ru)

        ! What kind of stabilisation are we?
        SELECT CASE (rafcstab%ctypeAFCstabilisation)
        CASE (AFCSTAB_UPWIND)

          ! Perform discrete upwinding without generating edge structure
          IF (bZeroRowsum) THEN
            
            SELECT CASE(ndim)
            CASE (NDIM1D)
              CALL doUpwindZRSMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_u, p_L)
            CASE (NDIM2D)
              CALL doUpwindZRSMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_u, p_L)
            CASE (NDIM3D)
              CALL doUpwindZRSMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_L)
            END SELECT

          ELSE

            SELECT CASE(ndim)
            CASE (NDIM1D)
              CALL doUpwindMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_u, p_L)
            CASE (NDIM2D)
              CALL doUpwindMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_u, p_L)
            CASE (NDIM3D)
              CALL doUpwindMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_L)
            END SELECT

          END IF
          

        CASE (AFCSTAB_FEMFCT)
          ! Set additional pointers
          CALL afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
          CALL afcstab_getbase_IverticesAtEdge(rafcstab,  p_IverticesAtEdge)
          CALL afcstab_getbase_DcoeffsAtEdge(rafcstab,    p_DcoefficientsAtEdge)
          
          ! Adopt no orientation convention and generate edge structure
          IF (bZeroRowsum) THEN
            
            SELECT CASE(ndim)
            CASE (NDIM1D)
              CALL doUpwindZRSMat9_AFC_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            CASE (NDIM2D)
              CALL doUpwindZRSMat9_AFC_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            CASE (NDIM3D)
              CALL doUpwindZRSMat9_AFC_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            END SELECT

          ELSE

            SELECT CASE(ndim)
            CASE (NDIM1D)
              CALL doUpwindMat9_AFC_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            CASE (NDIM2D)
              CALL doUpwindMat9_AFC_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            CASE (NDIM3D)
              CALL doUpwindMat9_AFC_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            END SELECT

          END IF
          

          ! Set state of stabilisation
          rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)
          rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_EDGEVALUES)

          
        CASE (AFCSTAB_FEMTVD, AFCSTAB_FEMGP)
          ! Set additional pointers
          CALL afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
          CALL afcstab_getbase_IverticesAtEdge(rafcstab,  p_IverticesAtEdge)
          CALL afcstab_getbase_DcoeffsAtEdge(rafcstab,   p_DcoefficientsAtEdge)
          
          ! Adopt orientation convention IJ, such that L_ij < L_ji and
          ! generate edge structure for the flux limiter
          IF (bZeroRowsum) THEN
            
            SELECT CASE(ndim)
            CASE (NDIM1D)
              CALL doUpwindZRSMat9_ordAFC_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            CASE (NDIM2D)
              CALL doUpwindZRSMat9_ordAFC_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            CASE (NDIM3D)
              CALL doUpwindZRSMat9_ordAFC_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            END SELECT

          ELSE

            SELECT CASE(ndim)
            CASE (NDIM1D)
              CALL doUpwindMat9_ordAFC_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            CASE (NDIM2D)
              CALL doUpwindMat9_ordAFC_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            CASE (NDIM3D)
              CALL doUpwindMat9_ordAFC_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                  rmatrixL%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_L, p_IsuperdiagonalEdgesIdx,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            END SELECT

          END IF
          
          
          ! Set state of stabilisation
          rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)
          rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_EDGEVALUES)
          rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_EDGEORIENTATION)

          
        CASE DEFAULT
          CALL output_line('Invalid type of AFC stabilisation!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvOperatorScalar')
          CALL sys_halt()
        END SELECT
      
      ELSEIF (bStabilize) THEN
        
        ! Perform discrete upwinding without generating edge structure
        IF (bZeroRowsum) THEN
          
          SELECT CASE(ndim)
          CASE (NDIM1D)
            CALL doUpwindZRSMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_u, p_L)
          CASE (NDIM2D)
            CALL doUpwindZRSMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_Cy, p_u, p_L)
          CASE (NDIM3D)
            CALL doUpwindZRSMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_L)
          END SELECT

        ELSE

          SELECT CASE(ndim)
          CASE (NDIM1D)
            CALL doUpwindMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_u, p_L)
          CASE (NDIM2D)
            CALL doUpwindMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_Cy, p_u, p_L)
          CASE (NDIM3D)
            CALL doUpwindMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_L)
          END SELECT

        END IF
        
      ELSE
        
        ! Apply standard Galerkin discretisation
        IF (bZeroRowsum) THEN
          
          SELECT CASE(ndim)
          CASE (NDIM1D)
            CALL doGalerkinZRSMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_u, p_L)
          CASE (NDIM2D)
            CALL doGalerkinZRSMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_Cy, p_u, p_L)
          CASE (NDIM3D)
            CALL doGalerkinZRSMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_L)
          END SELECT

        ELSE

          SELECT CASE(ndim)
          CASE (NDIM1D)
            CALL doGalerkinMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_u, p_L)
          CASE (NDIM2D)
            CALL doGalerkinMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_Cy, p_u, p_L)
          CASE (NDIM3D)
            CALL doGalerkinMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                rmatrixL%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_L)
          END SELECT

        END IF

      END IF

      ! Release diagonal separator
      CALL storage_free(h_Ksep)


    CASE DEFAULT
      CALL output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvOperatorScalar')
      CALL sys_halt()
    END SELECT

  CONTAINS

    ! Here, the working routine follow
    
    !**************************************************************
    ! Assemble high-order Galerkin operator K in 1D
    ! and assume zero row-sums.
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doGalerkinZRSMat7_1D(Kld, Kcol, Ksep, NEQ, Cx, u, K)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: K
      
      REAL(DP), DIMENSION(NDIM1D) :: v_ij,v_ji
      REAL(DP)                    :: k_ij,k_ji
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
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)
          k_ji = -v_ji(1)*Cx(ji)
          
          ! Assemble the global operator
          K(ij) = K(ij)+k_ij; K(ii) = K(ii)-k_ij
          K(ji) = K(ji)+k_ji; K(jj) = K(jj)-k_ji
        END DO
      END DO
    END SUBROUTINE doGalerkinZRSMat7_1D

    
    !**************************************************************
    ! Assemble high-order Galerkin operator K in 1D.
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doGalerkinMat7_1D(Kld, Kcol, Ksep, NEQ, Cx, u, K)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: K
      
      REAL(DP), DIMENSION(NDIM1D) :: v_ij,v_ji
      REAL(DP)                    :: k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji
      INTEGER(PREC_VECIDX)        :: i,j
            
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)

        ! Update the diagonal coefficient
        K(ii) = K(ii)+k_ii

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)
          k_ji = -v_ji(1)*Cx(ji)
          
          ! Assemble the global operator
          K(ij) = K(ij)+k_ij; K(ji) = K(ji)+k_ji
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7_1D

    
    !**************************************************************
    ! Assemble high-order Galerkin operator K in 2D
    ! and assume zero row-sums.
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doGalerkinZRSMat7_2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, K)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: K
      
      REAL(DP), DIMENSION(NDIM2D) :: v_ij,v_ji
      REAL(DP)                    :: k_ij,k_ji
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
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          
          ! Assemble the global operator
          K(ij) = K(ij)+k_ij; K(ii) = K(ii)-k_ij
          K(ji) = K(ji)+k_ji; K(jj) = K(jj)-k_ji
        END DO
      END DO
    END SUBROUTINE doGalerkinZRSMat7_2D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 2D.
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doGalerkinMat7_2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, K)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: K
      
      REAL(DP), DIMENSION(NDIM2D) :: v_ij,v_ji
      REAL(DP)                    :: k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji
      INTEGER(PREC_VECIDX)        :: i,j
            
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)-v_ij(2)*Cy(ii)

        ! Update the diagonal coefficient
        K(ii) = K(ii)+k_ii
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          
          ! Assemble the global operator
          K(ij) = K(ij)+k_ij; K(ji) = K(ji)+k_ji
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7_2D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 3D
    ! and assume zero row-sums
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doGalerkinZRSMat7_3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, K)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: K
      
      REAL(DP), DIMENSION(NDIM3D) :: v_ij,v_ji
      REAL(DP)                    :: k_ij,k_ji
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
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          
          ! Assemble the global operator
          K(ij) = K(ij)+k_ij; K(ii) = K(ii)-k_ij
          K(ji) = K(ji)+k_ji; K(jj) = K(jj)-k_ji
        END DO
      END DO
    END SUBROUTINE doGalerkinZRSMat7_3D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 3D.
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doGalerkinMat7_3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, K)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: K
      
      REAL(DP), DIMENSION(NDIM3D) :: v_ij,v_ji
      REAL(DP)                    :: k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji
      INTEGER(PREC_VECIDX)        :: i,j
            
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)-v_ij(2)*Cy(ii)-v_ij(3)*Cz(ii)

        ! Update the diagonal coefficient
        K(ii) = K(ii)+k_ii

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          
          ! Assemble the global operator
          K(ij) = K(ij)+k_ij; K(ji) = K(ji)+k_ji
        END DO
      END DO
    END SUBROUTINE doGalerkinMat7_3D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 1D
    ! and assume zero row-sums
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doGalerkinZRSMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, K)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: K
      
      REAL(DP), DIMENSION(NDIM1D) :: v_ij,v_ji
      REAL(DP)                    :: k_ij,k_ji
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
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)
          k_ji = -v_ji(1)*Cx(ji)
          
          ! Assemble the global operator
          K(ij) = K(ij)+k_ij; K(ii) = K(ii)-k_ij
          K(ji) = K(ji)+k_ji; K(jj) = K(jj)-k_ji
        END DO
      END DO
    END SUBROUTINE doGalerkinZRSMat9_1D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 1D.
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doGalerkinMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, K)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: K
      
      REAL(DP), DIMENSION(NDIM1D) :: v_ij,v_ji
      REAL(DP)                    :: k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji
      INTEGER(PREC_VECIDX)        :: i,j
            
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)
        
        ! Update the diagonal coefficient
        K(ii) = K(ii)+k_ii
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)
          k_ji = -v_ji(1)*Cx(ji)
          
          ! Assemble the global operator
          K(ij) = K(ij)+k_ij; K(ji) = K(ji)+k_ji
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9_1D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 2D
    ! and assume zero row-sums
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doGalerkinZRSMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, K)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: K
      
      REAL(DP), DIMENSION(NDIM2D) :: v_ij,v_ji
      REAL(DP)                    :: k_ij,k_ji
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
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          
          ! Assemble the global operator
          K(ij) = K(ij)+k_ij; K(ii) = K(ii)-k_ij
          K(ji) = K(ji)+k_ji; K(jj) = K(jj)-k_ji
        END DO
      END DO
    END SUBROUTINE doGalerkinZRSMat9_2D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 2D.
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doGalerkinMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, K)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: K
      
      REAL(DP), DIMENSION(NDIM2D) :: v_ij,v_ji
      REAL(DP)                    :: k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji
      INTEGER(PREC_VECIDX)        :: i,j
            
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)-v_ij(2)*Cy(ii)

        ! Update the diagonal coefficient
        K(ii) = K(ii)+k_ii

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          
          ! Assemble the global operator
          K(ij) = K(ij)+k_ij; K(ji) = K(ji)+k_ji
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9_2D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 3D
    ! and assume zero row-sums.
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doGalerkinZRSMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, K)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: K
      
      REAL(DP), DIMENSION(NDIM3D) :: v_ij,v_ji
      REAL(DP)                    :: k_ij,k_ji
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
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          
          ! Assemble the global operator
          K(ij) = K(ij)+k_ij; K(ii) = K(ii)-k_ij
          K(ji) = K(ji)+k_ji; K(jj) = K(jj)-k_ji
        END DO
      END DO
    END SUBROUTINE doGalerkinZRSMat9_3D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 3D.
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doGalerkinMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, K)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: K
      
      REAL(DP), DIMENSION(NDIM3D) :: v_ij,v_ji
      REAL(DP)                    :: k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji
      INTEGER(PREC_VECIDX)        :: i,j
            
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)-v_ij(2)*Cy(ii)-v_ij(3)*Cz(ii)

        ! Update the diagonal coefficient
        K(ii) = K(ii)+k_ii

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          
          ! Assemble the global operator
          K(ij) = K(ij)+k_ij; K(ji) = K(ji)+k_ji
        END DO
      END DO
    END SUBROUTINE doGalerkinMat9_3D


    !**************************************************************
    ! Assemble low-order operator L in 1D
    ! and assume zero row-sum.
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doUpwindZRSMat7_1D(Kld, Kcol, Ksep, NEQ, Cx, u, L)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      
      REAL(DP), DIMENSION(NDIM1D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ij,k_ji
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
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)
          k_ji = -v_ji(1)*Cx(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
        END DO
      END DO
    END SUBROUTINE doUpwindZRSMat7_1D


    !**************************************************************
    ! Assemble low-order operator L in 1D.
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doUpwindMat7_1D(Kld, Kcol, Ksep, NEQ, Cx, u, L)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      
      REAL(DP), DIMENSION(NDIM1D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j
            
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)

        ! Update the diagonal coefficient
        L(ii) = L(ii)+k_ii
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)
          k_ji = -v_ji(1)*Cx(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
        END DO
      END DO
    END SUBROUTINE doUpwindMat7_1D


    !**************************************************************
    ! Assemble low-order operator L in 2D
    ! and assume zero row-sums.
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doUpwindZRSMat7_2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, L)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      
      REAL(DP), DIMENSION(NDIM2D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ij,k_ji
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
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
        END DO
      END DO
    END SUBROUTINE doUpwindZRSMat7_2D


    !**************************************************************
    ! Assemble low-order operator L in 2D.
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doUpwindMat7_2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, L)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      
      REAL(DP), DIMENSION(NDIM2D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j
            
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)-v_ij(2)*Cy(ii)

        ! Update the diagonal coefficient
        L(ii) = L(ii)+k_ii

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
        END DO
      END DO
    END SUBROUTINE doUpwindMat7_2D


    !**************************************************************
    ! Assemble low-order operator L in 3D
    ! and assume zero row-sum.
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doUpwindZRSMat7_3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, L)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      
      REAL(DP), DIMENSION(NDIM3D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ij,k_ji
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
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
        END DO
      END DO
    END SUBROUTINE doUpwindZRSMat7_3D


    !**************************************************************
    ! Assemble low-order operator L in 3D.
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doUpwindMat7_3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, L)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      
      REAL(DP), DIMENSION(NDIM3D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j
            
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)-v_ij(2)*Cy(ii)-v_ij(3)*Cz(ii)

        ! Update the diagonal coefficient
        L(ii) = L(ii)+k_ii
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
        END DO
      END DO
    END SUBROUTINE doUpwindMat7_3D

    
    !**************************************************************
    ! Assemble low-order operator L in 1D
    ! and assume zero row-sum.
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doUpwindZRSMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, L)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      
      REAL(DP), DIMENSION(NDIM1D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ij,k_ji
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
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i),u(j),i,j,v_ij,v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)
          k_ji = -v_ji(1)*Cx(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
        END DO
      END DO
    END SUBROUTINE doUpwindZRSMat9_1D


    !**************************************************************
    ! Assemble low-order operator L in 1D.
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doUpwindMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, L)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      
      REAL(DP), DIMENSION(NDIM1D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j
            
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)

        ! Update the diagonal coefficient
        L(ii) = L(ii)+k_ii
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i),u(j),i,j,v_ij,v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)
          k_ji = -v_ji(1)*Cx(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
        END DO
      END DO
    END SUBROUTINE doUpwindMat9_1D


    !**************************************************************
    ! Assemble low-order operator L in 2D
    ! and assume zero row-sum.
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doUpwindZRSMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, L)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      
      REAL(DP), DIMENSION(NDIM2D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ij,k_ji
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
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i),u(j),i,j,v_ij,v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
        END DO
      END DO
    END SUBROUTINE doUpwindZRSMat9_2D


    !**************************************************************
    ! Assemble low-order operator L in 2D.
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doUpwindMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, L)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      
      REAL(DP), DIMENSION(NDIM2D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)-v_ij(2)*Cy(ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii)+k_ii
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i),u(j),i,j,v_ij,v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
        END DO
      END DO
    END SUBROUTINE doUpwindMat9_2D


    !**************************************************************
    ! Assemble low-order operator L in 3D
    ! and assume zero row-sums.
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doUpwindZRSMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, L)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      
      REAL(DP), DIMENSION(NDIM3D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ij,k_ji
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
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i),u(j),i,j,v_ij,v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
        END DO
      END DO
    END SUBROUTINE doUpwindZRSMat9_3D


    !**************************************************************
    ! Assemble low-order operator L in 3D.
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE doUpwindMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, L)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      
      REAL(DP), DIMENSION(NDIM3D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j
            
      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)-v_ij(2)*Cy(ii)-v_ij(3)*Cz(ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii)+k_ii
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i),u(j),i,j,v_ij,v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
        END DO
      END DO
    END SUBROUTINE doUpwindMat9_3D


    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 1D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    SUBROUTINE doUpwindZRSMat7_AFC_1D(Kld, Kcol, Ksep, NEQ, Cx, u,&
        L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM1D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)
          k_ji = -v_ji(1)*Cx(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindZRSMat7_AFC_1D


    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 1D.
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE doUpwindMat7_AFC_1D(Kld, Kcol, Ksep, NEQ, Cx, u,&
        L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM1D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii)+k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)
          k_ji = -v_ji(1)*Cx(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindMat7_AFC_1D


    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 2D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    SUBROUTINE doUpwindZRSMat7_AFC_2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u,&
        L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM2D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindZRSMat7_AFC_2D


    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 2D.
    ! All matrices are stored in matrix format 7

    SUBROUTINE doUpwindMat7_AFC_2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u,&
        L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM2D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)-v_ij(2)*Cy(ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii)+k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindMat7_AFC_2D

    
    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 3D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    SUBROUTINE doUpwindZRSMat7_AFC_3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u,&
        L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM3D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindZRSMat7_AFC_3D


    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 3D.
    ! All matrices are stored in matrix format 7

    SUBROUTINE doUpwindMat7_AFC_3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u,&
        L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM3D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)-v_ij(2)*Cy(ii)-v_ij(3)*Cz(ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii)+k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindMat7_AFC_3D

    
    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 1D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    SUBROUTINE doUpwindZRSMat9_AFC_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u,&
        L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM1D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)
          k_ji = -v_ji(1)*Cx(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindZRSMat9_AFC_1D


    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 1D.
    ! All matrices are stored in matrix format 9

    SUBROUTINE doUpwindMat9_AFC_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u,&
        L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM1D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii)+k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)
          k_ji = -v_ji(1)*Cx(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindMat9_AFC_1D


    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 2D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    SUBROUTINE doUpwindZRSMat9_AFC_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u,&
        L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM2D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindZRSMat9_AFC_2D


    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 2D.
    ! All matrices are stored in matrix format 9

    SUBROUTINE doUpwindMat9_AFC_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u,&
        L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM2D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)-v_ij(2)*Cy(ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii)+k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindMat9_AFC_2D


    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge 
    ! orientation in 3D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    SUBROUTINE doUpwindZRSMat9_AFC_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u,&
        L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM3D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindZRSMat9_AFC_3D

    
    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge 
    ! orientation in 3D.
    ! All matrices are stored in matrix format 9

    SUBROUTINE doUpwindMat9_AFC_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u,&
        L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM3D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)-v_ij(2)*Cy(ii)-v_ij(3)*Cz(ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii)+k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindMat9_AFC_3D


    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 1D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    SUBROUTINE doUpwindZRSMat7_ordAFC_1D(Kld, Kcol, Ksep, NEQ, Cx, &
        u, L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM1D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)
          k_ji = -v_ji(1)*Cx(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          IF (k_ij < k_ji) THEN
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          ELSE
            IverticesAtEdge(:,iedge)     =(/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) =(/d_ij, k_ji, k_ij/)
          END IF
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindZRSMat7_ordAFC_1D


    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 1D.
    ! All matrices are stored in matrix format 7

    SUBROUTINE doUpwindMat7_ordAFC_1D(Kld, Kcol, Ksep, NEQ, Cx, &
        u, L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM1D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii)+k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)
          k_ji = -v_ji(1)*Cx(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          IF (k_ij < k_ji) THEN
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          ELSE
            IverticesAtEdge(:,iedge)     =(/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) =(/d_ij, k_ji, k_ij/)
          END IF
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindMat7_ordAFC_1D
    
    
    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 2D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    SUBROUTINE doUpwindZRSMat7_ordAFC_2D(Kld, Kcol, Ksep, NEQ, Cx, Cy,&
        u, L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM2D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          IF (k_ij < k_ji) THEN
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          ELSE
            IverticesAtEdge(:,iedge)     =(/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) =(/d_ij, k_ji, k_ij/)
          END IF
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindZRSMat7_ordAFC_2D


    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 2D.
    ! All matrices are stored in matrix format 7

    SUBROUTINE doUpwindMat7_ordAFC_2D(Kld, Kcol, Ksep, NEQ, Cx, Cy,&
        u, L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM2D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)-v_ij(2)*Cy(ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii)+k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          IF (k_ij < k_ji) THEN
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          ELSE
            IverticesAtEdge(:,iedge)     =(/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) =(/d_ij, k_ji, k_ij/)
          END IF
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindMat7_ordAFC_2D
    
    
    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 3D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    SUBROUTINE doUpwindZRSMat7_ordAFC_3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz,&
        u, L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM3D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          IF (k_ij < k_ji) THEN
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          ELSE
            IverticesAtEdge(:,iedge)     =(/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) =(/d_ij, k_ji, k_ij/)
          END IF
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindZRSMat7_ordAFC_3D


    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 3D.
    ! All matrices are stored in matrix format 7

    SUBROUTINE doUpwindMat7_ordAFC_3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz,&
        u, L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM3D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)-v_ij(2)*Cy(ii)-v_ij(3)*Cz(ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii)+k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          IF (k_ij < k_ji) THEN
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          ELSE
            IverticesAtEdge(:,iedge)     =(/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) =(/d_ij, k_ji, k_ij/)
          END IF
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindMat7_ordAFC_3D
    
    
    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 1D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    SUBROUTINE doUpwindZRSMat9_ordAFC_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, &
        u, L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM1D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)
          k_ji = -v_ji(1)*Cx(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          IF (k_ij < k_ji) THEN
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          ELSE
            IverticesAtEdge(:,iedge)     = (/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ji, k_ij/)
          END IF
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindZRSMat9_ordAFC_1D


    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 1D.
    ! All matrices are stored in matrix format 9

    SUBROUTINE doUpwindMat9_ordAFC_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, &
        u, L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM1D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii)+k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)
          k_ji = -v_ji(1)*Cx(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          IF (k_ij < k_ji) THEN
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          ELSE
            IverticesAtEdge(:,iedge)     = (/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ji, k_ij/)
          END IF
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindMat9_ordAFC_1D


    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 2D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    SUBROUTINE doUpwindZRSMat9_ordAFC_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy,&
        u, L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM2D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          IF (k_ij < k_ji) THEN
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          ELSE
            IverticesAtEdge(:,iedge)     = (/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ji, k_ij/)
          END IF
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindZRSMat9_ordAFC_2D


    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 2D.
    ! All matrices are stored in matrix format 9

    SUBROUTINE doUpwindMat9_ordAFC_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy,&
        u, L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM2D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)-v_ij(2)*Cy(ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii)+k_ii

        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          IF (k_ij < k_ji) THEN
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          ELSE
            IverticesAtEdge(:,iedge)     = (/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ji, k_ij/)
          END IF
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindMat9_ordAFC_2D

    
    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 3D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    SUBROUTINE doUpwindZRSMat9_ordAFC_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz,&
        u, L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM3D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          IF (k_ij < k_ji) THEN
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          ELSE
            IverticesAtEdge(:,iedge)     = (/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ji, k_ij/)
          END IF
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindZRSMat9_ordAFC_3D
    
    
    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 3D.
    ! All matrices are stored in matrix format 9

    SUBROUTINE doUpwindMat9_ordAFC_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz,&
        u, L, IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz,u
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge
      
      ! local variables
      REAL(DP), DIMENSION(NDIM3D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,k_ii,k_ij,k_ji
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX)        :: i,j
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute local velocity coefficients for diagonal
        CALL fcb_getVelocity(u(i), u(i), i, i, v_ij, v_ji)
        
        ! Convective transport coefficients
        k_ii = -v_ij(1)*Cx(ii)-v_ij(2)*Cy(ii)-v_ij(3)*Cz(ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii)+k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute local velocity coefficients
          CALL fcb_getVelocity(u(i), u(j), i, j, v_ij, v_ji)
          
          ! Convective transport coefficients
          k_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          k_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          IF (k_ij < k_ji) THEN
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          ELSE
            IverticesAtEdge(:,iedge)     = (/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ji, k_ij/)
          END IF
        END DO
      END DO
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE doUpwindMat9_ordAFC_3D
  END SUBROUTINE gfsc_buildConvOperatorScalar

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_buildDiffusionOperator(rmatrixS, dscale,&
      bStabilise, bclear, rafcstab, rmatrixDest)

!<description>
    ! This subroutine assembles the diffusive part of the discrete
    ! transport operator which results from the discretisation of
    ! the scalar convection-diffusion-reaction equation.
    ! The matrix rmatrixS holds the unmodified diffusion operator
    ! which is possible modified by this routine. If the optional
    ! argument rmatrixDest is given, then the content of rmatrixS
    ! is first copied/added to the destination matrix prior to
    ! performing further modifications.
    ! If the parameter bclear is TRUE the destination matrix is cleared,
    ! that is, the content of rmatrixS is copied. Otherwise, the its
    ! content is combined linearly with that of the destination matrix.
    ! If the parameter bStabilse is TRUE, then symmetric stabilisation
    ! is applied so that the resulting diffusion operator is guaranted
    ! to ensure the discrete maximum principle. 
    ! If the optional parameter rafcstab is given, then the auxiliary
    ! data structures for of the stabilisation structure are updated.
    !
    ! At the moment, there is not yet any reference for this technique.
!</description>

!<input>
    ! scaling parameter by which the diffusion operator is scaled
    REAL(DP), INTENT(IN)                          :: dscale

    ! Switch for stabilisation
    ! TRUE  : perform stabilisation
    ! FALSE : perform no stabilisation
    LOGICAL, INTENT(IN)                           :: bStabilise

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    LOGICAL, INTENT(IN)                           :: bclear
!</input>

!<inputoutput>
    ! (anisotropic) diffusion operator
    TYPE(t_matrixScalar), INTENT(INOUT)           :: rmatrixS

    ! OPTIONAL: stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT), OPTIONAL      :: rafcstab

    ! OPTIONAL: destination matrix
    TYPE(t_matrixScalar), INTENT(INOUT), OPTIONAL :: rmatrixDest
!</inputoutput>
!</subroutine>


    ! local variables
    INTEGER(PREC_MATIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER   :: p_IsuperdiagonalEdgesIdx
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_Kld,p_Ksep
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_Kdiagonal
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER   :: p_Kcol
    REAL(DP), DIMENSION(:,:), POINTER             :: p_DcoefficientsAtEdge
    REAL(DP), DIMENSION(:), POINTER               :: p_S,p_L
    INTEGER :: h_Ksep

    
    ! Is destination matrix given?
    IF (PRESENT(rmatrixDest)) THEN

      ! Check if matrices are compatible
      CALL lsyssc_isMatrixCompatible(rmatrixS, rmatrixDest)

      ! Set pointers
      CALL lsyssc_getbase_double(rmatrixS,    p_S)
      CALL lsyssc_getbase_double(rmatrixDest, p_L)      

      ! Should matrix be cleared?
      IF (bclear) THEN
        CALL lalg_copyVectorDble (p_S, p_L)
        IF (dscale .NE. 1.0_DP) CALL lalg_scaleVectorDble(p_L, dscale)
      ELSE
        CALL lalg_vectorLinearCombDble(p_S, p_L, dscale, 1._DP)
      END IF

    ELSE

      ! Set pointers
      CALL lsyssc_getbase_double(rmatrixS, p_S)
      CALL lsyssc_getbase_double(rmatrixS, p_L)

      ! Scale matrix if required
      IF (dscale .NE. 1.0_DP) CALL lalg_scaleVectorDble(p_L, dscale)

    END IF

    ! Set pointers
    CALL lsyssc_getbase_Kld   (rmatrixS, p_Kld)
    CALL lsyssc_getbase_Kcol  (rmatrixS, p_Kcol)
    
    ! Do we have to stabilise?
    IF (PRESENT(rafcstab)) THEN

      ! Check if stabilisation has been initialised
      IF (IAND(rafcstab%iSpec, AFCSTAB_INITIALISED) .EQ. 0) THEN
        CALL output_line('Stabilisation has not been initialised',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildDiffOperatorScalar')
        CALL sys_halt()
      END IF

      ! Check if matrix and stabilisation
      ! structure are compatible to each other
      CALL gfsc_isMatrixCompatible(rafcstab, rmatrixS)

      ! What kind of stabilisation are we?
      SELECT CASE (rafcstab%ctypeAFCstabilisation)
      CASE (AFCSTAB_DMP)

        ! What kind of matrix are we?
        SELECT CASE(rmatrixS%cmatrixFormat)
        CASE(LSYSSC_MATRIX7)
          
          ! Create diagonal separator
          h_Ksep = ST_NOHANDLE
          CALL storage_copy(rmatrixS%h_Kld, h_Ksep)
          CALL storage_getbase_int(h_Ksep, p_Ksep, rmatrixS%NEQ+1)
          
          CALL do_loworderMat7(p_Kld, p_Kcol, p_Ksep, rmatrixS%NEQ, p_S, p_L)
          
          ! Release diagonal separator
          CALL storage_free(h_Ksep)
          
          
        CASE(LSYSSC_MATRIX9)
          
          ! Set pointers
          CALL lsyssc_getbase_Kdiagonal(rmatrixS, p_Kdiagonal)
          
          ! Create diagonal separator
          h_Ksep = ST_NOHANDLE
          CALL storage_copy(rmatrixS%h_Kld, h_Ksep)
          CALL storage_getbase_int(h_Ksep, p_Ksep, rmatrixS%NEQ+1)
          
          CALL do_loworderMat9(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                               rmatrixS%NEQ, p_S, p_L)
          
          ! Release diagonal separator
          CALL storage_free(h_Ksep)
          
        CASE DEFAULT
          CALL output_line('Unsupported matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildDiffOperatorScalar')
          CALL sys_halt()
        END SELECT


      CASE (AFCSTAB_SYMMETRIC)

        ! Set additional pointers
        CALL afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
        CALL afcstab_getbase_IverticesAtEdge(rafcstab,  p_IverticesAtEdge)
        CALL afcstab_getbase_DcoeffsAtEdge(rafcstab,    p_DcoefficientsAtEdge)
        
        ! What kind of matrix are we?
        SELECT CASE(rmatrixS%cmatrixFormat)
        CASE(LSYSSC_MATRIX7)
          
          ! Create diagonal separator
          h_Ksep = ST_NOHANDLE
          CALL storage_copy(rmatrixS%h_Kld, h_Ksep)
          CALL storage_getbase_int(h_Ksep, p_Ksep, rmatrixS%NEQ+1)
          
          CALL do_loworder_afcMat7(p_Kld, p_Kcol, p_Ksep, rmatrixS%NEQ,&
                                   p_S, p_L, p_IsuperdiagonalEdgesIdx,&
                                   p_IverticesAtEdge, p_DcoefficientsAtEdge)

          ! Set state of stabilisation
          rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)
          rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_EDGEVALUES)

          ! Release diagonal separator
          CALL storage_free(h_Ksep)
          
          
        CASE(LSYSSC_MATRIX9)

          ! Set pointers
          CALL lsyssc_getbase_Kdiagonal(rmatrixS, p_Kdiagonal)
          
          ! Create diagonal separator
          h_Ksep = ST_NOHANDLE
          CALL storage_copy(rmatrixS%h_Kld, h_Ksep)
          CALL storage_getbase_int(h_Ksep, p_Ksep, rmatrixS%NEQ+1)
          
          CALL do_loworder_afcMat9(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                   rmatrixS%NEQ, p_S, p_L, p_IsuperdiagonalEdgesIdx,&
                                   p_IverticesAtEdge, p_DcoefficientsAtEdge)

          ! Set state of stabilisation
          rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)
          rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_EDGEVALUES)

          ! Release diagonal separator
          CALL storage_free(h_Ksep)
          
          
        CASE DEFAULT
          CALL output_line('Unsupported matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildDiffOperatorScalar')
          CALL sys_halt()
        END SELECT
        
        
      CASE DEFAULT
        CALL output_line('Invalid type of AFC stabilisation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildDiffOperatorScalar')
        CALL sys_halt()
      END SELECT
      
    ELSEIF (bStabilise) THEN
      
      ! What kind of matrix are we?
      SELECT CASE(rmatrixS%cmatrixFormat)
      CASE(LSYSSC_MATRIX7)
        
        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        CALL storage_copy(rmatrixS%h_Kld, h_Ksep)
        CALL storage_getbase_int(h_Ksep, p_Ksep, rmatrixS%NEQ+1)
        
        CALL do_loworderMat7(p_Kld, p_Kcol, p_Ksep, rmatrixS%NEQ, p_S, p_L)
        
        ! Release diagonal separator
        CALL storage_free(h_Ksep)
        
        
      CASE(LSYSSC_MATRIX9)
        
        ! Set pointers
        CALL lsyssc_getbase_Kdiagonal(rmatrixS, p_Kdiagonal)
        
        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        CALL storage_copy(rmatrixS%h_Kld, h_Ksep)
        CALL storage_getbase_int(h_Ksep, p_Ksep, rmatrixS%NEQ+1)
        
        CALL do_loworderMat9(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                             rmatrixS%NEQ, p_S, p_L)
        
        ! Release diagonal separator
        CALL storage_free(h_Ksep)
        
      CASE DEFAULT
        CALL output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildDiffOperatorScalar')
        CALL sys_halt()
      END SELECT
    END IF
    
  CONTAINS
    
    ! Here, the working routine follow
    
    !**************************************************************
    ! Assemble low-order diffusion operator S.
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE do_loworderMat7(Kld, Kcol, Ksep, NEQ, S, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: S
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      
      REAL(DP)             :: d_ij
      INTEGER(PREC_MATIDX) :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX) :: i,j

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
          
          ! Artificial diffusion coefficient
          d_ij = MAX(0._DP, -S(ij)) 
          
          ! Assemble the global operator
          L(ii) = L(ii)-d_ij; L(ij) = L(ij)+d_ij
          L(ji) = L(ji)+d_ij; L(jj) = L(jj)-d_ij
        END DO
      END DO
    END SUBROUTINE do_loworderMat7

    
    !**************************************************************
    ! Assemble low-order diffusion operator S.
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE do_loworderMat9(Kld, Kcol, Kdiagonal, Ksep, NEQ, S, L)
      
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: S
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      
      REAL(DP)             :: d_ij
      INTEGER(PREC_MATIDX) :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX) :: i,j

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
          
          ! Artificial diffusion coefficient
          d_ij = MAX(0._DP, -S(ij)) 
          
          ! Assemble the global operator
          L(ii) = L(ii)-d_ij; L(ij) = L(ij)+d_ij
          L(ji) = L(ji)+d_ij; L(jj) = L(jj)-d_ij
        END DO
      END DO
    END SUBROUTINE do_loworderMat9


    !**************************************************************
    ! Assemble low-order diffusion operator S and AFC data.
    ! All matrices are stored in matrix format 7
    
    SUBROUTINE do_loworder_afcMat7(Kld, Kcol, Ksep, NEQ, S, L,&
        IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: S
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge

      REAL(DP)             :: d_ij,s_ij
      INTEGER(PREC_MATIDX) :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX) :: i,j

      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Artificial diffusion coefficient
          d_ij = MAX(0._DP, -S(ij))
          s_ij = MAX(0._DP,  S(ij))
          
          ! Assemble the global operator
          L(ii) = L(ii)-d_ij; L(ij) = L(ij)+d_ij
          L(ji) = L(ji)+d_ij; L(jj) = L(jj)-d_ij

          ! Increase edge counter
          iedge = iedge+1

          ! AFC
          IverticesAtEdge(:,iedge)     = (/i, j/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, s_ij/)
        END DO
      END DO

      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE do_loworder_afcMat7


    !**************************************************************
    ! Assemble low-order diffusion operator S and AFC data.
    ! All matrices are stored in matrix format 9
    
    SUBROUTINE do_loworder_afcMat9(Kld, Kcol, Kdiagonal, Ksep, NEQ, S, L,&
        IsuperdiagonalEdgesIdx, IverticesAtEdge, DcoefficientsAtEdge)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: S
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT)   :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: IverticesAtEdge
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: L
      REAL(DP), DIMENSION(:,:), INTENT(OUT)             :: DcoefficientsAtEdge

      REAL(DP)             :: d_ij,s_ij
      INTEGER(PREC_MATIDX) :: ii,ij,ji,jj,iedge
      INTEGER(PREC_VECIDX) :: i,j

      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      DO i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Artificial diffusion coefficient
          d_ij = MAX(0._DP, -S(ij))
          s_ij = MAX(0._DP,  S(ij))
          
          ! Assemble the global operator
          L(ii) = L(ii)-d_ij; L(ij) = L(ij)+d_ij
          L(ji) = L(ji)+d_ij; L(jj) = L(jj)-d_ij

          ! Increase edge counter
          iedge = iedge+1

          ! AFC
          IverticesAtEdge(:,iedge)     = (/i, j/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, s_ij/)
        END DO
      END DO

      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    END SUBROUTINE do_loworder_afcMat9    
  END SUBROUTINE gfsc_buildDiffusionOperator

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_buildResBlockFCT(rmatrixMC, rmatrixML, ru,&
      theta, tstep, binitResidual, rres, rafcstab)

!<description>
    ! This subroutine assembles the residual vector and applies
    ! stabilisation of FEM-FCT type.  Note that this routine serves as
    ! a wrapper for block vectors. If there is only one block, then
    ! the corresponding scalar routine is called.  Otherwise, an error
    ! is thrown.
!</description>

!<input>
    ! consistent mass matrix
    TYPE(t_matrixScalar), INTENT(IN)   :: rmatrixMC

    ! lumped mass matrix
    TYPE(t_matrixScalar), INTENT(IN)   :: rmatrixML

    ! solution vector
    TYPE(t_vectorBlock), INTENT(IN)    :: ru

    ! implicitness parameter
    REAL(DP), INTENT(IN)               :: theta

    ! time step size
    REAL(DP), INTENT(IN)               :: tstep
    
    ! Switch for residual
    ! TRUE  : build the initial residual
    ! FALSE : build an intermediate residual
    LOGICAL, INTENT(IN)                :: binitResidual
!</input>

!<inputoutput>
    ! residual vector
    TYPE(t_vectorBlock), INTENT(INOUT) :: rres

    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)     :: rafcstab
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    IF (ru%nblocks .NE. 1 .OR. rres%nblocks .NE. 1) THEN

      CALL output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResBlockFCT')
      CALL sys_halt()

    ELSE

      CALL gfsc_buildResScalarFCT(rmatrixMC, rmatrixML, ru%RvectorBlock(1),&
          theta, tstep, binitResidual, rres%RvectorBlock(1), rafcstab)

    END IF
  END SUBROUTINE gfsc_buildResBlockFCT

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_buildResScalarFCT(rmatrixMC, rmatrixML, ru,&
      theta, tstep, binitResidual, rres, rafcstab)

!<description>

    ! This subroutine assembles the residual vector and applies
    ! stabilisation of FEM-FCT type. The idea of flux corrected
    ! transport can be traced back to the early SHASTA algorithm by
    ! Boris and Bock in the early 1970s. Zalesak suggested a fully
    ! multi-dimensional generalisation of this approach and paved the
    ! way for a large family of FCT algorithms.
    !
    ! This subroutine provides different algorithms:
    !
    ! 1. Semi-explicit FEM-FCT algorithm
    !
    !    This is the classical algorithm which makes use of Zalesak's
    !    flux limiter and recomputes and auxiliary positivity-
    !    preserving solution in each iteration step. 
    !    The details of this method can be found in:
    !
    !    D. Kuzmin and M. Möller, "Algebraic flux correction I. Scalar
    !    conservation laws", Ergebnisberichte Angew. Math. 249,
    !    University of Dortmund, 2004.
    !
    ! 2. Iterative FEM-FCT algorithm
    !
    !    The main idea of the iterative algorithm is to reuse the
    !    amount of rejected antidiffusion in subsequent iterations so
    !    that more and more antidiffusion can be built into the
    !    residual as the iteration process continues. 
    !    The details of this method can be also found in the above reference.
    !
    ! 3. Semi-implicit FEM-FCT algorith
    !
    !    This is the FCT algorithm that should be used by default. It
    !    is quite efficient since the nodal correction factors are
    !    only computed in the first iteration and used to limit the
    !    antidiffusive flux from the first iteration. This explicit
    !    predictor is used in all subsequent iterations to constrain
    !    the actual target flux.
    !    The details of this method can be found in:
    !
    !    D. Kuzmin and D. Kourounis, " A semi-implicit FEM-FCT
    !    algorithm for efficient treatment of time-dependent
    !    problems", Ergebnisberichte Angew. Math. 302, University of
    !    Dortmund, 2005.
    !
    ! 4. Linearized FEM-FCT algorithms
    !
    !    A new trend in the development of FCT algorithms is to
    !    linearise the raw antidiffusive fluxes about an intermediate
    !    solution computed by a positivity-preserving low-order
    !    scheme. By virtue of this linearisation, the costly
    !    evaluation of correction factors needs to be performed just
    !    once per time step. Furthermore, no questionable
    !    `prelimiting' of antidiffusive fluxes is required, which
    !    eliminates the danger of artificial steepening.
    !    The details of this method can be found in:
    !
    !    D. Kuzmin, "Explicit and implicit FEM-FCT algorithms with
    !    flux linearization", Ergebnisberichte Angew. Math. 358,
    !    University of Dortmund, 2008.
!</description>

!<input>
    ! consistent mass matrix
    TYPE(t_matrixScalar), INTENT(IN)    :: rmatrixMC

    ! lumped mass matrix
    TYPE(t_matrixScalar), INTENT(IN)    :: rmatrixML

    ! solution vector
    TYPE(t_vectorScalar), INTENT(IN)    :: ru

    ! implicitness parameter
    REAL(DP), INTENT(IN)                :: theta

    ! time step size
    REAL(DP), INTENT(IN)                :: tstep
    
    ! Switch for residual
    ! TRUE  : build the initial residual
    ! FALSE : build an intermediate residual
    LOGICAL, INTENT(IN)                 :: binitResidual
!</input>

!<inputoutput>
    ! residual vector
    TYPE(t_vectorScalar), INTENT(INOUT) :: rres

    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)      :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    TYPE(t_vectorScalar), POINTER                 :: ruLow
    INTEGER(PREC_VECIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge
    REAL(DP), DIMENSION(:,:), POINTER             :: p_DcoefficientsAtEdge
    REAL(DP), DIMENSION(:), POINTER :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    REAL(DP), DIMENSION(:), POINTER :: p_flux,p_flux0
    REAL(DP), DIMENSION(:), POINTER :: p_u,p_ulow,p_res
    REAL(DP), DIMENSION(:), POINTER :: p_MC,p_ML

   
    ! Check if stabilisation is prepared
    IF (rafcstab%ctypeAFCstabilisation .NE. AFCSTAB_FEMFCT  .OR.&
        IAND(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE) .EQ. 0  .OR.&
        IAND(rafcstab%iSpec, AFCSTAB_EDGEVALUES)    .EQ. 0) THEN
      CALL output_line('Stabilisation does not provide required structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResScalarFCT')
      CALL sys_halt()
    END IF

    ! Check if vectors are compatible
    CALL lsyssc_isVectorCompatible(ru, rres)

    ! Set pointers
    CALL afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    CALL afcstab_getbase_DcoeffsAtEdge(rafcstab,   p_DcoefficientsAtEdge)
    CALL lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
    CALL lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
    CALL lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
    CALL lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
    CALL lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
    CALL lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
    CALL lsyssc_getbase_double(rafcstab%RedgeVectors(1),  p_flux)
    CALL lsyssc_getbase_double(rafcstab%RedgeVectors(2),  p_flux0)
    CALL lsyssc_getbase_double(rmatrixMC, p_MC)
    CALL lsyssc_getbase_double(rmatrixML, p_ML)
    CALL lsyssc_getbase_double(ru,   p_u)
    CALL lsyssc_getbase_double(rres, p_res)

    ! Should we build up the initial residual?
    IF (binitResidual) THEN
      
      ! Do we have a fully implicit time discretisation?
      IF (theta < 1.0_DP) THEN
        ruLow => rafcstab%RnodalVectors(5)
        CALL lsyssc_invertedDiagMatVec(rmatrixML, rres, 1._DP, ruLow)
        CALL lsyssc_vectorLinearComb(ru, ruLow, 1._DP, 1._DP-theta)
        CALL lsyssc_getbase_double(ruLow, p_ulow)
        
        CALL do_femfct_init(p_IverticesAtEdge, p_DcoefficientsAtEdge, p_MC, p_ML,&
            p_u, p_ulow, theta, tstep, rafcstab%NEDGE,&
            (rafcstab%imass .EQ. AFCSTAB_CONSISTENTMASS),&
            p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_flux, p_flux0)
        
      ELSE
        
        CALL do_femfct_init(p_IverticesAtEdge, p_DcoefficientsAtEdge, p_MC, p_ML,&
            p_u, p_u, theta, tstep, rafcstab%NEDGE,&
            (rafcstab%imass .EQ. AFCSTAB_CONSISTENTMASS),&
            p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_flux, p_flux0)
        
      END IF
      
      ! Set specifier
      rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_LIMITER)
      rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_FLUXES)        
    END IF
    
    ! Check if correction factors and fluxes are available
    IF (IAND(rafcstab%iSpec, AFCSTAB_LIMITER) .EQ. 0 .OR.&
        IAND(rafcstab%iSpec, AFCSTAB_FLUXES)  .EQ. 0) THEN
      CALL output_line('Stabilisation does not provide precomputed fluxes &
          &and/or nodal correction factors',OU_CLASS_ERROR,OU_MODE_STD,&
          'gfsc_buildResScalarFCT')
      CALL sys_halt()
    END IF
    
    ! Apply the limited
    CALL do_femfct_limit(p_IverticesAtEdge, p_DcoefficientsAtEdge, p_MC,&
        p_u, p_flux, p_flux0, theta, tstep, rafcstab%NEDGE,&
        (rafcstab%imass .EQ. AFCSTAB_CONSISTENTMASS), p_res)
    
  CONTAINS

    ! Here, the working routine follow
    
    !**************************************************************
    ! Initialisation of the semi-implicit FEM-FCT procedure
    
    SUBROUTINE do_femfct_init(IverticesAtEdge, DcoefficientsAtEdge, MC, ML,&
        u, ulow, theta, tstep, NEDGE, bmass, pp, pm, qp, qm, rp, rm, flux, flux0)
      
      INTEGER(PREC_VECIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtEdge
      REAL(DP), DIMENSION(:,:), INTENT(IN)             :: DcoefficientsAtEdge
      REAL(DP), DIMENSION(:), INTENT(IN)               :: MC,ML
      REAL(DP), DIMENSION(:), INTENT(IN)               :: u,ulow
      REAL(DP), INTENT(IN)                             :: theta,tstep
      INTEGER(PREC_MATIDX), INTENT(IN)                 :: NEDGE
      LOGICAL, INTENT(IN)                              :: bmass
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: pp,pm
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: qp,qm
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: rp,rm
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: flux,flux0
      
      INTEGER(PREC_MATIDX) :: iedge,ij
      INTEGER(PREC_VECIDX) :: i,j
      REAL(DP) :: d_ij,f_ij,m_ij,diff

      ! Clear nodal vectors
      CALL lalg_clearVectorDble(pp)
      CALL lalg_clearVectorDble(pm)
      CALL lalg_clearVectorDble(qp)
      CALL lalg_clearVectorDble(qm)

      ! Should we apply the consistent mass matrix?
      IF (bmass) THEN

        ! Loop over edges
        DO iedge = 1, NEDGE
          
          ! Determine indices
          i  = IverticesAtEdge(1,iedge)
          j  = IverticesAtEdge(2,iedge)
          ij = IverticesAtEdge(3,iedge)
          
          ! Determine coefficients
          d_ij = DcoefficientsAtEdge(1,iedge); m_ij = MC(ij)
          
          ! Determine fluxes
          diff = u(i)-u(j);   f_ij = tstep*d_ij*diff
          flux(iedge) = f_ij; flux0(iedge) = -m_ij*diff
          
          ! Sum of positive/negative fluxes
          pp(i) = pp(i)+MAX(0._DP,f_ij); pp(j) = pp(j)+MAX(0._DP,-f_ij)
          pm(i) = pm(i)+MIN(0._DP,f_ij); pm(j) = pm(j)+MIN(0._DP,-f_ij)
          
          ! Upper/lower bounds
          diff = ulow(j)-ulow(i)
          qp(i) = MAX(qp(i),diff); qp(j) = MAX(qp(j),-diff)
          qm(i) = MIN(qm(i),diff); qm(j) = MIN(qm(j),-diff)
        END DO

        ! Adopt the explicit part (if required)
        IF (theta < 1.0_DP) THEN
          CALL lalg_vectorLinearCombDble(flux, flux0, 1.0_DP-theta, 1.0_DP)
        END IF
        
      ELSE

        ! Loop over edges
        DO iedge = 1, NEDGE
          
          ! Determine indices
          i  = IverticesAtEdge(1,iedge)
          j  = IverticesAtEdge(2,iedge)
          ij = IverticesAtEdge(3,iedge)
          
          ! Determine coefficients
          d_ij = DcoefficientsAtEdge(1,iedge)
          
          ! Determine fluxes
          diff = u(i)-u(j);   f_ij = tstep*d_ij*diff
          flux(iedge) = f_ij
          
          ! Sum of positive/negative fluxes
          pp(i) = pp(i)+MAX(0._DP,f_ij); pp(j) = pp(j)+MAX(0._DP,-f_ij)
          pm(i) = pm(i)+MIN(0._DP,f_ij); pm(j) = pm(j)+MIN(0._DP,-f_ij)
          
          ! Upper/lower bounds
          diff = ulow(j)-ulow(i)
          qp(i) = MAX(qp(i),diff); qp(j) = MAX(qp(j),-diff)
          qm(i) = MIN(qm(i),diff); qm(j) = MIN(qm(j),-diff)
        END DO

        ! Adopt the explicit part (if required)
        IF (theta < 1.0_DP) THEN
          CALL lalg_copyVectorDble(flux, flux0)
          CALL lalg_scaleVectorDble(flux0, 1.0_DP-theta)
        ELSE
          CALL lalg_clearVectorDble(flux0)
        END IF

      END IF

      ! Apply the nodal limiter
      rp = ML*qp; rp = afcstab_limit( pp, rp, 0._DP)
      rm =-ML*qm; rm = afcstab_limit(-pm, rm, 0._DP)
      
      ! Limiting procedure
      DO iedge = 1, NEDGE

        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        IF (flux(iedge) > 0.0_DP) THEN
          flux(iedge) = MIN(rp(i), rm(j))*flux(iedge)
        ELSE
          flux(iedge) = MIN(rm(i), rp(j))*flux(iedge)
        END IF
      END DO
    END SUBROUTINE do_femfct_init
    

    !**************************************************************
    ! The semi-implicit FEM-FCT limiting procedure
    
    SUBROUTINE do_femfct_limit(IverticesAtEdge, DcoefficientsAtEdge, MC,&
        u, flux, flux0, theta, tstep, NEDGE, bmass, res)

      INTEGER(PREC_VECIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtEdge
      REAL(DP), DIMENSION(:,:), INTENT(IN)             :: DcoefficientsAtEdge
      REAL(DP), DIMENSION(:), INTENT(IN)               :: MC
      REAL(DP), DIMENSION(:), INTENT(IN)               :: u
      REAL(DP), DIMENSION(:), INTENT(IN)               :: flux,flux0
      REAL(DP), INTENT(IN)                             :: theta,tstep
      INTEGER(PREC_MATIDX), INTENT(IN)                 :: NEDGE
      LOGICAL, INTENT(IN)                              :: bmass
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: res
      
      INTEGER(PREC_MATIDX) :: iedge,ij
      INTEGER(PREC_VECIDX) :: i,j
      REAL(DP) :: d_ij,f_ij,m_ij

      ! Should we apply the consistent mass matrix
      IF (bmass) THEN
        
        ! Loop over edges
        DO iedge = 1, NEDGE
          
          ! Determine indices
          i  = IverticesAtEdge(1,iedge)
          j  = IverticesAtEdge(2,iedge)
          ij = IverticesAtEdge(3,iedge)
          
          ! Determine coefficients
          d_ij = DcoefficientsAtEdge(1,iedge); m_ij = MC(ij)
          
          ! Determine fluxes
          f_ij = (m_ij+theta*tstep*d_ij)*(u(i)-u(j))+flux0(iedge)
          
          IF (f_ij > 0.0_DP) THEN
            f_ij = MIN(f_ij,MAX(flux(iedge),0._DP))
          ELSE
            f_ij = MAX(f_ij,MIN(flux(iedge),0._DP))
          END IF
          
          ! Update the defect vector
          res(i) = res(i)+f_ij
          res(j) = res(j)-f_ij
        END DO
        
      ELSE
        
        ! Loop over edges
        DO iedge = 1, NEDGE
          
          ! Determine indices
          i  = IverticesAtEdge(1,iedge)
          j  = IverticesAtEdge(2,iedge)
          ij = IverticesAtEdge(3,iedge)
          
          ! Determine coefficients
          d_ij = DcoefficientsAtEdge(1,iedge)
          
          ! Determine fluxes
          f_ij = (theta*tstep*d_ij)*(u(i)-u(j))+flux0(iedge)
          
          IF (f_ij > 0.0_DP) THEN
            f_ij = MIN(f_ij,MAX(flux(iedge),0._DP))
          ELSE
            f_ij = MAX(f_ij,MIN(flux(iedge),0._DP))
          END IF
          
          ! Update the defect vector
          res(i) = res(i)+f_ij
          res(j) = res(j)-f_ij
        END DO
        
      END IF
    END SUBROUTINE do_femfct_limit
  END SUBROUTINE gfsc_buildResScalarFCT

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_buildResBlockTVD(rmatrixMC, ru, ru0,&
      theta, tstep, rres, rafcstab)

!<description>
    ! This subroutine assembles the residual vector
    ! and applies stabilisation of FEM-TVD type.
    ! Note that this routine serves as a wrapper for block vectors. If there
    ! is only one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
!</description>

!<input>
    ! consistent mass matrix
    TYPE(t_matrixScalar), INTENT(IN)   :: rmatrixMC

    ! solution vector
    TYPE(t_vectorBlock), INTENT(IN)    :: ru

    ! initial solution vector
    TYPE(t_vectorBlock), INTENT(IN)    :: ru0

    ! implicitness parameter
    REAL(DP), INTENT(IN)               :: theta

    ! time step size
    REAL(DP), INTENT(IN)               :: tstep
!</input>

!<inputoutput>
    ! residual vector
    TYPE(t_vectorBlock), INTENT(INOUT) :: rres

    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)     :: rafcstab
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    IF (ru%nblocks   .NE. 1 .OR.&
        ru0%nblocks  .NE. 1 .OR.&
        rres%nblocks .NE. 1) THEN

      CALL output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResBlockTVD')
      CALL sys_halt()

    ELSE

      CALL gfsc_buildResScalarTVD(rmatrixMC, ru%RvectorBlock(1),&
          ru0%RvectorBlock(1), theta, tstep, rres%RvectorBlock(1), rafcstab)
      
    END IF
  END SUBROUTINE gfsc_buildResBlockTVD
  
  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_buildResScalarTVD(rmatrixMC, ru, ru0,&
      theta, tstep, rres, rafcstab)

!<description>
    ! This subroutine assembles the residual vector and applies
    ! stabilisation of FEM-TVD type and/or employes the general
    ! purpose limiter.
    !
    ! A detailed description of the FEM-TVD limiter in general is given in:
    !
    !     D. Kuzmin and S. Turek, "Multidimensional FEM-TVD paradigm
    !     for convection-dominated flows" In:  Proceedings of the 
    !     IV European Congress on Computational Methods in Applied Sciences
    !     and Engineering (ECCOMAS 2004). Vol. II, ISBN 951-39-1869-6.
    !
    ! The method actually implemented in this routine is described in:
    !
    !     D. Kuzmin, "Algebraic flux correction for finite element
    !     discretizations of coupled systems" In: E. Onate,
    !     M. Papadrakakis and B. Schrefler (eds.) Computational
    !     Methods for Coupled Problems in Science and Engineering II,
    !     CIMNE, Barcelona, 2007, 653-656.
    !
    ! For time-dependent problems, the lack of a consistent mass
    ! matrix deteriorates the phase accuracy significnatly. It is thus
    ! that symmetric limiting strategies taylored to mass
    ! antidiffusion have been considered. If both convective and mass
    ! antidiffusion has to be limited then the so-called
    ! general-purpose (GP) flux limiter comes into play:
    !
    !     D. Kuzmin, "On the design of general-purpose flux 
    !     limiters for implicit FEM with a consistent mass matrix.
    !     I. Scalar convection."
    !     J. Comput. Phys.  219  (2006) 513-531.
    !
    ! Note however, that is is quite expensive and not recommended as
    ! a standard limiter. In fact, it is only implemented to
    ! demonstrate that the construction of general-purpose flux
    ! limiters is possible. If you want to recover the consistent mass
    ! matrix for time-dependent problems, you should apply flux
    ! correction of FCT type.
!</description>

!<input>
    ! consistent mass matrix
    TYPE(t_matrixScalar), INTENT(IN)    :: rmatrixMC

    ! solution vector
    TYPE(t_vectorScalar), INTENT(IN)    :: ru

    ! initial solution vector
    TYPE(t_vectorScalar), INTENT(IN)    :: ru0

    ! implicitness parameter
    REAL(DP), INTENT(IN)                :: theta

    ! time step size
    REAL(DP), INTENT(IN)                :: tstep
!</input>

!<inputoutput>
    ! residual vector
    TYPE(t_vectorScalar), INTENT(INOUT) :: rres

    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)      :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_VECIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge
    REAL(DP), DIMENSION(:,:), POINTER             :: p_DcoefficientsAtEdge
    REAL(DP), DIMENSION(:), POINTER :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    REAL(DP), DIMENSION(:), POINTER :: p_flux,p_flux0
    REAL(DP), DIMENSION(:), POINTER :: p_u,p_u0,p_res
    REAL(DP), DIMENSION(:), POINTER :: p_MC

    ! Check if vectors are compatible
    CALL lsyssc_isVectorCompatible(ru, rres)
    
    ! What kind of stabilisation are we?
    SELECT CASE(rafcstab%ctypeAFCstabilisation)
      
    CASE (AFCSTAB_FEMTVD)
      
      ! Check if stabilisation is prepared
      IF (IAND(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)  .EQ.0 .OR. &
          IAND(rafcstab%iSpec, AFCSTAB_EDGEORIENTATION).EQ.0 .OR. &
          IAND(rafcstab%iSpec, AFCSTAB_EDGEVALUES)     .EQ.0) THEN
        CALL output_line('Stabilisation does not provide required structures',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResScalarTVD')
        CALL sys_halt()
      END IF

      ! Set pointers
      CALL afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
      CALL afcstab_getbase_DcoeffsAtEdge(rafcstab,   p_DcoefficientsAtEdge)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
      CALL lsyssc_getbase_double(rafcstab%RedgeVectors(1),  p_flux)
      CALL lsyssc_getbase_double(ru,   p_u)
      CALL lsyssc_getbase_double(rres, p_res)

      CALL do_femtvd_limit(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
          p_u, tstep, rafcstab%NEDGE, p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_flux, p_res)

      ! Set specifier
      rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_BOUNDS)
      rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_ANTIDIFFUSION)
      rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_LIMITER)
      rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_FLUXES)
      

    CASE (AFCSTAB_FEMGP)
      
      ! Check if stabilisation is prepared
      IF (IAND(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)  .EQ.0 .OR. &
          IAND(rafcstab%iSpec, AFCSTAB_EDGEORIENTATION).EQ.0 .OR. &
          IAND(rafcstab%iSpec, AFCSTAB_EDGEVALUES)     .EQ.0) THEN
        CALL output_line('Stabilisation does not provide required structures',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResScalarTVD')
        CALL sys_halt()
      END IF

      ! Set pointers
      CALL afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
      CALL afcstab_getbase_DcoeffsAtEdge(rafcstab,   p_DcoefficientsAtEdge)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
      CALL lsyssc_getbase_double(rafcstab%RedgeVectors(1),  p_flux)
      CALL lsyssc_getbase_double(rafcstab%RedgeVectors(2),  p_flux0)
      CALL lsyssc_getbase_double(rmatrixMC, p_MC)
      CALL lsyssc_getbase_double(ru,   p_u)
      CALL lsyssc_getbase_double(ru0,  p_u0)
      CALL lsyssc_getbase_double(rres, p_res)

      IF (rafcstab%imass .EQ. AFCSTAB_CONSISTENTMASS) THEN
        CALL do_femgp_limit(p_IverticesAtEdge, p_DcoefficientsAtEdge, p_MC,&
            p_u, p_u0, theta, tstep, rafcstab%NEDGE,&
            p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_flux, p_flux0, p_res)
      ELSE
        CALL do_femtvd_limit(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
            p_u, tstep, rafcstab%NEDGE, p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_flux, p_res)
      END IF
      
      ! Set specifier
      rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_BOUNDS)
      rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_ANTIDIFFUSION)
      rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_LIMITER)
      rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_FLUXES)

    CASE DEFAULT
      CALL output_line('Invalid type of AFC stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResScalarTVD')
      CALL sys_halt()
    END SELECT

  CONTAINS

    ! Here, the working routine follow
    
    !**************************************************************
    ! The FEM-TVD limiting procedure
    
    SUBROUTINE do_femtvd_limit(IverticesAtEdge, DcoefficientsAtEdge,&
        u, tstep, NEDGE, pp, pm, qp, qm, rp, rm, flux, res)

      INTEGER(PREC_VECIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtEdge
      REAL(DP), DIMENSION(:,:), INTENT(IN)             :: DcoefficientsAtEdge
      REAL(DP), DIMENSION(:), INTENT(IN)               :: u
      REAL(DP), INTENT(IN)                             :: tstep
      INTEGER(PREC_MATIDX), INTENT(IN)                 :: NEDGE
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: pp,pm
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: qp,qm
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: rp,rm
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: flux
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: res
      
      INTEGER(PREC_MATIDX) :: iedge,ij
      INTEGER(PREC_VECIDX) :: i,j
      REAL(DP) :: d_ij,f_ij,l_ij,l_ji,diff

      ! Clear nodal vectors
      CALL lalg_clearVectorDble(pp)
      CALL lalg_clearVectorDble(pm)
      CALL lalg_clearVectorDble(qp)
      CALL lalg_clearVectorDble(qm)

      ! Assemble P's and Q's
      DO iedge = 1, NEDGE
        
        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Determine coefficients
        d_ij = DcoefficientsAtEdge(1,iedge)
        l_ij = DcoefficientsAtEdge(2,iedge)
        l_ji = DcoefficientsAtEdge(3,iedge)
        
        ! Determine solution difference
        diff = tstep*(u(i)-u(j))
        
        ! Prelimit the antidiffusive flux F'_IJ=MIN(-P_IJ,L_JI)(U_I-U_J)
        f_ij = MIN(d_ij,l_ji)*diff; flux(iedge) = f_ij
        
        ! Assemble P's accordingly
        pp(i) = pp(i)+MAX(0._DP,f_ij); pm(i) = pm(i)+MIN(0._DP,f_ij)
        
        ! Assemble Q's
        qp(i) = qp(i)+MAX(0._DP,-f_ij); qp(j) = qp(j)+MAX(0._DP, f_ij)
        qm(i) = qm(i)+MIN(0._DP,-f_ij); qm(j) = qm(j)+MIN(0._DP, f_ij)
      END DO
      
      ! Apply the nodal limiter
      rp = afcstab_limit( pp, qp, 0._DP, 1._DP)
      rm = afcstab_limit(-pm,-qm, 0._DP, 1._DP)

      ! Apply limiter
      DO iedge = 1, NEDGE
        
        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Get precomputed raw antidiffusive flux
        f_ij = flux(iedge)
        
        ! Apply correction factor and store limite flux
        f_ij = MERGE(rp(i), rm(i), f_ij > 0)*f_ij
        
        ! Update the defect vector
        res(i) = res(i)+f_ij
        res(j) = res(j)-f_ij
      END DO
    END SUBROUTINE do_femtvd_limit


    !**************************************************************
    ! The FEM-GP limiting procedure

    SUBROUTINE do_femgp_limit(IverticesAtEdge, DcoefficientsAtEdge, MC,&
        u, u0, theta, tstep, NEDGE, pp, pm, qp, qm, rp, rm, flux, flux0, res)

      INTEGER(PREC_VECIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtEdge
      REAL(DP), DIMENSION(:,:), INTENT(IN)             :: DcoefficientsAtEdge
      REAL(DP), DIMENSION(:), INTENT(IN)               :: MC
      REAL(DP), DIMENSION(:), INTENT(IN)               :: u,u0
      REAL(DP), INTENT(IN)                             :: theta,tstep
      INTEGER(PREC_MATIDX), INTENT(IN)                 :: NEDGE
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: pp,pm
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: qp,qm
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: rp,rm
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: flux,flux0
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: res
      
      INTEGER(PREC_MATIDX) :: iedge,ij
      INTEGER(PREC_VECIDX) :: i,j
      REAL(DP) :: d_ij,f_ij,l_ij,l_ji,m_ij,p_ij,pf_ij,df_ij,q_ij,q_ji
      REAL(DP) :: diff,diff0,diff1

      ! Clear nodal vectors
      CALL lalg_clearVectorDble(pp)
      CALL lalg_clearVectorDble(pm)
      CALL lalg_clearVectorDble(qp)
      CALL lalg_clearVectorDble(qm)

      ! Assemble P' and Q's
      DO iedge = 1, NEDGE
        
        ! Determine indices
        i  = IverticesAtEdge(1,iedge)
        j  = IverticesAtEdge(2,iedge)
        ij = IverticesAtEdge(3,iedge)
        
        ! Determine coefficients
        d_ij = DcoefficientsAtEdge(1,iedge)
        l_ij = DcoefficientsAtEdge(2,iedge)
        l_ji = DcoefficientsAtEdge(3,iedge)
        m_ij = MC(ij)
        
        ! Compute: diff1 = dt*theta*(u_i-u_j) + dt*(1-theta)*(u0_i-u0_j)
        diff1 = u(i)-u(j); diff0 = u0(i)-u0(j)
        diff  = tstep*(theta*diff1+(1._DP-theta)*diff0)
        
        ! Compute antidiffusive flux f_ij=min(0,p_ij)*(u_j-u_i)
        IF (ABS(diff) < SYS_EPSREAL) THEN
          p_ij = 0
          f_ij = 0
        ELSE
          p_ij = MAX(0._DP,m_ij*(diff1-diff0)/diff+d_ij)
          f_ij = p_ij*diff
        END IF
        
        ! Prelimit the antidiffusive flux F'_IJ=MIN(-P_IJ,L_JI)(U_I-U_J)
        pf_ij = MIN(p_ij,l_ji)*diff; flux0(iedge) = pf_ij
        
        ! Compute the remaining flux dF_IJ=F_IJ-F'_IJ
        df_ij = f_ij-pf_ij; flux(iedge) = df_ij
        
        ! Assemble P's accordingly
        pp(i) = pp(i)+MAX(0._DP,  f_ij); pm(i) = pm(i)+MIN(0._DP,  f_ij)
        pp(j) = pp(j)+MAX(0._DP,-df_ij); pm(j) = pm(j)+MIN(0._DP,-df_ij)
        
        q_ij = m_ij/tstep+l_ij; q_ji = m_ij/tstep+l_ji

        ! Assemble Q's
        qp(i) = qp(i)+q_ij*MAX(0._DP,-diff); qm(i) = qm(i)+q_ij*MIN(0._DP,-diff)
        qp(j) = qp(j)+q_ji*MAX(0._DP, diff); qm(j) = qm(j)+q_ji*MIN(0._DP, diff)
      END DO

      ! Apply nodal limiter
      rp = afcstab_limit( pp, qp, 0._DP, 1._DP)
      rm = afcstab_limit(-pm,-qm, 0._DP, 1._DP)

      ! Apply limiter
      DO iedge = 1, NEDGE

        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Get precomputed fluxes
        pf_ij = flux0(iedge); df_ij = flux(iedge)
        
        ! Limit upwind contribution
        IF (pf_ij > 0.0_DP) THEN
          pf_ij = rp(i)*pf_ij
        ELSE
          pf_ij = rm(i)*pf_ij
        END IF
        
        ! Limit symmetric contribution
        IF (df_ij > 0.0_DP) THEN
          df_ij = MIN(rp(i), rm(j))*df_ij
        ELSE
          df_ij = MIN(rm(i), rp(j))*df_ij
        END IF
        
        f_ij = pf_ij+df_ij
        
        ! Update the defect vector
        res(i) = res(i)+f_ij
        res(j) = res(j)-f_ij
      END DO
    END SUBROUTINE do_femgp_limit
  END SUBROUTINE gfsc_buildResScalarTVD

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_buildResBlockSymm(ru, dscale, rres, rafcstab)

!<description>
    ! This subroutine assembles the residual vector and applies stabilisation
    ! by means of symmetric flux limiting for diffusion operators.
    ! Note that this routine serves as a wrapper for block vectors. If there
    ! is only one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
!</description>

!<input>
    ! solution vector
    TYPE(t_vectorBlock), INTENT(IN)    :: ru

    ! scaling parameter
    REAL(DP), INTENT(IN)               :: dscale
!</input>

!<inputoutput>
    ! residual vector
    TYPE(t_vectorBlock), INTENT(INOUT) :: rres

    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)     :: rafcstab
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    IF (ru%nblocks .NE. 1 .OR. rres%nblocks .NE. 1) THEN

      CALL output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResBlockSymm')
      CALL sys_halt()

    ELSE

      CALL gfsc_buildResScalarSymm(ru%RvectorBlock(1), dscale,&
          rres%RvectorBlock(1), rafcstab)

    END IF
  END SUBROUTINE gfsc_buildResBlockSymm

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_buildResScalarSymm(ru, dscale, rres, rafcstab)

!<description>
    ! This subroutine assembles the residual vector and applies stabilisation
    ! by means of symmetric flux limiting for diffusion operators.
    !
    ! Yet, there is no publication available. This routine is based on
    ! private communication with D. Kuzmin.
    !
!</description>

!<input>
    ! solution vector
    TYPE(t_vectorScalar), INTENT(IN)    :: ru

    ! scaling parameter
    REAL(DP), INTENT(IN)                :: dscale
!</input>

!<inputoutput>
    ! residual vector
    TYPE(t_vectorScalar), INTENT(INOUT) :: rres

    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)      :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_VECIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge
    REAL(DP), DIMENSION(:,:), POINTER             :: p_DcoefficientsAtEdge
    REAL(DP), DIMENSION(:), POINTER :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    REAL(DP), DIMENSION(:), POINTER :: p_flux
    REAL(DP), DIMENSION(:), POINTER :: p_u,p_res
    

    ! Check if stabilisation is prepared
    IF (rafcstab%ctypeAFCstabilisation .NE. AFCSTAB_SYMMETRIC .OR.&
        IAND(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE) .EQ. 0    .OR.&
        IAND(rafcstab%iSpec, AFCSTAB_EDGEVALUES)    .EQ. 0) THEN
      CALL output_line('Stabilisation does not provide required structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResScalarSymm')
      CALL sys_halt()
    END IF

    ! Check if vectors are compatible
    CALL lsyssc_isVectorCompatible(ru, rres)

    ! Set pointers
    CALL afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    CALL afcstab_getbase_DcoeffsAtEdge(rafcstab,   p_DcoefficientsAtEdge)
    CALL lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
    CALL lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
    CALL lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
    CALL lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
    CALL lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
    CALL lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
    CALL lsyssc_getbase_double(rafcstab%RedgeVectors(1),  p_flux)
    CALL lsyssc_getbase_double(ru,   p_u)
    CALL lsyssc_getbase_double(rres, p_res)

    CALL do_symmetric_limit(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
        p_u, dscale, rafcstab%NEDGE, p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_flux, p_res)

    ! Set specifier
    rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_BOUNDS)
    rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_ANTIDIFFUSION)
    rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_LIMITER)
    rafcstab%iSpec = IOR(rafcstab%iSpec, AFCSTAB_FLUXES)

  CONTAINS
    
    ! Here, the working routine follow
    
    !**************************************************************
    ! Perform symmetric flux limiting
    
    SUBROUTINE do_symmetric_limit(IverticesAtEdge, DcoefficientsAtEdge,&
        u, dscale, NEDGE, pp, pm, qp, qm, rp, rm, flux, res)
      
      INTEGER(PREC_VECIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtEdge
      REAL(DP), DIMENSION(:,:), INTENT(IN)             :: DcoefficientsAtEdge
      REAL(DP), DIMENSION(:), INTENT(IN)               :: u
      REAL(DP), INTENT(IN)                             :: dscale
      INTEGER(PREC_MATIDX), INTENT(IN)                 :: NEDGE
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: pp,pm
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: qp,qm
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: rp,rm
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: flux
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: res
      
      INTEGER(PREC_MATIDX) :: iedge,ij
      INTEGER(PREC_VECIDX) :: i,j
      REAL(DP) :: d_ij,f_ij,s_ij,diff
      
      ! Clear nodal vectors
      CALL lalg_clearVectorDble(pp)
      CALL lalg_clearVectorDble(pm)
      CALL lalg_clearVectorDble(qp)
      CALL lalg_clearVectorDble(qm)
      
      ! Loop over edges
      DO iedge = 1, NEDGE
        
        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Determine coefficients
        d_ij = DcoefficientsAtEdge(1,iedge)
        s_ij = DcoefficientsAtEdge(2,iedge)
        
        ! Determine fluxes
        diff = u(i)-u(j); f_ij = d_ij*diff
        flux(iedge) = f_ij
        
        ! Sums of raw positive/negative fluxes
        pp(i) = pp(i)+MAX(0._DP, f_ij); pp(j) = pp(j)+MAX(0._DP,-f_ij)
        pm(i) = pm(i)+MIN(0._DP, f_ij); pm(j) = pm(j)+MIN(0._DP,-f_ij)
        
        ! Upper/lower bounds
        f_ij = -s_ij*diff
        qp(i) = qp(i)+MAX(0._DP, f_ij); qp(j) = qp(j)+MAX(0._DP,-f_ij)
        qm(i) = qm(i)+MIN(0._DP, f_ij); qm(j) = qm(j)+MIN(0._DP,-f_ij)
      END DO
      
      ! Apply the nodal limiter
      rp = afcstab_limit( pp, qp, 0._DP, 1._DP)
      rm = afcstab_limit(-pm,-qm, 0._DP, 1._DP)
      
      ! Apply limiter
      DO iedge = 1, NEDGE
        
        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Get precomputed raw antidiffusive flux
        f_ij = flux(iedge)
        
        IF (f_ij > 0.0_DP) THEN
          f_ij = dscale*MIN(rp(i), rm(j))*f_ij
        ELSE
          f_ij = dscale*MIN(rm(i), rp(j))*f_ij
        END IF
        
        ! Update the defect vector
        res(i) = res(i)+f_ij
        res(j) = res(j)-f_ij
      END DO
    END SUBROUTINE do_symmetric_limit
  END SUBROUTINE gfsc_buildResScalarSymm

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_buildConvJacobianBlock(RmatrixC, ru,&
      fcb_getVelocity, hstep, bclear, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the convective part
    ! of the discrete transport operator for a scalar convection equation.
    ! Note that this routine serves as a wrapper for block vectors. If there
    ! is only one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
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
    INCLUDE 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! Jacobian matrix
    TYPE(t_matrixScalar), INTENT(INOUT)            :: rmatrixJ
!</inputoutput>
!</subroutine>

    ! Check if block vector contains exactly one block
    IF (ru%nblocks .NE. 1) THEN
      
      CALL output_line('Solution vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvJacobianBlock')
      CALL sys_halt()

    ELSE
      
      CALL gfsc_buildConvJacobianScalar(RmatrixC, ru%RvectorBlock(1),&
          fcb_getVelocity, hstep, bclear, rmatrixJ)
      
    END IF
  END SUBROUTINE gfsc_buildConvJacobianBlock

   !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_buildConvJacobianScalar(RmatrixC, ru,&
      fcb_getVelocity, hstep, bclear, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the convective part
    ! of the discrete transport operator for a scalar convection equation.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: rmatrixC

    ! solution vector
    TYPE(t_vectorScalar), INTENT(IN)               :: ru
    
    ! perturbation parameter
    REAL(DP), INTENT(IN)                           :: hstep

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    LOGICAL, INTENT(IN)                            :: bclear

    ! callback functions to compute velocity
    INCLUDE 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! Jacobian matrix
    TYPE(t_matrixScalar), INTENT(INOUT)            :: rmatrixJ
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_Kld,p_Ksep
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_Kdiagonal
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER   :: p_Kcol
    REAL(DP), DIMENSION(:), POINTER               :: p_Cx,p_Cy,p_Cz,p_J,p_u
    INTEGER :: h_Ksep
    INTEGER :: idim,ndim
    
    ! Check if all matrices are compatible
    CALL lsyssc_isMatrixCompatible(ru, rmatrixJ, .FALSE.)
    ndim = SIZE(RmatrixC,1)
    DO idim = 1, ndim
      CALL lsyssc_isMatrixCompatible(RmatrixC(idim), rmatrixJ)
    END DO

    ! Clear matrix?
    IF (bclear) CALL lsyssc_clearMatrix(rmatrixJ)

    ! Set pointers
    CALL lsyssc_getbase_Kld(rmatrixJ,    p_Kld)
    CALL lsyssc_getbase_Kcol(rmatrixJ,   p_Kcol)
    CALL lsyssc_getbase_double(rmatrixJ, p_J)
    CALL lsyssc_getbase_double(ru,       p_u)

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
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvJacobianScalar')
      CALL sys_halt()
    END SELECT


    ! What kind of matrix are we?
    SELECT CASE(rmatrixJ%cmatrixFormat)
    CASE(LSYSSC_MATRIX7)

      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      CALL storage_copy(rmatrixJ%h_Kld, h_Ksep)
      CALL storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)

      ! Build Jacobian
      SELECT CASE(ndim)
      CASE (NDIM1D)
        CALL do_jacobianMat7_1D(p_Kld, p_Kcol, p_Ksep, rmatrixJ%NEQ, p_Cx, p_u, p_J)
      CASE (NDIM2D)
        CALL do_jacobianMat7_2D(p_Kld, p_Kcol, p_Ksep, rmatrixJ%NEQ, p_Cx, p_Cy, p_u, p_J)
      CASE (NDIM3D)
        CALL do_jacobianMat7_3D(p_Kld, p_Kcol, p_Ksep, rmatrixJ%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_J)
      END SELECT

      ! Release diagonal separator
      CALL storage_free(h_Ksep)

      
    CASE(LSYSSC_MATRIX9)
      
      ! Set pointers
      CALL lsyssc_getbase_Kdiagonal(rmatrixJ, p_Kdiagonal)
    
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      CALL storage_copy(rmatrixJ%h_Kld, h_Ksep)
      CALL storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)
      
      ! Build Jacobian
      SELECT CASE(ndim)
      CASE (NDIM1D)
        CALL do_jacobianMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, rmatrixJ%NEQ, p_Cx, p_u, p_J)
      CASE (NDIM2D)
        CALL do_jacobianMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, rmatrixJ%NEQ, p_Cx, p_Cy, p_u, p_J)
      CASE (NDIM3D)
        CALL do_jacobianMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, rmatrixJ%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_J)
      END SELECT
      
      ! Release diagonal separator
      CALL storage_free(h_Ksep)

    CASE DEFAULT
      CALL output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvJacobianScalar')
      CALL sys_halt()
    END SELECT

  CONTAINS

    ! Here, the working routine follow
    
    !**************************************************************
    ! Assemble Jacobian matrix for convective operator in 1D
    ! All matrices are stored in matrix format 7

    SUBROUTINE do_jacobianMat7_1D(Kld, Kcol, Ksep, NEQ, Cx, u, Jac)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: Jac
      
      REAL(DP), DIMENSION(NDIM1D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,l_ij,l_ji,xi_ij,xi_ji,diff
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j

      ! Loop over all rows I of Jacobian matrix
      DO i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients xi_ij and xi_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
                    
          ! solution difference u_j-u_i
          diff = u(j)-u(i)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients xi_ij^k which s defined as
          !   xi_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
          ! Compute perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i)+hstep, u(j), i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)
          l_ji = -v_ji(1)*Cx(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to xi_ij and xi_ji
          xi_ij = l_ij+d_ij; xi_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i)-hstep, u(j) ,i ,j ,v_ij ,v_ji)
          l_ij = -v_ij(1)*Cx(ij)
          l_ji = -v_ji(1)*Cx(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply the average of the perturbed coefficients
          Jac(ji) = Jac(ji)+0.5_DP*(xi_ji+l_ji+d_ij)
          Jac(jj) = Jac(jj)-0.5_DP*(xi_ji+l_ji+d_ij)
          
          ! Compute final coefficients xi_ij and xi_ji as the second
          ! order divided differences of the low-order coefficients
          xi_ij = 0.5_DP*(xi_ij-(l_ij+d_ij))/hstep
          xi_ji = 0.5_DP*(xi_ji-(l_ji+d_ij))/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+xi_ij*diff
          Jac(ji) = Jac(ji)-xi_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

          ! Compute perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i), u(j)+hstep, i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)
          l_ji = -v_ji(1)*Cx(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to xi_ij and xi_ji
          xi_ij = l_ij+d_ij; xi_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i), u(j)-hstep, i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)
          l_ji = -v_ji(1)*Cx(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply the average of the perturbed coefficients for J=K
          Jac(ij) = Jac(ij)+0.5_DP*(xi_ij+l_ij+d_ij)
          Jac(ii) = Jac(ii)-0.5_DP*(xi_ij+l_ij+d_ij)
          
          ! Compute final coefficients xi_ij and xi_ji as the second
          ! order divided differences of the low-order coefficients
          xi_ij = 0.5_DP*(xi_ij-(l_ij+d_ij))/hstep
          xi_ji = 0.5_DP*(xi_ji-(l_ji+d_ij))/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+xi_ij*diff
          Jac(jj) = Jac(jj)-xi_ji*diff
        END DO
      END DO
    END SUBROUTINE do_jacobianMat7_1D


    !**************************************************************
    ! Assemble Jacobian matrix for convective operator in 2D
    ! All matrices are stored in matrix format 7

    SUBROUTINE do_jacobianMat7_2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, Jac)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: Jac
      
      REAL(DP), DIMENSION(NDIM2D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,l_ij,l_ji,xi_ij,xi_ji,diff
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j

      ! Loop over all rows I of Jacobian matrix
      DO i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients xi_ij and xi_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
                    
          ! solution difference u_j-u_i
          diff = u(j)-u(i)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients xi_ij^k which s defined as
          !   xi_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
          ! Compute perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i)+hstep, u(j), i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to xi_ij and xi_ji
          xi_ij = l_ij+d_ij; xi_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i)-hstep, u(j) ,i ,j ,v_ij ,v_ji)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply the average of the perturbed coefficients
          Jac(ji) = Jac(ji)+0.5_DP*(xi_ji+l_ji+d_ij)
          Jac(jj) = Jac(jj)-0.5_DP*(xi_ji+l_ji+d_ij)
          
          ! Compute final coefficients xi_ij and xi_ji as the second
          ! order divided differences of the low-order coefficients
          xi_ij = 0.5_DP*(xi_ij-(l_ij+d_ij))/hstep
          xi_ji = 0.5_DP*(xi_ji-(l_ji+d_ij))/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+xi_ij*diff
          Jac(ji) = Jac(ji)-xi_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

          ! Compute perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i), u(j)+hstep, i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to xi_ij and xi_ji
          xi_ij = l_ij+d_ij; xi_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i), u(j)-hstep, i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply the average of the perturbed coefficients for J=K
          Jac(ij) = Jac(ij)+0.5_DP*(xi_ij+l_ij+d_ij)
          Jac(ii) = Jac(ii)-0.5_DP*(xi_ij+l_ij+d_ij)
          
          ! Compute final coefficients xi_ij and xi_ji as the second
          ! order divided differences of the low-order coefficients
          xi_ij = 0.5_DP*(xi_ij-(l_ij+d_ij))/hstep
          xi_ji = 0.5_DP*(xi_ji-(l_ji+d_ij))/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+xi_ij*diff
          Jac(jj) = Jac(jj)-xi_ji*diff
        END DO
      END DO
    END SUBROUTINE do_jacobianMat7_2D


    !**************************************************************
    ! Assemble Jacobian matrix for convective operator in 3D
    ! All matrices are stored in matrix format 7

    SUBROUTINE do_jacobianMat7_3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, Jac)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: Jac
      
      REAL(DP), DIMENSION(NDIM3D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,l_ij,l_ji,xi_ij,xi_ji,diff
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j

      ! Loop over all rows I of Jacobian matrix
      DO i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Ksep(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients xi_ij and xi_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
          
          ! solution difference u_j-u_i
          diff = u(j)-u(i)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients xi_ij^k which s defined as
          !   xi_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
          ! Compute perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i)+hstep, u(j), i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to xi_ij and xi_ji
          xi_ij = l_ij+d_ij; xi_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i)-hstep, u(j) ,i ,j ,v_ij ,v_ji)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply the average of the perturbed coefficients
          Jac(ji) = Jac(ji)+0.5_DP*(xi_ji+l_ji+d_ij)
          Jac(jj) = Jac(jj)-0.5_DP*(xi_ji+l_ji+d_ij)
          
          ! Compute final coefficients xi_ij and xi_ji as the second
          ! order divided differences of the low-order coefficients
          xi_ij = 0.5_DP*(xi_ij-(l_ij+d_ij))/hstep
          xi_ji = 0.5_DP*(xi_ji-(l_ji+d_ij))/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+xi_ij*diff
          Jac(ji) = Jac(ji)-xi_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

          ! Compute perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i), u(j)+hstep, i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to xi_ij and xi_ji
          xi_ij = l_ij+d_ij; xi_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i), u(j)-hstep, i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply the average of the perturbed coefficients for J=K
          Jac(ij) = Jac(ij)+0.5_DP*(xi_ij+l_ij+d_ij)
          Jac(ii) = Jac(ii)-0.5_DP*(xi_ij+l_ij+d_ij)
          
          ! Compute final coefficients xi_ij and xi_ji as the second
          ! order divided differences of the low-order coefficients
          xi_ij = 0.5_DP*(xi_ij-(l_ij+d_ij))/hstep
          xi_ji = 0.5_DP*(xi_ji-(l_ji+d_ij))/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+xi_ij*diff
          Jac(jj) = Jac(jj)-xi_ji*diff
        END DO
      END DO
    END SUBROUTINE do_jacobianMat7_3D
    
    
    !**************************************************************
    ! Assemble Jacobian matrix for convective operator in 1D
    ! All matrices are stored in matrix format 9

    SUBROUTINE do_jacobianMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, Jac)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: Jac
      
      REAL(DP), DIMENSION(NDIM1D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,l_ij,l_ji,xi_ij,xi_ji,diff
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j

      ! Loop over all rows I of Jacobian matrix
      DO i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients xi_ij and xi_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
          
          ! solution difference u_j-u_i
          diff = u(j)-u(i)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients xi_ij^k which s defined as
          !   xi_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
          ! Compute perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i)+hstep, u(j), i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)
          l_ji = -v_ji(1)*Cx(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to xi_ij and xi_ji
          xi_ij = l_ij+d_ij; xi_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i)-hstep, u(j), i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)
          l_ji = -v_ji(1)*Cx(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply the average of the perturbed coefficients
          Jac(ji) = Jac(ji)+0.5_DP*(xi_ji+l_ji+d_ij)
          Jac(jj) = Jac(jj)-0.5_DP*(xi_ji+l_ji+d_ij)
          
          ! Compute final coefficients xi_ij and xi_ji as the second
          ! order divided differences of the low-order coefficients
          xi_ij = 0.5_DP*(xi_ij-(l_ij+d_ij))/hstep
          xi_ji = 0.5_DP*(xi_ji-(l_ji+d_ij))/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+xi_ij*diff
          Jac(ji) = Jac(ji)-xi_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

          ! Compute perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i), u(j)+hstep, i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)
          l_ji = -v_ji(1)*Cx(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to xi_ij and xi_ji
          xi_ij = l_ij+d_ij; xi_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i), u(j)-hstep, i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)
          l_ji = -v_ji(1)*Cx(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply the average of the perturbed coefficients for J=K
          Jac(ij) = Jac(ij)+0.5_DP*(xi_ij+l_ij+d_ij)
          Jac(ii) = Jac(ii)-0.5_DP*(xi_ij+l_ij+d_ij)
          
          ! Compute final coefficients xi_ij and xi_ji as the second
          ! order divided differences of the low-order coefficients
          xi_ij = 0.5_DP*(xi_ij-(l_ij+d_ij))/hstep
          xi_ji = 0.5_DP*(xi_ji-(l_ji+d_ij))/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+xi_ij*diff
          Jac(jj) = Jac(jj)-xi_ji*diff
        END DO
      END DO
    END SUBROUTINE do_jacobianMat9_1D


    !**************************************************************
    ! Assemble Jacobian matrix for convective operator in 2D
    ! All matrices are stored in matrix format 9

    SUBROUTINE do_jacobianMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, Jac)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: Jac
      
      REAL(DP), DIMENSION(NDIM2D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,l_ij,l_ji,xi_ij,xi_ji,diff
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j

      ! Loop over all rows I of Jacobian matrix
      DO i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients xi_ij and xi_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
          
          ! solution difference u_j-u_i
          diff = u(j)-u(i)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients xi_ij^k which s defined as
          !   xi_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
          ! Compute perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i)+hstep, u(j), i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to xi_ij and xi_ji
          xi_ij = l_ij+d_ij; xi_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i)-hstep, u(j), i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply the average of the perturbed coefficients
          Jac(ji) = Jac(ji)+0.5_DP*(xi_ji+l_ji+d_ij)
          Jac(jj) = Jac(jj)-0.5_DP*(xi_ji+l_ji+d_ij)
          
          ! Compute final coefficients xi_ij and xi_ji as the second
          ! order divided differences of the low-order coefficients
          xi_ij = 0.5_DP*(xi_ij-(l_ij+d_ij))/hstep
          xi_ji = 0.5_DP*(xi_ji-(l_ji+d_ij))/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+xi_ij*diff
          Jac(ji) = Jac(ji)-xi_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

          ! Compute perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i), u(j)+hstep, i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to xi_ij and xi_ji
          xi_ij = l_ij+d_ij; xi_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i), u(j)-hstep, i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply the average of the perturbed coefficients for J=K
          Jac(ij) = Jac(ij)+0.5_DP*(xi_ij+l_ij+d_ij)
          Jac(ii) = Jac(ii)-0.5_DP*(xi_ij+l_ij+d_ij)
          
          ! Compute final coefficients xi_ij and xi_ji as the second
          ! order divided differences of the low-order coefficients
          xi_ij = 0.5_DP*(xi_ij-(l_ij+d_ij))/hstep
          xi_ji = 0.5_DP*(xi_ji-(l_ji+d_ij))/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+xi_ij*diff
          Jac(jj) = Jac(jj)-xi_ji*diff
        END DO
      END DO
    END SUBROUTINE do_jacobianMat9_2D


    !**************************************************************
    ! Assemble Jacobian matrix for convective operator in 3D
    ! All matrices are stored in matrix format 9

    SUBROUTINE do_jacobianMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, Jac)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz,u
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: Jac
      
      REAL(DP), DIMENSION(NDIM3D) :: v_ij,v_ji
      REAL(DP)                    :: d_ij,l_ij,l_ji,xi_ij,xi_ji,diff
      INTEGER(PREC_MATIDX)        :: ii,ij,ji,jj
      INTEGER(PREC_VECIDX)        :: i,j

      ! Loop over all rows I of Jacobian matrix
      DO i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        DO ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients xi_ij and xi_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
          
          ! solution difference u_j-u_i
          diff = u(j)-u(i)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients xi_ij^k which s defined as
          !   xi_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
          ! Compute perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i)+hstep, u(j), i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to xi_ij and xi_ji
          xi_ij = l_ij+d_ij; xi_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i)-hstep, u(j), i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply the average of the perturbed coefficients
          Jac(ji) = Jac(ji)+0.5_DP*(xi_ji+l_ji+d_ij)
          Jac(jj) = Jac(jj)-0.5_DP*(xi_ji+l_ji+d_ij)
          
          ! Compute final coefficients xi_ij and xi_ji as the second
          ! order divided differences of the low-order coefficients
          xi_ij = 0.5_DP*(xi_ij-(l_ij+d_ij))/hstep
          xi_ji = 0.5_DP*(xi_ji-(l_ji+d_ij))/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+xi_ij*diff
          Jac(ji) = Jac(ji)-xi_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

          ! Compute perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i), u(j)+hstep, i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to xi_ij and xi_ji
          xi_ij = l_ij+d_ij; xi_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          CALL fcb_getVelocity(u(i), u(j)-hstep, i, j, v_ij, v_ji)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          d_ij =  MAX(-l_ij, 0._DP, -l_ji)
          
          ! Apply the average of the perturbed coefficients for J=K
          Jac(ij) = Jac(ij)+0.5_DP*(xi_ij+l_ij+d_ij)
          Jac(ii) = Jac(ii)-0.5_DP*(xi_ij+l_ij+d_ij)
          
          ! Compute final coefficients xi_ij and xi_ji as the second
          ! order divided differences of the low-order coefficients
          xi_ij = 0.5_DP*(xi_ij-(l_ij+d_ij))/hstep
          xi_ji = 0.5_DP*(xi_ji-(l_ji+d_ij))/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+xi_ij*diff
          Jac(jj) = Jac(jj)-xi_ji*diff
        END DO
      END DO
    END SUBROUTINE do_jacobianMat9_3D
  END SUBROUTINE gfsc_buildConvJacobianScalar

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_buildStabJacLinearBlock_FCT(rmatrixMC, ru, &
      theta, tstep, hstep, bclear, rafcstab, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation part
    ! of the discrete transport operator for a scalar convection equation.
    ! Note that the velocity is assumed to be linear.
    ! Note that this routine serves as a wrapper for block vectors. If there
    ! is only one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
!</description>

!<input>
    ! consistent mass matrix
    TYPE(t_matrixScalar), INTENT(IN)    :: rmatrixMC

    ! solution vector
    TYPE(t_vectorBlock), INTENT(IN)     :: ru

    ! implicitness parameter
    REAL(DP), INTENT(IN)                :: theta

    ! time step size
    REAL(DP), INTENT(IN)                :: tstep

    ! perturbation parameter
    REAL(DP), INTENT(IN)                :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    LOGICAL, INTENT(IN)                 :: bclear
!</input>

!<inputoutput>
    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)      :: rafcstab

    ! Jacobian matrix
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrixJ   
!</inputoutput>
!</subroutine>

    IF (ru%nblocks  .NE. 1) THEN

      CALL output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacLinearBlock_FCT')
      CALL sys_halt()

    ELSE

      CALL gfsc_buildStabJacLinearScalar_FCT(rmatrixMC, ru%RvectorBlock(1),&
          theta, tstep, hstep, bclear, rafcstab, rmatrixJ)

    END IF
  END SUBROUTINE gfsc_buildStabJacLinearBlock_FCT

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_buildStabJacLinearScalar_FCT(rmatrixMC, ru, &
      theta, tstep, hstep, bclear, rafcstab, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation part
    ! of the discrete transport operator for a scalar convection equation.
    ! Note that the velocity is assumed to be linear.
!</description>

!<input>
    ! consistent mass matrix
    TYPE(t_matrixScalar), INTENT(IN)    :: rmatrixMC

    ! solution vector
    TYPE(t_vectorScalar), INTENT(IN)    :: ru

    ! implicitness parameter
    REAL(DP), INTENT(IN)                :: theta

    ! time step size
    REAL(DP), INTENT(IN)                :: tstep

    ! perturbation parameter
    REAL(DP), INTENT(IN)                :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    LOGICAL, INTENT(IN)                 :: bclear
!</input>

!<inputoutput>
    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)      :: rafcstab

    ! Jacobian matrix
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrixJ   
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_VECIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_Kld
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_Kdiagonal
    REAL(DP), DIMENSION(:,:), POINTER             :: p_DcoefficientsAtEdge
    REAL(DP), DIMENSION(:), POINTER               :: p_flux,p_flux0
    REAL(DP), DIMENSION(:), POINTER               :: p_u
    REAL(DP), DIMENSION(:), POINTER               :: p_MC,p_Jac
    
    
    ! Check if stabilisation is prepared
    IF (rafcstab%ctypeAFCstabilisation .NE. AFCSTAB_FEMFCT .OR.&
        IAND(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE) .EQ. 0 .OR.&
        IAND(rafcstab%iSpec, AFCSTAB_EDGEVALUES)    .EQ. 0 .OR.&
        IAND(rafcstab%iSpec, AFCSTAB_FLUXES)        .EQ. 0) THEN
      CALL output_line('Stabilisation does not provide required structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacLinearScalar_FCT')
      CALL sys_halt()
    END IF
    
    ! Check if matrices/vectors are compatible
    CALL lsyssc_isMatrixCompatible(ru, rmatrixMC, .FALSE.)
    CALL lsyssc_isMatrixCompatible(ru, rmatrixJ,  .FALSE.)
    CALL gfsc_isMatrixCompatible(rafcstab, rmatrixMC)
    CALL gfsc_isVectorCompatible(rafcstab, ru)
    
    ! Clear matrix?
    IF (bclear) CALL lsyssc_clearMatrix(rmatrixJ)

    ! Set pointers
    CALL afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    CALL afcstab_getbase_DcoeffsAtEdge(rafcstab,   p_DcoefficientsAtEdge)
    CALL lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_flux)
    CALL lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_flux0)
    CALL lsyssc_getbase_double(rmatrixMC, p_MC)
    CALL lsyssc_getbase_double(rmatrixJ,  p_Jac)
    CALL lsyssc_getbase_double(ru,        p_u)
    
    ! What kind of matrix are we?
    SELECT CASE(rmatrixJ%cmatrixFormat)
    CASE(LSYSSC_MATRIX7)
      CALL lsyssc_getbase_Kld(rmatrixJ, p_Kld)
      CALL do_femfct(p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kld, p_MC,&
          p_u, p_flux, p_flux0, theta, tstep, hstep, rafcstab%NEDGE,&
          (rafcstab%imass .EQ. AFCSTAB_CONSISTENTMASS), p_Jac)
      
    CASE(LSYSSC_MATRIX9)
      CALL lsyssc_getbase_Kdiagonal(rmatrixJ, p_Kdiagonal)
      CALL do_femfct(p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kdiagonal, p_MC,&
          p_u, p_flux, p_flux0, theta, tstep, hstep, rafcstab%NEDGE,&
          (rafcstab%imass .EQ. AFCSTAB_CONSISTENTMASS), p_Jac)
      
    CASE DEFAULT
      CALL output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacLinearScalar_FCT')
      CALL sys_halt()
    END SELECT
    
  CONTAINS
    
    ! Here, the working routine follow

    !**************************************************************
    ! Assemble the Jacobian matrix for semi-implicit FEM-FCT

    SUBROUTINE do_femfct(IverticesAtEdge, DcoefficientsAtEdge, Kdiagonal,&
        MC, u, flux, flux0, theta, tstep, hstep, NEDGE, bmass, Jac)

      INTEGER(PREC_VECIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtEdge
      REAL(DP), DIMENSION(:,:), INTENT(IN)             :: DcoefficientsAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)   :: Kdiagonal
      REAL(DP), DIMENSION(:), INTENT(IN)               :: MC
      REAL(DP), DIMENSION(:), INTENT(IN)               :: flux,flux0
      REAL(DP), DIMENSION(:), INTENT(IN)               :: u
      REAL(DP), INTENT(IN)                             :: theta,tstep,hstep
      INTEGER(PREC_MATIDX), INTENT(IN)                 :: NEDGE
      LOGICAL, INTENT(IN)                              :: bmass
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: Jac

      ! local variables
      INTEGER(PREC_MATIDX) :: iedge,ij,ji,ii,jj
      INTEGER(PREC_VECIDX) :: i,j
      REAL(DP)             :: f_i,f_j,f_ij,d_ij,a_ij
      REAL(DP)             :: diff,diff_i,diff_j

      ! Should we apply the consistent mass matrix?
      IF (bmass) THEN

        ! Loop over all edges
        DO iedge = 1, NEDGE
          
          ! Determine vertex numbers
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)
          
          ! Determine diagonal indices
          ii = Kdiagonal(i); jj = Kdiagonal(j)
          
          ! Determine coefficients
          d_ij = DcoefficientsAtEdge(1,iedge)
          a_ij = MC(ij)+theta*tstep*d_ij
          
          ! Compute solution difference
          diff = u(i)-u(j)
          
          ! Compute perturbed solution differences
          diff_i = diff+hstep
          diff_j = diff-hstep
          
          ! Compute limited antidiffusive flux f(u_ij+h*e_i)
          f_i = a_ij*diff_i+flux0(iedge)
          IF (f_i > 0.0_DP) THEN
            f_i = MIN(f_i, MAX(flux(iedge), 0._DP))
          ELSE
            f_i = MAX(f_i, MIN(flux(iedge), 0._DP))
          END IF
          
          ! Compute limited antidiffusive flux f(u_ij-h*e_j)
          f_j = a_ij*diff_j+flux0(iedge)
          IF (f_j > 0.0_DP) THEN
            f_j = MIN(f_j, MAX(flux(iedge), 0._DP))
          ELSE
            f_j = MAX(f_j, MIN(flux(iedge), 0._DP))
          END IF
          
          ! Compute divided differences of fluxes
          f_ij = 0.5_DP*(f_i-f_j)/hstep
          
          ! Apply i-th column
          Jac(ii) = Jac(ii)-f_ij
          Jac(ji) = Jac(ji)+f_ij
          
          ! Apply j-th column
          Jac(ij) = Jac(ij)+f_ij
          Jac(jj) = Jac(jj)-f_ij
        END DO

      ELSE

        ! Loop over all edges
        DO iedge = 1, NEDGE
          
          ! Determine vertex numbers
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)
          
          ! Determine diagonal indices
          ii = Kdiagonal(i); jj = Kdiagonal(j)
          
          ! Determine coefficients
          d_ij = DcoefficientsAtEdge(1,iedge)
          a_ij = theta*tstep*d_ij
          
          ! Compute solution difference
          diff = u(i)-u(j)
          
          ! Compute perturbed solution differences
          diff_i = diff+hstep
          diff_j = diff-hstep
          
          ! Compute limited antidiffusive flux f(u_ij+h*e_i) 
          f_i = a_ij*diff_i+flux0(iedge)
          IF (f_i > 0.0_DP) THEN
            f_i = MIN(f_i, MAX(flux(iedge), 0._DP))
          ELSE
            f_i = MAX(f_i, MIN(flux(iedge), 0._DP))
          END IF
          
          ! Compute limited antidiffusive flux f(u_ij-h*e_j)
          f_j = a_ij*diff_j+flux0(iedge)
          IF (f_j > 0.0_DP) THEN
            f_j = MIN(f_j, MAX(flux(iedge), 0._DP))
          ELSE
            f_j = MAX(f_j, MIN(flux(iedge), 0._DP))
          END IF
          
          ! Compute divided differences of fluxes
          f_ij = 0.5_DP*(f_i-f_j)/hstep
          
          ! Apply i-th column
          Jac(ii) = Jac(ii)-f_ij
          Jac(ji) = Jac(ji)+f_ij
          
          ! Apply j-th column
          Jac(ij) = Jac(ij)+f_ij
          Jac(jj) = Jac(jj)-f_ij
        END DO
      END IF     
    END SUBROUTINE do_femfct
  END SUBROUTINE gfsc_buildStabJacLinearScalar_FCT

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_buildStabJacLinearBlock_GPTVD(rmatrixMC, ru, ru0,&
      theta, tstep, hstep, bclear, rafcstab, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation part
    ! of the discrete transport operator for a scalar convection equation.
    ! Note that the velocity is assumed to be linear.
    ! Note that this routine serves as a wrapper for block vectors. If there
    ! is only one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
!</description>

!<input>
    ! consistent mass matrix
    TYPE(t_matrixScalar), INTENT(IN)    :: rmatrixMC

    ! solution vector
    TYPE(t_vectorBlock), INTENT(IN)     :: ru

    ! initial solution vector
    TYPE(t_vectorBlock), INTENT(IN)     :: ru0

    ! implicitness parameter
    REAL(DP), INTENT(IN)                :: theta

    ! time step size
    REAL(DP), INTENT(IN)                :: tstep

    ! perturbation parameter
    REAL(DP), INTENT(IN)                :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    LOGICAL, INTENT(IN)                 :: bclear
!</input>

!<inputoutput>
    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)      :: rafcstab

    ! Jacobian matrix
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrixJ   
!</inputoutput>
!</subroutine>

    IF (ru%nblocks  .NE. 1 .OR.&
        ru0%nblocks .NE. 1) THEN

      CALL output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacLinearBlock_GPTVD')
      CALL sys_halt()

    ELSE

      CALL gfsc_buildStabJacLinearScalar_GPTVD(rmatrixMC,&
          ru%RvectorBlock(1), ru0%RvectorBlock(1),&
          theta, tstep, hstep, bclear, rafcstab, rmatrixJ)

    END IF
  END SUBROUTINE gfsc_buildStabJacLinearBlock_GPTVD

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_buildStabJacLinearScalar_GPTVD(rmatrixMC, ru, ru0,&
      theta, tstep, hstep, bclear, rafcstab, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation part
    ! of the discrete transport operator for a scalar convection equation.
    ! Note that the velocity is assumed to be linear.
!</description>

!<input>
    ! consistent mass matrix
    TYPE(t_matrixScalar), INTENT(IN)    :: rmatrixMC

    ! solution vector
    TYPE(t_vectorScalar), INTENT(IN)    :: ru

    ! initial solution vector
    TYPE(t_vectorScalar), INTENT(IN)    :: ru0

    ! implicitness parameter
    REAL(DP), INTENT(IN)                :: theta

    ! time step size
    REAL(DP), INTENT(IN)                :: tstep

    ! perturbation parameter
    REAL(DP), INTENT(IN)                :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    LOGICAL, INTENT(IN)                 :: bclear
!</input>

!<inputoutput>
    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)      :: rafcstab

    ! Jacobian matrix
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrixJ   
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_VECIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER   :: p_IsuperdiagonalEdgesIdx
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_IsubdiagonalEdges
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_IsubdiagonalEdgesIdx
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_Kld,p_Ksep
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_Kdiagonal
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER   :: p_Kcol
    REAL(DP), DIMENSION(:,:), POINTER             :: p_DcoefficientsAtEdge
    REAL(DP), DIMENSION(:), POINTER               :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    REAL(DP), DIMENSION(:), POINTER               :: p_flux,p_flux0
    REAL(DP), DIMENSION(:), POINTER               :: p_u,p_u0
    REAL(DP), DIMENSION(:), POINTER               :: p_MC,p_Jac
    INTEGER :: h_Ksep
    LOGICAL :: bextend

    ! Clear matrix?
    IF (bclear) CALL lsyssc_clearMatrix(rmatrixJ)

    ! What kind of stabilisation are we?
    SELECT CASE(rafcstab%ctypeAFCstabilisation)

    CASE(AFCSTAB_FEMTVD)
      
      ! Check if stabilisation is prepared
      IF (IAND(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)   .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_EDGEORIENTATION) .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_EDGEVALUES)      .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_ANTIDIFFUSION)   .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_BOUNDS)          .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_FLUXES)          .EQ. 0) THEN
        CALL output_line('Stabilisation does not provide required structures',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacobianScalar_GPTVD')
        CALL sys_halt()
      END IF
      
      ! Check if matrices/vectors are compatible
      CALL lsyssc_isMatrixCompatible(ru, rmatrixJ, .FALSE.)
      CALL gfsc_isVectorCompatible(rafcstab, ru)

      ! Check if subdiagonal edges need to be generated
      IF (IAND(rafcstab%iSpec, AFCSTAB_SUBDIAGONALEDGES) .EQ. 0)&
          CALL afcstab_generateSubdiagEdges(rafcstab)

      ! Set pointers
      CALL afcstab_getbase_IverticesAtEdge(rafcstab,  p_IverticesAtEdge)
      CALL afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
      CALL afcstab_getbase_IsubdiagEdge(rafcstab,    p_IsubdiagonalEdges)
      CALL afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagonalEdgesIdx)
      CALL afcstab_getbase_DcoeffsAtEdge(rafcstab,    p_DcoefficientsAtEdge)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
      CALL lsyssc_getbase_double(rafcstab%RedgeVectors(1),  p_flux)
      CALL lsyssc_getbase_Kcol(rmatrixJ,   p_Kcol)
      CALL lsyssc_getbase_double(rmatrixJ, p_Jac)
      CALL lsyssc_getbase_double(ru,       p_u)

      ! What kind of matrix format are we?
      SELECT CASE(rmatrixJ%cmatrixFormat)
      CASE(LSYSSC_MATRIX7)
        
        ! Set pointers
        CALL lsyssc_getbase_Kld(rmatrixJ, p_Kld)

        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        CALL storage_copy(rmatrixJ%h_Kld, h_Ksep)
        CALL storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)
        CALL lalg_vectorAddScalarInt(p_Ksep, 1)

        ! Assembled extended Jacobian matrix
        bextend = (rafcstab%iextendedJacobian .NE. 0)
        CALL do_femtvd(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
            p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
            p_Kld, p_Kcol, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
            theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)
        
        ! Free storage
        CALL storage_free(h_Ksep)

      CASE(LSYSSC_MATRIX9)

        ! Set pointers
        CALL lsyssc_getbase_Kdiagonal(rmatrixJ, p_Kdiagonal)
        
        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        CALL storage_copy(rmatrixJ%h_Kld, h_Ksep)
        CALL storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)

        ! Assembled extended Jacobian matrix
        bextend = (rafcstab%iextendedJacobian .NE. 0)
        CALL do_femtvd(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
            p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
            p_Kdiagonal, p_Kcol, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
            theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)
        
        ! Free storage
        CALL storage_free(h_Ksep)
        
      CASE DEFAULT
        CALL output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacobianScalar_GPTVD')
        CALL sys_halt()
      END SELECT
      

    CASE(AFCSTAB_FEMGP)

       ! Check if stabilisation is prepared
      IF (IAND(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)   .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_EDGEORIENTATION) .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_EDGEVALUES)      .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_ANTIDIFFUSION)   .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_BOUNDS)          .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_FLUXES)          .EQ. 0) THEN
        CALL output_line('Stabilisation does not provide required structures',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacobianScalar')
        CALL sys_halt()
      END IF

      ! Check if matrices/vectors are compatible
      CALL lsyssc_isMatrixCompatible(ru, rmatrixMC, .FALSE.)
      CALL lsyssc_isMatrixCompatible(ru, rmatrixJ,  .FALSE.)
      CALL gfsc_isMatrixCompatible(rafcstab, rmatrixMC)
      CALL gfsc_isVectorCompatible(rafcstab, ru)
      CALL lsyssc_isVectorCompatible(ru, ru0)
      
      ! Check if subdiagonal edges need to be generated
      IF (IAND(rafcstab%iSpec, AFCSTAB_SUBDIAGONALEDGES) .EQ. 0)&
          CALL afcstab_generateSubdiagEdges(rafcstab)

      ! Set pointers
      CALL afcstab_getbase_IverticesAtEdge(rafcstab,  p_IverticesAtEdge)
      CALL afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
      CALL afcstab_getbase_IsubdiagEdge(rafcstab,    p_IsubdiagonalEdges)
      CALL afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagonalEdgesIdx)
      CALL afcstab_getbase_DcoeffsAtEdge(rafcstab,    p_DcoefficientsAtEdge)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
      CALL lsyssc_getbase_double(rafcstab%RedgeVectors(1),  p_flux)
      CALL lsyssc_getbase_double(rafcstab%RedgeVectors(2),  p_flux0)
      CALL lsyssc_getbase_Kcol(rmatrixJ,   p_Kcol)
      CALL lsyssc_getbase_double(rmatrixJ, p_Jac)
      CALL lsyssc_getbase_double(ru,       p_u)
      CALL lsyssc_getbase_double(ru0,      p_u0)

      ! What kind of matrix format are we?
      SELECT CASE(rmatrixJ%cmatrixFormat)
      CASE(LSYSSC_MATRIX7)
        
        ! Set pointers
        CALL lsyssc_getbase_Kld(rmatrixJ, p_Kld)

        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        CALL storage_copy(rmatrixJ%h_Kld, h_Ksep)
        CALL storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)
        CALL lalg_vectorAddScalarInt(p_Ksep, 1)

        ! Assembled extended Jacobian matrix
        bextend = (rafcstab%iextendedJacobian .NE. 0)
        IF (rafcstab%imass .EQ. AFCSTAB_CONSISTENTMASS) THEN
          ! Set pointers
          CALL lsyssc_getbase_double(rmatrixMC, p_MC)
          CALL lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
          CALL lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
          
          CALL do_femgp(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
              p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
              p_Kld, p_Kcol, p_MC, p_u, p_u0, p_flux, p_flux0, p_pp, p_pm, p_qp, p_qm,&
              p_rp, p_rm, theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
              rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)
          
        ELSE
          CALL do_femtvd(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
              p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
              p_Kld, p_Kcol, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
              theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
              rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)
        END IF

        ! Free storage
        CALL storage_free(h_Ksep)

        CASE(LSYSSC_MATRIX9)

        ! Set pointers
        CALL lsyssc_getbase_Kdiagonal(rmatrixJ, p_Kdiagonal)
        
        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        CALL storage_copy(rmatrixJ%h_Kld, h_Ksep)
        CALL storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)

        ! Assembled extended Jacobian matrix
        bextend = (rafcstab%iextendedJacobian .NE. 0)
        IF (rafcstab%imass .EQ. AFCSTAB_CONSISTENTMASS) THEN
          ! Set pointers
          CALL lsyssc_getbase_double(rmatrixMC, p_MC)
          CALL lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
          CALL lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)

          CALL do_femgp(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
              p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
              p_Kdiagonal, p_Kcol, p_MC, p_u, p_u0, p_flux, p_flux0,&
              p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, theta, tstep, hstep,&
              rafcstab%NEQ, rafcstab%NEDGE, rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)
        ELSE
          CALL do_femtvd(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
              p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
              p_Kdiagonal, p_Kcol, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
              theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
              rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)
        END IF

        ! Free storage
        CALL storage_free(h_Ksep)
        
      CASE DEFAULT
        CALL output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacobianScalar_GPTVD')
        CALL sys_halt()
      END SELECT
      
    CASE DEFAULT
      CALL output_line('Invalid type of AFC stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacLinearScalar_GPTVD')
      CALL sys_halt()
    END SELECT

  CONTAINS

    ! Here, the working routine follow
    
    !**************************************************************    
    ! Adjust the diagonal separator.
    ! The separator is initialied by the column separator (increased
    ! by one if the is necessary for matrix format 7).
    ! Based on the matric structure given by Kld/Kcol, the separator
    ! is "moved" to the given column "k". For efficiency reasons, only
    ! those entries are considered which are present in column "k".
    SUBROUTINE do_adjustKsep(Kld, Kcol, k, Ksep)
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: k
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      
      INTEGER(PREC_MATIDX) :: ild,isep
      INTEGER(PREC_VECIDX) :: l
      INTEGER :: iloc

      ! Loop over all entries of the k-th row
      DO ild = Kld(k), Kld(k+1)-1
        
        ! Get the column number and the position of the separator
        l = Kcol(ild); isep = Ksep(l)

        ! If the separator does not point to the k-th column
        ! it must be adjusted accordingly
        IF (Kcol(isep) < k) Ksep(l) = Ksep(l)+1
      END DO
    END SUBROUTINE do_adjustKsep


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-TVD
    SUBROUTINE do_femtvd(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
        IsubdiagonalEdgesIdx, IsubdiagonalEdges, DcoefficientsAtEdge,&
        Kld, Kcol, u, flux, pp, pm, qp, qm, theta, tstep, hstep,&
        NEQ, NEDGE, NNVEDGE, bextend, Ksep, Jac)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN)  :: IverticesAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: IsubdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: IsubdiagonalEdges
      REAL(DP), DIMENSION(:,:), INTENT(IN)              :: DcoefficientsAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      REAL(DP), DIMENSION(:), INTENT(IN)                :: u
      REAL(DP), DIMENSION(:), INTENT(IN)                :: flux
      REAL(DP), DIMENSION(:), INTENT(IN)                :: pp,pm,qp,qm
      REAL(DP), INTENT(IN)                              :: theta,tstep,hstep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER(PREC_MATIDX), INTENT(IN)                  :: NEDGE
      INTEGER, INTENT(IN)                               :: NNVEDGE
      LOGICAL, INTENT(IN)                               :: bextend
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: Jac
      
      ! local variables
      INTEGER(PREC_VECIDX), DIMENSION(5,NNVEDGE) :: Kloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: pploc,pmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: qploc,qmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: rploc,rmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: fluxloc
           
      INTEGER(PREC_MATIDX) :: ild,iedge
      INTEGER(PREC_VECIDX) :: k,l
      INTEGER :: iloc,nloc

      ! Loop over all columns of the Jacobian matrix
      DO k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        pploc(:,0) = 0; pmloc(:,0) = 0
        qploc(:,0) = 0; qmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        DO ild = IsubdiagonalEdgesIdx(k), IsubdiagonalEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagonalEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Update local coefficients
          CALL do_femtvd_update(IverticesAtEdge, DcoefficientsAtEdge,&
              u, pp, pm, qp, qm, tstep, hstep, iedge, iloc, k,&
              pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        END DO

        ! Loop over all superdiagonal edges
        DO iedge = IsuperdiagonalEdgesIdx(k), IsuperdiagonalEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                    
          ! Update local coefficients
          CALL do_femtvd_update(IverticesAtEdge, DcoefficientsAtEdge,&
              u, pp, pm, qp, qm, tstep, hstep, iedge, iloc, k,&
              pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        END DO
        
        ! Save total number of local neighbors
        nloc = iloc
        
        ! Adjust the diagonal separator
        CALL do_adjustKsep(Kld, Kcol, k, Ksep)

        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        rploc(:,0:nloc) = afcstab_limit( pploc(:,0:nloc), qploc(:,0:nloc), 0._DP, 1._DP)
        rmloc(:,0:nloc) = afcstab_limit(-pmloc(:,0:nloc),-qmloc(:,0:nloc), 0._DP, 1._DP)

        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        DO iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          DO ild = IsubdiagonalEdgesIdx(l), IsubdiagonalEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagonalEdges(ild)
            
            CALL do_femtvd_assemble(IverticesAtEdge, Kld, Kcol, flux, Kloc, rploc, rmloc,&
                fluxloc, hstep, iedge, iloc, k, l, bextend, Ksep, Jac)
          END DO

          ! Loop over all superdiagonal edges
          DO iedge = IsuperdiagonalEdgesIdx(l), IsuperdiagonalEdgesIdx(l+1)-1
            
            CALL do_femtvd_assemble(IverticesAtEdge, Kld, Kcol, flux, Kloc, rploc, rmloc,&
                fluxloc, hstep, iedge, iloc, k, l, bextend, Ksep, Jac)
          END DO
        END DO
      END DO   ! end-of k-loop
    END SUBROUTINE do_femtvd


    !**************************************************************
    ! Update the local coefficients for FEM-TVD
    SUBROUTINE do_femtvd_update(IverticesAtEdge, DcoefficientsAtEdge,&
        u, pp, pm, qp, qm, tstep, hstep, iedge, iloc, k,&
        pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
      
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN)     :: IverticesAtEdge
      REAL(DP), DIMENSION(:,:), INTENT(IN)                 :: DcoefficientsAtEdge
      REAL(DP), DIMENSION(:), INTENT(IN)                   :: u
      REAL(DP), DIMENSION(:), INTENT(IN)                   :: pp,pm,qp,qm
      REAL(DP), INTENT(IN)                                 :: tstep,hstep
      INTEGER(PREC_MATIDX), INTENT(IN)                     :: iedge
      INTEGER, INTENT(IN)                                  :: iloc
      INTEGER(PREC_VECIDX), INTENT(IN)                     :: k

      ! We actually know, that all local quantities start at index zero
      REAL(DP), DIMENSION(:,0:), INTENT(INOUT)             :: pploc,pmloc
      REAL(DP), DIMENSION(:,0:), INTENT(INOUT)             :: qploc,qmloc
      REAL(DP), DIMENSION(:,0:), INTENT(INOUT)             :: fluxloc
      INTEGER(PREC_VECIDX), DIMENSION(:,:), INTENT(INOUT)  :: Kloc

      ! local variables
      INTEGER(PREC_VECIDX) :: i,j
      REAL(DP) :: d_ij,f_ij,l_ij,l_ji,diff,hstep_ik,hstep_jk,dsign
      INTEGER  :: iperturb

      
      ! Determine indices. Obviously, either i or j must be equal
      ! to k. Otherwise, the edge ij would not be present in the
      ! list of incident edges for node k.
      i = IverticesAtEdge(1,iedge)
      j = IverticesAtEdge(2,iedge)
      
      !------------------------------------------------------------
      ! (1) unperturbed values: Retrieve the global Ps and Qs and
      !     copy their content to the local ones. Moreover,
      !     "eliminate" the contribution of the edge IJ for the
      !     unperturbed solution values u_i and u_j.
      !------------------------------------------------------------
      ! Determine coefficients
      d_ij = DcoefficientsAtEdge(1,iedge)
      l_ij = DcoefficientsAtEdge(2,iedge)
      l_ji = DcoefficientsAtEdge(3,iedge)
      
      ! Determine prelimited antidiffusive flux
      diff = tstep*(u(i)-u(j))
      f_ij = MIN(d_ij, l_ji)*diff
      
      IF (i .EQ. k) THEN

        ! Store global node number of the opposite node
        Kloc(1,iloc) = j

        ! Compute signed perturbation parameters
        hstep_ik = hstep; hstep_jk = 0._DP
        
        ! Update nodal coefficients for vertex j (!) which is the downwind node
        pploc(:,iloc) = pp(j)
        pmloc(:,iloc) = pm(j)
        qploc(:,iloc) = qp(j)-MAX(0._DP, f_ij)
        qmloc(:,iloc) = qm(j)-MIN(0._DP, f_ij)

      ELSE

        ! Store global node number of the opposite node
        Kloc(1,iloc) = i

        ! Compute signed perturbation parameters
        hstep_ik = 0._DP; hstep_jk = hstep
        
        ! Update nodal coefficients for vertex i (!) which is the upwind node
        pploc(:,iloc) = pp(i)-MAX(0._DP, f_ij)
        pmloc(:,iloc) = pm(i)-MIN(0._DP, f_ij)
        qploc(:,iloc) = qp(i)-MAX(0._DP,-f_ij)
        qmloc(:,iloc) = qm(i)-MIN(0._DP,-f_ij)
      END IF

      !------------------------------------------------------------
      ! (2) perturbed values: Now, the local Ps and Qs still
      !     require the contribution of the perturbed solution
      !     values u +/- h*e_k, whereby e_k denotes the k-th unit
      !     vector and h stands for the perturbation step length
      !------------------------------------------------------------

      !------------------------------------------------------------
      ! (3) perform the perturbation for "+/-h*e_k"
      !------------------------------------------------------------
        
      DO iperturb = 1, 2
        
        ! Compute correct sign of perturbation
        dsign = -2*iperturb+3
        
        ! The velocity is assumed to be linear.
        ! Hence the orientation convention for the edge ij is preserved
        
        ! Save oriented node numbers
        Kloc(2*iperturb:2*iperturb+1,iloc) = (/i,j/)
        
        ! In this case the orientation of edge ij remains unchanged
        f_ij = MIN(d_ij,l_ji)*(diff+tstep*dsign*(hstep_ik-hstep_jk))
        fluxloc(iperturb,iloc) = f_ij
        
        IF (i .EQ. k) THEN

          ! For node k which is the upwind node
          pploc(iperturb,0) = pploc(iperturb,0)+MAX(0._DP, f_ij)
          pmloc(iperturb,0) = pmloc(iperturb,0)+MIN(0._DP, f_ij)
          qploc(iperturb,0) = qploc(iperturb,0)+MAX(0._DP,-f_ij)
          qmloc(iperturb,0) = qmloc(iperturb,0)+MIN(0._DP,-f_ij)
          
          ! For node l opposite to k which is the downwind node
          qploc(iperturb,iloc) = qploc(iperturb,iloc)+MAX(0._DP,f_ij)
          qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+MIN(0._DP,f_ij)

        ELSE

          ! For node k which is the downwind node
          qploc(iperturb,0) = qploc(iperturb,0)+MAX(0._DP,f_ij)
          qmloc(iperturb,0) = qmloc(iperturb,0)+MIN(0._DP,f_ij)
          
          ! For node l opposite to k
          pploc(iperturb,iloc) = pploc(iperturb,iloc)+MAX(0._DP, f_ij)
          pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+MIN(0._DP, f_ij)
          qploc(iperturb,iloc) = qploc(iperturb,iloc)+MAX(0._DP,-f_ij)
          qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+MIN(0._DP,-f_ij)
        END IF
      END DO
    END SUBROUTINE do_femtvd_update


    !**************************************************************
    ! Assemble the given column of the Jacobian for FEM-TVD
    SUBROUTINE do_femtvd_assemble(IverticesAtEdge, Kdiagonal, Kcol, flux,&
        Kloc, rploc, rmloc, fluxloc, hstep, iedge, iloc, k, l, bextend, Ksep, Jac)

      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN)  :: IverticesAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      REAL(DP), DIMENSION(:), INTENT(IN)                :: flux
      INTEGER(PREC_VECIDX), DIMENSION(:,:), INTENT(IN)  :: Kloc
      REAL(DP), DIMENSION(:,0:), INTENT(IN)             :: rploc,rmloc
      REAL(DP), DIMENSION(:,0:), INTENT(IN)             :: fluxloc
      REAL(DP), INTENT(IN)                              :: hstep
      INTEGER(PREC_MATIDX), INTENT(IN)                  :: iedge
      INTEGER, INTENT(IN)                               :: iloc
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: k,l
      LOGICAL, INTENT(IN)                               :: bextend

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: Jac

      ! local variables
      INTEGER(PREC_MATIDX) :: ik,jk
      INTEGER(PREC_VECIDX) :: i,j,m
      REAL(DP)             :: f_ij
      INTEGER              :: iperturb
      
      ! Get global node number for edge IJ and the 
      ! number of the node m which is not l
      i = IverticesAtEdge(1,iedge)
      j = IverticesAtEdge(2,iedge)
      m = (i+j)-l
      
      ! We need to find out, which kind of edge is processed
      IF (m .EQ. k) THEN

        !------------------------------------------------------------
        ! 1. Case: primary edge
        !------------------------------------------------------------
        ! The current edge connects the perturbed node k with its
        ! direct neighbor l. Hence, all required information can be
        ! extracted from the local arrays and no global data
        ! retrieval has to be performed.
        !
        ! (a) The edge orientation needs to be adjusted for each
        !     perturbation direction
        
        DO iperturb = 1, 2
          
          ! Retrieve precomputed flux
          f_ij = fluxloc(iperturb,iloc)

          ! Adjust edge orientation
          i = Kloc(2*iperturb,iloc)
          j = Kloc(2*iperturb+1,iloc)
          
          ! Which node is located upwind?
          IF (i .EQ. k) THEN
            
            ! Get corresponding matrix indices
            ik = Kdiagonal(i); jk = Ksep(j)
            
            ! Limit flux 
            IF (f_ij > 0.0_DP) THEN
              f_ij = rploc(iperturb,0)*f_ij
            ELSE
              f_ij = rmloc(iperturb,0)*f_ij
            END IF
            
          ELSE
            
            ! Get corresponding matrix indices
            jk = Kdiagonal(j); ik = Ksep(i)
            
            ! Limit flux
            IF (f_ij > 0.0_DP) THEN
              f_ij = rploc(iperturb,iloc)*f_ij
            ELSE
              f_ij = rmloc(iperturb,iloc)*f_ij
            END IF
            
          END IF
          
          ! Adopt sign for perturbation direction
          f_ij = -(iperturb-1.5_DP)*f_ij/hstep

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij
        END DO
        
      ELSEIF (bextend) THEN
        
        !------------------------------------------------------------
        ! 2. Case: secondary edge
        !------------------------------------------------------------
        ! The current edge connects two nodes l and m which both are
        ! not equal to the perturbed vertex k. Thus, the influence of
        ! the solution perturbation can only be due to a change in
        ! the correction factors alpha_ij. Moreover, for upwind
        ! -biased flux limiting techniques only the nodal correction
        ! factors for the upwind node i is used. Hence, it suffices
        ! to check if node i corresponds to the modified vertex l.
        ! Interestingly enough, some edge LM which connects two
        ! direct neighbors of the perturbed vertex k is only
        ! processed once due to the fact that either l or (!) m
        ! corresponds to the upwind node.

        IF (i .EQ. l) THEN

          IF (flux(iedge) > 0.0_DP) THEN
            f_ij = 0.5_DP*(rploc(1,iloc)-rploc(2,iloc))*flux(iedge)/hstep
          ELSE
            f_ij = 0.5_DP*(rmloc(1,iloc)-rmloc(2,iloc))*flux(iedge)/hstep
          END IF
          
          ! Get corresponding matrix indices
          ik = Ksep(i); jk = Ksep(j)

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij

        END IF

      END IF
    END SUBROUTINE do_femtvd_assemble

    
    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-GP
    SUBROUTINE do_femgp(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
        IsubdiagonalEdgesIdx, IsubdiagonalEdges, DcoefficientsAtEdge,&
        Kld, Kcol, MC, u, u0, flux, flux0, pp, pm, qp, qm, rp, rm,&
        theta, tstep, hstep, NEQ, NEDGE, NNVEDGE, bextend, Ksep, Jac)
    
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN)  :: IverticesAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: IsubdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: IsubdiagonalEdges
      REAL(DP), DIMENSION(:,:), INTENT(IN)              :: DcoefficientsAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      REAL(DP), DIMENSION(:), INTENT(IN)                :: MC
      REAL(DP), DIMENSION(:), INTENT(IN)                :: u,u0
      REAL(DP), DIMENSION(:), INTENT(IN)                :: flux,flux0
      REAL(DP), DIMENSION(:), INTENT(IN)                :: pp,pm,qp,qm,rp,rm
      REAL(DP), INTENT(IN)                              :: theta,tstep,hstep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER(PREC_MATIDX), INTENT(IN)                  :: NEDGE
      INTEGER, INTENT(IN)                               :: NNVEDGE
      LOGICAL, INTENT(IN)                               :: bextend
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: Jac
      
      ! local variables
      INTEGER(PREC_VECIDX), DIMENSION(5,NNVEDGE) :: Kloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: pploc,pmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: qploc,qmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: rploc,rmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: fluxloc,fluxloc0

      INTEGER(PREC_MATIDX) :: ild,iedge
      INTEGER(PREC_VECIDX) :: k,l
      INTEGER :: iloc,nloc

      ! Loop over all columns of the Jacobian matrix
      DO k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        pploc(:,0) = 0; pmloc(:,0) = 0
        qploc(:,0) = 0; qmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        DO ild = IsubdiagonalEdgesIdx(k), IsubdiagonalEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagonalEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Update local coefficients
          CALL do_femgp_update(IverticesAtEdge, DcoefficientsAtEdge,&
              MC, u, u0, flux, flux0, pp, pm, qp, qm,&
              theta, tstep, hstep, iedge, iloc, k,&
              pploc, pmloc, qploc, qmloc, fluxloc, fluxloc0, Kloc)
        END DO

        ! Loop over all superdiagonal edges
        DO iedge = IsuperdiagonalEdgesIdx(k), IsuperdiagonalEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                    
          ! Update local coefficients
          CALL do_femgp_update(IverticesAtEdge, DcoefficientsAtEdge,&
              MC, u, u0, flux, flux0, pp, pm, qp, qm,&
              theta, tstep, hstep, iedge, iloc, k,&
              pploc, pmloc, qploc, qmloc, fluxloc, fluxloc0, Kloc)
        END DO
        
        ! Save total number of local neighbors
        nloc = iloc

        
        ! Adjust the diagonal separator
        CALL do_adjustKsep(Kld, Kcol, k, Ksep)

        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        rploc(:,0:nloc) = afcstab_limit( pploc(:,0:nloc), qploc(:,0:nloc), 0._DP, 1._DP)
        rmloc(:,0:nloc) = afcstab_limit(-pmloc(:,0:nloc),-qmloc(:,0:nloc), 0._DP, 1._DP)


        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        DO iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          DO ild = IsubdiagonalEdgesIdx(l), IsubdiagonalEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagonalEdges(ild)
            
            CALL do_femgp_assemble(IverticesAtEdge, Kld, Kcol, flux, flux0,&
                rp, rm, Kloc, rploc, rmloc, fluxloc, fluxloc0,&
                hstep, iedge, iloc, k, l, bextend, Ksep, Jac)
          END DO

          ! Loop over all superdiagonal edges
          DO iedge = IsuperdiagonalEdgesIdx(l), IsuperdiagonalEdgesIdx(l+1)-1
            
            CALL do_femgp_assemble(IverticesAtEdge, Kld,Kcol, flux, flux0,&
                rp, rm, Kloc, rploc, rmloc, fluxloc, fluxloc0,&
                hstep, iedge, iloc, k, l, bextend, Ksep, Jac)
          END DO
        END DO
      END DO   ! end-of k-loop
    END SUBROUTINE do_femgp

    
    !**************************************************************
    ! Update the local coefficients for FEM-GP
    SUBROUTINE do_femgp_update(IverticesAtEdge, DcoefficientsAtEdge,&
        MC, u, u0, flux, flux0, pp, pm, qp, qm, theta, tstep, hstep,&
        iedge, iloc, k, pploc, pmloc, qploc, qmloc, fluxloc, fluxloc0, Kloc)
      
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN)     :: IverticesAtEdge
      REAL(DP), DIMENSION(:,:), INTENT(IN)                 :: DcoefficientsAtEdge
      REAL(DP), DIMENSION(:), INTENT(IN)                   :: MC
      REAL(DP), DIMENSION(:), INTENT(IN)                   :: u,u0
      REAL(DP), DIMENSION(:), INTENT(IN)                   :: flux,flux0
      REAL(DP), DIMENSION(:), INTENT(IN)                   :: pp,pm,qp,qm
      REAL(DP), INTENT(IN)                                 :: theta,tstep,hstep
      INTEGER(PREC_MATIDX), INTENT(IN)                     :: iedge
      INTEGER, INTENT(IN)                                  :: iloc
      INTEGER(PREC_VECIDX), INTENT(IN)                     :: k

      ! We actually know, that all local quantities start at index zero
      REAL(DP), DIMENSION(:,0:), INTENT(INOUT)             :: pploc,pmloc
      REAL(DP), DIMENSION(:,0:), INTENT(INOUT)             :: qploc,qmloc
      REAL(DP), DIMENSION(:,0:), INTENT(INOUT)             :: fluxloc,fluxloc0
      INTEGER(PREC_VECIDX), DIMENSION(:,:), INTENT(INOUT)  :: Kloc

      ! local variables
      INTEGER(PREC_VECIDX)        :: i,j,ij
      REAL(DP) :: m_ij,d_ij,df_ij,f_ij,l_ij,l_ji,p_ij,pf_ij,q_ij,q_ji
      REAL(DP) :: diff,diff1,diff0,hstep_ik,hstep_jk,dsign
      INTEGER  :: iperturb

      
      ! Determine indices. Obviously, either i or j must be equal
      ! to k. Otherwise, the edge ij would not be present in the
      ! list of incident edges for node k.
      i  = IverticesAtEdge(1,iedge)
      j  = IverticesAtEdge(2,iedge)
      ij = IverticesAtEdge(3,iedge)
      
      !------------------------------------------------------------
      ! (1) unperturbed values: Retrieve the global Ps and Qs and
      !     copy their content to the local ones. Moreover,
      !     "eliminate" the contribution of the edge IJ for the
      !     unperturbed solution values u_i and u_j.
      !------------------------------------------------------------
      ! Determine coefficients
      d_ij = DcoefficientsAtEdge(1,iedge)
      l_ij = DcoefficientsAtEdge(2,iedge)
      l_ji = DcoefficientsAtEdge(3,iedge)
      
      ! Include consistent mass matrix
      m_ij = MC(ij)
      q_ij = m_ij/tstep+l_ij
      q_ji = m_ij/tstep+l_ji

      ! Determine solution differences
      diff1 = u(i)-u(j)
      diff0 = u0(i)-u0(j)

      ! Determine total solution difference
      diff = tstep*(theta*diff1+(1-theta)*diff0)

      ! Compute antidiffusive flux 
      IF(ABS(diff) < SYS_EPSREAL) THEN
        p_ij = 0
        f_ij = 0
      ELSE
        p_ij = MAX(0._DP, m_ij*(diff1-diff0)/diff+d_ij)
        f_ij = p_ij*diff
      END IF

      ! Prelimit the antidiffusive flux
      pf_ij = MIN(p_ij, l_ji)*diff
      
      ! Compute the remaining flux
      df_ij = f_ij-pf_ij

      IF (i .EQ. k) THEN

        ! Store global node number of the opposite node
        Kloc(1,iloc) = j

        ! Compute signed perturbation parameters
        hstep_ik = hstep; hstep_jk = 0._DP
        
        ! Update nodal coefficients for vertex j (!) which is the downwind node
        pploc(:,iloc) = pp(j)-MAX(0._DP,-df_ij)
        pmloc(:,iloc) = pm(j)-MIN(0._DP,-df_ij)
        qploc(:,iloc) = qp(j)-MAX(0._DP, diff)*q_ji
        qmloc(:,iloc) = qm(j)-MIN(0._DP, diff)*q_ji

      ELSE

        ! Store global node number of the opposite node
        Kloc(1,iloc) = i

        ! Compute signed perturbation parameters
        hstep_ik = 0._DP; hstep_jk = hstep
        
        ! Update nodal coefficients for vertex i (!) which is the upwind node
        pploc(:,iloc) = pp(i)-MAX(0._DP, f_ij)
        pmloc(:,iloc) = pm(i)-MIN(0._DP, f_ij)
        qploc(:,iloc) = qp(i)-MAX(0._DP,-diff)*q_ij
        qmloc(:,iloc) = qm(i)-MIN(0._DP,-diff)*q_ij

      END IF

      !------------------------------------------------------------
      ! (2) perturbed values: Now, the local Ps and Qs still
      !     require the contribution of the perturbed solution
      !     values u +/- h*e_k, whereby e_k denotes the k-th unit
      !     vector and h stands for the perturbation step length
      !------------------------------------------------------------
                  
      !------------------------------------------------------------
      ! (3) perform the perturbation for "+/-h*e_k"
      !------------------------------------------------------------
        
      DO iperturb = 1, 2
        
        ! Compute correct sign of perturbation
        dsign = -2*iperturb+3
        
        ! The velocity is assumed to be linear.
        ! Hence the orientation convention for the edge ij is preserved

        ! Save oriented node numbers
        Kloc(2*iperturb:2*iperturb+1,iloc) = (/i,j/)
        
        ! Update solution difference
        diff1 = u(i)-u(j)+dsign*(hstep_ik-hstep_jk)
        
        ! Update total solution difference
        diff = tstep*(theta*diff1+(1-theta)*diff0)
        
        ! Compute antidiffusive flux
        IF (ABS(diff) < SYS_EPSREAL) THEN
          p_ij = 0
          f_ij = 0
        ELSE
          p_ij = MAX(0._DP, m_ij*(diff1-diff0)/diff+d_ij)
          f_ij = p_ij*diff
        END IF
        
        ! Prelimit the antidiffusive flux
        pf_ij = MIN(p_ij, l_ji)*diff
        fluxloc0(iperturb,iloc) = pf_ij
        
        ! Compute the remaining flux
        df_ij = f_ij-pf_ij
        fluxloc(iperturb,iloc) = df_ij
        
        IF (i .EQ. k) THEN

          ! For node k which is the upwind node
          pploc(iperturb,0) = pploc(iperturb,0)+MAX(0._DP, f_ij)
          pmloc(iperturb,0) = pmloc(iperturb,0)+MIN(0._DP, f_ij)
          qploc(iperturb,0) = qploc(iperturb,0)+MAX(0._DP,-diff)*q_ij
          qmloc(iperturb,0) = qmloc(iperturb,0)+MIN(0._DP,-diff)*q_ij
          
          ! For node l opposite to k which is the downwind node
          pploc(iperturb,iloc) = pploc(iperturb,iloc)+MAX(0._DP,-df_ij)
          pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+MIN(0._DP,-df_ij)
          qploc(iperturb,iloc) = qploc(iperturb,iloc)+MAX(0._DP, diff)*q_ji
          qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+MIN(0._DP, diff)*q_ji

        ELSE

          ! For node k which is the downwind node
          pploc(iperturb,0) = pploc(iperturb,0)+MAX(0._DP,-df_ij)
          pmloc(iperturb,0) = pmloc(iperturb,0)+MIN(0._DP,-df_ij)
          qploc(iperturb,0) = qploc(iperturb,0)+MAX(0._DP, diff)*q_ji
          qmloc(iperturb,0) = qmloc(iperturb,0)+MIN(0._DP, diff)*q_ji
          
          ! For node l opposite to k
          pploc(iperturb,iloc) = pploc(iperturb,iloc)+MAX(0._DP, f_ij)
          pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+MIN(0._DP, f_ij)
          qploc(iperturb,iloc) = qploc(iperturb,iloc)+MAX(0._DP,-diff)*q_ij
          qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+MIN(0._DP,-diff)*q_ij

        END IF        
      END DO
    END SUBROUTINE do_femgp_update


    !**************************************************************
    ! Assemble the given column of the Jacobian for FEM-GP
    SUBROUTINE do_femgp_assemble(IverticesAtEdge, Kdiagonal, Kcol,&
        flux, flux0, rp, rm, Kloc, rploc, rmloc, fluxloc, fluxloc0,&
        hstep, iedge, iloc, k, l, bextend, Ksep, Jac)

      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN)  :: IverticesAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      REAL(DP), DIMENSION(:), INTENT(IN)                :: flux,flux0
      REAL(DP), DIMENSION(:), INTENT(IN)                :: rp,rm
      INTEGER(PREC_VECIDX), DIMENSION(:,:), INTENT(IN)  :: Kloc
      REAL(DP), DIMENSION(:,0:), INTENT(IN)             :: rploc,rmloc
      REAL(DP), DIMENSION(:,0:), INTENT(IN)             :: fluxloc,fluxloc0
      REAL(DP), INTENT(IN)                              :: hstep
      INTEGER(PREC_MATIDX), INTENT(IN)                  :: iedge
      INTEGER, INTENT(IN)                               :: iloc
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: k,l
      LOGICAL, INTENT(IN)                               :: bextend

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: Jac

      ! local variables
      INTEGER(PREC_MATIDX) :: ik,jk
      INTEGER(PREC_VECIDX) :: i,j,m
      REAL(DP)             :: f_ij,pf_ij,df_ij
      INTEGER              :: iperturb
      
      ! Get global node number for edge IJ and the 
      ! number of the node m which is not l
      i=IverticesAtEdge(1,iedge)
      j=IverticesAtEdge(2,iedge)
      m=(i+j)-l
      
      ! We need to find out, which kind of edge is processed
      IF (m .EQ. k) THEN

        !------------------------------------------------------------
        ! 1. Case: primary edge
        !------------------------------------------------------------
        ! The current edge connects the perturbed node k with its
        ! direct neighbor l. Hence, all required information can be
        ! extracted from the local arrays and no global data
        ! retrieval has to be performed.
        !
        ! (a) The edge orientation needs to be adjusted for each
        !     perturbation direction
        
        DO iperturb = 1, 2
          
          ! Retrieve precomputed fluxes
          df_ij = fluxloc(iperturb,iloc)
          pf_ij = fluxloc0(iperturb,iloc)

          ! Adjust edge orientation
          i = Kloc(2*iperturb,iloc)
          j = Kloc(2*iperturb+1,iloc)
          
          ! Which node is located upwind?
          IF (i .EQ. k) THEN
            
            ! Get corresponding matrix indices
            ik = Kdiagonal(i); jk = Ksep(j)
            
            ! Limit upwind contribution
            IF (pf_ij > 0.0_DP) THEN
              pf_ij = rploc(iperturb,0)*pf_ij
            ELSE
              pf_ij = rmloc(iperturb,0)*pf_ij
            END IF

            ! Limit symmetric contribution
            IF (df_ij > 0.0_DP) THEN
              df_ij = MIN(rploc(iperturb,0), rmloc(iperturb,iloc))*df_ij
            ELSE
              df_ij = MIN(rmloc(iperturb,0), rploc(iperturb,iloc))*df_ij
            END IF
            
          ELSE
            
            ! Get corresponding matrix indices
            jk = Kdiagonal(j); ik = Ksep(i)
            
            ! Limit upwind contribution
            IF (pf_ij > 0.0_DP) THEN
              pf_ij = rploc(iperturb,iloc)*pf_ij
            ELSE
              pf_ij = rmloc(iperturb,iloc)*pf_ij
            END IF

            ! Limit symmetric contribution
            IF (df_ij > 0.0_DP) THEN
              df_ij = MIN(rmloc(iperturb,0), rploc(iperturb,iloc))*df_ij
            ELSE
              df_ij = MIN(rploc(iperturb,0), rmloc(iperturb,iloc))*df_ij
            END IF
            
          END IF
          
          ! Combine both contributions and 
          ! adopt sign for perturbation direction
          f_ij = -(iperturb-1.5_DP)*(pf_ij+df_ij)/hstep

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij
        END DO
        
      ELSEIF (bextend) THEN
        
        !------------------------------------------------------------
        ! 2. Case: secondary edge
        !------------------------------------------------------------
        ! The current edge connects two nodes l and m which both are
        ! not equal to the perturbed vertex k. Thus, the influence of
        ! the solution perturbation can only be due to a change in
        ! the correction factors alpha_ij. Moreover, for upwind
        ! -biased flux limiting techniques only the nodal correction
        ! factors for the upwind node i is used. Hence, it suffices
        ! to check if node i corresponds to the modified vertex l.
        ! For the symmetric part, we must also check if node j
        ! correspond to the modified vertex l.
        
        IF (i .EQ. l) THEN

          ! Get precomputed fluxes
          pf_ij = flux0(iedge)
          df_ij = flux(iedge)

          ! Limit upwind contribution
          IF (pf_ij > 0.0_DP) THEN
            pf_ij = (rploc(1,iloc)-rploc(2,iloc))*pf_ij
          ELSE
            pf_ij = (rmloc(1,iloc)-rmloc(2,iloc))*pf_ij
          END IF

          ! Limit symmetric contribution
          IF (df_ij > 0.0_DP) THEN
            df_ij = (MIN(rploc(1,iloc), rm(j))-&
                     MIN(rploc(2,iloc), rm(j)))*df_ij
          ELSE
            df_ij = (MIN(rmloc(1,iloc), rp(j))-&
                     MIN(rmloc(2,iloc), rp(j)))*df_ij
          END IF

          ! Combine both contributions
          f_ij = 0.5_DP*(pf_ij+df_ij)/hstep

          ! Get corresponding matrix indices
          ik=Ksep(i); jk=Ksep(j)

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij

        ELSE

          ! Get precomputed flux (only symmetric part)
          df_ij = flux(iedge)

          ! Limit symmetric contribution
          IF (df_ij > 0.0_DP) THEN
            df_ij = (MIN(rp(i), rmloc(1,iloc))-&
                     MIN(rp(i), rmloc(2,iloc)))*df_ij
          ELSE
            df_ij = (MIN(rm(i), rploc(1,iloc))-&
                     MIN(rm(i), rploc(2,iloc)))*df_ij
          END IF

          ! Compute divided difference
          f_ij = 0.5_DP*df_ij/hstep

          ! Get corresponding matrix indices
          ik = Ksep(i); jk = Ksep(j)
          
          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij
          
        END IF
      END IF
    END SUBROUTINE do_femgp_assemble
  END SUBROUTINE gfsc_buildStabJacLinearScalar_GPTVD

  !*****************************************************************************
  
!<subroutine>

  SUBROUTINE gfsc_buildStabJacobianBlock_FCT(RmatrixC, rmatrixMC, ru, &
      fcb_getVelocity, theta, tstep, hstep, bclear, rafcstab, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation part
    ! of the discrete transport operator for a scalar convection equation.
    ! The velocity is assumed to be nonlinear/arbitrary.
    ! Note that this routine serves as a wrapper for block vectors. If there
    ! is only one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: RmatrixC

    ! consistent mass matrix
    TYPE(t_matrixScalar), INTENT(IN)               :: rmatrixMC

    ! solution vector
    TYPE(t_vectorBlock), INTENT(IN)                :: ru

    ! implicitness parameter
    REAL(DP), INTENT(IN)                           :: theta

    ! time step size
    REAL(DP), INTENT(IN)                           :: tstep

    ! perturbation parameter
    REAL(DP), INTENT(IN)                           :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    LOGICAL, INTENT(IN)                            :: bclear

     ! callback functions to compute velocity
    INCLUDE 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)                 :: rafcstab

    ! Jacobian matrix
    TYPE(t_matrixScalar), INTENT(INOUT)            :: rmatrixJ   
!</inputoutput>
!</subroutine>

    IF (ru%nblocks  .NE. 1) THEN

      CALL output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacobianBlock_FCT')
      CALL sys_halt()

    ELSE

      CALL gfsc_buildStabJacobianScalar_FCT(RmatrixC, rmatrixMC, ru%RvectorBlock(1),&
          fcb_getVelocity, theta, tstep, hstep, bclear, rafcstab, rmatrixJ)

    END IF
  END SUBROUTINE gfsc_buildStabJacobianBlock_FCT

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_buildStabJacobianScalar_FCT(RmatrixC, rmatrixMC, ru,&
      fcb_getVelocity, theta, tstep, hstep, bclear, rafcstab, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation
    ! part of the discrete transport operator for a scalar convection equation.
    ! The velocity is assumed to be nonlinear/arbitrary. 
    ! This routine will also work for linear velocities but then it is inefficient
    ! since the solution perturbation does not affect the velocity.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: RmatrixC

    ! consistent mass matrix
    TYPE(t_matrixScalar), INTENT(IN)               :: rmatrixMC

    ! solution vector
    TYPE(t_vectorScalar), INTENT(IN)               :: ru

    ! implicitness parameter
    REAL(DP), INTENT(IN)                           :: theta

    ! time step size
    REAL(DP), INTENT(IN)                           :: tstep

    ! perturbation parameter
    REAL(DP), INTENT(IN)                           :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    LOGICAL, INTENT(IN)                            :: bclear

     ! callback functions to compute velocity
    INCLUDE 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)                 :: rafcstab

    ! Jacobian matrix
    TYPE(t_matrixScalar), INTENT(INOUT)            :: rmatrixJ   
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_MATIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_Kld
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_Kdiagonal
    REAL(DP), DIMENSION(:,:), POINTER             :: p_DcoefficientsAtEdge
    REAL(DP), DIMENSION(:), POINTER               :: p_flux,p_flux0
    REAL(DP), DIMENSION(:), POINTER               :: p_u
    REAL(DP), DIMENSION(:), POINTER               :: p_Cx,p_Cy,p_Cz,p_MC,p_Jac
    INTEGER :: ndim

    
    ! Check if stabilisation is prepared
    IF (rafcstab%ctypeAFCstabilisation .NE. AFCSTAB_FEMFCT .OR.&
        IAND(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE) .EQ. 0 .OR.&
        IAND(rafcstab%iSpec, AFCSTAB_EDGEVALUES)    .EQ. 0 .OR.&
        IAND(rafcstab%iSpec, AFCSTAB_FLUXES)        .EQ. 0) THEN
      CALL output_line('Stabilisation does not provide required structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacobianScalar_FCT')
      CALL sys_halt()
    END IF
    
    ! Clear matrix?
    IF (bclear) CALL lsyssc_clearMatrix(rmatrixJ)

    ! Set pointers
    CALL afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    CALL afcstab_getbase_DcoeffsAtEdge(rafcstab,   p_DcoefficientsAtEdge)
    CALL lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_flux)
    CALL lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_flux0)
    CALL lsyssc_getbase_double(rmatrixMC, p_MC)
    CALL lsyssc_getbase_double(rmatrixJ,  p_Jac)
    CALL lsyssc_getbase_double(ru,        p_u)
    
    ! Set spatial dimensions
    ndim = SIZE(RmatrixC,1)

    ! What kind of matrix format are we?
    SELECT CASE(rmatrixJ%cmatrixFormat)
    CASE(LSYSSC_MATRIX7)
      CALL lsyssc_getbase_Kld(rmatrixJ, p_Kld)
      
      ! How many dimensions do we have?
      SELECT CASE(ndim)
      CASE (NDIM1D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        
        CALL do_femfct_1D(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
            p_Kld, p_Cx, p_MC, p_u, p_flux, p_flux0, theta, tstep, hstep,&
            rafcstab%NEDGE, (rafcstab%imass .EQ. AFCSTAB_CONSISTENTMASS), p_Jac)
        
      CASE (NDIM2D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
        
        CALL do_femfct_2D(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
            p_Kld, p_Cx, p_Cy, p_MC, p_u, p_flux, p_flux0, theta, tstep, hstep,&
            rafcstab%NEDGE, (rafcstab%imass .EQ. AFCSTAB_CONSISTENTMASS), p_Jac)
        
      CASE (NDIM3D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)        
        CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
        CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
        
        CALL do_femfct_3D(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
            p_Kld, p_Cx, p_Cy, p_Cz, p_MC, p_u, p_flux, p_flux0, theta, tstep, hstep,&
            rafcstab%NEDGE, (rafcstab%imass .EQ. AFCSTAB_CONSISTENTMASS), p_Jac)
      END SELECT
      
      
    CASE(LSYSSC_MATRIX9)
      CALL lsyssc_getbase_Kdiagonal(rmatrixJ, p_Kdiagonal)
      
      ! How many dimensions do we have?
      SELECT CASE(ndim)
      CASE (NDIM1D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        
        CALL do_femfct_1D(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
            p_Kdiagonal, p_Cx, p_MC, p_u, p_flux, p_flux0, theta, tstep, hstep,&
            rafcstab%NEDGE, (rafcstab%imass .EQ. AFCSTAB_CONSISTENTMASS), p_Jac)
        
      CASE (NDIM2D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)
        CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
        
        CALL do_femfct_2D(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
            p_Kdiagonal, p_Cx, p_Cy, p_MC, p_u, p_flux, p_flux0, theta, tstep, hstep,&
            rafcstab%NEDGE, (rafcstab%imass .EQ. AFCSTAB_CONSISTENTMASS), p_Jac)
        
      CASE (NDIM3D)
        CALL lsyssc_getbase_double(RmatrixC(1), p_Cx)        
        CALL lsyssc_getbase_double(RmatrixC(2), p_Cy)
        CALL lsyssc_getbase_double(RmatrixC(3), p_Cz)
        
        CALL do_femfct_3D(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
            p_Kdiagonal, p_Cx, p_Cy, p_Cz, p_MC, p_u, p_flux, p_flux0, theta, tstep, hstep,&
            rafcstab%NEDGE, (rafcstab%imass .EQ. AFCSTAB_CONSISTENTMASS), p_Jac)
      END SELECT
      
    CASE DEFAULT
      CALL output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacobianScalar_FCT')
      CALL sys_halt()
    END SELECT

  CONTAINS
    
    ! Here, the working routine follow

     !**************************************************************
    ! Assemble the Jacobian matrix for FEM-FCT in 1D
    ! All matrices can be stored in matrix format 7 or 9
    SUBROUTINE do_femfct_1D(IverticesAtEdge, DcoefficientsAtEdge, Kdiagonal,&
        Cx, MC, u, flux, flux0, theta, tstep, hstep, NEDGE, bmass, Jac)

      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtEdge
      REAL(DP), DIMENSION(:,:), INTENT(IN)             :: DcoefficientsAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)   :: Kdiagonal
      REAL(DP), DIMENSION(:), INTENT(IN)               :: Cx,MC
      REAL(DP), DIMENSION(:), INTENT(IN)               :: flux,flux0
      REAL(DP), DIMENSION(:), INTENT(IN)               :: u
      REAL(DP), INTENT(IN)                             :: theta,tstep,hstep
      INTEGER(PREC_MATIDX), INTENT(IN)                 :: NEDGE
      LOGICAL, INTENT(IN)                              :: bmass
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: Jac

      ! local variables
      REAL(DP), DIMENSION(NDIM1D) :: v_ij,v_ji
      INTEGER(PREC_MATIDX) :: iedge,ij,ji,ii,jj
      INTEGER(PREC_VECIDX) :: i,j
      REAL(DP) :: f_i,f_j,f_ij,d_ij,a_ij,b_ij,l_ij,l_ji
      REAL(DP) :: diff,diff_i,diff_j

      ! Should we apply the consistent mass matrix?
      IF (bmass) THEN

        ! Loop over all edges
        DO iedge = 1, NEDGE
          
          ! Determine vertex numbers
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)
          
          ! Determine diagonal indices
          ii = Kdiagonal(i); jj = Kdiagonal(j)
          
          ! Compute solution difference
          diff = u(i)-u(j)
          
          ! Determine perturbed solution differences
          diff_i = diff+hstep
          diff_j = diff-hstep
          
          !------------------------------------------------------------
          ! Compute flux for i-th column
          !------------------------------------------------------------
          
          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i)+hstep, u(j), i, j, v_ij, v_ji)

          ! Compute perturbed coefficient a_ij(u+hstep*e_i)
          l_ij = -v_ij(1)*Cx(ij)
          l_ji = -v_ji(1)*Cx(ji)
          a_ij = MC(ij)+theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_i)
          f_i = a_ij*diff_i+flux0(iedge)
          IF (f_i > 0.0_DP) THEN
            f_i = MIN(f_i, MAX(flux(iedge), 0._DP))
          ELSE
            f_i = MAX(f_i, MIN(flux(iedge), 0._DP))
          END IF

          
          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i)-hstep, u(j), i, j, v_ij, v_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_i)
          l_ij = -v_ij(1)*Cx(ij)
          l_ji = -v_ji(1)*Cx(ji)
          b_ij = MC(ij)+theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_j+flux0(iedge)
          IF (f_j > 0.0_DP) THEN
            f_j = MIN(f_j, MAX(flux(iedge), 0._DP))
          ELSE
            f_j = MAX(f_j, MIN(flux(iedge), 0._DP))
          END IF

          
          ! Compute divided differences of fluxes
          f_ij = 0.5_DP*(f_i-f_j)/hstep
          
          ! Apply i-th column
          Jac(ii) = Jac(ii)-f_ij
          Jac(ji) = Jac(ji)+f_ij

          !------------------------------------------------------------
          ! Compute flux for j-th column
          !------------------------------------------------------------
          
          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i), u(j)+hstep, i, j, v_ij, v_ji)

          ! Compute perturbed coefficient a_ij(u+hstep*e_j)
          l_ij = -v_ij(1)*Cx(ij)
          l_ji = -v_ji(1)*Cx(ji)
          a_ij = MC(ij)+theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_j) 
          f_i = a_ij*diff_j+flux0(iedge)
          IF (f_i > 0.0_DP) THEN
            f_i = MIN(f_i, MAX(flux(iedge), 0._DP))
          ELSE
            f_i = MAX(f_i, MIN(flux(iedge), 0._DP))
          END IF


          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i), u(j)-hstep, i, j, v_ij, v_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_j)
          l_ij = -v_ij(1)*Cx(ij)
          l_ji = -v_ji(1)*Cx(ji)
          b_ij = MC(ij)+theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_i+flux0(iedge)
          IF (f_j > 0.0_DP) THEN
            f_j = MIN(f_j, MAX(flux(iedge), 0._DP))
          ELSE
            f_j = MAX(f_j, MIN(flux(iedge), 0._DP))
          END IF
          

          ! Compute divided differences of fluxes
          f_ij = 0.5_DP*(f_i-f_j)/hstep
          
          ! Apply j-th column
          Jac(ij) = Jac(ij)-f_ij
          Jac(jj) = Jac(jj)+f_ij
        END DO

      ELSE

        ! Loop over all edges
        DO iedge = 1, NEDGE
          
          ! Determine vertex numbers
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)
          
          ! Determine diagonal indices
          ii = Kdiagonal(i); jj = Kdiagonal(j)
                    
          ! Compute solution difference
          diff = u(i)-u(j)
          
          ! Determine perturbed solution differences
          diff_i = diff+hstep
          diff_j = diff-hstep
          
          
          !------------------------------------------------------------
          ! Compute flux for i-th column
          !------------------------------------------------------------
          
          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i)+hstep, u(j), i, j, v_ij, v_ji)

          ! Compute perturbed coefficient a_ij(u+hstep*e_i)
          l_ij = -v_ij(1)*Cx(ij)
          l_ji = -v_ji(1)*Cx(ji)
          a_ij = theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_i)
          f_i = a_ij*diff_i+flux0(iedge)
          IF (f_i > 0.0_DP) THEN
            f_i = MIN(f_i, MAX(flux(iedge), 0._DP))
          ELSE
            f_i = MAX(f_i, MIN(flux(iedge), 0._DP))
          END IF


          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i)-hstep, u(j), i, j, v_ij, v_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_i)
          l_ij = -v_ij(1)*Cx(ij)
          l_ji = -v_ji(1)*Cx(ji)
          b_ij = theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
        
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_j+flux0(iedge)
          IF (f_j > 0.0_DP) THEN
            f_j = MIN(f_j, MAX(flux(iedge), 0._DP))
          ELSE
            f_j = MAX(f_j, MIN(flux(iedge), 0._DP))
          END IF
          

          ! Compute divided differences of fluxes
          f_ij = 0.5_DP*(f_i-f_j)/hstep
          
          ! Apply i-th column
          Jac(ii) = Jac(ii)-f_ij
          Jac(ji) = Jac(ji)+f_ij

          !------------------------------------------------------------
          ! Compute flux for j-th column
          !------------------------------------------------------------
          
          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i), u(j)+hstep, i, j, v_ij, v_ji)

          ! Compute perturbed coefficient a_ij(u+hstep*e_j)
          l_ij = -v_ij(1)*Cx(ij)
          l_ji = -v_ji(1)*Cx(ji)
          a_ij = theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_j)
          f_i = a_ij*diff_j+flux0(iedge)
          IF (f_i > 0.0_DP) THEN
            f_i = MIN(f_i, MAX(flux(iedge), 0._DP))
          ELSE
            f_i = MAX(f_i, MIN(flux(iedge), 0._DP))
          END IF


          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i), u(j)-hstep, i, j, v_ij, v_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_j)
          l_ij = -v_ij(1)*Cx(ij)
          l_ji = -v_ji(1)*Cx(ji)
          b_ij = theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_i+flux0(iedge)
          IF (f_j > 0.0_DP) THEN
            f_j = MIN(f_j, MAX(flux(iedge), 0._DP))
          ELSE
            f_j = MAX(f_j, MIN(flux(iedge), 0._DP))
          END IF
          

          ! Compute divided differences of fluxes
          f_ij = 0.5_DP*(f_i-f_j)/hstep
          
          ! Apply j-th column
          Jac(ij) = Jac(ij)-f_ij
          Jac(jj) = Jac(jj)+f_ij
        END DO
        
      END IF
    END SUBROUTINE do_femfct_1D


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-FCT in 2D
    ! All matrices can be stored in matrix format 7 or 9
    SUBROUTINE do_femfct_2D(IverticesAtEdge, DcoefficientsAtEdge, Kdiagonal,&
        Cx, Cy, MC, u, flux, flux0, theta, tstep, hstep, NEDGE, bmass, Jac)

      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtEdge
      REAL(DP), DIMENSION(:,:), INTENT(IN)             :: DcoefficientsAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)   :: Kdiagonal
      REAL(DP), DIMENSION(:), INTENT(IN)               :: Cx,Cy,MC
      REAL(DP), DIMENSION(:), INTENT(IN)               :: flux,flux0
      REAL(DP), DIMENSION(:), INTENT(IN)               :: u
      REAL(DP), INTENT(IN)                             :: theta,tstep,hstep
      INTEGER(PREC_MATIDX), INTENT(IN)                 :: NEDGE
      LOGICAL, INTENT(IN)                              :: bmass
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: Jac

      ! local variables
      REAL(DP), DIMENSION(NDIM2D) :: v_ij,v_ji
      INTEGER(PREC_MATIDX) :: iedge,ij,ji,ii,jj
      INTEGER(PREC_VECIDX) :: i,j
      REAL(DP) :: f_i,f_j,f_ij,d_ij,a_ij,b_ij,l_ij,l_ji
      REAL(DP) :: diff,diff_i,diff_j

      ! Should we apply the consistent mass matrix?
      IF (bmass) THEN

        ! Loop over all edges
        DO iedge = 1, NEDGE
          
          ! Determine vertex numbers
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)
          
          ! Determine diagonal indices
          ii = Kdiagonal(i); jj = Kdiagonal(j)
          
          ! Compute solution difference
          diff = u(i)-u(j)
          
          ! Determine perturbed solution differences
          diff_i = diff+hstep
          diff_j = diff-hstep
          
          !------------------------------------------------------------
          ! Compute flux for i-th column
          !------------------------------------------------------------
          
          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i)+hstep, u(j), i, j, v_ij, v_ji)

          ! Compute perturbed coefficient a_ij(u+hstep*e_i)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          a_ij = MC(ij)+theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_i)
          f_i = a_ij*diff_i+flux0(iedge)
          IF (f_i > 0.0_DP) THEN
            f_i = MIN(f_i, MAX(flux(iedge), 0._DP))
          ELSE
            f_i = MAX(f_i, MIN(flux(iedge), 0._DP))
          END IF

          
          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i)-hstep, u(j), i, j, v_ij, v_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_i)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          b_ij = MC(ij)+theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_j+flux0(iedge)
          IF (f_j > 0.0_DP) THEN
            f_j = MIN(f_j, MAX(flux(iedge), 0._DP))
          ELSE
            f_j = MAX(f_j, MIN(flux(iedge), 0._DP))
          END IF

          
          ! Compute divided differences of fluxes
          f_ij = 0.5_DP*(f_i-f_j)/hstep
          
          ! Apply i-th column
          Jac(ii) = Jac(ii)-f_ij
          Jac(ji) = Jac(ji)+f_ij

          !------------------------------------------------------------
          ! Compute flux for j-th column
          !------------------------------------------------------------
          
          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i), u(j)+hstep, i, j, v_ij, v_ji)

          ! Compute perturbed coefficient a_ij(u+hstep*e_j)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          a_ij = MC(ij)+theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_j) 
          f_i = a_ij*diff_j+flux0(iedge)
          IF (f_i > 0.0_DP) THEN
            f_i = MIN(f_i, MAX(flux(iedge), 0._DP))
          ELSE
            f_i = MAX(f_i, MIN(flux(iedge), 0._DP))
          END IF


          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i), u(j)-hstep, i, j, v_ij, v_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_j)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          b_ij = MC(ij)+theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_i+flux0(iedge)
          IF (f_j > 0.0_DP) THEN
            f_j = MIN(f_j, MAX(flux(iedge), 0._DP))
          ELSE
            f_j = MAX(f_j, MIN(flux(iedge), 0._DP))
          END IF
          

          ! Compute divided differences of fluxes
          f_ij = 0.5_DP*(f_i-f_j)/hstep
          
          ! Apply j-th column
          Jac(ij) = Jac(ij)-f_ij
          Jac(jj) = Jac(jj)+f_ij
        END DO

      ELSE

        ! Loop over all edges
        DO iedge = 1, NEDGE
          
          ! Determine vertex numbers
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)
          
          ! Determine diagonal indices
          ii = Kdiagonal(i); jj = Kdiagonal(j)
                    
          ! Compute solution difference
          diff = u(i)-u(j)
          
          ! Determine perturbed solution differences
          diff_i = diff+hstep
          diff_j = diff-hstep
          
          
          !------------------------------------------------------------
          ! Compute flux for i-th column
          !------------------------------------------------------------
          
          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i)+hstep, u(j), i, j, v_ij, v_ji)

          ! Compute perturbed coefficient a_ij(u+hstep*e_i)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          a_ij = theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_i)
          f_i = a_ij*diff_i+flux0(iedge)
          IF (f_i > 0.0_DP) THEN
            f_i = MIN(f_i, MAX(flux(iedge), 0._DP))
          ELSE
            f_i = MAX(f_i, MIN(flux(iedge), 0._DP))
          END IF


          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i)-hstep, u(j), i, j, v_ij, v_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_i)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          b_ij = theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
        
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_j+flux0(iedge)
          IF (f_j > 0.0_DP) THEN
            f_j = MIN(f_j, MAX(flux(iedge), 0._DP))
          ELSE
            f_j = MAX(f_j, MIN(flux(iedge), 0._DP))
          END IF
          

          ! Compute divided differences of fluxes
          f_ij = 0.5_DP*(f_i-f_j)/hstep
          
          ! Apply i-th column
          Jac(ii) = Jac(ii)-f_ij
          Jac(ji) = Jac(ji)+f_ij

          !------------------------------------------------------------
          ! Compute flux for j-th column
          !------------------------------------------------------------
          
          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i), u(j)+hstep, i, j, v_ij, v_ji)

          ! Compute perturbed coefficient a_ij(u+hstep*e_j)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          a_ij = theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_j)
          f_i = a_ij*diff_j+flux0(iedge)
          IF (f_i > 0.0_DP) THEN
            f_i = MIN(f_i, MAX(flux(iedge), 0._DP))
          ELSE
            f_i = MAX(f_i, MIN(flux(iedge), 0._DP))
          END IF


          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i), u(j)-hstep, i, j, v_ij, v_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_j)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)
          b_ij = theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_i+flux0(iedge)
          IF (f_j > 0.0_DP) THEN
            f_j = MIN(f_j, MAX(flux(iedge), 0._DP))
          ELSE
            f_j = MAX(f_j, MIN(flux(iedge), 0._DP))
          END IF
          

          ! Compute divided differences of fluxes
          f_ij = 0.5_DP*(f_i-f_j)/hstep
          
          ! Apply j-th column
          Jac(ij) = Jac(ij)-f_ij
          Jac(jj) = Jac(jj)+f_ij
        END DO
        
      END IF
    END SUBROUTINE do_femfct_2D


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-FCT in 3D
    ! All matrices can be stored in matrix format 7 or 9
    SUBROUTINE do_femfct_3D(IverticesAtEdge, DcoefficientsAtEdge, Kdiagonal,&
        Cx, Cy, Cz, MC, u, flux, flux0, theta, tstep, hstep, NEDGE, bmass, Jac)

      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtEdge
      REAL(DP), DIMENSION(:,:), INTENT(IN)             :: DcoefficientsAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)   :: Kdiagonal
      REAL(DP), DIMENSION(:), INTENT(IN)               :: Cx,Cy,Cz,MC
      REAL(DP), DIMENSION(:), INTENT(IN)               :: flux,flux0
      REAL(DP), DIMENSION(:), INTENT(IN)               :: u
      REAL(DP), INTENT(IN)                             :: theta,tstep,hstep
      INTEGER(PREC_MATIDX), INTENT(IN)                 :: NEDGE
      LOGICAL, INTENT(IN)                              :: bmass
      REAL(DP), DIMENSION(:), INTENT(INOUT)            :: Jac

      ! local variables
      REAL(DP), DIMENSION(NDIM3D) :: v_ij,v_ji
      INTEGER(PREC_MATIDX) :: iedge,ij,ji,ii,jj
      INTEGER(PREC_VECIDX) :: i,j
      REAL(DP) :: f_i,f_j,f_ij,d_ij,a_ij,b_ij,l_ij,l_ji
      REAL(DP) :: diff,diff_i,diff_j

      ! Should we apply the consistent mass matrix?
      IF (bmass) THEN

        ! Loop over all edges
        DO iedge = 1, NEDGE
          
          ! Determine vertex numbers
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)
          
          ! Determine diagonal indices
          ii = Kdiagonal(i); jj = Kdiagonal(j)
          
          ! Compute solution difference
          diff = u(i)-u(j)
          
          ! Determine perturbed solution differences
          diff_i = diff+hstep
          diff_j = diff-hstep
          
          !------------------------------------------------------------
          ! Compute flux for i-th column
          !------------------------------------------------------------
          
          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i)+hstep, u(j), i, j, v_ij, v_ji)

          ! Compute perturbed coefficient a_ij(u+hstep*e_i)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          a_ij = MC(ij)+theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_i)
          f_i = a_ij*diff_i+flux0(iedge)
          IF (f_i > 0.0_DP) THEN
            f_i = MIN(f_i, MAX(flux(iedge), 0._DP))
          ELSE
            f_i = MAX(f_i, MIN(flux(iedge), 0._DP))
          END IF

          
          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i)-hstep, u(j), i, j, v_ij, v_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_i)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          b_ij = MC(ij)+theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_j+flux0(iedge)
          IF (f_j > 0.0_DP) THEN
            f_j = MIN(f_j, MAX(flux(iedge), 0._DP))
          ELSE
            f_j = MAX(f_j, MIN(flux(iedge), 0._DP))
          END IF

          
          ! Compute divided differences of fluxes
          f_ij = 0.5_DP*(f_i-f_j)/hstep
          
          ! Apply i-th column
          Jac(ii) = Jac(ii)-f_ij
          Jac(ji) = Jac(ji)+f_ij

          !------------------------------------------------------------
          ! Compute flux for j-th column
          !------------------------------------------------------------
          
          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i), u(j)+hstep, i, j, v_ij, v_ji)

          ! Compute perturbed coefficient a_ij(u+hstep*e_j)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          a_ij = MC(ij)+theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_j) 
          f_i = a_ij*diff_j+flux0(iedge)
          IF (f_i > 0.0_DP) THEN
            f_i = MIN(f_i, MAX(flux(iedge), 0._DP))
          ELSE
            f_i = MAX(f_i, MIN(flux(iedge), 0._DP))
          END IF


          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i), u(j)-hstep, i, j, v_ij, v_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_j)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          b_ij = MC(ij)+theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_i+flux0(iedge)
          IF (f_j > 0.0_DP) THEN
            f_j = MIN(f_j, MAX(flux(iedge), 0._DP))
          ELSE
            f_j = MAX(f_j, MIN(flux(iedge), 0._DP))
          END IF
          

          ! Compute divided differences of fluxes
          f_ij = 0.5_DP*(f_i-f_j)/hstep
          
          ! Apply j-th column
          Jac(ij) = Jac(ij)-f_ij
          Jac(jj) = Jac(jj)+f_ij
        END DO

      ELSE

        ! Loop over all edges
        DO iedge = 1, NEDGE
          
          ! Determine vertex numbers
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)
          
          ! Determine diagonal indices
          ii = Kdiagonal(i); jj = Kdiagonal(j)
                    
          ! Compute solution difference
          diff = u(i)-u(j)
          
          ! Determine perturbed solution differences
          diff_i = diff+hstep
          diff_j = diff-hstep
          
          
          !------------------------------------------------------------
          ! Compute flux for i-th column
          !------------------------------------------------------------
          
          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i)+hstep, u(j), i, j, v_ij, v_ji)

          ! Compute perturbed coefficient a_ij(u+hstep*e_i)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          a_ij = theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_i)
          f_i = a_ij*diff_i+flux0(iedge)
          IF (f_i > 0.0_DP) THEN
            f_i = MIN(f_i, MAX(flux(iedge), 0._DP))
          ELSE
            f_i = MAX(f_i, MIN(flux(iedge), 0._DP))
          END IF


          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i)-hstep, u(j), i, j, v_ij, v_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_i)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          b_ij = theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
        
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_j+flux0(iedge)
          IF (f_j > 0.0_DP) THEN
            f_j = MIN(f_j, MAX(flux(iedge), 0._DP))
          ELSE
            f_j = MAX(f_j, MIN(flux(iedge), 0._DP))
          END IF
          

          ! Compute divided differences of fluxes
          f_ij = 0.5_DP*(f_i-f_j)/hstep
          
          ! Apply i-th column
          Jac(ii) = Jac(ii)-f_ij
          Jac(ji) = Jac(ji)+f_ij

          !------------------------------------------------------------
          ! Compute flux for j-th column
          !------------------------------------------------------------
          
          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i), u(j)+hstep, i, j, v_ij, v_ji)

          ! Compute perturbed coefficient a_ij(u+hstep*e_j)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          a_ij = theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_j)
          f_i = a_ij*diff_j+flux0(iedge)
          IF (f_i > 0.0_DP) THEN
            f_i = MIN(f_i, MAX(flux(iedge), 0._DP))
          ELSE
            f_i = MAX(f_i, MIN(flux(iedge), 0._DP))
          END IF


          ! Compute perturbed velocity
          CALL fcb_getVelocity(u(i), u(j)-hstep, i, j, v_ij, v_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_j)
          l_ij = -v_ij(1)*Cx(ij)-v_ij(2)*Cy(ij)-v_ij(3)*Cz(ij)
          l_ji = -v_ji(1)*Cx(ji)-v_ji(2)*Cy(ji)-v_ji(3)*Cz(ji)
          b_ij = theta*tstep*MAX(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_i+flux0(iedge)
          IF (f_j > 0.0_DP) THEN
            f_j = MIN(f_j, MAX(flux(iedge), 0._DP))
          ELSE
            f_j = MAX(f_j, MIN(flux(iedge), 0._DP))
          END IF
          

          ! Compute divided differences of fluxes
          f_ij = 0.5_DP*(f_i-f_j)/hstep
          
          ! Apply j-th column
          Jac(ij) = Jac(ij)-f_ij
          Jac(jj) = Jac(jj)+f_ij
        END DO
        
      END IF
    END SUBROUTINE do_femfct_3D
  END SUBROUTINE gfsc_buildStabJacobianScalar_FCT

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_buildStabJacobianBlock_GPTVD(RmatrixC, rmatrixMC, ru, ru0,&
      fcb_getVelocity, theta, tstep, hstep, bclear, rafcstab, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation part
    ! of the discrete transport operator for a scalar convection equation.
    ! The velocity is assumed to be nonlinear/arbitrary.
    ! Note that this routine serves as a wrapper for block vectors. If there
    ! is only one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: RmatrixC

    ! consistent mass matrix
    TYPE(t_matrixScalar), INTENT(IN)               :: rmatrixMC

    ! solution vector
    TYPE(t_vectorBlock), INTENT(IN)                :: ru

    ! initial solution vector
    TYPE(t_vectorBlock), INTENT(IN)                :: ru0

    ! implicitness parameter
    REAL(DP), INTENT(IN)                           :: theta

    ! time step size
    REAL(DP), INTENT(IN)                           :: tstep

    ! perturbation parameter
    REAL(DP), INTENT(IN)                           :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    LOGICAL, INTENT(IN)                            :: bclear

     ! callback functions to compute velocity
    INCLUDE 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)                 :: rafcstab

    ! Jacobian matrix
    TYPE(t_matrixScalar), INTENT(INOUT)            :: rmatrixJ   
!</inputoutput>
!</subroutine>

    IF (ru%nblocks  .NE. 1 .OR.&
        ru0%nblocks .NE. 1) THEN

      CALL output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacobianBlock_GPTVD')
      CALL sys_halt()

    ELSE

      CALL gfsc_buildStabJacobianScalar_GPTVD(RmatrixC, rmatrixMC, ru%RvectorBlock(1),&
          ru0%RvectorBlock(1), fcb_getVelocity, theta, tstep, hstep, bclear, rafcstab, rmatrixJ)

    END IF
  END SUBROUTINE gfsc_buildStabJacobianBlock_GPTVD
  
  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_buildStabJacobianScalar_GPTVD(RmatrixC, rmatrixMC, ru, ru0,&
      fcb_getVelocity, theta, tstep, hstep, bclear, rafcstab, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation
    ! part of the discrete transport operator for a scalar convection equation.
    ! The velocity is assumed to be nonlinear/arbitrary. 
    ! This routine will also work for linear velocities but then it is inefficient
    ! since the solution perturbation does not affect the velocity.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    TYPE(t_matrixScalar), DIMENSION(:), INTENT(IN) :: RmatrixC

    ! consistent mass matrix
    TYPE(t_matrixScalar), INTENT(IN)               :: rmatrixMC

    ! solution vector
    TYPE(t_vectorScalar), INTENT(IN)               :: ru

    ! initial solution vector
    TYPE(t_vectorScalar), INTENT(IN)               :: ru0

    ! implicitness parameter
    REAL(DP), INTENT(IN)                           :: theta

    ! time step size
    REAL(DP), INTENT(IN)                           :: tstep

    ! perturbation parameter
    REAL(DP), INTENT(IN)                           :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    LOGICAL, INTENT(IN)                            :: bclear

     ! callback functions to compute velocity
    INCLUDE 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)                 :: rafcstab

    ! Jacobian matrix
    TYPE(t_matrixScalar), INTENT(INOUT)            :: rmatrixJ   
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_MATIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER   :: p_IsuperdiagonalEdgesIdx
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_IsubdiagonalEdges
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_IsubdiagonalEdgesIdx
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_Kld,p_Ksep
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_Kdiagonal
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER   :: p_Kcol
    REAL(DP), DIMENSION(:,:), POINTER             :: p_DcoefficientsAtEdge
    REAL(DP), DIMENSION(:), POINTER               :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    REAL(DP), DIMENSION(:), POINTER               :: p_flux,p_flux0
    REAL(DP), DIMENSION(:), POINTER               :: p_u,p_u0
    REAL(DP), DIMENSION(:), POINTER               :: p_Cx,p_Cy,p_Cz,p_MC,p_Jac
    INTEGER :: h_Ksep
    INTEGER :: ndim
    LOGICAL :: bextend
    
    ! Set spatial dimensions
    ndim = SIZE(RmatrixC,1)
    
    ! Clear matrix?
    IF (bclear) CALL lsyssc_clearMatrix(rmatrixJ)
    
    ! What kind of stabilisation are we?
    SELECT CASE(rafcstab%ctypeAFCstabilisation)
      
    CASE(AFCSTAB_FEMTVD)
      
      ! Check if stabilisation is prepared
      IF (IAND(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)   .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_EDGEORIENTATION) .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_EDGEVALUES)      .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_ANTIDIFFUSION)   .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_BOUNDS)          .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_FLUXES)          .EQ. 0) THEN
        CALL output_line('Stabilisation does not provide required structures',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacobianScalar_GPTVD')
        CALL sys_halt()
      END IF

      ! Check if subdiagonal edges need to be generated
      IF (IAND(rafcstab%iSpec, AFCSTAB_SUBDIAGONALEDGES) .EQ. 0)&
          CALL afcstab_generateSubdiagEdges(rafcstab)

      ! Set pointers
      CALL afcstab_getbase_IverticesAtEdge(rafcstab,  p_IverticesAtEdge)
      CALL afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
      CALL afcstab_getbase_IsubdiagEdge(rafcstab,    p_IsubdiagonalEdges)
      CALL afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagonalEdgesIdx)
      CALL afcstab_getbase_DcoeffsAtEdge(rafcstab,    p_DcoefficientsAtEdge)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
      CALL lsyssc_getbase_double(rafcstab%RedgeVectors(1),  p_flux)
      CALL lsyssc_getbase_Kcol(rmatrixJ,   p_Kcol)
      CALL lsyssc_getbase_double(rmatrixJ, p_Jac)
      CALL lsyssc_getbase_double(ru,       p_u)

      ! What kind of matrix format are we?
      SELECT CASE(rmatrixJ%cmatrixFormat)
      CASE(LSYSSC_MATRIX7)
        
        ! Set pointers
        CALL lsyssc_getbase_Kld(rmatrixJ, p_Kld)

        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        CALL storage_copy(rmatrixJ%h_Kld, h_Ksep)
        CALL storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)
        CALL lalg_vectorAddScalarInt(p_Ksep, 1)

        ! Assembled extended Jacobian matrix
        bextend = (rafcstab%iextendedJacobian .NE. 0)

        ! How many dimensions do we have?
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)

          CALL do_femtvd_1D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
            p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
            p_Kld, p_Kcol, p_Cx, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
            theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)

        CASE (NDIM2D)
          CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(rmatrixC(2), p_Cy)

          CALL do_femtvd_2D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
            p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
            p_Kld, p_Kcol, p_Cx, p_Cy, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
            theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)
          
        CASE (NDIM3D)
          CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(rmatrixC(2), p_Cy)
          CALL lsyssc_getbase_double(rmatrixC(3), p_Cz)

          CALL do_femtvd_3D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
            p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
            p_Kld, p_Kcol, p_Cx, p_Cy, p_Cz, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
            theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)
        END SELECT
               
        ! Free storage
        CALL storage_free(h_Ksep)

        
      CASE(LSYSSC_MATRIX9)

        ! Set pointers
        CALL lsyssc_getbase_Kdiagonal(rmatrixJ, p_Kdiagonal)
        
        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        CALL storage_copy(rmatrixJ%h_Kld, h_Ksep)
        CALL storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)

        ! Assembled extended Jacobian matrix
        bextend = (rafcstab%iextendedJacobian .NE. 0)

        ! How many dimensions do we have?
        SELECT CASE(ndim)
        CASE (NDIM1D)
          CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)
          
          CALL do_femtvd_1D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
              p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
              p_Kdiagonal, p_Kcol, p_Cx, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
              theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
              rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)

        CASE (NDIM2D)
          CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(rmatrixC(2), p_Cy)

          CALL do_femtvd_2D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
              p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
              p_Kdiagonal, p_Kcol, p_Cx, p_Cy, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
              theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
              rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)

        CASE (NDIM3D)
          CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)
          CALL lsyssc_getbase_double(rmatrixC(2), p_Cy)
          CALL lsyssc_getbase_double(rmatrixC(3), p_Cz)

          CALL do_femtvd_3D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
              p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
              p_Kdiagonal, p_Kcol, p_Cx, p_Cy, p_Cz, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
              theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
              rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)
        END SELECT
        
        ! Free storage
        CALL storage_free(h_Ksep)

        
      CASE DEFAULT
        CALL output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacobianScalar_GPTVD')
        CALL sys_halt()
      END SELECT

      
    CASE(AFCSTAB_FEMGP)
      
      ! Check if stabilisation is prepared
      IF (IAND(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)   .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_EDGEORIENTATION) .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_EDGEVALUES)      .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_ANTIDIFFUSION)   .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_BOUNDS)          .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_FLUXES)          .EQ. 0) THEN
        CALL output_line('Stabilisation does not provide required structures',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacobianScalar_GPTVD')
        CALL sys_halt()
      END IF

      ! Check if subdiagonal edges need to be generated
      IF (IAND(rafcstab%iSpec, AFCSTAB_SUBDIAGONALEDGES) .EQ. 0)&
          CALL afcstab_generateSubdiagEdges(rafcstab)

      ! Set pointers
      CALL afcstab_getbase_IverticesAtEdge(rafcstab,  p_IverticesAtEdge)
      CALL afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
      CALL afcstab_getbase_IsubdiagEdge(rafcstab,    p_IsubdiagonalEdges)
      CALL afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagonalEdgesIdx)
      CALL afcstab_getbase_DcoeffsAtEdge(rafcstab,    p_DcoefficientsAtEdge)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
      CALL lsyssc_getbase_double(rafcstab%RedgeVectors(1),  p_flux)
      CALL lsyssc_getbase_double(rafcstab%RedgeVectors(2),  p_flux0)
      CALL lsyssc_getbase_Kcol(rmatrixJ,   p_Kcol)
      CALL lsyssc_getbase_double(rmatrixJ, p_Jac)
      CALL lsyssc_getbase_double(ru,       p_u)
      CALL lsyssc_getbase_double(ru0,      p_u0)

      ! What kind of matrix format are we?
      SELECT CASE(rmatrixJ%cmatrixFormat)
      CASE(LSYSSC_MATRIX7)
        
        ! Set pointers
        CALL lsyssc_getbase_Kld(rmatrixJ, p_Kld)

        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        CALL storage_copy(rmatrixJ%h_Kld, h_Ksep)
        CALL storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)
        CALL lalg_vectorAddScalarInt(p_Ksep, 1)

        ! Assembled extended Jacobian matrix
        bextend = (rafcstab%iextendedJacobian .NE. 0)
        IF (rafcstab%imass .EQ. AFCSTAB_CONSISTENTMASS) THEN
          
          ! Set pointers
          CALL lsyssc_getbase_double(rmatrixMC, p_MC)
          CALL lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
          CALL lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)

          ! How many dimensions do we have?
          SELECT CASE(ndim)
          CASE (NDIM1D)
            CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)

            CALL do_femgp_1D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
                p_Kld, p_Kcol, p_Cx, p_MC, p_u, p_u0, p_flux, p_flux0,&
                p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, theta, tstep, hstep,&
                rafcstab%NEQ, rafcstab%NEDGE, rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)

          CASE (NDIM2D)
            CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)
            CALL lsyssc_getbase_double(rmatrixC(2), p_Cy)

            CALL do_femgp_2D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
                p_Kld, p_Kcol, p_Cx, p_Cy, p_MC, p_u, p_u0, p_flux, p_flux0,&
                p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, theta, tstep, hstep,&
                rafcstab%NEQ, rafcstab%NEDGE, rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)
            
          CASE (NDIM3D)
            CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)
            CALL lsyssc_getbase_double(rmatrixC(2), p_Cy)
            CALL lsyssc_getbase_double(rmatrixC(3), p_Cz)

            CALL do_femgp_3D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
                p_Kld, p_Kcol, p_Cx, p_Cy, p_Cz, p_MC, p_u, p_u0, p_flux, p_flux0,&
                p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, theta, tstep, hstep,&
                rafcstab%NEQ, rafcstab%NEDGE, rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)
          END SELECT

        ELSE

          ! How many dimensions do we have?
          SELECT CASE(ndim)
          CASE (NDIM1D)
            CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)
            
            CALL do_femtvd_1D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
              p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
              p_Kld, p_Kcol, p_Cx, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
              theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
              rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)
            
          CASE (NDIM2D)
            CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)
            CALL lsyssc_getbase_double(rmatrixC(2), p_Cy)

            CALL do_femtvd_2D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
                p_Kld, p_Kcol, p_Cx, p_Cy, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)
            
          CASE (NDIM3D)
            CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)
            CALL lsyssc_getbase_double(rmatrixC(2), p_Cy)
            CALL lsyssc_getbase_double(rmatrixC(3), p_Cz)

            CALL do_femtvd_3D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
                p_Kld, p_Kcol, p_Cx, p_Cy, p_Cz, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)
          END SELECT

        END IF

        ! Free storage
        CALL storage_free(h_Ksep)
        
      CASE(LSYSSC_MATRIX9)

        ! Set pointers
        CALL lsyssc_getbase_Kdiagonal(rmatrixJ, p_Kdiagonal)
        
        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        CALL storage_copy(rmatrixJ%h_Kld, h_Ksep)
        CALL storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)

        ! Assembled extended Jacobian matrix
        bextend = (rafcstab%iextendedJacobian .NE. 0)
        IF (rafcstab%imass .EQ. AFCSTAB_CONSISTENTMASS) THEN
          
          ! Set pointers
          CALL lsyssc_getbase_double(rmatrixMC, p_MC)
          CALL lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
          CALL lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)

          ! How many dimensions do we have?
          SELECT CASE(ndim)
          CASE (NDIM1D)
            CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)
            
            CALL do_femgp_1D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
                p_Kdiagonal, p_Kcol, p_Cx, p_MC, p_u, p_u0, p_flux, p_flux0,&
                p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, theta, tstep, hstep,&
                rafcstab%NEQ, rafcstab%NEDGE, rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)

          CASE (NDIM2D)
            CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)
            CALL lsyssc_getbase_double(rmatrixC(2), p_Cy)

            CALL do_femgp_2D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
                p_Kdiagonal, p_Kcol, p_Cx, p_Cy, p_MC, p_u, p_u0, p_flux, p_flux0,&
                p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, theta, tstep, hstep,&
                rafcstab%NEQ, rafcstab%NEDGE, rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)
            
          CASE (NDIM3D)
            CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)
            CALL lsyssc_getbase_double(rmatrixC(2), p_Cy)
            CALL lsyssc_getbase_double(rmatrixC(3), p_Cz)

            CALL do_femgp_3D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
                p_Kdiagonal, p_Kcol, p_Cx, p_Cy, p_Cz, p_MC, p_u, p_u0, p_flux, p_flux0,&
                p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, theta, tstep, hstep,&
                rafcstab%NEQ, rafcstab%NEDGE, rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)
          END SELECT

        ELSE

          ! How many dimensions do we have?
          SELECT CASE(ndim)
          CASE (NDIM1D)
            CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)
            
            CALL do_femtvd_1D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
                p_Kdiagonal, p_Kcol, p_Cx, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)
            
          CASE (NDIM2D)
            CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)
            CALL lsyssc_getbase_double(rmatrixC(2), p_Cy)
            
            CALL do_femtvd_2D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
                p_Kdiagonal, p_Kcol, p_Cx, p_Cy, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)
            
          CASE (NDIM3D)
            CALL lsyssc_getbase_double(rmatrixC(1), p_Cx)
            CALL lsyssc_getbase_double(rmatrixC(2), p_Cy)
            CALL lsyssc_getbase_double(rmatrixC(3), p_Cz)

            CALL do_femtvd_3D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
                p_Kdiagonal, p_Kcol, p_Cx, p_Cy, p_Cz, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)
          END SELECT

        END IF

        ! Free storage
        CALL storage_free(h_Ksep)
        
      CASE DEFAULT
        CALL output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacobianScalar_GPTVD')
        CALL sys_halt()
      END SELECT

    CASE DEFAULT
      CALL output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacobianScalar_GPTVD')
      CALL sys_halt()
    END SELECT

  CONTAINS

    ! Here, the working routine follow
    
    !**************************************************************    
    ! Adjust the diagonal separator.
    ! The separator is initialied by the column separator (increased
    ! by one if the is necessary for matrix format 7).
    ! Based on the matric structure given by Kld/Kcol, the separator
    ! is "moved" to the given column "k". For efficiency reasons, only
    ! those entries are considered which are present in column "k".
    SUBROUTINE do_adjustKsep(Kld, Kcol, k, Ksep)
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: k
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      
      INTEGER(PREC_MATIDX) :: ild,isep
      INTEGER(PREC_VECIDX) :: l
      INTEGER :: iloc

      ! Loop over all entries of the k-th row
      DO ild = Kld(k), Kld(k+1)-1
        
        ! Get the column number and the position of the separator
        l = Kcol(ild); isep = Ksep(l)

        ! If the separator does not point to the k-th column
        ! it must be adjusted accordingly
        IF (Kcol(isep) < k) Ksep(l) = Ksep(l)+1
      END DO
    END SUBROUTINE do_adjustKsep


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-TVD in 1D
    SUBROUTINE do_femtvd_1D(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
        IsubdiagonalEdgesIdx, IsubdiagonalEdges, DcoefficientsAtEdge,&
        Kld, Kcol, Cx, u, flux, pp, pm, qp, qm, theta, tstep, hstep,&
        NEQ, NEDGE, NNVEDGE, bextend, Ksep, Jac)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN)  :: IverticesAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: IsubdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: IsubdiagonalEdges
      REAL(DP), DIMENSION(:,:), INTENT(IN)              :: DcoefficientsAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx
      REAL(DP), DIMENSION(:), INTENT(IN)                :: u
      REAL(DP), DIMENSION(:), INTENT(IN)                :: flux
      REAL(DP), DIMENSION(:), INTENT(IN)                :: pp,pm,qp,qm
      REAL(DP), INTENT(IN)                              :: theta,tstep,hstep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER(PREC_MATIDX), INTENT(IN)                  :: NEDGE
      INTEGER, INTENT(IN)                               :: NNVEDGE
      LOGICAL, INTENT(IN)                               :: bextend
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: Jac
      
      ! local variables
      INTEGER(PREC_VECIDX), DIMENSION(5,NNVEDGE) :: Kloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: pploc,pmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: qploc,qmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: rploc,rmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: fluxloc
      REAL(DP), DIMENSION(NDIM1D)      :: c_ij, c_ji

      INTEGER(PREC_MATIDX) :: ij,ji,ild,iedge
      INTEGER(PREC_VECIDX) :: i,j,k,l
      INTEGER :: iloc,nloc

      ! Loop over all columns of the Jacobian matrix
      DO k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        pploc(:,0) = 0; pmloc(:,0) = 0
        qploc(:,0) = 0; qmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        DO ild = IsubdiagonalEdgesIdx(k), IsubdiagonalEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagonalEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = Cx(ij)
          c_ji = Cx(ji)
          
          ! Update local coefficients
          CALL do_femtvd_update(DcoefficientsAtEdge, u, pp, pm, qp, qm,&
              c_ij, c_ji, tstep, hstep, iedge, i, j, ij, ji, NDIM1D,&
              iloc, k, pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        END DO

        ! Loop over all superdiagonal edges
        DO iedge = IsuperdiagonalEdgesIdx(k), IsuperdiagonalEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                    
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = Cx(ij)
          c_ji = Cx(ji)

          ! Update local coefficients
          CALL do_femtvd_update(DcoefficientsAtEdge, u, pp, pm, qp, qm,&
              c_ij, c_ji, tstep, hstep, iedge, i, j, ij, ji, NDIM1D,&
              iloc, k, pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        END DO
        
        ! Save total number of local neighbors
        nloc = iloc
        
        ! Adjust the diagonal separator
        CALL do_adjustKsep(Kld, Kcol, k, Ksep)

        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        rploc(:,0:nloc) = afcstab_limit( pploc(:,0:nloc), qploc(:,0:nloc), 0._DP, 1._DP)
        rmloc(:,0:nloc) = afcstab_limit(-pmloc(:,0:nloc),-qmloc(:,0:nloc), 0._DP, 1._DP)

        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        DO iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          DO ild = IsubdiagonalEdgesIdx(l), IsubdiagonalEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagonalEdges(ild)
            
            CALL do_femtvd_assemble(IverticesAtEdge, Kld, Kcol,&
                flux, Kloc, rploc, rmloc, fluxloc, hstep,&
                iedge, iloc, k, l, bextend, Ksep, Jac)
          END DO

          ! Loop over all superdiagonal edges
          DO iedge = IsuperdiagonalEdgesIdx(l), IsuperdiagonalEdgesIdx(l+1)-1
            
            CALL do_femtvd_assemble(IverticesAtEdge, Kld, Kcol,&
                flux, Kloc, rploc, rmloc, fluxloc, hstep,&
                iedge, iloc, k, l, bextend, Ksep, Jac)
          END DO
        END DO
      END DO   ! end-of k-loop
    END SUBROUTINE do_femtvd_1D


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-TVD in 2D
    SUBROUTINE do_femtvd_2D(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
        IsubdiagonalEdgesIdx, IsubdiagonalEdges, DcoefficientsAtEdge,&
        Kld, Kcol, Cx, Cy, u, flux, pp, pm, qp, qm, theta, tstep, hstep,&
        NEQ, NEDGE, NNVEDGE, bextend, Ksep, Jac)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN)  :: IverticesAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: IsubdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: IsubdiagonalEdges
      REAL(DP), DIMENSION(:,:), INTENT(IN)              :: DcoefficientsAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy
      REAL(DP), DIMENSION(:), INTENT(IN)                :: u
      REAL(DP), DIMENSION(:), INTENT(IN)                :: flux
      REAL(DP), DIMENSION(:), INTENT(IN)                :: pp,pm,qp,qm
      REAL(DP), INTENT(IN)                              :: theta,tstep,hstep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER(PREC_MATIDX), INTENT(IN)                  :: NEDGE
      INTEGER, INTENT(IN)                               :: NNVEDGE
      LOGICAL, INTENT(IN)                               :: bextend
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: Jac
      
      ! local variables
      INTEGER(PREC_VECIDX), DIMENSION(5,NNVEDGE) :: Kloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: pploc,pmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: qploc,qmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: rploc,rmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: fluxloc
      REAL(DP), DIMENSION(NDIM2D)      :: c_ij, c_ji

      INTEGER(PREC_MATIDX) :: ij,ji,ild,iedge
      INTEGER(PREC_VECIDX) :: i,j,k,l
      INTEGER :: iloc,nloc

      ! Loop over all columns of the Jacobian matrix
      DO k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        pploc(:,0) = 0; pmloc(:,0) = 0
        qploc(:,0) = 0; qmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        DO ild = IsubdiagonalEdgesIdx(k), IsubdiagonalEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagonalEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/Cx(ij),Cy(ij)/)
          c_ji = (/Cx(ji),Cy(ji)/)
          
          ! Update local coefficients
          CALL do_femtvd_update(DcoefficientsAtEdge, u, pp, pm, qp, qm,&
              c_ij, c_ji, tstep, hstep, iedge, i, j, ij, ji, NDIM2D,&
              iloc, k, pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        END DO

        ! Loop over all superdiagonal edges
        DO iedge = IsuperdiagonalEdgesIdx(k), IsuperdiagonalEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                    
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/Cx(ij),Cy(ij)/)
          c_ji = (/Cx(ji),Cy(ji)/)

          ! Update local coefficients
          CALL do_femtvd_update(DcoefficientsAtEdge, u, pp, pm, qp, qm,&
              c_ij, c_ji, tstep, hstep, iedge, i, j, ij, ji, NDIM2D,&
              iloc, k, pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        END DO
        
        ! Save total number of local neighbors
        nloc = iloc
        
        ! Adjust the diagonal separator
        CALL do_adjustKsep(Kld, Kcol, k, Ksep)

        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        rploc(:,0:nloc) = afcstab_limit( pploc(:,0:nloc), qploc(:,0:nloc), 0._DP, 1._DP)
        rmloc(:,0:nloc) = afcstab_limit(-pmloc(:,0:nloc),-qmloc(:,0:nloc), 0._DP, 1._DP)

        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        DO iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          DO ild = IsubdiagonalEdgesIdx(l), IsubdiagonalEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagonalEdges(ild)
            
            CALL do_femtvd_assemble(IverticesAtEdge, Kld, Kcol,&
                flux, Kloc, rploc, rmloc, fluxloc, hstep,&
                iedge, iloc, k, l, bextend, Ksep, Jac)
          END DO

          ! Loop over all superdiagonal edges
          DO iedge = IsuperdiagonalEdgesIdx(l), IsuperdiagonalEdgesIdx(l+1)-1
            
            CALL do_femtvd_assemble(IverticesAtEdge, Kld, Kcol,&
                flux, Kloc, rploc, rmloc, fluxloc, hstep,&
                iedge, iloc, k, l, bextend, Ksep, Jac)
          END DO
        END DO
      END DO   ! end-of k-loop
    END SUBROUTINE do_femtvd_2D


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-TVD in 3D
    SUBROUTINE do_femtvd_3D(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
        IsubdiagonalEdgesIdx, IsubdiagonalEdges, DcoefficientsAtEdge,&
        Kld, Kcol, Cx, Cy, Cz, u, flux, pp, pm, qp, qm, theta, tstep, hstep,&
        NEQ, NEDGE, NNVEDGE, bextend, Ksep, Jac)

      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN)  :: IverticesAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: IsubdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: IsubdiagonalEdges
      REAL(DP), DIMENSION(:,:), INTENT(IN)              :: DcoefficientsAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz
      REAL(DP), DIMENSION(:), INTENT(IN)                :: u
      REAL(DP), DIMENSION(:), INTENT(IN)                :: flux
      REAL(DP), DIMENSION(:), INTENT(IN)                :: pp,pm,qp,qm
      REAL(DP), INTENT(IN)                              :: theta,tstep,hstep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER(PREC_MATIDX), INTENT(IN)                  :: NEDGE
      INTEGER, INTENT(IN)                               :: NNVEDGE
      LOGICAL, INTENT(IN)                               :: bextend
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: Jac
      
      ! local variables
      INTEGER(PREC_VECIDX), DIMENSION(5,NNVEDGE) :: Kloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: pploc,pmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: qploc,qmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: rploc,rmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: fluxloc
      REAL(DP), DIMENSION(NDIM3D)      :: c_ij, c_ji

      INTEGER(PREC_MATIDX) :: ij,ji,ild,iedge
      INTEGER(PREC_VECIDX) :: i,j,k,l
      INTEGER :: iloc,nloc

      ! Loop over all columns of the Jacobian matrix
      DO k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        pploc(:,0) = 0; pmloc(:,0) = 0
        qploc(:,0) = 0; qmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        DO ild = IsubdiagonalEdgesIdx(k), IsubdiagonalEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagonalEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/Cx(ij),Cy(ij),Cz(ij)/)
          c_ji = (/Cx(ji),Cy(ji),Cz(ji)/)
          
          ! Update local coefficients
          CALL do_femtvd_update(DcoefficientsAtEdge, u, pp, pm, qp, qm,&
              c_ij, c_ji, tstep, hstep, iedge, i, j, ij, ji, NDIM3D,&
              iloc, k, pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        END DO

        ! Loop over all superdiagonal edges
        DO iedge = IsuperdiagonalEdgesIdx(k), IsuperdiagonalEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                    
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/Cx(ij),Cy(ij),Cz(ij)/)
          c_ji = (/Cx(ji),Cy(ji),Cz(ji)/)

          ! Update local coefficients
          CALL do_femtvd_update(DcoefficientsAtEdge, u, pp, pm, qp, qm,&
              c_ij, c_ji, tstep, hstep, iedge, i, j, ij, ji, NDIM3D,&
              iloc, k, pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        END DO
        
        ! Save total number of local neighbors
        nloc = iloc
        
        ! Adjust the diagonal separator
        CALL do_adjustKsep(Kld, Kcol, k, Ksep)

        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        rploc(:,0:nloc) = afcstab_limit( pploc(:,0:nloc), qploc(:,0:nloc), 0._DP, 1._DP)
        rmloc(:,0:nloc) = afcstab_limit(-pmloc(:,0:nloc),-qmloc(:,0:nloc), 0._DP, 1._DP)

        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        DO iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          DO ild = IsubdiagonalEdgesIdx(l), IsubdiagonalEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagonalEdges(ild)
            
            CALL do_femtvd_assemble(IverticesAtEdge, Kld, Kcol,&
                flux, Kloc, rploc, rmloc, fluxloc, hstep,&
                iedge, iloc, k, l, bextend, Ksep, Jac)
          END DO

          ! Loop over all superdiagonal edges
          DO iedge = IsuperdiagonalEdgesIdx(l), IsuperdiagonalEdgesIdx(l+1)-1
            
            CALL do_femtvd_assemble(IverticesAtEdge, Kld, Kcol,&
                flux, Kloc, rploc, rmloc, fluxloc, hstep,&
                iedge, iloc, k, l, bextend, Ksep, Jac)
          END DO
        END DO
      END DO   ! end-of k-loop
    END SUBROUTINE do_femtvd_3D


    !**************************************************************
    ! Update the local coefficients for FEM-TVD in arbitrary dimension
    SUBROUTINE do_femtvd_update(DcoefficientsAtEdge, u, pp, pm, qp, qm, &
        c_ij, c_ji, tstep, hstep, iedge, i, j, ij, ji, ndim, iloc, k,&
        pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
      
      REAL(DP), DIMENSION(:,:), INTENT(IN)                 :: DcoefficientsAtEdge
      REAL(DP), DIMENSION(:), INTENT(IN)                   :: u
      REAL(DP), DIMENSION(:), INTENT(IN)                   :: pp,pm,qp,qm
      REAL(DP), DIMENSION(:), INTENT(IN)                   :: c_ij,c_ji      
      REAL(DP), INTENT(IN)                                 :: tstep,hstep
      INTEGER(PREC_MATIDX), INTENT(IN)                     :: iedge
      INTEGER(PREC_VECIDX), INTENT(IN)                     :: i,j
      INTEGER(PREC_MATIDX), INTENT(IN)                     :: ij,ji
      INTEGER, INTENT(IN)                                  :: ndim
      INTEGER, INTENT(IN)                                  :: iloc
      INTEGER(PREC_VECIDX), INTENT(IN)                     :: k

      ! We actually know, that all local quantities start at index zero
      REAL(DP), DIMENSION(:,0:), INTENT(INOUT)             :: pploc,pmloc
      REAL(DP), DIMENSION(:,0:), INTENT(INOUT)             :: qploc,qmloc
      REAL(DP), DIMENSION(:,0:), INTENT(INOUT)             :: fluxloc
      INTEGER(PREC_VECIDX), DIMENSION(:,:), INTENT(INOUT)  :: Kloc

      ! local variables
      REAL(DP), DIMENSION(ndim) :: v_ij,v_ji
      REAL(DP) :: d_ij,f_ij,l_ij,l_ji,diff,hstep_ik,hstep_jk,dsign
      INTEGER  :: iperturb

      !------------------------------------------------------------
      ! (1) unperturbed values: Retrieve the global Ps and Qs and
      !     copy their content to the local ones. Moreover,
      !     "eliminate" the contribution of the edge IJ for the
      !     unperturbed solution values u_i and u_j.
      !------------------------------------------------------------
      ! Determine coefficients
      d_ij = DcoefficientsAtEdge(1,iedge)
      l_ij = DcoefficientsAtEdge(2,iedge)
      l_ji = DcoefficientsAtEdge(3,iedge)
      
      ! Determine prelimited antidiffusive flux
      diff = tstep*(u(i)-u(j))
      f_ij = MIN(d_ij,l_ji)*diff
      
      IF (i .EQ. k) THEN
        
        ! Store global node number of the opposite node
        Kloc(1,iloc) = j

        ! Compute signed perturbation parameters
        hstep_ik = hstep; hstep_jk = 0._DP
        
        ! Update nodal coefficients for vertex j (!) which is the downwind node
        pploc(:,iloc) = pp(j)
        pmloc(:,iloc) = pm(j)
        qploc(:,iloc) = qp(j)-MAX(0._DP, f_ij)
        qmloc(:,iloc) = qm(j)-MIN(0._DP, f_ij)

      ELSE

        ! Store global node number of the opposite node
        Kloc(1,iloc) = i

        ! Compute signed perturbation parameters
        hstep_ik = 0._DP; hstep_jk = hstep
        
        ! Update nodal coefficients for vertex i (!) which is the upwind node
        pploc(:,iloc) = pp(i)-MAX(0._DP, f_ij)
        pmloc(:,iloc) = pm(i)-MIN(0._DP, f_ij)
        qploc(:,iloc) = qp(i)-MAX(0._DP,-f_ij)
        qmloc(:,iloc) = qm(i)-MIN(0._DP,-f_ij)
      END IF

      !------------------------------------------------------------
      ! (2) perturbed values: Now, the local Ps and Qs still
      !     require the contribution of the perturbed solution
      !     values u +/- h*e_k, whereby e_k denotes the k-th unit
      !     vector and h stands for the perturbation step length
      !------------------------------------------------------------

      DO iperturb = 1, 2
        
        ! Compute correct sign of perturbation
        dsign = -2*iperturb+3
        
        ! Compute perturbed velocity
        CALL fcb_getVelocity(u(i)+dsign*hstep_ik, u(j)+dsign*hstep_jk, i, j, v_ij, v_ji)
        
        ! Compute perturbation coefficients k_ij and k_ji
        l_ij = -SUM(v_ij*c_ij)
        l_ji = -SUM(v_ji*c_ji)

        ! Compute diffusion coefficient
        d_ij = MAX(-l_ij, 0._DP, -l_ji)

        ! Apply discrete upwinding
        l_ij = l_ij+d_ij
        l_ji = l_ji+d_ij
        
        ! Due to the (possible) nonlinearity of the velocity vector
        ! the orientation convention for the edge ij may be violated,
        ! that is, the condition 0=l_ij < l_ji may not be valid. In this
        ! case the node number i and j must be swapped logically
        IF (l_ij .LE. l_ji) THEN
          
          ! Save oriented node numbers
          Kloc(2*iperturb:2*iperturb+1,iloc) = (/i,j/)
          
          ! In this case the orientation of edge ij remains unchanged
          f_ij = MIN(d_ij,l_ji)*(diff+tstep*dsign*(hstep_ik-hstep_jk))
          fluxloc(iperturb,iloc) = f_ij
        
          IF (i .EQ. k) THEN

            ! For node k which is the upwind node
            pploc(iperturb,0) = pploc(iperturb,0)+MAX(0._DP, f_ij)
            pmloc(iperturb,0) = pmloc(iperturb,0)+MIN(0._DP, f_ij)
            qploc(iperturb,0) = qploc(iperturb,0)+MAX(0._DP,-f_ij)
            qmloc(iperturb,0) = qmloc(iperturb,0)+MIN(0._DP,-f_ij)
            
            ! For node l opposite to k which is the downwind node
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+MAX(0._DP,f_ij)
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+MIN(0._DP,f_ij)

          ELSE

            ! For node k which is the downwind node
            qploc(iperturb,0) = qploc(iperturb,0)+MAX(0._DP,f_ij)
            qmloc(iperturb,0) = qmloc(iperturb,0)+MIN(0._DP,f_ij)
            
            ! For node l opposite to k
            pploc(iperturb,iloc) = pploc(iperturb,iloc)+MAX(0._DP, f_ij)
            pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+MIN(0._DP, f_ij)
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+MAX(0._DP,-f_ij)
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+MIN(0._DP,-f_ij)

          END IF
          
        ELSE
          
          ! Save oriented node numbers
          Kloc(2*iperturb:2*iperturb+1,iloc) = (/j,i/)
          
          ! In this case the orientation of edge ij needs to be
          ! reverted so as to let i denote the 'upwind' node
          f_ij = -MIN(d_ij,l_ij)*(diff+tstep*dsign*(hstep_ik-hstep_jk))
          fluxloc(iperturb,iloc) = f_ij
          
          IF (j .EQ. k) THEN

            ! For node k which is the upwind node
            pploc(iperturb,0) = pploc(iperturb,0)+MAX(0._DP, f_ij)
            pmloc(iperturb,0) = pmloc(iperturb,0)+MIN(0._DP, f_ij)
            qploc(iperturb,0) = qploc(iperturb,0)+MAX(0._DP,-f_ij)
            qmloc(iperturb,0) = qmloc(iperturb,0)+MIN(0._DP,-f_ij)
            
            ! For node "l" opposite to k which is the downwind node
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+MAX(0._DP,f_ij)
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+MIN(0._DP,f_ij)

          ELSE

            ! For node k which is the downwind node
            qploc(iperturb,0) = qploc(iperturb,0)+MAX(0._DP,f_ij)
            qmloc(iperturb,0) = qmloc(iperturb,0)+MIN(0._DP,f_ij)
            
            ! For node "l" opposite to k which is the upwind node
            pploc(iperturb,iloc) = pploc(iperturb,iloc)+MAX(0._DP, f_ij)
            pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+MIN(0._DP, f_ij)
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+MAX(0._DP,-f_ij)
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+MIN(0._DP,-f_ij)
          END IF
        END IF
      END DO
    END SUBROUTINE do_femtvd_update


    !**************************************************************
    ! Assemble the given column of the Jacobian for FEM-TVD in arbitrary dimension
    SUBROUTINE do_femtvd_assemble(IverticesAtEdge, Kdiagonal, Kcol, flux,&
        Kloc, rploc, rmloc, fluxloc, hstep, iedge, iloc, k, l, bextend, Ksep, Jac)

      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN)  :: IverticesAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      REAL(DP), DIMENSION(:), INTENT(IN)                :: flux
      INTEGER(PREC_VECIDX), DIMENSION(:,:), INTENT(IN)  :: Kloc
      REAL(DP), DIMENSION(:,0:), INTENT(IN)             :: rploc,rmloc
      REAL(DP), DIMENSION(:,0:), INTENT(IN)             :: fluxloc
      REAL(DP), INTENT(IN)                              :: hstep
      INTEGER(PREC_MATIDX), INTENT(IN)                  :: iedge
      INTEGER, INTENT(IN)                               :: iloc
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: k,l
      LOGICAL, INTENT(IN)                               :: bextend

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: Jac

      ! local variables
      INTEGER(PREC_MATIDX) :: ik,jk
      INTEGER(PREC_VECIDX) :: i,j,m
      REAL(DP)             :: f_ij
      INTEGER              :: iperturb
      
      ! Get global node number for edge IJ and the 
      ! number of the node m which is not l
      i = IverticesAtEdge(1,iedge)
      j = IverticesAtEdge(2,iedge)
      m = (i+j)-l
      
      ! We need to find out, which kind of edge is processed
      IF (m .EQ. k) THEN

        !------------------------------------------------------------
        ! 1. Case: primary edge
        !------------------------------------------------------------
        ! The current edge connects the perturbed node k with its
        ! direct neighbor l. Hence, all required information can be
        ! extracted from the local arrays and no global data
        ! retrieval has to be performed.
        !
        ! (a) The edge orientation needs to be adjusted for each
        !     perturbation direction
        
        DO iperturb = 1, 2
          
          ! Retrieve precomputed flux
          f_ij = fluxloc(iperturb,iloc)

          ! Adjust edge orientation
          i = Kloc(2*iperturb,iloc)
          j = Kloc(2*iperturb+1,iloc)
          
          ! Which node is located upwind?
          IF (i .EQ. k) THEN
            
            ! Get corresponding matrix indices
            ik = Kdiagonal(i); jk = Ksep(j)
            
            ! Limit flux 
            IF (f_ij > 0.0_DP) THEN
              f_ij = rploc(iperturb,0)*f_ij
            ELSE
              f_ij = rmloc(iperturb,0)*f_ij
            END IF
            
          ELSE
            
            ! Get corresponding matrix indices
            jk = Kdiagonal(j); ik = Ksep(i)
            
            ! Limit flux
            IF (f_ij > 0.0_DP) THEN
              f_ij = rploc(iperturb,iloc)*f_ij
            ELSE
              f_ij = rmloc(iperturb,iloc)*f_ij
            END IF
            
          END IF
          
          ! Adopt sign for perturbation direction
          f_ij = -(iperturb-1.5_DP)*f_ij/hstep

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij
        END DO
        
      ELSEIF (bextend) THEN
        
        !------------------------------------------------------------
        ! 2. Case: secondary edge
        !------------------------------------------------------------
        ! The current edge connects two nodes l and m which both are
        ! not equal to the perturbed vertex k. Thus, the influence of
        ! the solution perturbation can only be due to a change in
        ! the correction factors alpha_ij. Moreover, for upwind
        ! -biased flux limiting techniques only the nodal correction
        ! factors for the upwind node i is used. Hence, it suffices
        ! to check if node i corresponds to the modified vertex l.
        ! Interestingly enough, some edge LM which connects two
        ! direct neighbors of the perturbed vertex k is only
        ! processed once due to the fact that either l or (!) m
        ! corresponds to the upwind node.

        IF (i .EQ. l) THEN

          IF (flux(iedge) > 0.0_DP) THEN
            f_ij = 0.5_DP*(rploc(1,iloc)-rploc(2,iloc))*flux(iedge)/hstep
          ELSE
            f_ij = 0.5_DP*(rmloc(1,iloc)-rmloc(2,iloc))*flux(iedge)/hstep
          END IF
          
          ! Get corresponding matrix indices
          ik = Ksep(i); jk = Ksep(j)

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij
        END IF
      END IF
    END SUBROUTINE do_femtvd_assemble


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-GP in 1D
    SUBROUTINE do_femgp_1D(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
        IsubdiagonalEdgesIdx, IsubdiagonalEdges, DcoefficientsAtEdge,&
        Kld, Kcol, Cx, MC, u, u0, flux, flux0, pp, pm, qp, qm, rp, rm,&
        theta, tstep, hstep, NEQ, NEDGE, NNVEDGE, bextend, Ksep, Jac)
    
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN)  :: IverticesAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: IsubdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: IsubdiagonalEdges
      REAL(DP), DIMENSION(:,:), INTENT(IN)              :: DcoefficientsAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,MC
      REAL(DP), DIMENSION(:), INTENT(IN)                :: u,u0
      REAL(DP), DIMENSION(:), INTENT(IN)                :: flux,flux0
      REAL(DP), DIMENSION(:), INTENT(IN)                :: pp,pm,qp,qm,rp,rm
      REAL(DP), INTENT(IN)                              :: theta,tstep,hstep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER(PREC_MATIDX), INTENT(IN)                  :: NEDGE
      INTEGER, INTENT(IN)                               :: NNVEDGE
      LOGICAL, INTENT(IN)                               :: bextend
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: Jac
      
      ! local variables
      INTEGER(PREC_VECIDX), DIMENSION(5,NNVEDGE) :: Kloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: pploc,pmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: qploc,qmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: rploc,rmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: fluxloc,fluxloc0
      REAL(DP), DIMENSION(NDIM1D)      :: c_ij,c_ji
      
      INTEGER(PREC_MATIDX) :: ij,ji,ild,iedge
      INTEGER(PREC_VECIDX) :: i,j,k,l
      INTEGER :: iloc,nloc

      ! Loop over all columns of the Jacobian matrix
      DO k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        pploc(:,0) = 0; pmloc(:,0) = 0
        qploc(:,0) = 0; qmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        DO ild = IsubdiagonalEdgesIdx(k), IsubdiagonalEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagonalEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = Cx(ij)
          c_ji = Cx(ji)
          
          ! Update local coefficients
          CALL do_femgp_update(DcoefficientsAtEdge, MC, u, u0, flux, flux0,&
              pp, pm, qp, qm, c_ij, c_ji, theta, tstep, hstep,&
              iedge, i, j, ij, ji, NDIM1D, iloc, k, &
              pploc, pmloc, qploc, qmloc, fluxloc, fluxloc0, Kloc)
        END DO

        ! Loop over all superdiagonal edges
        DO iedge = IsuperdiagonalEdgesIdx(k), IsuperdiagonalEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                 
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = Cx(ij)
          c_ji = Cx(ji)
   
          ! Update local coefficients
          CALL do_femgp_update(DcoefficientsAtEdge, MC, u, u0, flux, flux0,&
              pp, pm, qp, qm, c_ij, c_ji, theta, tstep, hstep,&
              iedge, i, j, ij, ji, NDIM1D, iloc, k,&
              pploc, pmloc, qploc, qmloc, fluxloc, fluxloc0, Kloc)
        END DO
        
        ! Save total number of local neighbors
        nloc = iloc

        
        ! Adjust the diagonal separator
        CALL do_adjustKsep(Kld, Kcol, k, Ksep)

        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        rploc(:,0:nloc) = afcstab_limit( pploc(:,0:nloc), qploc(:,0:nloc), 0._DP, 1._DP)
        rmloc(:,0:nloc) = afcstab_limit(-pmloc(:,0:nloc),-qmloc(:,0:nloc), 0._DP, 1._DP)


        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        DO iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          DO ild = IsubdiagonalEdgesIdx(l), IsubdiagonalEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagonalEdges(ild)
            
            CALL do_femgp_assemble(IverticesAtEdge, Kld, Kcol, flux, flux0,&
                rp, rm, Kloc, rploc, rmloc, fluxloc, fluxloc0,&
                hstep, iedge, iloc, k, l, bextend, Ksep, Jac)
          END DO

          ! Loop over all superdiagonal edges
          DO iedge = IsuperdiagonalEdgesIdx(l), IsuperdiagonalEdgesIdx(l+1)-1
            
            CALL do_femgp_assemble(IverticesAtEdge, Kld, Kcol, flux, flux0,&
                rp, rm, Kloc, rploc, rmloc, fluxloc, fluxloc0,&
                hstep, iedge, iloc, k, l, bextend, Ksep, Jac)
          END DO
        END DO
      END DO   ! end-of k-loop
    END SUBROUTINE do_femgp_1D


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-GP in 2D
    SUBROUTINE do_femgp_2D(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
        IsubdiagonalEdgesIdx, IsubdiagonalEdges, DcoefficientsAtEdge,&
        Kld, Kcol, Cx, Cy, MC, u, u0, flux, flux0, pp, pm, qp, qm, rp, rm,&
        theta, tstep, hstep, NEQ, NEDGE, NNVEDGE, bextend, Ksep, Jac)
    
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN)  :: IverticesAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: IsubdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: IsubdiagonalEdges
      REAL(DP), DIMENSION(:,:), INTENT(IN)              :: DcoefficientsAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,MC
      REAL(DP), DIMENSION(:), INTENT(IN)                :: u,u0
      REAL(DP), DIMENSION(:), INTENT(IN)                :: flux,flux0
      REAL(DP), DIMENSION(:), INTENT(IN)                :: pp,pm,qp,qm,rp,rm
      REAL(DP), INTENT(IN)                              :: theta,tstep,hstep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER(PREC_MATIDX), INTENT(IN)                  :: NEDGE
      INTEGER, INTENT(IN)                               :: NNVEDGE
      LOGICAL, INTENT(IN)                               :: bextend
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: Jac
      
      ! local variables
      INTEGER(PREC_VECIDX), DIMENSION(5,NNVEDGE) :: Kloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: pploc,pmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: qploc,qmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: rploc,rmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: fluxloc,fluxloc0
      REAL(DP), DIMENSION(NDIM2D)      :: c_ij,c_ji
      
      INTEGER(PREC_MATIDX) :: ij,ji,ild,iedge
      INTEGER(PREC_VECIDX) :: i,j,k,l
      INTEGER :: iloc,nloc

      ! Loop over all columns of the Jacobian matrix
      DO k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        pploc(:,0) = 0; pmloc(:,0) = 0
        qploc(:,0) = 0; qmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        DO ild = IsubdiagonalEdgesIdx(k), IsubdiagonalEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagonalEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/Cx(ij),Cy(ij)/)
          c_ji = (/Cx(ji),Cy(ji)/)
          
          ! Update local coefficients
          CALL do_femgp_update(DcoefficientsAtEdge, MC, u, u0, flux, flux0,&
              pp, pm, qp, qm, c_ij, c_ji, theta, tstep, hstep,&
              iedge, i, j, ij, ji, NDIM2D, iloc, k, &
              pploc, pmloc, qploc, qmloc, fluxloc, fluxloc0, Kloc)
        END DO

        ! Loop over all superdiagonal edges
        DO iedge = IsuperdiagonalEdgesIdx(k), IsuperdiagonalEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                 
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/Cx(ij),Cy(ij)/)
          c_ji = (/Cx(ji),Cy(ji)/)
   
          ! Update local coefficients
          CALL do_femgp_update(DcoefficientsAtEdge, MC, u, u0, flux, flux0,&
              pp, pm, qp, qm, c_ij, c_ji, theta, tstep, hstep,&
              iedge, i, j, ij, ji, NDIM2D, iloc, k,&
              pploc, pmloc, qploc, qmloc, fluxloc, fluxloc0, Kloc)
        END DO
        
        ! Save total number of local neighbors
        nloc = iloc

        
        ! Adjust the diagonal separator
        CALL do_adjustKsep(Kld, Kcol, k, Ksep)

        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        rploc(:,0:nloc) = afcstab_limit( pploc(:,0:nloc), qploc(:,0:nloc), 0._DP, 1._DP)
        rmloc(:,0:nloc) = afcstab_limit(-pmloc(:,0:nloc),-qmloc(:,0:nloc), 0._DP, 1._DP)


        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        DO iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          DO ild = IsubdiagonalEdgesIdx(l), IsubdiagonalEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagonalEdges(ild)
            
            CALL do_femgp_assemble(IverticesAtEdge, Kld, Kcol, flux, flux0,&
                rp, rm, Kloc, rploc, rmloc, fluxloc, fluxloc0,&
                hstep, iedge, iloc, k, l, bextend, Ksep, Jac)
          END DO

          ! Loop over all superdiagonal edges
          DO iedge = IsuperdiagonalEdgesIdx(l), IsuperdiagonalEdgesIdx(l+1)-1
            
            CALL do_femgp_assemble(IverticesAtEdge, Kld, Kcol, flux, flux0,&
                rp, rm, Kloc, rploc, rmloc, fluxloc, fluxloc0,&
                hstep, iedge, iloc, k, l, bextend, Ksep, Jac)
          END DO
        END DO
      END DO   ! end-of k-loop
    END SUBROUTINE do_femgp_2D


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-GP in 3D
    SUBROUTINE do_femgp_3D(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
        IsubdiagonalEdgesIdx, IsubdiagonalEdges, DcoefficientsAtEdge,&
        Kld, Kcol, Cx, Cy, Cz, MC, u, u0, flux, flux0, pp, pm, qp, qm, rp, rm,&
        theta, tstep, hstep, NEQ, NEDGE, NNVEDGE, bextend, Ksep, Jac)
    
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN)  :: IverticesAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: IsubdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: IsubdiagonalEdges
      REAL(DP), DIMENSION(:,:), INTENT(IN)              :: DcoefficientsAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      REAL(DP), DIMENSION(:), INTENT(IN)                :: Cx,Cy,Cz,MC
      REAL(DP), DIMENSION(:), INTENT(IN)                :: u,u0
      REAL(DP), DIMENSION(:), INTENT(IN)                :: flux,flux0
      REAL(DP), DIMENSION(:), INTENT(IN)                :: pp,pm,qp,qm,rp,rm
      REAL(DP), INTENT(IN)                              :: theta,tstep,hstep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER(PREC_MATIDX), INTENT(IN)                  :: NEDGE
      INTEGER, INTENT(IN)                               :: NNVEDGE
      LOGICAL, INTENT(IN)                               :: bextend
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: Jac
      
      ! local variables
      INTEGER(PREC_VECIDX), DIMENSION(5,NNVEDGE) :: Kloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: pploc,pmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: qploc,qmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: rploc,rmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: fluxloc,fluxloc0
      REAL(DP), DIMENSION(NDIM3D)      :: c_ij,c_ji
      
      INTEGER(PREC_MATIDX) :: ij,ji,ild,iedge
      INTEGER(PREC_VECIDX) :: i,j,k,l
      INTEGER :: iloc,nloc

      ! Loop over all columns of the Jacobian matrix
      DO k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        pploc(:,0) = 0; pmloc(:,0) = 0
        qploc(:,0) = 0; qmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        DO ild = IsubdiagonalEdgesIdx(k), IsubdiagonalEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagonalEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/Cx(ij),Cy(ij),Cz(ij)/)
          c_ji = (/Cx(ji),Cy(ji),Cz(ji)/)
          
          ! Update local coefficients
          CALL do_femgp_update(DcoefficientsAtEdge, MC, u, u0, flux, flux0,&
              pp, pm, qp, qm, c_ij, c_ji, theta, tstep, hstep,&
              iedge, i, j, ij, ji, NDIM3D, iloc, k, &
              pploc, pmloc, qploc, qmloc, fluxloc, fluxloc0, Kloc)
        END DO

        ! Loop over all superdiagonal edges
        DO iedge = IsuperdiagonalEdgesIdx(k), IsuperdiagonalEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                 
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/Cx(ij),Cy(ij),Cz(ij)/)
          c_ji = (/Cx(ji),Cy(ji),Cz(ji)/)
   
          ! Update local coefficients
          CALL do_femgp_update(DcoefficientsAtEdge, MC, u, u0, flux, flux0,&
              pp, pm, qp, qm, c_ij, c_ji, theta, tstep, hstep,&
              iedge, i, j, ij, ji, NDIM3D, iloc, k,&
              pploc, pmloc, qploc, qmloc, fluxloc, fluxloc0, Kloc)
        END DO
        
        ! Save total number of local neighbors
        nloc = iloc

        
        ! Adjust the diagonal separator
        CALL do_adjustKsep(Kld, Kcol, k, Ksep)

        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        rploc(:,0:nloc) = afcstab_limit( pploc(:,0:nloc), qploc(:,0:nloc), 0._DP, 1._DP)
        rmloc(:,0:nloc) = afcstab_limit(-pmloc(:,0:nloc),-qmloc(:,0:nloc), 0._DP, 1._DP)


        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        DO iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          DO ild = IsubdiagonalEdgesIdx(l), IsubdiagonalEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagonalEdges(ild)
            
            CALL do_femgp_assemble(IverticesAtEdge, Kld, Kcol, flux, flux0,&
                rp, rm, Kloc, rploc, rmloc, fluxloc, fluxloc0,&
                hstep, iedge, iloc, k, l, bextend, Ksep, Jac)
          END DO

          ! Loop over all superdiagonal edges
          DO iedge = IsuperdiagonalEdgesIdx(l), IsuperdiagonalEdgesIdx(l+1)-1
            
            CALL do_femgp_assemble(IverticesAtEdge, Kld, Kcol, flux, flux0,&
                rp, rm, Kloc, rploc, rmloc, fluxloc, fluxloc0,&
                hstep, iedge, iloc, k, l, bextend, Ksep, Jac)
          END DO
        END DO
      END DO   ! end-of k-loop
    END SUBROUTINE do_femgp_3D

    
    !**************************************************************
    ! Update the local coefficients for FEM-GP in arbitrary dimensions
    SUBROUTINE do_femgp_update(DcoefficientsAtEdge, MC, u, u0, flux, flux0,&
        pp, pm, qp, qm, c_ij, c_ji, theta, tstep, hstep, &
        iedge, i, j, ij, ji, ndim, iloc, k,&
        pploc, pmloc, qploc, qmloc, fluxloc, fluxloc0, Kloc)
      
      REAL(DP), DIMENSION(:,:), INTENT(IN)                 :: DcoefficientsAtEdge
      REAL(DP), DIMENSION(:), INTENT(IN)                   :: MC
      REAL(DP), DIMENSION(:), INTENT(IN)                   :: u,u0
      REAL(DP), DIMENSION(:), INTENT(IN)                   :: flux,flux0
      REAL(DP), DIMENSION(:), INTENT(IN)                   :: pp,pm,qp,qm
      REAL(DP), DIMENSION(:), INTENT(IN)                   :: c_ij,c_ji
      REAL(DP), INTENT(IN)                                 :: theta,tstep,hstep
      INTEGER(PREC_MATIDX), INTENT(IN)                     :: iedge
      INTEGER(PREC_VECIDX), INTENT(IN)                     :: i,j
      INTEGER(PREC_MATIDX), INTENT(IN)                     :: ij,ji
      INTEGER, INTENT(IN)                                  :: ndim
      INTEGER, INTENT(IN)                                  :: iloc
      INTEGER(PREC_VECIDX), INTENT(IN)                     :: k

      ! We actually know, that all local quantities start at index zero
      REAL(DP), DIMENSION(:,0:), INTENT(INOUT)             :: pploc,pmloc
      REAL(DP), DIMENSION(:,0:), INTENT(INOUT)             :: qploc,qmloc
      REAL(DP), DIMENSION(:,0:), INTENT(INOUT)             :: fluxloc,fluxloc0
      INTEGER(PREC_VECIDX), DIMENSION(:,:), INTENT(INOUT)  :: Kloc

      ! local variables
      REAL(DP), DIMENSION(ndim) :: v_ij,v_ji
      REAL(DP) :: m_ij,d_ij,df_ij,f_ij,l_ij,l_ji,p_ij,pf_ij,q_ij,q_ji
      REAL(DP) :: diff,diff1,diff0,hstep_ik,hstep_jk,dsign
      INTEGER  :: iperturb

      
      !------------------------------------------------------------
      ! (1) unperturbed values: Retrieve the global Ps and Qs and
      !     copy their content to the local ones. Moreover,
      !     "eliminate" the contribution of the edge IJ for the
      !     unperturbed solution values u_i and u_j.
      !------------------------------------------------------------
      ! Determine coefficients
      d_ij = DcoefficientsAtEdge(1,iedge)
      l_ij = DcoefficientsAtEdge(2,iedge)
      l_ji = DcoefficientsAtEdge(3,iedge)
      
      ! Include consistent mass matrix
      m_ij = MC(ij)
      q_ij = m_ij/tstep+l_ij
      q_ji = m_ij/tstep+l_ji

      ! Determine solution differences
      diff1 = u(i)-u(j)
      diff0 = u0(i)-u0(j)

      ! Determine total solution difference
      diff = tstep*(theta*diff1+(1-theta)*diff0)

      ! Compute antidiffusive flux 
      IF(ABS(diff) < SYS_EPSREAL) THEN
        p_ij = 0
        f_ij = 0
      ELSE
        p_ij = MAX(0._DP, m_ij*(diff1-diff0)/diff+d_ij)
        f_ij = p_ij*diff
      END IF

      ! Prelimit the antidiffusive flux
      pf_ij = MIN(p_ij, l_ji)*diff
      
      ! Compute the remaining flux
      df_ij = f_ij-pf_ij

      IF (i .EQ. k) THEN
        
        ! Store global node number of the opposite node
        Kloc(1,iloc) = j

        ! Compute signed perturbation parameters
        hstep_ik = hstep; hstep_jk = 0._DP
        
        ! Update nodal coefficients for vertex j (!) which is the downwind node
        pploc(:,iloc) = pp(j)-MAX(0._DP,-df_ij)
        pmloc(:,iloc) = pm(j)-MIN(0._DP,-df_ij)
        qploc(:,iloc) = qp(j)-MAX(0._DP, diff)*q_ji
        qmloc(:,iloc) = qm(j)-MIN(0._DP, diff)*q_ji

      ELSE

        ! Store global node number of the opposite node
        Kloc(1,iloc) = i

        ! Compute signed perturbation parameters
        hstep_ik = 0._DP; hstep_jk = hstep
        
        ! Update nodal coefficients for vertex i (!) which is the upwind node
        pploc(:,iloc) = pp(i)-MAX(0._DP, f_ij)
        pmloc(:,iloc) = pm(i)-MIN(0._DP, f_ij)
        qploc(:,iloc) = qp(i)-MAX(0._DP,-diff)*q_ij
        qmloc(:,iloc) = qm(i)-MIN(0._DP,-diff)*q_ij
      END IF

      !------------------------------------------------------------
      ! (2) perturbed values: Now, the local Ps and Qs still
      !     require the contribution of the perturbed solution
      !     values u +/- h*e_k, whereby e_k denotes the k-th unit
      !     vector and h stands for the perturbation step length
      !------------------------------------------------------------
      
      DO iperturb = 1, 2
        
        ! Compute correct sign of perturbation
        dsign = -2*iperturb+3
        
        ! Compute perturbed velocity
        CALL fcb_getVelocity(u(i)+dsign*hstep_ik, u(j)+dsign*hstep_jk, i, j, v_ij, v_ji)

        ! Compute perturbation coefficients k_ij and k_ji
        l_ij = -SUM(v_ij*c_ij)
        l_ji = -SUM(v_ji*c_ji)
        
        ! Compute diffusion coefficient
        d_ij = MAX(-l_ij, 0._DP, -l_ji)

        ! Perform discrete upwinding
        l_ij = l_ij+d_ij
        l_ji = l_ji+d_ij

        q_ij = m_ij/tstep+l_ij
        q_ji = m_ij/tstep+l_ji

        ! Due to the (possible) nonlinearity of the velocity vector
        ! the orientation convention for the edge ij may be violated,
        ! that is, the condition 0=l_ij < l_ji may not be valid. In this
        ! case the node number i and j must be swapped logically
        IF (l_ij .LE. l_ji) THEN
          
          ! Save oriented node numbers
          Kloc(2*iperturb:2*iperturb+1,iloc) = (/i,j/)
          
          ! Update solution difference
          diff1 = u(i)-u(j)+dsign*(hstep_ik-hstep_jk)

          ! Update total solution difference
          diff = tstep*(theta*diff1+(1-theta)*diff0)

          ! Compute antidiffusive flux
          IF (ABS(diff) < SYS_EPSREAL) THEN
            p_ij = 0
            f_ij = 0
          ELSE
            p_ij = MAX(0._DP,m_ij*(diff1-diff0)/diff+d_ij)
            f_ij = p_ij*diff
          END IF
          
          ! Prelimit the antidiffusive flux
          pf_ij = MIN(p_ij,l_ji)*diff
          fluxloc0(iperturb,iloc) = pf_ij
          
          ! Compute the remaining flux
          df_ij = f_ij-pf_ij
          fluxloc(iperturb,iloc) = df_ij
        
          IF (i .EQ. k) THEN

            ! For node k which is the upwind node
            pploc(iperturb,0) = pploc(iperturb,0)+MAX(0._DP, f_ij)
            pmloc(iperturb,0) = pmloc(iperturb,0)+MIN(0._DP, f_ij)
            qploc(iperturb,0) = qploc(iperturb,0)+MAX(0._DP,-diff)*q_ij
            qmloc(iperturb,0) = qmloc(iperturb,0)+MIN(0._DP,-diff)*q_ij
            
            ! For node l opposite to k which is the downwind node
            pploc(iperturb,iloc) = pploc(iperturb,iloc)+MAX(0._DP,-df_ij)
            pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+MIN(0._DP,-df_ij)
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+MAX(0._DP, diff)*q_ji
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+MIN(0._DP, diff)*q_ji

          ELSE

            ! For node k which is the downwind node
            pploc(iperturb,0) = pploc(iperturb,0)+MAX(0._DP,-df_ij)
            pmloc(iperturb,0) = pmloc(iperturb,0)+MIN(0._DP,-df_ij)
            qploc(iperturb,0) = qploc(iperturb,0)+MAX(0._DP, diff)*q_ji
            qmloc(iperturb,0) = qmloc(iperturb,0)+MIN(0._DP, diff)*q_ji
            
            ! For node l opposite to k
            pploc(iperturb,iloc) = pploc(iperturb,iloc)+MAX(0._DP, f_ij)
            pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+MIN(0._DP, f_ij)
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+MAX(0._DP,-diff)*q_ij
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+MIN(0._DP,-diff)*q_ij

          END IF
          
        ELSE
          
          ! Save oriented node numbers
          Kloc(2*iperturb:2*iperturb+1,iloc) = (/j,i/)
          
          ! Update solution difference
          diff1 = u(i)-u(j)+dsign*(hstep_ik-hstep_jk)
          
          ! Update total solution difference
          diff = tstep*(theta*diff1+(1._DP-theta)*diff0)

          ! Compute antidiffusive flux
          IF (ABS(diff) < SYS_EPSREAL) THEN
            p_ij = 0
            f_ij = 0
          ELSE
            p_ij = MAX(0._DP, m_ij*(diff1-diff0)/diff+d_ij)
            f_ij = -p_ij*diff
          END IF

          ! Prelimit the antidiffusive flux
          pf_ij = -MIN(p_ij,l_ij)*diff
          fluxloc0(iperturb,iloc) = pf_ij

          ! Compute the remaining flux
          df_ij = f_ij-pf_ij
          fluxloc(iperturb,iloc) = df_ij
          
          IF (j .EQ. k) THEN
            
            ! For node k which is the upwind node
            pploc(iperturb,0) = pploc(iperturb,0)+MAX(0._DP, f_ij)
            pmloc(iperturb,0) = pmloc(iperturb,0)+MIN(0._DP, f_ij)
            qploc(iperturb,0) = qploc(iperturb,0)+MAX(0._DP, diff)*q_ij
            qmloc(iperturb,0) = qmloc(iperturb,0)+MIN(0._DP, diff)*q_ij
                       
            ! For node "l" opposite to k which is the downwind node
            pploc(iperturb,iloc) = pploc(iperturb,iloc)+MAX(0._DP,-df_ij)
            pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+MIN(0._DP,-df_ij)
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+MAX(0._DP,-diff)*q_ji
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+MIN(0._DP,-diff)*q_ji

          ELSE

            ! For node k which is the downwind node
            pploc(iperturb,0) = pploc(iperturb,0)+MAX(0._DP,-df_ij)
            pmloc(iperturb,0) = pmloc(iperturb,0)+MIN(0._DP,-df_ij)
            qploc(iperturb,0) = qploc(iperturb,0)+MAX(0._DP,-diff)*q_ji
            qmloc(iperturb,0) = qmloc(iperturb,0)+MIN(0._DP,-diff)*q_ji
            
            ! For node "l" opposite to k which is the upwind node
            pploc(iperturb,iloc) = pploc(iperturb,iloc)+MAX(0._DP, f_ij)
            pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+MIN(0._DP, f_ij)
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+MAX(0._DP, diff)*q_ij
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+MIN(0._DP, diff)*q_ij

          END IF
        END IF
      END DO
    END SUBROUTINE do_femgp_update


    !**************************************************************
    ! Assemble the given column of the Jacobian for FEM-GP for arbitrary dimensions
    SUBROUTINE do_femgp_assemble(IverticesAtEdge, Kdiagonal, Kcol,&
        flux, flux0, rp, rm, Kloc, rploc, rmloc, fluxloc, fluxloc0,&
        hstep, iedge, iloc, k, l, bextend, Ksep, Jac)

      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN)  :: IverticesAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      REAL(DP), DIMENSION(:), INTENT(IN)                :: flux,flux0
      REAL(DP), DIMENSION(:), INTENT(IN)                :: rp,rm
      INTEGER(PREC_VECIDX), DIMENSION(:,:), INTENT(IN)  :: Kloc
      REAL(DP), DIMENSION(:,0:), INTENT(IN)             :: rploc,rmloc
      REAL(DP), DIMENSION(:,0:), INTENT(IN)             :: fluxloc,fluxloc0
      REAL(DP), INTENT(IN)                              :: hstep
      INTEGER(PREC_MATIDX), INTENT(IN)                  :: iedge
      INTEGER, INTENT(IN)                               :: iloc
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: k,l
      LOGICAL, INTENT(IN)                               :: bextend

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: Jac

      ! local variables
      INTEGER(PREC_MATIDX) :: ik,jk
      INTEGER(PREC_VECIDX) :: i,j,m
      REAL(DP)             :: f_ij,pf_ij,df_ij
      INTEGER              :: iperturb
      
      ! Get global node number for edge IJ and the 
      ! number of the node m which is not l
      i = IverticesAtEdge(1,iedge)
      j = IverticesAtEdge(2,iedge)
      m = (i+j)-l
      
      ! We need to find out, which kind of edge is processed
      IF (m .EQ. k) THEN

        !------------------------------------------------------------
        ! 1. Case: primary edge
        !------------------------------------------------------------
        ! The current edge connects the perturbed node k with its
        ! direct neighbor l. Hence, all required information can be
        ! extracted from the local arrays and no global data
        ! retrieval has to be performed.
        !
        ! (a) The edge orientation needs to be adjusted for each
        !     perturbation direction
        
        DO iperturb = 1, 2
          
          ! Retrieve precomputed fluxes
          df_ij = fluxloc(iperturb,iloc)
          pf_ij = fluxloc0(iperturb,iloc)

          ! Adjust edge orientation
          i = Kloc(2*iperturb,iloc)
          j = Kloc(2*iperturb+1,iloc)
          
          ! Which node is located upwind?
          IF (i .EQ. k) THEN
            
            ! Get corresponding matrix indices
            ik = Kdiagonal(i); jk = Ksep(j)
            
            ! Limit upwind contribution
            IF (pf_ij > 0.0_DP) THEN
              pf_ij = rploc(iperturb,0)*pf_ij
            ELSE
              pf_ij = rmloc(iperturb,0)*pf_ij
            END IF

            ! Limit symmetric contribution
            IF (df_ij > 0.0_DP) THEN
              df_ij = MIN(rploc(iperturb,0), rmloc(iperturb,iloc))*df_ij
            ELSE
              df_ij = MIN(rmloc(iperturb,0), rploc(iperturb,iloc))*df_ij
            END IF
            
          ELSE
            
            ! Get corresponding matrix indices
            jk = Kdiagonal(j); ik = Ksep(i)
            
            ! Limit upwind contribution
            IF (pf_ij > 0.0_DP) THEN
              pf_ij = rploc(iperturb,iloc)*pf_ij
            ELSE
              pf_ij = rmloc(iperturb,iloc)*pf_ij
            END IF

            ! Limit symmetric contribution
            IF (df_ij > 0.0_DP) THEN
              df_ij = MIN(rmloc(iperturb,0), rploc(iperturb,iloc))*df_ij
            ELSE
              df_ij = MIN(rploc(iperturb,0), rmloc(iperturb,iloc))*df_ij
            END IF
            
          END IF
          
          ! Combine both contributions and 
          ! adopt sign for perturbation direction
          f_ij = -(iperturb-1.5_DP)*(pf_ij+df_ij)/hstep

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij
        END DO
        
      ELSEIF (bextend) THEN
        
        !------------------------------------------------------------
        ! 2. Case: secondary edge
        !------------------------------------------------------------
        ! The current edge connects two nodes l and m which both are
        ! not equal to the perturbed vertex k. Thus, the influence of
        ! the solution perturbation can only be due to a change in
        ! the correction factors alpha_ij. Moreover, for upwind
        ! -biased flux limiting techniques only the nodal correction
        ! factors for the upwind node i is used. Hence, it suffices
        ! to check if node i corresponds to the modified vertex l.
        ! For the symmetric part, we must also check if node j
        ! correspond to the modified vertex l.
        
        IF (i .EQ. l) THEN

          ! Get precomputed fluxes
          pf_ij = flux0(iedge)
          df_ij = flux(iedge)

          ! Limit upwind contribution
          IF (pf_ij > 0.0_DP) THEN
            pf_ij = (rploc(1,iloc)-rploc(2,iloc))*pf_ij
          ELSE
            pf_ij = (rmloc(1,iloc)-rmloc(2,iloc))*pf_ij
          END IF

          ! Limit symmetric contribution
          IF (df_ij > 0.0_DP) THEN
            df_ij = (MIN(rploc(1,iloc), rm(j))-&
                     MIN(rploc(2,iloc), rm(j)))*df_ij
          ELSE
            df_ij = (MIN(rmloc(1,iloc), rp(j))-&
                     MIN(rmloc(2,iloc), rp(j)))*df_ij
          END IF

          ! Combine both contributions
          f_ij = 0.5_DP*(pf_ij+df_ij)/hstep

          ! Get corresponding matrix indices
          ik = Ksep(i); jk = Ksep(j)

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij

        ELSE

          ! Get precomputed flux (only symmetric part)
          df_ij = flux(iedge)

          ! Limit symmetric contribution
          IF (df_ij > 0.0_DP) THEN
            df_ij = (MIN(rp(i), rmloc(1,iloc))-&
                     MIN(rp(i), rmloc(2,iloc)))*df_ij
          ELSE
            df_ij = (MIN(rm(i), rploc(1,iloc))-&
                     MIN(rm(i), rploc(2,iloc)))*df_ij
          END IF

          ! Compute divided difference
          f_ij = 0.5_DP*df_ij/hstep

          ! Get corresponding matrix indices
          ik = Ksep(i); jk = Ksep(j)
          
          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij
          
        END IF
      END IF
    END SUBROUTINE do_femgp_assemble
  END SUBROUTINE gfsc_buildStabJacobianScalar_GPTVD

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_buildStabJacobianBlock_Symm(ru, dscale, hstep, bclear, rafcstab, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation part
    ! of the discrete diffusion operator for a scalar convection equation.
    ! Note that this routine serves as a wrapper for block vectors. If there
    ! is only one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
!</description>

!<input>
    ! solution vector
    TYPE(t_vectorBlock), INTENT(IN)                :: ru

    ! scaling parameter
    REAL(DP), INTENT(IN)                           :: dscale

    ! perturbation parameter
    REAL(DP), INTENT(IN)                           :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    LOGICAL, INTENT(IN)                            :: bclear

     ! callback functions to compute velocity
    INCLUDE 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)                 :: rafcstab

    ! Jacobian matrix
    TYPE(t_matrixScalar), INTENT(INOUT)            :: rmatrixJ   
!</inputoutput>
!</subroutine>

    IF (ru%nblocks  .NE. 1) THEN

      CALL output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacobianBlock_Symm')
      CALL sys_halt()

    ELSE

      CALL gfsc_buildStabJacobianScalar_Symm(ru%RvectorBlock(1), dscale,&
          hstep, bclear, rafcstab, rmatrixJ)

    END IF
  END SUBROUTINE gfsc_buildStabJacobianBlock_Symm

  !*****************************************************************************

!<subroutine>

  SUBROUTINE gfsc_buildStabJacobianScalar_Symm(ru, dscale, hstep, bclear, rafcstab, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation
    ! part of the discrete diffusion operator for a scalar convection equation.
!</description>

!<input>
    ! solution vector
    TYPE(t_vectorScalar), INTENT(IN)               :: ru

    ! scaling parameter
    REAL(DP), INTENT(IN)                           :: dscale

    ! perturbation parameter
    REAL(DP), INTENT(IN)                           :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    LOGICAL, INTENT(IN)                            :: bclear

     ! callback functions to compute velocity
    INCLUDE 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)                 :: rafcstab

    ! Jacobian matrix
    TYPE(t_matrixScalar), INTENT(INOUT)            :: rmatrixJ   
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_MATIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER   :: p_IsuperdiagonalEdgesIdx
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_IsubdiagonalEdges
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_IsubdiagonalEdgesIdx
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_Kld,p_Ksep
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_Kdiagonal
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER   :: p_Kcol
    REAL(DP), DIMENSION(:,:), POINTER             :: p_DcoefficientsAtEdge
    REAL(DP), DIMENSION(:), POINTER               :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    REAL(DP), DIMENSION(:), POINTER               :: p_flux
    REAL(DP), DIMENSION(:), POINTER               :: p_u
    REAL(DP), DIMENSION(:), POINTER               :: p_Jac
    INTEGER :: h_Ksep
    LOGICAL :: bextend

    ! Clear matrix?
    IF (bclear) CALL lsyssc_clearMatrix(rmatrixJ)
    
    ! What kind of stabilisation are we?
    SELECT CASE(rafcstab%ctypeAFCstabilisation)

    CASE(AFCSTAB_SYMMETRIC)

      ! Check if stabilisation is prepared
      IF (IAND(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)   .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_EDGEVALUES)      .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_ANTIDIFFUSION)   .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_BOUNDS)          .EQ. 0 .OR.&
          IAND(rafcstab%iSpec, AFCSTAB_FLUXES)          .EQ. 0) THEN
        CALL output_line('Stabilisation does not provide required structures',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacobianScalar_Symm')
        CALL sys_halt()
      END IF
      
      ! Check if subdiagonal edges need to be generated
      IF (IAND(rafcstab%iSpec, AFCSTAB_SUBDIAGONALEDGES) .EQ. 0)&
          CALL afcstab_generateSubdiagEdges(rafcstab)
      
      ! Set pointers
      CALL afcstab_getbase_IverticesAtEdge(rafcstab,  p_IverticesAtEdge)
      CALL afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
      CALL afcstab_getbase_IsubdiagEdge(rafcstab,    p_IsubdiagonalEdges)
      CALL afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagonalEdgesIdx)
      CALL afcstab_getbase_DcoeffsAtEdge(rafcstab,    p_DcoefficientsAtEdge)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
      CALL lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
      CALL lsyssc_getbase_double(rafcstab%RedgeVectors(1),  p_flux)
      CALL lsyssc_getbase_Kcol(rmatrixJ,   p_Kcol)
      CALL lsyssc_getbase_double(rmatrixJ, p_Jac)
      CALL lsyssc_getbase_double(ru,       p_u)
      
      ! What kind of matrix format are we?
      SELECT CASE(rmatrixJ%cmatrixFormat)
      CASE(LSYSSC_MATRIX7)
        
        ! Set pointers
        CALL lsyssc_getbase_Kld(rmatrixJ, p_Kld)
        
        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        CALL storage_copy(rmatrixJ%h_Kld, h_Ksep)
        CALL storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)
        CALL lalg_vectorAddScalarInt(p_Ksep, 1)

        ! Assembled extended Jacobian matrix
        bextend = (rafcstab%iextendedJacobian .NE. 0)

        CALL do_symmetric(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
            p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
            p_Kld, p_Kcol, p_u, p_flux, p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, dscale,&
            hstep, rafcstab%NEQ, rafcstab%NEDGE, rafcstab%NNVEDGE, bextend, p_Ksep, p_Jac)

        ! Free storage
        CALL storage_free(h_Ksep)
        
        
      CASE DEFAULT
        CALL output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacobianScalar_Symm')
        CALL sys_halt()
      END SELECT
      

    CASE DEFAULT
      CALL output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildStabJacobianScalar_Symm')
      CALL sys_halt()
    END SELECT

  CONTAINS
    
    ! Here, the working routine follow
    
    !**************************************************************    
    ! Adjust the diagonal separator.
    ! The separator is initialied by the column separator (increased
    ! by one if the is necessary for matrix format 7).
    ! Based on the matric structure given by Kld/Kcol, the separator
    ! is "moved" to the given column "k". For efficiency reasons, only
    ! those entries are considered which are present in column "k".
    SUBROUTINE do_adjustKsep(Kld, Kcol, k, Ksep)
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: k
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      
      INTEGER(PREC_MATIDX) :: ild,isep
      INTEGER(PREC_VECIDX) :: l
      INTEGER :: iloc

      ! Loop over all entries of the k-th row
      DO ild = Kld(k), Kld(k+1)-1
        
        ! Get the column number and the position of the separator
        l = Kcol(ild); isep = Ksep(l)

        ! If the separator does not point to the k-th column
        ! it must be adjusted accordingly
        IF (Kcol(isep) < k) Ksep(l) = Ksep(l)+1
      END DO
    END SUBROUTINE do_adjustKsep


    !**************************************************************
    ! Assemble the Jacobian matrix for symmetric flux limiting
    SUBROUTINE do_symmetric(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
        IsubdiagonalEdgesIdx, IsubdiagonalEdges, DcoefficientsAtEdge,&
        Kld, Kcol, u, flux, pp, pm, qp, qm, rp, rm, dscale, hstep,&
        NEQ, NEDGE, NNVEDGE, bextend, Ksep, Jac)
      
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: IsuperdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN)  :: IverticesAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: IsubdiagonalEdgesIdx
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: IsubdiagonalEdges
      REAL(DP), DIMENSION(:,:), INTENT(IN)              :: DcoefficientsAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      REAL(DP), DIMENSION(:), INTENT(IN)                :: u
      REAL(DP), DIMENSION(:), INTENT(IN)                :: flux
      REAL(DP), DIMENSION(:), INTENT(IN)                :: pp,pm,qp,qm,rp,rm
      REAL(DP), INTENT(IN)                              :: dscale,hstep
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ
      INTEGER(PREC_MATIDX), INTENT(IN)                  :: NEDGE
      INTEGER, INTENT(IN)                               :: NNVEDGE
      LOGICAL, INTENT(IN)                               :: bextend
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: Jac

      ! local variables
      INTEGER(PREC_VECIDX), DIMENSION(5,NNVEDGE) :: Kloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: pploc,pmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: qploc,qmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: rploc,rmloc
      REAL(DP), DIMENSION(2,0:NNVEDGE) :: fluxloc

      INTEGER(PREC_MATIDX) :: ild,iedge
      INTEGER(PREC_VECIDX) :: k,l
      INTEGER :: iloc,nloc

      ! Loop over all columns of the Jacobian matrix
      DO k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        pploc(:,0) = 0; pmloc(:,0) = 0
        qploc(:,0) = 0; qmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        DO ild = IsubdiagonalEdgesIdx(k), IsubdiagonalEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagonalEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Update local coefficients
          CALL do_symmetric_update(IverticesAtEdge, DcoefficientsAtEdge,&
              u, pp, pm, qp, qm, hstep, iedge, iloc, k,&
              pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        END DO

        ! Loop over all superdiagonal edges
        DO iedge = IsuperdiagonalEdgesIdx(k), IsuperdiagonalEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                    
          ! Update local coefficients
          CALL do_symmetric_update(IverticesAtEdge, DcoefficientsAtEdge,&
              u, pp, pm, qp, qm, hstep, iedge, iloc, k,&
              pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        END DO

        ! Save total number of local neighbors
        nloc = iloc
        
        ! Adjust the diagonal separator
        CALL do_adjustKsep(Kld, Kcol, k, Ksep)

        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        rploc(:,0:nloc) = afcstab_limit( pploc(:,0:nloc), qploc(:,0:nloc), 0._DP, 1._DP)
        rmloc(:,0:nloc) = afcstab_limit(-pmloc(:,0:nloc),-qmloc(:,0:nloc), 0._DP, 1._DP)

        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        DO iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          DO ild = IsubdiagonalEdgesIdx(l), IsubdiagonalEdgesIdx(l+1)-1
            
            ! Get edge number
            iedge = IsubdiagonalEdges(ild)
            
            CALL do_symmetric_assemble(IverticesAtEdge, Kld, Kcol, flux, rp, rm,&
                Kloc, rploc, rmloc, fluxloc, dscale, hstep,&
                iedge, iloc, k, l, bextend, Ksep, Jac)
          END DO
          
          ! Loop over all superdiagonal edges
          DO iedge = IsuperdiagonalEdgesIdx(l), IsuperdiagonalEdgesIdx(l+1)-1
            
            CALL do_symmetric_assemble(IverticesAtEdge, Kld, Kcol, flux, rp, rm,&
                Kloc, rploc, rmloc, fluxloc, dscale, hstep,&
                iedge, iloc, k, l, bextend, Ksep, Jac)
          END DO
        END DO
      END DO   ! end-of k-loop
    END SUBROUTINE do_symmetric

    
    !**************************************************************
    ! Update the local coefficients for symmetric flux limiting
    SUBROUTINE do_symmetric_update(IverticesAtEdge, DcoefficientsAtEdge,&
        u, pp, pm, qp, qm, hstep, iedge, iloc, k,&
        pploc, pmloc, qploc, qmloc, fluxloc, Kloc)

      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN)     :: IverticesAtEdge
      REAL(DP), DIMENSION(:,:), INTENT(IN)                 :: DcoefficientsAtEdge
      REAL(DP), DIMENSION(:), INTENT(IN)                   :: u
      REAL(DP), DIMENSION(:), INTENT(IN)                   :: pp,pm,qp,qm
      REAL(DP), INTENT(IN)                                 :: hstep
      INTEGER(PREC_MATIDX), INTENT(IN)                     :: iedge
      INTEGER, INTENT(IN)                                  :: iloc
      INTEGER(PREC_VECIDX), INTENT(IN)                     :: k

      ! We actually know, that all local quantities start at index zero
      REAL(DP), DIMENSION(:,0:), INTENT(INOUT)             :: pploc,pmloc
      REAL(DP), DIMENSION(:,0:), INTENT(INOUT)             :: qploc,qmloc
      REAL(DP), DIMENSION(:,0:), INTENT(INOUT)             :: fluxloc
      INTEGER(PREC_VECIDX), DIMENSION(:,:), INTENT(INOUT)  :: Kloc

      ! local variables
      INTEGER(PREC_VECIDX) :: i,j
      REAL(DP) :: d_ij,f_ij,s_ij,diff,hstep_ik,hstep_jk,dsign
      INTEGER  :: iperturb


      ! Determine indices. Obviously, either i or j must be equal
      ! to k. Otherwise, the edge ij would not be present in the
      ! list of incident edges for node k.
      i = IverticesAtEdge(1,iedge)
      j = IverticesAtEdge(2,iedge)

      !------------------------------------------------------------
      ! (1) unperturbed values: Retrieve the global Ps and Qs and
      !     copy their content to the local ones. Moreover,
      !     "eliminate" the contribution of the edge IJ for the
      !     unperturbed solution values u_i and u_j.
      !------------------------------------------------------------
      ! Determine coefficients
      d_ij = DcoefficientsAtEdge(1,iedge)
      s_ij = DcoefficientsAtEdge(2,iedge)

      ! Determine solution difference
      diff = u(i)-u(j)

      IF (i .EQ. k) THEN
        
        ! Store global node number of the opposite node
        Kloc(1,iloc) = j

        ! Compute signed perturbation parameters
        hstep_ik = hstep; hstep_jk = 0._DP
        
        ! Compute raw antidiffusve flux
        f_ij = d_ij*diff
        
        ! Update sums of raw antidiffusive fluxes
        pploc(:,iloc) = pp(j)-MAX(0._DP, -f_ij)
        pmloc(:,iloc) = pm(j)-MAX(0._DP, -f_ij)
        
        ! Compute admissible edge contribution
        f_ij = -s_ij*diff
        
        ! Update upper/lower bounds
        qploc(:,iloc) = qp(j)-MAX(0._DP, -f_ij)
        qmloc(:,iloc) = qm(j)-MIN(0._DP, -f_ij)
        
      ELSE
        
        ! Store global node number of the opposite node
        Kloc(1,iloc) = i
        
        ! Compute signed perturbation parameters
        hstep_ik = 0._DP; hstep_jk = hstep
        
        ! Compute raw antidiffusve flux
        f_ij = d_ij*diff
        
        ! Update sums of raw antidiffusive fluxes
        pploc(:,iloc) = pp(i)-MAX(0._DP, f_ij)
        pmloc(:,iloc) = pm(i)-MIN(0._DP, f_ij)
        
        ! Compute admissible edge contribution
        f_ij = -s_ij*diff
        
        ! Update upper/lower bounds
        qploc(:,iloc) = qp(i)-MAX(0._DP, f_ij)
        qmloc(:,iloc) = qm(i)-MIN(0._DP, f_ij)
      END IF
      
      !------------------------------------------------------------
      ! (2) perturbed values: Now, the local Ps and Qs still
      !     require the contribution of the perturbed solution
      !     values u +/- h*e_k, whereby e_k denotes the k-th unit
      !     vector and h stands for the perturbation step length
      !------------------------------------------------------------
      
      !------------------------------------------------------------
      ! (3) perform the perturbation for "+/-h*e_k"
      !------------------------------------------------------------
        
      DO iperturb = 1, 2
        
        ! Compute correct sign of perturbation
        dsign = -2*iperturb+3  
        
        ! Save local node numbers
        Kloc(2*iperturb:2*iperturb+1,iloc) = (/i,j/)
        
        IF (i .EQ. k) THEN
          
          ! Compute raw antidiffusve flux
          f_ij = d_ij*(diff+dsign*(hstep_ik-hstep_jk))
          fluxloc(iperturb,iloc) = f_ij

          ! Update sums of raw antidiffusive fluxes
          pploc(iperturb,0)    = pploc(iperturb,0)+MAX(0._DP, f_ij)
          pmloc(iperturb,0)    = pmloc(iperturb,0)+MIN(0._DP, f_ij)
          pploc(iperturb,iloc) = pploc(iperturb,iloc)+MAX(0._DP, -f_ij)
          pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+MIN(0._DP, -f_ij)

          ! Compute admissible edge contribution
          f_ij = -s_ij*(diff+dsign*(hstep_ik-hstep_jk))

          ! Update upper/lower bounds
          qploc(iperturb,0)    = qploc(iperturb,0)+MAX(0._DP, f_ij)
          qmloc(iperturb,0)    = qmloc(iperturb,0)+MIN(0._DP, f_ij)
          qploc(iperturb,iloc) = qploc(iperturb,iloc)+MAX(0._DP, -f_ij)
          qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+MIN(0._DP, -f_ij)
          
        ELSE
          
          ! Compute raw antidiffusve flux
          f_ij = d_ij*(diff+dsign*(hstep_ik-hstep_jk))
          fluxloc(iperturb,iloc) = f_ij

          ! Update sums of raw antidiffusive fluxes
          pploc(iperturb,iloc) = pploc(iperturb,iloc)+MAX(0._DP, f_ij)
          pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+MIN(0._DP, f_ij)
          pploc(iperturb,0)    = pploc(iperturb,0)+MAX(0._DP, -f_ij)
          pmloc(iperturb,0)    = pmloc(iperturb,0)+MIN(0._DP, -f_ij)

          ! Compute admissible edge contribution
          f_ij = -s_ij*(diff+dsign*(hstep_ik-hstep_jk))
          
          ! Update upper/lower bounds
          qploc(iperturb,iloc) = qploc(iperturb,iloc)+MAX(0._DP, f_ij)
          qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+MIN(0._DP, f_ij)
          qploc(iperturb,0)    = qploc(iperturb,0)+MAX(0._DP, -f_ij)
          qmloc(iperturb,0)    = qmloc(iperturb,0)+MIN(0._DP, -f_ij)
        END IF
      END DO
    END SUBROUTINE do_symmetric_update

    
    !**************************************************************
    ! Assemble the given column of the Jacobian for symmetric flux limiting
    SUBROUTINE do_symmetric_assemble(IverticesAtEdge, Kdiagonal, Kcol,&
        flux, rp, rm, Kloc, rploc, rmloc, fluxloc, dscale, &
        hstep, iedge, iloc, k, l, bextend, Ksep, Jac)
      
      INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN)  :: IverticesAtEdge
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)    :: Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)    :: Kcol
      REAL(DP), DIMENSION(:), INTENT(IN)                :: flux
      REAL(DP), DIMENSION(:), INTENT(IN)                :: rp,rm
      INTEGER(PREC_VECIDX), DIMENSION(:,:), INTENT(IN)  :: Kloc
      REAL(DP), DIMENSION(:,0:), INTENT(IN)             :: rploc,rmloc
      REAL(DP), DIMENSION(:,0:), INTENT(IN)             :: fluxloc
      REAL(DP), INTENT(IN)                              :: dscale,hstep
      INTEGER(PREC_MATIDX), INTENT(IN)                  :: iedge
      INTEGER, INTENT(IN)                               :: iloc
      INTEGER(PREC_VECIDX), INTENT(IN)                  :: k,l
      LOGICAL, INTENT(IN)                               :: bextend

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Ksep
      REAL(DP), DIMENSION(:), INTENT(INOUT)             :: Jac
      
      ! local variables
      INTEGER(PREC_MATIDX) :: ik,jk
      INTEGER(PREC_VECIDX) :: i,j,m
      REAL(DP)             :: f_ij
      INTEGER              :: iperturb

      ! Get global node number for edge IJ and the 
      ! number of the node m which is not l
      i = IverticesAtEdge(1,iedge)
      j = IverticesAtEdge(2,iedge)
      m = (i+j)-l
      
      ! We need to find out, which kind of edge is processed
      IF (m .EQ. k) THEN

        !------------------------------------------------------------
        ! 1. Case: primary edge
        !------------------------------------------------------------
        ! The current edge connects the perturbed node k with its
        ! direct neighbor l. Hence, all required information can be
        ! extracted from the local arrays and no global data
        ! retrieval has to be performed.
        !
        ! (a) The edge orientation needs to be adjusted for each
        !     perturbation direction
        
        DO iperturb = 1, 2
          
          ! Retrieve precomputed flux
          f_ij = fluxloc(iperturb,iloc)

          ! Adjust edge orientation
          i = Kloc(2*iperturb,iloc)
          j = Kloc(2*iperturb+1,iloc)
          
          ! Which node is located upwind?
          IF (i .EQ. k) THEN
            
            ! Get corresponding matrix indices
            ik = Kdiagonal(i); jk = Ksep(j)
            
            ! Limit flux 
            IF (f_ij > 0.0_DP) THEN
              f_ij = dscale*MIN(rploc(iperturb,0), rmloc(iperturb,iloc))*f_ij
            ELSE
              f_ij = dscale*MIN(rmloc(iperturb,0), rploc(iperturb,iloc))*f_ij
            END IF
            
          ELSE
            
            ! Get corresponding matrix indices
            jk = Kdiagonal(j); ik = Ksep(i)
            
            ! Limit flux
            IF (f_ij > 0.0_DP) THEN
              f_ij = dscale*MIN(rploc(iperturb,iloc), rmloc(iperturb,0))*f_ij
            ELSE
              f_ij = dscale*MIN(rmloc(iperturb,iloc), rploc(iperturb,0))*f_ij
            END IF
            
          END IF
          
          ! Adopt sign for perturbation direction
          f_ij = -(iperturb-1.5_DP)*f_ij/hstep

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij
        END DO
        
      ELSEIF (bextend) THEN
        
        !------------------------------------------------------------
        ! 2. Case: secondary edge
        !------------------------------------------------------------
        ! The current edge connects two nodes l and m which both are
        ! not equal to the perturbed vertex k. Thus, the influence of
        ! the solution perturbation can only be due to a change in
        ! the correction factors alpha_ij. Moreover, for upwind
        ! -biased flux limiting techniques only the nodal correction
        ! factors for the upwind node i is used. Hence, it suffices
        ! to check if node i corresponds to the modified vertex l.
        ! Interestingly enough, some edge LM which connects two
        ! direct neighbors of the perturbed vertex k is only
        ! processed once due to the fact that either l or (!) m
        ! corresponds to the upwind node.

        IF (i .EQ. l) THEN

          IF (flux(iedge) > 0.0_DP) THEN
            f_ij = 0.5_DP*dscale*(MIN(rploc(1,iloc), rm(j))-&
                                  MIN(rploc(2,iloc), rm(j)))*flux(iedge)/hstep
          ELSE
            f_ij = 0.5_DP*dscale*(MIN(rmloc(1,iloc), rp(j))-&
                                  MIN(rmloc(2,iloc), rp(j)))*flux(iedge)/hstep
          END IF
          
          ! Get corresponding matrix indices
          ik = Ksep(i); jk = Ksep(j)

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij

        ELSE

          IF (flux(iedge) > 0.0_DP) THEN
            f_ij = 0.5_DP*dscale*(MIN(rp(i), rmloc(1,iloc))-&
                                  MIN(rp(i), rmloc(2,iloc)))*flux(iedge)/hstep
          ELSE
            f_ij = 0.5_DP*dscale*(MIN(rm(i), rploc(1,iloc))-&
                                  MIN(rm(i), rploc(2,iloc)))*flux(iedge)/hstep
          END IF
          
          ! Get corresponding matrix indices
          ik = Ksep(i); jk = Ksep(j)

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij

        END IF

      END IF
    END SUBROUTINE do_symmetric_assemble
  END SUBROUTINE gfsc_buildStabJacobianScalar_Symm
END MODULE groupfemscalar
