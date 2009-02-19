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
!# 6.) gfsys_buildResidualTVD = gfsys_buildResScalarTVD /
!#                              gfsys_buildResBlockTVD
!#     -> assemble the residual for FEM-TVD stabilisation
!#
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

module groupfemsystem

  use afcstabilisation
  use fsystem
  use genoutput
  use linearsystemblock
  use linearsystemscalar
  use storage

  implicit none
  
  private

  public :: gfsys_initStabilisation
  public :: gfsys_isMatrixCompatible
  public :: gfsys_isVectorCompatible
  public :: gfsys_buildDivOperator
  public :: gfsys_buildResidual
  public :: gfsys_buildResidualTVD

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<constants>
!<constantblock description="Global constants for directional splitting">

  ! unit vector in X-direction in 1D
  real(DP), dimension(NDIM1D) :: XDir1D = (/ 1.0_DP /)

  ! unit vector in X-direction in 2D
  real(DP), dimension(NDIM2D) :: XDir2D = (/ 1.0_DP, 0.0_DP /)

  ! unit vector in Y-direction in 2D
  real(DP), dimension(NDIM2D) :: YDir2D = (/ 0.0_DP, 1.0_DP /)

  ! unit vector in X-direction in 3D
  real(DP), dimension(NDIM3D) :: XDir3D = (/ 1.0_DP, 0.0_DP, 0.0_DP /)

  ! unit vector in Y-direction in 3D
  real(DP), dimension(NDIM3D) :: YDir3D = (/ 0.0_DP, 1.0_DP, 0.0_DP /)

  ! unit vector in Z-direction in 3D
  real(DP), dimension(NDIM3D) :: ZDir3D = (/ 0.0_DP, 0.0_DP, 1.0_DP /)

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
  private :: t_array
  type t_array

    ! Pointer to the double-valued matrix data
    real(DP), dimension(:), pointer :: Da

    ! Pointer to the single-valued matrix data
    real(SP), dimension(:), pointer :: Fa

    ! Pointer to the integer-valued matrix data
    integer, dimension(:), pointer :: Ia
  end type t_array
!</typeblock>
!</types>

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************
  
  interface gfsys_initStabilisation
    module procedure gfsys_initStabilisationScalar
    module procedure gfsys_initStabilisationBlock
  end interface

  interface gfsys_isMatrixCompatible
    module procedure gfsys_isMatrixCompatibleScalar
    module procedure gfsys_isMatrixCompatibleBlock
  end interface

  interface gfsys_isVectorCompatible
    module procedure gfsys_isVectorCompatibleScalar
    module procedure gfsys_isVectorCompatibleBlock
  end interface

  interface gfsys_buildDivOperator
     module procedure gfsys_buildDivOperatorScalar
     module procedure gfsys_buildDivOperatorBlock
  end interface

  interface gfsys_buildResidual
    module procedure gfsys_buildResScalar
    module procedure gfsys_buildResBlock
  end interface

  interface gfsys_buildResidualTVD
    module procedure gfsys_buildResScalarTVD
    module procedure gfsys_buildResBlockTVD
  end interface

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

contains

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_initStabilisationBlock(rmatrixBlockTemplate, rafcstab)

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
    type(t_matrixBlock), intent(IN) :: rmatrixBlockTemplate
!</input>

!<inputoutput>
    ! discrete operator structure
    type(t_afcstab), intent(INOUT) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i,j


    ! Check if block matrix has only one block
    if ((rmatrixBlockTemplate%nblocksPerCol .eq. 1) .and. &
        (rmatrixBlockTemplate%nblocksPerRow .eq. 1)) then
      call gfsys_initStabilisationScalar(&
          rmatrixBlockTemplate%RmatrixBlock(1,1), rafcstab)

      ! That's it
      return
    end if
    
    ! Check that number of columns equans number of rows
    if (rmatrixBlockTemplate%nblocksPerCol .ne. &
        rmatrixBlockTemplate%nblocksPerRow) then
      call output_line('Block matrix must have equal number of columns and rows!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_initStabilisationBlock')
      call sys_halt()
    end if
    
    ! Check if matrix exhibits group structure
    if (rmatrixBlockTemplate%imatrixSpec .ne. LSYSBS_MSPEC_GROUPMATRIX) then
      call output_line('Block matrix must have group structure!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_initStabilisationBlock')
      call sys_halt()
    end if
    
    ! Set atomic data from first block
    rafcstab%NVAR  = rmatrixBlockTemplate%nblocksPerCol
    rafcstab%NEQ   = rmatrixBlockTemplate%RmatrixBlock(1,1)%NEQ
    rafcstab%NEDGE = int(0.5*(rmatrixBlockTemplate%RmatrixBlock(1,1)%NA-&
                              rmatrixBlockTemplate%RmatrixBlock(1,1)%NEQ))

    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)
      
    case (AFCSTAB_GALERKIN,&
          AFCSTAB_UPWIND)
      ! do nothing


    case (AFCSTAB_FEMTVD)

      ! We need 6 nodal vectors for P's, Q's and R's
      allocate(rafcstab%RnodalVectors(6))
      do i = 1, 6
        call lsyssc_createVector(rafcstab%RnodalVectors(i), rafcstab%NEQ,&
                                 rafcstab%NVAR, .false., ST_DOUBLE)
      end do
      

    case (AFCSTAB_FEMFCT)

      ! We need 7 nodal vectors for P's, Q's and R's
      allocate(rafcstab%RnodalVectors(7))
      do i = 1, 7
        call lsyssc_createVector(rafcstab%RnodalVectors(i), rafcstab%NEQ,&
                                 rafcstab%NVAR, .false., ST_DOUBLE)
      end do

      ! We need 2 edgewise vectors for the fluxes
      allocate(rafcstab%RedgeVectors(4))
      do i = 1, 4
        call lsyssc_createVector(rafcstab%RedgeVectors(i), rafcstab%NEDGE,&
                                 rafcstab%NVAR, .false., ST_DOUBLE)
      end do


    case DEFAULT
      call output_line('Invalid type of stabilisation!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_initStabilisationBlock')
      call sys_halt()
    end select
    
    ! Set specifier
    rafcstab%iSpec = AFCSTAB_INITIALISED

  end subroutine gfsys_initStabilisationBlock

  ! *****************************************************************************

  subroutine gfsys_initStabilisationScalar(rmatrixTemplate, rafcstab)

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
    type(t_matrixScalar), intent(IN) :: rmatrixTemplate
!</input>

!<inputoutput>
    ! discrete operator structure
    type(t_afcstab), intent(INOUT) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i

  
    ! Set atomic data
    rafcstab%NVAR  = rmatrixTemplate%NVAR
    rafcstab%NEQ   = rmatrixTemplate%NEQ
    rafcstab%NEDGE = int(0.5*(rmatrixTemplate%NA-rmatrixTemplate%NEQ))

    
    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)
      
    case (AFCSTAB_GALERKIN,&
          AFCSTAB_UPWIND)
      ! do nothing

      
    case (AFCSTAB_FEMTVD)

      ! We need 6 nodal vectors for P's, Q's and R's
      allocate(rafcstab%RnodalVectors(6))
      do i = 1, 6
        call lsyssc_createVector(rafcstab%RnodalVectors(i), rafcstab%NEQ,&
                                 rafcstab%NVAR, .false., ST_DOUBLE)
      end do


    case (AFCSTAB_FEMFCT)

      ! We need 6 nodal vectors for P's, Q's and R's
      allocate(rafcstab%RnodalVectors(7))
      do i = 1, 7
        call lsyssc_createVector(rafcstab%RnodalVectors(i), rafcstab%NEQ,&
                                 rafcstab%NVAR, .false., ST_DOUBLE)
      end do

      ! We need 4 edgewise vectors for the fluxes
      allocate(rafcstab%RedgeVectors(4))
      do i = 1, 4
        call lsyssc_createVector(rafcstab%RedgeVectors(i), rafcstab%NEDGE,&
                                 rafcstab%NVAR, .false., ST_DOUBLE)
      end do


    case DEFAULT
      call output_line('Invalid type of stabilisation!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_initStabilisationScalar')
      call sys_halt()
    end select
    
    ! Set specifier
    rafcstab%iSpec = AFCSTAB_INITIALISED

  end subroutine gfsys_initStabilisationScalar

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_isMatrixCompatibleBlock(rafcstab, rmatrixBlock, bcompatible)

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
    type(t_matrixBlock), intent(IN) :: rmatrixBlock

    ! stabilisation structure
    type(t_afcstab), intent(IN) :: rafcstab
!</input>

!<output>
    ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
    ! whether matrix and stabilisation are compatible or not.
    ! If not given, an error will inform the user if the matrix/operator are
    ! not compatible and the program will halt.
    logical, intent(OUT), optional :: bcompatible
!</output>
!</subroutine>

    
    ! Check if matrix has only one block
    if ((rmatrixBlock%nblocksPerCol .eq. 1) .and. &
        (rmatrixBlock%nblocksPerRow .eq. 1)) then
      call gfsys_isMatrixCompatibleScalar(rafcstab,&
          rmatrixBlock%RmatrixBlock(1,1), bcompatible)
      return
    end if

    ! Check that number of columns equans number of rows
    if (rmatrixBlock%nblocksPerCol .ne. &
        rmatrixBlock%nblocksPerRow) then
      call output_line('Block matrix must have equal number of columns and rows!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_isMatrixCompatibleBlock')
      call sys_halt()
    end if

    ! Check if matrix exhibits group structure
    if (rmatrixBlock%imatrixSpec .ne. LSYSBS_MSPEC_GROUPMATRIX) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Block matrix must have group structure!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'gfsys_isMatrixCompatibleBlock')
        call sys_halt()
      end if
    end if

    ! Matrix/operator must have the same size
    if (rafcstab%NVAR  .ne. rmatrixBlock%nblocksPerCol .or.&
        rafcstab%NEQ   .ne. rmatrixBlock%RmatrixBlock(1,1)%NEQ  .or.&
        rafcstab%NEDGE .ne. int(0.5*(rmatrixBlock%RmatrixBlock(1,1)%NA-&
                                     rmatrixBlock%RmatrixBlock(1,1)%NEQ))) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Matrix/Operator not compatible, different structure!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'gfssy_isMatrixCompatibleBlock')
        call sys_halt()
      end if
    end if
  end subroutine gfsys_isMatrixCompatibleBlock

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_isMatrixCompatibleScalar(rafcstab, rmatrix, bcompatible)

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
    type(t_matrixScalar), intent(IN) :: rmatrix

    ! stabilisation structure
    type(t_afcstab), intent(IN) :: rafcstab
!</input>

!<output>
    ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
    ! whether matrix and stabilisation are compatible or not.
    ! If not given, an error will inform the user if the matrix/operator are
    ! not compatible and the program will halt.
    logical, intent(OUT), optional :: bcompatible
!</output>
!</subroutine>

    ! Matrix/operator must have the same size
    if (rafcstab%NEQ   .ne. rmatrix%NEQ  .or.&
        rafcstab%NVAR  .ne. rmatrix%NVAR .or.&
        rafcstab%NEDGE .ne. int(0.5*(rmatrix%NA-rmatrix%NEQ))) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Matrix/Operator not compatible, different structure!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'gfssy_isMatrixCompatibleScalar')
        call sys_halt()
      end if
    end if
  end subroutine gfsys_isMatrixCompatibleScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_isVectorCompatibleBlock(rafcstab, rvectorBlock, bcompatible)

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
    type(t_vectorBlock), intent(IN) :: rvectorBlock

    ! stabilisation structure
    type(t_afcstab), intent(IN) :: rafcstab
!</input>

!<output>
    ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
    ! whether matrix and stabilisation are compatible or not.
    ! If not given, an error will inform the user if the matrix/operator are
    ! not compatible and the program will halt.
    logical, intent(OUT), optional :: bcompatible
!</output>
!</subroutine>

    ! Check if block vectors has just one block
    if (rvectorBlock%nblocks .eq. 1) then
      call gfsys_isVectorCompatibleScalar(rafcstab,&
          rvectorBlock%RvectorBlock(1), bcompatible)
      return
    end if

    ! Vector/operator must have the same size
    if (rafcstab%NEQ*rafcstab%NVAR .ne. rvectorBlock%NEQ) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Vector/Operator not compatible, different structure!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'gfsys_isVectorCompatibleBlock')
        call sys_halt()
      end if
    end if
  end subroutine gfsys_isVectorCompatibleBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_isVectorCompatibleScalar(rafcstab, rvector, bcompatible)

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
    type(t_vectorScalar), intent(IN) :: rvector

    ! stabilisation structure
    type(t_afcstab), intent(IN) :: rafcstab
!</input>

!<output>
    ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
    ! whether matrix and stabilisation are compatible or not.
    ! If not given, an error will inform the user if the matrix/operator are
    ! not compatible and the program will halt.
    logical, intent(OUT), optional :: bcompatible
!</output>
!</subroutine>

    ! Matrix/operator must have the same size
    if (rafcstab%NEQ   .ne. rvector%NEQ .or.&
        rafcstab%NVAR  .ne. rvector%NVAR) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Vector/Operator not compatible, different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_isVectorCompatibleScalar')
        call sys_halt()
      end if
    end if
  end subroutine gfsys_isVectorCompatibleScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_buildDivOperatorBlock(RcoeffMatrices, ru, fcb_calcMatrix,&
                                         dscale, bclear, rdivMatrix)

!<description>
    ! This subroutine assembles the discrete transport operator $K$ so that
    !
    !   $$ (Ku)_i = -\sum_{j\ne i}\bf c}_{ij}\cdot({\bf F}_j-{\bf F}_i) $$
    !
    ! whereby the flux difference can be represented as
    !
    !   $$ {\bf F}_j-{\bf F}_i=\hat{\bf A}_{ij}(u_j-u_i) $$
    !
    ! and the matrix $\hat{\bf A}$ corresponds to the Jacobian tensor
    ! evaluated for some special set of averaged flow variables.
    !
    ! Note that this routine is designed for block matrices/vectors. 
    ! If there is only one block, then the corresponding scalar routine 
    ! is called. Otherwise, the global operator is treated as block matrix.
    ! This block matrix has to be in group structure, that is, the structure
    ! of subblock(1,1) will serve as template for all other submatrices.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(IN) :: RcoeffMatrices

    ! solution vector
    type(t_vectorBlock), intent(IN) :: ru

    ! scaling factor
    real(DP), intent(IN) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear

    ! callback functions to compute local Roe matrix and
    ! local dissipation matrix
    include 'intf_gfsyscallback.inc'
!</input>   

!<inputoutput>
    ! global transport operator
    type(t_matrixBlock), intent(INOUT) :: rdivMatrix
!</inputoutput>
!</subroutine>
    
    ! local variables
    type(t_array), dimension(ru%nblocks,ru%nblocks)  :: rarray
    real(DP), dimension(:), pointer :: p_Cx,p_Cy,p_Cz,p_u
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    integer :: h_Ksep,ndim
    logical :: bisFullMatrix

    ! Check if block vector contains only one block and if
    ! global operator is stored in interleave format.
    if ((ru%nblocks .eq. 1) .and.&
        (rdivMatrix%nblocksPerCol .eq. 1) .and. &
        (rdivMatrix%nblocksPerRow .eq. 1)) then
      call gfsys_buildDivOperatorScalar(RcoeffMatrices, ru%RvectorBlock(1), fcb_calcMatrix,&
                                        dscale, bclear, rdivMatrix%RmatrixBlock(1,1))
      return       
    end if

    ! Check if block matrix exhibits group structure
    if (rdivMatrix%imatrixSpec .ne. LSYSBS_MSPEC_GROUPMATRIX) then
      call output_line('Block matrix must have group structure!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
      call sys_halt()
    end if
    
    ! Clear matrix?
    if (bclear) call lsysbl_clearMatrix(rdivMatrix)
    
    ! Set pointers
    call gfsys_getbase_double(rdivMatrix, rarray, bisFullMatrix)
    call lsysbl_getbase_double(ru, p_u)
    
    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)

    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_Cz)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
      call sys_halt()
    end select
    
    ! What kind of matrix are we?
    select case(RcoeffMatrices(1)%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld(RcoeffMatrices(1), p_Kld)
      call lsyssc_getbase_Kcol(RcoeffMatrices(1), p_Kcol)

      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(RcoeffMatrices(1)%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, RcoeffMatrices(1)%NEQ+1)
      
      ! What type of matrix are we?
      if (bisFullMatrix) then
        
        select case(ndim)
        case (NDIM1D)
          call doOperatorMat7_1D(p_Kld, p_Kcol, p_Ksep, RcoeffMatrices(1)%NEQ,&
                                 ru%nblocks, p_Cx, p_u, dscale, rarray)
        case (NDIM2D)
          call doOperatorMat7_2D(p_Kld, p_Kcol, p_Ksep, RcoeffMatrices(1)%NEQ,&
                                 ru%nblocks, p_Cx, p_Cy, p_u, dscale, rarray)
        case (NDIM3D)
          call doOperatorMat7_3D(p_Kld, p_Kcol, p_Ksep, RcoeffMatrices(1)%NEQ,&
                                 ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, dscale, rarray)
        end select
        
      else   ! bisFullMatrix

        select case(ndim)
        case (NDIM1D)
          call doOperatorMat7Diag_1D(p_Kld, p_Kcol, p_Ksep, RcoeffMatrices(1)%NEQ,&
                                     ru%nblocks, p_Cx, p_u, dscale, rarray)
        case (NDIM2D)
          call doOperatorMat7Diag_2D(p_Kld, p_Kcol, p_Ksep, RcoeffMatrices(1)%NEQ,&
                                     ru%nblocks, p_Cx, p_Cy, p_u, dscale, rarray)
        case (NDIM3D)
          call doOperatorMat7Diag_3D(p_Kld, p_Kcol, p_Ksep, RcoeffMatrices(1)%NEQ,&
                                     ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, dscale, rarray)
        end select

      end if   ! bisFullMatrix
      
      ! Release diagonal separator
      call storage_free(h_Ksep)
      
      
    case(LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 9
      !-------------------------------------------------------------------------

      ! Set pointers
      call lsyssc_getbase_Kld(RcoeffMatrices(1), p_Kld)
      call lsyssc_getbase_Kcol(RcoeffMatrices(1), p_Kcol)
      call lsyssc_getbase_Kdiagonal(RcoeffMatrices(1), p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(RcoeffMatrices(1)%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, RcoeffMatrices(1)%NEQ+1)

      ! What type of matrix are we?
      if (bisFullMatrix) then
        
        select case(ndim)
        case (NDIM1D)
          call doOperatorMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                 RcoeffMatrices(1)%NEQ, ru%nblocks,&
                                 p_Cx, p_u, dscale, rarray)
        case (NDIM2D)
          call doOperatorMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                 RcoeffMatrices(1)%NEQ, ru%nblocks,&
                                 p_Cx, p_Cy, p_u, dscale, rarray)
        case (NDIM3D)
          call doOperatorMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                 RcoeffMatrices(1)%NEQ, ru%nblocks,&
                                 p_Cx, p_Cy, p_Cz, p_u, dscale, rarray)
        end select
        
      else   ! bisFullMatrix

        select case(ndim)
        case (NDIM1D)
          call doOperatorMat9Diag_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                     RcoeffMatrices(1)%NEQ, ru%nblocks,&
                                     p_Cx, p_u, dscale, rarray)
        case (NDIM2D)
          call doOperatorMat9Diag_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                     RcoeffMatrices(1)%NEQ, ru%nblocks,&
                                     p_Cx, p_Cy, p_u, dscale, rarray)
        case (NDIM3D)
          call doOperatorMat9Diag_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                     RcoeffMatrices(1)%NEQ, ru%nblocks,&
                                     p_Cx, p_Cy, p_Cz, p_u, dscale, rarray)
        end select
        
      end if   ! bisFullMatrix
      
      ! Release diagonal separator
      call storage_free(h_Ksep)
      

    case DEFAULT
      call output_line('Unsupported matrix format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the working routines follow
    
    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 1D
    ! All matrices are stored in matrix format 7
    
    subroutine doOperatorMat7Diag_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
                                     Cx, u, dscale, K)
      
      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NVAR

      type(t_array), dimension(:,:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: A_ij,S_ij,u_i,u_j
      integer :: ii,ij,ji,jj,i,j,ivar

      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute matrices
          call fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)

          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij
          
          ! Assemble the global operator
          do ivar = 1, NVAR
            K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii) + A_ij(ivar) + S_ij(ivar)
            K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij) - A_ij(ivar) - S_ij(ivar)
            K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji) + A_ij(ivar) - S_ij(ivar)
            K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj) - A_ij(ivar) + S_ij(ivar)
          end do
        end do
      end do
    end subroutine doOperatorMat7Diag_1D

    
    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 2D
    ! All matrices are stored in matrix format 7
    
    subroutine doOperatorMat7Diag_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
                                     Cx, Cy, u, dscale, K)
      
      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NVAR

      type(t_array), dimension(:,:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: A_ij,S_ij,u_i,u_j
      integer :: ii,ij,ji,jj,i,j,ivar

      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute matrices
          call fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)
          
          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij

          ! Assemble the global operator
          do ivar = 1, NVAR
            K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii) + A_ij(ivar) + S_ij(ivar)
            K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij) - A_ij(ivar) - S_ij(ivar)
            K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji) + A_ij(ivar) - S_ij(ivar)
            K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj) - A_ij(ivar) + S_ij(ivar)
          end do
        end do
      end do
    end subroutine doOperatorMat7Diag_2D


    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 3D
    ! All matrices are stored in matrix format 7
    
    subroutine doOperatorMat7Diag_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
                                     Cx, Cy, Cz, u, dscale, K)
      
      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NVAR

      type(t_array), dimension(:,:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: A_ij,S_ij,u_i,u_j
      integer :: ii,ij,ji,jj,i,j,ivar
      

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
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
          call fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)
          
          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij

          ! Assemble the global operator
          do ivar = 1, NVAR
            K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii) + A_ij(ivar) + S_ij(ivar)
            K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij) - A_ij(ivar) - S_ij(ivar)
            K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji) + A_ij(ivar) - S_ij(ivar)
            K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj) - A_ij(ivar) + S_ij(ivar)
          end do
        end do
      end do
    end subroutine doOperatorMat7Diag_3D


    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 1D
    ! All matrices are stored in matrix format 9
    
    subroutine doOperatorMat9Diag_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
                                     Cx, u, dscale, K)
      
      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR

      type(t_array), dimension(:,:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: A_ij,S_ij,u_i,u_j
      integer :: ii,ij,ji,jj,i,j,ivar

      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute matrices
          call fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)

          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij

          ! Assemble the global operator
          do ivar = 1, NVAR
            K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii) + A_ij(ivar) + S_ij(ivar)
            K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij) - A_ij(ivar) - S_ij(ivar)
            K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji) + A_ij(ivar) - S_ij(ivar)
            K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj) - A_ij(ivar) + S_ij(ivar)
          end do
        end do
      end do
    end subroutine doOperatorMat9Diag_1D


    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 2D
    ! All matrices are stored in matrix format 9
    
    subroutine doOperatorMat9Diag_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
                                     Cx, Cy, u, dscale, K)
      
      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR

      type(t_array), dimension(:,:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: A_ij,S_ij,u_i,u_j
      integer :: ii,ij,ji,jj,i,j,ivar


      ! Loop over all rows
      do i = 1, NEQ

        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute matrices
          call fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)

          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij
          
          ! Assemble the global operator
          do ivar = 1, NVAR
            K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii) + A_ij(ivar) + S_ij(ivar)
            K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij) - A_ij(ivar) - S_ij(ivar)
            K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji) + A_ij(ivar) - S_ij(ivar)
            K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj) - A_ij(ivar) + S_ij(ivar)
          end do
        end do
      end do
    end subroutine doOperatorMat9Diag_2D


    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 2D
    ! All matrices are stored in matrix format 9
    
    subroutine doOperatorMat9Diag_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
                                     Cx, Cy, Cz, u, dscale, K)
      
      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR

      type(t_array), dimension(:,:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: A_ij,S_ij,u_i,u_j
      integer :: ii,ij,ji,jj,i,j,ivar


      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
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
          call fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)

          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij

          ! Assemble the global operator
          do ivar = 1, NVAR
            K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii) + A_ij(ivar) + S_ij(ivar)
            K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij) - A_ij(ivar) - S_ij(ivar)
            K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji) + A_ij(ivar) - S_ij(ivar)
            K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj) - A_ij(ivar) + S_ij(ivar)
          end do
        end do
      end do
    end subroutine doOperatorMat9Diag_3D


    !**************************************************************
    ! Assemble divergence operator K in 1D
    ! All matrices are stored in matrix format 7
    
    subroutine doOperatorMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
                                 Cx, u, dscale, K)
      
      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NVAR

      type(t_array), dimension(:,:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,S_ij
      real(DP), dimension(NVAR) :: u_i,u_j
      integer :: ii,ij,ji,jj,i,j,ivar,jvar,idx

      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute matrices
          call fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)
          
          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij
          
          ! Assemble the global operator
          do ivar = 1, NVAR
            do jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii) + A_ij(idx) + S_ij(idx)
              K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij) - A_ij(idx) - S_ij(idx)
              K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji) + A_ij(idx) - S_ij(idx)
              K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj) - A_ij(idx) + S_ij(idx)
            end do
          end do
        end do
      end do
    end subroutine doOperatorMat7_1D

    
    !**************************************************************
    ! Assemble divergence operator K in 2D
    ! All matrices are stored in matrix format 7
    
    subroutine doOperatorMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
                                 Cx, Cy, u, dscale, K)
      
      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NVAR

      type(t_array), dimension(:,:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,S_ij
      real(DP), dimension(NVAR) :: u_i,u_j
      integer :: ii,ij,ji,jj,i,j,ivar,jvar,idx
      

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute matrices
          call fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)

          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij

          ! Assemble the global operator
          do ivar = 1, NVAR
            do jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii) + A_ij(idx) + S_ij(idx)
              K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij) - A_ij(idx) - S_ij(idx)
              K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji) + A_ij(idx) - S_ij(idx)
              K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj) - A_ij(idx) + S_ij(idx)
            end do
          end do
        end do
      end do
    end subroutine doOperatorMat7_2D


    !**************************************************************
    ! Assemble divergence operator K in 3D
    ! All matrices are stored in matrix format 7
    
    subroutine doOperatorMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
                                 Cx, Cy, Cz, u, dscale, K)
      
      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NVAR

      type(t_array), dimension(:,:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,S_ij
      real(DP), dimension(NVAR) :: u_i,u_j
      integer :: ii,ij,ji,jj,i,j,ivar,jvar,idx
      

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
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
          call fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)
          
          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij

          ! Assemble the global operator
          do ivar = 1, NVAR
            do jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii) + A_ij(idx) + S_ij(idx)
              K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij) - A_ij(idx) - S_ij(idx)
              K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji) + A_ij(idx) - S_ij(idx)
              K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj) - A_ij(idx) + S_ij(idx)
            end do
          end do
        end do
      end do
    end subroutine doOperatorMat7_3D

    
    !**************************************************************
    ! Assemble divergence operator K in 1D
    ! All matrices are stored in matrix format 9
    
    subroutine doOperatorMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
                                 Cx, u, dscale, K)
      
      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR

      type(t_array), dimension(:,:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,S_ij
      real(DP), dimension(NVAR) :: u_i,u_j
      integer :: ii,ij,ji,jj,i,j,ivar,jvar,idx


      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute matrices
          call fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)

          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij

          ! Assemble the global operator
          do ivar = 1, NVAR
            do jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii) + A_ij(idx) + S_ij(idx)
              K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij) - A_ij(idx) - S_ij(idx)
              K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji) + A_ij(idx) - S_ij(idx)
              K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj) - A_ij(idx) + S_ij(idx)
            end do
          end do
        end do
      end do
    end subroutine doOperatorMat9_1D


    !**************************************************************
    ! Assemble divergence operator K in 2D
    ! All matrices are stored in matrix format 9
    
    subroutine doOperatorMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
                                 Cx, Cy, u, dscale, K)

      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR

      type(t_array), dimension(:,:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,S_ij
      real(DP), dimension(NVAR) :: u_i,u_j
      integer :: ii,ij,ji,jj,i,j,ivar,jvar,idx
      

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute matrices
          call fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)
          
          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij

          ! Assemble the global operator
          do ivar = 1, NVAR
            do jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii) + A_ij(idx) + S_ij(idx)
              K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij) - A_ij(idx) - S_ij(idx)
              K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji) + A_ij(idx) - S_ij(idx)
              K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj) - A_ij(idx) + S_ij(idx)
            end do
          end do
        end do
      end do
    end subroutine doOperatorMat9_2D


    !**************************************************************
    ! Assemble divergence operator K in 3D
    ! All matrices are stored in matrix format 9
    
    subroutine doOperatorMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
                                 Cx, Cy, Cz, u, dscale, K)
      
      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR

      type(t_array), dimension(:,:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,S_ij
      real(DP), dimension(NVAR) :: u_i,u_j
      integer :: ii,ij,ji,jj,i,j,ivar,jvar,idx


      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
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
          call fcb_calcMatrix(u_i, u_j, C_ij, C_ji, A_ij, S_ij)

          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij
          
          ! Assemble the global operator
          do ivar = 1, NVAR
            do jvar = 1, NVAR
              idx = NVAR*(ivar-1)+jvar
              K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii) + A_ij(idx) + S_ij(idx)
              K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij) - A_ij(idx) - S_ij(idx)
              K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji) + A_ij(idx) - S_ij(idx)
              K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj) - A_ij(idx) + S_ij(idx)
            end do
          end do
        end do
      end do
    end subroutine doOperatorMat9_3D
  end subroutine gfsys_buildDivOperatorBlock

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildDivOperatorScalar(RcoeffMatrices, ru, fcb_calcMatrix,&
                                          dscale, bclear, rdivMatrix)

!<description>
    ! This subroutine assembles the discrete transport operator $K$ so that
    !
    !   $$ (Ku)_i=-\sum_{j\ne i}\bf c}_{ij}\cdot({\bf F}_j-{\bf F}_i) $$
    !
    ! whereby the flux difference can be represented as
    !
    !   $$ {\bf F}_j-{\bf F}_i=\hat{\bf A}_{ij}(u_j-u_i) $$
    !
    ! and the matrix $\hat{\bf A}$ corresponds to the Jacobian tensor
    ! evaluated for some special set of averaged variables.
    !
    ! Note that this routine requires the scalar matrices/vectors are
    ! stored in the interleave format.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(IN) :: RcoeffMatrices
    
    ! scalar solution vector
    type(t_vectorScalar), intent(IN) :: ru

    ! scaling factor
    real(DP), intent(IN) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear

    ! callback functions to compute local Roe matrix and
    ! local dissipation matrix
    include 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! scalar transport operator
    type(t_matrixScalar), intent(INOUT) :: rdivMatrix
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Cx,p_Cy,p_Cz,p_DivOp,p_u
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    integer :: h_Ksep,ndim
    
    
    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rdivMatrix)
    
    ! Set pointers
    call lsyssc_getbase_double(rdivMatrix, p_DivOp)
    call lsyssc_getbase_double(ru, p_u)

    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)

    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_Cz)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
      call sys_halt()
    end select
    
    ! What kind of matrix are we?
    select case(rdivMatrix%cmatrixFormat)
    case(LSYSSC_MATRIX7INTL)
      !-------------------------------------------------------------------------
      ! Matrix format 7 interleaved
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld(RcoeffMatrices(1), p_Kld)
      call lsyssc_getbase_Kcol(RcoeffMatrices(1), p_Kcol)

      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(RcoeffMatrices(1)%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, RcoeffMatrices(1)%NEQ+1)
      
      ! What type of matrix are we?
      select case(rdivMatrix%cinterleavematrixFormat)
        
      case (LSYSSC_MATRIX1)
        
        select case(ndim)
        case (NDIM1D)
          call doOperatorMat7_1D(p_Kld, p_Kcol, p_Ksep,&
                                 RcoeffMatrices(1)%NEQ, RcoeffMatrices(1)%NA,&
                                 ru%NVAR, p_Cx, p_u, dscale, p_DivOp)
        case (NDIM2D)
          call doOperatorMat7_2D(p_Kld, p_Kcol, p_Ksep,&
                                 RcoeffMatrices(1)%NEQ, RcoeffMatrices(1)%NA,&
                                 ru%NVAR, p_Cx, p_Cy, p_u, dscale, p_DivOp)
        case (NDIM3D)
          call doOperatorMat7_3D(p_Kld, p_Kcol, p_Ksep,&
                                 RcoeffMatrices(1)%NEQ, RcoeffMatrices(1)%NA,&
                                 ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, dscale, p_DivOp)
        end select
        

      case (LSYSSC_MATRIXD)
        
        select case(ndim)
        case (NDIM1D)
          call doOperatorMat7Diag_1D(p_Kld, p_Kcol, p_Ksep,&
                                     RcoeffMatrices(1)%NEQ, RcoeffMatrices(1)%NA,&
                                     ru%NVAR, p_Cx, p_u, dscale, p_DivOp)
        case (NDIM2D)
          call doOperatorMat7Diag_2D(p_Kld, p_Kcol, p_Ksep,&
                                     RcoeffMatrices(1)%NEQ, RcoeffMatrices(1)%NA,&
                                     ru%NVAR, p_Cx, p_Cy, p_u, dscale, p_DivOp)
        case (NDIM3D)
          call doOperatorMat7Diag_3D(p_Kld, p_Kcol, p_Ksep,&
                                     RcoeffMatrices(1)%NEQ, RcoeffMatrices(1)%NA,&
                                     ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, dscale, p_DivOp)
        end select

        
      case DEFAULT
        call output_line('Unsupported interleave matrix format!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
        call sys_halt()
      end select
      
      ! Release diagonal separator
      call storage_free(h_Ksep)
      
      
    case (LSYSSC_MATRIX9INTL)
      !-------------------------------------------------------------------------
      ! Matrix format 9 interleaved
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld(RcoeffMatrices(1), p_Kld)
      call lsyssc_getbase_Kcol(RcoeffMatrices(1), p_Kcol)
      call lsyssc_getbase_Kdiagonal(rdivMatrix, p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(RcoeffMatrices(1)%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, RcoeffMatrices(1)%NEQ+1)
      
      ! What type of matrix are we?
      select case(rdivMatrix%cinterleavematrixFormat)
        
      case (LSYSSC_MATRIX1)
        
        select case(ndim)
        case (NDIM1D)
          call doOperatorMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                 RcoeffMatrices(1)%NEQ, RcoeffMatrices(1)%NA,&
                                 ru%NVAR, p_Cx, p_u, dscale, p_DivOp)
        case (NDIM2D)
          call doOperatorMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                 RcoeffMatrices(1)%NEQ, RcoeffMatrices(1)%NA,&
                                 ru%NVAR, p_Cx, p_Cy, p_u, dscale, p_DivOp)
        case (NDIM3D)
          call doOperatorMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                 RcoeffMatrices(1)%NEQ, RcoeffMatrices(1)%NA,&
                                 ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, dscale, p_DivOp)
        end select

        
      case (LSYSSC_MATRIXD)
        
        select case(ndim)
        case (NDIM1D)
          call doOperatorMat9Diag_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                     RcoeffMatrices(1)%NEQ, RcoeffMatrices(1)%NA,&
                                     ru%NVAR, p_Cx, p_u, dscale, p_DivOp)
        case (NDIM2D)
          call doOperatorMat9Diag_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                     RcoeffMatrices(1)%NEQ, RcoeffMatrices(1)%NA,&
                                     ru%NVAR, p_Cx, p_Cy, p_u, dscale, p_DivOp)
        case (NDIM3D)
          call doOperatorMat9Diag_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                     RcoeffMatrices(1)%NEQ, RcoeffMatrices(1)%NA,&
                                     ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, dscale, p_DivOp)
        end select
        

      case DEFAULT
        call output_line('Unsupported interleave matrix format!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
        call sys_halt()
      end select
      
      ! Release diagonal separator
      call storage_free(h_Ksep)
      
      
    case DEFAULT
      call output_line('Unsupported matrix format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
      call sys_halt()
    end select

  contains
    
    ! Here, the working routines follow
    
    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 1D
    ! All matrices are stored in matrix format 7
    
    subroutine doOperatorMat7Diag_1D(Kld, Kcol, Ksep, NEQ, NA, NVAR, &
                                     Cx, u, dscale, K)
      
      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NA,NVAR

      real(DP), dimension(NVAR,NA), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: A_ij,S_ij
      integer :: ii,ij,ji,jj,i,j
      
      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Compute matrices
          call fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)
          
          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij

          ! Assemble the global operator
          K(:,ii) = K(:,ii) + A_ij + S_ij
          K(:,ij) = K(:,ij) - A_ij - S_ij
          K(:,ji) = K(:,ji) + A_ij - S_ij
          K(:,jj) = K(:,jj) - A_ij + S_ij
        end do
      end do
    end subroutine doOperatorMat7Diag_1D


    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 2D
    ! All matrices are stored in matrix format 7
    
    subroutine doOperatorMat7Diag_2D(Kld, Kcol, Ksep, NEQ, NA, NVAR,&
                                     Cx, Cy, u, dscale, K)
      
      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NA,NVAR

      real(DP), dimension(NVAR,NA), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: A_ij,S_ij
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute matrices
          call fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)

          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij

          ! Assemble the global operator
          K(:,ii) = K(:,ii) + A_ij + S_ij
          K(:,ij) = K(:,ij) - A_ij - S_ij
          K(:,ji) = K(:,ji) + A_ij - S_ij
          K(:,jj) = K(:,jj) - A_ij + S_ij
        end do
      end do
    end subroutine doOperatorMat7Diag_2D


    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 3D
    ! All matrices are stored in matrix format 7
    
    subroutine doOperatorMat7Diag_3D(Kld, Kcol, Ksep, NEQ, NA, NVAR,&
                                     Cx, Cy, Cz, u, dscale, K)
      
      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NA,NVAR

      real(DP), dimension(NVAR,NA), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: A_ij,S_ij
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)
          
          ! Compute matrices
          call fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)
          
          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij

          ! Assemble the global operator
          K(:,ii) = K(:,ii) + A_ij + S_ij
          K(:,ij) = K(:,ij) - A_ij - S_ij
          K(:,ji) = K(:,ji) + A_ij - S_ij
          K(:,jj) = K(:,jj) - A_ij + S_ij
        end do
      end do
    end subroutine doOperatorMat7Diag_3D

    
    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 1D
    ! All matrices are stored in matrix format 9
    
    subroutine doOperatorMat9Diag_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NA, NVAR,&
                                     Cx, u, dscale, K)
      
      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NA,NVAR

      real(DP), dimension(NVAR,NA), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: A_ij,S_ij
      integer :: ii,ij,ji,jj,i,j
      
      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute matrices
          call fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)

          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij

          ! Assemble the global operator
          K(:,ii) = K(:,ii) + A_ij + S_ij
          K(:,ij) = K(:,ij) - A_ij - S_ij
          K(:,ji) = K(:,ji) + A_ij - S_ij
          K(:,jj) = K(:,jj) - A_ij + S_ij
        end do
      end do
    end subroutine doOperatorMat9Diag_1D


    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 2D
    ! All matrices are stored in matrix format 9
    
    subroutine doOperatorMat9Diag_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NA, NVAR,&
                                     Cx, Cy, u, dscale, K)
      
      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NA,NVAR

      real(DP), dimension(NVAR,NA), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: A_ij,S_ij
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute matrices
          call fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)

          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij
          
          ! Assemble the global operator
          K(:,ii) = K(:,ii) + A_ij + S_ij
          K(:,ij) = K(:,ij) - A_ij - S_ij
          K(:,ji) = K(:,ji) + A_ij - S_ij
          K(:,jj) = K(:,jj) - A_ij + S_ij
        end do
      end do
    end subroutine doOperatorMat9Diag_2D


    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 3D
    ! All matrices are stored in matrix format 9
    
    subroutine doOperatorMat9Diag_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NA, NVAR,&
                                     Cx, Cy, Cz, u, dscale, K)
      
      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NA,NVAR

      real(DP), dimension(NVAR,NA), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: A_ij,S_ij
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)
          
          ! Compute matrices
          call fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)

          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij

          ! Assemble the global operator
          K(:,ii) = K(:,ii) + A_ij + S_ij
          K(:,ij) = K(:,ij) - A_ij - S_ij
          K(:,ji) = K(:,ji) + A_ij - S_ij
          K(:,jj) = K(:,jj) - A_ij + S_ij
        end do
      end do
    end subroutine doOperatorMat9Diag_3D

    
    !**************************************************************
    ! Assemble divergence operator K in 1D
    ! All matrices are stored in matrix format 7
    
    subroutine doOperatorMat7_1D(Kld, Kcol, Ksep, NEQ, NA, NVAR,&
                                 Cx, u, dscale, K)
      
      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NA,NVAR

      real(DP), dimension(NVAR*NVAR,NA), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,S_ij
      integer :: ii,ij,ji,jj,i,j
      
      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute matrices
          call fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)

          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij

          ! Assemble the global operator
          K(:,ii) = K(:,ii) + A_ij + S_ij
          K(:,ij) = K(:,ij) - A_ij - S_ij
          K(:,ji) = K(:,ji) + A_ij - S_ij
          K(:,jj) = K(:,jj) - A_ij + S_ij
        end do
      end do
    end subroutine doOperatorMat7_1D

    
    !**************************************************************
    ! Assemble divergence operator K in 2D
    ! All matrices are stored in matrix format 7
    
    subroutine doOperatorMat7_2D(Kld, Kcol, Ksep, NEQ, NA, NVAR,&
                                 Cx, Cy, u, dscale, K)
      
      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NA,NVAR

      real(DP), dimension(NVAR*NVAR,NA), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,S_ij
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute matrices
          call fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)

          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij

          ! Assemble the global operator
          K(:,ii) = K(:,ii) + A_ij + S_ij
          K(:,ij) = K(:,ij) - A_ij - S_ij
          K(:,ji) = K(:,ji) + A_ij - S_ij
          K(:,jj) = K(:,jj) - A_ij + S_ij
        end do
      end do
    end subroutine doOperatorMat7_2D


    !**************************************************************
    ! Assemble divergence operator K in 3D
    ! All matrices are stored in matrix format 7
    
    subroutine doOperatorMat7_3D(Kld, Kcol, Ksep, NEQ, NA, NVAR,&
                                 Cx, Cy, Cz, u, dscale, K)
      
      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NA,NVAR

      real(DP), dimension(NVAR*NVAR,NA), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,S_ij
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)
          
          ! Compute matrices
          call fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)

          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij

          ! Assemble the global operator
          K(:,ii) = K(:,ii) + A_ij + S_ij
          K(:,ij) = K(:,ij) - A_ij - S_ij
          K(:,ji) = K(:,ji) + A_ij - S_ij
          K(:,jj) = K(:,jj) - A_ij + S_ij
        end do
      end do
    end subroutine doOperatorMat7_3D


    !**************************************************************
    ! Assemble divergence operator K in 1D
    ! All matrices are stored in matrix format 9
    
    subroutine doOperatorMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NA, NVAR,&
                                 Cx, u, dscale, K)
      
      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NA,NVAR

      real(DP), dimension(NVAR*NVAR,NA), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,S_ij
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute matrices
          call fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)

          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij

          ! Assemble the global operator
          K(:,ii) = K(:,ii) + A_ij + S_ij
          K(:,ij) = K(:,ij) - A_ij - S_ij
          K(:,ji) = K(:,ji) + A_ij - S_ij
          K(:,jj) = K(:,jj) - A_ij + S_ij
        end do
      end do
    end subroutine doOperatorMat9_1D


    !**************************************************************
    ! Assemble divergence operator K in 2D
    ! All matrices are stored in matrix format 9
    
    subroutine doOperatorMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NA, NVAR,&
                                 Cx, Cy, u, dscale, K)
      
      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NA,NVAR

      real(DP), dimension(NVAR*NVAR,NA), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,S_ij
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute matrices
          call fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)

          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij

          ! Assemble the global operator
          K(:,ii) = K(:,ii) + A_ij + S_ij
          K(:,ij) = K(:,ij) - A_ij - S_ij
          K(:,ji) = K(:,ji) + A_ij - S_ij
          K(:,jj) = K(:,jj) - A_ij + S_ij
        end do
      end do
    end subroutine doOperatorMat9_2D


    !**************************************************************
    ! Assemble divergence operator K in 3D
    ! All matrices are stored in matrix format 9
    
    subroutine doOperatorMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NA, NVAR,&
                                 Cx, Cy, Cz, u, dscale, K)
      
      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NA,NVAR

      real(DP), dimension(NVAR*NVAR,NA), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,S_ij
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)
          
          ! Compute matrices
          call fcb_calcMatrix(u(:,i), u(:,j), C_ij, C_ji, A_ij, S_ij)

          ! Scale matrix blocks
          A_ij = dscale*A_ij; S_ij = dscale*S_ij

          ! Assemble the global operator
          K(:,ii) = K(:,ii) + A_ij + S_ij
          K(:,ij) = K(:,ij) - A_ij - S_ij
          K(:,ji) = K(:,ji) + A_ij - S_ij
          K(:,jj) = K(:,jj) - A_ij + S_ij
        end do
      end do
    end subroutine doOperatorMat9_3D
  end subroutine gfsys_buildDivOperatorScalar

  ! *****************************************************************************

!<subroutine>
  
  subroutine gfsys_buildResBlock(RcoeffMatrices, ru, fcb_calcFlux, dscale, bclear, rres)

!<description>
    ! This subroutine assembles the residual vector for block vectors.
    ! If the vector contains only one block, then the scalar counterpart
    ! of this routine is called with the scalar subvector.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(IN) :: RcoeffMatrices

    ! solution vector
    type(t_vectorBlock), intent(IN) :: ru

    ! scaling factor
    real(DP), intent(IN) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(IN) :: bclear

    ! callback functions to compute local matrices
    include 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! residual vector
    type(t_vectorBlock), intent(INOUT) :: rres
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Cx,p_Cy,p_Cz,p_u,p_res
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    integer :: h_Ksep,ndim

    
    ! Check if block vectors contain only one block.
    if ((ru%nblocks .eq. 1) .and. (rres%nblocks .eq. 1) ) then
      call gfsys_buildResScalar(RcoeffMatrices, ru%RvectorBlock(1),&
          fcb_calcFlux, dscale, bclear, rres%RvectorBlock(1))
      return       
    end if
    
    ! Clear vector?
    if (bclear) call lsysbl_clearVector(rres)

    ! Set pointers
    call lsysbl_getbase_double(ru, p_u)
    call lsysbl_getbase_double(rres, p_res)
    
    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)

    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_Cz)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlock')
      call sys_halt()
    end select

    ! What kind of matrix are we?
    select case(RcoeffMatrices(1)%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld (RcoeffMatrices(1), p_Kld)
      call lsyssc_getbase_Kcol (RcoeffMatrices(1), p_Kcol)

      ! Create diagonal Ksep=Kld
      h_Ksep = ST_NOHANDLE
      call storage_copy(RcoeffMatrices(1)%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, RcoeffMatrices(1)%NEQ+1)

      ! How many dimensions do we have?
      select case(ndim)
      case (NDIM1D)
        call doResidualMat7_1D(p_Kld, p_Kcol, p_Ksep, RcoeffMatrices(1)%NEQ,&
                               ru%nblocks, p_Cx, p_u, dscale, p_res)
      case (NDIM2D)
        call doResidualMat7_2D(p_Kld, p_Kcol, p_Ksep, RcoeffMatrices(1)%NEQ,&
                               ru%nblocks, p_Cx, p_Cy, p_u, dscale, p_res)
      case (NDIM3D)
        call doResidualMat7_3D(p_Kld, p_Kcol, p_Ksep, RcoeffMatrices(1)%NEQ,&
                               ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, dscale, p_res)
      end select

      ! Release Ksep
      call storage_free(h_Ksep)
      

    case (LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 9
      !-------------------------------------------------------------------------

      ! Set pointer
      call lsyssc_getbase_Kld (RcoeffMatrices(1), p_Kld)
      call lsyssc_getbase_Kcol (RcoeffMatrices(1), p_Kcol)
      call lsyssc_getbase_Kdiagonal(RcoeffMatrices(1), p_Kdiagonal)

      ! Create diagonal Ksep=Kld
      h_Ksep = ST_NOHANDLE
      call storage_copy(RcoeffMatrices(1)%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, RcoeffMatrices(1)%NEQ+1)

      ! How many dimensions do we have?
      select case(ndim)
      case (NDIM1D)
        call doResidualMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                               RcoeffMatrices(1)%NEQ, ru%nblocks,&
                               p_Cx, p_u, dscale, p_res)
      case (NDIM2D)
         call doResidualMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                RcoeffMatrices(1)%NEQ, ru%nblocks,&
                                p_Cx, p_Cy, p_u, dscale, p_res)
       case (NDIM3D)
        call doResidualMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                               RcoeffMatrices(1)%NEQ, ru%nblocks,&
                               p_Cx, p_Cy, p_Cz, p_u, dscale, p_res)
      end select

      ! Release Ksep
      call storage_free(h_Ksep)


    case DEFAULT
      call output_line('Unsupported matrix format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlock')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the working routines follow
    
    !**************************************************************
    ! Assemble residual vector in 1D
    ! All matrices are stored in matrix format 7
    
    subroutine doResidualMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
                                 Cx, u, dscale, res)

      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NVAR

      real(DP), dimension(NEQ,NVAR), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: u_i,u_j,F_ij,F_ji
      integer :: ij,ji,i,j
      
      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute the fluxes
          call fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
        end do
      end do
    end subroutine doResidualMat7_1D


    !**************************************************************
    ! Assemble residual vector in 2D
    ! All matrices are stored in matrix format 7

    subroutine doResidualMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
                                 Cx, Cy, u, dscale, res)
      
      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NVAR

      real(DP), dimension(NEQ,NVAR), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: u_i,u_j,F_ij,F_ji
      integer :: ij,ji,i,j
   

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute the fluxes
          call fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
        end do
      end do
    end subroutine doResidualMat7_2D

    
    !**************************************************************
    ! Assemble residual vector in 3D
    ! All matrices are stored in matrix format 7

    subroutine doResidualMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
                                 Cx, Cy, Cz, u, dscale, res)
      
      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NVAR

      real(DP), dimension(NEQ,NVAR), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: u_i,u_j,F_ij,F_ji
      integer :: ij,ji,i,j
   

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
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
          call fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)

          ! Assemble residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
        end do
      end do
    end subroutine doResidualMat7_3D

    
    !**************************************************************
    ! Assemble residual vector in 1D
    ! All matrices are stored in matrix format 9

    subroutine doResidualMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
                                 Cx, u, dscale, res)
      
      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR

      real(DP), dimension(NEQ,NVAR), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: u_i,u_j,F_ij,F_ji
      integer :: ij,ji,i,j

   
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute the fluxes
          call fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
        end do
      end do
    end subroutine doResidualMat9_1D


    !**************************************************************
    ! Assemble residual vector in 2D
    ! All matrices are stored in matrix format 9

    subroutine doResidualMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
                                 Cx, Cy, u, dscale, res)
      
      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR

      real(DP), dimension(NEQ,NVAR), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: u_i,u_j,F_ij,F_ji
      integer :: ij,ji,i,j

   
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute the fluxes
          call fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
        end do
      end do
    end subroutine doResidualMat9_2D


    !**************************************************************
    ! Assemble residual vector in 3D
    ! All matrices are stored in matrix format 9

    subroutine doResidualMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
                                 Cx, Cy, Cz, u, dscale, res)
      
      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR

      real(DP), dimension(NEQ,NVAR), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: u_i,u_j,F_ij,F_ji
      integer :: ij,ji,i,j
   

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
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
          call fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)

          ! Assemble residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
        end do
      end do
    end subroutine doResidualMat9_3D
  end subroutine gfsys_buildResBlock
  
  ! *****************************************************************************

!<subroutine>
  
  subroutine gfsys_buildResScalar(RcoeffMatrices, ru, fcb_calcFlux, dscale, bclear, rres)

!<description>
    ! This subroutine assembles the residual vector. Note that the vectors are
    ! required as scalar vectors which are stored in the interleave format.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(IN) :: RcoeffMatrices

    ! solution vector
    type(t_vectorScalar), intent(IN) :: ru

    ! scaling factor
    real(DP), intent(IN) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(IN) :: bclear

    ! callback functions to compute local matrices
    include 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! residual vector
    type(t_vectorScalar), intent(INOUT) :: rres
!</inputoutput>
!</subroutine>
  
    ! local variables
    real(DP), dimension(:), pointer :: p_Cx,p_Cy,p_Cz,p_u,p_res
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    integer :: h_Ksep,ndim

    
    ! Clear vector?
    if (bclear) call lsyssc_clearVector(rres)
    
    ! Set pointers
    call lsyssc_getbase_double(ru, p_u)
    call lsyssc_getbase_double(rres, p_res)

    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)

    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_Cz)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalar')
      call sys_halt()
    end select

    ! What kind of matrix are we?
    select case(RcoeffMatrices(1)%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------

      ! Set pointers
      call lsyssc_getbase_Kld (RcoeffMatrices(1), p_Kld)
      call lsyssc_getbase_Kcol (RcoeffMatrices(1), p_Kcol)

      ! Create diagonal Ksep=Kld
      h_Ksep = ST_NOHANDLE
      call storage_copy(RcoeffMatrices(1)%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, RcoeffMatrices(1)%NEQ+1)
      
      ! How many dimensions do we have?
      select case(ndim)
      case (NDIM1D)
        call doResidualMat7_1D(p_Kld, p_Kcol, p_Ksep, RcoeffMatrices(1)%NEQ,&
                               ru%NVAR, p_Cx, p_u, dscale, p_res)
      case (NDIM2D)
        call doResidualMat7_2D(p_Kld, p_Kcol, p_Ksep, RcoeffMatrices(1)%NEQ,&
                               ru%NVAR, p_Cx, p_Cy, p_u, dscale, p_res)
      case (NDIM3D)
        call doResidualMat7_3D(p_Kld, p_Kcol, p_Ksep, RcoeffMatrices(1)%NEQ,&
                               ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, dscale, p_res)
      end select

      ! Release Ksep
      call storage_free(h_Ksep)

      
    case (LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 9
      !-------------------------------------------------------------------------

      ! Set pointer
      call lsyssc_getbase_Kld (RcoeffMatrices(1), p_Kld)
      call lsyssc_getbase_Kcol (RcoeffMatrices(1), p_Kcol)
      call lsyssc_getbase_Kdiagonal(RcoeffMatrices(1), p_Kdiagonal)

      ! Create diagonal Ksep=Kld
      h_Ksep = ST_NOHANDLE
      call storage_copy(RcoeffMatrices(1)%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, RcoeffMatrices(1)%NEQ+1)

      ! How many dimensions do we have?
      select case(ndim)
      case (NDIM1D)
        call doResidualMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                               RcoeffMatrices(1)%NEQ, ru%NVAR,&
                               p_Cx, p_u, dscale, p_res)
      case (NDIM2D)
        call doResidualMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                               RcoeffMatrices(1)%NEQ, ru%NVAR,&
                               p_Cx, p_Cy, p_u, dscale, p_res)
      case (NDIM3D)
        call doResidualMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                               RcoeffMatrices(1)%NEQ, ru%NVAR,&
                               p_Cx, p_Cy, p_Cz, p_u, dscale, p_res)
      end select
      
      ! Release Ksep
      call storage_free(h_Ksep)
      
      
    case DEFAULT
      call output_line('Unsupported matrix format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalar')
      call sys_halt()
    end select
    
  contains

    ! Here, the working routines follow
    
    !**************************************************************
    ! Assemble residual vector in 1D
    ! All matrices are stored in matrix format 7

    subroutine doResidualMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR,&
                                 Cx, u, dscale, res)

      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx
      real(DP), intent(IN) :: dscale  
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij,F_ji
      integer :: ij,ji,i,j
      
      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Compute the fluxes
          call fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
        end do
      end do
    end subroutine doResidualMat7_1D


    !**************************************************************
    ! Assemble residual vector in 2D
    ! All matrices are stored in matrix format 7

    subroutine doResidualMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR,&
                                 Cx, Cy, u, dscale, res)
      
      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy
      real(DP), intent(IN) :: dscale  
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij,F_ji
      integer :: ij,ji,i,j
   

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute the fluxes
          call fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
        end do
      end do
    end subroutine doResidualMat7_2D

    
    !**************************************************************
    ! Assemble residual vector in 3D
    ! All matrices are stored in matrix format 7

    subroutine doResidualMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR,&
                                 Cx, Cy, Cz, u, dscale, res)
      
      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz
      real(DP), intent(IN) :: dscale  
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij,F_ji
      integer :: ij,ji,i,j
   

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)
          
          ! Compute the fluxes
          call fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
        end do
      end do
    end subroutine doResidualMat7_3D

    
    !**************************************************************
    ! Assemble residual vector in 1D
    ! All matrices are stored in matrix format 9

    subroutine doResidualMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
                                 Cx, u, dscale, res)
      
      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx
      real(DP), intent(IN) :: dscale  
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij,F_ji
      integer :: ij,ji,i,j

   
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute the fluxes
          call fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)

          ! Assemble residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
        end do
      end do
    end subroutine doResidualMat9_1D


    !**************************************************************
    ! Assemble residual vector in 2D
    ! All matrices are stored in matrix format 9

    subroutine doResidualMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
                                 Cx, Cy, u, dscale, res)
      
      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy
      real(DP), intent(IN) :: dscale  
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij,F_ji
      integer :: ij,ji,i,j

   
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Compute the fluxes
          call fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)

          ! Assemble residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
        end do
      end do
    end subroutine doResidualMat9_2D


    !**************************************************************
    ! Assemble residual vector in 3D
    ! All matrices are stored in matrix format 9

    subroutine doResidualMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
                                 Cx, Cy, Cz, u, dscale, res)
      
      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz
      real(DP), intent(IN) :: dscale  
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij,F_ji
      integer :: ij,ji,i,j
   

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute averaged coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Compute the fluxes
          call fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
        end do
      end do
    end subroutine doResidualMat9_3D
  end subroutine gfsys_buildResScalar

  ! *****************************************************************************
  
!<subroutine>
  
  subroutine gfsys_buildResBlockTVD(RcoeffMatrices, ru, fcb_calcFlux,&
                                    fcb_calcCharacteristics, rafcstab,&
                                    dscale, bclear, rres)

!<description>
    ! This subroutine assembles the residual vector for FEM-TVD schemes.
    ! If the vectors contain only one block, then the scalar counterpart
    ! of this routine is called with the scalar subvectors.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(IN) :: RcoeffMatrices

    ! solution vector
    type(t_vectorBlock), intent(IN) :: ru

    ! scaling factor
    real(DP), intent(IN) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(IN) :: bclear

    ! callback functions to compute local matrices
    include 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab

    ! residual vector
    type(t_vectorBlock), intent(INOUT) :: rres
!</inputoutput>
!</subroutine>
    
    ! local variables
    real(DP), dimension(:), pointer :: p_Cx,p_Cy,p_Cz,p_u,p_res
    real(DP), dimension(:), pointer :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    integer :: h_Ksep,ndim

    
    ! Check if block vectors contain only one block.
    if ((ru%nblocks .eq. 1) .and. (rres%nblocks .eq. 1) ) then
      call gfsys_buildResScalarTVD(RcoeffMatrices, ru%RvectorBlock(1),&
          fcb_calcFlux, fcb_calcCharacteristics, rafcstab, dscale,&
          bclear, rres%RvectorBlock(1))
      return
    end if
    
    ! Clear vector?
    if (bclear) call lsysbl_clearVector(rres)
    
    ! Check if stabilisation is prepeared
    if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlockTVD')
      call sys_halt()
    end if
    
    ! Set pointers
    call lsysbl_getbase_double(ru, p_u)
    call lsysbl_getbase_double(rres, p_res)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)

    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)

    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_Cz)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlockTVD')
      call sys_halt()
    end select
    
    ! What kind of matrix are we?
    select case(RcoeffMatrices(1)%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld (RcoeffMatrices(1), p_Kld)
      call lsyssc_getbase_Kcol (RcoeffMatrices(1), p_Kcol)
      
      ! Create diagonal Ksep=Kld
      h_Ksep = ST_NOHANDLE
      call storage_copy(RcoeffMatrices(1)%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, RcoeffMatrices(1)%NEQ+1)

      ! How many dimensions do we have?
      select case(ndim)
      case (NDIM1D)
        call doLimitTVDMat7_1D(p_Kld, p_Kcol, p_Ksep, RcoeffMatrices(1)%NEQ,&
                               ru%nblocks, p_Cx, p_u, dscale,&
                               p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
      case (NDIM2D)
        call doLimitTVDMat7_2D(p_Kld, p_Kcol, p_Ksep, RcoeffMatrices(1)%NEQ,&
                               ru%nblocks, p_Cx, p_Cy, p_u, dscale,&
                               p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
      case (NDIM3D)
        call doLimitTVDMat7_3D(p_Kld, p_Kcol, p_Ksep, RcoeffMatrices(1)%NEQ,&
                               ru%nblocks, p_Cx, p_Cy, p_Cz, p_u, dscale,&
                               p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
      end select

      ! Release Ksep
      call storage_free(h_Ksep)
      
      
    case (LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 9
      !-------------------------------------------------------------------------

      ! Set pointer
      call lsyssc_getbase_Kld (RcoeffMatrices(1), p_Kld)
      call lsyssc_getbase_Kcol (RcoeffMatrices(1), p_Kcol)
      call lsyssc_getbase_Kdiagonal(RcoeffMatrices(1), p_Kdiagonal)

      ! Create diagonal Ksep=Kld
      h_Ksep = ST_NOHANDLE
      call storage_copy(RcoeffMatrices(1)%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, RcoeffMatrices(1)%NEQ+1)
        
      ! How many dimensions do we have?
      select case(ndim)
      case (NDIM1D)
        call doLimitTVDMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                               RcoeffMatrices(1)%NEQ, ru%nblocks,&
                               p_Cx, p_u, dscale,&
                               p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
      case (NDIM2D)
        call doLimitTVDMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                               RcoeffMatrices(1)%NEQ, ru%nblocks,&
                               p_Cx, p_Cy, p_u, dscale,&
                               p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
      case (NDIM3D)
        call doLimitTVDMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                               RcoeffMatrices(1)%NEQ, ru%nblocks,&
                               p_Cx, p_Cy, p_Cz, p_u, dscale,&
                               p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
      end select
      
      ! Release Ksep
      call storage_free(h_Ksep)

      
    case DEFAULT
      call output_line('Unsupported matrix format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResBlockTVD')
      call sys_halt()
    end select
    
    ! Set specifiers for Ps, Qs and Rs
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_LIMITER)

  contains

    ! Here, the working routines follow
    
    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 1D
    ! All matrices are stored in matrix format 7

    subroutine doLimitTVDMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, u,&
                                 dscale, pp, pm, qp, qm, rp, rm, res)

      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NVAR
      
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: pp,pm,qp,qm,rp,rm
      real(DP), dimension(NEQ,NVAR), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,R_ij,L_ij
      real(DP), dimension(NVAR) :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij,u_i,u_j
      integer :: ij,ji,i,j,iloc,jloc,ivar

      
      ! Clear P's and Q's for X-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0
      
      ! Loop over all rows (forward)
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute the fluxes
          call fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)
           
          ! Assemble high-order residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u_i, u_j, XDir1D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR

            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do
          
        end do
      end do

      ! Compute nodal correction factors for X-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)
      
      ! Loop over all rows (backward)
      do i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u_i, u_j, XDir1D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij  =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij
        end do
      end do
    end subroutine doLimitTVDMat7_1D


    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 2D
    ! All matrices are stored in matrix format 7

    subroutine doLimitTVDMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, Cy, u,&
                                 dscale, pp, pm, qp, qm, rp, rm, res)

      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NVAR
      
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: pp,pm,qp,qm,rp,rm
      real(DP), dimension(NEQ,NVAR), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,R_ij,L_ij
      real(DP), dimension(NVAR) :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij,u_i,u_j
      integer :: ij,ji,i,j,iloc,jloc,ivar

      
      ! Clear P's and Q's for X-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0
      
      ! Loop over all rows (forward)
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute the fluxes
          call fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)
                   
          ! Assemble high-order residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u_i, u_j, XDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR

            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do
          
        end do
      end do
      
      ! Compute nodal correction factors for X-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)
      
      ! Clear P's and Q's for Y-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0

      ! Loop over all rows (backward)
      do i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u_i, u_j, XDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij   =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij

          ! Compute characteristic fluxes in Y-direction
          call fcb_calcCharacteristics(u_i, u_j, YDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ji)-Cy(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR

            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do
          
        end do
      end do

      ! Compute nodal correction factors for Y-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)
      
      ! Loop over all rows (forward)
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute characteristic fluxes in Y-direction
          call fcb_calcCharacteristics(u_i, u_j, YDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ji)-Cy(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij  =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar =1, NVAR

            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij
        end do
      end do
    end subroutine doLimitTVDMat7_2D

    
    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 3D
    ! All matrices are stored in matrix format 7

    subroutine doLimitTVDMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, Cy, Cz, u,&
                                 dscale, pp, pm, qp, qm, rp, rm, res)

      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NVAR
      
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: pp,pm,qp,qm,rp,rm
      real(DP), dimension(NEQ,NVAR), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,R_ij,L_ij
      real(DP), dimension(NVAR) :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij,u_i,u_j
      integer :: ij,ji,i,j,iloc,jloc,ivar

      
      ! Clear P's and Q's for X-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0

      ! Loop over all rows (forward)
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
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
          call fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble high-order residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u_i, u_j, XDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR

            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do
          
        end do
      end do

      ! Compute nodal correction factors for X-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's for Y-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0

      ! Loop over all rows (backward)
      do i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u_i, u_j, XDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij   =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij

          ! Compute characteristic fluxes in Y-direction
          call fcb_calcCharacteristics(u_i, u_j, YDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ji)-Cy(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR

            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do
          
        end do
      end do

      ! Compute nodal correction factors for Y-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's for Z-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0

      ! Loop over all rows (forward)
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute characteristic fluxes in Y-direction
          call fcb_calcCharacteristics(u_i, u_j, YDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ji)-Cy(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij  =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar =1, NVAR

            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij

          ! Compute characteristic fluxes in Z-direction
          call fcb_calcCharacteristics(u_i, u_j, ZDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cz(ji)-Cz(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cz(ij)+Cz(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij

          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR
            
            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do

        end do
      end do

      ! Compute nodal correction factors for Z-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)

      ! Loop over all rows (backward)
      do i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute characteristic fluxes in Z-direction
          call fcb_calcCharacteristics(u_i, u_j, ZDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cz(ji)-Cz(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cz(ji)-Cz(ij))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij  =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar =1, NVAR

            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij
        end do
      end do
    end subroutine doLimitTVDMat7_3D
    
    
    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 1D
    ! All matrices are stored in matrix format 9

    subroutine doLimitTVDMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
                                 Cx, u, dscale, pp, pm, qp, qm, rp, rm, res)

      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR
      
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: pp,pm,qp,qm,rp,rm
      real(DP), dimension(NEQ,NVAR), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,R_ij,L_ij
      real(DP), dimension(NVAR) :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij,u_i,u_j
      integer :: ij,ji,i,j,iloc,jloc,ivar
      

      ! Clear P's and Q's for X-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0
      
      ! Loop over all rows (forward)
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute the fluxes
          call fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble high-order residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u_i, u_j, XDir1D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR

            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do
          
        end do
      end do

      ! Compute nodal correction factors for X-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)
      
      ! Loop over all rows (backward)
      do i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); Ksep(j) = Ksep(j)-1; ji = Ksep(j)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u_i, u_j, XDir1D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij   =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij
        end do
      end do
    end subroutine doLimitTVDMat9_1D


    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 2D
    ! All matrices are stored in matrix format 9

    subroutine doLimitTVDMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
                                 Cx, Cy, u, dscale, pp, pm, qp, qm, rp, rm, res)

      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR
      
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: pp,pm,qp,qm,rp,rm
      real(DP), dimension(NEQ,NVAR), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,R_ij,L_ij
      real(DP), dimension(NVAR) :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij,u_i,u_j
      integer :: ij,ji,i,j,iloc,jloc,ivar

      ! Clear P's and Q's for X-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0
      
      ! Loop over all rows (forward)
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute the fluxes
          call fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble high-order residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u_i, u_j, XDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR

            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do
          
        end do
      end do
      
      ! Compute nodal correction factors for X-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)
      
      ! Clear P's and Q's for Y-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0

      ! Loop over all rows (backward)
      do i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); Ksep(j) = Ksep(j)-1; ji = Ksep(j)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u_i, u_j, XDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij   =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij

          ! Compute characteristic fluxes in Y-direction
          call fcb_calcCharacteristics(u_i, u_j, YDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ji)-Cy(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR

            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do
          
        end do
      end do

      ! Compute nodal correction factors for Y-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)
      
      ! Loop over all rows (forward)
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
          
          ! Compute characteristic fluxes in Y-direction
          call fcb_calcCharacteristics(u_i, u_j, YDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ji)-Cy(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij  =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar =1, NVAR

            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij
        end do
      end do
    end subroutine doLimitTVDMat9_2D

    
    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 3D
    ! All matrices are stored in matrix format 9

    subroutine doLimitTVDMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
                                 Cx, Cy, Cz, u, dscale, pp, pm, qp, qm, rp, rm, res)

      real(DP), dimension(NEQ,NVAR), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR
      
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: pp,pm,qp,qm,rp,rm
      real(DP), dimension(NEQ,NVAR), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,R_ij,L_ij
      real(DP), dimension(NVAR) :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij,u_i,u_j
      integer :: ij,ji,i,j,iloc,jloc,ivar
      
      
      ! Clear P's and Q's for X-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0
      
      ! Loop over all rows (forward)
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
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
          call fcb_calcFlux(u_i, u_j, C_ij, C_ji, dscale, F_ij, F_ji)
                 
          ! Assemble high-order residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u_i, u_j, XDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR

            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do
          
        end do
      end do

      ! Compute nodal correction factors for X-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's for Y-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0

      ! Loop over all rows (backward)
      do i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); Ksep(j) = Ksep(j)-1; ji = Ksep(j)
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)
                    
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u_i, u_j, XDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij   =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij

          ! Compute characteristic fluxes in Y-direction
          call fcb_calcCharacteristics(u_i, u_j, YDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ji)-Cy(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR

            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do
          
        end do
      end do

      ! Compute nodal correction factors for Y-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's for Z-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0

      ! Loop over all rows (forward)
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute characteristic fluxes in Y-direction
          call fcb_calcCharacteristics(u_i, u_j, YDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ji)-Cy(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij  =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar =1, NVAR

            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij

          ! Compute characteristic fluxes in Z-direction
          call fcb_calcCharacteristics(u_i, u_j, ZDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cz(ji)-Cz(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cz(ij)+Cz(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij

          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR
            
            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do

        end do
      end do

      ! Compute nodal correction factors for Z-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)

      ! Loop over all rows (backward)
      do i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); Ksep(j) = Ksep(j)-1; ji = Ksep(j)

          ! Get solution values at nodes
          u_i = u(i,:); u_j = u(j,:)

          ! Compute characteristic fluxes in Z-direction
          call fcb_calcCharacteristics(u_i, u_j, ZDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cz(ji)-Cz(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cz(ij)+Cz(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij  =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar =1, NVAR

            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(i,:) = res(i,:)+F_ij
          res(j,:) = res(j,:)-F_ij
        end do
      end do
    end subroutine doLimitTVDMat9_3D
  end subroutine gfsys_buildResBlockTVD

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildResScalarTVD(RcoeffMatrices, ru, fcb_calcFlux,&
                                     fcb_calcCharacteristics, rafcstab,&
                                     dscale, bclear, rres)
    
!<description>
    ! This subroutine assembles the residual vector for FEM-TVD schemes
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(IN) :: RcoeffMatrices

    ! solution vector
    type(t_vectorScalar), intent(IN) :: ru

    ! scaling factor
    real(DP), intent(IN) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(IN) :: bclear

    ! callback functions to compute local matrices
    include 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab
    
    ! residual vector
    type(t_vectorScalar), intent(INOUT) :: rres
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Cx,p_Cy,p_Cz,p_u,p_res
    real(DP), dimension(:), pointer :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    integer :: h_Ksep,ndim

    ! Clear vector?
    if (bclear) call lsyssc_clearVector(rres)
    
    ! Check if stabilisation is prepeared
    if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarTVD')
      call sys_halt()
    end if

    ! Set pointers
    call lsyssc_getbase_double(ru, p_u)
    call lsyssc_getbase_double(rres, p_res)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)

    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)

    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_Cz)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarTVD')
      call sys_halt()
    end select
      
    ! What kind of matrix are we?
    select case(RcoeffMatrices(1)%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld (RcoeffMatrices(1), p_Kld)
      call lsyssc_getbase_Kcol (RcoeffMatrices(1), p_Kcol)
      
      ! Create diagonal Ksep=Kld
      h_Ksep = ST_NOHANDLE
      call storage_copy(RcoeffMatrices(1)%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, RcoeffMatrices(1)%NEQ+1)
      
      ! How many dimensions do we have?
      select case(ndim)
      case (NDIM1D)
        call doLimitTVDMat7_1D(p_Kld, p_Kcol, p_Ksep, RcoeffMatrices(1)%NEQ,&
                               ru%NVAR, p_Cx, p_u, dscale,&
                               p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
      case (NDIM2D)
        call doLimitTVDMat7_2D(p_Kld, p_Kcol, p_Ksep, RcoeffMatrices(1)%NEQ,&
                               ru%NVAR, p_Cx, p_Cy, p_u, dscale,&
                               p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
      case (NDIM3D)
        call doLimitTVDMat7_3D(p_Kld, p_Kcol, p_Ksep, RcoeffMatrices(1)%NEQ,&
                               ru%NVAR, p_Cx, p_Cy, p_Cz, p_u, dscale,&
                               p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
      end select

      ! Release Ksep
      call storage_free(h_Ksep)

      
    case (LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 9
      !-------------------------------------------------------------------------

      ! Set pointer
      call lsyssc_getbase_Kld (RcoeffMatrices(1), p_Kld)
      call lsyssc_getbase_Kcol (RcoeffMatrices(1), p_Kcol)
      call lsyssc_getbase_Kdiagonal(RcoeffMatrices(1), p_Kdiagonal)

      ! Create diagonal Ksep=Kld
      h_Ksep = ST_NOHANDLE
      call storage_copy(RcoeffMatrices(1)%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, RcoeffMatrices(1)%NEQ+1)
      
      ! How many dimensions do we have?
      select case(ndim)
      case (NDIM1D)
        call doLimitTVDMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                               RcoeffMatrices(1)%NEQ, ru%NVAR,&
                               p_Cx, p_u, dscale,&
                               p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
      case (NDIM2D)
        call doLimitTVDMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                               RcoeffMatrices(1)%NEQ, ru%NVAR,&
                               p_Cx, p_Cy, p_u, dscale,&
                               p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
      case (NDIM3D)
        call doLimitTVDMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                               RcoeffMatrices(1)%NEQ, ru%NVAR,&
                               p_Cx, p_Cy, p_Cz, p_u, dscale,&
                               p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_res)
      end select
      
      ! Release Ksep
      call storage_free(h_Ksep)
    
  
    case DEFAULT
      call output_line('Unsupported matrix format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildResScalarTVD')
      call sys_halt()
    end select
    
    ! Set specifiers for Ps, Qs and Rs
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_LIMITER)
    
  contains

    ! Here, the working routines follow
    
    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 1D
    ! All matrices are stored in matrix format 7

    subroutine doLimitTVDMat7_1D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, u,&
                                 dscale, pp, pm, qp, qm, rp, rm, res)

      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NVAR
      
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: pp,pm,qp,qm,rp,rm
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,R_ij,L_ij
      real(DP), dimension(NVAR) :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij,u_i,u_j
      integer :: ij,ji,i,j,iloc,jloc,ivar
      
      
      ! Clear P's and Q's for X-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0
      
      ! Loop over all rows (forward)
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute the fluxes
          call fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble high-order residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), XDir1D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR

            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do
          
        end do
      end do

      ! Compute nodal correction factors for X-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)
      
      ! Loop over all rows (backward)
      do i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1
          
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), XDir1D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij   =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij
        end do
      end do
    end subroutine doLimitTVDMat7_1D


    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 2D
    ! All matrices are stored in matrix format 7

    subroutine doLimitTVDMat7_2D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, Cy, u,&
                                 dscale, pp, pm, qp, qm, rp, rm, res)

      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NVAR
      
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: pp,pm,qp,qm,rp,rm
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,R_ij,L_ij
      real(DP), dimension(NVAR) :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij,u_i,u_j
      integer :: ij,ji,i,j,iloc,jloc,ivar
      

      ! Clear P's and Q's
      pp = 0;  pm = 0;  qp = 0;  qm = 0

      ! Loop over all rows (forward)
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Compute the fluxes
          call fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble high-order residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), XDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR

            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do
          
        end do
      end do
      
      ! Compute nodal correction factors for X-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)
      
      ! Clear P's and Q's for Y-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0

      ! Loop over all rows (backward)
      do i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1
                    
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), XDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij   =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij

          ! Compute characteristic fluxes in Y-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), YDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ji)-Cy(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR

            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do
          
        end do
      end do

      ! Compute nodal correction factors for Y-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)
      
      ! Loop over all rows (forward)
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute characteristic fluxes in Y-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), YDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ji)-Cy(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij  =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar =1, NVAR

            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij
        end do
      end do
    end subroutine doLimitTVDMat7_2D


    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 3D
    ! All matrices are stored in matrix format 7

    subroutine doLimitTVDMat7_3D(Kld, Kcol, Ksep, NEQ, NVAR, Cx, Cy, Cz, u,&
                                 dscale, pp, pm, qp, qm, rp, rm, res)

      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ,NVAR
      
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: pp,pm,qp,qm,rp,rm
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,R_ij,L_ij
      real(DP), dimension(NVAR) :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij,u_i,u_j
      integer :: ij,ji,i,j,iloc,jloc,ivar
      

      ! Clear P's and Q's for X-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0
      
      ! Loop over all rows (forward)
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
              
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Compute the fluxes
          call fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
                    
          ! Assemble high-order residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), XDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR

            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do
          
        end do
      end do

      ! Compute nodal correction factors for X-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's for Y-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0

      ! Loop over all rows (backward)
      do i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1
          
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), XDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij   =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij

          ! Compute characteristic fluxes in Y-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), YDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ji)-Cy(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR

            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do
          
        end do
      end do

      ! Compute nodal correction factors for Y-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's for Z-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0

      ! Loop over all rows (forward)
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute characteristic fluxes in Y-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), YDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ji)-Cy(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij  =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar =1, NVAR

            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij

          ! Compute characteristic fluxes in Z-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), ZDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cz(ji)-Cz(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cz(ij)+Cz(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij

          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR
            
            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do

        end do
      end do

      ! Compute nodal correction factors for Z-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)

      ! Loop over all rows (backward)
      do i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1

          ! Compute characteristic fluxes in Z-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), ZDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cz(ji)-Cz(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cz(ij)+Cz(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij  =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar =1, NVAR

            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij
        end do
      end do
    end subroutine doLimitTVDMat7_3D


    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 1D
    ! All matrices are stored in matrix format 9

    subroutine doLimitTVDMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
                                 Cx, u, dscale, pp, pm, qp, qm, rp, rm, res)

      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR
      
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: pp,pm,qp,qm,rp,rm
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,R_ij,L_ij
      real(DP), dimension(NVAR) :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij,u_i,u_j
      integer :: ij,ji,i,j,iloc,jloc,ivar
      

      ! Clear P's and Q's for X-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0
      
      ! Loop over all rows (forward)
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
              
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Compute the fluxes
          call fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble high-order residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), XDir1D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR

            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do
          
        end do
      end do

      ! Compute nodal correction factors for X-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)
      
      ! Loop over all rows (backward)
      do i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); Ksep(j) = Ksep(j)-1; ji = Ksep(j)
          
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), XDir1D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij   =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij
        end do
      end do
    end subroutine doLimitTVDMat9_1D

    
    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 2D
    ! All matrices are stored in matrix format 9

    subroutine doLimitTVDMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
                                 Cx, Cy, u, dscale, pp, pm, qp, qm, rp, rm, res)

      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR
      
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: pp,pm,qp,qm,rp,rm
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,R_ij,L_ij
      real(DP), dimension(NVAR) :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij,u_i,u_j
      integer :: ij,ji,i,j,iloc,jloc,ivar


      ! Clear P's and Q's for X-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0      
      
      ! Loop over all rows (forward)
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Compute the fluxes
          call fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble high-order residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), XDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR

            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do
          
        end do
      end do
      
      ! Compute nodal correction factors for X-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's for Y-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0

      ! Loop over all rows (backward)
      do i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); Ksep(j) = Ksep(j)-1; ji = Ksep(j)
          
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), XDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij   =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar = 1, NVAR

            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij

          ! Compute characteristic fluxes in Y-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), YDir2D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ji)-Cy(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR

            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do
          
        end do
      end do

      ! Compute nodal correction factors for Y-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)

      ! Loop over all rows (forward)
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute characteristic fluxes in Y-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), YDir2D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ji)-Cy(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij  =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar =1, NVAR

            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij
        end do
      end do
    end subroutine doLimitTVDMat9_2D


    !**************************************************************
    ! Assemble residual for low-order operator plus
    ! algebraic flux correction of TVD-type in 3D
    ! All matrices are stored in matrix format 9

    subroutine doLimitTVDMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
                                 Cx, Cy, Cz, u, dscale, pp, pm, qp, qm, rp, rm, res)

      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR
      
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: pp,pm,qp,qm,rp,rm
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: res
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: A_ij,R_ij,L_ij
      real(DP), dimension(NVAR) :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij,u_i,u_j
      integer :: ij,ji,i,j,iloc,jloc,ivar
      

      ! Clear P's and Q's for X-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0
      
      ! Loop over all rows (forward)
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Compute the fluxes
          call fcb_calcFlux(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Assemble high-order residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)+F_ji
                   
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), XDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR

            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do
          
        end do
      end do

      ! Compute nodal correction factors for X-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's for Y-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0

      ! Loop over all rows (backward)
      do i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); Ksep(j) = Ksep(j)-1; ji = Ksep(j)
          
          ! Compute characteristic fluxes in X-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), XDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij  =  0.5_DP*(Cx(ji)-Cx(ij))*Lbd_ij
          ks_ij  =  0.5_DP*(Cx(ij)+Cx(ji))*Lbd_ij
          F_ij   = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij   =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar = 1, NVAR
            
            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij

          ! Compute characteristic fluxes in Y-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), YDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ji)-Cy(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          
          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR

            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do
          
        end do
      end do

      ! Compute nodal correction factors for Y-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's for Z-direction
      pp = 0;  pm = 0;  qp = 0;  qm = 0

      ! Loop over all rows (forward)
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute characteristic fluxes in Y-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), YDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cy(ji)-Cy(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cy(ij)+Cy(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij  =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar =1, NVAR

            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij

          ! Compute characteristic fluxes in Z-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), ZDir3D, W_ij, Lbd_ij)

          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cz(ji)-Cz(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cz(ij)+Cz(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij

          ! Update sums of downstream/upstream edge contributions
          do ivar = 1, NVAR
            
            ! Set node orientation
            if (ka_ij(ivar) > 0) then
              iloc = j; jloc = i
              F_ij(ivar) = -F_ij(ivar)
            else
              iloc = i; jloc = j
            end if
            
            ! Assemble P's and Q's
            if (F_ij(ivar) > 0) then
              pp(ivar,iloc) = pp(ivar,iloc)+F_ij(ivar)
              qm(ivar,iloc) = qm(ivar,iloc)-F_ij(ivar)
              qp(ivar,jloc) = qp(ivar,jloc)+F_ij(ivar)
            else
              pm(ivar,iloc) = pm(ivar,iloc)+F_ij(ivar)
              qp(ivar,iloc) = qp(ivar,iloc)-F_ij(ivar)
              qm(ivar,jloc) = qm(ivar,jloc)+F_ij(ivar)
            end if
          end do

        end do
      end do

      ! Compute nodal correction factors for Z-direction
      rp = afcstab_limit( pp, qp, 1.0_DP, 1.0_DP)
      rm = afcstab_limit(-pm,-qm, 1.0_DP, 1.0_DP)

      ! Loop over all rows (backward)
      do i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); Ksep(j) = Ksep(j)-1; ji = Ksep(j)

          ! Compute characteristic fluxes in Z-direction
          call fcb_calcCharacteristics(u(:,i), u(:,j), ZDir3D, W_ij, Lbd_ij, R_ij)
          
          ! Compute antidiffusive fluxes
          ka_ij =  0.5_DP*(Cz(ji)-Cz(ij))*Lbd_ij
          ks_ij =  0.5_DP*(Cz(ij)+Cz(ji))*Lbd_ij
          F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
          W_ij  =  abs(ka_ij)*W_ij
          
          ! Construct characteristic fluxes
          do ivar =1, NVAR

            ! Limit characteristic fluxes
            if (ka_ij(ivar) < 0) then
              if (F_ij(ivar) > 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,i)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,i)*F_ij(ivar)
              end if
            else
              if (F_ij(ivar) < 0) then
                W_ij(ivar) = W_ij(ivar)+rp(ivar,j)*F_ij(ivar)
              else
                W_ij(ivar) = W_ij(ivar)+rm(ivar,j)*F_ij(ivar)
              end if
            end if
          end do
          
          ! Transform back into conservative variables
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, W_ij, 1, 0.0_DP, F_ij, 1)
          
          ! Assemble high-resolution residual vector
          res(:,i) = res(:,i)+F_ij
          res(:,j) = res(:,j)-F_ij
        end do
      end do
    end subroutine doLimitTVDMat9_3D
  end subroutine gfsys_buildResScalarTVD

  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************

!<subroutine>

  subroutine gfsys_getbase_double(rmatrix, rarray, bisFullMatrix)

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
    type(t_matrixBlock), intent(IN) :: rmatrix
!</input>

!<output>
    ! The array
    type(t_array), dimension(:,:), intent(OUT) :: rarray

    ! OPTIONAL: indicator for full block matrix
    logical, intent(OUT), optional :: bisFullMatrix
!</output>
!</subroutine>
    
    ! local variables
    integer :: iblock,jblock
    logical :: bisFull
    
    ! Check if array is compatible
    if (rmatrix%nblocksPerCol .ne. size(rarray,1) .or.&
        rmatrix%nblocksPerRow .ne. size(rarray,2)) then
      call output_line('Block matrix and array are not compatible!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_getbase_double')
      call sys_halt()
    end if

    ! Assign pointers
    bisFull = .true.
    do iblock = 1, rmatrix%nblocksPerCol
      do jblock = 1, rmatrix%nblocksPerRow
        if (lsyssc_isExplicitMatrix1D(rmatrix%RmatrixBlock(iblock,jblock))) then
          call lsyssc_getbase_double(rmatrix%RmatrixBlock(iblock,jblock),&
                                     rarray(iblock,jblock)%Da)
        else
          nullify(rarray(iblock,jblock)%Da)
          bisFull = .false.
        end if
      end do
    end do

    if (present(bisFullMatrix)) bisFullMatrix = bisFull
  end subroutine gfsys_getbase_double

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_getbase_single(rmatrix, rarray, bisFullMatrix)

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
    type(t_matrixBlock), intent(IN) :: rmatrix
!</input>

!<output>
    ! The array
    type(t_array), dimension(:,:), intent(OUT) :: rarray

    ! OPTIONAL: indicator for full block matrix
    logical, intent(OUT), optional :: bisFullMatrix
!</output>
!</subroutine>
    
    ! local variables
    integer :: iblock,jblock
    logical :: bisFull
    
    ! Check if array is compatible
    if (rmatrix%nblocksPerCol .ne. size(rarray,1) .or.&
        rmatrix%nblocksPerRow .ne. size(rarray,2)) then
      call output_line('Block matrix and array are not compatible!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_getbase_single')
      call sys_halt()
    end if

    ! Assign pointers
    bisFull = .true.
    do iblock = 1, rmatrix%nblocksPerCol
      do jblock = 1, rmatrix%nblocksPerRow
        if (lsyssc_isExplicitMatrix1D(rmatrix%RmatrixBlock(iblock,jblock))) then
          call lsyssc_getbase_single(rmatrix%RmatrixBlock(iblock,jblock),&
                                     rarray(iblock,jblock)%Fa)
        else
          nullify(rarray(iblock,jblock)%Fa)
          bisFull = .false.
        end if
      end do
    end do

    if (present(bisFullMatrix)) bisFullMatrix = bisFull
  end subroutine gfsys_getbase_single

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_getbase_int(rmatrix, rarray, bisFullMatrix)

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
    type(t_matrixBlock), intent(IN) :: rmatrix
!</input>

!<output>
    ! The array
    type(t_array), dimension(:,:), intent(OUT) :: rarray

    ! OPTIONAL: indicator for full block matrix
    logical, intent(OUT), optional :: bisFullMatrix
!</output>
!</subroutine>
    
    ! local variables
    integer :: iblock,jblock
    logical :: bisFull
    
    ! Check if array is compatible
    if (rmatrix%nblocksPerCol .ne. size(rarray,1) .or.&
        rmatrix%nblocksPerRow .ne. size(rarray,2)) then
      call output_line('Block matrix and array are not compatible!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_getbase_int')
      call sys_halt()
    end if

    ! Assign pointers
    bisFull = .true.
    do iblock = 1, rmatrix%nblocksPerCol
      do jblock = 1, rmatrix%nblocksPerRow
        if (lsyssc_isExplicitMatrix1D(rmatrix%RmatrixBlock(iblock,jblock))) then
          call lsyssc_getbase_int(rmatrix%RmatrixBlock(iblock,jblock),&
                                  rarray(iblock,jblock)%Ia)
        else
          nullify(rarray(iblock,jblock)%Da)
          bisFull = .false.
        end if
      end do
    end do

    if (present(bisFullMatrix)) bisFullMatrix = bisFull
  end subroutine gfsys_getbase_int

end module groupfemsystem
