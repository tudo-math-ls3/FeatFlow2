!##############################################################################
!# ****************************************************************************
!# <name> blockmatassemblybase </name>
!# ****************************************************************************
!#
!# <purpose>
!# Contains the structures needed in the blockmatassembly routines.
!# </purpose>
!##############################################################################

module blockmatassemblybase

  use fsystem
  use derivatives
  use elementbase
  use cubature
  use transformation
  use boundary
  use triangulation
  use linearsystemscalar
  use perfconfig

  use feevaluation2

  implicit none

!<constants>

!<constantblock description="Constants for the initialisation">

  ! Do not compute extra information
  integer(I32), parameter :: BMA_CALC_NONE = 0

  ! Compute real world coordinates. Usually needed in callback routines
  integer(I32), parameter :: BMA_CALC_REALCOORDS = 2**0

  ! Default flags
  integer(I32), parameter :: BMA_CALC_STANDARD = BMA_CALC_REALCOORDS

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! Collects the information for the assembly of one
  ! matrix block in the block matrix.
  type t_bmaMatrixData

    !<!--
    ! ##############################################
    ! Data for the user. Needed during the assembly.
    ! ##############################################
    !-->

    ! Whether or not there is data in this block at all.
    ! If this is .false., all pointers in this structure are undefined!
    logical :: bhasData = .false.

    ! Defines whether or not the matrix data for this block is shared
    ! with another block in the block matrix. If this is the case,
    ! the data for this block is automatically calculated during the
    ! calculation of the other block and does not have to be assembled.
    logical :: bsharedMatrixData = .false.

    ! Arrays for the basis function values in the cubature points
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest => null()
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial => null()

    ! Arrays saving the DOF`s in the elements
    integer, dimension(:,:), pointer :: p_IdofsTest => null()
    integer, dimension(:,:), pointer :: p_IdofsTrial => null()

    ! Specifies whether the corresponding matrix is interleaved.
    ! =.FALSE.: Matrix is not interleaved. The local matrix entries
    !           must be saved to p_Dentry. p_DentryIntl is undefined.
    ! =.TRUE.:  Matrix is interleaved. The local matrix entries
    !           must be saved to p_DentryIntl. p_Dentry is undefined.
    logical :: bisInterleaved = .false.

    ! Arrays saving the entries of the local matrices.
    !
    !       p_DentryIntl(ndofTrial,ndofTest,nelements).
    !
    ! NOTE: Only for non-interleaved matrices; there is p_Dentry=>null for
    ! interleaved vectors. In this case, p_DentryIntl defines the matrix data.
    !
    ! NOTE: If the matrix data is shared with another block, p_Dentry=>null.
    !
    ! WARNING: All local matrices are saved transposed, i.e. we have
    !     p_Dentry(column,row,element) = p_Dentry(ndofTrial,ndofTest,nelements).
    ! This has to be taken care of during the assembly and was 
    ! introduced to gain higher speed as it avoids extensive 
    ! jumps in the memory access.
    real(DP), dimension(:,:,:), pointer :: p_Dentry => null()

    ! Arrays saving the entries of the local matrices.
    !
    !       p_DentryIntl(nvar,ndofTrial,ndofTest,nelements).
    !
    ! NOTE: Only for interleaved matrices; there is p_DentryIntl=>null for
    ! non-interleaved vectors. In this case, p_Dentry defines the matrix data.
    !
    ! NOTE: If the matrix data is shared with another block, p_DentryIntl=>null.
    !
    ! WARNING: All local matrices are saved transposed, i.e. we have
    !     p_DentryIntl(var,column,row,element) = p_Dentry(nvar,ndofTrial,ndofTest,nelements).
    ! This has to be taken care of during the assembly and was 
    ! introduced to gain higher speed as it avoids extensive 
    ! jumps in the memory access.
    real(DP), dimension(:,:,:,:), pointer :: p_DentryIntl => null()

    ! Number of local DOF`s in the trial/test space
    integer :: ndofTrial = 0
    integer :: ndofTest = 0

    ! Number of variables per matrix entry. Only for interleaved
    ! matrices. =1 for non-interleaved matrices.
    integer :: nvar = 1

    !<!--
    ! ##############################################
    ! Internal data from the assembly routines.
    ! ##############################################
    !-->

    ! Defines whether or not the matrix structure for this block is shared
    ! with another block in the block matrix. If this is the case,
    ! the data for this block (indices of matrix entries) is
    ! only calculated once and reused.
    logical :: bsharedMatrixStructure = .false.

    ! Type of element to evaluate in the trial and test space.
    integer(I32) :: celementTrial = 0
    integer(I32) :: celementTest = 0

    ! Whether trial and test space is identical
    logical :: bIdenticalTrialAndTest = .false.

    ! Arrays saving the indices of the local matrices
    ! => NULL() if the matrix data is shared with another block.
    integer, dimension(:,:,:), pointer :: p_Kentry => null()

    ! Reference to the matrix to be computed.
    type(t_matrixScalar), pointer :: p_rmatrix => null()

    ! Pointer to the matrix data.
    real(DP), dimension(:), pointer :: p_Da => null()

    ! Index of the corresponding FEM data structure for trial and test space
    integer :: iidxFemDataTrial = 0
    integer :: iidxFemDataTest = 0

  end type

!</typeblock>

  public :: t_bmaMatrixData

!<typeblock>

  ! Encapsules all dynamic data which is calculated during the
  ! matrix assembly. Roughly speaking, this is the data
  ! necessary for the assembly of a block of elements.
  ! Every block of elements may be assembled in parallel
  ! to other blocks.
  type t_bmaMatrixAssemblyData

    ! Element evaluation set encapsuling coordinates on of the cubature
    ! points on the reference element, real elements, Jacobian determinants etc.
    type(t_evalElementSet) :: revalElementSet

    ! The number of elements in revalElementSet of all submatrices
    ! whose cubature points on the reference element(s) have already 
    ! been initialised.
    integer :: ninitialisedElements = 0

    ! Number of elements in p_IelementList
    integer :: nelements = 0

    ! Pointer to the list of elements currently in progress
    integer, dimension(:), pointer :: p_IelementList => null()

    ! For every element and every cubature point in the element set,
    ! cubature weight to be used for calculating the matric entries.
    ! Is calculated from the actual cubature weight and the
    ! corresponding Jacobian determinant of the cubature point.
    !    p_DcubWeight(npoints,nelements)
    real(DP), dimension(:,:), pointer :: p_DcubWeight => null()

    ! Data of all involved FEM spaces.
    type(t_fev2FemData), dimension(:), pointer :: p_RfemData => null()

    ! Array of the blocks for all submatrices.
    ! Saves the data of all submatrices and contains temporary memory
    ! used during the assembly.
    type(t_bmaMatrixData), dimension(:,:), pointer :: p_RmatrixData => null()

  end type

!</typeblock>

  public :: t_bmaMatrixAssemblyData

!<typeblock>

  ! A matrix assembly structure that saves crucial data during the matrix assembly
  ! with bma_buildMatrix.
  type t_bmaMatrixAssembly

    ! Currently active element distribution
    integer :: ielementDistr = 0

    ! Maximum number of elements to handle simultaneously.
    integer :: nelementsPerBlock = 0

    ! Type of cubature formula to use.
    ! This cubature formula is simultaneously used for all blocks
    ! in the matrix.
    integer(I32) :: ccubType = CUB_UNDEFINED

    ! Type of transformation
    integer(I32) :: ctrafoType = TRAFO_ID_UNKNOWN

    ! Cubature weights
    real(DP), dimension(:), pointer :: p_Domega => null()

    ! An array that takes coordinates of the cubature formula on the reference element
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef => null()

    ! Number of cubature points per element
    integer :: ncubp = 0

    ! Basic evaluation tag; use the same for all the element spaces
    integer(I32) :: cevaluationTag = 0

    ! Underlying triangulation
    type(t_triangulation), pointer :: p_rtriangulation => null()

    ! Underlying boundary definition
    type(t_boundary), pointer :: p_rboundary => null()

    ! Template structure for the dynamic data that changes
    ! during the assembly of each element block.
    ! During the assembly, this is used to generate the
    ! actual dynamic data arrays that hold the values of the DOFs,...
    type(t_bmaMatrixAssemblyData) :: rassemblyDataTemplate

    ! Template structure for the dynamic evaluation of vectors.
    ! During the assembly, the vectors in this structure are
    ! automatically evaluated in the cubature points.
    type(t_fev2Vectors) :: revalVectorsTemplate

    ! Performance configuration. 
    type(t_perfconfig), pointer :: p_rperfconfig => null()

  end type

!</typeblock>

  public :: t_bmaMatrixAssembly

!<typeblock>

  ! Collects the information for the assembly of one
  ! vector block in the block vector.
  type t_bmaVectorData

    !<!--
    ! ##############################################
    ! Data for the user. Needed during the assembly.
    ! ##############################################
    !-->

    ! Arrays for the basis function values in the cubature points
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest => null()

    ! Arrays saving the DOF`s in the elements.
    !   p_IdofsTest(ndofTest,nelements)
    integer, dimension(:,:), pointer :: p_IdofsTest => null()

    ! Specifies whether the corresponding vector is interleaved.
    ! =.FALSE.: Vector is not interleaved. The local vector entries
    !           must be saved to p_Dentry. p_DentryIntl is undefined.
    ! =.TRUE.:  Vector is interleaved. The local vector entries
    !           must be saved to p_DentryIntl. p_Dentry is undefined.
    logical :: bisInterleaved = .false.

    ! Arrays saving the entries of the local vectors.
    !
    !     p_Dentry(ndofTest,nelements).
    !
    ! NOTE: Only for non-Interleaved vectors, there is p_Dentry=>null for
    ! interleaved vectors. In this case, p_DentryIntl defines the vector data.
    real(DP), dimension(:,:), pointer :: p_Dentry => null()

    ! Arrays saving the entries of the local vectors.
    !
    !     p_Dentry(nvar,ndofTest,nelements)
    !
    ! NOTE: Only for interleaved vectors, there is p_DentryIntl=>null for
    ! non-interleaved vectors. In this case, p_Dentry defines the vector data.
    real(DP), dimension(:,:,:), pointer :: p_DentryIntl => null()

    ! Number of local DOF`s in the test space
    integer :: ndofTest = 0

    ! Number of variables per vector entry. Only for interleaved
    ! vectors. =1 for non-interleaved vectors.
    integer :: nvar = 1

    !<!--
    ! ##############################################
    ! Internal data from the assembly routines.
    ! ##############################################
    !-->

    ! Type of element to evaluate in the trial and test space.
    integer(I32) :: celementTest = 0

    ! Reference to the vector to be computed.
    type(t_vectorScalar), pointer :: p_rvector => null()

    ! Pointer to the vector data.
    real(DP), dimension(:), pointer :: p_Ddata => null()

    ! Index of the corresponding FEM data structure for trial and test space
    integer :: iidxFemDataTest = 0

  end type

!</typeblock>

  public :: t_bmaVectorData

!<typeblock>

  ! Encapsules all dynamic data which is calculated during the
  ! vector assembly. Roughly speaking, this is the data
  ! necessary for the assembly of a block of elements.
  ! Every block of elements may be assembled in parallel
  ! to other blocks.
  type t_bmaVectorAssemblyData

    ! Element evaluation set encapsuling coordinates on of the cubature
    ! points on the reference element, real elements, Jacobian determinants etc.
    type(t_evalElementSet) :: revalElementSet

    ! The number of elements in revalElementSet of all submatrices
    ! whose cubature points on the reference element(s) have already 
    ! been initialised.
    integer :: ninitialisedElements = 0

    ! Number of elements in p_IelementList
    integer :: nelements = 0

    ! Pointer to the list of elements currently in progress
    integer, dimension(:), pointer :: p_IelementList => null()

    ! For every element and every cubature point in the element set,
    ! cubature weight to be used for calculating the matric entries.
    ! Is calculated from the actual cubature weight and the
    ! corresponding Jacobian determinant of the cubature point.
    !    p_DcubWeight(npoints,nelements)
    real(DP), dimension(:,:), pointer :: p_DcubWeight => null()

    ! Data of all involved FEM spaces.
    type(t_fev2FemData), dimension(:), pointer :: p_RfemData => null()

    ! Array of the blocks for all subvectors.
    ! Saves the data of all subvectors and contains temporary memory
    ! used during the assembly.
    type(t_bmaVectorData), dimension(:), pointer :: p_RvectorData => null()

  end type

!</typeblock>

  public :: t_bmaVectorAssemblyData

!<typeblock>

  ! A vector assembly structure that saves crucial data during the matrix assembly
  ! with bma_buildVector.
  type t_bmaVectorAssembly

    ! Currently active element distribution
    integer :: ielementDistr = 0

    ! Maximum number of elements to handle simultaneously.
    integer :: nelementsPerBlock = 0

    ! Type of cubature formula to use.
    ! This cubature formula is simultaneously used for all blocks
    ! in the matrix.
    integer(I32) :: ccubType = CUB_UNDEFINED

    ! Type of transformation
    integer(I32) :: ctrafoType = TRAFO_ID_UNKNOWN

    ! Cubature weights
    real(DP), dimension(:), pointer :: p_Domega => null()

    ! An array that takes coordinates of the cubature formula on the reference element
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef => null()

    ! Number of cubature points per element
    integer :: ncubp = 0

    ! Basic evaluation tag; use the same for all the element spaces
    integer(I32) :: cevaluationTag = 0

    ! Underlying triangulation
    type(t_triangulation), pointer :: p_rtriangulation => null()

    ! Underlying boundary definition
    type(t_boundary), pointer :: p_rboundary => null()

    ! Template structure for the dynamic data that changes
    ! during the assembly of each element block.
    ! During the assembly, this is used to generate the
    ! actual dynamic data arrays that hold the values of the DOFs,...
    type(t_bmaVectorAssemblyData) :: rassemblyDataTemplate

    ! Template structure for the dynamic evaluation of vectors.
    ! During the assembly, the vectors in this structure are
    ! automatically evaluated in the cubature points.
    type(t_fev2Vectors) :: revalVectorsTemplate

    ! Performance configuration. 
    type(t_perfconfig), pointer :: p_rperfconfig => null()

  end type

!</typeblock>

  public :: t_bmaVectorAssembly

!</types>

end module
