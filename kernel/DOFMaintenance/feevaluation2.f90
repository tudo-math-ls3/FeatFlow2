!##############################################################################
!# ****************************************************************************
!# <name> feevaluation2 </name>
!# ****************************************************************************
!#
!# <purpose>
!# Realises extended routines for the evaluation of FEM basis functions
!# and FEM functions.
!#
!# The following routines can be found here:
!#
!# 1.) fev2_getBderSize
!#     -> Determines the number of entries in a Bder array up to a given
!#        derivative
!#
!# 2.) fev2_prepareFemDataBMat
!#     -> Basic initialisation of a FEM data structure based on a block matrix
!#
!# 3.) fev2_prepareFemDataBVec
!#     -> Basic initialisation of a FEM data structure based on a block vector
!#
!# 4.) fev2_prepareFemDataSMat
!#     -> Basic initialisation of a FEM data structure based on a scalar matrix
!#
!# 5.) fev2_prepareFemDataSVec
!#     -> Basic initialisation of a FEM data structure based on a scalar vector
!#
!# 6.) fev2_cleanupFemData
!#     -> Cleans up a FEM data structure initialised by fev2_prepareFemDataXXXX.
!#
!# 7.) fev2_createFemData
!#     -> Full initialisation of a FEM data structure; allocate memory
!#
!# 8.) fev2_releaseFemData
!#     -> Cleanup of a FEM data structure
!#
!# 9.) fev2_evaluateFemData
!#     -> Save the values of FEM basis functions in a FEM data structure
!#
!# 10.) fev2_addVectorToEvalList
!#      -> Add a vector to a list of vectors to be evaluated simultaneously
!#         in a set of points on a set of elements
!#
!# 11.) fev2_addDummyVectorToEvalList
!#      -> Adds a scalar dummy entry to the list of vectors to be evaluated.
!#         Can be used as temporary memory during the evaluation.
!#
!# 12.) fev2_releaseVectorList
!#      -> Release a vector evaluation list
!#
!# 13.) fev2_initVectorEval
!#      -> Initialise a vector evaluation list for the evaluation
!#         of vectors
!#
!# 14.) fev2_prepareVectorEval
!#      -> (Re-)allocates temporary memory for the evaluation.
!#         Must be called after fev2_initVectorEval.
!#
!# 15.) fev2_doneVectorEval
!#      -> Release a vector evaluation list
!#
!# 16.) fev2_evaluateVectors
!#      -> Evaluate all vectors in a vector evaluation list in a set of points
!#         on a set of elements
!#
!# 17.) fev2_prepareFemDataVecEval
!#      -> Basic initialisation of a FEM data structure based on a
!#         vector evaluation structure
!#
!# 18.) fev2_calcDofMapping
!#      -> Calculates the DOF mapping 
!#
!# 19.) fev2_copyFemData
!#      -> Creates a copy of a FEM data structure.
!# </purpose>
!##############################################################################

module feevaluation2

  use fsystem
  use storage
  use basicgeometry
  use collection, only: t_collection
  use derivatives
  use element
  use elementpreprocessing
  use genoutput
  use scalarpde
  use dofmapping
  use spatialdiscretisation
  use transformation
  use triangulation
  use perfconfig
  
  use linearsystemscalar
  use linearsystemblock

  implicit none
  
  private

!<types>

!<typeblock>

  ! Holds all data for the evaluation of FEM basis functions.
  ! Use fev2_prepareFemData to prepare the structure for a FEM
  ! space and fev2_createFemData to allocate memory.
  ! Afterwards, use fev2_evaluateFemData to evaluate the FEM basis
  ! functions and save their values to the p_Dbas variable.
  type t_fev2FemData
  
    ! Type of element.
    integer(I32) :: celement = 0
    
    ! Number of local DOF`s in the trial/test space
    integer :: ndof = 0
    
    ! Dimension of the FE space.
    integer :: ndimfe = 0
    
    ! Maximum derivative to be computed.
    ! =-1: Automatic, evaluate all up to the highest possible derivative.
    integer :: nmaxDerivative = -1
    
    ! Last index in Bder which is set to TRUE.
    integer :: nmaxDerivativeIdx = 0
    
    ! Type of transformation
    integer(I32) :: ctrafoType = TRAFO_ID_UNKNOWN
    
    ! Array to tell the element which derivatives to calculate
    logical, dimension(EL_MAXNDER) :: Bder = .false.

    ! Arrays for the basis function values in the cubature points.
    !   Dbas(ndof*ndimfe,ideriv,ipt,iel)
    ! For vector valued basis functions of dimension d, there is
    !   Dbas(1..ndof,ideriv,ipt,iel) = first dimension
    !   Dbas((1..ndof) + ndof,ideriv,ipt,iel) = 2nd dimension
    !   ...
    !   Dbas((1..ndof) + ndimfe*ndof,ideriv,ipt,iel) = d'th dimension
    real(DP), dimension(:,:,:,:), pointer :: p_Dbas => null()
    
    ! Arrays saving the DOF`s in the elements
    !   Idofs(ndof,nelements)
    integer, dimension(:,:), pointer :: p_Idofs => null()

    ! Reference to the discretisation
    type(t_spatialDiscretisation), pointer :: p_rdiscr => null()

  end type

!</typeblock>

  public :: t_fev2FemData

!<typeblock>

  ! Data of all involved FEM spaces.
  type t_fev2FemDataBlocks
  
    ! Numer of FEM spaces.
    integer :: ncount = 0
    
    ! List of FEM space data arrays. May point to NULL
    ! if ncount=0.
    type(t_fev2FemData), dimension(:), pointer :: p_RfemData => null()

  end type

!</typeblock>

  public :: t_fev2FemDataBlocks

!<typeblock>

  ! This structure saves the values of a FEM function in a set
  ! of points on a set of elements. Used, e.g., for nonlinearities.
  type t_fev2VectorData

    ! Specifies whether the corresponding vector is interleaved.
    ! =.FALSE.: Vector is not interleaved. The local vector entries
    !           can be found in p_Ddata. p_DdataIntl is undefined.
    ! =.TRUE.:  Vector is interleaved. The local vector entries
    !           can be found in p_DdataIntl. p_Ddata is undefined.
    logical :: bisInterleaved = .false.
  
    ! Pointer to an array with values in all cubature points.
    !   Ddata(npointsPerElement,nelements,nmaxDerivativeIdx).
    !
    ! The last index specifies the type of the derivative,
    ! so Ddata(:,:,DER_DERIV2D_X) specifies the first X-derivative.
    !
    ! NOTE: For (non-interleaved) dummy vectors, the memory is not initialised.
    !
    ! NOTE: Only for non-Interleaved vectors, there is p_Ddata=>null for
    ! interleaved vectors. In this case, p_DentryIntl defines the vector data.
    !
    ! NOTE: For vector valued FE spaces, there is p_Ddata=>null(),
    ! p_DdataIntl=>null() and p_DdataVec specifying the values.
    real(DP), dimension(:,:,:), pointer :: p_Ddata => null()

    ! Pointer to an array with values in all cubature points.
    !   p_DdataIntl(nvar,npointsPerElement,nelements,nmaxDerivativeIdx).
    !
    ! The last index specifies the type of the derivative,
    ! so p_DdataIntl(:,:,:,DER_DERIV2D_X) specifies the first X-derivative.
    !
    ! NOTE: For (interleaved) dummy vectors, the memory is not initialised.
    !
    ! NOTE: Only for interleaved vectors, there is p_DdataIntl=>null for
    ! non-interleaved vectors. In this case, p_Dentry defines the vector data.
    !
    ! NOTE: For vector valued FE spaces, there is p_Ddata=>null(),
    real(DP), dimension(:,:,:,:), pointer :: p_DdataIntl => null()
    
    ! Pointer to an array with values in all cubature points.
    !   DdataVec(ndimfe,npointsPerElement,nelements,nmaxDerivativeIdx).
    !
    ! The last index specifies the type of the derivative,
    ! so Ddata(:,:,DER_DERIV2D_X) specifies the first X-derivative.
    !
    ! This only applies for vector valued FE spaces of dimension ndimfe.
    real(DP), dimension(:,:,:,:), pointer :: p_DdataVec => null()
    
    ! Reference to the vector or NULL, if there is no vector
    ! associated. The latter case appears for `dummy` vectors.
    ! Dummy vectors provide additional temporary memory which is
    ! preallocated and which can be arbitrarily used by the callback
    ! routines.
    type(t_vectorScalar), pointer :: p_rvector => null()
    
    ! Number of variables per vector entry. Only for interleaved
    ! vectors. =1 for non-interleaved vectors.
    integer :: nvar = 1
    
    ! Maximum derivative to be computed.
    ! If this is =-1, nmaxDerivativeIdx>0 and p_rvector=>null(), this vector
    ! is a dummy vector.
    integer :: nmaxDerivative = 0

    ! Last index of the corresponding Bder array which is set to TRUE.
    ! Specifies the last dimension in p_Ddata.
    integer :: nmaxDerivativeIdx = 0
    
    ! Dimension of the underlying FEM space (for vector valued FEM spaces.
    ! =1 for standard FEM spaces)
    integer :: ndimfe = 1
    
    ! Index of the corresponding FEM data structure for the discretisation
    ! in a corresponding arrays of FEM data structures.
    integer :: iidxFemData = 0

  end type

!</typeblock>

  public :: t_fev2VectorData

!<typeblock>

  ! Encapsules a list of vectors which is evaluated in all cubature
  ! points and passed to the callback routines.
  ! Vectors can be added to this list by fev2_addVectorToEvalList.
  ! With fev2_initVectorEval, memory for the evaluation can be allocated.
  ! A call to fev2_evaluateVectors evaluates all vectors in this list
  ! in a set of points on a set of elements (e.g., in all cubature points).
  ! Afterwards, the values of the FEM functions in the points
  ! can be obtained via p_RvectorData(.)%p_Ddata.
  type t_fev2Vectors
  
    ! Number of vectors in the list.
    integer :: ncount = 0
    
    ! Number of points per element to evaluate simultaneously
    integer :: npointsPerElement = 0
    
    ! Number of elements to evaluate simultaneously
    integer :: nelements = 0
    
    ! List of vectors.
    type(t_fev2VectorData), dimension(:), pointer :: p_RvectorData => null()

  end type

!</typeblock>

  public :: t_fev2Vectors

!</types>  

  public :: fev2_getBderSize
  public :: fev2_prepareFemDataBMat
  public :: fev2_prepareFemDataBVec
  public :: fev2_prepareFemDataSMat
  public :: fev2_prepareFemDataSVec
  public :: fev2_cleanupFemData
  public :: fev2_createFemData
  public :: fev2_calcDofMapping
  public :: fev2_releaseFemData
  public :: fev2_evaluateFemData
  public :: fev2_copyFemData
  
  public :: fev2_addVectorToEvalList
  public :: fev2_addDummyVectorToEvalList
  public :: fev2_releaseVectorList
  public :: fev2_initVectorEval
  public :: fev2_prepareVectorEval
  public :: fev2_doneVectorEval
  public :: fev2_evaluateVectors
  public :: fev2_prepareFemDataVecEval

contains

  !****************************************************************************

!<subroutine>

  subroutine fev2_getBderSize(ndim,imaxDerivative,imaxDerivativeIdx)

!<description>
  ! Calculates the last index in a BDER array which must be set to
  ! TRUE to evaluate up to the imaxDerivative'th derivative.
!</description>

!<input>
  ! Dimension of the space. NDIM1D, NDIM2D or NDIM3D
  integer, intent(in) :: ndim

  ! Maximum derivative.
  integer, intent(in) :: imaxDerivative
!</input>

!<output>
  ! Returns the maximum index in Bder which must be set = true.
  integer, intent(out) :: imaxDerivativeIdx
!</output>

!</subroutine>

    select case (ndim)
    case (NDIM1D)
      select case (imaxDerivative)
      case (0)
        imaxDerivativeIdx = DER_FUNC1D
      case (1)
        imaxDerivativeIdx = DER_DERIV1D_X
      case (2:)
        imaxDerivativeIdx = DER_DERIV1D_XX
      end select
    
    case (NDIM2D)
      select case (imaxDerivative)
      case (0)
        imaxDerivativeIdx = DER_FUNC2D
      case (1)
        imaxDerivativeIdx = DER_DERIV2D_Y
      case (2:)
        imaxDerivativeIdx = DER_DERIV2D_YY
      end select
      
    case (NDIM3D)
      select case (imaxDerivative)
      case (0)
        imaxDerivativeIdx = DER_FUNC3D
      case (1)
        imaxDerivativeIdx = DER_DERIV3D_Z
      case (2:)
        imaxDerivativeIdx = DER_DERIV3D_ZZ
      end select
    end select

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fev2_initBder(celement,imaxDerivative,imaxDerivativeIdx,Bder)

!<description>
  ! Initialises a BDER array (derivative evaluation flag) up to the
  ! imaxDerivative'th derivative
!</description>

!<input>
  ! Underlying finite element (for automatic detection
  integer(I32), intent(in) :: celement

  ! Maximum derivative.
  ! =-1: maximum available, detect automatically.
  integer, intent(in) :: imaxDerivative
!</input>

!<output>
  ! Returns the maximum index in Bder which is = true.
  integer, intent(out) :: imaxDerivativeIdx

  ! Derivative flags to initialise
  logical, dimension(:), intent(out) :: Bder
!</output>

!</subroutine>

    ! If it is =-1, take the default, i.e., the maximum available
    if (imaxDerivative .le. -1) then
      imaxDerivativeIdx = elem_getMaxDerivative(celement)
    else
      call fev2_getBderSize(elem_igetDimension(celement),&
          imaxDerivative,imaxDerivativeIdx)
    end if
    
    ! Fill BDER with TRUE up to element k.
    Bder(:) = .false.
    Bder(1:imaxDerivativeIdx) = .true.

  end subroutine

  !****************************************************************************

  subroutine addDiscr (rfemDataBlocks,rentry)
  ! Adds rentry to a list represented by an array
  type(t_fev2FemDataBlocks), intent(inout) :: rfemDataBlocks
  type(t_spatialDiscretisation), intent(in), target :: rentry
  
    type(t_fev2FemData), dimension(:), pointer :: p_RfemData
  
    ! Is data available?
    if (.not. associated(rfemDataBlocks%p_RfemData)) then
      allocate(rfemDataBlocks%p_RfemData(16))
    end if
    
    ! Array large enough?
    if (rfemDataBlocks%ncount .ge. ubound(rfemDataBlocks%p_RfemData,1)) then
      ! Reallocate
      allocate(p_RfemData(rfemDataBlocks%ncount+16))
      p_RfemData(1:rfemDataBlocks%ncount) = rfemDataBlocks%p_RfemData(:)
      deallocate (rfemDataBlocks%p_RfemData)
      rfemDataBlocks%p_RfemData => p_RfemData
    end if
  
    ! Append
    rfemDataBlocks%ncount = rfemDataBlocks%ncount + 1
    rfemDataBlocks%p_RfemData(rfemDataBlocks%ncount)%p_rdiscr => rentry
  end subroutine

  !****************************************************************************
  
  integer function containsDiscr (rfemDataBlocks,rentry)
  ! returns <> 0 (the index) if the list Rlist contains rentry
  type(t_fev2FemDataBlocks), intent(in) :: rfemDataBlocks
  type(t_spatialDiscretisation), intent(in), target :: rentry
  
    integer :: i
    
    do i=1,rfemDataBlocks%ncount
      if (associated(rfemDataBlocks%p_RfemData(i)%p_rdiscr,rentry)) then
        containsDiscr = i
        return
      end if
    end do
    
    containsDiscr = 0
    
  end function
    
  !****************************************************************************

  subroutine releaseDiscr(rfemDataBlocks)
  
  ! Cleans up a FEM evaluation structure which was prepared with
  ! fev2_prepareFemDataXXXX.
  type(t_fev2FemDataBlocks), intent(inout) :: rfemDataBlocks

    ! Deallocate, finish.
    if (associated(rfemDataBlocks%p_RfemData)) deallocate(rfemDataBlocks%p_RfemData)
    rfemDataBlocks%ncount = 0

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fev2_prepareFemDataBMat(rmatrix,rfemDataBlocks,ielementDistr, &
      RmaxDerivativeTest,RmaxDerivativeTrial)

!<description>
  ! Initialise a FEM structure based on all trial/test spaces
  ! appearing in a block matrix.
  ! No memory is allocated.
!</description>

!<input>
  ! The matrix which is going to be assembled.
  type(t_matrixBlock), intent(in) :: rmatrix
  
  ! Element distribution to use.
  integer, intent(in) :: ielementDistr
  
  ! OPTIONAL: For every block in the matrix, maximum
  ! derivative of the basis functions to be computed. If not
  ! specified or an entry is =-1, the maximum available derivative for 
  ! each FEM space is the default.
  integer, dimension(:,:), intent(in), optional :: RmaxDerivativeTest
  integer, dimension(:,:), intent(in), optional :: RmaxDerivativeTrial
!</input>

!<output>
  ! Structure for all involved FEM spaces.
  ! If the structure is empty, new memory is allocated.
  ! Otherwise, new FEM structures are appended.
  type(t_fev2FemDataBlocks), intent(inout) :: rfemDataBlocks
!</output>

!</subroutine>
  
    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rdiscrTrial,p_rdiscrTest
    integer :: i,j,imaxTest,imaxTrial
    
    ! Basic initialisation of the structure.
    ! Loop over all blocks. Figure out which blocks have data.
    ! Get the data arrays for these blocks.
    do j=1,rmatrix%nblocksPerRow
    
      do i=1,rmatrix%nblocksPerCol
        
        if (lsysbl_isSubmatrixPresent (rmatrix,i,j)) then
        
          ! Get the information about the trial and test spaces
          ! of this block.
          p_rdiscrTrial => rmatrix%RmatrixBlock(i,j)%p_rspatialDiscrTrial
          p_rdiscrTest  => rmatrix%RmatrixBlock(i,j)%p_rspatialDiscrTest
          
          ! Some basic checks...
          if (p_rdiscrTrial%inumFESpaces .ne. p_rdiscrTest%inumFESpaces) then
            call output_line ("Discretisation structures incompatible.",&
                OU_CLASS_ERROR,OU_MODE_STD,"fev2_prepareFemDataBMat")
            call sys_halt()
          end if
          
          ! Maximum derivative
          imaxTest = -1
          if (present(RmaxDerivativeTest)) imaxTest = RmaxDerivativeTest(i,j)

          imaxTrial = imaxTest
          if (present(RmaxDerivativeTrial)) imaxTrial = RmaxDerivativeTrial(i,j)

          ! Append the space and initialise if we do not have it.
          call fev2_prepareFemDataSMat(rmatrix%RmatrixBlock(i,j),rfemDataBlocks,&
              ielementDistr,imaxTest,imaxTrial)

        end if
      
      end do
    
    end do
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fev2_prepareFemDataBVec(rvector,rfemDataBlocks,ielementDistr, &
      RmaxDerivative)

!<description>
  ! Initialise a FEM structure based on all test spaces
  ! appearing in a block vector.
  ! No memory is allocated.
!</description>

!<input>
  ! The matrix which is going to be assembled.
  type(t_vectorBlock), intent(in) :: rvector
  
  ! Element distribution to use.
  integer, intent(in) :: ielementDistr
  
  ! OPTIONAL: For every block in the vector, maximum
  ! derivative of the basis functions to be computed. If not
  ! specified or an entry is =-1, the maximum available derivative for 
  ! each FEM space is the default.
  integer, dimension(:), intent(in), optional :: RmaxDerivative
!</input>

!<output>
  ! Structure for all involved FEM spaces.
  ! If the structure is empty, new memory is allocated.
  ! Otherwise, new FEM structures are appended.
  type(t_fev2FemDataBlocks), intent(inout) :: rfemDataBlocks
!</output>

!</subroutine>
  
    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rdiscrTest1,p_rdiscrTest
    integer :: i,imaxTest

    ! Basic initialisation of the structure.
    ! Loop over all blocks. Figure out which blocks have data.
    ! Get the data arrays for these blocks.
    do i=1,rvector%nblocks
      
      ! Get the information about the trial and test spaces
      ! of this block.
      p_rdiscrTest  => rvector%RvectorBlock(i)%p_rspatialDiscr
      if (i .eq. 1) then
        ! Discretisation of the first block
        p_rdiscrTest1  => rvector%RvectorBlock(i)%p_rspatialDiscr
      end if
      
      ! Some basic checks...
      if (p_rdiscrTest1%inumFESpaces .ne. p_rdiscrTest%inumFESpaces) then
        call output_line ("Discretisation structures incompatible.",&
            OU_CLASS_ERROR,OU_MODE_STD,"fev2_prepareFemDataBVec")
        call sys_halt()
      end if
      
      ! Maximum derivative
      imaxTest = -1
      if (present(RmaxDerivative)) imaxTest = RmaxDerivative(i)

      ! Append the FEM space and initialise if we do not have it.
      call fev2_prepareFemDataSVec(rvector%RvectorBlock(i),rfemDataBlocks, &
          ielementDistr,imaxTest)
      
    end do
    
  end subroutine

!****************************************************************************

!<subroutine>

  recursive subroutine fev2_prepareFemDataSVec(rvector,rfemDataBlocks,ielementDistr, &
      imaxDerivative)

!<description>
  ! Initialise a FEM structure based on the FEM spaces
  ! appearing in a scalar vector.
  ! No memory is allocated.
!</description>

!<input>
  ! The matrix which is going to be assembled.
  type(t_vectorScalar), intent(in) :: rvector
  
  ! Element distribution to use.
  integer, intent(in) :: ielementDistr
  
  ! OPTIONAL: Maximum derivative of the basis functions to be computed. If not
  ! specified or an entry is =-1, the maximum available derivative for 
  ! each FEM space is the default.
  integer, optional :: imaxDerivative
!</input>

!<output>
  ! Structure for all involved FEM spaces.
  ! If the structure is empty, new memory is allocated.
  ! Otherwise, new FEM structures are appended.
  type(t_fev2FemDataBlocks), intent(inout) :: rfemDataBlocks
!</output>

!</subroutine>
  
    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rdiscrTest
    integer :: ivecDiscrTest
    type(t_fev2FemData), dimension(:), pointer :: p_RfemData
    
    ! Get the information about the trial and test spaces
    ! of this block.
    p_rdiscrTest  => rvector%p_rspatialDiscr
    
    if (.not. associated(p_rdiscrTest)) return
    
    ! Do we have the test space already?
    ivecDiscrTest = containsDiscr (rfemDataBlocks,p_rdiscrTest)
    
    if (ivecDiscrTest .eq. 0) then
      
      ! Add the test space
      call addDiscr (rfemDataBlocks,p_rdiscrTest)

      ! Pointer to the data blocks
      p_RfemData => rfemDataBlocks%p_RfemData
    
      ! Initialise the content
      p_RfemData(rfemDataBlocks%ncount)%celement = &
          p_RfemData(rfemDataBlocks%ncount)%p_rdiscr%RelementDistr(ielementDistr)%celement
      p_RfemData(rfemDataBlocks%ncount)%ctrafoType = &
          p_RfemData(rfemDataBlocks%ncount)%p_rdiscr%RelementDistr(ielementDistr)%ctrafoType
      p_RfemData(rfemDataBlocks%ncount)%ndof = &
          elem_igetNDofLoc(p_RfemData(rfemDataBlocks%ncount)%celement)
      p_RfemData(rfemDataBlocks%ncount)%ndimfe = &
          elem_igetFeDimension(p_RfemData(rfemDataBlocks%ncount)%celement)

      ivecDiscrTest = rfemDataBlocks%ncount
    else
      ! Pointer to the data blocks
      p_RfemData => rfemDataBlocks%p_RfemData
    end if

    ! Maximum derivative for the test space.
    ! Exception: If nmaxDerivative = -1, the maximum derivative is taken
    ! anyway, so we must not taje the maximum.
    if (p_RfemData(ivecDiscrTest)%nmaxDerivative .ne. -1) then
      p_RfemData(ivecDiscrTest)%nmaxDerivative = &
          max(p_RfemData(ivecDiscrTest)%nmaxDerivative,imaxDerivative)
    end if
        
    call fev2_initBder(&
        p_RfemData(ivecDiscrTest)%celement,p_RfemData(ivecDiscrTest)%nmaxDerivative,&
        p_RfemData(ivecDiscrTest)%nmaxDerivativeIdx,p_RfemData(ivecDiscrTest)%Bder)

  end subroutine

!****************************************************************************

!<subroutine>

  subroutine fev2_prepareFemDataSMat(rmatrix,rfemDataBlocks,ielementDistr, &
      imaxDerivativeTest,imaxDerivativeTrial)

!<description>
  ! Initialise a FEM structure based on the FEM spaces
  ! appearing in a scalar matrix.
  ! No memory is allocated.
!</description>

!<input>
  ! The matrix which is going to be assembled.
  type(t_matrixScalar), intent(in) :: rmatrix
  
  ! Element distribution to use.
  integer, intent(in) :: ielementDistr
  
  ! OPTIONAL: Maximum derivative of the basis functions to be computed. If not
  ! specified or an entry is =-1, the maximum available derivative for 
  ! each FEM space is the default. Test space
  integer, optional :: imaxDerivativeTest

  ! OPTIONAL: Maximum derivative of the basis functions to be computed. If not
  ! specified or an entry is =-1, the maximum available derivative for 
  ! each FEM space is the default. Trial space
  integer, optional :: imaxDerivativeTrial
!</input>

!<output>
  ! Structure for all involved FEM spaces.
  ! If the structure is empty, new memory is allocated.
  ! Otherwise, new FEM structures are appended.
  type(t_fev2FemDataBlocks), intent(inout) :: rfemDataBlocks
!</output>

!</subroutine>
  
    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rdiscrTest,p_rdiscrTrial
    integer :: ivecDiscrTest,ivecDiscrTrial
    integer :: imaxTest, imaxTrial
    type(t_fev2FemData), dimension(:), pointer :: p_RfemData

    ! Maximum derivative
    imaxTest = -1
    if (present(imaxDerivativeTest)) imaxTest = imaxDerivativeTest
    
    imaxTrial = imaxTest
    if (present(imaxDerivativeTrial)) imaxTrial = imaxDerivativeTrial

    ! Get the information about the trial and test spaces
    ! of this block.
    p_rdiscrTest  => rmatrix%p_rspatialDiscrTest
    p_rdiscrTrial  => rmatrix%p_rspatialDiscrTest
    
    ! Do we have the test space already?
    if (associated(p_rdiscrTest)) then

      ivecDiscrTest = containsDiscr (rfemDataBlocks,p_rdiscrTest)
      
      if (ivecDiscrTest .eq. 0) then
        
        ! Add the test space
        call addDiscr (rfemDataBlocks,p_rdiscrTest)
        
        ! Pointer to the data blocks
        p_RfemData => rfemDataBlocks%p_RfemData
      
        ! Initialise the content
        p_RfemData(rfemDataBlocks%ncount)%celement = &
            p_RfemData(rfemDataBlocks%ncount)%p_rdiscr%RelementDistr(ielementDistr)%celement
        p_RfemData(rfemDataBlocks%ncount)%ctrafoType = &
            p_RfemData(rfemDataBlocks%ncount)%p_rdiscr%RelementDistr(ielementDistr)%ctrafoType
        p_RfemData(rfemDataBlocks%ncount)%ndof = &
            elem_igetNDofLoc(p_RfemData(rfemDataBlocks%ncount)%celement)
        p_RfemData(rfemDataBlocks%ncount)%ndimfe = &
            elem_igetFeDimension(p_RfemData(rfemDataBlocks%ncount)%celement)

        ivecDiscrTest = rfemDataBlocks%ncount
      else
        ! Pointer to the data blocks
        p_RfemData => rfemDataBlocks%p_RfemData
      end if

      ! Maximum derivative for the test space
      ! Exception: If nmaxDerivative = -1, the maximum derivative is taken
      ! anyway, so we must not taje the maximum.
      if (p_RfemData(ivecDiscrTest)%nmaxDerivative .ne. -1) then
        p_RfemData(ivecDiscrTest)%nmaxDerivative = &
            max(p_RfemData(ivecDiscrTest)%nmaxDerivative,imaxTest)
      end if

      call fev2_initBder(&
          p_RfemData(ivecDiscrTest)%celement,p_RfemData(ivecDiscrTest)%nmaxDerivative,&
          p_RfemData(ivecDiscrTest)%nmaxDerivativeIdx,p_RfemData(ivecDiscrTest)%Bder)

    end if
    
    ! Do we have the trial space already?
    if (associated(p_rdiscrTrial)) then
    
      ivecDiscrTrial = containsDiscr (rfemDataBlocks,p_rdiscrTrial)

      if (ivecDiscrTrial .eq. 0) then
        
        ! Add the trial space
        call addDiscr (rfemDataBlocks,p_rdiscrTrial)
      
        ! Pointer to the data blocks
        p_RfemData => rfemDataBlocks%p_RfemData

        ! Initialise the content
        p_RfemData(rfemDataBlocks%ncount)%celement = &
            p_RfemData(rfemDataBlocks%ncount)%p_rdiscr%RelementDistr(ielementDistr)%celement
        p_RfemData(rfemDataBlocks%ncount)%ctrafoType = &
            p_RfemData(rfemDataBlocks%ncount)%p_rdiscr%RelementDistr(ielementDistr)%ctrafoType
        p_RfemData(rfemDataBlocks%ncount)%ndof = &
            elem_igetNDofLoc(p_RfemData(rfemDataBlocks%ncount)%celement)
        p_RfemData(rfemDataBlocks%ncount)%ndimfe = &
            elem_igetFeDimension(p_RfemData(rfemDataBlocks%ncount)%celement)

        ivecDiscrTrial = rfemDataBlocks%ncount
      else
        ! Pointer to the data blocks
        p_RfemData => rfemDataBlocks%p_RfemData
      end if

      ! Maximum derivative for the trial space
      ! Exception: If nmaxDerivative = -1, the maximum derivative is taken
      ! anyway, so we must not taje the maximum.
      if (p_RfemData(ivecDiscrTrial)%nmaxDerivative .ne. -1) then
        p_RfemData(ivecDiscrTrial)%nmaxDerivative = &
            max(p_RfemData(ivecDiscrTrial)%nmaxDerivative,imaxTrial)
      end if
          
      call fev2_initBder(&
          p_RfemData(ivecDiscrTrial)%celement,p_RfemData(ivecDiscrTrial)%nmaxDerivative,&
          p_RfemData(ivecDiscrTrial)%nmaxDerivativeIdx,p_RfemData(ivecDiscrTrial)%Bder)

    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fev2_createFemData(rfemDataBlocks,npointsperelement,nelements)

!<description>
  ! Initialises a set of FEM data structurse which was prepared with
  ! fev2_prepareFemDataXXXX. Allocates memory such that the evaluation can
  ! be applied.
!</description>

!<input>
  ! Number of points per element
  integer, intent(in) :: npointsperelement

  ! Number of elements to simultaneously process
  integer, intent(in) :: nelements
!</input>

!<inputoutput>
  ! Structure for all involved FEM spaces.
  type(t_fev2FemDataBlocks), intent(inout) :: rfemDataBlocks
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_fev2FemData), pointer :: p_rfemData
    integer :: i

    ! Loop through all FEM data blocks.
    do i=1,rfemDataBlocks%ncount
    
      p_rfemData => rfemDataBlocks%p_RfemData(i)

      ! Allocate memory for the values of the basis functions    
      allocate(p_rfemData%p_Dbas(p_rfemData%ndof*p_rfemData%ndimfe,&
          p_rfemData%nmaxDerivativeIdx,npointsperelement,nelements))
      
      ! Allocate memory for the degrees of freedoms on the
      ! elements to be processed.
      allocate(p_rfemData%p_Idofs(p_rfemData%ndof,nelements))
    
    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fev2_releaseFemData(rfemDataBlocks)

!<description>
  ! Cleans up a set of FEM evaluation structures.
!</description>

!<inputoutput>
  ! Structure for all involved FEM spaces.
  type(t_fev2FemDataBlocks), intent(inout) :: rfemDataBlocks
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_fev2FemData), pointer :: p_rfemData
    integer :: i

    ! Release all allocated memory.
    do i=1,rfemDataBlocks%ncount
    
      p_rfemData => rfemDataBlocks%p_RfemData(i)
      deallocate(p_rfemData%p_Idofs)
      deallocate(p_rfemData%p_Dbas)
    
    end do
    
    deallocate(rfemDataBlocks%p_RfemData)
    rfemDataBlocks%ncount = 0

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fev2_cleanupFemData(rfemDataBlocks)

!<description>
  ! Cleans up a FEM evaluation structure which was prepared with
  ! fev2_prepareFemDataXXXX.
!</description>

!<inputoutput>
  ! Structure for all involved FEM spaces.
  type(t_fev2FemDataBlocks), intent(inout) :: rfemDataBlocks
!</inputoutput>

!</subroutine>

    call releaseDiscr(rfemDataBlocks)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fev2_copyFemData(rfemDataBlocksDest,rfemDataBlocksSource)

!<description>
  ! Creaates a copy of a FEM data block structure.
!</description>

!<input>
  ! Source structure
  type(t_fev2FemDataBlocks), intent(in) :: rfemDataBlocksSource
!</input>

!<output>
  ! Destination structure
  type(t_fev2FemDataBlocks), intent(inout) :: rfemDataBlocksDest
!</output>

!</subroutine>
    
    integer :: icount

    rfemDataBlocksDest%ncount = rfemDataBlocksSource%ncount
    
    icount = size(rfemDataBlocksSource%p_RfemData)
    allocate(rfemDataBlocksDest%p_RfemData(icount))
    rfemDataBlocksDest%p_RfemData(:) = rfemDataBlocksSource%p_RfemData(:)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fev2_calcDofMapping(rfemDataBlocks,IelementList)

!<description>
  ! Calculates the DOF mapping, i.e., the global degrees of freedom
  ! for all elements in the element list Ielement list.
!</description>

!<input>
  ! Element list
  integer, dimension(:), intent(in) :: IelementList
!</input>

!<inputoutput>
  ! All involved FEM spaces.
  type(t_fev2FemDataBlocks), intent(inout) :: rfemDataBlocks
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_fev2FemData), pointer :: p_rfemData
    integer :: i

    ! Loop through all FEM data blocks.
    do i=1,rfemDataBlocks%ncount
    
      p_rfemData => rfemDataBlocks%p_RfemData(i)

      ! Calculates the degrees of freedom on the elements
      if (associated(p_rfemData%p_rdiscr)) then
        call dof_locGlobMapping_mult(p_rfemData%p_rdiscr, &
            IelementList, p_rfemData%p_Idofs)
      end if
    
    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fev2_evaluateFemData(rfemDataBlocks,revalElementSet)

!<description>
  ! Evaluates the FEM basis functions according to the evaluation
  ! set revalElementSet. The values of the basis functions are saved
  ! in the p_Dbas arrays of the FEM data structures in rfemDataBlocks.
!</description>

!<input>
  ! Element evaluation set encapsuling coordinates on of the cubature
  ! points on the reference element, real elements, Jacobian determinants etc.
  type(t_evalElementSet), intent(in) :: revalElementSet
!</input>

!<inputoutput>
  ! All involved FEM spaces.
  ! The p_Dbas arrays in all the structures are filled with the
  ! values of the FEM basis functions.
  ! The structure must have been initialised by fev2_createFemData.
  type(t_fev2FemDataBlocks), intent(inout) :: rfemDataBlocks
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
    type(t_fev2FemData), pointer :: p_rfemData

    ! Loop through all FEM spaces to calculate the FEM basis functions
    do i=1,rfemDataBlocks%ncount
    
      p_rfemData => rfemDataBlocks%p_RfemData(i)

      ! Evaluate the basis functions.
      call elem_generic_sim2 (p_rfemData%celement, &
          revalElementSet, p_rfemData%Bder,p_rfemData%p_Dbas)

    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fev2_addVectorToEvalList(revalVectors,rvector,nmaxDerivative)

!<description>
  ! Adds a scalar vector to the list of vectors to be evaluated.
!</description>

!<input>
  ! Vector to be added to the list
  type(t_vectorScalar), intent(in), target :: rvector
  
  ! Maximum derivative to be calculated.
  ! =0: calculate the function value only,
  ! =1: Calculate the complete first derivative (X, Y and Z direction)
  ! etc.
  integer, intent(in) :: nmaxDerivative
!</input>

!<inputoutput>
  ! List of vectors to be automatically evaluated
  type(t_fev2Vectors), intent(inout) :: revalVectors
!</inputoutput>

!</subroutine>
  
    ! local variables
    type(t_fev2VectorData), dimension(:), pointer :: p_RvectorData 

    ! Add the vector
    if (revalVectors%ncount .eq. 0) then
      allocate(revalVectors%p_RvectorData(16))
    else
      if (revalVectors%ncount .ge. size(revalVectors%p_RvectorData)) then
        ! Reallocate
        allocate(p_RvectorData(size(revalVectors%p_RvectorData)+16))
        p_RvectorData(1:revalVectors%ncount) = &
            revalVectors%p_RvectorData(1:revalVectors%ncount)
        deallocate(revalVectors%p_RvectorData)
        revalVectors%p_RvectorData => p_RvectorData
      end if
    end if
    
    ! Append
    revalVectors%ncount = revalVectors%ncount + 1
    revalVectors%p_RvectorData(revalVectors%ncount)%p_rvector => rvector
    
    ! Maximum derivative to be calculated
    revalVectors%p_RvectorData(revalVectors%ncount)%nmaxDerivative = nmaxDerivative
    
    ! Variables per vector entry (for interleaved vectors)
    if (rvector%nvar .ne. 1) then
      revalVectors%p_RvectorData(revalVectors%ncount)%bisInterleaved = .true.
      revalVectors%p_RvectorData(revalVectors%ncount)%nvar = rvector%nvar
    end if
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fev2_addDummyVectorToEvalList(revalVectors,nsubarrays,nvar)

!<description>
  ! Adds a scalar dummy entry to the list of vectors to be evaluated.
  ! During the evaluation, memory is allocated for this dummy entry, which
  ! allows the assembly routines to compute intermediate data.
  ! However, since no actual vector is associated, no FEM function is
  ! evaluated in the cubature points, thus the `evaluation` does not need
  ! computational time.
  ! Note that the allocated memory stays uninitialised until given
  ! free.
!</description>

!<input>
  ! OPTIONAL: Number of subarrays.
  ! If not specified, there is exactly one subarray allocated.
  ! If specified, there are nsubarray memory blocks allocated in memory
  ! and associated to this dummy vector.
  integer, intent(in), optional :: nsubarrays
  
  ! OPTINOAL: Number of variables per vector entry.
  ! If set to > 1, a memory block for interleaved access is created.
  integer, intent(in), optional :: nvar
!</input>

!<inputoutput>
  ! List of vectors to be automatically evaluated
  type(t_fev2Vectors), intent(inout) :: revalVectors
!</inputoutput>

!</subroutine>
  
    ! local variables
    type(t_fev2VectorData), dimension(:), pointer :: p_RvectorData 

    ! Add a dummy entry
    if (revalVectors%ncount .eq. 0) then
      allocate(revalVectors%p_RvectorData(16))
    else
      if (revalVectors%ncount .ge. size(revalVectors%p_RvectorData)) then
        ! Reallocate
        allocate(p_RvectorData(size(revalVectors%p_RvectorData)+16))
        p_RvectorData(1:revalVectors%ncount) = &
            revalVectors%p_RvectorData(1:revalVectors%ncount)
        deallocate(revalVectors%p_RvectorData)
        revalVectors%p_RvectorData => p_RvectorData
      end if
    end if
    
    ! Append
    revalVectors%ncount = revalVectors%ncount + 1
    
    ! Nullify the pointer to mark it as dummy.
    nullify(revalVectors%p_RvectorData(revalVectors%ncount)%p_rvector)
    
    ! Maximum derivative to be calculated. Set to -1 here.
    revalVectors%p_RvectorData(revalVectors%ncount)%nmaxDerivative = -1
    
    ! Number of subarrays is saved to nmaxDerivativeIndex
    revalVectors%p_RvectorData(revalVectors%ncount)%nmaxDerivativeIdx = 1
    
    if (present(nsubarrays)) then
      revalVectors%p_RvectorData(revalVectors%ncount)%nmaxDerivativeIdx = nsubarrays
    end if
    
    ! Create an interleaved memory block?
    if (present(nvar)) then
      revalVectors%p_RvectorData(revalVectors%ncount)%bisInterleaved = .true.
      revalVectors%p_RvectorData(revalVectors%ncount)%nvar = nvar
    end if
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fev2_releaseVectorList(revalVectors)

!<description>
  ! Releases a list of vectors.
!</description>

!<inputoutput>
  ! List to be cleaned up
  type(t_fev2Vectors), intent(inout) :: revalVectors
!</inputoutput>

!</subroutine>
  
    if (associated(revalVectors%p_RvectorData)) then
      deallocate(revalVectors%p_RvectorData)
    end if
    
    revalVectors%ncount = 0

  end subroutine
  
  !****************************************************************************

!<subroutine>

  subroutine fev2_initVectorEval(revalVectors,rfemDataBlocks)

!<description>
  ! Auxiliary subroutine. 
  ! Initialises the evaluation of vectors representing nonlinearities.
!</description>

!<input>
  ! Structure for all involved FEM spaces.
  type(t_fev2FemDataBlocks), intent(in) :: rfemDataBlocks
!</input>

!<inputoutput>
  ! Vector structure to be initialised
  type(t_fev2Vectors), intent(inout) :: revalVectors
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i, nmaxDerivativeIdx
    type(t_fev2VectorData), pointer :: p_rvectorData
    
    ! For every vector to evaluate, allocate memory for the
    ! values in the points.
    
    do i=1,revalVectors%ncount
    
      p_rvectorData => revalVectors%p_RvectorData(i)

      if (p_rvectorData%nmaxDerivative .ge. 0) then
        ! Standard vector with data.
        
        ! Maximum index in Bder
        call fev2_getBderSize(&
            p_rvectorData%p_rvector%p_rspatialDiscr%p_rtriangulation%ndim,&
            p_rvectorData%nmaxDerivative,nmaxDerivativeIdx)
        p_rvectorData%nmaxDerivativeIdx = nmaxDerivativeIdx
      
        ! Remember the index of the discretisation in the list
        p_rvectorData%iidxFemData = &
            containsDiscr (rfemDataBlocks,&
                           p_rvectorData%p_rvector%p_rspatialDiscr)
                           
        ! Get the dimension of the underlying FEM space
        p_rvectorData%ndimfe = rfemDataBlocks%p_RfemData(p_rvectorData%iidxFemData)%ndimfe
      else
        ! Dummy vector, just allocate memory in advance.
        ! nmaxDerivativeIdx specifies the number of subarrays to allocate.
        nmaxDerivativeIdx = p_rvectorData%nmaxDerivativeIdx
      
      end if

    end do
        
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fev2_prepareVectorEval(revalVectors,revalElementSet)

!<description>
  ! Auxiliary subroutine. Called during the evaluaation.
  ! Prepares the evaluation of the vector data. Allocates/Deallocates/
  ! Reallocates temporary memory if necessary.
!</description>

!<input>
  ! Element evaluation set encapsuling coordinates on of the cubature
  ! points on the reference element, real elements, Jacobian determinants etc.
  type(t_evalElementSet), intent(in) :: revalElementSet
!</input>

!<inputoutput>
  ! Vector structure to be initialised
  type(t_fev2Vectors), intent(inout) :: revalVectors
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i, nmaxDerivativeIdx, npointsPerElement, nelements
    type(t_fev2VectorData), pointer :: p_rvectorData
    
    ! For every vector to evaluate, allocate memory for the
    ! values in the points.
    npointsPerElement = revalElementSet%npointsPerElement
    nelements = revalElementSet%nelements
    
    if ((revalVectors%npointsPerElement .ne. npointsPerElement) .or. &
        (revalVectors%nelements .ne. nelements)) then
      
      revalVectors%npointsPerElement = npointsPerElement
      revalVectors%nelements = nelements

      ! Allocate memory.
      do i=1,revalVectors%ncount

        p_rvectorData => revalVectors%p_RvectorData(i)

        nmaxDerivativeIdx = p_rvectorData%nmaxDerivativeIdx
      
        ! Deallocate memory if memory is allocated.
        if (associated(p_rvectorData%p_Ddata)) deallocate (p_rvectorData%p_Ddata)
        if (associated(p_rvectorData%p_DdataIntl)) deallocate (p_rvectorData%p_DdataIntl)
        if (associated(p_rvectorData%p_DdataVec)) deallocate (p_rvectorData%p_DdataVec)

        ! Allocate memory for the point values.
        if ((npointsPerElement .ne. 0) .and. (nelements .ne. 0)) then
          if (.not. p_rvectorData%bisInterleaved) then
          
            ! Get the dimension of the underlying FEM space.
            if (p_rvectorData%ndimfe .eq. 1) then
              ! Non-interleaved data
              allocate (p_rvectorData%p_Ddata(npointsPerElement,nelements,nmaxDerivativeIdx))
            else
              ! Non-interleaved data, vector valued basis functions.
              allocate (p_rvectorData%p_DdataVec(p_rvectorData%ndimfe,&
                  npointsPerElement,nelements,nmaxDerivativeIdx))
            end if
          else
            ! Interleaved data
            allocate (p_rvectorData%p_DdataIntl(&
                p_rvectorData%nvar,npointsPerElement,nelements,nmaxDerivativeIdx))
          end if
        end if
      
      end do
      
    end if
        
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fev2_doneVectorEval(revalVectors)

!<description>
  ! Cleans up an evaluation structure list.
!</description>

!<inputoutput>
  ! List to be cleaned up
  type(t_fev2Vectors), intent(inout) :: revalVectors
!</inputoutput>

!</subroutine>

    integer :: i
    
    do i=1,revalVectors%ncount
    
      ! Release memory.
      if (associated(revalVectors%p_RvectorData(i)%p_Ddata)) then
        deallocate(revalVectors%p_RvectorData(i)%p_Ddata)
      end if

      if (associated(revalVectors%p_RvectorData(i)%p_DdataIntl)) then
        deallocate(revalVectors%p_RvectorData(i)%p_DdataIntl)
      end if
      
      revalVectors%p_RvectorData(i)%iidxFemData = 0
    
    end do
        
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fev2_evaluateVectors(revalVectors,rfemDataBlocks)

!<description>
  ! Evaluates the vectors in revalVectors in all points
  ! specified in the element evaluation set.
  !
  ! The data for dummy vectors is not initialised.
!</description>

!<input>
  ! Structure for all involved FEM spaces.
  type(t_fev2FemDataBlocks), intent(in) :: rfemDataBlocks
!</input>

!<inputoutput>
  ! Vector structure.
  type(t_fev2Vectors), intent(inout) :: revalVectors
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ibas, ivector, npoints, nelements, iel, ipt, idimfe
    integer :: ndof, ideriv, nderiv, nvar, ivar, ndimfe
    real(DP), dimension(:,:,:,:), pointer :: p_Dbas
    real(DP) :: dval
    type(t_fev2FemData), pointer :: p_rfemData
    type(t_fev2VectorData), pointer :: p_rvectorData
    integer, dimension(:,:), pointer :: p_Idofs
    real(DP), dimension(:), pointer :: p_DvecData
    real(DP), dimension(:,:,:), pointer :: p_Ddata
    real(DP), dimension(:,:,:,:), pointer :: p_DdataIntl,p_DdataVec
    
    npoints = revalVectors%npointsPerElement
    nelements = revalVectors%nelements

    ! Loop through all vectors    
    do ivector=1,revalVectors%ncount
    
      ! Vector data
      p_rvectorData => revalVectors%p_RvectorData(ivector)
      
      ! Is this a vector to be evaluated or is it a
      ! dummy vector which data can be left uninitialised?
      if (associated(p_rvectorData%p_rvector)) then
    
        ! Corresponding FEM data
        p_rfemData => rfemDataBlocks%p_RfemData(p_rvectorData%iidxFemData)
        
        ! Get data arrays and array size
        ndof = p_rfemData%ndof
        ndimfe = p_rfemData%ndimfe
        p_Dbas => p_rfemData%p_Dbas
        p_Idofs => p_rfemData%p_Idofs
        nderiv = p_rvectorData%nmaxDerivativeIdx
        call lsyssc_getbase_double (p_rvectorData%p_rvector,p_DvecData)
        
        if (.not. p_rvectorData%bisInterleaved) then

          if (ndimfe .eq. 1) then
          
            p_Ddata => p_rvectorData%p_Ddata
        
            ! Loop over the derivatives, basis functions, sum up to the point value
            ! in every cubature point.
            do iel = 1,nelements
              do ideriv = 1,nderiv
                do ipt = 1,npoints
                
                  dval = 0.0_DP
                  do ibas = 1,ndof
                    dval = dval + &
                        p_DvecData(p_Idofs(ibas,iel))*p_Dbas(ibas,ideriv,ipt,iel)
                  end do
                  
                  ! Save the value
                  p_Ddata(ipt,iel,ideriv) = dval
                
                end do
              end do
            end do
            
          else
          
            p_DdataVec => p_rvectorData%p_DdataVec
          
            ! Loop over the derivatives, basis functions, sum up to the point value
            ! in every cubature point.
            do iel = 1,nelements
              do ideriv = 1,nderiv
                do ipt = 1,npoints
                
                  do idimfe = 0,ndimfe-1
                    dval = 0.0_DP
                    do ibas = 1,ndof
                      dval = dval + &
                          p_DvecData(p_Idofs(ibas,iel))*p_Dbas(ibas+idimfe*ndof,ideriv,ipt,iel)
                    end do
                    
                    ! Save the value
                    p_DdataVec(1+idimfe,ipt,iel,ideriv) = dval
                  end do
                
                end do
              end do
            end do
          end if            
        else
          
          ! Interleaved specification
          nvar = p_rvectorData%nvar

          p_DdataIntl => p_rvectorData%p_DdataIntl
        
          ! Loop over the derivatives, basis functions, sum up to the point value
          ! in every cubature point.
          do iel = 1,nelements
            do ideriv = 1,nderiv
              do ipt = 1,npoints
                do ivar = 1,nvar
              
                  dval = 0.0_DP
                  do ibas = 1,ndof
                    dval = dval + &
                        p_DvecData(nvar*(p_Idofs(ibas,iel)-1)+ivar)*p_Dbas(ibas,ideriv,ipt,iel)
                  end do
                  
                  ! Save the value
                  p_DdataIntl(ivar,ipt,iel,ideriv) = dval
                  
                end do
              
              end do
            end do
          end do
        
        end if
        
      end if
        
    end do
      
  end subroutine
  
!****************************************************************************

!<subroutine>

  subroutine fev2_prepareFemDataVecEval(revalVectors,rfemDataBlocks,ielementDistr)

!<description>
  ! Initialise a FEM structure based on the FEM spaces
  ! appearing in a t_fev2Vectors structure
  ! No memory is allocated.
!</description>

!<input>
  ! A t_fev2Vectors structure with a set of vectorsto be evaluated.
  type(t_fev2Vectors), intent(in) :: revalVectors
  
  ! Element distribution to use.
  integer, intent(in) :: ielementDistr
!</input>

!<output>
  ! Structure for all involved FEM spaces.
  ! If the structure is empty, new memory is allocated.
  ! Otherwise, new FEM structures are appended.
  type(t_fev2FemDataBlocks), intent(inout) :: rfemDataBlocks
!</output>

!</subroutine>

    ! local variables
    integer :: i

    ! Loop over all vectors to be evaluated  
    do i=1,revalVectors%ncount
      if (associated(revalVectors%p_RvectorData(i)%p_rvector)) then
        ! Add the FEM space if not done already.
        call fev2_prepareFemDataSVec(&
            revalVectors%p_RvectorData(i)%p_rvector,rfemDataBlocks,ielementDistr, &
            revalVectors%p_RvectorData(i)%nmaxDerivative)
      end if
    end do

  end subroutine

end module
