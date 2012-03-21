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
!# 6.) fev2_createFemData
!#     -> Full initialisation of a FEM data structure; allocate memory
!#
!# 7.) fev2_releaseFemData
!#     -> Cleanup of a FEM data structure
!#
!# 8.) fev2_evaluateFemData
!#     -> Save the values of FEM basis functions in a FEM data structure
!#
!# 9.) fev2_addVectorToEvalList
!#     -> Add a vector to a list of vectors to be evaluated simultaneously
!#        in a set of points on a set of elements
!#
!# 10.) fev2_releaseVectorList
!#      -> Release a vector evaluation list
!#
!# 11.) fev2_initVectorEval
!#      -> Initialise a vector evaluation list for the evaluation
!#         of vectors
!#
!# 12.) fev2_doneVectorEval
!#      -> Release a vector evaluation list
!#
!# 13.) fev2_evaluateVectors
!#      -> Evaluate all vectors in a vector evaluation list in a set of points
!#         on a set of elements
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
    
    ! Maximum derivative to be computed.
    ! =-1: Automatic, evaluate all up to the highest possible derivative.
    integer :: nmaxDerivative = -1
    
    ! Last index in Bder which is set to TRUE.
    integer :: nmaxDerivativeIdx = 0
    
    ! Type of transformation
    integer(I32) :: ctrafoType = TRAFO_ID_UNKNOWN
    
    ! Array to tell the element which derivatives to calculate
    logical, dimension(EL_MAXNDER) :: Bder = .false.

    ! Arrays for the basis function values in the cubature points
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

  ! This structure saves the values of a FEM function in a set
  ! of points on a set of elements. Used, e.g., for nonlinearities.
  type t_fev2VectorData
  
    ! Pointer to an array with values in all cubature points.
    !   Ddata(npointsPerElement,nelements,nmaxDerivativeIdx)
    ! The last index specifies the type of the derivative,
    ! so Ddata(:,:,DER_DERIV2D_X) specifies the first X-derivative.
    real(DP), dimension(:,:,:), pointer :: p_Ddata => null()
    
    ! Reference to the vector
    type(t_vectorScalar), pointer :: p_rvector => null()
    
    ! Maximum derivative to be computed
    integer :: nmaxDerivative = 0

    ! Last index of the corresponding Bder array which is set to TRUE.
    ! Specifies the last dimension in p_Ddata.
    integer :: nmaxDerivativeIdx = 0
    
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
  public :: fev2_createFemData
  public :: fev2_releaseFemData
  public :: fev2_evaluateFemData
  
  public :: fev2_addVectorToEvalList
  public :: fev2_releaseVectorList
  public :: fev2_initVectorEval
  public :: fev2_doneVectorEval
  public :: fev2_evaluateVectors

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
      select case (imaxDerivativeIdx)
      case (0)
        imaxDerivativeIdx = DER_FUNC2D
      case (1)
        imaxDerivativeIdx = DER_DERIV2D_Y
      case (2:)
        imaxDerivativeIdx = DER_DERIV2D_YY
      end select
      
    case (NDIM3D)
      select case (imaxDerivativeIdx)
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

  subroutine addDiscr (Rlist,nentries,rentry)
  ! Adds rentry to a list represented by an array
  type(t_fev2FemData), dimension(:), intent(inout) :: Rlist
  type(t_spatialDiscretisation), intent(in), target :: rentry
  integer, intent(inout) :: nentries
    nentries = nentries + 1
    Rlist(nentries)%p_rdiscr => rentry
  end subroutine

  !****************************************************************************
  
  integer function containsDiscr (Rlist,nentries,rentry)
  ! returns <> 0 (the index) if the list Rlist contains rentry
  type(t_fev2FemData), dimension(:), intent(in) :: Rlist
  type(t_spatialDiscretisation), intent(in), target :: rentry
  integer, intent(in) :: nentries
  
    integer :: i
    
    do i=1,nentries
      if (associated(Rlist(i)%p_rdiscr,rentry)) then
        containsDiscr = i
        return
      end if
    end do
    
    containsDiscr = 0
    
  end function
    
  !****************************************************************************

!<subroutine>

  subroutine fev2_prepareFemDataBMat(rmatrix,p_RfemData,ielementDistr, &
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
  ! Pointer to data of all involved FEM spaces.
  ! If this points to NULL, a new array is allocated.
  ! If this does not point to NULL, new FEM structures are appended.
  ! Memory is reallocated if necessary.
  type(t_fev2FemData), dimension(:), pointer :: p_RfemData
!</output>

!</subroutine>
  
    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rdiscrTrial,p_rdiscrTest
    integer :: i,j,k,ifirstNode
    
    integer, dimension(:,:), allocatable :: p_ImatDiscrTrial,p_ImatDiscrTest
    
    ! List of used discretisation structures
    type(t_fev2FemData), dimension(:), pointer :: p_RdiscrNodes
    integer :: ndiscrNodes
    
    ! Allocate memory for existence checks
    if (associated(p_RfemData)) then
    
      ifirstNode = size(p_RfemData)+1
      i = size(p_RfemData) + size(rmatrix%RmatrixBlock)
      allocate (p_RdiscrNodes(i))
      
      ! Copy the old content
      p_RdiscrNodes(1:size(p_RfemData)) = p_RfemData(:)

      ndiscrNodes = size(p_RfemData)
    
    else
    
      ! Only new content
      ifirstNode = 1
      i = size(rmatrix%RmatrixBlock)
      allocate (p_RdiscrNodes(i))
      ndiscrNodes = 0
      
    end if
    
    ! Arrays for remembering the discretisation structure index
    allocate(p_ImatDiscrTrial(rmatrix%nblocksPerCol,rmatrix%nblocksPerRow))
    allocate(p_ImatDiscrTest(rmatrix%nblocksPerCol,rmatrix%nblocksPerRow))

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

          ! Do we have the trial space already?
          p_ImatDiscrTrial(i,j) = containsDiscr (p_RdiscrNodes,ndiscrNodes,p_rdiscrTrial)
          if (p_ImatDiscrTrial(i,j) .eq. 0) then
            ! Remember the data for later checks
            call addDiscr (p_RdiscrNodes,ndiscrNodes,p_rdiscrTrial)
            p_ImatDiscrTrial(i,j) = ndiscrNodes
          end if

          ! What is with the test space? Shared? If not, do we have it?
          if (.not. associated(p_rdiscrTrial,p_rdiscrTest)) then
            p_ImatDiscrTest(i,j) = containsDiscr (p_RdiscrNodes,ndiscrNodes,p_rdiscrTest)
            if (p_ImatDiscrTest(i,j) .eq. 0) then
              ! Remember the data for later checks
              call addDiscr (p_RdiscrNodes,ndiscrNodes,p_rdiscrTest)
              p_ImatDiscrTest(i,j) = ndiscrNodes
            end if
          else
            p_ImatDiscrTest(i,j) = p_ImatDiscrTrial(i,j)
          end if

        end if
      
      end do
    
    end do
    
    if (ndiscrNodes .ne. 0) then

      ! Reallocate the FEM data structures and copy.
      if (associated(p_RfemData)) then
        deallocate(p_RfemData)
      end if
      
      allocate(p_RfemData(ndiscrNodes))
      p_RfemData(1:ndiscrNodes) = p_RdiscrNodes(1:ndiscrNodes)
      
      ! Initialise the actual data of the node
      do i = ifirstNode,ndiscrNodes
      
        p_RfemData(i)%celement = p_RfemData(i)%p_rdiscr%RelementDistr(ielementDistr)%celement
        p_RfemData(i)%ctrafoType = p_RfemData(i)%p_rdiscr%RelementDistr(ielementDistr)%ctrafoType
        p_RfemData(i)%ndof = elem_igetNDofLoc(p_RfemData(i)%celement)
        
      end do
      
      ! Initialise the maximum derivative to be used.
      do j=1,rmatrix%nblocksPerRow
      
        do i=1,rmatrix%nblocksPerCol
          if (lsysbl_isSubmatrixPresent (rmatrix,i,j)) then
          
            ! Maximum derivative for the test space
            if (present(RmaxDerivativeTest)) then
              k = p_ImatDiscrTrial(i,j)
              p_RfemData(k)%nmaxDerivative = &
                  max(p_RfemData(k)%nmaxDerivative,RmaxDerivativeTest(i,j))
            end if
            
            ! Maximum derivative for the trial space
            if (present(RmaxDerivativeTest)) then
              k = p_ImatDiscrTest(i,j)
              p_RfemData(k)%nmaxDerivative = &
                  max(p_RfemData(k)%nmaxDerivative,RmaxDerivativeTrial(i,j))
            end if
          
          end if
        end do
        
      end do
      
      ! Initialise the BDER array which is needed for the evaluation
      ! of the FEM functions.
      
      do i=ifirstNode,ndiscrNodes
      
        call fev2_initBder(&
            p_RfemData(i)%celement,p_RfemData(i)%nmaxDerivative,&
            p_RfemData(i)%nmaxDerivativeIdx,p_RfemData(i)%Bder)
        
      end do
      
    end if
    
    ! Release the temporary lists
    deallocate (p_RdiscrNodes)
    deallocate(p_ImatDiscrTrial)
    deallocate(p_ImatDiscrTest)
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fev2_prepareFemDataBVec(rvector,p_RfemData,ielementDistr, &
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
  ! Pointer to data of all involved FEM spaces.
  ! If this points to NULL, a new array is allocated.
  ! If this does not point to NULL, new FEM structures are appended.
  ! Memory is reallocated if necessary.
  type(t_fev2FemData), dimension(:), pointer :: p_RfemData
!</output>

!</subroutine>
  
    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rdiscrTest1,p_rdiscrTest
    integer :: i,k,ifirstnode
    
    integer, dimension(:), allocatable :: p_IvecDiscrTest

    ! List of used discretisation structures
    type(t_fev2FemData), dimension(:), pointer :: p_RdiscrNodes
    integer :: ndiscrNodes

    ! Allocate memory for existence checks
    if (associated(p_RfemData)) then
    
      ifirstNode = size(p_RfemData)+1
      i = size(p_RfemData) + rvector%nblocks
      allocate (p_RdiscrNodes(i))
      
      ! Copy the old content
      p_RdiscrNodes(1:size(p_RfemData)) = p_RfemData(:)

      ndiscrNodes = size(p_RfemData)
    
    else
    
      ! Only new content
      ifirstNode = 1
      i = rvector%nblocks
      allocate (p_RdiscrNodes(i))
      
    end if

    ! Arrays for remembering the discretisation structure index
    allocate(p_IvecDiscrTest(rvector%nblocks))

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

      ! Do we have the trial space already?
      p_IvecDiscrTest(i) = containsDiscr (p_RdiscrNodes,ndiscrNodes,p_rdiscrTest)
      if (p_IvecDiscrTest(i) .eq. 0) then
        ! Remember the data for later checks
        call addDiscr (p_RdiscrNodes,ndiscrNodes,p_rdiscrTest)
        p_IvecDiscrTest(i) = ndiscrNodes
      end if
    
    end do
    
    if (ndiscrNodes .ne. 0) then

      ! Reallocate the FEM data structures and copy.
      if (associated(p_RfemData)) then
        deallocate(p_RfemData)
      end if
      
      ! Prepare the FEM data structures. Copy the FEM data structures.
      allocate(p_RfemData(ndiscrNodes))
      p_RfemData(1:ndiscrNodes) = p_RdiscrNodes(1:ndiscrNodes)
      
      ! Initialise the actual data of the node
      do i = ifirstnode,ndiscrNodes
      
        p_RfemData(i)%celement = p_RfemData(i)%p_rdiscr%RelementDistr(ielementDistr)%celement
        p_RfemData(i)%ctrafoType = p_RfemData(i)%p_rdiscr%RelementDistr(ielementDistr)%ctrafoType
        p_RfemData(i)%ndof = elem_igetNDofLoc(p_RfemData(i)%celement)
        
      end do
      
      ! Initialise the maximum derivative to be used.
      do i=1,rvector%nblocks
        
        ! Maximum derivative for the test space
        if (present(RmaxDerivative)) then
          k = p_IvecDiscrTest(i)
          p_RfemData(k)%nmaxDerivative = &
              max(p_RfemData(k)%nmaxDerivative,RmaxDerivative(i))
        end if
        
      end do
        
      ! Initialise the BDER array which is needed for the evaluation
      ! of the FEM functions.
      
      do i=ifirstnode,ndiscrNodes
      
        call fev2_initBder(&
            p_RfemData(i)%celement,p_RfemData(i)%nmaxDerivative,&
            p_RfemData(i)%nmaxDerivativeIdx,p_RfemData(i)%Bder)
        
      end do
      
    end if
    
    ! Release the temporary lists
    deallocate (p_RdiscrNodes)
    
  contains
  
    subroutine addDiscr (Rlist,nentries,rentry)
    ! Adds rentry to a list represented by an array
    type(t_fev2FemData), dimension(:), intent(inout) :: Rlist
    type(t_spatialDiscretisation), intent(in), target :: rentry
    integer, intent(inout) :: nentries
      nentries = nentries + 1
      Rlist(nentries)%p_rdiscr => rentry
    end subroutine
    
    integer function containsDiscr (Rlist,nentries,rentry)
    ! returns <> 0 (the index) if the list Rlist contains rentry
    type(t_fev2FemData), dimension(:), intent(in) :: Rlist
    type(t_spatialDiscretisation), intent(in), target :: rentry
    integer, intent(in) :: nentries
    
      integer :: i
      
      do i=1,nentries
        if (associated(Rlist(nentries)%p_rdiscr,rentry)) then
          containsDiscr = i
          return
        end if
      end do
      
      containsDiscr = 0
      
    end function
    
  end subroutine

!****************************************************************************

!<subroutine>

  subroutine fev2_prepareFemDataSVec(rvector,p_RfemData,ielementDistr, &
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
  ! Pointer to data of all involved FEM spaces.
  ! If this points to NULL, a new array is allocated.
  ! If this does not point to NULL, new FEM structures are appended.
  ! Memory is reallocated if necessary.
  type(t_fev2FemData), dimension(:), pointer :: p_RfemData
!</output>

!</subroutine>
  
    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rdiscrTest
    integer :: i,ivecDiscrTest
    
    ! List of used discretisation structures
    type(t_fev2FemData), dimension(:), pointer :: p_RdiscrNodes
    integer :: ndiscrNodes

    ! Get the information about the trial and test spaces
    ! of this block.
    p_rdiscrTest  => rvector%p_rspatialDiscr
    
    ! Do we have the trial space already?
    ivecDiscrTest = containsDiscr (p_RdiscrNodes,ndiscrNodes,p_rdiscrTest)
    if (ivecDiscrTest .eq. 0) then
      
      ! Structure is new. Reallocate, add and initialise.
      if (associated(p_RfemData)) then
        ndiscrNodes = size(p_RfemData)+1
        allocate (p_RdiscrNodes(ndiscrNodes))
        p_RdiscrNodes(1:ndiscrNodes-1) = p_RfemData(:)
        deallocate(p_RfemData)
      else
        ndiscrNodes = 1
        allocate (p_RdiscrNodes(ndiscrNodes))
      end if
      
      p_RfemData => p_RdiscrNodes
      
      ndiscrNodes = ndiscrNodes-1
      call addDiscr (p_RdiscrNodes,ndiscrNodes,p_rdiscrTest)
    
      ! Initialise the content
      p_RfemData(ndiscrNodes)%celement = &
          p_RfemData(ndiscrNodes)%p_rdiscr%RelementDistr(ielementDistr)%celement
      p_RfemData(ndiscrNodes)%ctrafoType = &
          p_RfemData(ndiscrNodes)%p_rdiscr%RelementDistr(ielementDistr)%ctrafoType
      p_RfemData(ndiscrNodes)%ndof = elem_igetNDofLoc(p_RfemData(ndiscrNodes)%celement)

      ivecDiscrTest = ndiscrNodes
    end if

    ! Maximum derivative for the test space
    p_RfemData(ivecDiscrTest)%nmaxDerivative = &
        max(p_RfemData(ivecDiscrTest)%nmaxDerivative,imaxDerivative)
        
    call fev2_initBder(&
        p_RfemData(ivecDiscrTest)%celement,p_RfemData(i)%nmaxDerivative,&
        p_RfemData(ivecDiscrTest)%nmaxDerivativeIdx,p_RfemData(ivecDiscrTest)%Bder)

    ! Release the temporary lists
    deallocate (p_RdiscrNodes)
    
  contains
  
    subroutine addDiscr (Rlist,nentries,rentry)
    ! Adds rentry to a list represented by an array
    type(t_fev2FemData), dimension(:), intent(inout) :: Rlist
    type(t_spatialDiscretisation), intent(in), target :: rentry
    integer, intent(inout) :: nentries
      nentries = nentries + 1
      Rlist(nentries)%p_rdiscr => rentry
    end subroutine
    
    integer function containsDiscr (Rlist,nentries,rentry)
    ! returns <> 0 (the index) if the list Rlist contains rentry
    type(t_fev2FemData), dimension(:), intent(in) :: Rlist
    type(t_spatialDiscretisation), intent(in), target :: rentry
    integer, intent(in) :: nentries
    
      integer :: i
      
      do i=1,nentries
        if (associated(Rlist(nentries)%p_rdiscr,rentry)) then
          containsDiscr = i
          return
        end if
      end do
      
      containsDiscr = 0
      
    end function
    
  end subroutine

!****************************************************************************

!<subroutine>

  subroutine fev2_prepareFemDataSMat(rmatrix,p_RfemData,ielementDistr, &
      imaxDerivative)

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
  ! each FEM space is the default.
  integer, optional :: imaxDerivative
!</input>

!<output>
  ! Pointer to data of all involved FEM spaces.
  ! If this points to NULL, a new array is allocated.
  ! If this does not point to NULL, new FEM structures are appended.
  ! Memory is reallocated if necessary.
  type(t_fev2FemData), dimension(:), pointer :: p_RfemData
!</output>

!</subroutine>
  
    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rdiscrTest,p_rdiscrTrial
    integer :: ivecDiscrTest,ivecDiscrTrial
    
    ! List of used discretisation structures
    type(t_fev2FemData), dimension(:), pointer :: p_RdiscrNodes
    integer :: ndiscrNodes

    ! Get the information about the trial and test spaces
    ! of this block.
    p_rdiscrTest  => rmatrix%p_rspatialDiscrTest
    p_rdiscrTrial  => rmatrix%p_rspatialDiscrTest
    
    ! Do we have the trial space already?
    ivecDiscrTest = containsDiscr (p_RdiscrNodes,ndiscrNodes,p_rdiscrTest)
    
    if (ivecDiscrTest .eq. 0) then
      
      ! Structure is new. Reallocate, add and initialise.
      if (associated(p_RfemData)) then
        ndiscrNodes = size(p_RfemData)+1
        allocate (p_RdiscrNodes(ndiscrNodes))
        p_RdiscrNodes(1:ndiscrNodes-1) = p_RfemData(:)
        deallocate(p_RfemData)
      else
        ndiscrNodes = 1
        allocate (p_RdiscrNodes(ndiscrNodes))
      end if
      
      p_RfemData => p_RdiscrNodes
      
      ndiscrNodes = ndiscrNodes-1
      call addDiscr (p_RdiscrNodes,ndiscrNodes,p_rdiscrTest)
    
      ! Initialise the content
      p_RfemData(ndiscrNodes)%celement = &
          p_RfemData(ndiscrNodes)%p_rdiscr%RelementDistr(ielementDistr)%celement
      p_RfemData(ndiscrNodes)%ctrafoType = &
          p_RfemData(ndiscrNodes)%p_rdiscr%RelementDistr(ielementDistr)%ctrafoType
      p_RfemData(ndiscrNodes)%ndof = elem_igetNDofLoc(p_RfemData(ndiscrNodes)%celement)

      ivecDiscrTest = ndiscrNodes
    end if

    ivecDiscrTrial = containsDiscr (p_RdiscrNodes,ndiscrNodes,p_rdiscrTrial)

    if (ivecDiscrTrial .eq. 0) then
      
      ! Structure is new. Reallocate, add and initialise.
      if (associated(p_RfemData)) then
        ndiscrNodes = size(p_RfemData)+1
        allocate (p_RdiscrNodes(ndiscrNodes))
        p_RdiscrNodes(1:ndiscrNodes-1) = p_RfemData(:)
        deallocate(p_RfemData)
      else
        ndiscrNodes = 1
        allocate (p_RdiscrNodes(ndiscrNodes))
      end if
      
      p_RfemData => p_RdiscrNodes
      
      ndiscrNodes = ndiscrNodes-1
      call addDiscr (p_RdiscrNodes,ndiscrNodes,p_rdiscrTrial)
    
      ! Initialise the content
      p_RfemData(ndiscrNodes)%celement = &
          p_RfemData(ndiscrNodes)%p_rdiscr%RelementDistr(ielementDistr)%celement
      p_RfemData(ndiscrNodes)%ctrafoType = &
          p_RfemData(ndiscrNodes)%p_rdiscr%RelementDistr(ielementDistr)%ctrafoType
      p_RfemData(ndiscrNodes)%ndof = elem_igetNDofLoc(p_RfemData(ndiscrNodes)%celement)

      ivecDiscrTrial = ndiscrNodes
    end if

    ! Maximum derivative for the test space
    p_RfemData(ivecDiscrTest)%nmaxDerivative = &
        max(p_RfemData(ivecDiscrTest)%nmaxDerivative,imaxDerivative)

    p_RfemData(ivecDiscrTrial)%nmaxDerivative = &
        max(p_RfemData(ivecDiscrTrial)%nmaxDerivative,imaxDerivative)
        
    call fev2_initBder(&
        p_RfemData(ivecDiscrTest)%celement,p_RfemData(ivecDiscrTest)%nmaxDerivative,&
        p_RfemData(ivecDiscrTest)%nmaxDerivativeIdx,p_RfemData(ivecDiscrTest)%Bder)

    call fev2_initBder(&
        p_RfemData(ivecDiscrTrial)%celement,p_RfemData(ivecDiscrTrial)%nmaxDerivative,&
        p_RfemData(ivecDiscrTrial)%nmaxDerivativeIdx,p_RfemData(ivecDiscrTrial)%Bder)

    ! Release the temporary lists
    deallocate (p_RdiscrNodes)
    
  contains
  
    subroutine addDiscr (Rlist,nentries,rentry)
    ! Adds rentry to a list represented by an array
    type(t_fev2FemData), dimension(:), intent(inout) :: Rlist
    type(t_spatialDiscretisation), intent(in), target :: rentry
    integer, intent(inout) :: nentries
      nentries = nentries + 1
      Rlist(nentries)%p_rdiscr => rentry
    end subroutine
    
    integer function containsDiscr (Rlist,nentries,rentry)
    ! returns <> 0 (the index) if the list Rlist contains rentry
    type(t_fev2FemData), dimension(:), intent(in) :: Rlist
    type(t_spatialDiscretisation), intent(in), target :: rentry
    integer, intent(in) :: nentries
    
      integer :: i
      
      do i=1,nentries
        if (associated(Rlist(nentries)%p_rdiscr,rentry)) then
          containsDiscr = i
          return
        end if
      end do
      
      containsDiscr = 0
      
    end function
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fev2_createFemData(RfemData,npointsperelement,nelements)

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
  ! List all involved FEM spaces.
  type(t_fev2FemData), dimension(:), intent(inout), target :: RfemData
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_fev2FemData), pointer :: p_rfemData
    integer :: i

    ! Loop through all FEM data blocks.
    do i=1,size(RfemData)
    
      p_rfemData => RfemData(i)

      ! Allocate memory for the values of the basis functions    
      allocate(p_rfemData%p_Dbas(p_rfemData%ndof,&
          p_rfemData%nmaxDerivativeIdx,npointsperelement,nelements))
      
      ! Allocate memory for the degrees of freedoms on the
      ! elements to be processed.
      allocate(p_rfemData%p_Idofs(p_rfemData%ndof,nelements))
    
    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fev2_releaseFemData(RfemData)

!<description>
  ! Cleans up a set of FEM evaluation structures.
!</description>

!<inputoutput>
  ! Pointer to data of all involved FEM spaces.
  type(t_fev2FemData), dimension(:), intent(inout), target :: RfemData
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_fev2FemData), pointer :: p_rfemData
    integer :: i

    ! Release all allocated memory.
    do i=1,size(RfemData)
    
      p_rfemData => RfemData(i)
      deallocate(p_rfemData%p_Idofs)
      deallocate(p_rfemData%p_Dbas)
    
    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fev2_evaluateFemData(RfemData,revalElementSet)

!<description>
  ! Evaluates the FEM basis functions according to the evaluation
  ! set revalElementSet. The values of the basis functions are saved
  ! in the p_Dbas arrays of the FEM data structures in RfemData.
!</description>

!<input>
  ! Element evaluation set encapsuling coordinates on of the cubature
  ! points on the reference element, real elements, Jacobian determinants etc.
  type(t_evalElementSet), intent(in) :: revalElementSet
!</input>

!<inputoutput>
  ! Pointer to data of all involved FEM spaces.
  ! The p_Dbas arrays in all the structures are filled with the
  ! values of the FEM basis functions.
  ! The structure must have been initialised by fev2_createFemData.
  type(t_fev2FemData), dimension(:), intent(inout), target :: RfemData
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
    type(t_fev2FemData), pointer :: p_rfemData

    ! Loop through all FEM spaces to calculate the FEM basis functions
    do i=1,size(RfemData)
    
      p_rfemData => RfemData(i)

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
  type(t_vectorScalar), intent(inout), target :: rvector
  
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
        allocate(revalVectors%p_RvectorData(size(revalVectors%p_RvectorData)))
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

  subroutine fev2_initVectorEval(revalVectors,RfemData,revalElementSet)

!<description>
  ! Auxiliary subroutine. 
  ! Initialises the evaluation of vectors representing nonlinearities.
!</description>

!<input>
  ! List all involved FEM spaces that appear in the vectors
  ! collected in revalVectors
  type(t_fev2FemData), dimension(:), intent(in) :: RfemData
  
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
    integer :: i, nentries, nmaxDerivativeIdx, npointsPerElement, nelements
    
    ! For every vector to evaluate, allocate memory for the
    ! values in the points.
    
    nentries = size(RfemData)
    npointsPerElement = revalElementSet%npointsPerElement
    nelements = revalElementSet%nelements
    
    revalVectors%npointsPerElement = npointsPerElement
    revalVectors%nelements = nelements
    
    do i=1,revalVectors%ncount

      ! Maximum index in Bder
      call fev2_getBderSize(&
          revalVectors%p_RvectorData(i)%p_rvector%p_rspatialDiscr%p_rtriangulation%ndim,&
          revalVectors%p_RvectorData(i)%nmaxDerivative,nmaxDerivativeIdx)
      revalVectors%p_RvectorData%nmaxDerivativeIdx = nmaxDerivativeIdx
    
      ! Allocate memory for the point values
      allocate (revalVectors%p_RvectorData(i)%p_Ddata(&
          npointsPerElement,nelements,nmaxDerivativeIdx))
      
      ! Remember the index of the discretisation in the list
      revalVectors%p_RvectorData%iidxFemData = &
          containsDiscr (RfemData,nentries,&
                         revalVectors%p_RvectorData(i)%p_rvector%p_rspatialDiscr)
    
    end do
        
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
      
      revalVectors%p_RvectorData(i)%iidxFemData = 0
    
    end do
        
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fev2_evaluateVectors(revalVectors,RfemData)

!<description>
  ! Evaluates the vectors in revalVectors in all points
  ! specified in the element evaluation set.
!</description>

!<input>
  ! List all involved FEM spaces that appear in the vectors
  ! collected in revalVectors
  type(t_fev2FemData), dimension(:), intent(in), target :: RfemData
!</input>

!<inputoutput>
  ! Vector structure.
  type(t_fev2Vectors), intent(inout) :: revalVectors
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ibas, ivector, npoints, nelements, nentries, iel, ipt
    integer :: ndof, ideriv, nderiv
    real(DP), dimension(:,:,:,:), pointer :: p_Dbas
    real(DP) :: dval
    type(t_fev2FemData), pointer :: p_rfemData
    type(t_fev2VectorData), pointer :: p_rvectorData
    integer, dimension(:,:), pointer :: p_Idofs
    real(DP), dimension(:), pointer :: p_DvecData
    real(DP), dimension(:,:,:), pointer :: p_Ddata
    
    npoints = revalVectors%npointsPerElement
    nelements = revalVectors%nelements
    nentries = size(RfemData)

    ! Loop through all vectors    
    do ivector=1,revalVectors%ncount
    
      ! Vector data
      p_rvectorData => revalVectors%p_RvectorData(ivector)
    
      ! Corresponding FEM data
      p_rfemData => RfemData(p_rvectorData%iidxFemData)
      
      ! Get data arrays and array size
      ndof = p_rfemData%ndof
      p_Dbas => p_rfemData%p_Dbas
      p_Idofs => p_rfemData%p_Idofs
      p_Ddata => p_rvectorData%p_Ddata
      nderiv = p_rvectorData%nmaxDerivativeIdx
      call lsyssc_getbase_double (p_rvectorData%p_rvector,p_DvecData)
      
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
      
    end do
        
  end subroutine
  
end module
