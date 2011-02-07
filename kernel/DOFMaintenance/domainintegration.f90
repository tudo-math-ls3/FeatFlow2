!##############################################################################
!# ****************************************************************************
!# <name> domainintegration </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic structures and routines for integrating
!# a function over a domain.
!#
!# Routines in this module:
!#
!# 1.) domint_initIntegration
!#     -> Initialises a domain integration structure.
!#
!# 2.) domint_initIntegrationByEvalSet
!#     -> Initialises a domain integration structures based on an
!#        element evaluation set.
!#
!# 3.) domint_doneIntegration
!#     -> Releases a domain integration structure.
!# </purpose>
!##############################################################################

module domainintegration

  use collection
  use element
  use fsystem
  use genoutput
  use transformation
  use triangulation
  
  implicit none
  
  private

!<types>
  
!<typeblock>

  ! This structure is used during integration over a domain by the integration
  ! routines. It is passed to callback routines to inform them, which elements
  ! are currently in progress, which points are currently in progress,
  ! what are the properties of the mappings between the reference element
  ! and all real elements in progress, etc.
  type t_domainIntSubset
  
    ! Maximum number of elements in each element set.
    integer :: nelements = 0
    
    ! Number of (cubature) points per element.
    integer :: npointsPerElement = 0

    ! The currently active element distribution in the discretisation.
    ! Allows the routine to get the currently active element type for
    ! trial and test functions.
    integer :: ielementDistribution = 0
    
    ! Start index of the current element block in the current element 
    ! distribution ielementDistribution of the discretisation. 
    ! If this is =1, e.g., this is the very first element block
    ! that is currently being integrated.
    integer :: ielementStartIdx = 0

    ! The element set that is currently in progress by the integration 
    ! routine.
    integer, dimension(:), pointer :: p_Ielements => null()
    
    ! The orientation of the elements of the current element set.
    ! Typically, this pointer is not associated since the orientation of elements
    ! is not needed. At the boundary, it may be useful to know which edge/face
    ! of the element is located at the boundary
    integer, dimension(:), pointer :: p_IelementOrientation => null()

    ! If p_IdofsTrial is assigned, this is an element identifier that
    ! indicates the trial space.
    integer(I32) :: celement = 0
    
    ! An array containing the the degrees of freedom on all the
    ! elements. For multilinear forms (bilinear, trilinear), this is a pointer
    ! to the DOF`s of the trial space.
    ! DIMENSION(#dofPerElement,nelements)
    integer, dimension(:,:), pointer :: p_IdofsTrial => null()
    
    ! An element evaluation set structure that contains all information
    ! needed to evaluate the finite element on all elements in p_Ielements.
    type(t_evalElementSet), pointer :: p_revalElementSet => null()

    !<!--
    ! Information which may be shared with the element evaluation set
    ! -->
    
    ! A list of the corner vertices of all elements in progress.
    ! array [1..dimension,1..#vertices per element,1..Number of elements] of double
    real(DP), dimension(:,:,:), pointer :: p_Dcoords => null()
    
    ! A list of points in coordinates on the reference element.
    ! On each element in the current set of elements, this gives the
    ! coordinates of the cubature points on the reference element.
    ! Remark: As long as the same cubature formula is used on all
    ! elements, the coordinates here are the same for each element.
    ! array [1..dimension,1..npointsPerElement,1..Number of elements] of double
    real(DP), dimension(:,:,:), pointer :: p_DcubPtsRef => null()

    ! A list of points, corresponding to DcubPtsRef, in real coordinates.
    ! On each element in the current set of elements, this gives the
    ! coordinates of the cubature points on the real element.
    ! array [1..dimension,1..npointsPerElement,1..Number of elements] of double
    real(DP), dimension(:,:,:), pointer :: p_DcubPtsReal => null()

    ! The Jacobian matrix of the mapping between the reference and each
    ! real element, for all points on all elements in progress.
    ! array [1..dimension*dimension,1..npointsPerElement,1..Number of elements)
    real(DP), dimension(:,:,:),pointer :: p_Djac => null()
    
    ! The Jacobian determinant of the mapping of each point from the
    ! reference element to each real element in progress.
    ! array [1..npointsPerElement,1..Number of elements]
    real(DP), dimension(:,:), pointer :: p_Ddetj => null()
    
  end type
  
  public :: t_domainIntSubset
  
!</typeblock>

!</types>

  public :: domint_initIntegration
  public :: domint_initIntegrationByEvalSet
  public :: domint_doneIntegration
  
contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine domint_initIntegration (rintSubset,nelements,npointsPerElement,&
                                     icoordSystem,ndimSpace,nverticesPerElement)
  
!<description>
  ! This routine initialises a t_domainIntSubset structure. Memory is allocated
  ! for all entries in the structure.
  ! After the domain integration, the structure can be released with
  ! doneIntegration.
!</description>

!<input>
  ! Maximum number of elements in each element block.
  integer, intent(in) :: nelements
  
  ! Number of points per element
  integer, intent(in) :: npointsPerElement

  ! Coordinate system identifier. One of the TRAFO_CS_xxxx constants. Defines
  ! the type of the coordinate system that is used for specifying the coordinates
  ! on the reference element.
  integer(I32), intent(in) :: icoordSystem
  
  ! Number of space dimensions. Either NDIM2D or NDIM3D for 2D or 3D, 
  ! respectively.
  integer, intent(in) :: ndimSpace
  
  ! Number of vertices per element that are necessary to specify the 
  ! transformation from the reference to the real element
  integer, intent(in) :: nverticesPerElement
!</input>

!<output>
  ! An integration subset structure, initialised according to the parameters.
  type(t_domainIntSubset), intent(out) :: rintSubset
!</output>

!</subroutine>

    ! Initialise constants in the structure
    rintSubset%nelements            = nelements
    rintSubset%npointsPerElement    = npointsPerElement
    rintSubset%ielementDistribution = 1
    rintSubset%ielementStartIdx     = 1
    
    ! Allocate memory for the structures:
    !
    ! Allocate memory for corner coordinates.
    allocate(rintSubset%p_DCoords(ndimSpace,nverticesPerElement,nelements))

    ! Allocate arrays accepting cubature point coordinates.
    ! It is at most as large as number of elements or length
    ! of the element set.
    ! Check the coordinate system to find out what coordinate
    ! dimension to use
    select case (icoordSystem)
    case (TRAFO_CS_BARY2DTRI)
      allocate(rintSubset%p_DcubPtsRef(3,npointsPerElement,nelements))
      
    case (TRAFO_CS_BARY3DTETRA)
      allocate(rintSubset%p_DcubPtsRef(4,npointsPerElement,nelements))

    case (TRAFO_CS_REF2DTRI,TRAFO_CS_REF2DQUAD,&
          TRAFO_CS_REAL2DTRI,TRAFO_CS_REAL2DQUAD,TRAFO_CS_REF1D,&
          TRAFO_CS_REF3DTETRA,TRAFO_CS_REF3DHEXA,&
          TRAFO_CS_REAL3DTETRA,TRAFO_CS_REAL3DHEXA)
      allocate(rintSubset%p_DcubPtsRef(ndimSpace,npointsPerElement,nelements))
      
    case DEFAULT
      ! NULLIFY(rintSubset%p_DcubPtsRef)
      call output_line('Unknown coordinate system!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'domint_initIntegration')
      call sys_halt()
      
    end select

    ! The Jacobian matrix for the mapping from the reference to the real element
    ! has always length dimension^2
    allocate(rintSubset%p_Djac(ndimSpace*ndimSpace,npointsPerElement,nelements))
    
    ! The coordinate system on the real element has always the same dimension
    ! as the space.
    allocate(rintSubset%p_DcubPtsReal(ndimSpace,npointsPerElement,nelements))
    
    ! Allocate an array saving the coordinates of corner vertices of elements
    allocate(rintSubset%p_Ddetj(npointsPerElement,nelements))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine domint_initIntegrationByEvalSet (revalElementSet,rintSubset)
  
!<description>
  ! This routine initialises a t_domainIntSubset structure based on an
  ! element evaluation set. revalElementSet defines all arrays (Jacobian
  ! determinants, coordinates of points,...) where to evaluate.
  ! The routine will prepare the integration subset to fit to the
  ! element evaluation set. A pointer to the set is saved to rintSubset.
!</description>

!<input>
  ! Element evaluation set, the integration structure should be based on.
  type(t_evalElementSet), target :: revalElementSet
!</input>

!<output>
  ! An integration subset structure, initialised according to the parameters.
  type(t_domainIntSubset), intent(out) :: rintSubset
!</output>

!</subroutine>

    ! Initialise constants in the structure
    rintSubset%nelements            = revalElementSet%nelements
    rintSubset%npointsPerElement    = revalElementSet%npointsPerElement
    rintSubset%ielementDistribution = 1
    rintSubset%ielementStartIdx     = 1
    
    ! Let the pointers of the integration structure point to
    ! the structures of the evaluation set. We do not need the same arrays twice.
    rintSubset%p_Dcoords => revalElementSet%p_Dcoords
    rintSubset%p_DcubPtsRef => revalElementSet%p_DpointsRef
    rintSubset%p_DcubPtsReal => revalElementSet%p_DpointsReal
    rintSubset%p_Djac => revalElementSet%p_Djac
    rintSubset%p_Ddetj => revalElementSet%p_Ddetj
    
    ! Remember the evaluation set. This will prevent the DONE-routine from
    ! releasing the above pointers.
    rintSubset%p_revalElementSet => revalElementSet

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine domint_doneIntegration (rintSubset)
  
!<description>
  ! This routine reases memory allocated in domint_initIntegration
!</description>

!<inputoutput>
  ! The integation structure to be cleaned up.
  type(t_domainIntSubset), intent(inout) :: rintSubset
!</inputoutput>

!</subroutine>

    if (.not. associated(rintSubset%p_revalElementSet)) then
      ! Deallocate an array saving the coordinates of corner vertices of elements
      deallocate(rintSubset%p_Ddetj)
      if (associated(rintSubset%p_Djac)) deallocate(rintSubset%p_Djac)

      ! Deallocate arrays accepting cubature point coordinates.
      deallocate(rintSubset%p_DcubPtsReal)
      if (associated(rintSubset%p_DcubPtsRef)) deallocate(rintSubset%p_DcubPtsRef)

      ! Deallocate memory for corner coordinates
      deallocate(rintSubset%p_DCoords)
    else
      ! The pointers belong to the evaluation set, so we must not deallocate here!
      nullify(rintSubset%p_Ddetj)
      nullify(rintSubset%p_Djac)
      nullify(rintSubset%p_DcubPtsReal)
      nullify(rintSubset%p_DcubPtsRef)
      nullify(rintSubset%p_DCoords)
    end if

  end subroutine

end module
