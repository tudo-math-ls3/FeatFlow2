!##############################################################################
!# ****************************************************************************
!# <name> domainintegration </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic structures and routines for integrating
!# a function over a domain.
!# </purpose>
!##############################################################################

MODULE domainintegration

  USE fsystem
  USE triangulation
  USE collection
  
  IMPLICIT NONE

!<types>
  
!<typeblock>

  ! This structure is used during integration over a domain by the integration
  ! routines. It's passed to callback routines to inform them, which elements
  ! are currently in progress, which points are currently in progress,
  ! what are the properties of the mappings between the reference element
  ! and all real elements in progress, etc.
  TYPE t_domainIntSubset
  
    ! Maximum number of elements in each element set.
    INTEGER                                       :: nelements
    
    ! Number of (cubature) points per element.
    INTEGER                                       :: npointsPerElement

    ! The currently active element distribution in the discretisation.
    ! Allows the routine to get the currently active element type for
    ! trial and test functions.
    INTEGER                                       :: ielementDistribution
    
    ! Start index of the current element block in the current element 
    ! distribution ielementDistribution of the discretisation. 
    ! If this is =1, e.g., this is the very first element block
    ! that is currently being integrated.
    INTEGER(PREC_ELEMENTIDX)                      :: ielementStartIdx

    ! The element set that is currently in progress by the integration 
    ! routine.
    INTEGER(I32), DIMENSION(:), POINTER           :: p_Ielements
    
    ! A list of the corner vertices of all elements in progress.
    ! array [1..NDIM2D,1..TRIA_MAXNVE2D,1..Number of elements] of double
    REAL(DP), DIMENSION(:,:,:), POINTER           :: p_Dcoords
    
    ! A list of points in coordinates on the reference element.
    ! On each element in the current set of elements, this gives the
    ! coordinates of the cubature points on the reference element.
    ! Remark: As long as the same cubature formula is used on all
    !  elements, the coordinates here are the same for each element.
    ! array [1..dimension,1..npointsPerElement,1..Number of elements] of double
    REAL(DP), DIMENSION(:,:,:), POINTER           :: p_DcubPtsRef

    ! A list of points, corresponding to DcubPtsRef, in real coordinates.
    ! On each element in the current set of elements, this gives the
    ! coordinates of the cubature points on the real element.
    ! array [1..dimension,1..npointsPerElement,1..Number of elements] of double
    REAL(DP), DIMENSION(:,:,:), POINTER           :: p_DcubPtsReal

    ! The Jacobian matrix of the mapping between the reference and each
    ! real element, for all points on all elements in progress.
    ! array [1..dimension*dimension,1..npointsPerElement,1..Number of elements)
    REAL(DP), DIMENSION(:,:,:),POINTER            :: p_Djac
    
    ! The Jacobian determinant of the mapping of each point from the
    ! reference element to each real element in progress.
    ! array {1..npointsPerElement,1..Number of elements)
    REAL(DP), DIMENSION(:,:), POINTER             :: p_Ddetj

  END TYPE
  
!</typeblock>

!</types>

CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE domint_initIntegration (rintSubset,nelements,npointsPerElement,&
                                     icoordSystem,ndimSpace)
  
!<description>
  ! This routine initialises a t_domainIntSubset structure. Memory is allocated
  ! for all entries in the structure.
  ! After the domain integration, the structure can be released with
  ! doneIntegration.
!</description>

!<input>
  ! Maximum number of elements in each element block.
  INTEGER, INTENT(IN) :: nelements
  
  ! Number of points per element
  INTEGER, INTENT(IN) :: npointsPerElement

  ! Coordinate system identifier. One of the TRAFO_CS_xxxx constants. Defines
  ! the type of the coordinate system that is used for specifying the coordinates
  ! on the reference element.
  INTEGER, INTENT(IN) :: icoordSystem
  
  ! Number of space dimensions. Either NDIM2D or NDIM3D for 2D or 3D, 
  ! respectively.
  INTEGER, INTENT(IN) :: ndimSpace
!</input>

!<output>
  ! An integration subset structure, initialised according to the parameters.
  TYPE(t_domainIntSubset), INTENT(OUT) :: rintSubset
!</output>

    ! Initialise constants in the structure
    rintSubset%nelements            = nelements
    rintSubset%npointsPerElement    = npointsPerElement
    rintSubset%ielementDistribution = 1
    rintSubset%ielementStartIdx     = 1
    
    ! Allocate memory for the structures:
    !
    ! Allocate memory for corner coordinates.
    ALLOCATE(rintSubset%p_DCoords(ndimSpace,TRIA_MAXNVE2D,nelements))

    ! Allocate arrays accepting cubature point coordinates.
    ! It's at most as large as number of elements or length
    ! of the element set.
    ! Check the coordinate system to find out what coordinate
    ! dimension to use
    SELECT CASE (icoordSystem)
    CASE (TRAFO_CS_BARY2DTRI)
      ALLOCATE(rintSubset%p_DcubPtsRef(3,npointsPerElement,nelements))
      
    CASE (TRAFO_CS_REF2DTRI,TRAFO_CS_REF2DQUAD,&
          TRAFO_CS_REAL2DTRI,TRAFO_CS_REAL2DQUAD)
      ALLOCATE(rintSubset%p_DcubPtsRef(ndimSpace,npointsPerElement,nelements))
      
    CASE DEFAULT
      NULLIFY(rintSubset%p_DcubPtsRef)
      
    END SELECT

    ! The Jacobian matrix for the mapping from the reference to the real element
    ! has always length dimension^2
    ALLOCATE(rintSubset%p_Djac(ndimSpace*ndimSpace,npointsPerElement,nelements))
    
    ! The coordinate system on the real element has always the same dimension
    ! as the space.
    ALLOCATE(rintSubset%p_DcubPtsReal(ndimSpace,npointsPerElement,nelements))
    
    ! Allocate an array saving the coordinates of corner vertices of elements
    ALLOCATE(rintSubset%p_Ddetj(npointsPerElement,nelements))

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE domint_doneIntegration (rintSubset)
  
!<description>
  ! This routine reases memory allocated in domint_initIntegration
!</description>

!<inputoutput>
  ! The integation structure to be cleaned up.
  TYPE(t_domainIntSubset), INTENT(INOUT) :: rintSubset
!</inputoutput>

    ! Deallocate an array saving the coordinates of corner vertices of elements
    DEALLOCATE(rintSubset%p_Ddetj)
    IF (ASSOCIATED(rintSubset%p_Djac)) DEALLOCATE(rintSubset%p_Djac)

    ! Deallocate arrays accepting cubature point coordinates.
    DEALLOCATE(rintSubset%p_DcubPtsReal)
    IF (ASSOCIATED(rintSubset%p_DcubPtsRef)) DEALLOCATE(rintSubset%p_DcubPtsRef)

    ! Deallocate memory for corner coordinates
    DEALLOCATE(rintSubset%p_DCoords)

  END SUBROUTINE

END MODULE
