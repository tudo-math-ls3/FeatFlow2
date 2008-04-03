!##############################################################################
!# ****************************************************************************
!# <name> elementpreprocessing </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines that are usually executed to prepare the 
!# evaluation of a finite element in one or multiple points.
!#
!# The routines here will gather information from the triangulation, boundary
!# and whatever and set up evaluation structures. Calling the element routines
!# with these evaluation structures will then evaluate the element there.
!#
!# 1.) elprep_prepareSetForEvaluation
!#     -> Prepares an element for the evaluation on multiple points on multiple
!#        cells.
!#
!# 2.) elprep_releaseElementSet
!#     -> Releases an evaluation set for multiple points on multiple elements.
!#
!#  Frequently Asked Questions
!# ----------------------------
!#
!# 1.) How to evaluate a finite element? What are the routines here for?
!#
!#  Well, the routines in this module prepare the evaluation of basis
!#  functions. Imagine, you have
!#   - a mesh
!#   - a set of points
!#  and now you want to evaluate the finite element in these points. For that
!#  purpose, you have to call three functions:
!#
!#  a) cevaluationTag = elem_getEvaluationTag(ielement-type)
!#  b) elprep_prepareXXXX (eval-structure,mesh,points,...)
!#  c) elem_generic_XXXX (eval-structure,derivatives,values)
!#
!#  a) will at first grab from the element the so called 'evaluation tag'.
!#  This integer encodes a set of information that specifies 'what the element
!#  needs' to be evaluated; e.g. whether it needs coordinates on the reference
!#  element, on the real element, if mapping between reference and real
!#  element is needed etc.
!#
!#  b) uses this tag to create a so called 'evaluation structure' eval-structure.
!#  This structure contains all information that are necessary to evaluate
!#  the element -- e.g. it gathers cell information from the triangulation,
!#  writes coordinates of the evaluation points to the structure etc.
!#
!#  c) finally evaluates the basis functions using the prepared structure.
!#  The results are saved to 'values' in this example.
!#
!#  There's no harm in calling a elprep_prepareXXXX-routine more than once;
!#  each call overwrites previously generated information, amd memory is
!#  only reallocated if it's not enough. Nevertheless, when the evaluation
!#  is finished elprep_releaseXXXX shall be called to release allocated memory.
!#
!# </purpose>
!##############################################################################

MODULE elementpreprocessing

  USE fsystem
  USE triangulation
  USE element

CONTAINS

  !************************************************************************
  
!<subroutine>  

  SUBROUTINE elprep_prepareSetForEvaluation (revalElementSet, cevaluationTag, &
      rtriangulation, IelementList, ctrafoType, Dpoints, DpointsRef, DpointsReal)

!<description>
  ! This subroutine prepares a t_evalElementSet structure to be used for
  ! the evaluation of a finite element in a set of cells.
  ! Dpoints contains a list of coordinates on the reference element where to
  ! evaluate. These points are mapped onto all elements in the list
  ! IelementList from the triangulation rtriangulation.
  ! cevaluationTag specifies an 'evaluation tag' that defines which
  ! information must be prepared by this routine; that tag can be obtained
  ! by asking the finite element what it needs by calling 
  ! elem_getEvaluationTag.
!</description>

!<input>
  ! Evaluation tag. This is a bitfield that specifies which information is
  ! prepared in the element set; it's a combination of EL_EVLTAG_XXXX-constants.
  !
  ! Note: If EL_EVLTAG_REFCOORDS is not specified in this tag, the coordinates on
  ! the reference element are assumed to be initialised! Dpoints is ignored
  ! in that case. This can be used to calculate this information only once
  ! in a first call while then using the same set of reference coordinates
  ! for all subsequent calls.
  INTEGER(I32), INTENT(IN) :: cevaluationTag

  ! Underlying triangulation of the domain
  TYPE(t_triangulation), INTENT(IN) :: rtriangulation
  
  ! List of elements in the mesh where the integration is to be performed.
  INTEGER, DIMENSION(:), INTENT(IN) :: IelementList
  
  ! Type of transformation from the reference element to the real element.
  INTEGER, INTENT(IN) :: ctrafoType
  
  ! OPTIONAL: A set of npointsPerElement tuples (x,y) (or (x,y,z) in 3D) of the 
  ! points where to evaluate. These coordinates define on the reference the 
  ! coordinates of the cubature points where to evaluate the element.
  ! DIMENSION(ndimension,npointsPerElement)
  ! Ignored if EL_EVLTAG_REFCOORDS is not specified in cevaluationTag.
  REAL(DP), DIMENSION(:,:), OPTIONAL :: Dpoints
  
  ! OPTIONAL: An array with coordinates of the points where to evaluate,
  ! on the real element. For each element, a set of points can be
  ! specified here.
  ! A target to this array is saved in the revalElementSet structure. The array
  ! is assumed to be maintained by the caller.
  !
  ! If not specified, the routine will automatically set up that array
  ! if necessary using Dpoints.
  REAL(DP), DIMENSION(:,:,:), TARGET, OPTIONAL :: DpointsRef

  ! OPTIONAL: An array with coordinates of the points where to evaluate,
  ! on the reference element. For each element, a set of points can be
  ! specified here.
  ! A target to this array is saved in the revalElementSet structure. The array
  ! is assumed to be maintained by the caller.
  !
  ! If not specified, the routine will automatically set up that array
  ! if necessary using Dpoints.
  REAL(DP), DIMENSION(:,:,:), TARGET, OPTIONAL :: DpointsReal
!</input>
  
!<inputoutput>
  ! The element set that is to be initialised. If this is already initialised,
  ! previous information is overwritten.
  TYPE(t_evalElementSet), INTENT(INOUT) :: revalElementSet
!</inputoutput>

! </subroutine>

    INTEGER :: ndimRef,ndim,nverticesPerElement,nelements,ive,npointsPerElement
    INTEGER(I32), DIMENSION(:,:), POINTER :: p_ItwistIndexEdge
    INTEGER(PREC_ELEMENTIDX) :: iel,ielidx
    
    ! Fetch some information
    ndimRef = trafo_igetReferenceDimension(ctrafoType)
    ndim = trafo_igetDimension(ctrafoType)
    nverticesPerElement = trafo_igetNVE(ctrafoType)
    nelements = SIZE(IelementList)
    IF (revalElementSet%npointsPerElement .NE. 0) THEN 
      ! Structure already initialised
      npointsPerElement = revalElementSet%npointsPerElement
    ELSE IF (PRESENT(Dpoints)) THEN
      npointsPerElement = UBOUND(Dpoints,2)
    ELSE IF (PRESENT(DpointsRef)) THEN
      npointsPerElement = UBOUND(DpointsRef,2)
    ELSE IF (PRESENT(DpointsReal)) THEN
      npointsPerElement = UBOUND(DpointsReal,2)
    ELSE
      CALL output_line ('Cannot determine npointsPerElement!', &
          OU_CLASS_ERROR,OU_MODE_STD,'elprep_prepareSetForEvaluation')
    END IF
    
    NULLIFY(p_ItwistIndexEdge)
    IF (rtriangulation%h_ItwistIndex .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_int2d (rtriangulation%h_ItwistIndex,p_ItwistIndexEdge)
    END IF
  
    ! Probably save pointers to DpointsRef/DpointsReal.
    ! The default initialisation of the structure ensures that
    ! bforeignCoordsRef/bforeignCoordsReal is initialised, even in an
    ! uninitialised structure.
    IF (PRESENT(DpointsRef)) THEN
      ! Save a pointer to DpointRef
      IF (.NOT. revalElementSet%bforeignCoordsRef) THEN
        DEALLOCATE(revalElementSet%p_DpointsRef)
      END IF
      revalElementSet%p_DpointsRef => DpointsRef
      revalElementSet%bforeignCoordsRef = .TRUE.
    END IF 

    IF (PRESENT(DpointsReal)) THEN
      ! Save a pointer to DpointReal
      IF (.NOT. revalElementSet%bforeignCoordsReal) THEN
        DEALLOCATE(revalElementSet%p_DpointsReal)
      END IF
      revalElementSet%p_DpointsRef => DpointsReal
      revalElementSet%bforeignCoordsReal = .TRUE.
    END IF 
    
    ! Is the structure initialised?
    ! If yes, check the size of all the arrays. If it's different, deallocate
    ! and reallocate.
    IF (revalElementSet%npointsPerElement .NE. 0) THEN

      ! Corner coordinates for the transformation
      IF (ASSOCIATED(revalElementSet%p_Dcoords)) THEN
        IF (SIZE(revalElementSet%p_Dcoords) .LT. & 
            (ndim*nverticesPerElement*nelements)) THEN
          DEALLOCATE(revalElementSet%p_Dcoords)
        END IF
      END IF
    
      IF (IAND(cevaluationTag,EL_EVLTAG_REFCOORDS   ) .NE. 0) THEN
        IF (ASSOCIATED(revalElementSet%p_DpointsRef)) THEN
          IF (SIZE(revalElementSet%p_DpointsRef) .LT. & 
              (ndimRef*npointsPerElement*nelements)) THEN
            IF (.NOT. revalElementSet%bforeignCoordsRef) THEN
              DEALLOCATE(revalElementSet%p_DpointsRef)
            ELSE
              CALL output_line ('DpointsRef is too small!', &
                  OU_CLASS_ERROR,OU_MODE_STD,'elprep_prepareSetForEvaluation')
              CALL sys_halt()
            END IF
          END IF
        END IF
      END IF
      
      IF (IAND(cevaluationTag,EL_EVLTAG_REALCOORDS   ) .NE. 0) THEN
        IF (ASSOCIATED(revalElementSet%p_DpointsReal)) THEN
          IF (SIZE(revalElementSet%p_DpointsReal) .LT. &
              (ndim*npointsPerElement*nelements)) THEN
            IF (.NOT. revalElementSet%bforeignCoordsRef) THEN
              DEALLOCATE(revalElementSet%p_DpointsReal)
            ELSE
              CALL output_line ('DpointsReal is too small!', &
                  OU_CLASS_ERROR,OU_MODE_STD,'elprep_prepareSetForEvaluation')
              CALL sys_halt()
            END IF
          END IF
        END IF
      END IF

      IF (IAND(cevaluationTag,EL_EVLTAG_JAC   ) .NE. 0) THEN
        IF (ASSOCIATED(revalElementSet%p_Djac)) THEN
          IF (SIZE(revalElementSet%p_Djac) .LT. &
              (ndim*ndim*npointsPerElement*nelements)) THEN
            DEALLOCATE(revalElementSet%p_Djac)
          END IF
        END IF
      END IF
      
      IF (IAND(cevaluationTag,EL_EVLTAG_DETJ   ) .NE. 0) THEN
        IF (ASSOCIATED(revalElementSet%p_Ddetj)) THEN
          IF (SIZE(revalElementSet%p_Ddetj) .LT. (npointsPerElement*nelements)) THEN
            DEALLOCATE(revalElementSet%p_Ddetj)
          END IF
        END IF
      END IF
      
      IF (IAND(cevaluationTag,EL_EVLTAG_TWISTIDXEDGE   ) .NE. 0) THEN
        IF (ASSOCIATED(revalElementSet%p_ItwistIndexEdge)) THEN
          IF (SIZE(revalElementSet%p_ItwistIndexEdge) .LT. &
              (UBOUND(p_ItwistIndexEdge,1)*nelements)) THEN
            DEALLOCATE(revalElementSet%p_ItwistIndexEdge)
          END IF
        END IF
      END IF
      
      IF (IAND(cevaluationTag,EL_EVLTAG_TWISTIDXFACE   ) .NE. 0) THEN
      END IF
      
    END IF
    
    ! Set up basic parameters in the structure
    revalElementSet%npointsPerElement = npointsPerElement
    revalElementSet%nelements = nelements

    ! Allocate necessary data arrays
    
    ! Calculate the transformation
    IF (.NOT. ASSOCIATED(revalElementSet%p_Dcoords)) THEN
      ALLOCATE(revalElementSet%p_Dcoords(ndim,nverticesPerElement,nelements))
    END IF
    CALL trafo_getCoords_sim (ctrafoType,&
          rtriangulation,IelementList,revalElementSet%p_Dcoords)
          
    ! Coordinates on the reference element are always necessary.
    IF ((IAND(cevaluationTag,EL_EVLTAG_REFCOORDS   ) .NE. 0) .OR. &   
        (IAND(cevaluationTag,EL_EVLTAG_REALCOORDS   ) .NE. 0)) THEN
      IF (.NOT. ASSOCIATED(revalElementSet%p_DpointsRef)) THEN
        ALLOCATE(revalElementSet%p_DpointsRef(ndimRef,npointsPerElement,nelements))
      END IF
      
      ! Calculate the coordinates on the reference element
      IF (IAND(cevaluationTag,EL_EVLTAG_REFCOORDS   ) .NE. 0) THEN
        IF (.NOT. PRESENT(Dpoints)) THEN
          CALL output_line ('Dpoints not specified!', &
              OU_CLASS_ERROR,OU_MODE_STD,'elprep_prepareSetForEvaluation')
        END IF
        CALL calc_refCoords (&
            revalElementSet%p_DpointsRef,Dpoints,npointsPerElement,nelements)
      END IF
    END IF
    
    ! Prepare memory for the Jacobian matrix and/or determinant
    IF ((IAND(cevaluationTag,EL_EVLTAG_JAC   ) .NE. 0) .OR. &
        (IAND(cevaluationTag,EL_EVLTAG_DETJ  ) .NE. 0) .OR. &
        (IAND(cevaluationTag,EL_EVLTAG_REALCOORDS   ) .NE. 0)) THEN
        
      IF (.NOT. ASSOCIATED(revalElementSet%p_Djac)) THEN
        ALLOCATE(revalElementSet%p_Djac(ndim*ndim,npointsPerElement,nelements))
      END IF

      IF (.NOT. ASSOCIATED(revalElementSet%p_Ddetj)) THEN
        ALLOCATE(revalElementSet%p_Ddetj(npointsPerElement,nelements))
      END IF
      
      ! If real world coordinates must be calculated, we don't have to calculate
      ! the Jacobian stuff here; it's done below.
      IF (IAND(cevaluationTag,EL_EVLTAG_REALCOORDS   ) .EQ. 0) THEN
        CALL trafo_calctrafo_sim (ctrafoType,&
             nelements,npointsPerElement,revalElementSet%p_Dcoords,&
             revalElementSet%p_DpointsRef,revalElementSet%p_Djac,revalElementSet%p_Ddetj)
      END IF
      
    END IF
    
    ! Calculate the coordinates on the real element
    IF (IAND(cevaluationTag,EL_EVLTAG_REALCOORDS   ) .NE. 0) THEN
    
      IF (.NOT. ASSOCIATED(revalElementSet%p_DpointsReal)) THEN
        ALLOCATE(revalElementSet%p_DpointsReal(ndim,npointsPerElement,nelements))
      END IF
      
      ! Calculate the real world coordinates from the coordinates on the
      ! reference element. Calculate the Jacobial of the mapping from the
      ! reference element to the real elements.
      CALL trafo_calctrafo_sim (ctrafoType,&
          nelements,npointsPerElement,revalElementSet%p_Dcoords,&
          revalElementSet%p_DpointsRef,revalElementSet%p_Djac,revalElementSet%p_Ddetj,&
          revalElementSet%p_DpointsReal)
      
    END IF

    ! Get the twist indices or the orientation of edges.
    IF (IAND(cevaluationTag,EL_EVLTAG_TWISTIDXEDGE   ) .NE. 0) THEN
      
      ! Allocate memory
      IF (.NOT. ASSOCIATED(revalElementSet%p_ItwistIndexEdge)) THEN
        ALLOCATE(&
            revalElementSet%p_ItwistIndexEdge(UBOUND(p_ItwistIndexEdge,1),nelements))
      END IF
      
      ! Fetch the twist indices from the triangulation.
      DO ielidx = 1,nelements
        iel = IelementList(ielIdx)
        DO ive = 1,UBOUND(p_ItwistIndexEdge,1)
          revalElementSet%p_ItwistIndexEdge(ive,ielidx) = p_ItwistIndexEdge(ive,iel)
        END DO
      END DO
      
    END IF
    
    IF (IAND(cevaluationTag,EL_EVLTAG_TWISTIDXFACE   ) .NE. 0) THEN
      ! not yet defined
    END IF
    
  CONTAINS
  
    SUBROUTINE calc_refCoords (DpointsRef,Dpoints,npointsPerElement,nelements)
    
    ! Initialises the cubature point array on the reference element(s).
    
    ! Array with coordinates of points on the reference element.
    ! To be initialised.
    REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: DpointsRef
    
    ! A set of npointsPerElement tuples (x,y) (or (x,y,z) in 3D) of the points 
    ! where to evaluate. These coordinates define on the reference the coordinates
    ! of the cubature points where to evaluate the element.
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dpoints
    
    ! Number of points per element
    INTEGER, INTENT(IN) :: npointsPerElement
    
    ! Number of elements
    INTEGER, INTENT(IN) :: nelements
    
      ! local variables
      INTEGER :: iel,ipt,idim
    
      ! Copy the coordinates. They are the same for all elements.
      DO iel=1,nelements
        DO ipt=1,npointsPerElement
          DO idim=1,UBOUND(Dpoints,1)
            DpointsRef(idim,ipt,iel) = Dpoints(idim,ipt)
          END DO
        END DO
      END DO
    
    END SUBROUTINE
    
  END SUBROUTINE

  !************************************************************************
  
!<subroutine>  

  SUBROUTINE elprep_releaseElementSet (revalElementSet)

!<description>
  ! Releases all information in revalElementSet from memory. Cleans up the
  ! structure.
!</description>

!<inputoutput>
  ! The element set that is to be released.
  TYPE(t_evalElementSet), INTENT(INOUT) :: revalElementSet
!</inputoutput>

! </subroutine>

    ! Release all allocated information
    IF (ASSOCIATED(revalElementSet%p_Dcoords)) THEN
      DEALLOCATE(revalElementSet%p_Dcoords)
    END IF

    IF (ASSOCIATED(revalElementSet%p_DpointsRef) .AND. &
        (.NOT. revalElementSet%bforeignCoordsRef)) THEN
      DEALLOCATE(revalElementSet%p_DpointsRef)
    ELSE
      NULLIFY(revalElementSet%p_DpointsRef)
    END IF
  
    IF (ASSOCIATED(revalElementSet%p_DpointsReal) .AND. &
        (.NOT. revalElementSet%bforeignCoordsReal)) THEN
      DEALLOCATE(revalElementSet%p_DpointsReal)
    ELSE
      NULLIFY(revalElementSet%p_DpointsReal)
    END IF

    IF (ASSOCIATED(revalElementSet%p_Djac)) THEN
      DEALLOCATE(revalElementSet%p_Djac)
    END IF
  
    IF (ASSOCIATED(revalElementSet%p_Ddetj)) THEN
      DEALLOCATE(revalElementSet%p_Ddetj)
    END IF
  
    IF (ASSOCIATED(revalElementSet%p_ItwistIndexEdge)) THEN
      DEALLOCATE(revalElementSet%p_ItwistIndexEdge)
    END IF
  
    IF (ASSOCIATED(revalElementSet%p_ItwistIndexFace)) THEN
      DEALLOCATE(revalElementSet%p_ItwistIndexFace)
    END IF

    revalElementSet%npointsPerElement = 0
    revalElementSet%nelements = 0

  END SUBROUTINE

  
END MODULE 
