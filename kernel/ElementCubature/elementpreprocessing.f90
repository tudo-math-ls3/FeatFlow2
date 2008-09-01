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
!# 1.) elprep_init
!#     -> Initialises an element evaluation structure with default values.
!#
!# 2.) elprep_prepareForEvaluation
!#     -> Prepares an element for the evaluation on a single point in a 
!#        single cell.
!#
!# 3.) elprep_prepareSetForEvaluation
!#     -> Prepares an element for the evaluation on multiple points on multiple
!#        cells.
!#
!# 4.) elprep_releaseElementSet
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

module elementpreprocessing

  use fsystem
  use triangulation
  use element
  
  implicit none

contains

  !************************************************************************
  
!<subroutine>  

  subroutine elprep_init (revalElementSet)

!<description>
  ! Initialises an element set with default values.
!</description>

!<output>
  ! The element set that is to be initialised. If this is already initialised,
  ! previous information is overwritten.
  type(t_evalElementSet), intent(OUT) :: revalElementSet
!</output>

! </subroutine>

    ! No commands here; use default initialisation of Fortran 90.

  end subroutine

  !************************************************************************
  
!<subroutine>  

  subroutine elprep_prepareSetForEvaluation (revalElementSet, cevaluationTag, &
      rtriangulation, IelementList, ctrafoType, Dpoints, DpointsRef, DpointsReal, Dcoords)

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
  ! Note: If EL_EVLTAG_REFPOINTS is not specified in this tag, the coordinates on
  ! the reference element are assumed to be initialised! Dpoints is ignored
  ! in that case. This can be used to calculate this information only once
  ! in a first call while then using the same set of reference coordinates
  ! for all subsequent calls.
  integer(I32), intent(IN) :: cevaluationTag

  ! Underlying triangulation of the domain
  type(t_triangulation), intent(IN) :: rtriangulation
  
  ! List of elements in the mesh where the integration is to be performed.
  integer, dimension(:), intent(IN) :: IelementList
  
  ! Type of transformation from the reference element to the real element.
  integer, intent(IN) :: ctrafoType
  
  ! OPTIONAL: A set of npointsPerElement tuples (x,y) (or (x,y,z) in 3D) of the 
  ! points where to evaluate. These coordinates define on the reference the 
  ! coordinates of the cubature points where to evaluate the element.
  ! DIMENSION(ndimension,npointsPerElement)
  ! Ignored if EL_EVLTAG_REFPOINTS is not specified in cevaluationTag.
  real(DP), dimension(:,:), optional :: Dpoints
  
  ! OPTIONAL: An array with coordinates of the points where to evaluate,
  ! on the reference element. For each element, a set of points can be
  ! specified here.
  ! A pointer to this array is saved in the revalElementSet structure. The array
  ! is assumed to be maintained by the caller.
  !
  ! If not specified, the routine will automatically set up that array
  ! using Dpoints.
  real(DP), dimension(:,:,:), target, optional :: DpointsRef

  ! OPTIONAL: An array with coordinates of the points where to evaluate,
  ! on the real element. For each element, a set of points can be
  ! specified here.
  ! A pointer to this array is saved in the revalElementSet structure. The array
  ! is assumed to be maintained by the caller.
  !
  ! If not specified, the routine will automatically set up that array
  ! using Dpoints.
  real(DP), dimension(:,:,:), target, optional :: DpointsReal
  
  ! OPTIONAL: An array with vertices defining the elements.
  ! A pointer to this array is saved in the revalElementSet structure. The array
  ! is assumed to be maintained by the caller.
  !
  ! If not specified, the routine will automatically set up that array
  ! using rtriangulation.
  real(DP), dimension(:,:,:), target, optional :: Dcoords
!</input>
  
!<inputoutput>
  ! The element set that is to be initialised. If this is already initialised,
  ! previous information is overwritten.
  type(t_evalElementSet), intent(INOUT) :: revalElementSet
!</inputoutput>

! </subroutine>

    integer :: ndimRef,ndim,nverticesPerElement,nelements,npointsPerElement
    integer(I32), dimension(:), pointer :: p_ItwistIndexEdge
    integer(PREC_ELEMENTIDX) :: iel,ielidx
    
    ! Fetch some information
    ndimRef = trafo_igetReferenceDimension(ctrafoType)
    ndim = trafo_igetDimension(ctrafoType)
    nverticesPerElement = trafo_igetNVE(ctrafoType)
    nelements = size(IelementList)
    if (revalElementSet%npointsPerElement .ne. 0) then 
      ! Structure already initialised
      npointsPerElement = revalElementSet%npointsPerElement
    else if (present(Dpoints)) then
      npointsPerElement = ubound(Dpoints,2)
    else if (present(DpointsRef)) then
      npointsPerElement = ubound(DpointsRef,2)
    else if (present(DpointsReal)) then
      npointsPerElement = ubound(DpointsReal,2)
    else
      call output_line ('Cannot determine npointsPerElement!', &
          OU_CLASS_ERROR,OU_MODE_STD,'elprep_prepareSetForEvaluation')
    end if
    
    ! Probably save pointers to DpointsRef/DpointsReal.
    ! The default initialisation of the structure ensures that
    ! bforeignPointsRef/bforeignPointsReal is initialised, even in an
    ! uninitialised structure.
    if (present(DpointsRef)) then
      ! Save a pointer to DpointRef
      if (.not. revalElementSet%bforeignPointsRef) then
        deallocate(revalElementSet%p_DpointsRef)
      end if
      revalElementSet%p_DpointsRef => DpointsRef
      revalElementSet%bforeignPointsRef = .true.
    end if 

    if (present(DpointsReal)) then
      ! Save a pointer to DpointReal
      if (.not. revalElementSet%bforeignPointsReal) then
        deallocate(revalElementSet%p_DpointsReal)
      end if
      revalElementSet%p_DpointsReal => DpointsReal
      revalElementSet%bforeignPointsReal = .true.
    end if 

    ! Is the structure initialised?
    ! If yes, check the size of all the arrays. If it's different, deallocate
    ! and reallocate.
    if (revalElementSet%npointsPerElement .ne. 0) then

      ! Corner coordinates for the transformation
      if (iand(cevaluationTag,EL_EVLTAG_COORDS      ) .ne. 0) then
        if (size(revalElementSet%p_Dcoords) .lt. & 
            (ndim*nverticesPerElement*nelements)) then
          if (.not. revalElementSet%bforeignCoords) then
            ! Release the old Dcoords to create a new one.
            deallocate(revalElementSet%p_Dcoords)
          else
            call output_line ('Dcoords is too small!', &
                OU_CLASS_ERROR,OU_MODE_STD,'elprep_prepareSetForEvaluation')
            call sys_halt()
          end if
        end if 
      end if
    
      if (iand(cevaluationTag,EL_EVLTAG_REFPOINTS   ) .ne. 0) then
        if (associated(revalElementSet%p_DpointsRef)) then
          if (size(revalElementSet%p_DpointsRef) .lt. & 
              (ndimRef*npointsPerElement*nelements)) then
            if (.not. revalElementSet%bforeignPointsRef) then
              deallocate(revalElementSet%p_DpointsRef)
            else
              call output_line ('DpointsRef is too small!', &
                  OU_CLASS_ERROR,OU_MODE_STD,'elprep_prepareSetForEvaluation')
              call sys_halt()
            end if
          end if
        end if
      end if
      
      if (iand(cevaluationTag,EL_EVLTAG_REALPOINTS   ) .ne. 0) then
        if (associated(revalElementSet%p_DpointsReal)) then
          if (size(revalElementSet%p_DpointsReal) .lt. &
              (ndim*npointsPerElement*nelements)) then
            if (.not. revalElementSet%bforeignPointsRef) then
              deallocate(revalElementSet%p_DpointsReal)
            else
              call output_line ('DpointsReal is too small!', &
                  OU_CLASS_ERROR,OU_MODE_STD,'elprep_prepareSetForEvaluation')
              call sys_halt()
            end if
          end if
        end if
      end if

      if (iand(cevaluationTag,EL_EVLTAG_JAC   ) .ne. 0) then
        if (associated(revalElementSet%p_Djac)) then
          if (size(revalElementSet%p_Djac) .lt. &
              (ndim*ndim*npointsPerElement*nelements)) then
            deallocate(revalElementSet%p_Djac)
          end if
        end if
      end if
      
      if (iand(cevaluationTag,EL_EVLTAG_DETJ   ) .ne. 0) then
        if (associated(revalElementSet%p_Ddetj)) then
          if (size(revalElementSet%p_Ddetj) .lt. (npointsPerElement*nelements)) then
            deallocate(revalElementSet%p_Ddetj)
          end if
        end if
      end if
      
      if (iand(cevaluationTag,EL_EVLTAG_TWISTIDXEDGE   ) .ne. 0) then

        if (associated(revalElementSet%p_ItwistIndexEdge)) then
          if (size(revalElementSet%p_ItwistIndexEdge) .lt. nelements) then
            deallocate(revalElementSet%p_ItwistIndexEdge)
          end if
        end if
      end if
      
      if (iand(cevaluationTag,EL_EVLTAG_TWISTIDXFACE   ) .ne. 0) then
      end if
      
    end if
    
    ! Set up basic parameters in the structure
    revalElementSet%npointsPerElement = npointsPerElement
    revalElementSet%nelements = nelements

    ! Allocate necessary data arrays
    
    ! Calculate the transformation
    if (iand(cevaluationTag,EL_EVLTAG_COORDS      ) .ne. 0) then
      if (.not. present(Dcoords)) then
        if (.not. associated(revalElementSet%p_Dcoords)) then
          allocate(revalElementSet%p_Dcoords(ndim,nverticesPerElement,nelements))
          revalElementSet%bforeignCoords = .false.
        end if
        call trafo_getCoords_sim (ctrafoType,&
              rtriangulation,IelementList,revalElementSet%p_Dcoords)
      else
        revalElementSet%p_Dcoords => Dcoords
        revalElementSet%bforeignCoords = .true.
      end if
    end if
          
    ! Coordinates on the reference element are always necessary.
    if ((iand(cevaluationTag,EL_EVLTAG_REFPOINTS   ) .ne. 0) .or. &   
        (iand(cevaluationTag,EL_EVLTAG_REALPOINTS   ) .ne. 0)) then
      if (.not. associated(revalElementSet%p_DpointsRef)) then
        allocate(revalElementSet%p_DpointsRef(ndimRef,npointsPerElement,nelements))
      end if
      
      ! Calculate the coordinates on the reference element
      if (iand(cevaluationTag,EL_EVLTAG_REFPOINTS   ) .ne. 0) then
        if (.not. present(Dpoints)) then
          call output_line ('Dpoints not specified!', &
              OU_CLASS_ERROR,OU_MODE_STD,'elprep_prepareSetForEvaluation')
        end if
        call calc_refCoords (&
            revalElementSet%p_DpointsRef,Dpoints,npointsPerElement,nelements)
      end if
    end if
    
    ! Prepare memory for the Jacobian matrix and/or determinant
    if ((iand(cevaluationTag,EL_EVLTAG_JAC   ) .ne. 0) .or. &
        (iand(cevaluationTag,EL_EVLTAG_DETJ  ) .ne. 0) .or. &
        (iand(cevaluationTag,EL_EVLTAG_REALPOINTS   ) .ne. 0)) then
        
      if (.not. associated(revalElementSet%p_Djac)) then
        allocate(revalElementSet%p_Djac(ndim*ndim,npointsPerElement,nelements))
      end if

      if (.not. associated(revalElementSet%p_Ddetj)) then
        allocate(revalElementSet%p_Ddetj(npointsPerElement,nelements))
      end if
      
      if (.not. associated(revalElementSet%p_Dcoords)) then
        call output_line ('Dcoords not available!', &
            OU_CLASS_ERROR,OU_MODE_STD,'elprep_prepareSetForEvaluation')
        call sys_halt()
      end if
      
      ! If real world coordinates must be calculated, we don't have to calculate
      ! the Jacobian stuff here; it's done below.
      if (iand(cevaluationTag,EL_EVLTAG_REALPOINTS   ) .eq. 0) then
        call trafo_calctrafo_sim (ctrafoType,&
             nelements,npointsPerElement,revalElementSet%p_Dcoords,&
             revalElementSet%p_DpointsRef,revalElementSet%p_Djac,revalElementSet%p_Ddetj)
      end if
      
    end if
    
    ! Calculate the coordinates on the real element
    if (iand(cevaluationTag,EL_EVLTAG_REALPOINTS   ) .ne. 0) then
    
      if (.not. associated(revalElementSet%p_DpointsReal)) then
        allocate(revalElementSet%p_DpointsReal(ndim,npointsPerElement,nelements))
      end if
      
      if (.not. associated(revalElementSet%p_Dcoords)) then
        call output_line ('Dcoords not available!', &
            OU_CLASS_ERROR,OU_MODE_STD,'elprep_prepareSetForEvaluation')
        call sys_halt()
      end if

      ! Calculate the real world coordinates from the coordinates on the
      ! reference element. Calculate the Jacobial of the mapping from the
      ! reference element to the real elements.
      call trafo_calctrafo_sim (ctrafoType,&
          nelements,npointsPerElement,revalElementSet%p_Dcoords,&
          revalElementSet%p_DpointsRef,revalElementSet%p_Djac,revalElementSet%p_Ddetj,&
          revalElementSet%p_DpointsReal)
      
    end if

    ! Get the twist indices or the orientation of edges.
    if (iand(cevaluationTag,EL_EVLTAG_TWISTIDXEDGE   ) .ne. 0) then
      
      ! Get the twist index array
      if (rtriangulation%h_ItwistIndexEdges .ne. ST_NOHANDLE) then
        call storage_getbase_int (rtriangulation%h_ItwistIndexEdges,p_ItwistIndexEdge)
      else
        call output_line ('Twist indices not available!', &
            OU_CLASS_ERROR,OU_MODE_STD,'elprep_prepareSetForEvaluation')
        call sys_halt()
      end if
      
      ! Allocate memory
      if (.not. associated(revalElementSet%p_ItwistIndexEdge)) then
        allocate(&
            revalElementSet%p_ItwistIndexEdge(nelements))
      end if

      ! Fetch the twist indices from the triangulation.
      do ielidx = 1,nelements
        iel = IelementList(ielIdx)
        revalElementSet%p_ItwistIndexEdge(ielidx) = p_ItwistIndexEdge(iel)
      end do
      
    end if
    
    if (iand(cevaluationTag,EL_EVLTAG_TWISTIDXFACE   ) .ne. 0) then
      ! not yet defined
    end if
    
  contains
  
    subroutine calc_refCoords (DpointsRef,Dpoints,npointsPerElement,nelements)
    
    ! Initialises the cubature point array on the reference element(s).
    
    ! Array with coordinates of points on the reference element.
    ! To be initialised.
    real(DP), dimension(:,:,:), intent(INOUT) :: DpointsRef
    
    ! A set of npointsPerElement tuples (x,y) (or (x,y,z) in 3D) of the points 
    ! where to evaluate. These coordinates define on the reference the coordinates
    ! of the cubature points where to evaluate the element.
    real(DP), dimension(:,:), intent(IN) :: Dpoints
    
    ! Number of points per element
    integer, intent(IN) :: npointsPerElement
    
    ! Number of elements
    integer, intent(IN) :: nelements
    
      ! local variables
      integer :: iel,ipt,idim
    
      ! Copy the coordinates. They are the same for all elements.
      do iel=1,nelements
        do ipt=1,npointsPerElement
          do idim=1,ubound(Dpoints,1)
            DpointsRef(idim,ipt,iel) = Dpoints(idim,ipt)
          end do
        end do
      end do
    
    end subroutine
    
  end subroutine

  !************************************************************************
  
!<subroutine>  

  subroutine elprep_releaseElementSet (revalElementSet)

!<description>
  ! Releases all information in revalElementSet from memory. Cleans up the
  ! structure.
!</description>

!<inputoutput>
  ! The element set that is to be released.
  type(t_evalElementSet), intent(INOUT) :: revalElementSet
!</inputoutput>

! </subroutine>

    ! Release all allocated information
    if (associated(revalElementSet%p_Dcoords)) then
      deallocate(revalElementSet%p_Dcoords)
    end if

    if (associated(revalElementSet%p_DpointsRef) .and. &
        (.not. revalElementSet%bforeignPointsRef)) then
      deallocate(revalElementSet%p_DpointsRef)
    else
      nullify(revalElementSet%p_DpointsRef)
    end if
  
    if (associated(revalElementSet%p_DpointsReal) .and. &
        (.not. revalElementSet%bforeignPointsReal)) then
      deallocate(revalElementSet%p_DpointsReal)
    else
      nullify(revalElementSet%p_DpointsReal)
    end if

    if (associated(revalElementSet%p_Djac)) then
      deallocate(revalElementSet%p_Djac)
    end if
  
    if (associated(revalElementSet%p_Ddetj)) then
      deallocate(revalElementSet%p_Ddetj)
    end if
  
    if (associated(revalElementSet%p_ItwistIndexEdge)) then
      deallocate(revalElementSet%p_ItwistIndexEdge)
    end if
  
    if (associated(revalElementSet%p_ItwistIndexFace)) then
      deallocate(revalElementSet%p_ItwistIndexFace)
    end if

    revalElementSet%npointsPerElement = 0
    revalElementSet%nelements = 0

  end subroutine

  !************************************************************************
  
!<subroutine>  

  subroutine elprep_prepareForEvaluation (revalElement, cevaluationTag, &
      rtriangulation, ielement, ctrafoType, DpointRef, DpointReal, Dcoords)

!<description>
  ! This subroutine prepares a t_evalElement structure to be used for
  ! the evaluation of a finite element in a single point on a single cell.
  ! Dpoint contains the coordinates on the reference element where to
  ! evaluate. This point is mapped onto the element ielement from the 
  ! triangulation rtriangulation.
  ! cevaluationTag specifies an 'evaluation tag' that defines which
  ! information must be prepared by this routine; that tag can be obtained
  ! by asking the finite element what it needs by calling 
  ! elem_getEvaluationTag.
!</description>

!<input>
  ! Evaluation tag. This is a bitfield that specifies which information is
  ! prepared in the element set; it's a combination of EL_EVLTAG_XXXX-constants.
  !
  ! Note: If EL_EVLTAG_REFPOINTS is not specified in this tag, the coordinates on
  ! the reference element are assumed to be initialised! Dpoints is ignored
  ! in that case. This can be used to calculate this information only once
  ! in a first call while then using the same set of reference coordinates
  ! for all subsequent calls.
  integer(I32), intent(IN) :: cevaluationTag

  ! Underlying triangulation of the domain
  type(t_triangulation), intent(IN) :: rtriangulation
  
  ! Element where to evaluate the basis functions
  integer(PREC_ELEMENTIDX), intent(IN) :: ielement
  
  ! Type of transformation from the reference element to the real element.
  integer, intent(IN) :: ctrafoType
  
  ! OPTIONAL: A tuple (x,y) (or (x,y,z) in 3D) of the point where to evaluate.
  ! These coordinates define on the reference the coordinates of the cubature
  ! points where to evaluate the element.
  ! Ignored if EL_EVLTAG_REFPOINTS is not specified in cevaluationTag.
  real(DP), dimension(:), optional :: DpointRef
  
  ! OPTIONAL: A tuple (x,y) (or (x,y,z) in 3D) of the point where to evaluate
  ! on the real element. 
  !
  ! If not specified, the routine will automatically calculate its position
  ! if necessary.
  real(DP), dimension(:), target, optional :: DpointReal
  
  ! OPTIONAL: An array with vertices defining the element.
  !
  ! If not specified, the routine will automatically set up that array
  ! using rtriangulation.
  real(DP), dimension(:,:), target, optional :: Dcoords
!</input>
  
!<inputoutput>
  ! The element set that is to be initialised. If this is already initialised,
  ! previous information is overwritten.
  type(t_evalElement), intent(INOUT) :: revalElement
!</inputoutput>

! </subroutine>

    integer(I32), dimension(:), pointer :: p_ItwistIndexEdge
    
    ! Calculate the transformation
    if (iand(cevaluationTag,EL_EVLTAG_COORDS      ) .ne. 0) then
      if (.not. present(Dcoords)) then
        call trafo_getCoords (ctrafoType,&
              rtriangulation,ielement,revalElement%Dcoords)
      else
        revalElement%Dcoords = Dcoords
      end if
    end if
          
    ! Coordinates on the reference element are always necessary.
    if ((iand(cevaluationTag,EL_EVLTAG_REFPOINTS   ) .ne. 0) .or. &   
        (iand(cevaluationTag,EL_EVLTAG_REALPOINTS   ) .ne. 0)) then
      
      ! Calculate the coordinates on the reference element
      if (iand(cevaluationTag,EL_EVLTAG_REFPOINTS   ) .ne. 0) then
        if (.not. present(DpointRef)) then
          call output_line ('Dpoints not specified!', &
              OU_CLASS_ERROR,OU_MODE_STD,'elprep_prepareSetForEvaluation')
        end if
        revalElement%DpointRef = DpointRef
      end if
    end if
    
    ! Prepare memory for the Jacobian matrix and/or determinant
    if ((iand(cevaluationTag,EL_EVLTAG_JAC   ) .ne. 0) .or. &
        (iand(cevaluationTag,EL_EVLTAG_DETJ  ) .ne. 0)) then
        
      ! If real world coordinates must be calculated, we don't have to calculate
      ! the Jacobian stuff here; it's done below.
      if (iand(cevaluationTag,EL_EVLTAG_REALPOINTS   ) .eq. 0) then
        call trafo_calctrafo (ctrafoType,&
             revalElement%Dcoords,revalElement%DpointRef,revalElement%Djac,revalElement%ddetj)
      end if
      
    end if
    
    ! Calculate the coordinates on the real element
    if (iand(cevaluationTag,EL_EVLTAG_REALPOINTS   ) .ne. 0) then
    
      ! Calculate the real world coordinates from the coordinates on the
      ! reference element. Calculate the Jacobial of the mapping from the
      ! reference element to the real element.
      call trafo_calctrafo (ctrafoType,&
          revalElement%Dcoords,revalElement%DpointRef,revalElement%Djac,revalElement%ddetj,&
          revalElement%DpointReal)
      
    else
    
      ! Take the predefined coordinates -- if specified.
      if (present(DpointReal)) then
        ! Save a pointer to DpointReal
        revalElement%DpointReal = DpointReal
      end if 

    end if

    ! Get the twist indices or the orientation of edges.
    if (iand(cevaluationTag,EL_EVLTAG_TWISTIDXEDGE   ) .ne. 0) then
      
      ! Get the twist index array
      if (rtriangulation%h_ItwistIndexEdges .ne. ST_NOHANDLE) then
        call storage_getbase_int (rtriangulation%h_ItwistIndexEdges,p_ItwistIndexEdge)
      else
        call output_line ('Twist indices not available!', &
            OU_CLASS_ERROR,OU_MODE_STD,'elprep_prepareSetForEvaluation')
        call sys_halt()
      end if
      
      ! Fetch the twist indices from the triangulation.
      revalElement%itwistIndexEdges = p_ItwistIndexEdge(ielement)
      
    end if
    
    if (iand(cevaluationTag,EL_EVLTAG_TWISTIDXFACE   ) .ne. 0) then
      ! not yet defined
    end if
    
  end subroutine
  
end module 
