!##############################################################################
!# ****************************************************************************
!# <name> hadaptaux1d </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# WARNING: Do not USE this module in your applications unless you really
!#          know what you are doing. This module does no error checking!!!
!#
!# This module contains all auxiliary routines which are required for
!# performing h-adaptivity in 1D. Unlike other modules, all subroutines
!# are declared PUBLIC since they are used by module HADAPTIVITY.
!#
!# The following routines are available:
!#
!#  1.) hadapt_markRefinement1D
!#      -> Marks elements for refinement in 1D
!#
!#  2.) hadapt_markCoarsening1D
!#      -> Marks elements for recoarsening in 1D
!#
!#  3.) hadapt_calcNumberOfElements1D
!#      -> Calculates the number of elements after refinement
!#
!#  4.) hadapt_refine1D
!#      -> Performs refinement of marked elements in 1D
!#
!#  5.) hadapt_coarsen1D
!#      -> Performs re-coarsening of marked elements in 1D
!#
!#
!# The following auxiliary routines are available:
!#
!#  1.) add_vertex1D
!#      -> Adds a new vertex to the adaptation data structure in 1D
!#
!#  2.) remove_vertex1D
!#      -> Removes an existing vertex from the adaptation data structure in 1D
!#
!#  3.) replace_element1D
!#      -> Replaces an existing element by another element of he same type in 1D
!#
!#  4.) add_element1D
!#      -> Adds a new element to the adaptation data structure in 1D
!#
!#  5.) remove_element1D
!#      -> Removes an existing element from the adaptation data structure in 1D
!#
!#  6.) update_ElementNeighbors1D
!#      -> Updates the list of neighboring elements in 1D
!#
!#  7.) update_AllElementNeighbors1D
!#      -> Updates the lists of neighboring elements of ALL adjacent elements in 1D
!#
!#  8.) refine_Line2Line
!#      -> Refines a line by subdivision into two lines
!#
!#  9.) coarsen_2Line1Line
!#      -> Coarsens two neighbouring lines into a single line
!#
!#
!#  FAQ - Some explanation \\
!# ----------------------- \\
!# 1.) So, how does grid refinement wirk in detail?
!#
!#     In 1D, grid refinement is very easy. Each element (i.e. line)
!#     is subdivided into two similar lines, and a new vertex is
!#     inserted at the midpoint of the element. That is it.
!#
!# 2.)
!# </purpose>
!##############################################################################

module hadaptaux1d

!$use omp_lib
  use arraylistInt
  use collection
  use fsystem
  use genoutput
  use hadaptaux
  use linearsystemscalar
  use storage
  use triangulation

  implicit none

  private
  public :: hadapt_markRefinement1D
  public :: hadapt_markCoarsening1D
  public :: hadapt_calcNumberOfElements1D
  public :: hadapt_refine1D
  public :: hadapt_coarsen1D

!<constants>

!<constantblock description="Constants for element marker in 1D">

  ! Mark for keeping element 'as is'
  integer, parameter, public :: MARK_ASIS                   = 0

  ! Mark element for 1-line : 2-lines refinement
  integer, parameter, public :: MARK_REF_LINE2LINE          = 1

  ! Mark element for recoarsening with right neighbor
  integer, parameter, public :: MARK_CRS_LINE_RIGHT         = -1

  ! Mark element for recoarsening with left neighbor
  integer, parameter, public :: MARK_CRS_LINE_LEFT          = -2

!</constantblock>

!</constants>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_markRefinement1D(rhadapt, rindicator)

!<description>
    ! This subroutine marks all elements that should be refined
    ! due to accuracy reasons. The decision is based on some
    ! indicator vector which must be given element-wise.
!</description>

!<input>
    ! Indicator vector for refinement
    type(t_vectorScalar), intent(in) :: rindicator
!</input>

!<inputoutput>
    ! Adaptive data structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dindicator
    integer, dimension(:),  pointer :: p_Imarker
    integer, dimension(TRIA_MAXNVE1D) :: IverticesAtElement
    integer  :: ivt,iel,ive,isubdivide

    ! Check if dynamic data structures are generated and contain data
    if (iand(rhadapt%iSpec, HADAPT_HAS_DYNAMICDATA1D) .ne.&
                            HADAPT_HAS_DYNAMICDATA1D) then
      call output_line('Dynamic data structures are not generated!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_markRefinement1D')
      call sys_halt()
    end if

    ! Initialise marker structure for NEL0 elements
    if (rhadapt%h_Imarker .ne. ST_NOHANDLE) call storage_free(rhadapt%h_Imarker)
    call storage_new('hadapt_markRefinement1D', 'Imarker', rhadapt%NEL0,&
                     ST_INT, rhadapt%h_Imarker, ST_NEWBLOCK_ZERO)


    ! Set pointers
    call lsyssc_getbase_double(rindicator, p_Dindicator)
    call storage_getbase_int(rhadapt%h_Imarker, p_Imarker)

    ! Set state of all vertices to "free". Note that vertices of the
    ! initial triangulation are always "locked", i.e. have no positive age.
    do ivt = 1, size(rhadapt%p_IvertexAge, 1)
      rhadapt%p_IvertexAge(ivt) = abs(rhadapt%p_IvertexAge(ivt))
    end do

    ! Loop over all elements and mark those for which the
    ! indicator is greater than the prescribed treshold
    mark_elem: do iel = 1, rhadapt%NEL

      ! Check if element indicator exceeds tolerance
      if (p_Dindicator(iel) .gt. rhadapt%drefinementTolerance) then

        ! Retrieve local data
        IverticesAtElement(1:TRIA_NVELINE1D) =&
            rhadapt%p_IverticesAtElement(1:TRIA_NVELINE1D, iel)

        ! An element can only be refined, if all of its vertices do
        ! not exceed the number of admissible subdivision steps.
        isubdivide = maxval(abs(rhadapt%p_IvertexAge(&
                                IverticesAtElement(1:TRIA_NVELINE1D))))

        if (isubdivide .lt. rhadapt%nsubdividemax) then

          ! Mark element IEL for subdivision
          p_Imarker(iel) = MARK_REF_LINE2LINE

          ! Increase number of vertices by one
          rhadapt%increaseNVT = rhadapt%increaseNVT+1

          ! "Lock" all vertices connected to element IEL since they cannot be removed.
          do ive = 1, TRIA_NVELINE1D
            ivt = IverticesAtElement(ive)
            rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
          end do

        else

          ! Unmark element IEL for refinement
          p_Imarker(iel) = MARK_ASIS

          ! According to the indicator, this element should be refined. Since the
          ! maximum admissible refinement level has been reached no refinement
          ! was performed. At the same time, all vertices of the element should
          ! be "locked" to prevent this element from coarsening
          do ive = 1, TRIA_NVELINE1D
            ivt = IverticesAtElement(ive)
            rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
          end do

        end if

      else   ! p_Dindicator(IEL) <= drefinementTolerance

        ! Unmark element for refinement
        p_Imarker(iel) = MARK_ASIS

      end if
    end do mark_elem

    ! Set specifier to "marked for refinement"
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_MARKEDREFINE)

  end subroutine hadapt_markRefinement1D

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_markCoarsening1D(rhadapt, rindicator)

!<description>
    ! This routine marks all elements that sould be recoarsened due to accuracy
    ! reasons. The decision is based on some indicator vector which must be
    ! given element-wise. The recoarsening strategy is as follows: A subset of
    ! elements can only be coarsened if this subset results from a previous
    ! refinement step. In other words, the grid cannot become coarser than the
    ! initial triangulation. Moreover, one recoarsening step can only "undo"
    ! what one refinement step can "do". This is in contrast to other
    ! "node-removal" techniques, which remove all superficial vertices and
    ! retriangulate the generated "holes". However, such algorithms cannot
    ! guarantee a grid hierarchy between refined and recoarsened grids.
!</description>

!<input>
    ! Indicator vector for refinement
    type(t_vectorScalar), intent(in) :: rindicator
!</input>

!<inputoutput>
    ! Adaptive data structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dindicator
    integer, dimension(:), pointer :: p_Imarker
    integer :: ivt,jvt,iel,ive

    ! Check if dynamic data structures are generated and contain data
    if (iand(rhadapt%iSpec, HADAPT_HAS_DYNAMICDATA1D) .ne.&
                            HADAPT_HAS_DYNAMICDATA1D .or.&
        iand(rhadapt%iSpec, HADAPT_MARKEDREFINE) .ne.&
                            HADAPT_MARKEDREFINE) then
      call output_line('Dynamic data structures are not ' // &
                       'generated or no marker for grid refinement exists!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_markCoarsening1D')
      call sys_halt()
    end if

    ! Set pointers
    call storage_getbase_int(rhadapt%h_Imarker, p_Imarker)
    call lsyssc_getbase_double(rindicator, p_Dindicator)

    ! All nodes of the initial triangulation have age-0, and hence,
    ! they will never be deleted by the re-coarsening procedure.
    !
    ! Phase 1: Lock vertices which cannot be deleted at first sight.
    !
    ! In the first phase, all vertices which do belong to some
    ! macro-element are "locked" since they cannot be deleted.
    !
    ! For each element which is marked for refinement,
    ! the corner vertices are "locked" unconditionally.
    !
    ! For each element that is not marked for refinement, the
    ! indicator is considered and all vertices are "locked" if
    ! re-coarsening is not allowed due to accuracy reasons.
    do iel = 1, rhadapt%NEL

      ! Check if the current element is marked for refinement.
      if ((p_Imarker(iel) .eq. MARK_REF_LINE2LINE) .or.&
          lockElement(iel)) then

        ! The current element is marked for refinement
        ! or it should be locked to to accuracy reasons
        do ive = 1, TRIA_NVELINE1D
          ! "Lock" vertex at corner
          ivt = rhadapt%p_IverticesAtElement(ive, iel)
          rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
        end do

      else   ! no need to keep element

        ! Get global vertex numbers
        ivt = rhadapt%p_IverticesAtElement(1, iel)
        jvt = rhadapt%p_IverticesAtElement(2, iel)

        ! "Lock" the older endpoint
        if (abs(rhadapt%p_IvertexAge(ivt)) .lt.&
            abs(rhadapt%p_IvertexAge(jvt))) then
          rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
        else
          rhadapt%p_IvertexAge(jvt) = -abs(rhadapt%p_IvertexAge(jvt))
        end if

      end if
    end do

    ! Loop over all elements and determine those elements which
    ! can be combined together with its neighbouring cell
    do iel = 1, rhadapt%NEL

      ! Get global vertex numbers
      ivt = rhadapt%p_IverticesAtElement(1, iel)
      jvt = rhadapt%p_IverticesAtElement(2, iel)

      ! Check if left vertex is "free"
      if (rhadapt%p_IvertexAge(ivt) .gt. 0) then

        ! Set marker for left neighbour
        p_Imarker(iel) = MARK_CRS_LINE_LEFT

      ! Check if right vertex is "free"
      elseif (rhadapt%p_IvertexAge(jvt) .gt. 0) then

        ! Set marker for right neighbour
        p_Imarker(iel) = MARK_CRS_LINE_RIGHT
      end if
    end do

    ! Set specifier to "marked for coarsening"
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_MARKEDCOARSEN)

  contains

    ! Here, the real working routines follow.

    !**************************************************************
    ! For a given element IEL, check if all corners vertices need
    ! to be locked. In particular, this is the case if the element

    function lockElement(iel) result (blockElement)

      integer, intent(in) :: iel
      logical :: blockElement

      ! Check if element marker is available for this element
      blockElement = (iel .le. size(p_Dindicator, 1))

      if (blockElement) then
        blockElement = (p_Dindicator(iel) .ge. rhadapt%dcoarseningTolerance)
      end if
    end function lockElement
  end subroutine hadapt_markCoarsening1D

  ! ***************************************************************************

!<function>

  function hadapt_CalcNumberOfElements1D(rhadapt) result(nel)

!<description>
    ! This function calculates the number of elements present
    ! in the triangulation after refinement has been performed
!</description>

!<input>
    ! adaptivity structure
    type(t_hadapt), intent(in) :: rhadapt
!</input>

!<result>
    ! number of elements after refinement
    integer :: nel
!</result>
!</function>

    ! local variables
    integer, dimension(:), pointer :: p_Imarker
    integer :: iel

    call storage_getbase_int(rhadapt%h_Imarker, p_Imarker)

    ! Initialise number of elements by current number
    nel = rhadapt%NEL0

    ! Loop over all elements and check marker
    do iel = 1, rhadapt%NEL0
      select case(p_Imarker(iel))
      case(MARK_REF_LINE2LINE)
        nel = nel+1
      end select
    end do
  end function hadapt_CalcNumberOfElements1D

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_refine1D(rhadapt, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine refinement the elements according to the marker
!</description>

!<input>
    ! Callback function
    include 'intf_hadaptcallback.inc'
    optional :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! adaptive data structure
    type(t_hadapt), intent(inout) :: rhadapt

    ! OPTIONAL: Collection
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Imarker
    integer :: iel

    ! Check if dynamic data structures are o.k. and if
    ! cells are marked for refinement
    if (iand(rhadapt%iSpec, HADAPT_HAS_DYNAMICDATA1D) .ne.&
                            HADAPT_HAS_DYNAMICDATA1D .or.&
        iand(rhadapt%iSpec, HADAPT_MARKEDREFINE) .ne.&
                            HADAPT_MARKEDREFINE) then
      call output_line('Dynamic data structures are not generated ' // &
                       'or no marker for refinement is available!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_refine1D')
      call sys_halt()
    end if

    ! Set pointers
    call storage_getbase_int(rhadapt%h_Imarker, p_Imarker)

    ! Perform red-green refinement
    do iel = 1, size(p_Imarker, 1)

      select case(p_Imarker(iel))
      case(MARK_ASIS,&
           MARK_CRS_LINE_RIGHT,&
           MARK_CRS_LINE_LEFT)

        ! Do nothing for elements that should be kept 'as is'
        ! and those which are marked for re-coarsening

      case (MARK_REF_LINE2LINE)

        ! Regular refinement
        call refine_Line2Line(rhadapt, iel,&
                              rcollection, fcb_hadaptCallback)
        p_Imarker(iel) = MARK_ASIS

      case DEFAULT
        call output_line('Invalid element refinement marker!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_refine1D')
        call sys_halt()
      end select
    end do

    ! Increase the number of refinement steps by one
    rhadapt%nRefinementSteps = rhadapt%nRefinementSteps+1

    ! The markers are no longer valid
    rhadapt%iSpec = iand(rhadapt%iSpec, not(HADAPT_MARKEDREFINE))
  end subroutine hadapt_refine1D

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_coarsen1D(rhadapt, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine coarsens the elements according to the marker
!</description>

!<input>
    ! callback routines
    include 'intf_hadaptcallback.inc'
    optional :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! adativity structure
    type(t_hadapt), intent(inout) :: rhadapt

    ! OPTIONAL Collection
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(it_arraylistInt) :: ralstIter
    integer,  dimension(:), pointer :: p_Imarker
    integer :: iel,ive,ivt,ivtReplace,jel

    ! Check if dynamic data structures are o.k. and if
    ! cells are marked for refinement
    if (iand(rhadapt%iSpec, HADAPT_HAS_DYNAMICDATA1D) .ne.&
                            HADAPT_HAS_DYNAMICDATA1D .or.&
        iand(rhadapt%iSpec, HADAPT_MARKEDCOARSEN) .ne.&
                            HADAPT_MARKEDCOARSEN) then
      call output_line('Dynamic data structures are not generated ' // &
                       'or no marker for coarsening is available!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_coarsen1D')
      call sys_halt()
    end if

    ! Set pointers
    call storage_getbase_int(rhadapt%h_Imarker, p_Imarker)

    ! Perform hierarchical recoarsening
    element: do iel = size(p_Imarker,1), 1, -1
      select case(p_Imarker(iel))
      case(MARK_ASIS:)

        ! Do nothing for elements that should be kept 'as is'
        ! and those which are marked for refinement.

      case(MARK_CRS_LINE_LEFT)
        ! Do nothing since this element will be processed
        ! below initiated by the right element neighbour

      case(MARK_CRS_LINE_RIGHT)
        call coarsen_2Line1Line(rhadapt, iel, p_Imarker(iel),&
                                rcollection, fcb_hadaptCallback)

      case DEFAULT
        call output_line('Invalid recoarsening marker!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_coarsen1D')
        call sys_halt()
      end select
    end do element


    ! Loop over all vertices 1...NVT0 present in the triangulation
    ! before refinement and check if they are free for vertex removal.
    vertex: do ivt = rhadapt%NVT0, 1, -1

      ! If the vertex is locked, then skip this vertex
      if (rhadapt%p_IvertexAge(ivt) .le. 0) cycle vertex

      ! Remove vertex physically. Note that this vertex is no longer
      ! associated to any element. All associations have been removed
      ! in the above element coarsening/conversion step. In order to
      ! prevent "holes" in the vertex list, vertex IVT is replaced by
      ! the last vertex if it is not the last one itself.
      call remove_vertex1D(rhadapt, ivt, ivtReplace)

      ! If vertex IVT was not the last one, update the
      ! "elements-meeting-at-vertex" list
      if (ivtReplace .ne. 0) then

        ! Start with first element in "elements-meeting-at-vertex"
        ! list of the replaced vertex
        ralstIter = alst_begin(rhadapt%rElementsAtVertex, ivtReplace)

        update: do while(.not.alst_isNull(ralstIter))

          ! Get element number JEL from arraylist
          jel = alst_get(rhadapt%rElementsAtVertex, ralstIter)

          ! Proceed to next element
          call alst_next(ralstIter)

          ! Look for vertex ivtReplace in element JEL and replace it by IVT
          do ive = 1, TRIA_NVELINE1D
            if (rhadapt%p_IverticesAtElement(ive, jel) .eq. ivtReplace) then
              rhadapt%p_IverticesAtElement(ive, jel) = ivt
              cycle update
            end if
          end do

          ! If the replaced vertex ivtReplace could not be found in element JEL
          ! something is wrong and we stop the simulation
          call output_line('Unable to find replacement vertex in element',&
                           OU_CLASS_ERROR,OU_MODE_STD,'hadapt_coarsen1D')
          call sys_halt()
        end do update

        ! Swap tables IVT and ivtReplace in arraylist and release table ivtReplace
        call alst_swapTbl(rhadapt%rElementsAtVertex, ivt, ivtReplace)
        call alst_releaseTbl(rhadapt%rElementsAtVertex, ivtReplace)

      else

        ! Release table IVT
        call alst_releaseTbl(rhadapt%rElementsAtVertex, ivt)
      end if

      ! Optionally, invoke callback function
      if (present(fcb_hadaptCallback) .and. present(rcollection)) then
        rcollection%IquickAccess(1:2) = (/ivt, ivtReplace/)
        call fcb_hadaptCallback(HADAPT_OPR_REMOVEVERTEX, rcollection)
      end if
    end do vertex

    ! Increase the number of recoarsening steps by one
    rhadapt%nCoarseningSteps = rhadapt%nCoarseningSteps+1

    ! The markers are no longer valid
    rhadapt%iSpec = iand(rhadapt%iSpec, not(HADAPT_MARKEDCOARSEN))
  end subroutine hadapt_coarsen1D

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<subroutine>

  subroutine add_vertex1D(rhadapt, i1, i2, i12,&
                          rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine adds a new vertex at the midpoint of a given line.
!</description>

!<input>
    ! First point of edge on which new vertex will be added
    integer, intent(in) :: i1

    ! Second point of edge on which new vertex will be added
    integer, intent(in) :: i2

    ! Callback routines
    include 'intf_hadaptcallback.inc'
    optional :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! Adaptive data structure
    type(t_hadapt), intent(inout) :: rhadapt

    ! OPTIONAL: Collection
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Number of the new vertex located between i1 and i2
    integer, intent(out) :: i12
!</output>
!</subroutine>

    ! local variables
    real(DP) :: x1,x2,x12
    integer :: ivt

    ! Get coordinates of vertices
    x1 = rhadapt%p_DvertexCoords1D(1,i1)
    x2 = rhadapt%p_DvertexCoords1D(1,i2)

    ! Compute coordinates of new vertex
    x12 = 0.5_DP*(x1+x2)

    ! Search for vertex coordinate: Stupid linear search.
    do ivt = 1, rhadapt%NVT
      if (abs(x12 - rhadapt%p_DvertexCoords1D(1,ivt)) .lt. SYS_EPSREAL_DP) return
    end do

    ! Otherwise, update number of vertices
    rhadapt%NVT = rhadapt%NVT+1
    i12         = rhadapt%NVT

    ! Set age of vertex
    rhadapt%p_IvertexAge(i12) = &
        1+max(abs(rhadapt%p_IvertexAge(i1)),&
              abs(rhadapt%p_IvertexAge(i2)))

    ! Set nodal property
    rhadapt%p_InodalProperty(i12) = 0

    ! Add new entry to vertex coordinates
    rhadapt%p_DvertexCoords1D(1,i12) = x12

    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:3) = (/i12, i1, i2/)
      call fcb_hadaptCallback(HADAPT_OPR_INSERTVERTEXEDGE, rcollection)
    end if
  end subroutine add_vertex1D

  ! ***************************************************************************

!<subroutine>

  subroutine remove_vertex1D(rhadapt, ivt, ivtReplace)

!<description>
    ! This subroutine removes an existing vertex from the adaptivity structure
    ! and moves the last vertex at its position. The number of the replacement
    ! vertex is returned as ivtReplace. If the vertex to be replace is the last
    ! vertex then ivtReplace=0 is returned on output.
!</description>

!<input>
    ! Number of the vertex to be deleted
    integer, intent(in) :: ivt


!</input>

!<inputoutput>
    ! Adaptive data structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>

!<output>
    ! Number of the vertex to replace the deleted one
    integer, intent(out) :: ivtReplace
!</output>
!</subroutine>

    ! Check if we are not the last vertex
    if (ivt .lt. rhadapt%NVT) then

      ! Get last vertex of the adaptivity structure
      ivtReplace = rhadapt%NVT

      ! Copy data from last position to current position
      rhadapt%p_InodalProperty(ivt)    = rhadapt%p_InodalProperty(ivtReplace)
      rhadapt%p_IvertexAge(ivt)        = rhadapt%p_IvertexAge(ivtReplace)
      rhadapt%p_DvertexCoords1D(1,ivt) = rhadapt%p_DvertexCoords1D(1,ivtReplace)

      ! Clear data for node IVTREPLACE
      rhadapt%p_InodalProperty(ivtReplace)    = 0
      rhadapt%p_IvertexAge(ivtReplace)        = 0
      rhadapt%p_DvertexCoords1D(1,ivtReplace) = 0.0_DP

    else

      ! IVT is the last vertex of the adaptivity structure
      ivtReplace = 0

      ! Clear data for node IVT
      rhadapt%p_InodalProperty(ivt)    = 0
      rhadapt%p_IvertexAge(ivt)        = 0
      rhadapt%p_DvertexCoords1D(1,ivt) = 0.0_DP

    end if

    ! Decrease number of vertices by one
    rhadapt%NVT = rhadapt%NVT-1

  end subroutine remove_vertex1D

  ! ***************************************************************************

!<subroutine>

  subroutine replace_element1D(rhadapt, iel, i1, i2, e1, e2)

!<description>
    ! This subroutine replaces the vertices and
    ! adjacent elements for a given element in 1D
!</description>

!<input>
    ! number of the element
    integer, intent(in) :: iel

    ! numbers of the element nodes
    integer, intent(in) :: i1,i2

    ! numbers of the surrounding elements
    integer, intent(in) :: e1,e2
!</input>

!<inputoutput>
    ! Adaptive data structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    ! Replace triangular element
    rhadapt%p_IverticesAtElement(1:2,iel)   = (/i1,i2/)
    rhadapt%p_IneighboursAtElement(1:2,iel) = (/e1,e2/)
  end subroutine replace_element1D

  ! ***************************************************************************

!<subroutine>

  subroutine add_element1D(rhadapt, i1, i2, e1, e2)

!<description>
    ! This subroutine adds a new element connected to
    ! two vertices and surrounded by two adjacent elements
!</description>

!<input>
    ! numbers of the element nodes
    integer, intent(in) :: i1,i2

    ! numbers of the surrounding elements
    integer, intent(in) :: e1,e2
!</input>

!<inputoutput>
    ! Adaptive data structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    ! Increase number of elements
    rhadapt%NEL = rhadapt%NEL+1
    rhadapt%InelOfType(TRIA_NVELINE1D) = rhadapt%InelOfType(TRIA_NVELINE1D)+1

    rhadapt%p_IverticesAtElement(1:2,rhadapt%NEL)   = (/i1,i2/)
    rhadapt%p_IneighboursAtElement(1:2,rhadapt%NEL) = (/e1,e2/)
  end subroutine add_element1D

! ***************************************************************************

!<subroutine>

  subroutine remove_element1D(rhadapt, iel, ielReplace)

!<description>
    ! This subroutine removes an existing element and moves the last
    ! element of the adaptation data structure to its position.
    ! The routine returns the former number ielReplace of the last
    ! element. If iel is the last element, then ielReplace=0 is returned.
!</description>

!<input>
    ! Number of the element that should be removed
    integer, intent(in) :: iel
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>

!<output>
    ! Former number of the replacement element
    integer, intent(out) :: ielReplace
!</output>
!</subroutine>

    ! local variables
    type(it_arrayListInt) :: ralstIter
    integer, pointer :: p_iel
    integer :: ivt,jel,ive,jve

    ! Replace element by the last element and delete last element
    ielReplace = rhadapt%NEL

    ! Check if we are the last element
    if (iel .ne. ielReplace) then

      ! Element is not the last one. Then the element that should be removed must
      ! have a smaller element number. If this is not the case, something is wrong.
      if (iel > ielReplace) then
        call output_line('Number of replacement element must not be smaller than that of' // &
                         ' the removed elements!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'remove_element1D')
        call sys_halt()
      end if

      ! The element which formally was labeled ielReplace is now labeled IEL.
      ! This modification must be updated in the list of adjacent element
      ! neighbors of all surrounding elements. Moreover, the modified element
      ! number must be updated in the "elements-meeting-at-vertex" lists of the
      ! corner nodes of element IEL. Both operations are performed below.
      update: do ive = 1, TRIA_NVELINE1D

        ! Get vertex number of corner node
        ivt = rhadapt%p_IverticesAtElement(ive, ielReplace)

        ! Start with first element in "elements-meeting-at-vertex" list
        ralstIter = alst_begin(rhadapt%rElementsAtVertex, ivt)
        elements: do while(.not.alst_isNull(ralstIter))

          ! Check if element number corresponds to the replaced element
          call alst_getbase_key(rhadapt%rElementsAtVertex, ralstIter, p_iel)
          if (p_iel .eq. ielReplace) then
            p_iel = iel
            exit elements
          end if

          ! Proceed to next element in list
          call alst_next(ralstIter)
        end do elements


        ! Get element number of element JEL which
        ! is adjacent to element ielReplace
        jel = rhadapt%p_IneighboursAtElement(ive, ielReplace)

        ! If we are at the boundary or if element IEL to be removed is
        ! adjacent to element ielReplace, then skip this edge.
        if (jel .eq. 0 .or. jel .eq. iel) cycle update

        adjacent: do jve = 1, TRIA_NVELINE1D
          if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. ielReplace) then
              rhadapt%p_IneighboursAtElement(jve, jel) = iel
            exit adjacent
          end if
        end do adjacent

      end do update

      ! Copy data from element ielReplace to element IEL
      rhadapt%p_IverticesAtElement(:,iel) =&
          rhadapt%p_IverticesAtElement(:,ielReplace)
      rhadapt%p_IneighboursAtElement(:,iel) =&
          rhadapt%p_IneighboursAtElement(:,ielReplace)

    else

      ! Element iel is the last element
      ielReplace = 0
    end if

    ! Decrease number of elements
    rhadapt%NEL = rhadapt%NEL-1
  end subroutine remove_element1D

  ! ***************************************************************************

!<subroutine>

  subroutine update_ElementNeighbors1D(rhadapt, jel, iel0, iel)

!<description>
    ! This subroutine updates the list of elements
    ! adjacent to another element
!</description>

!<input>
    ! Number of the neighboring element
    integer, intent(in) :: jel

    ! Number of the updated macro-element
    integer, intent(in) :: iel0

    ! Number of the new neighboring element
    integer, intent(in) :: iel
!</input>

!<inputoutput>
    ! Adaptive data structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ive

    ! Do nothing for elements adjacent to the boundary
    if (jel .eq. 0) return

    do ive = 1, TRIA_NVELINE1D
      if (rhadapt%p_IneighboursAtElement(ive, jel) .eq. iel0) then
          rhadapt%p_IneighboursAtElement(ive, jel) = iel
        return
      end if
    end do

    ! If we end up here, the neighboring element could not be found
    call output_line('Inconsistent adjacency lists!',&
                     OU_CLASS_ERROR,OU_MODE_STD,'update_ElementNeighbors1D')
    call sys_halt()
  end subroutine update_ElementNeighbors1D

  ! ***************************************************************************

!<subroutine>

  subroutine update_AllElementNeighbors1D(rhadapt, iel0, iel)

!<description>
    ! This subroutine updates the list of elements adjacent to another elements.
    ! For all elements jel which are adjacent to the old element iel0 the new
    ! value iel is stored in the neighbours-at-element structure.
!</description>

!<input>
    ! Number of the element to be updated
    integer, intent(in) :: iel0

    ! New value of the element
    integer, intent(in) :: iel
!</input>

!<inputoutput>
    ! Adaptive data structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ive,jel,jve

    ! Check if the old element is still present in the triangulation
    if (iel0 .gt. rhadapt%NEL) return

    ! Loop over adjacent elements
    adjacent: do ive = 1, TRIA_NVELINE1D
      ! Get number of adjacent element
      jel = rhadapt%p_IneighboursAtElement(ive, iel0)

      ! Are we at the boundary?
      if (jel .eq. 0) cycle adjacent

      ! Find position of element IEL0 in adjacent element JEL
      do jve = 1, TRIA_NVELINE1D
        if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. iel0) then
          rhadapt%p_IneighboursAtElement(jve, jel) = iel
          cycle adjacent
        end if
      end do

      ! If the old element number was not found in adjacent element JEL
      ! then something went wrong and we should not proceed.
      call output_line('Inconsistent adjacency lists!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'update_AllElementNeighbors1D')
      call sys_halt()
    end do adjacent
  end subroutine update_AllElementNeighbors1D

  ! ***************************************************************************

!<subroutine>

  subroutine refine_Line2Line(rhadapt, iel, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one line into two lines
    ! by subdividion the element at the midpoint.
    !
    ! In the illustration, i1, i2 and i3 denote the vertices whereby i3
    ! is the new vertex. Moreover, (e1) and (e2) stand for the element numbers
    ! which are adjacent to the current elements. The new element is assigned
    ! the total number of elements currently present in the triangulation
    ! increased by one.
    ! <verb>
    !    initial line               subdivided line
    !
    ! (e1)    iel      (e2)      (e1) iel   nel+1 (e2)
    !     +-----------+              +-----+-----+
    !     i1          i2             i1    i3    i2
    ! </verb>
!</description>

!<input>
    ! Number of element to be refined
    integer, intent(in) :: iel

    ! Callback function
    include 'intf_hadaptcallback.inc'
    optional :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! adativity structure
    type(t_hadapt), intent(inout) :: rhadapt

    ! OPTIONAL: Collection
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(it_arrayListInt) :: ralstIter
    integer :: e1,e2,i1,i2,i3,nel0

    ! Store vertex- and element-values of the current element
    i1 = rhadapt%p_IverticesAtElement(1, iel)
    i2 = rhadapt%p_IverticesAtElement(2, iel)

    e1 = rhadapt%p_IneighboursAtElement(1, iel)
    e2 = rhadapt%p_IneighboursAtElement(2, iel)

    ! Store total number of elements before refinement
    nel0 = rhadapt%NEL

    ! Add new vertex I3 at the midpoint of edge (I1,I2)
    call add_vertex1D(rhadapt, i1, i2, i3,&
                      rcollection, fcb_hadaptCallback)

    ! Replace element IEL and add new element NEL0+1
    call replace_element1D(rhadapt, iel, i1, i3, e1, nel0+1)
    call add_element1D(rhadapt, i3, i2, iel, e2)

    ! Update list of neighboring elements
    call update_ElementNeighbors1D(rhadapt, e2, iel, nel0+1)

    ! Update list of elements meeting at vertices
    ralstIter = alst_find(rhadapt%relementsAtVertex, i2, iel)
    if (alst_isNull(ralstIter)) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'refine_Line2Line')
      call sys_halt()
    else
      ralstIter = alst_erase(rhadapt%relementsAtVertex, ralstIter)
    end if

    call alst_push_back(rhadapt%relementsAtVertex, i2, nel0+1)
    call alst_push_back(rhadapt%relementsAtVertex, i3, iel)
    call alst_push_back(rhadapt%relementsAtVertex, i3, nel0+1)

    ! Optionally, invoke callback routine
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:3) = (/i1, i2, i3/)
      call fcb_hadaptCallback(HADAPT_OPR_REF_LINE2LINE, rcollection)
    end if
  end subroutine refine_Line2Line

  ! ***************************************************************************

!<subroutine>

  subroutine coarsen_2Line1Line(rhadapt, iel, imarker,&
                                rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine combines two lines resulting from a 1-line : 2-line
    ! refinement into the original macro line. Depending on the value of
    ! IMARKER the given element IEL is combined with its left or right
    ! neighbouring element.
    !
    ! In the illustration, i1, i2 and i3 denote the vertices whereby i3
    ! is removed. Moreover, (e1) and (e2) stand for the element numbers
    ! which are adjacent to the element iel and iel1, respectively. The
    ! element with the larger element number JEL=MIN(IEL,IEL1) is preserved
    ! whereas the one with the larger element number is removed. The
    ! total number of elements in the triangulation is decreased by one.
    !
    !    initial lines              combined lines
    !
    ! (e1) iel   iel1  (e2)      (e1)     jel     (e2)
    !     +-----+-----+              +-----------+
    !     i1    i3    i2             i1          i2
    !
!</description>

!<input>
    ! Element number of one line
    integer, intent(in) :: iel

    ! Identifier for element marker
    integer, intent(in) :: imarker

    ! callback routines
    include 'intf_hadaptcallback.inc'
    optional :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt

    ! OPTIONAL: Collection
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(it_arraylistInt) :: ralstIter
    integer :: e1,e2,i1,i2,i3,iel1,ielRemove,ielReplace,jel

    select case(imarker)
    case (MARK_CRS_LINE_RIGHT)
      ! Get right-adjacent element number
      iel1 = rhadapt%p_IneighboursAtElement(2,iel)

      ! Store vertex- and element values of the two elements
      i1 = rhadapt%p_IverticesAtElement(1, iel)
      i3 = rhadapt%p_IverticesAtElement(2, iel)
      i2 = rhadapt%p_IverticesAtElement(2, iel1)

      e1 = rhadapt%p_IneighboursAtElement(1, iel)
      e2 = rhadapt%p_IneighboursAtElement(2, iel1)

      ! Which element has smaller element number?
      if (iel .lt. iel1) then
        ! Keep IEL, remove, IEL1
        jel=iel; ielRemove=iel1

        ! Update list of neighboring elements
        call update_ElementNeighbors1D(rhadapt, e2, iel1, iel)

      else
        ! Keep IEL1, remove IEL
        jel=iel1; ielRemove=iel

        ! Update list of neighboring elements
        call update_ElementNeighbors1D(rhadapt, e1, iel, iel1)
      end if

      ! Replace the element JEL
      call replace_element1D(rhadapt, jel, i1, i2, e1, e2)

      ! Remove the other element
      call remove_element1D(rhadapt, ielRemove, ielReplace)
      if (ielReplace .ne. 0)&
          call update_AllElementNeighbors1D(rhadapt, ielReplace, ielRemove)

      ! Update list of elements meeting at vertices
      if (ielRemove .eq. iel1) then

        ! Remove element IELREMOVE from right endpoint
        ralstIter = alst_find(rhadapt%relementsAtVertex, i2, ielRemove)
        if (alst_isNull(ralstIter)) then
          call output_line('Unable to delete element from vertex list!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Line1Line')
          call sys_halt()
        else
          ralstIter = alst_erase(rhadapt%relementsAtVertex, ralstIter)
        end if

        ! Add new element JEL to right endpoint
        call alst_push_back(rhadapt%relementsAtVertex, i2, jel)
      else

        ! Remove element IELREMOVE from left endpoint
        ralstIter = alst_find(rhadapt%relementsAtVertex, i1, ielRemove)
        if (alst_isNull(ralstIter)) then
          call output_line('Unable to delete element from vertex list!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Line1Line')
          call sys_halt()
        else
          ralstIter = alst_erase(rhadapt%relementsAtVertex, ralstIter)
        end if

        ! Add new element JEL to left endpoint
        call alst_push_back(rhadapt%relementsAtVertex, i1, jel)
      end if

      ! Optionally, invoke callback function
      if (present(fcb_hadaptCallback) .and. present(rcollection)) then
        rcollection%IquickAccess(1:3) = (/i1, i2, i3/)
        call fcb_hadaptCallback(HADAPT_OPR_CRS_2LINE1LINE, rcollection)
      end if

    case DEFAULT
      call output_line('Invalid element coarsening marker!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Line1Line')
      call sys_halt()
    end select
  end subroutine coarsen_2Line1Line

end module hadaptaux1d
