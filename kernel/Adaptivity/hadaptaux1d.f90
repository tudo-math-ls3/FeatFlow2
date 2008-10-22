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
!#  7.) refine_Line2Line
!#      -> Refines a line by subdivision into two lines
!#
!#
!#  FAQ - Some explanation
!# -----------------------
!# 1.) So, how dies grid refinement wirk in detail?
!#
!#     In 1D, grid refinement is very easy. Each element (i.e. line)
!#     is subdivided into two similar lines, and a new vertex is
!#     inserted at the midpoint of the element. That's it.
!#
!# 2.) 

!# </purpose>
!##############################################################################

module hadaptaux1d

  use collection
  use fsystem
  use hadaptaux
  
  implicit none
  
  private
  public :: hadapt_markRefinement1D
  public :: hadapt_markCoarsening1D
  public :: hadapt_calcNumberOfElements1D
  public :: hadapt_refine1D
!  public :: hadapt_coarsen1D

!<constants>

!<constantblock description="Constants for element marker in 1D">
  
  ! Mark for keeping element 'as is'
  integer, parameter :: MARK_ASIS                   = 0

  ! Mark element for 1-line : 2-lines refinement
  integer, parameter :: MARK_REF_LINE2LINE          = 1

  ! Mark element for recoarsening with right neighbor
  integer, parameter :: MARK_CRS_LINE_RIGHT         = -1

  ! Mark element for recoarsening with left neighbor
  integer, parameter :: MARK_CRS_LINE_LEFT          = -2

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
    type(t_vectorScalar), intent(IN) :: rindicator
!</input>

!<inputoutput>
    ! Adaptive data structure
    type(t_hadapt), intent(INOUT) :: rhadapt
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

    ! Initialize marker structure for NEL0 elements
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

        if (isubdivide .lt. rhadapt%NSUBDIVIDEMAX) then
          
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
    ! what one refinement step can "do". This is in contrast to other "node-
    ! removal" techniques, which remove all superficial vertices and retriangulate
    ! the generated "holes". However, such algorithms cannot guarantee a grid
    ! hierarchy between refined and recoarsened grids.
!</description>

!<input>
    ! Indicator vector for refinement
    type(t_vectorScalar), intent(IN) :: rindicator
!</input>

!<inputoutput>
    ! Adaptive data structure
    type(t_hadapt), intent(INOUT) :: rhadapt
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
      call output_line('Dynamic data structures are not &
                       & generated or no marker for grid refinement exists!',&
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
    
  contains
    
    ! Here, the real working routines follow.
    
    !**************************************************************
    ! For a given element IEL, check if all corners vertices need
    ! to be locked. In particular, this is the case if the element
    
    function lockElement(iel) result (blockElement)
      
      integer, intent(IN) :: iel
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
    type(t_hadapt), intent(IN) :: rhadapt
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
    
    ! Initialize number of elements by current number
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
    type(t_hadapt), intent(INOUT) :: rhadapt

    ! OPTIONAL: Collection
    type(t_collection), intent(INOUT), optional :: rcollection
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
      call output_line('Dynamic data structures are not generated &
                       &or no marker for refinement is available!',&
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
  ! ***************************************************************************
  ! ***************************************************************************

!<subroutine>

  subroutine add_vertex1D(rhadapt, i1, i2, e1, i12,&
                          rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine adds a new vertex at the midpoint of a given line.
!</description>

!<input>
    ! First point of edge on which new vertex will be added
    integer, intent(IN) :: i1

    ! Second point of edge on which new vertex will be added
    integer, intent(IN) :: i2

    ! Number of the right-adjacent element w.r.t. to the oriented edge (I1,I2)
    integer, intent(IN) :: e1

    ! Callback routines
    include 'intf_hadaptcallback.inc'
    optional :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! Adaptive data structure
    type(t_hadapt), intent(INOUT) :: rhadapt

    ! OPTIONAL: Collection
    type(t_collection), intent(INOUT), optional :: rcollection
!</inputoutput>

!<output>
    ! Number of the new vertex located between i1 and i2
    integer, intent(OUT) :: i12
!</output>
!</subroutine>

    ! local variables
    real(DP) :: x1,x2,x12
    integer :: ipos
    integer, dimension(3) :: Ivertices
    integer, dimension(1) :: Ielements

    ! Get coordinates of vertices
    x1 = rhadapt%rVertexCoordinates1D%p_Dkey(i1)
    x2 = rhadapt%rVertexCoordinates1D%p_Dkey(i2)
    
    ! Compute coordinates of new vertex
    x12 = 0.5_DP*(x1+x2)
    
print *, i1,i2,x1,x2,x12,rhadapt%NVT+1

    ! Search for vertec coordinate in binary tree:
    ! If the vertex already exists, e.g., it was added when the 
    ! adjacent element was refined, then nothing needs to be done 
    ! for this vertex
    if (btree_searchInTree(rhadapt%rVertexCoordinates1D, x12, ipos) .eq.&
                           BTREE_FOUND) return

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
    call btree_insertIntoTree(rhadapt%rVertexCoordinates1D, x12)

    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      Ivertices = (/i12, i1, i2/)
      Ielements = (/0/)
      call fcb_hadaptCallback(rcollection, HADAPT_OPR_INSERTVERTEXEDGE,&
                              Ivertices, Ielements)
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
    integer, intent(IN) :: ivt

    
!</input>

!<inputoutput>
    ! Adaptive data structure
    type(t_hadapt), intent(INOUT) :: rhadapt
!</inputoutput>

!<output>
    ! Number of the vertex to replace the deleted one
    integer, intent(OUT) :: ivtReplace
!</output>
!</subroutine>

    ! Check if we are the last 
    
  end subroutine remove_vertex1D

  ! ***************************************************************************

!<subroutine>

  subroutine replace_element1D(rhadapt, ipos, i1, i2, e1, e2)
  
!<description>
    ! This subroutine replaces the vertices and 
    ! adjacent elements for a given element in 1D
!</description>

!<input>
    ! position number of the element in dynamic data structure
    integer, intent(IN) :: ipos

    ! numbers of the element nodes
    integer, intent(IN) :: i1,i2

    ! numbers of the surrounding elements
    integer, intent(IN) :: e1,e2
!</input>

!<inputoutput>
    ! Adaptive data structure
    type(t_hadapt), intent(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    ! Replace triangular element
    rhadapt%p_IverticesAtElement(:,ipos)      = (/i1,i2/)
    rhadapt%p_IneighboursAtElement(:,ipos)    = (/e1,e2/)
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
    integer, intent(IN) :: i1,i2

    ! numbers of the surrounding elements
    integer, intent(IN) :: e1,e2
!</input>

!<inputoutput>
    ! Adaptive data structure
    type(t_hadapt), intent(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>
    
    ! Increase number of elements
    rhadapt%NEL = rhadapt%NEL+1
    rhadapt%InelOfType(TRIA_NVELINE1D) = rhadapt%InelOfType(TRIA_NVELINE1D)+1

    rhadapt%p_IverticesAtElement(:,rhadapt%NEL)      = (/i1,i2/)
    rhadapt%p_IneighboursAtElement(:,rhadapt%NEL)    = (/e1,e2/)
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
    integer, intent(IN) :: iel
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(INOUT) :: rhadapt
!</inputoutput>

!<output>
    ! Former number of the replacement element
    integer, intent(OUT) :: ielReplace
!</output>
!</subroutine>

    ! local variables
    integer :: ipos,ivt,jel,ive,jve

    ! Replace element by the last element and delete last element
    ielReplace = rhadapt%NEL

    ! Check if we are the last element
    if (iel .ne. ielReplace) then
      
      ! Element is not the last one. Then the element that should be removed must
      ! have a smaller element number. If this is not the case, something is wrong.
      if (iel > ielReplace) then
        call output_line('Number of replacement element must not be smaller than that of&
                         & the removed elements!',&
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
        ipos = arrlst_getNextInArraylist(rhadapt%rElementsAtVertex, ivt, .true.)
        elements: do while(ipos .gt. ARRLST_NULL)
          
          ! Check if element number corresponds to the replaced element
          if (rhadapt%rElementsAtVertex%p_IData(ipos) .eq. ielReplace) then
              rhadapt%rElementsAtVertex%p_IData(ipos) = iel
            exit elements
          end if
          
          ! Proceed to next element in list
          ipos = arrlst_getNextInArraylist(rhadapt%rElementsAtVertex, ivt, .false.)
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
    integer, intent(IN) :: jel

    ! Number of the updated macro-element
    integer, intent(IN) :: iel0

    ! Number of the new neighboring element
    integer, intent(IN) :: iel
!</input>

!<inputoutput>
    ! Adaptive data structure
    type(t_hadapt), intent(INOUT) :: rhadapt
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
    !
    !    initial line               subdivided line
    !
    ! (e1)    iel      (e2)      (e1) iel  nel+1  (e2)
    !     +-----------+              +-----+-----+
    !     i1          i2             i1    i3    i2
    !
!</description>

!<input>
    ! Number of element to be refined
    integer, intent(IN) :: iel

    ! Callback function
    include 'intf_hadaptcallback.inc'
    optional :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! adativity structure
    type(t_hadapt), intent(INOUT) :: rhadapt

    ! OPTIONAL: Collection
    type(t_collection), intent(INOUT), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(1) :: Ielements
    integer, dimension(3) :: Ivertices
    integer :: ipos
    integer :: nel0,e1,e2,i1,i2,i3

    ! Store vertex- and element-values of the current element
    i1 = rhadapt%p_IverticesAtElement(1, iel)
    i2 = rhadapt%p_IverticesAtElement(2, iel)

    e1 = rhadapt%p_IneighboursAtElement(1, iel)
    e2 = rhadapt%p_IneighboursAtElement(2, iel)

    ! Store total number of elements before refinement
    nel0 = rhadapt%NEL

    ! Add new vertex I3 at the midpoint of edge (I1,I2)
    call add_vertex1D(rhadapt, i1, i2, e1, i3,&
                      rcollection, fcb_hadaptCallback)

    ! Replace element IEL and add new element NEL0+1
    call replace_element1D(rhadapt, iel, i1, i3, e1, nel0+1)
    call add_element1D(rhadapt, i3, i2, iel, nel0+1)

    ! Update list of neighboring elements
    call update_ElementNeighbors1D(rhadapt, e2, iel, nel0+1)

    ! Update list of elements meeting at vertices
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i2, iel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'refine_Line2Line')
      call sys_halt()
    end if

    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i2, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i3, iel,    ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i3, nel0+1, ipos)
    
    ! Optionally, invoke callback routine
    if (present(fcb_hadaptCallback).and.present(rcollection)) then
      Ivertices = (/i1,i2,i3/)
      Ielements = (/e1/)
      call fcb_hadaptCallback(rcollection, HADAPT_OPR_REF_LINE2LINE,&
                              Ivertices, Ielements)
    end if
  end subroutine refine_Line2Line
  
end module hadaptaux1d
