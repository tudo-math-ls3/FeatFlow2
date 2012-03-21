!##############################################################################
!# ****************************************************************************
!# <name> bcassemblybase </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic routines use to discretise analytically
!# given boundary conditions.
!#
!# The following routines can be found here:
!#
!# 1.) bcasm_getVertInBdRegion
!#      -> Determine the vertices in a boundary region
!#
!# 2.) bcasm_getEdgesInBdRegion
!#     -> Determine the edges in a boundary region
!#
!# 3.) bcasm_getElementsInBdRegion
!#     -> Determines elements, edges and vertices in a boundary region
!#
!# 4.) bcasm_getDOFsOnBoundary
!#     -> Calculates the DOF's on the whole boundary.
!#
!# 5.) bcasm_getDOFsInBDRegion
!#     -> Calculates the DOF's in a boundary region
!#
!# </purpose>
!##############################################################################

module bcassemblybase

!$use omp_lib
  use basicgeometry
  use boundary
  use dofmapping
  use element
  use fsystem
  use genoutput
  use spatialdiscretisation
  use storage
  use triangulation

  implicit none
  
  private

  public :: bcasm_getVertInBdRegion
  public :: bcasm_getEdgesInBdRegion
  public :: bcasm_getElementsInBdRegion
  public :: bcasm_getDOFsOnBoundary
  public :: bcasm_getDOFsInBDRegion
  
  public :: bcasm_getVertInBCregion  ! Deprecated
  public :: bcasm_getEdgesInBCregion ! Deprecated
  
  ! Deprecated interface, routine was renamed!
  public :: bcasm_getElementsInBCregion
  interface bcasm_getElementsInBCregion
    module procedure bcasm_getElementsInBdRegion
  end interface

contains

  ! ***************************************************************************
  
!<subroutine>

  subroutine bcasm_getVertInBCregion (rtriangulation,rboundary,rregion, &
      IminIndex,ImaxIndex,icount)
  
!<description>
  ! DEPRECATED! Not used anymore, since multiply connected boundary regions
  ! must be supported. Use bcasm_getElementsInBdRegion!
  !
  ! This routine receives a boundary region rregion describing a part on the
  ! boundary. According to the triangulation, this boundary region contains
  ! some vertices and edges on the boundary. The routine now figures out,
  ! which index in the IverticesAtBoundary in the triangulation structure
  ! of the vertices that are inside of this boundary region.
  !
  ! Each boundary region is simply connected, but they may cross the
  ! maximum parameter value. Therefore, we find wither one or two 'sets'
  ! of indices in the IverticesAtBoundary array:
  ! 1.) One simly connected index set:
  !     ..... Vertex Vertex Vertex ......
  !     Here is icount=1 and IminIndex(1) / ImaxIndex(1) gives the
  !     first/last index in IverticesAtBoundary
  ! 2.) Two sets crossing the maximum parameter value:
  !     Vertex Vertex Vertex ........... Vertex Vertex Vertex
  !     Here is icount=2, IminIndex(2) / ImaxIndex(2) gives the
  !     first/last vertex of the first set (which was found at first when
  !     starting to search in the given interval) and IminIndex(1) / ImaxIndex(1)
  !     the first/last vertex if the 2nd set.
!</description>

!<input>
  ! The triangulation structure.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! The description of the domain boundary
  type(t_boundary), intent(in) :: rboundary
  
  ! A boundary region structure describing a part of a boundary
  type(t_boundaryRegion), intent(in) :: rregion
!</input>

!<output>
  ! The index of the first vertex in IverticesAtBoundary belonging to the
  ! bonudary region rregion. = 0, if no vertices belong to the given region.
  integer, dimension(2), intent(out) :: IminIndex

  ! The index of the last vertex in IverticesAtBoundary belonging to the
  ! bonudary region rregion. = -1, if no vertices belong to the given region.
  integer, dimension(2), intent(out) :: ImaxIndex
  
  ! Number of index sets in IverticesAtBoundary belonging
  ! to rregion. IminIndex(1)..ImaxIndex(1) is the first set,...
  ! IminIndex(icount)..ImaxIndex(icount) is the last set.
  integer, intent(out) :: icount
!</output>

!</subroutine>

  ! local variables
  logical :: binside,bfinish
  real(DP), dimension(:), pointer :: p_DvertexParameterValue
  real(DP), dimension(2) :: dbegin,dend
  real(DP) :: dmaxpar
  integer, dimension(:), pointer :: p_IboundaryCpIdx
  integer :: i, ifoundRegions
  integer, dimension(2) :: Iindex
  
#if WARN_DEPREC
  call output_line ("Using deprecated feature. Please update your code.", &
      OU_CLASS_WARNING,OU_MODE_STD,"bcasm_getVertInBCregion")
#endif
  
  ! Get the parameter value array from the triangulation
  call storage_getbase_double (rtriangulation%h_DvertexParameterValue, &
                               p_DvertexParameterValue)
                               
  ! Get the array describing the start/end of each boundary component
  call storage_getbase_int (rtriangulation%h_IboundaryCpIdx, &
                            p_IboundaryCpIdx)
  
  ! Get the real min/max parameter value of the boundary region
  dmaxpar = boundary_dgetMaxParVal(rboundary, rregion%iboundCompIdx)
  dbegin(1) = rregion%dminParam
  dend(1) = rregion%dmaxParam
  if (dend(1)-dbegin(1) .ge. dmaxpar) then
    ! Occupy everything
    icount = 1
    dbegin(1) = 0.0_DP
    dend(1) = dmaxpar
    dbegin(2) = 0.0_DP
    dend(2) = -1.0_DP
  else
    ! Calculate the real min/max
    dbegin(1) = mod(rregion%dminParam,dmaxpar)
    dend(1) = mod(rregion%dmaxParam,dmaxpar)
    icount = 1
    if (dend(1) .lt. dbegin(1)) then
      ! probably two sets
      dbegin(2) = dend(1)
      dend(2) = dmaxpar
      dbegin(1) = 0.0_DP
      icount = 2
    end if
  end if
  
  ! All regions in [dbegin(.),dend(.)] are now in asecnding order.
  ! icount indicates the maximum number of regions we expect to have
  ! vertices inside.
  ! Remark: dbegin/dend are not used anymore here - maybe in a later
  ! version if necessary for some reasons...
  !
  ! For now, use a simple linear search.
  ! The boundary region structure tells us which boundary component
  ! to use; loop through all vertices there.
  
  ifoundRegions = 0
  IminIndex = 0
  ImaxIndex = -1
  binside = .false.
  bfinish = .false.
  
  ! Loop through all vertices in the current boundary component
  
  do I = p_IboundaryCpIdx(rregion%iboundCompIdx), &
         p_IboundaryCpIdx(rregion%iboundCompIdx+1)-1
         
    ! Check if the vertex is inside the region.
    if (boundary_isInRegion (rregion,rregion%iboundCompIdx,p_DvertexParameterValue(I))) then
      ! We are inside
      if (.not. binside) then
        ! We are inside for the first time
        binside = .true.
        IminIndex(ifoundRegions+1) = I
      end if
    else
      ! We are outside - for the first time?
      ! If yes, quit the loop.
      if (binside) then
        binside = .false.
        ImaxIndex(ifoundRegions+1) = I-1
        
        ! We completed the current region successfully
        ifoundRegions = ifoundRegions + 1
        
        ! Decrement icount. If it is still > 0, there is another region
        ! in dbegin/dend that may contain points
        icount = icount - 1
        
        if (icount .le. 0) then
          ! Finish that, we quit the search here as no more regions
          ! are expected.
          icount = ifoundRegions
          return
        end if
      end if
    end if
         
  end do
  
  ! The loop is completed. Question: Are we still inside a region or not?
  ! If yes...
  
  if (binside) then
    ! Save the last vertex number
    ImaxIndex(ifoundRegions+1) = p_IboundaryCpIdx(rregion%iboundCompIdx+1)-1
    ! Complete the region and finish.
    icount = ifoundRegions + 1
  else
    ! No we are not. So we were awaiting points for another part of the region
    ! that never came! Reset the last index pair and quit.
    IminIndex(ifoundRegions+1) = 0
    ImaxIndex(ifoundRegions+1) = -1
    icount = ifoundRegions
  end if
  
  if (icount .eq. 2) then
    ! We found the intervals in the 'wrong order'. When going through the boundary
    ! starting at a parameter value x, we would first find the 2nd segment
    ! x <= y..TMAX, then the first one 0..z. on the other hand set up IminIndex/
    ! ImaxIndex in the other order. So exchange IxxxIndex(1) with IxxxIndex(2)
    ! to get the right order again.
    Iindex(1) = IminIndex(2)
    IminIndex(2) = IminIndex(1)
    IminIndex(1) = Iindex(1)

    Iindex(2) = ImaxIndex(2)
    ImaxIndex(2) = ImaxIndex(1)
    ImaxIndex(1) = Iindex(2)
  end if
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine bcasm_getEdgesInBCregion (rtriangulation,rboundary,rregion, &
      IminIndex,ImaxIndex,icount)
  
!<description>
  ! DEPRECATED! Not used anymore, since multiply connected boundary regions
  ! must be supported. Use bcasm_getElementsInBdRegion!
  !
  ! This routine receives a boundary region rregion describing a part on the
  ! boundary. According to the triangulation, this boundary region contains
  ! some vertices and edges on the boundary. The routine now figures out,
  ! which  index in the IverticesAtBoundary in the triangulation structure
  ! of the vertices that are inside of this boundary region.
  !
  ! Each boundary region is simply connected, but they may cross the
  ! maximum parameter value. Therefore, we find wither one or two 'sets'
  ! of indices in the IverticesAtBoundary array:
  ! 1.) One simly connected index set:
  !     ..... Edge Edge Edge......
  !     Here is icount=1 and IminIndex(1) / ImaxIndex(1) gives the
  !     first/last index in IverticesAtBoundary
  ! 2.) Two sets crossing the maximum parameter value:
  !     Edge Edge Edge ........... Edge Edge Edge
  !     Here is icount=2, IminIndex(2) / ImaxIndex(2) gives the
  !     first/last vertex of the first set (which was found at first when
  !     starting to search in the given interval) and IminIndex(1) / ImaxIndex(1)
  !     the first/last vertex if the 2nd set.
!</description>

!<input>
  ! The triangulation structure.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! The description of the domain boundary
  type(t_boundary), intent(in) :: rboundary
  
  ! A boundary region structure describing a part of a boundary
  type(t_boundaryRegion), intent(in) :: rregion
!</input>

!<output>
  ! The index of the first edge in IedgesAtBoundary belonging to the
  ! bonudary region rregion. = 0, if no edges belong to the given region.
  integer, dimension(2), intent(out) :: IminIndex

  ! The index of the last edges in IedgesAtBoundary belonging to the
  ! bonudary region rregion. = -1, if no edges belong to the given region.
  integer, dimension(2), intent(out) :: ImaxIndex
  
  ! Number of index sets in IedgesAtBoundary belonging
  ! to rregion. IminIndex(1)..ImaxIndex(1) is the first set,...
  ! IminIndex(icount)..ImaxIndex(icount) is the last set.
  integer, intent(out) :: icount
!</output>

!</subroutine>

  ! local variables
  logical binside,bfinish
  real(DP), dimension(:), pointer :: p_DedgeParameterValue
  real(DP), dimension(2) :: dbegin,dend
  real(DP) :: dmaxpar
  integer, dimension(:), pointer :: p_IboundaryCpIdx
  integer :: i, ifoundRegions
  integer, dimension(2) :: Iindex
  
#if WARN_DEPREC
  call output_line ("Using deprecated feature. Please update your code.", &
      OU_CLASS_WARNING,OU_MODE_STD,"bcasm_getEdgesInBCregion")
#endif

  ! Get the parameter value array from the triangulation
  call storage_getbase_double (rtriangulation%h_DedgeParameterValue, &
                               p_DedgeParameterValue)
                               
  ! Get the array describing the start/end of each boundary component
  call storage_getbase_int (rtriangulation%h_IboundaryCpIdx, &
                            p_IboundaryCpIdx)
  
  ! Get the real min/max parameter value of the boundary region
  dmaxpar = boundary_dgetMaxParVal(rboundary, rregion%iboundCompIdx)
  dbegin(1) = rregion%dminParam
  dend(1) = rregion%dmaxParam
  if (dend(1)-dbegin(1) .ge. dmaxpar) then
    ! Occupy everything
    icount = 1
    dbegin(1) = 0.0_DP
    dend(1) = dmaxpar
    dbegin(2) = 0.0_DP
    dend(2) = -1.0_DP
  else
    ! Calculate the real min/max
    dbegin(1) = mod(rregion%dminParam,dmaxpar)
    dend(1) = mod(rregion%dmaxParam,dmaxpar)
    icount = 1
    if (dend(1) .lt. dbegin(1)) then
      ! probably two sets
      dbegin(2) = dend(1)
      dend(2) = dmaxpar
      dbegin(1) = 0.0_DP
      icount = 2
    end if
  end if
  
  ! All regions in [dbegin(.),dend(.)] are now in asecnding order.
  ! icount indicates the maximum number of regions we expect to have
  ! vertices inside.
  ! Remark: dbegin/dend are not used anymore here - maybe in a later
  ! version if necessary for some reasons...
  !
  ! For now, use a simple linear search.
  ! The boundary region structure tells us which boundary component
  ! to use; loop through all vertices there.
  
  ifoundRegions = 0
  IminIndex = 0
  ImaxIndex = -1
  binside = .false.
  bfinish = .false.
  
  ! Loop through all vertices in the current boundary component
  
  do I = p_IboundaryCpIdx(rregion%iboundCompIdx), &
         p_IboundaryCpIdx(rregion%iboundCompIdx+1)-1
         
    ! Check if the edge is inside the region.
    if (boundary_isInRegion (rregion,rregion%iboundCompIdx,p_DedgeParameterValue(I))) then
      ! We are inside
      if (.not. binside) then
        ! We are inside for the first time
        binside = .true.
        IminIndex(ifoundRegions+1) = I
      end if
    else
      ! We are outside - for the first time?
      ! If yes, quit the loop.
      if (binside) then
        binside = .false.
        ImaxIndex(ifoundRegions+1) = I-1
        
        ! We completed the current region successfully
        ifoundRegions = ifoundRegions + 1
        
        ! Decrement icount. If it is still > 0, there is another region
        ! in dbegin/dend that may contain points
        icount = icount - 1
        
        if (icount .le. 0) then
          ! Finish that, we quit the search here as no more regions
          ! are expected.
          icount = ifoundRegions
          return
        end if
      end if
    end if
         
  end do
  
  ! The loop is completed. Question: Are we still inside a region or not?
  ! If yes...
  
  if (binside) then
    ! Save the last edge number
    ImaxIndex(ifoundRegions+1) = p_IboundaryCpIdx(rregion%iboundCompIdx+1)-1
    ! Complete the region and finish.
    icount = ifoundRegions + 1
  else
    ! No we are not. So we were awaiting points for another part of the region
    ! that never came! Reset the last index pair and quit.
    IminIndex(ifoundRegions+1) = 0
    ImaxIndex(ifoundRegions+1) = -1
    icount = ifoundRegions
  end if

  if (icount .eq. 2) then
    ! We found the intervals in the 'wrong order'. When going through the boundary
    ! starting at a parameter value x, we would first find the 2nd segment
    ! x <= y..TMAX, then the first one 0..z. on the other hand set up IminIndex/
    ! ImaxIndex in the other order. So exchange IxxxIndex(1) with IxxxIndex(2)
    ! to get the right order again.
    Iindex(1) = IminIndex(2)
    IminIndex(2) = IminIndex(1)
    IminIndex(1) = Iindex(1)

    Iindex(2) = ImaxIndex(2)
    ImaxIndex(2) = ImaxIndex(1)
    ImaxIndex(1) = Iindex(2)
  end if

  end subroutine
    
  ! ***************************************************************************
  
!<subroutine>

  subroutine bcasm_getVertInBdRegion(rtriangulation,rboundaryRegion, &
      ncount, Ivertices, Itemp)
  
!<description>
  ! Calculates the vertices in a boundary region.
!</description>

!<input>
  ! The triangulation structure.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! A boundary region structure describing a part of a boundary
  type(t_boundaryRegion), intent(in) :: rboundaryRegion

  ! Temporary array. Must be large enough to hold all vertices in the region.
  ! Must be specified if Iedges is present!
  integer, dimension(:), intent(inout), optional :: Itemp
!</input>

!<output>
  ! Number of vertices in the region.
  ! If Iedges is not present: Estimated number of vertices in the boundary region;
  ! the number is larger or equal than the actual numer of vertices in the region.
  integer, intent(out) :: ncount
  
  ! A list of edges in the region. The array must be large enough to hold
  ! all vertices in the region.
  ! If not specified, the routine only calculates ncount, the number
  ! of edges on the boundary.
  integer, dimension(:), intent(out), optional :: Ivertices
!</output>

!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer :: i
    
    if (present(Ivertices) .and. .not. present (Itemp)) then
      call output_line ('Itemp not specified!', &
          OU_CLASS_ERROR,OU_MODE_STD,'bcasm_getVertInBdRegion')
      call sys_halt()
    end if
    
    ! Get ncount
    call bcasm_getElementsInBdRegion (rtriangulation,rboundaryRegion, ncount, &
        Ivertices,IvtLocal=Itemp)

    ! Calculate the edge numbers.
    if (present(Ivertices)) then
    
      call storage_getbase_int2d(rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    
      do i=1,ncount
        Ivertices(i) = p_IverticesAtElement(Itemp(i),Ivertices(i))
      end do
      
    end if

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine bcasm_getEdgesInBdRegion(rtriangulation,rboundaryRegion, ncount, Iedges, Itemp)
  
!<description>
  ! Calculates the edges in a boundary region.
!</description>

!<input>
  ! The triangulation structure.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! A boundary region structure describing a part of a boundary
  type(t_boundaryRegion), intent(in) :: rboundaryRegion

  ! Temporary array. Must be large enough to hold all edges in the region.
  ! Must be specified if Iedges is present!
  integer, dimension(:), intent(inout), optional :: Itemp
!</input>

!<output>
  ! Number of edges in the region.
  ! If Iedges is not present: Estimated number of edges in the boundary region;
  ! the number is larger or equal than the actual numer of edges in the region.
  integer, intent(out) :: ncount
  
  ! A list of edges in the region. The array must be large enough to hold
  ! all edges in the region.
  ! If not specified, the routine only calculates ncount, the number
  ! of edges on the boundary.
  integer, dimension(:), intent(out), optional :: Iedges
!</output>

!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer :: i
    
    if (present(Iedges) .and. .not. present (Itemp)) then
      call output_line ('Itemp not specified!', &
          OU_CLASS_ERROR,OU_MODE_STD,'bcasm_getEdgesInBdRegion')
      call sys_halt()
    end if
    
    ! Get ncount
    call bcasm_getElementsInBdRegion (rtriangulation,rboundaryRegion, ncount, &
        Iedges,IedgeLocal=Itemp)

    ! Calculate the edge numbers.
    if (present(Iedges)) then
    
      call storage_getbase_int2d(rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
    
      do i=1,ncount
        Iedges(i) = p_IedgesAtElement(Itemp(i),Iedges(i))
      end do
      
    end if

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine bcasm_getElementsInBdRegion (rtriangulation,rregion,ncount, &
                                          IelList,IelListIdx,IvtLocal,IedgeLocal)
  
!<description>
  ! This routine receives a boundary region rregion describing a part on the
  ! boundary. According to the triangulation, this boundary region contains
  ! some vertices and edges on the boundary.
  !
  ! The routine will search all elements that touch (with a vertex or an edge)
  ! the boundary region. All these elements are written to IelList.
  ! If a vertex of the element touches it, the corresponding entry in IvtLocal
  ! is set to the local vertex number of that vertex, otherwise it is set to 0.
  ! If an edge of the element touches it, the corresponding entry in IedgeLocal
  ! is set to the local edge number of that edge, otherwise it is set to to 0.
  !
  ! If neighter IvtLocal nor IedgeLocal is present, the routine calculates
  ! the maximum number of elements that the boundary region covers, including
  ! all edges and vertices.
!</description>

!<input>
  ! The triangulation structure.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! A boundary region structure describing a part of a boundary
  type(t_boundaryRegion), intent(in) :: rregion
!</input>

!<output>
  ! Number of edges in the boundary region, including repetitions.
  integer, intent(out) :: ncount

  ! OPTIONAL: Array that receives a list of elements that touch the boundary region.
  ! Elements may be inside here more than once.
  integer, dimension(:), intent(out), optional :: IelList

  ! OPTIONAL: Index list. Receives for every element touching the BC region an index
  ! into the IelementsAtBoundary where the element can be found.
  integer, dimension(:), intent(out), optional :: IelListIdx
  
  ! OPTIONAL: For each element on the boundary, local index of a vertex touching the
  ! boundary region. If more than one vertex touches it, the element is
  ! repeated in IelList and IvertexList contains all vertices that touch.
  ! If not present, elements that touch the boundary region only with a point
  ! are ignored.
  ! If neighter IvtLocal nor IedgeLocal is present, the routine calculates
  ! the maximum number of elements that the boundary region covers, including
  ! all edges and vertices.
  integer, dimension(:), intent(out), optional :: IvtLocal

  ! OPTIONAL: For each element on the boundary, local index of an edge touching the
  ! boundary region. If more than one edge touches it, the element is
  ! repeated in IelList and IedgeList contains all edges that touch.
  ! If not present, elements that touch the boundary region only with an edge
  ! are ignored.
  ! If neighter IvtLocal nor IedgeLocal is present, the routine calculates
  ! the maximum number of elements that the boundary region covers, including
  ! all edges and vertices.
  integer, dimension(:), intent(out), optional :: IedgeLocal
  
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    real(DP), dimension(:), pointer :: p_DedgeParameterValue
    integer, dimension(:), pointer :: p_IelementsAtBoundary, &
               p_IverticesAtBoundary, p_IedgesAtBoundary
    integer, dimension(:,:), pointer :: p_IverticesAtElement,p_IedgesAtElement
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer :: i,iidx,iel,ivt,iedge
    logical :: bvertexInside, bedgeInside,bcheckall
    
    ! If both destination arrays are missing, calculate the maximum number
    ! of elements in the region.
    bcheckall = .false.
    if (.not. present(IvtLocal) .and. .not. present(IedgeLocal)) &
        bcheckall = .true.
    
    ! Get the parameter value array from the triangulation
    call storage_getbase_int (rtriangulation%h_IelementsAtBoundary, &
                              p_IelementsAtBoundary)

    ! What spatial dimension are we?
    select case(rtriangulation%ndim)
    case (NDIM1D)
      ! Get the parameter value array from the triangulation
      if (present(IvtLocal) .or. bcheckall) then
        call storage_getbase_int (rtriangulation%h_IverticesAtBoundary, &
                                  p_IverticesAtBoundary)
        call storage_getbase_int2d (rtriangulation%h_IverticesAtElement, &
                                    p_IverticesAtElement)
      else
        nullify(p_IverticesAtBoundary)
        nullify(p_IverticesAtElement)
      end if
      
      ! Get the array describing the start/end of each boundary component
      call storage_getbase_int (rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
      
      ncount = 0
      if (present(IedgeLocal)) IedgeLocal = 0
      
      ! Cancel if this region does not exist.
      if ((rregion%iboundCompIdx .lt. 1) .or. &
          (rregion%iboundCompIdx .gt. rtriangulation%nbct)) then
        return
      end if
      
      bvertexInside = .false.

      ! Loop through all elements on the boundary. If we find any in the boundary
      ! region, take it!
      elementloop1d: do i = p_IboundaryCpIdx(rregion%iboundCompIdx),&
                            p_IboundaryCpIdx(rregion%iboundCompIdx+1)-1

        ! Got one. Add the element
        ncount = ncount + 1
        iel    = p_IelementsAtBoundary(i)
        
        if (present(IelList))    IelList(ncount)    = iel
        if (present(IelListIdx)) IelListIdx(ncount) = i
        if (present(IvtLocal))   IvtLocal(ncount)   = 0
        
        ! Remember the vertex
        if (present(IvtLocal)) then
          
          ivt = p_IverticesAtBoundary(i)
          
          ! Search the local number of the vertex that we found
          do iidx = 1,ubound(p_IverticesAtElement,1)
            if (p_IverticesAtElement(iidx,iel) .eq. ivt) then
              
              IvtLocal(ncount) = iidx
              
              ! Next element
              cycle elementloop1d
            end if
          end do
          
          ! Oops, the triangulation is destroyed!
          call output_line ('Local vertex number not found!', &
              OU_CLASS_ERROR,OU_MODE_STD,'bcasm_getElementsInBdRegion')
          call sys_halt()
        end if

      end do elementloop1d
      

    case (NDIM2D)
      ! Get the parameter value array from the triangulation
      if (present(IvtLocal) .or. bcheckall) then
        call storage_getbase_double (rtriangulation%h_DvertexParameterValue, &
                                     p_DvertexParameterValue)
        call storage_getbase_int (rtriangulation%h_IverticesAtBoundary, &
                                  p_IverticesAtBoundary)
        call storage_getbase_int2d (rtriangulation%h_IverticesAtElement, &
                                    p_IverticesAtElement)
      else
        nullify(p_DvertexParameterValue)
        nullify(p_IverticesAtBoundary)
        nullify(p_IverticesAtElement)
      end if
      
      ! Get the parameter value array from the triangulation
      if (present(IedgeLocal) .or. bcheckall) then
        call storage_getbase_double (rtriangulation%h_DedgeParameterValue, &
                                     p_DedgeParameterValue)
        call storage_getbase_int (rtriangulation%h_IedgesAtBoundary, &
                                  p_IedgesAtBoundary)
        call storage_getbase_int2d (rtriangulation%h_IedgesAtElement, &
                                    p_IedgesAtElement)
      else
        nullify(p_DedgeParameterValue)
        nullify(p_IedgesAtBoundary)
        nullify(p_IedgesAtElement)
      end if
                                 
      ! Get the array describing the start/end of each boundary component
      call storage_getbase_int (rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
      
      ncount = 0
      
      ! Cancel if this region does not exist.
      if ((rregion%iboundCompIdx .lt. 1) .or. &
          (rregion%iboundCompIdx .gt. rtriangulation%nbct)) then
        return
      end if
      
      bvertexInside = .false.
      bedgeInside = .false.
      
      ! Loop through all elements on the boundary. If we find any in the boundary
      ! region, take it!
      elementloop2d: do i = p_IboundaryCpIdx(rregion%iboundCompIdx),&
                            p_IboundaryCpIdx(rregion%iboundCompIdx+1)-1
        
        ! Check if the vertex or edge touches the boundary region
        if (present(IvtLocal) .or. bcheckall) &
            bvertexInside = boundary_isInRegion (rregion, rregion%iboundCompIdx,&
                                                 p_DvertexParameterValue(i))
        
        if (present(IedgeLocal) .or. bcheckall) &
            bedgeInside = boundary_isInRegion (rregion, rregion%iboundCompIdx,&
                                               p_DedgeParameterValue(i))
      
        if (bvertexInside .or. bedgeInside) then
          
          ! Got one. Add the element
          ncount = ncount + 1
          iel    = p_IelementsAtBoundary(i)
          
          if (present(IelList))    IelList(ncount)    = iel
          if (present(IelListIdx)) IelListIdx(ncount) = i
          if (present(IvtLocal))   IvtLocal(ncount)   = 0
          if (present(IedgeLocal)) IedgeLocal(ncount) = 0
          
          ! Remember the vertex and/or edge
          if (bvertexInside .and. present(IvtLocal)) then
            
            ivt = p_IverticesAtBoundary(i)
            
            ! Search the local number of the vertex that we found
            do iidx = 1,ubound(p_IverticesAtElement,1)
              if (p_IverticesAtElement(iidx,iel) .eq. ivt) then
                
                IvtLocal(ncount) = iidx
                
                if (bedgeInside .and. present(IedgeLocal)) then
                  ! The edge is also inside -- and has the same local number.
                  IedgeLocal(ncount) = iidx
                end if
                
                ! Next element
                cycle elementloop2d
              end if
            end do
            
            ! Oops, the triangulation is destroyed!
            call output_line ('Local vertex number not found!', &
                OU_CLASS_ERROR,OU_MODE_STD,'bcasm_getElementsInBdRegion')
            call sys_halt()
            
          end if
          
          ! This point is only reached, if the vertex is not inside.
          ! Check the edge.
          if (bedgeInside .and. present(IedgeLocal)) then
            
            iedge = p_IedgesAtBoundary(i)
            
            ! Search the local number of the edge that we found
            do iidx = 1,ubound(p_IedgesAtElement,1)
              if (p_IedgesAtElement(iidx,iel) .eq. iedge) then
                
                IedgeLocal(ncount) = iidx
                
                ! Next element
                cycle elementloop2d
              end if
            end do
            
            ! Oops, the triangulation is destroyed!
            call output_line ('Local edge number not found!', &
                OU_CLASS_ERROR,OU_MODE_STD,'bcasm_getElementsInBdRegion')
            call sys_halt()
          end if
          
        end if
      end do elementloop2d

    case default
      call output_line ('Unsupported spatial dimension!', &
          OU_CLASS_ERROR,OU_MODE_STD,'bcasm_getElementsInBdRegion')
      call sys_halt()
    end select
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine bcasm_getDOFsOnBoundary (rspatialDiscr, h_Idofs, ndofs)
  
!<description>
  ! Calculates the DOF's on the whole boundary.
!</description>

!<input>
  ! The discretisation structure of the underlying discretisation.
  type(t_spatialDiscretisation), intent(in), target :: rspatialDiscr
!</input>

!<inputoutput>
  ! Handle to an array that contains the DOF's.
  ! If <> ST_NOHANDLE, the size must be large enough to hold all
  ! DOF's. If this is set to ST_NOHANDLE, memory is automatically allocated
  ! or left =ST_NOHANDLE if there are no DOF's on the region.
  ! The caller must deallocate the memory!
  integer, intent(out) :: h_Idofs
!</inputoutput>

!<output>
  ! OPTINOAL: Number of DOF's on the boundary.
  integer, intent(out), optional :: ndofs
!</output>

!</subroutine>

    integer :: icount,icounttotal,ibc,iseg,ielement,ielidx,nve,I,J
    integer :: ipoint,ipoint1,ipoint2,ilocalEdge
    integer(I32) :: celement
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_boundaryRegion) :: rboundaryRegion
    integer, dimension(:,:), allocatable :: Idofs
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:), pointer :: p_IelementDistr
    integer, dimension(:), pointer :: p_Idofs,p_IdofsLocal
    integer, dimension(:), pointer :: p_IelementsAtBoundary
    integer, dimension(:), pointer :: p_InodalProperty

    ! List of element distributions in the discretisation structure
    type(t_elementDistribution), dimension(:), pointer :: p_RelementDistribution
    
    ! Initialise the output
    icounttotal = 0
    
    ! Do we have a boundary structure?
    if (associated(rspatialDiscr%p_rboundary)) then
      
      ! This is the default case in 2D/3D since a boundary structure
      ! is typically attached to the spatial discretisation structure

      ! Loop through the boundary components and segments.
      do ibc = 1,boundary_igetNBoundComp(rspatialDiscr%p_rboundary)
        do iseg = 1,boundary_igetNsegments(rspatialDiscr%p_rboundary,ibc)
          
          ! Calculate the number of DOF's there.
          call boundary_createRegion (rspatialDiscr%p_rboundary, &
              ibc, iseg, rboundaryRegion)
          call bcasm_getDOFsInBDRegion (rspatialDiscr, &
              rboundaryRegion, ndofs=icount)
          
          ! Sum up.
          icounttotal = icounttotal+icount
          
        end do
      end do
    
      if (present(ndofs)) then
        ndofs = icounttotal
      end if
      
    elseif (associated(rspatialDiscr%p_rtriangulation)) then

      ! This is the default case in 1D since no boundary structure is
      ! attached to the spatial discretisation structure. It also
      ! works in 2D/3D but there, the above approach should be more efficient

      ! For easier access:
      p_rtriangulation => rspatialDiscr%p_rtriangulation
      p_RelementDistribution => rspatialDiscr%RelementDistr

      ! Set pointers
      call storage_getbase_int(&
          p_rtriangulation%h_IelementsAtBoundary, p_IelementsAtBoundary)
      call storage_getbase_int2d(&
          p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
      call storage_getbase_int(&
          p_rtriangulation%h_InodalProperty, p_InodalProperty)

      ! Set total number of DOFs at the boundary
      icounttotal = size(p_IelementsAtBoundary)

      ! Reserve some memory to save temporarily all DOF`s of all boundary
      ! elements.
      ! We handle all boundary elements simultaneously - let us hope that there are
      ! never so many elements on the boundary that our memory runs out :-)
      allocate (Idofs(EL_MAXNBAS,icounttotal))

      Idofs(:,:) = 0

      ! Now the elements with indices iminidx..imaxidx in the ItrialElements
      ! of the triangulation are on the boundary. Some elements may appear
      ! twice (on edges e.g.) but we do not care.
      !
      ! Ask the DOF-mapping routine to get us those DOF`s belonging to elements
      ! on the boundary.
      !
      ! The 'mult' call only works on uniform discretisations. We cannot assume
      ! that and have to call dof_locGlobMapping for every element separately.
      if (rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) then
        ! All elements are of the samne type. Get it in advance.
        celement = rspatialDiscr%RelementDistr(1)%celement
        nve = elem_igetNVE (celement)
        call dof_locGlobMapping_mult(rspatialDiscr, &
            p_IelementsAtBoundary, Idofs)
      else
        ! Every element can be of different type.
        call storage_getbase_int(rspatialDiscr%h_IelementDistr,&
            p_IelementDistr)
        do ielement = 1,icounttotal
          call dof_locGlobMapping(rspatialDiscr,&
              p_IelementsAtBoundary(ielement), Idofs(:,ielement))
        end do
      end if
      
      ! Loop through the elements
      do ielidx = 1,icounttotal
        
        ! Get the element and information about it.
        ielement = p_IelementsAtBoundary (ielidx)
        
        ! Get the element type in case we do not have a uniform triangulation.
        ! Otherwise, celement was set to the trial element type above.
        if (rspatialDiscr%ccomplexity .ne. SPDISC_UNIFORM) then
          celement = p_RelementDistribution(p_IelementDistr(ielement))%celement
          nve = elem_igetNVE (celement)
        end if
        
        ! Now the element-dependent part. For each element type, we
        ! have to figure out which DOF`s are on the boundary!
        !
        ! We proceed as follows: We figure out, which DOF is on the
        ! boundary. Then, we ask our computation routine to
        ! calculate the necessary value and translate them into a
        ! DOF value.  All DOF values are collected later.
        select case (elem_getPrimaryElement(celement))
          
        case (EL_P0_1D,EL_P0_2D,EL_Q0_2D,EL_P0_3D,EL_Q0_3D)
          
          ! These elements have no DOF's associated to boundary edges.
          
        case (EL_P1_1D)
          
          ! Left point at the boundary?
          ipoint = p_IverticesAtElement(1,ielement)
          if (p_InodalProperty(ipoint) .ne. 0) then
            ! Set the DOF number < 0 to indicate that this DOF is on the boundary
            Idofs(1,ielidx) = -abs(Idofs(1,ielidx))
          end if
          
          ! Right point at the boundary?
          ipoint = p_IverticesAtElement(2,ielement)
          if (p_InodalProperty(ipoint) .ne. 0) then
            ! Set the DOF number < 0 to indicate that this DOF is on the boundary
            Idofs(2,ielidx) = -abs(Idofs(2,ielidx))
          end if

        case (EL_P1_2D,EL_Q1_2D)
          
          ! Loop over all vertices
          do ilocalEdge = 1, nve
            
            ipoint1 = p_IverticesAtElement(ilocalEdge,ielement)
            ipoint2 = p_IverticesAtElement(mod(ilocalEdge,nve)+1,ielement)
            
            if ((p_InodalProperty(ipoint1) .eq. 0) .or.&
                (p_InodalProperty(ipoint2) .eq. 0)) then
              ipoint1 = 0
              ipoint2 = 0
            end if

            ! Left point inside? -> Corresponding DOF must be computed
            if ( ipoint1 .ne. 0 ) then
              ! Set the DOF number < 0 to indicate that this DOF is in the region.
              Idofs(ilocalEdge,ielidx) = -abs(Idofs(ilocalEdge,ielidx))
            end if
            
            ! Right point inside? -> Corresponding DOF must be computed
            if ( ipoint2 .ne. 0 ) then
              ! Set the DOF number < 0 to indicate that this DOF is in the region.
              Idofs(mod(ilocalEdge,nve)+1,ielidx) = -abs(Idofs(mod(ilocalEdge,nve)+1,ielidx))
            end if
          end do

        case default
          call output_line('Unsupported element!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bcasm_getDOFsOnBoundary')
          call sys_halt()
        end select

      end do
            
    else
      call output_line ('Neither boundary nor triangulation available!', &
          OU_CLASS_ERROR,OU_MODE_STD,'bcasm_getDOFsOnBoundary')
      call sys_halt()
    end if
    
  
    ! Cancel if there are no DOF's.
    if (icounttotal .eq. 0) return
    
    ! Allocate and fetch the DOF's.
    if (h_Idofs .eq. ST_NOHANDLE) then
      call storage_new('bcasm_getDOFsOnBoundary', 'h_Idofs', &
          icounttotal, ST_INT, h_Idofs, ST_NEWBLOCK_NOINIT)
      call storage_getbase_int (h_Idofs,p_Idofs)
    else
      call storage_getbase_int (h_Idofs,p_Idofs)
      if (size(p_Idofs) .lt. icount) then
        call output_line ('Output array not large enough!', &
            OU_CLASS_ERROR,OU_MODE_STD,'bcasm_getDOFsOnBoundary')
        call sys_halt()
      end if
    end if


    ! Initialise the output
    icounttotal = 0
      
    ! Do we have a boundary structure?
    if (associated(rspatialDiscr%p_rboundary)) then

      ! Loop through the boundary components and segments.
      do ibc = 1,boundary_igetNBoundComp(rspatialDiscr%p_rboundary)
        do iseg = 1,boundary_igetNsegments(rspatialDiscr%p_rboundary,ibc)
          
          ! Calculate the DOF's there.
          call boundary_createRegion (rspatialDiscr%p_rboundary, &
              ibc, iseg, rboundaryRegion)
          
          p_IdofsLocal => p_Idofs(icounttotal+1:)
          
          ! Store in a subarray of p_Idofs.
          call bcasm_getDOFsInBDRegion (rspatialDiscr, &
              rboundaryRegion,ndofs=icount,IdofsArray=p_IdofsLocal)
          
          icounttotal = icounttotal+icount
      
        end do
      end do

    else
      
      ! Transfer the DOF`s and their values to these arrays.
      do J=1,size(Idofs,2)
        do I=1,size(Idofs,1)
          if (Idofs(I,J) < 0) then
            icounttotal = icounttotal + 1
            p_Idofs(icounttotal) = abs(Idofs(I,J))
          end if
        end do
      end do

      ! Remove temporary memory, finish.
      deallocate (Idofs)

    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine bcasm_getDOFsInBDRegion (rspatialDiscr, &
      rboundaryRegion, h_Idofs, ndofs, IdofsArray)
  
!<description>
  ! Calculates all DOF's associated to edges/vertices on the boundary
  ! in the boundary region rboundaryRegion.
!</description>

!<input>
  ! The discretisation structure of the underlying discretisation.
  type(t_spatialDiscretisation), intent(in), target :: rspatialDiscr

  ! A boundary-condition-region object, describing the position on the
  ! boundary where boundary conditions should be imposed.
  type(t_boundaryRegion), intent(in) :: rboundaryRegion
!</input>

!<inputoutput>
  ! OPTIONAL: Handle to an array that contains the DOF's.
  ! If the handle is set to ST_NOHANDLE, memory is automatically allocated.
  ! Otherwise, the memory must be large enough to hold all DOF's.
  ! If the pointer is ST_NOHANDLE, it stays =ST_NOHANDLE if there are no DOF's on
  ! the region. The caller must deallocate the memory!
  integer, intent(out), optional :: h_Idofs
  
  ! OPTIONAL: Number of DOF's in the boundary region.
  integer, intent(out), optional :: ndofs
  
  ! OPTIONAL: Array where to save the DOF's to. If this is specified,
  ! it must be large enough and all DOF's are written to here. h_Idofs is ignored
  ! in this case.
  integer, dimension(:), intent(out), target, optional :: IdofsArray
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,j,ilocalEdge,icount,ielidx
    integer(I32) :: celement
    integer :: ielement
    integer :: iedge,ipoint1,ipoint2,NVT
    type(t_triangulation), pointer :: p_rtriangulation
    integer, dimension(:), pointer :: p_IelementDistr
    integer, dimension(:,:), allocatable :: Idofs
    integer, dimension(:), pointer :: p_Idofs
    real(DP), dimension(:), pointer :: p_DedgeParameterValue
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_InodalProperty
    integer, dimension(:), allocatable :: IverticesAtBoundaryIdx
    integer, dimension(:), allocatable :: IedgesAtBoundaryIdx
    integer, dimension(:), allocatable :: IelementsAtBoundary
    integer, dimension(:), allocatable :: IelementsAtBoundaryIdx
    
    integer :: nve,nnve
    real(DP) :: Dpar
    
    ! Position of cubature points for 2-point Gauss formula on an edge.
    ! Used for Q2T.
    real(DP), parameter :: Q2G1 = -0.577350269189626_DP !-SQRT(1.0_DP/3.0_DP)
    real(DP), parameter :: Q2G2 =  0.577350269189626_DP ! SQRT(1.0_DP/3.0_DP)
    
    ! List of element distributions in the discretisation structure
    type(t_elementDistribution), dimension(:), pointer :: p_RelementDistribution

    ! For easier access:
    p_rtriangulation => rspatialDiscr%p_rtriangulation
    p_RelementDistribution => rspatialDiscr%RelementDistr
    
    call storage_getbase_int (p_rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    call storage_getbase_int (p_rtriangulation%h_InodalProperty,&
        p_InodalProperty)
    call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    
    ! The edge-at-element array may not be initialised.
    if (p_rtriangulation%h_IedgesAtElement .ne. ST_NOHANDLE) then
      call storage_getbase_int2D(p_rtriangulation%h_IedgesAtElement,&
          p_IedgesAtElement)
    else
      nullify(p_IedgesAtElement)
    end if
    
    ! The parameter value arrays may not be initialised.
    if (p_rtriangulation%h_DedgeParameterValue .ne. ST_NOHANDLE) then
      call storage_getbase_double(p_rtriangulation%h_DedgeParameterValue,&
          p_DedgeParameterValue)
    else
      nullify(p_DedgeParameterValue)
    end if

    if (p_rtriangulation%h_DvertexParameterValue .ne. ST_NOHANDLE) then
      call storage_getbase_double(p_rtriangulation%h_DvertexParameterValue,&
          p_DvertexParameterValue)
    else
      nullify(p_DvertexParameterValue)
    end if

    NVT = p_rtriangulation%NVT
    nnve = p_rtriangulation%NNVE

    if (rspatialDiscr%ccomplexity .ne. SPDISC_UNIFORM) then
      ! Every element can be of different type.
      call storage_getbase_int(rspatialDiscr%h_IelementDistr,&
          p_IelementDistr)
    else
      ! All elements are of the samne type. Get it in advance.
      celement = rspatialDiscr%RelementDistr(1)%celement
      nve = elem_igetNVE (celement)
    end if
    
    ! We have to deal with all DOF`s on the boundary. This is highly element
    ! dependent and therefore a little bit tricky :(
    !
    ! As we are in 2D, we can use parameter values at first to figure out,
    ! which points and which edges are on the boundary.
    ! What we have is a boundary segment. Now ask the boundary-index routine
    ! to give us the vertices and edges on the boundary that belong to
    ! this boundary segment.
    
    allocate(IverticesAtBoundaryIdx(p_rtriangulation%NVBD))
    allocate(IedgesAtBoundaryIdx(p_rtriangulation%NVBD))
    allocate(IelementsAtBoundary(p_rtriangulation%NVBD))
    allocate(IelementsAtBoundaryIdx(p_rtriangulation%NVBD))
    
    IelementsAtBoundary = 0
    IelementsAtBoundaryIdx = 0

    if (associated(p_IedgesAtElement) .and. associated(p_DedgeParameterValue)) then
      call bcasm_getElementsInBdRegion (p_rtriangulation, rboundaryRegion, &
          icount, IelementsAtBoundary, IelementsAtBoundaryIdx, &
          IverticesAtBoundaryIdx, IedgesAtBoundaryIdx)
    else
      call bcasm_getElementsInBdRegion (p_rtriangulation, rboundaryRegion, &
          icount, IelementsAtBoundary, IelementsAtBoundaryIdx, &
          IverticesAtBoundaryIdx)
      IedgesAtBoundaryIdx = 0
    end if
                                   
    if (icount .eq. 0) then
      deallocate(IverticesAtBoundaryIdx)
      deallocate(IedgesAtBoundaryIdx)
      deallocate(IelementsAtBoundary)
      deallocate(IelementsAtBoundaryIdx)
      return
    end if
                                   
    ! Reserve some memory to save temporarily all DOF`s of all boundary
    ! elements.
    ! We handle all boundary elements simultaneously - let us hope that there are
    ! never so many elements on the boundary that our memory runs out :-)
    allocate (Idofs(EL_MAXNBAS,icount))
    
    Idofs(:,:) = 0

    ! Now the elements with indices iminidx..imaxidx in the ItrialElements
    ! of the triangulation are on the boundary. Some elements may appear
    ! twice (on edges e.g.) but we do not care.
    !
    ! Ask the DOF-mapping routine to get us those DOF`s belonging to elements
    ! on the boundary.
    !
    ! The 'mult' call only works on uniform discretisations. We cannot assume
    ! that and have to call dof_locGlobMapping for every element separately.
    if (rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) then
      call dof_locGlobMapping_mult(rspatialDiscr, &
                IelementsAtBoundary(1:icount), Idofs)
    else
      do ielement = 1,icount
        call dof_locGlobMapping(rspatialDiscr, IelementsAtBoundary(ielement),&
            Idofs(:,ielement))
      end do
    end if

    ! Loop through the elements
    do ielidx = 1,icount

      ! Get the element and information about it.
      ielement = IelementsAtBoundary (ielidx)
      
      ! Index in the boundary arrays.
      I = IelementsAtBoundaryIdx (ielidx)
      
      ! Get the element type in case we do not have a uniform triangulation.
      ! Otherwise, celement was set to the trial element type above.
      if (rspatialDiscr%ccomplexity .ne. SPDISC_UNIFORM) then
        celement = p_RelementDistribution(p_IelementDistr(ielement))%celement
        nve = elem_igetNVE (celement)
      end if
        
      ilocaledge = 0
      
      if (IverticesAtBoundaryIdx(ielidx) .ne. 0) then
        ! Get the local index of the edge -- it coincides with the local index
        ! of the vertex.
        ilocaledge = IverticesAtBoundaryIdx(ielidx)
        ipoint1 = p_IverticesAtElement(IverticesAtBoundaryIdx(ielidx),ielement)
        ipoint2 = p_IverticesAtElement(mod(IverticesAtBoundaryIdx(ielidx),nve)+1,ielement)
      else
        ipoint1 = 0
        ipoint2 = 0
      end if

      if (associated(p_IedgesAtElement) .and. (IedgesAtBoundaryIdx(ielidx) .ne. 0)) then
        ! Get the local index of the edge -- it coincides with the local index
        ! of the vertex.
        ilocaledge = IedgesAtBoundaryIdx(ielidx)

        ! Get the edge
        iedge = p_IedgesAtElement(IedgesAtBoundaryIdx(ielidx),ielement)
      else
        iedge = 0
      end if
      
      dpar = -1.0_DP
    
      ! Now the element-dependent part. For each element type, we have to
      ! figure out which DOF`s are on the boundary!
      !
      ! We proceed as follows: We figure out, which DOF is on the
      ! boundary. Then, we ask our computation routine to calculate
      ! the necessary value and translate them into a DOF value.
      ! All DOF values are collected later.
      select case (elem_getPrimaryElement(celement))
      
      case (EL_P0_1D,EL_P0_2D,EL_Q0_2D,EL_P0_3D,EL_Q0_3D)

        ! This element has no DOF's associated to boundary edges.
        
      case (EL_P1_1D)

        ! Left point inside? -> Corresponding DOF must be computed
        ! if it is lokcated at the boundary
        if ( (ipoint1 .ne. 0) .and. (p_InodalProperty(ipoint1) .ne. 0) ) then
          ! Set the DOF number < 0 to indicate that this DOF is in the region.
          Idofs(ilocalEdge,ielidx) = -abs(Idofs(ilocalEdge,ielidx))
        end if

        ! Right point inside? -> Corresponding DOF must be computed
        ! if it is lokcated at the boundary
        if ( (ipoint2 .ne. 0) .and. (p_InodalProperty(ipoint2) .ne. 0) ) then
          ! Set the DOF number < 0 to indicate that this DOF is in the region.
          Idofs(mod(ilocalEdge,2)+1,ielidx) = -abs(Idofs(mod(ilocalEdge,2)+1,ielidx))
        end if
        
      case (EL_P1_2D,EL_Q1_2D)

        ! Left point inside? -> Corresponding DOF must be computed
        if ( ipoint1 .ne. 0 ) then
          ! Set the DOF number < 0 to indicate that this DOF is in the region.
          Idofs(ilocalEdge,ielidx) = -abs(Idofs(ilocalEdge,ielidx))
        end if
        
        ! The right point does not have to be checked! It comes later
        ! with the next edge. The situation when an element crosses the
        ! maximum parameter value with its boundary is handled by the
        ! outer DO-LOOP:
        ! A boundary region with parameter value e.g. [3.0,TMAX]
        ! will produce two index sets: One index set for [0.0, 0.0]
        ! and one for [3.0, TMAX).
        
      case (EL_DG_P1_2D)

        ! Left point inside? -> Corresponding DOF must be computed
        if ( ipoint1 .ne. 0 ) then
          ! Set the DOF number < 0 to indicate that this DOF is in the region.
          Idofs(ilocalEdge,ielidx) = -abs(Idofs(ilocalEdge,ielidx))
        end if
        
        ! Right point inside? -> Corresponding DOF must be computed
        if ( ipoint2 .ne. 0 ) then
          ! Set the DOF number < 0 to indicate that this DOF is in the region.
          Idofs(ilocalEdge,ielidx) = -abs(Idofs(ilocalEdge,ielidx))
          Idofs(mod(ilocalEdge,3)+1,ielidx) = -abs(Idofs(mod(ilocalEdge,3)+1,ielidx))
        end if

      case (EL_P2_2D,EL_Q2_2D)

        ! Left point inside? -> Corresponding DOF must be computed
        if ( ipoint1 .ne. 0 ) then
          ! Set the DOF number < 0 to indicate that this DOF is in the region.
          Idofs(ilocalEdge,ielidx) = -abs(Idofs(ilocalEdge,ielidx))
        end if
        
        ! The right point does not have to be checked! It comes later
        ! with the next edge. The situation when an element crosses the
        ! maximum parameter value with its boundary is handled by the
        ! outer DO-LOOP:
        ! A boundary region with parameter value e.g. [3.0,TMAX]
        ! will produce two index sets: One index set for [0.0, 0.0]
        ! and one for [3.0, TMAX).
        !
        ! Edge inside? -> Calculate point value on midpoint of edge iedge
        if ( iedge .ne. 0 ) then
          
          ! Set the DOF number < 0 to indicate that this DOF is in the region.
          Idofs(ilocalEdge+nve,ielidx) = -abs(Idofs(ilocalEdge+nve,ielidx))
              
          ! The element midpoint does not have to be considered, as it cannot
          ! be on the boundary.
        end if

      case (EL_QP1_2D)
        ! Three DOF`s: Function value in the element midpoint
        ! and derivatives.
        ! No DOF is on a boundary edge.

      case (EL_P1T_2D)

        ! Edge midpoint based element.
        !
        ! Edge inside? -> Calculate point value on midpoint of edge iedge
        if ( iedge .ne. 0 ) then
          ! Set the DOF number < 0 to indicate that this DOF is in the region.
          Idofs(ilocalEdge,ielidx) = -abs(Idofs(ilocalEdge,ielidx))
        end if

      case (EL_Q1T_2D,EL_Q1TB_2D)
      
        ! The Q1T-element has different variants. Check which variant we have
        ! and choose the right way to calculate boundary values.
      
        if (iand(celement,int(2**16,I32)) .ne. 0) then
        
          ! Integral mean value based element.
          
          ! Edge inside? -> Calculate integral mean value over the edge
          if ( iedge .ne. 0 ) then
          
            ! Set the DOF number < 0 to indicate that this DOF is in the region.
            Idofs(ilocalEdge,ielidx) = -abs(Idofs(ilocalEdge,ielidx))
              
          end if
                                      
        else
          
          ! Edge midpoint based element.
          !
          ! Edge inside? -> Calculate point value on midpoint of edge iedge
          if ( iedge .ne. 0 ) then
            ! Set the DOF number < 0 to indicate that this DOF is in the region.
            Idofs(ilocalEdge,ielidx) = -abs(Idofs(ilocalEdge,ielidx))
          end if

        end if

      case (EL_Q2T_2D,EL_Q2TB_2D)
      
        ! The Q2T-element is only integral mean value based.
        ! On the one hand, we have integral mean values over the edges.
        ! On the other hand, we have integral mean values of function*parameter
        ! value on the edge.
        !
        ! Edge inside?
        if ( iedge .ne. 0 ) then
          
          ! Set the DOF number < 0 to indicate that this DOF is in the region.
          Idofs(ilocalEdge,ielidx) = -abs(Idofs(ilocalEdge,ielidx))

          ! Set the DOF number < 0 to indicate that this DOF is in the region.
          Idofs(ilocalEdge+nve,ielidx) = -abs(Idofs(ilocalEdge+nve,ielidx))
          
        end if

      case (EL_Q3T_2D)
      
        ! The Q3T-element is only integral mean value based.
        ! On the one hand, we have integral mean values over the edges.
        ! On the other hand, we have integral mean values of function*parameter
        ! value on the edge.
        !
        ! Edge inside?
        if ( iedge .ne. 0 ) then
          
          ! Set the DOF number < 0 to indicate that this DOF is in the region.
          Idofs(ilocalEdge,ielidx) = -abs(Idofs(ilocalEdge,ielidx))
          Idofs(ilocalEdge+nve,ielidx) = -abs(Idofs(ilocalEdge+nve,ielidx))
          Idofs(ilocalEdge+2*nve,ielidx) = -abs(Idofs(ilocalEdge+2*nve,ielidx))
          
        end if
      case default
        call output_line('Unsupported element!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bcasm_getDOFsInBDRegion')
        call sys_halt()
      
      end select
      
    end do

    ! Temp arrays no more necessary
    deallocate(IverticesAtBoundaryIdx)
    deallocate(IedgesAtBoundaryIdx)
    deallocate(IelementsAtBoundary)
    deallocate(IelementsAtBoundaryIdx)
    
    ! Now count how many values we actually have.
    icount = 0
    do J=1,size(Idofs,2)
      do I=1,size(Idofs,1)
        if (Idofs(I,J) < 0) icount = icount + 1
      end do
    end do

    nullify(p_Idofs)
        
    if ((icount .gt. 0) .and. (present(IdofsArray) .or. present(h_Idofs))) then
      if (present(IdofsArray)) then
        ! This array has priority
        p_Idofs => IdofsArray
      else
        ! Allocate arrays for storing these DOF`s and their values - if values are
        ! computed.
        if (h_Idofs .eq. ST_NOHANDLE) then
          call storage_new('bcasm_getDOFsInBDRegion', 'h_Idofs', &
              icount, ST_INT, h_Idofs, ST_NEWBLOCK_NOINIT)
          call storage_getbase_int (h_Idofs,p_Idofs)
        else
          call storage_getbase_int (h_Idofs,p_Idofs)
        end if
      end if

      if (associated(p_Idofs)) then
        if (size(p_Idofs) .lt. icount) then
          call output_line ('Output array not large enough!', &
              OU_CLASS_ERROR,OU_MODE_STD,'bcasm_getDOFsInBDRegion')
          call sys_halt()
        end if
    
        ! Transfer the DOF`s and their values to these arrays.
        icount = 0
        do J=1,size(Idofs,2)
          do I=1,size(Idofs,1)
            if (Idofs(I,J) < 0) then
              icount = icount + 1
              p_Idofs(icount) = abs(Idofs(I,J))
            end if
          end do
        end do
      end if
    
    end if
    
    if (present(ndofs)) ndofs=icount
    
    ! Remove temporary memory, finish.
    deallocate (Idofs)

  end subroutine

end module bcassemblybase
  

    
