!##############################################################################
!# ****************************************************************************
!# <name> graph </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic definitions and routines to maintain
!# graphs. A graph consists of a set of vertices V and a set of edge E
!# connecting vertices. Graphs are stored as adjacency matrices which
!# can be modified dynamically. A new vertex can be added or an
!# existing vertex can be removed. Moreover, edges can be inserted
!# between two vertices or removed.  This module provides data
!# structures for both directed and undirected graphs.  You should be
!# aware of the fact, that adjacency lists are perhaps not the optimal
!# choice for nearly complete graphs, that is, graphs which possess
!# the majority of all possible edges. In this case it would make
!# sense to adopt a static matrix data structure. Be warned, this is
!# not implemented in this module so that the performance of the
!# provided routines may deteriorate for nearly complete graphs.
!#
!# To keep things simple, graphs can be seen as sparce matrices stored
!# in the CSR format (either 7 or 9) which can be modified
!# dynamically. Hence, it makes sence to provide conversion routines
!# between a scalar matrix and a graph. If you are familiar with the
!# CSR format used for scalar matrices then this hint will be helpful:
!# The list of vertices is stored in the Kld array whereas the edges
!# are associated with the Kcol array. In short, an edge (I,J) exists
!# if and only if there is an entry J in the subarray
!# Kcol(Kld(i):Kld(i+1)-1).
!#
!#
!#
!# The following routines can be found in this module:
!#
!#  1.) grph_createGraph
!#      -> Create graph directly
!#
!#  2.) grph_releaseGraph
!#      -> Release an existing graph
!#
!#  3.) grph_printGraph
!#      -> Print the graph to screen
!#
!#  4.) grph_hasVertex
!#      -> Check if the graph has the given vertex
!#
!#  5.) grph_insertVertex
!#      -> Insert vertex into the graph. Do nothing if vertex already exists.
!#
!#  6.) grph_removeVertex
!#      -> Remove vertex from graph. Do nothing if vertex does not exist.
!#
!#  7.) grph_hasEdge
!#      -> Check if the graph has the given edge
!#
!#  8.) grph_insertEdge
!#      -> Insert edge into the graph. Do nothing if edge already exists.
!#
!#  9.) grph_removeEdge
!#      -> Remove edge from graph. Do nothing if edge does not exist.
!#
!# 10.) grph_infoGraph
!#      -> Print information about the graph.
!#
!# 11.) grph_duplicateGraph
!#      -> Create a duplicate / backup of a graph.
!#
!# 12.) grph_restoreGraph
!#      -> Restores a graph previously backed up with grph_duplicateGraph
!#
!# </purpose>
!##############################################################################

module graph

!$use omp_lib
  use arraylistInt
  use fsystem
  use genoutput
  use mapInt_Int
  use mapbase
  use storage

  implicit none

  private
  public :: t_graph
  public :: grph_createGraph
  public :: grph_releaseGraph
  public :: grph_printGraph
  public :: grph_hasVertex
  public :: grph_insertVertex
  public :: grph_removeVertex
  public :: grph_hasEdge
  public :: grph_insertEdge
  public :: grph_removeEdge
  public :: grph_infoGraph
  public :: grph_duplicateGraph
  public :: grph_restoreGraph

!<constants>
!<constantblock description="Global format flags for graphs">

  ! Unidentified graph format
  integer, parameter, public :: GRPH_GRAPHUNDEFINED = 0

  ! Identifier for directed graph with unordered entries
  integer, parameter, public :: GRPH_GRAPHUNORDERED_DIRECTED = 1

  ! Identifier for directed graph with ordered entries
  integer, parameter, public :: GRPH_GRAPHORDERED_DIRECTED = 2

  ! Identifier for undirected graph with (un)ordered entries
  integer, parameter, public :: GRPH_UNDIRECTED = 3

!</constantblock>
!</constants>

!<types>
!<typeblock>

  ! This type contains all data structures to handle graphs.
  type t_graph

    ! Format-tag. Identifies the format of the graph
    ! Can take one of the GRPH_GRAPHx format flags.
    integer :: cgraphFormat = GRPH_GRAPHUNDEFINED

    ! Number of vertices
    integer :: NVT = 0

    ! Number of edges
    integer :: NEDGE = 0

    ! Is set to true, if the graph is dense, that is, a graph with NVT
    ! vertices contains the vertices labeled 1,...,NVT. In this case
    ! the vertex number coincides with the table under which its
    ! adjacency list is stored.  If the graph is not marked as dense,
    ! then for each vertex IVT the table number is stored in the
    ! vertex list and must be looked-up each time, the adjacency
    ! information is required.  Note that dense graphs can only store
    ! positive vertex labels !!!
    logical :: bisDense = .true.

    ! Map for the list of vertices
    type(t_mapInt_Int) :: rVertices

    ! Array of lists for the edges
    type(t_arraylistInt) :: rEdges

  end type t_graph

!</typeblock>
!</types>

contains

!<subroutine>

  subroutine grph_createGraph(rgraph,cgraphFormat,nvtMax,nedgeMax,bisDense)

!<description>
    ! This routine creates a graph rgraph. The graph format must be
    ! one of the GRPH_GRAPHx format specifiers. The parameters nvtMax
    ! and nedgeMax define the maximum number of vertices/egdes that
    ! can be initially stored in the graph. If either more vertices or
    ! edges are inserted into the graph, then reallocation of memory
    ! will be performed.  The optional parameter bisDense defines, if
    ! the graph should be treated as dense graph or not. In a dense
    ! graph, the adjacency list of vertex ivt is stored in table
    ! itable=ivt. It makes sense to use dense graphs, e.g., if NVT
    ! vertices numbered by 1..NVT are present in the graph so that the
    ! table look-up can be performed in time O(1). If the graph
    ! consists of NVT vertices which can be labeled arbitrarily, e.g.,
    ! <tex>1.. ($1^k$)*NVT, k >> 1</tex>, then an enormous amount of
    ! memory would be wasted for the look-up tables. Then it makes
    ! sense to define the graph as "not dense". In this case, the
    ! tables are used sequentially if new vertices are inserted but
    ! finding the correct table number for vertex ivt requires a
    ! search in the map, i.e. time O(log NVT).
!</description>

!<input>
    ! Format of the graph
    integer, intent(in) :: cgraphFormat

    ! Maximum number of vertices
    integer, intent(in) :: nvtMax

    ! Maximum number of edges
    integer, intent(in) :: nedgeMax

    ! OPTIONAL: Indicates if the graph is dense or not
    logical, intent(in), optional :: bisDense
!</input>

!<output>
    ! The sparsity graph
    type(t_graph), intent(out) :: rgraph
!</output>
!</subroutine>

    ! Set dimensions
    rgraph%NVT    = 0
    rgraph%NEDGE  = 0

    if (present(bisDense)) rgraph%bisDense=bisDense
    rgraph%cgraphFormat = cgraphFormat
    
    ! Create map for the list of vertices
    if (rgraph%bisDense) then
      call map_create(rgraph%rVertices,rgraph%NVT+1,0)
    else
      call map_create(rgraph%rVertices,rgraph%NVT+1,1)
    end if
    
    ! Create array of lists for edges
    call alst_create(rgraph%rEdges,nvtMax,nedgeMax)
    
  end subroutine grph_createGraph

  ! ***************************************************************************

!<subroutine>

  subroutine grph_releaseGraph(rgraph)

!<description>
    ! This routine releases an existing graph rgraph.
!</description>

!<inputoutput>
    ! The graph that should be released
    type(t_graph), intent(inout) :: rgraph
!</inputoutput>
!</subroutine>

    ! Release map for the list of vertices
    call map_release(rgraph%rVertices)

    ! Release array of lists for edges
    call alst_release(rgraph%rEdges)

    ! Reset data
    rgraph%cgraphFormat = GRPH_GRAPHUNDEFINED
    rgraph%NVT          = 0
    rgraph%NEDGE        = 0
    rgraph%bisDense     = .true.

  end subroutine grph_releaseGraph

  ! ***************************************************************************

!<subroutine>

  subroutine grph_printGraph(rgraph)

!<description>
    ! This subroutine prints the content of the graph.
!</description>

!<inputoutput>
    ! The graph
    type(t_graph), intent(inout) :: rgraph
!</inputoutput>
!</subroutine>

    ! local variables
    type(it_mapInt_Int) :: rmapIter
    integer, dimension(:), pointer :: p_Idata
    integer :: ikey,itable

    ! Are we dense graph?
    if (rgraph%bisDense) then

      ! Output dense graph
      rmapIter = map_begin(rgraph%rVertices)
      do while (.not.map_isNull(rmapIter))

        ! Get key value
        ikey = map_get(rgraph%rVertices,rmapIter)

        ! Print neighbours of vertex IKEY
        call output_line('Vertex number: '//trim(sys_siL(ikey,15)))
        call alst_print(rgraph%rEdges,ikey)

        ! Proceed to next vertex
        call map_next(rmapIter)
      end do

    else

      ! Output sparse graph
      rmapIter = map_begin(rgraph%rVertices)
      do while (.not.map_isNull(rmapIter))

        ! Get key value
        ikey = map_get(rgraph%rVertices,rmapIter)

        ! Get auxiliary data value
        call map_getbase_data(rgraph%rVertices,rmapIter,p_Idata)
        itable = p_Idata(1)

        ! Print neighbours of vertex IKEY
        call output_line('Vertex number: '//trim(sys_siL(ikey,15)))
        call alst_print(rgraph%rEdges,itable)

        ! Proceed to next vertex
        call map_next(rmapIter)
      end do

    end if

  end subroutine grph_printGraph

  ! ***************************************************************************

!<function>

  function grph_hasVertex(rgraph,iVertex) result(bexists)

!<description>
    ! This function returns TRUE if the graph has the given vertex.
    ! Otherwise, it returns FALSE.
!</description>

!<input>
    ! Number of the vertex
    integer, intent(in) :: iVertex
!</input>

!<inputoutput>
    ! The graph
    type(t_graph), intent(inout) :: rgraph
!</inputoutput>

!<result>
    ! Flag for existence of the vertex
    logical :: bexists
!</result>
!</function>

    ! Search in the vertex list
    bexists  = .not.map_isNull(map_find(rgraph%rVertices,iVertex))

  end function grph_hasVertex

  ! ***************************************************************************

!<subroutine>

  subroutine grph_insertVertex(rgraph,iVertex)

!<description>
    ! This subroutine inserts the vertex with number iVertex into the
    ! graph. In addition, the trivial edge (iVertex,iVertex) is inserted.
!</description>

!<input>
    ! Number of the vertex
    integer, intent(in) :: iVertex
!</input>

!<inputoutput>
    ! The graph
    type(t_graph), intent(inout) :: rgraph
!</inputoutput>

!</subroutine>

    ! local variabes
    type(it_mapInt_Int) :: rmapIter
    integer :: itable

    ! Try to add entry with key iVertex
    if (rgraph%bisDense) then
      itable   = iVertex
      rmapIter = map_insert(rgraph%rVertices,iVertex)
    else
      itable   = rgraph%NVT+1
      rmapIter = map_insert(rgraph%rVertices,iVertex,(/itable/))
    end if

    ! Check for special flag
    if (map_hasSpec(rmapIter, MAP_MSPEC_EXISTS)) return

    ! Increase number of vertices by one
    rgraph%NVT = rgraph%NVT+1

    ! Add trivial edge (iVertex,iVertex) to the list of edges
    call alst_push_front(rgraph%rEdges,itable,iVertex)
    
    ! Increase number of edges by one
    rgraph%NEDGE = rgraph%NEDGE+1

  end subroutine grph_insertVertex

  ! ***************************************************************************

!<subroutine>

  subroutine grph_removeVertex(rgraph,iVertex,ireplacementVertex)

!<description>
    ! This subroutine removes the vertex with number iVertex from the
    ! graph and eliminates all of its incoming and outgoing edges. If
    ! the optional parameter ireplacemantVertex is present, then vertex
    ! iVertex is removed and the vertex with number ireplacementVertex
    ! is moved at its former position. This can be useful, if one
    ! needs a complete set of vertices, e.g., 1..NVT without "holes"
    ! at the positions of removed vertices.
!</description>

!<input>
    ! Number of the vertex
    integer, intent(in) :: iVertex

    ! OPTIONAL: Number of the replacement vertex
    integer, intent(in), optional :: ireplacementVertex
!</input>

!<inputoutput>
    ! The graph
    type(t_graph), intent(inout) :: rgraph
!</inputoutput>
!</subroutine>

    ! local variabes
    type(it_mapInt_Int) :: rmapIter
    type(it_arraylistInt) :: ralstIter,rposition
    integer, dimension(:), pointer :: p_Idata
    integer :: ireplaceTable,ireplaceVertex,itable,jVertex,jtable
    logical :: bdoReplace

    ! Check if vertex is present in map
    rmapIter = map_find(rgraph%rVertices,iVertex)
    if (map_isNull(rmapIter)) then
      call output_line('Vertex does not exist in graph!',&
          OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
      call sys_halt()
    end if

    ! Replace vertex or leave "holes"?
    bdoReplace = present(ireplacementVertex)
    if (bdoReplace) then
      ireplaceVertex = ireplacementVertex
    else
      ireplaceVertex = 0
    end if

    ! Are we dense graph?
    if (rgraph%bisDense) then

      ! What graph format are ?
      select case(rgraph%cgraphFormat)

      case(GRPH_GRAPHUNORDERED_DIRECTED,&
           GRPH_GRAPHORDERED_DIRECTED)
        call output_line('We cannot remove directed graphs at the moment!',&
            OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
        call sys_halt()

      case(GRPH_UNDIRECTED)
        ! We are in the lucky position, that the graph is undirected,
        ! that is, there exits an edge (i,j) if and only if there
        ! exists the edge (j,i).  Hence, we do not have to loop
        ! through the list of all vertices but it suffices to visit
        ! those which are present in the adjacency list of vertex
        ! iVertex.

        ! Step 1: Delete corresponding vertex from map. Note that for
        !         dense graphs the vertex map does not store further
        !         information and can be eliminated in the first
        !         step. If the vertex should be replaced by the last
        !         vertex, then the last vertex is removed from the map
        !         instead of vertex numbered iVertex.
        rmapIter = map_find(rgraph%rVertices,&
                            merge(ireplaceVertex,iVertex,bdoReplace))
        if (map_isNull(rmapIter)) then
          call output_line('Vertex does not exist in graph!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
          call sys_halt()
        else
          call map_erase(rgraph%rVertices, rmapIter)
        end if

        ! Step 2: Loop through adjacency list of vertex iVertex and
        !         delete all edges (jVertex,iVertex) from the
        !         adjacency list of jVertex.

        ! Find position of first entry in adjacency list
        ralstIter = alst_begin(rgraph%rEdges,iVertex)
        do while(.not.alst_isNull(ralstIter))

          ! Get number of adjacent vertex
          jVertex = alst_get(rgraph%rEdges,ralstIter)

          ! Get position of next entry in adjacency list
          call alst_next(ralstIter)

          ! Do nothing if both vertices are the same. The adjacency
          ! list of iVertex is removed/replaced anyway so we can leave
          ! it 'as is'
          if (iVertex .eq. jVertex) cycle

          ! Remove vertex iVertex from adjacency list of vertex jVertex
          rposition = alst_find(rgraph%rEdges,jVertex,iVertex)
          if (alst_isNull(rposition)) then
            call output_line('Unable to delete vertex from adjacency list!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            call sys_halt()
          else
            rposition = alst_erase(rgraph%rEdges,rposition)
          end if

          ! Decrease number of edges by two; for the edge
          ! (iVertex,jVertex) and for the edge (jVertex,iVertex) that
          ! exists in an undirected graph
          rgraph%NEDGE = rgraph%NEDGE-1
        end do

        ! Now, vertex iVertex does no longer exist in any adjacency list.
        ! Check if replacement vertex is different from iVertex.
        if (bdoReplace .and. (iVertex .ne. ireplaceVertex)) then

          ! Remove the trivial edge (ireplaceVertex,ireplaceVertex)
          rposition = alst_find(rgraph%rEdges,ireplaceVertex,ireplaceVertex)
          if (alst_isNull(rposition)) then
            call output_line('Unable to delete vertex from adjacency list!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            call sys_halt()
          else
            rposition = alst_erase(rgraph%rEdges,rposition)
          end if

          ! Swap adjacency list of vertices iVertex and ireplaceVertex
          ! From now onward, the adjacencey list of the replacement
          ! vertex is stored at position iVertex and vice versa. Keep
          ! this in mind!
          call alst_swapTbl(rgraph%rEdges,iVertex,ireplaceVertex)

          ! Release adjacency list of vertex ireplaceVertex
          call alst_releaseTbl(rgraph%rEdges,ireplaceVertex)

          ! Insert trivial edge (iVertex,iVertex) into adjacency list
          rposition = alst_find(rgraph%rEdges,iVertex,iVertex,.false.)
          if (alst_isNull(rposition)) then
            rposition = alst_insert(rgraph%rEdges,rposition,iVertex)
          else
            call output_line('Vertex already exists in adjacency list!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            call sys_halt()
          end if

          ! Step 3(a): Loop through adjacency list of vertex
          !            ireplaceVertex (which is already stored at its
          !            new position iVertex) and delete all edges
          !            (jVertex,ireplaceVertex) from the adjacency
          !            lists of jVertex. Afterwards, add the edge
          !            (jVertex,iVertex) to the adjacency list of
          !            jVertex.

          ! Find position of first entry in adjacency list
          ralstIter = alst_begin(rgraph%rEdges,iVertex)
          do while(.not.alst_isNull(ralstIter))

            ! Get number of adjacent vertex
            jVertex = alst_get(rgraph%rEdges,ralstIter)

            ! Get position of next entry in adjacency list
            call alst_next(ralstIter)

            ! Do nothing if jVertex is identical iVertex
            if (jVertex .eq. iVertex) cycle

            ! Remove vertex ireplaceVertex from adjacency list of
            ! vertex jVertex
            rposition = alst_find(rgraph%rEdges,jVertex,ireplaceVertex)
            if (alst_isNull(rposition)) then
              call output_line('Unable to delete vertex from adjacency list!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
              call sys_halt()
            else
              rposition = alst_erase(rgraph%rEdges,rposition)
            end if

            ! Insert edge (jVertex,iVertex) into adjacency list
            rposition = alst_find(rgraph%rEdges,jVertex,iVertex,.false.)
            if (alst_isNull(rposition)) then
              rposition = alst_insert(rgraph%rEdges,rposition,iVertex)
            else
              call output_line('Vertex already exists in adjacency list!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
              call sys_halt()
            end if

          end do

          ! Decrease number of vertices by one
          rgraph%NVT = rgraph%NVT-1

          ! Decrease number of edges by one; for the trivial edge
          ! (lastVertex,lastVertex)
          rgraph%NEDGE = rgraph%NEDGE-1

        else

          ! Step 3(b): Release adjacency list of vertex iVertex
          call alst_releaseTbl(rgraph%rEdges,iVertex)

          ! Decrease number of vertices by one
          rgraph%NVT = rgraph%NVT-1

          ! Decrease number of edges by one; for the trivial edge
          ! (iVertex,iVertex)
          rgraph%NEDGE = rgraph%NEDGE-1
        end if

      case default
        call output_line('Invalid graph format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
        call sys_halt()
      end select

    else   ! - - - - - - The following is for non-dense graphs - - - - - -

      ! What graph format are ?
      select case(rgraph%cgraphFormat)

      case(GRPH_GRAPHUNORDERED_DIRECTED,&
           GRPH_GRAPHORDERED_DIRECTED)
        call output_line('We cannot remove directed graphs at the moment!',&
            OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
        call sys_halt()

      case(GRPH_UNDIRECTED)
        ! We are in the lucky position, that the graph is undirected,
        ! that is, there exits an edge (i,j) if and only if there
        ! exists the edge (j,i).  Hence, we do not have to loop
        ! through the list of all vertices but it suffices to visit
        ! those which are present in the adjacency list of vertex
        ! iVertex.

        ! Get table for iVertex
        rmapIter = map_find(rgraph%rVertices,iVertex)
        if (map_isNull(rmapIter)) then
          call output_line('Vertex does not exist in graph!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
          call sys_halt()
        else
          call map_getbase_data(rgraph%rVertices,rmapIter,p_Idata)
          itable = p_Idata(1)
        end if

        ! Get table for ireplaceVertex (if required)
        if (bdoReplace) then
          rmapIter = map_find(rgraph%rVertices,ireplaceVertex)
          if (map_isNull(rmapIter)) then
            call output_line('Vertex does not exist in graph!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            call sys_halt()
          else
            call map_getbase_data(rgraph%rVertices,rmapIter,p_Idata)
            ireplaceTable = p_Idata(1)
          end if
        else
          ireplaceTable = 0
        end if

        ! Step 1: Delete corresponding vertex from map. Note that for
        !         dense graphs the vertex map does not store further
        !         information and can be eliminated in the first
        !         step. If the vertex should be replaced by the last
        !         vertex, then the last vertex is removed from the map
        !         instead of vertex numbered iVertex.
        rmapIter = map_find(rgraph%rVertices,&
                            merge(ireplaceTable,iVertex,bdoReplace))
        if (map_isNull(rmapIter)) then
          call output_line('Vertex does not exist in graph!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
          call sys_halt()
        else
          call map_erase(rgraph%rVertices,rmapIter)
        end if

        ! Step 2: Loop through adjacency list of vertex iVertex and
        !         delete all edges (jVertex,iVertex) from the
        !         adjacency lists of jVertex.

        ! Find position of first entry in adjacency list
        ralstIter = alst_begin(rgraph%rEdges,itable)
        do while(.not.alst_isNull(ralstIter))

          ! Get number of adjacent vertex
          jVertex = alst_get(rgraph%rEdges,ralstIter)

          ! Get position of next entry in adjacency list
          call alst_next(ralstIter)

          ! Do nothing if both vertices are the same
          if (iVertex .eq. jVertex) cycle

          ! In addition, do nothing if the current vertex is identical
          ! to the replacement vertex. Otherwise, we would have to
          ! re-insert it afterwards. Hence, it does not make sense to
          ! remove before.
          if (bdoReplace .and. (ireplaceVertex .eq. jVertex)) cycle

          ! Get table for jVertex
          rmapIter = map_find(rgraph%rVertices,jVertex)
          if (map_isNull(rmapIter)) then
            call output_line('Vertex does not exist in graph!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            call sys_halt()
          else
            call map_getbase_data(rgraph%rVertices,rmapIter,p_Idata)
            jtable = p_Idata(1)
          end if

          ! Remove vertex iVertex from adjacency list of vertex jVertex
          rposition = alst_find(rgraph%rEdges,jtable,iVertex)
          if (alst_isNull(rposition)) then
            call output_line('Unable to delete vertex from adjacency list!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            call sys_halt()
          else
            rposition = alst_erase(rgraph%rEdges, rposition)
          end if

          ! Decrease number of edges by two; for the edge
          ! (iVertex,jVertex) and for the edge (jVertex,iVertex) that
          ! exists in an undirected graph
          rgraph%NEDGE = rgraph%NEDGE-2
        end do

        ! Now, vertex iVertex does no longer exist in any adjacency
        ! list.  Check if replacement vertex needs to be moved to
        ! position iVertex.
        if (bdoReplace .and. (iVertex .ne. ireplaceVertex)) then

          ! Step 3(a): Loop through adjacency list of vertex
          !            ireplaceVertex and delete all edges
          !            (jVertex,ireplaceVertex) from the adjacency
          !            lists of jVertex. Afterwards, add the edge
          !            (jVertex,iVertex) to the adjacency list of
          !            jVertex. Finally, swap adjacency list of
          !            vertices iVertex and ireplaceVertex.

          ! Find position of first entry in adjacency list
          ralstIter = alst_begin(rgraph%rEdges,ireplaceTable)
          do while(.not.alst_isNull(ralstIter))

            ! Get number of adjacent vertex
            jVertex = alst_get(rgraph%rEdges,ralstIter)

            ! Get position of next entry in adjacency list
            call alst_next(ralstIter)

            ! Do nothing if both vertices are the same. This situation
            ! required special treatment (see below)
            if ((ireplaceVertex .eq. jVertex) .or. iVertex .eq. jVertex) cycle

            ! Get table for jVertex
            rmapIter = map_find(rgraph%rVertices,jVertex)
            if (map_isNull(rmapIter)) then
              call output_line('Vertex does not exist in graph!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
              call sys_halt()
            else
              call map_getbase_data(rgraph%rVertices,rmapIter,p_Idata)
              jtable = p_Idata(1)
            end if

            ! Remove vertex lastVertex from adjacency list of vertex jVertex
            rposition = alst_find(rgraph%rEdges,jtable,ireplaceVertex)
            if (alst_isNull(rposition)) then
              call output_line('Unable to update vertex in adjacency list!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
              call sys_halt()
            else
              rposition = alst_erase(rgraph%rEdges,rposition)
            end if

            ! Insert edge (jVertex,iVertex) into adjacency list
            rposition = alst_find(rgraph%rEdges,jtable,iVertex,.false.)
            if (alst_isNull(rposition)) then
              rposition = alst_insert(rgraph%rEdges,rposition,iVertex)
            else
              call output_line('Vertex replacement already exists in adjacency list!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
              call sys_halt()
            end if
          end do

          ! Remove the trivial edge (ireplaceVertex,ireplaceVertex)
          ! from the adjacency list. Note that the trivial edge
          ! (iVertex,iVertex) still exists and has not been removed in
          ! the upper removal loop
          rposition = alst_find(rgraph%rEdges,ireplaceTable,ireplaceVertex)
          if (alst_isNull(rposition)) then
            call output_line('Unable to update vertex in adjacency list!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            call sys_halt()
          else
            rposition = alst_erase(rgraph%rEdges,rposition)
          end if

          ! Swap adjacency list of vertices iVertex and ireplaceVertex
          call alst_swapTbl(rgraph%rEdges,itable,ireplaceTable)

          ! Release adjacency list of vertex ireplaceVertex
          call alst_releaseTbl(rgraph%rEdges,ireplaceTable)

          ! Decrease number of vertices by one
          rgraph%NVT = rgraph%NVT-1

          ! Decrease number of edges by one; for the trivial edge
          ! (lastVertex,lastVertex)
          rgraph%NEDGE = rgraph%NEDGE-1

        else

          ! Step 3(b): Release adjacency list of vertex iVertex
          call alst_releaseTbl(rgraph%rEdges,itable)

          ! Decrease number of vertices by one
          rgraph%NVT = rgraph%NVT-1

          ! Decrease number of edges by one; for the trivial edge
          ! (iVertex,iVertex)
          rgraph%NEDGE = rgraph%NEDGE-1
        end if

      case default
        call output_line('Invalid graph format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
        call sys_halt()
      end select
    end if

  end subroutine grph_removeVertex

  ! ***************************************************************************

!<function>

  function grph_hasEdge(rgraph,iFromVertex,iToVertex) result(bexists)

!<description>
    ! This function returns TRUE if the graph has the given edge from
    ! vertex iFromVertex to vertex iToVertex.  Otherwise, it returns
    ! FALSE. If the optional parameter iEdgePosition is present, then
    ! the function will also return the position of the edge
    ! (iFromVertex,iToVertex) in the array list.
!</description>

!<input>
    ! Number of the starting vertex
    integer, intent(in) :: iFromVertex

    ! Number of the ending vertex
    integer, intent(in) :: iToVertex
!</input>

!<inputoutput>
    ! The graph
    type(t_graph), intent(inout) :: rgraph
!</inputoutput>

!<result>
    ! Flag for existence of the vertex
    logical :: bexists
!</result>
!</function>

    ! local variables
    type(it_mapInt_Int) :: rmapIter
    integer, dimension(:), pointer :: p_Idata
    integer :: itable

    ! Determine table
    if (rgraph%bisDense) then
      itable = iFromVertex
    else
      rmapIter = map_find(rgraph%rVertices,iFromVertex)
      if (map_isNull(rmapIter)) then
        bexists = .false.
        return
      else
        call map_getbase_data(rgraph%rVertices,rmapIter,p_Idata)
        itable = p_Idata(1)
      end if
    end if

    ! Look for edge (iFromVertex,iToVertex)
    bexists = .not.alst_isNull(alst_find(rgraph%rEdges,itable,iToVertex))

  end function grph_hasEdge

  ! ***************************************************************************

!<subroutine>

  subroutine grph_insertEdge(rgraph,iFromVertex,iToVertex)

!<description>
    ! This subroutine inserts the edge between vertices iFromVertex
    ! and iToVertex. If the graph is undirected, then the opposite
    ! edge (iToVertex,iFromVertex) is also added to attain symmetry of
    ! the sparsity pattern.
!</description>

!<input>
    ! Number of the starting vertex
    integer, intent(in) :: iFromVertex

    ! Number of the ending vertex
    integer, intent(in) :: iToVertex
!</input>

!<inputoutput>
    ! The graph
    type(t_graph), intent(inout) :: rgraph
!</inputoutput>
!</subroutine>

    ! local variables
    type(it_mapInt_Int) :: rmapIter
    type(it_arraylistInt) :: ralstIter
    integer, dimension(:), pointer :: p_Idata
    integer :: itable

    ! What graph format are we?
    select case(rgraph%cgraphFormat)

    case(GRPH_GRAPHUNORDERED_DIRECTED)

      if (rgraph%bisDense) then

        ! Insert entry iToVertex into adjacency list of vertex iFromVertex
        ralstIter = alst_find(rgraph%rEdges,iToVertex,iFromVertex)
        if (alst_isNull(ralstIter)) then
          call alst_push_back(rgraph%rEdges,iToVertex,iFromVertex)
          rgraph%NEDGE = rgraph%NEDGE+1
        end if

      else

        ! Get table associated with vertex iToVertex
        rmapIter = map_find(rgraph%rVertices,iToVertex)
        if (map_isNull(rmapIter)) then
          call output_line('Vertex does not exists!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_insertEdge')
          call sys_halt()
        else
          call map_getbase_data(rgraph%rVertices,rmapIter,p_Idata)
          itable = p_Idata(1)
        end if

        ! Insert entry iToVertex into adjacency list of vertex iFromVertex
        ralstIter = alst_find(rgraph%rEdges,itable,iFromVertex)
        if (alst_isNull(ralstIter)) then
          call alst_push_back(rgraph%rEdges,itable,iFromVertex)
          rgraph%NEDGE = rgraph%NEDGE+1
        end if

      end if


    case(GRPH_GRAPHORDERED_DIRECTED)

      if (rgraph%bisDense) then

        ! Insert entry iToVertex into adjacency list of vertex iFromVertex
        ralstIter = alst_find(rgraph%rEdges,iToVertex,iFromVertex,.false.)
        if (alst_isNull(ralstIter)) then
          ralstIter = alst_insert(rgraph%rEdges,ralstIter,iFromVertex)
          rgraph%NEDGE = rgraph%NEDGE+1
        end if

      else

        ! Get table associated with vertex iToVertex
        rmapIter = map_find(rgraph%rVertices,iToVertex)
        if (map_isNull(rmapIter)) then
          call output_line('Vertex does not exists!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_insertEdge')
          call sys_halt()
        else
          call map_getbase_data(rgraph%rVertices,rmapIter,p_Idata)
          itable = p_Idata(1)
        end if

        ! Insert entry iToVertex into adjacency list of vertex iFromVertex
        ralstIter = alst_find(rgraph%rEdges,itable,iFromVertex,.false.)
        if (alst_isNull(ralstIter)) then
          ralstIter = alst_insert(rgraph%rEdges,ralstIter,iFromVertex)
          rgraph%NEDGE = rgraph%NEDGE+1
        end if
      end if


    case(GRPH_UNDIRECTED)

      if (rgraph%bisDense) then

        ! Insert entry iFromVertex into adjacency list of vertex iToVertex
        ralstIter = alst_find(rgraph%rEdges,iToVertex,iFromVertex,.false.)
        if (alst_isNull(ralstIter)) then
          ralstIter = alst_insert(rgraph%rEdges,ralstIter,iFromVertex)
          rgraph%NEDGE = rgraph%NEDGE+1
        end if

        ! Insert entry iToVertex into adjacency list of vertex iFromVertex
        ralstIter = alst_find(rgraph%rEdges,iFromVertex,iToVertex,.false.)
        if (alst_isNull(ralstIter)) then
          ralstIter = alst_insert(rgraph%rEdges,ralstIter,iToVertex)
          rgraph%NEDGE = rgraph%NEDGE+1
        end if

      else

        ! Get table associated with vertex iToVertex
        rmapIter = map_find(rgraph%rVertices,iToVertex)
        if (map_isNull(rmapIter)) then
          call output_line('Vertex does not exists!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_insertEdge')
          call sys_halt()
        else
          call map_getbase_data(rgraph%rVertices,rmapIter,p_Idata)
          itable = p_Idata(1)
        end if

        ! Insert entry iFromVertex into adjacency list of vertex iToVertex
        ralstIter = alst_find(rgraph%rEdges,itable,iFromVertex,.false.)
        if (alst_isNull(ralstIter)) then
          ralstIter = alst_insert(rgraph%rEdges,ralstIter,iFromVertex)
          rgraph%NEDGE = rgraph%NEDGE+1
        end if

        ! Get table associated with vertex iFromVertex
        rmapIter = map_find(rgraph%rVertices,iFromVertex)
        if (map_isNull(rmapIter)) then
          call output_line('Vertex does not exists!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_insertEdge')
          call sys_halt()
        else
          call map_getbase_data(rgraph%rVertices,rmapIter,p_Idata)
          itable = p_Idata(1)
        end if

        ! Insert entry iToVertex into adjacency list of vertex iFromVertex
        ralstIter = alst_find(rgraph%rEdges,itable,iToVertex,.false.)
        if (alst_isNull(ralstIter)) then
          ralstIter = alst_insert(rgraph%rEdges,ralstIter,iToVertex)
          rgraph%NEDGE = rgraph%NEDGE+1
        end if

      end if


    case default
      call output_line('Unsupported graph format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'grph_insertEdge')
      call sys_halt()
    end select

  end subroutine grph_insertEdge

  ! ***************************************************************************

!<subroutine>

  subroutine grph_removeEdge(rgraph,iFromVertex,iToVertex)

!<description>
    ! This subroutine removes the edge between vertices iFromVertex
    ! and iToVertex. If the graph is undirected, then the opposite
    ! edge (iToVertex,iFromVertex) is also removed to attain symmetry
    ! of the sparsity pattern.
!</description>

!<input>
    ! Number of the starting vertex
    integer, intent(in) :: iFromVertex

    ! Number of the ending vertex
    integer, intent(in) :: iToVertex
!</input>

!<inputoutput>
    ! The graph
    type(t_graph), intent(inout) :: rgraph
!</inputoutput>
!</subroutine>

    ! local variables
    type(it_mapInt_Int) :: rmapIter
    type(it_arraylistInt) :: ralstIter
    integer, dimension(:), pointer :: p_Idata
    integer :: itable

    ! What graph format are we=
    select case(rgraph%cgraphFormat)

    case(GRPH_GRAPHUNORDERED_DIRECTED,&
         GRPH_GRAPHORDERED_DIRECTED)

      if (rgraph%bisDense) then

        ! Remove vertex iToVertex from table iFromVertex
        ralstIter = alst_find(rgraph%rEdges,iFromVertex,iToVertex)
        if (.not.alst_isNull(ralstIter)) then
          ralstIter = alst_erase(rgraph%rEdges,ralstIter)
          rgraph%NEDGE = rgraph%NEDGE-1
        end if

      else

        ! Get table associated with vertex iFromVertex
        rmapIter = map_find(rgraph%rVertices,iFromVertex)
        if (map_isNull(rmapIter)) then
          call output_line('Vertex does not exist!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_removeEdge')
          call sys_halt()
        else
          call map_getbase_data(rgraph%rVertices,rmapIter,p_Idata)
          itable = p_Idata(1)
        end if

        ! Remove vertex iToVertex from table iFromVertex
        ralstIter = alst_find(rgraph%rEdges,itable,iToVertex)
        if (.not.alst_isNull(ralstIter)) then
          ralstIter = alst_erase(rgraph%rEdges,ralstIter)
          rgraph%NEDGE = rgraph%NEDGE-1
        end if
      end if


    case(GRPH_UNDIRECTED)

      if (rgraph%bisDense) then

        ! Remove vertex iToVertex from table iFromVertex
        ralstIter = alst_find(rgraph%rEdges,iFromVertex,iToVertex)
        if (.not.alst_isNull(ralstIter)) then
          ralstIter = alst_erase(rgraph%rEdges,ralstIter)
          rgraph%NEDGE = rgraph%NEDGE-1
        end if

        ! Remove vertex iFromVertex from table iToVertex
        ralstIter = alst_find(rgraph%rEdges,iToVertex,iFromVertex)
        if (.not.alst_isNull(ralstIter)) then
          ralstIter = alst_erase(rgraph%rEdges,ralstIter)
          rgraph%NEDGE = rgraph%NEDGE-1
        end if

      else

        ! Get table associated with vertex iFromVertex
        rmapIter = map_find(rgraph%rVertices,iFromVertex)
        if (map_isNull(rmapIter)) then
          call output_line('Vertex does not exist!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_removeEdge')
          call sys_halt()
        else
          call map_getbase_data(rgraph%rVertices,rmapIter,p_Idata)
          itable = p_Idata(1)
        end if

        ! Remove vertex iToVertex from table iFromVertex
        ralstIter = alst_find(rgraph%rEdges,itable,iToVertex)
        if (.not.alst_isNull(ralstIter)) then
          ralstIter = alst_erase(rgraph%rEdges,ralstIter)
          rgraph%NEDGE = rgraph%NEDGE-1
        end if

        ! Get table associated with vertex iToVertex
        rmapIter = map_find(rgraph%rVertices,iToVertex)
        if (map_isNull(rmapIter)) then
          call output_line('Vertex does not exist!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_removeEdge')
          call sys_halt()
        else
          call map_getbase_data(rgraph%rVertices,rmapIter,p_Idata)
          itable = p_Idata(1)
        end if

        ! Remove vertex iFromVertex from table iToVertex
        ralstIter = alst_find(rgraph%rEdges,itable,iFromVertex)
        if (.not.alst_isNull(ralstIter)) then
          ralstIter = alst_erase(rgraph%rEdges,ralstIter)
          rgraph%NEDGE = rgraph%NEDGE-1
        end if
      end if


    case default
      call output_line('Unsupported graph format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'grph_removeEdge')
      call sys_halt()
    end select

  end subroutine grph_removeEdge

  ! ***************************************************************************

!<subroutine>

  subroutine grph_infoGraph(rgraph)

!<description>
    ! This subroutine prints information about the grahp
!</description>

!<input>
    ! graph
    type(t_graph), intent(in) :: rgraph
!</input>
!</subroutine>

    call output_line('Graph statistics:')
    call output_line('-----------------')
    call output_line('cgraphFormat: '//trim(sys_siL(rgraph%cgraphFormat,2)))
    call output_line('bisDense:     '//merge('Yes','No ',rgraph%bisDense))
    call output_line('NVT:          '//trim(sys_siL(rgraph%NVT,15)))
    call output_line('NEDGE:        '//trim(sys_siL(rgraph%NEDGE,15)))
    call output_lbrk()

  end subroutine grph_infoGraph

  ! ***************************************************************************

!<subroutine>

  subroutine grph_duplicateGraph(rgraph,rgraphBackup)

!<description>
    ! This subroutine makes a copy of a graph in memory. It does not
    ! make sense to share some information between graphs, so each
    ! vectors is physically copied from the source graph to the
    ! destination graph.
!</description>

!<input>
    ! Source graph
    type(t_graph), intent(in) :: rgraph
!</input>

!<inputoutput>
    ! Destination graph
    type(t_graph), intent(inout) :: rgraphBackup
!</inputoutput>
!</subroutine>

    ! Release backup graph
    call grph_releaseGraph(rgraphBackup)

    ! Copy all data
    rgraphBackup = rgraph

    call map_duplicate(rgraph%rVertices,rgraphBackup%rVertices)
    call alst_duplicate(rgraph%rEdges,rgraphBackup%rEdges)

  end subroutine grph_duplicateGraph

  ! ***************************************************************************

!<subroutine>

  subroutine grph_restoreGraph(rgraphBackup,rgraph)

!<description>
    ! This subroutine restores a graph from a previous backup.
    ! The format of both graphs must be the same.
!</description>

!<input>
    ! Backup of a graph
    type(t_graph), intent(in) :: rgraphBackup
!</input>

!<inputoutput>
    ! Destination graph
    type(t_graph), intent(inout) :: rgraph
!</inputoutput>
!</subroutine>

    ! Check that both graphs are compatible
    if (rgraph%cgraphFormat .ne. rgraphBackup%cgraphFormat .or.&
        rgraph%bisDense     .neqv. rgraphBackup%bisDense) then
      call output_line('Incompatible graphs!',&
          OU_CLASS_ERROR,OU_MODE_STD,'grph_restoreGraph')
      call sys_halt()
    end if

    ! Release graph
    call grph_releaseGraph(rgraph)

    ! Duplicate the backup
    call grph_duplicateGraph(rgraphBackup,rgraph)
  end subroutine grph_restoreGraph

end module graph
