!##############################################################################
!# ****************************************************************************
!# <name> graph </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic definitions and routines to maintain
!# graphs. A graph consists of a set of vertices V and a set of edged E
!# connecting vertices. Graphs are stored as adjacency matrices which can be
!# modified dynamically. A new vertex can be added or an existing vertex can
!# be removed. Moreover, edges can be inserted between two vertices or removed.
!# This module provides data structures for both directed and undirected graphs.
!# You should be aware of the fact, that adjacency lists are perhaps not the
!# optimal choice for nearly complete graphs, that is, graphs which possess
!# the majority of all possible edges. In this case it would make sense to
!# adopt a static matrix data structure. Be warned, this is not implemented
!# in this module so that the performance of the provided routines may
!# deteriorate for nearly complete graphs.
!#
!# To keep things simple, graphs can be seen as sparce matrices stored in the
!# CSR format (either 7 or 9) which can be modified dynamically. Hence, it
!# makes sence to provide conversion routines between a scalar matrix and a
!# graph. If you are familiar with the CSR format used for scalar matrices
!# then this hint will be helpful: The list of vertices is stored in the Kld
!# array whereas the edges are associated with the Kcol array. In short, an
!# edge (I,J) exists if and only if there is an entry J in the subarray
!# Kcol(Kld(i):Kld(i+1)-1).
!#
!#
!#
!# The following routines can be found in this module:
!#
!#  1.) grph_createGraph
!#      -> Create graph directly
!#
!#  2.) grph_createGraphFromMatrix
!#      -> Create graph from scalar matrix
!#
!#  3.) grph_releaseGraph
!#      -> Release an existing graph
!#
!#  4.) grph_generateMatrix
!#      -> Generate scalar matrix from sparsity graph
!#
!#  5.) grph_printGraph
!#      -> Print the graph to screen
!#
!#  6.) grph_hasVertex
!#      -> Check if the graph has the given vertex
!#
!#  7.) grph_insertVertex
!#      -> Insert vertex into the graph. Do nothing if vertex already exists.
!#
!#  8.) grph_removeVertex
!#      -> Remove vertex from graph. Do nothing if vertex does not exist.
!#
!#  9.) grph_hasEdge
!#      -> Check if the graph has the given edge
!#
!# 10.) grph_insertEdge
!#      -> Insert edge into the graph. Do nothing if edge already exists.
!#
!# 11.) grph_removeEdge
!#      -> Remove edge from graph. Do nothing if edge does not exist.
!#
!# 12.) grph_infoGraph
!#      -> Print information about the graph.
!#
!# 13.) grph_duplicateGraph
!#      -> Create a duplicate / backup of a graph.
!#
!# 14.) grph_restoreGraph
!#      -> Restores a graph previously backed up with grph_duplicateGraph
!#
!# </purpose>
!##############################################################################

module graph

  use arraylist
  use binarytree
  use fsystem
  use genoutput
  use linearsystemscalar
  use storage

  implicit none

  private
  public :: t_graph
  public :: grph_createGraph
  public :: grph_createGraphFromMatrix
  public :: grph_releaseGraph
  public :: grph_generateMatrix
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

  ! Identifier for undirected graph with ordered entries in CSR format 9
  integer, parameter, public :: GRPH_GRAPH9 = 9

  ! Identifier for undirected graph with ordered entries in CSR format 7,
  ! that is, CSR with diagonal element in front
  integer, parameter, public :: GRPH_GRAPH7 = 7

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

    ! Is set to true, if the graph is dense, that is, a graph with NVT vertices
    ! contains the mainly the vertices labeled 1,...,NVT. In this case the vertex
    ! number coincides with the table under which its adjacency list is stored.
    ! If the graph is note marked as dense, then for each vertex IVT the table
    ! number is stored in the vertex list and must be looked-up each time, the
    ! adjacency information is required.
    ! Note that dense graphs can only store positive vertex labels !!!
    logical :: bisDense = .true.

    ! Search tree for the list of vertices
    type(t_btree) :: rVertices

    ! Array of lists for the edges
    type(t_arraylist) :: rEdges

  end type t_graph
    
!</typeblock>
!</types>

contains

!<subroutine>

  subroutine grph_createGraph(rgraph,cgraphFormat,nvtMax,nedgeMax,bisDense)

!<description>
    ! This routine creates a graph rgraph. The graph format must be one of
    ! the GRPH_GRAPHx format specifiers. The parameters nvtMax and
    ! nedgeMax define the maximum number of vertices/egdes that can be
    ! initially stored in the graph. If either more vertices or edges are
    ! inserted into the graph, then reallocation of memory will be performed.
    ! The optional parameter bisDense defines, if the graph should be treted
    ! as dense graph or not. In a dense graph, the adjacency list of vertex
    ! ivt is stored in table itable=ivt. It makes sense to use dense graphs,
    ! e.g., if NVT vertices numbered by 1..NVT are present in the graph so
    ! that the table look-up can be performed in time O(1). If the graph
    ! consists of NVT vertices which can be labeled arbitrarily, e.g.,
    ! <tex>1.. ($1^k$)*NVT, k >> 1</tex>, then an enormous amount of memory would be
    ! wasted for the look-up tables. Then it makes sense to define the graph as
    ! "not dense". In this case, the tables are used sequentially if new
    ! vertices are inserted but finding the correct table number for vertex
    ! ivt requires a search in the binary tree, i.e. time O(log NVT).
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

    ! What graph format are we?
    select case(cgraphFormat)

    case(GRPH_GRAPHORDERED_DIRECTED)
      rgraph%cgraphFormat = GRPH_GRAPHORDERED_DIRECTED

      ! Create search tree for the list of vertices
      if (rgraph%bisDense) then
        call btree_createTree(rgraph%rVertices,rgraph%NVT+1,ST_INT,0,0,0)
      else
        call btree_createTree(rgraph%rVertices,rgraph%NVT+1,ST_INT,0,0,1)
      end if
      
      ! Create array of lists for edges
      call arrlst_createArrayList(rgraph%rEdges,nvtMax,nedgeMax,ST_INT,ARRAYLIST_INCREASING)

    case(GRPH_GRAPHUNORDERED_DIRECTED)
      rgraph%cgraphFormat = GRPH_GRAPHUNORDERED_DIRECTED

      ! Create search tree for the list of vertices
      if (rgraph%bisDense) then
        call btree_createTree(rgraph%rVertices,rgraph%NVT+1,ST_INT,0,0,0)
      else
        call btree_createTree(rgraph%rVertices,rgraph%NVT+1,ST_INT,0,0,1)
      end if
      
      ! Create array of lists for edges
      call arrlst_createArrayList(rgraph%rEdges,nvtMax,nedgeMax,ST_INT,ARRAYLIST_UNORDERED)

    case(GRPH_GRAPH7)
      rgraph%cgraphFormat = GRPH_GRAPH7

      ! Create search tree for the list of vertices
      if (rgraph%bisDense) then
        call btree_createTree(rgraph%rVertices,rgraph%NVT+1,ST_INT,0,0,0)
      else
        call btree_createTree(rgraph%rVertices,rgraph%NVT+1,ST_INT,0,0,1)
      end if

      ! Create array of lists for edges
      call arrlst_createArrayList(rgraph%rEdges,nvtMax,nedgeMax,ST_INT,ARRAYLIST_CSR7)

    case(GRPH_GRAPH9)
      rgraph%cgraphFormat = GRPH_GRAPH9

      ! Create search tree for the list of vertices
      if (rgraph%bisDense) then
        call btree_createTree(rgraph%rVertices,rgraph%NVT+1,ST_INT,0,0,0)
      else
        call btree_createTree(rgraph%rVertices,rgraph%NVT+1,ST_INT,0,0,1)
      end if

      ! Create array of lists for edges
      call arrlst_createArrayList(rgraph%rEdges,nvtMax,nedgeMax,ST_INT,ARRAYLIST_INCREASING)
      
    case DEFAULT
      call output_line('Invalid matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'grph_createGraph')
      call sys_halt()
    end select
  end subroutine grph_createGraph

  ! ***************************************************************************

!<subroutine>

  subroutine grph_createGraphFromMatrix(rscalarMatrix,rgraph,nvtMax,nedgeMax)

!<description>
    ! This routine creates a graph rgraph from a scalar matrix rscalarMatrix.
    ! The matrix must be stored in CSR format 7 or 9. Note that a graph that
    ! is generated from a matrix is assumed to be dense automatically.
    ! If the optional parameters nvtMax or nedgeMax are present then
    ! their values will be adopted for the maximum number of vertices/edges
    ! than can be stored in the graph initially without reallocation.
!</description>

!<input>
    ! The scalar matrix which should be used to create the sparsity graph
    type(t_matrixScalar), intent(in) :: rscalarMatrix

    ! OPTIONAL: maximum number of vertices
    integer, intent(in), optional :: nvtMax

    ! OPTIONAL: maximum number of edges
    integer, intent(in), optional :: nedgeMax
!</input>

!<output>
    ! The sparsity graph
    type(t_graph), intent(out) :: rgraph
!</output>

!</subroutine>

    ! local variables
    integer :: h_Key
    integer :: nvt,nedge
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Key

    ! Set dimensions
    rgraph%NVT   = rscalarMatrix%NEQ
    rgraph%NEDGE = rscalarMatrix%NA

    ! Estimate maximum number of vertices
    if (present(nvtMax)) then
      nvt = max(nvtMax,rgraph%NVT)
    else
      nvt = int(1.5_DP*rgraph%NVT)
    end if

    ! Estimate maximum number of edges
    if (present(nedgeMax)) then
      nedge = max(nedgeMax,rgraph%NEDGE)
    else
      nedge = int(0.5_DP*(rgraph%NEDGE+rgraph%NVT))
    end if

    ! What matrix format are we?
    select case(rscalarMatrix%cmatrixFormat)

    case(LSYSSC_MATRIX7,&
         LSYSSC_MATRIX7INTL)
      rgraph%cgraphFormat = GRPH_GRAPH7

      ! Create search tree for the list of vertices
      call btree_createTree(rgraph%rVertices,rgraph%NVT+1,ST_INT,0,0,0)
      
      ! Create array of lists for edges
      call arrlst_createArrayList(rgraph%rEdges,nvt,nedge,ST_INT,ARRAYLIST_CSR7)

      ! Set pointers
      call lsyssc_getbase_Kld(rscalarMatrix,p_Kld)
      call lsyssc_getbase_Kcol(rscalarMatrix,p_Kcol)

      ! Fill list of edges
      call arrlst_copyArrayListTable(p_Kcol,rgraph%rEdges,p_Kld)

      ! Generate p_Key = array [1,2,3,...,NVT]
      call storage_new('grph_createGraphFromMatrix','p_Key',rgraph%NVT,ST_INT,&
          h_Key,ST_NEWBLOCK_ORDERED)
      call storage_getbase_int(h_Key,p_Key)
      call btree_copyToTree(p_Key,rgraph%rVertices)
      call storage_free(h_Key)
      

    case(LSYSSC_MATRIX9,&
         LSYSSC_MATRIX9INTL)
      rgraph%cgraphFormat = GRPH_GRAPH9
      
      ! Create search tree for the list of vertices
      call btree_createTree(rgraph%rVertices,rgraph%NVT+1,ST_INT,0,0,0)
      
      ! Create array of lists for edges
      call arrlst_createArrayList(rgraph%rEdges,nvt,nedge,ST_INT,ARRAYLIST_INCREASING)

      ! Set pointers
      call lsyssc_getbase_Kld(rscalarMatrix,p_Kld)
      call lsyssc_getbase_Kcol(rscalarMatrix,p_Kcol)

      ! Fill list of edges
      call arrlst_copyArrayListTable(p_Kcol,rgraph%rEdges,p_Kld)

      ! Generate p_Key = array [1,2,3,...,NVT]
      call storage_new('grph_createGraphFromMatrix','p_Key',rgraph%NVT,ST_INT,&
          h_Key,ST_NEWBLOCK_ORDERED)
      call storage_getbase_int(h_Key,p_Key)
      call btree_copyToTree(p_Key,rgraph%rVertices)
      call storage_free(h_Key)


    case DEFAULT
      call output_line('Invalid matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'grph_createGraphFromMatrix')
      call sys_halt()
    end select
  end subroutine grph_createGraphFromMatrix
  
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

    ! Release search tree for the list of vertices
    call btree_releaseTree(rgraph%rVertices)

    ! Release array of lists for edges
    call arrlst_releaseArrayList(rgraph%rEdges)

    ! Reset data
    rgraph%cgraphFormat = GRPH_GRAPHUNDEFINED
    rgraph%NVT          = 0
    rgraph%NEDGE        = 0
    rgraph%bisDense     = .true.
  end subroutine grph_releaseGraph

  ! ***************************************************************************

!<subroutine>

  subroutine grph_generateMatrix(rgraph,rscalarMatrix)

!<description>
    ! This subroutine generates a sparse matrix rscalarMatrix stored in
    ! format CSR 7 or 9 whereby the sparsity pattern is adopted from
    ! the graph rgraph.
!</description>

!<input>
    ! The graph
    type(t_graph), intent(inout) :: rgraph
!</input>

!<inputoutput>
    ! The matrix
    type(t_matrixScalar), intent(inout) :: rscalarMatrix
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Kdiagonal
    integer :: itable,ieq,ia,ncols,ipred,ipos
    integer :: isize
    
    ! Check that matrix and graph have the same format
    select case(rscalarMatrix%cmatrixFormat)
      
    case(LSYSSC_MATRIX7,&
         LSYSSC_MATRIX7INTL)
      if (rgraph%cgraphFormat .ne. GRPH_GRAPH7) then
        call output_line('Matrix/graph have incompatible format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'grph_generateMatrix')
        call sys_halt()
      end if

    case(LSYSSC_MATRIX9,&
         LSYSSC_MATRIX9INTL)
      if (rgraph%cgraphFormat .ne. GRPH_GRAPH9) then
        call output_line('Matrix/graph have incompatible format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'grph_generateMatrix')
        call sys_halt()
      end if

    case DEFAULT
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'grph_generateMatrix')
      call sys_halt()
    end select
 
    ! Set number of edges = number of nonzero matrix entries
    rscalarMatrix%NA = rgraph%NEDGE
    
    ! Set number of columns/rows. This is a little bit ugly because the
    ! number of vertices (NVT) may be different from the number of tables.
    ! If vertices are inserted as follows 1,2,3,...,NVT, then the number of
    ! tables is equal to the number of vertices. However, if only vertices
    ! 1,5,9,10 are present in the graph, then NVT=4 but the largest number in
    ! the graph is 10. Hence, the matrix has NEQ=NCOLS=10 and some rows, e.g.
    ! 2,3,4 have no entries. For a finite element matrix, this does not make
    ! sense but this graph module should be as general as possible.
    rscalarMatrix%NEQ   = max(rgraph%NVT,rgraph%rEdges%NTABLE)
    rscalarMatrix%NCOLS = rscalarMatrix%NEQ

    if (rgraph%bisDense) then
      
      ! Convert array list to matrix
      call arrlst_copyArrayListTable(rgraph%rEdges,rscalarMatrix%h_Kcol,rscalarMatrix%h_Kld)
      
    else

      ! Check if matrix is empty
      if ((rscalarMatrix%NEQ .eq. 0) .or. rscalarMatrix%NA .eq. 0) return

      ! Convert array list step-by-step
      if (rscalarMatrix%h_Kld .eq. ST_NOHANDLE) then
        call storage_new('grph_generateMatrix','p_Kld',&
            rscalarMatrix%NEQ+1,ST_INT,rscalarMatrix%h_Kld,ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(rscalarMatrix%h_Kld,isize)
        if (isize < rscalarMatrix%NEQ+1) then
          call storage_realloc('grph_generateMatrix',&
              rscalarMatrix%NEQ+1,rscalarMatrix%h_Kld,ST_NEWBLOCK_NOINIT,.false.)
        end if
      end if

      if (rscalarMatrix%h_Kcol .eq. ST_NOHANDLE) then
        call storage_new('grph_generateMatrix','p_Kcol',&
            rscalarMatrix%NA,ST_INT,rscalarMatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(rscalarMatrix%h_Kcol,isize)
        if (isize < rscalarMatrix%NA) then
          call storage_realloc('grph_generateMatrix',&
              rscalarMatrix%NA,rscalarMatrix%h_Kcol,ST_NEWBLOCK_NOINIT,.false.)
        end if
      end if

      ! Set pointers
      call lsyssc_getbase_Kld(rscalarMatrix,p_Kld)
      call lsyssc_getbase_Kcol(rscalarMatrix,p_Kcol)

      ! Initialization
      ia=1

      ! Loop over all equations
      do ieq=1,rscalarMatrix%NEQ

        p_Kld(ieq) = ia

        ! Check if vertex with number ieq exists
        if (btree_searchInTree(rgraph%rVertices,ieq,ipred).eq.BTREE_FOUND) then
          ipos = rgraph%rVertices%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))
          itable = rgraph%rVertices%p_IData(1,ipos)

          ! Restore row from table
          call arrlst_copyArrayList(rgraph%rEdges,itable,p_Kcol(ia:),ncols)
        end if

        ia = ia+ncols
      end do
      p_Kld(rscalarMatrix%NEQ+1)=ia
    end if

    ! Do we have to rebuild the diagonal?
    if (rscalarMatrix%cmatrixFormat .eq. LSYSSC_MATRIX9 .or.&
        rscalarMatrix%cmatrixFormat .eq. LSYSSC_MATRIX9INTL) then

      ! Create new memory or resize existing memory
      if (rscalarMatrix%h_Kdiagonal .eq. ST_NOHANDLE) then
        call storage_new('grph_generateMatrix','p_Kdiagonal',&
            rscalarMatrix%NEQ, ST_INT, rscalarMatrix%h_Kdiagonal,ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(rscalarMatrix%h_Kdiagonal,isize)
        if (isize < rscalarMatrix%NEQ) then
          call storage_realloc('grph_generateMatrix',&
              rscalarMatrix%NEQ, rscalarMatrix%h_Kdiagonal, ST_NEWBLOCK_NOINIT, .false.)
        end if
      end if

      ! Set pointers
      call lsyssc_getbase_Kld(rscalarMatrix, p_Kld)
      call lsyssc_getbase_Kcol(rscalarMatrix, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rscalarMatrix, p_Kdiagonal)

      ! Rebuild array
      call lsyssc_rebuildKdiagonal(p_Kcol, p_Kld, p_Kdiagonal, rscalarMatrix%NEQ)
    end if
  end subroutine grph_generateMatrix

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

    ! Are we dense graph?
    if (rgraph%bisDense) then
      call inorderDense(rgraph%rVertices%p_Kchild(TRIGHT,TROOT))
    else
      call inorderSparse(rgraph%rVertices%p_Kchild(TRIGHT,TROOT))
    end if

  contains
    
    ! Here, the real working routines follow.
    
    !**************************************************************
    ! Inorder traversal of the tree storing the vertex numbers
    recursive subroutine inorderDense(i)
      integer, intent(in) :: i

      ! local variables
      integer :: ikey

      ! Check if position is valid
      if (i .eq. TNULL) return

      ! Proceed with left child if it exists
      if (rgraph%rVertices%p_Kchild(TLEFT,i) .ne. TNULL)&
          call inorderDense(rgraph%rVertices%p_Kchild(TLEFT,i))

      ! We are in the lucky position that itable = ikey
      ikey   = rgraph%rVertices%p_IKey(i)
      
      call output_line('Vertex number: '//trim(sys_siL(ikey,15)))
      call arrlst_printArrayList(rgraph%rEdges,ikey)
        
      ! Proceed with right child if it exists
      if (rgraph%rVertices%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call inorderDense(rgraph%rVertices%p_Kchild(TRIGHT,i))
    end subroutine inorderDense

     !**************************************************************
    ! Inorder traversal of the tree storing the vertex numbers
    recursive subroutine inorderSparse(i)
      integer, intent(in) :: i

      ! local variables
      integer :: itable,ikey

      ! Check if position is valid
      if (i .eq. TNULL) return

      ! Proceed with left child if it exists
      if (rgraph%rVertices%p_Kchild(TLEFT,i) .ne. TNULL)&
          call inorderSparse(rgraph%rVertices%p_Kchild(TLEFT,i))

      ! In general ikey != itable
      ikey   = rgraph%rVertices%p_IKey(i)
      itable = rgraph%rVertices%p_IData(1,i)

      call output_line('Vertex number: '//trim(sys_siL(ikey,15)))
      call arrlst_printArrayList(rgraph%rEdges,itable)

      ! Proceed with right child if it exists
      if (rgraph%rVertices%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call inorderSparse(rgraph%rVertices%p_Kchild(TRIGHT,i))
    end subroutine inorderSparse
  end subroutine grph_printGraph

  ! ***************************************************************************

!<function>

  function grph_hasVertex(rgraph,iVertex,iVertexPosition) result(bexists)

!<description>
    ! This function returns TRUE if the graph has the given vertex.
    ! Otherwise, it returns FALSE. If the optional parameter
    ! iVertexPosition is present, then the function will also
    ! return the position of the vertex with number iVertex in the
    ! tree.
!</description>

!<input>
    ! Number of the vertex
    integer, intent(in) :: iVertex
!</input>

!<inputoutput>
    ! The graph
    type(t_graph), intent(inout) :: rgraph
!</inputoutput>

!<output>
    ! OPTIONAL: position of vertex in list
    integer, intent(out), optional :: iVertexPosition
!</output>

!<result>
    ! Flag for existence of the vertex
    logical :: bexists
!</result>
!</function>

    ! local variabes
    integer :: ipred

    ! Search in the vertex list
    bexists = (btree_searchInTree(rgraph%rVertices,iVertex,ipred).eq.BTREE_FOUND)

    ! Return vertex position if required
    if (bexists .and. present(iVertexPosition)) then
      iVertexPosition = rgraph%rVertices%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))
    end if
  end function grph_hasVertex

  ! ***************************************************************************

!<subroutine>

  subroutine grph_insertVertex(rgraph,iVertex,iVertexPosition,iEdgePosition)

!<description>
    ! This subroutine inserts the vertex with number iVertex into the graph.
    ! In addition, the trivial edge (iVertex,iVertex) is inserted.
    ! If the optional parameter iVertexPosition is present, then
    ! the posiion of vertex iVertex in the tree is returned.
    ! If the optional parameter iEdgePosition is present, then the position
    ! of the trivial edge (iVertex,iVertex) is returned
!</description>

!<input>
    ! Number of the vertex
    integer, intent(in) :: iVertex
!</input>

!<inputoutput>
    ! The graph
    type(t_graph), intent(inout) :: rgraph
!</inputoutput>

!<output>
    ! OPTIONAL: position of vertex in list
    integer, intent(out), optional :: iVertexPosition

    ! OPTIONAL: position of edge (iVertex,iVertex)
    integer, intent(out), optional :: iEdgePosition
!</output>
!</subroutine>

    ! local variabes
    integer :: iposVertex,iposEdge,itable
    
    ! Try to add entry with key iVertex
    if (rgraph%bisDense) then
      itable=iVertex
      call btree_insertIntoTree(rgraph%rVertices,iVertex,iposOpt=iposVertex)
    else
      itable=rgraph%NVT+1
      call btree_insertIntoTree(rgraph%rVertices,iVertex,IData=(/itable/),iposOpt=iposVertex)
    end if

    ! Return if vertex already exists
    if (iposVertex < 0) then
      if (present(iVertexPosition)) iVertexPosition=abs(iposVertex)
      return
    end if

    ! Increase number of vertices by one
    rgraph%NVT=rgraph%NVT+1
    
    ! Add trivial edge (iVertex,iVertex) to the list of edges
    select case(rgraph%cgraphFormat)

    case(GRPH_GRAPH7,&
         GRPH_GRAPH9,&
         GRPH_GRAPHORDERED_DIRECTED,&
         GRPH_GRAPHUNORDERED_DIRECTED)
      call arrlst_prependToArrayList(rgraph%rEdges,itable,iVertex,iposEdge)
      
    case DEFAULT
      call output_line('Unsupported graph format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'grph_insertVertex')
      call sys_halt()
    end select
    
    ! Increase number of edges by one
    rgraph%NEDGE=rgraph%NEDGE+1
    
    ! Return vertex/edge position?
    if (present(iVertexPosition)) iVertexPosition=iposVertex
    if (present(iEdgePosition))   iEdgePosition=iposEdge
  end subroutine grph_insertVertex

  ! ***************************************************************************

!<subroutine>

  subroutine grph_removeVertex(rgraph,iVertex,ireplacementVertex)

!<description>
    ! This subroutine removes the vertex with number iVertex from the graph
    ! and eliminates all of its incoming and outgoing edges.
    ! If the optional parameter ireplacemantVertex is true, then vertex iVertex
    ! is removed and the vertex with number ireplacementVertex is moved at its
    ! former position. This can be useful, if one needs a complete set of
    ! vertices, e.g., 1..NVT without "holes" at the positions of removed vertices.
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
    integer :: itable,jtable,ireplaceTable,jVertex,ireplaceVertex,ipred,ipos,iposVertex
    logical :: bdoReplace

    ! Check if vertex is present in tree?
    if (btree_searchInTree(rgraph%rVertices,iVertex,ipred).eq.BTREE_NOT_FOUND) then
      call output_line('Vertex does not exist in graph!',&
          OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
      call sys_halt()
    end if

    ! Replace vertex or leave "holes"?
    bdoReplace=present(ireplacementVertex)
    if (bdoReplace) then
      ireplaceVertex=ireplacementVertex
    else
      ireplaceVertex=0
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
      
      case(GRPH_GRAPH7,&
           GRPH_GRAPH9)
        ! We are in the lucky position, that the graph is undirected, that is,
        ! there exits an edge (i,j) if and only if there exists the edge (j,i).
        ! Hence, we do not have to loop through the list of all vertices but
        ! it suffices to visit those which are present in the adjacency list
        ! of vertex iVertex.
        
        ! Step 1: Delete corresponding vertex from tree. Note that for dense
        !         graphs the vertex tree does not store further information and
        !         can be eliminated in the first step. If the vertex should be
        !         replaced by the last vertex, then the last vertex is
        !         removed from the tree instead of vertex numbered iVertex.
        if (btree_deleteFromTree(rgraph%rVertices,&
            merge(ireplaceVertex,iVertex,bdoReplace)).eq.BTREE_NOT_FOUND) then
          call output_line('Vertex does not exist in graph!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
          call sys_halt()
        end if
               
        ! Step 2: Loop through adjacency list of vertex iVertex and delete all
        !         edges (jVertex,iVertex) from the adjacency list of jVertex.

        ! Find position of first entry in adjacency list
        iposVertex=arrlst_getNextInArrayList(rgraph%rEdges,iVertex,.true.)
        do while(iposVertex .ne. ARRLST_NULL)
          
          ! Get number of adjacent vertex
          jVertex=rgraph%rEdges%p_IData(iposVertex)

          ! Get position of next entry in adjacency list
          iposVertex=arrlst_getNextInArrayList(rgraph%rEdges,iVertex,.false.)
          
          ! Do nothing if both vertices are the same. The adjacency list of
          ! iVertex is removed/replaced anyway so we can leave it 'as is'
          if (iVertex .eq. jVertex) cycle

          ! Remove vertex iVertex from adjacency list of vertex jVertex
          if (arrlst_deleteFromArrayList(rgraph%rEdges,jVertex,iVertex).eq.&
              ARRAYLIST_NOT_FOUND) then
            call output_line('Unable to delete vertex from adjacency list!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            call sys_halt()
          end if
           
          ! Decrease number of edges by two; for the edge (iVertex,jVertex)
          ! and for the edge (jVertex,iVertex) that exists in an undirected graph
          rgraph%NEDGE = rgraph%NEDGE-1
        end do
        
        ! Now, vertex iVertex does no longer exist in any adjacency list.
        ! Check if replacement vertex is different from iVertex.
        if (bdoReplace .and. (iVertex .ne. ireplaceVertex)) then
          
          ! Remove the trivial edge (ireplaceVertex,ireplaceVertex)
          if (arrlst_deleteFromArrayList(rgraph%rEdges,ireplaceVertex,ireplaceVertex).eq.&
              ARRAYLIST_NOT_FOUND) then
            call output_line('Unable to delete vertex from adjacency list!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            call sys_halt()
          end if
          
          ! Swap adjacency list of vertices iVertex and ireplaceVertex
          ! From now onward, the adjacencey list of the replacement vertex is
          ! stored at position iVertex and vice versa. Keep this in mind!
          call arrlst_swapArrayList(rgraph%rEdges,iVertex,ireplaceVertex)

          ! Release adjacency list of vertex ireplaceVertex
          call arrlst_releaseArrayList(rgraph%rEdges,ireplaceVertex)

          ! Look for position of trivial edge (iVertex,iVertex)
          if (arrlst_searchInArrayList(rgraph%rEdges,iVertex,iVertex,ipred).eq.&
              ARRAYLIST_FOUND) then
            call output_line('Vertex already exists in adjacency list!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            call sys_halt()
          end if

          ! Insert trivial edge (iVertex,iVertex) into adjacency list
          call arrlst_insertIntoArrayList(rgraph%rEdges,iVertex,iVertex,ipred,ipos)

          ! Step 3(a): Loop through adjacency list of vertex ireplaceVertex
          !            (which is already stored at its new position iVertex) and
          !            delete all edges (jVertex,ireplaceVertex) from the adjacency
          !            lists of jVertex. Afterwards, add the edge (jVertex,iVertex)
          !            to the adjacency list of jVertex.

          ! Find position of first entry in adjacency list
          iposVertex=arrlst_getNextInArrayList(rgraph%rEdges,iVertex,.true.)
          do while(iposVertex .ne. ARRLST_NULL)
            
            ! Get number of adjacent vertex
            jVertex=rgraph%rEdges%p_IData(iposVertex)

            ! Get position of next entry in adjacency list
            iposVertex=arrlst_getNextInArrayList(rgraph%rEdges,iVertex,.false.)
            
            ! Do nothing if jVertex is identical iVertex
            if (jVertex .eq. iVertex) cycle
            
            ! Remove vertex ireplaceVertex from adjacency list of vertex jVertex
            if (arrlst_deleteFromArrayList(rgraph%rEdges,jVertex,ireplaceVertex).eq.&
                ARRAYLIST_NOT_FOUND) then
              call output_line('Unable to delete vertex from adjacency list!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
              call sys_halt()
            end if

            ! Look for position of edge (jVertex,iVertex)
            if (arrlst_searchInArrayList(rgraph%rEdges,jVertex,iVertex,ipred).eq.&
                ARRAYLIST_FOUND) then
              call output_line('Vertex already exists in adjacency list!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
              call sys_halt()
            end if
            
            ! Insert edge (jVertex,iVertex) into adjacency list
            call arrlst_insertIntoArrayList(rgraph%rEdges,jVertex,iVertex,ipred,ipos)
          end do
                  
          ! Decrease number of vertices by one
          rgraph%NVT = rgraph%NVT-1
          
          ! Decrease number of edges by one; for the trivial edge (lastVertex,lastVertex)
          rgraph%NEDGE = rgraph%NEDGE-1

        else
          
          ! Step 3(b): Release adjacency list of vertex iVertex
          call arrlst_releaseArrayList(rgraph%rEdges,iVertex)
          
          ! Decrease number of vertices by one
          rgraph%NVT = rgraph%NVT-1
          
          ! Decrease number of edges by one; for the trivial edge (iVertex,iVertex)
          rgraph%NEDGE = rgraph%NEDGE-1
        end if

      case DEFAULT
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

      case(GRPH_GRAPH7,&
           GRPH_GRAPH9)
        ! We are in the lucky position, that the graph is undirected, that is,
        ! there exits an edge (i,j) if and only if there exists the edge (j,i).
        ! Hence, we do not have to loop through the list of all vertices but
        ! it suffices to visit those which are present in the adjacency list
        ! of vertex iVertex.
        
        ! Get table for iVertex
        if (btree_searchInTree(rgraph%rVertices,iVertex,ipred).eq.BTREE_NOT_FOUND) then
          call output_line('Vertex does not exist in graph!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
          call sys_halt()
        end if
        ipos   = rgraph%rVertices%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))
        itable = rgraph%rVertices%p_IData(1,ipos)
        
        ! Get table for ireplaceVertex if required
        if (bdoReplace) then
          if (btree_searchInTree(rgraph%rVertices,ireplaceVertex,ipred).eq.BTREE_NOT_FOUND) then
            call output_line('Vertex does not exist in graph!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            call sys_halt()
          end if
          ipos          = rgraph%rVertices%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))
          ireplaceTable = rgraph%rVertices%p_IData(1,ipos)
        else
          ireplaceTable=0
        end if
        
        ! Step 1: Delete corresponding vertex from tree. Note that for dense
        !         graphs the vertex tree does not store further information and
        !         can be eliminated in the first step. If the vertex should be
        !         replaced by the last vertex, then the last vertex is
        !         removed from the tree instead of vertex numbered iVertex.
        
        if (btree_deleteFromTree(rgraph%rVertices,&
            merge(ireplaceTable,iVertex,bdoReplace)).eq.BTREE_NOT_FOUND) then
          call output_line('Vertex does not exist in graph!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
          call sys_halt()
        end if
        
        ! Step 2: Loop through adjacency list of vertex iVertex and delete all
        !         edges (jVertex,iVertex) from the adjacency lists of jVertex.
        
        ! Find position of first entry in adjacency list
        iposVertex=arrlst_getNextInArrayList(rgraph%rEdges,itable,.true.)
        do while(iposVertex .ne. ARRLST_NULL)
          
          ! Get number of adjacent vertex
          jVertex=rgraph%rEdges%p_IData(iposVertex)
          
          ! Get position of next entry in adjacency list
          iposVertex=arrlst_getNextInArrayList(rgraph%rEdges,itable,.false.)
          
          ! Do nothing if both vertices are the same
          if (iVertex .eq. jVertex) cycle

          ! In addition, do nothing if the current vertex is identical to
          ! the replacement vertex. Otherwise, we would have to re-insert
          ! it afterwards. Hence, it does not make sense to remove before.
          if (bdoReplace .and. (ireplaceVertex .eq. jVertex)) cycle

          ! Get table for jVertex
          if (btree_searchInTree(rgraph%rVertices,jVertex,ipred).eq.BTREE_NOT_FOUND) then
            call output_line('Vertex does not exist in graph!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            call sys_halt()
          end if
          ipos   = rgraph%rVertices%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))
          jtable = rgraph%rVertices%p_IData(1,ipos)

          ! Remove vertex iVertex from adjacency list of vertex jVertex
          if (arrlst_deleteFromArrayList(rgraph%rEdges,jtable,iVertex).eq.&
              ARRAYLIST_NOT_FOUND) then
            call output_line('Unable to delete vertex from adjacency list!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            call sys_halt()
          end if
          
          ! Decrease number of edges by two; for the edge (iVertex,jVertex)
          ! and for the edge (jVertex,iVertex) that exists in an undirected graph
          rgraph%NEDGE = rgraph%NEDGE-2
        end do

        ! Now, vertex iVertex does no longer exist in any adjacency list.
        ! Check if replacement vertex needs to be moved to position iVertex.
        if (bdoReplace .and. (iVertex .ne. ireplaceVertex)) then
          
          ! Step 3(a): Loop through adjacency list of vertex ireplaceVertex and
          !            delete all edges (jVertex,ireplaceVertex) from the adjacency
          !            lists of jVertex. Afterwards, add the edge (jVertex,iVertex)
          !            to the adjacency list of jVertex. Finally, swap adjacency
          !            list of vertices iVertex and ireplaceVertex.
          
          ! Find position of first entry in adjacency list
          iposVertex=arrlst_getNextInArrayList(rgraph%rEdges,ireplaceTable,.true.)
          do while(iposVertex .ne. ARRLST_NULL)

            ! Get number of adjacent vertex
            jVertex=rgraph%rEdges%p_IData(iposVertex)
            
            ! Get position of next entry in adjacency list
            iposVertex=arrlst_getNextInArrayList(rgraph%rEdges,ireplaceTable,.false.)
            
            ! Do nothing if both vertices are the same. This situation required
            ! special treatment (see below)
            if ((ireplaceVertex .eq. jVertex) .or. iVertex .eq. jVertex) cycle
            
            ! Get table for jVertex
            if (btree_searchInTree(rgraph%rVertices,jVertex,ipred).eq.BTREE_NOT_FOUND) then
              call output_line('Vertex does not exist in graph!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
              call sys_halt()
            end if
            ipos   = rgraph%rVertices%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))
            jtable = rgraph%rVertices%p_IData(1,ipos)
            
            ! Remove vertex lastVertex from adjacency list of vertex jVertex
            if (arrlst_deleteFromArrayList(rgraph%rEdges,jtable,ireplaceVertex).eq.&
                ARRAYLIST_NOT_FOUND) then
              call output_line('Unable to update vertex in adjacency list!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
              call sys_halt()
            end if
            
            ! Look for position of edge (jVertex,iVertex)
            if (arrlst_searchInArrayList(rgraph%rEdges,jtable,&
                iVertex,ipred).eq.ARRAYLIST_FOUND) then
              call output_line('Vertex replacement already exists in adjacency list!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
              call sys_halt()
            end if
            
            ! Insert edge (jVertex,iVertex) into adjacency list
            call arrlst_insertIntoArrayList(rgraph%rEdges,jtable,iVertex,ipred,ipos)
          end do

          ! Remove the trivial edge (ireplaceVertex,ireplaceVertex) from the adjacency
          ! list. Note that the trivial edge (iVertex,iVertex) still exists and has
          ! not been removed in the upper removal loop
          if (arrlst_deleteFromArrayList(rgraph%rEdges,ireplaceTable,ireplaceVertex).eq.&
              ARRAYLIST_NOT_FOUND) then
            call output_line('Unable to update vertex in adjacency list!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            call sys_halt()
          end if

          ! Swap adjacency list of vertices iVertex and ireplaceVertex
          call arrlst_swapArrayList(rgraph%rEdges,itable,ireplaceTable)

          ! Release adjacency list of vertex ireplaceVertex
          call arrlst_releaseArrayList(rgraph%rEdges,ireplaceTable)
          
          ! Decrease number of vertices by one
          rgraph%NVT = rgraph%NVT-1
          
          ! Decrease number of edges by one; for the trivial edge (lastVertex,lastVertex)
          rgraph%NEDGE = rgraph%NEDGE-1

        else
          
          ! Step 3(b): Release adjacency list of vertex iVertex
          call arrlst_releaseArrayList(rgraph%rEdges,itable)
          
          ! Decrease number of vertices by one
          rgraph%NVT = rgraph%NVT-1
          
          ! Decrease number of edges by one; for the trivial edge (iVertex,iVertex)
          rgraph%NEDGE = rgraph%NEDGE-1
        end if
        
      case DEFAULT
        call output_line('Invalid graph format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
        call sys_halt()
      end select
    end if
  end subroutine grph_removeVertex

  ! ***************************************************************************

!<function>

  function grph_hasEdge(rgraph,iFromVertex,iToVertex,iEdgePosition) result(bexists)

!<description>
    ! This function returns TRUE if the graph has the given edge from
    ! vertex iFromVertex to vertex iToVertex.
    ! Otherwise, it returns FALSE. If the optional parameter iEdgePosition
    ! is present, then the function will also return the position of the
    ! edge (iFromVertex,iToVertex) in the array list.
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

!<output>
    ! OPTIONAL: position of vertex in list
    integer, intent(out), optional :: iEdgePosition
!</output>

!<result>
    ! Flag for existence of the vertex
    logical :: bexists
!</result>
!</function>

    ! local variables
    integer :: itable,ipred,ipos

    ! Determine table
    if (rgraph%bisDense) then
      itable = iFromVertex
    else
      if (btree_searchInTree(rgraph%rVertices,iFromVertex,ipred).eq.BTREE_NOT_FOUND) then
        bexists = .false.
        return
      end if
      ipos = rgraph%rVertices%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))
      itable = rgraph%rVertices%p_IData(1,ipos)
    end if

    bexists = (arrlst_searchInArrayList(rgraph%rEdges,itable,iToVertex,ipos).eq.ARRAYLIST_FOUND)

    ! Return edge position if required
    if (bexists .and. present(iEdgePosition)) then
      iEdgePosition = ipos
    end if
  end function grph_hasEdge

  ! ***************************************************************************

!<subroutine>

  subroutine grph_insertEdge(rgraph,iFromVertex,iToVertex,iToEdgePosition,iFromEdgePosition)

!<description>
    ! This subroutine inserts the edge between vertices iFromVertex and
    ! iToVertex. If the graph is undirected, then the opposite edge
    ! (iToVertex,iFromVertex) is also added to attain symmetry of the
    ! sparsity pattern.
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

!<output>
    ! OPTIONAL: Position of edge (iToVertex,iFromVertex) in the array list
    integer, intent(out), optional :: iToEdgePosition

    ! OPTIONAL: Position of edge (iFromVertex,iToVertex) in the array list
    integer, intent(out), optional :: iFromEdgePosition
!</output>
!</subroutine>

    ! local variables
    integer :: itable,ipred,ipos
    
    ! What graph format are we?
    select case(rgraph%cgraphFormat)

    case(GRPH_GRAPHUNORDERED_DIRECTED)

      if (rgraph%bisDense) then
        
        ! Insert entry iToVertex into adjacency list of vertex iFromVertex
        if (arrlst_searchInArrayList(rgraph%rEdges,iToVertex,&
            iFromVertex,ipred).eq.ARRAYLIST_NOT_FOUND) then
          
          call arrlst_appendToArrayList(rgraph%rEdges,iToVertex,iFromVertex,ipos)
          rgraph%NEDGE = rgraph%NEDGE+1

          if (present(iToEdgePosition)) iToEdgePosition=ipos
        end if
      else

        ! Get table associated with vertex iToVertex
        if (btree_searchInTree(rgraph%rVertices,iToVertex,ipred).eq.BTREE_NOT_FOUND) then
          call output_line('Vertex does not exists!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_insertEdge')
          call sys_halt()
        end if
        ipos = rgraph%rVertices%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))
        itable = rgraph%rVertices%p_IData(1,ipos)

        ! Insert entry iToVertex into adjacency list of vertex iFromVertex
        if (arrlst_searchInArrayList(rgraph%rEdges,itable,&
            iFromVertex,ipred).eq.ARRAYLIST_NOT_FOUND) then
          
          call arrlst_appendToArrayList(rgraph%rEdges,itable,iFromVertex,ipos)
          rgraph%NEDGE = rgraph%NEDGE+1

          if (present(iToEdgePosition)) iToEdgePosition=ipos
        end if
      end if


    case(GRPH_GRAPHORDERED_DIRECTED)

      if (rgraph%bisDense) then
        
        ! Insert entry iToVertex into adjacency list of vertex iFromVertex
        if (arrlst_searchInArrayList(rgraph%rEdges,iToVertex,&
            iFromVertex,ipred).eq.ARRAYLIST_NOT_FOUND) then
          
          call arrlst_insertIntoArrayList(rgraph%rEdges,iToVertex,iFromVertex,ipred,ipos)
          rgraph%NEDGE = rgraph%NEDGE+1

          if (present(iToEdgePosition)) iToEdgePosition=ipos
        end if
      else

        ! Get table associated with vertex iToVertex
        if (btree_searchInTree(rgraph%rVertices,iToVertex,ipred).eq.BTREE_NOT_FOUND) then
          call output_line('Vertex does not exists!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_insertEdge')
          call sys_halt()
        end if
        ipos = rgraph%rVertices%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))
        itable = rgraph%rVertices%p_IData(1,ipos)

        ! Insert entry iToVertex into adjacency list of vertex iFromVertex
        if (arrlst_searchInArrayList(rgraph%rEdges,itable,&
            iFromVertex,ipred).eq.ARRAYLIST_NOT_FOUND) then
          
          call arrlst_insertIntoArrayList(rgraph%rEdges,itable,iFromVertex,ipred,ipos)
          rgraph%NEDGE = rgraph%NEDGE+1

          if (present(iToEdgePosition)) iToEdgePosition=ipos
        end if
      end if

      
    case(GRPH_GRAPH7,&
         GRPH_GRAPH9)
      
      if (rgraph%bisDense) then

        ! Insert entry iFromVertex into adjacency list of vertex iToVertex
        if (arrlst_searchInArrayList(rgraph%rEdges,iToVertex,&
            iFromVertex,ipred).ne.ARRAYLIST_FOUND) then

          call arrlst_insertIntoArrayList(rgraph%rEdges,iToVertex,iFromVertex,ipred,ipos)
          rgraph%NEDGE = rgraph%NEDGE+1

          if (present(iToEdgePosition)) iToEdgePosition=ipos
        end if
        
        ! Insert entry iToVertex into adjacency list of vertex iFromVertex
        if (arrlst_searchInArrayList(rgraph%rEdges,iFromVertex,&
            iToVertex,ipred).eq.ARRAYLIST_NOT_FOUND) then
          
          call arrlst_insertIntoArrayList(rgraph%rEdges,iFromVertex,iToVertex,ipred,ipos)
          rgraph%NEDGE = rgraph%NEDGE+1
          
          if (present(iFromEdgePosition)) iFromEdgePosition=ipos
        end if
      else
        
        ! Get table associated with vertex iToVertex
        if (btree_searchInTree(rgraph%rVertices,iToVertex,ipred).eq.BTREE_NOT_FOUND) then
          call output_line('Vertex does not exists!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_insertEdge')
          call sys_halt()
        end if
        ipos = rgraph%rVertices%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))
        itable = rgraph%rVertices%p_IData(1,ipos)

        ! Insert entry iFromVertex into adjacency list of vertex iToVertex
        if (arrlst_searchInArrayList(rgraph%rEdges,itable,&
            iFromVertex,ipred).ne.ARRAYLIST_FOUND) then
          
          call arrlst_insertIntoArrayList(rgraph%rEdges,itable,iFromVertex,ipred,ipos)
          rgraph%NEDGE = rgraph%NEDGE+1

          if (present(iToEdgePosition)) iToEdgePosition=ipos
        end if

        ! Get table associated with vertex iFromVertex
        if (btree_searchInTree(rgraph%rVertices,iFromVertex,ipred).eq.BTREE_NOT_FOUND) then
          call output_line('Vertex does not exists!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_insertEdge')
          call sys_halt()
        end if
        ipos = rgraph%rVertices%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))
        itable = rgraph%rVertices%p_IData(1,ipos)

        ! Insert entry iToVertex into adjacency list of vertex iFromVertex
        if (arrlst_searchInArrayList(rgraph%rEdges,itable,&
            iToVertex,ipred).eq.ARRAYLIST_NOT_FOUND) then
          
          call arrlst_insertIntoArrayList(rgraph%rEdges,itable,iToVertex,ipred,ipos)
          rgraph%NEDGE = rgraph%NEDGE+1

          if (present(iFromEdgePosition)) iFromEdgePosition=ipos
        end if
      end if

      
    case DEFAULT
      call output_line('Unsupported graph format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'grph_insertEdge')
      call sys_halt()
    end select
  end subroutine grph_insertEdge

  ! ***************************************************************************

!<subroutine>

  subroutine grph_removeEdge(rgraph,iFromVertex,iToVertex)

!<description>
    ! This subroutine removes the edge between vertices iFromVertex and
    ! iToVertex. If the graph is undirected, then the opposite edge
    ! (iToVertex,iFromVertex) is also removed to attain symmetry of the
    ! sparsity pattern.
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
    integer :: itable,ipred,ipos
    
    ! What graph format are we=
    select case(rgraph%cgraphFormat)

    case(GRPH_GRAPHUNORDERED_DIRECTED,&
         GRPH_GRAPHORDERED_DIRECTED)

      if (rgraph%bisDense) then

        ! Remove vertex iToVertex from table iFromVertex
        if (arrlst_deleteFromArrayList(rgraph%rEdges,iFromVertex,iToVertex)&
            .eq.ARRAYLIST_FOUND) rgraph%NEDGE = rgraph%NEDGE-1
      else

        ! Get table associated with vertex iFromVertex
        if (btree_searchInTree(rgraph%rVertices,iFromVertex,ipred).eq.BTREE_NOT_FOUND) then
          call output_line('Vertex does not exist!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_removeEdge')
          call sys_halt()
        end if
        ipos = rgraph%rVertices%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))
        itable = rgraph%rVertices%p_IData(1,ipos)
        
        ! Remove vertex iToVertex from table iFromVertex
        if (arrlst_deleteFromArrayList(rgraph%rEdges,itable,iToVertex)&
            .eq.ARRAYLIST_FOUND) rgraph%NEDGE = rgraph%NEDGE-1
      end if


    case(GRPH_GRAPH7,&
         GRPH_GRAPH9)

      if (rgraph%bisDense) then
        
        ! Remove vertex iToVertex from table iFromVertex
        if (arrlst_deleteFromArrayList(rgraph%rEdges,iFromVertex,iToVertex)&
            .eq.ARRAYLIST_FOUND) rgraph%NEDGE = rgraph%NEDGE-1

        ! Remove vertex iFromVertex from table iToVertex
        if (arrlst_deleteFromArrayList(rgraph%rEdges,iToVertex,iFromVertex)&
            .eq.ARRAYLIST_FOUND) rgraph%NEDGE = rgraph%NEDGE-1
      else

        ! Get table associated with vertex iFromVertex
        if (btree_searchInTree(rgraph%rVertices,iFromVertex,ipred).eq.BTREE_NOT_FOUND) then
          call output_line('Vertex does not exist!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_removeEdge')
          call sys_halt()
        end if
        ipos = rgraph%rVertices%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))
        itable = rgraph%rVertices%p_IData(1,ipos)
        
        ! Remove vertex iToVertex from table iFromVertex
        if (arrlst_deleteFromArrayList(rgraph%rEdges,itable,iToVertex)&
            .eq.ARRAYLIST_FOUND) rgraph%NEDGE = rgraph%NEDGE-1

        ! Get table associated with vertex iFromVertex
        if (btree_searchInTree(rgraph%rVertices,iToVertex,ipred).eq.BTREE_NOT_FOUND) then
          call output_line('Vertex does not exist!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_removeEdge')
          call sys_halt()
        end if
        ipos = rgraph%rVertices%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))
        itable = rgraph%rVertices%p_IData(1,ipos)

        ! Remove vertex iFromVertex from table iToVertex
        if (arrlst_deleteFromArrayList(rgraph%rEdges,itable,iFromVertex)&
            .eq.ARRAYLIST_FOUND) rgraph%NEDGE = rgraph%NEDGE-1
      end if


    case DEFAULT
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
    ! This subroutine makes a copy of a graph in memory.
    ! It does not make sense to share some information between graphs,
    ! so each vectors is physically copied from the source graph
    ! to the destination graph.
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

    call btree_duplicateTree(rgraph%rVertices,rgraphBackup%rVertices)
    call arrlst_duplicateArrayList(rgraph%rEdges,rgraphBackup%rEdges)
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
