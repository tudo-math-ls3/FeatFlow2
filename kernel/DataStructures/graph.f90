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

MODULE graph

  USE arraylist
  USE binarytree
  USE fsystem
  USE genoutput
  USE linearsystemscalar
  USE storage

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_graph
  PUBLIC :: grph_createGraph
  PUBLIC :: grph_createGraphFromMatrix
  PUBLIC :: grph_releaseGraph
  PUBLIC :: grph_generateMatrix
  PUBLIC :: grph_printGraph
  PUBLIC :: grph_hasVertex
  PUBLIC :: grph_insertVertex
  PUBLIC :: grph_removeVertex
  PUBLIC :: grph_hasEdge
  PUBLIC :: grph_insertEdge
  PUBLIC :: grph_removeEdge
  PUBLIC :: grph_infoGraph
  PUBLIC :: grph_duplicateGraph
  PUBLIC :: grph_restoreGraph

!<constants>
!<constantblock description="Global format flags for graphs">

  ! Unidentified graph format
  INTEGER, PARAMETER, PUBLIC :: GRPH_GRAPHUNDEFINED = 0

  ! Identifier for directed graph with unordered entries
  INTEGER, PARAMETER, PUBLIC :: GRPH_GRAPHUNORDERED_DIRECTED = 1

  ! Identifier for directed graph with ordered entries
  INTEGER, PARAMETER, PUBLIC :: GRPH_GRAPHORDERED_DIRECTED = 2

  ! Identifier for undirected graph with ordered entries in CSR format 9
  INTEGER, PARAMETER, PUBLIC :: GRPH_GRAPH9 = 9

  ! Identifier for undirected graph with ordered entries in CSR format 7,
  ! that is, CSR with diagonal element in front
  INTEGER, PARAMETER, PUBLIC :: GRPH_GRAPH7 = 7

!</constantblock>
!</constants>

!<types>
!<typeblock>

  ! This type contains all data structures to handle graphs.
  TYPE t_graph

    ! Format-tag. Identifies the format of the graph
    ! Can take one of the GRPH_GRAPHx format flags.
    INTEGER :: cgraphFormat = GRPH_GRAPHUNDEFINED

    ! Number of vertices
    INTEGER :: NVT = 0

    ! Number of edges
    INTEGER :: NEDGE = 0

    ! Is set to true, if the graph is dense, that is, a graph with NVT vertices
    ! contains the mainly the vertices labeled 1,...,NVT. In this case the vertex
    ! number coincides with the table under which its adjacency list is stored.
    ! If the graph is note marked as dense, then for each vertex IVT the table
    ! number is stored in the vertex list and must be looked-up each time, the
    ! adjacency information is required.
    ! Note that dense graphs can only store positive vertex labels !!!
    LOGICAL :: bisDense = .TRUE.

    ! Search tree for the list of vertices
    TYPE(t_btree) :: rVertices

    ! Array of lists for the edges
    TYPE(t_arraylist) :: rEdges

  END TYPE t_graph
    
!</typeblock>
!</types>

CONTAINS

!<subroutine>

  SUBROUTINE grph_createGraph(rgraph,cgraphFormat,nvtMax,nedgeMax,bisDense)

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
    ! 1.. (1^k)*NVT, k>> 1, then an enormous amount of memory would be wasted
    ! for the look-up tables. Then it makes sense to define the graph as
    ! "not dense". In this case, the tables are used sequentially if new
    ! vertices are inserted but finding the correct table number for vertex
    ! ivt requires a search in the binary tree, i.e. time O(log NVT).
!</description>

!<input>
    ! Format of the graph
    INTEGER, INTENT(IN) :: cgraphFormat

    ! Maximum number of vertices
    INTEGER(PREC_VECIDX), INTENT(IN) :: nvtMax

    ! Maximum number of edges
    INTEGER(PREC_MATIDX), INTENT(IN) :: nedgeMax

    ! OPTIONAL: Indicates if the graph is dense or not
    LOGICAL, INTENT(IN), OPTIONAL :: bisDense
!</input>

!<output>
    ! The sparsity graph
    TYPE(t_graph), INTENT(OUT) :: rgraph
!</output>
!</subroutine>

    ! Set dimensions
    rgraph%NVT    = 0
    rgraph%NEDGE  = 0
    
    IF (PRESENT(bisDense)) rgraph%bisDense=bisDense

    ! What graph format are we?
    SELECT CASE(cgraphFormat)

    CASE(GRPH_GRAPHORDERED_DIRECTED)
      rgraph%cgraphFormat = GRPH_GRAPHORDERED_DIRECTED

      ! Create search tree for the list of vertices
      IF (rgraph%bisDense) THEN
        CALL btree_createTree(rgraph%rVertices,rgraph%NVT+1,ST_INT,0,0,0)
      ELSE
        CALL btree_createTree(rgraph%rVertices,rgraph%NVT+1,ST_INT,0,0,1)
      END IF
      
      ! Create array of lists for edges
      CALL arrlst_createArrayList(rgraph%rEdges,nvtMax,nedgeMax,ST_INT,ARRAYLIST_INCREASING)

    CASE(GRPH_GRAPHUNORDERED_DIRECTED)
      rgraph%cgraphFormat = GRPH_GRAPHUNORDERED_DIRECTED

      ! Create search tree for the list of vertices
      IF (rgraph%bisDense) THEN
        CALL btree_createTree(rgraph%rVertices,rgraph%NVT+1,ST_INT,0,0,0)
      ELSE
        CALL btree_createTree(rgraph%rVertices,rgraph%NVT+1,ST_INT,0,0,1)
      END IF
      
      ! Create array of lists for edges
      CALL arrlst_createArrayList(rgraph%rEdges,nvtMax,nedgeMax,ST_INT,ARRAYLIST_UNORDERED)

    CASE(GRPH_GRAPH7)
      rgraph%cgraphFormat = GRPH_GRAPH7

      ! Create search tree for the list of vertices
      IF (rgraph%bisDense) THEN
        CALL btree_createTree(rgraph%rVertices,rgraph%NVT+1,ST_INT,0,0,0)
      ELSE
        CALL btree_createTree(rgraph%rVertices,rgraph%NVT+1,ST_INT,0,0,1)
      END IF

      ! Create array of lists for edges
      CALL arrlst_createArrayList(rgraph%rEdges,nvtMax,nedgeMax,ST_INT,ARRAYLIST_CSR7)

    CASE(GRPH_GRAPH9)
      rgraph%cgraphFormat = GRPH_GRAPH9

      ! Create search tree for the list of vertices
      IF (rgraph%bisDense) THEN
        CALL btree_createTree(rgraph%rVertices,rgraph%NVT+1,ST_INT,0,0,0)
      ELSE
        CALL btree_createTree(rgraph%rVertices,rgraph%NVT+1,ST_INT,0,0,1)
      END IF

      ! Create array of lists for edges
      CALL arrlst_createArrayList(rgraph%rEdges,nvtMax,nedgeMax,ST_INT,ARRAYLIST_INCREASING)
      
    CASE DEFAULT
      CALL output_line('Invalid matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'grph_createGraph')
      CALL sys_halt()
    END SELECT
  END SUBROUTINE grph_createGraph

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE grph_createGraphFromMatrix(rscalarMatrix,rgraph,nvtMax,nedgeMax)

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
    TYPE(t_matrixScalar), INTENT(IN) :: rscalarMatrix

    ! OPTIONAL: maximum number of vertices
    INTEGER(PREC_VECIDX), INTENT(IN), OPTIONAL :: nvtMax

    ! OPTIONAL: maximum number of edges
    INTEGER(PREC_MATIDX), INTENT(IN), OPTIONAL :: nedgeMax
!</input>

!<output>
    ! The sparsity graph
    TYPE(t_graph), INTENT(OUT) :: rgraph
!</output>

!</subroutine>

    ! local variables
    INTEGER :: h_Key
    INTEGER(PREC_VECIDX) :: nvt
    INTEGER(PREC_MATIDX) :: nedge
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol,p_Key

    ! Set dimensions
    rgraph%NVT   = rscalarMatrix%NEQ
    rgraph%NEDGE = rscalarMatrix%NA

    ! Estimate maximum number of vertices
    IF (PRESENT(nvtMax)) THEN
      nvt = MAX(nvtMax,rgraph%NVT)
    ELSE
      nvt = INT(1.5_DP*rgraph%NVT)
    END IF

    ! Estimate maximum number of edges
    IF (PRESENT(nedgeMax)) THEN
      nedge = MAX(nedgeMax,rgraph%NEDGE)
    ELSE
      nedge = INT(0.5_DP*(rgraph%NEDGE+rgraph%NVT))
    END IF

    ! What matrix format are we?
    SELECT CASE(rscalarMatrix%cmatrixFormat)

    CASE(LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)
      rgraph%cgraphFormat = GRPH_GRAPH7

      ! Create search tree for the list of vertices
      CALL btree_createTree(rgraph%rVertices,rgraph%NVT+1,ST_INT,0,0,0)
      
      ! Create array of lists for edges
      CALL arrlst_createArrayList(rgraph%rEdges,nvt,nedge,ST_INT,ARRAYLIST_CSR7)

      ! Set pointers
      CALL lsyssc_getbase_Kld(rscalarMatrix,p_Kld)
      CALL lsyssc_getbase_Kcol(rscalarMatrix,p_Kcol)

      ! Fill list of edges
      CALL arrlst_copyArrayListTable(p_Kcol,rgraph%rEdges,p_Kld)

      ! Generate p_Key = array [1,2,3,...,NVT]
      CALL storage_new('grph_createGraphFromMatrix','p_Key',rgraph%NVT,ST_INT,&
          h_Key,ST_NEWBLOCK_ORDERED)
      CALL storage_getbase_int(h_Key,p_Key)
      CALL btree_copyToTree(p_Key,rgraph%rVertices)
      CALL storage_free(h_Key)
      

    CASE(LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)
      rgraph%cgraphFormat = GRPH_GRAPH9
      
      ! Create search tree for the list of vertices
      CALL btree_createTree(rgraph%rVertices,rgraph%NVT+1,ST_INT,0,0,0)
      
      ! Create array of lists for edges
      CALL arrlst_createArrayList(rgraph%rEdges,nvt,nedge,ST_INT,ARRAYLIST_INCREASING)

      ! Set pointers
      CALL lsyssc_getbase_Kld(rscalarMatrix,p_Kld)
      CALL lsyssc_getbase_Kcol(rscalarMatrix,p_Kcol)

      ! Fill list of edges
      CALL arrlst_copyArrayListTable(p_Kcol,rgraph%rEdges,p_Kld)

      ! Generate p_Key = array [1,2,3,...,NVT]
      CALL storage_new('grph_createGraphFromMatrix','p_Key',rgraph%NVT,ST_INT,&
          h_Key,ST_NEWBLOCK_ORDERED)
      CALL storage_getbase_int(h_Key,p_Key)
      CALL btree_copyToTree(p_Key,rgraph%rVertices)
      CALL storage_free(h_Key)


    CASE DEFAULT
      CALL output_line('Invalid matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'grph_createGraphFromMatrix')
      CALL sys_halt()
    END SELECT
  END SUBROUTINE grph_createGraphFromMatrix
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE grph_releaseGraph(rgraph)

!<description>
    ! This routine releases an existing graph rgraph.
!</description>

!<inputoutput>
    ! The graph that should be released
    TYPE(t_graph), INTENT(INOUT) :: rgraph
!</inputoutput>
!</subroutine>

    ! Release search tree for the list of vertices
    CALL btree_releaseTree(rgraph%rVertices)

    ! Release array of lists for edges
    CALL arrlst_releaseArrayList(rgraph%rEdges)

    ! Reset data
    rgraph%cgraphFormat = GRPH_GRAPHUNDEFINED
    rgraph%NVT          = 0
    rgraph%NEDGE        = 0
    rgraph%bisDense     = .TRUE.
  END SUBROUTINE grph_releaseGraph

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE grph_generateMatrix(rgraph,rscalarMatrix)

!<description>
    ! This subroutine generates a sparse matrix rscalarMatrix stored in
    ! format CSR 7 or 9 whereby the sparsity pattern is adopted from
    ! the graph rgraph.
!</description>

!<input>
    ! The graph
    TYPE(t_graph), INTENT(INOUT) :: rgraph
!</input>

!<inputoutput>
    ! The matrix
    TYPE(t_matrixScalar), INTENT(INOUT) :: rscalarMatrix
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kdiagonal
    INTEGER(PREC_TABLEIDX) :: itable
    INTEGER(PREC_VECIDX) :: ieq
    INTEGER(PREC_MATIDX) :: ia,ncols
    INTEGER(PREC_TREEIDX) :: ipred,ipos
    INTEGER(I32) :: isize
    
    ! Check that matrix and graph have the same format
    SELECT CASE(rscalarMatrix%cmatrixFormat)
      
    CASE(LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)
      IF (rgraph%cgraphFormat .NE. GRPH_GRAPH7) THEN
        CALL output_line('Matrix/graph have incompatible format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'grph_generateMatrix')
        CALL sys_halt()
      END IF

    CASE(LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)
      IF (rgraph%cgraphFormat .NE. GRPH_GRAPH9) THEN
        CALL output_line('Matrix/graph have incompatible format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'grph_generateMatrix')
        CALL sys_halt()
      END IF

    CASE DEFAULT
      CALL output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'grph_generateMatrix')
      CALL sys_halt()
    END SELECT
 
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
    rscalarMatrix%NEQ   = MAX(rgraph%NVT,rgraph%rEdges%NTABLE)
    rscalarMatrix%NCOLS = rscalarMatrix%NEQ

    IF (rgraph%bisDense) THEN
      
      ! Convert array list to matrix
      CALL arrlst_copyArrayListTable(rgraph%rEdges,rscalarMatrix%h_Kcol,rscalarMatrix%h_Kld)
      
    ELSE

      ! Check if matrix is empty
      IF ((rscalarMatrix%NEQ .EQ. 0) .OR. rscalarMatrix%NA .EQ. 0) RETURN

      ! Convert array list step-by-step
      IF (rscalarMatrix%h_Kld .EQ. ST_NOHANDLE) THEN
        CALL storage_new('grph_generateMatrix','p_Kld',&
            rscalarMatrix%NEQ+1,ST_INT,rscalarMatrix%h_Kld,ST_NEWBLOCK_NOINIT)
      ELSE
        CALL storage_getsize(rscalarMatrix%h_Kld,isize)
        IF (isize < rscalarMatrix%NEQ+1) THEN
          CALL storage_realloc('grph_generateMatrix',&
              rscalarMatrix%NEQ+1,rscalarMatrix%h_Kld,ST_NEWBLOCK_NOINIT,.FALSE.)
        END IF
      END IF

      IF (rscalarMatrix%h_Kcol .EQ. ST_NOHANDLE) THEN
        CALL storage_new('grph_generateMatrix','p_Kcol',&
            rscalarMatrix%NA,ST_INT,rscalarMatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
      ELSE
        CALL storage_getsize(rscalarMatrix%h_Kcol,isize)
        IF (isize < rscalarMatrix%NA) THEN
          CALL storage_realloc('grph_generateMatrix',&
              rscalarMatrix%NA,rscalarMatrix%h_Kcol,ST_NEWBLOCK_NOINIT,.FALSE.)
        END IF
      END IF

      ! Set pointers
      CALL lsyssc_getbase_Kld(rscalarMatrix,p_Kld)
      CALL lsyssc_getbase_Kcol(rscalarMatrix,p_Kcol)

      ! Initialization
      ia=1

      ! Loop over all equations
      DO ieq=1,rscalarMatrix%NEQ

        p_Kld(ieq) = ia

        ! Check if vertex with number ieq exists
        IF (btree_searchInTree(rgraph%rVertices,ieq,ipred).EQ.BTREE_FOUND) THEN
          ipos = rgraph%rVertices%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
          itable = rgraph%rVertices%p_IData(1,ipos)

          ! Restore row from table
          CALL arrlst_copyArrayList(rgraph%rEdges,itable,p_Kcol(ia:),ncols)
        END IF

        ia = ia+ncols
      END DO
      p_Kld(rscalarMatrix%NEQ+1)=ia
    END IF

    ! Do we have to rebuild the diagonal?
    IF (rscalarMatrix%cmatrixFormat .EQ. LSYSSC_MATRIX9 .OR.&
        rscalarMatrix%cmatrixFormat .EQ. LSYSSC_MATRIX9INTL) THEN

      ! Create new memory or resize existing memory      
      IF (rscalarMatrix%h_Kdiagonal .EQ. ST_NOHANDLE) THEN
        CALL storage_new('grph_generateMatrix','p_Kdiagonal',&
            rscalarMatrix%NEQ, ST_INT, rscalarMatrix%h_Kdiagonal,ST_NEWBLOCK_NOINIT)
      ELSE
        CALL storage_getsize(rscalarMatrix%h_Kdiagonal,isize)
        IF (isize < rscalarMatrix%NEQ) THEN
          CALL storage_realloc('grph_generateMatrix',&
              rscalarMatrix%NEQ, rscalarMatrix%h_Kdiagonal, ST_NEWBLOCK_NOINIT, .FALSE.)
        END IF
      END IF

      ! Set pointers
      CALL lsyssc_getbase_Kld(rscalarMatrix, p_Kld)
      CALL lsyssc_getbase_Kcol(rscalarMatrix, p_Kcol)
      CALL lsyssc_getbase_Kdiagonal(rscalarMatrix, p_Kdiagonal)

      ! Rebuild array
      CALL lsyssc_rebuildKdiagonal(p_Kcol, p_Kld, p_Kdiagonal, rscalarMatrix%NEQ)
    END IF
  END SUBROUTINE grph_generateMatrix

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE grph_printGraph(rgraph)

!<description>
    ! This subroutine prints the content of the graph.
!</description>

!<inputoutput>
    ! The graph
    TYPE(t_graph), INTENT(INOUT) :: rgraph
!</inputoutput>
!</subroutine>

    ! Are we dense graph?
    IF (rgraph%bisDense) THEN
      CALL inorderDense(rgraph%rVertices%p_Kchild(TRIGHT,TROOT))
    ELSE
      CALL inorderSparse(rgraph%rVertices%p_Kchild(TRIGHT,TROOT))
    END IF

  CONTAINS
    
    ! Here, the real working routines follow.
    
    !**************************************************************
    ! Inorder traversal of the tree storing the vertex numbers
    RECURSIVE SUBROUTINE inorderDense(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i

      ! local variables
      INTEGER(PREC_VECIDX)   :: ikey

      ! Check if position is valid
      IF (i .EQ. TNULL) RETURN

      ! Proceed with left child if it exists
      IF (rgraph%rVertices%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderDense(rgraph%rVertices%p_Kchild(TLEFT,i))

      ! We are in the lucky position that itable = ikey
      ikey   = rgraph%rVertices%p_IKey(i)
      
      CALL output_line('Vertex number: '//TRIM(sys_siL(ikey,15)))
      CALL arrlst_printArrayList(rgraph%rEdges,ikey)
        
      ! Proceed with right child if it exists
      IF (rgraph%rVertices%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderDense(rgraph%rVertices%p_Kchild(TRIGHT,i))
    END SUBROUTINE inorderDense

     !**************************************************************
    ! Inorder traversal of the tree storing the vertex numbers
    RECURSIVE SUBROUTINE inorderSparse(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i

      ! local variables
      INTEGER(PREC_TABLEIDX) :: itable
      INTEGER(PREC_VECIDX)   :: ikey

      ! Check if position is valid
      IF (i .EQ. TNULL) RETURN

      ! Proceed with left child if it exists
      IF (rgraph%rVertices%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderSparse(rgraph%rVertices%p_Kchild(TLEFT,i))

      ! In general ikey != itable
      ikey   = rgraph%rVertices%p_IKey(i)
      itable = rgraph%rVertices%p_IData(1,i)

      CALL output_line('Vertex number: '//TRIM(sys_siL(ikey,15)))
      CALL arrlst_printArrayList(rgraph%rEdges,itable)

      ! Proceed with right child if it exists
      IF (rgraph%rVertices%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderSparse(rgraph%rVertices%p_Kchild(TRIGHT,i))
    END SUBROUTINE inorderSparse
  END SUBROUTINE grph_printGraph

  ! ***************************************************************************

!<function>

  FUNCTION grph_hasVertex(rgraph,iVertex,iVertexPosition) RESULT(bexists)

!<description>
    ! This function returns TRUE if the graph has the given vertex.
    ! Otherwise, it returns FALSE. If the optional parameter
    ! iVertexPosition is present, then the function will also
    ! return the position of the vertex with number iVertex in the
    ! tree.
!</description>

!<input>
    ! Number of the vertex
    INTEGER(PREC_VECIDX), INTENT(IN) :: iVertex
!</input>

!<inputoutput>
    ! The graph
    TYPE(t_graph), INTENT(INOUT) :: rgraph
!</inputoutput>

!<output>
    ! OPTIONAL: position of vertex in list
    INTEGER(PREC_TREEIDX), INTENT(OUT), OPTIONAL :: iVertexPosition
!</output>

!<result>
    ! Flag for existence of the vertex
    LOGICAL :: bexists
!</result>
!</function>

    ! local variabes
    INTEGER(PREC_TREEIDX) :: ipred

    ! Search in the vertex list
    bexists = (btree_searchInTree(rgraph%rVertices,iVertex,ipred).EQ.BTREE_FOUND)

    ! Return vertex position if required
    IF (bexists .AND. PRESENT(iVertexPosition)) THEN
      iVertexPosition = rgraph%rVertices%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
    END IF
  END FUNCTION grph_hasVertex

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE grph_insertVertex(rgraph,iVertex,iVertexPosition,iEdgePosition)

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
    INTEGER(PREC_VECIDX), INTENT(IN) :: iVertex
!</input>

!<inputoutput>
    ! The graph
    TYPE(t_graph), INTENT(INOUT) :: rgraph
!</inputoutput>

!<output>
    ! OPTIONAL: position of vertex in list
    INTEGER(PREC_TREEIDX), INTENT(OUT), OPTIONAL :: iVertexPosition

    ! OPTIONAL: position of edge (iVertex,iVertex)
    INTEGER(PREC_ARRAYLISTIDX), INTENT(OUT), OPTIONAL :: iEdgePosition
!</output>
!</subroutine>

    ! local variabes
    INTEGER(PREC_TREEIDX) :: iposVertex,iposEdge
    INTEGER(PREC_TABLEIDX) :: itable
    
    ! Try to add entry with key iVertex
    IF (rgraph%bisDense) THEN
      itable=iVertex
      CALL btree_insertIntoTree(rgraph%rVertices,iVertex,iposOpt=iposVertex)
    ELSE
      itable=rgraph%NVT+1
      CALL btree_insertIntoTree(rgraph%rVertices,iVertex,IData=(/itable/),iposOpt=iposVertex)
    END IF

    ! Return if vertex already exists
    IF (iposVertex < 0) THEN
      IF (PRESENT(iVertexPosition)) iVertexPosition=ABS(iposVertex)
      RETURN
    END IF

    ! Increase number of vertices by one
    rgraph%NVT=rgraph%NVT+1
    
    ! Add trivial edge (iVertex,iVertex) to the list of edges
    SELECT CASE(rgraph%cgraphFormat)

    CASE(GRPH_GRAPH7,GRPH_GRAPH9,&
        GRPH_GRAPHORDERED_DIRECTED,&
        GRPH_GRAPHUNORDERED_DIRECTED)
      CALL arrlst_prependToArrayList(rgraph%rEdges,itable,iVertex,iposEdge)
      
    CASE DEFAULT
      CALL output_line('Unsupported graph format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'grph_insertVertex')
      CALL sys_halt()
    END SELECT
    
    ! Increase number of edges by one
    rgraph%NEDGE=rgraph%NEDGE+1
    
    ! Return vertex/edge position?
    IF (PRESENT(iVertexPosition)) iVertexPosition=iposVertex
    IF (PRESENT(iEdgePosition))   iEdgePosition=iposEdge
  END SUBROUTINE grph_insertVertex

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE grph_removeVertex(rgraph,iVertex,ireplacementVertex)

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
    INTEGER(PREC_VECIDX), INTENT(IN) :: iVertex

    ! OPTIONAL: Number of the replacement vertex
    INTEGER(PREC_VECIDX), INTENT(IN), OPTIONAL :: ireplacementVertex
!</input>

!<inputoutput>
    ! The graph
    TYPE(t_graph), INTENT(INOUT) :: rgraph
!</inputoutput>
!</subroutine>

    ! local variabes
    INTEGER(PREC_TABLEIDX) :: itable,jtable,ireplaceTable
    INTEGER(PREC_VECIDX)   :: jVertex,ireplaceVertex
    INTEGER(PREC_TREEIDX)  :: ipred,ipos,iposVertex
    LOGICAL :: bdoReplace

    ! Check if vertex is present in tree?
    IF (btree_searchInTree(rgraph%rVertices,iVertex,ipred).EQ.BTREE_NOT_FOUND) THEN
      CALL output_line('Vertex does not exist in graph!',&
          OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
      CALL sys_halt()
    END IF

    ! Replace vertex or leave "holes"?
    bdoReplace=PRESENT(ireplacementVertex)
    IF (bdoReplace) THEN
      ireplaceVertex=ireplacementVertex
    ELSE
      ireplaceVertex=0
    END IF

    ! Are we dense graph?
    IF (rgraph%bisDense) THEN
      
      ! What graph format are ?
      SELECT CASE(rgraph%cgraphFormat)
        
      CASE(GRPH_GRAPHUNORDERED_DIRECTED,GRPH_GRAPHORDERED_DIRECTED)
        CALL output_line('We cannot remove directed graphs at the moment!',&
            OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
        CALL sys_halt()
      
      CASE(GRPH_GRAPH7,GRPH_GRAPH9)
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
        IF (btree_deleteFromTree(rgraph%rVertices,&
            MERGE(ireplaceVertex,iVertex,bdoReplace)).EQ.BTREE_NOT_FOUND) THEN
          CALL output_line('Vertex does not exist in graph!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
          CALL sys_halt()
        END IF
               
        ! Step 2: Loop through adjacency list of vertex iVertex and delete all 
        !         edges (jVertex,iVertex) from the adjacency list of jVertex.

        ! Find position of first entry in adjacency list
        iposVertex=arrlst_getNextInArrayList(rgraph%rEdges,iVertex,.TRUE.)
        DO WHILE(iposVertex .NE. ARRLST_NULL)
          
          ! Get number of adjacent vertex
          jVertex=rgraph%rEdges%p_IData(iposVertex)

          ! Get position of next entry in adjacency list
          iposVertex=arrlst_getNextInArrayList(rgraph%rEdges,iVertex,.FALSE.)
          
          ! Do nothing if both vertices are the same. The adjacency list of
          ! iVertex is removed/replaced anyway so we can leave it 'as is'
          IF (iVertex .EQ. jVertex) CYCLE

          ! Remove vertex iVertex from adjacency list of vertex jVertex
          IF (arrlst_deleteFromArrayList(rgraph%rEdges,jVertex,iVertex).EQ.&
              ARRAYLIST_NOT_FOUND) THEN
            CALL output_line('Unable to delete vertex from adjacency list!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            CALL sys_halt()
          END IF
           
          ! Decrease number of edges by two; for the edge (iVertex,jVertex)
          ! and for the edge (jVertex,iVertex) that exists in an undirected graph
          rgraph%NEDGE = rgraph%NEDGE-1
        END DO
        
        ! Now, vertex iVertex does no longer exist in any adjacency list.
        ! Check if replacement vertex is different from iVertex.
        IF (bdoReplace .AND. (iVertex .NE. ireplaceVertex)) THEN
          
          ! Remove the trivial edge (ireplaceVertex,ireplaceVertex)
          IF (arrlst_deleteFromArrayList(rgraph%rEdges,ireplaceVertex,ireplaceVertex).EQ.&
              ARRAYLIST_NOT_FOUND) THEN
            CALL output_line('Unable to delete vertex from adjacency list!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            CALL sys_halt()
          END IF
          
          ! Swap adjacency list of vertices iVertex and ireplaceVertex
          ! From now onward, the adjacencey list of the replacement vertex is
          ! stored at position iVertex and vice versa. Keep this in mind!
          CALL arrlst_swapArrayList(rgraph%rEdges,iVertex,ireplaceVertex)

          ! Release adjacency list of vertex ireplaceVertex
          CALL arrlst_releaseArrayList(rgraph%rEdges,ireplaceVertex)

          ! Look for position of trivial edge (iVertex,iVertex)
          IF (arrlst_searchInArrayList(rgraph%rEdges,iVertex,iVertex,ipred).EQ.&
              ARRAYLIST_FOUND) THEN             
            CALL output_line('Vertex already exists in adjacency list!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            CALL sys_halt()
          END IF

          ! Insert trivial edge (iVertex,iVertex) into adjacency list
          CALL arrlst_insertIntoArrayList(rgraph%rEdges,iVertex,iVertex,ipred,ipos)

          ! Step 3(a): Loop through adjacency list of vertex ireplaceVertex 
          !            (which is already stored at its new position iVertex) and 
          !            delete all edges (jVertex,ireplaceVertex) from the adjacency
          !            lists of jVertex. Afterwards, add the edge (jVertex,iVertex)
          !            to the adjacency list of jVertex.

          ! Find position of first entry in adjacency list
          iposVertex=arrlst_getNextInArrayList(rgraph%rEdges,iVertex,.TRUE.)
          DO WHILE(iposVertex .NE. ARRLST_NULL)
            
            ! Get number of adjacent vertex
            jVertex=rgraph%rEdges%p_IData(iposVertex)

            ! Get position of next entry in adjacency list
            iposVertex=arrlst_getNextInArrayList(rgraph%rEdges,iVertex,.FALSE.)
            
            ! Do nothing if jVertex is identical iVertex
            IF (jVertex .EQ. iVertex) CYCLE
            
            ! Remove vertex ireplaceVertex from adjacency list of vertex jVertex
            IF (arrlst_deleteFromArrayList(rgraph%rEdges,jVertex,ireplaceVertex).EQ.&
                ARRAYLIST_NOT_FOUND) THEN
              CALL output_line('Unable to delete vertex from adjacency list!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
              CALL sys_halt()
            END IF

            ! Look for position of edge (jVertex,iVertex)
            IF (arrlst_searchInArrayList(rgraph%rEdges,jVertex,iVertex,ipred).EQ.&
                ARRAYLIST_FOUND) THEN             
              CALL output_line('Vertex already exists in adjacency list!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
              CALL sys_halt()
            END IF
            
            ! Insert edge (jVertex,iVertex) into adjacency list
            CALL arrlst_insertIntoArrayList(rgraph%rEdges,jVertex,iVertex,ipred,ipos)            
          END DO
                  
          ! Decrease number of vertices by one
          rgraph%NVT = rgraph%NVT-1
          
          ! Decrease number of edges by one; for the trivial edge (lastVertex,lastVertex)
          rgraph%NEDGE = rgraph%NEDGE-1

        ELSE
          
          ! Step 3(b): Release adjacency list of vertex iVertex
          CALL arrlst_releaseArrayList(rgraph%rEdges,iVertex)
          
          ! Decrease number of vertices by one
          rgraph%NVT = rgraph%NVT-1
          
          ! Decrease number of edges by one; for the trivial edge (iVertex,iVertex)
          rgraph%NEDGE = rgraph%NEDGE-1
        END IF

      CASE DEFAULT
        CALL output_line('Invalid graph format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
        CALL sys_halt()
      END SELECT

    ELSE   ! - - - - - - The following is for non-dense graphs - - - - - -

      ! What graph format are ?
      SELECT CASE(rgraph%cgraphFormat)
        
      CASE(GRPH_GRAPHUNORDERED_DIRECTED,GRPH_GRAPHORDERED_DIRECTED)
        CALL output_line('We cannot remove directed graphs at the moment!',&
            OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
        CALL sys_halt()

      CASE(GRPH_GRAPH7,GRPH_GRAPH9)
        ! We are in the lucky position, that the graph is undirected, that is,
        ! there exits an edge (i,j) if and only if there exists the edge (j,i).
        ! Hence, we do not have to loop through the list of all vertices but
        ! it suffices to visit those which are present in the adjacency list
        ! of vertex iVertex.
        
        ! Get table for iVertex
        IF (btree_searchInTree(rgraph%rVertices,iVertex,ipred).EQ.BTREE_NOT_FOUND) THEN
          CALL output_line('Vertex does not exist in graph!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
          CALL sys_halt()
        END IF
        ipos   = rgraph%rVertices%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
        itable = rgraph%rVertices%p_IData(1,ipos)
        
        ! Get table for ireplaceVertex if required
        IF (bdoReplace) THEN
          IF (btree_searchInTree(rgraph%rVertices,ireplaceVertex,ipred).EQ.BTREE_NOT_FOUND) THEN
            CALL output_line('Vertex does not exist in graph!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            CALL sys_halt()
          END IF
          ipos          = rgraph%rVertices%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
          ireplaceTable = rgraph%rVertices%p_IData(1,ipos)
        ELSE
          ireplaceTable=0
        END IF
        
        ! Step 1: Delete corresponding vertex from tree. Note that for dense 
        !         graphs the vertex tree does not store further information and
        !         can be eliminated in the first step. If the vertex should be 
        !         replaced by the last vertex, then the last vertex is 
        !         removed from the tree instead of vertex numbered iVertex.
        
        IF (btree_deleteFromTree(rgraph%rVertices,&
            MERGE(ireplaceTable,iVertex,bdoReplace)).EQ.BTREE_NOT_FOUND) THEN
          CALL output_line('Vertex does not exist in graph!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
          CALL sys_halt()
        END IF
        
        ! Step 2: Loop through adjacency list of vertex iVertex and delete all 
        !         edges (jVertex,iVertex) from the adjacency lists of jVertex.
        
        ! Find position of first entry in adjacency list
        iposVertex=arrlst_getNextInArrayList(rgraph%rEdges,itable,.TRUE.)
        DO WHILE(iposVertex .NE. ARRLST_NULL)
          
          ! Get number of adjacent vertex
          jVertex=rgraph%rEdges%p_IData(iposVertex)
          
          ! Get position of next entry in adjacency list
          iposVertex=arrlst_getNextInArrayList(rgraph%rEdges,itable,.FALSE.)
          
          ! Do nothing if both vertices are the same
          IF (iVertex .EQ. jVertex) CYCLE

          ! In addition, do nothing if the current vertex is identical to 
          ! the replacement vertex. Otherwise, we would have to re-insert
          ! it afterwards. Hence, it does not make sense to remove before.
          IF (bdoReplace .AND. (ireplaceVertex .EQ. jVertex)) CYCLE

          ! Get table for jVertex
          IF (btree_searchInTree(rgraph%rVertices,jVertex,ipred).EQ.BTREE_NOT_FOUND) THEN
            CALL output_line('Vertex does not exist in graph!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            CALL sys_halt()
          END IF
          ipos   = rgraph%rVertices%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
          jtable = rgraph%rVertices%p_IData(1,ipos)

          ! Remove vertex iVertex from adjacency list of vertex jVertex
          IF (arrlst_deleteFromArrayList(rgraph%rEdges,jtable,iVertex).EQ.&
              ARRAYLIST_NOT_FOUND) THEN
            CALL output_line('Unable to delete vertex from adjacency list!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            CALL sys_halt()
          END IF
          
          ! Decrease number of edges by two; for the edge (iVertex,jVertex)
          ! and for the edge (jVertex,iVertex) that exists in an undirected graph
          rgraph%NEDGE = rgraph%NEDGE-2          
        END DO

        ! Now, vertex iVertex does no longer exist in any adjacency list.
        ! Check if replacement vertex needs to be moved to position iVertex.
        IF (bdoReplace .AND. (iVertex .NE. ireplaceVertex)) THEN
          
          ! Step 3(a): Loop through adjacency list of vertex ireplaceVertex and 
          !            delete all edges (jVertex,ireplaceVertex) from the adjacency
          !            lists of jVertex. Afterwards, add the edge (jVertex,iVertex)
          !            to the adjacency list of jVertex. Finally, swap adjacency 
          !            list of vertices iVertex and ireplaceVertex.
          
          ! Find position of first entry in adjacency list
          iposVertex=arrlst_getNextInArrayList(rgraph%rEdges,ireplaceTable,.TRUE.)
          DO WHILE(iposVertex .NE. ARRLST_NULL)

            ! Get number of adjacent vertex
            jVertex=rgraph%rEdges%p_IData(iposVertex)
            
            ! Get position of next entry in adjacency list
            iposVertex=arrlst_getNextInArrayList(rgraph%rEdges,ireplaceTable,.FALSE.)
            
            ! Do nothing if both vertices are the same. This situation required
            ! special treatment (see below)
            IF ((ireplaceVertex .EQ. jVertex) .OR. iVertex .EQ. jVertex) CYCLE
            
            ! Get table for jVertex
            IF (btree_searchInTree(rgraph%rVertices,jVertex,ipred).EQ.BTREE_NOT_FOUND) THEN
              CALL output_line('Vertex does not exist in graph!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
              CALL sys_halt()
            END IF
            ipos   = rgraph%rVertices%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
            jtable = rgraph%rVertices%p_IData(1,ipos)
            
            ! Remove vertex lastVertex from adjacency list of vertex jVertex
            IF (arrlst_deleteFromArrayList(rgraph%rEdges,jtable,ireplaceVertex).EQ.&
                ARRAYLIST_NOT_FOUND) THEN
              CALL output_line('Unable to update vertex in adjacency list!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
              CALL sys_halt()
            END IF
            
            ! Look for position of edge (jVertex,iVertex)
            IF (arrlst_searchInArrayList(rgraph%rEdges,jtable,&
                iVertex,ipred).EQ.ARRAYLIST_FOUND) THEN
              CALL output_line('Vertex replacement already exists in adjacency list!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
              CALL sys_halt()
            END IF
            
            ! Insert edge (jVertex,iVertex) into adjacency list
            CALL arrlst_insertIntoArrayList(rgraph%rEdges,jtable,iVertex,ipred,ipos)            
          END DO

          ! Remove the trivial edge (ireplaceVertex,ireplaceVertex) from the adjacency
          ! list. Note that the trivial edge (iVertex,iVertex) still exists and has
          ! not been removed in the upper removal loop
          IF (arrlst_deleteFromArrayList(rgraph%rEdges,ireplaceTable,ireplaceVertex).EQ.&
              ARRAYLIST_NOT_FOUND) THEN
            CALL output_line('Unable to update vertex in adjacency list!',&
                OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
            CALL sys_halt()
          END IF

          ! Swap adjacency list of vertices iVertex and ireplaceVertex
          CALL arrlst_swapArrayList(rgraph%rEdges,itable,ireplaceTable)

          ! Release adjacency list of vertex ireplaceVertex
          CALL arrlst_releaseArrayList(rgraph%rEdges,ireplaceTable)
          
          ! Decrease number of vertices by one
          rgraph%NVT = rgraph%NVT-1
          
          ! Decrease number of edges by one; for the trivial edge (lastVertex,lastVertex)
          rgraph%NEDGE = rgraph%NEDGE-1

        ELSE
          
          ! Step 3(b): Release adjacency list of vertex iVertex
          CALL arrlst_releaseArrayList(rgraph%rEdges,itable)
          
          ! Decrease number of vertices by one
          rgraph%NVT = rgraph%NVT-1
          
          ! Decrease number of edges by one; for the trivial edge (iVertex,iVertex)
          rgraph%NEDGE = rgraph%NEDGE-1
        END IF
        
      CASE DEFAULT
        CALL output_line('Invalid graph format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'grph_removeVertex')
        CALL sys_halt()
      END SELECT   
    END IF
  END SUBROUTINE grph_removeVertex

  ! ***************************************************************************

!<function>

  FUNCTION grph_hasEdge(rgraph,iFromVertex,iToVertex,iEdgePosition) RESULT(bexists)

!<description>
    ! This function returns TRUE if the graph has the given edge from
    ! vertex iFromVertex to vertex iToVertex.
    ! Otherwise, it returns FALSE. If the optional parameter iEdgePosition
    ! is present, then the function will also return the position of the
    ! edge (iFromVertex,iToVertex) in the array list.
!</description>

!<input>
    ! Number of the starting vertex
    INTEGER(PREC_VECIDX), INTENT(IN) :: iFromVertex

    ! Number of the ending vertex
    INTEGER(PREC_VECIDX), INTENT(IN) :: iToVertex
!</input>

!<inputoutput>
    ! The graph
    TYPE(t_graph), INTENT(INOUT) :: rgraph
!</inputoutput>

!<output>
    ! OPTIONAL: position of vertex in list
    INTEGER(PREC_TREEIDX), INTENT(OUT), OPTIONAL :: iEdgePosition
!</output>

!<result>
    ! Flag for existence of the vertex
    LOGICAL :: bexists
!</result>
!</function>

    ! local variables
    INTEGER(PREC_TABLEIDX) :: itable
    INTEGER(PREC_TREEIDX) :: ipred,ipos

    ! Determine table
    IF (rgraph%bisDense) THEN
      itable = iFromVertex
    ELSE
      IF (btree_searchInTree(rgraph%rVertices,iFromVertex,ipred).EQ.BTREE_NOT_FOUND) THEN
        bexists = .FALSE.
        RETURN
      END IF
      ipos = rgraph%rVertices%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
      itable = rgraph%rVertices%p_IData(1,ipos)
    END IF

    bexists = (arrlst_searchInArrayList(rgraph%rEdges,itable,iToVertex,ipos).EQ.ARRAYLIST_FOUND)

    ! Return edge position if required
    IF (bexists .AND. PRESENT(iEdgePosition)) THEN
      iEdgePosition = ipos
    END IF
  END FUNCTION grph_hasEdge

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE grph_insertEdge(rgraph,iFromVertex,iToVertex,iToEdgePosition,iFromEdgePosition)

!<description>
    ! This subroutine inserts the edge between vertices iFromVertex and
    ! iToVertex. If the graph is undirected, then the opposite edge
    ! (iToVertex,iFromVertex) is also added to attain symmetry of the
    ! sparsity pattern.
!</description>

!<input>
    ! Number of the starting vertex
    INTEGER(PREC_VECIDX), INTENT(IN) :: iFromVertex

    ! Number of the ending vertex
    INTEGER(PREC_VECIDX), INTENT(IN) :: iToVertex    
!</input>

!<inputoutput>
    ! The graph
    TYPE(t_graph), INTENT(INOUT) :: rgraph
!</inputoutput>

!<output>
    ! OPTIONAL: Position of edge (iToVertex,iFromVertex) in the array list
    INTEGER(PREC_ARRAYLISTIDX), INTENT(OUT), OPTIONAL :: iToEdgePosition

    ! OPTIONAL: Position of edge (iFromVertex,iToVertex) in the array list
    INTEGER(PREC_ARRAYLISTIDX), INTENT(OUT), OPTIONAL :: iFromEdgePosition
!</output>
!</subroutine>

    ! local variables
    INTEGER(PREC_TABLEIDX) :: itable
    INTEGER(PREC_TREEIDX) :: ipred,ipos
    
    ! What graph format are we?
    SELECT CASE(rgraph%cgraphFormat)

    CASE(GRPH_GRAPHUNORDERED_DIRECTED)

      IF (rgraph%bisDense) THEN
        
        ! Insert entry iToVertex into adjacency list of vertex iFromVertex
        IF (arrlst_searchInArrayList(rgraph%rEdges,iToVertex,&
            iFromVertex,ipred).EQ.ARRAYLIST_NOT_FOUND) THEN
          
          CALL arrlst_appendToArrayList(rgraph%rEdges,iToVertex,iFromVertex,ipos)
          rgraph%NEDGE = rgraph%NEDGE+1

          IF (PRESENT(iToEdgePosition)) iToEdgePosition=ipos
        END IF
      ELSE

        ! Get table associated with vertex iToVertex
        IF (btree_searchInTree(rgraph%rVertices,iToVertex,ipred).EQ.BTREE_NOT_FOUND) THEN
          CALL output_line('Vertex does not exists!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_insertEdge')
          CALL sys_halt()
        END IF
        ipos = rgraph%rVertices%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
        itable = rgraph%rVertices%p_IData(1,ipos)

        ! Insert entry iToVertex into adjacency list of vertex iFromVertex
        IF (arrlst_searchInArrayList(rgraph%rEdges,itable,&
            iFromVertex,ipred).EQ.ARRAYLIST_NOT_FOUND) THEN
          
          CALL arrlst_appendToArrayList(rgraph%rEdges,itable,iFromVertex,ipos)
          rgraph%NEDGE = rgraph%NEDGE+1

          IF (PRESENT(iToEdgePosition)) iToEdgePosition=ipos
        END IF
      END IF


    CASE(GRPH_GRAPHORDERED_DIRECTED)

      IF (rgraph%bisDense) THEN
        
        ! Insert entry iToVertex into adjacency list of vertex iFromVertex
        IF (arrlst_searchInArrayList(rgraph%rEdges,iToVertex,&
            iFromVertex,ipred).EQ.ARRAYLIST_NOT_FOUND) THEN
          
          CALL arrlst_insertIntoArrayList(rgraph%rEdges,iToVertex,iFromVertex,ipred,ipos)
          rgraph%NEDGE = rgraph%NEDGE+1

          IF (PRESENT(iToEdgePosition)) iToEdgePosition=ipos
        END IF
      ELSE

        ! Get table associated with vertex iToVertex
        IF (btree_searchInTree(rgraph%rVertices,iToVertex,ipred).EQ.BTREE_NOT_FOUND) THEN
          CALL output_line('Vertex does not exists!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_insertEdge')
          CALL sys_halt()
        END IF
        ipos = rgraph%rVertices%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
        itable = rgraph%rVertices%p_IData(1,ipos)

        ! Insert entry iToVertex into adjacency list of vertex iFromVertex
        IF (arrlst_searchInArrayList(rgraph%rEdges,itable,&
            iFromVertex,ipred).EQ.ARRAYLIST_NOT_FOUND) THEN
          
          CALL arrlst_insertIntoArrayList(rgraph%rEdges,itable,iFromVertex,ipred,ipos)
          rgraph%NEDGE = rgraph%NEDGE+1

          IF (PRESENT(iToEdgePosition)) iToEdgePosition=ipos
        END IF
      END IF

      
    CASE(GRPH_GRAPH7,GRPH_GRAPH9)
      
      IF (rgraph%bisDense) THEN

        ! Insert entry iFromVertex into adjacency list of vertex iToVertex
        IF (arrlst_searchInArrayList(rgraph%rEdges,iToVertex,&
            iFromVertex,ipred).NE.ARRAYLIST_FOUND) THEN

          CALL arrlst_insertIntoArrayList(rgraph%rEdges,iToVertex,iFromVertex,ipred,ipos)
          rgraph%NEDGE = rgraph%NEDGE+1

          IF (PRESENT(iToEdgePosition)) iToEdgePosition=ipos
        END IF
        
        ! Insert entry iToVertex into adjacency list of vertex iFromVertex
        IF (arrlst_searchInArrayList(rgraph%rEdges,iFromVertex,&
            iToVertex,ipred).EQ.ARRAYLIST_NOT_FOUND) THEN
          
          CALL arrlst_insertIntoArrayList(rgraph%rEdges,iFromVertex,iToVertex,ipred,ipos)
          rgraph%NEDGE = rgraph%NEDGE+1
          
          IF (PRESENT(iFromEdgePosition)) iFromEdgePosition=ipos
        END IF
      ELSE
        
        ! Get table associated with vertex iToVertex
        IF (btree_searchInTree(rgraph%rVertices,iToVertex,ipred).EQ.BTREE_NOT_FOUND) THEN
          CALL output_line('Vertex does not exists!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_insertEdge')
          CALL sys_halt()
        END IF
        ipos = rgraph%rVertices%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
        itable = rgraph%rVertices%p_IData(1,ipos)

        ! Insert entry iFromVertex into adjacency list of vertex iToVertex
        IF (arrlst_searchInArrayList(rgraph%rEdges,itable,&
            iFromVertex,ipred).NE.ARRAYLIST_FOUND) THEN
          
          CALL arrlst_insertIntoArrayList(rgraph%rEdges,itable,iFromVertex,ipred,ipos)
          rgraph%NEDGE = rgraph%NEDGE+1

          IF (PRESENT(iToEdgePosition)) iToEdgePosition=ipos
        END IF

        ! Get table associated with vertex iFromVertex
        IF (btree_searchInTree(rgraph%rVertices,iFromVertex,ipred).EQ.BTREE_NOT_FOUND) THEN
          CALL output_line('Vertex does not exists!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_insertEdge')
          CALL sys_halt()
        END IF
        ipos = rgraph%rVertices%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
        itable = rgraph%rVertices%p_IData(1,ipos)

        ! Insert entry iToVertex into adjacency list of vertex iFromVertex
        IF (arrlst_searchInArrayList(rgraph%rEdges,itable,&
            iToVertex,ipred).EQ.ARRAYLIST_NOT_FOUND) THEN
          
          CALL arrlst_insertIntoArrayList(rgraph%rEdges,itable,iToVertex,ipred,ipos)
          rgraph%NEDGE = rgraph%NEDGE+1

          IF (PRESENT(iFromEdgePosition)) iFromEdgePosition=ipos
        END IF    
      END IF

      
    CASE DEFAULT
      CALL output_line('Unsupported graph format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'grph_insertEdge')
      CALL sys_halt()
    END SELECT
  END SUBROUTINE grph_insertEdge

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE grph_removeEdge(rgraph,iFromVertex,iToVertex)

!<description>
    ! This subroutine removes the edge between vertices iFromVertex and
    ! iToVertex. If the graph is undirected, then the opposite edge
    ! (iToVertex,iFromVertex) is also removed to attain symmetry of the
    ! sparsity pattern.
!</description>

!<input>
    ! Number of the starting vertex
    INTEGER(PREC_VECIDX), INTENT(IN) :: iFromVertex

    ! Number of the ending vertex
    INTEGER(PREC_VECIDX), INTENT(IN) :: iToVertex    
!</input>

!<inputoutput>
    ! The graph
    TYPE(t_graph), INTENT(INOUT) :: rgraph
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_TABLEIDX) :: itable
    INTEGER(PREC_TREEIDX) :: ipred,ipos
    
    ! What graph format are we=
    SELECT CASE(rgraph%cgraphFormat)

    CASE(GRPH_GRAPHUNORDERED_DIRECTED,GRPH_GRAPHORDERED_DIRECTED)

      IF (rgraph%bisDense) THEN

        ! Remove vertex iToVertex from table iFromVertex
        IF (arrlst_deleteFromArrayList(rgraph%rEdges,iFromVertex,iToVertex)&
            .EQ.ARRAYLIST_FOUND) rgraph%NEDGE = rgraph%NEDGE-1
      ELSE

        ! Get table associated with vertex iFromVertex
        IF (btree_searchInTree(rgraph%rVertices,iFromVertex,ipred).EQ.BTREE_NOT_FOUND) THEN
          CALL output_line('Vertex does not exist!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_removeEdge')
          CALL sys_halt()
        END IF
        ipos = rgraph%rVertices%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
        itable = rgraph%rVertices%p_IData(1,ipos)
        
        ! Remove vertex iToVertex from table iFromVertex
        IF (arrlst_deleteFromArrayList(rgraph%rEdges,itable,iToVertex)&
            .EQ.ARRAYLIST_FOUND) rgraph%NEDGE = rgraph%NEDGE-1
      END IF


    CASE(GRPH_GRAPH7,GRPH_GRAPH9)

      IF (rgraph%bisDense) THEN
        
        ! Remove vertex iToVertex from table iFromVertex
        IF (arrlst_deleteFromArrayList(rgraph%rEdges,iFromVertex,iToVertex)&
            .EQ.ARRAYLIST_FOUND) rgraph%NEDGE = rgraph%NEDGE-1

        ! Remove vertex iFromVertex from table iToVertex
        IF (arrlst_deleteFromArrayList(rgraph%rEdges,iToVertex,iFromVertex)&
            .EQ.ARRAYLIST_FOUND) rgraph%NEDGE = rgraph%NEDGE-1
      ELSE

        ! Get table associated with vertex iFromVertex
        IF (btree_searchInTree(rgraph%rVertices,iFromVertex,ipred).EQ.BTREE_NOT_FOUND) THEN
          CALL output_line('Vertex does not exist!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_removeEdge')
          CALL sys_halt()
        END IF
        ipos = rgraph%rVertices%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
        itable = rgraph%rVertices%p_IData(1,ipos)
        
        ! Remove vertex iToVertex from table iFromVertex
        IF (arrlst_deleteFromArrayList(rgraph%rEdges,itable,iToVertex)&
            .EQ.ARRAYLIST_FOUND) rgraph%NEDGE = rgraph%NEDGE-1

        ! Get table associated with vertex iFromVertex
        IF (btree_searchInTree(rgraph%rVertices,iToVertex,ipred).EQ.BTREE_NOT_FOUND) THEN
          CALL output_line('Vertex does not exist!',&
              OU_CLASS_ERROR,OU_MODE_STD,'grph_removeEdge')
          CALL sys_halt()
        END IF
        ipos = rgraph%rVertices%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
        itable = rgraph%rVertices%p_IData(1,ipos)

        ! Remove vertex iFromVertex from table iToVertex
        IF (arrlst_deleteFromArrayList(rgraph%rEdges,itable,iFromVertex)&
            .EQ.ARRAYLIST_FOUND) rgraph%NEDGE = rgraph%NEDGE-1
      END IF


    CASE DEFAULT
      CALL output_line('Unsupported graph format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'grph_removeEdge')
      CALL sys_halt()
    END SELECT
  END SUBROUTINE grph_removeEdge
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE grph_infoGraph(rgraph)

!<description>
    ! This subroutine prints information about the grahp
!</description>

!<input>
    ! graph
    TYPE(t_graph), INTENT(IN) :: rgraph
!</input>
!</subroutine>

    CALL output_line('Graph statistics:')
    CALL output_line('-----------------')
    CALL output_line('cgraphFormat: '//TRIM(sys_siL(rgraph%cgraphFormat,2)))
    CALL output_line('bisDense:     '//MERGE('Yes','No ',rgraph%bisDense))
    CALL output_line('NVT:          '//TRIM(sys_siL(rgraph%NVT,15)))
    CALL output_line('NEDGE:        '//TRIM(sys_siL(rgraph%NEDGE,15)))
    CALL output_lbrk()    
  END SUBROUTINE grph_infoGraph

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE grph_duplicateGraph(rgraph,rgraphBackup)

!<description>
    ! This subroutine makes a copy of a graph in memory.
    ! It does not make sense to share some information between graphs,
    ! so each vectors is physically copied from the source graph
    ! to the destination graph.
!</description>

!<input>
    ! Source graph
    TYPE(t_graph), INTENT(IN) :: rgraph
!</input>

!<inputoutput>
    ! Destination graph
    TYPE(t_graph), INTENT(INOUT) :: rgraphBackup
!</inputoutput>
!</subroutine>

    ! Release backup graph
    CALL grph_releaseGraph(rgraphBackup)

    ! Copy all data
    rgraphBackup = rgraph

    CALL btree_duplicateTree(rgraph%rVertices,rgraphBackup%rVertices)
    CALL arrlst_duplicateArrayList(rgraph%rEdges,rgraphBackup%rEdges)
  END SUBROUTINE grph_duplicateGraph

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE grph_restoreGraph(rgraphBackup,rgraph)

!<description>
    ! This subroutine restores a graph from a previous backup.
    ! The format of both graphs must be the same.
!</description>

!<input>
    ! Backup of a graph
    TYPE(t_graph), INTENT(IN) :: rgraphBackup
!</input>

!<inputoutput>
    ! Destination graph
    TYPE(t_graph), INTENT(INOUT) :: rgraph
!</inputoutput>
!</subroutine>

    ! Check that both graphs are compatible
    IF (rgraph%cgraphFormat .NE. rgraphBackup%cgraphFormat .OR.&
        rgraph%bisDense     .NEQV. rgraphBackup%bisDense) THEN
      CALL output_line('Incompatible graphs!',&
          OU_CLASS_ERROR,OU_MODE_STD,'grph_restoreGraph')
      CALL sys_halt()
    END IF

    ! Release graph
    CALL grph_releaseGraph(rgraph)

    ! Duplicate the backup
    CALL grph_duplicateGraph(rgraphBackup,rgraph)
  END SUBROUTINE grph_restoreGraph
END MODULE graph
