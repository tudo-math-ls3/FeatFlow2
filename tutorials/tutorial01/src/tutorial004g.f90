!##############################################################################
!# Tutorial 004g: Datastructures - graph
!##############################################################################

module tutorial004g

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use graph
  use linearsystemscalar
  use matrixio

  implicit none
  private
  
  public :: start_tutorial004g

contains

  ! ***************************************************************************

  subroutine start_tutorial004g

    ! Declare some variables
    type(t_graph) :: rgraph
    type(t_matrixScalar) :: rmatrix
    integer :: i

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 004g")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Create undirected graph with space
    ! for 10 vertices and 100 edges 
    ! before new memory is allocated
    ! =================================

    call grph_createGraph(rgraph, GRPH_UNDIRECTED, 10, 100, .false.)

    ! =================================
    ! Insert 5 new vertices into the graph
    ! =================================

    do i=1,5
      call grph_insertVertex(rgraph, i)
    end do

    ! =================================
    ! Insert some edges into the graph
    ! =================================

    call grph_insertEdge(rgraph, 1,2)
    call grph_insertEdge(rgraph, 1,3)
    call grph_insertEdge(rgraph, 2,4)
    call grph_insertEdge(rgraph, 2,3)
    call grph_insertEdge(rgraph, 2,5)
    call grph_insertEdge(rgraph, 3,1)
    call grph_insertEdge(rgraph, 3,5)
    call grph_insertEdge(rgraph, 4,2)
    call grph_insertEdge(rgraph, 4,3)
    call grph_insertEdge(rgraph, 5,1)
    call grph_insertEdge(rgraph, 5,2)

    ! =================================
    ! Print graph
    ! =================================

    call grph_printGraph(rgraph)

    ! =================================
    ! Check if vertex 3 exists
    ! =================================

    if (grph_hasVertex(rgraph, 3)) then
      call output_line ("Graph has vertex 3.")
    else
      call output_line ("Graph does not have vertex 3.")
    end if

    ! =================================
    ! Check if edge (3,1) exists
    ! =================================

    if (grph_hasEdge(rgraph, 3,1)) then
      call output_line ("Graph has edge (3,1).")
    else
      call output_line ("Graph does not have edge (3,1).")
    end if

    ! =================================
    ! Create CSR9 matrix from graph
    ! =================================

    rmatrix%cmatrixFormat = LSYSSC_MATRIX9
    call lsyssc_createMatrixFromGraph(rgraph, rmatrix)

    ! =================================
    ! Write matrix to textfile
    ! =================================

    call matio_spyMatrix(&
        "post/tutorial004g_matrix","matrix",rmatrix,.false.)

    ! =================================
    ! Release graph
    ! =================================

    call grph_releaseGraph(rgraph)

    ! =================================
    ! Create graph from matrix
    ! =================================
    
    call lsyssc_createGraphFromMatrix(rmatrix, rgraph)

    ! =================================
    ! Print graph
    ! =================================

    call grph_printGraph(rgraph)

    ! =================================
    ! Release graph and matrix
    ! =================================

    call grph_releaseGraph(rgraph)
    call lsyssc_releaseMatrix(rmatrix)

  end subroutine

end module
