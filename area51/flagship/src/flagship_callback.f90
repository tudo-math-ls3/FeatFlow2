!##############################################################################
!# ****************************************************************************
!# <name> flagship_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains general callback functions which are useful
!# for all applications. In particular, the general tasks performed
!# for h-adaptivity is handled in thos module.
!#
!# The following callback functions are available:
!#
!# 1.) flagship_hadaptCallback1d
!#     -> Performs application specific tasks in the adaptation algorithm in 1D
!#
!# 2.) flagship_hadaptCallback2d
!#     -> Performs application specific tasks in the adaptation algorithm in 2D
!#
!# 3.) flagship_hadaptCallback3d
!#     -> Performs application specific tasks in the adaptation algorithm in 3D
!#
!# </purpose>
!##############################################################################

module flagship_callback

  use collection
  use fsystem
  use genoutput
  use graph
  use hadaptaux
  use linearsystemscalar
   use storage

  implicit none

  private
  public :: flagship_hadaptCallback1d
  public :: flagship_hadaptCallback2d
  public :: flagship_hadaptCallback3d

contains

  !*****************************************************************************

!<subroutine>

  subroutine flagship_hadaptCallback1D(rcollection, iOperation, Ivertices, Ielements)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 1D.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(in) :: iOperation

    ! Array of vertices involved in the adaptivity step
    integer, dimension(:), intent(in) :: Ivertices

    ! Array of elements involved in the adaptivity step
    integer, dimension(:), intent(in) :: Ielements
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_graph), pointer, save :: rgraph

    
    ! What operation should be performed?
    select case(iOperation)
      
    case(HADAPT_OPR_INITCALLBACK)
      ! This subroutine assumes that the name of the sparsity pattern
      ! is stored in the first quick access string.

      ! Retrieve sparsity pattern from collection and build sparsity-graph.
      rgraph => collct_getvalue_graph(rcollection,&
                                      trim(rcollection%SquickAccess(1)))

      
    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify sparsity graph
      nullify(rgraph)

      
    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into sparsity graph
      call grph_insertVertex(rgraph, Ivertices(1))

      
    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from sparsity graph
      if (Ivertices(2) .ne. 0) then
        call grph_removeVertex(rgraph, Ivertices(1), Ivertices(2))
      else
        call grph_removeVertex(rgraph, Ivertices(1))
      end if
      
      
    case(HADAPT_OPR_REF_LINE2LINE)
      ! Delete broken edge (I1,I2) and add new edges (I1,I3),(I2,I3)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(2))
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(3))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(3))
      
      
    case(HADAPT_OPR_CRS_2LINE1LINE)
      ! Delete broken edges (I1,I3) and (I2,I3) and add new edge (I1,I2)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(3))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(3))
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(2))
      

    case DEFAULT
      call output_line('Unsupported operation!',&
                       OU_CLASS_ERROR, OU_MODE_STD,'flagship_hadaptCallback1D')
      call sys_halt()
    end select
  end subroutine flagship_hadaptCallback1D

  !*****************************************************************************

!<subroutine>

  subroutine flagship_hadaptCallback2D(rcollection, iOperation, Ivertices, Ielements)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 2D.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(in) :: iOperation

    ! Array of vertices involved in the adaptivity step
    integer, dimension(:), intent(in) :: Ivertices

    ! Array of elements involved in the adaptivity step
    integer, dimension(:), intent(in) :: Ielements
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_graph), pointer, save :: rgraph

    
    ! What operation should be performed
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! This subroutine assumes that the name of the sparsity pattern
      ! is stored in the first quick access string.

      ! Retrieve sparsity pattern from collection and build sparsity-graph.
      rgraph => collct_getvalue_graph(rcollection,&
                                      trim(rcollection%SquickAccess(1)))


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify sparsity graph
      nullify(rgraph)
      
      
    case(HADAPT_OPR_INSERTVERTEXEDGE,&
         HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into sparsity graph
      call grph_insertVertex(rgraph, Ivertices(1))

   
    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from sparsity graph and solution
      if (Ivertices(2) .ne. 0) then
        call grph_removeVertex(rgraph,Ivertices(1), Ivertices(2))
      else
        call grph_removeVertex(rgraph, Ivertices(1))
      end if


    case(HADAPT_OPR_REF_TRIA2TRIA)
      ! Delete broken edge (I1,I2) and add three new edges 
      ! (I1,I4), (I2,I4), and (I3,I4) if this is necessary
      if (Ielements(1) .eq. Ielements(4)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(2))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(4))
      end if
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(4))


    case(HADAPT_OPR_REF_TRIA3TRIA12)
      ! Delete broken edge (I1,I2) and add new edges (I1,I4),(I2,I4)
      if (Ielements(1) .eq. Ielements(4)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(2))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(4))
      end if

      ! Delete broken edge (I2,I3) and add new edges (I2,I5),(I3,I5)
      if (Ielements(2) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(3))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(5))
      end if

      ! Add new edges (I3,I4),(I4,I5)
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(4))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(5))


    case(HADAPT_OPR_REF_TRIA3TRIA23)
      ! Delete broken edge (I1,I2) and add new edges (I1,I4),(I2,I4
      if (Ielements(1) .eq. Ielements(4)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(2))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(4))
      end if

      ! Delete broken edges (I2,I3) and add new edges (I2,I5),(I3,I5)
      if (Ielements(2) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(3))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(5))
      end if

      ! Add new edges (I1,I5),(I4,I5)
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(5))


    case(HADAPT_OPR_REF_TRIA4TRIA)
      ! Delete broken edge (I1,I2) and add new edges (I1,I4),(I2,I4)
      if (Ielements(1) .eq. Ielements(4)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(2))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(4))
      end if

      ! Delete broken edge (I2,I3) and add new edges (I2,I5),(I3,I5)
      if (Ielements(2) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(3))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(5))
      end if

      ! Delete broken edge (I1,I3) and add new edges (I3,I6),(I1,I6)
      if (Ielements(3) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(3))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(6))
      end if

      ! Add new edges I4,I5),(I4,I6), and (I5,I6)
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(6))


    case(HADAPT_OPR_REF_QUAD2QUAD)
      ! Delete broken edges (I1,I3),(I2,I4)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(3))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(4))

      ! Delete broken edge (I1,I2) and add new edges (I1,I5),(I2,I5)
      if (Ielements(1) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(2))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(5))
      end if

      ! Delete broken edge (I3,I4) and add new edges (I3,I6),(I4,I6)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(4), Ivertices(6))
      end if

      ! Add new edges (I1,I6),(I4,I5),(I2,I6),(I3,I5),(I5,I6)
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(6))


    case(HADAPT_OPR_REF_QUAD3TRIA)
      ! Delete broken edges (I1,I3),(I2,I4)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(3))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(4))

      ! Delete broken edge (I1,I2) and add new edges (I1,I5),(I2,I5)
      if (Ielements(1) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(2))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(5))
      end if

      ! Add new edges (I3,I5),(I4,I5)
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(5))


    case(HADAPT_OPR_REF_QUAD4TRIA)
      ! Delete broken edges (I1,I3),(I2,I4)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(3))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(4))

      ! Delete broken edge (I1,I2) and add new edges (I1,I5),(I2,I5)
      if (Ielements(1) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(2))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(5))
      end if

      ! Delete broken edge (I2,I3) and add new edges (I2,I6),(I3,I6)
      if (Ielements(2) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(3))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(6))
      end if

      ! Add new edges (I4,I5),(I4,I6), and (I5,I6)
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(6))


    case(HADAPT_OPR_REF_QUAD4QUAD)
      ! Delete broken edges (I1,I3) and (I2,I4)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(3))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(4))

      ! Delete broken edge (I1,I2) and add new edges (I1,I5),(I2,I5)
      if (Ielements(1) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(2))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(5))
      end if

      ! Delete broken edge (I2,I3) and add new edges (I2,I6),(I3,I6)
      if (Ielements(2) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(3))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(6))
      end if

      ! Delete broken edge (I3,I4) and add new edges (I3,I7),(I4,I7)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(7))
        call grph_insertEdge(rgraph, Ivertices(4), Ivertices(7))
      end if

      ! Delete broken edge (I1,I4) and add new edges (I4,I8),(I1,I8)
      if (Ielements(4) .eq. Ielements(8)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(4), Ivertices(8))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(8))
      end if

      ! Add new edges (I5,I9),(I6,I9),(I7,I9),(I8,I9),(I1,I9),(I5,I8),
      ! (I2,I9),(I5,I6),(I3,I9),(I6,I7),(I4,I9), and (I7,I8)
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(6), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(7), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(8), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(8))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(6), Ivertices(7))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(7), Ivertices(8))


    case(HADAPT_OPR_CVT_TRIA2TRIA)
      ! Delete broken edge (I2,I3) and add new edges (I2,I5),(I3,I5)
      if (Ielements(2) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(3))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(5))
      end if

      ! Delete broken edge (I1,I3) and add new edges (I1,I6),(I3,I6)
      if (Ielements(3) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(3))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(6))
      end if

      ! Delete broken edge (I3,I4) and add new edges (I5,I6),(I4,I6),(I4,I5)
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(4))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(6))


    case(HADAPT_OPR_CVT_QUAD2QUAD)

      ! Delete broken edge (I2,I3) and add new edges (I2,I6),(I3,I6)
      if (Ielements(2) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(3))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(6))
      end if

      ! Delete broken edge (I1,I4) and add new edges  (I1,I8),(I4,I8)
      if (Ielements(4) .eq. Ielements(8)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(8))
        call grph_insertEdge(rgraph, Ivertices(4), Ivertices(8))
      end if

      ! Delete broken edges (I5,I7),(I2,I7),(I3,I5),(I1,I7),(I4,I5) and
      ! add new edges (I5,I9),(I6,I9),(I7,I9),(I8,I9),(I1,I9),(I2,I9),
      ! (I3,I9),(I4,I9),(I5,I8),(I5,I6),(I6,I7),(I7,I8)
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(5))
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(6), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(7), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(8), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(8))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(6), Ivertices(7))
      call grph_insertEdge(rgraph, Ivertices(7), Ivertices(8))


    case(HADAPT_OPR_CVT_QUAD3TRIA)
      ! Delete broken edge (I2,I3) and add new edges (I2,I6),(I3,I6)
      if (Ielements(2) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(3))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(6))
      end if

      ! Delete broken edge (I3,I4) and add new edges (I3,I7),(I4,I7)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(7))
        call grph_insertEdge(rgraph, Ivertices(4), Ivertices(7))
      end if

      ! Delete broken edge (I1,I4) and add new edges (I4,I8),(I1,I8)
      if (Ielements(4) .eq. Ielements(8)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(4), Ivertices(8))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(8))
      end if

      ! Delete broken edges (I5,I3) and (I5,I4) and add new edges
      ! (I5,I9),(I6,I9),(I7,I9),(I8,I9),(I1,I9),(I5,I8),(I2,I9),(I5,I6)
      ! (I3,I9),I6,I7),(I4,I9), and (I7,I8)
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(5))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(5))
      
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(6), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(7), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(8), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(6), Ivertices(7))
      call grph_insertEdge(rgraph, Ivertices(7), Ivertices(8))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(8))

    case(HADAPT_OPR_CVT_QUAD4TRIA)
      ! Delete broken edge (I3,I4) and add new edges (I3,I7),(I4,I7)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(7))
        call grph_insertEdge(rgraph, Ivertices(4), Ivertices(7))
      end if

      ! Delete broken edge (I1,I4) and add new edges (I4,I8),(I1,I8)
      if (Ielements(4) .eq. Ielements(8)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(4), Ivertices(8))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(8))
      end if

      ! Delete broken edges (I4,I5),(I5,I6) and (I4,I6) and add new
      ! edges (I5,I9),(I6,I9),(I7,I9),(I8,I9),(I1,I9),(I5,I8),(I2,I9),
      ! (I5,I6),(I3,I9),I6,I7),(I4,I9), and (I7,I8)
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(6), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(7), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(8), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(8))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(6), Ivertices(7))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(9))
      call grph_insertEdge(rgraph, Ivertices(7), Ivertices(8))   


    case(HADAPT_OPR_CRS_2TRIA1TRIA)
      ! Delete broken edges (I1,I4), and (I2,I4) and add new edge (I1,I2)
      if (Ielements(1) .eq. Ielements(4)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(2))
      end if
      
      ! Delete broken edge (I3,I4)
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(4))
      

    case(HADAPT_OPR_CRS_4TRIA1TRIA)
      ! Delete broken edges (I4,I5),(I4,I6), and (I5,I6)
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(6))

      ! Delete broken edges (I1,I4), and (I2,I4) and add new edge (I1,I2)
      if (Ielements(1) .eq. Ielements(4)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(2))
      end if

      ! Delete broken edges (I2,I5), and (I3,I5) and add new edge (I2,I3)
      if (Ielements(2) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(5))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(3))
      end if

      ! Delete broken edges (I1,I6), and (I3,I6) and add new edge (I1,I3)
      if (Ielements(3) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(6))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(3))
      end if


    case(HADAPT_OPR_CRS_4TRIA2TRIA1)
      ! Delete broken edges (I4,I5),(I5,I6), and (I4,I6) and add new edge (I3,I4)
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(4))

      ! Delete broken edges (I2,I5), and (I3,I5) and add new edge (I2,I3)
      if (Ielements(2) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(5))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(3))
      end if

      ! Delete broken edges (I1,I6), and (I3,I6) and add new edge (I1,I3)
      if (Ielements(3) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(6))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(3))
      end if


    case(HADAPT_OPR_CRS_4TRIA2TRIA2)
      ! Delete broken edges (I4,I5),(I5,I6), and (I4,I6) and add new edge (I1,I5)
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(5))

      ! Delete broken edges (I1,I4), and (I2,I4) and add new edge (I1,I4)
      if (Ielements(1) .eq. Ielements(4)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(2))
      end if

      ! Delete broken edges (I1,I6), and (I3,I6) and add new edge (I1,I3)
      if (Ielements(3) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(6))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(3))
      end if


    case(HADAPT_OPR_CRS_4TRIA2TRIA3)
      ! Delete broken edges (I4,I5),(I5,I6), and (I4,I6) and add new edge (I2,I6)
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(6))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(6))

      ! Delete broken edges (I1,I4), and (I2,I4) and add new edge (I1,I4)
      if (Ielements(1) .eq. Ielements(4)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(4))
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(4))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(2))
      end if

      ! Delete broken edges (I2,I5), and (I3,I5) and add new edge (I2,I3)
      if (Ielements(2) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(5))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(3))
      end if


    case(HADAPT_OPR_CRS_4QUAD1QUAD)
      ! Delete broken edges (I1,I9),(I2,I9),(I3,I9),(I4,I9),(I5,I9),(I6,I9),
      ! (I7,I9),(I8,I9),(I5,I6),(I6,I7),(I7,I8), and (I5,I8)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(7), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(8), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(7), Ivertices(8))
      call grph_removeEdge(rgraph, Ivertices(8), Ivertices(5))
     
      ! Delete broken edges (I1,I5) and (I2,I5) and add new edge (I1,I2)
      if (Ielements(1) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(5))
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(2))
      end if

      ! Delete broken edges (I2,I6) and (I3,I6) and add new edge (I2,I3)
      if (Ielements(2) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(6))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(3))
      end if

      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(7))
        call grph_removeEdge(rgraph, Ivertices(4), Ivertices(7))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(4))
      end if

      ! Delete broken edges (I4,I8) and (I1,I8) and add new edge (I4,I1)
      if (Ielements(4) .eq. Ielements(8)) then
        call grph_removeEdge(rgraph, Ivertices(4), Ivertices(8))
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(8))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(4))
      end if
      
      ! Add new edges (I1,I3) and (I2,I4)
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(3))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(4))
 
      
    case(HADAPT_OPR_CRS_4QUAD2QUAD)
      ! Delete broken edges (I1,I9),(I2,I9),(I3,I9),(I4,I9),(I5,I9),(I6,I9),
      ! (I7,I9),(I8,I9),(I5,I6),(I6,I7),(I7,I8), and (I5,I8)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(7), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(8), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(7), Ivertices(8))
      call grph_removeEdge(rgraph, Ivertices(8), Ivertices(5))
      
      ! Delete broken edges (I2,I6) and (I3,I6) and add new edge (I2,I3)
      if (Ielements(2) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(6))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(3))
      end if
      
      ! Delete broken edges (I4,I8) and (I1,I8) and add new edge (I4,I1)
      if (Ielements(4) .eq. Ielements(8)) then
        call grph_removeEdge(rgraph, Ivertices(4), Ivertices(8))
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(8))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(4))
      end if
      
      ! Add new edges (I5,I7),(I1,I7),(I4,I5),(I2,I7) and (I3,I5)
      call grph_insertEdge(rgraph, Ivertices(5), Ivertices(7))
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(7))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(7))
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(5))
      
      
    case(HADAPT_OPR_CRS_4QUAD3TRIA)
      ! Delete broken edges (I1,I9),(I2,I9),(I3,I9),(I4,I9),(I5,I9),(I6,I9),
      ! (I7,I9),(I8,I9),(I5,I6),(I6,I7),(I7,I8), and (I5,I8)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(7), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(8), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(7), Ivertices(8))
      call grph_removeEdge(rgraph, Ivertices(8), Ivertices(5))
      
      ! Delete broken edges (I2,I6) and (I3,I6) and add new edge (I2,I3)
      if (Ielements(2) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(6))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(3))
      end if
      
      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(7))
        call grph_removeEdge(rgraph, Ivertices(4), Ivertices(7))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(4))
      end if
      
      ! Delete broken edges (I4,I8) and (I1,I8) and add new edge (I4,I1)
      if (Ielements(4) .eq. Ielements(8)) then
        call grph_removeEdge(rgraph, Ivertices(4), Ivertices(8))
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(8))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(4))
      end if
      
      ! Add new edges (I3,I5) and (I4,I5)
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(5))
      
      
    case(HADAPT_OPR_CRS_4QUAD4TRIA)
      ! Delete broken edges (I1,I9),(I2,I9),(I3,I9),(I4,I9),(I5,I9),(I6,I9),
      ! (I7,I9),(I8,I9),(I5,I6),(I6,I7), and (I7,I8)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(7), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(8), Ivertices(9))
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(7), Ivertices(8))
      
      ! Delete broken edges (I2,I6) and (I3,I6) and add new edge (I2,I3)
      if (Ielements(2) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(6))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(3))
      end if
      
      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(7))
        call grph_removeEdge(rgraph, Ivertices(4), Ivertices(7))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(4))
      end if
      
      ! Add new edges (I3,I5) and (I3,I8)
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(8))
      
      
    case(HADAPT_OPR_CRS_2QUAD1QUAD)
      ! Delete broken edges (I5,I7), (I1,I7),(I4,I5),(I2,(I7) and (I3,I5)
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(5))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(5))
      
      ! Delete broken edges (I1,I5) and (I2,I5) and add new edge (I1,I2)
      if (Ielements(1) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(5))
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(2))
      end if
      
      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(7))
        call grph_removeEdge(rgraph, Ivertices(4), Ivertices(7))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(4))
      end if
      
      ! Add new edges (I1,I3) and (I2,I4)
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(3))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(4))
      
      
    case(HADAPT_OPR_CRS_2QUAD3TRIA)
      ! Delete broken edges (I5,I7),(I1,I7),(I4,I5),(I2,(I7) and (I3,I5)
      call grph_removeEdge(rgraph, Ivertices(5), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(2), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(5))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(5))
      
      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(7))
        call grph_removeEdge(rgraph, Ivertices(4), Ivertices(7))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(4))
      end if
      
      ! Add new edges (I3,I5) and (I4,I5)
      call grph_insertEdge(rgraph, Ivertices(3), Ivertices(5))
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(5))
      
      
    case(HADAPT_OPR_CRS_3TRIA1QUAD)
      ! Delete broken edges (I3,I5) and (I4,I5)
      call grph_removeEdge(rgraph, Ivertices(3), Ivertices(5))
      call grph_removeEdge(rgraph, Ivertices(4), Ivertices(5))
      
      ! Delete broken edges (I1,I5) and (I2,I5) and add new edge (I1,I2)
      if (Ielements(1) .eq. Ielements(5)) then
        call grph_removeEdge(rgraph, Ivertices(1), Ivertices(5))
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(5))
        call grph_insertEdge(rgraph, Ivertices(1), Ivertices(2))
      end if
      
      ! Add new edges (I1,I3) and (I2,I4)
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(3))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(4))
      
      
    case(HADAPT_OPR_CRS_4TRIA1QUAD)
      ! Delete broken edges (I1,I6),(I1,I7) and (I6,I7)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(7))
      
      ! Delete broken edges (I2,I6) and (I3,I6) and add new edge (I2,I3)
      if (Ielements(2) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(6))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(3))
      end if
      
      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(7))
        call grph_removeEdge(rgraph, Ivertices(4), Ivertices(7))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(4))
      end if
      
      ! Add new edges (I1,I3) and (I2,I4)
      call grph_insertEdge(rgraph, Ivertices(1), Ivertices(3))
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(4))
      
      
    case(HADAPT_OPR_CRS_4TRIA3TRIA2)
      ! Delete broken edges (I1,I7) and (I6,I7)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(7))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(7))
      
      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (Ielements(3) .eq. Ielements(7)) then
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(7))
        call grph_removeEdge(rgraph, Ivertices(4), Ivertices(7))
        call grph_insertEdge(rgraph, Ivertices(3), Ivertices(4))
      end if
      
      ! Add new edge (I4,I6)
      call grph_insertEdge(rgraph, Ivertices(4), Ivertices(6))
      
      
    case(HADAPT_OPR_CRS_4TRIA3TRIA3)
      ! Delete broken edges (I1,I6) and (I6,I7)
      call grph_removeEdge(rgraph, Ivertices(1), Ivertices(6))
      call grph_removeEdge(rgraph, Ivertices(6), Ivertices(7))
      
      ! Delete broken edges (I2,I6) and (I3,I6) and add new edge (I2,I3)
      if (Ielements(2) .eq. Ielements(6)) then
        call grph_removeEdge(rgraph, Ivertices(2), Ivertices(6))
        call grph_removeEdge(rgraph, Ivertices(3), Ivertices(6))
        call grph_insertEdge(rgraph, Ivertices(2), Ivertices(3))
      end if
      
      ! Add new edge (I2,I7)
      call grph_insertEdge(rgraph, Ivertices(2), Ivertices(7))


    case DEFAULT
      call output_line('Unsupported operation!',&
                       OU_CLASS_ERROR, OU_MODE_STD,'flagship_hadaptCallback2D')
      call sys_halt()
    end select
  end subroutine flagship_hadaptCallback2D

  !*****************************************************************************

!<subroutine>

  subroutine flagship_hadaptCallback3D(rcollection, iOperation, Ivertices, Ielements)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 3D.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(in) :: iOperation

    ! Array of vertices involved in the adaptivity step
    integer, dimension(:), intent(in) :: Ivertices

    ! Array of elements involved in the adaptivity step
    integer, dimension(:), intent(in) :: Ielements
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_graph), pointer, save :: rgraph
    

    ! What operation should be performed?
    select case(iOperation)
      
    case(HADAPT_OPR_INITCALLBACK)
      ! This subroutine assumes that the name of the sparsity pattern
      ! is stored in the first quick access string.

      ! Retrieve sparsity pattern from collection and build sparsity-graph.
      rgraph => collct_getvalue_graph(rcollection,&
                                      trim(rcollection%SquickAccess(1)))

      
    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify sparsity graph
      nullify(rgraph)


    case(HADAPT_OPR_INSERTVERTEXEDGE,&
         HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into sparsity graph
      call grph_insertVertex(rgraph, Ivertices(1))

   
    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from sparsity graph and solution
      if (Ivertices(2) .ne. 0) then
        call grph_removeVertex(rgraph,Ivertices(1), Ivertices(2))
      else
        call grph_removeVertex(rgraph, Ivertices(1))
      end if

      
    case DEFAULT
      call output_line('Unsupported operation!',&
                       OU_CLASS_ERROR, OU_MODE_STD,'flagship_hadaptCallback3D')
      call sys_halt()
    end select
  end subroutine flagship_hadaptCallback3D

end module flagship_callback
