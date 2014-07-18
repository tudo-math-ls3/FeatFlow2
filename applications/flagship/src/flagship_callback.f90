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

#include "flagship.h"

!$ use omp_lib
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

  subroutine flagship_hadaptCallback1D(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 1D.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(in) :: iOperation
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_graph), pointer, save :: rgraph
    character(len=SYS_STRLEN) :: ssectionName
    integer :: i1,i2,i3


    ! What operation should be performed?
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! This subroutine assumes that the name of the sparsity pattern
      ! is stored in the first quick access string.

      ! Get section name
      call collct_getvalue_string(rcollection,&
          'ssectionname', ssectionName)

      ! Retrieve sparsity pattern from collection and build sparsity-graph.
      rgraph => collct_getvalue_graph(rcollection,&
          trim(rcollection%SquickAccess(1)), ssectionName=ssectionName)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify sparsity graph
      nullify(rgraph)


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into sparsity graph
      call grph_insertVertex(rgraph, rcollection%IquickAccess(1))


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from sparsity graph
      if (rcollection%IquickAccess(2) .ne. 0) then
        call grph_removeVertex(rgraph, rcollection%IquickAccess(1),&
                               rcollection%IquickAccess(2))
      else
        call grph_removeVertex(rgraph, rcollection%IquickAccess(1))
      end if


    case(HADAPT_OPR_REF_LINE2LINE)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)

      ! Delete broken edge (I1,I2) and add new edges (I1,I3),(I2,I3)
      call grph_removeEdge(rgraph, i1, i2)
      call grph_insertEdge(rgraph, i1, i3)
      call grph_insertEdge(rgraph, i2, i3)


    case(HADAPT_OPR_CRS_2LINE1LINE)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)

      ! Delete broken edges (I1,I3) and (I2,I3) and add new edge (I1,I2)
      call grph_removeEdge(rgraph, i1, i3)
      call grph_removeEdge(rgraph, i2, i3)
      call grph_insertEdge(rgraph, i1, i2)


    case default
      call output_line('Unsupported operation!',&
                       OU_CLASS_ERROR, OU_MODE_STD, 'flagship_hadaptCallback1D')
      call sys_halt()
    end select
  end subroutine flagship_hadaptCallback1D

  !*****************************************************************************

!<subroutine>

  subroutine flagship_hadaptCallback2D(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 2D.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(in) :: iOperation
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_graph), pointer, save :: rgraph => null()
    character(len=SYS_STRLEN) :: ssectionName
    integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9
    integer :: e1,e2,e3,e4,e5,e6,e7,e8


    ! What operation should be performed
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! This subroutine assumes that the name of the sparsity pattern
      ! is stored in the first quick access string.

      ! Get section name
      call collct_getvalue_string(rcollection,&
          'ssectionname', ssectionName)

      ! Retrieve sparsity pattern from collection and build sparsity-graph.
      rgraph => collct_getvalue_graph(rcollection,&
          trim(rcollection%SquickAccess(1)), ssectionName=ssectionName)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify sparsity graph
      nullify(rgraph)


    case(HADAPT_OPR_INSERTVERTEXEDGE,&
         HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into sparsity graph
      call grph_insertVertex(rgraph, rcollection%IquickAccess(1))


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from sparsity graph and solution
      if (rcollection%IquickAccess(2) .ne. 0) then
        call grph_removeVertex(rgraph, rcollection%IquickAccess(1),&
                               rcollection%IquickAccess(2))
      else
        call grph_removeVertex(rgraph, rcollection%IquickAccess(1))
      end if


    case(HADAPT_OPR_REF_TRIA2TRIA)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)

      e1 = rcollection%IquickAccess(5)
      e4 = rcollection%IquickAccess(8)

      ! Delete broken edge (I1,I2) and add three new edges
      ! (I1,I4), (I2,I4), and (I3,I4) if this is necessary
      if (e1 .eq. e4) then
        call grph_removeEdge(rgraph, i1, i2)
        call grph_insertEdge(rgraph, i1, i4)
        call grph_insertEdge(rgraph, i2, i4)
      end if
      call grph_insertEdge(rgraph, i3, i4)


    case(HADAPT_OPR_REF_TRIA3TRIA12)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)

      e1 = rcollection%IquickAccess(6)
      e2 = rcollection%IquickAccess(7)
      e3 = rcollection%IquickAccess(8)
      e4 = rcollection%IquickAccess(9)
      e5 = rcollection%IquickAccess(10)


      ! Delete broken edge (I1,I2) and add new edges (I1,I4),(I2,I4)
      if (e1 .eq. e4) then
        call grph_removeEdge(rgraph, i1, i2)
        call grph_insertEdge(rgraph, i1, i4)
        call grph_insertEdge(rgraph, i2, i4)
      end if

      ! Delete broken edge (I2,I3) and add new edges (I2,I5),(I3,I5)
      if (e2 .eq. e5) then
        call grph_removeEdge(rgraph, i2, i3)
        call grph_insertEdge(rgraph, i2, i5)
        call grph_insertEdge(rgraph, i3, i5)
      end if

      ! Add new edges (I3,I4),(I4,I5)
      call grph_insertEdge(rgraph, i3, i4)
      call grph_insertEdge(rgraph, i4, i5)


    case(HADAPT_OPR_REF_TRIA3TRIA23)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)

      e1 = rcollection%IquickAccess(6)
      e2 = rcollection%IquickAccess(7)
      e3 = rcollection%IquickAccess(8)
      e4 = rcollection%IquickAccess(9)
      e5 = rcollection%IquickAccess(10)

      ! Delete broken edge (I1,I2) and add new edges (I1,I4),(I2,I4
      if (e1 .eq. e4) then
        call grph_removeEdge(rgraph, i1, i2)
        call grph_insertEdge(rgraph, i1, i4)
        call grph_insertEdge(rgraph, i2, i4)
      end if

      ! Delete broken edges (I2,I3) and add new edges (I2,I5),(I3,I5)
      if (e2 .eq. e5) then
        call grph_removeEdge(rgraph, i2, i3)
        call grph_insertEdge(rgraph, i2, i5)
        call grph_insertEdge(rgraph, i3, i5)
      end if

      ! Add new edges (I1,I5),(I4,I5)
      call grph_insertEdge(rgraph, i4, i5)
      call grph_insertEdge(rgraph, i1, i5)


    case(HADAPT_OPR_REF_TRIA4TRIA)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)

      e1 = rcollection%IquickAccess(7)
      e2 = rcollection%IquickAccess(8)
      e3 = rcollection%IquickAccess(9)
      e4 = rcollection%IquickAccess(10)
      e5 = rcollection%IquickAccess(11)
      e6 = rcollection%IquickAccess(12)

      ! Delete broken edge (I1,I2) and add new edges (I1,I4),(I2,I4)
      if (e1 .eq. e4) then
        call grph_removeEdge(rgraph, i1, i2)
        call grph_insertEdge(rgraph, i1, i4)
        call grph_insertEdge(rgraph, i2, i4)
      end if

      ! Delete broken edge (I2,I3) and add new edges (I2,I5),(I3,I5)
      if (e2 .eq. e5) then
        call grph_removeEdge(rgraph, i2, i3)
        call grph_insertEdge(rgraph, i2, i5)
        call grph_insertEdge(rgraph, i3, i5)
      end if

      ! Delete broken edge (I1,I3) and add new edges (I3,I6),(I1,I6)
      if (e3 .eq. e6) then
        call grph_removeEdge(rgraph, i1, i3)
        call grph_insertEdge(rgraph, i3, i6)
        call grph_insertEdge(rgraph, i1, i6)
      end if

      ! Add new edges I4,I5),(I4,I6), and (I5,I6)
      call grph_insertEdge(rgraph, i4, i5)
      call grph_insertEdge(rgraph, i4, i6)
      call grph_insertEdge(rgraph, i5, i6)


    case(HADAPT_OPR_REF_QUAD2QUAD)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)

      e1 = rcollection%IquickAccess(7)
      e2 = rcollection%IquickAccess(8)
      e3 = rcollection%IquickAccess(9)
      e4 = rcollection%IquickAccess(10)
      e5 = rcollection%IquickAccess(11)
      e6 = rcollection%IquickAccess(12)
      e7 = rcollection%IquickAccess(13)
      e8 = rcollection%IquickAccess(14)

      ! Delete broken edges (I1,I3),(I2,I4)
      call grph_removeEdge(rgraph, i1, i3)
      call grph_removeEdge(rgraph, i2, i4)

      ! Delete broken edge (I1,I2) and add new edges (I1,I5),(I2,I5)
      if (e1 .eq. e5) then
        call grph_removeEdge(rgraph, i1, i2)
        call grph_insertEdge(rgraph, i1, i5)
        call grph_insertEdge(rgraph, i2, i5)
      end if

      ! Delete broken edge (I3,I4) and add new edges (I3,I6),(I4,I6)
      if (e3 .eq. e7) then
        call grph_removeEdge(rgraph, i3, i4)
        call grph_insertEdge(rgraph, i3, i6)
        call grph_insertEdge(rgraph, i4, i6)
      end if

      ! Add new edges (I1,I6),(I4,I5),(I2,I6),(I3,I5),(I5,I6)
      call grph_insertEdge(rgraph, i1, i6)
      call grph_insertEdge(rgraph, i4, i5)
      call grph_insertEdge(rgraph, i2, i6)
      call grph_insertEdge(rgraph, i3, i5)
      call grph_insertEdge(rgraph, i5, i6)


    case(HADAPT_OPR_REF_QUAD3TRIA)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)

      e1 = rcollection%IquickAccess(6)
      e2 = rcollection%IquickAccess(7)
      e3 = rcollection%IquickAccess(8)
      e4 = rcollection%IquickAccess(9)
      e5 = rcollection%IquickAccess(10)
      e6 = rcollection%IquickAccess(11)
      e7 = rcollection%IquickAccess(12)
      e8 = rcollection%IquickAccess(13)

      ! Delete broken edges (I1,I3),(I2,I4)
      call grph_removeEdge(rgraph, i1, i3)
      call grph_removeEdge(rgraph, i2, i4)

      ! Delete broken edge (I1,I2) and add new edges (I1,I5),(I2,I5)
      if (e1 .eq. e5) then
        call grph_removeEdge(rgraph, i1, i2)
        call grph_insertEdge(rgraph, i1, i5)
        call grph_insertEdge(rgraph, i2, i5)
      end if

      ! Add new edges (I3,I5),(I4,I5)
      call grph_insertEdge(rgraph, i3, i5)
      call grph_insertEdge(rgraph, i4, i5)


    case(HADAPT_OPR_REF_QUAD4TRIA)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)

      e1 = rcollection%IquickAccess(7)
      e2 = rcollection%IquickAccess(8)
      e3 = rcollection%IquickAccess(9)
      e4 = rcollection%IquickAccess(10)
      e5 = rcollection%IquickAccess(11)
      e6 = rcollection%IquickAccess(12)
      e7 = rcollection%IquickAccess(13)
      e8 = rcollection%IquickAccess(14)

      ! Delete broken edges (I1,I3),(I2,I4)
      call grph_removeEdge(rgraph, i1, i3)
      call grph_removeEdge(rgraph, i2, i4)

      ! Delete broken edge (I1,I2) and add new edges (I1,I5),(I2,I5)
      if (e1 .eq. e5) then
        call grph_removeEdge(rgraph, i1, i2)
        call grph_insertEdge(rgraph, i1, i5)
        call grph_insertEdge(rgraph, i2, i5)
      end if

      ! Delete broken edge (I2,I3) and add new edges (I2,I6),(I3,I6)
      if (e2 .eq. e6) then
        call grph_removeEdge(rgraph, i2, i3)
        call grph_insertEdge(rgraph, i2, i6)
        call grph_insertEdge(rgraph, i3, i6)
      end if

      ! Add new edges (I4,I5),(I4,I6), and (I5,I6)
      call grph_insertEdge(rgraph, i4, i5)
      call grph_insertEdge(rgraph, i4, i6)
      call grph_insertEdge(rgraph, i5, i6)


    case(HADAPT_OPR_REF_QUAD4QUAD)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)
      i7 = rcollection%IquickAccess(7)
      i8 = rcollection%IquickAccess(8)
      i9 = rcollection%IquickAccess(9)

      e1 = rcollection%IquickAccess(10)
      e2 = rcollection%IquickAccess(11)
      e3 = rcollection%IquickAccess(12)
      e4 = rcollection%IquickAccess(13)
      e5 = rcollection%IquickAccess(14)
      e6 = rcollection%IquickAccess(15)
      e7 = rcollection%IquickAccess(16)
      e8 = rcollection%IquickAccess(17)

      ! Delete broken edges (I1,I3) and (I2,I4)
      call grph_removeEdge(rgraph, i1, i3)
      call grph_removeEdge(rgraph, i2, i4)

      ! Delete broken edge (I1,I2) and add new edges (I1,I5),(I2,I5)
      if (e1 .eq. e5) then
        call grph_removeEdge(rgraph, i1, i2)
        call grph_insertEdge(rgraph, i1, i5)
        call grph_insertEdge(rgraph, i2, i5)
      end if

      ! Delete broken edge (I2,I3) and add new edges (I2,I6),(I3,I6)
      if (e2 .eq. e6) then
        call grph_removeEdge(rgraph, i2, i3)
        call grph_insertEdge(rgraph, i2, i6)
        call grph_insertEdge(rgraph, i3, i6)
      end if

      ! Delete broken edge (I3,I4) and add new edges (I3,I7),(I4,I7)
      if (e3 .eq. e7) then
        call grph_removeEdge(rgraph, i3, i4)
        call grph_insertEdge(rgraph, i3, i7)
        call grph_insertEdge(rgraph, i4, i7)
      end if

      ! Delete broken edge (I1,I4) and add new edges (I4,I8),(I1,I8)
      if (e4 .eq. e8) then
        call grph_removeEdge(rgraph, i1, i4)
        call grph_insertEdge(rgraph, i4, i8)
        call grph_insertEdge(rgraph, i1, i8)
      end if

      ! Add new edges (I5,I9),(I6,I9),(I7,I9),(I8,I9),(I1,I9),(I5,I8),
      ! (I2,I9),(I5,I6),(I3,I9),(I6,I7),(I4,I9), and (I7,I8)
      call grph_insertEdge(rgraph, i5, i9)
      call grph_insertEdge(rgraph, i6, i9)
      call grph_insertEdge(rgraph, i7, i9)
      call grph_insertEdge(rgraph, i8, i9)
      call grph_insertEdge(rgraph, i1, i9)
      call grph_insertEdge(rgraph, i5, i8)
      call grph_insertEdge(rgraph, i2, i9)
      call grph_insertEdge(rgraph, i5, i6)
      call grph_insertEdge(rgraph, i3, i9)
      call grph_insertEdge(rgraph, i6, i7)
      call grph_insertEdge(rgraph, i4, i9)
      call grph_insertEdge(rgraph, i7, i8)


    case(HADAPT_OPR_CVT_TRIA2TRIA)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)

      e1 = rcollection%IquickAccess(7)
      e2 = rcollection%IquickAccess(8)
      e3 = rcollection%IquickAccess(9)
      e4 = rcollection%IquickAccess(10)
      e5 = rcollection%IquickAccess(11)
      e6 = rcollection%IquickAccess(12)
      e7 = rcollection%IquickAccess(13)
      e8 = rcollection%IquickAccess(14)

      ! Delete broken edge (I2,I3) and add new edges (I2,I5),(I3,I5)
      if (e2 .eq. e5) then
        call grph_removeEdge(rgraph, i2, i3)
        call grph_insertEdge(rgraph, i2, i5)
        call grph_insertEdge(rgraph, i3, i5)
      end if

      ! Delete broken edge (I1,I3) and add new edges (I1,I6),(I3,I6)
      if (e3 .eq. e6) then
        call grph_removeEdge(rgraph, i1, i3)
        call grph_insertEdge(rgraph, i3, i6)
        call grph_insertEdge(rgraph, i1, i6)
      end if

      ! Delete broken edge (I3,I4) and add new edges (I5,I6),(I4,I6),(I4,I5)
      call grph_removeEdge(rgraph, i3, i4)
      call grph_insertEdge(rgraph, i4, i5)
      call grph_insertEdge(rgraph, i5, i6)
      call grph_insertEdge(rgraph, i4, i6)


    case(HADAPT_OPR_CVT_QUAD2QUAD)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)
      i7 = rcollection%IquickAccess(7)
      i8 = rcollection%IquickAccess(8)
      i9 = rcollection%IquickAccess(9)

      e1 = rcollection%IquickAccess(10)
      e2 = rcollection%IquickAccess(11)
      e3 = rcollection%IquickAccess(12)
      e4 = rcollection%IquickAccess(13)
      e5 = rcollection%IquickAccess(14)
      e6 = rcollection%IquickAccess(15)
      e7 = rcollection%IquickAccess(16)
      e8 = rcollection%IquickAccess(17)

      ! Delete broken edge (I2,I3) and add new edges (I2,I6),(I3,I6)
      if (e2 .eq. e6) then
        call grph_removeEdge(rgraph, i2, i3)
        call grph_insertEdge(rgraph, i2, i6)
        call grph_insertEdge(rgraph, i3, i6)
      end if

      ! Delete broken edge (I1,I4) and add new edges  (I1,I8),(I4,I8)
      if (e4 .eq. e8) then
        call grph_removeEdge(rgraph, i1, i4)
        call grph_insertEdge(rgraph, i1, i8)
        call grph_insertEdge(rgraph, i4, i8)
      end if

      ! Delete broken edges (I5,I7),(I2,I7),(I3,I5),(I1,I7),(I4,I5) and
      ! add new edges (I5,I9),(I6,I9),(I7,I9),(I8,I9),(I1,I9),(I2,I9),
      ! (I3,I9),(I4,I9),(I5,I8),(I5,I6),(I6,I7),(I7,I8)
      call grph_removeEdge(rgraph, i5, i7)
      call grph_removeEdge(rgraph, i2, i7)
      call grph_removeEdge(rgraph, i3, i5)
      call grph_removeEdge(rgraph, i1, i7)
      call grph_removeEdge(rgraph, i4, i5)
      call grph_insertEdge(rgraph, i5, i9)
      call grph_insertEdge(rgraph, i6, i9)
      call grph_insertEdge(rgraph, i7, i9)
      call grph_insertEdge(rgraph, i8, i9)
      call grph_insertEdge(rgraph, i1, i9)
      call grph_insertEdge(rgraph, i2, i9)
      call grph_insertEdge(rgraph, i3, i9)
      call grph_insertEdge(rgraph, i4, i9)
      call grph_insertEdge(rgraph, i5, i8)
      call grph_insertEdge(rgraph, i5, i6)
      call grph_insertEdge(rgraph, i6, i7)
      call grph_insertEdge(rgraph, i7, i8)


    case(HADAPT_OPR_CVT_QUAD3TRIA)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)
      i7 = rcollection%IquickAccess(7)
      i8 = rcollection%IquickAccess(8)
      i9 = rcollection%IquickAccess(9)

      e1 = rcollection%IquickAccess(10)
      e2 = rcollection%IquickAccess(11)
      e3 = rcollection%IquickAccess(12)
      e4 = rcollection%IquickAccess(13)
      e5 = rcollection%IquickAccess(14)
      e6 = rcollection%IquickAccess(15)
      e7 = rcollection%IquickAccess(16)
      e8 = rcollection%IquickAccess(17)

      ! Delete broken edge (I2,I3) and add new edges (I2,I6),(I3,I6)
      if (e2 .eq. e6) then
        call grph_removeEdge(rgraph, i2, i3)
        call grph_insertEdge(rgraph, i2, i6)
        call grph_insertEdge(rgraph, i3, i6)
      end if

      ! Delete broken edge (I3,I4) and add new edges (I3,I7),(I4,I7)
      if (e3 .eq. e7) then
        call grph_removeEdge(rgraph, i3, i4)
        call grph_insertEdge(rgraph, i3, i7)
        call grph_insertEdge(rgraph, i4, i7)
      end if

      ! Delete broken edge (I1,I4) and add new edges (I4,I8),(I1,I8)
      if (e4 .eq. e8) then
        call grph_removeEdge(rgraph, i1, i4)
        call grph_insertEdge(rgraph, i4, i8)
        call grph_insertEdge(rgraph, i1, i8)
      end if

      ! Delete broken edges (I5,I3) and (I5,I4) and add new edges
      ! (I5,I9),(I6,I9),(I7,I9),(I8,I9),(I1,I9),(I5,I8),(I2,I9),(I5,I6)
      ! (I3,I9),I6,I7),(I4,I9), and (I7,I8)
      call grph_removeEdge(rgraph, i3, i5)
      call grph_removeEdge(rgraph, i4, i5)

      call grph_insertEdge(rgraph, i1, i9)
      call grph_insertEdge(rgraph, i2, i9)
      call grph_insertEdge(rgraph, i3, i9)
      call grph_insertEdge(rgraph, i4, i9)
      call grph_insertEdge(rgraph, i5, i9)
      call grph_insertEdge(rgraph, i6, i9)
      call grph_insertEdge(rgraph, i7, i9)
      call grph_insertEdge(rgraph, i8, i9)
      call grph_insertEdge(rgraph, i5, i6)
      call grph_insertEdge(rgraph, i6, i7)
      call grph_insertEdge(rgraph, i7, i8)
      call grph_insertEdge(rgraph, i5, i8)


    case(HADAPT_OPR_CVT_QUAD4TRIA)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)
      i7 = rcollection%IquickAccess(7)
      i8 = rcollection%IquickAccess(8)
      i9 = rcollection%IquickAccess(9)

      e1 = rcollection%IquickAccess(10)
      e2 = rcollection%IquickAccess(11)
      e3 = rcollection%IquickAccess(12)
      e4 = rcollection%IquickAccess(13)
      e5 = rcollection%IquickAccess(14)
      e6 = rcollection%IquickAccess(15)
      e7 = rcollection%IquickAccess(16)
      e8 = rcollection%IquickAccess(17)

      ! Delete broken edge (I3,I4) and add new edges (I3,I7),(I4,I7)
      if (e3 .eq. e7) then
        call grph_removeEdge(rgraph, i3, i4)
        call grph_insertEdge(rgraph, i3, i7)
        call grph_insertEdge(rgraph, i4, i7)
      end if

      ! Delete broken edge (I1,I4) and add new edges (I4,I8),(I1,I8)
      if (e4 .eq. e8) then
        call grph_removeEdge(rgraph, i1, i4)
        call grph_insertEdge(rgraph, i4, i8)
        call grph_insertEdge(rgraph, i1, i8)
      end if

      ! Delete broken edges (I4,I5),(I5,I6) and (I4,I6) and add new
      ! edges (I5,I9),(I6,I9),(I7,I9),(I8,I9),(I1,I9),(I5,I8),(I2,I9),
      ! (I5,I6),(I3,I9),I6,I7),(I4,I9), and (I7,I8)
      call grph_removeEdge(rgraph, i4, i5)
      call grph_removeEdge(rgraph, i5, i6)
      call grph_removeEdge(rgraph, i4, i6)
      call grph_insertEdge(rgraph, i5, i9)
      call grph_insertEdge(rgraph, i6, i9)
      call grph_insertEdge(rgraph, i7, i9)
      call grph_insertEdge(rgraph, i8, i9)
      call grph_insertEdge(rgraph, i1, i9)
      call grph_insertEdge(rgraph, i5, i8)
      call grph_insertEdge(rgraph, i2, i9)
      call grph_insertEdge(rgraph, i5, i6)
      call grph_insertEdge(rgraph, i3, i9)
      call grph_insertEdge(rgraph, i6, i7)
      call grph_insertEdge(rgraph, i4, i9)
      call grph_insertEdge(rgraph, i7, i8)


    case(HADAPT_OPR_CRS_2TRIA1TRIA)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)

      e1 = rcollection%IquickAccess(5)
      e2 = rcollection%IquickAccess(6)
      e3 = rcollection%IquickAccess(7)
      e4 = rcollection%IquickAccess(8)

      ! Delete broken edges (I1,I4), and (I2,I4) and add new edge (I1,I2)
      if (e1 .eq. e4) then
        call grph_removeEdge(rgraph, i1, i4)
        call grph_removeEdge(rgraph, i2, i4)
        call grph_insertEdge(rgraph, i1, i2)
      end if

      ! Delete broken edge (I3,I4)
      call grph_removeEdge(rgraph, i3, i4)


    case(HADAPT_OPR_CRS_4TRIA1TRIA)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)

      e1 = rcollection%IquickAccess(7)
      e2 = rcollection%IquickAccess(8)
      e3 = rcollection%IquickAccess(9)
      e4 = rcollection%IquickAccess(10)
      e5 = rcollection%IquickAccess(11)
      e6 = rcollection%IquickAccess(12)

      ! Delete broken edges (I4,I5),(I4,I6), and (I5,I6)
      call grph_removeEdge(rgraph, i4, i5)
      call grph_removeEdge(rgraph, i5, i6)
      call grph_removeEdge(rgraph, i4, i6)

      ! Delete broken edges (I1,I4), and (I2,I4) and add new edge (I1,I2)
      if (e1 .eq. e4) then
        call grph_removeEdge(rgraph, i1, i4)
        call grph_removeEdge(rgraph, i2, i4)
        call grph_insertEdge(rgraph, i1, i2)
      end if

      ! Delete broken edges (I2,I5), and (I3,I5) and add new edge (I2,I3)
      if (e2 .eq. e5) then
        call grph_removeEdge(rgraph, i2, i5)
        call grph_removeEdge(rgraph, i3, i5)
        call grph_insertEdge(rgraph, i2, i3)
      end if

      ! Delete broken edges (I1,I6), and (I3,I6) and add new edge (I1,I3)
      if (e3 .eq. e6) then
        call grph_removeEdge(rgraph, i1, i6)
        call grph_removeEdge(rgraph, i3, i6)
        call grph_insertEdge(rgraph, i1, i3)
      end if


    case(HADAPT_OPR_CRS_4TRIA2TRIA1)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)

      e1 = rcollection%IquickAccess(7)
      e2 = rcollection%IquickAccess(8)
      e3 = rcollection%IquickAccess(9)
      e4 = rcollection%IquickAccess(10)
      e5 = rcollection%IquickAccess(11)
      e6 = rcollection%IquickAccess(12)

      ! Delete broken edges (I4,I5),(I5,I6), and (I4,I6) and add new edge (I3,I4)
      call grph_removeEdge(rgraph, i4, i5)
      call grph_removeEdge(rgraph, i5, i6)
      call grph_removeEdge(rgraph, i4, i6)
      call grph_insertEdge(rgraph, i3, i4)

      ! Delete broken edges (I2,I5), and (I3,I5) and add new edge (I2,I3)
      if (e2 .eq. e5) then
        call grph_removeEdge(rgraph, i2, i5)
        call grph_removeEdge(rgraph, i3, i5)
        call grph_insertEdge(rgraph, i2, i3)
      end if

      ! Delete broken edges (I1,I6), and (I3,I6) and add new edge (I1,I3)
      if (e3 .eq. e6) then
        call grph_removeEdge(rgraph, i1, i6)
        call grph_removeEdge(rgraph, i3, i6)
        call grph_insertEdge(rgraph, i1, i3)
      end if


    case(HADAPT_OPR_CRS_4TRIA2TRIA2)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)

      e1 = rcollection%IquickAccess(7)
      e2 = rcollection%IquickAccess(8)
      e3 = rcollection%IquickAccess(9)
      e4 = rcollection%IquickAccess(10)
      e5 = rcollection%IquickAccess(11)
      e6 = rcollection%IquickAccess(12)

      ! Delete broken edges (I4,I5),(I5,I6), and (I4,I6) and add new edge (I1,I5)
      call grph_removeEdge(rgraph, i4, i5)
      call grph_removeEdge(rgraph, i5, i6)
      call grph_removeEdge(rgraph, i4, i6)
      call grph_insertEdge(rgraph, i1, i5)

      ! Delete broken edges (I1,I4), and (I2,I4) and add new edge (I1,I4)
      if (e1 .eq. e4) then
        call grph_removeEdge(rgraph, i1, i4)
        call grph_removeEdge(rgraph, i2, i4)
        call grph_insertEdge(rgraph, i1, i2)
      end if

      ! Delete broken edges (I1,I6), and (I3,I6) and add new edge (I1,I3)
      if (e3 .eq. e6) then
        call grph_removeEdge(rgraph, i1, i6)
        call grph_removeEdge(rgraph, i3, i6)
        call grph_insertEdge(rgraph, i1, i3)
      end if


    case(HADAPT_OPR_CRS_4TRIA2TRIA3)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)

      e1 = rcollection%IquickAccess(7)
      e2 = rcollection%IquickAccess(8)
      e3 = rcollection%IquickAccess(9)
      e4 = rcollection%IquickAccess(10)
      e5 = rcollection%IquickAccess(11)
      e6 = rcollection%IquickAccess(12)

      ! Delete broken edges (I4,I5),(I5,I6), and (I4,I6) and add new edge (I2,I6)
      call grph_removeEdge(rgraph, i4, i5)
      call grph_removeEdge(rgraph, i5, i6)
      call grph_removeEdge(rgraph, i4, i6)
      call grph_insertEdge(rgraph, i2, i6)

      ! Delete broken edges (I1,I4), and (I2,I4) and add new edge (I1,I4)
      if (e1 .eq. e4) then
        call grph_removeEdge(rgraph, i1, i4)
        call grph_removeEdge(rgraph, i2, i4)
        call grph_insertEdge(rgraph, i1, i2)
      end if

      ! Delete broken edges (I2,I5), and (I3,I5) and add new edge (I2,I3)
      if (e2 .eq. e5) then
        call grph_removeEdge(rgraph, i2, i5)
        call grph_removeEdge(rgraph, i3, i5)
        call grph_insertEdge(rgraph, i2, i3)
      end if


    case(HADAPT_OPR_CRS_4QUAD1QUAD)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)
      i7 = rcollection%IquickAccess(7)
      i8 = rcollection%IquickAccess(8)
      i9 = rcollection%IquickAccess(9)

      e1 = rcollection%IquickAccess(10)
      e2 = rcollection%IquickAccess(11)
      e3 = rcollection%IquickAccess(12)
      e4 = rcollection%IquickAccess(13)
      e5 = rcollection%IquickAccess(14)
      e6 = rcollection%IquickAccess(15)
      e7 = rcollection%IquickAccess(16)
      e8 = rcollection%IquickAccess(17)

      ! Delete broken edges (I1,I9),(I2,I9),(I3,I9),(I4,I9),(I5,I9),(I6,I9),
      ! (I7,I9),(I8,I9),(I5,I6),(I6,I7),(I7,I8), and (I5,I8)
      call grph_removeEdge(rgraph, i1, i9)
      call grph_removeEdge(rgraph, i2, i9)
      call grph_removeEdge(rgraph, i3, i9)
      call grph_removeEdge(rgraph, i4, i9)
      call grph_removeEdge(rgraph, i5, i9)
      call grph_removeEdge(rgraph, i6, i9)
      call grph_removeEdge(rgraph, i7, i9)
      call grph_removeEdge(rgraph, i8, i9)
      call grph_removeEdge(rgraph, i5, i6)
      call grph_removeEdge(rgraph, i6, i7)
      call grph_removeEdge(rgraph, i7, i8)
      call grph_removeEdge(rgraph, i8, i5)

      ! Delete broken edges (I1,I5) and (I2,I5) and add new edge (I1,I2)
      if (e1 .eq. e5) then
        call grph_removeEdge(rgraph, i1, i5)
        call grph_removeEdge(rgraph, i2, i5)
        call grph_insertEdge(rgraph, i1, i2)
      end if

      ! Delete broken edges (I2,I6) and (I3,I6) and add new edge (I2,I3)
      if (e2 .eq. e6) then
        call grph_removeEdge(rgraph, i2, i6)
        call grph_removeEdge(rgraph, i3, i6)
        call grph_insertEdge(rgraph, i2, i3)
      end if

      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (e3 .eq. e7) then
        call grph_removeEdge(rgraph, i3, i7)
        call grph_removeEdge(rgraph, i4, i7)
        call grph_insertEdge(rgraph, i3, i4)
      end if

      ! Delete broken edges (I4,I8) and (I1,I8) and add new edge (I4,I1)
      if (e4 .eq. e8) then
        call grph_removeEdge(rgraph, i4, i8)
        call grph_removeEdge(rgraph, i1, i8)
        call grph_insertEdge(rgraph, i1, i4)
      end if

      ! Add new edges (I1,I3) and (I2,I4)
      call grph_insertEdge(rgraph, i1, i3)
      call grph_insertEdge(rgraph, i2, i4)


    case(HADAPT_OPR_CRS_4QUAD2QUAD)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)
      i7 = rcollection%IquickAccess(7)
      i8 = rcollection%IquickAccess(8)
      i9 = rcollection%IquickAccess(9)

      e1 = rcollection%IquickAccess(10)
      e2 = rcollection%IquickAccess(11)
      e3 = rcollection%IquickAccess(12)
      e4 = rcollection%IquickAccess(13)
      e5 = rcollection%IquickAccess(14)
      e6 = rcollection%IquickAccess(15)
      e7 = rcollection%IquickAccess(16)
      e8 = rcollection%IquickAccess(17)

      ! Delete broken edges (I1,I9),(I2,I9),(I3,I9),(I4,I9),(I5,I9),(I6,I9),
      ! (I7,I9),(I8,I9),(I5,I6),(I6,I7),(I7,I8), and (I5,I8)
      call grph_removeEdge(rgraph, i1, i9)
      call grph_removeEdge(rgraph, i2, i9)
      call grph_removeEdge(rgraph, i3, i9)
      call grph_removeEdge(rgraph, i4, i9)
      call grph_removeEdge(rgraph, i5, i9)
      call grph_removeEdge(rgraph, i6, i9)
      call grph_removeEdge(rgraph, i7, i9)
      call grph_removeEdge(rgraph, i8, i9)
      call grph_removeEdge(rgraph, i5, i6)
      call grph_removeEdge(rgraph, i6, i7)
      call grph_removeEdge(rgraph, i7, i8)
      call grph_removeEdge(rgraph, i8, i5)

      ! Delete broken edges (I2,I6) and (I3,I6) and add new edge (I2,I3)
      if (e2 .eq. e6) then
        call grph_removeEdge(rgraph, i2, i6)
        call grph_removeEdge(rgraph, i3, i6)
        call grph_insertEdge(rgraph, i2, i3)
      end if

      ! Delete broken edges (I4,I8) and (I1,I8) and add new edge (I4,I1)
      if (e4 .eq. e8) then
        call grph_removeEdge(rgraph, i4, i8)
        call grph_removeEdge(rgraph, i1, i8)
        call grph_insertEdge(rgraph, i1, i4)
      end if

      ! Add new edges (I5,I7),(I1,I7),(I4,I5),(I2,I7) and (I3,I5)
      call grph_insertEdge(rgraph, i5, i7)
      call grph_insertEdge(rgraph, i1, i7)
      call grph_insertEdge(rgraph, i4, i5)
      call grph_insertEdge(rgraph, i2, i7)
      call grph_insertEdge(rgraph, i3, i5)


    case(HADAPT_OPR_CRS_4QUAD3TRIA)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)
      i7 = rcollection%IquickAccess(7)
      i8 = rcollection%IquickAccess(8)
      i9 = rcollection%IquickAccess(9)

      e1 = rcollection%IquickAccess(10)
      e2 = rcollection%IquickAccess(11)
      e3 = rcollection%IquickAccess(12)
      e4 = rcollection%IquickAccess(13)
      e5 = rcollection%IquickAccess(14)
      e6 = rcollection%IquickAccess(15)
      e7 = rcollection%IquickAccess(16)
      e8 = rcollection%IquickAccess(17)

      ! Delete broken edges (I1,I9),(I2,I9),(I3,I9),(I4,I9),(I5,I9),(I6,I9),
      ! (I7,I9),(I8,I9),(I5,I6),(I6,I7),(I7,I8), and (I5,I8)
      call grph_removeEdge(rgraph, i1, i9)
      call grph_removeEdge(rgraph, i2, i9)
      call grph_removeEdge(rgraph, i3, i9)
      call grph_removeEdge(rgraph, i4, i9)
      call grph_removeEdge(rgraph, i5, i9)
      call grph_removeEdge(rgraph, i6, i9)
      call grph_removeEdge(rgraph, i7, i9)
      call grph_removeEdge(rgraph, i8, i9)
      call grph_removeEdge(rgraph, i5, i6)
      call grph_removeEdge(rgraph, i6, i7)
      call grph_removeEdge(rgraph, i7, i8)
      call grph_removeEdge(rgraph, i8, i5)

      ! Delete broken edges (I2,I6) and (I3,I6) and add new edge (I2,I3)
      if (e2 .eq. e6) then
        call grph_removeEdge(rgraph, i2, i6)
        call grph_removeEdge(rgraph, i3, i6)
        call grph_insertEdge(rgraph, i2, i3)
      end if

      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (e3 .eq. e7) then
        call grph_removeEdge(rgraph, i3, i7)
        call grph_removeEdge(rgraph, i4, i7)
        call grph_insertEdge(rgraph, i3, i4)
      end if

      ! Delete broken edges (I4,I8) and (I1,I8) and add new edge (I4,I1)
      if (e4 .eq. e8) then
        call grph_removeEdge(rgraph, i4, i8)
        call grph_removeEdge(rgraph, i1, i8)
        call grph_insertEdge(rgraph, i1, i4)
      end if

      ! Add new edges (I3,I5) and (I4,I5)
      call grph_insertEdge(rgraph, i3, i5)
      call grph_insertEdge(rgraph, i4, i5)


    case(HADAPT_OPR_CRS_4QUAD4TRIA)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)
      i7 = rcollection%IquickAccess(7)
      i8 = rcollection%IquickAccess(8)
      i9 = rcollection%IquickAccess(9)

      e1 = rcollection%IquickAccess(10)
      e2 = rcollection%IquickAccess(11)
      e3 = rcollection%IquickAccess(12)
      e4 = rcollection%IquickAccess(13)
      e5 = rcollection%IquickAccess(14)
      e6 = rcollection%IquickAccess(15)
      e7 = rcollection%IquickAccess(16)
      e8 = rcollection%IquickAccess(17)

      ! Delete broken edges (I1,I9),(I2,I9),(I3,I9),(I4,I9),(I5,I9),(I6,I9),
      ! (I7,I9),(I8,I9),(I5,I6),(I6,I7), and (I7,I8)
      call grph_removeEdge(rgraph, i1, i9)
      call grph_removeEdge(rgraph, i2, i9)
      call grph_removeEdge(rgraph, i3, i9)
      call grph_removeEdge(rgraph, i4, i9)
      call grph_removeEdge(rgraph, i5, i9)
      call grph_removeEdge(rgraph, i6, i9)
      call grph_removeEdge(rgraph, i7, i9)
      call grph_removeEdge(rgraph, i8, i9)
      call grph_removeEdge(rgraph, i5, i6)
      call grph_removeEdge(rgraph, i6, i7)
      call grph_removeEdge(rgraph, i7, i8)

      ! Delete broken edges (I2,I6) and (I3,I6) and add new edge (I2,I3)
      if (e2 .eq. e6) then
        call grph_removeEdge(rgraph, i2, i6)
        call grph_removeEdge(rgraph, i3, i6)
        call grph_insertEdge(rgraph, i2, i3)
      end if

      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (e3 .eq. e7) then
        call grph_removeEdge(rgraph, i3, i7)
        call grph_removeEdge(rgraph, i4, i7)
        call grph_insertEdge(rgraph, i3, i4)
      end if

      ! Add new edges (I3,I5) and (I3,I8)
      call grph_insertEdge(rgraph, i3, i5)
      call grph_insertEdge(rgraph, i3, i8)


    case(HADAPT_OPR_CRS_2QUAD1QUAD)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)
      i7 = rcollection%IquickAccess(7)
      i8 = rcollection%IquickAccess(8)

      e1 = rcollection%IquickAccess(9)
      e2 = rcollection%IquickAccess(10)
      e3 = rcollection%IquickAccess(11)
      e4 = rcollection%IquickAccess(12)
      e5 = rcollection%IquickAccess(13)
      e6 = rcollection%IquickAccess(14)
      e7 = rcollection%IquickAccess(15)
      e8 = rcollection%IquickAccess(16)

      ! Delete broken edges (I5,I7), (I1,I7),(I4,I5),(I2,I7) and (I3,I5)
      call grph_removeEdge(rgraph, i5, i7)
      call grph_removeEdge(rgraph, i1, i7)
      call grph_removeEdge(rgraph, i2, i7)
      call grph_removeEdge(rgraph, i3, i5)
      call grph_removeEdge(rgraph, i4, i5)

      ! Delete broken edges (I1,I5) and (I2,I5) and add new edge (I1,I2)
      if (e1 .eq. e5) then
        call grph_removeEdge(rgraph, i1, i5)
        call grph_removeEdge(rgraph, i2, i5)
        call grph_insertEdge(rgraph, i1, i2)
      end if

      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (e3 .eq. e7) then
        call grph_removeEdge(rgraph, i3, i7)
        call grph_removeEdge(rgraph, i4, i7)
        call grph_insertEdge(rgraph, i3, i4)
      end if

      ! Add new edges (I1,I3) and (I2,I4)
      call grph_insertEdge(rgraph, i1, i3)
      call grph_insertEdge(rgraph, i2, i4)


    case(HADAPT_OPR_CRS_2QUAD3TRIA)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)
      i7 = rcollection%IquickAccess(7)
      i8 = rcollection%IquickAccess(8)

      e1 = rcollection%IquickAccess(9)
      e2 = rcollection%IquickAccess(10)
      e3 = rcollection%IquickAccess(11)
      e4 = rcollection%IquickAccess(12)
      e5 = rcollection%IquickAccess(13)
      e6 = rcollection%IquickAccess(14)
      e7 = rcollection%IquickAccess(15)
      e8 = rcollection%IquickAccess(16)

      ! Delete broken edges (I5,I7),(I1,I7),(I4,I5),(I2,(I7) and (I3,I5)
      call grph_removeEdge(rgraph, i5, i7)
      call grph_removeEdge(rgraph, i1, i7)
      call grph_removeEdge(rgraph, i2, i7)
      call grph_removeEdge(rgraph, i3, i5)
      call grph_removeEdge(rgraph, i4, i5)

      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (e3 .eq. e7) then
        call grph_removeEdge(rgraph, i3, i7)
        call grph_removeEdge(rgraph, i4, i7)
        call grph_insertEdge(rgraph, i3, i4)
      end if

      ! Add new edges (I3,I5) and (I4,I5)
      call grph_insertEdge(rgraph, i3, i5)
      call grph_insertEdge(rgraph, i4, i5)


    case(HADAPT_OPR_CRS_3TRIA1QUAD)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)

      e1 = rcollection%IquickAccess(6)
      e2 = rcollection%IquickAccess(7)
      e3 = rcollection%IquickAccess(8)
      e4 = rcollection%IquickAccess(9)
      e5 = rcollection%IquickAccess(10)
      e6 = rcollection%IquickAccess(11)
      e7 = rcollection%IquickAccess(12)
      e8 = rcollection%IquickAccess(13)

      ! Delete broken edges (I3,I5) and (I4,I5)
      call grph_removeEdge(rgraph, i3, i5)
      call grph_removeEdge(rgraph, i4, i5)

      ! Delete broken edges (I1,I5) and (I2,I5) and add new edge (I1,I2)
      if (e1 .eq. e5) then
        call grph_removeEdge(rgraph, i1, i5)
        call grph_removeEdge(rgraph, i2, i5)
        call grph_insertEdge(rgraph, i1, i2)
      end if

      ! Add new edges (I1,I3) and (I2,I4)
      call grph_insertEdge(rgraph, i1, i3)
      call grph_insertEdge(rgraph, i2, i4)


    case(HADAPT_OPR_CRS_4TRIA1QUAD)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)
      i7 = rcollection%IquickAccess(7)
      i8 = rcollection%IquickAccess(8)

      e1 = rcollection%IquickAccess(9)
      e2 = rcollection%IquickAccess(10)
      e3 = rcollection%IquickAccess(11)
      e4 = rcollection%IquickAccess(12)
      e5 = rcollection%IquickAccess(13)
      e6 = rcollection%IquickAccess(14)
      e7 = rcollection%IquickAccess(15)
      e8 = rcollection%IquickAccess(16)

      ! Delete broken edges (I1,I6),(I1,I7) and (I6,I7)
      call grph_removeEdge(rgraph, i1, i6)
      call grph_removeEdge(rgraph, i1, i7)
      call grph_removeEdge(rgraph, i6, i7)

      ! Delete broken edges (I2,I6) and (I3,I6) and add new edge (I2,I3)
      if (e2 .eq. e6) then
        call grph_removeEdge(rgraph, i2, i6)
        call grph_removeEdge(rgraph, i3, i6)
        call grph_insertEdge(rgraph, i2, i3)
      end if

      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (e3 .eq. e7) then
        call grph_removeEdge(rgraph, i3, i7)
        call grph_removeEdge(rgraph, i4, i7)
        call grph_insertEdge(rgraph, i3, i4)
      end if

      ! Add new edges (I1,I3) and (I2,I4)
      call grph_insertEdge(rgraph, i1, i3)
      call grph_insertEdge(rgraph, i2, i4)


    case(HADAPT_OPR_CRS_4TRIA3TRIA2)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)
      i7 = rcollection%IquickAccess(7)
      i8 = rcollection%IquickAccess(8)

      e1 = rcollection%IquickAccess(9)
      e2 = rcollection%IquickAccess(10)
      e3 = rcollection%IquickAccess(11)
      e4 = rcollection%IquickAccess(12)
      e5 = rcollection%IquickAccess(13)
      e6 = rcollection%IquickAccess(14)
      e7 = rcollection%IquickAccess(15)
      e8 = rcollection%IquickAccess(16)

      ! Delete broken edges (I1,I7) and (I6,I7)
      call grph_removeEdge(rgraph, i1, i7)
      call grph_removeEdge(rgraph, i6, i7)

      ! Delete broken edges (I3,I7) and (I4,I7) and add new edge (I3,I4)
      if (e3 .eq. e7) then
        call grph_removeEdge(rgraph, i3, i7)
        call grph_removeEdge(rgraph, i4, i7)
        call grph_insertEdge(rgraph, i3, i4)
      end if

      ! Add new edge (I4,I6)
      call grph_insertEdge(rgraph, i4, i6)


    case(HADAPT_OPR_CRS_4TRIA3TRIA3)
      i1 = rcollection%IquickAccess(1)
      i2 = rcollection%IquickAccess(2)
      i3 = rcollection%IquickAccess(3)
      i4 = rcollection%IquickAccess(4)
      i5 = rcollection%IquickAccess(5)
      i6 = rcollection%IquickAccess(6)
      i7 = rcollection%IquickAccess(7)
      i8 = rcollection%IquickAccess(8)

      e1 = rcollection%IquickAccess(9)
      e2 = rcollection%IquickAccess(10)
      e3 = rcollection%IquickAccess(11)
      e4 = rcollection%IquickAccess(12)
      e5 = rcollection%IquickAccess(13)
      e6 = rcollection%IquickAccess(14)
      e7 = rcollection%IquickAccess(15)
      e8 = rcollection%IquickAccess(16)

      ! Delete broken edges (I1,I6) and (I6,I7)
      call grph_removeEdge(rgraph, i1, i6)
      call grph_removeEdge(rgraph, i6, i7)

      ! Delete broken edges (I2,I6) and (I3,I6) and add new edge (I2,I3)
      if (e2 .eq. e6) then
        call grph_removeEdge(rgraph, i2, i6)
        call grph_removeEdge(rgraph, i3, i6)
        call grph_insertEdge(rgraph, i2, i3)
      end if

      ! Add new edge (I2,I7)
      call grph_insertEdge(rgraph, i2, i7)


    case default
      call output_line('Unsupported operation!',&
                       OU_CLASS_ERROR, OU_MODE_STD, 'flagship_hadaptCallback2D')
      call sys_halt()
    end select
  end subroutine flagship_hadaptCallback2D

  !*****************************************************************************

!<subroutine>

  subroutine flagship_hadaptCallback3D(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 3D.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(in) :: iOperation
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_graph), pointer, save :: rgraph
    character(len=SYS_STRLEN) :: ssectionName

    ! What operation should be performed?
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! This subroutine assumes that the name of the sparsity pattern
      ! is stored in the first quick access string.

      ! Get section name
      call collct_getvalue_string(rcollection,&
          'ssectionname', ssectionName)

      ! Retrieve sparsity pattern from collection and build sparsity-graph.
      rgraph => collct_getvalue_graph(rcollection,&
          trim(rcollection%SquickAccess(1)), ssectionName=ssectionName)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify sparsity graph
      nullify(rgraph)


    case(HADAPT_OPR_INSERTVERTEXEDGE,&
         HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into sparsity graph
      call grph_insertVertex(rgraph, rcollection%IquickAccess(1))


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from sparsity graph and solution
      if (rcollection%IquickAccess(2) .ne. 0) then
        call grph_removeVertex(rgraph, rcollection%IquickAccess(1),&
                               rcollection%IquickAccess(2))
      else
        call grph_removeVertex(rgraph, rcollection%IquickAccess(1))
      end if


    case default
      call output_line('Unsupported operation!',&
                       OU_CLASS_ERROR, OU_MODE_STD, 'flagship_hadaptCallback3D')
      call sys_halt()
    end select
  end subroutine flagship_hadaptCallback3D

end module flagship_callback
