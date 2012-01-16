!##############################################################################
!# ****************************************************************************
!# <name> hadaptaux2d </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# WARNING: Do not USE this module in your applications unless you really
!#          know what you are doing. This module does no error checking!!!
!#
!# This module contains all auxiliary routines which are required for
!# performing h-adaptivity in 2D.
!#
!# The following routines are available:
!#
!#  1.) hadapt_markRefinement2D
!#      -> Marks elements for refinement in 2D
!#
!#  2.) hadapt_markRedgreenCoarsening2D
!#      -> Marks elements for coarsening in 2D
!#
!#  3.) hadapt_markRedgreenRefinement2D
!#      -> Marks elements for refinement due to Red-Green strategy in 2D
!#
!#  4.) hadapt_calcNumberOfElements2D
!#      -> Calculates the number of elements after refinement
!#
!#  5.) hadapt_refine2D
!#      -> Performs refinement of marked elements in 2D
!#
!#  6.) hadapt_coarsen2D
!#      -> Performs re-coarsening of marked elements in 2D
!#
!#  7.) hadapt_writeGridSVG2D
!#      -> Writes the adapted grid to file in SVG format in 2D
!#
!#
!# The following auxiliary routines are available:
!#
!#  1.) add_vertex2D = add_vertex_atEdgeMidpoint2D /
!#                     add_vertex_atElementCenter2D
!#      -> Adds a new vertex to the adaptation data structure in 2D
!#
!#  2.) remove_vertex2D
!#      -> Removes an existing vertex from the adaptation data structure in 2D
!#
!#  3.) replace_element2D = replace_elementTria /
!#                          replace_elementQuad
!#      -> Replaces an existing element by another element of he same type in 2D
!#
!#  4.) add_element2D = add_elementTria /
!#                      add_elementQuad
!#      -> Adds a new element to the adaptation data structure in 2D
!#
!#  5.) remove_element2D
!#      -> Removes an existing element from the adaptation data structure in 2D
!#
!#  6.) update_ElementNeighbors2D = update_ElemNeighb2D_1to2 /
!#                                  update_ElemNeighb2D_2to2
!#      -> Updates the list of neighboring elements in 2D
!#
!#  7.) update_AllElementNeighbors2D
!#      -> Updates the lists of neighboring elements of ALL adjacent elements in 2D
!#
!#  8.) refine_Tria2Tria
!#      -> Refines a triangle by subdivision into two triangles
!#
!#  9.) refine_Tria3Tria
!#      -> Refines a triangle by subdivision into three triangles
!#
!# 10.) refine_Tria4Tria
!#      -> Refines a triangle by subdivision into four triangles
!#
!# 11.) refine_Quad2Quad
!#      -> Refines a quadrilateral by subdivision into two quadrilaterals
!#
!# 12.) refine_Quad3Tria
!#      -> Refines a quadrilateral by subdivision into three triangles
!#
!# 13.) refine_Quad4Tria
!#      -> Refines a quadrilateral by subdivision into four triangles
!#
!# 14.) refine_Quad4Quad
!#      -> Refines a quadrilateral by subdivision into four quadrilaterals
!#
!# 15.) convert_Tria2Tria
!#      -> Converts two neighboring triangles into four similar triangle
!#
!# 16.) convert_Quad2Quad
!#      -> Converts two neighboring quadrilaterals into four similar quadrilaterals
!#
!# 17.) convert_Quad3Tria
!#      -> Converts three neighboring triangles into four similar quadrilaterals
!#
!# 18.) convert_Quad4Tria
!#      -> Converts four neighboring triangles into four similar quadrilaterals
!#
!# 19.) coarsen_2Tria1Tria
!#      -> Coarsens two green triangles by combining them to the macro element
!#
!# 20.) coarsen_4Tria1Tria
!#      -> Coarsens four red triangles by combining them to the macro element
!#
!# 21.) coarsen_4Tria2Tria
!#      -> Coarsens four red triangles by combining them to two green elements
!#
!# 22.) coarsen_4Quad1Quad
!#      -> Coarsens four red quadrilaterals by combining them to the macro element
!#
!# 23.) coarsen_4Quad2Quad
!#      -> Coarsens four red quadrilaterals by combining them to two green elements
!#
!# 24.) coarsen_4Quad3Tria
!#      -> Coarsens four red quadrilaterals by combining them to three green elements
!#
!# 25.) coarsen_4Quad4Tria
!#      -> Coarsens four red quadrilaterals by combining them to four green elements
!#
!# 26.) coarsen_2Quad1Quad
!#      -> Coarsens two green quadrilaterals by combining them to the macro element
!#
!# 27.) coarsen_2Quad3Tria
!#      -> Coarsens two green quadrilaterals by combining them to three green triangles
!#
!# 28.) coarsen_3Tria1Quad
!#      -> Coarsens three green triangles by combining them to the macro element
!#
!# 29.) coarsen_4Tria1Quad
!#      -> Coarsens four green triangles by combining them to the macro element
!#
!# 30.) coarsen_4Tria3Tria
!#      -> Coarsens four green triangles by combining them to three green triangles
!#
!# 31. redgreen_getState2D
!#      -> Returns the state of an element in 2D
!#
!# 32.) redgreen_getStateTria
!#      -> Returns the state of a triangle
!#
!# 33.) redgreen_getStateQuad
!#      -> Returns the state of a quadrilateral
!#
!# 34.) redgreen_rotateState2D
!#      -> Computes the state of an element after rotation in 2D
!#
!#
!#  FAQ - Some explanation  \\
!# ------------------------ \\
!# 1.) So, how does red-green grid refinement work in detail?
!#
!#     In general, a conforming mesh in two spatial dimensions consists of
!#     triangles and/or quadrilaterals which have common vertices, that is,
!#     there is no vertex located somewhere in the edge of one element.
!#     If you want to refine your mesh, a natural choice is to subdivide
!#     (a) one triangle into four similar triangles by placing new vertices
!#         at all midpoints of the three edges and connecting these new nodes
!#     (b) one quadrilateral into four similar quadrilaterals by placing
!#         four new vertices at the edge midpoints and one vertex at the
!#         center and connect all new nodes to the centered vertex.
!#     These two strategies are known as red-refinement.
!#
!#     In order to remove the hanging nodes in adjacent elements you have
!#     to refine these elements, too. Here, all all possibilities:
!#     (c) Refine one triangle into two triangles by placing one new node
!#         at the midpoint of one edge.
!#     (d) Refine one quadrilateral into two quadrilaterals by inserting two
!#         new vertices at the midpoints of opposite edges and connect them.
!#     (e) Refine one quadrilateral into three triangles by inserting one new
!#         vertex at the midpoint of one edge and connect it to the two nodes
!#         of the opposite edge.
!#     (f) Refine one quadrilateral into four triangles by inserting two nodes
!#         at the midpoint of adjacent edges, connect these vertices with
!#         oneanother and with the opposite corner vertex.
!#     These four refinement strategies are known as green refinement.
!#
!#     Keep in mind, that there may be two hanging nodes for triangles or
!#     three hanging nodes for quadrilaterals. If these elements were refined
!#     without introducing new vertices, we would perform blue refinement.
!#     However, this strategy should be avoided and an element marked for
!#     blue refinement should be subdivided into four similar elements
!#     introducing one additional vertex at the edge midpoint.
!#
!#     Red refinement does not deteriorate the grid quality but green does.
!#     If one would proceed like this, then green-refined element would tend
!#     to produce very sharp angles leading to nasty "needle" elements.
!#     In order to refine a green element, the original element from which
!#     the green element was created has to be restored. Then red refinement
!#     is performed for this element.
!#
!#
!# 2.) What are the different numbers for the element marker MARK_xxx?
!#
!#     In two space dimensions, the maximum number of vertices per element and
!#     consequently the number of edges surrounding an element is two. Consider
!#     a 5-bitfield [bit4,bit3,bit2,bit1,bit0]. The bit0 is used to distinguish
!#     between triangles (bit0=0) and quadrilaterals (bit1=1). Each of bitX is
!#     associated with the X-th edge. If you want to mark an edge for refinement,
!#     then you have to set its bit to one. If you want to unmark it, then set
!#     bitX=0. Finally, the integer associated to this 5-bitfiled is the element
!#     marker. Here, is the table of bitfields and markers
!#
!#     <verb>
!#     bit4  bit3  bit2  bit1  bit0   integer   description
!#      0     0     0     0     0        0       triangle, not marked
!#      0     0     0     0     1        1       quadrilateral, not marked
!#      0     0     0     1     0        2       triangle, marked for 1-tria : 2-tria
!#                                               refinement along the first edge
!#      0     0     0     1     1        3       quadtrilateral, marked for 1-quad :
!#                                               3-tria refinement along the first edge
!#      0     0     1     0     0        4       triangle, marked for 1-tria : 2-tria
!#                                               refinement along the second edge
!#      0     0     1     0     1        5       quadtrilateral, marked for 1-quad :
!#                                               3-tria refinement along the second edge
!#      0     0     1     1     0        6       triangle, marked for 1-tria : 3-tria
!#                                               refinement along first/second edge
!#      0     0     1     1     1        7       quadrilateral, marked for 1-quad :
!#                                               4-tria refinement along first/second edge
!#      0     1     0     0     0        8       triangle, marked for 1-tria : 2-tria
!#                                               refinement along the third edge
!#      0     1     0     0     1        9       quadtrilateral, marked for 1-quad :
!#                                               3-tria refinement along the third edge
!#      0     1     0     1     0       10       triangle, marked for 1-tria : 3-tria
!#                                               refinement along first/third edge
!#      0     1     0     1     1       11       quadrilateral, makred for 1-quad :
!#                                               2-quad along first/third edge
!#      0     1     1     0     0       12       triangle, marked for 1-tria : 3-tria
!#                                               refinement along second/third edge
!#      0     1     1     0     1       13       quadrilateral, marked for 1-quad :
!#                                               4-tria refinement along second/third edge
!#      0     1     1     1     0       14       triangle marked for 1:4 red refinement
!#      0     1     1     1     1       15       quadrilateral, marked for blue refinement
!#      1     0     0     0     0       16       n.a.
!#      1     0     0     0     1       17       quadtrilateral, marked for 1-quad :
!#                                               3-tria refinement along the fourth edge
!#      1     0     0     1     0       18       n.a
!#      1     0     0     1     1       19       quadrilateral, marked for 1-quad :
!#                                               4-tria refinement along first/fourth edge
!#      1     0     1     0     0       20       n.a.
!#      1     0     1     0     1       21       quadrilateral, makred for 1-quad :
!#                                               2-quad along second/fourth edge
!#      1     0     1     1     0       22       n.a.
!#      1     0     1     1     1       23       quadrilateral, marked for blue refinement
!#      1     1     0     0     0       24       n.a.
!#      1     1     0     0     1       25       quadrilateral, marked for 1-quad :
!#                                               4-tria refinement along third/fourth edge
!#      1     1     0     1     0       26       n.a.
!#      1     1     0     1     1       27       quadrilateral, marked for blue refinement
!#      1     1     1     0     0       28       n.a.
!#      1     1     1     0     1       29       quadrilateral, marked for blue refinement
!#      1     1     1     1     0       30       n.a.
!#      1     1     1     1     1       31       quadrilateral marked for 1:4 red refinement
!#     </verb>
!#
!#
!# 3.) Ok, and what does the state of an element mean?
!#
!#     The state of an element is computed from the age of its vertices.
!#     In 2D, an element can have four vertices at most
!#
!#     Typically, a more or less complicated
!#     data structure is required to maintain the father-son relationship.
!#     In this implementaiton, we use a new strategy based on element states
!#     which allows for an efficient computation of father-son relations.
!#     For each vertex, we save its creation age. All vertices of the initial
!#     grid have age 0 and cannot be removed. If a new vertex is introduced
!#     its age is computed as 1+MAX("Age of surrounding vertices"), that is,
!#     1+MAX("Age of left neighbor","Age of right neighbour") for edge midpoints
!#     and 1+MAX("Age of four corner vertices") for centered node.
!#
!#
!# 4.) And how does grid re-coarsening work?
!#
!#     In short, each refinement operation is reversed by a uniquely define
!#     re-coarsening operation. In practice, grid re-caorsening is slightly
!#     more complicated since the difficult task is to identify groups of
!#     elements that can be combined (i.e. glued together) to the original
!#     macro elements. The identification of element patches is based on the
!#     age of vertices and on the element state.
!# </purpose>
!##############################################################################

module hadaptaux2d

  use fsystem
  use hadaptaux
  use linearsystemscalar
  use io
  use collection

  implicit none

  private
  public :: hadapt_markRefinement2D
  public :: hadapt_markRedgreenCoarsening2D
  public :: hadapt_markRedgreenRefinement2D
  public :: hadapt_calcNumberOfElements2D
  public :: hadapt_refine2D
  public :: hadapt_coarsen2D
  public :: hadapt_writeGridSVG2D

!<constants>

!<constantblock description="Constants for element marker in 2D">

  ! Mark sub-element for generic recoarsening. Note that the detailed
  ! recoarsening strategy must be determined form the element state.
  integer, parameter, public :: MARK_CRS_GENERIC            = -1

  ! Mark inner triangle of a 1-tria : 4-tria refinement for
  ! recoarsening into the macro element
  integer, parameter, public :: MARK_CRS_4TRIA1TRIA         = -2

  ! Mark inner triangle of a 1-tria : 4-tria refinement for
  ! recoarsening into two green triangles, whereby the first vertex
  ! if the inner triangle is connected to the opposite vertex
  integer, parameter, public :: MARK_CRS_4TRIA2TRIA_1       = -3

  ! Mark inner triangle of a 1-tria : 4-tria refinement for
  ! recoarsening into two green triangles, whereby the second vertex
  ! if the inner triangle is connected to the opposite vertex
  integer, parameter, public :: MARK_CRS_4TRIA2TRIA_2       = -4

  ! Mark inner triangle of a 1-tria : 4-tria refinement for
  ! recoarsening into two green triangles, whereby the third vertex
  ! if the inner triangle is connected to the opposite vertex
  integer, parameter, public :: MARK_CRS_4TRIA2TRIA_3       = -5

  ! Mark inner green triangle of a 1-quad : 3-tria refinement for
  ! recoarsening into the macro quadrilateral together with its neighbors
  integer, parameter, public :: MARK_CRS_3TRIA1QUAD         = -6

  ! Mark (left) green triangle of a 1-tria : 2-tria refinement for
  ! recoarsening into the macro triangle together with its right neighbor
  integer, parameter, public :: MARK_CRS_2TRIA1TRIA         = -7

  ! Mark "most inner" green triangle of a 1-quad : 4-tria refinement for
  ! recoarsening into the macro quadrilateral together with its three neighbors
  integer, parameter, public :: MARK_CRS_4TRIA1QUAD         = -8

  ! Mark "most inner" green triangle of a 1-quad : 4-tria refinement for
  ! conversion into three triangles keeping the left neighboring triangle
  integer, parameter, public :: MARK_CRS_4TRIA3TRIA_LEFT    = -9

  ! Mark "most inner" green triangle of a 1-quad : 4-tria refinement for
  ! conversion into three triangles keeping the right neighboring triangle
  integer, parameter, public :: MARK_CRS_4TRIA3TRIA_RIGHT   = -10

  ! Mark green quadrilateral of a 1-quad : 2-quad refinement for
  ! recoarsening into the macro quadrilateral together with its neighbor
  integer, parameter, public :: MARK_CRS_2QUAD1QUAD         = -11
  
  ! Mark green quadrilateral of a 1-quad : 2-quad refinement for
  ! conversion into three triangles keeping the second vertex
  integer, parameter, public :: MARK_CRS_2QUAD3TRIA         = -12

  ! Mark red quadrilateral of a 1-quad : 4-quad refinement for
  ! recoarsening into the macro quadrilateral together with its neighbors
  integer, parameter, public :: MARK_CRS_4QUAD1QUAD         = -13

  ! Mark red quadrilateral of a 1-quad : 4-quad refinement for
  ! recoarsening into two quadrilaterals
  integer, parameter, public :: MARK_CRS_4QUAD2QUAD         = -14

  ! Mark red quadrilateral of a 1-quad : 4-quad refinement for
  ! recoarsening into three triangles
  integer, parameter, public :: MARK_CRS_4QUAD3TRIA         = -15

  ! Mark red quadrilateral of a 1-quad : 4-quad refinement for
  ! recoarsening into four triangles
  integer, parameter, public :: MARK_CRS_4QUAD4TRIA         = -16
  
  ! Mark for keeping element 'as is'
  integer, parameter, public :: MARK_ASIS                   = 0
  integer, parameter, public :: MARK_ASIS_TRIA              = 0
  integer, parameter, public :: MARK_ASIS_QUAD              = 1

  ! Mark element for 1-tria : 2-tria refinement along first edge
  integer, parameter, public :: MARK_REF_TRIA2TRIA_1        = 2

  ! Mark element for 1-tria : 2-tria refinement along second edge
  integer, parameter, public :: MARK_REF_TRIA2TRIA_2        = 4

  ! Mark element for 1-tria : 2-tria refinement along third edge
  integer, parameter, public :: MARK_REF_TRIA2TRIA_3        = 8

  ! Mark element for 1-tria : 3-tria refinement along first and second edge
  integer, parameter, public :: MARK_REF_TRIA3TRIA_12       = 6

  ! Mark element for 1-tria : 3-tria refinement along second and third edge
  integer, parameter, public :: MARK_REF_TRIA3TRIA_23       = 12

  ! Mark element for 1-tria : 3-tria refinement along first and third edge
  integer, parameter, public :: MARK_REF_TRIA3TRIA_13       = 10

  ! Mark element for 1-tria : 4-tria red refinement
  integer, parameter, public :: MARK_REF_TRIA4TRIA          = 14

  ! Mark element for 1-quad : 3-tria refinement along first edge
  integer, parameter, public :: MARK_REF_QUAD3TRIA_1        = 3
  
  ! Mark element for 1-quad : 3-tria refinement along second edge
  integer, parameter, public :: MARK_REF_QUAD3TRIA_2        = 5
  
  ! Mark element for 1-quad : 3-tria refinement along third edge
  integer, parameter, public :: MARK_REF_QUAD3TRIA_3        = 9
  
  ! Mark element for 1-quad : 3-tria refinement along fourth edge
  integer, parameter, public :: MARK_REF_QUAD3TRIA_4        = 17

  ! Mark element for 1-quad : 4-tria refinement along first and second edge
  integer, parameter, public :: MARK_REF_QUAD4TRIA_12       = 7

  ! Mark element for 1-quad : 4-tria refinement along second and third edge
  integer, parameter, public :: MARK_REF_QUAD4TRIA_23       = 13

  ! Mark element for 1-quad : 4-tria refinement along third and fourth edge
  integer, parameter, public :: MARK_REF_QUAD4TRIA_34       = 25

  ! Mark element for 1-quad : 4-tria refinement along first and fourth edge
  integer, parameter, public :: MARK_REF_QUAD4TRIA_14       = 19

  ! Mark element for 1-quad : 4-quad red refinement
  integer, parameter, public :: MARK_REF_QUAD4QUAD          = 31

  ! Mark element for 1-quad : 2-quad refinement along first and third edge
  integer, parameter, public :: MARK_REF_QUAD2QUAD_13       = 11

  ! Mark element for 1-quad : 2-quad refinement along second and fourth edge
  integer, parameter, public :: MARK_REF_QUAD2QUAD_24       = 21

  ! Mark element for 1-quad blue refinement
  integer, parameter, public :: MARK_REF_QUADBLUE_412       = 23
  integer, parameter, public :: MARK_REF_QUADBLUE_234       = 29
  integer, parameter, public :: MARK_REF_QUADBLUE_123       = 15
  integer, parameter, public :: MARK_REF_QUADBLUE_341       = 27
  
!</constantblock>

!<constantblock description="Constants for element states in 2D">

  ! Triangle from the root triangulation
  integer, parameter :: STATE_TRIA_ROOT             = 0

  ! Inner triangle of a 1-tria : 4-tria red refinement
  integer, parameter :: STATE_TRIA_REDINNER         = 14

  ! Outer triangle of a 1-tria : 4-tria red refinement or
  ! outer/inner triangle of 1-quad : 4-tria refinement
  integer, parameter :: STATE_TRIA_OUTERINNER       = 4
  integer, parameter :: STATE_TRIA_OUTERINNER1      = 2   ! only theoretically
  integer, parameter :: STATE_TRIA_OUTERINNER2      = 8   ! only theoretically

  ! Inner triangle of a 1-quad : 3-tria refinement
  integer, parameter :: STATE_TRIA_GREENINNER       = -2

  ! Outer trignale of a green refinement
  integer, parameter :: STATE_TRIA_GREENOUTER_LEFT  = -4
  integer, parameter :: STATE_TRIA_GREENOUTER_RIGHT = -8

  ! Quadrilateral form the root triangulation
  integer, parameter :: STATE_QUAD_ROOT             = 1

  ! Quadrilateral of a 1-quad : 4-quad red refinement
  integer, parameter :: STATE_QUAD_RED1             = 25
  integer, parameter :: STATE_QUAD_RED2             = 19
  integer, parameter :: STATE_QUAD_RED3             = 7
  integer, parameter :: STATE_QUAD_RED4             = 13

  ! Quadrilateral of a 1-quad : 2-quad refinement
  integer, parameter :: STATE_QUAD_HALF1            = 5
  integer, parameter :: STATE_QUAD_HALF2            = 21
  
!</constantblock>

!</constants>

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

  interface add_vertex2D
    module procedure add_vertex_atEdgeMidpoint2D
    module procedure add_vertex_atElementCenter2D
  end interface

  interface replace_element2D
    module procedure replace_elementTria
    module procedure replace_elementQuad
  end interface

  interface add_element2D
    module procedure add_elementTria
    module procedure add_elementQuad
  end interface
  
  interface update_ElementNeighbors2D
    module procedure update_ElemNeighb2D_1to2
    module procedure update_ElemNeighb2D_2to2
  end interface

contains

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_markRefinement2D(rhadapt, rindicator)

!<description>
    ! This subroutine marks all elements that should be refined
    ! due to accuracy reasons. The decision is based on some
    ! indicator vector which must be given element-wise in 2D.
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
    integer, dimension(TRIA_MAXNVE2D) :: IverticesAtElement,IneighboursAtElement
    integer, dimension(TRIA_MAXNVE) :: IvertexAge
    integer  :: ivt,iel,jel,ive,nve,istate,jstate,isubdivide


    ! Check if dynamic data structures are generated and contain data
    if (iand(rhadapt%iSpec, HADAPT_HAS_DYNAMICDATA2D) .ne.&
                            HADAPT_HAS_DYNAMICDATA2D) then
      call output_line('Dynamic data structures are not generated!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_markRefinement2D')
      call sys_halt()
    end if

    ! Initialise marker structure for NEL0 elements
    if (rhadapt%h_Imarker .ne. ST_NOHANDLE) call storage_free(rhadapt%h_Imarker)
    call storage_new('hadapt_markRefinement2D', 'Imarker', rhadapt%NEL0,&
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
        
        ! Get number of vertices per elements
        nve = hadapt_getNVE(rhadapt, iel)

        ! Retrieve local data
        IverticesAtElement(1:nve)   = rhadapt%p_IverticesAtElement(1:nve, iel)
        IneighboursAtElement(1:nve) = rhadapt%p_IneighboursAtElement(1:nve, iel)

        ! An element can only be refined, if all of its vertices do
        ! not exceed the number of admissible subdivision steps.
        ! However, there is an exception to this rule: green elements can
        ! be refined/converted to the corresponding red elements.
        ! So, check the element type and loop over the 3 or 4 vertices.
        select case(nve)
          
        case(TRIA_NVETRI2D)
          ! If triangle has reached maximum number of refinement levels,
          ! then enforce no further refinement of this element
          isubdivide = maxval(abs(rhadapt%p_IvertexAge(&
                                  IverticesAtElement(1:TRIA_NVETRI2D))))

          if (isubdivide .eq. rhadapt%nsubdividemax) then

            ! Check if triangle is an inner/outer red triangle
            do ive = 1, TRIA_NVETRI2D
              IvertexAge(ive) = rhadapt%p_IvertexAge(IverticesAtElement(ive))
            end do
            istate = redgreen_getStateTria(IvertexAge(1:TRIA_NVETRI2D))

            if ((istate .eq. STATE_TRIA_REDINNER) .or.&
                (istate .eq. STATE_TRIA_ROOT)) then

              ! Inner red triangle or root triangle
              ! => no further refinement possible
              p_Imarker(iel) = MARK_ASIS
              
              ! According to the indicator, this element should be refined. Since the
              ! maximum admissible refinement level has been reached no refinement
              ! was performed. At the same time, all vertices of the element should
              ! be "locked" to prevent this element from coarsening
              do ive = 1, TRIA_NVETRI2D
                ivt = IverticesAtElement(ive)
                rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              end do

              ! That is it for element IEL
              cycle mark_elem

            elseif (istate .eq. STATE_TRIA_OUTERINNER) then
              
              ! Possibly outer red triangle, get neighboring triangle
              jel = rhadapt%p_IneighboursAtElement(2, iel)

              ! Check state of neighboring triangle
              do ive = 1, TRIA_NVETRI2D
                IvertexAge(ive) = rhadapt%p_IvertexAge(rhadapt%p_IverticesAtElement(ive, jel))
              end do
              jstate = redgreen_getStateTria(IvertexAge(1:TRIA_NVETRI2D))

              if (jstate .eq. STATE_TRIA_REDINNER) then

                ! Outer red triangle => no further refinement
                p_Imarker(iel) = MARK_ASIS
                
                ! According to the indicator, this element should be refined. Since the
                ! maximum admissible refinement level has been reached no refinement
                ! was performed. At the same time, all vertices of the element should
                ! be "locked" to prevent this element from coarsening
                do ive = 1, TRIA_NVETRI2D
                  ivt = IverticesAtElement(ive)
                  rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
                end do

                ! That is it for element IEL
                cycle mark_elem
              end if
            end if
          end if   ! isubdivide = nsubdividemax

          
          ! Otherwise, we can mark the triangle for refinement
          p_Imarker(iel) = MARK_REF_TRIA4TRIA

          
          ! "Lock" all vertices connected to element IEL since they cannot be removed.
          do ive = 1, TRIA_NVETRI2D
            ivt = IverticesAtElement(ive)
            rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
          end do
          
          
          ! Update number of new vertices. In principle, we can increase the number of
          ! new vertices by 3, i.e., one for each edge. However, some care must be taken
          ! if the edge belongs to some adjacent element which has been marked for
          ! refinement previously. Hence, a new vertex is only added, if the edge is
          ! connected to the boundary or if the adjacent element has not been marked.
          do ive = 1, TRIA_NVETRI2D
            if (IneighboursAtElement(ive) .eq. 0 .or.&
                IneighboursAtElement(ive) .gt. rhadapt%NEL0) then

              ! Edge is adjacent to boundary or has not jet been processed
              rhadapt%increaseNVT = rhadapt%increaseNVT+1
              
            elseif((p_Imarker(IneighboursAtElement(ive)) .ne. MARK_REF_TRIA4TRIA) .and.&
                   (p_Imarker(IneighboursAtElement(ive)) .ne. MARK_REF_QUAD4QUAD)) then
              
              ! Edge has not been marked in previous steps
              rhadapt%increaseNVT = rhadapt%increaseNVT+1
            end if
          end do


        case(TRIA_NVEQUAD2D)
          ! If quadrilateral has reached maximum number of refinement levels,
          ! then enforce no further refinement of this element
          isubdivide = maxval(abs(rhadapt%p_IvertexAge(&
                                  IverticesAtElement(1:TRIA_NVEQUAD2D))))

          if (isubdivide .eq. rhadapt%nsubdividemax) then
            
            ! Check if quadrilateral is a red quadrilateral
            do ive = 1, TRIA_NVEQUAD2D
              IvertexAge(ive) = rhadapt%p_IvertexAge(IverticesAtElement(ive))
            end do
            istate = redgreen_getStateQuad(IvertexAge(1:TRIA_NVEQUAD2D))

            if ((istate .eq. STATE_QUAD_RED4) .or.&
                (istate .eq. STATE_QUAD_ROOT)) then
              
              ! Red quadrilateral => no further refinement
              p_Imarker(iel) = MARK_ASIS
              
              ! According to the indicator, this element should be refined. Since the
              ! maximum admissible refinement level has been reached no refinement
              ! was performed. At the same time, all vertices of the element should
              ! be "locked" to prevent this element from coarsening
              do ive = 1, TRIA_NVEQUAD2D
                ivt = IverticesAtElement(ive)
                rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              end do

              ! That is it for element IEL
              cycle mark_elem
            end if
          end if  ! isubdivide = nsubdividemax
          
          ! Otherwise, we can mark the quadrilateral for refinement
          p_Imarker(iel) = MARK_REF_QUAD4QUAD
          

          ! "Lock" all vertices connected to element IEL since the cannot be removed.
          do ive = 1, TRIA_NVEQUAD2D
            ivt = IverticesAtElement(ive)
            rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
          end do

          
          ! Update number of new vertices. In principle, we can increase the number of
          ! new vertices by 4, i.e., one for each edge. However, some care must be taken
          ! if the edge belongs to some adjacent element which has been marked for
          ! refinement previously. Hence, a new vertex is only added, if the edge is
          ! connected to the boundary or if the adjacent element has not been marked.

          do ive = 1, TRIA_NVEQUAD2D
            if (IneighboursAtElement(ive) .eq. 0 .or. &
                IneighboursAtElement(ive) .gt. rhadapt%NEL0) then

              ! Edge is adjacent to boundary or has not jet been processed
              rhadapt%increaseNVT = rhadapt%increaseNVT+1
              
            elseif((p_Imarker(IneighboursAtElement(ive)) .ne. MARK_REF_TRIA4TRIA) .and.&
                   (p_Imarker(IneighboursAtElement(ive)) .ne. MARK_REF_QUAD4QUAD)) then

              ! Edge has not been marked in previous steps
              rhadapt%increaseNVT = rhadapt%increaseNVT+1
            end if
          end do


        case DEFAULT
          call output_line('Invalid element type!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'hadapt_markRefinement2D')
          call sys_halt()
        end select
        
      else   ! p_Dindicator(IEL) <= drefinementTolerance

        ! Unmark element for refinement
        p_Imarker(iel) = MARK_ASIS

      end if
    end do mark_elem

    ! Set specifier to "marked for refinement"
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_MARKEDREFINE)

  end subroutine hadapt_markRefinement2D

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_markRedgreenCoarsening2D(rhadapt, rindicator)

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
    !
    ! In the context of hierarchical red-green mesh adaptivity, each recoarsening
    ! step is the inverse operation of a previous refinement step. In other words,
    ! the result of a complete recoarsening step could also be generated by locally
    ! refining the initial triangulation. The algorithm is as follows:
    !
    ! 1) All vertices present in the initial triangulation cannot be removed. A
    ! vertex belongs to the initial triangulation if and only if its age is zero.
    !
    ! 2) All other nodes which do not belong to elements that are marked for
    ! refinement are potential candidates for removal. In order to prevent the
    ! removel of such vertices, the nodes must be "locked" by some mechanism.
    !
    ! 3) All three/four nodes of a red triangle/quadrilateral are "locked" if
    ! the element should not be coarsened due to accuracy reasons.
    !
    ! 4) All vertices which belong to the macro element are "locked". Thus it
    ! suffices to consider all edges and "lock" one endpoint if its generation
    ! number is smaller than the generation number of the opposite endpoint.
    !
    ! It remains to prevent the generaiton of blue elements to to nodal removal.
    ! This must be done iteratively since "locking" of one vertex may alter the
    ! situation for all adjacent elements. The generation of blue elements can
    ! be prevented by applying the following rules:
    !
    ! A) If an inner red triangle has two locked vertices,
    ! then the third vertex is locked, too.
    !
    ! B) If the youngest vertex of an outer green triangle (i.e. the vertex which
    ! is responsible for the green refinement) is locked, then lock all nodes
    ! connected to the green triangle. There is one exception to this rule. For
    ! the "inner" green triangles resulting from a 1-quad : 4-tria refinement,
    ! this rule should not be applied since the four triangles can be converted
    ! into three triangles if only one of the two "green vertices" are locked.
    ! (apply 4-tria : 1-quad coarsening and eliminate the handing node by suitable
    ! 1-quad : 3-tria refinement = 4-tria : 3-tria conversion).
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
    integer, dimension(TRIA_MAXNVE2D) :: IverticesAtElement,ImacroElement
    integer, dimension(TRIA_MAXNVE) :: IvertexAge
    integer :: ivt,jvt,iel,jel,jelmid,kel
    integer :: istate,jstate,ive,jve,kve,nve,nlocked
    logical :: bisModified

    ! Check if dynamic data structures are generated and contain data
    if (iand(rhadapt%iSpec, HADAPT_HAS_DYNAMICDATA2D) .ne.&
                            HADAPT_HAS_DYNAMICDATA2D .or.&
        iand(rhadapt%iSpec, HADAPT_MARKEDREFINE) .ne.&
                            HADAPT_MARKEDREFINE) then
      call output_line('Dynamic data structures are not ' // &
                       ' generated or no marker for grid refinement exists!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_markRedgreenCoarsening2D')
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
    ! For each element which is marked for red refinement,
    ! the three/four corner vertices are "locked" unconditionally.
    !
    ! For each element that is not marked for red refinement, the
    ! indicator is considered and all vertices are "locked" if
    ! re-coarsening is not allowed due to accuracy reasons.
    !
    ! There is one special case which needs to be considered.
    ! If a green element is marked for refinement, then it is
    ! converted into the corresponding red element prior to
    ! performing further refinement. Therefore, the intermediate
    ! grid after (!) the "redgreen-marking" procedure may not
    ! satisfy the conformity property, e.g. a quadrilateral that
    ! is marked for red refinement may be adjacent to some former
    ! green element which is (!) already subdivided regularly.
    ! In this case, the "hanging" node at the edge midpoint
    ! needs to be "locked" since it cannot be deleted. Note that
    ! such "hanging" nodes are not "locked" if the adjacent (not
    ! yet subdivided) element is not marked for red refinement.
    phase1: do iel = 1, rhadapt%NEL
      
      ! Check if the current element is marked for red refinement.
      select case(p_Imarker(iel))
        
      case(MARK_REF_QUAD4QUAD,&
           MARK_REF_TRIA4TRIA)

        ! The current element is marked for red refinement.
        ! Thus, we have to "lock" all vertices connected to it.
        ! Moreover, we have to "lock" vertices located at edge
        ! midpoints which may be generated due to conversion of
        ! neighboring green elements into corresponding red ones.
        do ive = 1, hadapt_getNVE(rhadapt, iel)

          ! "Lock" vertex at corner
          ivt = rhadapt%p_IverticesAtElement(ive, iel)
          rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
          
          ! Get numbers of element(s) adjacent to the current edge
          jel    = rhadapt%p_IneighboursAtElement(ive, iel)
          jelmid = rhadapt%p_ImidneighboursAtElement(ive, iel)

          ! Check of neighboring element has been subdivided
          if (jel .ne. jelmid) then

            ! Find element IEL in adjacency list of element JEL
            ! and "lock" the corresponding midpoint vertex
            do jve = 1, hadapt_getNVE(rhadapt, jel)
              if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. iel) then
                jvt = rhadapt%p_IverticesAtElement(jve, jel)
                rhadapt%p_IvertexAge(jvt) = -abs(rhadapt%p_IvertexAge(jvt))
                exit
              end if
            end do

          elseif (jel .ne. 0) then   ! adjacent element has not been
                                     ! subdivided and is not boundary

            ! Lock all vertices in adjacent element
            do jve = 1, hadapt_getNVE(rhadapt, jel)
              jvt = rhadapt%p_IverticesAtElement(jve, jel)
              rhadapt%p_IvertexAge(jvt) = -abs(rhadapt%p_IvertexAge(jvt))
            end do
            
          end if
        end do
        
      case(:MARK_ASIS_QUAD)

        ! The current element is not marked for red refinement, and hence,
        ! it is a potential candidate for the re-coarsening algorithm.
        
        ! Get number of vertices per element
        nve = hadapt_getNVE(rhadapt, iel)

        ! Get local data for element iel
        do ive = 1, nve
          IverticesAtElement(ive) = rhadapt%p_IverticesAtElement(ive, iel)
        end do
        
        ! Check if element should be "locked" due to accuracy reasons.
        if (lockElement(iel)) then
          
          phase1_lock_all: do ive = 1, nve

            ! "Lock" vertex at corner
            ivt = IverticesAtElement(ive)
            rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))

          end do phase1_lock_all
          
        else   ! no need to keep element due to accuray reasons
          
          ! Loop over adjacent edges and check if the two endpoints have
          ! different age. If this is the case then "lock" the one with
          ! smaller generation number since it cannot be removed.
          phase1_lock_older: do ive = 1, nve
            
            ! Get global vertex numbers
            ivt = IverticesAtElement(ive)
            jvt = IverticesAtElement(mod(ive, nve)+1)

            ! Get numbers of element(s) adjacent to the current edge
            jel    = rhadapt%p_IneighboursAtElement(ive, iel)
            jelmid = rhadapt%p_ImidneighboursAtElement(ive, iel)
            
            ! Check if neighboring element has been subdivided
            if (jel .ne. jelmid) then

              ! "Lock" both endpoints of the current edge
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              rhadapt%p_IvertexAge(jvt) = -abs(rhadapt%p_IvertexAge(jvt))
              
            else   ! adjacent element has not bee subdivided

              ! Check if endpoints have different age and "lock" the older one
              if (abs(rhadapt%p_IvertexAge(ivt)) .lt.&
                  abs(rhadapt%p_IvertexAge(jvt))) then
                rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              elseif (abs(rhadapt%p_IvertexAge(jvt)) .lt.&
                  abs(rhadapt%p_IvertexAge(ivt))) then
                rhadapt%p_IvertexAge(jvt) = -abs(rhadapt%p_IvertexAge(jvt))
              end if
              
            end if
          end do phase1_lock_older

        end if   ! lock element due to accuracy reasons?

        
      case DEFAULT

        ! The current element is marked for some sort of green
        ! refinement, and hence, only the corner nodes have to be
        ! locked. Note that the rest is done elsewhere in this
        ! subroutine.
        
        do ive = 1, hadapt_getNVE(rhadapt, iel)
          
          ! "Lock" vertex at corner
          ivt = rhadapt%p_IverticesAtElement(ive, iel)
          rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
        end do
        
      end select
    end do phase1

    
    ! Phase 2: Prevent the generation of blue elements.
    !
    ! In the second phase, vertices are "locked" which would give rise to
    ! the creation of blue refinement pattern. In particular, there are
    ! three different situations which have to be considered:
    !
    ! (a) If two vertices of an inner red triangle are locked,
    !     then the third vertex is also locked
    !
    ! (b) If the "inner" node of a patch of four red quadrilaterals
    !     is locked, then lock all nine vertices of the patch
    !
    ! (c) If three midpoint vertices of a patch of four red quadrilaterals
    !     are locked, then lock all nine vertices of the patch
    !
    ! The locking of vertices may require the re-investigation of surrounding
    ! elements since they may give rise to the creation of blue elements.
    ! In order to circumvent this recursion, we make use of an infinite
    ! do-while loop which is terminated if no vertices have been "modified".
    
    ! Let us start
    bisModified = .true.
    phase2: do while(bisModified)

      ! Ok, hopefully no vertex will be modified
      bisModified = .false.

      ! Loop over all elements
      phase2_elem: do iel = 1, rhadapt%NEL

        ! Get number of vertices per element
        nve = hadapt_getNVE(rhadapt, iel)
        
        ! Get state of element IEL
        select case(nve)
        case(TRIA_NVETRI2D)
          do ive = 1, TRIA_NVETRI2D
            ivt                     = rhadapt%p_IverticesAtElement(ive, iel)
            IverticesAtElement(ive) = ivt
            IvertexAge(ive)         = rhadapt%p_IvertexAge(ivt)
          end do
          istate = redgreen_getStateTria(IvertexAge(1:TRIA_NVETRI2D))


        case(TRIA_NVEQUAD2D)
          do ive = 1, TRIA_NVEQUAD2D
            ivt                     = rhadapt%p_IverticesAtElement(ive, iel)
            IverticesAtElement(ive) = ivt
            IvertexAge(ive)         = rhadapt%p_IvertexAge(ivt)
          end do
          istate = redgreen_getStateQuad(IvertexAge(1:TRIA_NVEQUAD2D))


        case DEFAULT
          call output_line('Invalid number of vertices per element!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'hadapt_markRedgreenCoarsening2D')
          call sys_halt()
        end select


        ! What element state do we have?
        select case(istate)
          
        case(STATE_QUAD_RED4)
          ! Element IEL is one of four red quadrilaterals of a 1-quad : 4-quad refinement.
          
          ! Check if the interior vertex is locked, then lock all vertices in patch
          if (rhadapt%p_IvertexAge(IverticesAtElement(3)) .le. 0) then
            
            ! Lock all other vertices. Note that the first vertex and the third one
            ! are already locked, so that we have to consider vertices 2 and 4.
              
            ! Get global number of second vertex
            ivt = rhadapt%p_IverticesAtElement(2, iel)
            
            if (rhadapt%p_IvertexAge(ivt) .gt. 0) then
              ! Lock second vertex if required
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ! We modified some vertices in this iteration
              bisModified = .true.
            end if

          else   ! Interior vertex is not locked

            ! Count the number of locked midpoint vertices
            jel = iel; nlocked = 0
            do ive = 1, TRIA_NVEQUAD2D
              
              ! Get global number of second vertex
              ivt = rhadapt%p_IverticesAtElement(2, jel)

              ! Increase number of locked vertices if required
              if (rhadapt%p_IvertexAge(ivt) .le. 0) nlocked = nlocked+1

              ! Proceed to adjacent element along second edge
              jel = rhadapt%p_IneighboursAtElement(2, jel)
            end do

            ! Check if three or more vertices are locked, then lock all vertices
            if (nlocked .ge. 3) then
              
              ! Get global number of "inner" vertex
              ivt = rhadapt%p_IverticesAtElement(3, iel)
              
              if (rhadapt%p_IvertexAge(ivt) .gt. 0) then
                ! Lock inner vertex if required
                rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))

                ! We modified some vertices in this iteration
                bisModified = .true.
              end if
              
              ! Lock all midpoint vertices if required
              jel = iel
              do ive = 1, TRIA_NVEQUAD2D

                ! Get global number of second vertex
                ivt = rhadapt%p_IverticesAtElement(2, jel)

                if (rhadapt%p_IvertexAge(ivt) .gt. 0) then
                  ! Lock second vertex if required
                  rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))

                  ! We modified some vertices in this iteration
                  bisModified = .true.
                end if

                ! Proceed to adjacent element along second edge
                jel = rhadapt%p_IneighboursAtElement(2, jel)
              end do
              
            end if   ! nlocked >= 3
            
          end if   ! inner vertex is locked
          
          
        case(STATE_TRIA_REDINNER)
          ! Element IEL is inner red triangle of a 1-tria : 4-tria refinement.
          
          ! Determine number of locked vertices of the inner triangle.
          nlocked = 0
          do ive = 1, TRIA_NVETRI2D
            
            ! Get global number of vertex
            ivt = IverticesAtElement(ive)
            if (rhadapt%p_IvertexAge(ivt) .le. 0) nlocked = nlocked+1
          end do
          
          ! Check if two vertices are locked, then lock all vertices
          if (nlocked .eq. 2) then
            do ive = 1, TRIA_NVETRI2D
              ivt = IverticesAtElement(ive)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
            end do
            
            ! We definitively modified some vertex in this iteration
            bisModified = .true.
          end if


        case DEFAULT

          ! Get number of vertices per element
          nve = hadapt_getNVE(rhadapt, iel)
          
          phase2_edge: do ive = 1, nve
            
            ! Get numbers of element(s) adjacent to the current edge
            jel    = rhadapt%p_IneighboursAtElement(ive, iel)
            jelmid = rhadapt%p_ImidneighboursAtElement(ive, iel)
            
            ! Check if neighboring element has been subdivided
            if (jel .ne. jelmid) then
              
              ! Find element IEL in adjacency list of element JEL
              phase2_edge_midpoint: do jve = 1, hadapt_getNVE(rhadapt, jel)
                if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. iel) then

                  ! Get number of midpoint vertex
                  jvt = rhadapt%p_IverticesAtElement(jve, jel)
            
                  ! Check if midpoint vertex for the current edge is "locked",
                  ! then lock all corner vertices of element IEL if required
                  if (rhadapt%p_IvertexAge(jvt) .le. 0) then
                    phase2_lock_all: do kve = 1, nve
                      ivt = rhadapt%p_IverticesAtElement(kve, iel)
                      if (rhadapt%p_IvertexAge(ivt) .gt. 0) then
                        rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))

                        ! We definitively modified some vertex in this iteration
                        bisModified = .true.
                      end if
                    end do phase2_lock_all

                    exit phase2_edge
                  end if

                  exit phase2_edge_midpoint
                end if
              end do phase2_edge_midpoint
              
            end if
          end do phase2_edge

        end select
      end do phase2_elem
    end do phase2


    ! Phase 3: Determine re-coarsening rule based on "free" vertices.
    !
    ! In the third phase, all elements are considered and the
    ! corresponding re-coarsening rule is adopted based on the
    ! number and type of free vertices. Note that elements which
    ! have been marked for refinement may be modified.

    phase3: do iel = 1, rhadapt%NEL
      
      ! Check if the current element is marked for red refinement.
      select case(p_Imarker(iel))
        
      case(MARK_REF_QUAD4QUAD,&
           MARK_REF_TRIA4TRIA)

        ! The current element is marked for red refinement, and hence,
        ! it cannot be removed by the re-coarsening algorithm.


      case(MARK_REF_TRIA2TRIA_1,&
           MARK_REF_TRIA2TRIA_2,&
           MARK_REF_TRIA2TRIA_3)

        ! Get local edge number [2,4,8] -> [1,2,3]
        ive = ishft(p_Imarker(iel), -2)+1

        ! The current element is a triangle which is marker
        ! for green refinement along one single edge.
        jel    = rhadapt%p_IneighboursAtElement(ive, iel)
        jelmid = rhadapt%p_ImidneighboursAtElement(ive, iel)

        if (jel .ne. jelmid) then
          ! Find element IEL in adjacency list of element JEL
          do jve = 1, hadapt_getNVE(rhadapt, jel)
            if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. iel) then
              jvt = rhadapt%p_IverticesAtElement(jve, jel)
              if (rhadapt%p_IvertexAge(jvt) .gt. 0) then
                p_Imarker(iel) = ibclr(p_Imarker(iel), ive)
              end if
              exit
            end if
          end do
        end if


      case(MARK_REF_QUAD2QUAD_13,&
           MARK_REF_QUAD2QUAD_24)

        ! Get local edge number [11,21] -> [1,2]
        nve = ishft(p_Imarker(iel), -3)

        ! The current element is a quadrilateral which is marker
        ! for green refinement along two opposive edges.
        do ive = nve, nve+2, 2
          
          jel    = rhadapt%p_IneighboursAtElement(ive, iel)
          jelmid = rhadapt%p_ImidneighboursAtElement(ive, iel)
          
          if (jel .ne. jelmid) then
            ! Find element IEL in adjacency list of element JEL
            do jve = 1, hadapt_getNVE(rhadapt, jel)
              if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. iel) then
                jvt = rhadapt%p_IverticesAtElement(jve, jel)
                if (rhadapt%p_IvertexAge(jvt) .gt. 0) then
                  p_Imarker(iel) = ibclr(p_Imarker(iel), ive)
                end if
                exit
              end if
            end do
          end if
        end do
        
        
      case(MARK_REF_QUAD3TRIA_1,&
           MARK_REF_QUAD3TRIA_2,&
           MARK_REF_QUAD3TRIA_3,&
           MARK_REF_QUAD3TRIA_4)
        
        ! Get local edge number [3,5,9,17] -> [1,2,3,4]
        ive = min(4, ishft(p_Imarker(iel), -2)+1)

        ! The current element is a quadrilateral which is
        ! marker for green refinement along one single edge.
        jel    = rhadapt%p_IneighboursAtElement(ive, iel)
        jelmid = rhadapt%p_ImidneighboursAtElement(ive, iel)

        if (jel .ne. jelmid) then
          ! Find element IEL in adjacency list of element JEL
          do jve = 1, hadapt_getNVE(rhadapt, jel)
            if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. iel) then
              jvt = rhadapt%p_IverticesAtElement(jve, jel)
              if (rhadapt%p_IvertexAge(jvt) .gt. 0) then
                p_Imarker(iel) = ibclr(p_Imarker(iel), ive)
              end if
              exit
            end if
          end do
        end if

        
      case(MARK_REF_QUAD4TRIA_12,&
           MARK_REF_QUAD4TRIA_23,&
           MARK_REF_QUAD4TRIA_34,&
           MARK_REF_QUAD4TRIA_14)

        if (p_Imarker(iel) .eq. MARK_REF_QUAD4TRIA_12) then
          ive = 1
        elseif(p_Imarker(iel) .eq. MARK_REF_QUAD4TRIA_23) then
          ive = 2
        elseif(p_Imarker(iel) .eq. MARK_REF_QUAD4TRIA_34) then
          ive = 3
        else
          ive = 4
        end if

        ! The current element is a quadrilateral which is marker
        ! for green refinement along two connected edges.
          
        jel    = rhadapt%p_IneighboursAtElement(ive, iel)
        jelmid = rhadapt%p_ImidneighboursAtElement(ive, iel)
        
        if (jel .ne. jelmid) then
          ! Find element IEL in adjacency list of element JEL
          do jve = 1, hadapt_getNVE(rhadapt, jel)
            if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. iel) then
              jvt = rhadapt%p_IverticesAtElement(jve, jel)
              if (rhadapt%p_IvertexAge(jvt) .gt. 0) then
                p_Imarker(iel) = ibclr(p_Imarker(iel), ive)
              end if
              exit
            end if
          end do
        end if

        ive = mod(ive, 4)+1

        jel    = rhadapt%p_IneighboursAtElement(ive, iel)
        jelmid = rhadapt%p_ImidneighboursAtElement(ive, iel)
        
        if (jel .ne. jelmid) then
          ! Find element IEL in adjacency list of element JEL
          do jve = 1, hadapt_getNVE(rhadapt, jel)
            if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. iel) then
              jvt = rhadapt%p_IverticesAtElement(jve, jel)
              if (rhadapt%p_IvertexAge(jvt) .gt. 0) then
                p_Imarker(iel) = ibclr(p_Imarker(iel), ive)
              end if
              exit
            end if
          end do
        end if


      case DEFAULT
        
        ! The current element is not marked for red refinement, and hence,
        ! it is a potential candidate for the re-coarsening algorithm.
        
        ! Get number of vertices per element
        nve = hadapt_getNVE(rhadapt, iel)
        
        ! Get state of current element
        select case(nve)
        case(TRIA_NVETRI2D)
          do ive = 1, TRIA_NVETRI2D
            ivt                     = rhadapt%p_IverticesAtElement(ive, iel)
            IverticesAtElement(ive) = ivt
            IvertexAge(ive)         = rhadapt%p_IvertexAge(ivt)
          end do
          istate = redgreen_getStateTria(IvertexAge(1:TRIA_NVETRI2D))
          
          
        case(TRIA_NVEQUAD2D)
          do ive = 1, TRIA_NVEQUAD2D
            ivt                     = rhadapt%p_IverticesAtElement(ive, iel)
            IverticesAtElement(ive) = ivt
            IvertexAge(ive)         = rhadapt%p_IvertexAge(IverticesAtElement(ive))
          end do
          istate = redgreen_getStateQuad(IvertexAge(1:TRIA_NVEQUAD2D))
          
          
        case DEFAULT
          call output_line('Invalid number of vertices per element!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'hadapt_markRedgreenCoarsening2D')
          call sys_halt()
        end select
        
        
        ! What state are we?
        select case(istate)
        case(STATE_TRIA_ROOT,&
             STATE_QUAD_ROOT,&
             STATE_TRIA_OUTERINNER,&
             STATE_TRIA_GREENINNER,&
             STATE_TRIA_GREENOUTER_LEFT)
          ! Dummy states which are required to prevent
          ! immediate stop in the default branch.
          
          
        case(STATE_TRIA_REDINNER)
          ! If all vertices of the element are free, then the inner red
          ! triangle can be combined with its three neighboring red triangles.
          ! If one vertex of the inner triangle is locked, then the inner
          ! red tiangle together with its three neighboring elements can
          ! be converted into two green triangles.
          
          ! Convert position of locked vertices into integer
          nlocked = 0
          do ive = 1, TRIA_NVETRI2D
            ivt = IverticesAtElement(ive)
            if (rhadapt%p_IvertexAge(ivt) .le. 0) nlocked = ibset(nlocked, ive)
          end do
          
          ! How many vertices are locked?
          select case(nlocked)
          case(0)
            ! Mark element IEL for re-coarsening into macro element
            p_Imarker(iel) = MARK_CRS_4TRIA1TRIA
            
          case(2)
            ! Mark element IEL for re-coarsening into two triangles
            p_Imarker(iel) = MARK_CRS_4TRIA2TRIA_1
            
          case(4)
            ! Mark element IEL for re-coarsening into two triangles
            p_Imarker(iel) = MARK_CRS_4TRIA2TRIA_2
            
          case(8)
            ! Mark element IEL for re-coarsening into two triangles
            p_Imarker(iel) = MARK_CRS_4TRIA2TRIA_3
            
          case DEFAULT
            cycle phase3
          end select
          
          ! Mark all adjacent elements 'as is'
          do ive = 1, TRIA_NVETRI2D
            jel = rhadapt%p_IneighboursAtElement(ive, iel)
            p_Imarker(jel) = MARK_ASIS_TRIA
          end do
          
          
        case(STATE_QUAD_RED4)
          ! This is one of four quadrilaterals which result from a 1-quad : 4-quad
          ! refinement. Here, the situation is slightly more difficult because
          ! multiple rotationally-symmetric situations have to be considered. The
          ! four elements forming the macro element can be collected by visiting the
          ! neighboring element along the second edge starting at element IEL.
          
          ! Check if the "inner" vertex is locked, then no re-coarsening is possible.
          ivt = rhadapt%p_IverticesAtElement(3, iel)
          if (rhadapt%p_IvertexAge(ivt) .le. 0) cycle phase3
          
          ! Otherwise, determine the number of locked midpoints. This is done
          ! by visiting all four elements and checking the second vertex individually.
          ! Note that we do not only count the number of locked vertices but also
          ! remember its position by setting or clearing the corresponding bit of
          ! the integer nlocked. In the same loop, we determine the numbers of the
          ! four elements of the macro element. Of course, the first is element IEL
          ImacroElement(1) = iel
          
          ! Get global number of second vertex in element IEL
          ivt = rhadapt%p_IverticesAtElement(2, iel)
          nlocked = merge(2, 0, rhadapt%p_IvertexAge(ivt) .le. 0)
          do ive = 2, TRIA_NVEQUAD2D
            
            ! Get neighbouring element and store it
            jel                = rhadapt%p_IneighboursAtElement(2, ImacroElement(ive-1))
            ImacroElement(ive) = jel
            
            ! Get global number of second vertex in element JEL
            jvt = rhadapt%p_IverticesAtElement(2, jel)
            if (rhadapt%p_IvertexAge(jvt) .le. 0) nlocked = ibset(nlocked, ive)
          end do
          
          ! How many vertices are locked?
          select case(nlocked)
          case(0)
            ! Mark element IEL for re-coarsening into macro
            ! element and  keep all remaining elements 'as is'.
            p_Imarker(ImacroElement(1)) = MARK_CRS_4QUAD1QUAD
            p_Imarker(ImacroElement(2)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(3)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(4)) = MARK_ASIS_QUAD
            
            
          case(2)
            ! There is one vertex locked which is the second vertex
            ! of the first element. All other elements are kept 'as is'.
            p_Imarker(ImacroElement(1)) = MARK_CRS_4QUAD3TRIA
            p_Imarker(ImacroElement(2)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(3)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(4)) = MARK_ASIS_QUAD
            
            
          case(4)
            ! There is one vertex locked which is the second vertex
            ! of the second element. All other elements are kept 'as is'.
            p_Imarker(ImacroElement(1)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(2)) = MARK_CRS_4QUAD3TRIA
            p_Imarker(ImacroElement(3)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(4)) = MARK_ASIS_QUAD
            
            
          case(8)
            ! There is one vertex locked which is the second vertex
            ! of the third element. All other elements are kept 'as is'.
            p_Imarker(ImacroElement(1)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(2)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(3)) = MARK_CRS_4QUAD3TRIA
            p_Imarker(ImacroElement(4)) = MARK_ASIS_QUAD
            
            
          case(16)
            ! There is one vertex locked which is the second vertex
            ! of the fourth element. All other elements are kept 'as is'.
            p_Imarker(ImacroElement(1)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(2)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(3)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(4)) = MARK_CRS_4QUAD3TRIA
            
            
          case(10)
            ! There are two vertices locked which are the second vertices
            ! of the first and third elements. Mark the  first element
            ! for recoarsening. All other elements are kept 'as is'.
            p_Imarker(ImacroElement(1)) = MARK_CRS_4QUAD2QUAD
            p_Imarker(ImacroElement(2)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(3)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(4)) = MARK_ASIS_QUAD
            
            
          case(20)
            ! There are two vertices locked which are the second vertices
            ! of the second and fourth elements. Mark the second element
            ! for recoarsening. All other elements are kept 'as is'.
            p_Imarker(ImacroElement(1)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(2)) = MARK_CRS_4QUAD2QUAD
            p_Imarker(ImacroElement(3)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(4)) = MARK_ASIS_QUAD
            
            
          case(18)
            ! There are two vertices locked which are the second and
            ! fourth vertices of the first elements. Mark the firth element
            ! for recoarsening. All other elements are kept 'as is'.
            p_Imarker(ImacroElement(1)) = MARK_CRS_4QUAD4TRIA
            p_Imarker(ImacroElement(2)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(3)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(4)) = MARK_ASIS_QUAD
            
            
          case(6)
            ! There are two vertices locked which are the second and
            ! fourth vertices of the second elements. Mark the second element
            ! for recoarsening. All other elements are kept 'as is'.
            p_Imarker(ImacroElement(1)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(2)) = MARK_CRS_4QUAD4TRIA
            p_Imarker(ImacroElement(3)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(4)) = MARK_ASIS_QUAD
            
            
          case(12)
            ! There are two vertices locked which are the second and
            ! fourth vertices of the third elements. Mark the third element
            ! for recoarsening. All other elements are kept 'as is'.
            p_Imarker(ImacroElement(1)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(2)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(3)) = MARK_CRS_4QUAD4TRIA
            p_Imarker(ImacroElement(4)) = MARK_ASIS_QUAD
            
            
          case(24)
            ! There are two vertices locked which are the second and
            ! fourth vertices of the fourth elements. Mark the fourth element
            ! for recoarsening. All other elements are kept 'as is'.
            p_Imarker(ImacroElement(1)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(2)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(3)) = MARK_ASIS_QUAD
            p_Imarker(ImacroElement(4)) = MARK_CRS_4QUAD4TRIA
            
            
          case DEFAULT
            call output_line('Invalid number of locked vertices!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'hadapt_markRedgreenCoarsening2D')
            call sys_halt()
          end select
          
          
        case(STATE_TRIA_GREENOUTER_RIGHT)
          ! We have to considered several situations depending on the
          ! state of the adjacent element that shares the second edge.
          jel = rhadapt%p_IneighboursAtElement(2, iel)
          
          do ive = 1, TRIA_NVETRI2D
            jvt             = rhadapt%p_IverticesAtElement(ive, jel)
            IvertexAge(ive) = rhadapt%p_IvertexAge(jvt)
          end do
          jstate = redgreen_getStateTria(IvertexAge(1:TRIA_NVETRI2D))
          
          ! What is the state of the adjacent element?
          select case(jstate)
          case(STATE_TRIA_GREENOUTER_LEFT)
            ! If all vertices of element JEL are locked, then it is kept 'as is'.
            ! Since a left green triangle shares two vertices with its macro
            ! element, it suffices to check the state of the second vertex.
            ivt = rhadapt%p_IverticesAtElement(2, jel)
            if (rhadapt%p_IvertexAge(ivt) .gt. 0) then
              
              ! Mark element for re-coarsening.
              p_Imarker(jel) = MARK_CRS_2TRIA1TRIA
              
              ! Keep its neighbour, i.e. element IEL, 'as is'.
              p_Imarker(iel) = MARK_ASIS_TRIA
            end if
            
            
          case(STATE_TRIA_GREENINNER)
            ! If all vertices of element JEL are locked, then it is kept 'as is'.
            ! Since an inner green triangle shares two vertices with its macro
            ! element, it suffices to check the state of the first vertex.
            ivt = rhadapt%p_IverticesAtElement(1, jel)
            if (rhadapt%p_IvertexAge(ivt) .gt. 0) then
              
              ! Mark element for re-coarsening.
              p_Imarker(jel) = MARK_CRS_3TRIA1QUAD
              
              ! Keep its neighbour along first edge 'as is'.
              kel = rhadapt%p_IneighboursAtElement(1, jel)
              p_Imarker(kel) = MARK_ASIS_TRIA
              
              ! Keep its neighbour along third edge 'as is'.
              kel = rhadapt%p_IneighboursAtElement(3, jel)
              p_Imarker(kel) = MARK_ASIS_TRIA
            end if
            
            
          case(STATE_TRIA_OUTERINNER)
            ! If all vertices of element JEL are locked, then it is kept 'as is'.
            ! Since an inner-outer green triangle shares two vertices with its
            ! macro element, it suffices to check the state of the first vertex.
            ivt = rhadapt%p_IverticesAtElement(2, jel)
            jvt = rhadapt%p_IverticesAtElement(3, jel)
            
            if (rhadapt%p_IvertexAge(ivt) .gt. 0) then
              if (rhadapt%p_IvertexAge(jvt) .gt. 0) then
                
                ! Mark element for re-coarsening.
                p_Imarker(jel) = MARK_CRS_4TRIA1QUAD
                
                ! Keep its neighbour along first edge 'as is'
                kel = rhadapt%p_IneighboursAtElement(1, jel)
                p_Imarker(kel) = MARK_ASIS_TRIA
                
                ! Keep its neighbour along second edge 'as is'
                kel = rhadapt%p_IneighboursAtElement(2, jel)
                p_Imarker(kel) = MARK_ASIS_TRIA
                
                ! Keep its neighbour, i.e. element IEL, along third edge 'as is'
                p_Imarker(iel) = MARK_ASIS_TRIA
                
              else   ! third vertex is locked
                
                ! Mark element for re-coarsening, whereby the green
                ! outer triangle to the left is preserved.
                p_Imarker(jel) = MARK_CRS_4TRIA3TRIA_LEFT ! RIGHT
                
                ! Keep its neighbour along first edge 'as is'
                kel = rhadapt%p_IneighboursAtElement(1, jel)
                p_Imarker(kel) = MARK_ASIS_TRIA
                
                ! Keep its neighbour along second edge 'as is'
                kel = rhadapt%p_IneighboursAtElement(2, jel)
                p_Imarker(kel) = MARK_ASIS_TRIA
                
                ! Keep its neighbour, i.e. element IEL, along third edge 'as is'
                p_Imarker(iel) = MARK_ASIS_TRIA
                
              end if
              
            else    ! second vertex is locked
              
              if (rhadapt%p_IvertexAge(jvt) .gt. 0) then
                
                ! Mark element for re-coarsening, whereby the green
                ! outer triangle to the right is preserved.
                p_Imarker(jel) = MARK_CRS_4TRIA3TRIA_RIGHT ! LEFT
                
                ! Keep its neighbour along first edge 'as is'
                kel = rhadapt%p_IneighboursAtElement(1, jel)
                p_Imarker(kel) = MARK_ASIS_TRIA
                
                ! Keep its neighbour along second edge 'as is'
                kel = rhadapt%p_IneighboursAtElement(2, jel)
                p_Imarker(kel) = MARK_ASIS_TRIA
                
                ! Keep its neighbour, i.e. element IEL, along third edge 'as is'
                p_Imarker(iel)=MARK_ASIS_TRIA
                
              else   ! third vertex is locked
                
                ! Keep all elements in patch 'as is'.
                p_Imarker(iel) = MARK_ASIS_TRIA
                p_Imarker(jel) = MARK_ASIS_TRIA
                
                kel = rhadapt%p_IneighboursAtElement(1, jel)
                p_Imarker(kel) = MARK_ASIS
                
                kel = rhadapt%p_IneighboursAtElement(2, jel)
                p_Imarker(kel) = MARK_ASIS
                
              end if
              
            end if
            
            
          case DEFAULT
            call output_line('Invalid element state of adjacent element!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'hadapt_markRedgreenCoarsening2D')
            call sys_halt()
          end select
          

        case (STATE_QUAD_HALF1,STATE_QUAD_HALF2)
          ! This is one of two quadrilaterals which result from a
          ! 1-quad : 2-quad refinement. If both common vertices are
          ! locked, then the green quadrilaterals are kept 'as is'. If
          ! only one common vertex is locked, the two quadrilaterals
          ! can be coarsened into three triangles reducing the number
          ! of vertices by one. If both common vertices can be
          ! deleted, then both quadrilaterals are removed and the
          ! original macro element is restored.

          ! Get the number of the neighboring green element
          jel = rhadapt%p_IneighboursAtElement(2,iel)
          
          ! Check if the first common vertex of is locked
          if (rhadapt%p_IvertexAge(IverticesAtElement(2)) .le. 0) then

            ! Check if the second common vertex is also locked
            if (rhadapt%p_IvertexAge(IverticesAtElement(3)) .le. 0) then
              ! Keep both elements 'as is'
              p_Imarker(iel)=MARK_ASIS
              p_Imarker(jel)=MARK_ASIS
            else
              ! Mark element IEL for re-coarsening into three triangles,
              ! whereby the second vertex of element IEL is preserved
              p_Imarker(iel)=MARK_CRS_2QUAD3TRIA
              p_Imarker(jel)=MARK_ASIS
            end if

          else

            ! Check if the second common vertex is locked
            if (rhadapt%p_IvertexAge(IverticesAtElement(3)) .le. 0) then
              ! Mark element JEL for re-coarsening into three triangles,
              ! whereby the second vertex of element JEL is preserved
              p_Imarker(iel)=MARK_ASIS
              p_Imarker(jel)=MARK_CRS_2QUAD3TRIA
            else
              ! Mark element IEL for re-coarsening into the macro element
              p_Imarker(iel)=MARK_CRS_2QUAD1QUAD
              p_Imarker(jel)=MARK_ASIS
            end if
          end if
          
          
        case DEFAULT
          call output_line('Invalid element state!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'hadapt_markRedgreenCoarsening2D')
          call sys_halt()
        end select

      end select
    end do phase3
   
    ! Set specifier to "marked for coarsening"
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_MARKEDCOARSEN)

  contains

    ! Here, the real working routines follow.

    !**************************************************************
    ! For a given element IEL, check if all corners vertices need
    ! to be locked. In particular, this is the case if the element
    ! is a red one or belongs to the initial triangulation and, in
    ! addition, the element marker requires to keep this element.

    function lockElement(iel) result (blockElement)

      integer, intent(in) :: iel
      logical :: blockElement

      integer :: istate
      
      ! Check if element marker is available for this element
      blockElement = (iel .le. size(p_Dindicator, 1))
      
      if (blockElement) then
        
        ! Get element state
        istate = redgreen_getState2D(rhadapt, iel)
        
        select case(istate)

        case (STATE_TRIA_ROOT,&
              STATE_TRIA_REDINNER,&
              STATE_QUAD_ROOT,&
              STATE_QUAD_RED4)

          ! Element is a red element and/or belongs to the initial triangulation
          blockElement = (p_Dindicator(iel) .ge. rhadapt%dcoarseningTolerance)

        case DEFAULT
          blockElement = .false.
          
        end select
      end if

    end function lockElement

  end subroutine hadapt_markRedgreenCoarsening2D

  ! ***************************************************************************

!<subroutine>
  
  subroutine hadapt_markRedgreenRefinement2D(rhadapt, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine initialises the adaptive data structure for red-green
    ! refinement. Starting from the marker array the neighbors of elements
    ! which are marked for refinement are recursively marked for refinement
    ! until the resulting mesh satiesfies global conformity.
    ! Note that this subroutine is implemented in an iterative fashion
    ! rather than making use of recursive subroutine calls.
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
    integer, dimension(:), pointer :: p_Imarker,p_Imodifier
    integer, dimension(TRIA_MAXNVE2D) :: IverticesAtElement
    integer :: nvt,ivt,nel,iel,jel,kel,lel,iel1,iel2
    integer :: ive,jve,nve,mve,istate,jstate,kstate
    integer :: h_Imodifier,imodifier
    logical :: bisFirst,bisModified
    integer, dimension(TRIA_MAXNVE) :: IvertexAge
    
    !---------------------------------------------------------------------------
    ! At the moment, only those cells are marked for regular refinement
    ! for which the cell indicator does not satisfy some prescribed treshold.
    ! In order to globally restore the conformity of the final grid, we have
    ! to loop repeatedly over all elements and check if each of the three/four
    ! adjacent edges satisfies conformity.
    ! Some care must be taken for green elements. In order to retain the
    ! shape quality of the elements, green elements need to be converted into
    ! red ones prior to performing further refinement.
    !--------------------------------------------------------------------------
    !
    ! In the course of element marking, some green elements may have to be
    ! converted into regularly refined ones. In the worst case, i.e.
    ! Quad2Quad-refinement, each green element gives rise to 2.5 vertices
    ! and 2 new elements. Hence, the nodal arrays p_IvertexAge and
    ! p_InodalProperty as well as the element arrays p_IverticesAtElement,
    ! p_IneighboursAtElement, p_ImidneighboursAtElement and p_Imarker are
    ! adjusted to the new dimensione.
    !
    if (rhadapt%nGreenElements .gt. 0) then

      ! Predict new dimensions for the worst case
      nvt = rhadapt%NVT0 + ceiling(5*rhadapt%nGreenElements/2._DP)
      nel = rhadapt%NEL0 + rhadapt%nGreenElements

      ! Adjust nodal/elemental arrays
      call storage_realloc('hadapt_markRedgreenRefinement2D', nvt,&
                           rhadapt%h_IvertexAge, ST_NEWBLOCK_ZERO, .true.)
      call storage_realloc('hadapt_markRedgreenRefinement2D', nvt,&
                           rhadapt%h_InodalProperty, ST_NEWBLOCK_ZERO, .true.)
      call storage_realloc('hadapt_markRedgreenRefinement2D', nel,&
                           rhadapt%h_Imarker, ST_NEWBLOCK_ZERO, .true.)
      call storage_realloc('hadapt_markRedgreenRefinement2D', nel,&
                           rhadapt%h_IverticesAtElement, ST_NEWBLOCK_NOINIT, .true.)
      call storage_realloc('hadapt_markRedgreenRefinement2D', nel,&
                           rhadapt%h_IneighboursAtElement, ST_NEWBLOCK_NOINIT, .true.)
      call storage_realloc('hadapt_markRedgreenRefinement2D', nel,&
                           rhadapt%h_ImidneighboursAtElement, ST_NEWBLOCK_NOINIT, .true.)
      
      ! Reset pointers
      call storage_getbase_int(rhadapt%h_IvertexAge,&
                               rhadapt%p_IvertexAge)
      call storage_getbase_int(rhadapt%h_InodalProperty,&
                               rhadapt%p_InodalProperty)
      call storage_getbase_int2D(rhadapt%h_IverticesAtElement,&
                                 rhadapt%p_IverticesAtElement)
      call storage_getbase_int2D(rhadapt%h_IneighboursAtElement,&
                                 rhadapt%p_IneighboursAtElement)
      call storage_getbase_int2D(rhadapt%h_ImidneighboursAtElement,&
                                 rhadapt%p_ImidneighboursAtElement)

      ! Adjust dimension of solution vector
      if (present(fcb_hadaptCallback) .and. present(rcollection)) then
        rcollection%IquickAccess(1) = nvt
        call fcb_hadaptCallback(HADAPT_OPR_ADJUSTVERTEXDIM, rcollection)
      end if

      ! Create new array for modifier
      call storage_new('hadapt_markRedgreenRefinement2D', 'p_Imodifier', nel,&
                       ST_INT, h_Imodifier, ST_NEWBLOCK_ZERO)
      call storage_getbase_int(h_Imodifier, p_Imodifier)
      call storage_getbase_int(rhadapt%h_Imarker, p_Imarker)

    else
      
      ! No green elements have to be considered, hence use NEL0
      call storage_new('hadapt_markRedgreenRefinement2D', 'p_Imodifier', rhadapt%NEL0,&
                        ST_INT, h_Imodifier, ST_NEWBLOCK_ZERO)
      call storage_getbase_int(h_Imodifier, p_Imodifier)
      call storage_getbase_int(rhadapt%h_Imarker, p_Imarker)

    end if

    ! The task of the two arrays p_Imarker and p_Imodifier are as follows:
    !
    ! - p_Imarker stores for each element its individual state from which
    !   the task of refinement can be uniquely determined (see above;
    !   data type t_hadapt).
    !   If the state of one element is modified, then potentionally all
    !   neighboringelements have to be checked, if conformity is violated.
    !   In order to prevent this (expensive) check to be performed for all
    !   elements again and again, only those elements are investigated for
    !   which the modifier p_Imodifier is active, i.e.
    !   p_Imodifier(iel) .EQ. imodifier. In order to pre-select elements
    !   for investigation in the next loop, p_Imodifier(iel) is set to
    !   -imodifier and the modifier is reversed after each iteration.
    !
    !   This complicated step is necessary for the following reason:
    !
    !   If one element with mediate number initially is not marked but
    !   gets marked for subdivision at one edge, then it will be considered
    !   in the same iteration and green elements will eventually converted.
    !   However, it is possible that the same element would be marked for
    !   subdivision at another edge due to an adjacent element with large
    !   number. Hence, it is advisable to first process all elements which
    !   are activated and pre-select their neighbors for the next iteration.
    
    ! Ok, so let us start.
    imodifier = 1
    
    ! Initially, set the modified for all elements of the given triangulation
    ! which are  marked for refinement due to accuracy reasons.
    do iel = 1, rhadapt%NEL0
      if (p_Imarker(iel) .ne. MARK_ASIS) p_Imodifier(iel) = imodifier
    end do

    ! Initialization
    bisModified = .true.
    bisFirst    = .true.
    
    ! Iterate until global conformity is achieved
    conform: do while(bisModified)
      
      ! Ok, hopefully nothing gets modified
      bisModified = .false.
      
      ! Loop over all elements present in the initial grid (IEL <= NEL0)
      ! which are modified and mark their neighbors for further refinement
      ! if conformity is violated for some edge.
      element: do iel = 1, rhadapt%NEL0
        
        ! Skip those element which have not been modified
        if (p_Imodifier(iel) .ne. imodifier) cycle element


        ! Get number of vertices per element
        nve = hadapt_getNVE(rhadapt, iel)
        
        ! Get local data for element IEL
        do ive = 1, nve
          IverticesAtElement(ive) = rhadapt%p_IverticesAtElement(ive, iel)
        end do
      
        ! Are we triangular or quadrilateral element?
        select case(nve)
        case(TRIA_NVETRI2D)
          ! Triangular element: Clear bit0
          p_Imarker(iel) = ibclr(p_Imarker(iel), 0)
          
          ! Get element state
          do ive = 1, TRIA_NVETRI2D
            IvertexAge(ive) = rhadapt%p_IvertexAge(IverticesAtElement(ive))
          end do
          istate = redgreen_getStateTria(IvertexAge(1:TRIA_NVETRI2D))

        case(TRIA_NVEQUAD2D)
          ! Quadrilateral element: Set bit0
          p_Imarker(iel) = ibset(p_Imarker(iel), 0)

          ! Get element state
          do ive = 1, TRIA_NVEQUAD2D
            IvertexAge(ive) = rhadapt%p_IvertexAge(IverticesAtElement(ive))
          end do
          istate = redgreen_getStateQuad(IvertexAge(1:TRIA_NVEQUAD2D))

        case DEFAULT
          call output_line('Invalid number of vertices per element!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'hadapt_markRedgreenRefinement2D')
          call sys_halt()
        end select


        !-----------------------------------------------------------------------
        ! Check if the state of the current element IEL allows direct
        ! refinement or if the elements needs to be "converted" first
        !-----------------------------------------------------------------------
        select case(istate)
        case (STATE_TRIA_ROOT,&
              STATE_QUAD_ROOT,&
              STATE_TRIA_REDINNER,&
              STATE_QUAD_RED1,&
              STATE_QUAD_RED2,&
              STATE_QUAD_RED3,&
              STATE_QUAD_RED4)
          ! States which can be directly accepted for local refinement


        case(STATE_TRIA_OUTERINNER1,&
             STATE_TRIA_OUTERINNER2)
          ! Theoretically, these states may occure and have to be treated
          ! like CASE(4), see below. Due to the fact, that these states can
          ! only be generated by means of local red-green refinement and are
          ! not present in the initial grid we make use of additional know-
          ! ledge: Triangles which have too youngest vertices can either
          ! be outer red elements resulting from a Tria4Tria refinement
          ! or one of the opposite triangles resulting from Quad4Tria refinement.
          ! In all cases, the edge which connects the two nodes is opposite
          ! to the first local vertex. Hence, work must only be done for CASE(4)
          call output_line('These states must not occur!',&
                           OU_CLASS_ERROR, OU_MODE_STD, 'hadapt_markRedgreenRefinement2D')
          call sys_halt()

          
        case(STATE_TRIA_OUTERINNER)
          ! The triangle can either be an outer red element resulting from a Tria4Tria
          ! refinement or one of the two opposite triangles resulting from Quad4Tria
          ! refinement. This can be easily checked. If the opposite element is not
          ! an inner red triangle (14), then the element IEL and its opposite neighbor
          ! make up the inner "diamond" resulting from a Quad4Tria refinement.
          jel = rhadapt%p_IneighboursAtElement(2, iel)
       
          ! First, we need to check a special case. If the second edge of element IEL
          ! has two different neighbors, then the original neighbor along this edge
          ! was a green triangle that has been converted into an (outer) red one. In
          ! this case, the current element IEL does not need further modifications.
          if (jel .ne. rhadapt%p_ImidneighboursAtElement(2, iel)) goto 100
          
          ! Otherwise, determine the state of the edge neighbor JEL
          jstate = redgreen_getState2D(rhadapt, jel)
          
          if (jstate .eq. STATE_TRIA_OUTERINNER) then
            
            ! We know that element IEL and JEL make up the inner "diamond" resulting
            ! from a Quad4Tria refinement. It remains to find the two green triangle
            ! which make up the outer part of the macro element.
            ! Since we do not know if they are adjacent to element IEL or JEL we
            ! have to perform additional checks:
            !
            ! - The element along the first edge of the inner triangle
            !   must have state STATE_TRIA_GREENOUTER_LEFT
            !
            kel    = rhadapt%p_IneighboursAtElement(1, iel)
            kstate = redgreen_getState2D(rhadapt, kel)

            if (kstate .ne. STATE_TRIA_GREENOUTER_LEFT .or.&
                rhadapt%p_IneighboursAtElement(2, kel) .ne. iel) then
              ! At this stage we can be sure that element IEL is not (!)
              ! the inner triangle, hence, it must be element JEL
              
              ! We need to find the two missing triangles KEL and LEL
              ! which make up the original quadrilateral
              kel = rhadapt%p_IneighboursAtElement(1, jel)
              lel = rhadapt%p_IneighboursAtElement(3, jel)

              ! Mark the edge of the element adjacent to KEL for subdivision
              iel1 = rhadapt%p_IneighboursAtElement(3, kel)
              iel2 = rhadapt%p_ImidneighboursAtElement(3, kel)
              if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(kel, iel1)

              ! Mark the edge of the element adjacent to LEL for subdivision
              iel1 = rhadapt%p_IneighboursAtElement(1, lel)
              iel2 = rhadapt%p_ImidneighboursAtElement(1, lel)
              if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(lel, iel1)

              ! Physically convert the four triangles into four quadrilaterals
              call convert_Quad4Tria(rhadapt, kel, iel, lel, jel,&
                                     rcollection, fcb_hadaptCallback)

              ! Reduce the number of green elements by four
              rhadapt%nGreenElements = rhadapt%nGreenElements-4
              
              ! We modified some elements in this iteration
              bisModified = .true.


              ! All edges of element JEL must be unmarked.
              p_Imarker(jel) = ibset(0, 0)


              ! The second and third edge of element KEL must be unmarked.
              ! Moreover, the first edge of element KEL is marked for
              ! refinement if it is also marked from the adjacent element.
              iel1 = rhadapt%p_IneighboursAtElement(1, kel)
              iel2 = rhadapt%p_ImidneighboursAtElement(1, kel)
              if (ismarked_edge(iel1, iel2, kel)) then
                p_Imarker(kel) = ibset(ibset(0, 1), 0)
              else
                p_Imarker(kel) = ibset(0, 0)
              end if


              ! The second and third edge of element IEL must be unmarked.
              p_Imarker(iel) = ibset(0, 0)

              ! The fourth edge of element IEL is only marked if it
              ! is also marked from the adjacent element
              iel1 = rhadapt%p_IneighboursAtElement(4, iel)
              iel2 = rhadapt%p_ImidneighboursAtElement(4, iel)
              if (ismarked_edge(iel1, iel2, iel)) p_Imarker(iel) = ibset(p_Imarker(iel), 4)
              
              ! The first edge of element IEL is only marked if it
              ! is also marked from the adjacent element
              iel1 = rhadapt%p_IneighboursAtElement(1, iel)
              iel2 = rhadapt%p_ImidneighboursAtElement(1, iel)
              if (ismarked_edge(iel1, iel2, iel)) p_Imarker(iel) = ibset(p_Imarker(iel), 1)
                      

              ! The first, second and third edge of element LEL must be unmarked.
              p_Imarker(lel) = ibset(0,0)
              
              ! The fourth edge of element LEL is only marked if it
              ! is  also marked from the adjacent element
              iel1 = rhadapt%p_IneighboursAtElement(4, lel)
              iel2 = rhadapt%p_ImidneighboursAtElement(4, lel)
              if (ismarked_edge(iel1, iel2, lel)) p_Imarker(lel) = ibset(p_Imarker(lel), 4)
              
            else

              ! Otherwise, element IEL is (!) the inner triangle.
              
              ! We need to find the two missing triangles KEL and LEL
              ! which make up the original quadrilateral
              kel = rhadapt%p_IneighboursAtElement(1, iel)
              lel = rhadapt%p_IneighboursAtElement(3, iel)

              ! Mark the edge of the element adjacent to KEL for subdivision
              iel1 = rhadapt%p_IneighboursAtElement(3, kel)
              iel2 = rhadapt%p_ImidneighboursAtElement(3, kel)
              if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(kel, iel1)
              
              ! Mark the edge of the element adjacent to LEL for subdivision
              iel1 = rhadapt%p_IneighboursAtElement(1, lel)
              iel2 = rhadapt%p_ImidneighboursAtElement(1, lel)
              if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(lel, iel1)

              ! Physically convert the four triangles into four quadrilaterals
              call convert_Quad4Tria(rhadapt, kel, jel, lel, iel,&
                                     rcollection, fcb_hadaptCallback)

              ! Reduce the number of green elements by four
              rhadapt%nGreenElements = rhadapt%nGreenElements-4
              
              ! We modified some elements in this iteration
              bisModified = .true.

              
              ! All edges of element IEL must be unmarked.
              p_Imarker(iel) = ibset(0,0)
              
              ! The second and third edge of element KEL must be unmarked.
              ! Moreover, the first edge if element KEL is marked for
              ! refinement if it is also marked from the adjacent element
              iel1 = rhadapt%p_IneighboursAtElement(1, kel)
              iel2 = rhadapt%p_ImidneighboursAtElement(1, kel)
              if (ismarked_edge(iel1, iel2, kel)) then
                p_Imarker(kel) = ibset(ibset(0,1),0)
              else
                p_Imarker(kel) = ibset(0,0)
              end if


              ! The second and third edge of element JEL must be unmarked.
              p_Imarker(jel) = ibset(0,0)

              ! The fourth edge of element JEL is only marked if it
              ! is also marked from the adjacent element
              iel1 = rhadapt%p_IneighboursAtElement(4, jel)
              iel2 = rhadapt%p_ImidneighboursAtElement(4, jel)
              if (ismarked_edge(iel1, iel2, jel)) p_Imarker(jel) = ibset(p_Imarker(jel),4)

              ! The first edge of element JEL is only marked if it
              ! is also marked from the adjacent element
              iel1 = rhadapt%p_IneighboursAtElement(1, jel)
              iel2 = rhadapt%p_ImidneighboursAtElement(1, jel)
              if (ismarked_edge(iel1, iel2, jel)) p_Imarker(jel) = ibset(p_Imarker(jel),1)
              

              ! The first, second and third edge of element LEL must be unmarked.
              p_Imarker(lel) = ibset(0,0)

              ! The fourth edge of element LEL is only marked if it
              ! is also marked from the adjacent element
              iel1 = rhadapt%p_IneighboursAtElement(4, lel)
              iel2 = rhadapt%p_ImidneighboursAtElement(4, lel)
              if (ismarked_edge(iel1, iel2, lel)) p_Imarker(lel) = ibset(p_Imarker(lel),4)
              
            end if

            ! If we are in the first iteration which is caused by elements
            ! marked for refinement due to accuracy reasons, then "lock"
            ! the four midpoint vertices which must not be deleted.
            if (bisFirst) then
              ivt = rhadapt%p_IverticesAtElement(2, iel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ivt = rhadapt%p_IverticesAtElement(3, iel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))

              ivt = rhadapt%p_IverticesAtElement(2, jel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ivt = rhadapt%p_IverticesAtElement(2, kel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ivt = rhadapt%p_IverticesAtElement(2, lel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
            end if
            
          end if   ! jstate = STATE_TRIA_OUTERINNER

                    
        case(STATE_QUAD_HALF1,&
             STATE_QUAD_HALF2)
          ! Element is quadrilateral that results from a Quad2Quad
          ! refinement. Due to the orientation convention the other
          ! "halved" element is adjacent to the second edge.
          jel = rhadapt%p_IneighboursAtElement(2, iel)

          ! Mark the fourth edge of elements IEL for subdivision
          iel1 = rhadapt%p_IneighboursAtElement(4, iel)
          iel2 = rhadapt%p_ImidneighboursAtElement(4, iel)
          if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(iel, iel1)

          ! Mark the fourth edge of elements JEL for subdivision
          iel1 = rhadapt%p_IneighboursAtElement(4, jel)
          iel2 = rhadapt%p_ImidneighboursAtElement(4, jel)
          if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(jel, iel1)
          
          ! Physically convert the two quadrilaterals into four quadrilaterals
          call convert_Quad2Quad(rhadapt, iel, jel,&
                                 rcollection, fcb_hadaptCallback)

          ! Reduce the number of green elements by two
          rhadapt%nGreenElements = rhadapt%nGreenElements-2

          ! We modified some elements in this iteration
          bisModified = .true.


          ! The new elements NEL0+1 and NEL0+2 have zero markers by construction.
          kel= rhadapt%p_IneighboursAtElement(2, jel)
          lel= rhadapt%p_IneighboursAtElement(3, iel)
          
          ! As a first step, clear all four elements and mark them as quadrilaterals.
          p_Imarker(iel) = ibset(0,0)
          p_Imarker(jel) = ibset(0,0)
          p_Imarker(kel) = ibset(0,0)
          p_Imarker(lel) = ibset(0,0)
          
          ! In addition, we have to transfer some markers from element IEL and JEL
          ! to the new elements NEL0+1, NEL0+2.
          iel1 = rhadapt%p_IneighboursAtElement(1, iel)
          iel2 = rhadapt%p_ImidneighboursAtElement(1, iel)
          if (ismarked_edge(iel1, iel2, iel)) p_Imarker(iel) = ibset(p_Imarker(iel),1)
          
          iel1 = rhadapt%p_IneighboursAtElement(4, jel)
          iel2 = rhadapt%p_ImidneighboursAtElement(4, jel)
          if (ismarked_edge(iel1, iel2, jel)) p_Imarker(jel) = ibset(p_Imarker(jel),4)
          
          iel1 = rhadapt%p_IneighboursAtElement(1, kel)
          iel2 = rhadapt%p_ImidneighboursAtElement(1, kel)
          if (ismarked_edge(iel1, iel2, kel)) p_Imarker(kel) = ibset(p_Imarker(kel),1)
          
          iel1 = rhadapt%p_IneighboursAtElement(4, lel)
          iel2 = rhadapt%p_ImidneighboursAtElement(4, lel)
          if (ismarked_edge(iel1, iel2, lel)) p_Imarker(lel) = ibset(p_Imarker(lel),4)


          ! If we are in the first iteration which is caused by elements
          ! marked for refinement due to accuracy reasons, then "lock"
          ! the four midpoint vertices which must not be deleted.
          if (bisFirst) then
            ivt = rhadapt%p_IverticesAtElement(2, iel)
            rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
            
            ivt = rhadapt%p_IverticesAtElement(3, iel)
            rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))

            ivt = rhadapt%p_IverticesAtElement(2, jel)
            rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
            
            ivt = rhadapt%p_IverticesAtElement(2, kel)
            rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
            
            ivt = rhadapt%p_IverticesAtElement(2, lel)
            rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
          end if

         
        case(STATE_TRIA_GREENINNER)
          ! We are processing a green triangle. Due to our refinement
          ! convention, element IEL can only be the inner triangle
          ! resulting from a Quad3Tria refinement.
          jel = rhadapt%p_IneighboursAtElement(3, iel)
          kel = rhadapt%p_IneighboursAtElement(1, iel)

          ! We need to mark the edges of the adjacent elements for subdivision.
          iel1 = rhadapt%p_IneighboursAtElement(2, iel)
          iel2 = rhadapt%p_ImidneighboursAtElement(2, iel)
          if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(iel, iel1)
          
          ! Mark the edge of the element adjacent to JEL for subdivision
          iel1 = rhadapt%p_IneighboursAtElement(3, jel)
          iel2 = rhadapt%p_ImidneighboursAtElement(3, jel)
          if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(jel, iel1)

          ! Mark the edge of the element adjacent to KEL for subdivision
          iel1 = rhadapt%p_IneighboursAtElement(1, kel)
          iel2 = rhadapt%p_ImidneighboursAtElement(1, kel)
          if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(kel, iel1)

          ! Physically convert the three elements IEL, JEL and KEL
          ! into four similar quadrilaterals
          call convert_Quad3Tria(rhadapt, jel, kel, iel,&
                                 rcollection, fcb_hadaptCallback)

          ! Reduce the number of green elements by three
          rhadapt%nGreenElements = rhadapt%nGreenElements-3

          ! We modified some elements in this iteration
          bisModified = .true.
          

          ! The new element NEL0+1 has zero marker by construction
          ! but that of IEL must be nullified. The markers for the
          ! modified triangles also need to be adjusted.
          p_Imarker(iel) = ibset(0,0)
          p_Imarker(jel) = ibset(0,0)
          p_Imarker(kel) = ibset(0,0)

          ! The first edge is only marked if it is also marked from the adjacent element
          iel1 = rhadapt%p_IneighboursAtElement(1, jel)
          iel2 = rhadapt%p_ImidneighboursAtElement(1, jel)
          if (ismarked_edge(iel1, iel2, jel)) p_Imarker(jel) = ibset(p_Imarker(jel),1)
          
          ! The fourth edge is only marked if it is also marked from the adjacent element
          iel1 = rhadapt%p_IneighboursAtElement(4, kel)
          iel2 = rhadapt%p_ImidneighboursAtElement(4, kel)
          if (ismarked_edge(iel1, iel2, kel)) p_Imarker(kel) = ibset(p_Imarker(kel),4)

          
          ! If we are in the first iteration which is caused by elements
          ! marked for refinement due to accuracy reasons, then "lock"
          ! the four midpoint vertices which must not be deleted.
          if (bisFirst) then
            ivt = rhadapt%p_IverticesAtElement(2, iel)
            rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))

            ivt = rhadapt%p_IverticesAtElement(3, iel)
            rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
            
            ivt = rhadapt%p_IverticesAtElement(4, iel)
            rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))

            ivt = rhadapt%p_IverticesAtElement(2, jel)
            rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
            
            ivt = rhadapt%p_IverticesAtElement(4, jel)
            rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
          end if
          
          
        case(STATE_TRIA_GREENOUTER_LEFT)
          ! We are processing a green triangle. Here, we have to consider several
          ! cases. Thus, find out the state of the neighboring element JEL.
          jel    = rhadapt%p_IneighboursAtElement(2, iel)
          jstate = redgreen_getState2D(rhadapt, jel)

          ! What state is element JEL
          select case(jstate)
          case(STATE_TRIA_GREENOUTER_RIGHT)
            ! Element IEL and JEL are the result of a Tria2Tria refinement,
            ! whereby triangle IEL is located left to element JEL.
            ! We can safely convert both elements into the macro element
            ! and perform regular refinement afterwards.

            ! We need to mark the edges of the elements adjacent
            ! to IEL and JEL for subdivision.
            iel1 = rhadapt%p_IneighboursAtElement(3, iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(3, iel)
            if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(iel, iel1)
            
            ! The same procedure must be applied to the neighbor of element JEL
            iel1 = rhadapt%p_IneighboursAtElement(1, jel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1, jel)
            if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(jel, iel1)
            
            ! Physically convert the two elements IEL and JEL into four similar triangles
            call convert_Tria2Tria(rhadapt, iel, jel,&
                                   rcollection, fcb_hadaptCallback)

            ! Reduce the number of green elements by two
            rhadapt%nGreenElements = rhadapt%nGreenElements-2

            ! We modified some elements in this iteration
            bisModified = .true.
         
   
            ! The new elements NEL0+1 and NEL0+2 have zero markers by construction.
            ! The markers for the modified elements IEL and JEL need to be adjusted.
            p_Imarker(iel) = 0
            p_Imarker(jel) = 0

            iel1 = rhadapt%p_IneighboursAtElement(1, iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1, iel)
            if (ismarked_edge(iel1, iel2, iel)) p_Imarker(iel) = ibset(p_Imarker(iel),1)
            
            iel1 = rhadapt%p_IneighboursAtElement(3, jel)
            iel2 = rhadapt%p_ImidneighboursAtElement(3, jel)
            if (ismarked_edge(iel1, iel2, jel)) p_Imarker(jel) = ibset(p_Imarker(jel),3)

            ! If we are in the first iteration which is caused by elements
            ! marked for refinement due to accuracy reasons, then "lock"
            ! the four midpoint vertices which must not be deleted.
            if (bisFirst) then
              ivt = rhadapt%p_IverticesAtElement(2, iel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))

              ivt = rhadapt%p_IverticesAtElement(3, iel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))

              ivt = rhadapt%p_IverticesAtElement(2, jel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
            end if

          
          case(STATE_TRIA_GREENINNER)
            ! Element IEL and JEL are the result of a Quad3Tria refinement, whereby
            ! element JEL is the inner triangle and IEL is its left neighbor.

            ! We need to mark the edges of the elements adjacent to IEL,
            ! JEL and KEL for subdivision. Let us start with the third neighbor of IEL.
            iel1 = rhadapt%p_IneighboursAtElement(3, iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(3, iel)
            if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(iel, iel1)

            ! The same procedure must be applied to the second neighbor of element JEL
            iel1 = rhadapt%p_IneighboursAtElement(2, jel)
            iel2 = rhadapt%p_ImidneighboursAtElement(2, jel)
            if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(jel, iel1)
            
            ! And again, the same procedure must be applied to the first neighbor of
            ! element KEL which is the first neighbor of JEL, whereby KEL needs to
            ! be found in the dynamic data structure first
            kel = rhadapt%p_IneighboursAtElement(1, jel)
            
            ! Ok, now we can proceed to its first neighbor
            iel1 = rhadapt%p_IneighboursAtElement(1, kel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1, kel)
            if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(kel, iel1)
            
            ! Physically convert the three elements IEL,JEL and KEL
            ! into four similar quadrilaterals
            call convert_Quad3Tria(rhadapt, iel, kel, jel,&
                                   rcollection, fcb_hadaptCallback)

            ! Reduce the number of green elements by three
            rhadapt%nGreenElements = rhadapt%nGreenElements-3

            ! We modified some elements in this iteration
            bisModified = .true.
            

            ! The new element NEL0+1 has zero marker by construction
            ! but that of JEL must be nullified. The markers for the
            ! modified triangles also need to be adjusted.
            p_Imarker(iel) = ibset(0,0)
            p_Imarker(jel) = ibset(0,0)
            p_Imarker(kel) = ibset(0,0)

            ! The first edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(1, iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1, iel)
            if (ismarked_edge(iel1, iel2, iel)) p_Imarker(iel) = ibset(p_Imarker(iel),1)

            ! The fourth edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(4, kel)
            iel2 = rhadapt%p_ImidneighboursAtElement(4, kel)
            if (ismarked_edge(iel1, iel2, kel)) p_Imarker(kel) = ibset(p_Imarker(kel),4)

            ! If we are in the first iteration which is caused by elements
            ! marked for refinement due to accuracy reasons, then "lock"
            ! the four midpoint vertices which must not be deleted.
            if (bisFirst) then
              ivt = rhadapt%p_IverticesAtElement(2, iel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ivt = rhadapt%p_IverticesAtElement(4, iel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ivt = rhadapt%p_IverticesAtElement(2, jel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ivt = rhadapt%p_IverticesAtElement(3, jel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ivt = rhadapt%p_IverticesAtElement(4, jel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
            end if

           
          case(STATE_TRIA_OUTERINNER)
            ! Element IEL and JEL are the result of a Quad4Tria refinement, whereby
            ! element JEL is the inner triangles  and IEL is its right neighbor.
            kel = rhadapt%p_IneighboursAtElement(2, jel)
            lel = rhadapt%p_IneighboursAtElement(3, jel)

            ! Mark the edge of the element adjacent to IEL for subdivision
            iel1 = rhadapt%p_IneighboursAtElement(3, iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(3, iel)
            if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(iel, iel1)
            
            ! Mark the edge of the element adjacent to LEL for subdivision
            iel1 = rhadapt%p_IneighboursAtElement(1, lel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1, lel)
            if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(lel, iel1)

            ! Physically convert the four triangles into four quadrilaterals
            call convert_Quad4Tria(rhadapt, iel, kel, lel, jel,&
                                   rcollection, fcb_hadaptCallback)

            ! Reduce the number of green elements by four
            rhadapt%nGreenElements = rhadapt%nGreenElements-4
            
            ! We modified some elements in this iteration
            bisModified = .true.

            
            ! All four elements have to be converted from triangles to quadrilaterals.
            p_Imarker(iel) = ibset(0,0)
            p_Imarker(jel) = ibset(0,0)
            p_Imarker(kel) = ibset(0,0)
            p_Imarker(lel) = ibset(0,0)
            
            ! The first edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(1, iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1, iel)
            if (ismarked_edge(iel1, iel2, iel)) p_Imarker(iel) = ibset(p_Imarker(iel),1)

            ! The fourth edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(4, kel)
            iel2 = rhadapt%p_ImidneighboursAtElement(4, kel)
            if (ismarked_edge(iel1, iel2, kel)) p_Imarker(kel) = ibset(p_Imarker(kel),4)

            ! The first edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(1, kel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1, kel)
            if (ismarked_edge(iel1, iel2, kel)) p_Imarker(kel) = ibset(p_Imarker(kel),1)

            ! The fourth edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(4, lel)
            iel2 = rhadapt%p_ImidneighboursAtElement(4, lel)
            if (ismarked_edge(iel1, iel2, lel)) p_Imarker(lel) = ibset(p_Imarker(lel),4)

            
            ! If we are in the first iteration which is caused by elements
            ! marked for refinement due to accuracy reasons, then "lock"
            ! the four midpoint vertices which must not be deleted.
            if (bisFirst) then
              ivt = rhadapt%p_IverticesAtElement(2, iel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ivt = rhadapt%p_IverticesAtElement(3, iel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ivt = rhadapt%p_IverticesAtElement(2, jel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ivt = rhadapt%p_IverticesAtElement(2, kel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ivt = rhadapt%p_IverticesAtElement(2, lel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
            end if


          case DEFAULT
            call output_line('Invalid element state!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'hadapt_markRedgreenRefinement2D')
            call sys_halt()
          end select


        case(STATE_TRIA_GREENOUTER_RIGHT)
          ! We are processing a green triangle. Here, we have to consider several cases.
          ! First, let us find out the state of the neighboring element JEL.
          jel    = rhadapt%p_IneighboursAtElement(2, iel)
          jstate = redgreen_getState2D(rhadapt, jel)
          
          ! What state is element JEL
          select case(jstate)
          case(STATE_TRIA_GREENOUTER_LEFT)
            ! Element IEL and JEL are the result of a Tria2Tria refinement, whereby
            ! triangle IEL is located right to element JEL. We can safely convert both
            ! elements into one and perform regular refinement afterwards.
            
            ! We need to mark the edges of the elements
            ! adjacent to IEL and JEL for subdivision.
            iel1 = rhadapt%p_IneighboursAtElement(1, iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1, iel)
            if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(iel, iel1)
            
            ! The same procedure must be applied to the neighbor of element JEL
            iel1 = rhadapt%p_IneighboursAtElement(3, jel)
            iel2 = rhadapt%p_ImidneighboursAtElement(3, jel)
            if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(jel, iel1)
            
            ! Physically convert the two elements IEL and JEL into four similar triangles
            call convert_Tria2Tria(rhadapt, jel, iel,&
                                   rcollection, fcb_hadaptCallback)

            ! Reduce the number of green elements by two
            rhadapt%nGreenElements = rhadapt%nGreenElements-2
            
            ! We modified some elements in this iteration
            bisModified = .true.

            
            ! The new elements NEL0+1 and NEL0+2 have zero markers by construction.
            ! The markers for the modified elements IEL and JEL need to be adjusted.
            p_Imarker(iel) = 0
            p_Imarker(jel) = 0
            
            iel1 = rhadapt%p_IneighboursAtElement(1, jel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1, jel)
            if (ismarked_edge(iel1, iel2, jel)) p_Imarker(jel) = ibset(p_Imarker(jel),1)
            
            iel1 = rhadapt%p_IneighboursAtElement(3, iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(3, iel)
            if (ismarked_edge(iel1, iel2, iel)) p_Imarker(iel) = ibset(p_Imarker(iel),3)
            
            ! If we are in the first iteration which is caused by elements
            ! marked for refinement due to accuracy reasons, then "lock"
            ! the four midpoint vertices which must not be deleted.
            if (bisFirst) then
              ivt = rhadapt%p_IverticesAtElement(2, iel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))

              ivt = rhadapt%p_IverticesAtElement(3, iel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ivt = rhadapt%p_IverticesAtElement(3, jel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
            end if


          case(STATE_TRIA_GREENINNER)
            ! Element IEL and JEL are the result of a Quad3Tria refinement, whereby
            ! element JEL is the inner triangle and IEL is its left neighbor.

            ! We need to mark the edges of the elements adjacent to IEL, JEL and
            ! KEL for subdivision. Let us start with the firth neighbor of IEL.
            iel1 = rhadapt%p_IneighboursAtElement(1, iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1, iel)
            if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(iel, iel1)

            ! The same procedure must be applied to the second neighbor of element JEL
            iel1 = rhadapt%p_IneighboursAtElement(2, jel)
            iel2 = rhadapt%p_ImidneighboursAtElement(2, jel)
            if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(jel, iel1)

            ! And again, the same procedure must be applied to the third neighbor of
            ! element KEL which is the third neighbor of JEL.
            kel  = rhadapt%p_IneighboursAtElement(3, jel)
            iel1 = rhadapt%p_IneighboursAtElement(3, kel)
            iel2 = rhadapt%p_ImidneighboursAtElement(3, kel)
            if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(kel, iel1)

            ! Physically convert the three elements IEL,JEL and KEL
            ! into four similar quadrilaterals
            call convert_Quad3Tria(rhadapt, kel, iel, jel,&
                                   rcollection, fcb_hadaptCallback)

            ! Reduce the number of green elements by three
            rhadapt%nGreenElements = rhadapt%nGreenElements-3
            
            ! We modified some elements in this iteration
            bisModified = .true.
            
            
            ! The new element NEL0+1 has zero marker by construction
            ! but that of JEL must be nullified. The markers for the
            ! modified triangles also need to be adjusted.
            p_Imarker(iel) = ibset(0,0)
            p_Imarker(jel) = ibset(0,0)
            p_Imarker(kel) = ibset(0,0)

            ! The first edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(1, kel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1, kel)
            if (ismarked_edge(iel1, iel2, kel)) p_Imarker(kel) = ibset(p_Imarker(kel),1)

            ! The fourth edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(4, iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(4, iel)
            if (ismarked_edge(iel1, iel2, iel)) p_Imarker(iel) = ibset(p_Imarker(iel),4)

            ! If we are in the first iteration which is caused by elements
            ! marked for refinement due to accuracy reasons, then "lock"
            ! the four midpoint vertices which must not be deleted.
            if (bisFirst) then
              ivt = rhadapt%p_IverticesAtElement(2, jel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ivt = rhadapt%p_IverticesAtElement(3, jel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ivt = rhadapt%p_IverticesAtElement(4, jel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ivt = rhadapt%p_IverticesAtElement(2, kel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ivt = rhadapt%p_IverticesAtElement(4, kel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
            end if

            
          case(STATE_TRIA_OUTERINNER)
            ! Element IEL and JEL are the result of a Quad4Tria refinement, whereby
            ! element JEL is the inner triangles  and IEL is its left neighbor.
            kel = rhadapt%p_IneighboursAtElement(2, jel)
            lel = rhadapt%p_IneighboursAtElement(1, jel)

            ! Mark the edge of the element adjacent to IEL for subdivision
            iel1 = rhadapt%p_IneighboursAtElement(1, iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1, iel)
            if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(iel, iel1)

            ! Mark the edge of the element adjacent to LEL for subdivision
            iel1 = rhadapt%p_IneighboursAtElement(3, lel)
            iel2 = rhadapt%p_ImidneighboursAtElement(3, lel)
            if (iel1*iel2 .ne. 0 .and. iel1 .eq. iel2) call mark_edge(lel, iel1)

            ! Physically convert the four triangles into four quadrilaterals
            call convert_Quad4Tria(rhadapt, lel, kel, iel, jel,&
                                   rcollection, fcb_hadaptCallback)

            ! Reduce the number of green elements by four
            rhadapt%nGreenElements = rhadapt%nGreenElements-4

            ! We modified some elements in this iteration
            bisModified = .true.
            

            ! All four elements have to be converted from triangles to quadrilaterals.
            p_Imarker(iel) = ibset(0,0)
            p_Imarker(jel) = ibset(0,0)
            p_Imarker(kel) = ibset(0,0)
            p_Imarker(lel) = ibset(0,0)

            ! The first edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(1, lel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1, lel)
            if (ismarked_edge(iel1, iel2, lel)) p_Imarker(lel) = ibset(p_Imarker(lel),1)

            ! The fourth edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(4, kel)
            iel2 = rhadapt%p_ImidneighboursAtElement(4, kel)
            if (ismarked_edge(iel1, iel2, kel)) p_Imarker(kel) = ibset(p_Imarker(kel),4)

            ! The first edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(1, kel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1, kel)
            if (ismarked_edge(iel1, iel2, kel)) p_Imarker(kel) = ibset(p_Imarker(kel),1)

            ! The fourth edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(4, iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(4, iel)
            if (ismarked_edge(iel1, iel2, iel)) p_Imarker(iel) = ibset(p_Imarker(iel),4)

            ! If we are in the first iteration which is caused by elements
            ! marked for refinement due to accuracy reasons, then "lock"
            ! the four midpoint vertices which must not be deleted.
            if (bisFirst) then
              ivt = rhadapt%p_IverticesAtElement(2, iel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ivt = rhadapt%p_IverticesAtElement(3, iel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ivt = rhadapt%p_IverticesAtElement(2, jel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ivt = rhadapt%p_IverticesAtElement(2, kel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
              
              ivt = rhadapt%p_IverticesAtElement(2, lel)
              rhadapt%p_IvertexAge(ivt) = -abs(rhadapt%p_IvertexAge(ivt))
            end if


          case DEFAULT
            call output_line('Invalid element state!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'hadapt_markRedgreenRefinement2D')
            call sys_halt()
          end select

          
        case DEFAULT
          call output_line('Invalid element state!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'hadapt_markRedgreenRefinement2D')
          call sys_halt()
        end select
100     continue
        
        
        !-----------------------------------------------------------------------
        ! Loop over all adjacent cells of element IEL
        !-----------------------------------------------------------------------
        adjacent: do ive = 1, nve

          ! If element IEL is adjacent to two different elements along one edge
          ! then ignore this edge. This situation can arise, if the element IEL
          ! is adjacent to a green element JEL which has been converted previously,
          ! so that element IEL has two neighbors along the edge, i.e., IEL and one
          ! new element number >NEL0
          if (rhadapt%p_IneighboursAtElement(ive, iel) .ne. &
              rhadapt%p_ImidneighboursAtElement(ive, iel)) cycle adjacent
          
          ! If the edge shared by element IEL and its neighbor JEL has been
          ! marked for refinement and JEL is not outside of the domain, i.e.,
          ! element IEL is located at the boundary, then proceed to element JEL.
          if (btest(p_Imarker(iel), ive)) then
            
            ! Check if the current edge is located at the
            ! boundary, then nothing needs to be done
            jel = rhadapt%p_IneighboursAtElement(ive, iel)
            if (jel .eq. 0) cycle adjacent
            
            ! Get number of vertices per element
            mve = hadapt_getNVE(rhadapt, jel)
            
            ! Are we triangular or quadrilateral element?
            select case(mve)
            case(TRIA_NVETRI2D)
              p_Imarker(jel) = ibclr(p_Imarker(jel),0)

            case(TRIA_NVEQUAD2D)
              p_Imarker(jel) = ibset(p_Imarker(jel),0)

            case DEFAULT
              call output_line('Invalid number of vertices per element!',&
                               OU_CLASS_ERROR,OU_MODE_STD,'hadapt_markRedgreenRefinement2D')
              call sys_halt()
            end select

            ! Now, we need to find the local position of element IEL in the
            ! adjacency list of the nieghboring element JEL
            do jve = 1, mve
              if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. iel) exit
            end do

            if (jve .gt.mve) then
              call output_line('Unable to find element!',&
                               OU_CLASS_ERROR,OU_MODE_STD,'hadapt_markRedgreenRefinement2D')
              call sys_halt()
            end if

            ! If the edge is already marked for refinement then we can
            ! guarantee conformity for the edge shared by IEL and JEL.
            ! Otherwise, this edge needs to be marked "looking" from JEL
            if (.not.btest(p_Imarker(jel), jve)) then
              p_Imarker(jel) = ibset(p_Imarker(jel), jve)

              ! This element should not be considered in this iteration.
              ! It is therefore pre-marked for the next iteration.
              if (p_Imodifier(jel) .ne. imodifier) p_Imodifier(jel)= -imodifier

              ! We modified some elements in this iteration
              bisModified = .true.
            end if

          end if
        end do adjacent
      end do element
      
      
      ! Note that the number of elements (and vertices) may have
      ! changed due to conversion of green into red elements.
      
      ! Hence, adjust dimensions.
      rhadapt%NEL0 = rhadapt%NEL
      rhadapt%NVT0 = rhadapt%NVT

      !-------------------------------------------------------------------------
      ! We do not want to have elements with two and/or three divided edges,
      ! aka, blue refinement for triangles and quadrilaterals, respectively.
      ! Therefore, these elements are filtered and marked for red refinement.
      !-------------------------------------------------------------------------
      noblue: do iel = 1, rhadapt%NEL0
        if (p_Imodifier(iel)*imodifier .eq. 0) cycle noblue
        
        ! What type of element are we?
        select case(p_Imarker(iel))
        case(MARK_REF_TRIA3TRIA_12,&
             MARK_REF_TRIA3TRIA_13,&
             MARK_REF_TRIA3TRIA_23)
          ! Blue refinement for triangles is not allowed.
          ! Hence, mark triangle for red refinement
          p_Imarker(iel)      =  MARK_REF_TRIA4TRIA
          p_Imodifier(iel)    = -imodifier
          bisModified         =  .true.
          rhadapt%increaseNVT =  rhadapt%increaseNVT+1

          
        case(MARK_REF_QUADBLUE_123,&
             MARK_REF_QUADBLUE_412,&
             MARK_REF_QUADBLUE_341,&
             MARK_REF_QUADBLUE_234)
          ! Blue refinement for quadrilaterals is not allowed.
          ! Hence, mark quadrilateral for red refinement
          p_Imarker(iel)      =  MARK_REF_QUAD4QUAD
          p_Imodifier(iel)    = -imodifier
          bisModified         =  .true.
          rhadapt%increaseNVT =  rhadapt%increaseNVT+2


        case DEFAULT
          p_Imodifier(iel) = (p_Imodifier(iel)-imodifier)/2
        end select
      end do noblue
      
      ! Reverse modifier
      imodifier = -imodifier

      ! Once we arrived at this point, the first iteration is done
      bisFirst = .false.
    end do conform


    ! Loop over all elements and increase the number of vertices according
    ! to the number of quadrilaterals which are marked for red refinement.
    do iel = 1, size(p_Imarker, 1)
      if (p_Imarker(iel) .eq. MARK_REF_QUAD4QUAD) then
        rhadapt%increaseNVT = rhadapt%increaseNVT+1
      end if
    end do
    
    ! Free auxiliary storage
    call storage_free(h_Imodifier)

  contains

    ! Here, the real working routines follow.

    !**************************************************************
    ! For a given element IEL, mark the edge that connects IEL
    ! to its neighbor JEL in the marker array at position JEL
    
    subroutine mark_edge(iel, jel)

      integer, intent(in) :: iel
      integer, intent(in) :: jel

      ! local variables
      integer :: ive,nve

      ! Get number of vertices per element
      nve = hadapt_getNVE(rhadapt, jel)
      
      ! Find local position of element IEL in adjacency list of JEL
      select case(nve)
      case(TRIA_NVETRI2D)
        ! Triangular elements
        do ive = 1, nve
          if (rhadapt%p_IneighboursAtElement(ive, jel) .eq. iel) then
            p_Imarker(jel) = ibclr(ibset(p_Imarker(jel), ive),0)
            bisModified    =.true.
            if (p_Imodifier(jel) .ne. imodifier) p_Imodifier(jel) = -imodifier
            return
          end if
        end do

      case(TRIA_NVEQUAD2D)
        ! Quadrilateral element
        do ive = 1, nve
          if (rhadapt%p_IneighboursAtElement(ive, jel) .eq. iel) then
            p_Imarker(jel) = ibset(ibset(p_Imarker(jel), ive),0)
            bisModified    =.true.
            if (p_Imodifier(jel) .ne. imodifier) p_Imodifier(jel) = -imodifier
            return
          end if
        end do

      case DEFAULT
        call output_line('Invalid number of vertices per element!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'mark_edge')
        call sys_halt()
      end select

    end subroutine mark_edge

    !**************************************************************
    ! For a given element IEL, check if the edge that connects IEL
    ! with its adjacent element JEL is marked

    function ismarked_edge(iel, ielmid, jel) result (bismarked)

      integer, intent(in) :: iel,ielmid,jel
      logical :: bismarked

      ! local variables
      integer :: ive

      ! Are we at the boundary?
      if (iel*ielmid .eq. 0) then
        bismarked = .false.
        return
      elseif(iel .ne. ielmid) then
        bismarked = .true.
        return
      end if
      
      ! Loop over all edges of element IEL and try to find neighbor JEL
      do ive = 1, hadapt_getNVE(rhadapt, iel)
        if (rhadapt%p_IneighboursAtElement(ive, iel)    .eq. jel .or.&
            rhadapt%p_ImidneighboursAtElement(ive, iel) .eq. jel) then
          bismarked = btest(p_Imarker(iel), ive)
          return
        end if
      end do
      
      call output_line('Unable to find common egde!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'ismarked_edge')
      call sys_halt()

    end function ismarked_edge

  end subroutine hadapt_markRedgreenRefinement2D

  ! ***************************************************************************

!<function>

  function hadapt_CalcNumberOfElements2D(rhadapt) result(nel)

!<description>
    ! This function calculates the number of elements present
    ! in the triangulation after refinement has been performed
!</description>

!<input>
    ! Adaptivity structure
    type(t_hadapt), intent(in) :: rhadapt
!</input>

!<result>
    ! Number of elements after refinement
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
      case(MARK_REF_TRIA2TRIA_1,&
           MARK_REF_TRIA2TRIA_2,&
           MARK_REF_TRIA2TRIA_3,&
           MARK_REF_QUAD2QUAD_13,&
           MARK_REF_QUAD2QUAD_24,&
           MARK_CRS_2QUAD3TRIA)
        ! Interestingly enought, the coarsening of 2 quadrilaterals into three
        ! triangles which reduces the number of vertices by one also increases
        ! the number of elements by one which has to be taken into account here.
        nel = nel+1
        
      case(MARK_REF_QUAD3TRIA_1,&
           MARK_REF_QUAD3TRIA_2,&
           MARK_REF_QUAD3TRIA_3,&
           MARK_REF_QUAD3TRIA_4,&
           MARK_REF_TRIA3TRIA_12,&
           MARK_REF_TRIA3TRIA_23,&
           MARK_REF_TRIA3TRIA_13)
        nel = nel+2
        
      case(MARK_REF_TRIA4TRIA,&
           MARK_REF_QUAD4QUAD,&
           MARK_REF_QUAD4TRIA_12,&
           MARK_REF_QUAD4TRIA_23,&
           MARK_REF_QUAD4TRIA_34,&
           MARK_REF_QUAD4TRIA_14)
        nel = nel+3
      end select
    end do
  end function hadapt_CalcNumberOfElements2D

  ! ***************************************************************************
  
!<subroutine>

  subroutine hadapt_refine2D(rhadapt, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine refinement the elements according to the marker
!</description>

!<input>
    ! Callback function
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
    integer, dimension(:), pointer :: p_Imarker
    integer :: iel

    ! Check if dynamic data structures are o.k. and if
    ! cells are marked for refinement
    if (iand(rhadapt%iSpec, HADAPT_HAS_DYNAMICDATA2D) .ne.&
                            HADAPT_HAS_DYNAMICDATA2D .or.&
        iand(rhadapt%iSpec, HADAPT_MARKEDREFINE) .ne.&
                            HADAPT_MARKEDREFINE) then
      call output_line('Dynamic data structures are not generated ' // &
                       'or no marker for refinement is available!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_refine2D')
      call sys_halt()
    end if
    
    ! Set pointers
    call storage_getbase_int(rhadapt%h_Imarker, p_Imarker)

    ! Perform red-green refinement
    do iel = 1, size(p_Imarker, 1)
      
      select case(p_Imarker(iel))
      case(MARK_ASIS_TRIA,&
           MARK_ASIS_QUAD)

        ! Do nothing for elements that should be kept 'as is'

      case(:MARK_CRS_GENERIC)

        ! Do nothing for element that have been marked for coarsening

      case(MARK_REF_TRIA4TRIA)

        ! Red refinement triangle
        call refine_Tria4Tria(rhadapt, iel,&
                              rcollection, fcb_hadaptCallback)
        p_Imarker(iel) = MARK_ASIS

      case(MARK_REF_QUAD4QUAD)

        ! Red refinement quadrilateral
        call refine_Quad4Quad(rhadapt, iel,&
                              rcollection, fcb_hadaptCallback)
        p_Imarker(iel) = MARK_ASIS
        
      case(MARK_REF_TRIA2TRIA_1,&
           MARK_REF_TRIA2TRIA_2,&
           MARK_REF_TRIA2TRIA_3)

        ! Green refinement triangle
        call refine_Tria2Tria(rhadapt, iel, p_Imarker(iel),&
                              rcollection, fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements+2
        p_Imarker(iel)         = MARK_ASIS
        
      case(MARK_REF_QUAD3TRIA_1,&
           MARK_REF_QUAD3TRIA_2,&
           MARK_REF_QUAD3TRIA_3,&
           MARK_REF_QUAD3TRIA_4)

        ! Green refinement quadrilateral
        call refine_Quad3Tria(rhadapt, iel, p_Imarker(iel),&
                              rcollection, fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements+3
        p_Imarker(iel)         = MARK_ASIS
        
      case(MARK_REF_QUAD2QUAD_13,&
           MARK_REF_QUAD2QUAD_24)

        ! Green refinement quadrilateral
        call refine_Quad2Quad(rhadapt, iel, p_Imarker(iel),&
                              rcollection, fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements+2
        p_Imarker(iel)         = MARK_ASIS

      case(MARK_REF_QUAD4TRIA_12,&
           MARK_REF_QUAD4TRIA_23,&
           MARK_REF_QUAD4TRIA_34,&
           MARK_REF_QUAD4TRIA_14)

        ! Green refinement quadrilateral
        call refine_Quad4Tria(rhadapt, iel, p_Imarker(iel),&
                              rcollection, fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements+4
        p_Imarker(iel)         = MARK_ASIS

      case DEFAULT
        call output_line('Invalid element refinement marker!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_refine2D')
        call sys_halt()
      end select
    end do

    ! Increase the number of refinement steps by one
    rhadapt%nRefinementSteps = rhadapt%nRefinementSteps+1
    
    ! The markers are no longer valid
    rhadapt%iSpec = iand(rhadapt%iSpec, not(HADAPT_MARKEDREFINE))

  end subroutine hadapt_refine2D

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_coarsen2D(rhadapt, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine coarsens the elements according to the marker
!</description>

!<input>
    ! Callback function
    include 'intf_hadaptcallback.inc'
    optional :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
    
    ! OPTIONAL Collection
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer,  dimension(:), pointer :: p_Imarker
    integer :: ipos,iel,jel,ivt,ivtReplace,ive

    ! Check if dynamic data structures are o.k. and
    ! if  cells are marked for coarsening
    if (iand(rhadapt%iSpec, HADAPT_HAS_DYNAMICDATA2D) .ne.&
                            HADAPT_HAS_DYNAMICDATA2D .or.&
        iand(rhadapt%iSpec, HADAPT_MARKEDCOARSEN) .ne.&
                            HADAPT_MARKEDCOARSEN) then
      call output_line('Dynamic data structures are not generated ' // &
                       'or no marker for coarsening is available!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_coarsen2D')
      call sys_halt()
    end if

    ! Set pointers
    call storage_getbase_int(rhadapt%h_Imarker, p_Imarker)

    
    ! Perform hierarchical red-green recoarsening
    element: do iel = size(p_Imarker,1), 1, -1

      select case(p_Imarker(iel))
      case(MARK_CRS_GENERIC:)

        ! Do nothing for elements that should be kept 'as is'
        ! and those which are marked for refinement.
        
      case(MARK_CRS_2TRIA1TRIA)
        call coarsen_2Tria1Tria(rhadapt, p_Imarker, iel,&
                                rcollection, fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements-2
        p_Imarker(iel)         = MARK_ASIS
        
      case(MARK_CRS_4TRIA1TRIA)
        call coarsen_4Tria1Tria(rhadapt, p_Imarker, iel,&
                                rcollection, fcb_hadaptCallback)
        p_Imarker(iel) = MARK_ASIS

      case(MARK_CRS_4TRIA2TRIA_1,&
           MARK_CRS_4TRIA2TRIA_2,&
           MARK_CRS_4TRIA2TRIA_3)
        call coarsen_4Tria2Tria(rhadapt, p_Imarker, iel,&
                                rcollection, fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements+2
        p_Imarker(iel)         = MARK_ASIS

      case(MARK_CRS_3TRIA1QUAD)
        call coarsen_3Tria1Quad(rhadapt, p_Imarker, iel,&
                                rcollection, fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements-3
        p_Imarker(iel)         = MARK_ASIS

      case(MARK_CRS_4TRIA1QUAD)
        call coarsen_4Tria1Quad(rhadapt, p_Imarker, iel,&
                                rcollection, fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements-4
        p_Imarker(iel)         = MARK_ASIS
        
      case(MARK_CRS_4TRIA3TRIA_LEFT,&
           MARK_CRS_4TRIA3TRIA_RIGHT)
        call coarsen_4Tria3Tria(rhadapt, p_Imarker, iel,&
                                rcollection, fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements-1
        p_Imarker(iel)         = MARK_ASIS

      case(MARK_CRS_2QUAD1QUAD)
        call coarsen_2Quad1Quad(rhadapt, p_Imarker, iel,&
                                rcollection, fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements-2
        p_Imarker(iel)         = MARK_ASIS

      case(MARK_CRS_2QUAD3TRIA)
        call coarsen_2Quad3Tria(rhadapt, p_Imarker, iel,&
                                rcollection, fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements+1
        p_Imarker(iel)         = MARK_ASIS
        
      case(MARK_CRS_4QUAD1QUAD)
        call coarsen_4Quad1Quad(rhadapt, p_Imarker, iel,&
                                rcollection, fcb_hadaptCallback)
        p_Imarker(iel) = MARK_ASIS

      case(MARK_CRS_4QUAD2QUAD)
        call coarsen_4Quad2Quad(rhadapt, p_Imarker, iel,&
                                rcollection, fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements+2
        p_Imarker(iel)         = MARK_ASIS

      case(MARK_CRS_4QUAD3TRIA)
        call coarsen_4Quad3Tria(rhadapt, p_Imarker, iel,&
                                rcollection, fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements+3
        p_Imarker(iel)         = MARK_ASIS
        
      case(MARK_CRS_4QUAD4TRIA)
        call coarsen_4Quad4Tria(rhadapt, p_Imarker, iel,&
                                rcollection, fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements+4
        p_Imarker(iel)         = MARK_ASIS

      case DEFAULT
        call output_line('Invalid recoarsening marker!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_coarsen2D')
        call sys_halt()
      end select

    end do element


    ! Loop over all vertices 1...NVT0 present in the triangulation before
    ! refinement and check if they are free for vertex removal.
    vertex: do ivt = rhadapt%NVT0, 1, -1
      
      ! If the vertex is locked, then skip this vertex
      if (rhadapt%p_IvertexAge(ivt) .le. 0) cycle vertex
      
      ! Remove vertex physically. Note that this vertex is no longer associated
      ! to any element. All associations have been removed in the above element
      ! coarsening/conversion step. In order to prevent "holes" in the vertex list,
      ! vertex IVT is replaced by the last vertex if it is not the last one itself.
      call remove_vertex2D(rhadapt, ivt, ivtReplace)
      
      ! If vertex IVT was not the last one, update the "elements-meeting-at-vertex" list
      if (ivtReplace .ne. 0) then
        
        ! Start with first element in "elements-meeting-at-vertex" list of the replaced vertex
        ipos = arrlst_getNextInArraylist(rhadapt%rElementsAtVertex,&
                                         ivtReplace, .true.)

        update: do while(ipos .gt. ARRLST_NULL)
          
          ! Get element number JEL
          jel = rhadapt%rElementsAtVertex%p_IData(ipos)
          
          ! Proceed to next element
          ipos = arrlst_getNextInArraylist(rhadapt%rElementsAtVertex,&
                                           ivtReplace, .false.)
          
          ! Look for vertex ivtReplace in element JEL and replace it by IVT
          do ive = 1, hadapt_getNVE(rhadapt, jel)
            if (rhadapt%p_IverticesAtElement(ive, jel) .eq. ivtReplace) then
              rhadapt%p_IverticesAtElement(ive, jel) = ivt
              cycle update
            end if
          end do
          
          ! If the replaced vertex ivtReplace could not be found in element JEL
          ! something is wrong and we stop the simulation
          call output_line('Unable to find replacement vertex in element',&
                           OU_CLASS_ERROR,OU_MODE_STD,'hadapt_coarsen2D')
          call sys_halt()
        end do update
        
        ! Swap tables IVT and ivtReplace in arraylist and release table ivtReplace
        call arrlst_swapArrayList(rhadapt%rElementsAtVertex, ivt, ivtReplace)
        call arrlst_releaseArrayList(rhadapt%rElementsAtVertex, ivtReplace)
        
      else
        
        ! Release table IVT
        call arrlst_releaseArrayList(rhadapt%rElementsAtVertex, ivt)
      end if
            
      ! Optionally, invoke callback function
      if (present(fcb_hadaptCallback) .and. present(rcollection)) then
        rcollection%IquickAccess(1:2) = (/ivt,ivtReplace/)
        call fcb_hadaptCallback(HADAPT_OPR_REMOVEVERTEX, rcollection)
      end if
    end do vertex
        
    ! Increase the number of recoarsening steps by one
    rhadapt%nCoarseningSteps = rhadapt%nCoarseningSteps+1

    ! The markers are no longer valid
    rhadapt%iSpec = iand(rhadapt%iSpec, not(HADAPT_MARKEDCOARSEN))

  end subroutine hadapt_coarsen2D

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_writeGridSVG2D(rhadapt, coutputFile, width, height)

!<description>
    ! This subroutine outputs the current state of the adapted grid stored
    ! in the dynamic data structure to a given file in SVG format in 2D
!</description>

!<input>
    ! Output file name w/o suffix .svg
    character(LEN=*), intent(in) :: coutputFile

    ! OPTIONAL: Width of the generated file
    integer, intent(in), optional :: width

    ! OPTIONAL: Heigh of the generated file
    integer, intent(in), optional :: height
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local parameters
    integer, parameter :: defaultwidth  = 1280
    integer, parameter :: defaultheight = 1024
    integer, parameter :: xoffset       = 100
    integer, parameter :: yoffset       = 100
    integer, parameter :: font_size     = 18
    integer, parameter :: colsep        = 10
    integer, parameter :: linesep       = 32
    character(LEN=*), parameter :: font_family = "Arial"
    
    ! local variables
    real(DP), dimension(2*NDIM2D) :: bbox
    real(DP) :: x0,y0,xdim,ydim,dscale
    
    integer, dimension(:), pointer :: p_Imarker
    integer :: ivt,iel,ipos,xsize,ysize,iunit,ive,nve
    integer, save :: iout=0
    
    ! Check if dynamic data structures generated
    if (iand(rhadapt%iSpec, HADAPT_HAS_DYNAMICDATA2D) .ne.&
                            HADAPT_HAS_DYNAMICDATA2D) then
      call output_line('Dynamic data structures are not generated',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_writeGridSVG2D')
      call sys_halt()
    end if
    
    ! Increment the sample number
    iout=iout+1
    
    ! Open output file for writing
    call io_openFileForWriting(trim(adjustl(coutputFile))//'.'//&
        trim(sys_siL(iout,5))//'.svg', iunit, SYS_REPLACE, bformatted=.true.)
    
    ! Write prolog, XML-declaration, document type declaration)
    write(iunit,FMT='(A)') '<?xml version="1.0" encoding="utf-8" standalone="yes"?>'
    write(iunit,FMT='(A)') '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN"'
    write(iunit,FMT='(A)') ' "http://www.w3.org/Graphics/SVG/1.0/DTD/svg11.dtd">'
    write(iunit,FMT='(A)')
    write(iunit,FMT='(A)') '<!-- Created with Featflow2 (http://www.featflow.de/) -->'
    write(iunit,FMT='(A)')
!<!--
    ! Hide this block from the automatic documentation parser to prevent
    ! strings like "<polygon" being regarded as XML tags:
    write(iunit,FMT='(A)') '<svg version="1.0"'
    write(iunit,FMT='(A)') ' xmlns="http://www.w3.org/2000/svg"'
    write(iunit,FMT='(A)') ' xmlns:xlink="http://www.w3.org/1999/xlink"'

 
    if (present(width)) then
      xsize=width
    else
      xsize=defaultwidth
    end if

    if (present(height)) then
      ysize=height
    else
      ysize=defaultheight
    end if

    ! Determine bounding box
    bbox   = qtree_getBoundingBox(rhadapt%rVertexCoordinates2D)
    xdim   = bbox(3)-bbox(1); x0=bbox(1)
    ydim   = bbox(4)-bbox(2); y0=bbox(2)
    dscale = min(xsize/xdim,ysize/ydim)

    ! Set height and width of image
    write(iunit,FMT='(A)') ' width="100%" height="100%" xml:space="preserve"'
    write(iunit,FMT='(A)') ' viewBox="'//&
        trim(sys_siL(-xoffset,10))//' '//&
        trim(sys_siL(-yoffset,10))//' '//&
        trim(sys_siL(2*xoffset+xsize,10))//' '//&
        trim(sys_siL(2*yoffset+ysize,10))//'">'

    ! Write embedded java-scripts to file
    write(iunit,FMT='(A)') '<defs>'
    write(iunit,FMT='(A)') '<script type="text/javascript">'
    write(iunit,FMT='(A)') '<![CDATA['

    write(iunit,FMT='(A)') '  function ShowVertexInfo(evt,ivt,age,elems) {'
    write(iunit,FMT='(A)') '    var svgdoc=document.documentElement;'
    write(iunit,FMT='(A)') '    var x=evt.clientX; var y=evt.clientY;'
    write(iunit,FMT='(A)') '    var m=svgdoc.getScreenCTM();'
    write(iunit,FMT='(A)') '    var p=svgdoc.createSVGPoint();'
    write(iunit,FMT='(A)') '    p.x=evt.clientX; p.y=evt.clientY;'
    write(iunit,FMT='(A)') '    p=p.matrixTransform(m.inverse());'
    write(iunit,FMT='(A)') '    p.x=p.x+15; p.y=p.y+15;'

    write(iunit,FMT='(A)') '    var vinfo=svgdoc.getElementById("vinfo");'
    write(iunit,FMT='(A)') '    var vinfo_box=svgdoc.getElementById("vinfo_box");'
    write(iunit,FMT='(A)') '    var vinfo_text1=svgdoc.getElementById("vinfo_text1");'
    write(iunit,FMT='(A)') '    var vinfo_text2=svgdoc.getElementById("vinfo_text2");'
    write(iunit,FMT='(A)') '    var vinfo_text3=svgdoc.getElementById("vinfo_text3");'
    write(iunit,FMT='(A)') '    var vinfo_text4=svgdoc.getElementById("vinfo_text4");'
    write(iunit,FMT='(A)') '    var vinfo_text5=svgdoc.getElementById("vinfo_text5");'
    write(iunit,FMT='(A)') '    var vinfo_text6=svgdoc.getElementById("vinfo_text6");'
    
    write(iunit,FMT='(A)') '    vinfo_box.setAttribute("x",p.x);'
    write(iunit,FMT='(A)') '    vinfo_box.setAttribute("y",p.y);'
    write(iunit,FMT='(A)') '    vinfo_text1.setAttribute("x",p.x+'//trim(sys_siL(colsep,10))//');'
    write(iunit,FMT='(A)') '    vinfo_text1.setAttribute("y",p.y+'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    vinfo_text2.setAttribute("x",p.x+'//trim(sys_siL(colsep,10))//');'
    write(iunit,FMT='(A)') '    vinfo_text2.setAttribute("y",p.y+2*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    vinfo_text3.setAttribute("x",p.x+'//trim(sys_siL(colsep,10))//');'
    write(iunit,FMT='(A)') '    vinfo_text3.setAttribute("y",p.y+3*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    vinfo_text4.setAttribute("x",p.x+vinfo_text1.getComputedTextLength()+20);'
    write(iunit,FMT='(A)') '    vinfo_text4.setAttribute("y",p.y+'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    vinfo_text4.firstChild.nodeValue = ivt;'
    write(iunit,FMT='(A)') '    vinfo_text5.setAttribute("x",p.x+vinfo_text2.getComputedTextLength()+20);'
    write(iunit,FMT='(A)') '    vinfo_text5.setAttribute("y",p.y+2*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    vinfo_text5.firstChild.nodeValue = age;'
    write(iunit,FMT='(A)') '    vinfo_text6.setAttribute("x",p.x+vinfo_text3.getComputedTextLength()+20);'
    write(iunit,FMT='(A)') '    vinfo_text6.setAttribute("y",p.y+3*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    vinfo_text6.firstChild.nodeValue = elems;'

    write(iunit,FMT='(A)') '    var textlen = vinfo_text1.getComputedTextLength()+' // &
        'vinfo_text4.getComputedTextLength();'
    write(iunit,FMT='(A)') '    textlen = Math.max(textlen,vinfo_text2.getComputedTextLength()' // &
        '+vinfo_text5.getComputedTextLength());'
    write(iunit,FMT='(A)') '    textlen = Math.max(textlen,vinfo_text3.getComputedTextLength()' // &
        '+vinfo_text6.getComputedTextLength());'
    write(iunit,FMT='(A)') '    vinfo_box.setAttribute("width",textlen+30);'

    write(iunit,FMT='(A)') '    vinfo.setAttribute("style","visibility: visible");'
    write(iunit,FMT='(A)') '  }'
    
    write(iunit,FMT='(A)') '  function HideVertexInfo() {'
    write(iunit,FMT='(A)') '    var svgdoc=document.documentElement;'
    write(iunit,FMT='(A)') '    var vinfo=svgdoc.getElementById("vinfo");'
    write(iunit,FMT='(A)') '    vinfo.setAttribute("style","visibility: hidden");'
    write(iunit,FMT='(A)') '  }'

    write(iunit,FMT='(A)') '  function ShowElementInfo(evt,iel,state,marker,kadj,kmidadj,kvert) {'
    write(iunit,FMT='(A)') '    var svgdoc=document.documentElement;'
    write(iunit,FMT='(A)') '    var x=evt.clientX; var y=evt.clientY;'
    write(iunit,FMT='(A)') '    var m=svgdoc.getScreenCTM();'
    write(iunit,FMT='(A)') '    var p=svgdoc.createSVGPoint();'
    write(iunit,FMT='(A)') '    p.x=evt.clientX; p.y=evt.clientY;'
    write(iunit,FMT='(A)') '    p=p.matrixTransform(m.inverse());'
    write(iunit,FMT='(A)') '    p.x=p.x+15; p.y=p.y+15;'

    write(iunit,FMT='(A)') '    var einfo=svgdoc.getElementById("einfo");'
    write(iunit,FMT='(A)') '    var einfo_box=svgdoc.getElementById("einfo_box");'
    write(iunit,FMT='(A)') '    var einfo_text1=svgdoc.getElementById("einfo_text1");'
    write(iunit,FMT='(A)') '    var einfo_text2=svgdoc.getElementById("einfo_text2");'
    write(iunit,FMT='(A)') '    var einfo_text3=svgdoc.getElementById("einfo_text3");'
    write(iunit,FMT='(A)') '    var einfo_text4=svgdoc.getElementById("einfo_text4");'
    write(iunit,FMT='(A)') '    var einfo_text5=svgdoc.getElementById("einfo_text5");'
    write(iunit,FMT='(A)') '    var einfo_text6=svgdoc.getElementById("einfo_text6");'
    write(iunit,FMT='(A)') '    var einfo_text7=svgdoc.getElementById("einfo_text7");'
    write(iunit,FMT='(A)') '    var einfo_text8=svgdoc.getElementById("einfo_text8");'
    write(iunit,FMT='(A)') '    var einfo_text9=svgdoc.getElementById("einfo_text9");'
    write(iunit,FMT='(A)') '    var einfo_text10=svgdoc.getElementById("einfo_text10");'
    write(iunit,FMT='(A)') '    var einfo_text11=svgdoc.getElementById("einfo_text11");'
    write(iunit,FMT='(A)') '    var einfo_text12=svgdoc.getElementById("einfo_text12");'
    
    write(iunit,FMT='(A)') '    einfo_box.setAttribute("x",p.x);'
    write(iunit,FMT='(A)') '    einfo_box.setAttribute("y",p.y);'
    write(iunit,FMT='(A)') '    einfo_text1.setAttribute("x",p.x+'//trim(sys_siL(colsep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text1.setAttribute("y",p.y+'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text2.setAttribute("x",p.x+'//trim(sys_siL(colsep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text2.setAttribute("y",p.y+2*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text3.setAttribute("x",p.x+'//trim(sys_siL(colsep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text3.setAttribute("y",p.y+3*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text4.setAttribute("x",p.x+'//trim(sys_siL(colsep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text4.setAttribute("y",p.y+4*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text5.setAttribute("x",p.x+'//trim(sys_siL(colsep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text5.setAttribute("y",p.y+5*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text6.setAttribute("x",p.x+'//trim(sys_siL(colsep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text6.setAttribute("y",p.y+6*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text7.setAttribute("x",p.x+einfo_text1.getComputedTextLength()+20);'
    write(iunit,FMT='(A)') '    einfo_text7.setAttribute("y",p.y+'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text7.firstChild.nodeValue = iel;'
    write(iunit,FMT='(A)') '    einfo_text8.setAttribute("x",p.x+einfo_text2.getComputedTextLength()+20);'
    write(iunit,FMT='(A)') '    einfo_text8.setAttribute("y",p.y+2*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text8.firstChild.nodeValue = state;'
    write(iunit,FMT='(A)') '    einfo_text9.setAttribute("x",p.x+einfo_text3.getComputedTextLength()+20);'
    write(iunit,FMT='(A)') '    einfo_text9.setAttribute("y",p.y+3*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text9.firstChild.nodeValue = marker;'
    write(iunit,FMT='(A)') '    einfo_text10.setAttribute("x",p.x+einfo_text4.getComputedTextLength()+20);'
    write(iunit,FMT='(A)') '    einfo_text10.setAttribute("y",p.y+4*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text10.firstChild.nodeValue = kadj;'
    write(iunit,FMT='(A)') '    einfo_text11.setAttribute("x",p.x+einfo_text5.getComputedTextLength()+20);'
    write(iunit,FMT='(A)') '    einfo_text11.setAttribute("y",p.y+5*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text11.firstChild.nodeValue = kmidadj;'
    write(iunit,FMT='(A)') '    einfo_text12.setAttribute("x",p.x+einfo_text6.getComputedTextLength()+20);'
    write(iunit,FMT='(A)') '    einfo_text12.setAttribute("y",p.y+6*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text12.firstChild.nodeValue = kvert;'

    write(iunit,FMT='(A)') '    var textlen = einfo_text1.getComputedTextLength()+' // &
        'einfo_text7.getComputedTextLength();'
    write(iunit,FMT='(A)') '    textlen = Math.max(textlen,einfo_text2.getComputedTextLength()+' // &
        'einfo_text8.getComputedTextLength());'
    write(iunit,FMT='(A)') '    textlen = Math.max(textlen,einfo_text3.getComputedTextLength()+' // &
        'einfo_text9.getComputedTextLength());'
    write(iunit,FMT='(A)') '    textlen = Math.max(textlen,einfo_text4.getComputedTextLength()+' // &
        'einfo_text10.getComputedTextLength());'
    write(iunit,FMT='(A)') '    textlen = Math.max(textlen,einfo_text5.getComputedTextLength()+' // &
        'einfo_text11.getComputedTextLength());'
    write(iunit,FMT='(A)') '    textlen = Math.max(textlen,einfo_text6.getComputedTextLength()+' // &
        'einfo_text12.getComputedTextLength());'
    write(iunit,FMT='(A)') '    einfo_box.setAttribute("width",textlen+30);'

    write(iunit,FMT='(A)') '    einfo.setAttribute("style","visibility: visible");'
    write(iunit,FMT='(A)') '  }'

    write(iunit,FMT='(A)') '  function HideElementInfo() {'
    write(iunit,FMT='(A)') '    var svgdoc=document.documentElement;'
    write(iunit,FMT='(A)') '    var einfo=svgdoc.getElementById("einfo");'
    write(iunit,FMT='(A)') '    einfo.setAttribute("style","visibility: hidden");'
    write(iunit,FMT='(A)') '  }'
    write(iunit,FMT='(A)') ']]>'
    write(iunit,FMT='(A)') '</script>'
    write(iunit,FMT='(A)') '</defs>'

    !---------------------------------------------------------------------------
    ! Output all elements
    !---------------------------------------------------------------------------
    if (iand(rhadapt%iSpec, HADAPT_MARKEDREFINE)  .eq. HADAPT_MARKEDREFINE .or.&
        iand(rhadapt%iSpec, HADAPT_MARKEDCOARSEN) .eq. HADAPT_MARKEDCOARSEN) then

      ! Set pointer to marker
      call storage_getbase_int(rhadapt%h_Imarker, p_Imarker)

      ! Output elements and color those which are marked for refinement
      do iel = 1, rhadapt%NEL
        
        ! Get number of vertices per elements
        nve = hadapt_getNVE(rhadapt, iel)
        
        if (iel .le. rhadapt%NEL0) then
          
          ! What kind of element are we?
          select case(p_Imarker(iel))
            
          ! Element is neither marked for refinement nor coarsening
          case(MARK_ASIS_TRIA,&
               MARK_ASIS_QUAD)
            write(iunit,FMT='(A)') '<polygon id="el'//trim(sys_siL(iel,9))//&
                '" fill="white" stroke="black" stroke-width="1"'

          ! Element is marked for green refinement
          case(MARK_REF_TRIA2TRIA_1,&
               MARK_REF_TRIA2TRIA_2,&
               MARK_REF_TRIA2TRIA_3,&
               MARK_REF_QUAD3TRIA_1,&
               MARK_REF_QUAD3TRIA_2,&
               MARK_REF_QUAD3TRIA_3,&
               MARK_REF_QUAD3TRIA_4,&
               MARK_REF_QUAD4TRIA_12,&
               MARK_REF_QUAD4TRIA_23,&
               MARK_REF_QUAD4TRIA_34,&
               MARK_REF_QUAD4TRIA_14,&
               MARK_REF_QUAD2QUAD_13,&
               MARK_REF_QUAD2QUAD_24)
            write(iunit,FMT='(A)') '<polygon id="el'//trim(sys_siL(iel,9))//&
                '" fill="green" stroke="black" stroke-width="1"'

          ! Element is marked for blue refinement
          case(MARK_REF_TRIA3TRIA_12,&
               MARK_REF_TRIA3TRIA_23,&
               MARK_REF_TRIA3TRIA_13,&
               MARK_REF_QUADBLUE_412,&
               MARK_REF_QUADBLUE_234,&
               MARK_REF_QUADBLUE_123,&
               MARK_REF_QUADBLUE_341)
            write(iunit,FMT='(A)') '<polygon id="el'//trim(sys_siL(iel,9))//&
                '" fill="blue" stroke="black" stroke-width="1"'

          ! Element is marked for red refinement
          case(MARK_REF_TRIA4TRIA,&
               MARK_REF_QUAD4QUAD)
            write(iunit,FMT='(A)') '<polygon id="el'//trim(sys_siL(iel,9))//&
                '" fill="red" stroke="black" stroke-width="1"'

          ! Element is marked for generic coarsening
          case(MARK_CRS_GENERIC)
            write(iunit,FMT='(A)') '<polygon id="el'//trim(sys_siL(iel,9))//&
                '" fill="lightgray" stroke="black" stroke-width="1"'
            
          ! Element is marked for specific coarsening
          case(MARK_CRS_4QUAD4TRIA:MARK_CRS_4TRIA1TRIA)
            write(iunit,FMT='(A)') '<polygon id="el'//trim(sys_siL(iel,9))//&
                '" fill="yellow" stroke="black" stroke-width="1"'

          ! Unknown element marker
          case DEFAULT
            write(iunit,FMT='(A)') '<polygon id="el'//trim(sys_siL(iel,9))//&
                '" fill="hotpink" stroke="black" stroke-width="1"'
          end select
          
          ! Write data which is common to all elements
          select case(nve)
          case(TRIA_NVETRI2D)
            write(iunit,FMT='(A)') ' onmousemove="ShowElementInfo(evt,'''//&
                trim(sys_siL(iel,10))//''','''//&
                trim(sys_siL(redgreen_getState2D(rhadapt,iel),4))//''','''//&
                trim(sys_siL(p_Imarker(iel),5))//''','''//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(3,iel),10))//','//&
                trim(sys_siL(0,10))//''','''//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(3,iel),10))//','//&
                trim(sys_siL(0,10))//''','''//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(3,iel),10))//','//&
                trim(sys_siL(0,10))//''')"'
            write(iunit,FMT='(A)') ' onmouseout="HideElementInfo()" style="cursor: help"'

          case(TRIA_NVEQUAD2D)
            write(iunit,FMT='(A)') ' onmousemove="ShowElementInfo(evt,'''//&
                trim(sys_siL(iel,10))//''','''//&
                trim(sys_siL(redgreen_getState2D(rhadapt,iel),4))//''','''//&
                trim(sys_siL(p_Imarker(iel),5))//''','''//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(3,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(4,iel),10))//''','''//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(3,iel),10))//','//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(4,iel),10))//''','''//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(3,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(4,iel),10))//''')"'
            write(iunit,FMT='(A)') ' onmouseout="HideElementInfo()" style="cursor: help"'
          end select

        else
          
          ! For all new elements there si no marker available
          write(iunit,FMT='(A)') '<polygon id="el'//trim(sys_siL(iel,9))//&
              '" fill="white" stroke="black" stroke-width="1"'
          
          ! Write data which is common to all elements
          select case(nve)
          case(TRIA_NVETRI2D)
            write(iunit,FMT='(A)') ' onmousemove="ShowElementInfo(evt,'''//&
                trim(sys_siL(iel,10))//''','''//&
                trim(sys_siL(redgreen_getState2D(rhadapt,iel),4))//''','''//&
                trim(sys_siL(0,5))//''','''//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(3,iel),10))//','//&
                trim(sys_siL(0,10))//''','''//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(3,iel),10))//','//&
                trim(sys_siL(0,10))//''','''//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(3,iel),10))//','//&
                trim(sys_siL(0,10))//''')"'
            write(iunit,FMT='(A)') ' onmouseout="HideElementInfo()" style="cursor: help"'

          case(TRIA_NVEQUAD2D)
            write(iunit,FMT='(A)') ' onmousemove="ShowElementInfo(evt,'''//&
                trim(sys_siL(iel,10))//''','''//&
                trim(sys_siL(redgreen_getState2D(rhadapt,iel),4))//''','''//&
                trim(sys_siL(0,5))//''','''//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(3,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(4,iel),10))//''','''//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(3,iel),10))//','//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(4,iel),10))//''','''//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(3,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(4,iel),10))//''')"'
            write(iunit,FMT='(A)') ' onmouseout="HideElementInfo()" style="cursor: help"'
          end select
          
        end if
               
        ! Each element is a polygon made up from 3/4 points
        write(iunit,FMT='(A)',ADVANCE='NO') ' points="'
        do ive =1 , nve
          xdim = qtree_getX(rhadapt%rVertexCoordinates2D,&
                            rhadapt%p_IverticesAtElement(ive,iel))-x0
          ydim = qtree_getY(rhadapt%rVertexCoordinates2D,&
                            rhadapt%p_IverticesAtElement(ive,iel))-y0
          write(iunit,FMT='(A)',ADVANCE='NO') &
              trim(sys_siL(int(dscale*xdim),10))//','//&
              trim(sys_siL(ysize-int(dscale*ydim),10))//' '
        end do
        
        ! And of course, the polygone must be closed !
        xdim = qtree_getX(rhadapt%rVertexCoordinates2D,&
                          rhadapt%p_IverticesAtElement(1,iel))-x0
        ydim = qtree_getY(rhadapt%rVertexCoordinates2D,&
                          rhadapt%p_IverticesAtElement(1,iel))-y0
        write(iunit,FMT='(A)') &
            trim(sys_siL(int(dscale*xdim),10))//','//&
            trim(sys_siL(ysize-int(dscale*ydim),10))//'"/>'

      end do
      
    else
      
      ! Output only elements without individual coloring
      do iel = 1, rhadapt%NEL

        ! Get number of vertices per element
        nve = hadapt_getNVE(rhadapt, iel)
        
        ! For all new elements there si no marker available
        write(iunit,FMT='(A)') '<polygon id="el'//trim(sys_siL(iel,9))//&
            '" fill="white" stroke="black" stroke-width="1"'

        ! Write data which is common to all elements
        select case(nve)
        case(TRIA_NVETRI2D)
          write(iunit,FMT='(A)') ' onmousemove="ShowElementInfo(evt,'''//&
              trim(sys_siL(iel,10))//''','''//&
              trim(sys_siL(redgreen_getState2D(rhadapt,iel),4))//''','''//&
              trim(sys_siL(0,5))//''','''//&
              trim(sys_siL(rhadapt%p_IneighboursAtElement(1,iel),10))//','//&
              trim(sys_siL(rhadapt%p_IneighboursAtElement(2,iel),10))//','//&
              trim(sys_siL(rhadapt%p_IneighboursAtElement(3,iel),10))//','//&
              trim(sys_siL(0,10))//''','''//&
              trim(sys_siL(rhadapt%p_ImidneighboursAtElement(1,iel),10))//','//&
              trim(sys_siL(rhadapt%p_ImidneighboursAtElement(2,iel),10))//','//&
              trim(sys_siL(rhadapt%p_ImidneighboursAtElement(3,iel),10))//','//&
              trim(sys_siL(0,10))//''','''//&
              trim(sys_siL(rhadapt%p_IverticesAtElement(1,iel),10))//','//&
              trim(sys_siL(rhadapt%p_IverticesAtElement(2,iel),10))//','//&
              trim(sys_siL(rhadapt%p_IverticesAtElement(3,iel),10))//','//&
              trim(sys_siL(0,10))//''')"'
          write(iunit,FMT='(A)') ' onmouseout="HideElementInfo()" style="cursor: help"'
          
        case(TRIA_NVEQUAD2D)
          write(iunit,FMT='(A)') ' onmousemove="ShowElementInfo(evt,'''//&
              trim(sys_siL(iel,10))//''','''//&
              trim(sys_siL(redgreen_getState2D(rhadapt,iel),4))//''','''//&
              trim(sys_siL(0,5))//''','''//&
              trim(sys_siL(rhadapt%p_IneighboursAtElement(1,iel),10))//','//&
              trim(sys_siL(rhadapt%p_IneighboursAtElement(2,iel),10))//','//&
              trim(sys_siL(rhadapt%p_IneighboursAtElement(3,iel),10))//','//&
              trim(sys_siL(rhadapt%p_IneighboursAtElement(4,iel),10))//''','''//&
              trim(sys_siL(rhadapt%p_ImidneighboursAtElement(1,iel),10))//','//&
              trim(sys_siL(rhadapt%p_ImidneighboursAtElement(2,iel),10))//','//&
              trim(sys_siL(rhadapt%p_ImidneighboursAtElement(3,iel),10))//','//&
              trim(sys_siL(rhadapt%p_ImidneighboursAtElement(4,iel),10))//''','''//&
              trim(sys_siL(rhadapt%p_IverticesAtElement(1,iel),10))//','//&
              trim(sys_siL(rhadapt%p_IverticesAtElement(2,iel),10))//','//&
              trim(sys_siL(rhadapt%p_IverticesAtElement(3,iel),10))//','//&
              trim(sys_siL(rhadapt%p_IverticesAtElement(4,iel),10))//''')"'
          write(iunit,FMT='(A)') ' onmouseout="HideElementInfo()" style="cursor: help"'
        end select
        
        ! Each element is a polygon made up from 3/4 points
        write(iunit,FMT='(A)',ADVANCE='NO') ' points="'
        do ive = 1, nve
          xdim = qtree_getX(rhadapt%rVertexCoordinates2D,&
                            rhadapt%p_IverticesAtElement(ive,iel))-x0
          ydim = qtree_getY(rhadapt%rVertexCoordinates2D,&
                            rhadapt%p_IverticesAtElement(ive,iel))-y0
          write(iunit,FMT='(A)',ADVANCE='NO') &
              trim(sys_siL(int(dscale*xdim),10))//','//&
              trim(sys_siL(ysize-int(dscale*ydim),10))//' '
        end do
        
        ! And of course, the polygone must be closed !
        xdim = qtree_getX(rhadapt%rVertexCoordinates2D,&
                          rhadapt%p_IverticesAtElement(1,iel))-x0
        ydim = qtree_getY(rhadapt%rVertexCoordinates2D,&
                          rhadapt%p_IverticesAtElement(1,iel))-y0
        write(iunit,FMT='(A)') &
            trim(sys_siL(int(dscale*xdim),10))//','//&
            trim(sys_siL(ysize-int(dscale*ydim),10))//'"/>'
        
      end do
    end if
         
    ! Loop over all vertices
    do ivt = 1, qtree_getsize(rhadapt%rVertexCoordinates2D)
      xdim = qtree_getX(rhadapt%rVertexCoordinates2D,ivt)-x0
      ydim = qtree_getY(rhadapt%rVertexCoordinates2D,ivt)-y0
      
      ! Write vertices as points?
      if (rhadapt%p_IvertexAge(ivt) .gt.0) then
        write(iunit,FMT='(A)') '<circle id="vt'//trim(sys_siL(ivt,10))//'" cx="'//&
            trim(sys_siL(int(dscale*xdim),10))//'" cy="'//&
            trim(sys_siL(ysize-int(dscale*ydim),10))//&
            '" r="0pt" fill="white" stroke="black" stroke-width="1pt"'
      else
        write(iunit,FMT='(A)') '<circle id="vt'//trim(sys_siL(ivt,10))//'" cx="'//&
            trim(sys_siL(int(dscale*xdim),10))//'" cy="'//&
            trim(sys_siL(ysize-int(dscale*ydim),10))//&
            '" r="0pt" fill="black" stroke="black" stroke-width="1pt"'
      end if
      
      ! Write data which is common to all vertices
      write(iunit,FMT='(A)',ADVANCE='NO') ' onmousemove="ShowVertexInfo(evt,'''//&
          trim(sys_siL(ivt,10))//''','''//&
          trim(sys_siL(rhadapt%p_IvertexAge(ivt),5))//''','''
      
      ! Generate list of elements meeting at vertex
      ipos = arrlst_getNextInArrayList(rhadapt%rElementsAtVertex, ivt, .true.)
      do while(ipos .gt. ARRLST_NULL)
        ! Get element number IEL
        iel = rhadapt%rElementsAtVertex%p_IData(ipos)
        
        ! Proceed to next entry in array list
        ipos = arrlst_getNextInArraylist(rhadapt%rElementsAtVertex, ivt, .false.)

        if (ipos .gt. ARRLST_NULL) then
          write(iunit,FMT='(A)',ADVANCE='NO') trim(sys_siL(iel,10))//','
        else
          write(iunit,FMT='(A)',ADVANCE='NO') trim(sys_siL(iel,10))
        end if
      end do
      write(iunit,FMT='(A)') ''')"'
      write(iunit,FMT='(A)') ' onmouseout="HideVertexInfo()" style="cursor: crosshair"/>'
    end do


    ! Output info boxes
    write(iunit,FMT='(A)') '<g id="vinfo" style="visibility: hidden">'
    write(iunit,FMT='(A)') '  <rect id="vinfo_box"  x="0" y="0" width="200" height="110"'
    write(iunit,FMT='(A)') '        style="fill: #FFFFCC; stroke: #000000; stroke-width: 0.5px;"/>'
    write(iunit,FMT='(A)') '  <text id="vinfo_text1" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black">Vertex number:</text>'
    write(iunit,FMT='(A)') '  <text id="vinfo_text2" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black">Vertex age:</text>'
    write(iunit,FMT='(A)') '  <text id="vinfo_text3" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black">Elements meeting at vertex:</text>'
    write(iunit,FMT='(A)') '  <text id="vinfo_text4" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black">null</text>'
    write(iunit,FMT='(A)') '  <text id="vinfo_text5" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black">null</text>'
    write(iunit,FMT='(A)') '  <text id="vinfo_text6" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black">null</text>'
    write(iunit,FMT='(A)') '</g>'

    write(iunit,FMT='(A)') '<g id="einfo" style="visibility: hidden">'
    write(iunit,FMT='(A)') '  <rect id="einfo_box"  x="0" y="0" width="200" height="210"'
    write(iunit,FMT='(A)') '        style="fill: #FFFFCC; stroke: #000000; stroke-width: 0.5px;"/>'
    write(iunit,FMT='(A)') '  <text id="einfo_text1" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black">Element number:</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text2" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">Element state:</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text3" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">Element marker:</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text4" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">Element neighbours:</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text5" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">Element mid-neighbours:</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text6" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">Element vertices:</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text7" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">null</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text8" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">null</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text9" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">null</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text10" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">null</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text11" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">null</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text12" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">null</text>'
    write(iunit,FMT='(A)') '</g>'

    ! Close XML-file
    write(iunit,FMT='(A)') '</svg>'
    close(iunit)

    ! Unhide from here from the automatic documentation parser: -->

  end subroutine hadapt_writeGridSVG2D

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<subroutine>

  subroutine add_vertex_atEdgeMidpoint2D(rhadapt, i1, i2, e1, i12,&
                                         rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine adds a new vertex at the midpoint of a given egde.
    ! First, the coordinates of the new vertex are computed and added to
    ! the quadtree. If the new vertex is located at the boundary then its
    ! parametrization is computed and its normal vector is stored in
    ! the boundary data structure.
!</description>

!<input>
    ! First point of edge on which new vertex will be added
    integer, intent(in) :: i1

    ! Second point of edge on which new vertex will be added
    integer, intent(in) :: i2

    ! Number of the right-adjacent element w.r.t. to the oriented edge (I1,I2)
    integer, intent(in) :: e1

    ! Callback function
    include 'intf_hadaptcallback.inc'
    optional :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! Adaptivity structure
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
    real(DP), dimension(NDIM2D) :: Dcoord
    real(DP), dimension(1) :: Ddata
    real(DP) :: x1,y1,x2,y2,dvbdp1,dvbdp2
    integer, dimension(2) :: Idata
    integer :: inode,ipred,ipos,ibct
    
    ! Get coordinates of edge vertices
    x1 = qtree_getX(rhadapt%rVertexCoordinates2D, i1)
    y1 = qtree_getY(rhadapt%rVertexCoordinates2D, i1)
    x2 = qtree_getX(rhadapt%rVertexCoordinates2D, i2)
    y2 = qtree_getY(rhadapt%rVertexCoordinates2D, i2)

    ! Compute coordinates of new vertex
    Dcoord = 0.5_DP * (/x1+x2, y1+y2/)

    ! Search for vertex coordinates in quadtree: If the vertex already
    ! exists, e.g., it was added when the adjacent element was
    ! refined, then nothing needs to be done for this vertex

    if (qtree_searchInQuadtree(rhadapt%rVertexCoordinates2D,&
                               Dcoord, inode, ipos, i12) .eq. QTREE_NOT_FOUND) then
      
      ! Add new entry to vertex coordinates
      if (qtree_insertIntoQuadtree(rhadapt%rVertexCoordinates2D,&
                                   Dcoord, i12, inode) .ne. QTREE_INSERTED) then
        call output_line('An error occured while inserting coordinate to quadtree!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'add_vertex_atEdgeMidpoint2D')
        call sys_halt()
      end if

      ! Update number of vertices
      rhadapt%NVT = rhadapt%NVT+1

      if (i12 .ne. rhadapt%NVT) then
        call output_line('Inconsistent vertex number due to insertion!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'add_vertex_atEdgeMidpoint2D')
        call sys_halt()
      end if
      
      ! Set age of vertex
      rhadapt%p_IvertexAge(i12) = 1+max(abs(rhadapt%p_IvertexAge(i1)),&
                                        abs(rhadapt%p_IvertexAge(i2)))
      
      ! Are we at the boundary?
      if (e1 .eq. 0) then
        ! Set nodal property
        rhadapt%p_InodalProperty(i12) = rhadapt%p_InodalProperty(i1)

        ! Increment number of boundary nodes
        rhadapt%NVBD = rhadapt%NVBD+1
        
        ! Get number of boundary component
        ibct = rhadapt%p_InodalProperty(i1)
        
        ! Get parameter values of the boundary nodes
        if (btree_searchInTree(rhadapt%rBoundary(ibct), i1, ipred) .eq. BTREE_NOT_FOUND) then
          call output_line('Unable to find first vertex in boundary data structure!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'add_vertex_atEdgeMidpoint2D')
          call sys_halt()
        end if
        ipos   = rhadapt%rBoundary(ibct)%p_Kchild(merge(TLEFT, TRIGHT, ipred .lt. 0), abs(ipred))
        dvbdp1 = rhadapt%rBoundary(ibct)%p_DData(BdrValue, ipos)
        
        if (btree_searchInTree(rhadapt%rBoundary(ibct), i2, ipred) .eq. BTREE_NOT_FOUND) then
          call output_line('Unable to find second vertex in boundary data structure!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'add_vertex_atEdgeMidpoint2D')
          call sys_halt()
        end if
        ipos   = rhadapt%rBoundary(ibct)%p_Kchild(merge(TLEFT, TRIGHT, ipred .lt. 0), abs(ipred))
        dvbdp2 = rhadapt%rBoundary(ibct)%p_DData(BdrValue, ipos)
        
        ! If I2 is last(=first) node on boundary component IBCT round DVBDP2 to next integer
        if (dvbdp2 .le. dvbdp1) dvbdp2 = ceiling(dvbdp1)
        
        ! Add new entry to boundary structure
        Idata = (/i1, i2/)
        Ddata = (/0.5_DP*(dvbdp1+dvbdp2)/)
        call btree_insertIntoTree(rhadapt%rBoundary(ibct), i12,&
                                 Idata=Idata, Ddata=Ddata)
      else
        ! Set nodal property
        rhadapt%p_InodalProperty(i12) = 0
      end if

      ! Optionally, invoke callback function
      if (present(fcb_hadaptCallback) .and. present(rcollection)) then
        rcollection%IquickAccess(1:3) = (/i12, i1, i2/)
        call fcb_hadaptCallback(HADAPT_OPR_INSERTVERTEXEDGE, rcollection)
      end if

    end if

  end subroutine add_vertex_atEdgeMidpoint2D

  ! ***************************************************************************
  
!<subroutine>

  subroutine add_vertex_atElementCenter2D(rhadapt, i1, i2, i3, i4, i5,&
                                          rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine adds a new vertex at the center of a given quadtrilateral.
    ! First, the coordinates of the new vertex computed and added to the
    ! quadtree. The new vertex cannot be located at the boundary.
!</description>

!<input>
    ! Four corners of the quadrilateral
    integer, intent(in) :: i1,i2,i3,i4

    ! Callback function
    include 'intf_hadaptcallback.inc'
    optional :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt

    ! OPTIONAL: Collection
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Number of the new vertex
    integer, intent(out) :: i5
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(NDIM2D) :: Dcoord
    real(DP) :: x1,y1,x2,y2,x3,y3,x4,y4
    integer :: inode,ipos
    
    ! Compute coordinates of new vertex
    x1 = qtree_getX(rhadapt%rVertexCoordinates2D, i1)
    y1 = qtree_getY(rhadapt%rVertexCoordinates2D, i1)
    x2 = qtree_getX(rhadapt%rVertexCoordinates2D, i2)
    y2 = qtree_getY(rhadapt%rVertexCoordinates2D, i2)
    x3 = qtree_getX(rhadapt%rVertexCoordinates2D, i3)
    y3 = qtree_getY(rhadapt%rVertexCoordinates2D, i3)
    x4 = qtree_getX(rhadapt%rVertexCoordinates2D, i4)
    y4 = qtree_getY(rhadapt%rVertexCoordinates2D, i4)
    
    Dcoord = 0.25_DP * (/x1+x2+x3+x4, y1+y2+y3+y4/)

    ! Search for vertex coordinates in quadtree
    if (qtree_searchInQuadtree(rhadapt%rVertexCoordinates2D, Dcoord,&
                               inode, ipos, i5) .eq. QTREE_NOT_FOUND) then
      
      ! Add new entry to vertex coordinates
      if (qtree_insertIntoQuadtree(rhadapt%rVertexCoordinates2D,&
                                   Dcoord, i5, inode) .ne. QTREE_INSERTED) then
        call output_line('An error occured while inserting coordinate to quadtree!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'add_vertex_atElementCenter2D')
        call sys_halt()
      end if

      ! Update number of vertices
      rhadapt%NVT = rhadapt%NVT+1

      if (i5 .ne. rhadapt%NVT) then
        call output_line('Inconsistent vertex number during insertion!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'add_vertex_atElementCenter2D')
        call sys_halt()
      end if
      
      ! Set age of vertex
      rhadapt%p_IvertexAge(I5) = 1+max(abs(rhadapt%p_IvertexAge(i1)),&
                                       abs(rhadapt%p_IvertexAge(i2)),&
                                       abs(rhadapt%p_IvertexAge(i3)),&
                                       abs(rhadapt%p_IvertexAge(i4)))
      
      ! Set nodal property
      rhadapt%p_InodalProperty(i5) = 0
      
      ! Optionally, invoke callback function
      if (present(fcb_hadaptCallback) .and. present(rcollection)) then
        rcollection%IquickAccess(1:5) = (/i5, i1, i2, i3, i4/)
        call fcb_hadaptCallback(HADAPT_OPR_INSERTVERTEXCENTR, rcollection)
      end if

    end if

  end subroutine add_vertex_atElementCenter2D

  ! ***************************************************************************

!<subroutine>
  
  subroutine remove_vertex2D(rhadapt, ivt, ivtReplace)
  
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
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>

!<output>
    ! Number of the vertex to replace the deleted one
    integer, intent(out) :: ivtReplace
!</output>
!</subroutine>

    ! local variables
    integer :: i1,i2
    integer :: ipred,ipos
    integer :: ibct

    ! Remove vertex from coordinates and get number of replacement vertex
    if (qtree_deleteFromQuadtree(rhadapt%rVertexCoordinates2D,&
                                 ivt, ivtReplace) .ne. QTREE_DELETED) then
      call output_line('Unable to delete vertex coordinates!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'remove_vertex2D')
      call sys_halt()
    end if
    
    ! Decrease number of vertices by one
    rhadapt%NVT = rhadapt%NVT-1
    
    ! If IVT is a boundary node remove it from the boundary and
    ! connect its boundary neighbors with each other
    ibct = rhadapt%p_InodalProperty(ivt)
    
    ! Are we at the boundary?
    if (ibct .ne. 0) then
      
      ! Find position of vertex IVT in boundary array
      if (btree_searchInTree(rhadapt%rBoundary(ibct),&
                             ivt, ipred) .eq. BTREE_NOT_FOUND) then
        call output_line('Unable to find vertex in boundary data structure!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'remove_vertex2D')
        call sys_halt()
      end if
      ipos = rhadapt%rBoundary(ibct)%p_Kchild(merge(TLEFT,TRIGHT, ipred .lt. 0), abs(ipred))
      
      ! Get the two boundary neighbors: I1 <- IVT -> I2
      i1 = rhadapt%rBoundary(ibct)%p_IData(BdrPrev, ipos)
      i2 = rhadapt%rBoundary(ibct)%p_IData(BdrNext, ipos)
      
      ! Connect the boundary neighbors with each other: I1 <=> I2
      ! First, set I2 as next neighboring of I1
      if (btree_searchInTree(rhadapt%rBoundary(ibct),&
                             i1, ipred) .eq. BTREE_NOT_FOUND) then
        call output_line('Unable to find left neighboring vertex in boundary data structure!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'remove_vertex2D')
        call sys_halt()
      end if
      ipos = rhadapt%rBoundary(ibct)%p_Kchild(merge(TLEFT,TRIGHT, ipred .lt. 0), abs(ipred))
      rhadapt%rBoundary(ibct)%p_IData(BdrNext, ipos) = i2
      
      ! Second, set I1 as previous neighbor of I2
      if (btree_searchInTree(rhadapt%rBoundary(ibct),&
                             i2, ipred) .eq. BTREE_NOT_FOUND) then
        call output_line('Unable to find right neighboring vertex in boundary data structure!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'remove_vertex2D')
        call sys_halt()
      end if
      ipos = rhadapt%rBoundary(ibct)%p_Kchild(merge(TLEFT,TRIGHT, ipred .lt. 0), abs(ipred))
      rhadapt%rBoundary(ibct)%p_IData(BdrPrev, ipos) = i1
      
      ! And finally, delete IVT from the boundary
      if (btree_deleteFromTree(rhadapt%rBoundary(ibct),&
                               ivt) .eq. BTREE_NOT_FOUND) then
        call output_line('Unable to delete vertex from boundary data structure!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'remove_vertex2D')
        call sys_halt()
      end if

      ! Decrease number of boundary nodes
      rhadapt%NVBD = rhadapt%NVBD-1
    end if
    
    ! If IVT is not the last node then copy the data for vertex IVTREPLACE
    ! to IVT and prepare IVTREPLACE for elimination
    if (ivt .lt. ivtReplace) then
      
      ! If IVTREPLACE is a boundary node then remove IVTREPLACE from the boundary
      ! vector, insert IVT into the boundary vector instead, and connect the
      ! boundary neighbors of IVTREPLACE with IVT
      ibct = rhadapt%p_InodalProperty(ivtReplace)

      ! Are we at the boundary?
      if (ibct .ne. 0) then

        if (btree_searchInTree(rhadapt%rBoundary(ibct),&
                               ivtReplace, ipred) .eq. BTREE_NOT_FOUND) then
          call output_line('Unable to find replacement vertex in boundary data structure!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'remove_vertex2D')
          call sys_halt()
        end if
        ipos = rhadapt%rBoundary(ibct)%p_Kchild(merge(TLEFT,TRIGHT, ipred .lt. 0), abs(ipred))
        
        ! Insert IVT into the boundary vector
        call btree_insertIntoTree(rhadapt%rBoundary(ibct), ivt,&
                                  Idata=rhadapt%rBoundary(ibct)%p_IData(:, ipos),&
                                  Ddata=rhadapt%rBoundary(ibct)%p_DData(:, ipos))
        
        ! Get the two boundary neighbors: I1 <- IVTREPLACE -> I2
        i1 = rhadapt%rBoundary(ibct)%p_IData(BdrPrev, ipos)
        i2 = rhadapt%rBoundary(ibct)%p_IData(BdrNext, ipos)
        
        ! Connect the boundary neighbors with IVT: I1 <- IVT -> I2
        ! First, set IVT as next neighbor of I1
        if (btree_searchInTree(rhadapt%rBoundary(ibct),&
                               i1, ipred) .eq. BTREE_NOT_FOUND) then
          call output_line('Unable to find left neighboring vertex in boundary data structure!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'remove_vertex2D')
          call sys_halt()
        end if
        ipos = rhadapt%rBoundary(ibct)%p_Kchild(merge(TLEFT,TRIGHT, ipred .lt. 0), abs(ipred))
        rhadapt%rBoundary(ibct)%p_IData(BdrNext, ipos) = ivt
        
        ! Second, set IVT as previous neighbor of I2
        if (btree_searchInTree(rhadapt%rBoundary(ibct),&
                               i2, ipred) .eq. BTREE_NOT_FOUND) then
          call output_line('Unable to find right neighboring vertex in boundary data structure!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'remove_vertex2D')
          call sys_halt()
        end if
        ipos = rhadapt%rBoundary(ibct)%p_Kchild(merge(TLEFT,TRIGHT, ipred .lt. 0), abs(ipred))
        rhadapt%rBoundary(ibct)%p_IData(BdrPrev, ipos) = ivt
        
        ! Finally, delete IVTREPLACE from the boundary
        if (btree_deleteFromTree(rhadapt%rBoundary(ibct),&
                                 ivtReplace) .eq. BTREE_NOT_FOUND) then
          call output_line('Unable to delete vertex from the boundary data structure!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'remove_vertex2D')
          call sys_halt()
        end if
      end if
      
      ! Copy data from node IVTREPLACE to node IVT
      rhadapt%p_InodalProperty(ivt) = rhadapt%p_InodalProperty(ivtReplace)
      rhadapt%p_IvertexAge(ivt)     = rhadapt%p_IvertexAge(ivtReplace)

      ! Clear data for node IVTREPLACE
      rhadapt%p_InodalProperty(ivtReplace) = 0
      rhadapt%p_IvertexAge(ivtReplace)     = 0
      
    else

      ! IVT is the last vertex of the adaptivity structure
      ivtReplace = 0

      ! Clear data for node IVT
      rhadapt%p_InodalProperty(ivt) = 0
      rhadapt%p_IvertexAge(ivt)     = 0
      
    end if

  end subroutine remove_vertex2D

  ! ***************************************************************************

!<subroutine>

  subroutine replace_elementTria(rhadapt, ipos, i1, i2, i3,&
                                 e1, e2, e3, e4, e5, e6)
  
!<description>
    ! This subroutine replaces the vertices and adjacent elements for
    ! a given triangular element
!</description>

!<input>
    ! Position number of the element in dynamic data structure
    integer, intent(in) :: ipos

    ! Numbers of the element nodes
    integer, intent(in) :: i1,i2,i3

    ! Numbers of the surrounding elements
    integer, intent(in) :: e1,e2,e3

    ! Numbers of the surrounding mid-elements
    integer, intent(in) :: e4,e5,e6
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    ! Replace triangular element
    rhadapt%p_IverticesAtElement(:,ipos)      = (/i1,i2,i3,0/)
    rhadapt%p_IneighboursAtElement(:,ipos)    = (/e1,e2,e3,0/)
    rhadapt%p_ImidneighboursAtElement(:,ipos) = (/e4,e5,e6,0/)

  end subroutine replace_elementTria

  ! ***************************************************************************
  
!<subroutine>

  subroutine replace_elementQuad(rhadapt, ipos, i1, i2, i3, i4,&
                                 e1, e2, e3, e4, e5, e6, e7, e8)
  
!<description>
    ! This subroutine replaces the vertices and adjacent elements for
    ! a given quadrilateral element
!</description>

!<input>
    ! Position number of the element in dynamic data structure
    integer, intent(in) :: ipos

    ! Numbers of the element nodes
    integer, intent(in) :: i1,i2,i3,i4

    ! Numbers of the surrounding elements
    integer, intent(in) :: e1,e2,e3,e4

    ! Numbers of the surrounding mid-elements
    integer, intent(in) :: e5,e6,e7,e8
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    ! Replace quadrilateral element
    rhadapt%p_IverticesAtElement(:,ipos)      = (/i1,i2,i3,i4/)
    rhadapt%p_IneighboursAtElement(:,ipos)    = (/e1,e2,e3,e4/)
    rhadapt%p_ImidneighboursAtElement(:,ipos) = (/e5,e6,e7,e8/)

  end subroutine replace_elementQuad

  ! ***************************************************************************

!<subroutine>

  subroutine add_elementTria(rhadapt, i1, i2, i3, e1, e2, e3, e4, e5, e6)

!<description>
    ! This subroutine adds a new element connected to three vertices
    ! and surrounded by three adjacent elements
!</description>

!<input>
    ! Numbers of the element nodes
    integer, intent(in) :: i1,i2,i3

    ! Numbers of the surrounding elements
    integer, intent(in) :: e1,e2,e3

    ! Numbers of the surrounding mid-elements
    integer, intent(in) :: e4,e5,e6
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>
    
    ! Increase number of elements and number of triangles
    rhadapt%NEL = rhadapt%NEL+1
    rhadapt%InelOfType(TRIA_NVETRI2D) = rhadapt%InelOfType(TRIA_NVETRI2D)+1

    rhadapt%p_IverticesAtElement(:,rhadapt%NEL)      = (/i1,i2,i3,0/)
    rhadapt%p_IneighboursAtElement(:,rhadapt%NEL)    = (/e1,e2,e3,0/)
    rhadapt%p_ImidneighboursAtElement(:,rhadapt%NEL) = (/e4,e5,e6,0/)

  end subroutine add_elementTria

  ! ***************************************************************************

!<subroutine>
  
  subroutine add_elementQuad(rhadapt, i1, i2, i3, i4,&
                             e1, e2, e3, e4, e5, e6, e7, e8)

!<description>
    ! This subroutine adds a new element connected to four vertices
    ! and surrounded by four adjacent elements
!</description>

!<input>
    ! Numbers of the element nodes
    integer, intent(in) :: i1,i2,i3,i4

    ! Numbers of the surrounding elements
    integer, intent(in) :: e1,e2,e3,e4

    ! Numbers of the surrounding mid-elements
    integer, intent(in) :: e5,e6,e7,e8
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>
    
    ! Increase number of elements and number of quadrilaterals
    rhadapt%NEL = rhadapt%NEL+1
    rhadapt%InelOfType(TRIA_NVEQUAD2D) = rhadapt%InelOfType(TRIA_NVEQUAD2D)+1

    rhadapt%p_IverticesAtElement(:,rhadapt%NEL)      = (/i1,i2,i3,i4/)
    rhadapt%p_IneighboursAtElement(:,rhadapt%NEL)    = (/e1,e2,e3,e4/)
    rhadapt%p_ImidneighboursAtElement(:,rhadapt%NEL) = (/e5,e6,e7,e8/)

  end subroutine add_elementQuad

  ! ***************************************************************************

!<subroutine>

  subroutine remove_element2D(rhadapt, iel, ielReplace)
  
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
    integer :: ipos,ivt,jel,jelmid,ive,jve
    logical :: bfound

    ! Which kind of element are we?
    select case(hadapt_getNVE(rhadapt, iel))
    case(TRIA_NVETRI2D)
      rhadapt%InelOfType(TRIA_NVETRI2D) = rhadapt%InelOfType(TRIA_NVETRI2D)-1

    case(TRIA_NVEQUAD2D)
      rhadapt%InelOfType(TRIA_NVEQUAD2D) = rhadapt%InelOfType(TRIA_NVEQUAD2D)-1

    case DEFAULT
      call output_line('Invalid element type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'remove_element2D')
      call sys_halt()
    end select

    ! Replace element by the last element and delete last element
    ielReplace = rhadapt%NEL

    ! Check if we are the last element
    if (iel .ne. ielReplace) then
      
      ! Element is not the last one. Then the element that should be removed must
      ! have a smaller element number. If this is not the case, something is wrong.
      if (iel .gt.ielReplace) then
        call output_line('Number of replacement element must not be smaller than that of' // &
                         ' the removed elements!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'remove_element2D')
        call sys_halt()
      end if

      ! The element which formally was labeled ielReplace is now labeled IEL.
      ! This modification must be updated in the list of adjacent element
      ! neighbors of all surrounding elements. Moreover, the modified element
      ! number must be updated in the "elements-meeting-at-vertex" lists of the
      ! corner nodes of element IEL. Both operations are performed below.
      update: do ive = 1, hadapt_getNVE(rhadapt, ielReplace)
        
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

                
        ! Get element number of element JEL and JELMID which
        ! are (mid-)adjacent to element ielReplace
        jel    = rhadapt%p_IneighboursAtElement(ive, ielReplace)
        jelmid = rhadapt%p_ImidneighboursAtElement(ive, ielReplace)

        
        ! Do we have different neighbours along the first part and the
        ! second part of the common edge?
        if (jel .ne. jelmid) then
          ! There are two elements sharing the common edge with ielReplace.
          ! We have to find the position of element ielReplace in the lists
          ! of (mid-)adjacent elements for both element JEL and JELMID
          ! seperately. If we are at the boundary of ir element IEL to be
          ! removed is adjacent to element ielReplace, the skip this edge.
          bfound = .false.
          
          ! Find position of replacement element in adjacency list of
          ! element JEL and update the entry to new element number IEL
          if (jel .eq. 0 .or. jel .eq. iel) then
            bfound = .true.
          else
            adjacent1: do jve = 1, hadapt_getNVE(rhadapt, jel)
              if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. ielReplace) then
                  rhadapt%p_IneighboursAtElement(jve, jel)    = iel
                  rhadapt%p_ImidneighboursAtElement(jve, jel) = iel
                bfound = .true.
                exit adjacent1
              end if
            end do adjacent1
          end if
          
          ! Find position of replacement element in adjacentcy list of
          ! element JELMID and update the entry to new element number IEL
          if (jelmid .eq. 0 .or. jelmid .eq. iel) then
            bfound = bfound .and. .true.
          else
            adjacent2: do jve = 1, hadapt_getNVE(rhadapt, jelmid)
              if (rhadapt%p_IneighboursAtElement(jve, jelmid) .eq. ielReplace) then
                  rhadapt%p_IneighboursAtElement(jve, jelmid)    = iel
                  rhadapt%p_ImidneighboursAtElement(jve, jelmid) = iel
                bfound = bfound .and. .true.
                exit adjacent2
              end if
            end do adjacent2
          end if

        else
          ! There is only one element sharing the common edge with ielReplace.
          ! If we are at the boundary or if element IEL to be removed is
          ! adjacent to element ielReplace, then skip this edge.
          if (jel .eq. 0 .or. jel .eq. iel) cycle update

          ! We have to consider two possible situations. The neighbouring element
          ! JEL can be completely aligned with element ielReplace so that the
          ! element number must be updated in the list of adjacent and mid-
          ! adjacent elements. On the other hand, element JEL can share one
          ! half of the common edge with element ielReplace and the other half-
          ! edge with another element. In this case, the number ielReplace can
          ! only be found in either the list of adjacent or mid-adjacent elements.
          bfound = .false.
          adjacent3: do jve = 1, hadapt_getNVE(rhadapt, jel)
            if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. ielReplace) then
                rhadapt%p_IneighboursAtElement(jve, jel) = iel
              bfound = bfound .or. .true.
            end if

            if (rhadapt%p_ImidneighboursAtElement(jve, jel) .eq. ielReplace) then
                rhadapt%p_ImidneighboursAtElement(jve, jel) = iel
              bfound = bfound .or. .true.
            end if
            
            if (bfound) exit adjacent3
          end do adjacent3
          
        end if

        ! If the element could not be found, something is wrong
        if (.not.bfound) then
          call output_line('Unable to update element neighbor!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'remove_element2D')
          call sys_halt()
        end if
        
      end do update
      
      ! Copy data from element ielReplace to element IEL
      rhadapt%p_IverticesAtElement(:,iel) =&
          rhadapt%p_IverticesAtElement(:,ielReplace)
      rhadapt%p_IneighboursAtElement(:,iel) =&
          rhadapt%p_IneighboursAtElement(:,ielReplace)
      rhadapt%p_ImidneighboursAtElement(:,iel) =&
          rhadapt%p_ImidneighboursAtElement(:,ielReplace)
      
    else

      ! Element iel is the last element
      ielReplace = 0
    end if

    ! Decrease number of elements
    rhadapt%NEL = rhadapt%NEL-1

  end subroutine remove_element2D

  ! ***************************************************************************

!<subroutine>
  
  subroutine update_ElemNeighb2D_1to2(rhadapt, jel, jelmid,&
                                      iel0, iel, ielmid)

!<description>
    ! This subroutine updates the list of elements adjacent to another
    ! element and the list of elements mid-adjacent to another element.
    !
    ! The situation is as follows:
    !
    ! <verb>
    !  +---------------------+            +---------------------+
    !  |                     |            |          .          |
    !  |                     |            |          .          |
    !  |        IEL0         |            |  IELMID  .    IEL   |
    !  |                     |            |          .          |
    !  |                     |            |          .          |
    !  +---------------------+    --->    +----------+----------+
    !  |          .          |            |          .          |
    !  |          .          |            |          .          |
    !  |    JEL   .  JELMID  |            |   JEL    .  JELMID  |
    !  |          .          |            |          .          |
    !  |          .          |            |          .          |
    !  +---------------------+            +---------------------+
    ! </verb>
    !
    ! "Sitting" on element IEL0 we want to update the element lists:
    !
    ! 1) If IEL0 is located at the boundary then nothing needs to be done.
    ! 2) If the adjacent element has not been subdivided, that is,
    !    JEL = JELMID then it suffices to update its adjacency
    !    list for the entry IEL0 adopting the values IEL and IELMID.
    ! 3) If the adjacent element has been subdivided, that is,
    !    JEL != JELMID then the adjacency list if each of these two
    !    elements is updated by the value IEL and IELMID, respectively.
!</description>

!<input>
    ! Number of the neighboring element
    integer, intent(in) :: jel

    ! Number of the mid-neighboring element
    integer, intent(in) :: jelmid

    ! Number of the updated macro-element
    integer, intent(in) :: iel0

    ! Number of the new neighboring element
    integer, intent(in) :: iel

    ! Number of the new mid-neighboring element
    integer, intent(in) :: ielmid
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: ive,nve
    logical :: bfound1,bfound2,bfound

    ! Do nothing for elements adjacent to the boundary
    if (jel .eq. 0 .or. jelmid .eq. 0) return

    ! Check if adjacent and mid-adjacent elements are the same.
    if (jel .eq. jelmid) then

      ! Case 1: Adjacent element has not been subdivided.
      bfound1 = .false.; bfound2 = .false.

      ! What kind of element is neighboring element?
      nve = hadapt_getNVE(rhadapt, jel)
      
      ! Loop over all entries in the list of adjacent and/or mid-adjacent
      ! elements for element JEL and check if the value IEL0 is present.
      ! The situation may occur that element IEL0 is only adjacent to one
      ! "half" of the edge of element JEL. This may be the case, if element
      ! IEL0 was a green triangle which has already been subdivided in the
      ! marking procedure, whiereby element JEL is marked for red refinement.
      ! In order to consider this situation we are looking for element IEL0
      ! both in the list of adjacent and mid-adjacent elements of JEL.
      ! It suffices if IEL0 is found in one of these lists.
      do ive = 1, nve
        if (rhadapt%p_IneighboursAtElement(ive, jel) .eq. iel0) then
          rhadapt%p_IneighboursAtElement(ive, jel) = iel
          bfound1 = .true.
        end if

        if (rhadapt%p_ImidneighboursAtElement(ive, jel) .eq. iel0) then
          rhadapt%p_ImidneighboursAtElement(ive, jel) = ielmid
          bfound2 = .true.
        end if
        
        ! Exit if IEL0 has been found in the either the adjacency or the
        ! mid-adjacency list of element JEL.
        bfound = bfound1 .or. bfound2
        if (bfound) exit
      end do

    else
      
      ! Case 2: Adjacent element has already been subdivided.
      bfound1 = .false.; bfound2 = .false.
      
      ! What kind of element is neighboring element
      nve = hadapt_getNVE(rhadapt, jel)

      ! Loop over all entries in the list of adjacent elements for element
      ! JEL and check if the value IEL0 is present.
      ! If this is the case,  then update the corrseponding entries in the
      ! lists of (mid-)adjacent element neighbors.
      do ive = 1, nve
        if (rhadapt%p_IneighboursAtElement(ive, jel) .eq. iel0) then
          rhadapt%p_IneighboursAtElement(ive, jel)    = ielmid
          rhadapt%p_ImidneighboursAtElement(ive, jel) = ielmid
          bfound1 = .true.
          exit
        end if
      end do
      
      ! What kind of element is neighboring element
      nve = hadapt_getNVE(rhadapt, jelmid)
      
      ! Loop over all entries in the list of adjacent elements for element
      ! JELMID and check if the value IEL0 is present.
      ! If this is the case, then update the corrseponding entries in the
      ! lists of (mid-)adjacent element neighbors.
      do ive = 1, nve
        if (rhadapt%p_IneighboursAtElement(ive, jelmid) .eq. iel0) then
          rhadapt%p_IneighboursAtElement(ive, jelmid)    = iel
          rhadapt%p_ImidneighboursAtElement(ive, jelmid) = iel
          bfound2 = .true.
          exit
        end if
      end do

      ! Check success of both searches
      bfound = bfound1 .and. bfound2

    end if

    if (.not.bfound) then
      call output_line('Inconsistent adjacency lists!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'update_ElemNeighb2D_1to2')
      call sys_halt()
    end if

  end subroutine update_ElemNeighb2D_1to2

  ! ***************************************************************************

!<subroutine>
  
  subroutine update_ElemNeighb2D_2to2(rhadapt, jel, jelmid,&
                                      iel0, ielmid0, iel, ielmid)

!<description>
    ! This subroutine updates the list of elements adjacent to another
    ! element and the list of elements mid-adjacent to another element.
    !
    ! The situation is as follows:
    !
    ! <verb>
    !  +---------------------+            +---------------------+
    !  |          .          |            |          .          |
    !  |          .          |            |          .          |
    !  | IELMID0  .   IEL0   |            | IELMID   .   IEL    |
    !  |          .          |            |          .          |
    !  |          .          |            |          .          |
    !  +---------------------+    --->    +----------+----------+
    !  |          .          |            |          .          |
    !  |          .          |            |          .          |
    !  |    JEL   .  JELMID  |            |    JEL   .   JELMID |
    !  |          .          |            |          .          |
    !  |          .          |            |          .          |
    !  +---------------------+            +---------------------+
    ! </verb>
    !
    ! If IEL0 and IELMID0 are the same, then thins subroutine is identical
    ! to subroutine update_EkemNeighb2D_1to2 which is called in this case
!</description>

!<input>
    ! Number of the neighboring element
    integer, intent(in) :: jel

    ! Number of the mid-neighboring element
    integer, intent(in) :: jelmid

    ! Number of the updated macro-element
    integer, intent(in) :: iel0

    ! Number of the updated macro-element
    integer, intent(in) :: ielmid0

    ! Number of the new neighboring element
    integer, intent(in) :: iel

    ! Number of the new mid-neighboring element
    integer, intent(in) :: ielmid
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: ive,nve
    logical :: bfound1,bfound2,bfound
    
    ! Check if IEL0 and IELMID0 are the same.
    ! In this case call the corresponding subroutine
    if (iel0 .eq. ielmid0) then
      call update_ElemNeighb2D_1to2(rhadapt, jel, jelmid,&
                                    iel0, iel, ielmid)
      return
    end if

    ! Do nothing for elements adjacent to the boundary
    if (jel .eq. 0 .or. jelmid .eq. 0) return
    
    ! Check if adjacent and mid-adjacent elements are the same.
    if (jel .eq. jelmid) then

      ! Case 1: Adjacent element has not been subdivided, that is, the
      !         current edge contains a hanging node for the adjacent element.
      bfound = .false.

      ! What kind of element is neighboring element?
      nve = hadapt_getNVE(rhadapt, jel)
      
      ! Loop over all entries in the list of adjacent elements for element
      ! JEL and check if the value IEL0 or IELMID0 is present.
      ! If this is the case, then update the corrseponding entries in the
      ! lists of (mid-)adjacent element neighbors.
      do ive = 1, nve
        if (rhadapt%p_IneighboursAtElement(ive, jel)    .eq. iel0 .and.&
            rhadapt%p_ImidneighboursAtElement(ive, jel) .eq. ielmid0) then
          rhadapt%p_IneighboursAtElement(ive, jel)    = iel
          rhadapt%p_ImidneighboursAtElement(ive, jel) = ielmid
          bfound = .true.
          exit
        end if
      end do
                
    else
      
      ! Case 2: Adjacent element has already been subdivided.
      bfound1 = .false.; bfound2 = .false.
      
      ! What kind of element is neighboring element
      nve = hadapt_getNVE(rhadapt, jel)

      ! Loop over all entries in the list of adjacent elements for element
      ! JEL and check if the value IELMID0 is present.
      ! If this is the case,  then update the corrseponding entries in the
      ! lists of (mid-)adjacent element neighbors.
      do ive = 1, nve
        if (rhadapt%p_IneighboursAtElement(ive, jel) .eq. ielmid0 .and.&
            rhadapt%p_ImidneighboursAtElement(ive, jel) .eq. ielmid0) then
          rhadapt%p_IneighboursAtElement(ive, jel)    = ielmid
          rhadapt%p_ImidneighboursAtElement(ive, jel) = ielmid
          bfound1 = .true.
          exit
        end if
      end do
      
      ! What kind of element is neighboring element
      nve = hadapt_getNVE(rhadapt, jelmid)
      
      ! Loop over all entries in the list of adjacent elements for element
      ! JELMID and check if the value IEL0 is present.
      ! If this is the case, then update the corrseponding entries in the
      ! lists of (mid-)adjacent element neighbors.
      do ive = 1, nve
        if (rhadapt%p_IneighboursAtElement(ive, jelmid) .eq. iel0 .and.&
            rhadapt%p_ImidneighboursAtElement(ive, jelmid) .eq. iel0) then
          rhadapt%p_IneighboursAtElement(ive, jelmid)    = iel
          rhadapt%p_ImidneighboursAtElement(ive, jelmid) = iel
          bfound2 = .true.
          exit
        end if
      end do

      bfound = bfound1 .and. bfound2
    end if

    if (.not.bfound) then
      call output_line('Inconsistent adjacency lists!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'update_ElemNeighb2D_2to2')
      call sys_halt()
    end if

  end subroutine update_ElemNeighb2D_2to2

  ! ***************************************************************************

!<subroutine>

  subroutine update_AllElementNeighbors2D(rhadapt, iel0, iel)

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
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: jel,ive,jve
    logical :: bfound

    ! Check if the old element is still present in the triangulation
    if (iel0 .gt. rhadapt%NEL) return

    ! Loop over adjacent elements
    adjacent: do ive = 1, hadapt_getNVE(rhadapt, iel0)
      ! Get number of adjacent element
      jel = rhadapt%p_IneighboursAtElement(ive, iel0)

      ! Are we at the boundary?
      if (jel .eq. 0) cycle adjacent

      ! Initialise indicator
      bfound = .false.

      ! Find position of element IEL0 in adjacent element JEL
      do jve = 1, hadapt_getNVE(rhadapt, jel)
        if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. iel0) then
          rhadapt%p_IneighboursAtElement(jve, jel) = iel
          bfound = .true.
          exit
        end if
      end do

      ! Find position of element IEL0 in mid-adjacent element JEL
      do jve = 1, hadapt_getNVE(rhadapt, jel)
        if (rhadapt%p_ImidneighboursAtElement(jve, jel) .eq. iel0) then
          rhadapt%p_ImidneighboursAtElement(jve, jel) = iel
          bfound = (bfound .and. .true.)
          exit
        end if
      end do

      ! If the old element number was not found in adjacent element JEL
      ! then something went wrong and we should not proceed.
      if (.not.bfound) then
        call output_line('Inconsistent adjacency lists!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'update_AllElementNeighbors2D')
        call sys_halt()
      end if
    end do adjacent

  end subroutine update_AllElementNeighbors2D

  ! ***************************************************************************

!<subroutine>

  subroutine refine_Tria2Tria(rhadapt, iel, imarker,&
                              rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one triangular element into two triangular
    ! elements by subdividing one edge. The local number of the edge that is
    ! subdivided can be uniquely determined from the element marker.
    !
    ! In the illustration, i1,i2,i3 and i4 denote the vertices whereby i4
    ! is the new vertex. Moreover, (e1)-(e6) stand for the element numbers
    ! which are adjecent and mid-adjacent to the current element.
    ! The new element is assigned the total number of elements currently
    ! present in the triangulation increased by one.
    !
    ! <verb>
    !    initial triangle           subdivided triangle
    !
    !            i3                            i3
    !            +                             +
    !           / \                           /|\
    !     (e3) /   \ (e5)               (e3) / | \ (e5)
    !         /     \                       /  |  \
    !        /       \          ->         /   |   \
    !       /   iel   \                   /    |    \
    ! (e6) /           \ (e2)       (e6) / iel |nel+1\ (e2)
    !     /X            \               /X     |     X\
    !    +---------------+             +-------+-------+
    !   i1 (e1)     (e4) i2            i1 (e1) i4 (e4) i2
    ! </verb>
    !
!</description>

!<input>
    ! Number of element to be refined
    integer, intent(in) :: iel
    
    ! Identifier for element marker
    integer, intent(in) :: imarker

    ! Callback function
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
    integer :: ipos
    integer :: nel0,e1,e2,e3,e4,e5,e6
    integer :: i1,i2,i3,i4
    integer :: loc1,loc2,loc3
    
    ! Determine the local position of edge to
    ! be subdivided from the element marker
    select case(imarker)
    case(MARK_REF_TRIA2TRIA_1)
      loc1=1; loc2=2; loc3=3

    case(MARK_REF_TRIA2TRIA_2)
      loc1=2; loc2=3; loc3=1

    case(MARK_REF_TRIA2TRIA_3)
      loc1=3; loc2=1; loc3=2

    case DEFAULT
      call output_line('Invalid element marker!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'refine_Tria2Tria')
      call sys_halt()
    end select
    
    ! Store vertex- and element-values of the current element
    i1 = rhadapt%p_IverticesAtElement(loc1, iel)
    i2 = rhadapt%p_IverticesAtElement(loc2, iel)
    i3 = rhadapt%p_IverticesAtElement(loc3, iel)

    e1 = rhadapt%p_IneighboursAtElement(loc1, iel)
    e2 = rhadapt%p_IneighboursAtElement(loc2, iel)
    e3 = rhadapt%p_IneighboursAtElement(loc3, iel)
    
    e4 = rhadapt%p_ImidneighboursAtElement(loc1, iel)
    e5 = rhadapt%p_ImidneighboursAtElement(loc2, iel)
    e6 = rhadapt%p_ImidneighboursAtElement(loc3, iel)
    
    ! Store total number of elements before refinement
    nel0 = rhadapt%NEL

    
    ! Add one new vertex I4 at the midpoint of edge (I1,I2).
    call add_vertex2D(rhadapt, i1, i2, e1, i4,&
                      rcollection, fcb_hadaptCallback)
    
    ! Replace element IEL and add one new element numbered NEL0+1
    call replace_element2D(rhadapt, iel, i1, i4, i3,&
                           e1, nel0+1, e3, e1, nel0+1, e6)
    call add_element2D(rhadapt, i2, i3, i4,&
                       e2, iel, e4, e5, iel, e4)

    
    ! Update list of neighboring elements
    call update_ElementNeighbors2D(rhadapt, e1, e4, iel, nel0+1, iel)
    call update_ElementNeighbors2D(rhadapt, e2, e5, iel, nel0+1, nel0+1)

    
    ! Update list of elements meeting at vertices
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i2, iel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'refine_Tria2Tria')
      call sys_halt()
    end if

    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i2, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i3, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i4, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i4, iel,    ipos)


    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:10) = (/i1, i2, i3, i4, &
                                         e1, e2, e3, e4, e5, e6/)
      call fcb_hadaptCallback(HADAPT_OPR_REF_TRIA2TRIA, rcollection)
    end if

  end subroutine refine_Tria2Tria

  ! ***************************************************************************
  
!<subroutine>

  subroutine refine_Tria3Tria(rhadapt, iel, imarker,&
                              rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one triangular element into three triangular
    ! elements by subdividing the longest edge and connecting the new vertex
    ! with the opposite midpoint. The local numbers of the two edges that are
    ! subdivided can be uniquely determined from the element marker.
    !
    ! In the illustration, i1,i2,i3,i4 and i5 denote the vertices whereby i4 and
    ! i5 are the new vertices. Moreover, (e1)-(e6) stand for the element numbers
    ! which are adjecent and mid-adjacent to the current element.
    ! The new elements are assigned the total number of elements currently present
    ! in the triangulation increased by one and two, respectively.
    !
    ! <verb>
    !    initial triangle           subdivided triangle
    !
    !            i3                        i3
    !            +                         +
    !           / \                       /|\
    !     (e3) /   \ (e5)           (e3) / |X\ (e5)
    !         /     \                   /  |ne\
    !        /       \          ->     /   |l+2+i5
    !       /   iel   \               /    |  / \
    ! (e6) /           \ (e2)   (e6) / iel | /nel\ (e2)
    !     /X            \           /X     |/ +1 X\
    !    +---------------+         +-------+-------+
    !   i1 (e1)     (e4) i2       i1 (e1) i4  (e4) i2
    ! </verb>
    !
!</description>

!<input>
    ! Number of element to be refined
    integer, intent(in) :: iel
    
    ! Identifier for element marker
    integer, intent(in) :: imarker

    ! Callback function
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
    integer :: ipos
    integer :: nel0,e1,e2,e3,e4,e5,e6
    integer :: i1,i2,i3,i4,i5
    integer :: loc1,loc2,loc3
    real(DP) :: dlen12,dlen23,x,y
    
    ! Find corresponding edges according to convention to be subdivided
    select case(imarker)
    case(MARK_REF_TRIA3TRIA_12)
      loc1=1; loc2=2; loc3=3

    case(MARK_REF_TRIA3TRIA_23)
      loc1=2; loc2=3; loc3=1

    case(MARK_REF_TRIA3TRIA_13)
      loc1=3; loc2=1; loc3=2

    case DEFAULT
      call output_line('Invalid element marker!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'refine_Tria3Tria')
      call sys_halt()
    end select

    ! Store vertex- and element-values of the current element
    i1 = rhadapt%p_IverticesAtElement(loc1, iel)
    i2 = rhadapt%p_IverticesAtElement(loc2, iel)
    i3 = rhadapt%p_IverticesAtElement(loc3, iel)

    e1 = rhadapt%p_IneighboursAtElement(loc1, iel)
    e2 = rhadapt%p_IneighboursAtElement(loc2, iel)
    e3 = rhadapt%p_IneighboursAtElement(loc3, iel)

    e4 = rhadapt%p_ImidneighboursAtElement(loc1, iel)
    e5 = rhadapt%p_ImidneighboursAtElement(loc2, iel)
    e6 = rhadapt%p_ImidneighboursAtElement(loc3, iel)
    
    ! Store total number of elements before refinement
    nel0 = rhadapt%NEL
    

    ! Add two new vertices I4 and I5 at the midpoint of edges (I1,I2) and (I2,I3).
    call add_vertex2D(rhadapt, i1, i2, e1, i4,&
                      rcollection, fcb_hadaptCallback)
    call add_vertex2D(rhadapt, i2, i3, e2, i5,&
                      rcollection, fcb_hadaptCallback)
    
    ! Compute the length of edges (I1,I2) and (I2,I3)
    x = qtree_getX(rhadapt%rVertexCoordinates2D, i1)-&
        qtree_getX(rhadapt%rVertexCoordinates2D, i2)
    y = qtree_getY(rhadapt%rVertexCoordinates2D, i1)-&
        qtree_getY(rhadapt%rVertexCoordinates2D, i2)
    dlen12 = sqrt(x*x+y*y)

    x = qtree_getX(rhadapt%rVertexCoordinates2D, i2)-&
        qtree_getX(rhadapt%rVertexCoordinates2D, i3)
    y = qtree_getY(rhadapt%rVertexCoordinates2D, i2)-&
        qtree_getY(rhadapt%rVertexCoordinates2D, i3)
    dlen23 = sqrt(x*x+y*y)
    
    if (dlen12 .gt.dlen23) then
      
      ! 1st CASE: longest edge is (I1,I2)
      
      ! Replace element IEL and add two new elements numbered NEL0+1 and NEL0+2
      call replace_element2D(rhadapt, iel, i1, i4, i3,&
                             e1, nel0+2, e3, e1, nel0+2, e6)
      call add_element2D(rhadapt, i2, i5, i4,&
                         e2, nel0+2, e4, e2, nel0+2, e4)
      call add_element2D(rhadapt, i3, i4, i5,&
                         iel, nel0+1, e5, iel, nel0+1, e5)
      

      ! Update list of neighboring elements
      call update_ElementNeighbors2D(rhadapt, e1, e4, iel, nel0+1, iel)
      call update_ElementNeighbors2D(rhadapt, e2, e5, iel, nel0+2, nel0+1)

            
      ! Update list of elements meeting at vertices
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i2, iel) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'refine_Tria3Tria')
        call sys_halt()
      end if
      
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i2, iel,    ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i4, iel,    ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i4, nel0+1, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i4, nel0+2, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, nel0+1, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, nel0+2, ipos)


      ! Optionally, invoke callback function
      if (present(fcb_hadaptCallback) .and. present(rcollection)) then
        rcollection%IquickAccess(1:10) = (/i1, i2, i3, i4, i5, &
                                           e1, e2, e3, e4, e5/)
        call fcb_hadaptCallback(HADAPT_OPR_REF_TRIA3TRIA12, rcollection)
      end if
      
    else
      
      ! 2nd CASE: longest edge is (I2,I3)
      
      ! Replace element IEL and add two new elements numbered NEL0+1 and NEL0+2
      call replace_element2D(rhadapt, iel, i1, i5, i3,&
                             nel0+2, e5, e3, nel0+2, e5, e6)
      call add_element2D(rhadapt, i2, i5, i4,&
                         e2, nel0+2, e4, e2, nel0+2, e4)
      call add_element2D(rhadapt, i1, i4, i5,&
                         e1, nel0+1, iel, e1, nel0+1, iel)
      
      
      ! Update list of neighboring elements
      call update_ElementNeighbors2D(rhadapt, e1, e4, iel, nel0+1, nel0+2)
      call update_ElementNeighbors2D(rhadapt, e2, e5, iel, iel, nel0+1)
      
      
      ! Update list of elements meeting at vertices
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i2, iel) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'refine_Tria3Tria')
        call sys_halt()
      end if

      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i2, nel0+1, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i1, nel0+2, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i4, nel0+1, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i4, nel0+2, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, iel,    ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, nel0+1, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, nel0+2, ipos)


      ! Optionally, invoke callback function
      if (present(fcb_hadaptCallback) .and. present(rcollection)) then
        rcollection%IquickAccess(1:10) = (/i1, i2, i3, i4, i5,&
                                           e1, e2, e3, e4, e5/)
        call fcb_hadaptCallback(HADAPT_OPR_REF_TRIA3TRIA23, rcollection)
      end if
    end if

  end subroutine refine_Tria3Tria

  ! ***************************************************************************

!<subroutine>

  subroutine refine_Tria4Tria(rhadapt, iel, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one triangular element into four similar
    ! triangular elements by connecting the three edge midpoints.
    !
    ! In the illustration, i1-i6 denote the vertices whereby i4-i6 are the
    ! new vertices. Moreover, (e1)-(e6) stand for the element numbers which
    ! are adjecent and mid-adjacent to the current element.
    ! The new elements are assigned the total number of elements currently present
    ! in the triangulation increased by one, two and three, respectively.
    !
    ! <verb>
    !    initial triangle           subdivided triangle
    !
    !            i3                        i3
    !            +                         +
    !           / \                       /X\
    !     (e3) /   \ (e5)           (e3) /nel\ (e5)
    !         /     \                   / +3  \
    !        /       \          ->  i6 +-------+ i5
    !       /   iel   \               / \ iel / \
    ! (e6) /           \ (e2)   (e6) /nel\   /nel\ (e2)
    !     /X            \           /X +1 \X/ +2 X\
    !    +---------------+         +-------+-------+
    !   i1 (e1)     (e4) i2       i1 (e1)  i4 (e4) i2
    ! </verb>
    !
!</description>

!<input>
    ! Number of element to be refined
    integer, intent(in) :: iel

    ! Callback function
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
    integer :: ipos
    integer :: nel0,e1,e2,e3,e4,e5,e6
    integer :: i1,i2,i3,i4,i5,i6

    ! Store vertex- and element-values of the current element
    i1 = rhadapt%p_IverticesAtElement(1, iel)
    i2 = rhadapt%p_IverticesAtElement(2, iel)
    i3 = rhadapt%p_IverticesAtElement(3, iel)

    e1 = rhadapt%p_IneighboursAtElement(1, iel)
    e2 = rhadapt%p_IneighboursAtElement(2, iel)
    e3 = rhadapt%p_IneighboursAtElement(3, iel)

    e4 = rhadapt%p_ImidneighboursAtElement(1, iel)
    e5 = rhadapt%p_ImidneighboursAtElement(2, iel)
    e6 = rhadapt%p_ImidneighboursAtElement(3, iel)
        
    ! Store total number of elements before refinement
    nel0 = rhadapt%NEL


    ! Add three new vertices I4,I5 and I6 at the midpoint of
    ! edges (I1,I2), (I2,I3) and (I1,I3), respectively.
    call add_vertex2D(rhadapt, i1, i2, e1, i4,&
                      rcollection, fcb_hadaptCallback)
    call add_vertex2D(rhadapt, i2, i3, e2, i5,&
                      rcollection, fcb_hadaptCallback)
    call add_vertex2D(rhadapt, i3, i1, e3, i6,&
                      rcollection, fcb_hadaptCallback)
    

    ! Replace element IEL and add three new elements NEL0+1, NEL0+2 and NEL0+3
    call replace_element2D(rhadapt, iel, i4, i5, i6,&
                           nel0+2, nel0+3, nel0+1, nel0+2, nel0+3, nel0+1)
    call add_element2D(rhadapt, i1, i4, i6, e1, iel, e6, e1, iel, e6)
    call add_element2D(rhadapt, i2, i5, i4, e2, iel, e4, e2, iel, e4)
    call add_element2D(rhadapt, i3, i6, i5, e3, iel, e5, e3, iel, e5)

    
    ! Update list of neighboring elements
    call update_ElementNeighbors2D(rhadapt, e1, e4, iel, nel0+2, nel0+1)
    call update_ElementNeighbors2D(rhadapt, e2, e5, iel, nel0+3, nel0+2)
    call update_ElementNeighbors2D(rhadapt, e3, e6, iel, nel0+1, nel0+3)

    
    ! Update list of elements meeting at vertices
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i1, iel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'refine_Tria4Tria')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i2, iel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'refine_Tria4Tria')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i3, iel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'refine_Tria4Tria')
      call sys_halt()
    end if
   
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i1, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i2, nel0+2, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i3, nel0+3, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i4, iel   , ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i4, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i4, nel0+2, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, iel   , ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, nel0+2, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, nel0+3, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i6, iel   , ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i6, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i6, nel0+3, ipos)


    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:12) = (/i1, i2, i3, i4, i5, i6,&
                                         e1, e2, e3, e4, e5, e6/)
      call fcb_hadaptCallback(HADAPT_OPR_REF_TRIA4TRIA, rcollection)
    end if

  end subroutine refine_Tria4Tria

  ! ***************************************************************************

!<subroutine>
  
  subroutine refine_Quad2Quad(rhadapt, iel, imarker,&
                              rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one quadrilateral element into two quadrilateral
    ! elements. The local number of the edge that is subdivided can be uniquely
    ! determined from the element marker.
    !
    ! In the illustration, i1-i6 denote the vertices whereby i5 and i6
    ! are the new vertices. Moreover, (e1)-(e8) stand for the element numbers
    ! which are adjecent and mid-adjacent to the current element.
    ! The new element is assigned the total number of elements currently present
    ! in the triangulation increased by one.
    !
    ! <verb>
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7)     (e3) i3          i4 (e7)  i6 (e3) i3
    !      +---------------+            +-------+-------+
    !      |               |            |       |      *|
    ! (e4) |               | (e6)  (e4) |       |       | (e6)
    !      |               |            |       |       |
    !      |     iel       |     ->     |  iel  | nel+1 |
    !      |               |            |       |       |
    ! (e8) |               | (e2)  (e8) |       |       | (e2)
    !      |*              |            |*      |       |
    !      +---------------+            +-------+-------+
    !     i1 (e1)     (e5) i2          i1 (e1)  i5 (e5) i2
    ! </verb>
    !
!</description>

!<input>
    ! Number of element to be refined
    integer, intent(in) :: iel

    ! Identifier for element marker
    integer, intent(in) :: imarker

    ! Callback function
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
    integer :: ipos
    integer :: nel0,e1,e2,e3,e4,e5,e6,e7,e8
    integer :: i1,i2,i3,i4,i5,i6
    integer :: loc1,loc2,loc3,loc4
    
    ! Find local position of edge to be subdivided
    select case(imarker)
    case(MARK_REF_QUAD2QUAD_13)
      loc1=1; loc2=2; loc3=3; loc4=4

    case(MARK_REF_QUAD2QUAD_24)
      loc1=2; loc2=3; loc3=4; loc4=1

    case DEFAULT
      call output_line('Invalid element marker!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad2Quad')
      call sys_halt()
    end select
    
    ! Store vertex- and element-values of the current element
    i1 = rhadapt%p_IverticesAtElement(loc1, iel)
    i2 = rhadapt%p_IverticesAtElement(loc2, iel)
    i3 = rhadapt%p_IverticesAtElement(loc3, iel)
    i4 = rhadapt%p_IverticesAtElement(loc4, iel)

    e1 = rhadapt%p_IneighboursAtElement(loc1, iel)
    e2 = rhadapt%p_IneighboursAtElement(loc2, iel)
    e3 = rhadapt%p_IneighboursAtElement(loc3, iel)
    e4 = rhadapt%p_IneighboursAtElement(loc4, iel)

    e5 = rhadapt%p_ImidneighboursAtElement(loc1, iel)
    e6 = rhadapt%p_ImidneighboursAtElement(loc2, iel)
    e7 = rhadapt%p_ImidneighboursAtElement(loc3, iel)
    e8 = rhadapt%p_ImidneighboursAtElement(loc4, iel)
    
    ! Store total number of elements before refinement
    nel0 = rhadapt%NEL
    
    ! Add two new vertices I5 and I6 at the midpoint of edges (I1,I2) and
    ! (I3,I4), respectively.
    call add_vertex2D(rhadapt, i1, i2, e1, i5,&
                      rcollection, fcb_hadaptCallback)
    call add_vertex2D(rhadapt, i3, i4, e3, i6,&
                      rcollection, fcb_hadaptCallback)
    

    ! Replace element IEL and add one new element NEL0+1
    call replace_element2D(rhadapt, iel, i1, i5, i6, i4,&
                           e1, nel0+1, e7, e4, e1, nel0+1, e7, e8)
    call add_element2D(rhadapt, i3, i6, i5, i2,&
                       e3, iel, e5, e2, e3, iel, e5, e6)
    

    ! Update list of neighboring elements
    call update_ElementNeighbors2D(rhadapt, e1, e5, iel, nel0+1, iel)
    call update_ElementNeighbors2D(rhadapt, e2, e6, iel, nel0+1, nel0+1)
    call update_ElementNeighbors2D(rhadapt, e3, e7, iel, iel, nel0+1)

    
    ! Update list of elements meeting at vertices
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i2, iel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad2Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i3, iel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad2Quad')
      call sys_halt()
    end if

    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i2, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i3, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, iel,    ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i6, iel,    ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i6, nel0+1, ipos)
    

    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:14) = (/i1, i2, i3, i4, i5, i6,&
                                         e1, e2, e3, e4, e5, e6, e7, e8/)
      call fcb_hadaptCallback(HADAPT_OPR_REF_QUAD2QUAD, rcollection)
    end if

  end subroutine refine_Quad2Quad

  ! ***************************************************************************

!<subroutine>

  subroutine refine_Quad3Tria(rhadapt, iel, imarker,&
                              rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one quadrilateral element into three triangular
    ! elements. The local number of the edge that is
    ! subdivided can be uniquely determined from the element marker.
    !
    ! In the illustration, i1-i5 denote the vertices whereby i5 is the new
    ! vertex. Moreover, (e1)-(e8) stand for the element numbers which are
    ! adjecent and mid-adjacent to the current element.
    ! The new elements are assigned the total number of elements currently
    ! present in the triangulation increased by one and two, respectively.
    !
    ! <verb>
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7)      (e3) i3          i4 (e7)     (e3) i3
    !      +----------------+            +---------------+
    !      |                |            |\             /|
    ! (e4) |                | (e6)  (e4) | \           / | (e6)
    !      |                |            |  \  nel+2  /  |
    !      |      iel       |     ->     |   \       /   |
    !      |                |            |    \     /    |
    ! (e8) |                | (e2)  (e8) |     \   / nel | (e2)
    !      |*               |            |* iel \X/  +1 *|
    !      +----------------+            +-------+-------+
    !     i1 (e1)      (e5) i2          i1 (e1)  i5 (e5) i2
    ! </verb>
    !
!</description>

!<input>
    ! Number of element to be refined
    integer, intent(in) :: iel

    ! Identifier for element marker
    integer, intent(in) :: imarker

    ! Callback function
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
    integer :: ipos
    integer :: nel0,e1,e2,e3,e4,e5,e6,e7,e8
    integer :: i1,i2,i3,i4,i5
    integer :: loc1,loc2,loc3,loc4
    
    ! Find local position of edge to be subdivided
    select case(imarker)
    case(MARK_REF_QUAD3TRIA_1)
      loc1=1; loc2=2; loc3=3; loc4=4

    case(MARK_REF_QUAD3TRIA_2)
      loc1=2; loc2=3; loc3=4; loc4=1

    case(MARK_REF_QUAD3TRIA_3)
      loc1=3; loc2=4; loc3=1; loc4=2

    case(MARK_REF_QUAD3TRIA_4)
      loc1=4; loc2=1; loc3=2; loc4=3

    case DEFAULT
      call output_line('Invalid element marker!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad3Tria')
      call sys_halt()
    end select
    

    ! Store vertex- and element-values of the current element
    i1 = rhadapt%p_IverticesAtElement(loc1, iel)
    i2 = rhadapt%p_IverticesAtElement(loc2, iel)
    i3 = rhadapt%p_IverticesAtElement(loc3, iel)
    i4 = rhadapt%p_IverticesAtElement(loc4, iel)


    e1 = rhadapt%p_IneighboursAtElement(loc1, iel)
    e2 = rhadapt%p_IneighboursAtElement(loc2, iel)
    e3 = rhadapt%p_IneighboursAtElement(loc3, iel)
    e4 = rhadapt%p_IneighboursAtElement(loc4, iel)

    e5 = rhadapt%p_ImidneighboursAtElement(loc1, iel)
    e6 = rhadapt%p_ImidneighboursAtElement(loc2, iel)
    e7 = rhadapt%p_ImidneighboursAtElement(loc3, iel)
    e8 = rhadapt%p_ImidneighboursAtElement(loc4, iel)

    ! Store total number of elements before refinement
    nel0 = rhadapt%NEL
    

    ! Add one new vertex I5 it the midpoint of edge (I1,I2)
    call add_vertex2D(rhadapt, i1, i2, e1, i5,&
                      rcollection, fcb_hadaptCallback)

    
    ! Replace element IEL and add two new elements NEL0+1 and NEL0+2
    call replace_element2D(rhadapt, iel, i1, i5, i4,&
                           e1, nel0+2, e4, e1, nel0+2, e8)
    call add_element2D(rhadapt, i2, i3, i5,&
                       e2, nel0+2, e5, e6, nel0+2, e5)
    call add_element2D(rhadapt, i5, i3, i4,&
                       nel0+1, e3, iel, nel0+1, e7, iel)
    

    ! Update list of neighboring elements
    call update_ElementNeighbors2D(rhadapt, e1, e5, iel, nel0+1, iel)
    call update_ElementNeighbors2D(rhadapt, e2, e6, iel, nel0+1, nel0+1)
    call update_ElementNeighbors2D(rhadapt, e3, e7, iel, nel0+2, nel0+2)
    

    ! Update list of elements meeting at vertices
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i2, iel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad3Tria')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i3, iel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad3Tria')
      call sys_halt()
    end if

    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i2, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i3, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i3, nel0+2, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i4, nel0+2, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, iel,    ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, nel0+2, ipos)


    ! Adjust number of elements
    rhadapt%InelOfType(TRIA_NVETRI2D)  = rhadapt%InelOfType(TRIA_NVETRI2D)+1
    rhadapt%InelOfType(TRIA_NVEQUAD2D) = rhadapt%InelOfType(TRIA_NVEQUAD2D)-1
    

    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:13) = (/i1, i2, i3, i4, i5,&
                                         e1, e2, e3, e4, e5, e6, e7, e8/)
      call fcb_hadaptCallback(HADAPT_OPR_REF_QUAD3TRIA, rcollection)
    end if

  end subroutine refine_Quad3Tria
  
  ! ***************************************************************************

!<subroutine>

  subroutine refine_Quad4Tria(rhadapt, iel, imarker,&
                              rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one quadrilateral element into four triangular
    ! elements. The local numbers of the edges that are subdivided can be
    ! uniquely determined from the element marker.
    !
    ! In the illustration, i1-i6 denote the vertices whereby i5 and i6 are the
    ! new vertices. Moreover, (e1)-(e8) stand for the element numbers which are
    ! adjecent and mid-adjacent to the current element.
    ! The new elements are assigned the total number of elements currently
    ! present in the triangulation increased by one, two an three, respectively.
    !
    ! <verb>
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7)      (e3) i3          i4 (e7)     (e3) i3
    !      +----------------+            +---------------+
    !      |                |            |\\\\          *|
    ! (e4) |                | (e6)  (e4) | \*  \\\ nel+2 | (e6)
    !      |                |            |  \      \\\   |
    !      |      iel       |     ->     |   \         \\+i6
    !      |                |            |    \ nel+3  / |
    ! (e8) |                | (e2)  (e8) | iel \     /nel| (e2)
    !      |*               |            |*     \  / +1 *|
    !      +----------------+            +-------+-------+
    !     i1 (e1)      (e5) i2          i1 (e1)  i5 (e5) i2
    ! </verb>
    !
!</description>

!<input>
    ! Number of element to be refined
    integer, intent(in) :: iel

    ! Identifier for element marker
    integer, intent(in) :: imarker

    ! Callback function
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
    integer :: ipos
    integer :: nel0,e1,e2,e3,e4,e5,e6,e7,e8
    integer :: i1,i2,i3,i4,i5,i6
    integer :: loc1,loc2,loc3,loc4
    
    ! Find local position of edge to be subdivided
    select case(imarker)
    case(MARK_REF_QUAD4TRIA_12)
      loc1=1; loc2=2; loc3=3; loc4=4

    case(MARK_REF_QUAD4TRIA_23)
      loc1=2; loc2=3; loc3=4; loc4=1

    case(MARK_REF_QUAD4TRIA_34)
      loc1=3; loc2=4; loc3=1; loc4=2

    case(MARK_REF_QUAD4TRIA_14)
      loc1=4; loc2=1; loc3=2; loc4=3

    case DEFAULT
      call output_line('Invalid element marker!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad4Tria')
      call sys_halt()
    end select

    ! Store vertex- and element-values of the current element
    i1 = rhadapt%p_IverticesAtElement(loc1, iel)
    i2 = rhadapt%p_IverticesAtElement(loc2, iel)
    i3 = rhadapt%p_IverticesAtElement(loc3, iel)
    i4 = rhadapt%p_IverticesAtElement(loc4, iel)

    e1 = rhadapt%p_IneighboursAtElement(loc1, iel)
    e2 = rhadapt%p_IneighboursAtElement(loc2, iel)
    e3 = rhadapt%p_IneighboursAtElement(loc3, iel)
    e4 = rhadapt%p_IneighboursAtElement(loc4, iel)

    e5 = rhadapt%p_ImidneighboursAtElement(loc1, iel)
    e6 = rhadapt%p_ImidneighboursAtElement(loc2, iel)
    e7 = rhadapt%p_ImidneighboursAtElement(loc3, iel)
    e8 = rhadapt%p_ImidneighboursAtElement(loc4, iel)
    
    ! Store total number of elements before refinement
    nel0 = rhadapt%NEL
    

    ! Add new vertices I5 and I6 at the midpoint of edges (I1,I2) and
    ! (I2,I3), respectively.
    call add_vertex2D(rhadapt, i1, i2, e1, i5,&
                      rcollection, fcb_hadaptCallback)
    call add_vertex2D(rhadapt, i2, i3, e2, i6,&
                      rcollection, fcb_hadaptCallback)

    
    ! Replace element IEL and add three new elements NEL0+1,NEL0+2 and NEL0+3
    call replace_element2D(rhadapt, iel, i1, i5, i4,&
                           e1, nel0+3, e4, e1, nel0+3, e8)
    call add_element2D(rhadapt, i2, i6, i5,&
                       e2, nel0+3, e5, e2, nel0+3, e5)
    call add_element2D(rhadapt, i3, i4, i6,&
                       e3, nel0+3, e6, e7, nel0+3, e6)
    call add_element2D(rhadapt, i4, i5, i6,&
                       iel, nel0+1, nel0+2, iel, nel0+1, nel0+2)
    

    ! Update list of neighboring elements
    call update_ElementNeighbors2D(rhadapt, e1, e5, iel, nel0+1, iel)
    call update_ElementNeighbors2D(rhadapt, e2, e6, iel, nel0+2, nel0+1)
    call update_ElementNeighbors2D(rhadapt, e3, e7, iel, nel0+2, nel0+2)

    
    ! Update list of elements meeting at vertices
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i2, iel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad4Tria')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i3, iel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad4Tria')
      call sys_halt()
    end if
    
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i2, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i3, nel0+2, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i4, nel0+2, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i4, nel0+3, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, iel,    ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, nel0+3, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i6, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i6, nel0+2, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i6, nel0+3, ipos)


    ! Adjust number of elements
    rhadapt%InelOfType(TRIA_NVETRI2D)  = rhadapt%InelOfType(TRIA_NVETRI2D)+1
    rhadapt%InelOfType(TRIA_NVEQUAD2D) = rhadapt%InelOfType(TRIA_NVEQUAD2D)-1


    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:14) = (/i1, i2, i3, i4, i5, i6,&
                                         e1, e2, e3, e4, e5, e6, e7, e8/)
      call fcb_hadaptCallback(HADAPT_OPR_REF_QUAD4TRIA, rcollection)
    end if

  end subroutine refine_Quad4Tria
    
  ! ***************************************************************************

!<subroutine>
  
  subroutine refine_Quad4Quad(rhadapt, iel,&
                              rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one quadrilateral element four similar
    ! quadrilateral elements.
    !
    ! In the illustration, i1-i9 denote the vertices whereby i5-i9 are the new
    ! vertices. Moreover, (e1)-(e8) stand for the element numbers which are
    ! adjacent and mid-adjacent to the current element.
    ! The new elements are assigned the total number of elements currently present
    ! in the triangulation increased by one, two and three, respectively.
    !
    ! <verb>
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7)      (e3) i3          i4 (e7) i7  (e3) i3
    !      +----------------+            +-------+-------+
    !      |                |            |*      |      *|
    ! (e4) |                | (e6)  (e4) | nel+3 | nel+2 | (e6)
    !      |                |            |       |       |
    !      |       iel      |     ->   i8+-------+-------+i6
    !      |                |            |     i9|       |
    ! (e8) |                | (e2)  (e8) | iel   | nel+1 | (e2)
    !      |*               |            |*      |      *|
    !      +----------------+            +-------+-------+
    !     i1 (e1)      (e5) i2          i1 (e1) i5  (e5) i2
    ! </verb>
    !
!</description>

!<input>
    ! Number of element to be refined
    integer, intent(in) :: iel

    ! Callback function
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
    integer :: ipos
    integer :: nel0,e1,e2,e3,e4,e5,e6,e7,e8
    integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9
    
    ! Store vertex- and element-values of the current element
    i1 = rhadapt%p_IverticesAtElement(1, iel)
    i2 = rhadapt%p_IverticesAtElement(2, iel)
    i3 = rhadapt%p_IverticesAtElement(3, iel)
    i4 = rhadapt%p_IverticesAtElement(4, iel)

    e1 = rhadapt%p_IneighboursAtElement(1, iel)
    e2 = rhadapt%p_IneighboursAtElement(2, iel)
    e3 = rhadapt%p_IneighboursAtElement(3, iel)
    e4 = rhadapt%p_IneighboursAtElement(4, iel)

    e5 = rhadapt%p_ImidneighboursAtElement(1, iel)
    e6 = rhadapt%p_ImidneighboursAtElement(2, iel)
    e7 = rhadapt%p_ImidneighboursAtElement(3, iel)
    e8 = rhadapt%p_ImidneighboursAtElement(4, iel)

    ! Store total number of elements before refinement
    nel0 = rhadapt%NEL
    
    
    ! Add five new vertices I5,I6,I7,I8 and I9 at the midpoint of edges
    ! (I1,I2), (I2,I3), (I3,I4) and (I1,I4) and at the center of element IEL
    call add_vertex2D(rhadapt, i1, i2, e1, i5,&
                      rcollection, fcb_hadaptCallback)
    call add_vertex2D(rhadapt, i2, i3, e2, i6,&
                      rcollection, fcb_hadaptCallback)
    call add_vertex2D(rhadapt, i3, i4, e3, i7,&
                      rcollection, fcb_hadaptCallback)
    call add_vertex2D(rhadapt, i4, i1, e4, i8,&
                      rcollection, fcb_hadaptCallback)
    call add_vertex2D(rhadapt, i1, i2, i3, i4, i9,&
                      rcollection, fcb_hadaptCallback)
    

    ! Replace element IEL and add three new elements NEL0+1,NEL0+2 and NEL0+3
    call replace_element2D(rhadapt, iel, i1, i5, i9, i8,&
                           e1, nel0+1, nel0+3, e8, e1, nel0+1, nel0+3, e8)
    call add_element2D(rhadapt, i2, i6, i9, i5,&
                       e2, nel0+2, iel, e5, e2, nel0+2, iel, e5)
    call add_element2D(rhadapt, i3, i7, i9, i6,&
                       e3, nel0+3, nel0+1, e6, e3, nel0+3, nel0+1, e6)
    call add_element2D(rhadapt,i4, i8, i9, i7,&
                       e4, iel, nel0+2, e7, e4, iel, nel0+2, e7)

    
    ! Update list of neighboring elements
    call update_ElementNeighbors2D(rhadapt, e1, e5, iel, nel0+1, iel)
    call update_ElementNeighbors2D(rhadapt, e2, e6, iel, nel0+2, nel0+1)
    call update_ElementNeighbors2D(rhadapt, e3, e7, iel, nel0+3, nel0+2)
    call update_ElementNeighbors2D(rhadapt, e4, e8, iel, iel, nel0+3)

        
    ! Update list of elements meeting at vertices
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i2, iel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad4Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i3, iel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad4Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i4, iel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad4Quad')
      call sys_halt()
    end if
    
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i2, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i3, nel0+2, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i4, nel0+3, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, iel,    ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i6, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i6, nel0+2, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i7, nel0+2, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i7, nel0+3, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i8, iel,    ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i8, nel0+3, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i9, iel,    ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i9, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i9, nel0+2, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i9, nel0+3, ipos)


    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:17) = (/i1, i2, i3, i4, i5, i6, i7, i8, i9,&
                                         e1, e2, e3, e4, e5, e6, e7, e8/)
      call fcb_hadaptCallback(HADAPT_OPR_REF_QUAD4QUAD, rcollection)
    end if

  end subroutine refine_Quad4Quad
  
  ! ***************************************************************************

!<subroutine>

  subroutine convert_Tria2Tria(rhadapt, iel, jel,&
                               rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine combines two neighboring triangles into one triangle
    ! and performs regular refinement into four similar triangles afterwards.
    ! The local orientation of both elements can be uniquely determined
    ! from the elemental states so that the first node of each triangle
    ! is located at the midpoint of the bisected edge.
    !
    ! <verb>
    !    initial triangle           subdivided triangle
    !
    !            i3                         i3
    !            +                          +
    !           /|\                        /X\
    !     (e3) / | \ (e5)            (e3) /   \ (e5)
    !         /  |  \                    /nel+1\
    !        /   |   \        ->      i6+-------+i5
    !       /    |    \                / \nel+2/ \
    ! (e6) / iel | jel \ (e2)    (e6) /   \   /   \ (e2)
    !     /X     |     X\            /X iel\X/jel X\
    !    +-------+-------+          +-------+-------+
    !   i1(e1,e7)i4(e4,e8)i2       i1(e1,e7)i4(e4,e8)i2
    ! </verb>
    !
!</description>

!<input>
    ! Number of first (left) element
    integer, intent(in) :: iel

    ! Number of second (right) element
    integer, intent(in) :: jel

    ! Callback function
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
    integer :: ipos
    integer :: i1,i2,i3,i4,i5,i6
    integer :: nel0,e1,e2,e3,e4,e5,e6,e7,e8

    ! Find local positions of elements IEL and JEL
    i1 = rhadapt%p_IverticesAtElement(1, iel)
    i4 = rhadapt%p_IverticesAtElement(2, iel)
    i3 = rhadapt%p_IverticesAtElement(3, iel)
    i2 = rhadapt%p_IverticesAtElement(1, jel)
    
    e1 = rhadapt%p_IneighboursAtElement(1, iel)
    e3 = rhadapt%p_IneighboursAtElement(3, iel)
    e2 = rhadapt%p_IneighboursAtElement(1, jel)
    e4 = rhadapt%p_IneighboursAtElement(3, jel)
    
    e7 = rhadapt%p_ImidneighboursAtElement(1, iel)
    e6 = rhadapt%p_ImidneighboursAtElement(3, iel)
    e5 = rhadapt%p_ImidneighboursAtElement(1, jel)
    e8 = rhadapt%p_ImidneighboursAtElement(3, jel)
    
    ! Store total number of elements before conversion
    nel0 = rhadapt%NEL

    
    ! Add two new vertices I5 and I6 at the midpoint of edges (I2,I3)
    ! and (I1,I3), respectively.
    call add_vertex2D(rhadapt, i2, i3, e2, i5,&
                      rcollection, fcb_hadaptCallback)
    call add_vertex2D(rhadapt, i3, i1, e3, i6,&
                      rcollection, fcb_hadaptCallback)


    ! Replace elements IEL and JEL and add two new elements NEL0+1 and NEL0+2
    call replace_element2D(rhadapt, iel, i1, i4, i6,&
                           e1, nel0+2, e6, e7, nel0+2, e6)
    call replace_element2D(rhadapt, jel, i2, i5, i4,&
                           e2, nel0+2, e4, e2, nel0+2, e8)
    call add_element2D(rhadapt, i3, i6, i5,&
                       e3, nel0+2, e5, e3, nel0+2, e5)
    call add_element2D(rhadapt, i4, i5, i6,&
                       jel, nel0+1, iel, jel, nel0+1, iel)

    
    ! Update list of neighboring elements
    call update_ElementNeighbors2D(rhadapt, e2, e5, jel, nel0+1, jel)
    call update_ElementNeighbors2D(rhadapt, e3, e6, iel, iel, nel0+1)


    ! Update list of elements meeting at vertices
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i3, iel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'convert_Tria2Tria')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i3, jel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'convert_Tria2Tria')
      call sys_halt()
    end if

    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i3, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i4, nel0+2, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, jel,    ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, nel0+2, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i6, iel,    ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i6, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i6, nel0+2, ipos)


    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:14) = (/i1, i2, i3, i4, i5, i6,&
                                         e1, e2, e3, e4, e5, e6, e7, e8/)
      call fcb_hadaptCallback(HADAPT_OPR_CVT_TRIA2TRIA, rcollection)
    end if

  end subroutine convert_Tria2Tria

  ! ***************************************************************************

!<subroutine>

  subroutine convert_Quad2Quad(rhadapt, iel, jel,&
                               rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine combines two neighboring quadrilaterals into one
    ! and performs regular refinement into four similar quadrilaterals.
    ! The local orientation of both elements can be uniquely determined
    ! from the elemental states so that the first node of each
    ! Depending on the given state of the element, the corresponding
    ! neighboring element is given explicitly.
    !
    ! <verb>
    !     initial quadrilateral      subdivided quadrilateral
    !
    !     i4(e7,e3)i7(f5,f1)i3         i4(e7,e3)i7(f5,f1)i3
    !      +-------+-------+            +-------+-------+
    !      |       |      *|            |*      |      *|
    ! (e4) |       |       | (f8)  (e4) | nel+2 | nel+1 | (f8)
    !      |       |       |            |       |       |
    !      |  iel  | jel   |     ->   i8+-------+-------+i6
    !      |       |       |            |       |i9     |
    ! (e8) |       |       | (f4)  (e8) |  iel  |  jel  | (f4)
    !      |*      |       |            |*      |      *|
    !      +-------+-------+            +-------+-------+
    !     i1(e1,e5)i5(f3,f7)i2         i1(e1,e5)i5(f3,f7)i2
    ! </verb>
    !
!</description>

!<input>
    ! Number of first element
    integer, intent(in) :: iel

    ! Number of second element
    integer, intent(in) :: jel
    
    ! Callback function
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
    integer :: ipos
    integer :: nel0,e1,e3,e4,e5,e7,e8,f1,f3,f4,f5,f7,f8
    integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9

    ! Find local positions of elements IEL and JEL
    i1 = rhadapt%p_IverticesAtElement(1, iel)
    i5 = rhadapt%p_IverticesAtElement(2, iel)
    i7 = rhadapt%p_IverticesAtElement(3, iel)
    i4 = rhadapt%p_IverticesAtElement(4, iel)
    i3 = rhadapt%p_IverticesAtElement(1, jel)
    i2 = rhadapt%p_IverticesAtElement(4, jel)

    e1 = rhadapt%p_IneighboursAtElement(1, iel)
    e3 = rhadapt%p_IneighboursAtElement(3, iel)
    e4 = rhadapt%p_IneighboursAtElement(4, iel)
    f1 = rhadapt%p_IneighboursAtElement(1, jel)
    f3 = rhadapt%p_IneighboursAtElement(3, jel)
    f4 = rhadapt%p_IneighboursAtElement(4, jel)
    
    e5 = rhadapt%p_ImidneighboursAtElement(1, iel)
    e7 = rhadapt%p_ImidneighboursAtElement(3, iel)
    e8 = rhadapt%p_ImidneighboursAtElement(4, iel)
    f5 = rhadapt%p_ImidneighboursAtElement(1, jel)
    f7 = rhadapt%p_ImidneighboursAtElement(3, jel)
    f8 = rhadapt%p_ImidneighboursAtElement(4, jel)

    ! Store total number of elements before conversion
    nel0 = rhadapt%NEL

    
    ! Add two new vertices I6, I8, and I9 at the midpoint of edges (I2,I3),
    ! (I1,I4) and (I1,I2), respectively.
    call add_vertex2D(rhadapt, i2, i3, f4, i6,&
                      rcollection, fcb_hadaptCallback)
    call add_vertex2D(rhadapt, i4, i1, e4, i8,&
                      rcollection, fcb_hadaptCallback)
    call add_vertex2D(rhadapt, i1, i2, i3, i4, i9,&
                      rcollection, fcb_hadaptCallback)


    ! Replace element IEL and JEL and add two new elements NEL0+1 and NEL0+2
    call replace_element2D(rhadapt, iel, i1, i5, i9, i8,&
                           e1, jel, nel0+2, e8, e5, jel, nel0+2, e8)
    call replace_element2D(rhadapt, jel, i2, i6, i9, i5,&
                           f4, nel0+1, iel, f3, f4, nel0+1, iel, f7)
    call add_element2D(rhadapt, i3, i7, i9, i6,&
                       f1, nel0+2, jel, f8, f5, nel0+2, jel, f8)
    call add_element2D(rhadapt, i4, i8, i9, i7,&
                       e4, iel, nel0+1, e3, e4, iel, nel0+1, e7)


    ! Update list of neighboring elements
    call update_ElementNeighbors2D(rhadapt, f4, f8, jel, nel0+1, jel)
    call update_ElementNeighbors2D(rhadapt, e4, e8, iel, iel, nel0+2)
    call update_ElementNeighbors2D(rhadapt, f1, f5, jel, nel0+1, nel0+1)
    call update_ElementNeighbors2D(rhadapt, e3, e7, iel, nel0+2, nel0+2)


    ! Update list of elements meeting at vertices
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i3, jel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad2Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i4, iel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad2Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i7, iel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad2Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i7, jel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad2Quad')
      call sys_halt()
    end if

    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i3, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i4, nel0+2, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i6, jel,    ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i6, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i7, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i7, nel0+2, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i8, iel,    ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i8, nel0+2, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i9, iel,    ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i9, jel,    ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i9, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i9, nel0+2, ipos)


    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:17) = (/i1, i2, i3, i4, i5, i6, i7, i8, i9,&
                                         e1, f4, f1, e4, f3, f8, e3, e8/)
      call fcb_hadaptCallback(HADAPT_OPR_CVT_QUAD2QUAD, rcollection)
    end if

  end subroutine convert_Quad2Quad

  ! ***************************************************************************

!<subroutine>

  subroutine convert_Quad3Tria(rhadapt, iel1, iel2, iel3,&
                               rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine combines three neighboring triangles which result
    ! from a Quad4Tria refinement into one quadrilateral and performs
    ! regular refinement into four quadrilaterals afterwards.
    ! This subroutine is based on the convention that IEL1 denotes the
    ! left element, IEL2 denots the right element and IEL3 stands for
    ! the triangle which connects IEL1 and IEL2.
    !
    ! <verb>
    !   initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7)     (e3) i2          i4 (e7)  i7 (e3)  i3
    !      +---------------+            +-------+-------+
    !      |\             /|            |*      |      *|
    ! (e4) | \           / | (e6)  (e4) | nel+1 | iel3  | (e6)
    !      |  \   iel3  /  |            |       |       |
    !      |   \       /   |     ->   i8+-------+-------+i6
    !      |    \     /    |            |       |i9     |
    ! (e8) | iel1\   /iel2 | (e2)  (e8) | iel1  | iel2  | (e2)
    !      |*     \X/     *|            |*      |      *|
    !      +-------+-------+            +-------+-------+
    !     i1(e1,e9)i5(e5,e10)i2        i1(e1,e9)i5(e5,e10)i2
    ! </verb>
    !
!</description>

!<input>
    ! Number of first triangle
    integer, intent(in) :: iel1

    ! Number of second triangle
    integer, intent(in) :: iel2
    
    ! Number of third triangle
    integer, intent(in) :: iel3

    ! Callback function
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
    integer :: ipos
    integer :: nel0,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10
    integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9

    ! Get local data from elements IEL1, IEL2 and IEL3
    i1 = rhadapt%p_IverticesAtElement(1, iel1)
    i5 = rhadapt%p_IverticesAtElement(2, iel1)
    i4 = rhadapt%p_IverticesAtElement(3, iel1)
    i2 = rhadapt%p_IverticesAtElement(1, iel2)
    i3 = rhadapt%p_IverticesAtElement(2, iel2)

    e1 = rhadapt%p_IneighboursAtElement(1, iel1)
    e4 = rhadapt%p_IneighboursAtElement(3, iel1)
    e2 = rhadapt%p_IneighboursAtElement(1, iel2)
    e5 = rhadapt%p_IneighboursAtElement(3, iel2)
    e3 = rhadapt%p_IneighboursAtElement(2, iel3)
    
    e9  = rhadapt%p_ImidneighboursAtElement(1, iel1)
    e8  = rhadapt%p_ImidneighboursAtElement(3, iel1)
    e6  = rhadapt%p_ImidneighboursAtElement(1, iel2)
    e10 = rhadapt%p_ImidneighboursAtElement(3, iel2)
    e7  = rhadapt%p_ImidneighboursAtElement(2, iel3)
    
    ! Store total number of elements before conversion
    nel0 = rhadapt%NEL

    
    ! Add four new vertices I6,I7,I8 and I9 at the midpoint of edges
    ! (I2,I3), (I3,I4) and (I1,I4) and at the center of element IEL
    call add_vertex2D(rhadapt, i2, i3, e2, i6,&
                      rcollection, fcb_hadaptCallback)
    call add_vertex2D(rhadapt, i3, i4, e3, i7,&
                      rcollection, fcb_hadaptCallback)
    call add_vertex2D(rhadapt, i4, i1, e4, i8,&
                      rcollection, fcb_hadaptCallback)
    call add_vertex2D(rhadapt, i1, i2, i3, i4, i9,&
                      rcollection, fcb_hadaptCallback)

    
    ! Replace elements IEL1, IEL2 and IEL3 and add one new element
    call replace_element2D(rhadapt, iel1, i1, i5, i9, i8,&
                           e1, iel2, nel0+1, e8, e9, iel2, nel0+1, e8)
    call replace_element2D(rhadapt, iel2, i2, i6, i9, i5,&
                           e2, iel3, iel1, e5, e2, iel3, iel1, e10)
    call replace_element2D(rhadapt, iel3, i3, i7, i9, i6,&
                           e3, nel0+1, iel2, e6, e3, nel0+1, iel2, e6)
    call add_element2D(rhadapt, i4, i8, i9, i7,&
                       e4, iel1, iel3, e7, e4, iel1, iel3, e7)


    ! Update element neighbors
    call update_ElementNeighbors2D(rhadapt, e2, e6, iel2, iel3, iel2)
    call update_ElementNeighbors2D(rhadapt, e3, e7, iel3, nel0+1, iel3)
    call update_ElementNeighbors2D(rhadapt, e4, e8, iel1, iel1, nel0+1)


    ! Update list of elements meeting at vertices
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i3, iel2) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad3Tria')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i4, iel1) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad3Tria')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i4, iel3) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad3Tria')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i5, iel3) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad3Tria')
      call sys_halt()
    end if
    
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i4, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i6, iel2,   ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i6, iel3,   ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i7, iel3,   ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i7, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i8, nel0+1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i8, iel1,   ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i9, iel1,   ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i9, iel2,   ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i9, iel3,   ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i9, nel0+1, ipos)

    ! Finally, adjust numbers of triangles/quadrilaterals
    rhadapt%InelOfType(TRIA_NVETRI2D)  = rhadapt%InelOfType(TRIA_NVETRI2D)-3
    rhadapt%InelOfType(TRIA_NVEQUAD2D) = rhadapt%InelOfType(TRIA_NVEQUAD2D)+3


    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:17) = (/i1, i2, i3, i4, i5, i6, i7, i8, i9,&
                                         e1, e2, e3, e4, e5, e6, e7, e8/)
      call fcb_hadaptCallback(HADAPT_OPR_CVT_QUAD3TRIA, rcollection)
    end if

  end subroutine convert_Quad3Tria

  ! ***************************************************************************

!<subroutine>
  
  subroutine convert_Quad4Tria(rhadapt, iel1, iel2, iel3, iel4,&
                               rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine combines four neighboring triangles which result
    ! from a Quad4Tria refinement into one quadrilateral and performs
    ! regular refinement into four quadrilaterals afterwards.
    ! This subroutine is based on the convention, that all four elements
    ! are given in couterclockwise order. More precisely, IEL2 and IEL3
    ! make up the inner "diamond" of the refinement, whereas IEL1 and IEL4
    ! are the right and left outer triangles, respectively.
    !
    ! <verb>
    !    initial quadrilateral         subdivided quadrilateral
    !
    !     i4 (e7)     (e3) i3           i4 (e7)  i7 (e3) i3
    !      +---------------+             +-------+-------+
    !      |\\\\    iel3  *| (e12)       |*      |      *| (e12)
    ! (e4) | \*  \\\       | (e6)   (e4) | iel4  | iel3  | (e6)
    !      |  \      \\\   |             |       |       |
    !      |   \  iel4   \\+i6    ->   i8+-------+-------+i6
    !      |    \        / |             |       |i9     |
    ! (e8) | iel1\     /   | (e11)  (e8) | iel1  | iel2  | (e11)
    !      |*     \  /iel2*| (e2)        |*      |      *| (e2)
    !      +-------+-------+             +-------+-------+
    !     i1(e1,e9)i5(e5,e10)i2         i1(e1,e9)i5(e5,e10)i2
    ! </verb>
    !
!</description>

!<input>
    ! Number of first triangle
    integer, intent(in) :: iel1

    ! Number of second triangle
    integer, intent(in) :: iel2

    ! Number of third triangle
    integer, intent(in) :: iel3

    ! Number of fourth triangle
    integer, intent(in) :: iel4

    ! Callback function
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
    integer :: ipos
    integer :: e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12
    integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9
    
    ! Get local data from elements IEL1, IEL2 and IEL3
    i1 = rhadapt%p_IverticesAtElement(1, iel1)
    i5 = rhadapt%p_IverticesAtElement(2, iel1)
    i4 = rhadapt%p_IverticesAtElement(3, iel1)
    i2 = rhadapt%p_IverticesAtElement(1, iel2)
    i6 = rhadapt%p_IverticesAtElement(2, iel2)
    i3 = rhadapt%p_IverticesAtElement(1, iel3)

    e1 = rhadapt%p_IneighboursAtElement(1, iel1)
    e4 = rhadapt%p_IneighboursAtElement(3, iel1)
    e2 = rhadapt%p_IneighboursAtElement(1, iel2)
    e5 = rhadapt%p_IneighboursAtElement(3, iel2)
    e3 = rhadapt%p_IneighboursAtElement(1, iel3)
    e6 = rhadapt%p_IneighboursAtElement(3, iel3)

    e9  = rhadapt%p_ImidneighboursAtElement(1, iel1)
    e8  = rhadapt%p_ImidneighboursAtElement(3, iel1)
    e11 = rhadapt%p_ImidneighboursAtElement(1, iel2)
    e10 = rhadapt%p_ImidneighboursAtElement(3, iel2)
    e7  = rhadapt%p_ImidneighboursAtElement(1, iel3)
    e12 = rhadapt%p_ImidneighboursAtElement(3, iel3)
    

    ! Add three new vertices I7, I8 and I9 at the midpoint of edges
    ! (I3,I4) and (I1,I4) and at the center of element IEL
    call add_vertex2D(rhadapt, i3, i4, e3, i7,&
                      rcollection, fcb_hadaptCallback)
    call add_vertex2D(rhadapt, i4, i1, e4, i8,&
                      rcollection, fcb_hadaptCallback)
    call add_vertex2D(rhadapt, i1, i2, i3, i4, i9,&
                      rcollection, fcb_hadaptCallback)


    ! Replace all four elements
    call replace_element2D(rhadapt, iel1, i1, i5, i9, i8,&
                           e1, iel2, iel4, e8, e9, iel2, iel4, e8)
    call replace_element2D(rhadapt, iel2, i2, i6, i9, i5,&
                           e2, iel3, iel1, e5, e11, iel3, iel1, e10)
    call replace_element2D(rhadapt, iel3, i3, i7, i9, i6,&
                           e3, iel4, iel2, e6, e3, iel4, iel2, e12)
    call replace_element2D(rhadapt, iel4, i4, i8, i9, i7,&
                           e4, iel1, iel3, e7, e4, iel1, iel3, e7)


    ! Update element neighbors
    call update_ElementNeighbors2D(rhadapt, e3, e7, iel3, iel4, iel3)
    call update_ElementNeighbors2D(rhadapt, e4, e8, iel1, iel1, iel4)


    ! Update list of elements meeting at vertices
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i5, iel4) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad4Tria')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i6, iel4) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad4Tria')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i4, iel1) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad4Tria')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i4, iel3) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad4Tria')
      call sys_halt()
    end if
    
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i7, iel3, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i7, iel4, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i8, iel1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i8, iel4, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i9, iel1, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i9, iel2, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i9, iel3, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i9, iel4, ipos)


    ! Finally, adjust numbers of triangles/quadrilaterals
    rhadapt%InelOfType(TRIA_NVETRI2D)  = rhadapt%InelOfType(TRIA_NVETRI2D)-4
    rhadapt%InelOfType(TRIA_NVEQUAD2D) = rhadapt%InelOfType(TRIA_NVEQUAD2D)+4


    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:17) = (/i1, i2, i3, i4, i5, i6, i7, i8, i9,&
                                         e1, e2, e3, e4, e5, e6, e7, e8/)
      call fcb_hadaptCallback(HADAPT_OPR_CVT_QUAD4TRIA, rcollection)
    end if

  end subroutine convert_Quad4Tria

  ! ***************************************************************************

!<subroutine>

  subroutine coarsen_2Tria1Tria(rhadapt, Imarker, iel, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine combines two triangles resulting from a 1-tria : 2-tria
    ! green refinement into the original macro triangle. Note that element IEL
    ! must(!) be the left triangle which is combined with its right neighbour.
    ! This is not checked in the routine since this check is already performed
    ! in the marking routine.
    !
    ! Due to the numbering convention used in the refinement procedure
    ! the "second" element neighbor is the other green triangle that is removed.
    ! In order to restore the original orientation the outcome of the resulting
    ! triangle is considered. If the macro element belongs to the inital
    ! triangulation of if it is an inner red triangle, an arbitrary orientation
    ! is adopted. If the macro element is an outer red triangle, then the states
    ! 2 and 8 are trapped and converted into state 4.
    !
    ! <verb>
    !     initial triangle           subdivided triangle
    !
    !            i3                        i3
    !            +                         +
    !           /|\                       / \
    !     (e3) / | \ (e5)           (e3) /   \ (e5)
    !         /  |  \                   /     \
    !        /   |   \        ->       /       \
    !       /    |    \               /   jel   \
    ! (e6) / iel | iel1\ (e2)   (e6) /           \ (e2)
    !     /X     |     X\           /X            \
    !    +-------+-------+         +---------------+
    !    i1 (e1) i4 (e4) i2        i1 (e1)    (e4) i2
    ! </verb>
    !
!</description>

!<input>
    ! Element number of the inner red triangle
    integer, intent(in) :: iel

    ! Callback function
    include 'intf_hadaptcallback.inc'
    optional :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt

    ! Marker array
    integer, dimension(:), intent(inout) :: Imarker

    ! OPTIONAL: Collection
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(3) :: ImacroVertices
    integer, dimension(TRIA_NVETRI2D) :: IvertexAge
    integer :: ipos,istate
    integer :: iel1,e1,e2,e3,e4,e5,e6,jel,ielRemove,ielReplace
    integer :: i1,i2,i3,i4
    
          
    ! Get right-adjacent green elements
    iel1 = rhadapt%p_IneighboursAtElement(2, iel)
    
    ! Determine element with smaller element number
    if (iel .lt. iel1) then
      jel=iel; ielRemove=iel1
    else
      jel=iel1; ielRemove=iel
    end if
    
    ! Store vertex- and element values of the two elements
    i1 = rhadapt%p_IverticesAtElement(1, iel)
    i4 = rhadapt%p_IverticesAtElement(2, iel)
    i3 = rhadapt%p_IverticesAtElement(3, iel)
    i2 = rhadapt%p_IverticesAtElement(1, iel1)
    
    e1 = rhadapt%p_IneighboursAtElement(1, iel)
    e3 = rhadapt%p_IneighboursAtElement(3, iel)
    e2 = rhadapt%p_IneighboursAtElement(1, iel1)
    e4 = rhadapt%p_IneighboursAtElement(3, iel1)
    
    e5 = rhadapt%p_ImidneighboursAtElement(1, iel1)
    e6 = rhadapt%p_ImidneighboursAtElement(3, iel)

    
    ! Update list of neighboring elements
    call update_ElementNeighbors2D(rhadapt, e1, e4, iel1, iel, jel, jel)
    if (iel .lt. iel1) then
      call update_ElementNeighbors2D(rhadapt, e2, e5, iel1, jel, jel)
    else
      call update_ElementNeighbors2D(rhadapt, e3, e6, iel, jel, jel)
    end if

    
    ! The resulting triangle will possess one of the states STATE_TRIA_OUTERINNERx, whereby
    ! x can be blank, 1 or 2. Due to our refinement convention, the two states x=1 and x=2
    ! should not appear, that is, local numbering of the resulting triangle starts at the
    ! vertex which is opposite to the inner red triangle. To this end, we check the state of
    ! the provisional triangle (I1,I2,I3) and transform the orientation accordingly.
    ImacroVertices = (/i1, i2, i3/)
    IvertexAge = rhadapt%p_IvertexAge(ImacroVertices)
    istate     = redgreen_getstateTria(IvertexAge)
    
    select case(istate)
    case(STATE_TRIA_OUTERINNER,&
         STATE_TRIA_ROOT,&
         STATE_TRIA_REDINNER)
      ! Update element JEL = (I1,I2,I3)
      call replace_element2D(rhadapt, jel, i1, i2, i3,&
                             e1, e2, e3, e4, e5, e6)
      
    case(STATE_TRIA_OUTERINNER2)
      ! Update element JEL = (I2,I3,I1)
      call replace_element2D(rhadapt, jel, i2, i3, i1,&
                             e2, e3, e1, e5, e6, e4)
      
    case(STATE_TRIA_OUTERINNER1)
      ! Update element JEL = (I3,I1,I2)
      call replace_element2D(rhadapt,jel, i3, i1, i2,&
                             e3, e1, e2, e6, e4, e5)
      
    case DEFAULT
      call output_line('Invalid state of resulting triangle!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Tria1Tria')
      call sys_halt()
    end select
    

    ! Delete element IEL or IEL1 depending on which element has smaller
    ! element number, that is, is not equal to JEL
    call remove_element2D(rhadapt, ielRemove, ielReplace)
    if (ielReplace .ne. 0) then
      call update_AllElementNeighbors2D(rhadapt, ielReplace, ielRemove)
      if (ielReplace .lt. iel) then
        Imarker(ielRemove)  = Imarker(ielReplace)
        Imarker(ielReplace) = MARK_ASIS
      end if
    end if
    
    
    ! Update list of elements meeting at vertices
    if (ielRemove .eq. iel1) then
      
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i2, ielRemove) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Tria1Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i3, ielRemove) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Tria1Tria')
        call sys_halt()
      end if

      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i2, jel, ipos)
      
    else
      
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i1, ielRemove) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Tria1Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i3, ielRemove) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Tria1Tria')
        call sys_halt()
      end if

      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i1, jel, ipos)
      
    end if
    
    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:8) = (/i1, i2, i3, i4,&
                                        e1, e2, e3, e4/)
      call fcb_hadaptCallback(HADAPT_OPR_CRS_2TRIA1TRIA, rcollection)
    end if

  end subroutine coarsen_2Tria1Tria

  ! ***************************************************************************

!<subroutine>

  subroutine coarsen_4Tria1Tria(rhadapt, Imarker, iel, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine combines four triangles resulting from a
    ! 1-tria : 4-tria refinement into the original macro triangle.
    ! By definition, iel is the number of the inner red triangle.
    !
    ! <verb>
    !    initial triangle           subdivided triangle
    !
    !            i3                        i3
    !            +                         +
    !           /X\                       / \
    !     (e3) /   \ (e5)           (e3) /   \ (e5)
    !         / iel3\                   /     \
    !      i6+-------+i5      ->       /       \
    !       / \ iel / \               /   jel   \
    ! (e6) /   \   /   \ (e2)   (e6) /           \ (e2)
    !     /Xiel1\X/iel2X\           /             \
    !    +-------+-------+         +---------------+
    !   i1 (e1)  i4 (e4) i2       i1  (e1)   (e4)  i2
    ! </verb>
    !
!</description>

!<input>
    ! Element number of the inner red triangle
    integer, intent(in) :: iel

    ! Callback function
    include 'intf_hadaptcallback.inc'
    optional :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
    
    ! Marker array
    integer, dimension(:), intent(inout) :: Imarker

    ! OPTIONAL: Collection
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(4) :: IsortedElements
    integer, dimension(3) :: ImacroVertices
    integer, dimension(TRIA_NVETRI2D) :: IvertexAge
    integer :: ipos,istate
    integer :: iel1,iel2,iel3,e1,e2,e3,e4,e5,e6,jel,ielReplace
    integer :: i1,i2,i3,i4,i5,i6
    
    ! Retrieve patch of elements
    iel2 = rhadapt%p_IneighboursAtElement(1, iel)
    iel3 = rhadapt%p_IneighboursAtElement(2, iel)
    iel1 = rhadapt%p_IneighboursAtElement(3, iel)


    ! Store vertex- and element-values of the three neighboring elements
    i4 = rhadapt%p_IverticesAtElement(1, iel)
    i5 = rhadapt%p_IverticesAtElement(2, iel)
    i6 = rhadapt%p_IverticesAtElement(3, iel)
    i1 = rhadapt%p_IverticesAtElement(1, iel1)
    i2 = rhadapt%p_IverticesAtElement(1, iel2)
    i3 = rhadapt%p_IverticesAtElement(1, iel3)

    ! Store values of the elements adjacent to the resulting macro element
    e1 = rhadapt%p_IneighboursAtElement(1, iel1)
    e6 = rhadapt%p_IneighboursAtElement(3, iel1)
    e2 = rhadapt%p_IneighboursAtElement(1, iel2)
    e4 = rhadapt%p_IneighboursAtElement(3, iel2)
    e3 = rhadapt%p_IneighboursAtElement(1, iel3)
    e5 = rhadapt%p_IneighboursAtElement(3, iel3)


    ! Sort the four elements according to their number and
    ! determine the element with the smallest element number
    IsortedElements = (/iel,iel1,iel2,iel3/)
    call sort_I32(IsortedElements, SORT_INSERT)
    jel = IsortedElements(1)


    ! Update list of neighboring elements
    call update_ElementNeighbors2D(rhadapt, e1, e4, iel2, iel1, jel, jel)
    call update_ElementNeighbors2D(rhadapt, e2, e5, iel3, iel2, jel, jel)
    call update_ElementNeighbors2D(rhadapt, e3, e6, iel1, iel3, jel, jel)

    
    ! The resulting triangle will posses one of the states STATE_TRIA_REDINNER or
    ! STATE_TRIA_OUTERINNERx, whereby x can be blank, 1 or 2. Due to our refinement
    ! convention, the two states x=1 and x=2 should not appear, that is, local
    ! numbering of the resulting triangle starts at the vertex which is opposite
    ! to the inner red triangle. To this end, we check the state of the provisional
    ! triangle (I1,I2,I3) and transform the orientation accordingly.
    ImacroVertices = (/i1,i2,i3/)
    IvertexAge = rhadapt%p_IvertexAge(ImacroVertices)
    istate     = redgreen_getstateTria(IvertexAge)
    
    select case(istate)
    case(STATE_TRIA_OUTERINNER,&
         STATE_TRIA_ROOT,&
         STATE_TRIA_REDINNER)
      ! Update element JEL = (I1,I2,I3)
      call replace_element2D(rhadapt, jel, i1, i2, i3,&
                             e1, e2, e3, e4, e5, e6)
      
    case(STATE_TRIA_OUTERINNER1)
      ! Update element JEL = (I3,I1,I2)
      call replace_element2D(rhadapt, jel, i3, i1, i2,&
                             e3, e1, e2, e6, e4, e5)

    case(STATE_TRIA_OUTERINNER2)
      ! Update element JEL = (I2,I3,I1)
      call replace_element2D(rhadapt, jel, i2, i3, i1,&
                             e2, e3, e1, e5, e6, e4)

    case DEFAULT
      call output_line('Invalid state of resulting triangle!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Tria')
      call sys_halt()
    end select


    ! Delete elements IEL, IEL1, IEL2 and IEL3 depending on which
    ! element corresponds to element with minimum number JEL
    call remove_element2D(rhadapt, IsortedElements(4), ielReplace)
    if (ielReplace .ne. 0) then
      call update_AllElementNeighbors2D(rhadapt, ielReplace, IsortedElements(4))
      if (ielReplace .lt. iel) then
        Imarker(IsortedElements(4)) = Imarker(ielReplace)
        Imarker(ielReplace)         = MARK_ASIS
      end if
    end if

    call remove_element2D(rhadapt, IsortedElements(3), ielReplace)
    if (ielReplace .ne. 0) then
      call update_AllElementNeighbors2D(rhadapt, ielReplace, IsortedElements(3))
      if (ielReplace .lt. iel) then
        Imarker(IsortedElements(3)) = Imarker(ielReplace)
        Imarker(ielReplace)         = MARK_ASIS
      end if
    end if

    call remove_element2D(rhadapt, IsortedElements(2), ielReplace)
    if (ielReplace .ne. 0) then
      call update_AllElementNeighbors2D(rhadapt, ielReplace, IsortedElements(2))
      if (ielReplace .lt. iel) then
        Imarker(IsortedElements(2)) = Imarker(ielReplace)
        Imarker(ielReplace)         = MARK_ASIS
      end if
    end if


    ! Update list of elements meeting at vertices.
    ! Note that all elements are removed in the first step. Afterwards,
    ! element JEL is appended to the list of elements meeting at each vertex
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i1, iel1) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Tria')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i2, iel2) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Tria')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i3, iel3) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Tria')
      call sys_halt()
    end if

    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i1, jel, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i2, jel, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i3, jel, ipos)


    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:12) = (/i1, i2, i3, i4, i5, i6,&
                                         e1, e2, e3, e4, e5, e6/)
      call fcb_hadaptCallback(HADAPT_OPR_CRS_4TRIA1TRIA, rcollection)
    end if

  end subroutine coarsen_4Tria1Tria

  ! ***************************************************************************

!<subroutine>

  subroutine coarsen_4Tria2Tria(rhadapt, Imarker, iel,&
                                rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine combines four triangles resulting from a
    ! 1-tria : 4-tria refinement into two green triangle.
    ! The local position (node 1,2,3 of the interior element) of the midpoint
    ! vertex that should be kept is identified by the marker.
    ! By definition, iel is the number of the inner red triangle.
    ! Be warned, this numbering convention is not checked since the
    ! routine is only called internally and cannot be used from outside.
    !
    ! Due to the numbering convention used in the refinement procedure
    ! the "third" element neighbor of the inner triangle has the number
    ! of the original macro element which will be used for the first resulting
    ! triangle. The number of the second triangle will be the number of the
    ! "first" element neighbor of the inner triangle.
    !
    ! <verb>
    !    initial triangle           subdivided triangle
    !
    !            i3                        i3
    !            +                         +
    !           /X\                       /|\
    !    (e3)  /   \ (e5)           (e3) / | \ (e5)
    !         /iel3 \                   /  |  \
    !     i6 +-------+ i5     ->       /   |   \
    !       / \ iel / \               /    |    \
    ! (e6) /   \   /   \ (e2)   (e6) /     |     \ (e2)
    !     /Xiel1\X/iel2X\           /X jel1|jel2 X\
    !    +-------+-------+         +-------+-------+
    !    i1 (e1) i4 (e4) i2        i1 (e1) i4 (e4)  i2
    ! </verb>
    !
!</description>

!<input>
    ! Element number of the inner red triangle
    integer, intent(in) :: iel

    ! Callback function
    include 'intf_hadaptcallback.inc'
    optional :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt

    ! Marker array
    integer, dimension(:), intent(inout) :: Imarker

    ! OPTIONAL: Collection
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(4) :: IsortedElements
    integer :: ipos
    integer :: iel1,iel2,iel3,e1,e2,e3,e4,e5,e6,jel1,jel2,ielReplace
    integer :: i1,i2,i3,i4,i5,i6
    
    ! Retrieve patch of elements
    iel2 = rhadapt%p_IneighboursAtElement(1, iel)
    iel3 = rhadapt%p_IneighboursAtElement(2, iel)
    iel1 = rhadapt%p_IneighboursAtElement(3, iel)

    ! Store vertex- and element-values of the three neighboring elements
    i4 = rhadapt%p_IverticesAtElement(1, iel)
    i5 = rhadapt%p_IverticesAtElement(2, iel)
    i6 = rhadapt%p_IverticesAtElement(3, iel)
    i1 = rhadapt%p_IverticesAtElement(1, iel1)
    i2 = rhadapt%p_IverticesAtElement(1, iel2)
    i3 = rhadapt%p_IverticesAtElement(1, iel3)

    ! Store values of the elements adjacent to the resulting macro element
    e1 = rhadapt%p_IneighboursAtElement(1, iel1)
    e6 = rhadapt%p_IneighboursAtElement(3, iel1)
    e2 = rhadapt%p_IneighboursAtElement(1, iel2)
    e4 = rhadapt%p_IneighboursAtElement(3, iel2)
    e3 = rhadapt%p_IneighboursAtElement(1, iel3)
    e5 = rhadapt%p_IneighboursAtElement(3, iel3)


    ! Sort the four elements according to their number and
    ! determine the two elements withthe smallest element numbers
    IsortedElements=(/iel,iel1,iel2,iel3/)
    call sort_I32(IsortedElements, SORT_INSERT)
    jel1 = IsortedElements(1)
    jel2 = IsortedElements(2)


    ! Which midpoint vertex should be kept?
    select case(Imarker(iel))
    case(MARK_CRS_4TRIA2TRIA_1)
      ! Update list of neighboring elements
      call update_ElementNeighbors2D(rhadapt, e1, e4, iel2, iel1, jel2, jel1)
      call update_ElementNeighbors2D(rhadapt, e2, e5, iel3, iel2, jel2, jel2)
      call update_ElementNeighbors2D(rhadapt, e3, e6, iel1, iel3, jel1, jel1)


      ! Update elements JEL1 and JEL2
      call replace_element2D(rhadapt, jel1, i1, i4, i3,&
                             e1, jel2, e3, e1, jel2, e6)
      call replace_element2D(rhadapt, jel2, i2, i3, i4,&
                             e2, jel1, e4, e5, jel1, e4)
      

      ! Delete elements IEL, IEL1, IEL2 and IEL3 depending on which
      ! elements correspond to the two elements with smallest numbers
      call remove_element2D(rhadapt, IsortedElements(4), ielReplace)
      if (ielReplace .ne. 0) then
        call update_AllElementNeighbors2D(rhadapt, ielReplace, IsortedElements(4))
        if (ielReplace .lt. iel) then
          Imarker(IsortedElements(4)) = Imarker(ielReplace)
          Imarker(ielReplace)         = MARK_ASIS
        end if
      end if

      call remove_element2D(rhadapt, IsortedElements(3), ielReplace)
      if (ielReplace .ne. 0) then
        call update_AllElementNeighbors2D(rhadapt, ielReplace, IsortedElements(3))
        if (ielReplace .lt. iel) then
          Imarker(IsortedElements(3)) = Imarker(ielReplace)
          Imarker(ielReplace)         = MARK_ASIS
        end if
      end if

      
      ! Update list of elements meeting at vertices.
      ! Note that all elements are removed in the first step. Afterwards,
      ! element JEL is appended to the list of elements meeting at each vertex
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i1, iel1) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i2, iel2) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i3, iel3) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        call sys_halt()
      end if

      ! Note, this can be improved by checking against JEL1 and JEL2
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i4, iel) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i4, iel1) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i4, iel2) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        call sys_halt()
      end if

      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i1, jel1, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i2, jel2, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i3, jel1, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i3, jel2, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i4, jel1, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i4, jel2, ipos)


      ! Optionally, invoke callback function
      if (present(fcb_hadaptCallback) .and. present(rcollection)) then
        rcollection%IquickAccess(1:12) = (/i1, i2, i3, i4, i5, i6,&
                                           e1, e2, e3, e4, e5, e6/)
        call fcb_hadaptCallback(HADAPT_OPR_CRS_4TRIA2TRIA1, rcollection)
      end if


    case(MARK_CRS_4TRIA2TRIA_2)
      ! Update list of neighboring elements
      call update_ElementNeighbors2D(rhadapt, e1, e4, iel2, iel1, jel1, jel1)
      call update_ElementNeighbors2D(rhadapt, e2, e5, iel3, iel2, jel2, jel1)
      call update_ElementNeighbors2D(rhadapt, e3, e6, iel1, iel3, jel2, jel2)


      ! Update elements JEL1 and JEL2
      call replace_element2D(rhadapt, jel1, i2, i5, i1,&
                             e2, jel2, e1, e2, jel2, e4)
      call replace_element2D(rhadapt, jel2, i3, i1, i5,&
                             e3, jel1, e5, e6, jel1, e5)


      ! Delete elements IEL, IEL1, IEL2 and IEL3 depending on which
      ! elements correspond to the two elements with smallest numbers
      call remove_element2D(rhadapt, IsortedElements(4), ielReplace)
      if (ielReplace .ne. 0) then
        call update_AllElementNeighbors2D(rhadapt, ielReplace, IsortedElements(4))
        if (ielReplace .lt. iel) then
          Imarker(IsortedElements(4)) = Imarker(ielReplace)
          Imarker(ielReplace)         = MARK_ASIS
        end if
      end if

      call remove_element2D(rhadapt, IsortedElements(3), ielReplace)
      if (ielReplace .ne. 0) then
        call update_AllElementNeighbors2D(rhadapt, ielReplace, IsortedElements(3))
        if (ielReplace .lt. iel) then
          Imarker(IsortedElements(3)) = Imarker(ielReplace)
          Imarker(ielReplace)         = MARK_ASIS
        end if
      end if


      ! Update list of elements meeting at vertices
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i1, iel1) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i2, iel2) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i3, iel3) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        call sys_halt()
      end if

      ! Note, this can be improved by checking against JEL1 and JEL2
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i5, iel) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i5, iel2) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i5, iel3) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        call sys_halt()
      end if
      
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i1, jel1, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i1, jel2, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i2, jel1, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i3, jel2, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, jel1, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i5, jel2, ipos)


      ! Optionally, invoke callback function
      if (present(fcb_hadaptCallback) .and. present(rcollection)) then
        rcollection%IquickAccess(1:12) = (/i1, i2, i3, i4, i5, i6,&
                                           e1, e2, e3, e4, e5, e6/)
        call fcb_hadaptCallback(HADAPT_OPR_CRS_4TRIA2TRIA2, rcollection)
      end if


    case(MARK_CRS_4TRIA2TRIA_3)
      ! Update list of neighboring elements
      call update_ElementNeighbors2D(rhadapt, e1, e4, iel2, iel1, jel2, jel2)
      call update_ElementNeighbors2D(rhadapt, e2, e5, iel3, iel2, jel1, jel1)
      call update_ElementNeighbors2D(rhadapt, e3, e6, iel1, iel3, jel2, jel1)
      

      ! Update elements JEL1 and JEL2
      call replace_element2D(rhadapt, jel1, i3, i6, i2,&
                             e3, jel2, e2, e3, jel2, e5)
      call replace_element2D(rhadapt, jel2, i1, i2, i6,&
                             e1, jel1, e6, e4, jel1, e6)


      ! Delete elements IEL, IEL1, IEL2 and IEL3 depending on which
      ! elements correspond to the two elements with smallest numbers
      call remove_element2D(rhadapt, IsortedElements(4), ielReplace)
      if (ielReplace .ne. 0) then
        call update_AllElementNeighbors2D(rhadapt, ielReplace, IsortedElements(4))
        if (ielReplace .lt. iel) then
          Imarker(IsortedElements(4)) = Imarker(ielReplace)
          Imarker(ielReplace)         = MARK_ASIS
        end if
      end if

      call remove_element2D(rhadapt, IsortedElements(3), ielReplace)
      if (ielReplace .ne. 0) then
        call update_AllElementNeighbors2D(rhadapt, ielReplace, IsortedElements(3))
        if (ielReplace .lt. iel) then
          Imarker(IsortedElements(3)) = Imarker(ielReplace)
          Imarker(ielReplace)         = MARK_ASIS
        end if
      end if


      ! Update list of elements meeting at vertices.
      ! Note that all elements are removed in the first step. Afterwards,
      ! element JEL is appended to the list of elements meeting at each vertex
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i1, iel1) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i2, iel2) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i3, iel3) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        call sys_halt()
      end if

      ! Note, this can be improved by checking against JEL1 and JEL2
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i6, iel) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i6, iel1) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                     i6, iel3) .eq. ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        call sys_halt()
      end if

      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i1, jel2, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i2, jel1, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i2, jel2, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i3, jel1, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i6, jel1, ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i6, jel2, ipos)
      

      ! Optionally, invoke callback function
      if (present(fcb_hadaptCallback) .and. present(rcollection)) then
        rcollection%IquickAccess(1:12) = (/i1, i2, i3, i4, i5, i6,&
                                           e1, e2, e3, e4, e5, e6/)
        call fcb_hadaptCallback(HADAPT_OPR_CRS_4TRIA2TRIA3, rcollection)
      end if


    case DEFAULT
      call output_line('Invalid position of midpoint vertex!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
      call sys_halt()
    end select

  end subroutine coarsen_4Tria2Tria

  ! ***************************************************************************

!<subroutine>

  subroutine coarsen_4Quad1Quad(rhadapt, Imarker, iel, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine combines four quadrilaterals resulting from a
    ! 1-quad : 4-quad refinement into the original macro quadrilateral.
    ! The remaining element JEL is labeled by the smallest element number
    ! of the four neighboring elements. Due to the numbering convention
    ! all elements can be determined by visiting neighbors along the
    ! second edge starting at the given element IEL. The resulting element
    ! JEL is oriented such that local numbering starts at the vertex of the
    ! macro element, that is, such that the first vertex is the oldest one.
    !
    ! <verb>
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7)      (e3) i3          i4 (e7) i7  (e3) i3
    !      +-------+-------+             +---------------+
    !      |*      |      *|             |               |
    ! (e4) | iel3  |  iel2 | (e6)   (e4) |               | (e6)
    !      |       |       |             |               |
    !    i8+-------+-------+i6    ->     |      jel      |
    !      |     i9|       |             |               |
    ! (e8) | iel   |  iel1 | (e2)   (e8) |               | (e2)
    !      |*      |      *|             |*              |
    !      +-------+-------+             +---------------+
    !     i1 (e1) i5  (e5) i2           i1 (e1)     (e5) i2
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
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt

    ! Marker array
    integer, dimension(:), intent(inout) :: Imarker

    ! OPTIONAL: Collection
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(4) :: IsortedElements
    integer, dimension(TRIA_NVEQUAD2D) :: ImacroVertices,IvertexAge
    integer :: ipos,istate
    integer :: iel1,iel2,iel3,e1,e2,e3,e4,e5,e6,e7,e8,jel,ielReplace
    integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9


    ! Retrieve patch of elements
    iel1 = rhadapt%p_IneighboursAtElement(2, iel)
    iel2 = rhadapt%p_IneighboursAtElement(2, iel1)
    iel3 = rhadapt%p_IneighboursAtElement(2, iel2)

    ! Store vertex- and element-values of the four neighboring elements
    i1 = rhadapt%p_IverticesAtElement(1, iel)
    i5 = rhadapt%p_IverticesAtElement(2, iel)
    i9 = rhadapt%p_IverticesAtElement(3, iel)
    i8 = rhadapt%p_IverticesAtElement(4, iel)
    i2 = rhadapt%p_IverticesAtElement(1, iel1)
    i6 = rhadapt%p_IverticesAtElement(2, iel1)
    i3 = rhadapt%p_IverticesAtElement(1, iel2)
    i7 = rhadapt%p_IverticesAtElement(2, iel2)
    i4 = rhadapt%p_IverticesAtElement(1, iel3)

    ! Store values of the elements adjacent to the resulting macro element
    e1 = rhadapt%p_IneighboursAtElement(1, iel)
    e8 = rhadapt%p_IneighboursAtElement(4, iel)
    e2 = rhadapt%p_IneighboursAtElement(1, iel1)
    e5 = rhadapt%p_IneighboursAtElement(4, iel1)
    e3 = rhadapt%p_IneighboursAtElement(1, iel2)
    e6 = rhadapt%p_IneighboursAtElement(4, iel2)
    e4 = rhadapt%p_IneighboursAtElement(1, iel3)
    e7 = rhadapt%p_IneighboursAtElement(4, iel3)


    ! Sort the four elements according to their number and
    ! determine the element with the smallest element number
    IsortedElements = (/iel,iel1,iel2,iel3/)
    call sort_I32(IsortedElements, SORT_INSERT)
    jel = IsortedElements(1)


    ! Update list of neighboring elements
    call update_ElementNeighbors2D(rhadapt, e1, e5, iel1, iel, jel, jel)
    call update_ElementNeighbors2D(rhadapt, e2, e6, iel2, iel1, jel, jel)
    call update_ElementNeighbors2D(rhadapt, e3, e7, iel3, iel2, jel, jel)
    call update_ElementNeighbors2D(rhadapt, e4, e8, iel, iel3, jel, jel)


    ! The resulting quadrilateral will posses one of the states STATES_QUAD_REDx,
    ! whereby x can be 1,2,3 and 4. Die to our refinement convention, x=1,2,3
    ! should not appear, that is, local numbering of the resulting quadrilateral
    ! starts at the oldest vertex. to this end, we check the state of the
    ! provisional quadrilateral (I1,I2,I3,I4) and transform the orientation.
    ImacroVertices = (/i1,i2,i3,i4/)
    IvertexAge = rhadapt%p_IvertexAge(ImacroVertices)
    istate     = redgreen_getstateQuad(IvertexAge)

    select case(istate)
    case(STATE_QUAD_ROOT,&
         STATE_QUAD_RED4)
      ! Update element JEL = (I1,I2,I3,I4)
      call replace_element2D(rhadapt, jel, i1, i2, i3, i4,&
                             e1, e2, e3, e4, e5, e6, e7, e8)

    case(STATE_QUAD_RED1)
      ! Update element JEL = (I2,I3,I4,I1)
      call replace_element2D(rhadapt, jel, i2, i3, i4, i1,&
                             e2, e3, e4, e1, e6, e7, e8, e5)

    case(STATE_QUAD_RED2)
      ! Update element JEL = (I3,I4,I1,I2)
      call replace_element2D(rhadapt, jel, i3, i4, i1, i2,&
                             e3, e4, e1, e2, e7, e8, e5, e6)

    case(STATE_QUAD_RED3)
      ! Update element JEL = (I4,I1,I2,I3)
      call replace_element2D(rhadapt, jel, i4, i1, i2, i3,&
                             e4, e1, e2, e3, e8, e5, e6, e7)
      
    case DEFAULT
      call output_line('Invalid state of resulting quadrilateral!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad1Quad')
      call sys_halt()
    end select


    ! Delete elements IEL, IEL1, IEL2 and IEL3 depending on which
    ! element corresponds to element with minimum number JEL
    call remove_element2D(rhadapt,IsortedElements(4), ielReplace)
    if (ielReplace .ne. 0) then
      call update_AllElementNeighbors2D(rhadapt, ielReplace, IsortedElements(4))
      if (ielReplace .lt. iel) then
        Imarker(IsortedElements(4)) = Imarker(ielReplace)
        Imarker(ielReplace)         = MARK_ASIS
      end if
    end if
    
    call remove_element2D(rhadapt, IsortedElements(3), ielReplace)
    if (ielReplace .ne. 0) then
      call update_AllElementNeighbors2D(rhadapt, ielReplace, IsortedElements(3))
      if (ielReplace .lt. iel) then
        Imarker(IsortedElements(3)) = Imarker(ielReplace)
        Imarker(ielReplace)         = MARK_ASIS
      end if
    end if
    
    call remove_element2D(rhadapt, IsortedElements(2), ielReplace)
    if (ielReplace .ne. 0) then
      call update_AllElementNeighbors2D(rhadapt, ielReplace, IsortedElements(2))
      if (ielReplace .lt. iel) then
        Imarker(IsortedElements(2)) = Imarker(ielReplace)
        Imarker(ielReplace)         = MARK_ASIS
      end if
    end if


    ! Update list of elements meeting at vertices.
    ! Note that all elements are removed in the first step. Afterwards,
    ! element JEL is appended to the list of elements meeting at each vertex
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i1, iel) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad1Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i2, iel1) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad1Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i3, iel2) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad1Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,&
                                   i4, iel3) .eq. ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad1Quad')
      call sys_halt()
    end if

    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i1, jel, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i2, jel, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i3, jel, ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex, i4, jel, ipos)


    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:17) = (/i1, i2, i3, i4, i5, i6, i7, i8, i9,&
                                         e1, e2, e3, e4, e5, e6, e7, e8/)
      call fcb_hadaptCallback(HADAPT_OPR_CRS_4QUAD1QUAD, rcollection)
    end if

  end subroutine coarsen_4Quad1Quad

  ! ***************************************************************************

!<subroutine>

  subroutine coarsen_4Quad2Quad(rhadapt, Imarker, iel, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine combines four quadrilaterals resulting from a
    ! 1-quad : 4-quad refinement into two green quadrilaterals.
    ! The remaining elements JEL1 and JEL2 are labeled by the  two smallest
    ! element numbers of the four neighboring elements. Due to the numbering
    ! convention all elements can be determined by visiting neighbors along
    ! the second edge starting at the given element IEL. The resulting elements
    ! JEL1 and JEL2 are oriented such that local numbering starts at the vertices
    ! of the macro element, that is, such that the first vertex is the oldest one.
    !
    ! <verb>
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7) i7  (e3) i3          i4 (e7)  i7  (e3) i3
    !      +-------+-------+             +-------+-------+
    !      |*      |      *|             |       |      *|
    ! (e4) | iel3  |  iel2 | (e6)   (e4) |       |       | (e6)
    !      |       |       |             |       |       |
    !    i8+-------+-------+i6    ->     | jel1  | jel2  |
    !      |     i9|       |             |       |       |
    ! (e8) | iel   |  iel1 | (e2)   (e8) |       |       | (e2)
    !      |*      |      *|             |*      |       |
    !      +-------+-------+             +-------+-------+
    !     i1 (e1) i5  (e5) i2           i1 (e1) i5  (e5) i2
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
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt

    ! Marker array
    integer, dimension(:), intent(inout) :: Imarker

    ! OPTIONAL: Collection
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(4) :: IsortedElements
    integer :: ipos
    integer :: iel1,iel2,iel3,e1,e2,e3,e4,e5,e6,e7,e8,jel1,jel2,ielReplace
    integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9

    ! Retrieve patch of elements
    iel1 = rhadapt%p_IneighboursAtElement(2,iel)
    iel2 = rhadapt%p_IneighboursAtElement(2,iel1)
    iel3 = rhadapt%p_IneighboursAtElement(2,iel2)

    ! Store vertex- and element-values of the four neighboring elements
    i1 = rhadapt%p_IverticesAtElement(1,iel)
    i5 = rhadapt%p_IverticesAtElement(2,iel)
    i9 = rhadapt%p_IverticesAtElement(3,iel)
    i8 = rhadapt%p_IverticesAtElement(4,iel)
    i2 = rhadapt%p_IverticesAtElement(1,iel1)
    i6 = rhadapt%p_IverticesAtElement(2,iel1)
    i3 = rhadapt%p_IverticesAtElement(1,iel2)
    i7 = rhadapt%p_IverticesAtElement(2,iel2)
    i4 = rhadapt%p_IverticesAtElement(1,iel3)

    ! Store values of the elements adjacent to the resulting macro element
    e1 = rhadapt%p_IneighboursAtElement(1,iel)
    e8 = rhadapt%p_IneighboursAtElement(4,iel)
    e2 = rhadapt%p_IneighboursAtElement(1,iel1)
    e5 = rhadapt%p_IneighboursAtElement(4,iel1)
    e3 = rhadapt%p_IneighboursAtElement(1,iel2)
    e6 = rhadapt%p_IneighboursAtElement(4,iel2)
    e4 = rhadapt%p_IneighboursAtElement(1,iel3)
    e7 = rhadapt%p_IneighboursAtElement(4,iel3)


    ! Sort the four elements according to their number and
    ! determine the element with the smallest element number
    IsortedElements=(/iel,iel1,iel2,iel3/)
    call sort_I32(IsortedElements,SORT_INSERT)
    jel1=IsortedElements(1)
    jel2=IsortedElements(2)


    ! Update list of neighboring elements
    call update_ElementNeighbors2D(rhadapt,e1,e5,iel1,iel,jel2,jel1)
    call update_ElementNeighbors2D(rhadapt,e2,e6,iel2,iel1,jel2,jel2)
    call update_ElementNeighbors2D(rhadapt,e3,e7,iel3,iel2,jel1,jel2)
    call update_ElementNeighbors2D(rhadapt,e4,e8,iel,iel3,jel1,jel1)


    ! Update elements JEL1 = (I1,I5,I7,I4) and JEL2=(I3,I7,I5,I2)
    call replace_element2D(rhadapt,jel1,i1,i5,i7,i4,e1,jel2,e7,e4,e1,jel2,e7,e8)
    call replace_element2D(rhadapt,jel2,i3,i7,i5,i2,e3,jel1,e5,e2,e3,jel1,e5,e6)

    
    ! Delete elements IEL, IEL1, IEL2 and IEL3 depending on which
    ! elements corresponds to elements with minimum numbers JEL1 and JEL2
    call remove_element2D(rhadapt,IsortedElements(4),ielReplace)
    if (ielReplace.ne.0) then
      call update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(4))
      if (ielReplace .lt. iel) then
        Imarker(IsortedElements(4)) = Imarker(ielReplace)
        Imarker(ielReplace)         = MARK_ASIS
      end if
    end if
    
    call remove_element2D(rhadapt,IsortedElements(3),ielReplace)
    if (ielReplace.ne.0) then
      call update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(3))
      if (ielReplace .lt. iel) then
        Imarker(IsortedElements(3)) = Imarker(ielReplace)
        Imarker(ielReplace)         = MARK_ASIS
      end if
    end if
    

    ! Update list of elements meeting at vertices.
    ! Note that all elements are removed in the first step. Afterwards,
    ! element JEL is appended to the list of elements meeting at each vertex
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad2Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel1).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad2Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel2).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad2Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel3).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad2Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i5,iel).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad2Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i5,iel1).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad2Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i7,iel2).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad2Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i7,iel3).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad2Quad')
      call sys_halt()
    end if
   
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel1,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel2,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,jel2,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,jel1,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,jel1,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,jel2,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,jel1,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,jel2,ipos)


    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:17) = (/i1, i2, i3, i4, i5, i6, i7, i8, i9,&
                                         e1, e2, e3, e4, e5, e6, e7, e8/)
      call fcb_hadaptCallback(HADAPT_OPR_CRS_4QUAD2QUAD, rcollection)
    end if

  end subroutine coarsen_4Quad2Quad

  ! ***************************************************************************

!<subroutine>

  subroutine coarsen_4Quad3Tria(rhadapt, Imarker, iel, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine combines four quadrilaterals resulting from a
    ! 1-quad : 4-quad refinement into three green triangles.
    ! The elements remaining JEL1,JEL2, and JEL3 are the ones with
    ! the smallest element number, that is, only the largest element is
    ! removed. The numbering convention follows the strategy used for
    ! refinement, that is, the local numbering of the two outer triangles
    ! starts at the vertices of the macro element and the local numbering
    ! numbering of the interior triangle starts at the edge midpoint.
    !
    ! <verb>
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7) i7  (e3) i3          i4 (e7)      (e3) i3
    !      +-------+-------+             +---------------+
    !      |*      |      *|             |\             /|
    ! (e4) | iel3  |  iel2 | (e6)   (e4) | \           / | (e6)
    !      |       |       |             |  \         /  |
    !    i8+-------+-------+i6    ->     |   \  jel3 /   |
    !      |     i9|       |             |    \     /    |
    ! (e8) | iel   |  iel1 | (e2)   (e8) |jel1 \   / jel2| (e2)
    !      |*      |      *|             |*     \X/     *|
    !      +-------+-------+             +-------+-------+
    !     i1 (e1) i5  (e5) i2           i1 (e1) i5  (e5) i2
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
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt

    ! Marker array
    integer, dimension(:), intent(inout) :: Imarker

    ! OPTIONAL: Collection
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(4) :: IsortedElements
    integer :: ipos
    integer :: iel1,iel2,iel3,e1,e2,e3,e4,e5,e6,e7,e8
    integer :: jel1,jel2,jel3,ielReplace
    integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9

    ! Retrieve patch of elements
    iel1 = rhadapt%p_IneighboursAtElement(2,iel)
    iel2 = rhadapt%p_IneighboursAtElement(2,iel1)
    iel3 = rhadapt%p_IneighboursAtElement(2,iel2)

    ! Store vertex- and element-values of the four neighboring elements
    i1 = rhadapt%p_IverticesAtElement(1,iel)
    i5 = rhadapt%p_IverticesAtElement(2,iel)
    i9 = rhadapt%p_IverticesAtElement(3,iel)
    i8 = rhadapt%p_IverticesAtElement(4,iel)
    i2 = rhadapt%p_IverticesAtElement(1,iel1)
    i6 = rhadapt%p_IverticesAtElement(2,iel1)
    i3 = rhadapt%p_IverticesAtElement(1,iel2)
    i7 = rhadapt%p_IverticesAtElement(2,iel2)
    i4 = rhadapt%p_IverticesAtElement(1,iel3)

    ! Store values of the elements adjacent to the resulting macro element
    e1 = rhadapt%p_IneighboursAtElement(1,iel)
    e8 = rhadapt%p_IneighboursAtElement(4,iel)
    e2 = rhadapt%p_IneighboursAtElement(1,iel1)
    e5 = rhadapt%p_IneighboursAtElement(4,iel1)
    e3 = rhadapt%p_IneighboursAtElement(1,iel2)
    e6 = rhadapt%p_IneighboursAtElement(4,iel2)
    e4 = rhadapt%p_IneighboursAtElement(1,iel3)
    e7 = rhadapt%p_IneighboursAtElement(4,iel3)


    ! Sort the four elements according to their number and
    ! determine the elements with the smallest element numbers
    IsortedElements=(/iel,iel1,iel2,iel3/)
    call sort_I32(IsortedElements,SORT_INSERT)
    jel1=IsortedElements(1)
    jel2=IsortedElements(2)
    jel3=IsortedElements(3)


    ! Update list of neighboring elements
    call update_ElementNeighbors2D(rhadapt,e1,e5,iel1,iel,jel2,jel1)
    call update_ElementNeighbors2D(rhadapt,e2,e6,iel2,iel1,jel2,jel2)
    call update_ElementNeighbors2D(rhadapt,e3,e7,iel3,iel2,jel3,jel3)
    call update_ElementNeighbors2D(rhadapt,e4,e8,iel,iel3,jel1,jel1)


    ! Update elements JEL1 = (I1,I5,I4), JEL2 = (I2,I3,I5), and JEL3 = (I5,I3,I4)
    call replace_element2D(rhadapt,jel1,i1,i5,i4,e1,jel3,e4,e1,jel3,e8)
    call replace_element2D(rhadapt,jel2,i2,i3,i5,e2,jel3,e5,e6,jel3,e5)
    call replace_element2D(rhadapt,jel3,i5,i3,i4,jel2,e3,jel1,jel2,e7,jel1)


    ! Delete the element with the largest element number
    call remove_element2D(rhadapt,IsortedElements(4),ielReplace)
    if (ielReplace.ne.0) then
      call update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(4))
      if (ielReplace .lt. iel) then
        Imarker(IsortedElements(4)) = Imarker(ielReplace)
        Imarker(ielReplace)         = MARK_ASIS
      end if
    end if


    ! Update list of elements meeting at vertices.
    ! Note that all elements are removed in the first step. Afterwards,
    ! element JEL is appended to the list of elements meeting at each vertex
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad3Tria')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel1).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad3Tria')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel2).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad3Tria')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel3).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad3Tria')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i5,iel).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad3Tria')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i5,iel1).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad3Tria')
      call sys_halt()
    end if

    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel1,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel2,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,jel2,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,jel3,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,jel1,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,jel3,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,jel1,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,jel2,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,jel3,ipos)


    ! Adjust number of elements
    rhadapt%InelOfType(TRIA_NVETRI2D) = rhadapt%InelOfType(TRIA_NVETRI2D)+3
    rhadapt%InelOfType(TRIA_NVEQUAD2D) = rhadapt%InelOfType(TRIA_NVEQUAD2D)-3
    
    
    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:17) = (/i1, i2, i3, i4, i5, i6, i7, i8, i9,&
                                         e1, e2, e3, e4, e5, e6, e7, e8/)
      call fcb_hadaptCallback(HADAPT_OPR_CRS_4QUAD3TRIA, rcollection)
    end if

  end subroutine coarsen_4Quad3Tria

  ! ***************************************************************************

!<subroutine>

  subroutine coarsen_4Quad4Tria(rhadapt, Imarker, iel, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine combines four quadrilaterals resulting from a
    ! 1-quad : 4-quad refinement into four green triangles. At first glance
    ! this "coarsening" does not make sense since the number of elements
    ! is not reduced. However, the number of vertices is decreased by one.
    ! Moreover, the removal of this vertex may lead to further coarsening
    ! starting at the elements which meet at the deleted vertex.
    ! The numbering convention follows the strategy used for
    ! refinement, that is, the local numbering of the two outer triangles
    ! starts at the vertices of the macro element and the local numbering
    ! numbering of the interior triangle starts at the edge midpoint.
    !
    ! <verb>
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7) i7  (e3) i3          i4 (e7)      (e3) i3
    !      +-------+-------+             +---------------+
    !      |*      |      *|             |*          ---/|
    ! (e4) | iel3  |  iel2 | (e6)   (e4) |iel3   ---/ X/ | (e6)
    !      |       |       |             |   ---/     /  |
    !    i8+-------+-------+i6    ->   i8+--/  iel2  /   |
    !      |     i9|       |             |-\        /    |
    ! (e8) | iel   |  iel1 | (e2)   (e8) |  --\    /     | (e2)
    !      |*      |      *|             |*iel -\ / iel1*|
    !      +-------+-------+             +-------+-------+
    !     i1 (e1) i5  (e5) i2           i1 (e1) i5  (e5) i2
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
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt

    ! Marker array
    integer, dimension(:), intent(inout) :: Imarker

    ! OPTIONAL: Collection
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ipos
    integer :: iel1,iel2,iel3,e1,e2,e3,e4,e5,e6,e7,e8
    integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9
    
    ! Retrieve patch of elements
    iel1 = rhadapt%p_IneighboursAtElement(2,iel)
    iel2 = rhadapt%p_IneighboursAtElement(2,iel1)
    iel3 = rhadapt%p_IneighboursAtElement(2,iel2)

    ! Store vertex- and element-values of the four neighboring elements
    i1 = rhadapt%p_IverticesAtElement(1,iel)
    i5 = rhadapt%p_IverticesAtElement(2,iel)
    i9 = rhadapt%p_IverticesAtElement(3,iel)
    i8 = rhadapt%p_IverticesAtElement(4,iel)
    i2 = rhadapt%p_IverticesAtElement(1,iel1)
    i6 = rhadapt%p_IverticesAtElement(2,iel1)
    i3 = rhadapt%p_IverticesAtElement(1,iel2)
    i7 = rhadapt%p_IverticesAtElement(2,iel2)
    i4 = rhadapt%p_IverticesAtElement(1,iel3)

    ! Store values of the elements adjacent to the resulting macro element
    e1 = rhadapt%p_IneighboursAtElement(1,iel)
    e8 = rhadapt%p_IneighboursAtElement(4,iel)
    e2 = rhadapt%p_IneighboursAtElement(1,iel1)
    e5 = rhadapt%p_IneighboursAtElement(4,iel1)
    e3 = rhadapt%p_IneighboursAtElement(1,iel2)
    e6 = rhadapt%p_IneighboursAtElement(4,iel2)
    e4 = rhadapt%p_IneighboursAtElement(1,iel3)
    e7 = rhadapt%p_IneighboursAtElement(4,iel3)


    ! Update list of neighboring elements
    call update_ElementNeighbors2D(rhadapt,e2,e6,iel2,iel1,iel1,iel1)
    call update_ElementNeighbors2D(rhadapt,e3,e7,iel3,iel2,iel3,iel3)

    
    ! Update the four element
    call replace_element2D(rhadapt,iel,i1,i5,i8,e1,iel2,e8,e1,iel2,e8)
    call replace_element2D(rhadapt,iel1,i2,i3,i5,e2,iel2,e5,e6,iel2,e5)
    call replace_element2D(rhadapt,iel2,i3,i8,i5,iel3,iel,iel1,iel3,iel,iel1)
    call replace_element2D(rhadapt,iel3,i4,i8,i3,e4,iel2,e3,e4,iel2,e7)


    ! Update list of elements meeting at vertices. Note that we only have to
    ! add some elements to the vertices since all four elements are already
    ! "attached" to one of the four corner nodes
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,iel1,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,iel3,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,iel2,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i8,iel2,ipos)
    
    ! Adjust number of elements
    rhadapt%InelOfType(TRIA_NVETRI2D) = rhadapt%InelOfType(TRIA_NVETRI2D)+4
    rhadapt%InelOfType(TRIA_NVEQUAD2D) = rhadapt%InelOfType(TRIA_NVEQUAD2D)-4
    
    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:17) = (/i1, i2, i3, i4, i5, i6, i7, i8, i9,&
                                         e1, e2, e3, e4, e5, e6, e7, e8/)
      call fcb_hadaptCallback(HADAPT_OPR_CRS_4QUAD4TRIA, rcollection)
    end if

  end subroutine coarsen_4Quad4Tria
    
  ! ***************************************************************************

!<subroutine>

  subroutine coarsen_2Quad1Quad(rhadapt, Imarker, iel, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine combines two quadrilaterals resulting from a
    ! 1-quad : 2-quad refinement into the macro quadrilateral.
    ! The remaining element JEL is labeled by the smallest element number
    ! of the four neighboring elements. Due to the numbering convention
    ! all elements can be determined by visiting neighbors along the
    ! second edge starting at the given element IEL. The resulting element
    ! JEL is oriented such that local numbering starts at the vertex of the
    ! macro element, that is, such that the first vertex is the oldest one.
    !
    ! <verb>
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7) i7  (e3) i3          i4 (e7)      (e3) i3
    !      +-------+-------+             +---------------+
    !      |       |      *|             |               |
    ! (e4) |       |       | (e6)   (e4) |               | (e6)
    !      |       |       |             |               |
    !      |  iel  |  iel1 |      ->     |      jel      |
    !      |       |       |             |               |
    ! (e8) |       |       | (e2)   (e8) |               | (e2)
    !      |*      |      *|             |*              |
    !      +-------+-------+             +---------------+
    !     i1 (e1) i5  (e5) i2           i1 (e1)     (e5) i2
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
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt

    ! Marker array
    integer, dimension(:), intent(inout) :: Imarker

    ! OPTIONAL: Collection
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(TRIA_NVEQUAD2D) :: ImacroVertices, IvertexAge
    integer :: ipos,istate
    integer :: iel1,e1,e2,e3,e4,e5,e6,e7,e8,ielReplace
    integer :: i1,i2,i3,i4,i5,i7
    
    ! Retrieve neighboring element
    iel1 = rhadapt%p_IneighboursAtElement(2,iel)

    ! Store vertex- and element-values of the four neighboring elements
    i1 = rhadapt%p_IverticesAtElement(1,iel)
    i5 = rhadapt%p_IverticesAtElement(2,iel)
    i7 = rhadapt%p_IverticesAtElement(3,iel)
    i4 = rhadapt%p_IverticesAtElement(4,iel)
    i3 = rhadapt%p_IverticesAtElement(1,iel1)
    i2 = rhadapt%p_IverticesAtElement(4,iel1)

    ! Store values of the elements adjacent to the resulting macro element
    e1 = rhadapt%p_IneighboursAtElement(1,iel)
    e7 = rhadapt%p_IneighboursAtElement(3,iel)
    e4 = rhadapt%p_IneighboursAtElement(4,iel)
    e3 = rhadapt%p_IneighboursAtElement(1,iel1)
    e5 = rhadapt%p_IneighboursAtElement(3,iel1)
    e2 = rhadapt%p_IneighboursAtElement(4,iel1)

    e8 = rhadapt%p_ImidneighboursAtElement(4,iel)
    e6 = rhadapt%p_ImidneighboursAtElement(4,iel1)


    ! Which is the smaller element?
    if (iel .lt. iel1) then
      
      ! Update list of neighboring elements
      call update_ElementNeighbors2D(rhadapt,e1,e5,iel1,iel,iel,iel)
      call update_ElementNeighbors2D(rhadapt,e2,e6,iel1,iel,iel)
      call update_ElementNeighbors2D(rhadapt,e3,e7,iel,iel1,iel,iel)


      ! The resulting quadrilateral will posses one of the states STATES_QUAD_REDx,
      ! whereby x can be 1,2,3 and 4. Die to our refinement convention, x=1,2,3
      ! should not appear, that is, local numbering of the resulting quadrilateral
      ! starts at the oldest vertex. to this end, we check the state of the
      ! provisional quadrilateral (I1,I2,I3,I4) and transform the orientation.
      ImacroVertices = (/i1,i2,i3,i4/)
      IvertexAge = rhadapt%p_IvertexAge(ImacroVertices)
      istate = redgreen_getstateQuad(IvertexAge)
      
      select case(istate)
      case(STATE_QUAD_ROOT,STATE_QUAD_RED4)
        ! Update element IEL = (I1,I2,I3,I4)
        call replace_element2D(rhadapt,iel,i1,i2,i3,i4,e1,e2,e3,e4,e5,e6,e7,e8)
        
      case(STATE_QUAD_RED1)
        ! Update element IEL = (I2,I3,I4,I1)
        call replace_element2D(rhadapt,iel,i2,i3,i4,i1,e2,e3,e4,e1,e6,e7,e8,e5)

      case(STATE_QUAD_RED2)
        ! Update element IEL = (I3,I4,I1,I2)
        call replace_element2D(rhadapt,iel,i3,i4,i1,i2,e3,e4,e1,e2,e7,e8,e5,e6)

      case(STATE_QUAD_RED3)
        ! Update element IEL = (I4,I1,I2,I3)
        call replace_element2D(rhadapt,iel,i4,i1,i2,i3,e4,e1,e2,e3,e8,e5,e6,e7)
        
      case DEFAULT
        call output_line('Invalid state of resulting quadrilateral!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Quad1Quad')
        call sys_halt()
      end select


      ! Delete element IEL1
      call remove_element2D(rhadapt,iel1,ielReplace)
      if (ielReplace.ne.0) then
        call update_AllElementNeighbors2D(rhadapt,ielReplace,iel1)
        if (ielReplace .lt. iel) then
          Imarker(iel1)       = Imarker(ielReplace)
          Imarker(ielReplace) = MARK_ASIS
        end if
      end if


      ! Update list of elements meeting at vertices.
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel1).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Quad1Quad')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel1).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Quad1Quad')
        call sys_halt()
      end if

      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,iel,ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,iel,ipos)

    else

      ! Update list of neighboring elements
      call update_ElementNeighbors2D(rhadapt,e1,e5,iel1,iel,iel1,iel1)
      call update_ElementNeighbors2D(rhadapt,e3,e7,iel,iel1,iel1,iel1)
      call update_ElementNeighbors2D(rhadapt,e4,e8,iel,iel1,iel1)
      
      ! The resulting quadrilateral will posses one of the states STATES_QUAD_REDx,
      ! whereby x can be 1,2,3 and 4. Die to our refinement convention, x=1,2,3
      ! should not appear, that is, local numbering of the resulting quadrilateral
      ! starts at the oldest vertex. to this end, we check the state of the
      ! provisional quadrilateral (I1,I2,I3,I4) and transform the orientation.
      ImacroVertices = (/i1,i2,i3,i4/)
      IvertexAge = rhadapt%p_IvertexAge(ImacroVertices)
      istate = redgreen_getstateQuad(IvertexAge)
      
      select case(istate)
      case(STATE_QUAD_ROOT,STATE_QUAD_RED4)
        ! Update element IEL1 = (I1,I2,I3,I4)
        call replace_element2D(rhadapt,iel1,i1,i2,i3,i4,e1,e2,e3,e4,e5,e6,e7,e8)
        
      case(STATE_QUAD_RED1)
        ! Update element IEL1 = (I2,I3,I4,I1)
        call replace_element2D(rhadapt,iel1,i2,i3,i4,i1,e2,e3,e4,e1,e6,e7,e8,e5)

      case(STATE_QUAD_RED2)
        ! Update element IEL1 = (I3,I4,I1,I2)
        call replace_element2D(rhadapt,iel1,i3,i4,i1,i2,e3,e4,e1,e2,e7,e8,e5,e6)

      case(STATE_QUAD_RED3)
        ! Update element IEL1 = (I4,I1,I2,I3)
        call replace_element2D(rhadapt,iel1,i4,i1,i2,i3,e4,e1,e2,e3,e8,e5,e6,e7)
        
      case DEFAULT
        call output_line('Invalid state of resulting quadrilateral!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Quad1Quad')
        call sys_halt()
      end select


      ! Delete element IEL
      call remove_element2D(rhadapt,iel,ielReplace)
      if (ielReplace.ne.0) then
        call update_AllElementNeighbors2D(rhadapt,ielReplace,iel)
        if (ielReplace .lt. iel) then
          Imarker(iel)        = Imarker(ielReplace)
          Imarker(ielReplace) = MARK_ASIS
        end if
      end if
      

      ! Update list of elements meeting at vertices.
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Quad1Quad')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Quad1Quad')
        call sys_halt()
      end if

      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,iel1,ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,iel1,ipos)

    end if


    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:16) = (/i1, i2, i3, i4, i5, 0, i7, 0,&
                                         e1, e2, e3, e4, e5, e6, e7, e8/)
      call fcb_hadaptCallback(HADAPT_OPR_CRS_2QUAD1QUAD, rcollection)
    end if

  end subroutine coarsen_2Quad1Quad

  ! ***************************************************************************

!<subroutine>

  subroutine coarsen_2Quad3Tria(rhadapt, Imarker, iel, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine combines two quadrilaterals resulting from a
    ! 1-quad : 2-quad refinement into three green triangles. At first glance
    ! this "coarsening" does not make sense since the number of elements
    ! is even increased. However, the number of vertices is decreased by one.
    ! Moreover, the removal of this vertex may lead to further coarsening
    ! starting at the elements which meet at the deleted vertex.
    ! The numbering convention follows the strategy used for
    ! refinement, that is, the local numbering of the two outer triangles
    ! starts at the vertices of the macro element and the local numbering
    ! numbering of the interior triangle starts at the edge midpoint.
    !
    ! <verb>
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7) i7  (e3) i3          i4 (e7)      (e3) i3
    !      +-------+-------+             +---------------+
    !      |       |      *|             |\             /|
    ! (e4) |       |       | (e6)   (e4) | \           / | (e6)
    !      |       |       |             |  \  nel+1  /  |
    !      |  iel  |  iel1 |      ->     |   \       /   |
    !      |       |       |             |    \     /    |
    ! (e8) |       |       | (e2)   (e8) | iel \   / iel1| (e2)
    !      |*      |      *|             |*     \X/     *|
    !      +-------+-------+             +-------+-------+
    !     i1 (e1) i5  (e5) i2           i1 (e1) i5  (e5) i2
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
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt

    ! Marker array
    integer, dimension(:), intent(inout) :: Imarker

    ! OPTIONAL: Collection
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ipos
    integer :: iel1,nel0,e1,e2,e3,e4,e5,e6,e7,e8
    integer :: i1,i2,i3,i4,i5,i7

    ! Retrieve neighboring element
    iel1 = rhadapt%p_IneighboursAtElement(2,iel)

    ! Store vertex- and element-values of the four neighboring elements
    i1 = rhadapt%p_IverticesAtElement(1,iel)
    i5 = rhadapt%p_IverticesAtElement(2,iel)
    i7 = rhadapt%p_IverticesAtElement(3,iel)
    i4 = rhadapt%p_IverticesAtElement(4,iel)
    i3 = rhadapt%p_IverticesAtElement(1,iel1)
    i2 = rhadapt%p_IverticesAtElement(4,iel1)

    ! Store values of the elements adjacent to the resulting macro element
    e1 = rhadapt%p_IneighboursAtElement(1,iel)
    e7 = rhadapt%p_IneighboursAtElement(3,iel)
    e4 = rhadapt%p_IneighboursAtElement(4,iel)
    e3 = rhadapt%p_IneighboursAtElement(1,iel1)
    e5 = rhadapt%p_IneighboursAtElement(3,iel1)
    e2 = rhadapt%p_IneighboursAtElement(4,iel1)

    e8 = rhadapt%p_ImidneighboursAtElement(4,iel)
    e6 = rhadapt%p_ImidneighboursAtElement(4,iel1)

    ! Store total number of elements before refinement
    nel0 = rhadapt%NEL
    

    ! Replace elements IEL and IEL1 and add one new element JEL
    call replace_element2D(rhadapt,iel,i1,i5,i4,e1,nel0+1,e4,e1,nel0+1,e8)
    call replace_element2D(rhadapt,iel1,i2,i3,i5,e2,nel0+1,e5,e6,nel0+1,e5)
    call add_element2D(rhadapt,i5,i3,i4,iel1,e3,iel,iel1,e7,iel)


    ! Update list of neighboring elements
    call update_ElementNeighbors2D(rhadapt,e3,e7,iel,iel1,nel0+1,nel0+1)


    ! Update list of elements meeting at vertices.
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,nel0+1,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+1,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+1,ipos)


    ! Adjust number of elements
    rhadapt%InelOfType(TRIA_NVETRI2D) = rhadapt%InelOfType(TRIA_NVETRI2D)+2
    rhadapt%InelOfType(TRIA_NVEQUAD2D) = rhadapt%InelOfType(TRIA_NVEQUAD2D)-2

    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:16) = (/i1, i2, i3, i4, i5, 0, i7, 0,&
                                         e1, e2, e3, e4, e5, e6, e7, e8/)
      call fcb_hadaptCallback(HADAPT_OPR_CRS_2QUAD3TRIA, rcollection)
    end if

  end subroutine coarsen_2Quad3Tria

  ! ***************************************************************************

!<subroutine>

  subroutine coarsen_3Tria1Quad(rhadapt, Imarker, iel, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine combines three triangles resulting from a
    ! 1-quad : 3-tria refinement into the macro quadrilateral.
    !
    ! <verb>
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7)     (e3) i3          i4 (e7)      (e3) i3
    !      +---------------+             +---------------+
    !      |\             /|             |               |
    ! (e4) | \           / | (e6)   (e4) |               | (e6)
    !      |  \   iel   /  |             |               |
    !      |   \       /   |      ->     |      jel      |
    !      |    \     /    |             |               |
    ! (e8) | iel1\   /iel2 | (e2)   (e8) |               | (e2)
    !      |*     \X/     *|             |*              |
    !      +-------+-------+             +---------------+
    !     i1 (e1) i5  (e5) i2           i1 (e1)     (e5) i2
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
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt

    ! Marker array
    integer, dimension(:), intent(inout) :: Imarker

    ! OPTIONAL: Collection
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer, dimension(3) :: IsortedElements
    integer, dimension(TRIA_NVEQUAD2D) :: ImacroVertices,IvertexAge
    integer :: ipos,istate
    integer :: iel1,iel2,e1,e2,e3,e4,e5,e6,e7,e8,jel,ielReplace
    integer :: i1,i2,i3,i4,i5

    ! Retrieve neighboring element
    iel2 = rhadapt%p_IneighboursAtElement(1,iel)
    iel1 = rhadapt%p_IneighboursAtElement(3,iel)

    ! Store vertex- and element-values of the four neighboring elements
    i5 = rhadapt%p_IverticesAtElement(1,iel)
    i3 = rhadapt%p_IverticesAtElement(2,iel)
    i4 = rhadapt%p_IverticesAtElement(3,iel)
    i1 = rhadapt%p_IverticesAtElement(1,iel1)
    i2 = rhadapt%p_IverticesAtElement(1,iel2)

    ! Store values of the elements adjacent to the resulting macro element
    e3 = rhadapt%p_IneighboursAtElement(2,iel)
    e1 = rhadapt%p_IneighboursAtElement(1,iel1)
    e4 = rhadapt%p_IneighboursAtElement(3,iel1)
    e2 = rhadapt%p_IneighboursAtElement(1,iel2)
    e5 = rhadapt%p_IneighboursAtElement(3,iel2)

    e7 = rhadapt%p_ImidneighboursAtElement(2,iel)
    e8 = rhadapt%p_ImidneighboursAtElement(3,iel1)
    e6 = rhadapt%p_ImidneighboursAtElement(1,iel2)


    ! Sort the three elements according to their number and
    ! determine the element with the smallest element number
    IsortedElements=(/iel,iel1,iel2/)
    call sort_I32(IsortedElements,SORT_INSERT)
    jel=IsortedElements(1)

    ! Update list of neighboring elements
    call update_ElementNeighbors2D(rhadapt,e1,e5,iel2,iel1,jel,jel)
    call update_ElementNeighbors2D(rhadapt,e2,e6,iel2,jel,jel)
    call update_ElementNeighbors2D(rhadapt,e3,e7,iel,jel,jel)
    call update_ElementNeighbors2D(rhadapt,e4,e8,iel1,jel,jel)

    ! The resulting quadrilateral will posses one of the states STATES_QUAD_REDx,
    ! whereby x can be 1,2,3 and 4. Die to our refinement convention, x=1,2,3
    ! should not appear, that is, local numbering of the resulting quadrilateral
    ! starts at the oldest vertex. to this end, we check the state of the
    ! provisional quadrilateral (I1,I2,I3,I4) and transform the orientation.
    ImacroVertices = (/i1,i2,i3,i4/)
    IvertexAge = rhadapt%p_IvertexAge(ImacroVertices)
    istate = redgreen_getstateQuad(IvertexAge)

    select case(istate)
    case(STATE_QUAD_ROOT,STATE_QUAD_RED4)
      ! Update element JEL = (I1,I2,I3,I4)
      call replace_element2D(rhadapt,jel,i1,i2,i3,i4,e1,e2,e3,e4,e5,e6,e7,e8)
      
    case(STATE_QUAD_RED1)
      ! Update element JEL = (I2,I3,I4,I1)
      call replace_element2D(rhadapt,jel,i2,i3,i4,i1,e2,e3,e4,e1,e6,e7,e8,e5)
      
    case(STATE_QUAD_RED2)
      ! Update element JEL = (I3,I4,I1,I2)
      call replace_element2D(rhadapt,jel,i3,i4,i1,i2,e3,e4,e1,e2,e7,e8,e5,e6)
      
    case(STATE_QUAD_RED3)
      ! Update element JEL = (I4,I1,I2,I3)
      call replace_element2D(rhadapt,jel,i4,i1,i2,i3,e4,e1,e2,e3,e8,e5,e6,e7)
      
    case DEFAULT
      call output_line('Invalid state of resulting quadrilateral!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_3Tria1Quad')
      call sys_halt()
    end select


    ! Delete the two elements with the largest element numbers
    call remove_element2D(rhadapt,IsortedElements(3),ielReplace)
    if (ielReplace.ne.0) then
      call update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(3))
      if (ielReplace .lt. iel) then
        Imarker(IsortedElements(3)) = Imarker(ielReplace)
        Imarker(ielReplace)         = MARK_ASIS
      end if
    end if

    call remove_element2D(rhadapt,IsortedElements(2),ielReplace)
    if (ielReplace.ne.0) then
      call update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(2))
      if (ielReplace .lt. iel) then
        Imarker(IsortedElements(2)) = Imarker(ielReplace)
        Imarker(ielReplace)         = MARK_ASIS
      end if
    end if


    ! Update list of elements meeting at vertices.
    ! Note that all elements are removed in the first step. Afterwards,
    ! element JEL is appended to the list of elements meeting at each vertex
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel1).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_3Tria1Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel2).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_3Tria1Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_3Tria1Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel2).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_3Tria1Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_3Tria1Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel1).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_3Tria1Quad')
      call sys_halt()
    end if
    
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,jel,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,jel,ipos)


    ! Adjust number of elements
    rhadapt%InelOfType(TRIA_NVETRI2D) = rhadapt%InelOfType(TRIA_NVETRI2D)-1
    rhadapt%InelOfType(TRIA_NVEQUAD2D) = rhadapt%InelOfType(TRIA_NVEQUAD2D)+1


    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:13) = (/i1, i2, i3, i4, i5,&
                                         e1, e2, e3, e4, e5, e6, e7, e8/)
      call fcb_hadaptCallback(HADAPT_OPR_CRS_3TRIA1QUAD, rcollection)
    end if

  end subroutine coarsen_3Tria1Quad

  ! ***************************************************************************

!<subroutine>

  subroutine coarsen_4Tria1Quad(rhadapt, Imarker, iel, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine combines four triangles resulting from a
    ! 1-quad : 4-tria refinement into the macro quadrilateral.
    ! The remaining element JEL is the ones with the smallest element number.
    ! The numbering convention follows the strategy used for refinement.
    ! <verb>
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7)  i7 (e3) i3          i4 (e7)      (e3) i3
    !      +-------+-------+             +---------------+
    !      |*     / \-iel2*|             |               |
    ! (e4) |iel3 /    \--  | (e6)   (e4) |               | (e6)
    !      |    /        \-|             |               |
    !      |   /  iel    --+i6    ->     |      jel      |
    !      |  /      ---/  |             |               |
    ! (e8) | /X  ---/      | (e2)   (e8) |               | (e2)
    !      |/---/    iel1 *|             |*              |
    !      +---------------+             +---------------+
    !     i1 (e1)     (e5) i2           i1 (e1)     (e5) i2
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
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
    
    ! Marker array
    integer, dimension(:), intent(inout) :: Imarker

    ! OPTIONAL: Collection
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(4) :: IsortedElements
    integer, dimension(TRIA_NVEQUAD2D) :: ImacroVertices,IvertexAge
    integer :: ipos,istate
    integer :: iel1,iel2,iel3,e1,e2,e3,e4,e5,e6,e7,e8,jel,ielReplace
    integer :: i1,i2,i3,i4,i6,i7

    ! Retrieve patch of elements
    iel1 = rhadapt%p_IneighboursAtElement(1,iel)
    iel2 = rhadapt%p_IneighboursAtElement(2,iel)
    iel3 = rhadapt%p_IneighboursAtElement(3,iel)

    ! Store vertex- and element-values of the four neighboring elements
    i1 = rhadapt%p_IverticesAtElement(1,iel)
    i6 = rhadapt%p_IverticesAtElement(2,iel)
    i7 = rhadapt%p_IverticesAtElement(3,iel)
    i2 = rhadapt%p_IverticesAtElement(1,iel1)
    i3 = rhadapt%p_IverticesAtElement(1,iel2)
    i4 = rhadapt%p_IverticesAtElement(1,iel3)

    ! Store values of the elements adjacent to the resulting macro element
    e2 = rhadapt%p_IneighboursAtElement(1,iel1)
    e1 = rhadapt%p_IneighboursAtElement(3,iel1)
    e3 = rhadapt%p_IneighboursAtElement(1,iel2)
    e6 = rhadapt%p_IneighboursAtElement(3,iel2)
    e4 = rhadapt%p_IneighboursAtElement(1,iel3)
    e7 = rhadapt%p_IneighboursAtElement(3,iel3)

    e5 = rhadapt%p_ImidneighboursAtElement(3,iel1)
    e8 = rhadapt%p_ImidneighboursAtElement(1,iel3)

    
    ! Sort the four elements according to their number and
    ! determine the elements with the smallest element numbers
    IsortedElements=(/iel,iel1,iel2,iel3/)
    call sort_I32(IsortedElements,SORT_INSERT)
    jel=IsortedElements(1)

    
    ! Update list of neighboring elements
    call update_ElementNeighbors2D(rhadapt,e1,e5,iel1,jel,jel)
    call update_ElementNeighbors2D(rhadapt,e2,e6,iel2,iel1,jel,jel)
    call update_ElementNeighbors2D(rhadapt,e3,e7,iel3,iel2,jel,jel)
    call update_ElementNeighbors2D(rhadapt,e4,e8,iel3,jel,jel)


    ! The resulting quadrilateral will posses one of the states STATES_QUAD_REDx,
    ! whereby x can be 1,2,3 and 4. Die to our refinement convention, x=1,2,3
    ! should not appear, that is, local numbering of the resulting quadrilateral
    ! starts at the oldest vertex. to this end, we check the state of the
    ! provisional quadrilateral (I1,I2,I3,I4) and transform the orientation.
    ImacroVertices = (/i1,i2,i3,i4/)
    IvertexAge = rhadapt%p_IvertexAge(ImacroVertices)
    istate = redgreen_getstateQuad(IvertexAge)

    select case(istate)
    case(STATE_QUAD_ROOT,STATE_QUAD_RED4)
      ! Update element JEL = (I1,I2,I3,I4)
      call replace_element2D(rhadapt,jel,i1,i2,i3,i4,e1,e2,e3,e4,e5,e6,e7,e8)

    case(STATE_QUAD_RED1)
      ! Update element JEL = (I2,I3,I4,I1)
      call replace_element2D(rhadapt,jel,i2,i3,i4,i1,e2,e3,e4,e1,e6,e7,e8,e5)

    case(STATE_QUAD_RED2)
      ! Update element JEL = (I3,I4,I1,I2)
      call replace_element2D(rhadapt,jel,i3,i4,i1,i2,e3,e4,e1,e2,e7,e8,e5,e6)

    case(STATE_QUAD_RED3)
      ! Update element JEL = (I4,I1,I2,I3)
      call replace_element2D(rhadapt,jel,i4,i1,i2,i3,e4,e1,e2,e3,e8,e5,e6,e7)
      
    case DEFAULT
      call output_line('Invalid state of resulting quadrilateral!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Quad')
      call sys_halt()
    end select

    
    ! Delete elements IEL, IEL1, IEL2 and IEL3 depending on which
    ! element corresponds to element with minimum number JEL
    call remove_element2D(rhadapt,IsortedElements(4),ielReplace)
    if (ielReplace.ne.0) then
      call update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(4))
      if (ielReplace .lt. iel) then
        Imarker(IsortedElements(4)) = Imarker(ielReplace)
        Imarker(ielReplace)         = MARK_ASIS
      end if
    end if

    call remove_element2D(rhadapt,IsortedElements(3),ielReplace)
    if (ielReplace.ne.0) then
      call update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(3))
      if (ielReplace .lt. iel) then
        Imarker(IsortedElements(3)) = Imarker(ielReplace)
        Imarker(ielReplace)         = MARK_ASIS
      end if
    end if
    
    call remove_element2D(rhadapt,IsortedElements(2),ielReplace)
    if (ielReplace.ne.0) then
      call update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(2))
      if (ielReplace .lt. iel) then
        Imarker(IsortedElements(2)) = Imarker(ielReplace)
        Imarker(ielReplace)         = MARK_ASIS
      end if
    end if

    
    ! Update list of elements meeting at vertices.
    ! Note that all elements are removed in the first step. Afterwards,
    ! element JEL is appended to the list of elements meeting at each vertex
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel1).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel3).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel1).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel2).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Quad')
      call sys_halt()
    end if
    if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel3).eq.&
        ARRAYLIST_NOT_FOUND) then
      call output_line('Unable to delete element from vertex list!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Quad')
      call sys_halt()
    end if
        
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,jel,ipos)
    call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,jel,ipos)
    
    
    ! Adjust number of elements
    rhadapt%InelOfType(TRIA_NVETRI2D) = rhadapt%InelOfType(TRIA_NVETRI2D)-1
    rhadapt%InelOfType(TRIA_NVEQUAD2D) = rhadapt%InelOfType(TRIA_NVEQUAD2D)+1
    

    ! Optionally, invoke callback function
    if (present(fcb_hadaptCallback) .and. present(rcollection)) then
      rcollection%IquickAccess(1:16) = (/i1, i2, i3, i4, 0, i6, i7, 0,&
                                         e1, e2, e3, e4, e5, e6, e7, e8/)
      call fcb_hadaptCallback(HADAPT_OPR_CRS_4TRIA1QUAD, rcollection)
    end if

  end subroutine coarsen_4Tria1Quad

  ! ***************************************************************************

!<subroutine>

  subroutine coarsen_4Tria3Tria(rhadapt, Imarker, iel, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine combines four triangles resulting from a
    ! 1-quad : 4-tria refinement into three green triangles. The remaining
    ! elements JEL1, JEL2 and JEL3 are the ones with the smallest element number.
    ! The numbering convention follows the strategy used for refinement.
    ! <verb>
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7)  i7 (e3) i3          i4 (e7)      (e3) i3
    !      +-------+-------+             +---------------+
    !      |*     / \-iel2*|             |---\          *|
    ! (e4) |iel3 /    \--  | (e6)   (e4) |    ---\ jel2  | (e6)
    !      |    /        \-|             |        ---\   |
    !      |   /  iel    --+i6    ->     |  jel3     *-- +i6
    !      |  /      ---/  |             |        ---/   |
    ! (e8) | /X  ---/      | (e2)   (e8) |    ---/       | (e2)
    !      |/---/    iel1 *|             |---/     jel1 *|
    !      +---------------+             +---------------+
    !     i1 (e1)     (e5) i2           i1 (e1)     (e5) i2
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
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt

    ! Marker array
    integer, dimension(:), intent(inout) :: Imarker

    ! OPTIONAL: Collection
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(4) :: IsortedElements
    integer :: ipos
    integer :: iel1,iel2,iel3,e1,e2,e3,e4,e5,e6,e7,e8
    integer :: jel1,jel2,jel3,ielReplace
    integer :: i1,i2,i3,i4,i6,i7

    ! Retrieve patch of elements
    iel1 = rhadapt%p_IneighboursAtElement(1,iel)
    iel2 = rhadapt%p_IneighboursAtElement(2,iel)
    iel3 = rhadapt%p_IneighboursAtElement(3,iel)

    ! Store vertex- and element-values of the four neighboring elements
    i1 = rhadapt%p_IverticesAtElement(1,iel)
    i6 = rhadapt%p_IverticesAtElement(2,iel)
    i7 = rhadapt%p_IverticesAtElement(3,iel)
    i2 = rhadapt%p_IverticesAtElement(1,iel1)
    i3 = rhadapt%p_IverticesAtElement(1,iel2)
    i4 = rhadapt%p_IverticesAtElement(1,iel3)

    ! Store values of the elements adjacent to the resulting macro element
    e2 = rhadapt%p_IneighboursAtElement(1,iel1)
    e1 = rhadapt%p_IneighboursAtElement(3,iel1)
    e3 = rhadapt%p_IneighboursAtElement(1,iel2)
    e6 = rhadapt%p_IneighboursAtElement(3,iel2)
    e4 = rhadapt%p_IneighboursAtElement(1,iel3)
    e7 = rhadapt%p_IneighboursAtElement(3,iel3)

    e5 = rhadapt%p_ImidneighboursAtElement(3,iel1)
    e8 = rhadapt%p_ImidneighboursAtElement(1,iel3)


    ! Sort the four elements according to their number and
    ! determine the elements with the smallest element numbers
    IsortedElements=(/iel,iel1,iel2,iel3/)
    call sort_I32(IsortedElements,SORT_INSERT)
    jel1=IsortedElements(1)
    jel2=IsortedElements(2)
    jel3=IsortedElements(3)
    
    ! Which midpoint vertex should be kept?
    select case(Imarker(iel))
    case(MARK_CRS_4TRIA3TRIA_RIGHT)
      ! Update list of neighboring elements
      call update_ElementNeighbors2D(rhadapt,e1,e5,iel1,jel1,jel1)
      call update_ElementNeighbors2D(rhadapt,e2,e6,iel2,iel1,jel2,jel1)
      call update_ElementNeighbors2D(rhadapt,e3,e7,iel3,iel2,jel2,jel2)
      call update_ElementNeighbors2D(rhadapt,e4,e8,iel3,jel3,jel3)


      ! Update elements JEL1, JEL2 and JEL3
      call replace_element2D(rhadapt,jel1,i2,i6,i1,e2,jel3,e1,e2,jel3,e5)
      call replace_element2D(rhadapt,jel2,i3,i4,i6,e3,jel3,e6,e7,jel3,e6)
      call replace_element2D(rhadapt,jel3,i6,i4,i1,jel2,e4,jel1,jel2,e8,jel1)

      
      ! Delete elements IEL, IEL1, IEL2 and IEL3 depending on which
      ! elements correspond to the two elements with smallest numbers
      call remove_element2D(rhadapt,IsortedElements(4),ielReplace)
      if (ielReplace.ne.0) then
        call update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(4))
        if (ielReplace .lt. iel) then
          Imarker(IsortedElements(4)) = Imarker(ielReplace)
          Imarker(ielReplace)         = MARK_ASIS
        end if
      end if


      ! Update list of elements meeting at vertices.
      ! Note that all elements are removed in the first step. Afterwards,
      ! element JEL is appended to the list of elements meeting at each vertex
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel1).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel3).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel1).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel2).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel3).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i6,iel).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i6,iel1).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i6,iel2).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        call sys_halt()
      end if

      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel1,ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel3,ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel1,ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,jel2,ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,jel2,ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,jel3,ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,jel1,ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,jel2,ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,jel3,ipos)

      
      ! Optionally, invoke callback function
      if (present(fcb_hadaptCallback) .and. present(rcollection)) then
        rcollection%IquickAccess(1:16) = (/i1, i2, i3, i4, 0, i6, i7, 0,&
                                           e1, e2, e3, e4, e5, e6, e7, e8/)
        call fcb_hadaptCallback(HADAPT_OPR_CRS_4TRIA3TRIA2, rcollection)
      end if
      

    case(MARK_CRS_4TRIA3TRIA_LEFT)
      ! Update list of neighboring elements
      call update_ElementNeighbors2D(rhadapt,e1,e5,iel1,jel3,jel3)
      call update_ElementNeighbors2D(rhadapt,e2,e6,iel2,iel1,jel1,jel1)
      call update_ElementNeighbors2D(rhadapt,e3,e7,iel3,iel2,jel2,jel1)
      call update_ElementNeighbors2D(rhadapt,e4,e8,iel3,jel2,jel2)

      
      ! Update elements JEL1, JEL2 and JEL3
      call replace_element2D(rhadapt,jel1,i3,i7,i2,e3,jel3,e2,e3,jel3,e6)
      call replace_element2D(rhadapt,jel2,i4,i1,i7,e4,jel3,e7,e8,jel3,e7)
      call replace_element2D(rhadapt,jel3,i7,i1,i2,jel2,e1,jel1,jel2,e5,jel1)

      
      ! Delete elements IEL, IEL1, IEL2 and IEL3 depending on which
      ! elements correspond to the two elements with smallest numbers
      call remove_element2D(rhadapt,IsortedElements(4),ielReplace)
      if (ielReplace.ne.0) then
        call update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(4))
        if (ielReplace .lt. iel) then
          Imarker(IsortedElements(4)) = Imarker(ielReplace)
          Imarker(ielReplace)         = MARK_ASIS
        end if
      end if


      ! Update list of elements meeting at vertices.
      ! Note that all elements are removed in the first step. Afterwards,
      ! element JEL is appended to the list of elements meeting at each vertex
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel1).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel3).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel1).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel2).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel3).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i7,iel).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i7,iel2).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        call sys_halt()
      end if
      if (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i7,iel3).eq.&
          ARRAYLIST_NOT_FOUND) then
        call output_line('Unable to delete element from vertex list!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        call sys_halt()
      end if

      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel2,ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel3,ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel1,ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel3,ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,jel1,ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,jel2,ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,jel1,ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,jel2,ipos)
      call arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,jel3,ipos)

      ! Optionally, invoke callback function
      if (present(fcb_hadaptCallback) .and. present(rcollection)) then
        rcollection%IquickAccess(1:16) = (/i1, i2, i3, i4, 0, i6, i7, 0,&
                                           e1, e2, e3, e4, e5, e6, e7, e8/)
        call fcb_hadaptCallback(HADAPT_OPR_CRS_4TRIA3TRIA3, rcollection)
      end if

      
    case DEFAULT
      call output_line('Invalid position of midpoint vertex!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
      call sys_halt()
    end select

  end subroutine coarsen_4Tria3Tria

  ! ***************************************************************************

!<function>

  pure function redgreen_getState2D(rhadapt, iel) result(istate)

!<description>
    ! This function encodes the state of an element iel in 2D.
    ! Note: The state of an element is solely determined from the
    ! age information of its surrounding vertices.
!</description>

!<input>
    ! Adaptivity structure
    type(t_hadapt), intent(in) :: rhadapt

    ! Number of element for which state should be computed
    integer, intent(in) :: iel
!</input>

!<result>
    ! State of element
    integer :: istate
!</result>
!</function>

    integer, dimension(TRIA_MAXNVE) :: IvertexAge

    ! Do we have quadrilaterals in the triangulation?
    if (rhadapt%InelOfType(TRIA_NVEQUAD2D) .eq. 0) then
      
      ! There are no quadrilaterals in the current triangulation.
      IvertexAge(1:TRIA_NVETRI2D) = rhadapt%p_IvertexAge(&
          rhadapt%p_IverticesAtElement(1:TRIA_NVETRI2D,iel))
      istate = redgreen_getStateTria(IvertexAge(1:TRIA_NVETRI2D))
      
    else
      
      ! There are quadrilaterals and possible also triangles in
      ! the current triangulation. If the last entry of the vertices
      ! at element list is nonzero then the current element is a
      ! quadrilateral. Otherwise, we are dealing with  a triangle
      if (rhadapt%p_IverticesAtElement(TRIA_NVEQUAD2D,iel) .eq. 0) then
        IvertexAge(1:TRIA_NVETRI2D) = rhadapt%p_IvertexAge(&
            rhadapt%p_IverticesAtElement(1:TRIA_NVETRI2D,iel))
        istate = redgreen_getStateTria(IvertexAge(1:TRIA_NVETRI2D))
      else
        IvertexAge(1:TRIA_NVEQUAD2D) = rhadapt%p_IvertexAge(&
            rhadapt%p_IverticesAtElement(1:TRIA_NVEQUAD2D,iel))
        istate = redgreen_getStateQuad(IvertexAge(1:TRIA_NVEQUAD2D))
      end if
      
    end if
  end function redgreen_getState2D

  ! ***************************************************************************

!<function>

  pure function redgreen_getStateTria(IvertexAge) result(istate)

!<description>
    ! This pure functions encodes the state of the given triangle
    ! into a unique positive integer value. If the state cannot
    ! be determined uniquely, then the position of the youngest
    ! vertex is encoded but with negative sign.
!</description>

!<input>
    ! Age of the three vertices
    integer, dimension(TRIA_NVETRI2D), intent(in) :: IvertexAge
!</input>

!<result>
    ! State of the triangle
    integer :: istate
!</result>
!</function>

    ! local variables
    integer :: ive,ipos

    ! Reset state
    istate = 0
    
    ! Check if triangle belongs to the initial triangulation. Then nothing
    ! needs to be done. Otherwise, perform additional checks
    if (sum(abs(IvertexAge)) .ne. 0) then
      
      ! Check if any two vertices have the same age and mark that edge.
      ! In addition, determine the largest age and its position
      ipos = 1
      do ive = 1, TRIA_NVETRI2D
        if (abs(IvertexAge(ive)) .eq. abs(IvertexAge(mod(ive, TRIA_NVETRI2D)+1))) then
          ! Edge connects two nodes with the same age
          istate = ibset(istate, ive)
        elseif (abs(IvertexAge(ive)) .gt. abs(IvertexAge(ipos))) then
          ! The age of the "next" node is larger, so than it might be the largest
          ipos = ive
        end if
      end do
      
      ! If ISTATE=0, then there must be one youngest vertex. Its local number
      ! is returned but with negative sign. In this case the state of the
      ! triangle cannot be uniquely determined and additional checks are required.
      !
      ! Otherwise, either all nodes have the same age (1110 : inner red triangle)
      ! or exactly two nodes have the same age. In the latter case, we must check
      ! if the opposite vertex is the youngest one in the triple of nodes.
      select case(istate)
      case(STATE_TRIA_ROOT)
        istate = -ibset(0, ipos)

      case(STATE_TRIA_OUTERINNER1)
        if (ipos .eq. 3) istate = -ibset(0, ipos)

      case(STATE_TRIA_OUTERINNER)
        if (ipos .eq. 1) istate = -ibset(0, ipos)

      case(STATE_TRIA_OUTERINNER2)
        if (ipos .eq. 2) istate = -ibset(0, ipos)

      end select
    end if

    ! Mark state for triangle
    istate = ibclr(istate,0)
  end function redgreen_getStateTria

  ! ***************************************************************************

!<function>

  pure function redgreen_getStateQuad(IvertexAge) result(istate)

!<description>
    ! This pure functions encodes the state of the given quadrilateral
    ! into a unique integer value
!</description>

!<input>
    ! Age of the four vertices
    integer, dimension(TRIA_NVEQUAD2D), intent(in) :: IvertexAge
!</input>

!<result>
    ! State of the quadrilateral
    integer :: istate
!</result>
!</function>

    ! local variables
    integer :: ive

    ! Reset state
    istate = 0

    ! Check if quadrilateral belongs to the initial triangulation. Then nothing
    ! needs to be done. Otherwise, perform additional checks.
    if (sum(abs(IvertexAge)) .ne. 0) then
      
      ! Check if any two vertices have the same age and mark that edge.
      ! After this procedure, ISTATE must be different from 0 and a unique
      ! state for the quadrilateral has been determined.
      do ive = 1, TRIA_NVEQUAD2D
        if (abs(IvertexAge(ive)) .eq. abs(IvertexAge(mod(ive, TRIA_NVEQUAD2D)+1))) then
          ! Edge connects two nodes with the same age
          istate = ibset(istate,ive)
        end if
      end do
    end if

    ! Mark state for quadrilateral
    istate = ibset(istate,0)
  end function redgreen_getStateQuad

  ! ***************************************************************************

!<function>
  
  pure function redgreen_rotateState2D(istate, irotate) result(inewstate)
    
!<description>
    ! This pure function "rotates" a given state
!</description>

!<input>
    ! Current state
    integer, intent(in) :: istate

    ! Positive rotation
    integer, intent(in) :: irotate
!</input>

!<result>
    ! New state
    integer :: inewstate
!</result>
!</function>

    ! Are we triangular or quadrilateral?
    if (btest(istate,0)) then
      inewstate = ishft(istate, -1)
      inewstate = ishftc(inewstate, irotate, 4)
      inewstate = ishft(inewstate, 1)
      inewstate = ibset(inewstate, 0)
    else
      inewstate = ishft(istate, -1)
      inewstate = ishftc(inewstate, irotate, 3)
      inewstate = ishft(inewstate, 1)
      inewstate = ibclr(inewstate, 0)
    end if
  end function redgreen_rotateState2D

end module hadaptaux2d
