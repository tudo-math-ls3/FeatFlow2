!##############################################################################
!# ****************************************************************************
!# <name> hadaptaux </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# WARNING: Do not USE this module in your applications unless you really
!#          know what you are doing. This module does no error checking!!!
!#
!# This module contains the basic data structures and auxiliary routines
!# which are required for performing h-adaptivity in.
!#
!# The following routines are available:
!#
!#  1.) hadapt_getNVE
!#      -> return the number of vertices per element
!#
!#  2.) hadapt_setVertexCoords1D
!#      -> Set the coordinates of vertices to the adaptivity structure in 1D
!#
!#  3.) hadapt_getVertexCoords1D
!#      -> Get the coordinates of vertices from the adaptivity structure in 1D
!#
!#  4.) hadapt_setVertexCoords2D
!#      -> Set the coordinates of vertices to the adaptivity structure in 2D
!#
!#  5.) hadapt_getVertexCoords2D
!#      -> Get the coordinates of vertices from the adaptivity structure in 2D
!#
!#  6.) hadapt_setVertexCoords3D
!#      -> Set the coordinates of vertices to the adaptivity structure in 3D
!#
!#  7.) hadapt_getVertexCoords3D
!#      -> Get the coordinates of vertices from the adaptivity structure in 3D
!#
!#  8.) hadapt_setVerticesAtElement
!#      -> Set the "vertices-at-element" structure to the adaptivity structure
!#
!#  9.) hadapt_getVerticesAtElement
!#      -> Get the "vertices-at-element" structure from the adaptivity structure
!#
!# 10.) hadapt_setNeighboursAtElement
!#      -> Set the "neighbours-at-element" structure to the adaptivity structure
!#
!# 11.) hadapt_getNeighboursAtElement
!#      -> Get the "neighbours-at-element" structure from the adaptivity structure
!#
!# 12.) hadapt_setNelOfType
!#      -> Set the "InelOfType" array to the adaptivity structure
!#
!# 13.) hadapt_getNelOfType
!#      -> Get the "InelOfType" array from the adaptivity structure
!#
!# 14.) hadapt_setBoundary
!#      -> Set the boundary structure to the adaptivity structure
!#
!# 15.) hadapt_getBoundary
!#      -> Get the boundary structure form the adaptivity structure
!#
!# 16.) hadapt_setNodalProperty
!#      -> Set the "nodal property" list to the adaptivity structure
!#
!# 17.) hadapt_getNodalProperty
!#      -> Get the "nodal property" list from the adaptivity structure
!#
!# 18.) hadapt_genElementsAtVertex
!#      -> Generate the "elements-at-vertex" data structure
!#
!# </purpose>
!##############################################################################

module hadaptaux

!$use omp_lib
  use arraylistInt
  use basicgeometry
  use fsystem
  use genoutput
  use linearalgebra
  use mapInt_Dble
  use octreeDble
  use quadtreeDble
  use sort
  use storage
  use triangulation

  implicit none

  private
  public :: t_hadapt
  public :: hadapt_genElementsAtVertex
  public :: hadapt_getBoundary
  public :: hadapt_getNVE
  public :: hadapt_getNeighboursAtElement
  public :: hadapt_getNelOfType
  public :: hadapt_getNodalProperty
  public :: hadapt_getVertexCoords1D
  public :: hadapt_getVertexCoords2D
  public :: hadapt_getVertexCoords3D
  public :: hadapt_getVerticesAtElement
  public :: hadapt_setBoundary
  public :: hadapt_setNeighboursAtElement
  public :: hadapt_setNelOfType
  public :: hadapt_setNodalProperty
  public :: hadapt_setVertexCoords1D
  public :: hadapt_setVertexCoords2D
  public :: hadapt_setVertexCoords3D
  public :: hadapt_setVerticesAtElement

!<constants>

!<constantblock description="Global flags for mesh refinement/coarsening">

  ! No refinement and coarsening
  integer, parameter, public :: HADAPT_NOADAPTATION = 0

  ! Red-Green refinement and coarsening strategy (R. Banks et al.)
  integer, parameter, public :: HADAPT_REDGREEN     = 1

  ! Longest edge bisection strategy (M. Rivara)
  integer, parameter, public :: HADAPT_LONGESTEDGE  = 2

!</constantblock>


!<constantblock description="Bitfield identifiers for state of adaptation">

  ! Adaptation is undefined
  integer(I32), parameter, public :: HADAPT_UNDEFINED        = 2**0
  
  ! Parameters of adaptivity structure are initialised
  integer(I32), parameter, public :: HADAPT_HAS_PARAMETERS   = 2**1

  ! Map/quadtree/octree for vertex coordinates is generated
  integer(I32), parameter, public :: HADAPT_HAS_COORDS       = 2**2

  ! Array for IverticesAtElement is generated
  integer(I32), parameter, public :: HADAPT_HAS_VERTATELEM   = 2**3

  ! Array for IneighboursAtElement is generated
  integer(I32), parameter, public :: HADAPT_HAS_NEIGHATELEM  = 2**4

  ! Array for ImidneighboursAtElement is generated
  integer(I32), parameter, public :: HADAPT_HAS_MIDNEIGH     = 2**5

  ! Boundary data is generated
  integer(I32), parameter, public :: HADAPT_HAS_BOUNDARY     = 2**6

  ! Nodal property is generated
  integer(I32), parameter, public :: HADAPT_HAS_NODALPROP    = 2**7

  ! Number of elements for predefined type
  integer(I32), parameter, public :: HADAPT_HAS_NELOFTYPE    = 2**8

  ! Array for IelementsAtVertex is generated
  integer(I32), parameter, public :: HADAPT_HAS_ELEMATVERTEX = 2**9

  ! Dynamic data structures in 1D are all generated
  integer(I32), parameter, public :: HADAPT_HAS_DYNAMICDATA1D  = HADAPT_HAS_PARAMETERS+&
                                                                 HADAPT_HAS_COORDS+&
                                                                 HADAPT_HAS_VERTATELEM+&
                                                                 HADAPT_HAS_NEIGHATELEM+&
                                                                 HADAPT_HAS_NODALPROP+&
                                                                 HADAPT_HAS_NELOFTYPE+&
                                                                 HADAPT_HAS_ELEMATVERTEX

  ! Dynamic data structures in 2D are all generated
  integer(I32), parameter, public :: HADAPT_HAS_DYNAMICDATA2D  = HADAPT_HAS_PARAMETERS+&
                                                                 HADAPT_HAS_COORDS+&
                                                                 HADAPT_HAS_VERTATELEM+&
                                                                 HADAPT_HAS_NEIGHATELEM+&
                                                                 HADAPT_HAS_MIDNEIGH+&
                                                                 HADAPT_HAS_BOUNDARY+&
                                                                 HADAPT_HAS_NODALPROP+&
                                                                 HADAPT_HAS_NELOFTYPE+&
                                                                 HADAPT_HAS_ELEMATVERTEX
  
  ! Dynamic data structures in 3D are all generated
  integer(I32), parameter, public :: HADAPT_HAS_DYNAMICDATA3D  = HADAPT_HAS_PARAMETERS+&
                                                                 HADAPT_HAS_COORDS+&
                                                                 HADAPT_HAS_VERTATELEM+&
                                                                 HADAPT_HAS_NEIGHATELEM+&
                                                                 HADAPT_HAS_NODALPROP+&
                                                                 HADAPT_HAS_NELOFTYPE+&
                                                                 HADAPT_HAS_ELEMATVERTEX
  
  ! Cells are marked for refinement
  integer(I32), parameter, public :: HADAPT_MARKEDREFINE     = 2**10

  ! Cells are marked for coarsening
  integer(I32), parameter, public :: HADAPT_MARKEDCOARSEN    = 2**11

  ! Cells are marked
  integer(I32), parameter, public :: HADAPT_MARKED           = HADAPT_MARKEDREFINE+&
                                                               HADAPT_MARKEDCOARSEN
  
!</constantblock>


!<constantblock description="Global constants for grid modification operations">

  ! Operation identifier for initialization of callback function
  integer, parameter, public :: HADAPT_OPR_INITCALLBACK     = -1

  ! Operation identifier for finalization of callback function
  integer, parameter, public :: HADAPT_OPR_DONECALLBACK     = -2

  ! Operation identifier for adjustment of vertex dimension
  integer, parameter, public :: HADAPT_OPR_ADJUSTVERTEXDIM  = 1

  ! Operation identifier for vertex insertion at edge midpoint
  integer, parameter, public :: HADAPT_OPR_INSERTVERTEXEDGE = 2

  ! Operation identifier for vertex insertion at element centroid
  integer, parameter, public :: HADAPT_OPR_INSERTVERTEXCENTR= 3

  ! Operation identifier for vertex removal
  integer, parameter, public :: HADAPT_OPR_REMOVEVERTEX     = 4

  ! Operation identifier for refinment: 1-tria : 2-tria
  integer, parameter, public :: HADAPT_OPR_REF_TRIA2TRIA    = 5

  ! Operation identifier for refinment: 1-tria : 3-tria
  integer, parameter, public :: HADAPT_OPR_REF_TRIA3TRIA12  = 6
  integer, parameter, public :: HADAPT_OPR_REF_TRIA3TRIA23  = 7
 
  ! Operation identifier for refinment: 1-tria : 4-tria
  integer, parameter, public :: HADAPT_OPR_REF_TRIA4TRIA    = 8

  ! Operation identifier for refinment: 1-quad : 2-quad
  integer, parameter, public :: HADAPT_OPR_REF_QUAD2QUAD    = 9
  
  ! Operation identifier for refinment: 1-quad : 3-tria
  integer, parameter, public :: HADAPT_OPR_REF_QUAD3TRIA    = 10

  ! Operation identifier for refinment: 1-quad : 4-tria
  integer, parameter, public :: HADAPT_OPR_REF_QUAD4TRIA    = 11

  ! Operation identifier for refinment: 1-quad : 4-quad
  integer, parameter, public :: HADAPT_OPR_REF_QUAD4QUAD    = 12

  ! Operation identifier for conversion: 2-tria : 4-tria
  integer, parameter, public :: HADAPT_OPR_CVT_TRIA2TRIA    = 13

  ! Operation identifier for conversion: 2-quad : 4-tria
  integer, parameter, public :: HADAPT_OPR_CVT_QUAD2QUAD    = 14

  ! Operation identifier for conversion: 3-tria : 4-quad
  integer, parameter, public :: HADAPT_OPR_CVT_QUAD3TRIA    = 15

  ! Operation identifier for conversion: 4-tria : 4-quad
  integer, parameter, public :: HADAPT_OPR_CVT_QUAD4TRIA    = 16

  ! Operation identifier for coarsening: 2-tria : 1-tria
  integer, parameter, public :: HADAPT_OPR_CRS_2TRIA1TRIA   = 17

  ! Operation identifier for coarsening: 4-tria : 1-tria
  integer, parameter, public :: HADAPT_OPR_CRS_4TRIA1TRIA   = 18

  ! Operation identifier for coarsening: 4-tria : 2-tria
  integer, parameter, public :: HADAPT_OPR_CRS_4TRIA2TRIA1  = 19
  integer, parameter, public :: HADAPT_OPR_CRS_4TRIA2TRIA2  = 20
  integer, parameter, public :: HADAPT_OPR_CRS_4TRIA2TRIA3  = 21

  ! Operation identifier for coarsening: 4-quad : 1-quad
  integer, parameter, public :: HADAPT_OPR_CRS_4QUAD1QUAD   = 22

  ! Operation identifier for coarsening: 4-quad : 2-quad
  integer, parameter, public :: HADAPT_OPR_CRS_4QUAD2QUAD   = 23

  ! Operation identifier for coarsening: 4-quad : 3-tria
  integer, parameter, public :: HADAPT_OPR_CRS_4QUAD3TRIA   = 24

  ! Operation identifier for coarsening: 4-quad : 4-tria
  integer, parameter, public :: HADAPT_OPR_CRS_4QUAD4TRIA   = 25

  ! Operation identifier for coarsening: 2-quad : 1-quad
  integer, parameter, public :: HADAPT_OPR_CRS_2QUAD1QUAD   = 26

  ! Operation identifier for coarsening: 2-quad : 3-tria
  integer, parameter, public :: HADAPT_OPR_CRS_2QUAD3TRIA   = 27

  ! Operation identifier for coarsening: 3-tria : 1-quad
  integer, parameter, public :: HADAPT_OPR_CRS_3TRIA1QUAD   = 28

  ! Operation identifier for coarsening: 4-tria : 1-quad
  integer, parameter, public :: HADAPT_OPR_CRS_4TRIA1QUAD   = 29

  ! Operation identifier for coarsening: 4-tria : 3-tria
  integer, parameter, public :: HADAPT_OPR_CRS_4TRIA3TRIA2  = 30
  integer, parameter, public :: HADAPT_OPR_CRS_4TRIA3TRIA3  = 31

  ! Operation identifier for refinement: 1-line : 2-line
  integer, parameter, public :: HADAPT_OPR_REF_LINE2LINE    = 32
  
  ! Operation identifier for coarsening: 2-line : 1-line
  integer, parameter, public :: HADAPT_OPR_CRS_2LINE1LINE   = 33

!</constantblock>

!<constantblock
! description="Duplication flags. Specifies which information is shared between adaptivity structures">

  integer(I32), parameter, public :: HADAPT_SHARE_IMARKER            = 2** 0
  integer(I32), parameter, public :: HADAPT_SHARE_IVERTEXAGE         = 2** 1
  integer(I32), parameter, public :: HADAPT_SHARE_INODALPROPERTY     = 2** 2
  integer(I32), parameter, public :: HADAPT_SHARE_IVERTICESATELEMENT = 2** 3
  integer(I32), parameter, public :: HADAPT_SHARE_INEIGHATELEMENT    = 2** 4
  integer(I32), parameter, public :: HADAPT_SHARE_IMIDNEIGHATELEMENT = 2** 5
  integer(I32), parameter, public :: HADAPT_SHARE_RVERTEXCOORDINATES = 2** 6
  integer(I32), parameter, public :: HADAPT_SHARE_RBOUNDARY          = 2** 7
  integer(I32), parameter, public :: HADAPT_SHARE_RELEMENTSATVERTEX  = 2** 8

!</constantblock>

!</constants>

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<types>

  !<typeblock>
  
  ! This type contains all data structures to handle
  ! adaptive grid refinement and grid coarsening.
  type :: t_hadapt
    ! Format Tag: Specifies the state of adaptation
    integer(I32) :: iSpec = HADAPT_UNDEFINED

    ! Duplication flag. Bitfield that indicates which information is
    ! shared with another adaptivity structure.
    ! When a bit is set to 1, the corresponding array is
    ! maintained by another adaptivity structure and must
    ! not be deleted by hadapt_releaseAdaptation.
    ! When the bit is 0, the array is a real copy of another array
    ! and must be deleted in hadapt_releaseAdaptation.
    integer(I32) :: iduplicationFlag = 0

    ! Tag: Specifies the strategy for grid refinement and coarsening
    integer :: iadaptationStrategy = HADAPT_NOADAPTATION
    
    ! Maximum number of subdivisions from the original mesh
    integer :: nsubdividemax = 0

    ! Total number of grid refinement steps
    integer :: nRefinementSteps = 0

    ! Total number of grid coarsening steps
    integer :: nCoarseningSteps = 0

    ! Total number of grid smoothing steps
    integer :: nSmoothingSteps = 0

    ! Tolerance for refinement
    real(DP) :: drefinementTolerance = 0

    ! Tolerance for coarsening
    real(DP) :: dcoarseningTolerance = 0

    ! Dimension of the triangulation
    integer :: ndim = 0

    ! Total number of vertices (initially)
    integer :: NVT0 = 0
    
    ! Total number of vertices
    integer :: NVT = 0

    ! Increment of vertices
    integer :: increaseNVT = 0

    ! Total number of boundary vertives (initially)
    integer :: NVBD0 = 0

    ! Total number of boundary vertives
    integer :: NVBD = 0
    
    ! Total number of boundary components (should not change)
    integer :: NBCT = 0
    
    ! Total number of elements (initially)
    integer :: NEL0 = 0
    
    ! Total number of elements
    integer :: NEL = 0

    ! Maximum number of elements (before reallocation)
    integer :: NELMAX = 0
    
    ! Total number of green elements (required internally)
    integer :: nGreenElements = 0
    
    ! Nuber of elements with a defined number of vertices per element.
    ! InelOfType(TRIA_NVELINE1D) = number of lines in the mesh (1D).
    ! InelOfType(TRIA_NVETRI2D)  = number of triangles in the mesh (2D).
    ! InelOfType(TRIA_NVEQUAD2D) = number of quadrilaterals in the mesh (2D).
    integer, dimension(TRIA_MAXNVE) :: InelOfType = 0

    ! Same as InelOfType but this array stores the number of elements
    ! which are initially present in the mesh
    integer, dimension(TRIA_MAXNVE) :: InelOfType0 = 0
    
    ! Element marker array.
    ! Handle to
    !       p_Imarker = array [1..NEL] of integer.
    ! For each element there is one identifier [0,..,31] in the marker.
    ! BIT0: Is 0 for triangle and 1 for quadrilateral
    ! BIT1: Is 0 if first edge is not subdivided, 1 otherwise.
    ! BIT2: Is 0 if second edge is not subdivided, 1 otherwise.
    ! BIT3: Is 0 if third edge is not subdivided, 1 otherwise.
    ! BIT4: Is 0 if fourth edge is not subdivided, 1 otherwise.
    integer :: h_Imarker = ST_NOHANDLE
    
    ! Vertex age array.
    ! Handle to
    !       p_IvertexAge = array [1..NVT] of integer
    ! Each node has an individual age which is used to determine
    ! the type of the element (green/red triangle/quadrilateral).
    ! The nodes of the initial grid are assigned age 0 by default.
    ! Each node which is introduced at the midpoint of one edge (I,J)
    ! is given the age MAX(AGE(I),AGE(J))+1. For vertices which are
    ! inserted at the center of a quadrilateral, the maximum age of
    ! all four surrounding vertices is adopted and increased by one.
    ! A negative age means, the the vertex is lockes, that is, it
    ! cannot be removed. By definition, vertices of the initial
    ! triangulation cannot be removed and are always locked.
    integer :: h_IvertexAge = ST_NOHANDLE

    ! Pointer to h_IvertexAge.
    ! This array is introduced to increase performance and must
    ! not be modified by the user
    integer, dimension(:), pointer :: p_IvertexAge => null()

    ! Nodal property array.
    ! Handle to
    !       p_InodalProperty = array [1..NVT+NMT] of integer
    ! p_InodalProperty(i) defines for each vertex i=(1..NVT)
    ! and each edge i=(NVT+1..NVT+NMT) its function inside of the
    ! geometry. Generally said, the range of the p_InodalProperty-array
    ! characterises the type of the node (=vertex/edge):
    ! = 0    : The vertex/edge is an inner vertex/edge
    ! > 0    : The vertex/edge is a boundary vertex/edge on the real
    !           boundary. KNPR(.) defines the number of the boundary
    !           component.
    integer :: h_InodalProperty = ST_NOHANDLE

    ! Pointer to h_InodalProperty.
    ! This array is introduced to increase performance and must
    ! not be modified by the user
    integer, dimension(:), pointer :: p_InodalProperty => null()
    
    ! Vertices adjacent to an element.
    ! Handle to
    !       p_IverticesAtElement = array [1..TRIA_MAXNVE2D,1..NEL] of integer.
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    integer :: h_IverticesAtElement = ST_NOHANDLE

    ! Pointer to h_IverticesAtElement.
    ! This array is introduced to increase performance and must
    ! not be modified by the user
    integer, dimension(:,:), pointer :: p_IverticesAtElement => null ()

    ! Neighbour elements adjacent to an element.
    ! Handle to
    !       p_IneighboursAtElement = array [1..TRIA_MAXNME2D,1..NEL] of integer
    ! For each element, the numbers of adjacent elements in mathematically
    ! positive sense, metting the element in an edge.
    integer :: h_IneighboursAtElement = ST_NOHANDLE

    ! Pointer to h_IneighboursAtElement.
    ! This array is introduced to increase performance and must
    ! not be modified by the user
    integer, dimension(:,:), pointer :: p_IneighboursAtElement => null ()

    ! Midneighbour elements adjacent to an element.
    ! Handle to
    !       p_ImidneighboursAtElement = array [1..TRIA_MAXNME2D,1..NEL] of integer
    ! H-adaptivity is performed element-by-element. If one element is
    ! refined, then its neighbors need to be informed that there
    ! are possibly two elements adjacent along one edge. This is a
    ! nonconforming state. When time comes to process the other
    ! element (with hanging node), the element knows which elements
    ! are adjacent along the first and the second half of the edge.
    integer :: h_ImidneighboursAtElement = ST_NOHANDLE

    ! Pointer to h_ImidneighboursAtElement.
    ! This array is introduced to increase performance and must
    ! not be modified by the user
    integer, dimension(:,:), pointer :: p_ImidneighboursAtElement => null ()
    
    ! Handle to h_DVertexCoords in 1D
    integer :: h_DvertexCoords1D

    ! Pointer to handle h_DvertexCoords in 1D
    real(DP), dimension(:,:), pointer :: p_DvertexCoords1D => null()
    
    ! Quadtree storing the nodal coordinates in 2D
    type(t_quadtreeDble) :: rVertexCoordinates2D

    ! Octree storing the nodal coordinates in 2D
    type(t_octreeDble) :: rVertexCoordinates3D
    
    ! Array of maps storing the boundary data
    ! p_IboundaryCpIdx and p_IverticesAtBoundary
    type(t_mapInt_Dble), dimension(:), pointer :: rBoundary => null()

    ! Arraylist for elements-meeting-at-vertex structure
    type(t_arraylistInt) :: rElementsAtVertex
  end type t_hadapt

  !</typeblock>

!</types>



contains

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<function>

  pure function hadapt_getNVE(rhadapt, iel) result(nve)

!<description>
    ! This function returns the number of vertices present in the given element
!</description>

!<input>
    ! Adaptivity structure
    type(t_hadapt), intent(in) :: rhadapt

    ! Number of the element
    integer, intent(in) :: iel
!</input>

!<result>
    ! number of vertices per element
    integer :: nve
!</result>
!</function>

    ! Which spatial dimension do we have?
    select case(rhadapt%ndim)
    case (NDIM1D)
      nve = TRIA_NVELINE1D

    case (NDIM2D)
      
      ! Do we have quadrilaterals in the triangulation?
      if (rhadapt%InelOfType(TRIA_NVEQUAD2D) .eq. 0) then
        
        ! There are no quadrilaterals in the current triangulation.
        ! Hence, return TRIA_NVETRI2D by default.
        nve = TRIA_NVETRI2D
        
      else
        
        ! There are quadrilaterals and possible also triangles in
        ! the current triangulatin. If the last entry of the
        ! vertices-at-element list is nonzero then TRIA_NVEQUAD2D vertices
        ! are present in the current element. Otherwise return TRIA_NVETRI2D.
        if (rhadapt%p_IverticesAtElement(TRIA_NVEQUAD2D, iel) .eq. 0) then
          nve = TRIA_NVETRI2D
        else
          nve = TRIA_NVEQUAD2D
        end if
        
      end if

    case DEFAULT
      nve = 0
    end select
  end function hadapt_getNVE

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_setVertexCoords1D(rhadapt, h_DvertexCoords, nvt)

!<description>
    ! This subroutine sets the vertex coordinates given by the handle
    ! h_DvertexCoords that points to the two-dimensional array
    ! p_DvertexCoords and stores them in the map.
!</description>

!<input>
    ! Handle to the vertex coordinates
    integer, intent(in) :: h_DvertexCoords

    ! Number of vertices
    integer, intent(in) :: nvt
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    ! Check if handle is not empty
    if (h_DvertexCoords .eq. ST_NOHANDLE) then
      call output_line('Invalid handle!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setVertexCoords1D')
      call sys_halt()
    end if

    ! Check if structure is already generated, then remove old structure first
    if (iand(rhadapt%iSpec, HADAPT_HAS_COORDS) .eq.&
                            HADAPT_HAS_COORDS) then
      if (rhadapt%h_DvertexCoords1D .ne. ST_NOHANDLE)&
          call storage_free(rhadapt%h_DvertexCoords1D)
    end if

    ! Copy handle and set pointer
    call storage_copy(h_DvertexCoords, rhadapt%h_DvertexCoords1D)
    call storage_getbase_double2D(rhadapt%h_DvertexCoords1D,&
                                  rhadapt%p_DvertexCoords1D)

    ! Set specifier for structure
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_HAS_COORDS)

    ! Set dimensions
    rhadapt%ndim = NDIM1D
    rhadapt%NVT  = nvt
  end subroutine hadapt_setVertexCoords1D

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_getVertexCoords1D(rhadapt, h_DvertexCoords, nvt, ndim)

!<description>
    ! This subroutine gets the vertex coordinates from the structure
    ! and stores them in the two-dimensional array associated to the
    ! handle h_DvertexCoords. Note that the allocated memory will
    ! be reallocated if is does not provide the correct dimensions.
!</description>

!<input>
    ! Adaptivity structure
    type(t_hadapt), intent(in) :: rhadapt
!</input>

!<inputoutput>
    ! Handlt to the vertex coordinate vector
    integer, intent(inout) :: h_DvertexCoords
!</inputoutput>

!<output>
    ! OPTIONAL: number of vertices
    integer, intent(out), optional :: nvt

    ! OPTIONAL: number of spatial dimensions
    integer, intent(out), optional :: ndim
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:,:), pointer :: p_Ddata
    integer, dimension(2) :: Isize

    ! Check if coordinates exists
    if (iand(rhadapt%iSpec, HADAPT_HAS_COORDS) .ne.&
                            HADAPT_HAS_COORDS) then
      call output_line('Coords do not exist!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getVertexCoords1D')
      call sys_halt()
    end if

    ! Check if coordinates are given in 1D
    if (rhadapt%ndim .ne. NDIM1D) then
      call output_line('Invalid spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getVertexCoords1D')
      call sys_halt()
    end if

    ! Check if given handle is valid
    if (h_DvertexCoords .eq. ST_NOHANDLE) then
      Isize=(/1,rhadapt%NVT/)
      call storage_new('hadapt_getVertexCoords1D', 'DCORVG', Isize,&
                       ST_DOUBLE, h_DvertexCoords, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_DvertexCoords, Isize)
      if (Isize(1) .ne. 1) then
        call output_line('First component of DvertexCoords must be 1!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getVertexCoords1D')
        call sys_halt()
      end if

      if (Isize(2) .lt. rhadapt%NVT) then
        call storage_realloc('hadapt_getVertexCoords1D', rhadapt%NVT,&
                             h_DvertexCoords, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if

    ! Set pointers
    call storage_getbase_double2D(rhadapt%h_DvertexCoords1D,&
                                  p_Ddata, rhadapt%NVT)
    call storage_getbase_double2D(h_DvertexCoords, p_DvertexCoords)

    ! Copy data
    call lalg_copyVectorDble2D(p_Ddata, p_DvertexCoords, 1, rhadapt%NVT)

    ! Set dimension
    if (present(ndim)) ndim = rhadapt%ndim
    if (present(nvt))  nvt  = rhadapt%NVT
  end subroutine hadapt_getVertexCoords1D

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_setVertexCoords2D(rhadapt, h_DvertexCoords, nvt)

!<description>
    ! This subroutine sets the vertex coordinates given by the handle
    ! h_DvertexCoords that points to the two-dimensional array
    ! p_DvertexCoords and stores them in the quadtree.
!</description>

!<input>
    ! Handle to the vertex coordinates
    integer, intent(in) :: h_DvertexCoords

    ! Number of vertices
    integer, intent(in) :: nvt
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP) :: xmin,xmax,ymin,ymax
    integer  :: nnode

    ! Check if handle is not empty
    if (h_DvertexCoords .eq. ST_NOHANDLE) then
      call output_line('Invalid handle!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setVertexCoords2D')
      call sys_halt()
    end if

    ! Check if quadtree is already generated, then remove old quadtree/octree first
    if (iand(rhadapt%iSpec, HADAPT_HAS_COORDS) .eq.&
                            HADAPT_HAS_COORDS) then
      call qtree_release(rhadapt%rVertexCoordinates2D)
    end if
    
    ! Set pointer
    call storage_getbase_double2D(h_DvertexCoords, p_DvertexCoords, nvt)

    ! Get outer bounding-box of vertex coordinates
    xmin = minval(p_DvertexCoords(1,:))
    xmax = maxval(p_DvertexCoords(1,:))
    ymin = minval(p_DvertexCoords(2,:))
    ymax = maxval(p_DvertexCoords(2,:))
    
    ! Estimate number of initial quadrilaterals
    nnode = int(0.5_DP*nvt)

    ! Create quadtree for vertices
    call qtree_create(rhadapt%rVertexCoordinates2D, nvt,&
                      nnode, xmin, ymin, xmax, ymax)
    
    ! Copy vertex coordinates to quadtree
    call qtree_copy(p_DvertexCoords, rhadapt%rVertexCoordinates2D)

    ! Set specifier for quadtree
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_HAS_COORDS)

    ! Set dimensions
    rhadapt%ndim = NDIM2D
    rhadapt%NVT  = nvt
  end subroutine hadapt_setVertexCoords2D

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_getVertexCoords2D(rhadapt, h_DvertexCoords, nvt, ndim)

!<description>
    ! This subroutine gets the vertex coordinates from the quadtree
    ! and stores them in the two-dimensional array associated to the
    ! handle h_DvertexCoords. Note that the allocated memory will
    ! be reallocated if is does not provide the correct dimensions.
!</description>

!<input>
    ! Adaptivity structure
    type(t_hadapt), intent(in) :: rhadapt
!</input>

!<inputoutput>
    ! Handlt to the vertex coordinate vector
    integer, intent(inout) :: h_DvertexCoords
!</inputoutput>

!<output>
    ! OPTIONAL: number of vertices
    integer, intent(out), optional :: nvt

    ! OPTIONAL: number of spatial dimensions
    integer, intent(out), optional :: ndim
!</output>
!</subroutine>

    ! Check if coordinates exists
    if (iand(rhadapt%iSpec, HADAPT_HAS_COORDS) .ne.&
                            HADAPT_HAS_COORDS) then
      call output_line('Quadtree does not exist!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getVertexCoords2D')
      call sys_halt()
    end if

    ! Check if coordinates are given in 2D
    if (rhadapt%ndim .ne. NDIM2D) then
      call output_line('Invalid spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getVertexCoords2D')
      call sys_halt()
    end if

    ! Copy quadtree to handle h_DvertexCoords
    call qtree_copy(rhadapt%rVertexCoordinates2D, h_DvertexCoords)

    ! Set dimension
    if (present(ndim)) ndim = rhadapt%ndim
    if (present(nvt))  nvt  = rhadapt%NVT
  end subroutine hadapt_getVertexCoords2D
  
  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_setVertexCoords3D(rhadapt, h_DvertexCoords, nvt)

!<description>
    ! This subroutine sets the vertex coordinates given by the handle
    ! h_DvertexCoords that points to the three-dimensional array
    ! p_DvertexCoords and stores them in the octree.
!</description>

!<input>
    ! Handle to the vertex coordinates
    integer, intent(in) :: h_DvertexCoords

    ! Number of vertices
    integer, intent(in) :: nvt
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout)    :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP) :: xmin,xmax,ymin,ymax,zmin,zmax
    integer  :: nnode

    ! Check if handle is not empty
    if (h_DvertexCoords .eq. ST_NOHANDLE) then
      call output_line('Invalid handle!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setVertexCoords3D')
      call sys_halt()
    end if

    ! Check if quadtree is already generated, then remove old quadtree/octree first
    if (iand(rhadapt%iSpec, HADAPT_HAS_COORDS) .eq.&
                            HADAPT_HAS_COORDS) then
      call otree_release(rhadapt%rVertexCoordinates3D)
    end if
    
    ! Set pointer
    call storage_getbase_double2D(h_DvertexCoords, p_DvertexCoords, nvt)

    ! Get outer bounding-box of vertex coordinates
    xmin = minval(p_DvertexCoords(1,:))
    xmax = maxval(p_DvertexCoords(1,:))
    ymin = minval(p_DvertexCoords(2,:))
    ymax = maxval(p_DvertexCoords(2,:))
    zmin = minval(p_DvertexCoords(3,:))
    zmax = maxval(p_DvertexCoords(3,:))
    
    ! Estimate number of initial quadrilaterals
    nnode = int(0.5_DP*nvt)

    ! Create octree for vertices
    call otree_create(rhadapt%rVertexCoordinates3D, nvt,&
                            nnode, xmin, ymin, zmin, xmax, ymax, zmax)
    
    ! Copy vertex coordinates to octree
    call otree_copy(p_DvertexCoords, rhadapt%rVertexCoordinates3D)

    ! Set specifier for quadtree
    rhadapt%iSpec=ior(rhadapt%iSpec, HADAPT_HAS_COORDS)

    ! Set dimensions
    rhadapt%ndim = NDIM3D
    rhadapt%NVT  = nvt
  end subroutine hadapt_setVertexCoords3D

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_getVertexCoords3D(rhadapt, h_DvertexCoords, nvt, ndim)

!<description>
    ! This subroutine gets the vertex coordinates from the octree
    ! and stores them in the three-dimensional array associated to the
    ! handle h_DvertexCoords. Note that the allocated memory will
    ! be reallocated if is does not provide the correct dimensions.
!</description>

!<input>
    ! Adaptivity structure
    type(t_hadapt), intent(in) :: rhadapt
!</input>

!<inputoutput>
    ! Handlt to the vertex coordinate vector
    integer, intent(inout) :: h_DvertexCoords
!</inputoutput>

!<output>
    ! OPTIONAL: number of vertices
    integer, intent(out), optional :: nvt

    ! OPTIONAL: number of spatial dimensions
    integer, intent(out), optional :: ndim
!</output>
!</subroutine>

    ! Check if coordinates exists
    if (iand(rhadapt%iSpec, HADAPT_HAS_COORDS) .ne.&
                            HADAPT_HAS_COORDS) then
      call output_line('Octree does not exist!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getVertexCoords3D')
      call sys_halt()
    end if

    ! Check if coordinates are given in 2D
    if (rhadapt%ndim .ne. NDIM3D) then
      call output_line('Invalid spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getVertexCoords3D')
      call sys_halt()
    end if

    ! Copy octree to handle h_DvertexCoords
    call otree_copy(rhadapt%rVertexCoordinates3D, h_DvertexCoords)

    ! Set dimension
    if (present(ndim)) ndim = rhadapt%ndim
    if (present(nvt))  nvt  = rhadapt%NVT
  end subroutine hadapt_getVertexCoords3D

  ! ***************************************************************************

!<subroutine>
  
  subroutine hadapt_setVerticesAtElement(rhadapt, h_IverticesAtElement, nel)

!<description>
    ! This routine assigns the handle to the "vertices-at-element" array
    ! to the adaptivity structure and sets up internal links.
!</description>

!<input>
    ! Handle to p_IverticesAtElement
    integer, intent(in) :: h_IverticesAtElement

    ! Total number of elements
    integer, intent(in) :: nel
!</input>

!<inputoutput>
    ! adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    ! Check if handle is not empty
    if (h_IverticesAtElement .eq. ST_NOHANDLE) then
      call output_line('Invalid handle!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setVerticesAtElement')
      call sys_halt()
    end if

    rhadapt%h_IverticesAtElement = h_IverticesAtElement
    call storage_getbase_int2D(rhadapt%h_IverticesAtElement,&
                               rhadapt%p_IverticesAtElement)
    
    ! Set specifier for IverticesAtElement
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_HAS_VERTATELEM)

    ! Set dimensions
    rhadapt%NEL = nel
  end subroutine hadapt_setVerticesAtElement

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_getVerticesAtElement(rhadapt, h_IverticesAtElement, nel)

!<description>
    ! This routine assigns the "vertices-at-element" array from the
    ! adaptivity structure to the given handle.
!</description>

!<input>
    ! Adaptivity structure
    type(t_hadapt), intent(in) :: rhadapt
!</input>

!<inputoutput>
    ! Hande to p_IverticesAtElement
    integer, intent(inout) :: h_IverticesAtElement
!</inputoutput>

!<output>
    ! OPTIONAL: number of elements
    integer, intent(out), optional :: nel
!</output>
!</subroutine>

    ! Check if "vertices-at-element" array exists
    if (iand(rhadapt%iSpec, HADAPT_HAS_VERTATELEM) .ne.&
                            HADAPT_HAS_VERTATELEM) then
      call output_line('Structure does not exist!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getVerticesAtElement')
      call sys_halt()
    end if

    ! Check if handle needs to be freed first
    if ((h_IverticesAtElement .ne. ST_NOHANDLE) .and.&
        (h_IverticesAtElement .ne. rhadapt%h_IverticesAtElement))&
        call storage_free(h_IverticesAtElement)

    ! Assign handle
    h_IverticesAtElement = rhadapt%h_IverticesAtElement

    ! Set dimensions
    if(present(nel))   nel = rhadapt%NEL
  end subroutine hadapt_getVerticesAtElement
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine hadapt_setNeighboursAtElement(rhadapt, h_IneighboursAtElement)

!<description>
    ! This routine assigns the handle to the "neighbours-at-element" array
    ! to the adaptivity structure and sets up internal links.
!</description>

!<input>
    ! Handle to p_IneighboursAtElement
    integer, intent(in) :: h_IneighboursAtElement
!</input>

!<inputoutput>
    ! adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    ! Check if handle is not empty
    if (h_IneighboursAtElement .eq. ST_NOHANDLE) then
      call output_line('Invalid handle!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setNeighboursAtElement')
      call sys_halt()
    end if

    rhadapt%h_IneighboursAtElement = h_IneighboursAtElement
    call storage_getbase_int2D(rhadapt%h_IneighboursAtElement,&
                               rhadapt%p_IneighboursAtElement)
    
    ! What spatial dimension are we?
    select case(rhadapt%ndim)
    case (NDIM1D,&
          NDIM3D)
      ! Do nothing

    case(NDIM2D)
      ! Create structure of "mid-adjacent" neighbours
      call storage_copy(rhadapt%h_IneighboursAtElement,&
                        rhadapt%h_ImidneighboursAtElement)
      call storage_getbase_int2D(rhadapt%h_ImidneighboursAtElement,&
                                 rhadapt%p_ImidneighboursAtElement)

      ! Set specifier for ImidneighboursAtElement
      rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_HAS_MIDNEIGH)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setNeighboursAtElement')
      call sys_halt()
    end select
    
    ! Set specifier for IverticesAtElement
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_HAS_NEIGHATELEM)
  end subroutine hadapt_setNeighboursAtElement

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_getNeighboursAtElement(rhadapt,&
                                           h_IneighboursAtElement, nel)

!<description>
    ! This routine assigns the "neighbours-at-element" array from the
    ! adaptivity structure to the given handle
!</description>

!<input>
    ! Adaptivity structure
    type(t_hadapt), intent(in) :: rhadapt
!</input>

!<inputoutput>
    ! Hande to p_IneighboursAtElement
    integer, intent(inout) :: h_IneighboursAtElement
!</inputoutput>

!<output>
    ! OPTIONAL: number of elements
    integer, intent(out), optional :: nel
!</output>
!</subroutine>

    ! Check if "neighbours-at-element" array exists
    if (iand(rhadapt%iSpec, HADAPT_HAS_NEIGHATELEM) .ne.&
                            HADAPT_HAS_NEIGHATELEM) then
      call output_line('Structure does not exist!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getNeighboursAtElement')
      call sys_halt()
    end if

    ! Check if handle needs to be freed first
    if ((h_IneighboursAtElement .ne. ST_NOHANDLE) .and.&
        (h_IneighboursAtElement .ne. rhadapt%h_IneighboursAtElement))&
        call storage_free(h_IneighboursAtElement)

    ! Assign handle
    h_IneighboursAtElement = rhadapt%h_IneighboursAtElement

    ! Set dimensions
    if(present(nel))   nel = rhadapt%NEL
  end subroutine hadapt_getNeighboursAtElement

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_setNelOfType(rhadapt, InelOfType)

!<description>
    ! This subroutine sets the number of elements with a defined number
    ! of vertices per elements, that is, InelOfType(TRIA_MAXNVE)
!</description>

!<input>
    ! Number of elements with a defined number of vertices per element.
    integer, dimension(TRIA_MAXNVE), intent(in) :: InelOfType
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    rhadapt%InelOfType = InelOfType

    ! Set specifier for InelOfType
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_HAS_NELOFTYPE)
  end subroutine hadapt_setNelOfType

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_getNelOfType(rhadapt, InelOfType)

!<description>
    ! This subroutine returns the number of elements with a defined number
    ! of vertices per elements, that is, InelOfType(TRIA_MAXNVE)
!</description>

!<input>
    ! Adaptivity structure
    type(t_hadapt), intent(in) :: rhadapt
!</input>

!<output>
    ! Number of elements with a defined number of vertices per element.
    integer, dimension(TRIA_MAXNVE), intent(out) :: InelOfType
!</output>
!</subroutine>

    ! Check if "InelOfType" array exists
    if (iand(rhadapt%iSpec, HADAPT_HAS_NELOFTYPE) .ne.&
                            HADAPT_HAS_NELOFTYPE) then
      call output_line('Structure does not exist!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getNelOfType')
      call sys_halt()
    end if
    
    InelOfType = rhadapt%InelOfType
  end subroutine hadapt_getNelOfType

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_setBoundary(rhadapt, h_IboundaryCpIdx, h_IverticesAtBoundary,&
                                h_DvertexParameterValue, nbct, nvbd)

!<description>
    ! This subroutine sets the boundary structure to the
    ! adaptivity structure and initialises internal data
!</description>

!<input>
    ! Handle to p_IboundaryCpIdx
    integer, intent(in) :: h_IboundaryCpIdx

    ! Handle to p_IverticesAtBoundary
    integer, intent(in) :: h_IverticesAtBoundary

    ! Handle to p_DvertexParameterValue
    integer, intent(in) :: h_DvertexParameterValue

    ! Number of boundary components
    integer, intent(in) :: nbct

    ! Number of vertices at the boundary
    integer, intent(in) :: nvbd
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>
    
    ! local variables
    type(it_mapInt_Dble) :: rmapIter
    real(DP), dimension(:), pointer   :: p_DvertexParameterValue
    integer, dimension(:), pointer    :: p_IboundaryCpIdx
    integer, dimension(:), pointer    :: p_IverticesAtBoundary
    real(DP), dimension(3) :: Ddata
    integer :: ibct,ioff,ivbd,ivbdEnd,ivbdStart,lvbd
    
    ! Check if handle are not empty
    if ((h_IboundaryCpIdx .eq. ST_NOHANDLE) .or.&
        (h_IverticesAtBoundary .eq. ST_NOHANDLE) .or.&
        (h_DvertexParameterValue .eq. ST_NOHANDLE)) then
      call output_line('Invalid handle!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setBoundary')
      call sys_halt()
    end if

    ! Check if boundary structure exists and remove it
    if (associated(rhadapt%rBoundary)) then
      do ibct = 1, size(rhadapt%rBoundary,1)
        call map_release(rhadapt%rBoundary(ibct))
      end do
      deallocate(rhadapt%rBoundary)
    end if

    ! Set pointers
    call storage_getbase_double(h_DvertexParameterValue,&
                                p_DvertexParameterValue, nvbd)
    call storage_getbase_int(h_IboundaryCpIdx,&
                             p_IboundaryCpIdx, nbct+1)
    call storage_getbase_int(h_IverticesAtBoundary,&
                             p_IverticesAtBoundary, nvbd)
    
    ! Allocate array of maps
    allocate(rhadapt%rBoundary(nbct))
    
    ! Initialization
    ioff = 0

    do ibct = 1, nbct

      ! Create a separate map for each boundary component
      lvbd = p_IboundaryCpIdx(ibct+1)-p_IboundaryCpIdx(ibct)
      call map_create(rhadapt%rBoundary(ibct), lvbd, 3)

      ! Set subdimensions
      ivbdStart = ioff+1
      ivbdEnd   = ioff+lvbd
      ioff      = ioff+lvbd

      ! Insert first element into map
      Ddata(1) =      p_DvertexParameterValue(ivbdStart)
      Ddata(2) = real(p_IverticesAtBoundary(ivbdEnd),DP)
      Ddata(3) = real(p_IverticesAtBoundary(ivbdStart+1),DP)

      rmapIter = map_insert(rhadapt%rBoundary(ibct),&
                            p_IverticesAtBoundary(ivbdStart), Ddata)

      ! Insert second to last but one element into map
      do ivbd = ivbdStart+1, ivbdEnd-1
        Ddata(1) =      p_DvertexParameterValue(ivbd)
        Ddata(2) = real(p_IverticesAtBoundary(ivbd-1),DP)
        Ddata(3) = real(p_IverticesAtBoundary(ivbd+1),DP)
        
        rmapIter = map_insert(rhadapt%rBoundary(ibct),&
                              p_IverticesAtBoundary(ivbd), Ddata)
      end do

      ! Insert last element into map
      Ddata(1) =      p_DvertexParameterValue(ivbdEnd)
      Ddata(2) = real(p_IverticesAtBoundary(ivbdEnd-1),DP)
      Ddata(3) = real(p_IverticesAtBoundary(ivbdStart),DP)

      rmapIter = map_insert(rhadapt%rBoundary(ibct),&
                            p_IverticesAtBoundary(ivbdEnd), Ddata)
    end do
    
    ! Set dimensions
    rhadapt%NBCT  = nbct
    rhadapt%NVBD  = nvbd
    rhadapt%NVBD0 = nvbd

    ! Set specifier for boundary
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_HAS_BOUNDARY)

  end subroutine hadapt_setBoundary
  
  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_getBoundary(rhadapt, h_IboundaryCpIdx, h_IverticesAtBoundary,&
                                h_DvertexParameterValue, nvbd, nbct)

!<description>
    ! This subroutine extracts the boundary data from the adaptivity structure
    ! and generates the the arrays p_IboundaryCpIdx, p_IverticesAtBoundary and
    ! p_DvertexParameterValue. If the optional parameter nvbd is given
    ! then the number of boundary vertices is assigned.
!</description>

!<input>
    ! Adaptivity structure
    type(t_hadapt), intent(in) :: rhadapt
!</input>

!<inputoutput>
    ! Handle to p_IboundaryCpIdx
    integer, intent(inout) :: h_IboundaryCpIdx
    
    ! Handle to p_IverticesAtBoundary
    integer, intent(inout) :: h_IverticesAtBoundary

    ! Handle to p_DvertexParameterValue
    integer, intent(inout) :: h_DvertexParameterValue

    ! OPTIONAL: number of vertices at the boundary
    integer, intent(out), optional :: nvbd

    ! OPTIONAL: number of boundary components
    integer, intent(out), optional :: nbct
!</inputoutput>
!</subroutine>

    ! local variables
    type(it_mapInt_Dble) :: rmapIter
    real(DP), dimension(:), pointer :: p_DvertexParameterValue, p_Ddata
    integer, dimension(:), pointer  :: p_IboundaryCpIdx
    integer, dimension(:), pointer  :: p_IverticesAtBoundary
    integer :: ibct,ioff,isize,ivbd,ivbdEnd,ivbdStart,lvbd

    ! Check if boundary data exists
    if (iand(rhadapt%iSpec, HADAPT_HAS_BOUNDARY) .ne.&
                            HADAPT_HAS_BOUNDARY) then
      call output_line('Structure does not exist!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getBoundary')
      call sys_halt()
    end if

    ! Due to the fact, that the boundary data may be stored in
    ! multiple maps (one for each component) the static arrays are
    ! pre-allocated and filled step-by-step from the dynamic data in
    ! the boundary map(s).

    ! Check if the arrays p_IboundaryCpIdx, p_IverticesAtBoundary and
    ! p_DvertexParameterValue have correct dimension
    if (h_IboundaryCpIdx .eq. ST_NOHANDLE) then
      call storage_new('hadapt_getBoundary', 'p_IboundaryCpIdx',&
                       rhadapt%NBCT+1, ST_INT, h_IboundaryCpIdx,&
                       ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_IboundaryCpIdx, isize)
      if (isize < rhadapt%NBCT+1) then
        call storage_realloc('hadapt_getBoundary', rhadapt%NBCT+1,&
                             h_IboundaryCpIdx, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if

    if (h_IverticesAtBoundary .eq. ST_NOHANDLE) then
      call storage_new('hadapt_getBoundary', 'p_IverticesAtBoundary',&
                       rhadapt%NVBD, ST_INT, h_IverticesAtBoundary,&
                       ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_IverticesAtBoundary, isize)
      if (isize < rhadapt%NVBD) then
        call storage_realloc('hadapt_getBoundary', rhadapt%NVBD,&
                             h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if

    if (h_DvertexParameterValue .eq. ST_NOHANDLE) then
      call storage_new('hadapt_getBoundary', 'p_DvertexParameterValue',&
                       rhadapt%NVBD, ST_DOUBLE, h_DvertexParameterValue,&
                       ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_DvertexParameterValue, isize)
      if (isize < rhadapt%NVBD) then
        call storage_realloc('hadapt_getBoundary', rhadapt%NVBD,&
                             h_DvertexParameterValue, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if

    ! Set pointers
    call storage_getbase_int(h_IboundaryCpIdx,&
                             p_IboundaryCpIdx, rhadapt%NBCT+1)
    call storage_getbase_int(h_IverticesAtBoundary,&
                             p_IverticesAtBoundary, rhadapt%NVBD)
    call storage_getbase_double(h_DvertexParameterValue,&
                                p_DvertexParameterValue, rhadapt%NVBD)

    ! Initialization
    p_IboundaryCpIdx(1) = 1
    ioff = 0
    ivbd = 0

    do ibct = 1, rhadapt%NBCT
      
      ! Set subdimensions
      lvbd      = map_size(rhadapt%rBoundary(ibct))
      ivbdStart = ioff+1
      ivbdEnd   = ioff+lvbd
      ioff      = ioff+lvbd

      ! Set index for next boundary component
      p_IboundaryCpIdx(ibct+1) = ioff+1

      ! Get static boundary data from map
      rmapIter = map_begin(rhadapt%rBoundary(ibct))
      
      do while(.not.map_isNull(rmapIter))
        
        ! Increase counter
        ivbd = ivbd+1

        ! Get key value
        p_IverticesAtBoundary(ivbd) = map_get(rhadapt%rBoundary(ibct), rmapIter)

        ! Get auxiliary data
        call map_getbase_data(rhadapt%rBoundary(ibct), rmapIter, p_Ddata)
        p_DvertexParameterValue(ivbd) = p_Ddata(1)

        ! Increase iterators
        call map_next(rmapIter)        
      end do

      ! Consistency check
      if (ivbd .ne. ivbdEnd) then
        call output_line('An internal error occured!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getBoundary')
        call sys_halt()
      end if

      ! Sort array w.r.t. increasing parameter values
      call sort_dp(p_DvertexParameterValue(ivbdStart:ivbdEnd),&
                   SORT_QUICK, p_IverticesAtBoundary(ivbdStart:ivbdEnd))
    end do

    ! Set dimension
    if(present(nvbd)) nvbd = rhadapt%NVBD
    if(present(nbct)) nbct = rhadapt%NBCT
  end subroutine hadapt_getBoundary

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_setNodalProperty(rhadapt, h_InodalProperty)

!<description>
    ! This subroutine sets the nodal property list to the adaptivity structure
!</description>

!<input>
    ! Handle to p_InodalProperty
    integer, intent(in) :: h_InodalProperty
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    ! Check if handle is not empty
    if (h_InodalProperty .eq. ST_NOHANDLE) then
      call output_line('Invalid handle!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setNodalProperty')
      call sys_halt()
    end if

    rhadapt%h_InodalProperty = h_InodalProperty
    call storage_getbase_int(rhadapt%h_InodalProperty,&
                             rhadapt%p_InodalProperty)
    
    ! Set specifier for InodalProperty
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_HAS_NODALPROP)
  end subroutine hadapt_setNodalProperty

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_getNodalProperty(rhadapt, h_InodalProperty)

!<description>
    ! This subroutine assignes the "nodal property" array from the
    ! adaptivity structure to the given handle
!</description>

!<input>
    ! Adaptivity structure
    type(t_hadapt), intent(in) :: rhadapt
!</input>

!<inputoutput>
    ! Hande to p_InodalProperty
    integer, intent(inout) :: h_InodalProperty
!</inputoutput>
!</subroutine>

    ! Check if "nodal property" list exits
    if (iand(rhadapt%iSpec, HADAPT_HAS_NODALPROP) .ne.&
                            HADAPT_HAS_NODALPROP) then
      call output_line('Structure does not exist!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getNodalProperty')
      call sys_halt()
    end if

    ! Check if handle needs to be freed first
    if ((h_InodalProperty .ne. ST_NOHANDLE) .and.&
        (h_InodalProperty .ne. rhadapt%h_InodalProperty))&
        call storage_free(h_InodalProperty)

    ! Assign handle
    h_InodalProperty = rhadapt%h_InodalProperty
  end subroutine hadapt_getNodalProperty

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_genElementsAtVertex(rhadapt)

!<description>
    ! This subroutine generates a list of arrays for the "elements-at-vertex"
    ! structure. This is an internal structure of the adaptivity structure
    ! which is used, e.g., in the removel of vertices.
!</description>

!<inputoutput>
    ! Adaptivity structur
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: iel,ive

    ! Check if "vertices-at-element" list exists
    if (iand(rhadapt%iSpec, HADAPT_HAS_VERTATELEM) .ne.&
                            HADAPT_HAS_VERTATELEM) then
      call output_line('Structure does not exist!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_genElementsAtVertex')
      call sys_halt()
    end if

    ! Create arraylist
    call alst_create(rhadapt%relementsAtVertex, 2*rhadapt%NVT, 8*rhadapt%NEL)

    ! Fill arraylist
    do iel = 1, rhadapt%NEL
      do ive = 1, hadapt_getNVE(rhadapt,iel)
        call alst_push_back(rhadapt%relementsAtVertex,&
                            rhadapt%p_IverticesAtElement(ive,iel), iel)
      end do
    end do
    
    ! Set specifier for relementsAtVertex
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_HAS_ELEMATVERTEX)
  end subroutine hadapt_genElementsAtVertex

end module hadaptaux
