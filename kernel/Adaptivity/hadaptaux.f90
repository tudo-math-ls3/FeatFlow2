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
!# which are required for performing h-adaptivity in. Unlike other modules, 
!# all subroutines are declared PUBLIC since they are used by module HADAPTIVITY.
!#
!# The following routines are available:
!# 
!# 1.) hadapt_getNVE
!#     -> return the number of vertices per element
!#
!# </purpose>
!##############################################################################

module hadaptaux
  
  use arraylist
  use binarytree
  use fsystem
  use octree
  use quadtree
  use storage
  use triangulation

  implicit none

  public

!<constantblock description="Global flags for grid refinement/coarsening">

  ! No refinement and coarsening
  integer, parameter, public :: HADAPT_NOADAPTATION = 0

  ! Red-Green refinement and coarsening strategy (R. Banks et al.)
  integer, parameter, public :: HADAPT_REDGREEN     = 1

  ! Longest edge bisection strategy (M. Rivara)
  integer, parameter, public :: HADAPT_LONGESTEDGE  = 2

!</constantblock>

!<constantblock description="Bitfield identifiers for state of adaptation">

  ! Adaptation is undefined
  integer, parameter, public :: HADAPT_UNDEFINED = 2**0

  ! Parameters of adaptivity structure are initialized
  integer, parameter :: HADAPT_HAS_PARAMETERS   = 2**1

  ! Quadtree/octree for vertex coordinates is generated
  integer, parameter :: HADAPT_HAS_COORDS       = 2**2

  ! Array for IverticesAtElement is generated
  integer, parameter :: HADAPT_HAS_VERTATELEM   = 2**3

  ! Array for IneighboursAtElement is generated
  integer, parameter :: HADAPT_HAS_NEIGHATELEM  = 2**4

  ! Boundary data is generated
  integer, parameter :: HADAPT_HAS_BOUNDARY     = 2**5

  ! Nodal property is generated
  integer, parameter :: HADAPT_HAS_NODALPROP    = 2**6

  ! Number of elements for predefined type
  integer, parameter :: HADAPT_HAS_NELOFTYPE    = 2**7

  ! Array for IelementsAtVertex is generated
  integer, parameter :: HADAPT_HAS_ELEMATVERTEX = 2**8

  ! Dynamic data structures are all generated
  integer, parameter :: HADAPT_HAS_DYNAMICDATA  = HADAPT_HAS_PARAMETERS+&
                                                  HADAPT_HAS_COORDS+&
                                                  HADAPT_HAS_VERTATELEM+&
                                                  HADAPT_HAS_NEIGHATELEM+&
                                                  HADAPT_HAS_BOUNDARY+&
                                                  HADAPT_HAS_NODALPROP+&
                                                  HADAPT_HAS_NELOFTYPE+&
                                                  HADAPT_HAS_ELEMATVERTEX

  ! Cells are marked for refinement
  integer, parameter :: HADAPT_MARKEDREFINE     = 2**9

  ! Cells are marked for coarsening
  integer, parameter :: HADAPT_MARKEDCOARSEN    = 2**10

  ! Cells are marked
  integer, parameter :: HADAPT_MARKED           = HADAPT_MARKEDREFINE+&
                                                  HADAPT_MARKEDCOARSEN
  
  ! Grid has been refined
  integer, parameter :: HADAPT_REFINED          = 2**11
  
  ! Grid has been coarsened
  integer, parameter :: HADAPT_COARSENED        = 2**12

!</constantblock>

!<constantblock description="Constants for grid adaptation">
  
  ! Array position of the boundary
  integer, parameter :: BdrValue = 1
  
  ! Array position of the previous boundary vertex
  integer, parameter :: BdrPrev  = 1

  ! Array position of the next boundary vertex
  integer, parameter :: BdrNext  = 2

!</constantblock>

!<constantblock description="Duplication flags. Specifies which information is
!                            shared between adaptivity structures">

  integer(I32), parameter :: HADAPT_SHARE_IMARKER            = 2** 0
  integer(I32), parameter :: HADAPT_SHARE_IVERTEXAGE         = 2** 1
  integer(I32), parameter :: HADAPT_SHARE_INODALPROPERTY     = 2** 2
  integer(I32), parameter :: HADAPT_SHARE_IVERTICESATELEMENT = 2** 3
  integer(I32), parameter :: HADAPT_SHARE_INEIGHATELEMENT    = 2** 4
  integer(I32), parameter :: HADAPT_SHARE_IMIDNEIGHATELEMENT = 2** 5
  integer(I32), parameter :: HADAPT_SHARE_RVERTEXCOORDINATES = 2** 6
  integer(I32), parameter :: HADAPT_SHARE_RBOUNDARY          = 2** 7
  integer(I32), parameter :: HADAPT_SHARE_RELEMENTSATVERTEX  = 2** 8

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
    integer :: iSpec = HADAPT_UNDEFINED

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
    integer :: NSUBDIVIDEMAX = 0

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
    ! characterizes the type of the node (=vertex/edge):
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
    
    ! Binary tree storing the nodal coordinates in 1D
    type(t_btree) :: rVertexCoordinates1D
    
    ! Quadtree storing the nodal coordinates in 2D
    type(t_quadtree) :: rVertexCoordinates2D

    ! Octree storing the nodal coordinates in 2D
    type(t_octree) :: rVertexCoordinates3D
    
    ! Array of binary search trees storing the boundary data
    ! p_IboundaryCpIdx and p_IverticesAtBoundary
    type(t_btree), dimension(:), pointer :: rBoundary => null()

    ! Arraylist for elements-meeting-at-vertex structure
    type(t_arraylist) :: rElementsAtVertex
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
    type(t_hadapt), intent(IN) :: rhadapt

    ! Number of the element
    integer, intent(IN) :: iel
!</input>

!<result>
    ! number of vertices per element
    integer :: nve
!</result>
!</function

    ! Which spatial dimension do we have?
    select case(rhadapt%ndim)
    case (NDIM1D)
      nve = 2

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

end module hadaptaux
