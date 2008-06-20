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

MODULE hadaptaux
  
  USE arraylist
  USE binarytree
  USE fsystem
  USE octree
  USE quadtree
  USE storage
  USE triangulation

  IMPLICIT NONE

  PUBLIC

!<constantblock description="Global flags for grid refinement/coarsening">

  ! No refinement and coarsening
  INTEGER, PARAMETER, PUBLIC :: HADAPT_NOADAPTATION         = 0

  ! Red-Green refinement and coarsening strategy (R. Banks et al.)
  INTEGER, PARAMETER, PUBLIC :: HADAPT_REDGREEN             = 1

  ! Longest edge bisection strategy (M. Rivara)
  INTEGER, PARAMETER, PUBLIC :: HADAPT_LONGESTEDGE          = 2

!</constantblock>

!<constantblock description="Bitfield identifiers for state of adaptation">

  ! Adaptation is undefined
  INTEGER, PARAMETER, PUBLIC :: HADAPT_UNDEFINED     = 2**0

  ! Parameters of adaptivity structure are initialized
  INTEGER, PARAMETER :: HADAPT_HAS_PARAMETERS        = 2**1

  ! Quadtree/octree for vertex coordinates is generated
  INTEGER, PARAMETER :: HADAPT_HAS_COORDS            = 2**2

  ! Array for IverticesAtElement is generated
  INTEGER, PARAMETER :: HADAPT_HAS_VERTATELEM        = 2**3

  ! Array for IneighboursAtElement is generated
  INTEGER, PARAMETER :: HADAPT_HAS_NEIGHATELEM       = 2**4

  ! Boundary data is generated
  INTEGER, PARAMETER :: HADAPT_HAS_BOUNDARY          = 2**5

  ! Nodal property is generated
  INTEGER, PARAMETER :: HADAPT_HAS_NODALPROP         = 2**6

  ! Number of elements for predefined type
  INTEGER, PARAMETER :: HADAPT_HAS_NELOFTYPE         = 2**7

  ! Array for IelementsAtVertex is generated
  INTEGER, PARAMETER :: HADAPT_HAS_ELEMATVERTEX      = 2**8

  ! Dynamic data structures are all generated
  INTEGER, PARAMETER :: HADAPT_HAS_DYNAMICDATA       = HADAPT_HAS_PARAMETERS+&
                                                       HADAPT_HAS_COORDS+&
                                                       HADAPT_HAS_VERTATELEM+&
                                                       HADAPT_HAS_NEIGHATELEM+&
                                                       HADAPT_HAS_BOUNDARY+&
                                                       HADAPT_HAS_NODALPROP+&
                                                       HADAPT_HAS_NELOFTYPE+&
                                                       HADAPT_HAS_ELEMATVERTEX

  ! Cells are marked for refinement
  INTEGER, PARAMETER :: HADAPT_MARKEDREFINE          = 2**9

  ! Cells are marked for coarsening
  INTEGER, PARAMETER :: HADAPT_MARKEDCOARSEN         = 2**10

  ! Cells are marked
  INTEGER, PARAMETER :: HADAPT_MARKED                = HADAPT_MARKEDREFINE+&
                                                       HADAPT_MARKEDCOARSEN
  
  ! Grid has been refined
  INTEGER, PARAMETER :: HADAPT_REFINED               = 2**11
  
  ! Grid has been coarsened
  INTEGER, PARAMETER :: HADAPT_COARSENED             = 2**12

!</constantblock>

!<constantblock description="Constants for grid adaptation">
  
  ! Array position of the boundary
  INTEGER, PARAMETER :: BdrValue = 1
  
  ! Array position of the previous boundary vertex
  INTEGER, PARAMETER :: BdrPrev  = 1

  ! Array position of the next boundary vertex
  INTEGER, PARAMETER :: BdrNext  = 2

!</constantblock>

!<constantblock description="Duplication flags. Specifies which information is
!                            shared between adaptivity structures">

  INTEGER(I32), PARAMETER :: HADAPT_SHARE_IMARKER                 = 2** 0
  INTEGER(I32), PARAMETER :: HADAPT_SHARE_IVERTEXAGE              = 2** 1
  INTEGER(I32), PARAMETER :: HADAPT_SHARE_INODALPROPERTY          = 2** 2
  INTEGER(I32), PARAMETER :: HADAPT_SHARE_IVERTICESATELEMENT      = 2** 3
  INTEGER(I32), PARAMETER :: HADAPT_SHARE_INEIGHATELEMENT         = 2** 4
  INTEGER(I32), PARAMETER :: HADAPT_SHARE_IMIDNEIGHATELEMENT      = 2** 5
  INTEGER(I32), PARAMETER :: HADAPT_SHARE_RVERTEXCOORDINATES      = 2** 6
  INTEGER(I32), PARAMETER :: HADAPT_SHARE_RBOUNDARY               = 2** 7
  INTEGER(I32), PARAMETER :: HADAPT_SHARE_RELEMENTSATVERTEX       = 2** 8

!</constantblock>

!</constants>

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<types>

  !<typeblock>
  
  ! This type contains all data structures to handle 
  ! adaptive grid refinement and grid coarsening.
  TYPE :: t_hadapt
    ! Format Tag: Specifies the state of adaptation
    INTEGER :: iSpec                                 = HADAPT_UNDEFINED

    ! Duplication flag. Bitfield that indicates which information is
    ! shared with another adaptivity structure.
    ! When a bit is set to 1, the corresponding array is
    ! maintained by another adaptivity structure and must
    ! not be deleted by hadapt_releaseAdaptation. 
    ! When the bit is 0, the array is a real copy of another array 
    ! and must be deleted in hadapt_releaseAdaptation.
    INTEGER(I32) :: iduplicationFlag                 = 0

    ! Tag: Specified the strategy for grid refinement and coarsening
    INTEGER :: iadaptationStrategy                   = HADAPT_NOADAPTATION

    ! Maximum number of subdivisions from the original mesh
    INTEGER :: NSUBDIVIDEMAX                         = 0

    ! Total number of grid refinement steps
    INTEGER :: nRefinementSteps                      = 0

    ! Total number of grid coarsening steps
    INTEGER :: nCoarseningSteps                      = 0

    ! Total number of grid smoothing steps
    INTEGER :: nSmoothingSteps                       = 0

    ! Tolerance for refinement
    REAL(DP) :: drefinementTolerance                 = 0

    ! Tolerance for coarsening
    REAL(DP) :: dcoarseningTolerance                 = 0

    ! Dimension of the triangulation
    INTEGER :: ndim                                  = 0

    ! Total number of vertices (initially)
    INTEGER(PREC_VERTEXIDX) :: NVT0                  = 0
    
    ! Total number of vertices
    INTEGER(PREC_VERTEXIDX) :: NVT                   = 0

    ! Increment of vertices
    INTEGER(PREC_VERTEXIDX) :: increaseNVT           = 0

    ! Total number of boundary vertives (initially)
    INTEGER :: NVBD0                                 = 0

    ! Total number of boundary vertives
    INTEGER :: NVBD                                  = 0
    
    ! Total number of boundary components (should not change)
    INTEGER :: NBCT                                  = 0
    
    ! Total number of elements (initially)
    INTEGER(PREC_ELEMENTIDX) :: NEL0                 = 0
    
    ! Total number of elements
    INTEGER(PREC_ELEMENTIDX) :: NEL                  = 0

    ! Maximum number of elements (before reallocation)
    INTEGER(PREC_ELEMENTIDX) :: NELMAX               = 0
    
    ! Total number of green elements (required internally)
    INTEGER(PREC_ELEMENTIDX) :: nGreenElements       = 0
    
    ! Nuber of elements with a defined number of vertices per element.
    ! InelOfType(TRIA_NVETRI2D)  = number of triangles in the mesh (2D).
    ! InelOfType(TRIA_NVEQUAD2D) = number of quadrilaterals in the mesh (2D).
    INTEGER(PREC_ELEMENTIDX), DIMENSION(TRIA_MAXNVE) :: InelOfType = 0

    ! Same as InelOfType but this array stores the number of elements
    ! which are initially present in the mesh
    INTEGER(PREC_ELEMENTIDX), DIMENSION(TRIA_MAXNVE) :: InelOfType0 = 0
    
    ! Element marker array.
    ! Handle to
    !       p_Imarker = array [1..NEL] of integer.
    ! For each element there is one identifier [0,..,31] in the marker.
    ! BIT0: Is 0 for triangle and 1 for quadrilateral
    ! BIT1: Is 0 if first edge is not subdivided, 1 otherwise.
    ! BIT2: Is 0 if second edge is not subdivided, 1 otherwise.
    ! BIT3: Is 0 if third edge is not subdivided, 1 otherwise.
    ! BIT4: Is 0 if fourth edge is not subdivided, 1 otherwise.
    INTEGER :: h_Imarker = ST_NOHANDLE
    
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
    INTEGER :: h_IvertexAge = ST_NOHANDLE

    ! Pointer to h_IvertexAge.
    ! This array is introduced to increase performance and must
    ! not be modified by the user
    INTEGER, DIMENSION(:), POINTER :: p_IvertexAge => NULL()

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
    INTEGER :: h_InodalProperty = ST_NOHANDLE

    ! Pointer to h_InodalProperty.
    ! This array is introduced to increase performance and must
    ! not be modified by the user
    INTEGER, DIMENSION(:), POINTER :: p_InodalProperty => NULL()
    
    ! Vertices adjacent to an element.
    ! Handle to 
    !       p_IverticesAtElement = array [1..TRIA_MAXNVE2D,1..NEL] of integer.
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    INTEGER :: h_IverticesAtElement = ST_NOHANDLE

    ! Pointer to h_IverticesAtElement.
    ! This array is introduced to increase performance and must
    ! not be modified by the user
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement => NULL ()

    ! Neighbour elements adjacent to an element.
    ! Handle to
    !       p_IneighboursAtElement = array [1..TRIA_MAXNME2D,1..NEL] of integer
    ! For each element, the numbers of adjacent elements in mathematically
    ! positive sense, metting the element in an edge.
    INTEGER :: h_IneighboursAtElement = ST_NOHANDLE

    ! Pointer to h_IneighboursAtElement.
    ! This array is introduced to increase performance and must
    ! not be modified by the user
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElement => NULL ()

    ! Midneighbour elements adjacent to an element.
    ! Handle to
    !       p_ImidneighboursAtElement = array [1..TRIA_MAXNME2D,1..NEL] of integer
    ! H-adaptivity is performed element-by-element. If one element is
    ! refined, then its neighbors need to be informed that there
    ! are possibly two elements adjacent along one edge. This is a
    ! nonconforming state. When time comes to process the other
    ! element (with hanging node), the element knows which elements
    ! are adjacent along the first and the second half of the edge.
    INTEGER :: h_ImidneighboursAtElement = ST_NOHANDLE

    ! Pointer to h_ImidneighboursAtElement.
    ! This array is introduced to increase performance and must
    ! not be modified by the user
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_ImidneighboursAtElement => NULL ()
    
    ! Quadtree storing the nodal coordinates in 2D
    TYPE(t_quadtree) :: rVertexCoordinates2D

    ! Octree storing the nodal coordinates in 2D
    TYPE(t_octree) :: rVertexCoordinates3D
    
    ! Array of binary search trees storing the boundary data
    ! p_IboundaryCpIdx and p_IverticesAtBoundary
    TYPE(t_btree), DIMENSION(:), POINTER :: rBoundary => NULL()

    ! Arraylist for elements-meeting-at-vertex structure
    TYPE(t_arraylist) :: rElementsAtVertex
  END TYPE t_hadapt

  !</typeblock> 

!</types>

CONTAINS

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<function>

  PURE FUNCTION hadapt_getNVE(rhadapt, iel) RESULT(nve)

!<description>
    ! This function returns the number of vertices present in the given element
!</description>

!<input>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(IN)           :: rhadapt

    ! Number of the element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel
!</input>

!<result>
    ! number of vertices per element
    INTEGER :: nve
!</result>
!</function

    ! Which spatial dimension do we have?
    SELECT CASE(rhadapt%ndim)

    CASE (NDIM2D)
      
      ! Do we have quadrilaterals in the triangulation?
      IF (rhadapt%InelOfType(TRIA_NVEQUAD2D) .EQ. 0) THEN
        
        ! There are no quadrilaterals in the current triangulation.
        ! Hence, return TRIA_NVETRI2D by default.
        nve = TRIA_NVETRI2D
        
      ELSE
        
        ! There are quadrilaterals and possible also triangles in
        ! the current triangulatin. If the last entry of the
        ! vertices-at-element list is nonzero then TRIA_NVEQUAD2D vertices
        ! are present in the current element. Otherwise return TRIA_NVETRI2D.
        IF (rhadapt%p_IverticesAtElement(TRIA_NVEQUAD2D, iel) .EQ. 0) THEN
          nve = TRIA_NVETRI2D
        ELSE
          nve = TRIA_NVEQUAD2D
        END IF
        
      END IF

    CASE DEFAULT
      nve = 0
    END SELECT
  END FUNCTION hadapt_getNVE

END MODULE hadaptaux
