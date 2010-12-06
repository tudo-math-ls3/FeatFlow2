!##############################################################################
!# ****************************************************************************
!# <name> transformation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for the transformation between a reference
!# element and a real element.
!#
!# 1.) trafo_getCoords
!#     -> Determine the coordinates of those points on an element that
!#        specify a transformation.
!#
!# 2.) trafo_getCoords_sim
!#     -> Determine the coordinates of those points on a list of elements that
!#        specify the transformations.
!#
!# 3.) trafo_igetDimension
!#     -> Determine the underlying dimension of a transformation.
!#
!# 5.) trafo_igetNVE
!#     -> Determine the number of vertices on an element that are needed to
!#        specify the transformation from the reference to the real element
!#        (=number of degrees of freedom of the transformation formula)
!#
!# 6.) trafo_calctrafo / trafo_calctrafoabs
!#     -> Calculates the transformation of one point (from reference to the 
!#        real element).
!#        Supports triangular and quadrilateral mapping.
!#        trafo_calctrafoabs calculates the absolute value of the Jacobian 
!#        determinant.
!#
!# 7.) trafo_calctrafo_mult / trafo_calctrafoabs_mult
!#     -> Calculates the transformation for multiple points on one element.
!#        Supports triangular and quadrilateral mapping.
!#        trafo_calctrafoabs_mult calculates the absolute value of the Jacobian
!#        determinant.
!#
!# 8.) trafo_calctrafo_sim / trafo_calctrafoabs_sim
!#     -> Calculate the transformation for multiple points on multiple
!#        elements. Supports triangular and quadrilateral mapping.
!#        trafo_calctrafoabs_sim calculates the absolute value of the Jacobian
!#        determinant.
!#
!# 10.) trafo_prepJac
!#     -> calculates auxiliary Jacobian factors for the transformation
!#        from a reference quadrilateral to a real quadrilateral
!#
!# 12.) trafo_calcJac
!#     -> calculates the Jacobian matrix + Jacobian determinant of the mapping
!#        from  the reference to a real quadrilateral element
!#
!# 13.) trafo_calcRealCoords
!#      -> maps a point from the reference element to the real element
!#
!# 14.) trafo_calcTrafo_quad2d
!#      -> calculates the Jacobian matrix + Jacobian determinant of the mapping
!#         from  the reference to a real quadrilateral element
!#       -> maps a point from the reference element to the real element
!#      (so performing the same task as trafo_calcJac and trafo_calcRealCoords
!#       in one routine)
!#
!# 15.) trafo_mapCubPts1Dto2D
!#      -> Maps a set of 1D cubature point coordinates on the reference 
!#         interval to an edge of a reference element.
!#
!# 16.) trafo_mapCubPts1Dto2DRefQuad
!#      -> Maps a set of 1D cubature point coordinates on the reference 
!#         interval to an edge of the 2D reference quadrilateral.
!#
!# 17.) trafo_mapCubPts1Dto2DTriBary
!#      -> Maps a set of 1D cubature point coordinates on the reference 
!#         interval to an edge of the 2D reference triangle; barycentric
!#         coordinates.
!# 
!# 18.) trafo_calcBackTrafo_quad2d
!#      -> Find the coordinates on the reference element 
!#         for a given point in real world coordinates
!# 
!# 19.) trafo_calcRefCoords
!#      -> For a given point in world coordinates and an element, this 
!#         routine calculates the coordinates of the point on the reference
!#         element.
!#         Supports triangular and quadrilateral mapping.
!#
!# 20.) trafo_mapCubPts2Dto3DRefHexa
!#      -> Maps a set of 2D cubature point coordinates on the reference 
!#         quadrilateral to a face of the 3D reference hexahedron.
!#
!# 21.) trafo_mapCubPtsRef2LvlEdge1D
!#      -> Maps a set of 1D cubature points coordinate on the reference
!#         fine mesh edge to a coarse mesh edge according to the 2-level-
!#         ordering refinement algorithm.
!#
!# 22.) trafo_mapCubPtsRef2LvlTri2D
!#      -> Maps a set of 2D cubature points coordinate on the reference
!#         fine mesh triangle to a coarse mesh triangle according
!#         to the 2-level-ordering refinement algorithm.
!#
!# 23.) trafo_mapCubPtsRef2LvlQuad2D
!#      -> Maps a set of 2D cubature points coordinate on the reference
!#         fine mesh quadrilateral to a coarse mesh quadrilateral according
!#         to the 2-level-ordering refinement algorithm.
!#
!# 24.) trafo_mapCubPtsRef2LvlHexa3D
!#      -> Maps a set of 3D cubature points coordinate on the reference
!#         fine mesh hexahedron to a coarse mesh hexahedron according to
!#         the 2-level-ordering refinement algorithm.
!#
!# 25.) trafo_igetReferenceDimension
!#      -> Calculates the dimension of the coordinate system on the
!#         reference element
!#
!# 26.) trafo_calcBackTrafo_hexa
!#      -> Find the coordinates on the reference hexahedron
!#         for a given point in real world coordinates
!#
!#  FAQ - Some explainations  \\
!# -------------------------- \\
!# 1.) How is the ID code ctrafoType of the transformation defined?
!#
!#  The ID code of a transformation is a 32 bit integer field. The different
!#  bits in the transformation encode the type of transformation, the dimension
!#  and additional information, as far as necessary. The exact layout is
!#  defined as follows:
!#
!#  <verb>
!#     Bit | 31 ... 24 23 ... 16 | 15 ... 10 | 9   8 | 7 ............ 0|
!#     -----------------------------------------------------------------
!#         |         ****        | unused    |dimens |     trafo-ID    |
!#
!#   Bits 0..7   specify the type of transformation (triangular, quad, 
!#               linear, quadratic, ...). The different ID`s for this
!#               field are expressed in the TRAFO_ID_xxxx constants.
!#   Bits 8+ 9   encode the dimension of the reference element. =1: 1D, =2: 2D, =3: 3D.
!#   Bits 10-15  unused
!#   Bits 16..31 encode special variants of the transformation. Normally,these
!#               bits are all =0 except for those types of transformation that
!#               require special, additional information, like isoparametric elements.
!#         trafo-ID = ...
!#               TRAFO_ID_QUADTRI: 
!#                    Bit 16/17: These three bits encode how many edges of the tri
!#                               element are to be handled with an isoparametric mapping
!#                               to the reference element. 
!#                               =0: isoparametric mapping with one edge not linear
!#                               =1: isoparametric mapping with two edges not linear
!#                               =2: isoparametric mapping with three edges not linear
!#                    Bit 18/19/20: These three bits encode which of the edges are
!#                               to be mapped nonlinear from the reference to the real
!#                               element.
!#                               Bit 18:=1 1st edge maps nonlinear
!#                               Bit 19:=1 2nd edge maps nonlinear
!#                               Bit 20:=1 3rd edge maps nonlinear
!#               TRAFO_ID_MQUADQUAD: 
!#                    Bit 16/17: These three bits encode how many edges of the quad
!#                               element are to be handled with an isoparametric mapping
!#                               to the reference element. 
!#                               =0: isoparametric mapping with one edge not linear
!#                               =1: isoparametric mapping with two edges not linear
!#                               =2: isoparametric mapping with three edges not linear
!#                               =3: isoparametric mapping with all four edges not linear
!#                    Bit 18/19/20/21: These four bits encode which of the edges are
!#                               to be mapped nonlinear from the reference to the real
!#                               element.
!#                               Bit 18:=1 1st edge maps nonlinear
!#                               Bit 19:=1 2nd edge maps nonlinear
!#                               Bit 20:=1 3rd edge maps nonlinear
!#                               Bit 21:=1 4th edge maps nonlinear
!#   </verb>
!#
!#   Note: As long as the isoparametric mapping is not implemented, bits 16-31
!#   are unused!
!# </purpose>
!##############################################################################

module transformation

  use fsystem
  use basicgeometry
  use storage
  use triangulation
  use genoutput
  use geometryaux

  implicit none
  
  private
  
!<constants>
!<constantblock>
  ! Minimum number of elements for OpenMP parallelisation: If the number of
  ! elements is below this value, then no parallelisation is performed.
#ifndef TRAFO_NELEMMIN_OMP
  integer, parameter, public :: TRAFO_NELEMMIN_OMP = 1000
#endif
  
!</constantblock>


!<constantblock description="Constants for size of auxiliary arrays.">

  ! Number of entries in the array with the auxiliary Jacobian factors in 2D
  integer, parameter, public :: TRAFO_NAUXJAC2D = 8

  ! Number of entries in the array with the auxiliary Jacobian factors in 3D
  integer, parameter, public :: TRAFO_NAUXJAC3D = 24

  ! Max. Number of entries in the array, dimension independent
  integer, parameter, public :: TRAFO_NAUXJACMAX = TRAFO_NAUXJAC3D

!</constantblock>


!<constantblock description="Coordinate system constants.">

  ! Maximum number of dimensions on reference elements.
  ! =4 for Tetrahedra with barycentric coordinates in 3D
  integer, parameter, public :: TRAFO_MAXDIMREFCOORD = NDIM3D+1

!</constantblock>


!<constantblock description="Id values for bit 0..7 of the transformation ID.">

  ! Unspecified transformation
  integer, parameter, public :: TRAFO_ID_UNKNOWN       = 0

  ! Linear transformation for simplex-type elements (1D-lines, 2D-triangles, 
  ! 3D-tetrahedra)
  integer, parameter, public :: TRAFO_ID_LINSIMPLEX    = 1

  ! Multilinear (Bilinear/Trilinear) transformation for cubic-shaped elements 
  ! (2D-quadrilaterals, 3D-hexahedra)
  integer, parameter, public :: TRAFO_ID_MLINCUBE     = 2
  
  ! Multilinear transformation for pyramid shaped elements in 3D
  integer, parameter, public :: TRAFO_ID_MLINPYRAMID  = 3
  
  ! Multilinear transformation for prismic shaped elements in 3D
  integer, parameter, public :: TRAFO_ID_MLINPRISM    = 4
  
  ! Quadratic transformation for simplex-type elements (1D-lines, 2D-triangles, 
  ! 3D-tetrahedra)
  integer, parameter, public :: TRAFO_ID_QUADSIMPLEX  = 11

  ! Multiquadratic (Biquadratic/Triquadratic) quadrilateral transformation 
  ! for cubic-shaped elements (2D-quadrilaterals, 3D-hexahedra)
  integer, parameter, public :: TRAFO_ID_MQUADCUBE    = 12
!</constantblock>


!<constantblock description="Dimension constants for bit 8+9 of the transformation ID.">
  ! Bitmasks for dimension
  integer(I32), parameter, public :: TRAFO_DIM_DIMENSION = 2**8 + 2**9

  ! 1D element
  integer(I32), parameter, public :: TRAFO_DIM_1D = ishft(NDIM1D,8)

  ! 2D element
  integer(I32), parameter, public :: TRAFO_DIM_2D = ishft(NDIM2D,8)
  
  ! 3D element
  integer(I32), parameter, public :: TRAFO_DIM_3D = ishft(NDIM3D,8)

  ! Bitmask for transformation ID, without additional information
  integer(I32), parameter, public :: TRAFO_DIM_IDMASK = 255 

  ! Bitmask for transformation ID including dimension, without additional information
  integer(I32), parameter, public :: TRAFO_DIM_IDDIMMASK = TRAFO_DIM_IDMASK + TRAFO_DIM_DIMENSION
!</constantblock>


!<constantblock description="id values for coordinate systems">
  ! Undefined coordinate system or no coordinate system
  integer(I32), parameter, public :: TRAFO_CS_UNDEFINED   = 0

  ! Parameter value [-1..1] on reference interval in 1D
  integer(I32), parameter, public :: TRAFO_CS_REF1D       = TRAFO_DIM_1D + 1

  ! Barycentric coordinates on triangle
  integer(I32), parameter, public :: TRAFO_CS_BARY2DTRI   = TRAFO_DIM_2D + 2

  ! 2D coordinates on reference triangle
  integer(I32), parameter, public :: TRAFO_CS_REF2DTRI    = TRAFO_DIM_2D + 3

  ! 2D coordinates on real triangle
  integer(I32), parameter, public :: TRAFO_CS_REAL2DTRI   = TRAFO_DIM_2D + 4
  
  ! 2D coordinates on reference quadrilateral
  integer(I32), parameter, public :: TRAFO_CS_REF2DQUAD   = TRAFO_DIM_2D + 5

  ! 2D coordinates on real quadrilateral
  integer(I32), parameter, public :: TRAFO_CS_REAL2DQUAD  = TRAFO_DIM_2D + 6
  
  ! Barycentric coordinates on tetrahedron
  integer(I32), parameter, public :: TRAFO_CS_BARY3DTETRA = TRAFO_DIM_3D + 10
  
  ! 3D coordinates on reference tetrahedron
  integer(I32), parameter, public :: TRAFO_CS_REF3DTETRA  = TRAFO_DIM_3D + 11
  
  ! 3D coodinates on real tetrahedron
  integer(I32), parameter, public :: TRAFO_CS_REAL3DTETRA = TRAFO_DIM_3D + 12
  
  ! 3D coordinates on reference hexahedron
  integer(I32), parameter, public :: TRAFO_CS_REF3DHEXA   = TRAFO_DIM_3D + 13
  
  ! 3D coordinates on real hexahedron
  integer(I32), parameter, public :: TRAFO_CS_REAL3DHEXA  = TRAFO_DIM_3D + 14
  
  ! 3D coordinates on reference pyramid
  integer(I32), parameter, public :: TRAFO_CS_REF3DPYRA   = TRAFO_DIM_3D + 15
  
  ! 3D coordinates on real pyramid
  integer(I32), parameter, public :: TRAFO_CS_REAL3DPYRA  = TRAFO_DIM_3D + 16
  
  ! 3D coordinates on reference prism
  integer(I32), parameter, public :: TRAFO_CS_REF3DPRISM  = TRAFO_DIM_3D + 17
  
  ! 3D coordinates on real prism
  integer(I32), parameter, public :: TRAFO_CS_REAL3DPRISM = TRAFO_DIM_3D + 18
!</constantblock>

!</constants>

  interface trafo_calcRealCoords
    module procedure trafo_calcRealCoords2D
    module procedure trafo_calcRealCoords_general
  end interface
  
  public :: trafo_getCoords
  public :: trafo_getCoords_sim
  public :: trafo_igetDimension
  public :: trafo_igetNVE
  public :: trafo_calctrafo, trafo_calctrafoabs
  public :: trafo_calctrafo_mult, trafo_calctrafoabs_mult
  public :: trafo_calctrafo_sim, trafo_calctrafoabs_sim
  public :: trafo_calcRealCoords
  public :: trafo_mapCubPts1Dto2D
  public :: trafo_mapCubPts1Dto2DRefQuad
  public :: trafo_mapCubPts1Dto2DTriBary
  public :: trafo_calcBackTrafo_quad2d
  public :: trafo_calcRefCoords
  public :: trafo_mapCubPts2Dto3DRefHexa
  public :: trafo_mapCubPtsRef2LvlEdge1D
  public :: trafo_mapCubPtsRef2LvlTri2D
  public :: trafo_mapCubPtsRef2LvlQuad2D
  public :: trafo_mapCubPtsRef2LvlHexa3D
  public :: trafo_prepJac_quad2D
  public :: trafo_prepJac_hexa3D
  public :: trafo_prepJac_pyra3D
  public :: trafo_prepJac_prism3D
  public :: trafo_calcJac_quad2D
  public :: trafo_calcJac_hexa3D
  public :: trafo_calcJac_pyra3D
  public :: trafo_calcJac_prism3D
  public :: trafo_calcTrafo_quad2d 
  public :: trafo_calcTrafo_hexa3d 
  public :: trafo_calcTrafo_pyra3d 
  public :: trafo_calcTrafo_prism3d
  public :: trafo_igetReferenceDimension
  public :: trafo_calcBackTrafo_hexa

contains

  ! ***************************************************************************

!<function>  

  elemental integer function trafo_igetDimension(ctrafoType)

!<description>
  ! This function returns the dimensional constant that specifies which
  ! dimension (1D, 2D,...) of the space of the transformation.
  ! For a mapping from the reference to the real element, this is the
  ! number of dimensions of world coordinates.
!</description>

!<input>    
  ! ID of transformation to calculate
  integer(I32), intent(in) :: ctrafoType
!</input>

!<result>
  ! A constant that specifies the dimension of the transformation. NDIM2 for 2D,
  ! NDIM3D for 3D,...
!</result>

!</function>

    ! The dimension is encoded in two bits in the element quantifier!
    trafo_igetDimension = ishft(iand(ctrafoType,TRAFO_DIM_DIMENSION),-8)

  end function

  ! ***************************************************************************

!<function>  

  elemental integer function trafo_igetReferenceDimension(ctrafoType)

!<description>
  ! For a given transformation, this routine returns the dimension of
  ! the coordinate system on the reference element.
  ! E.g.: 3 for barycentric coordinates in 2D.
!</description>

!<input>    
  ! ID of transformation to calculate
  integer(I32), intent(in) :: ctrafoType
!</input>

!<result>
  ! Dimension of the coordinate system on the reference element.
!</result>

!</function>

    trafo_igetReferenceDimension = 0

    ! What type of transformation do we have? First decide on the dimension,
    ! then on the actual ID.
    select case (trafo_igetDimension(ctrafoType))

    case (NDIM1D)
      ! 1D elements. Lines. Check the actual transformation
      ! ID how to transform.
    
      select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
      
      case (TRAFO_ID_MLINCUBE)
        ! 1D linear line transformation. 
        trafo_igetReferenceDimension = 1
        
      end select
    
    case (NDIM2D)
      ! 2D elements. Triangles, Quadrilaterals. Check the actual transformation
      ! ID how to transform.
    
      select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
      
      case (TRAFO_ID_LINSIMPLEX)
        ! 2D simplex -> linear triangular transformation. 
        ! 3 DOF`s in the transformation (given by the corners of the element)
        trafo_igetReferenceDimension = 3
      
      case (TRAFO_ID_MLINCUBE)
        ! Bilinear transformation for cubic-shaped elements 
        ! -> Bilinear quadrilateral transformation.
        ! 4 DOF`s in the transformation (given by the corners of the element)
        trafo_igetReferenceDimension = 2
      
      end select
      
    case (NDIM3D)
      ! 3D elements. Tetrahedra, Hexahedra. Check the actual transformation
      ! ID how to transform.
    
      select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
      
      case (TRAFO_ID_LINSIMPLEX)
        ! 3D simplex -> 4 dimensions (1 more than space)
        trafo_igetReferenceDimension = 4
      
      case (TRAFO_ID_MLINCUBE)
        ! Trilinear hexahedral transformation -> dimension = space dimension
        trafo_igetReferenceDimension = 3
      
      case (TRAFO_ID_MLINPYRAMID)
        ! Bilinear pyramid transformation -> dimension = space dimension
        trafo_igetReferenceDimension = 3
      
      case (TRAFO_ID_MLINPRISM)
        ! Bilinear prism transformation -> dimension = space dimension
        trafo_igetReferenceDimension = 3
      
      end select
      
    end select

  end function

  ! ***************************************************************************

!<function>  

  elemental integer function trafo_igetNVE(ctrafoType)

!<description>
  ! This function returns for a given trynsformation ID the size NVE of the
  ! coordinate array which is used to hold all the coordinates of the vertices
  ! of an element (corners, edge midpoints,...).
!</description>

!<input>    
  ! The transformation code identifier.
  integer(I32), intent(in) :: ctrafoType
!</input>

!<result>
  ! The number vertices that are necessary to describe the transformation.
  ! (Note that this is the number of degrees of freedom of the transformation
  ! formula!)
!</result>

!</function>

    trafo_igetNVE = 0

    ! What type of transformation do we have? First decide on the dimension,
    ! then on the actual ID.
    select case (trafo_igetDimension(ctrafoType))

    case (NDIM1D)
      ! 1D elements. Lines. Check the actual transformation
      ! ID how to transform.
    
      select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
      
      case (TRAFO_ID_MLINCUBE)
        ! 1D linear line transformation. 
        ! 2 DOFs in the transformation (given by the corners of the element)
        trafo_igetNVE = 2
        
      end select
    
    case (NDIM2D)
      ! 2D elements. Triangles, Quadrilaterals. Check the actual transformation
      ! ID how to transform.
    
      select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
      
      case (TRAFO_ID_LINSIMPLEX)
        ! 2D simplex -> linear triangular transformation. 
        ! 3 DOFs in the transformation (given by the corners of the element)
        trafo_igetNVE = 3
      
      case (TRAFO_ID_MLINCUBE)
        ! Bilinear transformation for cubic-shaped elements 
        ! -> Bilinear quadrilateral transformation.
        ! 4 DOFs in the transformation (given by the corners of the element)
        trafo_igetNVE = 4
      
      end select
      
    case (NDIM3D)
      ! 3D elements. Tetrahedra, Hexahedra. Check the actual transformation
      ! ID how to transform.
    
      select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
      
      case (TRAFO_ID_LINSIMPLEX)
        ! 3D simplex -> linear triangular transformation. 
        ! 4 DOFs in the transformation (given by the corners of the element)
        trafo_igetNVE = 4
      
      case (TRAFO_ID_MLINCUBE)
        ! Trilinear transformation for cubic-shaped elements 
        ! -> Trilinear hexahedral transformation.
        ! 8 DOFs in the transformation (given by the corners of the element)
        trafo_igetNVE = 8
      
      case (TRAFO_ID_MLINPYRAMID)
        ! Bilinear transformation for pyramid shaped elements
        ! 5 DOFs in the transformation (given by the corners of the element)
        trafo_igetNVE = 5
      
      case (TRAFO_ID_MLINPRISM)
        ! Bilinear transformation for prismic shaped elements
        ! 6 DOFs in the transformation (given by the corners of the element)
        trafo_igetNVE = 6

      end select
      
    end select

  end function

! **********************************************************************

!<subroutine>

  subroutine trafo_getCoords (ctrafoType,rtriangulation,iel,Dcoords)

!<description>
  ! This routine filld the Dcoords array with the coordinates (of corners,
  ! edge midpoints,...) of element iel in the triangulation rtriangulation.
  ! Dcoords then defines the actual transformation from the reference
  ! to the real element.
!</description>

!<input>
  ! ID of transformation 
  integer(I32), intent(in) :: ctrafoType

  ! Triangulation structure that defines the elements.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! Number of the element whose information is to be fetched.
  integer, intent(in) :: iel
!</input>

!<output>
  ! Coordinates of the corners of the element
  !  Dcoord(1,i) = x-coordinates of corner i on an element, 
  !  Dcoord(2,i) = y-coordinates of corner i on an element.
  ! DIMENSION(#space dimensions,NVE).
  ! NVE depends on the type of transformation and can be determined with
  real(DP), dimension(:,:), intent(out) :: Dcoords
!</output>

!</subroutine>

    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! What type of transformation do we have? First decide on the dimension,
    ! then on the actual ID.
    select case (trafo_igetDimension(ctrafoType))
    
    case (NDIM1D)
      ! 1D elements. Lines
      ! We always need the corner-coordinate array, so get it 
      ! from the triangulation.
      call storage_getbase_double2d (rtriangulation%h_DvertexCoords,&  
                                     p_DvertexCoords)
      call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&  
                                  p_IverticesAtElement)
      
      ! Check the actual transformation ID how to transform.
    
      select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
      
      case (TRAFO_ID_MLINCUBE)
        ! 1D simplex -> linear line transformation. 
        ! Transfer the corners of the element.
        Dcoords (1,1:3) = p_DvertexCoords(1, p_IverticesAtElement(1:3, iel))
                                
      end select
    
    case (NDIM2D)
      ! 2D elements. Triangles, Quadrilaterals. 
      ! We always need the corner-coordinate array, so get it 
      ! from the triangulation.
      call storage_getbase_double2d (rtriangulation%h_DvertexCoords,&  
                                     p_DvertexCoords)
      call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&  
                                  p_IverticesAtElement)
      
      ! Check the actual transformation ID how to transform.
    
      select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
      
      case (TRAFO_ID_LINSIMPLEX)
        ! 2D simplex -> linear triangular transformation. 
        ! Transfer the corners of the element.
        Dcoords (1:NDIM2D,1:3) = p_DvertexCoords(1:NDIM2D,&
                                    p_IverticesAtElement(1:3,iel))
      
      case (TRAFO_ID_MLINCUBE)
        ! Bilinear transformation for cubic-shaped elements 
        ! -> Bilinear quadrilateral transformation.
        ! Transfer the corners of the element.
        Dcoords (1:NDIM2D,1:4) = p_DvertexCoords(1:NDIM2D,&
                                    p_IverticesAtElement(1:4,iel))
      
      end select

    case (NDIM3D)
      ! 3D elements. Tetrahedra, Hexahedra.
      ! We always need the corner-coordinate array, so get it 
      ! from the triangulation.
      call storage_getbase_double2d (rtriangulation%h_DvertexCoords,&  
                                     p_DvertexCoords)
      call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&  
                                  p_IverticesAtElement)
      
      ! Check the actual transformation ID how to transform.
    
      select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
      
      case (TRAFO_ID_LINSIMPLEX)
        ! 3D simplex -> linear tetrahedral transformation. 
        ! Transfer the corners of the element.
        Dcoords (1:NDIM3D,1:4) = p_DvertexCoords(1:NDIM3D,&
                                    p_IverticesAtElement(1:4,iel))
      
      case (TRAFO_ID_MLINCUBE)
        ! Trilinear transformation for cubic-shaped elements 
        ! -> Trilinear hexahedral transformation.
        ! Transfer the corners of the element.
        Dcoords (1:NDIM3D,1:8) = p_DvertexCoords(1:NDIM3D,&
                                    p_IverticesAtElement(1:8,iel))
      
      case (TRAFO_ID_MLINPYRAMID)
        ! Bilinear transformation for pyramid shaped elements 
        ! Transfer the corners of the element.
        Dcoords (1:NDIM3D,1:5) = p_DvertexCoords(1:NDIM3D,&
                                    p_IverticesAtElement(1:5,iel))

      case (TRAFO_ID_MLINPRISM)
        ! Bilinear transformation for prismic shaped elements 
        ! Transfer the corners of the element.
        Dcoords (1:NDIM3D,1:6) = p_DvertexCoords(1:NDIM3D,&
                                    p_IverticesAtElement(1:6,iel))

      end select
      
    end select
      
  end subroutine

! **********************************************************************

!<subroutine>

  subroutine trafo_getCoords_sim (ctrafoType,rtriangulation,Ielements,Dcoords)

!<description>
  ! This routine fills the Dcoords array with the coordinates (of corners,
  ! edge midpoints,...) of all the elements given in the element
  ! list Ielements, based on the triangulation rtriangulation.
  ! Dcoords then defines the actual transformation from the reference
  ! to the real elements.
!</description>

!<input>
  ! ID of transformation 
  integer(I32), intent(in) :: ctrafoType

  ! Triangulation structure that defines the elements.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! Array with element numbers whose coordinates should be extracted into
  ! Dcoords.
  integer, dimension(:), intent(in) :: Ielements
!</input>

!<output>
  ! Coordinates of the corners of all the elements
  !  Dcoord(1,i,.) = x-coordinates of corner i on an element, 
  !  Dcoord(2,i,.) = y-coordinates of corner i on an element.
  ! DIMENSION(#space dimensions,NVE,size of Ielements)
  ! NVE depends on the type of transformation and can be determined with
  real(DP), dimension(:,:,:), intent(out) :: Dcoords
!</output>

!</subroutine>

    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer :: ipoint
    integer :: iel

    ! What type of transformation do we have? First decide on the dimension,
    ! then on the actual ID.
    select case (trafo_igetDimension(ctrafoType))
    
    case (NDIM1D)
      ! 2D elements. Lines.
      ! We always need the corner-coordinate array, so get it 
      ! from the triangulation.
      call storage_getbase_double2d (rtriangulation%h_DvertexCoords,&  
                                     p_DvertexCoords)
      call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&  
                                  p_IverticesAtElement)
      
      ! Check the actual transformation ID how to transform.
    
      select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
      
      case (TRAFO_ID_MLINCUBE)
        ! 1D simplex -> linear line transformation. 
        ! Transfer the corners of the element.
        !$omp parallel do private(ipoint) &
        !$omp if(size(Ielements) > TRAFO_NELEMMIN_OMP)
        do iel=1,size(Ielements)
          do ipoint = 1,2
            Dcoords (1,ipoint,iel) = &
              p_DvertexCoords(1,p_IverticesAtElement(ipoint,Ielements(iel)))
          end do
        end do
        !$omp end parallel do
      
      end select

    case (NDIM2D)
      ! 2D elements. Triangles, Quadrilaterals. 
      ! We always need the corner-coordinate array, so get it 
      ! from the triangulation.
      call storage_getbase_double2d (rtriangulation%h_DvertexCoords,&  
                                     p_DvertexCoords)
      call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&  
                                  p_IverticesAtElement)
      
      ! Check the actual transformation ID how to transform.
    
      select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
      
      case (TRAFO_ID_LINSIMPLEX)
        ! 2D simplex -> linear triangular transformation. 
        ! Transfer the corners of the element.
        !$omp parallel do private(ipoint) &
        !$omp if(size(Ielements) > TRAFO_NELEMMIN_OMP)
        do iel=1,size(Ielements)
          do ipoint = 1,3
            Dcoords (1:NDIM2D,ipoint,iel) = &
              p_DvertexCoords(1:NDIM2D,p_IverticesAtElement(ipoint,Ielements(iel)))
          end do
        end do
        !$omp end parallel do
      
      case (TRAFO_ID_MLINCUBE)
        ! Bilinear transformation for cubic-shaped elements 
        ! -> Bilinear quadrilateral transformation.
        ! Transfer the corners of the element.
        !$omp parallel do private(ipoint) &
        !$omp if(size(Ielements) > TRAFO_NELEMMIN_OMP)
        do iel=1,size(Ielements)
          do ipoint = 1,4
            Dcoords (1:NDIM2D,ipoint,iel) = &
              p_DvertexCoords(1:NDIM2D,p_IverticesAtElement(ipoint,Ielements(iel)))
          end do
        end do
        !$omp end parallel do
      
      end select

    case (NDIM3D)
      ! 3D elements. Tetrahedra, Hexahedra.
      ! We always need the corner-coordinate array, so get it 
      ! from the triangulation.
      call storage_getbase_double2d (rtriangulation%h_DvertexCoords,&  
                                     p_DvertexCoords)
      call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&  
                                  p_IverticesAtElement)
      
      ! Check the actual transformation ID how to transform.
    
      select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
      
      case (TRAFO_ID_LINSIMPLEX)
        ! 3D simplex -> linear tetrahedral transformation. 
        ! Transfer the corners of the element.
        !$omp parallel do private(ipoint) &
        !$omp if(size(Ielements) > TRAFO_NELEMMIN_OMP)
        do iel=1,size(Ielements)
          do ipoint = 1,4
            Dcoords (1:NDIM3D,ipoint,iel) = &
              p_DvertexCoords(1:NDIM3D,p_IverticesAtElement(ipoint,Ielements(iel)))
          end do
        end do
        !$omp end parallel do
      
      case (TRAFO_ID_MLINCUBE)
        ! Trilinear transformation for cubic-shaped elements 
        ! -> Trilinear hexahedral transformation.
        ! Transfer the corners of the element.
        !$omp parallel do private(ipoint) &
        !$omp if(size(Ielements) > TRAFO_NELEMMIN_OMP)
        do iel=1,size(Ielements)
          do ipoint = 1,8
            Dcoords (1:NDIM3D,ipoint,iel) = &
              p_DvertexCoords(1:NDIM3D,p_IverticesAtElement(ipoint,Ielements(iel)))
          end do
        end do
        !$omp end parallel do
      
      case (TRAFO_ID_MLINPYRAMID)
        ! Bilinear transformation for pyramid shaped elements 
        ! Transfer the corners of the element.
        !$omp parallel do private(ipoint) &
        !$omp if(size(Ielements) > TRAFO_NELEMMIN_OMP)
        do iel=1,size(Ielements)
          do ipoint = 1,5
            Dcoords (1:NDIM3D,ipoint,iel) = &
              p_DvertexCoords(1:NDIM3D,p_IverticesAtElement(ipoint,Ielements(iel)))
          end do
        end do
        !$omp end parallel do
      
      
      case (TRAFO_ID_MLINPRISM)
        ! Bilinear transformation for prismic shaped elements 
        ! Transfer the corners of the element.
        !$omp parallel do private(ipoint) &
        !$omp if(size(Ielements) > TRAFO_NELEMMIN_OMP)
        do iel=1,size(Ielements)
          do ipoint = 1,6
            Dcoords (1:NDIM3D,ipoint,iel) = &
              p_DvertexCoords(1:NDIM3D,p_IverticesAtElement(ipoint,Ielements(iel)))
          end do
        end do
        !$omp end parallel do

      end select
      
    end select
      
  end subroutine

! **********************************************************************

!<subroutine>

  pure subroutine trafo_calctrafo (ctrafoType,Dcoords,&
                              DpointRef,Djac,ddetj,DpointReal)
!<description>
  ! General transformation.
  !
  ! The aim of this routine is to calculate the transformation between the
  ! reference element for one point. The element is given
  ! as a list of corner points in Dcoords.
  !
  ! For the point DpointsRef, the following information is calculated:
  ! 1.) Determinant of the mapping from the reference to the real element,
  ! 2.) the Jacobian matrix of the mapping,
  ! 3.) if the parameter DpointsReal is present: coordinates of the mapped
  !     points on the real element(s).
!</description>

!<input>
  ! ID of transformation to calculate
  integer(I32), intent(in) :: ctrafoType

  ! Coordinates of the corners of the element
  !  Dcoord(1,i) = x-coordinates of corner i on an element, 
  !  Dcoord(2,i) = y-coordinates of corner i on an element.
  ! DIMENSION(#space dimensions,NVE)
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Coordinates of the point on the reference element 
  !  DpointRef(1) = x-coordinates of point i on an element, 
  !  DpointRef(2) = y-coordinates of point i on an element.
  ! DIMENSION(#space dimension)
  real(DP), dimension(:), intent(in) :: DpointRef
!</input>

!<output>
  ! The Jacobian matrix of the mapping for each point.
  ! DIMENSION(number of entries in the matrix,npointsPerEl)
  real(DP), dimension(:), intent(out) :: Djac
  
  ! Jacobian determinant of the mapping.
  ! DIMENSION(npointsPerEl)
  real(DP), intent(out) :: ddetj
  
  ! OPTIONAL: Array receiving the coordinates of the points in DpointRef,
  ! mapped from the reference element to the real element.
  ! If not specified, they are not computed.
  ! DIMENSION(#space dimension,npointsPerEl)
  real(DP), dimension(:), intent(out), optional :: DpointReal
!</output>

!</subroutine>

  ! local variables
  real(DP) :: dax, day, dbx, dby, dcx, dcy
  
  ! auxiliary factors for the bilinear quad mapping
  real(DP), dimension(TRAFO_NAUXJACMAX) :: DjacPrep
  
  ! What type of transformation do we have? First decide on the dimension,
  ! then on the actual ID.
  select case (trafo_igetDimension(ctrafoType))
  
  case (NDIM1D)
    
    select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
    case (TRAFO_ID_MLINCUBE)
    
      ! The Jacobi matrix is simply the length of the interval multiplied by 0.5
      Djac(1) = 0.5_DP * (Dcoords(1,2) - Dcoords(1,1))   
      ddetj = Djac(1)
      
      ! Do we need to transform the reference point?
      if(present(DpointReal)) DpointReal(1) = &
          Dcoords(1,1) + (DpointRef(1)+1.0_DP) * ddetj
    
    end select
  
  case (NDIM2D)
    ! 2D elements. Triangles, Quadrilaterals. Check the actual transformation
    ! ID how to transform.
  
    select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
    case (TRAFO_ID_LINSIMPLEX)

      ! 2D simplex -> linear triangular transformation.
      !
      ! Calculate the Jacobian matrix and determinant.
      !
      ! Example where to find information about barycentric coordinates:
      ! http://mathworld.wolfram.com/BarycentricCoordinates.html
      !
      ! The determinant is simply the polygonal area of the parallelogram
      ! given by the two vectors (c-a,b-a); c.f.
      !
      ! http://mathworld.wolfram.com/Parallelogram.html
      
      dax = Dcoords(1, 1)
      day = Dcoords(2, 1)
      dbx = Dcoords(1, 2)
      dby = Dcoords(2, 2)
      dcx = Dcoords(1, 3)
      dcy = Dcoords(2, 3)
      
      Djac(1)=dbx-dax
      Djac(2)=dby-day
      Djac(3)=dcx-dax
      Djac(4)=dcy-day
      
      ddetj = Djac(1)*Djac(4) - Djac(2)*Djac(3)
        
      ! Calculate coordinates?
      if (present(DpointReal)) then
        ! Ok, that was easy. It is slightly more complicated to get
        ! the matrix...
        ! But as long as the matrix is not needed, we skip the calculation -
        ! this might be done in a future implementation!
        !
        ! Calculation of the real coordinates is also easy.
        DpointReal(1) = DpointRef(1)*dax &
                            + DpointRef(2)*dbx &
                            + DpointRef(3)*dcx 
        DpointReal(2) = DpointRef(1)*day &
                            + DpointRef(2)*dby &
                            + DpointRef(3)*dcy 
          
      end if
          
    case (TRAFO_ID_MLINCUBE)
    
      ! Prepare the calculation of the Jacobi determinants
      call trafo_prepJac_quad2d(Dcoords, DjacPrep)

      ! Calculate with or without coordinates?
      if (.not. present(DpointReal)) then
      
        ! Calculate the Jacobian matrix and determinant
        call trafo_calcJac_quad2d (DjacPrep,Djac(:),ddetj, &
                              DpointRef(1),DpointRef(2))
      
      else

        ! Calculate the Jacobian matrix and determinant
        call trafo_calcTrafo_quad2d (DjacPrep,Djac(:),ddetj, &
                                DpointRef(1),DpointRef(2), &
                                DpointReal(1),DpointReal(2))
      end if
      
    end select
  
  case (NDIM3D)
  
    ! 3D elements. Tetra-, Hexahedra. Check the actual transformation
    ! ID how to transform.
  
    select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
    case (TRAFO_ID_LINSIMPLEX)

      ! 3D simplex -> linear tetrahedral transformation.
      !
      ! Calculate the Jacobian matrix and determinant.
      Djac(1)=Dcoords(1,2)-Dcoords(1,1)
      Djac(2)=Dcoords(2,2)-Dcoords(2,1)
      Djac(3)=Dcoords(3,2)-Dcoords(3,1)
      Djac(4)=Dcoords(1,3)-Dcoords(1,1)
      Djac(5)=Dcoords(2,3)-Dcoords(2,1)
      Djac(6)=Dcoords(3,3)-Dcoords(3,1)
      Djac(7)=Dcoords(1,4)-Dcoords(1,1)
      Djac(8)=Dcoords(2,4)-Dcoords(2,1)
      Djac(9)=Dcoords(3,4)-Dcoords(3,1)
      
      ddetj = Djac(1)*(Djac(5)*Djac(9) - Djac(6)*Djac(8)) &
            + Djac(2)*(Djac(6)*Djac(7) - Djac(4)*Djac(9)) &
            + Djac(3)*(Djac(4)*Djac(8) - Djac(5)*Djac(7))
        
      ! Calculate coordinates?
      if (present(DpointReal)) then
        ! Calculation of the real coordinates is also easy.
        DpointReal(1) = DpointRef(1)*Dcoords(1,1) &
                            + DpointRef(2)*Dcoords(1,2) &
                            + DpointRef(3)*Dcoords(1,3) &
                            + DpointRef(4)*Dcoords(1,4)
        DpointReal(2) = DpointRef(1)*Dcoords(2,1) &
                            + DpointRef(2)*Dcoords(2,2) &
                            + DpointRef(3)*Dcoords(2,3) &
                            + DpointRef(4)*Dcoords(2,4)
        DpointReal(3) = DpointRef(1)*Dcoords(3,1) &
                            + DpointRef(2)*Dcoords(3,2) &
                            + DpointRef(3)*Dcoords(3,3) &
                            + DpointRef(4)*Dcoords(3,4)
          
      end if

    case (TRAFO_ID_MLINCUBE)
    
      ! Prepare the calculation of the Jacobi determinants
      call trafo_prepJac_hexa3d(Dcoords, DjacPrep)

      ! Calculate with or without coordinates?
      if (.not. present(DpointReal)) then
        
        ! Calculate the Jacobian matrix and determinant
        call trafo_calcJac_hexa3d (DjacPrep,Djac(:),ddetj, &
                              DpointRef(1),DpointRef(2),DpointRef(3))
      
      else
        
        ! Calculate the Jacobian matrix and determinant
        call trafo_calcTrafo_hexa3d (DjacPrep,Djac(:),ddetj, &
                                DpointRef(1),DpointRef(2),DpointRef(3),&
                                DpointReal(1),DpointReal(2),DpointReal(3))
      end if

    case (TRAFO_ID_MLINPYRAMID)
    
      ! Prepare the calculation of the Jacobi determinants
      call trafo_prepJac_pyra3d(Dcoords, DjacPrep)

      ! Calculate with or without coordinates?
      if (.not. present(DpointReal)) then
        
        ! Calculate the Jacobian matrix and determinant
        call trafo_calcJac_pyra3d (DjacPrep,Djac(:),ddetj, &
                              DpointRef(1),DpointRef(2),DpointRef(3))
      
      else
        
        ! Calculate the Jacobian matrix and determinant
        call trafo_calcTrafo_pyra3d (DjacPrep,Djac(:),ddetj, &
                                DpointRef(1),DpointRef(2),DpointRef(3),&
                                DpointReal(1),DpointReal(2),DpointReal(3))
      end if
      
    case (TRAFO_ID_MLINPRISM)
    
      ! Prepare the calculation of the Jacobi determinants
      call trafo_prepJac_prism3d(Dcoords, DjacPrep)

      ! Calculate with or without coordinates?
      if (.not. present(DpointReal)) then
        
        ! Calculate the Jacobian matrix and determinant
        call trafo_calcJac_prism3d (DjacPrep,Djac(:),ddetj, &
                              DpointRef(1),DpointRef(2),DpointRef(3))
      
      else
        
        ! Calculate the Jacobian matrix and determinant
        call trafo_calcTrafo_prism3d (DjacPrep,Djac(:),ddetj, &
                                DpointRef(1),DpointRef(2),DpointRef(3),&
                                DpointReal(1),DpointReal(2),DpointReal(3))
      end if
    
    end select
  
  end select

  end subroutine

! **********************************************************************

!<subroutine>

  subroutine trafo_calctrafo_mult (ctrafoType,npointsPerEl,Dcoords,&
                                   DpointsRef,Djac,Ddetj,DpointsReal)
!<description>
  ! General transformation support for multiple points on one element.
  !
  ! The aim of this routine is to calculate the transformation between the
  ! reference element for multiple points. The element is given
  ! as a list of corner points in Dcoords.
  !
  ! For every of these npointsPerEl points in the element specified 
  ! by DpointsRef, the following information is calculated:
  ! 1.) Determinant of the mapping from the reference to the real element,
  ! 2.) the Jacobian matrix of the mapping,
  ! 3.) if the parameter DpointsReal is present: coordinates of the mapped
  !     points on the real element(s).
!</description>

!<input>
  ! ID of transformation to calculate
  integer(I32), intent(in) :: ctrafoType

  ! Number of points in each element where to calculate the transformation
  integer, intent(in) :: npointsPerEl

  ! Coordinates of the corners of all the elements
  !  Dcoord(1,i) = x-coordinates of corner i on an element, 
  !  Dcoord(2,i) = y-coordinates of corner i on an element.
  ! DIMENSION(#space dimensions,NVE)
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Coordinates of the points on the reference element for each element 
  ! where to calculate the mapping.
  !  DpointsRef(1,i) = x-coordinates of point i on an element, 
  !  DpointsRef(2,i) = y-coordinates of point i on an element.
  ! ! DIMENSION(#space dimension,npointsPerEl)
  real(DP), dimension(:,:), intent(in) :: DpointsRef
!</input>

!<output>
  ! The Jacobian matrix of the mapping for each point.
  ! DIMENSION(number of entries in the matrix,npointsPerEl)
  real(DP), dimension(:,:), intent(out) :: Djac
  
  ! Jacobian determinants of the mapping for all the points from the
  ! reference element to the real element.
  ! DIMENSION(npointsPerEl)
  real(DP), dimension(:), intent(out) :: Ddetj
  
  ! OPTIONAL: Array receiving the coordinates of the points in DpointsRef,
  ! mapped from the reference element to the real element.
  ! If not specified, they are not computed.
  ! DIMENSION(#space dimension,npointsPerEl)
  real(DP), dimension(:,:), intent(out), optional :: DpointsReal
!</output>

!</subroutine>

  ! local variables
  integer :: ipt
  real(DP) :: dax, day, dbx, dby, dcx, dcy
  
  ! auxiliary factors for the bilinear quad mapping
  real(DP), dimension(TRAFO_NAUXJACMAX) :: DjacPrep
  
  ! What type of transformation do we have? First decide on the dimension,
  ! then on the actual ID.
  select case (trafo_igetDimension(ctrafoType))
  
  case (NDIM1D)
    ! 1D elements: Lines.
    
    select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
    case (TRAFO_ID_MLINCUBE)
    
      if(.not. present(DpointsReal)) then

        ! Loop over the points
        do ipt=1,npointsPerEl
        
          ! The Jacobi matrix is simply the length of the interval multiplied by 0.5
          Djac(1,ipt) = 0.5_DP * (Dcoords(1,2) - Dcoords(1,1))   
          ddetj(ipt) = Djac(1,ipt)
          
        end do
        
      else
      
        ! Loop over the points
        do ipt=1,npointsPerEl
        
          ! The Jacobi matrix is simply the length of the interval multiplied by 0.5
          Djac(1,ipt) = 0.5_DP * (Dcoords(1,2) - Dcoords(1,1))   
          ddetj(ipt) = Djac(1,ipt)
          
          ! Transform the reference point into real coordinates
          DpointsReal(1,ipt) = Dcoords(1,1) + (DpointsRef(1,ipt)+1.0_DP) * ddetj(ipt)
          
        end do

      end if
    
    end select
  
  case (NDIM2D)
    ! 2D elements. Triangles, Quadrilaterals. Check the actual transformation
    ! ID how to transform.
  
    select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
    
    case (TRAFO_ID_LINSIMPLEX)

      ! 2D simplex -> linear triangular transformation.
      !
      ! Calculate with or without coordinates?
      if (.not. present(DpointsReal)) then
      
        ! Loop over the points
        do ipt=1,npointsPerEl
        
          ! Calculate the Jacobian matrix and determinant.
          !
          ! Example where to find information about barycentric coordinates:
          ! http://mathworld.wolfram.com/BarycentricCoordinates.html
          !
          ! The determinant is simply the polygonal area of the parallelogram
          ! given by the two vectors (c-a,b-a); c.f.
          !
          ! http://mathworld.wolfram.com/Parallelogram.html
          
          dax = Dcoords(1, 1)
          day = Dcoords(2, 1)
          dbx = Dcoords(1, 2)
          dby = Dcoords(2, 2)
          dcx = Dcoords(1, 3)
          dcy = Dcoords(2, 3)
          
          Djac(1,ipt)=dbx-dax
          Djac(2,ipt)=dby-day
          Djac(3,ipt)=dcx-dax
          Djac(4,ipt)=dcy-day
          
          Ddetj(ipt) = Djac(1,ipt)*Djac(4,ipt) &
                     - Djac(2,ipt)*Djac(3,ipt)
          
        end do ! ipt
        
    else

        ! Loop over the points
        do ipt=1,npointsPerEl
        
          ! Calculate the Jacobian matrix and determinant.
          !
          ! Example where to find information about barycentric coordinates:
          ! http://mathworld.wolfram.com/BarycentricCoordinates.html
          !
          ! The determinant is simply the polygonal area of the parallelogram
          ! given by the two vectors (c-a,b-a); c.f.
          !
          ! http://mathworld.wolfram.com/Parallelogram.html
          
          dax = Dcoords(1, 1)
          day = Dcoords(2, 1)
          dbx = Dcoords(1, 2)
          dby = Dcoords(2, 2)
          dcx = Dcoords(1, 3)
          dcy = Dcoords(2, 3)
          
          Djac(1,ipt)=dbx-dax
          Djac(2,ipt)=dby-day
          Djac(3,ipt)=dcx-dax
          Djac(4,ipt)=dcy-day
          
          Ddetj(ipt) = Djac(1,ipt)*Djac(4,ipt) &
                     - Djac(2,ipt)*Djac(3,ipt)
          
          ! Ok, that was easy. It is slightly more complicated to get
          ! the matrix...
          ! But as long as the matrix is not needed, we skip the calculation -
          ! this might be done in a future implementation!
          !
          ! Calculation of the real coordinates is also easy.
          DpointsReal(1,ipt) = DpointsRef(1,ipt)*dax &
                             + DpointsRef(2,ipt)*dbx &
                             + DpointsRef(3,ipt)*dcx 
          DpointsReal(2,ipt) = DpointsRef(1,ipt)*day &
                             + DpointsRef(2,ipt)*dby &
                             + DpointsRef(3,ipt)*dcy 
          
        end do ! ipt
          
      end if
          
    case (TRAFO_ID_MLINCUBE)
    
      ! Prepare the calculation of the Jacobi determinants
      call trafo_prepJac_quad2D(Dcoords, DjacPrep)
        
      ! Calculate with or without coordinates?
      if (.not. present(DpointsReal)) then
      
        ! Loop over the points
        do ipt=1,npointsPerEl
          ! Calculate the Jacobian matrix and determinant
          call trafo_calcJac_quad2D (DjacPrep,Djac(:,ipt),Ddetj(ipt), &
                               DpointsRef(1,ipt),DpointsRef(2,ipt))
        end do ! ipt
      
      else

        ! Loop over the points
        do ipt=1,npointsPerEl
          ! Calculate the Jacobian matrix and determinant
          call trafo_calcTrafo_quad2d (DjacPrep,Djac(:,ipt),Ddetj(ipt), &
                                 DpointsRef(1,ipt),DpointsRef(2,ipt), &
                                 DpointsReal(1,ipt),DpointsReal(2,ipt))
        end do ! ipt

      end if
      
    end select
 
  case (NDIM3D)
  
    ! 3D elements. Tetra-, Hexahedra. Check the actual transformation
    ! ID how to transform.
  
    select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
    case (TRAFO_ID_LINSIMPLEX)

      ! 3D simplex -> linear tetrahedral transformation.
      !
      ! Calculate with or without coordinates?
      if (.not. present(DpointsReal)) then
      
        ! Loop over the points
        do ipt=1,npointsPerEl
        
          ! Calculate the Jacobian matrix and determinant.
          Djac(1,ipt)=Dcoords(1,2)-Dcoords(1,1)
          Djac(2,ipt)=Dcoords(2,2)-Dcoords(2,1)
          Djac(3,ipt)=Dcoords(3,2)-Dcoords(3,1)
          Djac(4,ipt)=Dcoords(1,3)-Dcoords(1,1)
          Djac(5,ipt)=Dcoords(2,3)-Dcoords(2,1)
          Djac(6,ipt)=Dcoords(3,3)-Dcoords(3,1)
          Djac(7,ipt)=Dcoords(1,4)-Dcoords(1,1)
          Djac(8,ipt)=Dcoords(2,4)-Dcoords(2,1)
          Djac(9,ipt)=Dcoords(3,4)-Dcoords(3,1)
          
          Ddetj(ipt) = Djac(1,ipt)*(Djac(5,ipt)*Djac(9,ipt) &
                              - Djac(6,ipt)*Djac(8,ipt)) &
                     + Djac(2,ipt)*(Djac(6,ipt)*Djac(7,ipt) &
                              - Djac(4,ipt)*Djac(9,ipt)) &
                     + Djac(3,ipt)*(Djac(4,ipt)*Djac(8,ipt) &
                              - Djac(5,ipt)*Djac(7,ipt))
        end do
        
      else    
          
        ! Loop over the points
        do ipt=1,npointsPerEl
        
          ! Calculate the Jacobian matrix and determinant.
          Djac(1,ipt)=Dcoords(1,2)-Dcoords(1,1)
          Djac(2,ipt)=Dcoords(2,2)-Dcoords(2,1)
          Djac(3,ipt)=Dcoords(3,2)-Dcoords(3,1)
          Djac(4,ipt)=Dcoords(1,3)-Dcoords(1,1)
          Djac(5,ipt)=Dcoords(2,3)-Dcoords(2,1)
          Djac(6,ipt)=Dcoords(3,3)-Dcoords(3,1)
          Djac(7,ipt)=Dcoords(1,4)-Dcoords(1,1)
          Djac(8,ipt)=Dcoords(2,4)-Dcoords(2,1)
          Djac(9,ipt)=Dcoords(3,4)-Dcoords(3,1)
          
          Ddetj(ipt) = Djac(1,ipt)*(Djac(5,ipt)*Djac(9,ipt) &
                              - Djac(6,ipt)*Djac(8,ipt)) &
                     + Djac(2,ipt)*(Djac(6,ipt)*Djac(7,ipt) &
                              - Djac(4,ipt)*Djac(9,ipt)) &
                     + Djac(3,ipt)*(Djac(4,ipt)*Djac(8,ipt) &
                              - Djac(5,ipt)*Djac(7,ipt))
            
          ! Calculation of the real coordinates is also easy.
          DpointsReal(1,ipt) = DpointsRef(1,ipt)*Dcoords(1,1) &
                             + DpointsRef(2,ipt)*Dcoords(1,2) &
                             + DpointsRef(3,ipt)*Dcoords(1,3) &
                             + DpointsRef(4,ipt)*Dcoords(1,4)
          DpointsReal(2,ipt) = DpointsRef(1,ipt)*Dcoords(2,1) &
                             + DpointsRef(2,ipt)*Dcoords(2,2) &
                             + DpointsRef(3,ipt)*Dcoords(2,3) &
                             + DpointsRef(4,ipt)*Dcoords(2,4)
          DpointsReal(3,ipt) = DpointsRef(1,ipt)*Dcoords(3,1) &
                             + DpointsRef(2,ipt)*Dcoords(3,2) &
                             + DpointsRef(3,ipt)*Dcoords(3,3) &
                             + DpointsRef(4,ipt)*Dcoords(3,4)
        end do
      end if

    case (TRAFO_ID_MLINCUBE)
    
      ! Prepare the calculation of the Jacobi determinants
      call trafo_prepJac_hexa3D(Dcoords, DjacPrep)

      ! Calculate with or without coordinates?
      if (.not. present(DpointsReal)) then
      
        ! Loop over the points
        do ipt=1,npointsPerEl
          ! Calculate the Jacobian matrix and determinant
          call trafo_calcJac_hexa3D (DjacPrep,Djac(:,ipt),Ddetj(ipt), &
              DpointsRef(1,ipt),DpointsRef(2,ipt),DpointsRef(3,ipt))
        end do ! ipt
        
      else
        
        ! Loop over the points
        do ipt=1,npointsPerEl
          ! Calculate the Jacobian matrix and determinant
          call trafo_calcTrafo_hexa3d(DjacPrep,Djac(:,ipt),Ddetj(ipt),&
              DpointsRef(1,ipt),DpointsRef(2,ipt),DpointsRef(3,ipt),&
              DpointsReal(1,ipt),DpointsReal(2,ipt),DpointsReal(3,ipt))
        end do ! ipt
        
      end if

    case (TRAFO_ID_MLINPYRAMID)
    
      ! Prepare the calculation of the Jacobi determinants
      call trafo_prepJac_pyra3D(Dcoords, DjacPrep)

      ! Calculate with or without coordinates?
      if (.not. present(DpointsReal)) then
      
        ! Loop over the points
        do ipt=1,npointsPerEl
          ! Calculate the Jacobian matrix and determinant
          call trafo_calcJac_pyra3D (DjacPrep,Djac(:,ipt),Ddetj(ipt), &
              DpointsRef(1,ipt),DpointsRef(2,ipt),DpointsRef(3,ipt))
        end do ! ipt
        
      else
        
        ! Loop over the points
        do ipt=1,npointsPerEl
          ! Calculate the Jacobian matrix and determinant
          call trafo_calcTrafo_pyra3d(DjacPrep,Djac(:,ipt),Ddetj(ipt),&
              DpointsRef(1,ipt),DpointsRef(2,ipt),DpointsRef(3,ipt),&
              DpointsReal(1,ipt),DpointsReal(2,ipt),DpointsReal(3,ipt))
        end do ! ipt
        
      end if

    case (TRAFO_ID_MLINPRISM)
    
      ! Prepare the calculation of the Jacobi determinants
      call trafo_prepJac_prism3D(Dcoords, DjacPrep)

      ! Calculate with or without coordinates?
      if (.not. present(DpointsReal)) then
      
        ! Loop over the points
        do ipt=1,npointsPerEl
          ! Calculate the Jacobian matrix and determinant
          call trafo_calcJac_prism3D (DjacPrep,Djac(:,ipt),Ddetj(ipt), &
              DpointsRef(1,ipt),DpointsRef(2,ipt),DpointsRef(3,ipt))
        end do ! ipt
        
      else
        
        ! Loop over the points
        do ipt=1,npointsPerEl
          ! Calculate the Jacobian matrix and determinant
          call trafo_calcTrafo_prism3d(DjacPrep,Djac(:,ipt),Ddetj(ipt),&
              DpointsRef(1,ipt),DpointsRef(2,ipt),DpointsRef(3,ipt),&
              DpointsReal(1,ipt),DpointsReal(2,ipt),DpointsReal(3,ipt))
        end do ! ipt
        
      end if
    
    end select
      
  end select

  end subroutine

! **********************************************************************

!<subroutine>

  subroutine trafo_calctrafo_sim (ctrafoType,nelements,npointsPerEl,Dcoords,&
                                  DpointsRef,Djac,Ddetj,DpointsReal)

!<description>
  ! General transformation support for multiple points on multiple
  ! elements. 
  !
  ! The aim of this routine is to calculate the transformation between the
  ! reference element and multiple real elements. The elements are given
  ! as a list of corner points in Dcoords.
  !
  ! On every of the nelements elements given in this list, there are
  ! npointsPerEl points inside the element given in reference coordinates.
  ! For every of these npointsPerEl*nelements points, the following 
  ! information is calculated:
  ! 1.) Determinant of the mapping from the reference to the real element,
  ! 2.) the Jacobian matrix of the mapping,
  ! 3.) if the parameter DpointsReal is present: coordinates of the mapped
  !     points on the real element(s).
!</description>

!<input>
  ! ID of transformation to calculate
  integer(I32), intent(in) :: ctrafoType

  ! Number of elements where to calculate the transformation
  integer, intent(in) :: nelements
  
  ! Number of points in each element where to calculate the transformation
  integer, intent(in) :: npointsPerEl

  ! Coordinates of the corners of all the elements
  !  Dcoord(1,i,.) = x-coordinates of corner i on an element, 
  !  Dcoord(2,i,.) = y-coordinates of corner i on an element.
  ! DIMENSION(#space dimensions,NVE,nelements)
  real(DP), dimension(:,:,:), intent(in) :: Dcoords

  ! Coordinates of the points on the reference element for each element 
  ! where to calculate the mapping.
  ! DIMENSION(#space dimensions,npointsPerEl,nelements) for quadrilateral elements and
  ! DIMENSION(#space dimensions+1,npointsPerEl,nelements) for triangular elements.
  !
  ! For QUAD elements:
  !  DpointsRef(1,i,.) = x-coordinates of point i on an element, 
  !  DpointsRef(2,i,.) = y-coordinates of point i on an element.
  !
  ! For triangular elements:
  !  DpointsRef(1,i,.) = First barycentric coordinate of point i on an element
  !  DpointsRef(2,i,.) = Second barycentric coordinate of point i on an element
  !  DpointsRef(3,i,.) = Third barycentric coordinate of point i on an element
  real(DP), dimension(:,:,:), intent(in) :: DpointsRef
!</input>

!<output>
  ! The Jacobian matrix of the mapping for each point.
  ! DIMENSION(number of entries in the matrix,npointsPerEl,nelements)
  real(DP), dimension(:,:,:), intent(out) :: Djac
  
  ! Jacobian determinants of the mapping for all the points from the
  ! reference element to the real element.
  ! DIMENSION(npointsPerEl,nelements)
  real(DP), dimension(:,:), intent(out) :: Ddetj
  
  ! OPTIONAL: Array receiving the coordinates of the points in DpointsRef,
  ! mapped from the reference element to the real element.
  ! If not specified, they are not computed.
  ! DIMENSION(#space dimensions,npointsPerEl,nelements)
  real(DP), dimension(:,:,:), intent(out), optional :: DpointsReal
!</output>

!</subroutine>

  ! local variables
  integer :: iel, ipt
  real(DP) :: dax, day, dbx, dby, dcx, dcy
  
  ! auxiliary factors for the bilinear quad mapping
  real(DP), dimension(TRAFO_NAUXJACMAX) :: DjacPrep
  
  ! What type of transformation do we have? First decide on the dimension,
  ! then on the actual ID.
  select case (trafo_igetDimension(ctrafoType))
  
  case (NDIM1D)
    ! 1D elements: Lines.
    
    select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
    case (TRAFO_ID_MLINCUBE)
    
      if(.not. present(DpointsReal)) then
        
        ! Loop over the elements
        !$omp parallel do private(ipt) &
        !$omp if(nelements > TRAFO_NELEMMIN_OMP)
        do iel=1,nelements
        
          ! Loop over the points
          do ipt=1,npointsPerEl
            
            ! The Jacobi matrix is simply the length of the interval multiplied by 0.5
            Djac(1,ipt,iel) = 0.5_DP * (Dcoords(1,2,iel) - Dcoords(1,1,iel))
            ddetj(ipt,iel) = Djac(1,ipt,iel)
              
          end do
        
        end do
        !$omp end parallel do
        
      else
      
        ! Loop over the elements
        !$omp parallel do private(ipt) &
        !$omp if(nelements > TRAFO_NELEMMIN_OMP)
        do iel=1,nelements
          
          ! Loop over the points
          do ipt=1,npointsPerEl
        
            ! The Jacobi matrix is simply the length of the interval multiplied by 0.5
            Djac(1,ipt,iel) = 0.5_DP * (Dcoords(1,2,iel) - Dcoords(1,1,iel))   
            ddetj(ipt,iel) = Djac(1,ipt,iel)
          
            ! Transform the reference point into real coordinates
            DpointsReal(1,ipt,iel) = Dcoords(1,1,iel) + &
                (DpointsRef(1,ipt,iel)+1.0_DP) * ddetj(ipt,iel)
          
          end do
          
        end do
        !$omp end parallel do

      end if
    
    end select
    
  case (NDIM2D)
    ! 2D elements. Triangles, Quadrilaterals. Check the actual transformation
    ! ID how to transform.
  
    select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
    
    case (TRAFO_ID_LINSIMPLEX)
    
      ! 2D simplex -> linear triangular transformation.
      !
      ! Calculate with or without coordinates?
      if (.not. present(DpointsReal)) then
      
        ! Loop over the elements
        !$omp parallel do private(ipt,dax,day,dbx,dby,dcx,dcy) &
        !$omp if(nelements > TRAFO_NELEMMIN_OMP)
        do iel = 1,nelements
          
          ! Loop over the points
          do ipt=1,npointsPerEl
          
            ! Calculate the Jacobian matrix and determinant.
            !
            ! Example where to find information about barycentric coordinates:
            ! http://mathworld.wolfram.com/BarycentricCoordinates.html
            !
            ! The determinant is simply the polygonal area of the parallelogram
            ! given by the two vectors (c-a,b-a); c.f.
            !
            ! http://mathworld.wolfram.com/Parallelogram.html
            
            dax = Dcoords(1, 1, iel)
            day = Dcoords(2, 1, iel)
            dbx = Dcoords(1, 2, iel)
            dby = Dcoords(2, 2, iel)
            dcx = Dcoords(1, 3, iel)
            dcy = Dcoords(2, 3, iel)
            
            Djac(1,ipt,iel)=dbx-dax
            Djac(2,ipt,iel)=dby-day
            Djac(3,ipt,iel)=dcx-dax
            Djac(4,ipt,iel)=dcy-day
            
            Ddetj(ipt,iel) = Djac(1,ipt,iel)*Djac(4,ipt,iel) &
                           - Djac(2,ipt,iel)*Djac(3,ipt,iel)
            
          end do ! ipt
          
        end do ! iel
        !$omp end parallel do
        
      else

        ! Loop over the elements
        !$omp parallel do private(ipt,dax,day,dbx,dby,dcx,dcy) &
        !$omp if(nelements > TRAFO_NELEMMIN_OMP)
        do iel = 1,nelements
          
          ! Loop over the points
          do ipt=1,npointsPerEl
          
            ! Calculate the Jacobian matrix and determinant.
            !
            ! Example where to find information about barycentric coordinates:
            ! http://mathworld.wolfram.com/BarycentricCoordinates.html
            !
            ! The determinant is simply the polygonal area of the parallelogram
            ! given by the two vectors (c-a,b-a); c.f.
            !
            ! http://mathworld.wolfram.com/Parallelogram.html
            
            dax = Dcoords(1, 1, iel)
            day = Dcoords(2, 1, iel)
            dbx = Dcoords(1, 2, iel)
            dby = Dcoords(2, 2, iel)
            dcx = Dcoords(1, 3, iel)
            dcy = Dcoords(2, 3, iel)
            
            Djac(1,ipt,iel)=dbx-dax
            Djac(2,ipt,iel)=dby-day
            Djac(3,ipt,iel)=dcx-dax
            Djac(4,ipt,iel)=dcy-day
            
            Ddetj(ipt,iel) = Djac(1,ipt,iel)*Djac(4,ipt,iel) &
                           - Djac(2,ipt,iel)*Djac(3,ipt,iel)
            
            ! Ok, that was easy. It is slightly more complicated to get
            ! the matrix...
            ! But as long as the matrix is not needed, we skip the calculation -
            ! this might be done in a future implementation!
            !
            ! Calculation of the real coordinates is also easy.
            DpointsReal(1,ipt,iel) = DpointsRef(1,ipt,iel)*dax &
                                   + DpointsRef(2,ipt,iel)*dbx &
                                   + DpointsRef(3,ipt,iel)*dcx 
            DpointsReal(2,ipt,iel) = DpointsRef(1,ipt,iel)*day &
                                   + DpointsRef(2,ipt,iel)*dby &
                                   + DpointsRef(3,ipt,iel)*dcy 
            
          end do ! ipt
          
        end do ! iel
        !$omp end parallel do
        
      end if
    
    case (TRAFO_ID_MLINCUBE)
    
      ! Bilinear transformation for cubic-shaped elements 
      ! -> Bilinear quadrilateral transformation.
      !
      ! Calculate with or without coordinates?
      if (.not. present(DpointsReal)) then
      
        ! Loop over the elements
        !$omp parallel do private(ipt,DjacPrep) &
        !$omp if(nelements > TRAFO_NELEMMIN_OMP)
        do iel = 1,nelements
          ! Prepare the calculation of the Jacobi determinants
          call trafo_prepJac_quad2D(Dcoords(:,:,iel), DjacPrep)
          
          ! Loop over the points
          do ipt=1,npointsPerEl
            ! Calculate the Jacobian matrix and determinant
            call trafo_calcJac_quad2D (DjacPrep,Djac(:,ipt,iel),Ddetj(ipt,iel), &
                                DpointsRef(1,ipt,iel),DpointsRef(2,ipt,iel))
          end do ! ipt
          
        end do ! iel
        !$omp end parallel do
      
      else

        ! Loop over the elements
        !$omp parallel do private(ipt,DjacPrep) &
        !$omp if(nelements > TRAFO_NELEMMIN_OMP)
        do iel = 1,nelements
          ! Prepare the calculation of the Jacobi determinants
          call trafo_prepJac_quad2D(Dcoords(:,:,iel), DjacPrep)
          
          ! Loop over the points
          do ipt=1,npointsPerEl
            ! Calculate the Jacobian matrix and determinant
            call trafo_calcTrafo_quad2d (DjacPrep,Djac(:,ipt,iel),Ddetj(ipt,iel), &
                                  DpointsRef(1,ipt,iel),DpointsRef(2,ipt,iel), &
                                  DpointsReal(1,ipt,iel),DpointsReal(2,ipt,iel))
          end do ! ipt
          
        end do ! iel
        !$omp end parallel do

      end if
      
    end select ! actual ID
 
  case (NDIM3D)
  
    ! 3D elements. Tetra-, Hexahedra. Check the actual transformation
    ! ID how to transform.
  
    select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
    case (TRAFO_ID_LINSIMPLEX)

      ! 3D simplex -> linear tetrahedral transformation.
      !
      ! Calculate with or without coordinates?
      if (.not. present(DpointsReal)) then
      
        ! Loop over the elements
        do iel=1, nelements

          ! Loop over the points
          do ipt=1,npointsPerEl
          
            ! Calculate the Jacobian matrix and determinant.
            Djac(1,ipt,iel)=Dcoords(1,2,iel)-Dcoords(1,1,iel)
            Djac(2,ipt,iel)=Dcoords(2,2,iel)-Dcoords(2,1,iel)
            Djac(3,ipt,iel)=Dcoords(3,2,iel)-Dcoords(3,1,iel)
            Djac(4,ipt,iel)=Dcoords(1,3,iel)-Dcoords(1,1,iel)
            Djac(5,ipt,iel)=Dcoords(2,3,iel)-Dcoords(2,1,iel)
            Djac(6,ipt,iel)=Dcoords(3,3,iel)-Dcoords(3,1,iel)
            Djac(7,ipt,iel)=Dcoords(1,4,iel)-Dcoords(1,1,iel)
            Djac(8,ipt,iel)=Dcoords(2,4,iel)-Dcoords(2,1,iel)
            Djac(9,ipt,iel)=Dcoords(3,4,iel)-Dcoords(3,1,iel)
            
            Ddetj(ipt,iel) = Djac(1,ipt,iel)*(Djac(5,ipt,iel)*Djac(9,ipt,iel) &
                                - Djac(6,ipt,iel)*Djac(8,ipt,iel)) &
                           + Djac(2,ipt,iel)*(Djac(6,ipt,iel)*Djac(7,ipt,iel) &
                                - Djac(4,ipt,iel)*Djac(9,ipt,iel)) &
                           + Djac(3,ipt,iel)*(Djac(4,ipt,iel)*Djac(8,ipt,iel) &
                                - Djac(5,ipt,iel)*Djac(7,ipt,iel))
           end do
           
        end do
        
      else    
        
        ! Loop over the elements
        do iel=1, nelements

          ! Loop over the points
          do ipt=1,npointsPerEl
          
            ! Calculate the Jacobian matrix and determinant.
            Djac(1,ipt,iel)=Dcoords(1,2,iel)-Dcoords(1,1,iel)
            Djac(2,ipt,iel)=Dcoords(2,2,iel)-Dcoords(2,1,iel)
            Djac(3,ipt,iel)=Dcoords(3,2,iel)-Dcoords(3,1,iel)
            Djac(4,ipt,iel)=Dcoords(1,3,iel)-Dcoords(1,1,iel)
            Djac(5,ipt,iel)=Dcoords(2,3,iel)-Dcoords(2,1,iel)
            Djac(6,ipt,iel)=Dcoords(3,3,iel)-Dcoords(3,1,iel)
            Djac(7,ipt,iel)=Dcoords(1,4,iel)-Dcoords(1,1,iel)
            Djac(8,ipt,iel)=Dcoords(2,4,iel)-Dcoords(2,1,iel)
            Djac(9,ipt,iel)=Dcoords(3,4,iel)-Dcoords(3,1,iel)
            
            Ddetj(ipt,iel) = Djac(1,ipt,iel)*(Djac(5,ipt,iel)*Djac(9,ipt,iel) &
                                - Djac(6,ipt,iel)*Djac(8,ipt,iel)) &
                           + Djac(2,ipt,iel)*(Djac(6,ipt,iel)*Djac(7,ipt,iel) &
                                - Djac(4,ipt,iel)*Djac(9,ipt,iel)) &
                           + Djac(3,ipt,iel)*(Djac(4,ipt,iel)*Djac(8,ipt,iel) &
                                - Djac(5,ipt,iel)*Djac(7,ipt,iel))
              
            ! Calculation of the real coordinates is also easy.
            DpointsReal(1,ipt,iel) = DpointsRef(1,ipt,iel)*Dcoords(1,1,iel) &
                                   + DpointsRef(2,ipt,iel)*Dcoords(1,2,iel) &
                                   + DpointsRef(3,ipt,iel)*Dcoords(1,3,iel) &
                                   + DpointsRef(4,ipt,iel)*Dcoords(1,4,iel)
            DpointsReal(2,ipt,iel) = DpointsRef(1,ipt,iel)*Dcoords(2,1,iel) &
                                   + DpointsRef(2,ipt,iel)*Dcoords(2,2,iel) &
                                   + DpointsRef(3,ipt,iel)*Dcoords(2,3,iel) &
                                   + DpointsRef(4,ipt,iel)*Dcoords(2,4,iel)
            DpointsReal(3,ipt,iel) = DpointsRef(1,ipt,iel)*Dcoords(3,1,iel) &
                                   + DpointsRef(2,ipt,iel)*Dcoords(3,2,iel) &
                                   + DpointsRef(3,ipt,iel)*Dcoords(3,3,iel) &
                                   + DpointsRef(4,ipt,iel)*Dcoords(3,4,iel)
          end do
          
        end do
        
      end if

    case (TRAFO_ID_MLINCUBE)
    
      ! Calculate with or without coordinates?
      if (.not. present(DpointsReal)) then
      
        !$omp parallel do private(ipt,DjacPrep) &
        !$omp if(nelements > TRAFO_NELEMMIN_OMP)
        do iel = 1,nelements
          ! Prepare the calculation of the Jacobi determinants
          call trafo_prepJac_hexa3D(Dcoords(:,:,iel), DjacPrep)
          
          ! Loop over the points
          do ipt=1,npointsPerEl
            ! Calculate the Jacobian matrix and determinant
            call trafo_calcJac_hexa3D (DjacPrep,Djac(:,ipt,iel),Ddetj(ipt,iel), &
                DpointsRef(1,ipt,iel),DpointsRef(2,ipt,iel),DpointsRef(3,ipt,iel))
          end do ! ipt
          
        end do !iel
        !$omp end parallel do
        
      else

        !$omp parallel do private(ipt,DjacPrep) &
        !$omp if(nelements > TRAFO_NELEMMIN_OMP)
        do iel = 1, nelements
          ! Prepare the calculation of the Jacobi determinants
          call trafo_prepJac_hexa3D(Dcoords(:,:,iel), DjacPrep)
          
          ! Loop over the points
          do ipt=1,npointsPerEl
            ! Calculate the Jacobian matrix and determinant
            call trafo_calcTrafo_hexa3d (DjacPrep,Djac(:,ipt,iel),Ddetj(ipt,iel), &
                DpointsRef(1,ipt,iel),DpointsRef(2,ipt,iel),DpointsRef(3,ipt,iel),&
                DpointsReal(1,ipt,iel),DpointsReal(2,ipt,iel),DpointsReal(3,ipt,iel))
          end do ! ipt
          
        end do ! iel
        !$omp end parallel do
        
      end if

    case (TRAFO_ID_MLINPYRAMID)
    
      ! Calculate with or without coordinates?
      if (.not. present(DpointsReal)) then
      
        !$omp parallel do private(ipt,DjacPrep) &
        !$omp if(nelements > TRAFO_NELEMMIN_OMP)
        do iel = 1,nelements
          ! Prepare the calculation of the Jacobi determinants
          call trafo_prepJac_pyra3D(Dcoords(:,:,iel), DjacPrep)
          
          ! Loop over the points
          do ipt=1,npointsPerEl
            ! Calculate the Jacobian matrix and determinant
            call trafo_calcJac_pyra3D (DjacPrep,Djac(:,ipt,iel),Ddetj(ipt,iel), &
                DpointsRef(1,ipt,iel),DpointsRef(2,ipt,iel),DpointsRef(3,ipt,iel))
          end do ! ipt
          
        end do !iel
        !$omp end parallel do
        
      else

        !$omp parallel do private(ipt,DjacPrep) &
        !$omp if(nelements > TRAFO_NELEMMIN_OMP)
        do iel = 1, nelements
          ! Prepare the calculation of the Jacobi determinants
          call trafo_prepJac_pyra3D(Dcoords(:,:,iel), DjacPrep)
          
          ! Loop over the points
          do ipt=1,npointsPerEl
            ! Calculate the Jacobian matrix and determinant
            call trafo_calcTrafo_pyra3d (DjacPrep,Djac(:,ipt,iel),Ddetj(ipt,iel), &
                DpointsRef(1,ipt,iel),DpointsRef(2,ipt,iel),DpointsRef(3,ipt,iel),&
                DpointsReal(1,ipt,iel),DpointsReal(2,ipt,iel),DpointsReal(3,ipt,iel))
          end do ! ipt
          
        end do ! iel
        !$omp end parallel do
        
      end if

    case (TRAFO_ID_MLINPRISM)
    
      ! Calculate with or without coordinates?
      if (.not. present(DpointsReal)) then
      
        !$omp parallel do private(ipt,DjacPrep) &
        !$omp if(nelements > TRAFO_NELEMMIN_OMP)
        do iel = 1,nelements
          ! Prepare the calculation of the Jacobi determinants
          call trafo_prepJac_prism3D(Dcoords(:,:,iel), DjacPrep)
          
          ! Loop over the points
          do ipt=1,npointsPerEl
            ! Calculate the Jacobian matrix and determinant
            call trafo_calcJac_prism3D (DjacPrep,Djac(:,ipt,iel),Ddetj(ipt,iel), &
                DpointsRef(1,ipt,iel),DpointsRef(2,ipt,iel),DpointsRef(3,ipt,iel))
          end do ! ipt
          
        end do !iel
        !$omp end parallel do
        
      else

        !$omp parallel do private(ipt,DjacPrep) &
        !$omp if(nelements > TRAFO_NELEMMIN_OMP)
        do iel = 1, nelements
          ! Prepare the calculation of the Jacobi determinants
          call trafo_prepJac_prism3D(Dcoords(:,:,iel), DjacPrep)
          
          ! Loop over the points
          do ipt=1,npointsPerEl
            ! Calculate the Jacobian matrix and determinant
            call trafo_calcTrafo_prism3d (DjacPrep,Djac(:,ipt,iel),Ddetj(ipt,iel), &
                DpointsRef(1,ipt,iel),DpointsRef(2,ipt,iel),DpointsRef(3,ipt,iel),&
                DpointsReal(1,ipt,iel),DpointsReal(2,ipt,iel),DpointsReal(3,ipt,iel))
          end do ! ipt
          
        end do ! iel
        !$omp end parallel do
        
      end if
    
    end select
        
  end select ! Dimension

  end subroutine
  
! **********************************************************************

!<subroutine>

  pure subroutine trafo_calctrafoabs (ctrafoType,Dcoords,&
                              DpointRef,Djac,ddetj,DpointReal)
!<description>
  ! General transformation.
  !
  ! The aim of this routine is to calculate the transformation between the
  ! reference element for one point. The element is given
  ! as a list of corner points in Dcoords.
  !
  ! For the point DpointsRef, the following information is calculated:
  ! 1.) Determinant of the mapping from the reference to the real element,
  ! 2.) the Jacobian matrix of the mapping,
  ! 3.) if the parameter DpointsReal is present: coordinates of the mapped
  !     points on the real element(s).
!</description>

!<input>
  ! ID of transformation to calculate
  integer(I32), intent(in) :: ctrafoType

  ! Coordinates of the corners of the element
  !  Dcoord(1,i) = x-coordinates of corner i on an element, 
  !  Dcoord(2,i) = y-coordinates of corner i on an element.
  ! DIMENSION(#space dimensions,NVE)
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Coordinates of the point on the reference element 
  !  DpointRef(1) = x-coordinates of point i on an element, 
  !  DpointRef(2) = y-coordinates of point i on an element.
  ! DIMENSION(#space dimension)
  real(DP), dimension(:), intent(in) :: DpointRef
!</input>

!<output>
  ! The Jacobian matrix of the mapping for each point.
  ! DIMENSION(number of entries in the matrix,npointsPerEl)
  real(DP), dimension(:), intent(out) :: Djac
  
  ! Absolute value of the Jacobian determinant of the mapping.
  ! DIMENSION(npointsPerEl)
  real(DP), intent(out) :: ddetj
  
  ! OPTIONAL: Array receiving the coordinates of the points in DpointRef,
  ! mapped from the reference element to the real element.
  ! If not specified, they are not computed.
  ! DIMENSION(#space dimension,npointsPerEl)
  real(DP), dimension(:), intent(out), optional :: DpointReal
!</output>

!</subroutine>

    call trafo_calctrafo (ctrafoType,Dcoords,&
                          DpointRef,Djac,ddetj,DpointReal)

    ! In 1D and 2D, the Jacobian determinant must always be positive.
    ! In 3D it can be negative.
  
    ddetj = abs(ddetj)

  end subroutine

! **********************************************************************

!<subroutine>

  subroutine trafo_calctrafoabs_mult (ctrafoType,npointsPerEl,Dcoords,&
                                      DpointsRef,Djac,Ddetj,DpointsReal)
!<description>
  ! General transformation support for multiple points on one element.
  !
  ! The aim of this routine is to calculate the transformation between the
  ! reference element for multiple points. The element is given
  ! as a list of corner points in Dcoords.
  !
  ! For every of these npointsPerEl points in the element specified 
  ! by DpointsRef, the following information is calculated:
  ! 1.) Determinant of the mapping from the reference to the real element,
  ! 2.) the Jacobian matrix of the mapping,
  ! 3.) if the parameter DpointsReal is present: coordinates of the mapped
  !     points on the real element(s).
!</description>

!<input>
  ! ID of transformation to calculate
  integer(I32), intent(in) :: ctrafoType

  ! Number of points in each element where to calculate the transformation
  integer, intent(in) :: npointsPerEl

  ! Coordinates of the corners of all the elements
  !  Dcoord(1,i) = x-coordinates of corner i on an element, 
  !  Dcoord(2,i) = y-coordinates of corner i on an element.
  ! DIMENSION(#space dimensions,NVE)
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Coordinates of the points on the reference element for each element 
  ! where to calculate the mapping.
  !  DpointsRef(1,i) = x-coordinates of point i on an element, 
  !  DpointsRef(2,i) = y-coordinates of point i on an element.
  ! ! DIMENSION(#space dimension,npointsPerEl)
  real(DP), dimension(:,:), intent(in) :: DpointsRef
!</input>

!<output>
  ! The Jacobian matrix of the mapping for each point.
  ! DIMENSION(number of entries in the matrix,npointsPerEl)
  real(DP), dimension(:,:), intent(out) :: Djac
  
  ! Absolute values of the Jacobian determinants of the mapping for all
  ! the points from the reference element to the real element.                                            
  ! DIMENSION(npointsPerEl)                                                           
  real(DP), dimension(:), intent(out) :: Ddetj                                        
  
  ! OPTIONAL: Array receiving the coordinates of the points in DpointsRef,
  ! mapped from the reference element to the real element.
  ! If not specified, they are not computed.
  ! DIMENSION(#space dimension,npointsPerEl)
  real(DP), dimension(:,:), intent(out), optional :: DpointsReal
!</output>

!</subroutine>

    call trafo_calctrafo_mult (ctrafoType,npointsPerEl,Dcoords,&
                               DpointsRef,Djac,Ddetj,DpointsReal)

    ! In 1D and 2D, the Jacobian determinant must always be positive.
    ! In 3D it can be negative.

    if (trafo_igetDimension(ctrafoType) .eq. NDIM3D) &
      Ddetj = abs(Ddetj)

  end subroutine

! **********************************************************************

!<subroutine>

  subroutine trafo_calctrafoabs_sim (ctrafoType,nelements,npointsPerEl,Dcoords,&
                                     DpointsRef,Djac,Ddetj,DpointsReal)

!<description>
  ! General transformation support for multiple points on multiple
  ! elements. 
  !
  ! The aim of this routine is to calculate the transformation between the
  ! reference element and multiple real elements. The elements are given
  ! as a list of corner points in Dcoords.
  !
  ! On every of the nelements elements given in this list, there are
  ! npointsPerEl points inside the element given in reference coordinates.
  ! For every of these npointsPerEl*nelements points, the following 
  ! information is calculated:
  ! 1.) Determinant of the mapping from the reference to the real element,
  ! 2.) the Jacobian matrix of the mapping,
  ! 3.) if the parameter DpointsReal is present: coordinates of the mapped
  !     points on the real element(s).
!</description>

!<input>
  ! ID of transformation to calculate
  integer(I32), intent(in) :: ctrafoType

  ! Number of elements where to calculate the transformation
  integer, intent(in) :: nelements
  
  ! Number of points in each element where to calculate the transformation
  integer, intent(in) :: npointsPerEl

  ! Coordinates of the corners of all the elements
  !  Dcoord(1,i,.) = x-coordinates of corner i on an element, 
  !  Dcoord(2,i,.) = y-coordinates of corner i on an element.
  ! DIMENSION(#space dimensions,NVE,nelements)
  real(DP), dimension(:,:,:), intent(in) :: Dcoords

  ! Coordinates of the points on the reference element for each element 
  ! where to calculate the mapping.
  ! DIMENSION(#space dimensions,npointsPerEl,nelements) for quadrilateral elements and
  ! DIMENSION(#space dimensions+1,npointsPerEl,nelements) for triangular elements.
  !
  ! For QUAD elements:
  !  DpointsRef(1,i,.) = x-coordinates of point i on an element, 
  !  DpointsRef(2,i,.) = y-coordinates of point i on an element.
  !
  ! For triangular elements:
  !  DpointsRef(1,i,.) = First barycentric coordinate of point i on an element
  !  DpointsRef(2,i,.) = Second barycentric coordinate of point i on an element
  !  DpointsRef(3,i,.) = Third barycentric coordinate of point i on an element
  real(DP), dimension(:,:,:), intent(in) :: DpointsRef
!</input>

!<output>
  ! The Jacobian matrix of the mapping for each point.
  ! DIMENSION(number of entries in the matrix,npointsPerEl,nelements)
  real(DP), dimension(:,:,:), intent(out) :: Djac
  
  ! Absolute values of the Jacobian determinants of the mapping for
  ! all the points from the reference element to the real element.
  ! DIMENSION(npointsPerEl,nelements)
  real(DP), dimension(:,:), intent(out) :: Ddetj
  
  ! OPTIONAL: Array receiving the coordinates of the points in DpointsRef,
  ! mapped from the reference element to the real element.
  ! If not specified, they are not computed.
  ! DIMENSION(#space dimensions,npointsPerEl,nelements)
  real(DP), dimension(:,:,:), intent(out), optional :: DpointsReal
!</output>

!</subroutine>

    call trafo_calctrafo_sim (ctrafoType,nelements,npointsPerEl,Dcoords,&
                              DpointsRef,Djac,Ddetj,DpointsReal)
            
    ! In 1D and 2D, the Jacobian determinant must always be positive.
    ! In 3D it can be negative.

    if (trafo_igetDimension(ctrafoType) .eq. NDIM3D) &
      Ddetj = abs(Ddetj)

  end subroutine
  

! ---------------------------------------------------------------------------
! Explaination of the quadrilateral transformation:
!
! We want to perform a transformation from the reference quadrilateral
! [-1,1]x[-1,1] onto a "real" quadrilateral:
!
!   (-1,1) ___ (1,1)         (x4,y4)   ___  (x3,y3)
!         |  |                        /  \
!         |__|           =>          /____\
!  (-1,-1)    (1,-1)         (x1,y1)        (x2,y2)   
!
! By theory this can by done with a bilinear mapping, i.e. a mapping
! Psi:R^2->R^2  of the form:
!
!      Psi (xi1) = ( a1 + a2*xi1 + a3*xi2 + a4*xi1*xi2 )
!          (xi2)   ( b1 + b2*xi1 + b3*xi2 + b4*xi1*xi2 )
!
! Our special transformation has to map:
!
!  Psi(-1,-1) = (x1,y1)
!  Psi( 1,-1) = (x2,y2)
!  Psi( 1, 1) = (x3,y3)
!  Psi(-1, 1) = (x4,y4)
!
! This gives the linear system:
!
!  a1 - a2 - a3 + a4 = x1       b1 - b2 - b3 + b4 = y1
!  a1 + a2 - a3 - a4 = x2       b1 + b2 - b3 - b4 = y2
!  a1 + a2 + a3 + a4 = x3       b1 + b2 + b3 + b4 = y3
!  a1 - a2 + a3 - a4 = x4       b1 - b2 + b3 - b4 = y4
!
! Reorder this to calculate the ai:
!
!  a1 = 1/4 * ( x1 + x2 + x3 + x4)      b1 = 1/4 * ( y1 + y2 + y3 + y4)
!  a2 = 1/4 * (-x1 + x2 + x3 - x4)      b2 = 1/4 * (-y1 + y2 + y3 - y4)
!  a3 = 1/4 * (-x1 - x2 + x3 + x4)      b3 = 1/4 * (-y1 - y2 + y3 + y4)
!  a4 = 1/4 * ( x1 - x2 + x3 - x4)      b4 = 1/4 * ( y1 - y2 + y3 - y4)
!
! The factors in the brackets in these equations are only dependent on 
! the corners of the quadrilateral, not of the current point. So they
! are constant for all points we want to map from the reference
! element to the real one. We call them here "auxiliary Jacobian factors".
! They can be calculated with trafo_calcJacPrepare in advance for all points
! that have to be mapped.
!
! The Jacobian matrix of the bilinear transformation is now 
! calculated as usual by partial differentiation. Thing above 
! coefficients ai and bi one can write:
!
! DPhi ( xi1 )  =  ( a2 + a4*xi2    a3 + a4*xi1 )
!      ( xi2 )     ( b2 + b4*xi2    b3 + b4*xi1 )
!
! which gives the Jacobian determinant of the 2x2-matrix:
!
!   det DPhi (xi1,xi2)  =  (DPhi[1,1]*DPhi[2,2] - DPhi[1,2]*DPhi[2,1]) (xi1,xi2)
!
! Using these information makes it possible to map a point (XI1,XI2)
! on the reference element to coordinates (XX,YY) on the real element by:
!
!   (XX,YY) := Psi(XI1,XI2) 
! ----------------------------------------------------------------------------

!<subroutine>

  pure subroutine trafo_prepJac_quad2D(Dcoords, DjacPrep)
  
!<description>
  ! Calculate auxiliary Jacobian factors for standard quadrilateral.
  ! This routine builds up constant factors that are later used during
  ! the transformation from the reference element to the real element.
!</description>

!<input>
  ! Coordinates of the four corners of the real quadrilateral.
  ! Dcoord(1,.) = x-coordinates, 
  ! Dcoord(2,.) = y-coordinates.
  ! DIMENSION(#space dimensions,NVE)
  real(DP), dimension(:,:), intent(in) :: Dcoords
!</input>
  
!<output>
  ! Arrays with auxiliary jacobian factors for later computation
  ! DIMENSION(TRAFO_NAUXJAC2D)
  real(DP), dimension(:), intent(out) :: DjacPrep
!</output>

!</subroutine>

  DjacPrep(1) = 0.25_DP*(Dcoords(1,1)+Dcoords(1,2)+Dcoords(1,3)+Dcoords(1,4))
  DjacPrep(2) = 0.25_DP*(-Dcoords(1,1)+Dcoords(1,2)+Dcoords(1,3)-Dcoords(1,4))
  DjacPrep(3) = 0.25_DP*(-Dcoords(1,1)-Dcoords(1,2)+Dcoords(1,3)+Dcoords(1,4))
  DjacPrep(4) = 0.25_DP*(Dcoords(1,1)-Dcoords(1,2)+Dcoords(1,3)-Dcoords(1,4))
  DjacPrep(5) = 0.25_DP*(Dcoords(2,1)+Dcoords(2,2)+Dcoords(2,3)+Dcoords(2,4))
  DjacPrep(6) = 0.25_DP*(-Dcoords(2,1)+Dcoords(2,2)+Dcoords(2,3)-Dcoords(2,4))
  DjacPrep(7) = 0.25_DP*(-Dcoords(2,1)-Dcoords(2,2)+Dcoords(2,3)+Dcoords(2,4))
  DjacPrep(8) = 0.25_DP*(Dcoords(2,1)-Dcoords(2,2)+Dcoords(2,3)-Dcoords(2,4))

  end subroutine 

  !************************************************************************

!<subroutine>

  pure subroutine trafo_calcTrafo_quad2d (DjacPrep,Djac,ddetj, &
                                   dparx,dpary,dxreal,dyreal)

!<description>
  ! This subroutine performs two tasks:
  ! -Initialisation of a a given 2x2 matrix with the
  !  mapping information from the reference element to the "real"
  !  quadrilateral. Calculation of the Jacobian determinant
  ! -Transformation of a given point on the reference element onto
  !  the "real" element
  ! Both things are performed simultaneously because the jacobian
  ! determinant is dependent of the point.
  !
  ! Before this routine can be called, the auxiliary factors DjacPrep
  ! have to be calculated with trafo_calcJacPrepare for the 
  ! considered element.
!</description>  

!<input>
  ! Auxiliary constants for the considered element with
  ! coordinates in Dcoord; have to be computed previously
  ! by trafo_calcJacPrepare2D.
  ! DIMENSION(TRAFO_NAUXJAC2D)
  real(DP), dimension(:), intent(in) :: DjacPrep
  
  ! Coordinates of a point on the reference element
  real(DP), intent(in) :: dparx,dpary
  
!</input>
  
!<output>
  ! The Jacobian matrix of the mapping from the reference to the 
  ! real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  real(DP), dimension(:), intent(out) :: Djac
  
  ! The determinant of the mapping.
  real(DP), intent(out) :: ddetj
  
  ! X-Coordinates of a point on the real element
  real(DP), intent(out) :: dxreal
  
  ! X-Coordinates of a point on the real element
  real(DP), intent(out) :: dyreal
!</output>
  
!</subroutine>

  ! Jacobian matrix of the mapping
  Djac(1) = DjacPrep(2) + DjacPrep(4)*dpary
  Djac(2) = DjacPrep(6) + DjacPrep(8)*dpary
  Djac(3) = DjacPrep(3) + DjacPrep(4)*dparx
  Djac(4) = DjacPrep(7) + DjacPrep(8)*dparx

  ! Determinant of the mapping
  ddetj = Djac(1)*Djac(4) - Djac(3)*Djac(2)
  
  ! Map the point to the real element
  dxreal = DjacPrep(1) + DjacPrep(2)*dparx + DjacPrep(3)*dpary &
         + DjacPrep(4)*dparx*dpary
  dyreal = DjacPrep(5) + DjacPrep(6)*dparx + DjacPrep(7)*dpary &
         + DjacPrep(8)*dparx*dpary
    
  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine trafo_calcJac_quad2D (DjacPrep,Djac,ddetj,dparx,dpary)

!<description>
  ! Calculate Jacobian determinant of mapping from reference- to
  ! real element.
  !
  ! This routine only calculates the Jacobian determinant of the
  ! mapping. This is on contrast to trafo_calcRealCoords, which not only 
  ! calculates this determinant but also maps the point. So this routine
  ! can be used to speed up the code if the coordinates of the
  ! mapped point already exist.
  !
  ! Before this routine can be called, the auxiliary factors DJF
  ! have to be calculated with trafo_calcJacPrepare for the considered element.
!</description>  

!<input>
  ! Auxiliary constants for the considered element with
  ! coordinates in Dcoord; have to be computed previously
  ! by trafo_calcJacPrepare2D.
  ! DIMENSION(TRAFO_NAUXJAC2D)
  real(DP), dimension(:), intent(in) :: DjacPrep
  
  ! Coordinates of a point on the reference element
  real(DP), intent(in) :: dparx,dpary
!</input>
  
!<output>
  ! The Jacobian matrix of the mapping from the reference to the 
  ! real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  real(DP), dimension(:), intent(out) :: Djac
  
  ! The determinant of the mapping.
  real(DP), intent(out) :: ddetj
!</output>
  
!</subroutine>

  ! Jacobian matrix of the mapping
  Djac(1) = DjacPrep(2) + DjacPrep(4)*dpary
  Djac(2) = DjacPrep(6) + DjacPrep(8)*dpary
  Djac(3) = DjacPrep(3) + DjacPrep(4)*dparx
  Djac(4) = DjacPrep(7) + DjacPrep(8)*dparx

  ! Determinant of the mapping
  ddetj = Djac(1)*Djac(4) - Djac(3)*Djac(2)
  
  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine trafo_calcRealCoords2D (DjacPrep,dparx,dpary,dxreal,dyreal)

!<description>
  ! This subroutine computes the real coordinates of a point which 
  ! is given by parameter values (dparx,dpary) on the reference element. 
  ! It is assumed that the array DjacPrep was initialised before using
  ! with trafo_calcJacPrepare.
!</description>  

!<input>
  ! Auxiliary constants for the considered element with
  ! coordinates in Dcoord; have to be computed previously
  ! by trafo_calcJacPrepare2D.
  ! DIMENSION(TRAFO_NAUXJAC2D)
  real(DP), dimension(:), intent(in) :: DjacPrep
  
  ! Coordinates of a point on the reference element
  real(DP), intent(in) :: dparx,dpary
!</input>
  
!<output>
  ! X-Coordinates of a point on the real element
  real(DP), intent(out) :: dxreal
  
  ! X-Coordinates of a point on the real element
  real(DP), intent(out) :: dyreal
!</output>
  
!</subroutine>

  ! Map the point to the real element
  dxreal = DjacPrep(1) + DjacPrep(2)*dparx + DjacPrep(3)*dpary &
         + DjacPrep(4)*dparx*dpary
  dyreal = DjacPrep(5) + DjacPrep(6)*dparx + DjacPrep(7)*dpary &
         + DjacPrep(8)*dparx*dpary
  
  end subroutine

  !****************************************************************************

!<subroutine>

  pure subroutine trafo_mapCubPts1Dto2DRefQuad(iedge, ncubp, Dxi1D, Dxi2D)

!<description>
  ! This routine maps the coordinates of the cubature points in 1D
  ! (given by Dxi1D(1..ncubp,1)) to an edge on the reference quadrilateral
  ! in 2D. iedge specifies the edge where to map the cubature points to.
!</description>

!<input>
  ! Number of the local edge of the element where to map the points to.
  integer, intent(in) :: iedge

  ! number of cubature points
  integer , intent(in) :: ncubp
  
  ! Cubature point coordinates on 1D reference interval [-1,1]
  !     Dxi(1..ncubp,1)=coordinates
  real(DP), dimension(:,:), intent(in) :: Dxi1D
!</input>
  
!<output>
  ! Coordinates of the cubature points on the edge in 2D.
  !        Dxi2D(1..ncubp,1)=x-coord, 
  !        Dxi2D(1..ncubp,2)=y-coord
  real(DP), dimension(:,:), intent(out) :: Dxi2D
!</output>

!</subroutine>    

  ! local variables
  integer :: ii

    select case(iedge)
    case (1)
      ! Edge 1 is on the bottom of the reference element
      do ii = 1,ncubp
        Dxi2D(ii,1) = Dxi1D(ii,1)
        Dxi2D(ii,2) = -1.0_DP
      end do
      
    case(2)
      ! Edge 2 is on the right of the reference element
      do ii = 1,ncubp
        Dxi2D(ii,1) = 1.0_DP
        Dxi2D(ii,2) = Dxi1D(ii,1)
      end do
      
    case(3)
      ! Edge 3 is on the top of the reference element
      do ii = 1,ncubp
        Dxi2D(ii,1) = -Dxi1D(ii,1)
        Dxi2D(ii,2) = 1.0_DP
      end do
      
    case(4)
      ! Edge 4 is on the left of the reference element
      do ii = 1,ncubp
        Dxi2D(ii,1) = -1.0_DP
        Dxi2D(ii,2) = -Dxi1D(ii,1)
      end do
      
    end select
  
  end subroutine

  !****************************************************************************

!<subroutine>

  pure subroutine trafo_mapCubPts1Dto2DTriBary(iedge, ncubp, Dxi1D, Dxi2D)

!<description>
  ! This routine maps the coordinates of the cubature points in 1D
  ! (given by Dxi1D(1..ncubp,1)) to an edge on a 2D triangle
  ! in barycentric coordinates. 
  ! iedge specifies the edge where to map the cubature points to.
!</description>

!<input>
  ! Number of the local edge of the element where to map the points to.
  integer, intent(in) :: iedge

  ! number of cubature points
  integer , intent(in) :: ncubp
  
  ! Cubature point coordinates on 1D reference interval [-1,1]
  !     Dxi(1..ncubp,1)=coordinates
  real(DP), dimension(:,:), intent(in) :: Dxi1D
!</input>
  
!<output>
  ! Coordinates of the cubature points on the edge in 2D.
  !        Dxi2D(1..ncubp,1)=1st barycentric coordinate, 
  !        Dxi2D(1..ncubp,2)=2nd barycentric coordinate,
  !        Dxi2D(1..ncubp,3)=3rd barycentric coordinate
  real(DP), dimension(:,:), intent(out) :: Dxi2D
!</output>

!</subroutine>    

  ! local variables
  integer :: ii

    select case(iedge)
    case(1)
      ! Edge 1 is between 1st and 2nd point
      do ii = 1,ncubp
        Dxi2D(ii,1) = 0.5_DP*(1.0_DP-Dxi1D(ii,1))
        Dxi2D(ii,2) = 0.5_DP*(1.0_DP+Dxi1D(ii,1))
        Dxi2D(ii,3) = 0.0_DP
      end do
      
    case(2)
      ! Edge 2 is between 2nd and 3rd point
      do ii = 1,ncubp
        Dxi2D(ii,1) = 0.0_DP
        Dxi2D(ii,2) = 0.5_DP*(1.0_DP-Dxi1D(ii,1))
        Dxi2D(ii,3) = 0.5_DP*(1.0_DP+Dxi1D(ii,1))
      end do
    
    case(3)
      ! Edge 3 is between 3rd and 1st point
      do ii = 1,ncubp
        Dxi2D(ii,1) = 0.5_DP*(1.0_DP+Dxi1D(ii,1))
        Dxi2D(ii,2) = 0.0_DP
        Dxi2D(ii,3) = 0.5_DP*(1.0_DP-Dxi1D(ii,1))
      end do
      
    end select
  
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine trafo_mapCubPts1Dto2D(icoordSystem, iedge, ncubp, Dxi1D, Dxi2D)

!<description>
  ! This routine maps the coordinates of the cubature points 
  ! (given by Dxi1D(1..ncubp,1)) on the 1D reference interval [-1,1] 
  ! to an edge on the 2D reference element.
  ! iedge specifies the edge where to map the cubature points to.
  !
  ! ctrafoType specifies the type of transformation that is necessary
  ! to transform coordinates from the 2D reference element to the 2D
  ! real element. This constant also specifies the type of reference element
  ! and the coordinate system there.
!</description>

!<input>
  ! Coordinate system identifier. One of the TRAFO_CS_xxxx constants. Defines
  ! the type of the coordinate system that is used for specifying the coordinates
  ! on the reference element.
  integer(I32), intent(in) :: icoordSystem
  
  ! Number of the local edge of the element where to map the points to.
  integer, intent(in) :: iedge

  ! number of cubature points
  integer , intent(in) :: ncubp
  
  ! Cubature point coordinates on 1D reference interval [-1,1]
  !     Dxi(1..ncubp,1)=coordinates
  real(DP), dimension(:,:), intent(in) :: Dxi1D
!</input>
  
!<output>
  ! Coordinates of the cubature points on the edge in 2D.
  ! For 2D triangles with barycentric coordinates:
  !        Dxi2D(1..ncubp,1)=1st barycentric coordinate, 
  !        Dxi2D(1..ncubp,2)=2nd barycentric coordinate,
  !        Dxi2D(1..ncubp,3)=3rd barycentric coordinate
  ! For 2D quads with standard coordinates:
  !        Dxi2D(1..ncubp,1)=x-coord, 
  !        Dxi2D(1..ncubp,2)=y-coord.
  real(DP), dimension(:,:), intent(out) :: Dxi2D
!</output>

!</subroutine>    

  ! What type of transformation do we have? First decide on the dimension,
  ! then on the actual ID.
  select case (icoordSystem)
  
  case (TRAFO_CS_BARY2DTRI)
    ! Triangle, barycentric coordinates
    call trafo_mapCubPts1Dto2DTriBary(iedge, ncubp, Dxi1D, Dxi2D)
  
  case (TRAFO_CS_REF2DQUAD,TRAFO_CS_REAL2DQUAD)
    ! Quadrilateral, [-1,1]^2
    call trafo_mapCubPts1Dto2DRefQuad(iedge, ncubp, Dxi1D, Dxi2D)
  
  case DEFAULT
    call output_line ('Unsupported coordinate system.', &
                      OU_CLASS_ERROR,OU_MODE_STD,'trafo_mapCubPts1Dto2D')  
  end select    

  end subroutine

  !****************************************************************************

!<subroutine>

  pure subroutine trafo_mapCubPts2Dto3DRefHexa(iface, ncubp, Dxi2D, Dxi3D)

!<description>
  ! This routine maps the coordinates of the cubature points in 2D
  ! (given by Dxi2D(1..ncubp,1..2)) to a face on the reference hexahedron
  ! in 3D. iface specifies the face where to map the cubature points to.
!</description>

!<input>
  ! Number of the local face of the element where to map the points to.
  integer, intent(in) :: iface

  ! number of cubature points
  integer , intent(in) :: ncubp
  
  ! Cubature point coordinates on 2D reference quadrilateral [-1,1]x[-1,1]
  ! Dxi2D(1..ncubp,1)= X-coordinates,
  ! Dxi2D(1..ncubp,2)= Y-coordinates
  real(DP), dimension(:,:), intent(in) :: Dxi2D
!</input>
  
!<output>
  ! Coordinates of the cubature points on the face in 3D.
  ! Dxi3D(1..ncubp,1)= X-coordinates,
  ! Dxi3D(1..ncubp,2)= Y-coordinates,
  ! Dxi3D(1..ncubp,3)= Z-coordinates
  real(DP), dimension(:,:), intent(out) :: Dxi3D
!</output>

!</subroutine>    

  ! local variables
  integer :: i

    select case(iface)
    case (1)
      ! Bottom face => Z = -1
      do i=1, ncubp
        Dxi3D(i,1) = Dxi2D(i,1)
        Dxi3D(i,2) = Dxi2D(i,2)
        Dxi3D(i,3) = -1.0_DP
      end do
      
    case (2)
      ! Front face => Y = -1
      do i=1, ncubp
        Dxi3D(i,1) = Dxi2D(i,1)
        Dxi3D(i,2) = -1.0_DP
        Dxi3D(i,3) = Dxi2D(i,2)
      end do
      
    case (3)
      ! Right face => X = 1
      do i=1, ncubp
        Dxi3D(i,1) = 1.0_DP
        Dxi3D(i,2) = Dxi2D(i,1)
        Dxi3D(i,3) = Dxi2D(i,2)
      end do
      
    case (4)
      ! Back face => Y = 1
      do i=1, ncubp
        Dxi3D(i,1) = -Dxi2D(i,1)
        Dxi3D(i,2) = 1.0_DP
        Dxi3D(i,3) = Dxi2D(i,2)
      end do
      
    case (5)
      ! Left face => X = -1
      do i=1, ncubp
        Dxi3D(i,1) = -1.0_DP
        Dxi3D(i,2) = -Dxi2D(i,1)
        Dxi3D(i,3) = Dxi2D(i,2)
      end do
      
    case (6)
      ! Top face => Z = 1
      do i=1, ncubp
        Dxi3D(i,1) = Dxi2D(i,1)
        Dxi3D(i,2) = Dxi2D(i,2)
        Dxi3D(i,3) = 1.0_DP
      end do
      
    end select
  
  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine trafo_calcBackTrafo_quad2d (Dcoord, &
     dxreal,dyreal,dparx,dpary)

!<description>
  ! This subroutine is to find the coordinates on the reference
  ! element (dparx,dpary) for a given point (dxreal,dyreal) 
  ! in real coordinates.
  !
  ! Remark: This is a difficult task, as usually in FEM codes the
  ! parameter values (dparx,dpary) are known and one wants to obtain 
  ! the real coordinates.
  !
  ! Inverting the bilinear trafo in a straightforward manner by using 
  ! pq-formula does not work very well, as it is numerically unstable.
  ! For parallelogram-shaped elements, one would have to introduce a
  ! special treatment. For nearly parallelogram-shaped elements, this 
  ! can cause a crash as the argument of the square root can become 
  ! negative due to rounding errors. In the case of points near
  ! the element borders, we divide nearly 0/0.
  !
  ! Therefore, we have implemented the algorithm described in
  !
  ! [Introduction to Finite Element Methods, Carlos Felippa,
  !  Department of Aerospace Engineering Sciences and Center for 
  !  Aerospace Structures, http://titan.colorado.edu/courses.d/IFEM.d/]
  !
  ! Note: For the transformation of coordinates to the barycentric
  ! coordinate system for triangles, one can use 
  ! gaux_getBarycentricCoords_tri2D.
!</description>  

!<input>
  ! Coordinates of the four points forming the element.
  ! Dcoord(1,.) = x-coordinates,
  ! Dcoord(2,.) = y-coordinates.
  ! DIMENSION(#space dimensions,NVE)
  real(DP), dimension(:,:), intent(in) :: Dcoord 
  
  ! Coordinates of a point on the real element
  real(DP), intent(in) :: dxreal,dyreal
!</input>
  
!<output>
  ! Coordinates of a point on the reference element
  real(DP), intent(out) :: dparx,dpary
!</output>
  
!</subroutine>

    ! local variables

    real(DP) :: X1,X2,X3,X4,Y1,Y2,Y3,Y4,XB,YB,XCX,YCX,XCE,YCE
    real(DP) :: A,J1,J2,X0,Y0,BXI,BETA,CXI,XP0,YP0,CETA
    real(DP) :: ROOT1,ROOT2,XIP1,XIP2,ETAP1,ETAP2,D1,D2

    ! Get nodal x-coordinates
    X1 = Dcoord(1,1)
    X2 = Dcoord(1,2)
    X3 = Dcoord(1,3)
    X4 = Dcoord(1,4)

    ! Get nodal y-coordinates
    Y1 = DCOORD(2,1)
    Y2 = DCOORD(2,2)
    Y3 = DCOORD(2,3)
    Y4 = DCOORD(2,4)

    XB = X1-X2+X3-X4
    YB = Y1-Y2+Y3-Y4

    XCX = X1+X2-X3-X4
    YCX = Y1+Y2-Y3-Y4

    XCE = X1-X2-X3+X4
    YCE = Y1-Y2-Y3+Y4

    A = 0.5*((X3-X1)*(Y4-Y2)-(X4-X2)*(Y3-Y1))

    J1 = (X3-X4)*(Y1-Y2)-(X1-X2)*(Y3-Y4)
    J2 = (X2-X3)*(Y1-Y4)-(X1-X4)*(Y2-Y3)

    X0 = 0.25*(X1+X2+X3+X4)
    Y0 = 0.25*(Y1+Y2+Y3+Y4)

    XP0 = dxreal-X0
    YP0 = dyreal-Y0

    BXI  =  A-XP0*YB+YP0*XB
    BETA = -A-XP0*YB+YP0*XB

    CXI  = XP0*YCX-YP0*XCX
    CETA = XP0*YCE-YP0*XCE

    ROOT1 = -sqrt(BXI**2-2.0*J1*CXI)-BXI
    ROOT2 =  sqrt(BXI**2-2.0*J1*CXI)-BXI
    if (ROOT1 .ne. 0.0_DP) then
      XIP1 = 2.0*CXI/ROOT1
    else
      XIP1 = 1E15_DP
    end if
    if (ROOT2 .ne. 0.0_DP) then
      XIP2 = 2.0_DP*CXI/ROOT2
    else
      XIP2 = 1E15_DP
    end if

    ROOT1 =  sqrt(BETA**2+2.0_DP*J2*CETA)-BETA
    ROOT2 = -sqrt(BETA**2+2.0_DP*J2*CETA)-BETA
    if (ROOT1 .ne. 0.0_DP) then
      ETAP1 = 2.0_DP*CETA/ROOT1
    else
      ETAP1 = 1E15_DP
    end if
    if (ROOT2 .ne. 0.0_DP) then
      ETAP2 = 2.0_DP*CETA/ROOT2
    else
      ETAP2 = 1E15_DP
    end if

    !D1 = SQRT(XIP1**2+ETAP1**2)
    !D2 = SQRT(XIP2**2+ETAP2**2)
    D1 = XIP1**2+ETAP1**2
    D2 = XIP2**2+ETAP2**2

    if (D1 .lt. D2) then
      dparx = XIP1
      dpary = ETAP1
    else
      dparx = XIP2
      dpary = ETAP2
    end if
  
  end subroutine

! ---------------------------------------------------------------------------
! Explaination of the hexahedral transformation:
!
! We want to perform a transformation from the reference hexahedron
! [-1,1]x[-1,1]x[-1,1] onto a "real" hexahedron:
!
!             H-------G              W---V
!            /|      /|             /|   |\
!           E-------F |            T-----U \
!           | |     | |     =>     | |    \ \
!           | D-----|-C            | S-----\-R
!           |/      |/             |/       \|
!           A-------B              P---------Q
!
! where:   x   y   z
!    A = (-1, -1, -1)         P = (x1, y1, z1)
!    B = ( 1, -1, -1)         Q = (x2, y2, z2)
!    C = ( 1,  1, -1)         R = (x3, y3, z3)
!    D = (-1,  1, -1)         S = (x4, y4, z4)
!    E = (-1, -1,  1)         T = (x5, y5, z5)
!    F = ( 1, -1,  1)         U = (x6, y6, z6)
!    G = ( 1,  1,  1)         V = (x7, y7, z7)
!    H = (-1,  1,  1)         W = (x8, y8, z8)
!
! By theory this can be done with a trilinear mapping, i.e. mapping
! Psi:R^3 -> R^3 of the form:
!
!    / yi1 \       / xi1 \
!    | yi2 | = Psi | xi2 |
!    \ yi3 /       \ xi3 /
!
! where:
!
!  yi1 = a1 + a2*xi1 + a3*xi2 + a4*xi3 + a5*xi1*xi2 
!      + a6*xi1*xi3 + a7*xi2*xi3 + a8*xi1*xi2*xi3
!  yi2 = b1 + b2*xi1 + b3*xi2 + b4*xi3 + b5*xi1*xi2
!      + b6*xi1*xi3 + b7*xi2*xi3 + b8*xi1*xi2*xi3
!  yi3 = c1 + c2*xi1 + c3*xi2 + c4*xi3 + c5*xi1*xi2
!      + c6*xi1*xi3 + c7*xi2*xi3 + c8*xi1*xi2*xi3
!
! Our special transformation has to map:
!
!  Psi(-1,-1,-1) = (x1, y1, z1)
!  Psi( 1,-1,-1) = (x2, y2, z2)
!  Psi( 1, 1,-1) = (x3, y3, z3)
!  Psi(-1, 1,-1) = (x4, y4, z4)
!  Psi(-1,-1, 1) = (x5, y5, z5)
!  Psi( 1,-1, 1) = (x6, y6, z6)
!  Psi( 1, 1, 1) = (x7, y7, z7)
!  Psi(-1, 1, 1) = (x8, y8, z8)
!
! This gives us 3 linear systems:
!
!  a1 - a2 - a3 - a4 + a5 + a6 + a7 - a8 = x1
!  a1 + a2 - a3 - a4 - a5 - a6 + a7 + a8 = x2
!  a1 + a2 + a3 - a4 + a5 - a6 - a7 - a8 = x3
!  a1 - a2 + a3 - a4 - a5 + a6 - a7 + a8 = x4
!  a1 - a2 - a3 + a4 + a5 - a6 - a7 + a8 = x5
!  a1 + a2 - a3 + a4 - a5 + a6 - a7 - a8 = x6
!  a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 = x7
!  a1 - a2 + a3 + a4 - a5 - a6 + a7 - a8 = x8
!
! (For the two other linear systems replace 'a' with 'b'/'c'
!  and 'x' with 'y'/'z' ...)
!
! Reorder the systems to calculate the ai, bi and ci:
!
!  a1 = 1/8 * ( x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8)
!  a2 = 1/8 * (-x1 + x2 + x3 - x4 - x5 + x6 + x7 - x8)
!  a3 = 1/8 * (-x1 - x2 + x3 + x4 - x5 - x6 + x7 + x8)
!  a4 = 1/8 * (-x1 - x2 - x3 - x4 + x5 + x6 + x7 + x8)
!  a5 = 1/8 * ( x1 - x2 + x3 - x4 + x5 - x6 + x7 - x8)
!  a6 = 1/8 * ( x1 - x2 - x3 + x4 - x5 + x6 + x7 - x8)
!  a7 = 1/8 * ( x1 + x2 - x3 - x4 - x5 - x6 + x7 + x8)
!  a8 = 1/8 * (-x1 + x2 - x3 + x4 + x5 - x6 + x7 - x8)
!
! The ai, bi and ci depend only on the corners of the hexahedron, but not on
! the current point. So they are constant for all points we want to map from
! the reference element to the real one. We call them `auxiliary Jacobian
! factors`. They are calculated with trafo_calcJacPrepare3D in advance for
! all points that have to be mapped.
! The Jacobian matrix of the trilinear transformation is now calculated
! as usual by partial differentiation:
!
! DPhi[1,1](xi1,xi2,xi3) = a2 + a5*xi2 + a6*xi3 + a8*xi2*xi3
! DPhi[1,2](xi1,xi2,xi3) = a3 + a5*xi1 + a7*xi3 + a8*xi1*xi3
! DPhi[1,3](xi1,xi2,xi3) = a4 + a6*xi1 + a7*xi2 + a8*xi1*xi2
! DPhi[1,1](xi1,xi2,xi3) = b2 + b5*xi2 + b6*xi3 + b8*xi2*xi3
! DPhi[1,2](xi1,xi2,xi3) = b3 + b5*xi1 + b7*xi3 + b8*xi1*xi3
! DPhi[1,3](xi1,xi2,xi3) = b4 + b6*xi1 + b7*xi2 + b8*xi1*xi2
! DPhi[1,1](xi1,xi2,xi3) = c2 + c5*xi2 + c6*xi3 + c8*xi2*xi3
! DPhi[1,2](xi1,xi2,xi3) = c3 + c5*xi1 + c7*xi3 + c8*xi1*xi3
! DPhi[1,3](xi1,xi2,xi3) = c4 + c6*xi1 + c7*xi2 + c8*xi1*xi2
!
! which gives the Jacobian determinant of the 3x3-matrix:
!
! det DPhi = DPhi[1,1]*(DPhi[2,2]*DPhi[3,3] - DPhi[2,3]*DPhi[3,2])
!          + DPhi[1,2]*(DPhi[2,3]*DPhi[3,1] - DPhi[2,1]*DPhi[3,3])
!          + DPhi[1,3]*(DPhi[2,1]*DPhi[3,2] - DPhi[2,2]*DPhi[3,1])
! ----------------------------------------------------------------------------

!<subroutine>

  pure subroutine trafo_prepJac_hexa3D(Dcoords, DjacPrep)
  
!<description>
  ! Calculate auxiliary Jacobian factors for standard hexahedron.
  ! This routine builds up constant factors that are later used during
  ! the transformation from the reference element to the real element.
!</description>

!<input>
  ! Coordinates of the four corners of the real hexahedron.
  ! Dcoord(1,.) = x-coordinates, 
  ! Dcoord(2,.) = y-coordinates,
  ! Dcoord(3,.) = z-coordinates.
  ! DIMENSION(#space dimensions,NVE)
  real(DP), dimension(:,:), intent(in) :: Dcoords
!</input>
  
!<output>
  ! Arrays with auxiliary jacobian factors for later computation
  ! DIMENSION(TRAFO_NAUXJAC3D)
  real(DP), dimension(:), intent(out) :: DjacPrep
!</output>

!</subroutine>

    DjacPrep( 1) = 0.125_DP *&
                 ( Dcoords(1,1)+Dcoords(1,2)+Dcoords(1,3)+Dcoords(1,4)&
                  +Dcoords(1,5)+Dcoords(1,6)+Dcoords(1,7)+Dcoords(1,8))
    DjacPrep( 2) = 0.125_DP *&
                 (-Dcoords(1,1)+Dcoords(1,2)+Dcoords(1,3)-Dcoords(1,4)&
                  -Dcoords(1,5)+Dcoords(1,6)+Dcoords(1,7)-Dcoords(1,8))
    DjacPrep( 3) = 0.125_DP *&
                 (-Dcoords(1,1)-Dcoords(1,2)+Dcoords(1,3)+Dcoords(1,4)&
                  -Dcoords(1,5)-Dcoords(1,6)+Dcoords(1,7)+Dcoords(1,8))
    DjacPrep( 4) = 0.125_DP *&
                 (-Dcoords(1,1)-Dcoords(1,2)-Dcoords(1,3)-Dcoords(1,4)&
                  +Dcoords(1,5)+Dcoords(1,6)+Dcoords(1,7)+Dcoords(1,8))
    DjacPrep( 5) = 0.125_DP *&
                 ( Dcoords(1,1)-Dcoords(1,2)+Dcoords(1,3)-Dcoords(1,4)&
                  +Dcoords(1,5)-Dcoords(1,6)+Dcoords(1,7)-Dcoords(1,8))
    DjacPrep( 6) = 0.125_DP *&
                 ( Dcoords(1,1)-Dcoords(1,2)-Dcoords(1,3)+Dcoords(1,4)&
                  -Dcoords(1,5)+Dcoords(1,6)+Dcoords(1,7)-Dcoords(1,8))
    DjacPrep( 7) = 0.125_DP *&
                 ( Dcoords(1,1)+Dcoords(1,2)-Dcoords(1,3)-Dcoords(1,4)&
                  -Dcoords(1,5)-Dcoords(1,6)+Dcoords(1,7)+Dcoords(1,8))
    DjacPrep( 8) = 0.125_DP *&
                 (-Dcoords(1,1)+Dcoords(1,2)-Dcoords(1,3)+Dcoords(1,4)&
                  +Dcoords(1,5)-Dcoords(1,6)+Dcoords(1,7)-Dcoords(1,8))
    DjacPrep( 9) = 0.125_DP *&
                 ( Dcoords(2,1)+Dcoords(2,2)+Dcoords(2,3)+Dcoords(2,4)&
                  +Dcoords(2,5)+Dcoords(2,6)+Dcoords(2,7)+Dcoords(2,8))
    DjacPrep(10) = 0.125_DP *&
                 (-Dcoords(2,1)+Dcoords(2,2)+Dcoords(2,3)-Dcoords(2,4)&
                  -Dcoords(2,5)+Dcoords(2,6)+Dcoords(2,7)-Dcoords(2,8))
    DjacPrep(11) = 0.125_DP *&
                 (-Dcoords(2,1)-Dcoords(2,2)+Dcoords(2,3)+Dcoords(2,4)&
                  -Dcoords(2,5)-Dcoords(2,6)+Dcoords(2,7)+Dcoords(2,8))
    DjacPrep(12) = 0.125_DP *&
                 (-Dcoords(2,1)-Dcoords(2,2)-Dcoords(2,3)-Dcoords(2,4)&
                  +Dcoords(2,5)+Dcoords(2,6)+Dcoords(2,7)+Dcoords(2,8))
    DjacPrep(13) = 0.125_DP *&
                 ( Dcoords(2,1)-Dcoords(2,2)+Dcoords(2,3)-Dcoords(2,4)&
                  +Dcoords(2,5)-Dcoords(2,6)+Dcoords(2,7)-Dcoords(2,8))
    DjacPrep(14) = 0.125_DP *&
                 ( Dcoords(2,1)-Dcoords(2,2)-Dcoords(2,3)+Dcoords(2,4)&
                  -Dcoords(2,5)+Dcoords(2,6)+Dcoords(2,7)-Dcoords(2,8))
    DjacPrep(15) = 0.125_DP *&
                 ( Dcoords(2,1)+Dcoords(2,2)-Dcoords(2,3)-Dcoords(2,4)&
                  -Dcoords(2,5)-Dcoords(2,6)+Dcoords(2,7)+Dcoords(2,8))
    DjacPrep(16) = 0.125_DP *&
                 (-Dcoords(2,1)+Dcoords(2,2)-Dcoords(2,3)+Dcoords(2,4)&
                  +Dcoords(2,5)-Dcoords(2,6)+Dcoords(2,7)-Dcoords(2,8))
    DjacPrep(17) = 0.125_DP *&
                 ( Dcoords(3,1)+Dcoords(3,2)+Dcoords(3,3)+Dcoords(3,4)&
                  +Dcoords(3,5)+Dcoords(3,6)+Dcoords(3,7)+Dcoords(3,8))
    DjacPrep(18) = 0.125_DP *&
                 (-Dcoords(3,1)+Dcoords(3,2)+Dcoords(3,3)-Dcoords(3,4)&
                  -Dcoords(3,5)+Dcoords(3,6)+Dcoords(3,7)-Dcoords(3,8))
    DjacPrep(19) = 0.125_DP *&
                 (-Dcoords(3,1)-Dcoords(3,2)+Dcoords(3,3)+Dcoords(3,4)&
                  -Dcoords(3,5)-Dcoords(3,6)+Dcoords(3,7)+Dcoords(3,8))
    DjacPrep(20) = 0.125_DP *&
                 (-Dcoords(3,1)-Dcoords(3,2)-Dcoords(3,3)-Dcoords(3,4)&
                  +Dcoords(3,5)+Dcoords(3,6)+Dcoords(3,7)+Dcoords(3,8))
    DjacPrep(21) = 0.125_DP *&
                 ( Dcoords(3,1)-Dcoords(3,2)+Dcoords(3,3)-Dcoords(3,4)&
                  +Dcoords(3,5)-Dcoords(3,6)+Dcoords(3,7)-Dcoords(3,8))
    DjacPrep(22) = 0.125_DP *&
                 ( Dcoords(3,1)-Dcoords(3,2)-Dcoords(3,3)+Dcoords(3,4)&
                  -Dcoords(3,5)+Dcoords(3,6)+Dcoords(3,7)-Dcoords(3,8))
    DjacPrep(23) = 0.125_DP *&
                 ( Dcoords(3,1)+Dcoords(3,2)-Dcoords(3,3)-Dcoords(3,4)&
                  -Dcoords(3,5)-Dcoords(3,6)+Dcoords(3,7)+Dcoords(3,8))
    DjacPrep(24) = 0.125_DP *&
                 (-Dcoords(3,1)+Dcoords(3,2)-Dcoords(3,3)+Dcoords(3,4)&
                  +Dcoords(3,5)-Dcoords(3,6)+Dcoords(3,7)-Dcoords(3,8))
  
  end subroutine 

  !************************************************************************

!<subroutine>

  pure subroutine trafo_calcTrafo_hexa3d (DjacPrep,Djac,ddetj, &
                                          dparx,dpary,dparz,&
                                          dxreal,dyreal,dzreal)

!<description>
  ! This subroutine performs two tasks:
  ! -Initialisation of a a given 3x3 matrix with the
  !  mapping information from the reference element to the "real"
  !  hexahedron. Calculation of the Jacobian determinant
  ! -Transformation of a given point on the reference element onto
  !  the "real" element
  ! Both things are performed simultaneously because the jacobian
  ! determinant is dependent of the point.
  !
  ! Before this routine can be called, the auxiliary factors DjacPrep
  ! have to be calculated with trafo_calcJacPrepare3D for the 
  ! considered element.
!</description>  

!<input>
  ! Auxiliary constants for the considered element with
  ! coordinates in Dcoord; have to be computed previously
  ! by trafo_calcJacPrepare3D.
  ! DIMENSION(TRAFO_NAUXJAC3D)
  real(DP), dimension(:), intent(in) :: DjacPrep
  
  ! Coordinates of a point on the reference element
  real(DP), intent(in) :: dparx,dpary,dparz
  
!</input>
  
!<output>
  ! The Jacobian matrix of the mapping from the reference to the 
  ! real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(3,1)
  !  Djac(4) = J(1,2)
  !  Djac(5) = J(2,2)
  !  Djac(6) = J(3,2)
  !  Djac(7) = J(1,3)
  !  Djac(8) = J(2,3)
  !  Djac(9) = J(3,3)
  real(DP), dimension(:), intent(out) :: Djac
  
  ! The determinant of the mapping.
  real(DP), intent(out) :: ddetj
  
  ! X-Coordinates of a point on the real element
  real(DP), intent(out) :: dxreal
  
  ! Y-Coordinates of a point on the real element
  real(DP), intent(out) :: dyreal

  ! Z-Coordinates of a point on the real element
  real(DP), intent(out) :: dzreal
  
!</output>
  
!</subroutine>



    ! Jacobian matrix of the mapping:
    Djac(1) = DjacPrep(2) + DjacPrep(5)*dpary + DjacPrep(6)*dparz &
            + DjacPrep(8)*dpary*dparz
    Djac(2) = DjacPrep(10) + DjacPrep(13)*dpary + DjacPrep(14)*dparz &
            + DjacPrep(16)*dpary*dparz
    Djac(3) = DjacPrep(18) + DjacPrep(21)*dpary + DjacPrep(22)*dparz &
            + DjacPrep(24)*dpary*dparz
    Djac(4) = DjacPrep(3) + DjacPrep(5)*dparx + DjacPrep(7)*dparz &
            + DjacPrep(8)*dparx*dparz
    Djac(5) = DjacPrep(11) + DjacPrep(13)*dparx + DjacPrep(15)*dparz &
            + DjacPrep(16)*dparx*dparz
    Djac(6) = DjacPrep(19) + DjacPrep(21)*dparx + DjacPrep(23)*dparz &
            + DjacPrep(24)*dparx*dparz
    Djac(7) = DjacPrep(4) + DjacPrep(6)*dparx + DjacPrep(7)*dpary &
            + DjacPrep(8)*dparx*dpary
    Djac(8) = DjacPrep(12) + DjacPrep(14)*dparx + DjacPrep(15)*dpary &
            + DjacPrep(16)*dparx*dpary
    Djac(9) = DjacPrep(20) + DjacPrep(22)*dparx + DjacPrep(23)*dpary &
            + DjacPrep(24)*dparx*dpary

    ! Determinant of the mapping
    ddetj = Djac(1)*(Djac(5)*Djac(9) - Djac(6)*Djac(8)) &
          + Djac(2)*(Djac(6)*Djac(7) - Djac(4)*Djac(9)) &
          + Djac(3)*(Djac(4)*Djac(8) - Djac(5)*Djac(7))
    
    ! Map the point to the real element
    dxreal = DjacPrep(1) + DjacPrep(2)*dparx + DjacPrep(3)*dpary &
           + DjacPrep(4)*dparz + DjacPrep(5)*dparx*dpary &
           + DjacPrep(6)*dparx*dparz + DjacPrep(7)*dpary*dparz &
           + DjacPrep(8)*dparx*dpary*dparz
    dyreal = DjacPrep(9) + DjacPrep(10)*dparx + DjacPrep(11)*dpary &
           + DjacPrep(12)*dparz + DjacPrep(13)*dparx*dpary &
           + DjacPrep(14)*dparx*dparz + DjacPrep(15)*dpary*dparz &
           + DjacPrep(16)*dparx*dpary*dparz
    dzreal = DjacPrep(17) + DjacPrep(18)*dparx + DjacPrep(19)*dpary &
           + DjacPrep(20)*dparz + DjacPrep(21)*dparx*dpary &
           + DjacPrep(22)*dparx*dparz + DjacPrep(23)*dpary*dparz &
           + DjacPrep(24)*dparx*dpary*dparz
  
  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine trafo_calcJac_hexa3D (DjacPrep,Djac,ddetj,dparx,dpary,dparz)

!<description>
  ! Calculate Jacobian determinant of mapping from reference- to
  ! real element.
  !
  ! This routine only calculates the Jacobian determinant of the
  ! mapping. This is on contrast to trafo_calcRealCoords, which not only 
  ! calculates this determinant but also maps the point. So this routine
  ! can be used to speed up the code if the coordinates of the
  ! mapped point already exist.
  !
  ! Before this routine can be called, the auxiliary factors DJF
  ! have to be calculated with trafo_calcJacPrepare for the considered element.
!</description>  

!<input>
  ! Auxiliary constants for the considered element with
  ! coordinates in Dcoord; have to be computed previously
  ! by trafo_calcJacPrepare.
  ! DIMENSION(TRAFO_NAUXJAC3D)
  real(DP), dimension(:), intent(in) :: DjacPrep
  
  ! Coordinates of a point on the reference element
  real(DP), intent(in) :: dparx,dpary,dparz
!</input>
  
!<output>
  ! The Jacobian matrix of the mapping from the reference to the 
  ! real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(3,1)
  !  Djac(4) = J(1,2)
  !  Djac(5) = J(2,2)
  !  Djac(6) = J(3,2)
  !  Djac(7) = J(1,3)
  !  Djac(8) = J(2,3)
  !  Djac(9) = J(3,3)
  real(DP), dimension(:), intent(out) :: Djac
  
  ! The determinant of the mapping.
  real(DP), intent(out) :: ddetj
!</output>
  
!</subroutine>

  
    ! Jacobian matrix of the mapping:
    Djac(1) = DjacPrep(2) + DjacPrep(5)*dpary + DjacPrep(6)*dparz &
            + DjacPrep(8)*dpary*dparz
    Djac(2) = DjacPrep(10) + DjacPrep(13)*dpary + DjacPrep(14)*dparz &
            + DjacPrep(16)*dpary*dparz
    Djac(3) = DjacPrep(18) + DjacPrep(21)*dpary + DjacPrep(22)*dparz &
            + DjacPrep(24)*dpary*dparz
    Djac(4) = DjacPrep(3) + DjacPrep(5)*dparx + DjacPrep(7)*dparz &
            + DjacPrep(8)*dparx*dparz
    Djac(5) = DjacPrep(11) + DjacPrep(13)*dparx + DjacPrep(15)*dparz &
            + DjacPrep(16)*dparx*dparz
    Djac(6) = DjacPrep(19) + DjacPrep(21)*dparx + DjacPrep(23)*dparz &
            + DjacPrep(24)*dparx*dparz
    Djac(7) = DjacPrep(4) + DjacPrep(6)*dparx + DjacPrep(7)*dpary &
            + DjacPrep(8)*dparx*dpary
    Djac(8) = DjacPrep(12) + DjacPrep(14)*dparx + DjacPrep(15)*dpary &
            + DjacPrep(16)*dparx*dpary
    Djac(9) = DjacPrep(20) + DjacPrep(22)*dparx + DjacPrep(23)*dpary &
            + DjacPrep(24)*dparx*dpary

    ! Determinant of the mapping
    ddetj = Djac(1)*(Djac(5)*Djac(9) - Djac(6)*Djac(8)) &
          + Djac(2)*(Djac(6)*Djac(7) - Djac(4)*Djac(9)) &
          + Djac(3)*(Djac(4)*Djac(8) - Djac(5)*Djac(7))
          
  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine trafo_prepJac_pyra3D(Dcoords, DjacPrep)
  
!<description>
  ! Calculate auxiliary Jacobian factors for standard pyramid.
  ! This routine builds up constant factors that are later used during
  ! the transformation from the reference element to the real element.
!</description>

!<input>
  ! Coordinates of the four corners of the real hexahedron.
  ! Dcoord(1,.) = x-coordinates, 
  ! Dcoord(2,.) = y-coordinates,
  ! Dcoord(3,.) = y-coordinates.
  ! DIMENSION(#space dimensions,NVE)
  real(DP), dimension(:,:), intent(in) :: Dcoords
!</input>
  
!<output>
  ! Arrays with auxiliary jacobian factors for later computation
  ! DIMENSION(TRAFO_NAUXJAC3D)
  real(DP), dimension(:), intent(out) :: DjacPrep
!</output>

!</subroutine>

    DjacPrep( 1) = 0.25_DP *&
                 ( Dcoords(1,1)+Dcoords(1,2)+Dcoords(1,3)+Dcoords(1,4))
    DjacPrep( 2) = 0.25_DP *&
                 (-Dcoords(1,1)+Dcoords(1,2)+Dcoords(1,3)-Dcoords(1,4))
    DjacPrep( 3) = 0.25_DP *&
                 (-Dcoords(1,1)-Dcoords(1,2)+Dcoords(1,3)+Dcoords(1,4))
    DjacPrep( 4) = 0.25_DP *&
                 (-Dcoords(1,1)-Dcoords(1,2)-Dcoords(1,3)-Dcoords(1,4))
    DjacPrep( 5) = Dcoords(1,5) - 0.25_DP * (&
                   Dcoords(1,1)+Dcoords(1,2)+Dcoords(1,3)+Dcoords(1,4))
    DjacPrep( 6) = 0.25_DP *&
                 ( Dcoords(2,1)+Dcoords(2,2)+Dcoords(2,3)+Dcoords(2,4))
    DjacPrep( 7) = 0.25_DP *&
                 (-Dcoords(2,1)+Dcoords(2,2)+Dcoords(2,3)-Dcoords(2,4))
    DjacPrep( 8) = 0.25_DP *&
                 (-Dcoords(2,1)-Dcoords(2,2)+Dcoords(2,3)+Dcoords(2,4))
    DjacPrep( 9) = 0.25_DP *&
                 (-Dcoords(2,1)-Dcoords(2,2)-Dcoords(2,3)-Dcoords(2,4))
    DjacPrep(10) = Dcoords(2,5) - 0.25_DP * (&
                   Dcoords(2,1)+Dcoords(2,2)+Dcoords(2,3)+Dcoords(2,4))
    DjacPrep(11) = 0.25_DP *&
                 ( Dcoords(3,1)+Dcoords(3,2)+Dcoords(3,3)+Dcoords(3,4))
    DjacPrep(12) = 0.25_DP *&
                 (-Dcoords(3,1)+Dcoords(3,2)+Dcoords(3,3)-Dcoords(3,4))
    DjacPrep(13) = 0.25_DP *&
                 (-Dcoords(3,1)-Dcoords(3,2)+Dcoords(3,3)+Dcoords(3,4))
    DjacPrep(14) = 0.25_DP *&
                 (-Dcoords(3,1)-Dcoords(3,2)-Dcoords(3,3)-Dcoords(3,4))
    DjacPrep(15) = Dcoords(3,5) - 0.25_DP * (&
                   Dcoords(3,1)+Dcoords(3,2)+Dcoords(3,3)+Dcoords(3,4))
  
  end subroutine 

  !************************************************************************

!<subroutine>

  pure subroutine trafo_calcTrafo_pyra3d (DjacPrep,Djac,ddetj, &
                                          dparx,dpary,dparz,&
                                          dxreal,dyreal,dzreal)

!<description>
  ! This subroutine performs two tasks:
  ! -Initialisation of a a given 3x3 matrix with the
  !  mapping information from the reference element to the "real"
  !  hexahedron. Calculation of the Jacobian determinant
  ! -Transformation of a given point on the reference element onto
  !  the "real" element
  ! Both things are performed simultaneously because the jacobian
  ! determinant is dependent of the point.
  !
  ! Before this routine can be called, the auxiliary factors DjacPrep
  ! have to be calculated with trafo_calcJacPrepare3D for the 
  ! considered element.
!</description>  

!<input>
  ! Auxiliary constants for the considered element with
  ! coordinates in Dcoord; have to be computed previously
  ! by trafo_calcJacPrepare3D.
  ! DIMENSION(TRAFO_NAUXJAC3D)
  real(DP), dimension(:), intent(in) :: DjacPrep
  
  ! Coordinates of a point on the reference element
  real(DP), intent(in) :: dparx,dpary,dparz
  
!</input>
  
!<output>
  ! The Jacobian matrix of the mapping from the reference to the 
  ! real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(3,1)
  !  Djac(4) = J(1,2)
  !  Djac(5) = J(2,2)
  !  Djac(6) = J(3,2)
  !  Djac(7) = J(1,3)
  !  Djac(8) = J(2,3)
  !  Djac(9) = J(3,3)
  real(DP), dimension(:), intent(out) :: Djac
  
  ! The determinant of the mapping.
  real(DP), intent(out) :: ddetj
  
  ! X-Coordinates of a point on the real element
  real(DP), intent(out) :: dxreal
  
  ! Y-Coordinates of a point on the real element
  real(DP), intent(out) :: dyreal

  ! Z-Coordinates of a point on the real element
  real(DP), intent(out) :: dzreal
  
!</output>
  
!</subroutine>

    ! Jacobian matrix of the mapping:
    Djac(1) = DjacPrep(2) + DjacPrep(4)*dpary
    Djac(2) = DjacPrep(7) + DjacPrep(9)*dpary
    Djac(3) = DjacPrep(12) + DjacPrep(14)*dpary
    Djac(4) = DjacPrep(3) + DjacPrep(4)*dparx
    Djac(5) = DjacPrep(8) + DjacPrep(9)*dparx
    Djac(6) = DjacPrep(13) + DjacPrep(14)*dparx
    Djac(7) = DjacPrep(5)
    Djac(8) = DjacPrep(10)
    Djac(9) = DjacPrep(15)

    ! Determinant of the mapping
    ddetj = Djac(1)*(Djac(5)*Djac(9) - Djac(6)*Djac(8)) &
          + Djac(2)*(Djac(6)*Djac(7) - Djac(4)*Djac(9)) &
          + Djac(3)*(Djac(4)*Djac(8) - Djac(5)*Djac(7))
    
    ! Map the point to the real element
    dxreal = DjacPrep(1) + DjacPrep(2)*dparx + DjacPrep(3)*dpary &
           + DjacPrep(4)*dparx*dpary + DjacPrep(5)*dparz
    dyreal = DjacPrep(6) + DjacPrep(7)*dparx + DjacPrep(8)*dpary &
           + DjacPrep(9)*dparx*dpary + DjacPrep(10)*dparz 
    dzreal = DjacPrep(11) + DjacPrep(12)*dparx + DjacPrep(13)*dpary &
           + DjacPrep(14)*dparx*dpary + DjacPrep(15)*dparz 
  
  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine trafo_calcJac_pyra3D (DjacPrep,Djac,ddetj,dparx,dpary,dparz)

!<description>
  ! Calculate Jacobian determinant of mapping from reference- to
  ! real element.
  !
  ! This routine only calculates the Jacobian determinant of the
  ! mapping. This is on contrast to trafo_calcRealCoords, which not only 
  ! calculates this determinant but also maps the point. So this routine
  ! can be used to speed up the code if the coordinates of the
  ! mapped point already exist.
  !
  ! Before this routine can be called, the auxiliary factors DJF
  ! have to be calculated with trafo_calcJacPrepare for the considered element.
!</description>  

!<input>
  ! Auxiliary constants for the considered element with
  ! coordinates in Dcoord; have to be computed previously
  ! by trafo_calcJacPrepare.
  ! DIMENSION(TRAFO_NAUXJAC3D)
  real(DP), dimension(:), intent(in) :: DjacPrep
  
  ! Coordinates of a point on the reference element
  real(DP), intent(in) :: dparx,dpary,dparz
!</input>
  
!<output>
  ! The Jacobian matrix of the mapping from the reference to the 
  ! real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(3,1)
  !  Djac(4) = J(1,2)
  !  Djac(5) = J(2,2)
  !  Djac(6) = J(3,2)
  !  Djac(7) = J(1,3)
  !  Djac(8) = J(2,3)
  !  Djac(9) = J(3,3)
  real(DP), dimension(:), intent(out) :: Djac
  
  ! The determinant of the mapping.
  real(DP), intent(out) :: ddetj
!</output>
  
!</subroutine>

    ! Jacobian matrix of the mapping:
    Djac(1) = DjacPrep(2) + DjacPrep(4)*dpary
    Djac(2) = DjacPrep(7) + DjacPrep(9)*dpary
    Djac(3) = DjacPrep(12) + DjacPrep(14)*dpary
    Djac(4) = DjacPrep(3) + DjacPrep(4)*dparx
    Djac(5) = DjacPrep(8) + DjacPrep(9)*dparx
    Djac(6) = DjacPrep(13) + DjacPrep(14)*dparx
    Djac(7) = DjacPrep(5)
    Djac(8) = DjacPrep(10)
    Djac(9) = DjacPrep(15)

    ! Determinant of the mapping
    ddetj = Djac(1)*(Djac(5)*Djac(9) - Djac(6)*Djac(8)) &
          + Djac(2)*(Djac(6)*Djac(7) - Djac(4)*Djac(9)) &
          + Djac(3)*(Djac(4)*Djac(8) - Djac(5)*Djac(7))
          
  end subroutine


  !************************************************************************

!<subroutine>

  pure subroutine trafo_prepJac_prism3D(Dcoords, DjacPrep)
  
!<description>
  ! Calculate auxiliary Jacobian factors for standard pyramid.
  ! This routine builds up constant factors that are later used during
  ! the transformation from the reference element to the real element.
!</description>

!<input>
  ! Coordinates of the four corners of the real hexahedron.
  ! Dcoord(1,.) = x-coordinates, 
  ! Dcoord(2,.) = y-coordinates,
  ! Dcoord(3,.) = y-coordinates.
  ! DIMENSION(#space dimensions,NVE)
  real(DP), dimension(:,:), intent(in) :: Dcoords
!</input>
  
!<output>
  ! Arrays with auxiliary jacobian factors for later computation
  ! DIMENSION(TRAFO_NAUXJAC3D)
  real(DP), dimension(:), intent(out) :: DjacPrep
!</output>

!</subroutine>

    DjacPrep( 1) = 0.5_DP *&
                 ( Dcoords(1,1)+Dcoords(1,4))
    DjacPrep( 2) = 0.5_DP *&
                 (-Dcoords(1,1)+Dcoords(1,2)-Dcoords(1,4)+Dcoords(1,5))
    DjacPrep( 3) = 0.5_DP *&
                 (-Dcoords(1,1)+Dcoords(1,3)-Dcoords(1,3)+Dcoords(1,6))
    DjacPrep( 4) = 0.5_DP *&
                 (-Dcoords(1,1)+Dcoords(1,4))
    DjacPrep( 5) = 0.5_DP *&
                 ( Dcoords(1,1)-Dcoords(1,2)-Dcoords(1,4)+Dcoords(1,5))
    DjacPrep( 6) = 0.5_DP *&
                 ( Dcoords(1,1)-Dcoords(1,3)-Dcoords(1,3)+Dcoords(1,6))
    DjacPrep( 7) = 0.5_DP *&
                 ( Dcoords(2,1)+Dcoords(2,4))
    DjacPrep( 8) = 0.5_DP *&
                 (-Dcoords(2,1)+Dcoords(2,2)-Dcoords(2,4)+Dcoords(2,5))
    DjacPrep( 9) = 0.5_DP *&
                 (-Dcoords(2,1)+Dcoords(2,3)-Dcoords(2,3)+Dcoords(2,6))
    DjacPrep(10) = 0.5_DP *&
                 (-Dcoords(2,1)+Dcoords(2,4))
    DjacPrep(11) = 0.5_DP *&
                 ( Dcoords(2,1)-Dcoords(2,2)-Dcoords(2,4)+Dcoords(2,5))
    DjacPrep(12) = 0.5_DP *&
                 ( Dcoords(2,1)-Dcoords(2,3)-Dcoords(2,3)+Dcoords(2,6))
    DjacPrep(13) = 0.5_DP *&
                 ( Dcoords(3,1)+Dcoords(3,4))
    DjacPrep(14) = 0.5_DP *&
                 (-Dcoords(3,1)+Dcoords(3,2)-Dcoords(3,4)+Dcoords(3,5))
    DjacPrep(15) = 0.5_DP *&
                 (-Dcoords(3,1)+Dcoords(3,3)-Dcoords(3,3)+Dcoords(3,6))
    DjacPrep(16) = 0.5_DP *&
                 (-Dcoords(3,1)+Dcoords(3,4))
    DjacPrep(17) = 0.5_DP *&
                 ( Dcoords(3,1)-Dcoords(3,2)-Dcoords(3,4)+Dcoords(3,5))
    DjacPrep(18) = 0.5_DP *&
                 ( Dcoords(3,1)-Dcoords(3,3)-Dcoords(3,3)+Dcoords(3,6))
  
  end subroutine 

  !************************************************************************

!<subroutine>

  pure subroutine trafo_calcTrafo_prism3d (DjacPrep,Djac,ddetj, &
                        dparx,dpary,dparz,dxreal,dyreal,dzreal)

!<description>
  ! This subroutine performs two tasks:
  ! -Initialisation of a a given 3x3 matrix with the
  !  mapping information from the reference element to the "real"
  !  hexahedron. Calculation of the Jacobian determinant
  ! -Transformation of a given point on the reference element onto
  !  the "real" element
  ! Both things are performed simultaneously because the jacobian
  ! determinant is dependent of the point.
  !
  ! Before this routine can be called, the auxiliary factors DjacPrep
  ! have to be calculated with trafo_calcJacPrepare3D for the 
  ! considered element.
!</description>  

!<input>
  ! Auxiliary constants for the considered element with
  ! coordinates in Dcoord; have to be computed previously
  ! by trafo_calcJacPrepare3D.
  ! DIMENSION(TRAFO_NAUXJAC3D)
  real(DP), dimension(:), intent(in) :: DjacPrep
  
  ! Coordinates of a point on the reference element
  real(DP), intent(in) :: dparx,dpary,dparz
  
!</input>
  
!<output>
  ! The Jacobian matrix of the mapping from the reference to the 
  ! real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(3,1)
  !  Djac(4) = J(1,2)
  !  Djac(5) = J(2,2)
  !  Djac(6) = J(3,2)
  !  Djac(7) = J(1,3)
  !  Djac(8) = J(2,3)
  !  Djac(9) = J(3,3)
  real(DP), dimension(:), intent(out) :: Djac
  
  ! The determinant of the mapping.
  real(DP), intent(out) :: ddetj
  
  ! X-Coordinates of a point on the real element
  real(DP), intent(out) :: dxreal
  
  ! Y-Coordinates of a point on the real element
  real(DP), intent(out) :: dyreal

  ! Z-Coordinates of a point on the real element
  real(DP), intent(out) :: dzreal
  
!</output>
  
!</subroutine>

    ! Jacobian matrix of the mapping:
    Djac(1) = DjacPrep(2) + DjacPrep(5)*dparz
    Djac(2) = DjacPrep(8) + DjacPrep(11)*dparz
    Djac(3) = DjacPrep(14) + DjacPrep(17)*dparz
    Djac(4) = DjacPrep(3) + DjacPrep(6)*dparz
    Djac(5) = DjacPrep(9) + DjacPrep(12)*dparz
    Djac(6) = DjacPrep(15) + DjacPrep(18)*dparz
    Djac(7) = DjacPrep(4) + DjacPrep(5)*dparx + DjacPrep(6)*dpary
    Djac(8) = DjacPrep(10) + DjacPrep(11)*dparx + DjacPrep(12)*dpary
    Djac(9) = DjacPrep(16) + DjacPrep(17)*dparx + DjacPrep(18)*dpary

    ! Determinant of the mapping
    ddetj = Djac(1)*(Djac(5)*Djac(9) - Djac(6)*Djac(8)) &
          + Djac(2)*(Djac(6)*Djac(7) - Djac(4)*Djac(9)) &
          + Djac(3)*(Djac(4)*Djac(8) - Djac(5)*Djac(7))
    
    ! Map the point to the real element
    dxreal = DjacPrep(1) + DjacPrep(2)*dparx + DjacPrep(3)*dpary &
           + dparz*(DjacPrep(4) + DjacPrep(5)*dparx + DjacPrep(6)*dpary)
    dyreal = DjacPrep(7) + DjacPrep(8)*dparx + DjacPrep(9)*dpary &
           + dparz*(DjacPrep(10) + DjacPrep(11)*dparx + DjacPrep(12)*dpary)
    dzreal = DjacPrep(13) + DjacPrep(14)*dparx + DjacPrep(15)*dpary &
           + dparz*(DjacPrep(16) + DjacPrep(17)*dparx + DjacPrep(18)*dpary)
  
  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine trafo_calcJac_prism3D (DjacPrep,Djac,ddetj,dparx,dpary,dparz)

!<description>
  ! Calculate Jacobian determinant of mapping from reference- to
  ! real element.
  !
  ! This routine only calculates the Jacobian determinant of the
  ! mapping. This is on contrast to trafo_calcRealCoords, which not only 
  ! calculates this determinant but also maps the point. So this routine
  ! can be used to speed up the code if the coordinates of the
  ! mapped point already exist.
  !
  ! Before this routine can be called, the auxiliary factors DJF
  ! have to be calculated with trafo_calcJacPrepare for the considered element.
!</description>  

!<input>
  ! Auxiliary constants for the considered element with
  ! coordinates in Dcoord; have to be computed previously
  ! by trafo_calcJacPrepare.
  ! DIMENSION(TRAFO_NAUXJAC3D)
  real(DP), dimension(:), intent(in) :: DjacPrep
  
  ! Coordinates of a point on the reference element
  real(DP), intent(in) :: dparx,dpary,dparz
!</input>
  
!<output>
  ! The Jacobian matrix of the mapping from the reference to the 
  ! real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(3,1)
  !  Djac(4) = J(1,2)
  !  Djac(5) = J(2,2)
  !  Djac(6) = J(3,2)
  !  Djac(7) = J(1,3)
  !  Djac(8) = J(2,3)
  !  Djac(9) = J(3,3)
  real(DP), dimension(:), intent(out) :: Djac
  
  ! The determinant of the mapping.
  real(DP), intent(out) :: ddetj
!</output>
  
!</subroutine>

    ! Jacobian matrix of the mapping:
    Djac(1) = DjacPrep(2) + DjacPrep(5)*dparz
    Djac(2) = DjacPrep(8) + DjacPrep(11)*dparz
    Djac(3) = DjacPrep(14) + DjacPrep(17)*dparz
    Djac(4) = DjacPrep(3) + DjacPrep(6)*dparz
    Djac(5) = DjacPrep(9) + DjacPrep(12)*dparz
    Djac(6) = DjacPrep(15) + DjacPrep(18)*dparz
    Djac(7) = DjacPrep(4) + DjacPrep(5)*dparx + DjacPrep(6)*dpary
    Djac(8) = DjacPrep(10) + DjacPrep(11)*dparx + DjacPrep(12)*dpary
    Djac(9) = DjacPrep(16) + DjacPrep(17)*dparx + DjacPrep(18)*dpary

    ! Determinant of the mapping
    ddetj = Djac(1)*(Djac(5)*Djac(9) - Djac(6)*Djac(8)) &
          + Djac(2)*(Djac(6)*Djac(7) - Djac(4)*Djac(9)) &
          + Djac(3)*(Djac(4)*Djac(8) - Djac(5)*Djac(7))
          
  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine trafo_calcRefCoords (ctrafoType,Dcoord,Dpoint,DparPoint)

!<description>
  ! This subroutine finds the coordinates on the reference
  ! element DparPoint for a given point Dpoint in real coordinates.
!</description>  

!<input>
  ! ID of transformation to calculate. This specifies the format of the
  ! destination array DparPoint (dimension, if to use barycentric 
  ! coordinates etc.)
  integer(I32), intent(in) :: ctrafoType
  
  ! Coordinates of the points forming the element.
  ! DIMENSION(1..ndim,1..nve)
  real(DP), dimension(:,:), intent(in) :: Dcoord

  ! Coordinates of the point in world coordinates.
  real(DP), dimension(:), intent(in) :: Dpoint
!</input>
  
!<output>
  ! Coordinates of the point in coordinates on there reference element.
  ! (The format depends on the type of coordinate on the reference element
  !  and is specified via ctrafoType.)
  real(DP), dimension(:), intent(out) :: DparPoint
!</output>
  
!</subroutine>
 integer :: iresult

    ! What type of transformation do we have? First decide on the dimension,
    ! then on the actual ID.
    select case (trafo_igetDimension(ctrafoType))
    
    case (NDIM1D)
      ! 1D elements: Lines.
      
      select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
      
      case (TRAFO_ID_MLINCUBE)
      
        ! Simple linear transformation of x \in [a,b] -> x` \in [-1,1]
        DparPoint(1) = -1.0_DP + 0.5_DP * (Dpoint(1) - Dcoord(1,1)) / (Dcoord(2,1)-Dcoord(1,1))
      
      end select
    
    case (NDIM2D)
      ! 2D elements. Triangles, Quadrilaterals. Check the actual transformation
      ! ID how to transform.
    
      select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
      
      case (TRAFO_ID_LINSIMPLEX)
      
        ! 2D simplex -> linear triangular transformation.
        !
        ! Call the corresponding routine in geometryaux to calculate
        ! the barycentric coordinates to Dpoint.
        call gaux_getBarycentricCoords_tri2D(&
          Dcoord,Dpoint(1),Dpoint(2),DparPoint(1),DparPoint(2),DparPoint(3))      
          
      case (TRAFO_ID_MLINCUBE)
      
        ! Bilinear transformation for cubic-shaped elements 
        ! -> Bilinear quadrilateral transformation.
        !
        ! Call the back transformation routine from above.
        call trafo_calcBackTrafo_quad2d (Dcoord, &
            Dpoint(1),Dpoint(2),DparPoint(1),DparPoint(2))
        
      end select ! actual ID
    case (NDIM3D)
      ! trilinear transformation for hexahedral elements 
      call trafo_calcBackTrafo_hexa(Dcoord, &
           Dpoint(1),Dpoint(2),Dpoint(3), &
           DparPoint(1),DparPoint(2),DparPoint(3),iresult)
      
    end select ! Dimension
  
  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine trafo_calcRealCoords_general (ctrafoType,Dcoord,&
    DpointRef,Dpoint)

!<description>
  ! This subroutine maps the coordinates of a point DpointRef on the
  ! reference element to a point Dpoint in world coordinates.
!</description>  

!<input>
  ! ID of transformation to calculate. This specifies the format of the
  ! destination array DparPoint (dimension, if to use barycentric 
  ! coordinates etc.)
  integer(I32), intent(in) :: ctrafoType
  
  ! Coordinates of the points forming the element.
  ! DIMENSION(1..ndim,1..nve)
  real(DP), dimension(:,:), intent(in) :: Dcoord

  ! Coordinates of the point in coordinates on there reference element.
  ! (The format depends on the type of coordinate on the reference element
  !  and is specified via ctrafoType.)
  real(DP), dimension(:), intent(in) :: DpointRef
!</input>
  
!<output>
  ! Coordinates of the point in world coordinates.
  real(DP), dimension(:), intent(out) :: Dpoint
!</output>
  
!</subroutine>

    ! local variables
    real(DP), dimension(TRAFO_NAUXJACMAX) :: Djac
    real(DP) :: ddetj

    ! What type of transformation do we have? First decide on the dimension,
    ! then on the actual ID.
    select case (trafo_igetDimension(ctrafoType))
    
    case (NDIM1D)
      ! 1D elements: Lines.
      
      select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
      
      case (TRAFO_ID_MLINCUBE)
      
        ! Simple linear transformation of x \in [-1,1] -> x` \in [a,b]
        Dpoint(1) = 0.5_DP*(Dcoord(2,1)-Dcoord(1,1))*(DpointRef(1)+1.0_DP)
      
      end select
    
    case (NDIM2D)
      ! 2D elements. Triangles, Quadrilaterals. Check the actual transformation
      ! ID how to transform.
    
      select case (iand(ctrafoType,TRAFO_DIM_IDMASK))
      
      case (TRAFO_ID_LINSIMPLEX)
      
        ! 2D simplex -> linear triangular transformation.
        !
        ! The world coordinates can be obtained by simple
        ! linear combination.
        Dpoint(1) = Dcoord(1,1)*DpointRef(1) + &
                    Dcoord(1,2)*DpointRef(2) + &
                    Dcoord(1,3)*DpointRef(3)
        Dpoint(2) = Dcoord(2,1)*DpointRef(1) + &
                    Dcoord(2,2)*DpointRef(2) + &
                    Dcoord(2,3)*DpointRef(3)
          
      case (TRAFO_ID_MLINCUBE)
      
        ! Bilinear transformation for cubic-shaped elements 
        ! -> Bilinear quadrilateral transformation.
        !
        ! This is a little bit more complicated, we need the whole 
        ! transformation. The Djac / Ddetj information is ignored.
        call trafo_calctrafo (ctrafoType,Dcoord,&
                              DpointRef,Djac,Ddetj,Dpoint)
        
      end select ! actual ID
    
    end select ! Dimension
  
  end subroutine

  !****************************************************************************
  
  ! The following routines (trafo_mapCubPtsRef2LvlXXXX) perform a
  ! '2-level-mapping' of cubature points on reference cells which is needed
  ! by the 2-Level-Mass matrix assembly in multileveloperators.f90.
  ! The 2-Level-Mass assembly needs to integrate over two FE basis functions,
  ! which, however, are defined on two different meshes: a coarse mesh
  ! and its 2-level-ordering-based refined mesh.
  ! We now need to specify a set of cubature points on both the coarse and
  ! fine mesh, such that the transformed 'real' points are equal. This means
  ! that if we have a set of cubature points given for the fine mesh, we need
  ! to transform these reference coordinates to reference coordinates on the
  ! coarse mesh.
  !
  ! Consider the following 1D example:
  !
  ! Let us assume that we have a set of 3 cubature points (a,b,c) on the
  ! 1D reference interval [-1,1]:
  !
  !   |----a----b----c----|
  !  -1                   1
  !
  ! Now during integration the trafo will be called to transform these 3
  ! points to 'real' coordinates on a fine mesh edge (X_i,X_i+1):
  !
  !    |------a------b------c------|
  !   X_i                        X_i+1
  !
  ! Now since this edge is on the fine mesh, it must have been refined from
  ! a coarse mesh edge (Y_i,Y_i+1):
  !
  !   Y_i                                                    Y_i+1
  !    |-------------------------------------------------------|
  !
  !                  2-Level-Ordering Refinement =>
  !
  !    |---------------------------|---------------------------|
  !   X_i                        X_i+1                       X_i+2
  !
  !                               OR
  !
  !    |---------------------------|---------------------------|
  !  X_i-1                        X_i                        X_i+1
  !
  ! (Note: It does not matter whether the edge (Y_i,Y_i+1) was refined into
  !  (X_i,X_i+1,X_i+2) or (X_i-1,X_i,X_i+1), so we will assume the first case)
  !
  ! Now we need to map the cubature points which lie on the fine mesh edges
  ! (X_i,X_i+1) and (X_i+1,X_i+2) onto the coarse mesh edge (Y_i,Y_i+1):
  !
  !    |------r------s------t-------------u------v------w------|
  !   Y_i     ^      ^      ^             ^      ^      ^    Y_i+1
  !           |      |      |             |      |      |
  !    |------a------b------c------|------a------b------c------|
  !   X_i                        X_i+1                       X_i+2
  !
  ! Afterwards we need to transform the 6 cubature points from the coarse
  ! mesh edge (Y_i,Y_i+1) to the reference interval [-1,1] as most of
  ! the FE basis functions need reference coordinates for evaluation:
  !
  !   |--r--s--t-----u--v--w--|
  !  -1                       1
  !
  ! But since we know (or better: we silently assume) that the mesh was
  ! refined using the 2-level-ordering algorithm, we can skip the whole
  ! transformations and directly map the points on the reference interval:
  !
  !   |----a----b----c----|  =>   |--r--s--t-----u--v--w--|
  !  -1                   1      -1                       1
  !
  ! And this is exactly what the trafo_mapCubPtsRef2LvlXXXX routines do.
  !
  ! The same thing holds for quadrilaterals in 2D and hexahedra in 3D.
  
!<subroutine>

  pure subroutine trafo_mapCubPtsRef2LvlEdge1D(ncubp, Dfine, Dcoarse)

!<description>
  ! This routine maps the coordinates of the cubature points of a fine mesh
  ! edge onto the coarse mesh edge.
!</description>

!<input>
  ! Number of cubature points for the edge
  integer , intent(in) :: ncubp
  
  ! Cubature point coordinates on 1D reference fine mesh edge [-1,1]
  ! Dfine(1,1:ncubp)= X-coordinates
  real(DP), dimension(:,:), intent(in) :: Dfine
!</input>
  
!<output>
  ! Cubature point coordinates on 1D reference coarse mesh edge [-1,1]
  ! Dcoarse(1,1:*ncubp)= X-coordinates
  real(DP), dimension(:,:), intent(out) :: Dcoarse
!</output>

!</subroutine>

  ! local variables
  integer :: i
  real(DP) :: dx
  
    do i = 1, ncubp
      dx = 0.5_DP * (Dfine(1,i) - 1.0_DP)
      Dcoarse(1,      i) = dx
      Dcoarse(1,ncubp+i) = dx + 1.0_DP
    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  pure subroutine trafo_mapCubPtsRef2LvlTri2D(ncubp, Dfine, Dcoarse)

!<description>
  ! This routine maps the coordinates of the cubature points of a fine mesh
  ! triangle onto the coarse mesh triangle.
!</description>

!<input>
  ! Number of cubature points for the triangle
  integer , intent(in) :: ncubp
  
  ! Cubature point coordinates for the fine mesh triangle given in 
  ! barycentric coordinates
  ! Dfine(1:3,1:ncubp) = barycentric coordinates
  real(DP), dimension(:,:), intent(in) :: Dfine
!</input>
  
!<output>
  ! Cubature point coordinates for the coarse mesh triangle given in
  ! barycentric coordinates.
  ! Dcoarse(1:3,1..4*ncubp) = barycentric coordinates
  real(DP), dimension(:,:), intent(out) :: Dcoarse
!</output>

!</subroutine>

  ! local variables
  integer :: i
  real(DP) :: d1,d2,d3
  
  ! Vertice / Edge-Midpoint coordinates
  real(DP), parameter :: V1x = 0.0_DP
  real(DP), parameter :: V1y = 0.0_DP
  real(DP), parameter :: V2x = 1.0_DP
  real(DP), parameter :: V2y = 0.0_DP
  real(DP), parameter :: V3x = 0.0_DP
  real(DP), parameter :: V3y = 1.0_DP
  real(DP), parameter :: E1x = 0.5_DP
  real(DP), parameter :: E1y = 0.0_DP
  real(DP), parameter :: E2x = 0.5_DP
  real(DP), parameter :: E2y = 0.5_DP
  real(DP), parameter :: E3x = 0.0_DP
  real(DP), parameter :: E3y = 0.5_DP
  
    ! Now as we do not use reference coordinates and a reference triangle
    ! but barycentric coordinates, the 2-Level mapping of cubature points
    ! on triangles gets a bit more tricky.
    ! The very first thing that one has to be aware of is that one of the
    ! barycentric coordinates is redundant as it can be calculated as 1
    ! minus the sum of the other two coordinates. So we can focus on
    ! calculating two barycentric coordinates to get the job done.
    !
    ! The trick we will use here can often be used when working with
    ! barycentric coordinates, as they are quite unhandy: We will simply
    ! define a reference triangle and a corresponding transformation from
    ! barycentric coordinates into reference coordinates and vice versa.
    ! Our reference triangle will be defined as the convex area with the
    ! three corner vertices V1 := (0,0), V2 := (1,0) and V3 := (0,1).
    ! Now if we have a point (p1,p2,p3) given in barycentric coordinates,
    ! we can transform the point into reference coordinates by using:
    !
    !                      (x,y) = p1*V1 + p2*V2 + p3*V3
    !
    ! Now inserting the definition of our corner vertices Vi we get:
    !
    !                      (x,y) = (p2, p3)
    !
    ! On the other hand, if we want to transform a point (x,y) given in
    ! reference coordinates into barycentric coordinates, we do this by:
    !
    !                 (p1,p2,p3) = (1-x-y, x, y)
    !
    ! Now in the following, we will also use the edge-midpoints given in
    ! reference coordinates to calculate the reference coordinates of the
    ! 2-Level mapping of the cubature points.
    !
    ! Take a look at the following ASCII-Art of a refined triangle:
    !
    !  y ^ 
    !    |
    !    1    V3
    !    |    |\__
    !    |    | 1 \__
    !    |    |      \__
    !    |    |         \__
    !    |    |     D      \__
    !    |    | 2           3 \
    !    |    E3---------------E2
    !    |    |\__3          2 |\__
    !    |    | 3 \__    A     | 2 \__
    !    |    |      \__       |      \__
    !    |    |         \__    |         \__
    !    |    |     B      \__1|     C      \__
    !    |    | 1           2 \| 3           1 \
    !    0    V1---------------E1---------------V2
    !
    !         0---------------------------------1--->
    !                                               x
    !
    ! Assume that we have a point (p1,p2,p3) given in barycentric coordinates
    ! on the triangle C. To map the point (p1,p2,p3) onto the coarse mesh
    ! triangle, all we have to calculate is:
    !
    !                 (x,y) = p1*V2 + p2*E2 + p3*E1
    !
    ! To figure out the rest, the reader is assumed to switch on his/her
    ! brain (or whatever organ the reader uses for thinking) and take
    ! a closer look at the following algorithm ^_^

    ! Remark:
    ! Please note that the following implementation focuses on being
    ! understandable rather than being efficient!
    do i = 1, ncubp
    
      ! Get fine mesh barycentric coordinates
      d1 = Dfine(1,i)
      d2 = Dfine(2,i)
      d3 = Dfine(3,i)
      
      ! Child A -> reference coordinates
      Dcoarse(2,        i) = d1*E1x + d2*E2x + d3*E3x
      Dcoarse(3,        i) = d1*E1y + d2*E2y + d3*E3y
      
      ! Child B -> reference coordinates
      Dcoarse(2,  ncubp+i) = d1*V1x + d2*E1x + d3*E3x
      Dcoarse(3,  ncubp+i) = d1*V1y + d2*E1y + d3*E3y
      
      ! Child C -> reference coordinates
      Dcoarse(2,2*ncubp+i) = d1*V2x + d2*E2x + d3*E1x
      Dcoarse(3,2*ncubp+i) = d1*V2y + d2*E2y + d3*E1y
      
      ! Child D -> reference coordinates
      Dcoarse(2,3*ncubp+i) = d1*V3x + d2*E3x + d3*E2x
      Dcoarse(3,3*ncubp+i) = d1*V3y + d2*E3y + d3*E2y
      
    end do
    
    ! Recalculate first barycentric coordinate
    do i = 1, 4*ncubp
      Dcoarse(1,i) = 1.0_DP - (Dcoarse(2,i) + Dcoarse(3,i))
    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  pure subroutine trafo_mapCubPtsRef2LvlQuad2D(ncubp, Dfine, Dcoarse)

!<description>
  ! This routine maps the coordinates of the cubature points of a fine mesh
  ! quadrilateral onto the coarse mesh quadrilateral.
!</description>

!<input>
  ! Number of cubature points for the quadrilateral
  integer , intent(in) :: ncubp
  
  ! Cubature point coordinates on 2D reference fine mesh quad [-1,1]x[-1,1]
  ! Dfine(1,1:ncubp) = X-coordinates,
  ! Dfine(2,1:ncubp) = Y-coordinates
  real(DP), dimension(:,:), intent(in) :: Dfine
!</input>
  
!<output>
  ! Cubature point coordinates on 2D reference coarse mesh quad [-1,1]x[-1,1]
  ! Dcoarse(1,1:4*ncubp) = X-coordinates,
  ! Dcoarse(2,1:4*ncubp) = Y-coordinates
  real(DP), dimension(:,:), intent(out) :: Dcoarse
!</output>

!</subroutine>

  ! local variables
  integer :: i
  real(DP) :: dx,dy
  
    do i = 1, ncubp
      dx = 0.5_DP * (Dfine(1,i) - 1.0_DP)
      dy = 0.5_DP * (Dfine(2,i) - 1.0_DP)
      Dcoarse(1:2,        i) = (/ dx,  dy/)
      Dcoarse(1:2,  ncubp+i) = (/-dy,  dx/)
      Dcoarse(1:2,2*ncubp+i) = (/-dx, -dy/)
      Dcoarse(1:2,3*ncubp+i) = (/ dy, -dx/)
    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  pure subroutine trafo_mapCubPtsRef2LvlHexa3D(ncubp, Dfine, Dcoarse)

!<description>
  ! This routine maps the coordinates of the cubature points of a fine mesh
  ! hexahedron onto the coarse mesh hexahedron.
!</description>

!<input>
  ! Number of cubature points for the hexahedron
  integer , intent(in) :: ncubp
  
  ! Cubature point coordinates on 3D reference fine mesh hexahedron [-1,1]^3
  ! Dfine(1,1:ncubp) = X-coordinates,
  ! Dfine(2,1:ncubp) = Y-coordinates,
  ! Dfine(3,1:ncubp) = Z-coordinates
  real(DP), dimension(:,:), intent(in) :: Dfine
!</input>
  
!<output>
  ! Cubature point coordinates on 3D reference coarse mesh hexahedron [-1,1]^3
  ! Dcoarse(1,1:8*ncubp) = X-coordinates,
  ! Dcoarse(2,1:8*ncubp) = Y-coordinates,
  ! Dcoarse(3,1:8*ncubp) = Z-coordinates
  real(DP), dimension(:,:), intent(out) :: Dcoarse
!</output>

!</subroutine>

  ! local variables
  integer :: i
  real(DP) :: dx,dy,dz
  
    do i = 1, ncubp
      dx = 0.5_DP * (Dfine(1,i) - 1.0_DP)
      dy = 0.5_DP * (Dfine(2,i) - 1.0_DP)
      dz = 0.5_DP * (Dfine(3,i) - 1.0_DP)
      Dcoarse(1:3,        i) = (/ dx,  dy,  dz/)
      Dcoarse(1:3,  ncubp+i) = (/-dy,  dx,  dz/)
      Dcoarse(1:3,2*ncubp+i) = (/-dx, -dy,  dz/)
      Dcoarse(1:3,3*ncubp+i) = (/ dy, -dx,  dz/)
      Dcoarse(1:3,4*ncubp+i) = (/ dy,  dx, -dz/)
      Dcoarse(1:3,5*ncubp+i) = (/-dx,  dy, -dz/)
      Dcoarse(1:3,6*ncubp+i) = (/-dy, -dx, -dz/)
      Dcoarse(1:3,7*ncubp+i) = (/ dx, -dy, -dz/)
    end do

  end subroutine
  
!************************************************************************  
  
!<subroutine>  
  pure subroutine trafo_calcBackTrafo_hexa(Dcoord, &
     dxreal,dyreal,dzreal,dparx,dpary,dparz,iresult)
!<description>
  ! This subroutine is to find the coordinates on the reference
  ! element (dparx,dpary,dparz) for a given point (dxreal,dyreal,dzreal) 
  ! in real coordinates.
  !
!<input>
  ! Coordinates of the eight points forming the element.
  ! Dcoord(1,.) = x-coordinates,
  ! Dcoord(2,.) = y-coordinates.
  ! Dcoord(3,.) = z-coordinates.
  ! DIMENSION(#space dimensions,NVE)
  real(DP), dimension(:,:), intent(in) :: Dcoord 
  
  ! Coordinates of a point on the real element
  real(DP), intent(in) :: dxreal,dyreal,dzreal
!</input>
  
!<output>
  ! Coordinates of a point on the reference element
  real(DP), intent(out) :: dparx,dpary,dparz
  integer,intent(out)   :: iresult
!</output>
  
!</subroutine>
  ! DIMENSION(TRAFO_NAUXJAC3D)
  real(DP), dimension(TRAFO_NAUXJAC3D) :: DjacPrep
  real(DP), dimension(9) :: Djac
  real(DP), dimension(3) :: DPoints
  real(DP)               :: ddetj,ddetjx,ddetjy,ddetjz
  real(dp) :: rhsx,rhsy,rhsz,eps
  real(DP) :: dparx_old,dpary_old,dparz_old,c1,c2,c3
  integer  :: i
  ! initialize
  dparx=0.0_dp
  dpary=0.0_dp
  dparz=0.0_dp
  eps  =1e-6
  
  !DEBUG
  !DPoints(:)=(/-1,-1,-1/)
  
  ! a1
  DjacPrep( 1) = 0.125_DP *&
               ( DCoord(1,1)+DCoord(1,2)+DCoord(1,3)+DCoord(1,4)&
                +DCoord(1,5)+DCoord(1,6)+DCoord(1,7)+DCoord(1,8))
  ! a2                
  DjacPrep( 2) = 0.125_DP *&
               (-DCoord(1,1)+DCoord(1,2)+DCoord(1,3)-DCoord(1,4)&
                -DCoord(1,5)+DCoord(1,6)+DCoord(1,7)-DCoord(1,8))
  ! a3                
  DjacPrep( 3) = 0.125_DP *&
               (-DCoord(1,1)-DCoord(1,2)+DCoord(1,3)+DCoord(1,4)&
                -DCoord(1,5)-DCoord(1,6)+DCoord(1,7)+DCoord(1,8))
  ! a4
  DjacPrep( 4) = 0.125_DP *&
               (-DCoord(1,1)-DCoord(1,2)-DCoord(1,3)-DCoord(1,4)&
                +DCoord(1,5)+DCoord(1,6)+DCoord(1,7)+DCoord(1,8))
  ! a5
  DjacPrep( 5) = 0.125_DP *&
               ( DCoord(1,1)-DCoord(1,2)+DCoord(1,3)-DCoord(1,4)&
                +DCoord(1,5)-DCoord(1,6)+DCoord(1,7)-DCoord(1,8))
  ! a6                
  DjacPrep( 6) = 0.125_DP *&
               ( DCoord(1,1)-DCoord(1,2)-DCoord(1,3)+DCoord(1,4)&
                -DCoord(1,5)+DCoord(1,6)+DCoord(1,7)-DCoord(1,8))
  ! a7
  DjacPrep( 7) = 0.125_DP *&
               ( DCoord(1,1)+DCoord(1,2)-DCoord(1,3)-DCoord(1,4)&
                -DCoord(1,5)-DCoord(1,6)+DCoord(1,7)+DCoord(1,8))
  ! a8
  DjacPrep( 8) = 0.125_DP *&
               (-DCoord(1,1)+DCoord(1,2)-DCoord(1,3)+DCoord(1,4)&
                +DCoord(1,5)-DCoord(1,6)+DCoord(1,7)-DCoord(1,8))
  DjacPrep( 9) = 0.125_DP *&
               ( DCoord(2,1)+DCoord(2,2)+DCoord(2,3)+DCoord(2,4)&
                +DCoord(2,5)+DCoord(2,6)+DCoord(2,7)+DCoord(2,8))
  DjacPrep(10) = 0.125_DP *&
               (-DCoord(2,1)+DCoord(2,2)+DCoord(2,3)-DCoord(2,4)&
                -DCoord(2,5)+DCoord(2,6)+DCoord(2,7)-DCoord(2,8))
  DjacPrep(11) = 0.125_DP *&
               (-DCoord(2,1)-DCoord(2,2)+DCoord(2,3)+DCoord(2,4)&
                -DCoord(2,5)-DCoord(2,6)+DCoord(2,7)+DCoord(2,8))
  DjacPrep(12) = 0.125_DP *&
               (-DCoord(2,1)-DCoord(2,2)-DCoord(2,3)-DCoord(2,4)&
                +DCoord(2,5)+DCoord(2,6)+DCoord(2,7)+DCoord(2,8))
  DjacPrep(13) = 0.125_DP *&
               ( DCoord(2,1)-DCoord(2,2)+DCoord(2,3)-DCoord(2,4)&
                +DCoord(2,5)-DCoord(2,6)+DCoord(2,7)-DCoord(2,8))
  DjacPrep(14) = 0.125_DP *&
               ( DCoord(2,1)-DCoord(2,2)-DCoord(2,3)+DCoord(2,4)&
                -DCoord(2,5)+DCoord(2,6)+DCoord(2,7)-DCoord(2,8))
  DjacPrep(15) = 0.125_DP *&
               ( DCoord(2,1)+DCoord(2,2)-DCoord(2,3)-DCoord(2,4)&
                -DCoord(2,5)-DCoord(2,6)+DCoord(2,7)+DCoord(2,8))
  DjacPrep(16) = 0.125_DP *&
               (-DCoord(2,1)+DCoord(2,2)-DCoord(2,3)+DCoord(2,4)&
                +DCoord(2,5)-DCoord(2,6)+DCoord(2,7)-DCoord(2,8))
  DjacPrep(17) = 0.125_DP *&
               ( DCoord(3,1)+DCoord(3,2)+DCoord(3,3)+DCoord(3,4)&
                +DCoord(3,5)+DCoord(3,6)+DCoord(3,7)+DCoord(3,8))
  DjacPrep(18) = 0.125_DP *&
               (-DCoord(3,1)+DCoord(3,2)+DCoord(3,3)-DCoord(3,4)&
                -DCoord(3,5)+DCoord(3,6)+DCoord(3,7)-DCoord(3,8))
  DjacPrep(19) = 0.125_DP *&
               (-DCoord(3,1)-DCoord(3,2)+DCoord(3,3)+DCoord(3,4)&
                -DCoord(3,5)-DCoord(3,6)+DCoord(3,7)+DCoord(3,8))
  DjacPrep(20) = 0.125_DP *&
               (-DCoord(3,1)-DCoord(3,2)-DCoord(3,3)-DCoord(3,4)&
                +DCoord(3,5)+DCoord(3,6)+DCoord(3,7)+DCoord(3,8))
  DjacPrep(21) = 0.125_DP *&
               ( DCoord(3,1)-DCoord(3,2)+DCoord(3,3)-DCoord(3,4)&
                +DCoord(3,5)-DCoord(3,6)+DCoord(3,7)-DCoord(3,8))
  DjacPrep(22) = 0.125_DP *&
               ( DCoord(3,1)-DCoord(3,2)-DCoord(3,3)+DCoord(3,4)&
                -DCoord(3,5)+DCoord(3,6)+DCoord(3,7)-DCoord(3,8))
  DjacPrep(23) = 0.125_DP *&
               ( DCoord(3,1)+DCoord(3,2)-DCoord(3,3)-DCoord(3,4)&
                -DCoord(3,5)-DCoord(3,6)+DCoord(3,7)+DCoord(3,8))
  DjacPrep(24) = 0.125_DP *&
               (-DCoord(3,1)+DCoord(3,2)-DCoord(3,3)+DCoord(3,4)&
                +DCoord(3,5)-DCoord(3,6)+DCoord(3,7)-DCoord(3,8))
  
  !  dparx = a1 + a2*xi1 + a3*xi2 + a4*xi3 + a5*xi1*xi2 
  !      + a6*xi1*xi3 + a7*xi2*xi3 + a8*xi1*xi2*xi3  dxreal,dyreal,dzreal
  ! DEBUG
!  DPoints(1)=-1.0_dp
!  DPoints(2)=-1.0_dp
!  DPoints(3)=-1.0_dp
!  
!  dparx_old = DjacPrep(1)+DjacPrep(2)*DPoints(1)+DjacPrep(3)*DPoints(2)+DjacPrep(4)*DPoints(3) &
!              +DjacPrep(5)*DPoints(1)*DPoints(2)+ DjacPrep(6)*DPoints(1)*DPoints(3) &
!              +DjacPrep(7)*DPoints(2)*DPoints(3)+DjacPrep(8)*DPoints(1)*DPoints(2)*DPoints(3)  
!
!  dpary_old = DjacPrep(9)+DjacPrep(10)*DPoints(1)+DjacPrep(11)*DPoints(2)+DjacPrep(12)*DPoints(3) &
!              +DjacPrep(13)*DPoints(1)*DPoints(2)+ DjacPrep(14)*DPoints(1)*DPoints(3) &
!              +DjacPrep(15)*DPoints(2)*DPoints(3)+DjacPrep(16)*DPoints(1)*DPoints(2)*DPoints(3)  
!
!  dparz_old = DjacPrep(17)+DjacPrep(18)*DPoints(1)+DjacPrep(19)*DPoints(2)+DjacPrep(20)*DPoints(3) &
!              +DjacPrep(21)*DPoints(1)*DPoints(2)+ DjacPrep(22)*DPoints(1)*DPoints(3) &
!              +DjacPrep(23)*DPoints(2)*DPoints(3)+DjacPrep(24)*DPoints(1)*DPoints(2)*DPoints(3)  
  
  do i=1,20
  
    dparx_old=dparx
    dpary_old=dpary
    dparz_old=dparz
    
    ! Jacobian matrix of the mapping:
    Djac(1) = DjacPrep(2) + DjacPrep(5)*dpary + DjacPrep(6)*dparz &
            + DjacPrep(8)*dpary*dparz
    Djac(2) = DjacPrep(10) + DjacPrep(13)*dpary + DjacPrep(14)*dparz &
            + DjacPrep(16)*dpary*dparz
    Djac(3) = DjacPrep(18) + DjacPrep(21)*dpary + DjacPrep(22)*dparz &
            + DjacPrep(24)*dpary*dparz
    Djac(4) = DjacPrep(3) + DjacPrep(5)*dparx + DjacPrep(7)*dparz &
            + DjacPrep(8)*dparx*dparz
    Djac(5) = DjacPrep(11) + DjacPrep(13)*dparx + DjacPrep(15)*dparz &
            + DjacPrep(16)*dparx*dparz
    Djac(6) = DjacPrep(19) + DjacPrep(21)*dparx + DjacPrep(23)*dparz &
            + DjacPrep(24)*dparx*dparz
    Djac(7) = DjacPrep(4) + DjacPrep(6)*dparx + DjacPrep(7)*dpary &
            + DjacPrep(8)*dparx*dpary
    Djac(8) = DjacPrep(12) + DjacPrep(14)*dparx + DjacPrep(15)*dpary &
            + DjacPrep(16)*dparx*dpary
    Djac(9) = DjacPrep(20) + DjacPrep(22)*dparx + DjacPrep(23)*dpary &
            + DjacPrep(24)*dparx*dpary

    ! build the RHS vector
    rhsx = DjacPrep(1) + DjacPrep(2)*dparx + DjacPrep(3)*dpary &
           + DjacPrep(4)*dparz + DjacPrep(5)*dparx*dpary &
           + DjacPrep(6)*dparx*dparz + DjacPrep(7)*dpary*dparz &
           + DjacPrep(8)*dparx*dpary*dparz - dxreal
    rhsy = DjacPrep(9) + DjacPrep(10)*dparx + DjacPrep(11)*dpary &
           + DjacPrep(12)*dparz + DjacPrep(13)*dparx*dpary &
           + DjacPrep(14)*dparx*dparz + DjacPrep(15)*dpary*dparz &
           + DjacPrep(16)*dparx*dpary*dparz - dyreal
    rhsz = DjacPrep(17) + DjacPrep(18)*dparx + DjacPrep(19)*dpary &
           + DjacPrep(20)*dparz + DjacPrep(21)*dparx*dpary &
           + DjacPrep(22)*dparx*dparz + DjacPrep(23)*dpary*dparz &
           + DjacPrep(24)*dparx*dpary*dparz - dzreal
    
    ! We have the rhs, now solve the system using cramer's rule
           
!    ! determinant 
!    ddetj = Djac(1)*(Djac(5)*Djac(9) - Djac(6)*Djac(8)) &
!          - Djac(2)*(Djac(4)*Djac(9) - Djac(6)*Djac(7)) &
!          + Djac(3)*(Djac(4)*Djac(8) - Djac(5)*Djac(7))
          
    ! Determinant of the mapping
    ddetj = Djac(1)*(Djac(5)*Djac(9) - Djac(6)*Djac(8)) &
          + Djac(2)*(Djac(6)*Djac(7) - Djac(4)*Djac(9)) &
          + Djac(3)*(Djac(4)*Djac(8) - Djac(5)*Djac(7))
          
          
    ! determinant substituting the 1st-column with the rhs
    ddetjx = rhsx*(Djac(5)*Djac(9) - Djac(6)*Djac(8)) &
           - Djac(4)*(rhsy*Djac(9) - Djac(8)*rhsz) &
           + Djac(7)*(rhsy*Djac(6) - Djac(5)*rhsz)
           
    ! determinant substituting the 2nd-column with the rhs
    ddetjy = Djac(1)*(rhsy*Djac(9) - Djac(8)*rhsz) &
           - rhsx*(Djac(2)*Djac(9) - Djac(3)*Djac(8)) &
           + Djac(3)*(Djac(2)*rhsz - rhsy*Djac(3))
           
    ! determinant substituting the 3rd-column with the rhs
    ddetjz = Djac(1)*(Djac(5)*rhsz - rhsy*Djac(6)) &
           - Djac(4)*(Djac(2)*rhsz - rhsy*Djac(3)) &
           + rhsx*(Djac(2)*Djac(6) - Djac(3)*Djac(5))

    ! calculate the solution
    dparx = ddetjx/ddetj
    dpary = ddetjy/ddetj
    dparz = ddetjz/ddetj
    
    ! DEBUG
!    c1=Djac(1)*dparx + Djac(2)*dpary + Djac(3)*dparz
!    
!    c2=Djac(4)*dparx + Djac(5)*dpary + Djac(6)*dparz
!    
!    c3=Djac(7)*dparx + Djac(8)*dpary + Djac(9)*dparz
    
    ! compute the correction
    dparx = dparx_old-dparx
    dpary = dpary_old-dpary
    dparz = dparz_old-dparz
    
    if((dparx .lt. -2.0_dp) .or. (dparx .gt. 2.0_dp))then
      dparx_old=dparx
      dpary_old=dpary
      dparz_old=dparz
      return
    end if

    if((dpary .lt. -2.0_dp) .or. (dpary .gt. 2.0_dp))then
      dparx_old=dparx
      dpary_old=dpary
      dparz_old=dparz
      return
    end if

    if((dparz .lt. -2.0_dp) .or. (dparz .gt. 2.0_dp))then
      dparx_old=dparx
      dpary_old=dpary
      dparz_old=dparz
      return
    end if


    if( (abs(dparx_old-dparx) .le. eps) .and. (abs(dpary_old-dpary) .le. eps) .and. &
        (abs(dparz_old-dparz) .le. eps) )then
        return
    end if
  
  end do
  
  if(i .gt. 20)then
    iresult=-1
  end if

  end subroutine trafo_calcBackTrafo_hexa
  
!************************************************************************
  

end module 
