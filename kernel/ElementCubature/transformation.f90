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
!# 6.) trafo_calctrafo_mult
!#     -> Calculates the transformation for multiple points on one element.
!#        Supports triangular and quadrilateral mapping.
!#
!# 7.) trafo_calctrafo_sim
!#     -> Calculate the transformation for multiple points on multiple
!#        elements. Supports triangular and quadrilateral mapping.
!#
!# 8.) trafo_calcJacPrepare
!#     -> calculates auxiliary Jacobian factors for the transformation
!#        from a reference quadrilateral to a real quadrilateral
!#
!# 9.) trafo_calcJac
!#     -> calculates the Jacobian matrix + Jacobian determinant of the mapping
!#        from  the reference to a real quadrilateral element
!#
!# 10.) trafo_calcRealCoords
!#      -> maps a point from the reference element to the real element
!#
!# 11.) trafo_calcTrafo
!#      -> calculates the Jacobian matrix + Jacobian determinant of the mapping
!#         from  the reference to a real quadrilateral element
!#       -> maps a point from the reference element to the real element
!#      (so performing the same task as elem_calcJac and elem_calcRealCoords
!#       in one routine)
!#
!# 12.) trafo_mapCubPts1Dto2DRefQuad
!#      -> Maps a set of 1D cubature point coordinates on the reference 
!#         interval to an edge of the 2D reference quadrilateral.
!#
!#  FAQ - Some explainations
!# --------------------------
!# 1.) How is the ID code ctrafoType of the transformation defined?
!#
!#  The ID code of a transformation is a 32 bit integer field. The different
!#  bits in the transformation encode the type of transformation, the dimension
!#  and additional information, as far as necessary. The exact layout is
!#  defined as follows:
!#
!#%     Bit | 31 ... 24 23 ... 16 | 15 ... 10 | 9   8 | 7 ............ 0|
!#%     -----------------------------------------------------------------
!#%         |         ****        | unused    |dimens |     trafo-ID    |
!#
!#   Bits 0..7   specify the type of transformation (triangular, quad, 
!#               linear, quadratic, ...). The different ID's for this
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
!#
!#   Note: As long as the isoparametric mapping is not implemented, bits 16-31
!#   are unused!
!# </purpose>
!##############################################################################

MODULE transformation

  USE fsystem
  USE basicgeometry
  USE triangulation

  IMPLICIT NONE
  
!<constants>
!<constantblock description="Constants for size of auxiliary arrays.">

  ! Number of entries in the array with the auxiliary Jacobian factors in 2D
  INTEGER, PARAMETER :: TRAFO_NAUXJAC2D = 4

!</constantblock>


!<constantblock description="Id values for bit 0..7 of the transformation ID.">

  ! Unspecified transformation
  INTEGER, PARAMETER :: TRAFO_ID_UNKNOWN       = 0

  ! Linear transformation for simplex-type elements (1D-lines, 2D-triangles, 
  ! 3D-tethrahedrons)
  INTEGER, PARAMETER :: TRAFO_ID_LINSIMPLEX    = 1

  ! Multilinear (Bilinear/Trilinear) transformation for cubic-shaped elements 
  ! (2D-quadrilaterals, 3D-hexahedrals)
  INTEGER, PARAMETER :: TRAFO_ID_MLINCUBE     = 2
  
  ! Quadratic transformation for simplex-type elements (1D-lines, 2D-triangles, 
  ! 3D-tethrahedrons)
  INTEGER, PARAMETER :: TRAFO_ID_QUADSIMPLEX   = 3

  ! Multiquadratic (Biquadratic/Triquadratic) quadrilateral transformation 
  ! for cubic-shaped elements (2D-quadrilaterals, 3D-hexahedrals)
  INTEGER, PARAMETER :: TRAFO_ID_MQUADCUBE    = 4
!</constantblock>


!<constantblock description="Dimension constants for bit 8+9 of the transformation ID.">
  ! Bitmasks for dimension
  INTEGER(I32), PARAMETER :: TRAFO_DIM_DIMENSION = 2**8 + 2**9

  ! 1D element
  INTEGER(I32), PARAMETER :: TRAFO_DIM_1D = ISHFT(NDIM1D,8)

  ! 2D element
  INTEGER(I32), PARAMETER :: TRAFO_DIM_2D = ISHFT(NDIM2D,8)
  
  ! 3D element
  INTEGER(I32), PARAMETER :: TRAFO_DIM_3D = ISHFT(NDIM3D,8)

  ! Bitmask for transformation ID, without additional information
  INTEGER(I32), PARAMETER :: TRAFO_DIM_IDMASK = 255 

  ! Bitmask for transformation ID including dimension, without additional information
  INTEGER(I32), PARAMETER :: TRAFO_DIM_IDDIMMASK = TRAFO_DIM_IDMASK + TRAFO_DIM_DIMENSION
!</constantblock>


!<constantblock description="id values for coordinate systems">
  ! Undefined coordinate system or no coordinate system
  INTEGER, PARAMETER :: TRAFO_CS_UNDEFINED   = 0

  ! Parameter value [-1..1] on reference interval in 1D
  INTEGER, PARAMETER :: TRAFO_CS_REF1D       = 1

  ! Barycentric coordinates on triangle
  INTEGER, PARAMETER :: TRAFO_CS_BARY2DTRI   = 2

  ! 2D coordinates on reference triangle
  INTEGER, PARAMETER :: TRAFO_CS_REF2DTRI    = 3

  ! 2D coordinates on real triangle
  INTEGER, PARAMETER :: TRAFO_CS_REAL2DTRI   = 4
  
  ! 2D coordinates on reference quadrilateral
  INTEGER, PARAMETER :: TRAFO_CS_REF2DQUAD   = 5

  ! 2D coordinates on real quadrilateral
  INTEGER, PARAMETER :: TRAFO_CS_REAL2DQUAD  = 6
!</constantblock>

!</constants>

CONTAINS

  ! ***************************************************************************

!<function>  

  ELEMENTAL INTEGER FUNCTION trafo_igetDimension(ctrafoType)

!<description>
  ! This function returns the dimensional constant that specifies which
  ! dimension (1D, 2D,...) an element uses.
!</description>

!<input>    
  ! ID of transformation to calculate
  INTEGER(I32), INTENT(IN) :: ctrafoType
!</input>

!<result>
  ! A constant that specifies the dimension of the transformation. NDIM2 for 2D,
  ! NDIM3D for 3D,...
!</result>

!</function>

    ! The dimension is encoded in two bits in the element quantifier!
    trafo_igetDimension = ISHFT(IAND(ctrafoType,TRAFO_DIM_DIMENSION),-8)

  END FUNCTION

  ! ***************************************************************************

!<function>  

  ELEMENTAL INTEGER FUNCTION trafo_igetNVE(ctrafoType)

!<description>
  ! This function returns for a given trynsformation ID the size NVE of the
  ! coordinate array which is used to hold all the coordinates of the vertices
  ! of an element (corners, edge midpoints,...).
!</description>

!<input>    
  ! The transformation code identifier.
  INTEGER(I32), INTENT(IN) :: ctrafoType
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
  SELECT CASE (trafo_igetDimension(ctrafoType))

  CASE (NDIM1D)
    ! 1D elements. Lines. Check the actual transformation
    ! ID how to transform.
  
    SELECT CASE (IAND(ctrafoType,TRAFO_DIM_IDMASK))
    
    CASE (TRAFO_ID_LINSIMPLEX)
      ! 1D simplex -> linear line transformation. 
      ! 2 DOF's in the transformation (given by the corners of the element)
      trafo_igetNVE = 2
      
    END SELECT
  
  CASE (NDIM2D)
    ! 2D elements. Triangles, Quadrilaterals. Check the actual transformation
    ! ID how to transform.
  
    SELECT CASE (IAND(ctrafoType,TRAFO_DIM_IDMASK))
    
    CASE (TRAFO_ID_LINSIMPLEX)
      ! 2D simplex -> linear triangular transformation. 
      ! 3 DOF's in the transformation (given by the corners of the element)
      trafo_igetNVE = 3
    
    CASE (TRAFO_ID_MLINCUBE)
      ! Bilinear transformation for cubic-shaped elements 
      ! -> Bilinear quadrilateral transformation.
      ! 4 DOF's in the transformation (given by the corners of the element)
      trafo_igetNVE = 4
    
    END SELECT
    
  CASE (NDIM3D)
    ! 3D elements. Tetrahedrals, Hexahedrals. Check the actual transformation
    ! ID how to transform.
  
    SELECT CASE (IAND(ctrafoType,TRAFO_DIM_IDMASK))
    
    CASE (TRAFO_ID_LINSIMPLEX)
      ! 3D simplex -> linear triangular transformation. 
      ! 4 DOF's in the transformation (given by the corners of the element)
      trafo_igetNVE = 4
    
    CASE (TRAFO_ID_MLINCUBE)
      ! Trilinear transformation for cubic-shaped elements 
      ! -> Trilinear hexahedral transformation.
      ! 8 DOF's in the transformation (given by the corners of the element)
      trafo_igetNVE = 8
    
    END SELECT
    
  END SELECT

  END FUNCTION

! **********************************************************************

!<subroutine>

  SUBROUTINE trafo_getCoords (ctrafoType,rtriangulation,iel,Dcoords)

!<description>
  ! This routine filld the Dcoords array with the coordinates (of corners,
  ! edge midpoints,...) of element iel in the triangulation rtriangulation.
  ! Dcoords then defines the actual transformation from the reference
  ! to the real element.
!</description>

!<input>
  ! ID of transformation 
  INTEGER(I32), INTENT(IN) :: ctrafoType

  ! Triangulation structure that defines the elements.
  TYPE(t_triangulation), INTENT(IN) :: rtriangulation
  
  ! Number of the element whose information is to be fetched.
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel
!</input>

!<output>
  ! Coordinates of the corners of the element
  !  Dcoord(1,i) = x-coordinates of corner i on an element, 
  !  Dcoord(2,i) = y-coordinates of corner i on an element.
  ! DIMENSION(#space dimensions,NVE).
  ! NVE depends on the type of transformation and can be determined with
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dcoords
!</output>

!</subroutine>

    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement

    ! What type of transformation do we have? First decide on the dimension,
    ! then on the actual ID.
    SELECT CASE (trafo_igetDimension(ctrafoType))
    
    CASE (NDIM1D)
      ! 1D elements. Lines
      ! We always need the corner-coordinate array, so get it 
      ! from the triangulation.
      CALL storage_getbase_double2d (rtriangulation%h_DvertexCoords,&  
                                     p_DvertexCoords)
      CALL storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&  
                                  p_IverticesAtElement)
      
      ! Check the actual transformation ID how to transform.
    
      SELECT CASE (IAND(ctrafoType,TRAFO_DIM_IDMASK))
      
      CASE (TRAFO_ID_LINSIMPLEX)
        ! 1D simplex -> linear line transformation. 
        ! Transfer the corners of the element.
        Dcoords (1,1:3) = p_DvertexCoords(1,&
                                p_IverticesAtElement(1:3,iel))
                                
      END SELECT
    
    CASE (NDIM2D)
      ! 2D elements. Triangles, Quadrilaterals. 
      ! We always need the corner-coordinate array, so get it 
      ! from the triangulation.
      CALL storage_getbase_double2d (rtriangulation%h_DvertexCoords,&  
                                     p_DvertexCoords)
      CALL storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&  
                                  p_IverticesAtElement)
      
      ! Check the actual transformation ID how to transform.
    
      SELECT CASE (IAND(ctrafoType,TRAFO_DIM_IDMASK))
      
      CASE (TRAFO_ID_LINSIMPLEX)
        ! 2D simplex -> linear triangular transformation. 
        ! Transfer the corners of the element.
        Dcoords (1:NDIM2D,1:3) = p_DvertexCoords(1:NDIM2D,&
                                    p_IverticesAtElement(1:3,iel))
      
      CASE (TRAFO_ID_MLINCUBE)
        ! Bilinear transformation for cubic-shaped elements 
        ! -> Bilinear quadrilateral transformation.
        ! Transfer the corners of the element.
        Dcoords (1:NDIM2D,1:4) = p_DvertexCoords(1:NDIM2D,&
                                    p_IverticesAtElement(1:4,iel))
      
      END SELECT

    CASE (NDIM3D)
      ! 3D elements. Tetrahedrals, Hexahedrals.
      ! We always need the corner-coordinate array, so get it 
      ! from the triangulation.
      CALL storage_getbase_double2d (rtriangulation%h_DvertexCoords,&  
                                     p_DvertexCoords)
      CALL storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&  
                                  p_IverticesAtElement)
      
      ! Check the actual transformation ID how to transform.
    
      SELECT CASE (IAND(ctrafoType,TRAFO_DIM_IDMASK))
      
      CASE (TRAFO_ID_LINSIMPLEX)
        ! 3D simplex -> linear tetrahedral transformation. 
        ! Transfer the corners of the element.
        Dcoords (1:NDIM3D,1:4) = p_DvertexCoords(1:NDIM3D,&
                                    p_IverticesAtElement(1:4,iel))
      
      CASE (TRAFO_ID_MLINCUBE)
        ! Trilinear transformation for cubic-shaped elements 
        ! -> Trilinear hexahedral transformation.
        ! Transfer the corners of the element.
        Dcoords (1:NDIM3D,1:8) = p_DvertexCoords(1:NDIM3D,&
                                    p_IverticesAtElement(1:8,iel))
      
      END SELECT
      
    END SELECT
      
  END SUBROUTINE

! **********************************************************************

!<subroutine>

  SUBROUTINE trafo_getCoords_sim (ctrafoType,rtriangulation,Ielements,Dcoords)

!<description>
  ! This routine fills the Dcoords array with the coordinates (of corners,
  ! edge midpoints,...) of all the elements given in the element
  ! list Ielements, based on the triangulation rtriangulation.
  ! Dcoords then defines the actual transformation from the reference
  ! to the real elements.
!</description>

!<input>
  ! ID of transformation 
  INTEGER(I32), INTENT(IN) :: ctrafoType

  ! Triangulation structure that defines the elements.
  TYPE(t_triangulation), INTENT(IN) :: rtriangulation
  
  ! Array with element numbers whose coordinates should be extracted into
  ! Dcoords.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: Ielements
!</input>

!<output>
  ! Coordinates of the corners of all the elements
  !  Dcoord(1,i,.) = x-coordinates of corner i on an element, 
  !  Dcoord(2,i,.) = y-coordinates of corner i on an element.
  ! DIMENSION(#space dimensions,NVE,size of Ielements)
  ! NVE depends on the type of transformation and can be determined with
  REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: Dcoords
!</output>

!</subroutine>

    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER :: ipoint
    INTEGER(PREC_ELEMENTIDX) :: iel

    ! What type of transformation do we have? First decide on the dimension,
    ! then on the actual ID.
    SELECT CASE (trafo_igetDimension(ctrafoType))
    
    CASE (NDIM1D)
      ! 2D elements. Lines.
      ! We always need the corner-coordinate array, so get it 
      ! from the triangulation.
      CALL storage_getbase_double2d (rtriangulation%h_DvertexCoords,&  
                                     p_DvertexCoords)
      CALL storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&  
                                  p_IverticesAtElement)
      
      ! Check the actual transformation ID how to transform.
    
      SELECT CASE (IAND(ctrafoType,TRAFO_DIM_IDMASK))
      
      CASE (TRAFO_ID_LINSIMPLEX)
        ! 1D simplex -> linear line transformation. 
        ! Transfer the corners of the element.
        DO iel=1,SIZE(Ielements)
          DO ipoint = 1,2
            Dcoords (1,ipoint,iel) = &
              p_DvertexCoords(1,p_IverticesAtElement(ipoint,Ielements(iel)))
          END DO
        END DO
      
      END SELECT

    CASE (NDIM2D)
      ! 2D elements. Triangles, Quadrilaterals. 
      ! We always need the corner-coordinate array, so get it 
      ! from the triangulation.
      CALL storage_getbase_double2d (rtriangulation%h_DvertexCoords,&  
                                     p_DvertexCoords)
      CALL storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&  
                                  p_IverticesAtElement)
      
      ! Check the actual transformation ID how to transform.
    
      SELECT CASE (IAND(ctrafoType,TRAFO_DIM_IDMASK))
      
      CASE (TRAFO_ID_LINSIMPLEX)
        ! 2D simplex -> linear triangular transformation. 
        ! Transfer the corners of the element.
        DO iel=1,SIZE(Ielements)
          DO ipoint = 1,3
            Dcoords (1:NDIM2D,ipoint,iel) = &
              p_DvertexCoords(1:NDIM2D,p_IverticesAtElement(ipoint,Ielements(iel)))
          END DO
        END DO
      
      CASE (TRAFO_ID_MLINCUBE)
        ! Bilinear transformation for cubic-shaped elements 
        ! -> Bilinear quadrilateral transformation.
        ! Transfer the corners of the element.
        DO iel=1,SIZE(Ielements)
          DO ipoint = 1,4
            Dcoords (1:NDIM2D,ipoint,iel) = &
              p_DvertexCoords(1:NDIM2D,p_IverticesAtElement(ipoint,Ielements(iel)))
          END DO
        END DO
      
      END SELECT

    CASE (NDIM3D)
      ! 3D elements. Tetrahedrals, Hexahedrals.
      ! We always need the corner-coordinate array, so get it 
      ! from the triangulation.
      CALL storage_getbase_double2d (rtriangulation%h_DvertexCoords,&  
                                     p_DvertexCoords)
      CALL storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&  
                                  p_IverticesAtElement)
      
      ! Check the actual transformation ID how to transform.
    
      SELECT CASE (IAND(ctrafoType,TRAFO_DIM_IDMASK))
      
      CASE (TRAFO_ID_LINSIMPLEX)
        ! 3D simplex -> linear tetrahedral transformation. 
        ! Transfer the corners of the element.
        DO iel=1,SIZE(Ielements)
          DO ipoint = 1,4
            Dcoords (1:NDIM3D,ipoint,iel) = &
              p_DvertexCoords(1:NDIM3D,p_IverticesAtElement(ipoint,Ielements(iel)))
          END DO
        END DO
      
      CASE (TRAFO_ID_MLINCUBE)
        ! Trilinear transformation for cubic-shaped elements 
        ! -> Trilinear hexahedral transformation.
        ! Transfer the corners of the element.
        DO iel=1,SIZE(Ielements)
          DO ipoint = 1,8
            Dcoords (1:NDIM3D,ipoint,iel) = &
              p_DvertexCoords(1:NDIM3D,p_IverticesAtElement(ipoint,Ielements(iel)))
          END DO
        END DO
      
      END SELECT
      
    END SELECT
      
  END SUBROUTINE

! **********************************************************************

!<subroutine>

  SUBROUTINE trafo_calctrafo_mult (ctrafoType,npointsPerEl,Dcoords,&
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
  INTEGER(I32), INTENT(IN) :: ctrafoType

  ! Number of points in each element where to calculate the transformation
  INTEGER, INTENT(IN) :: npointsPerEl

  ! Coordinates of the corners of all the elements
  !  Dcoord(1,i) = x-coordinates of corner i on an element, 
  !  Dcoord(2,i) = y-coordinates of corner i on an element.
  ! DIMENSION(#space dimensions,NVE)
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords

  ! Coordinates of the points on the reference element for each element 
  ! where to calculate the mapping.
  !  DpointsRef(1,i) = x-coordinates of point i on an element, 
  !  DpointsRef(2,i) = y-coordinates of point i on an element.
  ! ! DIMENSION(#space dimension,npointsPerEl)
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: DpointsRef
!</input>

!<output>
  ! The Jacobian matrix of the mapping for each point.
  ! DIMENSION(number of entries in the matrix,npointsPerEl)
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Djac
  
  ! Jacobian determinants of the mapping for all the points from the
  ! reference element to the real element.
  ! DIMENSION(npointsPerEl)
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Ddetj
  
  ! OPTIONAL: Array receiving the coordinates of the points in DpointsRef,
  ! mapped from the reference element to the real element.
  ! If not specified, they are not computed.
  ! DIMENSION(#space dimension,npointsPerEl)
  REAL(DP), DIMENSION(:,:), INTENT(OUT), OPTIONAL :: DpointsReal
!</output>

!</subroutine>

  ! local variables
  INTEGER :: ipt
  REAL(DP) :: dax, day, dbx, dby, dcx, dcy
  
  ! auxiliary factors for the bilinear quad mapping
  REAL(DP), DIMENSION(TRAFO_NAUXJAC2D) :: DjacPrep
  
  ! What type of transformation do we have? First decide on the dimension,
  ! then on the actual ID.
  SELECT CASE (trafo_igetDimension(ctrafoType))
  
  CASE (NDIM2D)
    ! 2D elements. Triangles, Quadrilaterals. Check the actual transformation
    ! ID how to transform.
  
    SELECT CASE (IAND(ctrafoType,TRAFO_DIM_IDMASK))
    
    CASE (TRAFO_ID_LINSIMPLEX)

      ! 2D simplex -> linear triangular transformation.
      !
      ! Calculate with or without coordinates?
      IF (.NOT. PRESENT(DpointsReal)) THEN
      
        ! Loop over the points
        DO ipt=1,npointsPerEl
        
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
          
        END DO ! ipt
        
    ELSE

        ! Loop over the points
        DO ipt=1,npointsPerEl
        
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
          
          ! Ok, that was easy. It's slightly more complicated to get
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
          
        END DO ! ipt
          
      END IF
          
    CASE (TRAFO_ID_MLINCUBE)
    
      ! Calculate with or without coordinates?
      IF (.NOT. PRESENT(DpointsReal)) THEN
      
        ! Prepare the calculation of the Jacobi determinants
        CALL trafo_calcJacPrepare(Dcoords, DjacPrep)
        
        ! Loop over the points
        DO ipt=1,npointsPerEl
          ! Calculate the Jacobian matrix and determinant
          CALL trafo_calcJac (Dcoords,DjacPrep,Djac(:,ipt),Ddetj(ipt), &
                               DpointsRef(1,ipt),DpointsRef(2,ipt))
        END DO ! ipt
      
      ELSE

        ! Prepare the calculation of the Jacobi determinants
        CALL trafo_calcJacPrepare(Dcoords, DjacPrep)
        
        ! Loop over the points
        DO ipt=1,npointsPerEl
          ! Calculate the Jacobian matrix and determinant
          CALL trafo_calcTrafo (Dcoords,DjacPrep,Djac(:,ipt),Ddetj(ipt), &
                                 DpointsRef(1,ipt),DpointsRef(2,ipt), &
                                 DpointsReal(1,ipt),DpointsReal(2,ipt))
        END DO ! ipt

      END IF
      
    END SELECT
    
  END SELECT

  END SUBROUTINE

! **********************************************************************

!<subroutine>

  SUBROUTINE trafo_calctrafo_sim (ctrafoType,nelements,npointsPerEl,Dcoords,&
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
  INTEGER(I32), INTENT(IN) :: ctrafoType

  ! Number of elements where to calculate the transformation
  INTEGER, INTENT(IN) :: nelements
  
  ! Number of points in each element where to calculate the transformation
  INTEGER, INTENT(IN) :: npointsPerEl

  ! Coordinates of the corners of all the elements
  !  Dcoord(1,i,.) = x-coordinates of corner i on an element, 
  !  Dcoord(2,i,.) = y-coordinates of corner i on an element.
  ! DIMENSION(#space dimensions,NVE,nelements)
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dcoords

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
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: DpointsRef
!</input>

!<output>
  ! The Jacobian matrix of the mapping for each point.
  ! DIMENSION(number of entries in the matrix,npointsPerEl,nelements)
  REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: Djac
  
  ! Jacobian determinants of the mapping for all the points from the
  ! reference element to the real element.
  ! DIMENSION(npointsPerEl,nelements)
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Ddetj
  
  ! OPTIONAL: Array receiving the coordinates of the points in DpointsRef,
  ! mapped from the reference element to the real element.
  ! If not specified, they are not computed.
  ! DIMENSION(#space dimensions,npointsPerEl,nelements)
  REAL(DP), DIMENSION(:,:,:), INTENT(OUT), OPTIONAL :: DpointsReal
!</output>

!</subroutine>

  ! local variables
  INTEGER :: iel, ipt
  REAL(DP) :: dax, day, dbx, dby, dcx, dcy
  
  ! auxiliary factors for the bilinear quad mapping
  REAL(DP), DIMENSION(TRAFO_NAUXJAC2D) :: DjacPrep
  
  ! What type of transformation do we have? First decide on the dimension,
  ! then on the actual ID.
  SELECT CASE (trafo_igetDimension(ctrafoType))
  
  CASE (NDIM2D)
    ! 2D elements. Triangles, Quadrilaterals. Check the actual transformation
    ! ID how to transform.
  
    SELECT CASE (IAND(ctrafoType,TRAFO_DIM_IDMASK))
    
    CASE (TRAFO_ID_LINSIMPLEX)
    
      ! 2D simplex -> linear triangular transformation.
      !
      ! Calculate with or without coordinates?
      IF (.NOT. PRESENT(DpointsReal)) THEN
      
        ! Loop over the elements
        DO iel = 1,nelements
          
          ! Loop over the points
          DO ipt=1,npointsPerEl
          
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
            
          END DO ! ipt
          
        END DO ! iel
        
      ELSE

        ! Loop over the elements
        DO iel = 1,nelements
          
          ! Loop over the points
          DO ipt=1,npointsPerEl
          
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
            
            ! Ok, that was easy. It's slightly more complicated to get
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
            
          END DO ! ipt
          
        END DO ! iel
        
      END IF
    
    CASE (TRAFO_ID_MLINCUBE)
    
      ! Bilinear transformation for cubic-shaped elements 
      ! -> Bilinear quadrilateral transformation.
      !
      ! Calculate with or without coordinates?
      IF (.NOT. PRESENT(DpointsReal)) THEN
      
        ! Loop over the elements
        DO iel = 1,nelements
          ! Prepare the calculation of the Jacobi determinants
          CALL trafo_calcJacPrepare(Dcoords(:,:,iel), DjacPrep)
          
          ! Loop over the points
          DO ipt=1,npointsPerEl
            ! Calculate the Jacobian matrix and determinant
            CALL trafo_calcJac (Dcoords(:,:,iel),DjacPrep,Djac(:,ipt,iel),Ddetj(ipt,iel), &
                                DpointsRef(1,ipt,iel),DpointsRef(2,ipt,iel))
          END DO ! ipt
          
        END DO ! iel
      
      ELSE

        ! Loop over the elements
        DO iel = 1,nelements
          ! Prepare the calculation of the Jacobi determinants
          CALL trafo_calcJacPrepare(Dcoords(:,:,iel), DjacPrep)
          
          ! Loop over the points
          DO ipt=1,npointsPerEl
            ! Calculate the Jacobian matrix and determinant
            CALL trafo_calcTrafo (Dcoords(:,:,iel),DjacPrep,Djac(:,ipt,iel),Ddetj(ipt,iel), &
                                  DpointsRef(1,ipt,iel),DpointsRef(2,ipt,iel), &
                                  DpointsReal(1,ipt,iel),DpointsReal(2,ipt,iel))
          END DO ! ipt
          
        END DO ! iel

      END IF
      
    END SELECT ! actual ID
  
  END SELECT ! Dimension

  END SUBROUTINE

! **********************************************************************
! Explaination of the quadrilateral transformation:
!
! We want to perform a transformation from the reference quadrilateral
! [-1,1]x[-1,1] onto a "real" quadrilaterl:
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
!  a1 - a2 - a3 - a4 = x1       b1 - b2 - b3 - b4 = y1
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
! They can be calculated with elem_calcJacPrepare in advance for all points that
! have to be mapped. To be more exact, elem_calcJacPrepare calculates:
!
!  J1 = 1/2 * (-x1 - x2 + x3 - x4)
!  J2 = 1/2 * ( x1 - x2 + x3 - x4)
!
!  J3 = 1/2 * (-y1 + y2 - y3 + y4)
!  J4 = 1/2 * (-y1 + y2 + y3 - y4)
!
! Using these factors, one can write:
!
!  a1  =  1/4 * ( x1 + x2 + x3 + x4)  =  1/2 * (x1 + x2 + J1)
!  a2  =  1/4 * (-x1 + x2 + x3 - x4)  =  1/2 * (x2 - x1 + J3)
!  a3  =  1/4 * (-x1 - x2 + x3 + x4)  =  1/2 * J1 
!  a4  =  1/4 * ( x1 - x2 + x3 - x4)  =  1/2 * J2
!
!  b1  =  1/4 * ( y1 + y2 + y3 + y4)  =  1/2 * (y1 + y3 + J2)
!  b2  =  1/4 * (-y1 + y2 + y3 - y4)  =  1/2 * J4
!  b3  =  1/4 * (-y1 - y2 + y3 + y4)  =  1/2 * (y3 - y1 - J4)
!  b4  =  1/4 * ( y1 - y2 + y3 - y4)  =  1/2 * J3
!
! The Jacobian matrix of the bilinear transformation is now 
! calculated as usual by partial differentiation. Thing above 
! coefficients ai and bi one can write:
!
! DPhi ( xi1 )  =  ( a2 + a4*xi2    a3 + a4*xi1 )
!      ( xi2 )     ( b2 + b4*xi2    b3 + b4*xi1 )
!
!               =  ( 1/2*(x2-x1+J2) + 1/2*J2*xi2           1/2*J1 + 1/2*J2*xi1 )
!                  (         1/2*J4 + 1/2*J3*xi2   1/2*(y3-y1-J4) + 1/2*J3*xi1 )
!
! which gives the Jacobian determinant of the 2x2-matrix:
!
!   det DPhi (xi1,xi2)  =  (DPhi[1,1]*DPhi[2,2] - DPhi[1,2]*DPhi[2,1]) (xi1,xi2)
!
! Using these information makes it possible to map a point (XI1,XI2)
! on the reference element to coordinates (XX,YY) on the real element by:
!
!   (XX,YY) := Psi(XI1,XI2) 
! **********************************************************************

!<subroutine>

  PURE SUBROUTINE trafo_calcJacPrepare(Dcoords, DjacPrep)
  
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
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
!</input>
  
!<output>
  ! Arrays with auxiliary jacobian factors for later computation
  ! DIMENSION(TRAFO_NAUXJAC2D)
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DjacPrep
!</output>

!</subroutine>

  DjacPrep(1) = 0.5_DP * (-Dcoords(1,1) - Dcoords(1,2) + Dcoords(1,3) + Dcoords(1,4))
  DjacPrep(2) = 0.5_DP * ( Dcoords(1,1) - Dcoords(1,2) + Dcoords(1,3) - Dcoords(1,4))
  DjacPrep(3) = 0.5_DP * (-Dcoords(2,1) + Dcoords(2,2) - Dcoords(2,3) + Dcoords(2,4))
  DjacPrep(4) = 0.5_DP * (-Dcoords(2,1) + Dcoords(2,2) + Dcoords(2,3) - Dcoords(2,4))

  END SUBROUTINE 

  !************************************************************************

!<subroutine>

  PURE SUBROUTINE trafo_calcTrafo (Dcoord,DjacPrep,Djac,ddetj, &
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
  ! have to be calculated with elem_calcJacPrepare for the 
  ! considered element.
!</description>  

!<input>
  
  ! Coordinates of the four points forming the element.
  ! Dcoord(1,.) = x-coordinates,
  ! Dcoord(2,.) = y-coordinates.
  ! DIMENSION(#space dimensions,NVE)
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoord 
  
  ! Auxiliary constants for the considered element with
  ! coordinates in Dcoord; have to be computed previously
  ! by elem_calcJacPrepare.
  ! DIMENSION(TRAFO_NAUXJAC2D)
  REAL(DP), DIMENSION(:), INTENT(IN) :: DjacPrep
  
  ! Coordinates of a point on the reference element
  REAL(DP), INTENT(IN) :: dparx,dpary
  
!</input>
  
!<output>
  ! The Jacobian matrix of the mapping from the reference to the 
  ! real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Djac
  
  ! The determinant of the mapping.
  REAL(DP), INTENT(OUT) :: ddetj
  
  ! X-Coordinates of a point on the real element
  REAL(DP), INTENT(OUT) :: dxreal
  
  ! X-Coordinates of a point on the real element
  REAL(DP), INTENT(OUT) :: dyreal
!</output>
  
!</subroutine>

  Djac(1) = 0.5E0_DP * ((Dcoord(1,2)-Dcoord(1,1)+DjacPrep(2)) &
                          + DjacPrep(2)*dpary )
  Djac(2) = 0.5E0_DP * ( DjacPrep(4) - DjacPrep(3)*dpary )
  Djac(3) = 0.5E0_DP * ( DjacPrep(1) + DjacPrep(2)*dparx )
  Djac(4) = 0.5E0_DP * ((Dcoord(2,3)-Dcoord(2,1)-DjacPrep(4)) &
                          - DjacPrep(3)*dparx )

  ! Determinant of the mapping

  ddetj = Djac(1)*Djac(4) - Djac(3)*Djac(2)
  
  ! Map the point to the real element
  
  dxreal = 0.5E0_DP*((Dcoord(1,1)+Dcoord(1,2)+DjacPrep(1)) + &
                     (Dcoord(1,2)-Dcoord(1,1)+DjacPrep(2))*dparx + &
                     DjacPrep(1)*dpary + &
                     DjacPrep(2)*dparx*dpary )
  dyreal = 0.5E0_DP*((Dcoord(2,1)+Dcoord(2,3)+DjacPrep(3)) + &
                     DjacPrep(4)*dparx + &
                     (Dcoord(2,3)-Dcoord(2,1)-DjacPrep(4))*dpary - &
                     DjacPrep(3)*dparx*dpary )
  
  END SUBROUTINE

  !************************************************************************

!<subroutine>

  PURE SUBROUTINE trafo_calcJac (Dcoord,DjacPrep,Djac,ddetj,dparx,dpary)

!<description>
  ! Calculate Jacobian determinant of mapping from reference- to
  ! real element.
  !
  ! This routine only calculates the Jacobian determinant of the
  ! mapping. This is on contrast to elem_calcRealCoords, which not only 
  ! calculates this determinant but also maps the point. So this routine
  ! can be used to speed up the code if the coordinates of the
  ! mapped point already exist.
  !
  ! Before this routine can be called, the auxiliary factors DJF
  ! have to be calculated with elem_calcJacPrepare for the considered element.
!</description>  

!<input>
  ! Coordinates of the four points forming the element.
  ! Dcoord(1,.) = x-coordinates,
  ! Dcoord(2,.) = y-coordinates.
  ! DIMENSION(#space dimensions,NVE)
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoord
  
  ! Auxiliary constants for the considered element with
  ! coordinates in Dcoord; have to be computed previously
  ! by elem_calcJacPrepare.
  ! DIMENSION(TRAFO_NAUXJAC2D)
  REAL(DP), DIMENSION(:), INTENT(IN) :: DjacPrep
  
  ! Coordinates of a point on the reference element
  REAL(DP), INTENT(IN) :: dparx,dpary
!</input>
  
!<output>
  ! The Jacobian matrix of the mapping from the reference to the 
  ! real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Djac
  
  ! The determinant of the mapping.
  REAL(DP), INTENT(OUT) :: ddetj
!</output>
  
!</subroutine>

  ! Jacobian matrix of the mapping:

  Djac(1) = 0.5E0_DP * (Dcoord(1,2)-Dcoord(1,1)+DjacPrep(2)) &
            + 0.5E0_DP * DjacPrep(2)*dpary
  Djac(2) = 0.5E0_DP * DjacPrep(4) - 0.5E0_DP * DjacPrep(3)*dpary
  Djac(3) = 0.5E0_DP * DjacPrep(1) + 0.5E0_DP * DjacPrep(2)*dparx
  Djac(4) = 0.5E0_DP * (Dcoord(2,3)-Dcoord(2,1)-DjacPrep(4)) &
            - 0.5E0_DP * DjacPrep(3)*dparx

  ! Determinant of the mapping

  ddetj = Djac(1)*Djac(4) - Djac(3)*Djac(2)
  
  END SUBROUTINE

  !************************************************************************

!<subroutine>

  PURE SUBROUTINE trafo_calcRealCoords2 (Dcoord,DjacPrep,dparx,dpary,dxreal,dyreal)

!<description>
  ! This subroutine computes the real coordinates of a point which 
  ! is given by parameter values (dparx,dpary) on the reference element. 
  ! It is assumed that the array DjacPrep was initialised before using
  ! with elem_calcJacPrepare.
!</description>  

!<input>
  ! Coordinates of the four points forming the element.
  ! Dcoord(1,.) = x-coordinates,
  ! Dcoord(2,.) = y-coordinates.
  ! DIMENSION(#space dimensions,NVE)
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoord
  
  ! Auxiliary constants for the considered element with
  ! coordinates in Dcoord; have to be computed previously
  ! by elem_calcJacPrepare.
  ! DIMENSION(TRAFO_NAUXJAC2D)
  REAL(DP), DIMENSION(:), INTENT(IN) :: DjacPrep
  
  ! Coordinates of a point on the reference element
  REAL(DP), INTENT(IN) :: dparx,dpary
!</input>
  
!<output>
  ! X-Coordinates of a point on the real element
  REAL(DP), INTENT(OUT) :: dxreal
  
  ! X-Coordinates of a point on the real element
  REAL(DP), INTENT(OUT) :: dyreal
!</output>
  
!</subroutine>

  ! Map the point to the real element
  
  dxreal = 0.5E0_DP*((Dcoord(1,1)+Dcoord(1,2)+DjacPrep(1)) + &
                     (Dcoord(1,2)-Dcoord(1,1)+DjacPrep(2))*dparx + &
                     DjacPrep(1)*dpary + &
                     DjacPrep(2)*dparx*dpary )
  dyreal = 0.5E0_DP*((Dcoord(2,1)+Dcoord(2,3)+DjacPrep(3)) + &
                     DjacPrep(4)*dparx + &
                     (Dcoord(2,3)-Dcoord(2,1)-DjacPrep(4))*dpary - &
                     DjacPrep(3)*dparx*dpary )
  
  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  PURE SUBROUTINE trafo_mapCubPts1Dto2DRefQuad(iedge, ncubp, Dxi1D, Dxi2D)

!<description>
  ! This routine maps the coordinates of the cubature points in 1D
  ! (given by Dxi1D(1..ncubp,1)) to an edge on the reference quadrilateral
  ! in 2D. iedge specifies the edge where to map the cubature points to.
!</description>

!<input>
  ! Number of the local edge of the element where to map the points to.
  INTEGER, INTENT(IN) :: iedge

  ! number of cubature points
  INTEGER , INTENT(IN) :: ncubp
  
  ! Cubature point coordinates on 1D reference interval [-1,1]
  !     Dxi(1..ncubp,1)=coordinates
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dxi1D
!</input>
  
!<output>
  ! Coordinates of the cubature points on the edge in 2D.
  !        Dxi2D(1..ncubp,1)=x-coord, 
  !        Dxi2D(1..ncubp,2)=y-coord
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dxi2D
!</output>

!</subroutine>    

    ! local variables
    INTEGER :: ii

    ! We have to transfer
    ! the coordinates of the cubature points from 1D to 2D depending
    ! on this edge.

    IF (iedge .EQ. 1) THEN
      ! Edge 1 is on the bottom of the reference element
      DO ii = 1,ncubp
        Dxi2D(ii,1) = Dxi1D(ii,1)
        Dxi2D(ii,2) = -1.0_DP
      END DO
    ELSE IF (iedge .EQ. 2) THEN
      ! Edge 2 is on the right of the reference element
      DO ii = 1,ncubp
        Dxi2D(ii,1) = 1.0_DP
        Dxi2D(ii,2) = Dxi1D(ii,1)
      END DO
    ELSE IF (iedge .EQ. 3) THEN
      ! Edge 3 is on the top of the reference element
      DO ii = 1,ncubp
        Dxi2D(ii,1) = -Dxi1D(ii,1)
        Dxi2D(ii,2) = 1.0_DP
      END DO
    ELSE 
      ! Edge 4 is on the left of the reference element
      DO ii = 1,ncubp
        Dxi2D(ii,1) = -1.0_DP
        Dxi2D(ii,2) = -Dxi1D(ii,1)
      END DO
    END IF
  
  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  PURE SUBROUTINE trafo_mapCubPts1Dto2DTriBary(iedge, ncubp, Dxi1D, Dxi2D)

!<description>
  ! This routine maps the coordinates of the cubature points in 1D
  ! (given by Dxi1D(1..ncubp,1)) to an edge on a 2D triangle
  ! in barycentric coordinates. 
  ! iedge specifies the edge where to map the cubature points to.
!</description>

!<input>
  ! Number of the local edge of the element where to map the points to.
  INTEGER, INTENT(IN) :: iedge

  ! number of cubature points
  INTEGER , INTENT(IN) :: ncubp
  
  ! Cubature point coordinates on 1D reference interval [-1,1]
  !     Dxi(1..ncubp,1)=coordinates
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dxi1D
!</input>
  
!<output>
  ! Coordinates of the cubature points on the edge in 2D.
  !        Dxi2D(1..ncubp,1)=1st barycentric coordinate, 
  !        Dxi2D(1..ncubp,2)=2nd barycentric coordinate,
  !        Dxi2D(1..ncubp,3)=3rd barycentric coordinate
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dxi2D
!</output>

!</subroutine>    

    ! local variables
    INTEGER :: ii

    ! We have to transfer
    ! the coordinates of the cubature points from 1D to 2D depending
    ! on this edge.

    IF (iedge .EQ. 1) THEN
      ! Edge 1 is between 1st and 2nd point
      DO ii = 1,ncubp
        Dxi2D(ii,1) = 0.5_DP*(1.0_DP-Dxi1D(ii,1))
        Dxi2D(ii,2) = 0.5_DP*(1.0_DP+Dxi1D(ii,1))
        Dxi2D(ii,3) = 0.0_DP
      END DO
    ELSE IF (iedge .EQ. 2) THEN
      ! Edge 2 is between 2nd and 2rd point
      DO ii = 1,ncubp
        Dxi2D(ii,1) = 0.0_DP
        Dxi2D(ii,2) = 0.5_DP*(1.0_DP-Dxi1D(ii,1))
        Dxi2D(ii,3) = 0.5_DP*(1.0_DP+Dxi1D(ii,1))
      END DO
    ELSE 
      ! Edge 3 is between 3rd and 1st point
      DO ii = 1,ncubp
        Dxi2D(ii,1) = 0.5_DP*(1.0_DP+Dxi1D(ii,1))
        Dxi2D(ii,2) = 0.0_DP
        Dxi2D(ii,3) = 0.5_DP*(1.0_DP-Dxi1D(ii,1))
      END DO
    END IF
  
  END SUBROUTINE

  !************************************************************************

!<subroutine>

  PURE SUBROUTINE trafo_calcBackTrafo (Dcoord, &
                                       dparx,dpary,dxreal,dyreal)

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
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoord 
  
  ! Coordinates of a point on the real element
  REAL(DP), INTENT(IN) :: dxreal,dyreal
!</input>
  
!<output>
  ! Coordinates of a point on the reference element
  REAL(DP), INTENT(OUT) :: dparx,dpary
!</output>
  
!</subroutine>

    ! local variables

    REAL(DP) :: X1,X2,X3,X4,Y1,Y2,Y3,Y4,XB,YB,XCX,YCX,XCE,YCE
    REAL(DP) :: A,J1,J2,X0,Y0,BXI,BETA,CXI,XP0,YP0,CETA
    REAL(DP) :: ROOT1,ROOT2,XIP1,XIP2,ETAP1,ETAP2,D1,D2

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

    ROOT1 = -SQRT(BXI**2-2.0*J1*CXI)-BXI
    ROOT2 =  SQRT(BXI**2-2.0*J1*CXI)-BXI
    IF (ROOT1 .NE. 0.0_DP) THEN
      XIP1 = 2.0*CXI/ROOT1
    ELSE
      XIP1 = 1E15_DP
    END IF
    IF (ROOT2 .NE. 0.0_DP) THEN
      XIP2 = 2.0_DP*CXI/ROOT2
    ELSE
      XIP2 = 1E15_DP
    END IF

    ROOT1 =  SQRT(BETA**2+2.0_DP*J2*CETA)-BETA
    ROOT2 = -SQRT(BETA**2+2.0_DP*J2*CETA)-BETA
    IF (ROOT1 .NE. 0.0_DP) THEN
      ETAP1 = 2.0_DP*CETA/ROOT1
    ELSE
      ETAP1 = 1E15_DP
    END IF
    IF (ROOT2 .NE. 0.0_DP) THEN
      ETAP2 = 2.0_DP*CETA/ROOT2
    ELSE
      ETAP2 = 1E15_DP
    END IF

    D1 = SQRT(XIP1**2+ETAP1**2)
    D2 = SQRT(XIP2**2+ETAP2**2)

    IF (D1 .LT. D2) THEN
      dparx = XIP1
      dpary = ETAP1
    ELSE
      dparx = XIP2
      dpary = ETAP2
    END IF
  
  END SUBROUTINE

END MODULE 
