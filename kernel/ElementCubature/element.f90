!##############################################################################
!# ****************************************************************************
!# <name> element </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to evaluate a finite element on an
!# element primitive.
!#
!# Each element is characterised by four things:\\
!# - an element type identifier (EL_xxxx constant)\\
!# - elem_XXXX
!#    -> an evaluation function for that element in one point\\
!# - elem_XXXX_mult
!#    -> an evaluation function for that element in multiple points
!#       (e.g. a set of cubature points).\\
!# - elem_XXXX_sim
!#    -> an evaluation function for that element in multiple points
!#       (e.g. a set of cubature points) on multiple elements simultaneously.\\
!#
!# The evaluation function for multiple points accepts a set of point 
!# coordinates where to evaluate and evaluates the finite element in 
!# all these points. 
!# The evaluation functions for multiple elements work like the variant
!# for multiple points, but also accept a list of elements where to
!# evaluate.
!# Depending on the type of the element, the coordinates must be given 
!# to the evaluation functions either on the reference element (for 
!# parametric elements) or on the real element (for nonparametric elements).\\
!#
!# There is also a couple of auxiliary routines which help to deal with
!# a finite element:\\
!#
!# 1.) elem_igetNDofLoc
!#     -> determines the number of local degrees of freedom for a finite 
!#        element\\
!#
!# 2.) elem_igetNVE
!#     -> get the number of vertices in the element primitive/element shape
!#        of the Finite Element (3=triangular, 4=quad)\\
!#
!# 3.) elem_igetCoordSystem
!#     -> get the type of coordinate system, a Finite Element uses\\
!#
!# 4.) elem_igetTrafoType
!#     -> Determine the type of transformation from the reference element
!#        to the real element.\\
!#
!# 5.) elem_isnonparametric
!#     -> Check whether an element is parametric or nonparametric\\
!#
!# 6.) elem_generic 
!#     -> Realises a generic element which can be used to evaluate a finite 
!#        element depending on its element identifier - in contrast to the 
!#        standard evaluation routines, which ignore the element quantifier 
!#        as they 'know' what they are...\\
!#
!# 7.) elem_generic_mult
!#     -> The multiple-point-evaluation routine for a generic element.\\
!#
!# 8.) elem_generic_sim
!#     -> The multiple-point/element-evaluation routine for a generic element.
!#
!#
!#
!#  FAQ - Some explainations
!# --------------------------
!# 1.) What is an element and what are Degrees Of Freedoms (DOF's)?
!#
!#   Imagine a triangulation of a domain with quadrilateral polygons, e.g.
!#   rectangles or even squares. Take for example four of these squares and
!#   form a large square, called a "patch", like
!#
!#       O-----O-----O
!#       |     |     |
!#       |     |     |
!#       O-----1-----O
!#       |     |     |
!#       |     |     |
!#       O-----O-----O
!#
!#   A point-based finite element basis function like $Q_1$ is a function
!#   that is =1 at the center of this patch and =0 in all the other corner
!#   points. On each of the squares, the function is a polynom, which
!#   coefficients are chosen in such a way, that it is =1 in the center and
!#   =0 in the other corner points. Outside of this patch, the function
!#   is identically =0.
!#
!#   If you call one such a function $phi_i$, the whole function 
!#   $u:R^2 \to R$ is then defined as a sum
!#
!#     $$ u(x,y) = \sum_i  u_i  phi_i(x,y) $$
!#
!#   where the coefficients $u_i$ can be modified to change the whole
!#   function. This set ${u_i}$ are called "Degrees of Freedom" of the FE
!#   function $u$. If the $\phi_i$ are chosen to be =1 in the center
!#   and =0 in the other corners (like in $P_1$ and $Q_1$, by "lucky 
!#   chance" the $u_i$ are exactly the values of $u$ in these points --
!#   $u_i = u(x_i,y_i)$.
!#
!#   An "element" in the sense of this library focuses on one single
!#   polygon, not on the whole patch.
!#
!# </purpose>
!##############################################################################

MODULE element

  USE fsystem
  USE basicgeometry
  USE derivatives
  USE transformation

  IMPLICIT NONE
  
!<constants>
!<constantblock description="Element identifiers for 2D elements">

  ! unspecified element
  INTEGER, PARAMETER :: EL_UNDEFINED = -1

  ! ID of bilinear conforming triangular FE, P0 (just for the FEAST-users...)
  INTEGER, PARAMETER :: EL_P0   = 0

  ! ID of bilinear conforming triangular FE, P0
  INTEGER, PARAMETER :: EL_E000 = 0

  ! ID of bilinear conforming triangular FE, P1 (just for the FEAST-users...)
  INTEGER, PARAMETER :: EL_P1   = 1

  ! ID of bilinear conforming triangular FE, P1 
  INTEGER, PARAMETER :: EL_E001 = 1

  ! ID of bilinear conforming triangular FE, P2 (just for the FEAST-users...)
  INTEGER, PARAMETER :: EL_P2   = 2

  ! ID of bilinear conforming triangular FE, P2 
  INTEGER, PARAMETER :: EL_E002 = 2

  ! ID of bilinear conforming triangular FE, P3 (just for the FEAST-users...)
  INTEGER, PARAMETER :: EL_P3   = 3

  ! ID of bilinear conforming triangular FE, P3
  INTEGER, PARAMETER :: EL_E003 = 3

  ! ID of bilinear conforming triangular FE, Q0 (just for the FEAST-users...)
  INTEGER, PARAMETER :: EL_Q0   = 10

  ! ID of bilinear conforming triangular FE, Q0
  INTEGER, PARAMETER :: EL_E010 = 10

  ! ID of bilinear conforming quadrilateral FE, Q1 (just for the FEAST-users...)
  INTEGER, PARAMETER :: EL_Q1   = 11

  ! ID of bilinear conforming quadrilateral FE, Q1 
  INTEGER, PARAMETER :: EL_E011 = 11 

  ! ID of biquadratic conforming quadrilateral FE, Q2 (just for the FEAST-users...)
  INTEGER, PARAMETER :: EL_Q2   = 13

  ! ID of biquadratic conforming quadrilateral FE, Q2 
  INTEGER, PARAMETER :: EL_E013 = 13

  ! ID of bicubic conforming quadrilateral FE, Q3 (just for the FEAT-users...)
  INTEGER, PARAMETER :: EL_Q3   = 14

  ! ID of biquadratic conforming quadrilateral FE, Q3 
  INTEGER, PARAMETER :: EL_E014 = 14

  ! ID of biquadratic conforming quadrilateral FE, Q1~, integral
  ! mean value based
  INTEGER, PARAMETER :: EL_E030 = 30

  ! ID of biquadratic conforming quadrilateral FE, Q1~, edge-midpoint based
  INTEGER, PARAMETER :: EL_E031 = 31

  ! ID of biquadratic nonconforming quadrilateral FE, Q1~, integral
  ! mean value based
  INTEGER, PARAMETER :: EL_EM30 = -30

  ! ID of biquadratic nonconforming quadrilateral FE, Q1~, edge-midpoint based
  INTEGER, PARAMETER :: EL_EM31 = -31

!</constantblock>

!<constantblock description="maximal values">

  ! Maximum number of basic functions = maximum number of
  ! local DOF's per element.
  INTEGER, PARAMETER :: EL_MAXNBAS = 21

  ! (maximum) number of vertices per element
  INTEGER, PARAMETER :: EL_MAXNVE = 4
  
  ! Maximum different derivative descriptor types supported by elements;
  ! corresponds to the maximum value of DER_xxxx constants in the 
  ! module 'derivatives'. Can be used as first quantifier in the
  ! array given to evaluation subroutines which receives function values/
  ! derivatives/...
  INTEGER, PARAMETER :: EL_MAXNDER = 6
  
  INTEGER, PARAMETER :: EL_MAXNCOF = 6
  
  ! Number of entries in the Jacobian matrix, defining the mapping between
  ! the reference element and the real element. 2x2-matrix=4 elements
  INTEGER, PARAMETER :: EL_NJACENTRIES2D = 4

  ! Number of entries in the array with the auxiliary Jacobian factors
  INTEGER, PARAMETER :: EL_NAUXJAC2D = 4

!</constantblock>

!</constants>

CONTAINS

  ! ***************************************************************************

!<function>  

  ELEMENTAL INTEGER(I32) FUNCTION elem_igetNDofLoc(ieltype)

!<description>
  ! This function returns for a given element type the number of local
  ! degrees of freedom.
!</description>

!<input>    

  ! The element type identifier.
  INTEGER, INTENT(IN) :: ieltype

!</input>

!<result>
  ! The number of local DOF's on this element.
!</result>

!</function>

  SELECT CASE (ieltype)
  CASE (EL_P0, EL_Q0)
    ! local DOF's for Q0
    elem_igetNDofLoc = 1
  CASE (EL_P1)
    ! local DOF's for P1
    elem_igetNDofLoc = 3
  CASE (EL_P2)
    ! local DOF's for P2
    elem_igetNDofLoc = 6
  CASE (EL_P3)
    ! local DOF's for P3
    elem_igetNDofLoc = 9
  CASE (EL_Q1)
    ! local DOF's for Q1
    elem_igetNDofLoc = 4
  CASE (EL_Q2)
    ! local DOF's for Q2
    elem_igetNDofLoc = 9
  CASE (EL_Q3)
    ! local DOF's for Q3
    elem_igetNDofLoc = 16
  CASE (EL_E030, EL_E031, EL_EM30, EL_EM31)
    ! local DOF's for Ex30
    elem_igetNDofLoc = 4
  CASE DEFAULT
    elem_igetNDofLoc = 0
  END SELECT

  END FUNCTION

  ! ***************************************************************************

!<function>  

  PURE SUBROUTINE elem_igetNDofLocAssignment(ieltype, &
      ndofAtVertices, ndofAtEdges, ndofAtElement)

!<description>
  ! This function returns for a given element type the number of local
  ! degrees of freedom that is assigned to vertices, edges and elements.
!</description>

!<input>    
  ! The element type identifier.
  INTEGER, INTENT(IN) :: ieltype
!</input>

!<output>
  ! Number of DOF's assigned to the vertices on one element.
  INTEGER, INTENT(OUT) :: ndofAtVertices
  
  ! Number of DOF's assigned to the edges on one element.
  INTEGER, INTENT(OUT) :: ndofAtEdges
  
  ! Number of DOF's assigned to one element, which do not belong to
  ! vertices or edges.
  INTEGER, INTENT(OUT) :: ndofAtElement

!</output>

!</function>

  SELECT CASE (ieltype)
  CASE (EL_P0, EL_Q0)
    ! local DOF's for Q0
    ndofAtVertices = 0
    ndofAtEdges    = 0
    ndofAtElement  = 1
  CASE (EL_P1)
    ! local DOF's for P1
    ndofAtVertices = 3
    ndofAtEdges    = 0
    ndofAtElement  = 0
  CASE (EL_P2)
    ! local DOF's for P2
    ndofAtVertices = 3
    ndofAtEdges    = 3
    ndofAtElement  = 0
  CASE (EL_P3)
    ! local DOF's for P3
    ndofAtVertices = 3
    ndofAtEdges    = 6
    ndofAtElement  = 0
  CASE (EL_Q1)
    ! local DOF's for Q1
    ndofAtVertices = 4
    ndofAtEdges    = 0
    ndofAtElement  = 0
  CASE (EL_Q2)
    ! local DOF's for Q2
    ndofAtVertices = 4
    ndofAtEdges    = 4
    ndofAtElement  = 1
  CASE (EL_Q3)
    ! local DOF's for Q3
    ndofAtVertices = 4
    ndofAtEdges    = 8
    ndofAtElement  = 4
  CASE (EL_E030, EL_E031, EL_EM30, EL_EM31)
    ! local DOF's for Ex30
    ndofAtVertices = 0
    ndofAtEdges    = 4
    ndofAtElement  = 0
  CASE DEFAULT
    ndofAtVertices = 0
    ndofAtEdges    = 0
    ndofAtElement  = 0
  END SELECT

  END SUBROUTINE

  ! ***************************************************************************

!<function>  

  ELEMENTAL INTEGER(I32) FUNCTION elem_igetNVE(ieltype)

!<description>
  ! This function returns for a given element type the number of vertices
  ! in the corresponding element primitive, i.e. 3 for triangular and
  ! 4 for quadrilateral elements.
!</description>

!<input>    

  ! The element type identifier.
  INTEGER, INTENT(IN) :: ieltype

!</input>

!<result>
  ! The number vertices in the element primitive/shape where the finite
  ! element is defined on.
!</result>

!</function>

  SELECT CASE (ieltype)
  CASE (EL_P0,EL_P1,EL_P2,EL_P3)
    elem_igetNVE = 3
  CASE (EL_Q0,EL_Q1,EL_Q2,EL_Q3,EL_E030, EL_E031, EL_EM30, EL_EM31)
    elem_igetNVE = 4
  CASE DEFAULT
    elem_igetNVE = 0
  END SELECT

  END FUNCTION

  ! ***************************************************************************

!<function>  

  ELEMENTAL INTEGER(I32) FUNCTION elem_igetCoordSystem(ieltype)

!<description>
  ! This function returns for a given element type the type of the coordinate
  ! system. Whenever an element of this type is evaluated, it expects the
  ! coordinates where to evaluate in the coordinates of this coordinate
  ! system. For example, triangular $P_1$ elements usually work in the
  ! barycentric triangular coordinate system, identified by the coordinate
  ! system identifier TRAFO_CS_BARY2DTRI.
!</description>

!<input>    

  ! The element type identifier.
  INTEGER, INTENT(IN) :: ieltype

!</input>

!<result>
  ! The type of the coordinate system the element ieltype acts on. One of
  ! the TRAFO_CS_xxxx constants from the transformation module.
!</result>

!</function>

  SELECT CASE (ieltype)
  CASE (EL_P0,EL_P1,EL_P2,EL_P3)
    ! Triangular elements work in barycentric coordinates
    elem_igetCoordSystem = TRAFO_CS_BARY2DTRI
  CASE (EL_Q0,EL_Q1,EL_Q2,EL_Q3,EL_E030, EL_E031)
    ! These work on the reference quadrilateral
    elem_igetCoordSystem = TRAFO_CS_REF2DQUAD
  CASE (EL_EM30, EL_EM31)
    ! These work in real coordinates
    elem_igetCoordSystem = TRAFO_CS_REAL2DQUAD
  CASE DEFAULT
    elem_igetCoordSystem = TRAFO_CS_UNDEFINED
  END SELECT
  
  END FUNCTION

  ! ***************************************************************************

!<function>  

  ELEMENTAL INTEGER FUNCTION elem_igetTrafoType(ieltype)

!<description>
  ! This function returns the typical type of transformation, a specific 
  ! element uses for transferring coordinates of the reference element to 
  ! the real element.
!</description>

!<input>    

  ! The element type identifier.
  INTEGER, INTENT(IN) :: ieltype

!</input>

!<result>
  ! The maximum derivative of that element. <= MAX_NDER.
!</result>

!</function>

  SELECT CASE (ieltype)
  CASE (EL_P0, EL_P1, EL_P2, EL_P3)
    ! Linear triangular transformation
    elem_igetTrafoType = TRAFO_IDLINTRI
  CASE (EL_Q0,EL_Q1,EL_Q2,EL_Q3,EL_E030, EL_E031, EL_EM30, EL_EM31)
    ! Bilinear quadrilateral transformation
    elem_igetTrafoType = TRAFO_IDBILINQUAD
  CASE DEFAULT
    elem_igetTrafoType = TRAFO_IDUNKNOWN
  END SELECT

  END FUNCTION

  ! ***************************************************************************

!<function>  

  ELEMENTAL INTEGER FUNCTION elem_getMaxDerivative(ieltype)

!<description>
  ! This function determines for a given element type the maximum derivative
  ! that the element can calculate (One of the DER_xxxx constants and 
  ! <= MAX_NDER). This can be used to specify the size of the array DBas.
!</description>

!<input>    

  ! The element type identifier.
  INTEGER, INTENT(IN) :: ieltype

!</input>

!<result>
  ! One of the TRAFO_IDxxxx constants of the module 'transformation' identifying
  ! the type of transformation.
!</result>

!</function>

  SELECT CASE (ieltype)
  CASE (EL_P0, EL_Q0)
    ! Function + 1st derivative
    elem_getMaxDerivative = 1
  CASE (EL_P1)
    ! Function + 1st derivative
    elem_getMaxDerivative = 3
  CASE (EL_P2)
    ! Function + 1st derivative
    elem_getMaxDerivative = 3
  CASE (EL_P3)
    ! Function + 1st derivative
    elem_getMaxDerivative = 3
  CASE (EL_Q1)
    ! Function + 1st derivative
    elem_getMaxDerivative = 3
  CASE (EL_Q2)
    ! Function + 1st derivative
    elem_getMaxDerivative = 3
  CASE (EL_Q3)
    ! Function + 1st derivative
    elem_getMaxDerivative = 3
  CASE (EL_E030, EL_E031, EL_EM30, EL_EM31)
    ! Function + 1st derivative
    elem_getMaxDerivative = 3
  CASE DEFAULT
    ! We don't know
    elem_getMaxDerivative = DER_MAXNDER
  END SELECT

  END FUNCTION

!**************************************************************************
! Generic element routines
! The following two routines define a generic element. Depending on
! the element type identifier ieltyp, the right element is called.
 
!<subroutine>  

  PURE SUBROUTINE elem_generic (ieltyp, Dcoords, Djac, ddetj, Bder, &
                                Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  ! ieltyp defines the element type that is used. Depending on
  ! the type of the element (parametric or nonparametric), dx
  ! and dy must be given either on the reference element or on the
  ! real element.
!</description>

!<input>
  ! Element type identifier.
  INTEGER, INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  REAL(DP), DIMENSION(EL_NJACENTRIES2D), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  REAL(DP), INTENT(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Coordinate of the point where to evaluate.
  ! The dimension depends on the coordinate system of the actual element.
  ! For triangular elements with barycentric 2D coordinates:
  !  Dcoord(1) = 1st barycentric coordinate
  !  Dcoord(2) = 2nd barycentric coordinate
  !  Dcoord(3) = 3rd barycentric coordinate
  ! For quadrilateral elements in real 2D coordinates
  !  Dcoord(1) = x-coordinate
  !  Dcoord(2) = y-coordinate
  ! For parametric elements, the point must be on the reference element, 
  ! for nonparametric elements on the real element.
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dpoint
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i'th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>

  ! Choose the right element subroutine to call.
  SELECT CASE (ieltyp)
  CASE (EL_P0)
    CALL elem_P0 (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_P1)
    CALL elem_P1 (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_P2)
    CALL elem_P2 (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_Q0)
    CALL elem_Q0 (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_Q1)
    CALL elem_Q1 (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_Q2)
    CALL elem_Q2 (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_EM30)
    CALL elem_EM30 (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE DEFAULT
    ! Element not implemened!
    ! Thow a floating point exception so that the program stops here!
    ! We cannot use "PRINT" here as the routine is PURE!
    CALL sys_throwFPE()
  END SELECT

  END SUBROUTINE 
  
  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_generic_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                     Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier
  INTEGER, INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  REAL(DP), DIMENSION(EL_NJACENTRIES2D,npoints), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points:
  REAL(DP), DIMENSION(npoints), INTENT(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! For parametric elements, the coordinates are expected on the 
  ! reference element. For nonparametric elements, the coordinates
  ! are expected on the real element!
  !   Dpoints(1,.)=x-coordinates,
  !   Dpoints(2,.)=y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,npoints), INTENT(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>

  ! local variables
  INTEGER :: i

  ! Choose the right element subroutine to call.
  SELECT CASE (ieltyp)
  CASE (EL_P0)
    CALL elem_P0_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_P1)
    CALL elem_P1_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_P2)
    CALL elem_P2_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_Q0)
    CALL elem_Q0_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_Q1)
    CALL elem_Q1_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_EM30)
    CALL elem_EM30_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE DEFAULT
    ! Compatibility handling: evaluate all points separately
    DO i=1,npoints
      CALL elem_generic (ieltyp, Dcoords, Djac(:,i), Ddetj(i), Bder, &
                         Dpoints(:,i), Dbas(:,:,i))
    END DO
  END SELECT

  END SUBROUTINE 
  
  !************************************************************************
  
!<subroutine>  

  SUBROUTINE elem_generic_sim (ieltyp, Dcoords, Djac, Ddetj, &
                               Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier
  INTEGER, INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE,nelements), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
  REAL(DP), DIMENSION(EL_NJACENTRIES2D,npoints,nelements), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  REAL(DP), DIMENSION(npoints,nelements), INTENT(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected 
  ! - on the reference element, if ieltyp identifies a parametric element
  ! - on the real element, if ieltyp identifies a nonparametric element
  ! It's assumed that:
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! array [1..EL_MAXNBAS,1..EL_MAXNDER,1..npoints,nelements] of double
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  REAL(DP), DIMENSION(:,:,:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>

  INTEGER :: i

  ! Choose the right element subroutine to call.
  SELECT CASE (ieltyp)
  CASE (EL_P0)
    CALL elem_P0_sim (ieltyp, Dcoords, Djac, Ddetj, &
                      Bder, Dbas, npoints, nelements, Dpoints)
  CASE (EL_P1)
    CALL elem_P1_sim (ieltyp, Dcoords, Djac, Ddetj, &
                      Bder, Dbas, npoints, nelements, Dpoints)
  CASE (EL_P2)
    CALL elem_P2_sim (ieltyp, Dcoords, Djac, Ddetj, &
                      Bder, Dbas, npoints, nelements, Dpoints)
  CASE (EL_Q0)
    CALL elem_Q0_sim (ieltyp, Dcoords, Djac, Ddetj, &
                      Bder, Dbas, npoints, nelements, Dpoints)
  CASE (EL_Q1)
    CALL elem_Q1_sim (ieltyp, Dcoords, Djac, Ddetj, &
                      Bder, Dbas, npoints, nelements, Dpoints)
  CASE (EL_Q2)
    CALL elem_Q2_sim (ieltyp, Dcoords, Djac, Ddetj, &
                      Bder, Dbas, npoints, nelements, Dpoints)
  CASE (EL_EM30)
    CALL elem_EM30_sim (ieltyp, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
  CASE DEFAULT
    ! Compatibility handling: evaluate on all elements separately
    DO i=1,nelements
      CALL elem_generic_mult (ieltyp, Dcoords(:,:,i), Djac(:,:,i), Ddetj(:,i), &
                              Bder, Dbas(:,:,:,i), npoints, Dpoints(:,:,i))
    END DO
  END SELECT

  END SUBROUTINE 
  
  !************************************************************************
  
!<function>  

  ELEMENTAL LOGICAL FUNCTION elem_isnonparametric (ieltyp) RESULT (inonpar)

  !<description>
  
  ! Determines whether an element ieltyp is parametric or nonparametric.
  
  !</description>
  
  !<result>
  ! =true, if the element is nonparametric. 
  ! =false, if the element is parametric.
  !</result>
  
  !<input>
  
  ! Element type qualifier.
  INTEGER, INTENT(IN) :: ieltyp
  
  !</input>
 
!</function>
 
  inonpar = ieltyp .LT. -1
  
  END FUNCTION

!**************************************************************************
! Element subroutines for parametric P0 element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  PURE SUBROUTINE elem_P0 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                           Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q0.
  INTEGER, INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(EL_NJACENTRIES2D), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! REMARK: Not used by this special type of element!
  REAL(DP), INTENT(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Barycentric coordinates of the point where to evaluate
  REAL(DP), DIMENSION(3), INTENT(IN) :: Dpoint
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i'th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>

  ! Clear the output array
  !Dbas = 0.0_DP
    
  ! Q0 is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC) = 1.0_DP

  END SUBROUTINE 
  
  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_P0_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_P0.
  INTEGER, INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(EL_NJACENTRIES2D,npoints), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(npoints), INTENT(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(3,npoints)
  !  Dpoints(1,.) = 1st barycentric coordinate
  !  Dpoints(2,.) = 2nd barycentric coordinate
  !  Dpoints(3,.) = 3rd barycentric coordinate
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>

  ! Clear the output array
  !Dbas = 0.0_DP
  
  ! Q0 is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC,:) = 1.0_DP

  END SUBROUTINE 

  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_P0_sim (ieltyp, Dcoords, Djac, Ddetj, &
                               Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>

  ! Element type identifier. Must be =EL_P0.
  INTEGER, INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(EL_NJACENTRIES2D,npoints), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(npoints), INTENT(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(3,npoints,nelements)
  !  Dpoints(1,.) = 1st barycentric coordinate
  !  Dpoints(2,.) = 2nd barycentric coordinate
  !  Dpoints(3,.) = 3rd barycentric coordinate
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,EL_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  REAL(DP), DIMENSION(:,:,:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>

  ! Clear the output array
  !Dbas = 0.0_DP
  
  ! Q0 is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC,:,:) = 1.0_DP

  END SUBROUTINE 

!**************************************************************************
! Element subroutines for parametric P1 element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  PURE SUBROUTINE elem_P1 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                           Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_P1.
  INTEGER, INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  REAL(DP), DIMENSION(EL_NJACENTRIES2D), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  REAL(DP), INTENT(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Barycentric coordinates of the point where to evaluate
  REAL(DP), DIMENSION(3), INTENT(IN) :: Dpoint
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i'th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>

  REAL(DP) :: dxj !auxiliary variable
  
  ! Clear the output array
  !Dbas = 0.0_DP
    
  ! Remark: The P1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
!  if (el_bder(DER_FUNC)) then
    Dbas(1:3,DER_FUNC) = Dpoint(1:3)
!  endif
  
  !if x-or y-derivatives are desired
!  if ((el_bder(DER_DERIV_X)) .or. (el_bder(DER_DERIV_Y))) then
    dxj = 1E0_DP / ddetj
    
    !x-derivatives on current element
!    if (el_bder(DER_DERIV_X)) then
      Dbas(1,DER_DERIV_X) = -(Djac(4)-Djac(2))*dxj
      Dbas(2,DER_DERIV_X) =  Djac(4)*dxj
      Dbas(3,DER_DERIV_X) = -Djac(2)*dxj
!    endif
    
    !y-derivatives on current element
!    if (el_bder(DER_DERIV_Y)) then
      Dbas(1,DER_DERIV_Y) = (Djac(3)-Djac(1))*dxj
      Dbas(2,DER_DERIV_Y) = -Djac(3)*dxj
      Dbas(3,DER_DERIV_Y) =  Djac(1)*dxj
!    endif
!  endif
    
  END SUBROUTINE 
  
  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_P1_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

  !<input>

  ! Element type identifier. Must be =EL_P1.
  INTEGER, INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  REAL(DP), DIMENSION(EL_NJACENTRIES2D,npoints), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  REAL(DP), DIMENSION(npoints), INTENT(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(3,npoints).
  !  Dpoints(1,.)=1st barycentric coordinate
  !  Dpoints(2,.)=2nd barycentric coordinate
  !  Dpoints(3,.)=3rd barycentric coordinate
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dpoints

  !</input>
  
  !<output>
  
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: Dbas
  
  !</output>

! </subroutine>

  REAL(DP),DIMENSION(npoints) :: dxj ! auxiliary variable
  
  INTEGER :: i   ! point counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The P1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  !IF (Bder(DER_FUNC)) THEN
    DO i=1,npoints
      Dbas(1,DER_FUNC,i) = Dpoints(1,i)
      Dbas(2,DER_FUNC,i) = Dpoints(2,i)
      Dbas(3,DER_FUNC,i) = Dpoints(3,i)
    END DO
  !ENDIF
  
  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN
    dxj = 1E0_DP / Ddetj
    
    !x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      DO i=1,npoints
        Dbas(1,DER_DERIV_X,i) = -(Djac(4,i)-Djac(2,i))*dxj(i)
        Dbas(2,DER_DERIV_X,i) =  Djac(4,i)*dxj(i)
        Dbas(3,DER_DERIV_X,i) = -Djac(2,i)*dxj(i)
!      END DO
!    ENDIF
    
    !y-derivatives on current element
!    IF (Bder(DER_DERIV_Y)) THEN
!      DO i=1,npoints
        Dbas(1,DER_DERIV_Y,i) = (Djac(3,i)-Djac(1,i))*dxj(i)
        Dbas(2,DER_DERIV_Y,i) = -Djac(3,i)*dxj(i)
        Dbas(3,DER_DERIV_Y,i) =  Djac(1,i)*dxj(i)
      END DO
!    ENDIF
!  ENDIF
    
  END SUBROUTINE 

  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_P1_sim (ieltyp, Dcoords, Djac, Ddetj, &
                               Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_P1.
  INTEGER, INTENT(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE,nelements), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
  REAL(DP), DIMENSION(EL_NJACENTRIES2D,npoints,nelements), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  REAL(DP), DIMENSION(npoints,nelements), INTENT(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(3,npoints,nelements)
  !  Dpoints(1,.) = 1st barycentric coordinate
  !  Dpoints(2,.) = 2nd barycentric coordinate
  !  Dpoints(3,.) = 3rd barycentric coordinate
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i'th 
  !   basis function of the finite element k in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,EL_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  REAL(DP), DIMENSION(:,:,:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>

  REAL(DP),DIMENSION(npoints) :: dxj !auxiliary variable
  
  INTEGER :: i   ! point counter
  INTEGER :: j   ! element counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  !if function values are desired
  IF (Bder(DER_FUNC)) THEN
  
    DO j=1,nelements
    
      DO i=1,npoints
        Dbas(1,DER_FUNC,i,j) = Dpoints(1,i,j)
        Dbas(2,DER_FUNC,i,j) = Dpoints(2,i,j)
        Dbas(3,DER_FUNC,i,j) = Dpoints(3,i,j)
      END DO
      
    END DO
    
  END IF
    
  !if x-or y-derivatives are desired
  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN
  
    DO j=1,nelements
      dxj = 1E0_DP / Ddetj(:,j)
      
      !x-derivatives on current element
!      IF (Bder(DER_DERIV_X)) THEN
        DO i=1,npoints
          Dbas(1,DER_DERIV_X,i,j) = -(Djac(4,i,j)-Djac(2,i,j))*dxj(i)
          Dbas(2,DER_DERIV_X,i,j) =  Djac(4,i,j)*dxj(i)
          Dbas(3,DER_DERIV_X,i,j) = -Djac(2,i,j)*dxj(i)
        END DO
!      ENDIF
      
      !y-derivatives on current element
!      IF (Bder(DER_DERIV_Y)) THEN
        DO i=1,npoints
          Dbas(1,DER_DERIV_Y,i,j) = (Djac(3,i,j)-Djac(1,i,j))*dxj(i)
          Dbas(2,DER_DERIV_Y,i,j) = -Djac(3,i,j)*dxj(i)
          Dbas(3,DER_DERIV_Y,i,j) =  Djac(1,i,j)*dxj(i)
        END DO
!      ENDIF

    END DO
      
  END IF
    
  END SUBROUTINE 

!**************************************************************************
! Element subroutines for parametric P2 element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  PURE SUBROUTINE elem_P2 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                           Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_P1.
  INTEGER, INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  REAL(DP), DIMENSION(EL_NJACENTRIES2D), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  REAL(DP), INTENT(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Barycentric coordinates of the point where to evaluate
  REAL(DP), DIMENSION(3), INTENT(IN) :: Dpoint
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i'th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>

  REAL(DP) :: dxj,dp1,dp2,dp3 !auxiliary variables
  
  ! Clear the output array
  !Dbas = 0.0_DP
    
  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  if (Bder(DER_FUNC)) then
    Dbas(1,DER_FUNC)= Dpoint(1)*(Dpoint(1)-Dpoint(2)-Dpoint(3))
    Dbas(2,DER_FUNC)= Dpoint(2)*(Dpoint(2)-Dpoint(1)-Dpoint(3))
    Dbas(3,DER_FUNC)= Dpoint(3)*(Dpoint(3)-Dpoint(1)-Dpoint(2))
    Dbas(4,DER_FUNC)= 4.0_DP*Dpoint(1)*Dpoint(2)
    Dbas(5,DER_FUNC)= 4.0_DP*Dpoint(2)*Dpoint(3)
    Dbas(6,DER_FUNC)= 4.0_DP*Dpoint(1)*Dpoint(3)
  endif
  
  !if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then
    dxj = 1E0_DP / ddetj
    dp1=(1.0_DP-4.0_DP*Dpoint(1))*dxj
    dp2=(4.0_DP*Dpoint(2)-1.0_DP)*dxj
    dp3=(4.0_DP*Dpoint(3)-1.0_DP)*dxj
    
    !x-derivatives on current element
    if (Bder(DER_DERIV_X)) then
      Dbas(1,DER_DERIV_X)= (Djac(4)-Djac(2))*dp1
      Dbas(2,DER_DERIV_X)= Djac(4)*dp2
      Dbas(3,DER_DERIV_X)=-Djac(2)*dp3
      Dbas(4,DER_DERIV_X)= &
          4.0_DP*(Dpoint(1)*Djac(4)-Dpoint(2)*(Djac(4)-Djac(2)))*dxj
      Dbas(5,DER_DERIV_X)= &
          4.0_DP*(Dpoint(3)*Djac(4)-Dpoint(2)*Djac(2))*dxj
      Dbas(6,DER_DERIV_X)= &
          4.0_DP*(-Dpoint(1)*Djac(2)-Dpoint(3)*(Djac(4)-Djac(2)))*dxj
    endif
    
    !y-derivatives on current element
    if (Bder(DER_DERIV_Y)) then
      Dbas(1,DER_DERIV_Y)=-(Djac(3)-Djac(1))*dp1
      Dbas(2,DER_DERIV_Y)=-Djac(3)*dp2
      Dbas(3,DER_DERIV_Y)= Djac(1)*dp3
      Dbas(4,DER_DERIV_Y)= &
          4.0_DP*(-Dpoint(1)*Djac(3)+Dpoint(2)*(Djac(3)-Djac(1)))*dxj
      Dbas(5,DER_DERIV_Y)= &
          4.0_DP*(-Dpoint(3)*Djac(3)+Dpoint(2)*Djac(1))*dxj
      Dbas(6,DER_DERIV_Y)= &
          4.0_DP*(Dpoint(1)*Djac(1)+Dpoint(3)*(Djac(3)-Djac(1)))*dxj
    endif
  endif
    
  END SUBROUTINE 

  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_P2_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

  !<input>

  ! Element type identifier. Must be =EL_P2.
  INTEGER, INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  REAL(DP), DIMENSION(EL_NJACENTRIES2D,npoints), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  REAL(DP), DIMENSION(npoints), INTENT(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(3,npoints).
  !  Dpoints(1,.)=1st barycentric coordinate
  !  Dpoints(2,.)=2nd barycentric coordinate
  !  Dpoints(3,.)=3rd barycentric coordinate
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dpoints

  !</input>
  
  !<output>
  
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: Dbas
  
  !</output>

! </subroutine>

  REAL(DP),DIMENSION(npoints) :: Dxj,Dp1,Dp2,Dp3 ! auxiliary variable
  
  INTEGER :: i   ! point counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The P1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  IF (Bder(DER_FUNC)) THEN
    DO i=1,npoints
      Dbas(1,DER_FUNC,i)= Dpoints(1,i)*(Dpoints(1,i)-Dpoints(2,i)-Dpoints(3,i))
      Dbas(2,DER_FUNC,i)= Dpoints(2,i)*(Dpoints(2,i)-Dpoints(1,i)-Dpoints(3,i))
      Dbas(3,DER_FUNC,i)= Dpoints(3,i)*(Dpoints(3,i)-Dpoints(1,i)-Dpoints(2,i))
      Dbas(4,DER_FUNC,i)= 4.0_DP*Dpoints(1,i)*Dpoints(2,i)
      Dbas(5,DER_FUNC,i)= 4.0_DP*Dpoints(2,i)*Dpoints(3,i)
      Dbas(6,DER_FUNC,i)= 4.0_DP*Dpoints(1,i)*Dpoints(3,i)
    END DO
  ENDIF
  
  !if x-or y-derivatives are desired
  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN
    Dxj = 1E0_DP / Ddetj
    DO i=1,npoints
      Dp1(i)=(1.0_DP-4.0_DP*Dpoints(1,i))*Dxj(i)
      Dp2(i)=(4.0_DP*Dpoints(2,i)-1.0_DP)*Dxj(i)
      Dp3(i)=(4.0_DP*Dpoints(3,i)-1.0_DP)*Dxj(i)
    END DO
    
    !x-derivatives on current element
    IF (Bder(DER_DERIV_X)) THEN
      DO i=1,npoints
        Dbas(1,DER_DERIV_X,i)= (Djac(4,i)-Djac(2,i))*Dp1(i)
        Dbas(2,DER_DERIV_X,i)= Djac(4,i)*Dp2(i)
        Dbas(3,DER_DERIV_X,i)=-Djac(2,i)*Dp3(i)
        Dbas(4,DER_DERIV_X,i)= &
            4.0_DP*(Dpoints(1,i)*Djac(4,i) &
            -Dpoints(2,i)*(Djac(4,i)-Djac(2,i)))*Dxj(i)
        Dbas(5,DER_DERIV_X,i)= &
            4.0_DP*(Dpoints(3,i)*Djac(4,i) &
            -Dpoints(2,i)*Djac(2,i))*Dxj(i)
        Dbas(6,DER_DERIV_X,i)= &
            4.0_DP*(-Dpoints(1,i)*Djac(2,i)&
            -Dpoints(3,i)*(Djac(4,i)-Djac(2,i)))*Dxj(i)
      END DO
    ENDIF
    
    !y-derivatives on current element
    IF (Bder(DER_DERIV_Y)) THEN
      DO i=1,npoints
        Dbas(1,DER_DERIV_Y,i)=-(Djac(3,i)-Djac(1,i))*Dp1(i)
        Dbas(2,DER_DERIV_Y,i)=-Djac(3,i)*Dp2(i)
        Dbas(3,DER_DERIV_Y,i)= Djac(1,i)*Dp3(i)
        Dbas(4,DER_DERIV_Y,i)= &
            4.0_DP*(-Dpoints(1,i)*Djac(3,i) &
                    +Dpoints(2,i)*(Djac(3,i)-Djac(1,i)))*Dxj(i)
        Dbas(5,DER_DERIV_Y,i)= &
            4.0_DP*(-Dpoints(3,i)*Djac(3,i) &
                    +Dpoints(2,i)*Djac(1,i))*Dxj(i)
        Dbas(6,DER_DERIV_Y,i)= &
            4.0_DP*(Dpoints(1,i)*Djac(1,i) &
                    +Dpoints(3,i)*(Djac(3,i)-Djac(1,i)))*Dxj(i)
      END DO
    ENDIF
  ENDIF
    
  END SUBROUTINE 

  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_P2_sim (ieltyp, Dcoords, Djac, Ddetj, &
                               Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_P2.
  INTEGER, INTENT(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE,nelements), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
  REAL(DP), DIMENSION(EL_NJACENTRIES2D,npoints,nelements), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  REAL(DP), DIMENSION(npoints,nelements), INTENT(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(3,npoints,nelements)
  !  Dpoints(1,.) = 1st barycentric coordinate
  !  Dpoints(2,.) = 2nd barycentric coordinate
  !  Dpoints(3,.) = 3rd barycentric coordinate
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i'th 
  !   basis function of the finite element k in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,EL_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  REAL(DP), DIMENSION(:,:,:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>

  REAL(DP),DIMENSION(npoints) :: dxj,Dp1,Dp2,Dp3 !auxiliary variable
  
  INTEGER :: i   ! point counter
  INTEGER :: j   ! element counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  !if function values are desired
  IF (Bder(DER_FUNC)) THEN
  
    DO j=1,nelements
    
      DO i=1,npoints
        Dbas(1,DER_FUNC,i,j)= &
            Dpoints(1,i,j)*(Dpoints(1,i,j)-Dpoints(2,i,j)-Dpoints(3,i,j))
        Dbas(2,DER_FUNC,i,j)= &
            Dpoints(2,i,j)*(Dpoints(2,i,j)-Dpoints(1,i,j)-Dpoints(3,i,j))
        Dbas(3,DER_FUNC,i,j)= &
            Dpoints(3,i,j)*(Dpoints(3,i,j)-Dpoints(1,i,j)-Dpoints(2,i,j))
        Dbas(4,DER_FUNC,i,j)= &
            4.0_DP*Dpoints(1,i,j)*Dpoints(2,i,j)
        Dbas(5,DER_FUNC,i,j)= &
            4.0_DP*Dpoints(2,i,j)*Dpoints(3,i,j)
        Dbas(6,DER_FUNC,i,j)= &
            4.0_DP*Dpoints(1,i,j)*Dpoints(3,i,j)
      END DO
      
    END DO
    
  END IF
    
  !if x-or y-derivatives are desired
  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN
  
    DO j=1,nelements
      dxj = 1E0_DP / Ddetj(:,j)
      
      DO i=1,npoints
        Dp1(i)=(1.0_DP-4.0_DP*Dpoints(1,i,j))*Dxj(i)
        Dp2(i)=(4.0_DP*Dpoints(2,i,j)-1.0_DP)*Dxj(i)
        Dp3(i)=(4.0_DP*Dpoints(3,i,j)-1.0_DP)*Dxj(i)
      END DO
      
      
      !x-derivatives on current element
!      IF (Bder(DER_DERIV_X)) THEN
        DO i=1,npoints
          Dbas(1,DER_DERIV_X,i,j)= (Djac(4,i,j)-Djac(2,i,j))*Dp1(i)
          Dbas(2,DER_DERIV_X,i,j)= Djac(4,i,j)*Dp2(i)
          Dbas(3,DER_DERIV_X,i,j)=-Djac(2,i,j)*Dp3(i)
          Dbas(4,DER_DERIV_X,i,j)= &
              4.0_DP*(Dpoints(1,i,j)*Djac(4,i,j) &
              -Dpoints(2,i,j)*(Djac(4,i,j)-Djac(2,i,j)))*Dxj(i)
          Dbas(5,DER_DERIV_X,i,j)= &
              4.0_DP*(Dpoints(3,i,j)*Djac(4,i,j) &
              -Dpoints(2,i,j)*Djac(2,i,j))*Dxj(i)
          Dbas(6,DER_DERIV_X,i,j)= &
              4.0_DP*(-Dpoints(1,i,j)*Djac(2,i,j)&
              -Dpoints(3,i,j)*(Djac(4,i,j)-Djac(2,i,j)))*Dxj(i)
        END DO
!      ENDIF
      
      !y-derivatives on current element
!      IF (Bder(DER_DERIV_Y)) THEN
        DO i=1,npoints
          Dbas(1,DER_DERIV_Y,i,j)=-(Djac(3,i,j)-Djac(1,i,j))*Dp1(i)
          Dbas(2,DER_DERIV_Y,i,j)=-Djac(3,i,j)*Dp2(i)
          Dbas(3,DER_DERIV_Y,i,j)= Djac(1,i,j)*Dp3(i)
          Dbas(4,DER_DERIV_Y,i,j)= &
              4.0_DP*(-Dpoints(1,i,j)*Djac(3,i,j) &
                      +Dpoints(2,i,j)*(Djac(3,i,j)-Djac(1,i,j)))*Dxj(i)
          Dbas(5,DER_DERIV_Y,i,j)= &
              4.0_DP*(-Dpoints(3,i,j)*Djac(3,i,j) &
                      +Dpoints(2,i,j)*Djac(1,i,j))*Dxj(i)
          Dbas(6,DER_DERIV_Y,i,j)= &
              4.0_DP*(Dpoints(1,i,j)*Djac(1,i,j) &
                      +Dpoints(3,i,j)*(Djac(3,i,j)-Djac(1,i,j)))*Dxj(i)
        END DO
!      ENDIF

    END DO
      
  END IF
    
  END SUBROUTINE 

!**************************************************************************
! Element subroutines for parametric Q0 element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  PURE SUBROUTINE elem_Q0 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                           Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q0.
  INTEGER, INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(EL_NJACENTRIES2D), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! REMARK: Not used by this special type of element!
  REAL(DP), INTENT(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  REAL(DP), DIMENSION(2), INTENT(IN) :: Dpoint
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i'th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>

  ! Clear the output array
  !Dbas = 0.0_DP
    
  ! Q0 is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC) = 1.0_DP

  END SUBROUTINE 
  
  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_Q0_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, Dpoints)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q0.
  INTEGER, INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(EL_NJACENTRIES2D,npoints), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(npoints), INTENT(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! Dpoints(1,.)=x-coordinates,
  ! Dpoints(2,.)=y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,npoints), INTENT(IN) :: Dpoints
  
  !</input>
  
  !<output>
  
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: Dbas
  
  !</output>

! </subroutine>

  ! Clear the output array
  !Dbas = 0.0_DP
  
  ! Q0 is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC,:) = 1.0_DP

  END SUBROUTINE 

  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_Q0_sim (ieltyp, Dcoords, Djac, Ddetj, &
                               Bder, Dbas, npoints, nelements, Dpoints)

  !<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q0.
  INTEGER, INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(EL_NJACENTRIES2D,npoints), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(npoints), INTENT(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  REAL(DP), DIMENSION(NDIM2D,npoints), INTENT(IN) :: Dpoints
  
  !</input>
  
  !<output>
  
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,EL_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  REAL(DP), DIMENSION(:,:,:,:), INTENT(OUT) :: Dbas
  
  !</output>

! </subroutine>

  ! Clear the output array
  !Dbas = 0.0_DP
  
  ! Q0 is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC,:,:) = 1.0_DP

  END SUBROUTINE 

!**************************************************************************
! Element subroutines for parametric Q1 element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  PURE SUBROUTINE elem_Q1 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                           Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q1.
  INTEGER, INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  REAL(DP), DIMENSION(EL_NJACENTRIES2D), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  REAL(DP), INTENT(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  REAL(DP), DIMENSION(2), INTENT(IN) :: Dpoint
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i'th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>

  !auxiliary vector containing the first derivatives on the reference element
  REAL(DP), DIMENSION(EL_MAXNVE,NDIM2D) :: Dhelp
  REAL(DP) :: dx,dy

  REAL(DP) :: dxj !auxiliary variable
  
  ! Clear the output array
  !Dbas = 0.0_DP
  
  dx = Dpoint(1)
  dy = Dpoint(2)
    
  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
!  if (el_bder(DER_FUNC)) then
    Dbas(1,DER_FUNC) = 0.25E0_DP*(1E0_DP-dx)*(1E0_DP-dy)
    Dbas(2,DER_FUNC) = 0.25E0_DP*(1E0_DP+dx)*(1E0_DP-dy)
    Dbas(3,DER_FUNC) = 0.25E0_DP*(1E0_DP+dx)*(1E0_DP+dy)
    Dbas(4,DER_FUNC) = 0.25E0_DP*(1E0_DP-dx)*(1E0_DP+dy)
!  endif
  
  !if x-or y-derivatives are desired
!  if ((el_bder(DER_DERIV_X)) .or. (el_bder(DER_DERIV_Y))) then
    dxj = 0.25E0_DP / ddetj
    
    !x- and y-derivatives on reference element
    Dhelp(1,1) =-(1E0_DP-dy)
    Dhelp(2,1) = (1E0_DP-dy)
    Dhelp(3,1) = (1E0_DP+dy)
    Dhelp(4,1) =-(1E0_DP+dy)
    Dhelp(1,2) =-(1E0_DP-dx)
    Dhelp(2,2) =-(1E0_DP+dx)
    Dhelp(3,2) = (1E0_DP+dx)
    Dhelp(4,2) = (1E0_DP-dx)
      
    !x-derivatives on current element
!    if (el_bder(DER_DERIV_X)) then
      Dbas(1,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(1,1) - Djac(2) * Dhelp(1,2))
      Dbas(2,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(2,1) - Djac(2) * Dhelp(2,2))
      Dbas(3,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(3,1) - Djac(2) * Dhelp(3,2))
      Dbas(4,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(4,1) - Djac(2) * Dhelp(4,2))
!    endif
    
    !y-derivatives on current element
!    if (el_bder(DER_DERIV_Y)) then
      Dbas(1,DER_DERIV_Y) = -dxj * (Djac(3) * Dhelp(1,1) - Djac(1) * Dhelp(1,2))
      Dbas(2,DER_DERIV_Y) = -dxj * (Djac(3) * Dhelp(2,1) - Djac(1) * Dhelp(2,2))
      Dbas(3,DER_DERIV_Y) = -dxj * (Djac(3) * Dhelp(3,1) - Djac(1) * Dhelp(3,2))
      Dbas(4,DER_DERIV_Y) = -dxj * (Djac(3) * Dhelp(4,1) - Djac(1) * Dhelp(4,2))
!    endif
!  endif
    
  END SUBROUTINE 
  
  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_Q1_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1.
  INTEGER, INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  REAL(DP), DIMENSION(EL_NJACENTRIES2D,npoints), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  REAL(DP), DIMENSION(npoints), INTENT(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(NDIM2D,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  REAL(DP), DIMENSION(EL_MAXNVE,NDIM2D,npoints) :: Dhelp

  REAL(DP),DIMENSION(npoints) :: dxj !auxiliary variable
  
  INTEGER :: i   ! point counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  !IF (Bder(DER_FUNC)) THEN
    DO i=1,npoints
      Dbas(1,DER_FUNC,i) = 0.25E0_DP*(1E0_DP-Dpoints(1,i))*(1E0_DP-Dpoints(2,i))
      Dbas(2,DER_FUNC,i) = 0.25E0_DP*(1E0_DP+Dpoints(1,i))*(1E0_DP-Dpoints(2,i))
      Dbas(3,DER_FUNC,i) = 0.25E0_DP*(1E0_DP+Dpoints(1,i))*(1E0_DP+Dpoints(2,i))
      Dbas(4,DER_FUNC,i) = 0.25E0_DP*(1E0_DP-Dpoints(1,i))*(1E0_DP+Dpoints(2,i))
    END DO
  !ENDIF
  
  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN
    dxj = 0.25E0_DP / Ddetj
    
    !x- and y-derivatives on reference element
    DO i=1,npoints
      Dhelp(1,1,i) =-(1E0_DP-Dpoints(2,i))
      Dhelp(2,1,i) = (1E0_DP-Dpoints(2,i))
      Dhelp(3,1,i) = (1E0_DP+Dpoints(2,i))
      Dhelp(4,1,i) =-(1E0_DP+Dpoints(2,i))
      Dhelp(1,2,i) =-(1E0_DP-Dpoints(1,i))
      Dhelp(2,2,i) =-(1E0_DP+Dpoints(1,i))
      Dhelp(3,2,i) = (1E0_DP+Dpoints(1,i))
      Dhelp(4,2,i) = (1E0_DP-Dpoints(1,i))
    END DO
      
    !x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      DO i=1,npoints
        Dbas(1,DER_DERIV_X,i) = dxj(i) * (Djac(4,i) * Dhelp(1,1,i) - Djac(2,i) * Dhelp(1,2,i))
        Dbas(2,DER_DERIV_X,i) = dxj(i) * (Djac(4,i) * Dhelp(2,1,i) - Djac(2,i) * Dhelp(2,2,i))
        Dbas(3,DER_DERIV_X,i) = dxj(i) * (Djac(4,i) * Dhelp(3,1,i) - Djac(2,i) * Dhelp(3,2,i))
        Dbas(4,DER_DERIV_X,i) = dxj(i) * (Djac(4,i) * Dhelp(4,1,i) - Djac(2,i) * Dhelp(4,2,i))
!      END DO
!    ENDIF
    
    !y-derivatives on current element
!    IF (Bder(DER_DERIV_Y)) THEN
!      DO i=1,npoints
        Dbas(1,DER_DERIV_Y,i) = -dxj(i) * (Djac(3,i) * Dhelp(1,1,i) - Djac(1,i) * Dhelp(1,2,i))
        Dbas(2,DER_DERIV_Y,i) = -dxj(i) * (Djac(3,i) * Dhelp(2,1,i) - Djac(1,i) * Dhelp(2,2,i))
        Dbas(3,DER_DERIV_Y,i) = -dxj(i) * (Djac(3,i) * Dhelp(3,1,i) - Djac(1,i) * Dhelp(3,2,i))
        Dbas(4,DER_DERIV_Y,i) = -dxj(i) * (Djac(3,i) * Dhelp(4,1,i) - Djac(1,i) * Dhelp(4,2,i))
      END DO
!    ENDIF
!  ENDIF
    
  END SUBROUTINE 

  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_Q1_sim (ieltyp, Dcoords, Djac, Ddetj, &
                               Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1.
  INTEGER, INTENT(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE,nelements), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
  REAL(DP), DIMENSION(EL_NJACENTRIES2D,npoints,nelements), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  REAL(DP), DIMENSION(npoints,nelements), INTENT(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(NDIM2D,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i'th 
  !   basis function of the finite element k in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,EL_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  REAL(DP), DIMENSION(:,:,:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  REAL(DP), DIMENSION(EL_MAXNVE,NDIM2D,npoints) :: Dhelp

  REAL(DP),DIMENSION(npoints) :: dxj !auxiliary variable
  
  INTEGER :: i   ! point counter
  INTEGER :: j   ! element counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  IF (Bder(DER_FUNC)) THEN
  
    DO j=1,nelements
    
      DO i=1,npoints
        Dbas(1,DER_FUNC,i,j) = 0.25_DP*(1.0_DP-Dpoints(1,i,j))*(1.0_DP-Dpoints(2,i,j))
        Dbas(2,DER_FUNC,i,j) = 0.25_DP*(1.0_DP+Dpoints(1,i,j))*(1.0_DP-Dpoints(2,i,j))
        Dbas(3,DER_FUNC,i,j) = 0.25_DP*(1.0_DP+Dpoints(1,i,j))*(1.0_DP+Dpoints(2,i,j))
        Dbas(4,DER_FUNC,i,j) = 0.25_DP*(1.0_DP-Dpoints(1,i,j))*(1.0_DP+Dpoints(2,i,j))
      END DO
      
    END DO
    
  END IF
    
  !if x-or y-derivatives are desired
  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN
  
    DO j=1,nelements
      dxj = 0.25E0_DP / Ddetj(:,j)
      
      !x- and y-derivatives on reference element
      DO i=1,npoints
        Dhelp(1,1,i) =-(1E0_DP-Dpoints(2,i,j))
        Dhelp(2,1,i) = (1E0_DP-Dpoints(2,i,j))
        Dhelp(3,1,i) = (1E0_DP+Dpoints(2,i,j))
        Dhelp(4,1,i) =-(1E0_DP+Dpoints(2,i,j))
        Dhelp(1,2,i) =-(1E0_DP-Dpoints(1,i,j))
        Dhelp(2,2,i) =-(1E0_DP+Dpoints(1,i,j))
        Dhelp(3,2,i) = (1E0_DP+Dpoints(1,i,j))
        Dhelp(4,2,i) = (1E0_DP-Dpoints(1,i,j))
      END DO
        
      !x-derivatives on current element
!      IF (Bder(DER_DERIV_X)) THEN
        DO i=1,npoints
          Dbas(1,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(1,1,i) &
                                    - Djac(2,i,j) * Dhelp(1,2,i))
          Dbas(2,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(2,1,i) &
                                    - Djac(2,i,j) * Dhelp(2,2,i))
          Dbas(3,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(3,1,i) &
                                    - Djac(2,i,j) * Dhelp(3,2,i))
          Dbas(4,DER_DERIV_X,i,j) = dxj(i) * (Djac(4,i,j) * Dhelp(4,1,i) &
                                    - Djac(2,i,j) * Dhelp(4,2,i))
        END DO
!      ENDIF
      
      !y-derivatives on current element
!      IF (Bder(DER_DERIV_Y)) THEN
        DO i=1,npoints
          Dbas(1,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(1,1,i) &
                                    - Djac(1,i,j) * Dhelp(1,2,i))
          Dbas(2,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(2,1,i) &
                                    - Djac(1,i,j) * Dhelp(2,2,i))
          Dbas(3,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(3,1,i) &
                                    - Djac(1,i,j) * Dhelp(3,2,i))
          Dbas(4,DER_DERIV_Y,i,j) = -dxj(i) * (Djac(3,i,j) * Dhelp(4,1,i) &
                                    - Djac(1,i,j) * Dhelp(4,2,i))
        END DO
!      ENDIF

    END DO
      
  END IF
    
  END SUBROUTINE 

!**************************************************************************
! Element subroutines for parametric Q2 element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  PURE SUBROUTINE elem_Q2 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                           Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q2.
  INTEGER, INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  REAL(DP), DIMENSION(EL_NJACENTRIES2D), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  REAL(DP), INTENT(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Cartesian coordinates of the evaluation point on reference element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  REAL(DP), DIMENSION(2), INTENT(IN) :: Dpoint
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i'th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>

  !auxiliary vector containing the first derivatives on the reference element
  REAL(DP), DIMENSION(9,NDIM2D) :: Dhelp
  REAL(DP) :: dx,dy

  INTEGER :: idof
  REAL(DP) :: dxj !auxiliary variable
  
  REAL(DP), PARAMETER :: Q2 = 0.5_DP
  REAL(DP), PARAMETER :: Q4 = 0.25_DP
  
  ! Clear the output array
  !Dbas = 0.0_DP
  
  dx = Dpoint(1)
  dy = Dpoint(2)
    
  !if function values are desired
  if (Bder(DER_FUNC)) then
    Dbas(1,DER_FUNC)= Q4*(1.0_DP-dx)*(1.0_DP-dy)*dx*dy
    Dbas(2,DER_FUNC)=-Q4*(1.0_DP+dx)*(1.0_DP-dy)*dx*dy
    Dbas(3,DER_FUNC)= Q4*(1.0_DP+dx)*(1.0_DP+dy)*dx*dy
    Dbas(4,DER_FUNC)=-Q4*(1.0_DP-dx)*(1.0_DP+dy)*dx*dy
    Dbas(5,DER_FUNC)=-Q2*(1.0_DP-dx*dx)*(1.0_DP-dy)*dy
    Dbas(6,DER_FUNC)= Q2*(1.0_DP+dx)*(1.0_DP-dy*dy)*dx
    Dbas(7,DER_FUNC)= Q2*(1.0_DP-dx*dx)*(1.0_DP+dy)*dy
    Dbas(8,DER_FUNC)=-Q2*(1.0_DP-dx)*(1.0_DP-dy*dy)*dx
    Dbas(9,DER_FUNC)= (1.0_DP-dx*dx)*(1.0_DP-dy*dy)
  endif
  
  ! if x-or y-derivatives are desired
  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then
    dxj = 1.0E0_DP / ddetj
    
    !x- and y-derivatives on reference element
    Dhelp(1,1)= Q4*(1.0_DP-2.0_DP*dx)*(1.0_DP-dy)*dy
    Dhelp(2,1)=-Q4*(1.0_DP+2.0_DP*dx)*(1.0_DP-dy)*dy
    Dhelp(3,1)= Q4*(1.0_DP+2.0_DP*dx)*(1.0_DP+dy)*dy
    Dhelp(4,1)=-Q4*(1.0_DP-2.0_DP*dx)*(1.0_DP+dy)*dy
    Dhelp(5,1)= (1.0_DP-dy)*dx*dy
    Dhelp(6,1)= Q2*(1.0_DP+2.0_DP*dx)*(1.0_DP-dy*dy)
    Dhelp(7,1)=-(1.0_DP+dy)*dx*dy
    Dhelp(8,1)=-Q2*(1.0_DP-2.0_DP*dx)*(1.0_DP-dy*dy)
    Dhelp(9,1)=-2.0_DP*(1.0_DP-dy*dy)*dx

    Dhelp(1,2)= Q4*(1.0_DP-dx)*(1.0_DP-2.0_DP*dy)*dx
    Dhelp(2,2)=-Q4*(1.0_DP+dx)*(1.0_DP-2.0_DP*dy)*dx
    Dhelp(3,2)= Q4*(1.0_DP+dx)*(1.0_DP+2.0_DP*dy)*dx
    Dhelp(4,2)=-Q4*(1.0_DP-dx)*(1.0_DP+2.0_DP*dy)*dx
    Dhelp(5,2)=-Q2*(1.0_DP-dx*dx)*(1.0_DP-2.0_DP*dy)
    Dhelp(6,2)=-(1.0_DP+dx)*dx*dy
    Dhelp(7,2)= Q2*(1.0_DP-dx*dx)*(1.0_DP+2.0_DP*dy)
    Dhelp(8,2)= (1.0_DP-dx)*dx*dy
    Dhelp(9,2)=-2.0_DP*(1.0_DP-dx*dx)*dy

    ! x-derivatives on current element
    if (Bder(DER_DERIV_X)) then
      DO idof = 1,9
        Dbas(idof,DER_DERIV_X) = &
            dxj * (Djac(4) * Dhelp(idof,1) - Djac(2) * Dhelp(idof,2))
      END DO
    endif
    
    ! y-derivatives on current element
    if (Bder(DER_DERIV_Y)) then
      DO idof = 1,9
        Dbas(idof,DER_DERIV_Y) = &
            -dxj * (Djac(3) * Dhelp(idof,1) - Djac(1) * Dhelp(idof,2))
      END DO
    endif
  endif
    
  END SUBROUTINE 
  
  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_Q2_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q2.
  INTEGER, INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  REAL(DP), DIMENSION(EL_NJACENTRIES2D,npoints), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  REAL(DP), DIMENSION(npoints), INTENT(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(NDIM2D,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  REAL(DP), DIMENSION(9,NDIM2D,npoints) :: Dhelp
  REAL(DP) :: dx,dy
  REAL(DP),DIMENSION(npoints) :: Dxj !auxiliary variable
  
  INTEGER :: i,idof   ! point counter
    
  REAL(DP), PARAMETER :: Q2 = 0.5_DP
  REAL(DP), PARAMETER :: Q4 = 0.25_DP
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  !IF (Bder(DER_FUNC)) THEN
    DO i=1,npoints
      dx = Dpoints(1,i)
      dy = Dpoints(2,i)
      Dbas(1,DER_FUNC,i)= Q4*(1.0_DP-dx)*(1.0_DP-dy)*dx*dy
      Dbas(2,DER_FUNC,i)=-Q4*(1.0_DP+dx)*(1.0_DP-dy)*dx*dy
      Dbas(3,DER_FUNC,i)= Q4*(1.0_DP+dx)*(1.0_DP+dy)*dx*dy
      Dbas(4,DER_FUNC,i)=-Q4*(1.0_DP-dx)*(1.0_DP+dy)*dx*dy
      Dbas(5,DER_FUNC,i)=-Q2*(1.0_DP-dx*dx)*(1.0_DP-dy)*dy
      Dbas(6,DER_FUNC,i)= Q2*(1.0_DP+dx)*(1.0_DP-dy*dy)*dx
      Dbas(7,DER_FUNC,i)= Q2*(1.0_DP-dx*dx)*(1.0_DP+dy)*dy
      Dbas(8,DER_FUNC,i)=-Q2*(1.0_DP-dx)*(1.0_DP-dy*dy)*dx
      Dbas(9,DER_FUNC,i)= (1.0_DP-dx*dx)*(1.0_DP-dy*dy)
    END DO
  !ENDIF
  
  !if x-or y-derivatives are desired
  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN
    Dxj = 1.0E0_DP / Ddetj
    
    !x- and y-derivatives on reference element
    DO i=1,npoints
      dx = Dpoints(1,i)
      dy = Dpoints(2,i)

      !x- and y-derivatives on reference element
      Dhelp(1,1,i)= Q4*(1.0_DP-2.0_DP*dx)*(1.0_DP-dy)*dy
      Dhelp(2,1,i)=-Q4*(1.0_DP+2.0_DP*dx)*(1.0_DP-dy)*dy
      Dhelp(3,1,i)= Q4*(1.0_DP+2.0_DP*dx)*(1.0_DP+dy)*dy
      Dhelp(4,1,i)=-Q4*(1.0_DP-2.0_DP*dx)*(1.0_DP+dy)*dy
      Dhelp(5,1,i)= (1.0_DP-dy)*dx*dy
      Dhelp(6,1,i)= Q2*(1.0_DP+2.0_DP*dx)*(1.0_DP-dy*dy)
      Dhelp(7,1,i)=-(1.0_DP+dy)*dx*dy
      Dhelp(8,1,i)=-Q2*(1.0_DP-2.0_DP*dx)*(1.0_DP-dy*dy)
      Dhelp(9,1,i)=-2.0_DP*(1.0_DP-dy*dy)*dx

      Dhelp(1,2,i)= Q4*(1.0_DP-dx)*(1.0_DP-2.0_DP*dy)*dx
      Dhelp(2,2,i)=-Q4*(1.0_DP+dx)*(1.0_DP-2.0_DP*dy)*dx
      Dhelp(3,2,i)= Q4*(1.0_DP+dx)*(1.0_DP+2.0_DP*dy)*dx
      Dhelp(4,2,i)=-Q4*(1.0_DP-dx)*(1.0_DP+2.0_DP*dy)*dx
      Dhelp(5,2,i)=-Q2*(1.0_DP-dx*dx)*(1.0_DP-2.0_DP*dy)
      Dhelp(6,2,i)=-(1.0_DP+dx)*dx*dy
      Dhelp(7,2,i)= Q2*(1.0_DP-dx*dx)*(1.0_DP+2.0_DP*dy)
      Dhelp(8,2,i)= (1.0_DP-dx)*dx*dy
      Dhelp(9,2,i)=-2.0_DP*(1.0_DP-dx*dx)*dy
    END DO
      
    ! x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      DO i=1,npoints
        DO idof = 1,9
          Dbas(idof,DER_DERIV_X,i) = &
              Dxj(i) * (Djac(4,i) * Dhelp(idof,1,i) - Djac(2,i) * Dhelp(idof,2,i))
        END DO              
!      END DO
!    ENDIF
    
    ! y-derivatives on current element
!    IF (Bder(DER_DERIV_Y)) THEN
!      DO i=1,npoints
        DO idof = 1,9
          Dbas(idof,DER_DERIV_Y,i) = &
              -Dxj(i) * (Djac(3,i) * Dhelp(idof,1,i) - Djac(1,i) * Dhelp(idof,2,i))
        END DO
      END DO
!    ENDIF
  ENDIF
    
  END SUBROUTINE 

  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_Q2_sim (ieltyp, Dcoords, Djac, Ddetj, &
                               Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q2.
  INTEGER, INTENT(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE,nelements), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
  REAL(DP), DIMENSION(EL_NJACENTRIES2D,npoints,nelements), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  REAL(DP), DIMENSION(npoints,nelements), INTENT(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(NDIM2D,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i'th 
  !   basis function of the finite element k in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,EL_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  REAL(DP), DIMENSION(:,:,:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>

  ! auxiliary vector containing the first derivatives on the reference element
  REAL(DP), DIMENSION(9,NDIM2D,npoints) :: Dhelp
  REAL(DP) :: dx,dy
  REAL(DP),DIMENSION(npoints) :: Dxj !auxiliary variable
  
  INTEGER :: i   ! point counter
  INTEGER :: j   ! element counter
  INTEGER :: idof

  REAL(DP), PARAMETER :: Q2 = 0.5_DP
  REAL(DP), PARAMETER :: Q4 = 0.25_DP
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  IF (Bder(DER_FUNC)) THEN
  
    DO j=1,nelements
    
      DO i=1,npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)
        Dbas(1,DER_FUNC,i,j)= Q4*(1.0_DP-dx)*(1.0_DP-dy)*dx*dy
        Dbas(2,DER_FUNC,i,j)=-Q4*(1.0_DP+dx)*(1.0_DP-dy)*dx*dy
        Dbas(3,DER_FUNC,i,j)= Q4*(1.0_DP+dx)*(1.0_DP+dy)*dx*dy
        Dbas(4,DER_FUNC,i,j)=-Q4*(1.0_DP-dx)*(1.0_DP+dy)*dx*dy
        Dbas(5,DER_FUNC,i,j)=-Q2*(1.0_DP-dx*dx)*(1.0_DP-dy)*dy
        Dbas(6,DER_FUNC,i,j)= Q2*(1.0_DP+dx)*(1.0_DP-dy*dy)*dx
        Dbas(7,DER_FUNC,i,j)= Q2*(1.0_DP-dx*dx)*(1.0_DP+dy)*dy
        Dbas(8,DER_FUNC,i,j)=-Q2*(1.0_DP-dx)*(1.0_DP-dy*dy)*dx
        Dbas(9,DER_FUNC,i,j)= (1.0_DP-dx*dx)*(1.0_DP-dy*dy)
      END DO
      
    END DO
    
  END IF
    
  !if x-or y-derivatives are desired
  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN
  
    DO j=1,nelements
      Dxj = 1.0E0_DP / Ddetj(:,j)
      
      !x- and y-derivatives on reference element
      DO i=1,npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)

        !x- and y-derivatives on reference element
        Dhelp(1,1,i)= Q4*(1.0_DP-2.0_DP*dx)*(1.0_DP-dy)*dy
        Dhelp(2,1,i)=-Q4*(1.0_DP+2.0_DP*dx)*(1.0_DP-dy)*dy
        Dhelp(3,1,i)= Q4*(1.0_DP+2.0_DP*dx)*(1.0_DP+dy)*dy
        Dhelp(4,1,i)=-Q4*(1.0_DP-2.0_DP*dx)*(1.0_DP+dy)*dy
        Dhelp(5,1,i)= (1.0_DP-dy)*dx*dy
        Dhelp(6,1,i)= Q2*(1.0_DP+2.0_DP*dx)*(1.0_DP-dy*dy)
        Dhelp(7,1,i)=-(1.0_DP+dy)*dx*dy
        Dhelp(8,1,i)=-Q2*(1.0_DP-2.0_DP*dx)*(1.0_DP-dy*dy)
        Dhelp(9,1,i)=-2.0_DP*(1.0_DP-dy*dy)*dx

        Dhelp(1,2,i)= Q4*(1.0_DP-dx)*(1.0_DP-2.0_DP*dy)*dx
        Dhelp(2,2,i)=-Q4*(1.0_DP+dx)*(1.0_DP-2.0_DP*dy)*dx
        Dhelp(3,2,i)= Q4*(1.0_DP+dx)*(1.0_DP+2.0_DP*dy)*dx
        Dhelp(4,2,i)=-Q4*(1.0_DP-dx)*(1.0_DP+2.0_DP*dy)*dx
        Dhelp(5,2,i)=-Q2*(1.0_DP-dx*dx)*(1.0_DP-2.0_DP*dy)
        Dhelp(6,2,i)=-(1.0_DP+dx)*dx*dy
        Dhelp(7,2,i)= Q2*(1.0_DP-dx*dx)*(1.0_DP+2.0_DP*dy)
        Dhelp(8,2,i)= (1.0_DP-dx)*dx*dy
        Dhelp(9,2,i)=-2.0_DP*(1.0_DP-dx*dx)*dy
      END DO
        
      !x-derivatives on current element
!      IF (Bder(DER_DERIV_X)) THEN
        DO i=1,npoints
          DO idof = 1,9
            Dbas(idof,DER_DERIV_X,i,j) = &
                Dxj(i) * (Djac(4,i,j) * Dhelp(idof,1,i) &
                          - Djac(2,i,j) * Dhelp(idof,2,i))
          END DO              
        END DO
!      ENDIF
      
      !y-derivatives on current element
!      IF (Bder(DER_DERIV_Y)) THEN
        DO i=1,npoints
          DO idof = 1,9
            Dbas(idof,DER_DERIV_Y,i,j) = &
                -Dxj(i) * (Djac(3,i,j) * Dhelp(idof,1,i) &
                - Djac(1,i,j) * Dhelp(idof,2,i))
          END DO
        END DO
!      ENDIF

    END DO
      
  END IF
    
  END SUBROUTINE 

!**************************************************************************
! Element subroutines for nonparametric Q1~ element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  PURE SUBROUTINE elem_EM30 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                             Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point. The coordinates are expected
  ! on the real element!
!</description>

  !<input>

  ! Element type identifier. Must be =EL_EM30.
  INTEGER, INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(EL_NJACENTRIES2D), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  ! REMARK: Not used by this special type of element!
  REAL(DP), INTENT(IN) :: ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Cartesian coordinates of the evaluation point on the real element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
  REAL(DP), DIMENSION(2), INTENT(IN) :: Dpoint
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i'th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>

  ! This element clearly works only with standard quadrilaterals
  INTEGER, PARAMETER :: NVE = 4

  ! auxiliary variables  
  INTEGER :: IVE,IA,IK
  REAL(DP) :: PXL,PYL,PXU,PYU, D1,D2
  REAL(DP),DIMENSION(4) :: DXM,DYM,DLX,DLY,CKH,F
  REAL(DP),DIMENSION(4,4) :: A       ! local 4x4 system
  REAL(DP),DIMENSION(4,4) :: CK
  REAL(DP),DIMENSION(256) :: Dlapack
  REAL(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC3
  REAL(DP),DIMENSION(EL_MAXNBAS,EL_MAXNCOF) :: COB
  INTEGER :: INFO
  INTEGER, DIMENSION(NVE) :: IPIV
  REAL(DP) :: dx, dy
  
  INTERFACE
   
    ! Declare the two LAPACK subroutines as PURE using an interface.
    ! Then we can use them here.
  
    PURE SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
    USE fsystem
    INTEGER, INTENT(IN) :: LDA, M, N
    INTEGER, DIMENSION(N), INTENT(OUT) :: IPIV
    INTEGER, INTENT(OUT) :: INFO
    REAL(DP),INTENT(INOUT) :: A(LDA,N)
    END SUBROUTINE
    
    PURE SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
    USE fsystem
    INTEGER, INTENT(IN) :: LDA, N, LWORK
    INTEGER, DIMENSION(N), INTENT(OUT) :: IPIV
    INTEGER, INTENT(OUT) :: INFO
    REAL(DP),INTENT(INOUT) :: A(LDA,N),WORK(256)
    END SUBROUTINE

    PURE SUBROUTINE INVERT(A,F,X,IPAR)
    USE fsystem
    REAL(DP), DIMENSION(4,4), INTENT(INOUT) :: A
    REAL(DP), DIMENSION(4), INTENT(INOUT) :: F,X
    INTEGER, INTENT(IN) :: IPAR
    END SUBROUTINE

  END INTERFACE
  
  ! Clear the output array
  !Dbas = 0.0_DP
  
  dx = Dpoint(1)
  dy = Dpoint(2)
  
  DO IVE=1,NVE
    DXM(IVE)=0.5_DP*(Dcoords(1,IVE)+Dcoords(1,MOD(IVE,4)+1))
    DYM(IVE)=0.5_DP*(Dcoords(2,IVE)+Dcoords(2,MOD(IVE,4)+1))
    DLX(IVE)=0.5_DP*(Dcoords(1,MOD(IVE,4)+1)-Dcoords(1,IVE))
    DLY(IVE)=0.5_DP*(Dcoords(2,MOD(IVE,4)+1)-Dcoords(2,IVE))
  END DO

  D1 = 1.0_DP / SQRT((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
  D2 = 1.0_DP / SQRT((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)
  CA1 = (DXM(2)-DXM(4)) * D1
  CB1 = (DYM(2)-DYM(4)) * D1
  CA2 = (DXM(3)-DXM(1)) * D2
  CB2 = (DYM(3)-DYM(1)) * D2
  CA3 = CA1**2-CA2**2
  CB3 = 2.0_DP*(CA1*CB1-CA2*CB2)
  CC3 = CB1**2-CB2**2
  
  DO IA = 1,NVE
    PXL = DXM(IA)-SQRT(1.0_DP/3.0_DP)*DLX(IA)
    PYL = DYM(IA)-SQRT(1.0_DP/3.0_DP)*DLY(IA)
    PXU = DXM(IA)+SQRT(1.0_DP/3.0_DP)*DLX(IA)
    PYU = DYM(IA)+SQRT(1.0_DP/3.0_DP)*DLY(IA)
    A(IA,1)=0.5_DP*( F1(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                  +F1(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
    A(IA,2)=0.5_DP*( F2(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                  +F2(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
    A(IA,3)=0.5_DP*( F3(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                  +F3(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
    A(IA,4)=0.5_DP*( F4(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                  +F4(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
  END DO
  
  CALL INVERT(A,F,CKH,0)

  ! Invert the matrix with LAPACK
  !CALL DGETRF (4,4,A,4,IPIV,INFO)
  !CALL DGETRI (4,A,4,IPIV, Dlapack, 256, INFO )
  
!  DO IK1=1,4
!    DO IK2=1,4
!      CK(IK1,IK2)=A(IK2,IK1)
!    END DO
!  END DO
  CK = TRANSPOSE(A)

  DO IK=1,4
    COB(IK,1) = CK(IK,4)*CA3
    COB(IK,2) = CK(IK,4)*CC3
    COB(IK,3) = CK(IK,4)*CB3
    COB(IK,4) = CK(IK,2)*CA1+CK(IK,3)*CA2
    COB(IK,5) = CK(IK,2)*CB1+CK(IK,3)*CB2
    COB(IK,6) = CK(IK,1)
  END DO

  ! Function values

  IF (BDER(DER_FUNC)) THEN
    Dbas(1,DER_FUNC)= COB(1,1)*dx**2+COB(1,2)*dy**2+COB(1,3)*dx*dy &
                    +COB(1,4)*dx   +COB(1,5)*dy   +COB(1,6)
    Dbas(2,DER_FUNC)= COB(2,1)*dx**2+COB(2,2)*dy**2+COB(2,3)*dx*dy &
                    +COB(2,4)*dx   +COB(2,5)*dy   +COB(2,6)
    Dbas(3,DER_FUNC)= COB(3,1)*dx**2+COB(3,2)*dy**2+COB(3,3)*dx*dy &
                    +COB(3,4)*dx   +COB(3,5)*dy   +COB(3,6)
    Dbas(4,DER_FUNC)= COB(4,1)*dx**2+COB(4,2)*dy**2+COB(4,3)*dx*dy &
                    +COB(4,4)*dx   +COB(4,5)*dy   +COB(4,6)
  END IF

  ! Derivatives:
  
  IF (BDER(DER_DERIV_X) .OR. BDER(DER_DERIV_Y)) THEN
    ! x-derivatives
               
    Dbas(1,DER_DERIV_X) = 2.0_DP*COB(1,1)*dx+COB(1,3)*dy+COB(1,4)
    Dbas(2,DER_DERIV_X) = 2.0_DP*COB(2,1)*dx+COB(2,3)*dy+COB(2,4)
    Dbas(3,DER_DERIV_X) = 2.0_DP*COB(3,1)*dx+COB(3,3)*dy+COB(3,4)
    Dbas(4,DER_DERIV_X) = 2.0_DP*COB(4,1)*dx+COB(4,3)*dy+COB(4,4)

    ! y-derivatives
          
    Dbas(1,DER_DERIV_Y) = 2.0_DP*COB(1,2)*dy+COB(1,3)*dx+COB(1,5)
    Dbas(2,DER_DERIV_Y) = 2.0_DP*COB(2,2)*dy+COB(2,3)*dx+COB(2,5)
    Dbas(3,DER_DERIV_Y) = 2.0_DP*COB(3,2)*dy+COB(3,3)*dx+COB(3,5)
    Dbas(4,DER_DERIV_Y) = 2.0_DP*COB(4,2)*dy+COB(4,3)*dx+COB(4,5)
  END IF
             
  CONTAINS

    ! Auxiliary functions
      
    ELEMENTAL REAL(DP) FUNCTION F1(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    REAL(DP), INTENT(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F1 = 1.0_DP
    END FUNCTION
    
    ELEMENTAL REAL(DP) FUNCTION F2(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    REAL(DP), INTENT(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F2 = CA1*X  +CB1*Y
    END FUNCTION
    
    ELEMENTAL REAL(DP) FUNCTION F3(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    REAL(DP), INTENT(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F3=CA2*X  +CB2*Y
    END FUNCTION
    
    ELEMENTAL REAL(DP) FUNCTION F4(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    REAL(DP), INTENT(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F4 = CA3*X*X+CB3*X*Y+CC3*Y*Y
    END FUNCTION

  END SUBROUTINE 
  
  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_EM30_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given points. The coordinates are expected
  ! on the real element!
!</description>

!<input>
  ! Element type identifier. Must be =EL_EM30.
  INTEGER, INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(EL_NJACENTRIES2D,npoints), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points.
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(npoints), INTENT(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the reference element.
  ! DIMENSION(NDIM2D,npoints)
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>

  ! This element clearly works only with standard quadrilaterals
  INTEGER, PARAMETER :: NVE = 4

  ! auxiliary variables  
  INTEGER :: IVE,IA,IK,i
  REAL(DP) :: PXL,PYL,PXU,PYU, D1,D2
  REAL(DP),DIMENSION(4) :: DXM,DYM,DLX,DLY,CKH,F
  REAL(DP),DIMENSION(4,4) :: A       ! local 4x4 system
  REAL(DP),DIMENSION(4,4) :: CK
  REAL(DP),DIMENSION(256) :: Dlapack
  REAL(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC3
  REAL(DP),DIMENSION(EL_MAXNBAS,EL_MAXNCOF) :: COB
  INTEGER :: INFO
  INTEGER, DIMENSION(NVE) :: IPIV
  
  INTERFACE
   
    ! Declare the two LAPACK subroutines as PURE using an interface.
    ! Then we can use them here.
  
    PURE SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
    INTEGER, INTENT(IN) :: LDA, M, N
    INTEGER, DIMENSION(N), INTENT(OUT) :: IPIV
    INTEGER, INTENT(OUT) :: INFO
    DOUBLE PRECISION,INTENT(INOUT) :: A(LDA,N)
    END SUBROUTINE
    
    PURE SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
    INTEGER, INTENT(IN) :: LDA, N, LWORK
    INTEGER, DIMENSION(N), INTENT(OUT) :: IPIV
    INTEGER, INTENT(OUT) :: INFO
    DOUBLE PRECISION,INTENT(INOUT) :: A(LDA,N),WORK(256)
    END SUBROUTINE

    PURE SUBROUTINE INVERT(A,F,X,IPAR)
    USE fsystem
    REAL(DP), DIMENSION(4,4), INTENT(INOUT) :: A
    REAL(DP), DIMENSION(4), INTENT(INOUT) :: F,X
    INTEGER, INTENT(IN) :: IPAR
    END SUBROUTINE

  END INTERFACE
  
  ! Clear the output array
  !Dbas = 0.0_DP
  
  DO IVE=1,NVE
    DXM(IVE)=0.5_DP*(Dcoords(1,IVE)+Dcoords(1,MOD(IVE,4)+1))
    DYM(IVE)=0.5_DP*(Dcoords(2,IVE)+Dcoords(2,MOD(IVE,4)+1))
    DLX(IVE)=0.5_DP*(Dcoords(1,MOD(IVE,4)+1)-Dcoords(1,IVE))
    DLY(IVE)=0.5_DP*(Dcoords(2,MOD(IVE,4)+1)-Dcoords(2,IVE))
  END DO

  D1 = 1.0_DP / SQRT((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
  D2 = 1.0_DP / SQRT((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)
  CA1 = (DXM(2)-DXM(4)) * D1
  CB1 = (DYM(2)-DYM(4)) * D1
  CA2 = (DXM(3)-DXM(1)) * D2
  CB2 = (DYM(3)-DYM(1)) * D2
  CA3 = CA1**2-CA2**2
  CB3 = 2.0_DP*(CA1*CB1-CA2*CB2)
  CC3 = CB1**2-CB2**2
  
  DO IA = 1,4
    PXL = DXM(IA)-SQRT(1.0_DP/3.0_DP)*DLX(IA)
    PYL = DYM(IA)-SQRT(1.0_DP/3.0_DP)*DLY(IA)
    PXU = DXM(IA)+SQRT(1.0_DP/3.0_DP)*DLX(IA)
    PYU = DYM(IA)+SQRT(1.0_DP/3.0_DP)*DLY(IA)
    A(IA,1)=0.5_DP*( F1(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                  +F1(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
    A(IA,2)=0.5_DP*( F2(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                  +F2(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
    A(IA,3)=0.5_DP*( F3(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                  +F3(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
    A(IA,4)=0.5_DP*( F4(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                  +F4(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
  END DO
  
  CALL INVERT(A,F,CKH,0)

  ! Invert the matrix with LAPACK
  !CALL DGETRF (4,4,A,4,IPIV,INFO)
  !CALL DGETRI (4,A,4,IPIV, Dlapack, 256, INFO )
  
!  DO IK1=1,4
!    DO IK2=1,4
!      CK(IK1,IK2)=A(IK2,IK1)
!    END DO
!  END DO
  CK = TRANSPOSE(A)

  DO IK=1,4
    COB(IK,1) = CK(IK,4)*CA3
    COB(IK,2) = CK(IK,4)*CC3
    COB(IK,3) = CK(IK,4)*CB3
    COB(IK,4) = CK(IK,2)*CA1+CK(IK,3)*CA2
    COB(IK,5) = CK(IK,2)*CB1+CK(IK,3)*CB2
    COB(IK,6) = CK(IK,1)
  END DO

  ! Function values
  
  IF (BDER(DER_FUNC)) THEN
    DO i=1,npoints
      Dbas(1,DER_FUNC,i)= COB(1,1)*Dpoints(1,i)**2 &
                +COB(1,2)*Dpoints(2,i)**2 &
                +COB(1,3)*Dpoints(1,i)*Dpoints(2,i) &
                +COB(1,4)*Dpoints(1,i)   &
                +COB(1,5)*Dpoints(2,i)   &
                +COB(1,6)
      Dbas(2,DER_FUNC,i)= COB(2,1)*Dpoints(1,i)**2 &
                +COB(2,2)*Dpoints(2,i)**2 &
                +COB(2,3)*Dpoints(1,i)*Dpoints(2,i) &
                +COB(2,4)*Dpoints(1,i)   &
                +COB(2,5)*Dpoints(2,i)   &
                +COB(2,6)
      Dbas(3,DER_FUNC,i)= COB(3,1)*Dpoints(1,i)**2 &
                +COB(3,2)*Dpoints(2,i)**2 &
                +COB(3,3)*Dpoints(1,i)*Dpoints(2,i) &
                +COB(3,4)*Dpoints(1,i)   &
                +COB(3,5)*Dpoints(2,i)   &
                +COB(3,6)
      Dbas(4,DER_FUNC,i)= COB(4,1)*Dpoints(1,i)**2 &
                +COB(4,2)*Dpoints(2,i)**2 &
                +COB(4,3)*Dpoints(1,i)*Dpoints(2,i) &
                +COB(4,4)*Dpoints(1,i)    &
                +COB(4,5)*Dpoints(2,i)    &
                +COB(4,6)
    END DO
  END IF

  ! Derivatives:
  
  IF (BDER(DER_DERIV_X) .OR. BDER(DER_DERIV_Y)) THEN
    ! x-derivatives
         
    DO i=1,npoints      
      Dbas(1,DER_DERIV_X,i) = 2.0_DP*COB(1,1)*Dpoints(1,i)+COB(1,3)*Dpoints(2,i)+COB(1,4)
      Dbas(2,DER_DERIV_X,i) = 2.0_DP*COB(2,1)*Dpoints(1,i)+COB(2,3)*Dpoints(2,i)+COB(2,4)
      Dbas(3,DER_DERIV_X,i) = 2.0_DP*COB(3,1)*Dpoints(1,i)+COB(3,3)*Dpoints(2,i)+COB(3,4)
      Dbas(4,DER_DERIV_X,i) = 2.0_DP*COB(4,1)*Dpoints(1,i)+COB(4,3)*Dpoints(2,i)+COB(4,4)
  !  END DO

    ! y-derivatives
          
  !  DO i=1,npoints
      Dbas(1,DER_DERIV_Y,i) = 2.0_DP*COB(1,2)*Dpoints(2,i)+COB(1,3)*Dpoints(1,i)+COB(1,5)
      Dbas(2,DER_DERIV_Y,i) = 2.0_DP*COB(2,2)*Dpoints(2,i)+COB(2,3)*Dpoints(1,i)+COB(2,5)
      Dbas(3,DER_DERIV_Y,i) = 2.0_DP*COB(3,2)*Dpoints(2,i)+COB(3,3)*Dpoints(1,i)+COB(3,5)
      Dbas(4,DER_DERIV_Y,i) = 2.0_DP*COB(4,2)*Dpoints(2,i)+COB(4,3)*Dpoints(1,i)+COB(4,5)
    END DO
  END IF
             
  CONTAINS

    ! Auxiliary functions
      
    ELEMENTAL REAL(DP) FUNCTION F1(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    REAL(DP), INTENT(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F1 = 1.0_DP
    END FUNCTION
    
    ELEMENTAL REAL(DP) FUNCTION F2(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    REAL(DP), INTENT(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F2 = CA1*X  +CB1*Y
    END FUNCTION
    
    ELEMENTAL REAL(DP) FUNCTION F3(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    REAL(DP), INTENT(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F3=CA2*X  +CB2*Y
    END FUNCTION
    
    ELEMENTAL REAL(DP) FUNCTION F4(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    REAL(DP), INTENT(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F4 = CA3*X*X+CB3*X*Y+CC3*Y*Y
    END FUNCTION

  END SUBROUTINE 

  !************************************************************************
  
!<subroutine>  

  SUBROUTINE elem_EM30_sim (ieltyp, Dcoords, Djac, Ddetj, &
                            Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_EM30.
  INTEGER, INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE,nelements), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
  REAL(DP), DIMENSION(EL_NJACENTRIES2D,npoints,nelements), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  REAL(DP), DIMENSION(npoints,nelements), INTENT(IN) :: Ddetj
  
  ! Derivative quantifier array. array [1..EL_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  LOGICAL, DIMENSION(EL_MAXNDER), INTENT(IN) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected on the real element.
  ! DIMENSION(NDIM2D,npoints,nelements)
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i'th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,EL_MAXNDER,npoints,nelements), INTENT(OUT) :: Dbas
  REAL(DP), DIMENSION(:,:,:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>

  ! This element clearly works only with standard quadrilaterals
  INTEGER, PARAMETER :: NVE = 4

  ! auxiliary variables  
  INTEGER :: IVE,IA,IK,i,j
  REAL(DP) :: PXL,PYL,PXU,PYU, D1,D2
  REAL(DP),DIMENSION(4) :: DXM,DYM,DLX,DLY , CKH,F
  REAL(DP), DIMENSION(4,4) :: A       ! local 4x4 system
  REAL(DP), DIMENSION(4,4) :: CK
  REAL(DP) :: Dlapack(256)
  REAL(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC3,X,Y
  REAL(DP) :: COB(EL_MAXNBAS,EL_MAXNCOF)
  INTEGER :: INFO
  INTEGER, DIMENSION(NVE) :: IPIV
  INTEGER :: npointsfunc,npointsderx,npointsdery
  
  !INCLUDE 'cbasictria.inc'
  !INCLUDE 'cbasicelem.inc'
  !INCLUDE 'celem.inc'
  
  INTERFACE
   
    ! Declare the two LAPACK subroutines as PURE using an interface.
    ! Then we can use them here.
  
    PURE SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
    INTEGER, INTENT(IN) :: LDA, M, N
    INTEGER, DIMENSION(N), INTENT(OUT) :: IPIV
    INTEGER, INTENT(OUT) :: INFO
    DOUBLE PRECISION,INTENT(INOUT) :: A(LDA,N)
    END SUBROUTINE
    
    PURE SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
    INTEGER, INTENT(IN) :: LDA, N, LWORK
    INTEGER, DIMENSION(N), INTENT(OUT) :: IPIV
    INTEGER, INTENT(OUT) :: INFO
    DOUBLE PRECISION,INTENT(INOUT) :: A(LDA,N),WORK(256)
    END SUBROUTINE
    
    PURE SUBROUTINE INVERT(A,F,X,IPAR)
    USE fsystem
    REAL(DP), DIMENSION(4,4), INTENT(INOUT) :: A
    REAL(DP), DIMENSION(4), INTENT(INOUT) :: F,X
    INTEGER, INTENT(IN) :: IPAR
    END SUBROUTINE

  END INTERFACE
  
  ! Calculate the loop counters in advance. Help us to get rid
  ! of any if-commands in the element loop below.
  npointsfunc = 0
  npointsderx = 0
  npointsdery = 0
  IF (Bder(DER_FUNC)) npointsfunc = npoints
  IF (Bder(DER_DERIV_X)) npointsderx = npoints
  IF (Bder(DER_DERIV_Y)) npointsdery = npoints
  
  ! Clear the output array
  !Dbas = 0.0_DP
  
  ! Loop over the elements
  
  DO j=1,nelements
    
!    DO ive=1,4
!      DX(ive) = Dcoords(1,IVE,j)
!      DY(ive) = Dcoords(2,IVE,j)
!    END DO
!    BDER = Bder1
!    
!    CALL EM30(0.0_DP,0.0_DP,-2)
!    
!    DO ive=1,npoints
!      CALL EM30(Dpoints(1,ive,j),Dpoints(2,ive,j),0)
!      DBas1(:,1:3,1,1) = DBAS(:,1:3)
!      !DBas1(:,1:3,ive,j) = DBAS(:,1:3)
!    END DO
!    
!  
    DO IVE=1,NVE
      DXM(IVE)=0.5_DP*(Dcoords(1,IVE,j)+Dcoords(1,MOD(IVE,4)+1,j))
      DYM(IVE)=0.5_DP*(Dcoords(2,IVE,j)+Dcoords(2,MOD(IVE,4)+1,j))
      DLX(IVE)=0.5_DP*(Dcoords(1,MOD(IVE,4)+1,j)-Dcoords(1,IVE,j))
      DLY(IVE)=0.5_DP*(Dcoords(2,MOD(IVE,4)+1,j)-Dcoords(2,IVE,j))
    END DO

    D1 = 1.0_DP / SQRT((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
    D2 = 1.0_DP / SQRT((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)
    CA1 = (DXM(2)-DXM(4)) * D1
    CB1 = (DYM(2)-DYM(4)) * D1
    CA2 = (DXM(3)-DXM(1)) * D2
    CB2 = (DYM(3)-DYM(1)) * D2
    CA3 = CA1**2-CA2**2
    CB3 = 2.0_DP*(CA1*CB1-CA2*CB2)
    CC3 = CB1**2-CB2**2
    
    DO IA = 1,4
      PXL = DXM(IA)-SQRT(1.0_DP/3.0_DP)*DLX(IA)
      PYL = DYM(IA)-SQRT(1.0_DP/3.0_DP)*DLY(IA)
      PXU = DXM(IA)+SQRT(1.0_DP/3.0_DP)*DLX(IA)
      PYU = DYM(IA)+SQRT(1.0_DP/3.0_DP)*DLY(IA)
      A(IA,1)=0.5_DP*( F1(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                    +F1(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
      A(IA,2)=0.5_DP*( F2(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                    +F2(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
      A(IA,3)=0.5_DP*( F3(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                    +F3(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
      A(IA,4)=0.5_DP*( F4(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                    +F4(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
    END DO
    
    CALL INVERT(A,F,CKH,0)

    ! Invert the matrix with LAPACK
    !CALL DGETRF (4,4,A,4,IPIV,INFO)
    !CALL DGETRI (4,A,4,IPIV, Dlapack, 256, INFO )
    
!    DO IK1=1,4
!      DO IK2=1,4
!        CK(IK1,IK2)=A(IK2,IK1)
!      END DO
!    END DO
    CK = TRANSPOSE(A)

    DO IK=1,4
      COB(IK,1) = CK(IK,4)*CA3
      COB(IK,2) = CK(IK,4)*CC3
      COB(IK,3) = CK(IK,4)*CB3
      COB(IK,4) = CK(IK,2)*CA1+CK(IK,3)*CA2
      COB(IK,5) = CK(IK,2)*CB1+CK(IK,3)*CB2
      COB(IK,6) = CK(IK,1)
    END DO

    ! Function values
    DO i=1,npointsfunc   ! either 0 or npoints

      Dbas(1,DER_FUNC,i,j)= COB(1,1)*Dpoints(1,i,j)**2 &
                +COB(1,2)*Dpoints(2,i,j)**2 &
                +COB(1,3)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                +COB(1,4)*Dpoints(1,i,j)   &
                +COB(1,5)*Dpoints(2,i,j)   &
                +COB(1,6)
      Dbas(2,DER_FUNC,i,j)= COB(2,1)*Dpoints(1,i,j)**2 &
                +COB(2,2)*Dpoints(2,i,j)**2 &
                +COB(2,3)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                +COB(2,4)*Dpoints(1,i,j)   &
                +COB(2,5)*Dpoints(2,i,j)   &
                +COB(2,6)
      Dbas(3,DER_FUNC,i,j)= COB(3,1)*Dpoints(1,i,j)**2 &
                +COB(3,2)*Dpoints(2,i,j)**2 &
                +COB(3,3)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                +COB(3,4)*Dpoints(1,i,j)   &
                +COB(3,5)*Dpoints(2,i,j)   &
                +COB(3,6)
      Dbas(4,DER_FUNC,i,j)= COB(4,1)*Dpoints(1,i,j)**2 &
                +COB(4,2)*Dpoints(2,i,j)**2 &
                +COB(4,3)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                +COB(4,4)*Dpoints(1,i,j)    &
                +COB(4,5)*Dpoints(2,i,j)    &
                +COB(4,6)
    END DO
    
    ! x-derivatives
          
    DO i=1,npointsderx   ! either 0 or npoints
      Dbas(1,DER_DERIV_X,i,j) = 2.0_DP*COB(1,1)*Dpoints(1,i,j) &
                                      +COB(1,3)*Dpoints(2,i,j)+COB(1,4)
      Dbas(2,DER_DERIV_X,i,j) = 2.0_DP*COB(2,1)*Dpoints(1,i,j) &
                                      +COB(2,3)*Dpoints(2,i,j)+COB(2,4)
      Dbas(3,DER_DERIV_X,i,j) = 2.0_DP*COB(3,1)*Dpoints(1,i,j) &
                                      +COB(3,3)*Dpoints(2,i,j)+COB(3,4)
      Dbas(4,DER_DERIV_X,i,j) = 2.0_DP*COB(4,1)*Dpoints(1,i,j) &
                                      +COB(4,3)*Dpoints(2,i,j)+COB(4,4)
    END DO

    ! y-derivatives
          
    DO i=1,npointsdery   ! either 0 or npoints
      Dbas(1,DER_DERIV_Y,i,j) = 2.0_DP*COB(1,2)*Dpoints(2,i,j) &
                                      +COB(1,3)*Dpoints(1,i,j)+COB(1,5)
      Dbas(2,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,2)*Dpoints(2,i,j) &
                                      +COB(2,3)*Dpoints(1,i,j)+COB(2,5)
      Dbas(3,DER_DERIV_Y,i,j) = 2.0_DP*COB(3,2)*Dpoints(2,i,j) &
                                      +COB(3,3)*Dpoints(1,i,j)+COB(3,5)
      Dbas(4,DER_DERIV_Y,i,j) = 2.0_DP*COB(4,2)*Dpoints(2,i,j) &
                                      +COB(4,3)*Dpoints(1,i,j)+COB(4,5)

    END DO
    
  END DO ! ielement
             
  CONTAINS

    ! Auxiliary functions
      
    ELEMENTAL REAL(DP) FUNCTION F1(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    REAL(DP), INTENT(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F1 = 1.0_DP
    END FUNCTION
    
    ELEMENTAL REAL(DP) FUNCTION F2(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    REAL(DP), INTENT(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F2 = CA1*X  +CB1*Y
    END FUNCTION
    
    ELEMENTAL REAL(DP) FUNCTION F3(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    REAL(DP), INTENT(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F3=CA2*X  +CB2*Y
    END FUNCTION
    
    ELEMENTAL REAL(DP) FUNCTION F4(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    REAL(DP), INTENT(IN) :: X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3
    F4 = CA3*X*X+CB3*X*Y+CC3*Y*Y
    END FUNCTION

  END SUBROUTINE 

END MODULE 
