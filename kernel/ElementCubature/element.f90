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
!# 5.) elem_igetDimension
!#     -> Determine the dimension of an element; 1D, 2D, 3D,...\\
!#
!# 6.) elem_isNonparametric
!#     -> Check whether an element is parametric or nonparametric\\
!#
!# 7.) elem_getPrimaryElement
!#     -> Returns the element identifier of the 'primary' element without
!#        any subtype information, which identifies an element family.\\
!#
!# 8.) elem_generic 
!#     -> Realises a generic element which can be used to evaluate a finite 
!#        element depending on its element identifier - in contrast to the 
!#        standard evaluation routines, which ignore the element quantifier 
!#        as they 'know' what they are...\\
!#
!# 9.) elem_generic_mult
!#     -> The multiple-point-evaluation routine for a generic element.\\
!#
!# 10.) elem_generic_sim
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
!# 2.) What does the variable 'ieltyp' mean, which must be passed to
!#   the element subroutines?
!#
!#   This variable identifies the type of the element. For example
!#   ieltyp=EL_Q1 identifies a Q1-element.
!#
!# 3.) I wrote a subroutine which works for a special element family. More 
!#   precisely for example, my routine works for all variants of $\tilde Q_1$.
!#   Do I really have to check all element types at the beginning of the
!#   routine? That would be legthy!?! Like in
!#
!#    SUBROUTINE abc (ieltype)
!#    ...
!#      IF ( (ieltype .NE. EL_E030) .AND. (ieltype .NE. EL_EM30) ...) THEN
!#        CALL sys_halt()
!#      END IF
!#
!#   No you don't need to. The simplest method for such element families:
!#   Use elem_getPrimaryElement! This returns an element identifier
!#   identifying all elements 'of the same type'. The above can be
!#   shortened this way to:
!#
!#    SUBROUTINE abc (ieltype)
!#    ...
!#      IF (elem_getPrimaryElement(ieltype) .NE. EL_Q1T) THEN
!#        CALL sys_halt()
!#      END IF
!#
!#   Which directly describes all $\tilde Q_1$ elements. If you additionally
!#   want to check it whether it is parametric/nonparametric, use an
!#   additional check with elem_isNonparametric!
!#
!#   Of course, for 'standard' elements, elem_getPrimaryElement does nothing.
!#   So writing for $Q_1$ for example,
!#
!#      IF (elem_getPrimaryElement(ieltype) .NE. EL_Q1) THEN
!#
!#   is exactly the same as
!#
!#      IF (ieltype .NE. EL_Q1) THEN
!#
!#   but the first version is cleaner :-)
!#
!# 4.) How is this thing with the 'primary' element realised? I mean, for example, 
!#   EL_EM30 is defined as EL_EM30=EL_Q1T+2**8!?! What does this mean?
!#
!#   The element constants follow a special pattern oriented on a bitfield to 
!#   encode as much information into an integer as possible.
!#   Every element identifier consists of a 32 bit integer, which is coded as
!#   follows:
!#
!#% Bit | 31 ... 24 23 ... 16 | 15                | 14 ... 10 | 9   8 | 7 ............ 0|
!#% -------------------------------------------------------------------------------------
!#%     |         ****        | =1: nonparametric | unused    |dimens |  Element number |
!#
!#   Bits 0..7   specifies the element number/family (1=P1, 11=Q1,...).
!#   Bits 8+ 9   encodes the dimension of the element. =1: 1D, =2: 2D, =3: 3D.
!#   Bit    15   specifies whether an element is nonparametric (1) or parametric (0).
!#   Bits 16..31 encodes special variants of elements. This is used only for some special
!#               type elements:
!#               Q1T: Bit 16 =1: element nonconformal, integral mean value based, 
!#                           =0: element conformal, midpoint value based
!#                    Bit 17:=0: Standard handling
!#                           =1: For nonparametric element: Don't use pivoting when 
!#                               solving local 4x4 systems; faster but less stable on 
!#                               cruel elements.
!#                    Bit 18:=0: Standard handling
!#                           =1: No scaling of the local coodinate system for 
!#                               nonparametric $\tilde Q_1$ element. Slightly faster but 
!#                               less stable on cruel elements.
!#               Q1HN: Bit 16/17: These three bits encode the local number of the edge
!#                               where there is a hanging node.
!#                               =0: local edge 1, =1: local edge 2,...
!#               P2ISO: Bit 16/17: These three bits encode how many edges of the $P_2$
!#                               element are to be handled with an isoparametric mapping
!#                               to the reference element. 
!#                               =0: isoparametric mapping with one edge not linear
!#                               =1: isoparametric mapping with two edges not linear
!#                               =2: isoparametric mapping with three edges not linear
!#                    Bit 18/19/20/21: These four bits encode which of the edges are
!#                               to be mapped nonlinear from the reference to the real
!#                               element.
!#                               Bit 18:=1 1st edge maps nonlinear
!#                               Bit 19:=1 2nd edge maps nonlinear
!#                               Bit 20:=1 3rd edge maps nonlinear
!#               Q2ISO: Bit 16/17: These three bits encode how many edges of the $Q_2$
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
!#   To obtain the actual element identifier, one must mask out all bits 
!#   except for bit 0..7, i.e. to check whether an element is Q1~, one can use
!#   elem_getPrimaryElement, which masks out all unimportant bits with IAND:
!#
!#     if ( elem_getPrimaryElement(ieltype) .EQ. EL_Q1T ) then ...
!#
!#   or to check all variants
!#
!#     if ( (ieltype .eq. EL_E030) .OR. (ieltype .eq. EL_EM30). OR. ...) then ...
!# 
!#   When it's clear that it's a $\tilde Q_1$ element, one can have a closer
!#   look which variant it is -- if this is necessary.
!#
!#   Final note: For implementational reasons, the bits 8+9+16..31 coincide --
!#   depending on the element -- with the ID of the transformation formula
!#   from the reference to the real element! This allows an easier and more
!#   transparent implementation of isoparametric elements!
!#
!# </purpose>
!##############################################################################

MODULE element

  USE fsystem
  USE basicgeometry
  USE derivatives
  USE transformation
  USE mprimitives

  IMPLICIT NONE
!<constants>

!<constantblock description="Internal constants for element ID bitfield.">
  ! Bitmasks for dimension; coincides on purpose with TRAFO_DIM_DIMENSION!
  INTEGER(I32), PARAMETER :: EL_DIMENSION = 2**8 + 2**9
  
  ! 1D element; coincides on purpose with TRAFO_DIM_1D!
  INTEGER(I32), PARAMETER :: EL_1D = ISHFT(NDIM1D,8)

  ! 2D element; coincides on purpose with TRAFO_DIM_2D!
  INTEGER(I32), PARAMETER :: EL_2D = ISHFT(NDIM2D,8)
  
  ! 3D element; coincides on purpose with TRAFO_DIM_3D!
  INTEGER(I32), PARAMETER :: EL_3D = ISHFT(NDIM3D,8)
  
  ! Bitmask for element number including dimension
  INTEGER(I32), PARAMETER :: EL_ELNRMASK = 255 + EL_DIMENSION
  
  ! Bitmask specifying a nonparametric element
  INTEGER(I32), PARAMETER :: EL_NONPARAMETRIC = 2**15
!</constantblock>

!<constantblock description="Element identifiers for 1D elements">

  ! ID of constant conforming line FE, P0
  INTEGER(I32), PARAMETER :: EL_P0_1D = EL_1D + 0
  
  ! ID of linear conforming line FE, P1
  INTEGER(I32), PARAMETER :: EL_P1_1D = EL_1D + 1
  
  ! ID of quadratic conforming line FE, P2
  INTEGER(I32), PARAMETER :: EL_P2_1D = EL_1D + 2
  
!</constantblock>
  
!<constantblock description="Element identifiers for 2D elements">

  ! unspecified element
  INTEGER(I32), PARAMETER :: EL_UNDEFINED = -1

  ! ID of constant conforming triangular FE, P0 (just for the FEAST-users...)
  INTEGER(I32), PARAMETER :: EL_P0   = EL_2D + 0

  ! ID of constant conforming triangular FE, P0
  INTEGER(I32), PARAMETER :: EL_E000 = EL_P0

  ! ID of linear conforming triangular FE, P1 (just for the FEAST-users...)
  INTEGER(I32), PARAMETER :: EL_P1   = EL_2D + 1

  ! ID of linear conforming triangular FE, P1 
  INTEGER(I32), PARAMETER :: EL_E001 = EL_P1

  ! ID of quadratic conforming triangular FE, P2 (just for the FEAST-users...)
  INTEGER(I32), PARAMETER :: EL_P2   = EL_2D + 2

  ! ID of quadratic conforming triangular FE, P2 
  INTEGER(I32), PARAMETER :: EL_E002 = EL_P2

  ! ID of cubic conforming triangular FE, P3 (just for the FEAST-users...)
  INTEGER(I32), PARAMETER :: EL_P3   = EL_2D + 3

  ! ID of cubic conforming triangular FE, P3
  INTEGER(I32), PARAMETER :: EL_E003 = EL_P3

  ! ID of constant conforming quadrilateral FE, Q0 (just for the FEAST-users...)
  INTEGER(I32), PARAMETER :: EL_Q0   = EL_2D + 10

  ! ID of constant conforming quadrilateral FE, Q0
  INTEGER(I32), PARAMETER :: EL_E010 = EL_Q0

  ! ID of bilinear conforming quadrilateral FE, Q1 (just for the FEAST-users...)
  INTEGER(I32), PARAMETER :: EL_Q1   = EL_2D + 11

  ! ID of bilinear conforming quadrilateral FE, Q1 
  INTEGER(I32), PARAMETER :: EL_E011 = EL_Q1 

  ! ID of biquadratic conforming quadrilateral FE, Q2 (just for the FEAST-users...)
  INTEGER(I32), PARAMETER :: EL_Q2   = EL_2D + 13

  ! ID of biquadratic conforming quadrilateral FE, Q2 
  INTEGER(I32), PARAMETER :: EL_E013 = EL_Q2

  ! ID of bicubic conforming quadrilateral FE, Q3 (just for the FEAST-users...)
  INTEGER(I32), PARAMETER :: EL_Q3   = EL_2D + 14

  ! ID of bicubic conforming quadrilateral FE, Q3 
  INTEGER(I32), PARAMETER :: EL_E014 = EL_Q3
  
  ! ID of nonconforming parametric linear P1 element on a quadrilareral
  ! element, given by function value in the midpoint and the two
  ! derivatives.
  INTEGER(I32), PARAMETER :: EL_QP1  = EL_2D + 21

  ! General rotated bilinear $\tilde Q1$ element, all variants (conformal, 
  ! nonconformal, parametric, nonparametric).
  ! Simplest variant is: parametric, edge midpoint-value based.
  INTEGER(I32), PARAMETER :: EL_Q1T  = EL_2D + 30
  
  ! Quadrilateral $Q_1$ element with one hanging node. In the property
  ! bitfield of the element identifier, the local number of the edge where there
  ! is a hanging node must be encoded! (i.e. one must add a corresponding
  ! value to the constant EL_Q1HN1 to get the actual element identifier!)
  INTEGER(I32), PARAMETER :: EL_Q1HN1 = EL_2D + 40
  
  ! Quadrilateral $Q_2$ element with isoparametric mapping from the reference
  ! to the real element. In the property bitfield, one must set the corresponding
  ! bits to identify the edge that should map isoparametric!
  INTEGER(I32), PARAMETER :: EL_Q2ISO = EL_2D + 50

!</constantblock>

!<constantblock description="Element identifiers for 3D elements">

  ! ID of constant conforming tetrahedral FE, P0
  INTEGER(I32), PARAMETER :: EL_P0_3D = EL_3D + 0
  
  ! ID of linear conforming tetrahedral FE, P1
  INTEGER(I32), PARAMETER :: EL_P1_3D = EL_3D + 1

  ! ID of constant conforming hexahedral FE, Q0
  INTEGER(I32), PARAMETER :: EL_Q0_3D = EL_3D + 10

  ! ID of trilinear conforming hexahedral FE, Q1
  INTEGER(I32), PARAMETER :: EL_Q1_3D = EL_3D + 11

!</constantblock>

!<constantblock description="Special element variants.">

  ! ID of rotated bilinear conforming quadrilateral FE, Q1~, integral
  ! mean value based
  INTEGER(I32), PARAMETER :: EL_E030 = EL_Q1T + 2**16

  ! ID of rotated bilinear conforming quadrilateral FE, Q1~, edge-midpoint based
  INTEGER(I32), PARAMETER :: EL_E031 = EL_Q1T

  ! ID of rotated bilinear nonconforming quadrilateral FE, Q1~, integral
  ! mean value based
  INTEGER(I32), PARAMETER :: EL_EM30 = EL_Q1T + EL_NONPARAMETRIC+ 2**16

  ! ID of rotated bilinear nonconforming quadrilateral FE, Q1~, integral
  ! mean value based; 'unpivoted' variant, solving local 4x4 systems directly 
  ! without pivoting. Faster but less stable.
  INTEGER(I32), PARAMETER :: EL_EM30_UNPIVOTED = EL_Q1T + EL_NONPARAMETRIC + 2**17 + 2**16
  
  ! ID of rotated bilinear nonconforming quadrilateral FE, Q1~, integral
  ! mean value based; 'unscaled' variant, which does not scale the local
  ! coordinate system on every element. Faster but less stable.
  INTEGER(I32), PARAMETER :: EL_EM30_UNSCALED = EL_Q1T + EL_NONPARAMETRIC + 2**18 + 2**16

  ! ID of rotated bilinear nonconforming quadrilateral FE, Q1~, edge-midpoint based
  INTEGER(I32), PARAMETER :: EL_EM31 = EL_Q1T + EL_NONPARAMETRIC

  ! Isoparametric $Q_2$ element with one edge mapped nonlinear from the reference
  ! to the real element. Additionally, one bit in the property bitfield must
  ! be set to identify the edge.
  INTEGER(I32), PARAMETER :: EL_Q2ISO1 = EL_Q2ISO 

  ! Isoparametric $Q_2$ element with two edges mapped nonlinear from the reference
  ! to the real element. Additionally, two bits in the property bitfield must
  ! be set to identify the edges.
  INTEGER(I32), PARAMETER :: EL_Q2ISO2 = EL_Q2ISO + 2**16

  ! Isoparametric $Q_2$ element with three edges mapped nonlinear from the reference
  ! to the real element. Additionally, three bits in the property bitfield must
  ! be set to identify the edges.
  INTEGER(I32), PARAMETER :: EL_Q2ISO3 = EL_Q2ISO + 2**17

  ! Isoparametric $Q_2$ element with four edges mapped nonlinear from the reference
  ! to the real element. Additionally, four bits in the property bitfield must
  ! be set to identify the edges.
  INTEGER(I32), PARAMETER :: EL_Q2ISO4 = EL_Q2ISO + 2**16 + 2**17

!</constantblock>

!<constantblock description="maximal values">

  ! Maximum number of basic functions = maximum number of
  ! local DOF's per element.
  INTEGER, PARAMETER :: EL_MAXNBAS = 21

  ! Maximum different derivative descriptor types supported by elements;
  ! corresponds to the maximum value of DER_xxxx constants in the 
  ! module 'derivatives'. Can be used as first quantifier in the
  ! array given to evaluation subroutines which receives function values/
  ! derivatives/...
  INTEGER, PARAMETER :: EL_MAXNDER = 6
  
  INTEGER, PARAMETER :: EL_MAXNCOF = 6
  
  ! Number of entries in the Jacobian matrix, defining the mapping between
  ! the reference element and the real element. 1x1-matrix=1 elements
  INTEGER, PARAMETER :: EL_NJACENTRIES1D = 1

  ! Number of entries in the Jacobian matrix, defining the mapping between
  ! the reference element and the real element. 2x2-matrix=4 elements
  INTEGER, PARAMETER :: EL_NJACENTRIES2D = 4

  ! Number of entries in the Jacobian matrix, defining the mapping between
  ! the reference element and the real element. 3x3-matrix=9 elements
  INTEGER, PARAMETER :: EL_NJACENTRIES3D = 9

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
  INTEGER(I32), INTENT(IN) :: ieltype

!</input>

!<result>
  ! The number of local DOF's on this element.
!</result>

!</function>

  SELECT CASE (elem_getPrimaryElement(ieltype))
  
  ! -= 1D element types =-
  CASE (EL_P0_1D)
    ! local DOF's for 1D P0
    elem_igetNDofLoc = 1
  CASE (EL_P1_1D)
    ! local DOF's for 1D P1
    elem_igetNDofLoc = 2
  CASE (EL_P2_1D)
    ! local DOF's for 1D P2
    elem_igetNDofLoc = 3
  
  ! -= 2D element types =-
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
  CASE (EL_QP1)
    ! local DOF's for QP1
    elem_igetNDofLoc = 3
  CASE (EL_Q1T)
    ! local DOF's for Ex30
    elem_igetNDofLoc = 4
  
  ! -= 3D element types =-
  CASE (EL_P0_3D, EL_Q0_3D)
    ! local DOF's for 3D P0, Q0
    elem_igetNDofLoc = 1
  CASE (EL_P1_3D)
    ! local DOF's for 3D P1
    elem_igetNDofLoc = 4
  CASE (EL_Q1_3D)
    ! local DOF's for 3D Q1
    elem_igetNDofLoc = 8
    
  CASE DEFAULT
    elem_igetNDofLoc = 0
  END SELECT

  END FUNCTION

  ! ***************************************************************************

!<subroutine>  

  PURE SUBROUTINE elem_igetNDofLocAssignment(ieltype, &
      ndofAtVertices, ndofAtEdges, ndofAtFaces, ndofAtElement)

!<description>
  ! This function returns for a given element type the number of local
  ! degrees of freedom that is assigned to vertices, edges, faces and elements.
!</description>

!<input>    
  ! The element type identifier.
  INTEGER(I32), INTENT(IN) :: ieltype
!</input>

!<output>
  ! Number of DOF's assigned to the vertices on one element.
  INTEGER, INTENT(OUT) :: ndofAtVertices
  
  ! Number of DOF's assigned to the edges on one element.
  ! Is always 0 for 1D element types.
  INTEGER, INTENT(OUT) :: ndofAtEdges
  
  ! Number of DOF's assigned to the faces on one element.
  ! Is always 0 for 1D/2D element types.
  INTEGER, INTENT(OUT) :: ndofAtFaces
  
  ! Number of DOF's assigned to one element, which do not belong to
  ! vertices or edges.
  INTEGER, INTENT(OUT) :: ndofAtElement

!</output>

!</subroutine>

  ! Default setup
  ndofAtVertices = 0
  ndofAtEdges = 0
  ndofAtFaces = 0
  ndofAtElement = 0

  SELECT CASE (elem_getPrimaryElement(ieltype))
  
  ! -= 1D element types =-
  CASE (EL_P0_1D)
    ! local DOF's for P0
    ndofAtElement  = 1
  CASE (EL_P1_1D)
    ! local DOF's for P1
    ndofAtVertices = 2
  CASE (EL_P2_1D)
    ! local DOF's for P2
    ndofAtVertices = 2
    ndofAtElement  = 1
    
  ! -= 2D element types =-
  CASE (EL_P0, EL_Q0)
    ! local DOF's for Q0
    ndofAtElement  = 1
  CASE (EL_P1)
    ! local DOF's for P1
    ndofAtVertices = 3
  CASE (EL_P2)
    ! local DOF's for P2
    ndofAtVertices = 3
    ndofAtEdges    = 3
  CASE (EL_P3)
    ! local DOF's for P3
    ndofAtVertices = 3
    ndofAtEdges    = 6
  CASE (EL_Q1)
    ! local DOF's for Q1
    ndofAtVertices = 4
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
  CASE (EL_QP1)
    ! local DOF's for QP1
    ndofAtElement  = 3
  CASE (EL_Q1T)
    ! local DOF's for Ex30
    ndofAtEdges    = 4
  
  ! -= 3D element types =-
  CASE (EL_P0_3D, EL_Q0_3D)
    ! local DOF's for P0
    ndofAtElement  = 1
  CASE (EL_P1_3D)
    ! local DOF's for P1
    ndofAtVertices = 4
  CASE (EL_Q1_3D)
    ! local DOF's for Q1
    ndofAtVertices = 8
  END SELECT

  END SUBROUTINE

  ! ***************************************************************************

!<function>  

  ELEMENTAL INTEGER FUNCTION elem_igetNVE(ieltype)

!<description>
  ! This function returns for a given element type the number of vertices
  ! in the corresponding element primitive, i.e. 3 for triangular and
  ! 4 for quadrilateral elements.
!</description>

!<input>    

  ! The element type identifier.
  INTEGER(I32), INTENT(IN) :: ieltype

!</input>

!<result>
  ! The number vertices in the element primitive/shape where the finite
  ! element is defined on.
!</result>

!</function>

  ! NVE coincides with that NVE from the transformation.
  ! Use the transformation routine to determine that value!
  
  elem_igetNVE = trafo_igetNVE(elem_igetTrafoType(ieltype))

  !SELECT CASE (elem_getPrimaryElement(ieltype))
  !CASE (EL_P0,EL_P1,EL_P2,EL_P3)
  !  elem_igetNVE = 3
  !CASE (EL_Q0, EL_Q1, EL_Q2, EL_Q3, EL_QP1, EL_Q1T)
  !  elem_igetNVE = 4
  !CASE DEFAULT
  !  elem_igetNVE = 0
  !END SELECT

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
  INTEGER(I32), INTENT(IN) :: ieltype

!</input>

!<result>
  ! The type of the coordinate system the element ieltype acts on. One of
  ! the TRAFO_CS_xxxx constants from the transformation module.
!</result>

!</function>

  SELECT CASE (IAND(ieltype,EL_ELNRMASK+EL_NONPARAMETRIC))
  CASE (EL_P0_1D,EL_P1_1D,EL_P2_1D)
    ! Line elements
    elem_igetCoordSystem = TRAFO_CS_REF1D
  CASE (EL_P0,EL_P1,EL_P2,EL_P3)
    ! Triangular elements work in barycentric coordinates
    elem_igetCoordSystem = TRAFO_CS_BARY2DTRI
  CASE (EL_Q0,EL_Q1,EL_Q2,EL_Q3,EL_QP1,EL_Q1T)
    ! These work on the reference quadrilateral
    elem_igetCoordSystem = TRAFO_CS_REF2DQUAD
  CASE (EL_Q1T+EL_NONPARAMETRIC)
    ! EM30, EM31; these work in real coordinates
    elem_igetCoordSystem = TRAFO_CS_REAL2DQUAD
  CASE (EL_P0_3D, EL_P1_3D)
    ! Tetrahedral elements work in barycentric coordinates
    elem_igetCoordSystem = TRAFO_CS_BARY3DTETRA
  CASE (EL_Q0_3D, EL_Q1_3D)
    ! These work on the reference hexahedron
    elem_igetCoordSystem = TRAFO_CS_REF3DHEXA
  CASE DEFAULT
    elem_igetCoordSystem = TRAFO_CS_UNDEFINED
  END SELECT
  
  END FUNCTION

  ! ***************************************************************************

!<function>  

  ELEMENTAL INTEGER(I32) FUNCTION elem_igetTrafoType(ieltype)

!<description>
  ! This function returns the typical type of transformation, a specific 
  ! element uses for transferring coordinates of the reference element to 
  ! the real element.
!</description>

!<input>    
  ! The element type identifier.
  INTEGER(I32), INTENT(IN) :: ieltype
!</input>

!<result>
  ! The maximum derivative of that element. <= MAX_NDER.
!</result>

!</function>

  SELECT CASE (elem_getPrimaryElement(ieltype))
  
  CASE (EL_P0_1D, EL_P1_1D, EL_P2_1D)
    ! Linear line transformation, 1D
    elem_igetTrafoType = TRAFO_ID_LINSIMPLEX + TRAFO_DIM_1D
  
  CASE (EL_P0, EL_P1, EL_P2, EL_P3)
    ! Linear triangular transformation, 2D
    elem_igetTrafoType = TRAFO_ID_LINSIMPLEX + TRAFO_DIM_2D
    
  CASE (EL_Q0,EL_Q1,EL_Q2,EL_Q3,EL_QP1,EL_Q1T)
    ! Bilinear quadrilateral transformation, 2D.
    elem_igetTrafoType = TRAFO_ID_MLINCUBE + TRAFO_DIM_2D
  
  CASE (EL_P0_3D, EL_P1_3D)
    ! Linear tetrahedral transrormation, 3D
    elem_igetTrafoType = TRAFO_ID_LINSIMPLEX + TRAFO_DIM_3D
  
  CASE (EL_Q0_3D, EL_Q1_3D)
    ! Trilinear hexahedral transformation, 3D
    elem_igetTrafoType = TRAFO_ID_MLINCUBE + TRAFO_DIM_3D
  
  CASE DEFAULT
    elem_igetTrafoType = TRAFO_ID_UNKNOWN
    
  END SELECT

  END FUNCTION

  ! ***************************************************************************

!<function>  

  ELEMENTAL INTEGER FUNCTION elem_igetDimension(ieltype)

!<description>
  ! This function returns the dimensional constant that specifies which
  ! dimension (1D, 2D,...) an element uses.
!</description>

!<input>    
  ! The element type identifier.
  INTEGER(I32), INTENT(IN) :: ieltype
!</input>

!<result>
  ! A constant that specifies the dimension of an element. NDIM2 for 2D,
  ! NDIM3D for 3D,...
!</result>

!</function>

    ! The dimension is encoded in two bits in the element quantifier!
    elem_igetDimension = ISHFT(IAND(ieltype,EL_DIMENSION),-8)

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
  INTEGER(I32), INTENT(IN) :: ieltype

!</input>

!<result>
  ! One of the TRAFO_IDxxxx constants of the module 'transformation' identifying
  ! the type of transformation.
!</result>

!</function>

  SELECT CASE (elem_getPrimaryElement(ieltype))
  CASE (EL_P0_1D)
    ! Function
    elem_getMaxDerivative = 1
  CASE (EL_P1_1D)
    ! Function + 1st derivative
    elem_getMaxDerivative = 2
  CASE (EL_P2_1D)
    ! Function + 1st + 2nd derivative
    elem_getMaxDerivative = 3
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
  CASE (EL_QP1)
    ! Function + 1st derivative
    elem_getMaxDerivative = 3
  CASE (EL_Q1T)
    ! Function + 1st derivative
    elem_getMaxDerivative = 3
  CASE (EL_P0_3D, EL_Q0_3D)
    ! Function + 1st derivative
    elem_getMaxDerivative = 1
  CASE (EL_P1_3D)
    ! Function + 1st derivative
    elem_getMaxDerivative = 4
  CASE (EL_Q1_3D)
    ! Function + 1st derivative
    elem_getMaxDerivative = 4
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
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates,
  ! Dcoords(3,.)=z-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. The size of this array depends
  ! on the dimension of the space.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  REAL(DP), DIMENSION(:), INTENT(IN) :: Djac
  
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
  CASE (EL_P0_1D)
    CALL elem_P0_1D (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_P1_1D)
    CALL elem_P1_1D (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_P2_1D)
    CALL elem_P2_1D (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
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
  CASE (EL_QP1)
    CALL elem_QP1 (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_EM30,EL_EM30_UNPIVOTED,EL_EM30_UNSCALED)
    CALL elem_EM30 (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_E030)
    CALL elem_E030 (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_EM31)
    CALL elem_EM31 (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_E031)
    CALL elem_E031 (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  !CASE (EL_P0_3D)
  !  CALL elem_P0_3D (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  !CASE (EL_P1_3D)
  !  CALL elem_P1_3D (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_Q0_3D)
    CALL elem_Q0_3D (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_Q1_3D)
    CALL elem_Q1_3D (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
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
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  !REAL(DP), DIMENSION(:,npoints), INTENT(IN) :: Djac
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Djac
  
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
  ! DIMENSION(#space dimensions,npoints)
  !   Dpoints(1,.)=x-coordinates,
  !   Dpoints(2,.)=y-coordinates.
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

  ! local variables
  INTEGER :: i

  ! Choose the right element subroutine to call.
  SELECT CASE (ieltyp)
  CASE (EL_P0_1D)
    CALL elem_P0_1D_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_P1_1D)
    CALL elem_P1_1D_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_P2_1D)
    CALL elem_P2_1D_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
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
  CASE (EL_EM30,EL_EM30_UNPIVOTED,EL_EM30_UNSCALED)
    CALL elem_EM30_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_E030)
    CALL elem_E030_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_EM31)
    CALL elem_EM31_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_E031)
    CALL elem_E031_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  !CASE (EL_P0_3D)
  !  CALL elem_P0_3D_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  !CASE (EL_P1_3D)
  !  CALL elem_P1_3D_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_Q0_3D)
    CALL elem_Q0_3D_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_Q1_3D)
    CALL elem_Q1_3D_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
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
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
  !REAL(DP), DIMENSION(:,npoints,nelements), INTENT(IN) :: Djac
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Djac
  
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
  CASE (EL_P0_1D)
    CALL elem_P0_1D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                         Bder, Dbas, npoints, nelements, Dpoints)
  CASE (EL_P1_1D)
    CALL elem_P1_1D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                         Bder, Dbas, npoints, nelements, Dpoints)
  CASE (EL_P2_1D)
    CALL elem_P2_1D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                         Bder, Dbas, npoints, nelements, Dpoints)
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
  CASE (EL_EM30,EL_EM30_UNPIVOTED,EL_EM30_UNSCALED)
    CALL elem_EM30_sim (ieltyp, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
  CASE (EL_E030)
    CALL elem_E030_sim (ieltyp, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
  CASE (EL_EM31)
    CALL elem_EM31_sim (ieltyp, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
  CASE (EL_E031)
    CALL elem_E031_sim (ieltyp, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
  !CASE (EL_P0_3D)
  !  CALL elem_P0_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
  !                    Bder, Dbas, npoints, nelements, Dpoints)
  !CASE (EL_P1_3D)
  !  CALL elem_P1_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
  !                    Bder, Dbas, npoints, nelements, Dpoints)
  CASE (EL_Q0_3D)
    CALL elem_Q0_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                      Bder, Dbas, npoints, nelements, Dpoints)
  CASE (EL_Q1_3D)
    CALL elem_Q1_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
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

  ELEMENTAL LOGICAL FUNCTION elem_isNonparametric (ieltyp) RESULT (inonpar)

  !<description>
  
  ! Determines whether an element ieltyp is parametric or nonparametric.
  
  !</description>
  
  !<result>
  ! =true, if the element is nonparametric. 
  ! =false, if the element is parametric.
  !</result>
  
  !<input>
  
  ! Element type qualifier.
  INTEGER(I32), INTENT(IN) :: ieltyp
  
  !</input>
 
!</function>
 
    inonpar = IAND(ieltyp,EL_NONPARAMETRIC) .NE. 0
  
  END FUNCTION

  !************************************************************************
  
!<function>  

  ELEMENTAL INTEGER(I32) FUNCTION elem_getPrimaryElement (ieltyp) RESULT (iresult)

!<description>
  ! Determines the 'primary' element type identifier. For standard elements
  ! like $Q_1$, there is usually ieltyp=elem_getPrimaryElement(ieltyp).
  ! But some elements like $\tilde Q_1$ have different modifications,
  ! which are encoded in ieltyp. In this case, elem_getPrimaryElement
  ! returns the element identifier of the element without any
  ! modification; this can be seen as 'standard' or 'primary' element
  ! type.
  ! In case of $\tilde Q_1$ for example, there is
  !    elem_getPrimaryElement(EL_E030) = elem_getPrimaryElement(EL_EM30)
  !  = elem_getPrimaryElement(EL_E031) = elem_getPrimaryElement(EL_EM31)
  !  = EL_Q1T,
  ! which identifies the 'general' $\tilde Q_1$ element.
!</description>
  
!<result>
  ! The identifier of the 'standard' element, ieltyp refers to.
!</result>
  
  !<input>
  
  ! Element type qualifier.
   INTEGER(I32), INTENT(IN) :: ieltyp
  
  !</input>
 
!</function>
 
    ! To get the standard identifier, we just have to mask out all bits
    ! except for bit 0..9. These 10 bits encode the standard element
    ! identifier plus dimension.
    iresult = IAND(ieltyp,EL_ELNRMASK)
  
  END FUNCTION

!**************************************************************************
! Element subroutines for parametric P0_1D element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  PURE SUBROUTINE elem_P0_1D (ieltyp, Dcoords, Djac, ddetj, Bder, &
                              Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_P0_1D.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(EL_NJACENTRIES1D), INTENT(IN) :: Djac
  
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
  REAL(DP), DIMENSION(1), INTENT(IN) :: Dpoint
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
    
  ! P0_1D is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC) = 1.0_DP

  END SUBROUTINE 
  
  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_P0_1D_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                   Bder, Dbas, npoints, Dpoints)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_P0_1D.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(EL_NJACENTRIES1D,npoints), INTENT(IN) :: Djac
  
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
  ! DIMENSION(#space dimensions,npoints)
  ! Dpoints(1,.)=x-coordinates,
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
  
  ! P0_1D is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC,:) = 1.0_DP

  END SUBROUTINE 

  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_P0_1D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, nelements, Dpoints)

  !<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_P0_1D.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE^,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(EL_NJACENTRIES1D,npoints), INTENT(IN) :: Djac
  
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
  ! DIMENSION(#space dimensions,npoints,nelements)
  !  Dpoints(1,.)=x-coordinates,
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
  
  ! P0_1D is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC,:,:) = 1.0_DP

  END SUBROUTINE 

!**************************************************************************
! Element subroutines for parametric P1_1D element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  PURE SUBROUTINE elem_P1_1D (ieltyp, Dcoords, Djac, ddetj, Bder, &
                             Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_P1_1D.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  REAL(DP), DIMENSION(EL_NJACENTRIES1D), INTENT(IN) :: Djac
  
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
  REAL(DP), DIMENSION(1), INTENT(IN) :: Dpoint
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

  ! The P1_1D element is specified by two polynomials on the reference element.
  ! These two polynomials are:
  !
  !  P1(X,Y) = 1/2 (1-x)
  !  P2(X,Y) = 1/2 (1+x)
  !
  ! Each of them calculated that way that Pi(Xj)=delta_ij (Kronecker)
  ! for X1,X2 the two corners of the reference element [-1,1].
  
  ! Clear the output array
  !Dbas = 0.0_DP
  
  ! Remark: The P1_1D-element always computes function value and 1st derivative.
  ! That's even faster than when using two IF commands for preventing
  ! the computation of one of the values!
      
  ! If function values are desired, calculate them.
!  if (el_bder(DER_FUNC)) then
    Dbas(1,DER_FUNC) = 0.5_DP * (1.0_DP - Dpoint(1))
    Dbas(2,DER_FUNC) = 0.5_DP * (1.0_DP + Dpoint(1))
!  endif
  
  ! If x-derivatives are desired, calculate them.
  ! Since P1_1D is linear, the first derivative is constant!
!  if (Bder(DER_DERIV_X)) then
    Dbas(1,DER_DERIV_X) = -0.5_DP / ddetj
    Dbas(2,DER_DERIV_X) = 0.5_DP / ddetj
!  endif
    
  END SUBROUTINE 
  
  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_P1_1D_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                   Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_P1_1D.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  REAL(DP), DIMENSION(EL_NJACENTRIES1D,npoints), INTENT(IN) :: Djac
  
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
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
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

  INTEGER :: i   ! point counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The P1_1D-element always computes function value and 1st derivatives
  ! That's even faster than when using two IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  !IF (Bder(DER_FUNC)) THEN
    DO i=1,npoints
      Dbas(1,DER_FUNC,i) = 0.5_DP * (1_DP - Dpoints(1,i))
      Dbas(2,DER_FUNC,i) = 0.5_DP * (1_DP + Dpoints(1,i))
    END DO
  !ENDIF
  
  !if x-derivatives are desired
!  IF (Bder(DER_DERIV_X)) THEN
    DO i=1,npoints
      Dbas(1,DER_DERIV_X,i) = -0.5_DP / Ddetj(i)
      Dbas(2,DER_DERIV_X,i) = 0.5_DP / Ddetj(i)
    END DO
!  ENDIF
    
  END SUBROUTINE 

  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_P1_1D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_P1_1D.
  INTEGER(I32), INTENT(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements).
  !  Dcoords(1,.,.)=x-coordinates,
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
  REAL(DP), DIMENSION(EL_NJACENTRIES1D,npoints,nelements), INTENT(IN) :: Djac
  
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
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
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

  INTEGER :: i   ! point counter
  INTEGER :: j   ! element counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  !if function values are desired
  IF (Bder(DER_FUNC)) THEN
  
    DO j=1,nelements
    
      DO i=1,npoints
        Dbas(1,DER_FUNC,i,j) = 0.5_DP * (1.0_DP - Dpoints(1,i,j))
        Dbas(2,DER_FUNC,i,j) = 0.5_DP * (1.0_DP + Dpoints(1,i,j))
      END DO
      
    END DO
    
  END IF
    
  !if x-derivatives are desired
  IF (Bder(DER_DERIV_X)) THEN
  
    DO j=1,nelements
    
      DO i=1,npoints
        Dbas(1,DER_DERIV_X,i,j) = -0.5 / Ddetj(i,j)
        Dbas(2,DER_DERIV_X,i,j) = 0.5 / Ddetj(i,j)
      END DO

    END DO
      
  END IF
    
  END SUBROUTINE 
  
!**************************************************************************
! Element subroutines for parametric P2_1D element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  PURE SUBROUTINE elem_P2_1D (ieltyp, Dcoords, Djac, ddetj, Bder, &
                             Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_P2_1D.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  REAL(DP), DIMENSION(EL_NJACENTRIES1D), INTENT(IN) :: Djac
  
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
  REAL(DP), DIMENSION(1), INTENT(IN) :: Dpoint
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

  REAL(DP) :: d

  ! The P2_1D element is specified by three polynomials on the reference element.
  ! These three polynomials are:
  !
  !  P1(x) = 1/2*x*(x-1)
  !  P2(x) = 1/2*x*(x+1)
  !  P3(x) = 1-x*x
  !
  
  ! Clear the output array
  !Dbas = 0.0_DP
  
  ! Remark: The P2_1D-element always computes function value and 1st derivative.
  ! That's even faster than when using two IF commands for preventing
  ! the computation of one of the values!
      
  ! If function values are desired, calculate them.
!  if (el_bder(DER_FUNC)) then
    Dbas(1,DER_FUNC) = 0.5_DP * Dpoint(1) * (Dpoint(1) - 1.0_DP)
    Dbas(2,DER_FUNC) = 0.5_DP * Dpoint(1) * (Dpoint(1) + 1.0_DP)
    Dbas(3,DER_FUNC) = 1.0_DP - Dpoint(1)**2
!  endif
  
  ! If x-derivatives are desired, calculate them.
!  if (Bder(DER_DERIV_X)) then
    d = 1.0_DP / ddetj
    Dbas(1,DER_DERIV_X) = (Dpoint(1) - 0.5_DP) * d
    Dbas(2,DER_DERIV_X) = (Dpoint(1) + 0.5_DP) * d
    Dbas(3,DER_DERIV_X) = -2.0_DP * Dpoint(1) * d
!  endif
    
  END SUBROUTINE 
  
  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_P2_1D_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                   Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_P2_1D.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  REAL(DP), DIMENSION(EL_NJACENTRIES1D,npoints), INTENT(IN) :: Djac
  
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
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
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

  REAL(DP) :: d
  INTEGER :: i   ! point counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The P2_1D-element always computes function value and 1st derivatives
  ! That's even faster than when using two IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  !IF (Bder(DER_FUNC)) THEN
    DO i=1,npoints
      Dbas(1,DER_FUNC,i) = 0.5_DP * Dpoints(1,i) * (Dpoints(1,i) - 1.0_DP)
      Dbas(2,DER_FUNC,i) = 0.5_DP * Dpoints(1,i) * (Dpoints(1,i) + 1.0_DP)
      Dbas(3,DER_FUNC,i) = 1.0_DP - Dpoints(1,i)**2
    END DO
  !ENDIF
 
  !if x-derivatives are desired
!  IF (Bder(DER_DERIV_X)) THEN
    DO i=1,npoints
      d = 1.0_DP / Ddetj(i)
      Dbas(1,DER_DERIV_X,i) = (Dpoints(1,i) - 0.5_DP) * d
      Dbas(2,DER_DERIV_X,i) = (Dpoints(1,i) + 0.5_DP) * d
      Dbas(3,DER_DERIV_X,i) = -2.0_DP * Dpoints(1,i) * d
    END DO
!  ENDIF
    
  END SUBROUTINE 

  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_P2_1D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_P2_1D.
  INTEGER(I32), INTENT(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements).
  !  Dcoords(1,.,.)=x-coordinates,
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
  REAL(DP), DIMENSION(EL_NJACENTRIES1D,npoints,nelements), INTENT(IN) :: Djac
  
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
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
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

  REAL(DP) :: d
  INTEGER :: i   ! point counter
  INTEGER :: j   ! element counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  !if function values are desired
  IF (Bder(DER_FUNC)) THEN
  
    DO j=1,nelements
    
      DO i=1,npoints
        Dbas(1,DER_FUNC,i,j) = 0.5_DP*Dpoints(1,i,j)*(Dpoints(1,i,j) - 1.0_DP)
        Dbas(2,DER_FUNC,i,j) = 0.5_DP*Dpoints(1,i,j)*(Dpoints(1,i,j) + 1.0_DP)
        Dbas(3,DER_FUNC,i,j) = 1.0_DP - Dpoints(1,i,j)**2
      END DO
      
    END DO
    
  END IF
    
  !if x-derivatives are desired
  IF (Bder(DER_DERIV_X)) THEN
  
    DO j=1,nelements
    
      DO i=1,npoints
        d = 1.0_DP / Ddetj(i,j)
        Dbas(1,DER_DERIV_X,i,j) = (Dpoints(1,i,j) - 0.5_DP) * d
        Dbas(2,DER_DERIV_X,i,j) = (Dpoints(1,i,j) + 0.5_DP) * d
        Dbas(3,DER_DERIV_X,i,j) = -2.0_DP * Dpoints(1,i,j) * d
      END DO

    END DO
      
  END IF
    
  END SUBROUTINE 

!**************************************************************************
! General information: Function values and derivatives of 
!                      triangular elements
!                      with linear transformation from the reference
!                      to the real element.
!
! The element subroutines return
! - the function value and
! - the X- and Y-derivatives
! of the basis function in a (cubature) point (x,y) on the real mesh!
! The coordinates of a (cubature) point is given
! - as coordinate triple (xi1, xi2, xi3) on the reference element, if the
!   the element is parametric; the actual cubature point is then at
!   (x,y) = s(t(xi1,xi2,xi3))
! - as coordinate pair (x,y) on the real element, if the element is
!   nonparametric.
! The mapping  s=(s1,s2):R^2->R^2  is the bilinear mapping "sigma" from 
! transformation.f90, that maps the reference element T^ to the real
! element T; its shape is of no importance here.
! The transformation t:R^3->R^2 maps the coordinates from the barycenric
! coordinate system on the reference element to the standard space
! (to get the (X,Y) coordinate on the reference element), so
!
!    t(xi1,xi2,xi3) = xi2 * [1,0]  +  xi3 * [0,1]
!
! The linear transformation s(.) from the reference to the real element
! has the form
!
!    s(X,Y) = c1  +  c2 X  +  c3 Y  =:  (x,y)
!
! Let u be an arbitrary FE (basis) function on an element T and
! p be the associated polynomial on the reference element T^. Then,
! we can obtain the function value of u at (x,y) easily:
!
!   u(x,y) = u(s(t(xi1,xi2,xi3)) = p(t^-1(s^-1(x,y))) = p(xi1,xi2,xi3)
!
!   [0,1]
!   | \               s(t(.))                C
!   |   \           --------->              / \
!   |  T^ \                               /     \
!   |       \                           /    T    \
!   [0,0]----[1,0]                     A-----------B
!
! The derivative is a little bit harder to calculate because of the
! mapping. We have:
!
!    grad(u)(x,y) = grad( p(t^-1(s^-1(x,y))) )
!
! Let's use the notation 
!
!    P(X,Y)  :=  p( t^-1 (X,Y) )  =  p (xi1,xi2,xi3)
!
! which is the 'standard' form of the polynomials on the
! reference element without barycentric coordinates involved (i.e.
! for the P1-element it's of the form P(x,y) = c1 x + c2 y + c3).
!
! Because of the chain rule, we have:
!
!    grad( p( t^-1 (s^-1) ) ) = (DP)( (s^-1) * D(s^-1) )
!
!       = ( P_X(s^-1)  P_Y(s^-1)) * ( (s1^-1)_x   (s1^-1)_y )
!                                   ( (s2^-1)_x   (s2^-1)_y )
!
!      =: ( P_X(s^-1)  P_Y(s^-1)) * ( e f )
!                                   ( g h )
!
! With s^-1(x,y)=(X,Y), we therefore have:
!
!    grad(u)(x,y) = ( P_X(X,Y) * e  +  P_Y(X,Y) * g )
!                   ( P_X(X,Y) * f  +  P_Y(X,Y) * h )
!
! Now, from e.g. http://mathworld.wolfram.com/MatrixInverse.html we know,
! that:
! 
!     A = ( a b )    =>   A^-1  =  1/det(A) (  d -b )
!         ( c d )                           ( -c  a )
!
! Therefore:
!
!    ( e f )  =  D(s^-1)  =  (Ds)^-1  =  1/(ad-bc) (  d -b )
!    ( g h )                                       ( -c  a )
!
! with
!
!    A = ( a b )  =  ( s1_X   s1_Y )  =  ( B-A  C-A )
!        ( c d )     ( s2_X   s2_Y )
!
! being the matrix from the transformation (according to
! http://mathworld.wolfram.com/BarycentricCoordinates.html).
!
!**************************************************************************

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
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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
  
  ! We have no derivatives.

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
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelemens)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dcoords
  
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
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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
  
  ! The P1 space consists of 'linear' finite elements. We have three basis 
  ! functions on the reference element, which can be written down in
  ! standard coordinates (-> P(.)) as well as in barycentric coordinates
  ! (-> p(.)). These are:
  !
  !   p1(xi1,xi2,xi3) = xi1 =  1 - X - Y  = P1(X,Y)
  !   p2(xi1,xi2,xi3) = xi2 =  X          = P2(X,Y)
  !   p3(xi1,xi2,xi3) = xi3 =  Y          = P3(X,Y)
  
  ! Clear the output array
  !Dbas = 0.0_DP
    
  ! Remark: The P1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  ! If function values are desired, calculate them.
  ! Use the p(.) representation in barycentric coordinates to calculate the
  ! function values.
!  if (el_bder(DER_FUNC)) then
    Dbas(1:3,DER_FUNC) = Dpoint(1:3)
!  endif
  
  ! If x-or y-derivatives are desired, calculate them.
  ! Here, we use the P(.) representation to get P_X and P_Y (which are
  ! only 0, 1 or -1)!
  ! These are then multiplied with the inverse of the transformation
  ! as described above to get the actual values of the derivatives.
  
!  if ((el_bder(DER_DERIV_X)) .or. (el_bder(DER_DERIV_Y))) then
    dxj = 1E0_DP / ddetj
    
    ! x-derivatives on current element.
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
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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
  INTEGER(I32), INTENT(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dcoords
  
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
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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
  INTEGER(I32), INTENT(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dcoords
  
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
! General information: Function values and derivatives of 
!                      quadrilateral elements
!                      with bilinear transformation between the
!                      reference and the real element.
!
! The element subroutines return
! - the function value and
! - the X- and Y-derivatives
! of the basis function in a (cubature) point (x,y) on the real mesh!
! The coordinates of a (cubature) point is given
! - as coordinate pair (xi1, xi2) on the reference element, if the
!   the element is parametric; the actual cubature point is then at
!   (x,y) = s(xi1,xi2)
! - as coordinate pair (x,y) on the real element, if the element is
!   nonparametric.
! The mapping  s=(s1,s2):R^2->R^2  is the bilinear mapping "sigma" from 
! transformation.f90, that maps the reference element T^ to the real
! element T; its shape is of no importance here.
!
! Let u be an arbitrary FE (basis) function on an element T and
! p be the associated polynomial on the reference element T^. Then,
! we can obtain the function value of u at (x,y) easily:
!
!   u(x,y) = u(s(xi1,xi2)) = p(s^-1(x,y)) = p(xi1,xi2)
!
! The derivative is a little bit harder to calculate because of the
! mapping. We have:
!
!    grad(u)(x,y) = grad( p(s^-1(x,y)) )
!
! Because of the chain rule, we have:
!
!    grad( p(s^-1) ) = (Dp)(s^-1) * D(s^-1)
!
!       = ( p_xi1(s^-1)  p_xi2(s^-1)) * ( (s1^-1)_x   (s1^-1)_y )
!                                       ( (s2^-1)_x   (s2^-1)_y )
!
!      =: ( p_xi1(s^-1)  p_xi2(s^-1)) * ( e f )
!                                       ( g h )
!
! With s^-1(x,y)=(xi1,xi2), we therefore have:
!
!    grad(u)(x,y) = ( p_xi1(xi1,xi2) * e  +  p_xi2(xi1,xi2) * g )
!                   ( p_xi1(xi1,xi2) * f  +  p_xi2(xi1,xi2) * h )
!
! Now, from e.g. http://mathworld.wolfram.com/MatrixInverse.html we know,
! that:
! 
!     A = ( a b )    =>   A^-1  =  1/det(A) (  d -b )
!         ( c d )                           ( -c  a )
!
! Therefore:
!
!    ( e f )  =  D(s^-1)  =  (Ds)^-1  =  1/(ad-bc) (  d -b )
!    ( g h )                                       ( -c  a )
!
! with
!
!    A = ( a b )  =  ( s1_x   s1_y )
!        ( c d )     ( s2_x   s2_y )
!
! being the matrix from the transformation.
!**************************************************************************

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
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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
  ! DIMENSION(#space dimensions,npoints)
  ! Dpoints(1,.)=x-coordinates,
  ! Dpoints(2,.)=y-coordinates.
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

  PURE SUBROUTINE elem_Q0_sim (ieltyp, Dcoords, Djac, Ddetj, &
                               Bder, Dbas, npoints, nelements, Dpoints)

  !<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q0.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE^,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dcoords
  
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
  ! DIMENSION(#space dimensions,npoints,nelements)
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
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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
  REAL(DP), DIMENSION(4,NDIM2D) :: Dhelp
  REAL(DP) :: dx,dy

  REAL(DP) :: dxj !auxiliary variable
  
  ! The Q1 element is specified by four polynomials on the reference element.
  ! These four polynomials are:
  !
  !  P1(X,Y) = 1/4 (1-x) (1-y)
  !  P2(X,Y) = 1/4 (1+x) (1-y)
  !  P3(X,Y) = 1/4 (1+x) (1+y)
  !  P4(X,Y) = 1/4 (1-x) (1+y)
  !
  ! Each of them calculated that way that Pi(Xj)=delta_ij (Kronecker)
  ! for X1,X2,X3,X4 the four corners of the reference element [-1,1]x[-1,1].
  
  ! Clear the output array
  !Dbas = 0.0_DP
  
  dx = Dpoint(1)
  dy = Dpoint(2)
    
  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  ! If function values are desired, calculate them.
!  if (el_bder(DER_FUNC)) then
    Dbas(1,DER_FUNC) = 0.25E0_DP*(1E0_DP-dx)*(1E0_DP-dy)
    Dbas(2,DER_FUNC) = 0.25E0_DP*(1E0_DP+dx)*(1E0_DP-dy)
    Dbas(3,DER_FUNC) = 0.25E0_DP*(1E0_DP+dx)*(1E0_DP+dy)
    Dbas(4,DER_FUNC) = 0.25E0_DP*(1E0_DP-dx)*(1E0_DP+dy)
!  endif
  
  ! If x-or y-derivatives are desired, calculate them.
  ! The values of the derivatives are calculated by taking the
  ! derivative of the polynomials and multiplying them with the
  ! inverse of the transformation matrix (in each point) as
  ! stated above.
!  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then
    dxj = 0.25E0_DP / ddetj
    
    ! x- and y-derivatives on reference element
    Dhelp(1,1) =-(1E0_DP-dy)
    Dhelp(2,1) = (1E0_DP-dy)
    Dhelp(3,1) = (1E0_DP+dy)
    Dhelp(4,1) =-(1E0_DP+dy)
    Dhelp(1,2) =-(1E0_DP-dx)
    Dhelp(2,2) =-(1E0_DP+dx)
    Dhelp(3,2) = (1E0_DP+dx)
    Dhelp(4,2) = (1E0_DP-dx)
      
    ! x-derivatives on current element
!    if (Bder(DER_DERIV_X)) then
      Dbas(1,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(1,1) - Djac(2) * Dhelp(1,2))
      Dbas(2,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(2,1) - Djac(2) * Dhelp(2,2))
      Dbas(3,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(3,1) - Djac(2) * Dhelp(3,2))
      Dbas(4,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(4,1) - Djac(2) * Dhelp(4,2))
!    endif
    
    ! y-derivatives on current element
!    if (Bder(DER_DERIV_Y)) then
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
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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
  ! DIMENSION(#space dimensions,npoints).
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
  REAL(DP), DIMENSION(4,NDIM2D,npoints) :: Dhelp

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
  INTEGER(I32), INTENT(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements).
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dcoords
  
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
  ! DIMENSION(#space dimensions,npoints,nelements).
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
  REAL(DP), DIMENSION(4,NDIM2D,npoints) :: Dhelp

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
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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
  ! DIMENSION(#space dimensions,npoints).
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
  INTEGER(I32), INTENT(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dcoords
  
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
  ! DIMENSION(#space dimensions,npoints,nelements).
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
! Element subroutines for parametric QP1 element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  PURE SUBROUTINE elem_QP1 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                            Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_QP1.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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

  ! local variables
  ! REAL(DP), DIMENSION(2,2) :: Dhelp
  REAL(DP) :: dxj

    ! The element is given the function value and the X- and Y-derivative
    ! in the midpoint of the reference element.
    ! That means: p(x,y) = a*1 + b*x + c*y, coefficients (a,b,c).
    !
    ! So the value of the first basis function is 1 on the whole
    ! element.
    ! The value of the basis function "x" is obviously x
    ! and the value of the basis function "y" is y.
    
    Dbas(1,DER_FUNC) = 1
    Dbas(2,DER_FUNC) = Dpoint(1)
    Dbas(3,DER_FUNC) = Dpoint(2)
    
    ! 1st X- and Y-derivative on the reference element are obvious...
    
    ! Dhelp(1,1) = 1.0_DP   ! dx/dx = 1
    ! Dhelp(2,1) = 0.0_DP   ! dy/dx = 0

    ! Dhelp(1,2) = 0.0_DP   ! dx/dy = 0
    ! Dhelp(2,2) = 1.0_DP   ! dy/dy = 1
  
    ! Use them to calculate the derivative in the cubature point
    ! on the real element. 
  
    dxj = 0.25_DP/ddetj
    
    ! X-derivatives on current element
    ! Dbas(1,DER_DERIV_X) = 0
    ! Dbas(2,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(1,1) - Djac(2) * Dhelp(1,2))
    ! Dbas(3,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(2,1) - Djac(2) * Dhelp(2,2))
    Dbas(1,DER_DERIV_X) = 0.0_DP
    Dbas(2,DER_DERIV_X) = dxj * Djac(4)
    Dbas(3,DER_DERIV_X) = -dxj * Djac(2) 
    
    ! y-derivatives on current element
    ! Dbas(1,DER_DERIV_Y) = 0
    ! Dbas(2,DER_DERIV_Y) = dxj * (-Djac(3) * Dhelp(1,1) + Djac(1) * Dhelp(1,2))
    ! Dbas(3,DER_DERIV_Y) = dxj * (-Djac(3) * Dhelp(2,1) + Djac(1) * Dhelp(2,2))
    Dbas(1,DER_DERIV_Y) = 0.0_DP
    Dbas(2,DER_DERIV_Y) = -dxj * Djac(3)
    Dbas(3,DER_DERIV_Y) = dxj * Djac(1) 
  
  END SUBROUTINE 

!**************************************************************************
! Element subroutines for Q1~ element, integral mean value based.
!**************************************************************************
  
! The standard integral-based Q1~-element looks locally as follows:
!
!                 phi_3
!           +-----X-----+                            +-----e3----+
!           |           |                            |           |
!           |           |                            |           |
!     phi_4 X           X phi_2    for the edges     e4          e2
!           |           |                            |           |
!           |           |                            |           |
!           +-----X-----+                            +-----e1----+
!                 phi_1
! 
! on the reference element [-1,1] x [-1,1].
!
! On the element, we can see the four basis functions phi_1, ..., phi_4.
! Correspondiong to these, we have the following four local basis functions:
!
!  p_1(x,y) = a1 (x^2-y^2)  +  b1 x  +  c1 y  +  d1
!  p_2(x,y) = a2 (x^2-y^2)  +  b2 x  +  c2 y  +  d2
!  p_3(x,y) = a3 (x^2-y^2)  +  b3 x  +  c3 y  +  d3
!  p_4(x,y) = a4 (x^2-y^2)  +  b4 x  +  c4 y  +  d4
!
! each of them designed (with ai, bi, ci and di) in such a way, such that
!
!      1/|ei| int_ei p_j(x,y) ds = delta_ij
!
! Solving this 4x4 system given by this integral set gives for the standard 
! parametric integral mean based Q1~ element the following local polynomials:
!
!  p_1(x,y) = -3/8 (x^2-y^2)  -  1/2 y  +  1/4
!  p_2(x,y) =  3/8 (x^2-y^2)  +  1/2 x  +  1/4
!  p_3(x,y) = -3/8 (x^2-y^2)  +  1/2 y  +  1/4
!  p_4(x,y) =  3/8 (x^2-y^2)  -  1/2 x  +  1/4
!
! with x=-1..1, y=-1..1.
!
! The nonparametric variant of the integral-mean based Q1~ extends this.
! We have the usual mapping between the reference element and the real element:
!
!   +-----e3----+                        +--------E3--------+
!   |           |                       /                    \
!   |           |       sigma          /                      \
!   e4          e2     ------->       E4                       \
!   |           |                    /                          E2
!   |           |                   /                            \ 
!   +-----e1----+                  /                              \
!                                 +----_____                       \
!                                           -----__E1_              \
!                                                      -----_____    \
!                                                                -----+
!
! 
! with the bilinear mapping
!
!    sigma: [-1,1]^2 -> T
!
!    sigma(x,y) = s1 x^2  +  s2 y^2  +  s3 x y  +  s4 x  +  s5 y  + s6
!
! On the real element T, we fix the midpoints of the edges on E1, ..., E4 and
! define a normalised local coordinate system in terms of xi and eta by taking
! the opposite midpoints as vectors of the coordinate system. 
!                                           
!           +---------X--------+            
!          /          ^         \           
!         /           |vec_2     \          
!        X--------____|___        \               
!       /             |   -------->X    
!      /              |     vec_1   \       
!     /               |              \      
!    +----_____       |               \     
!              -----__X__              \    
!                         -----_____    \   
!                                   -----+
!
! We shift both of these vectors vec_1 and vec_2 into the origin (0,0)
! and normalise them to have length=1 (otherwise, the LBB condition might
! get violated!). So the picture we are looking at here is the following:
!
!   ^ xi             +---------X--------+          
!   |               /          ^         \         
!   |              /           |vec_2     \        
!   |             X--------____|___        \       
!   |            /             |   -------->X    
!   |           /              |     vec_1   \     
!   |          /               |              \    
!   |         +----_____       |               \   
!   |                   -----__X__              \  
!   |                              -----_____    \ 
!   |                                        -----+
!   |
!   |
!   X--------________          
! (0,0)              --------________            
!                                    ---> eta
!
! Every point (x,y) in the 'old' standard coordinate system has a
! representation (z1,z2) in the new coordinate system. To get this
! representation, we need a linear mapping. We define it as:
!
!   ( z1 ) := r(x,y) := ( k11 k12 ) ( x ) 
!   ( z2 )              ( k21 k22 ) ( y )
!
! This mapping should fulfill:
!
!   ( eta_1 ) = ( k11 k12 ) ( 1 ) 
!   ( eta_2 )   ( k21 k22 ) ( 0 )
!  
!   ( xi_1 ) = ( k11 k12 ) ( 0 ) 
!   ( xi_2 )   ( k21 k22 ) ( 1 )
!
! so that the vector eta has the coordinates (1,0) and xi the coordinates
! (0,1). This simply means:
!
!   r(x,y) = ( eta_1 xi_1 ) ( x )  =  ( eta_1 x  +  xi_1 y ) 
!            ( eta_2 xi_2 ) ( y )     ( eta_2 x  +  xi_2 y )
!
! Then, we set up the local basis functions in terms of xi and eta,
! i.e. in the coordinate space of the new coordinate system, 
! with a new set of (unknown) coefficents ai, bi, ci and di:
!
!  P1(z1,z2)  :=  d1 (z1^2 - z2^2)  +  c1 z2  +  b1 z1  +  a1
!  P2(z1,z2)  :=  d2 (z1^2 - z2^2)  +  c2 z2  +  b2 z1  +  a2
!  P3(z1,z2)  :=  d3 (z1^2 - z2^2)  +  c3 z2  +  b3 z1  +  a3
!  P4(z1,z2)  :=  d4 (z1^2 - z2^2)  +  c4 z2  +  b4 z1  +  a4
!
! Later, we want to evaluate these P_i. Each P_i consists of a linear 
! combination of some coefficients ai, bi, ci, di with the four monoms
!
!  m1(xi,eta) := 1
!  m2(xi,eta) := z1
!  m3(xi,eta) := z2
!  m4(xi,eta) := (z1^2-z2^2)
!
! To evaluate these mi's in the new coordinate system, we concatenate them
! with the mapping r(.,.). As result, we get the four functions F1,F2,F3,F4
! which are defined as functions at the bottom of this routine:
!
!  F1(x,y) := m1(r(x,y)) = 1
!  F2(x,y) := m2(r(x,y)) = eta_1 x  +  xi_1 y
!  F3(x,y) := m3(r(x,y)) = eta_2 x  +  xi_2 y
!  F4(x,y) := m4(r(x,y)) = ( eta_1 x + xi_1 y )^2 - ( eta_2 x + xi_2 y )^2
!                        =           ( eta_1^2 - eta_2^2 ) x^2 
!                          + 2 ( eta_1 xi_1 - eta_2 xi_2 ) x y
!                          +           ( xi_1^2 - xi_2^2 ) y^2
!
! So the polynomials have now the form:
!
!  P1(r(x,y)) = a1 F1(x,y)  +  b1 F2(x,y)  +  c1 F3(x,y)  +  d1 F4(x,y)
!  P2(r(x,y)) = a2 F1(x,y)  +  b2 F2(x,y)  +  c2 F3(x,y)  +  d2 F4(x,y)
!  P3(r(x,y)) = a3 F1(x,y)  +  b3 F2(x,y)  +  c3 F3(x,y)  +  d3 F4(x,y)
!  P4(r(x,y)) = a4 F1(x,y)  +  b4 F2(x,y)  +  c4 F3(x,y)  +  d4 F4(x,y)
!
! It doesn't matter whether the local coordinate system starts in (0,0) or in
! the midpoint of the element or whereever. As the rotation of the element
! coincides with the rotation of the new coordinate system, the polynomial 
! space is unisolvent and therefore exist the above local basis functions uniquely.
!
! The coefficients ai, bi, ci, di are to be calculated in such a way, that
!
!    1/|ei| int_ei Pj(r(.)) ds = delta_ij
!
! holds. The integral "1/|ei| int_ei ... ds" over the edge Ei is approximated
! by a 2-point gauss integral.
! This gives four 4x4 systems for the computation of ai, bi, ci and di, which
! can be written as:
!
!   ( . . . . ) ( a1 a2 a3 a4 ) = ( 1 0 0 0 )
!   ( . .V. . ) ( b1 b2 b3 b4 )   ( 0 1 0 0 )
!   ( . . . . ) ( c1 c2 c3 c4 )   ( 0 0 1 0 )
!   ( . . . . ) ( d1 d2 d3 d4 )   ( 0 0 0 1 )
!
! So to get all the coefficients, one has to calculate V^-1 !
! The entries of the matrix V = {v_ij} are defined (because of the linearity
! of the integral) as
!
!         vij = 1/|ei| int_ei Fj(x,y) ds
!
! Now let's go...
 
!**************************************************************************
! Element subroutines for nonparametric Q1~ element, integral mean value
! based.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
!**************************************************************************

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
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE),
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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
  REAL(DP),DIMENSION(4) :: DXM,DYM,DLX,DLY
  REAL(DP),DIMENSION(4,4) :: A       ! local 4x4 system
  REAL(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC3
  REAL(DP),DIMENSION(EL_MAXNBAS,EL_MAXNCOF) :: COB
  REAL(DP) :: dx, dy
  
  ! Clear the output array
  !Dbas = 0.0_DP
  
  ! Where to evaluate? On the real element T...
  
  dx = Dpoint(1)
  dy = Dpoint(2)
  
  ! Calculate the edge midpoints and length of edges: 
  !  DXM(:) := X-coordinates of midpoints
  !  DYM(:) := Y-coordinates of midpoints
  !  DLX(:) := length of each edge in X-direction
  !  DLY(:) := length of each edge in Y-direction
  ! So SQRT(DLX(:)**2+DLY(:)**2) = length of the edges.
  
  DO IVE=1,NVE
    DXM(IVE)=0.5_DP*(Dcoords(1,IVE)+Dcoords(1,MOD(IVE,4)+1))
    DYM(IVE)=0.5_DP*(Dcoords(2,IVE)+Dcoords(2,MOD(IVE,4)+1))
    DLX(IVE)=0.5_DP*(Dcoords(1,MOD(IVE,4)+1)-Dcoords(1,IVE))
    DLY(IVE)=0.5_DP*(Dcoords(2,MOD(IVE,4)+1)-Dcoords(2,IVE))
  END DO

  ! Calculate the scaling factors for the local coordinate system.
  !  D1 := 1 / ||xi||_2
  !  D2 := 1 / ||eta||_2
  
  D1 = 1.0_DP / SQRT((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
  D2 = 1.0_DP / SQRT((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)
  
  ! Calculate the vector eta = (CA1,CB1); these numbers coincide
  ! with the coefficients of the polynomial F2(x,y) := m2(r(x,y))
  CA1 = (DXM(2)-DXM(4)) * D1
  CB1 = (DYM(2)-DYM(4)) * D1
  
  ! Calculate the vector xi = (CA2,CB2); these numbers coincide
  ! with the coefficients of the polynomial F3(x,y) := m3(r(x,y))
  CA2 = (DXM(3)-DXM(1)) * D2
  CB2 = (DYM(3)-DYM(1)) * D2
  
  ! Calculate the coefficients of the polynomial F4(x,y) := m4(r(x,y))
  CA3 = CA1**2-CA2**2
  CB3 = 2.0_DP*(CA1*CB1-CA2*CB2)
  CC3 = CB1**2-CB2**2
  
  ! Calculate the matrix V (=A) with vij = int_ei Fj (x,y) d(x,y).
  ! Loop over the edges.
  DO IA = 1,NVE
  
    ! Calculate the X- and Y-coordinates of the two Gauss points on the
    ! current edge. (PXL,PYL) = 1st Gauss point, (PXU,PYU) = 2nd Gauss point.
    
    PXL = DXM(IA)-SQRT(1.0_DP/3.0_DP)*DLX(IA)
    PYL = DYM(IA)-SQRT(1.0_DP/3.0_DP)*DLY(IA)
    PXU = DXM(IA)+SQRT(1.0_DP/3.0_DP)*DLX(IA)
    PYU = DYM(IA)+SQRT(1.0_DP/3.0_DP)*DLY(IA)
    
    ! Set up the coefficients of the linear system to calculate ai, bi, ci and di;
    ! i.e. "1/|ei| int_ei Fj(x,y) ds"
    A(1,IA)=0.5_DP*( F1(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                    +F1(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
    A(2,IA)=0.5_DP*( F2(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                    +F2(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
    A(3,IA)=0.5_DP*( F3(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                    +F3(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
    A(4,IA)=0.5_DP*( F4(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                    +F4(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
  END DO
  
  ! Invert that matrix V to get the matrix of the coefficients of the
  ! four polynomials. The matix A (=V) is replaced by the inverse.
  !CALL INVERT(A,F,CKH,0)
  CALL mprim_invertMatrixPivotDble(A,4)

  ! Ok, the coefficients ai, bi, ci, di are calculated.
  ! The next point is: We want to evaluate the polynoms Pi(r(.))
  ! in the point (x,y) which is specified as input parameter to this routine!
  !
  ! For this purpose, we first transform the polynom Pi(r(.)) into the
  ! monomial representation:
  !
  !  Pi(r(x,y)) = ai F4(x,y)  +  bi F3(x,y)  +  ci F2(x,y)  +  di F1(x,y)
  !             =   COB(i,1) x^2  +  COB(i,2) y^2  +  COB(i,3) x y
  !               + COB(i,4) x    +  COB(i,5) y
  !               + COB(i,6)

  DO IK=1,4
    COB(IK,1) = A(IK,4)*CA3
    COB(IK,2) = A(IK,4)*CC3
    COB(IK,3) = A(IK,4)*CB3
    COB(IK,4) = A(IK,2)*CA1+A(IK,3)*CA2
    COB(IK,5) = A(IK,2)*CB1+A(IK,3)*CB2
    COB(IK,6) = A(IK,1)
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
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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
  ! DIMENSION(#space dimensions,npoints)
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
  REAL(DP),DIMENSION(4) :: DXM,DYM,DLX,DLY
  REAL(DP),DIMENSION(4,4) :: A       ! local 4x4 system
  REAL(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC3
  REAL(DP),DIMENSION(EL_MAXNBAS,EL_MAXNCOF) :: COB
  
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
    A(1,IA)=0.5_DP*( F1(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                  +F1(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
    A(2,IA)=0.5_DP*( F2(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                  +F2(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
    A(3,IA)=0.5_DP*( F3(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                  +F3(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
    A(4,IA)=0.5_DP*( F4(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                  +F4(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
  END DO
  
  !CALL INVERT(A,F,CKH,0)
  CALL mprim_invertMatrixPivotDble(A,4)

  DO IK=1,4
    COB(IK,1) = A(IK,4)*CA3
    COB(IK,2) = A(IK,4)*CC3
    COB(IK,3) = A(IK,4)*CB3
    COB(IK,4) = A(IK,2)*CA1+A(IK,3)*CA2
    COB(IK,5) = A(IK,2)*CB1+A(IK,3)*CB2
    COB(IK,6) = A(IK,1)
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
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dcoords
  
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
  ! DIMENSION(#space dimensions,npoints,nelements)
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
  REAL(DP),DIMENSION(4) :: DXM,DYM,DLX,DLY
  REAL(DP), DIMENSION(4,4) :: A       ! local 4x4 system
  REAL(DP), DIMENSION(4,4) :: CK
  REAL(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC3
  REAL(DP) :: COB(EL_MAXNBAS,EL_MAXNCOF)
  INTEGER :: npointsfunc,npointsderx,npointsdery
  
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
  
  ! Check which element variant we have; with or without pivoting...
  IF (IAND(ieltyp,INT(2**17,I32)) .EQ. 0) THEN

    ! Check whether to scale the local coordinate system or not.
  
    IF (IAND(ieltyp,INT(2**18,I32)) .EQ. 0) THEN
  
      ! Use pivoting and scaled local coordinate system for 
      ! increased numerical stability.
    
      ! Loop over the elements
      
      DO j=1,nelements
        
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
          
          A(1,IA)=0.5_DP*( F1(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F1(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(2,IA)=0.5_DP*( F2(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F2(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(3,IA)=0.5_DP*( F3(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F3(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(4,IA)=0.5_DP*( F4(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F4(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
        END DO
        
        ! Invert the matrix in-place.
        !CALL INVERT(A,F,CKH,0)
        CALL mprim_invertMatrixPivotDble(A,4)

        ! In comparison to the standard EM30 routine above, we use the
        ! COB-array transposed!:
        !
        !  Pi(r(x,y)) = ai F4(x,y)  +  bi F3(x,y)  +  ci F2(x,y)  +  di F1(x,y)
        !             =   COB(1,i) x^2  +  COB(2,i) y^2  +  COB(3,i) x y
        !               + COB(4,i) x    +  COB(5,i) y
        !               + COB(6,i)
        !
        ! This gives easier array access to the processor and gains a little
        ! bit speed!

        DO IK=1,4
          COB(1,IK) = A(IK,4)*CA3
          COB(2,IK) = A(IK,4)*CC3
          COB(3,IK) = A(IK,4)*CB3
          COB(4,IK) = A(IK,2)*CA1+A(IK,3)*CA2
          COB(5,IK) = A(IK,2)*CB1+A(IK,3)*CB2
          COB(6,IK) = A(IK,1)
        END DO

        ! Function values
        DO i=1,npointsfunc   ! either 0 or npoints

          Dbas(1,DER_FUNC,i,j)=  COB(1,1)*Dpoints(1,i,j)**2 &
                                +COB(2,1)*Dpoints(2,i,j)**2 &
                                +COB(3,1)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                                +COB(4,1)*Dpoints(1,i,j)   &
                                +COB(5,1)*Dpoints(2,i,j)   &
                                +COB(6,1)
          Dbas(2,DER_FUNC,i,j)=  COB(1,2)*Dpoints(1,i,j)**2 &
                                +COB(2,2)*Dpoints(2,i,j)**2 &
                                +COB(3,2)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                                +COB(4,2)*Dpoints(1,i,j)   &
                                +COB(5,2)*Dpoints(2,i,j)   &
                                +COB(6,2)
          Dbas(3,DER_FUNC,i,j)=  COB(1,3)*Dpoints(1,i,j)**2 &
                                +COB(2,3)*Dpoints(2,i,j)**2 &
                                +COB(3,3)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                                +COB(4,3)*Dpoints(1,i,j)   &
                                +COB(5,3)*Dpoints(2,i,j)   &
                                +COB(6,3)
          Dbas(4,DER_FUNC,i,j)=  COB(1,4)*Dpoints(1,i,j)**2 &
                                +COB(2,4)*Dpoints(2,i,j)**2 &
                                +COB(3,4)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                                +COB(4,4)*Dpoints(1,i,j)    &
                                +COB(5,4)*Dpoints(2,i,j)    &
                                +COB(6,4)
        END DO
        
        ! x-derivatives
              
        DO i=1,npointsderx   ! either 0 or npoints
          Dbas(1,DER_DERIV_X,i,j) = 2.0_DP*COB(1,1)*Dpoints(1,i,j) &
                                          +COB(3,1)*Dpoints(2,i,j) &
                                          +COB(4,1)
          Dbas(2,DER_DERIV_X,i,j) = 2.0_DP*COB(1,2)*Dpoints(1,i,j) &
                                          +COB(3,2)*Dpoints(2,i,j) &
                                          +COB(4,2)
          Dbas(3,DER_DERIV_X,i,j) = 2.0_DP*COB(1,3)*Dpoints(1,i,j) &
                                          +COB(3,3)*Dpoints(2,i,j) &
                                          +COB(4,3)
          Dbas(4,DER_DERIV_X,i,j) = 2.0_DP*COB(1,4)*Dpoints(1,i,j) &
                                          +COB(3,4)*Dpoints(2,i,j) &
                                          +COB(4,4)
        END DO

        ! y-derivatives
              
        DO i=1,npointsdery   ! either 0 or npoints
          Dbas(1,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,1)*Dpoints(2,i,j) &
                                          +COB(3,1)*Dpoints(1,i,j) &
                                          +COB(5,1)
          Dbas(2,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,2)*Dpoints(2,i,j) &
                                          +COB(3,2)*Dpoints(1,i,j) &
                                          +COB(5,2)
          Dbas(3,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,3)*Dpoints(2,i,j) &
                                          +COB(3,3)*Dpoints(1,i,j) &
                                          +COB(5,3)
          Dbas(4,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,4)*Dpoints(2,i,j) &
                                          +COB(3,4)*Dpoints(1,i,j) &
                                          +COB(5,4)

        END DO
        
      END DO ! ielement
      
    ELSE
    
      ! Use pivoting for increased numerical stability.
      ! Don't scaled local coordinate system 
    
      ! Loop over the elements
      
      DO j=1,nelements
        
        DO IVE=1,NVE
          DXM(IVE)=0.5_DP*(Dcoords(1,IVE,j)+Dcoords(1,MOD(IVE,4)+1,j))
          DYM(IVE)=0.5_DP*(Dcoords(2,IVE,j)+Dcoords(2,MOD(IVE,4)+1,j))
          DLX(IVE)=0.5_DP*(Dcoords(1,MOD(IVE,4)+1,j)-Dcoords(1,IVE,j))
          DLY(IVE)=0.5_DP*(Dcoords(2,MOD(IVE,4)+1,j)-Dcoords(2,IVE,j))
        END DO

        ! Don't scale the local coordinate system in this approach.
        ! Saves a little it numerical effort.

        CA1 = (DXM(2)-DXM(4))
        CB1 = (DYM(2)-DYM(4))
        CA2 = (DXM(3)-DXM(1))
        CB2 = (DYM(3)-DYM(1))
        CA3 = CA1**2-CA2**2
        CB3 = 2.0_DP*(CA1*CB1-CA2*CB2)
        CC3 = CB1**2-CB2**2
        
        DO IA = 1,4
          PXL = DXM(IA)-SQRT(1.0_DP/3.0_DP)*DLX(IA)
          PYL = DYM(IA)-SQRT(1.0_DP/3.0_DP)*DLY(IA)
          PXU = DXM(IA)+SQRT(1.0_DP/3.0_DP)*DLX(IA)
          PYU = DYM(IA)+SQRT(1.0_DP/3.0_DP)*DLY(IA)
          
          A(1,IA)=0.5_DP*( F1(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F1(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(2,IA)=0.5_DP*( F2(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F2(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(3,IA)=0.5_DP*( F3(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F3(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(4,IA)=0.5_DP*( F4(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F4(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
        END DO
        
        ! Invert the matrix in-place.
        !CALL INVERT(A,F,CKH,0)
        CALL mprim_invertMatrixPivotDble(A,4)

        ! In comparison to the standard EM30 routine above, we use the
        ! COB-array transposed!:
        !
        !  Pi(r(x,y)) = ai F4(x,y)  +  bi F3(x,y)  +  ci F2(x,y)  +  di F1(x,y)
        !             =   COB(1,i) x^2  +  COB(2,i) y^2  +  COB(3,i) x y
        !               + COB(4,i) x    +  COB(5,i) y
        !               + COB(6,i)
        !
        ! This gives easier array access to the processor and gains a little
        ! bit speed!

        DO IK=1,4
          COB(1,IK) = A(IK,4)*CA3
          COB(2,IK) = A(IK,4)*CC3
          COB(3,IK) = A(IK,4)*CB3
          COB(4,IK) = A(IK,2)*CA1+A(IK,3)*CA2
          COB(5,IK) = A(IK,2)*CB1+A(IK,3)*CB2
          COB(6,IK) = A(IK,1)
        END DO

        ! Function values
        DO i=1,npointsfunc   ! either 0 or npoints

          Dbas(1,DER_FUNC,i,j)=  COB(1,1)*Dpoints(1,i,j)**2 &
                                +COB(2,1)*Dpoints(2,i,j)**2 &
                                +COB(3,1)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                                +COB(4,1)*Dpoints(1,i,j)   &
                                +COB(5,1)*Dpoints(2,i,j)   &
                                +COB(6,1)
          Dbas(2,DER_FUNC,i,j)=  COB(1,2)*Dpoints(1,i,j)**2 &
                                +COB(2,2)*Dpoints(2,i,j)**2 &
                                +COB(3,2)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                                +COB(4,2)*Dpoints(1,i,j)   &
                                +COB(5,2)*Dpoints(2,i,j)   &
                                +COB(6,2)
          Dbas(3,DER_FUNC,i,j)=  COB(1,3)*Dpoints(1,i,j)**2 &
                                +COB(2,3)*Dpoints(2,i,j)**2 &
                                +COB(3,3)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                                +COB(4,3)*Dpoints(1,i,j)   &
                                +COB(5,3)*Dpoints(2,i,j)   &
                                +COB(6,3)
          Dbas(4,DER_FUNC,i,j)=  COB(1,4)*Dpoints(1,i,j)**2 &
                                +COB(2,4)*Dpoints(2,i,j)**2 &
                                +COB(3,4)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                                +COB(4,4)*Dpoints(1,i,j)    &
                                +COB(5,4)*Dpoints(2,i,j)    &
                                +COB(6,4)
        END DO
        
        ! x-derivatives
              
        DO i=1,npointsderx   ! either 0 or npoints
          Dbas(1,DER_DERIV_X,i,j) = 2.0_DP*COB(1,1)*Dpoints(1,i,j) &
                                          +COB(3,1)*Dpoints(2,i,j) &
                                          +COB(4,1)
          Dbas(2,DER_DERIV_X,i,j) = 2.0_DP*COB(1,2)*Dpoints(1,i,j) &
                                          +COB(3,2)*Dpoints(2,i,j) &
                                          +COB(4,2)
          Dbas(3,DER_DERIV_X,i,j) = 2.0_DP*COB(1,3)*Dpoints(1,i,j) &
                                          +COB(3,3)*Dpoints(2,i,j) &
                                          +COB(4,3)
          Dbas(4,DER_DERIV_X,i,j) = 2.0_DP*COB(1,4)*Dpoints(1,i,j) &
                                          +COB(3,4)*Dpoints(2,i,j) &
                                          +COB(4,4)
        END DO

        ! y-derivatives
              
        DO i=1,npointsdery   ! either 0 or npoints
          Dbas(1,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,1)*Dpoints(2,i,j) &
                                          +COB(3,1)*Dpoints(1,i,j) &
                                          +COB(5,1)
          Dbas(2,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,2)*Dpoints(2,i,j) &
                                          +COB(3,2)*Dpoints(1,i,j) &
                                          +COB(5,2)
          Dbas(3,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,3)*Dpoints(2,i,j) &
                                          +COB(3,3)*Dpoints(1,i,j) &
                                          +COB(5,3)
          Dbas(4,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,4)*Dpoints(2,i,j) &
                                          +COB(3,4)*Dpoints(1,i,j) &
                                          +COB(5,4)

        END DO
        
      END DO ! ielement

    END IF
    
  ELSE
  
    ! Don't use pivoting.
    
    ! Check whether to scae the local coordinate system or not.
  
    IF (IAND(ieltyp,INT(2**18,I32)) .EQ. 0) THEN
  
      ! Loop over the elements
      
      DO j=1,nelements
        
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
          A(1,IA)=0.5_DP*( F1(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F1(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(2,IA)=0.5_DP*( F2(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F2(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(3,IA)=0.5_DP*( F3(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F3(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(4,IA)=0.5_DP*( F4(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F4(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
        END DO
        
        ! Invert the matrix to get the coefficients.
        ! Use direct inversion and save the result to CK directly.
        CALL mprim_invert4x4MatrixDirectDble(A,CK)

        ! In comparison to the standard EM30 routine above, we use the
        ! COB-array transposed!:
        !
        !  Pi(r(x,y)) = ai F4(x,y)  +  bi F3(x,y)  +  ci F2(x,y)  +  di F1(x,y)
        !             =   COB(1,i) x^2  +  COB(2,i) y^2  +  COB(3,i) x y
        !               + COB(4,i) x    +  COB(5,i) y
        !               + COB(6,i)
        !
        ! This gives easier array access to the processor and gains a little
        ! bit speed!

        ! Calculate the coefficients of the monoms.
        DO IK=1,4
          COB(1,IK) = CK(IK,4)*CA3
          COB(2,IK) = CK(IK,4)*CC3
          COB(3,IK) = CK(IK,4)*CB3
          COB(4,IK) = CK(IK,2)*CA1+CK(IK,3)*CA2
          COB(5,IK) = CK(IK,2)*CB1+CK(IK,3)*CB2
          COB(6,IK) = CK(IK,1)
        END DO

        ! Function values
        DO i=1,npointsfunc   ! either 0 or npoints

          Dbas(1,DER_FUNC,i,j)=  COB(1,1)*Dpoints(1,i,j)**2 &
                                +COB(2,1)*Dpoints(2,i,j)**2 &
                                +COB(3,1)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                                +COB(4,1)*Dpoints(1,i,j)   &
                                +COB(5,1)*Dpoints(2,i,j)   &
                                +COB(6,1)
          Dbas(2,DER_FUNC,i,j)=  COB(1,2)*Dpoints(1,i,j)**2 &
                                +COB(2,2)*Dpoints(2,i,j)**2 &
                                +COB(3,2)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                                +COB(4,2)*Dpoints(1,i,j)   &
                                +COB(5,2)*Dpoints(2,i,j)   &
                                +COB(6,2)
          Dbas(3,DER_FUNC,i,j)=  COB(1,3)*Dpoints(1,i,j)**2 &
                                +COB(2,3)*Dpoints(2,i,j)**2 &
                                +COB(3,3)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                                +COB(4,3)*Dpoints(1,i,j)   &
                                +COB(5,3)*Dpoints(2,i,j)   &
                                +COB(6,3)
          Dbas(4,DER_FUNC,i,j)=  COB(1,4)*Dpoints(1,i,j)**2 &
                                +COB(2,4)*Dpoints(2,i,j)**2 &
                                +COB(3,4)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                                +COB(4,4)*Dpoints(1,i,j)    &
                                +COB(5,4)*Dpoints(2,i,j)    &
                                +COB(6,4)
        END DO
        
        ! x-derivatives
              
        DO i=1,npointsderx   ! either 0 or npoints
          Dbas(1,DER_DERIV_X,i,j) = 2.0_DP*COB(1,1)*Dpoints(1,i,j) &
                                          +COB(3,1)*Dpoints(2,i,j) &
                                          +COB(4,1)
          Dbas(2,DER_DERIV_X,i,j) = 2.0_DP*COB(1,2)*Dpoints(1,i,j) &
                                          +COB(3,2)*Dpoints(2,i,j) &
                                          +COB(4,2)
          Dbas(3,DER_DERIV_X,i,j) = 2.0_DP*COB(1,3)*Dpoints(1,i,j) &
                                          +COB(3,3)*Dpoints(2,i,j) &
                                          +COB(4,3)
          Dbas(4,DER_DERIV_X,i,j) = 2.0_DP*COB(1,4)*Dpoints(1,i,j) &
                                          +COB(3,4)*Dpoints(2,i,j) &
                                          +COB(4,4)
        END DO

        ! y-derivatives
              
        DO i=1,npointsdery   ! either 0 or npoints
          Dbas(1,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,1)*Dpoints(2,i,j) &
                                          +COB(3,1)*Dpoints(1,i,j) &
                                          +COB(5,1)
          Dbas(2,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,2)*Dpoints(2,i,j) &
                                          +COB(3,2)*Dpoints(1,i,j) &
                                          +COB(5,2)
          Dbas(3,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,3)*Dpoints(2,i,j) &
                                          +COB(3,3)*Dpoints(1,i,j) &
                                          +COB(5,3)
          Dbas(4,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,4)*Dpoints(2,i,j) &
                                          +COB(3,4)*Dpoints(1,i,j) &
                                          +COB(5,4)

        END DO
        
      END DO ! ielement
    
    ELSE
    
      ! Loop over the elements
      
      DO j=1,nelements
        
        DO IVE=1,NVE
          DXM(IVE)=0.5_DP*(Dcoords(1,IVE,j)+Dcoords(1,MOD(IVE,4)+1,j))
          DYM(IVE)=0.5_DP*(Dcoords(2,IVE,j)+Dcoords(2,MOD(IVE,4)+1,j))
          DLX(IVE)=0.5_DP*(Dcoords(1,MOD(IVE,4)+1,j)-Dcoords(1,IVE,j))
          DLY(IVE)=0.5_DP*(Dcoords(2,MOD(IVE,4)+1,j)-Dcoords(2,IVE,j))
        END DO

        ! Don't scale the local coordinate in this approach.
        ! Saves a little bit numerical effort.

        CA1 = (DXM(2)-DXM(4))
        CB1 = (DYM(2)-DYM(4))
        CA2 = (DXM(3)-DXM(1))
        CB2 = (DYM(3)-DYM(1))
        CA3 = CA1**2-CA2**2
        CB3 = 2.0_DP*(CA1*CB1-CA2*CB2)
        CC3 = CB1**2-CB2**2
        
        DO IA = 1,4
          PXL = DXM(IA)-SQRT(1.0_DP/3.0_DP)*DLX(IA)
          PYL = DYM(IA)-SQRT(1.0_DP/3.0_DP)*DLY(IA)
          PXU = DXM(IA)+SQRT(1.0_DP/3.0_DP)*DLX(IA)
          PYU = DYM(IA)+SQRT(1.0_DP/3.0_DP)*DLY(IA)
          A(1,IA)=0.5_DP*( F1(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F1(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(2,IA)=0.5_DP*( F2(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F2(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(3,IA)=0.5_DP*( F3(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F3(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
          A(4,IA)=0.5_DP*( F4(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3) &
                          +F4(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
        END DO
        
        ! Invert the matrix to get the coefficients.
        ! Use direct inversion and save the result to CK directly.
        CALL mprim_invert4x4MatrixDirectDble(A,CK)

        ! In comparison to the standard EM30 routine above, we use the
        ! COB-array transposed!:
        !
        !  Pi(r(x,y)) = ai F4(x,y)  +  bi F3(x,y)  +  ci F2(x,y)  +  di F1(x,y)
        !             =   COB(1,i) x^2  +  COB(2,i) y^2  +  COB(3,i) x y
        !               + COB(4,i) x    +  COB(5,i) y
        !               + COB(6,i)
        !
        ! This gives easier array access to the processor and gains a little
        ! bit speed!

        ! Calculate the coefficients of the monoms.
        DO IK=1,4
          COB(1,IK) = CK(IK,4)*CA3
          COB(2,IK) = CK(IK,4)*CC3
          COB(3,IK) = CK(IK,4)*CB3
          COB(4,IK) = CK(IK,2)*CA1+CK(IK,3)*CA2
          COB(5,IK) = CK(IK,2)*CB1+CK(IK,3)*CB2
          COB(6,IK) = CK(IK,1)
        END DO

        ! Function values
        DO i=1,npointsfunc   ! either 0 or npoints

          Dbas(1,DER_FUNC,i,j)=  COB(1,1)*Dpoints(1,i,j)**2 &
                                +COB(2,1)*Dpoints(2,i,j)**2 &
                                +COB(3,1)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                                +COB(4,1)*Dpoints(1,i,j)   &
                                +COB(5,1)*Dpoints(2,i,j)   &
                                +COB(6,1)
          Dbas(2,DER_FUNC,i,j)=  COB(1,2)*Dpoints(1,i,j)**2 &
                                +COB(2,2)*Dpoints(2,i,j)**2 &
                                +COB(3,2)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                                +COB(4,2)*Dpoints(1,i,j)   &
                                +COB(5,2)*Dpoints(2,i,j)   &
                                +COB(6,2)
          Dbas(3,DER_FUNC,i,j)=  COB(1,3)*Dpoints(1,i,j)**2 &
                                +COB(2,3)*Dpoints(2,i,j)**2 &
                                +COB(3,3)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                                +COB(4,3)*Dpoints(1,i,j)   &
                                +COB(5,3)*Dpoints(2,i,j)   &
                                +COB(6,3)
          Dbas(4,DER_FUNC,i,j)=  COB(1,4)*Dpoints(1,i,j)**2 &
                                +COB(2,4)*Dpoints(2,i,j)**2 &
                                +COB(3,4)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                                +COB(4,4)*Dpoints(1,i,j)    &
                                +COB(5,4)*Dpoints(2,i,j)    &
                                +COB(6,4)
        END DO
        
        ! x-derivatives
              
        DO i=1,npointsderx   ! either 0 or npoints
          Dbas(1,DER_DERIV_X,i,j) = 2.0_DP*COB(1,1)*Dpoints(1,i,j) &
                                          +COB(3,1)*Dpoints(2,i,j) &
                                          +COB(4,1)
          Dbas(2,DER_DERIV_X,i,j) = 2.0_DP*COB(1,2)*Dpoints(1,i,j) &
                                          +COB(3,2)*Dpoints(2,i,j) &
                                          +COB(4,2)
          Dbas(3,DER_DERIV_X,i,j) = 2.0_DP*COB(1,3)*Dpoints(1,i,j) &
                                          +COB(3,3)*Dpoints(2,i,j) &
                                          +COB(4,3)
          Dbas(4,DER_DERIV_X,i,j) = 2.0_DP*COB(1,4)*Dpoints(1,i,j) &
                                          +COB(3,4)*Dpoints(2,i,j) &
                                          +COB(4,4)
        END DO

        ! y-derivatives
              
        DO i=1,npointsdery   ! either 0 or npoints
          Dbas(1,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,1)*Dpoints(2,i,j) &
                                          +COB(3,1)*Dpoints(1,i,j) &
                                          +COB(5,1)
          Dbas(2,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,2)*Dpoints(2,i,j) &
                                          +COB(3,2)*Dpoints(1,i,j) &
                                          +COB(5,2)
          Dbas(3,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,3)*Dpoints(2,i,j) &
                                          +COB(3,3)*Dpoints(1,i,j) &
                                          +COB(5,3)
          Dbas(4,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,4)*Dpoints(2,i,j) &
                                          +COB(3,4)*Dpoints(1,i,j) &
                                          +COB(5,4)

        END DO
        
      END DO ! ielement

    END IF
  
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

!**************************************************************************
! Element subroutines for parametric Q1~ element, integral mean value
! based.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
!**************************************************************************
 
!<subroutine>  

  PURE SUBROUTINE elem_E030 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                             Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q1.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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

  ! auxiliary vector containing the first derivatives on the reference element
  REAL(DP), DIMENSION(4,NDIM2D) :: Dhelp
  REAL(DP) :: dx,dy

  REAL(DP) :: dxj !auxiliary variable
  
  ! The Q1 element is specified by four polynomials on the reference element.
  ! These four polynomials are:
  !
  !  p_1(x,y) = -3/8 (x^2-y^2)  -  1/2 y  +  1/4
  !  p_2(x,y) =  3/8 (x^2-y^2)  +  1/2 x  +  1/4
  !  p_3(x,y) = -3/8 (x^2-y^2)  +  1/2 y  +  1/4
  !  p_4(x,y) =  3/8 (x^2-y^2)  -  1/2 x  +  1/4
  !
  ! Each of them calculated that way that Pi(Xj)=delta_ij (Kronecker)
  ! for X1,X2,X3,X4 the ingegral over the four edges of the reference 
  ! element [-1,1]x[-1,1].
  
  ! Clear the output array
  !Dbas = 0.0_DP
  
  dx = Dpoint(1)
  dy = Dpoint(2)
    
  ! Remark: The Q1~-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  ! If function values are desired, calculate them.
!  if (el_bder(DER_FUNC)) then
    Dbas(1,DER_FUNC) = 0.125E0_DP*(-3.0_DP*(dx**2-dy**2)-4.0_DP*dy+2.0_DP) 
    Dbas(2,DER_FUNC) = 0.125E0_DP*( 3.0_DP*(dx**2-dy**2)+4.0_DP*dx+2.0_DP)
    Dbas(3,DER_FUNC) = 0.125E0_DP*( 3.0_DP*(dx**2-dy**2)-4.0_DP*dy-2.0_DP)
    Dbas(4,DER_FUNC) = 0.125E0_DP*(-3.0_DP*(dx**2-dy**2)+4.0_DP*dx-2.0_DP)
!  endif
  
  ! If x-or y-derivatives are desired, calculate them.
  ! The values of the derivatives are calculated by taking the
  ! derivative of the polynomials and multiplying them with the
  ! inverse of the transformation matrix (in each point) as
  ! stated above.
!  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then
    dxj = 0.125E0_DP / ddetj
    
    ! x- and y-derivatives on reference element
    Dhelp(1,1) = -6.0_DP*dx
    Dhelp(2,1) =  6.0_DP*dx+4.0_DP
    Dhelp(3,1) = -6.0_DP*dx
    Dhelp(4,1) =  6.0_DP*dx-4.0_DP
    Dhelp(1,2) =  6.0_DP*dy-4.0_DP
    Dhelp(2,2) = -6.0_DP*dy
    Dhelp(3,2) =  6.0_DP*dy+4.0_DP
    Dhelp(4,2) = -6.0_DP*dy
      
    ! x-derivatives on current element
!    if (Bder(DER_DERIV_X)) then
      Dbas(1,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(1,1) - Djac(2) * Dhelp(1,2))
      Dbas(2,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(2,1) - Djac(2) * Dhelp(2,2))
      Dbas(3,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(3,1) - Djac(2) * Dhelp(3,2))
      Dbas(4,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(4,1) - Djac(2) * Dhelp(4,2))
!    endif
    
    ! y-derivatives on current element
!    if (Bder(DER_DERIV_Y)) then
      Dbas(1,DER_DERIV_Y) = -dxj * (Djac(3) * Dhelp(1,1) - Djac(1) * Dhelp(1,2))
      Dbas(2,DER_DERIV_Y) = -dxj * (Djac(3) * Dhelp(2,1) - Djac(1) * Dhelp(2,2))
      Dbas(3,DER_DERIV_Y) = -dxj * (Djac(3) * Dhelp(3,1) - Djac(1) * Dhelp(3,2))
      Dbas(4,DER_DERIV_Y) = -dxj * (Djac(3) * Dhelp(4,1) - Djac(1) * Dhelp(4,2))
!    endif
!  endif
    
  END SUBROUTINE 
  
  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_E030_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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
  ! DIMENSION(#space dimensions,npoints).
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
  REAL(DP), DIMENSION(4,NDIM2D,npoints) :: Dhelp

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
      Dbas(1,DER_FUNC,i) = 0.125E0_DP* &
          (-3.0_DP*(Dpoints(1,i)**2-Dpoints(2,i)**2)-4.0_DP*Dpoints(2,i)+2.0_DP) 
      Dbas(2,DER_FUNC,i) = 0.125E0_DP* &
          ( 3.0_DP*(Dpoints(1,i)**2-Dpoints(2,i)**2)+4.0_DP*Dpoints(1,i)+2.0_DP)
      Dbas(3,DER_FUNC,i) = 0.125E0_DP* &
          ( 3.0_DP*(Dpoints(1,i)**2-Dpoints(2,i)**2)-4.0_DP*Dpoints(2,i)-2.0_DP)
      Dbas(4,DER_FUNC,i) = 0.125E0_DP* &
          (-3.0_DP*(Dpoints(1,i)**2-Dpoints(2,i)**2)+4.0_DP*Dpoints(1,i)-2.0_DP)
    END DO
  !ENDIF
  
  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN
    dxj = 0.125E0_DP / Ddetj
    
    !x- and y-derivatives on reference element
    DO i=1,npoints
      Dhelp(1,1,i) = -6.0_DP*Dpoints(1,i)
      Dhelp(2,1,i) =  6.0_DP*Dpoints(1,i)+4.0_DP
      Dhelp(3,1,i) = -6.0_DP*Dpoints(1,i)
      Dhelp(4,1,i) =  6.0_DP*Dpoints(1,i)-4.0_DP
      Dhelp(1,2,i) =  6.0_DP*Dpoints(2,i)-4.0_DP
      Dhelp(2,2,i) = -6.0_DP*Dpoints(2,i)
      Dhelp(3,2,i) =  6.0_DP*Dpoints(2,i)+4.0_DP
      Dhelp(4,2,i) = -6.0_DP*Dpoints(2,i)
    END DO
      
    !x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      DO i=1,npoints
        Dbas(1,DER_DERIV_X,i) = &
            dxj(i) * (Djac(4,i) * Dhelp(1,1,i) - Djac(2,i) * Dhelp(1,2,i))
        Dbas(2,DER_DERIV_X,i) = &
            dxj(i) * (Djac(4,i) * Dhelp(2,1,i) - Djac(2,i) * Dhelp(2,2,i))
        Dbas(3,DER_DERIV_X,i) = &
            dxj(i) * (Djac(4,i) * Dhelp(3,1,i) - Djac(2,i) * Dhelp(3,2,i))
        Dbas(4,DER_DERIV_X,i) = &
            dxj(i) * (Djac(4,i) * Dhelp(4,1,i) - Djac(2,i) * Dhelp(4,2,i))
!      END DO
!    ENDIF
    
    !y-derivatives on current element
!    IF (Bder(DER_DERIV_Y)) THEN
!      DO i=1,npoints
        Dbas(1,DER_DERIV_Y,i) = &
            -dxj(i) * (Djac(3,i) * Dhelp(1,1,i) - Djac(1,i) * Dhelp(1,2,i))
        Dbas(2,DER_DERIV_Y,i) = &
            -dxj(i) * (Djac(3,i) * Dhelp(2,1,i) - Djac(1,i) * Dhelp(2,2,i))
        Dbas(3,DER_DERIV_Y,i) = &
            -dxj(i) * (Djac(3,i) * Dhelp(3,1,i) - Djac(1,i) * Dhelp(3,2,i))
        Dbas(4,DER_DERIV_Y,i) = &
            -dxj(i) * (Djac(3,i) * Dhelp(4,1,i) - Djac(1,i) * Dhelp(4,2,i))
      END DO
!    ENDIF
!  ENDIF
    
  END SUBROUTINE 

  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_E030_sim (ieltyp, Dcoords, Djac, Ddetj, &
                                 Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1.
  INTEGER(I32), INTENT(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dcoords
  
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
  ! DIMENSION(#space dimensions,npoints,nelements).
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
  REAL(DP), DIMENSION(4,NDIM2D,npoints) :: Dhelp

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
        Dbas(1,DER_FUNC,i,j) = 0.125E0_DP* &
            (-3.0_DP*(Dpoints(1,i,j)**2-Dpoints(2,i,j)**2)-4.0_DP*Dpoints(2,i,j)+2.0_DP) 
        Dbas(2,DER_FUNC,i,j) = 0.125E0_DP* &
            ( 3.0_DP*(Dpoints(1,i,j)**2-Dpoints(2,i,j)**2)+4.0_DP*Dpoints(1,i,j)+2.0_DP)
        Dbas(3,DER_FUNC,i,j) = 0.125E0_DP* &
            ( 3.0_DP*(Dpoints(1,i,j)**2-Dpoints(2,i,j)**2)-4.0_DP*Dpoints(2,i,j)-2.0_DP)
        Dbas(4,DER_FUNC,i,j) = 0.125E0_DP* &
            (-3.0_DP*(Dpoints(1,i,j)**2-Dpoints(2,i,j)**2)+4.0_DP*Dpoints(1,i,j)-2.0_DP)
      END DO
      
    END DO
    
  END IF
    
  !if x-or y-derivatives are desired
  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN
  
    DO j=1,nelements
      dxj = 0.125E0_DP / Ddetj(:,j)
      
      !x- and y-derivatives on reference element
      DO i=1,npoints
        Dhelp(1,1,i) = -6.0_DP*Dpoints(1,i,j)
        Dhelp(2,1,i) =  6.0_DP*Dpoints(1,i,j)+4.0_DP
        Dhelp(3,1,i) = -6.0_DP*Dpoints(1,i,j)
        Dhelp(4,1,i) =  6.0_DP*Dpoints(1,i,j)-4.0_DP
        Dhelp(1,2,i) =  6.0_DP*Dpoints(2,i,j)-4.0_DP
        Dhelp(2,2,i) = -6.0_DP*Dpoints(2,i,j)
        Dhelp(3,2,i) =  6.0_DP*Dpoints(2,i,j)+4.0_DP
        Dhelp(4,2,i) = -6.0_DP*Dpoints(2,i,j)
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
! Element subroutines for Q1~ element, midpoint value based.
!**************************************************************************
  
! The standard integral-based Q1~-element looks locally as follows:
!
!                 phi_3
!           +-----X-----+                            +-----m3----+
!           |           |                            |           |
!           |           |                            |           |
!     phi_4 X           X phi_2  for the midpoints   m4          m2
!           |           |                            |           |
!           |           |                            |           |
!           +-----X-----+                            +-----m1----+
!                 phi_1
! 
! on the reference element [-1,1] x [-1,1].
!
! On the element, we can see the four basis functions phi_1, ..., phi_4.
! Correspondiong to these, we have the following four local basis functions:
!
!  p_1(x,y) = a1 (x^2-y^2)  +  b1 x  +  c1 y  +  d1
!  p_2(x,y) = a2 (x^2-y^2)  +  b2 x  +  c2 y  +  d2
!  p_3(x,y) = a3 (x^2-y^2)  +  b3 x  +  c3 y  +  d3
!  p_4(x,y) = a4 (x^2-y^2)  +  b4 x  +  c4 y  +  d4
!
! each of them designed (with ai, bi, ci and di) in such a way, such that
!
!      p_i(mj) = delta_ij
!
! Solving this 4x4 system given by this integral set gives for the standard 
! parametric integral mean based Q1~ element the following local polynomials:
!
!  p_1(x,y) = -1/4 (x^2-y^2)  -  1/2 y  +  1/4
!  p_2(x,y) =  1/4 (x^2-y^2)  +  1/2 x  +  1/4
!  p_3(x,y) = -1/4 (x^2-y^2)  +  1/2 y  +  1/4
!  p_4(x,y) =  1/4 (x^2-y^2)  -  1/2 x  +  1/4
!
! with x=-1..1, y=-1..1.
!
! The extension in the nonparametric case is done as above for the
! integral ean value based element.

!**************************************************************************
! Element subroutines for nonparametric Q1~ element, midpoint value based.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
!**************************************************************************

!<subroutine>  

  PURE SUBROUTINE elem_EM31 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                             Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point. The coordinates are expected
  ! on the real element!
!</description>

  !<input>

  ! Element type identifier. Must be =EL_EM30.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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
  REAL(DP) :: D1,D2
  REAL(DP),DIMENSION(4) :: DXM,DYM,DLX,DLY
  REAL(DP),DIMENSION(4,4) :: A       ! local 4x4 system
  REAL(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC3
  REAL(DP),DIMENSION(EL_MAXNBAS,EL_MAXNCOF) :: COB
  REAL(DP) :: dx, dy
  
  ! Clear the output array
  !Dbas = 0.0_DP
  
  ! Where to evaluate? On the real element T...
  
  dx = Dpoint(1)
  dy = Dpoint(2)
  
  ! Calculate the edge midpoints and length of edges: 
  !  DXM(:) := X-coordinates of midpoints
  !  DYM(:) := Y-coordinates of midpoints
  !  DLX(:) := length of each edge in X-direction
  !  DLY(:) := length of each edge in Y-direction
  ! So SQRT(DLX(:)**2+DLY(:)**2) = length of the edges.
  
  DO IVE=1,NVE
    DXM(IVE)=0.5_DP*(Dcoords(1,IVE)+Dcoords(1,MOD(IVE,4)+1))
    DYM(IVE)=0.5_DP*(Dcoords(2,IVE)+Dcoords(2,MOD(IVE,4)+1))
    DLX(IVE)=0.5_DP*(Dcoords(1,MOD(IVE,4)+1)-Dcoords(1,IVE))
    DLY(IVE)=0.5_DP*(Dcoords(2,MOD(IVE,4)+1)-Dcoords(2,IVE))
  END DO

  ! Calculate the scaling factors for the local coordinate system.
  !  D1 := 1 / ||xi||_2
  !  D2 := 1 / ||eta||_2
  
  D1 = 1.0_DP / SQRT((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
  D2 = 1.0_DP / SQRT((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)
  
  ! Calculate the vector eta = (CA1,CB1); these numbers coincide
  ! with the coefficients of the polynomial F2(x,y) := m2(r(x,y))
  CA1 = (DXM(2)-DXM(4)) * D1
  CB1 = (DYM(2)-DYM(4)) * D1
  
  ! Calculate the vector xi = (CA2,CB2); these numbers coincide
  ! with the coefficients of the polynomial F3(x,y) := m3(r(x,y))
  CA2 = (DXM(3)-DXM(1)) * D2
  CB2 = (DYM(3)-DYM(1)) * D2
  
  ! Calculate the coefficients of the polynomial F4(x,y) := m4(r(x,y))
  CA3 = CA1**2-CA2**2
  CB3 = 2.0_DP*(CA1*CB1-CA2*CB2)
  CC3 = CB1**2-CB2**2
  
  ! Calculate the matrix V (=A) with vij = int_ei Fj (x,y) d(x,y).
  ! Loop over the edges.
  DO IA = 1,NVE
    ! Set up the coefficients of the linear system to calculate ai, bi, ci and di.
    ! Use the X- and Y-coordinates of the midpointof every edge to evaluate
    ! the Fi.
    A(1,IA)=F1(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    A(2,IA)=F2(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    A(3,IA)=F3(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    A(4,IA)=F4(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
  END DO
  
  ! Invert that matrix V to get the matrix of the coefficients of the
  ! four polynomials. The matix A (=V) is replaced by the inverse.
  !CALL INVERT(A,F,CKH,0)
  CALL mprim_invertMatrixPivotDble(A,4)

  ! Ok, the coefficients ai, bi, ci, di are calculated.
  ! The next point is: We want to evaluate the polynoms Pi(r(.))
  ! in the point (x,y) which is specified as input parameter to this routine!
  !
  ! For this purpose, we first transform the polynom Pi(r(.)) into the
  ! monomial representation:
  !
  !  Pi(r(x,y)) = ai F4(x,y)  +  bi F3(x,y)  +  ci F2(x,y)  +  di F1(x,y)
  !             =   COB(i,1) x^2  +  COB(i,2) y^2  +  COB(i,3) x y
  !               + COB(i,4) x    +  COB(i,5) y
  !               + COB(i,6)

  DO IK=1,4
    COB(IK,1) = A(IK,4)*CA3
    COB(IK,2) = A(IK,4)*CC3
    COB(IK,3) = A(IK,4)*CB3
    COB(IK,4) = A(IK,2)*CA1+A(IK,3)*CA2
    COB(IK,5) = A(IK,2)*CB1+A(IK,3)*CB2
    COB(IK,6) = A(IK,1)
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

  PURE SUBROUTINE elem_EM31_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given points. The coordinates are expected
  ! on the real element!
!</description>

!<input>
  ! Element type identifier. Must be =EL_EM30.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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
  ! DIMENSION(#space dimensions,npoints)
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
  REAL(DP) :: D1,D2
  REAL(DP),DIMENSION(4) :: DXM,DYM,DLX,DLY
  REAL(DP),DIMENSION(4,4) :: A       ! local 4x4 system
  REAL(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC3
  REAL(DP),DIMENSION(EL_MAXNBAS,EL_MAXNCOF) :: COB
  
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
    A(1,IA)=F1(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    A(2,IA)=F2(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    A(3,IA)=F3(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
    A(4,IA)=F4(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
  END DO
  
  !CALL INVERT(A,F,CKH,0)
  CALL mprim_invertMatrixPivotDble(A,4)

  DO IK=1,4
    COB(IK,1) = A(IK,4)*CA3
    COB(IK,2) = A(IK,4)*CC3
    COB(IK,3) = A(IK,4)*CB3
    COB(IK,4) = A(IK,2)*CA1+A(IK,3)*CA2
    COB(IK,5) = A(IK,2)*CB1+A(IK,3)*CB2
    COB(IK,6) = A(IK,1)
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

  SUBROUTINE elem_EM31_sim (ieltyp, Dcoords, Djac, Ddetj, &
                            Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_EM30.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dcoords
  
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
  ! DIMENSION(#space dimensions,npoints,nelements)
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
  REAL(DP) :: D1,D2
  REAL(DP),DIMENSION(4) :: DXM,DYM,DLX,DLY
  REAL(DP), DIMENSION(4,4) :: A       ! local 4x4 system
  REAL(DP), DIMENSION(4,4) :: CK
  REAL(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC3
  REAL(DP) :: COB(EL_MAXNBAS,EL_MAXNCOF)
  INTEGER :: npointsfunc,npointsderx,npointsdery
  
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
  
  ! Check which element variant we have; with or without pivoting...
  IF (IAND(ieltyp,INT(2**17,I32)) .EQ. 0) THEN
  
    ! Use pivoting for increased numerical stability.
  
    ! Loop over the elements
    
    DO j=1,nelements
      
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
        A(1,IA)=F1(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
        A(2,IA)=F2(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
        A(3,IA)=F3(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
        A(4,IA)=F4(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
      END DO
      
      ! Invert the matrix in-place.
      !CALL INVERT(A,F,CKH,0)
      CALL mprim_invertMatrixPivotDble(A,4)

      ! In comparison to the standard EM30 routine above, we use the
      ! COB-array transposed!:
      !
      !  Pi(r(x,y)) = ai F4(x,y)  +  bi F3(x,y)  +  ci F2(x,y)  +  di F1(x,y)
      !             =   COB(1,i) x^2  +  COB(2,i) y^2  +  COB(3,i) x y
      !               + COB(4,i) x    +  COB(5,i) y
      !               + COB(6,i)
      !
      ! This gives easier array access to the processor and gains a little
      ! bit speed!

      DO IK=1,4
        COB(1,IK) = A(IK,4)*CA3
        COB(2,IK) = A(IK,4)*CC3
        COB(3,IK) = A(IK,4)*CB3
        COB(4,IK) = A(IK,2)*CA1+A(IK,3)*CA2
        COB(5,IK) = A(IK,2)*CB1+A(IK,3)*CB2
        COB(6,IK) = A(IK,1)
      END DO

      ! Function values
      DO i=1,npointsfunc   ! either 0 or npoints

        Dbas(1,DER_FUNC,i,j)=  COB(1,1)*Dpoints(1,i,j)**2 &
                              +COB(2,1)*Dpoints(2,i,j)**2 &
                              +COB(3,1)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                              +COB(4,1)*Dpoints(1,i,j)   &
                              +COB(5,1)*Dpoints(2,i,j)   &
                              +COB(6,1)
        Dbas(2,DER_FUNC,i,j)=  COB(1,2)*Dpoints(1,i,j)**2 &
                              +COB(2,2)*Dpoints(2,i,j)**2 &
                              +COB(3,2)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                              +COB(4,2)*Dpoints(1,i,j)   &
                              +COB(5,2)*Dpoints(2,i,j)   &
                              +COB(6,2)
        Dbas(3,DER_FUNC,i,j)=  COB(1,3)*Dpoints(1,i,j)**2 &
                              +COB(2,3)*Dpoints(2,i,j)**2 &
                              +COB(3,3)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                              +COB(4,3)*Dpoints(1,i,j)   &
                              +COB(5,3)*Dpoints(2,i,j)   &
                              +COB(6,3)
        Dbas(4,DER_FUNC,i,j)=  COB(1,4)*Dpoints(1,i,j)**2 &
                              +COB(2,4)*Dpoints(2,i,j)**2 &
                              +COB(3,4)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                              +COB(4,4)*Dpoints(1,i,j)    &
                              +COB(5,4)*Dpoints(2,i,j)    &
                              +COB(6,4)
      END DO
      
      ! x-derivatives
            
      DO i=1,npointsderx   ! either 0 or npoints
        Dbas(1,DER_DERIV_X,i,j) = 2.0_DP*COB(1,1)*Dpoints(1,i,j) &
                                        +COB(3,1)*Dpoints(2,i,j) &
                                        +COB(4,1)
        Dbas(2,DER_DERIV_X,i,j) = 2.0_DP*COB(1,2)*Dpoints(1,i,j) &
                                        +COB(3,2)*Dpoints(2,i,j) &
                                        +COB(4,2)
        Dbas(3,DER_DERIV_X,i,j) = 2.0_DP*COB(1,3)*Dpoints(1,i,j) &
                                        +COB(3,3)*Dpoints(2,i,j) &
                                        +COB(4,3)
        Dbas(4,DER_DERIV_X,i,j) = 2.0_DP*COB(1,4)*Dpoints(1,i,j) &
                                        +COB(3,4)*Dpoints(2,i,j) &
                                        +COB(4,4)
      END DO

      ! y-derivatives
            
      DO i=1,npointsdery   ! either 0 or npoints
        Dbas(1,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,1)*Dpoints(2,i,j) &
                                        +COB(3,1)*Dpoints(1,i,j) &
                                        +COB(5,1)
        Dbas(2,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,2)*Dpoints(2,i,j) &
                                        +COB(3,2)*Dpoints(1,i,j) &
                                        +COB(5,2)
        Dbas(3,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,3)*Dpoints(2,i,j) &
                                        +COB(3,3)*Dpoints(1,i,j) &
                                        +COB(5,3)
        Dbas(4,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,4)*Dpoints(2,i,j) &
                                        +COB(3,4)*Dpoints(1,i,j) &
                                        +COB(5,4)

      END DO
      
    END DO ! ielement
    
  ELSE
  
    ! Don't use pivoting.
    
    ! Loop over the elements
    
    DO j=1,nelements
      
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
        A(1,IA)=F1(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
        A(2,IA)=F2(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
        A(3,IA)=F3(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
        A(4,IA)=F4(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
      END DO
      
      ! Invert the matrix to get the coefficients.
      ! Use direct inversion and save the result to CK directly.
      CALL mprim_invert4x4MatrixDirectDble(A,CK)

      ! In comparison to the standard EM31 routine above, we use the
      ! COB-array transposed!:
      !
      !  Pi(r(x,y)) = ai F4(x,y)  +  bi F3(x,y)  +  ci F2(x,y)  +  di F1(x,y)
      !             =   COB(1,i) x^2  +  COB(2,i) y^2  +  COB(3,i) x y
      !               + COB(4,i) x    +  COB(5,i) y
      !               + COB(6,i)
      !
      ! This gives easier array access to the processor and gains a little
      ! bit speed!

      ! Calculate the coefficients of the monoms.
      DO IK=1,4
        COB(1,IK) = CK(IK,4)*CA3
        COB(2,IK) = CK(IK,4)*CC3
        COB(3,IK) = CK(IK,4)*CB3
        COB(4,IK) = CK(IK,2)*CA1+CK(IK,3)*CA2
        COB(5,IK) = CK(IK,2)*CB1+CK(IK,3)*CB2
        COB(6,IK) = CK(IK,1)
      END DO

      ! Function values
      DO i=1,npointsfunc   ! either 0 or npoints

        Dbas(1,DER_FUNC,i,j)=  COB(1,1)*Dpoints(1,i,j)**2 &
                              +COB(2,1)*Dpoints(2,i,j)**2 &
                              +COB(3,1)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                              +COB(4,1)*Dpoints(1,i,j)   &
                              +COB(5,1)*Dpoints(2,i,j)   &
                              +COB(6,1)
        Dbas(2,DER_FUNC,i,j)=  COB(1,2)*Dpoints(1,i,j)**2 &
                              +COB(2,2)*Dpoints(2,i,j)**2 &
                              +COB(3,2)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                              +COB(4,2)*Dpoints(1,i,j)   &
                              +COB(5,2)*Dpoints(2,i,j)   &
                              +COB(6,2)
        Dbas(3,DER_FUNC,i,j)=  COB(1,3)*Dpoints(1,i,j)**2 &
                              +COB(2,3)*Dpoints(2,i,j)**2 &
                              +COB(3,3)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                              +COB(4,3)*Dpoints(1,i,j)   &
                              +COB(5,3)*Dpoints(2,i,j)   &
                              +COB(6,3)
        Dbas(4,DER_FUNC,i,j)=  COB(1,4)*Dpoints(1,i,j)**2 &
                              +COB(2,4)*Dpoints(2,i,j)**2 &
                              +COB(3,4)*Dpoints(1,i,j)*Dpoints(2,i,j) &
                              +COB(4,4)*Dpoints(1,i,j)    &
                              +COB(5,4)*Dpoints(2,i,j)    &
                              +COB(6,4)
      END DO
      
      ! x-derivatives
            
      DO i=1,npointsderx   ! either 0 or npoints
        Dbas(1,DER_DERIV_X,i,j) = 2.0_DP*COB(1,1)*Dpoints(1,i,j) &
                                        +COB(3,1)*Dpoints(2,i,j) &
                                        +COB(4,1)
        Dbas(2,DER_DERIV_X,i,j) = 2.0_DP*COB(1,2)*Dpoints(1,i,j) &
                                        +COB(3,2)*Dpoints(2,i,j) &
                                        +COB(4,2)
        Dbas(3,DER_DERIV_X,i,j) = 2.0_DP*COB(1,3)*Dpoints(1,i,j) &
                                        +COB(3,3)*Dpoints(2,i,j) &
                                        +COB(4,3)
        Dbas(4,DER_DERIV_X,i,j) = 2.0_DP*COB(1,4)*Dpoints(1,i,j) &
                                        +COB(3,4)*Dpoints(2,i,j) &
                                        +COB(4,4)
      END DO

      ! y-derivatives
            
      DO i=1,npointsdery   ! either 0 or npoints
        Dbas(1,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,1)*Dpoints(2,i,j) &
                                        +COB(3,1)*Dpoints(1,i,j) &
                                        +COB(5,1)
        Dbas(2,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,2)*Dpoints(2,i,j) &
                                        +COB(3,2)*Dpoints(1,i,j) &
                                        +COB(5,2)
        Dbas(3,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,3)*Dpoints(2,i,j) &
                                        +COB(3,3)*Dpoints(1,i,j) &
                                        +COB(5,3)
        Dbas(4,DER_DERIV_Y,i,j) = 2.0_DP*COB(2,4)*Dpoints(2,i,j) &
                                        +COB(3,4)*Dpoints(1,i,j) &
                                        +COB(5,4)

      END DO
      
    END DO ! ielement
  
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

!**************************************************************************
! Element subroutines for parametric Q1~ element, midpoint value based.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
!**************************************************************************
 
!<subroutine>  

  PURE SUBROUTINE elem_E031 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                             Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q1.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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
  REAL(DP), DIMENSION(4,NDIM2D) :: Dhelp
  REAL(DP) :: dx,dy

  REAL(DP) :: dxj !auxiliary variable
  
  ! The Q1 element is specified by four polynomials on the reference element.
  ! These four polynomials are:
  !
  !  p_1(x,y) = -1/4 (x^2-y^2)  -  1/2 y  +  1/4
  !  p_2(x,y) =  1/4 (x^2-y^2)  +  1/2 x  +  1/4
  !  p_3(x,y) = -1/4 (x^2-y^2)  +  1/2 y  +  1/4
  !  p_4(x,y) =  1/4 (x^2-y^2)  -  1/2 x  +  1/4
  !
  ! Each of them calculated that way that Pi(Xj)=delta_ij (Kronecker)
  ! for X1,X2,X3,X4 the four midpoints of the reference element [-1,1]x[-1,1].
  
  ! Clear the output array
  !Dbas = 0.0_DP
  
  dx = Dpoint(1)
  dy = Dpoint(2)
    
  ! Remark: The Q1~-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  ! If function values are desired, calculate them.
!  if (el_bder(DER_FUNC)) then
    Dbas(1,DER_FUNC) = 0.25E0_DP*(-(dx**2-dy**2)-2.0_DP*dy+1.0_DP) 
    Dbas(2,DER_FUNC) = 0.25E0_DP*( (dx**2-dy**2)+2.0_DP*dx+1.0_DP)
    Dbas(3,DER_FUNC) = 0.25E0_DP*(-(dx**2-dy**2)+2.0_DP*dy+1.0_DP)
    Dbas(4,DER_FUNC) = 0.25E0_DP*( (dx**2-dy**2)-2.0_DP*dx+1.0_DP)
!  endif
  
  ! If x-or y-derivatives are desired, calculate them.
  ! The values of the derivatives are calculated by taking the
  ! derivative of the polynomials and multiplying them with the
  ! inverse of the transformation matrix (in each point) as
  ! stated above.
!  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then
    dxj = 0.5E0_DP / ddetj
    
    ! x- and y-derivatives on reference element
    Dhelp(1,1) = -dx
    Dhelp(2,1) =  dx+1.0_DP
    Dhelp(3,1) = -dx
    Dhelp(4,1) =  dx-1.0_DP
    Dhelp(1,2) =  dy-1.0_DP
    Dhelp(2,2) = -dy
    Dhelp(3,2) =  dy+1.0_DP
    Dhelp(4,2) = -dy
      
    ! x-derivatives on current element
!    if (Bder(DER_DERIV_X)) then
      Dbas(1,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(1,1) - Djac(2) * Dhelp(1,2))
      Dbas(2,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(2,1) - Djac(2) * Dhelp(2,2))
      Dbas(3,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(3,1) - Djac(2) * Dhelp(3,2))
      Dbas(4,DER_DERIV_X) = dxj * (Djac(4) * Dhelp(4,1) - Djac(2) * Dhelp(4,2))
!    endif
    
    ! y-derivatives on current element
!    if (Bder(DER_DERIV_Y)) then
      Dbas(1,DER_DERIV_Y) = -dxj * (Djac(3) * Dhelp(1,1) - Djac(1) * Dhelp(1,2))
      Dbas(2,DER_DERIV_Y) = -dxj * (Djac(3) * Dhelp(2,1) - Djac(1) * Dhelp(2,2))
      Dbas(3,DER_DERIV_Y) = -dxj * (Djac(3) * Dhelp(3,1) - Djac(1) * Dhelp(3,2))
      Dbas(4,DER_DERIV_Y) = -dxj * (Djac(3) * Dhelp(4,1) - Djac(1) * Dhelp(4,2))
!    endif
!  endif
    
  END SUBROUTINE 
  
  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_E031_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
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
  ! DIMENSION(#space dimensions,npoints).
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
  REAL(DP), DIMENSION(4,NDIM2D,npoints) :: Dhelp

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
      Dbas(1,DER_FUNC,i) = 0.25E0_DP*&
          (-(Dpoints(1,i)**2-Dpoints(2,i)**2)-2.0_DP*Dpoints(2,i)+1.0_DP) 
      Dbas(2,DER_FUNC,i) = 0.25E0_DP*&
          ( (Dpoints(1,i)**2-Dpoints(2,i)**2)+2.0_DP*Dpoints(1,i)+1.0_DP)
      Dbas(3,DER_FUNC,i) = 0.25E0_DP*&
          (-(Dpoints(1,i)**2-Dpoints(2,i)**2)+2.0_DP*Dpoints(2,i)+1.0_DP)
      Dbas(4,DER_FUNC,i) = 0.25E0_DP*&
          ( (Dpoints(1,i)**2-Dpoints(2,i)**2)-2.0_DP*Dpoints(1,i)+1.0_DP)
    END DO
  !ENDIF
  
  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN
    dxj = 0.5E0_DP / Ddetj
    
    !x- and y-derivatives on reference element
    DO i=1,npoints
      Dhelp(1,1,i) = -Dpoints(1,i)
      Dhelp(2,1,i) =  Dpoints(1,i)+1.0_DP
      Dhelp(3,1,i) = -Dpoints(1,i)
      Dhelp(4,1,i) =  Dpoints(1,i)-1.0_DP
      Dhelp(1,2,i) =  Dpoints(2,i)-1.0_DP
      Dhelp(2,2,i) = -Dpoints(2,i)
      Dhelp(3,2,i) =  Dpoints(2,i)+1.0_DP
      Dhelp(4,2,i) = -Dpoints(2,i)
    END DO
      
    !x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      DO i=1,npoints
        Dbas(1,DER_DERIV_X,i) = &
            dxj(i) * (Djac(4,i) * Dhelp(1,1,i) - Djac(2,i) * Dhelp(1,2,i))
        Dbas(2,DER_DERIV_X,i) = & 
            dxj(i) * (Djac(4,i) * Dhelp(2,1,i) - Djac(2,i) * Dhelp(2,2,i))
        Dbas(3,DER_DERIV_X,i) = &
            dxj(i) * (Djac(4,i) * Dhelp(3,1,i) - Djac(2,i) * Dhelp(3,2,i))
        Dbas(4,DER_DERIV_X,i) = &
            dxj(i) * (Djac(4,i) * Dhelp(4,1,i) - Djac(2,i) * Dhelp(4,2,i))
!      END DO
!    ENDIF
    
    !y-derivatives on current element
!    IF (Bder(DER_DERIV_Y)) THEN
!      DO i=1,npoints
        Dbas(1,DER_DERIV_Y,i) = &
            -dxj(i) * (Djac(3,i) * Dhelp(1,1,i) - Djac(1,i) * Dhelp(1,2,i))
        Dbas(2,DER_DERIV_Y,i) = &
            -dxj(i) * (Djac(3,i) * Dhelp(2,1,i) - Djac(1,i) * Dhelp(2,2,i))
        Dbas(3,DER_DERIV_Y,i) = &
            -dxj(i) * (Djac(3,i) * Dhelp(3,1,i) - Djac(1,i) * Dhelp(3,2,i))
        Dbas(4,DER_DERIV_Y,i) = &
            -dxj(i) * (Djac(3,i) * Dhelp(4,1,i) - Djac(1,i) * Dhelp(4,2,i))
      END DO
!    ENDIF
!  ENDIF
    
  END SUBROUTINE 

  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_E031_sim (ieltyp, Dcoords, Djac, Ddetj, &
                                 Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1.
  INTEGER(I32), INTENT(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dcoords
  
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
  ! DIMENSION(#space dimensions,npoints,nelements).
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
  REAL(DP), DIMENSION(4,NDIM2D,npoints) :: Dhelp

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
        Dbas(1,DER_FUNC,i,j) = 0.25E0_DP*&
            (-(Dpoints(1,i,j)**2-Dpoints(2,i,j)**2)-2.0_DP*Dpoints(2,i,j)+1.0_DP) 
        Dbas(2,DER_FUNC,i,j) = 0.25E0_DP*&
            ( (Dpoints(1,i,j)**2-Dpoints(2,i,j)**2)+2.0_DP*Dpoints(1,i,j)+1.0_DP)
        Dbas(3,DER_FUNC,i,j) = 0.25E0_DP*&
            (-(Dpoints(1,i,j)**2-Dpoints(2,i,j)**2)+2.0_DP*Dpoints(2,i,j)+1.0_DP)
        Dbas(4,DER_FUNC,i,j) = 0.25E0_DP*&
            ( (Dpoints(1,i,j)**2-Dpoints(2,i,j)**2)-2.0_DP*Dpoints(1,i,j)+1.0_DP)
      END DO
      
    END DO
    
  END IF
    
  !if x-or y-derivatives are desired
  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN
  
    DO j=1,nelements
      dxj = 0.5E0_DP / Ddetj(:,j)
      
      !x- and y-derivatives on reference element
      DO i=1,npoints
        Dhelp(1,1,i) = -Dpoints(1,i,j)
        Dhelp(2,1,i) =  Dpoints(1,i,j)+1.0_DP
        Dhelp(3,1,i) = -Dpoints(1,i,j)
        Dhelp(4,1,i) =  Dpoints(1,i,j)-1.0_DP
        Dhelp(1,2,i) =  Dpoints(2,i,j)-1.0_DP
        Dhelp(2,2,i) = -Dpoints(2,i,j)
        Dhelp(3,2,i) =  Dpoints(2,i,j)+1.0_DP
        Dhelp(4,2,i) = -Dpoints(2,i,j)
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
! Element subroutines for parametric P0_3D element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  PURE SUBROUTINE elem_P0_3D (ieltyp, Dcoords, Djac, ddetj, Bder, &
                              Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_P0_3D.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates,
  ! Dcoords(3,.)=z-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(EL_NJACENTRIES3D), INTENT(IN) :: Djac
  
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
  REAL(DP), DIMENSION(4), INTENT(IN) :: Dpoint
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
  
  ! We have no derivatives.

  END SUBROUTINE 
  
  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_P0_3D_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                   Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_P0.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(EL_NJACENTRIES3D,npoints), INTENT(IN) :: Djac
  
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
  ! DIMENSION(4,npoints)
  !  Dpoints(1,.) = 1st barycentric coordinate
  !  Dpoints(2,.) = 2nd barycentric coordinate
  !  Dpoints(3,.) = 3rd barycentric coordinate
  !  Dpoints(4,.) = 4th barycentric coordinate
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

  PURE SUBROUTINE elem_P0_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>

  ! Element type identifier. Must be =EL_P0.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelemens)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(EL_NJACENTRIES3D,npoints), INTENT(IN) :: Djac
  
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
  ! DIMENSION(4,npoints,nelements)
  !  Dpoints(1,.) = 1st barycentric coordinate
  !  Dpoints(2,.) = 2nd barycentric coordinate
  !  Dpoints(3,.) = 3rd barycentric coordinate
  !  Dpoints(4,.) = 4th barycentric coordinate
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
! Element subroutines for parametric 3D Q0 element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  PURE SUBROUTINE elem_Q0_3D (ieltyp, Dcoords, Djac, ddetj, Bder, &
                              Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q0_3D.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates,
  ! Dcoords(3,.)=z-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(3,1)
  !  Djac(4) = J(1,2)
  !  Djac(5) = J(2,2)
  !  Djac(6) = J(3,2)
  !  Djac(7) = J(1,3)
  !  Djac(8) = J(2,3)
  !  Djac(9) = J(3,3)
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(EL_NJACENTRIES3D), INTENT(IN) :: Djac
  
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
  ! Dpoint(3) = z-coordinate
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

  PURE SUBROUTINE elem_Q0_3D_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                   Bder, Dbas, npoints, Dpoints)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q0_3D.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates,
  ! Dcoords(3,.)=z-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(3,1)
  !  Djac(4,i) = J_i(1,2)
  !  Djac(5,i) = J_i(2,2)
  !  Djac(6,i) = J_i(3,2)
  !  Djac(7,i) = J_i(1,3)
  !  Djac(8,i) = J_i(2,3)
  !  Djac(9,i) = J_I(3,3)
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(EL_NJACENTRIES3D,npoints), INTENT(IN) :: Djac
  
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
  ! DIMENSION(#space dimensions,npoints)
  ! Dpoints(1,.)=x-coordinates,
  ! Dpoints(2,.)=y-coordinates,
  ! Dpoints(3,.)=z-coordinates.
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

  PURE SUBROUTINE elem_Q0_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, nelements, Dpoints)

  !<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q0_3D.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements

  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE^,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates,
  !  Dcoords(3,.,.)=z-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(3,1,.)
  !  Djac(4,i,.) = J_i(1,2,.)
  !  Djac(5,i,.) = J_i(2,2,.)
  !  Djac(6,i,.) = J_i(3,2,.)
  !  Djac(7,i,.) = J_i(1,3,.)
  !  Djac(8,i,.) = J_i(2,3,.)
  !  Djac(9,i,.) = J_i(3,3,.)
  ! REMARK: Not used by this special type of element!
  REAL(DP), DIMENSION(EL_NJACENTRIES3D,npoints), INTENT(IN) :: Djac
  
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
  ! DIMENSION(#space dimensions,npoints,nelements)
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates,
  !  Dpoints(3,.)=z-coordinates.
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
! Element subroutines for parametric Q1 element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  PURE SUBROUTINE elem_Q1_3D (ieltyp, Dcoords, Djac, ddetj, Bder, &
                              Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q1_3D.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates,
  ! Dcoords(3,.)=z-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(3,1)
  !  Djac(4) = J(1,2)
  !  Djac(5) = J(2,2)
  !  Djac(6) = J(3,2)
  !  Djac(7) = J(1,3)
  !  Djac(8) = J(2,3)
  !  Djac(9) = J(3,3)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  REAL(DP), DIMENSION(EL_NJACENTRIES3D), INTENT(IN) :: Djac
  
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
  ! Dpoint(2) = y-coordinate,
  ! Dpoint(3) = z-coordinate
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

  !auxiliary vector containing the first derivatives on the reference element
  REAL(DP), DIMENSION(8,NDIM3D) :: Dhelp
  REAL(DP) :: dx,dy,dz, djx, djy, djz

  REAL(DP) :: dxj !auxiliary variable
  
  ! The Q1 element is specified by eight polynomials on the reference element.
  ! These eight polynomials are:
  !
  !  P1(x,y,z) = 1/8 (1-x) (1-y) (1-z)
  !  P2(x,y,z) = 1/8 (1+x) (1-y) (1-z)
  !  P3(x,y,z) = 1/8 (1+x) (1+y) (1-z)
  !  P4(x,y,z) = 1/8 (1-x) (1+y) (1-z)
  !  P5(x,y,z) = 1/8 (1-x) (1-y) (1+z)
  !  P6(x,y,z) = 1/8 (1+x) (1-y) (1+z)
  !  P7(x,y,z) = 1/8 (1+x) (1+y) (1+z)
  !  P8(x,y,z) = 1/8 (1-x) (1+y) (1+z)
  !
  ! Each of them calculated that way that Pi(Xj)=delta_ij (Kronecker)
  ! for X1,...,X8 the eight corners of the reference element
  ! [-1,1]x[-1,1]x[-1,1].
  
  ! Clear the output array
  !Dbas = 0.0_DP
  
  dx = Dpoint(1)
  dy = Dpoint(2)
  dz = Dpoint(3)
    
  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  ! If function values are desired, calculate them.
!  if (el_bder(DER_FUNC3D)) then
    Dbas(1,DER_FUNC3D) = 0.125_DP*(1.0_DP-dx)*(1.0_DP-dy)*(1.0_DP-dz)
    Dbas(2,DER_FUNC3D) = 0.125_DP*(1.0_DP+dx)*(1.0_DP-dy)*(1.0_DP-dz)
    Dbas(3,DER_FUNC3D) = 0.125_DP*(1.0_DP+dx)*(1.0_DP+dy)*(1.0_DP-dz)
    Dbas(4,DER_FUNC3D) = 0.125_DP*(1.0_DP-dx)*(1.0_DP+dy)*(1.0_DP-dz)
    Dbas(5,DER_FUNC3D) = 0.125_DP*(1.0_DP-dx)*(1.0_DP-dy)*(1.0_DP+dz)
    Dbas(6,DER_FUNC3D) = 0.125_DP*(1.0_DP+dx)*(1.0_DP-dy)*(1.0_DP+dz)
    Dbas(7,DER_FUNC3D) = 0.125_DP*(1.0_DP+dx)*(1.0_DP+dy)*(1.0_DP+dz)
    Dbas(8,DER_FUNC3D) = 0.125_DP*(1.0_DP-dx)*(1.0_DP+dy)*(1.0_DP+dz)
!  endif
  
  ! If x-, y- or z-derivatives are desired, calculate them.
  ! The values of the derivatives are calculated by taking the
  ! derivative of the polynomials and multiplying them with the
  ! inverse of the transformation matrix (in each point) as
  ! stated above.
!  if ((Bder(DER_DERIV3D_X)) .or. (Bder(DER_DERIV3D_Y)) .or. &
!      (Bder(DER_DERIV3D_Z))) then
    dxj = 0.125_DP / ddetj
    
    ! x-, y- and z-derivatives on reference element
    Dhelp(1,1) =-(1.0_DP-dy)*(1.0_DP-dz)
    Dhelp(2,1) = (1.0_DP-dy)*(1.0_DP-dz)
    Dhelp(3,1) = (1.0_DP+dy)*(1.0_DP-dz)
    Dhelp(4,1) =-(1.0_DP+dy)*(1.0_DP-dz)
    Dhelp(5,1) =-(1.0_DP-dy)*(1.0_DP+dz)
    Dhelp(6,1) = (1.0_DP-dy)*(1.0_DP+dz)
    Dhelp(7,1) = (1.0_DP+dy)*(1.0_DP+dz)
    Dhelp(8,1) =-(1.0_DP+dy)*(1.0_DP+dz)
    Dhelp(1,2) =-(1.0_DP-dx)*(1.0_DP-dz)
    Dhelp(2,2) =-(1.0_DP+dx)*(1.0_DP-dz)
    Dhelp(3,2) = (1.0_DP+dx)*(1.0_DP-dz)
    Dhelp(4,2) = (1.0_DP-dx)*(1.0_DP-dz)
    Dhelp(5,2) =-(1.0_DP-dx)*(1.0_DP+dz)
    Dhelp(6,2) =-(1.0_DP+dx)*(1.0_DP+dz)
    Dhelp(7,2) = (1.0_DP+dx)*(1.0_DP+dz)
    Dhelp(8,2) = (1.0_DP-dx)*(1.0_DP+dz)
    Dhelp(1,3) =-(1.0_DP-dx)*(1.0_DP-dy)
    Dhelp(2,3) =-(1.0_DP+dx)*(1.0_DP-dy)
    Dhelp(3,3) =-(1.0_DP+dx)*(1.0_DP+dy)
    Dhelp(4,3) =-(1.0_DP-dx)*(1.0_DP+dy)
    Dhelp(5,3) = (1.0_DP-dx)*(1.0_DP-dy)
    Dhelp(6,3) = (1.0_DP+dx)*(1.0_DP-dy)
    Dhelp(7,3) = (1.0_DP+dx)*(1.0_DP+dy)
    Dhelp(8,3) = (1.0_DP-dx)*(1.0_DP+dy)
      
    ! x-derivatives on current element
!    if (Bder(DER_DERIV3D_X)) then
      djx = Djac(5)*Djac(9) - Djac(6)*Djac(8)
      djy = Djac(8)*Djac(3) - Djac(2)*Djac(9)
      djz = Djac(2)*Djac(6) - Djac(5)*Djac(3)
      Dbas(1,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(1,1) + djy*Dhelp(1,2) + djz*Dhelp(1,3))
      Dbas(2,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(2,1) + djy*Dhelp(2,2) + djz*Dhelp(2,3))
      Dbas(3,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(3,1) + djy*Dhelp(3,2) + djz*Dhelp(3,3))
      Dbas(4,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(4,1) + djy*Dhelp(4,2) + djz*Dhelp(4,3))
      Dbas(5,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(5,1) + djy*Dhelp(5,2) + djz*Dhelp(5,3))
      Dbas(6,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(6,1) + djy*Dhelp(6,2) + djz*Dhelp(6,3))
      Dbas(7,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(7,1) + djy*Dhelp(7,2) + djz*Dhelp(7,3))
      Dbas(8,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(8,1) + djy*Dhelp(8,2) + djz*Dhelp(8,3))
!    endif
    
    ! y-derivatives on current element
!    if (Bder(DER_DERIV3D_Y)) then
      djx = Djac(7)*Djac(6) - Djac(4)*Djac(9)
      djy = Djac(1)*Djac(9) - Djac(7)*Djac(3)
      djz = Djac(4)*Djac(3) - Djac(1)*Djac(6)
      Dbas(1,DER_DERIV3D_Y) = dxj * &
          (djx*Dhelp(1,1) + djy*Dhelp(1,2) + djz*Dhelp(1,3))
      Dbas(2,DER_DERIV3D_Y) = dxj * &
          (djx*Dhelp(2,1) + djy*Dhelp(2,2) + djz*Dhelp(2,3))
      Dbas(3,DER_DERIV3D_Y) = dxj * &
          (djx*Dhelp(3,1) + djy*Dhelp(3,2) + djz*Dhelp(3,3))
      Dbas(4,DER_DERIV3D_Y) = dxj * &
          (djx*Dhelp(4,1) + djy*Dhelp(4,2) + djz*Dhelp(4,3))
      Dbas(5,DER_DERIV3D_Y) = dxj * &
          (djx*Dhelp(5,1) + djy*Dhelp(5,2) + djz*Dhelp(5,3))
      Dbas(6,DER_DERIV3D_Y) = dxj * &
          (djx*Dhelp(6,1) + djy*Dhelp(6,2) + djz*Dhelp(6,3))
      Dbas(7,DER_DERIV3D_Y) = dxj * &
          (djx*Dhelp(7,1) + djy*Dhelp(7,2) + djz*Dhelp(7,3))
      Dbas(8,DER_DERIV3D_Y) = dxj * &
          (djx*Dhelp(8,1) + djy*Dhelp(8,2) + djz*Dhelp(8,3))
!    endif

    ! z-derivatives on current element
!    if (Bder(DER_DERIV3D_Z)) then
      djx = Djac(4)*Djac(8) - Djac(7)*Djac(5)
      djy = Djac(7)*Djac(2) - Djac(1)*Djac(8)
      djz = Djac(1)*Djac(5) - Djac(4)*Djac(2)
      Dbas(1,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(1,1) + djy*Dhelp(1,2) + djz*Dhelp(1,3))
      Dbas(2,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(2,1) + djy*Dhelp(2,2) + djz*Dhelp(2,3))
      Dbas(3,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(3,1) + djy*Dhelp(3,2) + djz*Dhelp(3,3))
      Dbas(4,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(4,1) + djy*Dhelp(4,2) + djz*Dhelp(4,3))
      Dbas(5,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(5,1) + djy*Dhelp(5,2) + djz*Dhelp(5,3))
      Dbas(6,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(6,1) + djy*Dhelp(6,2) + djz*Dhelp(6,3))
      Dbas(7,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(7,1) + djy*Dhelp(7,2) + djz*Dhelp(7,3))
      Dbas(8,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(8,1) + djy*Dhelp(8,2) + djz*Dhelp(8,3))
!    endif
!  endif
    
  END SUBROUTINE 
  
  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_Q1_3D_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                   Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1_3D.
  INTEGER(I32), INTENT(IN)  :: ieltyp
  
  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates,
  ! Dcoords(3,.)=z-coordinates.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(3,1)
  !  Djac(4,i) = J_i(1,2)
  !  Djac(5,i) = J_i(2,2)
  !  Djac(6,i) = J_i(3,2)
  !  Djac(7,i) = J_i(1,3)
  !  Djac(8,i) = J_i(2,3)
  !  Djac(9,i) = J_I(3,3)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  REAL(DP), DIMENSION(EL_NJACENTRIES3D,npoints), INTENT(IN) :: Djac
  
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
  ! DIMENSION(#space dimensions,npoints).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates,
  !  Dpoints(3,.)=z-coordinates.
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
  REAL(DP), DIMENSION(8,NDIM3D,npoints) :: Dhelp

  REAL(DP),DIMENSION(npoints) :: Dxj !auxiliary variable
  REAL(DP) :: djx, djy, djz
  INTEGER :: i   ! point counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  !IF (Bder(DER_FUNC3D)) THEN
    DO i=1,npoints
      Dbas(1,DER_FUNC3D,i) = 0.125_DP*(1.0_DP-Dpoints(1,i))*&
          (1.0_DP-Dpoints(2,i))*(1.0_DP-Dpoints(3,i))
      Dbas(2,DER_FUNC3D,i) = 0.125_DP*(1.0_DP+Dpoints(1,i))*&
          (1.0_DP-Dpoints(2,i))*(1.0_DP-Dpoints(3,i))
      Dbas(3,DER_FUNC3D,i) = 0.125_DP*(1.0_DP+Dpoints(1,i))*&
          (1.0_DP+Dpoints(2,i))*(1.0_DP-Dpoints(3,i))
      Dbas(4,DER_FUNC3D,i) = 0.125_DP*(1.0_DP-Dpoints(1,i))*&
          (1.0_DP+Dpoints(2,i))*(1.0_DP-Dpoints(3,i))
      Dbas(5,DER_FUNC3D,i) = 0.125_DP*(1.0_DP-Dpoints(1,i))*&
          (1.0_DP-Dpoints(2,i))*(1.0_DP+Dpoints(3,i))
      Dbas(6,DER_FUNC3D,i) = 0.125_DP*(1.0_DP+Dpoints(1,i))*&
          (1.0_DP-Dpoints(2,i))*(1.0_DP+Dpoints(3,i))
      Dbas(7,DER_FUNC3D,i) = 0.125_DP*(1.0_DP+Dpoints(1,i))*&
          (1.0_DP+Dpoints(2,i))*(1.0_DP+Dpoints(3,i))
      Dbas(8,DER_FUNC3D,i) = 0.125_DP*(1.0_DP-Dpoints(1,i))*&
          (1.0_DP+Dpoints(2,i))*(1.0_DP+Dpoints(3,i))
    END DO
  !ENDIF
  
  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV3D_X)) .OR. (Bder(DER_DERIV3D_Y)) .OR.&
!      (Bder(DER_DERIV3D_Z))) THEN
    Dxj = 0.125E0_DP / Ddetj
    
    !x-, y- and z-derivatives on reference element
    DO i=1,npoints
      Dhelp(1,1,i) =-(1.0_DP-Dpoints(2,i))*(1.0_DP-Dpoints(3,i))
      Dhelp(2,1,i) = (1.0_DP-Dpoints(2,i))*(1.0_DP-Dpoints(3,i))
      Dhelp(3,1,i) = (1.0_DP+Dpoints(2,i))*(1.0_DP-Dpoints(3,i))
      Dhelp(4,1,i) =-(1.0_DP+Dpoints(2,i))*(1.0_DP-Dpoints(3,i))
      Dhelp(5,1,i) =-(1.0_DP-Dpoints(2,i))*(1.0_DP+Dpoints(3,i))
      Dhelp(6,1,i) = (1.0_DP-Dpoints(2,i))*(1.0_DP+Dpoints(3,i))
      Dhelp(7,1,i) = (1.0_DP+Dpoints(2,i))*(1.0_DP+Dpoints(3,i))
      Dhelp(8,1,i) =-(1.0_DP+Dpoints(2,i))*(1.0_DP+Dpoints(3,i))
      Dhelp(1,2,i) =-(1.0_DP-Dpoints(1,i))*(1.0_DP-Dpoints(3,i))
      Dhelp(2,2,i) =-(1.0_DP+Dpoints(1,i))*(1.0_DP-Dpoints(3,i))
      Dhelp(3,2,i) = (1.0_DP+Dpoints(1,i))*(1.0_DP-Dpoints(3,i))
      Dhelp(4,2,i) = (1.0_DP-Dpoints(1,i))*(1.0_DP-Dpoints(3,i))
      Dhelp(5,2,i) =-(1.0_DP-Dpoints(1,i))*(1.0_DP+Dpoints(3,i))
      Dhelp(6,2,i) =-(1.0_DP+Dpoints(1,i))*(1.0_DP+Dpoints(3,i))
      Dhelp(7,2,i) = (1.0_DP+Dpoints(1,i))*(1.0_DP+Dpoints(3,i))
      Dhelp(8,2,i) = (1.0_DP-Dpoints(1,i))*(1.0_DP+Dpoints(3,i))
      Dhelp(1,3,i) =-(1.0_DP-Dpoints(1,i))*(1.0_DP-Dpoints(2,i))
      Dhelp(2,3,i) =-(1.0_DP+Dpoints(1,i))*(1.0_DP-Dpoints(2,i))
      Dhelp(3,3,i) =-(1.0_DP+Dpoints(1,i))*(1.0_DP+Dpoints(2,i))
      Dhelp(4,3,i) =-(1.0_DP-Dpoints(1,i))*(1.0_DP+Dpoints(2,i))
      Dhelp(5,3,i) = (1.0_DP-Dpoints(1,i))*(1.0_DP-Dpoints(2,i))
      Dhelp(6,3,i) = (1.0_DP+Dpoints(1,i))*(1.0_DP-Dpoints(2,i))
      Dhelp(7,3,i) = (1.0_DP+Dpoints(1,i))*(1.0_DP+Dpoints(2,i))
      Dhelp(8,3,i) = (1.0_DP-Dpoints(1,i))*(1.0_DP+Dpoints(2,i))
    END DO
      
    !x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      DO i=1,npoints
        djx = Djac(5,i)*Djac(9,i) - Djac(6,i)*Djac(8,i)
        djy = Djac(8,i)*Djac(3,i) - Djac(2,i)*Djac(9,i)
        djz = Djac(2,i)*Djac(6,i) - Djac(5,i)*Djac(3,i)
        Dbas(1,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(1,1,i) + djy*Dhelp(1,2,i) + djz*Dhelp(1,3,i))
        Dbas(2,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(2,1,i) + djy*Dhelp(2,2,i) + djz*Dhelp(2,3,i))
        Dbas(3,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(3,1,i) + djy*Dhelp(3,2,i) + djz*Dhelp(3,3,i))
        Dbas(4,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(4,1,i) + djy*Dhelp(4,2,i) + djz*Dhelp(4,3,i))
        Dbas(5,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(5,1,i) + djy*Dhelp(5,2,i) + djz*Dhelp(5,3,i))
        Dbas(6,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(6,1,i) + djy*Dhelp(6,2,i) + djz*Dhelp(6,3,i))
        Dbas(7,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(7,1,i) + djy*Dhelp(7,2,i) + djz*Dhelp(7,3,i))
        Dbas(8,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(8,1,i) + djy*Dhelp(8,2,i) + djz*Dhelp(8,3,i))
!      END DO
!    ENDIF
    
    !y-derivatives on current element
!    IF (Bder(DER_DERIV_Y)) THEN
!      DO i=1,npoints
        djx = Djac(7,i)*Djac(6,i) - Djac(4,i)*Djac(9,i)
        djy = Djac(1,i)*Djac(9,i) - Djac(7,i)*Djac(3,i)
        djz = Djac(4,i)*Djac(3,i) - Djac(1,i)*Djac(6,i)
        Dbas(1,DER_DERIV3D_Y,i) = Dxj(i) * &
            (djx*Dhelp(1,1,i) + djy*Dhelp(1,2,i) + djz*Dhelp(1,3,i))
        Dbas(2,DER_DERIV3D_Y,i) = Dxj(i) * &
            (djx*Dhelp(2,1,i) + djy*Dhelp(2,2,i) + djz*Dhelp(2,3,i))
        Dbas(3,DER_DERIV3D_Y,i) = Dxj(i) * &
            (djx*Dhelp(3,1,i) + djy*Dhelp(3,2,i) + djz*Dhelp(3,3,i))
        Dbas(4,DER_DERIV3D_Y,i) = Dxj(i) * &
            (djx*Dhelp(4,1,i) + djy*Dhelp(4,2,i) + djz*Dhelp(4,3,i))
        Dbas(5,DER_DERIV3D_Y,i) = Dxj(i) * &
            (djx*Dhelp(5,1,i) + djy*Dhelp(5,2,i) + djz*Dhelp(5,3,i))
        Dbas(6,DER_DERIV3D_Y,i) = Dxj(i) * &
            (djx*Dhelp(6,1,i) + djy*Dhelp(6,2,i) + djz*Dhelp(6,3,i))
        Dbas(7,DER_DERIV3D_Y,i) = Dxj(i) * &
            (djx*Dhelp(7,1,i) + djy*Dhelp(7,2,i) + djz*Dhelp(7,3,i))
        Dbas(8,DER_DERIV3D_Y,i) = Dxj(i) * &
            (djx*Dhelp(8,1,i) + djy*Dhelp(8,2,i) + djz*Dhelp(8,3,i))
!      END DO
!    ENDIF

    !z-derivatives on current element
!    IF (Bder(DER_DERIV3D_Z)) THEN
!      DO i=1,npoints
        djx = Djac(4,i)*Djac(8,i) - Djac(7,i)*Djac(5,i)
        djy = Djac(7,i)*Djac(2,i) - Djac(1,i)*Djac(8,i)
        djz = Djac(1,i)*Djac(5,i) - Djac(4,i)*Djac(2,i)
        Dbas(1,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(1,1,i) + djy*Dhelp(1,2,i) + djz*Dhelp(1,3,i))
        Dbas(2,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(2,1,i) + djy*Dhelp(2,2,i) + djz*Dhelp(2,3,i))
        Dbas(3,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(3,1,i) + djy*Dhelp(3,2,i) + djz*Dhelp(3,3,i))
        Dbas(4,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(4,1,i) + djy*Dhelp(4,2,i) + djz*Dhelp(4,3,i))
        Dbas(5,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(5,1,i) + djy*Dhelp(5,2,i) + djz*Dhelp(5,3,i))
        Dbas(6,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(6,1,i) + djy*Dhelp(6,2,i) + djz*Dhelp(6,3,i))
        Dbas(7,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(7,1,i) + djy*Dhelp(7,2,i) + djz*Dhelp(7,3,i))
        Dbas(8,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(8,1,i) + djy*Dhelp(8,2,i) + djz*Dhelp(8,3,i))
      END DO
!    ENDIF
!  ENDIF
    
  END SUBROUTINE 

  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_Q1_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                                  Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1_3D.
  INTEGER(I32), INTENT(IN)  :: ieltyp

  ! Number of points on every element where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  INTEGER, INTENT(IN)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements).
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates,
  !  Dcoords(3,.,.)=z-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(3,1,.)
  !  Djac(4,i,.) = J_i(1,2,.)
  !  Djac(5,i,.) = J_i(2,2,.)
  !  Djac(6,i,.) = J_i(3,2,.)
  !  Djac(7,i,.) = J_i(1,3,.)
  !  Djac(8,i,.) = J_i(2,3,.)
  !  Djac(9,i,.) = J_i(3,3,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
  REAL(DP), DIMENSION(EL_NJACENTRIES3D,npoints,nelements), INTENT(IN) :: Djac
  
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
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates,
  !  Dpoints(3,.)=z-coordinates.
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
  REAL(DP), DIMENSION(8,NDIM3D,npoints) :: Dhelp
  REAL(DP) :: djx, djy, djz
  REAL(DP),DIMENSION(npoints) :: Dxj !auxiliary variable
  
  INTEGER :: i   ! point counter
  INTEGER :: j   ! element counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  IF (Bder(DER_FUNC3D)) THEN
  
    DO j=1,nelements
    
      DO i=1,npoints
        Dbas(1,DER_FUNC3D,i,j) = 0.125_DP*(1.0_DP-Dpoints(1,i,j))*&
            (1.0_DP-Dpoints(2,i,j))*(1.0_DP-Dpoints(3,i,j))
        Dbas(2,DER_FUNC3D,i,j) = 0.125_DP*(1.0_DP+Dpoints(1,i,j))*&
            (1.0_DP-Dpoints(2,i,j))*(1.0_DP-Dpoints(3,i,j))
        Dbas(3,DER_FUNC3D,i,j) = 0.125_DP*(1.0_DP+Dpoints(1,i,j))*&
            (1.0_DP+Dpoints(2,i,j))*(1.0_DP-Dpoints(3,i,j))
        Dbas(4,DER_FUNC3D,i,j) = 0.125_DP*(1.0_DP-Dpoints(1,i,j))*&
            (1.0_DP+Dpoints(2,i,j))*(1.0_DP-Dpoints(3,i,j))
        Dbas(5,DER_FUNC3D,i,j) = 0.125_DP*(1.0_DP-Dpoints(1,i,j))*&
            (1.0_DP-Dpoints(2,i,j))*(1.0_DP+Dpoints(3,i,j))
        Dbas(6,DER_FUNC3D,i,j) = 0.125_DP*(1.0_DP+Dpoints(1,i,j))*&
            (1.0_DP-Dpoints(2,i,j))*(1.0_DP+Dpoints(3,i,j))
        Dbas(7,DER_FUNC3D,i,j) = 0.125_DP*(1.0_DP+Dpoints(1,i,j))*&
            (1.0_DP+Dpoints(2,i,j))*(1.0_DP+Dpoints(3,i,j))
        Dbas(8,DER_FUNC3D,i,j) = 0.125_DP*(1.0_DP-Dpoints(1,i,j))*&
            (1.0_DP+Dpoints(2,i,j))*(1.0_DP+Dpoints(3,i,j))
      END DO
      
    END DO
    
  END IF
    
  !if x-, y- or z-derivatives are desired
  IF ((Bder(DER_DERIV3D_X)) .OR. (Bder(DER_DERIV3D_Y)) .OR. &
      (Bder(DER_DERIV3D_Z))) THEN
  
    DO j=1,nelements
      Dxj = 0.125_DP / Ddetj(:,j)
      
      !x-, y- and z-derivatives on reference element
      DO i=1,npoints
        Dhelp(1,1,i) =-(1.0_DP-Dpoints(2,i,j))*(1.0_DP-Dpoints(3,i,j))
        Dhelp(2,1,i) = (1.0_DP-Dpoints(2,i,j))*(1.0_DP-Dpoints(3,i,j))
        Dhelp(3,1,i) = (1.0_DP+Dpoints(2,i,j))*(1.0_DP-Dpoints(3,i,j))
        Dhelp(4,1,i) =-(1.0_DP+Dpoints(2,i,j))*(1.0_DP-Dpoints(3,i,j))
        Dhelp(5,1,i) =-(1.0_DP-Dpoints(2,i,j))*(1.0_DP+Dpoints(3,i,j))
        Dhelp(6,1,i) = (1.0_DP-Dpoints(2,i,j))*(1.0_DP+Dpoints(3,i,j))
        Dhelp(7,1,i) = (1.0_DP+Dpoints(2,i,j))*(1.0_DP+Dpoints(3,i,j))
        Dhelp(8,1,i) =-(1.0_DP+Dpoints(2,i,j))*(1.0_DP+Dpoints(3,i,j))
        Dhelp(1,2,i) =-(1.0_DP-Dpoints(1,i,j))*(1.0_DP-Dpoints(3,i,j))
        Dhelp(2,2,i) =-(1.0_DP+Dpoints(1,i,j))*(1.0_DP-Dpoints(3,i,j))
        Dhelp(3,2,i) = (1.0_DP+Dpoints(1,i,j))*(1.0_DP-Dpoints(3,i,j))
        Dhelp(4,2,i) = (1.0_DP-Dpoints(1,i,j))*(1.0_DP-Dpoints(3,i,j))
        Dhelp(5,2,i) =-(1.0_DP-Dpoints(1,i,j))*(1.0_DP+Dpoints(3,i,j))
        Dhelp(6,2,i) =-(1.0_DP+Dpoints(1,i,j))*(1.0_DP+Dpoints(3,i,j))
        Dhelp(7,2,i) = (1.0_DP+Dpoints(1,i,j))*(1.0_DP+Dpoints(3,i,j))
        Dhelp(8,2,i) = (1.0_DP-Dpoints(1,i,j))*(1.0_DP+Dpoints(3,i,j))
        Dhelp(1,3,i) =-(1.0_DP-Dpoints(1,i,j))*(1.0_DP-Dpoints(2,i,j))
        Dhelp(2,3,i) =-(1.0_DP+Dpoints(1,i,j))*(1.0_DP-Dpoints(2,i,j))
        Dhelp(3,3,i) =-(1.0_DP+Dpoints(1,i,j))*(1.0_DP+Dpoints(2,i,j))
        Dhelp(4,3,i) =-(1.0_DP-Dpoints(1,i,j))*(1.0_DP+Dpoints(2,i,j))
        Dhelp(5,3,i) = (1.0_DP-Dpoints(1,i,j))*(1.0_DP-Dpoints(2,i,j))
        Dhelp(6,3,i) = (1.0_DP+Dpoints(1,i,j))*(1.0_DP-Dpoints(2,i,j))
        Dhelp(7,3,i) = (1.0_DP+Dpoints(1,i,j))*(1.0_DP+Dpoints(2,i,j))
        Dhelp(8,3,i) = (1.0_DP-Dpoints(1,i,j))*(1.0_DP+Dpoints(2,i,j))
      END DO
        
      !x-derivatives on current element
!      IF (Bder(DER_DERIV3D_X)) THEN
        DO i=1,npoints
          djx = Djac(5,i,j)*Djac(9,i,j) - Djac(6,i,j)*Djac(8,i,j)
          djy = Djac(8,i,j)*Djac(3,i,j) - Djac(2,i,j)*Djac(9,i,j)
          djz = Djac(2,i,j)*Djac(6,i,j) - Djac(5,i,j)*Djac(3,i,j)
          Dbas(1,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(1,1,i) + djy*Dhelp(1,2,i) + djz*Dhelp(1,3,i))
          Dbas(2,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(2,1,i) + djy*Dhelp(2,2,i) + djz*Dhelp(2,3,i))
          Dbas(3,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(3,1,i) + djy*Dhelp(3,2,i) + djz*Dhelp(3,3,i))
          Dbas(4,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(4,1,i) + djy*Dhelp(4,2,i) + djz*Dhelp(4,3,i))
          Dbas(5,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(5,1,i) + djy*Dhelp(5,2,i) + djz*Dhelp(5,3,i))
          Dbas(6,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(6,1,i) + djy*Dhelp(6,2,i) + djz*Dhelp(6,3,i))
          Dbas(7,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(7,1,i) + djy*Dhelp(7,2,i) + djz*Dhelp(7,3,i))
          Dbas(8,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(8,1,i) + djy*Dhelp(8,2,i) + djz*Dhelp(8,3,i))
        END DO
!      ENDIF
      
      !y-derivatives on current element
!      IF (Bder(DER_DERIV3D_Y)) THEN
        DO i=1,npoints
          djx = Djac(7,i,j)*Djac(6,i,j) - Djac(4,i,j)*Djac(9,i,j)
          djy = Djac(1,i,j)*Djac(9,i,j) - Djac(7,i,j)*Djac(3,i,j)
          djz = Djac(4,i,j)*Djac(3,i,j) - Djac(1,i,j)*Djac(6,i,j)
          Dbas(1,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (djx*Dhelp(1,1,i) + djy*Dhelp(1,2,i) + djz*Dhelp(1,3,i))
          Dbas(2,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (djx*Dhelp(2,1,i) + djy*Dhelp(2,2,i) + djz*Dhelp(2,3,i))
          Dbas(3,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (djx*Dhelp(3,1,i) + djy*Dhelp(3,2,i) + djz*Dhelp(3,3,i))
          Dbas(4,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (djx*Dhelp(4,1,i) + djy*Dhelp(4,2,i) + djz*Dhelp(4,3,i))
          Dbas(5,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (djx*Dhelp(5,1,i) + djy*Dhelp(5,2,i) + djz*Dhelp(5,3,i))
          Dbas(6,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (djx*Dhelp(6,1,i) + djy*Dhelp(6,2,i) + djz*Dhelp(6,3,i))
          Dbas(7,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (djx*Dhelp(7,1,i) + djy*Dhelp(7,2,i) + djz*Dhelp(7,3,i))
          Dbas(8,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (djx*Dhelp(8,1,i) + djy*Dhelp(8,2,i) + djz*Dhelp(8,3,i))
        END DO
!      ENDIF

      !z-derivatives on current element
!      IF (Bder(DER_DERIV3D_Z)) THEN
        DO i=1,npoints
          djx = Djac(4,i,j)*Djac(8,i,j) - Djac(7,i,j)*Djac(5,i,j)
          djy = Djac(7,i,j)*Djac(2,i,j) - Djac(1,i,j)*Djac(8,i,j)
          djz = Djac(1,i,j)*Djac(5,i,j) - Djac(4,i,j)*Djac(2,i,j)
          Dbas(1,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(1,1,i) + djy*Dhelp(1,2,i) + djz*Dhelp(1,3,i))
          Dbas(2,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(2,1,i) + djy*Dhelp(2,2,i) + djz*Dhelp(2,3,i))
          Dbas(3,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(3,1,i) + djy*Dhelp(3,2,i) + djz*Dhelp(3,3,i))
          Dbas(4,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(4,1,i) + djy*Dhelp(4,2,i) + djz*Dhelp(4,3,i))
          Dbas(5,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(5,1,i) + djy*Dhelp(5,2,i) + djz*Dhelp(5,3,i))
          Dbas(6,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(6,1,i) + djy*Dhelp(6,2,i) + djz*Dhelp(6,3,i))
          Dbas(7,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(7,1,i) + djy*Dhelp(7,2,i) + djz*Dhelp(7,3,i))
          Dbas(8,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(8,1,i) + djy*Dhelp(8,2,i) + djz*Dhelp(8,3,i))
        END DO
!      ENDIF
    END DO
      
  END IF
    
  END SUBROUTINE 

END MODULE 
