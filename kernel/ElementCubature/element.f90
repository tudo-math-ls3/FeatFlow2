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
!#                               solving local 4x4 / 6x6 systems; faster but less
!#                               stable on cruel elements.
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
  
  ! ID of cubic conforming line FE, 3,1-Spline
  INTEGER(I32), PARAMETER :: EL_S31_1D = EL_1D + 17
  
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

  ! General rotated linear $\tilde P1$ element (Crouzeix-Raviart)
  INTEGER(I32), PARAMETER :: EL_P1T  = EL_2D + 20

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

  ! General rotated biquadratic $\tilde Q2$ element, all variants.
  INTEGER(I32), PARAMETER :: EL_Q2T  = EL_2D + 35
  
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

  ! ID of rotated trilinear non-conforming hexahedral FE, Q1~
  ! parametric, face-midpoint-value based
  INTEGER(I32), PARAMETER :: EL_Q1T_3D = EL_3D + 30
  INTEGER(I32), PARAMETER :: EL_E031_3D = EL_Q1T_3D
  
  ! ID of rotated trilinear non-conforming hexahedral FE, Q1~
  ! parametric, integral mean value based
  INTEGER(I32), PARAMETER :: EL_E030_3D = EL_Q1T_3D + 2**16

  ! ID of rotated trilinear non-conforming hexahedral FE, Q1~
  ! non-parametric, integral mean value based
  INTEGER(I32), PARAMETER :: EL_EM30_3D = EL_Q1T_3D + EL_NONPARAMETRIC + 2**16

!</constantblock>

!<constantblock description="Special 2D element variants.">
  
  ! ID of rotated linear nonconforming triangular FE, P1~, edge-midpoint based
  INTEGER(I32), PARAMETER :: EL_E020 = EL_P1T + 2**16

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

  ! ID of rotated biquadratic nonconforming quadrilateral FE, Q2~.
  INTEGER(I32), PARAMETER :: EL_E035 = EL_Q2T
  
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
  CASE (EL_S31_1D)
    ! local DOF's for 1D S31
    elem_igetNDofLoc = 4
  
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
  CASE (EL_P1T)
    ! local DOF's for Ex20
    elem_igetNDofLoc = 3
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
  CASE (EL_Q2T)
    ! local DOF's for Ex35
    elem_igetNDofLoc = 8
  
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
  CASE (EL_Q1T_3D)
    ! local DOF's for 3D Ex30
    elem_igetNDofLoc = 6
    
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
  CASE (EL_S31_1D)
    ! local DOF's for S31
    ndofAtVertices = 4
    
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
  CASE (EL_P1T)
    ! local DOF's for Ex20
    ndofAtEdges    = 3
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
  CASE (EL_Q2T)
    ! local DOF's for Ex35
    ndofAtEdges    = 8
  
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
  CASE (EL_Q1T_3D)
    ! local DOF's for Ex30
    ndofAtFaces = 6
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
  ! 1D Element types
  CASE (EL_P0_1D,EL_P1_1D,EL_P2_1D,EL_S31_1D)
    ! Line elements
    elem_igetCoordSystem = TRAFO_CS_REF1D
  
  ! 2D Element Types
  CASE (EL_P0,EL_P1,EL_P2,EL_P3,EL_P1T)
    ! Triangular elements work in barycentric coordinates
    elem_igetCoordSystem = TRAFO_CS_BARY2DTRI
  CASE (EL_Q0,EL_Q1,EL_Q2,EL_Q3,EL_QP1,EL_Q1T,EL_Q2T)
    ! These work on the reference quadrilateral
    elem_igetCoordSystem = TRAFO_CS_REF2DQUAD
  CASE (EL_Q1T+EL_NONPARAMETRIC)
    ! EM30, EM31; these work in real coordinates
    elem_igetCoordSystem = TRAFO_CS_REAL2DQUAD
  
  ! 3D Element types
  CASE (EL_P0_3D, EL_P1_3D)
    ! Tetrahedral elements work in barycentric coordinates
    elem_igetCoordSystem = TRAFO_CS_BARY3DTETRA
  CASE (EL_Q0_3D, EL_Q1_3D, EL_Q1T_3D)
    ! These work on the reference hexahedron
    elem_igetCoordSystem = TRAFO_CS_REF3DHEXA
  CASE (EL_Q1T_3D+EL_NONPARAMETRIC)
    ! EM30; these work in real coordinates
    elem_igetCoordSystem = TRAFO_CS_REAL3DHEXA
    
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
  
  CASE (EL_P0_1D, EL_P1_1D, EL_P2_1D, EL_S31_1D)
    ! Linear line transformation, 1D
    elem_igetTrafoType = TRAFO_ID_LINSIMPLEX + TRAFO_DIM_1D
  
  CASE (EL_P0, EL_P1, EL_P2, EL_P3, EL_P1T)
    ! Linear triangular transformation, 2D
    elem_igetTrafoType = TRAFO_ID_LINSIMPLEX + TRAFO_DIM_2D
    
  CASE (EL_Q0,EL_Q1,EL_Q2,EL_Q3,EL_QP1,EL_Q1T,EL_Q2T)
    ! Bilinear quadrilateral transformation, 2D.
    elem_igetTrafoType = TRAFO_ID_MLINCUBE + TRAFO_DIM_2D
  
  CASE (EL_P0_3D, EL_P1_3D)
    ! Linear tetrahedral transrormation, 3D
    elem_igetTrafoType = TRAFO_ID_LINSIMPLEX + TRAFO_DIM_3D
  
  CASE (EL_Q0_3D, EL_Q1_3D, EL_Q1T_3D)
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
    ! Function + 1st derivative
    elem_getMaxDerivative = 2
  CASE (EL_S31_1D)
    ! Function + 1st derivative
    elem_getMaxDerivative = 2
  CASE (EL_P0, EL_Q0)
    ! Function
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
  CASE (EL_P1T)
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
  CASE (EL_Q2T)
    ! Function + 1st derivative
    elem_getMaxDerivative = 3
  CASE (EL_P0_3D, EL_Q0_3D)
    ! Function
    elem_getMaxDerivative = 1
  CASE (EL_P1_3D)
    ! Function + 1st derivative
    elem_getMaxDerivative = 4
  CASE (EL_Q1_3D)
    ! Function + 1st derivative
    elem_getMaxDerivative = 4
  CASE (EL_Q1T_3D)
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
  CASE (EL_S31_1D)
    CALL elem_S31_1D (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  ! 2D elements
  CASE (EL_P0)
    CALL elem_P0 (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_P1)
    CALL elem_P1 (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_P2)
    CALL elem_P2 (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_P1T)
    CALL elem_P1T (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
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
  CASE (EL_E035)
    CALL elem_E035 (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  ! 3D elements
  CASE (EL_P0_3D)
    CALL elem_P0_3D (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_P1_3D)
    CALL elem_P1_3D (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_Q0_3D)
    CALL elem_Q0_3D (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_Q1_3D)
    CALL elem_Q1_3D (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_EM30_3D)
    CALL elem_EM30_3D (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_E030_3D)
    CALL elem_E030_3D (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
  CASE (EL_E031_3D)
    CALL elem_E031_3D (ieltyp, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
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
  CASE (EL_S31_1D)
    CALL elem_S31_1D_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_P0)
    CALL elem_P0_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_P1)
    CALL elem_P1_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_P2)
    CALL elem_P2_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_P1T)
    CALL elem_P1T_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
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
  CASE (EL_P0_3D)
    CALL elem_P0_3D_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_P1_3D)
    CALL elem_P1_3D_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_Q0_3D)
    CALL elem_Q0_3D_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_Q1_3D)
    CALL elem_Q1_3D_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_EM30_3D)
    CALL elem_EM30_3D_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_E030_3D)
    CALL elem_E030_3D_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
  CASE (EL_E031_3D)
    CALL elem_E031_3D_mult (ieltyp, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
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
  CASE (EL_S31_1D)
    CALL elem_S31_1D_sim (ieltyp, Dcoords, Djac, Ddetj, &
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
  CASE (EL_P1T)
    CALL elem_P1T_sim (ieltyp, Dcoords, Djac, Ddetj, &
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
  CASE (EL_P0_3D)
    CALL elem_P0_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                      Bder, Dbas, npoints, nelements, Dpoints)
  CASE (EL_P1_3D)
    CALL elem_P1_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                      Bder, Dbas, npoints, nelements, Dpoints)
  CASE (EL_Q0_3D)
    CALL elem_Q0_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                      Bder, Dbas, npoints, nelements, Dpoints)
  CASE (EL_Q1_3D)
    CALL elem_Q1_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                      Bder, Dbas, npoints, nelements, Dpoints)
  CASE (EL_EM30_3D)
    CALL elem_EM30_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
  CASE (EL_E030_3D)
    CALL elem_E030_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                      Bder, Dbas, npoints, nelements, Dpoints)
  CASE (EL_E031_3D)
    CALL elem_E031_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
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
  DBas(1,DER_FUNC1D) = 1.0_DP

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
  DBas(1,DER_FUNC1D,:) = 1.0_DP

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
  DBas(1,DER_FUNC1D,:,:) = 1.0_DP

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
  !  P1(x) = 1/2 (1-x)
  !  P2(x) = 1/2 (1+x)
  !
  ! Each of them calculated that way that Pi(Xj)=delta_ij (Kronecker)
  ! for X1,X2 the two corners of the reference element [-1,1].
  
  ! Clear the output array
  !Dbas = 0.0_DP
  
  ! Remark: The P1_1D-element always computes function value and 1st derivative.
  ! That's even faster than when using two IF commands for preventing
  ! the computation of one of the values!
      
  ! If function values are desired, calculate them.
!  if (el_bder(DER_FUNC1D)) then
    Dbas(1,DER_FUNC1D) = 0.5_DP * (1.0_DP - Dpoint(1))
    Dbas(2,DER_FUNC1D) = 0.5_DP * (1.0_DP + Dpoint(1))
!  endif
  
  ! If x-derivatives are desired, calculate them.
  ! Since P1_1D is linear, the first derivative is constant!
!  if (Bder(DER_DERIV1D_X)) then
    Dbas(1,DER_DERIV1D_X) = -0.5_DP / ddetj
    Dbas(2,DER_DERIV1D_X) = 0.5_DP / ddetj
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
  !IF (Bder(DER_FUNC1D)) THEN
    DO i=1,npoints
      Dbas(1,DER_FUNC1D,i) = 0.5_DP * (1_DP - Dpoints(1,i))
      Dbas(2,DER_FUNC1D,i) = 0.5_DP * (1_DP + Dpoints(1,i))
    END DO
  !ENDIF
  
  !if x-derivatives are desired
!  IF (Bder(DER_DERIV1D_X)) THEN
    DO i=1,npoints
      Dbas(1,DER_DERIV1D_X,i) = -0.5_DP / Ddetj(i)
      Dbas(2,DER_DERIV1D_X,i) = 0.5_DP / Ddetj(i)
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
  IF (Bder(DER_FUNC1D)) THEN
  
    DO j=1,nelements
    
      DO i=1,npoints
        Dbas(1,DER_FUNC1D,i,j) = 0.5_DP * (1.0_DP - Dpoints(1,i,j))
        Dbas(2,DER_FUNC1D,i,j) = 0.5_DP * (1.0_DP + Dpoints(1,i,j))
      END DO
      
    END DO
    
  END IF
    
  !if x-derivatives are desired
  IF (Bder(DER_DERIV1D_X)) THEN
  
    DO j=1,nelements
    
      DO i=1,npoints
        Dbas(1,DER_DERIV1D_X,i,j) = -0.5 / Ddetj(i,j)
        Dbas(2,DER_DERIV1D_X,i,j) = 0.5 / Ddetj(i,j)
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
!  if (el_bder(DER_FUNC1D)) then
    Dbas(1,DER_FUNC1D) = 0.5_DP * Dpoint(1) * (Dpoint(1) - 1.0_DP)
    Dbas(2,DER_FUNC1D) = 0.5_DP * Dpoint(1) * (Dpoint(1) + 1.0_DP)
    Dbas(3,DER_FUNC1D) = 1.0_DP - Dpoint(1)**2
!  endif
  
  ! If x-derivatives are desired, calculate them.
!  if (Bder(DER_DERIV1D_X)) then
    d = 1.0_DP / ddetj
    Dbas(1,DER_DERIV1D_X) = (Dpoint(1) - 0.5_DP) * d
    Dbas(2,DER_DERIV1D_X) = (Dpoint(1) + 0.5_DP) * d
    Dbas(3,DER_DERIV1D_X) = -2.0_DP * Dpoint(1) * d
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
  !IF (Bder(DER_FUNC1D)) THEN
    DO i=1,npoints
      Dbas(1,DER_FUNC1D,i) = 0.5_DP * Dpoints(1,i) * (Dpoints(1,i) - 1.0_DP)
      Dbas(2,DER_FUNC1D,i) = 0.5_DP * Dpoints(1,i) * (Dpoints(1,i) + 1.0_DP)
      Dbas(3,DER_FUNC1D,i) = 1.0_DP - Dpoints(1,i)**2
    END DO
  !ENDIF
 
  !if x-derivatives are desired
!  IF (Bder(DER_DERIV1D_X)) THEN
    DO i=1,npoints
      d = 1.0_DP / Ddetj(i)
      Dbas(1,DER_DERIV1D_X,i) = (Dpoints(1,i) - 0.5_DP) * d
      Dbas(2,DER_DERIV1D_X,i) = (Dpoints(1,i) + 0.5_DP) * d
      Dbas(3,DER_DERIV1D_X,i) = -2.0_DP * Dpoints(1,i) * d
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
  IF (Bder(DER_FUNC1D)) THEN
  
    DO j=1,nelements
    
      DO i=1,npoints
        Dbas(1,DER_FUNC1D,i,j) = 0.5_DP*Dpoints(1,i,j)*(Dpoints(1,i,j) - 1.0_DP)
        Dbas(2,DER_FUNC1D,i,j) = 0.5_DP*Dpoints(1,i,j)*(Dpoints(1,i,j) + 1.0_DP)
        Dbas(3,DER_FUNC1D,i,j) = 1.0_DP - Dpoints(1,i,j)**2
      END DO
      
    END DO
    
  END IF
    
  !if x-derivatives are desired
  IF (Bder(DER_DERIV1D_X)) THEN
  
    DO j=1,nelements
    
      DO i=1,npoints
        d = 1.0_DP / Ddetj(i,j)
        Dbas(1,DER_DERIV1D_X,i,j) = (Dpoints(1,i,j) - 0.5_DP) * d
        Dbas(2,DER_DERIV1D_X,i,j) = (Dpoints(1,i,j) + 0.5_DP) * d
        Dbas(3,DER_DERIV1D_X,i,j) = -2.0_DP * Dpoints(1,i,j) * d
      END DO

    END DO
      
  END IF
    
  END SUBROUTINE 

!**************************************************************************
! Element subroutines for parametric S31_1D element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  PURE SUBROUTINE elem_S31_1D (ieltyp, Dcoords, Djac, ddetj, Bder, &
                               Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_S31_1D.
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

  ! auxiliary variable
  REAL(DP) :: dxj

  ! The S31_1D element is specified by four polynomials on the reference element.
  ! These four polynomials are:
  !
  ! P1(x) = 1/4 * x * (x^2 - 3) + 1/2
  ! P2(x) = 1/4 * x * (3 - x^2) + 1/2
  ! Q1(x) = 1/4 * (x + 1) * (x - 1)^2
  ! Q2(x) = 1/4 * (x - 1) * (x + 1)^2
  !
  ! The basis functions are constructed that way that they fulfill the
  ! following conditions on the reference line [-1,1]:
  ! P1(-1) = 1    P1(1) = 0    P1_x(-1) = 0    P1_x(1) = 0
  ! P2(-1) = 0    P2(1) = 1    P2_x(-1) = 0    P2_x(1) = 0
  ! Q1(-1) = 0    Q1(1) = 0    Q1_x(-1) = 1    Q1_x(1) = 0
  ! Q2(-1) = 0    Q2(1) = 0    Q2_x(-1) = 0    Q2_x(1) = 1
  !
  ! Now if the FEM function is transformed onto a "real" line [a,b], the
  ! transformed basis functions pi, qi fulfill the following conditions:
  ! p1(a) = 1    p1(b) = 0    p1_x(a) = 0    p1_x(b) = 0
  ! p2(a) = 0    p2(b) = 1    p2_x(a) = 0    p2_x(b) = 0
  ! q1(a) = 0    q1(b) = 0    q1_x(a) = L    q1_x(b) = 0
  ! q2(a) = 0    q2(b) = 0    q2_x(a) = 0    q2_x(b) = L
  !
  ! where L = 2/(b-a)
  !
  ! We now want to enforce that q1_x(a)=1 and q2_x(b)=1.
  ! To do this, we need multiply the reference basis functions Q1 and Q2
  ! by 1/L. Now L is exactly the inverse of the determinant of the Jacobian
  ! matrix of the linear line transformation (see transformation.f90), i.e.
  ! 1/L = ddetj. So all we need to do is to multiply Q1 and Q2 by ddetj.
      
  ! If function values are desired, calculate them.
!  if (el_bder(DER_FUNC1D)) then
    ! P1, P2
    Dbas(1,DER_FUNC1D) = 0.25_DP*Dpoint(1)*(Dpoint(1)**2 - 3.0_DP) + 0.5_DP
    Dbas(2,DER_FUNC1D) = 0.25_DP*Dpoint(1)*(3.0_DP - Dpoint(1)**2) + 0.5_DP
    ! Q1, Q2
    Dbas(3,DER_FUNC1D) = ddetj*0.25*(Dpoint(1) + 1.0_DP)*(Dpoint(1) - 1.0_DP)**2
    Dbas(4,DER_FUNC1D) = ddetj*0.25*(Dpoint(1) - 1.0_DP)*(Dpoint(1) + 1.0_DP)**2
!  endif

  ! If x-derivatives are desired, calculate them.
!  if (Bder(DER_DERIV1D_X)) then
    dxj = 1.0_DP / ddetj
    ! P1, P2
    Dbas(1,DER_DERIV1D_X) = 0.75_DP*(Dpoint(1)**2 - 1.0_DP)*dxj
    Dbas(2,DER_DERIV1D_X) = 0.75_DP*(1.0_DP - Dpoint(1)**2)*dxj
    ! Q1, Q2
    Dbas(3,DER_DERIV1D_X) = 0.25_DP*(3.0_DP*Dpoint(1)**2 - Dpoint(1) - 1.0_DP)
    Dbas(4,DER_DERIV1D_X) = 0.25_DP*(3.0_DP*Dpoint(1)**2 + Dpoint(1) - 1.0_DP)
!  endif
    
  END SUBROUTINE 
  
  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_S31_1D_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                    Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_S31_1D.
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
  REAL(DP) :: dxj
    
  ! Clear the output array
  !Dbas = 0.0_DP

  !if function values are desired
  IF (Bder(DER_FUNC1D)) THEN
    DO i=1,npoints
      Dbas(1,DER_FUNC1D,i) = 0.25_DP*Dpoints(1,i)*&
                            (Dpoints(1,i)**2 - 3.0_DP) + 0.5_DP
      Dbas(2,DER_FUNC1D,i) = 0.25_DP*Dpoints(1,i)*&
                            (3.0_DP - Dpoints(1,i)**2) + 0.5_DP
      Dbas(3,DER_FUNC1D,i) = Ddetj(i)*0.25*(Dpoints(1,i) + 1.0_DP)*&
                            (Dpoints(1,i) - 1.0_DP)**2
      Dbas(4,DER_FUNC1D,i) = Ddetj(i)*0.25*(Dpoints(1,i) - 1.0_DP)*&
                            (Dpoints(1,i) + 1.0_DP)**2
    END DO
  ENDIF
  
  !if x-derivatives are desired
  IF (Bder(DER_DERIV1D_X)) THEN
    DO i=1,npoints
      dxj = 1.0_DP / Ddetj(i)
      Dbas(1,DER_DERIV1D_X,i) = 0.75_DP*(Dpoints(1,i)**2 - 1.0_DP)*dxj
      Dbas(2,DER_DERIV1D_X,i) = 0.75_DP*(1.0_DP - Dpoints(1,i)**2)*dxj
      Dbas(3,DER_DERIV1D_X,i) = 0.25_DP*(3.0_DP*Dpoints(1,i)**2 - &
                                         Dpoints(1,i) - 1.0_DP)
      Dbas(4,DER_DERIV1D_X,i) = 0.25_DP*(3.0_DP*Dpoints(1,i)**2 + &
                                         Dpoints(1,i) - 1.0_DP)
    END DO
  ENDIF
    
  END SUBROUTINE 

  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_S31_1D_sim (ieltyp, Dcoords, Djac, Ddetj, &
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
  REAL(DP) :: dxj
    
  ! Clear the output array
  !Dbas = 0.0_DP

  !if function values are desired
  IF (Bder(DER_FUNC1D)) THEN
  
    DO j=1,nelements

      DO i=1,npoints
        Dbas(1,DER_FUNC1D,i,j) = 0.25_DP*Dpoints(1,i,j)*&
                                (Dpoints(1,i,j)**2 - 3.0_DP) + 0.5_DP
        Dbas(2,DER_FUNC1D,i,j) = 0.25_DP*Dpoints(1,i,j)*&
                                (3.0_DP - Dpoints(1,i,j)**2) + 0.5_DP
        Dbas(3,DER_FUNC1D,i,j) = Ddetj(i,j)*0.25*(Dpoints(1,i,j) + 1.0_DP)*&
                                (Dpoints(1,i,j) - 1.0_DP)**2
        Dbas(4,DER_FUNC1D,i,j) = Ddetj(i,j)*0.25*(Dpoints(1,i,j) - 1.0_DP)*&
                                (Dpoints(1,i,j) + 1.0_DP)**2
      END DO
      
    END DO
    
  END IF
    
  !if x-derivatives are desired
  IF (Bder(DER_DERIV1D_X)) THEN
  
    DO j=1,nelements
    
      DO i=1,npoints
        dxj = 1.0_DP / Ddetj(i,j)
        Dbas(1,DER_DERIV1D_X,i,j) = 0.75_DP*(Dpoints(1,i,j)**2 - 1.0_DP)*dxj
        Dbas(2,DER_DERIV1D_X,i,j) = 0.75_DP*(1.0_DP - Dpoints(1,i,j)**2)*dxj
        Dbas(3,DER_DERIV1D_X,i,j) = 0.25_DP*(3.0_DP*Dpoints(1,i,j)**2 - &
                                             Dpoints(1,i,j) - 1.0_DP)
        Dbas(4,DER_DERIV1D_X,i,j) = 0.25_DP*(3.0_DP*Dpoints(1,i,j)**2 + &
                                             Dpoints(1,i,j) - 1.0_DP)
      END DO

    END DO
      
  END IF
    
  END SUBROUTINE 

! ----------------------------------------------------------------------------
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
! -----------------------------------------------------------------------------

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
! Element subroutines for parametric P1~ element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  PURE SUBROUTINE elem_P1T (ieltyp, Dcoords, Djac, ddetj, Bder, &
                            Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_P1T.
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
  
  ! The P1~ space consists of 'linear' finite elements. We have three basis 
  ! functions on the reference element, which can be written down in
  ! standard coordinates (-> P(.)) as well as in barycentric coordinates
  ! (-> p(.)). These are:
  !
  !   p1(xi1,xi2,xi3) = 1 - 2*xi3 =  1 - 2*Y       = P1(X,Y)
  !   p2(xi1,xi2,xi3) = 1 - 2*xi1 = -1 + 2*X + 2*Y = P2(X,Y)
  !   p3(xi1,xi2,xi3) = 1 - 2*xi2 =  1 - 2*X       = P3(X,Y)
  
  ! Clear the output array
  !Dbas = 0.0_DP
    
  ! Remark: The P1~-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  ! If function values are desired, calculate them.
  ! Use the p(.) representation in barycentric coordinates to calculate the
  ! function values.
!  if (el_bder(DER_FUNC)) then
    Dbas(1,DER_FUNC) = 1._DP -2E0_DP*Dpoint(3)
    Dbas(2,DER_FUNC) = 1._DP -2E0_DP*Dpoint(1)
    Dbas(3,DER_FUNC) = 1._DP -2E0_DP*Dpoint(2)
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
      Dbas(1,DER_DERIV_X) =  2E0_DP*Djac(2)*dxj
      Dbas(2,DER_DERIV_X) =  2E0_DP*(Djac(4)-Djac(2))*dxj
      Dbas(3,DER_DERIV_X) = -2E0_DP*Djac(2)*dxj
!    endif
    
    !y-derivatives on current element
!    if (el_bder(DER_DERIV_Y)) then
      Dbas(1,DER_DERIV_Y) = -2E0_DP*Djac(1)*dxj
      Dbas(2,DER_DERIV_Y) = -2E0_DP*(Djac(3)- Djac(1))*dxj
      Dbas(3,DER_DERIV_Y) =  2E0_DP*Djac(3)*dxj
!    endif
!  endif
    
  END SUBROUTINE 

  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_P1T_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                 Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

  !<input>

  ! Element type identifier. Must be =EL_P1T.
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

  ! Remark: The P1~-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  !IF (Bder(DER_FUNC)) THEN
    DO i=1,npoints
      Dbas(1,DER_FUNC,i) = 1._DP -2E0_DP*Dpoints(3,i)
      Dbas(2,DER_FUNC,i) = 1._DP -2E0_DP*Dpoints(1,i)
      Dbas(3,DER_FUNC,i) = 1._DP -2E0_DP*Dpoints(2,i)
    END DO
  !ENDIF
  
  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV_X)) .OR. (Bder(DER_DERIV_Y))) THEN
    dxj = 1E0_DP / Ddetj
    
    !x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      DO i=1,npoints
        Dbas(1,DER_DERIV_X,i) =  2E0_DP*Djac(2,i)*dxj(i)
        Dbas(2,DER_DERIV_X,i) =  2E0_DP*(Djac(4,i)-Djac(2,i))*dxj(i)
        Dbas(3,DER_DERIV_X,i) = -2E0_DP*Djac(2,i)*dxj(i)
!      END DO
!    ENDIF
    
    !y-derivatives on current element
!    IF (Bder(DER_DERIV_Y)) THEN
!      DO i=1,npoints
        Dbas(1,DER_DERIV_Y,i) = -2E0_DP*Djac(1,i)*dxj(i)
        Dbas(2,DER_DERIV_Y,i) = -2E0_DP*(Djac(3,i)- Djac(1,i))*dxj(i)
        Dbas(3,DER_DERIV_Y,i) =  2E0_DP*Djac(3,i)*dxj(i)
      END DO
!    ENDIF
!  ENDIF
    
  END SUBROUTINE

  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_P1T_sim (ieltyp, Dcoords, Djac, Ddetj, &
                                Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_P1T.
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
        Dbas(1,DER_FUNC,i,j) = 1._DP -2E0_DP*Dpoints(3,i,j)
        Dbas(2,DER_FUNC,i,j) = 1._DP -2E0_DP*Dpoints(1,i,j)
        Dbas(3,DER_FUNC,i,j) = 1._DP -2E0_DP*Dpoints(2,i,j)
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
          Dbas(1,DER_DERIV_X,i,j) =  2E0_DP*Djac(2,i,j)*dxj(i)
          Dbas(2,DER_DERIV_X,i,j) =  2E0_DP*(Djac(4,i,j)-Djac(2,i,j))*dxj(i)
          Dbas(3,DER_DERIV_X,i,j) = -2E0_DP*Djac(2,i,j)*dxj(i)
        END DO
!      ENDIF
      
      !y-derivatives on current element
!      IF (Bder(DER_DERIV_Y)) THEN
        DO i=1,npoints
          Dbas(1,DER_DERIV_Y,i,j) = -2E0_DP*Djac(1,i,j)*dxj(i)
          Dbas(2,DER_DERIV_Y,i,j) = -2E0_DP*(Djac(3,i,j)- Djac(1,i,j))*dxj(i)
          Dbas(3,DER_DERIV_Y,i,j) =  2E0_DP*Djac(3,i,j)*dxj(i)
        END DO
!      ENDIF

    END DO
      
  END IF
    
  END SUBROUTINE 
  
! ----------------------------------------------------------------------------
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
! -----------------------------------------------------------------------------

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
    Dbas(3,DER_FUNC) = -0.125E0_DP*( 3.0_DP*(dx**2-dy**2)-4.0_DP*dy-2.0_DP)
    Dbas(4,DER_FUNC) = -0.125E0_DP*(-3.0_DP*(dx**2-dy**2)+4.0_DP*dx-2.0_DP)
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
      Dbas(3,DER_FUNC,i) = -0.125E0_DP* &
          ( 3.0_DP*(Dpoints(1,i)**2-Dpoints(2,i)**2)-4.0_DP*Dpoints(2,i)-2.0_DP)
      Dbas(4,DER_FUNC,i) = -0.125E0_DP* &
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
        Dbas(3,DER_FUNC,i,j) = -0.125E0_DP* &
            ( 3.0_DP*(Dpoints(1,i,j)**2-Dpoints(2,i,j)**2)-4.0_DP*Dpoints(2,i,j)-2.0_DP)
        Dbas(4,DER_FUNC,i,j) = -0.125E0_DP* &
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
! Element subroutines for parametric Q2~ element, midpoint value based.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
!**************************************************************************
 
!<subroutine>  

  PURE SUBROUTINE elem_E035 (ieltyp, Dcoords, Djac, ddetj, Bder, &
                             Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_Q2T.
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
  REAL(DP), DIMENSION(8,NDIM2D) :: Dhelp
  REAL(DP) :: dx,dy

  REAL(DP) :: dxj !auxiliary variable

  REAL(DP), PARAMETER :: Q1 = 2.25_DP
  REAL(DP), PARAMETER :: Q2 = 1.875_DP
  REAL(DP), PARAMETER :: Q3 = 0.75_DP
  REAL(DP), PARAMETER :: Q4 = 0.25_DP
  REAL(DP), PARAMETER :: Q5 = 1.5_DP
  REAL(DP), PARAMETER :: Q6 = 5.625_DP
  REAL(DP), PARAMETER :: Q7 = 4.5_DP
  
  ! The Q1 element is specified by four polynomials on the reference element.
  ! These four polynomials are:
  !
  !  p_1(x,y) = -1/4 - 3/4*y + 3/4*y^2 + 3/4*y*x^2
  !  p_2(x,y) = -1/4 + 3/4*x + 3/4*x^2 - 3/4*x*y^2
  !  p_3(x,y) = ...
  !  p_4(x,y) = 
  !  p_5(x,y) = 
  !  p_6(x,y) = 
  !  p_7(x,y) = 
  !  p_8(x,y) = 
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
    Dbas(1,DER_FUNC) = -Q4 - Q3*dy + Q3*dy**2 + Q3*dx**2*dy
    Dbas(2,DER_FUNC) = -Q4 + Q3*dx + Q3*dx**2 - Q3*dx*dy**2
    Dbas(3,DER_FUNC) = -Q4 + Q3*dy + Q3*dy**2 - Q3*dx**2*dy
    Dbas(4,DER_FUNC) = -Q4 - Q3*dx + Q3*dx**2 + Q3*dx*dy**2
    Dbas(5,DER_FUNC) = -Q3*dx - Q3*dx*dy + Q1*dx*dy**2 + Q2*dx**3*dy - Q2*dx*dy**3
    Dbas(6,DER_FUNC) = -Q3*dy + Q3*dx*dy + Q1*dx**2*dy + Q2*dx**3*dy - Q2*dx*dy**3
    Dbas(7,DER_FUNC) =  Q3*dx - Q3*dx*dy - Q1*dx*dy**2 + Q2*dx**3*dy - Q2*dx*dy**3
    Dbas(8,DER_FUNC) =  Q3*dy + Q3*dx*dy - Q1*dx**2*dy + Q2*dx**3*dy - Q2*dx*dy**3
    
!  endif
  
  ! If x-or y-derivatives are desired, calculate them.
  ! The values of the derivatives are calculated by taking the
  ! derivative of the polynomials and multiplying them with the
  ! inverse of the transformation matrix (in each point) as
  ! stated above.
!  if ((Bder(DER_DERIV_X)) .or. (Bder(DER_DERIV_Y))) then
    dxj = 0.5_DP / ddetj
    
    ! x- and y-derivatives on reference element
    Dhelp(1,1) =  Q5*dx*dy
    Dhelp(2,1) =  Q3 + Q5*dx - Q3*dy**2
    Dhelp(3,1) = -Q5*dx*dy
    Dhelp(4,1) = -Q3 + Q5*dx + Q3*dy**2
    Dhelp(5,1) = -Q3 - Q3*dy + Q1*dy**2 + Q6*dx**2*dy - Q2*dy**3
    Dhelp(6,1) =  Q3*dy + Q7*dx*dy + Q6*dx**2*dy - Q2*dy**3
    Dhelp(7,1) =  Q3 - Q3*dy - Q1*dy**2 + Q6*dx**2*dy - Q2*dy**3
    Dhelp(8,1) =  Q3*dy - Q7*dx*dy + Q6*dx**2*dy - Q2*dy**3
    Dhelp(1,2) = -Q3 + Q5*dy + Q3*dx**2
    Dhelp(2,2) = -Q5*dx*dy
    Dhelp(3,2) =  Q3 + Q5*dy - Q3*dx**2
    Dhelp(4,2) =  Q5*dx*dy
    Dhelp(5,2) = -Q3*dx + Q7*dx*dy + Q2*dx**3 - Q6*dy**2*dx
    Dhelp(6,2) = -Q3 + Q3*dx + Q1*dx**2 + Q2*dx**3 - Q6*dy**2*dx
    Dhelp(7,2) = -Q3*dx - Q7*dx*dy + Q2*dx**3 - Q6*dy**2*dx
    Dhelp(8,2) =  Q3 + Q3*dx - Q1*dx**2 + Q2*dx**3 - Q6*dy**2*dx
      
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

!</subroutine>

  ! Clear the output array
  !Dbas = 0.0_DP
  
  ! Q0 is a single basis function, constant in the element.
  ! The function value of the basis function is =1, the derivatives are all 0!
  DBas(1,DER_FUNC,:,:) = 1.0_DP

  END SUBROUTINE 

! ----------------------------------------------------------------------------
! General information: Function values and derivatives of 
!                      tetrahedral elements
!                      with linear transformation from the reference
!                      to the real element.
!
! The element subroutines return
! - the function value and
! - the X-, Y- and Z-derivatives
! of the basis function in a (cubature) point (x,y,z) on the real mesh!
! The coordinates of a (cubature) point is given
! - as barycentric coordinate tuple (xi1, xi2, xi3, xi4), if the
!   the element is parametric; the actual cubature point is then at
!   (x,y,z) = s(t(xi1,xi2,xi3,xi4))
! - as coordinate pair (x,y,z) on the real element, if the element is
!   nonparametric.
! The mapping  s=(s1,s2,s3):R^3->R^3  is the linear mapping "sigma" from 
! transformation.f90, that maps the reference element T^ to the real
! element T; its shape is of no importance here.
! The transformation t:R^4->R^3 maps the coordinates from the barycenric
! coordinate system on the reference element to the standard space
! (to get the (X,Y,Z) coordinate on the reference element), so
!
!    t(xi1,xi2,xi3,xi4) = xi2*[1,0,0] + xi3*[0,1,0] + xi4*[0,0,1]
!
! The linear transformation s(.) from the reference to the real element
! has the form
!
!    s(X,Y,Z) = c1 + c2*X + c3*Y + c4*Z  =:  (x,y,z)
!
! Let u be an arbitrary FE (basis) function on an element T and
! p be the associated polynomial on the reference element T^. Then,
! we can obtain the function value of u at (x,y,z) easily:
!
!   u(x,y,z) = u(s(t(xi1,xi2,xi3,xi4)) = p(t^-1(s^-1(x,y,z)))
!            = p(xi1,xi2,xi3,xi4)
!
! The derivative is a little bit harder to calculate because of the
! mapping. We have:
!
!    grad(u)(x,y,z) = grad( p(t^-1(s^-1(x,y,z))) )
!
! Let's use the notation 
!
!    P(X,Y,Z)  :=  p( t^-1 (X,Y,Z) )  =  p (xi1,xi2,xi3,xi4)
!
! which is the 'standard' form of the polynomials on the
! reference element without barycentric coordinates involved (i.e.
! for the P1-element it's of the form P(x,y,z) = c1*x + c2*y + c3*z + c4).
!
! Because of the chain rule, we have:
!
!    grad( p( t^-1 (s^-1) ) ) = (DP)( (s^-1) * D(s^-1) )
!
!                                           / (s1^-1)_x (s1^-1)_y (s1^-1)_z \
!       = (P_X(s^-1) P_Y(s^-1) P_Z(s^-1)) * | (s2^-1)_x (s2^-1)_y (s1^-1)_z |
!                                           \ (s3^-1)_x (s3^-1)_y (s3^-1)_z /
!
! With s^-1(x,y,z)=(X,Y,Z), we therefore have:
!                     
!    grad(u)(x,y,z) = 
!
! / P_X(X,Y,Z)*(s1^-1)_x  +  P_Y(X,Y,Z)*(s1^-1)_y  +  P_Z(X,Y,Z)*(s1^-1)_z \
! | P_X(X,Y,Z)*(s2^-1)_x  +  P_Y(X,Y,Z)*(s2^-1)_y  +  P_Z(X,Y,Z)*(s2^-1)_z |
! \ P_X(X,Y,Z)*(s3^-1)_x  +  P_Y(X,Y,Z)*(s3^-1)_y  +  P_Z(X,Y,Z)*(s3^-1)_z /
!
! Now, from e.g. http://mathworld.wolfram.com/MatrixInverse.html we know,
! that:
! 
!         / a b c \                       / ei-fh  ch-bi  bf-ce \
!     A = | d e f |  =>   A^-1 = 1/det(A) | fg-di  ai-cg  cd-af |
!         \ g h i /                       \ dh-eg  bg-ah  ae-bd /
!
! So we have:
! (s1^-1)_x = s2_Y * s3_Z - s2_Z * s3_Y
! (s1^-1)_y = s1_Z * s3_Y - s1_Y * s3_Z
! (s1^-1)_z = s1_Y * s2_Z - s1_Z * s2_Y
! (s2^-1)_x = s2_Z * s3_X - s2_X * s3_Z
! (s2^-1)_y = s1_X * s3_Z - s1_Z * s3_X
! (s2^-1)_z = s1_Z * s2_X - s1_X * s2_Z
! (s3^-1)_x = s2_X * s3_Y - s2_Y * s3_X
! (s3^-1)_y = s1_Y * s3_X - s1_X * s3_Y
! (s3^-1)_z = s1_X * s2_Y - s1_Y * s2_X
!
! being the matrix from the transformation (according to
! http://mathworld.wolfram.com/BarycentricCoordinates.html).
! ----------------------------------------------------------------------------

!**************************************************************************
! Element subroutines for parametric P1_3D element.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  PURE SUBROUTINE elem_P1_3D (ieltyp, Dcoords, Djac, ddetj, Bder, &
                              Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_P1_3D.
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

  REAL(DP) :: dxj !auxiliary variable
  
  ! A matrix for the inverse jacobian matrix
  REAL(DP), DIMENSION(9) :: Dinv
  
  ! The P1 space consists of 'linear' finite elements. We have four basis 
  ! functions on the reference element, which can be written down in
  ! standard coordinates (-> P(.)) as well as in barycentric coordinates
  ! (-> p(.)). These are:
  !
  !   p1(xi1,xi2,xi3,xi4) = xi1 =  1 - X - Y - Z = P1(X,Y,Z)
  !   p2(xi1,xi2,xi3,xi4) = xi2 =  X             = P2(X,Y,Z)
  !   p3(xi1,xi2,xi3,xi4) = xi3 =  Y             = P3(X,Y,Z)
  !   p4(xi1,xi2,xi3,xi4) = xi4 =  Z             = P4(X,Y,Z)
  
  ! Clear the output array
  !Dbas = 0.0_DP
    
  ! Remark: The P1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  ! If function values are desired, calculate them.
  ! Use the p(.) representation in barycentric coordinates to calculate the
  ! function values.
!  if (el_bder(DER_FUNC)) then
    Dbas(1:4,DER_FUNC) = Dpoint(1:4)
!  endif
  
  ! If x-, y- or z-derivatives are desired, calculate them.
  ! Here, we use the P(.) representation to get P_X, P_Y and P_Z (which are
  ! only 0, 1 or -1)!
  ! These are then multiplied with the inverse of the transformation
  ! as described above to get the actual values of the derivatives.
  
!  if ((el_bder(DER_DERIV3D_X)) .or. (el_bder(DER_DERIV3D_Y)) .or. &
!      (el_bder(DER_DERIV3D_Z))) then
    dxj = 1.0_DP / ddetj
    
    ! Invert the jacobian matrix
    Dinv(1)=(Djac(5)*Djac(9)-Djac(8)*Djac(6))*dxj
    Dinv(2)=(Djac(8)*Djac(3)-Djac(2)*Djac(9))*dxj
    Dinv(3)=(Djac(2)*Djac(6)-Djac(5)*Djac(3))*dxj
    Dinv(4)=(Djac(7)*Djac(6)-Djac(4)*Djac(9))*dxj
    Dinv(5)=(Djac(1)*Djac(9)-Djac(7)*Djac(3))*dxj
    Dinv(6)=(Djac(4)*Djac(3)-Djac(1)*Djac(6))*dxj
    Dinv(7)=(Djac(4)*Djac(8)-Djac(7)*Djac(5))*dxj
    Dinv(8)=(Djac(7)*Djac(2)-Djac(1)*Djac(8))*dxj
    Dinv(9)=(Djac(1)*Djac(5)-Djac(4)*Djac(2))*dxj
    
    ! x-derivatives on current element.
!    if (el_bder(DER_DERIV3D_X)) then
      Dbas(1,DER_DERIV3D_X) = -(Dinv(1)+Dinv(2)+Dinv(3))
      Dbas(2,DER_DERIV3D_X) = Dinv(1)
      Dbas(3,DER_DERIV3D_X) = Dinv(2)
      Dbas(4,DER_DERIV3D_X) = Dinv(3)
!    endif
    
    !y-derivatives on current element
!    if (el_bder(DER_DERIV3D_Y)) then
      Dbas(1,DER_DERIV3D_Y) = -(Dinv(4)+Dinv(5)+Dinv(6))
      Dbas(2,DER_DERIV3D_Y) = Dinv(4)
      Dbas(3,DER_DERIV3D_Y) = Dinv(5)
      Dbas(4,DER_DERIV3D_Y) = Dinv(6)
!    endif

    !z-derivatives on current element
!    if (el_bder(DER_DERIV3D_Z)) then
      Dbas(1,DER_DERIV3D_Z) = -(Dinv(7)+Dinv(8)+Dinv(9))
      Dbas(2,DER_DERIV3D_Z) = Dinv(7)
      Dbas(3,DER_DERIV3D_Z) = Dinv(8)
      Dbas(4,DER_DERIV3D_Z) = Dinv(9)
!    endif
!  endif
    
  END SUBROUTINE 
  
  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_P1_3D_mult (ieltyp, Dcoords, Djac, Ddetj, &
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

  REAL(DP) :: dxj ! auxiliary variable
  
  INTEGER :: i   ! point counter

  ! A matrix for the inverse jacobian matrix
  REAL(DP), DIMENSION(9) :: Dinv
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The P1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  !IF (Bder(DER_FUNC3D)) THEN
    DO i=1,npoints
      Dbas(1,DER_FUNC3D,i) = Dpoints(1,i)
      Dbas(2,DER_FUNC3D,i) = Dpoints(2,i)
      Dbas(3,DER_FUNC3D,i) = Dpoints(3,i)
      Dbas(4,DER_FUNC3D,i) = Dpoints(4,i)
    END DO
  !ENDIF
  
  !if x-or y-derivatives are desired
!  IF ((el_bder(DER_DERIV3D_X)) .OR. (el_bder(DER_DERIV3D_Y)) .OR. &
!      (el_bder(DER_DERIV3D_Z))) THEN

    ! Since the jacobian matrix (and therefore also its determinant) is
    ! constant for all points, we only need to invert the matrix once.
    dxj = 1.0_DP / Ddetj(1)
    Dinv(1)=(Djac(5,1)*Djac(9,1)-Djac(8,1)*Djac(6,1))*dxj
    Dinv(2)=(Djac(8,1)*Djac(3,1)-Djac(2,1)*Djac(9,1))*dxj
    Dinv(3)=(Djac(2,1)*Djac(6,1)-Djac(5,1)*Djac(3,1))*dxj
    Dinv(4)=(Djac(7,1)*Djac(6,1)-Djac(4,1)*Djac(9,1))*dxj
    Dinv(5)=(Djac(1,1)*Djac(9,1)-Djac(7,1)*Djac(3,1))*dxj
    Dinv(6)=(Djac(4,1)*Djac(3,1)-Djac(1,1)*Djac(6,1))*dxj
    Dinv(7)=(Djac(4,1)*Djac(8,1)-Djac(7,1)*Djac(5,1))*dxj
    Dinv(8)=(Djac(7,1)*Djac(2,1)-Djac(1,1)*Djac(8,1))*dxj
    Dinv(9)=(Djac(1,1)*Djac(5,1)-Djac(4,1)*Djac(2,1))*dxj

    !x-derivatives on current element
!    IF (Bder(DER_DERIV3D_X)) THEN
      DO i=1,npoints
        Dbas(1,DER_DERIV3D_X,i) = -(Dinv(1)+Dinv(2)+Dinv(3))
        Dbas(2,DER_DERIV3D_X,i) = Dinv(1)
        Dbas(3,DER_DERIV3D_X,i) = Dinv(2)
        Dbas(4,DER_DERIV3D_X,i) = Dinv(3)
!      END DO
!    ENDIF
    
    !y-derivatives on current element
!    IF (Bder(DER_DERIV3D_Y)) THEN
!      DO i=1,npoints
        Dbas(1,DER_DERIV3D_Y,i) = -(Dinv(4)+Dinv(5)+Dinv(6))
        Dbas(2,DER_DERIV3D_Y,i) = Dinv(4)
        Dbas(3,DER_DERIV3D_Y,i) = Dinv(5)
        Dbas(4,DER_DERIV3D_Y,i) = Dinv(6)
!      END DO
!    ENDIF

    !z-derivatives on current element
!    IF (el_bder(DER_DERIV3D_Z)) THEN
!      DO i=1,npoints
        Dbas(1,DER_DERIV3D_Z,i) = -(Dinv(7)+Dinv(8)+Dinv(9))
        Dbas(2,DER_DERIV3D_Z,i) = Dinv(7)
        Dbas(3,DER_DERIV3D_Z,i) = Dinv(8)
        Dbas(4,DER_DERIV3D_Z,i) = Dinv(9)
      END DO
!    ENDIF
!  ENDIF
    
  END SUBROUTINE 

  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_P1_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
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

  REAL(DP) :: dxj !auxiliary variable
  
  INTEGER :: i   ! point counter
  INTEGER :: j   ! element counter

  ! A matrix for the inverse jacobian matrix
  REAL(DP), DIMENSION(9) :: Dinv
    
  ! Clear the output array
  !Dbas = 0.0_DP

  !if function values are desired
  IF (Bder(DER_FUNC3D)) THEN
  
    DO j=1,nelements
    
      DO i=1,npoints
        Dbas(1,DER_FUNC3D,i,j) = Dpoints(1,i,j)
        Dbas(2,DER_FUNC3D,i,j) = Dpoints(2,i,j)
        Dbas(3,DER_FUNC3D,i,j) = Dpoints(3,i,j)
        Dbas(4,DER_FUNC3D,i,j) = Dpoints(4,i,j)
      END DO
      
    END DO
    
  END IF
    
  !if x-or y-derivatives are desired
  IF ((Bder(DER_DERIV3D_X)) .OR. (Bder(DER_DERIV3D_Y)) .OR. &
      (Bder(DER_DERIV3D_Z))) THEN
  
    DO j=1,nelements

      ! Since the jacobian matrix (and therefore also its determinant) is
      ! constant for all points on one element, we only need to invert
      ! the matrix once per element.
      dxj = 1.0_DP / Ddetj(1,j)
      Dinv(1)=(Djac(5,1,j)*Djac(9,1,j)-Djac(8,1,j)*Djac(6,1,j))*dxj
      Dinv(2)=(Djac(8,1,j)*Djac(3,1,j)-Djac(2,1,j)*Djac(9,1,j))*dxj
      Dinv(3)=(Djac(2,1,j)*Djac(6,1,j)-Djac(5,1,j)*Djac(3,1,j))*dxj
      Dinv(4)=(Djac(7,1,j)*Djac(6,1,j)-Djac(4,1,j)*Djac(9,1,j))*dxj
      Dinv(5)=(Djac(1,1,j)*Djac(9,1,j)-Djac(7,1,j)*Djac(3,1,j))*dxj
      Dinv(6)=(Djac(4,1,j)*Djac(3,1,j)-Djac(1,1,j)*Djac(6,1,j))*dxj
      Dinv(7)=(Djac(4,1,j)*Djac(8,1,j)-Djac(7,1,j)*Djac(5,1,j))*dxj
      Dinv(8)=(Djac(7,1,j)*Djac(2,1,j)-Djac(1,1,j)*Djac(8,1,j))*dxj
      Dinv(9)=(Djac(1,1,j)*Djac(5,1,j)-Djac(4,1,j)*Djac(2,1,j))*dxj
     
      !x-derivatives on current element
!      IF (Bder(DER_DERIV3D_X)) THEN
        DO i=1,npoints
          Dbas(1,DER_DERIV3D_X,i,j) = -(Dinv(1)+Dinv(2)+Dinv(3))
          Dbas(2,DER_DERIV3D_X,i,j) = Dinv(1)
          Dbas(3,DER_DERIV3D_X,i,j) = Dinv(2)
          Dbas(4,DER_DERIV3D_X,i,j) = Dinv(3)
!        END DO
!      ENDIF
    
      !y-derivatives on current element
!      IF (Bder(DER_DERIV3D_Y)) THEN
!        DO i=1,npoints
          Dbas(1,DER_DERIV3D_Y,i,j) = -(Dinv(4)+Dinv(5)+Dinv(6))
          Dbas(2,DER_DERIV3D_Y,i,j) = Dinv(4)
          Dbas(3,DER_DERIV3D_Y,i,j) = Dinv(5)
          Dbas(4,DER_DERIV3D_Y,i,j) = Dinv(6)
!        END DO
!      ENDIF

      !z-derivatives on current element
!      IF (el_bder(DER_DERIV3D_Z)) THEN
!        DO i=1,npoints
          Dbas(1,DER_DERIV3D_Z,i,j) = -(Dinv(7)+Dinv(8)+Dinv(9))
          Dbas(2,DER_DERIV3D_Z,i,j) = Dinv(7)
          Dbas(3,DER_DERIV3D_Z,i,j) = Dinv(8)
          Dbas(4,DER_DERIV3D_Z,i,j) = Dinv(9)
        END DO
!      ENDIF
    END DO
      
  END IF
    
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
          (djx*Dhelp(1,1) - djy*Dhelp(1,2) + djz*Dhelp(1,3))
      Dbas(2,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(2,1) - djy*Dhelp(2,2) + djz*Dhelp(2,3))
      Dbas(3,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(3,1) - djy*Dhelp(3,2) + djz*Dhelp(3,3))
      Dbas(4,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(4,1) - djy*Dhelp(4,2) + djz*Dhelp(4,3))
      Dbas(5,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(5,1) - djy*Dhelp(5,2) + djz*Dhelp(5,3))
      Dbas(6,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(6,1) - djy*Dhelp(6,2) + djz*Dhelp(6,3))
      Dbas(7,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(7,1) - djy*Dhelp(7,2) + djz*Dhelp(7,3))
      Dbas(8,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(8,1) - djy*Dhelp(8,2) + djz*Dhelp(8,3))
!    endif
    
    ! y-derivatives on current element
!    if (Bder(DER_DERIV3D_Y)) then
      djx = Djac(7)*Djac(6) - Djac(4)*Djac(9)
      djy = Djac(1)*Djac(9) - Djac(7)*Djac(3)
      djz = Djac(4)*Djac(3) - Djac(1)*Djac(6)
      Dbas(1,DER_DERIV3D_Y) = dxj * &
          (-djx*Dhelp(1,1) + djy*Dhelp(1,2) - djz*Dhelp(1,3))
      Dbas(2,DER_DERIV3D_Y) = dxj * &
          (-djx*Dhelp(2,1) + djy*Dhelp(2,2) - djz*Dhelp(2,3))
      Dbas(3,DER_DERIV3D_Y) = dxj * &
          (-djx*Dhelp(3,1) + djy*Dhelp(3,2) - djz*Dhelp(3,3))
      Dbas(4,DER_DERIV3D_Y) = dxj * &
          (-djx*Dhelp(4,1) + djy*Dhelp(4,2) - djz*Dhelp(4,3))
      Dbas(5,DER_DERIV3D_Y) = dxj * &
          (-djx*Dhelp(5,1) + djy*Dhelp(5,2) - djz*Dhelp(5,3))
      Dbas(6,DER_DERIV3D_Y) = dxj * &
          (-djx*Dhelp(6,1) + djy*Dhelp(6,2) - djz*Dhelp(6,3))
      Dbas(7,DER_DERIV3D_Y) = dxj * &
          (-djx*Dhelp(7,1) + djy*Dhelp(7,2) - djz*Dhelp(7,3))
      Dbas(8,DER_DERIV3D_Y) = dxj * &
          (-djx*Dhelp(8,1) + djy*Dhelp(8,2) - djz*Dhelp(8,3))
!    endif

    ! z-derivatives on current element
!    if (Bder(DER_DERIV3D_Z)) then
      djx = Djac(4)*Djac(8) - Djac(7)*Djac(5)
      djy = Djac(7)*Djac(2) - Djac(1)*Djac(8)
      djz = Djac(1)*Djac(5) - Djac(4)*Djac(2)
      Dbas(1,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(1,1) - djy*Dhelp(1,2) + djz*Dhelp(1,3))
      Dbas(2,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(2,1) - djy*Dhelp(2,2) + djz*Dhelp(2,3))
      Dbas(3,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(3,1) - djy*Dhelp(3,2) + djz*Dhelp(3,3))
      Dbas(4,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(4,1) - djy*Dhelp(4,2) + djz*Dhelp(4,3))
      Dbas(5,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(5,1) - djy*Dhelp(5,2) + djz*Dhelp(5,3))
      Dbas(6,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(6,1) - djy*Dhelp(6,2) + djz*Dhelp(6,3))
      Dbas(7,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(7,1) - djy*Dhelp(7,2) + djz*Dhelp(7,3))
      Dbas(8,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(8,1) - djy*Dhelp(8,2) + djz*Dhelp(8,3))
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
            (djx*Dhelp(1,1,i) - djy*Dhelp(1,2,i) + djz*Dhelp(1,3,i))
        Dbas(2,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(2,1,i) - djy*Dhelp(2,2,i) + djz*Dhelp(2,3,i))
        Dbas(3,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(3,1,i) - djy*Dhelp(3,2,i) + djz*Dhelp(3,3,i))
        Dbas(4,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(4,1,i) - djy*Dhelp(4,2,i) + djz*Dhelp(4,3,i))
        Dbas(5,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(5,1,i) - djy*Dhelp(5,2,i) + djz*Dhelp(5,3,i))
        Dbas(6,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(6,1,i) - djy*Dhelp(6,2,i) + djz*Dhelp(6,3,i))
        Dbas(7,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(7,1,i) - djy*Dhelp(7,2,i) + djz*Dhelp(7,3,i))
        Dbas(8,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(8,1,i) - djy*Dhelp(8,2,i) + djz*Dhelp(8,3,i))
!      END DO
!    ENDIF
    
    !y-derivatives on current element
!    IF (Bder(DER_DERIV_Y)) THEN
!      DO i=1,npoints
        djx = Djac(7,i)*Djac(6,i) - Djac(4,i)*Djac(9,i)
        djy = Djac(1,i)*Djac(9,i) - Djac(7,i)*Djac(3,i)
        djz = Djac(4,i)*Djac(3,i) - Djac(1,i)*Djac(6,i)
        Dbas(1,DER_DERIV3D_Y,i) = Dxj(i) * &
            (-djx*Dhelp(1,1,i) + djy*Dhelp(1,2,i) - djz*Dhelp(1,3,i))
        Dbas(2,DER_DERIV3D_Y,i) = Dxj(i) * &
            (-djx*Dhelp(2,1,i) + djy*Dhelp(2,2,i) - djz*Dhelp(2,3,i))
        Dbas(3,DER_DERIV3D_Y,i) = Dxj(i) * &
            (-djx*Dhelp(3,1,i) + djy*Dhelp(3,2,i) - djz*Dhelp(3,3,i))
        Dbas(4,DER_DERIV3D_Y,i) = Dxj(i) * &
            (-djx*Dhelp(4,1,i) + djy*Dhelp(4,2,i) - djz*Dhelp(4,3,i))
        Dbas(5,DER_DERIV3D_Y,i) = Dxj(i) * &
            (-djx*Dhelp(5,1,i) + djy*Dhelp(5,2,i) - djz*Dhelp(5,3,i))
        Dbas(6,DER_DERIV3D_Y,i) = Dxj(i) * &
            (-djx*Dhelp(6,1,i) + djy*Dhelp(6,2,i) - djz*Dhelp(6,3,i))
        Dbas(7,DER_DERIV3D_Y,i) = Dxj(i) * &
            (-djx*Dhelp(7,1,i) + djy*Dhelp(7,2,i) - djz*Dhelp(7,3,i))
        Dbas(8,DER_DERIV3D_Y,i) = Dxj(i) * &
            (-djx*Dhelp(8,1,i) + djy*Dhelp(8,2,i) - djz*Dhelp(8,3,i))
!      END DO
!    ENDIF

    !z-derivatives on current element
!    IF (Bder(DER_DERIV3D_Z)) THEN
!      DO i=1,npoints
        djx = Djac(4,i)*Djac(8,i) - Djac(7,i)*Djac(5,i)
        djy = Djac(7,i)*Djac(2,i) - Djac(1,i)*Djac(8,i)
        djz = Djac(1,i)*Djac(5,i) - Djac(4,i)*Djac(2,i)
        Dbas(1,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(1,1,i) - djy*Dhelp(1,2,i) + djz*Dhelp(1,3,i))
        Dbas(2,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(2,1,i) - djy*Dhelp(2,2,i) + djz*Dhelp(2,3,i))
        Dbas(3,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(3,1,i) - djy*Dhelp(3,2,i) + djz*Dhelp(3,3,i))
        Dbas(4,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(4,1,i) - djy*Dhelp(4,2,i) + djz*Dhelp(4,3,i))
        Dbas(5,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(5,1,i) - djy*Dhelp(5,2,i) + djz*Dhelp(5,3,i))
        Dbas(6,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(6,1,i) - djy*Dhelp(6,2,i) + djz*Dhelp(6,3,i))
        Dbas(7,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(7,1,i) - djy*Dhelp(7,2,i) + djz*Dhelp(7,3,i))
        Dbas(8,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(8,1,i) - djy*Dhelp(8,2,i) + djz*Dhelp(8,3,i))
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
              (djx*Dhelp(1,1,i) - djy*Dhelp(1,2,i) + djz*Dhelp(1,3,i))
          Dbas(2,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(2,1,i) - djy*Dhelp(2,2,i) + djz*Dhelp(2,3,i))
          Dbas(3,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(3,1,i) - djy*Dhelp(3,2,i) + djz*Dhelp(3,3,i))
          Dbas(4,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(4,1,i) - djy*Dhelp(4,2,i) + djz*Dhelp(4,3,i))
          Dbas(5,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(5,1,i) - djy*Dhelp(5,2,i) + djz*Dhelp(5,3,i))
          Dbas(6,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(6,1,i) - djy*Dhelp(6,2,i) + djz*Dhelp(6,3,i))
          Dbas(7,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(7,1,i) - djy*Dhelp(7,2,i) + djz*Dhelp(7,3,i))
          Dbas(8,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(8,1,i) - djy*Dhelp(8,2,i) + djz*Dhelp(8,3,i))
        END DO
!      ENDIF
      
      !y-derivatives on current element
!      IF (Bder(DER_DERIV3D_Y)) THEN
        DO i=1,npoints
          djx = Djac(7,i,j)*Djac(6,i,j) - Djac(4,i,j)*Djac(9,i,j)
          djy = Djac(1,i,j)*Djac(9,i,j) - Djac(7,i,j)*Djac(3,i,j)
          djz = Djac(4,i,j)*Djac(3,i,j) - Djac(1,i,j)*Djac(6,i,j)
          Dbas(1,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (-djx*Dhelp(1,1,i) + djy*Dhelp(1,2,i) - djz*Dhelp(1,3,i))
          Dbas(2,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (-djx*Dhelp(2,1,i) + djy*Dhelp(2,2,i) - djz*Dhelp(2,3,i))
          Dbas(3,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (-djx*Dhelp(3,1,i) + djy*Dhelp(3,2,i) - djz*Dhelp(3,3,i))
          Dbas(4,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (-djx*Dhelp(4,1,i) + djy*Dhelp(4,2,i) - djz*Dhelp(4,3,i))
          Dbas(5,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (-djx*Dhelp(5,1,i) + djy*Dhelp(5,2,i) - djz*Dhelp(5,3,i))
          Dbas(6,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (-djx*Dhelp(6,1,i) + djy*Dhelp(6,2,i) - djz*Dhelp(6,3,i))
          Dbas(7,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (-djx*Dhelp(7,1,i) + djy*Dhelp(7,2,i) - djz*Dhelp(7,3,i))
          Dbas(8,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (-djx*Dhelp(8,1,i) + djy*Dhelp(8,2,i) - djz*Dhelp(8,3,i))
        END DO
!      ENDIF

      !z-derivatives on current element
!      IF (Bder(DER_DERIV3D_Z)) THEN
        DO i=1,npoints
          djx = Djac(4,i,j)*Djac(8,i,j) - Djac(7,i,j)*Djac(5,i,j)
          djy = Djac(7,i,j)*Djac(2,i,j) - Djac(1,i,j)*Djac(8,i,j)
          djz = Djac(1,i,j)*Djac(5,i,j) - Djac(4,i,j)*Djac(2,i,j)
          Dbas(1,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(1,1,i) - djy*Dhelp(1,2,i) + djz*Dhelp(1,3,i))
          Dbas(2,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(2,1,i) - djy*Dhelp(2,2,i) + djz*Dhelp(2,3,i))
          Dbas(3,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(3,1,i) - djy*Dhelp(3,2,i) + djz*Dhelp(3,3,i))
          Dbas(4,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(4,1,i) - djy*Dhelp(4,2,i) + djz*Dhelp(4,3,i))
          Dbas(5,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(5,1,i) - djy*Dhelp(5,2,i) + djz*Dhelp(5,3,i))
          Dbas(6,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(6,1,i) - djy*Dhelp(6,2,i) + djz*Dhelp(6,3,i))
          Dbas(7,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(7,1,i) - djy*Dhelp(7,2,i) + djz*Dhelp(7,3,i))
          Dbas(8,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(8,1,i) - djy*Dhelp(8,2,i) + djz*Dhelp(8,3,i))
        END DO
!      ENDIF
    END DO
      
  END IF
    
  END SUBROUTINE 

!**************************************************************************
! Element subroutines for parametric Q1~ element, integral mean value based.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  PURE SUBROUTINE elem_E030_3D (ieltyp, Dcoords, Djac, ddetj, Bder, &
                                Dpoint, Dbas)

  !<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  !</description>

  !<input>

  ! Element type identifier. Must be =EL_E030_3D.
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

  ! auxiliary vector containing the first derivatives on the reference element
  REAL(DP), DIMENSION(6,NDIM3D) :: Dhelp
  ! auxiliary variables
  REAL(DP) :: dx,dy,dz,dxy,dyz, djx, djy, djz
  REAL(DP) :: dxj
  REAL(DP), PARAMETER :: R16 = 1.0_DP/6.0_DP

  ! The Q1~ element is specified by six polynomials on the reference element.
  ! These six polynomials are:
  ! P1(x,y,z) = 1/6 - 1/2*z - 1/4*(x^2 - y^2) - 1/2*(y^2 - z^2)
  ! P2(x,y,z) = 1/6 - 1/2*y - 1/4*(x^2 - y^2) + 1/4*(y^2 - z^2)
  ! P3(x,y,z) = 1/6 + 1/2*x + 1/2*(x^2 - y^2) + 1/4*(y^2 - z^2)
  ! P4(x,y,z) = 1/6 + 1/2*y - 1/4*(x^2 - y^2) + 1/4*(y^2 - z^2)
  ! P5(x,y,z) = 1/6 - 1/2*x + 1/2*(x^2 - y^2) + 1/4*(y^2 - z^2)
  ! P6(x,y,z) = 1/6 + 1/2*z - 1/4*(x^2 - y^2) - 1/2*(y^2 - z^2)
  !
  ! Each of them calculated that way that Pi(Xj)=delta_ij (Kronecker)
  ! for X1,...,X6 the integrals over the six faces of the reference element
  ! [-1,1]x[-1,1]x[-1,1].
  
  ! Clear the output array
  !Dbas = 0.0_DP
  
  dx = Dpoint(1)
  dy = Dpoint(2)
  dz = Dpoint(3)
  dxy = dx**2 - dy**2
  dyz = dy**2 - dz**2
    
  ! Remark: The Q1~-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  ! If function values are desired, calculate them.
!  if (el_bder(DER_FUNC3D)) then
    Dbas(1,DER_FUNC3D) = R16 - 0.25_DP*dxy - 0.5_DP*(dz + dyz)
    Dbas(2,DER_FUNC3D) = R16 - 0.25_DP*(dxy - dyz) - 0.5_DP*dy
    Dbas(3,DER_FUNC3D) = R16 + 0.25_DP*dyz + 0.5_DP*(dx + dxy)
    Dbas(4,DER_FUNC3D) = R16 - 0.25_DP*(dxy - dyz) + 0.5_DP*dy
    Dbas(5,DER_FUNC3D) = R16 + 0.25_DP*dyz - 0.5_DP*(dx - dxy)
    Dbas(6,DER_FUNC3D) = R16 - 0.25_DP*dxy + 0.5_DP*(dz - dyz)
!  endif
  
  ! If x-, y- or z-derivatives are desired, calculate them.
  ! The values of the derivatives are calculated by taking the
  ! derivative of the polynomials and multiplying them with the
  ! inverse of the transformation matrix (in each point) as
  ! stated above.
!  if ((Bder(DER_DERIV3D_X)) .or. (Bder(DER_DERIV3D_Y)) .or. &
!      (Bder(DER_DERIV3D_Z))) then
    dxj = 1.0_DP / ddetj
    
    ! x-, y- and z-derivatives on reference element
    djx = -0.5_DP * dx
    Dhelp(1,1) = djx
    Dhelp(2,1) = djx
    Dhelp(3,1) = dx + 0.5_DP
    Dhelp(4,1) = djx
    Dhelp(5,1) = dx - 0.5_DP
    Dhelp(6,1) = djx
    
    djy = -0.5_DP * dy
    Dhelp(1,2) = djy
    Dhelp(2,2) = dy - 0.5_DP
    Dhelp(3,2) = djy
    Dhelp(4,2) = dy + 0.5_DP
    Dhelp(5,2) = djy
    Dhelp(6,2) = djy
    
    djz = -0.5_DP * dz
    Dhelp(1,3) = dz - 0.5_DP
    Dhelp(2,3) = djz
    Dhelp(3,3) = djz
    Dhelp(4,3) = djz
    Dhelp(5,3) = djz
    Dhelp(6,3) = dz + 0.5_DP
      
    ! x-derivatives on current element
!    if (Bder(DER_DERIV3D_X)) then
      djx = Djac(5)*Djac(9) - Djac(6)*Djac(8)
      djy = Djac(8)*Djac(3) - Djac(2)*Djac(9)
      djz = Djac(2)*Djac(6) - Djac(5)*Djac(3)
      Dbas(1,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(1,1) - djy*Dhelp(1,2) + djz*Dhelp(1,3))
      Dbas(2,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(2,1) - djy*Dhelp(2,2) + djz*Dhelp(2,3))
      Dbas(3,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(3,1) - djy*Dhelp(3,2) + djz*Dhelp(3,3))
      Dbas(4,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(4,1) - djy*Dhelp(4,2) + djz*Dhelp(4,3))
      Dbas(5,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(5,1) - djy*Dhelp(5,2) + djz*Dhelp(5,3))
      Dbas(6,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(6,1) - djy*Dhelp(6,2) + djz*Dhelp(6,3))
!    endif
    
    ! y-derivatives on current element
!    if (Bder(DER_DERIV3D_Y)) then
      djx = Djac(7)*Djac(6) - Djac(4)*Djac(9)
      djy = Djac(1)*Djac(9) - Djac(7)*Djac(3)
      djz = Djac(4)*Djac(3) - Djac(1)*Djac(6)
      Dbas(1,DER_DERIV3D_Y) = dxj * &
          (-djx*Dhelp(1,1) + djy*Dhelp(1,2) - djz*Dhelp(1,3))
      Dbas(2,DER_DERIV3D_Y) = dxj * &
          (-djx*Dhelp(2,1) + djy*Dhelp(2,2) - djz*Dhelp(2,3))
      Dbas(3,DER_DERIV3D_Y) = dxj * &
          (-djx*Dhelp(3,1) + djy*Dhelp(3,2) - djz*Dhelp(3,3))
      Dbas(4,DER_DERIV3D_Y) = dxj * &
          (-djx*Dhelp(4,1) + djy*Dhelp(4,2) - djz*Dhelp(4,3))
      Dbas(5,DER_DERIV3D_Y) = dxj * &
          (-djx*Dhelp(5,1) + djy*Dhelp(5,2) - djz*Dhelp(5,3))
      Dbas(6,DER_DERIV3D_Y) = dxj * &
          (-djx*Dhelp(6,1) + djy*Dhelp(6,2) - djz*Dhelp(6,3))
!    endif

    ! z-derivatives on current element
!    if (Bder(DER_DERIV3D_Z)) then
      djx = Djac(4)*Djac(8) - Djac(7)*Djac(5)
      djy = Djac(7)*Djac(2) - Djac(1)*Djac(8)
      djz = Djac(1)*Djac(5) - Djac(4)*Djac(2)
      Dbas(1,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(1,1) - djy*Dhelp(1,2) + djz*Dhelp(1,3))
      Dbas(2,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(2,1) - djy*Dhelp(2,2) + djz*Dhelp(2,3))
      Dbas(3,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(3,1) - djy*Dhelp(3,2) + djz*Dhelp(3,3))
      Dbas(4,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(4,1) - djy*Dhelp(4,2) + djz*Dhelp(4,3))
      Dbas(5,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(5,1) - djy*Dhelp(5,2) + djz*Dhelp(5,3))
      Dbas(6,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(6,1) - djy*Dhelp(6,2) + djz*Dhelp(6,3))
!    endif
!  endif
    
  END SUBROUTINE 

  
  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_E030_3D_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                     Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1T_3D.
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
  REAL(DP), DIMENSION(6,NDIM3D,npoints) :: Dhelp
  ! auxiliary variables
  REAL(DP) :: dx,dy,dz,dxy,dyz, djx, djy, djz
  REAL(DP), DIMENSION(npoints) :: Dxj
  REAL(DP), PARAMETER :: R16 = 1.0_DP/6.0_DP
  INTEGER :: i   ! point counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  !IF (Bder(DER_FUNC3D)) THEN
    DO i=1,npoints
      dx = Dpoints(1,i)
      dy = Dpoints(2,i)
      dz = Dpoints(3,i)
      dxy = dx**2 - dy**2
      dyz = dy**2 - dz**2
      Dbas(1,DER_FUNC3D,i) = R16 - 0.25_DP*dxy - 0.5_DP*(dz + dyz)
      Dbas(2,DER_FUNC3D,i) = R16 - 0.25_DP*(dxy - dyz) - 0.5_DP*dy
      Dbas(3,DER_FUNC3D,i) = R16 + 0.25_DP*dyz + 0.5_DP*(dx + dxy)
      Dbas(4,DER_FUNC3D,i) = R16 - 0.25_DP*(dxy - dyz) + 0.5_DP*dy
      Dbas(5,DER_FUNC3D,i) = R16 + 0.25_DP*dyz - 0.5_DP*(dx - dxy)
      Dbas(6,DER_FUNC3D,i) = R16 - 0.25_DP*dxy + 0.5_DP*(dz - dyz)
    END DO
  !ENDIF
  
  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV3D_X)) .OR. (Bder(DER_DERIV3D_Y)) .OR.&
!      (Bder(DER_DERIV3D_Z))) THEN
    Dxj = 1.0_DP / Ddetj
    
    !x-, y- and z-derivatives on reference element
    DO i=1,npoints
      dx = Dpoints(1,i)
      dy = Dpoints(2,i)
      dz = Dpoints(3,i)
      djx = -0.5_DP * dx
      Dhelp(1,1,i) = djx
      Dhelp(2,1,i) = djx
      Dhelp(3,1,i) = dx + 0.5_DP
      Dhelp(4,1,i) = djx
      Dhelp(5,1,i) = dx - 0.5_DP
      Dhelp(6,1,i) = djx
      djy = -0.5_DP * dy
      Dhelp(1,2,i) = djy
      Dhelp(2,2,i) = dy - 0.5_DP
      Dhelp(3,2,i) = djy
      Dhelp(4,2,i) = dy + 0.5_DP
      Dhelp(5,2,i) = djy
      Dhelp(6,2,i) = djy
      djz = -0.5_DP * dz
      Dhelp(1,3,i) = dz - 0.5_DP
      Dhelp(2,3,i) = djz
      Dhelp(3,3,i) = djz
      Dhelp(4,3,i) = djz
      Dhelp(5,3,i) = djz
      Dhelp(6,3,i) = dz + 0.5_DP
    END DO
      
    !x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      DO i=1,npoints
        djx = Djac(5,i)*Djac(9,i) - Djac(6,i)*Djac(8,i)
        djy = Djac(8,i)*Djac(3,i) - Djac(2,i)*Djac(9,i)
        djz = Djac(2,i)*Djac(6,i) - Djac(5,i)*Djac(3,i)
        Dbas(1,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(1,1,i) - djy*Dhelp(1,2,i) + djz*Dhelp(1,3,i))
        Dbas(2,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(2,1,i) - djy*Dhelp(2,2,i) + djz*Dhelp(2,3,i))
        Dbas(3,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(3,1,i) - djy*Dhelp(3,2,i) + djz*Dhelp(3,3,i))
        Dbas(4,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(4,1,i) - djy*Dhelp(4,2,i) + djz*Dhelp(4,3,i))
        Dbas(5,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(5,1,i) - djy*Dhelp(5,2,i) + djz*Dhelp(5,3,i))
        Dbas(6,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(6,1,i) - djy*Dhelp(6,2,i) + djz*Dhelp(6,3,i))
!      END DO
!    ENDIF
    
    !y-derivatives on current element
!    IF (Bder(DER_DERIV_Y)) THEN
!      DO i=1,npoints
        djx = Djac(7,i)*Djac(6,i) - Djac(4,i)*Djac(9,i)
        djy = Djac(1,i)*Djac(9,i) - Djac(7,i)*Djac(3,i)
        djz = Djac(4,i)*Djac(3,i) - Djac(1,i)*Djac(6,i)
        Dbas(1,DER_DERIV3D_Y,i) = Dxj(i) * &
            (-djx*Dhelp(1,1,i) + djy*Dhelp(1,2,i) - djz*Dhelp(1,3,i))
        Dbas(2,DER_DERIV3D_Y,i) = Dxj(i) * &
            (-djx*Dhelp(2,1,i) + djy*Dhelp(2,2,i) - djz*Dhelp(2,3,i))
        Dbas(3,DER_DERIV3D_Y,i) = Dxj(i) * &
            (-djx*Dhelp(3,1,i) + djy*Dhelp(3,2,i) - djz*Dhelp(3,3,i))
        Dbas(4,DER_DERIV3D_Y,i) = Dxj(i) * &
            (-djx*Dhelp(4,1,i) + djy*Dhelp(4,2,i) - djz*Dhelp(4,3,i))
        Dbas(5,DER_DERIV3D_Y,i) = Dxj(i) * &
            (-djx*Dhelp(5,1,i) + djy*Dhelp(5,2,i) - djz*Dhelp(5,3,i))
        Dbas(6,DER_DERIV3D_Y,i) = Dxj(i) * &
            (-djx*Dhelp(6,1,i) + djy*Dhelp(6,2,i) - djz*Dhelp(6,3,i))
!      END DO
!    ENDIF

    !z-derivatives on current element
!    IF (Bder(DER_DERIV3D_Z)) THEN
!      DO i=1,npoints
        djx = Djac(4,i)*Djac(8,i) - Djac(7,i)*Djac(5,i)
        djy = Djac(7,i)*Djac(2,i) - Djac(1,i)*Djac(8,i)
        djz = Djac(1,i)*Djac(5,i) - Djac(4,i)*Djac(2,i)
        Dbas(1,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(1,1,i) - djy*Dhelp(1,2,i) + djz*Dhelp(1,3,i))
        Dbas(2,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(2,1,i) - djy*Dhelp(2,2,i) + djz*Dhelp(2,3,i))
        Dbas(3,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(3,1,i) - djy*Dhelp(3,2,i) + djz*Dhelp(3,3,i))
        Dbas(4,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(4,1,i) - djy*Dhelp(4,2,i) + djz*Dhelp(4,3,i))
        Dbas(5,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(5,1,i) - djy*Dhelp(5,2,i) + djz*Dhelp(5,3,i))
        Dbas(6,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(6,1,i) - djy*Dhelp(6,2,i) + djz*Dhelp(6,3,i))
      END DO
!    ENDIF
!  ENDIF
    
  END SUBROUTINE 

  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_E030_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
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
  REAL(DP), DIMENSION(6,NDIM3D,npoints) :: Dhelp
  ! auxiliary variables
  REAL(DP) :: dx,dy,dz,dxy,dyz, djx, djy, djz
  REAL(DP),DIMENSION(npoints) :: Dxj
  REAL(DP), PARAMETER :: R16 = 1.0_DP/6.0_DP
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
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)
        dz = Dpoints(3,i,j)
        dxy = dx**2 - dy**2
        dyz = dy**2 - dz**2
        Dbas(1,DER_FUNC3D,i,j) = R16 - 0.25_DP*dxy - 0.5_DP*(dz + dyz)
        Dbas(2,DER_FUNC3D,i,j) = R16 - 0.25_DP*(dxy - dyz) - 0.5_DP*dy
        Dbas(3,DER_FUNC3D,i,j) = R16 + 0.25_DP*dyz + 0.5_DP*(dx + dxy)
        Dbas(4,DER_FUNC3D,i,j) = R16 - 0.25_DP*(dxy - dyz) + 0.5_DP*dy
        Dbas(5,DER_FUNC3D,i,j) = R16 + 0.25_DP*dyz - 0.5_DP*(dx - dxy)
        Dbas(6,DER_FUNC3D,i,j) = R16 - 0.25_DP*dxy + 0.5_DP*(dz - dyz)
      END DO
      
    END DO
    
  END IF
    
  !if x-, y- or z-derivatives are desired
  IF ((Bder(DER_DERIV3D_X)) .OR. (Bder(DER_DERIV3D_Y)) .OR. &
      (Bder(DER_DERIV3D_Z))) THEN
  
    DO j=1,nelements
      Dxj = 1.0_DP / Ddetj(:,j)
      
      !x-, y- and z-derivatives on reference element
      DO i=1,npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)
        dz = Dpoints(3,i,j)
        djx = -0.5_DP * dx
        Dhelp(1,1,i) = djx
        Dhelp(2,1,i) = djx
        Dhelp(3,1,i) = dx + 0.5_DP
        Dhelp(4,1,i) = djx
        Dhelp(5,1,i) = dx - 0.5_DP
        Dhelp(6,1,i) = djx
        djy = -0.5_DP * dy
        Dhelp(1,2,i) = djy
        Dhelp(2,2,i) = dy - 0.5_DP
        Dhelp(3,2,i) = djy
        Dhelp(4,2,i) = dy + 0.5_DP
        Dhelp(5,2,i) = djy
        Dhelp(6,2,i) = djy
        djz = -0.5_DP * dz
        Dhelp(1,3,i) = dz - 0.5_DP
        Dhelp(2,3,i) = djz
        Dhelp(3,3,i) = djz
        Dhelp(4,3,i) = djz
        Dhelp(5,3,i) = djz
        Dhelp(6,3,i) = dz + 0.5_DP
      END DO
        
      !x-derivatives on current element
!      IF (Bder(DER_DERIV3D_X)) THEN
        DO i=1,npoints
          djx = Djac(5,i,j)*Djac(9,i,j) - Djac(6,i,j)*Djac(8,i,j)
          djy = Djac(8,i,j)*Djac(3,i,j) - Djac(2,i,j)*Djac(9,i,j)
          djz = Djac(2,i,j)*Djac(6,i,j) - Djac(5,i,j)*Djac(3,i,j)
          Dbas(1,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(1,1,i) - djy*Dhelp(1,2,i) + djz*Dhelp(1,3,i))
          Dbas(2,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(2,1,i) - djy*Dhelp(2,2,i) + djz*Dhelp(2,3,i))
          Dbas(3,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(3,1,i) - djy*Dhelp(3,2,i) + djz*Dhelp(3,3,i))
          Dbas(4,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(4,1,i) - djy*Dhelp(4,2,i) + djz*Dhelp(4,3,i))
          Dbas(5,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(5,1,i) - djy*Dhelp(5,2,i) + djz*Dhelp(5,3,i))
          Dbas(6,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(6,1,i) - djy*Dhelp(6,2,i) + djz*Dhelp(6,3,i))
        END DO
!      ENDIF
      
      !y-derivatives on current element
!      IF (Bder(DER_DERIV3D_Y)) THEN
        DO i=1,npoints
          djx = Djac(7,i,j)*Djac(6,i,j) - Djac(4,i,j)*Djac(9,i,j)
          djy = Djac(1,i,j)*Djac(9,i,j) - Djac(7,i,j)*Djac(3,i,j)
          djz = Djac(4,i,j)*Djac(3,i,j) - Djac(1,i,j)*Djac(6,i,j)
          Dbas(1,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (-djx*Dhelp(1,1,i) + djy*Dhelp(1,2,i) - djz*Dhelp(1,3,i))
          Dbas(2,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (-djx*Dhelp(2,1,i) + djy*Dhelp(2,2,i) - djz*Dhelp(2,3,i))
          Dbas(3,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (-djx*Dhelp(3,1,i) + djy*Dhelp(3,2,i) - djz*Dhelp(3,3,i))
          Dbas(4,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (-djx*Dhelp(4,1,i) + djy*Dhelp(4,2,i) - djz*Dhelp(4,3,i))
          Dbas(5,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (-djx*Dhelp(5,1,i) + djy*Dhelp(5,2,i) - djz*Dhelp(5,3,i))
          Dbas(6,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (-djx*Dhelp(6,1,i) + djy*Dhelp(6,2,i) - djz*Dhelp(6,3,i))
        END DO
!      ENDIF

      !z-derivatives on current element
!      IF (Bder(DER_DERIV3D_Z)) THEN
        DO i=1,npoints
          djx = Djac(4,i,j)*Djac(8,i,j) - Djac(7,i,j)*Djac(5,i,j)
          djy = Djac(7,i,j)*Djac(2,i,j) - Djac(1,i,j)*Djac(8,i,j)
          djz = Djac(1,i,j)*Djac(5,i,j) - Djac(4,i,j)*Djac(2,i,j)
          Dbas(1,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(1,1,i) - djy*Dhelp(1,2,i) + djz*Dhelp(1,3,i))
          Dbas(2,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(2,1,i) - djy*Dhelp(2,2,i) + djz*Dhelp(2,3,i))
          Dbas(3,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(3,1,i) - djy*Dhelp(3,2,i) + djz*Dhelp(3,3,i))
          Dbas(4,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(4,1,i) - djy*Dhelp(4,2,i) + djz*Dhelp(4,3,i))
          Dbas(5,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(5,1,i) - djy*Dhelp(5,2,i) + djz*Dhelp(5,3,i))
          Dbas(6,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(6,1,i) - djy*Dhelp(6,2,i) + djz*Dhelp(6,3,i))
        END DO
!      ENDIF
    END DO
      
  END IF
    
  END SUBROUTINE 

!**************************************************************************
! Element subroutines for parametric Q1~ element, face-midpoint-value based.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
 
!<subroutine>  

  PURE SUBROUTINE elem_E031_3D (ieltyp, Dcoords, Djac, ddetj, Bder, &
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
  REAL(DP), DIMENSION(6,NDIM3D) :: Dhelp
  REAL(DP) :: dx,dy,dz,dxy,dyz, djx, djy, djz
  REAL(DP), PARAMETER :: R12 = 0.5_DP
  REAL(DP), PARAMETER :: R13 = 1.0_DP/3.0_DP
  REAL(DP), PARAMETER :: R16 = 1.0_DP/6.0_DP

  REAL(DP) :: dxj !auxiliary variable

  ! The Q1~ element is specified by six polynomials on the reference element.
  ! These six polynomials are:
  ! P1(x,y,z) = 1/6 - 1/2*z - 1/6*(x^2 - y^2) - 1/3*(y^2 - z^2)
  ! P2(x,y,z) = 1/6 - 1/2*y - 1/6*(x^2 - y^2) + 1/6*(y^2 - z^2)
  ! P3(x,y,z) = 1/6 + 1/2*x + 1/3*(x^2 - y^2) + 1/6*(y^2 - z^2)
  ! P4(x,y,z) = 1/6 + 1/2*y - 1/6*(x^2 - y^2) + 1/6*(y^2 - z^2)
  ! P5(x,y,z) = 1/6 - 1/2*x + 1/3*(x^2 - y^2) + 1/6*(y^2 - z^2)
  ! P6(x,y,z) = 1/6 + 1/2*z - 1/6*(x^2 - y^2) - 1/3*(y^2 - z^2)
  !
  ! Each of them calculated that way that Pi(Xj)=delta_ij (Kronecker)
  ! for X1,...,X6 the six face midpoints of the reference element
  ! [-1,1]x[-1,1]x[-1,1].
  
  ! Clear the output array
  !Dbas = 0.0_DP
  
  dx = Dpoint(1)
  dy = Dpoint(2)
  dz = Dpoint(3)
  dxy = dx**2 - dy**2
  dyz = dy**2 - dz**2
    
  ! Remark: The Q1~-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  ! If function values are desired, calculate them.
!  if (el_bder(DER_FUNC3D)) then
    Dbas(1,DER_FUNC3D) = R16 - R12*dz - R16*dxy - R13*dyz
    Dbas(2,DER_FUNC3D) = R16 - R12*dy - R16*dxy + R16*dyz
    Dbas(3,DER_FUNC3D) = R16 + R12*dx + R13*dxy + R16*dyz
    Dbas(4,DER_FUNC3D) = R16 + R12*dy - R16*dxy + R16*dyz
    Dbas(5,DER_FUNC3D) = R16 - R12*dx + R13*dxy + R16*dyz
    Dbas(6,DER_FUNC3D) = R16 + R12*dz - R16*dxy - R13*dyz
!  endif
  
  ! If x-, y- or z-derivatives are desired, calculate them.
  ! The values of the derivatives are calculated by taking the
  ! derivative of the polynomials and multiplying them with the
  ! inverse of the transformation matrix (in each point) as
  ! stated above.
!  if ((Bder(DER_DERIV3D_X)) .or. (Bder(DER_DERIV3D_Y)) .or. &
!      (Bder(DER_DERIV3D_Z))) then
    dxj = 1.0_DP / ddetj
    
    ! x-, y- and z-derivatives on reference element
    djx = -R13*dx
    Dhelp(1,1) = djx
    Dhelp(2,1) = djx
    Dhelp(3,1) = R16*(4.0_DP*dx + 3.0_DP)
    Dhelp(4,1) = djx
    Dhelp(5,1) = R16*(4.0_DP*dx - 3.0_DP)
    Dhelp(6,1) = djx
    
    djy = -R13*dy
    Dhelp(1,2) = djy
    Dhelp(2,2) = R16*(4.0_DP*dy - 3.0_DP)
    Dhelp(3,2) = djy
    Dhelp(4,2) = R16*(4.0_DP*dy + 3.0_DP)
    Dhelp(5,2) = djy
    Dhelp(6,2) = djy
    
    djz = -R13*dz
    Dhelp(1,3) = R16*(4.0_DP*dz - 3.0_DP)
    Dhelp(2,3) = djz
    Dhelp(3,3) = djz
    Dhelp(4,3) = djz
    Dhelp(5,3) = djz
    Dhelp(6,3) = R16*(4.0_DP*dz + 3.0_DP)
      
    ! x-derivatives on current element
!    if (Bder(DER_DERIV3D_X)) then
      djx = Djac(5)*Djac(9) - Djac(6)*Djac(8)
      djy = Djac(8)*Djac(3) - Djac(2)*Djac(9)
      djz = Djac(2)*Djac(6) - Djac(5)*Djac(3)
      Dbas(1,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(1,1) - djy*Dhelp(1,2) + djz*Dhelp(1,3))
      Dbas(2,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(2,1) - djy*Dhelp(2,2) + djz*Dhelp(2,3))
      Dbas(3,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(3,1) - djy*Dhelp(3,2) + djz*Dhelp(3,3))
      Dbas(4,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(4,1) - djy*Dhelp(4,2) + djz*Dhelp(4,3))
      Dbas(5,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(5,1) - djy*Dhelp(5,2) + djz*Dhelp(5,3))
      Dbas(6,DER_DERIV3D_X) = dxj * &
          (djx*Dhelp(6,1) - djy*Dhelp(6,2) + djz*Dhelp(6,3))
!    endif
    
    ! y-derivatives on current element
!    if (Bder(DER_DERIV3D_Y)) then
      djx = Djac(7)*Djac(6) - Djac(4)*Djac(9)
      djy = Djac(1)*Djac(9) - Djac(7)*Djac(3)
      djz = Djac(4)*Djac(3) - Djac(1)*Djac(6)
      Dbas(1,DER_DERIV3D_Y) = dxj * &
          (-djx*Dhelp(1,1) + djy*Dhelp(1,2) - djz*Dhelp(1,3))
      Dbas(2,DER_DERIV3D_Y) = dxj * &
          (-djx*Dhelp(2,1) + djy*Dhelp(2,2) - djz*Dhelp(2,3))
      Dbas(3,DER_DERIV3D_Y) = dxj * &
          (-djx*Dhelp(3,1) + djy*Dhelp(3,2) - djz*Dhelp(3,3))
      Dbas(4,DER_DERIV3D_Y) = dxj * &
          (-djx*Dhelp(4,1) + djy*Dhelp(4,2) - djz*Dhelp(4,3))
      Dbas(5,DER_DERIV3D_Y) = dxj * &
          (-djx*Dhelp(5,1) + djy*Dhelp(5,2) - djz*Dhelp(5,3))
      Dbas(6,DER_DERIV3D_Y) = dxj * &
          (-djx*Dhelp(6,1) + djy*Dhelp(6,2) - djz*Dhelp(6,3))
!    endif

    ! z-derivatives on current element
!    if (Bder(DER_DERIV3D_Z)) then
      djx = Djac(4)*Djac(8) - Djac(7)*Djac(5)
      djy = Djac(7)*Djac(2) - Djac(1)*Djac(8)
      djz = Djac(1)*Djac(5) - Djac(4)*Djac(2)
      Dbas(1,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(1,1) - djy*Dhelp(1,2) + djz*Dhelp(1,3))
      Dbas(2,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(2,1) - djy*Dhelp(2,2) + djz*Dhelp(2,3))
      Dbas(3,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(3,1) - djy*Dhelp(3,2) + djz*Dhelp(3,3))
      Dbas(4,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(4,1) - djy*Dhelp(4,2) + djz*Dhelp(4,3))
      Dbas(5,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(5,1) - djy*Dhelp(5,2) + djz*Dhelp(5,3))
      Dbas(6,DER_DERIV3D_Z) = dxj * &
          (djx*Dhelp(6,1) - djy*Dhelp(6,2) + djz*Dhelp(6,3))
!    endif
!  endif
    
  END SUBROUTINE 

  
  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_E031_3D_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                     Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier. Must be =EL_Q1T_3D.
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
  REAL(DP), DIMENSION(6,NDIM3D,npoints) :: Dhelp

  REAL(DP),DIMENSION(npoints) :: Dxj !auxiliary variable
  REAL(DP) :: dx,dy,dz,dxy,dyz, djx, djy, djz
  REAL(DP), PARAMETER :: R12 = 0.5_DP
  REAL(DP), PARAMETER :: R13 = 1.0_DP/3.0_DP
  REAL(DP), PARAMETER :: R16 = 1.0_DP/6.0_DP
  INTEGER :: i   ! point counter
    
  ! Clear the output array
  !Dbas = 0.0_DP

  ! Remark: The Q1-element always computes function value and 1st derivatives.
  ! That's even faster than when using three IF commands for preventing
  ! the computation of one of the values!
      
  !if function values are desired
  !IF (Bder(DER_FUNC3D)) THEN
    DO i=1,npoints
      dx = Dpoints(1,i)
      dy = Dpoints(2,i)
      dz = Dpoints(3,i)
      dxy = dx**2 - dy**2
      dyz = dy**2 - dz**2
      Dbas(1,DER_FUNC3D,i) = R16 - R12*dz - R16*dxy - R13*dyz
      Dbas(2,DER_FUNC3D,i) = R16 - R12*dy - R16*dxy + R16*dyz
      Dbas(3,DER_FUNC3D,i) = R16 + R12*dx + R13*dxy + R16*dyz
      Dbas(4,DER_FUNC3D,i) = R16 + R12*dy - R16*dxy + R16*dyz
      Dbas(5,DER_FUNC3D,i) = R16 - R12*dx + R13*dxy + R16*dyz
      Dbas(6,DER_FUNC3D,i) = R16 + R12*dz - R16*dxy - R13*dyz
    END DO
  !ENDIF
  
  !if x-or y-derivatives are desired
!  IF ((Bder(DER_DERIV3D_X)) .OR. (Bder(DER_DERIV3D_Y)) .OR.&
!      (Bder(DER_DERIV3D_Z))) THEN
    Dxj = 1.0_DP / Ddetj
    
    !x-, y- and z-derivatives on reference element
    DO i=1,npoints
      dx = Dpoints(1,i)
      dy = Dpoints(2,i)
      dz = Dpoints(3,i)
      djx = -R13*dx
      Dhelp(1,1,i) = djx
      Dhelp(2,1,i) = djx
      Dhelp(3,1,i) = R16*(4.0_DP*dx + 3.0_DP)
      Dhelp(4,1,i) = djx
      Dhelp(5,1,i) = R16*(4.0_DP*dx - 3.0_DP)
      Dhelp(6,1,i) = djx
      djy = -R13*dy
      Dhelp(1,2,i) = djy
      Dhelp(2,2,i) = R16*(4.0_DP*dy - 3.0_DP)
      Dhelp(3,2,i) = djy
      Dhelp(4,2,i) = R16*(4.0_DP*dy + 3.0_DP)
      Dhelp(5,2,i) = djy
      Dhelp(6,2,i) = djy
      djz = -R13*dz
      Dhelp(1,3,i) = R16*(4.0_DP*dz - 3.0_DP)
      Dhelp(2,3,i) = djz
      Dhelp(3,3,i) = djz
      Dhelp(4,3,i) = djz
      Dhelp(5,3,i) = djz
      Dhelp(6,3,i) = R16*(4.0_DP*dz + 3.0_DP)
    END DO
      
    !x-derivatives on current element
!    IF (Bder(DER_DERIV_X)) THEN
      DO i=1,npoints
        djx = Djac(5,i)*Djac(9,i) - Djac(6,i)*Djac(8,i)
        djy = Djac(8,i)*Djac(3,i) - Djac(2,i)*Djac(9,i)
        djz = Djac(2,i)*Djac(6,i) - Djac(5,i)*Djac(3,i)
        Dbas(1,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(1,1,i) - djy*Dhelp(1,2,i) + djz*Dhelp(1,3,i))
        Dbas(2,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(2,1,i) - djy*Dhelp(2,2,i) + djz*Dhelp(2,3,i))
        Dbas(3,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(3,1,i) - djy*Dhelp(3,2,i) + djz*Dhelp(3,3,i))
        Dbas(4,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(4,1,i) - djy*Dhelp(4,2,i) + djz*Dhelp(4,3,i))
        Dbas(5,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(5,1,i) - djy*Dhelp(5,2,i) + djz*Dhelp(5,3,i))
        Dbas(6,DER_DERIV3D_X,i) = Dxj(i) * &
            (djx*Dhelp(6,1,i) - djy*Dhelp(6,2,i) + djz*Dhelp(6,3,i))
!      END DO
!    ENDIF
    
    !y-derivatives on current element
!    IF (Bder(DER_DERIV_Y)) THEN
!      DO i=1,npoints
        djx = Djac(7,i)*Djac(6,i) - Djac(4,i)*Djac(9,i)
        djy = Djac(1,i)*Djac(9,i) - Djac(7,i)*Djac(3,i)
        djz = Djac(4,i)*Djac(3,i) - Djac(1,i)*Djac(6,i)
        Dbas(1,DER_DERIV3D_Y,i) = Dxj(i) * &
            (-djx*Dhelp(1,1,i) + djy*Dhelp(1,2,i) - djz*Dhelp(1,3,i))
        Dbas(2,DER_DERIV3D_Y,i) = Dxj(i) * &
            (-djx*Dhelp(2,1,i) + djy*Dhelp(2,2,i) - djz*Dhelp(2,3,i))
        Dbas(3,DER_DERIV3D_Y,i) = Dxj(i) * &
            (-djx*Dhelp(3,1,i) + djy*Dhelp(3,2,i) - djz*Dhelp(3,3,i))
        Dbas(4,DER_DERIV3D_Y,i) = Dxj(i) * &
            (-djx*Dhelp(4,1,i) + djy*Dhelp(4,2,i) - djz*Dhelp(4,3,i))
        Dbas(5,DER_DERIV3D_Y,i) = Dxj(i) * &
            (-djx*Dhelp(5,1,i) + djy*Dhelp(5,2,i) - djz*Dhelp(5,3,i))
        Dbas(6,DER_DERIV3D_Y,i) = Dxj(i) * &
            (-djx*Dhelp(6,1,i) + djy*Dhelp(6,2,i) - djz*Dhelp(6,3,i))
!      END DO
!    ENDIF

    !z-derivatives on current element
!    IF (Bder(DER_DERIV3D_Z)) THEN
!      DO i=1,npoints
        djx = Djac(4,i)*Djac(8,i) - Djac(7,i)*Djac(5,i)
        djy = Djac(7,i)*Djac(2,i) - Djac(1,i)*Djac(8,i)
        djz = Djac(1,i)*Djac(5,i) - Djac(4,i)*Djac(2,i)
        Dbas(1,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(1,1,i) - djy*Dhelp(1,2,i) + djz*Dhelp(1,3,i))
        Dbas(2,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(2,1,i) - djy*Dhelp(2,2,i) + djz*Dhelp(2,3,i))
        Dbas(3,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(3,1,i) - djy*Dhelp(3,2,i) + djz*Dhelp(3,3,i))
        Dbas(4,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(4,1,i) - djy*Dhelp(4,2,i) + djz*Dhelp(4,3,i))
        Dbas(5,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(5,1,i) - djy*Dhelp(5,2,i) + djz*Dhelp(5,3,i))
        Dbas(6,DER_DERIV3D_Z,i) = Dxj(i) * &
            (djx*Dhelp(6,1,i) - djy*Dhelp(6,2,i) + djz*Dhelp(6,3,i))
      END DO
!    ENDIF
!  ENDIF
    
  END SUBROUTINE 

  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_E031_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
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
  REAL(DP), DIMENSION(6,NDIM3D,npoints) :: Dhelp
  REAL(DP) :: dx,dy,dz,dxy,dyz, djx, djy, djz
  REAL(DP), PARAMETER :: R12 = 0.5_DP
  REAL(DP), PARAMETER :: R13 = 1.0_DP/3.0_DP
  REAL(DP), PARAMETER :: R16 = 1.0_DP/6.0_DP
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
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)
        dz = Dpoints(3,i,j)
        dxy = dx**2 - dy**2
        dyz = dy**2 - dz**2
        Dbas(1,DER_FUNC3D,i,j) = R16 - R12*dz - R16*dxy - R13*dyz
        Dbas(2,DER_FUNC3D,i,j) = R16 - R12*dy - R16*dxy + R16*dyz
        Dbas(3,DER_FUNC3D,i,j) = R16 + R12*dx + R13*dxy + R16*dyz
        Dbas(4,DER_FUNC3D,i,j) = R16 + R12*dy - R16*dxy + R16*dyz
        Dbas(5,DER_FUNC3D,i,j) = R16 - R12*dx + R13*dxy + R16*dyz
        Dbas(6,DER_FUNC3D,i,j) = R16 + R12*dz - R16*dxy - R13*dyz
      END DO
      
    END DO
    
  END IF
    
  !if x-, y- or z-derivatives are desired
  IF ((Bder(DER_DERIV3D_X)) .OR. (Bder(DER_DERIV3D_Y)) .OR. &
      (Bder(DER_DERIV3D_Z))) THEN
  
    DO j=1,nelements
      Dxj = 1.0_DP / Ddetj(:,j)
      
      !x-, y- and z-derivatives on reference element
      DO i=1,npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)
        dz = Dpoints(3,i,j)
        djx = -R13*dx
        Dhelp(1,1,i) = djx
        Dhelp(2,1,i) = djx
        Dhelp(3,1,i) = R16*(4.0_DP*dx + 3.0_DP)
        Dhelp(4,1,i) = djx
        Dhelp(5,1,i) = R16*(4.0_DP*dx - 3.0_DP)
        Dhelp(6,1,i) = djx
        djy = -R13*dy
        Dhelp(1,2,i) = djy
        Dhelp(2,2,i) = R16*(4.0_DP*dy - 3.0_DP)
        Dhelp(3,2,i) = djy
        Dhelp(4,2,i) = R16*(4.0_DP*dy + 3.0_DP)
        Dhelp(5,2,i) = djy
        Dhelp(6,2,i) = djy
        djz = -R13*dz
        Dhelp(1,3,i) = R16*(4.0_DP*dz - 3.0_DP)
        Dhelp(2,3,i) = djz
        Dhelp(3,3,i) = djz
        Dhelp(4,3,i) = djz
        Dhelp(5,3,i) = djz
        Dhelp(6,3,i) = R16*(4.0_DP*dz + 3.0_DP)
      END DO
        
      !x-derivatives on current element
!      IF (Bder(DER_DERIV3D_X)) THEN
        DO i=1,npoints
          djx = Djac(5,i,j)*Djac(9,i,j) - Djac(6,i,j)*Djac(8,i,j)
          djy = Djac(8,i,j)*Djac(3,i,j) - Djac(2,i,j)*Djac(9,i,j)
          djz = Djac(2,i,j)*Djac(6,i,j) - Djac(5,i,j)*Djac(3,i,j)
          Dbas(1,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(1,1,i) - djy*Dhelp(1,2,i) + djz*Dhelp(1,3,i))
          Dbas(2,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(2,1,i) - djy*Dhelp(2,2,i) + djz*Dhelp(2,3,i))
          Dbas(3,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(3,1,i) - djy*Dhelp(3,2,i) + djz*Dhelp(3,3,i))
          Dbas(4,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(4,1,i) - djy*Dhelp(4,2,i) + djz*Dhelp(4,3,i))
          Dbas(5,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(5,1,i) - djy*Dhelp(5,2,i) + djz*Dhelp(5,3,i))
          Dbas(6,DER_DERIV3D_X,i,j) = Dxj(i) * &
              (djx*Dhelp(6,1,i) - djy*Dhelp(6,2,i) + djz*Dhelp(6,3,i))
        END DO
!      ENDIF
      
      !y-derivatives on current element
!      IF (Bder(DER_DERIV3D_Y)) THEN
        DO i=1,npoints
          djx = Djac(7,i,j)*Djac(6,i,j) - Djac(4,i,j)*Djac(9,i,j)
          djy = Djac(1,i,j)*Djac(9,i,j) - Djac(7,i,j)*Djac(3,i,j)
          djz = Djac(4,i,j)*Djac(3,i,j) - Djac(1,i,j)*Djac(6,i,j)
          Dbas(1,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (-djx*Dhelp(1,1,i) + djy*Dhelp(1,2,i) - djz*Dhelp(1,3,i))
          Dbas(2,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (-djx*Dhelp(2,1,i) + djy*Dhelp(2,2,i) - djz*Dhelp(2,3,i))
          Dbas(3,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (-djx*Dhelp(3,1,i) + djy*Dhelp(3,2,i) - djz*Dhelp(3,3,i))
          Dbas(4,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (-djx*Dhelp(4,1,i) + djy*Dhelp(4,2,i) - djz*Dhelp(4,3,i))
          Dbas(5,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (-djx*Dhelp(5,1,i) + djy*Dhelp(5,2,i) - djz*Dhelp(5,3,i))
          Dbas(6,DER_DERIV3D_Y,i,j) = Dxj(i) * &
              (-djx*Dhelp(6,1,i) + djy*Dhelp(6,2,i) - djz*Dhelp(6,3,i))
        END DO
!      ENDIF

      !z-derivatives on current element
!      IF (Bder(DER_DERIV3D_Z)) THEN
        DO i=1,npoints
          djx = Djac(4,i,j)*Djac(8,i,j) - Djac(7,i,j)*Djac(5,i,j)
          djy = Djac(7,i,j)*Djac(2,i,j) - Djac(1,i,j)*Djac(8,i,j)
          djz = Djac(1,i,j)*Djac(5,i,j) - Djac(4,i,j)*Djac(2,i,j)
          Dbas(1,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(1,1,i) - djy*Dhelp(1,2,i) + djz*Dhelp(1,3,i))
          Dbas(2,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(2,1,i) - djy*Dhelp(2,2,i) + djz*Dhelp(2,3,i))
          Dbas(3,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(3,1,i) - djy*Dhelp(3,2,i) + djz*Dhelp(3,3,i))
          Dbas(4,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(4,1,i) - djy*Dhelp(4,2,i) + djz*Dhelp(4,3,i))
          Dbas(5,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(5,1,i) - djy*Dhelp(5,2,i) + djz*Dhelp(5,3,i))
          Dbas(6,DER_DERIV3D_Z,i,j) = Dxj(i) * &
              (djx*Dhelp(6,1,i) - djy*Dhelp(6,2,i) + djz*Dhelp(6,3,i))
        END DO
!      ENDIF
    END DO
      
  END IF
    
  END SUBROUTINE 

!**************************************************************************
! Element subroutines for Q1~ element, integral mean value based.
!**************************************************************************
! It is recommended that you have read (and understood) the documentation
! of the 2D EM30 element as this documentation of the 3D EM30 is quite
! similar, although not as detailed and without that fancy ASCII-art...
!
! The standard integral-based Q1~-element has the following six local
! basis functions:
!
!  p_1(x,y,z) = a1 + b1*x + c1*y + d1*z + e1*(x^2 - y^2) + f1*(y^2 - z^2)
!  p_2(x,y,z) = a2 + b2*x + c2*y + d2*z + e2*(x^2 - y^2) + f2*(y^2 - z^2)
!  p_3(x,y,z) = a3 + b3*x + c3*y + d3*z + e3*(x^2 - y^2) + f3*(y^2 - z^2)
!  p_4(x,y,z) = a4 + b4*x + c4*y + d4*z + e4*(x^2 - y^2) + f4*(y^2 - z^2)
!  p_5(x,y,z) = a5 + b5*x + c5*y + d5*z + e5*(x^2 - y^2) + f5*(y^2 - z^2)
!  p_6(x,y,z) = a6 + b6*x + c6*y + d6*z + e6*(x^2 - y^2) + f6*(y^2 - z^2)
!
! each of them designed in such a way, such that
!
!      1/|Xi| int_Xi p_j(x,y,z) d(x,y,z) = delta_ij
!
! with Xi being the i-th local face of a hexahedron.
!
! We will now use linear mapping between the face midpoints of our
! reference hexahedron [-1,1]^3 and the real hexahedron. So our
! mapping t:[-1,1]^3 -> R^3 will have the following form:
!
!   / X \                / t11 t12 t13 \   / x \ 
!   | Y | := t(x,y,z) := | t21 t22 t23 | * | y |
!   \ Z /                \ t31 t32 t33 /   \ z /
!
! This mapping should fulfill:
!
!   / eta_1 \   / t11 t12 t13 \   / 1 \ 
!   | eta_2 | = | t21 t22 t23 | * | 0 |                       (1)
!   \ eta_3 /   \ t31 t32 t33 /   \ 0 /
!  
!   / xi_1 \    / t11 t12 t13 \   / 0 \ 
!   | xi_2 |  = | t21 t22 t23 | * | 1 |                       (2)
!   \ xi_3 /    \ t31 t32 t33 /   \ 0 /
!
!   / rho_1 \   / t11 t12 t13 \   / 0 \ 
!   | rho_2 | = | t21 t22 t23 | * | 0 |                       (3)
!   \ rho_3 /   \ t31 t32 t33 /   \ 1 /
!
! where eta is the vector from the hexahedron midpoint to the midpoint of
! the third face, xi being the vector from the hexahedron midpoint to the
! midpoint of the sixth face and rho being the vector from the hexahedron
! midpoint to the midpoint of the second face.
! Note: Check the 2D EM30 documentation for an ASCII-art.
!
! The linear system given by the equations (1),(2) and (3) can easily
! be solved to calculate the unknown transformation coefficients tXY:
!
!   / X \                / eta_1 xi_1 rho_1 \   / x \ 
!   | Y | := t(x,y,z) := | eta_2 xi_2 rho_2 | * | y |
!   \ Z /                \ eta_3 xi_3 rho_3 /   \ z /
!
! Then, we set up the local basis functions in terms of xi, eta and rho,
! i.e. in the coordinate space of the new coordinate system, 
! with a new set of (unknown) coefficents ai, bi, etc::
!
!  P1(X,Y,Z) = a1 + b1*X + c1*Y + d1*Z + e1*(X^2 - Y^2) + f1*(Y^2 - Z^2)
!  P2(X,Y,Z) = a2 + b2*X + c2*Y + d2*Z + e2*(X^2 - Y^2) + f2*(Y^2 - Z^2)
!  P3(X,Y,Z) = a3 + b3*X + c3*Y + d3*Z + e3*(X^2 - Y^2) + f3*(Y^2 - Z^2)
!  P4(X,Y,Z) = a4 + b4*X + c4*Y + d4*Z + e4*(X^2 - Y^2) + f4*(Y^2 - Z^2)
!  P5(X,Y,Z) = a5 + b5*X + c5*Y + d5*Z + e5*(X^2 - Y^2) + f5*(Y^2 - Z^2)
!  P6(X,Y,Z) = a6 + b6*X + c6*Y + d6*Z + e6*(X^2 - Y^2) + f6*(Y^2 - Z^2)
!
! Later, we want to evaluate these Pi. Each Pi consists of a linear 
! combination of some coefficients ai, bi, etc. with the six monoms:
!
!  M1(X,Y,Z) := 1
!  M2(X,Y,Z) := X
!  M3(X,Y,Z) := Y
!  M4(X,Y,Z) := Z
!  M5(X,Y,Z) := X^2 - Y^2
!  M6(X,Y,Z) := Y^2 - Z^2
!
! To evaluate these Mi's in the new coordinate system, we concatenate them
! with the mapping t(.,.,.). As result, we get the six functions Fi,
! which are defined as functions at the bottom of this routine:
!
!  F1(x,y,z) := M1(t(x,y,z)) = 1
!  F2(x,y,z) := M2(t(x,y,z)) = eta_1*x + xi_1*y + rho_1*z
!  F3(x,y,z) := M3(t(x,y,z)) = eta_2*x + xi_2*y + rho_2*z
!  F4(x,y,z) := M4(t(x,y,z)) = eta_3*x + xi_3*y + rho_3*z
!
!  F5(x,y,z) := M5(t(x,y,z)) =   (eta_1*x + xi_1*y + rho_1*z)^2
!                              - (eta_2*x + xi_2*y + rho_2*z)^2
!                            =   (eta_1^2 - eta_2^2) * x^2 
!                              + (xi_1^2  - xi_2^2 ) * y^2
!                              + (rho1_^2 - rho_2^2) * z^2
!                              + 2*(eta_1*xi_1  - eta_2*xi_2 ) * x*y
!                              + 2*(eta_1*rho_1 - eta_2*rho_2) * x*z
!                              + 2*( xi_1*rho_1 -  xi_2*rho_2) * y*z
!
!  F6(x,y,z) := M6(t(x,y,z)) =   (eta_2*x + xi_2*y + rho_2*z)^2
!                              - (eta_3*x + xi_3*y + rho_3*z)^2
!                            =   (eta_2^2 - eta_3^2) * x^2 
!                              + (xi_2^2  - xi_3^2 ) * y^2
!                              + (rho1_^2 - rho_3^2) * z^2
!                              + 2*(eta_2*xi_2  - eta_3*xi_3 ) * x*y
!                              + 2*(eta_2*rho_2 - eta_3*rho_3) * x*z
!                              + 2*( xi_2*rho_2 -  xi_3*rho_3) * y*z
!
! So the polynomials have now the form:
!
!  Pi(t(x,y,z)) = ai*F1(x,y,z) + bi*F2(x,y,z) + ci*F3(x,y,z)
!               + di*F4(x,y,z) + ei*F5(x,y,z) + fi*F6(x,y,z)
!
! It doesn't matter whether the local coordinate system starts in (0,0,0) or in
! the midpoint of the element or whereever. As the rotation of the element
! coincides with the rotation of the new coordinate system, the polynomial 
! space is unisolvent and therefore exist the above local basis functions uniquely.
!
! The coefficients ai, bi, ci, di are to be calculated in such a way, that
!
!    1/|Xi| int_Xi Pj(t(x,y,z)) d(x,y,z) = delta_ij
!
! holds. The integral "1/|Xi| int_Xi ... d(x,y,z)" over the face Xi is approximated
! by a cubature formula with 8 points.
! This gives us six 6x6 systems for the computation of ai, bi, etc.
! The system matix A is given as A = {a_ij}, where 
!             a_ij := 1/|Xi| int_Xi Pj(t(x,y,z)) d(x,y,z)
!
! Once the system matrix is inverted, the A^-1 will hold the coefficients
! ai, bi, etc. which define the basis polynomials Pi with the desired property.
!
! Let's go for it...

!**************************************************************************
! Element subroutines for nonparametric 3D Q1~ element, integral mean value
! based.
! The routines are defines with the F95 PURE statement as they work 
! only on the parameters; helps some compilers in optimisation.
!**************************************************************************

!<subroutine>  

  PURE SUBROUTINE elem_EM30_3D (ieltyp, Dcoords, Djac, ddetj, Bder, &
                                Dpoint, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point. The coordinates are expected
  ! on the real element!
!</description>

  !<input>

  ! Element type identifier. Must be =EL_EM30_3D.
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
  
  ! Cartesian coordinates of the evaluation point on the real element.
  ! Dpoint(1) = x-coordinate,
  ! Dpoint(2) = y-coordinate
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

  REAL(DP), PARAMETER :: Q1 =   0.25_DP
  REAL(DP), PARAMETER :: Q2 =   1.0_DP /   3.0_DP
  REAL(DP), PARAMETER :: Q3 =   0.2_DP
  REAL(DP), PARAMETER :: Q4 =   0.6_DP
  REAL(DP), PARAMETER :: Q8 =  -9.0_DP /  16.0_DP
  REAL(DP), PARAMETER :: Q9 = 100.0_DP / 192.0_DP
  
  REAL(DP), DIMENSION(3,8) :: P                       ! cubature points
  REAL(DP) :: PP1,PP2,PP3,PP4,PP5,PS1,PS2,PQ1,PQ2,PQ3 ! cubature weights
  REAL(DP), DIMENSION(6,6) :: A, B                    ! coefficient matrices
  REAL(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC1,&           ! auxiliary coefficients
              CC2,CC3,CD1,CD2,CD3,CD4,CD5,&
              CD6,CE1,CE2,CE3,CE4,CE5,CE6
  REAL(DP), DIMENSION(6,10) :: COB                    ! monomial coefficients
  REAL(DP), DIMENSION(3,4,6) :: V                     ! face corner vertices
  REAL(DP), DIMENSION(2:6,8) :: F                     ! transformed monomials
  REAL(DP) :: dx,dy,dz
  INTEGER :: k,l
  
    ! Initialise the vertex-coordinates - we do this so we can re-use
    ! the code for the calculation of the cubature points, otherwise we
    ! could not loop over all faces of the hexahedron, and we would have
    ! to handle each one seperately...
    V(1:3,1,1) = Dcoords(1:3,1)
    V(1:3,2,1) = Dcoords(1:3,2)
    V(1:3,3,1) = Dcoords(1:3,3)
    V(1:3,4,1) = Dcoords(1:3,4)
    V(1:3,1,2) = Dcoords(1:3,1)
    V(1:3,2,2) = Dcoords(1:3,2)
    V(1:3,3,2) = Dcoords(1:3,6)
    V(1:3,4,2) = Dcoords(1:3,5)
    V(1:3,1,3) = Dcoords(1:3,2)
    V(1:3,2,3) = Dcoords(1:3,3)
    V(1:3,3,3) = Dcoords(1:3,7)
    V(1:3,4,3) = Dcoords(1:3,6)
    V(1:3,1,4) = Dcoords(1:3,3)
    V(1:3,2,4) = Dcoords(1:3,4)
    V(1:3,3,4) = Dcoords(1:3,8)
    V(1:3,4,4) = Dcoords(1:3,7)
    V(1:3,1,5) = Dcoords(1:3,4)
    V(1:3,2,5) = Dcoords(1:3,1)
    V(1:3,3,5) = Dcoords(1:3,5)
    V(1:3,4,5) = Dcoords(1:3,8)
    V(1:3,1,6) = Dcoords(1:3,5)
    V(1:3,2,6) = Dcoords(1:3,6)
    V(1:3,3,6) = Dcoords(1:3,7)
    V(1:3,4,6) = Dcoords(1:3,8)

    ! Calculate the coefficients of the transformed monomials F2..F6:
    ! Note: the vectors eta, xi and rho are:
    !  eta = (CA1, CB1, CC1)
    !  xi  = (CA2, CB2, CC2)
    !  rho = (CA3, CB3, CC3)
    CA1 = Q1*((Dcoords(1,1)+Dcoords(1,2)+Dcoords(1,3)+Dcoords(1,4)) &
             -(Dcoords(1,5)+Dcoords(1,6)+Dcoords(1,7)+Dcoords(1,8)))
    CB1 = Q1*((Dcoords(2,1)+Dcoords(2,2)+Dcoords(2,3)+Dcoords(2,4)) &
             -(Dcoords(2,5)+Dcoords(2,6)+Dcoords(2,7)+Dcoords(2,8)))
    CC1 = Q1*((Dcoords(3,1)+Dcoords(3,2)+Dcoords(3,3)+Dcoords(3,4)) &
             -(Dcoords(3,5)+Dcoords(3,6)+Dcoords(3,7)+Dcoords(3,8)))
    CA2 = Q1*((Dcoords(1,1)+Dcoords(1,2)+Dcoords(1,6)+Dcoords(1,5)) &
             -(Dcoords(1,3)+Dcoords(1,7)+Dcoords(1,8)+Dcoords(1,4)))
    CB2 = Q1*((Dcoords(2,1)+Dcoords(2,2)+Dcoords(2,6)+Dcoords(2,5)) &
             -(Dcoords(2,3)+Dcoords(2,7)+Dcoords(2,8)+Dcoords(2,4)))
    CC2 = Q1*((Dcoords(3,1)+Dcoords(3,2)+Dcoords(3,6)+Dcoords(3,5)) &
             -(Dcoords(3,3)+Dcoords(3,7)+Dcoords(3,8)+Dcoords(3,4)))
    CA3 = Q1*((Dcoords(1,2)+Dcoords(1,3)+Dcoords(1,7)+Dcoords(1,6)) &
             -(Dcoords(1,1)+Dcoords(1,4)+Dcoords(1,8)+Dcoords(1,5)))
    CB3 = Q1*((Dcoords(2,2)+Dcoords(2,3)+Dcoords(2,7)+Dcoords(2,6)) &
             -(Dcoords(2,1)+Dcoords(2,4)+Dcoords(2,8)+Dcoords(2,5)))
    CC3 = Q1*((Dcoords(3,2)+Dcoords(3,3)+Dcoords(3,7)+Dcoords(3,6)) &
             -(Dcoords(3,1)+Dcoords(3,4)+Dcoords(3,8)+Dcoords(3,5)))
    ! Coefficients for F5
    CD1 = CA1**2 - CA2**2
    CD2 = CB1**2 - CB2**2
    CD3 = CC1**2 - CC2**2
    CD4 = 2.0_DP*(CA1*CB1 - CA2*CB2)
    CD5 = 2.0_DP*(CA1*CC1 - CA2*CC2)
    CD6 = 2.0_DP*(CB1*CC1 - CB2*CC2)
    ! Coefficients for F6
    CE1 = CA2**2 - CA3**2
    CE2 = CB2**2 - CB3**2
    CE3 = CC2**2 - CC3**2
    CE4 = 2.0_DP*(CA2*CB2 - CA3*CB3)
    CE5 = 2.0_DP*(CA2*CC2 - CA3*CC3)
    CE6 = 2.0_DP*(CB2*CC2 - CB3*CC3)

    ! We now have the coefficients of the transformed monomials in our
    ! new coordinate system - so we now can start with the assembly of
    ! our system matrix A
    
    ! Loop over all faces of the hexahedron
    DO k=1,6
    
      ! Calculate the eight cubature points for this face
      P(1,1)=Q2*(V(1,1,k)+V(1,2,k)+V(1,3,k))
      P(2,1)=Q2*(V(2,1,k)+V(2,2,k)+V(2,3,k))
      P(3,1)=Q2*(V(3,1,k)+V(3,2,k)+V(3,3,k))
      P(1,2)=Q4*V(1,1,k)+Q3*V(1,2,k)+Q3*V(1,3,k)
      P(2,2)=Q4*V(2,1,k)+Q3*V(2,2,k)+Q3*V(2,3,k)
      P(3,2)=Q4*V(3,1,k)+Q3*V(3,2,k)+Q3*V(3,3,k)
      P(1,3)=Q3*V(1,1,k)+Q4*V(1,2,k)+Q3*V(1,3,k)
      P(2,3)=Q3*V(2,1,k)+Q4*V(2,2,k)+Q3*V(2,3,k)
      P(3,3)=Q3*V(3,1,k)+Q4*V(3,2,k)+Q3*V(3,3,k)
      P(1,4)=Q3*V(1,1,k)+Q3*V(1,2,k)+Q4*V(1,3,k)
      P(2,4)=Q3*V(2,1,k)+Q3*V(2,2,k)+Q4*V(2,3,k)
      P(3,4)=Q3*V(3,1,k)+Q3*V(3,2,k)+Q4*V(3,3,k)
      P(1,5)=Q2*(V(1,1,k)+V(1,3,k)+V(1,4,k))
      P(2,5)=Q2*(V(2,1,k)+V(2,3,k)+V(2,4,k))
      P(3,5)=Q2*(V(3,1,k)+V(3,3,k)+V(3,4,k))
      P(1,6)=Q4*V(1,1,k)+Q3*V(1,3,k)+Q3*V(1,4,k)
      P(2,6)=Q4*V(2,1,k)+Q3*V(2,3,k)+Q3*V(2,4,k)
      P(3,6)=Q4*V(3,1,k)+Q3*V(3,3,k)+Q3*V(3,4,k)
      P(1,7)=Q3*V(1,1,k)+Q4*V(1,3,k)+Q3*V(1,4,k)
      P(2,7)=Q3*V(2,1,k)+Q4*V(2,3,k)+Q3*V(2,4,k)
      P(3,7)=Q3*V(3,1,k)+Q4*V(3,3,k)+Q3*V(3,4,k)
      P(1,8)=Q3*V(1,1,k)+Q3*V(1,3,k)+Q4*V(1,4,k)
      P(2,8)=Q3*V(2,1,k)+Q3*V(2,3,k)+Q4*V(2,4,k)
      P(3,8)=Q3*V(3,1,k)+Q3*V(3,3,k)+Q4*V(3,4,k)
      
      ! And calculate the weights for the cubature formula
      PP1=SQRT((V(1,1,k)-V(1,2,k))**2+(V(2,1,k)-V(2,2,k))**2+&
          (V(3,1,k)-V(3,2,k))**2)
      PP2=SQRT((V(1,2,k)-V(1,3,k))**2+(V(2,2,k)-V(2,3,k))**2+&
          (V(3,2,k)-V(3,3,k))**2)
      PP3=SQRT((V(1,1,k)-V(1,3,k))**2+(V(2,1,k)-V(2,3,k))**2+&
          (V(3,1,k)-V(3,3,k))**2)
      PP4=SQRT((V(1,3,k)-V(1,4,k))**2+(V(2,3,k)-V(2,4,k))**2+&
          (V(3,3,k)-V(3,4,k))**2)
      PP5=SQRT((V(1,1,k)-V(1,4,k))**2+(V(2,1,k)-V(2,4,k))**2+&
          (V(3,1,k)-V(3,4,k))**2)
      PS1=(PP1+PP2+PP3)*0.5_DP
      PS2=(PP3+PP4+PP5)*0.5_DP
      PQ1=SQRT(PS1*(PS1-PP1)*(PS1-PP2)*(PS1-PP3))
      PQ2=SQRT(PS2*(PS2-PP3)*(PS2-PP4)*(PS2-PP5))
      PQ3=1.0_DP / (PQ1 + PQ2)
      
      ! Evalute the F2..F6 in the eight cubature points
      ! Remember: F1 = 1
      DO l=1,8
        F(2,l) = CA1*P(1,l) + CB1*P(2,l) + CC1*P(3,l)
        F(3,l) = CA2*P(1,l) + CB2*P(2,l) + CC2*P(3,l)
        F(4,l) = CA3*P(1,l) + CB3*P(2,l) + CC3*P(3,l)
        F(5,l) = CD1*P(1,l)**2 + CD2*P(2,l)**2 + CD3*P(3,l)**2 &
               + CD4*P(1,l)*P(2,l) + CD5*P(1,l)*P(3,l) + CD6*P(2,l)*P(3,l)
        F(6,l) = CE1*P(1,l)**2 + CE2*P(2,l)**2 + CE3*P(3,l)**2 &
               + CE4*P(1,l)*P(2,l) + CE5*P(1,l)*P(3,l) + CE6*P(2,l)*P(3,l)
      END DO

      ! Build up the matrix entries
      A(k,1)= 1.0_DP
      DO l=2,6
        A(k,l) = PQ3*(PQ1*(Q8*F(l,1)+Q9*(F(l,2)+F(l,3)+F(l,4))) &
                     +PQ2*(Q8*F(l,5)+Q9*(F(l,6)+F(l,7)+F(l,8))))
      END DO
               
    END DO
    
    ! Now we have the coeffienct matrix - so invert it to get the
    ! coefficients of our basis polynomials.
    CALL mprim_invert6x6MatrixDirectDble(A,B)

    ! Basically, we now have everything we need to start with the evaluation
    ! of the basis polynomial, however, our polynomials are currently given as:
    !
    !  Pi(t(x,y,z)) = ai*F1(x,y,z) + bi*F2(x,y,z) + ci*F3(x,y,z)
    !               + di*F4(x,y,z) + ei*F5(x,y,z) + fi*F6(x,y,z)
    !
    ! We would like to transfom the Pi back to monomial base, so that the
    ! Pi can be written as:
    !
    ! Pi(t(x,y,z)) = COB(i,1)*x^2 + COB(i,2)*y^2 + COB(i,3)*z^2
    !              + COB(i,4)*x*y + COB(i,5)*x*z + COB(i,6)*y*z
    !              + COB(i,7)*x + COB(i,8)*y + COB(i,9)*z + COB(1,10)
    !
    ! Spo transform the coefficiencts into monomial base
    DO k=1,6
      COB(k,1) = B(6,k)*CE1 + B(5,k)*CD1
      COB(k,2) = B(6,k)*CE2 + B(5,k)*CD2
      COB(k,3) = B(6,k)*CE3 + B(5,k)*CD3
      COB(k,4) = B(6,k)*CE4 + B(5,k)*CD4
      COB(k,5) = B(6,k)*CE5 + B(5,k)*CD5
      COB(k,6) = B(6,k)*CE6 + B(5,k)*CD6
      COB(k,7) = B(4,k)*CA3 + B(3,k)*CA2 + B(2,k)*CA1
      COB(k,8) = B(4,k)*CB3 + B(3,k)*CB2 + B(2,k)*CB1
      COB(k,9) = B(4,k)*CC3 + B(3,k)*CC2 + B(2,k)*CC1
      COB(k,10)= B(1,k)
    END DO
    
    ! Now we can finally start to evaluate...
    dx = Dpoint(1)
    dy = Dpoint(2)
    dz = Dpoint(3)

    ! Function values
    IF (BDER(DER_FUNC3D)) THEN
      Dbas(1,DER_FUNC3D)= COB(1,1)*dx**2+COB(1,2)*dy**2+COB(1,3)*dz**2&
               +COB(1,4)*dx*dy+COB(1,5)*dx*dz+COB(1,6)*dy*dz&
               +COB(1,7)*dx+COB(1,8)*dy+COB(1,9)*dz+COB(1,10)
      Dbas(2,DER_FUNC3D)= COB(2,1)*dx**2+COB(2,2)*dy**2+COB(2,3)*dz**2&
               +COB(2,4)*dx*dy+COB(2,5)*dx*dz+COB(2,6)*dy*dz&
               +COB(2,7)*dx+COB(2,8)*dy+COB(2,9)*dz+COB(2,10)
      Dbas(3,DER_FUNC3D)= COB(3,1)*dx**2+COB(3,2)*dy**2+COB(3,3)*dz**2&
               +COB(3,4)*dx*dy+COB(3,5)*dx*dz+COB(3,6)*dy*dz&
               +COB(3,7)*dx+COB(3,8)*dy+COB(3,9)*dz+COB(3,10)
      Dbas(4,DER_FUNC3D)= COB(4,1)*dx**2+COB(4,2)*dy**2+COB(4,3)*dz**2&
               +COB(4,4)*dx*dy+COB(4,5)*dx*dz+COB(4,6)*dy*dz&
               +COB(4,7)*dx+COB(4,8)*dy+COB(4,9)*dz+COB(4,10)
      Dbas(5,DER_FUNC3D)= COB(5,1)*dx**2+COB(5,2)*dy**2+COB(5,3)*dz**2&
               +COB(5,4)*dx*dy+COB(5,5)*dx*dz+COB(5,6)*dy*dz&
               +COB(5,7)*dx+COB(5,8)*dy+COB(5,9)*dz+COB(5,10)
      Dbas(6,DER_FUNC3D)= COB(6,1)*dx**2+COB(6,2)*dy**2+COB(6,3)*dz**2&
               +COB(6,4)*dx*dy+COB(6,5)*dx*dz+COB(6,6)*dy*dz&
               +COB(6,7)*dx+COB(6,8)*dy+COB(6,9)*dz+COB(6,10)
    END IF

    ! X-derivatives
    IF(BDER(DER_DERIV3D_X)) THEN
      Dbas(1,DER_DERIV3D_X)=2.0_DP*COB(1,1)*dx+COB(1,4)*dy+COB(1,5)*dz+COB(1,7)
      Dbas(2,DER_DERIV3D_X)=2.0_DP*COB(2,1)*dx+COB(2,4)*dy+COB(2,5)*dz+COB(2,7)
      Dbas(3,DER_DERIV3D_X)=2.0_DP*COB(3,1)*dx+COB(3,4)*dy+COB(3,5)*dz+COB(3,7)
      Dbas(4,DER_DERIV3D_X)=2.0_DP*COB(4,1)*dx+COB(4,4)*dy+COB(4,5)*dz+COB(4,7)
      Dbas(5,DER_DERIV3D_X)=2.0_DP*COB(5,1)*dx+COB(5,4)*dy+COB(5,5)*dz+COB(5,7)
      Dbas(6,DER_DERIV3D_X)=2.0_DP*COB(6,1)*dx+COB(6,4)*dy+COB(6,5)*dz+COB(6,7)
    END IF

    ! Y-derivatives
    IF(BDER(DER_DERIV3D_Y)) THEN
      Dbas(1,DER_DERIV3D_Y)=2.0_DP*COB(1,2)*dy+COB(1,4)*dx+COB(1,6)*dz+COB(1,8)
      Dbas(2,DER_DERIV3D_Y)=2.0_DP*COB(2,2)*dy+COB(2,4)*dx+COB(2,6)*dz+COB(2,8)
      Dbas(3,DER_DERIV3D_Y)=2.0_DP*COB(3,2)*dy+COB(3,4)*dx+COB(3,6)*dz+COB(3,8)
      Dbas(4,DER_DERIV3D_Y)=2.0_DP*COB(4,2)*dy+COB(4,4)*dx+COB(4,6)*dz+COB(4,8)
      Dbas(5,DER_DERIV3D_Y)=2.0_DP*COB(5,2)*dy+COB(5,4)*dx+COB(5,6)*dz+COB(5,8)
      Dbas(6,DER_DERIV3D_Y)=2.0_DP*COB(6,2)*dy+COB(6,4)*dx+COB(6,6)*dz+COB(6,8)
    END IF
    
    ! Z-derivatives
    IF(BDER(DER_DERIV3D_Z)) THEN
      Dbas(1,DER_DERIV3D_Z)=2.0_DP*COB(1,3)*dz+COB(1,5)*dx+COB(1,6)*dy+COB(1,9)
      Dbas(2,DER_DERIV3D_Z)=2.0_DP*COB(2,3)*dz+COB(2,5)*dx+COB(2,6)*dy+COB(2,9)
      Dbas(3,DER_DERIV3D_Z)=2.0_DP*COB(3,3)*dz+COB(3,5)*dx+COB(3,6)*dy+COB(3,9)
      Dbas(4,DER_DERIV3D_Z)=2.0_DP*COB(4,3)*dz+COB(4,5)*dx+COB(4,6)*dy+COB(4,9)
      Dbas(5,DER_DERIV3D_Z)=2.0_DP*COB(5,3)*dz+COB(5,5)*dx+COB(5,6)*dy+COB(5,9)
      Dbas(6,DER_DERIV3D_Z)=2.0_DP*COB(6,3)*dz+COB(6,5)*dx+COB(6,6)*dy+COB(6,9)
    END IF

    ! That's it
    
  END SUBROUTINE


  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_EM30_3D_mult (ieltyp, Dcoords, Djac, Ddetj, &
                                    Bder, Dbas, npoints, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element.
!</description>

!<input>
  ! Element type identifier. Must be =EL_EM30_3D.
  INTEGER(I32), INTENT(IN)  :: ieltyp

  ! Number of points where to evalate the basis functions.
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE).
  !  Dcoords(1,.)=x-coordinates,
  !  Dcoords(2,.)=y-coordinates,
  !  Dcoords(3,.)=z-coordinates.
  ! furthermore:
  !  Dcoords(:,i) = Coordinates of vertex i
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(3,1)
  !  Djac(4,i) = J_i(1,2)
  !  Djac(5,i) = J_i(2,2)
  !  Djac(6,i) = J_i(3,2)
  !  Djac(7,i) = J_i(1,3)
  !  Djac(8,i) = J_i(2,3)
  !  Djac(9,i) = J_i(3,3)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  REAL(DP), DIMENSION(EL_NJACENTRIES3D,npoints), INTENT(IN) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i) = Determinant of point i
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
  ! DIMENSION(#space dimensions,npoints,nelements).
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates,
  !  Dpoints(3,.)=z-coordinates.
  ! furthermore:
  !  Dpoints(:,i) = Coordinates of point i
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dpoints
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j,k) defines the value of the i'th 
  !   basis function of the finite element k in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i'th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.,.) is undefined.
  !REAL(DP), DIMENSION(EL_MAXNBAS,EL_MAXNDER,npoints), INTENT(OUT) :: Dbas
  REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: Dbas
!</output>

! </subroutine>
  REAL(DP), PARAMETER :: Q1 =   0.25_DP
  REAL(DP), PARAMETER :: Q2 =   1.0_DP /   3.0_DP
  REAL(DP), PARAMETER :: Q3 =   0.2_DP
  REAL(DP), PARAMETER :: Q4 =   0.6_DP
  REAL(DP), PARAMETER :: Q8 =  -9.0_DP /  16.0_DP
  REAL(DP), PARAMETER :: Q9 = 100.0_DP / 192.0_DP
  
  REAL(DP), DIMENSION(3,8) :: P                       ! cubature points
  REAL(DP) :: PP1,PP2,PP3,PP4,PP5,PS1,PS2,PQ1,PQ2,PQ3 ! cubature weights
  REAL(DP), DIMENSION(6,6) :: A, B                    ! coefficient matrices
  REAL(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC1,&           ! auxiliary coefficients
              CC2,CC3,CD1,CD2,CD3,CD4,CD5,&
              CD6,CE1,CE2,CE3,CE4,CE5,CE6
  REAL(DP), DIMENSION(6,10) :: COB                    ! monomial coefficients
  REAL(DP), DIMENSION(3,4,6) :: V                     ! face corner vertices
  REAL(DP), DIMENSION(2:6,8) :: F                     ! transformed monomials
  REAL(DP) :: dx,dy,dz
  INTEGER :: i,k,l
  
    ! Initialise the vertex-coordinates
    V(1:3,1,1) = Dcoords(1:3,1)
    V(1:3,2,1) = Dcoords(1:3,2)
    V(1:3,3,1) = Dcoords(1:3,3)
    V(1:3,4,1) = Dcoords(1:3,4)
    V(1:3,1,2) = Dcoords(1:3,1)
    V(1:3,2,2) = Dcoords(1:3,2)
    V(1:3,3,2) = Dcoords(1:3,6)
    V(1:3,4,2) = Dcoords(1:3,5)
    V(1:3,1,3) = Dcoords(1:3,2)
    V(1:3,2,3) = Dcoords(1:3,3)
    V(1:3,3,3) = Dcoords(1:3,7)
    V(1:3,4,3) = Dcoords(1:3,6)
    V(1:3,1,4) = Dcoords(1:3,3)
    V(1:3,2,4) = Dcoords(1:3,4)
    V(1:3,3,4) = Dcoords(1:3,8)
    V(1:3,4,4) = Dcoords(1:3,7)
    V(1:3,1,5) = Dcoords(1:3,4)
    V(1:3,2,5) = Dcoords(1:3,1)
    V(1:3,3,5) = Dcoords(1:3,5)
    V(1:3,4,5) = Dcoords(1:3,8)
    V(1:3,1,6) = Dcoords(1:3,5)
    V(1:3,2,6) = Dcoords(1:3,6)
    V(1:3,3,6) = Dcoords(1:3,7)
    V(1:3,4,6) = Dcoords(1:3,8)

    ! Calculate the coefficients of the Fi
    CA1 = Q1*((Dcoords(1,1)+Dcoords(1,2)+Dcoords(1,3)+Dcoords(1,4)) &
             -(Dcoords(1,5)+Dcoords(1,6)+Dcoords(1,7)+Dcoords(1,8)))
    CB1 = Q1*((Dcoords(2,1)+Dcoords(2,2)+Dcoords(2,3)+Dcoords(2,4)) &
             -(Dcoords(2,5)+Dcoords(2,6)+Dcoords(2,7)+Dcoords(2,8)))
    CC1 = Q1*((Dcoords(3,1)+Dcoords(3,2)+Dcoords(3,3)+Dcoords(3,4)) &
             -(Dcoords(3,5)+Dcoords(3,6)+Dcoords(3,7)+Dcoords(3,8)))
    CA2 = Q1*((Dcoords(1,1)+Dcoords(1,2)+Dcoords(1,6)+Dcoords(1,5)) &
             -(Dcoords(1,3)+Dcoords(1,7)+Dcoords(1,8)+Dcoords(1,4)))
    CB2 = Q1*((Dcoords(2,1)+Dcoords(2,2)+Dcoords(2,6)+Dcoords(2,5)) &
             -(Dcoords(2,3)+Dcoords(2,7)+Dcoords(2,8)+Dcoords(2,4)))
    CC2 = Q1*((Dcoords(3,1)+Dcoords(3,2)+Dcoords(3,6)+Dcoords(3,5)) &
             -(Dcoords(3,3)+Dcoords(3,7)+Dcoords(3,8)+Dcoords(3,4)))
    CA3 = Q1*((Dcoords(1,2)+Dcoords(1,3)+Dcoords(1,7)+Dcoords(1,6)) &
             -(Dcoords(1,1)+Dcoords(1,4)+Dcoords(1,8)+Dcoords(1,5)))
    CB3 = Q1*((Dcoords(2,2)+Dcoords(2,3)+Dcoords(2,7)+Dcoords(2,6)) &
             -(Dcoords(2,1)+Dcoords(2,4)+Dcoords(2,8)+Dcoords(2,5)))
    CC3 = Q1*((Dcoords(3,2)+Dcoords(3,3)+Dcoords(3,7)+Dcoords(3,6)) &
             -(Dcoords(3,1)+Dcoords(3,4)+Dcoords(3,8)+Dcoords(3,5)))
    CD1 = CA1**2 - CA2**2
    CD2 = CB1**2 - CB2**2
    CD3 = CC1**2 - CC2**2
    CD4 = 2.0_DP*(CA1*CB1 - CA2*CB2)
    CD5 = 2.0_DP*(CA1*CC1 - CA2*CC2)
    CD6 = 2.0_DP*(CB1*CC1 - CB2*CC2)
    CE1 = CA2**2 - CA3**2
    CE2 = CB2**2 - CB3**2
    CE3 = CC2**2 - CC3**2
    CE4 = 2.0_DP*(CA2*CB2 - CA3*CB3)
    CE5 = 2.0_DP*(CA2*CC2 - CA3*CC3)
    CE6 = 2.0_DP*(CB2*CC2 - CB3*CC3)

    DO k=1,6
      ! Calculate the eight cubature poknts
      P(1,1)=Q2*(V(1,1,k)+V(1,2,k)+V(1,3,k))
      P(2,1)=Q2*(V(2,1,k)+V(2,2,k)+V(2,3,k))
      P(3,1)=Q2*(V(3,1,k)+V(3,2,k)+V(3,3,k))
      P(1,2)=Q4*V(1,1,k)+Q3*V(1,2,k)+Q3*V(1,3,k)
      P(2,2)=Q4*V(2,1,k)+Q3*V(2,2,k)+Q3*V(2,3,k)
      P(3,2)=Q4*V(3,1,k)+Q3*V(3,2,k)+Q3*V(3,3,k)
      P(1,3)=Q3*V(1,1,k)+Q4*V(1,2,k)+Q3*V(1,3,k)
      P(2,3)=Q3*V(2,1,k)+Q4*V(2,2,k)+Q3*V(2,3,k)
      P(3,3)=Q3*V(3,1,k)+Q4*V(3,2,k)+Q3*V(3,3,k)
      P(1,4)=Q3*V(1,1,k)+Q3*V(1,2,k)+Q4*V(1,3,k)
      P(2,4)=Q3*V(2,1,k)+Q3*V(2,2,k)+Q4*V(2,3,k)
      P(3,4)=Q3*V(3,1,k)+Q3*V(3,2,k)+Q4*V(3,3,k)
      P(1,5)=Q2*(V(1,1,k)+V(1,3,k)+V(1,4,k))
      P(2,5)=Q2*(V(2,1,k)+V(2,3,k)+V(2,4,k))
      P(3,5)=Q2*(V(3,1,k)+V(3,3,k)+V(3,4,k))
      P(1,6)=Q4*V(1,1,k)+Q3*V(1,3,k)+Q3*V(1,4,k)
      P(2,6)=Q4*V(2,1,k)+Q3*V(2,3,k)+Q3*V(2,4,k)
      P(3,6)=Q4*V(3,1,k)+Q3*V(3,3,k)+Q3*V(3,4,k)
      P(1,7)=Q3*V(1,1,k)+Q4*V(1,3,k)+Q3*V(1,4,k)
      P(2,7)=Q3*V(2,1,k)+Q4*V(2,3,k)+Q3*V(2,4,k)
      P(3,7)=Q3*V(3,1,k)+Q4*V(3,3,k)+Q3*V(3,4,k)
      P(1,8)=Q3*V(1,1,k)+Q3*V(1,3,k)+Q4*V(1,4,k)
      P(2,8)=Q3*V(2,1,k)+Q3*V(2,3,k)+Q4*V(2,4,k)
      P(3,8)=Q3*V(3,1,k)+Q3*V(3,3,k)+Q4*V(3,4,k)
      
      ! And calculate the weights for the cubature formula
      PP1=SQRT((V(1,1,k)-V(1,2,k))**2+(V(2,1,k)-V(2,2,k))**2+&
          (V(3,1,k)-V(3,2,k))**2)
      PP2=SQRT((V(1,2,k)-V(1,3,k))**2+(V(2,2,k)-V(2,3,k))**2+&
          (V(3,2,k)-V(3,3,k))**2)
      PP3=SQRT((V(1,1,k)-V(1,3,k))**2+(V(2,1,k)-V(2,3,k))**2+&
          (V(3,1,k)-V(3,3,k))**2)
      PP4=SQRT((V(1,3,k)-V(1,4,k))**2+(V(2,3,k)-V(2,4,k))**2+&
          (V(3,3,k)-V(3,4,k))**2)
      PP5=SQRT((V(1,1,k)-V(1,4,k))**2+(V(2,1,k)-V(2,4,k))**2+&
          (V(3,1,k)-V(3,4,k))**2)
      PS1=(PP1+PP2+PP3)*0.5_DP
      PS2=(PP3+PP4+PP5)*0.5_DP
      PQ1=SQRT(PS1*(PS1-PP1)*(PS1-PP2)*(PS1-PP3))
      PQ2=SQRT(PS2*(PS2-PP3)*(PS2-PP4)*(PS2-PP5))
      PQ3=1.0_DP / (PQ1 + PQ2)
      
      ! Evalute the F2..F6 in the eight cubature points
      ! Remember: F1 = 1
      DO l=1,8
        F(2,l) = CA1*P(1,l) + CB1*P(2,l) + CC1*P(3,l)
        F(3,l) = CA2*P(1,l) + CB2*P(2,l) + CC2*P(3,l)
        F(4,l) = CA3*P(1,l) + CB3*P(2,l) + CC3*P(3,l)
        F(5,l) = CD1*P(1,l)**2 + CD2*P(2,l)**2 + CD3*P(3,l)**2 &
               + CD4*P(1,l)*P(2,l) + CD5*P(1,l)*P(3,l) + CD6*P(2,l)*P(3,l)
        F(6,l) = CE1*P(1,l)**2 + CE2*P(2,l)**2 + CE3*P(3,l)**2 &
               + CE4*P(1,l)*P(2,l) + CE5*P(1,l)*P(3,l) + CE6*P(2,l)*P(3,l)
      END DO

      ! Build up the matrix entries
      A(k,1)= 1.0_DP
      DO l=2,6
        A(k,l) = PQ3*(PQ1*(Q8*F(l,1)+Q9*(F(l,2)+F(l,3)+F(l,4))) &
                     +PQ2*(Q8*F(l,5)+Q9*(F(l,6)+F(l,7)+F(l,8))))
      END DO
               
    END DO
    
    ! Now we have the coeffienct matrix - so invert it
    CALL mprim_invert6x6MatrixDirectDble(A,B)

    ! Transform coefficiencts into monomial base
    DO k=1,6
      COB(k,1) = B(6,k)*CE1 + B(5,k)*CD1
      COB(k,2) = B(6,k)*CE2 + B(5,k)*CD2
      COB(k,3) = B(6,k)*CE3 + B(5,k)*CD3
      COB(k,4) = B(6,k)*CE4 + B(5,k)*CD4
      COB(k,5) = B(6,k)*CE5 + B(5,k)*CD5
      COB(k,6) = B(6,k)*CE6 + B(5,k)*CD6
      COB(k,7) = B(4,k)*CA3 + B(3,k)*CA2 + B(2,k)*CA1
      COB(k,8) = B(4,k)*CB3 + B(3,k)*CB2 + B(2,k)*CB1
      COB(k,9) = B(4,k)*CC3 + B(3,k)*CC2 + B(2,k)*CC1
      COB(k,10)= B(1,k)
    END DO
    
    ! Function values
    IF (BDER(DER_FUNC3D)) THEN
      DO i=1, npoints
        dx = Dpoints(1,i)
        dy = Dpoints(2,i)
        dz = Dpoints(3,i)
        Dbas(1,DER_FUNC3D,i)= COB(1,1)*dx**2+COB(1,2)*dy**2+COB(1,3)*dz**2&
                 +COB(1,4)*dx*dy+COB(1,5)*dx*dz+COB(1,6)*dy*dz&
                 +COB(1,7)*dx+COB(1,8)*dy+COB(1,9)*dz+COB(1,10)
        Dbas(2,DER_FUNC3D,i)= COB(2,1)*dx**2+COB(2,2)*dy**2+COB(2,3)*dz**2&
                 +COB(2,4)*dx*dy+COB(2,5)*dx*dz+COB(2,6)*dy*dz&
                 +COB(2,7)*dx+COB(2,8)*dy+COB(2,9)*dz+COB(2,10)
        Dbas(3,DER_FUNC3D,i)= COB(3,1)*dx**2+COB(3,2)*dy**2+COB(3,3)*dz**2&
                 +COB(3,4)*dx*dy+COB(3,5)*dx*dz+COB(3,6)*dy*dz&
                 +COB(3,7)*dx+COB(3,8)*dy+COB(3,9)*dz+COB(3,10)
        Dbas(4,DER_FUNC3D,i)= COB(4,1)*dx**2+COB(4,2)*dy**2+COB(4,3)*dz**2&
                 +COB(4,4)*dx*dy+COB(4,5)*dx*dz+COB(4,6)*dy*dz&
                 +COB(4,7)*dx+COB(4,8)*dy+COB(4,9)*dz+COB(4,10)
        Dbas(5,DER_FUNC3D,i)= COB(5,1)*dx**2+COB(5,2)*dy**2+COB(5,3)*dz**2&
                 +COB(5,4)*dx*dy+COB(5,5)*dx*dz+COB(5,6)*dy*dz&
                 +COB(5,7)*dx+COB(5,8)*dy+COB(5,9)*dz+COB(5,10)
        Dbas(6,DER_FUNC3D,i)= COB(6,1)*dx**2+COB(6,2)*dy**2+COB(6,3)*dz**2&
                 +COB(6,4)*dx*dy+COB(6,5)*dx*dz+COB(6,6)*dy*dz&
                 +COB(6,7)*dx+COB(6,8)*dy+COB(6,9)*dz+COB(6,10)
      END DO
    END IF

    ! X-derivatives
    IF(BDER(DER_DERIV3D_X)) THEN
      DO i=1, npoints
        dx = Dpoints(1,i)
        dy = Dpoints(2,i)
        dz = Dpoints(3,i)
        Dbas(1,DER_DERIV3D_X,i)=2.0_DP*COB(1,1)*dx+COB(1,4)*dy+COB(1,5)*dz+COB(1,7)
        Dbas(2,DER_DERIV3D_X,i)=2.0_DP*COB(2,1)*dx+COB(2,4)*dy+COB(2,5)*dz+COB(2,7)
        Dbas(3,DER_DERIV3D_X,i)=2.0_DP*COB(3,1)*dx+COB(3,4)*dy+COB(3,5)*dz+COB(3,7)
        Dbas(4,DER_DERIV3D_X,i)=2.0_DP*COB(4,1)*dx+COB(4,4)*dy+COB(4,5)*dz+COB(4,7)
        Dbas(5,DER_DERIV3D_X,i)=2.0_DP*COB(5,1)*dx+COB(5,4)*dy+COB(5,5)*dz+COB(5,7)
        Dbas(6,DER_DERIV3D_X,i)=2.0_DP*COB(6,1)*dx+COB(6,4)*dy+COB(6,5)*dz+COB(6,7)
      END DO
    END IF

    ! Y-derivatives
    IF(BDER(DER_DERIV3D_Y)) THEN
      DO i=1, npoints
        dx = Dpoints(1,i)
        dy = Dpoints(2,i)
        dz = Dpoints(3,i)
        Dbas(1,DER_DERIV3D_Y,i)=2.0_DP*COB(1,2)*dy+COB(1,4)*dx+COB(1,6)*dz+COB(1,8)
        Dbas(2,DER_DERIV3D_Y,i)=2.0_DP*COB(2,2)*dy+COB(2,4)*dx+COB(2,6)*dz+COB(2,8)
        Dbas(3,DER_DERIV3D_Y,i)=2.0_DP*COB(3,2)*dy+COB(3,4)*dx+COB(3,6)*dz+COB(3,8)
        Dbas(4,DER_DERIV3D_Y,i)=2.0_DP*COB(4,2)*dy+COB(4,4)*dx+COB(4,6)*dz+COB(4,8)
        Dbas(5,DER_DERIV3D_Y,i)=2.0_DP*COB(5,2)*dy+COB(5,4)*dx+COB(5,6)*dz+COB(5,8)
        Dbas(6,DER_DERIV3D_Y,i)=2.0_DP*COB(6,2)*dy+COB(6,4)*dx+COB(6,6)*dz+COB(6,8)
      END DO
    END IF
    
    ! Z-derivatives
    IF(BDER(DER_DERIV3D_Z)) THEN
      DO i=1, npoints
        dx = Dpoints(1,i)
        dy = Dpoints(2,i)
        dz = Dpoints(3,i)
        Dbas(1,DER_DERIV3D_Z,i)=2.0_DP*COB(1,3)*dz+COB(1,5)*dx+COB(1,6)*dy+COB(1,9)
        Dbas(2,DER_DERIV3D_Z,i)=2.0_DP*COB(2,3)*dz+COB(2,5)*dx+COB(2,6)*dy+COB(2,9)
        Dbas(3,DER_DERIV3D_Z,i)=2.0_DP*COB(3,3)*dz+COB(3,5)*dx+COB(3,6)*dy+COB(3,9)
        Dbas(4,DER_DERIV3D_Z,i)=2.0_DP*COB(4,3)*dz+COB(4,5)*dx+COB(4,6)*dy+COB(4,9)
        Dbas(5,DER_DERIV3D_Z,i)=2.0_DP*COB(5,3)*dz+COB(5,5)*dx+COB(5,6)*dy+COB(5,9)
        Dbas(6,DER_DERIV3D_Z,i)=2.0_DP*COB(6,3)*dz+COB(6,5)*dx+COB(6,6)*dy+COB(6,9)
      END DO
    END IF
  
  ! That's it
  
  END SUBROUTINE

  !************************************************************************
  
!<subroutine>  

  PURE SUBROUTINE elem_EM30_3D_sim (ieltyp, Dcoords, Djac, Ddetj, &
                                    Bder, Dbas, npoints, nelements, Dpoints)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier. Must be =EL_EM30_3D.
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
  REAL(DP), PARAMETER :: Q1 =   0.25_DP
  REAL(DP), PARAMETER :: Q2 =   1.0_DP /   3.0_DP
  REAL(DP), PARAMETER :: Q3 =   0.2_DP
  REAL(DP), PARAMETER :: Q4 =   0.6_DP
  REAL(DP), PARAMETER :: Q8 =  -9.0_DP /  16.0_DP
  REAL(DP), PARAMETER :: Q9 = 100.0_DP / 192.0_DP
  
  REAL(DP), DIMENSION(3,8) :: P                       ! cubature points
  REAL(DP) :: PP1,PP2,PP3,PP4,PP5,PS1,PS2,PQ1,PQ2,PQ3 ! cubature weights
  REAL(DP), DIMENSION(6,6) :: A, B                    ! coefficient matrices
  REAL(DP) :: CA1,CA2,CA3,CB1,CB2,CB3,CC1,&           ! auxiliary coefficients
              CC2,CC3,CD1,CD2,CD3,CD4,CD5,&
              CD6,CE1,CE2,CE3,CE4,CE5,CE6
  REAL(DP), DIMENSION(6,10) :: COB                    ! monomial coefficients
  REAL(DP), DIMENSION(3,4,6) :: V                     ! face corner vertices
  REAL(DP), DIMENSION(2:6,8) :: F                     ! transformed monomials
  REAL(DP) :: dx,dy,dz
  INTEGER :: i,j,k,l
  
  ! Loop over all elements
  DO j=1, nelements
  
    ! Initialise the vertex-coordinates
    V(1:3,1,1) = Dcoords(1:3,1,j)
    V(1:3,2,1) = Dcoords(1:3,2,j)
    V(1:3,3,1) = Dcoords(1:3,3,j)
    V(1:3,4,1) = Dcoords(1:3,4,j)
    V(1:3,1,2) = Dcoords(1:3,1,j)
    V(1:3,2,2) = Dcoords(1:3,2,j)
    V(1:3,3,2) = Dcoords(1:3,6,j)
    V(1:3,4,2) = Dcoords(1:3,5,j)
    V(1:3,1,3) = Dcoords(1:3,2,j)
    V(1:3,2,3) = Dcoords(1:3,3,j)
    V(1:3,3,3) = Dcoords(1:3,7,j)
    V(1:3,4,3) = Dcoords(1:3,6,j)
    V(1:3,1,4) = Dcoords(1:3,3,j)
    V(1:3,2,4) = Dcoords(1:3,4,j)
    V(1:3,3,4) = Dcoords(1:3,8,j)
    V(1:3,4,4) = Dcoords(1:3,7,j)
    V(1:3,1,5) = Dcoords(1:3,4,j)
    V(1:3,2,5) = Dcoords(1:3,1,j)
    V(1:3,3,5) = Dcoords(1:3,5,j)
    V(1:3,4,5) = Dcoords(1:3,8,j)
    V(1:3,1,6) = Dcoords(1:3,5,j)
    V(1:3,2,6) = Dcoords(1:3,6,j)
    V(1:3,3,6) = Dcoords(1:3,7,j)
    V(1:3,4,6) = Dcoords(1:3,8,j)

    ! Calculate the coefficients of the Fi
    CA1 = Q1*((Dcoords(1,1,j)+Dcoords(1,2,j)+Dcoords(1,3,j)+Dcoords(1,4,j)) &
             -(Dcoords(1,5,j)+Dcoords(1,6,j)+Dcoords(1,7,j)+Dcoords(1,8,j)))
    CB1 = Q1*((Dcoords(2,1,j)+Dcoords(2,2,j)+Dcoords(2,3,j)+Dcoords(2,4,j)) &
             -(Dcoords(2,5,j)+Dcoords(2,6,j)+Dcoords(2,7,j)+Dcoords(2,8,j)))
    CC1 = Q1*((Dcoords(3,1,j)+Dcoords(3,2,j)+Dcoords(3,3,j)+Dcoords(3,4,j)) &
             -(Dcoords(3,5,j)+Dcoords(3,6,j)+Dcoords(3,7,j)+Dcoords(3,8,j)))
    CA2 = Q1*((Dcoords(1,1,j)+Dcoords(1,2,j)+Dcoords(1,6,j)+Dcoords(1,5,j)) &
             -(Dcoords(1,3,j)+Dcoords(1,7,j)+Dcoords(1,8,j)+Dcoords(1,4,j)))
    CB2 = Q1*((Dcoords(2,1,j)+Dcoords(2,2,j)+Dcoords(2,6,j)+Dcoords(2,5,j)) &
             -(Dcoords(2,3,j)+Dcoords(2,7,j)+Dcoords(2,8,j)+Dcoords(2,4,j)))
    CC2 = Q1*((Dcoords(3,1,j)+Dcoords(3,2,j)+Dcoords(3,6,j)+Dcoords(3,5,j)) &
             -(Dcoords(3,3,j)+Dcoords(3,7,j)+Dcoords(3,8,j)+Dcoords(3,4,j)))
    CA3 = Q1*((Dcoords(1,2,j)+Dcoords(1,3,j)+Dcoords(1,7,j)+Dcoords(1,6,j)) &
             -(Dcoords(1,1,j)+Dcoords(1,4,j)+Dcoords(1,8,j)+Dcoords(1,5,j)))
    CB3 = Q1*((Dcoords(2,2,j)+Dcoords(2,3,j)+Dcoords(2,7,j)+Dcoords(2,6,j)) &
             -(Dcoords(2,1,j)+Dcoords(2,4,j)+Dcoords(2,8,j)+Dcoords(2,5,j)))
    CC3 = Q1*((Dcoords(3,2,j)+Dcoords(3,3,j)+Dcoords(3,7,j)+Dcoords(3,6,j)) &
             -(Dcoords(3,1,j)+Dcoords(3,4,j)+Dcoords(3,8,j)+Dcoords(3,5,j)))
    CD1 = CA1**2 - CA2**2
    CD2 = CB1**2 - CB2**2
    CD3 = CC1**2 - CC2**2
    CD4 = 2.0_DP*(CA1*CB1 - CA2*CB2)
    CD5 = 2.0_DP*(CA1*CC1 - CA2*CC2)
    CD6 = 2.0_DP*(CB1*CC1 - CB2*CC2)
    CE1 = CA2**2 - CA3**2
    CE2 = CB2**2 - CB3**2
    CE3 = CC2**2 - CC3**2
    CE4 = 2.0_DP*(CA2*CB2 - CA3*CB3)
    CE5 = 2.0_DP*(CA2*CC2 - CA3*CC3)
    CE6 = 2.0_DP*(CB2*CC2 - CB3*CC3)

    DO k=1,6
      ! Calculate the eight cubature points
      P(1,1)=Q2*(V(1,1,k)+V(1,2,k)+V(1,3,k))
      P(2,1)=Q2*(V(2,1,k)+V(2,2,k)+V(2,3,k))
      P(3,1)=Q2*(V(3,1,k)+V(3,2,k)+V(3,3,k))
      P(1,2)=Q4*V(1,1,k)+Q3*V(1,2,k)+Q3*V(1,3,k)
      P(2,2)=Q4*V(2,1,k)+Q3*V(2,2,k)+Q3*V(2,3,k)
      P(3,2)=Q4*V(3,1,k)+Q3*V(3,2,k)+Q3*V(3,3,k)
      P(1,3)=Q3*V(1,1,k)+Q4*V(1,2,k)+Q3*V(1,3,k)
      P(2,3)=Q3*V(2,1,k)+Q4*V(2,2,k)+Q3*V(2,3,k)
      P(3,3)=Q3*V(3,1,k)+Q4*V(3,2,k)+Q3*V(3,3,k)
      P(1,4)=Q3*V(1,1,k)+Q3*V(1,2,k)+Q4*V(1,3,k)
      P(2,4)=Q3*V(2,1,k)+Q3*V(2,2,k)+Q4*V(2,3,k)
      P(3,4)=Q3*V(3,1,k)+Q3*V(3,2,k)+Q4*V(3,3,k)
      P(1,5)=Q2*(V(1,1,k)+V(1,3,k)+V(1,4,k))
      P(2,5)=Q2*(V(2,1,k)+V(2,3,k)+V(2,4,k))
      P(3,5)=Q2*(V(3,1,k)+V(3,3,k)+V(3,4,k))
      P(1,6)=Q4*V(1,1,k)+Q3*V(1,3,k)+Q3*V(1,4,k)
      P(2,6)=Q4*V(2,1,k)+Q3*V(2,3,k)+Q3*V(2,4,k)
      P(3,6)=Q4*V(3,1,k)+Q3*V(3,3,k)+Q3*V(3,4,k)
      P(1,7)=Q3*V(1,1,k)+Q4*V(1,3,k)+Q3*V(1,4,k)
      P(2,7)=Q3*V(2,1,k)+Q4*V(2,3,k)+Q3*V(2,4,k)
      P(3,7)=Q3*V(3,1,k)+Q4*V(3,3,k)+Q3*V(3,4,k)
      P(1,8)=Q3*V(1,1,k)+Q3*V(1,3,k)+Q4*V(1,4,k)
      P(2,8)=Q3*V(2,1,k)+Q3*V(2,3,k)+Q4*V(2,4,k)
      P(3,8)=Q3*V(3,1,k)+Q3*V(3,3,k)+Q4*V(3,4,k)
      
      ! And calculate the weights for the cubature formula
      PP1=SQRT((V(1,1,k)-V(1,2,k))**2+(V(2,1,k)-V(2,2,k))**2+&
          (V(3,1,k)-V(3,2,k))**2)
      PP2=SQRT((V(1,2,k)-V(1,3,k))**2+(V(2,2,k)-V(2,3,k))**2+&
          (V(3,2,k)-V(3,3,k))**2)
      PP3=SQRT((V(1,1,k)-V(1,3,k))**2+(V(2,1,k)-V(2,3,k))**2+&
          (V(3,1,k)-V(3,3,k))**2)
      PP4=SQRT((V(1,3,k)-V(1,4,k))**2+(V(2,3,k)-V(2,4,k))**2+&
          (V(3,3,k)-V(3,4,k))**2)
      PP5=SQRT((V(1,1,k)-V(1,4,k))**2+(V(2,1,k)-V(2,4,k))**2+&
          (V(3,1,k)-V(3,4,k))**2)
      PS1=(PP1+PP2+PP3)*0.5_DP
      PS2=(PP3+PP4+PP5)*0.5_DP
      PQ1=SQRT(PS1*(PS1-PP1)*(PS1-PP2)*(PS1-PP3))
      PQ2=SQRT(PS2*(PS2-PP3)*(PS2-PP4)*(PS2-PP5))
      PQ3=1.0_DP / (PQ1 + PQ2)
      
      ! Evalute the F2..F6 in the eight cubature points
      ! Remember: F1 = 1
      DO l=1,8
        F(2,l) = CA1*P(1,l) + CB1*P(2,l) + CC1*P(3,l)
        F(3,l) = CA2*P(1,l) + CB2*P(2,l) + CC2*P(3,l)
        F(4,l) = CA3*P(1,l) + CB3*P(2,l) + CC3*P(3,l)
        F(5,l) = CD1*P(1,l)**2 + CD2*P(2,l)**2 + CD3*P(3,l)**2 &
               + CD4*P(1,l)*P(2,l) + CD5*P(1,l)*P(3,l) + CD6*P(2,l)*P(3,l)
        F(6,l) = CE1*P(1,l)**2 + CE2*P(2,l)**2 + CE3*P(3,l)**2 &
               + CE4*P(1,l)*P(2,l) + CE5*P(1,l)*P(3,l) + CE6*P(2,l)*P(3,l)
      END DO

      ! Build up the matrix entries
      A(k,1)= 1.0_DP
      DO l=2,6
        A(k,l) = PQ3*(PQ1*(Q8*F(l,1)+Q9*(F(l,2)+F(l,3)+F(l,4))) &
                     +PQ2*(Q8*F(l,5)+Q9*(F(l,6)+F(l,7)+F(l,8))))
      END DO
               
    END DO
    
    ! Now we have the coeffienct matrix - so invert it
    CALL mprim_invert6x6MatrixDirectDble(A,B)

    ! Transform coefficiencts into monomial base
    DO k=1,6
      COB(k,1) = B(6,k)*CE1 + B(5,k)*CD1
      COB(k,2) = B(6,k)*CE2 + B(5,k)*CD2
      COB(k,3) = B(6,k)*CE3 + B(5,k)*CD3
      COB(k,4) = B(6,k)*CE4 + B(5,k)*CD4
      COB(k,5) = B(6,k)*CE5 + B(5,k)*CD5
      COB(k,6) = B(6,k)*CE6 + B(5,k)*CD6
      COB(k,7) = B(4,k)*CA3 + B(3,k)*CA2 + B(2,k)*CA1
      COB(k,8) = B(4,k)*CB3 + B(3,k)*CB2 + B(2,k)*CB1
      COB(k,9) = B(4,k)*CC3 + B(3,k)*CC2 + B(2,k)*CC1
      COB(k,10)= B(1,k)
    END DO
    
    ! Function values
    IF (BDER(DER_FUNC3D)) THEN
      DO i=1, npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)
        dz = Dpoints(3,i,j)
        Dbas(1,DER_FUNC3D,i,j)= COB(1,1)*dx**2+COB(1,2)*dy**2+COB(1,3)*dz**2&
                 +COB(1,4)*dx*dy+COB(1,5)*dx*dz+COB(1,6)*dy*dz&
                 +COB(1,7)*dx+COB(1,8)*dy+COB(1,9)*dz+COB(1,10)
        Dbas(2,DER_FUNC3D,i,j)= COB(2,1)*dx**2+COB(2,2)*dy**2+COB(2,3)*dz**2&
                 +COB(2,4)*dx*dy+COB(2,5)*dx*dz+COB(2,6)*dy*dz&
                 +COB(2,7)*dx+COB(2,8)*dy+COB(2,9)*dz+COB(2,10)
        Dbas(3,DER_FUNC3D,i,j)= COB(3,1)*dx**2+COB(3,2)*dy**2+COB(3,3)*dz**2&
                 +COB(3,4)*dx*dy+COB(3,5)*dx*dz+COB(3,6)*dy*dz&
                 +COB(3,7)*dx+COB(3,8)*dy+COB(3,9)*dz+COB(3,10)
        Dbas(4,DER_FUNC3D,i,j)= COB(4,1)*dx**2+COB(4,2)*dy**2+COB(4,3)*dz**2&
                 +COB(4,4)*dx*dy+COB(4,5)*dx*dz+COB(4,6)*dy*dz&
                 +COB(4,7)*dx+COB(4,8)*dy+COB(4,9)*dz+COB(4,10)
        Dbas(5,DER_FUNC3D,i,j)= COB(5,1)*dx**2+COB(5,2)*dy**2+COB(5,3)*dz**2&
                 +COB(5,4)*dx*dy+COB(5,5)*dx*dz+COB(5,6)*dy*dz&
                 +COB(5,7)*dx+COB(5,8)*dy+COB(5,9)*dz+COB(5,10)
        Dbas(6,DER_FUNC3D,i,j)= COB(6,1)*dx**2+COB(6,2)*dy**2+COB(6,3)*dz**2&
                 +COB(6,4)*dx*dy+COB(6,5)*dx*dz+COB(6,6)*dy*dz&
                 +COB(6,7)*dx+COB(6,8)*dy+COB(6,9)*dz+COB(6,10)
      END DO
    END IF

    ! X-derivatives
    IF(BDER(DER_DERIV3D_X)) THEN
      DO i=1, npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)
        dz = Dpoints(3,i,j)
        Dbas(1,DER_DERIV3D_X,i,j)=2.0_DP*COB(1,1)*dx+COB(1,4)*dy+COB(1,5)*dz+COB(1,7)
        Dbas(2,DER_DERIV3D_X,i,j)=2.0_DP*COB(2,1)*dx+COB(2,4)*dy+COB(2,5)*dz+COB(2,7)
        Dbas(3,DER_DERIV3D_X,i,j)=2.0_DP*COB(3,1)*dx+COB(3,4)*dy+COB(3,5)*dz+COB(3,7)
        Dbas(4,DER_DERIV3D_X,i,j)=2.0_DP*COB(4,1)*dx+COB(4,4)*dy+COB(4,5)*dz+COB(4,7)
        Dbas(5,DER_DERIV3D_X,i,j)=2.0_DP*COB(5,1)*dx+COB(5,4)*dy+COB(5,5)*dz+COB(5,7)
        Dbas(6,DER_DERIV3D_X,i,j)=2.0_DP*COB(6,1)*dx+COB(6,4)*dy+COB(6,5)*dz+COB(6,7)
      END DO
    END IF

    ! Y-derivatives
    IF(BDER(DER_DERIV3D_Y)) THEN
      DO i=1, npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)
        dz = Dpoints(3,i,j)
        Dbas(1,DER_DERIV3D_Y,i,j)=2.0_DP*COB(1,2)*dy+COB(1,4)*dx+COB(1,6)*dz+COB(1,8)
        Dbas(2,DER_DERIV3D_Y,i,j)=2.0_DP*COB(2,2)*dy+COB(2,4)*dx+COB(2,6)*dz+COB(2,8)
        Dbas(3,DER_DERIV3D_Y,i,j)=2.0_DP*COB(3,2)*dy+COB(3,4)*dx+COB(3,6)*dz+COB(3,8)
        Dbas(4,DER_DERIV3D_Y,i,j)=2.0_DP*COB(4,2)*dy+COB(4,4)*dx+COB(4,6)*dz+COB(4,8)
        Dbas(5,DER_DERIV3D_Y,i,j)=2.0_DP*COB(5,2)*dy+COB(5,4)*dx+COB(5,6)*dz+COB(5,8)
        Dbas(6,DER_DERIV3D_Y,i,j)=2.0_DP*COB(6,2)*dy+COB(6,4)*dx+COB(6,6)*dz+COB(6,8)
      END DO
    END IF
    
    ! Z-derivatives
    IF(BDER(DER_DERIV3D_Z)) THEN
      DO i=1, npoints
        dx = Dpoints(1,i,j)
        dy = Dpoints(2,i,j)
        dz = Dpoints(3,i,j)
        Dbas(1,DER_DERIV3D_Z,i,j)=2.0_DP*COB(1,3)*dz+COB(1,5)*dx+COB(1,6)*dy+COB(1,9)
        Dbas(2,DER_DERIV3D_Z,i,j)=2.0_DP*COB(2,3)*dz+COB(2,5)*dx+COB(2,6)*dy+COB(2,9)
        Dbas(3,DER_DERIV3D_Z,i,j)=2.0_DP*COB(3,3)*dz+COB(3,5)*dx+COB(3,6)*dy+COB(3,9)
        Dbas(4,DER_DERIV3D_Z,i,j)=2.0_DP*COB(4,3)*dz+COB(4,5)*dx+COB(4,6)*dy+COB(4,9)
        Dbas(5,DER_DERIV3D_Z,i,j)=2.0_DP*COB(5,3)*dz+COB(5,5)*dx+COB(5,6)*dy+COB(5,9)
        Dbas(6,DER_DERIV3D_Z,i,j)=2.0_DP*COB(6,3)*dz+COB(6,5)*dx+COB(6,6)*dy+COB(6,9)
      END DO
    END IF

  END DO ! element loop
  
  ! That's it
  END SUBROUTINE
    
END MODULE 
