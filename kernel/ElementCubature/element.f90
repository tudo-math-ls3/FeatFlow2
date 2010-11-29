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
!# a finite element:
!#
!# 1.) elem_igetNDofLoc
!#     -> determines the number of local degrees of freedom for a finite 
!#        element
!#
!# 2.) elem_igetNVE
!#     -> get the number of vertices in the element primitive/element shape
!#        of the Finite Element (3=triangular, 4=quad)
!#
!# 3.) elem_igetCoordSystem
!#     -> get the type of coordinate system, a Finite Element uses
!#
!# 4.) elem_igetTrafoType
!#     -> Determine the type of transformation from the reference element
!#        to the real element.
!#
!# 5.) elem_igetDimension
!#     -> Determine the dimension of an element; 1D, 2D, 3D,...
!#
!# 6.) elem_isNonparametric
!#     -> Check whether an element is parametric or nonparametric
!#
!# 7.) elem_getPrimaryElement
!#     -> Returns the element identifier of the 'primary' element without
!#        any subtype information, which identifies an element family.
!#
!# 8.) elem_generic 
!#     -> Realises a generic element which can be used to evaluate a finite 
!#        element depending on its element identifier - in contrast to the 
!#        standard evaluation routines, which ignore the element quantifier 
!#        as they 'know' what they are...
!#
!# 9.) elem_generic_mult
!#     -> The multiple-point-evaluation routine for a generic element.
!#
!# 10.) elem_generic_sim
!#     -> The multiple-point/element-evaluation routine for a generic element.
!#
!# 11.) elem_getEvaluationTag
!#     -> Returns for an element it is 'evaluation tag'. This is a bitfield
!#        which encodes what information an element needs if one wants
!#        to evaluate it in a point.
!#
!# 12.) elem_igetID
!#      -> Convert a string with an element identifier to the corresponding
!#         element id.
!#
!#  FAQ - Some explainations  \\
!# -------------------------- \\
!# 1.) What is an element and what are Degrees Of Freedoms (DOF`s)?
!#
!#   Imagine a triangulation of a domain with quadrilateral polygons, e.g.
!#   rectangles or even squares. Take for example four of these squares and
!#   form a large square, called a "patch", like
!#
!#   <verb>
!#       O-----O-----O
!#       |     |     |
!#       |     |     |
!#       O-----1-----O
!#       |     |     |
!#       |     |     |
!#       O-----O-----O
!#   </verb>
!#
!#   A point-based finite element basis function like <tex>$Q_1$</tex> is a function
!#   that is =1 at the center of this patch and =0 in all the other corner
!#   points. On each of the squares, the function is a polynom, which
!#   coefficients are chosen in such a way, that it is =1 in the center and
!#   =0 in the other corner points. Outside of this patch, the function
!#   is identically =0.
!#
!#   <tex>
!#   If you call one such a function $\phi_i$, the whole function 
!#   <tex>$u:R^2 \to R$</tex> is then defined as a sum
!#
!#     $$ u(x,y) = \sum_i  u_i  \phi_i(x,y) $$
!#
!#   where the coefficients $u_i$ can be modified to change the whole
!#   function. This set ${u_i}$ are called "Degrees of Freedom" of the FE
!#   function $u$. If the $\phi_i$ are chosen to be =1 in the center
!#   and =0 in the other corners (like in $P_1$ and $Q_1$, by `lucky 
!#   chance` the $u_i$ are exactly the values of $u$ in these points --
!#   $u_i = u(x_i,y_i)$.
!#   </tex>
!#
!#   An "element" in the sense of this library focuses on one single
!#   polygon, not on the whole patch.
!#
!# 2.) What does the variable 'celement' mean, which must be passed to
!#   the element subroutines?
!#
!#   This variable identifies the type of the element. For example
!#   celement=EL_Q1 identifies a Q1-element.
!#
!# 3.) I wrote a subroutine which works for a special element family. More 
!#   precisely for example, my routine works for all variants of <tex>$\tilde Q_1$</tex>.
!#   Do I really have to check all element types at the beginning of the
!#   routine? That would be legthy!?! Like in
!#
!# <code>
!#    subroutine abc (celement)
!#    ...
!#      if ( (celement .ne. EL_E030) .and. (celement .ne. EL_EM30) ...) then
!#        call sys_halt()
!#      end if
!# </code>
!#
!#   No you do not need to. The simplest method for such element families:
!#   Use elem_getPrimaryElement! This returns an element identifier
!#   identifying all elements 'of the same type'. The above can be
!#   shortened this way to:
!#
!# <code>
!#    subroutine abc (celement)
!#    ...
!#      if (elem_getPrimaryElement(celement) .ne. EL_Q1T) then
!#        call sys_halt()
!#      end if
!# </code>
!#
!#   Which directly describes all <tex>$\tilde Q_1$</tex> elements. If you additionally
!#   want to check it whether it is parametric/nonparametric, use an
!#   additional check with elem_isNonparametric!
!#
!#   Of course, for 'standard' elements, elem_getPrimaryElement does nothing.
!#   So writing for <tex>$Q_1$</tex> for example,
!#
!# <code>
!#      if (elem_getPrimaryElement(celement) .ne. EL_Q1) then
!# </code>
!#
!#   is exactly the same as
!#
!# <code>
!#      if (celement .ne. EL_Q1) then
!# </code>
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
!# <verb>
!# Bit | 31 ... 24 23 ... 16 | 15                | 14 ... 10 | 9   8 | 7 ............ 0|
!# -------------------------------------------------------------------------------------
!#     |         ****        | =1: nonparametric | unused    |dimens |  Element number |
!#
!#   Bits 0..7   specifies the element number/family (1=P1, 11=Q1,...).
!#   Bits 8+ 9   encodes the dimension of the element. =1: 1D, =2: 2D, =3: 3D.
!#   Bit    15   specifies whether an element is nonparametric (1) or parametric (0).
!#   Bits 16..31 encodes special variants of elements. This is used only for some special
!#               type elements:
!#               Q1T: Bit 16 =1: element nonconformal, integral mean value based, 
!#                           =0: element conformal, midpoint value based
!#                    Bit 17:=0: Standard handling
!#                           =1: For nonparametric element: Do not use pivoting when 
!#                               solving local 4x4 / 6x6 systems; faster but less
!#                               stable on cruel elements.
!#                    Bit 18:=0: Standard handling
!#                           =1: No scaling of the local coodinate system for 
!#                               nonparametric <tex>$\tilde Q_1$</tex> element. Slightly faster but 
!#                               less stable on cruel elements.
!#               Q1HN: Bit 16/17: These three bits encode the local number of the edge
!#                               where there is a hanging node.
!#                               =0: local edge 1, =1: local edge 2,...
!#               P2ISO: Bit 16/17: These three bits encode how many edges of the <tex>$P_2$</tex>
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
!#               Q2ISO: Bit 16/17: These three bits encode how many edges of the <tex>$Q_2$</tex>
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
!# </verb>
!#
!#   To obtain the actual element identifier, one must mask out all bits 
!#   except for bit 0..7, i.e. to check whether an element is Q1~, one can use
!#   elem_getPrimaryElement, which masks out all unimportant bits with IAND:
!#
!# <code>
!#     if ( elem_getPrimaryElement(celement) .eq. EL_Q1T ) then ...
!# </code>
!#
!#   or to check all variants
!#
!# <code>
!#     if ( (celement .eq. EL_E030) .or. (celement .eq. EL_EM30) .or. ...) then ...
!# </code>
!# 
!#   When it is clear that it is a <tex>$\tilde Q_1$</tex> element, one can have a closer
!#   look which variant it is -- if this is necessary.
!#
!#   Final note: For implementational reasons, the bits 8+9+16..31 coincide --
!#   depending on the element -- with the ID of the transformation formula
!#   from the reference to the real element! This allows an easier and more
!#   transparent implementation of isoparametric elements!
!#
!# </purpose>
!##############################################################################

module element

  use fsystem
  use genoutput
  use basicgeometry
  use derivatives
  use transformation
  use mprimitives
  
  use elementbase
  use element_line1d
  use element_tri2d
  use element_quad2d
  use element_tetra3d
  use element_hexa3d
  use element_pyra3d
  use element_prism3d

  implicit none
  
  private
  
  public :: elem_igetID
  public :: elem_igetNDofLoc
  public :: elem_igetNDofLocAssignment
  public :: elem_igetNVE
  public :: elem_igetCoordSystem
  public :: elem_igetTrafoType
  public :: elem_igetDimension
  public :: elem_getMaxDerivative
  public :: elem_getEvaluationTag
  public :: elem_isNonparametric
  public :: elem_getPrimaryElement
  public :: elem_igetShape
  public :: elem_generic1
  public :: elem_generic2
  public :: elem_generic_mult
  public :: elem_generic_sim1
  public :: elem_generic_sim2
  public :: elem_generic
  public :: elem_generic_sim
  
  ! Public entities imported from elementbase.f90
  public :: t_evalElement
  public :: t_evalElementSet
  public :: EL_EVLTAG_COORDS, EL_EVLTAG_REFPOINTS, EL_EVLTAG_REALPOINTS, &
            EL_EVLTAG_JAC, EL_EVLTAG_DETJ, EL_EVLTAG_TWISTIDX
  
!<constants>

!<constantblock description="Internal constants for element ID bitfield.">
  ! Bitmasks for dimension; coincides on purpose with TRAFO_DIM_DIMENSION!
  integer(I32), parameter, public :: EL_DIMENSION = 2**8 + 2**9
  
  ! 1D element; coincides on purpose with TRAFO_DIM_1D!
  integer(I32), parameter, public :: EL_1D = ishft(NDIM1D,8)

  ! 2D element; coincides on purpose with TRAFO_DIM_2D!
  integer(I32), parameter, public :: EL_2D = ishft(NDIM2D,8)
  
  ! 3D element; coincides on purpose with TRAFO_DIM_3D!
  integer(I32), parameter, public :: EL_3D = ishft(NDIM3D,8)
  
  ! Bitmask for element number including dimension
  integer(I32), parameter, public :: EL_ELNRMASK = 255 + EL_DIMENSION
  
  ! Bitmask specifying a nonparametric element
  integer(I32), parameter, public :: EL_NONPARAMETRIC = 2**15
!</constantblock>

!<constantblock description="Element identifiers for 1D elements">

  ! ID of constant conforming line FE, P0
  integer(I32), parameter, public :: EL_P0_1D   = EL_1D + 0
  integer(I32), parameter, public :: EL_E000_1D = EL_P0_1D
  
  ! ID of linear conforming line FE, P1
  integer(I32), parameter, public :: EL_P1_1D   = EL_1D + 1
  integer(I32), parameter, public :: EL_E001_1D = EL_P1_1D
  
  ! ID of quadratic conforming line FE, P2
  integer(I32), parameter, public :: EL_P2_1D   = EL_1D + 2
  integer(I32), parameter, public :: EL_E002_1D = EL_P2_1D
  
  ! ID of cubic conforming line FE, 3,1-Spline
  integer(I32), parameter, public :: EL_S31_1D  = EL_1D + 17
  
  ! ID of conforming line FE, Pn, 1 <= n <= 256
  integer(I32), parameter, public :: EL_PN_1D   = EL_1D + 128
  
  ! Discontinuous Galerkin taylor basis element - constant.
  integer(I32), parameter, public :: EL_DG_T0_1D   = EL_1D + 30
  
  ! Discontinuous Galerkin taylor basis element - linear.
  integer(I32), parameter, public :: EL_DG_T1_1D   = EL_1D + 31

  ! Discontinuous Galerkin taylor basis element - quadratic.
  integer(I32), parameter, public :: EL_DG_T2_1D   = EL_1D + 32
  
!</constantblock>
  
!<constantblock description="Element identifiers for 2D elements">

  ! unspecified element
  integer(I32), parameter, public :: EL_UNDEFINED = -1

  ! ID of constant conforming triangular FE, P0
  integer(I32), parameter, public :: EL_P0      = EL_2D + 0
  integer(I32), parameter, public :: EL_E000    = EL_P0
  integer(I32), parameter, public :: EL_P0_2D   = EL_P0
  integer(I32), parameter, public :: EL_E000_2D = EL_P0

  ! ID of linear conforming triangular FE, P1
  integer(I32), parameter, public :: EL_P1      = EL_2D + 1
  integer(I32), parameter, public :: EL_E001    = EL_P1
  integer(I32), parameter, public :: EL_P1_2D   = EL_P1
  integer(I32), parameter, public :: EL_E001_2D = EL_P1

  ! ID of quadratic conforming triangular FE, P2
  integer(I32), parameter, public :: EL_P2      = EL_2D + 2
  integer(I32), parameter, public :: EL_E002    = EL_P2
  integer(I32), parameter, public :: EL_P2_2D   = EL_P2
  integer(I32), parameter, public :: EL_E002_2D = EL_P2

  ! ID of cubic conforming triangular FE, P3
  integer(I32), parameter, public :: EL_P3      = EL_2D + 3
  integer(I32), parameter, public :: EL_E003    = EL_P3
  integer(I32), parameter, public :: EL_P3_2D   = EL_P3
  integer(I32), parameter, public :: EL_E003_2D = EL_P3

  ! ID for rotated linear <tex>$\tilde P_1$</tex> element (Crouzeix-Raviart)
  integer(I32), parameter, public :: EL_P1T     = EL_2D + 20
  integer(I32), parameter, public :: EL_P1T_2D  = EL_P1T

  ! ID of constant conforming quadrilateral FE, Q0
  integer(I32), parameter, public :: EL_Q0      = EL_2D + 10
  integer(I32), parameter, public :: EL_E010    = EL_Q0
  integer(I32), parameter, public :: EL_Q0_2D   = EL_Q0
  integer(I32), parameter, public :: EL_E010_2D = EL_Q0

  ! ID of bilinear conforming quadrilateral FE, Q1
  integer(I32), parameter, public :: EL_Q1      = EL_2D + 11
  integer(I32), parameter, public :: EL_E011    = EL_Q1 
  integer(I32), parameter, public :: EL_Q1_2D   = EL_Q1 
  integer(I32), parameter, public :: EL_E011_2D = EL_Q1 

  ! ID of biquadratic conforming quadrilateral FE, Q2
  integer(I32), parameter, public :: EL_Q2      = EL_2D + 13
  integer(I32), parameter, public :: EL_E013    = EL_Q2
  integer(I32), parameter, public :: EL_Q2_2D   = EL_Q2

  ! ID of bicubic conforming quadrilateral FE, Q3
  integer(I32), parameter, public :: EL_Q3      = EL_2D + 14
  integer(I32), parameter, public :: EL_E014    = EL_Q3
  integer(I32), parameter, public :: EL_Q3_2D   = EL_Q3
  
  ! ID of nonconforming parametric linear P1 element on a quadrilareral
  ! element, given by function value in the midpoint and the two
  ! derivatives.
  integer(I32), parameter, public :: EL_QP1     = EL_2D + 21
  integer(I32), parameter, public :: EL_QP1_2D  = EL_QP1
  
  ! QP1-element, nonparametric
  integer(I32), parameter, public :: EL_QP1NP    = EL_QP1 + EL_NONPARAMETRIC
  integer(I32), parameter, public :: EL_QP1NP_2D = EL_QP1NP

  ! QP1-element, nonparametric, direct on element
  integer(I32), parameter, public :: EL_QP1NPD   = EL_QP1 + EL_NONPARAMETRIC + 2**17
  integer(I32), parameter, public :: EL_QP1NPD_2D = EL_QP1NPD

  ! General rotated bilinear <tex>$\tilde Q_1$</tex> element, all variants (conformal, 
  ! nonconformal, parametric, nonparametric).
  ! Simplest variant is: parametric, edge midpoint-value based.
  integer(I32), parameter, public :: EL_Q1T     = EL_2D + 30
  integer(I32), parameter, public :: EL_Q1T_2D  = EL_Q1T

  ! General rotated bilinear <tex>$\tilde Q_1$</tex> element with bubble, 
  ! all variants (conformal, nonconformal, parametric, nonparametric).
  integer(I32), parameter, public :: EL_Q1TB    = EL_2D + 31
  integer(I32), parameter, public :: EL_Q1TB_2D = EL_Q1TB

  ! General rotated biquadratic <tex>$\tilde Q_2$</tex> element, all variants.
  integer(I32), parameter, public :: EL_Q2T     = EL_2D + 35
  integer(I32), parameter, public :: EL_Q2T_2D  = EL_Q2T

  ! General rotated biquadratic <tex>$\tilde Q_2$</tex> element with bubble, all variants.
  integer(I32), parameter, public :: EL_Q2TB    = EL_2D + 37
  integer(I32), parameter, public :: EL_Q2TB_2D = EL_Q2TB
  
  ! Quadrilateral <tex>$Q_1$</tex> element with one hanging node. In the property
  ! bitfield of the element identifier, the local number of the edge where there
  ! is a hanging node must be encoded! (i.e. one must add a corresponding
  ! value to the constant EL_Q1HN1 to get the actual element identifier!)
  integer(I32), parameter, public :: EL_Q1HN1 = EL_2D + 40
  
  ! Quadrilateral <tex>$Q_2$</tex> element with isoparametric mapping from the reference
  ! to the real element. In the property bitfield, one must set the corresponding
  ! bits to identify the edge that should map isoparametric!
  integer(I32), parameter, public :: EL_Q2ISO = EL_2D + 50
  
  ! Discontinuous Galerkin taylor basis element - constant.
  integer(I32), parameter, public :: EL_DG_T0_2D   = EL_2D + 60
  
  ! Discontinuous Galerkin taylor basis element - linear.
  integer(I32), parameter, public :: EL_DG_T1_2D   = EL_2D + 61

  ! Discontinuous Galerkin taylor basis element - quadratic.
  integer(I32), parameter, public :: EL_DG_T2_2D   = EL_2D + 62

!</constantblock>

!<constantblock description="Element identifiers for 3D elements">

  ! ID of constant conforming tetrahedral FE, P0
  integer(I32), parameter, public :: EL_P0_3D   = EL_3D + 0
  integer(I32), parameter, public :: EL_E000_3D = EL_P0_3D
  
  ! ID of linear conforming tetrahedral FE, P1
  integer(I32), parameter, public :: EL_P1_3D   = EL_3D + 1
  integer(I32), parameter, public :: EL_E001_3D = EL_P1_3D

  ! ID of quadratic conforming tetrahedral FE, P2
  integer(I32), parameter, public :: EL_P2_3D   = EL_3D + 2
  integer(I32), parameter, public :: EL_E002_3D = EL_P2_3D

  ! ID of constant conforming hexahedral FE, Q0
  integer(I32), parameter, public :: EL_Q0_3D   = EL_3D + 10
  integer(I32), parameter, public :: EL_E010_3D = EL_Q0_3D

  ! ID of trilinear conforming hexahedral FE, Q1
  integer(I32), parameter, public :: EL_Q1_3D   = EL_3D + 11
  integer(I32), parameter, public :: EL_E011_3D = EL_Q1_3D
  
  ! ID of triquadratic conforming hexahedral FE, Q2
  integer(I32), parameter, public :: EL_Q2_3D   = EL_3D + 13
  integer(I32), parameter, public :: EL_E013_3D = EL_Q2_3D

  ! ID of constant conforming pyramid FE, Y0 = Q0
  integer(I32), parameter, public :: EL_Y0_3D   = EL_3D + 60
  
  ! ID of sub-trilinear conforming pyramid FE, Y1 \subset Q1
  integer(I32), parameter, public :: EL_Y1_3D   = EL_3D + 61
  
  ! ID of constant conforming prism FE, R0 = Q0
  integer(I32), parameter, public :: EL_R0_3D   = EL_3D + 70
  
  ! ID of sub-trilinear conforming prism FE, R1 \subset Q1
  integer(I32), parameter, public :: EL_R1_3D   = EL_3D + 71

  ! ID of rotated trilinear non-conforming hexahedral FE, Q1~
  ! parametric, face-midpoint-value based
  integer(I32), parameter, public :: EL_Q1T_3D  = EL_3D + 30
  integer(I32), parameter, public :: EL_E031_3D = EL_Q1T_3D
  
  ! ID of rotated trilinear non-conforming hexahedral FE, Q1~
  ! parametric, integral mean value based
  integer(I32), parameter, public :: EL_E030_3D = EL_Q1T_3D + 2**16

  ! ID of rotated trilinear non-conforming hexahedral FE, Q1~
  ! non-parametric, integral mean value based
  integer(I32), parameter, public :: EL_EM30_3D = EL_Q1T_3D + EL_NONPARAMETRIC + 2**16
  
  ! ID of rotated trilinear non-conforming hexahedral FE, Q1~
  ! non-parametric, integral mean value based, new implementation
  integer(I32), parameter, public :: EL_EM30_NEW_3D = EL_Q1T_3D + EL_NONPARAMETRIC &
                                                      + 2**16 + 2**19

  ! ID of nonconforming quadrilateral FE, Q2~.
  integer(I32), parameter, public :: EL_Q2T_3D  = EL_3D + 50
  integer(I32), parameter, public :: EL_E050_3D = EL_Q2T_3D

  ! ID of nonconforming quadrilateral FE, Q2~, non-parametric
  integer(I32), parameter, public :: EL_EM50_3D = EL_Q2T_3D + EL_NONPARAMETRIC

  ! ID of discontinous parametric linear hexahedron FE, P1
  integer(I32), parameter, public :: EL_QP1_3D  = EL_3D + 21

  ! ID of discontinous parametric linear hexahedron FE, P1, nonparametric
  integer(I32), parameter, public :: EL_QP1NP_3D  = EL_QP1_3D + EL_NONPARAMETRIC

!</constantblock>

!<constantblock description="Special 2D element variants.">

  ! ID of bilinear conforming quadrilateral FE, Q1, non-parametric
  integer(I32), parameter, public :: EL_EM11    = EL_Q1 + EL_NONPARAMETRIC
  integer(I32), parameter, public :: EL_EM11_2D = EL_EM11
  
  ! ID of rotated linear nonconforming triangular FE, P1~, edge-midpoint based
  integer(I32), parameter, public :: EL_E020    = EL_P1T
  integer(I32), parameter, public :: EL_E020_2D = EL_E020

  ! ID of rotated bilinear conforming quadrilateral FE, Q1~, integral
  ! mean value based
  integer(I32), parameter, public :: EL_E030    = EL_Q1T + 2**16
  integer(I32), parameter, public :: EL_E030_2D = EL_E030

  ! ID of rotated bilinear conforming quadrilateral FE, Q1~ with bubble, integral
  ! mean value based
  integer(I32), parameter, public :: EL_EB30    = EL_Q1TB
  integer(I32), parameter, public :: EL_EB30_2D = EL_EB30

  ! ID of rotated bilinear conforming quadrilateral FE, Q1~, edge-midpoint based
  integer(I32), parameter, public :: EL_E031    = EL_Q1T
  integer(I32), parameter, public :: EL_E031_2D = EL_E031

  ! ID of rotated bilinear nonconforming quadrilateral FE, Q1~, integral
  ! mean value based
  integer(I32), parameter, public :: EL_EM30    = EL_Q1T + EL_NONPARAMETRIC+ 2**16
  integer(I32), parameter, public :: EL_EM30_2D = EL_EM30

  ! ID of rotated bilinear nonconforming quadrilateral FE, Q1~, integral
  ! mean value based; 'unpivoted' variant, solving local 4x4 systems directly 
  ! without pivoting. Faster but less stable.
  integer(I32), parameter, public :: EL_EM30_UNPIVOTED = EL_Q1T + EL_NONPARAMETRIC + 2**17 + 2**16
  
  ! ID of rotated bilinear nonconforming quadrilateral FE, Q1~, integral
  ! mean value based; 'unscaled' variant, which does not scale the local
  ! coordinate system on every element. Faster but less stable.
  integer(I32), parameter, public :: EL_EM30_UNSCALED = EL_Q1T + EL_NONPARAMETRIC + 2**18 + 2**16
  
  ! ID of rotated bilinear nonconforming quadrilateral FE, Q1~, integral
  ! mean value based; new interface implementations
  integer(I32), parameter, public :: EL_EM30_NEW = EL_Q1T + EL_NONPARAMETRIC + 2**19 + 2**16

  ! ID of rotated bilinear nonconforming quadrilateral FE, Q1~, edge-midpoint based
  integer(I32), parameter, public :: EL_EM31    = EL_Q1T + EL_NONPARAMETRIC
  integer(I32), parameter, public :: EL_EM31_2D = EL_EM31
  
  ! ID of rotated bilinear nonconforming quadrilateral FE, Q1~, edge-midpoint based
  ! new interface implementations
  integer(I32), parameter, public :: EL_EM31_NEW = EL_Q1T + EL_NONPARAMETRIC + 2**19
  
  ! ID of 'rotated bilinear enhanced' nonconforming quarilateral FE
  integer(I32), parameter, public :: EL_E032    = EL_Q1TB + 2**17
  integer(I32), parameter, public :: EL_E032_2D = EL_E032

  ! ID of rotated biquadratic nonconforming quadrilateral FE, Q2~.
  integer(I32), parameter, public :: EL_E050    = EL_Q2T
  integer(I32), parameter, public :: EL_E050_2D = EL_E050
  
  ! ID of rotated biquadratic nonconforming quadrilateral FE, Q2~,
  ! non-parametric version
  integer(I32), parameter, public :: EL_EM50    = EL_Q2T + EL_NONPARAMETRIC
  integer(I32), parameter, public :: EL_EM50_2D = EL_EM50

  ! ID of rotated biquadratic nonconforming quadrilateral FE, Q2~ with bubble.
  integer(I32), parameter, public :: EL_EB50    = EL_Q2TB
  integer(I32), parameter, public :: EL_EB50_2D = EL_EB50
  
  ! ID of <tex>$Q_2$</tex> element with hierarchical basis functions.
  ! WARNING: Do not use this element, as it is highly experimental and is not
  ! yet supported by the majority of the kernel routines!
  integer(I32), parameter, public :: EL_Q2H_2D = EL_Q2_2D + 2**16
  
  ! Isoparametric <tex>$Q_2$</tex> element with one edge mapped nonlinear from the reference
  ! to the real element. Additionally, one bit in the property bitfield must
  ! be set to identify the edge.
  integer(I32), parameter, public :: EL_Q2ISO1 = EL_Q2ISO 

  ! Isoparametric <tex>$Q_2$</tex> element with two edges mapped nonlinear from the reference
  ! to the real element. Additionally, two bits in the property bitfield must
  ! be set to identify the edges.
  integer(I32), parameter, public :: EL_Q2ISO2 = EL_Q2ISO + 2**16

  ! Isoparametric <tex>$Q_2$</tex> element with three edges mapped nonlinear from the reference
  ! to the real element. Additionally, three bits in the property bitfield must
  ! be set to identify the edges.
  integer(I32), parameter, public :: EL_Q2ISO3 = EL_Q2ISO + 2**17

  ! Isoparametric <tex>$Q_2$</tex> element with four edges mapped nonlinear from the reference
  ! to the real element. Additionally, four bits in the property bitfield must
  ! be set to identify the edges.
  integer(I32), parameter, public :: EL_Q2ISO4 = EL_Q2ISO + 2**16 + 2**17


  ! Discontinuous Galerkin element based on P1, 2D.
  integer(I32), parameter, public :: EL_DG_P1_2D   = EL_P1 + 101

!</constantblock>

!<constantblock description="maximal values">

  ! DEPRECATED: Maximum number of basic functions = maximum number of
  ! local DOF`s per element.
  ! Do not use this constant anymore - determine the number of local basis
  ! functions dynamically using the 'elem_igetNDofLoc' routine!
  integer, parameter, public :: EL_MAXNBAS = 27
  
  ! Maximum number of derivatives. Corresponds to DER_MAXNDER.
  integer, parameter, public :: EL_MAXNDER = DER_MAXNDER

  ! Maximum number of DOFs per vertice.
  integer, parameter, public :: EL_MAXNDOF_PER_VERT = 2
  
  ! Maximum number of DOFs per edge.
  integer, parameter, public :: EL_MAXNDOF_PER_EDGE = 2
  
  ! Maximum number of DOFs per face.
  integer, parameter, public :: EL_MAXNDOF_PER_FACE = 3
  
  ! Maximum number of DOFs per element.
  integer, parameter, public :: EL_MAXNDOF_PER_ELEM = 6

  ! Number of entries in the Jacobian matrix, defining the mapping between
  ! the reference element and the real element. 1x1-matrix=1 elements
  integer, parameter, public :: EL_NJACENTRIES1D = 1

  ! Number of entries in the Jacobian matrix, defining the mapping between
  ! the reference element and the real element. 2x2-matrix=4 elements
  integer, parameter, public :: EL_NJACENTRIES2D = 4

  ! Number of entries in the Jacobian matrix, defining the mapping between
  ! the reference element and the real element. 3x3-matrix=9 elements
  integer, parameter, public :: EL_NJACENTRIES3D = 9
  
!</constantblock>

!</constants>

  interface elem_generic
    module procedure elem_generic1
    module procedure elem_generic2
  end interface

  interface elem_generic_sim
    module procedure elem_generic_sim1
    module procedure elem_generic_sim2
  end interface

contains

  ! ***************************************************************************

!<function>  

  integer(I32) function elem_igetID(selemName, bcheck)
  
!<description>
  ! This routine returns the element id to a given element name.
  ! It is  case-insensitive. 
!</description>

!<result>
  ! id of the element
!</result>

  !<input>

  ! Element name - one of the EL_xxxx constants.
  character (LEN=*) :: selemName
  
  ! Check the element type. If set to TRUE and the element
  ! name is invalid, the program is not stopped, but 0 is returned.
  logical, intent(in), optional :: bcheck

  !</input>
  
!</function>

    character(len=len(selemName)+1) :: selem
    logical :: bchk
    
    ! SELECT CASE is not allowed for strings (although supported by a majority
    ! of compilers), therefore we have to use a couple of if-commands :(
    ! select case(trim(sys_upcase(scubName)))

    selem = trim(sys_upcase(selemName))

    ! -= 1D Line Elements =-
    if (selem .eq. "EL_P0_1D" .or. selem .eq. "EL_E000_1D") then
      elem_igetID = EL_P0_1D
    else if (selem .eq. "EL_P1_1D" .or. selem .eq. "EL_E001_1D") then
      elem_igetID = EL_P1_1D
    else if (selem .eq. "EL_P2_1D" .or. selem .eq. "EL_E002_1D") then
      elem_igetID = EL_P2_1D
    else if (selem .eq. "EL_S31_1D") then
      elem_igetID = EL_S31_1D
    else if (selem .eq. "EL_DG_T0_1D") then
      elem_igetID = EL_DG_T0_1D
    else if (selem .eq. "EL_DG_T1_1D") then
      elem_igetID = EL_DG_T1_1D
    else if (selem .eq. "EL_DG_T2_1D") then
      elem_igetID = EL_DG_T2_1D            
    
    ! -= 2D Triangle Elements =-
    else if (selem .eq. "EL_P0" .or. selem .eq. "EL_P0_2D" .or. &
             selem .eq. "EL_E000" .or. selem .eq. "EL_E000_2D") then
      elem_igetID = EL_P0_2D
    else if (selem .eq. "EL_P1" .or. selem .eq. "EL_P1_2D" .or. &
             selem .eq. "EL_E001" .or. selem .eq. "EL_E001_2D") then
      elem_igetID = EL_P1_2D
    else if (selem .eq. "EL_P2" .or. selem .eq. "EL_P2_2D" .or. &
             selem .eq. "EL_E002" .or. selem .eq. "EL_E002_2D") then
      elem_igetID = EL_P2_2D
    else if (selem .eq. "EL_P3" .or. selem .eq. "EL_P3_2D" .or. &
             selem .eq. "EL_E003" .or. selem .eq. "EL_E003_2D") then
      elem_igetID = EL_P3_2D
    else if (selem .eq. "EL_P1T" .or. selem .eq. "EL_P1T_2D" .or. &
             selem .eq. "EL_E020" .or. selem .eq. "EL_E020_2D") then
      elem_igetID = EL_P1T_2D
    
    ! -= 2D Quadrilateral Elements =-
    else if (selem .eq. "EL_Q0" .or. selem .eq. "EL_Q0_2D" .or. &
             selem .eq. "EL_E010" .or. selem .eq. "EL_E010_2D") then
      elem_igetID = EL_Q0_2D
    else if (selem .eq. "EL_Q1" .or. selem .eq. "EL_Q1_2D" .or. &
             selem .eq. "EL_E011" .or. selem .eq. "EL_E011_2D") then
      elem_igetID = EL_Q1_2D
    else if (selem .eq. "EL_EM11" .or. selem .eq. "EL_EM11_2D") then
      elem_igetID = EL_EM11_2D
    else if (selem .eq. "EL_Q2" .or. selem .eq. "EL_Q2_2D") then
      elem_igetID = EL_Q2_2D
    else if (selem .eq. "EL_Q2H_2D") then
      elem_igetID = EL_Q2H_2D
    else if (selem .eq. "EL_Q3" .or. selem .eq. "EL_Q3_2D") then
      elem_igetID = EL_Q3_2D
    else if (selem .eq. "EL_QP1" .or. selem .eq. "EL_QP1_2D") then
      elem_igetID = EL_QP1_2D
    else if (selem .eq. "EL_QP1NP" .or. selem .eq. "EL_QP1NP_2D") then
      elem_igetID = EL_QP1NP_2D
    else if (selem .eq. "EL_QP1NPD" .or. selem .eq. "EL_QP1NPD_2D") then
      elem_igetID = EL_QP1NPD_2D
    else if (selem .eq. "EL_Q1T" .or. selem .eq. "EL_Q1T_2D" .or. &
             selem .eq. "EL_E030" .or. selem .eq. "EL_E030_2D") then
      elem_igetID = EL_E030_2D
    else if (selem .eq. "EL_EM30" .or. selem .eq. "EL_EM30_2D") then
      elem_igetID = EL_EM30_2D
    else if (selem .eq. "EL_E031" .or. selem .eq. "EL_E031_2D") then
      elem_igetID = EL_E031_2D
    else if (selem .eq. "EL_EM31" .or. selem .eq. "EL_EM31_2D") then
      elem_igetID = EL_EM31_2D
    else if (selem .eq. "EL_EM30_UNPIVOTED") then
      elem_igetID = EL_EM30_UNPIVOTED
    else if (selem .eq. "EL_EM30_UNSCALED") then
      elem_igetID = EL_EM30_UNSCALED
    else if (selem .eq. "EL_EM30_NEW") then
      elem_igetID = EL_EM30_NEW
    else if (selem .eq. "EL_EM31_NEW") then
      elem_igetID = EL_EM31_NEW
    else if (selem .eq. "EL_Q1TB" .or. selem .eq. "EL_Q1TB_2D" .or. &
             selem .eq. "EL_EB30" .or. selem .eq. "EL_EB30_2D") then
      elem_igetID = EL_Q1TB_2D
    else if (selem .eq. "EL_E032" .or. selem .eq. "EL_E032_2D") then
      elem_igetID = EL_E032_2D
    else if (selem .eq. "EL_Q2T" .or. selem .eq. "EL_Q2T_2D" .or. &
             selem .eq. "EL_E050" .or. selem .eq. "EL_E050_2D") then
      elem_igetID = EL_E050_2D
    else if (selem .eq. "EL_Q2TB" .or. selem .eq. "EL_Q2TB_2D" .or. &
             selem .eq. "EL_EB50" .or. selem .eq. "EL_EB50_2D") then
      elem_igetID = EL_EB50_2D
    else if (selem .eq. "EL_EM50" .or. selem .eq. "EL_EM50_2D") then
      elem_igetID = EL_EM50_2D
    else if (selem .eq. "EL_DG_T0_2D") then
      elem_igetID = EL_DG_T0_2D
    else if (selem .eq. "EL_DG_T1_2D") then
      elem_igetID = EL_DG_T1_2D
    else if (selem .eq. "EL_DG_T2_2D") then
      elem_igetID = EL_DG_T2_2D
    
    ! -= 3D Tetrahedron Elements =-
    else if (selem .eq. "EL_P0_3D" .or. selem .eq. "EL_E000_3D") then
      elem_igetID = EL_P0_3D
    else if (selem .eq. "EL_P1_3D" .or. selem .eq. "EL_E001_3D") then
      elem_igetID = EL_P1_3D
    
    ! -= 3D Hexahedron Elements =-
    else if (selem .eq. "EL_Q0_3D" .or. selem .eq. "EL_E010_3D") then
      elem_igetID = EL_Q0_3D
    else if (selem .eq. "EL_Q1_3D" .or. selem .eq. "EL_E011_3D") then
      elem_igetID = EL_Q1_3D
    else if (selem .eq. "EL_Q2_3D" .or. selem .eq. "EL_E013_3D") then
      elem_igetID = EL_Q2_3D
    else if (selem .eq. "EL_QP1_3D") then
      elem_igetID = EL_QP1_3D
    else if (selem .eq. "EL_QP1NP_3D") then
      elem_igetID = EL_QP1NP_3D
    else if (selem .eq. "EL_Q1T_3D" .or. selem .eq. "EL_E031_3D") then
      elem_igetID = EL_E031_3D
    else if (selem .eq. "EL_E030_3D") then
      elem_igetID = EL_E030_3D
    else if (selem .eq. "EL_EM30_3D") then
      elem_igetID = EL_EM30_3D
    else if (selem .eq. "EL_EM30_NEW_3D") then
      elem_igetID = EL_EM30_NEW_3D
    else if (selem .eq. "EL_Q2T_3D" .or. selem .eq. "EL_E050_3D") then
      elem_igetID = EL_E050_3D
    else if (selem .eq. "EL_EM50_3D") then
      elem_igetID = EL_EM50_3D
    
    ! -= 3D Pyramid Elements =-
    else if (selem .eq. "EL_Y0_3D") then
      elem_igetID = EL_Y0_3D
    else if (selem .eq. "EL_Y1_3D") then
      elem_igetID = EL_Y1_3D
    
    ! -= 3D Prism Elements =-
    else if (selem .eq. "EL_R0_3D") then
      elem_igetID = EL_R0_3D
    else if (selem .eq. "EL_R1_3D") then
      elem_igetID = EL_R1_3D

    ! Unknown element
    else
      bchk = .false.
      if (present(bcheck)) bchk = bcheck
      if (.not. bcheck) then
        call output_line('Error: Unknown element: ' // selemName, &
                        OU_CLASS_ERROR,OU_MODE_STD,'elem_igetID')
        call sys_halt()
      else
        elem_igetID = 0
      end if
      
    end if
  
  end function
  
  ! ***************************************************************************

!<function>  

  elemental integer function elem_igetNDofLoc(celement)

!<description>
  ! This function returns for a given element type the number of local
  ! degrees of freedom.
!</description>

!<input>    

  ! The element type identifier.
  integer(I32), intent(in) :: celement

!</input>

!<result>
  ! The number of local DOF`s on this element.
!</result>

!</function>

    select case (elem_getPrimaryElement(celement))
    
    ! -= 1D Line Elements =-
    case (EL_P0_1D)
      ! local DOFs for 1D P0
      elem_igetNDofLoc = 1
    case (EL_P1_1D)
      ! local DOFs for 1D P1
      elem_igetNDofLoc = 2
    case (EL_P2_1D)
      ! local DOFs for 1D P2
      elem_igetNDofLoc = 3
    case (EL_S31_1D)
      ! local DOFs for 1D S31
      elem_igetNDofLoc = 4
    case (EL_PN_1D)
      ! local DOFs for 1D Pn
      elem_igetNDofLoc = 1 + iand(ishft(celement,-16),255_I32)
    case (EL_DG_T0_1D)
      ! local DOFs for 1D DG Taylor constant
      elem_igetNDofLoc = 1
    case (EL_DG_T1_1D)
      ! local DOFs for 1D DG Taylor linear
      elem_igetNDofLoc = 2
    case (EL_DG_T2_1D)
      ! local DOFs for 1D DG Taylor quadratic
      elem_igetNDofLoc = 3
    
    ! -= 2D Triangle Elements =-
    case (EL_P0)
      ! local DOFs for P0
      elem_igetNDofLoc = 1
    case (EL_P1,EL_DG_P1_2D)
      ! local DOFs for P1
      elem_igetNDofLoc = 3
    case (EL_P2)
      ! local DOFs for P2
      elem_igetNDofLoc = 6
    case (EL_P3)
      ! local DOFs for P3
      elem_igetNDofLoc = 9
    case (EL_P1T)
      ! local DOFs for Ex20
      elem_igetNDofLoc = 3

    ! -= 2D Quadrilateral Elements =-
    case (EL_Q0)
      ! local DOFs for Q0
      elem_igetNDofLoc = 1
    case (EL_Q1)
      ! local DOFs for Q1
      elem_igetNDofLoc = 4
    case (EL_Q2)
      ! local DOFs for Q2
      elem_igetNDofLoc = 9
    case (EL_Q3)
      ! local DOFs for Q3
      elem_igetNDofLoc = 16
    case (EL_QP1)
      ! local DOFs for QP1
      elem_igetNDofLoc = 3
    case (EL_Q1T)
      ! local DOFs for Ex30, Ex31
      elem_igetNDofLoc = 4
    case (EL_Q1TB)
      ! local DOFs for EB30 / E032
      elem_igetNDofLoc = 5
    case (EL_Q2T)
      ! local DOFs for Ex50 
      elem_igetNDofLoc = 9
    case (EL_Q2TB)
      ! local DOFs for EB50 
      elem_igetNDofLoc = 10
    case (EL_DG_T0_2D)
      ! local DOFs for 1D DG Taylor constant
      elem_igetNDofLoc = 1
    case (EL_DG_T1_2D)
      ! local DOFs for 1D DG Taylor linear
      elem_igetNDofLoc = 3
    case (EL_DG_T2_2D)
      ! local DOFs for 1D DG Taylor quadratic
      elem_igetNDofLoc = 6
    
    ! -= 3D Tetrahedron Elements =-
    case (EL_P0_3D)
      ! local DOFs for 3D P0, Q0, Y0, R0
      elem_igetNDofLoc = 1
    case (EL_P1_3D)
      ! local DOFs for 3D P1
      elem_igetNDofLoc = 4
    
    ! -= 3D Hexahedron Elements =-
    case (EL_Q0_3D)
      ! local DOFs for 3D Q0
      elem_igetNDofLoc = 1
    case (EL_Q1_3D)
      ! local DOFs for 3D Q1
      elem_igetNDofLoc = 8
    case (EL_Q2_3D)
      ! local DOFs for 3D Q2
      elem_igetNDofLoc = 27
    case (EL_QP1_3d)
      ! local DOFs for 3D QP1
      elem_igetNDofLoc = 4
    case (EL_Q1T_3D)
      ! local DOFs for 3D Ex30
      elem_igetNDofLoc = 6
    case (EL_Q2T_3D)
      ! local DOFs for 3D Ex50
      elem_igetNDofLoc = 19

    ! -= 3D Pyramid Elements =-
    case (EL_Y0_3D)
      ! local DOFs for 3D Y0
      elem_igetNDofLoc = 1
    case (EL_Y1_3D)
      ! local DOFs for 3D Y1
      elem_igetNDofLoc = 5
    
    ! -= 3D Prism Elements =-
    case (EL_R0_3D)
      ! local DOFs for 3D R0
      elem_igetNDofLoc = 1
    case (EL_R1_3D)
      ! local DOFs for 3D R1
      elem_igetNDofLoc = 6
      
    case default
      elem_igetNDofLoc = 0
    end select

  end function

  ! ***************************************************************************

!<subroutine>  

  pure subroutine elem_igetNDofLocAssignment(celement, &
      ndofAtVertices, ndofAtEdges, ndofAtFaces, ndofAtElement)

!<description>
  ! This function returns for a given element type the number of local
  ! degrees of freedom that is assigned to vertices, edges, faces and elements.
!</description>

!<input>    
  ! The element type identifier.
  integer(I32), intent(in) :: celement
!</input>

!<output>
  ! Number of DOF`s assigned to the vertices on one element.
  integer, intent(out) :: ndofAtVertices
  
  ! Number of DOF`s assigned to the edges on one element.
  ! Is always 0 for 1D element types.
  integer, intent(out) :: ndofAtEdges
  
  ! Number of DOF`s assigned to the faces on one element.
  ! Is always 0 for 1D/2D element types.
  integer, intent(out) :: ndofAtFaces
  
  ! Number of DOF`s assigned to one element, which do not belong to
  ! vertices or edges.
  integer, intent(out) :: ndofAtElement

!</output>

!</subroutine>

    ! Default setup
    ndofAtVertices = 0
    ndofAtEdges    = 0
    ndofAtFaces    = 0
    ndofAtElement  = 0

    select case (elem_getPrimaryElement(celement))
    
    ! -= 1D Line Elements =-
    case (EL_P0_1D)
      ! local DOFs for P0
      ndofAtElement  = 1
    case (EL_P1_1D)
      ! local DOFs for P1
      ndofAtVertices = 2
    case (EL_P2_1D)
      ! local DOFs for P2
      ndofAtVertices = 2
      ndofAtElement  = 1
    case (EL_S31_1D)
      ! local DOFs for S31
      ndofAtVertices = 4
    case (EL_PN_1D)
      ! local DOFs for Pn
      ndofAtVertices = 2
      ndofAtElement = iand(ishft(celement,-16),255_I32)-1
    case (EL_DG_T0_1D)
      ! local DOFs for 1D DG Taylor constant
      ndofAtElement = 1
    case (EL_DG_T1_1D)
      ! local DOFs for 1D DG Taylor linear
      ndofAtElement = 2
    case (EL_DG_T2_1D)
      ! local DOFs for 1D DG Taylor quadratic
      ndofAtElement = 3
    
    ! -= 2D Triangle Elements =-
    case (EL_P0)
      ! local DOFs for Q0
      ndofAtElement  = 1
    case (EL_P1)
      ! local DOFs for P1
      ndofAtVertices = 3
    case (EL_P2)
      ! local DOFs for P2
      ndofAtVertices = 3
      ndofAtEdges    = 3
    case (EL_P3)
      ! local DOFs for P3
      ndofAtVertices = 3
      ndofAtEdges    = 6
    case (EL_P1T)
      ! local DOFs for P1~
      ndofAtEdges    = 3

    case (EL_DG_P1_2D)
      ! local DOFs for P1
      ndofAtElement = 3
    
    ! -= 2D Quadrilateral Elements =-
    case (EL_Q0)
      ! local DOFs for Q0
      ndofAtElement  = 1
    case (EL_Q1)
      ! local DOFs for Q1
      ndofAtVertices = 4
    case (EL_Q2)
      ! local DOFs for Q2
      ndofAtVertices = 4
      ndofAtEdges    = 4
      ndofAtElement  = 1
    case (EL_Q3)
      ! local DOFs for Q3
      ndofAtVertices = 4
      ndofAtEdges    = 8
      ndofAtElement  = 4
    case (EL_QP1)
      ! local DOFs for QP1
      ndofAtElement  = 3
    case (EL_Q1T)
      ! local DOFs for Ex30
      ndofAtEdges    = 4
    case (EL_Q1TB)
      ! local DOFs for EB30 / E032
      ndofAtEdges    = 4
      ndofAtElement  = 1
    case (EL_Q2T)
      ! local DOFs for Ex50
      ndofAtEdges    = 8
      ndofAtElement  = 1
    case (EL_Q2TB)
      ! local DOFs for EB50
      ndofAtEdges    = 8
      ndofAtElement  = 2
    case (EL_DG_T0_2D)
      ! local DOFs for 1D DG Taylor constant
      ndofAtElement = 1
    case (EL_DG_T1_2D)
      ! local DOFs for 1D DG Taylor linear
      ndofAtElement = 3
    case (EL_DG_T2_2D)
      ! local DOFs for 1D DG Taylor quadratic
      ndofAtElement = 6
    
      
    ! -= 3D Tetrahedron Elements =-
    case (EL_P0_3D)
      ! local DOFs for P0,Q0,Y0,R0
      ndofAtElement  = 1
    case (EL_P1_3D)
      ! local DOFs for P1
      ndofAtVertices = 4
    
    ! -= 3D Hexahedron Elements =-
    case (EL_Q0_3D)
      ! local DOFs for Q0
      ndofAtElement  = 1
    case (EL_Q1_3D)
      ! local DOFs for Q1
      ndofAtVertices = 8
    case (EL_Q2_3D)
      ! local DOFs for Q2
      ndofAtVertices = 8
      ndofAtEdges = 12
      ndofAtFaces = 6
      ndofAtElement = 1
    case (EL_QP1_3D)
      ! local DOFs for QP1
      ndofAtElement  = 4
    case (EL_Q1T_3D)
      ! local DOFs for Ex30
      ndofAtFaces = 6
    case (EL_Q2T_3D)
      ! local DOFs for Ex50
      ndofAtFaces = 18
      ndofAtElement = 1
      
    ! -= 3D Pyramid Elements =-
    case (EL_Y0_3D)
      ! local DOFs for Y0
      ndofAtElement  = 1
    case (EL_Y1_3D)
      ! local DOFs for Y1
      ndofAtVertices = 5

    ! -= 3D Pyramid Elements =-
    case (EL_R0_3D)
      ! local DOFs for R0
      ndofAtElement  = 1
    case (EL_R1_3D)
      ! local DOFs for R1
      ndofAtVertices = 6
    end select

  end subroutine

  ! ***************************************************************************

!<function>  

  elemental integer function elem_igetNVE(celement)

!<description>
  ! This function returns for a given element type the number of vertices
  ! in the corresponding element primitive, i.e. 3 for triangular and
  ! 4 for quadrilateral elements.
!</description>

!<input>    

  ! The element type identifier.
  integer(I32), intent(in) :: celement

!</input>

!<result>
  ! The number vertices in the element primitive/shape where the finite
  ! element is defined on.
!</result>

!</function>

    ! NVE coincides with that NVE from the transformation.
    ! Use the transformation routine to determine that value!
    elem_igetNVE = trafo_igetNVE(elem_igetTrafoType(celement))

  end function

  ! ***************************************************************************

!<function>  

  elemental integer(I32) function elem_igetCoordSystem(celement)

!<description>
  ! This function returns for a given element type the type of the coordinate
  ! system. Whenever an element of this type is evaluated, it expects the
  ! coordinates where to evaluate in the coordinates of this coordinate
  ! system. For example, triangular <tex>$P_1$</tex> elements usually work in the
  ! barycentric triangular coordinate system, identified by the coordinate
  ! system identifier TRAFO_CS_BARY2DTRI.
!</description>

!<input>    

  ! The element type identifier.
  integer(I32), intent(in) :: celement

!</input>

!<result>
  ! The type of the coordinate system the element celement acts on. One of
  ! the TRAFO_CS_xxxx constants from the transformation module.
!</result>

!</function>

    select case (iand(celement,EL_ELNRMASK+EL_NONPARAMETRIC))
    ! 1D Element types
    case (EL_P0_1D, EL_P1_1D, EL_P2_1D, EL_S31_1D, EL_PN_1D,&
          EL_DG_T0_1D, EL_DG_T1_1D, EL_DG_T2_1D)
      ! Line elements
      elem_igetCoordSystem = TRAFO_CS_REF1D
    
    ! 2D Element Types
    case (EL_P0, EL_P1, EL_P2, EL_P3, EL_P1T, EL_DG_P1_2D)
      ! Triangular elements work in barycentric coordinates
      elem_igetCoordSystem = TRAFO_CS_BARY2DTRI
    case (EL_Q0, EL_Q1, EL_Q2, EL_Q3, EL_QP1,&
          EL_Q1T, EL_Q1TB, EL_Q2T, EL_Q2TB,&
          EL_DG_T0_2D, EL_DG_T1_2D, EL_DG_T2_2D)
      ! These work on the reference quadrilateral
      elem_igetCoordSystem = TRAFO_CS_REF2DQUAD
    case (EL_QP1NP,EL_QP1NPD)
      ! EL_QP1NP works in real coordinates
      elem_igetCoordSystem = TRAFO_CS_REAL2DQUAD
    case (EL_Q1 +EL_NONPARAMETRIC, &
          EL_Q1T+EL_NONPARAMETRIC, &
          EL_Q2T+EL_NONPARAMETRIC)
      ! EM11, EM30, EM31, EM50; these work in real coordinates
      elem_igetCoordSystem = TRAFO_CS_REAL2DQUAD
    
    ! 3D Element types
    case (EL_P0_3D, EL_P1_3D)
      ! Tetrahedral elements work in barycentric coordinates
      elem_igetCoordSystem = TRAFO_CS_BARY3DTETRA
    case (EL_Q0_3D, EL_Q1_3D, EL_Q2_3D, EL_QP1_3D, &
          EL_Q1T_3D, EL_Q2T_3D)
      ! These work on the reference hexahedron
      elem_igetCoordSystem = TRAFO_CS_REF3DHEXA
    case (EL_Y0_3D, EL_Y1_3D)
      ! These work on the reference pyramid
      elem_igetCoordSystem = TRAFO_CS_REF3DPYRA
    case (EL_R0_3D, EL_R1_3D)
      ! These work on the reference prism
      elem_igetCoordSystem = TRAFO_CS_REF3DPRISM
    case (EL_Q1T_3D+EL_NONPARAMETRIC,&
          EL_Q2T_3D+EL_NONPARAMETRIC)
      ! EM30, EM50; these work in real coordinates
      elem_igetCoordSystem = TRAFO_CS_REAL3DHEXA
    case (EL_QP1NP_3D)
      ! EM30, EM50; these work in real coordinates
      elem_igetCoordSystem = TRAFO_CS_REAL3DHEXA
      
    case default
      elem_igetCoordSystem = TRAFO_CS_UNDEFINED
    end select
  
  end function

  ! ***************************************************************************

!<function>  

  elemental integer(I32) function elem_igetTrafoType(celement)

!<description>
  ! This function returns the typical type of transformation, a specific 
  ! element uses for transferring coordinates of the reference element to 
  ! the real element.
!</description>

!<input>    
  ! The element type identifier.
  integer(I32), intent(in) :: celement
!</input>

!<result>
  ! The maximum derivative of that element. <= MAX_NDER.
!</result>

!</function>

    select case (elem_getPrimaryElement(celement))
    
    case (EL_P0_1D, EL_P1_1D, EL_P2_1D, EL_S31_1D, EL_PN_1D,&
          EL_DG_T0_1D, EL_DG_T1_1D, EL_DG_T2_1D)
      ! Linear line transformation, 1D
      elem_igetTrafoType = TRAFO_ID_MLINCUBE + TRAFO_DIM_1D
    
    case (EL_P0, EL_P1, EL_P2, EL_P3, EL_P1T, EL_DG_P1_2D)
      ! Linear triangular transformation, 2D
      elem_igetTrafoType = TRAFO_ID_LINSIMPLEX + TRAFO_DIM_2D
      
    case (EL_Q0, EL_Q1, EL_Q2, EL_Q3, EL_QP1,&
          EL_Q1T, EL_Q1TB, EL_Q2T, EL_Q2TB,&
          EL_DG_T0_2D, EL_DG_T1_2D, EL_DG_T2_2D)
      ! Bilinear quadrilateral transformation, 2D.
      elem_igetTrafoType = TRAFO_ID_MLINCUBE + TRAFO_DIM_2D
    
    case (EL_P0_3D, EL_P1_3D)
      ! Linear tetrahedral transrormation, 3D
      elem_igetTrafoType = TRAFO_ID_LINSIMPLEX + TRAFO_DIM_3D
    
    case (EL_Q0_3D, EL_Q1_3D, EL_Q2_3D, EL_QP1_3D, &
          EL_Q1T_3D, EL_Q2T_3D)
      ! Trilinear hexahedral transformation, 3D
      elem_igetTrafoType = TRAFO_ID_MLINCUBE + TRAFO_DIM_3D
    
    case (EL_Y0_3D, EL_Y1_3D)
      ! Sub-Trilinear pyramid transformation, 3D
      elem_igetTrafoType = TRAFO_ID_MLINPYRAMID + TRAFO_DIM_3D
    
    case (EL_R0_3D, EL_R1_3D)
      ! Sub-Trilinear prism transformation, 3D
      elem_igetTrafoType = TRAFO_ID_MLINPRISM + TRAFO_DIM_3D
    
    case default
      elem_igetTrafoType = TRAFO_ID_UNKNOWN
      
    end select

  end function

  ! ***************************************************************************

!<function>  

  elemental integer function elem_igetDimension(celement)

!<description>
  ! This function returns the dimensional constant that specifies which
  ! dimension (1D, 2D,...) an element uses.
!</description>

!<input>    
  ! The element type identifier.
  integer(I32), intent(in) :: celement
!</input>

!<result>
  ! A constant that specifies the dimension of an element. NDIM2 for 2D,
  ! NDIM3D for 3D,...
!</result>

!</function>

    ! The dimension is encoded in two bits in the element quantifier!
    elem_igetDimension = ishft(iand(celement,EL_DIMENSION),-8)

  end function

  ! ***************************************************************************

!<function>  

  elemental integer function elem_getMaxDerivative(celement)

!<description>
  ! This function determines for a given element type the maximum derivative
  ! that the element can calculate (One of the DER_xxxx constants and 
  ! <= MAX_NDER). This can be used to specify the size of the array DBas.
!</description>

!<input>    

  ! The element type identifier.
  integer(I32), intent(in) :: celement

!</input>

!<result>
  ! One of the TRAFO_IDxxxx constants of the module 'transformation' identifying
  ! the type of transformation.
!</result>

!</function>

    select case (elem_getPrimaryElement(celement))
    
    ! -= 1D elements =-
    case (EL_P0_1D)
      ! Function
      elem_getMaxDerivative = 1
    case (EL_P1_1D)
      ! Function + 1st derivative
      elem_getMaxDerivative = 2
    case (EL_P2_1D)
      ! Function + 1st derivative
      elem_getMaxDerivative = 2
    case (EL_S31_1D)
      ! Function + 1st derivative
      elem_getMaxDerivative = 2
    case (EL_PN_1D)
      ! Function + 1st derivative
      elem_getMaxDerivative = 2
    case (EL_DG_T0_1D)
      ! Function 
      elem_getMaxDerivative = 1
    case (EL_DG_T1_1D)
      ! Function + 1st derivative
      elem_getMaxDerivative = 2
    case (EL_DG_T2_1D)
      ! Function + 1st derivative + 2nd derivative
      elem_getMaxDerivative = 3
    
    ! -= 2D elements =-
    case (EL_P0, EL_Q0)
      ! Function
      elem_getMaxDerivative = 1
    case (EL_P1,EL_DG_P1_2D)
      ! Function + 1st derivative
      elem_getMaxDerivative = 3
    case (EL_P2)
      ! Function + 1st derivative
      elem_getMaxDerivative = 3
    case (EL_P3)
      ! Function + 1st derivative
      elem_getMaxDerivative = 3
    case (EL_P1T)
      ! Function + 1st derivative
      elem_getMaxDerivative = 3
    case (EL_Q1)
      ! Function + 1st derivative
      elem_getMaxDerivative = 3
    case (EL_Q2)
      ! Function + 1st derivative + 2nd derivative
      elem_getMaxDerivative = 6
    case (EL_Q3)
      ! Function + 1st derivative
      elem_getMaxDerivative = 3
    case (EL_QP1)
      ! Function + 1st derivative
      elem_getMaxDerivative = 3
    case (EL_Q1T,EL_Q1TB)
      ! Function + 1st derivative
      elem_getMaxDerivative = 3
    case (EL_Q2T,EL_Q2TB)
      ! Function + 1st derivative
      elem_getMaxDerivative = 3
    case (EL_DG_T0_2D)
      ! Function 
      elem_getMaxDerivative = 1
    case (EL_DG_T1_2D)
      ! Function + 1st derivative
      elem_getMaxDerivative = 3
    case (EL_DG_T2_2D)
      ! Function + 1st derivative + 2nd derivative
      elem_getMaxDerivative = 6
    
    ! -= 3D elements =-
    case (EL_P0_3D, EL_Q0_3D, EL_Y0_3D, EL_R0_3D)
      ! Function
      elem_getMaxDerivative = 1
    case (EL_P1_3D)
      ! Function + 1st derivative
      elem_getMaxDerivative = 4
    case (EL_Q1_3D)
      ! Function + 1st derivative
      elem_getMaxDerivative = 4
    case (EL_Q2_3D)
      ! Function + 1st derivative
      elem_getMaxDerivative = 4
    case (EL_QP1_3D)
      ! Function + 1st derivative
      elem_getMaxDerivative = 4
    case (EL_Y1_3D)
      ! Function + 1st derivative
      elem_getMaxDerivative = 4
    case (EL_R1_3D)
      ! Function + 1st derivative
      elem_getMaxDerivative = 4
    case (EL_Q1T_3D)
      ! Function + 1st derivative
      elem_getMaxDerivative = 4
    case (EL_Q2T_3D)
      ! Function + 1st derivative
      elem_getMaxDerivative = 4
    case default
      ! We do not know
      elem_getMaxDerivative = DER_MAXNDER
    end select

  end function

  ! ***************************************************************************

!<function>  

  elemental integer(I32) function elem_getEvaluationTag(celement)

!<description>
  ! Returns for a given element its 'evaluation tag'. This tag is a bitfield
  ! that defines for element celement which information have to be prepared
  ! to evaluate it in a set of points.
  !
  ! If more than one finite element has to be evaluated in the same points,
  ! the evaluation tags of all elements under consideration can be combined
  ! using OR. With the resulting tag, a t_evalElementXXXX structure can be 
  ! set up. This structure allows then to evaluate the element(s).
!</description>

!<input>    
  ! The element type identifier.
  integer(I32), intent(in) :: celement
!</input>

!<result>
  ! The 'evaluation tag' of element type celement.
!</result>

!</function>

    select case (elem_getPrimaryElement(celement))

      ! NOTE: Information about the Jacobian may not be necessary.
      !       However, some routines simply crash, if, e.g., the
      !       Jacobian matrix and/or its determinant are not there.

!!$    case (EL_P0_1D, EL_P0, EL_P0_3D)
!!$      ! No information about Jacobian necessary
!!$      elem_getEvaluationTag = 0
    case (EL_Q2T, EL_Q2TB, EL_Q2T_3D)
      ! We need the twist indices.
      elem_getEvaluationTag = EL_EVLTAG_REFPOINTS + &
        EL_EVLTAG_JAC + EL_EVLTAG_DETJ + EL_EVLTAG_TWISTIDX
    case default
      ! Standard evaluation tag. Evaluate reference coordinates 
      ! + Jac + Determinant of Jac for the transformation
      elem_getEvaluationTag = EL_EVLTAG_REFPOINTS + &
          EL_EVLTAG_JAC + EL_EVLTAG_DETJ
    end select
    
    ! Always evaluate the corners of the cells that define the element shape.
    elem_getEvaluationTag = &
      ior(elem_getEvaluationTag,EL_EVLTAG_COORDS)

    if (elem_isNonparametric(celement)) then
      ! Standard nonparametric element.
      ! Calculate real coordinates (for the element)
      elem_getEvaluationTag = &
        ior(elem_getEvaluationTag,EL_EVLTAG_REALPOINTS)
    else
      ! Standard parametric element
      elem_getEvaluationTag = &
        ior(elem_getEvaluationTag,EL_EVLTAG_REFPOINTS)
    end if

  end function

  !************************************************************************
  
!<function>  

  elemental logical function elem_isNonparametric (celement) result (inonpar)

  !<description>
  
  ! Determines whether an element celement is parametric or nonparametric.
  
  !</description>
  
  !<result>
  ! =true, if the element is nonparametric. 
  ! =false, if the element is parametric.
  !</result>
  
  !<input>
  
  ! Element type qualifier.
  integer(I32), intent(in) :: celement
  
  !</input>
 
!</function>
 
    inonpar = iand(celement,EL_NONPARAMETRIC) .ne. 0
  
  end function

  !************************************************************************
  
!<function>  

  elemental integer(I32) function elem_getPrimaryElement (celement) result (iresult)

!<description>
  ! Determines the 'primary' element type identifier. For standard elements
  ! like <tex>$Q_1$</tex>, there is usually celement=elem_getPrimaryElement(celement).
  ! But some elements like <tex>$\tilde Q_1$</tex> have different modifications,
  ! which are encoded in celement. In this case, elem_getPrimaryElement
  ! returns the element identifier of the element without any
  ! modification; this can be seen as 'standard' or 'primary' element
  ! type.
  ! In case of <tex>$\tilde Q_1$</tex> for example, there is
  !    elem_getPrimaryElement(EL_E030) = elem_getPrimaryElement(EL_EM30)
  !  = elem_getPrimaryElement(EL_E031) = elem_getPrimaryElement(EL_EM31)
  !  = EL_Q1T,
  ! which identifies the 'general' <tex>$\tilde Q_1$</tex> element.
!</description>
  
!<result>
  ! The identifier of the 'standard' element, celement refers to.
!</result>
  
  !<input>
  
  ! Element type qualifier.
   integer(I32), intent(in) :: celement
  
  !</input>
 
!</function>
 
    ! To get the standard identifier, we just have to mask out all bits
    ! except for bit 0..9. These 10 bits encode the standard element
    ! identifier plus dimension.
    iresult = iand(celement,EL_ELNRMASK)
  
  end function

  !************************************************************************
  
!<function>  

  integer(I32) function elem_igetShape(celement) result(ishp)
  
!<description>
  ! This function returns the element shape identifier for a given element
  ! id. The shape identifier is one of the BGEOM_SHAPE_XXXX constants
  ! defined in basicgeometry.f90.
!</description>


!<result>
  ! One of the BGEOM_SHAPE_XXXX shape identifiers.
!</result>

!<input>
  ! Element identifier
  integer(I32), intent(in) :: celement
!</input>
  
!</function>

    select case (elem_getPrimaryElement(celement))
    case (EL_P0_1D, EL_P1_1D, EL_P2_1D, EL_S31_1D, EL_PN_1D,&
          EL_DG_T0_1D, EL_DG_T1_1D, EL_DG_T2_1D)
      ! 1D Line
      ishp = BGEOM_SHAPE_LINE
    
    case (EL_P0, EL_P1, EL_P2, EL_P3, EL_P1T, EL_DG_P1_2D)
      ! 2D Triangle
      ishp = BGEOM_SHAPE_TRIA
      
    case (EL_Q0, EL_Q1, EL_Q2, EL_Q3, EL_QP1,&
          EL_Q1T, EL_Q1TB, EL_Q2T, EL_Q2TB,&
          EL_DG_T0_2D, EL_DG_T1_2D, EL_DG_T2_2D)
      ! 2D Quadrilateral
      ishp = BGEOM_SHAPE_QUAD
    
    case (EL_P0_3D, EL_P1_3D, EL_P2_3D)
      ! 3D Tetrahedron
      ishp = BGEOM_SHAPE_TETRA
    
    case (EL_Q0_3D, EL_Q1_3D, EL_Q2_3D, EL_QP1_3D, &
          EL_Q1T_3D, EL_Q2T_3D)
      ! 3D Hexahedron
      ishp = BGEOM_SHAPE_HEXA
    
    case (EL_Y0_3D, EL_Y1_3D)
      ! 3D Pyramid
      ishp = BGEOM_SHAPE_PYRA
    
    case (EL_R0_3D, EL_R1_3D)
      ! 3D Prism
      ishp = BGEOM_SHAPE_PRISM
    
    case default
      ! Unknown
      ishp = BGEOM_SHAPE_UNKNOWN
      
    end select

  end function

!**************************************************************************
! Generic element routines
! The following two routines define a generic element. Depending on
! the element type identifier celement, the right element is called.
 
!<subroutine>  

  subroutine elem_generic1 (celement, Dcoords, Djac, ddetj, Bder, &
                                Dpoint, Dbas, ItwistIndex)

!<description>
  ! DEPRECATED!!!
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  ! celement defines the element type that is used. Depending on
  ! the type of the element (parametric or nonparametric), dx
  ! and dy must be given either on the reference element or on the
  ! real element.
!</description>

!<input>
  ! Element type identifier.
  integer(I32), intent(in)  :: celement
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates,
  ! Dcoords(3,.)=z-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords
  
  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. The size of this array depends
  ! on the dimension of the space.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  real(DP), dimension(:), intent(in) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element.
  real(DP), intent(in) :: ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder
  
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
  real(DP), dimension(:), intent(in) :: Dpoint

  ! OPTIONAL: Twist index bitfield that defines the orientation of the edges
  ! of the element.
  ! Can be omitted if the element does not use this information.
  integer(I32), intent(in), optional :: ItwistIndex
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i-th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i-th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  real(DP), dimension(:,:), intent(out) :: Dbas
!</output>

! </subroutine>
  
  ! local variables for the sim2-wrapper
  logical :: bwrapSim2
  real(DP), dimension(size(Dcoords,1),size(Dcoords,2),1), target :: Dcoords2
  real(DP), dimension(size(Djac,1),1,1), target :: Djac2
  real(DP), dimension(1,1), target :: Ddetj2
  real(DP), dimension(size(Dpoint,1),1,1), target :: Dpoints2
  integer(I32), dimension(1), target :: ItwistIndex2
  real(DP), dimension(size(Dbas,1),size(Dbas,2),1,1) :: Dbas2
  type(t_evalElementSet) :: reval

    ! Take care of the 1D PN element
    if(elem_getPrimaryElement(celement) .eq. EL_PN_1D) then
      
      bwrapSim2 = .true.
    
    else
    
      ! Choose the right element subroutine to call.
      bwrapSim2 = .false.
      select case (celement)
      ! 1D elements
      case (EL_P0_1D)
        call elem_P0_1D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_P1_1D)
        call elem_P1_1D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_P2_1D)
        call elem_P2_1D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_S31_1D)
        call elem_S31_1D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_DG_T0_1D)
        call elem_DG_T0_1D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_DG_T1_1D)
        call elem_DG_T1_1D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_DG_T2_1D)
        call elem_DG_T2_1D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)        

      ! 2D elements
      case (EL_P0)
        call elem_P0 (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_P1,EL_DG_P1_2D)
        call elem_P1 (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_P2)
        call elem_P2 (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_P1T)
        call elem_P1T (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_Q0)
        call elem_Q0 (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_Q1)
        call elem_Q1 (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_EM11_2D)
        bwrapSim2 = .true.
      case (EL_Q2)
        call elem_Q2 (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_Q2H_2D)
        bwrapSim2 = .true.
      case (EL_QP1)
        call elem_QP1 (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_QP1NP,EL_QP1NPD)
        call elem_QP1NP (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_EM30,EL_EM30_UNPIVOTED,EL_EM30_UNSCALED)
        call elem_EM30 (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_EM30_NEW)
        bwrapSim2 = .true.
      case (EL_E030)
        call elem_E030 (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_EB30)
        call elem_EB30 (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_EM31)
        call elem_EM31 (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_EM31_NEW)
        bwrapSim2 = .true.
      case (EL_E031)
        call elem_E031 (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_E032)
        bwrapSim2 = .true.
      case (EL_E050)
        call elem_E050 (celement, Dcoords, ItwistIndex, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_EB50)
        call elem_EB50 (celement, Dcoords, ItwistIndex, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_EM50)
        bwrapSim2 = .true.
      case (EL_DG_T0_2D)
        call elem_DG_T0_2D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_DG_T1_2D)
        call elem_DG_T1_2D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_DG_T2_2D)
        call elem_DG_T2_2D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas) 

      ! 3D elements
      case (EL_P0_3D)
        call elem_P0_3D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_P1_3D)
        call elem_P1_3D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_Q0_3D)
        call elem_Q0_3D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_Q1_3D)
        call elem_Q1_3D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_Q2_3D)
        bwrapSim2 = .true.
      case (EL_QP1_3D,EL_QP1NP_3D)
        bwrapSim2 = .true.
      case (EL_Y0_3D)
        call elem_Y0_3D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_Y1_3D)
        call elem_Y1_3D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_R0_3D)
        call elem_R0_3D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_R1_3D)
        call elem_R1_3D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_EM30_3D)
        call elem_EM30_3D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_EM30_NEW_3D)
        bwrapSim2 = .true.
      case (EL_E030_3D)
        call elem_E030_3D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_E031_3D)
        call elem_E031_3D (celement, Dcoords, Djac, ddetj, Bder, Dpoint, Dbas)
      case (EL_EM50_3D)
        bwrapSim2 = .true.
      case (EL_E050_3D)
        bwrapSim2 = .true.

      case default
        ! Element not implemened!
        ! Throw a floating point exception so that the program stops here!
        ! We cannot use "PRINT" here as the routine is PURE!
        call sys_throwFPE()
        
      end select
    
    end if
    
    if(.not. bwrapSim2) return
    
    ! Prepare the variables for the sim2-wrapper
    Dcoords2(:,:,1) = Dcoords(:,:)
    Djac2(:,1,1) = Djac(:)
    Ddetj2(1,1) = ddetj
    Dpoints2(:,1,1) = Dpoint
    if(present(itwistIndex)) then
      ItwistIndex2(1) = itwistIndex
    else
      ItwistIndex2(1) = 0_I32
    end if
    
    ! Set up the structure
    reval%npointsPerElement = 1
    reval%nelements = 1
    reval%p_Dcoords => Dcoords2
    reval%p_Djac => Djac2
    reval%p_Ddetj => Ddetj2
    if(iand(celement,EL_NONPARAMETRIC) .ne. 0) then
      reval%p_DpointsReal => Dpoints2
    else
      reval%p_DpointsRef => Dpoints2
    end if
    reval%p_ItwistIndex => ItwistIndex2
    
    ! Call sim2-wrapper
    call elem_generic_sim2(celement, reval, Bder, Dbas2)
    
    ! Copy results to Dbas
    Dbas(:,:) = Dbas2(:,:,1,1)

  end subroutine 
  
!**************************************************************************
! Generic element routines
! The following two routines define a generic element. Depending on
! the element type identifier celement, the right element is called.
 
!<subroutine>  

  subroutine elem_generic2 (celement, revalElement, Bder, Dbas)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at the given point on the reference element. 
  ! celement defines the element type that is used. Depending on
  ! the type of the element (parametric or nonparametric), dx
  ! and dy must be given either on the reference element or on the
  ! real element.
!</description>

!<input>
  ! Element type identifier.
  integer(I32), intent(in)  :: celement
  
  ! t_evalElement-structure that contains cell-specific information and
  ! coordinates of the evaluation point. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElement), intent(in) :: revalElement

  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder
  
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC) defines the value of the i-th 
  !   basis function of the finite element in the point (dx,dy) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i-th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx) is undefined.
  real(DP), dimension(:,:), intent(out) :: Dbas
!</output>

! </subroutine>

    ! Choose the right element subroutine to call.
    select case (celement)
    ! 1D elements
    case (EL_P0_1D)
      call elem_P0_1D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_P1_1D)
      call elem_P1_1D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_P2_1D)
      call elem_P2_1D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_S31_1D)
      call elem_S31_1D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_DG_T0_1D)
      call elem_DG_T0_1D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_DG_T1_1D)
      call elem_DG_T1_1D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_DG_T2_1D)
      call elem_DG_T2_1D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)

    ! 2D elements
    case (EL_P0)
      call elem_P0 (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_P1,EL_DG_P1_2D)
      call elem_P1 (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_P2)
      call elem_P2 (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_P1T)
      call elem_P1T (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_Q0)
      call elem_Q0 (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_Q1)
      call elem_Q1 (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_Q2)
      call elem_Q2 (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_QP1)
      call elem_QP1 (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_QP1NP,EL_QP1NPD)
      call elem_QP1NP (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointReal, Dbas)
    case (EL_EM30, EL_EM30_UNPIVOTED, EL_EM30_UNSCALED)
      call elem_EM30 (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointReal, Dbas)
    case (EL_E030)
      call elem_E030 (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_EB30)
      call elem_EB30 (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_EM31)
      call elem_EM31 (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointReal, Dbas)
    case (EL_E031)
      call elem_E031 (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_E050)
      call elem_E050 (celement, revalElement%Dcoords, revalElement%itwistIndex, &
          revalElement%Djac, revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_EB50)
      call elem_EB50 (celement, revalElement%Dcoords, revalElement%itwistIndex, &
          revalElement%Djac, revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_DG_T0_2D)
      call elem_DG_T0_2D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_DG_T1_2D)
      call elem_DG_T1_2D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_DG_T2_2D)
      call elem_DG_T2_2D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)

    ! 3D elements
    case (EL_P0_3D)
      call elem_P0_3D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_P1_3D)
      call elem_P1_3D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_Q0_3D)
      call elem_Q0_3D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_Q1_3D)
      call elem_Q1_3D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_Y0_3D)
      call elem_Y0_3D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_Y1_3D)
      call elem_Y1_3D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_R0_3D)
      call elem_R0_3D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_R1_3D)
      call elem_R1_3D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_EM30_3D)
      call elem_EM30_3D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_E030_3D)
      call elem_E030_3D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)
    case (EL_E031_3D)
      call elem_E031_3D (celement, revalElement%Dcoords, revalElement%Djac, &
          revalElement%ddetj, Bder, revalElement%DpointRef, Dbas)

    case default
      ! Element not implemened!
      ! Throw a floating point exception so that the program stops here!
      ! We cannot use "PRINT" here as the routine is PURE!
      call sys_throwFPE()
    end select

  end subroutine 
  
  !************************************************************************
  
!<subroutine>  

  subroutine elem_generic_mult (celement, Dcoords, Djac, Ddetj, &
                                     Bder, Dbas, npoints, Dpoints, itwistIndex)

!<description>
  ! This subroutine calculates the values of the basic functions of the
  ! finite element at multiple given points on the reference element. 
!</description>

!<input>
  ! Element type identifier
  integer(I32), intent(in)  :: celement
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE)
  ! Dcoords(1,.)=x-coordinates,
  ! Dcoords(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real element. For every point i:
  !  Djac(1,i) = J_i(1,1)
  !  Djac(2,i) = J_i(2,1)
  !  Djac(3,i) = J_i(1,2)
  !  Djac(4,i) = J_i(2,2)
  !REAL(DP), DIMENSION(:,npoints), INTENT(in) :: Djac
  real(DP), dimension(:,:), intent(in) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! element for every of the npoints points:
  real(DP), dimension(:), intent(in) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! For parametric elements, the coordinates are expected on the 
  ! reference element. For nonparametric elements, the coordinates
  ! are expected on the real element!
  ! DIMENSION(#space dimensions,npoints)
  !   Dpoints(1,.)=x-coordinates,
  !   Dpoints(2,.)=y-coordinates.
  real(DP), dimension(:,:), intent(in) :: Dpoints
  
  ! OPTIONAL: Twist index bitfield that defines the orientation
  ! of the edges in the element.
  ! Can be omitted if the element does not use this information.
  integer(I32), intent(in), optional :: itwistIndex
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i-th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i-th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:), intent(out) :: Dbas
!</output>

! </subroutine>

  ! local variables
  integer :: i

    ! Choose the right element subroutine to call.
    select case (celement)
    case (EL_P0_1D)
      call elem_P0_1D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_P1_1D)
      call elem_P1_1D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_P2_1D)
      call elem_P2_1D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_S31_1D)
      call elem_S31_1D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_DG_T0_1D)
      call elem_DG_T0_1D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_DG_T1_1D)
      call elem_DG_T1_1D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_DG_T2_1D)
      call elem_DG_T2_1D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
      
    ! 2D elements
    case (EL_P0)
      call elem_P0_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_P1,EL_DG_P1_2D)
      call elem_P1_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_P2)
      call elem_P2_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_P1T)
      call elem_P1T_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_Q0)
      call elem_Q0_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_Q1)
      call elem_Q1_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_EM30, EL_EM30_UNPIVOTED ,EL_EM30_UNSCALED)
      call elem_EM30_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_E030)
      call elem_E030_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_EB30)
      call elem_EB30_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_EM31)
      call elem_EM31_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_E031)
      call elem_E031_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_E050)
      call elem_E050_mult (celement, Dcoords, itwistIndex,Djac, Ddetj, Bder, Dbas, &
                           npoints, Dpoints)
    case (EL_EB50)
      call elem_EB50_mult (celement, Dcoords, itwistIndex,Djac, Ddetj, Bder, Dbas, &
                           npoints, Dpoints)
    case (EL_DG_T0_2D)
      call elem_DG_T0_2D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, &
                               npoints, Dpoints)
    case (EL_DG_T1_2D)
      call elem_DG_T1_2D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, &
                               npoints, Dpoints)
    case (EL_DG_T2_2D)
      call elem_DG_T2_2D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, &
                               npoints, Dpoints)

    ! 3D elements
    case (EL_P0_3D)
      call elem_P0_3D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_P1_3D)
      call elem_P1_3D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_Q0_3D)
      call elem_Q0_3D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_Q1_3D)
      call elem_Q1_3D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_Y0_3D)
      call elem_Y0_3D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_Y1_3D)
      call elem_Y1_3D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_R0_3D)
      call elem_R0_3D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_R1_3D)
      call elem_R1_3D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_EM30_3D)
      call elem_EM30_3D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_E030_3D)
      call elem_E030_3D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)
    case (EL_E031_3D)
      call elem_E031_3D_mult (celement, Dcoords, Djac, Ddetj, Bder, Dbas, npoints, Dpoints)

    case default
      ! Compatibility handling: evaluate all points separately
      do i=1,npoints
        call elem_generic (celement, Dcoords, Djac(:,i), Ddetj(i), Bder, &
                           Dpoints(:,i), Dbas(:,:,i), ItwistIndex)
      end do
    end select

  end subroutine 
  
  !************************************************************************
  
!<subroutine>  

  subroutine elem_generic_sim1 (celement, Dcoords, Djac, Ddetj, &
                               Bder, Dbas, npoints, nelements, Dpoints, ItwistIndex)

!<description>
  ! DEPRECATED:
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! Element type identifier
  integer(I32), intent(in)  :: celement
  
  ! Number of points on every element where to evalate the basis functions.
  integer, intent(in) :: npoints
  
  ! Number of elements, the basis functions are evaluated at
  integer, intent(in)  :: nelements
  
  ! Array with coordinates of the corners that form the real element.
  ! DIMENSION(#space dimensions,NVE,nelements)
  !  Dcoords(1,.,.)=x-coordinates,
  !  Dcoords(2,.,.)=y-coordinates.
  ! furthermore:
  !  Dcoords(:,i,.) = Coordinates of vertex i
  ! furthermore:
  !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
  real(DP), dimension(:,:,:), intent(in) :: Dcoords

  ! Values of the Jacobian matrix that defines the mapping between the
  ! reference element and the real elements. For every point i:
  !  Djac(1,i,.) = J_i(1,1,.)
  !  Djac(2,i,.) = J_i(2,1,.)
  !  Djac(3,i,.) = J_i(1,2,.)
  !  Djac(4,i,.) = J_i(2,2,.)
  ! Remark: Only used for calculating derivatives; can be set to 0.0
  ! when derivatives are not used.
  !  Djac(:,:,j) refers to the determinants of the points of element j.
  !REAL(DP), DIMENSION(:,npoints,nelements), INTENT(in) :: Djac
  real(DP), dimension(:,:,:), intent(in) :: Djac
  
  ! Determinant of the mapping from the reference element to the real
  ! elements for every of the npoints points on all the elements.
  !  Ddetj(i,.) = Determinant of point i
  !  Ddetj(:,j) = determinants of all points on element j
  ! Remark: Only used for calculating derivatives; can be set to 1.0
  ! when derivatives are not needed. Must not be set to 0.0!
  real(DP), dimension(:,:), intent(in) :: Ddetj
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder
  
  ! Array with coordinates of the points where to evaluate.
  ! The coordinates are expected 
  ! - on the reference element, if celement identifies a parametric element
  ! - on the real element, if celement identifies a nonparametric element
  ! It is assumed that:
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(in) :: Dpoints
  
  ! OPTIONAL: List of twist indices. For every element, the corresponding
  ! entry is a bitfield that defines the orientation of the edge.
  ! Can be omitted if the element does not need it.
  ! Array with DIMENSION(nelements)
  integer(I32), dimension(:), intent(in), optional :: ItwistIndex
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! array [1..EL_MAXNBAS,1..DER_MAXNDER,1..npoints,nelements] of double
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i-th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i-th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:,:), intent(out) :: Dbas
!</output>

! </subroutine>

  integer :: i

    ! Choose the right element subroutine to call.
    select case (celement)
    ! 1D elements
    case (EL_P0_1D)
      call elem_P0_1D_sim (celement, Dcoords, Djac, Ddetj, &
                           Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_P1_1D)
      call elem_P1_1D_sim (celement, Dcoords, Djac, Ddetj, &
                           Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_P2_1D)
      call elem_P2_1D_sim (celement, Dcoords, Djac, Ddetj, &
                           Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_S31_1D)
      call elem_S31_1D_sim (celement, Dcoords, Djac, Ddetj, &
                            Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_DG_T0_1D)
      call elem_DG_T0_1D_sim (celement, Dcoords, Djac, Ddetj, &
                              Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_DG_T1_1D)
      call elem_DG_T1_1D_sim (celement, Dcoords, Djac, Ddetj, &
                              Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_DG_T2_1D)
      call elem_DG_T2_1D_sim (celement, Dcoords, Djac, Ddetj, &
                              Bder, Dbas, npoints, nelements, Dpoints)

    ! 2D elements
    case (EL_P0)
      call elem_P0_sim (celement, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_P1,EL_DG_P1_2D)
      call elem_P1_sim (celement, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_P2)
      call elem_P2_sim (celement, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_P1T)
      call elem_P1T_sim (celement, Dcoords, Djac, Ddetj, &
                         Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_Q0)
      call elem_Q0_sim (celement, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_Q1)
      call elem_Q1_sim (celement, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_Q2)
      call elem_Q2_sim (celement, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_EM30, EL_EM30_UNPIVOTED, EL_EM30_UNSCALED)
      call elem_EM30_sim (celement, Dcoords, Djac, Ddetj, &
                          Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_E030)
      call elem_E030_sim (celement, Dcoords, Djac, Ddetj, &
                          Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_EB30)
      call elem_EB30_sim (celement, Dcoords, Djac, Ddetj, &
                          Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_EM31)
      call elem_EM31_sim (celement, Dcoords, Djac, Ddetj, &
                          Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_E031)
      call elem_E031_sim (celement, Dcoords, Djac, Ddetj, &
                          Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_E050)
      call elem_E050_sim (celement, Dcoords, ItwistIndex, Djac, Ddetj, &
                          Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_EB50)
      call elem_EB50_sim (celement, Dcoords, ItwistIndex, Djac, Ddetj, &
                          Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_DG_T0_2D)
      call elem_DG_T0_2D_sim (celement, Dcoords, Djac, Ddetj, &
                              Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_DG_T1_2D)
      call elem_DG_T1_2D_sim (celement, Dcoords, Djac, Ddetj, &
                              Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_DG_T2_2D)
      call elem_DG_T2_2D_sim (celement, Dcoords, Djac, Ddetj, &
                              Bder, Dbas, npoints, nelements, Dpoints)

    ! 3D elements
    case (EL_P0_3D)
      call elem_P0_3D_sim (celement, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_P1_3D)
      call elem_P1_3D_sim (celement, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_Q0_3D)
      call elem_Q0_3D_sim (celement, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_Q1_3D)
      call elem_Q1_3D_sim (celement, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_Y0_3D)
      call elem_Y0_3D_sim (celement, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_Y1_3D)
      call elem_Y1_3D_sim (celement, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_R0_3D)
      call elem_R0_3D_sim (celement, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_R1_3D)
      call elem_R1_3D_sim (celement, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_EM30_3D)
      call elem_EM30_3D_sim (celement, Dcoords, Djac, Ddetj, &
                          Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_E030_3D)
      call elem_E030_3D_sim (celement, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)
    case (EL_E031_3D)
      call elem_E031_3D_sim (celement, Dcoords, Djac, Ddetj, &
                        Bder, Dbas, npoints, nelements, Dpoints)

    case default
      ! Compatibility handling: evaluate on all elements separately
      if (present(ItwistIndex)) then
        do i=1,nelements
          call elem_generic_mult (celement, Dcoords(:,:,i),&
              Djac(:,:,i), Ddetj(:,i), &
              Bder, Dbas(:,:,:,i), npoints, Dpoints(:,:,i), ItwistIndex(i))
        end do
      else
        do i=1,nelements
          call elem_generic_mult (celement, Dcoords(:,:,i),&
              Djac(:,:,i), Ddetj(:,i), &
              Bder, Dbas(:,:,:,i), npoints, Dpoints(:,:,i))
        end do
      end if
    end select

  end subroutine 
  
  !************************************************************************
  
!<subroutine>  

  subroutine elem_generic_sim2 (celement, revalElementSet, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic 
  ! functions of the finite element at multiple given points on the reference 
  ! element for multiple given elements.
!</description>

!<input>
  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in) :: revalElementSet

  ! Element type identifier
  integer(I32), intent(in)  :: celement
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder  
!</input>
  
!<output>
  ! Value/derivatives of basis functions. 
  ! array [1..EL_MAXNBAS,1..DER_MAXNDER,1..npointsPerElement,nelements] of double
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i-th 
  !   basis function of the finite element in the point Dcoords(j) on the 
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i-th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:,:), intent(out) :: Dbas
!</output>

! </subroutine>

  integer :: i
  
    ! Take care of the 1D PN element
    if(elem_getPrimaryElement(celement) .eq. EL_PN_1D) then
      call elem_eval_PN_1D(celement, revalElementSet, Bder, Dbas)
      return
    end if

    ! Choose the right element subroutine to call.
    select case (celement)
    
    ! *****************************************************
    ! 1D line elements
    case (EL_P0_1D)
      call elem_P0_1D_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
    
    case (EL_P1_1D)
      call elem_P1_1D_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
    
    case (EL_P2_1D)
      call elem_P2_1D_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
    
    case (EL_S31_1D)
      call elem_S31_1D_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
        
    case (EL_DG_T0_1D)
      call elem_DG_T0_1D_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
                              
    case (EL_DG_T1_1D)
      call elem_DG_T1_1D_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
                              
    case (EL_DG_T2_1D)
      call elem_DG_T2_1D_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
    
    ! *****************************************************
    ! 2D triangle elements
    case (EL_P0)
      call elem_P0_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
    
    case (EL_P1,EL_DG_P1_2D)
      call elem_P1_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
    
    case (EL_P2)
      call elem_P2_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
    
    case (EL_P1T)
      call elem_P1T_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
    
    ! *****************************************************
    ! 2D quadrilateral elements
    case (EL_Q0)
      call elem_Q0_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
    
    case (EL_Q1)
      call elem_Q1_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
      
      ! New implementation
      !call elem_eval_Q1_2D(celement, revalElementSet, Bder, Dbas)
    
    case (EL_EM11_2D)
      call elem_eval_EM11_2D(celement, revalElementSet, Bder, Dbas)
    
    case (EL_Q2)
      call elem_Q2_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)

    case (EL_Q2H_2D)
      call elem_eval_Q2H_2D(celement, revalElementSet, Bder, Dbas)
    
    case (EL_EM30, EL_EM30_UNPIVOTED, EL_EM30_UNSCALED)
      call elem_EM30_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsReal)
    
    case (EL_EM30_NEW)
      ! new implementation of 2D EM30
      call elem_eval_EM30_2D(celement, revalElementSet, Bder, Dbas)
    
    case (EL_E030)
      call elem_E030_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
    
    case (EL_EB30)
      call elem_EB30_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)

    case (EL_EM31)
      call elem_EM31_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsReal)
    
    case (EL_EM31_NEW)
      ! new implementation of 2D EM31
      call elem_eval_EM31_2D(celement, revalElementSet, Bder, Dbas)
    
    case (EL_E031)
      call elem_E031_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
    
    case (EL_E032)
      call elem_eval_E032_2D(celement, revalElementSet, Bder, Dbas)
    
    case (EL_E050)
      call elem_E050_sim (celement, revalElementSet%p_Dcoords,&
        revalElementSet%p_ItwistIndex, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
    
    case (EL_EB50)
      !call elem_EB50_sim (celement, revalElementSet%p_Dcoords,&
      !  revalElementSet%p_ItwistIndex, &
      !  revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
      !  Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
      !  revalElementSet%p_DpointsRef)

      ! New implementation
      call elem_eval_EB50_2D(celement, revalElementSet, Bder, Dbas)
    
    case (EL_EM50)
      call elem_eval_EM50_2D(celement, revalElementSet, Bder, Dbas)
         
    case (EL_DG_T0_2D)
      call elem_DG_T0_2D_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
                              
    case (EL_DG_T1_2D)
      call elem_DG_T1_2D_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
                              
    case (EL_DG_T2_2D)
      call elem_DG_T2_2D_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)

    ! *****************************************************
    ! 3D tetrahedron elements
    case (EL_P0_3D)
      call elem_P0_3D_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
    
    case (EL_P1_3D)
      call elem_P1_3D_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)

    ! *****************************************************
    ! 3D hexahedron elements
    case (EL_Q0_3D)
      call elem_Q0_3D_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
    
    case (EL_Q1_3D)
      call elem_Q1_3D_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)

    case (EL_Q2_3D)
      call elem_eval_Q2_3D(celement, revalElementSet, Bder, Dbas)
    
    case (EL_QP1_3D)
      call elem_eval_QP1_3D(celement, revalElementSet, Bder, Dbas)

    case (EL_QP1NP_3D)
      call elem_eval_QP1NP_3D(celement, revalElementSet, Bder, Dbas)

    case (EL_EM30_3D)
      call elem_EM30_3D_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsReal)

    case (EL_EM30_NEW_3D)
      ! new implementation of 3D EM30
      call elem_eval_EM30_3D(celement, revalElementSet, Bder, Dbas)
    
    case (EL_E030_3D)
      call elem_E030_3D_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
    
    case (EL_E031_3D)
      call elem_E031_3D_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
        
    case (EL_E050_3D)
      call elem_eval_E050_3D(celement, revalElementSet, Bder, Dbas)
        
    case (EL_EM50_3D)
      call elem_eval_EM50_3D(celement, revalElementSet, Bder, Dbas)

    ! *****************************************************
    ! 3D pyramid elements
    case (EL_Y0_3D)
      call elem_Y0_3D_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
    
    case (EL_Y1_3D)
      !call elem_Y1_3D_sim (celement, revalElementSet%p_Dcoords, &
      !  revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
      !  Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
      !  revalElementSet%p_DpointsRef)

      ! New implementation
      call elem_eval_Y1_3D(celement, revalElementSet, Bder, Dbas)

    ! *****************************************************
    ! 3D prism elements
    case (EL_R0_3D)
      call elem_R0_3D_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
    
    case (EL_R1_3D)
      call elem_R1_3D_sim (celement, revalElementSet%p_Dcoords, &
        revalElementSet%p_Djac, revalElementSet%p_Ddetj, &
        Bder, Dbas, revalElementSet%npointsPerElement, revalElementSet%nelements, &
        revalElementSet%p_DpointsRef)
    
    
          
    case default
      ! Compatibility handling: evaluate on all elements separately
      if (associated(revalElementSet%p_ItwistIndex)) then
        if (elem_isNonparametric(celement)) then
          do i=1,revalElementSet%nelements
            call elem_generic_mult (celement, revalElementSet%p_Dcoords(:,:,i),&
                revalElementSet%p_Djac(:,:,i), revalElementSet%p_Ddetj(:,i), &
                Bder, Dbas(:,:,:,i), revalElementSet%npointsPerElement, &
                revalElementSet%p_DpointsReal(:,:,i), revalElementSet%p_ItwistIndex(i))
          end do
        else
          do i=1,revalElementSet%nelements
            call elem_generic_mult (celement, revalElementSet%p_Dcoords(:,:,i),&
                revalElementSet%p_Djac(:,:,i), revalElementSet%p_Ddetj(:,i), &
                Bder, Dbas(:,:,:,i), revalElementSet%npointsPerElement, &
                revalElementSet%p_DpointsRef(:,:,i), revalElementSet%p_ItwistIndex(i))
          end do
        end if
      else
        if (elem_isNonparametric(celement)) then
          do i=1,revalElementSet%nelements
            call elem_generic_mult (celement, revalElementSet%p_Dcoords(:,:,i),&
                revalElementSet%p_Djac(:,:,i), revalElementSet%p_Ddetj(:,i), &
                Bder, Dbas(:,:,:,i), revalElementSet%npointsPerElement, &
                revalElementSet%p_DpointsReal(:,:,i))
          end do
        else
          do i=1,revalElementSet%nelements
            call elem_generic_mult (celement, revalElementSet%p_Dcoords(:,:,i),&
                revalElementSet%p_Djac(:,:,i), revalElementSet%p_Ddetj(:,i), &
                Bder, Dbas(:,:,:,i), revalElementSet%npointsPerElement, &
                revalElementSet%p_DpointsRef(:,:,i))
          end do
        end if
      end if
    end select

  end subroutine 
   
end module 
